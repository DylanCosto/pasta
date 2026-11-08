mod block;
pub(crate) mod compression_header;
pub(crate) mod embedded_reference;
mod header;
pub(crate) mod slice;

use std::{
    cmp,
    io::{self, Write},
    mem,
    sync::Arc,
};

use crate::sam;

pub use self::{
    block::{Block, write_block},
    header::write_header,
};
use self::{
    compression_header::{build_compression_header, write_compression_header},
    slice::{ExternalBlock, Slice, build_slice},
};
use super::{
    Context, Record, ReferenceBases,
    compression::{Choice, Compressor, Method, QualityLayout, Trial},
};
use crate::cram::container::{Header, ReferenceSequenceContext, block::ContentType};

/// Removes the MD and NM tags that readers regenerate from the external reference, as HTSlib
/// does (see [`Record::remove_regenerated_md_nm`]), unless the options keep them.
fn remove_regenerated_md_nm(
    ctx: &Context,
    header: &sam::Header,
    records: &mut [Record],
) -> io::Result<()> {
    let (store_md, store_nm) = ctx.stores_md_nm();

    if store_md && store_nm {
        return Ok(());
    }

    let mut reference_sequence = None;
    // The MD string of each record in turn.
    let mut md = Vec::new();

    for record in records
        .iter_mut()
        .filter(|record| record.reference_bases == ReferenceBases::External)
    {
        let Some(id) = record.reference_sequence_id else {
            continue;
        };

        if reference_sequence
            .as_ref()
            .is_none_or(|(prev_id, _, _)| *prev_id != id)
        {
            let (name, map) = header.reference_sequences().get_index(id).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
            })?;

            let Some(sequence) = ctx.reference_sequence_repository.get(name).transpose()? else {
                continue;
            };

            reference_sequence = Some((id, sequence, usize::from(map.length())));
        }

        if let Some((_, sequence, len)) = &reference_sequence {
            record.remove_regenerated_md_nm(
                (**sequence).as_ref(),
                *len,
                store_md,
                store_nm,
                &mut md,
            );
        }
    }

    Ok(())
}

/// Writes a container of `records`, if any, encoding and compressing it on this thread.
pub fn write_container<W>(
    writer: &mut W,
    ctx: &Context,
    compressor: Option<&mut Compressor>,
    header: &sam::Header,
    record_counter: u64,
    records: &mut [Record],
) -> io::Result<()>
where
    W: Write,
{
    if records.is_empty() {
        return Ok(());
    }

    let container = encode_container(ctx, header, record_counter, records)?;

    let container = match compressor {
        Some(compressor) => {
            let mut choosing = Choosing::new(container);

            while let Some(jobs) = choosing.step(compressor) {
                let results = jobs.iter().map(TrialJob::run).collect();
                choosing.end_trials(compressor, results);
            }

            choosing.into_container()
        }
        None => container,
    };

    let buf = finish_container(container)?;
    writer.write_all(&buf)
}

/// A container with its records encoded, whose external blocks may not be compressed yet.
///
/// Writing a container takes three steps: [`encode_container`] and [`finish_container`] work on
/// one container alone, so they can run on any thread, but choosing the compression methods
/// ([`Choosing`]) must see the containers in order.
pub struct EncodedContainer {
    /// The header, without the landmarks and block count.
    header: Header,
    compression_header_block: Block,
    /// The slices with their header blocks.
    slices: Vec<(Block, Slice)>,
}

/// Encodes a container of `records` (not empty), the first of which is record `record_counter`
/// of the stream.
pub fn encode_container(
    ctx: &Context,
    header: &sam::Header,
    record_counter: u64,
    records: &mut [Record],
) -> io::Result<EncodedContainer> {
    let mut slices = Vec::new();
    let mut slice_record_counter = record_counter;

    // A container either requires the external reference or not. With records on sequences the
    // reference lacks, every record stores its aligned bases verbatim, as in HTSlib's
    // non-reference mode, and slices on one reference sequence embed a consensus of their records
    // to store the bases against, as HTSlib's `embed_ref=2` does.
    let mut embedded_reference_sequences = Vec::new();

    if records
        .iter()
        .any(|record| record.reference_bases == ReferenceBases::Verbatim)
    {
        for record in records.iter_mut() {
            record.store_bases_verbatim();
        }

        for chunk in records.chunks_mut(ctx.records_per_slice) {
            let reference_sequence_context = slice::get_reference_sequence_context(chunk);
            embedded_reference_sequences.push(embedded_reference::embed(
                reference_sequence_context,
                chunk,
            )?);
        }
    } else {
        remove_regenerated_md_nm(ctx, header, records)?;
    }

    let compression_header = build_compression_header(ctx, records);

    for (i, chunk) in records.chunks_mut(ctx.records_per_slice).enumerate() {
        let embedded_reference_sequence = embedded_reference_sequences
            .get(i)
            .and_then(|sequence| sequence.as_deref());

        let slice = build_slice(
            ctx,
            header,
            slice_record_counter,
            &compression_header,
            chunk,
            embedded_reference_sequence,
        )?;

        slices.push(slice);

        let record_count = u64::try_from(chunk.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        slice_record_counter += record_count;
    }

    let reference_sequence_context = get_container_reference_sequence_context(&slices)?;
    let record_count = records.len();
    let base_count = calculate_base_count(records)?;

    let mut buf = Vec::new();

    write_compression_header(&mut buf, &compression_header)?;
    let compression_header_block = build_compression_header_block(&buf)?;

    let slices = slices
        .into_iter()
        .map(|slice| {
            buf.clear();
            slice::write_header(&mut buf, &slice.header)?;
            Ok((build_slice_header_block(&buf)?, slice))
        })
        .collect::<io::Result<_>>()?;

    let header = Header {
        reference_sequence_context,
        record_count,
        record_counter,
        base_count,
        block_count: 0,
        landmarks: Vec::new(),
    };

    Ok(EncodedContainer {
        header,
        compression_header_block,
        slices,
    })
}

/// Compressing a block with one of a trial's methods (see [`choose_compression`]).
pub struct TrialJob {
    method: Method,
    src: Arc<Vec<u8>>,
    read_lengths: Option<Arc<[usize]>>,
}

impl TrialJob {
    /// Compresses the block, or returns `None` if the method cannot.
    pub fn run(&self) -> Option<Vec<u8>> {
        let quality_layout = QualityLayout {
            read_lengths: self.read_lengths.as_deref(),
        };

        self.method.encode(&self.src, &quality_layout)
    }
}

/// Choosing the compression method of each external block of a container, block by block, as the
/// compressor learns from them. The container must follow the last container chosen for.
///
/// Some blocks are trials, compressed with every candidate method; the compressor must learn from
/// them before the blocks after them are chosen for. A content ID appears once in a slice, so the
/// trials of a slice can run together.
pub struct Choosing {
    container: EncodedContainer,
    next_slice: usize,
    /// The trials of the last slice chosen for: the block's index, the trial, and the block's
    /// data.
    trials: Vec<(usize, Trial, Arc<Vec<u8>>)>,
}

impl Choosing {
    pub fn new(container: EncodedContainer) -> Self {
        Self {
            container,
            next_slice: 0,
            trials: Vec::new(),
        }
    }

    /// Chooses for the slices not chosen for yet, up to and including one with trials, and
    /// returns the jobs of those trials, whose results go to [`Self::end_trials`]. Returns `None`
    /// once every block is chosen for.
    pub fn step(&mut self, compressor: &mut Compressor) -> Option<Vec<TrialJob>> {
        assert!(self.trials.is_empty(), "the last trials have ended");

        while let Some((_, slice)) = self.container.slices.get_mut(self.next_slice) {
            self.next_slice += 1;

            let mut jobs = Vec::new();

            for (i, block) in slice.external_data_blocks.iter_mut().enumerate() {
                let ExternalBlock::Unchosen { content_id, src } = block else {
                    continue;
                };

                match compressor.choose(*content_id, src.len()) {
                    Choice::Method(method) => {
                        let src = mem::take(src);
                        *block = ExternalBlock::Chosen { method, src };
                    }
                    Choice::Trial(trial) => {
                        let src = Arc::new(mem::take(src));

                        jobs.extend(trial.methods().map(|method| TrialJob {
                            method,
                            src: Arc::clone(&src),
                            read_lengths: slice.read_lengths.clone(),
                        }));

                        self.trials.push((i, trial, src));
                    }
                }
            }

            if !self.trials.is_empty() {
                return Some(jobs);
            }
        }

        None
    }

    /// Learns from the results of the jobs [`Self::step`] returned, in order. Does nothing when no
    /// trials are pending.
    pub fn end_trials(&mut self, compressor: &mut Compressor, results: Vec<Option<Vec<u8>>>) {
        if self.trials.is_empty() {
            return;
        }

        let (_, slice) = &mut self.container.slices[self.next_slice - 1];
        let mut results = results.into_iter();

        for (i, trial, src) in self.trials.drain(..) {
            let n = trial.methods().count();
            let (block, sizes) = trial.pick(&src, results.by_ref().take(n).collect());
            compressor.end_trial(&trial, sizes);
            slice.external_data_blocks[i] = ExternalBlock::Compressed(block);
        }
    }

    /// Returns the container, once [`Self::step`] returned `None`.
    pub fn into_container(self) -> EncodedContainer {
        self.container
    }
}

/// Compresses the external blocks of `container` with their chosen methods and returns the
/// container as written.
pub fn finish_container(container: EncodedContainer) -> io::Result<Vec<u8>> {
    let EncodedContainer {
        mut header,
        compression_header_block,
        slices,
    } = container;

    let mut container_size = compression_header_block.size()?;
    let mut blocks = vec![compression_header_block];

    let slice_count = slices.len();

    // § 7 "Container header structure" (2024-09-04): "...landmarks[0] is equal to the byte length
    // of the compression header."
    let mut landmarks = vec![container_size];

    for (i, (slice_header_block, slice)) in slices.into_iter().enumerate() {
        let Slice {
            core_data_block,
            external_data_blocks,
            read_lengths,
            ..
        } = slice;

        let quality_layout = QualityLayout {
            read_lengths: read_lengths.as_deref(),
        };

        let mut slice_size = slice_header_block.size()?;
        blocks.push(slice_header_block);

        slice_size += core_data_block.size()?;
        blocks.push(core_data_block);

        for block in external_data_blocks {
            let block = match block {
                ExternalBlock::Compressed(block) => block,
                ExternalBlock::Chosen { method, src } => method.compress(&src, &quality_layout),
                ExternalBlock::Unchosen { .. } => {
                    unreachable!("compression is chosen before a container is finished")
                }
            };

            slice_size += block.size()?;
            blocks.push(block);
        }

        container_size += slice_size;

        if i < slice_count - 1 {
            landmarks.push(container_size);
        }
    }

    header.block_count = blocks.len();
    header.landmarks = landmarks;

    let mut buf = Vec::new();
    write_header(&mut buf, &header, container_size)?;

    for block in &blocks {
        write_block(&mut buf, block)?;
    }

    Ok(buf)
}

fn get_container_reference_sequence_context(
    slices: &[Slice],
) -> io::Result<ReferenceSequenceContext> {
    assert!(!slices.is_empty());

    let first_slice = slices.first().expect("slices cannot be empty");
    let mut container_reference_sequence_context = first_slice.header.reference_sequence_context();

    for slice in slices.iter().skip(1) {
        let slice_reference_sequence_context = slice.header.reference_sequence_context();

        match (
            container_reference_sequence_context,
            slice_reference_sequence_context,
        ) {
            (
                ReferenceSequenceContext::Some(container_context),
                ReferenceSequenceContext::Some(slice_context),
            ) if container_context.reference_sequence_id()
                == slice_context.reference_sequence_id() =>
            {
                let alignment_start = cmp::min(
                    container_context.alignment_start(),
                    slice_context.alignment_start(),
                );

                let alignment_end = cmp::max(
                    container_context.alignment_end(),
                    slice_context.alignment_end(),
                );

                container_reference_sequence_context = ReferenceSequenceContext::some(
                    container_context.reference_sequence_id(),
                    alignment_start,
                    alignment_end,
                );
            }
            (ReferenceSequenceContext::None, ReferenceSequenceContext::None) => {}
            (ReferenceSequenceContext::Many, ReferenceSequenceContext::Many) => {}
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "invalid slice reference sequence context: expected {container_reference_sequence_context:?}, got {slice_reference_sequence_context:?}"
                    ),
                ));
            }
        }
    }

    Ok(container_reference_sequence_context)
}

fn build_compression_header_block(src: &[u8]) -> io::Result<Block> {
    Block::encode(ContentType::CompressionHeader, 0, None, src)
}

fn calculate_base_count(records: &[Record]) -> io::Result<u64> {
    let n: usize = records.iter().map(|record| record.read_length).sum();
    u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn build_slice_header_block(src: &[u8]) -> io::Result<Block> {
    Block::encode(ContentType::SliceHeader, 0, None, src)
}

// § 9 "End of file container" (2022-04-12)
pub static EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

pub fn write_eof_container<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&EOF)
}
