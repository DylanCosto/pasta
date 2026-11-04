//! `ntools index`: build BAI, CSI, or CRAI indexes (samtools index).

use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use pasta::bam::{self as bam, bai, index::Indexer};
use pasta::bgzf;
use pasta::core::Position;
use pasta::cram::{self as cram, crai};
use pasta::csi::{
    self as csi,
    binning_index::index::reference_sequence::{bin::Chunk, index::BinnedIndex},
};
use pasta::sam::{self as sam, alignment::Record as _};

use crate::{
    bgzf_stream::BgzfReader,
    source::{Format, Input, Stream},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Coordinate-sorted BAM or CRAM files to index (several with `-M`); without `-M`, a second
    /// argument that does not exist or is an index names the index to write.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Write the index to this path instead of next to the input.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Build a BAI index (the default for BAM).
    #[arg(short = 'b', long = "bai", overrides_with_all = ["csi", "min_shift"])]
    bai: bool,

    /// Build a CSI index.
    #[arg(short = 'c', long = "csi", overrides_with = "bai")]
    csi: bool,

    /// Build a CSI index with intervals of 2^INT bases at its finest level [14].
    #[arg(
        short = 'm',
        long = "min-shift",
        overrides_with = "bai",
        allow_negative_numbers = true
    )]
    min_shift: Option<i32>,

    /// Index every file argument.
    #[arg(short = 'M')]
    multiple: bool,

    /// Additional threads for BGZF decompression.
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,
}

/// HTSlib's default CSI minimum shift (`BAM_LIDX_SHIFT`).
const DEFAULT_MIN_SHIFT: i32 = 14;

pub fn run(args: Args) -> anyhow::Result<()> {
    let csi = args.csi || args.min_shift.is_some();
    let min_shift = args.min_shift.unwrap_or(DEFAULT_MIN_SHIFT);

    // As in HTSlib, a minimum shift of 0 or less means BAI.
    let kind = if csi && min_shift > 0 {
        IndexKind::CsiWithMinShift(u8::try_from(min_shift).context("invalid min_shift")?)
    } else {
        IndexKind::Default
    };

    let mut inputs = args.inputs;
    let mut output = args.output;

    // The legacy synopsis: `index in.bam out.index`.
    if inputs.len() == 2 && output.is_none() && is_missing_or_index(&inputs[1]) {
        output = inputs.pop();
    }

    if inputs.len() > 1 && !args.multiple {
        bail!("use -M to enable indexing more than one alignment file");
    }

    if output.is_some() && inputs.len() > 1 {
        bail!("can't use -o with multiple input alignment files");
    }

    for input in &inputs {
        index_file(input, output.clone(), kind, args.threads)?;
    }

    Ok(())
}

/// Whether a path does not exist or holds an index (BAI, CSI, TBI, or CRAI), as HTSlib's
/// `nonexistent_or_index` decides whether a second argument names the index to write.
fn is_missing_or_index(path: &Path) -> bool {
    use std::io::Read as _;

    let Ok(mut file) = File::open(path) else {
        return true;
    };

    let mut buf = Vec::new();

    if (&mut file).take(4).read_to_end(&mut buf).is_err() {
        return false;
    }

    if buf == b"BAI\x01" {
        return true;
    }

    if !buf.starts_with(&[0x1f, 0x8b]) {
        return false;
    }

    // Compressed: CSI and TBI have magic numbers; a CRAI is lines of six integers.
    let Ok(file) = File::open(path) else {
        return false;
    };

    let mut text = Vec::new();
    let mut decoder = flate2::read::MultiGzDecoder::new(file).take(256);

    if decoder.read_to_end(&mut text).is_err() && text.is_empty() {
        return false;
    }

    text.starts_with(b"CSI\x01") || text.starts_with(b"TBI\x01") || is_crai_text(&text)
}

/// Whether text starts like a CRAI index: tab-separated integers on each line.
fn is_crai_text(text: &[u8]) -> bool {
    let Some(line) = text.split(|&b| b == b'\n').next() else {
        return false;
    };

    let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();

    fields.len() == 6
        && fields.iter().all(|field| {
            let digits = field.strip_prefix(b"-").unwrap_or(field);
            !digits.is_empty() && digits.iter().all(u8::is_ascii_digit)
        })
}

/// Which index to build for BAM; CRAM always gets CRAI.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IndexKind {
    /// BAI for BAM.
    Default,
    /// CSI with HTSlib's default minimum shift.
    Csi,
    /// CSI with intervals of 2^n bases at its finest level.
    CsiWithMinShift(u8),
}

/// Indexes a BAM or CRAM file, writing next to it (`in.bam.bai`) unless `output` is given.
pub fn index_file(
    input: &Path,
    output: Option<PathBuf>,
    kind: IndexKind,
    threads: usize,
) -> anyhow::Result<()> {
    let src = input.to_str().context("input path is not valid UTF-8")?;
    let opened = Input::open(src)?;

    match opened.format {
        Format::Bam => {
            let (min_shift, extension) = match kind {
                IndexKind::Default => (None, "bai"),
                IndexKind::Csi => (Some(DEFAULT_MIN_SHIFT as u8), "csi"),
                IndexKind::CsiWithMinShift(min_shift) => (Some(min_shift), "csi"),
            };

            let dst = output.unwrap_or_else(|| with_appended_extension(input, extension));

            let index = index_bam(opened.reader, min_shift, threads)
                .with_context(|| format!("indexing {}", input.display()))?;

            write_bam_index(&dst, &index)
        }
        // As in samtools, CSI options do not apply to CRAM.
        Format::Cram => {
            let dst = output.unwrap_or_else(|| with_appended_extension(input, "crai"));

            let index =
                cram::fs::index(input).with_context(|| format!("indexing {}", input.display()))?;

            let file = File::create(&dst).with_context(|| format!("creating {}", dst.display()))?;
            let mut writer = crai::io::Writer::new(file);
            writer.write_index(&index)?;
            writer.finish()?;

            Ok(())
        }
        Format::Sam => bail!("only BAM and CRAM files can be indexed"),
    }
}

/// HTSlib's CSI layout for a longest reference sequence (`hts_adjust_csi_settings`): enough
/// levels above `min_shift` to cover it, or, past 9 levels, a larger minimum shift. Returns the
/// minimum shift and depth.
fn csi_layout(max_length: usize, min_shift: u8) -> (u8, u8) {
    const MAX_DEPTH: u8 = 9;

    let max_length = u128::try_from(max_length).unwrap_or(u128::MAX) + 256;
    let max_position = |min_shift: u8, depth: u8| {
        1u128
            .checked_shl(u32::from(min_shift) + 3 * u32::from(depth))
            .unwrap_or(u128::MAX)
    };

    if max_length <= max_position(min_shift, MAX_DEPTH) {
        let mut depth = 0;

        while max_length > max_position(min_shift, depth) {
            depth += 1;
        }

        (min_shift, depth)
    } else {
        let mut shift = min_shift;

        while max_length > max_position(shift, MAX_DEPTH) {
            shift += 1;
        }

        eprintln!(
            "[W::hts_adjust_csi_settings] Adjusted min_shift from {min_shift} to {shift} due to \
             longest reference of {} bases.",
            max_length - 256
        );

        (shift, MAX_DEPTH)
    }
}

/// A BAI or CSI indexer.
enum BamIndexer {
    Bai(Indexer),
    Csi(csi::binning_index::Indexer<BinnedIndex>),
}

impl BamIndexer {
    fn add_record(
        &mut self,
        context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> anyhow::Result<()> {
        match self {
            Self::Bai(indexer) => indexer.add_record(context, chunk)?,
            Self::Csi(indexer) => indexer.add_record(context, chunk)?,
        }

        Ok(())
    }

    fn build(self, reference_sequence_count: usize) -> bam::Index {
        match self {
            Self::Bai(indexer) => indexer.build(reference_sequence_count),
            Self::Csi(indexer) => bam::Index::Csi(indexer.build(reference_sequence_count)),
        }
    }
}

/// Builds a BAI, or with a minimum shift, a CSI.
fn index_bam(
    reader: BufReader<Stream>,
    min_shift: Option<u8>,
    threads: usize,
) -> anyhow::Result<bam::Index> {
    let mut reader = bam::io::Reader::from(BgzfReader::new(reader, threads));
    let header = reader.read_header()?;

    let max_reference_sequence_length = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length())
        .max()
        .unwrap_or(NonZero::<usize>::MIN);

    let max_position = Position::new(max_reference_sequence_length.get())
        .context("reference sequence length overflows a position")?;

    let mut indexer = match min_shift {
        Some(min_shift) => {
            let (min_shift, depth) = csi_layout(max_reference_sequence_length.get(), min_shift);
            let indexer = csi::binning_index::Indexer::new(min_shift, depth)
                .with_context(|| format!("invalid CSI layout: min_shift {min_shift}"))?;
            BamIndexer::Csi(indexer)
        }
        None => BamIndexer::Bai(
            Indexer::builder()
                .set_format(bam::index::Format::Bai)
                .set_max_position_hint(max_position)
                .build()?,
        ),
    };

    let mut record = bam::Record::default();
    let mut start = bgzf::io::Read::virtual_position(reader.get_ref());

    // As in samtools, the records' order is checked, not the header's sort order.
    let mut order = RecordOrder::new(header.reference_sequences().len());

    while reader.read_record(&mut record)? != 0 {
        let end = bgzf::io::Read::virtual_position(reader.get_ref());

        order
            .push(&record)
            .map_err(|message| unindexable(&header, &record, &message))?;

        let context = alignment_context(&record)?
            .map(|(id, start, end)| (id, start, end, !record.flags().is_unmapped()));

        indexer
            .add_record(context, Chunk::new(start, end))
            .map_err(|e| unindexable(&header, &record, &e.to_string()))?;

        start = end;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

/// The record order checks of HTSlib's `hts_idx_push`: each reference sequence's records form
/// one block, sorted by position, and unplaced records come last.
struct RecordOrder {
    /// The reference sequence of the previous record (`None` for unplaced), if any.
    previous: Option<Option<usize>>,
    /// The 0-based position of the previous record (-1 without one).
    previous_position: i64,
    /// Whether each reference sequence's block has started.
    started: Vec<bool>,
    unplaced: bool,
}

impl RecordOrder {
    fn new(reference_sequence_count: usize) -> Self {
        Self {
            previous: None,
            previous_position: -1,
            started: vec![false; reference_sequence_count],
            unplaced: false,
        }
    }

    /// Checks that a record may follow the previous ones, with HTSlib's message if not.
    fn push(&mut self, record: &bam::Record) -> Result<(), String> {
        let id = record
            .reference_sequence_id()
            .transpose()
            .map_err(|e| e.to_string())?;

        let position = match (id, record.alignment_start().transpose()) {
            (None, _) => -1,
            (Some(_), Ok(start)) => start.map_or(-1, |start| usize::from(start) as i64 - 1),
            (Some(_), Err(e)) => return Err(e.to_string()),
        };

        if self.previous != Some(id) {
            if let Some(id) = id {
                if self.unplaced {
                    return Err(format!(
                        "NO_COOR reads not in a single block at the end {id} {}",
                        self.previous.flatten().map_or(-1, |id| id as i64)
                    ));
                }

                match self.started.get_mut(id) {
                    Some(true) => return Err("Chromosome blocks not continuous".into()),
                    Some(started) => *started = true,
                    None => return Err(format!("invalid reference sequence ID {id}")),
                }
            }

            self.previous = Some(id);
        } else if let Some(id) = id
            && self.previous_position > position
        {
            return Err(format!(
                "Unsorted positions on sequence #{}: {} followed by {}",
                id + 1,
                self.previous_position + 1,
                position + 1
            ));
        }

        if id.is_none() {
            self.unplaced = true;
        }

        self.previous_position = position;

        Ok(())
    }
}

/// The error for a record out of order, with the details HTSlib's `sam_index` gives.
fn unindexable(header: &sam::Header, record: &bam::Record, message: &str) -> anyhow::Error {
    let reference_sequence = record
        .reference_sequence_id()
        .and_then(Result::ok)
        .and_then(|id| header.reference_sequences().get_index(id));

    let (name, length) = reference_sequence.map_or_else(
        || ("(null)".to_owned(), 0),
        |(name, rs)| (name.to_string(), rs.length().get()),
    );

    let position = record
        .alignment_start()
        .and_then(Result::ok)
        .map_or(0, usize::from);

    anyhow::anyhow!(
        "{message}: read '{}' with ref_name='{name}', ref_length={length}, flags={}, pos={position} \
         cannot be indexed",
        record.name().map(ToString::to_string).unwrap_or_default(),
        record.flags().bits(),
    )
}

fn write_bam_index(dst: &Path, index: &bam::Index) -> anyhow::Result<()> {
    let file = File::create(dst).with_context(|| format!("creating {}", dst.display()))?;

    match index {
        bam::Index::Bai(index) => {
            let mut writer = bai::io::Writer::new(BufWriter::new(file));
            writer.write_index(index)?;
            writer.get_mut().flush()?;
        }
        bam::Index::Csi(index) => {
            let mut writer = csi::io::Writer::new(file);
            writer.write_index(index)?;
            writer.get_mut().try_finish()?;
        }
    }

    Ok(())
}

/// The reference sequence, start, and end a record is indexed under.
fn alignment_context(record: &bam::Record) -> anyhow::Result<Option<(usize, Position, Position)>> {
    let (Some(id), Some(start)) = (
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };

    // As HTSlib's `bam_endpos`: an unmapped record ends where it starts, whatever its CIGAR.
    let end = if record.flags().is_unmapped() {
        Some(start)
    } else {
        record.alignment_end().transpose()?
    };

    Ok(end.map(|end| (id, start, end)))
}

/// Appends an extension, as samtools names indexes: `in.bam` becomes `in.bam.bai`.
fn with_appended_extension(path: &Path, extension: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(extension);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use bstr::BString;
    use pasta::sam::{
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags,
                cigar::{Op, op::Kind},
            },
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    use super::*;

    #[test]
    fn test_alignment_context_of_an_unmapped_record_with_a_cigar() -> anyhow::Result<()> {
        let header = sam::Header::builder()
            .set_reference_sequences(
                [(
                    BString::from("sq0"),
                    Map::<ReferenceSequence>::new(const { NonZero::new(100_000).unwrap() }),
                )]
                .into_iter()
                .collect(),
            )
            .build();

        let start = Position::try_from(16_380)?;
        let mut data = Vec::new();

        {
            let mut writer = bam::io::Writer::from(&mut data);
            writer.write_header(&header)?;

            for flags in [Flags::empty(), Flags::UNMAPPED] {
                let record = RecordBuf::builder()
                    .set_flags(flags)
                    .set_reference_sequence_id(0)
                    .set_alignment_start(start)
                    .set_cigar([Op::new(Kind::Match, 50)].into_iter().collect())
                    .set_sequence(vec![b'A'; 50].into())
                    .build();

                writer.write_alignment_record(&header, &record)?;
            }
        }

        let mut reader = bam::io::Reader::from(&data[..]);
        reader.read_header()?;
        let mut record = bam::Record::default();

        // An unmapped record is indexed under its start alone, as samtools does.
        reader.read_record(&mut record)?;
        assert_eq!(
            alignment_context(&record)?,
            Some((0, start, Position::try_from(16_429)?))
        );

        reader.read_record(&mut record)?;
        assert_eq!(alignment_context(&record)?, Some((0, start, start)));

        Ok(())
    }
}
