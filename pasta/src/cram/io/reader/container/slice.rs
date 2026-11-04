mod header;
pub mod records;

use std::{borrow::Cow, io, sync::Arc};

use crate::core::Position;
use crate::fasta;
use crate::sam::{self as sam, alignment::Record as _};

pub(crate) use self::records::MdNmTags;
use self::{
    header::read_header,
    records::{ExternalDataReaders, Records, Selection},
};

use super::read_block_as;
use crate::cram::{
    Record, calculate_normalized_sequence_digest,
    container::{
        CompressionHeader, ReferenceSequenceContext,
        block::{self, ContentType},
        slice::Header,
    },
    io::{
        BitReader,
        reader::{Fields, MdNm, Tags},
    },
    record::{Feature, Flags},
};

/// A container slice.
///
/// A slice contains a header, a core data block, and one or more external blocks. This is where
/// the CRAM records are stored.
pub struct Slice<'c> {
    header: Header,
    src: &'c [u8],
}

impl<'c> Slice<'c> {
    pub(crate) fn header(&self) -> &Header {
        &self.header
    }

    #[allow(clippy::type_complexity)]
    pub fn decode_blocks(
        &self,
    ) -> io::Result<(Cow<'c, [u8]>, Vec<(block::ContentId, Cow<'c, [u8]>)>)> {
        self.decode_selected_blocks(|_| true)
    }

    /// Decodes the core data block and the external blocks that decoding `fields` and `tags`
    /// reads.
    ///
    /// The other external blocks are left out. Pass the result to [`Self::records_with_fields`]
    /// with the same `fields` and `tags`.
    ///
    /// # Errors
    ///
    /// Returns an error if a block cannot be read or decompressed.
    #[allow(clippy::type_complexity)]
    pub fn decode_blocks_for(
        &self,
        compression_header: &CompressionHeader,
        fields: Fields,
        tags: &Tags,
    ) -> io::Result<(Cow<'c, [u8]>, Vec<(block::ContentId, Cow<'c, [u8]>)>)> {
        let selection = Selection::new(compression_header, fields, tags, MdNm::AsStored);
        let embedded_reference_id = self
            .header
            .embedded_reference_bases_block_content_id()
            .filter(|_| selection.reads_reference_sequences());

        self.decode_selected_blocks(|id| {
            selection.reads_block(id) || embedded_reference_id == Some(id)
        })
    }

    #[allow(clippy::type_complexity)]
    fn decode_selected_blocks<F>(
        &self,
        mut is_selected: F,
    ) -> io::Result<(Cow<'c, [u8]>, Vec<(block::ContentId, Cow<'c, [u8]>)>)>
    where
        F: FnMut(block::ContentId) -> bool,
    {
        let mut src = self.src;

        let block = read_block_as(&mut src, ContentType::CoreData)?;
        let core_data_src = block.decode()?;

        let external_data_block_count = self.header.block_count() - 1;
        let mut external_data_srcs = Vec::with_capacity(external_data_block_count);

        for _ in 0..external_data_block_count {
            let block = read_block_as(&mut src, ContentType::ExternalData)?;

            if is_selected(block.content_id) {
                external_data_srcs.push((block.content_id, block.decode()?));
            }
        }

        Ok((core_data_src, external_data_srcs))
    }

    /// Reads and returns a list of raw records in this slice.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use pasta::cram::{self as cram, io::reader::Container};
    /// use pasta::fasta;
    ///
    /// let data = [];
    /// let mut reader = cram::io::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     let compression_header = container.compression_header()?;
    ///
    ///     for result in container.slices() {
    ///         let slice = result?;
    ///
    ///         let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
    ///
    ///         let records = slice.records(
    ///             fasta::Repository::default(),
    ///             &header,
    ///             &compression_header,
    ///             &core_data_src,
    ///             &external_data_srcs,
    ///         )?;
    ///
    ///         // ...
    ///     }
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        self.records_with_md_nm(
            reference_sequence_repository,
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            MdNm::Regenerate,
        )
    }

    /// Reads and returns a list of raw records in this slice, adding the MD and NM tags the
    /// slice left out only if `md_nm` is [`MdNm::Regenerate`].
    ///
    /// [`Self::records`] always regenerates them.
    pub fn records_with_md_nm<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
        md_nm: MdNm,
    ) -> io::Result<Vec<Record<'c>>> {
        self.records_with_fields(
            reference_sequence_repository,
            header,
            compression_header,
            core_data_src,
            external_data_srcs,
            md_nm,
            Fields::all(),
            &Tags::All,
        )
    }

    /// Reads and returns a list of records in this slice, decoding only `fields` (and the fields
    /// that share their blocks) and, with [`Fields::DATA`], only `tags`; see [`Fields`] and
    /// [`Tags`].
    ///
    /// `core_data_src` and `external_data_srcs` are from [`Self::decode_blocks_for`] with the
    /// same `fields` and `tags`, or from [`Self::decode_blocks`].
    ///
    /// # Errors
    ///
    /// Returns an error if a record cannot be decoded or a reference sequence cannot be read.
    #[expect(
        clippy::too_many_arguments,
        reason = "`records_with_md_nm` with the fields to decode"
    )]
    pub fn records_with_fields<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
        md_nm: MdNm,
        fields: Fields,
        tags: &Tags,
    ) -> io::Result<Vec<Record<'c>>> {
        let selection = Selection::new(compression_header, fields, tags, md_nm);
        let md_nm = selection.md_nm();

        let core_data_reader = BitReader::new(core_data_src);

        let mut external_data_readers = ExternalDataReaders::new();

        for (block_content_id, src) in external_data_srcs {
            external_data_readers.insert(*block_content_id, src);
        }

        let reference_sequence_context = self.header.reference_sequence_context();
        let initial_id = self.header.record_counter();

        let mut reader = Records::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
            initial_id,
            &selection,
        )?;

        let reads_reference_sequences = selection.reads_reference_sequences();

        let slice_reference_sequence = if reads_reference_sequences {
            get_slice_reference_sequence(
                &reference_sequence_repository.clone(),
                header,
                compression_header,
                &self.header,
                external_data_srcs,
                true,
            )?
        } else {
            None
        };

        // Loaded for records that read outside a partial slice reference sequence.
        let mut whole_slice_reference_sequence = None;

        let substitution_matrix = compression_header.preservation_map().substitution_matrix();

        let mut records: Vec<Record<'c>> = vec![Record::default(); self.header.record_count()];

        for record in &mut records {
            reader.read_record(record)?;

            record.header = Some(header);

            if !reads_reference_sequences {
                // The sequence is not decoded, so it reads as missing rather than being rebuilt
                // from read features without a reference sequence.
                record.cram_flags.insert(Flags::SEQUENCE_IS_MISSING);
            } else if needs_reference_sequence(record) {
                // A slice with many reference sequences has no slice reference sequence. Its
                // records get theirs once all are read, so that each reference sequence can be
                // read once for the span of its records.
                if let Some(reference_sequence) = &slice_reference_sequence {
                    record.reference_sequence =
                        select_reference_sequence(reference_sequence, record, || {
                            if whole_slice_reference_sequence.is_none() {
                                whole_slice_reference_sequence =
                                    Some(get_slice_reference_sequence(
                                        &reference_sequence_repository,
                                        header,
                                        compression_header,
                                        &self.header,
                                        external_data_srcs,
                                        false,
                                    )?);
                            }

                            Ok(whole_slice_reference_sequence.clone().flatten())
                        })?;
                }

                record.substitution_matrix = substitution_matrix.clone();
            }
        }

        // Without a required external reference, records read their bases from features (and the
        // embedded reference, if any).
        if reads_reference_sequences
            && reference_sequence_context.is_many()
            && compression_header
                .preservation_map()
                .external_reference_sequence_is_required()
        {
            set_record_reference_sequences(&reference_sequence_repository, header, &mut records)?;
        }

        for record in &mut records {
            crate::cram::record::md_nm::add_missing_md_nm(record, md_nm);
        }

        let names = if compression_header.preservation_map().records_have_names() {
            Names::Stored
        } else {
            Names::Generated
        };

        resolve_mates(&mut records, names)?;

        Ok(records)
    }
}

impl<'c> Slice<'c> {
    /// Decodes the records without reference sequences, mates, or `MD`/`NM` tags.
    ///
    /// Their positions, read lengths, and features are set, which is enough to compute
    /// alignment spans (e.g., for indexing) without a reference sequence repository.
    pub(crate) fn records_without_reference_sequences<'ch: 'c>(
        &self,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        let core_data_reader = BitReader::new(core_data_src);

        let mut external_data_readers = ExternalDataReaders::new();

        for (block_content_id, src) in external_data_srcs {
            external_data_readers.insert(*block_content_id, src);
        }

        let mut reader = Records::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            self.header.reference_sequence_context(),
            self.header.record_counter(),
            &Selection::all(),
        )?;

        let mut records: Vec<Record<'c>> = vec![Record::default(); self.header.record_count()];

        for record in &mut records {
            reader.read_record(record)?;
        }

        Ok(records)
    }
}

pub fn read_slice<'c>(src: &mut &'c [u8]) -> io::Result<Slice<'c>> {
    let header = read_header(src)?;
    Ok(Slice { header, src })
}

/// Where read names come from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Names {
    /// The slice stores them; a missing one (`*`) stays missing.
    Stored,
    /// The slice has none (lossy read names): each record is named by its record counter, and
    /// mates take the name of the first record of their template.
    Generated,
}

fn resolve_mates(records: &mut [Record], names: Names) -> io::Result<()> {
    let mut mate_indices: Vec<_> = records
        .iter()
        .enumerate()
        .map(|(i, record)| record.mate_distance.map(|len| i + len + 1))
        .collect();

    if mate_indices.iter().flatten().any(|&j| j >= records.len()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "mate record is outside the slice",
        ));
    }

    for i in 0..records.len() {
        let record = &mut records[i];

        if names == Names::Generated && record.name.is_none() {
            let name = record.id.to_string().into_bytes();
            record.name = Some(Cow::from(name));
        }

        if mate_indices[i].is_none() {
            continue;
        }

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let mid = j + 1;
            let (left, right) = records.split_at_mut(mid);

            let record = &mut left[j];
            let mate = &mut right[mate_index - mid];
            set_mate(record, mate);

            if names == Names::Generated && mate.name.is_none() {
                mate.name = record.name.clone();
            }

            j = mate_index;
        }

        let (left, right) = records.split_at_mut(j);
        let record = &mut right[0];
        let mate = &mut left[i];
        set_mate(record, mate);

        set_template_lengths(records, &mut mate_indices, i);
    }

    Ok(())
}

fn set_mate(record: &mut Record, mate: &mut Record) {
    set_mate_chunk(
        &mut record.bam_flags,
        &mut record.mate_reference_sequence_id,
        &mut record.mate_alignment_start,
        mate.bam_flags,
        mate.reference_sequence_id,
        mate.alignment_start,
    );
}

fn set_mate_chunk(
    record_bam_flags: &mut sam::alignment::record::Flags,
    record_mate_reference_sequence_id: &mut Option<usize>,
    record_mate_alignment_start: &mut Option<Position>,
    mate_bam_flags: sam::alignment::record::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
) {
    if mate_bam_flags.is_reverse_complemented() {
        *record_bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_flags.is_unmapped() {
        *record_bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
    }

    *record_mate_reference_sequence_id = mate_reference_sequence_id;
    *record_mate_alignment_start = mate_alignment_start;
}

/// Returns the alignment start and end of a record, as HTSlib uses them for template lengths:
/// 0 for an unplaced record, and the start for an unmapped one or one that covers no reference
/// bases.
fn template_extent(record: &Record) -> (usize, usize) {
    use crate::cram::record::calculate_alignment_span;

    let start = record.alignment_start.map(usize::from).unwrap_or_default();

    if record.bam_flags.is_unmapped() {
        return (start, start);
    }

    let span = calculate_alignment_span(record.read_length, &record.features);
    let end = (start + span).saturating_sub(1).max(start);

    (start, end)
}

/// Sets the template lengths of the attached template that starts at record `i`, as HTSlib does
/// (`cram_decode_slice_xref`), and detaches its records in `mate_indices`.
///
/// The length spans all segments when they are on the same reference sequence, and is 0
/// otherwise. The first record's length is positive when it is leftmost (and not also rightmost
/// with another segment), with ties of identical extents broken by the first-segment flag, and the
/// other records get the opposite sign. Unmapped records and records with unmapped mates get 0.
fn set_template_lengths(records: &mut [Record], mate_indices: &mut [Option<usize>], i: usize) {
    let (first_start, first_end) = template_extent(&records[i]);
    let reference_sequence_id = records[i].reference_sequence_id;

    let (mut left, mut right) = (first_start, first_end);
    let (mut left_count, mut right_count) = (0, 0);
    let mut same_reference_sequence = reference_sequence_id.is_some();

    let mut j = Some(i);

    while let Some(k) = j {
        let (start, end) = template_extent(&records[k]);

        if start < left {
            (left, left_count) = (start, 1);
        } else if start == left {
            left_count += 1;
        }

        if end > right {
            (right, right_count) = (end, 1);
        } else if end == right {
            right_count += 1;
        }

        same_reference_sequence &= records[k].reference_sequence_id == reference_sequence_id;

        j = mate_indices[k];
    }

    let len = if same_reference_sequence {
        i32::try_from(right - left + 1).unwrap_or_default()
    } else {
        0
    };

    let first_is_leftmost = first_start == left && (first_end < right || left_count <= 1);
    let first_is_tied =
        first_start == left && first_end == right && left_count > 1 && right_count > 1;

    let first_len =
        if first_is_leftmost || (first_is_tied && records[i].bam_flags.is_first_segment()) {
            len
        } else {
            -len
        };

    let mut prev = i;
    let mut j = mate_indices[i].take();
    records[i].template_length = first_len;

    while let Some(k) = j {
        records[k].template_length = -first_len;
        clear_unmapped_template_length(records, prev, k);
        prev = k;
        j = mate_indices[k].take();
    }

    clear_unmapped_template_length(records, prev, i);
}

/// Sets the template length of record `i` to 0 when it or its mate, record `mate`, is unmapped.
fn clear_unmapped_template_length(records: &mut [Record], i: usize, mate: usize) {
    if records[i].bam_flags.is_unmapped() || records[mate].bam_flags.is_unmapped() {
        records[i].template_length = 0;
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum ReferenceSequence<'c> {
    Embedded {
        reference_start: Position,
        sequence: &'c [u8],
    },
    External {
        sequence: Arc<fasta::record::Sequence>,
    },
    /// The part of an external reference sequence that starts at `reference_start`.
    ExternalRange {
        reference_start: Position,
        sequence: Arc<fasta::record::Sequence>,
    },
}

fn get_slice_reference_sequence<'c>(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    compression_header: &CompressionHeader,
    slice_header: &Header,
    external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    allow_range: bool,
) -> io::Result<Option<ReferenceSequence<'c>>> {
    let reference_sequence_context = slice_header.reference_sequence_context();

    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok(None);
    };

    let external_reference_sequence_is_required = compression_header
        .preservation_map()
        .external_reference_sequence_is_required();

    let embedded_reference_bases_block_content_id =
        slice_header.embedded_reference_bases_block_content_id();

    if external_reference_sequence_is_required {
        let (reference_sequence_name, reference_sequence) = header
            .reference_sequences()
            .get_index(context.reference_sequence_id())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid slice reference sequence ID",
                )
            })?;

        // § 8.5 "Slice header block" (2024-09-04): "MD5sums should not be validated if the stored
        // checksum is all-zero." The reference sequence is then checked against the header.
        if slice_header.reference_md5().is_none() {
            validate_reference_sequence(
                reference_sequence_repository,
                reference_sequence_name,
                reference_sequence,
            )?;
        }

        if allow_range
            && let Some(reference_sequence) = get_slice_reference_sequence_range(
                reference_sequence_repository,
                reference_sequence_name,
                reference_sequence.length().get(),
                context.alignment_start(),
                context.alignment_end(),
                slice_header.reference_md5(),
            )
        {
            return Ok(Some(reference_sequence));
        }

        let sequence = reference_sequence_repository
            .get(reference_sequence_name)
            .transpose()?
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("missing reference sequence: {reference_sequence_name}"),
                )
            })?;

        // § 8.5 "Slice header block" (2024-09-04): "MD5sums should not be validated if the stored
        // checksum is all-zero."
        if let Some(expected_md5) = slice_header.reference_md5() {
            let interval = context.alignment_start()..=context.alignment_end();
            let subsequence = sequence.get(interval).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "slice alignment span is outside the reference sequence",
                )
            })?;
            validate_sequence(
                reference_sequence_name,
                context.alignment_start(),
                context.alignment_end(),
                subsequence,
                expected_md5,
            )?;
        }

        Ok(Some(ReferenceSequence::External { sequence }))
    } else if let Some(block_content_id) = embedded_reference_bases_block_content_id {
        let sequence = external_data_srcs
            .iter()
            .find(|(id, _)| *id == block_content_id)
            .map(|(_, src)| src)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "invalid block content ID")
            })?;

        Ok(Some(ReferenceSequence::Embedded {
            reference_start: context.alignment_start(),
            sequence,
        }))
    } else {
        Ok(None)
    }
}

/// Reads only the slice's span of a reference sequence, as HTSlib does when the span is less than
/// half of the reference sequence.
///
/// This returns `None` to use the whole reference sequence instead: when the span is too long, the
/// repository cannot read ranges, the read fails or is short, or the checksum does not match. The
/// whole sequence then gives the same result (or error) as without ranges.
fn get_slice_reference_sequence_range<'c>(
    reference_sequence_repository: &fasta::Repository,
    reference_sequence_name: &[u8],
    reference_sequence_length: usize,
    alignment_start: Position,
    alignment_end: Position,
    reference_md5: Option<&[u8; 16]>,
) -> Option<ReferenceSequence<'c>> {
    let len = usize::from(alignment_end)
        .checked_sub(usize::from(alignment_start))?
        .checked_add(1)?;

    if (len - 1).saturating_mul(2) >= reference_sequence_length {
        return None;
    }

    let sequence = reference_sequence_repository
        .get_range(
            reference_sequence_name,
            (alignment_start..=alignment_end).into(),
        )?
        .ok()?;

    if sequence.len() != len {
        return None;
    }

    if let Some(expected_md5) = reference_md5 {
        validate_sequence(
            reference_sequence_name,
            alignment_start,
            alignment_end,
            (*sequence).as_ref(),
            expected_md5,
        )
        .ok()?;
    }

    Some(ReferenceSequence::ExternalRange {
        reference_start: alignment_start,
        sequence,
    })
}

fn needs_reference_sequence(record: &Record<'_>) -> bool {
    !record.bam_flags.is_unmapped() && !record.cram_flags.sequence_is_missing()
}

/// Sets the reference sequences of the records of a slice with many reference sequences.
///
/// Each reference sequence is read once for the span of its records in the slice, as for a
/// single-reference slice. Records outside that span, and reference sequences that cannot be read
/// in part, use the whole reference sequence.
fn set_record_reference_sequences(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    records: &mut [Record<'_>],
) -> io::Result<()> {
    use std::collections::{HashMap, HashSet};

    // The (start, end) span of the records that need a reference sequence, by reference
    // sequence ID. Maps keep this linear in the number of reference sequences, which can be
    // as many as the records (e.g., a transcriptome).
    let mut spans: HashMap<usize, (Position, Position)> = HashMap::new();
    // The reference sequences that records need, in order, checked once each.
    let mut ids: Vec<usize> = Vec::new();
    let mut seen: HashSet<usize> = HashSet::new();

    for record in records
        .iter()
        .filter(|record| needs_reference_sequence(record))
    {
        let Some(id) = record.reference_sequence_id else {
            continue;
        };

        if seen.insert(id) {
            ids.push(id);
        }

        let Some((start, len)) = reference_span(record) else {
            continue;
        };

        let Some(end) = start.checked_add(len.max(1) - 1) else {
            continue;
        };

        spans
            .entry(id)
            .and_modify(|(span_start, span_end)| {
                *span_start = (*span_start).min(start);
                *span_end = (*span_end).max(end);
            })
            .or_insert((start, end));
    }

    // Slices with many reference sequences have no reference checksum, so the reference sequences
    // are checked against the header instead.
    for id in ids {
        let (name, reference_sequence) =
            header.reference_sequences().get_index(id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid record reference sequence ID",
                )
            })?;

        validate_reference_sequence(reference_sequence_repository, name, reference_sequence)?;
    }

    let ranges: HashMap<_, _> = spans
        .into_iter()
        .filter_map(|(id, (start, end))| {
            let (name, reference_sequence) = header.reference_sequences().get_index(id)?;

            // Slices with many reference sequences have no reference checksum.
            get_slice_reference_sequence_range(
                reference_sequence_repository,
                name,
                reference_sequence.length().get(),
                start,
                end,
                None,
            )
            .map(|reference_sequence| (id, reference_sequence))
        })
        .collect();

    for record in records
        .iter_mut()
        .filter(|record| needs_reference_sequence(record))
    {
        let range = record.reference_sequence_id.and_then(|id| ranges.get(&id));

        let get_whole =
            || get_record_reference_sequence(reference_sequence_repository, header, record);

        record.reference_sequence = match range {
            Some(range) => select_reference_sequence(range, record, get_whole)?,
            None => get_whole()?,
        };
    }

    Ok(())
}

/// Returns the given reference sequence for the record, or the whole reference sequence from
/// `get_whole` if the given one is a part that the record reads outside of.
fn select_reference_sequence<'c, F>(
    reference_sequence: &ReferenceSequence<'c>,
    record: &Record<'_>,
    get_whole: F,
) -> io::Result<Option<ReferenceSequence<'c>>>
where
    F: FnOnce() -> io::Result<Option<ReferenceSequence<'c>>>,
{
    match reference_sequence {
        ReferenceSequence::ExternalRange {
            reference_start,
            sequence,
        } if !reads_reference_within(record, *reference_start, sequence.len()) => get_whole(),
        _ => Ok(Some(reference_sequence.clone())),
    }
}

/// Returns the alignment start and the number of reference positions that reconstructing the
/// record's sequence reads from there.
///
/// Reference positions advance with the features and end at the alignment end, provided the
/// features are in order and do not run past the end of the read. This returns `None` for records
/// that are not.
fn reference_span(record: &Record<'_>) -> Option<(Position, usize)> {
    let alignment_start = record.alignment_start?;

    let mut read_end = 1;
    let mut alignment_span = record.read_length;

    for feature in &record.features {
        let read_len = match feature {
            Feature::Scores { .. } | Feature::QualityScore { .. } => continue,
            Feature::Bases { bases, .. } => bases.len(),
            Feature::ReadBase { .. } | Feature::Substitution { .. } => 1,
            Feature::Insertion { bases, .. } | Feature::SoftClip { bases, .. } => {
                alignment_span = alignment_span.checked_sub(bases.len())?;
                bases.len()
            }
            Feature::InsertBase { .. } => {
                alignment_span = alignment_span.checked_sub(1)?;
                1
            }
            Feature::Deletion { len, .. } | Feature::ReferenceSkip { len, .. } => {
                alignment_span = alignment_span.checked_add(*len)?;
                0
            }
            Feature::Padding { .. } | Feature::HardClip { .. } => 0,
        };

        let position = usize::from(feature.position());

        if position < read_end {
            return None;
        }

        read_end = position.checked_add(read_len)?;
    }

    (read_end <= record.read_length + 1).then_some((alignment_start, alignment_span))
}

/// Returns whether reconstructing the record's sequence reads only reference positions in the
/// `len` positions from `reference_start`.
fn reads_reference_within(record: &Record<'_>, reference_start: Position, len: usize) -> bool {
    reference_span(record).is_some_and(|(alignment_start, alignment_span)| {
        usize::from(alignment_start)
            .checked_sub(usize::from(reference_start))
            .and_then(|offset| offset.checked_add(alignment_span))
            .is_some_and(|end| end <= len)
    })
}

fn get_record_reference_sequence<'c>(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    record: &Record<'c>,
) -> io::Result<Option<ReferenceSequence<'c>>> {
    if record.bam_flags.is_unmapped() {
        return Ok(None);
    }

    let reference_sequence_name = record
        .reference_sequence(header)
        .transpose()?
        .map(|(name, _)| name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid record reference sequence ID",
            )
        })?;

    let sequence = reference_sequence_repository
        .get(reference_sequence_name)
        .transpose()?
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("missing reference sequence: {reference_sequence_name}"),
            )
        })?;

    Ok(Some(ReferenceSequence::External { sequence }))
}

/// Checks a reference sequence against its `@SQ` line, for records in slices without a reference
/// checksum (e.g., slices with many reference sequences).
///
/// Such records are otherwise decoded against whatever the reference holds: a wrong or truncated
/// one (e.g., another assembly's) would silently give wrong bases, or `N`s past its end. HTSlib
/// does not check them. Here, if the reference's sequence dictionary (see
/// [`fasta::Repository::load_dictionary`]) and the `@SQ` line both have a checksum, they must
/// match, and a match settles it, whatever `LN` says (pasta's CRAM writer, like HTSlib's, keeps
/// the input's `LN` beside the reference's `M5`). Otherwise the length is checked, and only a
/// length that differs from `LN` has the reference's checksum calculated, to settle it by `M5` as
/// the dictionary would; without `M5`, a reference shorter than `LN` is an error, and a longer
/// one gets HTSlib's warning. (HTSlib warns for both and decodes against the reference's length,
/// giving `N`s past its end.)
///
/// A matching length is not followed by calculating the checksum: hashing a whole chromosome on
/// first use took longer than decoding a small CRAM file.
fn validate_reference_sequence(
    reference_sequence_repository: &fasta::Repository,
    name: &[u8],
    reference_sequence: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> io::Result<()> {
    use crate::sam::header::record::value::map::reference_sequence::{Md5Checksum, tag};
    use bstr::ByteSlice;

    // A malformed checksum cannot be checked; the length is then.
    let expected_md5 = reference_sequence
        .other_fields()
        .get(&tag::MD5_CHECKSUM)
        .and_then(|value| value.to_str().ok()?.parse::<Md5Checksum>().ok());

    if let Some(expected_md5) = expected_md5
        && let Some(actual_md5) = reference_sequence_repository.dictionary_md5(name)
    {
        if actual_md5 != *expected_md5 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference sequence checksum mismatch for {}: the header (@SQ M5) has \
                     {expected_md5}, but the reference's sequence dictionary has {}; is it the \
                     wrong reference?",
                    name.as_bstr(),
                    Md5Checksum::from(actual_md5)
                ),
            ));
        }

        return Ok(());
    }

    let actual_len = reference_sequence_repository
        .sequence_len(name)
        .transpose()?
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("missing reference sequence: {}", name.as_bstr()),
            )
        })?;
    let expected_len = reference_sequence.length().get();

    // Only a length mismatch, which is rare, is worth calculating the checksum for.
    if actual_len != expected_len
        && let Some(expected_md5) = expected_md5
    {
        let actual_md5 = reference_sequence_repository
            .normalized_md5(name)
            .transpose()?
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("missing reference sequence: {}", name.as_bstr()),
                )
            })?;

        if actual_md5 != *expected_md5 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "reference sequence checksum mismatch for {}: the header (@SQ M5) has \
                     {expected_md5}, but the reference has {}; is it the wrong reference?",
                    name.as_bstr(),
                    Md5Checksum::from(actual_md5)
                ),
            ));
        }

        return Ok(());
    }

    if actual_len < expected_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "reference sequence length mismatch for {}: the header (@SQ LN) has \
                 {expected_len}, but the reference has {actual_len}; is it the wrong reference?",
                name.as_bstr()
            ),
        ));
    }

    if actual_len != expected_len {
        crate::variant::log::warn(
            "sanitise_SQ_lines",
            format_args!(
                "Header @SQ length mismatch for ref {}, {expected_len} vs {actual_len}",
                name.as_bstr()
            ),
        );
    }

    Ok(())
}

/// Checks the `start..=end` span of the reference sequence `name` (`sequence`) against the
/// slice's reference checksum, naming the region as HTSlib does ("MD5 checksum reference
/// mismatch at c1:10-39").
fn validate_sequence(
    name: &[u8],
    start: Position,
    end: Position,
    sequence: &[u8],
    expected_checksum: &[u8; 16],
) -> io::Result<()> {
    use crate::sam::header::record::value::map::reference_sequence::Md5Checksum;
    use bstr::ByteSlice;

    let actual_checksum = calculate_normalized_sequence_digest(sequence);

    if &actual_checksum == expected_checksum {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "reference sequence checksum mismatch at {}:{start}-{end}: the slice has {}, but \
                 the reference has {}; is it the wrong reference?",
                name.as_bstr(),
                Md5Checksum::from(*expected_checksum),
                Md5Checksum::from(actual_checksum),
            ),
        ))
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;
    use crate::cram::record::Flags;

    #[test]
    fn test_resolve_mates() -> io::Result<()> {
        let mut records = vec![
            Record {
                id: 1,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                bam_flags: sam::alignment::record::Flags::empty(),
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(5),
                mate_distance: Some(0),
                ..Default::default()
            },
            Record {
                id: 2,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                bam_flags: sam::alignment::record::Flags::empty(),
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(8),
                mate_distance: Some(1),
                ..Default::default()
            },
            Record {
                id: 3,
                ..Default::default()
            },
            Record {
                id: 4,
                bam_flags: sam::alignment::record::Flags::empty(),
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(13),
                ..Default::default()
            },
        ];

        resolve_mates(&mut records, Names::Generated)?;

        let name_1 = b"1".as_bstr();

        assert_eq!(records[0].name(), Some(name_1));
        assert_eq!(
            records[0].mate_reference_sequence_id,
            records[1].reference_sequence_id
        );
        assert_eq!(records[0].mate_alignment_start, records[1].alignment_start);
        assert_eq!(records[0].template_length, 12);

        assert_eq!(records[1].name(), Some(name_1));
        assert_eq!(
            records[1].mate_reference_sequence_id,
            records[3].reference_sequence_id
        );
        assert_eq!(records[1].mate_alignment_start, records[3].alignment_start);
        assert_eq!(records[1].template_length, -12);

        let name_3 = b"3".as_bstr();
        assert_eq!(records[2].name(), Some(name_3));

        assert_eq!(records[3].name(), Some(name_1));
        assert_eq!(
            records[3].mate_reference_sequence_id,
            records[0].reference_sequence_id
        );
        assert_eq!(records[3].mate_alignment_start, records[0].alignment_start);
        assert_eq!(records[3].template_length, -12);

        Ok(())
    }

    fn attached_pair<'c>(first: Record<'c>, second: Record<'c>) -> io::Result<[Record<'c>; 2]> {
        let mut records = [
            Record {
                id: 1,
                name: Some(Cow::from(b"r".to_vec())),
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                mate_distance: Some(0),
                ..first
            },
            Record { id: 2, ..second },
        ];

        resolve_mates(&mut records, Names::Generated)?;

        Ok(records)
    }

    #[test]
    fn test_resolve_mates_template_lengths() -> io::Result<()> {
        use sam::alignment::record::Flags as BamFlags;

        let mapped = |flags: BamFlags, start: usize, read_length: usize| Record {
            bam_flags: flags,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(start),
            read_length,
            ..Default::default()
        };

        // The leftmost segment gets the positive length.
        let [a, b] = attached_pair(
            mapped(BamFlags::LAST_SEGMENT, 100, 50),
            mapped(BamFlags::FIRST_SEGMENT, 200, 50),
        )?;
        assert_eq!((a.template_length, b.template_length), (150, -150));

        // Identical extents: the first segment gets the positive length, whatever the order.
        let [a, b] = attached_pair(
            mapped(BamFlags::LAST_SEGMENT, 100, 50),
            mapped(BamFlags::FIRST_SEGMENT, 100, 50),
        )?;
        assert_eq!((a.template_length, b.template_length), (-50, 50));

        // Same start, the first record is shorter: it is leftmost and not rightmost.
        let [a, b] = attached_pair(
            mapped(BamFlags::LAST_SEGMENT, 100, 10),
            mapped(BamFlags::FIRST_SEGMENT, 100, 50),
        )?;
        assert_eq!((a.template_length, b.template_length), (50, -50));

        // An unmapped segment gives both records 0.
        let unmapped = Record {
            bam_flags: BamFlags::UNMAPPED,
            ..mapped(BamFlags::LAST_SEGMENT, 100, 50)
        };
        let [a, b] = attached_pair(mapped(BamFlags::FIRST_SEGMENT, 100, 50), unmapped)?;
        assert_eq!((a.template_length, b.template_length), (0, 0));
        assert!(a.bam_flags.is_mate_unmapped());

        // Segments on different reference sequences give 0.
        let other = Record {
            reference_sequence_id: Some(1),
            ..mapped(BamFlags::LAST_SEGMENT, 300, 50)
        };
        let [a, b] = attached_pair(mapped(BamFlags::FIRST_SEGMENT, 100, 50), other)?;
        assert_eq!((a.template_length, b.template_length), (0, 0));

        // A mate past the end of the slice is an error.
        let mut records = [Record {
            mate_distance: Some(3),
            ..Default::default()
        }];
        assert!(matches!(
            resolve_mates(&mut records, Names::Generated),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_records_with_many_reference_sequences_and_no_required_reference()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::sam::alignment::{
            RecordBuf,
            io::Write as _,
            record::cigar::{Op, op::Kind},
        };

        use crate::cram::io::reader::Container;

        let header = "@SQ\tSN:sq0\tLN:4\n@SQ\tSN:sq1\tLN:4\n".parse()?;
        let repository = fasta::Repository::new(
            [b"sq0", b"sq1"]
                .map(|name| {
                    let definition = fasta::record::Definition::new(name, None);
                    fasta::Record::new(definition, b"ACGT".to_vec().into())
                })
                .to_vec(),
        );

        // A slice with a record on each reference sequence.
        let mut writer = crate::cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.many_references = true;
        writer.write_header(&header)?;
        for reference_sequence_id in [0, 1] {
            let record = RecordBuf::builder()
                .set_flags(sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(Position::MIN)
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(b"ACGT".to_vec().into())
                .set_quality_scores(vec![30; 4].into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }
        writer.try_finish(&header)?;
        let src = writer.into_inner();

        let mut reader = crate::cram::io::Reader::new(&src[..]);
        reader.read_header()?;
        let mut container = Container::default();
        reader.read_container(&mut container)?;
        let mut compression_header = container.compression_header()?;
        compression_header
            .preservation_map
            .external_reference_sequence_is_required = false;
        let slice = container.slices().next().unwrap()?;
        let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
        let records = slice.records(
            repository,
            &header,
            &compression_header,
            &core_data_src,
            &external_data_srcs,
        )?;

        // As in slices with one reference sequence (and HTSlib), records get no external
        // reference sequence when it is not required: their bases are stored verbatim.
        assert_eq!(records.len(), 2);
        for record in records {
            assert!(record.reference_sequence.is_none());
        }

        Ok(())
    }

    #[test]
    fn test_reads_reference_within() {
        let p = |n| Position::new(n).unwrap();
        let bases = |s: &'static [u8]| Cow::Borrowed(s);
        // 6 read bases: a soft clip, a 2-base deletion and a 2-base insertion span 5 positions.
        let features = vec![
            Feature::SoftClip {
                position: p(1),
                bases: bases(b"A"),
            },
            Feature::Deletion {
                position: p(3),
                len: 2,
            },
            Feature::Insertion {
                position: p(4),
                bases: bases(b"GG"),
            },
        ];
        let record = Record {
            read_length: 6,
            alignment_start: Some(p(10)),
            features,
            ..Default::default()
        };
        let span = crate::cram::record::calculate_alignment_span(6, &record.features);
        assert_eq!(reference_span(&record), Some((p(10), span)));
        for (start, len, expected) in [(10, 5, true), (8, 7, true), (10, 4, false), (11, 9, false)]
        {
            assert_eq!(reads_reference_within(&record, p(start), len), expected);
        }
        // Features out of order or past the end of the read are not trusted.
        let past_end = Feature::Bases {
            position: p(5),
            bases: bases(b"ACG"),
        };
        let unordered = vec![record.features[1].clone(), record.features[0].clone()];
        for features in [unordered, vec![past_end]] {
            let record = Record {
                features,
                ..record.clone()
            };
            assert!(!reads_reference_within(&record, p(1), 100));
        }
    }

    #[test]
    fn test_select_reference_sequence() -> io::Result<()> {
        let sequence = |s: &[u8]| Arc::new(fasta::record::Sequence::from(s.to_vec()));

        let range = ReferenceSequence::ExternalRange {
            reference_start: Position::MIN,
            sequence: sequence(b"ACGT"),
        };
        let whole = ReferenceSequence::External {
            sequence: sequence(b"ACGTACGT"),
        };
        let get_whole = || Ok(Some(whole.clone()));

        // The record reads reference positions 2-3.
        let record = Record {
            read_length: 2,
            alignment_start: Position::new(2),
            ..Default::default()
        };
        let actual = select_reference_sequence(&range, &record, get_whole)?;
        assert_eq!(actual, Some(range.clone()));

        // With a 2-base deletion, the record reads reference positions 2-5, past the range.
        let record = Record {
            features: vec![Feature::Deletion {
                position: Position::new(2).unwrap(),
                len: 2,
            }],
            ..record
        };
        let actual = select_reference_sequence(&range, &record, get_whole)?;
        assert_eq!(actual, Some(whole.clone()));

        // A whole reference sequence is used as is.
        let actual = select_reference_sequence(&whole, &record, || unreachable!())?;
        assert_eq!(actual, Some(whole));

        Ok(())
    }

    /// Slices with many reference sequences have no reference checksum, so their
    /// reference sequences are checked against `@SQ` `M5` (by the sequence dictionary's
    /// checksums, or by calculating them for a length that differs) or `LN`, instead of decoding
    /// against a wrong or truncated reference (as `N`s past its end).
    #[test]
    fn test_many_reference_sequences_check_the_reference_against_the_header()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::sam::alignment::{
            RecordBuf,
            io::Write as _,
            record::cigar::{Op, op::Kind},
        };

        let repository = |sq1: &[u8]| {
            fasta::Repository::new(vec![
                fasta::Record::new(
                    fasta::record::Definition::new("sq0", None),
                    b"ACGTACGT".to_vec().into(),
                ),
                fasta::Record::new(
                    fasta::record::Definition::new("sq1", None),
                    sq1.to_vec().into(),
                ),
            ])
        };

        let header: sam::Header = "@SQ\tSN:sq0\tLN:8\n@SQ\tSN:sq1\tLN:8\n".parse()?;

        // A slice with records on both reference sequences, the second reading to the end of sq1.
        let mut writer = crate::cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository(b"TTGCAAGG"))
            .build_from_writer(Vec::new());
        writer.many_references = true;
        writer.write_header(&header)?;
        for (reference_sequence_id, start, bases) in [(0, 1, b"ACGT"), (1, 5, b"AAGG")] {
            let record = RecordBuf::builder()
                .set_flags(sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(Position::new(start).unwrap())
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(bases.to_vec().into())
                .set_quality_scores(vec![30; 4].into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }
        writer.try_finish(&header)?;
        let src = writer.into_inner();

        // A sequence dictionary of the reference given, as `ntools dict` or Picard would write it.
        let dir = std::env::temp_dir().join(format!("pasta-slice-dict-{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let fasta_path = dir.join("ref.fa");

        let read_with = |sq1: &[u8],
                         strip_checksums: bool,
                         sq1_len: usize,
                         dictionary: bool|
         -> io::Result<Vec<Vec<u8>>> {
            let mut reader = crate::cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository(sq1))
                .build_from_reader(&src[..]);
            let mut header = reader.read_header()?;

            let mut slices = 0;
            let mut container = crate::cram::io::reader::Container::default();
            while reader.read_container(&mut container)? != 0 {
                for slice in container.slices() {
                    assert!(slice?.header().reference_sequence_context().is_many());
                    slices += 1;
                }
            }
            assert_eq!(slices, 1);

            if strip_checksums {
                for (_, reference_sequence) in header.reference_sequences_mut() {
                    reference_sequence.other_fields_mut().clear();
                }
            }

            if let Some(reference_sequence) = header.reference_sequences_mut().get_mut(&b"sq1"[..])
            {
                *reference_sequence.length_mut() = std::num::NonZero::new(sq1_len).unwrap();
            }

            let repository = repository(sq1);
            if dictionary {
                let md5 = |s: &[u8]| {
                    crate::sam::header::record::value::map::reference_sequence::Md5Checksum::from(
                        crate::cram::calculate_normalized_sequence_digest(s),
                    )
                };
                std::fs::write(
                    dir.join("ref.dict"),
                    format!(
                        "@HD\tVN:1.0\n@SQ\tSN:sq0\tLN:8\tM5:{}\n@SQ\tSN:sq1\tLN:{}\tM5:{}\n",
                        md5(b"ACGTACGT"),
                        sq1.len(),
                        md5(sq1)
                    ),
                )?;
                assert!(repository.load_dictionary(&fasta_path)?.is_some());
            }

            let mut reader = crate::cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_reader(&src[..]);
            reader.read_header()?;
            reader
                .records(&header)
                .map(|record| {
                    let record = record?;
                    Ok(record.sequence().as_ref().to_vec())
                })
                .collect()
        };

        let read = |sq1: &[u8], strip_checksums: bool, sq1_len: usize| {
            read_with(sq1, strip_checksums, sq1_len, false)
        };

        let expected = [b"ACGT".to_vec(), b"AAGG".to_vec()];
        assert_eq!(read(b"TTGCAAGG", false, 8)?, expected);
        assert_eq!(read_with(b"TTGCAAGG", false, 8, true)?, expected);
        assert_eq!(
            read_with(b"ttgcaagg", false, 8, true)?.len(),
            2,
            "M5 is case-insensitive"
        );

        // A matching M5 settles an LN that differs from the reference, as samtools reads it.
        assert_eq!(read(b"TTGCAAGG", false, 9)?, expected);
        // Without M5, a reference longer than LN is read, with HTSlib's warning.
        assert_eq!(read(b"TTGCAAGG", true, 7)?, expected);

        // A reference shorter than LN (formerly read as AANN): its checksum differs, and without
        // M5, its length does.
        let e = read(b"TTGCAA", false, 8).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
        assert!(e.to_string().contains("checksum mismatch for sq1"), "{e}");
        let e = read(b"TTGCAA", true, 8).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
        assert!(e.to_string().contains("length mismatch for sq1"), "{e}");

        // A reference of the right length but other bases: the sequence dictionary's checksum
        // tells; without one, it is read (as substitutions), as HTSlib reads it.
        let e = read_with(b"TTGCAAGC", false, 8, true).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
        assert!(e.to_string().contains("sequence dictionary has"), "{e}");
        assert_eq!(read(b"TTGCAAGC", false, 8)?.len(), 2);

        std::fs::remove_dir_all(&dir)?;

        Ok(())
    }

    #[test]
    fn test_slice_checksum_mismatch_names_the_region() -> Result<(), Box<dyn std::error::Error>> {
        use crate::sam::alignment::{
            RecordBuf,
            io::Write as _,
            record::cigar::{Op, op::Kind},
        };

        let repository = |sq0: &[u8]| {
            fasta::Repository::new(vec![fasta::Record::new(
                fasta::record::Definition::new("sq0", None),
                sq0.to_vec().into(),
            )])
        };

        let header: sam::Header = "@SQ\tSN:sq0\tLN:16\n".parse()?;

        let mut writer = crate::cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repository(b"ACGTACGTACGTACGT"))
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;
        let record = RecordBuf::builder()
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(3).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
            .set_sequence(b"GTAC".to_vec().into())
            .set_quality_scores(vec![30; 4].into())
            .build();
        writer.write_alignment_record(&header, &record)?;
        writer.try_finish(&header)?;
        let src = writer.into_inner();

        // Another reference of the same length: the slice's checksum of sq0:3-6 does not match.
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository(b"ACCAACGTACGTACGT"))
            .build_from_reader(&src[..]);
        // The header the records are read with has no M5, so the slice checksum is what fails.
        reader.read_header()?;
        let e = reader
            .records(&header)
            .collect::<io::Result<Vec<_>>>()
            .unwrap_err();

        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
        let message = e.to_string();
        assert!(message.contains("mismatch at sq0:3-6"), "{message}");
        assert!(message.contains("is it the wrong reference?"), "{message}");
        let expected = crate::cram::calculate_normalized_sequence_digest(b"GTAC");
        let expected =
            sam::header::record::value::map::reference_sequence::Md5Checksum::from(expected);
        assert!(
            message.contains(&format!("the slice has {expected}")),
            "{message}"
        );

        Ok(())
    }
}
