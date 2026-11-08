use std::{io, ops::Deref};

use bstr::BStr;

use crate::core::Position;
use crate::sam::{
    self as sam,
    alignment::record::{
        MappingQuality,
        data::field::{Tag, Type, Value},
    },
};
use crate::{bam, fasta};

use super::{Data, Feature, Record, ReferenceBases};
use crate::cram::{
    container::compression_header::preservation_map::substitution_matrix::Base,
    record::{Flags, MateFlags},
};

impl Record {
    /// Converts an alignment record, comparing its bases with the sequence that
    /// `get_reference_sequence` returns for its reference sequence name. When that is `None`, the
    /// bases are stored verbatim.
    pub fn try_from_alignment_record<F, R>(
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
        get_reference_sequence: F,
    ) -> io::Result<Self>
    where
        F: FnOnce(&[u8]) -> io::Result<Option<R>>,
        R: Deref<Target = fasta::record::Sequence>,
    {
        let cigar = record.cigar();
        let sequence = record.sequence();
        let quality_scores = record.quality_scores();
        let (data, read_group_id) =
            get_filtered_data(header, record.data().as_ref(), &mut ReadGroups::default())?;

        let parts = Parts {
            bam_flags: record.flags()?,
            reference_sequence_id: record.reference_sequence_id(header).transpose()?,
            alignment_start: record.alignment_start().transpose()?,
            mapping_quality: record.mapping_quality().transpose()?,
            mate_reference_sequence_id: record.mate_reference_sequence_id(header).transpose()?,
            mate_alignment_start: record.mate_alignment_start().transpose()?,
            template_length: record.template_length()?,
            name: record.name(),
            sequence: sequence.as_ref(),
            quality_scores: quality_scores.as_ref(),
            data,
            read_group_id,
        };

        parts.into_record(
            header,
            cigar.iter(),
            || cigar.read_length(),
            get_reference_sequence,
        )
    }

    /// Converts a BAM record, as [`Self::try_from_alignment_record`] does, reading its fields in
    /// place and taking its tags as encoded. `read_groups` remembers the last read group looked
    /// up.
    pub(crate) fn try_from_bam_record<F, R>(
        header: &sam::Header,
        record: &bam::Record,
        read_groups: &mut ReadGroups,
        get_reference_sequence: F,
    ) -> io::Result<Self>
    where
        F: FnOnce(&[u8]) -> io::Result<Option<R>>,
        R: Deref<Target = fasta::record::Sequence>,
    {
        let cigar = record.cigar();
        let sequence = record.sequence();
        let quality_scores = record.quality_scores();

        let mut read_group_id = None;

        let data = Data::from_bam(record.data().as_ref(), |field| {
            if field.tag == Tag::READ_GROUP
                && read_group_id.is_none()
                && let Some(name) = field.as_str().filter(|_| field.ty == Type::String)
                && let Some(id) = read_groups.get(header, name)
            {
                read_group_id = Some(id);
                return false;
            }

            true
        })?;

        let parts = Parts {
            bam_flags: record.flags(),
            reference_sequence_id: record.reference_sequence_id().transpose()?,
            alignment_start: record.alignment_start().transpose()?,
            mapping_quality: record.mapping_quality(),
            mate_reference_sequence_id: record.mate_reference_sequence_id().transpose()?,
            mate_alignment_start: record.mate_alignment_start().transpose()?,
            template_length: record.template_length(),
            name: record.name(),
            sequence: &sequence,
            quality_scores: &quality_scores,
            data,
            read_group_id,
        };

        parts.into_record(
            header,
            cigar.iter(),
            || sam::alignment::record::Cigar::read_length(&cigar),
            get_reference_sequence,
        )
    }
}

/// The fields of an alignment record that a CRAM record is converted from.
struct Parts<'a> {
    bam_flags: sam::alignment::record::Flags,
    reference_sequence_id: Option<usize>,
    alignment_start: Option<Position>,
    mapping_quality: Option<MappingQuality>,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
    template_length: i32,
    name: Option<&'a BStr>,
    sequence: &'a dyn sam::alignment::record::Sequence,
    quality_scores: &'a dyn sam::alignment::record::QualityScores,
    data: Data,
    read_group_id: Option<usize>,
}

impl Parts<'_> {
    /// Builds the record, with the CIGAR operations `cigar` (whose read length
    /// `cigar_read_length` gives).
    fn into_record<I, L, F, R>(
        self,
        header: &sam::Header,
        cigar: I,
        cigar_read_length: L,
        get_reference_sequence: F,
    ) -> io::Result<Record>
    where
        I: IntoIterator<Item = io::Result<sam::alignment::record::cigar::Op>>,
        L: FnOnce() -> io::Result<usize>,
        F: FnOnce(&[u8]) -> io::Result<Option<R>>,
        R: Deref<Target = fasta::record::Sequence>,
    {
        let bam_flags = self.bam_flags;
        let mut cram_flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;

        let record_sequence = self.sequence;
        let record_quality_scores = self.quality_scores;

        if record_sequence.is_empty() {
            cram_flags.insert(Flags::SEQUENCE_IS_MISSING);
        }

        // As in HTSlib, a mapped record without a sequence (`*`) has the read length of its
        // CIGAR, so decoders can rebuild the CIGAR from the features.
        let read_length = if record_sequence.is_empty() && !bam_flags.is_unmapped() {
            cigar_read_length()?
        } else {
            record_sequence.len()
        };

        if !record_quality_scores.is_empty() && record_quality_scores.len() != record_sequence.len()
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "sequence-quality scores length mismatch: expected {}, got {}",
                    record_sequence.len(),
                    record_quality_scores.len()
                ),
            ));
        }

        // The bases, the quality scores, then the name, in one buffer. Quality scores are stored
        // as an array of the read length; missing ones (`*`) are stored as 0xff, as in HTSlib.
        let name_len = self.name.map(|name| name.len());
        let mut bases =
            Vec::with_capacity(record_sequence.len() + read_length + name_len.unwrap_or(0));
        record_sequence.append_to(&mut bases);
        let sequence_len = bases.len();

        if record_quality_scores.is_empty() {
            bases.resize(sequence_len + read_length, 0xff);
        } else {
            record_quality_scores.append_to(&mut bases)?;
        }

        let (sequence, quality_scores) = bases.split_at(sequence_len);

        let (features, reference_bases) =
            if let (Some(id), Some(start)) = (self.reference_sequence_id, self.alignment_start) {
                let (reference_sequence_name, _) =
                    header.reference_sequences().get_index(id).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
                    })?;

                let reference_sequence = get_reference_sequence(reference_sequence_name)?;

                let reference_bases = if reference_sequence.is_some() {
                    ReferenceBases::External
                } else {
                    ReferenceBases::Verbatim
                };

                let features = cigar_to_features(
                    cigar,
                    reference_sequence.as_deref(),
                    cram_flags,
                    start,
                    sequence,
                    quality_scores,
                )?;

                (features, reference_bases)
            } else {
                (Vec::new(), ReferenceBases::Unused)
            };

        if let Some(name) = self.name {
            bases.extend_from_slice(name);
        }

        Ok(Record {
            bam_flags,
            cram_flags,
            reference_sequence_id: self.reference_sequence_id,
            read_length,
            alignment_start: self.alignment_start,
            read_group_id: self.read_group_id,
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: self.mate_reference_sequence_id,
            mate_alignment_start: self.mate_alignment_start,
            template_length: self.template_length,
            mate_distance: None,
            data: self.data,
            features,
            reference_bases,
            mapping_quality: self.mapping_quality,
            bases,
            sequence_len,
            name_len,
        })
    }
}

/// The index of the last read group found by name, as consecutive records usually share it.
#[derive(Debug, Default)]
pub(crate) struct ReadGroups {
    last: Option<usize>,
}

impl ReadGroups {
    /// The index in `header` of the read group named `name`.
    ///
    /// The last index found is checked against `header` rather than trusted, as records may come
    /// with another header.
    fn get(&mut self, header: &sam::Header, name: &[u8]) -> Option<usize> {
        let read_groups = header.read_groups();

        if let Some(id) = self.last
            && read_groups
                .get_index(id)
                .is_some_and(|(last_name, _)| last_name[..] == name[..])
        {
            return Some(id);
        }

        let id = read_groups.get_index_of(name);

        if id.is_some() {
            self.last = id;
        }

        id
    }
}

/// Returns the read length a record converts with (see [`Record::try_from_alignment_record`])
/// without converting it.
pub(crate) fn read_length(record: &bam::Record) -> io::Result<usize> {
    let sequence = record.sequence();

    if sequence.is_empty() && !record.flags().is_unmapped() {
        sam::alignment::record::Cigar::read_length(&record.cigar())
    } else {
        Ok(sequence.len())
    }
}

/// Builds the features of a record's CIGAR. Without a reference sequence, aligned bases are stored
/// verbatim, as stretches of bases.
fn cigar_to_features<I>(
    cigar: I,
    reference_sequence: Option<&fasta::record::Sequence>,
    flags: Flags,
    alignment_start: Position,
    bases: &[u8],
    scores: &[u8],
) -> io::Result<Vec<Feature>>
where
    I: IntoIterator<Item = io::Result<sam::alignment::record::cigar::Op>>,
{
    use crate::sam::alignment::record::cigar::op::Kind;

    let reference_bases: Option<&[u8]> = reference_sequence.map(|sequence| sequence.as_ref());
    let stores_scores = !flags.quality_scores_are_stored_as_array();

    // Many reads match the reference throughout and have none.
    let mut features = Vec::new();

    // 0-based offsets into the reference sequence and the read.
    let mut reference_offset = usize::from(alignment_start) - 1;
    let mut read_offset: usize = 0;

    for result in cigar {
        let op = result?;
        let len = op.len();
        let kind = op.kind();

        let read_end = read_offset
            .checked_add(len)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "CIGAR is too long"))?;

        // A missing sequence (`*`) has no bases to compare or store; HTSlib writes `N` for its
        // inserted and soft-clipped bases.
        if kind.consumes_read() && !bases.is_empty() && read_end > bases.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "CIGAR and query sequence are of different length",
            ));
        }

        let position = read_position(read_offset)?;

        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if !bases.is_empty() {
                    if stores_scores {
                        push_scores(&mut features, position, scores, read_offset..read_end);
                    }

                    if let Some(reference_bases) = reference_bases {
                        push_mismatches(
                            &mut features,
                            reference_bases,
                            reference_offset,
                            &bases[read_offset..read_end],
                            scores,
                            read_offset,
                        )?;
                    } else if len > 0 {
                        features.push(Feature::Bases {
                            position,
                            bases: bases[read_offset..read_end].to_vec(),
                        });
                    }
                }
            }
            Kind::Insertion => {
                let inserted = read_bases(bases, read_offset..read_end);

                if let [base] = inserted[..] {
                    features.push(Feature::InsertBase { position, base });
                } else {
                    features.push(Feature::Insertion {
                        position,
                        bases: inserted,
                    });
                }

                if stores_scores {
                    push_scores(&mut features, position, scores, read_offset..read_end);
                }
            }
            Kind::Deletion => features.push(Feature::Deletion { position, len }),
            Kind::Skip => features.push(Feature::ReferenceSkip { position, len }),
            Kind::SoftClip => {
                features.push(Feature::SoftClip {
                    position,
                    bases: read_bases(bases, read_offset..read_end),
                });

                if stores_scores {
                    push_scores(&mut features, position, scores, read_offset..read_end);
                }
            }
            Kind::HardClip => features.push(Feature::HardClip { position, len }),
            Kind::Pad => features.push(Feature::Padding { position, len }),
        }

        if kind.consumes_reference() {
            reference_offset = reference_offset
                .checked_add(len)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "CIGAR is too long"))?;
        }

        if kind.consumes_read() {
            read_offset = read_end;
        }
    }

    Ok(features)
}

fn read_position(offset: usize) -> io::Result<Position> {
    offset
        .checked_add(1)
        .and_then(Position::new)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read position"))
}

/// Returns the read bases in `range`, or `N`s when the sequence is missing.
fn read_bases(bases: &[u8], range: std::ops::Range<usize>) -> Vec<u8> {
    if bases.is_empty() {
        vec![b'N'; range.len()]
    } else {
        bases[range].to_vec()
    }
}

/// Pushes a quality score feature for the scores in `range`, if the read has them.
fn push_scores(
    features: &mut Vec<Feature>,
    position: Position,
    scores: &[u8],
    range: std::ops::Range<usize>,
) {
    let Some(scores) = scores.get(range) else {
        return;
    };

    if let [quality_score] = *scores {
        features.push(Feature::QualityScore {
            position,
            quality_score,
        });
    } else {
        features.push(Feature::Scores {
            position,
            quality_scores: scores.to_vec(),
        });
    }
}

/// Pushes substitutions (or verbatim read bases) for the aligned bases `read_bases`, which start
/// at `read_offset` in the read and `reference_offset` in the reference.
///
/// As in HTSlib, bases past the end of the reference are stored verbatim.
pub(crate) fn push_mismatches(
    features: &mut Vec<Feature>,
    reference_bases: &[u8],
    reference_offset: usize,
    read_bases: &[u8],
    scores: &[u8],
    read_offset: usize,
) -> io::Result<()> {
    let reference_bases = reference_bases.get(reference_offset..).unwrap_or_default();
    let n = reference_bases.len().min(read_bases.len());

    // Most bases match: compare chunks first (whole arrays, which compare inline), and bases only
    // in chunks that differ.
    const CHUNK_LEN: usize = 16;

    let mut push_chunk = |start: usize, reference_chunk: &[u8], read_chunk: &[u8]| {
        for (j, (&reference_base, &read_base)) in reference_chunk.iter().zip(read_chunk).enumerate()
        {
            if !reference_base.eq_ignore_ascii_case(&read_base) {
                let i = start + j;
                let position = read_position(read_offset + i)?;
                features.push(mismatch(position, reference_base, read_base, || {
                    scores.get(read_offset + i).copied()
                }));
            }
        }

        Ok::<_, io::Error>(())
    };

    let (reference_chunks, reference_rest) = reference_bases[..n].as_chunks::<CHUNK_LEN>();
    let (read_chunks, read_rest) = read_bases[..n].as_chunks::<CHUNK_LEN>();

    for (chunk_index, (reference_chunk, read_chunk)) in
        reference_chunks.iter().zip(read_chunks).enumerate()
    {
        if reference_chunk != read_chunk {
            push_chunk(chunk_index * CHUNK_LEN, reference_chunk, read_chunk)?;
        }
    }

    push_chunk(
        reference_chunks.len() * CHUNK_LEN,
        reference_rest,
        read_rest,
    )?;

    // Bases past the end of the reference.
    for (i, &base) in read_bases.iter().enumerate().skip(n) {
        features.push(Feature::ReadBase {
            position: read_position(read_offset + i)?,
            base,
            quality_score: scores.get(read_offset + i).copied().unwrap_or(0xff),
        });
    }

    Ok(())
}

/// A substitution, or the read base verbatim if either base cannot be substituted.
fn mismatch<F>(position: Position, reference_base: u8, read_base: u8, score: F) -> Feature
where
    F: FnOnce() -> Option<u8>,
{
    match (Base::try_from(reference_base), Base::try_from(read_base)) {
        (Ok(reference_base), Ok(read_base)) => Feature::Substitution {
            position,
            reference_base,
            read_base,
        },
        // HTSlib stores missing quality scores (`*`) as 0xff.
        _ => Feature::ReadBase {
            position,
            base: read_base,
            quality_score: score().unwrap_or(0xff),
        },
    }
}

/// Returns the record's tags as encoded, but for an `RG` tag naming a read group in the header,
/// whose index is returned instead.
///
/// As in HTSlib, an `RG` field naming a read group that the header lacks (or that is not a
/// string) is stored verbatim as a tag.
fn get_filtered_data(
    header: &sam::Header,
    data: &dyn sam::alignment::record::Data<'_>,
    read_groups: &mut ReadGroups,
) -> io::Result<(Data, Option<usize>)> {
    let mut data_buf = Data::default();
    let mut read_group_id = None;

    for result in data.iter() {
        let (tag, value) = result?;

        if tag == Tag::READ_GROUP
            && read_group_id.is_none()
            && let Value::String(name) = &value
            && let Some(id) = read_groups.get(header, name)
        {
            read_group_id = Some(id);
            continue;
        }

        data_buf.push(tag, &value.try_into()?)?;
    }

    Ok((data_buf, read_group_id))
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::sam::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::{Cigar, data::field::Value as ValueBuf},
    };

    use super::*;
    use crate::sam::alignment::record_buf::{QualityScores, Sequence};

    fn mapped_record(
        cigar: &str,
        sequence: &[u8],
        quality_scores: &[u8],
    ) -> sam::alignment::RecordBuf {
        use sam::alignment::record_buf::Cigar;

        let cigar: Cigar = cigar
            .split_inclusive(|c: char| c.is_ascii_alphabetic() || c == '=')
            .map(|s| {
                let (n, op) = s.split_at(s.len() - 1);
                let kind = match op {
                    "M" => Kind::Match,
                    "I" => Kind::Insertion,
                    "D" => Kind::Deletion,
                    "S" => Kind::SoftClip,
                    _ => unreachable!(),
                };
                Op::new(kind, n.parse().unwrap())
            })
            .collect();

        sam::alignment::RecordBuf::builder()
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_cigar(cigar)
            .set_sequence(Sequence::from(sequence.to_vec()))
            .set_quality_scores(QualityScores::from(quality_scores.to_vec()))
            .build()
    }

    fn convert(record: &sam::alignment::RecordBuf) -> io::Result<Record> {
        let header =
            sam::Header::builder()
                .add_reference_sequence(
                    "sq0",
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(std::num::NonZero::new(4).unwrap()),
                )
                .build();

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGT".to_vec()),
        )]);

        convert_both(&header, record, |name| repository.get(name).transpose())
    }

    /// Converts `record` as an alignment record and as the BAM record it encodes to, which must
    /// give the same record.
    fn convert_both<F>(
        header: &sam::Header,
        record: &sam::alignment::RecordBuf,
        get_reference_sequence: F,
    ) -> io::Result<Record>
    where
        F: Fn(&[u8]) -> io::Result<Option<Arc<fasta::record::Sequence>>>,
    {
        let actual = Record::try_from_alignment_record(header, record, &get_reference_sequence)?;

        let mut buf = Vec::new();
        crate::bam::record::codec::encoder::encode(&mut buf, header, record)?;
        let bam_record = bam::Record::try_from(buf)?;
        let from_bam = Record::try_from_bam_record(
            header,
            &bam_record,
            &mut ReadGroups::default(),
            &get_reference_sequence,
        )?;
        assert_eq!(from_bam, actual);

        Ok(actual)
    }

    fn data_of(fields: &[(Tag, ValueBuf)]) -> Data {
        let mut data = Data::default();

        for (tag, value) in fields {
            data.push(*tag, value).unwrap();
        }

        data
    }

    #[test]
    fn test_read_groups_follow_the_header() -> Result<(), Box<dyn std::error::Error>> {
        let a: sam::Header = "@RG\tID:x\n@RG\tID:y\n".parse()?;
        let b: sam::Header = "@RG\tID:y\n@RG\tID:x\n".parse()?;
        let mut read_groups = ReadGroups::default();

        assert_eq!(read_groups.get(&a, b"y"), Some(1));
        // The same name in another header, at another index.
        assert_eq!(read_groups.get(&b, b"y"), Some(0));
        assert_eq!(read_groups.get(&b, b"x"), Some(1));
        assert_eq!(read_groups.get(&b, b"z"), None);
        assert_eq!(read_groups.get(&a, b"x"), Some(0));

        Ok(())
    }

    #[test]
    fn missing_sequence_takes_the_cigar_read_length() -> io::Result<()> {
        // SEQ and QUAL `*`, as HTSlib stores them: the read length is the CIGAR's, qualities
        // are 0xff, and inserted and soft-clipped bases are `N`.
        let record = convert(&mapped_record("1S2M2I", b"", b""))?;

        assert!(record.cram_flags.contains(Flags::SEQUENCE_IS_MISSING));
        assert_eq!(record.read_length, 5);
        assert_eq!(record.quality_scores(), [0xff; 5]);
        assert_eq!(
            record.features,
            [
                Feature::SoftClip {
                    position: Position::MIN,
                    bases: b"N".to_vec(),
                },
                Feature::Insertion {
                    position: Position::try_from(4).unwrap(),
                    bases: b"NN".to_vec(),
                },
            ]
        );

        Ok(())
    }

    #[test]
    fn missing_quality_scores_are_stored_as_0xff() -> io::Result<()> {
        let record = convert(&mapped_record("2M", b"AC", b""))?;
        assert_eq!(record.read_length, 2);
        assert_eq!(record.quality_scores(), [0xff; 2]);
        assert!(record.features.is_empty());
        Ok(())
    }

    #[test]
    fn bases_past_the_reference_end_are_stored_verbatim() -> io::Result<()> {
        // The reference is ACGT; the read covers positions 3 to 6.
        let mut record = mapped_record("4M", b"GTCA", b"ABCD");
        *record.alignment_start_mut() = Position::new(3);
        let record = convert(&record)?;

        assert_eq!(
            record.features,
            [
                Feature::ReadBase {
                    position: Position::try_from(3).unwrap(),
                    base: b'C',
                    quality_score: b'C',
                },
                Feature::ReadBase {
                    position: Position::try_from(4).unwrap(),
                    base: b'A',
                    quality_score: b'D',
                },
            ]
        );

        Ok(())
    }

    #[test]
    fn quality_scores_shorter_than_the_sequence_are_an_error() {
        assert!(matches!(
            convert(&mapped_record("4M", b"ACGT", b"AB")),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }

    #[test]
    fn cigar_longer_than_the_sequence_is_an_error() {
        assert!(matches!(
            convert(&mapped_record("4M", b"AC", b"")),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }

    #[test]
    fn read_group_uses_dedicated_series_without_duplicate_auxiliary_tag() -> io::Result<()> {
        let header = sam::Header::builder()
            .add_read_group("rg0", Default::default())
            .build();
        let data = [(Tag::READ_GROUP, "rg0"), (Tag::COMMENT, "note")];
        let data = data
            .into_iter()
            .map(|(tag, s)| (tag, ValueBuf::from(s)))
            .collect();
        let record = sam::alignment::RecordBuf::builder().set_data(data).build();
        let actual = convert_both(&header, &record, |_| Ok(None))?;
        assert_eq!(actual.read_group_id, Some(0));
        assert_eq!(
            actual.data,
            data_of(&[(Tag::COMMENT, ValueBuf::from("note"))])
        );
        Ok(())
    }

    #[test]
    fn read_group_missing_from_the_header_is_stored_as_a_tag() -> io::Result<()> {
        let header = sam::Header::builder()
            .add_read_group("rg0", Default::default())
            .build();
        let data = [(Tag::READ_GROUP, "nosuch"), (Tag::COMMENT, "note")];
        let data: sam::alignment::record_buf::Data = data
            .into_iter()
            .map(|(tag, s)| (tag, ValueBuf::from(s)))
            .collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_data(data.clone())
            .build();
        let actual = convert_both(&header, &record, |_| Ok(None))?;
        assert_eq!(actual.read_group_id, None);
        let expected: Vec<_> = data
            .iter()
            .map(|(tag, value)| (tag, value.clone()))
            .collect();
        assert_eq!(actual.data, data_of(&expected));
        Ok(())
    }

    #[test]
    fn test_cigar_to_features() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGT".to_vec());

        let flags = Flags::default();
        let alignment_start = Position::MIN;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![Feature::QualityScore {
            position: Position::try_from(1)?,
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"C");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![Feature::Scores {
            position: Position::try_from(1)?,
            quality_scores: vec![45, 35],
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AT");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::C,
                read_base: Base::T,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::QualityScore {
                position: Position::try_from(2)?,
                quality_score: 35,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores: QualityScores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 43,
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::G,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::QualityScore {
                position: Position::try_from(2)?,
                quality_score: 35,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores: QualityScores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 43,
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::HardClip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::Padding {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_cigar_to_features_with_quality_scores_stored_as_array()
    -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGT".to_vec());

        let flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;
        let alignment_start = Position::MIN;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"C");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![Feature::Substitution {
            position: Position::try_from(1)?,
            reference_base: Base::A,
            read_base: Base::C,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AT");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![Feature::Substitution {
            position: Position::try_from(2)?,
            reference_base: Base::C,
            read_base: Base::T,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores: QualityScores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::G,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores: QualityScores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores: QualityScores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![Feature::HardClip {
            position: Position::try_from(1)?,
            len: 1,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores: QualityScores = [45].into_iter().collect();
        let actual = cigar_to_features(
            cigar.as_ref().iter().copied().map(Ok),
            Some(&reference_sequence),
            flags,
            alignment_start,
            sequence.as_ref(),
            quality_scores.as_ref(),
        )?;
        let expected = vec![Feature::Padding {
            position: Position::try_from(1)?,
            len: 1,
        }];
        assert_eq!(actual, expected);

        Ok(())
    }
}
