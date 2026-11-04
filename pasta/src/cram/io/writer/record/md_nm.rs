//! MD and NM tags left out because readers regenerate them, as HTSlib's writer does.

use crate::core::Position;
use crate::sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};

use super::{Feature, Record, ReferenceBases};
use crate::cram::record::md_nm::{Reference, Step, calculate_from_steps};

impl Record {
    /// Removes the MD and NM tags that readers regenerate with the same values, given the
    /// record's external `reference_sequence` and its length in the header.
    ///
    /// As in HTSlib, both are kept for alignments that reach past the end of the reference
    /// sequence or align an N in the read to an N in the reference, where readers disagree.
    /// `store_md` and `store_nm` keep them always (HTSlib's `store_md` and `store_nm`).
    pub(crate) fn remove_regenerated_md_nm(
        &mut self,
        reference_sequence: &[u8],
        reference_sequence_length: usize,
        store_md: bool,
        store_nm: bool,
    ) {
        let has = |tag| self.data.iter().filter(|(t, _)| *t == tag).count() == 1;
        let has_md = !store_md && has(Tag::MISMATCHED_POSITIONS);
        let has_nm = !store_nm && has(Tag::EDIT_DISTANCE);

        if !(has_md || has_nm)
            || self.reference_bases != ReferenceBases::External
            || self.bam_flags.is_unmapped()
            || self.sequence.is_empty()
        {
            return;
        }

        let (Some(start), Some(end)) = (self.alignment_start, self.alignment_end()) else {
            return;
        };

        if usize::from(end) > reference_sequence_length.min(reference_sequence.len())
            || self.aligns_n_to_n(reference_sequence, start)
        {
            return;
        }

        let reference = Reference {
            bases: reference_sequence,
            offset: 0,
            len: reference_sequence_length,
        };

        let steps = self
            .features
            .iter()
            .map(|feature| (usize::from(feature.position()), Step::from(feature)));

        let Some((md, nm)) =
            calculate_from_steps(&reference, start, self.read_length, steps, has_md)
        else {
            return;
        };

        self.data.retain(|(tag, value)| match *tag {
            Tag::MISMATCHED_POSITIONS => {
                store_md || !matches!(value, Value::String(s) if s[..] == md[..])
            }
            Tag::EDIT_DISTANCE => store_nm || value.as_int() != Some(i64::from(nm)),
            _ => true,
        });
    }
}

impl Record {
    /// Whether an aligned N in the read is on an N in `reference_sequence`, which holds the whole
    /// alignment from `start`.
    fn aligns_n_to_n(&self, reference_sequence: &[u8], start: Position) -> bool {
        let is_n = |b: &u8| b & !0x20 == b'N';

        let both_n = |read_pos: usize, ref_pos: usize, len: usize| {
            let read_bases = self
                .sequence
                .get(read_pos..read_pos + len)
                .unwrap_or_default();
            let ref_bases = &reference_sequence[ref_pos..ref_pos + len];
            read_bases
                .iter()
                .zip(ref_bases)
                .any(|(a, b)| is_n(a) && is_n(b))
        };

        // 0-based positions.
        let mut read_pos = 0;
        let mut ref_pos = usize::from(start) - 1;

        for feature in &self.features {
            let pos = usize::from(feature.position()) - 1;

            if pos > read_pos {
                let len = pos - read_pos;

                if both_n(read_pos, ref_pos, len) {
                    return true;
                }

                read_pos = pos;
                ref_pos += len;
            }

            match Step::from(feature) {
                Step::Substitution | Step::ReadBase => {
                    if both_n(read_pos, ref_pos, 1) {
                        return true;
                    }

                    read_pos += 1;
                    ref_pos += 1;
                }
                Step::Bases(len) => {
                    if both_n(read_pos, ref_pos, len) {
                        return true;
                    }

                    read_pos += len;
                    ref_pos += len;
                }
                Step::SoftClip(len) | Step::Insertion(len) => read_pos += len,
                Step::InsertBase => read_pos += 1,
                Step::Deletion(len) | Step::ReferenceSkip(len) => ref_pos += len,
                Step::None => {}
            }
        }

        self.read_length > read_pos && both_n(read_pos, ref_pos, self.read_length - read_pos)
    }
}

impl From<&Feature> for Step {
    fn from(feature: &Feature) -> Self {
        match feature {
            Feature::SoftClip { bases, .. } => Self::SoftClip(bases.len()),
            Feature::Substitution { .. } => Self::Substitution,
            Feature::Deletion { len, .. } => Self::Deletion(*len),
            Feature::Insertion { bases, .. } => Self::Insertion(bases.len()),
            Feature::InsertBase { .. } => Self::InsertBase,
            Feature::Bases { bases, .. } => Self::Bases(bases.len()),
            Feature::ReadBase { .. } => Self::ReadBase,
            Feature::ReferenceSkip { len, .. } => Self::ReferenceSkip(*len),
            Feature::Scores { .. }
            | Feature::QualityScore { .. }
            | Feature::Padding { .. }
            | Feature::HardClip { .. } => Self::None,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{io, num::NonZero};

    use crate::core::Position;
    use crate::fasta;
    use crate::sam::{
        self as sam,
        alignment::{
            RecordBuf,
            record::{
                Flags,
                cigar::{Op, op::Kind},
            },
            record_buf::Cigar,
        },
        header::record::value::{Map, map::ReferenceSequence},
    };

    use super::*;

    const REFERENCE_SEQUENCE: &[u8] = b"ACGTACGTACGTACGTNNNN";
    const XS: Tag = Tag::new(b'X', b'S');

    fn convert(
        start: usize,
        cigar: &[Op],
        sequence: &[u8],
        md: &str,
        nm: u8,
        has_reference: bool,
    ) -> io::Result<Record> {
        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(NonZero::new(REFERENCE_SEQUENCE.len()).unwrap()),
            )
            .build();

        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(start).unwrap())
            .set_cigar(cigar.iter().copied().collect::<Cigar>())
            .set_sequence(sequence.to_vec().into())
            .set_quality_scores(vec![30; sequence.len()].into())
            .set_data(
                [
                    (Tag::EDIT_DISTANCE, Value::from(nm)),
                    (Tag::MISMATCHED_POSITIONS, Value::from(md)),
                    (XS, Value::from(7u8)),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        let sequence = fasta::record::Sequence::from(REFERENCE_SEQUENCE.to_vec());

        let mut record = Record::try_from_alignment_record(&header, &record, |_| {
            Ok(has_reference.then(|| std::sync::Arc::new(sequence)))
        })?;

        record.remove_regenerated_md_nm(REFERENCE_SEQUENCE, REFERENCE_SEQUENCE.len(), false, false);

        Ok(record)
    }

    fn tags(record: &Record) -> Vec<Tag> {
        record.data.iter().map(|(tag, _)| *tag).collect()
    }

    #[test]
    fn test_remove_regenerated_md_nm() -> io::Result<()> {
        // A substitution at the third base and a 2-base deletion.
        let cigar = [
            Op::new(Kind::Match, 4),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Match, 4),
        ];

        let record = convert(1, &cigar, b"ACCTGTAC", "2G1^AC4", 3, true)?;
        assert_eq!(tags(&record), [XS]);

        // Only the matching one is removed.
        let record = convert(1, &cigar, b"ACCTGTAC", "8", 3, true)?;
        assert_eq!(tags(&record), [Tag::MISMATCHED_POSITIONS, XS]);

        let record = convert(1, &cigar, b"ACCTGTAC", "2G1^AC4", 2, true)?;
        assert_eq!(tags(&record), [Tag::EDIT_DISTANCE, XS]);

        // Bases stored verbatim are not compared with a reference.
        let record = convert(1, &cigar, b"ACCTGTAC", "2G1^AC4", 3, false)?;
        assert_eq!(
            tags(&record),
            [Tag::EDIT_DISTANCE, Tag::MISMATCHED_POSITIONS, XS]
        );

        Ok(())
    }

    #[test]
    fn test_remove_regenerated_md_nm_with_reference_ns_and_past_the_end() -> io::Result<()> {
        let cigar = [Op::new(Kind::Match, 4)];

        // Reference positions 17 and 18 are N: they are mismatches unless the read has Ns there.
        let record = convert(15, &cigar, b"GTAA", "2N0N0", 2, true)?;
        assert_eq!(tags(&record), [XS]);

        let record = convert(15, &cigar, b"GTNA", "2N0N0", 2, true)?;
        assert_eq!(
            tags(&record),
            [Tag::EDIT_DISTANCE, Tag::MISMATCHED_POSITIONS, XS]
        );

        let record = convert(18, &cigar, b"NNNN", "3", 0, true)?;
        assert_eq!(
            tags(&record),
            [Tag::EDIT_DISTANCE, Tag::MISMATCHED_POSITIONS, XS]
        );

        Ok(())
    }
}
