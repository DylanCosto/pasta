//! MD and NM tags regenerated from read features, as HTSlib does when decoding.
//!
//! HTSlib's CRAM writer drops MD and NM tags that match the alignment (unless asked to keep
//! them) and its reader rebuilds them for mapped records with a reference, appending MD and then
//! NM after the stored tags. This ports `cram_decode_seq`'s MD/NM generation (HTSlib 1.24), so
//! records decode with the same tags. A `cF` tag (CRAM flags) is removed; its bits 1 and 2 stop
//! MD and NM from being generated.
//!
//! Placement follows HTSlib: a regenerated MD goes after every stored tag, and a regenerated NM
//! after that (so a stored NM stays where it was, before a regenerated MD). The read group, a
//! separate data series, comes last. HTSlib's decoder also accepts `MD*` and `NM*` placeholder
//! tags that mark where regenerated tags go, but no HTSlib release writes them; noodles does not
//! read them.

use std::borrow::Cow;

use crate::core::Position;
use crate::sam::alignment::record::data::field::Tag;

use super::{Feature, Record, data::field::Value};
use crate::cram::io::reader::container::slice::{MdNmTags, ReferenceSequence};

const CRAM_FLAGS: Tag = Tag::new(b'c', b'F');

/// Removes a `cF` tag and adds those of the MD and NM tags HTSlib generates for this record that
/// `md_nm` names.
pub(crate) fn add_missing_md_nm(record: &mut Record<'_>, md_nm: MdNmTags) {
    let mut cram_flags = 0;

    if let Some(i) = record.data.iter().position(|(tag, _)| *tag == CRAM_FLAGS)
        && let Value::UInt8(flags) = record.data[i].1
    {
        cram_flags = flags;
        record.data.remove(i);
    }

    if !md_nm.any() || record.bam_flags.is_unmapped() || record.cram_flags.sequence_is_missing() {
        return;
    }

    let has_md = !md_nm.md
        || cram_flags & 1 != 0
        || record
            .data
            .iter()
            .any(|(tag, _)| *tag == Tag::MISMATCHED_POSITIONS);
    let has_nm = !md_nm.nm
        || cram_flags & 2 != 0
        || record
            .data
            .iter()
            .any(|(tag, _)| *tag == Tag::EDIT_DISTANCE);

    if has_md && has_nm {
        return;
    }

    let Some((md, nm)) = calculate(record, !has_md) else {
        return;
    };

    if !has_md {
        record.data.push((
            Tag::MISMATCHED_POSITIONS,
            Value::String(Cow::Owned(md.into())),
        ));
    }

    if !has_nm {
        let value = if let Ok(n) = u8::try_from(nm) {
            Value::UInt8(n)
        } else if let Ok(n) = u16::try_from(nm) {
            Value::UInt16(n)
        } else {
            Value::UInt32(nm)
        };

        record.data.push((Tag::EDIT_DISTANCE, value));
    }
}

/// A record's reference bases: the loaded bases, the 0-based reference position of the first,
/// and the reference sequence's length.
pub(crate) struct Reference<'a> {
    pub(crate) bases: &'a [u8],
    pub(crate) offset: usize,
    pub(crate) len: usize,
}

impl Reference<'_> {
    /// The uppercased base at a 0-based reference position, or `None` outside the loaded bases.
    fn get(&self, pos: usize) -> Option<u8> {
        pos.checked_sub(self.offset)
            .and_then(|i| self.bases.get(i))
            .map(|b| b & !0x20)
    }
}

fn reference<'a>(record: &'a Record<'_>) -> Option<Reference<'a>> {
    let reference_sequence_id = record.reference_sequence_id?;
    let len = record
        .header?
        .reference_sequences()
        .get_index(reference_sequence_id)
        .map(|(_, rs)| usize::from(rs.length()))?;

    let (bases, offset) = match record.reference_sequence.as_ref()? {
        ReferenceSequence::Embedded {
            reference_start,
            sequence,
        } => (*sequence, usize::from(*reference_start) - 1),
        ReferenceSequence::ExternalRange {
            reference_start,
            sequence,
        } => ((**sequence).as_ref(), usize::from(*reference_start) - 1),
        ReferenceSequence::External { sequence } => ((**sequence).as_ref(), 0),
    };

    Some(Reference { bases, offset, len })
}

/// The MD string being built, or nothing when only NM is needed.
struct Md<'a> {
    buf: Option<&'a mut Vec<u8>>,
}

impl Md<'_> {
    fn push(&mut self, b: u8) {
        if let Some(buf) = &mut self.buf {
            buf.push(b);
        }
    }

    fn push_number(&mut self, n: i64) {
        use lexical_core::FormattedSize;

        if let Some(buf) = &mut self.buf {
            let mut dst = [0; i64::FORMATTED_SIZE_DECIMAL];

            // Most are a digit or two, pushed rather than copied.
            for &b in lexical_core::write(n, &mut dst).iter() {
                buf.push(b);
            }
        }
    }
}

/// What a read feature contributes to MD and NM, with the feature's read position.
pub(crate) enum Step {
    SoftClip(usize),
    Substitution,
    Deletion(usize),
    Insertion(usize),
    InsertBase,
    Bases(usize),
    ReadBase,
    ReferenceSkip(usize),
    /// Scores, padding and hard clips.
    None,
}

impl From<&Feature<'_>> for Step {
    fn from(feature: &Feature<'_>) -> Self {
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

/// The MD string (empty unless `with_md`) and NM count, or `None` if the record has no reference
/// or reads outside the loaded reference bases (where HTSlib fails to decode the slice).
fn calculate(record: &Record<'_>, with_md: bool) -> Option<(Vec<u8>, u32)> {
    let reference = reference(record)?;
    let steps = record
        .features
        .iter()
        .map(|feature| (usize::from(feature.position()), Step::from(feature)));

    let mut md = Vec::new();

    let nm = calculate_from_steps(
        &reference,
        record.alignment_start?,
        record.read_length,
        steps,
        with_md.then_some(&mut md),
    )?;

    Some((md, nm))
}

/// The NM count of a record from its features' `steps`, writing its MD string to `md`, if any
/// (cleared first), or `None` if it reads outside the loaded reference bases.
pub(crate) fn calculate_from_steps<I>(
    reference: &Reference<'_>,
    alignment_start: Position,
    read_length: usize,
    steps: I,
    md: Option<&mut Vec<u8>>,
) -> Option<u32>
where
    I: IntoIterator<Item = (usize, Step)>,
{
    let mut md = Md { buf: md };

    if let Some(buf) = &mut md.buf {
        buf.clear();
    }
    let mut md_dist: i64 = 0;
    let mut nm: u32 = 0;
    let mut seq_pos = 1;
    let mut ref_pos = usize::from(alignment_start) - 1;

    // `add_md_char`
    let push_mismatch = |md: &mut Md, md_dist: &mut i64, base: u8| {
        md.push_number(*md_dist);
        md.push(base);
        *md_dist = 0;
    };

    // Matching bases, where reference Ns are mismatches.
    let matches = |md: &mut Md, md_dist: &mut i64, nm: &mut u32, start: usize, len: usize| {
        if start + len > reference.len {
            if *md_dist >= 0 {
                *md_dist += len as i64;
            }

            return Some(());
        }

        // Most spans hold no N: they are all matches.
        let i = start.checked_sub(reference.offset)?;
        let bases = reference.bases.get(i..i + len)?;

        if memchr::memchr2(b'N', b'n', bases).is_none() {
            *md_dist += len as i64;
            return Some(());
        }

        for i in 0..len {
            if reference.get(start + i)? == b'N' {
                push_mismatch(md, md_dist, b'N');
                *nm += 1;
            } else {
                *md_dist += 1;
            }
        }

        Some(())
    };

    let push_dist = |md: &mut Md, md_dist: i64| {
        if md_dist >= 0 {
            md.push_number(md_dist);
        }
    };

    for (pos, step) in steps {
        if pos > seq_pos {
            let len = pos - seq_pos;
            matches(&mut md, &mut md_dist, &mut nm, ref_pos, len)?;
            ref_pos += len;
            seq_pos = pos;
        }

        match step {
            Step::SoftClip(len) => seq_pos += len,
            Step::Substitution => {
                if ref_pos >= reference.len {
                    push_dist(&mut md, md_dist);
                    md_dist = -1;
                    nm = nm.wrapping_sub(1);
                } else {
                    let base = reference.get(ref_pos)?;
                    push_mismatch(&mut md, &mut md_dist, base);
                }

                nm = nm.wrapping_add(1);
                seq_pos += 1;
                ref_pos += 1;
            }
            Step::Deletion(len) => {
                push_dist(&mut md, md_dist);

                if ref_pos + len <= reference.len {
                    md.push(b'^');

                    for i in 0..len {
                        md.push(reference.get(ref_pos + i)?);
                    }

                    md_dist = 0;
                    nm += len as u32;
                } else {
                    if reference.len >= ref_pos {
                        md.push(b'^');

                        for i in ref_pos..reference.len {
                            md.push(reference.get(i)?);
                        }

                        md.push(b'0');
                        nm += (reference.len - ref_pos) as u32;
                    }

                    md_dist = -1;
                }

                ref_pos += len;
            }
            Step::Insertion(len) => {
                seq_pos += len;
                nm += len as u32;
            }
            Step::InsertBase => {
                seq_pos += 1;
                nm += 1;
            }
            Step::Bases(len) => {
                push_dist(&mut md, md_dist);

                let mut x = 0;

                while x < len {
                    if x > 0 {
                        md.push(b'0');
                    }

                    if ref_pos + x >= reference.len {
                        break;
                    }

                    md.push(reference.get(ref_pos + x)?);
                    x += 1;
                }

                nm += x as u32;
                md_dist = 0;
                seq_pos += len;
                ref_pos += len;
            }
            Step::ReadBase => {
                push_dist(&mut md, md_dist);

                if ref_pos >= reference.len {
                    md_dist = -1;
                } else {
                    md.push(reference.get(ref_pos)?);
                    nm += 1;
                    md_dist = 0;
                }

                seq_pos += 1;
                ref_pos += 1;
            }
            Step::ReferenceSkip(len) => ref_pos += len,
            Step::None => {}
        }
    }

    if read_length >= seq_pos {
        let len = read_length - seq_pos + 1;
        matches(&mut md, &mut md_dist, &mut nm, ref_pos, len)?;
    }

    push_dist(&mut md, md_dist);

    Some(nm)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::core::Position;
    use crate::fasta;
    use crate::sam::{self as sam, header::record::value::map::ReferenceSequence as Sq};

    use bstr::BStr;

    use super::*;

    fn header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                sam::header::record::value::Map::<Sq>::new(
                    std::num::NonZero::new(12).expect("non-zero"),
                ),
            )
            .build()
    }

    fn record<'c>(header: &'c sam::Header, features: Vec<Feature<'c>>) -> Record<'c> {
        Record {
            header: Some(header),
            bam_flags: sam::alignment::record::Flags::empty(),
            reference_sequence_id: Some(0),
            alignment_start: Position::new(2),
            read_length: 8,
            reference_sequence: Some(ReferenceSequence::External {
                sequence: Arc::new(fasta::record::Sequence::from(b"ACGTNacgtACG".to_vec())),
            }),
            features,
            ..Record::default()
        }
    }

    fn md_nm(record: &Record<'_>) -> (String, u32) {
        let (md, nm) = calculate(record, true).expect("calculated");
        (String::from_utf8(md).expect("ASCII"), nm)
    }

    #[test]
    fn test_calculate() {
        let header = header();
        let position = |n| Position::new(n).expect("non-zero");

        // Reference CGTNacgt: the N is a mismatch.
        assert_eq!(md_nm(&record(&header, Vec::new())), ("3N4".into(), 1));

        let features = vec![
            Feature::Substitution {
                position: position(2),
                code: 0,
            },
            Feature::Deletion {
                position: position(4),
                len: 2,
            },
            Feature::InsertBase {
                position: position(6),
                base: b'A',
            },
        ];

        // C, G>?, T, ^NA, then A (inserted), CGTA...
        let r = record(&header, features);
        assert_eq!(md_nm(&r), ("1G1^NA4".into(), 4));

        // Without MD, only NM is counted.
        assert_eq!(calculate(&r, false), Some((Vec::new(), 4)));
    }

    #[test]
    fn test_add_missing_md_nm() {
        let header = header();
        let mut r = record(&header, Vec::new());
        r.data.push((CRAM_FLAGS, Value::UInt8(2)));
        add_missing_md_nm(&mut r, MdNmTags::ALL);

        assert_eq!(r.data.len(), 1);
        assert_eq!(r.data[0].0, Tag::MISMATCHED_POSITIONS);
        assert!(matches!(&r.data[0].1, Value::String(s) if s.as_ref() == BStr::new("3N4")));

        let mut r = record(&header, Vec::new());
        r.data.push((CRAM_FLAGS, Value::UInt8(2)));
        add_missing_md_nm(&mut r, MdNmTags::default());
        assert!(r.data.is_empty());

        // Only NM, without MD.
        let mut r = record(&header, Vec::new());
        add_missing_md_nm(
            &mut r,
            MdNmTags {
                md: false,
                nm: true,
            },
        );
        assert_eq!(r.data.len(), 1);
        assert_eq!(r.data[0].0, Tag::EDIT_DISTANCE);
        assert!(matches!(r.data[0].1, Value::UInt8(1)));
    }
}
