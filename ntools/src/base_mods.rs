//! Base modification tags (`MM`, `ML`, and `MN`; SAM tags specification § 1.7) of reads whose
//! bases are removed, such as soft clips.

use std::{io, sync::Once};

use crate::bam_raw;

const MM: [u8; 2] = *b"MM";
const ML: [u8; 2] = *b"ML";
const MN: [u8; 2] = *b"MN";

/// Rewrites the `MM`, `ML`, and `MN` fields of the encoded BAM record `record`, whose sequence is
/// `original` with only the bases `kept` selects remaining, to describe the remaining bases.
/// Calls on removed bases are dropped with their `ML` probabilities, and the skip counts are
/// recounted over the remaining bases. `original` and `kept` are in the read's original
/// orientation (reverse complemented for reverse-strand reads), the orientation `MM` counts in.
///
/// Tags that do not describe `original` (malformed, calls past its end, an `ML` of another
/// length, or an `MN` other than its length) cannot be remapped, and are removed with a warning.
///
/// # Errors
///
/// Returns an error if the record's auxiliary data is malformed.
pub fn trim(record: &mut Vec<u8>, original: &[u8], kept: &[bool]) -> io::Result<()> {
    let [mm_range, ml_range, mn_range] = bam_raw::find_fields(record, [MM, ML, MN])?;

    let Some(mm_range) = mm_range else {
        return Ok(());
    };

    let mm =
        (record[mm_range.start + 2] == b'Z').then(|| &record[mm_range.start + 3..mm_range.end - 1]);

    // `ML` is a `B:C` array.
    let ml = match &ml_range {
        Some(range) if record[range.start + 2..range.start + 4] == *b"BC" => {
            Some(Some(&record[range.start + 8..range.end]))
        }
        Some(_) => None,
        None => Some(None),
    };

    let mn_matches = mn_range
        .as_ref()
        .is_none_or(|range| integer_value(&record[range.clone()]) == Some(original.len() as i64));

    let remapped = match (mm, ml) {
        (Some(mm), Some(ml)) if mn_matches => remap(mm, ml, original, kept),
        _ => None,
    };

    let Some((new_mm, new_ml)) = remapped else {
        static WARNING: Once = Once::new();
        WARNING.call_once(|| {
            eprintln!(
                "[W::remove_soft_clips] Removing MM/ML/MN tags that do not match the read sequence"
            );
        });

        let mut dst = Vec::with_capacity(record.len());
        bam_raw::strip_fields(record, &mut dst, |tag| [MM, ML, MN].contains(&tag))?;
        *record = dst;
        return Ok(());
    };

    let new_len = kept.iter().filter(|&&k| k).count();

    let mut replacements = vec![(mm_range, string_field(MM, &new_mm))];

    if let (Some(range), Some(ml)) = (ml_range, new_ml) {
        replacements.push((range, byte_array_field(ML, &ml)));
    }

    if let Some(range) = mn_range {
        let mut field = MN.to_vec();
        field.push(b'i');
        let n =
            i32::try_from(new_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        field.extend_from_slice(&n.to_le_bytes());
        replacements.push((range, field));
    }

    // From the last field back, so the earlier ranges stay valid.
    replacements.sort_by_key(|(range, _)| std::cmp::Reverse(range.start));

    for (range, field) in replacements {
        record.splice(range, field);
    }

    Ok(())
}

/// Returns the value of an encoded integer field (its tag, type, and value).
fn integer_value(field: &[u8]) -> Option<i64> {
    let value = field.get(3..)?;

    Some(match field[2] {
        b'c' => i64::from(i8::from_le_bytes([*value.first()?])),
        b'C' => i64::from(*value.first()?),
        b's' => i64::from(i16::from_le_bytes(value.get(..2)?.try_into().ok()?)),
        b'S' => i64::from(u16::from_le_bytes(value.get(..2)?.try_into().ok()?)),
        b'i' => i64::from(i32::from_le_bytes(value.get(..4)?.try_into().ok()?)),
        b'I' => i64::from(u32::from_le_bytes(value.get(..4)?.try_into().ok()?)),
        _ => return None,
    })
}

fn string_field(tag: [u8; 2], value: &[u8]) -> Vec<u8> {
    let mut field = Vec::with_capacity(value.len() + 4);
    field.extend_from_slice(&tag);
    field.push(b'Z');
    field.extend_from_slice(value);
    field.push(0);
    field
}

fn byte_array_field(tag: [u8; 2], values: &[u8]) -> Vec<u8> {
    let count = u32::try_from(values.len()).expect("no more values than before");
    let mut field = Vec::with_capacity(values.len() + 8);
    field.extend_from_slice(&tag);
    field.extend_from_slice(b"BC");
    field.extend_from_slice(&count.to_le_bytes());
    field.extend_from_slice(values);
    field
}

/// Remaps an `MM` value and its `ML` probabilities (if any) from `original` to the bases of it
/// that `kept` selects, both in the original orientation. Returns `None` if they do not describe
/// `original`.
fn remap(
    mm: &[u8],
    ml: Option<&[u8]>,
    original: &[u8],
    kept: &[bool],
) -> Option<(Vec<u8>, Option<Vec<u8>>)> {
    let mut new_mm = Vec::with_capacity(mm.len());
    let mut new_ml = ml.map(|ml| Vec::with_capacity(ml.len()));
    let mut ml_pos = 0;

    let mut groups = mm.split(|&b| b == b';');
    // A value ends with `;`, so the last piece is empty.
    let last = groups.next_back()?;

    if !last.is_empty() {
        return None;
    }

    for group in groups {
        let mut parts = group.split(|&b| b == b',');
        let head = parts.next()?;
        let deltas: Vec<u64> = parts
            .map(|delta| std::str::from_utf8(delta).ok()?.parse().ok())
            .collect::<Option<_>>()?;

        // `C+m`, `C+mh.`, `N-12345?`: the base, the strand, and the modification codes, with an
        // optional implicit (`.`) or unknown (`?`) mode.
        let (&base, rest) = head.split_first()?;
        let (&strand, codes) = rest.split_first()?;
        let codes = codes
            .strip_suffix(b".")
            .or_else(|| codes.strip_suffix(b"?"))
            .unwrap_or(codes);

        if !b"ACGTUN".contains(&base) || !b"+-".contains(&strand) || codes.is_empty() {
            return None;
        }

        let code_count = if codes.iter().all(u8::is_ascii_digit) {
            1
        } else if codes.iter().all(u8::is_ascii_lowercase) {
            codes.len()
        } else {
            return None;
        };

        let matches = |b: u8| match base {
            b'N' => true,
            b'U' => b == b'T',
            _ => b == base,
        };

        new_mm.extend_from_slice(head);

        let mut calls = deltas.iter();
        let mut skip = calls.next().copied();
        // Remaining matching bases since the last remaining call.
        let mut pending = 0;
        let mut call = 0;

        for (&b, &is_kept) in original.iter().zip(kept) {
            let Some(n) = skip else {
                break;
            };

            if !matches(b.to_ascii_uppercase()) {
                continue;
            }

            if n > 0 {
                skip = Some(n - 1);
                pending += u64::from(is_kept);
                continue;
            }

            if is_kept {
                new_mm.push(b',');
                new_mm.extend_from_slice(pending.to_string().as_bytes());
                pending = 0;

                if let (Some(ml), Some(new_ml)) = (ml, new_ml.as_mut()) {
                    let start = ml_pos + call * code_count;
                    new_ml.extend_from_slice(ml.get(start..start + code_count)?);
                }
            }

            call += 1;
            skip = calls.next().copied();
        }

        // A call past the end of the sequence.
        if skip.is_some() {
            return None;
        }

        new_mm.push(b';');
        ml_pos += deltas.len() * code_count;
    }

    if ml.is_some_and(|ml| ml.len() != ml_pos) {
        return None;
    }

    Some((new_mm, new_ml))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mask(s: &str) -> Vec<bool> {
        s.bytes().map(|b| b == b'1').collect()
    }

    #[test]
    fn test_remap() {
        // `CACGCC` with `2S4M` keeps `CGCC`; the first call (base 0) is clipped, and the
        // second (base 2) becomes the first C of the remaining bases.
        assert_eq!(
            remap(b"C+m,0,0;", Some(&[10, 200]), b"CACGCC", &mask("001111")),
            Some((b"C+m,0;".to_vec(), Some(vec![200])))
        );

        // Skipped Cs in the clip are not counted; two codes per call; several groups; a group
        // without calls; the unknown mode.
        assert_eq!(
            remap(
                b"C+mh?,1,1;A+a,0;G-m;",
                Some(&[1, 2, 3, 4, 5]),
                b"CCACCA",
                &mask("011110"),
            ),
            Some((b"C+mh?,0,1;A+a,0;G-m;".to_vec(), Some(vec![1, 2, 3, 4, 5])))
        );

        assert_eq!(
            remap(b"N+n,1;", None, b"ACGT", &mask("0111")),
            Some((b"N+n,0;".to_vec(), None))
        );

        // Calls past the end, and an `ML` of the wrong length.
        assert_eq!(remap(b"C+m,4;", None, b"CACGCC", &mask("111111")), None);
        assert_eq!(remap(b"C+m,0;", Some(&[1, 2]), b"C", &mask("1")), None);
        assert_eq!(remap(b"C+m,0", None, b"C", &mask("1")), None);
    }
}
