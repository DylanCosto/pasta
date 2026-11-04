//! samtools' record sanitizer (`bam_sanitize` in `bam_mate.c`, 1.24), for `fixmate -z` and
//! `checksum -z`/`-a`.

use std::io;

use anyhow::bail;

use crate::bam_raw;

/// Tags removed from unmapped records (`FIX_AUX`).
const UNMAPPED_REMOVED_TAGS: [[u8; 2]; 4] = [*b"NM", *b"MD", *b"CG", *b"SM"];

const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAPPED: u16 = 0x4;

const CIGAR_MATCH: u32 = 0;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;
const CIGAR_SEQUENCE_MATCH: u32 = 7;
const CIGAR_SEQUENCE_MISMATCH: u32 = 8;

/// The sanitizer steps (samtools' `FIX_*` bits).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Sanitize(pub u8);

impl Sanitize {
    pub const POS: u8 = 2;
    pub const MQUAL: u8 = 4;
    pub const UNMAP: u8 = 8;
    pub const CIGAR: u8 = 16;
    pub const AUX: u8 = 32;
    pub const CIGDUP: u8 = 64;
    pub const CIGARX: u8 = 128;
    pub const ALL: u8 = 127;
    pub const ON: u8 = Self::MQUAL | Self::UNMAP | Self::CIGAR | Self::AUX;

    pub fn has(self, flag: u8) -> bool {
        self.0 & flag != 0
    }

    /// `bam_sanitize_options`: samtools' comma-separated keywords.
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        let mut opt = 0;

        for keyword in s.split(',').filter(|k| !k.is_empty()) {
            opt = match keyword {
                k if k.starts_with("all") || k.starts_with('*') => Self::ALL,
                k if k.starts_with("none") || k.starts_with("off") => 0,
                k if k.starts_with("on") => Self::ON,
                k if k.starts_with("pos") => opt | Self::POS,
                k if k.starts_with("mqual") => opt | Self::MQUAL,
                k if k.starts_with("unmap") => opt | Self::UNMAP,
                k if k.starts_with("cigdup") => opt | Self::CIGDUP,
                k if k.starts_with("cigarx") => opt | Self::CIGARX | Self::CIGDUP,
                k if k.starts_with("cigar") => opt | Self::CIGAR,
                k if k.starts_with("aux") => opt | Self::AUX,
                k => bail!("unrecognised sanitize keyword: {k}"),
            };
        }

        Ok(Self(opt))
    }
}

/// `bam_sanitize` on an encoded BAM record (without `block_size`). `lengths` are the reference
/// sequence lengths by ID.
pub fn sanitize_raw(record: &mut Vec<u8>, lengths: &[usize], sanitize: Sanitize) -> io::Result<()> {
    if sanitize.has(Sanitize::POS) && i32_at(record, 0) < 0 {
        set_i32(record, 4, -1);

        if sanitize.has(Sanitize::UNMAP) {
            add_flags(record, FLAG_UNMAPPED);
        }
    }

    if sanitize.has(Sanitize::CIGAR) && flags(record) & FLAG_UNMAPPED == 0 {
        let pos = i64::from(i32_at(record, 4));

        if pos < 0 && sanitize.has(Sanitize::UNMAP) {
            add_flags(record, FLAG_UNMAPPED);
        } else {
            let length = usize::try_from(i32_at(record, 0))
                .ok()
                .and_then(|id| lengths.get(id))
                .map_or(0, |&n| n as i64);

            if pos >= length && sanitize.has(Sanitize::UNMAP) {
                add_flags(record, FLAG_UNMAPPED);

                if sanitize.has(Sanitize::POS) {
                    set_i32(record, 0, -1);
                    set_i32(record, 4, -1);
                }
            } else if end_position(record)? > length {
                trim(record, length)?;
            }
        }
    }

    if flags(record) & FLAG_UNMAPPED != 0 {
        if sanitize.has(Sanitize::CIGAR) && n_cigar(record) > 0 {
            bam_raw::set_cigar(record, &[])?;
        }

        if sanitize.has(Sanitize::MQUAL) {
            record[9] = 0;
        }

        if sanitize.has(Sanitize::AUX) {
            // A long CIGAR's `CG` field is the CIGAR (HTSlib has moved it there), not a tag.
            let long_cigar = match bam_raw::Record::new(record)?.long_cigar_field()? {
                Some(_) => Some(cigar(record)?),
                None => None,
            };

            let mut kept = Vec::with_capacity(record.len());
            bam_raw::strip_fields(record, &mut kept, |tag| {
                UNMAPPED_REMOVED_TAGS.contains(&tag)
            })?;
            *record = kept;

            if let Some(ops) = long_cigar {
                bam_raw::set_cigar(record, &ops)?;
            }
        }
    }

    if sanitize.has(Sanitize::CIGARX) && flags(record) & FLAG_UNMAPPED == 0 {
        let cigar: Vec<u32> = cigar(record)?
            .into_iter()
            .map(|op| match op & 0xf {
                CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => (op & !0xf) | CIGAR_MATCH,
                _ => op,
            })
            .collect();
        bam_raw::set_cigar(record, &cigar)?;
    }

    if sanitize.has(Sanitize::CIGDUP) && flags(record) & FLAG_UNMAPPED == 0 {
        let cigar = merge_cigar_ops(&cigar(record)?);
        bam_raw::set_cigar(record, &cigar)?;
    }

    Ok(())
}

fn i32_at(record: &[u8], at: usize) -> i32 {
    i32::from_le_bytes([record[at], record[at + 1], record[at + 2], record[at + 3]])
}

fn set_i32(record: &mut [u8], at: usize, n: i32) {
    record[at..at + 4].copy_from_slice(&n.to_le_bytes());
}

fn flags(record: &[u8]) -> u16 {
    u16::from_le_bytes([record[14], record[15]])
}

fn add_flags(record: &mut [u8], flags_to_add: u16) {
    let f = flags(record) | flags_to_add;
    record[14..16].copy_from_slice(&f.to_le_bytes());
}

fn n_cigar(record: &[u8]) -> usize {
    usize::from(u16::from_le_bytes([record[12], record[13]]))
}

/// The CIGAR operations, a long CIGAR's from its `CG` field, as HTSlib reads them.
fn cigar(record: &[u8]) -> io::Result<Vec<u32>> {
    let mut ops = Vec::new();
    bam_raw::Record::new(record)?.read_cigar(&mut ops)?;
    Ok(ops)
}

/// `bam_cigar_type` bit 2: whether an operation consumes the reference.
fn consumes_reference(op: u32) -> bool {
    (0x3c1a7 >> ((op & 0xf) << 1)) & 2 != 0
}

/// `bam_endpos`
fn end_position(record: &[u8]) -> io::Result<i64> {
    let pos = i64::from(i32::from_le_bytes([
        record[4], record[5], record[6], record[7],
    ]));
    let flags = u16::from_le_bytes([record[14], record[15]]);

    let len: i64 = if flags & FLAG_UNMAPPED != 0 {
        0
    } else {
        cigar(record)?
            .iter()
            .filter(|&&op| consumes_reference(op))
            .map(|&op| i64::from(op >> 4))
            .sum()
    };

    Ok(pos + len.max(1))
}

/// `bam_trim`: soft clips the alignment past reference position `end`. An alignment entirely
/// past it is unmapped (its CIGAR is then cleared by the unmapped step).
fn trim(record: &mut Vec<u8>, end: i64) -> io::Result<()> {
    let ops = cigar(record)?;
    let mut pos = i64::from(i32::from_le_bytes([
        record[4], record[5], record[6], record[7],
    ]));
    let (mut op, mut len) = (0, 0i64);
    let mut i = 0;

    while i < ops.len() {
        op = ops[i] & 0xf;
        len = i64::from(ops[i] >> 4);

        if consumes_reference(ops[i]) {
            pos += len;

            if pos > end {
                break;
            }
        }

        i += 1;
    }

    if i == ops.len() {
        return Ok(());
    }

    let mut kept = ops[..i].to_vec();
    let mut trailing: Vec<u32> = Vec::new();

    if pos - len < end {
        kept.push((((end - (pos - len)) as u32) << 4) | op);
        trailing.push((((pos - end) as u32) << 4) | CIGAR_SOFT_CLIP);
    } else if pos > end {
        let flags = u16::from_le_bytes([record[14], record[15]]);
        let flags = (flags | FLAG_UNMAPPED) & !FLAG_PROPER_PAIR;
        record[14..16].copy_from_slice(&flags.to_le_bytes());
        return Ok(());
    } else {
        trailing.push(((len as u32) << 4) | CIGAR_SOFT_CLIP);
    }

    for &next in &ops[i + 1..] {
        if next & 0xf == CIGAR_HARD_CLIP {
            trailing.push(next);
        } else if let Some(last) = trailing.last_mut() {
            *last = (((*last >> 4) + (next >> 4)) << 4) | CIGAR_SOFT_CLIP;
        }
    }

    kept.extend(trailing);
    bam_raw::set_cigar(record, &kept)
}

/// `FIX_CIGDUP`: merges adjacent operations of the same kind and drops empty ones.
fn merge_cigar_ops(ops: &[u32]) -> Vec<u32> {
    const MAX_LEN: u32 = (1 << 28) - 1;

    let mut merged: Vec<u32> = Vec::with_capacity(ops.len());

    for &op in ops {
        let len = op >> 4;

        if len == 0 {
            continue;
        }

        match merged.last_mut() {
            Some(last) if *last & 0xf == op & 0xf => {
                let total = u64::from(*last >> 4) + u64::from(len);

                if total > u64::from(MAX_LEN) {
                    *last = (MAX_LEN << 4) | (op & 0xf);
                    merged.push((((total - u64::from(MAX_LEN)) as u32) << 4) | (op & 0xf));
                } else {
                    *last = ((total as u32) << 4) | (op & 0xf);
                }
            }
            _ => merged.push(op),
        }
    }

    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> anyhow::Result<()> {
        assert_eq!(Sanitize::parse("all")?, Sanitize(Sanitize::ALL));
        assert_eq!(
            Sanitize::parse("off,pos,mqual")?,
            Sanitize(Sanitize::POS | Sanitize::MQUAL)
        );
        assert!(Sanitize::parse("bogus").is_err());
        Ok(())
    }

    fn sam_record(line: &str) -> anyhow::Result<Vec<u8>> {
        use pasta::{bam, sam::alignment::io::Write as _};

        let header: pasta::sam::Header = "@SQ\tSN:c\tLN:200000\n".parse()?;
        let record = pasta::sam::Record::try_from(line.as_bytes())?;
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&header, &record)?;
        Ok(writer.get_ref()[4..].to_vec())
    }

    fn read_cigar(record: &[u8]) -> io::Result<Vec<u32>> {
        let mut ops = Vec::new();
        bam_raw::Record::new(record)?.read_cigar(&mut ops)?;
        Ok(ops)
    }

    #[test]
    fn test_sanitize_raw_with_a_long_cigar() -> anyhow::Result<()> {
        const LENGTHS: [usize; 1] = [200_000];

        let canonical = Sanitize(Sanitize::ALL | Sanitize::CIGARX);
        let long = "1=1X".repeat(35_000);
        let seq = "A".repeat(70_000);

        // The `=` and `X` operations in the `CG` field become one M, as in a SAM record.
        let mut record = sam_record(&format!("r\t0\tc\t1\t60\t{long}\t*\t0\t0\t{seq}\t*"))?;
        sanitize_raw(&mut record, &LENGTHS, canonical)?;
        assert_eq!(read_cigar(&record)?, [70_000 << 4]);
        assert!(bam_raw::find_field(&record, *b"CG")?.is_none());

        // A long CIGAR past the reference end is soft clipped there.
        let mut record = sam_record(&format!(
            "r\t0\tc\t150002\t60\t{}\t*\t0\t0\t{seq}\t*",
            "2=2X".repeat(17_500)
        ))?;
        sanitize_raw(&mut record, &LENGTHS, canonical)?;
        assert_eq!(read_cigar(&record)?, [49_999 << 4, (20_001 << 4) | 4]);

        // An unmapped record keeps a long CIGAR when tags are removed but the CIGAR is not.
        let ops: String = (0..35_000).map(|i| ["1M1I", "2M1D"][i % 2]).collect();
        let tags = "NM:i:1\tXA:i:2";
        let mut record = sam_record(&format!(
            "r\t4\tc\t1\t60\t{ops}\t*\t0\t0\t{}\t*\t{tags}",
            "A".repeat(70_000)
        ))?;
        let expected = read_cigar(&record)?;
        sanitize_raw(&mut record, &LENGTHS, Sanitize(Sanitize::AUX))?;
        assert_eq!(read_cigar(&record)?, expected);
        assert!(bam_raw::find_field(&record, *b"NM")?.is_none());
        assert!(bam_raw::find_field(&record, *b"XA")?.is_some());

        Ok(())
    }

    #[test]
    fn test_merge_cigar_ops() {
        let m = |len: u32| len << 4;
        assert_eq!(
            merge_cigar_ops(&[m(10), m(5), (3 << 4) | 1, 0]),
            [m(15), (3 << 4) | 1]
        );
    }
}
