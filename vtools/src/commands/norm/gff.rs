//! The transcripts of a GFF3 file, for `-g/--gff-annot`: their strands and extents, read
//! with the shared GFF reader ([`crate::gff`], bcftools' `gff.c`).

use std::{collections::HashMap, path::Path};

use crate::gff::{Gff, GffOptions, Strand};

/// A transcript's extent (0-based, inclusive), strand, and the largest end of the transcripts
/// sorted before it, for overlap queries.
#[derive(Clone, Copy, Debug)]
struct Span {
    beg: u32,
    end: u32,
    strand: Strand,
    max_end: u32,
}

/// Transcript extents and strands by sequence, for the alignment direction.
pub(super) struct Transcripts {
    by_seq: HashMap<Vec<u8>, Vec<Span>>,
}

impl Transcripts {
    /// Reads the transcripts of a GFF3 file (bcftools' `gff_parse`).
    pub(super) fn parse(path: &Path, verbosity: i32) -> anyhow::Result<Self> {
        let mut gff = Gff::parse(
            path,
            GffOptions {
                verbosity,
                ..GffOptions::default()
            },
        )?;
        let mut by_seq: HashMap<Vec<u8>, Vec<Span>> = HashMap::new();
        let regions: Vec<(Vec<u8>, u32)> = gff
            .idx_tscript
            .iter()
            .map(|r| (r.seq.to_vec(), *r.payload))
            .collect();
        for (seq, tr) in regions {
            let t = &gff.transcripts[tr as usize];
            by_seq.entry(seq).or_default().push(Span {
                beg: t.beg,
                end: t.end,
                strand: t.strand,
                max_end: 0,
            });
        }
        for list in by_seq.values_mut() {
            list.sort_by_key(|t| (t.beg, t.end));
            let mut max_end = 0;
            for t in list.iter_mut() {
                max_end = max_end.max(t.end);
                t.max_end = max_end;
            }
        }
        Ok(Self { by_seq })
    }

    /// Whether to left-align a variant over `beg..=end` (0-based) of a sequence: unless all
    /// transcripts it overlaps are on the forward strand (bcftools' `is_left_align`).
    pub(super) fn left_align(&self, chr: &[u8], beg: i64, end: i64) -> bool {
        let Some(list) = self.by_seq.get(chr) else {
            return true;
        };
        let upto = list.partition_point(|t| i64::from(t.beg) <= end);
        let mut has_fwd = false;
        for t in list[..upto].iter().rev() {
            if i64::from(t.max_end) < beg {
                break;
            }
            if i64::from(t.end) < beg {
                continue;
            }
            match t.strand {
                Strand::Rev => return true,
                Strand::Fwd => has_fwd = true,
                Strand::Unknown => {}
            }
        }
        !has_fwd
    }
}
