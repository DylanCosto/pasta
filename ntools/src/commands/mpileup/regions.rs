//! Positions and regions to pile up: BED and position files (`-l`) and `-r` regions.

use std::path::Path;

use anyhow::Context;

use crate::{hts_region, regions};

/// `-l` intervals by sequence name: 0-based, half-open, sorted by start (`bed_read`), with the
/// running maximum of their ends.
#[derive(Debug, Default)]
pub(super) struct Bed {
    pub(super) by_name: std::collections::HashMap<Vec<u8>, Intervals>,
}

/// Intervals sorted by start, and the running maximum of their ends.
pub(super) type Intervals = (Vec<(i64, i64)>, Vec<i64>);

impl Bed {
    /// `bed_read` (see [`regions::read_bed`]).
    pub(super) fn read(path: &Path) -> anyhow::Result<Self> {
        let intervals = regions::read_bed(path)
            .with_context(|| format!("Could not read file \"{}\"", path.display()))?;
        let mut bed = Self::default();

        for interval in intervals {
            bed.by_name
                .entry(interval.name)
                .or_default()
                .0
                .push((interval.start as i64, interval.end as i64));
        }

        bed.index();
        Ok(bed)
    }

    /// Sorts the intervals and computes the running maximum of their ends.
    pub(super) fn index(&mut self) {
        for (intervals, max_ends) in self.by_name.values_mut() {
            intervals.sort_unstable();

            *max_ends = intervals
                .iter()
                .scan(i64::MIN, |max, &(_, end)| {
                    *max = (*max).max(end);
                    Some(*max)
                })
                .collect();
        }
    }

    /// `bed_overlap`: whether `[beg, end)` overlaps an interval on `name`.
    pub(super) fn overlaps(&self, name: &[u8], beg: i64, end: i64) -> bool {
        let Some((intervals, max_ends)) = self.by_name.get(name) else {
            return false;
        };

        // Of the intervals starting before `end`, whether any ends after `beg`.
        let n = intervals.partition_point(|&(b, _)| b < end);
        n > 0 && max_ends[n - 1] > beg
    }
}

/// Parses `-r` as `sam_itr_querys` does (`hts_parse_region`): 0-based start and exclusive end.
pub(super) fn parse_region(raw: &str, names: &[Vec<u8>]) -> anyhow::Result<(i32, i64, i64)> {
    let (tid, beg, end) =
        hts_region::parse_region(raw, |name| names.iter().position(|n| n.as_slice() == name))
            .with_context(|| format!("[E::mpileup] fail to parse region '{raw}'"))?;
    let end = if end == hts_region::POS_MAX {
        i64::MAX
    } else {
        end
    };
    Ok((i32::try_from(tid)?, beg, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed_overlaps() {
        let mut bed = Bed::default();
        bed.by_name.insert(
            b"chr1".to_vec(),
            (vec![(200, 300), (10, 20), (0, 100_000)], Vec::new()),
        );
        bed.index();

        assert!(bed.overlaps(b"chr1", 15, 16));
        assert!(bed.overlaps(b"chr1", 99_999, 100_000));
        assert!(!bed.overlaps(b"chr1", 100_000, 100_001));
        assert!(!bed.overlaps(b"chr2", 15, 16));
    }
}
