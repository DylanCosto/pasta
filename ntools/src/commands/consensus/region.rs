//! Region strings (`hts_parse_region`) and BED region lists (`bed_read`, `bed_reglist`) as
//! samtools consensus reads them.

use std::{collections::HashMap, path::Path};

use anyhow::{Context, bail};

use crate::{hts_region, khash::StringSet, regions};

/// `hts_parse_region` over the header's names (see [`hts_region::parse_region`]).
pub fn parse_region(s: &str, names: &[Vec<u8>]) -> Option<(usize, i64, i64)> {
    hts_region::parse_region(s, |name| names.iter().position(|n| n.as_slice() == name))
}

/// A `bed_read` region list: intervals by name, with names in khash order and each name's
/// intervals sorted (not merged).
#[derive(Debug, Default)]
pub struct Bed {
    order: StringSet,
    intervals: HashMap<Vec<u8>, Vec<(i64, i64)>>,
}

impl Bed {
    /// `bed_read` (see [`regions::read_bed`]).
    pub fn read(path: &Path) -> anyhow::Result<Self> {
        let mut bed = Self::default();

        for interval in regions::read_bed(path)? {
            bed.order.insert(&interval.name);
            bed.intervals
                .entry(interval.name)
                .or_default()
                .push((interval.start as i64, interval.end as i64));
        }

        for intervals in bed.intervals.values_mut() {
            intervals.sort_unstable();
        }

        Ok(bed)
    }

    /// `bed_reglist` and `fill_reglist_tid`: the intervals of the names in the header, in
    /// khash order, with their reference sequence IDs. Unknown names are skipped with a warning.
    pub fn reglist(&self, names: &[Vec<u8>]) -> anyhow::Result<Vec<(usize, i64, i64)>> {
        if self.order.is_empty() {
            bail!("the regions file has no regions");
        }

        let mut list = Vec::new();
        let mut known = 0;

        for name in self.order.iter() {
            match names.iter().position(|n| n.as_slice() == name) {
                Some(tid) => {
                    known += 1;
                    let intervals = self.intervals.get(name).map_or(&[][..], Vec::as_slice);
                    list.extend(intervals.iter().map(|&(beg, end)| (tid, beg, end)));
                }
                None => eprintln!(
                    "[W::fill_reglist_tid] Region '{}' specifies an unknown reference name",
                    String::from_utf8_lossy(name)
                ),
            }
        }

        if known == 0 {
            bail!("no regions file region names are in the header");
        }

        Ok(list)
    }
}

/// Reads a BED file, with the error message samtools gives.
pub fn read_bed(path: &Path) -> anyhow::Result<Bed> {
    Bed::read(path).with_context(|| format!("Could not read file \"{}\"", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hts_region::POS_MAX;

    fn names() -> Vec<Vec<u8>> {
        vec![b"chr1".to_vec(), b"chr2".to_vec(), b"HLA:1".to_vec()]
    }

    #[test]
    fn test_parse_region() {
        let names = names();

        assert_eq!(parse_region("chr1", &names), Some((0, 0, POS_MAX)));
        assert_eq!(parse_region("chr2:100", &names), Some((1, 99, POS_MAX)));
        assert_eq!(parse_region("chr2:1,000-2k", &names), Some((1, 999, 2000)));
        assert_eq!(parse_region("chr1:-100", &names), Some((0, 0, 100)));
        assert_eq!(parse_region("chr1:100-", &names), Some((0, 99, POS_MAX)));
        assert_eq!(parse_region("HLA:1", &names), Some((2, 0, POS_MAX)));
        assert_eq!(parse_region("{HLA:1}:5-6", &names), Some((2, 4, 6)));
        assert_eq!(parse_region("chr1:200-100", &names), None);
        assert_eq!(parse_region("chr1:1-2x", &names), None);
        assert_eq!(parse_region("chr3", &names), None);
    }
}
