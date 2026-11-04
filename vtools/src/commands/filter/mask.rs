//! Regions for `--mask` and `--mask-file` (bcftools' `regidx` with its line parsers).

use std::{collections::HashMap, path::Path};

use anyhow::bail;
use pasta::variant::{index::split_idx_name, number, readlist};

/// The largest coordinate (0-based) a region may have (bcftools' `REGIDX_MAX`).
const MAX_COORD: i64 = 2_147_483_646;

/// How region lines are read (bcftools' `regidx_parse_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum LineFormat {
    /// `CHR:BEG-END`, 1-based (`regidx_parse_reg`).
    Reg,
    /// `CHR BEG END`, 1-based (`regidx_parse_tab`).
    Tab,
    /// `CHR BEG END`, 0-based half-open (`regidx_parse_bed`).
    Bed,
    /// `CHR POS ...`, a single position (`regidx_parse_vcf`).
    Vcf,
}

/// Sorted regions per sequence, queried for overlaps (bcftools' `regidx_t`).
#[derive(Debug, Default)]
pub struct Mask {
    seqs: HashMap<Vec<u8>, Vec<(i64, i64)>>,
}

/// The outcome of parsing one line.
enum Parsed<'a> {
    Skip,
    Region(&'a [u8], i64, i64),
    Error,
}

#[inline]
fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

fn skip_space(s: &[u8], mut i: usize) -> usize {
    while at(s, i).is_ascii_whitespace() {
        i += 1;
    }
    i
}

/// C's `(uint32_t) strtod(...)`.
fn coord(d: f64) -> i64 {
    if d <= 0.0 {
        0
    } else {
        d.min(u32::MAX as f64) as i64
    }
}

fn parse_line(line: &[u8], format: LineFormat) -> Parsed<'_> {
    let ss = skip_space(line, 0);
    if ss >= line.len() || line[ss] == b'#' {
        return Parsed::Skip;
    }
    let name_end = if format == LineFormat::Reg {
        ss + line[ss..]
            .iter()
            .position(|&c| c == b':')
            .unwrap_or(line.len() - ss)
    } else {
        ss + line[ss..]
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(line.len() - ss)
    };
    let chr = &line[ss..name_end];
    if name_end >= line.len() {
        return Parsed::Region(chr, 0, MAX_COORD);
    }
    let text = String::from_utf8_lossy(line);
    let start = name_end + 1;
    let (b, se) = number::strtod(line, start);
    if se == start {
        let kind = match format {
            LineFormat::Reg => "reg",
            LineFormat::Bed => "bed",
            LineFormat::Tab | LineFormat::Vcf => "tab",
        };
        eprintln!("Could not parse {kind} line: {text}");
        return Parsed::Error;
    }
    let beg = coord(b);
    if format == LineFormat::Bed {
        let s2 = se + 1;
        let (e, se2) = number::strtod(line, s2);
        if se2 == s2 {
            eprintln!("Could not parse bed line: {text}");
            return Parsed::Error;
        }
        return Parsed::Region(chr, beg, coord(e) - 1);
    }
    if beg == 0 {
        let kind = if format == LineFormat::Reg {
            "reg"
        } else {
            "tab"
        };
        eprintln!("Could not parse {kind} line, expected 1-based coordinate: {text}");
        return Parsed::Error;
    }
    let beg = beg - 1;
    if format == LineFormat::Vcf {
        return Parsed::Region(chr, beg, beg);
    }
    let end = if at(line, se) == 0 || at(line, se + 1) == 0 {
        if format == LineFormat::Reg && at(line, se) == b'-' {
            MAX_COORD
        } else {
            beg
        }
    } else {
        let s2 = se + 1;
        let (e, se2) = number::strtod(line, s2);
        let bad_tail =
            format == LineFormat::Tab && at(line, se2) != 0 && !at(line, se2).is_ascii_whitespace();
        if se2 == s2 || bad_tail {
            beg
        } else if coord(e) == 0 {
            let kind = if format == LineFormat::Reg {
                "reg"
            } else {
                "tab"
            };
            eprintln!("Could not parse {kind} line, expected 1-based coordinate: {text}");
            return Parsed::Error;
        } else {
            coord(e) - 1
        }
    };
    Parsed::Region(chr, beg, end)
}

impl Mask {
    fn push(&mut self, chr: &[u8], beg: i64, end: i64) {
        let (mut beg, mut end) = (beg.min(MAX_COORD), end.min(MAX_COORD));
        if beg > end {
            std::mem::swap(&mut beg, &mut end);
        }
        self.seqs.entry(chr.to_vec()).or_default().push((beg, end));
    }

    fn finish(mut self) -> Self {
        for regs in self.seqs.values_mut() {
            regs.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
            // Keep each region's end as the running maximum, so a binary search on the start
            // finds every overlap.
            let mut max_end = i64::MIN;
            for r in regs.iter_mut() {
                max_end = max_end.max(r.1);
                r.1 = max_end;
            }
        }
        self
    }

    /// Reads comma-separated `CHR[:BEG[-END]]` regions; malformed ones are reported and skipped
    /// (bcftools' `regidx_init_string` with `regidx_parse_reg`).
    pub fn from_regions(list: &[u8]) -> Self {
        let mut mask = Self::default();
        for line in list.split(|&c| c == b',' || c == b'\n' || c == b'\r') {
            if let Parsed::Region(chr, beg, end) = parse_line(line, LineFormat::Reg) {
                mask.push(chr, beg, end);
            }
        }
        mask.finish()
    }

    /// Reads a regions file, as BED for `.bed[.gz]`, a VCF for `.vcf[.gz]`, else tab-delimited
    /// (bcftools' `regidx_init`). `FILE##idx##INDEX` reads `FILE`, as `hts_open` does, and
    /// its format follows the name of `FILE`, where bcftools looks at the whole name.
    ///
    /// # Errors
    ///
    /// Fails when the file cannot be read or a line cannot be parsed.
    pub fn from_file(path: &str) -> anyhow::Result<Self> {
        let data = split_idx_name(path).0;
        let lower = data.to_ascii_lowercase();
        let format =
            if lower.ends_with(".bed.gz") || lower.ends_with(".bed.bgz") || lower.ends_with(".bed")
            {
                LineFormat::Bed
            } else if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
                LineFormat::Vcf
            } else {
                LineFormat::Tab
            };
        let lines = readlist::read_all_lines(Path::new(data))
            .inspect_err(|e| crate::output::log_open_failure(data, e))?;
        let mut mask = Self::default();
        for line in &lines {
            match parse_line(line, format) {
                Parsed::Skip => {}
                Parsed::Region(chr, beg, end) => mask.push(chr, beg, end),
                Parsed::Error => bail!("could not parse {path}"),
            }
        }
        Ok(mask.finish())
    }

    /// Whether any region of `chr` overlaps `[beg, end]` (0-based, inclusive)
    /// (`regidx_overlap`).
    pub fn overlaps(&self, chr: &[u8], beg: i64, end: i64) -> bool {
        let Some(regs) = self.seqs.get(chr) else {
            return false;
        };
        // Regions starting at or before `end`; the last of them has the largest end so far.
        let n = regs.partition_point(|r| r.0 <= end);
        n > 0 && regs[n - 1].1 >= beg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_strings() {
        let m = Mask::from_regions(b"chr1:100-200,chr2,chr3:50,chr4:10-");
        assert!(m.overlaps(b"chr1", 99, 99));
        assert!(!m.overlaps(b"chr1", 98, 98));
        assert!(m.overlaps(b"chr1", 150, 300));
        assert!(!m.overlaps(b"chr1", 200, 300));
        assert!(m.overlaps(b"chr2", 1_000_000, 1_000_000));
        assert!(m.overlaps(b"chr3", 49, 49));
        assert!(!m.overlaps(b"chr3", 50, 50));
        assert!(m.overlaps(b"chr4", 1 << 30, 1 << 30));
        assert!(!m.overlaps(b"chr5", 1, 1));
    }

    #[test]
    fn idx_names_read_the_data_file_by_its_name() {
        let (dir, name) =
            crate::raw_input::write_idx_named("mask-idx", "m.bed", "chr1\t100\t200\n");
        let m = Mask::from_file(&name);
        let _ = std::fs::remove_dir_all(&dir);
        let m = m.expect("the data file is read");
        assert!(m.overlaps(b"chr1", 100, 100));
        assert!(!m.overlaps(b"chr1", 99, 99));
    }

    #[test]
    fn nested_regions() {
        let m = Mask::from_regions(b"c:1-1000,c:10-20,c:2000-2001");
        assert!(m.overlaps(b"c", 500, 500));
        assert!(!m.overlaps(b"c", 1500, 1500));
    }
}
