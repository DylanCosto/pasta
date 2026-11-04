//! bcftools' region index (`regidx.c`) as mpileup uses it: the `-r/-R` regions it visits
//! through the BAM index and the `-t/-T` targets it filters reads and positions by.

use std::{collections::HashMap, fs::File, io::Read as _, path::Path};

use anyhow::Context;
use pasta::variant::{index::split_idx_name, number::strtod};

/// `REGIDX_MAX`: the largest 0-based coordinate.
pub const REGIDX_MAX: u32 = 2_147_483_646;

/// `iBIN`: the index bin of a coordinate.
fn ibin(x: u32) -> usize {
    (x >> 13) as usize
}

/// How a line of a regions file is read (`regidx_parse_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LineFormat {
    /// `CHR:BEG-END` (`regidx_parse_reg`), as in `-r` and `-t`.
    Region,
    /// `CHR BEG END`, 1-based inclusive (`regidx_parse_tab`).
    Tab,
    /// BED, 0-based half-open (`regidx_parse_bed`).
    Bed,
    /// `CHR POS ...` of a VCF (`regidx_parse_vcf`).
    Vcf,
}

impl LineFormat {
    /// The parser `regidx_init` picks for a file name: that of `FILE` for
    /// `FILE##idx##INDEX`, where `regidx_init` looks at the whole name.
    pub fn for_file(name: &str) -> Self {
        let lower = split_idx_name(name).0.to_ascii_lowercase();
        if lower.ends_with(".bed.gz") || lower.ends_with(".bed.bgz") || lower.ends_with(".bed") {
            Self::Bed
        } else if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
            Self::Vcf
        } else {
            Self::Tab
        }
    }
}

/// `reglist_t`: the regions of one sequence.
#[derive(Debug, Default)]
struct RegList {
    regs: Vec<(u32, u32)>,
    idx: Vec<u32>,
    unsorted: bool,
    merged: bool,
    indexed: bool,
}

/// `cmp_regs`: by start, longer intervals first.
fn cmp_regs(a: &(u32, u32), b: &(u32, u32)) -> std::cmp::Ordering {
    a.0.cmp(&b.0).then(b.1.cmp(&a.1))
}

impl RegList {
    fn sort(&mut self) {
        if self.unsorted {
            // A stable sort; with equal regions the order does not matter.
            self.regs.sort_by(cmp_regs);
            self.unsorted = false;
        }
    }

    fn merge(&mut self) {
        if self.merged {
            return;
        }
        let mut out: Vec<(u32, u32)> = Vec::with_capacity(self.regs.len());
        for &r in &self.regs {
            match out.last_mut() {
                Some(last) if last.1 >= r.0 => {
                    if last.1 < r.1 {
                        last.1 = r.1;
                    }
                }
                _ => out.push(r),
            }
        }
        self.regs = out;
        self.merged = true;
    }

    /// `_reglist_build_index`, after [`Self::merge`] when overlaps are merged.
    fn build_index(&mut self) {
        self.sort();
        self.idx.clear();
        for (j, &(beg, end)) in self.regs.iter().enumerate() {
            let (ib, ie) = (ibin(beg), ibin(end));
            if self.idx.len() <= ie {
                self.idx.resize(ie + 1, 0);
            }
            for k in ib..=ie {
                if self.idx[k] == 0 {
                    self.idx[k] = j as u32 + 1;
                }
            }
        }
        self.indexed = true;
    }
}

/// `regidx_t`
#[derive(Debug, Default)]
pub struct RegIdx {
    names: Vec<Vec<u8>>,
    by_name: HashMap<Vec<u8>, usize>,
    lists: Vec<RegList>,
    merge_overlaps: bool,
}

impl RegIdx {
    /// `regidx_insert_list`: parses a `delim`-separated list.
    ///
    /// # Errors
    ///
    /// Returns an error for a line the parser rejects.
    pub fn from_list(list: &str, format: LineFormat) -> anyhow::Result<Self> {
        let mut idx = Self::default();
        for item in list.as_bytes().split(|&c| c == b',') {
            idx.insert(item, format)?;
        }
        Ok(idx)
    }

    /// `regidx_init` of a file; `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read or a line cannot be parsed.
    pub fn from_file(path: &str, format: LineFormat) -> anyhow::Result<Self> {
        let text = read_maybe_compressed(Path::new(split_idx_name(path).0))?;
        let mut idx = Self::default();
        for line in text.split(|&c| c == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            // HTSlib's `hts_getline` loop stops at the first empty line.
            if line.is_empty() {
                break;
            }
            idx.insert(line, format)?;
        }
        Ok(idx)
    }

    /// `regidx_set(merge_overlaps, 1)`
    pub fn set_merge_overlaps(&mut self) {
        self.merge_overlaps = true;
    }

    /// `regidx_insert`
    fn insert(&mut self, line: &[u8], format: LineFormat) -> anyhow::Result<()> {
        let Some((chr, beg, end)) = parse_line(line, format)? else {
            return Ok(());
        };
        self.push(chr, beg, end);
        Ok(())
    }

    /// `regidx_push`
    fn push(&mut self, chr: &[u8], beg: u32, end: u32) {
        let (mut beg, mut end) = (beg.min(REGIDX_MAX), end.min(REGIDX_MAX));
        if beg > end {
            std::mem::swap(&mut beg, &mut end);
        }
        let rid = match self.by_name.get(chr) {
            Some(&rid) => rid,
            None => {
                self.names.push(chr.to_vec());
                self.lists.push(RegList::default());
                self.by_name.insert(chr.to_vec(), self.names.len() - 1);
                self.names.len() - 1
            }
        };
        let list = &mut self.lists[rid];
        list.regs.push((beg, end));
        let n = list.regs.len();
        if !list.unsorted && n > 1 && cmp_regs(&list.regs[n - 2], &list.regs[n - 1]).is_gt() {
            list.unsorted = true;
        }
    }

    /// All regions in `regitr_loop` order: sequences in insertion order, each sorted (and merged
    /// with merge_overlaps). Returns (sequence name, 0-based beg, 0-based inclusive end).
    pub fn all_regions(&mut self) -> Vec<(Vec<u8>, u32, u32)> {
        let mut out = Vec::new();
        for (name, list) in self.names.iter().zip(&mut self.lists) {
            if self.merge_overlaps {
                list.sort();
                list.merge();
            }
            out.extend(list.regs.iter().map(|&(b, e)| (name.clone(), b, e)));
        }
        out
    }

    /// Prepares the per-sequence lookup indexes so [`Self::overlaps`] can take `&self`.
    pub fn build(&mut self) {
        let merge = self.merge_overlaps;
        for list in &mut self.lists {
            if !list.regs.is_empty() {
                if merge {
                    list.sort();
                    list.merge();
                }
                list.build_index();
            }
        }
    }

    /// `regidx_overlap`: the regions overlapping `[beg, end]` (0-based, inclusive), in order.
    /// [`Self::build`] must have been called.
    pub fn overlaps(&self, chr: &[u8], beg: u32, end: u32) -> Overlaps<'_> {
        let empty = Overlaps {
            regs: &[],
            i: 0,
            beg,
            end,
        };
        let Some(&iseq) = self.by_name.get(chr) else {
            return empty;
        };
        let list = &self.lists[iseq];
        if list.regs.is_empty() {
            return empty;
        }
        if list.regs.len() == 1 {
            let r = list.regs[0];
            if beg > r.1 || end < r.0 {
                return empty;
            }
            return Overlaps {
                regs: &list.regs,
                i: 0,
                beg,
                end,
            };
        }
        debug_assert!(list.indexed, "RegIdx::build was not called");
        let ib = ibin(beg);
        let nidx = list.idx.len();
        if ib >= nidx {
            return empty;
        }
        let mut i = list.idx[ib];
        if i == 0 {
            let ie = ibin(end).min(nidx);
            let mut k = ib;
            while k <= ie && k < nidx && list.idx[k] == 0 {
                k += 1;
            }
            if k > ie || k >= nidx {
                return empty;
            }
            i = list.idx[k];
        }
        Overlaps {
            regs: &list.regs,
            i: i as usize - 1,
            beg,
            end,
        }
    }

    /// Whether any region overlaps `[beg, end]`.
    pub fn overlap(&self, chr: &[u8], beg: u32, end: u32) -> bool {
        self.overlaps(chr, beg, end).next().is_some()
    }
}

/// `regitr_overlap`: iterates over the regions overlapping a query.
pub struct Overlaps<'a> {
    regs: &'a [(u32, u32)],
    i: usize,
    beg: u32,
    end: u32,
}

impl Iterator for Overlaps<'_> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<(u32, u32)> {
        while self.i < self.regs.len() {
            let r = self.regs[self.i];
            if r.0 > self.end {
                self.i = self.regs.len();
                return None;
            }
            self.i += 1;
            if r.1 >= self.beg && r.0 <= self.end {
                return Some(r);
            }
        }
        None
    }
}

fn is_space(c: u8) -> bool {
    c.is_ascii_whitespace() || c == 0x0b
}

/// C's `strtod` on the rest of a line, as a `uint32_t` assignment truncates it.
fn strtod_u32(s: &[u8], start: usize) -> (u32, usize) {
    let (v, end) = strtod(s, start);
    (double_to_u32(v), end)
}

/// A C `double` to `uint32_t` conversion as x86-64 builds do it (through a 64-bit integer):
/// negative values wrap instead of saturating at 0.
fn double_to_u32(v: f64) -> u32 {
    v as i64 as u32
}

/// The `regidx_parse_*` functions: `None` for a skipped line.
fn parse_line(line: &[u8], format: LineFormat) -> anyhow::Result<Option<(&[u8], u32, u32)>> {
    let text = || String::from_utf8_lossy(line).into_owned();
    let mut ss = 0;
    while ss < line.len() && is_space(line[ss]) {
        ss += 1;
    }
    if ss == line.len() || line[ss] == b'#' {
        return Ok(None);
    }
    let mut se = ss;
    match format {
        LineFormat::Region => {
            while se < line.len() && line[se] != b':' {
                se += 1;
            }
        }
        _ => {
            while se < line.len() && !is_space(line[se]) {
                se += 1;
            }
        }
    }
    let chr = &line[ss..se];
    if se == line.len() {
        return Ok(Some((chr, 0, REGIDX_MAX)));
    }
    let what = match format {
        LineFormat::Region => "reg",
        LineFormat::Bed => "bed",
        _ => "tab",
    };
    let ss = se + 1;
    let (beg, se) = strtod_u32(line, ss);
    if ss == se {
        anyhow::bail!("Could not parse {what} line: {}", text());
    }
    if format == LineFormat::Bed {
        let ss = se + 1;
        let (end, se2) = strtod(line, ss.min(line.len()));
        if ss >= line.len() || ss == se2 {
            anyhow::bail!("Could not parse bed line: {}", text());
        }
        return Ok(Some((chr, beg, double_to_u32(end - 1.0))));
    }
    if beg == 0 {
        anyhow::bail!(
            "Could not parse {what} line, expected 1-based coordinate: {}",
            text()
        );
    }
    let beg = beg - 1;
    let at = |i: usize| line.get(i).copied().unwrap_or(0);
    let end = if at(se) == 0 || at(se + 1) == 0 {
        if format == LineFormat::Region && at(se) == b'-' {
            REGIDX_MAX
        } else {
            beg
        }
    } else {
        let ss = se + 1;
        let (end, se) = strtod_u32(line, ss);
        let bad_tail = format != LineFormat::Region && at(se) != 0 && !is_space(at(se));
        if ss == se || bad_tail {
            beg
        } else if end == 0 {
            anyhow::bail!(
                "Could not parse {what} line, expected 1-based coordinate: {}",
                text()
            );
        } else {
            end - 1
        }
    };
    let end = if format == LineFormat::Vcf { beg } else { end };
    Ok(Some((chr, beg, end)))
}

/// Whether HTSlib's `find_scheme_handler` sees a URL scheme: 2-11 characters of
/// `[A-Za-z0-9+.-]` before a colon.
fn has_url_scheme(name: &str) -> bool {
    name.find(':').is_some_and(|i| {
        (2..=11).contains(&i)
            && name[..i]
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"+-.".contains(&c))
    })
}

/// Reads a text file, decompressing it if it is gzip/BGZF compressed.
fn read_maybe_compressed(path: &Path) -> anyhow::Result<Vec<u8>> {
    let mut raw = Vec::new();
    File::open(path)
        .and_then(|mut f| f.read_to_end(&mut raw))
        .map_err(|e| {
            let name = path.display().to_string();
            // HTSlib reads `scheme:rest` as a URL, and has no handler for most schemes.
            let why = if has_url_scheme(&name) && e.kind() == std::io::ErrorKind::NotFound {
                "Protocol not supported".to_owned()
            } else {
                crate::raw_input::strerror(&e)
            };
            anyhow::anyhow!("[E::hts_open_format] Failed to open file \"{name}\" : {why}")
        })?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        pasta::bgzf::io::Reader::new(&raw[..])
            .read_to_end(&mut out)
            .with_context(|| format!("Could not read {}", path.display()))?;
        return Ok(out);
    }
    Ok(raw)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_regions() -> anyhow::Result<()> {
        let mut idx = RegIdx::from_list("chr1:10-20,chr1:15-30,chr2,chr1:5", LineFormat::Region)?;
        idx.set_merge_overlaps();
        assert_eq!(
            idx.all_regions(),
            [
                (b"chr1".to_vec(), 4, 4),
                (b"chr1".to_vec(), 9, 29),
                (b"chr2".to_vec(), 0, REGIDX_MAX)
            ]
        );
        Ok(())
    }

    #[test]
    fn converts_negative_coordinates_as_x86_does() -> anyhow::Result<()> {
        // strtod gives -100, which an x86-64 build converts to 2^32 - 100.
        assert_eq!(
            parse_line(b"c1:-100", LineFormat::Region)?,
            Some((&b"c1"[..], 4_294_967_195, 4_294_967_195))
        );
        assert_eq!(
            parse_line(b"c1\t5\t0", LineFormat::Bed)?,
            Some((&b"c1"[..], 5, u32::MAX))
        );
        Ok(())
    }

    #[test]
    fn stops_at_an_empty_line() -> anyhow::Result<()> {
        let dir = std::env::temp_dir().join(format!("vtools-regidx-{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let path = dir.join("regions.txt");
        std::fs::write(&path, "c1\t1\t10\n\nc1\t20\t30\n")?;
        let mut idx = RegIdx::from_file(path.to_str().expect("UTF-8 path"), LineFormat::Tab)?;
        std::fs::remove_file(&path)?;
        std::fs::remove_dir(&dir)?;
        assert_eq!(idx.all_regions(), [(b"c1".to_vec(), 0, 9)]);
        Ok(())
    }

    #[test]
    fn finds_overlaps() -> anyhow::Result<()> {
        let mut idx = RegIdx::from_list("c:100-200,c:300-400,c:100000-100010", LineFormat::Region)?;
        idx.build();
        assert!(idx.overlap(b"c", 150, 150));
        assert!(!idx.overlap(b"c", 250, 250));
        assert!(idx.overlap(b"c", 99_999, 100_000));
        assert!(!idx.overlap(b"d", 150, 150));
        Ok(())
    }
}
