//! An in-memory index of regions with payloads, for overlap queries in any order: a port of
//! bcftools' `regidx.c`.
//!
//! Unlike [`Regions`](super::synced::Regions), which is consumed in genome order as records
//! stream past, a [`RegIdx`] answers arbitrary overlap queries and carries a payload per
//! region (an annotation line, a name). Coordinates are 0-based and inclusive, as `u32` up to
//! [`REGIDX_MAX`].
//!
//! ```
//! use pasta::variant::regidx::{self, RegIdx};
//!
//! # fn main() -> Result<(), pasta::variant::regidx::RegIdxError> {
//! let mut idx: RegIdx<String> = RegIdx::new();
//! idx.insert_text(b"chr1\t100\t200\nchr1\t150\t300\n", |line| {
//!     Ok(regidx::parse_tab(line)?.map(|r| r.with_payload(String::from_utf8_lossy(line).into())))
//! })?;
//! // bcftools: regidx_overlap(idx, "chr1", 179, 179, itr) and regitr_overlap(itr).
//! let hits: Vec<(u32, u32)> = idx
//!     .overlap(b"chr1", 179, 179)
//!     .map(|r| (r.beg, r.end))
//!     .collect();
//! assert_eq!(hits, [(99, 199), (149, 299)]);
//! # Ok(())
//! # }
//! ```

use std::{ops::Range, path::Path};

use rustc_hash::FxHashMap;

use super::{index::data_file, readlist};

/// The largest coordinate (0-based) a region can have (`REGIDX_MAX`).
pub const REGIDX_MAX: u32 = 2_147_483_646;

/// Bins of 8 kb for the lookup table (`iBIN`).
const BIN_SHIFT: u32 = 13;

/// An error parsing or reading regions.
#[derive(Debug)]
pub struct RegIdxError(String);

impl RegIdxError {
    /// Creates an error with a message, for custom line parsers.
    #[must_use]
    pub fn new(msg: impl Into<String>) -> Self {
        Self(msg.into())
    }
}

impl std::fmt::Display for RegIdxError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RegIdxError {}

/// A parsed line: the sequence name's byte range within the line, the 0-based inclusive
/// coordinates, and the payload (bcftools' `regidx_parse_f` outputs).
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct Parsed<P = ()> {
    /// Where the sequence name is in the line.
    pub chr: Range<usize>,
    /// The 0-based start.
    pub beg: u32,
    /// The 0-based inclusive end.
    pub end: u32,
    /// The payload.
    pub payload: P,
}

impl Parsed<()> {
    /// Attaches a payload.
    #[must_use]
    pub fn with_payload<P>(self, payload: P) -> Parsed<P> {
        Parsed {
            chr: self.chr,
            beg: self.beg,
            end: self.end,
            payload,
        }
    }
}

/// The result of a line parser: a region, `None` to skip the line, or an error.
pub type ParseResult<P> = Result<Option<Parsed<P>>, RegIdxError>;

/// One region in a query result.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct Region<'a, P> {
    /// The sequence name.
    pub seq: &'a [u8],
    /// The 0-based start.
    pub beg: u32,
    /// The 0-based inclusive end.
    pub end: u32,
    /// The payload.
    pub payload: &'a P,
}

/// The regions of one sequence (`reglist_t`).
#[derive(Clone, Debug)]
struct RegList<P> {
    regs: Vec<(u32, u32)>,
    payloads: Vec<P>,
    /// For each 8 kb bin, 1 + the first region overlapping it, or 0 (`idx`).
    idx: Vec<u32>,
    nidx: usize,
    indexed: bool,
    unsorted: bool,
    merged: bool,
}

impl<P> Default for RegList<P> {
    fn default() -> Self {
        Self {
            regs: Vec::new(),
            payloads: Vec::new(),
            idx: Vec::new(),
            nidx: 0,
            indexed: false,
            unsorted: false,
            merged: false,
        }
    }
}

/// Orders regions by start, longer first (`cmp_regs`).
fn cmp_regs(a: (u32, u32), b: (u32, u32)) -> std::cmp::Ordering {
    a.0.cmp(&b.0).then(b.1.cmp(&a.1))
}

impl<P> RegList<P> {
    /// Sorts the regions, keeping the order of equal ones (`reglist_sort_`).
    fn sort(&mut self) {
        if !self.unsorted {
            return;
        }
        let mut order: Vec<usize> = (0..self.regs.len()).collect();
        order.sort_by(|&a, &b| cmp_regs(self.regs[a], self.regs[b]));
        let regs: Vec<(u32, u32)> = order.iter().map(|&i| self.regs[i]).collect();
        let mut slots: Vec<Option<P>> = self.payloads.drain(..).map(Some).collect();
        self.payloads = order
            .iter()
            .map(|&i| slots[i].take().expect("each region is moved once"))
            .collect();
        self.regs = regs;
        self.unsorted = false;
    }

    /// Merges overlapping regions, keeping the first one's payload (`reglist_merge_`).
    fn merge(&mut self) {
        if self.merged {
            return;
        }
        let mut j = 1;
        while j < self.regs.len() {
            if self.regs[j - 1].1 < self.regs[j].0 {
                j += 1;
                continue;
            }
            if self.regs[j - 1].1 < self.regs[j].1 {
                self.regs[j - 1].1 = self.regs[j].1;
            }
            self.regs.remove(j);
            self.payloads.remove(j);
        }
        self.merged = true;
    }

    /// Builds the bin lookup table (`_reglist_build_index`).
    fn build_index(&mut self, merge_overlaps: MergeOverlaps) {
        self.sort();
        if merge_overlaps == MergeOverlaps::Yes {
            self.merge();
        }
        self.idx.clear();
        self.nidx = 0;
        for (j, &(beg, end)) in self.regs.iter().enumerate() {
            let ibeg = (beg >> BIN_SHIFT) as usize;
            let iend = (end >> BIN_SHIFT) as usize;
            if self.idx.len() <= iend {
                self.idx.resize((iend + 1).next_power_of_two(), 0);
            }
            for slot in &mut self.idx[ibeg..=iend] {
                if *slot == 0 {
                    *slot = j as u32 + 1;
                }
            }
            self.nidx = self.nidx.max(iend + 1);
        }
        self.indexed = true;
    }

    /// Finds the first region overlapping `[beg, end]` (`regidx_overlap`).
    fn first_overlap(
        &mut self,
        beg: u32,
        end: u32,
        merge_overlaps: MergeOverlaps,
    ) -> Option<usize> {
        if self.regs.is_empty() {
            return None;
        }
        if self.regs.len() == 1 {
            let r = self.regs[0];
            return (beg <= r.1 && end >= r.0).then_some(0);
        }
        if !self.indexed {
            self.build_index(merge_overlaps);
        }
        let ibeg = (beg >> BIN_SHIFT) as usize;
        if ibeg >= self.nidx {
            return None;
        }
        let mut i = self.idx[ibeg];
        if i == 0 {
            let iend = ((end >> BIN_SHIFT) as usize).min(self.nidx);
            i = (ibeg..=iend)
                .map(|b| self.idx.get(b).copied().unwrap_or(0))
                .find(|&v| v != 0)?;
        }
        self.next_overlap(i as usize - 1, beg, end)
    }

    /// Finds the next region from `from` overlapping `[beg, end]` (`regitr_overlap`).
    fn next_overlap(&self, from: usize, beg: u32, end: u32) -> Option<usize> {
        for (k, &(rb, re)) in self.regs.iter().enumerate().skip(from) {
            if rb > end {
                return None;
            }
            if re >= beg && rb <= end {
                return Some(k);
            }
        }
        None
    }
}

/// Whether overlapping regions are merged into one (`regidx_set(idx, merge_overlaps, ..)`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub enum MergeOverlaps {
    /// Keep every region (the default).
    #[default]
    No,
    /// Merge overlapping regions, keeping the first one's payload.
    Yes,
}

/// An index of regions with payloads (bcftools' `regidx_t`).
#[derive(Clone, Debug)]
pub struct RegIdx<P = ()> {
    seqs: Vec<RegList<P>>,
    seq_names: Vec<Vec<u8>>,
    seq_map: FxHashMap<Vec<u8>, usize>,
    merge_overlaps: MergeOverlaps,
}

impl<P> Default for RegIdx<P> {
    fn default() -> Self {
        Self::new()
    }
}

impl RegIdx<()> {
    /// Reads regions from a file without payloads, choosing the parser by the name as
    /// bcftools' `regidx_init` does with no parser: `.bed`, `.bed.gz`, `.bed.bgz` are BED;
    /// `.vcf`, `.vcf.gz` are VCF (one base at POS); anything else is CHROM, POS[, END]. The
    /// name of a `FILE##idx##INDEX` path is that of `FILE` (HTSlib looks at the whole name).
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be read or a line cannot be parsed.
    pub fn from_file(path: &Path) -> Result<Self, RegIdxError> {
        let name = data_file(path).to_string_lossy().to_ascii_lowercase();
        let parser: fn(&[u8]) -> ParseResult<()> =
            if name.ends_with(".bed.gz") || name.ends_with(".bed.bgz") || name.ends_with(".bed") {
                parse_bed
            } else if name.ends_with(".vcf") || name.ends_with(".vcf.gz") {
                parse_vcf
            } else {
                parse_tab
            };
        Self::from_file_with(path, parser)
    }
}

impl<P> RegIdx<P> {
    /// Creates an empty index (`regidx_init(NULL, ...)`).
    #[must_use]
    pub fn new() -> Self {
        Self {
            seqs: Vec::new(),
            seq_names: Vec::new(),
            seq_map: FxHashMap::default(),
            merge_overlaps: MergeOverlaps::No,
        }
    }

    /// Sets whether overlapping regions are merged (`regidx_set(idx, merge_overlaps, ..)`).
    /// Takes effect when the index is first queried.
    pub fn set_merge_overlaps(&mut self, merge: MergeOverlaps) {
        self.merge_overlaps = merge;
    }

    /// Reads regions from a plain or gzip-compressed file, one per line, with `parser`
    /// (`regidx_init`). Reading stops at the first empty line, as in bcftools. A
    /// `FILE##idx##INDEX` path reads `FILE`, as `hts_open` does.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be read or the parser fails.
    pub fn from_file_with<F>(path: &Path, mut parser: F) -> Result<Self, RegIdxError>
    where
        F: FnMut(&[u8]) -> ParseResult<P>,
    {
        let lines = readlist::read_all_lines(data_file(path))
            .map_err(|e| RegIdxError(format!("{}: {e}", path.display())))?;
        let mut idx = Self::new();
        for line in lines.iter().take_while(|l| !l.is_empty()) {
            idx.insert(line, &mut parser)?;
        }
        Ok(idx)
    }

    /// Adds the regions of a text, one per line (`regidx_init_string`); unparsable lines are
    /// skipped, as bcftools ignores them.
    ///
    /// # Errors
    ///
    /// Never fails; the `Result` is for parsers that want to stop early.
    pub fn insert_text<F>(&mut self, text: &[u8], mut parser: F) -> Result<(), RegIdxError>
    where
        F: FnMut(&[u8]) -> ParseResult<P>,
    {
        let mut ss = 0;
        while ss < text.len() {
            while ss < text.len() && is_space(text[ss]) {
                ss += 1;
            }
            let mut se = ss;
            while se < text.len() && text[se] != b'\r' && text[se] != b'\n' {
                se += 1;
            }
            // bcftools ignores the result of each insertion here.
            let _ = self.insert(&text[ss..se], &mut parser);
            while se < text.len() && is_space(text[se]) {
                se += 1;
            }
            ss = se;
        }
        Ok(())
    }

    /// Parses one line and adds its region (`regidx_insert`).
    ///
    /// # Errors
    ///
    /// Returns the parser's error.
    pub fn insert<F>(&mut self, line: &[u8], parser: &mut F) -> Result<(), RegIdxError>
    where
        F: FnMut(&[u8]) -> ParseResult<P>,
    {
        if let Some(p) = parser(line)? {
            self.push(&line[p.chr], p.beg, p.end, p.payload);
        }
        Ok(())
    }

    /// Adds the regions of a `delim`-separated list (`regidx_insert_list`).
    ///
    /// # Errors
    ///
    /// Returns the parser's error.
    pub fn insert_list<F>(
        &mut self,
        list: &[u8],
        delim: u8,
        mut parser: F,
    ) -> Result<(), RegIdxError>
    where
        F: FnMut(&[u8]) -> ParseResult<P>,
    {
        // As in bcftools, a trailing delimiter does not add an empty item.
        let list = list.strip_suffix(&[delim]).unwrap_or(list);
        if list.is_empty() {
            return Ok(());
        }
        for item in list.split(|&c| c == delim) {
            self.insert(item, &mut parser)?;
        }
        Ok(())
    }

    /// Adds a region with 0-based inclusive coordinates (`regidx_push`); coordinates are
    /// capped at [`REGIDX_MAX`] and swapped when reversed.
    pub fn push(&mut self, chr: &[u8], beg: u32, end: u32, payload: P) {
        let (mut beg, mut end) = (beg.min(REGIDX_MAX), end.min(REGIDX_MAX));
        if beg > end {
            std::mem::swap(&mut beg, &mut end);
        }
        let rid = match self.seq_map.get(chr) {
            Some(&r) => r,
            None => {
                let r = self.seqs.len();
                self.seq_names.push(chr.to_vec());
                self.seq_map.insert(chr.to_vec(), r);
                self.seqs.push(RegList::default());
                r
            }
        };
        let list = &mut self.seqs[rid];
        if let Some(&last) = list.regs.last()
            && cmp_regs(last, (beg, end)).is_gt()
        {
            list.unsorted = true;
        }
        list.regs.push((beg, end));
        list.payloads.push(payload);
        list.indexed = false;
        list.merged = false;
    }

    /// Returns the sequence names in the order first seen (`regidx_seq_names`).
    #[must_use]
    pub fn seq_names(&self) -> &[Vec<u8>] {
        &self.seq_names
    }

    /// Returns the number of regions of a sequence (`regidx_seq_nregs`).
    #[must_use]
    pub fn seq_region_count(&self, seq: &[u8]) -> usize {
        self.seq_map
            .get(seq)
            .map_or(0, |&i| self.seqs[i].regs.len())
    }

    /// Returns the total number of regions (`regidx_nregs`).
    #[must_use]
    pub fn region_count(&self) -> usize {
        self.seqs.iter().map(|s| s.regs.len()).sum()
    }

    /// Returns whether `chr:[beg, end]` (0-based inclusive) overlaps a region
    /// (`regidx_overlap` without an iterator).
    pub fn overlaps(&mut self, chr: &[u8], beg: u32, end: u32) -> bool {
        self.overlap(chr, beg, end).next().is_some()
    }

    /// Returns the regions overlapping `chr:[beg, end]` (0-based inclusive), sorted by start
    /// (`regidx_overlap` then `regitr_overlap`).
    pub fn overlap(&mut self, chr: &[u8], beg: u32, end: u32) -> Overlaps<'_, P> {
        let merge = self.merge_overlaps;
        let found = self.seq_map.get(chr).copied().and_then(|iseq| {
            let first = self.seqs[iseq].first_overlap(beg, end, merge)?;
            Some((iseq, first))
        });
        Overlaps {
            idx: self,
            found,
            beg,
            end,
        }
    }

    /// Returns every region, sequence by sequence (`regitr_loop`). Regions are in insertion
    /// order, or sorted and merged with [`MergeOverlaps::Yes`].
    pub fn iter(&mut self) -> impl Iterator<Item = Region<'_, P>> + '_ {
        if self.merge_overlaps == MergeOverlaps::Yes {
            for list in &mut self.seqs {
                list.sort();
                list.merge();
            }
        }
        self.seqs
            .iter()
            .zip(&self.seq_names)
            .flat_map(|(list, name)| {
                list.regs
                    .iter()
                    .zip(&list.payloads)
                    .map(move |(&(beg, end), payload)| Region {
                        seq: name,
                        beg,
                        end,
                        payload,
                    })
            })
    }
}

/// The regions overlapping a query (bcftools' `regitr_t` after `regidx_overlap`).
#[derive(Debug)]
pub struct Overlaps<'a, P> {
    idx: &'a RegIdx<P>,
    /// The sequence and the next overlapping region.
    found: Option<(usize, usize)>,
    beg: u32,
    end: u32,
}

impl<'a, P> Iterator for Overlaps<'a, P> {
    type Item = Region<'a, P>;

    fn next(&mut self) -> Option<Self::Item> {
        let (iseq, k) = self.found?;
        let list = &self.idx.seqs[iseq];
        self.found = list
            .next_overlap(k + 1, self.beg, self.end)
            .map(|n| (iseq, n));
        Some(Region {
            seq: &self.idx.seq_names[iseq],
            beg: list.regs[k].0,
            end: list.regs[k].1,
            payload: &list.payloads[k],
        })
    }
}

/// C's `isspace`.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// C's `strtod` on a prefix: the value and the index after it (`start` when nothing parses).
fn strtod(s: &[u8], start: usize) -> (f64, usize) {
    if start >= s.len() {
        return (0.0, start);
    }
    let mut i = start;
    while i < s.len() && is_space(s[i]) {
        i += 1;
    }
    let num_start = i;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        i += 1;
    }
    let rest = &s[i..];
    for (word, value) in [
        (&b"infinity"[..], f64::INFINITY),
        (b"inf", f64::INFINITY),
        (b"nan", f64::NAN),
    ] {
        if rest.len() >= word.len() && rest[..word.len()].eq_ignore_ascii_case(word) {
            let neg = s[num_start] == b'-';
            let end = i + word.len();
            return (if neg { -value } else { value }, end);
        }
    }
    let digits_start = i;
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
    }
    let mut ndigits = i - digits_start;
    if i < s.len() && s[i] == b'.' {
        let frac = i + 1;
        let mut j = frac;
        while j < s.len() && s[j].is_ascii_digit() {
            j += 1;
        }
        ndigits += j - frac;
        if ndigits > 0 {
            i = j;
        }
    }
    if ndigits == 0 {
        return (0.0, start);
    }
    if i < s.len() && (s[i] == b'e' || s[i] == b'E') {
        let mut j = i + 1;
        if j < s.len() && (s[j] == b'+' || s[j] == b'-') {
            j += 1;
        }
        let exp_digits = j;
        while j < s.len() && s[j].is_ascii_digit() {
            j += 1;
        }
        if j > exp_digits {
            i = j;
        }
    }
    let text = std::str::from_utf8(&s[num_start..i]).unwrap_or("0");
    (text.parse().unwrap_or(0.0), i)
}

/// The sequence name at the start of a line (after spaces), up to `stop`; `None` for blank
/// lines and comments.
fn chrom_span(line: &[u8], stop: impl Fn(u8) -> bool) -> Option<(usize, usize)> {
    let mut ss = 0;
    while ss < line.len() && is_space(line[ss]) {
        ss += 1;
    }
    if ss == line.len() || line[ss] == b'#' {
        return None;
    }
    let mut se = ss;
    while se < line.len() && !stop(line[se]) {
        se += 1;
    }
    Some((ss, se))
}

/// Converts a C double to `uint32_t` (saturating where C is undefined).
fn to_u32(v: f64) -> u32 {
    v as u32
}

fn whole_sequence(ss: usize, se: usize) -> Parsed {
    Parsed {
        chr: ss..se,
        beg: 0,
        end: REGIDX_MAX,
        payload: (),
    }
}

/// Parses CHROM, or whitespace-separated CHROM, FROM, TO (0-based, half-open)
/// (`regidx_parse_bed`).
///
/// # Errors
///
/// Returns an error when a coordinate is missing.
pub fn parse_bed(line: &[u8]) -> ParseResult<()> {
    let Some((ss, se)) = chrom_span(line, is_space) else {
        return Ok(None);
    };
    if se == line.len() {
        return Ok(Some(whole_sequence(ss, se)));
    }
    let fail = || {
        RegIdxError(format!(
            "Could not parse bed line: {}",
            String::from_utf8_lossy(line)
        ))
    };
    let (beg, e) = strtod(line, se + 1);
    if e == se + 1 {
        return Err(fail());
    }
    let (end, e2) = strtod(line, e + 1);
    if e2 == e + 1 {
        return Err(fail());
    }
    Ok(Some(Parsed {
        chr: ss..se,
        beg: to_u32(beg),
        end: to_u32(end - 1.0),
        payload: (),
    }))
}

/// Parses CHROM, or whitespace-separated CHROM, POS[, END] (1-based, inclusive)
/// (`regidx_parse_tab`).
///
/// # Errors
///
/// Returns an error when a coordinate is missing or 0.
pub fn parse_tab(line: &[u8]) -> ParseResult<()> {
    let Some((ss, se)) = chrom_span(line, is_space) else {
        return Ok(None);
    };
    if se == line.len() {
        return Ok(Some(whole_sequence(ss, se)));
    }
    let text = String::from_utf8_lossy(line);
    let (beg, e) = strtod(line, se + 1);
    if e == se + 1 {
        return Err(RegIdxError(format!("Could not parse tab line: {text}")));
    }
    let beg = to_u32(beg);
    if beg == 0 {
        return Err(RegIdxError(format!(
            "Could not parse tab line, expected 1-based coordinate: {text}"
        )));
    }
    let beg = beg - 1;
    let end = if e >= line.len() || e + 1 >= line.len() {
        beg
    } else {
        let (v, e2) = strtod(line, e + 1);
        if e2 == e + 1 || (e2 < line.len() && !is_space(line[e2])) {
            beg
        } else {
            let v = to_u32(v);
            if v == 0 {
                return Err(RegIdxError(format!(
                    "Could not parse tab line, expected 1-based coordinate: {text}"
                )));
            }
            v - 1
        }
    };
    Ok(Some(Parsed {
        chr: ss..se,
        beg,
        end,
        payload: (),
    }))
}

/// Parses the CHROM and POS of a VCF line as a one-base region (`regidx_parse_vcf`).
///
/// # Errors
///
/// Returns an error when POS is missing or 0.
pub fn parse_vcf(line: &[u8]) -> ParseResult<()> {
    Ok(parse_tab(line)?.map(|mut p| {
        p.end = p.beg;
        p
    }))
}

/// Parses CHROM, CHROM:POS, CHROM:FROM-TO, or CHROM:FROM- (1-based, inclusive)
/// (`regidx_parse_reg`).
///
/// # Errors
///
/// Returns an error when a coordinate is malformed or 0.
pub fn parse_reg(line: &[u8]) -> ParseResult<()> {
    let Some((ss, se)) = chrom_span(line, |c| c == b':') else {
        return Ok(None);
    };
    if se == line.len() {
        return Ok(Some(whole_sequence(ss, se)));
    }
    let text = String::from_utf8_lossy(line);
    let (beg, e) = strtod(line, se + 1);
    if e == se + 1 {
        return Err(RegIdxError(format!("Could not parse reg line: {text}")));
    }
    let beg = to_u32(beg);
    if beg == 0 {
        return Err(RegIdxError(format!(
            "Could not parse reg line, expected 1-based coordinate: {text}"
        )));
    }
    let beg = beg - 1;
    let end = if e >= line.len() || e + 1 >= line.len() {
        if line.get(e) == Some(&b'-') {
            REGIDX_MAX
        } else {
            beg
        }
    } else {
        let (v, e2) = strtod(line, e + 1);
        if e2 == e + 1 {
            beg
        } else {
            let v = to_u32(v);
            if v == 0 {
                return Err(RegIdxError(format!(
                    "Could not parse reg line, expected 1-based coordinate: {text}"
                )));
            }
            v - 1
        }
    };
    Ok(Some(Parsed {
        chr: ss..se,
        beg,
        end,
        payload: (),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(text: &str, parser: fn(&[u8]) -> ParseResult<()>) -> RegIdx {
        let mut idx = RegIdx::new();
        idx.insert_text(text.as_bytes(), parser)
            .expect("the text parses");
        idx
    }

    fn hits(idx: &mut RegIdx, chr: &str, beg: u32, end: u32) -> Vec<(u32, u32)> {
        idx.overlap(chr.as_bytes(), beg, end)
            .map(|r| (r.beg, r.end))
            .collect()
    }

    #[test]
    fn idx_names_read_the_data_file_by_its_name() {
        let dir = std::env::temp_dir().join(format!("pasta-regidx-idx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let bed = dir.join("r.bed");
        std::fs::write(&bed, "chr1\t100\t200\n").expect("write the test file");
        let mut name = bed.into_os_string();
        name.push("##idx##other.tbi");
        let idx = RegIdx::from_file(Path::new(&name));
        let _ = std::fs::remove_dir_all(&dir);
        let mut idx = idx.expect("the data file is read");
        assert_eq!(hits(&mut idx, "chr1", 0, 1000), [(100, 199)]);
    }

    #[test]
    fn parsers_follow_bcftools() {
        let p = |r: ParseResult<()>| r.expect("parses").map(|p| (p.chr, p.beg, p.end));
        assert_eq!(p(parse_tab(b"chr1\t100")), Some((0..4, 99, 99)));
        assert_eq!(p(parse_tab(b"chr1\t100\t200")), Some((0..4, 99, 199)));
        assert_eq!(p(parse_tab(b"chr1\t100\tx")), Some((0..4, 99, 99)));
        assert_eq!(p(parse_tab(b"  chr1")), Some((2..6, 0, REGIDX_MAX)));
        assert_eq!(p(parse_tab(b"# comment")), None);
        assert!(parse_tab(b"chr1\t0").is_err());
        assert_eq!(p(parse_bed(b"chr1 99 200")), Some((0..4, 99, 199)));
        assert_eq!(p(parse_vcf(b"chr1\t100\t.\tA\tC")), Some((0..4, 99, 99)));
        assert_eq!(p(parse_reg(b"chr1:100-200")), Some((0..4, 99, 199)));
        assert_eq!(p(parse_reg(b"chr1:100-")), Some((0..4, 99, REGIDX_MAX)));
        assert_eq!(p(parse_reg(b"chr1:1e3")), Some((0..4, 999, 999)));
    }

    #[test]
    fn overlaps_across_bins() {
        let mut idx = build(
            "1\t100\t200\n1\t50000\t60000\n1\t150\t100000\n2\t10\n1\t5\t5\n",
            parse_tab,
        );
        assert_eq!(hits(&mut idx, "1", 160, 170), [(99, 199), (149, 99_999)]);
        assert_eq!(
            hits(&mut idx, "1", 55000, 55000),
            [(149, 99_999), (49_999, 59_999)]
        );
        assert_eq!(hits(&mut idx, "1", 4, 4), [(4, 4)]);
        assert!(hits(&mut idx, "1", 200_000, 300_000).is_empty());
        assert_eq!(hits(&mut idx, "2", 0, 100), [(9, 9)]);
        assert!(hits(&mut idx, "3", 0, 100).is_empty());
        assert_eq!(idx.region_count(), 5);
        assert_eq!(idx.seq_region_count(b"1"), 4);
    }

    #[test]
    fn payloads_follow_their_regions() {
        let mut idx: RegIdx<&str> = RegIdx::new();
        idx.push(b"1", 300, 400, "c");
        idx.push(b"1", 100, 200, "a");
        idx.push(b"1", 100, 250, "b");
        let got: Vec<&str> = idx.overlap(b"1", 150, 350).map(|r| *r.payload).collect();
        assert_eq!(got, ["b", "a", "c"]);
        idx.set_merge_overlaps(MergeOverlaps::Yes);
        let mut merged = idx.clone();
        let all: Vec<(u32, u32, &str)> =
            merged.iter().map(|r| (r.beg, r.end, *r.payload)).collect();
        assert_eq!(all, [(100, 250, "b"), (300, 400, "c")]);
    }
}
