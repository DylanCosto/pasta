//! A region index with payloads (bcftools' `regidx.c`), as `call` uses it for ploidy regions
//! and targets: regions per sequence, sorted by start (longer first on ties, otherwise in
//! insertion order), queried for overlaps in that order.

use std::collections::HashMap;

/// The largest coordinate (0-based) bcftools keeps (`REGIDX_MAX`).
pub const REGIDX_MAX: u32 = 2_147_483_646;

#[derive(Clone, Debug)]
struct RegList<T> {
    regs: Vec<(u32, u32)>,
    payload: Vec<T>,
    /// The largest end among the regions up to each one, for finding overlaps.
    max_end: Vec<u32>,
}

/// bcftools' `regidx_t`.
#[derive(Clone, Debug)]
pub struct RegIdx<T> {
    names: Vec<Vec<u8>>,
    map: HashMap<Vec<u8>, usize>,
    lists: Vec<RegList<T>>,
}

impl<T> Default for RegIdx<T> {
    fn default() -> Self {
        Self {
            names: Vec::new(),
            map: HashMap::new(),
            lists: Vec::new(),
        }
    }
}

/// A region found by [`RegIdx::overlaps`]: its sequence, coordinates (0-based, inclusive),
/// and position in the index (for [`RegIdx::payload_mut`]).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Hit {
    pub seq: usize,
    pub beg: u32,
    pub end: u32,
    pub index: usize,
}

impl<T> RegIdx<T> {
    /// Adds a region (`regidx_push`); [`Self::finish`] must be called before queries.
    pub fn push(&mut self, chr: &[u8], beg: u32, end: u32, payload: T) {
        let (mut beg, mut end) = (beg.min(REGIDX_MAX), end.min(REGIDX_MAX));
        if beg > end {
            std::mem::swap(&mut beg, &mut end);
        }
        let seq = match self.map.get(chr) {
            Some(&i) => i,
            None => {
                self.names.push(chr.to_vec());
                self.map.insert(chr.to_vec(), self.lists.len());
                self.lists.push(RegList {
                    regs: Vec::new(),
                    payload: Vec::new(),
                    max_end: Vec::new(),
                });
                self.lists.len() - 1
            }
        };
        let list = &mut self.lists[seq];
        list.regs.push((beg, end));
        list.payload.push(payload);
    }

    /// Sorts the regions (`reglist_sort_`, a stable sort).
    pub fn finish(&mut self) {
        for list in &mut self.lists {
            let mut order: Vec<usize> = (0..list.regs.len()).collect();
            order.sort_by(|&a, &b| {
                let (ra, rb) = (list.regs[a], list.regs[b]);
                ra.0.cmp(&rb.0).then(rb.1.cmp(&ra.1))
            });
            let regs: Vec<(u32, u32)> = order.iter().map(|&i| list.regs[i]).collect();
            let mut slots: Vec<Option<T>> = std::mem::take(&mut list.payload)
                .into_iter()
                .map(Some)
                .collect();
            list.payload = order
                .iter()
                .map(|&i| slots[i].take().expect("each region is moved once"))
                .collect();
            list.regs = regs;
            let mut m = 0;
            list.max_end = list
                .regs
                .iter()
                .map(|r| {
                    m = m.max(r.1);
                    m
                })
                .collect();
        }
    }

    /// Returns the sequence names in the order they were first seen (`regidx_seq_names`).
    pub fn seq_names(&self) -> &[Vec<u8>] {
        &self.names
    }

    /// Returns the regions of `chr` overlapping `[beg, end]` in index order
    /// (`regidx_overlap` and `regitr_overlap`).
    pub fn overlaps(&self, chr: &[u8], beg: u32, end: u32) -> impl Iterator<Item = Hit> + '_ {
        let seq = self.map.get(chr).copied();
        let (lo, hi) = match seq {
            Some(s) => {
                let list = &self.lists[s];
                let lo = list.max_end.partition_point(|&m| m < beg);
                let hi = list.regs.partition_point(|r| r.0 <= end);
                (lo, hi.max(lo))
            }
            None => (0, 0),
        };
        let regs: &[(u32, u32)] = seq.map_or(&[], |s| &self.lists[s].regs[lo..hi]);
        regs.iter()
            .enumerate()
            .filter(move |(_, r)| r.1 >= beg && r.0 <= end)
            .map(move |(i, r)| Hit {
                seq: seq.unwrap_or(0),
                beg: r.0,
                end: r.1,
                index: lo + i,
            })
    }

    /// Returns the payload of a region found by [`Self::overlaps`].
    pub fn payload(&self, hit: &Hit) -> &T {
        &self.lists[hit.seq].payload[hit.index]
    }

    /// Returns the payload of a region found by [`Self::overlaps`] for changes.
    pub fn payload_mut(&mut self, hit: &Hit) -> &mut T {
        &mut self.lists[hit.seq].payload[hit.index]
    }
}

/// A sequence name and a 0-based inclusive interval.
pub type TabRegion<'a> = (&'a [u8], u32, u32);

/// Parses `CHROM BEG [END]` as `regidx_parse_tab` does: `Ok(None)` for blank and comment
/// lines, and an error message for malformed ones. Returns the name and the 0-based
/// interval, and where the parsing stopped.
pub fn parse_tab(line: &[u8]) -> Result<Option<TabRegion<'_>>, String> {
    let ws = |c: u8| c.is_ascii_whitespace();
    let mut i = line.iter().position(|&c| !ws(c)).unwrap_or(line.len());
    if i == line.len() || line[i] == b'#' {
        return Ok(None);
    }
    let start = i;
    while i < line.len() && !ws(line[i]) {
        i += 1;
    }
    let chr = &line[start..i];
    if i == line.len() {
        return Ok(Some((chr, 0, REGIDX_MAX)));
    }
    let text = String::from_utf8_lossy(line);
    let (beg, used) = strtod_u32(&line[i + 1..]);
    if used == 0 {
        return Err(format!("Could not parse tab line: {text}"));
    }
    if beg == 0 {
        return Err(format!(
            "Could not parse tab line, expected 1-based coordinate: {text}"
        ));
    }
    let beg = beg - 1;
    let se = i + 1 + used;
    let end = if se >= line.len() || se + 1 >= line.len() {
        beg
    } else {
        let (end, used2) = strtod_u32(&line[se + 1..]);
        let after = se + 1 + used2;
        if used2 == 0 || (after < line.len() && !ws(line[after])) {
            beg
        } else if end == 0 {
            return Err(format!(
                "Could not parse tab line, expected 1-based coordinate: {text}"
            ));
        } else {
            end - 1
        }
    };
    Ok(Some((chr, beg, end)))
}

/// C's `strtod` stored into a `uint32_t`: the value and the number of bytes used.
pub fn strtod_u32(s: &[u8]) -> (u32, usize) {
    let (v, end) = pasta::variant::number::strtod(s, 0);
    if end == 0 {
        return (0, 0);
    }
    // C converts the double to uint32_t; out-of-range values are undefined there, so clamp.
    let v = if v.is_nan() || v <= 0.0 {
        0
    } else if v >= f64::from(u32::MAX) {
        u32::MAX
    } else {
        v as u32
    };
    (v, end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_overlaps_in_order() {
        let mut idx = RegIdx::default();
        idx.push(b"1", 10, 20, 'a');
        idx.push(b"1", 5, 7, 'b');
        idx.push(b"1", 10, 30, 'c');
        idx.push(b"2", 0, 0, 'd');
        idx.finish();
        let got: Vec<char> = idx
            .overlaps(b"1", 15, 15)
            .map(|h| *idx.payload(&h))
            .collect();
        assert_eq!(got, ['c', 'a']);
        assert_eq!(idx.overlaps(b"1", 8, 9).count(), 0);
        assert_eq!(idx.overlaps(b"3", 0, 10).count(), 0);
        assert_eq!(idx.seq_names(), [b"1".to_vec(), b"2".to_vec()]);
    }

    #[test]
    fn parses_tab_lines() {
        assert_eq!(
            parse_tab(b"chr1 100 200"),
            Ok(Some((b"chr1".as_slice(), 99, 199)))
        );
        assert_eq!(
            parse_tab(b"chr1\t100"),
            Ok(Some((b"chr1".as_slice(), 99, 99)))
        );
        assert_eq!(parse_tab(b"# x"), Ok(None));
        assert!(parse_tab(b"chr1 0").is_err());
    }
}
