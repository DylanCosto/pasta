//! Tabix: line intervals of tab-delimited files (HTSlib's `tbx.c`), used to index and query
//! bgzipped VCF and region files.

use std::{io, path::Path, sync::atomic::AtomicBool};

use rustc_hash::FxHashMap;

use super::{
    Index, IndexFormat, PushError, TBI_MIN_SHIFT, TBI_N_LVLS, TBX_CONF_LEN, adjust_csi_settings,
};
use crate::variant::{
    log,
    number::{at, strtol},
    record::rlen::svlen_on_ref_for_alt,
};

/// Tabix's generic preset (`TBX_GENERIC`).
pub const PRESET_GENERIC: i32 = 0;
/// Tabix's SAM preset (`TBX_SAM`).
pub const PRESET_SAM: i32 = 1;
/// Tabix's VCF preset (`TBX_VCF`).
pub const PRESET_VCF: i32 = 2;
/// Tabix's flag for 0-based, half-open coordinates (`TBX_UCSC`).
pub const PRESET_UCSC: i32 = 0x10000;

/// The largest shift a tabix-built CSI uses (HTSlib's `TBX_MAX_SHIFT`).
const MAX_SHIFT: u32 = 31;

/// Which columns hold the name and coordinates (HTSlib's `tbx_conf_t`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct TabixConfig {
    /// The preset (`TBX_*`, with [`PRESET_UCSC`] for 0-based coordinates).
    pub preset: i32,
    /// The 1-based sequence-name column.
    pub sc: i32,
    /// The 1-based begin column.
    pub bc: i32,
    /// The 1-based end column (0 for none).
    pub ec: i32,
    /// Lines starting with this byte are comments.
    pub meta_char: i32,
    /// The number of leading lines to skip.
    pub line_skip: i32,
}

impl TabixConfig {
    /// VCF (HTSlib's `tbx_conf_vcf`).
    pub const VCF: Self = Self {
        preset: PRESET_VCF,
        sc: 1,
        bc: 2,
        ec: 0,
        meta_char: b'#' as i32,
        line_skip: 0,
    };

    /// BED (HTSlib's `tbx_conf_bed`).
    pub const BED: Self = Self {
        preset: PRESET_UCSC,
        sc: 1,
        bc: 2,
        ec: 3,
        meta_char: b'#' as i32,
        line_skip: 0,
    };

    fn to_bytes(self) -> [u8; 24] {
        let mut out = [0u8; 24];
        for (i, v) in [
            self.preset,
            self.sc,
            self.bc,
            self.ec,
            self.meta_char,
            self.line_skip,
        ]
        .into_iter()
        .enumerate()
        {
            out[i * 4..i * 4 + 4].copy_from_slice(&v.to_le_bytes());
        }
        out
    }

    /// The preset's name in `get_intv`'s error. HTSlib compares `preset & 0xffff` with
    /// `TBX_UCSC` (0x10000), which never matches, so 0-based presets are named `TBX_GENERIC`.
    fn type_name(&self) -> &'static str {
        match self.preset & 0xffff {
            PRESET_SAM => "TBX_SAM",
            PRESET_VCF => "TBX_VCF",
            3 => "TBX_GAF",
            _ => "TBX_GENERIC",
        }
    }
}

/// A line's interval (HTSlib's `tbx_intv_t`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct Interval {
    /// The start of the name in the line.
    pub name_start: usize,
    /// The end of the name in the line.
    pub name_end: usize,
    /// The 0-based start.
    pub beg: i64,
    /// The 0-based exclusive end.
    pub end: i64,
}

/// C's `strtoll(s, &end, 0)`: decimal, `0x` hex, or `0` octal. Returns the value and end.
fn strtoll0(s: &[u8], start: usize) -> (i64, usize) {
    let mut i = start;
    while crate::variant::number::is_space(at(s, i)) {
        i += 1;
    }
    let neg = match at(s, i) {
        b'-' => {
            i += 1;
            true
        }
        b'+' => {
            i += 1;
            false
        }
        _ => false,
    };
    let (radix, digits_start) = if at(s, i) == b'0'
        && matches!(at(s, i + 1), b'x' | b'X')
        && (at(s, i + 2) as char).is_ascii_hexdigit()
    {
        (16, i + 2)
    } else if at(s, i) == b'0' {
        (8, i)
    } else {
        (10, i)
    };
    let mut j = digits_start;
    let mut n: i128 = 0;
    while let Some(d) = (at(s, j) as char).to_digit(radix) {
        n = (n * i128::from(radix) + i128::from(d)).min(i128::from(i64::MAX) + 1);
        j += 1;
    }
    if j == digits_start {
        return (0, start);
    }
    let v = if neg { -n } else { n };
    (
        v.clamp(i128::from(i64::MIN), i128::from(i64::MAX)) as i64,
        j,
    )
}

/// C's `atoll`.
fn atoll(s: &[u8], start: usize) -> i64 {
    strtol(s, start).0
}

static END_REPORTED: AtomicBool = AtomicBool::new(false);

/// Parses a line's interval (HTSlib's `tbx_parse1`), or `None` when it has none.
#[must_use]
pub fn parse_interval(conf: &TabixConfig, line: &[u8]) -> Option<Interval> {
    let mut intv = Interval {
        name_start: usize::MAX,
        name_end: usize::MAX,
        beg: -1,
        end: -1,
    };
    let len = line.len();
    let preset = conf.preset & 0xffff;
    let mut b = 0usize;
    let mut id = 1;
    let mut getlen = false;
    let mut alcnt = 0usize;
    let mut use_svlen = false;
    let mut lenpos: i64 = -1;
    // The ALT column; HTSlib keeps a bit per allele for SVLEN-on-REF alleles, which
    // `parse_vcf_info` recomputes from the column instead.
    let mut alts = 0..0;
    let (mut reflen, mut svlen, mut fmtlen) = (0i64, 0i64, 0i64);
    let mut i = 0usize;
    'columns: while i <= len {
        // A tab or a NUL ends a column.
        i += memchr::memchr2(b'\t', 0, &line[i..]).unwrap_or(len - i);
        if id == conf.sc {
            intv.name_start = b;
            intv.name_end = i;
        } else if id == conf.bc {
            let (v, s) = strtoll0(line, b);
            intv.beg = v;
            if conf.bc <= conf.ec {
                intv.end = intv.beg;
            }
            if s == b {
                return None;
            }
            if conf.preset & PRESET_UCSC == 0 {
                intv.beg -= 1;
            } else if conf.bc <= conf.ec {
                intv.end += 1;
            }
            if intv.beg < 0 {
                log::warn(
                    "tbx_parse1",
                    format_args!("Coordinate <= 0 detected. Did you forget to use the -0 option?"),
                );
                intv.beg = 0;
            }
            if intv.end < 1 {
                intv.end = 1;
            }
        } else if preset == PRESET_GENERIC {
            if id == conf.ec {
                let (v, s) = strtoll0(line, b);
                if s == b {
                    return None;
                }
                intv.end = v;
            }
        } else if preset == PRESET_SAM {
            if id == 6 {
                let mut l: i64 = 0;
                let mut s = b;
                while s < i {
                    let (x, t) = strtol(&line[..i], s);
                    let op = at(line, t).to_ascii_uppercase();
                    if matches!(op, b'M' | b'D' | b'N') {
                        l += x;
                    }
                    s = t + 1;
                }
                if l == 0 {
                    l = 1;
                }
                intv.end = intv.beg + l;
            }
        } else if preset == PRESET_VCF {
            let field = &line[b..i];
            if id == 4 {
                if b < i {
                    intv.end = intv.beg + (i - b) as i64;
                }
                alcnt += 1;
                reflen = (i - b) as i64;
            }
            if id == 5 {
                alts = b..i;
                for alt in field.split(|&c| c == b',') {
                    alcnt += 1;
                    if svlen_on_ref_for_alt(alt) {
                        use_svlen = true;
                    } else if alt == b"<*>" || alt == b"<NON_REF>" {
                        getlen = true;
                    }
                    if alcnt >= 65536 {
                        break;
                    }
                }
            } else if id == 8 {
                let svlen_alts = if use_svlen { &line[alts.clone()] } else { &[] };
                parse_vcf_info(line, b, i, &mut intv, alcnt, svlen_alts, &mut svlen);
                if !getlen {
                    // Nothing after INFO matters.
                    break 'columns;
                }
            } else if getlen && id == 9 {
                let mut pos = -1i64;
                for key in field.split(|&c| c == b':') {
                    pos += 1;
                    if key == b"LEN" {
                        lenpos = pos;
                        break;
                    }
                }
                if lenpos == -1 {
                    break 'columns;
                }
            } else if id > 9 && getlen && lenpos != -1 {
                let mut tmp = 0i64;
                let mut s = b;
                for d in 0..=lenpos {
                    if d == lenpos {
                        tmp = atoll(&line[..i], s);
                        break;
                    }
                    match line[s..i].iter().position(|&c| c == b':') {
                        Some(t) => s += t + 1,
                        None => break,
                    }
                }
                fmtlen = fmtlen.max(tmp);
            }
        }
        b = i + 1;
        id += 1;
        i += 1;
    }
    if preset == PRESET_VCF {
        let tmp = reflen.max(svlen).max(fmtlen) + intv.beg;
        intv.end = intv.end.max(tmp);
    }
    if intv.name_start == usize::MAX || intv.beg < 0 || intv.end < 0 {
        return None;
    }
    Some(intv)
}

/// The INFO part of [`parse_interval`] (HTSlib's `tbx_parse1`): END and SVLEN.
fn parse_vcf_info(
    line: &[u8],
    b: usize,
    i: usize,
    intv: &mut Interval,
    alcnt: usize,
    svlen_alts: &[u8],
    svlen: &mut i64,
) {
    let info = &line[..i];
    let find =
        |pat: &[u8]| -> Option<usize> { memchr::memmem::find(&info[b..], pat).map(|p| p + b) };
    let mut s = find(b"END=");
    if s == Some(b) {
        s = Some(b + 4);
    } else if s.is_some() {
        s = find(b";END=").map(|p| p + 5);
    }
    if let Some(p) = s
        && at(info, p) != b'.'
    {
        let (end, _) = strtoll0(info, p);
        if end <= intv.beg {
            if !log::fired(&END_REPORTED) {
                let name = if intv.name_start == usize::MAX {
                    String::new()
                } else {
                    String::from_utf8_lossy(&line[intv.name_start..intv.name_end]).into_owned()
                };
                log::warn_once(
                    &END_REPORTED,
                    "tbx_parse1",
                    format_args!(
                        "VCF INFO/END={end} is smaller than POS at {name}:{}\nThis tag will be ignored. Note: only one invalid END tag will be reported.",
                        intv.beg
                    ),
                );
            }
        } else {
            intv.end = end;
        }
    }
    let mut s = find(b"SVLEN=");
    if s == Some(b) {
        s = Some(b + 6);
    } else if s.is_some() {
        s = find(b";SVLEN=").map(|p| p + 7);
    }
    // Whether each ALT allele takes its length from SVLEN (empty when none does).
    let mut on_ref = svlen_alts.split(|&c| c == b',').map(svlen_on_ref_for_alt);
    let mut d = 1usize;
    while let Some(p) = s {
        if d >= alcnt {
            break;
        }
        let comma = memchr::memchr(b',', &info[p..]).map(|t| p + t);
        let tmp = if on_ref.next() == Some(true) {
            atoll(info, p).abs()
        } else {
            1
        };
        *svlen = (*svlen).max(tmp);
        s = comma.map(|t| t + 1);
        d += 1;
    }
}

/// A tabix index with its configuration and sequence names (HTSlib's `tbx_t`).
#[derive(Clone, Debug)]
pub struct Tabix {
    /// The column configuration.
    pub conf: TabixConfig,
    names: Vec<Vec<u8>>,
    dict: FxHashMap<Vec<u8>, u32>,
    /// The binning index.
    pub index: Index,
}

impl Tabix {
    /// Loads a tabix-style index (TBI, or CSI with tabix metadata) from `path` (HTSlib's
    /// `tbx_index_load`).
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be read or has no tabix metadata.
    pub fn load(path: &Path) -> io::Result<Self> {
        let index = Index::read(path)?;
        Self::from_index(index).ok_or_else(|| {
            log::error(
                "index_load",
                format_args!("Invalid index header for {}", path.display()),
            );
            io::Error::new(io::ErrorKind::InvalidData, "invalid index header")
        })
    }

    /// Wraps an index with tabix metadata, or `None` when the metadata is missing.
    #[must_use]
    pub fn from_index(index: Index) -> Option<Self> {
        let meta = index.meta();
        if meta.len() < TBX_CONF_LEN {
            return None;
        }
        let int =
            |k: usize| i32::from_le_bytes(meta[k * 4..k * 4 + 4].try_into().expect("4 bytes"));
        let conf = TabixConfig {
            preset: int(0),
            sc: int(1),
            bc: int(2),
            ec: int(3),
            meta_char: int(4),
            line_skip: int(5),
        };
        let l_nm = int(6) as u32 as usize;
        if l_nm > meta.len() - TBX_CONF_LEN {
            return None;
        }
        let mut tbx = Self {
            conf,
            names: Vec::new(),
            dict: FxHashMap::default(),
            index: Index::new(0, IndexFormat::Csi, 0, TBI_MIN_SHIFT, TBI_N_LVLS),
        };
        let names = &meta[TBX_CONF_LEN..TBX_CONF_LEN + l_nm];
        let mut p = 0;
        while p < names.len() {
            let end = names[p..]
                .iter()
                .position(|&c| c == 0)
                .map_or(names.len(), |e| p + e);
            tbx.add_name(&names[p..end]);
            p = end + 1;
        }
        tbx.index = index;
        Some(tbx)
    }

    /// Returns the ID of a sequence name, adding the name when it is new (HTSlib's `get_tid`
    /// with `is_add`).
    fn add_name(&mut self, name: &[u8]) -> u32 {
        if let Some(&t) = self.dict.get(name) {
            return t;
        }
        let t = self.names.len() as u32;
        self.dict.insert(name.to_vec(), t);
        self.names.push(name.to_vec());
        t
    }

    /// Returns the ID of a sequence name (HTSlib's `tbx_name2id`).
    #[must_use]
    pub fn name_to_id(&self, name: &[u8]) -> Option<u32> {
        self.dict.get(name).copied()
    }

    /// Returns the sequence names in ID order (HTSlib's `tbx_seqnames`).
    #[must_use]
    pub fn names(&self) -> &[Vec<u8>] {
        &self.names
    }

    /// Parses a line's interval and target (HTSlib's `get_intv`), reporting unparsable lines
    /// as HTSlib does.
    ///
    /// # Errors
    ///
    /// Returns an error for an unparsable line or an unknown sequence.
    pub fn line_interval(&self, line: &[u8]) -> io::Result<(u32, Interval)> {
        let intv = self.parse_line(line)?;
        let tid = self
            .name_to_id(&line[intv.name_start..intv.name_end])
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "unknown sequence"))?;
        Ok((tid, intv))
    }

    /// [`Self::line_interval`], adding the sequence name when it is new (`get_intv` with
    /// `is_add`).
    fn line_interval_adding(&mut self, line: &[u8]) -> io::Result<(u32, Interval)> {
        let intv = self.parse_line(line)?;
        Ok((self.add_name(&line[intv.name_start..intv.name_end]), intv))
    }

    /// Parses a line's interval, reporting an unparsable line as HTSlib does.
    fn parse_line(&self, line: &[u8]) -> io::Result<Interval> {
        let Some(intv) = parse_interval(&self.conf, line) else {
            log::error(
                "get_intv",
                format_args!(
                    "Failed to parse {}: was wrong -p [type] used?\nThe offending line was: \"{}\"",
                    self.conf.type_name(),
                    // `%s` stops at a NUL.
                    String::from_utf8_lossy(&line[..memchr::memchr(0, line).unwrap_or(line.len())])
                ),
            );
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "unparsable line",
            ));
        };
        Ok(intv)
    }
}

/// Builds a tabix-style index from lines and their end offsets (HTSlib's `tbx_index`).
#[derive(Debug)]
pub struct TabixBuilder {
    tbx: Tabix,
    fmt: IndexFormat,
    min_shift: u32,
    n_lvls: u32,
    lineno: i64,
    last_off: u64,
    max_ref_len: i64,
    started: bool,
}

impl TabixBuilder {
    /// Starts an index: TBI when `min_shift` is 0, otherwise CSI with that minimum shift.
    #[must_use]
    pub fn new(conf: TabixConfig, min_shift: u32) -> Self {
        let (fmt, min_shift, n_lvls) = if min_shift > 0 {
            (
                IndexFormat::Csi,
                min_shift,
                (MAX_SHIFT - min_shift).div_ceil(3),
            )
        } else {
            (IndexFormat::Tbi, TBI_MIN_SHIFT, TBI_N_LVLS)
        };
        Self {
            tbx: Tabix {
                conf,
                names: Vec::new(),
                dict: FxHashMap::default(),
                index: Index::new(0, fmt, 0, min_shift, n_lvls),
            },
            fmt,
            min_shift,
            n_lvls,
            lineno: 0,
            last_off: 0,
            max_ref_len: 0,
            started: false,
        }
    }

    /// Returns the index format being built.
    #[must_use]
    pub fn format(&self) -> IndexFormat {
        self.fmt
    }

    /// Adds a line (without its newline) and the virtual offset just past it.
    ///
    /// # Errors
    ///
    /// Returns an error for an unparsable line or unsorted records.
    pub fn push_line(&mut self, line: &[u8], end_offset: u64) -> io::Result<()> {
        self.lineno += 1;
        let first = line.first().copied().map_or(0, i32::from);
        let is_meta = first == self.tbx.conf.meta_char;
        if is_meta && self.fmt == IndexFormat::Csi && self.tbx.conf.preset == PRESET_VCF {
            adjust_max_ref_len_vcf(line, &mut self.max_ref_len);
        }
        if self.lineno <= i64::from(self.tbx.conf.line_skip) || is_meta {
            self.last_off = end_offset;
            return Ok(());
        }
        if !self.started {
            if self.fmt == IndexFormat::Csi {
                if self.max_ref_len > 0 {
                    adjust_csi_settings(self.max_ref_len, &mut self.min_shift, &mut self.n_lvls);
                } else {
                    self.n_lvls = if self.min_shift < 10 {
                        9
                    } else if self.min_shift < 25 {
                        9 - (self.min_shift - 10) / 3
                    } else {
                        4
                    };
                }
            }
            self.tbx.index = Index::new(0, self.fmt, self.last_off, self.min_shift, self.n_lvls);
            self.started = true;
        }
        let (tid, intv) = self.tbx.line_interval_adding(line)?;
        self.tbx
            .index
            .push(tid as i32, intv.beg, intv.end, end_offset)
            .map_err(|e: PushError| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    /// Finishes the index at the final offset and sets the tabix metadata.
    #[must_use]
    pub fn finish(mut self, final_offset: u64) -> Tabix {
        if !self.started {
            self.tbx.index = Index::new(0, self.fmt, self.last_off, self.min_shift, self.n_lvls);
        }
        self.tbx.index.finish(final_offset);
        let mut meta = self.tbx.conf.to_bytes().to_vec();
        let l_nm: usize = self.tbx.names.iter().map(|n| n.len() + 1).sum();
        meta.extend_from_slice(&(l_nm as u32).to_le_bytes());
        for n in &self.tbx.names {
            meta.extend_from_slice(n);
            meta.push(0);
        }
        self.tbx.index.set_meta(meta);
        self.tbx
    }
}

/// Reads `##contig=<...length=N...>` lengths (HTSlib's `adjust_max_ref_len_vcf`).
fn adjust_max_ref_len_vcf(line: &[u8], max_ref_len: &mut i64) {
    if !line.starts_with(b"##contig") {
        return;
    }
    let Some(p) = line[8..].windows(6).position(|w| w == b"length") else {
        return;
    };
    let mut p = p + 8 + 6;
    while matches!(at(line, p), b' ' | b'=') {
        p += 1;
    }
    let (len, _) = strtol(line, p);
    if *max_ref_len < len {
        *max_ref_len = len;
    }
}

/// The tabix metadata of a VCF index built while writing (HTSlib's `vcf_idx_init`): the VCF
/// configuration with no names yet.
#[must_use]
pub fn vcf_meta() -> Vec<u8> {
    let mut meta = TabixConfig::VCF.to_bytes().to_vec();
    meta.extend_from_slice(&0u32.to_le_bytes());
    meta
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vcf_intervals() {
        let conf = TabixConfig::VCF;
        let i = parse_interval(&conf, b"1\t100\t.\tACG\tA\t.\t.\t.").unwrap();
        assert_eq!((i.beg, i.end), (99, 102));
        let i = parse_interval(&conf, b"1\t100\t.\tA\t<DEL>\t.\t.\tSVLEN=-50;END=150").unwrap();
        assert_eq!((i.beg, i.end), (99, 150));
        let i = parse_interval(&conf, b"chr2\t5\t.\tA\tC").unwrap();
        assert_eq!(&b"chr2\t5"[i.name_start..i.name_end], b"chr2");
    }

    #[test]
    fn bed_intervals() {
        let i = parse_interval(&TabixConfig::BED, b"1\t10\t20").unwrap();
        assert_eq!((i.beg, i.end), (10, 20));
    }
}
