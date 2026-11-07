//! Binning indexes for VCF and BCF (HTSlib's `hts_idx_t`, `hts_itr_query`, and tabix's
//! `tbx_t`): building them as HTSlib does (`hts_idx_push`/`hts_idx_finish`), reading and
//! writing CSI and TBI files, and planning region queries.
//!
//! Indexes built here have the same bins, chunks, linear indexes, and per-contig counts as
//! HTSlib's, and list each target's bins in the order of HTSlib's hash table, and are saved in
//! HTSlib's BGZF blocks, so an index of the same data file has the same data as HTSlib's, block
//! for block. Its compressed bytes are zflate's, at pasta's default level, not libdeflate's.

mod khash;
pub mod tabix;

use std::{
    fs::File,
    io::{self, BufWriter, Read, Write},
    ops::Range,
    path::Path,
};

use indexmap::IndexMap;

use self::khash::BinTable;
use crate::bgzf;
use crate::variant::{
    io::{bgzf_out::BlockWriter, writer::htslib_level},
    log,
};

/// The index flavor (HTSlib's `HTS_FMT_CSI`/`HTS_FMT_TBI`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum IndexFormat {
    /// Coordinate-sorted index (`.csi`).
    Csi,
    /// Tabix index (`.tbi`), with a fixed 14-bit minimum shift and 5 levels.
    Tbi,
}

/// The minimum distance between chunks merged into a parent bin (HTSlib's
/// `HTS_MIN_MARKER_DIST`).
const MIN_MARKER_DIST: u64 = 0x10000;

/// The minimum shift of every TBI index.
pub(crate) const TBI_MIN_SHIFT: u32 = 14;

/// The number of levels of every TBI index.
pub(crate) const TBI_N_LVLS: u32 = 5;

/// The length of the tabix configuration that starts a TBI index's metadata: seven
/// little-endian 32-bit integers, the last one the length of the sequence names that follow.
pub(crate) const TBX_CONF_LEN: usize = 28;

/// Where the sequence names' length sits in the tabix configuration.
const TBX_L_NM: std::ops::Range<usize> = 24..TBX_CONF_LEN;

/// The separator of an explicit index name in an input name, `FILE##idx##INDEX` (HTSlib's
/// `HTS_IDX_DELIM`).
const IDX_DELIM: &str = "##idx##";

/// Splits an input name at `##idx##` into the data path and the explicit index path, as
/// `hts_open` and `hts_idx_load` do; a name without `##idx##` is all data path.
#[must_use]
pub fn split_idx_name(name: &str) -> (&str, Option<&str>) {
    match name.find(IDX_DELIM) {
        Some(i) => (&name[..i], Some(&name[i + IDX_DELIM.len()..])),
        None => (name, None),
    }
}

/// Returns the data file of a path that may carry an `##idx##` index name, the file
/// `hts_open` opens. A path that is not UTF-8 is returned whole.
#[must_use]
pub fn data_file(path: &Path) -> &Path {
    path.to_str()
        .map_or(path, |s| Path::new(split_idx_name(s).0))
}

/// The first bin of level `l` (HTSlib's `hts_bin_first`).
#[inline]
#[must_use]
pub fn bin_first(l: u32) -> u32 {
    ((1u32 << (l * 3)) - 1) / 7
}

/// The parent of a bin (HTSlib's `hts_bin_parent`).
#[inline]
#[must_use]
pub fn bin_parent(b: u32) -> u32 {
    (b.wrapping_sub(1)) >> 3
}

/// The level of a bin (HTSlib's `hts_bin_level`).
#[inline]
#[must_use]
pub fn bin_level(bin: u32) -> u32 {
    let mut l = 0;
    let mut b = bin;
    while b != 0 {
        l += 1;
        b = bin_parent(b);
    }
    l
}

/// The first linear-index window of a bin (HTSlib's `hts_bin_bot`).
#[inline]
#[must_use]
pub fn bin_bot(bin: u32, n_lvls: u32) -> u64 {
    let l = bin_level(bin);
    u64::from(bin - bin_first(l)) << ((n_lvls - l) * 3)
}

/// The exclusive maximum position a binning scheme covers (HTSlib's `hts_bin_maxpos`).
#[inline]
#[must_use]
pub fn bin_maxpos(min_shift: u32, n_lvls: u32) -> i64 {
    1i64 << (min_shift + n_lvls * 3)
}

/// The smallest bin containing `[beg, end)` (HTSlib's `hts_reg2bin`).
#[inline]
#[must_use]
pub fn reg2bin(beg: i64, end: i64, min_shift: u32, n_lvls: u32) -> u32 {
    let end = end - 1;
    let mut s = min_shift;
    let mut t = bin_first(n_lvls);
    let mut l = n_lvls;
    while l > 0 {
        if beg >> s == end >> s {
            // As C's int arithmetic: a record without coordinates (-1..0) is in bin `t - 1`.
            return t.wrapping_add((beg >> s) as u32);
        }
        l -= 1;
        s += 3;
        t -= 1 << (l * 3);
    }
    0
}

/// Adjusts the CSI levels (or min shift) so the longest contig fits (HTSlib's
/// `hts_adjust_csi_settings`).
pub fn adjust_csi_settings(max_len_in: i64, min_shift: &mut u32, n_lvls: &mut u32) {
    const MAX_N_LVLS: u32 = 9;
    let max_len = max_len_in + 256;
    if max_len <= bin_maxpos(*min_shift, MAX_N_LVLS) {
        let mut maxpos = bin_maxpos(*min_shift, *n_lvls);
        while max_len > maxpos {
            *n_lvls += 1;
            maxpos *= 8;
        }
    } else {
        let mut shift = *min_shift;
        let mut maxpos = bin_maxpos(shift, MAX_N_LVLS);
        while max_len > maxpos {
            shift += 1;
            maxpos *= 2;
        }
        log::warn(
            "hts_adjust_csi_settings",
            format_args!(
                "Adjusted min_shift from {} to {shift} due to longest reference of {max_len_in} bases.",
                *min_shift
            ),
        );
        *n_lvls = MAX_N_LVLS;
        *min_shift = shift;
    }
}

/// A bin's chunks and linear offset (HTSlib's `bins_t`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Bin {
    /// The smallest file offset of records in or below this bin (CSI).
    pub loff: u64,
    /// The `(begin, end)` virtual offsets of the chunks.
    pub chunks: Vec<(u64, u64)>,
}

/// The state of an index being built (HTSlib's `idx->z`).
#[derive(Clone, Debug)]
struct BuildState {
    last_bin: u32,
    save_bin: u32,
    last_coor: i64,
    last_tid: i32,
    save_tid: i32,
    finished: bool,
    last_off: u64,
    save_off: u64,
    off_beg: u64,
    off_end: u64,
    /// The records of the current target so far, for the pseudo-bin; the unmapped count
    /// written beside it is always 0, as VCF and BCF records are all mapped.
    n_mapped: u64,
}

/// An error building an index from records.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum PushError {
    /// The record does not fit the binning scheme.
    OutOfRange,
    /// Records are not sorted by position.
    Unsorted,
    /// A contig's records are not contiguous.
    NotContiguous,
    /// A record ends before it begins.
    InvalidInterval,
    /// Records without coordinates are followed by records with them.
    NoCoorNotLast,
}

impl std::fmt::Display for PushError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::OutOfRange => "region cannot be stored in the index",
            Self::Unsorted => "unsorted positions",
            Self::NotContiguous => "chromosome blocks not continuous",
            Self::InvalidInterval => "invalid record interval",
            Self::NoCoorNotLast => "records without coordinates are not at the end",
        })
    }
}

impl std::error::Error for PushError {}

/// The last record pushed, for the order checks of `hts_idx_push`.
#[derive(Clone, Copy, Debug)]
struct Last {
    tid: i32,
    coor: i64,
    /// Whether records without coordinates were pushed.
    no_coor: bool,
}

impl Last {
    /// The state of a new index.
    const NEW: Self = Self {
        tid: -1,
        coor: 0xffff_ffff,
        no_coor: false,
    };

    /// Checks that a record (with `beg` -1 and `end` 0 when it has no coordinates) may follow
    /// this one, as `hts_idx_push` does, and reports an error as it does. `seen` says whether
    /// records of `tid` came before. Returns whether the record starts a new target.
    #[inline]
    fn check(
        self,
        tid: i32,
        beg: i64,
        end: i64,
        seen: impl FnOnce() -> bool,
    ) -> Result<bool, PushError> {
        let new_target = self.tid != tid || (self.tid >= 0 && tid < 0);
        if new_target {
            if tid >= 0 && self.no_coor {
                log::error(
                    "hts_idx_push",
                    format_args!(
                        "NO_COOR reads not in a single block at the end {tid} {}",
                        self.tid
                    ),
                );
                return Err(PushError::NoCoorNotLast);
            }
            if tid >= 0 && seen() {
                log::error(
                    "hts_idx_push",
                    format_args!("Chromosome blocks not continuous"),
                );
                return Err(PushError::NotContiguous);
            }
        } else if tid >= 0 && self.coor > beg {
            log::error(
                "hts_idx_push",
                format_args!(
                    "Unsorted positions on sequence #{}: {} followed by {}",
                    tid + 1,
                    self.coor + 1,
                    beg + 1
                ),
            );
            return Err(PushError::Unsorted);
        }
        if end < beg {
            log::error(
                "hts_idx_push",
                format_args!(
                    "Invalid record on sequence #{}: end {end} < begin {}",
                    tid + 1,
                    beg + 1
                ),
            );
            return Err(PushError::InvalidInterval);
        }
        Ok(new_target)
    }
}

/// The checks [`Index::push`] makes on a record, made when the record is written rather than
/// when it is pushed. An index built while writing with compression threads pushes a record
/// once the address of its block is known, which can be many blocks later; checking it as it
/// is written stops the output at the record that fails, as a serial run and HTSlib's
/// threads (which push before writing a block) do.
#[derive(Clone, Debug)]
pub(crate) struct PushCheck {
    last: Last,
    /// The targets of the records checked.
    seen: Vec<bool>,
}

impl PushCheck {
    /// Starts checking the records for a new index.
    pub(crate) fn new() -> Self {
        Self {
            last: Last::NEW,
            seen: Vec::new(),
        }
    }

    /// Checks the next record to push to `index`, which must be new and get the records
    /// checked here in the same order, as [`Index::push`] will check it.
    ///
    /// # Errors
    ///
    /// Returns the error [`Index::push`] will return (and reports it as HTSlib does).
    pub(crate) fn check(
        &mut self,
        index: &Index,
        tid: i32,
        beg: i64,
        end: i64,
    ) -> Result<(), PushError> {
        let (beg, end) = if tid < 0 { (-1, 0) } else { (beg, end) };
        index.check_range(tid, beg, end)?;
        let seen = &self.seen;
        if self.last.check(tid, beg, end, || {
            usize::try_from(tid).is_ok_and(|t| seen.get(t).copied().unwrap_or(false))
        })? {
            self.last.tid = tid;
        }
        if let Ok(t) = usize::try_from(tid) {
            if t >= self.seen.len() {
                self.seen.resize(t + 1, false);
            }
            self.seen[t] = true;
            // As `push` clamps it.
            self.last.coor = beg.max(0);
        } else {
            self.last.no_coor = true;
            self.last.coor = beg;
        }
        Ok(())
    }
}

/// A CSI or TBI index (HTSlib's `hts_idx_t`).
#[derive(Clone, Debug)]
pub struct Index {
    fmt: IndexFormat,
    min_shift: u32,
    n_lvls: u32,
    n_bins: u32,
    meta: Vec<u8>,
    bidx: Vec<Option<IndexMap<u32, Bin>>>,
    lidx: Vec<Vec<u64>>,
    n_no_coor: u64,
    z: BuildState,
    tbi_n: i32,
    last_tbi_tid: i32,
}

impl Index {
    /// Starts an index (HTSlib's `hts_idx_init`) with `n` targets, the offset of the first
    /// record, and the binning scheme.
    #[must_use]
    pub fn new(n: usize, fmt: IndexFormat, offset0: u64, min_shift: u32, n_lvls: u32) -> Self {
        Self {
            fmt,
            min_shift,
            n_lvls,
            n_bins: ((1u32 << (3 * n_lvls + 3)) - 1) / 7,
            meta: Vec::new(),
            bidx: vec![None; n],
            lidx: vec![Vec::new(); n],
            n_no_coor: 0,
            z: BuildState {
                last_bin: u32::MAX,
                save_bin: u32::MAX,
                last_coor: Last::NEW.coor,
                last_tid: Last::NEW.tid,
                save_tid: -1,
                finished: false,
                last_off: offset0,
                save_off: offset0,
                off_beg: offset0,
                off_end: offset0,
                n_mapped: 0,
            },
            tbi_n: -1,
            last_tbi_tid: -1,
        }
    }

    /// Returns the index flavor.
    #[inline]
    #[must_use]
    pub fn format(&self) -> IndexFormat {
        self.fmt
    }

    /// Returns the minimum shift.
    #[inline]
    #[must_use]
    pub fn min_shift(&self) -> u32 {
        self.min_shift
    }

    /// Returns the number of levels.
    #[inline]
    #[must_use]
    pub fn n_lvls(&self) -> u32 {
        self.n_lvls
    }

    fn meta_bin(&self) -> u32 {
        self.n_bins + 1
    }

    /// Returns the exclusive maximum position this index covers.
    #[inline]
    #[must_use]
    pub fn max_pos(&self) -> i64 {
        bin_maxpos(self.min_shift, self.n_lvls)
    }

    /// Sets the metadata (the tabix configuration and names).
    pub fn set_meta(&mut self, meta: Vec<u8>) {
        self.meta = meta;
    }

    /// Returns the metadata.
    #[must_use]
    pub fn meta(&self) -> &[u8] {
        &self.meta
    }

    /// Returns the number of targets (HTSlib's `hts_idx_nseq`).
    #[must_use]
    pub fn nseq(&self) -> usize {
        self.bidx.len()
    }

    /// Returns whether a target has any bins.
    #[must_use]
    pub fn has_target(&self, tid: usize) -> bool {
        self.bidx.get(tid).is_some_and(Option::is_some)
    }

    /// Returns the bins of a target: in HTSlib's hash table order once built, in file order
    /// when read.
    #[must_use]
    pub fn bins(&self, tid: usize) -> Option<&IndexMap<u32, Bin>> {
        self.bidx.get(tid)?.as_ref()
    }

    /// Returns the number of records without coordinates.
    #[must_use]
    pub fn n_no_coor(&self) -> u64 {
        self.n_no_coor
    }

    /// Returns a target's mapped and unmapped record counts from the metadata pseudo-bin
    /// (HTSlib's `hts_idx_get_stat`), or `None` when it has none.
    #[must_use]
    pub fn stat(&self, tid: usize) -> Option<(u64, u64)> {
        let bins = self.bidx.get(tid)?.as_ref()?;
        let meta = bins.get(&self.meta_bin())?;
        let c = meta.chunks.get(1)?;
        Some((c.0, c.1))
    }

    /// Maps a VCF contig to its tabix target, adding its name to the metadata on first use
    /// (HTSlib's `hts_idx_tbi_name`).
    pub fn tbi_name(&mut self, tid: i32, name: &[u8]) -> i32 {
        if tid == self.last_tbi_tid || tid < 0 {
            return self.tbi_n;
        }
        self.meta.extend_from_slice(name);
        self.meta.push(0);
        if self.meta.len() >= TBX_CONF_LEN {
            let l = u32::from_le_bytes(self.meta[TBX_L_NM].try_into().expect("4 bytes"));
            let l = l + name.len() as u32 + 1;
            self.meta[TBX_L_NM].copy_from_slice(&l.to_le_bytes());
        }
        self.last_tbi_tid = tid;
        self.tbi_n += 1;
        self.tbi_n
    }

    /// Moves the start of the next record (HTSlib's `hts_idx_amend_last`).
    pub fn amend_last(&mut self, offset: u64) {
        self.z.last_off = offset;
    }

    /// Checks that a record fits the binning scheme (HTSlib's `hts_idx_check_range`), reporting
    /// an error as HTSlib does.
    fn check_range(&self, tid: i32, beg: i64, end: i64) -> Result<(), PushError> {
        let maxpos = self.max_pos();
        if tid < 0 || (beg <= maxpos && end <= maxpos) {
            return Ok(());
        }
        match self.fmt {
            IndexFormat::Csi => log::error(
                "hts_idx_check_range",
                format_args!(
                    "Region {beg}..{end} cannot be stored in a csi index with these parameters. Please use a larger min_shift or depth"
                ),
            ),
            IndexFormat::Tbi => log::error(
                "hts_idx_check_range",
                format_args!(
                    "Region {beg}..{end} cannot be stored in a tbi index. Try using a csi index"
                ),
            ),
        }
        Err(PushError::OutOfRange)
    }

    /// Adds a record (HTSlib's `hts_idx_push`): its target, 0-based half-open interval, and
    /// the file offset just past it. The previous push's offset is taken as its start.
    ///
    /// # Errors
    ///
    /// Returns an error (also reported as HTSlib reports it) for unsorted or invalid records.
    pub fn push(
        &mut self,
        tid: i32,
        mut beg: i64,
        mut end: i64,
        offset: u64,
    ) -> Result<(), PushError> {
        if tid < 0 {
            beg = -1;
            end = 0;
        }
        self.check_range(tid, beg, end)?;
        if tid >= 0 && tid as usize >= self.bidx.len() {
            self.bidx.resize(tid as usize + 1, None);
            self.lidx.resize(tid as usize + 1, Vec::new());
        }
        if self.z.finished {
            return Ok(());
        }
        let last = Last {
            tid: self.z.last_tid,
            coor: self.z.last_coor,
            no_coor: self.n_no_coor > 0,
        };
        if last.check(tid, beg, end, || self.bidx[tid as usize].is_some())? {
            self.z.last_tid = tid;
            self.z.last_bin = u32::MAX;
        }
        if tid >= 0 {
            let t = tid as usize;
            if self.bidx[t].is_none() {
                self.bidx[t] = Some(IndexMap::new());
            }
            if beg < 0 {
                beg = 0;
            }
            if end <= 0 {
                end = 1;
            }
            insert_to_l(&mut self.lidx[t], beg, end, self.z.last_off, self.min_shift);
        } else {
            self.n_no_coor += 1;
        }
        let bin = reg2bin(beg, end, self.min_shift, self.n_lvls);
        if self.z.last_bin != bin {
            if self.z.save_bin != u32::MAX {
                let st = self.z.save_tid as usize;
                insert_to_b(
                    &mut self.bidx[st],
                    self.z.save_bin,
                    self.z.save_off,
                    self.z.last_off,
                );
            }
            if self.z.last_bin == u32::MAX && self.z.save_bin != u32::MAX {
                self.z.off_end = self.z.last_off;
                let st = self.z.save_tid as usize;
                let mb = self.meta_bin();
                insert_to_b(&mut self.bidx[st], mb, self.z.off_beg, self.z.off_end);
                insert_to_b(&mut self.bidx[st], mb, self.z.n_mapped, 0);
                self.z.n_mapped = 0;
                self.z.off_beg = self.z.off_end;
            }
            self.z.save_off = self.z.last_off;
            self.z.save_bin = bin;
            self.z.last_bin = bin;
            self.z.save_tid = tid;
        }
        self.z.n_mapped += 1;
        self.z.last_off = offset;
        self.z.last_coor = beg;
        Ok(())
    }

    /// Finishes building (HTSlib's `hts_idx_finish`): closes the last bin, fills the linear
    /// index, sets bin offsets, and merges small bins into their parents.
    pub fn finish(&mut self, final_offset: u64) {
        if self.z.finished {
            return;
        }
        if self.z.save_tid >= 0 {
            let st = self.z.save_tid as usize;
            let mb = self.meta_bin();
            insert_to_b(
                &mut self.bidx[st],
                self.z.save_bin,
                self.z.save_off,
                final_offset,
            );
            insert_to_b(&mut self.bidx[st], mb, self.z.off_beg, final_offset);
            insert_to_b(&mut self.bidx[st], mb, self.z.n_mapped, 0);
        }
        for i in 0..self.bidx.len() {
            self.update_loff(i);
            if self.fmt == IndexFormat::Csi {
                // CSI files keep no linear index.
                self.lidx[i] = Vec::new();
            }
            let table = self.bidx[i].as_ref().map(bin_table);
            self.compress_binning(i);
            if let (Some(bins), Some(table)) = (&mut self.bidx[i], table) {
                sort_by_table(bins, &table);
            }
        }
        self.z.finished = true;
    }

    /// Fills unset linear-index entries and sets each bin's `loff` (HTSlib's `update_loff`).
    fn update_loff(&mut self, i: usize) {
        let lidx = &mut self.lidx[i];
        for l in (0..lidx.len().saturating_sub(1)).rev() {
            if lidx[l] == u64::MAX {
                lidx[l] = lidx[l + 1];
            }
        }
        let n_bins = self.n_bins;
        let n_lvls = self.n_lvls;
        if let Some(bins) = &mut self.bidx[i] {
            for (&key, bin) in bins.iter_mut() {
                bin.loff = if key < n_bins {
                    let bot = bin_bot(key, n_lvls);
                    if (bot as usize) < lidx.len() {
                        lidx[bot as usize]
                    } else {
                        0
                    }
                } else {
                    0
                };
            }
        }
    }

    /// Merges bins whose chunks span less than a BGZF block into their parents, and merges
    /// chunks starting in the same block (HTSlib's `compress_binning`).
    fn compress_binning(&mut self, i: usize) {
        let n_bins = self.n_bins;
        let n_lvls = self.n_lvls;
        let Some(bins) = &mut self.bidx[i] else {
            return;
        };
        for l in (1..=n_lvls).rev() {
            let start = bin_first(l);
            let keys: Vec<u32> = bins
                .keys()
                .copied()
                .filter(|&k| k < n_bins && k >= start)
                .collect();
            for k in keys {
                let Some(p) = bins.get_mut(&k) else { continue };
                if l < n_lvls && p.chunks.len() > 1 {
                    p.chunks.sort_by_key(|c| c.0);
                }
                let first_u = p.chunks[0].0;
                let last_v = p.chunks[p.chunks.len() - 1].1;
                if (last_v >> 16).wrapping_sub(first_u >> 16) < MIN_MARKER_DIST {
                    let parent = bin_parent(k);
                    if !bins.contains_key(&parent) {
                        continue;
                    }
                    let moved = bins.shift_remove(&k).expect("bin exists").chunks;
                    bins.get_mut(&parent)
                        .expect("parent exists")
                        .chunks
                        .extend(moved);
                }
            }
        }
        if let Some(b0) = bins.get_mut(&0) {
            b0.chunks.sort_by_key(|c| c.0);
        }
        for (&k, p) in bins.iter_mut() {
            if k >= n_bins {
                continue;
            }
            let mut m = 0;
            for l in 1..p.chunks.len() {
                if p.chunks[m].1 >> 16 >= p.chunks[l].0 >> 16 {
                    if p.chunks[m].1 < p.chunks[l].1 {
                        p.chunks[m].1 = p.chunks[l].1;
                    }
                } else {
                    m += 1;
                    p.chunks[m] = p.chunks[l];
                }
            }
            p.chunks.truncate(m + 1);
        }
    }

    // ---- reading and writing ----

    /// Writes the index as `fmt` to a BGZF-compressed file (HTSlib's `hts_idx_save_as`), in
    /// HTSlib's blocks at pasta's default level: the same data as HTSlib's file, block for
    /// block, in zflate's compressed bytes.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn save(&self, path: &Path, fmt: IndexFormat) -> io::Result<()> {
        let file = File::create(path)?;
        let mut data = Vec::new();
        self.write_to(&mut data, fmt)?;
        let mut w = BlockWriter::new(BufWriter::new(file), htslib_level(None), 0);
        w.write_all(&data)?;
        let (mut inner, _) = w.finish()?;
        inner.flush()
    }

    /// Writes the uncompressed index data (HTSlib's `hts_idx_write_out`).
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn write_to<W: Write>(&self, w: &mut W, fmt: IndexFormat) -> io::Result<()> {
        match fmt {
            IndexFormat::Csi => {
                w.write_all(b"CSI\x01")?;
                w.write_all(&(self.min_shift as i32).to_le_bytes())?;
                w.write_all(&(self.n_lvls as i32).to_le_bytes())?;
                w.write_all(&(self.meta.len() as u32).to_le_bytes())?;
                w.write_all(&self.meta)?;
            }
            IndexFormat::Tbi => w.write_all(b"TBI\x01")?,
        }
        let is_vcf = self.meta.len() >= 4 && self.meta[..4] == 2u32.to_le_bytes();
        let n = self.bidx.len();
        let nids = if is_vcf {
            self.bidx.iter().filter(|b| b.is_some()).count()
        } else {
            n
        };
        w.write_all(&(nids as i32).to_le_bytes())?;
        if fmt == IndexFormat::Tbi && !self.meta.is_empty() {
            w.write_all(&self.meta)?;
        }
        for i in 0..n {
            let bins = self.bidx[i].as_ref();
            if nids == n || bins.is_some() {
                w.write_all(&(bins.map_or(0, IndexMap::len) as i32).to_le_bytes())?;
            }
            if let Some(bins) = bins {
                for (&key, p) in bins {
                    w.write_all(&key.to_le_bytes())?;
                    if fmt == IndexFormat::Csi {
                        w.write_all(&p.loff.to_le_bytes())?;
                    }
                    w.write_all(&(p.chunks.len() as i32).to_le_bytes())?;
                    for c in &p.chunks {
                        w.write_all(&c.0.to_le_bytes())?;
                        w.write_all(&c.1.to_le_bytes())?;
                    }
                }
            }
            if fmt != IndexFormat::Csi {
                let lidx = &self.lidx[i];
                w.write_all(&(lidx.len() as i32).to_le_bytes())?;
                for o in lidx {
                    w.write_all(&o.to_le_bytes())?;
                }
            }
        }
        w.write_all(&self.n_no_coor.to_le_bytes())
    }

    /// Reads a CSI or TBI index file (HTSlib's `idx_read`).
    ///
    /// # Errors
    ///
    /// Returns an error for an unreadable or malformed file.
    pub fn read(path: &Path) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut r = bgzf::io::Reader::new(io::BufReader::new(file));
        let mut data = Vec::new();
        r.read_to_end(&mut data)?;
        Self::parse(&data)
    }

    /// Parses uncompressed index data.
    ///
    /// # Errors
    ///
    /// Returns an error for malformed data.
    pub fn parse(data: &[u8]) -> io::Result<Self> {
        let bad = || io::Error::new(io::ErrorKind::InvalidData, "malformed index");
        let mut c = Cursor { data, pos: 0 };
        let magic = c.take(4).ok_or_else(bad)?;
        let mut idx = match magic {
            b"CSI\x01" => {
                let min_shift = c.u32().ok_or_else(bad)?;
                let n_lvls = c.u32().ok_or_else(bad)?;
                let l_meta = c.u32().ok_or_else(bad)? as usize;
                let meta = c.take(l_meta).ok_or_else(bad)?.to_vec();
                let n = c.u32().ok_or_else(bad)?;
                if n > i32::MAX as u32 {
                    return Err(bad());
                }
                let mut idx = Self::new(n as usize, IndexFormat::Csi, 0, min_shift, n_lvls);
                idx.meta = meta;
                idx
            }
            b"TBI\x01" => {
                let n = c.u32().ok_or_else(bad)?;
                if n > i32::MAX as u32 {
                    return Err(bad());
                }
                let conf = c.take(TBX_CONF_LEN).ok_or_else(bad)?.to_vec();
                let l_nm = u32::from_le_bytes(conf[TBX_L_NM].try_into().expect("4 bytes")) as usize;
                let names = c.take(l_nm).ok_or_else(bad)?;
                let mut idx = Self::new(n as usize, IndexFormat::Tbi, 0, TBI_MIN_SHIFT, TBI_N_LVLS);
                idx.meta = conf;
                idx.meta.extend_from_slice(names);
                idx
            }
            _ => return Err(bad()),
        };
        let fmt = idx.fmt;
        for i in 0..idx.bidx.len() {
            let n = c.i32().ok_or_else(bad)?;
            if n < 0 {
                return Err(bad());
            }
            let mut bins = IndexMap::with_capacity(n as usize);
            for _ in 0..n {
                let key = c.u32().ok_or_else(bad)?;
                let loff = if fmt == IndexFormat::Csi {
                    c.u64().ok_or_else(bad)?
                } else {
                    0
                };
                let nc = c.i32().ok_or_else(bad)?;
                if nc < 0 {
                    return Err(bad());
                }
                let mut chunks = Vec::with_capacity(nc as usize);
                for _ in 0..nc {
                    let u = c.u64().ok_or_else(bad)?;
                    let v = c.u64().ok_or_else(bad)?;
                    chunks.push((u, v));
                }
                if bins.insert(key, Bin { loff, chunks }).is_some() {
                    return Err(bad());
                }
            }
            idx.bidx[i] = Some(bins);
            if fmt != IndexFormat::Csi {
                let nl = c.i32().ok_or_else(bad)?;
                if nl < 0 {
                    return Err(bad());
                }
                let mut l = Vec::with_capacity(nl as usize);
                for _ in 0..nl {
                    l.push(c.u64().ok_or_else(bad)?);
                }
                let mut k = 0;
                while k < l.len() && l[k] == 0 {
                    k += 1;
                }
                let mut j = l.len().saturating_sub(1);
                while j > k {
                    if l[j - 1] == 0 {
                        l[j - 1] = l[j];
                    }
                    j -= 1;
                }
                idx.lidx[i] = l;
                idx.update_loff(i);
            }
        }
        idx.n_no_coor = c.u64().unwrap_or(0);
        idx.z.finished = true;
        Ok(idx)
    }

    // ---- queries ----

    /// Plans a query of `[beg, end)` on target `tid` (HTSlib's `hts_itr_query`): the merged
    /// chunks to read. `None` when `end < beg`.
    #[must_use]
    pub fn query(&self, tid: usize, beg: i64, end: i64) -> Option<QueryPlan> {
        let mut plan = QueryPlan {
            tid,
            beg,
            end,
            chunks: Vec::new(),
        };
        let Some(Some(bins)) = self.bidx.get(tid) else {
            return Some(plan);
        };
        let beg = beg.max(0);
        plan.beg = beg;
        if end < beg {
            return None;
        }
        if bins.is_empty() || beg >= self.max_pos() {
            return Some(plan);
        }
        let min_off = self.min_offset(tid, bins, beg);
        let max_off = self.max_offset(bins, end);
        let off = self.overlapping_chunks(tid, bins, beg, end, min_off..max_off);
        if !off.is_empty() {
            plan.chunks = merge_chunks(off);
        }
        Some(plan)
    }

    /// The smallest file offset a query starting at `beg` can need, from the bins and the
    /// linear index (part of HTSlib's `hts_itr_query`).
    fn min_offset(&self, tid: usize, bins: &IndexMap<u32, Bin>, beg: i64) -> u64 {
        let meta = bins.get(&self.meta_bin());
        let unmapped = meta.and_then(|m| m.chunks.get(1)).map_or(1, |c| c.1);
        let rel_off = (beg >> self.min_shift) as usize;
        let mut bin = bin_first(self.n_lvls) + rel_off as u32;
        let mut found = None;
        loop {
            if let Some(b) = bins.get(&bin) {
                found = Some(b);
                break;
            }
            let first = (bin_parent(bin) << 3) + 1;
            if bin > first {
                bin -= 1;
            } else {
                bin = bin_parent(bin);
            }
            if bin == 0 {
                break;
            }
        }
        if bin == 0 {
            found = bins.get(&0);
        }
        let mut min_off = found.map_or(0, |b| b.loff);
        let lidx = &self.lidx[tid];
        if !lidx.is_empty() && rel_off < lidx.len() {
            if min_off < lidx[rel_off] {
                min_off = lidx[rel_off];
            }
            if unmapped != 0 {
                let mut tmp_off: isize = rel_off as isize - 1;
                while tmp_off >= 0 {
                    if lidx[tmp_off as usize] < min_off {
                        min_off = lidx[tmp_off as usize];
                        break;
                    }
                    tmp_off -= 1;
                }
                if let Some(b) = found
                    && (min_off < b.loff || tmp_off < 0)
                {
                    min_off = b.loff;
                }
            }
        } else if unmapped != 0
            && let Some(b) = found
        {
            min_off = b.loff;
        }
        min_off
    }

    /// The file offset a query ending at `end` need not read past: the start of the first
    /// bin after it (part of HTSlib's `hts_itr_query`).
    fn max_offset(&self, bins: &IndexMap<u32, Bin>, end: i64) -> u64 {
        if end > self.max_pos() {
            return u64::MAX;
        }
        let mut bin = bin_first(self.n_lvls) + ((end - 1) >> self.min_shift) as u32 + 1;
        if bin >= self.n_bins {
            bin = 0;
        }
        loop {
            while bin % 8 == 1 {
                bin = bin_parent(bin);
            }
            if bin == 0 {
                return u64::MAX;
            }
            if let Some(b) = bins.get(&bin)
                && !b.chunks.is_empty()
            {
                return b.chunks[0].0;
            }
            bin += 1;
        }
    }

    /// The chunks of the bins overlapping `[beg, end)`, clipped to `offsets`, as
    /// `(begin, end, order)` with `order` keeping equal begins in bin order.
    fn overlapping_chunks(
        &self,
        tid: usize,
        bins: &IndexMap<u32, Bin>,
        beg: i64,
        end: i64,
        offsets: Range<u64>,
    ) -> Vec<(u64, u64, u64)> {
        let Range {
            start: min_off,
            end: max_off,
        } = offsets;
        let mut off: Vec<(u64, u64, u64)> = Vec::new();
        let s0 = self.min_shift + self.n_lvls * 3;
        let end_cl = end.min(1i64 << s0);
        if beg >= end_cl {
            return off;
        }
        let e1 = end_cl - 1;
        let mut add_bin = |bin: u32| {
            if let Some(p) = bins.get(&bin) {
                for (j, c) in p.chunks.iter().enumerate() {
                    if c.1 > min_off && c.0 < max_off {
                        off.push((
                            c.0.max(min_off),
                            c.1.min(max_off),
                            ((tid as u64) << 32) | j as u64,
                        ));
                    }
                }
            }
        };
        // The bin range of each level (HTSlib's `reg2bins`).
        let ranges = (0..=self.n_lvls).map(|l| {
            let s = s0.saturating_sub(3 * l);
            let t = i64::from(bin_first(l));
            (t + (beg >> s), t + (e1 >> s))
        });
        let reg_bin_count: i64 = ranges.clone().map(|(b, e)| e - b + 1).sum();
        if reg_bin_count > bins.len() as i64 {
            // A region wider than the index: test each bin instead (`reg2bins_wide`), in
            // the same ascending order.
            let mut wanted: Vec<u32> = bins
                .keys()
                .copied()
                .filter(|&bin| {
                    let l = bin_level(bin);
                    l <= self.n_lvls && {
                        let s = s0 - 3 * l;
                        let t = i64::from(bin_first(l));
                        (t + (beg >> s)..=t + (e1 >> s)).contains(&i64::from(bin))
                    }
                })
                .collect();
            wanted.sort_unstable();
            for bin in wanted {
                add_bin(bin);
            }
        } else {
            for (b, e) in ranges {
                for i in b..=e {
                    add_bin(i as u32);
                }
            }
        }
        off
    }

    /// Returns the offset of the first record (HTSlib's `hts_itr_off(HTS_IDX_START)`), or
    /// `None` when there are no records.
    #[must_use]
    pub fn first_offset(&self) -> Option<u64> {
        let mb = self.meta_bin();
        let off = self
            .bidx
            .iter()
            .flatten()
            .filter_map(|b| b.get(&mb).and_then(|m| m.chunks.first()).map(|c| c.0))
            .min();
        off.or((self.n_no_coor > 0).then_some(0))
    }
}

/// Sorts a query's chunks and merges them (the end of HTSlib's `hts_itr_query`): drops the
/// chunks contained in others, trims overlaps, and joins the chunks that meet in one BGZF
/// block.
fn merge_chunks(mut off: Vec<(u64, u64, u64)>) -> Vec<(u64, u64)> {
    off.sort_by_key(|a| (a.0, a.2));
    let mut l = 0;
    for i in 1..off.len() {
        if off[l].1 < off[i].1 {
            l += 1;
            off[l] = off[i];
        }
    }
    off.truncate(l + 1);
    for i in 1..off.len() {
        if off[i - 1].1 >= off[i].0 {
            off[i - 1].1 = off[i].0;
        }
    }
    let mut l = 0;
    for i in 1..off.len() {
        if off[l].1 >> 16 == off[i].0 >> 16 {
            off[l].1 = off[i].1;
        } else {
            l += 1;
            off[l] = off[i];
        }
    }
    off.truncate(l + 1);
    off.into_iter().map(|(u, v, _)| (u, v)).collect()
}

/// The chunks to read for a region query (HTSlib's `hts_itr_t` offsets).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QueryPlan {
    /// The target.
    pub tid: usize,
    /// The 0-based start of the region.
    pub beg: i64,
    /// The 0-based exclusive end of the region.
    pub end: i64,
    /// The `(begin, end)` virtual offsets to read, in order.
    pub chunks: Vec<(u64, u64)>,
}

impl QueryPlan {
    /// Returns whether a record with `[beg, end)` on `tid` overlaps the region (the record
    /// filter of HTSlib's `hts_itr_next`).
    #[inline]
    #[must_use]
    pub fn overlaps(&self, beg: i64, end: i64) -> bool {
        end > self.beg && self.end > beg
    }

    /// Returns whether a record at `beg` on `tid` is past the region, ending the query.
    #[inline]
    #[must_use]
    pub fn is_past(&self, tid: i64, beg: i64) -> bool {
        tid != self.tid as i64 || beg >= self.end
    }
}

/// HTSlib's hash table of a target's bins once every record is in. A bin is put each time a
/// run of records in it ends, so the bins enter the table in the order they were first
/// inserted here; the puts of bins already present add nothing, and the table grows at the
/// next put after it fills either way. The last bin to enter is the pseudo-bin, put twice when
/// the target ends: that second put can still grow the table.
fn bin_table(bins: &IndexMap<u32, Bin>) -> BinTable {
    let mut table = BinTable::default();
    for &key in bins.keys() {
        table.put(key);
    }
    if let Some(&last) = bins.last().map(|(key, _)| key) {
        table.put(last);
    }
    table
}

/// Puts `bins` in the bucket order of `table` (which holds every one of them, and those that
/// `compress_binning` merged into their parents as `kh_del` leaves them, as deleted buckets).
fn sort_by_table(bins: &mut IndexMap<u32, Bin>, table: &BinTable) {
    let mut sorted = IndexMap::with_capacity(bins.len());
    for key in table.keys() {
        if let Some(bin) = bins.swap_remove(&key) {
            sorted.insert(key, bin);
        }
    }
    debug_assert!(bins.is_empty(), "every bin is in the table");
    // Keep any bin the table lacks rather than drop it: the order would differ, the index not.
    sorted.extend(bins.drain(..));
    *bins = sorted;
}

fn insert_to_b(bins: &mut Option<IndexMap<u32, Bin>>, bin: u32, beg: u64, end: u64) {
    bins.get_or_insert_with(IndexMap::new)
        .entry(bin)
        .or_default()
        .chunks
        .push((beg, end));
}

fn insert_to_l(l: &mut Vec<u64>, beg: i64, end: i64, offset: u64, min_shift: u32) {
    let b = (beg >> min_shift) as usize;
    let e = ((end - 1) >> min_shift) as usize;
    if l.len() < e + 1 {
        l.resize(e + 1, u64::MAX);
    }
    for o in &mut l[b..=e] {
        if *o == u64::MAX {
            *o = offset;
        }
    }
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        let s = self.data.get(self.pos..self.pos.checked_add(n)?)?;
        self.pos += n;
        Some(s)
    }

    fn u32(&mut self) -> Option<u32> {
        self.take(4)
            .map(|b| u32::from_le_bytes(b.try_into().expect("4 bytes")))
    }

    fn i32(&mut self) -> Option<i32> {
        self.take(4)
            .map(|b| i32::from_le_bytes(b.try_into().expect("4 bytes")))
    }

    fn u64(&mut self) -> Option<u64> {
        self.take(8)
            .map(|b| u64::from_le_bytes(b.try_into().expect("8 bytes")))
    }
}

/// Warns, as HTSlib's `hts_idx_load3` does before loading an index, when the index file is
/// older than its data file (by whole seconds of modification time).
pub fn warn_if_older(data: &Path, index: &Path) {
    let mtime = |p: &Path| {
        std::fs::metadata(p)
            .and_then(|m| m.modified())
            .ok()
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|d| d.as_secs())
    };
    if let (Some(data_time), Some(index_time)) = (mtime(data), mtime(index))
        && index_time < data_time
    {
        log::warn(
            "hts_idx_load3",
            format_args!(
                "The index file is older than the data file: {}",
                index.display()
            ),
        );
    }
}

/// Reports an index that could not be loaded, as HTSlib's `hts_idx_load3` does: with the
/// system error, or `Invalid argument` for a file that is not an index.
pub fn report_load_error(index: &Path, e: &io::Error) {
    let reason = if e.raw_os_error().is_some() {
        let msg = e.to_string();
        msg.split(" (os error").next().unwrap_or(&msg).to_owned()
    } else {
        "Invalid argument".to_owned()
    };
    log::error(
        "hts_idx_load3",
        format_args!(
            "Could not load local index file '{}' : {reason}",
            index.display()
        ),
    );
}

/// Finds the index file of a data file as HTSlib does (`hts_idx_check_local`): `FILE.csi`,
/// then the extension replaced by `.csi`, then (for tabix) `FILE.tbi` and the replaced form.
/// A name containing `##idx##` gives the index explicitly.
#[must_use]
pub fn locate(path: &str, fmt: IndexFormat) -> Option<std::path::PathBuf> {
    if let (_, Some(index)) = split_idx_name(path) {
        return Some(index.into());
    }
    let try_ext = |ext: &str| -> Option<std::path::PathBuf> {
        let full = format!("{path}{ext}");
        if Path::new(&full).exists() {
            return Some(full.into());
        }
        if let Some(dot) = path.rfind('.').filter(|&d| d > 0) {
            let replaced = format!("{}{ext}", &path[..dot]);
            if Path::new(&replaced).exists() {
                return Some(replaced.into());
            }
        }
        None
    };
    try_ext(".csi").or_else(|| (fmt == IndexFormat::Tbi).then(|| try_ext(".tbi")).flatten())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checking_ahead_fails_on_the_record_push_fails_on() {
        let max = bin_maxpos(14, 5);
        // (tid, beg, end) sequences, each ending in the record push fails on, if any.
        type Records<'a> = &'a [(i32, i64, i64)];
        let cases: &[(Records<'_>, Option<PushError>)] = &[
            (&[(0, 5, 6), (0, 5, 9), (1, 0, 1), (-1, 7, 3)], None),
            (&[(0, 10, 11), (0, 9, 10)], Some(PushError::Unsorted)),
            (&[(0, -5, 1), (0, 0, 1)], None),
            (
                &[(0, -5, 1), (0, 0, 1), (0, -1, 1)],
                Some(PushError::Unsorted),
            ),
            (
                &[(0, 1, 2), (1, 1, 2), (0, 3, 4)],
                Some(PushError::NotContiguous),
            ),
            (
                &[(0, 1, 2), (-1, 0, 0), (1, 1, 2)],
                Some(PushError::NoCoorNotLast),
            ),
            (&[(0, 1, 2), (0, 5, 4)], Some(PushError::InvalidInterval)),
            (&[(0, 1, 2), (0, max, max + 1)], Some(PushError::OutOfRange)),
            (
                &[(2, 1, 2), (5, 1, 2), (2, 9, 10)],
                Some(PushError::NotContiguous),
            ),
        ];
        for (records, want) in cases {
            let mut pushed = Index::new(1, IndexFormat::Csi, 0, 14, 5);
            let mut checked = PushCheck::new();
            for (i, &(tid, beg, end)) in records.iter().enumerate() {
                let got = checked.check(&pushed, tid, beg, end);
                let push = pushed.push(tid, beg, end, (i as u64 + 1) << 16);
                assert_eq!(got, push, "{records:?}, record {i}");
                if got.is_err() {
                    assert_eq!(i + 1, records.len(), "{records:?}: failed early");
                }
            }
            // Checked ahead of every push, as with compression threads.
            let index = Index::new(1, IndexFormat::Csi, 0, 14, 5);
            let mut check = PushCheck::new();
            let result = records
                .iter()
                .try_for_each(|&(tid, beg, end)| check.check(&index, tid, beg, end));
            assert_eq!(result.err(), *want, "{records:?}");
        }
    }

    #[test]
    fn reg2bin_matches_htslib() {
        assert_eq!(reg2bin(0, 1, 14, 5), 4681);
        assert_eq!(reg2bin(0, 1 << 29, 14, 5), 0);
        assert_eq!(reg2bin(16384, 16385, 14, 5), 4682);
    }

    #[test]
    fn build_and_query() {
        let mut idx = Index::new(1, IndexFormat::Csi, 100, 14, 5);
        idx.push(0, 10, 20, 200).unwrap();
        idx.push(0, 100_000, 100_010, 300 << 16).unwrap();
        idx.finish(400 << 16);
        assert_eq!(idx.stat(0), Some((2, 0)));
        let plan = idx.query(0, 0, 50).unwrap();
        assert!(!plan.chunks.is_empty());
        assert_eq!(plan.chunks[0].0, 100);
        let mut buf = Vec::new();
        idx.write_to(&mut buf, IndexFormat::Csi).unwrap();
        let back = Index::parse(&buf).unwrap();
        assert_eq!(back.stat(0), Some((2, 0)));
        assert_eq!(back.query(0, 0, 50).unwrap().chunks, plan.chunks);
    }

    #[test]
    fn bins_are_in_htslib_hash_order() {
        // One record in each of four bottom-level bins of a second target, a block (64 KiB of
        // file) apart so no bin merges into its parent: bcftools writes their bins as
        // 4688, 37450 (the pseudo-bin), 4685, 4686, 4687.
        let mut idx = Index::new(2, IndexFormat::Tbi, 0, 14, 5);
        idx.push(0, 0, 1, 1 << 32).unwrap();
        for (i, bin) in (4685..=4688).enumerate() {
            let beg = i64::from(bin - 4681) << 14;
            idx.push(1, beg, beg + 1, (i as u64 + 2) << 32).unwrap();
        }
        idx.finish(6 << 32);
        let keys: Vec<u32> = idx.bins(1).unwrap().keys().copied().collect();
        assert_eq!(keys, [4688, 37450, 4685, 4686, 4687]);
        let keys: Vec<u32> = idx.bins(0).unwrap().keys().copied().collect();
        assert_eq!(keys, [4681, 37450]);
    }
}
