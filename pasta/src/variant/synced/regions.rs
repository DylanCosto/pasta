//! Region lists for the synced reader (HTSlib's `bcf_sr_regions_t`).

use std::{
    io,
    path::{Path, PathBuf},
};

use rustc_hash::FxHashMap;

use super::SyncedError;
use crate::variant::{
    Record, Unpack, VariantType,
    index::{
        Index, IndexFormat, data_file, locate, split_idx_name,
        tabix::{PRESET_UCSC, Tabix},
    },
    io::{Query, ReadError, Reader},
    log,
    number::{ParseDecimalFlags, parse_decimal},
    readlist,
};

/// The largest coordinate a default CSI index can hold (HTSlib's `MAX_CSI_COOR`).
pub const MAX_CSI_COOR: i64 = (1i64 << (14 + 30)) - 1;

/// The end of a whole-sequence tabix query (HTSlib's `HTS_POS_MAX`).
const WHOLE_SEQUENCE: i64 = i64::MAX;

/// When a record counts as overlapping a region (HTSlib's `BCF_SR_*_OVERLAP` values).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub enum Overlap {
    /// The record's POS is in the region (0).
    #[default]
    Pos,
    /// Any base of the record (POS to POS+rlen-1) is in the region (1).
    Record,
    /// Any base of the variant, after the shared prefix of REF and ALT, is in the region (2).
    Variant,
}

impl Overlap {
    /// Parses bcftools' `--regions-overlap`/`--targets-overlap` (`parse_overlap_option`).
    #[must_use]
    pub fn parse(s: &str) -> Option<Self> {
        if s.eq_ignore_ascii_case("pos") || s == "0" {
            Some(Self::Pos)
        } else if s.eq_ignore_ascii_case("record") || s == "1" {
            Some(Self::Record)
        } else if s.eq_ignore_ascii_case("variant") || s == "2" {
            Some(Self::Variant)
        } else {
            None
        }
    }

    /// Returns the 0-based inclusive interval a record covers under this rule.
    pub fn interval(self, rec: &mut Record) -> (i64, i64) {
        match self {
            Self::Pos => (rec.pos(), rec.pos()),
            Self::Record => (rec.pos(), rec.pos() + rec.rlen() - 1),
            Self::Variant => variant_boundaries(rec),
        }
    }
}

/// The variant's bases after the prefix all ALT alleles share with REF (HTSlib's
/// `_set_variant_boundaries`).
fn variant_boundaries(rec: &mut Record) -> (i64, i64) {
    let mut off = 0i64;
    if rec.allele_count() > 0 {
        off = rec.rlen();
        rec.unpack(Unpack::STR);
        let r = rec.allele(0);
        for i in 1..rec.allele_count() {
            let a = rec.allele(i);
            let j = r.iter().zip(a).take_while(|(x, y)| x == y).count() as i64;
            if off > j {
                off = j;
            }
            if off == 0 {
                break;
            }
        }
    }
    (rec.pos() + off, rec.pos() + rec.rlen() - 1)
}

#[derive(Clone, Debug, Default)]
struct SeqRegions {
    /// 0-based inclusive regions; merged-away ones have `start > end`.
    regs: Vec<(i64, i64)>,
    creg: isize,
}

/// Where reading a tabix-indexed region file goes on (HTSlib's `itr`, or none).
#[derive(Debug)]
enum Cursor {
    /// Back at the first line, to be read in file order (after [`Regions::rewind`]).
    Start,
    /// Reading the file in order, as no sequence was sought yet.
    InOrder,
    /// Reading one sequence through the index, after a seek.
    Query(Query),
}

/// Regions streamed from a tabix-indexed file (HTSlib's `bcf_sr_regions_t` with `tbx`): read
/// in file order until a seek, then one sequence at a time through the index. Only the last
/// line read is kept.
struct FileRegions {
    name: String,
    file: Reader,
    tbx: Tabix,
    cursor: Cursor,
    /// The last line read (HTSlib's `line`).
    line: Vec<u8>,
    /// Whether `line` holds the current region, which is not so after a seek or a failed read.
    has_line: bool,
    ichr: usize,
    ifrom: usize,
    ito: usize,
    /// Whether the index has the BED preset (0-based starts).
    bed: bool,
}

impl std::fmt::Debug for FileRegions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FileRegions")
            .field("name", &self.name)
            .field("cursor", &self.cursor)
            .field("has_line", &self.has_line)
            .finish_non_exhaustive()
    }
}

impl FileRegions {
    /// Reads the next line into `line`, up to any NUL byte (HTSlib's lines are C strings);
    /// `false` at the end of the file or of the sequence sought.
    fn read_line(&mut self) -> io::Result<bool> {
        if matches!(self.cursor, Cursor::Start) {
            self.file.seek(0)?;
            self.cursor = Cursor::InOrder;
        }
        let got = match &mut self.cursor {
            Cursor::Query(q) => self
                .file
                .query_next_line(q, &mut self.tbx, &mut self.line)?,
            Cursor::Start | Cursor::InOrder => self.file.read_line(&mut self.line)?,
        };
        if let Some(nul) = memchr::memchr(0, &self.line) {
            self.line.truncate(nul);
        }
        Ok(got)
    }
}

/// The type of the alleles listed for the current target (HTSlib's `als_type`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
enum AllelesType {
    #[default]
    Snp,
    Indel,
}

/// A region passed over without any record in it, reported when missed-region tracking is on
/// (HTSlib's `missed_reg_handler`).
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct MissedRegion {
    /// The sequence name.
    pub seq: Vec<u8>,
    /// The 0-based inclusive start.
    pub start: i64,
    /// The 0-based inclusive end.
    pub end: i64,
    /// The region's line, for regions read from a tabix-indexed file.
    pub line: Option<Vec<u8>>,
}

/// A list of regions (HTSlib's `bcf_sr_regions_t`).
///
/// Regions come from a string (`chr`, `chr:pos`, `chr:beg-end`, `chr:beg-`, `{name}:...`,
/// comma-separated) or a file (CHROM, BEG, END columns, or CHROM, POS; BED by extension). They
/// are sorted and merged per sequence, and sequences keep the order in which they were first
/// listed. A file with a tabix index is instead streamed as HTSlib streams it: in file order
/// without merging until a seek, then a sequence at a time through the index, with its
/// sequences in index order. Reading such a file can fail at any line, so the methods that
/// read return a [`RegionsError`], reported as HTSlib reports it.
///
/// Used as the synced reader's regions, they are visited in order by [`Self::next_region`];
/// used as targets, records are tested in order by [`Self::overlaps`], which advances through
/// the regions.
#[derive(Debug)]
pub struct Regions {
    seq_names: Vec<Vec<u8>>,
    seq_hash: FxHashMap<Vec<u8>, usize>,
    seqs: Vec<SeqRegions>,
    file: Option<Box<FileRegions>>,
    pub(super) iseq: isize,
    pub(super) start: i64,
    pub(super) end: i64,
    prev_seq: isize,
    prev_start: i64,
    pub(super) prev_end: i64,
    als: Vec<Vec<u8>>,
    als_parsed: bool,
    als_type: AllelesType,
    missed: Option<Vec<MissedRegion>>,
    /// The overlap rule for records. The synced reader uses its own setting
    /// ([`super::SyncedReaders::set_regions_overlap`]), which it copies here, as HTSlib does.
    pub overlap: Overlap,
}

/// An error reading a region list.
#[derive(Debug)]
pub struct RegionsError {
    msg: String,
    /// Whether a line of a tabix-indexed file could not be read, rather than being malformed.
    read_failed: bool,
}

impl RegionsError {
    fn new(msg: String) -> Self {
        Self {
            msg,
            read_failed: false,
        }
    }

    /// Returns whether reading a line of a tabix-indexed file failed (a truncated or corrupt
    /// BGZF block, say), rather than the list or a line being malformed.
    #[must_use]
    pub fn is_read_failure(&self) -> bool {
        self.read_failed
    }
}

impl std::fmt::Display for RegionsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.msg)
    }
}

impl std::error::Error for RegionsError {}

impl Default for Regions {
    fn default() -> Self {
        Self::empty()
    }
}

impl Regions {
    /// Returns an empty list (HTSlib's `bcf_sr_regions_alloc`).
    #[must_use]
    pub fn empty() -> Self {
        Self {
            seq_names: Vec::new(),
            seq_hash: FxHashMap::default(),
            seqs: Vec::new(),
            file: None,
            iseq: 0,
            start: -1,
            end: -1,
            prev_seq: -1,
            prev_start: -1,
            prev_end: -1,
            als: Vec::new(),
            als_parsed: false,
            als_type: AllelesType::Snp,
            missed: None,
            overlap: Overlap::Pos,
        }
    }

    /// Adds a region with 1-based inclusive coordinates, or the whole sequence for
    /// `(-1, -1)` (HTSlib's `_regions_add`). Call [`Self::sort_and_merge`] afterwards.
    pub fn add(&mut self, chr: &[u8], start: i64, end: i64) {
        let (start, end) = if start == -1 && end == -1 {
            (0, MAX_CSI_COOR - 1)
        } else {
            (start - 1, end - 1)
        };
        let iseq = self.seq_index_or_add(chr);
        self.seqs[iseq].regs.push((start, end));
    }

    fn seq_index_or_add(&mut self, chr: &[u8]) -> usize {
        if let Some(&i) = self.seq_hash.get(chr) {
            return i;
        }
        let i = self.seq_names.len();
        self.seq_names.push(chr.to_vec());
        self.seq_hash.insert(chr.to_vec(), i);
        self.seqs.push(SeqRegions {
            regs: Vec::new(),
            creg: -1,
        });
        i
    }

    /// Sorts each sequence's regions and merges overlapping ones (HTSlib's
    /// `_regions_sort_and_merge`).
    pub fn sort_and_merge(&mut self) {
        for s in &mut self.seqs {
            s.regs.sort_unstable();
            let n = s.regs.len();
            let mut i = 0;
            while i < n {
                let mut j = i + 1;
                while j < n && s.regs[i].1 >= s.regs[j].0 {
                    if s.regs[i].1 < s.regs[j].1 {
                        s.regs[i].1 = s.regs[j].1;
                    }
                    s.regs[j] = (1, 0);
                    j += 1;
                }
                i = j;
            }
        }
    }

    /// Parses a comma-separated region list (HTSlib's `_regions_init_string`), then sorts and
    /// merges it.
    ///
    /// # Errors
    ///
    /// Returns an error (also reported as HTSlib reports it) for a malformed list.
    pub fn from_string(s: &[u8]) -> Result<Self, RegionsError> {
        let mut reg = Self::empty();
        let fail = |msg: String| {
            log::error("_regions_init_string", format_args!("{msg}"));
            RegionsError::new(msg)
        };
        let full = String::from_utf8_lossy(s).into_owned();
        let at = |i: usize| s.get(i).copied().unwrap_or(0);
        let mut sp = 0usize;
        let mut ep = 0usize;
        loop {
            let name: Vec<u8> = if at(ep) == b'{' {
                while at(ep) != 0 && at(ep) != b'}' {
                    ep += 1;
                }
                if at(ep) == 0 {
                    return Err(fail(format!(
                        "Could not parse the region, mismatching braces in: \"{full}\""
                    )));
                }
                ep += 1;
                s[sp + 1..ep - 1].to_vec()
            } else {
                while at(ep) != 0 && at(ep) != b',' && at(ep) != b':' {
                    ep += 1;
                }
                s[sp..ep].to_vec()
            };
            if at(ep) == b':' {
                sp = ep + 1;
                let (from, e) = parse_decimal(s, sp, ParseDecimalFlags::default());
                ep = e;
                if sp == ep {
                    return Err(fail(format!("Could not parse the region(s): {full}")));
                }
                if at(ep) == 0 || at(ep) == b',' {
                    reg.add(&name, from, from);
                    if at(ep) == 0 {
                        break;
                    }
                    // HTSlib leaves sp at the comma; the next name starts after it.
                    ep += 1;
                    sp = ep;
                    continue;
                }
                if at(ep) != b'-' {
                    return Err(fail(format!("Could not parse the region(s): {full}")));
                }
                ep += 1;
                sp = ep;
                let (mut to, e) = parse_decimal(s, sp, ParseDecimalFlags::default());
                ep = e;
                if at(ep) != 0 && at(ep) != b',' {
                    return Err(fail(format!("Could not parse the region(s): {full}")));
                }
                if sp == ep {
                    to = MAX_CSI_COOR - 1;
                }
                reg.add(&name, from, to);
                if at(ep) == 0 {
                    break;
                }
                ep += 1;
                sp = ep;
            } else if at(ep) == 0 || at(ep) == b',' {
                if !name.is_empty() {
                    reg.add(&name, -1, -1);
                }
                if at(ep) == 0 {
                    break;
                }
                ep += 1;
                sp = ep;
            } else {
                return Err(fail(format!("Could not parse the region(s): {full}")));
            }
        }
        reg.sort_and_merge();
        Ok(reg)
    }

    /// Reads regions from a file (HTSlib's `bcf_sr_regions_init` with a file and the default
    /// columns): tab-separated CHROM, BEG, END (1-based inclusive; 0-based half-open for
    /// `.bed`/`.bed.gz`), or CHROM and POS. A file with a tabix index is streamed in file order
    /// without merging, as HTSlib streams it.
    ///
    /// # Errors
    ///
    /// Returns an error for an unreadable, malformed, or empty file.
    pub fn from_file(path: &Path) -> Result<Self, RegionsError> {
        Self::from_file_columns(path, 0, 1, -2)
    }

    /// Reads regions from a file using the given 0-based columns (HTSlib's
    /// `bcf_sr_regions_init(fname, 1, ichr, ifrom, ito)`). A negative `to_col` is optional:
    /// column `-to_col` is used when present on the first line, else the region is one base.
    /// VCF files use CHROM and POS. A tabix-indexed file uses its index's columns and is only
    /// opened here: its lines are read as the regions are visited. `FILE##idx##INDEX` names
    /// the index; when it cannot be loaded, `FILE` is read whole.
    ///
    /// # Errors
    ///
    /// Returns an error for an unreadable file, or a malformed or empty file without an index.
    pub fn from_file_columns(
        path: &Path,
        chr_col: usize,
        from_col: usize,
        to_col: isize,
    ) -> Result<Self, RegionsError> {
        let name = path.to_string_lossy().into_owned();
        // HTSlib ignores an index it cannot load (`HTS_IDX_SILENT_FAIL`).
        let tbx = locate(&name, IndexFormat::Tbi)
            .and_then(|p| Index::read(&p).ok())
            .and_then(Tabix::from_index);
        match tbx {
            Some(tbx) => Self::from_tabix_file(name, tbx),
            None => Self::from_plain_file(path, &name, chr_col, from_col, to_col),
        }
    }

    /// Opens a tabix-indexed region file for streaming, with the sequences of its index.
    fn from_tabix_file(name: String, tbx: Tabix) -> Result<Self, RegionsError> {
        let file = Reader::open_text(&name, 0).map_err(|e| match e {
            ReadError::Io(e) => open_failed(&name, &e),
            e => open_failed(
                &name,
                &io::Error::new(io::ErrorKind::InvalidData, e.to_string()),
            ),
        })?;
        let mut reg = Self::empty();
        for seq in tbx.names() {
            reg.seq_index_or_add(seq);
        }
        let ifrom = usize::try_from(tbx.conf.bc - 1).unwrap_or(0);
        reg.file = Some(Box::new(FileRegions {
            ichr: usize::try_from(tbx.conf.sc - 1).unwrap_or(0),
            ifrom,
            ito: usize::try_from(tbx.conf.ec - 1).unwrap_or(ifrom),
            bed: tbx.conf.preset == PRESET_UCSC,
            name,
            file,
            tbx,
            cursor: Cursor::InOrder,
            line: Vec::new(),
            has_line: false,
        }));
        Ok(reg)
    }

    /// Reads a region file without an index line by line, keeping only its regions. The
    /// `.bed` test looks at the data file's name (HTSlib's at the whole `##idx##` name).
    fn from_plain_file(
        path: &Path,
        name: &str,
        chr_col: usize,
        from_col: usize,
        to_col: isize,
    ) -> Result<Self, RegionsError> {
        let data = data_file(path);
        let lower = data.to_string_lossy().to_ascii_lowercase();
        let is_bed = lower.ends_with(".bed") || lower.ends_with(".bed.gz");
        let mut reader = readlist::open_text(data).map_err(|e| open_failed(name, &e))?;
        let mut reg = Self::empty();
        let (ichr, ifrom) = (chr_col, from_col);
        let mut ito = to_col;
        let mut line = Vec::new();
        let mut iline = 0usize;
        // HTSlib stops reading a region file at its first empty line.
        while readlist::read_line(&mut reader, &mut line).map_err(|e| open_failed(name, &e))?
            && !line.is_empty()
        {
            iline += 1;
            if iline == 1 && line.starts_with(b"##fileformat=VCF") {
                ito = 1;
            }
            let mut parsed = parse_region_line(&line, ichr, ifrom, ito.unsigned_abs());
            if parsed.is_none() {
                if ito < 0 {
                    parsed = parse_region_line(&line, ichr, ifrom, ifrom);
                }
                if parsed.is_none() {
                    log::error(
                        "bcf_sr_regions_init",
                        format_args!(
                            "Could not parse {iline}-th line of file {name}, using the columns {},{}[,{}]",
                            ichr + 1,
                            ifrom + 1,
                            ito + 1
                        ),
                    );
                    return Err(RegionsError::new(format!("could not parse {name}")));
                }
                ito = ifrom as isize;
            } else if ito < 0 {
                ito = ito.abs();
            }
            let Some(Some((chr, mut from, to))) = parsed else {
                continue;
            };
            if is_bed {
                from += 1;
            }
            reg.add(chr, from, to);
        }
        if reg.seqs.is_empty() {
            return Err(RegionsError::new(format!("no regions in {name}")));
        }
        reg.sort_and_merge();
        Ok(reg)
    }

    /// Reads regions from a string or, with [`readlist::ListSource::File`], a file.
    ///
    /// # Errors
    ///
    /// Returns an error for a malformed list or file.
    pub fn new(spec: &[u8], source: readlist::ListSource) -> Result<Self, RegionsError> {
        match source {
            readlist::ListSource::Inline => Self::from_string(spec),
            readlist::ListSource::File => {
                Self::from_file(&PathBuf::from(String::from_utf8_lossy(spec).into_owned()))
            }
        }
    }

    /// Returns the sequence names in the order they were first seen.
    #[must_use]
    pub fn seq_names(&self) -> &[Vec<u8>] {
        &self.seq_names
    }

    /// Returns whether the regions are streamed from a tabix-indexed file rather than held in
    /// memory (HTSlib's `tbx` rather than `regs`).
    #[must_use]
    pub fn is_file_backed(&self) -> bool {
        self.file.is_some()
    }

    /// Returns the current region's sequence index, start, and end (0-based inclusive).
    #[must_use]
    pub fn current(&self) -> Option<(usize, i64, i64)> {
        (self.iseq >= 0).then_some((self.iseq as usize, self.start, self.end))
    }

    /// Returns the current region's sequence name.
    #[must_use]
    pub fn current_seq(&self) -> Option<&[u8]> {
        let i = usize::try_from(self.iseq).ok()?;
        self.seq_names.get(i).map(Vec::as_slice)
    }

    /// Returns the line of the current region, for regions read from a tabix-indexed file
    /// (HTSlib's `line`: the last line read, so after a seek still the line before it).
    #[must_use]
    pub fn current_line(&self) -> Option<&[u8]> {
        let file = self.file.as_ref()?;
        (self.iseq >= 0).then_some(file.line.as_slice())
    }

    /// Returns the alleles of the current target, parsed by allele-aware target matching
    /// (HTSlib's `als`); empty until then.
    #[must_use]
    pub fn alleles(&self) -> &[Vec<u8>] {
        if self.als_parsed { &self.als } else { &[] }
    }

    /// Starts recording the regions passed over without a record (HTSlib's
    /// `missed_reg_handler`), which [`Self::take_missed`] returns.
    pub fn track_missed(&mut self) {
        self.missed.get_or_insert_with(Vec::new);
    }

    /// Returns and clears the regions passed over since the last call.
    pub fn take_missed(&mut self) -> Vec<MissedRegion> {
        self.missed.as_mut().map(std::mem::take).unwrap_or_default()
    }

    /// Restarts iteration (HTSlib's `bcf_sr_seek_start`).
    pub fn rewind(&mut self) {
        for s in &mut self.seqs {
            s.creg = -1;
        }
        if let Some(f) = &mut self.file {
            f.cursor = Cursor::Start;
            f.has_line = false;
        }
        self.iseq = 0;
        self.start = -1;
        self.end = -1;
        self.prev_seq = -1;
        self.prev_start = -1;
        self.prev_end = -1;
    }

    /// Positions at a sequence's first region (HTSlib's `bcf_sr_regions_seek`); `false` when
    /// the sequence has no regions.
    pub fn seek(&mut self, seq: &[u8]) -> bool {
        self.start = -1;
        self.end = -1;
        match self.seq_hash.get(seq) {
            Some(&i) => {
                self.iseq = i as isize;
                if let Some(file) = &mut self.file {
                    // Sequence `i` is the index's sequence `i`. A whole-sequence query always
                    // has a plan; an empty one would end the sequence at once.
                    let plan = file
                        .tbx
                        .index
                        .query(i, 0, WHOLE_SEQUENCE)
                        .unwrap_or_default();
                    file.cursor = Cursor::Query(Query::new(plan));
                    file.has_line = false;
                } else {
                    self.seqs[i].creg = -1;
                }
                true
            }
            None => {
                self.iseq = -1;
                false
            }
        }
    }

    /// Positions the iteration as if the regions were read up to `seq:pos` (as HTSlib's
    /// `bcf_sr_seek` does with `_bcf_sr_regions_overlap`).
    ///
    /// # Errors
    ///
    /// Returns an error for an unreadable or malformed line of a tabix-indexed file.
    pub(super) fn seek_to(&mut self, seq: &[u8], pos: i64) -> Result<(), RegionsError> {
        if let Some(&i) = self.seq_hash.get(seq) {
            self.iseq = i as isize;
        }
        self.overlap_impl(seq, pos, pos, MissedTracking::Off, LineRead::AsNeeded)?;
        Ok(())
    }

    /// Moves to the next region (HTSlib's `bcf_sr_regions_next`); `false` when done.
    ///
    /// # Errors
    ///
    /// Returns an error, reported as HTSlib reports it, for an unreadable or malformed line of
    /// a tabix-indexed file. As in HTSlib, the regions can still be read after it: the next
    /// call reads the next line.
    pub fn next_region(&mut self) -> Result<bool, RegionsError> {
        if self.iseq < 0 {
            return Ok(false);
        }
        self.start = -1;
        self.end = -1;
        self.als_parsed = false;
        if self.file.is_some() {
            return self.next_file_region();
        }
        while (self.iseq as usize) < self.seqs.len() {
            let s = &mut self.seqs[self.iseq as usize];
            let mut i = s.creg + 1;
            while (i as usize) < s.regs.len() && s.regs[i as usize].0 > s.regs[i as usize].1 {
                i += 1;
            }
            s.creg = i;
            if (i as usize) < s.regs.len() {
                break;
            }
            self.iseq += 1;
        }
        if self.iseq as usize >= self.seqs.len() {
            self.iseq = -1;
            return Ok(false);
        }
        let s = &self.seqs[self.iseq as usize];
        let r = s.regs[s.creg as usize];
        self.start = r.0;
        self.end = r.1;
        Ok(true)
    }

    /// [`Self::next_region`] for a tabix-indexed file: reads lines until one with a region.
    fn next_file_region(&mut self) -> Result<bool, RegionsError> {
        let file = self.file.as_mut().expect("called for a tabix-indexed file");
        file.has_line = false;
        loop {
            let got = match file.read_line() {
                Ok(got) => got,
                Err(e) => {
                    // A failed query ends the sequence; a failed sequential read leaves the
                    // state as it was, as in HTSlib.
                    if matches!(file.cursor, Cursor::Query(_)) {
                        self.iseq = -1;
                    } else {
                        log::error(
                            "bcf_sr_regions_next",
                            format_args!("Error reading \"{}\": {e}", file.name),
                        );
                    }
                    return Err(RegionsError {
                        msg: format!("{}: {e}", file.name),
                        read_failed: true,
                    });
                }
            };
            if !got {
                self.iseq = -1;
                return Ok(false);
            }
            let Some(parsed) = parse_region_line(&file.line, file.ichr, file.ifrom, file.ito)
            else {
                log::error(
                    "bcf_sr_regions_next",
                    format_args!(
                        "Could not parse the file {}, using the columns {},{},{}",
                        file.name,
                        file.ichr + 1,
                        file.ifrom + 1,
                        file.ito + 1
                    ),
                );
                return Err(RegionsError::new(format!("could not parse {}", file.name)));
            };
            let Some((chr, from, to)) = parsed else {
                continue;
            };
            let Some(&iseq) = self.seq_hash.get(chr) else {
                log::error(
                    "bcf_sr_regions_next",
                    format_args!(
                        "Broken tabix index? The sequence \"{}\" not in dictionary [{}]",
                        String::from_utf8_lossy(chr),
                        String::from_utf8_lossy(&file.line)
                    ),
                );
                return Err(RegionsError::new(format!(
                    "broken tabix index of {}",
                    file.name
                )));
            };
            let from = if file.bed { from + 1 } else { from };
            self.iseq = iseq as isize;
            self.start = from - 1;
            self.end = to - 1;
            file.has_line = true;
            return Ok(true);
        }
    }

    /// Tests whether `[start, end]` (0-based inclusive) on `seq` overlaps a region, advancing
    /// through the regions as records arrive in order (HTSlib's `bcf_sr_regions_overlap`).
    /// With missed-region tracking on, regions passed over are recorded.
    ///
    /// # Errors
    ///
    /// Returns an error, reported as HTSlib reports it, for an unreadable or malformed line of
    /// a tabix-indexed file. HTSlib's callers count the interval as outside the regions then,
    /// and the next call goes on from the next line.
    pub fn overlaps(&mut self, seq: &[u8], start: i64, end: i64) -> Result<bool, RegionsError> {
        let tracking = if self.missed.is_some() {
            MissedTracking::On
        } else {
            MissedTracking::Off
        };
        let overlap = self.overlap_impl(seq, start, end, tracking, LineRead::AsNeeded)?;
        Ok(overlap == RegionOverlap::Overlap)
    }

    /// [`Self::overlaps`] for regions whose current line is used (bcftools annotate's
    /// `annotate_from_tab`): after a seek, a line is read even for an interval starting at -1
    /// (a record at POS 0), where HTSlib reports an overlap with no line read and bcftools
    /// then parses an empty line and crashes.
    ///
    /// # Errors
    ///
    /// As [`Self::overlaps`].
    pub fn overlaps_line(
        &mut self,
        seq: &[u8],
        start: i64,
        end: i64,
    ) -> Result<bool, RegionsError> {
        let overlap = self.overlap_impl(seq, start, end, MissedTracking::Off, LineRead::Always)?;
        Ok(overlap == RegionOverlap::Overlap)
    }

    /// HTSlib's `_bcf_sr_regions_overlap`.
    fn overlap_impl(
        &mut self,
        seq: &[u8],
        start: i64,
        end: i64,
        tracking: MissedTracking,
        line_read: LineRead,
    ) -> Result<RegionOverlap, RegionsError> {
        let Some(&iseq) = self.seq_hash.get(seq) else {
            return Ok(RegionOverlap::None);
        };
        let iseq = iseq as isize;
        if self.prev_seq == -1 || iseq != self.prev_seq || self.prev_start > start {
            if tracking == MissedTracking::On && self.prev_seq != -1 && self.iseq != -1 {
                self.flush_missed()?;
            }
            if !self.seek(seq) {
                return Ok(RegionOverlap::None);
            }
            self.start = -1;
            self.end = -1;
        }
        if self.prev_seq == iseq && self.iseq != iseq {
            return Ok(RegionOverlap::Done);
        }
        self.prev_seq = self.iseq;
        self.prev_start = start;
        while iseq == self.iseq && (self.end < start || self.needs_line(line_read)) {
            if !self.next_region()? {
                return Ok(RegionOverlap::Done);
            }
            if self.iseq != iseq {
                return Ok(RegionOverlap::None);
            }
            if tracking == MissedTracking::On && self.end < start {
                self.record_missed();
            }
        }
        Ok(if self.start <= end {
            RegionOverlap::Overlap
        } else {
            RegionOverlap::None
        })
    }

    /// Whether a line must be read before testing an overlap under `line_read`.
    fn needs_line(&self, line_read: LineRead) -> bool {
        line_read == LineRead::Always && self.file.as_ref().is_some_and(|f| !f.has_line)
    }

    /// Records every region not yet visited as missed (HTSlib's `bcf_sr_regions_flush`).
    ///
    /// # Errors
    ///
    /// Returns an error for an unreadable or malformed line of a tabix-indexed file.
    pub fn flush_missed(&mut self) -> Result<(), RegionsError> {
        if self.missed.is_none() || self.prev_seq == -1 {
            return Ok(());
        }
        while self.next_region()? {
            self.record_missed();
        }
        Ok(())
    }

    fn record_missed(&mut self) {
        let region = MissedRegion {
            seq: self.current_seq().unwrap_or_default().to_vec(),
            start: self.start,
            end: self.end,
            line: self.current_line().map(<[u8]>::to_vec),
        };
        if let Some(m) = &mut self.missed {
            m.push(region);
        }
    }

    /// Tests whether the current target's alleles (column `als_idx`, 0-based, of its line)
    /// have the type of the record's variants (HTSlib's `_regions_match_alleles`): indels if
    /// any listed allele is longer than one base, SNPs otherwise.
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::ApiUsage`] for targets held in memory, which have no lines.
    pub fn match_alleles(&mut self, als_idx: usize, rec: &mut Record) -> Result<bool, SyncedError> {
        if self.file.is_none() {
            log::error(
                "_regions_match_alleles",
                format_args!("Compressed and indexed targets file is required"),
            );
            return Err(SyncedError::ApiUsage);
        }
        if !self.als_parsed {
            let line = self.current_line().unwrap_or_default();
            let col = line
                .split(|&c| c == b'\t')
                .nth(als_idx)
                .unwrap_or_default()
                .to_vec();
            self.als.clear();
            let mut max_len = 0;
            for a in col.split(|&c| c == b',') {
                max_len = max_len.max(a.len());
                self.als.push(a.to_vec());
            }
            self.als_type = if max_len > 1 {
                AllelesType::Indel
            } else {
                AllelesType::Snp
            };
            self.als_parsed = true;
        }
        let is_indel = rec.variant_types().contains(VariantType::INDEL);
        Ok(match self.als_type {
            AllelesType::Indel => is_indel,
            AllelesType::Snp => !is_indel,
        })
    }
}

/// Whether `_bcf_sr_regions_overlap` reports missed regions.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MissedTracking {
    On,
    Off,
}

/// When `_bcf_sr_regions_overlap` reads the next line.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum LineRead {
    /// While the current region ends before the interval (HTSlib).
    AsNeeded,
    /// Also while a tabix-indexed file has no current line (see [`Regions::overlaps_line`]).
    Always,
}

/// The outcome of `_bcf_sr_regions_overlap`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RegionOverlap {
    /// 0: the interval overlaps a region.
    Overlap,
    /// -1: no overlap.
    None,
    /// -2: no regions are left on the sequence.
    Done,
}

/// Reports a region file that cannot be opened or read, as `hts_open` (naming the data file)
/// and `bcf_sr_regions_init` (the whole name) do.
fn open_failed(name: &str, e: &io::Error) -> RegionsError {
    if let Some(code) = e.raw_os_error() {
        log::error(
            "hts_open_format",
            format_args!(
                "Failed to open file \"{}\" : {}",
                split_idx_name(name).0,
                SyncedError::Open(io::Error::from_raw_os_error(code))
            ),
        );
    }
    log::error(
        "bcf_sr_regions_init",
        format_args!("Could not open file: {name}"),
    );
    RegionsError::new(format!("{name}: {e}"))
}

/// Parses a region file line (HTSlib's `_regions_parse_line`): `Some(None)` for a comment,
/// `Some(Some((chr, from, to)))` with the coordinates as written, or `None` when malformed.
fn parse_region_line(
    line: &[u8],
    ichr: usize,
    ifrom: usize,
    ito: usize,
) -> Option<Option<(&[u8], i64, i64)>> {
    if line.first() == Some(&b'#') {
        return Some(None);
    }
    // The columns needed, without collecting them all (this runs for every line of a file).
    let last = ichr.max(ifrom).max(ito);
    let mut picked: [&[u8]; 3] = [&[]; 3];
    let mut ncols = 0;
    for (i, col) in line.split(|&c| c == b'\t').enumerate().take(last + 1) {
        for (j, want) in [ichr, ifrom, ito].into_iter().enumerate() {
            if i == want {
                picked[j] = col;
            }
        }
        ncols = i + 1;
    }
    let col = |i: usize| -> &[u8] {
        if i == ichr {
            picked[0]
        } else if i == ifrom {
            picked[1]
        } else {
            picked[2]
        }
    };
    let (k, l) = if ifrom <= ito {
        (ifrom, ito)
    } else {
        (ito, ifrom)
    };
    let num = |s: &[u8]| -> Option<i64> {
        let (v, end) = parse_decimal(s, 0, ParseDecimalFlags::default());
        if end == 0 || end != s.len() {
            return None;
        }
        Some(v)
    };
    if ncols <= k || ncols <= l {
        return None;
    }
    let (from, to) = if k == l {
        let v = num(col(k))?;
        (v, v)
    } else {
        let a = num(col(k))?;
        let b = num(col(l))?;
        if k == ifrom { (a, b) } else { (b, a) }
    };
    if ncols <= ichr {
        return None;
    }
    Some(Some((picked[0], from, to)))
}

#[cfg(test)]
mod tests {
    use std::io::Write as _;

    use super::*;
    use crate::{
        bgzf,
        variant::index::tabix::{TabixBuilder, TabixConfig},
    };

    /// A file in a test directory removed on drop.
    struct TestFile {
        dir: PathBuf,
        path: PathBuf,
    }

    impl TestFile {
        /// Writes `text` as `name` in a new directory for test `test`.
        fn plain(test: &str, name: &str, text: &str) -> Self {
            let dir =
                std::env::temp_dir().join(format!("pasta-regions-{}-{test}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("create the test directory");
            let path = dir.join(name);
            std::fs::write(&path, text).expect("write the test file");
            Self { dir, path }
        }

        /// Writes `text` BGZF-compressed as `name`, with a tabix index of `conf`.
        fn tabix(test: &str, name: &str, text: &str, conf: TabixConfig) -> Self {
            let file = Self::plain(test, name, "");
            let mut w = bgzf::io::Writer::new(
                std::fs::File::create(&file.path).expect("create the test file"),
            );
            let mut builder = TabixBuilder::new(conf, 0);
            for line in text.lines() {
                writeln!(w, "{line}").expect("write a line");
                builder
                    .push_line(line.as_bytes(), u64::from(w.virtual_position()))
                    .expect("index a line");
            }
            let tbx = builder.finish(u64::from(w.virtual_position()));
            w.finish().expect("finish the BGZF file");
            let mut idx = file.path.clone().into_os_string();
            idx.push(".tbi");
            tbx.index
                .save(Path::new(&idx), IndexFormat::Tbi)
                .expect("save the index");
            file
        }
    }

    impl Drop for TestFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.dir);
        }
    }

    /// CHROM, BEG, END columns, 1-based.
    const GENERIC: TabixConfig = TabixConfig {
        preset: crate::variant::index::tabix::PRESET_GENERIC,
        sc: 1,
        bc: 2,
        ec: 3,
        meta_char: b'#' as i32,
        line_skip: 0,
    };

    /// Tests an overlap with regions that read nothing, so cannot fail.
    fn overlaps(r: &mut Regions, seq: &[u8], start: i64, end: i64) -> bool {
        r.overlaps(seq, start, end)
            .expect("in-memory regions read nothing")
    }

    /// The current region and its line.
    fn current(r: &Regions) -> Option<(usize, i64, i64, String)> {
        let (iseq, start, end) = r.current()?;
        let line = String::from_utf8_lossy(r.current_line().unwrap_or_default()).into_owned();
        Some((iseq, start, end, line))
    }

    /// A sequence name and its regions as inclusive 0-based intervals.
    type SeqList = (Vec<u8>, Vec<(i64, i64)>);

    fn regs(s: &str) -> Vec<SeqList> {
        let r = Regions::from_string(s.as_bytes()).expect("the regions parse");
        r.seq_names
            .iter()
            .zip(&r.seqs)
            .map(|(n, s)| {
                (
                    n.clone(),
                    s.regs.iter().copied().filter(|(a, b)| a <= b).collect(),
                )
            })
            .collect()
    }

    #[test]
    fn region_strings_parse_as_htslib() {
        assert_eq!(
            regs("chr1"),
            [(b"chr1".to_vec(), vec![(0, MAX_CSI_COOR - 1)])]
        );
        assert_eq!(regs("chr1:100"), [(b"chr1".to_vec(), vec![(99, 99)])]);
        assert_eq!(
            regs("chr1:100-"),
            [(b"chr1".to_vec(), vec![(99, MAX_CSI_COOR - 2)])]
        );
        assert_eq!(regs("chr1:1k-2k"), [(b"chr1".to_vec(), vec![(999, 1999)])]);
        assert_eq!(
            regs("{a:b}:5-6,chr2,chr1:10-20,chr1:15-30"),
            [
                (b"a:b".to_vec(), vec![(4, 5)]),
                (b"chr2".to_vec(), vec![(0, MAX_CSI_COOR - 1)]),
                (b"chr1".to_vec(), vec![(9, 29)]),
            ]
        );
        assert!(Regions::from_string(b"chr1:x").is_err());
        assert!(Regions::from_string(b"{chr1").is_err());
    }

    #[test]
    fn overlaps_advance_in_order() {
        let mut r = Regions::from_string(b"1:10-20,1:30-40,2:5-5").expect("the regions parse");
        assert!(!overlaps(&mut r, b"1", 5, 5));
        assert!(overlaps(&mut r, b"1", 12, 12));
        assert!(!overlaps(&mut r, b"1", 25, 25));
        assert!(overlaps(&mut r, b"1", 29, 35));
        assert!(!overlaps(&mut r, b"1", 50, 50));
        assert!(overlaps(&mut r, b"2", 4, 4));
        assert!(!overlaps(&mut r, b"3", 4, 4));
    }

    #[test]
    fn missed_regions_are_recorded() {
        let mut r = Regions::from_string(b"1:10-20,1:30-40,1:50-60").expect("the regions parse");
        r.track_missed();
        assert!(overlaps(&mut r, b"1", 54, 54));
        let missed = r.take_missed();
        assert_eq!(
            missed.iter().map(|m| (m.start, m.end)).collect::<Vec<_>>(),
            [(9, 19), (29, 39)]
        );
        r.flush_missed().expect("in-memory regions read nothing");
        assert!(r.take_missed().is_empty());
    }

    #[test]
    fn lines_parse_as_htslib() {
        assert_eq!(
            parse_region_line(b"chr1\t10\t20\tx", 0, 1, 2),
            Some(Some((b"chr1".as_slice(), 10, 20)))
        );
        assert_eq!(
            parse_region_line(b"chr1\t10\tx", 0, 1, 1),
            Some(Some((b"chr1".as_slice(), 10, 10)))
        );
        assert_eq!(parse_region_line(b"#chr1\t10", 0, 1, 1), Some(None));
        assert_eq!(parse_region_line(b"chr1\t1x", 0, 1, 1), None);
        assert_eq!(parse_region_line(b"chr1", 0, 1, 1), None);
        assert_eq!(parse_region_line(b"", 0, 1, 1), None);
    }

    #[test]
    fn plain_files_are_read_up_to_the_first_empty_line() {
        let f = TestFile::plain("plain", "r.txt", "1\t10\t20\n2\t5\t5\n\n1\t30\t40\n");
        let r = Regions::from_file(&f.path).expect("the regions parse");
        assert!(!r.is_file_backed());
        assert_eq!(r.seq_names(), [b"1".to_vec(), b"2".to_vec()]);
        assert_eq!(r.seqs[0].regs, [(9, 19)]);
        assert_eq!(r.seqs[1].regs, [(4, 4)]);

        // Without an END column on the first line, every region is one base.
        let f = TestFile::plain("one", "r.txt", "1\t10\n1\t30\t40\n");
        let r = Regions::from_file(&f.path).expect("the regions parse");
        assert_eq!(r.seqs[0].regs, [(9, 9), (29, 29)]);

        let f = TestFile::plain("bed", "r.bed", "1\t10\t20\n");
        let r = Regions::from_file(&f.path).expect("the regions parse");
        assert_eq!(r.seqs[0].regs, [(10, 19)]);

        let f = TestFile::plain("vcf", "r.vcf", "##fileformat=VCFv4.2\n1\t10\t20\n");
        let r = Regions::from_file(&f.path).expect("the regions parse");
        assert_eq!(r.seqs[0].regs, [(9, 9)]);

        let f = TestFile::plain("bad", "r.txt", "1\t10\t20\n1\tx\n");
        assert!(Regions::from_file(&f.path).is_err());
        let f = TestFile::plain("empty", "r.txt", "\n1\t10\t20\n");
        assert!(Regions::from_file(&f.path).is_err());
    }

    #[test]
    fn idx_names_read_the_data_file() {
        let path = |f: &TestFile, index: &str| {
            let mut name = f.path.clone().into_os_string();
            name.push(format!("##idx##{}", f.dir.join(index).display()));
            PathBuf::from(name)
        };
        // An index that cannot be loaded leaves the file to be read whole, by its own name.
        let f = TestFile::plain("idxplain", "r.bed", "1\t10\t20\n");
        let r = Regions::from_file(&path(&f, "missing.tbi")).expect("the regions parse");
        assert!(!r.is_file_backed());
        assert_eq!(r.seqs[0].regs, [(10, 19)]);

        // A named index streams the file through it.
        let f = TestFile::tabix("idxtbx", "r.txt.gz", "1\t10\t20\n", GENERIC);
        let tbi = f.dir.join("other.tbi");
        let mut sibling = f.path.clone().into_os_string();
        sibling.push(".tbi");
        std::fs::rename(&sibling, &tbi).expect("move the index");
        assert!(
            !Regions::from_file(&f.path)
                .expect("the file reads")
                .is_file_backed()
        );
        let mut r = Regions::from_file(&path(&f, "other.tbi")).expect("the file opens");
        assert!(r.is_file_backed());
        assert!(r.next_region().expect("the line parses"));
        assert_eq!(current(&r), Some((0, 9, 19, "1\t10\t20".to_owned())));
    }

    #[test]
    fn tabix_files_stream_in_file_order_then_by_sequence() {
        let text = "#c\n1\t10\t20\tA\n1\t30\t40\tB\n2\t5\t5\tC\n";
        let f = TestFile::tabix("stream", "r.txt.gz", text, GENERIC);
        let mut r = Regions::from_file(&f.path).expect("the file opens");
        assert!(r.is_file_backed());
        assert_eq!(r.seq_names(), [b"1".to_vec(), b"2".to_vec()]);
        let mut got = Vec::new();
        while r.next_region().expect("the lines parse") {
            got.push(current(&r).expect("a region is current"));
        }
        assert_eq!(
            got,
            [
                (0, 9, 19, "1\t10\t20\tA".to_owned()),
                (0, 29, 39, "1\t30\t40\tB".to_owned()),
                (1, 4, 4, "2\t5\t5\tC".to_owned()),
            ]
        );

        // Overlap tests query a sequence at a time.
        r.rewind();
        assert!(r.overlaps(b"2", 4, 4).expect("the lines parse"));
        assert_eq!(current(&r), Some((1, 4, 4, "2\t5\t5\tC".to_owned())));
        assert!(r.overlaps(b"1", 35, 35).expect("the lines parse"));
        assert_eq!(current(&r), Some((0, 29, 39, "1\t30\t40\tB".to_owned())));
        assert!(!r.overlaps(b"1", 45, 45).expect("the lines parse"));
        assert!(
            !r.overlaps(b"3", 1, 1)
                .expect("an unknown sequence reads nothing")
        );

        // A rewind goes back to the first line.
        r.rewind();
        assert!(r.next_region().expect("the lines parse"));
        assert_eq!(current(&r), Some((0, 9, 19, "1\t10\t20\tA".to_owned())));
    }

    #[test]
    fn a_malformed_tabix_line_fails_when_reached() {
        let text = "1\t10\t10\n1\t25x\t35\n1\t40\t40\n";
        let f = TestFile::tabix("malformed", "r.txt.gz", text, GENERIC);
        let mut r = Regions::from_file(&f.path).expect("the file opens");
        assert!(r.overlaps(b"1", 9, 9).expect("the first line parses"));
        let e = r
            .overlaps(b"1", 19, 19)
            .expect_err("the second line is malformed");
        assert!(!e.is_read_failure(), "{e}");
        // As in HTSlib, the next test goes on from the next line.
        assert!(!r.overlaps(b"1", 29, 29).expect("the next line parses"));
        assert!(r.overlaps(b"1", 39, 39).expect("the next line parses"));

        let mut r = Regions::from_file(&f.path).expect("the file opens");
        assert!(r.next_region().expect("the first line parses"));
        assert!(r.next_region().is_err());
        assert!(r.next_region().expect("the next line parses"));
        assert_eq!(current(&r), Some((0, 39, 39, "1\t40\t40".to_owned())));
    }

    #[test]
    fn a_corrupt_block_of_a_tabix_file_is_a_read_failure() {
        let text: String = (1..20_000)
            .map(|i| format!("1\t{}\t{}\t{i:040}\n", i * 10, i * 10 + 5))
            .collect();
        let f = TestFile::tabix("corrupt", "r.txt.gz", &text, GENERIC);
        // A wrong CRC32 in the middle block.
        let mut data = std::fs::read(&f.path).expect("read the file");
        let mut blocks = Vec::new();
        let mut p = 0;
        while p < data.len() {
            let size = usize::from(u16::from_le_bytes([data[p + 16], data[p + 17]])) + 1;
            blocks.push((p, size));
            p += size;
        }
        assert!(blocks.len() > 3, "{} blocks", blocks.len());
        let (start, size) = blocks[blocks.len() / 2];
        data[start + size - 8] ^= 0xff;
        std::fs::write(&f.path, data).expect("write the file");

        for in_order in [true, false] {
            let mut r = Regions::from_file(&f.path).expect("the file opens");
            let mut tested = 0;
            let e = loop {
                tested += 1;
                let result = if in_order {
                    r.next_region()
                } else {
                    r.overlaps(b"1", tested * 10 - 1, tested * 10 - 1)
                };
                match result {
                    Ok(found) => assert!(found, "line {tested} is there"),
                    Err(e) => break e,
                }
            };
            assert!(e.is_read_failure(), "{e}");
            assert!(
                tested > 1000,
                "the lines before the block are read: {tested}"
            );
        }
    }

    #[test]
    fn overlaps_line_reads_a_line_for_pos_zero() {
        let text = "1\t10\t20\n";
        let f = TestFile::tabix("pos0", "r.txt.gz", text, GENERIC);
        let mut r = Regions::from_file(&f.path).expect("the file opens");
        // HTSlib counts POS 0 (start -1) as overlapping before reading any line.
        assert!(r.overlaps(b"1", -1, -1).expect("nothing is read"));
        let mut r = Regions::from_file(&f.path).expect("the file opens");
        assert!(!r.overlaps_line(b"1", -1, -1).expect("the line parses"));
        assert_eq!(current(&r), Some((0, 9, 19, "1\t10\t20".to_owned())));
    }

    #[test]
    fn tabix_bed_files_have_zero_based_starts() {
        let f = TestFile::tabix("tbxbed", "r.bed.gz", "1\t10\t20\n", TabixConfig::BED);
        let mut r = Regions::from_file(&f.path).expect("the file opens");
        assert!(r.next_region().expect("the line parses"));
        assert_eq!(r.current(), Some((0, 10, 19)));
    }
}
