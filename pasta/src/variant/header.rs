//! The VCF/BCF header (HTSlib's `bcf_hdr_t`): header lines in order, the shared
//! FILTER/INFO/FORMAT ID dictionary, the contig dictionary, and the sample list.
//!
//! Dictionary IDs follow HTSlib exactly: FILTER, INFO, and FORMAT share one index space, `PASS`
//! is always ID 0, `IDX=` fields in BCF headers fix IDs (leaving gaps), and removing a line
//! keeps its dictionary entry (records referring to it still format), as in HTSlib.

mod ops;
mod record;

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use rustc_hash::FxHashMap;

pub use self::record::{HeaderRecord, HeaderRecordType, HeaderTextFlavor, Quoting};
pub(crate) use self::record::{ParsedLine, parse_line};
use crate::variant::{log, number};

/// The VCF version 4.2, as HTSlib encodes versions (`major * 1000000 + minor * 1000`).
pub const VCF_DEFAULT_VERSION: i32 = 4_002_000;
/// The VCF version 4.4.
pub const VCF44: i32 = 4_004_000;
/// The VCF version 4.5.
pub const VCF45: i32 = 4_005_000;

/// The value type of an INFO or FORMAT field (HTSlib's `BCF_HT_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[repr(u8)]
pub enum ValueType {
    /// `Type=Flag` (`BCF_HT_FLAG`).
    Flag = 0,
    /// `Type=Integer` (`BCF_HT_INT`).
    Integer = 1,
    /// `Type=Float` (`BCF_HT_REAL`).
    Float = 2,
    /// `Type=String` or `Type=Character` (`BCF_HT_STR`).
    String = 3,
}

impl ValueType {
    /// Converts a raw `BCF_HT_*` code; 15 (undefined) and others give `None`.
    #[inline]
    #[must_use]
    pub fn from_code(code: u32) -> Option<Self> {
        Some(match code {
            0 => Self::Flag,
            1 => Self::Integer,
            2 => Self::Float,
            3 => Self::String,
            _ => return None,
        })
    }
}

/// The `Number=` of an INFO or FORMAT definition (HTSlib's `BCF_VL_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Number {
    /// A fixed count (`BCF_VL_FIXED`).
    Fixed(u32),
    /// `.`, any number (`BCF_VL_VAR`).
    Variable,
    /// `A`, one per ALT allele.
    A,
    /// `G`, one per genotype.
    G,
    /// `R`, one per allele including REF.
    R,
    /// `P`, one per allele in GT (FORMAT only, VCF 4.4).
    P,
    /// `LA` (FORMAT only, VCF 4.5).
    LocalA,
    /// `LG` (FORMAT only, VCF 4.5).
    LocalG,
    /// `LR` (FORMAT only, VCF 4.5).
    LocalR,
    /// `M`, one per base modification (FORMAT only, VCF 4.5).
    M,
}

impl Number {
    /// Returns the `BCF_VL_*` code.
    #[must_use]
    pub fn code(self) -> u32 {
        match self {
            Self::Fixed(_) => 0,
            Self::Variable => 1,
            Self::A => 2,
            Self::G => 3,
            Self::R => 4,
            Self::P => 5,
            Self::LocalA => 6,
            Self::LocalG => 7,
            Self::LocalR => 8,
            Self::M => 9,
        }
    }

    fn from_info(info: u32) -> Option<Self> {
        Some(match (info >> 8) & 0xf {
            0 => Self::Fixed(info >> 12),
            1 => Self::Variable,
            2 => Self::A,
            3 => Self::G,
            4 => Self::R,
            5 => Self::P,
            6 => Self::LocalA,
            7 => Self::LocalG,
            8 => Self::LocalR,
            9 => Self::M,
            _ => return None,
        })
    }
}

/// An error building or modifying a header.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum HeaderError {
    /// A header line could not be parsed.
    InvalidLine(String),
    /// Two lines claim the same `IDX`.
    ConflictingIdx(i32, String),
    /// The `#CHROM` line is malformed.
    InvalidSampleLine,
    /// The header has no `#CHROM` line.
    MissingSampleLine,
    /// A sample name is empty.
    EmptySampleName,
    /// A sample name appears twice.
    DuplicateSample(String),
    /// A sample list names a sample that is not in the header (1-based position in the list).
    UnknownSample(usize),
    /// A sample list file could not be read.
    SampleList(String),
}

impl std::fmt::Display for HeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidLine(s) => write!(f, "could not parse the header line: {s}"),
            Self::ConflictingIdx(i, t) => {
                write!(
                    f,
                    "conflicting IDX={i} lines in the header dictionary at {t}"
                )
            }
            Self::InvalidSampleLine => f.write_str("could not parse the #CHROM line"),
            Self::MissingSampleLine => f.write_str("sample line not found"),
            Self::EmptySampleName => f.write_str("empty sample name"),
            Self::DuplicateSample(s) => write!(f, "duplicated sample name '{s}'"),
            Self::UnknownSample(i) => write!(f, "sample {i} of the list is not in the header"),
            Self::SampleList(s) => write!(f, "could not read the sample list {s}"),
        }
    }
}

impl std::error::Error for HeaderError {}

/// An entry of the shared FILTER/INFO/FORMAT dictionary (HTSlib's `bcf_idinfo_t`).
#[derive(Clone, Debug)]
pub(crate) struct IdEntry {
    pub(crate) name: Box<[u8]>,
    /// Packed `Number:20 | var:4 | Type:4 | ColType:4` per FILTER/INFO/FORMAT, 15 if undefined.
    pub(crate) info: [u32; 3],
    /// The defining line of each type, by serial.
    pub(crate) hrec: [Option<u32>; 3],
}

/// An entry of the contig dictionary.
#[derive(Clone, Debug)]
pub(crate) struct ContigEntry {
    pub(crate) name: Box<[u8]>,
    pub(crate) length: i64,
    pub(crate) hrec: Option<u32>,
}

/// Sample selection applied while reading (HTSlib's `bcf_hdr_set_samples`).
#[derive(Clone, Debug)]
pub(crate) struct KeepSamples {
    pub(crate) bits: Vec<u8>,
    pub(crate) n_original: usize,
}

impl KeepSamples {
    #[inline]
    pub(crate) fn contains(&self, i: usize) -> bool {
        self.bits[i / 8] & (1 << (i % 8)) != 0
    }
}

/// The translation tables built by [`crate::variant::Record::translate`] (`bcf_translate`),
/// for the destination whose dictionaries had the [`Dictionaries`] stamp `dst`.
#[derive(Clone, Debug, Default)]
pub(crate) enum Translation {
    /// Not computed yet.
    #[default]
    Unknown,
    /// No ID differs between the headers.
    Identity { dst: u64 },
    /// Per source ID (ID dictionary, contig dictionary): the destination ID or -1.
    Map {
        dst: u64,
        ids: Vec<i32>,
        ctgs: Vec<i32>,
    },
}

impl Translation {
    /// Whether the tables were built for the destination with the stamp `dst`.
    pub(crate) fn is_for(&self, dst: u64) -> bool {
        matches!(self, Self::Identity { dst: d } | Self::Map { dst: d, .. } if *d == dst)
    }
}

/// A stamp of a header's ID and contig dictionaries, new for every header, every clone, and
/// every change to them, so that translation tables are never used for another destination.
/// HTSlib keeps the tables of the first destination for good.
#[derive(Debug)]
pub(crate) struct Dictionaries(u64);

impl Dictionaries {
    fn next() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        Self(NEXT.fetch_add(1, Ordering::Relaxed))
    }

    pub(crate) fn stamp(&self) -> u64 {
        self.0
    }
}

impl Default for Dictionaries {
    fn default() -> Self {
        Self::next()
    }
}

impl Clone for Dictionaries {
    fn clone(&self) -> Self {
        Self::next()
    }
}

/// A VCF/BCF header (HTSlib's `bcf_hdr_t`).
#[derive(Clone, Debug, Default)]
pub struct Header {
    /// Header lines by serial; `None` once removed.
    slab: Vec<Option<HeaderRecord>>,
    /// Serials of the lines in header order.
    order: Vec<u32>,
    pub(crate) ids: Vec<Option<IdEntry>>,
    id_map: FxHashMap<Box<[u8]>, u32>,
    pub(crate) contigs: Vec<Option<ContigEntry>>,
    contig_map: FxHashMap<Box<[u8]>, u32>,
    pub(crate) samples: Vec<Box<[u8]>>,
    sample_map: FxHashMap<Box<[u8]>, u32>,
    /// `##key=value` and `##key=<ID=id>` lines, for duplicate detection and lookups.
    generic: FxHashMap<Vec<u8>, u32>,
    /// The cached version (0 when unknown).
    version: i32,
    pub(crate) keep_samples: Option<KeepSamples>,
    pub(crate) translation: Translation,
    pub(crate) dictionaries: Dictionaries,
    /// The IDs of tags HTSlib looks up for every record: END, SVLEN, LEN, and GT.
    special: [Option<u32>; 4],
}

/// The tags cached in [`Header::special_id`], in order.
const SPECIAL_TAGS: [&[u8]; 4] = [b"END", b"SVLEN", b"LEN", b"GT"];

/// A tag whose ID the header caches for per-record lookups.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SpecialTag {
    End = 0,
    Svlen = 1,
    Len = 2,
    Gt = 3,
}

const PASS_LINE: &[u8] = b"##FILTER=<ID=PASS,Description=\"All filters passed\">";

impl Header {
    /// Creates an empty header (HTSlib's `bcf_hdr_init("r")`).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a header for writing, with `##fileformat=VCFv4.2` and the `PASS` filter
    /// (HTSlib's `bcf_hdr_init("w")`).
    #[must_use]
    pub fn new_for_writing() -> Self {
        let mut h = Self::new();
        h.append(b"##fileformat=VCFv4.2")
            .expect("the default fileformat line is valid");
        h.append(PASS_LINE).expect("the PASS line is valid");
        h.version = VCF_DEFAULT_VERSION;
        h
    }

    /// Parses header text (HTSlib's `bcf_hdr_parse`), ending with the `#CHROM` line.
    ///
    /// # Errors
    ///
    /// Returns an error when the sample line is missing or malformed, or IDs conflict.
    pub fn parse(text: &[u8]) -> Result<Self, HeaderError> {
        let mut h = Self::new();
        h.parse_into(text)?;
        Ok(h)
    }

    /// Parses header text into this (empty) header.
    pub(crate) fn parse_into(&mut self, text: &[u8]) -> Result<(), HeaderError> {
        // Check sanity: "fileformat" string must come as first. HTSlib parses the first line
        // here and again in the loop below (where it is a duplicate).
        let first = match parse_line(text) {
            ParsedLine::Record(r, _) => Some(r),
            _ => None,
        };
        if first
            .as_ref()
            .is_none_or(|r| !r.key.eq_ignore_ascii_case(b"fileformat"))
        {
            log::warn(
                "bcf_hdr_parse",
                format_args!(
                    "The first line should be ##fileformat; is the VCF/BCF header broken?"
                ),
            );
        }
        if let Some(r) = first {
            self.add_record(r)?;
        }
        // The filter PASS must appear first in the dictionary
        self.append(PASS_LINE)?;

        let mut p = 0;
        loop {
            loop {
                match parse_line(&text[p..]) {
                    ParsedLine::Record(r, len) => {
                        self.add_record(r)?;
                        p += len;
                    }
                    ParsedLine::Malformed(len) => p += len,
                    ParsedLine::NotMeta => break,
                }
            }
            let rest = &text[p..];
            if rest.starts_with(b"#CHROM\t") || rest.starts_with(b"#CHROM ") {
                break;
            }
            let eol = rest.iter().position(|&c| c == b'\n');
            if number::at(rest, 0) != 0 {
                let line = &rest[..eol.unwrap_or(rest.len())];
                let line = &line[..line.iter().position(|&c| c == 0).unwrap_or(line.len())];
                log::warn(
                    "bcf_hdr_parse",
                    format_args!(
                        "Could not parse header line: {}",
                        log::strprint(line, b'"', 320)
                    ),
                );
            }
            match eol {
                Some(e) => p += e + 1,
                None => {
                    log::error(
                        "bcf_hdr_parse",
                        format_args!("Could not parse the header, sample line not found"),
                    );
                    return Err(HeaderError::MissingSampleLine);
                }
            }
        }
        self.parse_sample_line(&text[p..])?;
        self.check_sanity();
        Ok(())
    }

    /// Parses the `#CHROM` line and adds its samples (HTSlib's `bcf_hdr_parse_sample_line`).
    fn parse_sample_line(&mut self, s: &[u8]) -> Result<(), HeaderError> {
        const MANDATORY: &[u8] = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let line_for_msg = || {
            let end = s.iter().position(|&c| c == 0).unwrap_or(s.len());
            String::from_utf8_lossy(&s[..end]).into_owned()
        };
        if !s.starts_with(MANDATORY) {
            log::error(
                "bcf_hdr_parse_sample_line",
                format_args!(
                    "Could not parse the \"#CHROM..\" line, either the fields are incorrect or spaces are present instead of tabs:\n\t{}",
                    line_for_msg()
                ),
            );
            return Err(HeaderError::InvalidSampleLine);
        }
        let mut beg = MANDATORY.len();
        let c = number::at(s, beg);
        if c == 0 || c == b'\n' {
            return Ok(());
        }
        if !s[beg..].starts_with(b"\tFORMAT\t") {
            log::error(
                "bcf_hdr_parse_sample_line",
                format_args!(
                    "Could not parse the \"#CHROM..\" line, either FORMAT is missing or spaces are present instead of tabs:\n\t{}",
                    line_for_msg()
                ),
            );
            return Err(HeaderError::InvalidSampleLine);
        }
        beg += 8;
        while number::at(s, beg) != 0 {
            let mut end = beg;
            while !matches!(number::at(s, end), 0 | b'\t' | b'\n') {
                end += 1;
            }
            self.add_sample(&s[beg..end])?;
            if matches!(number::at(s, end), 0 | b'\n') {
                break;
            }
            beg = end + 1;
        }
        Ok(())
    }

    /// Adds a sample (HTSlib's `bcf_hdr_add_sample`).
    ///
    /// # Errors
    ///
    /// Returns an error for an empty (or all-whitespace) or duplicate name.
    pub fn add_sample(&mut self, name: &[u8]) -> Result<(), HeaderError> {
        if name.iter().all(|&c| number::is_space(c)) || name.first() == Some(&0) {
            log::error(
                "bcf_hdr_add_sample_len",
                format_args!("Empty sample name: trailing spaces/tabs in the header line?"),
            );
            return Err(HeaderError::EmptySampleName);
        }
        if self.sample_map.contains_key(name) {
            let s = String::from_utf8_lossy(name).into_owned();
            log::error(
                "bcf_hdr_add_sample_len",
                format_args!("Duplicated sample name '{s}'"),
            );
            return Err(HeaderError::DuplicateSample(s));
        }
        let id = self.samples.len() as u32;
        self.sample_map.insert(name.into(), id);
        self.samples.push(name.into());
        Ok(())
    }

    /// Parses a `##` line and adds it (HTSlib's `bcf_hdr_append`).
    ///
    /// # Errors
    ///
    /// Returns an error when the line cannot be parsed or its IDX conflicts.
    pub fn append(&mut self, line: &[u8]) -> Result<bool, HeaderError> {
        match parse_line(line) {
            ParsedLine::Record(r, _) => self.add_record(r),
            _ => Err(HeaderError::InvalidLine(
                String::from_utf8_lossy(line).into_owned(),
            )),
        }
    }

    /// Adds a header line (HTSlib's `bcf_hdr_add_hrec`), registering FILTER/INFO/FORMAT/contig
    /// definitions in the dictionaries. Returns whether the line was added: duplicates of
    /// existing lines are dropped.
    ///
    /// # Errors
    ///
    /// Returns an error when the line's `IDX` conflicts with another line's.
    pub fn add_record(&mut self, mut hrec: HeaderRecord) -> Result<bool, HeaderError> {
        hrec.set_type();
        check_record(&hrec);
        let serial = self.slab.len() as u32;
        let registered = self.register(&mut hrec, serial)?;
        let mut gen_key: Option<Vec<u8>> = None;
        if !registered {
            if hrec.kind != HeaderRecordType::Generic {
                return Ok(false);
            }
            let key = gen_key_generic(&hrec);
            if self.generic.contains_key(&key) {
                return Ok(false);
            }
            if hrec.key == b"fileformat" {
                self.version = parse_version_string(hrec.value.as_deref().unwrap_or(b""));
            }
            gen_key = Some(key);
        }
        if hrec.kind == HeaderRecordType::Structured
            && let Some(i) = hrec.find_key(b"ID")
        {
            let key = gen_key_structured(&hrec.key, &hrec.vals[i]);
            if self.generic.contains_key(&key) {
                return Ok(false);
            }
            gen_key = Some(key);
        }
        if let Some(key) = gen_key {
            self.generic.insert(key, serial);
        }
        self.slab.push(Some(hrec));
        self.order.push(serial);
        Ok(true)
    }

    /// Registers a line in the dictionaries (HTSlib's `bcf_hdr_register_hrec`); returns
    /// whether it is to be added (a new structured line or definition).
    fn register(&mut self, hrec: &mut HeaderRecord, serial: u32) -> Result<bool, HeaderError> {
        self.dictionaries_changed();
        match hrec.kind {
            HeaderRecordType::Contig => self.register_contig(hrec, serial),
            HeaderRecordType::Structured => Ok(true),
            HeaderRecordType::Generic => Ok(false),
            _ => self.register_id(hrec, serial),
        }
    }

    fn register_contig(
        &mut self,
        hrec: &mut HeaderRecord,
        serial: u32,
    ) -> Result<bool, HeaderError> {
        let len = match hrec.find_key(b"length") {
            None => 0,
            Some(i) => {
                let (len, end) = number::strtol(&hrec.vals[i], 0);
                if end == 0 || len < 0 {
                    return Ok(false);
                }
                len
            }
        };
        let Some(i) = hrec.find_key(b"ID") else {
            return Ok(false);
        };
        let name: Box<[u8]> = hrec.vals[i].as_slice().into();
        let existing = self.contig_map.get(&name).copied();
        if let Some(id) = existing
            && self.contigs[id as usize]
                .as_ref()
                .is_some_and(|e| e.hrec.is_some())
        {
            return Ok(false);
        }
        let mut idx: i64 = -1;
        if let Some(k) = hrec.find_key(b"IDX") {
            match parse_idx(&hrec.vals[k]) {
                Some(v) => idx = v,
                None => {
                    log::warn(
                        "bcf_hdr_register_hrec",
                        format_args!("Error parsing the IDX tag, skipping"),
                    );
                    return Ok(false);
                }
            }
        }
        let id = self.set_idx(Dict::Contig, &name, idx)?;
        self.contig_map.insert(name.clone(), id);
        self.contigs[id as usize] = Some(ContigEntry {
            name,
            length: len,
            hrec: Some(serial),
        });
        if idx == -1 {
            add_idx(hrec, id);
        }
        Ok(true)
    }

    fn register_id(&mut self, hrec: &mut HeaderRecord, serial: u32) -> Result<bool, HeaderError> {
        let coltype = hrec.kind as u32;
        let mut id: Option<Vec<u8>> = None;
        let mut ty = u32::MAX;
        let mut var = u32::MAX;
        let mut num: i32 = -1;
        let mut idx: i64 = -1;
        let is_fmt = hrec.kind == HeaderRecordType::Format;
        for i in 0..hrec.keys.len() {
            let (k, v) = (&hrec.keys[i], &hrec.vals[i]);
            match k.as_slice() {
                b"ID" => id = Some(v.clone()),
                b"IDX" => match parse_idx(v) {
                    Some(x) => idx = x,
                    None => {
                        log::warn(
                            "bcf_hdr_register_hrec",
                            format_args!("Error parsing the IDX tag, skipping"),
                        );
                        return Ok(false);
                    }
                },
                b"Type" => {
                    ty = match v.as_slice() {
                        b"Integer" => 1,
                        b"Float" => 2,
                        b"String" | b"Character" => 3,
                        b"Flag" => 0,
                        _ => {
                            log::warn(
                                "bcf_hdr_register_hrec",
                                format_args!(
                                    "The type \"{}\" is not supported, assuming \"String\"",
                                    String::from_utf8_lossy(v)
                                ),
                            );
                            3
                        }
                    }
                }
                b"Number" => {
                    var = match v.as_slice() {
                        b"A" => 2,
                        b"R" => 4,
                        b"G" => 3,
                        b"." => 1,
                        b"P" if is_fmt => 5,
                        b"LA" if is_fmt => 6,
                        b"LR" if is_fmt => 8,
                        b"LG" if is_fmt => 7,
                        b"M" if is_fmt => 9,
                        _ => match scanf_int(v) {
                            Some(n) => {
                                num = n;
                                0
                            }
                            None => var,
                        },
                    };
                    if var != 0 {
                        num = 0xfffff;
                    }
                }
                _ => {}
            }
        }
        if matches!(hrec.kind, HeaderRecordType::Info | HeaderRecordType::Format) {
            let article = if hrec.key.first() == Some(&b'I') {
                "An"
            } else {
                "A"
            };
            let key = String::from_utf8_lossy(&hrec.key).into_owned();
            if ty == u32::MAX {
                log::warn(
                    "bcf_hdr_register_hrec",
                    format_args!("{article} {key} field has no Type defined. Assuming String"),
                );
                ty = 3;
            }
            if var == u32::MAX {
                log::warn(
                    "bcf_hdr_register_hrec",
                    format_args!("{article} {key} field has no Number defined. Assuming '.'"),
                );
                var = 1;
            }
            if ty == 0 && (var != 0 || num != 0) {
                log::warn(
                    "bcf_hdr_register_hrec",
                    format_args!(
                        "The definition of Flag \"{key}/{}\" is invalid, forcing Number=0",
                        String::from_utf8_lossy(id.as_deref().unwrap_or(b"(null)"))
                    ),
                );
                var = 0;
                num = 0;
            }
        }
        let info = ((num as u32) & 0xfffff) << 12 | (var & 0xf) << 8 | (ty & 0xf) << 4 | coltype;
        let Some(name) = id else {
            return Ok(false);
        };
        let slot = coltype as usize;
        if let Some(&existing) = self.id_map.get(name.as_slice()) {
            let entry = self.ids[existing as usize]
                .as_mut()
                .expect("dictionary entries exist for mapped names");
            if entry.hrec[slot].is_some() {
                return Ok(false);
            }
            entry.info[slot] = info;
            entry.hrec[slot] = Some(serial);
            if idx == -1 {
                add_idx(hrec, existing);
            }
            return Ok(true);
        }
        let name: Box<[u8]> = name.into();
        let id = self.set_idx(Dict::Id, &name, idx)?;
        let mut entry = IdEntry {
            name: name.clone(),
            info: [15, 15, 15],
            hrec: [None, None, None],
        };
        entry.info[slot] = info;
        entry.hrec[slot] = Some(serial);
        if let Some(k) = SPECIAL_TAGS.iter().position(|t| **t == *name) {
            self.special[k] = Some(id);
        }
        self.ids[id as usize] = Some(entry);
        self.id_map.insert(name, id);
        if idx == -1 {
            add_idx(hrec, id);
        }
        Ok(true)
    }

    /// Picks the dictionary slot for a new entry (HTSlib's `bcf_hdr_set_idx`).
    fn set_idx(&mut self, dict: Dict, tag: &[u8], idx: i64) -> Result<u32, HeaderError> {
        let n = match dict {
            Dict::Id => self.ids.len(),
            Dict::Contig => self.contigs.len(),
        };
        let occupied = |h: &Self, i: usize| match dict {
            Dict::Id => h.ids.get(i).is_some_and(Option::is_some),
            Dict::Contig => h.contigs.get(i).is_some_and(Option::is_some),
        };
        let id = if idx == -1 {
            n
        } else {
            let i = idx as usize;
            if i < n && occupied(self, i) {
                let t = String::from_utf8_lossy(tag).into_owned();
                log::error(
                    "bcf_hdr_set_idx",
                    format_args!(
                        "Conflicting IDX={idx} lines in the header dictionary, the new tag is {t}"
                    ),
                );
                return Err(HeaderError::ConflictingIdx(idx as i32, t));
            }
            i
        };
        let new_n = (id + 1).max(n);
        match dict {
            Dict::Id => self.ids.resize(new_n, None),
            Dict::Contig => self.contigs.resize(new_n, None),
        }
        Ok(id as u32)
    }

    /// Removes lines of a kind (HTSlib's `bcf_hdr_remove`): all of them with `key` `None`,
    /// otherwise those with that `ID` (or, for generic lines, that key).
    pub fn remove(&mut self, kind: HeaderRecordType, key: Option<&[u8]>) {
        self.dictionaries_changed();
        let Some(key) = key else {
            let mut i = 0;
            while i < self.order.len() {
                let serial = self.order[i];
                let rec_kind = self.record_by_serial(serial).map(HeaderRecord::kind);
                if rec_kind != Some(kind) {
                    i += 1;
                    continue;
                }
                self.unregister(serial);
                self.remove_from_gen(serial);
                self.order.remove(i);
                self.slab[serial as usize] = None;
            }
            return;
        };
        loop {
            let serial = match kind {
                HeaderRecordType::Filter
                | HeaderRecordType::Info
                | HeaderRecordType::Format
                | HeaderRecordType::Contig => {
                    let Some(serial) = self.dict_hrec(kind, key) else {
                        return;
                    };
                    self.unregister(serial);
                    serial
                }
                _ => {
                    let found = self.order.iter().copied().find(|&s| {
                        let r = self.record_by_serial(s).expect("ordered lines exist");
                        if r.kind != kind {
                            return false;
                        }
                        if kind == HeaderRecordType::Generic {
                            r.key == key
                        } else {
                            r.find_key(b"ID").is_some_and(|j| r.vals[j] == key)
                        }
                    });
                    let Some(serial) = found else { return };
                    self.remove_from_gen(serial);
                    serial
                }
            };
            if let Some(pos) = self.order.iter().position(|&s| s == serial) {
                self.order.remove(pos);
            }
            self.slab[serial as usize] = None;
        }
    }

    /// Marks a FILTER/INFO/FORMAT ID as undefined for `kind` after its line was removed, so
    /// [`Self::id_exists`] is false for it until it is defined again (what bcftools annotate
    /// does with `info[type] |= 0xf` when it removes a class of lines).
    pub fn undefine_id(&mut self, kind: HeaderRecordType, name: &[u8]) {
        let Some(slot) = kind.id_slot() else { return };
        if let Some(&id) = self.id_map.get(name)
            && let Some(e) = self.ids[id as usize].as_mut()
        {
            e.info[slot] |= 0xf;
        }
    }

    /// Renames a FILTER/INFO/FORMAT ID in place, keeping its numeric ID: the dictionary name
    /// and the `ID` of its definition lines change (as bcftools annotate `--rename-annots`
    /// edits the dictionary). Returns `false` when `id` is not in the dictionary or
    /// `new_name` already is.
    pub fn rename_id(&mut self, id: u32, new_name: &[u8]) -> bool {
        self.dictionaries_changed();
        if self.id_map.contains_key(new_name) {
            return false;
        }
        let Some(entry) = self.ids.get_mut(id as usize).and_then(Option::as_mut) else {
            return false;
        };
        let old = std::mem::replace(&mut entry.name, new_name.into());
        let serials = entry.hrec;
        self.id_map.remove(&old);
        self.id_map.insert(new_name.into(), id);
        for serial in serials.into_iter().flatten() {
            if let Some(rec) = self.slab[serial as usize].as_mut()
                && let Some(i) = rec.find_key(b"ID")
            {
                rec.vals[i] = new_name.to_vec();
            }
        }
        for (k, tag) in SPECIAL_TAGS.iter().enumerate() {
            if **tag == *old {
                self.special[k] = None;
            }
            if **tag == *new_name {
                self.special[k] = Some(id);
            }
        }
        true
    }

    /// Renames a contig in place, keeping its numeric ID: the dictionary name and the `ID` of
    /// its `##contig` line change (as bcftools annotate `--rename-chrs` edits the header).
    /// Returns `false` when `rid` has no `##contig` line.
    pub fn rename_contig(&mut self, rid: u32, new_name: &[u8]) -> bool {
        self.dictionaries_changed();
        let Some(entry) = self.contigs.get_mut(rid as usize).and_then(Option::as_mut) else {
            return false;
        };
        let Some(serial) = entry.hrec else {
            return false;
        };
        let old = std::mem::replace(&mut entry.name, new_name.into());
        if self.contig_map.get(&old) == Some(&rid) {
            self.contig_map.remove(&old);
        }
        self.contig_map.insert(new_name.into(), rid);
        if let Some(rec) = self.slab[serial as usize].as_mut()
            && let Some(i) = rec.find_key(b"ID")
        {
            rec.vals[i] = new_name.to_vec();
        }
        true
    }

    /// Clears the dictionary's reference to a line (HTSlib's `bcf_hdr_unregister_hrec`).
    fn unregister(&mut self, serial: u32) {
        let Some(rec) = self.record_by_serial(serial) else {
            return;
        };
        let kind = rec.kind;
        let Some(id) = rec.find_key(b"ID").map(|i| rec.vals[i].clone()) else {
            return;
        };
        match kind {
            HeaderRecordType::Contig => {
                if let Some(&c) = self.contig_map.get(id.as_slice())
                    && let Some(e) = self.contigs[c as usize].as_mut()
                    && e.hrec == Some(serial)
                {
                    e.hrec = None;
                }
            }
            HeaderRecordType::Filter | HeaderRecordType::Info | HeaderRecordType::Format => {
                let slot = kind as usize;
                if let Some(&c) = self.id_map.get(id.as_slice())
                    && let Some(e) = self.ids[c as usize].as_mut()
                    && e.hrec[slot] == Some(serial)
                {
                    e.hrec[slot] = None;
                }
            }
            _ => {}
        }
    }

    /// Forgets the translation tables built from this header's dictionaries, and stamps them
    /// anew for the tables built into this header.
    fn dictionaries_changed(&mut self) {
        self.dictionaries = Dictionaries::next();
        self.translation = Translation::Unknown;
    }

    fn remove_from_gen(&mut self, serial: u32) {
        self.generic.retain(|_, s| *s != serial);
    }

    fn record_by_serial(&self, serial: u32) -> Option<&HeaderRecord> {
        self.slab.get(serial as usize).and_then(Option::as_ref)
    }

    /// Returns the dictionary's defining line of a FILTER/INFO/FORMAT/contig ID.
    fn dict_hrec(&self, kind: HeaderRecordType, name: &[u8]) -> Option<u32> {
        match kind {
            HeaderRecordType::Contig => {
                let &c = self.contig_map.get(name)?;
                self.contigs[c as usize].as_ref()?.hrec
            }
            HeaderRecordType::Filter | HeaderRecordType::Info | HeaderRecordType::Format => {
                let &c = self.id_map.get(name)?;
                self.ids[c as usize].as_ref()?.hrec[kind as usize]
            }
            _ => None,
        }
    }

    /// Finds a header line (HTSlib's `bcf_hdr_get_hrec`):
    ///
    /// - generic lines by key and, optionally, value;
    /// - structured lines of class `class` whose field `key` has `value`;
    /// - FILTER/INFO/FORMAT/contig lines by ID (`value`; `key` is ignored).
    #[must_use]
    pub fn get_record(
        &self,
        kind: HeaderRecordType,
        key: &[u8],
        value: Option<&[u8]>,
        class: Option<&[u8]>,
    ) -> Option<&HeaderRecord> {
        match kind {
            HeaderRecordType::Generic => {
                if let Some(value) = value {
                    let mut k = b"##".to_vec();
                    k.extend_from_slice(key);
                    k.push(b'=');
                    k.extend_from_slice(value);
                    return self.generic.get(&k).and_then(|&s| self.record_by_serial(s));
                }
                self.records()
                    .find(|r| r.kind == HeaderRecordType::Generic && r.key == key)
            }
            HeaderRecordType::Structured => {
                let class = class?;
                let value = value?;
                if key == b"ID" {
                    let k = gen_key_structured(class, value);
                    return self.generic.get(&k).and_then(|&s| self.record_by_serial(s));
                }
                self.records().find(|r| {
                    r.kind == HeaderRecordType::Structured
                        && r.key == class
                        && r.find_key(key).is_some_and(|j| r.vals[j] == value)
                })
            }
            _ => {
                let serial = self.dict_hrec(kind, value?)?;
                self.record_by_serial(serial)
            }
        }
    }

    /// Returns the definition of an INFO, FILTER, or FORMAT ID, or a contig.
    #[must_use]
    pub fn definition(&self, kind: HeaderRecordType, id: &[u8]) -> Option<&HeaderRecord> {
        self.get_record(kind, b"ID", Some(id), None)
    }

    /// Iterates over the header lines in order.
    pub fn records(&self) -> impl Iterator<Item = &HeaderRecord> {
        self.order
            .iter()
            .filter_map(|&s| self.slab[s as usize].as_ref())
    }

    /// Returns the number of header lines.
    #[inline]
    #[must_use]
    pub fn record_count(&self) -> usize {
        self.order.len()
    }

    /// Returns a mutable reference to the `i`-th header line. Changing the ID, type, or
    /// number of a definition through it does not update the dictionaries, as in HTSlib.
    pub fn record_mut(&mut self, i: usize) -> Option<&mut HeaderRecord> {
        let serial = *self.order.get(i)?;
        self.slab[serial as usize].as_mut()
    }

    /// Formats the header as text, ending with the `#CHROM` line (HTSlib's `bcf_hdr_format`).
    pub fn format_into(&self, out: &mut Vec<u8>, flavor: HeaderTextFlavor) {
        for r in self.records() {
            r.format_into(out, flavor);
        }
        out.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            out.extend_from_slice(b"\tFORMAT");
            for s in &self.samples {
                out.push(b'\t');
                out.extend_from_slice(s);
            }
        }
        out.push(b'\n');
    }

    /// Returns the header as VCF text.
    #[must_use]
    pub fn to_vcf_text(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.format_into(&mut out, HeaderTextFlavor::Vcf);
        out
    }

    /// Updates internal tables after edits (HTSlib's `bcf_hdr_sync`). The dictionaries are
    /// always kept current, so this exists for API parity and does nothing.
    pub fn sync(&mut self) {}

    // ---- version ----

    /// Returns the `##fileformat` value, warning and assuming `VCFv4.2` when absent
    /// (HTSlib's `bcf_hdr_get_version`).
    #[must_use]
    pub fn version_string(&self) -> Vec<u8> {
        match self.get_record(HeaderRecordType::Generic, b"fileformat", None, None) {
            Some(r) => r.value.clone().unwrap_or_default(),
            None => {
                log::warn(
                    "bcf_hdr_get_version",
                    format_args!("No version string found, assuming VCFv4.2"),
                );
                b"VCFv4.2".to_vec()
            }
        }
    }

    /// Returns the VCF version as HTSlib encodes it (`major * 1000000 + minor * 1000`), using
    /// the cached value when there is one (HTSlib's `bcf_get_version`).
    #[must_use]
    pub fn version(&self) -> i32 {
        if self.version != 0 {
            return self.version;
        }
        parse_version_string(&self.version_string())
    }

    /// Returns the cached version, 0 when no `##fileformat` line was seen (HTSlib's
    /// `aux->version`).
    #[inline]
    #[must_use]
    pub fn cached_version(&self) -> i32 {
        self.version
    }

    /// Sets the `##fileformat` value (HTSlib's `bcf_hdr_set_version`).
    pub fn set_version(&mut self, version: &[u8]) {
        let pos = self.order.iter().position(|&s| {
            self.record_by_serial(s)
                .is_some_and(|r| r.kind == HeaderRecordType::Generic && r.key == b"fileformat")
        });
        match pos {
            None => {
                // HTSlib parses the line but does not add it.
                let mut line = b"##fileformat=".to_vec();
                line.extend_from_slice(version);
                self.version = parse_version_string(version);
            }
            Some(pos) => {
                let serial = self.order[pos];
                self.remove_from_gen(serial);
                let rec = self.slab[serial as usize]
                    .as_mut()
                    .expect("ordered lines exist");
                rec.value = Some(version.to_vec());
                let key = gen_key_generic(rec);
                self.generic.insert(key, serial);
                self.version = parse_version_string(version);
            }
        }
    }

    // ---- dictionaries ----

    /// Returns the ID of a FILTER/INFO/FORMAT name (HTSlib's `bcf_hdr_id2int(BCF_DT_ID)`).
    #[inline]
    #[must_use]
    pub fn id(&self, name: &[u8]) -> Option<u32> {
        self.id_map.get(name).copied()
    }

    /// Returns the ID of END, SVLEN, LEN, or GT without a hash lookup.
    #[inline]
    pub(crate) fn special_id(&self, tag: SpecialTag) -> Option<u32> {
        self.special[tag as usize]
    }

    /// Returns the number of slots in the ID dictionary (HTSlib's `bcf_hdr_nids`).
    #[inline]
    #[must_use]
    pub fn id_count(&self) -> usize {
        self.ids.len()
    }

    /// Returns the name of a FILTER/INFO/FORMAT ID (HTSlib's `bcf_hdr_int2id(BCF_DT_ID)`).
    #[inline]
    #[must_use]
    pub fn id_name(&self, id: u32) -> Option<&[u8]> {
        self.ids.get(id as usize)?.as_ref().map(|e| &*e.name)
    }

    /// Returns the packed definition word of an ID for a kind (`info[type]`), or `None` when
    /// the slot is empty.
    #[inline]
    #[must_use]
    pub fn id_info(&self, kind: HeaderRecordType, id: u32) -> Option<u32> {
        let slot = kind.id_slot()?;
        Some(self.ids.get(id as usize)?.as_ref()?.info[slot])
    }

    /// Returns whether an ID is defined for a kind (HTSlib's `bcf_hdr_idinfo_exists`). Lines
    /// removed with [`Self::remove`] still count, as in HTSlib.
    #[inline]
    #[must_use]
    pub fn id_exists(&self, kind: HeaderRecordType, id: u32) -> bool {
        self.id_info(kind, id).is_some_and(|i| i & 0xf != 0xf)
    }

    /// Returns the value type of an ID for a kind (HTSlib's `bcf_hdr_id2type`).
    #[inline]
    #[must_use]
    pub fn id_type(&self, kind: HeaderRecordType, id: u32) -> Option<ValueType> {
        ValueType::from_code((self.id_info(kind, id)? >> 4) & 0xf)
    }

    /// Returns the raw `BCF_HT_*` type code of an ID (15 when undefined).
    #[inline]
    #[must_use]
    pub fn id_type_code(&self, kind: HeaderRecordType, id: u32) -> u32 {
        self.id_info(kind, id).map_or(15, |i| (i >> 4) & 0xf)
    }

    /// Returns the `Number=` of an ID for a kind (HTSlib's `bcf_hdr_id2length` and
    /// `bcf_hdr_id2number`).
    #[inline]
    #[must_use]
    pub fn id_number(&self, kind: HeaderRecordType, id: u32) -> Option<Number> {
        Number::from_info(self.id_info(kind, id)?)
    }

    /// Returns the definition line of an ID for a kind (HTSlib's `bcf_hdr_id2hrec`).
    #[must_use]
    pub fn id_record(&self, kind: HeaderRecordType, id: u32) -> Option<&HeaderRecord> {
        let serial = match kind {
            HeaderRecordType::Contig => self.contigs.get(id as usize)?.as_ref()?.hrec?,
            _ => self.ids.get(id as usize)?.as_ref()?.hrec[kind.id_slot()?]?,
        };
        self.record_by_serial(serial)
    }

    /// Returns the ID of a contig (HTSlib's `bcf_hdr_name2id`).
    #[inline]
    #[must_use]
    pub fn contig_id(&self, name: &[u8]) -> Option<u32> {
        self.contig_map.get(name).copied()
    }

    /// Returns the name of a contig ID (HTSlib's `bcf_hdr_id2name`), including removed ones.
    #[inline]
    #[must_use]
    pub fn contig_name(&self, rid: i32) -> Option<&[u8]> {
        let rid = usize::try_from(rid).ok()?;
        self.contigs.get(rid)?.as_ref().map(|e| &*e.name)
    }

    /// Returns the length of a contig (0 when not given).
    #[inline]
    #[must_use]
    pub fn contig_length(&self, rid: u32) -> Option<i64> {
        self.contigs.get(rid as usize)?.as_ref().map(|e| e.length)
    }

    /// Returns the number of slots in the contig dictionary.
    #[inline]
    #[must_use]
    pub fn contig_count(&self) -> usize {
        self.contigs.len()
    }

    /// Returns the defined contig names in ID order, skipping removed ones (HTSlib's
    /// `bcf_hdr_seqnames`).
    #[must_use]
    pub fn contig_names(&self) -> Vec<&[u8]> {
        self.contig_map
            .values()
            .copied()
            .filter_map(|id| {
                let e = self.contigs[id as usize].as_ref()?;
                e.hrec?;
                Some((id, &*e.name))
            })
            .collect::<std::collections::BTreeMap<_, _>>()
            .into_values()
            .collect()
    }

    /// Returns the sample names.
    #[inline]
    #[must_use]
    pub fn samples(&self) -> &[Box<[u8]>] {
        &self.samples
    }

    /// Returns the number of samples (HTSlib's `bcf_hdr_nsamples`).
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// Returns the index of a sample (HTSlib's `bcf_hdr_id2int(BCF_DT_SAMPLE)`).
    #[inline]
    #[must_use]
    pub fn sample_id(&self, name: &[u8]) -> Option<u32> {
        self.sample_map.get(name).copied()
    }

    /// Returns a copy made by formatting and reparsing (HTSlib's `bcf_hdr_dup`), which keeps
    /// `IDX` fields and drops removed dictionary entries.
    ///
    /// # Errors
    ///
    /// Returns an error when the text does not reparse.
    pub fn dup(&self) -> Result<Self, HeaderError> {
        let mut text = Vec::new();
        self.format_into(&mut text, HeaderTextFlavor::Bcf);
        Self::parse(&text)
    }

    /// Checks standard INFO and FORMAT definitions, warning once per tag about a wrong Number
    /// or Type (HTSlib's `bcf_hdr_check_sanity`).
    fn check_sanity(&self) {
        let version = self.version();
        for (i, t) in SANITY_INFO.iter().enumerate() {
            self.check_sanity_tag(HeaderRecordType::Info, t, &INFO_WARNED[i], version);
        }
        for (i, t) in SANITY_FMT.iter().enumerate() {
            self.check_sanity_tag(HeaderRecordType::Format, t, &FMT_WARNED[i], version);
        }
    }

    fn check_sanity_tag(
        &self,
        kind: HeaderRecordType,
        t: &SanityTag,
        warned: &'static AtomicBool,
        version: i32,
    ) {
        use std::sync::atomic::Ordering;
        if warned.load(Ordering::Relaxed) {
            return;
        }
        let Some(id) = self.id(t.name.as_bytes()) else {
            return;
        };
        if !self.id_exists(kind, id) {
            return;
        }
        let info = self.id_info(kind, id).expect("existing IDs have info");
        let length = (info >> 8) & 0xf;
        let number = info >> 12;
        let mut flag = false;
        if kind == HeaderRecordType::Info {
            flag = (length != t.number && length != 1)
                || (length == 0 && i64::from(number) != number::atoi(t.number_str.as_bytes()));
        } else if length != t.number {
            if (version < VCF44 && length != 1) || (version >= VCF44 && version >= t.version) {
                flag = true;
            }
        } else if length == 0 && i64::from(number) != number::atoi(t.number_str.as_bytes()) {
            flag = true;
        }
        if flag {
            log::warn(
                "bcf_hdr_check_sanity",
                format_args!("{} should be declared as Number={}", t.name, t.number_str),
            );
        }
        let ty = (info >> 4) & 0xf;
        if ty != t.ty {
            log::warn(
                "bcf_hdr_check_sanity",
                format_args!(
                    "{} should be declared as Type={}",
                    t.name,
                    ["Flag", "Integer", "Float", "String"][t.ty as usize]
                ),
            );
            flag = true;
        }
        if flag {
            warned.store(true, Ordering::Relaxed);
        }
    }

    /// Returns whether `id` is defined as a FILTER/INFO/FORMAT of `kind` with the dictionary's
    /// `info` slot (for internal checks like `info[BCF_HL_INFO] == 15`).
    #[inline]
    pub(crate) fn id_slot_undefined(&self, id: u32, slot: usize) -> bool {
        self.ids
            .get(id as usize)
            .and_then(Option::as_ref)
            .is_none_or(|e| e.info[slot] == 15)
    }
}

#[derive(Clone, Copy)]
enum Dict {
    Id,
    Contig,
}

/// Parses an `IDX=` value as HTSlib does (`strtol` with nothing after the digits).
fn parse_idx(v: &[u8]) -> Option<i64> {
    let (idx, end) = number::strtol(v, 0);
    if end != v.len() || idx < 0 || idx >= i64::from(i32::MAX) - 1 {
        return None;
    }
    Some(idx)
}

/// C's `sscanf(s, "%d", &n) == 1`.
fn scanf_int(v: &[u8]) -> Option<i32> {
    let (n, end) = number::strtol(v, 0);
    if end == 0 {
        return None;
    }
    Some(n as i32)
}

/// Appends `IDX=id` (HTSlib's `hrec_add_idx`).
fn add_idx(hrec: &mut HeaderRecord, id: u32) {
    hrec.keys.push(b"IDX".to_vec());
    hrec.vals.push(id.to_string().into_bytes());
}

fn gen_key_generic(hrec: &HeaderRecord) -> Vec<u8> {
    let mut k = b"##".to_vec();
    k.extend_from_slice(&hrec.key);
    k.push(b'=');
    if let Some(v) = &hrec.value {
        k.extend_from_slice(v);
    } else {
        k.extend_from_slice(b"(null)");
    }
    k
}

fn gen_key_structured(class: &[u8], id: &[u8]) -> Vec<u8> {
    let mut k = b"##".to_vec();
    k.extend_from_slice(class);
    k.extend_from_slice(b"=<ID=");
    k.extend_from_slice(id);
    k.push(b'>');
    k
}

/// Parses `VCFv4.3` into 4003000 (HTSlib's `bcf_get_version`), warning and returning 4.2 for a
/// malformed string.
pub(crate) fn parse_version_string(s: &[u8]) -> i32 {
    let fail = || {
        log::warn(
            "bcf_get_version",
            format_args!("Couldn't get VCF version, considering as 4.2"),
        );
        VCF_DEFAULT_VERSION
    };
    let Some(pos) = s.windows(4).position(|w| w == b"VCFv") else {
        return fail();
    };
    let major = &s[pos + 4..];
    let Some(dot) = major.iter().position(|&c| c == b'.') else {
        return fail();
    };
    let (maj, _) = number::strtol(major, 0);
    let (min, _) = number::strtol(&major[dot + 1..], 0);
    (maj as i32)
        .wrapping_mul(100 * 10000)
        .wrapping_add((min as i32).wrapping_mul(1000))
}

/// Warns about invalid IDs in contig, INFO and FORMAT lines (HTSlib's `bcf_hrec_check`).
fn check_record(hrec: &HeaderRecord) {
    let kind = hrec.kind;
    if !matches!(
        kind,
        HeaderRecordType::Contig | HeaderRecordType::Info | HeaderRecordType::Format
    ) {
        return;
    }
    let Some(i) = hrec.find_key(b"ID") else {
        log::warn(
            "bcf_hrec_check",
            format_args!("Missing ID attribute in one or more header lines"),
        );
        return;
    };
    let val = &hrec.vals[i];
    let first = val.first().copied().unwrap_or(0);
    let bad = if kind == HeaderRecordType::Contig {
        first == b'*'
            || first == b'='
            || !valid_ctg(first)
            || !val[1..].iter().all(|&c| valid_ctg(c))
    } else {
        if kind == HeaderRecordType::Info && val.as_slice() == b"1000G" {
            return;
        }
        first == b'.'
            || first.is_ascii_digit()
            || !valid_tag(first)
            || !val[1..].iter().all(|&c| valid_tag(c))
    };
    if bad {
        let what = if kind == HeaderRecordType::Contig {
            "contig"
        } else {
            "tag"
        };
        log::warn(
            "bcf_hrec_check",
            format_args!(
                "Invalid {what} name: \"{}\"",
                String::from_utf8_lossy(
                    &val[..val.iter().position(|&c| c == 0).unwrap_or(val.len())]
                )
            ),
        );
    }
}

fn valid_ctg(c: u8) -> bool {
    c.is_ascii_alphanumeric() || b"!#$%&*+./:;=?@^_|~-".contains(&c)
}

fn valid_tag(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

struct SanityTag {
    name: &'static str,
    number_str: &'static str,
    number: u32,
    version: i32,
    ty: u32,
}

const fn tag(
    name: &'static str,
    number_str: &'static str,
    number: u32,
    version: i32,
    ty: u32,
) -> SanityTag {
    SanityTag {
        name,
        number_str,
        number,
        version,
        ty,
    }
}

const V42: i32 = VCF_DEFAULT_VERSION;
static SANITY_INFO: [SanityTag; 21] = [
    tag("AD", "R", 4, V42, 1),
    tag("ADF", "R", 4, V42, 1),
    tag("ADR", "R", 4, V42, 1),
    tag("AC", "A", 2, V42, 1),
    tag("AF", "A", 2, V42, 2),
    tag("CIGAR", "A", 2, V42, 3),
    tag("AA", "1", 0, V42, 3),
    tag("AN", "1", 0, V42, 1),
    tag("BQ", "1", 0, V42, 2),
    tag("DB", "0", 0, V42, 0),
    tag("DP", "1", 0, V42, 1),
    tag("END", "1", 0, V42, 1),
    tag("H2", "0", 0, V42, 0),
    tag("H3", "0", 0, V42, 0),
    tag("MQ", "1", 0, V42, 2),
    tag("MQ0", "1", 0, V42, 1),
    tag("NS", "1", 0, V42, 1),
    tag("SB", "4", 0, V42, 1),
    tag("SOMATIC", "0", 0, V42, 0),
    tag("VALIDATED", "0", 0, V42, 0),
    tag("1000G", "0", 0, V42, 0),
];
static SANITY_FMT: [SanityTag; 28] = [
    tag("AD", "R", 4, V42, 1),
    tag("ADF", "R", 4, V42, 1),
    tag("ADR", "R", 4, V42, 1),
    tag("EC", "A", 2, V42, 1),
    tag("GL", "G", 3, V42, 2),
    tag("GP", "G", 3, V42, 2),
    tag("PL", "G", 3, V42, 1),
    tag("PP", "G", 3, V42, 1),
    tag("DP", "1", 0, V42, 1),
    tag("LEN", "1", 0, V42, 1),
    tag("FT", "1", 0, V42, 3),
    tag("GQ", "1", 0, V42, 1),
    tag("GT", "1", 0, V42, 3),
    tag("HQ", "2", 0, V42, 1),
    tag("MQ", "1", 0, V42, 1),
    tag("PQ", "1", 0, V42, 1),
    tag("PS", "1", 0, V42, 1),
    tag("PSL", "P", 5, VCF44, 3),
    tag("PSO", "P", 5, VCF44, 1),
    tag("PSQ", "P", 5, VCF44, 1),
    tag("LGL", "LG", 7, VCF45, 1),
    tag("LGP", "LG", 7, VCF45, 1),
    tag("LPL", "LG", 7, VCF45, 1),
    tag("LPP", "LG", 7, VCF45, 1),
    tag("LEC", "LA", 6, VCF45, 1),
    tag("LAD", "LR", 8, VCF45, 1),
    tag("LADF", "LR", 8, VCF45, 1),
    tag("LADR", "LR", 8, VCF45, 1),
];
#[expect(
    clippy::declare_interior_mutable_const,
    reason = "array initializer for statics"
)]
const UNWARNED: AtomicBool = AtomicBool::new(false);
static INFO_WARNED: [AtomicBool; 21] = [UNWARNED; 21];
static FMT_WARNED: [AtomicBool; 28] = [UNWARNED; 28];

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=q10,Description=\"Low quality\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n";

    #[test]
    fn parse_and_format() {
        let h = Header::parse(TEXT.as_bytes()).unwrap();
        let mut want = TEXT.replace(
            "##FILTER=<ID=q10",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">\n##FILTER=<ID=q10",
        );
        assert_eq!(String::from_utf8(h.to_vcf_text()).unwrap(), want);
        assert_eq!(h.id(b"PASS"), Some(0));
        assert_eq!(h.id(b"q10"), Some(1));
        assert_eq!(h.id(b"DP"), Some(2));
        assert_eq!(h.id(b"GT"), Some(3));
        assert!(h.id_exists(HeaderRecordType::Info, 2));
        assert!(h.id_exists(HeaderRecordType::Format, 2));
        assert!(!h.id_exists(HeaderRecordType::Info, 3));
        assert_eq!(h.contig_id(b"1"), Some(0));
        assert_eq!(h.contig_length(0), Some(1000));
        assert_eq!(h.sample_count(), 2);
        let mut bcf = Vec::new();
        h.format_into(&mut bcf, HeaderTextFlavor::Bcf);
        assert!(String::from_utf8(bcf).unwrap().contains("Depth\",IDX=2>"));
        want.clear();
    }

    #[test]
    fn renames_keep_ids() {
        let mut h = Header::parse(TEXT.as_bytes()).unwrap();
        assert!(h.rename_id(1, b"lowq"));
        assert_eq!(h.id(b"lowq"), Some(1));
        assert_eq!(h.id(b"q10"), None);
        assert!(!h.rename_id(1, b"DP"), "an existing name is refused");
        assert!(h.rename_contig(0, b"chr1"));
        assert_eq!(h.contig_id(b"chr1"), Some(0));
        assert_eq!(h.contig_name(0), Some(b"chr1".as_slice()));
        let text = String::from_utf8(h.to_vcf_text()).unwrap();
        assert!(text.contains("##FILTER=<ID=lowq,"));
        assert!(text.contains("##contig=<ID=chr1,length=1000>"));
        h.remove(HeaderRecordType::Info, Some(b"DP"));
        assert!(h.id_exists(HeaderRecordType::Info, 2));
        h.undefine_id(HeaderRecordType::Info, b"DP");
        assert!(!h.id_exists(HeaderRecordType::Info, 2));
        assert!(h.id_exists(HeaderRecordType::Format, 2));
    }

    #[test]
    fn idx_is_kept() {
        let text = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=5>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let h = Header::parse(text.as_bytes()).unwrap();
        assert_eq!(h.id(b"DP"), Some(5));
        assert_eq!(h.id_count(), 6);
        assert_eq!(h.id_name(3), None);
    }

    #[test]
    fn remove_keeps_dictionary() {
        let mut h = Header::parse(TEXT.as_bytes()).unwrap();
        h.remove(HeaderRecordType::Info, Some(b"DP"));
        assert!(h.definition(HeaderRecordType::Info, b"DP").is_none());
        assert!(h.id_exists(HeaderRecordType::Info, 2));
        h.remove(HeaderRecordType::Format, None);
        assert!(!h.to_vcf_text().windows(8).any(|w| w == b"##FORMAT"));
    }
}
