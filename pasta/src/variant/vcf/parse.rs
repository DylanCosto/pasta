//! Parsing VCF lines into records (HTSlib's `vcf_parse` and its helpers), with every warning,
//! error, and quirk of HTSlib 1.24.

use std::{cell::RefCell, sync::atomic::AtomicBool};

use crate::variant::{
    Header, Record,
    encoding::{
        self, BcfType, INT32_MISSING, MAX_INT32, MIN_INT32, enc_int1, enc_size, enc_vchar,
        enc_vfloat, enc_vint,
    },
    header::{ParsedLine, parse_line as parse_header_line},
    log, number,
    record::{IS_64BIT, RecordError, rlen::get_rlen},
};

/// An error parsing a VCF line.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
#[non_exhaustive]
pub enum ParseError {
    /// The line is malformed (the reason was reported on stderr as HTSlib reports it).
    Invalid,
    /// The line uses a contig or tag missing from the header, and the header could not be
    /// updated because it was shared ([`parse_shared`]).
    HeaderUpdateNeeded,
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Invalid => "VCF parse error",
            Self::HeaderUpdateNeeded => "the line needs a header update",
        })
    }
}

impl std::error::Error for ParseError {}

/// Header access during parsing: mutable headers gain dummy definitions for undefined tags.
pub(super) trait HeaderAccess {
    fn get(&self) -> &Header;
    /// Returns whether undefined tags can be added.
    fn can_update(&self) -> bool;
    /// Adds a dummy header line; `Err` when the header cannot change.
    fn add_dummy(&mut self, line: &[u8]) -> Result<(), ParseError>;
}

impl HeaderAccess for &mut Header {
    #[inline]
    fn get(&self) -> &Header {
        self
    }

    #[inline]
    fn can_update(&self) -> bool {
        true
    }

    fn add_dummy(&mut self, line: &[u8]) -> Result<(), ParseError> {
        if let ParsedLine::Record(r, _) = parse_header_line(line) {
            // Errors (conflicting IDX) leave the tag undefined; the caller reports it.
            let _ = self.add_record(r);
        }
        Ok(())
    }
}

struct Shared<'h>(&'h Header);

impl HeaderAccess for Shared<'_> {
    #[inline]
    fn get(&self) -> &Header {
        self.0
    }

    #[inline]
    fn can_update(&self) -> bool {
        false
    }

    fn add_dummy(&mut self, _line: &[u8]) -> Result<(), ParseError> {
        Err(ParseError::HeaderUpdateNeeded)
    }
}

/// Parses a VCF line (without its newline) into `rec` (HTSlib's `vcf_parse`). Contigs and
/// tags missing from the header are added to it with HTSlib's dummy definitions and
/// warnings, and recorded in [`Record::errors`].
///
/// The line buffer is used as scratch space: its contents are changed.
///
/// # Errors
///
/// Returns [`ParseError::Invalid`] for a malformed line, after reporting why on stderr.
pub fn parse(header: &mut Header, line: &mut Vec<u8>, rec: &mut Record) -> Result<(), ParseError> {
    let mut access = header;
    parse_with(&mut access, line, rec)
}

/// Parses a VCF line like [`parse`] without changing the header, for parsing on worker
/// threads: a line that needs a header update fails with [`ParseError::HeaderUpdateNeeded`]
/// (before reporting anything about it), so it can be reparsed with [`parse`] in order.
///
/// # Errors
///
/// Returns [`ParseError::HeaderUpdateNeeded`] or [`ParseError::Invalid`].
pub fn parse_shared(
    header: &Header,
    line: &mut Vec<u8>,
    rec: &mut Record,
) -> Result<(), ParseError> {
    parse_with(&mut Shared(header), line, rec)
}

/// Position of a tab-separated token, as HTSlib's `kstrtok` finds it.
struct Tokens {
    next: usize,
    finished: bool,
}

impl Tokens {
    /// Returns the next token's start and end; the end byte is set to NUL.
    #[inline]
    fn next(&mut self, s: &mut [u8]) -> Option<(usize, usize)> {
        if self.finished {
            return None;
        }
        let start = self.next;
        let end = start + memchr::memchr2(b'\t', 0, &s[start..]).unwrap_or(s.len() - start);
        if s[end] == 0 {
            self.finished = true;
        }
        s[end] = 0;
        self.next = end + 1;
        Some((start, end))
    }
}

#[inline]
fn not_dot(s: &[u8], p: usize) -> bool {
    !(s[p] == b'.' && s[p + 1] == 0)
}

/// The C string starting at `p`.
#[inline]
pub(super) fn cstr(s: &[u8], p: usize) -> &[u8] {
    let rest = &s[p..];
    &rest[..memchr::memchr(0, rest).unwrap_or(rest.len())]
}

pub(super) fn lossy(s: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(s)
}

static MISSING_HEADERS: &str = "Missing headers may cause later processing to fail";

fn parse_with<H: HeaderAccess>(
    h: &mut H,
    line: &mut Vec<u8>,
    v: &mut Record,
) -> Result<(), ParseError> {
    let len = line.len();
    // Room for overruns, as HTSlib's 4 NUL bytes past the line.
    line.extend_from_slice(&[0, 0, 0, 0]);
    let s = &mut line[..];
    let result = parse_fields(h, s, len, v);
    line.truncate(len);
    result
}

fn parse_fields<H: HeaderAccess>(
    h: &mut H,
    s: &mut [u8],
    len: usize,
    v: &mut Record,
) -> Result<(), ParseError> {
    v.clear();
    let mut tok = Tokens {
        next: 0,
        finished: false,
    };

    // CHROM
    let (p, _) = tok.next(s).ok_or(ParseError::Invalid)?;
    let chrom = cstr(s, p);
    v.rid = match cached_contig(h.get(), chrom) {
        Some(id) => id as i32,
        None => {
            if !h.can_update() {
                return Err(ParseError::HeaderUpdateNeeded);
            }
            log::warn(
                "vcf_parse",
                format_args!(
                    "Contig '{}' is not defined in the header. (Quick workaround: index the file with tabix.)",
                    lossy(chrom)
                ),
            );
            log::warn("vcf_parse", format_args!("{MISSING_HEADERS}"));
            v.errcode = RecordError::CTG_UNDEF;
            let mut dummy = b"##contig=<ID=".to_vec();
            dummy.extend_from_slice(chrom);
            dummy.push(b'>');
            h.add_dummy(&dummy)?;
            match h.get().contig_id(chrom) {
                Some(id) => id as i32,
                None => {
                    log::error(
                        "vcf_parse",
                        format_args!("Could not add dummy header for contig '{}'", lossy(chrom)),
                    );
                    v.errcode |= RecordError::CTG_INVALID;
                    return Err(ParseError::Invalid);
                }
            }
        }
    };

    // POS
    let (p, _) = tok.next(s).ok_or(ParseError::Invalid)?;
    let (pos, end, overflow) = number::str2uint(s, p, 62);
    if overflow {
        log::error(
            "vcf_parse",
            format_args!("Position value '{}' is too large", lossy(cstr(s, p))),
        );
        return Err(ParseError::Invalid);
    } else if s[end] != 0 {
        log::error(
            "vcf_parse",
            format_args!("Could not parse the position '{}'", lossy(cstr(s, p))),
        );
        return Err(ParseError::Invalid);
    }
    v.pos = pos as i64 - 1;
    if v.pos >= i64::from(i32::MAX) {
        v.unpacked |= IS_64BIT;
    }

    // ID
    let (p, q) = tok.next(s).ok_or(ParseError::Invalid)?;
    if not_dot(s, p) {
        enc_vchar(&mut v.shared, &s[p..q]);
    } else {
        enc_size(&mut v.shared, 0, BcfType::Char as u8);
    }

    // REF
    let (p, q) = tok.next(s).ok_or(ParseError::Invalid)?;
    enc_vchar(&mut v.shared, &s[p..q]);
    v.n_allele = 1;
    v.rlen = (q - p) as i64;

    // ALT
    let (p, q) = tok.next(s).ok_or(ParseError::Invalid)?;
    if not_dot(s, p) {
        let mut t = p;
        let mut r = p;
        loop {
            if s[r] == b',' || s[r] == 0 {
                if v.n_allele == u32::from(u16::MAX) {
                    log::error(
                        "vcf_parse",
                        format_args!(
                            "Too many ALT alleles at {}:{}",
                            lossy(h.get().contig_name(v.rid).unwrap_or(b"(unknown)")),
                            v.pos + 1
                        ),
                    );
                    v.errcode |= RecordError::LIMITS;
                    return Err(ParseError::Invalid);
                }
                enc_vchar(&mut v.shared, &s[t..r]);
                t = r + 1;
                v.n_allele += 1;
            }
            if r == q {
                break;
            }
            r += 1;
        }
    }

    // QUAL
    let (p, _) = tok.next(s).ok_or(ParseError::Invalid)?;
    if not_dot(s, p) {
        v.qual = number::strtod(s, p).0 as f32;
    } else {
        v.qual = encoding::float_missing();
    }
    let max_unpack = v.max_unpack.bits();
    if max_unpack != 0 && max_unpack >> 1 == 0 {
        return finish(h, v);
    }

    // FILTER
    let (p, q) = tok.next(s).ok_or(ParseError::Invalid)?;
    if not_dot(s, p) {
        parse_filter(h, s, v, p, q)?;
    } else {
        enc_vint(&mut v.shared, &[], None);
    }
    if max_unpack != 0 && max_unpack >> 2 == 0 {
        return finish(h, v);
    }

    // INFO
    let (p, q) = tok.next(s).ok_or(ParseError::Invalid)?;
    if not_dot(s, p) {
        parse_info(h, s, v, p, q)?;
    }
    if max_unpack != 0 && max_unpack >> 3 == 0 {
        return finish(h, v);
    }

    // FORMAT; optional
    if let Some((p, q)) = tok.next(s) {
        super::format_fields::parse_format(h, s, len, v, p, q)?;
    }
    finish(h, v)
}

#[inline]
fn finish<H: HeaderAccess>(h: &H, v: &mut Record) -> Result<(), ParseError> {
    v.rlen = get_rlen(h.get(), v);
    Ok(())
}

pub(super) fn record_location(h: &Header, v: &Record) -> String {
    format!(
        "{}:{}",
        lossy(h.contig_name(v.rid).unwrap_or(b"(unknown)")),
        v.pos + 1
    )
}

/// Warns about an undefined tag and adds its dummy definition (shared by FILTER, INFO, and
/// FORMAT). Returns the tag's ID, or an error.
pub(super) fn define_dummy<H: HeaderAccess>(
    h: &mut H,
    v: &mut Record,
    kind: &str,
    name: &[u8],
    warning_suffix: &str,
    dummy: &[u8],
    slot: usize,
) -> Result<u32, ParseError> {
    if !h.can_update() {
        return Err(ParseError::HeaderUpdateNeeded);
    }
    let func = match kind {
        "FILTER" => "vcf_parse_filter",
        "INFO" => "vcf_parse_info",
        _ => "vcf_parse_format_dict2",
    };
    log::warn(
        func,
        format_args!(
            "{kind} '{}' at {} is not defined in the header{warning_suffix}",
            lossy(name),
            record_location(h.get(), v)
        ),
    );
    if !v
        .errcode
        .intersects(RecordError::TAG_UNDEF | RecordError::CTG_UNDEF)
    {
        log::warn(func, format_args!("{MISSING_HEADERS}"));
    }
    h.add_dummy(dummy)?;
    v.errcode |= RecordError::TAG_UNDEF;
    let id = h.get().id(name);
    let defined = id.is_some_and(|id| slot == usize::MAX || !h.get().id_slot_undefined(id, slot));
    match id {
        Some(id) if defined => Ok(id),
        _ => {
            log::error(
                func,
                format_args!(
                    "Could not add dummy header for {kind} '{}' at {}",
                    lossy(name),
                    record_location(h.get(), v)
                ),
            );
            v.errcode |= RecordError::TAG_INVALID;
            Err(ParseError::Invalid)
        }
    }
}

fn parse_filter<H: HeaderAccess>(
    h: &mut H,
    s: &mut [u8],
    v: &mut Record,
    p: usize,
    q: usize,
) -> Result<(), ParseError> {
    if q > 0 && s[q - 1] == b';' {
        s[q - 1] = 0;
    }
    SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        let flt = &mut scratch.flt;
        flt.clear();
        let mut t = p;
        loop {
            let end = t + memchr::memchr2(b';', 0, &s[t..]).unwrap_or(s.len() - t);
            let last = s[end] == 0;
            s[end] = 0;
            let name = &s[t..end];
            let id = match h.get().id(name) {
                Some(id) => id,
                None => {
                    let mut dummy = b"##FILTER=<ID=".to_vec();
                    dummy.extend_from_slice(name);
                    dummy.extend_from_slice(b",Description=\"Dummy\">");
                    let name = name.to_vec();
                    define_dummy(h, v, "FILTER", &name, "", &dummy, usize::MAX)?
                }
            };
            flt.push(id as i32);
            if last {
                break;
            }
            t = end + 1;
        }
        enc_vint(&mut v.shared, flt, None);
        Ok(())
    })
}

static EXTREME_INT_WARNED: AtomicBool = AtomicBool::new(false);
static NEGATIVE_RLEN_WARNED: AtomicBool = AtomicBool::new(false);

fn parse_info<H: HeaderAccess>(
    h: &mut H,
    s: &mut [u8],
    v: &mut Record,
    p: usize,
    q: usize,
) -> Result<(), ParseError> {
    v.n_info = 0;
    if q > 0 && s[q - 1] == b';' {
        s[q - 1] = 0;
    }
    let mut key = p;
    let mut r = p;
    loop {
        r += memchr::memchr3(b';', b'=', 0, &s[r..]).unwrap_or(s.len() - r);
        if v.n_info == u32::from(u16::MAX) {
            log::error(
                "vcf_parse_info",
                format_args!("Too many INFO entries at {}", record_location(h.get(), v)),
            );
            v.errcode |= RecordError::LIMITS;
            return Err(ParseError::Invalid);
        }
        let key_end = r;
        let mut c = s[r];
        s[r] = 0;
        let (val, end) = if c == b'=' {
            let vstart = r + 1;
            let e = vstart + memchr::memchr2(b';', 0, &s[vstart..]).unwrap_or(s.len() - vstart);
            c = s[e];
            s[e] = 0;
            (Some(vstart), e)
        } else {
            (None, r)
        };
        if s[key] == 0 {
            // Faulty VCF, ";;" in the INFO
            if c == 0 {
                break;
            }
            r = end;
            key = r + 1;
            r += 1;
            continue;
        }
        let key_bytes = &s[key..key_end];
        let id = match cached_id(h.get(), CacheSlot::Info, v.n_info as usize, key_bytes) {
            Some(id) if !h.get().id_slot_undefined(id, 1) => id,
            _ => {
                let name = key_bytes.to_vec();
                let mut dummy = b"##INFO=<ID=".to_vec();
                dummy.extend_from_slice(&name);
                dummy.extend_from_slice(b",Number=1,Type=String,Description=\"Dummy\">");
                define_dummy(h, v, "INFO", &name, ", assuming Type=String", &dummy, 1)?
            }
        };
        let y = h
            .get()
            .id_info(crate::variant::HeaderRecordType::Info, id)
            .unwrap_or(15);
        v.n_info += 1;
        enc_int1(&mut v.shared, id as i32);
        match val {
            None => enc_size(&mut v.shared, 0, BcfType::Null as u8),
            Some(vs) => {
                let ty = (y >> 4) & 0xf;
                if ty == 0 || ty == 3 {
                    enc_vchar(&mut v.shared, &s[vs..end]);
                } else {
                    parse_info_numbers(h.get(), s, v, key, vs, end, ty)?;
                }
            }
        }
        if c == 0 {
            break;
        }
        r = end;
        key = r + 1;
        r += 1;
    }
    Ok(())
}

fn parse_info_numbers(
    h: &Header,
    s: &[u8],
    v: &mut Record,
    key: usize,
    val: usize,
    end: usize,
    ty: u32,
) -> Result<(), ParseError> {
    let n_val = 1 + memchr::memchr_iter(b',', &s[val..end]).count();
    SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        if ty == 1 {
            let a = &mut scratch.ints;
            a.clear();
            let mut t = val;
            for _ in 0..n_val {
                let (mut tmp, te, overflow) = number::str2int(s, t, 64);
                if te == t {
                    tmp = i64::from(INT32_MISSING);
                } else if overflow || tmp < i64::from(MIN_INT32) || tmp > i64::from(MAX_INT32) {
                    if !log::fired(&EXTREME_INT_WARNED) {
                        log::warn_once(
                            &EXTREME_INT_WARNED,
                            "vcf_parse_info",
                            format_args!(
                                "Extreme INFO/{} value encountered and set to missing at {}",
                                lossy(cstr(s, key)),
                                record_location(h, v)
                            ),
                        );
                    }
                    tmp = i64::from(INT32_MISSING);
                }
                a.push(tmp as i32);
                t = te;
                while s[t] != 0 && s[t] != b',' {
                    t += 1;
                }
                t += 1;
            }
            if n_val == 1 {
                let val1 = a[0];
                enc_int1(&mut v.shared, val1);
                if val1 != INT32_MISSING
                    && cstr(s, key) == b"END"
                    && i64::from(val1) <= v.pos
                    && !log::fired(&NEGATIVE_RLEN_WARNED)
                {
                    {
                        log::warn_once(
                            &NEGATIVE_RLEN_WARNED,
                            "vcf_parse_info",
                            format_args!(
                                "INFO/END={val1} is smaller than POS at {}",
                                record_location(h, v)
                            ),
                        );
                    }
                }
            } else {
                enc_vint(&mut v.shared, a, None);
            }
        } else if ty == 2 {
            let f = &mut scratch.floats;
            f.clear();
            let mut t = val;
            for _ in 0..n_val {
                let (d, te, overflow) = number::str2dbl(s, t);
                f.push(if te == t || overflow {
                    encoding::float_missing()
                } else {
                    d as f32
                });
                t = te;
                while s[t] != 0 && s[t] != b',' {
                    t += 1;
                }
                t += 1;
            }
            enc_vfloat(&mut v.shared, f);
        }
        Ok(())
    })
}

/// Reusable per-thread parsing buffers (HTSlib keeps them in the header's `mem`).
#[derive(Default)]
pub(super) struct Scratch {
    pub(super) flt: Vec<i32>,
    pub(super) ints: Vec<i32>,
    pub(super) floats: Vec<f32>,
    pub(super) fmt: Vec<super::format_fields::FmtAux>,
}

thread_local! {
    pub(super) static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
    static KEY_CACHE: RefCell<KeyCache> = RefCell::new(KeyCache::default());
}

/// The IDs found for the previous line's keys, by position, to skip hash lookups when lines
/// repeat their keys. A cached ID is used only when the header maps it back to the same name.
#[derive(Default)]
struct KeyCache {
    chrom: Option<u32>,
    info: Vec<u32>,
    format: Vec<u32>,
}

/// Looks up a FILTER/INFO/FORMAT name, trying the ID cached for this position first.
#[inline]
fn cached_id(header: &Header, which: CacheSlot, pos: usize, name: &[u8]) -> Option<u32> {
    KEY_CACHE.with(|c| {
        let mut c = c.borrow_mut();
        let cache = match which {
            CacheSlot::Info => &mut c.info,
            CacheSlot::Format => &mut c.format,
        };
        if let Some(&id) = cache.get(pos)
            && header.id_name(id) == Some(name)
        {
            return Some(id);
        }
        let id = header.id(name)?;
        if cache.len() <= pos {
            cache.resize(pos + 1, u32::MAX);
        }
        cache[pos] = id;
        Some(id)
    })
}

/// Which key list [`cached_id`] uses.
#[derive(Clone, Copy)]
pub(super) enum CacheSlot {
    Info,
    Format,
}

/// Looks up a contig name, trying the previous line's contig first.
#[inline]
fn cached_contig(header: &Header, name: &[u8]) -> Option<u32> {
    KEY_CACHE.with(|c| {
        let mut c = c.borrow_mut();
        if let Some(id) = c.chrom
            && header.contig_name(id as i32) == Some(name)
            && header.contigs[id as usize]
                .as_ref()
                .is_some_and(|e| e.hrec.is_some())
        {
            return Some(id);
        }
        let id = header.contig_id(name)?;
        c.chrom = Some(id);
        Some(id)
    })
}

/// [`cached_id`] for the FORMAT parser.
#[inline]
pub(super) fn cached_format_id(header: &Header, pos: usize, name: &[u8]) -> Option<u32> {
    cached_id(header, CacheSlot::Format, pos, name)
}
