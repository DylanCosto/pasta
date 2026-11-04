//! Parsing the FORMAT and sample columns (HTSlib's `vcf_parse_format` and its numbered steps:
//! `empty1`, `dict2`, `max3`, `alloc4`, `fill5`, `gt6`, `check7`).

use std::sync::atomic::AtomicBool;

use super::parse::{HeaderAccess, ParseError, SCRATCH, define_dummy, lossy, record_location};
use crate::variant::{
    HeaderRecordType, Record, VCF44,
    encoding::{
        self, BcfType, INT32_MISSING, INT32_VECTOR_END, MAX_INT32, MIN_INT32, enc_int1, enc_size,
        enc_vint,
    },
    log,
    number::{self, at},
    record::RecordError,
};

/// The most FORMAT fields a record can have (HTSlib's `MAX_N_FMT`).
const MAX_N_FMT: usize = 255;

/// The per-field state of FORMAT parsing (HTSlib's `fmt_aux_t`), with reusable buffers.
#[derive(Debug, Default)]
pub(crate) struct FmtAux {
    key: i32,
    max_m: usize,
    /// Bytes per sample, or `None` for a field to be dropped (duplicate or too large).
    size: Option<usize>,
    is_gt: bool,
    max_g: usize,
    max_l: usize,
    y: u32,
    ints: Vec<i32>,
    floats: Vec<f32>,
    chars: Vec<u8>,
}

impl FmtAux {
    #[inline]
    fn htype(&self) -> u32 {
        (self.y >> 4) & 0xf
    }
}

static EXCESSIVE_WARNED: AtomicBool = AtomicBool::new(false);
static DUPLICATE_WARNED: AtomicBool = AtomicBool::new(false);
static EXTREME_VAL_WARNED: AtomicBool = AtomicBool::new(false);

/// Parses the FORMAT column `s[p..q]` and the sample columns after it, up to `end`.
pub(super) fn parse_format<H: HeaderAccess>(
    h: &mut H,
    s: &mut [u8],
    end: usize,
    v: &mut Record,
    p: usize,
    q: usize,
) -> Result<(), ParseError> {
    let nsamples = h.get().sample_count();
    if nsamples == 0 {
        return Ok(());
    }

    // empty1: detect FORMAT "."
    if q >= end {
        log::error(
            "vcf_parse_format_empty1",
            format_args!(
                "FORMAT column with no sample columns starting at {}",
                record_location(h.get(), v)
            ),
        );
        v.errcode |= RecordError::NCOLS;
        // HTSlib's caller maps this failure to success.
        return Ok(());
    }
    v.n_fmt = 0;
    if s[p] == b'.' && s[p + 1] == 0 {
        v.n_sample = nsamples as u32;
        return Ok(());
    }

    SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        let fmt = &mut scratch.fmt;
        let n_fmt = dict2(h, s, v, p, fmt)?;
        let fmt = &mut fmt[..n_fmt];
        if let [gt] = &*fmt
            && gt.is_gt
            && gt.htype() == 3
            && h.get().keep_samples.is_none()
            && diploid_gt_only(&s[q + 1..end], nsamples, gt.key, v)
        {
            return Ok(());
        }
        max3(h, s, end, v, q, fmt)?;
        alloc4(h, v, fmt)?;
        fill5(h, s, end, v, q, fmt)?;
        gt6(v, fmt);
        check7(h, v)
    })
}

/// The common case of a FORMAT of only GT whose samples are all unprefixed diploid genotypes of
/// single digits (`0|1`, `1/0`): encodes the sample block in one pass, exactly as the general
/// steps (`max3` to `gt6`) would, as two int8 values per sample. Returns `false`, leaving the
/// record as it was, for anything else.
fn diploid_gt_only(samples: &[u8], nsamples: usize, key: i32, v: &mut Record) -> bool {
    // "a|b" per sample with a tab between samples; `alloc4`'s memory limit.
    if samples.len() + 1 != 4 * nsamples || 8 * nsamples as u64 > i32::MAX as u64 {
        return false;
    }
    let start = v.indiv.len();
    enc_int1(&mut v.indiv, key);
    enc_size(&mut v.indiv, 2, BcfType::Int8 as u8);
    v.indiv.reserve(2 * nsamples);
    for gt in samples.chunks(4) {
        let (a, sep, b) = (gt[0], gt[1], gt[2]);
        let phased = u8::from(sep == b'|');
        if !(a.is_ascii_digit()
            && b.is_ascii_digit()
            && (phased == 1 || sep == b'/')
            && gt.get(3).is_none_or(|&c| c == b'\t'))
        {
            v.indiv.truncate(start);
            return false;
        }
        v.indiv.push(((a - b'0' + 1) << 1) | phased);
        v.indiv.push(((b - b'0' + 1) << 1) | phased);
    }
    v.n_sample = nsamples as u32;
    true
}

/// Looks up the FORMAT keys, defining missing ones (HTSlib's `vcf_parse_format_dict2`).
/// Returns the number of fields.
fn dict2<H: HeaderAccess>(
    h: &mut H,
    s: &mut [u8],
    v: &mut Record,
    p: usize,
    fmt: &mut Vec<FmtAux>,
) -> Result<usize, ParseError> {
    let mut j = 0usize;
    let mut t = p;
    loop {
        let e = t + memchr::memchr2(b':', 0, &s[t..]).unwrap_or(s.len() - t);
        let last = s[e] == 0;
        if j >= MAX_N_FMT {
            v.errcode |= RecordError::LIMITS;
            log::error(
                "vcf_parse_format_dict2",
                format_args!(
                    "FORMAT column at {} lists more identifiers than htslib can handle",
                    record_location(h.get(), v)
                ),
            );
            return Err(ParseError::Invalid);
        }
        s[e] = 0;
        let name = &s[t..e];
        let id = match super::parse::cached_format_id(h.get(), j, name) {
            Some(id) if !h.get().id_slot_undefined(id, 2) => id,
            _ => {
                if name == b"." {
                    if !h.can_update() {
                        return Err(ParseError::HeaderUpdateNeeded);
                    }
                    log::error(
                        "vcf_parse_format_dict2",
                        format_args!(
                            "Invalid FORMAT tag name '.' at {}",
                            record_location(h.get(), v)
                        ),
                    );
                    v.errcode |= RecordError::TAG_INVALID;
                    return Err(ParseError::Invalid);
                }
                let name = name.to_vec();
                let mut dummy = b"##FORMAT=<ID=".to_vec();
                dummy.extend_from_slice(&name);
                dummy.extend_from_slice(b",Number=1,Type=String,Description=\"Dummy\">");
                define_dummy(h, v, "FORMAT", &name, ", assuming Type=String", &dummy, 2)?
            }
        };
        if fmt.len() <= j {
            fmt.push(FmtAux::default());
        }
        let f = &mut fmt[j];
        f.max_l = 0;
        f.max_m = 0;
        f.max_g = 0;
        f.key = id as i32;
        f.is_gt = &s[t..e] == b"GT";
        f.y = h.get().id_info(HeaderRecordType::Format, id).unwrap_or(15);
        v.n_fmt += 1;
        j += 1;
        if last {
            break;
        }
        t = e + 1;
    }
    Ok(j)
}

/// Measures the largest vector, string, and ploidy of each field over the samples (HTSlib's
/// `vcf_parse_format_max3`), turning sample-separating tabs into NULs.
fn max3<H: HeaderAccess>(
    h: &H,
    s: &mut [u8],
    end: usize,
    v: &mut Record,
    q: usize,
    fmt: &mut [FmtAux],
) -> Result<(), ParseError> {
    let header = h.get();
    let nsamples = header.sample_count() as u32;
    let keep = header.keep_samples.as_ref();
    let mut n_sample_ori: isize = -1;
    let mut r = q + 1;
    v.n_sample = 0;
    while r < end {
        if let Some(keep) = keep {
            n_sample_ori += 1;
            if !keep.contains(n_sample_ori as usize) {
                while s[r] != b'\t' && r < end {
                    r += 1;
                }
                if s[r] == b'\t' {
                    s[r] = 0;
                    r += 1;
                }
                continue;
            }
        }
        // One sample: its fields up to the tab (or the end of the line).
        let mut j = 0usize;
        let mut r_start = 0usize;
        let (mut m, mut g) = (1usize, 1usize);
        let mut seg_end = s.len();
        for (k, &c) in s[r..].iter().enumerate() {
            match c {
                b',' => m += 1,
                b'|' | b'/' => g += 1,
                b':' | b'\t' | 0 => {
                    let f = &mut fmt[j];
                    f.max_m = f.max_m.max(m);
                    f.max_l = f.max_l.max(k - r_start);
                    if f.is_gt {
                        f.max_g = f.max_g.max(g);
                    }
                    r_start = k;
                    m = 1;
                    g = 1;
                    if c != b':' {
                        seg_end = r + k;
                        break;
                    }
                    j += 1;
                    if j >= fmt.len() {
                        log::error(
                            "vcf_parse_format_max3",
                            format_args!(
                                "Incorrect number of FORMAT fields at {}:{}",
                                lossy(header.contig_name(v.rid).unwrap_or(b"(null)")),
                                v.pos + 1
                            ),
                        );
                        v.errcode |= RecordError::NCOLS;
                        return Err(ParseError::Invalid);
                    }
                }
                _ => {}
            }
        }
        r = seg_end;
        if s[r] == b'\t' {
            s[r] = 0;
        }
        v.n_sample += 1;
        if v.n_sample == nsamples {
            break;
        }
        r += 1;
    }
    Ok(())
}

/// Sizes the per-field buffers and drops duplicate or oversized fields (HTSlib's
/// `vcf_parse_format_alloc4`).
fn alloc4<H: HeaderAccess>(h: &H, v: &mut Record, fmt: &mut [FmtAux]) -> Result<(), ParseError> {
    let n_sample = v.n_sample as usize;
    let mut mem: u64 = 0;
    for f in fmt.iter_mut() {
        if f.max_m == 0 {
            f.max_m = 1;
        }
        let size = match f.htype() {
            3 => {
                if f.is_gt {
                    f.max_g << 2
                } else {
                    f.max_l
                }
            }
            1 | 2 => f.max_m << 2,
            ty => {
                log::error(
                    "vcf_parse_format_alloc4",
                    format_args!(
                        "The format type {ty} at {} is currently not supported",
                        record_location(h.get(), v)
                    ),
                );
                v.errcode |= RecordError::TAG_INVALID;
                return Err(ParseError::Invalid);
            }
        };
        mem = (mem + 7) & !7;
        if mem + n_sample as u64 * size as u64 > i32::MAX as u64 {
            if !log::fired(&EXCESSIVE_WARNED) {
                log::warn_once(
                    &EXCESSIVE_WARNED,
                    "vcf_parse_format_alloc4",
                    format_args!(
                        "Excessive memory required by FORMAT fields at {}",
                        record_location(h.get(), v)
                    ),
                );
            }
            v.errcode |= RecordError::LIMITS;
            f.size = None;
            continue;
        }
        mem += n_sample as u64 * size as u64;
        f.size = Some(size);
    }

    // check for duplicate tags
    for i in 1..fmt.len() {
        if fmt[i].size.is_none() {
            continue;
        }
        for j in 0..i {
            if fmt[j].size.is_none() || fmt[i].key != fmt[j].key {
                continue;
            }
            if !log::fired(&DUPLICATE_WARNED) {
                let header = h.get();
                log::warn_once(
                    &DUPLICATE_WARNED,
                    "vcf_parse_format_alloc4",
                    format_args!(
                        "Duplicate FORMAT tag {} at {}",
                        lossy(header.id_name(fmt[i].key as u32).unwrap_or(b"")),
                        record_location(header, v)
                    ),
                );
            }
            v.errcode |= RecordError::TAG_INVALID;
            fmt[i].size = None;
            break;
        }
    }

    for f in fmt.iter_mut() {
        let Some(size) = f.size else { continue };
        match f.htype() {
            3 if !f.is_gt => {
                f.chars.clear();
                f.chars.resize(n_sample * size, 0);
            }
            2 => {
                f.floats.clear();
                f.floats.resize(n_sample * (size >> 2), 0.0);
            }
            _ => {
                f.ints.clear();
                f.ints.resize(n_sample * (size >> 2), 0);
            }
        }
    }
    Ok(())
}

/// Fills the per-field buffers from the sample columns (HTSlib's `vcf_parse_format_fill5`).
fn fill5<H: HeaderAccess>(
    h: &H,
    s: &[u8],
    end: usize,
    v: &mut Record,
    q: usize,
    fmt: &mut [FmtAux],
) -> Result<(), ParseError> {
    let header = h.get();
    let keep = header.keep_samples.as_ref();
    let nsamples = header.sample_count();
    let ver = header.version();
    let mut n_sample_ori: isize = -1;
    let mut t = q + 1;
    let mut m = 0usize;
    while t < end {
        if let Some(keep) = keep {
            n_sample_ori += 1;
            if !keep.contains(n_sample_ori as usize) {
                while s[t] != 0 && t < end {
                    t += 1;
                }
                t += 1;
                continue;
            }
        }
        if m == nsamples {
            break;
        }
        let mut j = 0usize;
        while t < end {
            let z = &mut fmt[j];
            j += 1;
            match z.size {
                None => {
                    while s[t] != b':' && s[t] != 0 {
                        t += 1;
                    }
                }
                Some(size) => match z.htype() {
                    3 if z.is_gt => {
                        t = match fill_gt_fast(s, t, z, size, m) {
                            Some(t) => t,
                            None => fill_gt(header, s, t, v, z, size, m, ver)?,
                        }
                    }
                    3 => {
                        let x = &mut z.chars[size * m..size * (m + 1)];
                        let mut l = 0;
                        while s[t] != b':' && s[t] != 0 {
                            if l < x.len() {
                                x[l] = s[t];
                            }
                            l += 1;
                            t += 1;
                        }
                        if l < x.len() {
                            x[l..].fill(0);
                        }
                    }
                    1 => t = fill_int(header, s, t, v, z, size, m),
                    2 => t = fill_float(header, s, t, v, z, size, m),
                    ty => {
                        log::error(
                            "vcf_parse_format_fill5",
                            format_args!(
                                "Unknown FORMAT field type {ty} at {}",
                                record_location(header, v)
                            ),
                        );
                        v.errcode |= RecordError::TAG_INVALID;
                        return Err(ParseError::Invalid);
                    }
                },
            }
            match s[t] {
                0 => break,
                b':' => t += 1,
                c => {
                    log::error(
                        "vcf_parse_format_fill5",
                        format_args!(
                            "Invalid character {} in '{}' FORMAT field at {}",
                            log::strprint(&[c], b'\'', 8),
                            lossy(header.id_name(fmt[j - 1].key as u32).unwrap_or(b"")),
                            record_location(header, v)
                        ),
                    );
                    v.errcode |= RecordError::CHAR;
                    return Err(ParseError::Invalid);
                }
            }
        }

        // fill end-of-vector values
        for z in &mut fmt[j..] {
            let Some(size) = z.size else { continue };
            let nv = size >> 2;
            match z.htype() {
                3 if z.is_gt => {
                    let x = &mut z.ints[m * nv..(m + 1) * nv];
                    if size > 0 {
                        x[0] = INT32_MISSING;
                    }
                    if nv > 1 {
                        x[1..].fill(INT32_VECTOR_END);
                    }
                }
                3 => {
                    let x = &mut z.chars[m * size..(m + 1) * size];
                    if size > 0 {
                        x[0] = b'.';
                        x[1..].fill(0);
                    }
                }
                1 => {
                    let x = &mut z.ints[m * nv..(m + 1) * nv];
                    if let Some(first) = x.first_mut() {
                        *first = INT32_MISSING;
                    }
                    if nv > 1 {
                        x[1..].fill(INT32_VECTOR_END);
                    }
                }
                2 => {
                    let x = &mut z.floats[m * nv..(m + 1) * nv];
                    if let Some(first) = x.first_mut() {
                        *first = encoding::float_missing();
                    }
                    if nv > 1 {
                        x[1..].fill(encoding::float_vector_end());
                    }
                }
                _ => {}
            }
        }
        m += 1;
        t += 1;
    }
    Ok(())
}

/// Parses an unprefixed diploid genotype of single digits, such as `0|1`, the common case of
/// [`fill_gt`]; `None` for anything else.
#[inline(always)]
fn fill_gt_fast(s: &[u8], t: usize, z: &mut FmtAux, size: usize, m: usize) -> Option<usize> {
    let nv = size >> 2;
    let (c0, sep, c1, c3) = (s[t], at(s, t + 1), at(s, t + 2), at(s, t + 3));
    if nv >= 2
        && c0.is_ascii_digit()
        && (sep == b'|' || sep == b'/')
        && c1.is_ascii_digit()
        && !c3.is_ascii_digit()
        && c3 != b'|'
        && c3 != b'/'
    {
        let x = &mut z.ints[m * nv..(m + 1) * nv];
        let phased = i32::from(sep == b'|');
        x[0] = ((i32::from(c0 - b'0') + 1) << 1) | phased;
        x[1] = ((i32::from(c1 - b'0') + 1) << 1) | phased;
        x[2..].fill(INT32_VECTOR_END);
        return Some(t + 3);
    }
    None
}

/// Parses one genotype (the GT part of `vcf_parse_format_fill5`), returning the position after
/// it.
#[expect(
    clippy::too_many_arguments,
    reason = "a line-for-line port of HTSlib's GT parser"
)]
#[inline(never)]
fn fill_gt(
    header: &crate::variant::Header,
    s: &[u8],
    mut t: usize,
    v: &Record,
    z: &mut FmtAux,
    size: usize,
    m: usize,
    ver: i32,
) -> Result<usize, ParseError> {
    let nv = size >> 2;
    let x = &mut z.ints[m * nv..(m + 1) * nv];

    let mut is_phased: i32 = 0;
    let mut unreadable = false;
    let mut max: u32 = 0;
    let mut overflow = false;
    let mut ploidy = 0;
    let mut anyunphased = false;
    let mut phasingprfx = false;
    let mut unknown1 = false;

    // With prefixed phasing, the first allele's phasing is explicit.
    if ver >= VCF44 && (s[t] == b'|' || s[t] == b'/') {
        is_phased = i32::from(s[t] == b'|');
        t += 1;
        phasingprfx = true;
    }
    let mut l = 0usize;
    loop {
        ploidy += 1;
        if s[t] == b'.' {
            t += 1;
            if l < x.len() {
                x[l] = is_phased;
            }
            l += 1;
            if l == 1 {
                unknown1 = true;
            }
        } else {
            let tt = t;
            let val: u32;
            if s[t].is_ascii_digit() && !at(s, t + 1).is_ascii_digit() {
                val = u32::from(s[t] - b'0');
                t += 1;
            } else {
                let (n, e, of) = number::str2uint(s, t, 4 * 127 - 2);
                val = n as u32;
                overflow |= of;
                t = e;
                unreadable |= tt == t;
            }
            if max < val {
                max = val;
            }
            if l < x.len() {
                x[l] = ((val.wrapping_add(1)) << 1) as i32 | is_phased;
            }
            l += 1;
        }
        anyunphased |= ploidy != 1 && is_phased == 0;
        is_phased = i32::from(s[t] == b'|');
        if s[t] != b'|' && s[t] != b'/' {
            break;
        }
        t += 1;
    }
    if !phasingprfx && !x.is_empty() {
        if ploidy == 1 {
            if !unknown1 {
                x[0] |= 1;
            }
        } else if !anyunphased {
            x[0] |= 1;
        }
    }
    if overflow || max > (i32::MAX as u32 >> 1) - 1 {
        log::error(
            "vcf_parse_format_fill5",
            format_args!(
                "Couldn't read GT data: value too large at {}",
                record_location(header, v)
            ),
        );
        return Err(ParseError::Invalid);
    }
    if unreadable {
        log::error(
            "vcf_parse_format_fill5",
            format_args!(
                "Couldn't read GT data: value not a number or '.' at {}",
                record_location(header, v)
            ),
        );
        return Err(ParseError::Invalid);
    }
    if l == 0 && !x.is_empty() {
        x[0] = 0;
        l = 1;
    }
    if l < x.len() {
        x[l..].fill(INT32_VECTOR_END);
    }
    Ok(t)
}

fn fill_int(
    header: &crate::variant::Header,
    s: &[u8],
    mut t: usize,
    v: &Record,
    z: &mut FmtAux,
    size: usize,
    m: usize,
) -> usize {
    let nv = size >> 2;
    let key = z.key;
    let x = &mut z.ints[m * nv..(m + 1) * nv];
    let mut l = 0usize;
    loop {
        let value = if s[t] == b'.' {
            t += 1;
            INT32_MISSING
        } else {
            let (tmp, te, overflow) = number::str2int(s, t, 64);
            let value = if te == t
                || overflow
                || tmp < i64::from(MIN_INT32)
                || tmp > i64::from(MAX_INT32)
            {
                if !log::fired(&EXTREME_VAL_WARNED) {
                    log::warn_once(
                        &EXTREME_VAL_WARNED,
                        "vcf_parse_format_fill5",
                        format_args!(
                            "Extreme FORMAT/{} value encountered and set to missing at {}",
                            lossy(header.id_name(key as u32).unwrap_or(b"")),
                            record_location(header, v)
                        ),
                    );
                }
                INT32_MISSING
            } else {
                tmp as i32
            };
            t = te;
            value
        };
        if l < x.len() {
            x[l] = value;
        }
        l += 1;
        if s[t] != b',' {
            break;
        }
        t += 1;
    }
    if l == 0 && !x.is_empty() {
        x[0] = INT32_MISSING;
        l = 1;
    }
    if l < x.len() {
        x[l..].fill(INT32_VECTOR_END);
    }
    t
}

fn fill_float(
    header: &crate::variant::Header,
    s: &[u8],
    mut t: usize,
    v: &Record,
    z: &mut FmtAux,
    size: usize,
    m: usize,
) -> usize {
    let nv = size >> 2;
    let key = z.key;
    let x = &mut z.floats[m * nv..(m + 1) * nv];
    let mut l = 0usize;
    loop {
        let value = if s[t] == b'.' && !at(s, t + 1).is_ascii_digit() {
            t += 1;
            encoding::float_missing()
        } else {
            let (d, te, overflow) = number::str2dbl(s, t);
            if (te == t || overflow) && !log::fired(&EXTREME_VAL_WARNED) {
                let chrom = header
                    .contig_name(v.rid)
                    .map_or_else(|| "(null)".into(), lossy);
                log::warn_once(
                    &EXTREME_VAL_WARNED,
                    "vcf_parse_format_fill5",
                    format_args!(
                        "Extreme FORMAT/{} value encountered at {chrom}:{}",
                        lossy(header.id_name(key as u32).unwrap_or(b"")),
                        v.pos + 1
                    ),
                );
            }
            t = te;
            d as f32
        };
        if l < x.len() {
            x[l] = value;
        }
        l += 1;
        if s[t] != b',' {
            break;
        }
        t += 1;
    }
    if l == 0 && !x.is_empty() {
        x[0] = encoding::float_missing();
        l = 1;
    }
    if l < x.len() {
        x[l..].fill(encoding::float_vector_end());
    }
    t
}

/// Encodes the per-field buffers into the sample block (HTSlib's `vcf_parse_format_gt6`).
fn gt6(v: &mut Record, fmt: &[FmtAux]) {
    let n_sample = v.n_sample as usize;
    let str = &mut v.indiv;
    let mut dropped = 0u32;
    if n_sample > 0 {
        for z in fmt {
            let Some(size) = z.size else {
                // HTSlib drops these fields only when there are samples.
                dropped += 1;
                continue;
            };
            enc_int1(str, z.key);
            match z.htype() {
                3 if !z.is_gt => {
                    enc_size(str, size, BcfType::Char as u8);
                    str.extend_from_slice(&z.chars[..size * n_sample]);
                }
                2 => {
                    enc_size(str, size >> 2, BcfType::Float as u8);
                    encoding::push_floats(str, &z.floats[..(size >> 2) * n_sample]);
                }
                _ => enc_vint(str, &z.ints[..(size >> 2) * n_sample], Some(size >> 2)),
            }
        }
    }
    v.n_fmt -= dropped;
}

/// Checks the sample count (HTSlib's `vcf_parse_format_check7`).
fn check7<H: HeaderAccess>(h: &H, v: &mut Record) -> Result<(), ParseError> {
    let header = h.get();
    if v.n_sample as usize != header.sample_count() {
        log::error(
            "vcf_parse_format_check7",
            format_args!(
                "Number of columns at {} does not match the number of samples ({} vs {})",
                record_location(header, v),
                v.n_sample,
                header.sample_count()
            ),
        );
        v.errcode |= RecordError::NCOLS;
        return Err(ParseError::Invalid);
    }
    if v.indiv.len() > u32::MAX as usize {
        log::error(
            "vcf_parse_format_check7",
            format_args!("The FORMAT at {} is too long", record_location(header, v)),
        );
        v.errcode |= RecordError::LIMITS;
        v.n_fmt = 0;
        return Err(ParseError::Invalid);
    }
    Ok(())
}
