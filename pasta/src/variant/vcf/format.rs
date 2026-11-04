//! Formatting records as VCF text (HTSlib's `vcf_format`, `bcf_fmt_array`,
//! `bcf_format_gt_v2`).

use crate::variant::{
    Header, Record, VCF44,
    encoding::{
        self, FLOAT_MISSING, FLOAT_VECTOR_END, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING,
        INT16_VECTOR_END, INT32_MISSING, INT32_VECTOR_END, STR_MISSING, type_shift,
    },
    log,
    number::{push_double, push_int},
    record::{FieldData, Unpack},
};

/// An error formatting a record: it refers to an ID the header does not define.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum FormatError {
    /// The contig ID is not in the header.
    InvalidContig(i32),
    /// A FILTER ID is not in the header.
    InvalidFilter(i32),
    /// An INFO ID is not in the header.
    InvalidInfo(i32),
    /// A FORMAT ID is not in the header.
    InvalidFormat(i32),
    /// A value has a type that cannot be formatted.
    UnexpectedType(u8),
}

impl std::fmt::Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidContig(id) => {
                write!(f, "invalid BCF, CONTIG id={id} not present in the header")
            }
            Self::InvalidFilter(id) => {
                write!(f, "invalid BCF, FILTER id={id} not present in the header")
            }
            Self::InvalidInfo(id) => {
                write!(f, "invalid BCF, INFO id={id} not present in the header")
            }
            Self::InvalidFormat(id) => {
                write!(f, "invalid BCF, FORMAT id={id} not present in the header")
            }
            Self::UnexpectedType(t) => write!(f, "unexpected BCF type {t}"),
        }
    }
}

impl std::error::Error for FormatError {}

/// Appends one value of a FORMAT field with one value per sample (HTSlib's `bcf_fmt_array1`).
#[inline]
fn fmt_array1(s: &mut Vec<u8>, ty: u8, p: &[u8]) -> Result<(), FormatError> {
    match ty {
        7 => s.push(if p[0] == STR_MISSING { b'.' } else { p[0] }),
        1 => {
            let v = p[0] as i8;
            if v != INT8_VECTOR_END {
                if v == INT8_MISSING {
                    s.push(b'.');
                } else {
                    push_int(s, i64::from(v));
                }
            }
        }
        2 => {
            let v = encoding::le_i16(p);
            if v != INT16_VECTOR_END {
                if v == INT16_MISSING {
                    s.push(b'.');
                } else {
                    push_int(s, i64::from(v));
                }
            }
        }
        3 => {
            let v = encoding::le_i32(p);
            if v != INT32_VECTOR_END {
                if v == INT32_MISSING {
                    s.push(b'.');
                } else {
                    push_int(s, i64::from(v));
                }
            }
        }
        5 => {
            let v = encoding::le_u32(p);
            if v != FLOAT_VECTOR_END {
                if v == FLOAT_MISSING {
                    s.push(b'.');
                } else {
                    push_double(s, f64::from(f32::from_bits(v)));
                }
            }
        }
        _ => {
            log::error("bcf_fmt_array1", format_args!("Unexpected type {ty}"));
            return Err(FormatError::UnexpectedType(ty));
        }
    }
    Ok(())
}

/// Appends `n` values of raw type `ty` from `data`, comma-separated, stopping at a vector end
/// (HTSlib's `bcf_fmt_array`). Zero values print as `.`.
pub(crate) fn fmt_array(s: &mut Vec<u8>, n: usize, ty: u8, data: &[u8]) -> bool {
    if n == 0 {
        s.push(b'.');
        return true;
    }
    match ty {
        7 => {
            let d = &data[..n.min(data.len())];
            let end = if n >= 8 {
                memchr::memchr(0, d).unwrap_or(d.len())
            } else {
                d.iter().position(|&c| c == 0).unwrap_or(d.len())
            };
            s.extend_from_slice(&d[..end]);
        }
        1 => {
            for (j, &b) in data[..n].iter().enumerate() {
                let v = b as i8;
                if v == INT8_VECTOR_END {
                    break;
                }
                if j > 0 {
                    s.push(b',');
                }
                if v == INT8_MISSING {
                    s.push(b'.');
                } else {
                    push_int(s, i64::from(v));
                }
            }
        }
        2 => {
            for j in 0..n {
                let v = encoding::le_i16(&data[j * 2..]);
                if v == INT16_VECTOR_END {
                    break;
                }
                if j > 0 {
                    s.push(b',');
                }
                if v == INT16_MISSING {
                    s.push(b'.');
                } else {
                    push_int(s, i64::from(v));
                }
            }
        }
        3 => {
            for j in 0..n {
                let v = encoding::le_i32(&data[j * 4..]);
                if v == INT32_VECTOR_END {
                    break;
                }
                if j > 0 {
                    s.push(b',');
                }
                if v == INT32_MISSING {
                    s.push(b'.');
                } else {
                    push_int(s, i64::from(v));
                }
            }
        }
        5 => {
            for j in 0..n {
                let v = encoding::le_u32(&data[j * 4..]);
                if v == FLOAT_VECTOR_END {
                    break;
                }
                if j > 0 {
                    s.push(b',');
                }
                if v == FLOAT_MISSING {
                    s.push(b'.');
                } else {
                    push_double(s, f64::from(f32::from_bits(v)));
                }
            }
        }
        _ => {
            log::error("bcf_fmt_array", format_args!("Unexpected type {ty}"));
            return false;
        }
    }
    true
}

/// Appends one sample's genotype (HTSlib's `bcf_format_gt_v2`). With a VCF 4.4+ header, the
/// first allele's phasing is written as a prefix when it cannot be inferred.
pub(crate) fn format_gt(
    s: &mut Vec<u8>,
    version: i32,
    ty: u8,
    n: usize,
    p: &[u8],
) -> Result<(), FormatError> {
    // Fast path: a diploid 8-bit genotype of alleles 0-9 or missing, before VCF 4.4.
    if ty == 1 && n >= 2 && version < VCF44 {
        let (v0, v1) = (p[0] as i8, p[1] as i8);
        if (0..22).contains(&v0)
            && (0..22).contains(&v1)
            && (n == 2 || p[2] as i8 == INT8_VECTOR_END)
        {
            let digit = |v: i8| {
                if v >> 1 == 0 {
                    b'.'
                } else {
                    b'0' + ((v >> 1) - 1) as u8
                }
            };
            s.extend_from_slice(&[digit(v0), if v1 & 1 != 0 { b'|' } else { b'/' }, digit(v1)]);
            return Ok(());
        }
    }
    let pos = s.len();
    let mut anyunphased = false;
    let mut val0: i32 = 0;
    let mut emit = |i: usize, val: i32, s: &mut Vec<u8>| {
        if i == 0 {
            val0 = val;
        } else {
            s.push(if val & 1 != 0 { b'|' } else { b'/' });
            anyunphased |= val & 1 == 0;
        }
        let allele = (val >> 1) - 1;
        if val >> 1 == 0 {
            s.push(b'.');
        } else if (0..10).contains(&allele) {
            s.push(b'0' + allele as u8);
        } else {
            push_int(s, i64::from(allele));
        }
    };
    let ploidy = match ty {
        1 => {
            let mut i = 0;
            while i < n {
                let v = p[i] as i8;
                if v == INT8_VECTOR_END {
                    break;
                }
                emit(i, i32::from(v), s);
                i += 1;
            }
            if i == 0 {
                s.push(b'.');
            }
            i
        }
        2 => {
            let mut i = 0;
            while i < n {
                let v = encoding::le_i16(&p[i * 2..]);
                if v == INT16_VECTOR_END {
                    break;
                }
                emit(i, i32::from(v), s);
                i += 1;
            }
            if i == 0 {
                s.push(b'.');
            }
            i
        }
        3 => {
            let mut i = 0;
            while i < n {
                let v = encoding::le_i32(&p[i * 4..]);
                if v == INT32_VECTOR_END {
                    break;
                }
                emit(i, v, s);
                i += 1;
            }
            if i == 0 {
                s.push(b'.');
            }
            i
        }
        0 => {
            s.push(b'.');
            1
        }
        _ => {
            log::error("bcf_format_gt_v2", format_args!("Unexpected type {ty}"));
            return Err(FormatError::UnexpectedType(ty));
        }
    };
    if version >= VCF44 {
        let insert = if val0 & 1 != 0 {
            if (ploidy > 1 && anyunphased) || (ploidy <= 1 && val0 >> 1 == 0) {
                Some(b'|')
            } else {
                None
            }
        } else if (ploidy <= 1 && val0 != 0) || (ploidy > 1 && !anyunphased) {
            Some(b'/')
        } else {
            None
        };
        if let Some(c) = insert {
            s.insert(pos, c);
        }
    }
    Ok(())
}

/// Appends `n` values of raw BCF type `ty` as VCF text: comma-separated, `.` for missing,
/// stopping at a vector end, and `.` for no values (HTSlib's `bcf_fmt_array`). Returns `false`
/// for a type that cannot be formatted.
pub fn format_array(out: &mut Vec<u8>, n: usize, ty: u8, data: &[u8]) -> bool {
    fmt_array(out, n, ty, data)
}

/// Appends one sample's genotype as VCF text (HTSlib's `bcf_format_gt_v2`): `data` holds the
/// sample's `n` values of raw type `ty`; with a VCF 4.4+ header (`header.cached_version()`),
/// the first allele's phasing is written when it cannot be inferred.
///
/// # Errors
///
/// Returns an error for a non-integer type.
pub fn format_genotype(
    out: &mut Vec<u8>,
    version: i32,
    ty: u8,
    n: usize,
    data: &[u8],
) -> Result<(), FormatError> {
    format_gt(out, version, ty, n, data)
}

/// Appends a record as a VCF line, with its newline (HTSlib's `vcf_format`).
///
/// The record need not be unpacked: undecoded parts are formatted straight from the encoded
/// blocks. Decoded parts (and edits to them) are used when present.
///
/// # Errors
///
/// Returns an error when the record refers to a contig or tag missing from the header (also
/// reported on stderr as HTSlib does).
pub fn format_record(header: &Header, v: &Record, s: &mut Vec<u8>) -> Result<(), FormatError> {
    let Some(chrom) = header.contig_name(v.rid) else {
        log::error(
            "vcf_format",
            format_args!("Invalid BCF, CONTIG id={} not present in the header", v.rid),
        );
        return Err(FormatError::InvalidContig(v.rid));
    };
    let shared = &v.shared;
    let str_unpacked = v.unpacked & Unpack::STR.bits() != 0;
    let flt_unpacked = v.unpacked & Unpack::FLT.bits() != 0;
    s.reserve(256 + shared.len() * 2 + v.indiv.len() * 3);
    s.extend_from_slice(chrom);
    s.push(b'\t');
    push_int(s, v.pos + 1);
    s.push(b'\t');

    // ID, REF, ALT
    let mut p; // offset of FILTER in the packed block
    // A record built through the API (HTSlib prints its `d` fields, with a NULL ID as `.`).
    let built = !str_unpacked && shared.is_empty() && v.n_allele > 0;
    if str_unpacked || built {
        if built && v.d.id.is_empty() {
            s.push(b'.');
        } else {
            s.extend_from_slice(&v.d.id);
        }
        s.push(b'\t');
        if v.n_allele > 0 {
            s.extend_from_slice(v.allele(0));
        } else {
            s.push(b'.');
        }
        s.push(b'\t');
        if v.n_allele > 1 {
            for i in 1..v.n_allele as usize {
                if i > 1 {
                    s.push(b',');
                }
                s.extend_from_slice(v.allele(i));
            }
        } else {
            s.push(b'.');
        }
        p = if str_unpacked {
            v.unpack_size[0] + v.unpack_size[1]
        } else {
            0
        };
    } else if shared.is_empty() {
        s.extend_from_slice(b".\t.\t.");
        p = 0;
    } else {
        p = sized_text(s, shared, 0);
        s.push(b'\t');
        if v.n_allele > 0 {
            p = sized_text(s, shared, p);
        } else {
            s.push(b'.');
        }
        s.push(b'\t');
        if v.n_allele > 1 {
            for i in 1..v.n_allele {
                if i > 1 {
                    s.push(b',');
                }
                p = sized_text(s, shared, p);
            }
        } else {
            s.push(b'.');
        }
    }

    // QUAL
    s.push(b'\t');
    if encoding::float_is_missing(v.qual) {
        s.push(b'.');
    } else {
        push_double(s, f64::from(v.qual));
    }

    // FILTER
    s.push(b'\t');
    let max_id = header.id_count() as i64;
    let write_filter = |s: &mut Vec<u8>, i: usize, idx: i32| -> Result<(), FormatError> {
        let name = (idx >= 0 && i64::from(idx) < max_id)
            .then(|| header.id_name(idx as u32))
            .flatten();
        let Some(name) = name else {
            log::error(
                "vcf_format",
                format_args!(
                    "Invalid BCF, the FILTER tag id={idx} at {}:{} not present in the header",
                    String::from_utf8_lossy(chrom),
                    v.pos + 1
                ),
            );
            return Err(FormatError::InvalidFilter(idx));
        };
        if i > 0 {
            s.push(b';');
        }
        s.extend_from_slice(name);
        Ok(())
    };
    if flt_unpacked {
        if v.d.flt.is_empty() {
            s.push(b'.');
        }
        for (i, &idx) in v.d.flt.iter().enumerate() {
            write_filter(s, i, idx)?;
        }
        p = v.unpack_size[0] + v.unpack_size[1] + v.unpack_size[2];
    } else if p < shared.len() {
        let (n, ty, used) = encoding::dec_size(&shared[p..]);
        p += used;
        let n = if shared[p - used] >> 4 != 0 {
            n.max(0) as usize
        } else {
            0
        };
        if n == 0 {
            s.push(b'.');
        }
        for i in 0..n {
            let (idx, used) = encoding::dec_int1(&shared[p..], ty);
            p += used;
            write_filter(s, i, idx as i32)?;
        }
        if n == 0 && shared[p - 1] >> 4 == 0 {
            // A zero-length FILTER descriptor; nothing more to skip.
        }
    } else {
        s.push(b'.');
    }

    // INFO
    s.push(b'\t');
    format_info(header, v, s, p, chrom)?;

    // FORMAT and samples
    if v.n_sample > 0 {
        if v.n_fmt > 0 {
            format_samples(header, v, s, chrom)?;
        } else {
            for _ in 0..=v.n_sample {
                s.extend_from_slice(b"\t.");
            }
        }
    }
    s.push(b'\n');
    Ok(())
}

/// Appends a sized character array as `bcf_fmt_sized_array` does, returning the next offset.
#[inline]
fn sized_text(s: &mut Vec<u8>, data: &[u8], p: usize) -> usize {
    let (n, ty, used) = encoding::dec_size(&data[p..]);
    let start = p + used;
    let n = n.max(0) as usize;
    fmt_array(s, n, ty, &data[start..]);
    start + (n << type_shift(ty))
}

fn info_error(header: &Header, key: i32, chrom: &[u8], pos: i64) -> FormatError {
    let why = if key < 0 {
        "negative"
    } else if key as usize >= header.id_count() {
        "too large"
    } else {
        "not present in the header"
    };
    log::error(
        "vcf_format",
        format_args!(
            "Invalid BCF, the INFO tag id={key} is {why} at {}:{}",
            String::from_utf8_lossy(chrom),
            pos + 1
        ),
    );
    FormatError::InvalidInfo(key)
}

fn format_info(
    header: &Header,
    v: &Record,
    s: &mut Vec<u8>,
    mut p: usize,
    chrom: &[u8],
) -> Result<(), FormatError> {
    if v.n_info == 0 {
        s.push(b'.');
        return Ok(());
    }
    let shared = &v.shared;
    let info_packed = v.unpacked & Unpack::INFO.bits() == 0 && !shared.is_empty();
    let mut first = true;
    if info_packed {
        for _ in 0..v.n_info {
            let (key, used) = encoding::dec_typed_int1(&shared[p..]);
            p += used;
            let (len, ty, used) = encoding::dec_size(&shared[p..]);
            p += used;
            let key = key as i32;
            let Some(name) = (key >= 0).then(|| header.id_name(key as u32)).flatten() else {
                return Err(info_error(header, key, chrom, v.pos));
            };
            if !first {
                s.push(b';');
            }
            first = false;
            s.extend_from_slice(name);
            let len = len.max(0) as usize;
            let data = &shared[p..];
            p += len << type_shift(ty);
            if len == 0 {
                continue;
            }
            s.push(b'=');
            if !fmt_array(s, len, ty, data) {
                return Err(FormatError::UnexpectedType(ty));
            }
        }
    } else {
        for z in &v.d.info[..v.n_info as usize] {
            if z.data == FieldData::Removed {
                continue;
            }
            let Some(name) = (z.key >= 0).then(|| header.id_name(z.key as u32)).flatten() else {
                return Err(info_error(header, z.key, chrom, v.pos));
            };
            if !first {
                s.push(b';');
            }
            first = false;
            s.extend_from_slice(name);
            if z.len <= 0 {
                continue;
            }
            s.push(b'=');
            if z.len != 1 {
                if !fmt_array(s, z.len as usize, z.ty, v.info_data(z)) {
                    return Err(FormatError::UnexpectedType(z.ty));
                }
                continue;
            }
            match z.ty {
                5 => {
                    let f = z.v1 as u32;
                    if f == FLOAT_MISSING {
                        s.push(b'.');
                    } else {
                        push_double(s, f64::from(f32::from_bits(f)));
                    }
                }
                7 => s.push(z.v1 as u8),
                0..=3 => {
                    let missing = [0i64, -128, -32768, i64::from(i32::MIN)][z.ty as usize];
                    if z.v1 == missing {
                        s.push(b'.');
                    } else {
                        push_int(s, i64::from(z.v1 as i32));
                    }
                }
                4 => {
                    if z.v1 == encoding::INT64_MISSING {
                        s.push(b'.');
                    } else {
                        push_int(s, z.v1);
                    }
                }
                ty => {
                    log::error(
                        "vcf_format",
                        format_args!(
                            "Unexpected type {ty} at {}:{}",
                            String::from_utf8_lossy(chrom),
                            v.pos + 1
                        ),
                    );
                    return Err(FormatError::UnexpectedType(ty));
                }
            }
        }
    }
    if first {
        s.push(b'.');
    }
    Ok(())
}

/// A FORMAT field located for formatting.
#[derive(Clone, Copy, Debug, Default)]
struct FmtDesc {
    id: i32,
    n: usize,
    size: usize,
    ty: u8,
    removed: bool,
    /// Offset of the values in the sample block (packed records).
    start: usize,
    /// Index of the decoded field (unpacked records).
    idx: usize,
}

thread_local! {
    static FMT_DESCS: std::cell::RefCell<Vec<FmtDesc>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn format_samples(
    header: &Header,
    v: &Record,
    s: &mut Vec<u8>,
    chrom: &[u8],
) -> Result<(), FormatError> {
    FMT_DESCS.with(|descs| {
        let mut descs = descs.borrow_mut();
        format_samples_with(header, v, s, chrom, &mut descs)
    })
}

fn format_samples_with(
    header: &Header,
    v: &Record,
    s: &mut Vec<u8>,
    chrom: &[u8],
    descs: &mut Vec<FmtDesc>,
) -> Result<(), FormatError> {
    let n_sample = v.n_sample as usize;
    let fmt_packed = v.unpacked & Unpack::FMT.bits() == 0;
    descs.clear();
    if fmt_packed {
        let indiv = &v.indiv;
        let mut p = 0;
        for _ in 0..v.n_fmt {
            let (id, used) = encoding::dec_typed_int1(&indiv[p..]);
            p += used;
            let (n, ty, used) = encoding::dec_size(&indiv[p..]);
            p += used;
            let n = n.max(0) as usize;
            let size = n << type_shift(ty);
            descs.push(FmtDesc {
                id: id as i32,
                n,
                size,
                ty,
                removed: false,
                start: p,
                idx: 0,
            });
            p += n_sample * size;
        }
    } else {
        for (idx, f) in v.d.fmt[..(v.n_fmt as usize).min(v.d.fmt.len())]
            .iter()
            .enumerate()
        {
            descs.push(FmtDesc {
                id: f.id,
                n: f.n(),
                size: f.size(),
                ty: f.ty,
                removed: f.is_removed(),
                start: 0,
                idx,
            });
        }
    }
    let data_of = |d: &FmtDesc| -> &[u8] {
        if fmt_packed {
            &v.indiv[d.start..d.start + n_sample * d.size]
        } else {
            v.format_data(&v.d.fmt[d.idx])
        }
    };

    // KEYS
    let mut gt_i: Option<usize> = None;
    let mut first = true;
    for (i, z) in descs.iter().enumerate() {
        if z.removed {
            continue;
        }
        s.push(if first { b'\t' } else { b':' });
        first = false;
        let Some(name) = (z.id >= 0).then(|| header.id_name(z.id as u32)).flatten() else {
            log::error(
                "vcf_format",
                format_args!(
                    "Invalid BCF, the FORMAT tag id={} at {}:{} not present in the header",
                    z.id,
                    String::from_utf8_lossy(chrom),
                    v.pos + 1
                ),
            );
            return Err(FormatError::InvalidFormat(z.id));
        };
        s.extend_from_slice(name);
        if name == b"GT" {
            gt_i = Some(i);
        }
    }
    if first {
        s.extend_from_slice(b"\t.");
    }

    // Single GT field: the common case, formatted without the per-field loop.
    let version = header.cached_version();
    if descs.len() == 1 && gt_i == Some(0) {
        let f = descs[0];
        let data = data_of(&f);
        for j in 0..n_sample {
            s.push(b'\t');
            format_gt(s, version, f.ty, f.n, &data[j * f.size..]).inspect_err(|_| gt_error(0))?;
        }
        return Ok(());
    }

    // VALUES per sample
    for j in 0..n_sample {
        s.push(b'\t');
        let mut first = true;
        for (i, f) in descs.iter().enumerate() {
            if f.removed {
                continue;
            }
            if !first {
                s.push(b':');
            }
            first = false;
            let p = &data_of(f)[j * f.size..];
            if gt_i == Some(i) {
                format_gt(s, version, f.ty, f.n, p).inspect_err(|_| gt_error(i))?;
            } else if f.n == 1 {
                fmt_array1(s, f.ty, p)?;
            } else if !fmt_array(s, f.n, f.ty, p) {
                return Err(FormatError::UnexpectedType(f.ty));
            }
        }
        if first {
            s.push(b'.');
        }
    }
    Ok(())
}

fn gt_error(i: usize) {
    log::error(
        "vcf_format",
        format_args!("Failed to format GT value for sample {i}, returned -2"),
    );
}
