//! BCF records and headers on the wire (HTSlib's `bcf_read1_core`, `bcf_record_check`,
//! `bcf_write`, `bcf_hdr_read`, `bcf_hdr_write`).

use std::{
    io::{self, BufRead, Read},
    sync::atomic::AtomicBool,
};

use crate::variant::{
    Header, HeaderTextFlavor, Record, VCF44,
    encoding::{self, type_shift},
    log,
    record::{IS_64BIT, RecordError, rlen::get_rlen},
};

/// The BCF 2.2 magic bytes.
pub const MAGIC: &[u8; 5] = b"BCF\x02\x02";

/// Reads the BCF magic and header text and parses it (HTSlib's `bcf_hdr_read` for BCF).
///
/// # Errors
///
/// Returns an error for a wrong magic number, a short read, or an unparsable header.
pub fn read_header<R: Read>(reader: &mut R) -> io::Result<Header> {
    let mut magic = [0u8; 5];
    if read_full(reader, &mut magic)? != 5 {
        log::error(
            "bcf_hdr_read",
            format_args!("Failed to read the header (reading BCF in text mode?)"),
        );
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "failed to read the BCF header",
        ));
    }
    if &magic != MAGIC {
        if &magic[..3] == b"BCF" {
            log::error(
                "bcf_hdr_read",
                format_args!("Invalid BCF2 magic string: only BCFv2.2 is supported"),
            );
        } else {
            log::error("bcf_hdr_read", format_args!("Invalid BCF2 magic string"));
        }
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BCF magic string",
        ));
    }
    let fail = || {
        log::error("bcf_hdr_read", format_args!("Failed to read BCF header"));
        io::Error::new(io::ErrorKind::InvalidData, "failed to read the BCF header")
    };
    let mut len = [0u8; 4];
    if read_full(reader, &mut len)? != 4 {
        return Err(fail());
    }
    let hlen = u32::from_le_bytes(len) as usize;
    let mut text = vec![0u8; hlen];
    if read_full(reader, &mut text)? != hlen {
        return Err(fail());
    }
    // The text is NUL-terminated; HTSlib parses it as a C string.
    if let Some(nul) = memchr::memchr(0, &text) {
        text.truncate(nul);
    }
    Header::parse(&text).map_err(|_| fail())
}

/// Returns the BCF header block: magic, length, and the text with IDX fields and a NUL
/// (HTSlib's `bcf_hdr_write` for BCF).
#[must_use]
pub fn encode_header(header: &Header) -> Vec<u8> {
    let mut text = Vec::new();
    header.format_into(&mut text, HeaderTextFlavor::Bcf);
    text.push(0);
    let mut out = Vec::with_capacity(text.len() + 9);
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(&text);
    out
}

fn read_full<R: Read>(reader: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match reader.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// The outcome of reading one BCF record.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ReadOutcome {
    /// A record was read.
    Record,
    /// The stream ended cleanly.
    Eof,
}

/// Replaces `out` with the next `n` bytes of `reader`, copied straight from its buffer (no
/// zero-filling, one copy). Returns `false` when the input ends first.
fn read_exact_into<R: BufRead>(reader: &mut R, out: &mut Vec<u8>, n: usize) -> io::Result<bool> {
    out.clear();
    out.reserve(n);
    while out.len() < n {
        let buf = match reader.fill_buf() {
            Ok(buf) => buf,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        };
        if buf.is_empty() {
            return Ok(false);
        }
        let k = buf.len().min(n - out.len());
        out.extend_from_slice(&buf[..k]);
        reader.consume(k);
    }
    Ok(true)
}

/// Reads the 32-byte fixed part of a record, returning how many bytes were read (fewer at
/// the end of the input).
#[inline]
fn read_fixed<R: BufRead>(reader: &mut R, x: &mut [u8; 32]) -> io::Result<usize> {
    // The common case: the whole fixed part is in the current block. A read error is the
    // reader's to report now; trying again could report a different one.
    let buf = loop {
        match reader.fill_buf() {
            Ok(buf) => break buf,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    };
    if let Some(head) = buf.first_chunk::<32>() {
        *x = *head;
        reader.consume(32);
        return Ok(32);
    }
    read_full(reader, x)
}

/// Reads one record's blocks (HTSlib's `bcf_read1_core`), without checking them.
///
/// # Errors
///
/// Returns an error for a truncated or malformed record.
pub fn read_record_raw<R: BufRead>(reader: &mut R, v: &mut Record) -> io::Result<ReadOutcome> {
    let mut x = [0u8; 32];
    let n = read_fixed(reader, &mut x)?;
    if n == 0 {
        return Ok(ReadOutcome::Eof);
    }
    if n != 32 {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated BCF record",
        ));
    }
    v.clear();
    let shared_len = encoding::le_u32(&x[0..]);
    if shared_len < 24 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "malformed BCF record",
        ));
    }
    let shared_len = (shared_len - 24) as usize;
    let indiv_len = encoding::le_u32(&x[4..]) as usize;
    v.rid = encoding::le_i32(&x[8..]);
    let pos = encoding::le_u32(&x[12..]);
    v.pos = if pos == u32::MAX { -1 } else { i64::from(pos) };
    v.rlen = i64::from(encoding::le_i32(&x[16..]));
    v.qual = f32::from_bits(encoding::le_u32(&x[20..]));
    v.n_info = u32::from(u16::from_le_bytes([x[24], x[25]]));
    v.n_allele = u32::from(u16::from_le_bytes([x[26], x[27]]));
    v.n_sample = encoding::le_u32(&x[28..]) & 0xff_ffff;
    v.n_fmt = u32::from(x[31]);
    // Silent fix of broken BCFs produced by bcf_subset prior to HTSlib commit bd6ed8b4.
    if (indiv_len == 0 || v.n_sample == 0) && v.n_fmt != 0 {
        v.n_fmt = 0;
    }
    if !read_exact_into(reader, &mut v.shared, shared_len)?
        || !read_exact_into(reader, &mut v.indiv, indiv_len)?
    {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated BCF record",
        ));
    }
    Ok(ReadOutcome::Record)
}

/// The size of the length that starts a [`pack_record`] entry.
pub const PACKED_HEADER_LEN: usize = 4;

/// Appends a record in a lossless in-memory form (after bcftools sort's `pack_bcf_data`): the
/// length of the rest, the contig, the 64-bit position and QUAL, the other fixed fields and the
/// length of the shared block as varints, then the shared and per-sample blocks. Unlike the
/// BCF encoding it keeps positions and values that need 64 bits. As in bcftools, the 64-bit
/// flag is not kept, so an unpacked record is refused as BCF only for what BCF cannot store.
/// Edits are synced first.
pub fn pack_record(rec: &mut Record, out: &mut Vec<u8>) {
    rec.sync();
    let start = out.len();
    out.reserve(PACKED_HEADER_LEN + 16 + 6 * 10 + rec.shared.len() + rec.indiv.len());
    out.extend_from_slice(&[0; PACKED_HEADER_LEN]);
    out.extend_from_slice(&rec.rid.to_le_bytes());
    out.extend_from_slice(&rec.pos.to_le_bytes());
    out.extend_from_slice(&rec.qual.to_bits().to_le_bytes());
    push_varint(out, ((rec.rlen << 1) ^ (rec.rlen >> 63)).cast_unsigned());
    for n in [rec.n_info, rec.n_allele, rec.n_fmt, rec.n_sample] {
        push_varint(out, u64::from(n));
    }
    push_varint(out, rec.shared.len() as u64);
    out.extend_from_slice(&rec.shared);
    out.extend_from_slice(&rec.indiv);
    let len = (out.len() - start - PACKED_HEADER_LEN) as u32;
    out[start..start + PACKED_HEADER_LEN].copy_from_slice(&len.to_le_bytes());
}

/// Appends `v` as an unsigned LEB128 varint.
fn push_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push(v as u8 | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// Reads an unsigned LEB128 varint at `*p`, moving `*p` past it; `None` when `data` ends first
/// or the value does not fit 64 bits.
fn read_varint(data: &[u8], p: &mut usize) -> Option<u64> {
    let mut v = 0u64;
    for shift in (0..64).step_by(7) {
        let b = *data.get(*p)?;
        *p += 1;
        v |= u64::from(b & 0x7f) << shift;
        if b & 0x80 == 0 {
            return Some(v);
        }
    }
    None
}

/// Returns the total length of the [`pack_record`] entry starting at `data`, or `None` when
/// `data` is shorter than its length.
#[must_use]
#[inline]
pub fn packed_len(data: &[u8]) -> Option<usize> {
    if data.len() < PACKED_HEADER_LEN {
        return None;
    }
    Some(PACKED_HEADER_LEN + encoding::le_u32(data) as usize)
}

/// Reads a [`pack_record`] entry into `rec` (bcftools sort's `read_packed_bcf`), returning the
/// number of bytes it took.
///
/// # Errors
///
/// Returns an error when `data` is shorter than the entry or the entry is malformed.
pub fn unpack_record(data: &[u8], rec: &mut Record) -> io::Result<usize> {
    let truncated = || io::Error::new(io::ErrorKind::UnexpectedEof, "truncated packed record");
    let len = packed_len(data).ok_or_else(truncated)?;
    let body = data.get(PACKED_HEADER_LEN..len).ok_or_else(truncated)?;
    let malformed = || io::Error::new(io::ErrorKind::InvalidData, "malformed packed record");
    if body.len() < 16 {
        return Err(malformed());
    }
    rec.clear();
    rec.rid = encoding::le_i32(body);
    rec.pos = encoding::le_i64(&body[4..]);
    rec.qual = f32::from_bits(encoding::le_u32(&body[12..]));
    let mut p = 16;
    let mut next = || read_varint(body, &mut p).ok_or_else(malformed);
    let rlen = next()?;
    rec.rlen = (rlen >> 1).cast_signed() ^ -((rlen & 1).cast_signed());
    let mut count = || u32::try_from(next()?).map_err(|_| malformed());
    rec.n_info = count()?;
    rec.n_allele = count()?;
    rec.n_fmt = count()?;
    rec.n_sample = count()?;
    let shared = usize::try_from(next()?).map_err(|_| malformed())?;
    let blocks = &body[p..];
    if shared > blocks.len() {
        return Err(malformed());
    }
    rec.shared.extend_from_slice(&blocks[..shared]);
    rec.indiv.extend_from_slice(&blocks[shared..]);
    Ok(len)
}

fn type_name(t: u8) -> &'static str {
    match t {
        0 => "null",
        1 => "int (8-bit)",
        2 => "int (16 bit)",
        3 => "int (32 bit)",
        5 => "float",
        7 => "char",
        _ => "unknown",
    }
}

static RLEN_WARNED: AtomicBool = AtomicBool::new(false);

/// Validates a record read from BCF against the header (HTSlib's `bcf_record_check`),
/// warning about bad IDs and types, fixing GT phasing for pre-4.4 headers, and repairing a
/// negative `rlen`.
///
/// # Errors
///
/// Returns an error when the record is malformed or refers to undefined IDs (the problems are
/// also recorded in [`Record::errors`]).
pub fn check_record(header: Option<&Header>, rec: &mut Record) -> Result<(), RecordError> {
    let mut err = RecordError::empty();
    let max_id = header.map_or(0, |h| h.id_count() as i64);
    let ids = header.map(|h| &h.ids[..]);
    // Whether a FILTER/INFO/FORMAT ID is unknown (always false without a header).
    let bad_id =
        |key: i64| key < 0 || ids.is_some_and(|ids| key >= max_id || ids[key as usize].is_none());
    let idgt = header
        .filter(|h| h.version() < VCF44)
        .and_then(|h| h.special_id(crate::variant::header::SpecialTag::Gt))
        .map_or(-1, i64::from);
    let loc = |rec: &Record| {
        format!(
            "{}:{}",
            String::from_utf8_lossy(
                header
                    .and_then(|h| h.contig_name(rec.rid))
                    .unwrap_or(b"(unknown)")
            ),
            rec.pos + 1
        )
    };
    let bad_shared = |rec: &Record| {
        log::error(
            "bcf_record_check",
            format_args!(
                "Bad BCF record at {} - shared section malformed or too short",
                loc(rec)
            ),
        );
        RecordError::empty()
    };
    let bad_indiv = |rec: &Record| {
        log::error(
            "bcf_record_check",
            format_args!(
                "Bad BCF record at {} - individuals section malformed or too short",
                loc(rec)
            ),
        );
        RecordError::empty()
    };

    if rec.rid < 0
        || header.is_some_and(|h| {
            rec.rid as usize >= h.contig_count() || h.contig_name(rec.rid).is_none()
        })
    {
        log::warn(
            "bcf_record_check",
            format_args!(
                "Bad BCF record at {}: Invalid CONTIG id {}",
                rec.pos + 1,
                rec.rid
            ),
        );
        err |= RecordError::CTG_INVALID;
    }

    let shared = &rec.shared;
    let mut p = 0usize;
    // ID
    let Some((num, ty, used)) = encoding::dec_size_safe(&shared[p..]) else {
        return Err(bad_shared(rec));
    };
    p += used;
    if ty != 7 {
        log::warn(
            "bcf_record_check",
            format_args!(
                "Bad BCF record at {}: Invalid ID type {ty} ({})",
                loc(rec),
                type_name(ty)
            ),
        );
        err |= RecordError::TAG_INVALID;
    }
    let bytes = (num as usize) << type_shift(ty);
    if shared.len() - p < bytes {
        return Err(bad_shared(rec));
    }
    p += bytes;

    // REF and ALT
    if rec.n_allele < 1 {
        log::warn(
            "bcf_record_check",
            format_args!("Bad BCF record at {}: No REF allele", loc(rec)),
        );
        err |= RecordError::TAG_UNDEF;
    }
    let mut reports = 0;
    for _ in 0..rec.n_allele {
        let Some((num, ty, used)) = encoding::dec_size_safe(&shared[p..]) else {
            return Err(bad_shared(rec));
        };
        p += used;
        if ty != 7 {
            if reports == 0 || log::verbosity() >= 5 {
                log::warn(
                    "bcf_record_check",
                    format_args!(
                        "Bad BCF record at {}: Invalid REF/ALT type {ty} ({})",
                        loc(rec),
                        type_name(ty)
                    ),
                );
            }
            reports += 1;
            err |= RecordError::CHAR;
        }
        let bytes = (num as usize) << type_shift(ty);
        if shared.len() - p < bytes {
            return Err(bad_shared(rec));
        }
        p += bytes;
    }

    // FILTER
    let mut reports = 0;
    let Some((num, ty, used)) = encoding::dec_size_safe(&shared[p..]) else {
        return Err(bad_shared(rec));
    };
    p += used;
    if num > 0 {
        let bytes = (num as usize) << type_shift(ty);
        if !matches!(ty, 1..=3) {
            log::warn(
                "bcf_record_check",
                format_args!(
                    "Bad BCF record at {}: Invalid FILTER type {ty} ({})",
                    loc(rec),
                    type_name(ty)
                ),
            );
            err |= RecordError::TAG_INVALID;
            if shared.len() - p < bytes {
                return Err(bad_shared(rec));
            }
            p += bytes;
        } else {
            if shared.len() - p < bytes {
                return Err(bad_shared(rec));
            }
            for _ in 0..num {
                let (key, used) = encoding::dec_int1(&shared[p..], ty);
                p += used;
                if bad_id(key) {
                    if reports == 0 || log::verbosity() >= 5 {
                        log::warn(
                            "bcf_record_check",
                            format_args!("Bad BCF record at {}: Invalid FILTER id {key}", loc(rec)),
                        );
                    }
                    reports += 1;
                    err |= RecordError::TAG_UNDEF;
                }
            }
        }
    }

    // INFO
    let mut reports = 0;
    for _ in 0..rec.n_info {
        let Some((key, used)) = encoding::dec_typed_int1_safe(&shared[p..]) else {
            return Err(bad_shared(rec));
        };
        p += used;
        let key = i64::from(key);
        if bad_id(key) {
            if reports == 0 || log::verbosity() >= 5 {
                log::warn(
                    "bcf_record_check",
                    format_args!("Bad BCF record at {}: Invalid INFO id {key}", loc(rec)),
                );
            }
            reports += 1;
            err |= RecordError::TAG_UNDEF;
        }
        let Some((num, ty, used)) = encoding::dec_size_safe(&shared[p..]) else {
            return Err(bad_shared(rec));
        };
        p += used;
        if !matches!(ty, 0..=3 | 5 | 7) || (ty == 0 && num > 0) {
            if reports == 0 || log::verbosity() >= 5 {
                log::warn(
                    "bcf_record_check",
                    format_args!(
                        "Bad BCF record at {}: Invalid INFO type {ty} ({})",
                        loc(rec),
                        type_name(ty)
                    ),
                );
            }
            reports += 1;
            err |= RecordError::TAG_INVALID;
        }
        let bytes = (num as usize) << type_shift(ty);
        if shared.len() - p < bytes {
            return Err(bad_shared(rec));
        }
        p += bytes;
    }

    // FORMAT and individual information
    let n_sample = rec.n_sample as usize;
    let mut p = 0usize;
    let mut reports = 0u32;
    let fmt_err = |rec: &Record, what: &str, i: i64, reports: &mut u32| {
        if *reports == 0 || log::verbosity() >= 5 {
            log::warn(
                "bcf_record_check_err",
                format_args!("Bad BCF record at {}: Invalid FORMAT {what} {i}", loc(rec)),
            );
        }
        *reports += 1;
    };
    for _ in 0..rec.n_fmt {
        let indiv = &rec.indiv;
        let Some((key, used)) = encoding::dec_typed_int1_safe(&indiv[p..]) else {
            return Err(bad_indiv(rec));
        };
        p += used;
        let key = i64::from(key);
        if bad_id(key) {
            fmt_err(rec, "id", key, &mut reports);
            err |= RecordError::TAG_UNDEF;
        }
        let Some((num, ty, used)) = encoding::dec_size_safe(&rec.indiv[p..]) else {
            return Err(bad_indiv(rec));
        };
        p += used;
        if !matches!(ty, 0..=3 | 5 | 7) || (ty == 0 && num > 0) {
            fmt_err(rec, "type", i64::from(ty), &mut reports);
            err |= RecordError::TAG_INVALID;
        }
        let ndata = num as u64 * n_sample as u64;
        if ndata > (i32::MAX as u64 >> type_shift(ty)) {
            log::error(
                "bcf_record_check",
                format_args!(
                    "Bad BCF record at {} - individuals section data too large",
                    loc(rec)
                ),
            );
            return Err(RecordError::empty());
        }
        let bytes = (ndata as usize) << type_shift(ty);
        if rec.indiv.len() - p < bytes {
            return Err(bad_indiv(rec));
        }
        if idgt >= 0
            && idgt == key
            && !update_phasing(&mut rec.indiv[p..p + bytes], n_sample, num as usize, ty)
        {
            err |= RecordError::TAG_INVALID;
        }
        p += bytes;
    }

    if err.is_empty() && rec.rlen < 0 {
        if !log::fired(&RLEN_WARNED) {
            log::warn_once(
                &RLEN_WARNED,
                "bcf_record_check",
                format_args!(
                    "BCF record at {} has invalid RLEN ({}). Only one invalid RLEN will be reported.",
                    loc(rec),
                    rec.rlen
                ),
            );
        }
        let len = header.map_or(0, |h| get_rlen(h, rec));
        rec.rlen = len.max(0);
    }
    rec.errcode |= err;
    if err.is_empty() { Ok(()) } else { Err(err) }
}

/// Sets the first allele's phase bit from the others, as VCF 4.4 implies it (HTSlib's
/// `updatephasing`), working on the low byte of each value.
fn update_phasing(p: &mut [u8], samples: usize, ploidy: usize, ty: u8) -> bool {
    let inc = 1usize << type_shift(ty);
    if p.len() < samples * ploidy * inc {
        return false;
    }
    let mut off = 0;
    match ploidy {
        0 => {}
        1 => {
            for _ in 0..samples {
                if p[off] != 0 {
                    p[off] |= 1;
                }
                off += inc;
            }
        }
        2 if inc == 1 => {
            // Four int8 diploid genotypes per word: each first allele takes the second's
            // phase bit.
            let (words, rest) = p[..samples * 2].as_chunks_mut::<8>();
            for w in words {
                let x = u64::from_le_bytes(*w);
                *w = (x | ((x >> 8) & 0x0001_0001_0001_0001)).to_le_bytes();
            }
            for c in rest.as_chunks_mut::<2>().0 {
                c[0] |= c[1] & 1;
            }
        }
        2 => {
            for c in p[..samples * 2 * inc].chunks_exact_mut(2 * inc) {
                c[0] |= c[inc] & 1;
            }
        }
        _ => {
            for _ in 0..samples {
                let mut allphased: u8 = 1;
                for k in 1..ploidy {
                    allphased &= p[off + inc * k];
                }
                p[off] |= allphased;
                off += ploidy * inc;
            }
        }
    }
    true
}

/// Reads and checks one record (HTSlib's `bcf_read`): raw read, [`check_record`], and the
/// header's read-time sample selection.
///
/// # Errors
///
/// Returns an error for a truncated or malformed record.
pub fn read_record<R: BufRead>(
    reader: &mut R,
    header: &Header,
    rec: &mut Record,
) -> io::Result<ReadOutcome> {
    if read_record_raw(reader, rec)? == ReadOutcome::Eof {
        return Ok(ReadOutcome::Eof);
    }
    if check_record(Some(header), rec).is_err() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "BCF read error"));
    }
    if header.keep_samples.is_some() {
        rec.subset_format(header);
    }
    Ok(ReadOutcome::Record)
}

/// An error encoding a record as BCF.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum EncodeError {
    /// The record's sample count differs from the header's.
    SampleCount,
    /// The record has unchecked parse errors (undefined contigs or tags).
    Unchecked(RecordError),
    /// The record holds 64-bit values.
    Requires64Bit,
    /// A count or length does not fit its BCF field.
    TooLarge(Limit),
}

/// A BCF field that a record's count or length can overflow.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Limit {
    /// More than 65535 alleles (`n_allele` has 16 bits).
    Alleles,
    /// More than 65535 INFO fields (`n_info` has 16 bits).
    Info,
    /// More than 255 FORMAT fields (`n_fmt` has 8 bits).
    Format,
    /// More than 2^24 - 1 samples (`n_sample` has 24 bits).
    Samples,
    /// A shared block of 2^31 bytes or more, with the 24 bytes of fixed fields its length
    /// counts.
    Shared,
    /// A per-sample block of 2^31 bytes or more.
    Indiv,
}

/// The counts and block lengths of a record, narrowed to the widths of their BCF fields.
#[derive(Debug, Eq, PartialEq)]
struct FixedCounts {
    l_shared: u32,
    l_indiv: u32,
    n_info: u16,
    n_allele: u16,
    /// `n_fmt` in the top 8 bits, `n_sample` in the other 24.
    n_fmt_sample: u32,
}

impl FixedCounts {
    /// Narrows the counts and block lengths, or names the first that does not fit.
    ///
    /// The blocks are held below 2^31 bytes rather than at their 32-bit fields: every typed
    /// value in a block is then shorter than 2^31 elements, the most its size (a signed 32-bit
    /// integer) holds.
    fn new(
        n_allele: u32,
        n_info: u32,
        n_fmt: u32,
        n_sample: u32,
        shared_len: usize,
        indiv_len: usize,
    ) -> Result<Self, Limit> {
        let n_allele = u16::try_from(n_allele).map_err(|_| Limit::Alleles)?;
        let n_info = u16::try_from(n_info).map_err(|_| Limit::Info)?;
        let n_fmt = u8::try_from(n_fmt).map_err(|_| Limit::Format)?;
        if n_sample > 0xff_ffff {
            return Err(Limit::Samples);
        }
        let l_shared = i32::try_from(shared_len + 24).map_err(|_| Limit::Shared)?;
        let l_indiv = i32::try_from(indiv_len).map_err(|_| Limit::Indiv)?;
        Ok(Self {
            l_shared: l_shared.cast_unsigned(),
            l_indiv: l_indiv.cast_unsigned(),
            n_info,
            n_allele,
            n_fmt_sample: (u32::from(n_fmt) << 24) | n_sample,
        })
    }
}

impl std::fmt::Display for EncodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SampleCount => f.write_str("the number of samples does not match the header"),
            Self::Unchecked(e) => write!(f, "unchecked error ({})", e.describe(1024)),
            Self::Requires64Bit => f.write_str("64-bit values cannot be written as BCF"),
            Self::TooLarge(limit) => write!(f, "too large for BCF ({limit:?})"),
        }
    }
}

impl std::error::Error for EncodeError {}

/// Checks that a record's sample count matches the output header (the first check of
/// HTSlib's `bcf_write`, for both VCF and BCF).
///
/// # Errors
///
/// Returns [`EncodeError::SampleCount`] (after reporting it) on a mismatch.
pub fn check_sample_count(header: &Header, rec: &Record) -> Result<(), EncodeError> {
    if header.sample_count() != rec.n_sample as usize {
        log::error(
            "bcf_write",
            format_args!(
                "Broken VCF record, the number of columns at {}:{} does not match the number of samples ({} vs {})",
                String::from_utf8_lossy(rec.chrom_safe(header)),
                rec.pos + 1,
                rec.n_sample,
                header.sample_count()
            ),
        );
        return Err(EncodeError::SampleCount);
    }
    Ok(())
}

/// Returns whether `bcf_write` refuses a record as needing 64 bits. HTSlib refuses only records
/// flagged `BCF_IS_64BIT`, which the VCF parser sets from POS 2^31 on, and writes any other
/// record with POS as an unsigned and rlen as a signed 32-bit integer. Those fields hold POS
/// up to 2^32 - 1 (0-based 2^32 - 2: `bcf_read` takes `u32::MAX` for -1) and rlen up to
/// `i32::MAX`, so larger values, which HTSlib truncates silently in records that `norm`,
/// `merge`, or `sort` rebuilt, are refused here as well.
fn needs_64bit(v: &Record) -> bool {
    v.unpacked & IS_64BIT != 0 || v.pos >= i64::from(u32::MAX) || v.rlen > i64::from(i32::MAX)
}

/// Appends a record as BCF (the BCF part of HTSlib's `bcf_write`), syncing edits first.
///
/// # Errors
///
/// Returns an error (after reporting it) for a sample-count mismatch, unchecked parse errors,
/// 64-bit values, or counts and lengths that do not fit their BCF fields.
pub fn encode_record(
    header: &Header,
    v: &mut Record,
    out: &mut Vec<u8>,
) -> Result<(), EncodeError> {
    check_sample_count(header, v)?;
    let unchecked = v.errcode - RecordError::LIMITS;
    if !unchecked.is_empty() {
        log::error(
            "bcf_write",
            format_args!(
                "Unchecked error ({} {}) at {}:{}",
                v.errcode.bits(),
                v.errcode.describe(1024),
                String::from_utf8_lossy(v.chrom_safe(header)),
                v.pos + 1
            ),
        );
        return Err(EncodeError::Unchecked(v.errcode));
    }
    v.sync();
    if needs_64bit(v) {
        log::error(
            "bcf_write",
            format_args!(
                "Data at {}:{} contains 64-bit values not representable in BCF. Please use VCF instead",
                String::from_utf8_lossy(v.chrom_safe(header)),
                v.pos + 1
            ),
        );
        return Err(EncodeError::Requires64Bit);
    }
    // HTSlib keeps the counts in bit fields of these widths and writes the lengths as 32-bit
    // integers, so larger values (which edits such as `norm -m+` or `annotate` can make) wrap
    // into a corrupt record, or crash bcftools first.
    let counts = match FixedCounts::new(
        v.n_allele,
        v.n_info,
        v.n_fmt,
        v.n_sample,
        v.shared.len(),
        v.indiv.len(),
    ) {
        Ok(counts) => counts,
        Err(limit) => {
            let at = format!(
                "{}:{}",
                String::from_utf8_lossy(v.chrom_safe(header)),
                v.pos + 1
            );
            let message = match limit {
                Limit::Alleles => format!("Too many ALT alleles at {at}"),
                Limit::Info => format!("Too many INFO entries at {at}"),
                Limit::Format => format!("Too many FORMAT fields at {at}"),
                Limit::Samples => format!("Too many samples at {at}"),
                Limit::Shared => format!("The record at {at} is too long"),
                Limit::Indiv => format!("The FORMAT at {at} is too long"),
            };
            log::error("bcf_write", format_args!("{message}"));
            return Err(EncodeError::TooLarge(limit));
        }
    };
    out.reserve(32 + v.shared.len() + v.indiv.len());
    out.extend_from_slice(&counts.l_shared.to_le_bytes());
    out.extend_from_slice(&counts.l_indiv.to_le_bytes());
    out.extend_from_slice(&v.rid.to_le_bytes());
    // POS -1 is written as u32::MAX, as HTSlib does; `needs_64bit` refused larger values.
    out.extend_from_slice(&(v.pos as u32).to_le_bytes());
    out.extend_from_slice(&(v.rlen as u32).to_le_bytes());
    out.extend_from_slice(&v.qual.to_bits().to_le_bytes());
    out.extend_from_slice(&counts.n_info.to_le_bytes());
    out.extend_from_slice(&counts.n_allele.to_le_bytes());
    out.extend_from_slice(&counts.n_fmt_sample.to_le_bytes());
    out.extend_from_slice(&v.shared);
    out.extend_from_slice(&v.indiv);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_records_round_trip() {
        let text = b"##fileformat=VCFv4.2\n##contig=<ID=1>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let mut header = Header::parse(text).expect("the header parses");
        let mut line = b"1\t3000000000\trs1\tAC\tA,T\t5\t.\tDP=7\tGT\t0|1".to_vec();
        let mut rec = Record::new();
        crate::variant::vcf::parse(&mut header, &mut line, &mut rec).expect("the line parses");
        let mut packed = vec![0xaa];
        pack_record(&mut rec, &mut packed);
        assert_eq!(packed_len(&packed[1..]), Some(packed.len() - 1));
        let mut back = Record::new();
        assert_eq!(
            unpack_record(&packed[1..], &mut back).unwrap(),
            packed.len() - 1
        );
        let (mut a, mut b) = (Vec::new(), Vec::new());
        crate::variant::vcf::format_record(&header, &rec, &mut a).unwrap();
        crate::variant::vcf::format_record(&header, &back, &mut b).unwrap();
        assert_eq!(a, b);
        assert_eq!(back.pos(), 2_999_999_999);
        assert!(unpack_record(&packed[1..packed.len() - 1], &mut back).is_err());
    }

    #[test]
    fn packed_records_keep_counts_and_lengths_that_take_long_varints() {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
             ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        (0..200).for_each(|i| text.push_str(&format!("\ts{i}")));
        text.push('\n');
        let mut header = Header::parse(text.as_bytes()).expect("the header parses");
        let alts = (0..130)
            .map(|i| format!("C{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let mut line = format!(
            "1\t5\t.\tA\t{alts}\t.\t.\tEND=20000\tDP{}",
            "\t300".repeat(200)
        )
        .into_bytes();
        let mut rec = Record::new();
        crate::variant::vcf::parse(&mut header, &mut line, &mut rec).expect("the line parses");
        let mut packed = Vec::new();
        pack_record(&mut rec, &mut packed);
        let mut back = Record::new();
        assert_eq!(unpack_record(&packed, &mut back).unwrap(), packed.len());
        assert_eq!(
            (back.rlen(), back.allele_count(), back.sample_count()),
            (19996, 131, 200)
        );
        assert_eq!(back.shared_bytes(), rec.shared_bytes());
        assert_eq!(back.indiv_bytes(), rec.indiv_bytes());
        // Entries too short for the fixed fields, with a varint that does not end, or with a
        // shared block longer than the rest are refused.
        let mut short = packed.clone();
        short[..4].copy_from_slice(&10u32.to_le_bytes());
        let mut endless = packed.clone();
        endless[20..31].fill(0xff);
        let mut long_shared = 32u32.to_le_bytes().to_vec();
        long_shared.extend_from_slice(&[0; 16]);
        long_shared.extend_from_slice(&[2, 0, 1, 0, 0, 50]);
        long_shared.extend_from_slice(&[0; 10]);
        for bad in [short, endless, long_shared] {
            assert!(unpack_record(&bad, &mut back).is_err());
        }
    }

    #[test]
    fn positions_beyond_32_bits_are_refused_after_copies_and_edits() {
        let text = b"##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let mut header = Header::parse(text).expect("the header parses");
        let mut line = b"1\t3000000000\t.\tA\tC\t.\t.\t.".to_vec();
        let mut flagged = Record::new();
        crate::variant::vcf::parse(&mut header, &mut line, &mut flagged).expect("the line parses");
        assert!(flagged.is_64bit());
        let mut line = b"1\t5\t.\tA\tC\t.\t.\t.".to_vec();
        let mut small = Record::new();
        crate::variant::vcf::parse(&mut header, &mut line, &mut small).expect("the line parses");

        // Records without the flag are written while BCF's fields read them back unchanged:
        // POS as an unsigned 32-bit integer short of u32::MAX (-1), rlen as a signed one.
        let last_pos = i64::from(u32::MAX) - 1;
        for pos in [
            i64::from(i32::MAX) - 1,
            i64::from(i32::MAX),
            2_999_999_999,
            last_pos,
        ] {
            let mut rec = small.clone();
            rec.set_pos(pos);
            let mut out = Vec::new();
            encode_record(&header, &mut rec, &mut out).expect("a 32-bit POS encodes");
            let mut back = Record::new();
            read_record_raw(&mut out.as_slice(), &mut back).expect("the record reads back");
            assert_eq!(back.pos(), pos);
        }
        let mut long = small.clone();
        long.set_rlen(i64::from(i32::MAX));
        encode_record(&header, &mut long, &mut Vec::new()).expect("rlen i32::MAX encodes");

        let mut copy = Record::new();
        copy.copy_from(&mut flagged);
        assert!(copy.is_64bit());
        let mut wrapped = small.clone();
        wrapped.set_pos(i64::from(u32::MAX));
        let mut truncated = small.clone();
        truncated.set_pos(4_294_967_299);
        let mut longer = small.clone();
        longer.set_rlen(i64::from(i32::MAX) + 1);
        for mut rec in [flagged, copy, wrapped, truncated, longer] {
            let mut out = Vec::new();
            assert_eq!(
                encode_record(&header, &mut rec, &mut out),
                Err(EncodeError::Requires64Bit),
                "pos {} rlen {}",
                rec.pos(),
                rec.rlen()
            );
            assert!(out.is_empty());
        }
    }

    #[test]
    fn packed_records_drop_the_64bit_flag() {
        let text = b"##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let mut header = Header::parse(text).expect("the header parses");
        for (pos, encodes) in [(3_000_000_000_i64, true), (4_294_967_300, false)] {
            let mut line = format!("1\t{pos}\t.\tA\tC\t.\t.\t.").into_bytes();
            let mut rec = Record::new();
            crate::variant::vcf::parse(&mut header, &mut line, &mut rec).expect("parses");
            let mut packed = Vec::new();
            pack_record(&mut rec, &mut packed);
            let mut back = Record::new();
            unpack_record(&packed, &mut back).expect("unpacks");
            assert!(!back.is_64bit());
            let encoded = encode_record(&header, &mut back, &mut Vec::new());
            assert_eq!(encoded.is_ok(), encodes, "POS {pos}");
        }
    }

    #[test]
    fn counts_and_lengths_beyond_their_fields_are_refused() {
        let fits = |a, i, f, s, ls, li| FixedCounts::new(a, i, f, s, ls, li).err();
        let block = i32::MAX as usize;
        assert_eq!(fits(65535, 65535, 255, 0xff_ffff, block - 24, block), None);
        assert_eq!(fits(65536, 0, 0, 0, 0, 0), Some(Limit::Alleles));
        assert_eq!(fits(1, 65536, 0, 0, 0, 0), Some(Limit::Info));
        assert_eq!(fits(1, 0, 256, 1, 0, 0), Some(Limit::Format));
        assert_eq!(fits(1, 0, 1, 0x100_0000, 0, 0), Some(Limit::Samples));
        assert_eq!(fits(1, 0, 0, 0, block - 23, 0), Some(Limit::Shared));
        assert_eq!(fits(1, 0, 1, 1, 0, block + 1), Some(Limit::Indiv));
        let packed = FixedCounts::new(3, 2, 255, 0xff_ffff, 10, 20).expect("fits");
        assert_eq!(
            (packed.l_shared, packed.l_indiv, packed.n_fmt_sample),
            (34, 20, u32::MAX)
        );

        // Edits that pass a field's width (bcftools wraps or crashes) are not written.
        let mut text = String::from("##fileformat=VCFv4.2\n##contig=<ID=1>\n");
        for i in 0..65536 {
            text.push_str(&format!(
                "##INFO=<ID=F{i},Number=0,Type=Flag,Description=\"f\">\n"
            ));
        }
        for i in 0..256 {
            text.push_str(&format!(
                "##FORMAT=<ID=X{i},Number=1,Type=Integer,Description=\"x\">\n"
            ));
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n");
        let mut header = Header::parse(text.as_bytes()).expect("the header parses");
        let flags = (0..65535).map(|i| format!("F{i}")).collect::<Vec<_>>();
        let keys = (0..255).map(|i| format!("X{i}")).collect::<Vec<_>>();
        let mut line = format!(
            "1\t5\t.\tA\tC\t.\t.\t{}\t{}\t{}",
            flags.join(";"),
            keys.join(":"),
            vec!["1"; 255].join(":")
        )
        .into_bytes();
        let mut base = Record::new();
        crate::variant::vcf::parse(&mut header, &mut line, &mut base).expect("the line parses");
        encode_record(&header, &mut base.clone(), &mut Vec::new())
            .expect("65535 INFO and 255 FORMAT fields fit");

        let mut info = base.clone();
        info.update_info_flag(&header, b"F65535", true)
            .expect("the flag is added");
        let mut alleles = base.clone();
        let als = (0..65536).map(|i| format!("A{i:b}")).collect::<Vec<_>>();
        alleles.update_alleles(&header, &als);
        let mut format = base.clone();
        format
            .update_format_i32(&header, b"X255", &[1])
            .expect("the field is added");
        for (mut rec, limit) in [
            (info, Limit::Info),
            (alleles, Limit::Alleles),
            (format, Limit::Format),
        ] {
            let mut out = Vec::new();
            assert_eq!(
                encode_record(&header, &mut rec, &mut out),
                Err(EncodeError::TooLarge(limit))
            );
            assert!(out.is_empty());
        }
    }

    #[test]
    fn diploid_phasing_matches_the_scalar_rule() {
        // 11 int8 samples: two full words and a remainder of three.
        let mut gts: Vec<u8> = (0..22u8).map(|i| i.wrapping_mul(37) ^ 0x55).collect();
        let want: Vec<u8> = gts
            .chunks(2)
            .flat_map(|c| [c[0] | (c[1] & 1), c[1]])
            .collect();
        assert!(update_phasing(&mut gts, 11, 2, 1));
        assert_eq!(gts, want);

        // int16 values change only their low byte.
        let mut wide = vec![2u8, 0, 5, 0, 4, 0, 4, 0];
        assert!(update_phasing(&mut wide, 2, 2, 2));
        assert_eq!(wide, [3, 0, 5, 0, 4, 0, 4, 0]);
    }
}
