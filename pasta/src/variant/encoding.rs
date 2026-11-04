//! BCF typed-value encoding: type codes, missing and vector-end sentinels, and the
//! encoders/decoders HTSlib uses (`bcf_enc_*`, `bcf_dec_*`).

/// A BCF value type code (HTSlib's `BCF_BT_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[repr(u8)]
pub enum BcfType {
    /// No value (`BCF_BT_NULL`).
    Null = 0,
    /// 8-bit integers (`BCF_BT_INT8`).
    Int8 = 1,
    /// 16-bit integers (`BCF_BT_INT16`).
    Int16 = 2,
    /// 32-bit integers (`BCF_BT_INT32`).
    Int32 = 3,
    /// 64-bit integers (`BCF_BT_INT64`), internal to HTSlib.
    Int64 = 4,
    /// 32-bit floats (`BCF_BT_FLOAT`).
    Float = 5,
    /// Characters (`BCF_BT_CHAR`).
    Char = 7,
}

impl BcfType {
    /// Converts a raw 4-bit type code.
    #[inline]
    #[must_use]
    pub fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            0 => Self::Null,
            1 => Self::Int8,
            2 => Self::Int16,
            3 => Self::Int32,
            4 => Self::Int64,
            5 => Self::Float,
            7 => Self::Char,
            _ => return None,
        })
    }

    /// Returns the size in bytes of one value (`1 << bcf_type_shift[type]`).
    #[inline]
    #[must_use]
    pub fn size(self) -> usize {
        1 << type_shift(self as u8)
    }

    /// Returns whether this is an integer type.
    #[inline]
    #[must_use]
    pub fn is_int(self) -> bool {
        matches!(self, Self::Int8 | Self::Int16 | Self::Int32 | Self::Int64)
    }
}

/// HTSlib's `bcf_type_shift` table, indexed by raw type code.
pub const TYPE_SHIFT: [u8; 16] = [0, 0, 1, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

/// Returns `bcf_type_shift[code & 0xf]`.
#[inline]
#[must_use]
pub fn type_shift(code: u8) -> u32 {
    u32::from(TYPE_SHIFT[usize::from(code & 0xf)])
}

/// The 8-bit integer missing value.
pub const INT8_MISSING: i8 = i8::MIN;
/// The 8-bit integer vector-end value.
pub const INT8_VECTOR_END: i8 = i8::MIN + 1;
/// The 16-bit integer missing value.
pub const INT16_MISSING: i16 = i16::MIN;
/// The 16-bit integer vector-end value.
pub const INT16_VECTOR_END: i16 = i16::MIN + 1;
/// The 32-bit integer missing value (`bcf_int32_missing`).
pub const INT32_MISSING: i32 = i32::MIN;
/// The 32-bit integer vector-end value (`bcf_int32_vector_end`).
pub const INT32_VECTOR_END: i32 = i32::MIN + 1;
/// The 64-bit integer missing value.
pub const INT64_MISSING: i64 = i64::MIN;
/// The 64-bit integer vector-end value.
pub const INT64_VECTOR_END: i64 = i64::MIN + 1;
/// The bits of the float missing value (`bcf_float_missing`), a NaN payload.
pub const FLOAT_MISSING: u32 = 0x7F80_0001;
/// The bits of the float vector-end value (`bcf_float_vector_end`), a NaN payload.
pub const FLOAT_VECTOR_END: u32 = 0x7F80_0002;
/// The string missing value (`bcf_str_missing`).
pub const STR_MISSING: u8 = 7;

/// The largest 8-bit value that is not a sentinel.
pub const MAX_INT8: i32 = 0x7f;
/// The smallest 8-bit value that is not a sentinel.
pub const MIN_INT8: i32 = -120;
/// The largest 16-bit value that is not a sentinel.
pub const MAX_INT16: i32 = 0x7fff;
/// The smallest 16-bit value that is not a sentinel.
pub const MIN_INT16: i32 = -32760;
/// The largest 32-bit value that is not a sentinel.
pub const MAX_INT32: i32 = 0x7fff_ffff;
/// The smallest 32-bit value that is not a sentinel.
pub const MIN_INT32: i32 = -2_147_483_640;

/// Returns the float missing value.
#[inline]
#[must_use]
pub fn float_missing() -> f32 {
    f32::from_bits(FLOAT_MISSING)
}

/// Returns the float vector-end value.
#[inline]
#[must_use]
pub fn float_vector_end() -> f32 {
    f32::from_bits(FLOAT_VECTOR_END)
}

/// Returns whether `f` is the float missing value (compared bitwise, as HTSlib does).
#[inline]
#[must_use]
pub fn float_is_missing(f: f32) -> bool {
    f.to_bits() == FLOAT_MISSING
}

/// Returns whether `f` is the float vector-end value.
#[inline]
#[must_use]
pub fn float_is_vector_end(f: f32) -> bool {
    f.to_bits() == FLOAT_VECTOR_END
}

/// Reads a little-endian `i16`.
#[inline(always)]
#[must_use]
pub fn le_i16(p: &[u8]) -> i16 {
    i16::from_le_bytes([p[0], p[1]])
}

/// Reads a little-endian `i32`.
#[inline(always)]
#[must_use]
pub fn le_i32(p: &[u8]) -> i32 {
    i32::from_le_bytes([p[0], p[1], p[2], p[3]])
}

/// Reads a little-endian `u32`.
#[inline(always)]
#[must_use]
pub fn le_u32(p: &[u8]) -> u32 {
    u32::from_le_bytes([p[0], p[1], p[2], p[3]])
}

/// Reads a little-endian `i64`.
#[inline(always)]
#[must_use]
pub fn le_i64(p: &[u8]) -> i64 {
    i64::from_le_bytes([p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]])
}

/// Appends a type descriptor (HTSlib's `bcf_enc_size`).
#[inline]
pub fn enc_size(s: &mut Vec<u8>, size: usize, ty: u8) {
    if size < 15 {
        s.push(((size as u8) << 4) | ty);
        return;
    }
    s.push((15 << 4) | ty);
    if size < 128 {
        s.push((1 << 4) | BcfType::Int8 as u8);
        s.push(size as u8);
    } else if size < 32768 {
        s.push((1 << 4) | BcfType::Int16 as u8);
        s.extend_from_slice(&(size as i16).to_le_bytes());
    } else {
        // A size of 2^31 or more does not fit: it is written negative, which BCF readers
        // refuse, and `bcf::encode_record` refuses the blocks long enough to hold such values.
        s.push((1 << 4) | BcfType::Int32 as u8);
        s.extend_from_slice(&i32::try_from(size).unwrap_or(i32::MIN).to_le_bytes());
    }
}

/// Returns the smallest integer type that holds `x` (HTSlib's `bcf_enc_inttype`).
#[inline]
#[must_use]
pub fn enc_inttype(x: i64) -> BcfType {
    if x <= i64::from(MAX_INT8) && x >= i64::from(MIN_INT8) {
        BcfType::Int8
    } else if x <= i64::from(MAX_INT16) && x >= i64::from(MIN_INT16) {
        BcfType::Int16
    } else {
        BcfType::Int32
    }
}

/// Appends a single typed integer (HTSlib's `bcf_enc_int1`).
#[inline]
pub fn enc_int1(s: &mut Vec<u8>, x: i32) {
    if x == INT32_VECTOR_END {
        s.push((1 << 4) | BcfType::Int8 as u8);
        s.push(INT8_VECTOR_END as u8);
    } else if x == INT32_MISSING {
        s.push((1 << 4) | BcfType::Int8 as u8);
        s.push(INT8_MISSING as u8);
    } else if (MIN_INT8..=MAX_INT8).contains(&x) {
        s.push((1 << 4) | BcfType::Int8 as u8);
        s.push(x as u8);
    } else if (MIN_INT16..=MAX_INT16).contains(&x) {
        s.push((1 << 4) | BcfType::Int16 as u8);
        s.extend_from_slice(&(x as i16).to_le_bytes());
    } else {
        s.push((1 << 4) | BcfType::Int32 as u8);
        s.extend_from_slice(&x.to_le_bytes());
    }
}

/// Appends a typed integer vector (HTSlib's `bcf_enc_vint`). `wsize` is the vector size stored
/// in the descriptor (the values per sample for FORMAT data), or `None` for `a.len()`.
pub fn enc_vint(s: &mut Vec<u8>, a: &[i32], wsize: Option<usize>) {
    let n = a.len();
    if n == 0 {
        enc_size(s, 0, BcfType::Null as u8);
        return;
    }
    if n == 1 {
        enc_int1(s, a[0]);
        return;
    }
    let wsize = wsize.filter(|&w| w > 0).unwrap_or(n);
    let max = a.iter().copied().fold(i32::MIN, i32::max);
    let min = a
        .iter()
        .copied()
        .fold(i32::MAX, |m, x| if x > i32::MIN + 1 { m.min(x) } else { m });
    if max <= MAX_INT8 && min >= MIN_INT8 {
        enc_size(s, wsize, BcfType::Int8 as u8);
        s.extend(a.iter().map(|&x| {
            if x == INT32_VECTOR_END {
                INT8_VECTOR_END as u8
            } else if x == INT32_MISSING {
                INT8_MISSING as u8
            } else {
                x as u8
            }
        }));
    } else if max <= MAX_INT16 && min >= MIN_INT16 {
        enc_size(s, wsize, BcfType::Int16 as u8);
        s.reserve(n * 2);
        for &x in a {
            let v = if x == INT32_VECTOR_END {
                INT16_VECTOR_END
            } else if x == INT32_MISSING {
                INT16_MISSING
            } else {
                x as i16
            };
            s.extend_from_slice(&v.to_le_bytes());
        }
    } else {
        enc_size(s, wsize, BcfType::Int32 as u8);
        s.reserve(n * 4);
        for &x in a {
            s.extend_from_slice(&x.to_le_bytes());
        }
    }
}

/// Appends a typed float vector (HTSlib's `bcf_enc_vfloat`).
pub fn enc_vfloat(s: &mut Vec<u8>, a: &[f32]) {
    enc_size(s, a.len(), BcfType::Float as u8);
    push_floats(s, a);
}

/// Appends floats in little-endian byte order (HTSlib's `serialize_float_array`).
#[inline]
pub fn push_floats(s: &mut Vec<u8>, a: &[f32]) {
    s.reserve(a.len() * 4);
    for &f in a {
        s.extend_from_slice(&f.to_bits().to_le_bytes());
    }
}

/// Appends a typed character vector (HTSlib's `bcf_enc_vchar`).
#[inline]
pub fn enc_vchar(s: &mut Vec<u8>, a: &[u8]) {
    enc_size(s, a.len(), BcfType::Char as u8);
    s.extend_from_slice(a);
}

/// Decodes one integer of raw type `ty` at `p[0..]` (HTSlib's `bcf_dec_int1`), returning the
/// value and its size, or `(0, 0)` for a non-integer type.
#[inline(always)]
#[must_use]
pub fn dec_int1(p: &[u8], ty: u8) -> (i64, usize) {
    match ty {
        1 => (i64::from(p[0] as i8), 1),
        2 => (i64::from(le_i16(p)), 2),
        3 => (i64::from(le_i32(p)), 4),
        4 => (le_i64(p), 8),
        _ => (0, 0),
    }
}

/// Decodes a typed integer (HTSlib's `bcf_dec_typed_int1`), returning the value and the number
/// of bytes consumed (1 for an invalid type, which consumes only the type byte).
#[inline(always)]
#[must_use]
pub fn dec_typed_int1(p: &[u8]) -> (i64, usize) {
    let (v, n) = dec_int1(&p[1..], p[0] & 0xf);
    (v, 1 + n)
}

/// Decodes a type descriptor (HTSlib's `bcf_dec_size`), returning the vector size, the raw type
/// code, and the number of bytes consumed.
#[inline(always)]
#[must_use]
pub fn dec_size(p: &[u8]) -> (i64, u8, usize) {
    let ty = p[0] & 0xf;
    if p[0] >> 4 != 15 {
        (i64::from(p[0] >> 4), ty, 1)
    } else {
        let (v, n) = dec_typed_int1(&p[1..]);
        (v, ty, 1 + n)
    }
}

/// Bounds-checked [`dec_typed_int1`] (HTSlib's `bcf_dec_typed_int1_safe`): the value and bytes
/// consumed, or `None` when the data is short or the type is not an 8-, 16- or 32-bit integer.
#[inline(always)]
#[must_use]
pub fn dec_typed_int1_safe(p: &[u8]) -> Option<(i32, usize)> {
    if p.len() < 2 {
        return None;
    }
    let t = p[0] & 0xf;
    if t == 1 {
        return Some((i32::from(p[1] as i8), 2));
    }
    let need = 1usize << type_shift(t);
    if p.len() - 1 < need {
        return None;
    }
    match t {
        2 => Some((i32::from(le_i16(&p[1..])), 3)),
        3 => Some((le_i32(&p[1..]), 5)),
        _ => None,
    }
}

/// Bounds-checked [`dec_size`] (HTSlib's `bcf_dec_size_safe`): the size, raw type, and bytes
/// consumed, or `None` for short data or a negative size.
#[inline(always)]
#[must_use]
pub fn dec_size_safe(p: &[u8]) -> Option<(i32, u8, usize)> {
    let first = *p.first()?;
    let ty = first & 0xf;
    if first >> 4 != 15 {
        return Some((i32::from(first >> 4), ty, 1));
    }
    let (num, n) = dec_typed_int1_safe(&p[1..])?;
    if num < 0 {
        return None;
    }
    Some((num, ty, 1 + n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn size_round_trip() {
        for n in [0usize, 1, 14, 15, 127, 128, 32767, 32768, 100_000] {
            let mut s = Vec::new();
            enc_size(&mut s, n, BcfType::Char as u8);
            let (got, ty, used) = dec_size(&s);
            assert_eq!((got as usize, ty, used), (n, 7, s.len()), "size {n}");
        }
    }

    #[test]
    fn vint_picks_smallest_type() {
        let mut s = Vec::new();
        enc_vint(&mut s, &[1, INT32_MISSING, INT32_VECTOR_END], None);
        assert_eq!(s, [0x31, 1, 0x80, 0x81]);
        s.clear();
        enc_vint(&mut s, &[1, 300], None);
        assert_eq!(s[0], 0x22);
        s.clear();
        enc_vint(&mut s, &[-121, 0], None);
        assert_eq!(s[0], 0x22);
    }
}
