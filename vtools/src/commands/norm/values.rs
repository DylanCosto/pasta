//! Integer and float INFO/FORMAT values handled alike, as bcftools' `BRANCH_NUMERIC` macros
//! handle `int32_t` and `float`.

use pasta::variant::{GetError, Header, Record, UpdateError, encoding};

use super::split::FormatStrings;

/// Reusable buffers of one value type.
#[derive(Debug, Default)]
pub(super) struct Buffers<T> {
    pub(super) a: Vec<T>,
    pub(super) b: Vec<T>,
    pub(super) out: Vec<T>,
}

/// Reusable buffers for INFO/FORMAT values.
#[derive(Debug, Default)]
pub(super) struct Scratch {
    pub(super) i: Buffers<i32>,
    pub(super) f: Buffers<f32>,
    pub(super) bytes: Vec<u8>,
    pub(super) bytes2: Vec<u8>,
    pub(super) strings: Vec<Vec<u8>>,
    pub(super) info_keys: Vec<i32>,
    pub(super) fmt_ids: Vec<i32>,
    /// Per sample: whether a Number=G vector is haploid.
    pub(super) haploid: Vec<bool>,
    pub(super) format_strings: FormatStrings,
}

/// A numeric INFO/FORMAT value: `i32` or `f32`, with HTSlib's missing and vector-end
/// sentinels.
pub(super) trait Value: Copy + std::fmt::Debug {
    fn missing() -> Self;
    fn vector_end() -> Self;
    fn is_missing(self) -> bool;
    fn is_vector_end(self) -> bool;
    /// `a += b` as C does it (wrapping for integers).
    fn add(self, other: Self) -> Self;
    fn buffers(s: &mut Scratch) -> &mut Buffers<Self>;
    fn get_info(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        out: &mut Vec<Self>,
    ) -> Result<usize, GetError>;
    fn get_format(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        out: &mut Vec<Self>,
    ) -> Result<usize, GetError>;
    fn update_info(rec: &mut Record, h: &Header, tag: &[u8], v: &[Self])
    -> Result<(), UpdateError>;
    fn update_format(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        v: &[Self],
    ) -> Result<(), UpdateError>;
}

impl Value for i32 {
    #[inline]
    fn missing() -> Self {
        encoding::INT32_MISSING
    }
    #[inline]
    fn vector_end() -> Self {
        encoding::INT32_VECTOR_END
    }
    #[inline]
    fn is_missing(self) -> bool {
        self == encoding::INT32_MISSING
    }
    #[inline]
    fn is_vector_end(self) -> bool {
        self == encoding::INT32_VECTOR_END
    }
    #[inline]
    fn add(self, other: Self) -> Self {
        self.wrapping_add(other)
    }
    fn buffers(s: &mut Scratch) -> &mut Buffers<Self> {
        &mut s.i
    }
    fn get_info(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        out: &mut Vec<Self>,
    ) -> Result<usize, GetError> {
        rec.get_info_i32(h, tag, out)
    }
    fn get_format(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        out: &mut Vec<Self>,
    ) -> Result<usize, GetError> {
        rec.get_format_i32(h, tag, out)
    }
    fn update_info(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        v: &[Self],
    ) -> Result<(), UpdateError> {
        rec.update_info_i32(h, tag, v)
    }
    fn update_format(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        v: &[Self],
    ) -> Result<(), UpdateError> {
        rec.update_format_i32(h, tag, v)
    }
}

impl Value for f32 {
    #[inline]
    fn missing() -> Self {
        encoding::float_missing()
    }
    #[inline]
    fn vector_end() -> Self {
        encoding::float_vector_end()
    }
    #[inline]
    fn is_missing(self) -> bool {
        self.to_bits() == encoding::FLOAT_MISSING
    }
    #[inline]
    fn is_vector_end(self) -> bool {
        self.to_bits() == encoding::FLOAT_VECTOR_END
    }
    #[inline]
    fn add(self, other: Self) -> Self {
        self + other
    }
    fn buffers(s: &mut Scratch) -> &mut Buffers<Self> {
        &mut s.f
    }
    fn get_info(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        out: &mut Vec<Self>,
    ) -> Result<usize, GetError> {
        rec.get_info_f32(h, tag, out)
    }
    fn get_format(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        out: &mut Vec<Self>,
    ) -> Result<usize, GetError> {
        rec.get_format_f32(h, tag, out)
    }
    fn update_info(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        v: &[Self],
    ) -> Result<(), UpdateError> {
        rec.update_info_f32(h, tag, v)
    }
    fn update_format(
        rec: &mut Record,
        h: &Header,
        tag: &[u8],
        v: &[Self],
    ) -> Result<(), UpdateError> {
        rec.update_format_f32(h, tag, v)
    }
}

/// Finds the `nth` comma-separated field of `s[start..]`, stopping at a NUL, as bcftools'
/// `string_move_nth` does. Returns the field's start and end (the end is where the scan
/// stopped: a comma, a NUL, or the end), or `None` when the field is missing or empty at the
/// end of the string.
pub(super) fn nth_field(s: &[u8], start: usize, nth: usize) -> Option<(usize, usize)> {
    if start >= s.len() {
        return None;
    }
    let (mut ss, mut se) = (start, start);
    let mut j = 0;
    while se < s.len() && s[se] != 0 {
        if s[se] == b',' {
            if j == nth {
                break;
            }
            j += 1;
            ss = se + 1;
        }
        se += 1;
    }
    if j != nth || ss >= s.len() || s[ss] == 0 {
        return None;
    }
    Some((ss, se))
}

/// Copies field `isrc` of the comma-separated `src` into field `idst` of `dst` when that is
/// still missing (`.`), as bcftools' `copy_string_field` (vcfmerge.c) does. Returns `false`
/// when either field does not exist.
pub(super) fn copy_string_field(src: &[u8], isrc: usize, dst: &mut Vec<u8>, idst: usize) -> bool {
    let (mut ith, mut start) = (0, 0);
    while ith < isrc && start < src.len() {
        if src[start] == b',' {
            ith += 1;
        }
        start += 1;
    }
    if ith != isrc {
        return false;
    }
    let mut end = start;
    while end < src.len() && src[end] != 0 && src[end] != b',' {
        end += 1;
    }
    if end - start == 1 && src[start] == b'.' {
        return true;
    }
    let (mut ith, mut sd) = (0, 0);
    while ith < idst && sd < dst.len() {
        if dst[sd] == b',' {
            ith += 1;
        }
        sd += 1;
    }
    if ith != idst {
        return false;
    }
    let mut ed = sd;
    while ed < dst.len() && dst[ed] != b',' {
        ed += 1;
    }
    if ed - sd > 1 || dst.get(sd) != Some(&b'.') {
        return true;
    }
    dst.splice(sd..ed, src[start..end].iter().copied());
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nth_field_matches_string_move_nth() {
        let s = b"a,bb,,c";
        assert_eq!(nth_field(s, 0, 0), Some((0, 1)));
        assert_eq!(nth_field(s, 0, 1), Some((2, 4)));
        assert_eq!(nth_field(s, 0, 2), Some((5, 5)));
        assert_eq!(nth_field(s, 0, 3), Some((6, 7)));
        assert_eq!(nth_field(s, 0, 4), None);
        assert_eq!(nth_field(b"a,", 0, 1), None);
        assert_eq!(nth_field(b"a\0,b", 0, 1), None);
    }

    #[test]
    fn copy_string_field_fills_missing_only() {
        let mut dst = b".,.,.".to_vec();
        assert!(copy_string_field(b"x,yy", 1, &mut dst, 2));
        assert_eq!(dst, b".,.,yy");
        assert!(copy_string_field(b"z", 0, &mut dst, 2));
        assert_eq!(dst, b".,.,yy");
        assert!(copy_string_field(b".", 0, &mut dst, 0));
        assert_eq!(dst, b".,.,yy");
        assert!(!copy_string_field(b"a", 3, &mut dst, 0));
    }
}
