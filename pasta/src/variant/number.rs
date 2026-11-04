//! Number parsing and formatting with HTSlib's exact semantics.
//!
//! VCF output must reproduce HTSlib byte for byte, so floats are formatted with `kputd`'s rounding
//! rules (six significant digits, `%g` outside `[0.0001, 999999]`) and integers and floats are
//! parsed with HTSlib's fast parsers (`hts_str2int`, `hts_str2uint`, `hts_str2dbl`), including
//! their quirks (`-` alone parses as 0, for example).
//!
//! Parsers take a byte slice and a start index and return the parsed value with the index just
//! past it; bytes past the end of the slice read as NUL, like the terminator of a C string.

/// Returns the byte at `i`, or NUL past the end (a C string terminator).
#[inline]
pub(crate) fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// Two-digit lookup table for integer formatting.
const DIGITS2: &[u8; 200] = b"\
0001020304050607080910111213141516171819\
2021222324252627282930313233343536373839\
4041424344454647484950515253545556575859\
6061626364656667686970717273747576777879\
8081828384858687888990919293949596979899";

/// Appends the decimal representation of `v` (HTSlib's `kputw`/`kputll`).
#[inline]
pub fn push_int(out: &mut Vec<u8>, v: i64) {
    if v < 0 {
        out.push(b'-');
    }
    push_uint(out, v.unsigned_abs());
}

/// Appends the decimal representation of `v` (HTSlib's `kputuw`/`kputull`).
#[inline(always)]
pub fn push_uint(out: &mut Vec<u8>, v: u64) {
    // Most values in VCF (depths, counts, qualities) are below 100: one or two bytes.
    if v < 10 {
        out.push(b'0' + v as u8);
    } else if v < 100 {
        let d = v as usize * 2;
        out.extend_from_slice(&DIGITS2[d..d + 2]);
    } else {
        push_uint_long(out, v);
    }
}

/// [`push_uint`] for values from 100 on. Values that fit 32 bits (positions, most others) are
/// divided in 32 bits, which is cheaper.
#[inline(never)]
fn push_uint_long(out: &mut Vec<u8>, v: u64) {
    let mut buf = [0u8; 20];
    let i = match u32::try_from(v) {
        Ok(v) => digits_u32(&mut buf, v),
        Err(_) => digits_u64(&mut buf, v),
    };
    out.extend_from_slice(&buf[i..]);
}

/// Writes the digits of `v` at the end of `buf`, returning where they start.
#[inline]
fn digits_u32(buf: &mut [u8], mut v: u32) -> usize {
    let mut i = buf.len();
    while v >= 100 {
        let d = (v % 100) as usize * 2;
        v /= 100;
        i -= 2;
        buf[i..i + 2].copy_from_slice(&DIGITS2[d..d + 2]);
    }
    if v >= 10 {
        let d = v as usize * 2;
        i -= 2;
        buf[i..i + 2].copy_from_slice(&DIGITS2[d..d + 2]);
    } else {
        i -= 1;
        buf[i] = b'0' + v as u8;
    }
    i
}

/// [`digits_u32`] for values past 32 bits.
#[cold]
fn digits_u64(buf: &mut [u8; 20], mut v: u64) -> usize {
    let mut i = buf.len();
    while v > u64::from(u32::MAX) {
        let d = (v % 100) as usize * 2;
        v /= 100;
        i -= 2;
        buf[i..i + 2].copy_from_slice(&DIGITS2[d..d + 2]);
    }
    let mut tail = [0u8; 20];
    let j = digits_u32(&mut tail, v as u32);
    let n = tail.len() - j;
    i -= n;
    buf[i..i + n].copy_from_slice(&tail[j..]);
    i
}

/// Appends a float as HTSlib's `kputd` does: up to six significant digits without trailing zeros
/// for magnitudes in `[0.0001, 999999]`, and C's `%g` otherwise.
pub fn push_double(out: &mut Vec<u8>, mut d: f64) {
    if d == 0.0 {
        if d.is_sign_negative() {
            out.extend_from_slice(b"-0");
        } else {
            out.push(b'0');
        }
        return;
    }
    if d < 0.0 {
        out.push(b'-');
        d = -d;
    }
    if !(0.0001..=999_999.0).contains(&d) {
        push_g(out, d);
        return;
    }

    // Correction for rounding, as kputd does: scale to a six-digit integer.
    let mut buf = [0u8; 21];
    let mut cp: usize = 20;
    let rint = |x: f64| x.round_ties_even() as u32;
    let mut i: u32 = if d < 0.001 {
        cp -= 1;
        rint(d * 1_000_000_000.0)
    } else if d < 0.01 {
        cp -= 2;
        rint(d * 100_000_000.0)
    } else if d < 0.1 {
        cp -= 3;
        rint(d * 10_000_000.0)
    } else if d < 1.0 {
        cp -= 4;
        rint(d * 1_000_000.0)
    } else if d < 10.0 {
        cp -= 5;
        rint(d * 100_000.0)
    } else if d < 100.0 {
        cp -= 6;
        rint(d * 10_000.0)
    } else if d < 1000.0 {
        cp -= 7;
        rint(d * 1000.0)
    } else if d < 10000.0 {
        cp -= 8;
        rint(d * 100.0)
    } else if d < 100_000.0 {
        cp -= 9;
        rint(d * 10.0)
    } else {
        cp -= 10;
        rint(d)
    };

    for k in 0..3 {
        let pair = (i % 100) as usize * 2;
        cp -= 2;
        buf[cp..cp + 2].copy_from_slice(&DIGITS2[pair..pair + 2]);
        if k < 2 {
            i /= 100;
        }
    }
    // Except when it rounds up (d=0.009999999 is i=1000000).
    if i >= 100 {
        cp -= 1;
        buf[cp] = b'0' + (i / 100) as u8;
    }

    let mut p = 20 - cp;
    let mut ep;
    if p <= 10 {
        // 0.00123 is 123, so add leading zeros and 0.
        ep = cp + 5;
        while p < 10 {
            cp -= 1;
            buf[cp] = b'0';
            p += 1;
        }
        cp -= 1;
        buf[cp] = b'.';
        cp -= 1;
        buf[cp] = b'0';
    } else {
        // 123.001 is 123001 with p==13, so move 123 down and add ".".
        cp -= 1;
        let mut xp = cp;
        ep = cp + 6;
        while p > 10 {
            buf[xp] = buf[xp + 1];
            xp += 1;
            p -= 1;
        }
        buf[xp] = b'.';
    }

    // Cull trailing zeros.
    while buf[ep] == b'0' && ep > cp {
        ep -= 1;
    }
    // End can be 1 out due to the mostly-6 but occasionally 7 (i==1) case, and "123." is "123".
    if buf[ep] != 0 && buf[ep] != b'.' {
        ep += 1;
    }
    out.extend_from_slice(&buf[cp..ep]);
}

/// Appends a non-negative float (or NaN) as C's `printf("%g")` does.
fn push_g(out: &mut Vec<u8>, d: f64) {
    use std::io::Write;

    if d.is_nan() {
        out.extend_from_slice(if d.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
        return;
    }
    if d.is_infinite() {
        out.extend_from_slice(b"inf");
        return;
    }

    // Six significant digits, rounded exactly (half to even), as glibc does.
    let mut sci = Vec::with_capacity(16);
    write!(sci, "{d:.5e}").expect("writing to a Vec cannot fail");
    let epos = sci
        .iter()
        .position(|&c| c == b'e')
        .expect("scientific formatting has an exponent");
    let exponent: i32 = std::str::from_utf8(&sci[epos + 1..])
        .ok()
        .and_then(|s| s.parse().ok())
        .expect("scientific formatting has an integer exponent");
    let digits: Vec<u8> = sci[..epos].iter().copied().filter(|c| *c != b'.').collect();

    if !(-4..6).contains(&exponent) {
        let mut mantissa = Vec::with_capacity(8);
        mantissa.push(digits[0]);
        mantissa.push(b'.');
        mantissa.extend_from_slice(&digits[1..]);
        strip_fraction_zeros(&mut mantissa);
        out.extend_from_slice(&mantissa);
        out.push(b'e');
        out.push(if exponent < 0 { b'-' } else { b'+' });
        let e = exponent.unsigned_abs();
        if e < 10 {
            out.push(b'0');
        }
        push_uint(out, u64::from(e));
    } else {
        let mut fixed = Vec::with_capacity(16);
        if exponent < 0 {
            fixed.extend_from_slice(b"0.");
            fixed.extend(std::iter::repeat_n(b'0', (-exponent - 1) as usize));
            fixed.extend_from_slice(&digits);
        } else {
            let int_digits = exponent as usize + 1;
            fixed.extend_from_slice(&digits[..int_digits]);
            fixed.push(b'.');
            fixed.extend_from_slice(&digits[int_digits..]);
        }
        strip_fraction_zeros(&mut fixed);
        out.extend_from_slice(&fixed);
    }
}

/// Removes trailing zeros after a decimal point, and the point itself if nothing follows it.
fn strip_fraction_zeros(s: &mut Vec<u8>) {
    if !s.contains(&b'.') {
        return;
    }
    while s.last() == Some(&b'0') {
        s.pop();
    }
    if s.last() == Some(&b'.') {
        s.pop();
    }
}

/// Parses a signed decimal integer as HTSlib's `hts_str2int` does, limited to `bits` bits.
///
/// Returns the value, the index after the parsed text, and whether the value overflowed (in
/// which case it is clamped). A lone sign parses as 0 and is consumed.
pub fn str2int(s: &[u8], start: usize, bits: u32) -> (i64, usize, bool) {
    let mut limit: u64 = (1u64 << (bits - 1)) - 1;
    let mut fast: u32 = (bits - 1) * 1000 / 3322 + 1;
    let mut n: u64 = 0;
    let mut v = start;
    let mut failed = false;

    let neg = match at(s, v) {
        b'-' => {
            limit += 1;
            v += 1;
            true
        }
        b'+' => {
            v += 1;
            false
        }
        _ => false,
    };
    loop {
        fast -= 1;
        let c = at(s, v);
        if fast == 0 || !c.is_ascii_digit() {
            break;
        }
        n = n.wrapping_mul(10).wrapping_add(u64::from(c - b'0'));
        v += 1;
    }
    if at(s, v) >= b'0' && fast == 0 {
        let limit_d_10 = limit / 10;
        let limit_m_10 = limit - 10 * limit_d_10;
        loop {
            let d = at(s, v).wrapping_sub(b'0');
            if d >= 10 {
                break;
            }
            let d = u64::from(d);
            if n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10) {
                n = n * 10 + d;
                v += 1;
            } else {
                loop {
                    v += 1;
                    if at(s, v).wrapping_sub(b'0') >= 10 {
                        break;
                    }
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    let value = if neg {
        (n as i64).wrapping_neg()
    } else {
        n as i64
    };
    (value, v, failed)
}

/// Parses an unsigned decimal integer as HTSlib's `hts_str2uint` does, limited to `bits` bits.
///
/// Returns the value, the index after the parsed text, and whether the value overflowed.
pub fn str2uint(s: &[u8], start: usize, bits: u32) -> (u64, usize, bool) {
    let limit: u64 = if bits < 64 {
        (1u64 << bits) - 1
    } else {
        u64::MAX
    };
    let mut fast: u32 = bits * 1000 / 3322 + 1;
    let mut n: u64 = 0;
    let mut v = start;
    let mut failed = false;
    if at(s, v) == b'+' {
        v += 1;
    }
    loop {
        fast -= 1;
        let c = at(s, v);
        if fast == 0 || !c.is_ascii_digit() {
            break;
        }
        n = n.wrapping_mul(10).wrapping_add(u64::from(c - b'0'));
        v += 1;
    }
    if at(s, v).wrapping_sub(b'0') < 10 && fast == 0 {
        let limit_d_10 = limit / 10;
        let limit_m_10 = limit - 10 * limit_d_10;
        loop {
            let d = at(s, v).wrapping_sub(b'0');
            if d >= 10 {
                break;
            }
            let d = u64::from(d);
            if n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10) {
                n = n * 10 + d;
                v += 1;
            } else {
                loop {
                    v += 1;
                    if at(s, v).wrapping_sub(b'0') >= 10 {
                        break;
                    }
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    (n, v, failed)
}

/// Powers of ten for [`str2dbl`]'s fast path (HTSlib's table, which starts with 1, 1).
const POW10: [f64; 22] = [
    1.0, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
    1e16, 1e17, 1e18, 1e19, 1e20,
];

/// Returns whether `c` is C's `isspace`.
#[inline]
pub(crate) fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Parses a float as HTSlib's `hts_str2dbl` does: a fast path for plain decimals of up to 14
/// digits, and C's `strtod` for anything else (exponents, `inf`, `nan`, hex).
///
/// Returns the value, the index after the parsed text, and whether parsing failed (nothing was
/// consumed).
pub fn str2dbl(s: &[u8], start: usize) -> (f64, usize, bool) {
    let mut v = start;
    while is_space(at(s, v)) {
        v += 1;
    }
    let mut neg = false;
    match at(s, v) {
        b'-' => {
            neg = true;
            v += 1;
        }
        b'+' => v += 1,
        _ => {}
    }
    match at(s, v) {
        b'1'..=b'9' => {}
        b'0' if !matches!(at(s, v + 1), b'x' | b'X') => {}
        _ => {
            let (d, end) = strtod(s, start);
            return (d, end, end == start);
        }
    }
    while at(s, v) == b'0' {
        v += 1;
    }
    let digits_start = v;
    let mut max_len: i32 = 15;
    let mut n: u64 = 0;
    let mut point: isize = -1;
    loop {
        max_len -= 1;
        let c = at(s, v);
        if max_len == 0 || !c.is_ascii_digit() {
            break;
        }
        n = n * 10 + u64::from(c - b'0');
        v += 1;
    }
    if max_len != 0 && at(s, v) == b'.' {
        point = (v - digits_start) as isize;
        v += 1;
        loop {
            max_len -= 1;
            let c = at(s, v);
            if max_len == 0 || !c.is_ascii_digit() {
                break;
            }
            n = n * 10 + u64::from(c - b'0');
            v += 1;
        }
    }
    if point < 0 {
        point = (v - digits_start) as isize;
    }
    if max_len == 0 || matches!(at(s, v), b'e' | b'E') {
        let (d, end) = strtod(s, start);
        return (d, end, end == start);
    }
    let d = n as f64 / POW10[((v - digits_start) as isize - point) as usize];
    (if neg { -d } else { d }, v, false)
}

/// C's `atof`: [`strtod`] without the end position (0 when nothing parses).
pub fn atof(s: &[u8]) -> f64 {
    strtod(s, 0).0
}

/// Parses a float as glibc's `strtod` does, returning the value and the index after the longest
/// valid prefix (`start` when nothing parses, with a value of 0).
pub fn strtod(s: &[u8], start: usize) -> (f64, usize) {
    let mut i = start;
    while is_space(at(s, i)) {
        i += 1;
    }
    let sign_pos = i;
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
    let apply = |d: f64| if neg { -d } else { d };
    let lower = |k: usize| at(s, k).to_ascii_lowercase();

    // inf, infinity, nan, nan(chars)
    if lower(i) == b'i' && lower(i + 1) == b'n' && lower(i + 2) == b'f' {
        let mut end = i + 3;
        if (0..5).all(|k| lower(end + k) == b"inity"[k]) {
            end += 5;
        }
        return (apply(f64::INFINITY), end);
    }
    if lower(i) == b'n' && lower(i + 1) == b'a' && lower(i + 2) == b'n' {
        let mut end = i + 3;
        if at(s, end) == b'(' {
            let mut k = end + 1;
            while at(s, k).is_ascii_alphanumeric() || at(s, k) == b'_' {
                k += 1;
            }
            if at(s, k) == b')' {
                end = k + 1;
            }
        }
        return (apply(f64::NAN), end);
    }
    if at(s, i) == b'0' && lower(i + 1) == b'x' {
        if let Some((d, end)) = parse_hex_float(s, i + 2) {
            return (apply(d), end);
        }
        // "0x" without hex digits parses as 0, consuming the "0".
        return (apply(0.0), i + 1);
    }

    let mant_start = i;
    let mut ndigits = 0;
    while at(s, i).is_ascii_digit() {
        i += 1;
        ndigits += 1;
    }
    if at(s, i) == b'.' {
        i += 1;
        while at(s, i).is_ascii_digit() {
            i += 1;
            ndigits += 1;
        }
    }
    if ndigits == 0 {
        let _ = sign_pos;
        return (0.0, start);
    }
    let mant_end = i;
    let mut end = i;
    if matches!(at(s, i), b'e' | b'E') {
        let mut k = i + 1;
        if matches!(at(s, k), b'+' | b'-') {
            k += 1;
        }
        if at(s, k).is_ascii_digit() {
            while at(s, k).is_ascii_digit() {
                k += 1;
            }
            end = k;
        }
    }
    let text = &s[mant_start..end.min(s.len())];
    let parsed = std::str::from_utf8(text)
        .ok()
        .and_then(|t| t.parse::<f64>().ok())
        .or_else(|| {
            std::str::from_utf8(&s[mant_start..mant_end])
                .ok()
                .and_then(|t| t.parse::<f64>().ok())
        })
        .unwrap_or(0.0);
    (apply(parsed), end)
}

/// Parses the digits of a C99 hex float after `0x`, returning the value and the end index.
fn parse_hex_float(s: &[u8], start: usize) -> Option<(f64, usize)> {
    let mut i = start;
    let mut mantissa: u64 = 0;
    let mut exp: i64 = 0;
    let mut ndigits = 0;
    let hexval = |c: u8| (c as char).to_digit(16);
    while let Some(h) = hexval(at(s, i)) {
        if mantissa >> 60 == 0 {
            mantissa = mantissa * 16 + u64::from(h);
        } else {
            exp += 4;
        }
        i += 1;
        ndigits += 1;
    }
    if at(s, i) == b'.' {
        i += 1;
        while let Some(h) = hexval(at(s, i)) {
            if mantissa >> 60 == 0 {
                mantissa = mantissa * 16 + u64::from(h);
                exp -= 4;
            }
            i += 1;
            ndigits += 1;
        }
    }
    if ndigits == 0 {
        return None;
    }
    if matches!(at(s, i), b'p' | b'P') {
        let mut k = i + 1;
        let eneg = match at(s, k) {
            b'-' => {
                k += 1;
                true
            }
            b'+' => {
                k += 1;
                false
            }
            _ => false,
        };
        if at(s, k).is_ascii_digit() {
            let mut e: i64 = 0;
            while at(s, k).is_ascii_digit() {
                e = (e * 10 + i64::from(at(s, k) - b'0')).min(100_000);
                k += 1;
            }
            exp += if eneg { -e } else { e };
            i = k;
        }
    }
    let exp = i32::try_from(exp.clamp(-100_000, 100_000)).unwrap_or(0);
    Some(((mantissa as f64) * 2f64.powi(exp), i))
}

/// Flags for [`parse_decimal`].
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ParseDecimalFlags {
    /// Accept commas between digits (HTSlib's `HTS_PARSE_THOUSANDS_SEP`).
    pub thousands_separator: bool,
}

/// Parses a genomic coordinate as HTSlib's `hts_parse_decimal` does: an optional sign, digits
/// (optionally with comma separators), an optional fraction, and an `e` exponent or a `k`/`m`/`g`
/// suffix.
///
/// Returns the value and the index after the parsed text, or `start` when there are no digits.
pub fn parse_decimal(s: &[u8], start: usize, flags: ParseDecimalFlags) -> (i64, usize) {
    let mut str_start = start;
    while is_space(at(s, str_start)) {
        str_start += 1;
    }
    let mut p = str_start;
    let mut n: u64 = 0;
    let mut digits = 0;
    let mut decimals: i32 = 0;
    let mut e: i32 = 0;
    let mut lost = 0u64;
    let mut sign = b'+';
    if matches!(at(s, p), b'+' | b'-') {
        sign = at(s, p);
        p += 1;
    }
    loop {
        let c = at(s, p);
        if c.is_ascii_digit() {
            digits += 1;
            n = n.wrapping_mul(10).wrapping_add(u64::from(c - b'0'));
            p += 1;
        } else if c == b',' && flags.thousands_separator {
            p += 1;
        } else {
            break;
        }
    }
    if at(s, p) == b'.' {
        p += 1;
        while at(s, p).is_ascii_digit() {
            decimals += 1;
            digits += 1;
            n = n.wrapping_mul(10).wrapping_add(u64::from(at(s, p) - b'0'));
            p += 1;
        }
    }
    match at(s, p) {
        b'e' | b'E' => {
            p += 1;
            let mut esign = b'+';
            if matches!(at(s, p), b'+' | b'-') {
                esign = at(s, p);
                p += 1;
            }
            while at(s, p).is_ascii_digit() {
                e = e.wrapping_mul(10).wrapping_add(i32::from(at(s, p) - b'0'));
                p += 1;
            }
            if esign == b'-' {
                e = -e;
            }
        }
        b'k' | b'K' => {
            e += 3;
            p += 1;
        }
        b'm' | b'M' => {
            e += 6;
            p += 1;
        }
        b'g' | b'G' => {
            e += 9;
            p += 1;
        }
        _ => {}
    }
    e -= decimals;
    while e > 0 && n != 0 {
        n = n.wrapping_mul(10);
        e -= 1;
    }
    while e < 0 && n != 0 {
        lost += n % 10;
        n /= 10;
        e += 1;
    }
    if lost > 0 {
        crate::variant::log::warn(
            "hts_parse_decimal",
            format_args!(
                "Discarding fractional part of {}",
                String::from_utf8_lossy(&s[str_start..p.min(s.len())])
            ),
        );
    }
    let end = if digits > 0 { p } else { start };
    let value = if sign == b'+' {
        n as i64
    } else {
        (n as i64).wrapping_neg()
    };
    (value, end)
}

/// C's `strtol(s, &end, 10)` clamped to `i64`: returns the value and the end index (`start` when
/// no digits were found).
pub fn strtol(s: &[u8], start: usize) -> (i64, usize) {
    let mut i = start;
    while is_space(at(s, i)) {
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
    if !at(s, i).is_ascii_digit() {
        return (0, start);
    }
    let mut n: i128 = 0;
    while at(s, i).is_ascii_digit() {
        n = (n * 10 + i128::from(at(s, i) - b'0')).min(i128::from(i64::MAX) + 1);
        i += 1;
    }
    let v = if neg { -n } else { n };
    (
        v.clamp(i128::from(i64::MIN), i128::from(i64::MAX)) as i64,
        i,
    )
}

/// C's `atoi`/`atol`: [`strtol`] without the end position.
pub fn atoi(s: &[u8]) -> i64 {
    strtol(s, 0).0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fmt(d: f64) -> String {
        let mut v = Vec::new();
        push_double(&mut v, d);
        String::from_utf8(v).unwrap()
    }

    #[test]
    fn kputd_matches_htslib() {
        for (d, want) in [
            (0.0, "0"),
            (-0.0, "-0"),
            (1.5, "1.5"),
            (0.5, "0.5"),
            (0.001234, "0.001234"),
            (0.0099999999, "0.01"),
            (123.456789, "123.457"),
            (1.0, "1"),
            (100.0, "100"),
            (999999.0, "999999"),
            (1234565.0, "1.23456e+06"),
            (0.00001, "1e-05"),
            (0.000099999995, "0.0001"),
            (1e300, "1e+300"),
            (-2.5, "-2.5"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (f64::NAN, "nan"),
            (f64::from(0.1f32), "0.1"),
            (f64::from(1.23457f32), "1.23457"),
        ] {
            assert_eq!(fmt(d), want, "formatting {d:e}");
        }
    }

    #[test]
    fn ints() {
        let mut v = Vec::new();
        for x in [0i64, 7, 10, 99, 100, 12345, -1, -128, i64::MIN, i64::MAX] {
            v.clear();
            push_int(&mut v, x);
            assert_eq!(v, x.to_string().as_bytes());
        }
        let edge = u64::from(u32::MAX);
        for x in [
            edge - 1,
            edge,
            edge + 1,
            edge * 10,
            10_000_000_000,
            4_294_967_300,
            u64::MAX,
        ] {
            v.clear();
            push_uint(&mut v, x);
            assert_eq!(v, x.to_string().as_bytes());
        }
    }

    #[test]
    fn str2int_quirks() {
        assert_eq!(str2int(b"-", 0, 64), (0, 1, false));
        assert_eq!(str2int(b"12,3", 0, 64), (12, 2, false));
        assert_eq!(str2int(b"x", 0, 64), (0, 0, false));
        let (v, end, failed) = str2int(b"99999999999999999999999", 0, 64);
        assert_eq!((v, end, failed), (i64::MAX, 23, true));
        assert_eq!(str2int(b"-9223372036854775808", 0, 64).0, i64::MIN);
    }

    #[test]
    fn str2dbl_paths() {
        assert_eq!(str2dbl(b"1.5", 0), (1.5, 3, false));
        assert_eq!(str2dbl(b"0.05,", 0), (0.05, 4, false));
        assert_eq!(str2dbl(b"1e3", 0), (1000.0, 3, false));
        assert_eq!(str2dbl(b"abc", 0), (0.0, 0, true));
        let (d, end, failed) = str2dbl(b"nan", 0);
        assert!(d.is_nan() && end == 3 && !failed);
        assert_eq!(str2dbl(b".5", 0), (0.5, 2, false));
    }

    #[test]
    fn parse_decimal_suffixes() {
        let f = ParseDecimalFlags::default();
        assert_eq!(parse_decimal(b"1k", 0, f), (1000, 2));
        assert_eq!(parse_decimal(b"1.5M", 0, f), (1_500_000, 4));
        assert_eq!(parse_decimal(b"1,000", 0, f), (1, 1));
        let t = ParseDecimalFlags {
            thousands_separator: true,
        };
        assert_eq!(parse_decimal(b"1,000-", 0, t), (1000, 5));
        assert_eq!(parse_decimal(b"-", 0, f), (0, 0));
    }
}
