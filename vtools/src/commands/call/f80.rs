//! The x87 80-bit `long double` arithmetic that bcftools' `prob1.c` sums in.
//!
//! On x86-64 Linux, C's `long double` keeps a 64-bit significand, and `prob1.c` accumulates
//! allele-frequency sums in it before rounding the results to `double`. Summing in `f64`
//! instead changes the last bits of QUAL, AF1, and FQ now and then. On x86-64 the operations
//! run on the x87 unit itself, as the C code does; elsewhere a software emulation covers the
//! few operations used: conversion from `f64`, addition of non-negative values,
//! multiplication, division, and rounding back to `f64`, all rounding to nearest with ties
//! to even as the x87 unit does.

#[cfg(target_arch = "x86_64")]
pub use self::x87::F80;

#[cfg(not(target_arch = "x86_64"))]
pub use self::soft::F80;

/// The extended format on the x87 unit.
#[cfg(target_arch = "x86_64")]
mod x87 {
    use std::arch::asm;

    /// An 80-bit extended-precision value as the x87 unit stores it (`fstp tbyte`), in 16
    /// bytes with zero padding.
    #[derive(Clone, Copy, Debug, PartialEq)]
    #[repr(C, align(16))]
    pub struct F80([u8; 16]);

    /// One half, for `x / 2.` (exact in the extended format).
    static HALF: f64 = 0.5;

    // SAFETY (all asm blocks below): every load reads 8 bytes from a live `f64` (or from
    // `len` elements of a slice) or 10 bytes from an `F80`'s 16-byte buffer, and every store
    // writes 8 bytes into a live `f64` or 10 bytes into an `F80`'s buffer. The x87 register
    // stack is empty on entry and is left empty (each push is matched by a pop), and all
    // eight registers are declared clobbered. The control word is the process default
    // (64-bit significand, round to nearest), which C's long double arithmetic uses too.
    impl F80 {
        /// Converts exactly from a double (C's `(long double)x`).
        #[cfg(test)]
        #[must_use]
        pub fn from_f64(x: f64) -> Self {
            let mut out = Self([0; 16]);
            // SAFETY: see the impl block.
            unsafe {
                asm!(
                    "fld qword ptr [{x}]",
                    "fstp tbyte ptr [{o}]",
                    x = in(reg) &raw const x,
                    o = in(reg) out.0.as_mut_ptr(),
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }

        /// Converts from an integer (C's `(long double)k`).
        #[must_use]
        #[inline]
        pub fn from_i64(k: i64) -> Self {
            let mut out = Self([0; 16]);
            // SAFETY: see the impl block.
            unsafe {
                asm!(
                    "fild qword ptr [{k}]",
                    "fstp tbyte ptr [{o}]",
                    k = in(reg) &raw const k,
                    o = in(reg) out.0.as_mut_ptr(),
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }

        /// `sum += (long double)a[i] * b[i]` over the slices, from zero.
        ///
        /// # Panics
        ///
        /// Panics when the slices differ in length.
        #[must_use]
        pub fn dot(a: &[f64], b: &[f64]) -> Self {
            Self::dot_scaled(a, b, 1.0)
        }

        /// `sum += (long double)a[i] / 2. * b[i]` over the slices, from zero.
        ///
        /// # Panics
        ///
        /// Panics when the slices differ in length.
        #[must_use]
        pub fn dot_half(a: &[f64], b: &[f64]) -> Self {
            Self::dot_scaled(a, b, HALF)
        }

        /// The sum of `(long double)a[i] * scale * b[i]`, with `scale` 1 or 1/2 (both exact).
        fn dot_scaled(a: &[f64], b: &[f64], scale: f64) -> Self {
            assert_eq!(
                a.len(),
                b.len(),
                "dot product of slices of different lengths"
            );
            let mut out = Self([0; 16]);
            // SAFETY: see the impl block; the loop reads `a.len()` elements of each slice.
            unsafe {
                asm!(
                    "fldz",
                    "test {n}, {n}",
                    "jz 3f",
                    "2:",
                    "fld qword ptr [{a}]",
                    "fmul qword ptr [{s}]",
                    "fmul qword ptr [{b}]",
                    "faddp st(1), st",
                    "add {a}, 8",
                    "add {b}, 8",
                    "dec {n}",
                    "jnz 2b",
                    "3:",
                    "fstp tbyte ptr [{o}]",
                    a = inout(reg) a.as_ptr() => _,
                    b = inout(reg) b.as_ptr() => _,
                    n = inout(reg) a.len() => _,
                    s = in(reg) &raw const scale,
                    o = in(reg) out.0.as_mut_ptr(),
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }

        /// `sum += (long double)a[i]` over the slice, from zero.
        #[must_use]
        pub fn sum(a: &[f64]) -> Self {
            let mut out = Self([0; 16]);
            // SAFETY: see the impl block; the loop reads `a.len()` elements.
            unsafe {
                asm!(
                    "fldz",
                    "test {n}, {n}",
                    "jz 3f",
                    "2:",
                    "fadd qword ptr [{a}]",
                    "add {a}, 8",
                    "dec {n}",
                    "jnz 2b",
                    "3:",
                    "fstp tbyte ptr [{o}]",
                    a = inout(reg) a.as_ptr() => _,
                    n = inout(reg) a.len() => _,
                    o = in(reg) out.0.as_mut_ptr(),
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }

        /// `self / den` in extended precision, rounded to a double.
        #[must_use]
        #[inline]
        pub fn ratio(self, den: Self) -> f64 {
            let mut out = 0.0f64;
            // SAFETY: see the impl block.
            unsafe {
                asm!(
                    "fld tbyte ptr [{a}]",
                    "fld tbyte ptr [{b}]",
                    "fdivp st(1), st",
                    "fstp qword ptr [{o}]",
                    a = in(reg) self.0.as_ptr(),
                    b = in(reg) den.0.as_ptr(),
                    o = in(reg) &raw mut out,
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }

        /// `(long double)x / den`, rounded to a double.
        #[must_use]
        #[inline]
        pub fn div_into(x: f64, den: Self) -> f64 {
            let mut out = 0.0f64;
            // SAFETY: see the impl block.
            unsafe {
                asm!(
                    "fld qword ptr [{x}]",
                    "fld tbyte ptr [{b}]",
                    "fdivp st(1), st",
                    "fstp qword ptr [{o}]",
                    x = in(reg) &raw const x,
                    b = in(reg) den.0.as_ptr(),
                    o = in(reg) &raw mut out,
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }

        /// The stored bytes: the significand, then the sign and exponent.
        #[cfg(test)]
        pub fn to_bytes(self) -> [u8; 16] {
            self.0
        }

        /// Rounds to a double (C's conversion of `long double` to `double`).
        #[must_use]
        #[inline]
        pub fn to_f64(self) -> f64 {
            let mut out = 0.0f64;
            // SAFETY: see the impl block.
            unsafe {
                asm!(
                    "fld tbyte ptr [{a}]",
                    "fstp qword ptr [{o}]",
                    a = in(reg) self.0.as_ptr(),
                    o = in(reg) &raw mut out,
                    out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
                    out("st(4)") _, out("st(5)") _, out("st(6)") _, out("st(7)") _,
                    options(nostack),
                );
            }
            out
        }
    }
}

/// The software emulation, used where there is no x87 unit (and in tests, against it).
#[cfg_attr(target_arch = "x86_64", cfg(test))]
mod soft {
    /// A non-negative extended-precision value, or a value outside the emulated range (negative,
    /// infinite, or NaN) carried as an `f64`.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub enum F80 {
        /// `mant * 2^exp`, with the top bit of `mant` set, or zero when `mant` is 0.
        Finite { mant: u64, exp: i32 },
        /// A value the emulation does not cover; arithmetic falls back to `f64`.
        Special(f64),
    }

    impl F80 {
        /// Zero.
        pub const ZERO: Self = Self::Finite { mant: 0, exp: 0 };

        /// `sum += (long double)a[i] * b[i]` over the slices, from zero.
        #[must_use]
        pub fn dot(a: &[f64], b: &[f64]) -> Self {
            assert_eq!(
                a.len(),
                b.len(),
                "dot product of slices of different lengths"
            );
            a.iter().zip(b).fold(Self::ZERO, |sum, (&x, &y)| {
                sum.add(Self::from_f64(x).mul(Self::from_f64(y)))
            })
        }

        /// `sum += (long double)a[i] / 2. * b[i]` over the slices, from zero.
        #[must_use]
        pub fn dot_half(a: &[f64], b: &[f64]) -> Self {
            assert_eq!(
                a.len(),
                b.len(),
                "dot product of slices of different lengths"
            );
            a.iter().zip(b).fold(Self::ZERO, |sum, (&x, &y)| {
                sum.add(Self::from_f64(x).half().mul(Self::from_f64(y)))
            })
        }

        /// `sum += (long double)a[i]` over the slice, from zero.
        #[must_use]
        pub fn sum(a: &[f64]) -> Self {
            a.iter()
                .fold(Self::ZERO, |sum, &x| sum.add(Self::from_f64(x)))
        }

        /// `self / den` in extended precision, rounded to a double.
        #[must_use]
        pub fn ratio(self, den: Self) -> f64 {
            self.div(den).to_f64()
        }

        /// `(long double)x / den`, rounded to a double.
        #[must_use]
        pub fn div_into(x: f64, den: Self) -> f64 {
            Self::from_f64(x).ratio(den)
        }

        /// Converts exactly from a double (C's `(long double)x`).
        #[must_use]
        pub fn from_f64(x: f64) -> Self {
            if x == 0.0 {
                return Self::ZERO;
            }
            if !x.is_finite() || x < 0.0 {
                return Self::Special(x);
            }
            let bits = x.to_bits();
            let e = ((bits >> 52) & 0x7ff) as i32;
            let frac = bits & ((1u64 << 52) - 1);
            if e == 0 {
                let lz = frac.leading_zeros() as i32;
                Self::Finite {
                    mant: frac << lz,
                    exp: -1074 - lz,
                }
            } else {
                Self::Finite {
                    mant: (frac | (1u64 << 52)) << 11,
                    exp: e - 1075 - 11,
                }
            }
        }

        /// Converts from an integer (C's `(long double)k`).
        #[must_use]
        pub fn from_i64(k: i64) -> Self {
            Self::from_f64(k as f64)
        }

        /// Rounds `m * 2^exp` (plus a nonzero tail below `m` when `sticky`) to 64 bits.
        fn round(m: u128, exp: i32, sticky: bool) -> Self {
            if m == 0 {
                return Self::ZERO;
            }
            let bits = 128 - m.leading_zeros() as i32;
            if bits <= 64 {
                debug_assert!(!sticky, "a tail needs guard bits");
                let shift = 64 - bits;
                return Self::Finite {
                    mant: (m << shift) as u64,
                    exp: exp - shift,
                };
            }
            let shift = bits - 64;
            let rem = m & ((1u128 << shift) - 1);
            let half = 1u128 << (shift - 1);
            let mut mant = m >> shift;
            let mut exp = exp + shift;
            if rem > half || (rem == half && (sticky || mant & 1 == 1)) {
                mant += 1;
                if mant == 1u128 << 64 {
                    mant >>= 1;
                    exp += 1;
                }
            }
            Self::Finite {
                mant: mant as u64,
                exp,
            }
        }

        /// Adds two values (C's `+` on `long double`).
        #[must_use]
        pub fn add(self, other: Self) -> Self {
            let (Self::Finite { mant: ma, exp: ea }, Self::Finite { mant: mb, exp: eb }) =
                (self, other)
            else {
                return Self::Special(self.to_f64() + other.to_f64());
            };
            if ma == 0 {
                return other;
            }
            if mb == 0 {
                return self;
            }
            let ((mbig, ebig), (msmall, esmall)) = if (ea, ma) >= (eb, mb) {
                ((ma, ea), (mb, eb))
            } else {
                ((mb, eb), (ma, ea))
            };
            // 63 guard bits keep the sum below 2^128.
            let big = u128::from(mbig) << 63;
            let small = u128::from(msmall) << 63;
            let d = (ebig - esmall) as u32;
            let (shifted, sticky) = if d >= 127 {
                (0, true)
            } else {
                (small >> d, small & ((1u128 << d) - 1) != 0)
            };
            Self::round(big + shifted, ebig - 63, sticky)
        }

        /// Halves a value (C's `x / 2.` on `long double`, which is exact in the extended range).
        #[must_use]
        pub fn half(self) -> Self {
            match self {
                Self::Finite { mant: 0, .. } => self,
                Self::Finite { mant, exp } => Self::Finite { mant, exp: exp - 1 },
                Self::Special(x) => Self::Special(x / 2.0),
            }
        }

        /// Multiplies two values (C's `*` on `long double`).
        #[must_use]
        pub fn mul(self, other: Self) -> Self {
            let (Self::Finite { mant: ma, exp: ea }, Self::Finite { mant: mb, exp: eb }) =
                (self, other)
            else {
                return Self::Special(self.to_f64() * other.to_f64());
            };
            if ma == 0 || mb == 0 {
                return Self::ZERO;
            }
            Self::round(u128::from(ma) * u128::from(mb), ea + eb, false)
        }

        /// Divides two values (C's `/` on `long double`).
        #[must_use]
        pub fn div(self, other: Self) -> Self {
            let (Self::Finite { mant: ma, exp: ea }, Self::Finite { mant: mb, exp: eb }) =
                (self, other)
            else {
                return Self::Special(self.to_f64() / other.to_f64());
            };
            if mb == 0 {
                return Self::Special(if ma == 0 { f64::NAN } else { f64::INFINITY });
            }
            if ma == 0 {
                return Self::ZERO;
            }
            let n = u128::from(ma) << 64;
            let d = u128::from(mb);
            let q = n / d;
            let r = n % d;
            let exp = ea - 64 - eb;
            if q >> 64 != 0 {
                Self::round(q, exp, r != 0)
            } else {
                // One more quotient bit, so that rounding sees a guard bit.
                let r2 = r << 1;
                let (q2, r2) = if r2 >= d {
                    ((q << 1) | 1, r2 - d)
                } else {
                    (q << 1, r2)
                };
                Self::round(q2, exp - 1, r2 != 0)
            }
        }

        /// Rounds to a double (C's conversion of `long double` to `double`).
        #[must_use]
        pub fn to_f64(self) -> f64 {
            let (mant, exp) = match self {
                Self::Special(x) => return x,
                Self::Finite { mant: 0, .. } => return 0.0,
                Self::Finite { mant, exp } => (mant, exp),
            };
            // The value is in [2^e, 2^(e+1)).
            let e = exp + 63;
            if e > 1023 {
                return f64::INFINITY;
            }
            let m = u128::from(mant);
            // Bits to drop: 11 for a normal double, more for a subnormal one.
            let shift = if e >= -1022 { 11 } else { -(exp + 1074) };
            if shift >= 128 {
                return 0.0;
            }
            let shift = shift as u32;
            let rem = m & ((1u128 << shift) - 1);
            let half = 1u128 << (shift - 1);
            let mut r = m >> shift;
            if rem > half || (rem == half && r & 1 == 1) {
                r += 1;
            }
            if e >= -1022 {
                let mut e = e;
                if r == 1u128 << 53 {
                    r >>= 1;
                    e += 1;
                    if e > 1023 {
                        return f64::INFINITY;
                    }
                }
                f64::from_bits((((e + 1023) as u64) << 52) | (r as u64 & ((1u64 << 52) - 1)))
            } else {
                // A subnormal (or, after rounding up, the smallest normal) double.
                f64::from_bits(r as u64)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_doubles() {
        for x in [1.0, 0.1, 3.5e-300, 5e-324, 1.7e308, 2.2250738585072014e-308] {
            assert_eq!(F80::from_f64(x).to_f64(), x, "{x}");
        }
    }

    #[test]
    fn keeps_extra_precision() {
        let one = F80::from_f64(1.0);
        // 1 + 2^-60 is exact in extended precision but not in double.
        let s = F80::sum(&[1.0, 2f64.powi(-60)]);
        assert_ne!(s, one);
        assert_eq!(s.to_f64(), 1.0);
        assert_eq!(s.ratio(one), 1.0);
        // 2^-53 is half an ulp of 1.0: the tie rounds to even in double, but the extra bits
        // of 1 + 2^-53 + 2^-60 round up.
        assert_eq!(F80::sum(&[1.0, 2f64.powi(-53)]).to_f64(), 1.0);
        assert_eq!(
            F80::sum(&[1.0, 2f64.powi(-53), 2f64.powi(-60)]).to_f64(),
            1.0 + f64::EPSILON
        );
    }

    #[test]
    fn divides_and_multiplies() {
        let three = F80::from_f64(3.0);
        assert_eq!(F80::from_f64(1.0).ratio(three), 1.0 / 3.0);
        assert_eq!(F80::from_i64(3), three);
        assert_eq!(F80::dot(&[1.0, 2.0], &[3.0, 4.0]).to_f64(), 11.0);
        assert_eq!(F80::dot_half(&[1.0, 2.0], &[3.0, 4.0]).to_f64(), 5.5);
        assert_eq!(F80::dot(&[], &[]), F80::sum(&[]));
        let tiny = 5e-324;
        assert_eq!(
            F80::dot_half(&[tiny], &[4.0]).ratio(F80::from_f64(2.0)),
            tiny
        );
    }

    /// The emulation matches the x87 unit on the kernels over a wide range of magnitudes.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn emulation_matches_x87() {
        use super::soft;

        // The x87 value as (significand, exponent of its lowest bit), or None for zero.
        fn parts(x: F80) -> Option<(u64, i32)> {
            let bytes = x.to_bytes();
            let mant = u64::from_le_bytes(bytes[..8].try_into().expect("8 bytes"));
            let e = i32::from(u16::from_le_bytes([bytes[8], bytes[9]]) & 0x7fff);
            (mant != 0).then_some((mant, e - 16383 - 63))
        }
        fn soft_parts(x: soft::F80) -> Option<(u64, i32)> {
            match x {
                soft::F80::Finite { mant: 0, .. } => None,
                soft::F80::Finite { mant, exp } => Some((mant, exp)),
                soft::F80::Special(v) => panic!("unexpected special value {v}"),
            }
        }
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let m = (state >> 11) as f64 / (1u64 << 53) as f64 + 0.5;
            m * 2f64.powi((state % 200) as i32 - 100)
        };
        for len in 0..2000 {
            let a: Vec<f64> = (0..len % 9).map(|_| next()).collect();
            let b: Vec<f64> = (0..len % 9).map(|_| next()).collect();
            let den = next();
            for (hw, sw) in [
                (F80::dot(&a, &b), soft::F80::dot(&a, &b)),
                (F80::dot_half(&a, &b), soft::F80::dot_half(&a, &b)),
                (F80::sum(&a), soft::F80::sum(&a)),
            ] {
                assert_eq!(parts(hw), soft_parts(sw), "{a:?} {b:?}");
                assert_eq!(hw.to_f64(), sw.to_f64(), "{a:?} {b:?}");
                assert_eq!(
                    hw.ratio(F80::from_f64(den)),
                    sw.ratio(soft::F80::from_f64(den)),
                    "{a:?} {b:?} {den}"
                );
                assert_eq!(
                    F80::from_f64(den).ratio(hw),
                    soft::F80::from_f64(den).ratio(sw),
                    "{a:?} {b:?} {den}"
                );
                assert_eq!(
                    F80::div_into(den, hw),
                    soft::F80::div_into(den, sw),
                    "{a:?} {b:?} {den}"
                );
            }
        }
        for k in [1, 2, 3, 1000, 123_456_789] {
            assert_eq!(
                parts(F80::from_i64(k)),
                soft_parts(soft::F80::from_i64(k)),
                "{k}"
            );
        }
    }
}
