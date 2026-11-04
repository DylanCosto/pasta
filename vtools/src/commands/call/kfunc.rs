//! Special functions from HTSlib's `kfunc.c` and bcftools' `kmin.c`, with the same
//! floating-point operations in the same order.

/// The C library's `lgamma` (glibc's, as bcftools links it).
pub fn lgamma(x: f64) -> f64 {
    unsafe extern "C" {
        fn lgamma_r(x: f64, sign: *mut libc::c_int) -> f64;
    }
    let mut sign: libc::c_int = 0;
    // SAFETY: `lgamma_r` reads `x` and writes the sign of Γ(x) to the valid `sign`.
    unsafe { lgamma_r(x, &raw mut sign) }
}

/// Log gamma (HTSlib's `kf_lgamma`, AS245).
pub fn kf_lgamma(z: f64) -> f64 {
    let mut x = 0.0;
    x += 0.1659470187408462e-06 / (z + 7.0);
    x += 0.9934937113930748e-05 / (z + 6.0);
    x -= 0.1385710331296526 / (z + 5.0);
    x += 12.50734324009056 / (z + 4.0);
    x -= 176.6150291498386 / (z + 3.0);
    x += 771.3234287757674 / (z + 2.0);
    x -= 1259.139216722289 / (z + 1.0);
    x += 676.5203681218835 / z;
    x += 0.9999999999995183;
    x.ln() - 5.581_061_466_795_328 - z + (z - 0.5) * (z + 6.5).ln()
}

const KF_GAMMA_EPS: f64 = 1e-14;
const KF_TINY: f64 = 1e-290;

/// Regularized lower incomplete gamma function by series expansion.
fn kf_gammap_series(s: f64, z: f64) -> f64 {
    let mut sum = 1.0;
    let mut x = 1.0;
    for k in 1..100 {
        x *= z / (s + f64::from(k));
        sum += x;
        if x / sum < KF_GAMMA_EPS {
            break;
        }
    }
    (s * z.ln() - z - kf_lgamma(s + 1.0) + sum.ln()).exp()
}

/// Regularized upper incomplete gamma function by continued fraction.
fn kf_gammaq_cf(s: f64, z: f64) -> f64 {
    let mut f = 1.0 + z - s;
    let mut c = f;
    let mut d = 0.0;
    for j in 1..100 {
        let a = f64::from(j) * (s - f64::from(j));
        let b = f64::from((j << 1) + 1) + z - s;
        d = b + a * d;
        if d < KF_TINY {
            d = KF_TINY;
        }
        c = b + a / c;
        if c < KF_TINY {
            c = KF_TINY;
        }
        d = 1.0 / d;
        let dd = c * d;
        f *= dd;
        if (dd - 1.0).abs() < KF_GAMMA_EPS {
            break;
        }
    }
    (s * z.ln() - z - kf_lgamma(s) - f.ln()).exp()
}

/// HTSlib's `kf_gammaq`.
pub fn kf_gammaq(s: f64, z: f64) -> f64 {
    if z <= 1.0 || z < s {
        1.0 - kf_gammap_series(s, z)
    } else {
        kf_gammaq_cf(s, z)
    }
}

fn kf_betai_aux(a: f64, b: f64, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    let mut f = 1.0;
    let mut c = f;
    let mut d = 0.0;
    for j in 1..200 {
        let m = f64::from(j >> 1);
        let aa = if j & 1 != 0 {
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        } else {
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        };
        d = 1.0 + aa * d;
        if d < KF_TINY {
            d = KF_TINY;
        }
        c = 1.0 + aa / c;
        if c < KF_TINY {
            c = KF_TINY;
        }
        d = 1.0 / d;
        let dd = c * d;
        f *= dd;
        if (dd - 1.0).abs() < KF_GAMMA_EPS {
            break;
        }
    }
    (kf_lgamma(a + b) - kf_lgamma(a) - kf_lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp() / a / f
}

/// Regularized incomplete beta function (HTSlib's `kf_betai`).
pub fn kf_betai(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        kf_betai_aux(a, b, x)
    } else {
        1.0 - kf_betai_aux(b, a, 1.0 - x)
    }
}

/// log(n choose k) with the C library's `lgamma`.
fn lbinom(n: i32, k: i32) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    lgamma(f64::from(n.wrapping_add(1)))
        - lgamma(f64::from(k.wrapping_add(1)))
        - lgamma(f64::from(n.wrapping_sub(k).wrapping_add(1)))
}

fn hypergeo(n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
    (lbinom(n1_, n11) + lbinom(n.wrapping_sub(n1_), n_1.wrapping_sub(n11)) - lbinom(n, n_1)).exp()
}

#[derive(Default)]
struct HgAcc {
    n11: i32,
    n1_: i32,
    n_1: i32,
    n: i32,
    p: f64,
}

impl HgAcc {
    /// HTSlib's `hypergeo_acc`: the probability of a table, incrementally when only `n11`
    /// changes by one.
    fn next(&mut self, n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
        if n1_ != 0 || n_1 != 0 || n != 0 {
            self.n11 = n11;
            self.n1_ = n1_;
            self.n_1 = n_1;
            self.n = n;
        } else {
            // The counts wrap as C's `int` does on x86, for tables made of wrapped I16 values.
            let n22_less_n11 = self.n.wrapping_sub(self.n1_).wrapping_sub(self.n_1);
            if n11 % 11 != 0 && n11.wrapping_add(n22_less_n11) != 0 {
                if n11 == self.n11.wrapping_add(1) {
                    self.p *= f64::from(self.n1_.wrapping_sub(self.n11)) / f64::from(n11)
                        * f64::from(self.n_1.wrapping_sub(self.n11))
                        / f64::from(n11.wrapping_add(n22_less_n11));
                    self.n11 = n11;
                    return self.p;
                }
                if n11 == self.n11.wrapping_sub(1) {
                    self.p *= f64::from(self.n11) / f64::from(self.n1_.wrapping_sub(n11))
                        * f64::from(self.n11.wrapping_add(n22_less_n11))
                        / f64::from(self.n_1.wrapping_sub(n11));
                    self.n11 = n11;
                    return self.p;
                }
            }
            self.n11 = n11;
        }
        self.p = hypergeo(self.n11, self.n1_, self.n_1, self.n);
        self.p
    }
}

/// Fisher's exact test (HTSlib's `kt_fisher_exact`); returns the two-tailed P-value.
///
/// The sums wrap as C's `int` does on x86: `call` passes I16 counts of 2^31 or more as
/// `INT_MIN` (see `ccall::c_int`), and bcftools then works on the wrapped sums.
pub fn kt_fisher_exact_two(n11: i32, n12: i32, n21: i32, n22: i32) -> f64 {
    let n1_ = n11.wrapping_add(n12);
    let n_1 = n11.wrapping_add(n21);
    let n = n1_.wrapping_add(n21).wrapping_add(n22);
    let max = n_1.min(n1_);
    let min = n1_.wrapping_add(n_1).wrapping_sub(n).max(0);
    if min == max {
        return 1.0;
    }
    let mut aux = HgAcc::default();
    let q = aux.next(n11, n1_, n_1, n);
    if q == 0.0 {
        return 0.0;
    }
    let mut p = aux.next(min, 0, 0, 0);
    let mut left = 0.0;
    let mut i = min.wrapping_add(1);
    while p < 0.99999999 * q && i <= max {
        left += p;
        p = aux.next(i, 0, 0, 0);
        i = i.wrapping_add(1);
    }
    if p < 1.00000001 * q {
        left += p;
    }
    p = aux.next(max, 0, 0, 0);
    let mut right = 0.0;
    let mut j = max.wrapping_sub(1);
    while p < 0.99999999 * q && j >= 0 {
        right += p;
        p = aux.next(j, 0, 0, 0);
        j -= 1;
    }
    if p < 1.00000001 * q {
        right += p;
    }
    let two = left + right;
    if two > 1.0 { 1.0 } else { two }
}

/// Brent's method for a one-dimensional minimum (bcftools' `kmin_brent`); returns the
/// minimum point.
pub fn kmin_brent<F: FnMut(f64) -> f64>(mut func: F, mut a: f64, mut b: f64, tol: f64) -> f64 {
    const GOLD1: f64 = 1.6180339887;
    const GOLD2: f64 = 0.3819660113;
    const TINY: f64 = 1e-20;
    const MAX_ITER: usize = 100;

    let mut fa = func(a);
    let mut fb = func(b);
    if fb > fa {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = b + GOLD1 * (b - a);
    let mut fc = func(c);
    let mut u;
    let mut fu;
    while fb > fc {
        let bound = b + 100.0 * (c - b);
        let r = (b - a) * (fb - fc);
        let q = (b - c) * (fb - fa);
        let tmp = if (q - r).abs() < TINY {
            if q > r { TINY } else { 0.0 - TINY }
        } else {
            q - r
        };
        u = b - ((b - c) * q - (b - a) * r) / (2.0 * tmp);
        if (b > u && u > c) || (b < u && u < c) {
            fu = func(u);
            if fu < fc {
                a = b;
                b = u;
                fa = fb;
                fb = fu;
                break;
            } else if fu > fb {
                c = u;
                break;
            }
            u = c + GOLD1 * (c - b);
            fu = func(u);
        } else if (c > u && u > bound) || (c < u && u < bound) {
            fu = func(u);
            if fu < fc {
                b = c;
                c = u;
                u = c + GOLD1 * (c - b);
                fb = fc;
                fc = fu;
                fu = func(u);
            } else {
                a = b;
                b = c;
                c = u;
                fa = fb;
                fb = fc;
                break;
            }
        } else if (u > bound && bound > c) || (u < bound && bound < c) {
            u = bound;
            fu = func(u);
        } else {
            u = c + GOLD1 * (c - b);
            fu = func(u);
        }
        a = b;
        b = c;
        c = u;
        fa = fb;
        fb = fc;
        fc = fu;
    }
    let _ = fa;
    if a > c {
        std::mem::swap(&mut a, &mut c);
    }

    let mut e: f64 = 0.0;
    let mut d: f64 = 0.0;
    let mut w = b;
    let mut v = b;
    let mut fv = fb;
    let mut fw = fb;
    for _ in 0..MAX_ITER {
        let mid = 0.5 * (a + c);
        let tol1 = tol * b.abs() + TINY;
        let tol2 = 2.0 * tol1;
        if (b - mid).abs() <= (tol2 - 0.5 * (c - a)) {
            return b;
        }
        if e.abs() > tol1 {
            let r = (b - w) * (fb - fv);
            let mut q = (b - v) * (fb - fw);
            let mut p = (b - v) * q - (b - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = 0.0 - p;
            } else {
                q = 0.0 - q;
            }
            let eold = e;
            e = d;
            if p.abs() >= (0.5 * q * eold).abs() || p <= q * (a - b) || p >= q * (c - b) {
                e = if b >= mid { a - b } else { c - b };
                d = GOLD2 * e;
            } else {
                d = p / q;
                u = b + d;
                if u - a < tol2 || c - u < tol2 {
                    d = if mid > b { tol1 } else { 0.0 - tol1 };
                }
            }
        } else {
            e = if b >= mid { a - b } else { c - b };
            d = GOLD2 * e;
        }
        u = if d.abs() >= tol1 {
            b + d
        } else {
            b + if d > 0.0 { tol1 } else { -tol1 }
        };
        fu = func(u);
        if fu <= fb {
            if u >= b {
                a = b;
            } else {
                c = b;
            }
            v = w;
            w = b;
            b = u;
            fv = fw;
            fw = fb;
            fb = fu;
        } else {
            if u < b {
                a = u;
            } else {
                c = u;
            }
            if fu <= fw || w == b {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if fu <= fv || v == b || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fisher_matches_known_values() {
        assert_eq!(kt_fisher_exact_two(0, 0, 0, 0), 1.0);
        let p = kt_fisher_exact_two(3, 1, 1, 3);
        assert!((p - 0.48571428571428565).abs() < 1e-12, "{p}");
    }

    #[test]
    fn fisher_sums_wrap_like_c() {
        let m = i32::MIN;
        assert_eq!(kt_fisher_exact_two(m, m, m, m), 1.0);
        assert_eq!(kt_fisher_exact_two(m, 0, 0, 5), 0.0);
        assert_eq!(kt_fisher_exact_two(5, 6, m, m), 0.0);
        // The largest count is INT_MIN, so the right tail starts from INT_MIN - 1.
        assert_eq!(kt_fisher_exact_two(m, 0, 0, 0), 1.0);
    }

    #[test]
    fn gamma_functions_are_sane() {
        assert!((kf_lgamma(5.0) - 24f64.ln()).abs() < 1e-10);
        assert!((kf_gammaq(0.5, 0.0) - 1.0).abs() < 1e-12);
        assert!((kf_betai(1.0, 1.0, 0.25) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn brent_finds_a_parabola_minimum() {
        let x = kmin_brent(|x| (x - 0.3) * (x - 0.3), 0.1, 0.9, 1e-5);
        assert!((x - 0.3).abs() < 1e-4, "{x}");
    }
}
