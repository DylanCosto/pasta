//! Statistical tests from HTSlib's `kfunc.c` and bcftools' `calc_binom_two_sided`.

unsafe extern "C" {
    /// The C library's `lgamma`, which HTSlib's Fisher test uses; libm is linked by std.
    safe fn lgamma(x: f64) -> f64;
}

const KF_GAMMA_EPS: f64 = 1e-14;
const KF_TINY: f64 = 1e-290;

/// HTSlib's `kf_lgamma`: the log gamma function by the Lanczos approximation.
#[expect(
    clippy::excessive_precision,
    reason = "HTSlib's constants, kept digit for digit"
)]
fn kf_lgamma(z: f64) -> f64 {
    let mut x = 0.0;
    x += 0.165_947_018_740_846_2e-06 / (z + 7.0);
    x += 0.993_493_711_393_074_8e-05 / (z + 6.0);
    x -= 0.138_571_033_129_652_6 / (z + 5.0);
    x += 12.507_343_240_090_56 / (z + 4.0);
    x -= 176.615_029_149_838_6 / (z + 3.0);
    x += 771.323_428_775_767_4 / (z + 2.0);
    x -= 1_259.139_216_722_289 / (z + 1.0);
    x += 676.520_368_121_883_5 / z;
    x += 0.999_999_999_999_518_3;
    x.ln() - 5.581_061_466_795_327_77 - z + (z - 0.5) * (z + 6.5).ln()
}

/// The continued fraction of the regularized incomplete beta function (`kf_betai_aux`).
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
    // Modified Lentz's algorithm for computing the continued fraction.
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

/// The regularized incomplete beta function (HTSlib's `kf_betai`).
fn kf_betai(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        kf_betai_aux(a, b, x)
    } else {
        1.0 - kf_betai_aux(b, a, 1.0 - x)
    }
}

/// The two-sided binomial test of `na` and `nb` observations with probability `aprob`, or -1
/// when both are zero (bcftools' `calc_binom_two_sided`).
pub(crate) fn binom_two_sided(na: i32, nb: i32, aprob: f64) -> f64 {
    if na == 0 && nb == 0 {
        return -1.0;
    }
    if na == nb {
        return 1.0;
    }
    let prob = if na > nb {
        2.0 * kf_betai(f64::from(na), f64::from(nb) + 1.0, aprob)
    } else {
        2.0 * kf_betai(f64::from(nb), f64::from(na) + 1.0, aprob)
    };
    // Not `min`: a NaN (from the missing-value sentinel of an INFO vector) stays NaN.
    if prob > 1.0 { 1.0 } else { prob }
}

/// `log(n choose k)`.
fn lbinom(n: i32, k: i32) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    lgamma(f64::from(n) + 1.0) - lgamma(f64::from(k) + 1.0) - lgamma(f64::from(n - k) + 1.0)
}

/// The hypergeometric probability of a 2x2 table.
fn hypergeo(n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
    (lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1)).exp()
}

/// The state of [`hypergeo_acc`].
#[derive(Default)]
struct HgAcc {
    n11: i32,
    n1_: i32,
    n_1: i32,
    n: i32,
    p: f64,
}

/// The incremental hypergeometric probability (HTSlib's `hypergeo_acc`): with the margins 0,
/// only `n11` changed. Inlined into the tail loops of [`fisher_exact_two`] so that the
/// probability stays in a register (the loops run hundreds of times per test on deep data).
#[inline(always)]
fn hypergeo_acc(n11: i32, n1_: i32, n_1: i32, n: i32, aux: &mut HgAcc) -> f64 {
    if n1_ != 0 || n_1 != 0 || n != 0 {
        aux.n11 = n11;
        aux.n1_ = n1_;
        aux.n_1 = n_1;
        aux.n = n;
    } else {
        if n11 % 11 != 0 && n11 + aux.n - aux.n1_ - aux.n_1 != 0 {
            if n11 == aux.n11 + 1 {
                aux.p *= f64::from(aux.n1_ - aux.n11) / f64::from(n11)
                    * f64::from(aux.n_1 - aux.n11)
                    / f64::from(n11 + aux.n - aux.n1_ - aux.n_1);
                aux.n11 = n11;
                return aux.p;
            }
            if n11 == aux.n11 - 1 {
                aux.p *= f64::from(aux.n11) / f64::from(aux.n1_ - n11)
                    * f64::from(aux.n11 + aux.n - aux.n1_ - aux.n_1)
                    / f64::from(aux.n_1 - n11);
                aux.n11 = n11;
                return aux.p;
            }
        }
        aux.n11 = n11;
    }
    aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
    aux.p
}

/// Fisher's exact test of a 2x2 table, returning the two-tailed p-value (HTSlib's
/// `kt_fisher_exact`).
pub(super) fn fisher_exact_two(n11: i32, n12: i32, n21: i32, n22: i32) -> f64 {
    // C's int arithmetic: the missing-value sentinel of an INFO vector wraps around.
    let n1_ = n11.wrapping_add(n12);
    let n_1 = n11.wrapping_add(n21);
    let n = n1_.wrapping_add(n21).wrapping_add(n22);
    let max = n_1.min(n1_);
    let min = n1_.wrapping_add(n_1).wrapping_sub(n).max(0);
    if min == max {
        return 1.0;
    }
    let mut aux = HgAcc::default();
    let q = hypergeo_acc(n11, n1_, n_1, n, &mut aux);
    if q == 0.0 {
        return 0.0;
    }
    // left tail
    let mut p = hypergeo_acc(min, 0, 0, 0, &mut aux);
    let mut left = 0.0;
    let mut i = min + 1;
    while p < 0.999_999_99 * q && i <= max {
        left += p;
        p = hypergeo_acc(i, 0, 0, 0, &mut aux);
        i += 1;
    }
    if p < 1.000_000_01 * q {
        left += p;
    }
    // right tail
    p = hypergeo_acc(max, 0, 0, 0, &mut aux);
    let mut right = 0.0;
    let mut j = max - 1;
    while p < 0.999_999_99 * q && j >= 0 {
        right += p;
        p = hypergeo_acc(j, 0, 0, 0, &mut aux);
        j -= 1;
    }
    if p < 1.000_000_01 * q {
        right += p;
    }
    // Not `min`: a NaN stays NaN, as in HTSlib.
    let two = left + right;
    if two > 1.0 { 1.0 } else { two }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fisher_matches_known_values() {
        // R: fisher.test(matrix(c(3,1,1,3),2))$p.value = 0.4857143
        let p = fisher_exact_two(3, 1, 1, 3);
        assert!((p - 0.485_714_3).abs() < 1e-6, "{p}");
        assert_eq!(fisher_exact_two(0, 0, 0, 0), 1.0);
    }

    #[test]
    fn binomial_two_sided() {
        assert_eq!(binom_two_sided(0, 0, 0.5), -1.0);
        assert_eq!(binom_two_sided(5, 5, 0.5), 1.0);
        // R: binom.test(10, 12)$p.value = 0.03857422
        let p = binom_two_sided(10, 2, 0.5);
        assert!((p - 0.038_574_22).abs() < 1e-6, "{p}");
        // The missing-value sentinel of an INFO vector (`binom(INFO/TAG)` of `.,6`) gives NaN,
        // which compares false, as in bcftools.
        assert!(binom_two_sided(i32::MIN, 6, 0.5).is_nan());
    }
}
