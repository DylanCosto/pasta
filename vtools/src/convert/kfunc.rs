//! The two-sided binomial test of `%PBINOM`, from HTSlib's `kfunc.c` (`kf_lgamma`,
//! `kf_betai`) and bcftools' `calc_binom_two_sided`.

/// Relative tolerance of the continued fraction (HTSlib's `KF_GAMMA_EPS`).
const GAMMA_EPS: f64 = 1e-14;
/// Floor that keeps the continued fraction from dividing by zero (HTSlib's `KF_TINY`).
const TINY: f64 = 1e-290;

/// The log gamma function (HTSlib's `kf_lgamma`, Lanczos approximation).
fn lgamma(z: f64) -> f64 {
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

/// HTSlib's `kf_betai_aux`: the continued fraction of the incomplete beta function.
fn betai_aux(a: f64, b: f64, x: f64) -> f64 {
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
        if d < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < GAMMA_EPS {
            break;
        }
    }
    (lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp() / a / f
}

/// The regularized incomplete beta function (HTSlib's `kf_betai`).
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        betai_aux(a, b, x)
    } else {
        1.0 - betai_aux(b, a, 1.0 - x)
    }
}

/// The two-sided binomial test of `na` against `nb` successes (bcftools'
/// `calc_binom_two_sided`): -1 when both are zero.
pub(super) fn binom_two_sided(na: i32, nb: i32, aprob: f64) -> f64 {
    if na == 0 && nb == 0 {
        return -1.0;
    }
    if na == nb {
        return 1.0;
    }
    let prob = if na > nb {
        2.0 * betai(f64::from(na), f64::from(nb) + 1.0, aprob)
    } else {
        2.0 * betai(f64::from(nb), f64::from(na) + 1.0, aprob)
    };
    // Machine precision can push it over one (kf_betai(1, 0, 0.5), for example).
    if prob > 1.0 { 1.0 } else { prob }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binomial_test_is_symmetric_and_bounded() {
        assert_eq!(binom_two_sided(0, 0, 0.5), -1.0);
        assert_eq!(binom_two_sided(3, 3, 0.5), 1.0);
        let p = binom_two_sided(10, 0, 0.5);
        assert!((p - 2.0 * 0.5f64.powi(10)).abs() < 1e-9, "p = {p}");
        assert_eq!(binom_two_sided(2, 9, 0.5), binom_two_sided(9, 2, 0.5));
    }
}
