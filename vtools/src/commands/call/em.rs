//! Single-locus EM estimates of the consensus caller (bcftools' `em.c`, `bcf_em1` as
//! `bcftools call -c` runs it: one group, frequency and genotype-frequency estimates and the
//! HWE test).

use super::kfunc::{kf_gammaq, kmin_brent};

const ITER_MAX: usize = 50;
const ITER_TRY: usize = 10;
const EPS: f64 = 1e-5;

/// A rough estimate of the alternate allele frequency (`est_freq`), or -1 without data.
fn est_freq(pdg: &[f64]) -> f64 {
    let mut gcnt = [0i32; 3];
    for p in pdg.as_chunks::<3>().0 {
        if p[0] != 1. || p[1] != 1. || p[2] != 1. {
            let mut which = if p[0] > p[1] { 0 } else { 1 };
            which = if p[which] > p[2] { which } else { 2 };
            gcnt[which] += 1;
        }
    }
    let tmp1 = gcnt[0] + gcnt[1] + gcnt[2];
    if tmp1 == 0 {
        -1.0
    } else {
        (0.5 * f64::from(gcnt[1]) + f64::from(gcnt[2])) / f64::from(tmp1)
    }
}

/// The negative log-likelihood of frequency `f` (`prob1` in `em.c`).
fn prob1(f: f64, pdg: &[f64]) -> f64 {
    if !(0.0..=1.0).contains(&f) {
        return 1e300;
    }
    let f3 = [(1. - f) * (1. - f), 2. * f * (1. - f), f * f];
    let mut p = 1.0;
    let mut l = 0.0;
    for d in pdg.as_chunks::<3>().0 {
        p *= d[0] * f3[0] + d[1] * f3[1] + d[2] * f3[2];
        if p < 1e-200 {
            l -= p.ln();
            p = 1.0;
        }
    }
    l - p.ln()
}

/// One EM iteration of the frequency estimate (`freq_iter`); returns the change.
fn freq_iter(f: &mut f64, pdg: &[f64]) -> f64 {
    let f0 = *f;
    let f3 = [(1. - f0) * (1. - f0), 2. * f0 * (1. - f0), f0 * f0];
    let mut acc = 0.0;
    for d in pdg.as_chunks::<3>().0 {
        acc += (d[1] * f3[1] + 2. * d[2] * f3[2]) / (d[0] * f3[0] + d[1] * f3[1] + d[2] * f3[2]);
    }
    let n = pdg.len() / 3;
    acc /= (n * 2) as f64;
    let err = (acc - *f).abs();
    *f = acc;
    err
}

/// EM, switching to Brent's method when EM converges slowly (`freqml`).
fn freqml(f0: f64, pdg: &[f64]) -> f64 {
    let mut f = f0;
    let mut i = 0;
    while i < ITER_TRY {
        if freq_iter(&mut f, pdg) < EPS {
            break;
        }
        i += 1;
    }
    if i == ITER_TRY {
        let a = if f0 == f { 0.5 * f0 } else { f0 };
        f = kmin_brent(|x| prob1(x, pdg), a, f, EPS);
    }
    f
}

/// One EM iteration of the genotype-frequency estimate (`g3_iter`).
fn g3_iter(g: &mut [f64; 3], pdg: &[f64]) -> f64 {
    let n = (pdg.len() / 3) as f64;
    let mut gg = [0.0; 3];
    for d in pdg.as_chunks::<3>().0 {
        let tmp = [d[0] * g[0], d[1] * g[1], d[2] * g[2]];
        let sum = (tmp[0] + tmp[1] + tmp[2]) * n;
        gg[0] += tmp[0] / sum;
        gg[1] += tmp[1] / sum;
        gg[2] += tmp[2] / sum;
    }
    let d0 = (gg[0] - g[0]).abs();
    let d1 = (gg[1] - g[1]).abs();
    let d2 = (gg[2] - g[2]).abs();
    let mut err = if d0 > d1 { d0 } else { d1 };
    err = if err > d2 { err } else { d2 };
    *g = gg;
    err
}

/// `bcf_em1` with flag 0x1ff and no second group: `x[0]` is the reference allele
/// frequency, `x[1..4]` the genotype frequencies, `x[4]` the HWE P-value; -1 when not
/// computed.
pub fn em1(pdg: &[f64], n_allele: usize) -> [f64; 10] {
    let mut x = [-1.0; 10];
    if n_allele < 2 {
        return x;
    }
    x[0] = est_freq(pdg);
    if x[0] < 0.0 {
        return x;
    }
    x[0] = freqml(x[0], pdg);
    let f3 = [
        (1. - x[0]) * (1. - x[0]),
        2. * x[0] * (1. - x[0]),
        x[0] * x[0],
    ];
    let mut g = f3;
    for _ in 0..ITER_MAX {
        if g3_iter(&mut g, pdg) < EPS {
            break;
        }
    }
    let mut r = 1.0;
    for p in pdg.as_chunks::<3>().0 {
        r *= (p[0] * g[0] + p[1] * g[1] + p[2] * g[2])
            / (p[0] * f3[0] + p[1] * f3[1] + p[2] * f3[2]);
    }
    x[1] = g[0];
    x[2] = g[1];
    x[3] = g[2];
    x[4] = kf_gammaq(0.5, r.ln());
    x
}
