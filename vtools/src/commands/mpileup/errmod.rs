//! HTSlib's error model for genotype likelihoods (`errmod.c`), the `drand48` generator it
//! shuffles deep pileups with, and the `kfunc.c` statistics bcftools mpileup uses.
//!
//! The floating-point operations follow the C code, in single or double precision as it uses
//! them, and call the same libm functions, so the likelihoods match bcftools bit for bit.

use std::sync::LazyLock;

/// bcftools' `CALL_DEFTHETA`: errmod is set up with a dependency of `1 - theta`.
const THETA: f64 = 0.83;
/// The `eta` HTSlib's `errmod_init` uses.
const ETA: f64 = 0.03;

unsafe extern "C" {
    /// libm's reentrant `lgamma`.
    fn lgamma_r(x: f64, sign: *mut libc::c_int) -> f64;
}

/// libm's `lgamma`, as the C code calls it (through the reentrant variant, which gives the
/// same value without the global sign).
fn lgamma(x: f64) -> f64 {
    let mut sign: libc::c_int = 0;
    // SAFETY: `lgamma_r` only reads `x` and writes the sign through the valid pointer.
    unsafe { lgamma_r(x, &raw mut sign) }
}

/// `errmod_t`: the tables of `cal_coef`.
pub struct ErrorModel {
    fk: Vec<f64>,
    beta: Vec<f64>,
    lhet: Vec<f64>,
}

/// The one error model bcftools mpileup uses (`errmod_init(1 - 0.83)`), built once.
pub static ERROR_MODEL: LazyLock<ErrorModel> = LazyLock::new(|| ErrorModel::new(1.0 - THETA));

impl ErrorModel {
    /// `errmod_init` and `cal_coef`.
    fn new(depcorr: f64) -> Self {
        let mut fk = vec![0.0; 256];
        fk[0] = 1.0;
        for (n, v) in fk.iter_mut().enumerate().skip(1) {
            *v = (1.0 - depcorr).powf(n as f64) * (1.0 - ETA) + ETA;
        }

        // logbinomial_table(256)
        let mut lc = vec![0.0f64; 256 * 256];
        for n in 1..256usize {
            let lfn = lgamma(n as f64 + 1.0);
            for k in 1..=n {
                lc[n << 8 | k] = lfn - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0);
            }
        }

        let mut beta = vec![0.0f64; 256 * 256 * 64];
        for q in 1..64usize {
            let e = 10f64.powf(-(q as f64) / 10.0);
            let le = e.ln();
            let le1 = (1.0 - e).ln();
            for n in 1..=255usize {
                // Laid out by depth first (HTSlib uses q<<16|n<<8|k), so that the bases of one
                // pileup, which share their depth, read nearby entries.
                let base = n << 14 | q << 8;
                let mut sum1 = lc[n << 8 | n] + n as f64 * le;
                beta[base + n] = f64::INFINITY;
                for k in (0..n).rev() {
                    let sum = sum1
                        + (lc[n << 8 | k] + k as f64 * le + (n - k) as f64 * le1 - sum1)
                            .exp()
                            .ln_1p();
                    beta[base + k] = -10.0 / std::f64::consts::LN_10 * (sum1 - sum);
                    sum1 = sum;
                }
            }
        }

        let mut lhet = vec![0.0f64; 256 * 256];
        for n in 0..256usize {
            for k in 0..256usize {
                lhet[n << 8 | k] = lc[n << 8 | k] - std::f64::consts::LN_2 * n as f64;
            }
        }

        Self { fk, beta, lhet }
    }

    /// `errmod_cal`: the phred-scaled likelihood `q[i*m+j]` of each genotype from the bases of
    /// one sample (6-bit quality, strand, 4-bit base each). More than 255 bases are shuffled
    /// with `rng` and the first 255 kept, as HTSlib does.
    pub fn cal(&self, bases: &mut [u16], m: usize, q: &mut [f32], rng: &mut Rng48) {
        q[..m * m].fill(0.0);
        let mut n = bases.len();
        if n == 0 {
            return;
        }
        if n > 255 {
            rng.shuffle(bases);
            n = 255;
        }
        let bases = &mut bases[..n];
        bases.sort_unstable();

        let mut fsum = [0f64; 16];
        let mut bsum = [0f64; 16];
        let mut c = [0u32; 16];
        let mut w = [0usize; 32];
        for &b in bases.iter().rev() {
            let qual = ((b >> 5) as usize).clamp(4, 63);
            let basestrand = (b & 0x1f) as usize;
            let base = (b & 0xf) as usize;
            let fk = self.fk[w[basestrand]];
            fsum[base] += fk;
            bsum[base] += fk * self.beta[n << 14 | qual << 8 | c[base] as usize];
            c[base] += 1;
            w[basestrand] += 1;
        }

        for j in 0..m {
            // homozygous
            let mut tmp1: f32 = 0.0;
            let mut tmp2: u32 = 0;
            for k in 0..m {
                if k == j {
                    continue;
                }
                tmp1 = (f64::from(tmp1) + bsum[k]) as f32;
                tmp2 += c[k];
            }
            if tmp2 != 0 {
                q[j * m + j] = tmp1;
            }
            // heterozygous
            for k in j + 1..m {
                let cjk = (c[j] + c[k]) as usize;
                let mut tmp1: f32 = 0.0;
                let mut tmp2: u32 = 0;
                for i in 0..m {
                    if i == j || i == k {
                        continue;
                    }
                    tmp1 = (f64::from(tmp1) + bsum[i]) as f32;
                    tmp2 += c[i];
                }
                let lhet = -4.343 * self.lhet[cjk << 8 | c[k] as usize];
                let v = if tmp2 != 0 {
                    (lhet + f64::from(tmp1)) as f32
                } else {
                    lhet as f32
                };
                q[j * m + k] = v;
                q[k * m + j] = v;
            }
            for k in 0..m {
                if q[j * m + k] < 0.0 {
                    q[j * m + k] = 0.0;
                }
            }
        }
    }
}

/// glibc's `drand48` (a 48-bit linear congruential generator), which HTSlib's `ks_shuffle`
/// draws from; bcftools mpileup seeds it with `srand48(--seed)`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Rng48 {
    x: u64,
}

const LCG_A: u64 = 0x5_DEEC_E66D;
const LCG_C: u64 = 0xB;
const MASK48: u64 = (1 << 48) - 1;

impl Rng48 {
    /// `srand48(seed)`.
    pub fn new(seed: i64) -> Self {
        Self {
            x: ((seed as u64 & 0xffff_ffff) << 16 | 0x330e) & MASK48,
        }
    }

    /// `drand48()`.
    pub fn next_f64(&mut self) -> f64 {
        self.x = (LCG_A.wrapping_mul(self.x).wrapping_add(LCG_C)) & MASK48;
        self.x as f64 / (1u64 << 48) as f64
    }

    /// The state after `n` more draws, computed in O(log n).
    pub fn advanced(self, mut n: u64) -> Self {
        // Compose x -> a*x + c with itself n times.
        let (mut acc_a, mut acc_c) = (1u64, 0u64);
        let (mut a, mut c) = (LCG_A, LCG_C);
        while n > 0 {
            if n & 1 != 0 {
                acc_a = acc_a.wrapping_mul(a) & MASK48;
                acc_c = (acc_c.wrapping_mul(a).wrapping_add(c)) & MASK48;
            }
            c = (c.wrapping_mul(a).wrapping_add(c)) & MASK48;
            a = a.wrapping_mul(a) & MASK48;
            n >>= 1;
        }
        Self {
            x: (acc_a.wrapping_mul(self.x).wrapping_add(acc_c)) & MASK48,
        }
    }

    /// `ks_shuffle`: a Fisher-Yates shuffle drawing `len - 1` numbers.
    fn shuffle(&mut self, a: &mut [u16]) {
        let mut i = a.len();
        while i > 1 {
            let j = (self.next_f64() * i as f64) as usize;
            a.swap(j, i - 1);
            i -= 1;
        }
    }
}

/// `kf_erfc`: the complementary error function (AS66).
pub fn kf_erfc(x: f64) -> f64 {
    const P0: f64 = 220.2068679123761;
    const P1: f64 = 221.2135961699311;
    const P2: f64 = 112.0792914978709;
    const P3: f64 = 33.912866078383;
    const P4: f64 = 6.37396220353165;
    const P5: f64 = 0.7003830644436881;
    const P6: f64 = 0.03526249659989109;
    const Q0: f64 = 440.4137358247522;
    const Q1: f64 = 793.8265125199484;
    const Q2: f64 = 637.3336333788311;
    const Q3: f64 = 296.5642487796737;
    const Q4: f64 = 86.78073220294608;
    const Q5: f64 = 16.06417757920695;
    const Q6: f64 = 1.755667163182642;
    const Q7: f64 = 0.08838834764831844;
    let z = x.abs() * std::f64::consts::SQRT_2;
    if z > 37.0 {
        return if x > 0.0 { 0.0 } else { 2.0 };
    }
    let expntl = (z * z * -0.5).exp();
    let p = if z < 10.0 / std::f64::consts::SQRT_2 {
        expntl * ((((((P6 * z + P5) * z + P4) * z + P3) * z + P2) * z + P1) * z + P0)
            / (((((((Q7 * z + Q6) * z + Q5) * z + Q4) * z + Q3) * z + Q2) * z + Q1) * z + Q0)
    } else {
        expntl / 2.506628274631001 / (z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 0.65)))))
    };
    if x > 0.0 { 2.0 * p } else { 2.0 * (1.0 - p) }
}

/// `lbinom`: log of the binomial coefficient.
fn lbinom(n: i32, k: i32) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    lgamma(f64::from(n) + 1.0) - lgamma(f64::from(k) + 1.0) - lgamma(f64::from(n - k) + 1.0)
}

fn hypergeo(n11: i32, n1_: i32, n_1: i32, n: i32) -> f64 {
    (lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1)).exp()
}

/// `hgacc_t`
#[derive(Default)]
struct HgAcc {
    n11: i32,
    n1_: i32,
    n_1: i32,
    n: i32,
    p: f64,
}

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

/// `kt_fisher_exact`: returns the probability of the table and the two-tailed p-value.
pub fn kt_fisher_exact(n11: i32, n12: i32, n21: i32, n22: i32) -> (f64, f64) {
    let n1_ = n11 + n12;
    let n_1 = n11 + n21;
    let n = n11 + n12 + n21 + n22;
    let max = n_1.min(n1_);
    let min = (n1_ + n_1 - n).max(0);
    if min == max {
        return (1.0, 1.0);
    }
    let mut aux = HgAcc::default();
    let q = hypergeo_acc(n11, n1_, n_1, n, &mut aux);
    if q == 0.0 {
        return (0.0, 0.0);
    }
    let mut p = hypergeo_acc(min, 0, 0, 0, &mut aux);
    let mut left = 0.0;
    let mut i = min + 1;
    while p < 0.99999999 * q && i <= max {
        left += p;
        p = hypergeo_acc(i, 0, 0, 0, &mut aux);
        i += 1;
    }
    if p < 1.00000001 * q {
        left += p;
    }
    p = hypergeo_acc(max, 0, 0, 0, &mut aux);
    let mut right = 0.0;
    let mut j = max - 1;
    while p < 0.99999999 * q && j >= 0 {
        right += p;
        p = hypergeo_acc(j, 0, 0, 0, &mut aux);
        j -= 1;
    }
    if p < 1.00000001 * q {
        right += p;
    }
    let two = (left + right).min(1.0);
    (q, two)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drand48_matches_glibc() {
        // glibc: srand48(0) then drand48() three times; srand48(5) then drand48().
        let mut rng = Rng48::new(0);
        let got: Vec<f64> = (0..3).map(|_| rng.next_f64()).collect();
        assert_eq!(
            got,
            [
                0.170_828_036_106_289_72,
                0.749_901_980_484_963_8,
                0.096_371_655_623_567_42
            ]
        );
        assert_eq!(Rng48::new(5).next_f64(), 0.524_839_579_434_232_2);
    }

    #[test]
    fn advance_matches_stepping() {
        let mut a = Rng48::new(42);
        let b = a.advanced(1000);
        for _ in 0..1000 {
            a.next_f64();
        }
        assert_eq!(a, b);
    }

    #[test]
    fn fisher_of_a_balanced_table() {
        let (q, two) = kt_fisher_exact(10, 10, 10, 10);
        assert!(q > 0.0 && q < 1.0, "{q}");
        assert!((two - 1.0).abs() < 1e-9, "{two}");
    }
}
