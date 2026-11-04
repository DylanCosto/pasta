//! The allele-frequency model of the consensus caller (bcftools' `prob1.c`), for one group
//! of samples (the `-1` two-group contrast was never exposed by `bcftools call`).

use super::f80::F80;

/// Likelihoods below this are squashed to 0 (`TINY` in `prob1.c`).
const TINY: f64 = 1e-20;
/// The default indel prior scale (`MC_DEF_INDEL`).
const MC_DEF_INDEL: f64 = 0.15;

/// The type of a site, which picks the prior of [`P1Aux::cal`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum VariantKind {
    /// A SNP (`bcf_is_snp`).
    Snp,
    /// Anything else.
    Indel,
}

/// What [`P1Aux::call_gt`] calls: the most likely genotype, or the reference genotype.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SiteCall {
    /// A variant site: the most likely genotype.
    Variant,
    /// A reference site: the reference genotype, with its quality.
    Reference,
}

/// bcftools' `bcf_p1aux_t`.
#[derive(Clone, Debug)]
pub struct P1Aux {
    /// Number of samples.
    n: usize,
    /// Number of chromosomes.
    m: usize,
    kind: VariantKind,
    /// Per-sample ploidy, or `None` when all samples are diploid.
    ploidy: Option<Vec<u8>>,
    q2p: [f64; 256],
    /// P(D|g) per sample, in reverse genotype order (alt/alt first).
    pub pdg: Vec<f64>,
    phi: Vec<f64>,
    phi_indel: Vec<f64>,
    z: Vec<f64>,
    zswap: Vec<f64>,
    afs1: Vec<f64>,
    /// Scratch space for the terms of `mc_cal_afs`'s sums.
    scratch: Vec<f64>,
    t: f64,
}

/// bcftools' `bcf_p1rst_t`, the fields `bcftools call -c` uses.
#[derive(Clone, Copy, Debug, Default)]
pub struct P1Result {
    pub rank0: i32,
    pub ac: i32,
    pub f_exp: f64,
    pub p_ref_folded: f64,
    pub p_ref: f64,
    pub p_var_folded: f64,
    pub p_var: f64,
}

impl P1Aux {
    /// bcftools' `bcf_p1_init` with the default prior (`MC_PTYPE_FULL`, theta 1e-3).
    pub fn new(n_smpl: usize, ploidy: Option<&[u8]>) -> Self {
        let mut m = 2 * n_smpl;
        let mut pl = None;
        if let Some(p) = ploidy {
            m = p.iter().map(|&x| usize::from(x)).sum();
            if m != 2 * n_smpl {
                pl = Some(p.to_vec());
            }
        }
        let mut q2p = [0.0; 256];
        for (i, q) in q2p.iter_mut().enumerate() {
            *q = 10f64.powf(-(i as f64) / 10.0);
        }
        let mut aux = Self {
            n: n_smpl,
            m,
            kind: VariantKind::Snp,
            ploidy: pl,
            q2p,
            pdg: vec![0.0; 3 * n_smpl],
            phi: vec![0.0; m + 1],
            phi_indel: vec![0.0; m + 1],
            z: vec![0.0; m + 1],
            zswap: vec![0.0; m + 1],
            afs1: vec![0.0; m + 1],
            scratch: Vec::with_capacity(m + 1),
            t: 0.0,
        };
        aux.init_prior(1e-3);
        aux
    }

    /// `init_prior` with `MC_PTYPE_FULL`, then `bcf_p1_indel_prior`.
    fn init_prior(&mut self, theta: f64) {
        let m = self.m;
        let mut sum = 0.0;
        for i in 0..m {
            self.phi[i] = theta / (m - i) as f64;
            sum += self.phi[i];
        }
        self.phi[m] = 1.0 - sum;
        for i in 0..m {
            self.phi_indel[i] = self.phi[i] * MC_DEF_INDEL;
        }
        self.phi_indel[m] = 1.0 - self.phi[m] * MC_DEF_INDEL;
    }

    fn q2p(&self, pl: i32) -> f64 {
        match usize::try_from(pl) {
            Ok(i) if i < 256 => self.q2p[i],
            _ => 10f64.powf(-f64::from(pl) / 10.0),
        }
    }

    /// `cal_pdg`: sets P(D|g) and returns the rank of the least likely allele.
    fn cal_pdg(&mut self, pls: &[i32], pl_len: usize, n_allele: usize) -> i32 {
        let mut p = vec![0i64; n_allele];
        for j in 0..self.n {
            let pi = &pls[j * pl_len..];
            let pl = |k: usize| pi.get(k).copied().unwrap_or(0);
            self.pdg[j * 3] = self.q2p(pl(2));
            self.pdg[j * 3 + 1] = self.q2p(pl(1));
            self.pdg[j * 3 + 2] = self.q2p(pl(0));
            for (i, pp) in p.iter_mut().enumerate() {
                *pp += i64::from(pl((i + 1) * (i + 2) / 2 - 1));
            }
        }
        for (i, pp) in p.iter_mut().enumerate() {
            *pp = (*pp << 4) | i as i64;
        }
        for i in 1..n_allele {
            let mut j = i;
            while j > 0 && p[j] < p[j - 1] {
                p.swap(j, j - 1);
                j -= 1;
            }
        }
        let mut i = n_allele as i32 - 1;
        while i >= 0 {
            if p[i as usize] & 0xf == 0 {
                break;
            }
            i -= 1;
        }
        i
    }

    /// `bcf_p1_call_gt`: the genotype (`q<<2 | index`) of sample `k`; `f0` is the frequency
    /// of the reference allele.
    pub fn call_gt(&self, f0: f64, k: usize, site: SiteCall) -> i32 {
        let pdg = &self.pdg[k * 3..k * 3 + 3];
        let ploidy = self.ploidy.as_ref().map_or(2, |p| p[k]);
        let f3 = if ploidy == 2 {
            [(1. - f0) * (1. - f0), 2. * f0 * (1. - f0), f0 * f0]
        } else {
            [1. - f0, 0., f0]
        };
        let mut g = [0.0; 3];
        let mut sum = 0.0;
        for i in 0..3 {
            g[i] = pdg[i] * f3[i];
            sum += g[i];
        }
        let mut max = -1.0;
        let mut max_i = 0;
        for (i, gi) in g.iter_mut().enumerate() {
            *gi /= sum;
            if *gi > max {
                max = *gi;
                max_i = i as i32;
            }
        }
        if site == SiteCall::Reference {
            max_i = 2;
            max = g[2];
        }
        max = 1. - max;
        if max < 1e-308 {
            max = 1e-308;
        }
        let q = ((-4.343 * max.ln() + 0.499) as i32).min(99);
        (q << 2) | max_i
    }

    /// `mc_cal_y_core` with `beg` 0.
    fn cal_y(&mut self) {
        let n = self.n;
        let mut z0 = std::mem::take(&mut self.z);
        let mut z1 = std::mem::take(&mut self.zswap);
        z0.fill(0.0);
        z1.fill(0.0);
        z0[0] = 1.0;
        let mut last_min = 0usize;
        let mut last_max = 0usize;
        self.t = 0.0;
        let mut m_acc = 0usize;
        for j in 0..n {
            let ploidy = match &self.ploidy {
                None => 2,
                Some(p) => p[j],
            };
            let mut min = last_min;
            let mut max = last_max;
            let pdg = &self.pdg[j * 3..j * 3 + 3];
            while min < max && z0[min] < TINY {
                z0[min] = 0.0;
                z1[min] = 0.0;
                min += 1;
            }
            while max > min && z0[max] < TINY {
                z0[max] = 0.0;
                z1[max] = 0.0;
                max -= 1;
            }
            let m0 = m_acc as f64;
            m_acc += usize::from(ploidy);
            if ploidy == 1 {
                let p = [pdg[0], pdg[2]];
                max += 1;
                if min == 0 {
                    z1[0] = (m0 + 1.0) * p[0] * z0[0];
                }
                for k in min.max(1)..=max {
                    let kf = k as f64;
                    z1[k] = (m0 + 1.0 - kf) * p[0] * z0[k] + kf * p[1] * z0[k - 1];
                }
                let mut sum = 0.0;
                for &v in &z1[min..=max] {
                    sum += v;
                }
                self.t += (sum / m_acc as f64).ln();
                for v in &mut z1[min..=max] {
                    *v /= sum;
                }
                if min >= 1 {
                    z1[min - 1] = 0.0;
                }
                if j < n - 1 {
                    z1[max + 1] = 0.0;
                }
            } else if ploidy == 2 {
                let p = [pdg[0], 2.0 * pdg[1], pdg[2]];
                max += 2;
                if min == 0 {
                    z1[0] = (m0 + 1.0) * (m0 + 2.0) * p[0] * z0[0];
                }
                if min <= 1 {
                    z1[1] = m0 * (m0 + 1.0) * p[0] * z0[1] + (m0 + 1.0) * p[1] * z0[0];
                }
                for k in min.max(2)..=max {
                    let kf = k as f64;
                    z1[k] = (m0 - kf + 1.0) * (m0 - kf + 2.0) * p[0] * z0[k]
                        + kf * (m0 - kf + 2.0) * p[1] * z0[k - 1]
                        + kf * (kf - 1.0) * p[2] * z0[k - 2];
                }
                let mut sum = 0.0;
                for &v in &z1[min..=max] {
                    sum += v;
                }
                let mf = m_acc as f64;
                self.t += (sum / (mf * (mf - 1.0))).ln();
                for v in &mut z1[min..=max] {
                    *v /= sum;
                }
                if min >= 1 {
                    z1[min - 1] = 0.0;
                }
                if min >= 2 {
                    z1[min - 2] = 0.0;
                }
                if j < n - 1 {
                    z1[max + 1] = 0.0;
                    z1[max + 2] = 0.0;
                }
            } else {
                // Other ploidies leave the arrays as they are, swapped (as bcftools does).
            }
            std::mem::swap(&mut z0, &mut z1);
            last_min = min;
            last_max = max;
        }
        self.z = z0;
        self.zswap = z1;
    }

    /// `mc_cal_afs`: returns the expected frequency and the folded probabilities.
    fn cal_afs(&mut self) -> (f64, f64, f64) {
        let m = self.m;
        self.afs1.fill(0.0);
        self.cal_y();
        let phi = if self.kind == VariantKind::Indel {
            &self.phi_indel
        } else {
            &self.phi
        };
        let z = &self.z;
        let sum = F80::dot(&phi[..=m], &z[..=m]);
        for k in 0..=m {
            let v = F80::div_into(phi[k] * z[k], sum);
            self.afs1[k] = v;
            if v.is_nan() || v.is_infinite() {
                return (-1.0, 0.0, 0.0);
            }
        }
        // The folded sums: `(long double)(phi[k] + phi[M - k]) / 2. * z[k]` over all k, and
        // over 1..M.
        let folded = &mut self.scratch;
        folded.clear();
        folded.extend((0..=m).map(|k| phi[k] + phi[m - k]));
        let sum = F80::dot_half(folded, &z[..=m]);
        let sum2 = if m >= 2 {
            F80::dot_half(&folded[1..m], &z[1..m])
        } else {
            F80::dot_half(&[], &[])
        };
        let p_var_folded = sum2.ratio(sum);
        // bcftools reads `phi[k] + phi[M - k]` with `k` left at M by the loop above.
        let a = (phi[m] + phi[0]) / 2.0 * (z[m] + z[0]);
        let p_ref_folded = F80::div_into(a, sum);
        // `sum += k * afs1[k]`, the product in double.
        folded.clear();
        folded.extend(self.afs1.iter().enumerate().map(|(k, &a)| k as f64 * a));
        let f_exp = F80::sum(folded).ratio(F80::from_i64(m as i64));
        (f_exp, p_ref_folded, p_var_folded)
    }

    /// `bcf_p1_cal`: `None` for a site with fewer than two alleles.
    pub fn cal(
        &mut self,
        pls: &[i32],
        pl_len: usize,
        n_allele: usize,
        kind: VariantKind,
    ) -> Option<P1Result> {
        self.kind = kind;
        if n_allele < 2 {
            return None;
        }
        let m = self.m;
        let rank0 = self.cal_pdg(pls, pl_len, n_allele);
        let (f_exp, p_ref_folded, p_var_folded) = self.cal_afs();
        let p_ref = self.afs1[m];
        let p_var = F80::sum(&self.afs1[..m]).to_f64();
        let mut max = -1.0;
        let mut ac = -1i32;
        for (k, &z) in self.z.iter().enumerate() {
            if max < z {
                max = z;
                ac = k as i32;
            }
        }
        Some(P1Result {
            rank0,
            ac: m as i32 - ac,
            f_exp,
            p_ref_folded,
            p_ref,
            p_var_folded,
            p_var,
        })
    }
}
