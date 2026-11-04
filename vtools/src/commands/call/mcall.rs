//! The multiallelic and rare-variant caller, `bcftools call -m` (bcftools' `mcall.c`).
//!
//! The functions keep the names and the order of floating-point operations of the C code,
//! including its `float` intermediates, so that QUAL, GQ, and GP come out the same.

use anyhow::{Context, bail};
use pasta::variant::{
    Header, HeaderRecordType, Number, Record, ValueType,
    encoding::{INT32_MISSING, INT32_VECTOR_END, float_missing, float_vector_end},
    log,
};

use super::{
    CallFlags, NumberR, OutputTags, Shared, Tag,
    ccall::{c_int, test16},
    get_code,
    vcmp::Vcmp,
};

/// The tags `-m` reads and writes at every site.
#[derive(Clone, Copy, Debug)]
struct MTags {
    pl: Tag,
    qs: Tag,
    i16: Tag,
    ac: Tag,
    an: Tag,
    dp4: Tag,
    mq: Tag,
    gq: Tag,
    gp: Tag,
    pv4: Tag,
}

impl MTags {
    fn new(hdr: &Header) -> Self {
        Self {
            pl: Tag::new(hdr, b"PL"),
            qs: Tag::new(hdr, b"QS"),
            i16: Tag::new(hdr, b"I16"),
            ac: Tag::new(hdr, b"AC"),
            an: Tag::new(hdr, b"AN"),
            dp4: Tag::new(hdr, b"DP4"),
            mq: Tag::new(hdr, b"MQ"),
            gq: Tag::new(hdr, b"GQ"),
            gp: Tag::new(hdr, b"GP"),
            pv4: Tag::new(hdr, b"PV4"),
        }
    }
}

/// The most alleles a site can have to be called: allele sets are bit masks of a `u32`.
const MAX_ALLELES: usize = 32;

/// The likelihoods of the group whose call has the highest quality.
#[derive(Clone, Copy, Debug)]
struct SiteLikelihoods {
    ref_lk: f64,
    lk_sum: f64,
    max_qual: f64,
}

/// `bcf_gt_missing`.
const GT_MISSING: i32 = 0;

/// `bcf_gt_unphased`.
#[inline]
fn gt_unphased(a: i32) -> i32 {
    (a + 1) << 1
}

/// `bcf_gt_allele`.
#[inline]
fn gt_allele(g: i32) -> i32 {
    (g >> 1) - 1
}

/// `bcf_alleles2gt`.
#[inline]
fn alleles2gt(a: i32, b: i32) -> i32 {
    if a > b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

/// `bcf_gt2alleles`.
fn gt2alleles(igt: i32) -> (i32, i32) {
    let mut k = 0;
    let mut dk = 1;
    while k < igt {
        dk += 1;
        k += dk;
    }
    let b = dk - 1;
    let a = igt - k + b;
    (a, b)
}

/// log(exp(a)+exp(b)) (`logsumexp2`).
#[inline]
fn logsumexp2(a: f64, b: f64) -> f64 {
    if a > b {
        (1.0 + (b - a).exp()).ln() + a
    } else {
        (1.0 + (a - b).exp()).ln() + b
    }
}

/// A group of samples called together (bcftools' `smpl_grp_t`).
#[derive(Clone, Debug, Default)]
pub struct SampleGroup {
    pub smpl: Vec<u32>,
    qsum: Vec<f32>,
    ref_lk: f64,
    max_lk: f64,
    lk_sum: f64,
    als: u32,
    nals: u32,
}

impl SampleGroup {
    /// A group of the given samples.
    pub fn new(smpl: Vec<u32>) -> Self {
        Self {
            smpl,
            ..Self::default()
        }
    }
}

/// The settings and scratch buffers of `-m` calling (the mcall part of bcftools' `call_t`).
#[derive(Clone, Debug)]
pub struct MCall {
    groups: Vec<SampleGroup>,
    /// The FORMAT tag used with `-G` (QS or AD).
    group_tag: Option<Vec<u8>>,
    pl2p: [f64; 256],
    theta: f64,
    pls: Vec<i32>,
    pdg: Vec<f64>,
    gts: Vec<i32>,
    ac: Vec<i32>,
    gps: Vec<f32>,
    gqs: Vec<i32>,
    ads: Vec<i32>,
    als_map: Vec<i32>,
    pl_map: Vec<usize>,
    anno16: Vec<f32>,
    itmp: Vec<i32>,
    ftmp: Vec<f32>,
    alleles: Vec<u8>,
    vcmp: Vcmp,
    number_r: NumberR,
    tags: MTags,
    nals_new: usize,
    als_new: u32,
}

/// What a site turned out to be.
pub type CallResult = i32;

/// Stores a genotype's likelihood in a sample's GP values (`gps[igt] = lk`); a write past the
/// end of the buffer is dropped.
fn set_gp(gps: &mut [f32], igt: i32, lk: f64) {
    if let Some(g) = usize::try_from(igt).ok().and_then(|i| gps.get_mut(i)) {
        *g = lk as f32;
    }
}

impl MCall {
    /// bcftools' `mcall_init` after the header lines are added: the groups, the PL table,
    /// and the prior (`theta`, scaled by Watterson's factor and logged).
    pub fn new(
        hdr: &Header,
        groups: Vec<SampleGroup>,
        group_tag: Option<Vec<u8>>,
        theta: f64,
        ploidy: &[u8],
        nsmpl: usize,
    ) -> Self {
        let mut pl2p = [0.0; 256];
        for (i, p) in pl2p.iter_mut().enumerate() {
            *p = 10f64.powf(-(i as f64) / 10.0);
        }
        let mut theta = theta;
        if theta > 0.0 {
            let n: i64 = if ploidy.is_empty() {
                2 * nsmpl as i64
            } else {
                ploidy.iter().map(|&p| i64::from(p)).sum()
            };
            let mut am = 1.0;
            for i in 2..n {
                am += 1.0 / i as f64;
            }
            theta *= am;
            if theta >= 1.0 {
                eprintln!("The prior is too big (theta*aM={theta:.2}), going with 0.99");
                theta = 0.99;
            }
            theta = theta.ln();
        }
        Self {
            groups,
            group_tag,
            pl2p,
            theta,
            pls: Vec::new(),
            pdg: Vec::new(),
            gts: vec![0; nsmpl * 2],
            ac: Vec::new(),
            gps: Vec::new(),
            gqs: vec![0; nsmpl],
            ads: Vec::new(),
            als_map: Vec::new(),
            pl_map: Vec::new(),
            anno16: Vec::new(),
            itmp: Vec::new(),
            ftmp: Vec::new(),
            alleles: Vec::new(),
            vcmp: Vcmp::default(),
            number_r: NumberR::new(hdr),
            tags: MTags::new(hdr),
            nals_new: 0,
            als_new: 0,
        }
    }

    #[inline]
    fn pl2p(&self, pl: i32) -> f64 {
        if (0..256).contains(&pl) {
            self.pl2p[pl as usize]
        } else {
            10f64.powf(-f64::from(pl) / 10.0)
        }
    }

    /// `set_pdg`: converts PLs to normalized probabilities, filling missing values from the
    /// unseen allele (and changing the PLs as bcftools does).
    fn set_pdg(&mut self, n_smpl: usize, n_gt: usize, unseen: i32) {
        let nals = {
            let (a, b) = gt2alleles(n_gt as i32 - 1);
            debug_assert_eq!(a, b);
            b as usize + 1
        };
        self.pdg.clear();
        self.pdg.resize(n_smpl * n_gt, 0.0);
        for i in 0..n_smpl {
            let base = i * n_gt;
            let mut sum = 0.0;
            let mut j = 0;
            while j < n_gt {
                let pl = self.pls[base + j];
                if pl == INT32_VECTOR_END {
                    j = 0;
                    break;
                }
                if pl == INT32_MISSING {
                    break;
                }
                self.pdg[base + j] = self.pl2p(pl);
                sum += self.pdg[base + j];
                j += 1;
            }
            if j == 0 {
                j = n_gt;
                sum = n_gt as f64;
            } else if j < n_gt && unseen < 0 {
                sum = 0.0;
                for j in 0..n_gt {
                    if self.pls[base + j] == INT32_MISSING {
                        self.pls[base + j] = 255;
                    }
                    self.pdg[base + j] = self.pl2p(self.pls[base + j]);
                    sum += self.pdg[base + j];
                }
                j = n_gt;
            }
            if j < n_gt {
                let mut jj = 0;
                sum = 0.0;
                for ia in 0..nals as i32 {
                    for ib in 0..=ia {
                        if self.pls[base + jj] == INT32_MISSING {
                            let pl_at = |k: i32| self.pls.get(base + k as usize).copied();
                            let mut k = alleles2gt(ia, unseen);
                            if pl_at(k) == Some(INT32_MISSING) {
                                k = alleles2gt(ib, unseen);
                            }
                            if pl_at(k) == Some(INT32_MISSING) {
                                k = alleles2gt(unseen, unseen);
                            }
                            self.pls[base + jj] = match pl_at(k) {
                                Some(INT32_MISSING) | None => 255,
                                Some(v) => v,
                            };
                        }
                        self.pdg[base + jj] = self.pl2p(self.pls[base + jj]);
                        sum += self.pdg[base + jj];
                        jj += 1;
                    }
                }
            }
            let pdg = &mut self.pdg[base..base + n_gt];
            if sum == n_gt as f64 {
                pdg.fill(0.0);
            } else {
                for p in pdg {
                    *p /= sum;
                }
            }
        }
    }

    /// `init_allele_trimming_maps`.
    fn init_allele_trimming_maps(&mut self, nals_ori: usize, als_out: u32) {
        self.als_map.clear();
        let mut nout = 0;
        for i in 0..nals_ori {
            if als_out & (1 << i) != 0 {
                self.als_map.push(nout);
                nout += 1;
            } else {
                self.als_map.push(-1);
            }
        }
        self.pl_map.clear();
        let mut l = 0;
        for i in 0..nals_ori {
            for j in 0..=i {
                if als_out & (1 << i) != 0 && als_out & (1 << j) != 0 {
                    self.pl_map.push(l);
                }
                l += 1;
            }
        }
    }

    /// `mcall_find_best_alleles`: the most likely combination of up to three alleles.
    fn find_best_alleles(&mut self, nals: usize, ig: usize, ploidy: &[u8]) {
        let theta = self.theta;
        let ngts = nals * (nals + 1) / 2;
        let pdg = &self.pdg;
        let grp = &mut self.groups[ig];
        let mut max_als = 0u32;
        let mut ref_lk = f64::NEG_INFINITY;
        let mut max_lk = f64::NEG_INFINITY;
        let mut lk_sum = f64::NEG_INFINITY;
        let ploidy_of = |s: usize| ploidy.get(s).copied().unwrap_or(2);

        macro_rules! update_max {
            ($lk_tot:expr, $set:expr, $als:expr, $sum:expr) => {
                if max_lk < $lk_tot && $set {
                    max_lk = $lk_tot;
                    max_als = $als;
                }
                if $sum {
                    lk_sum = logsumexp2($lk_tot, lk_sum);
                }
            };
        }

        for ia in 0..nals {
            let mut lk_tot = 0.0;
            let mut lk_tot_set = false;
            let iaa = (ia + 1) * (ia + 2) / 2 - 1;
            for &s in &grp.smpl {
                let p = pdg[s as usize * ngts + iaa];
                if p != 0.0 {
                    lk_tot += p.ln();
                    lk_tot_set = true;
                }
            }
            if ia == 0 {
                ref_lk = lk_tot;
            } else {
                lk_tot += theta;
            }
            update_max!(lk_tot, lk_tot_set, 1 << ia, ia > 0 && lk_tot_set);
        }

        let q = &grp.qsum;
        if nals > 1 {
            for ia in 0..nals {
                if q[ia] == 0.0 {
                    continue;
                }
                let iaa = (ia + 1) * (ia + 2) / 2 - 1;
                for ib in 0..ia {
                    if q[ib] == 0.0 {
                        continue;
                    }
                    let mut lk_tot = 0.0;
                    let mut lk_tot_set = false;
                    let fa = f64::from(q[ia] / (q[ia] + q[ib]));
                    let fb = f64::from(q[ib] / (q[ia] + q[ib]));
                    let fa2 = fa * fa;
                    let fb2 = fb * fb;
                    let fab = 2.0 * fa * fb;
                    let ibb = (ib + 1) * (ib + 2) / 2 - 1;
                    let iab = iaa - ia + ib;
                    for &s in &grp.smpl {
                        let p = &pdg[s as usize * ngts..];
                        let pl = ploidy_of(s as usize);
                        let mut val = 0.0;
                        if pl == 2 {
                            val = fa2 * p[iaa] + fb2 * p[ibb] + fab * p[iab];
                        } else if pl == 1 {
                            val = fa * p[iaa] + fb * p[ibb];
                        }
                        if val != 0.0 {
                            lk_tot += val.ln();
                            lk_tot_set = true;
                        }
                    }
                    if ia != 0 {
                        lk_tot += theta;
                    }
                    if ib != 0 {
                        lk_tot += theta;
                    }
                    update_max!(lk_tot, lk_tot_set, 1 << ia | 1 << ib, lk_tot_set);
                }
            }
        }

        if nals > 2 {
            for ia in 0..nals {
                if q[ia] == 0.0 {
                    continue;
                }
                let iaa = (ia + 1) * (ia + 2) / 2 - 1;
                for ib in 0..ia {
                    if q[ib] == 0.0 {
                        continue;
                    }
                    let ibb = (ib + 1) * (ib + 2) / 2 - 1;
                    let iab = iaa - ia + ib;
                    for ic in 0..ib {
                        if q[ic] == 0.0 {
                            continue;
                        }
                        let mut lk_tot = 0.0;
                        let mut lk_tot_set = false;
                        let fa = f64::from(q[ia] / (q[ia] + q[ib] + q[ic]));
                        let fb = f64::from(q[ib] / (q[ia] + q[ib] + q[ic]));
                        let fc = f64::from(q[ic] / (q[ia] + q[ib] + q[ic]));
                        let fa2 = fa * fa;
                        let fb2 = fb * fb;
                        let fc2 = fc * fc;
                        let fab = 2.0 * fa * fb;
                        let fac = 2.0 * fa * fc;
                        let fbc = 2.0 * fb * fc;
                        let icc = (ic + 1) * (ic + 2) / 2 - 1;
                        let iac = iaa - ia + ic;
                        let ibc = ibb - ib + ic;
                        for &s in &grp.smpl {
                            let p = &pdg[s as usize * ngts..];
                            let pl = ploidy_of(s as usize);
                            let mut val = 0.0;
                            if pl == 2 {
                                val = fa2 * p[iaa]
                                    + fb2 * p[ibb]
                                    + fc2 * p[icc]
                                    + fab * p[iab]
                                    + fac * p[iac]
                                    + fbc * p[ibc];
                            } else if pl == 1 {
                                val = fa * p[iaa] + fb * p[ibb] + fc * p[icc];
                            }
                            if val != 0.0 {
                                lk_tot += val.ln();
                                lk_tot_set = true;
                            }
                        }
                        if ia != 0 {
                            lk_tot += theta;
                        }
                        if ib != 0 {
                            lk_tot += theta;
                        }
                        if ic != 0 {
                            lk_tot += theta;
                        }
                        update_max!(lk_tot, lk_tot_set, 1 << ia | 1 << ib | 1 << ic, lk_tot_set);
                    }
                }
            }
        }

        grp.max_lk = max_lk;
        grp.ref_lk = ref_lk;
        grp.lk_sum = lk_sum;
        grp.als = max_als;
        grp.nals = max_als.count_ones();
    }

    /// `mcall_set_ref_genotypes`: 0/0, or missing for samples without data.
    fn set_ref_genotypes(&mut self, nals_ori: usize, nsmpl: usize, ploidy: &[u8]) {
        let ngts = nals_ori * (nals_ori + 1) / 2;
        for a in &mut self.ac[..nals_ori] {
            *a = 0;
        }
        for s in 0..nsmpl {
            let pl = ploidy.get(s).copied().unwrap_or(2);
            let pdg = &self.pdg[s * ngts..(s + 1) * ngts];
            let gts = &mut self.gts[s * 2..s * 2 + 2];
            if pl == 0 || pdg.iter().all(|&p| p == 0.0) {
                gts[0] = GT_MISSING;
                gts[1] = if pl == 2 {
                    GT_MISSING
                } else {
                    INT32_VECTOR_END
                };
            } else {
                gts[0] = gt_unphased(0);
                gts[1] = if pl == 2 {
                    gt_unphased(0)
                } else {
                    INT32_VECTOR_END
                };
                self.ac[0] += i32::from(pl);
            }
        }
    }

    /// `mcall_call_genotypes`: the most likely genotypes of a group's samples, and their
    /// GQ and GP.
    fn call_genotypes(&mut self, nals_ori: usize, ig: usize, ploidy: &[u8], tags: OutputTags) {
        let ngts_ori = nals_ori * (nals_ori + 1) / 2;
        let ngts_new = self.nals_new * (self.nals_new + 1) / 2;
        let grp = &self.groups[ig];
        for &s in &grp.smpl {
            let s = s as usize;
            let pdg = &self.pdg[s * ngts_ori..(s + 1) * ngts_ori];
            // The rest of the flat buffer, as bcftools' `call->GPs + ismpl*ngts_new`: when the
            // unseen allele is among the best alleles but not on output, its genotypes index
            // past this sample's values (into the next sample's, which overwrites them, or
            // past the end, where nothing reads them).
            let gps = &mut self.gps[s * ngts_new..];
            let gts = &mut self.gts[s * 2..s * 2 + 2];
            let pl = ploidy.get(s).copied().unwrap_or(2);
            assert!(pl <= 2, "ploidy {pl} is not supported by -m");
            if pl == 0 {
                gts[0] = GT_MISSING;
                gts[1] = INT32_VECTOR_END;
                gps[0] = -1.0;
                continue;
            }
            if pdg.iter().all(|&p| p == 0.0) {
                gts[0] = GT_MISSING;
                gts[1] = if pl == 2 {
                    GT_MISSING
                } else {
                    INT32_VECTOR_END
                };
                gps[0] = -1.0;
                continue;
            }
            gts[0] = gt_unphased(0);
            gts[1] = if pl == 2 {
                gt_unphased(0)
            } else {
                INT32_VECTOR_END
            };
            let mut best_lk = 0.0;
            for ia in 0..nals_ori {
                if grp.als & (1 << ia) == 0 {
                    continue;
                }
                let iaa = (ia + 1) * (ia + 2) / 2 - 1;
                let qa = f64::from(grp.qsum[ia]);
                let lk = if pl == 2 {
                    pdg[iaa] * qa * qa
                } else {
                    pdg[iaa] * qa
                };
                let m = self.als_map[ia];
                let igt = if pl == 2 { alleles2gt(m, m) } else { m };
                set_gp(gps, igt, lk);
                if best_lk < lk {
                    best_lk = lk;
                    gts[0] = gt_unphased(m);
                }
            }
            if pl == 2 {
                gts[1] = gts[0];
                for ia in 0..nals_ori {
                    if grp.als & (1 << ia) == 0 {
                        continue;
                    }
                    let iaa = (ia + 1) * (ia + 2) / 2 - 1;
                    for ib in 0..ia {
                        if grp.als & (1 << ib) == 0 {
                            continue;
                        }
                        let iab = iaa - ia + ib;
                        let lk = 2.0 * pdg[iab] * f64::from(grp.qsum[ia]) * f64::from(grp.qsum[ib]);
                        let igt = alleles2gt(self.als_map[ia], self.als_map[ib]);
                        set_gp(gps, igt, lk);
                        if best_lk < lk {
                            best_lk = lk;
                            gts[0] = gt_unphased(self.als_map[ib]);
                            gts[1] = gt_unphased(self.als_map[ia]);
                        }
                    }
                }
            } else {
                gts[1] = INT32_VECTOR_END;
            }
            // The best genotype names an allele that is not on output (the unseen allele):
            // bcftools writes it anyway, pointing past the ALTs, with GP from values it never
            // set (NaN). That allele is missing here (`0/.`), and so are GQ and GP.
            let mut dropped = false;
            for g in gts.iter_mut() {
                if *g != INT32_VECTOR_END
                    && *g != GT_MISSING
                    && gt_allele(*g) as usize >= self.nals_new
                {
                    *g = GT_MISSING;
                    dropped = true;
                }
            }
            if dropped {
                for &g in gts.iter() {
                    if g != INT32_VECTOR_END && g != GT_MISSING {
                        self.ac[gt_allele(g) as usize] += 1;
                    }
                }
                gps[0] = -1.0;
                continue;
            }
            self.ac[gt_allele(gts[0]) as usize] += 1;
            if gts[1] != INT32_VECTOR_END {
                self.ac[gt_allele(gts[1]) as usize] += 1;
            }
        }
        if !(tags.gq || tags.gp) {
            return;
        }
        for &s in &grp.smpl {
            let s = s as usize;
            let gps = &mut self.gps[s * ngts_new..];
            // A haploid sample's count can exceed ngts_new in the same case (see above).
            let nmax = match ploidy.get(s).copied().unwrap_or(2) {
                2 => ngts_new,
                1 => grp.nals as usize,
                _ => 0,
            }
            .min(gps.len());
            let mut max = f64::from(gps[0]);
            if max < 0.0 || nmax == 0 {
                if tags.gp {
                    let mut nmax = nmax;
                    for g in &mut gps[..nmax] {
                        *g = 0.0;
                    }
                    if nmax == 0 {
                        gps[0] = float_missing();
                        nmax += 1;
                    }
                    if nmax < ngts_new {
                        gps[nmax] = float_vector_end();
                    }
                }
                self.gqs[s] = 0;
                continue;
            }
            let mut sum = f64::from(gps[0]);
            for &g in &gps[1..nmax] {
                let g = f64::from(g);
                if max < g {
                    max = g;
                }
                sum += g;
            }
            max = -4.34294 * (1.0 - max / sum).ln();
            self.gqs[s] = if max <= f64::from(i8::MAX) {
                max as i32
            } else {
                i32::from(i8::MAX)
            };
            if tags.gp {
                for g in &mut gps[..nmax] {
                    *g = (f64::from(*g) / sum) as f32;
                }
                for g in gps.iter_mut().take(ngts_new).skip(nmax) {
                    *g = float_vector_end();
                }
            }
        }
    }

    /// `mcall_trim_and_update_PLs`.
    fn trim_and_update_pls(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        nals_ori: usize,
        nals_new: usize,
        nsmpl: usize,
        ploidy: &[u8],
    ) -> anyhow::Result<()> {
        let npls_src = nals_ori * (nals_ori + 1) / 2;
        let npls_dst = nals_new * (nals_new + 1) / 2;
        self.itmp.clear();
        for s in 0..nsmpl {
            let src = &self.pls[s * npls_src..(s + 1) * npls_src];
            let start = self.itmp.len();
            match ploidy.get(s).copied().unwrap_or(2) {
                2 => self
                    .itmp
                    .extend(self.pl_map[..npls_dst].iter().map(|&k| src[k])),
                1 => {
                    for ia in 0..nals_new {
                        let isrc = (ia + 1) * (ia + 2) / 2 - 1;
                        self.itmp.push(src[self.pl_map[isrc]]);
                    }
                    if nals_new < npls_dst {
                        self.itmp.push(INT32_VECTOR_END);
                    }
                }
                _ => {
                    self.itmp.push(INT32_MISSING);
                    self.itmp.push(INT32_VECTOR_END);
                }
            }
            // bcftools trims in place and leaves older values after the vector end; they
            // are never read.
            self.itmp.resize(start + npls_dst, INT32_VECTOR_END);
        }
        self.tags
            .pl
            .update_format_i32(hdr, rec, &self.itmp)
            .context("Could not update FORMAT/PL")
    }

    /// `mcall_constrain_alleles`: replaces the alleles with those of the `-C alleles` target,
    /// estimating PLs of new alleles from the unseen allele.
    fn constrain_alleles(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        tgt: &[Vec<u8>],
        unseen: &mut i32,
    ) -> anyhow::Result<()> {
        assert!(!tgt.is_empty(), "targets have alleles");
        let n_allele = rec.allele_count();
        if !self.vcmp.set_ref(rec.allele(0), &tgt[0]) {
            bail!(
                "The reference alleles are not compatible at {}:{} .. {} vs {}",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1,
                String::from_utf8_lossy(&tgt[0]),
                String::from_utf8_lossy(rec.allele(0))
            );
        }
        let mut als: Vec<Vec<u8>> = Vec::with_capacity(tgt.len() + 1);
        self.als_map.clear();
        als.push(tgt[0].clone());
        self.als_map.push(0);
        let mut has_new = false;
        for t in &tgt[1..] {
            als.push(t.clone());
            let rec_alts: Vec<&[u8]> = (1..n_allele).map(|i| rec.allele(i)).collect();
            match self.vcmp.find_allele(&rec_alts, t) {
                Some(j) => self.als_map.push(j as i32 + 1),
                None => {
                    self.als_map.push(if *unseen >= 0 {
                        *unseen
                    } else {
                        n_allele as i32 - 1
                    });
                    has_new = true;
                }
            }
        }
        if *unseen != 0 {
            self.als_map.push(*unseen);
            als.push(rec.allele(*unseen as usize).to_vec());
        }
        let nals = als.len();
        if !has_new && nals == n_allele {
            return Ok(());
        }
        rec.update_alleles(hdr, &als);

        let mut pl_map: Vec<usize> = Vec::with_capacity(nals * (nals + 1) / 2);
        for i in 0..nals {
            for j in 0..=i {
                let (a, b) = (self.als_map[i], self.als_map[j]);
                pl_map.push(alleles2gt(a, b) as usize);
            }
        }

        let nsmpl = hdr.sample_count();
        let npls = get_code(rec.get_format_i32(hdr, b"PL", &mut self.pls));
        let npls_ori = if nsmpl > 0 && npls > 0 {
            npls as usize / nsmpl
        } else {
            0
        };
        self.itmp.clear();
        for i in 0..nsmpl {
            let ori = &self.pls[i * npls_ori..(i + 1) * npls_ori];
            let get = |k: i32| ori.get(k as usize).copied().unwrap_or(INT32_VECTOR_END);
            for (k, &k_ori) in pl_map.iter().enumerate() {
                let mut v = get(k_ori as i32);
                if v == INT32_MISSING && *unseen >= 0 {
                    let (ia, ib) = gt2alleles(k_ori as i32);
                    let mut k2 = alleles2gt(ia, *unseen);
                    if get(k2) == INT32_MISSING {
                        k2 = alleles2gt(ib, *unseen);
                    }
                    if get(k2) == INT32_MISSING {
                        k2 = alleles2gt(*unseen, *unseen);
                    }
                    v = get(k2);
                }
                if k == 0 && v == INT32_VECTOR_END {
                    v = INT32_MISSING;
                }
                self.itmp.push(v);
            }
        }
        if npls > 0 {
            rec.update_format_i32(hdr, b"PL", &self.itmp)
                .context("Could not update FORMAT/PL")?;
        }

        // QS
        let nqs = get_code(rec.get_info_f32(hdr, b"QS", &mut self.ftmp));
        let qs: Vec<f32> = (0..nals)
            .map(|i| {
                let m = self.als_map[i];
                if m < nqs { self.ftmp[m as usize] } else { 0.0 }
            })
            .collect();
        let _ = rec.update_info_f32(hdr, b"QS", &qs);

        self.constrain_number_r(hdr, rec, nals)?;
        if *unseen != 0 {
            *unseen = nals as i32 - 1;
        }
        Ok(())
    }

    /// Reorders the FORMAT Number=R fields to the `-C alleles` alleles (the end of
    /// `mcall_constrain_alleles`).
    fn constrain_number_r(
        &self,
        hdr: &Header,
        rec: &mut Record,
        nals: usize,
    ) -> anyhow::Result<()> {
        let nsmpl = hdr.sample_count();
        rec.unpack(pasta::variant::Unpack::FMT);
        let fmts: Vec<(Vec<u8>, ValueType, usize)> = rec
            .format_fields()
            .iter()
            .filter(|f| !f.is_removed())
            .filter_map(|f| {
                let id = f.id() as u32;
                if hdr.id_number(HeaderRecordType::Format, id) != Some(Number::R) {
                    return None;
                }
                let ty = hdr.id_type(HeaderRecordType::Format, id)?;
                Some((hdr.id_name(id)?.to_vec(), ty, f.n()))
            })
            .collect();
        for (key, ty, n) in fmts {
            match ty {
                ValueType::Integer => {
                    let mut ori = Vec::new();
                    if get_code(rec.get_format_i32(hdr, &key, &mut ori)) <= 0 {
                        continue;
                    }
                    let mut out = Vec::with_capacity(nsmpl * nals);
                    for j in 0..nsmpl {
                        for k in 0..nals {
                            let m = self.als_map[k] as usize;
                            out.push(ori.get(j * n + m).copied().unwrap_or(INT32_VECTOR_END));
                        }
                    }
                    rec.update_format_i32(hdr, &key, &out)
                        .context("Could not update a Number=R tag")?;
                }
                ValueType::Float => {
                    let mut ori = Vec::new();
                    if get_code(rec.get_format_f32(hdr, &key, &mut ori)) <= 0 {
                        continue;
                    }
                    let mut out = Vec::with_capacity(nsmpl * nals);
                    for j in 0..nsmpl {
                        for k in 0..nals {
                            let m = self.als_map[k] as usize;
                            out.push(ori.get(j * n + m).copied().unwrap_or(float_vector_end()));
                        }
                    }
                    rec.update_format_f32(hdr, &key, &out)
                        .context("Could not update a Number=R tag")?;
                }
                _ => bail!(
                    "Only Integer and Float Number=R tags are supported: FORMAT/{}",
                    String::from_utf8_lossy(&key)
                ),
            }
        }
        Ok(())
    }

    /// bcftools' `mcall`: calls one site. Returns 0 when the site is not a variant (or was
    /// skipped), 1 for a reference site, or the number of alleles of a variant site.
    pub fn call(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        ploidy: &[u8],
        unseen: usize,
        tgt_als: Option<&[Vec<u8>]>,
    ) -> anyhow::Result<CallResult> {
        let mut unseen = unseen as i32;
        if let Some(tgt) = tgt_als {
            self.constrain_alleles(hdr, rec, tgt, &mut unseen)?;
        }
        let nsmpl = hdr.sample_count();
        let nals_ori = rec.allele_count();
        self.read_likelihoods(hdr, rec, shared, nals_ori, nsmpl, unseen)?;

        if nals_ori > MAX_ALLELES {
            log::line(format_args!(
                "Too many alleles at {}:{}, skipping.",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1
            ));
            return Ok(0);
        }

        let lks = self.find_site_alleles(nals_ori, ploidy);
        let is_variant = self.als_new != 1;
        let flags = shared.flags;
        if flags.var_only && !is_variant {
            return Ok(0);
        }
        self.select_output_alleles(nals_ori, unseen, flags);
        let nals_new = self.nals_new;
        self.init_allele_trimming_maps(nals_ori, self.als_new);

        let mut n_ac = 0;
        if self.als_new == 1 {
            self.set_ref_genotypes(nals_ori, nsmpl, ploidy);
            self.tags.pl.remove_format(hdr, rec);
        } else if !is_variant {
            self.set_ref_genotypes(nals_ori, nsmpl, ploidy);
            self.trim_and_update_pls(hdr, rec, nals_ori, nals_new, nsmpl, ploidy)?;
        } else {
            let Some(n) = self.call_variant(hdr, rec, shared, nals_ori, nsmpl, ploidy)? else {
                return Ok(0);
            };
            n_ac = n;
        }
        if nals_ori != nals_new {
            self.number_r
                .trim(hdr, rec, &self.als_map, nals_ori, nals_new)?;
        }

        rec.set_qual(self.site_qual(n_ac, lks));
        if nals_new > 1 {
            self.tags
                .ac
                .update_info_i32(hdr, rec, &self.ac[1..nals_new])
                .context("Could not update INFO/AC")?;
        }
        n_ac += self.ac[0];
        self.tags
            .an
            .update_info_i32(hdr, rec, &[n_ac])
            .context("Could not update INFO/AN")?;

        self.update_alleles_and_genotypes(hdr, rec, nals_ori, nsmpl)?;
        self.update_i16_annotations(hdr, rec, shared.tags)?;

        Ok(if is_variant { nals_new as i32 } else { 1 })
    }

    /// The start of `mcall`: reads the PLs into genotype probabilities, and the allele
    /// frequencies of each group (QS or the `-G` tag, mixed with `-F` priors) normalized.
    #[inline]
    fn read_likelihoods(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        nals_ori: usize,
        nsmpl: usize,
        unseen: i32,
    ) -> anyhow::Result<()> {
        let ngts_ori = nals_ori * (nals_ori + 1) / 2;
        if self.ac.len() < nals_ori.max(1) {
            self.ac.resize(nals_ori.max(1), 0);
        }

        let npls = self.tags.pl.get_format_i32(hdr, rec, &mut self.pls);
        if npls != (nsmpl * ngts_ori) as i32 && npls != (nsmpl * nals_ori) as i32 {
            bail!("Wrong number of PL fields? nals={nals_ori} npl={npls}");
        }
        // Haploid-only PLs are read with the diploid stride, as bcftools does; pad them.
        self.pls.resize(nsmpl * ngts_ori, INT32_VECTOR_END);
        self.set_pdg(nsmpl, ngts_ori, unseen);

        if self.groups.len() == 1 {
            let qsum = &mut self.groups[0].qsum;
            let nqs = self.tags.qs.get_info_f32(hdr, rec, qsum);
            if nqs <= 0 {
                bail!(
                    "The QS annotation not present at {}:{}",
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1
                );
            }
            if (nqs as usize) < nals_ori {
                qsum.resize(nals_ori, 0.0);
            }
        } else {
            self.groups_qsum(hdr, rec, nals_ori, nsmpl)?;
        }

        if let (Some(an_tag), Some(ac_tag)) = (&shared.prior_an, &shared.prior_ac) {
            self.apply_prior(hdr, rec, an_tag, ac_tag, nals_ori)?;
        }

        for grp in &mut self.groups {
            let mut sum: f32 = 0.0;
            for &q in &grp.qsum[..nals_ori] {
                sum += q;
            }
            if sum != 0.0 {
                for q in &mut grp.qsum[..nals_ori] {
                    *q /= sum;
                }
            }
        }

        self.tags.qs.remove_info(hdr, rec);
        Ok(())
    }

    /// Finds the best alleles of each group, and their union in `als_new` (always with
    /// REF); returns the likelihoods of the group with the highest quality.
    #[inline]
    fn find_site_alleles(&mut self, nals_ori: usize, ploidy: &[u8]) -> SiteLikelihoods {
        self.als_new = 0;
        let mut lks = SiteLikelihoods {
            ref_lk: f64::NEG_INFINITY,
            lk_sum: f64::NEG_INFINITY,
            max_qual: f64::NEG_INFINITY,
        };
        for ig in 0..self.groups.len() {
            self.find_best_alleles(nals_ori, ig, ploidy);
            let grp = &self.groups[ig];
            self.als_new |= grp.als;
            if grp.max_lk == f64::NEG_INFINITY {
                continue;
            }
            let qual = -4.343 * (grp.ref_lk - logsumexp2(grp.lk_sum, grp.ref_lk));
            if lks.max_qual < qual {
                lks.max_qual = qual;
                lks.lk_sum = grp.lk_sum;
                lks.ref_lk = grp.ref_lk;
            }
        }
        if self.als_new & 1 == 0 {
            self.als_new |= 1;
        }
        lks
    }

    /// Adds the alleles `-A` and `-*` keep to `als_new`, and counts them in `nals_new`.
    #[inline]
    fn select_output_alleles(&mut self, nals_ori: usize, unseen: i32, flags: CallFlags) {
        self.nals_new = 0;
        for i in 0..nals_ori {
            if flags.keep_unseen && i as i32 == unseen && self.nals_new == 1 {
                self.nals_new += 1;
                self.als_new |= 1 << i;
            }
            if i > 0 && i as i32 == unseen {
                continue;
            }
            if flags.keep_alts {
                self.als_new |= 1 << i;
            }
            if self.als_new & (1 << i) != 0 {
                self.nals_new += 1;
            }
        }
    }

    /// Calls the genotypes of a variant site and sets FORMAT/GP, FORMAT/GQ, and the PLs;
    /// returns the number of ALT alleles called, or `None` when the site is skipped.
    #[inline]
    fn call_variant(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        nals_ori: usize,
        nsmpl: usize,
        ploidy: &[u8],
    ) -> anyhow::Result<Option<i32>> {
        let flags = shared.flags;
        let nals_new = self.nals_new;
        let ngts_new = nals_new * (nals_new + 1) / 2;
        // Room past the last sample for the genotypes of an unseen allele that is not
        // kept (see `call_genotypes`): bcftools' buffer, rounded up by `hts_expand`,
        // takes such writes, and a haploid sample reads them back.
        let len = ngts_new * nsmpl + nals_ori * (nals_ori + 1) / 2;
        if self.gps.len() < len {
            self.gps.resize(len, 0.0);
        }
        if self.ac.len() < nals_new {
            self.ac.resize(nals_new, 0);
        }
        for a in &mut self.ac[..nals_new] {
            *a = 0;
        }
        if flags.constrain_trio && nals_new > 4 {
            log::line(format_args!(
                "Too many alleles at {}:{}, skipping.",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1
            ));
            return Ok(None);
        }
        if shared.tags.gq || shared.tags.gp {
            self.gps[..nsmpl * ngts_new].fill(0.0);
            self.gqs.fill(0);
        }
        for ig in 0..self.groups.len() {
            if flags.constrain_trio {
                bail!("todo: constrained trio calling temporarily disabled");
            }
            self.call_genotypes(nals_ori, ig, ploidy, shared.tags);
        }
        let mut n_ac = 0;
        for &a in &self.ac[1..nals_new] {
            n_ac += a;
        }
        if n_ac == 0 && flags.var_only {
            return Ok(None);
        }
        if shared.tags.gp {
            self.tags
                .gp
                .update_format_f32(hdr, rec, &self.gps[..nsmpl * ngts_new])
                .context("Could not update FORMAT/GP")?;
        }
        if shared.tags.gq {
            self.tags
                .gq
                .update_format_i32(hdr, rec, &self.gqs)
                .context("Could not update FORMAT/GQ")?;
        }
        self.trim_and_update_pls(hdr, rec, nals_ori, nals_new, nsmpl, ploidy)?;
        Ok(Some(n_ac))
    }

    /// The site's QUAL: from the best group when ALT alleles were called, else from the
    /// likelihood of a variant, else from the prior when REF was called.
    #[inline]
    fn site_qual(&self, n_ac: i32, lks: SiteLikelihoods) -> Option<f32> {
        if n_ac != 0 {
            Some(lks.max_qual as f32)
        } else if lks.lk_sum != f64::NEG_INFINITY {
            Some((-4.343 * (lks.lk_sum - logsumexp2(lks.lk_sum, lks.ref_lk))) as f32)
        } else if self.ac[0] != 0 {
            let q = if self.theta != 0.0 {
                -4.343 * self.theta
            } else {
                0.0
            };
            Some(q as f32)
        } else {
            None
        }
    }

    /// Writes the kept alleles and the genotypes.
    #[inline]
    fn update_alleles_and_genotypes(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        nals_ori: usize,
        nsmpl: usize,
    ) -> anyhow::Result<()> {
        // The allele map keeps the order, so the kept alleles are listed in order; an index
        // past nals_new (the unseen allele among the best alleles but not kept) is not listed.
        self.alleles.clear();
        let kept = 0..self.nals_new as i32;
        for i in (0..nals_ori).filter(|&i| kept.contains(&self.als_map[i])) {
            if !self.alleles.is_empty() {
                self.alleles.push(b',');
            }
            self.alleles.extend_from_slice(rec.allele(i));
        }
        let alleles = std::mem::take(&mut self.alleles);
        rec.update_alleles_str(hdr, &alleles);
        self.alleles = alleles;
        rec.update_genotypes(hdr, &self.gts[..nsmpl * 2])
            .context("Could not update FORMAT/GT")?;
        Ok(())
    }

    /// Replaces INFO/I16 with DP4, MQ, and with `-a PV4` the bias tests.
    #[inline]
    fn update_i16_annotations(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        tags: OutputTags,
    ) -> anyhow::Result<()> {
        if self.tags.i16.get_info_f32(hdr, rec, &mut self.anno16) == 16 {
            let a = &self.anno16;
            // A missing I16 value converts to INT_MIN, written as missing (see c_int).
            let dp = [c_int(a[0]), c_int(a[1]), c_int(a[2]), c_int(a[3])];
            self.tags
                .dp4
                .update_info_i32(hdr, rec, &dp)
                .context("Could not update INFO/DP4")?;
            let dp4_sum = c_int(a[0] + a[1] + a[2] + a[3]);
            let mq = if dp4_sum > 0 {
                c_int((a[8] + a[10]) / dp4_sum as f32)
            } else {
                0
            };
            self.tags
                .mq
                .update_info_i32(hdr, rec, &[mq])
                .context("Could not update INFO/MQ")?;
            if tags.pv4
                && let Some(t) = test16(a)
                && t.is_tested
            {
                let p = [t.p[0] as f32, t.p[1] as f32, t.p[2] as f32, t.p[3] as f32];
                self.tags
                    .pv4
                    .update_info_f32(hdr, rec, &p)
                    .context("Could not update INFO/PV4")?;
            }
        }
        self.tags.i16.remove_info(hdr, rec);
        Ok(())
    }

    /// The allele frequencies of `-G` groups from FORMAT/QS or FORMAT/AD.
    fn groups_qsum(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        nals_ori: usize,
        nsmpl: usize,
    ) -> anyhow::Result<()> {
        let tag = self.group_tag.as_deref().unwrap_or(b"AD");
        let nad = get_code(rec.get_format_i32(hdr, tag, &mut self.ads));
        if nad < 1 {
            bail!(
                "Error: FORMAT/{} is required with the -G option, mpileup must be run with \"-a AD\" or \"-a QS\"",
                String::from_utf8_lossy(tag)
            );
        }
        let nad = nad as usize / nsmpl;
        for grp in &mut self.groups {
            grp.qsum.clear();
            grp.qsum.resize(nals_ori.max(nad), 0.0);
            for &s in &grp.smpl {
                let ptr = &self.ads[s as usize * nad..(s as usize + 1) * nad];
                let mut sum: f32 = 0.0;
                for &v in ptr {
                    if v == INT32_VECTOR_END {
                        break;
                    }
                    if v != INT32_MISSING {
                        sum += v as f32;
                    }
                }
                if sum != 0.0 {
                    for (j, &v) in ptr.iter().enumerate() {
                        if v == INT32_VECTOR_END {
                            break;
                        }
                        if v != INT32_MISSING {
                            grp.qsum[j] += v as f32 / sum;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `-F AN,AC`: mixes reference-panel allele frequencies into the QS estimates.
    fn apply_prior(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        an_tag: &[u8],
        ac_tag: &[u8],
        nals_ori: usize,
    ) -> anyhow::Result<()> {
        let mut vals = std::mem::take(&mut self.itmp);
        if get_code(rec.get_info_i32(hdr, an_tag, &mut vals)) != 1 {
            self.itmp = vals;
            return Ok(());
        }
        let an = vals[0];
        if an > 0 && get_code(rec.get_info_i32(hdr, ac_tag, &mut vals)) == nals_ori as i32 - 1 {
            let mut ac0 = an;
            for (i, &ac) in vals.iter().enumerate().take(nals_ori - 1) {
                if ac == INT32_VECTOR_END {
                    break;
                }
                if ac == INT32_MISSING {
                    continue;
                }
                ac0 -= ac;
                for grp in &mut self.groups {
                    let q = (f64::from(grp.qsum[i + 1]) + 0.5 * f64::from(ac))
                        / (f64::from(grp.smpl.len() as u32) + 0.5 * f64::from(an));
                    grp.qsum[i + 1] = q as f32;
                }
            }
            if ac0 < 0 {
                bail!(
                    "Incorrect {},{} values at {}:{}",
                    String::from_utf8_lossy(an_tag),
                    String::from_utf8_lossy(ac_tag),
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1
                );
            }
            for grp in &mut self.groups {
                let q = (f64::from(grp.qsum[0]) + 0.5 * f64::from(ac0))
                    / (f64::from(grp.smpl.len() as u32) + 0.5 * f64::from(an));
                grp.qsum[0] = q as f32;
            }
        }
        self.itmp = vals;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_indices_round_trip() {
        for a in 0..5 {
            for b in 0..=a {
                assert_eq!(gt2alleles(alleles2gt(a, b)), (b, a));
            }
        }
    }
}
