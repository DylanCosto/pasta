//! The consensus caller, `bcftools call -c` (bcftools' `ccall.c`, with `prob1.c` and
//! `em.c`), and the I16 bias tests shared with `-m`.

use anyhow::{Context, bail};
use pasta::variant::{Header, Record, encoding::INT32_MISSING, encoding::INT32_VECTOR_END};

use super::{
    NumberR, Shared, Tag, em,
    kfunc::{kf_betai, kt_fisher_exact_two},
    prob1::{P1Aux, P1Result, SiteCall, VariantKind},
};

/// bcftools' `anno16_t`: the bias tests of the I16 annotation.
#[derive(Clone, Copy, Debug, Default)]
pub struct Anno16 {
    /// P-values of strand, base quality, mapping quality, and tail distance bias.
    pub p: [f64; 4],
    /// Root-mean-square mapping quality.
    pub mq: i32,
    pub is_tested: bool,
}

/// C's conversion of a float or double to `int` as x86 performs it (`cvttss2si`,
/// `cvttsd2si`): truncation, with NaN and values out of range giving `INT_MIN`. bcftools
/// converts I16 values this way, so a missing I16 value (a NaN) gives `INT_MIN`, which HTSlib
/// writes as a missing integer.
pub(super) fn c_int(x: impl Into<f64>) -> i32 {
    let x: f64 = x.into();
    if x.is_nan() || !(-2_147_483_648.0..2_147_483_648.0).contains(&x) {
        i32::MIN
    } else {
        x as i32
    }
}

/// The t-test of `test16_core` (`ttest` in `ccall.c`).
///
/// `n1 + n2` wraps as C's `int` does on x86, for counts made of wrapped I16 values.
fn ttest(n1: i32, n2: i32, a: &[f32]) -> f64 {
    let n = n1.wrapping_add(n2);
    if n1 == 0 || n2 == 0 || n < 3 {
        return 1.0;
    }
    let (n1f, n2f) = (f64::from(n1), f64::from(n2));
    let v = f64::from(n - 2);
    let u1 = f64::from(a[0]) / n1f;
    let u2 = f64::from(a[2]) / n2f;
    if u1 <= u2 {
        return 1.0;
    }
    let t = (u1 - u2)
        / (((f64::from(a[1]) - n1f * u1 * u1) + (f64::from(a[3]) - n2f * u2 * u2)) / v
            * (1. / n1f + 1. / n2f))
            .sqrt();
    if t < 0.0 {
        1.0
    } else {
        0.5 * kf_betai(0.5 * v, 0.5, v / (v + t * t))
    }
}

/// bcftools' `test16`: `None` when the depth is zero. A missing value among the first four
/// gives a depth of `INT_MIN` and leaves the site untested; one among the mapping qualities
/// makes MQ `INT_MIN` (missing). The counts convert to `int` as bcftools converts them (see
/// [`c_int`]), so a count of 2^31 or more takes part in the tests as `INT_MIN`.
pub fn test16(anno: &[f32]) -> Option<Anno16> {
    let mut a = Anno16 {
        p: [1.0; 4],
        ..Anno16::default()
    };
    let depth = c_int(anno[0] + anno[1] + anno[2] + anno[3]);
    a.is_tested = anno[0] + anno[1] > 0.0 && anno[2] + anno[3] > 0.0;
    if depth == 0 {
        return None;
    }
    a.mq = c_int(f64::from((anno[9] + anno[11]) / depth as f32).sqrt() + 0.499);
    a.p[0] = kt_fisher_exact_two(
        c_int(anno[0]),
        c_int(anno[1]),
        c_int(anno[2]),
        c_int(anno[3]),
    );
    let n1 = c_int(anno[0] + anno[1]);
    let n2 = c_int(anno[2] + anno[3]);
    for i in 1..4 {
        a.p[i] = ttest(n1, n2, &anno[4 * i..4 * i + 4]);
    }
    Some(a)
}

/// The tags `-c` reads and writes at every site.
#[derive(Clone, Copy, Debug)]
struct CTags {
    pl: Tag,
    i16: Tag,
    qs: Tag,
    af1: Tag,
    g3: Tag,
    hwe: Tag,
    ac1: Tag,
    dp4: Tag,
    mq: Tag,
    fq: Tag,
    pv4: Tag,
    gq: Tag,
}

impl CTags {
    fn new(hdr: &Header) -> Self {
        Self {
            pl: Tag::new(hdr, b"PL"),
            i16: Tag::new(hdr, b"I16"),
            qs: Tag::new(hdr, b"QS"),
            af1: Tag::new(hdr, b"AF1"),
            g3: Tag::new(hdr, b"G3"),
            hwe: Tag::new(hdr, b"HWE"),
            ac1: Tag::new(hdr, b"AC1"),
            dp4: Tag::new(hdr, b"DP4"),
            mq: Tag::new(hdr, b"MQ"),
            fq: Tag::new(hdr, b"FQ"),
            pv4: Tag::new(hdr, b"PV4"),
            gq: Tag::new(hdr, b"GQ"),
        }
    }
}

/// The PL value `-m` gives a missing likelihood when the unseen allele's is not available.
const MISSING_PL: i32 = 255;

/// Copies the PLs of every sample (`stride` values each) for the likelihoods of `-c`.
///
/// bcftools converts the first three values of each sample through a 256-entry table, so a
/// missing or vector-end value there makes it read far outside the table and crash. Those
/// samples are read as `-m` reads them (bcftools' `set_pdg`): a sample whose first value is
/// missing, or whose values end early, has no data (all 0, which `mpileup` writes for a sample
/// without reads), and another missing value is [`MISSING_PL`]. Every other sample is copied
/// as it is.
fn pls_for_calling(pls: &[i32], stride: usize, out: &mut Vec<i32>) {
    out.clear();
    out.extend_from_slice(pls);
    if stride == 0 {
        return;
    }
    for sample in out.chunks_mut(stride) {
        let head = &sample[..sample.len().min(3)];
        if head.first() == Some(&INT32_MISSING) || head.contains(&INT32_VECTOR_END) {
            sample.fill(0);
        } else {
            for pl in sample.iter_mut().take(3) {
                if *pl == INT32_MISSING {
                    *pl = MISSING_PL;
                }
            }
        }
    }
}

/// The number of INFO/I16 values `-c` uses.
const N_I16: usize = 16;

/// What a batch of records called on a worker thread (`--threads`) did with the INFO/I16
/// values that `-c` gives a record without I16: those of the last record called with I16.
#[derive(Clone, Debug, Default)]
pub struct I16Batch {
    /// Whether a record used any of the values the batch started from.
    pub took_start: bool,
    /// The values the batch leaves, when a record wrote any. They depend on the values the
    /// batch started from only when `took_start` is set.
    pub last: Option<Vec<f32>>,
}

/// What a record did with the I16 values held by [`CCall`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum I16Read {
    /// It has no I16 (or one that cannot be read): it uses the values held.
    Absent,
    /// It has fewer than 16 values, which replace the first ones held.
    Partial,
    /// It has 16 values or more, which replace all those held.
    Full,
}

/// Where the I16 values held by [`CCall`] come from, in the current batch.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum I16Source {
    /// The values the batch started from, some of them replaced when `changed`; `taken`
    /// tells whether a record used any of them.
    Start { taken: bool, changed: bool },
    /// Values all written by records of the batch; `took_start` tells whether a record
    /// before them used the start values.
    Batch { took_start: bool },
}

impl I16Source {
    /// The source after a record does `read`.
    fn after(self, read: I16Read) -> Self {
        match (self, read) {
            (Self::Start { taken, .. }, I16Read::Full) => Self::Batch { took_start: taken },
            (Self::Start { changed, .. }, I16Read::Absent) => Self::Start {
                taken: true,
                changed,
            },
            (Self::Start { .. }, I16Read::Partial) => Self::Start {
                taken: true,
                changed: true,
            },
            (batch @ Self::Batch { .. }, _) => batch,
        }
    }
}

/// The settings and scratch buffers of `-c` calling.
#[derive(Clone, Debug)]
pub struct CCall {
    p1: P1Aux,
    pl2p: [f64; 256],
    pls: Vec<i32>,
    /// The PLs the likelihoods are computed from (see [`pls_for_calling`]).
    pls_calc: Vec<i32>,
    pdg: Vec<f64>,
    gts: Vec<i32>,
    gqs: Vec<i32>,
    /// The 16 I16 values of the last records called with I16 (bcftools' `anno16` buffer).
    anno16: Vec<f32>,
    anno16_source: I16Source,
    /// The I16 values of the record being called, before they go into `anno16`.
    i16_read: Vec<f32>,
    als: Vec<u8>,
    als_map: Vec<i32>,
    pls_out: Vec<i32>,
    number_r: NumberR,
    tags: CTags,
}

impl CCall {
    /// bcftools' `ccall_init` (after the header lines are added).
    pub fn new(hdr: &Header, nsmpl: usize, ploidy: &[u8]) -> Self {
        let mut pl2p = [0.0; 256];
        for (i, p) in pl2p.iter_mut().enumerate() {
            *p = 10f64.powf(-(i as f64) / 10.0);
        }
        Self {
            p1: P1Aux::new(
                nsmpl,
                if ploidy.is_empty() {
                    None
                } else {
                    Some(ploidy)
                },
            ),
            pl2p,
            pls: Vec::new(),
            pls_calc: Vec::new(),
            pdg: Vec::new(),
            gts: vec![0; nsmpl * 2],
            gqs: vec![0; nsmpl],
            anno16: vec![0.0; N_I16],
            anno16_source: I16Source::Start {
                taken: false,
                changed: false,
            },
            i16_read: Vec::new(),
            als: Vec::new(),
            als_map: Vec::new(),
            pls_out: Vec::new(),
            number_r: NumberR::new(hdr),
            tags: CTags::new(hdr),
        }
    }

    /// Starts a batch of records on a worker thread: records without I16 take `anno16` until
    /// one reads its own.
    pub fn start_batch(&mut self, anno16: &[f32]) {
        self.anno16.clear();
        self.anno16.extend_from_slice(anno16);
        self.anno16.resize(N_I16, 0.0);
        self.anno16_source = I16Source::Start {
            taken: false,
            changed: false,
        };
    }

    /// What the batch since [`Self::start_batch`] did with the I16 values.
    pub fn end_batch(&self) -> I16Batch {
        match self.anno16_source {
            I16Source::Start { taken, changed } => I16Batch {
                took_start: taken,
                last: changed.then(|| self.anno16.clone()),
            },
            I16Source::Batch { took_start } => I16Batch {
                took_start,
                last: Some(self.anno16.clone()),
            },
        }
    }

    /// Reads the record's INFO/I16 into the values held, as bcftools'
    /// `bcf_get_info_float(..., &call->anno16, ...)`: without I16 (or with one that cannot
    /// be read) the values are those of the previous records, and fewer than 16 values
    /// replace only the first ones. Values no record has set are zeros, where bcftools crashes
    /// or reads memory it never set.
    fn read_i16(&mut self, hdr: &Header, rec: &mut Record) {
        let read = if self.tags.i16.get_info_f32(hdr, rec, &mut self.i16_read) >= 0 {
            let n = self.i16_read.len().min(N_I16);
            self.anno16[..n].copy_from_slice(&self.i16_read[..n]);
            if n == N_I16 {
                I16Read::Full
            } else {
                I16Read::Partial
            }
        } else {
            I16Read::Absent
        };
        self.anno16_source = self.anno16_source.after(read);
    }

    fn pl2p(&self, pl: i32) -> f64 {
        match usize::try_from(pl) {
            Ok(i) if i < 256 => self.pl2p[i],
            _ => 10f64.powf(-f64::from(pl) / 10.0),
        }
    }

    /// bcftools' `ccall`: returns whether the site is a variant (0 also when it is skipped
    /// with `-v`).
    pub fn call(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        ploidy: &[u8],
        unseen: usize,
    ) -> anyhow::Result<i32> {
        let nsmpl = hdr.sample_count();
        let nals = rec.allele_count();
        let npls = self.tags.pl.get_format_i32(hdr, rec, &mut self.pls);
        let ngts = nals * (nals + 1) / 2;
        if npls != (nsmpl * ngts) as i32 && npls != (nsmpl * nals) as i32 {
            bail!("Wrong number of PL fields? nals={nals} npl={npls}");
        }
        let stride = (npls as usize).checked_div(nsmpl).unwrap_or(0);
        pls_for_calling(&self.pls, stride, &mut self.pls_calc);
        // set_pdg3: only the first two alleles, in reverse order.
        self.pdg.clear();
        for i in 0..nsmpl {
            let pl = |k: usize| self.pls_calc.get(i * ngts + k).copied().unwrap_or(0);
            self.pdg.push(self.pl2p(pl(2)));
            self.pdg.push(self.pl2p(pl(1)));
            self.pdg.push(self.pl2p(pl(0)));
        }
        let em = em::em1(&self.pdg, nals);
        let kind = if super::is_snp(rec) {
            VariantKind::Snp
        } else {
            VariantKind::Indel
        };
        let Some(pr) = self.p1.cal(&self.pls_calc, stride, nals, kind) else {
            bail!("Something is wrong");
        };
        if pr.p_ref >= shared.pref && shared.flags.var_only {
            return Ok(0);
        }
        self.update_bcf1(hdr, rec, shared, &pr, &em, ploidy, unseen, stride)
    }

    /// bcftools' `update_bcf1`: writes the call into the record.
    #[expect(
        clippy::too_many_arguments,
        reason = "the state of bcftools' update_bcf1"
    )]
    fn update_bcf1(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        pr: &P1Result,
        em: &[f64; 10],
        ploidy: &[u8],
        unseen: usize,
        pl_stride: usize,
    ) -> anyhow::Result<i32> {
        self.read_i16(hdr, rec);
        let test = test16(&self.anno16);
        let a = test.unwrap_or(Anno16 {
            p: [1.0; 4],
            ..Anno16::default()
        });
        let has_i16 = test.is_some();

        if em[0] >= 0.0 {
            self.tags
                .af1
                .update_info_f32(hdr, rec, &[(1.0 - em[0]) as f32])
                .context("Could not update INFO/AF1")?;
        }
        if em[4] >= 0.0 && em[4] <= 0.05 {
            let g3 = [em[3] as f32, em[2] as f32, em[1] as f32];
            self.tags
                .g3
                .update_info_f32(hdr, rec, &g3)
                .context("Could not update INFO/G3")?;
            self.tags
                .hwe
                .update_info_f32(hdr, rec, &[em[4] as f32])
                .context("Could not update INFO/HWE")?;
        }

        let site = if pr.p_ref < shared.pref {
            SiteCall::Variant
        } else {
            SiteCall::Reference
        };
        let is_var = site == SiteCall::Variant;
        let r = if is_var { pr.p_ref } else { pr.p_var } as f32;

        self.tags
            .ac1
            .update_info_i32(hdr, rec, &[pr.ac])
            .context("Could not update INFO/AC1")?;
        let dp = [
            c_int(self.anno16[0]),
            c_int(self.anno16[1]),
            c_int(self.anno16[2]),
            c_int(self.anno16[3]),
        ];
        self.tags
            .dp4
            .update_info_i32(hdr, rec, &dp)
            .context("Could not update INFO/DP4")?;
        self.tags
            .mq
            .update_info_i32(hdr, rec, &[a.mq])
            .context("Could not update INFO/MQ")?;

        let mut fq = if pr.p_ref_folded < 0.5 {
            (-4.343 * pr.p_ref_folded.ln()) as f32
        } else {
            (4.343 * pr.p_var_folded.ln()) as f32
        };
        fq = fq.clamp(-999.0, 999.0);
        self.tags
            .fq
            .update_info_f32(hdr, rec, &[fq])
            .context("Could not update INFO/FQ")?;

        if has_i16 && a.is_tested {
            let p = [a.p[0] as f32, a.p[1] as f32, a.p[2] as f32, a.p[3] as f32];
            self.tags
                .pv4
                .update_info_f32(hdr, rec, &p)
                .context("Could not update INFO/PV4")?;
        }
        self.tags.i16.remove_info(hdr, rec);
        self.tags.qs.remove_info(hdr, rec);

        let mut qual = if f64::from(r) < 1e-100 {
            999.0
        } else {
            (-4.343 * f64::from(r).ln()) as f32
        };
        if qual > 999.0 {
            qual = 999.0;
        }
        rec.set_qual(Some(qual));

        let nals_ori = rec.allele_count();
        let mut nals = if !is_var && !shared.flags.keep_alts {
            1
        } else if pr.rank0 < 2 {
            2
        } else {
            pr.rank0 as usize + 1
        };
        if shared.flags.keep_alts && unseen == nals - 1 {
            nals -= 1;
        }
        if nals < nals_ori {
            self.trim_alleles(hdr, rec, nals, ploidy, pl_stride)?;
        }
        self.set_genotypes(hdr, rec, shared, pr, site, ploidy)?;

        self.als_map.clear();
        self.als_map
            .extend((0..nals_ori).map(|i| if i < nals { i as i32 } else { -1 }));
        if nals_ori != nals {
            self.number_r
                .trim(hdr, rec, &self.als_map, nals_ori, nals)?;
        }
        Ok(i32::from(is_var))
    }

    /// Keeps the first `nals` alleles and their PLs (part of `update_bcf1`).
    fn trim_alleles(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        nals: usize,
        ploidy: &[u8],
        pl_stride: usize,
    ) -> anyhow::Result<()> {
        self.als.clear();
        for i in 0..nals {
            if i > 0 {
                self.als.push(b',');
            }
            self.als.extend_from_slice(rec.allele(i));
        }
        let als = std::mem::take(&mut self.als);
        rec.update_alleles_str(hdr, &als);
        self.als = als;

        let nsmpl = rec.sample_count();
        let npls_src = pl_stride;
        let npls_dst = nals * (nals + 1) / 2;
        let mut out = std::mem::take(&mut self.pls_out);
        out.clear();
        for s in 0..nsmpl {
            let src = &self.pls[s * npls_src..];
            let get = |k: usize| src.get(k).copied().unwrap_or(INT32_VECTOR_END);
            let start = out.len();
            if ploidy.get(s).is_none_or(|&p| p == 2) {
                out.extend((0..npls_dst).map(get));
            } else {
                for i in 0..nals {
                    out.push(get((i + 1) * (i + 2) / 2 - 1));
                }
            }
            out.resize(start + npls_dst, INT32_VECTOR_END);
        }
        let updated = self.tags.pl.update_format_i32(hdr, rec, &out);
        self.pls_out = out;
        updated.context("Could not update FORMAT/PL")?;
        Ok(())
    }

    /// Sets FORMAT/GT, and FORMAT/GQ with `-a GQ` (part of `update_bcf1`).
    fn set_genotypes(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        pr: &P1Result,
        site: SiteCall,
        ploidy: &[u8],
    ) -> anyhow::Result<()> {
        let nsmpl = rec.sample_count();
        for i in 0..nsmpl {
            let x = if site == SiteCall::Variant || shared.tags.gq {
                self.p1.call_gt(pr.f_exp, i, site)
            } else {
                2
            };
            let gt = x & 3;
            if ploidy.get(i).is_none_or(|&p| p == 2) {
                let (a, b) = match gt {
                    1 => (0, 1),
                    0 => (1, 1),
                    _ => (0, 0),
                };
                self.gts[2 * i] = (a + 1) << 1;
                self.gts[2 * i + 1] = (b + 1) << 1;
                self.gqs[i] = x >> 2;
            } else {
                self.gts[2 * i] = if gt == 0 { 2 << 1 } else { 1 << 1 };
                self.gts[2 * i + 1] = INT32_VECTOR_END;
                self.gqs[i] = INT32_MISSING;
            }
        }
        rec.update_genotypes(hdr, &self.gts[..nsmpl * 2])
            .context("Could not update FORMAT/GT")?;
        if shared.tags.gq {
            self.tags
                .gq
                .update_format_i32(hdr, rec, &self.gqs[..nsmpl])
                .context("Could not update FORMAT/GQ")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_missing_pls_as_mcall_does() {
        let m = INT32_MISSING;
        let e = INT32_VECTOR_END;
        let mut out = Vec::new();
        pls_for_calling(&[0, 3, 30, m, 5, 6, 0, m, 9, 4, e, e], 3, &mut out);
        assert_eq!(out, [0, 3, 30, 0, 0, 0, 0, 255, 9, 0, 0, 0]);
        // Values past the first three are left as they are.
        pls_for_calling(&[0, 3, 30, e, e, e, 0, 1, 2, m, m, m], 6, &mut out);
        assert_eq!(out, [0, 3, 30, e, e, e, 0, 1, 2, m, m, m]);
        pls_for_calling(&[], 0, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn missing_i16_values_convert_to_int_min() {
        let m = pasta::variant::encoding::float_missing();
        assert_eq!(c_int(m), INT32_MISSING);
        assert_eq!(c_int(3e9_f32), i32::MIN);
        assert_eq!(c_int(-7.9_f32), -7);
        assert_eq!(c_int(59.999_999_9_f64), 59);
        let mut i16 = [
            5.0, 6.0, 7.0, 8.0, 330.0, 9950.0, 525.0, 18500.0, 660.0, 39700.0, 750.0, 37600.0,
            220.0, 4500.0, 375.0, 9500.0,
        ];
        let a = test16(&i16).expect("the depth is not zero");
        assert_eq!(a.mq, 55);
        assert!(a.is_tested);
        // A missing count makes the depth INT_MIN: the site has MQ missing and is not tested.
        i16[2] = m;
        let a = test16(&i16).expect("the depth is INT_MIN");
        assert_eq!(a.mq, INT32_MISSING);
        assert!(!a.is_tested);
        // A missing mapping-quality sum makes only MQ missing.
        i16[2] = 7.0;
        i16[11] = m;
        let a = test16(&i16).expect("the depth is not zero");
        assert_eq!(a.mq, INT32_MISSING);
        assert!(a.is_tested);
        // A negative sum of squares has no square root: MQ is missing too.
        i16[11] = -1e6;
        let a = test16(&i16).expect("the depth is not zero");
        assert_eq!(a.mq, INT32_MISSING);
    }

    #[test]
    fn huge_i16_counts_take_part_in_the_tests_as_int_min() {
        let mut i16 = [
            3e9, 3e9, 3e9, 3e9, 400.0, 16000.0, 200.0, 8000.0, 600.0, 36000.0, 300.0, 18000.0,
            100.0, 2500.0, 50.0, 1250.0,
        ];
        // Four INT_MIN counts sum to 0 in C's `int`, which leaves nothing to test.
        let a = test16(&i16).expect("the depth is INT_MIN");
        assert_eq!(a.p, [1.0; 4]);
        assert_eq!(a.mq, INT32_MISSING);
        assert!(a.is_tested);
        // bcftools gives a strand-bias P-value of 0 for one or two INT_MIN counts.
        i16[2] = 5.0;
        i16[3] = 6.0;
        let a = test16(&i16).expect("the depth is INT_MIN");
        assert_eq!(a.p, [0.0, 1.0, 1.0, 1.0]);
    }
}
