//! Genotype dosages and the per-pair sums of bcftools' `process_line`.
//!
//! A genotype is reduced to a dosage bitmask (`dsg`): bit 0, 1, or 2 for zero, one, or two
//! alternate alleles, 0 for missing. A PL triple sets every bit whose likelihood is the
//! smallest, so a pair of samples matches when their masks share a bit. With
//! `--error-probability`, each sample also gets the negative log probabilities of the three
//! genotypes, and a pair adds the smallest sum of its two samples' values.
//!
//! The sums over sample pairs are the hot loop. They are written without branches on the
//! data (a missing genotype adds zero), which keeps every pair's sums bit-identical to
//! bcftools' (adding `+0.0` never changes a sum that starts at `+0.0`) while letting the
//! compiler vectorize them.

use std::ops::Range;

use pasta::variant::encoding::{INT32_MISSING, INT32_VECTOR_END};

/// Which FORMAT field a file's genotypes are read from (bcftools' `use_GT`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(super) enum Tag {
    /// FORMAT/PL (`use_GT` 0).
    Pl,
    /// FORMAT/GT (`use_GT` 1).
    Gt,
}

impl Tag {
    /// bcftools' `use_GT` value, the index of the `nused` counters.
    pub(super) fn index(self) -> usize {
        match self {
            Self::Pl => 0,
            Self::Gt => 1,
        }
    }

    /// The number of values per sample of a diploid biallelic site.
    pub(super) fn values(self) -> usize {
        match self {
            Self::Pl => 3,
            Self::Gt => 2,
        }
    }
}

/// The dosage bits of the homozygous genotypes (zero or two alternate alleles).
const HOM_DSG: u8 = 0b101;

/// `--homs-only`: which genotypes of the `-g` side are compared.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Genotypes {
    All,
    /// Only homozygous ones (bcftools' `hom_only`).
    HomsOnly,
}

impl Genotypes {
    /// Whether a sample with dosage bits `d` is left out of the comparison.
    #[inline]
    pub(super) fn skips(self, d: u8) -> bool {
        self == Self::HomsOnly && d & HOM_DSG == 0
    }
}

/// Whether the sums keep the HWE scores of matching genotypes (bcftools' `calc_hwe_prob`,
/// off with `--no-HWE-prob`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum HweSums {
    Off,
    On,
}

/// HTSlib's `bcf_gt_is_missing`.
#[inline]
fn gt_is_missing(v: i32) -> bool {
    v >> 1 == 0
}

/// HTSlib's `bcf_gt_allele`.
#[inline]
fn gt_allele(v: i32) -> i32 {
    (v >> 1) - 1
}

/// bcftools' `gt_to_dsg`: the dosage bit of a diploid genotype, 0 when missing.
#[inline]
fn gt_to_dsg(p: &[i32]) -> u8 {
    if gt_is_missing(p[0]) || gt_is_missing(p[1]) || p[1] == INT32_VECTOR_END {
        return 0;
    }
    let dsg = u8::from(gt_allele(p[0]) != 0) + u8::from(gt_allele(p[1]) != 0);
    1 << dsg
}

/// bcftools' `pl_to_dsg`: the dosages with the smallest PL, 0 when missing.
#[inline]
fn pl_to_dsg(p: &[i32]) -> u8 {
    if p[0] == INT32_MISSING || p[1] == INT32_MISSING || p[2] == INT32_MISSING {
        return 0;
    }
    if p[1] == INT32_VECTOR_END || p[2] == INT32_VECTOR_END {
        return 0;
    }
    let min = p[0].min(p[1]).min(p[2]);
    u8::from(p[0] == min) | (u8::from(p[1] == min) << 1) | (u8::from(p[2] == min) << 2)
}

/// The dosage of one sample's values.
#[inline]
pub(super) fn dsg(tag: Tag, p: &[i32]) -> u8 {
    match tag {
        Tag::Gt => gt_to_dsg(p),
        Tag::Pl => pl_to_dsg(p),
    }
}

/// The genotype probability tables of `--error-probability` (bcftools' `dsg2prob` and
/// `pl2prob`).
#[derive(Clone, Debug)]
pub(super) struct Probs {
    dsg2prob: [[f64; 3]; 8],
    pl2prob: [f64; 256],
}

impl Probs {
    /// Builds the tables for a phred-scaled genotyping error `gt_err` (bcftools'
    /// `init_data`).
    pub(super) fn new(gt_err: i32) -> Self {
        let mut dsg2prob = [[f64::INFINITY; 3]; 8];
        let eprob = 10f64.powf(-0.1 * f64::from(gt_err));
        let e = -eprob.ln();
        let e2 = -2.0 * eprob.ln();
        dsg2prob[1] = [0.0, e, e2];
        dsg2prob[2] = [e, 0.0, e];
        dsg2prob[4] = [e2, e, 0.0];
        let mut pl2prob = [0.0; 256];
        for (i, p) in pl2prob.iter_mut().enumerate() {
            *p = 10f64.powf(-0.1 * i as f64);
        }
        Self { dsg2prob, pl2prob }
    }

    /// bcftools' `gt_to_prob` and `pl_to_prob`: the dosage of one sample's values and, when
    /// not missing, the negative log probabilities of the three genotypes.
    #[inline]
    pub(super) fn dsg_prob(&self, tag: Tag, p: &[i32], prob: &mut [f64; 3]) -> u8 {
        match tag {
            Tag::Gt => {
                let d = gt_to_dsg(p);
                if d != 0 {
                    *prob = self.dsg2prob[usize::from(d)];
                }
                d
            }
            Tag::Pl => {
                let d = pl_to_dsg(p);
                if d != 0 {
                    let pl = |v: i32| match usize::try_from(v) {
                        Ok(i) if i < 255 => self.pl2prob[i],
                        _ => self.pl2prob[255],
                    };
                    let (a, b, c) = (pl(p[0]), pl(p[1]), pl(p[2]));
                    let sum = a + b + c;
                    *prob = [-(a / sum).ln(), -(b / sum).ln(), -(c / sum).ln()];
                }
                d
            }
        }
    }
}

/// The negative log HWE probabilities of each match bitmask at a site (bcftools'
/// `hwe_dsg`), from the alternate allele counts.
pub(super) fn hwe_dsg(ac: [i32; 2]) -> [f64; 8] {
    let af = if ac[0] + ac[1] != 0 {
        f64::from(ac[1]) / f64::from(ac[0] + ac[1])
    } else {
        1e-6
    };
    let hwe = [
        -((1.0 - af) * (1.0 - af)).ln(),
        -(2.0 * af * (1.0 - af)).ln(),
        -(af * af).ln(),
    ];
    let mut out = [0.0; 8];
    for (i, v) in out.iter_mut().enumerate().skip(1) {
        *v = f64::INFINITY;
        for (k, h) in hwe.iter().enumerate() {
            if (1 << k) & i != 0 && *v > *h {
                *v = *h;
            }
        }
    }
    out
}

/// The discordance of each pair: mismatch counts (`ndiff`), or scores with
/// `--error-probability` (`pdiff`).
#[derive(Debug)]
pub(super) enum Discordance {
    Count(Vec<u32>),
    Score(Vec<f64>),
}

/// The sums of every pair of samples (bcftools' `ndiff`/`pdiff`, `ncnt`, `hwe_prob`,
/// `nmatch`).
#[derive(Debug)]
pub(super) struct Sums {
    pub(super) disc: Discordance,
    pub(super) ncnt: Vec<u32>,
    /// `hwe_prob` and `nmatch`, without `--no-HWE-prob`.
    pub(super) hwe: Option<(Vec<f64>, Vec<u32>)>,
}

impl Sums {
    /// Empty sums of `npairs` pairs: scores with `--error-probability` (`probs`), else
    /// mismatch counts.
    pub(super) fn new(npairs: usize, probs: Option<&Probs>, hwe: HweSums) -> Self {
        Self {
            disc: if probs.is_some() {
                Discordance::Score(vec![0.0; npairs])
            } else {
                Discordance::Count(vec![0; npairs])
            },
            ncnt: vec![0; npairs],
            hwe: (hwe == HweSums::On).then(|| (vec![0.0; npairs], vec![0; npairs])),
        }
    }

    fn slices(&mut self) -> SumSlices<'_> {
        SumSlices {
            disc: match &mut self.disc {
                Discordance::Count(v) => DiscSlice::Count(v),
                Discordance::Score(v) => DiscSlice::Score(v),
            },
            ncnt: &mut self.ncnt,
            hwe: self
                .hwe
                .as_mut()
                .map(|(h, m)| (h.as_mut_slice(), m.as_mut_slice())),
        }
    }

    /// The number of matching genotypes of pair `i`, 0 without HWE sums.
    pub(super) fn nmatch(&self, i: usize) -> u32 {
        self.hwe.as_ref().map_or(0, |(_, m)| m[i])
    }

    /// The average HWE score of pair `i` (0 when nothing matched or without HWE sums).
    pub(super) fn hwe_avg(&self, i: usize) -> f64 {
        match &self.hwe {
            Some((h, m)) if m[i] != 0 => h[i] / f64::from(m[i]),
            _ => 0.0,
        }
    }
}

enum DiscSlice<'a> {
    Count(&'a mut [u32]),
    Score(&'a mut [f64]),
}

/// Mutable views of a range of pairs' sums.
struct SumSlices<'a> {
    disc: DiscSlice<'a>,
    ncnt: &'a mut [u32],
    hwe: Option<(&'a mut [f64], &'a mut [u32])>,
}

impl<'a> SumSlices<'a> {
    /// Splits off the first `n` pairs.
    fn split_at(self, n: usize) -> (Self, Self) {
        let (d0, d1) = match self.disc {
            DiscSlice::Count(v) => {
                let (a, b) = v.split_at_mut(n);
                (DiscSlice::Count(a), DiscSlice::Count(b))
            }
            DiscSlice::Score(v) => {
                let (a, b) = v.split_at_mut(n);
                (DiscSlice::Score(a), DiscSlice::Score(b))
            }
        };
        let (c0, c1) = self.ncnt.split_at_mut(n);
        let (h0, h1) = match self.hwe {
            Some((h, m)) => {
                let (ha, hb) = h.split_at_mut(n);
                let (ma, mb) = m.split_at_mut(n);
                (Some((ha, ma)), Some((hb, mb)))
            }
            None => (None, None),
        };
        (
            Self {
                disc: d0,
                ncnt: c0,
                hwe: h0,
            },
            Self {
                disc: d1,
                ncnt: c1,
                hwe: h1,
            },
        )
    }
}

/// The dosages (and probabilities) of one side of the comparison at a batch of sites.
#[derive(Debug, Default)]
pub(super) struct Side {
    n: usize,
    dsg: Vec<u8>,
    /// The three genotype scores, one array each, with `--error-probability`.
    prob: [Vec<f64>; 3],
}

impl Side {
    fn new(n: usize) -> Self {
        Self {
            n,
            ..Self::default()
        }
    }

    fn clear(&mut self) {
        self.dsg.clear();
        for p in &mut self.prob {
            p.clear();
        }
    }

    /// Appends one site's values: `get(k)` gives sample k's values.
    #[inline]
    fn push_site<'a>(
        &mut self,
        tag: Tag,
        probs: Option<&Probs>,
        genotypes: Genotypes,
        get: impl Fn(usize) -> &'a [i32],
    ) {
        match probs {
            None => {
                for k in 0..self.n {
                    let mut d = dsg(tag, get(k));
                    if genotypes.skips(d) {
                        d = 0;
                    }
                    self.dsg.push(d);
                }
            }
            Some(probs) => {
                for k in 0..self.n {
                    let mut prob = [0.0; 3];
                    let mut d = probs.dsg_prob(tag, get(k), &mut prob);
                    if genotypes.skips(d) {
                        d = 0;
                    }
                    self.dsg.push(d);
                    for (dst, v) in self.prob.iter_mut().zip(prob) {
                        dst.push(v);
                    }
                }
            }
        }
    }

    fn site_dsg(&self, s: usize) -> &[u8] {
        &self.dsg[s * self.n..(s + 1) * self.n]
    }

    fn site_prob(&self, s: usize) -> [&[f64]; 3] {
        let r = s * self.n..(s + 1) * self.n;
        [
            &self.prob[0][r.clone()],
            &self.prob[1][r.clone()],
            &self.prob[2][r],
        ]
    }
}

/// How the pairs of the sample matrix are laid out.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Layout {
    /// Every query sample against every earlier query sample (`cross_check`): row `i` has
    /// `i` pairs, starting at `i*(i-1)/2`.
    Triangle,
    /// Every row sample against every column sample: rows of `ncols` pairs.
    Rect { ncols: usize },
}

impl Layout {
    /// The index of the first pair of row `i`.
    fn row_start(self, i: usize) -> usize {
        match self {
            Self::Triangle => i * i.saturating_sub(1) / 2,
            Self::Rect { ncols } => i * ncols,
        }
    }

    /// The number of pairs of row `i`.
    fn row_len(self, i: usize) -> usize {
        match self {
            Self::Triangle => i,
            Self::Rect { ncols } => ncols,
        }
    }
}

/// Sites collected for the sample-matrix sums, which are added in batches so that the
/// pairs' sums stay in cache and can be shared out between threads.
///
/// The sums of a query sample against `-g` samples form a row, so that the innermost loop
/// runs over many pairs; when there are fewer `-g` samples than query samples, the matrix
/// is stored transposed (a row per `-g` sample).
#[derive(Debug)]
pub(super) struct Batch {
    layout: Layout,
    /// Whether the rows are the `-g` samples.
    transposed: bool,
    rows: Side,
    cols: Side,
    hwe: Vec<[f64; 8]>,
    len: usize,
}

/// The number of sites summed together.
pub(super) const BATCH_SITES: usize = 512;

/// About how many pairs' sums are updated together for a whole batch (to fit in L2 cache).
const CHUNK_PAIRS: usize = 8192;

/// The fewest pair updates (pairs times sites) a batch gives each thread: a smaller batch is
/// summed on fewer threads, as starting a thread costs about as much as this many updates.
const MIN_THREAD_WORK: usize = 1 << 20;

impl Batch {
    /// A batch for the cross-check of `n` query samples.
    pub(super) fn triangle(n: usize) -> Self {
        Self {
            layout: Layout::Triangle,
            transposed: false,
            rows: Side::new(n),
            cols: Side::new(0),
            hwe: Vec::new(),
            len: 0,
        }
    }

    /// A batch for `nqry` query samples against `ngt` genotype samples.
    pub(super) fn rect(nqry: usize, ngt: usize) -> Self {
        let transposed = ngt < nqry;
        let (nrows, ncols) = if transposed { (ngt, nqry) } else { (nqry, ngt) };
        Self {
            layout: Layout::Rect { ncols },
            transposed,
            rows: Side::new(nrows),
            cols: Side::new(ncols),
            hwe: Vec::new(),
            len: 0,
        }
    }

    /// The index in the sums of query sample `i` against genotype sample `j` (of a
    /// [`Self::rect`] batch).
    pub(super) fn pair_index(&self, i: usize, j: usize) -> usize {
        if self.transposed {
            j * self.cols.n + i
        } else {
            i * self.cols.n + j
        }
    }

    pub(super) fn len(&self) -> usize {
        self.len
    }

    /// Adds a site: the query and genotype values per sample, and the HWE table.
    #[expect(
        clippy::too_many_arguments,
        reason = "the two sides' values and settings of one site"
    )]
    #[inline]
    pub(super) fn push<'q, 'g>(
        &mut self,
        probs: Option<&Probs>,
        qry_tag: Tag,
        qry: impl Fn(usize) -> &'q [i32],
        gt_tag: Tag,
        gt: impl Fn(usize) -> &'g [i32],
        genotypes: Genotypes,
        hwe: [f64; 8],
    ) {
        let (qry_side, gt_side) = if self.transposed {
            (&mut self.cols, &mut self.rows)
        } else {
            (&mut self.rows, &mut self.cols)
        };
        qry_side.push_site(qry_tag, probs, Genotypes::All, qry);
        if let Layout::Rect { .. } = self.layout {
            gt_side.push_site(gt_tag, probs, genotypes, gt);
        }
        self.hwe.push(hwe);
        self.len += 1;
    }

    /// Adds the batch to the sums, on `threads` threads (1 for the calling thread only), and
    /// empties it.
    pub(super) fn flush(&mut self, sums: &mut Sums, threads: usize) {
        // A thread is started only for enough work to pay for starting it.
        let total = self.layout.row_start(self.rows.n);
        let worth = (total * self.len / MIN_THREAD_WORK).max(1);
        self.flush_parts(sums, threads.min(worth));
    }

    /// [`Self::flush`] with the rows split between `parts` threads (at most one per row).
    fn flush_parts(&mut self, sums: &mut Sums, parts: usize) {
        if self.len == 0 {
            return;
        }
        let nrows = self.rows.n;
        let total = self.layout.row_start(nrows);
        let parts = parts.clamp(1, nrows.max(1));
        let this = &*self;
        if parts == 1 {
            this.add_rows(0..nrows, sums.slices());
        } else {
            // Rows split so that each thread gets about the same number of pairs.
            let mut bounds = vec![0];
            let mut row = 0;
            for p in 1..parts {
                let target = total * p / parts;
                while row < nrows && this.layout.row_start(row) < target {
                    row += 1;
                }
                bounds.push(row);
            }
            bounds.push(nrows);
            std::thread::scope(|scope| {
                let mut rest = sums.slices();
                for w in bounds.windows(2) {
                    let (a, b) = (w[0], w[1]);
                    let n = this.layout.row_start(b) - this.layout.row_start(a);
                    let (mine, others) = rest.split_at(n);
                    rest = others;
                    // Rows whose thread cannot be started are added here.
                    if a < b {
                        pasta::util::threads::spawn_scoped_or_run(scope, move || {
                            this.add_rows(a..b, mine);
                        });
                    }
                }
            });
        }
        self.rows.clear();
        self.cols.clear();
        self.hwe.clear();
        self.len = 0;
    }

    /// Adds every site of the batch to the pairs of `rows`, whose sums `out` holds, with the
    /// widest vector instructions the CPU has.
    fn add_rows(&self, rows: Range<usize>, out: SumSlices<'_>) {
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2, checked just above.
            unsafe { self.add_rows_avx2(rows, out) };
            return;
        }
        self.add_rows_portable(rows, out);
    }

    /// [`Self::add_rows`] compiled for AVX2 (four doubles per instruction).
    ///
    /// # Safety
    ///
    /// The CPU must support AVX2.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn add_rows_avx2(&self, rows: Range<usize>, out: SumSlices<'_>) {
        self.add_rows_portable(rows, out);
    }

    /// [`Self::add_rows`], a chunk of rows at a time.
    #[inline(always)]
    fn add_rows_portable(&self, rows: Range<usize>, mut out: SumSlices<'_>) {
        let mut i = rows.start;
        while i < rows.end {
            let mut j = i + 1;
            while j < rows.end
                && self.layout.row_start(j + 1) - self.layout.row_start(i) <= CHUNK_PAIRS
            {
                j += 1;
            }
            let n = self.layout.row_start(j) - self.layout.row_start(i);
            let (mut chunk, rest) = out.split_at(n);
            out = rest;
            for s in 0..self.len {
                self.add_site_rows(s, i..j, &mut chunk);
            }
            i = j;
        }
    }

    /// Adds site `s` to the pairs of `rows`.
    #[inline(always)]
    fn add_site_rows(&self, s: usize, rows: Range<usize>, out: &mut SumSlices<'_>) {
        let base = self.layout.row_start(rows.start);
        let qd = self.rows.site_dsg(s);
        let gt_side = match self.layout {
            Layout::Triangle => &self.rows,
            Layout::Rect { .. } => &self.cols,
        };
        let gd = gt_side.site_dsg(s);
        let hwe = &self.hwe[s];
        let none: [&[f64]; 3] = [&[], &[], &[]];
        let (qp, gp) = match out.disc {
            DiscSlice::Score(_) => (self.rows.site_prob(s), gt_side.site_prob(s)),
            DiscSlice::Count(_) => (none, none),
        };
        for i in rows {
            let q = qd[i];
            if q == 0 {
                continue;
            }
            let ng = self.layout.row_len(i);
            let r = self.layout.row_start(i) - base..self.layout.row_start(i) - base + ng;
            let gd = &gd[..ng];
            let ncnt = &mut out.ncnt[r.clone()];
            match &mut out.disc {
                DiscSlice::Count(ndiff) => {
                    let ndiff = &mut ndiff[r.clone()];
                    match &mut out.hwe {
                        Some((h, m)) => {
                            row_count_hwe(q, gd, hwe, ndiff, ncnt, &mut h[r.clone()], &mut m[r]);
                        }
                        None => row_count(q, gd, ndiff, ncnt),
                    }
                }
                DiscSlice::Score(pdiff) => {
                    let q3 = [qp[0][i], qp[1][i], qp[2][i]];
                    let g3 = [&gp[0][..ng], &gp[1][..ng], &gp[2][..ng]];
                    let pdiff = &mut pdiff[r.clone()];
                    match &mut out.hwe {
                        Some((h, m)) => row_score_hwe(
                            q,
                            q3,
                            gd,
                            g3,
                            hwe,
                            pdiff,
                            ncnt,
                            &mut h[r.clone()],
                            &mut m[r],
                        ),
                        None => row_score(q3, gd, g3, pdiff, ncnt),
                    }
                }
            }
        }
    }
}

/// Mismatch counts of one query sample (dosage `q`) against the genotype samples.
#[inline(always)]
fn row_count(q: u8, gd: &[u8], ndiff: &mut [u32], ncnt: &mut [u32]) {
    for ((&g, nd), nc) in gd.iter().zip(ndiff).zip(ncnt) {
        let valid = g != 0;
        *nd += u32::from(valid & (q & g == 0));
        *nc += u32::from(valid);
    }
}

/// [`row_count`] with the HWE sums of the matches.
#[inline(always)]
fn row_count_hwe(
    q: u8,
    gd: &[u8],
    hwe: &[f64; 8],
    ndiff: &mut [u32],
    ncnt: &mut [u32],
    hwe_prob: &mut [f64],
    nmatch: &mut [u32],
) {
    let n = gd.len();
    let (ndiff, ncnt, hwe_prob, nmatch) = (
        &mut ndiff[..n],
        &mut ncnt[..n],
        &mut hwe_prob[..n],
        &mut nmatch[..n],
    );
    if q.is_power_of_two() {
        // One dosage (always so with GT): a match adds hwe[q], anything else hwe[0] = 0.
        let hq = hwe[usize::from(q & 7)];
        for j in 0..n {
            let g = gd[j];
            let valid = g != 0;
            let hit = q & g != 0;
            ndiff[j] += u32::from(valid & !hit);
            ncnt[j] += u32::from(valid);
            hwe_prob[j] += if hit { hq } else { 0.0 };
            nmatch[j] += u32::from(hit);
        }
        return;
    }
    for j in 0..n {
        let g = gd[j];
        let valid = g != 0;
        let m = q & g & 7;
        ndiff[j] += u32::from(valid & (m == 0));
        ncnt[j] += u32::from(valid);
        // A missing genotype has no match bit: it adds hwe[0], which is 0.
        hwe_prob[j] += hwe[usize::from(m)];
        nmatch[j] += u32::from(m != 0);
    }
}

/// The smallest of the three genotype score sums, compared as bcftools compares them.
#[inline(always)]
pub(super) fn min_sum(q: [f64; 3], g0: f64, g1: f64, g2: f64) -> f64 {
    let mut m = q[0] + g0;
    let b = q[1] + g1;
    if m > b {
        m = b;
    }
    let c = q[2] + g2;
    if m > c {
        m = c;
    }
    m
}

/// Discordance scores of one query sample (scores `q`) against the genotype samples.
#[inline(always)]
fn row_score(q: [f64; 3], gd: &[u8], g: [&[f64]; 3], pdiff: &mut [f64], ncnt: &mut [u32]) {
    let n = gd.len();
    let (g0, g1, g2) = (&g[0][..n], &g[1][..n], &g[2][..n]);
    let (pdiff, ncnt) = (&mut pdiff[..n], &mut ncnt[..n]);
    for j in 0..n {
        let valid = gd[j] != 0;
        let m = min_sum(q, g0[j], g1[j], g2[j]);
        pdiff[j] += if valid { m } else { 0.0 };
        ncnt[j] += u32::from(valid);
    }
}

/// [`row_score`] with the HWE sums.
#[expect(clippy::too_many_arguments, reason = "the slices of one row's sums")]
#[inline(always)]
fn row_score_hwe(
    qd: u8,
    q: [f64; 3],
    gd: &[u8],
    g: [&[f64]; 3],
    hwe: &[f64; 8],
    pdiff: &mut [f64],
    ncnt: &mut [u32],
    hwe_prob: &mut [f64],
    nmatch: &mut [u32],
) {
    let n = gd.len();
    let (g0, g1, g2) = (&g[0][..n], &g[1][..n], &g[2][..n]);
    let (pdiff, ncnt, hwe_prob, nmatch) = (
        &mut pdiff[..n],
        &mut ncnt[..n],
        &mut hwe_prob[..n],
        &mut nmatch[..n],
    );
    if qd.is_power_of_two() {
        let hq = hwe[usize::from(qd & 7)];
        for j in 0..n {
            let gdj = gd[j];
            let valid = gdj != 0;
            let hit = qd & gdj != 0;
            let m = min_sum(q, g0[j], g1[j], g2[j]);
            pdiff[j] += if valid { m } else { 0.0 };
            ncnt[j] += u32::from(valid);
            hwe_prob[j] += if hit { hq } else { 0.0 };
            nmatch[j] += u32::from(hit);
        }
        return;
    }
    for j in 0..n {
        let gdj = gd[j];
        let valid = gdj != 0;
        let m = min_sum(q, g0[j], g1[j], g2[j]);
        pdiff[j] += if valid { m } else { 0.0 };
        ncnt[j] += u32::from(valid);
        let mt = qd & gdj & 7;
        hwe_prob[j] += hwe[usize::from(mt)];
        nmatch[j] += u32::from(mt != 0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const fn gt(a: i32) -> i32 {
        (a + 1) << 1
    }

    #[test]
    fn gt_dosages() {
        assert_eq!(gt_to_dsg(&[gt(0), gt(0)]), 1);
        assert_eq!(gt_to_dsg(&[gt(0), gt(1)]), 2);
        assert_eq!(gt_to_dsg(&[gt(1), gt(1) | 1]), 4);
        assert_eq!(gt_to_dsg(&[0, gt(1)]), 0);
        assert_eq!(gt_to_dsg(&[gt(1), INT32_VECTOR_END]), 0);
    }

    #[test]
    fn pl_dosages() {
        assert_eq!(pl_to_dsg(&[0, 10, 100]), 1);
        assert_eq!(pl_to_dsg(&[0, 0, 10]), 3);
        assert_eq!(pl_to_dsg(&[5, 5, 5]), 7);
        assert_eq!(pl_to_dsg(&[INT32_MISSING, 0, 0]), 0);
        assert_eq!(pl_to_dsg(&[0, 0, INT32_VECTOR_END]), 0);
    }

    #[test]
    fn hwe_takes_the_smallest_score() {
        let h = hwe_dsg([3, 1]);
        assert_eq!(h[0], 0.0);
        assert_eq!(h[1], -(0.75f64 * 0.75).ln());
        assert_eq!(h[3], h[1].min(h[2]));
    }

    #[test]
    fn triangle_rows_tile_the_pairs() {
        let l = Layout::Triangle;
        assert_eq!(l.row_start(0), 0);
        assert_eq!(l.row_start(1), 0);
        assert_eq!(l.row_start(2), 1);
        assert_eq!(l.row_start(4), 6);
    }

    /// The batched, threaded sums equal a direct pair-by-pair port of bcftools' loop.
    #[test]
    fn batch_matches_direct_sums() {
        let n = 37;
        let probs = Probs::new(40);
        let mut sites: Vec<Vec<i32>> = Vec::new();
        let mut seed = 12345u64;
        let mut next = || {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            (seed >> 33) as i32
        };
        for _ in 0..50 {
            let mut v = Vec::new();
            for _ in 0..n {
                let k = next() % 5;
                v.extend_from_slice(&match k {
                    0 => [0, 30, 300],
                    1 => [20, 0, 50],
                    2 => [255, 40, 0],
                    3 => [0, 0, 7],
                    _ => [INT32_MISSING, INT32_MISSING, INT32_MISSING],
                });
            }
            sites.push(v);
        }
        for score in [false, true] {
            let npairs = n * (n - 1) / 2;
            let mut want = Sums::new(npairs, score.then_some(&probs), HweSums::On);
            for v in &sites {
                let hwe = hwe_dsg([10, 3]);
                let mut idx = 0;
                for i in 0..n {
                    let mut qp = [0.0; 3];
                    let qd = probs.dsg_prob(Tag::Pl, &v[i * 3..], &mut qp);
                    if qd == 0 {
                        idx += i;
                        continue;
                    }
                    for j in 0..i {
                        let mut gp = [0.0; 3];
                        let gd = probs.dsg_prob(Tag::Pl, &v[j * 3..], &mut gp);
                        if gd == 0 {
                            idx += 1;
                            continue;
                        }
                        let m = qd & gd;
                        match &mut want.disc {
                            Discordance::Count(d) => d[idx] += u32::from(m == 0),
                            Discordance::Score(d) => d[idx] += min_sum(qp, gp[0], gp[1], gp[2]),
                        }
                        let (h, nm) = want.hwe.as_mut().expect("hwe");
                        if m != 0 || score {
                            h[idx] += hwe[usize::from(m)];
                        }
                        nm[idx] += u32::from(m != 0);
                        want.ncnt[idx] += 1;
                        idx += 1;
                    }
                }
            }
            for threads in [1, 3] {
                let mut got = Sums::new(npairs, score.then_some(&probs), HweSums::On);
                let mut batch = Batch::triangle(n);
                for v in &sites {
                    batch.push(
                        score.then_some(&probs),
                        Tag::Pl,
                        |k| &v[k * 3..k * 3 + 3],
                        Tag::Pl,
                        |k| &v[k * 3..k * 3 + 3],
                        Genotypes::All,
                        hwe_dsg([10, 3]),
                    );
                    if batch.len() == 7 {
                        batch.flush_parts(&mut got, threads);
                    }
                }
                batch.flush_parts(&mut got, threads);
                assert_eq!(got.ncnt, want.ncnt);
                assert_eq!(got.hwe, want.hwe);
                match (&got.disc, &want.disc) {
                    (Discordance::Count(a), Discordance::Count(b)) => assert_eq!(a, b),
                    (Discordance::Score(a), Discordance::Score(b)) => assert_eq!(a, b),
                    _ => panic!("mismatched discordance kinds"),
                }
            }
        }
    }

    /// PL triples of `n` samples at `nsites` sites, some missing or tied.
    fn random_sites(n: usize, nsites: usize, mut seed: u64) -> Vec<Vec<i32>> {
        let mut next = || {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            (seed >> 33) as i32
        };
        (0..nsites)
            .map(|_| {
                (0..n)
                    .flat_map(|_| match next() % 6 {
                        0 => [0, 30, 300],
                        1 => [20, 0, 50],
                        2 => [255, 40, 0],
                        3 => [0, 0, 7],
                        4 => [9, 0, 0],
                        _ => [INT32_MISSING, INT32_MISSING, INT32_MISSING],
                    })
                    .collect()
            })
            .collect()
    }

    /// Query samples against genotype samples, stored either way round, equal the direct
    /// loop, with `-H` applied to the genotype side.
    #[test]
    fn rect_batch_matches_direct_sums() {
        let probs = Probs::new(25);
        for (nqry, ngt) in [(9, 4), (4, 9), (5, 5)] {
            let sites = random_sites(nqry + ngt, 40, 99);
            for score in [false, true] {
                for genotypes in [Genotypes::All, Genotypes::HomsOnly] {
                    let npairs = nqry * ngt;
                    let mut want = Sums::new(npairs, score.then_some(&probs), HweSums::On);
                    let hwe = hwe_dsg([7, 5]);
                    for v in &sites {
                        for i in 0..nqry {
                            let mut qp = [0.0; 3];
                            let qd = probs.dsg_prob(Tag::Pl, &v[i * 3..], &mut qp);
                            for j in 0..ngt {
                                let mut gp = [0.0; 3];
                                let k = (nqry + j) * 3;
                                let mut gd = probs.dsg_prob(Tag::Pl, &v[k..], &mut gp);
                                if genotypes == Genotypes::HomsOnly && gd & 5 == 0 {
                                    gd = 0;
                                }
                                if qd == 0 || gd == 0 {
                                    continue;
                                }
                                let idx = i * ngt + j;
                                let m = qd & gd;
                                match &mut want.disc {
                                    Discordance::Count(d) => d[idx] += u32::from(m == 0),
                                    Discordance::Score(d) => {
                                        d[idx] += min_sum(qp, gp[0], gp[1], gp[2]);
                                    }
                                }
                                let (h, nm) = want.hwe.as_mut().expect("hwe");
                                h[idx] += hwe[usize::from(m)];
                                nm[idx] += u32::from(m != 0);
                                want.ncnt[idx] += 1;
                            }
                        }
                    }
                    let mut got = Sums::new(npairs, score.then_some(&probs), HweSums::On);
                    let mut batch = Batch::rect(nqry, ngt);
                    for v in &sites {
                        batch.push(
                            score.then_some(&probs),
                            Tag::Pl,
                            |k| &v[k * 3..k * 3 + 3],
                            Tag::Pl,
                            |k| &v[(nqry + k) * 3..(nqry + k) * 3 + 3],
                            genotypes,
                            hwe,
                        );
                    }
                    batch.flush_parts(&mut got, 2);
                    for i in 0..nqry {
                        for j in 0..ngt {
                            let (a, b) = (batch.pair_index(i, j), i * ngt + j);
                            assert_eq!(got.ncnt[a], want.ncnt[b], "{nqry}x{ngt} pair {i},{j}");
                            assert_eq!(got.nmatch(a), want.nmatch(b));
                            assert_eq!(got.hwe_avg(a).to_bits(), want.hwe_avg(b).to_bits());
                            match (&got.disc, &want.disc) {
                                (Discordance::Count(x), Discordance::Count(y)) => {
                                    assert_eq!(x[a], y[b]);
                                }
                                (Discordance::Score(x), Discordance::Score(y)) => {
                                    assert_eq!(x[a].to_bits(), y[b].to_bits());
                                }
                                _ => panic!("mismatched discordance kinds"),
                            }
                        }
                    }
                }
            }
        }
    }
}
