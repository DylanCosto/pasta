//! The per-site and per-sample work of `roh`: allele frequencies, emission probabilities, and
//! the HMM runs that print the states (`vcfroh`, `process_line`, and `flush_viterbi` in
//! `vcfroh.c`).

use std::{fmt, io::Write, sync::LazyLock, thread};

use anyhow::bail;
use pasta::variant::{Header, Record, Unpack, encoding::INT8_VECTOR_END, number, record::GetError};

use super::{
    genmap::{self, AZ, HW, MapCursor, MapPoint, Transitions},
    hmm::{Hmm, Matrix, Snapshot},
};
use crate::commands::call::is_snp;

/// The number of precalculated transition matrix powers (`hmm_init(2, tprob, 10000)`).
const NTPROB: usize = 10000;

/// BCF type codes of FORMAT values.
const BCF_BT_INT8: u8 = 1;
const BCF_BT_INT16: u8 = 2;
const BCF_BT_INT32: u8 = 3;

/// Genotype likelihoods of phred-scaled values 0 to 255 (`args->pl2p`).
static PL2P: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut pl2p = [0.0; 256];
    for (i, p) in pl2p.iter_mut().enumerate() {
        *p = 10f64.powf(-(i as f64) / 10.0);
    }
    pl2p
});

/// Bytes of output collected before they are written.
const OUT_CHUNK: usize = 1 << 16;

/// The number of sites the flush workers run in one wave, which bounds the output they hold
/// (about 30 bytes a site).
const WAVE_SITES: usize = 1 << 20;

/// Where allele frequencies come from (in `process_line`'s order of preference).
#[derive(Clone, Debug)]
pub enum AfSource {
    /// `--AF-tag`: an INFO tag.
    Tag(String),
    /// `--AF-file`: the fourth column of the targets file.
    File,
    /// `--AF-dflt` alone: the default frequency.
    Default,
    /// `-e`: estimated from GT or PL of some samples (indices) or all (`None`).
    Estimate {
        from: EstimateFrom,
        samples: Option<Vec<usize>>,
    },
    /// INFO/AC and INFO/AN.
    AcAn,
}

/// The FORMAT field `-e` estimates allele frequencies from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EstimateFrom {
    Gt,
    Pl,
}

/// A PL field of a type other than integer, on which bcftools exits with status 1.
#[derive(Debug)]
pub struct PlTypeError {
    /// The line of `vcfroh.c` that bcftools names.
    line: u32,
    /// The BCF type code.
    ty: u8,
}

impl fmt::Display for PlTypeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Unknown format type for PL: vcfroh.c:{} .. fmt->type={}",
            self.line, self.ty
        )
    }
}

impl std::error::Error for PlTypeError {}

/// What `-O` asks for.
#[derive(Clone, Copy, Debug, Default)]
pub struct Output {
    /// Per-site states (`ST` lines).
    pub sites: bool,
    /// Autozygous regions (`RG` lines).
    pub regions: bool,
}

/// `-V`: Viterbi training, with its convergence threshold.
#[derive(Clone, Copy, Debug)]
pub struct Training {
    pub threshold: f64,
}

/// The settings of a run (the options part of `args_t`).
#[derive(Debug)]
pub struct Config {
    /// The target samples, as header indices.
    pub roh_samples: Vec<usize>,
    /// The number of `-e` samples, as reported.
    pub naf_samples: usize,
    pub af_source: AfSource,
    /// `--AF-dflt`, 0 when not given.
    pub af_dflt: f64,
    /// `-G`: the genotype likelihood of unseen genotypes; `None` to use PL.
    pub unseen_pl: Option<f64>,
    pub pl_id: Option<u32>,
    pub gt_id: Option<u32>,
    pub skip_homref: bool,
    pub snps_only: bool,
    pub include_noalt: bool,
    pub t2az: f64,
    pub t2hw: f64,
    pub training: Option<Training>,
    /// Sites per buffer, 0 for unlimited.
    pub nbuf_max: usize,
    pub nbuf_olap: usize,
    pub output: Output,
    /// `-m`: the genetic map file name, possibly with `{CHROM}`.
    pub genmap: Option<String>,
    pub rec_rate: f64,
}

/// Line counts for the summary.
#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub ntot: u64,
    pub nfiltered: u64,
    pub nno_af: u64,
    pub nno_alt: u64,
    pub nmultiallelic: u64,
    pub ndup: u64,
}

/// The `--AF-file` line overlapping the current record.
#[derive(Clone, Copy, Debug)]
pub struct AfTarget<'a> {
    pub fname: &'a str,
    pub line: &'a [u8],
    pub alleles: &'a [Vec<u8>],
}

/// An open autozygous region of a sample (`smpl_t.rg`).
#[derive(Clone, Copy, Debug, Default)]
struct Region {
    beg: u32,
    end: u32,
    nqual: u32,
    qual: f64,
    rid: i32,
    open: bool,
}

/// The HMM data of one sample (`smpl_t`).
#[derive(Clone, Debug, Default)]
struct Sample {
    /// Emission probabilities, HW then AZ for each site.
    eprob: Vec<[f64; 2]>,
    sites: Vec<u32>,
    igenmap: usize,
    nused: u64,
    /// With Viterbi training: each chromosome and the index of its first site.
    rids: Vec<(i32, usize)>,
    snapshot: Option<Snapshot<2>>,
    rg: Region,
}

/// What `process_line` found for the allele frequency.
enum Af {
    Value(f64),
    /// The frequency is unknown (`ret < 0`).
    Missing,
    /// The site is skipped without counting it.
    Skip,
}

/// The state of a run (`args_t`).
pub struct Roh {
    cfg: Config,
    hmm: Hmm<2>,
    samples: Vec<Sample>,
    /// The genetic map of the last sequence loaded; empty when its file was missing, or moved
    /// to `seq_maps`.
    genmap: Vec<MapPoint>,
    /// With Viterbi training and a genetic map: the map of each sequence, for the flush at the
    /// end.
    seq_maps: Vec<(i32, Vec<MapPoint>)>,
    skip_rid: i32,
    prev_rid: i32,
    prev_pos: i64,
    pub stats: Stats,
    afs: Vec<f32>,
    ints: Vec<i32>,
    threads: usize,
    /// The models and output buffers of the flush workers, one per thread.
    workers: Vec<Worker>,
    /// Per target sample: whether its buffer is full at the current site.
    full: Vec<bool>,
    out_buf: Vec<u8>,
}

/// The model and output buffer a flush worker owns while it runs one sample.
struct Worker {
    hmm: Hmm<2>,
    out: Vec<u8>,
}

/// C's `%.1f`, with glibc's spelling of values that are not finite.
fn push_f1(out: &mut Vec<u8>, x: f64) {
    if x.is_nan() {
        out.extend_from_slice(if x.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
    } else if x.is_infinite() {
        out.extend_from_slice(if x < 0.0 { b"-inf" } else { b"inf" });
    } else if x.is_sign_positive() && x < 1e8 {
        // Rounding x * 10 is exact unless it is within its rounding error of a tie, which
        // the exact formatter decides.
        let y = x * 10.0;
        let floor = y.floor();
        let frac = y - floor;
        if (frac - 0.5).abs() > 1e-6 {
            let n = floor as u64 + u64::from(frac > 0.5);
            number::push_uint(out, n / 10);
            out.push(b'.');
            out.push(b'0' + (n % 10) as u8);
        } else {
            // Writing to a Vec cannot fail.
            let _ = write!(out, "{x:.1}");
        }
    } else {
        let _ = write!(out, "{x:.1}");
    }
}

/// C's `%e`: one digit, six decimals, and an exponent of at least two digits.
pub fn fmt_e(x: f64) -> String {
    if x.is_nan() {
        return if x.is_sign_negative() { "-nan" } else { "nan" }.to_owned();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    let s = format!("{x:.6e}");
    match s.split_once('e') {
        Some((mantissa, exp)) => {
            let (sign, digits) = match exp.strip_prefix('-') {
                Some(d) => ('-', d),
                None => ('+', exp),
            };
            format!("{mantissa}e{sign}{digits:0>2}")
        }
        None => s,
    }
}

/// bcftools' `phred_score`.
fn phred_score(prob: f64) -> f64 {
    if prob == 0.0 {
        return 99.0;
    }
    let p = -4.3429 * compat_log(prob);
    if p > 99.0 { 99.0 } else { p }
}

/// `log` as bcftools links it (glibc's `log@GLIBC_2.2.5` wrapper): a negative argument gives a
/// positive NaN (printed `nan`), where the current `log` gives a negative one (`-nan`). A
/// posterior just above 1 (from an `--AF-dflt` above 1, say) reaches it.
fn compat_log(x: f64) -> f64 {
    if x < 0.0 { f64::NAN } else { x.ln() }
}

/// HTSlib's `bcf_gt_is_missing` on a raw int8 GT value.
#[inline]
fn gt_is_missing(v: i8) -> bool {
    (i32::from(v) >> 1) == 0
}

/// HTSlib's `bcf_gt_allele` on a raw int8 GT value.
#[inline]
fn gt_allele(v: i8) -> i32 {
    (i32::from(v) >> 1) - 1
}

/// Value `i` of an integer FORMAT array of BCF type `ty`.
#[inline]
fn int_value(data: &[u8], ty: u8, i: usize) -> i32 {
    match ty {
        BCF_BT_INT8 => i32::from(data[i] as i8),
        BCF_BT_INT16 => i32::from(i16::from_le_bytes([data[2 * i], data[2 * i + 1]])),
        _ => i32::from_le_bytes([
            data[4 * i],
            data[4 * i + 1],
            data[4 * i + 2],
            data[4 * i + 3],
        ]),
    }
}

/// A record's FORMAT field: its BCF type, values per sample, and data.
struct FmtData<'a> {
    ty: u8,
    n: usize,
    data: &'a [u8],
}

/// Finds a FORMAT field by dictionary ID (FORMAT must be unpacked).
fn fmt_data(rec: &Record, id: Option<u32>) -> Option<FmtData<'_>> {
    let id = i32::try_from(id?).ok()?;
    let f = rec.format_fields().iter().find(|f| f.id() == id)?;
    Some(FmtData {
        ty: f.type_code(),
        n: f.n(),
        data: rec.format_data(f),
    })
}

/// The diploid GTs of a record (`get_GT`), `None` when absent or not diploid.
fn get_gt(rec: &Record, id: Option<u32>) -> anyhow::Result<Option<&[u8]>> {
    let Some(f) = fmt_data(rec, id) else {
        return Ok(None);
    };
    if f.n != 2 {
        return Ok(None);
    }
    if f.ty != BCF_BT_INT8 {
        bail!("This is unexpected, GT type is {}", f.ty);
    }
    Ok(Some(f.data))
}

/// Checks the type of a PL field as the switches of `vcfroh.c` do, at its `line`.
fn check_pl_type(ty: u8, line: u32) -> Result<(), PlTypeError> {
    if matches!(ty, BCF_BT_INT8 | BCF_BT_INT16 | BCF_BT_INT32) {
        Ok(())
    } else {
        Err(PlTypeError { line, ty })
    }
}

impl Roh {
    pub fn new(cfg: Config, threads: usize) -> Self {
        let mut tprob: Matrix<2> = [[0.0; 2]; 2];
        tprob[HW][HW] = 1.0 - cfg.t2az;
        tprob[HW][AZ] = cfg.t2hw;
        tprob[AZ][HW] = cfg.t2az;
        tprob[AZ][AZ] = 1.0 - cfg.t2hw;
        let nsamples = cfg.roh_samples.len();
        Self {
            hmm: Hmm::new(&tprob, NTPROB),
            samples: vec![Sample::default(); nsamples],
            cfg,
            genmap: Vec::new(),
            seq_maps: Vec::new(),
            skip_rid: -1,
            prev_rid: -1,
            prev_pos: 0,
            stats: Stats::default(),
            afs: Vec::new(),
            ints: Vec::new(),
            threads,
            workers: Vec::new(),
            full: vec![false; nsamples],
            out_buf: Vec::with_capacity(OUT_CHUNK * 2),
        }
    }

    /// The fewest sites any target sample used (`nmin`).
    pub fn min_used(&self) -> u64 {
        self.samples.iter().map(|s| s.nused).min().unwrap_or(0)
    }

    /// Processes one record that passed the filters (`vcfroh`).
    ///
    /// # Errors
    ///
    /// Returns bcftools' errors for unsorted input, a bad genetic map or AF file, and output
    /// errors.
    pub fn site(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        target: Option<AfTarget<'_>>,
        out: &mut impl Write,
    ) -> anyhow::Result<()> {
        if rec.rid() == self.skip_rid {
            return Ok(());
        }
        rec.unpack(Unpack::STR);
        // Raw calls may have the symbolic unseen allele <*>.
        let mut ial = 0;
        let mut nalt = rec.allele_count() as i64 - 1;
        for i in 1..rec.allele_count() {
            let a = rec.allele(i);
            if a == b"<*>" || a == b"<NON_REF>" {
                nalt -= 1;
            } else if ial == 0 {
                ial = i;
            }
        }
        if nalt == 0 {
            self.stats.nno_alt += 1;
            if !self.cfg.include_noalt {
                return Ok(());
            }
        } else if nalt > 1 {
            self.stats.nmultiallelic += 1;
            return Ok(());
        }
        if self.cfg.snps_only && !is_snp(rec) {
            return Ok(());
        }

        let mut skip = false;
        if self.prev_rid != rec.rid() {
            if self.cfg.training.is_none() {
                self.flush_samples(hdr, out, Flush::NewSequence)?;
            }
            if self.prev_rid < 0 {
                // bcftools loads the first sequence's map twice.
                self.load_genmap(rec.chrom_safe(hdr))?;
            }
            self.prev_rid = rec.rid();
            self.prev_pos = rec.pos();
            skip = !self.load_genmap(rec.chrom_safe(hdr))?;
            // bcftools trains and decodes every sequence with the map loaded last, as it
            // flushes only at the end with -V; each sequence keeps its own map here.
            if !skip && self.cfg.training.is_some() && self.cfg.genmap.is_some() {
                let map = std::mem::take(&mut self.genmap);
                self.seq_maps.push((rec.rid(), map));
            }
        } else if self.prev_pos == rec.pos() {
            self.stats.ndup += 1;
            return Ok(());
        }
        if skip {
            eprintln!(
                "Skipping the sequence, no genmap for {}",
                String::from_utf8_lossy(rec.chrom_safe(hdr))
            );
            self.skip_rid = rec.rid();
            return Ok(());
        }
        if self.prev_pos > rec.pos() {
            bail!("The file is not sorted?!");
        }
        self.prev_rid = rec.rid();
        self.prev_pos = rec.pos();
        self.process_line(hdr, rec, ial, target, out)
    }

    /// Loads the genetic map of a sequence; `false` when its file is missing (`load_genmap`).
    fn load_genmap(&mut self, chr: &[u8]) -> anyhow::Result<bool> {
        let Some(pattern) = &self.cfg.genmap else {
            return Ok(true);
        };
        match genmap::load(pattern, chr)? {
            Some(map) => {
                self.genmap = map;
                Ok(true)
            }
            None => {
                self.genmap.clear();
                Ok(false)
            }
        }
    }

    /// Flushes every sample at the end of a sequence or of the input.
    ///
    /// # Errors
    ///
    /// Returns output errors.
    pub fn finish(&mut self, hdr: &Header, out: &mut impl Write) -> std::io::Result<()> {
        self.flush_samples(hdr, out, Flush::End)?;
        out.write_all(&self.out_buf)?;
        self.out_buf.clear();
        Ok(())
    }

    /// The allele frequency of a site, by `--AF-tag`, `--AF-file`, `--AF-dflt`, `-e`, or AC/AN.
    fn allele_frequency(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        ial: usize,
        target: Option<AfTarget<'_>>,
    ) -> anyhow::Result<Af> {
        Ok(match &self.cfg.af_source {
            // A missing value (bcftools uses its NaN) leaves the frequency unknown.
            AfSource::Tag(tag) => match rec.get_info_f32(hdr, tag.as_bytes(), &mut self.afs) {
                Ok(n) if n > 0 && ial >= 1 && ial <= n && !self.afs[ial - 1].is_nan() => {
                    Af::Value(f64::from(self.afs[ial - 1]))
                }
                Err(GetError::WrongType) => bail!(
                    "Type mismatch for INFO/{tag} tag at {}:{}",
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1
                ),
                _ => Af::Missing,
            },
            AfSource::File => match target {
                Some(t) => read_af(t, rec)?,
                None => Af::Missing,
            },
            AfSource::Default => Af::Value(self.cfg.af_dflt),
            AfSource::Estimate { from, samples } => {
                rec.unpack(Unpack::FMT);
                let nsamples = hdr.sample_count();
                if *from == EstimateFrom::Pl {
                    let Some(pl) = fmt_data(rec, self.cfg.pl_id) else {
                        return Ok(Af::Skip);
                    };
                    estimate_af_from_pl(&pl, ial, samples.as_deref(), nsamples, &PL2P)?
                } else {
                    let Some(gt) = get_gt(rec, self.cfg.gt_id)? else {
                        return Ok(Af::Skip);
                    };
                    estimate_af_from_gt(gt, samples.as_deref(), nsamples)
                }
            }
            AfSource::AcAn => {
                let mut ac = -1;
                let mut an = 0;
                if let Ok(1) = rec.get_info_i32(hdr, b"AN", &mut self.ints) {
                    an = self.ints[0];
                    if let Ok(n) = rec.get_info_i32(hdr, b"AC", &mut self.ints)
                        && n > 0
                    {
                        ac = self.ints[0];
                    }
                }
                if an <= 0 || ac < 0 {
                    Af::Missing
                } else {
                    Af::Value(f64::from(ac) / f64::from(an))
                }
            }
        })
    }

    /// Adds a site to the target samples (`process_line`).
    fn process_line(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        ial: usize,
        target: Option<AfTarget<'_>>,
        out: &mut impl Write,
    ) -> anyhow::Result<()> {
        let af = self.allele_frequency(hdr, rec, ial, target)?;
        let dflt = self.cfg.af_dflt;
        let alt_freq = match af {
            Af::Skip => return Ok(()),
            Af::Value(v) if dflt > 0.0 && v == 0.0 => dflt,
            Af::Missing if dflt > 0.0 => dflt,
            Af::Value(v) if v != 0.0 => v,
            Af::Value(_) | Af::Missing => {
                self.stats.nno_af += 1;
                return Ok(());
            }
        };

        rec.unpack(Unpack::FMT);
        let ira = ial * (ial + 1) / 2;
        let iaa = ira + ial;
        let likelihoods = match self.cfg.unseen_pl {
            Some(unseen) => {
                // bcftools would crash on a site without diploid GTs.
                let Some(gt) = get_gt(rec, self.cfg.gt_id)? else {
                    return Ok(());
                };
                Likelihoods::Gt { gt, unseen }
            }
            None => {
                let Some(pl) = fmt_data(rec, self.cfg.pl_id) else {
                    return Ok(());
                };
                if iaa >= pl.n {
                    return Ok(());
                }
                // bcftools checks the type in its loop over the samples.
                if !self.cfg.roh_samples.is_empty() {
                    check_pl_type(pl.ty, 950)?;
                }
                Likelihoods::Pl {
                    pl,
                    ira,
                    iaa,
                    pl2p: &PL2P,
                }
            }
        };

        let pos = rec.pos() as u32;
        let rid = rec.rid();
        let mut nfull = 0;
        for i in 0..self.cfg.roh_samples.len() {
            let ismpl = self.cfg.roh_samples[i];
            let Some(mut pdg) = likelihoods.pdg(ismpl) else {
                continue;
            };
            let sum = pdg[0] + pdg[1] + pdg[2];
            if sum == 0.0 {
                continue;
            }
            for p in &mut pdg {
                *p /= sum;
            }
            if self.cfg.skip_homref && pdg[0] > 0.99 {
                continue;
            }

            let smpl = &mut self.samples[i];
            smpl.nused += 1;
            // P(D|HW) and P(D|AZ).
            let f = alt_freq;
            let az = pdg[0] * (1.0 - f) + pdg[2] * f;
            let hw = pdg[0] * (1.0 - f) * (1.0 - f) + 2.0 * pdg[1] * (1.0 - f) * f + pdg[2] * f * f;
            smpl.eprob.push([hw, az]);
            smpl.sites.push(pos);

            if self.cfg.training.is_some() {
                if smpl.rids.last().is_none_or(|&(r, _)| r != rid) {
                    smpl.rids.push((rid, smpl.sites.len() - 1));
                }
            } else if self.cfg.nbuf_max > 0 && smpl.sites.len() >= self.cfg.nbuf_max {
                self.full[i] = true;
                nfull += 1;
            }
        }
        // The full buffers are flushed in sample order, as bcftools flushes each in its loop
        // (the samples' runs are independent).
        if nfull > 0 {
            self.flush_samples(hdr, out, Flush::Full)?;
        }
        Ok(())
    }

    /// Runs the HMM of some samples (`flush_viterbi`) and writes their output in sample order:
    /// every sample at the end of a sequence or of the input, or those whose buffers are full.
    /// With threads, the samples run on workers, each with its own copy of the model.
    fn flush_samples(
        &mut self,
        hdr: &Header,
        out: &mut impl Write,
        why: Flush,
    ) -> std::io::Result<()> {
        let selected = |full: &mut [bool], i: usize| match why {
            Flush::Full => std::mem::take(&mut full[i]),
            Flush::NewSequence | Flush::End => true,
        };
        let ctx = FlushCtx {
            hdr,
            rid: self.prev_rid,
            cfg: &self.cfg,
            transitions: transitions(&self.cfg, &self.genmap),
            seq_maps: &self.seq_maps,
        };
        let names = hdr.samples();
        let roh_samples = &self.cfg.roh_samples;
        let nworkers = self.threads.min(self.samples.len());
        if nworkers <= 1 {
            for (i, smpl) in self.samples.iter_mut().enumerate() {
                if !selected(&mut self.full, i) {
                    continue;
                }
                flush_sample(
                    &mut self.hmm,
                    smpl,
                    &names[roh_samples[i]],
                    &ctx,
                    &mut self.out_buf,
                );
                if self.out_buf.len() >= OUT_CHUNK {
                    out.write_all(&self.out_buf)?;
                    self.out_buf.clear();
                }
            }
        } else {
            out.write_all(&self.out_buf)?;
            self.out_buf.clear();
            while self.workers.len() < nworkers {
                self.workers.push(Worker {
                    hmm: self.hmm.clone(),
                    out: Vec::new(),
                });
            }
            let mut picked: Vec<(usize, &mut Sample)> = self
                .samples
                .iter_mut()
                .enumerate()
                .filter(|&(i, _)| selected(&mut self.full, i))
                .collect();
            // A wave gives each worker a run of consecutive samples, so that the outputs come
            // in sample order; small buffers (`-b`) are batched, long sequences run one sample
            // per worker so that the outputs held stay small.
            let nsites: usize = picked.iter().map(|(_, s)| s.sites.len()).sum();
            let per_sample = (nsites / picked.len().max(1)).max(1);
            let wave_len = nworkers * (WAVE_SITES / nworkers / per_sample).max(1);
            let mut done = Vec::with_capacity(nworkers);
            for wave in picked.chunks_mut(wave_len) {
                let part_len = wave.len().div_ceil(nworkers);
                let ctx = &ctx;
                let idle = self.workers.drain(..wave.len().div_ceil(part_len));
                thread::scope(|s| {
                    // Each worker owns its model and buffer while it runs, so that no two
                    // threads write to the same cache lines. A part whose thread cannot be
                    // started runs here.
                    let handles: Vec<_> = wave
                        .chunks_mut(part_len)
                        .zip(idle)
                        .map(|(part, mut w)| {
                            pasta::util::threads::spawn_scoped_or_run(s, move || {
                                w.out.clear();
                                for (i, smpl) in part {
                                    let name = &names[roh_samples[*i]];
                                    flush_sample(&mut w.hmm, smpl, name, ctx, &mut w.out);
                                }
                                w
                            })
                        })
                        .collect();
                    done.extend(
                        handles
                            .into_iter()
                            .map(|h| h.join().expect("a roh flush worker panicked")),
                    );
                });
                for w in &done {
                    out.write_all(&w.out)?;
                }
                self.workers.append(&mut done);
            }
        }
        if why == Flush::NewSequence {
            for smpl in &mut self.samples {
                if let Some(s) = &mut smpl.snapshot {
                    s.deactivate();
                }
            }
        }
        Ok(())
    }
}

/// The transitions for the loaded genetic map or the recombination rate.
fn transitions<'a>(cfg: &Config, genmap: &'a [MapPoint]) -> Transitions<'a> {
    if cfg.genmap.is_some() {
        Transitions::GenMap {
            cursor: MapCursor::new(genmap, 0),
            rec_rate: cfg.rec_rate,
        }
    } else if cfg.rec_rate > 0.0 {
        Transitions::RecRate(cfg.rec_rate)
    } else {
        Transitions::Constant
    }
}

/// Why samples are flushed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Flush {
    /// The buffers of the samples marked in `Roh::full` are full (`-b`).
    Full,
    /// A new sequence starts: every sample, then their snapshots are deactivated.
    NewSequence,
    /// The input ended: every sample.
    End,
}

/// The genotype likelihoods of a site: from PL, or made up from GT with `-G`.
enum Likelihoods<'a> {
    Gt {
        gt: &'a [u8],
        unseen: f64,
    },
    Pl {
        pl: FmtData<'a>,
        ira: usize,
        iaa: usize,
        pl2p: &'a [f64; 256],
    },
}

impl Likelihoods<'_> {
    /// P(D|G) of a sample for RR, RA, and AA, or `None` when missing or uninformative.
    #[inline]
    fn pdg(&self, ismpl: usize) -> Option<[f64; 3]> {
        match self {
            Self::Gt { gt, unseen } => {
                let (g0, g1) = (gt[2 * ismpl] as i8, gt[2 * ismpl + 1] as i8);
                // A haploid GT next to diploid ones ends with the vector end, which bcftools
                // reads as an allele: the call becomes a het. A haploid call carries no
                // information on autozygosity, so the sample is skipped, as it is without -G
                // (its PL has no AA value).
                if gt_is_missing(g0) || gt_is_missing(g1) || g1 == INT8_VECTOR_END {
                    return None;
                }
                let (a, b) = (gt_allele(g0), gt_allele(g1));
                let u = *unseen;
                Some(if a != b {
                    [u, 1.0 - 2.0 * u, u]
                } else if a == 0 {
                    [1.0 - u - u * u, u, u * u]
                } else {
                    [u * u, u, 1.0 - u - u * u]
                })
            }
            Self::Pl { pl, ira, iaa, pl2p } => {
                let base = pl.n * ismpl;
                let rr = int_value(pl.data, pl.ty, base);
                let ra = int_value(pl.data, pl.ty, base + ira);
                let aa = int_value(pl.data, pl.ty, base + iaa);
                if rr < 0 || ra < 0 || aa < 0 || (rr == ra && rr == aa) {
                    return None;
                }
                Some([
                    pl2p[rr.min(255) as usize],
                    pl2p[ra.min(255) as usize],
                    pl2p[aa.min(255) as usize],
                ])
            }
        }
    }
}

/// The frequency of the `--AF-file` line matching a record's alleles (`read_AF`).
fn read_af(t: AfTarget<'_>, rec: &Record) -> anyhow::Result<Af> {
    if t.alleles.len() < 2 {
        bail!(
            "Expected two comma-separated alleles (REF,ALT) in the third column of {}, found:\n\t{}",
            t.fname,
            String::from_utf8_lossy(t.line)
        );
    }
    if t.alleles.len() != rec.allele_count()
        || t.alleles.iter().zip(rec.alleles()).any(|(a, b)| a != b)
    {
        return Ok(Af::Missing);
    }
    let line = t.line;
    let mut i = 0;
    let mut tabs = 0;
    while i < line.len() && tabs < 3 {
        if line[i] == b'\t' {
            tabs += 1;
        }
        i += 1;
    }
    let (af, end) = number::strtod(line, i);
    if end < line.len() && !line[end].is_ascii_whitespace() {
        let next = line.get(i + 1);
        if line.get(i) == Some(&b'.') && next.is_none_or(u8::is_ascii_whitespace) {
            return Ok(Af::Missing);
        }
        bail!("Could not parse: [{}]", String::from_utf8_lossy(line));
    }
    if !(0.0..=1.0).contains(&af) {
        bail!("Could not parse AF: [{}]", String::from_utf8_lossy(line));
    }
    Ok(Af::Value(af))
}

/// The frequency of non-reference alleles in the GTs (`estimate_AF_from_GT`).
fn estimate_af_from_gt(gt: &[u8], samples: Option<&[usize]>, nsamples: usize) -> Af {
    let (mut nalt, mut nref) = (0u32, 0u32);
    let mut count = |i: usize| {
        let (a, b) = (gt[2 * i] as i8, gt[2 * i + 1] as i8);
        if gt_is_missing(a) || gt_is_missing(b) {
            return;
        }
        // A haploid GT next to diploid ones has one allele: bcftools counts the vector end
        // padding it as an ALT.
        let alleles: &[i8] = if b == INT8_VECTOR_END { &[a] } else { &[a, b] };
        for &v in alleles {
            if gt_allele(v) != 0 {
                nalt += 1;
            } else {
                nref += 1;
            }
        }
    };
    match samples {
        Some(list) => list.iter().for_each(|&i| count(i)),
        None => (0..nsamples).for_each(&mut count),
    }
    if nalt == 0 && nref == 0 {
        return Af::Missing;
    }
    Af::Value(f64::from(nalt) / f64::from(nalt + nref))
}

/// The expected frequency of allele `ial` from the PLs (`estimate_AF_from_PL`).
fn estimate_af_from_pl(
    pl: &FmtData<'_>,
    ial: usize,
    samples: Option<&[usize]>,
    nsamples: usize,
    pl2p: &[f64; 256],
) -> Result<Af, PlTypeError> {
    let ira = ial * (ial + 1) / 2;
    let iaa = ira + ial;
    if iaa >= pl.n {
        return Ok(Af::Missing);
    }
    check_pl_type(pl.ty, if samples.is_some() { 778 } else { 808 })?;
    let mut af = 0.0;
    let mut naf = 0u32;
    let mut add = |i: usize| {
        let base = pl.n * i;
        let rr = int_value(pl.data, pl.ty, base);
        let ra = int_value(pl.data, pl.ty, base + ira);
        let aa = int_value(pl.data, pl.ty, base + iaa);
        if rr < 0 || ra < 0 || aa < 0 || (rr == ra && rr == aa) {
            return;
        }
        let mut prob = [
            pl2p[rr.min(255) as usize],
            pl2p[ra.min(255) as usize],
            pl2p[aa.min(255) as usize],
        ];
        let mut norm = 0.0;
        for p in &prob {
            norm += p;
        }
        for p in &mut prob {
            *p /= norm;
        }
        af += 0.5 * prob[1] + prob[2];
        naf += 1;
    };
    match samples {
        Some(list) => list.iter().for_each(|&i| add(i)),
        None => (0..nsamples).for_each(&mut add),
    }
    if naf == 0 {
        return Ok(Af::Missing);
    }
    Ok(Af::Value(af / f64::from(naf)))
}

/// What one sample's HMM run needs besides the model and the sample.
struct FlushCtx<'a> {
    hdr: &'a Header,
    /// The sequence of the buffered sites (`args->prev_rid`).
    rid: i32,
    cfg: &'a Config,
    transitions: Transitions<'a>,
    /// With Viterbi training: the genetic map of each sequence.
    seq_maps: &'a [(i32, Vec<MapPoint>)],
}

impl<'a> FlushCtx<'a> {
    /// The transitions of a sequence's sites with Viterbi training, from its own genetic map.
    fn seq_transitions(&self, rid: i32) -> Transitions<'a> {
        let map = self
            .seq_maps
            .iter()
            .find(|(r, _)| *r == rid)
            .map_or(&[][..], |(_, m)| m.as_slice());
        transitions(self.cfg, map)
    }
}

/// Appends `NAME\tCHR\t` to a line.
fn push_prefix(out: &mut Vec<u8>, tag: &[u8], name: &[u8], chr: &[u8]) {
    out.extend_from_slice(tag);
    out.push(b'\t');
    out.extend_from_slice(name);
    out.push(b'\t');
    out.extend_from_slice(chr);
    out.push(b'\t');
}

/// Prints a 1-based position or a length. Positions are `uint32_t` as in bcftools, so a record
/// at position 0 prints as 0; bcftools prints them with `%d`, as negative numbers from 2^31
/// on.
fn push_pos(out: &mut Vec<u8>, v: u32) {
    number::push_uint(out, u64::from(v));
}

/// Prints a closed autozygous region.
fn push_region(out: &mut Vec<u8>, name: &[u8], hdr: &Header, rg: &Region) {
    let chr = hdr.contig_name(rg.rid).unwrap_or_default();
    push_prefix(out, b"RG", name, chr);
    push_pos(out, rg.beg.wrapping_add(1));
    out.push(b'\t');
    push_pos(out, rg.end.wrapping_add(1));
    out.push(b'\t');
    push_pos(out, rg.end.wrapping_sub(rg.beg).wrapping_add(1));
    out.push(b'\t');
    number::push_int(out, i64::from(rg.nqual));
    out.push(b'\t');
    push_f1(out, rg.qual / f64::from(rg.nqual));
    out.push(b'\n');
}

/// Runs the HMM over a sample's buffered sites and prints the results (`flush_viterbi`).
fn flush_sample(
    hmm: &mut Hmm<2>,
    smpl: &mut Sample,
    name: &[u8],
    ctx: &FlushCtx<'_>,
    out: &mut Vec<u8>,
) {
    if smpl.sites.is_empty() {
        return;
    }
    if let Some(training) = ctx.cfg.training {
        train_and_print(hmm, smpl, name, ctx, training.threshold, out);
        return;
    }
    let cfg = ctx.cfg;
    hmm.restore(smpl.snapshot.as_ref());
    let nsites = smpl.sites.len();
    let olap = cfg.nbuf_olap;
    let end = if cfg.nbuf_max > 0 && nsites >= cfg.nbuf_max && nsites > olap {
        nsites - olap
    } else {
        nsites
    };
    // The runs fill in the sample's own snapshot where the next buffer will start. bcftools
    // has one snapshot pointer for all samples, so a sample's last run could overwrite the
    // snapshot of the sample that last took one.
    let mut snapshot = None;
    if end < nsites {
        let s = smpl.snapshot.get_or_insert_with(Snapshot::new);
        s.take_at(smpl.sites[nsites - olap - 1]);
        snapshot = Some(s);
    }
    let mut tr = ctx.transitions;
    tr.set_map_index(smpl.igenmap);
    hmm.run_viterbi(&smpl.eprob, &smpl.sites, &mut tr, snapshot.as_deref_mut());
    hmm.run_fwd_bwd(&smpl.eprob, &smpl.sites, &mut tr, snapshot);

    let chr = ctx.hdr.contig_name(ctx.rid).unwrap_or_default();
    // A local copy: with threads, the samples of a worker's neighbors are next to this one.
    let mut rg = smpl.rg;
    for i in 0..end {
        let state = usize::from(hmm.viterbi_state(i) == AZ as u8);
        let qual = phred_score(1.0 - hmm.posterior(i, state));
        if cfg.output.sites {
            push_prefix(out, b"ST", name, chr);
            push_pos(out, smpl.sites[i].wrapping_add(1));
            out.push(b'\t');
            number::push_int(out, state as i64);
            out.push(b'\t');
            push_f1(out, qual);
            out.push(b'\n');
        }
        if cfg.output.regions {
            if (state == 1) != rg.open {
                if state == 0 {
                    // The region ends.
                    push_region(out, name, ctx.hdr, &rg);
                    rg.open = false;
                } else {
                    rg = Region {
                        beg: smpl.sites[i],
                        end: smpl.sites[i],
                        nqual: 1,
                        qual,
                        rid: ctx.rid,
                        open: true,
                    };
                }
            } else if state == 1 {
                rg.nqual += 1;
                rg.qual += qual;
                rg.end = smpl.sites[i];
            }
        }
    }

    if end < nsites {
        let keep = nsites - olap;
        smpl.sites.drain(..keep);
        smpl.eprob.drain(..keep);
        smpl.igenmap = tr.map_index();
    } else {
        smpl.sites.clear();
        smpl.eprob.clear();
        smpl.igenmap = 0;
        if rg.open {
            push_region(out, name, ctx.hdr, &rg);
            rg.open = false;
        }
    }
    smpl.rg = rg;
}

/// Viterbi training over all of a sample's sequences, then the states (`flush_viterbi` with
/// `-V`).
fn train_and_print(
    hmm: &mut Hmm<2>,
    smpl: &Sample,
    name: &[u8],
    ctx: &FlushCtx<'_>,
    threshold: f64,
    out: &mut Vec<u8>,
) {
    let cfg = ctx.cfg;
    let mut tprob = *hmm.tprob();
    tprob[HW][HW] = 1.0 - cfg.t2az;
    tprob[HW][AZ] = cfg.t2hw;
    tprob[AZ][HW] = cfg.t2az;
    tprob[AZ][AZ] = 1.0 - cfg.t2hw;
    hmm.set_tprob(&tprob, NTPROB);

    let chunks: Vec<(i32, std::ops::Range<usize>)> = smpl
        .rids
        .iter()
        .enumerate()
        .map(|(i, &(rid, off))| {
            let end = smpl.rids.get(i + 1).map_or(smpl.sites.len(), |&(_, o)| o);
            (rid, off..end)
        })
        .collect();
    let mut niter = 0;
    loop {
        let t2az_prev = hmm.tprob()[AZ][HW];
        let t2hw_prev = hmm.tprob()[HW][AZ];
        let mut tprob_new: Matrix<2> = [[0.0; 2]; 2];
        for (rid, r) in &chunks {
            let mut tr = ctx.seq_transitions(*rid);
            let t = hmm.run_baum_welch(&smpl.eprob[r.clone()], &smpl.sites[r.clone()], &mut tr);
            for (n, v) in tprob_new.iter_mut().flatten().zip(t.iter().flatten()) {
                *n += v;
            }
        }
        for v in tprob_new.iter_mut().flatten() {
            *v /= chunks.len() as f64;
        }
        hmm.set_tprob(&tprob_new, NTPROB);

        let deltaz = (tprob_new[AZ][HW] - t2az_prev).abs();
        let delthw = (tprob_new[HW][AZ] - t2hw_prev).abs();
        niter += 1;
        let _ = writeln!(
            out,
            "VT\t{}\t{niter}\t{}\t{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(name),
            fmt_e(deltaz),
            fmt_e(delthw),
            fmt_e(1.0 - tprob_new[HW][HW]),
            fmt_e(tprob_new[AZ][HW]),
            fmt_e(1.0 - tprob_new[AZ][AZ]),
            fmt_e(tprob_new[HW][AZ])
        );
        if !(deltaz > threshold || delthw > threshold) {
            break;
        }
    }

    for (rid, r) in &chunks {
        let mut tr = ctx.seq_transitions(*rid);
        let eprob = &smpl.eprob[r.clone()];
        let sites = &smpl.sites[r.clone()];
        hmm.run_viterbi(eprob, sites, &mut tr, None);
        hmm.run_fwd_bwd(eprob, sites, &mut tr, None);
        let chr = ctx.hdr.contig_name(*rid).unwrap_or_default();
        for (j, &site) in sites.iter().enumerate() {
            let state = usize::from(hmm.viterbi_state(j) == AZ as u8);
            push_prefix(out, b"ROH", name, chr);
            push_pos(out, site.wrapping_add(1));
            out.push(b'\t');
            number::push_int(out, state as i64);
            out.push(b'\t');
            push_f1(out, phred_score(1.0 - hmm.posterior(j, state)));
            out.push(b'\n');
        }
    }
}

#[cfg(test)]
mod tests {
    use pasta::variant::vcf;

    use super::*;

    #[test]
    fn c_number_formats() {
        assert_eq!(fmt_e(6.7e-8), "6.700000e-08");
        assert_eq!(fmt_e(0.0), "0.000000e+00");
        let mut v = Vec::new();
        push_f1(&mut v, 0.25);
        push_f1(&mut v, f64::NAN);
        assert_eq!(v, b"0.2nan");
    }

    #[test]
    fn fast_f1_matches_exact_formatting() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut values = vec![0.0, -0.0, 0.05, 0.15, 0.25, 0.75, 2.5, 99.0, 99.95, 1e-300];
        for _ in 0..200_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            values.push((state >> 11) as f64 / (1u64 << 53) as f64 * 120.0);
        }
        let mut got = Vec::new();
        for x in values {
            got.clear();
            push_f1(&mut got, x);
            assert_eq!(got, format!("{x:.1}").as_bytes(), "{x:e}");
        }
    }

    #[test]
    fn phred_caps_at_99() {
        assert_eq!(phred_score(0.0), 99.0);
        assert_eq!(phred_score(1e-30), 99.0);
        assert!((phred_score(0.1) - 4.3429 * 10f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn phred_of_negative_is_positive_nan() {
        let q = phred_score(-1e-9);
        assert!(q.is_nan() && q.is_sign_positive(), "{q}");
    }

    /// GT values as raw int8 BCF bytes: (allele + 1) << 1, and the vector end.
    const REF: u8 = 2;
    const ALT: u8 = 4;
    const END: u8 = INT8_VECTOR_END as u8;

    #[test]
    fn gt_likelihoods_skip_haploid_calls() {
        // A haploid 0 padded with the vector end next to a diploid 0/0.
        let gt = [REF, END, REF, REF];
        let l = Likelihoods::Gt {
            gt: &gt,
            unseen: 1e-3,
        };
        assert_eq!(l.pdg(0), None);
        assert!(l.pdg(1).is_some_and(|p| p[0] > 0.99));
    }

    #[test]
    fn gt_af_counts_one_allele_of_haploid_calls() {
        // 0, 0/0, 1, ./. (missing, not counted).
        let gt = [REF, END, REF, REF, ALT, END, 0, 0];
        let Af::Value(af) = estimate_af_from_gt(&gt, None, 4) else {
            panic!("no AF");
        };
        assert_eq!(af, 1.0 / 4.0);
        let Af::Value(af) = estimate_af_from_gt(&gt, Some(&[0, 1]), 4) else {
            panic!("no AF");
        };
        assert_eq!(af, 0.0);
    }

    /// Runs `roh -G30 --AF-dflt 0.4 -O s` (changed by `edit`) over records given as
    /// `CHROM POS GT...` (an A>C SNP on contig 1 or 2) and returns the output.
    fn run_roh(samples: &[&str], recs: &[String], edit: impl FnOnce(&mut Config)) -> String {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        for s in samples {
            text += &format!("\t{s}");
        }
        text.push('\n');
        let mut hdr = Header::parse(text.as_bytes()).expect("the header parses");
        let mut cfg = Config {
            roh_samples: (0..samples.len()).collect(),
            naf_samples: 0,
            af_source: AfSource::Default,
            af_dflt: 0.4,
            unseen_pl: Some(1e-3),
            pl_id: None,
            gt_id: hdr.id(b"GT"),
            skip_homref: false,
            snps_only: false,
            include_noalt: false,
            t2az: 6.7e-8,
            t2hw: 5e-9,
            training: None,
            nbuf_max: 0,
            nbuf_olap: 0,
            output: Output {
                sites: true,
                regions: false,
            },
            genmap: None,
            rec_rate: 0.0,
        };
        edit(&mut cfg);
        let mut roh = Roh::new(cfg, 0);
        let mut out = Vec::new();
        for r in recs {
            let fields: Vec<&str> = r.split(' ').collect();
            let [chrom, pos, gts @ ..] = &fields[..] else {
                panic!("bad record {r:?}");
            };
            let mut line =
                format!("{chrom}\t{pos}\t.\tA\tC\t.\t.\t.\tGT\t{}", gts.join("\t")).into_bytes();
            let mut rec = Record::new();
            vcf::parse(&mut hdr, &mut line, &mut rec).expect("the line parses");
            roh.site(&hdr, &mut rec, None, &mut out)
                .expect("the site is processed");
        }
        roh.finish(&hdr, &mut out).expect("the output is written");
        String::from_utf8(out).expect("the output is text")
    }

    #[test]
    fn gts_only_skips_haploid_samples_at_mixed_ploidy_sites() {
        let recs: Vec<String> = (0..30)
            .map(|i| format!("1 {} 0 0/0", 1000 + 10_000 * i))
            .collect();
        let out = run_roh(&["male", "female"], &recs, |_| {});
        assert!(!out.contains("\tmale\t"), "{out}");
        assert_eq!(out.matches("ST\tfemale\t").count(), 30, "{out}");
    }

    /// A scratch directory of genetic maps, removed when dropped.
    struct Maps(std::path::PathBuf);

    impl Maps {
        fn new(name: &str) -> Self {
            let dir =
                std::env::temp_dir().join(format!("vtools-roh-{name}-{}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("the map directory is created");
            Self(dir)
        }

        /// Writes `NAME.txt`, a map of `cm_per_mb` up to 8 Mb, and returns the `-m` pattern of
        /// the directory.
        fn write(&self, name: &str, cm_per_mb: f64) -> String {
            let mut text = String::from("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n");
            for i in 0..=16 {
                text += &format!("{} 1 {}\n", i * 500_000, f64::from(i) * 0.5 * cm_per_mb);
            }
            std::fs::write(self.0.join(format!("{name}.txt")), text).expect("the map is written");
            self.pattern()
        }

        fn pattern(&self) -> String {
            format!("{}/{{CHROM}}.txt", self.0.display())
        }
    }

    impl Drop for Maps {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Sites on `chroms` every 20 kb, het outside an autozygous block at 1.5 to 3.5 Mb.
    fn block_sites(chroms: &[&str]) -> Vec<String> {
        let mut recs = Vec::new();
        for chrom in chroms {
            for i in 0..300 {
                let pos = 1000 + 20_000 * i;
                let gt = if (1_500_000..3_500_000).contains(&pos) || i % 3 == 0 {
                    "1/1"
                } else {
                    "0/1"
                };
                recs.push(format!("{chrom} {pos} {gt}"));
            }
        }
        recs
    }

    /// `roh -V` with the maps in `maps`.
    fn train(maps: &Maps, recs: &[String]) -> String {
        run_roh(&["S"], recs, |cfg| {
            cfg.training = Some(Training { threshold: 1e-10 });
            cfg.genmap = Some(maps.pattern());
        })
    }

    #[test]
    fn viterbi_training_uses_each_sequences_map() {
        let recs = block_sites(&["1", "2"]);
        let (a, b) = (Maps::new("train-a"), Maps::new("train-b"));
        a.write("1", 1e-4);
        b.write("1", 1e3);
        for m in [&a, &b] {
            m.write("2", 1.0);
        }
        // Only chromosome 1's map differs; bcftools used the map of 2 for both.
        assert_ne!(train(&a, &recs), train(&b, &recs));

        // Without a map for the last sequence, that sequence is skipped and the others keep
        // their maps; bcftools trained them with no map at all.
        let c = Maps::new("train-c");
        c.write("1", 1e3);
        assert_eq!(train(&c, &recs), train(&c, &block_sites(&["1"])));
    }

    #[test]
    fn gt_af_counts_alleles() {
        // 0/1, 1/1, ./., 0/0 as raw int8 BCF values (allele + 1) << 1.
        let gt = [2u8, 4, 4, 4, 0, 0, 2, 2];
        let Af::Value(af) = estimate_af_from_gt(&gt, None, 4) else {
            panic!("no AF");
        };
        assert_eq!(af, 3.0 / 6.0);
    }
}
