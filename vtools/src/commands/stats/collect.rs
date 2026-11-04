//! Collecting the statistics of bcftools stats, record by record (`do_vcf_stats` and its
//! helpers in `vcfstats.c`).
//!
//! A [`Collector`] holds one [`Stats`] per set (one for a single file; known and novel with
//! `-I`; A only, B only, and both for two files) and the scratch space of a record. Collectors
//! of worker threads are merged at the end; every statistic is a sum, so the merged result is
//! the serial one.

use pasta::variant::{
    FormatField, Header, HeaderRecordType, Record, Unpack, ValueType, VariantType,
    encoding::{INT8_VECTOR_END, INT16_VECTOR_END, INT32_MISSING, INT32_VECTOR_END},
    log,
    record::TypedValues,
    synced::Regions,
    utils::{self, AcSource, GenotypeType},
};

use super::{
    bins::{AfBins, DepthBins, Dist, IDist},
    indel_context::{IndelContext, REP_LEN},
};
use crate::zeroed;

/// The number of heterozygosity bins of the HWE statistics (bcftools' `naf_hwe`).
pub const NAF_HWE: usize = 100;
/// The longest indel length counted separately (bcftools' `m_indel`).
pub const M_INDEL: usize = 60;
/// The bins of the VAF distribution: 0 to 1 in steps of 0.05.
pub const VAF_BINS: usize = 21;

/// An error that ends the program with `code`, after printing `message` (when given) to
/// stderr. Messages of HTSlib functions are printed through [`log`] where they occur.
#[derive(Debug)]
pub struct Fatal {
    /// The exit status.
    pub code: u8,
    /// The message to print.
    pub message: Option<String>,
}

impl Fatal {
    /// An HTSlib `exit(1)` whose message was already logged.
    fn exit1() -> Self {
        Self {
            code: 1,
            message: None,
        }
    }

    /// A bcftools `error()`: the message, then exit status 255.
    pub fn error(message: String) -> Self {
        Self {
            code: 255,
            message: Some(message),
        }
    }
}

/// A `-u/--user-tstv` definition: `TAG[idx]:min:max:nbins`.
#[derive(Clone, Debug)]
pub struct UserStatDef {
    /// The INFO tag.
    pub tag: Vec<u8>,
    /// The value index.
    pub idx: usize,
    /// The lower bound of the first bin.
    pub min: f32,
    /// The upper bound of the last bin.
    pub max: f32,
    /// The number of bins.
    pub nbins: usize,
}

/// The settings shared by all collectors.
#[derive(Debug)]
pub struct Config {
    /// The number of allele-frequency bins (bcftools' `m_af`).
    pub m_af: usize,
    /// `--af-bins`.
    pub af_bins: Option<AfBins>,
    /// `--af-tag`.
    pub af_tag: Option<Vec<u8>>,
    /// `-1/--1st-allele-only`.
    pub first_allele_only: bool,
    /// `-I/--split-by-ID`.
    pub split_by_id: bool,
    /// The header columns of the samples of `-s/-S` in each input (bcftools'
    /// `reader->samples`); every input has the same number of them.
    pub samples: Vec<Vec<usize>>,
    /// Whether the header defines FORMAT/AD.
    pub has_fmt_ad: bool,
    /// Where allele counts come from (INFO only unless samples are given).
    pub ac_source: AcSource,
    /// `-d/--depth`.
    pub depth: DepthBins,
    /// `-u/--user-tstv`.
    pub user: Vec<UserStatDef>,
    /// Whether `-E/--exons` was given.
    pub exons: bool,
}

/// Per-sample counts (bcftools' `smpl_*` arrays).
#[derive(Clone, Copy, Debug, Default)]
pub struct SampleCounts {
    pub hets: u64,
    pub hom_rr: u64,
    pub hom_aa: u64,
    pub ts: u64,
    pub tv: u64,
    pub indels: u64,
    pub ndp: u64,
    pub dp: u64,
    pub sngl: u64,
    pub hap_ref: u64,
    pub hap_alt: u64,
    pub missing: u64,
    pub ins_hets: u64,
    pub del_hets: u64,
    pub ins_homs: u64,
    pub del_homs: u64,
    /// Indel genotypes not applicable to, in, and out of frame (with `-E`).
    pub frm_shifts: [u64; 3],
}

impl SampleCounts {
    fn merge(&mut self, o: &Self) {
        self.hets += o.hets;
        self.hom_rr += o.hom_rr;
        self.hom_aa += o.hom_aa;
        self.ts += o.ts;
        self.tv += o.tv;
        self.indels += o.indels;
        self.ndp += o.ndp;
        self.dp += o.dp;
        self.sngl += o.sngl;
        self.hap_ref += o.hap_ref;
        self.hap_alt += o.hap_alt;
        self.missing += o.missing;
        self.ins_hets += o.ins_hets;
        self.del_hets += o.del_hets;
        self.ins_homs += o.ins_homs;
        self.del_homs += o.del_homs;
        for (a, b) in self.frm_shifts.iter_mut().zip(o.frm_shifts) {
            *a += b;
        }
    }
}

/// The VAF distributions of one sample (bcftools' `vaf_t`).
#[derive(Clone, Copy, Debug)]
pub struct Vaf {
    pub snv: [u64; VAF_BINS],
    pub indel: [u64; VAF_BINS],
}

impl Default for Vaf {
    fn default() -> Self {
        Self {
            snv: [0; VAF_BINS],
            indel: [0; VAF_BINS],
        }
    }
}

/// The counts of one `-u/--user-tstv` tag.
#[derive(Clone, Debug)]
pub struct UserStats {
    /// Whether the tag is a Float (else an Integer).
    pub is_float: bool,
    pub vals_ts: Vec<u64>,
    pub vals_tv: Vec<u64>,
}

/// The statistics of one set (bcftools' `stats_t`).
#[derive(Clone, Debug)]
pub struct Stats {
    pub n_snps: u64,
    pub n_indels: u64,
    pub n_mnps: u64,
    pub n_others: u64,
    pub n_mals: u64,
    pub n_snp_mals: u64,
    pub n_records: u64,
    pub n_noalts: u64,
    /// By allele-frequency bin; bin 0 holds singletons.
    pub af_ts: Vec<u64>,
    pub af_tv: Vec<u64>,
    pub af_snps: Vec<u64>,
    /// Heterozygosity by allele frequency: `[iaf * NAF_HWE + ihet]` (with samples).
    pub af_hwe: Vec<u64>,
    /// Repeat-consistent and inconsistent deletions and insertions by repeat length.
    pub n_repeat: [[u64; 4]; REP_LEN],
    pub n_repeat_na: u64,
    /// Repeat-consistent, inconsistent, and not applicable indels by allele frequency.
    pub af_repeats: [Vec<u64>; 3],
    pub ts_alt1: u64,
    pub tv_alt1: u64,
    pub qual_ts: Dist,
    pub qual_tv: Dist,
    pub qual_indels: Dist,
    pub insertions: Vec<u64>,
    pub deletions: Vec<u64>,
    pub in_frame: i64,
    pub out_frame: i64,
    pub na_frame: i64,
    pub in_frame_alt1: i64,
    pub out_frame_alt1: i64,
    pub na_frame_alt1: i64,
    pub subst: [u64; 15],
    pub smpl: Vec<SampleCounts>,
    /// Per-sample VAF distributions (with samples, when the header defines FORMAT/AD).
    pub smpl_vaf: Option<Vec<Vaf>>,
    pub dp: IDist,
    pub dp_sites: IDist,
    pub usr: Vec<UserStats>,
    /// The sum and number of VAFs by allele length, `-M_INDEL..=M_INDEL` (with samples).
    pub dvaf: Vec<f64>,
    pub nvaf: Vec<u64>,
}

impl Config {
    /// The number of samples of `-s/-S` (bcftools' `files->n_smpl`).
    #[inline]
    pub fn n_smpl(&self) -> usize {
        self.samples.first().map_or(0, Vec::len)
    }
}

impl Stats {
    /// Creates empty statistics (`init_stats`); `user_is_float` gives the type of each `-u`
    /// tag in this set's header.
    ///
    /// # Errors
    ///
    /// Fails, with status 255, when the system cannot allocate the bins of `-d` or `-u`.
    pub fn try_new(cfg: &Config, user_is_float: &[bool]) -> Result<Self, Fatal> {
        let n_smpl = cfg.n_smpl();
        let with_samples = n_smpl > 0;
        let depth = || {
            IDist::try_new(cfg.depth).map_err(|e| {
                let DepthBins { min, max, step } = cfg.depth;
                Fatal::error(format!("Error: {e} for --depth {min},{max},{step}"))
            })
        };
        let user_bins = |u: &UserStatDef| {
            zeroed::counts(u.nbins).map_err(|e| {
                let tag = String::from_utf8_lossy(&u.tag);
                Fatal::error(format!(
                    "Error: {e} for --user-tstv {tag} ({} bins)",
                    u.nbins
                ))
            })
        };
        let mut usr = Vec::with_capacity(cfg.user.len());
        for (u, &is_float) in cfg.user.iter().zip(user_is_float) {
            usr.push(UserStats {
                is_float,
                vals_ts: user_bins(u)?,
                vals_tv: user_bins(u)?,
            });
        }
        Ok(Self {
            n_snps: 0,
            n_indels: 0,
            n_mnps: 0,
            n_others: 0,
            n_mals: 0,
            n_snp_mals: 0,
            n_records: 0,
            n_noalts: 0,
            af_ts: vec![0; cfg.m_af],
            af_tv: vec![0; cfg.m_af],
            af_snps: vec![0; cfg.m_af],
            af_hwe: if with_samples {
                vec![0; cfg.m_af * NAF_HWE]
            } else {
                Vec::new()
            },
            n_repeat: [[0; 4]; REP_LEN],
            n_repeat_na: 0,
            af_repeats: [vec![0; cfg.m_af], vec![0; cfg.m_af], vec![0; cfg.m_af]],
            ts_alt1: 0,
            tv_alt1: 0,
            qual_ts: Dist::new(5),
            qual_tv: Dist::new(5),
            qual_indels: Dist::new(5),
            insertions: vec![0; M_INDEL],
            deletions: vec![0; M_INDEL],
            in_frame: 0,
            out_frame: 0,
            na_frame: 0,
            in_frame_alt1: 0,
            out_frame_alt1: 0,
            na_frame_alt1: 0,
            subst: [0; 15],
            smpl: vec![SampleCounts::default(); n_smpl],
            smpl_vaf: (with_samples && cfg.has_fmt_ad).then(|| vec![Vaf::default(); n_smpl]),
            dp: depth()?,
            dp_sites: depth()?,
            usr,
            dvaf: if with_samples {
                vec![0.0; M_INDEL * 2 + 1]
            } else {
                Vec::new()
            },
            nvaf: if with_samples {
                vec![0; M_INDEL * 2 + 1]
            } else {
                Vec::new()
            },
        })
    }

    /// Adds the counts of `o`.
    pub fn merge(&mut self, o: &Self) {
        fn add(a: &mut [u64], b: &[u64]) {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
        self.n_snps += o.n_snps;
        self.n_indels += o.n_indels;
        self.n_mnps += o.n_mnps;
        self.n_others += o.n_others;
        self.n_mals += o.n_mals;
        self.n_snp_mals += o.n_snp_mals;
        self.n_records += o.n_records;
        self.n_noalts += o.n_noalts;
        add(&mut self.af_ts, &o.af_ts);
        add(&mut self.af_tv, &o.af_tv);
        add(&mut self.af_snps, &o.af_snps);
        add(&mut self.af_hwe, &o.af_hwe);
        for (a, b) in self.n_repeat.iter_mut().zip(&o.n_repeat) {
            add(a, b);
        }
        self.n_repeat_na += o.n_repeat_na;
        for (a, b) in self.af_repeats.iter_mut().zip(&o.af_repeats) {
            add(a, b);
        }
        self.ts_alt1 += o.ts_alt1;
        self.tv_alt1 += o.tv_alt1;
        self.qual_ts.merge(&o.qual_ts);
        self.qual_tv.merge(&o.qual_tv);
        self.qual_indels.merge(&o.qual_indels);
        add(&mut self.insertions, &o.insertions);
        add(&mut self.deletions, &o.deletions);
        self.in_frame += o.in_frame;
        self.out_frame += o.out_frame;
        self.na_frame += o.na_frame;
        self.in_frame_alt1 += o.in_frame_alt1;
        self.out_frame_alt1 += o.out_frame_alt1;
        self.na_frame_alt1 += o.na_frame_alt1;
        add(&mut self.subst, &o.subst);
        for (a, b) in self.smpl.iter_mut().zip(&o.smpl) {
            a.merge(b);
        }
        if let (Some(a), Some(b)) = (&mut self.smpl_vaf, &o.smpl_vaf) {
            for (x, y) in a.iter_mut().zip(b) {
                add(&mut x.snv, &y.snv);
                add(&mut x.indel, &y.indel);
            }
        }
        self.dp.merge(&o.dp);
        self.dp_sites.merge(&o.dp_sites);
        for (a, b) in self.usr.iter_mut().zip(&o.usr) {
            add(&mut a.vals_ts, &b.vals_ts);
            add(&mut a.vals_tv, &b.vals_tv);
        }
        for (a, b) in self.dvaf.iter_mut().zip(&o.dvaf) {
            *a += b;
        }
        add(&mut self.nvaf, &o.nvaf);
    }
}

/// What a record's alleles are, computed once per record.
#[derive(Clone, Copy, Debug, Default)]
struct AlleleInfo {
    /// The variant type without INS/DEL (`bcf_get_variant_type`).
    ty: VariantType,
    /// The variant length (`bcf_variant_length`).
    n: i32,
    /// The first base as 0-3 for A, C, G, T, else -1 (`bcf_acgt2int`).
    acgt: i32,
}

/// The header IDs of the tags read from every record, looked up again only when the header
/// grows (a VCF line can add tags, or define an existing ID in another section).
#[derive(Clone, Copy, Debug, Default)]
struct TagIds {
    /// The header's ID and line counts when these were looked up (zeros: not yet).
    version: (usize, usize),
    gt: Option<i32>,
    ad: Option<i32>,
    dp_fmt: Option<i32>,
    /// INFO/DP when it is an Integer tag (`bcf_get_info_int32` fails otherwise).
    dp_info: Option<i32>,
    /// Whether the header has both AN and AC, which `calc_ac` needs to read allele counts
    /// from INFO.
    an_ac: bool,
}

impl TagIds {
    #[inline]
    fn refresh(&mut self, header: &Header) {
        let version = (header.id_count(), header.record_count());
        if self.version == version {
            return;
        }
        let id = |tag: &[u8], kind: HeaderRecordType| {
            header
                .id(tag)
                .filter(|&id| header.id_exists(kind, id))
                .map(|id| id as i32)
        };
        self.version = version;
        self.gt = id(b"GT", HeaderRecordType::Format);
        self.ad = id(b"AD", HeaderRecordType::Format);
        self.dp_fmt = id(b"DP", HeaderRecordType::Format);
        self.dp_info = id(b"DP", HeaderRecordType::Info).filter(|&i| {
            header.id_type(HeaderRecordType::Info, i as u32) == Some(ValueType::Integer)
        });
        self.an_ac = header.id(b"AN").is_some() && header.id(b"AC").is_some();
    }
}

/// The value of an integer INFO field that has exactly one value
/// (`bcf_get_info_int32(...) == 1`).
fn info_int_single(rec: &mut Record, id: i32) -> Option<i32> {
    rec.unpack(Unpack::INFO);
    let f = rec.info_fields().iter().find(|f| f.key() == id)?;
    if f.is_removed() {
        return None;
    }
    let values = TypedValues::new(f.type_code(), rec.info_data(f));
    let (mut n, mut first) = (0, 0);
    for j in 0..f.len() {
        let v = values.int(j)?;
        if v == INT32_VECTOR_END {
            break;
        }
        if n == 0 {
            first = v;
        }
        n += 1;
    }
    (n == 1).then_some(first)
}

/// The per-record scratch space (bcftools' `tmp_iaf`, `tmp_frm`, and value buffers).
#[derive(Debug, Default)]
struct Scratch {
    /// The tag IDs of each input's header.
    ids: [TagIds; 2],
    /// The allele-frequency bin of each allele; 0 for singletons.
    tmp_iaf: Vec<i32>,
    /// The bin of the first ALT allele of the last record that had one: bcftools reads
    /// `tmp_iaf[1]` for a REF-only record without clearing it, so it keeps this value.
    iaf1: Option<i32>,
    /// The frameshift class of each allele (0 n/a, 1 in frame, 2 out of frame).
    tmp_frm: Vec<u8>,
    alleles: Vec<AlleleInfo>,
    farr: Vec<f32>,
    iarr: Vec<i32>,
}

/// One thread's statistics and scratch space.
#[derive(Debug)]
pub struct Collector {
    /// The statistics of each set.
    pub sets: Vec<Stats>,
    scratch: Scratch,
    indel_ctx: Option<IndelContext>,
    exons: Option<Regions>,
}

/// bcftools' `bcf_acgt2int`.
#[inline]
fn acgt2int(c: u8) -> i32 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => -1,
    }
}

/// C's conversion of a float to `int` on x86-64: out-of-range values become `INT_MIN`.
#[inline]
fn c_f32_to_i32(x: f32) -> i32 {
    if x.is_nan() || !(-2_147_483_648.0..2_147_483_648.0).contains(&x) {
        i32::MIN
    } else {
        x as i32
    }
}

/// The QUAL histogram value: `1 + QUAL*10`, or 0 for missing and negative values.
#[inline]
fn qual_bin(qual: f32) -> u32 {
    if qual.is_nan() || qual < 0.0 {
        0
    } else {
        1i32.wrapping_add(c_f32_to_i32(qual * 10.0)) as u32
    }
}

/// The byte width of an integer BCF type, `None` for other types.
#[inline]
fn int_width(ty: u8) -> Option<usize> {
    match ty {
        1 => Some(1),
        2 => Some(2),
        3 => Some(4),
        _ => None,
    }
}

/// Reads integer `i` of `data` (of integer type `ty`), with sentinels widened to 32 bits;
/// past the end of `data` reads as a vector end.
#[inline]
fn read_int(ty: u8, data: &[u8], i: usize) -> i32 {
    let w = int_width(ty).unwrap_or(4);
    if (i + 1) * w > data.len() {
        return INT32_VECTOR_END;
    }
    TypedValues::new(ty, data)
        .int(i)
        .unwrap_or(INT32_VECTOR_END)
}

fn todo(func: &str, ty: u8) -> Fatal {
    log::error(func, format_args!("todo: {ty}"));
    Fatal::exit1()
}

impl Collector {
    /// Creates a collector with `nstats` sets; `user_types[set]` gives the `-u` tag types of
    /// each set's header.
    ///
    /// # Errors
    ///
    /// Fails, with status 255, when the system cannot allocate the bins of `-d` or `-u`.
    pub fn try_new(
        cfg: &Config,
        user_types: &[Vec<bool>],
        indel_ctx: Option<IndelContext>,
        exons: Option<Regions>,
    ) -> Result<Self, Fatal> {
        Ok(Self {
            sets: user_types
                .iter()
                .map(|u| Stats::try_new(cfg, u))
                .collect::<Result<_, _>>()?,
            scratch: Scratch::default(),
            indel_ctx,
            exons,
        })
    }

    /// The frequency bin of the first ALT allele of the last record that had one, which
    /// bcftools reads as `tmp_iaf[1]` when comparing genotypes; `None` before any such
    /// record (bcftools' `ntmp_iaf <= 1`).
    #[inline]
    pub fn first_alt_bin(&self) -> Option<usize> {
        self.scratch.iaf1.map(|b| b as usize)
    }

    /// Adds the counts of another collector.
    pub fn merge(&mut self, other: &Self) {
        for (a, b) in self.sets.iter_mut().zip(&other.sets) {
            a.merge(b);
        }
    }

    /// Collects the statistics of a record of input `reader` (0 or 1) for set `set` (the body
    /// of `do_vcf_stats`).
    ///
    /// # Errors
    ///
    /// Returns the errors bcftools exits on.
    pub fn add(
        &mut self,
        cfg: &Config,
        header: &Header,
        rec: &mut Record,
        set: usize,
        reader: usize,
    ) -> Result<(), Fatal> {
        let line_type = rec.variant_types();
        let sc = &mut self.scratch;
        sc.alleles.clear();
        for i in 0..rec.allele_count() {
            let v = rec.variant(i).unwrap_or_default();
            sc.alleles.push(AlleleInfo {
                ty: v.ty & VariantType::ORIGINAL,
                n: v.n,
                acgt: acgt2int(rec.allele(i).first().copied().unwrap_or(0)),
            });
        }
        sc.ids[reader].refresh(header);
        init_iaf(cfg, sc, reader, header, rec)?;

        let set = if cfg.split_by_id && rec.id() == b"." {
            1
        } else {
            set
        };
        let stats = &mut self.sets[set];
        stats.n_records += 1;

        if line_type.is_empty() {
            stats.n_noalts += 1;
        }
        if line_type.contains(VariantType::SNP) {
            do_snp_stats(cfg, stats, sc, header, rec);
        }
        if line_type.contains(VariantType::INDEL) {
            do_indel_stats(
                cfg,
                stats,
                sc,
                header,
                rec,
                self.indel_ctx.as_mut(),
                self.exons.as_mut(),
            )?;
        }
        if line_type.contains(VariantType::MNP) {
            stats.n_mnps += 1;
        }
        if line_type.contains(VariantType::OTHER) {
            stats.n_others += 1;
        }
        if rec.allele_count() > 2 {
            stats.n_mals += 1;
            if line_type == VariantType::SNP {
                stats.n_snp_mals += 1;
            }
        }
        if cfg.n_smpl() > 0 {
            do_sample_stats(cfg, stats, sc, reader, rec)?;
        }
        if let Some(id) = sc.ids[reader].dp_info
            && let Some(dp) = info_int_single(rec, id)
        {
            stats.dp_sites.add(dp);
        }
        Ok(())
    }
}

/// Bins each allele by frequency (`init_iaf`), and remembers the first ALT allele's bin.
fn init_iaf(
    cfg: &Config,
    sc: &mut Scratch,
    reader: usize,
    header: &Header,
    rec: &mut Record,
) -> Result<(), Fatal> {
    bin_alleles(cfg, sc, reader, header, rec)?;
    if rec.allele_count() > 1 {
        sc.iaf1 = Some(sc.tmp_iaf[1]);
    }
    Ok(())
}

/// Bins each allele by frequency: `tmp_iaf[i]` is 0 for singletons and unknown frequencies,
/// else 1 + the frequency bin.
fn bin_alleles(
    cfg: &Config,
    sc: &mut Scratch,
    reader: usize,
    header: &Header,
    rec: &mut Record,
) -> Result<(), Fatal> {
    let n = rec.allele_count();
    sc.tmp_iaf.clear();
    sc.tmp_iaf.resize(n.max(1), 0);
    let to_bin = |af: f32| -> i32 {
        let af = af.clamp(0.0, 1.0);
        match &cfg.af_bins {
            Some(b) => b.index(af) as i32,
            None => c_f32_to_i32(af * (cfg.m_af as f32 - 2.0)),
        }
    };
    if let Some(tag) = &cfg.af_tag {
        let ok =
            matches!(rec.get_info_f32(header, tag, &mut sc.farr), Ok(k) if k > 0 && k + 1 == n);
        if !ok {
            sc.tmp_iaf.fill(0);
            return Ok(());
        }
        sc.tmp_iaf[0] = 0;
        for i in 1..n {
            let af = sc.farr[i - 1];
            // bcftools crashes on a missing value; it goes to the unknown bin here.
            sc.tmp_iaf[i] = if af.is_nan() { 0 } else { to_bin(af) + 1 };
        }
        return Ok(());
    }
    // Without AN and AC in the header, calc_ac goes straight to the genotypes (if allowed);
    // deciding that here saves its tag lookups.
    let source = match (sc.ids[reader].an_ac, cfg.ac_source) {
        (true, s) => Some(s),
        (false, AcSource::Info) => None,
        (false, _) => Some(AcSource::Genotypes),
    };
    let found = match source {
        Some(s) => utils::calc_ac(header, rec, &mut sc.tmp_iaf, s).map_err(|_| Fatal::exit1())?,
        None => false,
    };
    if !found {
        sc.tmp_iaf.fill(0);
        return Ok(());
    }
    let an: i32 = sc.tmp_iaf[..n].iter().fold(0i32, |a, &b| a.wrapping_add(b));
    if n > 0 {
        sc.tmp_iaf[0] = 0;
    }
    for i in 1..n {
        let ac = sc.tmp_iaf[i];
        sc.tmp_iaf[i] = if ac == 1 {
            0
        } else if an == 0 {
            1
        } else {
            to_bin(ac as f32 / an as f32) + 1
        };
    }
    Ok(())
}

fn do_snp_stats(
    cfg: &Config,
    stats: &mut Stats,
    sc: &mut Scratch,
    header: &Header,
    rec: &mut Record,
) {
    stats.n_snps += 1;
    let r = sc.alleles[0].acgt;
    if r < 0 {
        return;
    }
    let iqual = qual_bin(rec.qual_raw());
    for i in 1..rec.allele_count() {
        if cfg.first_allele_only && i > 1 {
            break;
        }
        if !sc.alleles[i].ty.contains(VariantType::SNP) {
            continue;
        }
        let a = sc.alleles[i].acgt;
        if a < 0 || r == a {
            continue;
        }
        stats.subst[(r << 2 | a) as usize] += 1;
        let iaf = sc.tmp_iaf[i] as usize;
        stats.af_snps[iaf] += 1;
        let is_ts = (r - a).abs() == 2;
        if i == 1 {
            if is_ts {
                stats.ts_alt1 += 1;
                stats.qual_ts.insert(iqual);
            } else {
                stats.tv_alt1 += 1;
                stats.qual_tv.insert(iqual);
            }
            do_user_stats(cfg, stats, sc, header, rec, is_ts);
        }
        if is_ts {
            stats.af_ts[iaf] += 1;
        } else {
            stats.af_tv[iaf] += 1;
        }
    }
}

fn do_user_stats(
    cfg: &Config,
    stats: &mut Stats,
    sc: &mut Scratch,
    header: &Header,
    rec: &mut Record,
    is_ts: bool,
) {
    for (def, usr) in cfg.user.iter().zip(&mut stats.usr) {
        let val = if usr.is_float {
            match rec.get_info_f32(header, &def.tag, &mut sc.farr) {
                Ok(n) if n > 0 && def.idx < n => sc.farr[def.idx],
                _ => continue,
            }
        } else {
            match rec.get_info_i32(header, &def.tag, &mut sc.iarr) {
                Ok(n) if n > 0 && def.idx < n => sc.iarr[def.idx] as f32,
                _ => continue,
            }
        };
        // bcftools crashes on a missing Float value; it is skipped here.
        if val.is_nan() {
            continue;
        }
        let idx = if val <= def.min {
            0
        } else if val >= def.max {
            def.nbins - 1
        } else {
            c_f32_to_i32((val - def.min) / (def.max - def.min) * (def.nbins as f32 - 1.0)) as usize
        };
        if is_ts {
            usr.vals_ts[idx] += 1;
        } else {
            usr.vals_tv[idx] += 1;
        }
    }
}

fn do_indel_stats(
    cfg: &Config,
    stats: &mut Stats,
    sc: &mut Scratch,
    header: &Header,
    rec: &mut Record,
    mut indel_ctx: Option<&mut IndelContext>,
    exons: Option<&mut Regions>,
) -> Result<(), Fatal> {
    stats.n_indels += 1;
    stats.qual_indels.insert(qual_bin(rec.qual_raw()));

    let n = rec.allele_count();
    let pos = rec.pos();
    let chr = header.contig_name(rec.rid()).unwrap_or(b"");
    let mut exon = None;
    if let Some(exons) = exons {
        // bcftools counts a site as outside the exons when an exon line is malformed (the
        // error was reported). It does the same for a line that cannot be read, which loses
        // the exons after it: that is an error here.
        let overlap = match exons.overlaps(chr, pos, pos) {
            Ok(overlap) => overlap,
            Err(e) if e.is_read_failure() => {
                return Err(Fatal::error(format!(
                    "Error: failed to read the exons: {e}"
                )));
            }
            Err(_) => false,
        };
        if overlap {
            exon = exons.current().map(|(_, s, e)| (s, e));
        }
        sc.tmp_frm.clear();
        sc.tmp_frm.resize(n, 0);
    }
    for i in 1..n {
        if cfg.first_allele_only && i > 1 {
            break;
        }
        if !sc.alleles[i].ty.contains(VariantType::INDEL) {
            continue;
        }
        let len = sc.alleles[i].n;
        let iaf = sc.tmp_iaf[i] as usize;

        if let Some(ctx) = indel_ctx.as_deref_mut() {
            let c = ctx
                .classify(chr, pos + 1, rec.allele(0), rec.allele(i))
                .map_err(|e| Fatal::error(e.0))?;
            if c.nlen <= 1 || c.nrep <= 1 {
                stats.n_repeat_na += 1;
                stats.af_repeats[2][iaf] += 1;
            } else if !(c.ndel.unsigned_abs() as usize).is_multiple_of(c.nlen) {
                stats.n_repeat[c.nlen - 1][if c.ndel < 0 { 1 } else { 3 }] += 1;
                stats.af_repeats[1][iaf] += 1;
            } else {
                stats.n_repeat[c.nlen - 1][if c.ndel < 0 { 0 } else { 2 }] += 1;
                stats.af_repeats[0][iaf] += 1;
            }
        } else {
            stats.af_repeats[2][iaf] += 1;
        }

        // Frameshifts.
        let mut tlen: i64 = 0;
        if let Some((estart, eend)) = exon {
            let alen = i64::from(len.unsigned_abs());
            if len > 0 {
                if estart <= pos && eend > pos {
                    tlen = alen;
                }
            } else if estart <= pos + alen {
                tlen = alen;
                if pos < estart {
                    tlen -= estart - pos + 1;
                }
                if eend < pos + alen {
                    tlen -= pos + alen - eend;
                }
            }
        }
        if tlen != 0 {
            if tlen % 3 != 0 {
                stats.out_frame += 1;
                sc.tmp_frm[i] = 2;
            } else {
                stats.in_frame += 1;
                sc.tmp_frm[i] = 1;
            }
            if i == 1 {
                if tlen % 3 != 0 {
                    stats.out_frame_alt1 += 1;
                } else {
                    stats.in_frame_alt1 += 1;
                }
            }
        } else {
            if i == 1 {
                stats.na_frame_alt1 += 1;
            }
            stats.na_frame += 1;
        }

        // The indel length distribution.
        let counts = if len < 0 {
            &mut stats.deletions
        } else {
            &mut stats.insertions
        };
        let l = (len.unsigned_abs() as usize)
            .saturating_sub(1)
            .min(M_INDEL - 1);
        counts[l] += 1;
    }
    Ok(())
}

/// The per-record genotype tallies for the HWE and singleton statistics.
#[derive(Default)]
struct GtTally {
    nref_tot: u32,
    nhet_tot: u32,
    nalt_tot: u32,
    n_nref: u32,
    i_nref: usize,
}

/// Returns allele `i`'s type, or HTSlib's error for an allele out of range
/// (`bcf_get_variant_type`).
#[inline]
fn allele_type(alleles: &[AlleleInfo], i: i32) -> Result<VariantType, Fatal> {
    match alleles.get(i as usize) {
        Some(a) if i >= 0 => Ok(a.ty),
        _ => {
            log::error(
                "bcf_get_variant_type",
                format_args!("Requested allele outside valid range"),
            );
            Err(Fatal::exit1())
        }
    }
}

/// Counts one sample's genotype (`sample_gt_stats`).
#[expect(
    clippy::too_many_arguments,
    reason = "a line-for-line port of sample_gt_stats"
)]
#[inline]
fn sample_gt_stats(
    sc: &Scratch,
    smpl: &mut SampleCounts,
    tally: &mut GtTally,
    is: usize,
    gt: GenotypeType,
    ial: i32,
    jal: i32,
    frm_shifts: bool,
) -> Result<(), Fatal> {
    if gt == GenotypeType::Unknown {
        smpl.missing += 1;
        return Ok(());
    }
    let mut var_type = VariantType::empty();
    if ial > 0 {
        var_type |= allele_type(&sc.alleles, ial)?;
    }
    if jal > 0 {
        var_type |= allele_type(&sc.alleles, jal)?;
    }
    if matches!(gt, GenotypeType::HaploidRef | GenotypeType::HaploidAlt) {
        if var_type.contains(VariantType::INDEL) && frm_shifts {
            smpl.frm_shifts[usize::from(sc.tmp_frm[ial as usize])] += 1;
        }
        if gt == GenotypeType::HaploidRef {
            smpl.hap_ref += 1;
        } else {
            smpl.hap_alt += 1;
        }
        return Ok(());
    }
    if gt != GenotypeType::HomRef {
        tally.n_nref += 1;
        tally.i_nref = is;
    }
    match gt {
        GenotypeType::HomRef => tally.nref_tot += 1,
        GenotypeType::HetRefAlt => tally.nhet_tot += 1,
        GenotypeType::HetAltAlt | GenotypeType::HomAlt => tally.nalt_tot += 1,
        _ => {}
    }

    if var_type.contains(VariantType::SNP) || var_type.is_empty() {
        match gt {
            GenotypeType::HetRefAlt | GenotypeType::HetAltAlt => smpl.hets += 1,
            GenotypeType::HomRef => smpl.hom_rr += 1,
            GenotypeType::HomAlt => smpl.hom_aa += 1,
            _ => {}
        }
        if gt != GenotypeType::HomRef {
            for (al, check) in [(ial, true), (jal, ial != jal)] {
                if !check || !sc.alleles[al as usize].ty.contains(VariantType::SNP) {
                    continue;
                }
                let r = sc.alleles[0].acgt;
                let a = sc.alleles[al as usize].acgt;
                if a < 0 {
                    return Ok(());
                }
                if (r - a).abs() == 2 {
                    smpl.ts += 1;
                } else {
                    smpl.tv += 1;
                }
            }
        }
    }
    if var_type.contains(VariantType::INDEL) {
        if gt != GenotypeType::HomRef {
            smpl.indels += 1;
            if matches!(gt, GenotypeType::HetRefAlt | GenotypeType::HetAltAlt) {
                let (mut is_ins, mut is_del) = (false, false);
                for al in [ial, jal] {
                    let a = sc.alleles[al as usize];
                    if a.ty.contains(VariantType::INDEL) {
                        if a.n < 0 {
                            is_del = true;
                        } else {
                            is_ins = true;
                        }
                    }
                }
                // Alt-het genotypes with both an insertion and a deletion count twice.
                if is_del {
                    smpl.del_hets += 1;
                }
                if is_ins {
                    smpl.ins_hets += 1;
                }
            } else if gt == GenotypeType::HomAlt {
                if sc.alleles[ial as usize].n < 0 {
                    smpl.del_homs += 1;
                } else {
                    smpl.ins_homs += 1;
                }
            }
        }
        if frm_shifts {
            smpl.frm_shifts[usize::from(sc.tmp_frm[ial as usize])] += 1;
            smpl.frm_shifts[usize::from(sc.tmp_frm[jal as usize])] += 1;
        }
    }
    Ok(())
}

/// Classifies a genotype from its raw BCF values (HTSlib's `bcf_gt_type`), returning the
/// type and the smallest and largest alternate allele (0-based; 0 when none).
#[inline]
fn gt_type(ty: u8, p: &[u8]) -> (GenotypeType, i32, i32) {
    match ty {
        1 => gt_type_values(p.iter().map(|&b| {
            let v = b as i8;
            if v == INT8_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(v)
            }
        })),
        2 => gt_type_values(p.as_chunks::<2>().0.iter().map(|&c| {
            let v = i16::from_le_bytes(c);
            if v == INT16_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(v)
            }
        })),
        _ => gt_type_values(p.as_chunks::<4>().0.iter().map(|&c| i32::from_le_bytes(c))),
    }
}

#[inline(always)]
fn gt_type_values(values: impl Iterator<Item = i32>) -> (GenotypeType, i32, i32) {
    let (mut nals, mut has_ref, mut has_alt) = (0u32, false, 0u8);
    let (mut ial, mut jal) = (0i32, 0i32);
    for val in values {
        if val == INT32_VECTOR_END {
            break;
        }
        let tmp = val >> 1;
        if tmp == 0 {
            return (GenotypeType::Unknown, 0, 0);
        }
        if tmp > 1 {
            if ial == 0 {
                ial = tmp;
                has_alt = 1;
            } else if tmp != ial {
                if tmp < ial {
                    jal = ial;
                    ial = tmp;
                } else {
                    jal = tmp;
                }
                has_alt = 2;
            }
        } else {
            has_ref = true;
        }
        nals += 1;
    }
    let ial = if ial > 0 { ial - 1 } else { ial };
    let jal = if jal > 0 { jal - 1 } else { jal };
    let ty = match (nals, has_ref, has_alt) {
        (0, _, _) => GenotypeType::Unknown,
        (1, true, _) => GenotypeType::HaploidRef,
        (1, false, _) => GenotypeType::HaploidAlt,
        (_, false, 1) => GenotypeType::HomAlt,
        (_, false, _) => GenotypeType::HetAltAlt,
        (_, true, 0) => GenotypeType::HomRef,
        _ => GenotypeType::HetRefAlt,
    };
    (ty, ial, jal)
}

/// A sample's depth: FORMAT/DP, else the sum of FORMAT/AD, else -1 (`calc_sample_depth`).
#[inline]
fn calc_sample_depth(
    rec: &Record,
    ismpl: usize,
    ad: Option<&FormatField>,
    dp: Option<&FormatField>,
) -> Result<i32, Fatal> {
    if let Some(f) = dp {
        let ty = f.type_code();
        if int_width(ty).is_none() {
            return Err(todo("calc_sample_depth", ty));
        }
        let v = read_int(ty, rec.format_sample_data(f, ismpl), 0);
        return Ok(if v == INT32_MISSING || v == INT32_VECTOR_END {
            -1
        } else {
            v
        });
    }
    if let Some(f) = ad {
        let ty = f.type_code();
        if int_width(ty).is_none() {
            return Err(todo("calc_sample_depth", ty));
        }
        let data = rec.format_sample_data(f, ismpl);
        let (mut dp, mut has_value) = (0i32, false);
        for iv in 0..f.n() {
            let v = read_int(ty, data, iv);
            if v == INT32_VECTOR_END {
                break;
            }
            if v == INT32_MISSING {
                continue;
            }
            has_value = true;
            dp = dp.wrapping_add(v);
        }
        return Ok(if has_value { dp } else { -1 });
    }
    Ok(-1)
}

/// The largest non-reference AD of a sample and its allele (`get_ad`).
fn get_ad(rec: &Record, f: &FormatField, ismpl: usize) -> Result<(i32, i32), Fatal> {
    let ty = f.type_code();
    if int_width(ty).is_none() {
        return Err(todo("get_ad", ty));
    }
    let data = rec.format_sample_data(f, ismpl);
    let (mut ad, mut ial) = (0i32, 0i32);
    let mut iv = 1;
    while iv < f.n() && iv < rec.allele_count() {
        let v = read_int(ty, data, iv);
        if v == INT32_VECTOR_END {
            break;
        }
        if v != INT32_MISSING && ad < v {
            ad = v;
            ial = iv as i32;
        }
        iv += 1;
    }
    Ok((ad, ial))
}

/// The AD of allele `ial` of a sample (`get_iad`); bcftools reads past the sample's values
/// when `ial` is beyond them, which reads the next sample's.
fn get_iad(rec: &Record, f: &FormatField, ismpl: usize, ial: i32) -> Result<i32, Fatal> {
    let ty = f.type_code();
    let Some(width) = int_width(ty) else {
        return Err(todo("get_iad", ty));
    };
    let off = f.size() * ismpl + ial as usize * width;
    let v = match rec.format_data(f).get(off..off + width) {
        Some(bytes) => read_int(ty, bytes, 0),
        None => INT32_MISSING,
    };
    Ok(if v == INT32_VECTOR_END || v == INT32_MISSING {
        0
    } else {
        v
    })
}

fn update_dvaf(stats: &mut Stats, sc: &Scratch, ial: i32, vaf: f32) {
    let len = sc.alleles[ial as usize]
        .n
        .clamp(-(M_INDEL as i32), M_INDEL as i32);
    let bin = (M_INDEL as i32 + len) as usize;
    stats.nvaf[bin] += 1;
    stats.dvaf[bin] += f64::from(vaf);
}

/// Rounds a small non-negative value to the nearest integer, ties to even (C's `nearbyintf`),
/// without a libm call.
#[inline]
fn round_half_even(x: f32) -> usize {
    let t = x as usize;
    let frac = x - t as f32;
    if frac > 0.5 || (frac == 0.5 && t % 2 == 1) {
        t + 1
    } else {
        t
    }
}

fn update_vaf(vaf_counts: Option<&mut Vaf>, sc: &Scratch, ial: i32, vaf: f32) {
    let Some(v) = vaf_counts else { return };
    if (0.0..=1.0).contains(&vaf) {
        let idx = round_half_even((f64::from(vaf) / 0.05) as f32);
        if sc.alleles[ial as usize].ty == VariantType::SNP {
            v.snv[idx] += 1;
        } else {
            v.indel[idx] += 1;
        }
    }
}

fn do_sample_stats(
    cfg: &Config,
    stats: &mut Stats,
    sc: &mut Scratch,
    reader: usize,
    rec: &mut Record,
) -> Result<(), Fatal> {
    rec.unpack(Unpack::FMT);
    let rec = &*rec;
    let find =
        |id: Option<i32>| id.and_then(|id| rec.format_fields().iter().find(|f| f.id() == id));
    let ids = sc.ids[reader];
    let gt_fmt = find(ids.gt);
    let ad_fmt = find(ids.ad);
    let dp_fmt = find(ids.dp_fmt);
    let frm_shifts = cfg.exons;
    let mut tally = GtTally::default();
    let gt_data = gt_fmt.map_or(&[][..], |f| rec.format_data(f));

    if let Some(f) = gt_fmt
        && !matches!(f.type_code(), 1..=3)
    {
        log::error(
            "bcf_gt_type",
            format_args!("Unexpected type {}", f.type_code()),
        );
        return Err(Fatal::exit1());
    }
    let has_depth = ad_fmt.is_some() || dp_fmt.is_some();
    // Diploid 8-bit genotypes, the common case, have a shortcut for homozygous reference.
    let diploid8 = gt_fmt.is_some_and(|f| f.type_code() == 1 && f.size() == 2);
    for (is, &ismpl) in cfg.samples[reader].iter().enumerate() {
        let dp = if has_depth {
            calc_sample_depth(rec, ismpl, ad_fmt, dp_fmt)?
        } else {
            -1
        };
        if dp > 0 {
            stats.dp.add(dp);
            stats.smpl[is].ndp += 1;
            stats.smpl[is].dp += dp as u64;
        }

        let (mut gt, mut ial, mut jal) = (GenotypeType::Unknown, 0, 0);
        if diploid8 && gt_data[ismpl * 2] & 0xfe == 2 && gt_data[ismpl * 2 + 1] & 0xfe == 2 {
            // 0/0: what sample_gt_stats counts for it.
            gt = GenotypeType::HomRef;
            stats.smpl[is].hom_rr += 1;
            tally.nref_tot += 1;
        } else if let Some(f) = gt_fmt {
            let size = f.size();
            let p = &gt_data[ismpl * size..(ismpl + 1) * size];
            (gt, ial, jal) = gt_type(f.type_code(), p);
            sample_gt_stats(
                sc,
                &mut stats.smpl[is],
                &mut tally,
                is,
                gt,
                ial,
                jal,
                frm_shifts,
            )?;
        }

        if dp > 0
            && let Some(f) = ad_fmt
        {
            let (mut iad, mut jad) = (0f32, 0f32);
            if gt == GenotypeType::Unknown {
                let (ad, al) = get_ad(rec, f, ismpl)?;
                iad = ad as f32;
                ial = al;
            } else {
                if ial != 0 {
                    iad = get_iad(rec, f, ismpl, ial)? as f32;
                }
                if jal != 0 {
                    jad = get_iad(rec, f, ismpl, jal)? as f32;
                }
            }
            if iad != 0.0 {
                let vaf = iad / dp as f32;
                update_dvaf(stats, sc, ial, vaf);
                update_vaf(stats.smpl_vaf.as_mut().map(|v| &mut v[is]), sc, ial, vaf);
            }
            if jad != 0.0 && iad != jad {
                let vaf = jad / dp as f32;
                update_dvaf(stats, sc, jal, vaf);
                update_vaf(stats.smpl_vaf.as_mut().map(|v| &mut v[is]), sc, jal, vaf);
            }
        }
    }
    if tally.n_nref == 1 {
        stats.smpl[tally.i_nref].sngl += 1;
    }

    if gt_fmt.is_some()
        && rec.allele_count() > 1
        && (tally.nref_tot != 0 || tally.nhet_tot != 0 || tally.nalt_tot != 0)
    {
        let het_frac =
            tally.nhet_tot as f32 / (tally.nref_tot + tally.nhet_tot + tally.nalt_tot) as f32;
        let ihet = c_f32_to_i32(het_frac * (NAF_HWE as f32 - 1.0));
        let idx = ihet as usize + sc.tmp_iaf[1] as usize * NAF_HWE;
        stats.af_hwe[idx] += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rounding_ties_to_even() {
        assert_eq!(round_half_even(0.0), 0);
        assert_eq!(round_half_even(0.5), 0);
        assert_eq!(round_half_even(1.5), 2);
        assert_eq!(round_half_even(2.5), 2);
        assert_eq!(round_half_even(2.6), 3);
        assert_eq!(round_half_even(19.999), 20);
    }

    #[test]
    fn qual_bins() {
        assert_eq!(qual_bin(f32::NAN), 0);
        assert_eq!(qual_bin(-1.0), 0);
        assert_eq!(qual_bin(0.0), 1);
        assert_eq!(qual_bin(12.34), 124);
        assert_eq!(qual_bin(1e30), 2_147_483_649);
    }
}
