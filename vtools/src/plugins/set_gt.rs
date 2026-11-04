//! `+setGT`: set chosen genotypes (missing, partially missing, all, those passing `-i/-e`, a
//! binomial test, or a random fraction) to missing, reference, major or minor allele, the
//! allele with the most reads, a custom genotype, or change their phase (bcftools'
//! `plugins/setGT.c`).

use anyhow::bail;
use clap::{ArgMatches, CommandFactory, FromArgMatches, Parser};
use pasta::variant::{
    Header, HeaderRecordType, Record, Unpack,
    encoding::INT32_VECTOR_END,
    utils::{self, AcSource},
};

use super::{Exit, InitHeaders, Loaded, Output, Plugin};
use crate::{
    filter::{Filter, FilterLogic, kfunc::binom_two_sided},
    shortopt,
};

/// The description `plugin -lv` prints.
pub const ABOUT: &str =
    "Set genotypes: partially missing to missing, missing to ref/major allele, etc.\n";

/// The plugin's usage text.
pub fn usage() -> String {
    "About: Sets genotypes. The target genotypes can be specified as:\n\
\x20          ./.     .. completely missing (\".\" or \"./.\", depending on ploidy)\n\
\x20          ./x     .. partially missing (e.g., \"./0\" or \".|1\" but not \"./.\")\n\
\x20          .       .. partially or completely missing\n\
\x20          a       .. all genotypes\n\
\x20          b       .. heterozygous genotypes failing two-tailed binomial test (example below)\n\
\x20          q       .. select genotypes using -i/-e options\n\
\x20          r:FLOAT .. select randomly a proportion of FLOAT genotypes (can be combined with other modes)\n\
\x20      and the new genotype can be one of:\n\
\x20          .       .. missing (\".\" or \"./.\", keeps ploidy)\n\
\x20          0       .. reference allele (e.g. 0/0 or 0, keeps ploidy)\n\
\x20          c:GT    .. custom genotype (e.g. 0/0, 0, 0/1, m/M, 0/X, ./., .;  overrides ploidy)\n\
\x20          i       .. invert the genotype phase (0|1 becomes 1|0)\n\
\x20          m       .. minor (the second most common) allele as determined from INFO/AC or FMT/GT (e.g. 1/1 or 1, keeps ploidy)\n\
\x20          M       .. major allele as determined from INFO/AC or FMT/GT (e.g. 1/1 or 1, keeps ploidy)\n\
\x20          X       .. allele with bigger read depth as determined from FMT/AD\n\
\x20          p       .. phase genotype (0/1 becomes 0|1)\n\
\x20          u       .. unphase genotype and sort by allele (1|0 becomes 0/1)\n\
Usage: vtools +setGT [General Options] -- [Plugin Options]\n\
Options:\n\
\x20  run \"vtools plugin\" for a list of common options\n\
\n\
Plugin options:\n\
\x20  -e, --exclude EXPR        Exclude a genotype if true (requires -t q)\n\
\x20  -i, --include EXPR        include a genotype if true (requires -t q)\n\
\x20  -n, --new-gt TYPE         Genotypes to set, see above\n\
\x20  -s, --seed INT            Random seed to use with -t r [0]\n\
\x20  -t, --target-gt TYPE      Genotypes to change, see above\n\
\n\
Example:\n\
\x20  # set missing genotypes (\"./.\") to phased ref genotypes (\"0|0\")\n\
\x20  vtools +setGT in.vcf -- -t . -n 0p\n\
\n\
\x20  # set missing genotypes with DP>0 and GQ>20 to ref genotypes (\"0/0\")\n\
\x20  vtools +setGT in.vcf -- -t q -n 0 -i 'GT=\".\" && FMT/DP>0 && GQ>20'\n\
\n\
\x20  # set partially missing genotypes to completely missing\n\
\x20  vtools +setGT in.vcf -- -t ./x -n .\n\
\n\
\x20  # set heterozygous genotypes to 0/0 if binom.test(nAlt,nRef+nAlt,0.5)<1e-3\n\
\x20  vtools +setGT in.vcf -- -t \"b:AD<1e-3\" -n 0\n\
\n\
\x20  # force unphased heterozygous genotype if binom.test(nAlt,nRef+nAlt,0.5)>0.1\n\
\x20  vtools +setGT in.vcf -- -t ./x -n c:'m/M'\n\
\n"
    .to_owned()
}

/// The plugin options (bcftools' getopt string `?hn:t:i:e:s:`).
#[derive(Debug, Parser)]
#[command(
    no_binary_name = true,
    disable_help_flag = true,
    disable_version_flag = true
)]
struct Args {
    /// Arguments left over (ignored, as getopt leaves them).
    #[arg(hide = true)]
    _rest: Vec<String>,
    /// Exclude a genotype if true (requires -t q).
    #[arg(short = 'e', long = "exclude", allow_hyphen_values = true)]
    exclude: Vec<String>,
    /// Include a genotype if true (requires -t q).
    #[arg(short = 'i', long = "include", allow_hyphen_values = true)]
    include: Vec<String>,
    /// Genotypes to set.
    #[arg(short = 'n', long = "new-gt", allow_hyphen_values = true)]
    new_gt: Vec<String>,
    /// Random seed to use with -t r [0].
    #[arg(short = 's', long = "seed", allow_hyphen_values = true)]
    seed: Vec<String>,
    /// Genotypes to change.
    #[arg(short = 't', long = "target-gt", allow_hyphen_values = true)]
    target_gt: Vec<String>,
    /// Print the usage.
    #[arg(short = 'h', long = "help", short_alias = '?', action = clap::ArgAction::Count)]
    help: u8,
}

/// An option occurrence, in command-line order: bcftools checks each as getopt reads it.
#[derive(Clone, Copy, Debug)]
enum Opt<'a> {
    Seed(&'a str),
    Filter,
    NewGt(&'a str),
    TargetGt(&'a str),
    Help,
}

/// The masks and values of `-n`, `-t`, and `-s` (bcftools' option loop of `init`).
#[derive(Debug, Default)]
struct Parsed {
    rand_seed: i64,
    new_mask: u32,
    custom_str: Option<String>,
    tgt_mask: u32,
    rand_frac: f64,
    binom: Option<(Vec<u8>, Cmp, f64)>,
}

/// Parses the options in command-line order, failing at the first bad one (or printing the
/// usage at `-h`) as bcftools' getopt loop does.
fn parse_options(
    args: &Args,
    matches: &ArgMatches,
    hdr_in: Option<&Header>,
) -> anyhow::Result<Parsed> {
    let events = super::in_order(
        matches,
        &[
            ("seed", &args.seed, Opt::Seed),
            ("include", &args.include, |_| Opt::Filter),
            ("exclude", &args.exclude, |_| Opt::Filter),
            ("new_gt", &args.new_gt, Opt::NewGt),
            ("target_gt", &args.target_gt, Opt::TargetGt),
        ],
        &[("help", args.help, Opt::Help)],
    );
    let mut p = Parsed::default();
    let mut filters = 0;
    for event in events {
        match event {
            Opt::Help => return Err(Exit::usage(usage())),
            Opt::Seed(s) => {
                let (v, end) = pasta::variant::number::strtol(s.as_bytes(), 0);
                if end != s.len() {
                    bail!("Could not parse: -s {s}");
                }
                p.rand_seed = v;
            }
            Opt::Filter => {
                filters += 1;
                if filters > 1 {
                    bail!(
                        "Error: only one -i or -e expression can be given, and they cannot be combined"
                    );
                }
            }
            Opt::NewGt(n) => parse_new_gt(&mut p, n)?,
            Opt::TargetGt(t) => parse_target_gt(&mut p, t, hdr_in)?,
        }
    }
    Ok(p)
}

/// One `-n` value (each replaces the earlier ones).
fn parse_new_gt(p: &mut Parsed, n: &str) -> anyhow::Result<()> {
    p.new_mask = 0;
    for (c, bit) in [
        ('.', gt::MISSING),
        ('0', gt::REF),
        ('m', gt::MINOR),
        ('M', gt::MAJOR),
        ('X', gt::X_VAF),
        ('p', gt::PHASED),
        ('u', gt::UNPHASED),
        ('i', gt::INVPHASE),
    ] {
        if n.contains(c) {
            p.new_mask |= bit;
        }
    }
    if n.starts_with("c:") {
        p.new_mask |= gt::CUSTOM;
        p.custom_str = Some(n.to_owned());
    }
    if p.new_mask == 0 {
        bail!("Unknown parameter to --new-gt: {n}");
    }
    Ok(())
}

/// One `-t` value (they add up).
fn parse_target_gt(p: &mut Parsed, t: &str, hdr_in: Option<&Header>) -> anyhow::Result<()> {
    match t {
        "." => p.tgt_mask |= gt::MISSING | gt::PARTIAL,
        "./x" => p.tgt_mask |= gt::PARTIAL,
        "./." => p.tgt_mask |= gt::MISSING,
        "a" => p.tgt_mask |= gt::ALL,
        "q" | "?" => p.tgt_mask |= gt::QUERY,
        _ => {}
    }
    if let Some(frac) = t.strip_prefix("r:") {
        let (v, end) = pasta::variant::number::strtod(frac.as_bytes(), 0);
        if end != frac.len() {
            bail!("Could not parse: -t {t}");
        }
        if v <= 0.0 || v >= 1.0 {
            bail!("Expected value between 0 and 1 with -t");
        }
        p.rand_frac = v;
        p.tgt_mask |= gt::RAND;
    }
    if let Some(i) = t.find('b') {
        p.binom = Some(parse_binom_expr(hdr_in, &t[i..])?);
        p.tgt_mask |= gt::BINOM;
    }
    if p.tgt_mask == 0 {
        bail!("Unknown parameter to --target-gt: {t}");
    }
    Ok(())
}

/// Bits of the target (`-t`) and new genotype (`-n`) masks (bcftools' `GT_*`).
mod gt {
    pub const MISSING: u32 = 1;
    pub const PARTIAL: u32 = 1 << 1;
    pub const REF: u32 = 1 << 2;
    pub const MAJOR: u32 = 1 << 3;
    pub const PHASED: u32 = 1 << 4;
    pub const UNPHASED: u32 = 1 << 5;
    pub const ALL: u32 = 1 << 6;
    pub const QUERY: u32 = 1 << 7;
    pub const BINOM: u32 = 1 << 8;
    pub const MINOR: u32 = 1 << 9;
    pub const CUSTOM: u32 = 1 << 10;
    pub const X_VAF: u32 = 1 << 11;
    pub const RAND: u32 = 1 << 12;
    pub const INVPHASE: u32 = 1 << 13;
}

/// An allele of a custom genotype (`-n c:GT`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CustomAllele {
    Index(i64),
    Minor,
    Major,
    XVaf,
    Missing,
}

/// The comparison of a binomial test (`-t b:TAG<VAL`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Cmp {
    Eq,
    Le,
    Ge,
    Lt,
    Gt,
}

impl Cmp {
    fn test(self, a: f64, b: f64) -> bool {
        match self {
            Self::Eq => a == b,
            Self::Le => a <= b,
            Self::Ge => a >= b,
            Self::Lt => a < b,
            Self::Gt => a > b,
        }
    }
}

/// glibc's `drand48`, which HTSlib's `hts_drand48` is on Linux.
#[derive(Clone, Copy, Debug)]
struct Drand48(u64);

impl Drand48 {
    /// `srand48(seed)`.
    fn new(seed: i64) -> Self {
        Self(((seed as u64 & 0xffff_ffff) << 16) | 0x330e)
    }

    /// `drand48()`.
    fn next(&mut self) -> f64 {
        self.0 = 0x5_DEEC_E66D_u64.wrapping_mul(self.0).wrapping_add(0xB) & ((1 << 48) - 1);
        self.0 as f64 / (1u64 << 48) as f64
    }
}

/// The plugin state (bcftools' `args_t`).
struct SetGt {
    tgt_mask: u32,
    new_mask: u32,
    new_gt: i32,
    filter: Option<(Filter, FilterLogic)>,
    custom: Vec<(CustomAllele, bool)>,
    m_allele: i32,
    big_m_allele: i32,
    x_vaf_allele: i32,
    rand_frac: f64,
    rng: Drand48,
    binom: Option<(Vec<u8>, Cmp, f64)>,
    nchanged: u64,
    gts: Vec<i32>,
    /// The genotypes widened to a custom genotype's larger ploidy.
    expanded: Vec<i32>,
    arr: Vec<i32>,
    iarr: Vec<i32>,
    xarr: Vec<i32>,
    smpl_pass: Vec<bool>,
}

/// bcftools' `_parse_binom_expr_error`.
fn binom_expr_error(s: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "Error parsing the expression: {s}\n\
Expected TAG CMP VAL, where\n\
\x20  TAG .. one of the format tags\n\
\x20  CMP .. operator, one of <, <=, >, >=\n\
\x20  VAL .. value\n\
For example:\n\
\x20  vtools +setGT in.vcf -- -t \"b:AD>1e-3\" -n 0\n\n"
    )
}

/// C's `isspace`.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Parses `b:TAG CMP VAL` (bcftools' `parse_binom_expr`), `s` starting at the `b`.
fn parse_binom_expr(hdr: Option<&Header>, s: &str) -> anyhow::Result<(Vec<u8>, Cmp, f64)> {
    let b = s.as_bytes();
    if b.get(1) != Some(&b':') {
        return Err(binom_expr_error(s));
    }
    let mut beg = 2;
    while beg < b.len() && is_space(b[beg]) {
        beg += 1;
    }
    if beg == b.len() {
        return Err(binom_expr_error(s));
    }
    let mut end = beg;
    while end < b.len() && !(is_space(b[end]) || matches!(b[end], b'<' | b'=' | b'>')) {
        end += 1;
    }
    if end == b.len() {
        return Err(binom_expr_error(s));
    }
    let tag = b[beg..end].to_vec();
    if let Some(hdr) = hdr {
        let ok = hdr
            .id(&tag)
            .is_some_and(|id| hdr.id_exists(HeaderRecordType::Format, id));
        if !ok {
            bail!(
                "The FORMAT tag \"{}\" is not present in the VCF",
                String::from_utf8_lossy(&tag)
            );
        }
    }
    while end < b.len() && is_space(b[end]) {
        end += 1;
    }
    if end == b.len() {
        return Err(binom_expr_error(s));
    }
    let rest = &b[end..];
    let (cmp, n) = if rest.starts_with(b"<=") {
        (Cmp::Le, 2)
    } else if rest.starts_with(b">=") {
        (Cmp::Ge, 2)
    } else if rest.starts_with(b"==") {
        (Cmp::Eq, 2)
    } else if rest.starts_with(b"<") {
        (Cmp::Lt, 1)
    } else if rest.starts_with(b">") {
        (Cmp::Gt, 1)
    } else if rest.starts_with(b"=") {
        (Cmp::Eq, 1)
    } else {
        return Err(binom_expr_error(s));
    };
    let mut beg = end + n;
    while beg < b.len() && is_space(b[beg]) {
        beg += 1;
    }
    if beg == b.len() {
        return Err(binom_expr_error(s));
    }
    let (val, mut vend) = pasta::variant::number::strtod(b, beg);
    while vend < b.len() && is_space(b[vend]) {
        vend += 1;
    }
    if vend != b.len() {
        return Err(binom_expr_error(s));
    }
    Ok((tag, cmp, val))
}

/// Parses a custom genotype, `c:GT` (bcftools' `init`), adding the m/M/X bits it uses.
fn parse_custom(gt_str: &str, new_mask: &mut u32) -> anyhow::Result<Vec<(CustomAllele, bool)>> {
    let bad = || anyhow::anyhow!("Could not parse the genotype: {gt_str}");
    let b = gt_str.as_bytes();
    let mut out: Vec<(CustomAllele, bool)> = Vec::new();
    let mut p = 2;
    while p < b.len() {
        let allele = match b[p] {
            b'm' => {
                *new_mask |= gt::MINOR;
                CustomAllele::Minor
            }
            b'M' => {
                *new_mask |= gt::MAJOR;
                CustomAllele::Major
            }
            b'X' => {
                *new_mask |= gt::X_VAF;
                CustomAllele::XVaf
            }
            b'.' => CustomAllele::Missing,
            b'/' | b'|' => {
                let Some(last) = out.last_mut() else {
                    return Err(bad());
                };
                last.1 = b[p] == b'|';
                p += 1;
                continue;
            }
            _ => {
                let (v, end) = pasta::variant::number::strtol(b, p);
                if end == p || (end < b.len() && b[end] != b'/' && b[end] != b'|') {
                    return Err(bad());
                }
                p = end - 1;
                CustomAllele::Index(v)
            }
        };
        out.push((allele, false));
        p += 1;
    }
    // The phasing sign comes before the allele: move it to the allele that follows.
    for i in (1..out.len()).rev() {
        out[i].1 = out[i - 1].1;
    }
    Ok(out)
}

/// Initializes the plugin (bcftools' `init`).
pub fn init(argv: &[String], hdrs: InitHeaders<'_>) -> anyhow::Result<Loaded> {
    // argv[0] is the input file, which getopt takes as the program name (for its messages).
    let prog = argv.first().map_or("setGT", String::as_str);
    let argv = shortopt::rewrite_strings(&Args::command(), prog, argv.get(1..).unwrap_or_default())
        .map_err(|msg| Exit::usage(format!("{msg}{}", usage())))?;
    let matches = Args::command()
        .try_get_matches_from(argv)
        .map_err(|_| Exit::usage(usage()))?;
    let args = Args::from_arg_matches(&matches).map_err(|_| Exit::usage(usage()))?;
    let hdr_in = hdrs.as_ref().map(|(h, _)| &**h);
    let Parsed {
        rand_seed,
        mut new_mask,
        custom_str,
        mut tgt_mask,
        rand_frac,
        binom,
    } = parse_options(&args, &matches, hdr_in)?;
    if new_mask == 0 {
        bail!("Expected -n option");
    }
    if tgt_mask == 0 {
        bail!("Expected -t option");
    }
    if tgt_mask & gt::RAND != 0 && tgt_mask == gt::RAND {
        tgt_mask |= gt::ALL;
    }
    let mut new_gt = 0;
    if new_mask & gt::MISSING != 0 {
        new_gt = 0;
    }
    if new_mask & gt::REF != 0 {
        new_gt = if new_mask & gt::PHASED != 0 { 3 } else { 2 };
    }
    let mut custom = Vec::new();
    if new_mask & gt::CUSTOM != 0 {
        let s = custom_str.as_deref().expect("c: sets the custom genotype");
        custom = parse_custom(s, &mut new_mask)?;
    }
    let filter_expr = args
        .include
        .last()
        .map(|e| (e, FilterLogic::Include))
        .or_else(|| args.exclude.last().map(|e| (e, FilterLogic::Exclude)));
    if filter_expr.is_some() && tgt_mask & gt::QUERY == 0 {
        bail!("Expected -tq with -i/-e");
    }
    if filter_expr.is_none() && tgt_mask & gt::QUERY != 0 {
        bail!("Expected -i/-e with -tq");
    }
    let Some((hdr, out)) = hdrs else {
        // bcftools goes on without a header and crashes.
        return Err(super::general_usage());
    };
    let filter = match filter_expr {
        Some((e, logic)) => Some((Filter::init(hdr, e)?, logic)),
        None => None,
    };
    let fmt_exists = |name: &[u8]| {
        hdr.id(name)
            .is_some_and(|id| hdr.id_exists(HeaderRecordType::Format, id))
    };
    if !fmt_exists(b"GT") {
        let _ = out.append(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    }
    if new_mask & gt::X_VAF != 0 && !fmt_exists(b"AD") {
        bail!(
            "Error: the FORMAT/AD annotation does exist, cannot run with --new-gt {}",
            custom_str.as_deref().unwrap_or("(null)")
        );
    }
    Ok(Loaded {
        plugin: Box::new(SetGt {
            tgt_mask,
            new_mask,
            new_gt,
            filter,
            custom,
            m_allele: 0,
            big_m_allele: 0,
            x_vaf_allele: 0,
            rand_frac,
            rng: Drand48::new(rand_seed),
            binom,
            nchanged: 0,
            gts: Vec::new(),
            expanded: Vec::new(),
            arr: Vec::new(),
            iarr: Vec::new(),
            xarr: Vec::new(),
            smpl_pass: Vec::new(),
        }),
        output: Output::Records,
        max_unpack: Unpack::empty(),
    })
}

/// `bcf_gt_unphased`.
fn unphased(a: i32) -> i32 {
    (a + 1) << 1
}

/// `bcf_gt_phased`.
fn phased(a: i32) -> i32 {
    ((a + 1) << 1) | 1
}

/// `bcf_gt_is_missing`.
fn is_missing(v: i32) -> bool {
    v >> 1 == 0
}

/// `bcf_gt_allele`.
fn allele(v: i32) -> i32 {
    (v >> 1) - 1
}

/// Phases a genotype (bcftools' `phase_gt`); returns the number of changed alleles.
fn phase_gt(ptr: &mut [i32]) -> u64 {
    let mut changed = 0;
    for v in ptr {
        if *v == INT32_VECTOR_END {
            break;
        }
        if *v & 1 != 0 {
            continue;
        }
        *v = phased(allele(*v));
        changed += 1;
    }
    changed
}

/// Unphases a genotype and sorts its alleles (bcftools' `unphase_gt`).
fn unphase_gt(ptr: &mut [i32]) -> u64 {
    let mut changed = 0;
    let mut j = 0;
    while j < ptr.len() {
        if ptr[j] == INT32_VECTOR_END {
            break;
        }
        if ptr[j] & 1 != 0 {
            ptr[j] = unphased(allele(ptr[j]));
            changed += 1;
        }
        j += 1;
    }
    // Insertion sort, as bcftools sorts the few alleles.
    for k in 1..j {
        let x = ptr[k];
        let mut l = k;
        while l > 0 && ptr[l - 1] > x {
            ptr[l] = ptr[l - 1];
            l -= 1;
        }
        ptr[l] = x;
    }
    changed
}

/// Swaps the alleles of a diploid genotype (bcftools' `invert_phase_gt`); `ngts` is the
/// record's largest ploidy.
fn invert_phase_gt(ptr: &mut [i32], ngts: usize) -> u64 {
    if ngts != 2 || ptr[0] == INT32_VECTOR_END || ptr[1] == INT32_VECTOR_END {
        return 0;
    }
    let gt = allele(ptr[1]);
    ptr[1] = if ptr[1] & 1 != 0 {
        phased(allele(ptr[0]))
    } else {
        unphased(allele(ptr[0]))
    };
    ptr[0] = unphased(gt);
    2
}

/// Sets every allele of a genotype (bcftools' `set_gt`).
fn set_gt(ptr: &mut [i32], value: i32) -> u64 {
    let mut changed = 0;
    for v in ptr {
        if *v == INT32_VECTOR_END {
            break;
        }
        if *v != value {
            changed += 1;
        }
        *v = value;
    }
    changed
}

/// Which samples a filter selects after `-i`/`-e` (bcftools' `GT_QUERY` branch).
enum Selection {
    /// Leave the record as it is.
    Nothing,
    /// Every sample.
    All,
    /// The samples marked in `smpl_pass`.
    Marked,
}

/// The alleles `m`, `M`, and `X` of a custom genotype stand for at the current record and
/// sample.
#[derive(Clone, Copy, Debug)]
struct CustomAlleles {
    minor: i32,
    major: i32,
    x_vaf: i32,
}

/// Sets a genotype to the custom one (bcftools' `set_gt_custom`).
fn set_gt_custom(
    custom: &[(CustomAllele, bool)],
    alleles: CustomAlleles,
    ptr: &mut [i32],
    nals: usize,
) -> u64 {
    let nals = nals as i64;
    let mut changed = false;
    for (i, &(a, is_phased)) in custom.iter().enumerate() {
        let new = match a {
            CustomAllele::Minor => i64::from(alleles.minor),
            CustomAllele::Major => i64::from(alleles.major),
            CustomAllele::XVaf => {
                if is_missing(alleles.x_vaf) {
                    nals
                } else {
                    i64::from(allele(alleles.x_vaf))
                }
            }
            CustomAllele::Missing => nals,
            CustomAllele::Index(n) => n,
        };
        let new = if new >= nals {
            0
        } else if is_phased {
            phased(new as i32)
        } else {
            unphased(new as i32)
        };
        if ptr[i] == new {
            continue;
        }
        ptr[i] = new;
        changed = true;
    }
    for v in &mut ptr[custom.len()..] {
        if !changed && *v != INT32_VECTOR_END {
            changed = true;
        }
        *v = INT32_VECTOR_END;
    }
    u64::from(changed)
}

impl SetGt {
    /// `random_draw`: whether the genotype is skipped.
    fn skip_random(&mut self) -> bool {
        self.tgt_mask & gt::RAND != 0 && self.rng.next() > self.rand_frac
    }

    /// Changes the genotype of sample `i` as `-n` asks.
    #[inline]
    fn apply(&mut self, i: usize, ngts: usize, nals: usize) -> u64 {
        let m = self.new_mask;
        let gt = &mut self.gts[i * ngts..(i + 1) * ngts];
        if m & gt::UNPHASED != 0 {
            unphase_gt(gt)
        } else if m == gt::PHASED {
            phase_gt(gt)
        } else if m & gt::CUSTOM != 0 {
            if m & gt::X_VAF != 0 {
                self.x_vaf_allele = self.xarr[i];
            }
            let alleles = CustomAlleles {
                minor: self.m_allele,
                major: self.big_m_allele,
                x_vaf: self.x_vaf_allele,
            };
            set_gt_custom(&self.custom, alleles, gt, nals)
        } else if m & gt::X_VAF != 0 {
            set_gt(gt, self.xarr[i])
        } else if m & gt::INVPHASE != 0 {
            invert_phase_gt(gt, ngts)
        } else {
            set_gt(gt, self.new_gt)
        }
    }

    fn calc_ac(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<()> {
        // HTSlib logs a bad allele and exits with status 1.
        let ok = utils::calc_ac(hdr, rec, &mut self.arr, AcSource::Genotypes)
            .map_err(|_| Exit::status(1))?;
        if !ok {
            bail!(
                "Could not calculate allele count at {}:{}",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1
            );
        }
        Ok(())
    }

    /// Finds the major and minor alleles and the alleles with the most reads, as `-n` needs.
    /// Returns `false` when FORMAT/AD is unusable for `X` (the record is left as it is).
    fn prepare_alleles(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        let nals = rec.allele_count();
        if self.new_mask & gt::MAJOR != 0 {
            self.calc_ac(hdr, rec)?;
            let mut max_ac = -1;
            let mut major = -1;
            for (i, &ac) in self.arr.iter().enumerate().take(nals) {
                if ac > max_ac {
                    max_ac = ac;
                    major = i as i32;
                }
            }
            self.new_gt = if self.new_mask & gt::PHASED != 0 {
                phased(major)
            } else {
                unphased(major)
            };
            if self.new_mask & gt::CUSTOM != 0 {
                self.big_m_allele = major;
            }
        }
        if self.new_mask & gt::MINOR != 0 {
            self.calc_ac(hdr, rec)?;
            let arr = &self.arr;
            let mut imax = 0;
            for i in 1..nals {
                if arr[imax] < arr[i] {
                    imax = i;
                }
            }
            let mut imax2 = if imax > 0 { 0 } else { usize::from(nals > 1) };
            for i in 0..nals {
                if i != imax && arr[imax2] < arr[i] {
                    imax2 = i;
                }
            }
            self.new_gt = if self.new_mask & gt::PHASED != 0 {
                phased(imax2 as i32)
            } else {
                unphased(imax2 as i32)
            };
            if self.new_mask & gt::CUSTOM != 0 {
                self.m_allele = imax2 as i32;
            }
        }
        if self.new_mask & gt::X_VAF != 0 {
            let nsmpl = rec.sample_count();
            match rec.get_format_i32(hdr, b"AD", &mut self.xarr) {
                Ok(n) if n == nals * nsmpl => {}
                _ => return Ok(false),
            }
            for i in 0..nsmpl {
                let src = &self.xarr[i * nals..(i + 1) * nals];
                let mut jmax: Option<usize> = None;
                for (j, &v) in src.iter().enumerate() {
                    if v == INT32_VECTOR_END {
                        break;
                    }
                    if v == pasta::variant::encoding::INT32_MISSING {
                        continue;
                    }
                    if jmax.is_none_or(|m| src[m] < v) {
                        jmax = Some(j);
                    }
                }
                self.xarr[i] = jmax.map_or(0, |j| unphased(j as i32));
            }
        }
        Ok(true)
    }

    /// The samples `-i/-e` select (bcftools' `GT_QUERY` branch), in `smpl_pass`.
    fn query(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<Selection> {
        let (filter, logic) = self.filter.as_mut().expect("-t q requires -i/-e");
        let mut pass_site = filter.test(hdr, rec)?;
        let samples = filter.sample_pass();
        self.smpl_pass.clear();
        if let Some(s) = samples {
            self.smpl_pass.extend_from_slice(s);
        }
        let has_samples = samples.is_some();
        if *logic == FilterLogic::Exclude {
            // -i can include a site but exclude a sample, -e exclude a site but include a sample.
            if pass_site {
                if !has_samples {
                    return Ok(Selection::Nothing);
                }
                pass_site = false;
                for p in &mut self.smpl_pass {
                    *p = !*p;
                    pass_site |= *p;
                }
                if !pass_site {
                    return Ok(Selection::Nothing);
                }
            } else if has_samples {
                self.smpl_pass.fill(true);
            }
        } else if !pass_site {
            return Ok(Selection::Nothing);
        }
        Ok(if has_samples {
            Selection::Marked
        } else {
            Selection::All
        })
    }

    /// The binomial-test branch: heterozygous diploid genotypes whose allele counts in the
    /// tag pass the test.
    fn binom_branch(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        ngts: usize,
        nbinom: usize,
    ) -> anyhow::Result<u64> {
        let mut sample_mask = None;
        if let Some((filter, logic)) = &mut self.filter {
            filter.test(hdr, rec)?;
            if let Some(s) = filter.sample_pass() {
                self.smpl_pass.clear();
                self.smpl_pass.extend_from_slice(s);
                sample_mask = Some(*logic);
            }
        }
        let (cmp, val) = self
            .binom
            .as_ref()
            .map(|(_, cmp, val)| (*cmp, *val))
            .expect("the binomial test is set");
        let nals = rec.allele_count();
        let mut changed = 0;
        for i in 0..rec.sample_count() {
            if let Some(logic) = sample_mask
                && self.smpl_pass[i] != (logic == FilterLogic::Include)
            {
                continue;
            }
            let (a, b) = (self.gts[i * ngts], self.gts[i * ngts + 1]);
            if is_missing(a) || is_missing(b) || b == INT32_VECTOR_END || a == b {
                continue;
            }
            let (ia, ib) = (allele(a), allele(b));
            if ia as usize >= nbinom || ib as usize >= nbinom {
                let tag = self.binom.as_ref().map_or(&[][..], |(tag, _, _)| tag);
                bail!(
                    "The sample {} has incorrect number of {} fields at {}:{}",
                    String::from_utf8_lossy(&hdr.samples()[i]),
                    String::from_utf8_lossy(tag),
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1
                );
            }
            let prob = binom_two_sided(
                self.iarr[i * nbinom + ia as usize],
                self.iarr[i * nbinom + ib as usize],
                0.5,
            );
            if !cmp.test(prob, val) || self.skip_random() {
                continue;
            }
            changed += self.apply(i, ngts, nals);
        }
        Ok(changed)
    }
}

impl Plugin for SetGt {
    fn process(&mut self, hdr: &Header, out: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        let nsmpl = rec.sample_count();
        if nsmpl == 0 {
            return Ok(true);
        }
        let mut ngts = rec.get_genotypes(hdr, &mut self.gts).unwrap_or_else(|_| {
            self.gts.clear();
            0
        }) / nsmpl;
        if self.new_mask & gt::CUSTOM != 0 && self.custom.len() > ngts {
            // The custom genotype increases the ploidy: expand the array.
            let p = self.custom.len();
            self.expanded.clear();
            self.expanded.resize(p * nsmpl, INT32_VECTOR_END);
            for i in 0..nsmpl {
                self.expanded[i * p..i * p + ngts]
                    .copy_from_slice(&self.gts[i * ngts..(i + 1) * ngts]);
            }
            std::mem::swap(&mut self.gts, &mut self.expanded);
            ngts = p;
        }
        if !self.prepare_alleles(hdr, rec)? {
            return Ok(true);
        }
        let mut nbinom = 0;
        if let Some((tag, _, _)) = &self.binom {
            nbinom = rec.get_format_i32(hdr, tag, &mut self.iarr).unwrap_or(0) / nsmpl;
        }
        let nals = rec.allele_count();
        let mut changed = 0;
        if nbinom > 0 && ngts >= 2 {
            changed = self.binom_branch(hdr, rec, ngts, nbinom)?;
        } else if self.tgt_mask & gt::QUERY != 0 {
            let selection = self.query(hdr, rec)?;
            if let Selection::Nothing = selection {
                return Ok(true);
            }
            for i in 0..nsmpl {
                if matches!(selection, Selection::Marked) && !self.smpl_pass[i] {
                    continue;
                }
                if self.skip_random() {
                    continue;
                }
                changed += self.apply(i, ngts, nals);
            }
        } else {
            for i in 0..nsmpl {
                let do_set = self.tgt_mask & gt::ALL != 0 || {
                    let gt = &self.gts[i * ngts..(i + 1) * ngts];
                    let ploidy = gt.iter().take_while(|&&v| v != INT32_VECTOR_END).count();
                    let nmiss = gt[..ploidy].iter().filter(|&&v| is_missing(v)).count();
                    (self.tgt_mask & gt::PARTIAL != 0 && nmiss > 0)
                        || (self.tgt_mask & gt::MISSING != 0 && ploidy == nmiss)
                };
                if !do_set || self.skip_random() {
                    continue;
                }
                changed += self.apply(i, ngts, nals);
            }
        }
        self.nchanged += changed;
        if changed > 0 {
            rec.update_genotypes(out, &self.gts[..ngts * nsmpl])
                .map_err(|_| {
                    anyhow::anyhow!(
                        "Error: failed to update genotypes at {}:{}",
                        String::from_utf8_lossy(rec.chrom_safe(hdr)),
                        rec.pos() + 1
                    )
                })?;
        }
        Ok(true)
    }

    fn destroy(&mut self) -> anyhow::Result<()> {
        eprintln!("Filled {} alleles", self.nchanged);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drand48_matches_glibc() {
        // glibc: srand48(0) then drand48() twice.
        let mut r = Drand48::new(0);
        assert_eq!(r.next(), 0.170_828_036_106_289_72);
        assert_eq!(r.next(), 0.749_901_980_484_963_8);
    }

    #[test]
    fn custom_genotypes_take_the_phase_of_their_separator() {
        let mut mask = 0;
        let c = parse_custom("c:0|1/2", &mut mask).expect("valid");
        assert_eq!(
            c,
            [
                (CustomAllele::Index(0), true),
                (CustomAllele::Index(1), true),
                (CustomAllele::Index(2), false)
            ]
        );
        let c = parse_custom("c:m/M", &mut mask).expect("valid");
        assert_eq!(c.len(), 2);
        assert_eq!(mask, gt::MINOR | gt::MAJOR);
        assert!(parse_custom("c:/1", &mut mask).is_err());
        assert!(parse_custom("c:1x", &mut mask).is_err());
    }

    #[test]
    fn binomial_expressions() {
        let (tag, cmp, val) = parse_binom_expr(None, "b:AD<=1e-3").expect("valid");
        assert_eq!((tag.as_slice(), cmp, val), (&b"AD"[..], Cmp::Le, 1e-3));
        assert!(parse_binom_expr(None, "b:AD").is_err());
        assert!(parse_binom_expr(None, "b:AD<1x").is_err());
    }
}
