//! `vtools merge`: merge VCF/BCF files of different samples into one multi-sample file
//! (bcftools merge).
//!
//! The files are read in step through the synced reader. At each position, the records of all
//! files are grouped into output lines by their alleles and the `-m` rules (`maux`), then each
//! line's fields are merged: alleles normalized across files, Number=A/G/R values remapped,
//! INFO combined by the `-i` rules (`site`), and FORMAT values laid out sample by sample
//! (`format`). gVCF blocks are cut where records of other files interrupt them (`gvcf`).

mod alleles;
mod format;
mod gvcf;
mod header;
mod maux;
mod rules;
mod site;

use std::{
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::bail;
use clap::ArgMatches;
use pasta::variant::{
    Header, Record, Unpack,
    io::{OutputType, WriteError, Writer},
    readlist::{self, ListSource},
    regidx::{self, RegIdx},
    synced::{IndexRequirement, Overlap, SyncedReaders},
};

use self::{
    gvcf::Flush,
    header::DuplicateSamples,
    maux::{Inputs, Maux},
    rules::{Defaults, InfoRule, MissingRule},
};
use super::{
    call::vcmp::Vcmp,
    head::apply_verbosity,
    norm::{emit::Emitter, reference::Reference},
};
use crate::options::later;
use crate::{
    output::{self, OutputSetting},
    raw_input::strerror,
    version,
};

const USAGE: &str = "
About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.
         Note that only records from different files can be merged, never from the same file. For
         \"vertical\" merge take a look at \"vtools norm\" instead.
Usage:   vtools merge [options] <A.vcf.gz> <B.vcf.gz> [...]

Options:
        --force-no-index              Merge unindexed files, synonymous to --no-index
        --force-samples               Resolve duplicate sample names
        --force-single                Run even if there is only one file on input
        --print-header                Print only the merged header and exit
        --use-header FILE             Use the provided header
    -0  --missing-to-ref              Assume genotypes at missing sites are 0/0
    -f, --apply-filters LIST          Require at least one of the listed FILTER strings (e.g. \"PASS,.\")
    -F, --filter-logic x|+            Remove filters if some input is PASS (\"x\"), or apply all filters (\"+\") [+]
    -g, --gvcf -|REF.FA               Merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,MIN_DP:min,I16:sum,IDV:max,IMF:max -M PL:max,AD:0
    -i, --info-rules TAG:METHOD,..    Rules for merging INFO fields (method is one of sum,avg,min,max,join) or \"-\" to turn off the default [DP:sum,DP4:sum]
    -l, --file-list FILE              Read file names from the file
    -L, --local-alleles INT           If more than INT alt alleles are encountered, drop FMT/PL and output LAA+LPL instead; 0=unlimited [0]
    -m, --merge STRING[*|**]          Allow multiallelic records for snps,indels,both,snp-ins-del,all,exact,none,id,*,**; see man page for details [both]
    -M, --missing-rules TAG:METHOD    Rules for replacing missing values in numeric vectors (.,0,max) when unknown allele <*> is not present [.]
        --no-index                    Merge unindexed files, the same chromosomal order is required and -r/-R are not allowed
        --no-version                  Do not append version and command line to the header
    -o, --output FILE                 Write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9]    u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -r, --regions REGION              Restrict to comma-separated list of regions
    -R, --regions-file FILE           Restrict to regions listed in a file
        --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
        --threads INT                 Use multithreading with INT worker threads [0]
    -v, --verbosity INT               Verbosity level
    -W, --write-index[=FMT]           Automatically index the output files [off]

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF files, bgzipped and indexed unless --no-index is given.
    inputs: Vec<String>,

    /// Merge unindexed files, synonymous to --no-index.
    #[arg(long = "force-no-index")]
    force_no_index: bool,

    /// Resolve duplicate sample names.
    #[arg(long = "force-samples")]
    force_samples: bool,

    /// Run even if there is only one file on input.
    #[arg(long = "force-single")]
    force_single: bool,

    /// Print only the merged header and exit.
    #[arg(long = "print-header")]
    print_header: bool,

    /// Use the provided header.
    #[arg(long = "use-header", value_name = "FILE", allow_hyphen_values = true)]
    use_header: Option<String>,

    /// Assume genotypes at missing sites are 0/0.
    #[arg(short = '0', long = "missing-to-ref")]
    missing_to_ref: bool,

    /// Require at least one of the listed FILTER strings (e.g. "PASS,.").
    #[arg(
        short = 'f',
        long = "apply-filters",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    apply_filters: Option<String>,

    /// Remove filters if some input is PASS ("x"), or apply all filters ("+") [+].
    #[arg(
        short = 'F',
        long = "filter-logic",
        value_name = "x|+",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    filter_logic: Vec<String>,

    /// Merge gVCF blocks, INFO/END tag is expected; the reference fills in REF.
    #[arg(
        short = 'g',
        long = "gvcf",
        value_name = "-|REF.FA",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    gvcf: Vec<String>,

    /// Rules for merging INFO fields (sum,avg,min,max,join) or "-" to turn off the default.
    #[arg(
        short = 'i',
        long = "info-rules",
        value_name = "TAG:METHOD,..",
        allow_hyphen_values = true
    )]
    info_rules: Option<String>,

    /// Read file names from the file.
    #[arg(
        short = 'l',
        long = "file-list",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    file_list: Option<String>,

    /// If more than INT alt alleles are encountered, drop FMT/PL and output LAA+LPL instead.
    #[arg(
        short = 'L',
        long = "local-alleles",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    local_alleles: Vec<String>,

    /// Allow multiallelic records for snps,indels,both,snp-ins-del,all,exact,none,id,*,** [both].
    #[arg(
        short = 'm',
        long = "merge",
        value_name = "STRING",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    merge: Vec<String>,

    /// Rules for replacing missing values in numeric vectors (.,0,max).
    #[arg(
        short = 'M',
        long = "missing-rules",
        value_name = "TAG:METHOD",
        allow_hyphen_values = true
    )]
    missing_rules: Option<String>,

    /// Merge unindexed files; the same chromosomal order is required.
    #[arg(long = "no-index")]
    no_index: bool,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Write output to a file [standard output].
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "u|b|v|z[0-9]",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    output_type: Vec<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(
        short = 'r',
        long = "regions",
        value_name = "REGION",
        allow_hyphen_values = true
    )]
    regions: Option<String>,

    /// Restrict to regions listed in a file.
    #[arg(
        short = 'R',
        long = "regions-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    regions_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(
        long = "regions-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    regions_overlap: Vec<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,

    /// Automatically index the output files [off]: csi (default) or tbi.
    #[arg(
        short = 'W',
        long = "write-index",
        value_name = "FMT",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "",
        action = clap::ArgAction::Append
    )]
    write_index: Vec<String>,

    /// Print help.
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Count)]
    help: u8,
}

/// bcftools merge's `-m` modes (`args->collapse`, HTSlib's `COLLAPSE_*` bits). bcftools
/// compares the bits both as sets and for equality, so they are kept as bits.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Collapse(u32);

impl Collapse {
    /// SNPs with SNPs (`COLLAPSE_SNPS`).
    const SNPS: Self = Self(1);
    /// Indels with indels (`COLLAPSE_INDELS`).
    const INDELS: Self = Self(1 << 1);
    /// Anything with anything (`COLLAPSE_ANY`).
    const ANY: Self = Self(1 << 2);
    /// Only identical allele sets (`BCF_SR_PAIR_EXACT`).
    const EXACT: Self = Self(1 << 6);
    /// SNPs, insertions, and deletions separately (`COLLAPSE_SNP_INS_DEL`).
    const SNP_INS_DEL: Self = Self(1 << 10);
    /// SNPs with SNPs and indels with indels (`COLLAPSE_BOTH`).
    const BOTH: Self = Self(1 | 1 << 1);

    /// No merging of different alleles (`COLLAPSE_NONE`).
    const fn empty() -> Self {
        Self(0)
    }

    /// Whether all bits of `other` are set.
    const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }

    /// Whether any bit of `other` is set.
    const fn intersects(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

impl std::ops::BitOr for Collapse {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

impl std::ops::BitOrAssign for Collapse {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

/// `-F`: how the FILTERs of the merged records combine.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum FilterLogic {
    /// Apply all filters (`+`).
    #[default]
    Add,
    /// PASS when some record passes (`x`).
    Remove,
}

/// The options the merging needs (from bcftools' `args_t`).
#[derive(Debug, Default)]
struct Opts {
    collapse: Collapse,
    merge_by_id: bool,
    do_gvcf: bool,
    filter_logic: FilterLogic,
    missing_to_ref: bool,
    /// `-L`: the most ALT alleles output without localizing; 0 for no limit.
    local_alleles: usize,
    /// AC/AN have INFO rules, so they are merged rather than recomputed.
    keep_ac_an: bool,
    /// `-m *`: trim an unobserved allele at variant sites (1) or at all sites (2).
    trim_star_allele: u8,
}

/// An error already reported, which ends the program with status 1 (bcftools' `exit(1)`).
#[derive(Debug)]
struct ExitOne;

impl std::fmt::Display for ExitOne {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("exit 1")
    }
}

impl std::error::Error for ExitOne {}

/// Marks dictionary IDs at one site without clearing a table per site.
#[derive(Debug, Default)]
struct Seen {
    stamp: Vec<u32>,
    val: Vec<i64>,
    now: u32,
}

impl Seen {
    /// Forgets everything.
    fn clear(&mut self) {
        self.now = self.now.wrapping_add(1);
        if self.now == 0 {
            self.stamp.fill(0);
            self.now = 1;
        }
    }

    /// Returns the value stored for an ID at this site.
    #[inline]
    fn get(&self, id: u32) -> Option<i64> {
        let i = id as usize;
        (self.stamp.get(i) == Some(&self.now)).then(|| self.val[i])
    }

    /// Stores a value for an ID.
    #[inline]
    fn set(&mut self, id: u32, v: i64) {
        let i = id as usize;
        if i >= self.stamp.len() {
            self.stamp.resize(i + 1, 0);
            self.val.resize(i + 1, 0);
        }
        self.stamp[i] = self.now;
        self.val[i] = v;
    }
}

/// An output dictionary ID looked up for a reader's ID.
#[derive(Clone, Copy, Debug, Default)]
enum OutId {
    /// Not looked up yet.
    #[default]
    Unknown,
    /// The output header lacks the name.
    Absent,
    Id(u32),
}

/// Translates each reader's dictionary IDs to the output header's, by name.
#[derive(Debug, Default)]
struct IdMaps {
    /// Per reader, by dictionary ID.
    maps: Vec<Vec<OutId>>,
    /// The output IDs of INFO/END, AN, and AC.
    end: Option<u32>,
    an: Option<u32>,
    ac: Option<u32>,
}

impl IdMaps {
    /// Returns the output ID of reader `i`'s dictionary ID `id`.
    #[inline]
    fn out_id(&mut self, i: usize, id: i32, hdr: &Header, out: &Header) -> Option<u32> {
        let map = &mut self.maps[i];
        let k = id as usize;
        if k >= map.len() {
            map.resize(k + 1, OutId::Unknown);
        }
        if let OutId::Unknown = map[k] {
            map[k] = hdr
                .id_name(id as u32)
                .and_then(|name| out.id(name))
                .map_or(OutId::Absent, OutId::Id);
        }
        match map[k] {
            OutId::Id(id) => Some(id),
            OutId::Unknown | OutId::Absent => None,
        }
    }
}

/// The merge (bcftools' `args_t` while merging).
struct Merger {
    inp: Inputs,
    out_hdr: Header,
    out: Record,
    ma: Maux,
    opts: Opts,
    rules: Vec<InfoRule>,
    missing_rules: Vec<(Vec<u8>, MissingRule)>,
    regs: Option<RegIdx>,
    fai: Option<Reference>,
    writer: Writer,
    /// Writes the output lines, formatting or encoding them on worker threads with `--threads`.
    emitter: Option<Emitter>,
    output_name: String,
    vcmp: Vcmp,
    ids: IdMaps,
    seen: Seen,
    fmt: format::FormatState,
    site: site::SiteState,
    scratch_i: Vec<i32>,
    scratch_f: Vec<f32>,
    /// The alleles `-m *` removes.
    trim_mask: Vec<bool>,
}

/// One `-m` value applied to the mode (bcftools' option loop): a prefix of a mode name, with
/// `*` or `**` to trim unobserved alleles.
fn parse_merge_mode(value: &str, opts: &mut Opts) -> anyhow::Result<()> {
    let mut s = value.as_bytes();
    for _ in 0..2 {
        if let Some((b'*', rest)) = s.split_last() {
            opts.trim_star_allele += 1;
            s = rest;
        }
    }
    if let Some((b',', rest)) = s.split_last() {
        s = rest;
    }
    let is = |name: &str| name.as_bytes().starts_with(s);
    opts.collapse = Collapse::empty();
    if is("snp-ins-del") {
        opts.collapse = Collapse::SNP_INS_DEL | Collapse::SNPS;
    } else if is("snps") {
        opts.collapse |= Collapse::SNPS;
    } else if is("indels") {
        opts.collapse |= Collapse::INDELS;
    } else if is("id") {
        opts.collapse = Collapse::empty();
        opts.merge_by_id = true;
    } else if is("any") || is("all") {
        opts.collapse |= Collapse::ANY;
    } else if is("both") {
        opts.collapse |= Collapse::BOTH;
    } else if is("none") {
        opts.collapse = Collapse::empty();
    } else if is("exact") {
        opts.collapse = Collapse::EXACT;
    } else {
        bail!("The -m type \"{value}\" is not recognised.");
    }
    Ok(())
}

/// C's `strtol(s, &end, 10)` requiring the whole string to be read.
fn strtol_full(s: &str) -> Option<i64> {
    let t = s.trim_start();
    if t.is_empty() {
        return s.is_empty().then_some(0);
    }
    t.parse::<i64>().ok()
}

pub fn run(mut args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    // Every -v is checked, in order, as getopt hands each one to bcftools.
    for v in &args.verbosity {
        apply_verbosity(v)?;
    }
    if args.help > 0 {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    let mut opts = Opts {
        collapse: Collapse::BOTH,
        ..Opts::default()
    };
    for l in &args.local_alleles {
        let Some(n) = strtol_full(l) else {
            bail!("Could not parse argument: --local-alleles {l}");
        };
        if n < 1 {
            bail!(
                "Error: \"--local-alleles {l}\" makes no sense, expected value bigger or equal than 1"
            );
        }
        opts.local_alleles = n as usize;
    }
    for f in &args.filter_logic {
        opts.filter_logic = match f.as_str() {
            "+" => FilterLogic::Add,
            "x" => FilterLogic::Remove,
            _ => bail!("Filter logic not recognised: {f}"),
        };
    }
    opts.missing_to_ref = args.missing_to_ref;
    let mut fai = None;
    for g in &args.gvcf {
        opts.do_gvcf = true;
        if g != "-" {
            match Reference::open(Path::new(g)) {
                Ok(r) => fai = Some(r),
                Err(_) => bail!("Failed to load the fai index: {g}"),
            }
        }
    }
    for m in &args.merge {
        parse_merge_mode(m, &mut opts)?;
    }
    let mut settings = Vec::new();
    if let Some(idx) = matches.indices_of("output_type") {
        let mut ordered: Vec<(usize, &String)> = idx.zip(&args.output_type).collect();
        ordered.sort_by_key(|(i, _)| *i);
        for (_, t) in ordered {
            settings.push(OutputSetting::Type(t.clone()));
            output::requested_type(&settings)?;
        }
    }
    let regions_overlap = args
        .regions_overlap
        .iter()
        .try_fold(Overlap::Record, |_, s| match Overlap::parse(s) {
            Some(o) => Ok(o),
            None => bail!("Could not parse: --regions-overlap {s}"),
        })?;
    // strtol: a negative count fails when the thread pool is made.
    let threads = args
        .threads
        .as_deref()
        .map_or(Some(0), |t| usize::try_from(output::thread_count(t)).ok());
    let write_index = args
        .write_index
        .iter()
        .try_fold(None, |_, w| output::parse_write_index(w).map(Some))?;

    if args.inputs.is_empty() && args.file_list.is_none() {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    // -R makes the regions a file even when a later -r gives the list (bcftools' flag).
    let regions_source = if args.regions_file.is_some() {
        ListSource::File
    } else {
        ListSource::Inline
    };
    let regions = later(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    )
    .map(|(spec, _)| spec);
    let no_index = args.no_index || args.force_no_index;

    let mut sr = SyncedReaders::new();
    sr.set_apply_filters(args.apply_filters.as_deref().map(str::as_bytes));
    let mut regs = None;
    if no_index {
        if regions.is_some() {
            bail!("Error: cannot combine --no-index with -r/-R");
        }
        sr.set_index_requirement(IndexRequirement::AllowMissing);
    } else {
        sr.set_index_requirement(IndexRequirement::Required);
    }
    if let Some(spec) = regions {
        sr.set_regions_overlap(regions_overlap);
        if sr
            .set_regions_spec(spec.as_bytes(), regions_source)
            .is_err()
        {
            bail!("Failed to read the regions: {spec}");
        }
        regs = Some(load_regidx(spec, regions_source)?);
    }
    let Some(threads) = threads.filter(|&t| output::threads_start(t)) else {
        bail!("Failed to create threads");
    };
    let mut files = std::mem::take(&mut args.inputs);
    if let Some(list) = &args.file_list {
        let Ok(names) = readlist::read_lines_or_inline(list.as_bytes()) else {
            bail!("Failed to read from {list}");
        };
        files.extend(
            names
                .iter()
                .map(|n| String::from_utf8_lossy(n).into_owned()),
        );
    }
    sr.set_shared_threads(threads, files.len());
    for f in &files {
        if let Err(e) = sr.add_reader(f) {
            bail!("Failed to open {f}: {e}");
        }
    }
    if sr.reader_count() == 0 {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    if sr.reader_count() == 1 && !args.force_single {
        bail!(
            "Expected two or more files to merge, got only one. Use --force-single to proceed anyway"
        );
    }
    merge_vcf(&args, sr, opts, regs, fai, &settings, threads, write_index)
}

/// The regions as an index, for trimming gVCF blocks and checking merged lines (bcftools'
/// `regidx_init` of `-r`/`-R`).
fn load_regidx(spec: &str, source: ListSource) -> anyhow::Result<RegIdx> {
    let idx = match source {
        ListSource::File => RegIdx::from_file(Path::new(spec)).ok(),
        ListSource::Inline => {
            let mut idx = RegIdx::new();
            idx.insert_list(spec.as_bytes(), b',', regidx::parse_reg)
                .ok()
                .map(|()| idx)
        }
    };
    match idx {
        Some(idx) => Ok(idx),
        None => bail!("Could not parse the regions: {spec}"),
    }
}

/// bcftools' `merge_vcf`: the output header, then the merge position by position.
#[expect(
    clippy::too_many_arguments,
    reason = "the parsed options of main_vcfmerge"
)]
fn merge_vcf(
    args: &Args,
    sr: SyncedReaders,
    mut opts: Opts,
    regs: Option<RegIdx>,
    fai: Option<Reference>,
    settings: &[OutputSetting],
    threads: usize,
    write_index: Option<output::WriteIndex>,
) -> anyhow::Result<ExitCode> {
    let output_name = args
        .output
        .as_deref()
        .map_or_else(|| "-".to_owned(), |p| p.to_string_lossy().into_owned());
    let (out_type, level) = output::resolve(settings, args.output.as_deref())?;
    let writer = create_writer(
        args.output.as_deref(),
        out_type,
        level,
        threads,
        &output_name,
    )?;
    let (mut out_hdr, mut writer) = output_header(args, &sr, &opts, writer)?;
    let defaults = if opts.do_gvcf {
        Defaults::Gvcf
    } else {
        Defaults::Vcf
    };
    let (rules, keep_ac_an) =
        rules::init_info_rules(args.info_rules.as_deref(), &mut out_hdr, defaults)?;
    opts.keep_ac_an = keep_ac_an;
    let missing_rules =
        rules::init_missing_rules(args.missing_rules.as_deref(), &out_hdr, defaults)?;

    if writer.write_header(&out_hdr).is_err() {
        writer.abandon();
        bail!("[merge_vcf] Error: cannot write to {output_name}");
    }
    if args.print_header {
        return match writer.finish() {
            Ok(()) => Ok(ExitCode::SUCCESS),
            Err(e) => Err(output::write_failure(e, || {
                anyhow::anyhow!("[merge_vcf] Error: close failed .. {output_name}")
            })),
        };
    }
    if let Some(fmt) = write_index {
        let ok = match output::index_request(fmt, out_type, args.output.as_deref()) {
            Some(req) => writer.request_index(req).is_ok(),
            None => false,
        };
        if !ok {
            writer.abandon();
            bail!("Error: failed to initialise index for {output_name}");
        }
        writer.start_index(&out_hdr)?;
    }

    let n = sr.reader_count();
    let ma = Maux::new(&sr);
    if ma.nout_smpl != out_hdr.sample_count() {
        // bcftools fails an assertion here.
        writer.abandon();
        bail!(
            "Error: the header has {} samples, the files {}",
            out_hdr.sample_count(),
            ma.nout_smpl
        );
    }
    // Worker threads format VCF lines; BCF records are already encoded by then.
    let emitter = (threads > 0 && writer.output_type().is_vcf()).then(|| {
        Emitter::new(&out_hdr, &writer, threads, output_name.clone())
            .with_write_error(format!("[merge_line] Error: cannot write to {output_name}"))
    });
    let out_hdr_ids = (out_hdr.id(b"END"), out_hdr.id(b"AN"), out_hdr.id(b"AC"));
    let mut merger = Merger {
        inp: Inputs {
            sr,
            gvcf_lines: (0..n).map(|_| Record::new()).collect(),
            gvcf_ref: vec![None; n],
        },
        out_hdr,
        out: Record::new(),
        ma,
        opts,
        rules,
        missing_rules,
        regs,
        fai,
        emitter,
        writer,
        output_name,
        vcmp: Vcmp::default(),
        ids: IdMaps {
            maps: vec![Vec::new(); n],
            end: out_hdr_ids.0,
            an: out_hdr_ids.1,
            ac: out_hdr_ids.2,
        },
        seen: Seen::default(),
        fmt: format::FormatState::default(),
        site: site::SiteState::default(),
        scratch_i: Vec::new(),
        scratch_f: Vec::new(),
        trim_mask: Vec::new(),
    };
    if merger.opts.do_gvcf {
        merger.enable_gvcf();
    }
    merger.out.clear();
    merger.run_and_finish()
}

/// Opens the output (bcftools' `hts_open` of `-o`).
fn create_writer(
    output: Option<&Path>,
    out_type: OutputType,
    level: Option<u8>,
    threads: usize,
    output_name: &str,
) -> anyhow::Result<Writer> {
    let out_path = output.map(output::data_path);
    match Writer::create(
        out_path.as_deref().filter(|p| p.as_os_str() != "-"),
        &output::write_options(out_type, level, threads),
    ) {
        Ok(w) => Ok(w),
        Err(e) => bail!("Can't write to \"{output_name}\": {}", strerror(&e)),
    }
}

/// The output header: the `--use-header` file, or the inputs' headers merged. The writer is
/// abandoned where bcftools exits before closing it.
fn output_header(
    args: &Args,
    sr: &SyncedReaders,
    opts: &Opts,
    writer: Writer,
) -> anyhow::Result<(Header, Writer)> {
    let mut hdr = Header::new_for_writing();
    hdr.set_version(&sr.header(0).version_string());
    if let Some(f) = &args.use_header {
        if header::read_header_file(&mut hdr, f).is_err() {
            writer.abandon();
            bail!("Could not read/parse the header: {f}");
        }
        return Ok((hdr, writer));
    }
    let duplicates = if args.force_samples {
        DuplicateSamples::Rename
    } else {
        DuplicateSamples::Fail
    };
    for i in 0..sr.reader_count() {
        if let Err(e) =
            header::merge_headers(&mut hdr, sr.header(i), &(i + 1).to_string(), duplicates)
        {
            writer.abandon();
            return Err(e);
        }
    }
    if opts.local_alleles > 0 {
        header::add_localized_tags(&mut hdr)?;
    }
    if !args.no_version {
        version::append_version(&mut hdr, "merge")?;
    }
    Ok((hdr, writer))
}

impl Merger {
    /// Merges the inputs and closes the output, abandoning it on an error (the end of
    /// bcftools' `merge_vcf`).
    fn run_and_finish(mut self) -> anyhow::Result<ExitCode> {
        let result = self.run().and_then(|()| {
            let Self {
                emitter,
                out_hdr,
                writer,
                ..
            } = &mut self;
            emitter
                .as_mut()
                .map_or(Ok(()), |e| e.finish(out_hdr, writer))
        });
        if result.is_err() {
            // Write what a serial run would have written before the error.
            let Self {
                emitter,
                out_hdr,
                writer,
                ..
            } = &mut self;
            if let Some(e) = emitter {
                e.abort(out_hdr, writer);
            }
        }
        let Self {
            writer,
            output_name,
            ..
        } = self;
        if let Err(e) = result {
            writer.abandon();
            if e.downcast_ref::<ExitOne>().is_some() {
                return Ok(ExitCode::from(1));
            }
            return Err(e);
        }
        match writer.finish() {
            Ok(()) => Ok(ExitCode::SUCCESS),
            Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
            Err(e @ WriteError::Deferred(_)) => Err(output::write_failure(e, || {
                anyhow::anyhow!("[merge_line] Error: cannot write to {output_name}")
            })),
            Err(e) => Err(output::write_failure(e, || {
                anyhow::anyhow!("[merge_vcf] Error: close failed .. {output_name}")
            })),
        }
    }

    /// The main loop of bcftools' `merge_vcf`.
    fn run(&mut self) -> anyhow::Result<()> {
        let mut rid_tab = vec![0i32; self.ma.n];
        loop {
            match self.inp.sr.next_line() {
                Ok(0) => break,
                Ok(_) => {}
                Err(e) => bail!("Error: {e}"),
            }
            // Output the cached gVCF blocks that end before the new record.
            if self.opts.do_gvcf {
                self.gvcf_flush(Flush::BeforeNextRecord)?;
            }
            self.ma.reset(&self.inp.sr, &mut rid_tab);
            self.unpack_position();
            if self.opts.do_gvcf {
                self.gvcf_stage(self.ma.pos)?;
            }
            while self.can_merge()? {
                self.stage_line();
                self.merge_line()?;
            }
            self.clean_buffer();
        }
        if let Some(e) = self.inp.sr.take_pending_error() {
            bail!("Error: {e}");
        }
        if self.opts.do_gvcf {
            self.gvcf_flush(Flush::All)?;
        }
        Ok(())
    }

    /// Decodes the records at this position.
    fn unpack_position(&mut self) {
        for i in 0..self.ma.n {
            let (beg, end, src) = {
                let b = &self.ma.buf[i];
                (b.beg, b.end, b.src)
            };
            for j in beg..end {
                let line = self.inp.line_mut(src, i, j);
                line.unpack(Unpack::ALL);
                line.variant_types_all();
            }
        }
    }

    /// Merges the staged records into one output line and writes it (bcftools' `merge_line`).
    fn merge_line(&mut self) -> anyhow::Result<()> {
        self.merge_chrom2qual()?;
        if let Some(regs) = &mut self.regs {
            let chr = self.ma.chr.as_deref().unwrap_or_default();
            let beg = self.out.pos();
            let end = beg + self.out.rlen() - 1;
            if !regs.overlaps(chr, to_u32(beg), to_u32(end)) {
                self.out.clear();
                return Ok(());
            }
        }
        self.merge_filter()?;
        self.merge_info()?;
        if self.opts.do_gvcf
            && let Some(id) = self.ids.end
        {
            let _ = self.out.remove_info_id(&self.out_hdr, id);
        }
        self.merge_format()?;
        self.trim_star_allele()?;
        self.write_out()
    }

    /// Trims an unobserved allele, at variant sites or everywhere (`-m *` or `**`).
    fn trim_star_allele(&mut self) -> anyhow::Result<()> {
        let trim = self.opts.trim_star_allele;
        if trim == 0 || (self.out.allele_count() <= 2 && trim <= 1) {
            return Ok(());
        }
        let iunseen = alleles::unseen_allele(self.out.alleles());
        if iunseen == 0 {
            return Ok(());
        }
        let rm = &mut self.trim_mask;
        rm.clear();
        rm.extend((0..self.out.allele_count()).map(|i| i == iunseen));
        if pasta::variant::utils::remove_allele_set(&self.out_hdr, &mut self.out, rm).is_err() {
            bail!(
                "[merge_line] Error: failed to trim the unobserved allele at {}:{}",
                String::from_utf8_lossy(self.out.chrom_safe(&self.out_hdr)),
                self.out.pos() + 1
            );
        }
        Ok(())
    }

    /// Writes the output line and clears it.
    fn write_out(&mut self) -> anyhow::Result<()> {
        // Encode the line built field by field, so that it is written from its blocks.
        self.out.sync();
        match &mut self.emitter {
            None => {
                if let Err(e) = self.writer.write_record(&self.out_hdr, &mut self.out) {
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[merge_line] Error: cannot write to {}", self.output_name)
                    }));
                }
            }
            Some(emitter) => {
                let next = emitter.spare().unwrap_or_default();
                let rec = std::mem::replace(&mut self.out, next);
                if let Some(rec) = emitter.write(&self.out_hdr, &mut self.writer, rec)? {
                    self.out = rec;
                }
            }
        }
        self.out.clear();
        Ok(())
    }
}

/// A 0-based coordinate for the region index, clamped to its range.
#[inline]
fn to_u32(x: i64) -> u32 {
    x.clamp(0, i64::from(regidx::REGIDX_MAX)) as u32
}
