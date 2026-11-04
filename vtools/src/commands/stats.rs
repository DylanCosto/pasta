//! `vtools stats`: statistics of a VCF/BCF file for plot-vcfstats (bcftools stats).
//!
//! With one input, records may be collected on worker threads ([`parallel`]). With two, the
//! synced reader pairs their records ([`pair`]): records of only one input and records in both
//! are counted in separate sets, and with `-s/-S` the genotypes of pairs are compared
//! ([`compare`]).

mod bins;
mod collect;
mod compare;
mod indel_context;
mod pair;
mod parallel;
mod report;

use std::{
    collections::HashSet,
    io::{self, IsTerminal, Write},
    process::ExitCode,
};

use anyhow::bail;
use clap::ArgMatches;
use pasta::variant::{
    Header, HeaderRecordType, Record, Unpack,
    io::ReadError,
    log,
    number::strtod,
    readlist::{self, ListSource},
    synced::{Overlap, Regions, SyncedError, SyncedOptions, SyncedReader},
    utils::AcSource,
};

use self::{
    bins::{AfBins, DepthBins},
    collect::{Collector, Config, Fatal, UserStatDef},
    indel_context::IndelContext,
    report::ReportInfo,
};
use crate::options::later;
use crate::{
    commands::head::apply_verbosity,
    filter::{self, FilterLogic, RecordFilter},
};

const USAGE: &str = "
About:   Parses VCF or BCF and produces stats which can be plotted using plot-vcfstats.
         When two files are given, the program generates separate stats for intersection
         and the complements. By default only sites are compared, -s/-S must given to include
         also sample columns.
Usage:   vtools stats [options] <A.vcf.gz> [<B.vcf.gz>]

Options:
        --af-bins LIST               Allele frequency bins, a list (0.1,0.5,1) or a file (0.1\\n0.5\\n1)
        --af-tag STRING              Allele frequency tag to use, by default estimated from AN,AC or GT
    -1, --1st-allele-only            Include only 1st allele at multiallelic sites
    -c, --collapse STRING            Treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]
    -d, --depth INT,INT,INT          Depth distribution: min,max,bin size [0,500,1]
    -e, --exclude EXPR               Exclude sites for which the expression is true (see man page for details)
    -E, --exons FILE.gz              Tab-delimited file with exons for indel frameshifts (chr,beg,end; 1-based, inclusive, bgzip compressed)
    -f, --apply-filters LIST         Require at least one of the listed FILTER strings (e.g. \"PASS,.\")
    -F, --fasta-ref FILE             Faidx indexed reference sequence file to determine INDEL context
    -i, --include EXPR               Select sites for which the expression is true (see man page for details)
    -I, --split-by-ID                Collect stats for sites with ID separately (known vs novel)
    -r, --regions REGION             Restrict to comma-separated list of regions
    -R, --regions-file FILE          Restrict to regions listed in a file
        --regions-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --samples LIST               List of samples for sample stats, \"-\" to include all samples
    -S, --samples-file FILE          File of samples to include
    -t, --targets REGION             Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE          Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
    -u, --user-tstv TAG[:min:max:n]  Collect Ts/Tv stats for any tag using the given binning [0:1:100]
                                       A subfield can be selected as e.g. 'PV4[0]', here the first value of the PV4 tag
        --threads INT                Use multithreading with <int> worker threads [0]
    -v, --verbosity INT              Verbosity level

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin), and optionally a second file to compare.
    inputs: Vec<String>,

    /// Allele frequency bins, a list (0.1,0.5,1) or a file (0.1\n0.5\n1).
    #[arg(long = "af-bins", value_name = "LIST", allow_hyphen_values = true)]
    af_bins: Option<String>,

    /// Allele frequency tag to use, by default estimated from AN,AC or GT.
    #[arg(long = "af-tag", value_name = "STRING", allow_hyphen_values = true)]
    af_tag: Option<String>,

    /// Include only 1st allele at multiallelic sites.
    #[arg(short = '1', long = "1st-allele-only")]
    first_allele_only: bool,

    /// Treat as identical records with <snps|indels|both|all|some|none> [none].
    #[arg(short = 'c', long = "collapse", value_name = "STRING", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    collapse: Vec<String>,

    /// Depth distribution: min,max,bin size [0,500,1].
    #[arg(short = 'd', long = "depth", value_name = "INT,INT,INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    depth: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Tab-delimited file with exons for indel frameshifts (chr,beg,end; 1-based, inclusive).
    #[arg(
        short = 'E',
        long = "exons",
        value_name = "FILE.gz",
        allow_hyphen_values = true
    )]
    exons: Option<String>,

    /// Require at least one of the listed FILTER strings (e.g. "PASS,.").
    #[arg(
        short = 'f',
        long = "apply-filters",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    apply_filters: Option<String>,

    /// Faidx indexed reference sequence file to determine INDEL context.
    #[arg(
        short = 'F',
        long = "fasta-ref",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    fasta_ref: Option<String>,

    /// Select sites for which the expression is true.
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

    /// Collect stats for sites with ID separately (known vs novel).
    #[arg(short = 'I', long = "split-by-ID")]
    split_by_id: bool,

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
    #[arg(long = "regions-overlap", value_name = "0|1|2", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions_overlap: Vec<String>,

    /// List of samples for sample stats, "-" to include all samples.
    #[arg(
        short = 's',
        long = "samples",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to include.
    #[arg(
        short = 'S',
        long = "samples-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    samples_file: Option<String>,

    /// Similar to -r but streams rather than index-jumps.
    #[arg(
        short = 't',
        long = "targets",
        value_name = "REGION",
        allow_hyphen_values = true
    )]
    targets: Option<String>,

    /// Similar to -R but streams rather than index-jumps.
    #[arg(
        short = 'T',
        long = "targets-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    targets_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0].
    #[arg(long = "targets-overlap", value_name = "0|1|2", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    targets_overlap: Vec<String>,

    /// Collect Ts/Tv stats for any tag using the given binning [0:1:100].
    #[arg(short = 'u', long = "user-tstv", value_name = "TAG[:min:max:n]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    user_tstv: Vec<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Verbosity level (`-v` alone: verbose per-site output with two files).
    #[arg(
        short = 'v',
        long = "verbosity",
        visible_alias = "verbose",
        value_name = "INT",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = NO_VALUE,
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,

    /// Print the usage (to standard error, exiting with 1, as bcftools does).
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Count)]
    help: u8,
}

/// The value clap gives `-v` without a value: a NUL, which no command-line word can hold, so that
/// `--verbosity=` (an empty value, level 0) stays apart from `-v` (verbose output).
const NO_VALUE: &str = "\0";

/// The size of the report's output buffer.
const OUT_BUFFER: usize = 1 << 16;

/// An option that bcftools' option loop validates when it reads it.
enum OptionValue<'a> {
    Collapse(&'a str),
    Verbosity(&'a str),
    Depth(&'a str),
    User(&'a str),
    Filter(&'a str, FilterLogic),
    RegionsOverlap(&'a str),
    TargetsOverlap(&'a str),
    Help,
}

/// The occurrences of the options bcftools validates as it reads them, in command-line order.
fn ordered_options<'a>(args: &'a Args, matches: &ArgMatches) -> Vec<(usize, OptionValue<'a>)> {
    let mut events: Vec<(usize, OptionValue<'a>)> = Vec::new();
    let mut add = |id: &str, values: &'a [String], make: fn(&'a str) -> OptionValue<'a>| {
        if let Some(idx) = matches.indices_of(id) {
            events.extend(idx.zip(values).map(|(i, v)| (i, make(v))));
        }
    };
    add("collapse", &args.collapse, OptionValue::Collapse);
    add("verbosity", &args.verbosity, OptionValue::Verbosity);
    add("depth", &args.depth, OptionValue::Depth);
    add("user_tstv", &args.user_tstv, OptionValue::User);
    add("include", &args.include, |v| {
        OptionValue::Filter(v, FilterLogic::Include)
    });
    add("exclude", &args.exclude, |v| {
        OptionValue::Filter(v, FilterLogic::Exclude)
    });
    add(
        "regions_overlap",
        &args.regions_overlap,
        OptionValue::RegionsOverlap,
    );
    add(
        "targets_overlap",
        &args.targets_overlap,
        OptionValue::TargetsOverlap,
    );
    // A count has a default value, with an index of its own: only real occurrences count.
    if args.help > 0
        && let Some(idx) = matches.indices_of("help")
    {
        events.extend(idx.map(|i| (i, OptionValue::Help)));
    }
    events.sort_unstable_by_key(|e| e.0);
    events
}

/// The options of bcftools' option loop, read in command-line order.
struct LoopOptions {
    verbose_sites: bool,
    depth: DepthBins,
    filter: Option<(String, FilterLogic)>,
    user: Vec<UserStatDef>,
    regions_overlap: Overlap,
    targets_overlap: Overlap,
    /// The other options are not read: `-h` prints the usage.
    help: bool,
}

impl LoopOptions {
    /// Reads the options in order, failing at the first invalid one as bcftools does.
    fn parse(args: &Args, matches: &ArgMatches) -> anyhow::Result<Self> {
        let mut o = Self {
            verbose_sites: false,
            depth: DepthBins::default(),
            filter: None,
            user: Vec::new(),
            regions_overlap: Overlap::Record,
            targets_overlap: Overlap::Pos,
            help: false,
        };
        for (_, opt) in ordered_options(args, matches) {
            match opt {
                OptionValue::Collapse(c) => {
                    if !matches!(
                        c,
                        "snps" | "indels" | "both" | "any" | "all" | "some" | "none"
                    ) {
                        bail!("The --collapse string \"{c}\" not recognised.");
                    }
                }
                OptionValue::Verbosity(NO_VALUE) => o.verbose_sites = true,
                OptionValue::Verbosity(v) => {
                    apply_verbosity(v)?;
                    // bcftools keeps the level in an int.
                    let level = pasta::variant::number::strtol(v.as_bytes(), 0).0 as i32;
                    o.verbose_sites = level != 0;
                }
                OptionValue::Depth(d) => o.depth = parse_depth(d)?,
                OptionValue::User(u) => o.user.push(parse_user_stat(u)?),
                OptionValue::Filter(expr, logic) => {
                    if o.filter.is_some() {
                        bail!(
                            "Error: only one -i or -e expression can be given, and they cannot be combined"
                        );
                    }
                    o.filter = Some((expr.to_owned(), logic));
                }
                OptionValue::RegionsOverlap(v) => {
                    o.regions_overlap = parse_overlap(v, "regions-overlap")?;
                }
                OptionValue::TargetsOverlap(v) => {
                    o.targets_overlap = parse_overlap(v, "targets-overlap")?;
                }
                OptionValue::Help => {
                    o.help = true;
                    break;
                }
            }
        }
        Ok(o)
    }
}

/// Parses `-d min,max,step` as `sscanf("%d,%d,%d")` does.
fn parse_depth(s: &str) -> anyhow::Result<DepthBins> {
    let b = s.as_bytes();
    let mut pos = 0;
    let mut vals = [0i32; 3];
    for (k, v) in vals.iter_mut().enumerate() {
        if k > 0 {
            if b.get(pos) != Some(&b',') {
                bail!("Could not parse --depth {s}");
            }
            pos += 1;
        }
        while b.get(pos).is_some_and(u8::is_ascii_whitespace) {
            pos += 1;
        }
        let (n, end) = pasta::variant::number::strtol(b, pos);
        if end == pos {
            bail!("Could not parse --depth {s}");
        }
        *v = n as i32;
        pos = end;
    }
    let [min, max, step] = vals;
    // bcftools divides by a zero step; it is rejected here.
    if min < 0 || min >= max || step > max - min + 1 || step <= 0 {
        bail!("Is this a typo? --depth {s}");
    }
    // bcftools' bin count overflows an int (and it crashes) for huge ranges.
    if (i64::from(max) - i64::from(min)) / i64::from(step) + 4 > i64::from(i32::MAX) {
        bail!("Is this a typo? --depth {s}");
    }
    Ok(DepthBins { min, max, step })
}

/// Parses a `-u TAG[idx][:min:max:nbins]` definition (`add_user_stats`). bcftools reads past
/// the end of the argument when it stops before `nbins`; here the fields given are used.
fn parse_user_stat(s: &str) -> anyhow::Result<UserStatDef> {
    let b = s.as_bytes();
    let colon = s.find(':').unwrap_or(s.len());
    let mut tag_end = colon;
    let mut idx = 0usize;
    if colon > 0
        && b[colon - 1] == b']'
        && let Some(open) = s[..colon].rfind('[')
    {
        let (v, end) = pasta::variant::number::strtol(b, open + 1);
        if end == open + 1 || end != colon - 1 {
            bail!("Could not parse the index in \"{s}\"");
        }
        if v < 0 {
            bail!("Error: negative index is not allowed: \"{s}\"");
        }
        idx = v as usize;
        tag_end = open;
    }
    // bcftools cuts the argument at '[' and names what is left in its messages.
    let shown = if tag_end < colon { &s[..tag_end] } else { s };
    let tag = &s[..tag_end];
    let (mut min, mut max, mut nbins) = (0f32, 1f32, 100usize);
    let mut t = colon;
    if t < b.len() {
        t += 1;
        let (v, end) = strtod(b, t);
        if end == t {
            bail!("Could not parse {shown}");
        }
        min = v as f32;
        t = end + 1;
    }
    if t < b.len() {
        let (v, end) = strtod(b, t);
        if end == t {
            bail!("Could not parse {shown}");
        }
        max = v as f32;
        t = end + 1;
    }
    if t < b.len() {
        let (v, end) = pasta::variant::number::strtol(b, t);
        if end == t {
            bail!("Could not parse {shown}");
        }
        if v <= 0 {
            bail!("Number of bins does not make sense ({v}): {shown}.");
        }
        nbins = v as usize;
    }
    // bcftools crashes on bins with an infinite or NaN bound.
    if !min.is_finite() || !max.is_finite() {
        bail!("Could not parse {shown}");
    }
    Ok(UserStatDef {
        tag: tag.as_bytes().to_vec(),
        idx,
        min,
        max,
        nbins,
    })
}

fn parse_overlap(s: &str, what: &str) -> anyhow::Result<Overlap> {
    Overlap::parse(s).ok_or_else(|| anyhow::anyhow!("Could not parse: --{what} {s}"))
}

/// Resolves `-s/-S` to header columns and names (HTSlib's `bcf_sr_set_samples`).
fn select_samples(
    header: &Header,
    fname: &str,
    spec: &str,
    source: ListSource,
) -> anyhow::Result<Vec<(usize, Vec<u8>)>> {
    let fail = || -> anyhow::Error {
        if header.sample_count() == 0 {
            anyhow::anyhow!("No sample columns in {fname}")
        } else {
            anyhow::anyhow!("Unable to parse the samples: \"{spec}\"")
        }
    };
    let all = || -> Vec<Vec<u8>> { header.samples().iter().map(|s| s.to_vec()).collect() };
    let names: Vec<Vec<u8>> = if spec == "-" {
        all()
    } else if let Some(body) = spec.strip_prefix('^') {
        // HTSlib keeps the '^' on the first name and then excludes every listed name, so
        // nothing is left; here the listed samples are excluded from the header's.
        let Ok(list) = readlist::read_list(body.as_bytes(), source) else {
            log::error(
                "bcf_sr_set_samples",
                format_args!("Could not read the file: \"{spec}\""),
            );
            return Err(fail());
        };
        let excluded: HashSet<Vec<u8>> = list.into_iter().collect();
        all()
            .into_iter()
            .filter(|s| !excluded.contains(s))
            .collect()
    } else {
        let Ok(list) = readlist::read_list(spec.as_bytes(), source) else {
            log::error(
                "bcf_sr_set_samples",
                format_args!("Could not read the file: \"{spec}\""),
            );
            return Err(fail());
        };
        list
    };
    let mut out = Vec::new();
    for name in names {
        match header.sample_id(&name) {
            Some(i) => out.push((i as usize, name)),
            None => log::warn(
                "bcf_sr_set_samples",
                format_args!(
                    "The sample \"{}\" was not found in {fname}, skipping",
                    String::from_utf8_lossy(&name)
                ),
            ),
        }
    }
    if out.is_empty() {
        return Err(fail());
    }
    Ok(out)
}

/// The settings of a run, besides the shared [`Config`].
struct Setup {
    cfg: Config,
    /// The `-u` tag types of each set's header.
    user_types: Vec<Vec<bool>>,
    filter: Option<(String, FilterLogic)>,
    fasta_ref: Option<String>,
    sample_names: Vec<Vec<u8>>,
    /// Whether VCF FORMAT columns are parsed (with samples, or a filter that reads them).
    parse_format: bool,
}

impl Setup {
    /// How much of each VCF line is parsed: FORMAT only when samples are analysed (or the
    /// filter reads it), as bcftools sets `max_unpack`.
    fn max_unpack(&self) -> Unpack {
        if self.parse_format {
            Unpack::empty()
        } else {
            Unpack::INFO
        }
    }
}

/// Prints the usage and fails, as bcftools' `usage()` does.
fn usage() -> ExitCode {
    eprint!("{USAGE}");
    ExitCode::FAILURE
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let LoopOptions {
        verbose_sites,
        depth,
        filter,
        user,
        regions_overlap,
        targets_overlap,
        help,
    } = LoopOptions::parse(&args, matches)?;
    if help {
        return Ok(usage());
    }
    // strtol: a negative count fails to create the thread pool.
    let threads = args
        .threads
        .as_deref()
        .map_or(0, |t| i64::from(crate::output::thread_count(t)));

    let mut files = args.inputs.clone();
    if files.is_empty() {
        if io::stdin().is_terminal() {
            return Ok(usage());
        }
        files.push("-".into());
    }
    if files.len() > 2 {
        return Ok(usage());
    }
    if files.len() == 2 && args.split_by_id {
        bail!("Only one file can be given with -i.");
    }

    let sample_spec = later(
        matches,
        ("samples", args.samples.as_deref()),
        ("samples_file", args.samples_file.as_deref()),
    );
    let options = synced_options(&args, matches, regions_overlap, targets_overlap, threads)?;
    let run = Run {
        files: &files,
        sample_spec,
        filter,
        depth,
        user,
        verbose_sites,
    };
    if files.len() == 2 {
        return pair::run(&args, run, options);
    }
    run_single(&args, run, options)
}

/// The reader's targets (`-t/-T`) and regions (`-r/-R`), the later of each pair winning.
fn synced_options(
    args: &Args,
    matches: &ArgMatches,
    regions_overlap: Overlap,
    targets_overlap: Overlap,
    threads: i64,
) -> anyhow::Result<SyncedOptions> {
    let target_spec = later(
        matches,
        ("targets", args.targets.as_deref()),
        ("targets_file", args.targets_file.as_deref()),
    );
    let region_spec = later(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    );
    let mut targets = None;
    let mut targets_exclude = false;
    if let Some((spec, source)) = target_spec {
        let body = match spec.strip_prefix('^') {
            Some(rest) => {
                targets_exclude = true;
                rest
            }
            None => spec,
        };
        let mut t = Regions::new(body.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the targets: {spec}"))?;
        t.overlap = targets_overlap;
        targets = Some(t);
    }
    let mut regions = None;
    if let Some((spec, source)) = region_spec {
        let mut r = Regions::new(spec.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = regions_overlap;
        regions = Some(r);
    }
    if threads < 0 || !crate::output::threads_start(threads as usize) {
        bail!("Failed to create threads");
    }
    Ok(SyncedOptions {
        regions,
        targets,
        targets_exclude,
        apply_filters: args.apply_filters.as_ref().map(|f| f.as_bytes().to_vec()),
        threads: threads as usize,
    })
}

/// Runs `stats A`: collects the records on worker threads with `--threads` (unless `-E` is
/// given), else on this one.
fn run_single(args: &Args, run: Run<'_>, options: SyncedOptions) -> anyhow::Result<ExitCode> {
    let files = run.files;
    let verbose_sites = run.verbose_sites;
    let threads = options.threads;
    let fname = files[0].as_str();
    let mut reader = match SyncedReader::open(fname, options) {
        Ok(r) => r,
        Err(e) => return Err(open_error(fname, &e)),
    };
    let setup = single_setup(args, run, fname, reader.header())?;

    let exons = open_exons(args)?;
    reader.set_max_unpack(setup.max_unpack());
    let indel_ctx = setup.fasta_ref.as_deref().and_then(IndelContext::open);
    let has_indel_ctx = indel_ctx.is_some();

    let nsamples = [reader.header().sample_count()];
    let mut out = io::BufWriter::with_capacity(OUT_BUFFER, crate::output::stdout());
    let info = ReportInfo {
        files,
        nsamples: &nsamples,
        sample_names: &setup.sample_names,
        split_by_id: args.split_by_id,
        verbose_sites,
        exons: args.exons.is_some(),
        indel_ctx: has_indel_ctx,
    };
    report::print_header(&mut out, &info)?;

    let parallel = if threads > 0 && exons.is_none() {
        parallel::collect(&mut reader, &setup, indel_ctx.as_ref(), threads)
    } else {
        Ok(None)
    };
    // Without threads, or when no worker thread could start or get its tables.
    let collected = match parallel {
        Ok(None) => collect_serial(&mut reader, &setup, indel_ctx, exons),
        Ok(Some(c)) => Ok(c),
        Err(f) => Err(f),
    };
    let mut collector = match collected {
        Ok(c) => c,
        Err(f) => return fatal_exit(&mut out, f),
    };
    if let Some(e) = reader.take_pending_error() {
        out.flush()?;
        bail!("Error: {e}");
    }

    report::print_stats(&mut out, &setup.cfg, &info, &mut collector.sets, None)?;
    out.flush()?;
    Ok(ExitCode::SUCCESS)
}

/// The setup for one input `fname` with this header.
fn single_setup(args: &Args, run: Run<'_>, fname: &str, header: &Header) -> anyhow::Result<Setup> {
    let bins = check_options(args, &[header], run.filter.as_ref())?;
    let (columns, names): (Vec<usize>, Vec<Vec<u8>>) = match run.sample_spec {
        Some((spec, source)) => select_samples(header, fname, spec, source)?
            .into_iter()
            .unzip(),
        None => (Vec::new(), Vec::new()),
    };
    let samples = Samples {
        columns: vec![columns],
        names,
        given: run.sample_spec.is_some(),
    };
    init_stats(args, &[header], samples, bins, run)
}

/// The parsed options that both modes share.
struct Run<'a> {
    files: &'a [String],
    sample_spec: Option<(&'a str, ListSource)>,
    filter: Option<(String, FilterLogic)>,
    depth: DepthBins,
    user: Vec<UserStatDef>,
    verbose_sites: bool,
}

/// The samples of `-s/-S`.
struct Samples {
    /// Their header columns in each input.
    columns: Vec<Vec<usize>>,
    names: Vec<Vec<u8>>,
    /// Whether `-s/-S` was given.
    given: bool,
}

/// The error of an input that cannot be opened (`bcf_sr_add_reader` failing; the synced
/// reader has logged HTSlib's message).
fn open_error(fname: &str, e: &SyncedError) -> anyhow::Error {
    let name = if fname == "-" {
        "standard input"
    } else {
        fname
    };
    anyhow::anyhow!("Failed to read from {name}: {e}")
}

/// Reads `-E/--exons`.
fn open_exons(args: &Args) -> anyhow::Result<Option<Regions>> {
    match &args.exons {
        Some(path) => match Regions::from_file(std::path::Path::new(path)) {
            Ok(r) => Ok(Some(r)),
            Err(_) => {
                bail!("Error occurred while reading, was the file compressed with bgzip: {path}?")
            }
        },
        None => Ok(None),
    }
}

/// Ends the program after a [`Fatal`] error, keeping what was printed so far.
fn fatal_exit(out: &mut impl Write, f: Fatal) -> anyhow::Result<ExitCode> {
    out.flush()?;
    if let Some(m) = f.message {
        eprint!("{m}");
    }
    Ok(ExitCode::from(f.code))
}

/// The first checks of `init_stats`: the filter compiles for every input, the frequency bins,
/// and `--af-tag` (in the first input). Returns the number of frequency bins and the bins of
/// `--af-bins`.
fn check_options(
    args: &Args,
    headers: &[&Header],
    filter: Option<&(String, FilterLogic)>,
) -> anyhow::Result<(usize, Option<AfBins>)> {
    if let Some((expr, _)) = filter {
        // Compile once here to report errors before anything is printed.
        for header in headers {
            filter::init(header, expr)?;
        }
    }
    let bins = match &args.af_bins {
        None => {
            let most = headers.iter().map(|h| h.sample_count() + 1).max();
            (most.unwrap_or(0).max(101), None)
        }
        Some(list) => {
            let b = AfBins::parse(list)?;
            (b.len(), Some(b))
        }
    };
    if let Some(tag) = &args.af_tag {
        let header = headers[0];
        let ok = header
            .id(tag.as_bytes())
            .is_some_and(|id| header.id_exists(HeaderRecordType::Info, id));
        if !ok {
            bail!("No such INFO tag: {tag}");
        }
    }
    Ok(bins)
}

/// The types of the `-u` tags in a header (`init_user_stats`): whether each is a Float.
fn user_types(header: &Header, user: &[UserStatDef]) -> anyhow::Result<Vec<bool>> {
    let mut types = Vec::with_capacity(user.len());
    for u in user {
        let tag = String::from_utf8_lossy(&u.tag);
        let Some(id) = header
            .id(&u.tag)
            .filter(|&id| header.id_exists(HeaderRecordType::Info, id))
        else {
            bail!("The INFO tag \"{tag}\" is not defined in the header");
        };
        match header.id_type_code(HeaderRecordType::Info, id) {
            1 => types.push(false),
            2 => types.push(true),
            t => bail!("The INFO tag \"{tag}\" is not of Float or Integer type ({t})"),
        }
    }
    Ok(types)
}

/// Sets up the statistics (`init_stats`) once [`check_options`] passed and the samples are
/// known: one set for a single input (or known and novel with `-I`), and for two inputs one
/// set for the records of each alone and one for records in both.
fn init_stats(
    args: &Args,
    headers: &[&Header],
    samples: Samples,
    (m_af, af_bins): (usize, Option<AfBins>),
    run: Run<'_>,
) -> anyhow::Result<Setup> {
    let header = headers[0];
    let has_fmt_ad = header
        .id(b"AD")
        .is_some_and(|id| header.id_exists(HeaderRecordType::Format, id));
    let nstats = match (args.split_by_id, headers.len()) {
        (true, _) => 2,
        (false, 1) => 1,
        (false, _) => 3,
    };
    // bcftools adds what the filter reads to `max_unpack` (`filter_max_unpack`), so sample
    // columns are parsed only when samples are analysed or the expression reads them.
    let filter_unpack = match &run.filter {
        Some((expr, _)) => filter::init(header, expr)?.max_unpack(),
        None => Unpack::empty(),
    };
    let parse_format = samples.given || filter_unpack.contains(Unpack::FMT);
    // Set 1 holds the second input's records; the others are read with the first's header.
    let mut types = Vec::with_capacity(nstats);
    for set in 0..nstats {
        let h = headers.get(1).filter(|_| set == 1).unwrap_or(&header);
        types.push(user_types(h, &run.user)?);
    }
    let cfg = Config {
        m_af,
        af_bins,
        af_tag: args.af_tag.as_ref().map(|t| t.as_bytes().to_vec()),
        first_allele_only: args.first_allele_only,
        split_by_id: args.split_by_id,
        samples: samples.columns,
        has_fmt_ad,
        ac_source: if samples.given {
            AcSource::InfoThenGenotypes
        } else {
            AcSource::Info
        },
        depth: run.depth,
        user: run.user,
        exons: args.exons.is_some(),
    };
    Ok(Setup {
        cfg,
        user_types: types,
        filter: run.filter,
        fasta_ref: args.fasta_ref.clone(),
        sample_names: samples.names,
        parse_format,
    })
}

/// A compiled `-i/-e` filter.
type Filter = (Box<dyn RecordFilter + Send>, FilterLogic);

fn compile_filter(setup: &Setup, header: &Header) -> Result<Option<Filter>, Fatal> {
    match &setup.filter {
        Some((expr, logic)) => match filter::init(header, expr) {
            Ok(f) => Ok(Some((f, *logic))),
            Err(e) => Err(Fatal::error(format!("{e:#}\n"))),
        },
        None => Ok(None),
    }
}

/// Applies the filter and collects one record.
#[inline]
fn process(
    collector: &mut Collector,
    cfg: &Config,
    filter: &mut Option<Filter>,
    header: &Header,
    rec: &mut Record,
) -> Result<(), Fatal> {
    if let Some((f, logic)) = filter {
        let pass = f
            .pass(header, rec)
            .map_err(|e| Fatal::error(format!("{e:#}\n")))?;
        if pass != (*logic == FilterLogic::Include) {
            return Ok(());
        }
    }
    collector.add(cfg, header, rec, 0, 0)
}

fn read_error(e: &SyncedError) -> Fatal {
    Fatal::error(format!("Error: {e}\n"))
}

/// Collects the statistics on this thread (`do_vcf_stats`).
fn collect_serial(
    reader: &mut SyncedReader,
    setup: &Setup,
    indel_ctx: Option<IndelContext>,
    exons: Option<Regions>,
) -> Result<Collector, Fatal> {
    let cfg = &setup.cfg;
    let mut collector = Collector::try_new(cfg, &setup.user_types, indel_ctx, exons)?;
    let mut filter = compile_filter(setup, reader.header())?;
    let mut rec = Record::new();
    if reader.is_plain_stream() {
        // Records are read straight into `rec`, without the synced reader's buffer (which
        // moves every record several times). The buffer only matters on a read error after
        // records at one position, which drops them; the report is not printed then anyway.
        rec.set_max_unpack(setup.max_unpack());
        let mut header = reader.header().clone();
        let mut line = Vec::new();
        loop {
            match reader
                .reader_mut()
                .read_record(&mut header, &mut line, &mut rec)
            {
                Ok(true) => {}
                Ok(false) => break,
                Err(ReadError::Parse(_)) => return Err(read_error(&SyncedError::VcfParse)),
                Err(_) => return Err(read_error(&SyncedError::BcfRead)),
            }
            process(&mut collector, cfg, &mut filter, &header, &mut rec)?;
        }
        return Ok(collector);
    }
    loop {
        match reader.next_record(&mut rec) {
            Ok(true) => {}
            Ok(false) => break,
            Err(e) => return Err(read_error(&e)),
        }
        process(&mut collector, cfg, &mut filter, reader.header(), &mut rec)?;
    }
    Ok(collector)
}
