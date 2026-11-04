//! `vtools roh`: an HMM for detecting runs of autozygosity (bcftools roh, `vcfroh.c`).
//!
//! Each target sample's sites are collected with their emission probabilities (from PL, or from
//! GT with `-G`) and the allele frequency of the site, then a two-state HMM (Hardy-Weinberg and
//! autozygous) is run over each chromosome: the Viterbi path gives the state of each site and
//! the forward-backward posteriors its quality. The arithmetic follows bcftools step for step,
//! so the printed numbers match.

mod genmap;
mod hmm;
mod output;
mod samples;
mod sites;

use std::{
    io::{self, IsTerminal, Write},
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, HeaderRecordType, Record, ValueType,
    io::InputFormat,
    number,
    readlist::ListSource,
    synced::{Overlap, SyncedReaders},
};

use self::{
    output::{Compression, OutFile},
    samples::{SampleList, Verbosity},
    sites::{AfSource, Config, EstimateFrom, Output, PlTypeError, Roh, Training},
};
use crate::options::later;
use crate::{
    commands::head::apply_verbosity,
    filter::{self, FilterLogic},
    raw_input::strerror,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin).
    inputs: Vec<String>,

    /// If AF is not known, use this allele frequency [skip].
    #[arg(long = "AF-dflt", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    af_dflt: Vec<String>,

    /// Use TAG for allele frequency.
    #[arg(long = "AF-tag", value_name = "TAG", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    af_tag: Vec<String>,

    /// Read allele frequencies from a tabix-indexed file (CHR\tPOS\tREF,ALT\tAF).
    #[arg(long = "AF-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    af_file: Vec<String>,

    /// Buffer size and the number of overlapping sites, 0 for unlimited [0].
    #[arg(
        short = 'b',
        long = "buffer-size",
        value_name = "INT[,INT]",
        allow_hyphen_values = true
    )]
    buffer_size: Option<String>,

    /// Estimate AF from FORMAT/TAG (GT or PL) of all samples ("-") or samples listed in FILE.
    #[arg(short = 'e', long = "estimate-AF", value_name = "[TAG],FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    estimate_af: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(long = "exclude", value_name = "EXPR", allow_hyphen_values = true)]
    exclude: Vec<String>,

    /// Use GTs and ignore PLs, using FLOAT for PL of the two least likely genotypes.
    #[arg(short = 'G', long = "GTs-only", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    gts_only: Vec<String>,

    /// Select sites for which the expression is true.
    #[arg(long = "include", value_name = "EXPR", allow_hyphen_values = true)]
    include: Vec<String>,

    /// Skip hom-ref genotypes (0/0).
    #[arg(short = 'i', long = "ignore-homref")]
    ignore_homref: bool,

    /// Include sites with no ALT allele (ignored by default).
    #[arg(long = "include-noalt")]
    include_noalt: bool,

    /// Skip indels as their genotypes are enriched for errors.
    #[arg(short = 'I', long = "skip-indels")]
    skip_indels: bool,

    /// Genetic map in IMPUTE2 format, single file or mask with "{CHROM}".
    #[arg(
        short = 'm',
        long = "genetic-map",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    genetic_map: Option<String>,

    /// Constant recombination rate per bp.
    #[arg(short = 'M', long = "rec-rate", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    rec_rate: Vec<String>,

    /// Write output to a file [standard output].
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<String>,

    /// Output s:per-site, r:regions, z:compressed [sr].
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "[srz]",
        allow_hyphen_values = true
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
    #[arg(long = "regions-overlap", value_name = "0|1|2", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions_overlap: Vec<String>,

    /// List of samples to analyze [all samples].
    #[arg(
        short = 's',
        long = "samples",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to analyze [all samples].
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

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Verbosity level.
    #[arg(short = 'v', long = "verbosity", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    verbosity: Vec<String>,

    /// P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state
    /// [6.7e-8].
    #[arg(short = 'a', long = "hw-to-az", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    hw_to_az: Vec<String>,

    /// P(HW|AZ) transition probability from AZ to HW state [5e-9].
    #[arg(short = 'H', long = "az-to-hw", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    az_to_hw: Vec<String>,

    /// Estimate HMM parameters, FLOAT is the convergence threshold, e.g. 1e-10 (experimental).
    #[arg(short = 'V', long = "viterbi-training", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    viterbi_training: Vec<String>,

    /// Print help.
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Count)]
    help: u8,

    /// Print help (getopt's `-?`).
    #[arg(short = '?', hide = true, action = clap::ArgAction::Count)]
    help_q: u8,
}

const USAGE: &str = "
About:   HMM model for detecting runs of autozygosity.
Usage:   vtools roh [options] <in.vcf.gz>

General Options:
        --AF-dflt FLOAT              If AF is not known, use this allele frequency [skip]
        --AF-tag TAG                 Use TAG for allele frequency
        --AF-file FILE               read allele frequencies from file (CHR\\tPOS\\tREF,ALT\\tAF)
    -b  --buffer-size INT[,INT]      Buffer size and the number of overlapping sites, 0 for unlimited [0]
                                       If the first number is negative, it is interpreted as the maximum memory to
                                       use, in MB. The default overlap is set to roughly 1% of the buffer size.
    -e, --estimate-AF [TAG],FILE     Estimate AF from FORMAT/TAG (GT or PL) of all samples (\"-\") or samples listed
                                       in FILE. If TAG is not given, the frequency is estimated from GT by default
        --exclude EXPR               Exclude sites for which the expression is true
    -G, --GTs-only FLOAT             Use GTs and ignore PLs, instead using FLOAT for PL of the two least likely genotypes.
                                       Safe value to use is 30 to account for GT errors.
        --include EXPR               Select sites for which the expression is true
    -i, --ignore-homref              Skip hom-ref genotypes (0/0)
        --include-noalt              Include sites with no ALT allele (ignored by default)
    -I, --skip-indels                Skip indels as their genotypes are enriched for errors
    -m, --genetic-map FILE           Genetic map in IMPUTE2 format, single file or mask, where string \"{CHROM}\"
                                       is replaced with chromosome name
    -M, --rec-rate FLOAT             Constant recombination rate per bp
    -o, --output FILE                Write output to a file [standard output]
    -O, --output-type [srz]          Output s:per-site, r:regions, z:compressed [sr]
    -r, --regions REGION             Restrict to comma-separated list of regions
    -R, --regions-file FILE          Restrict to regions listed in a file
        --regions-overlap 0|1|2      include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --samples LIST               List of samples to analyze [all samples]
    -S, --samples-file FILE          File of samples to analyze [all samples]
    -t, --targets REGION             Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE          Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT                Use multithreading with <int> worker threads [0]
    -v, --verbosity INT              Verbosity level

HMM Options:
    -a, --hw-to-az FLOAT             P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state [6.7e-8]
    -H, --az-to-hw FLOAT             P(HW|AZ) transition probability from AZ to HW state [5e-9]
    -V, --viterbi-training FLOAT     Estimate HMM parameters, FLOAT is the convergence threshold, e.g. 1e-10 (experimental)

Example:
   # Find RoH regions assuming default allele frequency 0.4
   vtools roh -G30 --AF-dflt 0.4 test.vcf -o out.txt
";

/// Prints the usage and returns bcftools' exit status for it.
fn usage() -> ExitCode {
    eprintln!("{USAGE}");
    ExitCode::from(1)
}

/// C's `strtod` over the whole option value, or bcftools' `Could not parse` error.
fn parse_float(value: &str, option: &str) -> anyhow::Result<f64> {
    let (v, end) = number::strtod(value.as_bytes(), 0);
    if end != value.len() {
        bail!("Could not parse: {option} {value}");
    }
    Ok(v)
}

fn parse_overlap(s: &str, what: &str) -> anyhow::Result<Overlap> {
    Overlap::parse(s).with_context(|| format!("Could not parse: --{what} {s}"))
}

/// The output bits of `-O` (`OUTPUT_ST`, `OUTPUT_RG`, `OUTPUT_GZ`).
fn output_kind(types: &[String]) -> (Output, Compression) {
    let (mut sites, mut regions, mut gz) = (false, false, false);
    for t in types {
        sites |= t.contains(['s', 'S']);
        regions |= t.contains(['r', 'R']);
        gz |= t.contains('z');
    }
    if !sites && !regions {
        (sites, regions) = (true, true);
    }
    let compression = if gz {
        Compression::Gzip
    } else {
        Compression::Plain
    };
    (Output { sites, regions }, compression)
}

/// An option that bcftools' option loop validates (or counts) when it reads it.
enum OptionValue<'a> {
    Verbosity(&'a str),
    /// `--AF-tag`, `--AF-file`, or `-e`: each occurrence counts towards their exclusivity.
    AfOption,
    AfDflt(&'a str),
    Filter(&'a str, FilterLogic),
    GtsOnly(&'a str),
    RecRate(&'a str),
    HwToAz(&'a str),
    AzToHw(&'a str),
    TargetsOverlap(&'a str),
    RegionsOverlap(&'a str),
    ViterbiTraining(&'a str),
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
    add("verbosity", &args.verbosity, OptionValue::Verbosity);
    add("af_tag", &args.af_tag, |_| OptionValue::AfOption);
    add("af_file", &args.af_file, |_| OptionValue::AfOption);
    add("estimate_af", &args.estimate_af, |_| OptionValue::AfOption);
    add("af_dflt", &args.af_dflt, OptionValue::AfDflt);
    add("include", &args.include, |v| {
        OptionValue::Filter(v, FilterLogic::Include)
    });
    add("exclude", &args.exclude, |v| {
        OptionValue::Filter(v, FilterLogic::Exclude)
    });
    add("gts_only", &args.gts_only, OptionValue::GtsOnly);
    add("rec_rate", &args.rec_rate, OptionValue::RecRate);
    add("hw_to_az", &args.hw_to_az, OptionValue::HwToAz);
    add("az_to_hw", &args.az_to_hw, OptionValue::AzToHw);
    add(
        "targets_overlap",
        &args.targets_overlap,
        OptionValue::TargetsOverlap,
    );
    add(
        "regions_overlap",
        &args.regions_overlap,
        OptionValue::RegionsOverlap,
    );
    add(
        "viterbi_training",
        &args.viterbi_training,
        OptionValue::ViterbiTraining,
    );
    // A count has a default value, with an index of its own: only real occurrences count.
    for (id, n) in [("help", args.help), ("help_q", args.help_q)] {
        if n > 0
            && let Some(idx) = matches.indices_of(id)
        {
            events.extend(idx.map(|i| (i, OptionValue::Help)));
        }
    }
    events.sort_unstable_by_key(|e| e.0);
    events
}

/// The options of bcftools' option loop that are validated as they are read.
struct LoopOptions {
    af_dflt: f64,
    unseen_pl: Option<f64>,
    rec_rate: f64,
    t2az: f64,
    t2hw: f64,
    training: Option<f64>,
    regions_overlap: Overlap,
    targets_overlap: Overlap,
    filter: Option<(String, FilterLogic)>,
    /// How many times `--AF-tag`, `--AF-file`, and `-e` were given.
    naf_opts: usize,
    /// `-h` was reached: the usage is printed.
    help: bool,
}

impl LoopOptions {
    /// Reads the options in command-line order, failing at the first invalid one as bcftools
    /// does.
    fn parse(args: &Args, matches: &ArgMatches) -> anyhow::Result<Self> {
        let mut n = Self {
            af_dflt: 0.0,
            unseen_pl: None,
            rec_rate: 0.0,
            t2az: 6.7e-8,
            t2hw: 5e-9,
            training: None,
            regions_overlap: Overlap::Record,
            targets_overlap: Overlap::Pos,
            filter: None,
            naf_opts: 0,
            help: false,
        };
        for (_, opt) in ordered_options(args, matches) {
            match opt {
                OptionValue::Verbosity(v) => apply_verbosity(v)?,
                OptionValue::AfOption => n.naf_opts += 1,
                OptionValue::AfDflt(v) => n.af_dflt = parse_float(v, "--AF-dflt")?,
                OptionValue::Filter(expr, logic) => {
                    if n.filter.is_some() {
                        bail!(
                            "Error: only one --include or --exclude expression can be given, and they cannot be combined"
                        );
                    }
                    n.filter = Some((expr.to_owned(), logic));
                }
                OptionValue::GtsOnly(v) => {
                    n.unseen_pl = Some(10f64.powf(-parse_float(v, "-G")? / 10.0));
                }
                OptionValue::RecRate(v) => n.rec_rate = parse_float(v, "-M")?,
                OptionValue::HwToAz(v) => n.t2az = parse_float(v, "-a")?,
                OptionValue::AzToHw(v) => n.t2hw = parse_float(v, "-H")?,
                OptionValue::TargetsOverlap(v) => {
                    n.targets_overlap = parse_overlap(v, "targets-overlap")?;
                }
                OptionValue::RegionsOverlap(v) => {
                    n.regions_overlap = parse_overlap(v, "regions-overlap")?;
                }
                OptionValue::ViterbiTraining(v) => {
                    n.training = Some(parse_float(v, "--viterbi-training")?);
                }
                OptionValue::Help => {
                    n.help = true;
                    break;
                }
            }
        }
        Ok(n)
    }
}

/// The value of the last occurrence of an option.
fn last(values: &[String]) -> Option<&str> {
    values.last().map(String::as_str)
}

/// Parses `-b`: the buffer size (sites, or with a negative number, megabytes) and overlap.
fn parse_buffer_size(spec: &str, nsamples: usize) -> anyhow::Result<(usize, usize)> {
    let bad = || anyhow::anyhow!("Could not parse: --buffer-size {spec}");
    let s = spec.as_bytes();
    let (size, end) = number::strtod(s, 0);
    let mut olap: Option<i64> = None;
    if end < s.len() {
        if s[end] != b',' {
            return Err(bad());
        }
        let (v, e) = number::strtol(s, end + 1);
        if e < s.len() || v < 0 {
            return Err(bad());
        }
        olap = Some(v);
    }
    let nbuf_max = if size < 0.0 {
        (size.abs() * 1e6 / f64::from(4 + 8 * 2) / nsamples as f64) as i32
    } else {
        size as i32
    };
    let nbuf_olap = olap.unwrap_or((f64::from(nbuf_max) * 0.01) as i64);
    Ok((nbuf_max.max(0) as usize, nbuf_olap as usize))
}

/// The header of the output (`init_data`).
fn header_text(output: Output, training: Option<Training>) -> String {
    use std::fmt::Write as _;

    let mut s = format!(
        "# This file was produced by: vtools roh({}+pasta-{})\n# The command line was:\tvtools roh",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION")
    );
    for a in std::env::args().skip(2) {
        s.push(' ');
        s.push_str(&a);
    }
    s.push_str("\n#\n");
    let columns = |s: &mut String, prefix: &str, names: &[&str]| {
        s.push_str(prefix);
        for (i, n) in names.iter().enumerate() {
            // Writing to a String cannot fail.
            let _ = write!(s, "\t[{}]{n}", i + 2);
        }
        s.push('\n');
    };
    if output.regions {
        columns(
            &mut s,
            "# RG",
            &[
                "Sample",
                "Chromosome",
                "Start",
                "End",
                "Length (bp)",
                "Number of markers",
                "Quality (average fwd-bwd phred score)",
            ],
        );
    }
    if output.sites {
        columns(
            &mut s,
            "# ST",
            &[
                "Sample",
                "Chromosome",
                "Position",
                "State (0:HW, 1:AZ)",
                "Quality (fwd-bwd phred score)",
            ],
        );
    }
    if training.is_some() {
        columns(
            &mut s,
            "# VT, Viterbi Training",
            &[
                "Sample",
                "Iteration",
                "dAZ",
                "dHW",
                "1 - P(HW|HW)",
                "P(AZ|HW)",
                "1 - P(AZ|AZ)",
                "P(HW|AZ)",
            ],
        );
    }
    s
}

fn display_name(fname: &str) -> &str {
    if fname == "-" {
        "standard input"
    } else {
        fname
    }
}

/// Sets up the reader with the regions, targets, and `--AF-file`, and opens the input.
fn open_input(
    args: &Args,
    matches: &ArgMatches,
    numbers: &LoopOptions,
    fname: &str,
    threads: i64,
) -> anyhow::Result<SyncedReaders> {
    let LoopOptions {
        regions_overlap,
        targets_overlap,
        ..
    } = *numbers;
    let targets = later(
        matches,
        ("targets", args.targets.as_deref()),
        ("targets_file", args.targets_file.as_deref()),
    );
    let af_file = last(&args.af_file);
    if af_file.is_some() && targets.is_some() {
        bail!("Error: The options --AF-file and -t are mutually exclusive");
    }

    let mut sr = SyncedReaders::new();
    if let Some((spec, source)) = later(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    ) {
        sr.set_regions_overlap(regions_overlap);
        sr.set_regions_spec(spec.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
    }
    if let Some((spec, source)) = targets {
        sr.set_targets_overlap(targets_overlap);
        sr.set_targets_spec(spec.as_bytes(), source, None)
            .map_err(|_| anyhow::anyhow!("Failed to read the targets: {spec}"))?;
    }
    if let Some(af_file) = af_file {
        // The REF,ALT alleles are in the third column.
        sr.set_targets_spec(af_file.as_bytes(), ListSource::File, Some(2))
            .map_err(|_| anyhow::anyhow!("Failed to read the targets: {af_file}"))?;
    }
    if threads < 0 || !crate::output::threads_start(threads as usize) {
        bail!("Failed to create threads");
    }
    sr.set_threads(threads as usize);
    sr.add_reader(fname)
        .map_err(|e| anyhow::anyhow!("Failed to read from {}: {e}", display_name(fname)))?;
    Ok(sr)
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let numbers = LoopOptions::parse(&args, matches)?;
    if numbers.help {
        return Ok(usage());
    }
    let filter_expr = numbers
        .filter
        .as_ref()
        .map(|(expr, logic)| (expr.as_str(), *logic));
    // strtol: a negative count fails to create the thread pool.
    let threads = args
        .threads
        .as_deref()
        .map_or(0, |t| i64::from(crate::output::thread_count(t)));
    let naf_opts = numbers.naf_opts;

    // bcftools writes standard output for "stdout" (its default name) and "-".
    let output_name = args.output.as_deref().unwrap_or("stdout");
    let (output, compression) = output_kind(&args.output_type);
    let fname = match args.inputs.first() {
        Some(f) => f.clone(),
        None if !io::stdin().is_terminal() => "-".to_owned(),
        None => return Ok(usage()),
    };

    if numbers.training.is_some() && args.buffer_size.is_some() {
        bail!("Error: cannot use -b with -V");
    }
    if !(0.0..=1.0).contains(&numbers.t2az) {
        bail!(
            "Error: The parameter --hw-to-az is not in [0,1] .. {}",
            sites::fmt_e(numbers.t2az)
        );
    }
    if !(0.0..=1.0).contains(&numbers.t2hw) {
        bail!(
            "Error: The parameter --az-to-hw is not in [0,1] .. {}",
            sites::fmt_e(numbers.t2hw)
        );
    }
    if naf_opts > 1 {
        bail!("Error: The options --AF-tag, --AF-file and -e are mutually exclusive");
    }
    let mut sr = open_input(&args, matches, &numbers, &fname, threads)?;

    let samples = later(
        matches,
        ("samples", args.samples.as_deref()),
        ("samples_file", args.samples_file.as_deref()),
    );
    let setup = Setup {
        args: &args,
        numbers: &numbers,
        filter_expr,
        samples,
        output,
    };
    let (config, mut filter) = setup.init(&mut sr)?;

    eprintln!("Number of target samples: {}", config.roh_samples.len());
    let naf_samples = config.naf_samples;
    eprintln!("Number of --estimate-AF samples: {naf_samples}");
    if config.nbuf_max > 0 {
        eprintln!(
            "Number of sites in the buffer/overlap: {}/{}",
            config.nbuf_max, config.nbuf_olap
        );
    } else {
        eprintln!("Number of sites in the buffer/overlap: unlimited");
    }

    let mut out = OutFile::open(output_name, compression)?;
    // A closed pipe passes as is, for `main` to die of SIGPIPE.
    let write_err = |e: io::Error| {
        if e.kind() == io::ErrorKind::BrokenPipe {
            anyhow::Error::from(e)
        } else {
            anyhow::anyhow!("Error writing {output_name}: {}", strerror(&e))
        }
    };
    out.write_all(header_text(output, config.training).as_bytes())
        .map_err(write_err)?;

    let mut roh = Roh::new(config, usize::try_from(threads).unwrap_or(0));
    let af_file = last(&args.af_file);
    if let Err(e) = read_records(&mut sr, &mut roh, filter.as_mut(), af_file, &mut out) {
        // bcftools exits with 1 on a PL field that is not an integer, and with -1 otherwise.
        if let Some(pl) = e.downcast_ref::<PlTypeError>() {
            eprintln!("{pl}");
            return Ok(ExitCode::from(1));
        }
        return Err(match e.downcast::<io::Error>() {
            Ok(io) => write_err(io),
            Err(e) => e,
        });
    }
    roh.finish(sr.header(0), &mut out).map_err(write_err)?;

    // The summary comes before the output is closed (bcftools' `destroy_data`), so it is
    // printed also when the close fails or the reader of the output went away.
    let stats = &roh.stats;
    let nmin = roh.min_used();
    if af_file.is_some() {
        eprintln!(
            "Number of lines overlapping with --AF-file/processed: {}/{nmin}",
            stats.ntot
        );
    } else {
        eprintln!("Number of lines total/processed: {}/{nmin}", stats.ntot);
    }
    eprintln!(
        "Number of lines filtered/no AF/no alt/multiallelic/dup: {}/{}/{}/{}/{}",
        stats.nfiltered, stats.nno_af, stats.nno_alt, stats.nmultiallelic, stats.ndup
    );
    if nmin == 0 {
        eprintln!("No usable sites were found.");
        if naf_opts == 0 && numbers.af_dflt == 0.0 {
            eprintln!(" Consider using one of the AF options.");
        }
    }
    out.close().map_err(|e| {
        if e.kind() == io::ErrorKind::BrokenPipe {
            anyhow::Error::from(e)
        } else {
            anyhow::anyhow!("Error: close failed .. {output_name}")
        }
    })?;
    Ok(ExitCode::SUCCESS)
}

/// Reads the records, filters them, and hands them to the model (`main_vcfroh`'s loop).
fn read_records(
    sr: &mut SyncedReaders,
    roh: &mut Roh,
    mut filter: Option<&mut (BoxedFilter, FilterLogic)>,
    af_file: Option<&str>,
    out: &mut OutFile,
) -> anyhow::Result<()> {
    let mut rec = Record::new();
    loop {
        match sr.next_line() {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => bail!("Error: {e}"),
        }
        sr.swap_line(0, &mut rec);
        roh.stats.ntot += 1;
        let header = sr.header(0);
        if let Some((f, logic)) = filter.as_deref_mut() {
            let mut pass = f.pass(header, &mut rec)?;
            if *logic == FilterLogic::Exclude {
                pass = !pass;
            }
            if !pass {
                roh.stats.nfiltered += 1;
                continue;
            }
        }
        let target = af_file.zip(sr.targets()).map(|(fname, t)| sites::AfTarget {
            fname,
            line: t.current_line().unwrap_or_default(),
            alleles: t.alleles(),
        });
        roh.site(header, &mut rec, target, out)?;
    }
    Ok(())
}

/// Everything `init_data` needs.
struct Setup<'a> {
    args: &'a Args,
    numbers: &'a LoopOptions,
    filter_expr: Option<(&'a str, FilterLogic)>,
    samples: Option<(&'a str, ListSource)>,
    output: Output,
}

type BoxedFilter = Box<dyn filter::RecordFilter + Send>;

impl Setup<'_> {
    /// Checks the header, reads the sample lists, and subsets the samples of VCF input
    /// (`init_data`).
    fn init(
        &self,
        sr: &mut SyncedReaders,
    ) -> anyhow::Result<(Config, Option<(BoxedFilter, FilterLogic)>)> {
        let args = self.args;
        let hdr = sr.header(0);
        if hdr.sample_count() == 0 {
            bail!("No samples in the VCF?");
        }
        let filter = match self.filter_expr {
            Some((expr, logic)) => Some((filter::init(hdr, expr)?, logic)),
            None => None,
        };
        let fake_pls = self.numbers.unseen_pl.is_some();
        let mut pl_id = None;
        if !fake_pls {
            let id = format_id(hdr, b"PL", "PL", " consider running with -G")?;
            if hdr.id_type(HeaderRecordType::Format, id) != Some(ValueType::Integer) {
                bail!("Error: The FORMAT/PL tag not defined as Integer in the header");
            }
            pl_id = Some(id);
        }
        let mut gt_id = None;

        // -e [GT,|PL,]FILE or -.
        let mut estimate: Option<(EstimateFrom, Option<&str>)> = None;
        if let Some(e) = last(&args.estimate_af) {
            let (from, list) = if let Some(rest) = e.strip_prefix("GT,") {
                (EstimateFrom::Gt, rest)
            } else if let Some(rest) = e.strip_prefix("PL,") {
                (EstimateFrom::Pl, rest)
            } else {
                (EstimateFrom::Gt, e)
            };
            estimate = Some((from, (list != "-").then_some(list)));
        }
        let af_list = estimate.and_then(|(_, list)| list);
        let read_af = |hdr: &Header, verbosity: Verbosity| {
            af_list
                .map(|l| SampleList::new(hdr, Some((l, ListSource::File)), verbosity))
                .transpose()
        };
        let read_roh =
            |hdr: &Header, verbosity: Verbosity| SampleList::new(hdr, self.samples, verbosity);
        let mut af_smpl = read_af(hdr, Verbosity::Verbose)?;
        let af_from_pl = estimate.is_some_and(|(from, _)| from == EstimateFrom::Pl);
        if estimate.is_some() || fake_pls {
            if af_from_pl {
                pl_id = Some(format_id(hdr, b"PL", "PL", "")?);
            } else {
                gt_id = Some(format_id(hdr, b"GT", "GT", "")?);
            }
        }
        if fake_pls {
            gt_id = Some(format_id(hdr, b"GT", "GT", "")?);
        }

        let mut roh_smpl = read_roh(hdr, Verbosity::Verbose)?;
        // A VCF is parsed faster when only the samples used are kept. With -e -, AF comes
        // from all samples, so they are all kept; a filter on FORMAT fields was compiled for
        // all samples.
        let filter_reads_fmt = filter
            .as_ref()
            .is_some_and(|(f, _)| f.max_unpack().contains(pasta::variant::Unpack::FMT));
        let all_for_af = estimate.is_some() && af_list.is_none();
        if self.samples.is_some()
            && sr.file_mut(0).format() == InputFormat::Vcf
            && !all_for_af
            && !filter_reads_fmt
            && !roh_smpl.is_empty()
            && subset_samples(sr, &roh_smpl, af_smpl.as_ref())?
        {
            roh_smpl = read_roh(sr.header(0), Verbosity::Quiet)?;
            af_smpl = read_af(sr.header(0), Verbosity::Quiet)?;
        }
        let hdr = sr.header(0);
        // When all samples are listed, no lookup is needed.
        if af_smpl
            .as_ref()
            .is_some_and(|l| l.len() == hdr.sample_count())
        {
            af_smpl = None;
        }

        let (nbuf_max, nbuf_olap) = match &args.buffer_size {
            Some(spec) => parse_buffer_size(spec, roh_smpl.len())?,
            None => (0, 0),
        };

        let naf_samples = match (&af_smpl, estimate) {
            (Some(l), _) => l.len(),
            (None, Some(_)) => hdr.sample_count(),
            (None, None) => 0,
        };
        let numbers = self.numbers;
        let af_source = if let Some(tag) = last(&args.af_tag) {
            AfSource::Tag(tag.to_owned())
        } else if !args.af_file.is_empty() {
            AfSource::File
        } else if numbers.af_dflt > 0.0 {
            AfSource::Default
        } else if let Some((from, _)) = estimate {
            AfSource::Estimate {
                from,
                samples: af_smpl.map(SampleList::into_vec),
            }
        } else {
            AfSource::AcAn
        };
        let config = Config {
            roh_samples: roh_smpl.into_vec(),
            naf_samples,
            af_source,
            af_dflt: numbers.af_dflt,
            unseen_pl: numbers.unseen_pl,
            pl_id,
            gt_id,
            skip_homref: args.ignore_homref,
            snps_only: args.skip_indels,
            include_noalt: args.include_noalt,
            t2az: numbers.t2az,
            t2hw: numbers.t2hw,
            training: numbers.training.map(|threshold| Training { threshold }),
            nbuf_max,
            nbuf_olap,
            output: self.output,
            genmap: args.genetic_map.clone(),
            rec_rate: numbers.rec_rate,
        };
        Ok((config, filter))
    }
}

/// The dictionary ID of a FORMAT tag, or bcftools' error when the header does not define it.
fn format_id(hdr: &Header, tag: &[u8], name: &str, hint: &str) -> anyhow::Result<u32> {
    match hdr
        .id(tag)
        .filter(|&id| hdr.id_exists(HeaderRecordType::Format, id))
    {
        Some(id) => Ok(id),
        None if hint.is_empty() => bail!("Error: The FORMAT/{name} tag not found in the header"),
        None => bail!("Error: The FORMAT/{name} tag not found in the header,{hint}"),
    }
}

/// Keeps only the samples of the two lists in a VCF input's records, which makes parsing
/// faster; `false` when they are all the samples.
fn subset_samples(
    sr: &mut SyncedReaders,
    roh_smpl: &SampleList,
    af_smpl: Option<&SampleList>,
) -> anyhow::Result<bool> {
    let hdr = sr.header(0);
    let mut keep = vec![false; hdr.sample_count()];
    for &i in roh_smpl.iter().chain(af_smpl.iter().flat_map(|l| l.iter())) {
        keep[i] = true;
    }
    if keep.iter().all(|&k| k) {
        return Ok(false);
    }
    let names: Vec<&[u8]> = hdr
        .samples()
        .iter()
        .zip(&keep)
        .filter(|&(_, &k)| k)
        .map(|(s, _)| &s[..])
        .collect();
    let list = names.join(&b","[..]);
    let ret = sr
        .header_mut(0)
        .set_samples(Some(&list), ListSource::Inline)
        .map_err(|_| {
            anyhow::anyhow!(
                "Error parsing the list of samples: {}",
                String::from_utf8_lossy(&list)
            )
        })?;
    if ret > 0 {
        bail!(
            "The {ret}-th sample not found in the VCF: {}",
            String::from_utf8_lossy(&list)
        );
    }
    Ok(true)
}
