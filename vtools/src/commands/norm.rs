//! `vtools norm`: left-align and normalize indels, check REF alleles, split multiallelic sites
//! and join biallelic ones, atomize complex variants, and remove duplicates (bcftools norm,
//! `vcfnorm.c`).
//!
//! The structure follows bcftools: records are normalized as they are read and kept in a
//! buffer sorted by position (realignment can move them), which is flushed once records are
//! far enough behind the newest one. Flushing writes them, joining (`-m+`) or deduplicating
//! (`-d`) records at the same position on the way.

mod atomize;
pub(crate) mod emit;
mod gff;
pub(crate) mod input;
mod join;
mod ploidy;
mod realign;
pub(crate) mod reference;
mod split;
mod values;

use std::{
    collections::VecDeque,
    io::{self, IsTerminal},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, HeaderRecordType, Number, Record, Unpack, VariantType,
    io::{OutputType, WriteError, Writer},
    number,
    synced::{Overlap, Regions, SyncedOptions, SyncedReader},
};

use self::{
    atomize::{Atomized, Atomizer},
    emit::Emitter,
    gff::Transcripts,
    join::Joiner,
    reference::Reference,
    values::Scratch,
};
use crate::options::{last_index, later_file_sticky};
use crate::{
    commands::head::parse_verbosity,
    filter::{self, FilterLogic, RecordFilter},
    output::{self, OutputSetting},
    version,
};

/// bcftools' usage text (`usage()`), printed for `-h`, `--help`, and no input on a terminal.
const USAGE: &str = "
About:   Left-align and normalize indels; check if REF alleles match the reference;
         split multiallelic sites into multiple rows; recover multiallelics from
         multiple rows.
Usage:   vtools norm [options] <in.vcf.gz>

Options:
    -a, --atomize                   Decompose complex variants (e.g. MNVs become consecutive SNVs)
        --atom-overlaps '*'|.       Use the star allele (*) for overlapping alleles or set to missing (.) [*]
    -c, --check-ref e|w|x|s         Check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [e]
    -D, --remove-duplicates         Remove duplicate lines of the same type.
    -d, --rm-dup TYPE               Remove duplicate snps|indels|both|all|exact
    -e, --exclude EXPR              Do not normalize records for which the expression is true (see man page for details)
    -f, --fasta-ref FILE            Reference sequence
        --force                     Try to proceed even if malformed tags are encountered. Experimental, use at your own risk
    -g, --gff-annot FILE            Follow HGVS 3'rule and right-align variants in transcripts on the forward strand
    -i, --include EXPR              Normalize only records for which the expression is true (see man page for details)
        --keep-sum TAG,..           Keep vector sum constant when splitting multiallelics (see github issue #360)
    -m, --multiallelics -|+TYPE     Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]
        --multi-overlaps 0|.        Fill in the reference (0) or missing (.) allele when splitting multiallelics [0]
        --no-version                Do not append version and command line to the header
    -N, --no-realign [NUM]          Skip realignment; optionally only for events longer than NUM bp (e.g. -N1000, no space)
        --old-rec-tag STR           Annotate modified records with INFO/STR indicating the original variant
    -o, --output FILE               Write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -r, --regions REGION            Restrict to comma-separated list of regions
    -R, --regions-file FILE         Restrict to regions listed in a file
        --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --strict-filter             When merging (-m+), merged site is PASS only if all sites being merged PASS
    -S, --sort METHOD               Sort order: chr_pos,lex [chr_pos]
    -t, --targets REGION            Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE         Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT               Use multithreading with INT worker threads [0]
    -v, --verbosity INT             Verbosity level
    -w, --site-win INT              Buffer for sorting lines which changed position during realignment [1000]
    -W, --write-index[=FMT]         Automatically index the output files [off]

Examples:
   # normalize and left-align indels
   vtools norm -f ref.fa in.vcf

   # split multi-allelic sites
   vtools norm -m- in.vcf

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin); further arguments are ignored, as in
    /// bcftools.
    inputs: Vec<String>,

    /// Decompose complex variants (e.g. MNVs become consecutive SNVs).
    #[arg(short = 'a', long = "atomize")]
    atomize: bool,

    /// Use the star allele (*) for overlapping alleles or set to missing (.) [*].
    #[arg(long = "atom-overlaps", value_name = "'*'|.", action = clap::ArgAction::Append)]
    atom_overlaps: Vec<String>,

    /// Check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [e].
    #[arg(short = 'c', long = "check-ref", value_name = "e|w|x|s", action = clap::ArgAction::Append)]
    check_ref: Vec<String>,

    /// Remove duplicate lines of the same type (deprecated alias of `-d exact`).
    #[arg(short = 'D', long = "remove-duplicates", action = clap::ArgAction::Count)]
    remove_duplicates: u8,

    /// Remove duplicate snps|indels|both|all|exact.
    #[arg(short = 'd', long = "rm-dup", value_name = "TYPE", action = clap::ArgAction::Append)]
    rm_dup: Vec<String>,

    /// Do not normalize records for which the expression is true.
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Reference sequence.
    #[arg(short = 'f', long = "fasta-ref", value_name = "FILE")]
    fasta_ref: Option<PathBuf>,

    /// Try to proceed even if malformed tags are encountered.
    #[arg(long = "force")]
    force: bool,

    /// Follow HGVS 3'rule and right-align variants in transcripts on the forward strand.
    #[arg(short = 'g', long = "gff-annot", value_name = "FILE")]
    gff_annot: Option<PathBuf>,

    /// Normalize only records for which the expression is true.
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

    /// Keep vector sum constant when splitting multiallelics.
    #[arg(long = "keep-sum", value_name = "TAG,..", action = clap::ArgAction::Append)]
    keep_sum: Vec<String>,

    /// Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both].
    #[arg(
        short = 'm',
        long = "multiallelics",
        value_name = "-|+TYPE",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    multiallelics: Vec<String>,

    /// Fill in the reference (0) or missing (.) allele when splitting multiallelics [0].
    #[arg(long = "multi-overlaps", value_name = "0|.", action = clap::ArgAction::Append)]
    multi_overlaps: Vec<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Skip realignment; optionally only for events longer than NUM bp (e.g. -N1000).
    #[arg(
        short = 'N',
        long = "no-realign",
        value_name = "NUM",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "",
        action = clap::ArgAction::Append
    )]
    no_realign: Vec<String>,

    /// Skip realignment (as `-N`).
    #[arg(long = "do-not-normalize")]
    do_not_normalize: bool,

    /// Right-align instead of left-align (for debugging).
    #[arg(long = "right-align", hide = true)]
    right_align: bool,

    /// Annotate modified records with INFO/STR indicating the original variant.
    #[arg(long = "old-rec-tag", value_name = "STR")]
    old_rec_tag: Option<String>,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "u|b|v|z[0-9]", action = clap::ArgAction::Append)]
    output_type: Vec<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(short = 'r', long = "regions", value_name = "REGION")]
    regions: Option<String>,

    /// Restrict to regions listed in a file.
    #[arg(short = 'R', long = "regions-file", value_name = "FILE")]
    regions_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(long = "regions-overlap", value_name = "0|1|2", action = clap::ArgAction::Append)]
    regions_overlap: Vec<String>,

    /// When merging (-m+), merged site is PASS only if all sites being merged PASS.
    #[arg(short = 's', long = "strict-filter")]
    strict_filter: bool,

    /// Sort order: pos,lex [pos].
    #[arg(short = 'S', long = "sort", value_name = "METHOD", action = clap::ArgAction::Append)]
    sort: Vec<String>,

    /// Similar to -r but streams rather than index-jumps.
    #[arg(
        short = 't',
        long = "targets",
        value_name = "REGION",
        allow_hyphen_values = true
    )]
    targets: Option<String>,

    /// Similar to -R but streams rather than index-jumps.
    #[arg(short = 'T', long = "targets-file", value_name = "FILE")]
    targets_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0].
    #[arg(long = "targets-overlap", value_name = "0|1|2", action = clap::ArgAction::Append)]
    targets_overlap: Vec<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        visible_alias = "verbose",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,

    /// Buffer for sorting lines which changed position during realignment [1000].
    #[arg(
        short = 'w',
        long = "site-win",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    site_win: Vec<String>,

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

/// `-c` bits (bcftools' `CHECK_REF_*`).
const CHECK_REF_EXIT: u8 = 1;
const CHECK_REF_WARN: u8 = 2;
const CHECK_REF_SKIP: u8 = 4;
const CHECK_REF_FIX: u8 = 8;

/// The most records kept for reuse; more are dropped.
const MAX_POOLED_RECORDS: usize = 4096;

/// `-d` bits (HTSlib's `BCF_SR_PAIR_*`).
const PAIR_SNPS: u8 = 1;
const PAIR_INDELS: u8 = 2;
const PAIR_ANY: u8 = 4;
const PAIR_EXACT: u8 = 64;

/// What `-m` does with multiallelic sites.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MultiRows {
    /// `-m-`: split them.
    Split,
    /// `-m+`: join records at one position.
    Join,
}

/// Which variant types `-m` splits or joins (bcftools' `COLLAPSE_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Collapse {
    Snps,
    Indels,
    Both,
    Any,
}

/// The order of records at one position (`-S`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SortOrder {
    /// By position only (`cmp_bcf_pos`).
    Pos,
    /// By position, then REF and ALT (`cmp_bcf_pos_ref_alt`).
    Lex,
}

/// The options that drive normalization.
#[derive(Debug)]
struct Options {
    check_ref: u8,
    multi_rows: Option<MultiRows>,
    collapse: Collapse,
    rmdup: u8,
    atomize: bool,
    use_star_allele: bool,
    ma_use_ref_allele: bool,
    old_rec_tag: Option<Vec<u8>>,
    /// The header ID of AD with `--keep-sum AD`.
    keep_sum_ad: Option<i32>,
    strict_filter: bool,
    force: bool,
    /// Realign records with `rlen` below this (`-N`).
    do_realign: u32,
    /// The realignment window (bcftools' `aln_win`).
    aln_win: i64,
    /// The sorting window (`-w`).
    buf_win: i64,
    right_align: bool,
    sort: SortOrder,
}

/// The counts printed at the end.
#[derive(Debug, Default)]
struct Stats {
    total: u64,
    split: u64,
    joined: u64,
    changed: u64,
    skipped: u64,
    rmdup: u64,
    filtered: u64,
    ref_total: u64,
    ref_set: u64,
    ref_swap: u64,
}

/// A buffered record and whether it passed `-i/-e`.
struct Line {
    rec: Record,
    pass: bool,
}

/// The alleles of a record written at the current position, for `-d exact` (bcftools'
/// `cmpals1_t`).
#[derive(Default)]
struct CmpAls {
    n: usize,
    ref_: Vec<u8>,
    alt: Vec<u8>,
    alts: Vec<Vec<u8>>,
}

/// The state of `-d` within one flush (bcftools' locals of `flush_buffer`).
struct DupState {
    prev: Option<(i32, i64)>,
    prev_type: VariantType,
}

/// The normalization state (bcftools' `args_t`, minus the reader and writer).
struct Normalizer {
    opts: Options,
    stats: Stats,
    fai: Option<Reference>,
    gff: Option<Transcripts>,
    filter: Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    filter_pass: bool,
    /// The `--old-rec-tag` value being built (bcftools' `old_rec_tag_kstr`).
    old_rec: Vec<u8>,
    force_warned: bool,
    scratch: Scratch,
    seqs: realign::Seqs,
    /// The records waiting to be written, sorted by position (bcftools' `rbuf`/`lines`).
    lines: VecDeque<Line>,
    /// Records to reuse.
    pool: Vec<Record>,
    /// The biallelic records of a split site.
    split_lines: Vec<Record>,
    joiner: Joiner,
    atomizer: Atomizer,
    cmpals: Vec<CmpAls>,
    ncmpals: usize,
    emitter: Emitter,
    /// The parallel VCF parser, with `--threads` on a plain VCF stream.
    input: Option<input::ParallelVcf>,
    svlen_warned: bool,
}

/// How a record's normalization ended (bcftools' `ERR_*` codes of `realign`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Realigned {
    DupAllele,
    RefMismatch,
    Ok,
    Symbolic,
    SpanningDeletion,
}

/// Whether `normalize_line` kept the record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Kept {
    Yes,
    /// Dropped by `-c x`.
    Skipped,
}

/// Collects `-O` occurrences in command-line order.
fn output_settings(args: &Args) -> Vec<OutputSetting> {
    args.output_type
        .iter()
        .map(|v| OutputSetting::Type(v.clone()))
        .collect()
}

/// C's `strtol(s, &end, 10)` requiring the whole string to be a number.
fn strtol_strict(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    let (v, end) = pasta::variant::number::strtol(b, 0);
    (end == b.len()).then_some(v)
}

/// Parses every `--regions-overlap` or `--targets-overlap` given; the last one counts.
fn parse_overlap(values: &[String], default: Overlap, what: &str) -> anyhow::Result<Overlap> {
    values
        .iter()
        .try_fold(default, |_, s| match Overlap::parse(s) {
            Some(o) => Ok(o),
            None => bail!("Could not parse: --{what} {s}"),
        })
}

/// Parses a `-N` value (empty for a bare `-N`) as bcftools does, with `strtod`: realign only
/// records shorter than the value. NaN passes bcftools' range check and converts to 0.
fn parse_no_realign(n: &str) -> anyhow::Result<u32> {
    if n.is_empty() {
        return Ok(0);
    }
    let b = n.as_bytes();
    let (val, end) = pasta::variant::number::strtod(b, 0);
    if val.is_nan() && end == b.len() {
        return Ok(0);
    }
    if end == 0 || end != b.len() || !(0.0..=f64::from(u32::MAX)).contains(&val) {
        bail!("Could not parse argument: -N, --no-realign {n}");
    }
    Ok(val as u32)
}

/// What bcftools' `init_data` sets up besides the reader and the output: the reference, the
/// transcripts, the `-i/-e` filter, and the output header.
struct Init {
    fai: Option<Reference>,
    gff: Option<Transcripts>,
    filter: Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    out_hdr: Header,
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    // Every occurrence of an option is checked, as getopt hands each one to bcftools; the
    // last one wins.
    let mut verbose = 1;
    for v in &args.verbosity {
        verbose = parse_verbosity(v)?;
    }
    if args.help > 0 {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    // bcftools checks the -W format while parsing the options.
    let write_index = args
        .write_index
        .iter()
        .try_fold(None, |_, w| output::parse_write_index(w).map(Some))?;
    let keep_sum = args.keep_sum.iter().try_fold(false, |_, k| {
        if k.eq_ignore_ascii_case("ad") {
            Ok(true)
        } else {
            Err(anyhow::anyhow!(
                "Error: only --keep-sum AD is currently supported. See https://github.com/samtools/bcftools/issues/360 for more."
            ))
        }
    })?;
    let filter_expr =
        FilterLogic::single(&args.include, &args.exclude)?.map(|(e, logic)| (e.to_owned(), logic));
    let mut opts = parse_options(&args, matches)?;
    // strtol: a negative count fails when the thread pool is made.
    let threads = args
        .threads
        .as_deref()
        .map_or(Some(0), |t| usize::try_from(output::thread_count(t)).ok());
    let regions_overlap = parse_overlap(&args.regions_overlap, Overlap::Record, "regions-overlap")?;
    let targets_overlap = parse_overlap(&args.targets_overlap, Overlap::Pos, "targets-overlap")?;

    let fname = match args.inputs.first() {
        Some(f) => f.clone(),
        None if io::stdin().is_terminal() => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
        None => "-".to_owned(),
    };
    if args.fasta_ref.is_none() && opts.multi_rows.is_none() && opts.rmdup == 0 && !opts.atomize {
        bail!("Expected -a, -f, -m, -D or -d option");
    }
    if opts.check_ref == 0 && args.fasta_ref.is_some() {
        opts.check_ref = CHECK_REF_EXIT;
    }
    if opts.check_ref != 0 && args.fasta_ref.is_none() {
        bail!("Expected --fasta-ref with --check-ref");
    }

    let synced = synced_options(&args, matches, regions_overlap, targets_overlap, threads)?;
    let threads = synced.threads;
    let mut reader = match SyncedReader::open(&fname, synced) {
        Ok(r) => r,
        Err(e) => {
            let name = if fname == "-" {
                "standard input"
            } else {
                fname.as_str()
            };
            bail!("Failed to read from {name}: {e}");
        }
    };

    let Init {
        fai,
        gff,
        filter,
        mut out_hdr,
    } = init_data(
        &args,
        &mut opts,
        &mut reader,
        keep_sum,
        verbose,
        filter_expr.as_ref(),
    )?;
    let output_name = args
        .output
        .as_deref()
        .map_or_else(|| "-".to_owned(), |p| p.display().to_string());
    let (mut writer, out_type) =
        open_writer(&args, &mut out_hdr, threads, write_index, &output_name)?;

    let mut norm = Normalizer {
        opts,
        stats: Stats::default(),
        fai,
        gff,
        filter,
        filter_pass: true,
        old_rec: Vec::new(),
        force_warned: false,
        scratch: Scratch::default(),
        seqs: realign::Seqs::default(),
        lines: VecDeque::new(),
        pool: Vec::new(),
        split_lines: Vec::new(),
        joiner: Joiner::default(),
        atomizer: Atomizer::default(),
        cmpals: Vec::new(),
        ncmpals: 0,
        emitter: Emitter::new(&out_hdr, &writer, threads, output_name.clone()),
        input: input::ParallelVcf::new(&mut reader, threads),
        svlen_warned: false,
    };
    let result = norm.normalize_vcf(&mut reader, &out_hdr, &mut writer);
    if let Err(e) = result {
        // bcftools exits from error() without flushing its output buffers; with threads,
        // write first what a serial run would have written.
        norm.emitter.abort(&out_hdr, &mut writer);
        writer.abandon();
        return Err(e);
    }
    close_writer(writer, &args, out_type, &output_name)?;
    norm.print_stats();
    Ok(ExitCode::SUCCESS)
}

/// Parses the options that drive normalization, each occurrence in turn as getopt hands them
/// to bcftools.
fn parse_options(args: &Args, matches: &ArgMatches) -> anyhow::Result<Options> {
    let mut opts = Options {
        check_ref: 0,
        multi_rows: None,
        collapse: Collapse::Both,
        rmdup: 0,
        atomize: args.atomize,
        use_star_allele: true,
        ma_use_ref_allele: true,
        old_rec_tag: args.old_rec_tag.as_ref().map(|s| s.as_bytes().to_vec()),
        keep_sum_ad: None,
        strict_filter: args.strict_filter,
        force: args.force,
        do_realign: u32::MAX,
        aln_win: 100,
        buf_win: 1000,
        right_align: args.right_align,
        sort: SortOrder::Pos,
    };
    for s in &args.sort {
        opts.sort = if s.eq_ignore_ascii_case("pos") {
            SortOrder::Pos
        } else if s.eq_ignore_ascii_case("lex") {
            SortOrder::Lex
        } else {
            bail!("Error: the sort order --sort {s} is not recognised");
        };
    }
    for a in &args.atom_overlaps {
        opts.use_star_allele = match a.as_bytes().first() {
            Some(b'*') => true,
            Some(b'.') => false,
            _ => bail!(
                "Invalid argument to --atom-overlaps. Perhaps you wanted: \"--atom-overlaps '*'\"?"
            ),
        };
    }
    for m in &args.multi_overlaps {
        opts.ma_use_ref_allele = match m.as_bytes().first() {
            Some(b'0') => true,
            Some(b'.') => false,
            _ => bail!("Invalid argument to --multi-overlaps"),
        };
    }
    for n in &args.no_realign {
        opts.do_realign = parse_no_realign(n)?;
    }
    // --do-not-normalize after the last -N turns realignment off again.
    let last = |id: &str| last_index(matches, id);
    if args.do_not_normalize && last("do_not_normalize") > last("no_realign") {
        opts.do_realign = 0;
    }
    opts.rmdup = parse_rm_dup(args, matches)?;
    // A later -m without a type keeps the type of an earlier one.
    for m in &args.multiallelics {
        let b = m.as_bytes();
        opts.multi_rows = Some(match b.first() {
            Some(b'-') => MultiRows::Split,
            Some(b'+') => MultiRows::Join,
            _ => bail!("Expected '+' or '-' with -m"),
        });
        if b.len() > 1 {
            opts.collapse = match &m[1..] {
                "snps" => Collapse::Snps,
                "indels" => Collapse::Indels,
                "both" => Collapse::Both,
                "any" => Collapse::Any,
                _ => bail!("The argument to -m not recognised: {m}"),
            };
        }
    }
    for c in &args.check_ref {
        if c.contains('w') {
            opts.check_ref |= CHECK_REF_WARN;
        }
        if c.contains('x') {
            opts.check_ref |= CHECK_REF_SKIP;
        }
        if c.contains('s') {
            opts.check_ref |= CHECK_REF_FIX;
        }
        if c.contains('e') {
            opts.check_ref = CHECK_REF_EXIT;
        }
    }
    for w in &args.site_win {
        let Some(v) = strtol_strict(w) else {
            bail!("Could not parse argument: --site-win {w}");
        };
        opts.buf_win = v;
    }
    Ok(opts)
}

/// Parses `-d` and `-D` in command-line order into `PAIR_*` bits; the last one wins (0 when
/// neither is given).
fn parse_rm_dup(args: &Args, matches: &ArgMatches) -> anyhow::Result<u8> {
    let mut dups: Vec<(usize, Option<&str>)> = Vec::new();
    if let Some(i) = matches.indices_of("rm_dup") {
        dups.extend(i.zip(args.rm_dup.iter().map(|d| Some(d.as_str()))));
    }
    if args.remove_duplicates > 0
        && let Some(i) = matches.indices_of("remove_duplicates")
    {
        dups.extend(i.map(|i| (i, None)));
    }
    dups.sort_unstable();
    let mut rmdup = 0;
    for (_, d) in dups {
        rmdup = match d {
            None => {
                eprintln!(
                    "Warning: `-D` is functional but deprecated, replaced by and alias of `-d exact`."
                );
                PAIR_EXACT
            }
            Some("snps") => PAIR_SNPS,
            Some("indels") => PAIR_INDELS,
            Some("both") => PAIR_SNPS | PAIR_INDELS,
            Some("all" | "any") => PAIR_ANY,
            Some("none" | "exact") => PAIR_EXACT,
            Some(d) => bail!("The argument to -d not recognised: {d}"),
        };
    }
    Ok(rmdup)
}

/// Reads `-r/-R` and `-t/-T` into the reader's options, and checks the thread count
/// (`None`: one that cannot make a thread pool).
fn synced_options(
    args: &Args,
    matches: &ArgMatches,
    regions_overlap: Overlap,
    targets_overlap: Overlap,
    threads: Option<usize>,
) -> anyhow::Result<SyncedOptions> {
    let mut regions = None;
    if let Some((spec, source)) = later_file_sticky(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    ) {
        let mut r = Regions::new(spec.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = regions_overlap;
        regions = Some(r);
    }
    let mut targets = None;
    let mut targets_exclude = false;
    if let Some((spec, source)) = later_file_sticky(
        matches,
        ("targets", args.targets.as_deref()),
        ("targets_file", args.targets_file.as_deref()),
    ) {
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
    let Some(threads) = threads.filter(|&t| output::threads_start(t)) else {
        bail!("Failed to create threads");
    };
    Ok(SyncedOptions {
        regions,
        targets,
        targets_exclude,
        apply_filters: None,
        threads,
    })
}

/// bcftools' `init_data`, but for opening the output: adds `--old-rec-tag` to the header, and
/// opens the reference, the transcripts, and the filter.
fn init_data(
    args: &Args,
    opts: &mut Options,
    reader: &mut SyncedReader,
    keep_sum: bool,
    verbose: i32,
    filter_expr: Option<&(String, FilterLogic)>,
) -> anyhow::Result<Init> {
    if keep_sum {
        let Some(id) = reader.header().id(b"AD") else {
            bail!("Error: --keep-sum-ad requested but the tag AD is not present");
        };
        opts.keep_sum_ad = Some(id as i32);
    }
    if let Some(tag) = &opts.old_rec_tag {
        let line = format!(
            "##INFO=<ID={},Number=1,Type=String,Description=\"Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX\">",
            String::from_utf8_lossy(tag)
        );
        reader
            .header_mut()
            .append(line.as_bytes())
            .with_context(|| format!("Could not add the header line: {line}"))?;
    }
    let fai = match &args.fasta_ref {
        Some(path) => Some(Reference::open(path)?),
        None => None,
    };
    if opts.multi_rows == Some(MultiRows::Split) {
        let h = reader.header();
        if let Some(id) = h.id(b"SVLEN")
            && h.id_number(HeaderRecordType::Info, id) != Some(Number::A)
        {
            eprintln!(
                "Warning: the tag INFO/SVLEN must be defined as Number=A in order for the field to be split\n         (the command `vtools reheader` can be used to fix the header)"
            );
        }
    }
    let gff = match &args.gff_annot {
        Some(path) => Some(Transcripts::parse(path, verbose)?),
        None => None,
    };
    let filter = match filter_expr {
        Some((expr, logic)) => Some((filter::init(reader.header(), expr)?, *logic)),
        None => None,
    };
    let out_hdr = reader
        .header()
        .dup()
        .context("Failed to duplicate the header")?;
    Ok(Init {
        fai,
        gff,
        filter,
        out_hdr,
    })
}

/// Opens the output, adds the version lines to its header, sets up `-W`, and writes the
/// header. Returns the writer and the output type.
fn open_writer(
    args: &Args,
    out_hdr: &mut Header,
    threads: usize,
    write_index: Option<output::WriteIndex>,
    output_name: &str,
) -> anyhow::Result<(Writer, OutputType)> {
    let (out_type, level) = output::resolve(&output_settings(args), args.output.as_deref())?;
    let out_path = args.output.as_deref().map(output::data_path);
    let mut writer = match Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, threads),
    ) {
        Ok(w) => w,
        Err(e) => bail!("Can't write to \"{output_name}\": {e}"),
    };
    if !args.no_version {
        version::append_version(out_hdr, "norm")?;
    }
    if let Some(fmt) = write_index {
        let ok = output::index_request(fmt, out_type, args.output.as_deref())
            .is_some_and(|req| writer.request_index(req).is_ok());
        if !ok {
            if writer.write_header(out_hdr).is_err() {
                bail!("[normalize_vcf] Error: cannot write to {output_name}");
            }
            writer.abandon();
            bail!("Error: failed to initialise index for {output_name}");
        }
    }
    if writer.write_header(out_hdr).is_err() {
        bail!("[normalize_vcf] Error: cannot write to {output_name}");
    }
    Ok((writer, out_type))
}

/// Closes the output and its index.
fn close_writer(
    writer: Writer,
    args: &Args,
    out_type: OutputType,
    output_name: &str,
) -> anyhow::Result<()> {
    let write_index = !args.write_index.is_empty();
    match writer.finish() {
        Ok(()) => Ok(()),
        Err(WriteError::Index(_)) if write_index => {
            bail!(
                "Error: cannot write to index {}",
                index_name(args, out_type)
            )
        }
        Err(WriteError::Deferred(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("[flush_buffer] Error: cannot write to {output_name}")
        })),
        Err(WriteError::Io(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("[normalize_vcf] Error: close failed .. {output_name}")
        })),
        Err(e) => Err(anyhow::Error::new(e).context(format!(
            "[normalize_vcf] Error: close failed .. {output_name}"
        ))),
    }
}

/// The index file name, for messages.
fn index_name(args: &Args, ty: OutputType) -> String {
    let fmt = args
        .write_index
        .last()
        .and_then(|w| output::parse_write_index(w).ok())
        .unwrap_or(output::WriteIndex::Csi);
    output::index_request(fmt, ty, args.output.as_deref())
        .map_or_else(|| "(null)".to_owned(), |r| r.path.display().to_string())
}

/// Compares two records by position (and alleles for `-S lex`); bcftools' `cmp_bcf_pos` and
/// `cmp_bcf_pos_ref_alt`.
fn compare(order: SortOrder, a: &Record, b: &Record) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let o = (a.rid(), a.pos()).cmp(&(b.rid(), b.pos()));
    if o != Ordering::Equal || order == SortOrder::Pos {
        return o;
    }
    for i in 0..a.allele_count() {
        if i >= b.allele_count() {
            return Ordering::Greater;
        }
        let o = strcasecmp(a.allele(i), b.allele(i));
        if o != Ordering::Equal {
            return o;
        }
    }
    a.allele_count().cmp(&b.allele_count())
}

/// C's `strcasecmp` ordering.
fn strcasecmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    let la = a.iter().map(u8::to_ascii_lowercase);
    let lb = b.iter().map(u8::to_ascii_lowercase);
    la.cmp(lb)
}

/// Whether two alleles are equal ignoring case (`!strcasecmp`).
fn eq_nocase(a: &[u8], b: &[u8]) -> bool {
    a.eq_ignore_ascii_case(b)
}

impl Normalizer {
    /// Takes a cleared record from the pool.
    fn take_record(&mut self) -> Record {
        let mut r = self
            .pool
            .pop()
            .or_else(|| self.emitter.spare())
            .unwrap_or_default();
        r.clear();
        r
    }

    /// Returns a record to the pool.
    fn recycle(&mut self, rec: Record) {
        if self.pool.len() < MAX_POOLED_RECORDS {
            self.pool.push(rec);
        }
    }

    /// Prints the counts of `main_vcfnorm`'s summary to stderr.
    fn print_stats(&self) {
        let s = &self.stats;
        eprintln!(
            "Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\t{}/{}/{}/{}/{}/{}/{}",
            s.total, s.split, s.joined, s.changed, s.skipped, s.rmdup, s.filtered
        );
        if self.opts.check_ref & CHECK_REF_FIX != 0 {
            eprintln!(
                "REF/ALT total/modified/added:  \t{}/{}/{}",
                s.ref_total, s.ref_swap, s.ref_set
            );
        }
    }

    /// The main loop (bcftools' `normalize_vcf` after the header is written).
    fn normalize_vcf(
        &mut self,
        reader: &mut SyncedReader,
        out_hdr: &Header,
        writer: &mut Writer,
    ) -> anyhow::Result<()> {
        let mut read_error = None;
        loop {
            // Buffer lines until a different position is reached.
            let mut done = false;
            loop {
                loop {
                    match self.split_and_normalize(reader, out_hdr) {
                        Ok(Some(Kept::Skipped)) => continue,
                        Ok(Some(Kept::Yes)) => {}
                        Ok(None) => done = true,
                        Err(ReadOrFatal::Read(e)) => {
                            read_error = Some(e);
                            done = true;
                        }
                        Err(ReadOrFatal::Fatal(e)) => return Err(e),
                    }
                    break;
                }
                if done {
                    break;
                }
                let (Some(first), Some(last)) = (self.lines.front(), self.lines.back()) else {
                    continue;
                };
                if first.rec.rid() != last.rec.rid() || first.rec.pos() != last.rec.pos() {
                    break;
                }
            }
            if done {
                break;
            }
            // How many sites to flush.
            let (Some(first), Some(last)) = (self.lines.front(), self.lines.back()) else {
                continue;
            };
            let (first_rid, last_rid, last_pos) = (first.rec.rid(), last.rec.rid(), last.rec.pos());
            let mut n = 0;
            for line in &self.lines {
                if first_rid != last_rid {
                    // Two chromosomes in the buffer: flush the first one.
                    if line.rec.rid() != first_rid {
                        break;
                    }
                    n += 1;
                    continue;
                }
                // One chromosome: flush only lines unlikely to change order on realigning.
                if last_pos - line.rec.pos() < self.opts.buf_win {
                    break;
                }
                n += 1;
            }
            if n > 0 {
                self.flush_buffer(reader.header(), out_hdr, writer, n)?;
            }
        }
        let n = self.lines.len();
        self.flush_buffer(reader.header(), out_hdr, writer, n)?;
        self.emitter.finish(out_hdr, writer)?;
        if let Some(e) = read_error.or_else(|| reader.take_pending_error()) {
            bail!("Error: {e}");
        }
        Ok(())
    }

    /// Reads and normalizes the next record (bcftools' `split_and_normalize`). Returns `None`
    /// at the end of the input.
    fn split_and_normalize(
        &mut self,
        reader: &mut SyncedReader,
        out_hdr: &Header,
    ) -> Result<Option<Kept>, ReadOrFatal> {
        let mut line = self.take_record();
        let got = match &mut self.input {
            Some(p) => p.next(reader, &mut line, &mut self.pool),
            None => reader.next_record(&mut line),
        };
        match got {
            Ok(true) => {}
            Ok(false) => {
                self.recycle(line);
                return Ok(None);
            }
            Err(e) => {
                self.recycle(line);
                return Err(ReadOrFatal::Read(e));
            }
        }
        let hdr = reader.header();
        self.stats.total += 1;
        if let Some((filter, logic)) = &mut self.filter {
            let pass = filter.pass(hdr, &mut line)?;
            self.filter_pass = match logic {
                FilterLogic::Include => pass,
                FilterLogic::Exclude => !pass,
            };
            if !self.filter_pass {
                self.stats.filtered += 1;
            }
        }
        if self.opts.multi_rows != Some(MultiRows::Split)
            || line.allele_count() <= 2
            || !self.filter_pass
        {
            return Ok(Some(self.normalize_line(hdr, out_hdr, line)?));
        }
        // Any restrictions on variant types to split?
        if self.opts.collapse != Collapse::Both && self.opts.collapse != Collapse::Any {
            let ty = if self.opts.collapse == Collapse::Snps {
                VariantType::SNP
            } else {
                VariantType::INDEL
            };
            if !line.variant_types().intersects(ty) {
                return Ok(Some(self.normalize_line(hdr, out_hdr, line)?));
            }
        }
        self.stats.split += 1;
        self.split_multiallelic(hdr, out_hdr, &mut line)?;
        self.recycle(line);
        let lines = std::mem::take(&mut self.split_lines);
        let mut result = Kept::Yes;
        let mut iter = lines.into_iter();
        for rec in iter.by_ref() {
            let kept = self.normalize_line(hdr, out_hdr, rec)?;
            if kept != Kept::Yes {
                result = kept;
                break;
            }
        }
        for rec in iter {
            self.recycle(rec);
        }
        Ok(Some(result))
    }

    /// Normalizes one record and inserts it, or its atoms, into the sorted buffer (bcftools'
    /// `normalize_line`).
    fn normalize_line(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        mut line: Record,
    ) -> anyhow::Result<Kept> {
        if self.fai.is_some() {
            if self.filter_pass && self.opts.check_ref & CHECK_REF_FIX != 0 {
                self.old_rec_tag_init(hdr, &mut line);
                if self.fix_ref(hdr, out_hdr, &mut line)? {
                    self.old_rec_tag_set(out_hdr, &mut line, 0)?;
                }
            }
            if line.rlen() < i64::from(self.opts.do_realign) {
                let ret = if self.filter_pass {
                    self.realign(hdr, out_hdr, &mut line)?
                } else {
                    Realigned::Ok
                };
                if ret == Realigned::RefMismatch && self.opts.check_ref & CHECK_REF_SKIP != 0 {
                    self.stats.skipped += 1;
                    self.recycle(line);
                    return Ok(Kept::Skipped);
                }
                if ret == Realigned::DupAllele {
                    if self.opts.check_ref & CHECK_REF_FIX != 0 {
                        self.fix_dup_alt(hdr, out_hdr, &mut line);
                    } else if self.opts.check_ref == CHECK_REF_EXIT {
                        bail!(
                            "Duplicate alleles at {}:{}; run with -cw to turn the error into warning or with -cs to fix.",
                            String::from_utf8_lossy(line.chrom_safe(hdr)),
                            line.pos() + 1
                        );
                    } else if self.opts.check_ref & CHECK_REF_WARN != 0 {
                        eprintln!(
                            "ALT_DUP\t{}\t{}",
                            String::from_utf8_lossy(line.chrom_safe(hdr)),
                            line.pos() + 1
                        );
                    }
                }
            }
        }

        if self.filter_pass && self.opts.atomize {
            let mut atoms = std::mem::take(&mut self.atomizer.out);
            let atomized = self.atomize(hdr, &mut line, &mut atoms);
            match atomized {
                Ok(Atomized::Unchanged) => self.insert(line, true),
                Ok(Atomized::Split) => {
                    self.recycle(line);
                    for rec in atoms.drain(..) {
                        self.insert(rec, true);
                    }
                }
                Err(_) => self.recycle(line),
            }
            self.atomizer.out = atoms;
            atomized?;
        } else {
            let pass = self.filter_pass;
            self.insert(line, pass);
        }
        Ok(Kept::Yes)
    }

    /// Appends a record to the buffer, moving it up to its sorted place. bcftools buffers a
    /// `bcf_dup` copy, which loses the 64-bit flag: BCF output then fails only for positions
    /// that BCF cannot store.
    fn insert(&mut self, mut rec: Record, pass: bool) {
        rec.clear_64bit();
        let order = self.opts.sort;
        if order == SortOrder::Lex {
            rec.unpack(Unpack::STR);
        }
        self.lines.push_back(Line { rec, pass });
        let mut k = self.lines.len() - 1;
        while k > 0 {
            let (a, b) = (&self.lines[k - 1].rec, &self.lines[k].rec);
            if a.rid() != b.rid() || compare(order, a, b) != std::cmp::Ordering::Greater {
                // The lines before are sorted already.
                break;
            }
            self.lines.swap(k - 1, k);
            k -= 1;
        }
    }

    /// Writes the first `n` buffered lines, joining or deduplicating them (bcftools'
    /// `flush_buffer`).
    fn flush_buffer(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        writer: &mut Writer,
        n: usize,
    ) -> anyhow::Result<()> {
        let mut dup = DupState {
            prev: None,
            prev_type: VariantType::empty(),
        };
        for _ in 0..n {
            let Some(Line { mut rec, pass }) = self.lines.pop_front() else {
                break;
            };
            if !pass {
                self.write(out_hdr, writer, rec)?;
                continue;
            }
            if self.opts.multi_rows == Some(MultiRows::Join) {
                if self.joiner.can_flush(&rec) {
                    self.flush_joined(hdr, out_hdr, writer)?;
                }
                let old = self.joiner.push(rec, self.opts.collapse);
                if let Some(old) = old {
                    self.recycle(old);
                }
                continue;
            } else if self.opts.rmdup != 0 && self.is_duplicate(hdr, &mut dup, &mut rec) {
                self.stats.rmdup += 1;
                self.recycle(rec);
                continue;
            }
            self.write(out_hdr, writer, rec)?;
        }
        if self.opts.multi_rows == Some(MultiRows::Join) && self.lines.is_empty() {
            self.flush_joined(hdr, out_hdr, writer)?;
        }
        Ok(())
    }

    /// Writes all joined records waiting (bcftools' `while (line = mrows_flush())` loop).
    fn flush_joined(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        writer: &mut Writer,
    ) -> anyhow::Result<()> {
        loop {
            let collapse = self.opts.collapse;
            let Some(which) = self.joiner.flush(collapse) else {
                return Ok(());
            };
            match which {
                join::Flushed::Single(i) => {
                    let spare = self.take_record();
                    let rec = std::mem::replace(&mut self.joiner.rows[i], spare);
                    self.write(out_hdr, writer, rec)?;
                }
                join::Flushed::Merge(beg, n) => {
                    let mut out = self.take_record();
                    let mut rows = std::mem::take(&mut self.joiner.rows);
                    let r = self.merge_biallelics(hdr, out_hdr, &mut out, &mut rows[beg..beg + n]);
                    self.joiner.rows = rows;
                    r?;
                    self.write(out_hdr, writer, out)?;
                }
            }
        }
    }

    /// Whether a record duplicates one already written at its position (the `-d` part of
    /// bcftools' `flush_buffer`); records it otherwise.
    fn is_duplicate(&mut self, hdr: &Header, st: &mut DupState, rec: &mut Record) -> bool {
        let line_type = rec.variant_types();
        let snp_mnp = VariantType::SNP | VariantType::MNP;
        let rmdup = self.opts.rmdup;
        if st.prev == Some((rec.rid(), rec.pos())) {
            if rmdup & PAIR_ANY != 0 {
                return true;
            }
            if rmdup & PAIR_SNPS != 0
                && line_type.intersects(snp_mnp)
                && st.prev_type.intersects(snp_mnp)
            {
                return true;
            }
            if rmdup & PAIR_INDELS != 0
                && line_type.intersects(VariantType::INDEL)
                && st.prev_type.intersects(VariantType::INDEL)
            {
                return true;
            }
            if rmdup & PAIR_EXACT != 0 && self.cmpals_match(hdr, rec) {
                return true;
            }
        } else {
            st.prev = Some((rec.rid(), rec.pos()));
            st.prev_type = VariantType::empty();
            if rmdup & PAIR_EXACT != 0 {
                self.ncmpals = 0;
            }
        }
        st.prev_type |= line_type;
        if rmdup & PAIR_EXACT != 0 {
            self.cmpals_add(hdr, rec);
        }
        false
    }

    /// The ALT allele for duplicate comparison: symbolic alleles get their SVLEN (bcftools'
    /// `strdup_alt_svlen`).
    fn alt_svlen(&mut self, hdr: &Header, rec: &mut Record, ial: usize, out: &mut Vec<u8>) {
        out.clear();
        rec.unpack(Unpack::STR);
        out.extend_from_slice(rec.allele(ial));
        if out.first() != Some(&b'<') {
            return;
        }
        let svlen = &mut self.scratch.i.a;
        let n = rec.get_info_i32(hdr, b"SVLEN", svlen).unwrap_or(0);
        if n == 0 {
            return;
        }
        if n + 1 != rec.allele_count() && !self.svlen_warned {
            eprintln!(
                "TODO: different number of ALT alleles and SVLEN fields {}:{}",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1
            );
            self.svlen_warned = true;
        }
        let v = svlen.get(ial - 1).copied().unwrap_or(0);
        out.push(b'.');
        number::push_int(out, i64::from(v));
    }

    /// bcftools' `cmpals_add`.
    fn cmpals_add(&mut self, hdr: &Header, rec: &mut Record) {
        if self.ncmpals == self.cmpals.len() {
            self.cmpals.push(CmpAls::default());
        }
        rec.unpack(Unpack::STR);
        let mut alt = std::mem::take(&mut self.cmpals[self.ncmpals].alt);
        if rec.allele_count() == 2 {
            self.alt_svlen(hdr, rec, 1, &mut alt);
        }
        let c = &mut self.cmpals[self.ncmpals];
        c.alt = alt;
        c.ref_.clear();
        c.ref_.extend_from_slice(rec.allele(0));
        c.n = rec.allele_count();
        if c.n != 2 {
            c.alts.clear();
            c.alts.extend(rec.alleles().skip(1).map(<[u8]>::to_vec));
        }
        self.ncmpals += 1;
    }

    /// bcftools' `cmpals_match`.
    fn cmpals_match(&mut self, hdr: &Header, rec: &mut Record) -> bool {
        rec.unpack(Unpack::STR);
        let mut alt = std::mem::take(&mut self.scratch.bytes);
        alt.clear();
        if rec.allele_count() == 2 {
            self.alt_svlen(hdr, rec, 1, &mut alt);
        }
        let found = self.cmpals[..self.ncmpals].iter().any(|c| {
            if rec.allele_count() != c.n || !eq_nocase(rec.allele(0), &c.ref_) {
                return false;
            }
            if c.n == 2 {
                return eq_nocase(&alt, &c.alt);
            }
            rec.alleles()
                .skip(1)
                .all(|a| c.alts.iter().any(|b| b.as_slice() == a))
        });
        self.scratch.bytes = alt;
        found
    }

    /// Writes a record (or queues it for the output workers).
    fn write(&mut self, out_hdr: &Header, writer: &mut Writer, rec: Record) -> anyhow::Result<()> {
        if let Some(rec) = self.emitter.write(out_hdr, writer, rec)? {
            self.recycle(rec);
        }
        Ok(())
    }

    /// Starts the `--old-rec-tag` value from a record (bcftools' `old_rec_tag_init`).
    fn old_rec_tag_init(&mut self, hdr: &Header, line: &mut Record) {
        if self.opts.old_rec_tag.is_none() {
            return;
        }
        line.unpack(Unpack::STR);
        let s = &mut self.old_rec;
        s.clear();
        s.extend_from_slice(line.chrom_safe(hdr));
        s.push(b'|');
        number::push_int(s, line.pos() + 1);
        s.push(b'|');
        s.extend_from_slice(line.allele(0));
        s.push(b'|');
        for i in 1..line.allele_count() {
            s.extend_from_slice(line.allele(i));
            if i + 1 < line.allele_count() {
                s.push(b',');
            }
        }
    }

    /// Writes the `--old-rec-tag` value to a record unless it has the tag (bcftools'
    /// `old_rec_tag_set`).
    fn old_rec_tag_set(
        &mut self,
        out_hdr: &Header,
        line: &mut Record,
        ialt: usize,
    ) -> anyhow::Result<()> {
        let Some(tag) = &self.opts.old_rec_tag else {
            return Ok(());
        };
        if self.old_rec.is_empty() {
            return Ok(());
        }
        let id = out_hdr.id(tag).map(|i| i as i32);
        line.unpack(Unpack::INFO);
        if line.info_fields().iter().any(|f| Some(f.key()) == id) {
            return Ok(());
        }
        if ialt > 0 {
            self.old_rec.push(b'|');
            number::push_uint(&mut self.old_rec, ialt as u64);
        }
        if line
            .update_info_string(out_hdr, tag, Some(&self.old_rec))
            .is_err()
        {
            bail!(
                "An error occurred while updating INFO/{}",
                String::from_utf8_lossy(tag)
            );
        }
        self.old_rec.clear();
        Ok(())
    }
}

/// An error while reading (reported after flushing what was read) or a fatal one.
enum ReadOrFatal {
    Read(pasta::variant::synced::SyncedError),
    Fatal(anyhow::Error),
}

impl From<anyhow::Error> for ReadOrFatal {
    fn from(e: anyhow::Error) -> Self {
        Self::Fatal(e)
    }
}

/// HTSlib's `bcf_alleles2gt`: the index of genotype a/b in a Number=G list.
#[inline]
fn alleles2gt(a: usize, b: usize) -> usize {
    if a > b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

/// Whether an allele is a sequence (not symbolic, missing, `*`, or a breakend); bcftools'
/// `is_sequence_allele`.
fn is_sequence_allele(a: &[u8]) -> bool {
    match a.first() {
        None | Some(b'<') => false,
        Some(b'.' | b'*') if a.len() == 1 => false,
        _ => !a.iter().any(|&c| c == b'[' || c == b']'),
    }
}

#[cfg(test)]
mod tests {
    use std::{fs, path::Path};

    use clap::{Args as _, Command, FromArgMatches as _};
    use pasta::variant::io::Reader;

    use super::*;

    /// Runs `norm` with BCF output on `input`; returns the positions read back.
    fn norm_to_bcf(input: &Path, options: &[&str]) -> anyhow::Result<Vec<i64>> {
        let out = input.with_extension(format!("{}.bcf", options.join("")));
        let out_name = out.to_str().expect("UTF-8 path");
        let argv = ["norm", "--no-version", "-Ou", "-o", out_name];
        let input_name = input.to_str().expect("UTF-8 path");
        let matches = Args::augment_args(Command::new("norm"))
            .try_get_matches_from(argv.iter().chain(options).chain([&input_name]))?;
        let ran = run(Args::from_arg_matches(&matches)?, &matches);
        let positions = ran.and_then(|_| {
            let mut reader = Reader::open(out_name, 0)?;
            let mut header = reader.read_header()?;
            let (mut line, mut rec, mut positions) = (Vec::new(), Record::new(), Vec::new());
            while reader.read_record(&mut header, &mut line, &mut rec)? {
                positions.push(rec.pos() + 1);
            }
            Ok(positions)
        });
        if out.exists() {
            fs::remove_file(&out)?;
        }
        positions
    }

    /// Runs `norm` with VCF output on `input`; returns the output after the `##` lines.
    fn norm_to_vcf(input: &Path, options: &[&str]) -> anyhow::Result<String> {
        let out = input.with_extension(format!("{}.out.vcf", options.join("")));
        let out_name = out.to_str().expect("UTF-8 path");
        let argv = ["norm", "--no-version", "-Ov", "-o", out_name];
        let input_name = input.to_str().expect("UTF-8 path");
        let matches = Args::augment_args(Command::new("norm"))
            .try_get_matches_from(argv.iter().chain(options).chain([&input_name]))?;
        let ran = run(Args::from_arg_matches(&matches)?, &matches);
        let text = ran.and_then(|_| Ok(fs::read_to_string(&out)?));
        if out.exists() {
            fs::remove_file(&out)?;
        }
        Ok(text?
            .lines()
            .filter(|l| !l.starts_with("##"))
            .flat_map(|l| [l, "\n"])
            .collect())
    }

    #[test]
    fn split_writes_every_field_kind_to_each_alt() -> anyhow::Result<()> {
        // INFO and FORMAT values of Number A, R, G, and fixed counts, numeric and string, with
        // haploid, diploid, and missing samples in one record. The first record's output is
        // bcftools'. In the second, bcftools does not move past the missing Number=G string of
        // S1, so it leaves the later samples unsplit and drops XG; vtools keeps XG and splits
        // every sample (README).
        let dir = std::env::temp_dir().join(format!("vtools-norm-split-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let input = dir.join("split.vcf");
        let header = "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
            ##INFO=<ID=IA,Number=A,Type=Integer,Description=\"a\">\n\
            ##INFO=<ID=FR,Number=R,Type=Float,Description=\"r\">\n\
            ##INFO=<ID=IG,Number=G,Type=Integer,Description=\"g\">\n\
            ##INFO=<ID=SG,Number=G,Type=String,Description=\"g\">\n\
            ##INFO=<ID=SA,Number=A,Type=String,Description=\"a\">\n\
            ##INFO=<ID=FL,Number=0,Type=Flag,Description=\"f\">\n\
            ##INFO=<ID=I2,Number=2,Type=Integer,Description=\"2\">\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"r\">\n\
            ##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"a\">\n\
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"g\">\n\
            ##FORMAT=<ID=XG,Number=G,Type=String,Description=\"g\">\n\
            ##FORMAT=<ID=XR,Number=R,Type=String,Description=\"r\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n";
        fs::write(
            &input,
            format!(
                "{header}\
                 1\t10\trs1\tA\tC,G,T\t30\tPASS\t\
                 IA=1,2,3;FR=0.5,1.5,2.5,3.5;IG=0,1,2,3,4,5,6,7,8,9;SG=a,b,c,d,e,f,g,h,i,j;SA=x,,z;FL;I2=7,8\t\
                 GT:AD:AF:PL:XG:XR\t\
                 0/2:1,2,3,4:0.1,.,0.3:0,1,2,3,4,5,6,7,8,9:p0,p1,p2,p3,p4,p5,p6,p7,p8,p9:r0,r1,r2,r3\t\
                 3:5,6,.,8:.,0.2,0.3:10,11,12,13:h0,h1,,h3:s0,,s2,s3\t\
                 1|.:.:.:20,21,22,23,24,25,26,27,28,29:.:.\n\
                 1\t20\t.\tA\tC,G\t.\t.\t.\tGT:XG\t0/1:.\t1/2:a,b,c,d,e,f\t2:x,y,z\n"
            ),
        )?;
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n";
        let second = "1\t20\t.\tA\tC\t.\t.\t.\tGT:XG\t0/1:.\t1/0:a,b,c\t0:x,y\n\
                      1\t20\t.\tA\tG\t.\t.\t.\tGT:XG\t0/0:.\t0/1:a,d,f\t1:x,z\n";
        assert_eq!(
            norm_to_vcf(&input, &["-m-"])?,
            format!(
                "{columns}\
                 1\t10\trs1\tA\tC\t30\tPASS\tIA=1;FR=0.5,1.5;IG=0,1,2;SG=a,b,c;SA=x;FL;I2=7,8\t\
                 GT:AD:AF:PL:XG:XR\t0/0:1,2:0.1:0,1,2:p0,p1,p2:r0,r1\t0:5,6:.:10,11:h0,h1:s0,\t\
                 1|.:.,.:.:20,21,22:.:.\n\
                 1\t10\trs1\tA\tG\t30\tPASS\tIA=2;FR=0.5,2.5;IG=0,3,5;SG=a,d,f;SA;FL;I2=7,8\t\
                 GT:AD:AF:PL:XG:XR\t0/1:1,3:.:0,3,5:p0,p3,p5:r0,r2\t0:5,.:0.2:10,12:h0,:s0,s2\t\
                 0|.:.,.:.:20,23,25:.:.\n\
                 1\t10\trs1\tA\tT\t30\tPASS\tIA=3;FR=0.5,3.5;IG=0,6,9;SG=a,g,j;SA=z;FL;I2=7,8\t\
                 GT:AD:AF:PL:XG:XR\t0/0:1,4:0.3:0,6,9:p0,p6,p9:r0,r3\t1:5,8:0.3:10,13:h0,h3:s0,s3\t\
                 0|.:.,.:.:20,26,29:.:.\n{second}"
            )
        );
        assert_eq!(
            norm_to_vcf(&input, &["-m-", "--keep-sum", "AD", "--old-rec-tag", "OLD"])?,
            format!(
                "{columns}\
                 1\t10\trs1\tA\tC\t30\tPASS\tIA=1;FR=0.5,1.5;IG=0,1,2;SG=a,b,c;SA=x;FL;I2=7,8;\
                 OLD=1|10|A|C,G,T|1\tGT:AD:AF:PL:XG:XR\t0/0:8,2:0.1:0,1,2:p0,p1,p2:r0,r1\t\
                 0:13,6:.:10,11:h0,h1:s0,\t1|.:.,.:.:20,21,22:.:.\n\
                 1\t10\trs1\tA\tG\t30\tPASS\tIA=2;FR=0.5,2.5;IG=0,3,5;SG=a,d,f;SA;FL;I2=7,8;\
                 OLD=1|10|A|C,G,T|2\tGT:AD:AF:PL:XG:XR\t0/1:7,3:.:0,3,5:p0,p3,p5:r0,r2\t\
                 0:19,.:0.2:10,12:h0,:s0,s2\t0|.:.,.:.:20,23,25:.:.\n\
                 1\t10\trs1\tA\tT\t30\tPASS\tIA=3;FR=0.5,3.5;IG=0,6,9;SG=a,g,j;SA=z;FL;I2=7,8;\
                 OLD=1|10|A|C,G,T|3\tGT:AD:AF:PL:XG:XR\t0/0:6,4:0.3:0,6,9:p0,p6,p9:r0,r3\t\
                 1:11,8:0.3:10,13:h0,h3:s0,s3\t0|.:.,.:.:20,26,29:.:.\n\
                 1\t20\t.\tA\tC\t.\t.\tOLD=1|20|A|C,G|1\tGT:XG\t0/1:.\t1/0:a,b,c\t0:x,y\n\
                 1\t20\t.\tA\tG\t.\t.\tOLD=1|20|A|C,G|2\tGT:XG\t0/0:.\t0/1:a,d,f\t1:x,z\n"
            )
        );
        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    #[test]
    fn split_time_grows_with_the_record_not_with_alts_times_record() -> anyhow::Result<()> {
        // One indel record with 250 ALTs and full diploid PL (31,626 values per sample): -m-snps
        // reads and writes it whole, and -m-indels splits it. Decoding PL again for each ALT
        // made the split about 45 times slower than the copy; decoded once, it takes about
        // as long. The ratio leaves room for a loaded machine.
        let dir = std::env::temp_dir().join(format!("vtools-norm-cubic-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let input = dir.join("alts.vcf");
        let (nalt, nsmpl) = (250, 8);
        let ng = (nalt + 1) * (nalt + 2) / 2;
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        for s in 0..nsmpl {
            vcf.push_str(&format!("\tS{s}"));
        }
        let alts: Vec<String> = (1..=nalt).map(|i| format!("A{}", "C".repeat(i))).collect();
        vcf.push_str(&format!(
            "\n1\t100\t.\tA\t{}\t.\t.\t.\tGT:PL",
            alts.join(",")
        ));
        let pl: Vec<String> = (0..ng).map(|i| (i % 100).to_string()).collect();
        let pl = pl.join(",");
        for s in 0..nsmpl {
            vcf.push_str(&format!("\t0/{}:{pl}", s + 1));
        }
        vcf.push('\n');
        fs::write(&input, vcf)?;
        let time = |option: &str| -> anyhow::Result<std::time::Duration> {
            let start = std::time::Instant::now();
            norm_to_bcf(&input, &[option])?;
            Ok(start.elapsed())
        };
        // The faster of two runs each, against noise from other processes.
        let copy = time("-m-snps")?.min(time("-m-snps")?);
        let split = time("-m-indels")?.min(time("-m-indels")?);
        fs::remove_dir_all(&dir)?;
        assert!(
            split < copy * 10,
            "splitting took {split:?}, copying the record {copy:?}"
        );
        Ok(())
    }

    #[test]
    fn bcf_output_holds_every_position_below_2_32() -> anyhow::Result<()> {
        // bcftools buffers bcf_dup copies, which lose the VCF parser's 64-bit flag, so its BCF
        // output holds the records up to POS 2^32 - 1 (where BCF's unsigned POS ends), split
        // or not. From 2^32 on, bcftools truncates POS and vtools fails as bcf_write does.
        let dir = std::env::temp_dir().join(format!("vtools-norm-pos64-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        for pos in [
            2_147_483_648_i64,
            3_000_000_000,
            4_294_967_295,
            4_294_967_296,
        ] {
            let input = dir.join(format!("{pos}.vcf"));
            fs::write(
                &input,
                format!(
                    "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     1\t{pos}\t.\tA\tC,G\t.\t.\t.\n1\t{pos}\t.\tA\tT\t.\t.\t.\n"
                ),
            )?;
            for (options, records) in [(&["-m-"][..], 3), (&["-m+"], 1), (&["-d", "exact"], 2)] {
                let positions = norm_to_bcf(&input, options);
                if pos < 1 << 32 {
                    assert_eq!(positions?, vec![pos; records], "POS {pos} {options:?}");
                } else {
                    assert!(positions.is_err(), "POS {pos} {options:?}");
                }
            }
        }
        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    #[test]
    fn a_panicking_parser_worker_ends_the_run() -> anyhow::Result<()> {
        // One worker panics on its batch of lines while the others wait for jobs, holding the
        // results channel open: the main thread must stop waiting and re-raise the panic.
        let dir = crate::raw_input::TestDir::new("norm-panic");
        let input = dir.join("panic.vcf");
        fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 1\t1\t{}\tA\tG,T\t.\t.\t.\n",
                crate::workers::PANIC_MARKER
            ),
        )?;
        let panicked =
            crate::raw_input::within_seconds("norm with a panicking worker", move || {
                std::panic::catch_unwind(|| norm_to_vcf(&input, &["-m-", "--threads", "3"]))
                    .is_err()
            });
        assert!(panicked, "the worker's panic is re-raised");
        Ok(())
    }

    #[test]
    fn records_with_fewer_samples_than_the_header_fail_without_a_panic() -> anyhow::Result<()> {
        // BCF records of three samples under a header of four. Splitting the Number=G string
        // read past the end of the record's FORMAT data and panicked; now the split records
        // fail to write, as bcf_write fails on such records. (bcftools reads past the record
        // too, and fails on whatever it finds there.)
        let text = "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
            ##FORMAT=<ID=XG,Number=G,Type=String,Description=\"g\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3";
        let mut three = Header::parse(format!("{text}\n").as_bytes())?;
        let four = Header::parse(format!("{text}\tS4\n").as_bytes())?;
        let mut data = pasta::variant::bcf::encode_header(&four);
        let mut line =
            b"1\t10\t.\tA\tC,G\t.\t.\t.\tGT:XG\t0/1:a,b,c,d,e,f\t1:x,y,z\t2:p,q,r".to_vec();
        let mut rec = Record::new();
        pasta::variant::vcf::parse(&mut three, &mut line, &mut rec)?;
        pasta::variant::bcf::encode_record(&three, &mut rec, &mut data)?;
        let dir = crate::raw_input::TestDir::new("norm-fewer-samples");
        let input = dir.join("fewer.bcf");
        fs::write(&input, data)?;
        for threads in ["0", "2"] {
            let options = ["-m-", "--threads", threads];
            let e = norm_to_vcf(&input, &options).expect_err("the split records are not written");
            let out = input.with_extension(format!("{}.out.vcf", options.join("")));
            let want = format!("[flush_buffer] Error: cannot write to {}", out.display());
            assert_eq!(e.to_string(), want);
        }
        Ok(())
    }

    #[test]
    fn number_g_fields_of_any_ploidy_join_and_split() -> anyhow::Result<()> {
        // Haploid to tetraploid samples next to each other, with missing values. bcftools
        // stops on the joins ("could not merge FORMAT tag") or writes the values of every
        // ploidy as haploid ones, and refuses the split ("wrong number of fields") (README).
        let dir = std::env::temp_dir().join(format!("vtools-norm-ploidy-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let input = dir.join("ploidy.vcf");
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n";
        fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
                 ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"g\">\n\
                 ##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"g\">\n\
                 ##FORMAT=<ID=XG,Number=G,Type=String,Description=\"g\">\n{columns}\
                 1\t10\t.\tA\tC\t.\t.\t.\tPL:GP\t0,1,2,3,4:.\t.:.\t5,6:0.5,0.6\t7,8,9,10:.\n\
                 1\t10\t.\tA\tG\t.\t.\t.\tPL\t10,11,12,13,14\t20,21,22\t.\t.\n\
                 1\t10\t.\tA\tT\t.\t.\t.\tGP\t1,2,3,4,5\t.\t.\t0.1,0.2,0.3,0.4\n\
                 1\t20\t.\tA\tC,G\t.\t.\t.\tPL\t0,1,2,3,4,5,6,7,8,9,10,11,12,13,14\t.\t5,6,7\t1,2,3,4,5,6\n\
                 1\t30\t.\tA\tC\t.\t.\t.\tXG\ta,b,c,d\te,f,g\th,i\tj,k,l,m,n\n\
                 1\t30\t.\tA\tG\t.\t.\t.\tXG\tp,q,r,s\tt,u,v\tw,x\t.\n\
                 1\t40\t.\tA\tC,G\t.\t.\t.\tXG\ta,b,c,d,e,f,g,h,i,j\t.\tx,y,z\tp,q,r,s,t,u\n"
            ),
        )?;
        // Tetraploid with four alleles: 35 values, of which line A,G gives 0000, 0002, 0022,
        // 0222, and 2222 (0, 5, 9, 12, and 14) and line A,T the ones of allele 3.
        let m = |n: usize| vec!["."; n].join(",");
        let joined = format!(
            "{columns}\
             1\t10\t.\tA\tC,G,T\t.\t.\t.\tPL:GP\t\
             10,1,2,3,4,11,.,.,.,12,.,.,13,.,14,{}:1,{},2,{},3,{},4,.,.,5\t\
             20,.,.,21,.,22,{}:.\t5,6,.,.:0.5,0.6,.,.\t\
             7,8,9,10,{}:0.1,{},0.2,{},0.3,.,.,0.4\n\
             1\t20\t.\tA\tC,G\t.\t.\t.\tPL\t0,1,2,3,4,5,6,7,8,9,10,11,12,13,14\t.\t5,6,7\t1,2,3,4,5,6\n\
             1\t30\t.\tA\tC,G\t.\t.\t.\tXG\ta,b,c,d,q,.,.,r,.,s\te,f,g,u,.,v\th,i,x\tj,k,l,m,n,{}\n\
             1\t40\t.\tA\tC,G\t.\t.\t.\tXG\ta,b,c,d,e,f,g,h,i,j\t.\tx,y,z\tp,q,r,s,t,u\n",
            m(20),
            m(14),
            m(9),
            m(5),
            m(4),
            m(16),
            m(9),
            m(5),
            m(10),
        );
        let split = format!(
            "{columns}\
             1\t10\t.\tA\tC\t.\t.\t.\tPL:GP\t0,1,2,3,4:.\t.:.\t5,6:0.5,0.6\t7,8,9,10:.\n\
             1\t10\t.\tA\tG\t.\t.\t.\tPL\t10,11,12,13,14\t20,21,22\t.\t.\n\
             1\t10\t.\tA\tT\t.\t.\t.\tGP\t1,2,3,4,5\t.\t.\t0.1,0.2,0.3,0.4\n\
             1\t20\t.\tA\tC\t.\t.\t.\tPL\t0,1,2,3,4\t.\t5,6\t1,2,3\n\
             1\t20\t.\tA\tG\t.\t.\t.\tPL\t0,5,9,12,14\t.\t5,7\t1,4,6\n\
             1\t30\t.\tA\tC\t.\t.\t.\tXG\ta,b,c,d\te,f,g\th,i\tj,k,l,m,n\n\
             1\t30\t.\tA\tG\t.\t.\t.\tXG\tp,q,r,s\tt,u,v\tw,x\t.\n\
             1\t40\t.\tA\tC\t.\t.\t.\tXG\ta,b,c,d\t.\tx,y\tp,q,r\n\
             1\t40\t.\tA\tG\t.\t.\t.\tXG\ta,e,h,j\t.\tx,z\tp,s,u\n"
        );
        let bcf = dir.join("ploidy.bcf");
        let argv = [
            "view",
            "--no-version",
            "-Ob",
            "-o",
            bcf.to_str().expect("UTF-8 path"),
        ];
        let input_name = input.to_str().expect("UTF-8 path");
        let matches = crate::commands::view::Args::augment_args(Command::new("view"))
            .try_get_matches_from(argv.iter().chain([&input_name]))?;
        crate::commands::view::run(
            crate::commands::view::Args::from_arg_matches(&matches)?,
            &matches,
        )?;
        for file in [&input, &bcf] {
            assert_eq!(norm_to_vcf(file, &["-m+"])?, joined, "{file:?}");
            assert_eq!(norm_to_vcf(file, &["-m-"])?, split, "{file:?}");
        }
        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    #[test]
    fn polyploid_values_that_cannot_be_joined_stop() -> anyhow::Result<()> {
        let dir =
            std::env::temp_dir().join(format!("vtools-norm-ploidy-bad-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let header = "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
                      ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"g\">\n\
                      ##FORMAT=<ID=XG,Number=G,Type=String,Description=\"g\">\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
        // Ploidy 150 over six alleles: 7.5e8 values a sample, past the 1 GiB limit.
        let wide: Vec<String> = (0..=150).map(|v| v.to_string()).collect();
        let wide = wide.join(",");
        let wide = ["C", "G", "T", "N", "M"]
            .map(|a| format!("1\t10\t.\tA\t{a}\t.\t.\t.\tPL\t{wide}\t.\n"))
            .concat();
        // Lines A,C and A,G of one field with the values of S1 and S2: a sample whose ploidy
        // changes between the lines of a field with a triploid sample, in numbers and in
        // strings. From 3 to 1, from 1 or 2 to 3, from 2 to 1 next to a triploid sample, and
        // from 2 to 3 where the first line has no triploid sample.
        let pairs = [
            ("pl", "PL", ["1,2,3,4", "1,2"], ["5,6", "."]),
            ("xg", "XG", ["a,b,c,d", "x,y"], ["A,B", "X,Y"]),
            ("pl_dip", "PL", ["1,2,3", "1,2,3,4"], ["5,6,7,8", "5,6,7,8"]),
            ("xg_dip", "XG", ["a,b,c", "a,b,c,d"], ["p,q,r,s", "p,q,r,s"]),
            ("pl_hap", "PL", ["1,2", "1,2,3,4"], ["5,6,7,8", "5,6,7,8"]),
            ("xg_hap", "XG", ["a,b", "a,b,c,d"], ["p,q,r,s", "p,q,r,s"]),
            ("pl_next", "PL", ["1,2,3,4", "1,2,3"], ["5,6,7,8", "5,6"]),
            ("xg_next", "XG", ["a,b,c,d", "a,b,c"], ["p,q,r,s", "p,q"]),
            ("pl_later", "PL", ["1,2,3", "1,2,3"], ["5,6,7,8", "5,6,7"]),
            ("xg_later", "XG", ["a,b,c", "a,b,c"], ["p,q,r,s", "p,q,r"]),
        ];
        let pairs = pairs.map(|(name, tag, c, g)| {
            let records = format!(
                "1\t10\t.\tA\tC\t.\t.\t.\t{tag}\t{}\t{}\n\
                 1\t10\t.\tA\tG\t.\t.\t.\t{tag}\t{}\t{}\n",
                c[0], c[1], g[0], g[1]
            );
            (name, records)
        });
        let cases = std::iter::once(("wide", wide)).chain(pairs);
        for (name, records) in cases {
            let input = dir.join(format!("{name}.vcf"));
            fs::write(&input, format!("{header}{records}"))?;
            let err = norm_to_vcf(&input, &["-m+"]).expect_err(name);
            assert!(
                err.to_string().contains("could not merge FORMAT tag"),
                "{name}: {err}"
            );
        }
        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    #[test]
    fn no_realign_is_read_as_strtod_reads_it() {
        for (value, limit) in [
            ("", 0),
            ("5", 5),
            ("1e1", 10),
            ("0.5", 0),
            ("-0", 0),
            ("nan", 0),
        ] {
            assert_eq!(parse_no_realign(value).ok(), Some(limit), "value {value:?}");
        }
        for value in ["x", "=5", "-1", "inf", "4294967296", "5x"] {
            assert!(parse_no_realign(value).is_err(), "value {value:?}");
        }
    }
}
