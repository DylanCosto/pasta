//! `vtools concat`: concatenate or combine VCF/BCF files with the same samples (bcftools
//! concat).
//!
//! Plain concatenation streams the files one after another, checking that each chromosome is
//! one sorted block. `-a` merges overlapping files position by position through the synced
//! reader, `-l` ligates overlapping phased chunks, and `-n` copies the compressed blocks
//! without decompressing the records.

mod ligate;
mod naive;

use std::{cell::Cell, path::PathBuf, process::ExitCode, time::Instant};

use anyhow::bail;
use clap::ArgMatches;
use pasta::util::inputs;
use pasta::variant::{
    Header, HeaderRecordType, Record,
    index::split_idx_name,
    io::{
        Compression as ReadCompression, InputFormat, OutputType, ReadError, Reader, WriteError,
        Writer, reader::is_stream,
    },
    log,
    readlist::{self, ListSource},
    synced::{IndexRequirement, Overlap, Pairing, Regions, SyncedReaders},
};

use crate::options::last_index;
use crate::{
    commands::norm::emit::Emitter,
    output::{self, OutputSetting},
    raw_input::{
        Compression, Format, FrameCopier, InflatingReader, RawInput, report_exec_format_error,
        strerror,
    },
    version,
};

const USAGE: &str = "
About:   Concatenate or combine VCF/BCF files. All source files must have the same sample
         columns appearing in the same order. The program can be used, for example, to
         concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel
         VCF into one. The input files must be sorted by chr and position. The files
         must be given in the correct order to produce sorted VCF on output unless
         the -a, --allow-overlaps option is specified. With the --naive option, the files
         are concatenated without being recompressed, which is very fast.
Usage:   vtools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]

Options:
   -a, --allow-overlaps           First coordinate of the next file can precede last record of the current file.
   -c, --compact-PS               Do not output PS tag at each site, only at the start of a new phase set block.
   -d, --rm-dups STRING           Output duplicate records present in multiple files only once: <snps|indels|both|all|exact>
   -D, --remove-duplicates        Alias for -d exact
   -f, --file-list FILE           Read the list of files from a file.
   -G, --drop-genotypes           Drop individual genotype information.
   -l, --ligate                   Ligate phased VCFs by matching phase at overlapping haplotypes
       --ligate-force             Ligate even non-overlapping chunks, keep all sites
       --ligate-warn              Drop sites in imperfect overlaps
       --no-version               Do not append version and command line to the header
   -n, --naive                    Concatenate files without recompression, a header check compatibility is performed
       --naive-force              Same as --naive, but header compatibility is not checked. Dangerous, use with caution.
   -o, --output FILE              Write output to a file [standard output]
   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
   -q, --min-PQ INT               Break phase set if phasing quality is lower than <int> [30]
   -r, --regions REGION           Restrict to comma-separated list of regions
   -R, --regions-file FILE        Restrict to regions listed in a file
       --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
       --threads INT              Use multithreading with <int> worker threads [0]
   -v, --verbosity INT            Set verbosity level
   -W, --write-index[=FMT]        Automatically index the output files [off]

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF files.
    inputs: Vec<String>,

    /// First coordinate of the next file can precede last record of the current file.
    #[arg(short = 'a', long = "allow-overlaps")]
    allow_overlaps: bool,

    /// Do not output PS tag at each site, only at the start of a new phase set block.
    #[arg(short = 'c', long = "compact-PS")]
    compact_ps: bool,

    /// Output duplicate records present in multiple files only once:
    /// snps|indels|both|all|exact.
    #[arg(
        short = 'd',
        long = "rm-dups",
        value_name = "STRING",
        allow_hyphen_values = true
    )]
    rm_dups: Option<String>,

    /// Alias for -d exact.
    #[arg(short = 'D', long = "remove-duplicates")]
    remove_duplicates: bool,

    /// Read the list of files from a file.
    #[arg(
        short = 'f',
        long = "file-list",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    file_list: Option<String>,

    /// Drop individual genotype information.
    #[arg(short = 'G', long = "drop-genotypes")]
    drop_genotypes: bool,

    /// Ligate phased VCFs by matching phase at overlapping haplotypes.
    #[arg(short = 'l', long = "ligate")]
    ligate: bool,

    /// Ligate even non-overlapping chunks, keep all sites.
    #[arg(long = "ligate-force")]
    ligate_force: bool,

    /// Drop sites in imperfect overlaps.
    #[arg(long = "ligate-warn")]
    ligate_warn: bool,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Concatenate files without recompression, after checking that the headers are
    /// compatible.
    #[arg(short = 'n', long = "naive")]
    naive: bool,

    /// Same as --naive, but header compatibility is not checked.
    #[arg(long = "naive-force")]
    naive_force: bool,

    /// Write output to a file [standard output].
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "u|b|v|z[0-9]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output_type: Vec<String>,

    /// Break phase set if phasing quality is lower than INT [30].
    #[arg(
        short = 'q',
        long = "min-PQ",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    min_pq: Option<String>,

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

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2).
    #[arg(
        long = "regions-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true
    )]
    regions_overlap: Option<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Set verbosity level (0 or 1).
    #[arg(
        short = 'v',
        long = "verbosity",
        visible_alias = "verbose",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    verbosity: Option<String>,

    /// Automatically index the output files [off]: csi (default) or tbi.
    #[arg(
        short = 'W',
        long = "write-index",
        value_name = "FMT",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = ""
    )]
    write_index: Option<String>,

    /// Print help (`-h` takes an argument in bcftools and prints the usage).
    #[arg(
        short = 'h',
        long = "help",
        num_args = 0..=1,
        default_missing_value = ""
    )]
    help: Option<String>,

    /// Print the usage (`-?`, which takes no argument).
    #[arg(short = '?', hide = true, action = clap::ArgAction::Count)]
    question: u8,
}

/// The state of bcftools concat (`args_t`).
pub struct Concat {
    fnames: Vec<String>,
    out_hdr: Header,
    seen_seq: Vec<bool>,
    sites_only: bool,
    verbose: bool,
    threads: usize,
    output_name: String,
    /// The output type the writer was opened with.
    out_type: OutputType,
    has_index: bool,
    /// How the last record was written, to word a failure that the compression threads
    /// report only when the output is closed.
    last_write: Cell<LastWrite>,
}

/// How a record was written: bcftools words a failed write by the way it wrote.
#[derive(Clone, Copy, Debug, Default)]
enum LastWrite {
    /// Nothing was written.
    #[default]
    Nothing,
    /// A VCF line of this many bytes, copied as it is.
    Line(usize),
    /// A record of a plain concatenation.
    Record,
    /// A record written by this bcftools function (`-a` and `-l`).
    Named(&'static str),
}

impl LastWrite {
    /// bcftools' error for a failed write of this kind to `output_name`.
    fn failure(self, output_name: &str) -> anyhow::Error {
        match self {
            // Only a written block can fail before closing.
            Self::Nothing => anyhow::anyhow!("Error: close failed .. {output_name}"),
            Self::Line(n) => anyhow::anyhow!("\nFailed to write {n} bytes"),
            Self::Record => anyhow::anyhow!("\nFailed to write"),
            Self::Named(func) => {
                anyhow::anyhow!("[{func}] Error: cannot write to {output_name}")
            }
        }
    }
}

/// An error that bcftools reports with `exit(1)` rather than `error()`.
#[derive(Debug)]
pub struct ExitOne;

impl std::fmt::Display for ExitOne {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("exit 1")
    }
}

impl std::error::Error for ExitOne {}

/// Collects `-O` in command-line order.
fn output_settings(args: &Args) -> Vec<OutputSetting> {
    args.output_type
        .iter()
        .map(|t| OutputSetting::Type(t.clone()))
        .collect()
}

/// C's `strtol(s, &end, 10)` followed by bcftools' `*end` check: the whole string must parse,
/// and an empty one is 0.
fn strtol_whole(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    let (v, end) = pasta::variant::number::strtol(b, 0);
    (end == b.len()).then_some(v)
}

/// C's `strtol(s, &end, 0)` (decimal, `0x` hexadecimal, or `0` octal) followed by bcftools'
/// `*end` check: the whole string must parse, and an empty one is 0.
fn strtol_base0_whole(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    if b.is_empty() {
        return Some(0);
    }
    let mut i = b
        .iter()
        .take_while(|c| matches!(c, b' ' | b'\t'..=b'\r'))
        .count();
    let neg = match b.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let hex = b.get(i) == Some(&b'0')
        && matches!(b.get(i + 1), Some(b'x' | b'X'))
        && b.get(i + 2).is_some_and(u8::is_ascii_hexdigit);
    let (radix, start) = match (hex, b.get(i)) {
        (true, _) => (16, i + 2),
        (false, Some(b'0')) => (8, i),
        _ => (10, i),
    };
    let digits = &b[start.min(b.len())..];
    let n = digits
        .iter()
        .take_while(|&&c| char::from(c).is_digit(radix))
        .count();
    if n == 0 || start + n != b.len() {
        return None;
    }
    let magnitude = digits.iter().try_fold(0i64, |acc, &c| {
        let d = i64::from(char::from(c).to_digit(radix)?);
        acc.checked_mul(i64::from(radix))?.checked_add(d)
    });
    // An out-of-range value saturates in C; any such value is rejected by the callers' ranges.
    let magnitude = magnitude.unwrap_or(i64::MAX);
    Some(if neg { -magnitude } else { magnitude })
}

/// `-d/-D`: which of the records at one position `-a` writes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Duplicates {
    /// All of them.
    Keep,
    /// Only the first file's.
    Remove,
}

/// The options parsed as bcftools' option loop parses them.
struct Parsed {
    verbose: bool,
    min_pq: i32,
    /// `--threads` as given, negative when no thread pool can be made.
    n_threads: i64,
    threads: usize,
    regions_overlap: Overlap,
    write_index: Option<output::WriteIndex>,
    settings: Vec<OutputSetting>,
    /// `-d` (or `exact` for `-D`), whichever was given last.
    remove_dups: Option<String>,
    /// `-r` or `-R`, whichever was given last.
    regions: Option<(String, ListSource)>,
}

/// Parses the option values, in bcftools' order of checks.
fn parse_options(args: &Args, matches: &ArgMatches) -> anyhow::Result<Parsed> {
    let mut verbose = true;
    if let Some(v) = &args.verbosity {
        match strtol_base0_whole(v) {
            Some(n) if (0..=1).contains(&n) => verbose = n == 1,
            _ => bail!("Error: currently only --verbose 0 or --verbose 1 is supported"),
        }
    }
    let min_pq = match &args.min_pq {
        Some(q) => match strtol_whole(q) {
            Some(v) => v as i32,
            None => bail!("Could not parse argument: --min-PQ {q}"),
        },
        None => 30,
    };
    // `int n_threads = strtol(optarg, 0, 0)`.
    let n_threads = args
        .threads
        .as_deref()
        .map_or(0, |t| i64::from(output::thread_count(t)));
    let threads = usize::try_from(n_threads).unwrap_or(0);
    let regions_overlap = match &args.regions_overlap {
        Some(s) => match Overlap::parse(s) {
            Some(o) => o,
            None => bail!("Could not parse: --regions-overlap {s}"),
        },
        None => Overlap::Pos,
    };
    let write_index = args
        .write_index
        .as_deref()
        .map(output::parse_write_index)
        .transpose()?;
    let settings = output_settings(args);
    // -O is checked as it is parsed.
    output::resolve(&settings, None)?;
    // -d and -D: the last one wins.
    let remove_dups = match (&args.rm_dups, args.remove_duplicates) {
        (Some(d), true) => {
            if last_index(matches, "rm_dups") > last_index(matches, "remove_duplicates") {
                Some(d.clone())
            } else {
                Some("exact".to_owned())
            }
        }
        (Some(d), false) => Some(d.clone()),
        (None, true) => Some("exact".to_owned()),
        (None, false) => None,
    };
    // -r and -R: the last one wins (-R also marks the list as a file for a later -r).
    let regions = match (&args.regions, &args.regions_file) {
        (Some(r), Some(f)) => Some(
            if last_index(matches, "regions") > last_index(matches, "regions_file") {
                (r.clone(), ListSource::File)
            } else {
                (f.clone(), ListSource::File)
            },
        ),
        (Some(r), None) => Some((r.clone(), ListSource::Inline)),
        (None, Some(f)) => Some((f.clone(), ListSource::File)),
        (None, None) => None,
    };
    Ok(Parsed {
        verbose,
        min_pq,
        n_threads,
        threads,
        regions_overlap,
        write_index,
        settings,
        remove_dups,
        regions,
    })
}

/// The options that cannot be combined.
fn check_combinations(args: &Args) -> anyhow::Result<()> {
    if args.ligate_force && args.ligate_warn {
        bail!("The options cannot be combined: --ligate-force and --ligate-warn");
    }
    if args.allow_overlaps && args.ligate {
        bail!("The options -a and -l should not be combined. Please run with -l only.");
    }
    if args.drop_genotypes && args.ligate {
        bail!("The options --drop-genotypes and --ligate cannot be combined");
    }
    if args.compact_ps && !args.ligate {
        bail!("The -c option is intended only with -l");
    }
    Ok(())
}

/// The input file names: the arguments, or the `-f` list.
fn input_names(args: &mut Args) -> anyhow::Result<Vec<String>> {
    let mut fnames = std::mem::take(&mut args.inputs);
    if let Some(list) = &args.file_list {
        if !fnames.is_empty() {
            bail!("Cannot combine -l with file names on command line.");
        }
        // hts_readlines: a file, or a `:a,b` list when no such file opens.
        match readlist::read_lines_or_inline(list.as_bytes()) {
            Ok(lines) => {
                fnames = lines
                    .into_iter()
                    .map(|l| String::from_utf8_lossy(&l).into_owned())
                    .collect();
            }
            Err(_) => bail!("Could not read the file: {list}"),
        }
    }
    Ok(fnames)
}

pub fn run(mut args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    if args.help.is_some() || args.question > 0 {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    let opts = parse_options(&args, matches)?;
    check_combinations(&args)?;
    let fnames = input_names(&mut args)?;
    if fnames.is_empty() {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    // The inputs are closed when the output is opened, and read again after.
    if let Some(out) = &args.output {
        let data: Vec<&str> = fnames.iter().map(|f| split_idx_name(f).0).collect();
        if let Err(e) = inputs::check_output_against(&output::data_path(out), &data) {
            bail!("Can't write to \"{}\": {e}", out.display());
        }
    }
    if opts.remove_dups.is_some() && !args.allow_overlaps {
        bail!("The -D option is supported only with -a");
    }
    if opts.regions.is_some() && !args.allow_overlaps {
        bail!("The -r/-R option is supported only with -a");
    }
    if args.naive || args.naive_force {
        run_naive(&args, &fnames, &opts)?;
        return Ok(ExitCode::SUCCESS);
    }
    let output_name = args
        .output
        .as_deref()
        .map_or_else(|| "-".to_owned(), |p| p.display().to_string());

    if opts.verbose {
        eprintln!(
            "Checking the headers and starting positions of {} files",
            fnames.len()
        );
    }
    let (out_hdr, start_pos) = init_header(&args, &fnames)?;
    let seen_seq = vec![false; out_hdr.contig_count()];
    let (out_type, level) = output::resolve(&opts.settings, args.output.as_deref())?;
    let out_path = args.output.as_deref().map(output::data_path);
    // bcftools opens the output and then starts its thread pool. The output's compression
    // threads start with it here, and may take the room the check needs, so the check comes
    // first and fails after the output's error.
    let threads_start = output::threads_start(opts.threads);
    let mut writer = match Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, opts.threads),
    ) {
        Ok(w) => w,
        Err(e) => {
            output::log_open_failure(&output_name, &e);
            bail!("Can't write to \"{output_name}\": {}", strerror(&e))
        }
    };
    let mut files = (args.allow_overlaps || args.ligate).then(|| {
        let mut sr = SyncedReaders::new();
        sr.set_index_requirement(IndexRequirement::Required);
        // `-a` reads every file in step; `-l` reads two chunks at a time.
        let in_step = if args.ligate {
            fnames.len().min(2)
        } else {
            fnames.len()
        };
        sr.set_shared_threads(opts.threads, in_step);
        sr
    });
    let n_threads = opts.n_threads;
    if n_threads < 0 || !threads_start {
        writer.abandon();
        if files.is_some() {
            bail!("Failed to create threads");
        }
        bail!("Failed to initialize {n_threads} threads");
    }
    if writer.write_header(&out_hdr).is_err() {
        writer.abandon();
        bail!("[init_data] Error: cannot write the header to {output_name}");
    }
    let mut index_name = String::new();
    if let Some(fmt) = opts.write_index {
        let ok = match output::index_request(fmt, out_type, args.output.as_deref()) {
            Some(req) => {
                index_name = req.path.display().to_string();
                writer.request_index(req).is_ok()
            }
            None => false,
        };
        if !ok {
            writer.abandon();
            bail!("Error: failed to initialise index for {output_name}");
        }
        writer.start_index(&out_hdr)?;
    }
    let mut state = Concat {
        fnames,
        out_hdr,
        seen_seq,
        sites_only: args.drop_genotypes,
        verbose: opts.verbose,
        threads: opts.threads,
        output_name,
        out_type,
        has_index: opts.write_index.is_some(),
        last_write: Cell::default(),
    };
    // A failure to open the inputs abandons the output too, rather than dropping it: with
    // threads, a dropped output is closed by the compression threads.
    let result = if let Some(files) = files.as_mut() {
        if args.allow_overlaps {
            let remove_dups = opts.remove_dups.as_deref();
            let duplicates = if remove_dups.is_some() {
                Duplicates::Remove
            } else {
                Duplicates::Keep
            };
            state
                .init_overlaps(files, opts.regions, opts.regions_overlap, remove_dups)
                .and_then(|()| state.concat_overlaps(files, &mut writer, duplicates))
        } else {
            ligate::Ligate::new(
                &mut state,
                start_pos,
                ligate::Options {
                    compact_ps: args.compact_ps,
                    min_pq: opts.min_pq,
                    force: args.ligate_force,
                    warn: args.ligate_warn,
                },
            )
            .and_then(|mut lig| lig.run(&mut state, files, &mut writer))
        }
    } else {
        state.concat_plain(&mut writer)
    };
    state.finish(writer, result, &index_name)
}

/// `-n` and `--naive-force`: the files' compressed blocks copied one after another.
fn run_naive(args: &Args, fnames: &[String], opts: &Parsed) -> anyhow::Result<()> {
    let mode = if args.naive_force {
        "--naive-force"
    } else {
        "--naive"
    };
    if opts.write_index.is_some() {
        bail!("Error: cannot --write-index in the {mode} mode");
    }
    if args.allow_overlaps {
        bail!("The option --naive cannot be combined with --allow-overlaps");
    }
    if args.ligate {
        bail!("The option --naive cannot be combined with --ligate");
    }
    if args.drop_genotypes {
        bail!("The option --naive cannot be combined with --drop-genotypes");
    }
    naive::naive_concat(&naive::Options {
        fnames,
        output: args.output.as_deref(),
        settings: &opts.settings,
        explicit_output_type: !args.output_type.is_empty(),
        check_headers: !args.naive_force,
        verbose: opts.verbose,
    })
}

/// bcftools' `init_data`: the output header, merged from every file's, and with `-l` the
/// position of each file's first record (-1 on a new chromosome, -2 for no records).
fn init_header(args: &Args, fnames: &[String]) -> anyhow::Result<(Header, Vec<i64>)> {
    let mut out_hdr: Option<Header> = None;
    let mut start_pos: Vec<i64> = Vec::new();
    let mut prev_chrid: Option<u32> = None;
    let mut line = Vec::new();
    let mut rec = Record::new();
    for f in fnames {
        let (mut reader, hdr) = open_input(f, 0, Reading::First)?;
        let mut hdr = if args.drop_genotypes {
            drop_genotypes(&hdr)?
        } else {
            hdr
        };
        let merged = match out_hdr.take() {
            None => Header::merge_new(&hdr)?,
            Some(mut h) => {
                h.merge(&hdr)?;
                h
            }
        };
        if hdr.sample_count() != merged.sample_count() {
            bail!(
                "Different number of samples in {f}. Perhaps \"vtools merge\" is what you are looking for?"
            );
        }
        if hdr.samples() != merged.samples() {
            bail!(
                "Different sample names in {f}. Perhaps \"vtools merge\" is what you are looking for?"
            );
        }
        if args.ligate {
            match reader.read_record(&mut hdr, &mut line, &mut rec) {
                Ok(true) => {
                    let chrid = merged.contig_id(rec.chrom_safe(&hdr));
                    start_pos.push(if chrid == prev_chrid { rec.pos() } else { -1 });
                    prev_chrid = chrid;
                }
                // A record that fails to parse ends the chunk, as it ends plain concatenation
                // and as HTSlib's reader takes it.
                Ok(false) | Err(ReadError::Parse(_)) => start_pos.push(-2),
                // A damaged or truncated chunk must not be dropped from the output quietly.
                Err(_) => return Err(bcf_read_error(f)),
            }
        }
        out_hdr = Some(merged);
    }
    let mut out_hdr = out_hdr.expect("there is at least one file");
    if args.ligate {
        out_hdr.append(
            b"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phasing Quality (bigger is better)\">",
        )?;
        out_hdr.append(b"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">")?;
    }
    if !args.no_version {
        version::append_version(&mut out_hdr, "concat")?;
    }
    Ok((out_hdr, start_pos))
}

/// Whether [`open_input`] reads an input for the first time.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Reading {
    /// The pass over the headers ([`init_header`]).
    First,
    /// The pass over the records, after the headers.
    Again,
}

/// Opens an input and reads its header, with bcftools' messages: `hts_open` fails on an I/O
/// error or on binary data it does not know, and `bcf_hdr_read` on other data that is neither
/// VCF nor BCF. Read again, a stream goes on from where the first reading left it, as in
/// bcftools, without waiting for the writer of a FIFO (bcftools waits forever).
fn open_input(f: &str, threads: usize, reading: Reading) -> anyhow::Result<(Reader, Header)> {
    // The pass over the records follows "Concatenating FILE" on the same line.
    let (lead, opened) = match reading {
        Reading::First => ("", Reader::open(f, threads)),
        Reading::Again => ("\n", Reader::open_again(f, threads)),
    };
    let mut reader = match opened {
        Ok(reader) => reader,
        Err(ReadError::UnknownBinary) => {
            report_exec_format_error(f);
            bail!("{lead}Failed to open: {f}");
        }
        Err(ReadError::UnknownFormat) => {
            log::error(
                "bcf_hdr_read",
                format_args!("Input is not detected as bcf or vcf format"),
            );
            bail!("{lead}Failed to parse header: {f}");
        }
        Err(ReadError::Io(e)) => {
            output::log_open_failure(f, &e);
            bail!("{lead}Failed to open: {f}");
        }
        Err(_) => bail!("{lead}Failed to open: {f}"),
    };
    let Ok(hdr) = reader.read_header() else {
        bail!("{lead}Failed to parse header: {f}");
    };
    Ok((reader, hdr))
}

/// The error for a truncated or malformed BCF record, or a truncated or corrupt BGZF block (of
/// BCF or VCF), in plain concatenation. bcftools' `bcf_read` loop takes it for the end of the
/// file and goes on with the next one, so the records after it are lost with status 0.
fn bcf_read_error(f: &str) -> anyhow::Error {
    output::InputReadFailure(format!("\nError: BCF read error: {f}")).into()
}

/// HTSlib's warning when reading reaches the end of a BGZF file whose last block is not the
/// EOF marker, as a file cut at a block boundary is. (A pipe is checked as it is read.)
fn warn_if_no_eof(reader: &Reader, f: &str) {
    if reader.compression() == ReadCompression::Bgzf {
        super::index::warn_if_no_eof(f);
    }
}

/// bcftools' `drop_hdr_genotypes`: the header without samples or FORMAT definitions.
fn drop_genotypes(hdr: &Header) -> anyhow::Result<Header> {
    let (mut h, _) = hdr.subset(&[])?;
    h.remove(HeaderRecordType::Format, None);
    Ok(h)
}

/// C's `strtol(s + 1, &end, 10)` on the text after the CHROM column: the position and whether
/// any digits were read.
fn parse_pos(line: &[u8], tab: usize) -> Option<i64> {
    let start = tab + 1;
    if start > line.len() {
        return None;
    }
    let (v, end) = pasta::variant::number::strtol(line, start);
    (end != start).then_some(v)
}

impl Concat {
    /// Closes the output after concatenating (`result`), or abandons it after an error. After
    /// a failed read the output is closed, keeping the records before it as bcftools does, and
    /// the read error is returned unless closing fails.
    fn finish(
        &self,
        writer: Writer,
        result: anyhow::Result<()>,
        index_name: &str,
    ) -> anyhow::Result<ExitCode> {
        match result {
            Ok(()) => self.close(writer, index_name),
            Err(e) if e.is::<output::InputReadFailure>() => {
                self.close(writer, index_name)?;
                Err(e)
            }
            Err(e) => {
                writer.abandon();
                if e.is::<ExitOne>() {
                    return Ok(ExitCode::from(1));
                }
                Err(e)
            }
        }
    }

    /// Flushes and closes the output, saving its index.
    fn close(&self, writer: Writer, index_name: &str) -> anyhow::Result<ExitCode> {
        match writer.finish() {
            Ok(()) => Ok(ExitCode::SUCCESS),
            Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
            Err(WriteError::SaveIndex(_)) => bail!("Error: cannot write to index {index_name}"),
            Err(e @ WriteError::Deferred(_)) => Err(output::write_failure(e, || {
                self.last_write.get().failure(&self.output_name)
            })),
            Err(e) => Err(output::write_failure(e, || {
                anyhow::anyhow!("Error: close failed .. {}", self.output_name)
            })),
        }
    }

    /// Writes a record, translated to the output header.
    fn translate(&self, rec: &mut Record, hdr: &mut Header) -> anyhow::Result<()> {
        if self.sites_only {
            rec.subset(&[]);
        }
        if rec.translate(&self.out_hdr, hdr).is_err() {
            return Err(ExitOne.into());
        }
        Ok(())
    }

    /// Plain concatenation (bcftools' `concat` without `-a` or `-l`). With `--threads`,
    /// records are formatted or encoded on worker threads.
    fn concat_plain(&mut self, writer: &mut Writer) -> anyhow::Result<()> {
        let mut emitter = Emitter::new(&self.out_hdr, writer, self.threads, String::new())
            .with_write_error(LastWrite::Record.failure("").to_string());
        let fnames = std::mem::take(&mut self.fnames);
        let result = self.concat_files(&fnames, writer, &mut emitter);
        self.fnames = fnames;
        match result {
            Ok(()) => emitter.finish(&self.out_hdr, writer),
            // The records before a failed read are written, as without threads.
            Err(e) if e.is::<output::InputReadFailure>() => {
                emitter.finish(&self.out_hdr, writer)?;
                Err(e)
            }
            Err(e) => {
                emitter.abort(&self.out_hdr, writer);
                Err(e)
            }
        }
    }

    fn concat_files(
        &mut self,
        fnames: &[String],
        writer: &mut Writer,
        emitter: &mut Emitter,
    ) -> anyhow::Result<()> {
        let mut order = BlockOrder {
            prev_chr_id: -1,
            prev_pos: -1,
        };
        for f in fnames {
            let t0 = Instant::now();
            if self.verbose {
                eprint!("Concatenating {f}");
            }
            // Copied blocks and lines are written directly, after the records queued before
            // them.
            emitter.finish(&self.out_hdr, writer)?;
            if !self.copy_blocks(f, writer, &mut order)? {
                self.concat_file(f, writer, &mut order, emitter)?;
            }
            if self.verbose {
                eprintln!("\t{:.6} seconds", t0.elapsed().as_secs_f64());
            }
        }
        Ok(())
    }

    /// Concatenates one file, reading and writing its records or lines.
    fn concat_file(
        &mut self,
        f: &str,
        writer: &mut Writer,
        order: &mut BlockOrder,
        emitter: &mut Emitter,
    ) -> anyhow::Result<()> {
        let (mut reader, hdr) = open_input(f, self.threads, Reading::Again)?;
        let mut hdr = if self.sites_only {
            drop_genotypes(&hdr)?
        } else {
            hdr
        };
        let mut line = Vec::new();
        if reader.format() == InputFormat::Vcf && self.out_type.is_vcf() && !self.has_index {
            // VCF in and out: lines are copied as they are.
            while reader.read_line(&mut line).map_err(|_| bcf_read_error(f))? {
                if self.sites_only {
                    let mut ntab = 0;
                    if let Some(i) = line.iter().position(|&c| {
                        if c == b'\t' {
                            ntab += 1;
                        }
                        c == 0 || (c == b'\t' && ntab == 8)
                    }) {
                        line.truncate(i);
                    }
                }
                order.check_line(&self.out_hdr, &mut self.seen_seq, &line, f)?;
                line.push(b'\n');
                let written = LastWrite::Line(line.len());
                self.last_write.set(written);
                if let Err(e) = writer.write_raw(&line) {
                    return Err(output::write_failure(e, || written.failure("")));
                }
            }
            warn_if_no_eof(&reader, f);
            return Ok(());
        }
        let mut rec = emitter.spare().unwrap_or_default();
        loop {
            match reader.read_record(&mut hdr, &mut line, &mut rec) {
                Ok(true) => {}
                Ok(false) => {
                    warn_if_no_eof(&reader, f);
                    break;
                }
                // An unparsable VCF line ends the file quietly, as in bcftools (HTSlib reports
                // it); a failed read does not.
                Err(ReadError::Parse(_)) => break,
                Err(_) => return Err(bcf_read_error(f)),
            }
            self.translate(&mut rec, &mut hdr)?;
            order.check_record(&self.out_hdr, &mut self.seen_seq, &rec, f)?;
            self.last_write.set(LastWrite::Record);
            let next = emitter.spare().unwrap_or_default();
            if let Some(back) =
                emitter.write(&self.out_hdr, writer, std::mem::replace(&mut rec, next))?
            {
                rec = back;
            }
        }
        Ok(())
    }

    /// Concatenates one bgzipped file by copying its compressed blocks, recompressing only
    /// what reading changes. Returns `false` when the file needs [`Self::concat_file`]: other
    /// input or output types, a stream, `-G`, or an index to write.
    fn copy_blocks(
        &mut self,
        f: &str,
        writer: &mut Writer,
        order: &mut BlockOrder,
    ) -> anyhow::Result<bool> {
        if is_stream(f) || self.sites_only || self.has_index {
            return Ok(false);
        }
        let data = split_idx_name(f).0;
        let Ok(input) = RawInput::open(data) else {
            return Ok(false);
        };
        let is_vcf = match (input.format, input.compression, self.out_type) {
            (Format::Vcf, Compression::Bgzf, OutputType::VcfBgzf) => true,
            (Format::Bcf, Compression::Bgzf, OutputType::Bcf) => false,
            _ => return Ok(false),
        };
        let mut copier = FrameCopier::new(InflatingReader::new(input.source, data, self.threads));
        if !is_vcf {
            let Some((hdr, start)) = copier.read_bcf_header() else {
                bail!("\nFailed to parse header: {f}");
            };
            let mut thdr = hdr.clone();
            let (out_hdr, seen) = (&self.out_hdr, &mut self.seen_seq);
            self.last_write.set(LastWrite::Record);
            let write_err = || LastWrite::Record.failure("");
            copier.copy_bcf_records(start, &hdr, out_hdr, writer, &write_err, |rec| {
                if rec.translate(out_hdr, &mut thdr).is_err() {
                    return Err(ExitOne.into());
                }
                order.check_record(out_hdr, seen, rec, f)
            })?;
            if copier.read_failed() {
                return Err(bcf_read_error(f));
            }
            return Ok(true);
        }
        // The header: lines up to the first that starts with one '#'.
        let mut p = 0u64;
        while let Some(end) = copier.line_end(p) {
            let line = copier.slice(p, end);
            p = end;
            if line.len() > 1 && line[0] == b'#' && line[1] != b'#' {
                break;
            }
        }
        copier.skip_to(p);
        let mut line = Vec::new();
        // The blocks copied as they are count as a line of 0 bytes.
        self.last_write.set(LastWrite::Line(0));
        let write_err = |n: usize| LastWrite::Line(n).failure("");
        while let Some(end) = copier.line_end(p) {
            line.clear();
            line.extend_from_slice(copier.slice(p, end));
            let had_newline = line.last() == Some(&b'\n');
            if !had_newline && copier.read_failed() {
                // A line cut off by a truncated or corrupt block is not read (`hts_getline`
                // fails).
                break;
            }
            if had_newline {
                line.pop();
            }
            let had_cr = line.last() == Some(&b'\r');
            if had_cr {
                line.pop();
            }
            order.check_line(&self.out_hdr, &mut self.seen_seq, &line, f)?;
            if had_cr || !had_newline {
                line.push(b'\n');
                copier
                    .replace(writer, p, end, &line)
                    .map_err(|e| output::write_failure(e, || write_err(line.len())))?;
            }
            p = end;
            copier
                .advance(writer, p)
                .map_err(|e| output::write_failure(e, || write_err(0)))?;
        }
        copier
            .emit_original(writer, p)
            .map_err(|e| output::write_failure(e, || write_err(0)))?;
        if copier.read_failed() {
            return Err(bcf_read_error(f));
        }
        Ok(true)
    }

    /// Sets up the synced reader for `-a` (bcftools' `init_data`).
    fn init_overlaps(
        &mut self,
        files: &mut SyncedReaders,
        regions: Option<(String, ListSource)>,
        overlap: Overlap,
        remove_dups: Option<&str>,
    ) -> anyhow::Result<()> {
        if let Some((spec, source)) = regions {
            let Ok(r) = Regions::new(spec.as_bytes(), source) else {
                bail!("Failed to read the regions: {spec}");
            };
            files.set_regions_overlap(overlap);
            files.set_regions(r);
        }
        if let Some(d) = remove_dups {
            files.set_pairing(match d {
                "snps" => Pairing::SNPS,
                "indels" => Pairing::INDELS,
                "both" => Pairing::BOTH,
                "any" | "all" => Pairing::ANY,
                "none" | "exact" => Pairing::empty(),
                _ => bail!("The -D string \"{d}\" not recognised."),
            });
        }
        for f in &self.fnames {
            if let Err(e) = files.add_reader_again(f) {
                bail!("Failed to open {f}: {e}");
            }
        }
        Ok(())
    }

    /// `-a`: records of all files in position order (bcftools' `concat` with the synced
    /// reader).
    fn concat_overlaps(
        &mut self,
        files: &mut SyncedReaders,
        writer: &mut Writer,
        duplicates: Duplicates,
    ) -> anyhow::Result<()> {
        while next_line(files)? {
            for i in 0..files.reader_count() {
                let Some((line, hdr)) = files.line_with_header_mut(i) else {
                    continue;
                };
                self.translate(line, hdr)?;
                let written = LastWrite::Named("concat");
                self.last_write.set(written);
                if let Err(e) = writer.write_record(&self.out_hdr, line) {
                    return Err(output::write_failure(e, || {
                        written.failure(&self.output_name)
                    }));
                }
                if duplicates == Duplicates::Remove {
                    break;
                }
            }
        }
        if let Some(e) = files.take_pending_error() {
            bail!("Error: {e}");
        }
        Ok(())
    }
}

/// Moves the synced reader to the next position (`bcf_sr_next_line`): `false` at the end.
fn next_line(files: &mut SyncedReaders) -> anyhow::Result<bool> {
    match files.next_line() {
        Ok(n) => Ok(n > 0),
        Err(e) => bail!("Error: {e}"),
    }
}

/// The first reader with a record at the current position (bcftools' `_get_active_index`).
fn active_index(files: &SyncedReaders) -> Option<usize> {
    (0..files.reader_count()).find(|&i| files.has_line(i))
}

/// The checks of plain concatenation that each chromosome is one block in order (bcftools'
/// `prev_chr_id` and `prev_pos`, kept across files).
struct BlockOrder {
    prev_chr_id: i64,
    prev_pos: i64,
}

impl BlockOrder {
    fn check(
        &mut self,
        chr_id: i64,
        pos: i64,
        seen_seq: &mut [bool],
        chr: &dyn Fn() -> String,
    ) -> anyhow::Result<()> {
        if self.prev_chr_id != chr_id {
            self.prev_pos = -1;
            if seen_seq[chr_id as usize] {
                bail!(
                    "\nThe chromosome block {} is not contiguous, consider running with -a.",
                    chr()
                );
            }
        }
        if self.prev_pos > pos {
            bail!(
                "\nThe chromosome block {} is not sorted, consider running with -a.",
                chr()
            );
        }
        seen_seq[chr_id as usize] = true;
        self.prev_chr_id = chr_id;
        self.prev_pos = pos;
        Ok(())
    }

    /// Checks a VCF line (without its newline) copied as it is.
    fn check_line(
        &mut self,
        out_hdr: &Header,
        seen_seq: &mut [bool],
        line: &[u8],
        f: &str,
    ) -> anyhow::Result<()> {
        let tab = line
            .iter()
            .position(|&c| c == b'\t' || c == 0)
            .unwrap_or(line.len());
        let chr = &line[..tab];
        let name = || String::from_utf8_lossy(chr).into_owned();
        let Some(chr_id) = out_hdr.contig_id(chr) else {
            bail!(
                "\nThe sequence \"{}\" not defined in the header: {f}\n(Quick workaround: index the file.)",
                name()
            );
        };
        let Some(pos) = parse_pos(line, tab) else {
            bail!("Could not parse line: {}", String::from_utf8_lossy(line));
        };
        self.check(i64::from(chr_id), pos - 1, seen_seq, &name)
    }

    /// Checks a record translated to the output header.
    fn check_record(
        &mut self,
        out_hdr: &Header,
        seen_seq: &mut [bool],
        rec: &Record,
        f: &str,
    ) -> anyhow::Result<()> {
        let name = || String::from_utf8_lossy(rec.chrom_safe(out_hdr)).into_owned();
        let rid = i64::from(rec.rid());
        if usize::try_from(rid)
            .ok()
            .filter(|&r| r < seen_seq.len())
            .is_none()
        {
            bail!(
                "\nThe sequence \"{}\" not defined in the header: {f}",
                name()
            );
        }
        self.check(rid, rec.pos(), seen_seq, &name)
    }
}

#[cfg(test)]
mod tests {
    use pasta::variant::synced::SyncedError;

    use super::*;
    use crate::raw_input::{make_fifo, within_seconds};

    #[test]
    fn a_fifo_read_again_does_not_wait_for_a_writer() {
        let (dir, fifo) = make_fifo("concat-again");
        let name = fifo.to_string_lossy().into_owned();
        let again = {
            let name = name.clone();
            within_seconds("opening a FIFO again", move || {
                open_input(&name, 0, Reading::Again).map(|_| ())
            })
        };
        let synced = {
            let name = name.clone();
            within_seconds("adding a FIFO again", move || {
                SyncedReaders::new().add_reader_again(&name)
            })
        };
        let _ = std::fs::remove_dir_all(&dir);
        let msg = again.map_err(|e| e.to_string()).err();
        assert_eq!(msg, Some(format!("\nFailed to parse header: {name}")));
        assert!(matches!(synced, Err(SyncedError::FileType)), "{synced:?}");
    }

    #[test]
    fn min_pq_is_read_as_strtol_base_10() {
        assert_eq!(strtol_whole(""), Some(0));
        assert_eq!(strtol_whole(" +40"), Some(40));
        assert_eq!(strtol_whole("-5"), Some(-5));
        assert_eq!(strtol_whole("5x"), None);
        assert_eq!(strtol_whole(" "), None);
        assert_eq!(strtol_whole("4.5"), None);
    }

    #[test]
    fn verbose_is_read_as_strtol_base_0() {
        assert_eq!(strtol_base0_whole(""), Some(0));
        assert_eq!(strtol_base0_whole("0x1"), Some(1));
        assert_eq!(strtol_base0_whole("01"), Some(1));
        assert_eq!(strtol_base0_whole(" -0"), Some(0));
        assert_eq!(strtol_base0_whole("010"), Some(8));
        assert_eq!(strtol_base0_whole("08"), None);
        assert_eq!(strtol_base0_whole("0x"), None);
        assert_eq!(strtol_base0_whole("1 "), None);
        assert_eq!(strtol_base0_whole("99999999999999999999"), Some(i64::MAX));
    }
}
