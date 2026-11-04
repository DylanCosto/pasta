//! `ntools`: a samtools-compatible command-line toolkit built on pasta.

mod alignment;
mod bam_blocks;
mod bam_raw;
mod bam_sam;
mod base_mods;
mod bgzf_stream;
mod commands;
mod cram_file;
mod cram_raw;
mod depth;
mod expr;
mod file_list;
mod filter;
mod flags;
mod getopt;
mod global_args;
mod hts_region;
mod khash;
mod output;
mod parallel_cram;
mod parallel_sam;
mod program;
mod reference;
mod regions;
mod runs;
mod sanitize;
mod sequence;
mod sigpipe;
mod sort_order;
mod source;
mod template_coordinate;

use std::{
    io::{self, IsTerminal, Write},
    process::ExitCode,
};

use clap::{CommandFactory, FromArgMatches, Parser, Subcommand};

#[derive(Debug, Parser)]
// A repeated option takes its last value, as with samtools' getopt, rather than being an error.
// Options whose repeats samtools combines (view -f/-F/-G) are lists.
#[command(
    name = "ntools",
    version,
    long_version = concat!(
        env!("CARGO_PKG_VERSION"),
        "\nUsing pasta ",
        env!("CARGO_PKG_VERSION"),
        "\nCompatible with samtools 1.24"
    ),
    about,
    args_override_self = true,
    max_term_width = 80
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Convert, filter, count, and query SAM/BAM/CRAM.
    View(Box<commands::view::Args>),
    /// Count records by flag category.
    Flagstat(commands::flagstat::Args),
    /// Report per-reference mapped and unmapped counts.
    Idxstats(commands::idxstats::Args),
    /// Build a BAI, CSI, or CRAI index.
    Index(commands::index::Args),
    /// Index a FASTA file or extract regions from it.
    Faidx(commands::faidx::Args),
    /// Index a FASTQ file or extract regions from it.
    Fqidx(commands::faidx::Args),
    /// Sort alignments by coordinate or read name.
    Sort(commands::sort::Args),
    /// Merge sorted SAM, BAM, or CRAM files.
    Merge(commands::merge::Args),
    /// Convert alignments to FASTQ.
    #[command(visible_alias = "bam2fq")]
    Fastq(commands::fastq::Args),
    /// Convert alignments to FASTA.
    Fasta(commands::fastq::Args),
    /// Check that files look intact.
    Quickcheck(commands::quickcheck::Args),
    /// Replace a BAM or CRAM header without recompressing records.
    Reheader(commands::reheader::Args),
    /// Concatenate BAM or CRAM files without recompressing records.
    Cat(commands::cat::Args),
    /// Report per-position read depth.
    Depth(commands::depth::Args),
    /// Summarize coverage per reference sequence.
    Coverage(commands::coverage::Args),
    /// Write one file per read group or tag value.
    Split(commands::split::Args),
    /// Group records by read name.
    Collate(commands::collate::Args),
    /// Fill in mate information from name-grouped alignments.
    Fixmate(commands::fixmate::Args),
    /// Mark duplicate alignments.
    Markdup(commands::markdup::Args),
    /// Recalculate MD and NM tags.
    Calmd(commands::calmd::Args),
    /// Text pileup of alignments.
    Mpileup(commands::mpileup::Args),
    /// Alignment statistics for plot-bamstats.
    Stats(commands::stats::Args),
    /// Write a SAM sequence dictionary for a FASTA file.
    Dict(commands::dict::Args),
    /// Sum read depth over BED regions.
    Bedcov(commands::bedcov::Args),
    /// Add or replace read group tags and header lines.
    Addreplacerg(commands::addreplacerg::Args),
    /// Clip primer regions from reads.
    Ampliconclip(commands::ampliconclip::Args),
    /// Print the header and first records as SAM.
    Head(commands::head::Args),
    /// Convert FASTQ or FASTA to unaligned SAM/BAM/CRAM.
    Import(commands::import::Args),
    /// Turn alignments back into unaligned reads.
    Reset(commands::reset::Args),
    /// Order-independent checksums of reads, qualities, and tags.
    Checksum(commands::checksum::Args),
    /// Convert between numeric and textual flags.
    Flags(commands::flags::Args),
    /// Report CRAM block sizes by content ID and data series.
    #[command(alias = "cram_size")]
    CramSize(commands::cram_size::Args),
    /// Rebuild a reference from MD tags or embedded CRAM references.
    Reference(commands::reference::Args),
    /// Remove PCR duplicates (superseded by markdup).
    Rmdup(commands::rmdup::Args),
    /// List the samples (read group SM values) of alignment files.
    Samples(commands::samples::Args),
    /// Generate a consensus sequence from aligned reads.
    #[command(disable_help_flag = true)]
    Consensus(commands::consensus::Args),
    /// Compress or decompress files in BGZF format (HTSlib's bgzip).
    #[command(disable_help_flag = true, disable_version_flag = true)]
    Bgzip(commands::bgzip::Args),
}

/// Prints a command-line error and exits as samtools does on a usage error, with status 1 (clap
/// uses 2). Help and version requests print to stdout and exit 0, as before.
fn exit_with_usage_error(e: clap::Error) -> ! {
    if !e.use_stderr() {
        e.exit();
    }

    let _ = e.print();
    std::process::exit(1)
}

/// When samtools prints a command's usage instead of running it.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum UsageWhen {
    /// Run without arguments.
    Bare,
    /// Run without arguments, with standard input a terminal.
    BareOnTerminal,
    /// Run without an input file, with standard input a terminal (and other options allowed).
    NoInputOnTerminal,
}

/// Where a usage goes, and the exit status after it.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum UsageExit {
    /// Standard output, exit 0.
    Success,
    /// Standard output, exit 1.
    StdoutFailure,
    /// Standard error, exit 1.
    Failure,
}

/// A command that prints its usage as samtools does, with the line samtools prints before it.
struct Usage {
    command: &'static str,
    when: UsageWhen,
    exit: UsageExit,
    message: Option<&'static str>,
}

impl Usage {
    const fn new(command: &'static str, when: UsageWhen, exit: UsageExit) -> Self {
        Self {
            command,
            when,
            exit,
            message: None,
        }
    }
}

/// The commands that print their usage when samtools does (rather than failing, or reading
/// standard input). The first entry that applies is used.
const USAGES: &[Usage] = &[
    Usage::new("coverage", UsageWhen::Bare, UsageExit::Success),
    Usage::new("depth", UsageWhen::Bare, UsageExit::Success),
    Usage::new("faidx", UsageWhen::Bare, UsageExit::Success),
    Usage::new("fixmate", UsageWhen::Bare, UsageExit::Success),
    Usage::new("flagstat", UsageWhen::Bare, UsageExit::Success),
    Usage::new("fqidx", UsageWhen::Bare, UsageExit::Success),
    Usage::new("idxstats", UsageWhen::Bare, UsageExit::Success),
    Usage::new("import", UsageWhen::Bare, UsageExit::Success),
    Usage::new("index", UsageWhen::Bare, UsageExit::Success),
    Usage::new("merge", UsageWhen::Bare, UsageExit::Success),
    Usage::new("mpileup", UsageWhen::Bare, UsageExit::Failure),
    Usage::new("reset", UsageWhen::BareOnTerminal, UsageExit::Success),
    Usage::new("checksum", UsageWhen::BareOnTerminal, UsageExit::Success),
    Usage::new("view", UsageWhen::BareOnTerminal, UsageExit::Success),
    Usage {
        message: Some("ntools view: No input provided or missing option argument.\n\n"),
        ..Usage::new("view", UsageWhen::NoInputOnTerminal, UsageExit::Failure)
    },
    Usage::new("sort", UsageWhen::NoInputOnTerminal, UsageExit::Success),
    Usage::new("fastq", UsageWhen::NoInputOnTerminal, UsageExit::Success),
    Usage::new("fasta", UsageWhen::NoInputOnTerminal, UsageExit::Success),
    Usage::new("head", UsageWhen::NoInputOnTerminal, UsageExit::Success),
    Usage::new(
        "cram-size",
        UsageWhen::NoInputOnTerminal,
        UsageExit::Success,
    ),
    Usage::new(
        "reference",
        UsageWhen::NoInputOnTerminal,
        UsageExit::Success,
    ),
    Usage::new(
        "stats",
        UsageWhen::NoInputOnTerminal,
        UsageExit::StdoutFailure,
    ),
    Usage::new("dict", UsageWhen::NoInputOnTerminal, UsageExit::Failure),
    Usage::new("samples", UsageWhen::NoInputOnTerminal, UsageExit::Failure),
];

/// Prints the usage of `ntools NAME` if samtools would print its usage for these arguments
/// (see [`USAGES`]), and returns the exit status. `matches` are the parsed arguments, or `None`
/// before parsing, when only the bare cases are checked.
fn print_usage(
    args: &[std::ffi::OsString],
    matches: Option<&clap::ArgMatches>,
) -> Option<ExitCode> {
    let name = args.get(1)?.to_str()?;
    let bare = args.len() == 2;

    let mut cli = Cli::command();
    cli.build();
    let command = cli.find_subcommand(name)?;

    let applies = |when| match (when, matches) {
        (UsageWhen::Bare, None) => bare,
        (UsageWhen::BareOnTerminal, None) => bare && io::stdin().is_terminal(),
        (UsageWhen::NoInputOnTerminal, Some(matches)) => {
            let Some((_, matches)) = matches.subcommand() else {
                return false;
            };

            io::stdin().is_terminal()
                && command.get_positionals().all(|arg| {
                    matches.value_source(arg.get_id().as_str())
                        != Some(clap::parser::ValueSource::CommandLine)
                })
        }
        _ => false,
    };

    let usage = USAGES
        .iter()
        .find(|usage| usage.command == command.get_name() && applies(usage.when))?;

    let mut text = usage.message.unwrap_or_default().to_owned();
    text.push_str(&command.clone().render_help().to_string());

    let (result, status) = match usage.exit {
        UsageExit::Success => (io::stdout().write_all(text.as_bytes()), ExitCode::SUCCESS),
        UsageExit::StdoutFailure => (io::stdout().write_all(text.as_bytes()), ExitCode::FAILURE),
        UsageExit::Failure => (io::stderr().write_all(text.as_bytes()), ExitCode::FAILURE),
    };

    if let Err(e) = result
        && e.kind() == io::ErrorKind::BrokenPipe
    {
        sigpipe::die();
    }

    Some(status)
}

/// Lets every option that requires a value take a negative number, as getopt does
/// (`view -q -5`), rather than reading it as an option.
fn allow_negative_values(cli: clap::Command) -> clap::Command {
    cli.mut_subcommands(|command| {
        command.mut_args(|arg| {
            let optional = arg.get_num_args().is_some_and(|n| n.min_values() == 0);

            if arg.is_positional() || optional || !arg.get_action().takes_values() {
                arg
            } else {
                arg.allow_negative_numbers(true)
            }
        })
    })
}

fn main() -> ExitCode {
    let args = expand_output_format(std::env::args_os());

    if let Some(status) = print_usage(&args, None) {
        return status;
    }

    let matches = allow_negative_values(Cli::command())
        .try_get_matches_from(&args)
        .unwrap_or_else(|e| exit_with_usage_error(e));

    if let Some(status) = print_usage(&args, Some(&matches)) {
        return status;
    }

    let mut cli = Cli::from_arg_matches(&matches).unwrap_or_else(|e| exit_with_usage_error(e));

    // Some options apply in command-line order, which the parsed arguments do not keep.
    if let (Command::View(args), Some(("view", view_matches))) =
        (&mut cli.command, matches.subcommand())
    {
        args.record_order(view_matches);
    }

    let result = match cli.command {
        Command::View(args) => commands::view::run(*args).map(|()| ExitCode::SUCCESS),
        Command::Flagstat(args) => commands::flagstat::run(args).map(|()| ExitCode::SUCCESS),
        Command::Idxstats(args) => commands::idxstats::run(args).map(|()| ExitCode::SUCCESS),
        Command::Index(args) => commands::index::run(args).map(|()| ExitCode::SUCCESS),
        Command::Faidx(args) => commands::faidx::run(args),
        Command::Fqidx(args) => commands::faidx::run_fqidx(args),
        Command::Sort(args) => commands::sort::run(args).map(|()| ExitCode::SUCCESS),
        Command::Merge(args) => commands::merge::run(args).map(|()| ExitCode::SUCCESS),
        Command::Fastq(args) => commands::fastq::run(args, commands::fastq::ReadFormat::Fastq)
            .map(|()| ExitCode::SUCCESS),
        Command::Fasta(args) => commands::fastq::run(args, commands::fastq::ReadFormat::Fasta)
            .map(|()| ExitCode::SUCCESS),
        Command::Quickcheck(args) => commands::quickcheck::run(args),
        Command::Reheader(args) => commands::reheader::run(args).map(|()| ExitCode::SUCCESS),
        Command::Cat(args) => commands::cat::run(args).map(|()| ExitCode::SUCCESS),
        Command::Depth(args) => commands::depth::run(args).map(|()| ExitCode::SUCCESS),
        Command::Coverage(args) => commands::coverage::run(args).map(|()| ExitCode::SUCCESS),
        Command::Split(args) => commands::split::run(args).map(|()| ExitCode::SUCCESS),
        Command::Collate(args) => commands::collate::run(args).map(|()| ExitCode::SUCCESS),
        Command::Fixmate(args) => commands::fixmate::run(args).map(|()| ExitCode::SUCCESS),
        Command::Markdup(args) => commands::markdup::run(args).map(|()| ExitCode::SUCCESS),
        Command::Calmd(args) => commands::calmd::run(args).map(|()| ExitCode::SUCCESS),
        Command::Mpileup(args) => commands::mpileup::run(args).map(|()| ExitCode::SUCCESS),
        Command::Stats(args) => commands::stats::run(args).map(|()| ExitCode::SUCCESS),
        Command::Dict(args) => commands::dict::run(args).map(|()| ExitCode::SUCCESS),
        Command::Bedcov(args) => commands::bedcov::run(args),
        Command::Addreplacerg(args) => {
            commands::addreplacerg::run(args).map(|()| ExitCode::SUCCESS)
        }
        Command::Ampliconclip(args) => {
            commands::ampliconclip::run(args).map(|()| ExitCode::SUCCESS)
        }
        Command::Reset(args) => commands::reset::run(args).map(|()| ExitCode::SUCCESS),
        Command::Checksum(args) => commands::checksum::run(args).map(|()| ExitCode::SUCCESS),
        Command::Import(args) => commands::import::run(args).map(|()| ExitCode::SUCCESS),
        Command::Head(args) => commands::head::run(args).map(|()| ExitCode::SUCCESS),
        Command::Flags(args) => commands::flags::run(args),
        Command::CramSize(args) => commands::cram_size::run(args).map(|()| ExitCode::SUCCESS),
        Command::Reference(args) => commands::reference::run(args).map(|()| ExitCode::SUCCESS),
        Command::Rmdup(args) => commands::rmdup::run(args).map(|()| ExitCode::SUCCESS),
        Command::Samples(args) => commands::samples::run(args),
        Command::Consensus(args) => commands::consensus::run(args),
        Command::Bgzip(args) => commands::bgzip::run(args),
    };

    match result {
        Ok(code) => code,
        // A closed downstream pipe (`| head`): die of SIGPIPE, with no message, as samtools
        // does. The command's temporary files are gone by now (see `sigpipe`).
        Err(e) if is_broken_pipe(&e) => sigpipe::die(),
        Err(e) => {
            // As samtools, die of SIGPIPE when standard error is a closed pipe.
            if let Err(e) = writeln!(io::stderr(), "ntools: {e:#}")
                && e.kind() == io::ErrorKind::BrokenPipe
            {
                sigpipe::die();
            }
            ExitCode::FAILURE
        }
    }
}

fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|cause| cause.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::BrokenPipe)
}

/// Splits HTSlib's `-O FORMAT,OPT,...` (or `--output-fmt`) into `-O FORMAT` and one
/// `--output-fmt-option OPT` per option, as HTSlib parses it. Only values that start with an
/// alignment format and a comma are split, so other `-O` options (such as `fastq -O`) are left
/// alone.
fn expand_output_format<I>(args: I) -> Vec<std::ffi::OsString>
where
    I: IntoIterator<Item = std::ffi::OsString>,
{
    use std::ffi::OsString;

    fn split(value: &str) -> Option<(&str, std::str::Split<'_, char>)> {
        let (format, options) = value.split_once(',')?;
        ["sam", "bam", "cram"]
            .iter()
            .any(|name| format.eq_ignore_ascii_case(name))
            .then(|| (format, options.split(',')))
    }

    let mut expanded = Vec::new();
    let mut args = args.into_iter();

    while let Some(arg) = args.next() {
        let Some(s) = arg.to_str() else {
            expanded.push(arg);
            continue;
        };

        if s == "--" {
            expanded.push(arg);
            expanded.extend(args.by_ref());
            break;
        }

        let (flag, value) = match s {
            "-O" | "--output-fmt" => match args.next() {
                Some(value) => (s.to_owned(), value),
                None => {
                    expanded.push(arg);
                    continue;
                }
            },
            _ if s.starts_with("--output-fmt=") => (
                "--output-fmt".to_owned(),
                OsString::from(&s["--output-fmt=".len()..]),
            ),
            _ if s.starts_with("-O") && s.len() > 2 => ("-O".to_owned(), OsString::from(&s[2..])),
            _ => {
                expanded.push(arg);
                continue;
            }
        };

        match value.to_str().and_then(split) {
            Some((format, options)) => {
                expanded.push(OsString::from(flag));
                expanded.push(OsString::from(format));

                for option in options.filter(|option| !option.is_empty()) {
                    expanded.push(OsString::from("--output-fmt-option"));
                    expanded.push(OsString::from(option));
                }
            }
            None if s.starts_with("--output-fmt=") || (s.starts_with("-O") && s.len() > 2) => {
                expanded.push(arg);
            }
            None => {
                expanded.push(OsString::from(flag));
                expanded.push(value);
            }
        }
    }

    expanded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cli() {
        allow_negative_values(Cli::command()).debug_assert();
    }

    #[test]
    fn test_expand_output_format() {
        let expand = |args: &[&str]| -> Vec<String> {
            expand_output_format(args.iter().map(std::ffi::OsString::from))
                .into_iter()
                .map(|arg| arg.into_string().unwrap())
                .collect()
        };

        assert_eq!(
            expand(&["ntools", "view", "-O", "cram,version=3.0,small", "x.bam"]),
            [
                "ntools",
                "view",
                "-O",
                "cram",
                "--output-fmt-option",
                "version=3.0",
                "--output-fmt-option",
                "small",
                "x.bam"
            ]
        );
        assert_eq!(
            expand(&["ntools", "sort", "--output-fmt=BAM,level=1"]),
            [
                "ntools",
                "sort",
                "--output-fmt",
                "BAM",
                "--output-fmt-option",
                "level=1"
            ]
        );
        assert_eq!(
            expand(&["ntools", "view", "-Ocram,level=9"]),
            [
                "ntools",
                "view",
                "-O",
                "cram",
                "--output-fmt-option",
                "level=9"
            ]
        );

        // Left alone: plain formats, other -O options, and anything after `--`.
        for args in [
            &["ntools", "view", "-O", "bam", "x.bam"][..],
            &["ntools", "fastq", "-O", "x,y.bam"],
            &["ntools", "flagstat", "-O", "json", "x.bam"],
            &["ntools", "view", "--", "-O", "cram,level=1"],
        ] {
            assert_eq!(expand(args), args);
        }
    }
}
