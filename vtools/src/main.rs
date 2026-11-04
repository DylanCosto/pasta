//! `vtools`: a bcftools-compatible command-line toolkit built on pasta.

// std's print macros panic when a write fails with anything but EBADF (standard output or
// error on /dev/full, say), where bcftools' printf and fprintf carry on. These replace them in
// every module declared below; unit tests keep std's, which the test harness captures.
#[cfg(not(test))]
macro_rules! print {
    ($($arg:tt)*) => {
        $crate::print_to(&mut ::std::io::stdout().lock(), ::std::format_args!($($arg)*), "")
    };
}
#[cfg(not(test))]
macro_rules! println {
    () => {
        print!("\n")
    };
    ($($arg:tt)*) => {
        $crate::print_to(&mut ::std::io::stdout().lock(), ::std::format_args!($($arg)*), "\n")
    };
}
#[cfg(not(test))]
macro_rules! eprint {
    ($($arg:tt)*) => {
        $crate::print_to(&mut ::std::io::stderr().lock(), ::std::format_args!($($arg)*), "")
    };
}
#[cfg(not(test))]
macro_rules! eprintln {
    () => {
        eprint!("\n")
    };
    ($($arg:tt)*) => {
        $crate::print_to(&mut ::std::io::stderr().lock(), ::std::format_args!($($arg)*), "\n")
    };
}

mod closed_stdout;
mod commands;
mod convert;
mod faidx;
mod filter;
mod gff;
mod limits;
mod options;
mod output;
mod pipeline;
mod plugins;
mod raw_input;
mod shortopt;
mod sigpipe;
mod version;
mod workers;
mod zeroed;

use std::{
    ffi::OsString,
    io::{self, Write},
    process::ExitCode,
};

use clap::{
    ArgMatches, Command as ClapCommand, CommandFactory, FromArgMatches, Parser, Subcommand,
    error::ErrorKind,
};
use pasta::variant::{
    io::{reader::take_missing_eof, writer::dropped_unfinished},
    utils::CalcAcError,
};

use crate::shortopt::{OnError, Parsed, Style};

#[derive(Debug, Parser)]
#[command(name = "vtools", about, max_term_width = 80)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Convert, subset, and filter VCF/BCF files.
    View(Box<commands::view::Args>),
    /// Print the header and first records of a VCF/BCF file.
    Head(commands::head::Args),
    /// Index bgzipped VCF and BCF files, or report index statistics.
    Index(commands::index::Args),
    /// Apply fixed-threshold filters.
    Filter(Box<commands::filter::Args>),
    /// Create intersections, unions and complements of VCF files.
    Isec(Box<commands::isec::Args>),
    /// Print fields of VCF/BCF records in a user-defined format.
    Query(Box<commands::query::Args>),
    /// Convert VCF/BCF to other formats and back.
    Convert(Box<commands::convert::Args>),
    /// Concatenate or combine VCF/BCF files.
    Concat(Box<commands::concat::Args>),
    /// Create consensus sequence by applying VCF variants to a reference fasta file.
    Consensus(Box<commands::consensus::Args>),
    /// Modify the header, change sample names.
    Reheader(commands::reheader::Args),
    /// Sort VCF/BCF files.
    Sort(commands::sort::Args),
    /// Index or query bgzipped tab-delimited files (bcftools' tabix subcommand).
    Tabix(commands::tabix::Args),
    /// Left-align and normalize indels, split or join multiallelic sites.
    Norm(Box<commands::norm::Args>),
    /// Annotate and edit VCF/BCF files.
    Annotate(Box<commands::annotate::Args>),
    /// Produce VCF/BCF stats for plot-vcfstats.
    Stats(Box<commands::stats::Args>),
    /// SNP/indel variant calling from mpileup output.
    Call(Box<commands::call::Args>),
    /// Genotype likelihoods from alignments (SAM/BAM/CRAM).
    Mpileup(Box<commands::mpileup::Args>),
    /// Check sample identity by genotype discordance.
    Gtcheck(Box<commands::gtcheck::Args>),
    /// Merge VCF/BCF files from non-overlapping sample sets.
    Merge(Box<commands::merge::Args>),
    /// Identify runs of autozygosity (HMM).
    Roh(Box<commands::roh::Args>),
    /// Haplotype-aware consequence caller.
    Csq(Box<commands::csq::Args>),
    /// Run a built-in plugin (also `vtools +NAME`): counts, fill-tags, setGT, split-vep.
    Plugin(plugins::Args),
}

/// Writes `args` and then `end` to `out` for the print macros, in one write (standard error is
/// unbuffered), ignoring errors.
#[cfg(not(test))]
fn print_to(out: &mut impl io::Write, args: std::fmt::Arguments<'_>, end: &str) {
    let mut text = args.to_string();
    text.push_str(end);
    let _ = out.write_all(text.as_bytes());
}

fn main() -> ExitCode {
    sigpipe::restore();
    limits::fit_arenas();
    let args = match top_level(std::env::args_os().collect()) {
        TopLevel::Run(args) => args,
        TopLevel::Plugin(args) => return exit_status(plugins::run(&args)),
        TopLevel::Exit(code) => return code,
    };
    let name = args[1].to_string_lossy().into_owned();
    let matches = match Cli::command().try_get_matches_from(&args) {
        Ok(m) => m,
        Err(e) => return clap_error(&e, &name),
    };
    let cli = match Cli::from_arg_matches(&matches) {
        Ok(cli) => cli,
        Err(e) => return clap_error(&e, &name),
    };

    let sub = matches
        .subcommand()
        .map(|(_, sub)| sub)
        .expect("clap requires a subcommand");
    exit_status(run(cli.command, sub))
}

/// Runs a subcommand with its arguments and their matches.
fn run(command: Command, sub: &ArgMatches) -> anyhow::Result<ExitCode> {
    match command {
        Command::View(args) => commands::view::run(*args, sub),
        Command::Head(args) => commands::head::run(args, sub),
        Command::Index(args) => commands::index::run(args, sub),
        Command::Filter(args) => commands::filter::run(*args, sub),
        Command::Isec(args) => commands::isec::run(*args, sub),
        Command::Query(args) => commands::query::run(*args, sub),
        Command::Convert(args) => commands::convert::run(*args, sub),
        Command::Concat(args) => commands::concat::run(*args, sub),
        Command::Consensus(args) => commands::consensus::run(*args, sub),
        Command::Reheader(args) => commands::reheader::run(args, sub),
        Command::Sort(args) => commands::sort::run(args),
        Command::Tabix(args) => commands::tabix::run(args, sub),
        Command::Norm(args) => commands::norm::run(*args, sub),
        Command::Annotate(args) => commands::annotate::run(*args, sub),
        Command::Stats(args) => commands::stats::run(*args, sub),
        Command::Call(args) => commands::call::run(*args, sub),
        Command::Mpileup(args) => commands::mpileup::run(*args, sub),
        Command::Gtcheck(args) => commands::gtcheck::run(*args, sub),
        Command::Merge(args) => commands::merge::run(*args, sub),
        Command::Roh(args) => commands::roh::run(*args, sub),
        Command::Csq(args) => commands::csq::run(*args, sub),
        // `top_level` runs the plugin interface itself, with the arguments as given.
        Command::Plugin(_) => unreachable!("top_level runs `plugin` without clap"),
    }
}

/// The exit status of a subcommand's result, after printing its error. A command that
/// succeeded fails after all when it read a BGZF input to its end without finding the EOF
/// marker block (HTSlib only warns; README Damaged input).
fn exit_status(result: anyhow::Result<ExitCode>) -> ExitCode {
    // Only `finish` completes an output; a writer dropped without it is abandoned, which a
    // command that did not fail never does.
    debug_assert!(
        result.is_err() || dropped_unfinished() == 0,
        "an output was dropped without being finished"
    );
    match result {
        Ok(code) if code == ExitCode::SUCCESS => match take_missing_eof() {
            Some(name) => {
                let name = if name == "-" { "standard input" } else { &name };
                // Nothing more can be reported when standard error is closed.
                let _ = writeln!(
                    io::stderr(),
                    "Error: no BGZF EOF marker; the input may be truncated: {name}"
                );
                ExitCode::from(255)
            }
            None => code,
        },
        Ok(code) => code,
        // A closed downstream pipe (`| head`) when SIGPIPE was blocked: die of SIGPIPE, as
        // bcftools does when SIGPIPE is not blocked (see `sigpipe`).
        Err(e) if is_broken_pipe(&e) => sigpipe::die(),
        // HTSlib's `bcf_calc_ac` logs bad AN/AC counts and exits with status 1 (reached
        // through the AC, AF, MAC, and MAF of filter expressions).
        Err(e) if e.chain().any(|cause| cause.is::<CalcAcError>()) => ExitCode::from(1),
        Err(e) => {
            let msg = format!("{e:#}");
            if msg.ends_with('\n') {
                eprint!("{msg}");
            } else {
                eprintln!("{msg}");
            }
            // bcftools' `error()` exits with -1.
            ExitCode::from(255)
        }
    }
}

fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|cause| cause.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::BrokenPipe)
}

/// What the program's own arguments (before a subcommand's) ask for.
enum TopLevel {
    /// Run the subcommand: `vtools NAME` and the subcommand's arguments as clap reads them.
    Run(Vec<OsString>),
    /// Run the plugin interface with the arguments after `vtools plugin`, as given (clap would
    /// drop a leading `--`, which ends bcftools' options before any plugin name).
    Plugin(Vec<OsString>),
    /// Exit now with this status, after printing what was asked for.
    Exit(ExitCode),
}

/// The commands in `vtools`' usage, by group (bcftools' `cmds`, without the dropped commands
/// and without `tabix`, which bcftools does not list either).
const GROUPS: &[(&str, &[(&str, &str)])] = &[
    ("Indexing", &[("index", "index VCF/BCF files")]),
    (
        "VCF/BCF manipulation",
        &[
            ("annotate", "annotate and edit VCF/BCF files"),
            (
                "concat",
                "concatenate VCF/BCF files from the same set of samples",
            ),
            (
                "convert",
                "convert VCF/BCF files to different formats and back",
            ),
            ("head", "view VCF/BCF file headers"),
            ("isec", "intersections of VCF/BCF files"),
            (
                "merge",
                "merge VCF/BCF files files from non-overlapping sample sets",
            ),
            ("norm", "left-align and normalize indels"),
            ("plugin", "built-in plugins"),
            ("query", "transform VCF/BCF into user-defined formats"),
            ("reheader", "modify VCF/BCF header, change sample names"),
            ("sort", "sort VCF/BCF file"),
            (
                "view",
                "VCF/BCF conversion, view, subset and filter VCF/BCF files",
            ),
        ],
    ),
    (
        "VCF/BCF analysis",
        &[
            ("call", "SNP/indel calling"),
            (
                "consensus",
                "create consensus sequence by applying VCF variants",
            ),
            ("csq", "call variation consequences"),
            ("filter", "filter VCF/BCF files using fixed thresholds"),
            (
                "gtcheck",
                "check sample concordance, detect sample swaps and contamination",
            ),
            ("mpileup", "multi-way pileup producing genotype likelihoods"),
            ("roh", "identify runs of autozygosity (HMM)"),
            ("stats", "produce VCF/BCF stats"),
        ],
    ),
];

/// `vtools`' usage (bcftools' `usage`).
fn usage() -> String {
    use std::fmt::Write;

    let version = env!("CARGO_PKG_VERSION");
    let mut s = format!(
        "\nProgram: vtools (Tools for variant calling and manipulating VCFs and BCFs)\n\
         Version: {version} (using pasta {version}; compatible with bcftools 1.24)\n\n\
         Usage:   vtools [--version|--version-only] [--help] <command> <argument>\n\n\
         Commands:\n"
    );
    // Writing to a String cannot fail.
    for (group, commands) in GROUPS {
        let _ = write!(s, "\n -- {group}\n");
        for (name, help) in *commands {
            let _ = writeln!(s, "    {name:<12} {help}");
        }
    }
    let _ = write!(
        s,
        "\n -- Plugins (collection of programs for calling, file manipulation & analysis)\n    \
         {} plugins available, run \"vtools plugin -lv\" to see a complete list\n\n \
         Most commands accept VCF, bgzipped VCF, and BCF with the file type detected\n \
         automatically even when streaming from a pipe. Indexed VCF and BCF will work\n \
         in all situations. Un-indexed VCF and BCF and streams will work in most but\n \
         not all situations.\n\n",
        plugins::count()
    );
    s
}

/// Handles the program's own arguments as bcftools' `main` does: no arguments, `--version`,
/// `--version-only`, `help [COMMAND]`, `+PLUGIN`, and an unknown command. A subcommand's
/// arguments are then parsed as getopt does (see [`shortopt`]).
fn top_level(mut argv: Vec<OsString>) -> TopLevel {
    let Some(first) = argv.get(1) else {
        eprint!("{}", usage());
        return TopLevel::Exit(ExitCode::from(1));
    };
    let version = env!("CARGO_PKG_VERSION");
    let mut plus_name = true;
    match first.to_str() {
        Some("version" | "--version" | "-v") => {
            print!(
                "vtools {version}\nUsing pasta {version}\nCompatible with bcftools 1.24\n\
                 License Expat: The MIT/Expat license\n\
                 This is free software: you are free to change and redistribute it.\n\
                 There is NO WARRANTY, to the extent permitted by law.\n"
            );
            return TopLevel::Exit(ExitCode::SUCCESS);
        }
        Some("--version-only") => {
            println!("{version}+pasta-{version}");
            return TopLevel::Exit(ExitCode::SUCCESS);
        }
        Some("help" | "--help" | "-h") => {
            if argv.len() == 2 {
                print!("{}", usage());
                return TopLevel::Exit(ExitCode::SUCCESS);
            }
            // `vtools help COMMAND [...]` runs `vtools COMMAND`, which prints its usage
            // (`help +NAME` is an unknown command in bcftools).
            argv.truncate(3);
            argv.remove(1);
            plus_name = false;
        }
        _ => {}
    }
    let argv = if plus_name {
        plugins::expand_plus_name(argv)
    } else {
        argv
    };
    let cli = Cli::command();
    let sub = argv[1]
        .to_str()
        .filter(|name| *name != "help")
        .and_then(|name| cli.find_subcommand(name));
    let Some(sub) = sub else {
        let name = argv[1].to_string_lossy();
        eprintln!("[E::main] unrecognized command '{name}'");
        if let Some(similar) = similar_command(&cli, &name) {
            eprintln!("  tip: a similar command exists: '{similar}'");
        }
        return TopLevel::Exit(ExitCode::from(1));
    };
    let name = sub.get_name().to_owned();
    if name == "plugin" {
        // The plugin interface parses its own arguments.
        return TopLevel::Plugin(argv[2..].to_vec());
    }
    let mut sub = sub.clone().bin_name(format!("vtools {name}"));
    sub.build();
    // bcftools' tabix uses plain getopt and only warns about a bad option.
    let (style, on_error) = if name == "tabix" {
        (Style::Short, OnError::Ignore)
    } else {
        (Style::Long, OnError::Usage)
    };
    let run = |args: Vec<OsString>| {
        let mut all = Vec::with_capacity(args.len() + 2);
        all.extend(argv[..2].iter().cloned());
        all.extend(args);
        TopLevel::Run(all)
    };
    match shortopt::rewrite(&sub, &name, &argv[2..], style, on_error) {
        Parsed::Args { args, warnings } => {
            for w in warnings {
                eprint!("{w}");
            }
            run(args)
        }
        Parsed::Help => help(&sub, &name),
        Parsed::Error(msg) => {
            eprint!("{msg}");
            // bcftools then prints the command's usage, as its `-h`/`--help` does.
            let help_flag = [OsString::from("--help")];
            match shortopt::rewrite(&sub, &name, &help_flag, style, on_error) {
                Parsed::Args { args, .. } => run(args),
                Parsed::Help | Parsed::Error(_) => help(&sub, &name),
            }
        }
    }
}

/// The exit status of subcommand `name`'s usage errors: 255 where bcftools reports them with
/// `error()` or `exit(-1)` (call and csq), else 1.
fn usage_status(name: &str) -> u8 {
    match name {
        "call" | "csq" => 255,
        _ => 1,
    }
}

/// Prints the subcommand's options to standard error and exits as bcftools' `usage()` does.
fn help(sub: &ClapCommand, name: &str) -> TopLevel {
    eprint!("{}", sub.clone().render_long_help());
    TopLevel::Exit(ExitCode::from(usage_status(name)))
}

/// The visible command closest to `name`, if one is within a few edits of it (a third of its
/// length, at least 1), for a hint after an unrecognized command.
fn similar_command<'a>(cli: &'a ClapCommand, name: &str) -> Option<&'a str> {
    cli.get_subcommands()
        .filter(|sub| !sub.is_hide_set())
        .map(|sub| (edit_distance(name, sub.get_name()), sub.get_name()))
        .filter(|&(distance, candidate)| distance <= (candidate.len() / 3).max(1))
        .min_by_key(|&(distance, _)| distance)
        .map(|(_, candidate)| candidate)
}

/// The number of insertions, deletions, substitutions, and swaps of adjacent characters that
/// turn `a` into `b` (optimal string alignment distance).
fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let mut d = vec![vec![0; b.len() + 1]; a.len() + 1];

    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i;
    }
    for (j, cell) in d[0].iter_mut().enumerate() {
        *cell = j;
    }

    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            let mut best = (d[i - 1][j] + 1)
                .min(d[i][j - 1] + 1)
                .min(d[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                best = best.min(d[i - 2][j - 2] + 1);
            }
            d[i][j] = best;
        }
    }

    d[a.len()][b.len()]
}

/// Reports an error clap found in the rewritten arguments (a value it cannot parse, or a
/// missing or conflicting argument) with bcftools' exit statuses rather than clap's 2.
fn clap_error(e: &clap::Error, name: &str) -> ExitCode {
    // Nothing more can be reported when standard error is closed.
    let _ = e.print();
    match e.kind() {
        ErrorKind::DisplayHelp | ErrorKind::DisplayVersion => ExitCode::SUCCESS,
        ErrorKind::InvalidValue | ErrorKind::ValueValidation | ErrorKind::InvalidUtf8 => {
            ExitCode::from(255)
        }
        _ => ExitCode::from(usage_status(name)),
    }
}

#[cfg(test)]
mod tests {
    use pasta::variant::io::{OutputType, WriteOptions, Writer};

    use super::*;

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "an output was dropped without being finished")]
    fn a_command_that_succeeds_without_finishing_its_output_is_caught() {
        let options = WriteOptions {
            output_type: OutputType::VcfBgzf,
            level: None,
            threads: 2,
        };
        let mut writer = Writer::from_sink(Box::new(io::sink()), &options);
        writer.write_raw(b"a record\n").expect("write to a sink");
        // A command that forgot `finish` and then succeeded.
        drop(writer);
        exit_status(Ok(ExitCode::SUCCESS));
    }
}
