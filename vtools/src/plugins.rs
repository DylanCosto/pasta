//! bcftools' plugin interface (`vcfplugin.c`): `vtools +NAME [OPTIONS] FILE [-- PLUGIN_OPTIONS]`
//! and `vtools plugin NAME ...`, with the plugins built in rather than loaded from shared
//! libraries.
//!
//! A plugin either takes over the command line (bcftools' `run`, as `+split-vep` does) or is
//! driven by this module (`init`, `process`, `destroy`): the input is read with the shared
//! options (`-i/-e`, `-r/-R/-t/-T`), each record goes through the plugin, and the records it
//! keeps are written with `-o/-O/--threads/-W`.

mod counts;
mod fill_tags;
mod set_gt;
mod split_vep;

use std::{
    ffi::OsString,
    fmt,
    io::{self, IsTerminal},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::bail;
use clap::{ArgMatches, CommandFactory, FromArgMatches, Parser, builder::TypedValueParser as _};
use pasta::variant::{
    Header, Record, Unpack,
    io::{WriteError, Writer},
    log,
    readlist::ListSource,
    synced::{Overlap, Regions, SyncedError, SyncedOptions, SyncedReader},
};

use crate::options::later;
use crate::{
    commands::{
        head::apply_verbosity,
        norm::{emit::Emitter, input::ParallelVcf},
    },
    filter::{Filter, FilterLogic},
    output::{self, OutputSetting},
    raw_input::strerror,
    shortopt::{self, OnError, Parsed, Style},
    version,
};

/// The arguments of `vtools plugin`, for clap's list of commands only: `top_level` hands the
/// arguments to [`run`] as given (clap would drop a leading `--`), and each plugin parses its own.
#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// The plugin name, the shared options, the input file, and `--` with the plugin's options.
    #[arg(trailing_var_arg = true, allow_hyphen_values = true, num_args = 0..)]
    pub(crate) args: Vec<OsString>,
}

/// Rewrites `vtools +NAME ...` as `vtools plugin NAME ...` (bcftools' `main`).
pub fn expand_plus_name(mut args: Vec<OsString>) -> Vec<OsString> {
    let name = args
        .get(1)
        .and_then(|a| a.to_str())
        .and_then(|a| a.strip_prefix('+'))
        .map(str::to_owned);
    if let Some(name) = name {
        args[1] = OsString::from(name);
        args.insert(1, OsString::from("plugin"));
    }
    args
}

/// An error that ends the program with a given exit status after printing its message, as
/// bcftools' `usage()` functions (status 1) and plugin usage errors do.
#[derive(Debug)]
pub(crate) struct Exit {
    /// The exit status.
    pub code: u8,
    /// The text printed to standard error, as is.
    pub message: String,
}

impl Exit {
    /// Exits with status 1 after printing `message`.
    pub fn usage(message: impl Into<String>) -> anyhow::Error {
        anyhow::Error::new(Self {
            code: 1,
            message: message.into(),
        })
    }

    /// Exits with `code` without printing anything, after a message HTSlib already logged
    /// (as its `exit(1)` calls do).
    pub fn status(code: u8) -> anyhow::Error {
        anyhow::Error::new(Self {
            code,
            message: String::new(),
        })
    }

    /// Exits with status 255 (bcftools' `error()`) after printing `message` as is.
    pub fn error(message: impl Into<String>) -> anyhow::Error {
        anyhow::Error::new(Self {
            code: 255,
            message: message.into(),
        })
    }
}

impl fmt::Display for Exit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Exit {}

/// A plugin driven record by record (bcftools' `init`, `process`, and `destroy`).
pub(crate) trait Plugin {
    /// Processes a record read with `hdr`; values are set with the output header `out`.
    /// Returns whether the record is written (bcftools' `process` returning the record).
    ///
    /// # Errors
    ///
    /// Returns the plugin's fatal errors.
    fn process(&mut self, hdr: &Header, out: &Header, rec: &mut Record) -> anyhow::Result<bool>;

    /// Finishes after the last record (bcftools' `destroy`), printing any summary.
    ///
    /// # Errors
    ///
    /// Returns errors writing the summary.
    fn destroy(&mut self) -> anyhow::Result<()> {
        Ok(())
    }
}

/// Whether the plugin's records are written (bcftools' `init` returning 0) or not (1).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Output {
    /// The header and the records the plugin keeps are written.
    Records,
    /// Nothing is written; the plugin prints its own results.
    None,
}

/// An initialized plugin.
pub(crate) struct Loaded {
    pub plugin: Box<dyn Plugin>,
    pub output: Output,
    /// How much of each record the plugin reads; empty for everything. A limit lets the VCF
    /// parser skip the rest of the line.
    pub max_unpack: Unpack,
}

/// The input and output headers given to a plugin's `init`, or none when bcftools calls it
/// without an input (from a terminal), where plugins only parse their options.
pub(crate) type InitHeaders<'a> = Option<(&'a Header, &'a mut Header)>;

type InitFn = fn(&[String], InitHeaders<'_>) -> anyhow::Result<Loaded>;
type RunFn = fn(&[String]) -> anyhow::Result<ExitCode>;

/// How a plugin is started.
#[derive(Clone, Copy)]
enum Entry {
    /// Driven record by record.
    Init(InitFn),
    /// Takes over the command line (bcftools' `run`).
    Run(RunFn),
}

/// A built-in plugin.
struct Info {
    name: &'static str,
    about: &'static str,
    usage: fn() -> String,
    entry: Entry,
}

/// The built-in plugins, sorted by name as `plugin -l` lists them.
const PLUGINS: [Info; 4] = [
    Info {
        name: "counts",
        about: counts::ABOUT,
        usage: counts::usage,
        entry: Entry::Init(counts::init),
    },
    Info {
        name: "fill-tags",
        about: fill_tags::ABOUT,
        usage: fill_tags::usage,
        entry: Entry::Init(fill_tags::init),
    },
    Info {
        name: "setGT",
        about: set_gt::ABOUT,
        usage: set_gt::usage,
        entry: Entry::Init(set_gt::init),
    },
    Info {
        name: "split-vep",
        about: split_vep::ABOUT,
        usage: split_vep::usage,
        entry: Entry::Run(split_vep::run),
    },
];

/// Reads the shared options as bcftools' getopt loop does, printing getopt's message for each
/// bad option: the arguments as clap is to read them, and whether there was a bad option.
fn general_args(prog: &str, args: &[String]) -> (Vec<OsString>, bool) {
    let mut cmd = General::command();
    cmd.build();
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    match shortopt::rewrite(&cmd, prog, &args, Style::Long, OnError::Ignore) {
        Parsed::Args { args, warnings } => {
            for w in &warnings {
                eprint!("{w}");
            }
            (args, !warnings.is_empty())
        }
        Parsed::Help | Parsed::Error(_) => {
            unreachable!("-h is a plain flag here and getopt errors are warnings")
        }
    }
}

/// The shared options (bcftools' `main_plugin` getopt string `h?o:O:r:R:t:T:li:e:v::VW::`).
#[derive(Debug, Parser)]
#[command(
    no_binary_name = true,
    disable_help_flag = true,
    disable_version_flag = true,
    args_override_self = true
)]
struct General {
    /// The input file, then arguments passed on to the plugin.
    inputs: Vec<String>,

    /// The plugin's options, after `--`.
    #[arg(last = true)]
    plugin_args: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(short = 'e', long = "exclude", allow_hyphen_values = true)]
    exclude: Vec<String>,

    /// Select sites for which the expression is true.
    #[arg(short = 'i', long = "include", allow_hyphen_values = true)]
    include: Vec<String>,

    /// List the plugin's options.
    #[arg(short = 'h', long = "help", short_alias = '?')]
    help: bool,

    /// List the available plugins.
    #[arg(short = 'l', long = "list-plugins")]
    list: bool,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_parser = clap::builder::OsStringValueParser::new().map(PathBuf::from))]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", action = clap::ArgAction::Append)]
    output_type: Vec<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(short = 'r', long = "regions")]
    regions: Option<String>,

    /// Restrict to regions listed in a file.
    #[arg(short = 'R', long = "regions-file")]
    regions_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(long = "regions-overlap", action = clap::ArgAction::Append)]
    regions_overlap: Vec<String>,

    /// Similar to -r but streams rather than index-jumps.
    #[arg(short = 't', long = "targets", allow_hyphen_values = true)]
    targets: Option<String>,

    /// Similar to -R but streams rather than index-jumps.
    #[arg(short = 'T', long = "targets-file")]
    targets_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0].
    #[arg(long = "targets-overlap", action = clap::ArgAction::Append)]
    targets_overlap: Vec<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Print version string and exit.
    #[arg(short = 'V', long = "version")]
    version: bool,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        visible_alias = "verbose",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "",
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,

    /// Automatically index the output files [off].
    #[arg(
        short = 'W',
        long = "write-index",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "",
        action = clap::ArgAction::Append
    )]
    write_index: Vec<String>,
}

/// bcftools' `usage()` of `main_plugin`.
fn general_usage() -> anyhow::Error {
    Exit::usage(
        "\n\
About:   Run user defined plugin\n\
Usage:   vtools plugin <name> [OPTIONS] <file> [-- PLUGIN_OPTIONS]\n\
\x20        vtools +name [OPTIONS] <file>  [-- PLUGIN_OPTIONS]\n\
\n\
VCF input options:\n\
\x20  -e, --exclude EXPR             Exclude sites for which the expression is true\n\
\x20  -i, --include EXPR             Select sites for which the expression is true\n\
\x20  -r, --regions REGION           Restrict to comma-separated list of regions\n\
\x20  -R, --regions-file FILE        Restrict to regions listed in a file\n\
\x20      --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n\
\x20  -t, --targets REGION           Similar to -r but streams rather than index-jumps\n\
\x20  -T, --targets-file FILE        Similar to -R but streams rather than index-jumps\n\
\x20      --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n\
VCF output options:\n\
\x20      --no-version               Do not append version and command line to the header\n\
\x20  -o, --output FILE              Write output to a file [standard output]\n\
\x20  -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n\
\x20      --threads INT              Use multithreading with <int> worker threads [0]\n\
Plugin options:\n\
\x20  -h, --help                     List plugin's options\n\
\x20  -l, --list-plugins             List available plugins (counts, fill-tags, setGT, split-vep)\n\
\x20  -v, --verbosity INT            Verbosity level\n\
\x20  -V, --version                  Print version string and exit\n\
\x20  -W, --write-index[=FMT]        Automatically index the output files [off]\n\
\n",
    )
}

/// Runs `vtools plugin ...` (bcftools' `main_plugin`).
///
/// # Errors
///
/// Returns bcftools' fatal errors; usage errors become their exit status.
pub fn run(args: &[OsString]) -> anyhow::Result<ExitCode> {
    let args: Vec<String> = args
        .iter()
        .map(|a| a.to_string_lossy().into_owned())
        .collect();
    match main_plugin(&args) {
        Ok(code) => Ok(code),
        Err(e) => match e.downcast::<Exit>() {
            Ok(exit) => {
                eprint!("{}", exit.message);
                Ok(ExitCode::from(exit.code))
            }
            Err(e) => Err(e),
        },
    }
}

/// The number of built-in plugins (for `vtools`' usage).
pub(crate) fn count() -> usize {
    PLUGINS.len()
}

/// Prints what glibc's `getopt` prints, with `prog` as `argv[0]`, for an option a plugin does
/// not know, before the plugin prints its usage.
pub(crate) fn getopt_message(prog: &str, e: &clap::Error) {
    use clap::error::{ContextKind, ContextValue, ErrorKind};
    if e.kind() != ErrorKind::UnknownArgument {
        return;
    }
    let Some(ContextValue::String(arg)) = e.get(ContextKind::InvalidArg) else {
        return;
    };
    if arg.starts_with("--") {
        eprintln!("{prog}: unrecognized option '{arg}'");
    } else if let Some(c) = arg.chars().nth(1) {
        eprintln!("{prog}: invalid option -- '{c}'");
    }
}

/// An option with values for [`in_order`]: its argument ID, its values (an `Append`
/// argument), and the event of each value.
pub(crate) type Valued<'a, T> = (&'a str, &'a [String], fn(&'a str) -> T);

/// The occurrences of the given options in command-line order, as getopt reads them: one
/// event per value of each [`Valued`] option, and one per occurrence of each
/// `(id, count, event)` (arguments with `Count`).
pub(crate) fn in_order<'a, T: Copy>(
    matches: &ArgMatches,
    valued: &[Valued<'a, T>],
    flags: &[(&str, u8, T)],
) -> Vec<T> {
    let mut events: Vec<(usize, T)> = Vec::new();
    for &(id, values, make) in valued {
        if let Some(idx) = matches.indices_of(id) {
            events.extend(idx.zip(values).map(|(i, v)| (i, make(v))));
        }
    }
    // A count has a default value, with an index of its own: only real occurrences count.
    for &(id, n, event) in flags {
        if n > 0
            && let Some(idx) = matches.indices_of(id)
        {
            events.extend(idx.map(|i| (i, event)));
        }
    }
    events.sort_by_key(|e| e.0);
    events.into_iter().map(|e| e.1).collect()
}

/// A shared option that bcftools' `main_plugin` and `+split-vep` check as getopt reads it.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Checked<'a> {
    /// `-O`.
    OutputType(&'a str),
    /// `-i` or `-e`: a second one is an error.
    Filter,
    /// `--regions-overlap`.
    RegionsOverlap(&'a str),
    /// `--targets-overlap`.
    TargetsOverlap(&'a str),
    /// `-W`.
    WriteIndex(&'a str),
    /// `-v`, checked by the caller.
    Verbosity(&'a str),
}

/// Checks one shared option the way bcftools' option loops do; `filters` counts `-i/-e`.
pub(crate) fn check_shared(event: Checked<'_>, filters: &mut usize) -> anyhow::Result<()> {
    match event {
        Checked::OutputType(t) => {
            output::requested_type(&[OutputSetting::Type(t.to_owned())])?;
        }
        Checked::Filter => {
            *filters += 1;
            if *filters > 1 {
                bail!(
                    "Error: only one -i or -e expression can be given, and they cannot be combined"
                );
            }
        }
        Checked::RegionsOverlap(o) => {
            parse_overlap(Some(o), Overlap::Record, "regions-overlap")?;
        }
        Checked::TargetsOverlap(o) => {
            parse_overlap(Some(o), Overlap::Pos, "targets-overlap")?;
        }
        Checked::WriteIndex(w) => {
            output::parse_write_index(w)?;
        }
        Checked::Verbosity(_) => {}
    }
    Ok(())
}

/// Finds a built-in plugin by name.
fn find(name: &str) -> Option<&'static Info> {
    PLUGINS.iter().find(|p| p.name == name)
}

/// What `-l` lists.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Listing {
    /// The plugin names.
    Names,
    /// The names and descriptions (`-l -v`).
    Descriptions,
}

/// Prints the plugin names, or their descriptions with `-v` (bcftools' `list_plugins`).
fn list_plugins(listing: Listing) -> ExitCode {
    let mut out = String::new();
    for p in &PLUGINS {
        match listing {
            Listing::Descriptions => out.push_str(&format!("\n-- {} --\n{}", p.name, p.about)),
            Listing::Names => {
                out.push_str(p.name);
                out.push('\n');
            }
        }
    }
    if listing == Listing::Descriptions {
        out.push('\n');
    }
    print!("{out}");
    ExitCode::SUCCESS
}

/// Checks the shared options in command-line order, as `main_plugin`'s getopt loop does, and
/// returns the verbosity from `-v` (one level each) and `--verbosity=INT`.
fn check_general(g: &General, matches: &ArgMatches) -> anyhow::Result<i32> {
    let events = in_order(
        matches,
        &[
            ("verbosity", &g.verbosity, Checked::Verbosity),
            ("output_type", &g.output_type, Checked::OutputType),
            ("include", &g.include, |_| Checked::Filter),
            ("exclude", &g.exclude, |_| Checked::Filter),
            (
                "regions_overlap",
                &g.regions_overlap,
                Checked::RegionsOverlap,
            ),
            (
                "targets_overlap",
                &g.targets_overlap,
                Checked::TargetsOverlap,
            ),
            ("write_index", &g.write_index, Checked::WriteIndex),
        ],
        &[],
    );
    let mut verbose = 0;
    let mut filters = 0;
    for event in events {
        match event {
            Checked::Verbosity("") => verbose += 1,
            Checked::Verbosity(v) => {
                verbose = match v.parse::<i32>() {
                    Ok(n) if n >= 0 => n,
                    _ => bail!("Could not parse argument: --verbosity {v}"),
                };
                apply_verbosity(v)?;
            }
            other => check_shared(other, &mut filters)?,
        }
    }
    Ok(verbose)
}

fn main_plugin(args: &[String]) -> anyhow::Result<ExitCode> {
    let Some(first) = args.first() else {
        return Err(general_usage());
    };
    let mut info = None;
    let mut rest = args;
    if !first.starts_with('-') {
        let Some(p) = find(first) else {
            return Err(Exit::error(format!(
                "\nThe vtools plugin \"{first}\" was not found or is not functional. The built-in plugins are counts, fill-tags, setGT, and split-vep.\n\nCould not load \"{first}\".\n\n"
            )));
        };
        if let Entry::Run(run) = p.entry {
            return run(&args[1..]);
        }
        info = Some(p);
        rest = &args[1..];
    }

    // bcftools' loop handles a bad option as `-h` (`usage_only`) and goes on reading the
    // options, so `-l` and `-V` still win over it; getopt names the plugin in its messages.
    let prog = info.map_or("plugin", |p| p.name);
    let (rest, bad_option) = general_args(prog, rest);
    let parsed = General::command()
        .try_get_matches_from(rest)
        .and_then(|m| Ok((General::from_arg_matches(&m)?, m)));
    let general = match parsed {
        Ok(g) => Some(g),
        Err(e) => {
            // getopt reports a bad option and bcftools prints the plugin's usage.
            let Some(p) = info else {
                return Err(general_usage());
            };
            getopt_message(p.name, &e);
            None
        }
    };
    let Some((mut general, matches)) = general else {
        let p = info.expect("a plugin was named");
        return Err(Exit {
            code: 0,
            message: (p.usage)(),
        }
        .into());
    };
    general.help |= bad_option;
    let verbose = check_general(&general, &matches)?;
    if general.list {
        let listing = if verbose > 0 {
            Listing::Descriptions
        } else {
            Listing::Names
        };
        return Ok(list_plugins(listing));
    }
    let Some(info) = info else {
        return Err(general_usage());
    };
    let Entry::Init(init) = info.entry else {
        unreachable!("plugins that take over the command line return above");
    };
    if general.version {
        let v = env!("CARGO_PKG_VERSION");
        println!("vtools  {v} using pasta {v}");
        println!("plugin at {v} using pasta {v}\n");
        return Ok(ExitCode::SUCCESS);
    }
    if general.help {
        return Err(Exit {
            code: 0,
            message: (info.usage)(),
        }
        .into());
    }
    let lists = Lists::new(&matches, &general);
    Driver::new(general, lists, init)?.run()
}

/// The regions (`-r`/`-R`) and targets (`-t`/`-T`) to read, each the later of its inline list
/// and its file.
pub(crate) struct Lists {
    pub regions: Option<(String, ListSource)>,
    pub targets: Option<(String, ListSource)>,
}

impl Lists {
    fn new(matches: &ArgMatches, g: &General) -> Self {
        Self {
            regions: later(
                matches,
                ("regions", g.regions.as_deref()),
                ("regions_file", g.regions_file.as_deref()),
            )
            .map(|(spec, source)| (spec.to_owned(), source)),
            targets: later(
                matches,
                ("targets", g.targets.as_deref()),
                ("targets_file", g.targets_file.as_deref()),
            )
            .map(|(spec, source)| (spec.to_owned(), source)),
        }
    }
}

/// Reads records through a plugin and writes the ones it keeps (the rest of `main_plugin`).
struct Driver {
    general: General,
    lists: Lists,
    init: InitFn,
    /// The input, `-` for standard input.
    fname: String,
    /// The arguments for the plugin's `init`.
    plugin_args: Vec<String>,
}

impl Driver {
    /// Splits the input name from the plugin's arguments as getopt leaves them: the first
    /// non-option argument is the input, unless the plugin's options come first.
    fn new(mut general: General, lists: Lists, init: InitFn) -> anyhow::Result<Self> {
        let mut rest = std::mem::take(&mut general.inputs);
        rest.append(&mut general.plugin_args);
        let starts_with_option = rest
            .first()
            .is_none_or(|a| a.starts_with('-') && a.len() > 1);
        // The plugin's arguments start with getopt's `argv[0]`: the input, or the `--` before
        // the plugin's options.
        let fname = if starts_with_option {
            if io::stdin().is_terminal() {
                if rest.is_empty() {
                    return Err(general_usage());
                }
                // bcftools calls the plugin's init without an input, for `-- -h` and
                // `+fill-tags -- -l`.
                rest.insert(0, "--".to_owned());
                init(&rest, None)?;
                return Err(general_usage());
            }
            rest.insert(0, "--".to_owned());
            "-".to_owned()
        } else {
            rest[0].clone()
        };
        Ok(Self {
            general,
            lists,
            init,
            fname,
            plugin_args: rest,
        })
    }

    fn run(self) -> anyhow::Result<ExitCode> {
        let g = &self.general;
        // `strtol(optarg, 0, 0)`; HTSlib runs a count below one without threads.
        let threads = output::threads(g.threads.as_deref()).unwrap_or(0);
        let options = synced_options(
            &self.lists,
            g.regions_overlap.last().map(String::as_str),
            g.targets_overlap.last().map(String::as_str),
            threads,
        )?;
        let mut reader = SyncedReader::open(&self.fname, options).map_err(|e| {
            let name = if self.fname == "-" {
                "standard input"
            } else {
                &self.fname
            };
            anyhow::anyhow!("Failed to read from {name}: {e}")
        })?;

        // init_data
        let mut hdr_out = reader.header().dup()?;
        let Loaded {
            mut plugin,
            output,
            mut max_unpack,
        } = (self.init)(&self.plugin_args, Some((reader.header(), &mut hdr_out)))?;
        let filter = match (g.include.first(), g.exclude.first()) {
            (Some(e), _) => Some((Filter::init(reader.header(), e)?, FilterLogic::Include)),
            (None, Some(e)) => Some((Filter::init(reader.header(), e)?, FilterLogic::Exclude)),
            (None, None) => None,
        };
        if let Some((f, _)) = &filter
            && !max_unpack.is_empty()
        {
            let needs = f.max_unpack();
            max_unpack = if needs.is_empty() {
                Unpack::empty()
            } else {
                max_unpack | needs
            };
        }
        if !g.no_version {
            version::append_version(&mut hdr_out, "plugin")?;
        }
        let out_name = g
            .output
            .as_deref()
            .map_or_else(|| "-".to_owned(), |p| p.display().to_string());
        let mut sink = match output {
            Output::Records => Some(Sink::open(g, &hdr_out, &out_name, threads)?),
            Output::None => None,
        };
        // The reader's lookahead records are parsed as far as the plugin needs, with or
        // without worker threads.
        reader.set_max_unpack(max_unpack);
        let mut source = Source::new(&mut reader, threads, max_unpack);
        let result = process_records(
            &mut reader,
            &mut source,
            filter,
            plugin.as_mut(),
            &hdr_out,
            sink.as_mut(),
        );
        let error = match result {
            Ok(()) => reader
                .take_pending_error()
                .map(|e| anyhow::anyhow!("Error: {e}")),
            Err(Fatal::Read(e)) => Some(anyhow::anyhow!("Error: {e}")),
            Err(Fatal::Other(e)) => Some(e),
        };
        if let Some(e) = error {
            if let Some(s) = sink {
                s.abandon(&hdr_out);
            }
            return Err(e);
        }
        plugin.destroy()?;
        if let Some(s) = sink {
            s.finish(&hdr_out, &out_name)?;
        }
        Ok(ExitCode::SUCCESS)
    }
}

fn parse_overlap(s: Option<&str>, default: Overlap, what: &str) -> anyhow::Result<Overlap> {
    match s {
        None => Ok(default),
        Some(s) => {
            Overlap::parse(s).ok_or_else(|| anyhow::anyhow!("Could not parse: --{what} {s}"))
        }
    }
}

/// The reader options of `-r/-R`, `-t/-T` (with `^`), and their overlap modes.
pub(crate) fn synced_options(
    lists: &Lists,
    regions_overlap: Option<&str>,
    targets_overlap: Option<&str>,
    threads: usize,
) -> anyhow::Result<SyncedOptions> {
    let regions_overlap = parse_overlap(regions_overlap, Overlap::Record, "regions-overlap")?;
    let targets_overlap = parse_overlap(targets_overlap, Overlap::Pos, "targets-overlap")?;
    let mut options = SyncedOptions {
        threads,
        ..SyncedOptions::default()
    };
    if let Some((spec, source)) = &lists.regions {
        let mut r = Regions::new(spec.as_bytes(), *source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = regions_overlap;
        options.regions = Some(r);
    }
    if let Some((spec, source)) = &lists.targets {
        let body = match spec.strip_prefix('^') {
            Some(rest) => {
                options.targets_exclude = true;
                rest
            }
            None => spec.as_str(),
        };
        let mut t = Regions::new(body.as_bytes(), *source)
            .map_err(|_| anyhow::anyhow!("Failed to read the targets: {spec}"))?;
        t.overlap = targets_overlap;
        options.targets = Some(t);
    }
    Ok(options)
}

/// A read error (reported as `Error: ...`) or any other fatal error.
enum Fatal {
    Read(SyncedError),
    Other(anyhow::Error),
}

impl From<anyhow::Error> for Fatal {
    fn from(e: anyhow::Error) -> Self {
        Self::Other(e)
    }
}

/// The records of the input, parsed on worker threads with `--threads` when the input is a
/// plain VCF stream.
pub(crate) struct Source {
    parallel: Option<ParallelVcf>,
    pool: Vec<Record>,
    max_unpack: Unpack,
}

impl Source {
    /// Starts reading, with `threads` parsing workers when they help; `max_unpack` limits
    /// how much of each line they parse, as the reader's own limit does.
    pub fn new(reader: &mut SyncedReader, threads: usize, max_unpack: Unpack) -> Self {
        let mut parallel = ParallelVcf::new(reader, threads);
        if let Some(p) = &mut parallel {
            p.set_max_unpack(max_unpack);
        }
        Self {
            parallel,
            pool: Vec::new(),
            max_unpack,
        }
    }

    /// Reads the next record into `rec`; `false` at the end (`bcf_sr_next_line`).
    ///
    /// # Errors
    ///
    /// Returns the reader's errors.
    pub fn next(
        &mut self,
        reader: &mut SyncedReader,
        rec: &mut Record,
    ) -> Result<bool, SyncedError> {
        match &mut self.parallel {
            Some(p) => p.next(reader, rec, &mut self.pool),
            None => {
                // The record goes into the reader's lookahead buffer, to be parsed later.
                rec.set_max_unpack(self.max_unpack);
                reader.next_record(rec)
            }
        }
    }

    /// Returns a record for reuse.
    pub fn recycle(&mut self, rec: Record) {
        if self.pool.len() < 4096 {
            self.pool.push(rec);
        }
    }

    /// Takes a spare record.
    pub fn take(&mut self) -> Record {
        self.pool.pop().unwrap_or_default()
    }
}

/// The VCF/BCF output: the writer and, with `--threads`, the formatting workers.
pub(crate) struct Sink {
    pub writer: Writer,
    emitter: Emitter,
}

impl Sink {
    /// Opens the output of `-o/-O` and writes the header (bcftools' `init_data`).
    fn open(g: &General, hdr_out: &Header, out_name: &str, threads: usize) -> anyhow::Result<Self> {
        let settings: Vec<OutputSetting> = g
            .output_type
            .iter()
            .map(|t| OutputSetting::Type(t.clone()))
            .collect();
        let sink = Self::create(
            &settings,
            g.output.as_deref(),
            g.write_index.last().map(String::as_str),
            hdr_out,
            out_name,
            threads,
        )?;
        Ok(sink.with_write_error(format!("[main_plugin] Error: cannot write to {out_name}")))
    }

    /// Opens an output and writes the header.
    pub fn create(
        settings: &[OutputSetting],
        output: Option<&std::path::Path>,
        write_index: Option<&str>,
        hdr_out: &Header,
        out_name: &str,
        threads: usize,
    ) -> anyhow::Result<Self> {
        // The plugins' output name defaults to "-", which the level error names.
        let (out_type, level) =
            output::resolve(settings, Some(output.unwrap_or(std::path::Path::new("-"))))?;
        let out_path = output.map(output::data_path);
        let mut writer = match Writer::create(
            out_path.as_deref(),
            &output::write_options(out_type, level, threads),
        ) {
            Ok(w) => w,
            Err(e) => {
                log::error(
                    "hts_open_format",
                    format_args!("Failed to open file \"{out_name}\" : {}", strerror(&e)),
                );
                bail!("Can't write to \"{out_name}\": {}", strerror(&e))
            }
        };
        if writer.write_header(hdr_out).is_err() {
            bail!("[init_data] Error: cannot write to {out_name}");
        }
        if let Some(wi) = write_index {
            let fmt = output::parse_write_index(wi)?;
            let ok = output::index_request(fmt, out_type, output)
                .is_some_and(|req| writer.request_index(req).is_ok());
            if !ok {
                bail!("Error: failed to initialise index for {out_name}");
            }
        }
        writer.start_index(hdr_out)?;
        let emitter = Emitter::new(hdr_out, &writer, threads, out_name.to_owned());
        Ok(Self { writer, emitter })
    }

    /// Reports `message` when a record cannot be written.
    #[must_use]
    pub fn with_write_error(mut self, message: String) -> Self {
        self.emitter = self.emitter.with_write_error(message);
        self
    }

    /// Writes a record (or queues it for the workers); returns it when it can be reused.
    pub fn write(&mut self, hdr_out: &Header, rec: Record) -> anyhow::Result<Option<Record>> {
        self.emitter.write(hdr_out, &mut self.writer, rec)
    }

    /// Takes a record the workers are done with.
    pub fn spare(&mut self) -> Option<Record> {
        self.emitter.spare()
    }

    /// Drops the output after a fatal error without flushing it, as bcftools' `error()` exits;
    /// with threads, what a serial run would have written is written first.
    pub fn abandon(mut self, hdr_out: &Header) {
        self.emitter.abort(hdr_out, &mut self.writer);
        self.writer.abandon();
    }

    /// Writes everything queued, then closes the output and writes its index.
    pub fn finish(mut self, hdr_out: &Header, out_name: &str) -> anyhow::Result<()> {
        self.emitter.finish(hdr_out, &mut self.writer)?;
        match self.writer.finish() {
            Ok(()) => Ok(()),
            Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
            Err(e @ WriteError::Deferred(_)) => Err(self.emitter.write_failure(e)),
            Err(WriteError::Io(e)) => Err(output::write_failure(e, || {
                anyhow::anyhow!("[destroy_data] Error: close failed .. {out_name}")
            })),
            Err(e) => Err(anyhow::Error::new(e)
                .context(format!("[destroy_data] Error: close failed .. {out_name}"))),
        }
    }
}

/// The main loop of `main_plugin`.
fn process_records(
    reader: &mut SyncedReader,
    source: &mut Source,
    mut filter: Option<(Filter, FilterLogic)>,
    plugin: &mut dyn Plugin,
    hdr_out: &Header,
    mut sink: Option<&mut Sink>,
) -> Result<(), Fatal> {
    let mut rec = source.take();
    loop {
        match source.next(reader, &mut rec) {
            Ok(true) => {}
            Ok(false) => break,
            Err(e) => return Err(Fatal::Read(e)),
        }
        let hdr = reader.header();
        if let Some((f, logic)) = &mut filter {
            let pass = f.test(hdr, &mut rec)?;
            if pass != (*logic == FilterLogic::Include) {
                continue;
            }
        }
        if !plugin.process(hdr, hdr_out, &mut rec)? {
            continue;
        }
        let Some(sink) = sink.as_deref_mut() else {
            continue;
        };
        let errors = rec.errors();
        if !errors.is_empty() {
            return Err(Fatal::Other(anyhow::anyhow!(
                "[E::main_plugin] Unchecked error ({}), exiting",
                errors.bits()
            )));
        }
        let next = sink.spare().unwrap_or_else(|| source.take());
        let done = std::mem::replace(&mut rec, next);
        if let Some(back) = sink.write(hdr_out, done)? {
            source.recycle(back);
        }
    }
    if let Some(sink) = sink {
        sink.emitter.finish(hdr_out, &mut sink.writer)?;
    }
    Ok(())
}
