//! A `getopt_long` front end for the subcommands' clap arguments.
//!
//! bcftools parses its options with glibc's `getopt_long`, which differs from clap in several
//! ways that users rely on:
//!
//! - a long option may be abbreviated to any unique prefix (`--mask-w` for `--mask-with`);
//! - a required value is always the next argument, even when it starts with `-`
//!   (`--threads -1`, `-s -x`);
//! - an optional value (`N::`, `W::`) is taken only when attached (`-N1000`, `--verbosity=4`);
//! - short options cluster (`-Hv3`), and options and operands may come in any order (unless
//!   `POSIXLY_CORRECT` is set: then the first operand ends the options);
//! - an unknown option, a missing value, or an ambiguous prefix is reported as
//!   `PROG: unrecognized option '--x'` (and similar) and then handled by the command, which
//!   prints its usage.
//!
//! [`rewrite`] parses the command line of one subcommand that way, using the clap definition of
//! its arguments, and rewrites it to a form clap reads unambiguously: every option with its
//! canonical name and any value attached with `=`, then the operands (those after a `--` still
//! after `--`, where the plugin interface finds the plugin's own options). Clap's derive
//! arguments stay the single definition of each command's options.

use std::{
    ffi::{OsStr, OsString},
    os::unix::ffi::OsStrExt,
};

use clap::{ArgAction, Command};

/// What happens after a getopt error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OnError {
    /// The error ends the parse; the command prints its usage (bcftools' `case '?'`).
    Usage,
    /// The error is a warning and the option is ignored (bcftools' `tabix`, whose loop has no
    /// `case '?'`).
    Ignore,
}

/// Whether long options are recognized.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Style {
    /// `getopt_long`: `--name` is a long option.
    Long,
    /// Plain `getopt`: `--name` is a cluster of short options starting with `-`.
    Short,
}

/// The result of parsing a subcommand's command line.
#[derive(Debug, Eq, PartialEq)]
pub enum Parsed {
    /// The command line rewritten for clap, with the getopt warnings already met (only with
    /// [`OnError::Ignore`]).
    Args {
        /// The rewritten arguments, without the program and subcommand names.
        args: Vec<OsString>,
        /// The getopt messages, each ending with a newline.
        warnings: Vec<String>,
    },
    /// A help option (one with clap's `ArgAction::Help`): the command prints its usage.
    Help,
    /// A getopt error (with [`OnError::Usage`]); the message ends with a newline.
    Error(String),
}

/// What an option takes, as in getopt's option string (`x`, `x:`, `x::`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Takes {
    Nothing,
    Required,
    Optional,
}

/// One option of the subcommand.
#[derive(Debug)]
struct Opt {
    /// The name clap knows it by: `--long`, or `-c` without a long name.
    canonical: String,
    takes: Takes,
    help: bool,
}

/// The bcftools commands with a `--help` long option. vtools has one on every command; on the
/// others it is matched only in full, so that `--he` still finds `--header-only` in `view`.
const LONG_HELP: &[&str] = &[
    "call", "csq", "gtcheck", "isec", "merge", "norm", "plugin", "query", "sort", "stats",
];

/// The subcommand's options, by short and long name (aliases included).
#[derive(Debug)]
struct Spec {
    opts: Vec<Opt>,
    shorts: Vec<(char, usize)>,
    longs: Vec<(String, usize)>,
    /// Whether `--help` may be abbreviated (when bcftools has it).
    help_prefixes: bool,
}

/// A long name looked up as getopt does.
enum Long<'a> {
    Found(usize, &'a str),
    Ambiguous(Vec<&'a str>),
    Unknown,
}

impl Spec {
    /// Reads the options of subcommand `name`, `cmd`, which must be built (`Command::build`)
    /// so that each argument's number of values is known.
    fn new(cmd: &Command, name: &str) -> Self {
        let mut spec = Self {
            opts: Vec::new(),
            shorts: Vec::new(),
            longs: Vec::new(),
            help_prefixes: LONG_HELP.contains(&name),
        };
        for arg in cmd.get_arguments().filter(|a| !a.is_positional()) {
            let action = arg.get_action();
            let takes = if !action.takes_values() {
                Takes::Nothing
            } else if arg.is_require_equals_set()
                || arg.get_num_args().is_some_and(|n| n.min_values() == 0)
            {
                Takes::Optional
            } else {
                Takes::Required
            };
            let canonical = match (arg.get_long(), arg.get_short()) {
                (Some(l), _) => format!("--{l}"),
                (None, Some(c)) => format!("-{c}"),
                (None, None) => continue,
            };
            let i = spec.opts.len();
            spec.opts.push(Opt {
                canonical,
                takes,
                help: matches!(
                    action,
                    ArgAction::Help | ArgAction::HelpShort | ArgAction::HelpLong
                ),
            });
            let shorts = arg.get_short().into_iter();
            spec.shorts.extend(
                shorts
                    .chain(arg.get_all_short_aliases().unwrap_or_default())
                    .map(|c| (c, i)),
            );
            let longs = arg.get_long().into_iter();
            spec.longs.extend(
                longs
                    .chain(arg.get_all_aliases().unwrap_or_default())
                    .map(|l| (l.to_owned(), i)),
            );
        }
        spec
    }

    fn short(&self, c: char) -> Option<usize> {
        self.shorts.iter().find(|(s, _)| *s == c).map(|&(_, i)| i)
    }

    /// Finds a long option by its exact name or a unique prefix (`getopt_long`). Prefixes of
    /// several names of one option are not ambiguous.
    fn long(&self, name: &str) -> Long<'_> {
        if let Some((l, i)) = self.longs.iter().find(|(l, _)| l == name) {
            return Long::Found(*i, l);
        }
        let mut found: Option<(usize, &str)> = None;
        let mut names = Vec::new();
        let mut ambiguous = false;
        let candidates = self
            .longs
            .iter()
            .filter(|(l, _)| (self.help_prefixes || l != "help") && l.starts_with(name));
        for (l, i) in candidates {
            names.push(l.as_str());
            match found {
                None => found = Some((*i, l)),
                Some((j, _)) if j != *i => ambiguous = true,
                Some(_) => {}
            }
        }
        match found {
            _ if ambiguous => Long::Ambiguous(names),
            Some((i, l)) => Long::Found(i, l),
            None => Long::Unknown,
        }
    }
}

/// Parses the arguments of subcommand `cmd` (built, see [`Spec::new`]) as bcftools' getopt
/// loop does. `prog` is the subcommand name, as in bcftools' messages (its `argv[0]`).
pub fn rewrite(
    cmd: &Command,
    prog: &str,
    args: &[OsString],
    style: Style,
    on_error: OnError,
) -> Parsed {
    // glibc stops at the first operand when POSIXLY_CORRECT is set.
    let order = if std::env::var_os("POSIXLY_CORRECT").is_some() {
        Order::Posix
    } else {
        Order::Permute
    };
    parse(cmd, prog, args, style, on_error, order)
}

/// Where getopt looks for options.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Order {
    /// Anywhere: operands are moved after the options (glibc's default).
    Permute,
    /// Up to the first operand, which ends the options (`POSIXLY_CORRECT`).
    Posix,
}

/// [`rewrite`] with the order getopt reads options in.
fn parse(
    cmd: &Command,
    prog: &str,
    args: &[OsString],
    style: Style,
    on_error: OnError,
    order: Order,
) -> Parsed {
    let spec = Spec::new(cmd, prog);
    let mut out = Vec::with_capacity(args.len() + 1);
    let mut operands = Vec::new();
    // The arguments after `--`, or from the first operand on with `Order::Posix`.
    let mut rest: Option<&[OsString]> = None;
    let mut warnings = Vec::new();
    let mut i = 0;
    // Reports a getopt error: the parse ends, or the option is ignored.
    macro_rules! fail {
        ($($msg:tt)*) => {{
            let msg = format!("{prog}: {}\n", format_args!($($msg)*));
            match on_error {
                OnError::Usage => return Parsed::Error(msg),
                OnError::Ignore => warnings.push(msg),
            }
        }};
    }
    while i < args.len() {
        let arg = &args[i];
        i += 1;
        let bytes = arg.as_encoded_bytes();
        if bytes == b"--" {
            rest = Some(&args[i..]);
            break;
        }
        if bytes.len() < 2 || bytes[0] != b'-' {
            if order == Order::Posix {
                rest = Some(&args[i - 1..]);
                break;
            }
            operands.push(arg.clone());
            continue;
        }
        if style == Style::Long
            && let Some(body) = bytes.strip_prefix(b"--")
        {
            let (name, value) = match body.iter().position(|&c| c == b'=') {
                Some(eq) => (&body[..eq], Some(OsStr::from_bytes(&body[eq + 1..]))),
                None => (body, None),
            };
            let name = String::from_utf8_lossy(name);
            let (idx, full) = match spec.long(&name) {
                Long::Found(idx, full) => (idx, full),
                Long::Ambiguous(names) => {
                    let list: String = names.iter().map(|n| format!(" '--{n}'")).collect();
                    // getopt names the whole argument, `=VALUE` included.
                    fail!(
                        "option '{}' is ambiguous; possibilities:{list}",
                        arg.to_string_lossy()
                    );
                    continue;
                }
                Long::Unknown => {
                    fail!("unrecognized option '{}'", arg.to_string_lossy());
                    continue;
                }
            };
            let opt = &spec.opts[idx];
            match (opt.takes, value) {
                (Takes::Nothing, Some(_)) => {
                    fail!("option '--{full}' doesn't allow an argument");
                }
                _ if opt.help => return Parsed::Help,
                (Takes::Nothing, None) | (Takes::Optional, None) => {
                    out.push(OsString::from(&opt.canonical));
                }
                (_, Some(v)) => out.push(with_value(&opt.canonical, v)),
                (Takes::Required, None) => match args.get(i) {
                    Some(v) => {
                        i += 1;
                        out.push(with_value(&opt.canonical, v));
                    }
                    None => fail!("option '--{full}' requires an argument"),
                },
            }
            continue;
        }
        // A cluster of short options, one byte each as getopt reads them.
        let cluster = &bytes[1..];
        for (pos, &c) in cluster.iter().enumerate() {
            let Some(idx) = spec.short(char::from(c)).filter(|_| c.is_ascii()) else {
                fail!("invalid option -- '{}'", String::from_utf8_lossy(&[c]));
                continue;
            };
            let opt = &spec.opts[idx];
            if opt.help {
                return Parsed::Help;
            }
            let rest = &cluster[pos + 1..];
            match opt.takes {
                Takes::Nothing => {
                    out.push(OsString::from(&opt.canonical));
                    continue;
                }
                Takes::Optional if rest.is_empty() => out.push(OsString::from(&opt.canonical)),
                _ if !rest.is_empty() => {
                    out.push(with_value(&opt.canonical, OsStr::from_bytes(rest)));
                }
                _ => match args.get(i) {
                    Some(v) => {
                        i += 1;
                        out.push(with_value(&opt.canonical, v));
                    }
                    None => fail!("option requires an argument -- '{}'", char::from(c)),
                },
            }
            break;
        }
    }
    // Before `--`, an operand is `-` or does not start with `-`, so clap reads it as one.
    out.append(&mut operands);
    if let Some(rest) = rest {
        out.push(OsString::from("--"));
        out.extend(rest.iter().cloned());
    }
    Parsed::Args {
        args: out,
        warnings,
    }
}

/// [`rewrite`] for the plugin interface and the plugins, whose parsers take `String`s without a
/// program name: the arguments as clap is to read them, or the getopt message on an error (empty
/// for a help option). `prog` is the name in the messages.
pub fn rewrite_strings(cmd: &Command, prog: &str, args: &[String]) -> Result<Vec<String>, String> {
    let mut cmd = cmd.clone();
    cmd.build();
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    match rewrite(&cmd, prog, &args, Style::Long, OnError::Usage) {
        // Every argument came from a `String`.
        Parsed::Args { args, .. } => Ok(args
            .into_iter()
            .map(|a| a.to_string_lossy().into_owned())
            .collect()),
        Parsed::Help => Err(String::new()),
        Parsed::Error(msg) => Err(msg),
    }
}

/// `--long=VALUE` or `-c=VALUE`: clap takes an `=`-attached value as is, even when it starts
/// with `-` (and drops just the one `=` after a short option).
fn with_value(canonical: &str, value: &OsStr) -> OsString {
    let mut s = OsString::with_capacity(canonical.len() + 1 + value.len());
    s.push(canonical);
    s.push("=");
    s.push(value);
    s
}

#[cfg(test)]
mod tests {
    use clap::Arg;

    use super::*;

    fn norm() -> Command {
        let mut cmd = Command::new("norm")
            .disable_help_flag(true)
            .arg(
                Arg::new("no_realign")
                    .short('N')
                    .long("no-realign")
                    .num_args(0..=1)
                    .require_equals(true),
            )
            .arg(Arg::new("m").short('m').long("multiallelics"))
            .arg(Arg::new("mask").long("mask"))
            .arg(Arg::new("mask_with").long("mask-with"))
            .arg(Arg::new("o").short('o'))
            .arg(Arg::new("v").short('v').action(ArgAction::SetTrue))
            .arg(Arg::new("H").short('H').action(ArgAction::Count))
            .arg(
                Arg::new("min_mq")
                    .long("min-MQ")
                    .visible_alias("min-mq")
                    .short('q'),
            )
            .arg(Arg::new("help").long("help").action(ArgAction::Help))
            .arg(Arg::new("files").num_args(0..));
        cmd.build();
        cmd
    }

    fn args(v: &[&str]) -> Vec<OsString> {
        v.iter().map(OsString::from).collect()
    }

    fn rewritten(v: &[&str]) -> Vec<OsString> {
        match rewrite(&norm(), "norm", &args(v), Style::Long, OnError::Usage) {
            Parsed::Args { args, warnings } => {
                assert!(warnings.is_empty(), "{warnings:?}");
                args
            }
            other => panic!("{v:?} parsed as {other:?}"),
        }
    }

    fn error(v: &[&str]) -> String {
        match rewrite(&norm(), "norm", &args(v), Style::Long, OnError::Usage) {
            Parsed::Error(msg) => msg,
            other => panic!("{v:?} parsed as {other:?}"),
        }
    }

    #[test]
    fn optional_values_are_attached_only() {
        assert_eq!(
            rewritten(&["-N1000", "-N", "x", "--no-realign=5", "--no-realign"]),
            args(&[
                "--no-realign=1000",
                "--no-realign",
                "--no-realign=5",
                "--no-realign",
                "x"
            ])
        );
    }

    #[test]
    fn required_values_take_the_next_argument() {
        assert_eq!(
            rewritten(&["-m", "-", "-o", "-x", "--mask", "--x", "-m-any", "-o="]),
            args(&[
                "--multiallelics=-",
                "-o=-x",
                "--mask=--x",
                "--multiallelics=-any",
                "-o=="
            ])
        );
    }

    #[test]
    fn clusters_and_operands() {
        assert_eq!(
            rewritten(&["a.vcf", "-vHHo", "out", "-", "-Hm+", "--", "-v", "b"]),
            args(&[
                "-v",
                "-H",
                "-H",
                "-o=out",
                "-H",
                "--multiallelics=+",
                "a.vcf",
                "-",
                "--",
                "-v",
                "b"
            ])
        );
    }

    #[test]
    fn long_prefixes() {
        assert_eq!(
            rewritten(&[
                "--mask-w",
                "f",
                "--multi=both",
                "--min-mq",
                "5",
                "--min-M=6"
            ]),
            args(&[
                "--mask-with=f",
                "--multiallelics=both",
                "--min-MQ=5",
                "--min-MQ=6"
            ])
        );
        // An exact name wins over a longer one it prefixes.
        assert_eq!(rewritten(&["--mask", "f"]), args(&["--mask=f"]));
        // Both names of --min-MQ: not ambiguous.
        assert_eq!(rewritten(&["--min", "5"]), args(&["--min-MQ=5"]));
    }

    #[test]
    fn errors() {
        assert_eq!(
            error(&["--ma", "x"]),
            "norm: option '--ma' is ambiguous; possibilities: '--mask' '--mask-with'\n"
        );
        assert_eq!(
            error(&["--nosuch=1"]),
            "norm: unrecognized option '--nosuch=1'\n"
        );
        assert_eq!(error(&["-vZ"]), "norm: invalid option -- 'Z'\n");
        assert_eq!(error(&["-o"]), "norm: option requires an argument -- 'o'\n");
        assert_eq!(
            error(&["--mask"]),
            "norm: option '--mask' requires an argument\n"
        );
        assert_eq!(
            error(&["--help=x"]),
            "norm: option '--help' doesn't allow an argument\n"
        );
    }

    #[test]
    fn values_keep_their_bytes() {
        use std::os::unix::ffi::OsStringExt;

        let raw = |b: &[u8]| OsString::from_vec(b.to_vec());
        let mut v = vec![
            raw(b"-o\xffx"),
            raw(b"--mask=\xfe"),
            raw(b"-m"),
            raw(b"\xfd"),
        ];
        v.push(raw(b"f\xfc"));
        let Parsed::Args { args, .. } = rewrite(&norm(), "norm", &v, Style::Long, OnError::Usage)
        else {
            panic!("{v:?} did not parse");
        };
        let want = [
            &b"-o=\xffx"[..],
            b"--mask=\xfe",
            b"--multiallelics=\xfd",
            b"f\xfc",
        ];
        let got: Vec<&[u8]> = args.iter().map(|a| a.as_encoded_bytes()).collect();
        assert_eq!(got, want);
        // A byte that is not ASCII is not an option.
        assert_eq!(error(&["-v\u{e9}"]), "norm: invalid option -- '\u{fffd}'\n");
    }

    #[test]
    fn posixly_correct_stops_at_the_first_operand() {
        let parsed = parse(
            &norm(),
            "norm",
            &args(&["-v", "a.vcf", "-H", "--", "b"]),
            Style::Long,
            OnError::Usage,
            Order::Posix,
        );
        assert_eq!(
            parsed,
            Parsed::Args {
                args: args(&["-v", "--", "a.vcf", "-H", "--", "b"]),
                warnings: Vec::new(),
            }
        );
    }

    #[test]
    fn errors_as_warnings_go_on() {
        // bcftools' plugin interface and tabix read on after a bad option.
        let parsed = rewrite(
            &norm(),
            "plugin",
            &args(&["--nosuch", "-v", "--ma=x", "-Zo", "out", "--mask"]),
            Style::Long,
            OnError::Ignore,
        );
        assert_eq!(
            parsed,
            Parsed::Args {
                args: args(&["-v", "-o=out"]),
                warnings: vec![
                    "plugin: unrecognized option '--nosuch'\n".to_owned(),
                    "plugin: option '--ma=x' is ambiguous; possibilities: '--mask' '--mask-with'\n"
                        .to_owned(),
                    "plugin: invalid option -- 'Z'\n".to_owned(),
                    "plugin: option '--mask' requires an argument\n".to_owned(),
                ],
            }
        );
    }

    #[test]
    fn help() {
        let parse = |v: &[&str]| rewrite(&norm(), "norm", &args(v), Style::Long, OnError::Usage);
        assert_eq!(parse(&["-v", "--help"]), Parsed::Help);
        assert_eq!(parse(&["--he"]), Parsed::Help);
        // Where bcftools has no --help, it is matched only in full.
        assert_eq!(
            rewrite(
                &norm(),
                "view",
                &args(&["--he"]),
                Style::Long,
                OnError::Usage
            ),
            Parsed::Error("view: unrecognized option '--he'\n".to_owned())
        );
        // The first problem wins.
        assert!(matches!(parse(&["-Z", "--help"]), Parsed::Error(_)));
    }

    /// One command's options in bcftools 1.24 (`bcftools-1.24.tsv`, made from the C sources).
    struct Bcftools {
        name: &'static str,
        /// getopt's option string.
        optstring: &'static str,
        /// The long options: name, argument, and the value getopt returns (options with the
        /// same argument and value are one option to getopt).
        longs: Vec<(&'static str, Takes, &'static str)>,
    }

    fn bcftools() -> Vec<Bcftools> {
        let takes = |k: &str| match k {
            "0" => Takes::Nothing,
            "1" => Takes::Required,
            _ => Takes::Optional,
        };
        include_str!("shortopt/bcftools-1.24.tsv")
            .lines()
            .map(|line| {
                let mut fields = line.split('\t');
                let name = fields.next().unwrap_or_default();
                let optstring = fields.next().unwrap_or_default();
                let longs = fields
                    .next()
                    .unwrap_or_default()
                    .split(',')
                    .filter(|e| !e.is_empty())
                    .map(|e| {
                        let mut parts = e.split('=');
                        let (long, k, val) = (parts.next(), parts.next(), parts.next());
                        let (Some(long), Some(k), Some(val)) = (long, k, val) else {
                            panic!("{name}: bad long option {e:?}");
                        };
                        (long, takes(k), val)
                    })
                    .collect();
                Bcftools {
                    name,
                    optstring,
                    longs,
                }
            })
            .collect()
    }

    /// The vtools subcommand's options.
    fn vtools_spec(name: &str) -> Spec {
        use clap::CommandFactory;

        let mut sub = crate::Cli::command()
            .find_subcommand(name)
            .unwrap_or_else(|| panic!("no subcommand {name}"))
            .clone();
        sub.build();
        Spec::new(&sub, name)
    }

    /// Lists how each bcftools option differs in vtools: missing, or taking a different kind of
    /// argument.
    fn audit() -> Vec<String> {
        let mut found = Vec::new();
        for bt in bcftools() {
            let name = bt.name;
            let spec = vtools_spec(name);
            let chars: Vec<char> = bt.optstring.chars().collect();
            let mut i = 0;
            while i < chars.len() {
                let c = chars[i];
                let colons = chars[i + 1..].iter().take_while(|&&x| x == ':').count();
                i += 1 + colons;
                let want = [Takes::Nothing, Takes::Required, Takes::Optional][colons.min(2)];
                match spec.short(c) {
                    None => found.push(format!("{name} -{c}: missing")),
                    Some(j) if spec.opts[j].takes != want => found.push(format!(
                        "{name} -{c}: {:?}, bcftools {want:?}",
                        spec.opts[j].takes
                    )),
                    Some(_) => {}
                }
            }
            for &(long, want, _) in &bt.longs {
                match spec.longs.iter().find(|(l, _)| l == long) {
                    None => found.push(format!("{name} --{long}: missing")),
                    Some(&(_, j)) if spec.opts[j].takes != want => found.push(format!(
                        "{name} --{long}: {:?}, bcftools {want:?}",
                        spec.opts[j].takes
                    )),
                    Some(_) => {}
                }
            }
        }
        found
    }

    #[test]
    fn every_bcftools_option_is_defined() {
        // Differences that end the same way. getopt reports an unknown option as '?', which
        // these commands' `case '?'` (or `default`) handles as their usage, so a missing `-?`
        // or a dead `-1` makes no difference; concat's `-h` prints the usage with or without a
        // value. query's --disable-automatic-newline takes a value in bcftools by mistake (its
        // -N takes none); vtools takes none either way (README query).
        const SAME: &[&str] = &[
            "filter -?: missing",
            "sort -?: missing",
            "concat -?: missing",
            "consensus -?: missing",
            "csq -?: missing",
            "consensus -1: missing",
            "concat -h: Optional, bcftools Required",
            "query --disable-automatic-newline: Nothing, bcftools Required",
        ];
        let found: Vec<String> = audit()
            .into_iter()
            .filter(|f| !SAME.contains(&f.as_str()))
            .collect();
        assert!(found.is_empty(), "{}", found.join("\n"));
    }

    /// glibc's `getopt_long` lookup of `--name` in a bcftools table: the long option found, or
    /// `Err(true)` when ambiguous and `Err(false)` when unknown.
    fn bcftools_long<'a>(bt: &'a Bcftools, name: &str) -> Result<&'a str, bool> {
        if let Some(&(l, _, _)) = bt.longs.iter().find(|(l, _, _)| *l == name) {
            return Ok(l);
        }
        let mut found: Option<(&str, Takes, &str)> = None;
        for &(l, takes, val) in bt.longs.iter().filter(|(l, _, _)| l.starts_with(name)) {
            match found {
                None => found = Some((l, takes, val)),
                Some((_, t, v)) if t != takes || v != val => return Err(true),
                Some(_) => {}
            }
        }
        found.map(|(l, _, _)| l).ok_or(false)
    }

    #[test]
    fn long_prefixes_resolve_as_in_bcftools() {
        let mut found = Vec::new();
        for bt in bcftools().iter().filter(|bt| bt.name != "tabix") {
            let spec = vtools_spec(bt.name);
            let names = bt
                .longs
                .iter()
                .map(|(l, _, _)| *l)
                .chain(spec.longs.iter().map(|(l, _)| l.as_str()));
            let mut prefixes: Vec<&str> = names
                .flat_map(|l| l.char_indices().map(move |(i, c)| &l[..i + c.len_utf8()]))
                .collect();
            prefixes.sort_unstable();
            prefixes.dedup();
            for p in prefixes {
                let vt = spec.long(p);
                let same = match (bcftools_long(bt, p), &vt) {
                    // Each ends in the usage.
                    (Ok("help") | Err(_), Long::Unknown | Long::Ambiguous(_)) => true,
                    (Ok("help"), Long::Found(i, l)) => spec.opts[*i].help || *l == "help",
                    (Ok(l), Long::Found(i, _)) => spec.longs.contains(&(l.to_owned(), *i)),
                    // A vtools extension.
                    (Err(false), Long::Found(_, l)) => !bt.longs.iter().any(|(b, _, _)| b == l),
                    _ => false,
                };
                // vtools' concat also takes --output-file (README command lines).
                if !(same || bt.name == "concat" && p == "output-") {
                    let vt = match vt {
                        Long::Found(_, l) => format!("--{l}"),
                        Long::Ambiguous(_) => "ambiguous".to_owned(),
                        Long::Unknown => "unknown".to_owned(),
                    };
                    found.push(format!(
                        "{} --{p}: bcftools {:?}, vtools {vt}",
                        bt.name,
                        bcftools_long(bt, p)
                    ));
                }
            }
        }
        assert!(found.is_empty(), "{}", found.join("\n"));
    }

    #[test]
    fn plain_getopt_ignores_errors() {
        let parsed = rewrite(
            &norm(),
            "tabix",
            &args(&["--no-such", "f", "-Z"]),
            Style::Short,
            OnError::Ignore,
        );
        assert_eq!(
            parsed,
            Parsed::Args {
                args: args(&["-o=-such", "f"]),
                warnings: vec![
                    "tabix: invalid option -- '-'\n".to_owned(),
                    "tabix: invalid option -- 'n'\n".to_owned(),
                    "tabix: invalid option -- 'Z'\n".to_owned(),
                ],
            }
        );
    }
}
