//! `vtools isec`: intersections, unions, and complements of VCF/BCF files (bcftools isec).

use std::{
    env,
    fs::{self, File},
    io::{self, BufWriter, Write},
    os::unix::fs::DirBuilderExt,
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::bail;
use clap::ArgMatches;
use pasta::variant::{
    Header, Record, Unpack,
    io::{IndexRequest, OutputType, WriteError, Writer},
    readlist::{self, ListSource},
    synced::{IndexRequirement, Overlap, Pairing, SyncedReaders},
};

use crate::options::later;
use crate::{
    commands::head::apply_verbosity,
    filter::{self, FilterLogic, RecordFilter},
    output::{self, OutputSetting, WriteIndex},
    raw_input::strerror,
    version,
};

const USAGE: &str = "
About:   Create intersections, unions and complements of VCF files.
Usage:   vtools isec [options] <A.vcf.gz> <B.vcf.gz> [...]

Options:
    -c, --collapse STRING          Treat as identical records with <snps|indels|both|all|some|none|id>, see man page for details [none]
    -C, --complement               Output positions present only in the first file but missing in the others
    -e, --exclude EXPR             Exclude sites for which the expression is true
    -f, --apply-filters LIST       Require at least one of the listed FILTER strings (e.g. \"PASS,.\")
    -i, --include EXPR             Include only sites for which the expression is true
    -l, --file-list FILE           Read the input file names from the file
        --no-version               Do not append version and command line to the header
    -n, --nfiles [+-=~]INT         Output positions present in this many (=), this many or more (+), this many or fewer (-), the exact (~) files
    -o, --output FILE              Write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -p, --prefix DIR               If given, subset each of the input files accordingly, see also -w
    -r, --regions REGION           Restrict to comma-separated list of regions
    -R, --regions-file FILE        Restrict to regions listed in a file
        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -t, --targets REGION           Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT              Use multithreading with INT worker threads [0]
    -v, --verbosity INT            Verbosity level
    -w, --write LIST               List of files to write with -p given as 1-based indexes. By default, all files are written
    -W, --write-index[=FMT]        Automatically index the output files [off]

Examples:
   # Create intersection and complements of two sets saving the output in dir/*
   vtools isec A.vcf.gz B.vcf.gz -p dir

   # Filter sites in A and B (but not in C) and create intersection
   vtools isec -e'MAF<0.01' -i'dbSNP=1' -e - A.vcf.gz B.vcf.gz C.vcf.gz -p dir

   # Extract and write records from A shared by both A and B using exact allele match
   vtools isec A.vcf.gz B.vcf.gz -p dir -n =2 -w 1

   # Extract and write records from C found in A and C but not in B
   vtools isec A.vcf.gz B.vcf.gz C.vcf.gz -p dir -n~101 -w 3

   # Extract records private to A or B comparing by position only
   vtools isec A.vcf.gz B.vcf.gz -p dir -n -1 -c all

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF files, bgzipped and indexed.
    inputs: Vec<String>,

    /// Treat as identical records with <snps|indels|both|all|some|none|id> [none].
    #[arg(
        short = 'c',
        long = "collapse",
        value_name = "STRING",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    collapse: Vec<String>,

    /// Output positions present only in the first file but missing in the others.
    #[arg(short = 'C', long = "complement", action = clap::ArgAction::Count)]
    complement: u8,

    /// Exclude sites for which the expression is true (one, or one per file; - for none).
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    exclude: Vec<String>,

    /// Require at least one of the listed FILTER strings (e.g. "PASS,.").
    #[arg(
        short = 'f',
        long = "apply-filters",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    apply_filters: Option<String>,

    /// Include only sites for which the expression is true (one, or one per file; - for none).
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    include: Vec<String>,

    /// Read the input file names from the file.
    #[arg(
        short = 'l',
        long = "file-list",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    file_list: Option<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Output positions present in this many (=), this many or more (+), this many or fewer
    /// (-), the exact (~) files.
    #[arg(
        short = 'n',
        long = "nfiles",
        value_name = "[+-=~]INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    nfiles: Vec<String>,

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

    /// If given, subset each of the input files accordingly, see also -w.
    #[arg(
        short = 'p',
        long = "prefix",
        value_name = "DIR",
        allow_hyphen_values = true
    )]
    prefix: Option<String>,

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
        allow_hyphen_values = true
    )]
    regions_overlap: Option<String>,

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
    #[arg(
        long = "targets-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true
    )]
    targets_overlap: Option<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    verbosity: Option<String>,

    /// List of files to write with -p given as 1-based indexes. By default, all files are
    /// written.
    #[arg(
        short = 'w',
        long = "write",
        value_name = "LIST",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    write: Vec<String>,

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

    /// Print help.
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Count)]
    help: u8,
}

/// Which positions are output (bcftools' `OP_*`).
#[derive(Clone, Debug, Eq, PartialEq)]
enum Op {
    /// `-n +N`: present in N or more files.
    Plus(i64),
    /// `-n -N`: present in N or fewer files.
    Minus(i64),
    /// `-n =N` or `-n N`: present in exactly N files.
    Equal(i64),
    /// Two files without `-n`: the four-way split written with `-p`.
    Venn,
    /// `-C`: present only in the first file.
    Complement,
    /// `-n ~MASK`: present in exactly the files of the bitmask.
    Exact(Vec<bool>),
}

/// glibc's `sscanf(s, "%d")`: leading white space, a sign, and digits read as `strtol` reads
/// them (saturating at the `long` range), then truncated to `int`; `None` when no digits.
fn scan_int(s: &str) -> Option<i32> {
    let t = s.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
    let (neg, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let end = digits.bytes().take_while(u8::is_ascii_digit).count();
    if end == 0 {
        return None;
    }
    let v = digits[..end].bytes().fold(0i64, |v, d| {
        let d = i64::from(d - b'0');
        if neg {
            v.saturating_mul(10).saturating_sub(d)
        } else {
            v.saturating_mul(10).saturating_add(d)
        }
    });
    #[expect(
        clippy::cast_possible_truncation,
        reason = "glibc stores the long in an int without a range check"
    )]
    Some(v as i32)
}

/// Parses `-n` (bcftools' option loop).
fn parse_nfiles(s: &str) -> anyhow::Result<Op> {
    let bad = || anyhow::anyhow!("Could not parse --nfiles {s}");
    let (kind, rest) = match s.as_bytes().first() {
        Some(b'-') => ('-', &s[1..]),
        Some(b'+') => ('+', &s[1..]),
        Some(b'=') => ('=', &s[1..]),
        Some(b'~') => ('~', &s[1..]),
        Some(c) if c.is_ascii_digit() => ('=', s),
        _ => return Err(bad()),
    };
    if kind == '~' {
        return Ok(Op::Exact(rest.bytes().map(|b| b == b'1').collect()));
    }
    let n = i64::from(scan_int(rest).ok_or_else(bad)?);
    Ok(match kind {
        '-' => Op::Minus(n),
        '+' => Op::Plus(n),
        _ => Op::Equal(n),
    })
}

/// Returns the values of an option with their command-line positions.
fn indexed<'a>(matches: &ArgMatches, id: &str, values: &'a [String]) -> Vec<(usize, &'a str)> {
    matches
        .indices_of(id)
        .map(|idx| idx.zip(values.iter().map(String::as_str)).collect())
        .unwrap_or_default()
}

fn parse_overlap(s: Option<&str>, default: Overlap, what: &str) -> anyhow::Result<Overlap> {
    match s {
        None => Ok(default),
        Some(s) => {
            Overlap::parse(s).ok_or_else(|| anyhow::anyhow!("Could not parse: --{what} {s}"))
        }
    }
}

/// A file being written: its name, writer, and the index it should get.
struct OutFile {
    name: String,
    writer: Writer,
    index: IndexPlan,
}

/// The index bcftools builds for a `-p` output after closing it.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum IndexPlan {
    /// None.
    None,
    /// A TBI index built while writing (bcftools' `tbx_index_build`, whose failure it ignores).
    Tbi,
    /// A CSI index built while writing (bcftools' `bcf_index_build`: "Could not index" when
    /// it fails).
    Csi,
    /// Requested for an output that cannot be indexed: an error after closing.
    Impossible,
}

/// A file's filter expression and whether it includes or excludes.
type BoxedFilter = (Box<dyn RecordFilter + Send>, FilterLogic);

/// The state of bcftools' `args_t` that the main loop needs.
struct Isec {
    sr: SyncedReaders,
    op: Op,
    filters: Vec<Option<BoxedFilter>>,
    /// The files to write with `-p` (all when `None`).
    write: Option<Vec<bool>>,
    iwrite: usize,
    nwrite: usize,
    prefix: Option<String>,
    /// Whether the records of one file go to standard output or `-o` (`out_std`).
    out_std: bool,
    std_writer: Option<Writer>,
    out_name: String,
    /// The index file of the standard output or `-o` output, when `-W` asks for one.
    out_index_name: String,
    sites: Option<Box<dyn Write>>,
    sites_name: String,
    /// The `-p` outputs, by bcftools' file number (four for the two-file split).
    outs: Vec<Option<OutFile>>,
}

/// `-C` and `-n`, checked in command-line order as bcftools' option loop does.
fn parse_op(matches: &ArgMatches, args: &Args) -> anyhow::Result<Option<Op>> {
    let mut op: Option<Op> = None;
    let mut ops: Vec<(usize, Option<&str>)> = indexed(matches, "nfiles", &args.nfiles)
        .into_iter()
        .map(|(i, v)| (i, Some(v)))
        .collect();
    if args.complement > 0
        && let Some(idx) = matches.indices_of("complement")
    {
        ops.extend(idx.map(|i| (i, None)));
    }
    ops.sort_unstable_by_key(|(i, _)| *i);
    for (_, v) in ops {
        match v {
            None => {
                if op.as_ref().is_some_and(|o| *o != Op::Complement) {
                    bail!("Error: either -C or -n should be given, not both.");
                }
                op = Some(Op::Complement);
            }
            Some(v) => {
                if op == Some(Op::Complement) {
                    bail!("Error: either -C or -n should be given, not both.");
                }
                if op.is_some() {
                    bail!("Error: -n should be given only once.");
                }
                op = Some(parse_nfiles(v)?);
            }
        }
    }
    Ok(op)
}

/// `-c`: which records count as the same.
fn parse_collapse(values: &[String]) -> anyhow::Result<Pairing> {
    let mut pairing = Pairing::empty();
    for c in values {
        pairing = match c.as_str() {
            "snps" => pairing | Pairing::SNPS,
            "indels" => pairing | Pairing::INDELS,
            "both" => pairing | Pairing::BOTH,
            "any" | "all" => pairing | Pairing::ANY,
            "some" => pairing | Pairing::SOME,
            "none" => Pairing::empty(),
            "id" => pairing | Pairing::ID,
            _ => bail!("The --collapse string \"{c}\" not recognised."),
        };
    }
    Ok(pairing)
}

/// `-i/-e` in command-line order; "-" stands for no filter on that file.
fn ordered_filters(matches: &ArgMatches, args: &Args) -> Vec<(Option<String>, FilterLogic)> {
    let mut flt: Vec<(usize, Option<String>, FilterLogic)> = Vec::new();
    for (i, e) in indexed(matches, "include", &args.include) {
        flt.push((i, (e != "-").then(|| e.to_owned()), FilterLogic::Include));
    }
    for (i, e) in indexed(matches, "exclude", &args.exclude) {
        flt.push((i, (e != "-").then(|| e.to_owned()), FilterLogic::Exclude));
    }
    flt.sort_by_key(|(i, ..)| *i);
    flt.into_iter().map(|(_, e, l)| (e, l)).collect()
}

pub fn run(mut args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    if args.help > 0 {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }

    let mut op = parse_op(matches, &args)?;
    if args.write.len() > 1 {
        bail!("The option -w accepts a list of indices and can be given only once");
    }
    let pairing = parse_collapse(&args.collapse)?;
    let flt = ordered_filters(matches, &args);

    let regions_overlap = parse_overlap(
        args.regions_overlap.as_deref(),
        Overlap::Record,
        "regions-overlap",
    )?;
    let targets_overlap = parse_overlap(
        args.targets_overlap.as_deref(),
        Overlap::Pos,
        "targets-overlap",
    )?;
    let threads = isec_threads(&args);
    let write_index = args
        .write_index
        .as_deref()
        .map(output::parse_write_index)
        .transpose()?;
    let mut settings = Vec::new();
    for (_, t) in indexed(matches, "output_type", &args.output_type) {
        settings.push(OutputSetting::Type(t.to_owned()));
        // Validate each -O as bcftools does while parsing.
        output::requested_type(&settings)?;
    }

    if args.inputs.is_empty() && args.file_list.is_none() {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
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

    let mut sr = SyncedReaders::new();
    sr.set_pairing(pairing);
    sr.set_apply_filters(args.apply_filters.as_deref().map(str::as_bytes));
    sr.set_threads(threads);
    let targets = later(
        matches,
        ("targets", args.targets.as_deref()),
        ("targets_file", args.targets_file.as_deref()),
    );
    if let Some((spec, source)) = targets {
        sr.set_targets_overlap(targets_overlap);
        if sr.set_targets_spec(spec.as_bytes(), source, None).is_err() {
            bail!("Failed to read the targets: {spec}");
        }
    }
    let regions = later(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    );
    if let Some((spec, source)) = regions {
        sr.set_regions_overlap(regions_overlap);
        if sr.set_regions_spec(spec.as_bytes(), source).is_err() {
            bail!("Failed to read the regions: {spec}");
        }
    }
    if files.len() == 2 && op.is_none() {
        op = Some(Op::Venn);
        if args.prefix.is_none() {
            bail!("Expected the -p option");
        }
    }
    let op = op.unwrap_or(Op::Plus(1));
    sr.set_index_requirement(IndexRequirement::Required);
    sr.set_shared_threads(threads, files.len());
    for f in &files {
        if let Err(e) = sr.add_reader(f) {
            bail!("Failed to open {f}: {e}");
        }
    }

    let mut isec = init(&args, sr, op, flt, write_index, &settings, targets)?;
    isec.run(&args, write_index, &settings)?;
    Ok(ExitCode::SUCCESS)
}

/// The filters of the files: one `-i/-e` for all of them, or one each.
fn init_filters(
    sr: &SyncedReaders,
    flt: Vec<(Option<String>, FilterLogic)>,
) -> anyhow::Result<Vec<Option<BoxedFilter>>> {
    let nreaders = sr.reader_count();
    let mut filters = Vec::new();
    if flt.is_empty() {
        return Ok(filters);
    }
    if flt.len() > 1 && flt.len() != nreaders {
        bail!("Error: expected either one -i/-e option or as many as there are input files");
    }
    let exprs: Vec<(Option<String>, FilterLogic)> = if flt.len() < nreaders {
        if flt[0].0.is_none() {
            bail!("Error: useless use of -i/-e");
        }
        vec![flt[0].clone(); nreaders]
    } else {
        flt
    };
    for (i, (expr, logic)) in exprs.into_iter().enumerate() {
        filters.push(match expr {
            Some(e) => Some((filter::init(sr.header(i), &e)?, logic)),
            None => None,
        });
    }
    Ok(filters)
}

/// `-w`: the files to write (all when `None`), the last one listed, and how many are listed.
fn parse_write(list: &str, nreaders: usize) -> anyhow::Result<(Option<Vec<bool>>, usize, usize)> {
    let mut write: Option<Vec<bool>> = None;
    let mut iwrite = 0;
    let mut nwrite = 0;
    let mut p = list;
    while !p.is_empty() {
        let w = write.get_or_insert_with(|| vec![false; nreaders]);
        let Some(i) = scan_int(p) else {
            bail!("Could not parse --write {list}");
        };
        let Some(k) = usize::try_from(i)
            .ok()
            .filter(|&k| (1..=nreaders).contains(&k))
        else {
            bail!("The index is out of range: {i} (-w {list})");
        };
        w[k - 1] = true;
        iwrite = k - 1;
        nwrite += 1;
        p = p.find(',').map_or("", |c| &p[c + 1..]);
    }
    Ok((write, iwrite, nwrite))
}

/// bcftools' `init_data`: filters, the bitmask, `-w`, and the `-p` directory with its README.
/// `targets` are the `-t/-T` targets, with which the records of a single file go to the
/// output.
fn init(
    args: &Args,
    sr: SyncedReaders,
    op: Op,
    flt: Vec<(Option<String>, FilterLogic)>,
    write_index: Option<WriteIndex>,
    settings: &[OutputSetting],
    targets: Option<(&str, ListSource)>,
) -> anyhow::Result<Isec> {
    let nreaders = sr.reader_count();
    let filters = init_filters(&sr, flt)?;

    if let Op::Exact(mask) = &op {
        let text = args
            .nfiles
            .first()
            .map(|s| s[1..].to_owned())
            .unwrap_or_default();
        if mask.len() != nreaders {
            bail!("The number of files does not match the bitmask: {nreaders} vs {text}");
        }
        if text.bytes().any(|b| b != b'0' && b != b'1') {
            bail!("Unexpected bitmask: {text}");
        }
    }

    let (write, iwrite, nwrite) = match args.write.first() {
        Some(list) => parse_write(list, nreaders)?,
        None => (None, 0, 0),
    };
    if nwrite > 1 && args.prefix.is_none() {
        bail!(
            "Expected -p when multiple output files given: --write {}",
            args.write[0]
        );
    }
    if op == Op::Complement && nwrite > 0 {
        if nwrite > 1 {
            bail!("Multiple files to -w make no sense with -C");
        }
        if !write.as_ref().is_some_and(|w| w[0]) {
            bail!("Only -w1 makes sense with -C");
        }
    }

    let mut isec = Isec {
        sr,
        op,
        filters,
        write,
        iwrite,
        nwrite,
        prefix: args.prefix.clone(),
        out_std: (nwrite == 1 && args.prefix.is_none()) || (targets.is_some() && nreaders == 1),
        std_writer: None,
        out_name: String::new(),
        out_index_name: String::new(),
        sites: None,
        sites_name: String::new(),
        outs: Vec::new(),
    };

    if let Some(prefix) = &args.prefix {
        isec.open_prefix_outputs(args, prefix, write_index, settings)?;
    } else {
        // bcftools opens -o for the sites list even when records go there instead.
        let sites: Box<dyn Write> = match &args.output {
            Some(path) => match File::create(path) {
                Ok(f) => Box::new(BufWriter::new(f)),
                Err(e) => bail!("Can't write to \"{}\": {}", path.display(), strerror(&e)),
            },
            None => Box::new(BufWriter::new(output::stdout())),
        };
        if !isec.out_std {
            isec.sites = Some(sites);
            isec.sites_name = args
                .output
                .as_deref()
                .map_or_else(|| "standard output".to_owned(), |p| p.display().to_string());
        }
    }
    Ok(isec)
}

impl Isec {
    /// Opens the `-p` outputs, the README listing them, and the sites list (bcftools'
    /// `init_data`).
    fn open_prefix_outputs(
        &mut self,
        args: &Args,
        prefix: &str,
        write_index: Option<WriteIndex>,
        settings: &[OutputSetting],
    ) -> anyhow::Result<()> {
        let nreaders = self.sr.reader_count();
        let readme_name = format!("{prefix}/README.txt");
        mkdir_p(&readme_name)?;
        let mut log = match File::create(&readme_name) {
            Ok(f) => BufWriter::new(f),
            Err(e) => bail!("{readme_name}: {}", strerror(&e)),
        };
        let mut cmd = String::from("vtools isec ");
        for a in env::args().skip(2) {
            cmd.push(' ');
            cmd.push_str(&a);
        }
        write!(
            log,
            "This file was produced by vcfisec.\nThe command line was:\t{cmd}\n\nUsing the following file names:\n"
        )?;
        let base_type = output::requested_type(settings)?;
        let suffix = match base_type {
            OutputType::Bcf | OutputType::BcfUncompressed => "bcf",
            OutputType::VcfBgzf => "vcf.gz",
            OutputType::Vcf => "vcf",
        };
        let plan = |ty: OutputType| {
            let is_tbi = write_index.is_none_or(|w| w == WriteIndex::Tbi);
            let wants = ty == OutputType::VcfBgzf || ty == OutputType::Bcf || !is_tbi;
            match (wants, ty.is_bgzf()) {
                (false, _) => (IndexPlan::None, 0),
                (true, true) if ty == OutputType::VcfBgzf && is_tbi => (IndexPlan::Tbi, 0),
                (true, true) => (IndexPlan::Csi, 14),
                (true, false) => (IndexPlan::Impossible, 0),
            }
        };
        let nfiles = if self.op == Op::Venn { 4 } else { nreaders };
        self.outs = (0..nfiles).map(|_| None).collect();
        for (i, j, legend) in self.prefix_outputs() {
            let name = format!("{prefix}/{i:04}.{suffix}");
            mkdir_p(&name)?;
            let (ty, level) = output::resolve(settings, Some(Path::new(&name)))?;
            let mut writer = match Writer::create(
                Some(Path::new(&name)),
                &output::write_options(ty, level, isec_threads(args)),
            ) {
                Ok(w) => w,
                Err(_) => bail!("Could not open {name}"),
            };
            if !args.no_version {
                version::append_version(self.sr.header_mut(j), "isec")?;
            }
            let header = self.sr.header(j);
            if writer.write_header(header).is_err() {
                bail!("[init_data] Error: cannot write to {name}");
            }
            let (index, min_shift) = plan(base_type);
            if matches!(index, IndexPlan::Tbi | IndexPlan::Csi) {
                let ext = if min_shift == 0 { "tbi" } else { "csi" };
                writer.request_index(IndexRequest {
                    path: PathBuf::from(format!("{name}.{ext}")),
                    min_shift,
                })?;
                writer.start_index(header)?;
            }
            writeln!(log, "{name}\t{legend}")?;
            self.outs[i] = Some(OutFile {
                name,
                writer,
                index,
            });
        }
        log.flush()?;
        let sites_name = format!("{prefix}/sites.txt");
        mkdir_p(&sites_name)?;
        match File::create(&sites_name) {
            Ok(f) => self.sites = Some(Box::new(BufWriter::new(f))),
            Err(e) => bail!("{sites_name}: {}", strerror(&e)),
        }
        self.sites_name = sites_name;
        Ok(())
    }

    /// The `-p` outputs to open: their number, the reader whose header they take, and their
    /// README legend.
    fn prefix_outputs(&self) -> Vec<(usize, usize, String)> {
        let fname = |i: usize| self.sr.file_name(i).to_owned();
        let wants = |i: usize| self.write.as_ref().is_none_or(|w| w[i]);
        let mut open: Vec<(usize, usize, String)> = Vec::new();
        if self.op == Op::Venn {
            let (a, b) = (fname(0), fname(1));
            if wants(0) {
                open.push((0, 0, format!("for records private to\t{a}")));
            }
            if wants(1) {
                open.push((1, 1, format!("for records private to\t{b}")));
            }
            if wants(0) {
                open.push((
                    2,
                    0,
                    format!("for records from {a} shared by both\t{a} {b}"),
                ));
            }
            if wants(1) {
                open.push((
                    3,
                    1,
                    format!("for records from {b} shared by both\t{a} {b}"),
                ));
            }
        } else {
            for i in 0..self.sr.reader_count() {
                if !wants(i) {
                    continue;
                }
                if self.op == Op::Complement && i > 0 {
                    break;
                }
                open.push((i, i, format!("for stripped\t{}", fname(i))));
            }
        }
        open
    }
}

/// `--threads` as bcftools reads it (`strtol(optarg, 0, 0)`); a negative count means none.
fn isec_threads(args: &Args) -> usize {
    crate::output::threads(args.threads.as_deref()).unwrap_or(0)
}

/// Creates the directories of a file's path, one level at a time (bcftools' `mkdir_p`).
fn mkdir_p(path: &str) -> anyhow::Result<()> {
    let bytes = path.as_bytes();
    let mut p = 1;
    while p < bytes.len() {
        while p < bytes.len() && bytes[p] != b'/' {
            p += 1;
        }
        if p == bytes.len() {
            break;
        }
        let dir = &path[..p];
        if let Err(e) = fs::DirBuilder::new().mode(0o775).create(dir)
            && e.kind() != io::ErrorKind::AlreadyExists
        {
            bail!("Error creating directory {dir}: {}", strerror(&e));
        }
        while p < bytes.len() && bytes[p] == b'/' {
            p += 1;
        }
    }
    Ok(())
}

impl Isec {
    /// bcftools' `isec_vcf`, then the closing part of `destroy_data`.
    fn run(
        &mut self,
        args: &Args,
        write_index: Option<WriteIndex>,
        settings: &[OutputSetting],
    ) -> anyhow::Result<()> {
        let nreaders = self.sr.reader_count();
        let out_std = self.out_std;
        if out_std {
            self.open_out_std(args, write_index, settings)?;
        }
        if self.nwrite == 0 && !out_std && self.prefix.is_none() {
            eprintln!("Note: -w option not given, printing list of sites...");
        }
        // VCF lines of records that are never written need no FORMAT parsing (pairing reads
        // up to INFO/END), unless a filter may look at them.
        for i in 0..nreaders {
            let written = if out_std {
                i == self.iwrite
            } else {
                self.prefix.is_some()
                    && self.write.as_ref().is_none_or(|w| w[i])
                    && (self.op != Op::Complement || i == 0)
            };
            if !written && self.filters.get(i).is_none_or(Option::is_none) {
                self.sr.set_reader_max_unpack(i, Unpack::SHR);
            }
        }

        let result = self.main_loop(nreaders);
        if let Err(e) = result {
            self.abandon();
            return Err(e);
        }
        if let Some(e) = self.sr.take_pending_error() {
            self.abandon();
            bail!("Error: {e}");
        }
        if let Some(w) = self.std_writer.take() {
            match w.finish() {
                Ok(()) => {}
                Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
                Err(WriteError::SaveIndex(_)) => {
                    bail!("Error: cannot write to index {}", self.out_index_name)
                }
                Err(e @ WriteError::Deferred(_)) => {
                    let name = &self.out_name;
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[isec_vcf] Error: cannot write to {name}")
                    }));
                }
                Err(e) => {
                    // bcftools names standard output "-" here.
                    let name = if self.out_name == "standard output" {
                        "-"
                    } else {
                        &self.out_name
                    };
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[isec_vcf] Error: close failed .. {name}")
                    }));
                }
            }
        }
        // destroy_data: close (and index) each -p output in turn.
        for k in 0..self.outs.len() {
            let Some(OutFile {
                name,
                writer,
                index,
            }) = self.outs[k].take()
            else {
                continue;
            };
            match (writer.finish(), index) {
                (Ok(()), IndexPlan::Impossible)
                | (Err(WriteError::SaveIndex(_)), IndexPlan::Csi) => {
                    self.abandon();
                    bail!("Could not index {name}");
                }
                (Ok(()), _) | (Err(WriteError::SaveIndex(_)), IndexPlan::Tbi) => {}
                (Err(WriteError::Deferred(_)), _) => {
                    self.abandon();
                    bail!("[isec_vcf] Error: cannot write");
                }
                (Err(_), _) => {
                    self.abandon();
                    bail!("[destroy_data] Error: close failed .. {name}");
                }
            }
        }
        if let Some(s) = &mut self.sites
            && let Err(e) = s.flush()
        {
            let name = &self.sites_name;
            return Err(output::write_failure(e, || {
                anyhow::anyhow!("[isec_vcf] Error: close failed .. {name}")
            }));
        }
        Ok(())
    }

    /// Opens the single output written to standard output or `-o` (`out_std`).
    fn open_out_std(
        &mut self,
        args: &Args,
        write_index: Option<WriteIndex>,
        settings: &[OutputSetting],
    ) -> anyhow::Result<()> {
        let (ty, level) = output::resolve(settings, args.output.as_deref())?;
        let out_path = args.output.as_deref().map(output::data_path);
        self.out_name = args
            .output
            .as_deref()
            .map_or_else(|| "standard output".to_owned(), |p| p.display().to_string());
        let mut writer = match Writer::create(
            out_path.as_deref(),
            &output::write_options(ty, level, isec_threads(args)),
        ) {
            Ok(w) => w,
            Err(e) => bail!("Can't write to {}: {}", self.out_name, strerror(&e)),
        };
        if !args.no_version {
            version::append_version(self.sr.header_mut(self.iwrite), "isec")?;
        }
        let header = self.sr.header(self.iwrite);
        if writer.write_header(header).is_err() {
            bail!("[isec_vcf] Error: cannot write to {}", self.out_name);
        }
        if let Some(fmt) = write_index {
            let ok = match output::index_request(fmt, ty, args.output.as_deref()) {
                Some(req) => {
                    self.out_index_name = req.path.display().to_string();
                    writer.request_index(req).is_ok()
                }
                None => false,
            };
            if !ok {
                bail!("Error: failed to initialise index for {}", self.out_name);
            }
            writer.start_index(header)?;
        }
        self.std_writer = Some(writer);
        Ok(())
    }

    /// Drops all output without flushing, as bcftools' `error()` exits.
    fn abandon(&mut self) {
        if let Some(w) = self.std_writer.take() {
            w.abandon();
        }
        for out in self.outs.iter_mut().filter_map(Option::take) {
            out.writer.abandon();
        }
        // stdio streams are flushed at exit.
        if let Some(s) = &mut self.sites {
            let _ = s.flush();
        }
    }

    /// Tests the per-file filters; returns how many files keep their record.
    fn apply_filters(&mut self, mut n: usize) -> anyhow::Result<usize> {
        for i in 0..self.sr.reader_count() {
            if !self.sr.has_line(i) {
                continue;
            }
            let Some(Some((flt, logic))) = self.filters.get_mut(i) else {
                continue;
            };
            let logic = *logic;
            let (rec, header) = self.sr.line_with_header(i).expect("the reader has a line");
            let mut pass = flt.pass(header, rec)?;
            if logic == FilterLogic::Exclude {
                pass = !pass;
            }
            if !pass {
                self.sr.clear_line(i);
                n -= 1;
            }
        }
        Ok(n)
    }

    fn main_loop(&mut self, nreaders: usize) -> anyhow::Result<()> {
        let mut site = Vec::new();
        loop {
            let n = match self.sr.next_line() {
                Ok(0) => break,
                Ok(n) => n,
                Err(e) => bail!("Error: {e}"),
            };
            let n = self.apply_filters(n)?;
            let Some(first) = (0..nreaders).find(|&i| self.sr.has_line(i)) else {
                continue; // filtered in all files
            };
            let n = n as i64;
            let keep = match &self.op {
                Op::Complement => n == 1 && self.sr.has_line(0),
                Op::Equal(k) => n == *k,
                Op::Plus(k) => n >= *k,
                Op::Minus(k) => n <= *k,
                Op::Exact(mask) => (0..nreaders).all(|i| self.sr.has_line(i) == mask[i]),
                Op::Venn => true,
            };
            if !keep {
                continue;
            }
            if let Some(w) = &mut self.std_writer {
                if let Some((rec, header)) = self.sr.line_with_header(self.iwrite)
                    && let Err(e) = w.write_record(header, rec)
                {
                    let name = &self.out_name;
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[isec_vcf] Error: cannot write to {name}")
                    }));
                }
                continue;
            }
            if let Some(sites) = &mut self.sites {
                site.clear();
                let (rec, header) = self
                    .sr
                    .line_with_header(first)
                    .expect("the reader has a line");
                site_line(&mut site, header, rec);
                site.push(b'\t');
                for i in 0..nreaders {
                    site.push(if self.sr.has_line(i) { b'1' } else { b'0' });
                }
                site.push(b'\n');
                if let Err(e) = sites.write_all(&site) {
                    let (len, name) = (site.len(), &self.sites_name);
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[isec_vcf] Error: failed to write {len} bytes to {name}")
                    }));
                }
            }
            if self.prefix.is_some() {
                self.write_prefix(nreaders)?;
            }
        }
        Ok(())
    }

    /// Writes the current records to the `-p` outputs.
    fn write_prefix(&mut self, nreaders: usize) -> anyhow::Result<()> {
        if self.op == Op::Venn && self.sr.has_line(0) && self.sr.has_line(1) {
            // Shared by both: each file's record goes to its own "shared" output.
            self.write_to(2, 0)?;
            return self.write_to(3, 1);
        }
        for i in 0..nreaders {
            if self.sr.has_line(i) {
                self.write_to(i, i)?;
            }
        }
        Ok(())
    }

    /// Writes reader `i`'s record to `-p` output `k`, unless `-w` leaves file `i` out.
    fn write_to(&mut self, k: usize, i: usize) -> anyhow::Result<()> {
        if self.write.as_ref().is_some_and(|w| !w[i]) {
            return Ok(());
        }
        let Some(out) = self.outs.get_mut(k).and_then(Option::as_mut) else {
            return Ok(());
        };
        let (rec, header) = self.sr.line_with_header(i).expect("the reader has a line");
        if out.writer.write_record(header, rec).is_err() {
            bail!("[isec_vcf] Error: cannot write");
        }
        Ok(())
    }
}

/// Formats CHROM, POS, REF, and ALT of a sites line.
fn site_line(out: &mut Vec<u8>, header: &Header, rec: &mut Record) {
    rec.unpack(Unpack::STR);
    out.extend_from_slice(header.contig_name(rec.rid()).unwrap_or_default());
    out.push(b'\t');
    write!(out, "{}", rec.pos() + 1).expect("writing to a Vec cannot fail");
    out.push(b'\t');
    let n = rec.allele_count();
    if n > 0 {
        out.extend_from_slice(rec.allele(0));
    } else {
        out.push(b'.');
    }
    out.push(b'\t');
    if n > 1 {
        out.extend_from_slice(rec.allele(1));
    } else {
        out.push(b'.');
    }
    for i in 2..n {
        out.push(b',');
        out.extend_from_slice(rec.allele(i));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nfiles_parse_as_bcftools() {
        assert_eq!(parse_nfiles("2").ok(), Some(Op::Equal(2)));
        assert_eq!(parse_nfiles("=2").ok(), Some(Op::Equal(2)));
        assert_eq!(parse_nfiles("+1").ok(), Some(Op::Plus(1)));
        assert_eq!(parse_nfiles("-1").ok(), Some(Op::Minus(1)));
        assert_eq!(
            parse_nfiles("~101").ok(),
            Some(Op::Exact(vec![true, false, true]))
        );
        assert!(parse_nfiles("x").is_err());
        assert!(parse_nfiles("+x").is_err());
        assert_eq!(scan_int(" 12abc"), Some(12));
    }
}
