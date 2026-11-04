//! `vtools annotate`: add or remove annotations (bcftools annotate, `vcfannotate.c`).
//!
//! Annotations come from a second VCF/BCF read in step with the input (`-a file.vcf.gz`), from
//! a tab-delimited file (`-a file.tab.gz -c CHROM,POS,...`, streamed with its tabix index, or
//! held in memory when it has END columns), or from the record itself (`-c NEW:=OLD`). The
//! columns (`-c`) decide what is transferred and how ([`columns`]); the setters
//! ([`setters`], [`format`]) copy values into each record. `-x` removes annotations, `-I`
//! sets IDs, `-m` marks sites, and `--rename-*` rename sequences and tags.

mod columns;
mod format;
mod line;
mod remove;
mod rename;
mod samples;
mod setters;
mod tab;
mod vcmp;

use std::{
    io::{self, IsTerminal},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, Record,
    index::split_idx_name,
    io::{WriteError, WriteOptions, Writer, reader::is_stream},
    readlist,
    regidx::{self, RegIdx},
    synced::{
        IndexRequirement, Overlap, Pairing, Regions, SyncedError, SyncedReader, SyncedReaders,
    },
};

use self::{
    columns::{AnnotKind, Column, ColumnsSpec},
    line::AnnotLine,
    remove::{Removals, RemoveSpec},
    rename::Rename,
    setters::{AlleleColumns, CopiedAlleles, Ctx, Scratch, Set, Src, lossy, site},
    tab::{Annotate, Ext, ExtFilter, TabAnnot, TabSource},
    vcmp::Vcmp,
};
use crate::options::{last_index, later_file_sticky};
use crate::{
    commands::{
        head::apply_verbosity,
        norm::{emit::Emitter, input::ParallelVcf},
    },
    convert::Convert,
    filter::{Filter, FilterLogic},
    output::{self, OutputSetting},
    raw_input::{Format, RawInput, report_exec_format_error},
    version,
};

/// bcftools' usage text (`usage()`), printed for `-?`, `--help`, and no input on a terminal.
const USAGE: &str = "
About:   Annotate and edit VCF/BCF files.
Usage:   vtools annotate [options] VCF

Options:
   -a, --annotations FILE          VCF file or tabix-indexed FILE with annotations: CHR\\tPOS[\\tVALUE]+
   -c, --columns LIST              List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details
   -C, --columns-file FILE         Read -c columns from FILE, one name per row, with optional --merge-logic TYPE: NAME[ TYPE]
   -e, --exclude EXPR              Exclude sites for which the expression is true (see man page for details)
       --force                     Continue despite parsing error (at your own risk!)
   -H, --header-line STR           Header line which should be appended to the VCF header, can be given multiple times
   -h, --header-lines FILE         Lines which should be appended to the VCF header
   -I, --set-id [+]FORMAT          Set ID column using a `vtools query`-like expression, see man page for details
   -i, --include EXPR              Select sites for which the expression is true (see man page for details)
   -k, --keep-sites                Leave -i/-e sites unchanged instead of discarding them
   -l, --merge-logic TAG:TYPE      Merge logic for multiple overlapping regions (see man page for details), EXPERIMENTAL
   -m, --mark-sites [+-]TAG        Add INFO/TAG flag to sites which are (\"+\") or are not (\"-\") listed in the -a file
       --min-overlap ANN:VCF       Required overlap as a fraction of variant in the -a file (ANN), the VCF (:VCF), or reciprocal (ANN:VCF)
       --no-version                Do not append version and command line to the header
   -o, --output FILE               Write output to a file [standard output]
   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
       --pair-logic STR            Matching records by <snps|indels|both|all|some|exact|id>, see man page for details [some]
   -r, --regions REGION            Restrict to comma-separated list of regions
   -R, --regions-file FILE         Restrict to regions listed in FILE
       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
       --rename-annots FILE        Rename annotations: TYPE/old\\tnew, where TYPE is one of FILTER,INFO,FORMAT
       --rename-chrs FILE          Rename sequences according to the mapping: old\\tnew
   -s, --samples [^]LIST           Comma separated list of samples to annotate (or exclude with \"^\" prefix)
   -S, --samples-file [^]FILE      File of samples to annotate (or exclude with \"^\" prefix)
       --single-overlaps           Keep memory low by avoiding complexities arising from handling multiple overlapping intervals
   -x, --remove LIST               List of annotations (e.g. ID,INFO/DP,FORMAT/DP,FILTER) to remove (or keep with \"^\" prefix). See man page for details
       --threads INT               Number of extra output compression threads [0]
   -v, --verbosity INT             Verbosity level
   -W, --write-index[=FMT]         Automatically index the output files [off]

Examples:
   http://samtools.github.io/bcftools/howtos/annotate.html

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin); further arguments are ignored, as in
    /// bcftools.
    inputs: Vec<String>,

    /// VCF file or tabix-indexed FILE with annotations: CHR\tPOS[\tVALUE]+.
    #[arg(short = 'a', long = "annotations", value_name = "FILE")]
    annotations: Option<String>,

    /// List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG.
    #[arg(short = 'c', long = "columns", value_name = "LIST", allow_hyphen_values = true, action = clap::ArgAction::Append)]
    columns: Vec<String>,

    /// Read -c columns from FILE, one name per row, with optional --merge-logic TYPE: NAME[ TYPE].
    #[arg(short = 'C', long = "columns-file", value_name = "FILE", action = clap::ArgAction::Append)]
    columns_file: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(short = 'e', long = "exclude", value_name = "EXPR", allow_hyphen_values = true, action = clap::ArgAction::Append)]
    exclude: Vec<String>,

    /// Continue despite parsing error (at your own risk!).
    #[arg(short = 'f', long = "force")]
    force: bool,

    /// Header line which should be appended to the VCF header, can be given multiple times.
    #[arg(short = 'H', long = "header-line", value_name = "STR", allow_hyphen_values = true, action = clap::ArgAction::Append)]
    header_line: Vec<String>,

    /// Lines which should be appended to the VCF header.
    #[arg(short = 'h', long = "header-lines", value_name = "FILE")]
    header_lines: Option<String>,

    /// Set ID column using a `bcftools query`-like expression.
    #[arg(
        short = 'I',
        long = "set-id",
        value_name = "[+]FORMAT",
        allow_hyphen_values = true
    )]
    set_id: Option<String>,

    /// Select sites for which the expression is true.
    #[arg(short = 'i', long = "include", value_name = "EXPR", allow_hyphen_values = true, action = clap::ArgAction::Append)]
    include: Vec<String>,

    /// Leave -i/-e sites unchanged instead of discarding them.
    #[arg(short = 'k', long = "keep-sites")]
    keep_sites: bool,

    /// Merge logic for multiple overlapping regions, TAG:TYPE.
    #[arg(short = 'l', long = "merge-logic", value_name = "TAG:TYPE", action = clap::ArgAction::Append)]
    merge_logic: Vec<String>,

    /// Add INFO/TAG flag to sites which are ("+") or are not ("-") listed in the -a file.
    #[arg(
        short = 'm',
        long = "mark-sites",
        value_name = "[+-]TAG",
        allow_hyphen_values = true
    )]
    mark_sites: Option<String>,

    /// Required overlap as a fraction of variant in the -a file (ANN), the VCF (:VCF), or both.
    #[arg(long = "min-overlap", value_name = "ANN:VCF")]
    min_overlap: Option<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "u|b|v|z[0-9]", action = clap::ArgAction::Append)]
    output_type: Vec<String>,

    /// Matching records by <snps|indels|both|all|some|exact|id> [some].
    #[arg(long = "pair-logic", visible_alias = "collapse", value_name = "STR", action = clap::ArgAction::Append)]
    pair_logic: Vec<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(short = 'r', long = "regions", value_name = "REGION")]
    regions: Option<String>,

    /// Restrict to regions listed in FILE.
    #[arg(short = 'R', long = "regions-file", value_name = "FILE")]
    regions_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(long = "regions-overlap", value_name = "0|1|2", action = clap::ArgAction::Append)]
    regions_overlap: Vec<String>,

    /// Rename annotations: TYPE/old\tnew, where TYPE is one of FILTER,INFO,FORMAT.
    #[arg(long = "rename-annots", value_name = "FILE")]
    rename_annots: Option<String>,

    /// Rename sequences according to the mapping: old\tnew.
    #[arg(long = "rename-chrs", value_name = "FILE")]
    rename_chrs: Option<String>,

    /// Comma separated list of samples to annotate (or exclude with "^" prefix).
    #[arg(
        short = 's',
        long = "samples",
        value_name = "[^]LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to annotate (or exclude with "^" prefix).
    #[arg(short = 'S', long = "samples-file", value_name = "[^]FILE")]
    samples_file: Option<String>,

    /// Keep memory low by avoiding complexities arising from handling multiple overlapping
    /// intervals.
    #[arg(long = "single-overlaps")]
    single_overlaps: bool,

    /// List of annotations (e.g. ID,INFO/DP,FORMAT/DP,FILTER) to remove (or keep with "^").
    #[arg(
        short = 'x',
        long = "remove",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    remove: Option<String>,

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
    #[arg(short = '?', long = "help", action = clap::ArgAction::Count)]
    help: u8,
}

/// `-m` logic (bcftools' `MARK_LISTED`, `MARK_UNLISTED`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MarkLogic {
    Listed,
    Unlisted,
}

/// Where annotations come from.
enum Source {
    /// Nothing to transfer (only `-x`, `-I`, `-m`, or renames).
    None,
    /// From the record itself (`-c` without `-a`).
    SelfColumns,
    /// A tab-delimited file.
    Tab(Box<TabAnnot>),
    /// A VCF read as the second reader.
    Vcf,
}

/// The input: one reader (with parallel parsing) or the input and the annotation VCF.
enum Input {
    Single(Box<SyncedReader>, Option<Box<ParallelVcf>>),
    Pair(Box<SyncedReaders>),
}

/// The annotation state (bcftools' `args_t`).
struct Annotator {
    ctx: Ctx,
    cols: Vec<Column>,
    source: Source,
    rm: Option<Removals>,
    renames: Vec<Rename>,
    filter: Option<(Filter, FilterLogic)>,
    filter_ext: Option<ExtFilter>,
    keep_sites: bool,
    set_ids: Option<(Convert, bool)>,
    mark_sites: Option<(Vec<u8>, MarkLogic)>,
    has_targets: bool,
    force: bool,
    id_buf: Vec<u8>,
    force_warned: bool,
    /// Whether `--force` went past a reader error (bcftools' `vcf_parse_error_warned`).
    errnum_warned: bool,
}

/// C's `strtod(s, &end)` over a whole string.
fn strtod_prefix(s: &str) -> (f64, usize) {
    pasta::variant::number::strtod(s.as_bytes(), 0)
}

/// The last value of a repeated option across two option ids, by command-line position.
fn last_vec<'a>(
    matches: &ArgMatches,
    a: (&str, &'a [String]),
    b: (&str, &'a [String]),
) -> Option<&'a str> {
    let pos = |id: &str| last_index(matches, id);
    match (a.1.last(), b.1.last()) {
        (Some(x), Some(y)) => Some(if pos(a.0) > pos(b.0) { x } else { y }),
        (Some(x), None) => Some(x),
        (None, Some(y)) => Some(y),
        (None, None) => None,
    }
}

/// Parses `--pair-logic` values in order (bcftools' option 2).
fn parse_pair_logic(values: &[String]) -> anyhow::Result<Option<Pairing>> {
    let mut p: Option<Pairing> = None;
    for v in values {
        let cur = p.unwrap_or_default();
        p = Some(match v.as_str() {
            "snps" => cur | Pairing::SNP_REF,
            "indels" => cur | Pairing::INDEL_REF,
            "both" => cur | Pairing::BOTH_REF,
            "any" | "all" => cur | Pairing::ANY,
            "some" => cur | Pairing::SOME,
            "none" | "exact" => Pairing::EXACT,
            "id" => cur | Pairing::ID,
            _ => bail!("The --pair-logic string \"{v}\" not recognised."),
        });
    }
    Ok(p)
}

/// Whether a file is VCF or BCF, as `hts_get_format` decides (other text is tab-delimited);
/// `hts_open` refuses binary data it does not know. `FILE##idx##INDEX` opens `FILE`.
fn is_vcf_or_bcf(path: &str) -> anyhow::Result<bool> {
    let (data, _) = split_idx_name(path);
    let input = match RawInput::open(data) {
        Ok(input) => input,
        Err(e) => {
            output::log_open_failure(data, &e);
            bail!("Failed to open {path}");
        }
    };
    if input.is_unknown_binary() {
        report_exec_format_error(data);
        bail!("Failed to open {path}");
    }
    Ok(input.format != Format::Other)
}

/// Parses the `{}` external values of a `-i/-e` expression (bcftools' `init_filters`):
/// returns the expression with each `{[type:]NAME}` replaced by `{}` (or `{str}`, `{int}`,
/// `{float}`) and the names. Braces inside quotes are left alone.
fn parse_ext_values(expr: &str) -> anyhow::Result<(String, Vec<Vec<u8>>)> {
    let s = expr.as_bytes();
    let mut out = Vec::with_capacity(s.len());
    let mut names = Vec::new();
    let mut i = 0;
    let mut quote: Option<u8> = None;
    while i < s.len() {
        let c = s[i];
        match quote {
            Some(q) => {
                if c == q {
                    quote = None;
                }
            }
            None if c == b'"' || c == b'\'' => quote = Some(c),
            None if c == b'{' => {
                let Some(close) = s[i + 1..].iter().position(|&c| c == b'}') else {
                    bail!("Could not parse the expression: {expr}");
                };
                let inner = &s[i + 1..i + 1 + close];
                let (name, token): (&[u8], &[u8]) = match inner.iter().position(|&c| c == b':') {
                    Some(colon) => {
                        let ty = &inner[..colon];
                        let token: Option<&[u8]> = if ty.eq_ignore_ascii_case(b"str") {
                            Some(b"{str}")
                        } else if ty.eq_ignore_ascii_case(b"int") {
                            Some(b"{int}")
                        } else if ty.eq_ignore_ascii_case(b"float") {
                            Some(b"{float}")
                        } else {
                            None
                        };
                        // An unrecognised type names the value after the text before the colon,
                        // as in bcftools.
                        match token {
                            Some(t) => (&inner[colon + 1..], t),
                            None => (ty, b"{}"),
                        }
                    }
                    None => (inner, b"{}"),
                };
                out.extend_from_slice(token);
                names.push(name.to_vec());
                i += close + 2;
                continue;
            }
            None => {}
        }
        out.push(c);
        i += 1;
    }
    if quote.is_some() && !names.is_empty() {
        bail!("Could not parse the expression: {expr}");
    }
    Ok((String::from_utf8_lossy(&out).into_owned(), names))
}

/// Reads `-h FILE` lines, up to the first empty line (`hts_getline` returning 0).
/// `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
fn read_header_file(path: &str) -> anyhow::Result<Vec<Vec<u8>>> {
    let data = split_idx_name(path).0;
    let lines = match readlist::read_all_lines(Path::new(data)) {
        Ok(lines) => lines,
        Err(e) => {
            output::log_open_failure(data, &e);
            bail!("Error reading {path}");
        }
    };
    Ok(lines.into_iter().take_while(|l| !l.is_empty()).collect())
}

/// Parses `--min-overlap ANN:VCF`.
fn parse_min_overlap(s: &str) -> anyhow::Result<(f32, f32)> {
    let fail =
        || anyhow::anyhow!("Could not parse \"--min-overlap {s}\", expected value(s) between 0-1");
    let (mut ann, mut vcf) = (0.0f32, 0.0f32);
    let mut rest = s;
    if !s.starts_with(':') {
        let (v, end) = strtod_prefix(s);
        ann = v as f32;
        rest = &s[end..];
        if !(0.0..=1.0).contains(&ann) || (!rest.is_empty() && !rest.starts_with(':')) {
            return Err(fail());
        }
    }
    if let Some(r) = rest.strip_prefix(':') {
        let (v, end) = strtod_prefix(r);
        vcf = v as f32;
        if !(0.0..=1.0).contains(&vcf) || end != r.len() {
            return Err(fail());
        }
    }
    Ok((ann, vcf))
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    // Every -v is checked, in order, as getopt hands each one to bcftools.
    for v in &args.verbosity {
        apply_verbosity(v)?;
    }
    if args.help > 0 {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    let filter_expr =
        FilterLogic::single(&args.include, &args.exclude)?.map(|(e, logic)| (e.to_owned(), logic));
    let mut pair_logic = parse_pair_logic(&args.pair_logic)?;
    let regions_overlap = args
        .regions_overlap
        .iter()
        .try_fold(Overlap::Record, |_, s| match Overlap::parse(s) {
            Some(o) => Ok(o),
            None => bail!("Could not parse: --regions-overlap {s}"),
        })?;
    // C's strtol(optarg, 0, 0): a negative count fails when the thread pool is created.
    let n_threads = args
        .threads
        .as_deref()
        .map_or(0, |t| i64::from(output::thread_count(t)));
    let threads = usize::try_from(n_threads).unwrap_or(0);
    let write_index = args
        .write_index
        .iter()
        .try_fold(None, |_, w| output::parse_write_index(w).map(Some))?;
    let settings: Vec<OutputSetting> = args
        .output_type
        .iter()
        .map(|v| OutputSetting::Type(v.clone()))
        .collect();
    let (out_type, level) = output::resolve(&settings, args.output.as_deref())?;

    let fname = match args.inputs.first() {
        Some(f) => f.clone(),
        None if io::stdin().is_terminal() => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
        None => "-".to_owned(),
    };

    let mut sr = SyncedReaders::new();
    if let Some((spec, source)) = later_file_sticky(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    ) {
        sr.set_regions_overlap(regions_overlap);
        if sr.set_regions_spec(spec.as_bytes(), source).is_err() {
            bail!("Failed to read the regions: {spec}");
        }
    }
    let kind = match &args.annotations {
        Some(a) if is_vcf_or_bcf(a)? => AnnotKind::Vcf,
        Some(_) => AnnotKind::Tab,
        None => AnnotKind::NoFile,
    };
    if kind == AnnotKind::Vcf {
        sr.set_index_requirement(IndexRequirement::Required);
        if args.min_overlap.is_some() {
            bail!("The --min-overlap option cannot be used when annotating from a VCF");
        }
    }
    if args.min_overlap.is_some() && args.single_overlaps {
        bail!("The options --single-overlaps and --min-overlap cannot be combined");
    }
    if n_threads < 0 || !output::threads_start(threads) {
        bail!("Failed to create threads");
    }
    // A VCF annotation file is read in step with the input.
    sr.set_shared_threads(threads, if kind == AnnotKind::Vcf { 2 } else { 1 });
    if let Err(e) = sr.add_reader(&fname) {
        let name = if fname == "-" {
            "standard input"
        } else {
            fname.as_str()
        };
        bail!("Failed to read from {name}: {e}");
    }

    let output_name = args
        .output
        .as_deref()
        .map_or_else(|| "-".to_owned(), |p| p.display().to_string());
    let (mut annotator, mut input, hdr_out) =
        init_data(&args, matches, sr, kind, &mut pair_logic, filter_expr)?;

    let options = output::write_options(out_type, level, threads);
    let mut writer = open_writer(&args, &hdr_out, &options, write_index, &output_name)?;
    annotator.ctx.hdr_out = hdr_out;
    let mut emitter = Emitter::new(
        &annotator.ctx.hdr_out,
        &writer,
        threads,
        output_name.clone(),
    )
    .with_write_error(format!(
        "[main_vcfannotate] Error: failed to write to {output_name}"
    ));
    if let Input::Single(r, par) = &mut input {
        *par = ParallelVcf::new(r, threads).map(Box::new);
    }
    let result = annotator.main_loop(&mut input, &mut writer, &mut emitter);
    if let Err(e) = result {
        // With threads, write first what a serial run would have written.
        emitter.abort(&annotator.ctx.hdr_out, &mut writer);
        writer.abandon();
        return Err(e);
    }
    close_writer(writer, write_index, &output_name)?;
    Ok(ExitCode::SUCCESS)
}

/// Opens the output, sets up `-W`, and writes the header (the output part of bcftools'
/// `init_data`).
fn open_writer(
    args: &Args,
    hdr_out: &Header,
    options: &WriteOptions,
    write_index: Option<output::WriteIndex>,
    output_name: &str,
) -> anyhow::Result<Writer> {
    let out_path = args.output.as_deref().map(output::data_path);
    let mut writer = match Writer::create(out_path.as_deref(), options) {
        Ok(w) => w,
        Err(e) => bail!(
            "[init_data] Error: cannot write to \"{output_name}\": {}",
            crate::raw_input::strerror(&e)
        ),
    };
    if let Some(fmt) = write_index {
        let ok = output::index_request(fmt, options.output_type, args.output.as_deref())
            .is_some_and(|req| writer.request_index(req).is_ok());
        if !ok {
            if writer.write_header(hdr_out).is_err() {
                bail!("[init_data] Error: failed to write the header to {output_name}");
            }
            writer.abandon();
            bail!("Error: failed to initialise index for {output_name}");
        }
    }
    if writer.write_header(hdr_out).is_err() {
        bail!("[init_data] Error: failed to write the header to {output_name}");
    }
    Ok(writer)
}

/// Closes the output and its index.
fn close_writer(
    writer: Writer,
    write_index: Option<output::WriteIndex>,
    output_name: &str,
) -> anyhow::Result<()> {
    match writer.finish() {
        Ok(()) => Ok(()),
        Err(WriteError::Index(_)) if write_index.is_some() => {
            bail!("Error: cannot write to index {output_name}");
        }
        Err(WriteError::Deferred(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("[main_vcfannotate] Error: failed to write to {output_name}")
        })),
        Err(WriteError::Io(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("Error: close failed .. {output_name}")
        })),
        Err(e) => {
            Err(anyhow::Error::new(e).context(format!("Error: close failed .. {output_name}")))
        }
    }
}

/// bcftools' `init_data`: everything but opening the output. Returns the annotator, the
/// input, and the output header.
fn init_data(
    args: &Args,
    matches: &ArgMatches,
    mut sr: SyncedReaders,
    kind: AnnotKind,
    pair_logic: &mut Option<Pairing>,
    filter_expr: Option<(String, FilterLogic)>,
) -> anyhow::Result<(Annotator, Input, Header)> {
    let tgts_is_vcf = kind == AnnotKind::Vcf;
    let mut hdr_out = sr
        .header(0)
        .dup()
        .context("Failed to duplicate the header")?;

    let mut set_ids = None;
    let mut set_ids_fmt = String::new();
    if let Some(fmt) = &args.set_id {
        let (replace, fmt) = match fmt.strip_prefix('+') {
            Some(rest) => (false, rest),
            None => (true, fmt.as_str()),
        };
        set_ids_fmt = fmt.to_owned();
        set_ids = Some((Convert::new(&hdr_out, None, fmt)?, replace));
    }
    let mut rm = None;
    if let Some(list) = &args.remove {
        let spec = RemoveSpec {
            list,
            keep_sites: args.keep_sites,
            set_id_reads_info: set_ids.as_ref().is_some_and(|(c, _)| c.uses_info_tags()),
            set_id_fmt: &set_ids_fmt,
        };
        rm = Some(remove::init_remove_annots(
            &spec,
            sr.header(0),
            &mut hdr_out,
        )?);
    }
    append_header_lines(args, &mut sr, &mut hdr_out)?;
    let targets_fname = args.annotations.as_deref();
    // Detecting the annotation file's type read it once already.
    if let (Some(a), true) = (targets_fname, tgts_is_vcf)
        && let Err(e) = sr.add_reader_again(a)
    {
        bail!("Failed to open {a}: {e}");
    }

    let mut merge_logic: Vec<u8> = args.merge_logic.join(",").into_bytes();
    let parsed = parse_columns(
        args,
        matches,
        &sr,
        kind,
        &mut hdr_out,
        pair_logic,
        &mut merge_logic,
    )?;
    let columns_given = parsed.is_some();
    let mut parsed = parsed.unwrap_or_else(|| columns::Columns {
        match_id: -1,
        match_end: -1,
        ..columns::Columns::default()
    });

    let mut source = if tgts_is_vcf {
        Source::Vcf
    } else if parsed.cols.is_empty() {
        Source::None
    } else {
        Source::SelfColumns
    };
    if let (Some(a), false) = (targets_fname, tgts_is_vcf) {
        if !columns_given {
            bail!("The -c option not given");
        }
        let tab = init_tab_source(args, a, &mut parsed, &merge_logic)?;
        source = Source::Tab(Box::new(tab));
    }
    init_merge_method(&mut parsed, &merge_logic, &mut source)?;
    let (filter, filter_ext) =
        init_filters(filter_expr.as_ref(), sr.header(0), &source, &parsed.cols)?;
    let mark_sites = match &args.mark_sites {
        Some(m) => Some(init_mark_sites(m, targets_fname, &mut hdr_out)?),
        None => None,
    };

    if !args.no_version {
        version::append_version(&mut hdr_out, "annotate")?;
    }
    if let Some(f) = &args.rename_chrs {
        rename::rename_chrs(f, sr.header(0), &mut hdr_out)?;
    }
    let mut renames = Vec::new();
    if args.rename_annots.is_some() || !parsed.renames.is_empty() {
        if let Some(f) = &args.rename_annots {
            rename::read_rename_file(f, &mut parsed.renames)?;
        }
        renames = rename::init_renames(&parsed.renames, &mut hdr_out)?;
    }
    if tgts_is_vcf {
        sr.set_pairing(pair_logic.unwrap_or(Pairing::SOME));
    }

    let ctx = Ctx {
        hdr_out: Header::new(),
        vcmp: Vcmp::default(),
        allele_cols: AlleleColumns::new(&parsed.cols),
        copied_alleles: CopiedAlleles::None,
        smpl: std::mem::take(&mut parsed.samples),
        targets_fname: targets_fname.unwrap_or("(null)").to_owned(),
        s: Scratch::default(),
    };
    let annotator = Annotator {
        ctx,
        cols: parsed.cols,
        source,
        rm,
        renames,
        filter,
        filter_ext,
        keep_sites: args.keep_sites,
        set_ids,
        mark_sites,
        has_targets: targets_fname.is_some(),
        force: args.force,
        id_buf: Vec::new(),
        force_warned: false,
        errnum_warned: false,
    };
    let input = if tgts_is_vcf {
        Input::Pair(Box::new(sr))
    } else {
        Input::Single(Box::new(SyncedReader::from(sr)), None)
    };
    Ok((annotator, input, hdr_out))
}

/// Appends `-h FILE` and `-H STR` lines to the output header and to the input's.
fn append_header_lines(
    args: &Args,
    sr: &mut SyncedReaders,
    hdr_out: &mut Header,
) -> anyhow::Result<()> {
    if let Some(path) = &args.header_lines {
        for l in read_header_file(path)? {
            if hdr_out.append(&l).is_err() {
                bail!("Could not parse {path}: {}", lossy(&l));
            }
            let _ = sr.header_mut(0).append(&l);
        }
    }
    for l in &args.header_line {
        if hdr_out.append(l.as_bytes()).is_err() {
            bail!("Could not parse the header line: {l}");
        }
        let _ = sr.header_mut(0).append(l.as_bytes());
    }
    Ok(())
}

/// Parses the last `-c` or `-C` with `init_columns`; `None` when neither is given.
fn parse_columns(
    args: &Args,
    matches: &ArgMatches,
    sr: &SyncedReaders,
    kind: AnnotKind,
    hdr_out: &mut Header,
    pair_logic: &mut Option<Pairing>,
    merge_logic: &mut Vec<u8>,
) -> anyhow::Result<Option<columns::Columns>> {
    let columns_arg = last_vec(
        matches,
        ("columns", &args.columns),
        ("columns_file", &args.columns_file),
    );
    let samples_arg = later_file_sticky(
        matches,
        ("samples", args.samples.as_deref()),
        ("samples_file", args.samples_file.as_deref()),
    );
    let Some(cols) = columns_arg else {
        return Ok(None);
    };
    let spec = ColumnsSpec {
        columns: cols,
        columns_is_file: !args.columns_file.is_empty(),
        targets_fname: args.annotations.as_deref(),
        tgts_hdr: (kind == AnnotKind::Vcf).then(|| sr.header(1)),
        samples: samples_arg,
    };
    let parsed = columns::init_columns(&spec, sr.header(0), hdr_out, pair_logic, merge_logic)?;
    Ok(Some(parsed))
}

/// Opens the tab-delimited annotation file `a` (the part of bcftools' `init_data` for
/// `-a` without a VCF): streamed through its tabix index, or read into memory when the
/// columns give the end of each region.
fn init_tab_source(
    args: &Args,
    a: &str,
    parsed: &mut columns::Columns,
    merge_logic: &[u8],
) -> anyhow::Result<TabAnnot> {
    let tc = &mut parsed.tab;
    if tc.chr == -1 {
        bail!("The -c CHROM option not given");
    }
    if tc.beg == -1 {
        bail!("The -c POS option not given");
    }
    if args.single_overlaps && !merge_logic.is_empty() {
        bail!("The options --merge-logic and --single-overlaps cannot be combined");
    }
    // Detecting the file's type read it once already. A stream has no index, and what it
    // holds then is gone (bcftools reads the rest of a pipe, which a small input leaves
    // empty, and waits forever on a FIFO); it is not opened again.
    let stream = is_stream(a);
    let tab_source = if tc.end == -1 || (args.single_overlaps && merge_logic.is_empty()) {
        tc.end = -tc.beg - 1;
        if stream {
            bail!("Could not initialize the annotation file: {a}");
        }
        // Without an index, bcftools reads the whole file to check it, then stops.
        let regions = Regions::from_file_columns(
            Path::new(a),
            tc.chr as usize,
            tc.beg as usize,
            tc.end as isize,
        )
        .map_err(|_| anyhow::anyhow!("Could not initialize the annotation file: {a}"))?;
        if !regions.is_file_backed() {
            bail!("Expected tabix-indexed annotation file: {a}");
        }
        TabSource::Tabix(Box::new(regions))
    } else {
        if tc.ref_ != -1 {
            bail!(
                "Error: the REF columns will be ignored when BEG,END (or FROM,TO) is present. Replace END (or TO) with \"-\"."
            );
        }
        // The data file's name, where bcftools tests the whole `##idx##` name and so reads
        // a BED file with a named index as 1-based.
        let lower = split_idx_name(a).0.to_ascii_lowercase();
        let is_bed =
            lower.ends_with(".bed.gz") || lower.ends_with(".bed.bgz") || lower.ends_with(".bed");
        let parse = |l: &[u8]| -> regidx::ParseResult<Vec<u8>> {
            let p = if is_bed {
                regidx::parse_bed(l)?
            } else {
                regidx::parse_tab(l)?
            };
            Ok(p.map(|p| p.with_payload(l.to_vec())))
        };
        let idx = if stream {
            RegIdx::new()
        } else {
            RegIdx::from_file_with(Path::new(a), parse).map_err(|e| {
                eprintln!("{e}");
                anyhow::anyhow!("Failed to parse: {a}")
            })?
        };
        TabSource::Idx(idx)
    };
    let mut tab = TabAnnot::new(tab_source, parsed.tab);
    tab.match_id = parsed.match_id;
    tab.match_end = parsed.match_end;
    if let Some(mo) = &args.min_overlap {
        let (ann, vcf) = parse_min_overlap(mo)?;
        tab.min_overlap_ann = ann;
        tab.min_overlap_vcf = vcf;
    }
    Ok(tab)
}

/// bcftools' `init_merge_method`, and the line of missing values that an appending tab
/// source fills unmatched columns with.
fn init_merge_method(
    parsed: &mut columns::Columns,
    merge_logic: &[u8],
    source: &mut Source,
) -> anyhow::Result<()> {
    let ref_idx = parsed.tab.ref_;
    let kind = match source {
        Source::Vcf => AnnotKind::Vcf,
        Source::Tab(_) => AnnotKind::Tab,
        Source::None | Source::SelfColumns => AnnotKind::NoFile,
    };
    let has_append_mode = columns::init_merge_method(&mut parsed.cols, merge_logic, kind, ref_idx)?;
    if let Source::Tab(tab) = source {
        tab.has_append_mode = has_append_mode;
        if has_append_mode {
            let ncol = parsed
                .cols
                .iter()
                .map(|c| c.icol + 1)
                .chain([ref_idx + 1])
                .max()
                .unwrap_or(0)
                .max(0) as usize;
            tab.aline_missing = AnnotLine::missing(ncol);
        }
    }
    Ok(())
}

/// The `-i/-e` filter, plain or with values from annotation columns.
type Filters = (Option<(Filter, FilterLogic)>, Option<ExtFilter>);

/// bcftools' `init_filters`: compiles `-i/-e`, whose `{NAME}` values come from the columns
/// of a tab-delimited annotation file.
fn init_filters(
    filter_expr: Option<&(String, FilterLogic)>,
    hdr: &Header,
    source: &Source,
    cols: &[Column],
) -> anyhow::Result<Filters> {
    let Some((expr, logic)) = filter_expr else {
        return Ok((None, None));
    };
    let (expr2, ext) = parse_ext_values(expr)?;
    let f = Filter::init(hdr, &expr2)?;
    let types = f.ext_types().to_vec();
    if types.len() != ext.len() {
        bail!(
            "Failed to parse the expression, unexpected number of dynamic variables ({} vs {}): {expr}",
            types.len(),
            ext.len()
        );
    }
    if ext.is_empty() {
        return Ok((Some((f, *logic)), None));
    }
    if !matches!(source, Source::Tab(_)) {
        bail!(
            "Error: dynamic variables in -i/-e expressions can be currently used only with tab-delimited file, not with VCF (todo)"
        );
    }
    let mut exts = Vec::with_capacity(ext.len());
    for (name, ty) in ext.into_iter().zip(types) {
        // bcftools leaves the column unset (uninitialised) when -c names no annotation
        // columns; that is an error here.
        let Some(c) = cols.iter().find(|c| c.key_src == name) else {
            bail!("No such column: {}", lossy(&name));
        };
        exts.push(Ext { ty, icol: c.icol });
    }
    let filter_ext = ExtFilter {
        filter: f,
        ext: exts,
        logic: *logic,
    };
    Ok((None, Some(filter_ext)))
}

/// Parses `-m [+-]TAG` and adds the definition of the flag to `hdr_out`.
fn init_mark_sites(
    m: &str,
    targets_fname: Option<&str>,
    hdr_out: &mut Header,
) -> anyhow::Result<(Vec<u8>, MarkLogic)> {
    let (tag, logic) = match m.as_bytes().first() {
        Some(b'+') => (&m[1..], MarkLogic::Listed),
        Some(b'-') => (&m[1..], MarkLogic::Unlisted),
        _ => (m, MarkLogic::Listed),
    };
    let line = match targets_fname {
        None => {
            if logic != MarkLogic::Listed {
                bail!("The -a option not given but -{tag} logic was requested");
            }
            eprintln!("Note: The -a option not given, all sites will be annotated with INFO/{tag}");
            format!(
                "##INFO=<ID={tag},Number=0,Type=Flag,Description=\"Sites marked with `bcftools annotate -m {tag}`\">"
            )
        }
        Some(_) => format!(
            "##INFO=<ID={tag},Number=0,Type=Flag,Description=\"Sites {}listed in {tag}\">",
            if logic == MarkLogic::Listed {
                ""
            } else {
                "not "
            }
        ),
    };
    let _ = hdr_out.append(line.as_bytes());
    Ok((tag.as_bytes().to_vec(), logic))
}

impl Annotator {
    /// The main loop of `main_vcfannotate`.
    fn main_loop(
        &mut self,
        input: &mut Input,
        writer: &mut Writer,
        emitter: &mut Emitter,
    ) -> anyhow::Result<()> {
        let mut rec = Record::new();
        let mut arec = Record::new();
        let mut spare: Vec<Record> = Vec::new();
        let mut read_error = None;
        loop {
            let has_annot = match input {
                Input::Single(r, par) => {
                    let got = match par {
                        Some(p) => p.next(r, &mut rec, &mut spare),
                        None => r.next_record(&mut rec),
                    };
                    match got {
                        Ok(true) => {
                            self.check_errnum(r.pending_error())?;
                            false
                        }
                        Ok(false) => break,
                        Err(e) => {
                            read_error = Some(e);
                            break;
                        }
                    }
                }
                Input::Pair(sr) => {
                    match sr.next_line() {
                        Ok(0) => break,
                        Ok(_) => self.check_errnum(sr.pending_error())?,
                        Err(e) => {
                            read_error = Some(e);
                            break;
                        }
                    }
                    if !sr.has_line(0) {
                        continue;
                    }
                    sr.swap_line(0, &mut rec);
                    let has = sr.has_line(1);
                    if has {
                        sr.swap_line(1, &mut arec);
                    }
                    has
                }
            };
            if !rec.errors().is_empty() {
                if !self.force {
                    bail!(
                        "Encountered an error, cannot proceed. Please check the error output above.\nIf feeling adventurous, use the --force option. (At your own risk!)"
                    );
                }
                if !self.force_warned {
                    eprintln!(
                        "Warning: Encountered an error, proceeding only because --force was given.\n         Note that this can result in a segfault or a silent corruption of the output file!"
                    );
                    self.force_warned = true;
                }
                rec.clear_errors();
            }
            let (hdr, ahdr) = match input {
                Input::Single(r, _) => (r.header(), None),
                Input::Pair(sr) => (sr.header(0), Some(sr.header(1))),
            };
            let mut write = true;
            if let Some((f, logic)) = &mut self.filter {
                let mut pass = f.test(hdr, &mut rec)?;
                if *logic == FilterLogic::Exclude {
                    pass = !pass;
                }
                if !pass {
                    write = self.keep_sites;
                    if write {
                        rec = self.emit(emitter, writer, rec)?;
                    }
                    continue;
                }
            }
            let annot = if has_annot {
                ahdr.map(|h| (&mut arec, h))
            } else {
                None
            };
            let keep = self.annotate_line(hdr, &mut rec, annot)?;
            if self.filter_ext.is_some() && !self.keep_sites && !keep {
                write = false;
            }
            if write {
                rec = self.emit(emitter, writer, rec)?;
            }
        }
        emitter.finish(&self.ctx.hdr_out, writer)?;
        let pending = match input {
            Input::Single(r, _) => r.take_pending_error(),
            Input::Pair(sr) => sr.take_pending_error(),
        };
        if let Some(e) = read_error.or(pending) {
            bail!("Error: {e}");
        }
        Ok(())
    }

    /// The check of the synced reader's `errnum` after each line: a missing BGZF EOF marker
    /// of the last input opened stops at the first record unless `--force` is given.
    fn check_errnum(&mut self, pending: Option<&SyncedError>) -> anyhow::Result<()> {
        let Some(e) = pending else {
            return Ok(());
        };
        if !self.force {
            bail!("Error: {e}");
        }
        if !self.errnum_warned {
            eprintln!("Warning: Encountered an error, proceeding only because --force was given.");
            self.errnum_warned = true;
        }
        Ok(())
    }

    /// Writes a record and returns one to read the next into.
    fn emit(
        &mut self,
        emitter: &mut Emitter,
        writer: &mut Writer,
        rec: Record,
    ) -> anyhow::Result<Record> {
        let back = emitter.write(&self.ctx.hdr_out, writer, rec)?;
        Ok(back.or_else(|| emitter.spare()).unwrap_or_default())
    }

    /// bcftools' `annotate_line`; returns whether an annotation overlapped.
    fn annotate_line(
        &mut self,
        hdr: &Header,
        line: &mut Record,
        annot: Option<(&mut Record, &Header)>,
    ) -> anyhow::Result<bool> {
        if let Some(rm) = &self.rm {
            rm.apply(hdr, line)?;
        }
        if !self.renames.is_empty() {
            rename::apply(&self.renames, line);
        }
        let mut has_overlap = false;
        match &mut self.source {
            Source::Tab(tab) => {
                let mut a = Annotate {
                    ctx: &mut self.ctx,
                    cols: &mut self.cols,
                    hdr,
                    ext: self.filter_ext.as_mut(),
                    keep_sites: self.keep_sites,
                    mark_sites: self.mark_sites.is_some(),
                };
                has_overlap = tab.annotate(&mut a, line)?;
            }
            Source::Vcf => {
                if let Some((arec, ahdr)) = annot {
                    self.ctx.copied_alleles = self.ctx.allele_cols.copied(line);
                    for j in 0..self.cols.len() {
                        if self.cols[j].setter.is_none() {
                            continue;
                        }
                        let r = self.ctx.set(
                            &mut self.cols[j],
                            hdr,
                            line,
                            Src::Vcf(&mut *arec, ahdr),
                        )?;
                        if r == Set::Failed {
                            bail!(
                                "fixme: Could not set {} at {}",
                                lossy(&self.cols[j].key_src),
                                site(hdr, line)
                            );
                        }
                    }
                    has_overlap = true;
                }
            }
            Source::SelfColumns => {
                for j in 0..self.cols.len() {
                    if self.cols[j].setter.is_none() {
                        continue;
                    }
                    if self.ctx.set(&mut self.cols[j], hdr, line, Src::Flush)? == Set::Failed {
                        bail!(
                            "fixme: Could not set {} at {}",
                            lossy(&self.cols[j].key_src),
                            site(hdr, line)
                        );
                    }
                }
            }
            Source::None => {}
        }
        if let Some((conv, replace)) = &mut self.set_ids {
            conv.line(&self.ctx.hdr_out, line, None, &mut self.id_buf)?;
            if !self.id_buf.is_empty() {
                line.unpack(pasta::variant::Unpack::STR);
                if *replace || line.id().is_empty() || line.id() == b"." {
                    line.update_id(Some(&self.id_buf));
                }
            }
        }
        if let Some((tag, logic)) = &self.mark_sites {
            if !self.has_targets {
                has_overlap = true;
            }
            let set = (*logic == MarkLogic::Listed) == has_overlap;
            let _ = line.update_info_flag(&self.ctx.hdr_out, tag, set);
        }
        Ok(has_overlap)
    }
}

#[cfg(test)]
mod tests {
    use pasta::variant::{
        index::{
            IndexFormat,
            tabix::{PRESET_GENERIC, TabixBuilder, TabixConfig},
        },
        io::reader::open_again,
    };

    use super::*;
    use crate::raw_input::{make_fifo, within_seconds};

    #[test]
    fn ext_values_are_rewritten() {
        let (e, ext) = parse_ext_values("INFO/DP>{SCORE} && ID=\"{x}\"").expect("parses");
        assert_eq!(e, "INFO/DP>{} && ID=\"{x}\"");
        assert_eq!(ext.len(), 1);
        assert_eq!(ext[0], b"SCORE");
        let (e, ext) = parse_ext_values("INFO/DP>{int:SCORE}").expect("parses");
        assert_eq!(e, "INFO/DP>{int}");
        assert_eq!(ext[0], b"SCORE");
    }

    #[test]
    fn min_overlap_parses() {
        assert_eq!(parse_min_overlap("0.5").expect("ok"), (0.5, 0.0));
        assert_eq!(parse_min_overlap(":0.25").expect("ok"), (0.0, 0.25));
        assert_eq!(parse_min_overlap("1:1").expect("ok"), (1.0, 1.0));
        assert!(parse_min_overlap("2").is_err());
    }

    /// Runs `annotate` with the arguments after the command name.
    fn run_annotate(argv: &[String]) -> anyhow::Result<()> {
        use clap::{Args as _, FromArgMatches as _};

        let argv = std::iter::once("annotate".to_owned()).chain(argv.iter().cloned());
        let matches =
            Args::augment_args(clap::Command::new("annotate")).try_get_matches_from(argv)?;
        let args = Args::from_arg_matches(&matches)?;
        run(args, &matches).map(drop)
    }

    #[test]
    fn a_tab_delimited_fifo_is_not_opened_again() {
        let (dir, fifo) = make_fifo("annotate-tab");
        let vcf = dir.join("in.vcf");
        let hdr = dir.join("hdr.txt");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t10\t.\tA\tC\t.\t.\t.\n",
        )
        .expect("write the VCF");
        std::fs::write(
            &hdr,
            "##INFO=<ID=TAG,Number=1,Type=Integer,Description=\"Tag\">\n",
        )
        .expect("write the header lines");
        let path = |p: &std::path::Path| p.to_string_lossy().into_owned();
        let mut results = Vec::new();
        // Without END the file needs an index (which a stream has not); with it, it is read
        // into memory. Either way, detecting the type read the FIFO once.
        for columns in ["CHROM,POS,TAG", "CHROM,FROM,TO,TAG"] {
            let feeder = {
                let fifo = fifo.clone();
                std::thread::spawn(move || std::fs::write(&fifo, "chr1\t10\t12\t5\n"))
            };
            let argv: Vec<String> = vec![
                "-a".into(),
                path(&fifo),
                "-h".into(),
                path(&hdr),
                "-c".into(),
                columns.into(),
                "-o".into(),
                path(&dir.join("out.vcf")),
                path(&vcf),
            ];
            let result = within_seconds("annotate reading a FIFO", move || {
                run_annotate(&argv).map_err(|e| e.to_string())
            });
            // A reader for the writer, should the FIFO not have been opened.
            let partner = open_again(&fifo);
            let _ = feeder.join().expect("the writer thread does not panic");
            drop(partner);
            results.push(result);
        }
        let _ = std::fs::remove_dir_all(&dir);
        let fifo = path(&fifo);
        assert_eq!(
            results,
            [
                Err(format!("Could not initialize the annotation file: {fifo}")),
                Ok(())
            ]
        );
    }

    /// Writes `text` BGZF-compressed as `path`, with a tabix index of `conf` as `index`.
    fn write_tabix(path: &Path, text: &str, conf: TabixConfig, index: &Path, fmt: IndexFormat) {
        use std::io::Write as _;

        let mut w = pasta::bgzf::io::Writer::new(
            std::fs::File::create(path).expect("create the annotation file"),
        );
        let mut builder = TabixBuilder::new(conf, if fmt == IndexFormat::Csi { 14 } else { 0 });
        for line in text.lines() {
            writeln!(w, "{line}").expect("write a line");
            builder
                .push_line(line.as_bytes(), u64::from(w.virtual_position()))
                .expect("index a line");
        }
        let tbx = builder.finish(u64::from(w.virtual_position()));
        w.finish().expect("finish the BGZF file");
        tbx.index.save(index, fmt).expect("save the index");
    }

    #[test]
    fn annotation_files_take_an_index_named_with_idx() {
        use clap::{Args as _, FromArgMatches as _};

        use crate::commands::view;

        let dir = std::env::temp_dir().join(format!("vtools-annotate-idx-{}", std::process::id()));
        let idx = dir.join("idx");
        std::fs::create_dir_all(&idx).expect("create the test directories");
        let path = |p: &Path| p.to_string_lossy().into_owned();
        let named = |data: &str, index: &str| {
            format!("{}##idx##{}", path(&dir.join(data)), path(&idx.join(index)))
        };
        let dst = dir.join("dst.vcf");
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        );
        for pos in [4, 5, 6, 10, 15, 16, 20] {
            text.push_str(&format!("1\t{pos}\t.\tA\tC\t.\t.\t.\n"));
        }
        std::fs::write(&dst, text).expect("write the VCF");
        let hdr = dir.join("hdr.txt");
        std::fs::write(
            &hdr,
            "##INFO=<ID=YY,Number=1,Type=Integer,Description=\"y\">\n",
        )
        .expect("write the header lines");
        let ann = dir.join("ann.vcf");
        std::fs::write(
            &ann,
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
             ##INFO=<ID=YY,Number=1,Type=Integer,Description=\"y\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t10\t.\tA\tC\t.\t.\tYY=42\n1\t20\t.\tA\tC\t.\t.\tYY=43\n",
        )
        .expect("write the annotation VCF");
        // Annotating from a VCF needs an indexed input too.
        let dst_gz = path(&dir.join("dst.vcf.gz"));
        for (src, out, ot, fmt) in [
            (&dst, dst_gz.clone(), "-Oz", "--write-index=tbi"),
            (
                &ann,
                named("ann.vcf.gz", "a.tbi"),
                "-Oz",
                "--write-index=tbi",
            ),
            (&ann, named("ann.bcf", "b.csi"), "-Ob", "--write-index=csi"),
        ] {
            let argv = ["view", "--no-version", ot, fmt, "-o", &out, &path(src)];
            let matches = view::Args::augment_args(clap::Command::new("view"))
                .try_get_matches_from(argv)
                .expect("the view options parse");
            let args = view::Args::from_arg_matches(&matches).expect("the view options parse");
            view::run(args, &matches).expect("view writes the annotation file");
        }
        let generic = |ec| TabixConfig {
            preset: PRESET_GENERIC,
            sc: 1,
            bc: 2,
            ec,
            meta_char: i32::from(b'#'),
            line_skip: 0,
        };
        let tab = |name: &str, text: &str, conf, index: &str, fmt| {
            write_tabix(&dir.join(name), text, conf, &idx.join(index), fmt);
        };
        tab(
            "pos.tab.gz",
            "1\t10\t42\n1\t20\t43\n",
            generic(2),
            "pos.tbi",
            IndexFormat::Tbi,
        );
        tab(
            "rng.tab.gz",
            "1\t5\t15\t42\n1\t20\t20\t43\n",
            generic(3),
            "rng.csi",
            IndexFormat::Csi,
        );
        tab(
            "reg.bed.gz",
            "1\t5\t15\t42\n1\t19\t20\t43\n",
            TabixConfig::BED,
            "reg.tbi",
            IndexFormat::Tbi,
        );

        // The YY value written at each position, or '.'.
        let annotate = |a: String, columns: &str, extra: &[&str]| -> Result<Vec<String>, String> {
            let out = dir.join("out.vcf");
            let mut argv: Vec<String> = vec!["--no-version".into(), "-a".into(), a, "-c".into()];
            argv.push(columns.into());
            argv.extend(extra.iter().map(|&s| s.to_owned()));
            argv.extend([
                "-h".into(),
                path(&hdr),
                "-o".into(),
                path(&out),
                dst_gz.clone(),
            ]);
            run_annotate(&argv).map_err(|e| e.to_string())?;
            let text = std::fs::read_to_string(&out).expect("read the output");
            Ok(text
                .lines()
                .filter(|l| !l.starts_with('#'))
                .map(|l| l.split('\t').nth(7).unwrap_or_default().to_owned())
                .collect())
        };
        let results = [
            annotate(named("ann.vcf.gz", "a.tbi"), "YY", &[]),
            annotate(named("ann.bcf", "b.csi"), "YY", &[]),
            annotate(named("pos.tab.gz", "pos.tbi"), "CHROM,POS,YY", &[]),
            annotate(named("rng.tab.gz", "rng.csi"), "CHROM,FROM,TO,YY", &[]),
            annotate(
                named("rng.tab.gz", "rng.csi"),
                "CHROM,BEG,END,YY",
                &["--single-overlaps"],
            ),
            annotate(named("reg.bed.gz", "reg.tbi"), "CHROM,FROM,TO,YY", &[]),
            annotate(
                named("reg.bed.gz", "reg.tbi"),
                "CHROM,FROM,TO,YY",
                &["--single-overlaps"],
            ),
            // Without its index, a file of positions cannot be read; one of ranges can.
            annotate(named("pos.tab.gz", "missing.tbi"), "CHROM,POS,YY", &[]),
            annotate(named("rng.tab.gz", "missing.tbi"), "CHROM,FROM,TO,YY", &[]),
        ];
        let _ = std::fs::remove_dir_all(&dir);
        let at = |v: [&str; 7]| Ok(v.map(str::to_owned).to_vec());
        let points = at([".", ".", ".", "YY=42", ".", ".", "YY=43"]);
        let ranges = at([".", "YY=42", "YY=42", "YY=42", "YY=42", ".", "YY=43"]);
        let bed = at([".", ".", "YY=42", "YY=42", "YY=42", ".", "YY=43"]);
        let missing = Err(format!(
            "Expected tabix-indexed annotation file: {}",
            named("pos.tab.gz", "missing.tbi")
        ));
        assert_eq!(
            results,
            [
                points.clone(),
                points.clone(),
                points,
                ranges.clone(),
                ranges.clone(),
                bed.clone(),
                bed,
                missing,
                ranges
            ]
        );
    }

    #[test]
    fn idx_names_read_the_data_file_of_header_lines() {
        let (dir, name) =
            crate::raw_input::write_idx_named("annotate-hdr-idx", "h.txt", "##x=1\n\n##y=2\n");
        let lines = read_header_file(&name);
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(lines.expect("the data file is read"), [b"##x=1".to_vec()]);
    }

    #[test]
    fn pair_logic_accumulates() {
        let p = parse_pair_logic(&["snps".to_owned(), "id".to_owned()]).expect("parses");
        assert_eq!(p, Some(Pairing::SNP_REF | Pairing::ID));
        let p = parse_pair_logic(&["snps".to_owned(), "exact".to_owned()]).expect("parses");
        assert_eq!(p, Some(Pairing::EXACT));
        assert!(parse_pair_logic(&["bogus".to_owned()]).is_err());
    }

    #[test]
    fn a_lone_missing_string_fills_every_allele() {
        use clap::{Args as _, FromArgMatches as _};

        use crate::commands::view;

        let dir = std::env::temp_dir().join(format!("vtools-annotate-ar-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let path = |p: &Path| p.to_string_lossy().into_owned();
        let defs = "##INFO=<ID=SA,Number=A,Type=String,Description=\"a\">\n\
                    ##INFO=<ID=SR,Number=R,Type=String,Description=\"r\">\n";
        let head = "##fileformat=VCFv4.2\n##contig=<ID=1>\n";
        let columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let dst = dir.join("dst.vcf");
        std::fs::write(
            &dst,
            format!(
                "{head}{defs}{columns}1\t10\t.\tA\tC,G\t.\t.\t.\n1\t20\t.\tA\tG,C,T\t.\t.\tSA=p,.,.\n"
            ),
        )
        .expect("write the VCF");
        let hdr = dir.join("hdr.txt");
        std::fs::write(&hdr, defs).expect("write the header lines");
        let ann = dir.join("ann.vcf");
        std::fs::write(
            &ann,
            format!(
                "{head}{defs}{columns}1\t10\t.\tA\tC,G\t.\t.\tSA=.;SR=.\n\
                 1\t20\t.\tA\tC,G\t.\t.\tSA=.;SR=.\n"
            ),
        )
        .expect("write the annotation VCF");
        let tab = dir.join("ann.tab.gz");
        let conf = TabixConfig {
            preset: PRESET_GENERIC,
            sc: 1,
            bc: 2,
            ec: 2,
            meta_char: i32::from(b'#'),
            line_skip: 0,
        };
        write_tabix(
            &tab,
            "1\t10\tA\tC,G\t.\t.\n1\t20\tA\tC,G\t.\t.\n",
            conf,
            &dir.join("ann.tab.gz.tbi"),
            IndexFormat::Tbi,
        );
        // Annotating from a VCF needs it indexed, and the input bgzipped.
        let (ann_bcf, dst_gz) = (path(&dir.join("ann.bcf")), path(&dir.join("dst.vcf.gz")));
        for (src, out, ot) in [(&ann, &ann_bcf, "-Ob"), (&dst, &dst_gz, "-Oz")] {
            let argv = ["view", "--no-version", ot, "-W", "-o", out, &path(src)];
            let matches = view::Args::augment_args(clap::Command::new("view"))
                .try_get_matches_from(argv)
                .expect("the view options parse");
            let args = view::Args::from_arg_matches(&matches).expect("the view options parse");
            view::run(args, &matches).expect("view writes the file");
        }
        // The INFO column of each output record.
        let annotate = |argv: &[&str]| -> Result<Vec<String>, String> {
            let out = dir.join("out.vcf");
            let mut argv: Vec<String> = argv.iter().map(|&s| s.to_owned()).collect();
            argv.extend([
                "--no-version".into(),
                "-o".into(),
                path(&out),
                dst_gz.clone(),
            ]);
            run_annotate(&argv).map_err(|e| e.to_string())?;
            let text = std::fs::read_to_string(&out).expect("read the output");
            Ok(text
                .lines()
                .filter(|l| !l.starts_with('#'))
                .map(|l| l.split('\t').nth(7).unwrap_or_default().to_owned())
                .collect())
        };
        let (tab, hdr) = (path(&tab), path(&hdr));
        let results = [
            annotate(&["-a", &ann_bcf, "-c", "INFO/SA,INFO/SR"]),
            annotate(&["-a", &tab, "-h", &hdr, "-c", "CHROM,POS,REF,ALT,.SA,.SR"]),
        ];
        let _ = std::fs::remove_dir_all(&dir);
        let want = Ok(vec![
            "SA=.,.;SR=.,.,.".to_owned(),
            "SA=p,.,.;SR=.,.,.,.".to_owned(),
        ]);
        assert_eq!(results, [want.clone(), want]);
    }
}
