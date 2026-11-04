//! `vtools filter`: apply fixed-threshold filters (bcftools filter).

mod gaps;
mod mask;

use std::{
    io::{self, IsTerminal},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::bail;
use clap::ArgMatches;
use pasta::variant::{
    GetError, Header, HeaderRecordType, Record, Unpack, VariantType,
    encoding::{INT32_MISSING, INT32_VECTOR_END},
    io::{WriteError, Writer},
    readlist::{self, ListSource},
    synced::{Overlap, Regions, SyncedError, SyncedOptions, SyncedReader},
};

use self::{
    gaps::{GapBuffer, GapFilters},
    mask::Mask,
};
use crate::options::{last_index, later};
use crate::{
    commands::head::apply_verbosity,
    filter::{Filter, FilterLogic},
    output::{self, OutputSetting},
    raw_input::strerror,
    version,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin), then regions to filter (needs an index).
    inputs: Vec<String>,

    /// Exclude sites for which the expression is true (see man page for details).
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Filter SNPs within INT base pairs of an indel (the default) or any combination of
    /// indel,mnp,bnd,other,overlap.
    #[arg(
        short = 'g',
        long = "SnpGap",
        value_name = "INT[:TYPE]",
        allow_hyphen_values = true
    )]
    snp_gap: Option<String>,

    /// Filter clusters of indels separated by INT or fewer base pairs allowing only one to pass.
    #[arg(
        short = 'G',
        long = "IndelGap",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    indel_gap: Option<String>,

    /// Include only sites for which the expression is true (see man page for details).
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

    /// Soft filter regions, "^" to negate.
    #[arg(long = "mask", value_name = "[^]REGION", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    mask: Vec<String>,

    /// Soft filter regions listed in a file, "^" to negate.
    #[arg(short = 'M', long = "mask-file", value_name = "[^]FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    mask_file: Vec<String>,

    /// Mask if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(
        long = "mask-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true
    )]
    mask_overlap: Option<String>,

    /// "+": do not replace but add to existing FILTER; "x": reset filters at sites which pass.
    #[arg(short = 'm', long = "mode", value_name = "[+x]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    mode: Vec<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

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

    /// Annotate FILTER column with STRING or unique filter name ("Filter%d") made up by the
    /// program ("+").
    #[arg(
        short = 's',
        long = "soft-filter",
        value_name = "STRING",
        allow_hyphen_values = true
    )]
    soft_filter: Option<String>,

    /// Set genotypes of failed samples to missing (.) or ref (0).
    #[arg(
        short = 'S',
        long = "set-GTs",
        value_name = ".|0",
        allow_hyphen_values = true
    )]
    set_gts: Option<String>,

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
    #[arg(short = 'h', short_alias = '?')]
    help: bool,

    /// Not an option of bcftools filter: reported as unrecognized, with the usage.
    #[arg(long = "help", hide = true)]
    long_help: bool,
}

const USAGE: &str = "
About:   Apply fixed-threshold filters.
Usage:   vtools filter [options] <in.vcf.gz>

Options:
    -e, --exclude EXPR             Exclude sites for which the expression is true (see man page for details)
    -g, --SnpGap INT[:TYPE]        Filter SNPs within <int> base pairs of an indel (the default) or any combination of indel,mnp,bnd,other,overlap
    -G, --IndelGap INT             Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass
    -i, --include EXPR             Include only sites for which the expression is true (see man page for details
        --mask [^]REGION           Soft filter regions, \"^\" to negate
    -M, --mask-file [^]FILE        Soft filter regions listed in a file, \"^\" to negate
        --mask-overlap 0|1|2       Mask if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -m, --mode [+x]                \"+\": do not replace but add to existing FILTER; \"x\": reset filters at sites which pass
        --no-version               Do not append version and command line to the header
    -o, --output FILE              Write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -r, --regions REGION           Restrict to comma-separated list of regions
    -R, --regions-file FILE        Restrict to regions listed in a file
        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --soft-filter STRING       Annotate FILTER column with <string> or unique filter name (\"Filter%d\") made up by the program (\"+\")
    -S, --set-GTs .|0              Set genotypes of failed samples to missing (.) or ref (0)
    -t, --targets REGION           Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT              Use multithreading with <int> worker threads [0]
    -v, --verbosity INT            Verbosity level
    -W, --write-index[=FMT]        Automatically index the output files [off]
";

/// What `-S/--set-GTs` sets the genotypes of failing samples to.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SetGts {
    /// `.`
    Missing,
    /// `0`
    Ref,
}

/// How `-m/--mode` treats the FILTER column.
#[derive(Clone, Copy, Debug, Default)]
struct AnnotMode {
    /// `+`: add the soft filter to the existing ones instead of replacing them.
    add: bool,
    /// `x`: reset FILTER to PASS at sites which pass.
    reset: bool,
}

/// The `--mask` setting.
struct MaskSetting {
    mask: Mask,
    negate: bool,
    overlap: Overlap,
}

/// The settings of the main loop.
struct Settings {
    filter: Option<(Filter, FilterLogic)>,
    mask: Option<MaskSetting>,
    soft_filter: bool,
    annot: AnnotMode,
    set_gts: Option<SetGts>,
    flt_pass: i32,
    flt_fail: i32,
}

fn output_settings(args: &Args) -> Vec<OutputSetting> {
    args.output_type
        .iter()
        .map(|t| OutputSetting::Type(t.clone()))
        .collect()
}

fn parse_overlap(s: Option<&str>, default: Overlap, what: &str) -> anyhow::Result<Overlap> {
    match s {
        None => Ok(default),
        Some(s) => {
            Overlap::parse(s).ok_or_else(|| anyhow::anyhow!("Could not parse: --{what} {s}"))
        }
    }
}

/// Parses `-g INT[:TYPE,...]` into the gap and the variant types near which SNPs fail.
fn parse_snp_gap(s: &str) -> anyhow::Result<(i64, u32, String)> {
    let (gap, end) = pasta::variant::number::strtol(s.as_bytes(), 0);
    let rest = &s[end..];
    if !rest.is_empty() && !rest.starts_with(':') {
        bail!("Could not parse argument: --SnpGap {s}");
    }
    let Some(list) = rest.strip_prefix(':') else {
        return Ok((gap, VariantType::INDEL.bits(), "indel".into()));
    };
    let mut types = 0;
    for key in readlist::read_list(list.as_bytes(), ListSource::Inline)? {
        let key = String::from_utf8_lossy(&key).to_ascii_lowercase();
        types |= match key.as_str() {
            "indel" => VariantType::INDEL,
            "mnp" => VariantType::MNP,
            "bnd" => VariantType::BND,
            "other" => VariantType::OTHER,
            "overlap" => VariantType::OVERLAP,
            _ => bail!("Could not parse \"{key}\" in \"--SnpGap {s}\""),
        }
        .bits();
    }
    Ok((gap, types, list.to_owned()))
}

/// The `-g/--SnpGap` and `-G/--IndelGap` settings.
struct GapOptions {
    snp_gap: i64,
    snp_gap_type: u32,
    /// The types as given, for the FILTER description.
    snp_gap_str: String,
    indel_gap: i64,
}

impl GapOptions {
    fn parse(args: &Args) -> anyhow::Result<Self> {
        let (snp_gap, snp_gap_type, snp_gap_str) = match &args.snp_gap {
            Some(s) => parse_snp_gap(s)?,
            None => (0, 0, String::new()),
        };
        let indel_gap = match &args.indel_gap {
            Some(s) => {
                let (v, end) = pasta::variant::number::strtol(s.as_bytes(), 0);
                if end != s.len() {
                    bail!("Could not parse argument: --IndelGap {s}");
                }
                v
            }
            None => 0,
        };
        Ok(Self {
            snp_gap,
            snp_gap_type,
            snp_gap_str,
            indel_gap,
        })
    }

    fn is_set(&self) -> bool {
        self.snp_gap != 0 || self.indel_gap != 0
    }

    /// Adds the `SnpGap` and `IndelGap` FILTER lines, returning their IDs, and warns when
    /// `soft`, the soft filter name when no expression takes it, is another name.
    fn add_filters(&self, header: &mut Header, soft: Option<&str>) -> (i32, i32) {
        if let Some(soft) = soft
            && soft != "+"
        {
            let mut names = String::new();
            if self.snp_gap != 0 {
                names.push_str("\"SnpGap\"");
            }
            if self.indel_gap != 0 {
                if !names.is_empty() {
                    names.push_str(" and ");
                }
                names.push_str("\"IndelGap\"");
            }
            // bcftools compares the soft filter with the first name, up to its length.
            let n = names.len() - 2;
            let first = &names.as_bytes()[1..=n];
            if (0..n).any(|k| soft.as_bytes().get(k) != Some(&first[k])) {
                eprintln!("Warning: using {names} filter name instead of \"{soft}\"");
            }
        }
        let mut ids = (-1, -1);
        if self.snp_gap != 0 {
            let line = format!(
                "##FILTER=<ID=SnpGap,Description=\"SNP within {} bp of {}\">",
                self.snp_gap, self.snp_gap_str
            );
            let _ = header.append(line.as_bytes());
            ids.0 = header.id(b"SnpGap").map_or(-1, |id| id as i32);
        }
        if self.indel_gap != 0 {
            let line = format!(
                "##FILTER=<ID=IndelGap,Description=\"Indel within {} bp of an indel\">",
                self.indel_gap
            );
            let _ = header.append(line.as_bytes());
            ids.1 = header.id(b"IndelGap").map_or(-1, |id| id as i32);
        }
        ids
    }
}

/// The `--mask`/`-M` regions given (the later one), and whether they name a file (any `-M`,
/// as bcftools remembers it).
fn mask_list<'a>(args: &'a Args, matches: &ArgMatches) -> Option<(&'a str, ListSource)> {
    let list = match (args.mask.last(), args.mask_file.last()) {
        (Some(m), Some(f)) => Some(
            if last_index(matches, "mask") > last_index(matches, "mask_file") {
                m
            } else {
                f
            },
        ),
        (m, f) => m.or(f),
    };
    let source = if args.mask_file.is_empty() {
        ListSource::Inline
    } else {
        ListSource::File
    };
    list.map(|l| (l.as_str(), source))
}

/// Reads the mask regions (bcftools' `regidx_init` and `regidx_init_string`).
fn load_mask(list: &str, source: ListSource, overlap: Overlap) -> anyhow::Result<MaskSetting> {
    let (negate, body) = match list.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, list),
    };
    let mask = match source {
        ListSource::File => Mask::from_file(body)
            .map_err(|_| anyhow::anyhow!("Could not initialize the mask: {list}"))?,
        ListSource::Inline => Mask::from_regions(body.as_bytes()),
    };
    Ok(MaskSetting {
        mask,
        negate,
        overlap,
    })
}

/// Opens the input with the regions (`-r/-R` or positional) and targets (`-t/-T`).
fn open_reader(
    args: &Args,
    matches: &ArgMatches,
    fname: &str,
    positional_regions: &[String],
    threads: usize,
) -> anyhow::Result<SyncedReader> {
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
    let mut regions = None;
    if let Some((spec, source)) = later(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    ) {
        let mut r = Regions::new(spec.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = regions_overlap;
        regions = Some(r);
    } else if !positional_regions.is_empty() {
        let joined = positional_regions.join(",");
        let mut r = Regions::new(joined.as_bytes(), ListSource::Inline)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: (null)"))?;
        r.overlap = regions_overlap;
        regions = Some(r);
    }
    let mut targets = None;
    let mut targets_exclude = false;
    if let Some((spec, source)) = later(
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
    let options = SyncedOptions {
        regions,
        targets,
        targets_exclude,
        apply_filters: None,
        threads,
    };
    SyncedReader::open(fname, options).map_err(|e| {
        let name = if fname == "-" {
            "standard input"
        } else {
            fname
        };
        anyhow::anyhow!("Failed to read from {name}: {e}")
    })
}

/// The options of `filter` that are checked before the input is opened.
struct Options<'a> {
    gap: GapOptions,
    annot: AnnotMode,
    filter_expr: Option<(&'a str, FilterLogic)>,
    set_gts: Option<SetGts>,
    threads: usize,
    mask_overlap: Overlap,
    mask_list: Option<(&'a str, ListSource)>,
    write_index: Option<output::WriteIndex>,
}

impl<'a> Options<'a> {
    /// Parses the options in bcftools' order, so that the first bad one is reported.
    fn parse(args: &'a Args, matches: &ArgMatches) -> anyhow::Result<Self> {
        let gap = GapOptions::parse(args)?;
        let mut annot = AnnotMode::default();
        for m in &args.mode {
            annot.reset |= m.contains('x');
            annot.add |= m.contains('+');
        }
        let filter_expr = FilterLogic::single(&args.include, &args.exclude)?;
        let set_gts = match args.set_gts.as_deref() {
            None => None,
            Some(".") => Some(SetGts::Missing),
            Some("0") => Some(SetGts::Ref),
            Some(s) => bail!("The argument to -S not recognised: {s}"),
        };
        // hts_set_threads ignores a negative count.
        let threads = output::threads(args.threads.as_deref()).unwrap_or(0);
        let mask_overlap = match args.mask_overlap.as_deref() {
            None | Some("1") => Overlap::Record,
            Some("0") => Overlap::Pos,
            Some("2") => Overlap::Variant,
            Some(s) => bail!("Could not parse: --mask-overlap {s}"),
        };
        let mask_list = mask_list(args, matches);
        // Checked in bcftools' option loop, before the input is opened.
        output::requested_type(&output_settings(args))?;
        let write_index = args
            .write_index
            .as_deref()
            .map(output::parse_write_index)
            .transpose()?;
        Ok(Self {
            gap,
            annot,
            filter_expr,
            set_gts,
            threads,
            mask_overlap,
            mask_list,
            write_index,
        })
    }
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    if args.long_help {
        eprintln!("filter: unrecognized option '--help'");
    }
    if args.help || args.long_help {
        return Ok(usage());
    }
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    let opts = Options::parse(&args, matches)?;

    let (fname, positional_regions) = match args.inputs.split_first() {
        Some((f, rest)) => (f.as_str(), rest),
        None => {
            if io::stdin().is_terminal() {
                return Ok(usage());
            }
            ("-", &[][..])
        }
    };
    if opts.mask_list.is_some() && args.soft_filter.is_none() {
        bail!("The option --soft-filter is required with --mask and --mask-file options");
    }
    let mut reader = open_reader(&args, matches, fname, positional_regions, opts.threads)?;

    // init_data
    let out_name = args
        .output
        .as_deref()
        .map_or_else(|| "-".to_owned(), |p| p.display().to_string());
    let (out_type, level) = output::resolve(&output_settings(&args), Some(Path::new(&out_name)))?;
    let out_path = args.output.as_deref().map(output::data_path);
    let mut writer = match Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, opts.threads),
    ) {
        Ok(w) => w,
        Err(e) => bail!("Can't write to \"{out_name}\": {}", strerror(&e)),
    };
    let mask = match opts.mask_list {
        Some((list, source)) => Some(load_mask(list, source, opts.mask_overlap)?),
        None => None,
    };
    let ids = add_filter_lines(&args, &opts, reader.header_mut())?;
    if !args.no_version {
        version::append_version(reader.header_mut(), "filter")?;
    }
    let filter = match opts.filter_expr {
        Some((expr, logic)) => Some((Filter::init(reader.header(), expr)?, logic)),
        None => None,
    };

    if let Some(fmt) = opts.write_index {
        let ok = output::index_request(fmt, out_type, args.output.as_deref())
            .is_some_and(|req| writer.request_index(req).is_ok());
        if !ok {
            bail!("Error: failed to initialise index for {out_name}");
        }
    }
    if writer.write_header(reader.header()).is_err() {
        bail!("[main_vcffilter] Error: cannot write the header to {out_name}");
    }
    writer.start_index(reader.header())?;

    let mut settings = Settings {
        filter,
        mask,
        soft_filter: args.soft_filter.is_some(),
        annot: opts.annot,
        set_gts: opts.set_gts,
        flt_pass: ids.pass,
        flt_fail: ids.fail,
    };
    let gap = &opts.gap;
    let gaps = gap.is_set().then(|| {
        GapBuffer::new(GapFilters {
            snp_gap: gap.snp_gap,
            snp_gap_type: gap.snp_gap_type,
            snp_gap_id: ids.snp_gap,
            indel_gap: gap.indel_gap,
            indel_gap_id: ids.indel_gap,
            soft: args.soft_filter.is_some(),
        })
    });

    let result = filter_records(&mut reader, &mut writer, &mut settings, gaps, &out_name);
    finish(&mut reader, writer, result, &out_name)?;
    Ok(ExitCode::SUCCESS)
}

/// The IDs of the FILTER values that `filter` sets, -1 for those it does not.
struct FilterIds {
    pass: i32,
    /// The soft filter of the expression or the mask.
    fail: i32,
    snp_gap: i32,
    indel_gap: i32,
}

/// Adds the `##FILTER` lines of the soft filter and the gap filters to the header.
fn add_filter_lines(
    args: &Args,
    opts: &Options<'_>,
    header: &mut Header,
) -> anyhow::Result<FilterIds> {
    let pass = header.id(b"PASS").map_or(-1, |id| id as i32);
    debug_assert_eq!(pass, 0, "PASS is the first filter in BCF");
    let mut fail = -1;
    if let Some(soft) = &args.soft_filter
        && (opts.filter_expr.is_some() || opts.mask_list.is_some())
    {
        fail = add_soft_filter(header, soft, opts.filter_expr)?;
    }
    let (snp_gap, indel_gap) = if opts.gap.is_set() {
        // Without an expression, the soft filter name is the gap filters'.
        let soft = args
            .soft_filter
            .as_deref()
            .filter(|_| opts.filter_expr.is_none());
        opts.gap.add_filters(header, soft)
    } else {
        (-1, -1)
    };
    Ok(FilterIds {
        pass,
        fail,
        snp_gap,
        indel_gap,
    })
}

/// Reports the result of the main loop and closes the output.
fn finish(
    reader: &mut SyncedReader,
    writer: Writer,
    result: Result<(), FilterError>,
    out_name: &str,
) -> anyhow::Result<()> {
    // bcftools exits from error() on any error, without flushing its output buffers; a
    // missing BGZF EOF marker is reported after reading everything.
    let error = match result {
        Ok(()) => reader
            .take_pending_error()
            .map(|e| anyhow::anyhow!("Error: {e}")),
        Err(FilterError::Synced(e)) => Some(anyhow::anyhow!("Error: {e}")),
        Err(FilterError::Fatal(e)) => Some(e),
    };
    if let Some(e) = error {
        writer.abandon();
        return Err(e);
    }
    match writer.finish() {
        Ok(()) => Ok(()),
        Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
        Err(WriteError::Deferred(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("[main_vcffilter] Error: cannot write to {out_name}")
        })),
        Err(WriteError::Io(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("[main_vcffilter] Error: close failed .. {out_name}")
        })),
        Err(e) => Err(anyhow::Error::new(e).context(format!(
            "[main_vcffilter] Error: close failed .. {out_name}"
        ))),
    }
}

/// Prints the usage and returns bcftools' exit status for it.
fn usage() -> ExitCode {
    eprintln!("{USAGE}");
    ExitCode::from(1)
}

/// Adds the `##FILTER` line of the soft filter, inventing a `FilterN` name for `+`; returns
/// its ID.
fn add_soft_filter(
    header: &mut Header,
    soft: &str,
    expr: Option<(&str, FilterLogic)>,
) -> anyhow::Result<i32> {
    let name = if soft == "+" {
        let mut i = 0;
        loop {
            i += 1;
            let name = format!("Filter{i}");
            let taken = header
                .id(name.as_bytes())
                .is_some_and(|id| header.id_exists(HeaderRecordType::Filter, id));
            if !taken {
                break name;
            }
        }
    } else {
        soft.to_owned()
    };
    let description = match expr {
        Some((expr, logic)) => {
            let mut d = String::from(if logic == FilterLogic::Include {
                "Set if not true: "
            } else {
                "Set if true: "
            });
            for c in expr.chars() {
                if c == '"' {
                    d.push('\\');
                }
                d.push(c);
            }
            d
        }
        None => "Record masked by region".to_owned(),
    };
    let line = format!("##FILTER=<ID={name},Description=\"{description}\">");
    if header.append(line.as_bytes()).is_err() {
        bail!("Failed to append header line: {line}");
    }
    Ok(header.id(name.as_bytes()).map_or(-1, |id| id as i32))
}

/// A reader error (reported with its exit status) or a fatal error.
enum FilterError {
    Synced(SyncedError),
    Fatal(anyhow::Error),
}

impl From<anyhow::Error> for FilterError {
    fn from(e: anyhow::Error) -> Self {
        Self::Fatal(e)
    }
}

fn write_one(
    writer: &mut Writer,
    header: &Header,
    rec: &mut Record,
    out_name: &str,
) -> anyhow::Result<()> {
    match writer.write_record(header, rec) {
        Ok(()) => Ok(()),
        Err(WriteError::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => {
            Err(anyhow::Error::new(e))
        }
        Err(_) => bail!("[main_vcffilter] Error: cannot write to {out_name}"),
    }
}

/// The main loop of bcftools filter.
fn filter_records(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    s: &mut Settings,
    mut gaps: Option<GapBuffer>,
    out_name: &str,
) -> Result<(), FilterError> {
    let mut rec = Box::new(Record::new());
    let mut gts = GtBuffers::default();
    loop {
        match reader.next_record_boxed(&mut rec) {
            Ok(true) => {}
            Ok(false) => break,
            Err(e) => return Err(FilterError::Synced(e)),
        }
        let header = reader.header();
        let logic = s.filter.as_ref().map(|(_, l)| *l);
        let mut pass = true;
        let mut smpl_pass: Option<&[bool]> = None;
        if let Some((filter, logic)) = &mut s.filter {
            pass = filter.test(header, &mut rec)?;
            if *logic == FilterLogic::Exclude {
                pass = !pass;
            }
            smpl_pass = filter.sample_pass();
        }
        if let Some(m) = &s.mask {
            let (beg, end) = m.overlap.interval(&mut rec);
            let chrom = header.contig_name(rec.rid()).unwrap_or(b"(null)");
            let mpass = m.mask.overlaps(chrom, beg, end) == m.negate;
            pass &= mpass;
        }
        if !(s.soft_filter || s.set_gts.is_some() || pass) {
            continue;
        }
        if pass {
            rec.unpack(Unpack::FLT);
            if s.annot.reset || rec.filters().is_empty() {
                rec.add_filter(s.flt_pass);
            }
        } else if s.soft_filter {
            if s.annot.add {
                rec.add_filter(s.flt_fail);
            } else {
                rec.update_filter(&[s.flt_fail]);
            }
        }
        if let Some(mode) = s.set_gts {
            set_genotypes(header, &mut rec, pass, smpl_pass, logic, mode, &mut gts)?;
        }
        match &mut gaps {
            None => write_one(writer, header, &mut rec, out_name)?,
            Some(g) => {
                let line = std::mem::replace(&mut *rec, g.spare());
                g.push(header, Some(line), &mut |r| {
                    write_one(writer, header, r, out_name)
                })?;
            }
        }
    }
    // bcftools reports a reading error (a missing BGZF EOF marker) before it flushes the
    // records held back for the gap filters.
    if let Some(e) = reader.take_pending_error() {
        return Err(FilterError::Synced(e));
    }
    if let Some(g) = &mut gaps {
        let header = reader.header();
        g.push(header, None, &mut |r| {
            write_one(writer, header, r, out_name)
        })?;
    }
    Ok(())
}

/// Buffers for [`set_genotypes`].
#[derive(Default)]
struct GtBuffers {
    gt: Vec<i32>,
    ac: Vec<i32>,
    /// The indexes of the missing AC values.
    ac_missing: Vec<usize>,
}

/// Sets the genotypes of failing samples to missing or reference, updating INFO/AC and AN
/// (bcftools' `set_genotypes`). A missing AN, or a missing AC value, stays missing; bcftools
/// counts down from the missing value's integer, writing numbers such as 2147483646.
fn set_genotypes(
    header: &Header,
    rec: &mut Record,
    pass_site: bool,
    smpl_pass: Option<&[bool]>,
    logic: Option<FilterLogic>,
    mode: SetGts,
    buf: &mut GtBuffers,
) -> anyhow::Result<()> {
    let nsmpl = header.sample_count();
    if nsmpl == 0 {
        return Ok(());
    }
    let exclude = logic == Some(FilterLogic::Exclude);
    if let Some(sp) = smpl_pass {
        let npass = sp.iter().filter(|&&p| p).count();
        if npass == nsmpl && logic == Some(FilterLogic::Include) {
            return Ok(());
        }
        if npass == 0 && exclude {
            return Ok(());
        }
    } else if pass_site {
        return Ok(());
    }

    let mut an = 0;
    let has_an =
        matches!(rec.get_info_i32(header, b"AN", &mut buf.ac), Ok(1)) && buf.ac[0] != INT32_MISSING;
    if has_an {
        an = buf.ac[0];
    }
    let n_allele = rec.allele_count();
    let has_ac = matches!(rec.get_info_i32(header, b"AC", &mut buf.ac), Ok(n) if n + 1 == n_allele);
    buf.ac_missing.clear();
    if has_ac && buf.ac[..n_allele - 1].contains(&INT32_MISSING) {
        set_aside_missing(&mut buf.ac[..n_allele - 1], &mut buf.ac_missing);
    }
    match rec.get_genotypes(header, &mut buf.gt) {
        Ok(n) => {
            let ngts = n / nsmpl;
            // bcf_gt_missing, or bcf_gt_unphased(0).
            let new_gt = match mode {
                SetGts::Missing => 0,
                SetGts::Ref => 2,
            };
            let ac: &mut [i32] = if has_ac {
                &mut buf.ac[..n_allele - 1]
            } else {
                &mut []
            };
            for i in 0..nsmpl {
                if let Some(sp) = smpl_pass
                    && sp[i] != exclude
                {
                    continue;
                }
                for gt in &mut buf.gt[i * ngts..(i + 1) * ngts] {
                    if *gt == INT32_VECTOR_END {
                        break;
                    }
                    // The allele's index plus one, 0 for a missing allele.
                    let code = *gt >> 1;
                    if code == 0 {
                        if mode == SetGts::Ref {
                            an = an.wrapping_add(1);
                        }
                    } else {
                        if mode == SetGts::Missing {
                            an = an.wrapping_sub(1);
                        }
                        // An ALT allele of the record leaves AC (REF, code 1, and alleles
                        // past the record's or negative codes index nothing).
                        if let Some(c) = ac.get_mut((code as u32).wrapping_sub(2) as usize) {
                            *c = c.wrapping_sub(1);
                        }
                    }
                    *gt = new_gt;
                }
            }
            rec.update_genotypes(header, &buf.gt)?;
        }
        Err(e) => set_unreadable_genotypes(header, rec, e)?,
    }
    if has_an {
        rec.update_info_i32(header, b"AN", &[an])?;
    }
    if has_ac {
        for &i in &buf.ac_missing {
            buf.ac[i] = INT32_MISSING;
        }
        rec.update_info_i32(header, b"AC", &buf.ac[..n_allele - 1])?;
    }
    Ok(())
}

/// Sets the missing values of `ac` to 0, for [`set_genotypes`] to count down from, and lists
/// their indexes in `missing` to be put back. Kept out of line: AC is rarely missing.
#[cold]
#[inline(never)]
fn set_aside_missing(ac: &mut [i32], missing: &mut Vec<usize>) {
    for (i, v) in ac.iter_mut().enumerate() {
        if *v == INT32_MISSING {
            missing.push(i);
            *v = 0;
        }
    }
}

/// Passes on bcftools' handling of a record whose GT cannot be read: it divides HTSlib's
/// negative return code by the number of samples and sets GT from that many values per
/// sample. A count of 0 removes GT; a negative one (with one to three samples) sets an empty
/// GT, printed as `.`, when the header defines GT. A record with fewer samples than the header
/// is an error.
fn set_unreadable_genotypes(header: &Header, rec: &mut Record, e: GetError) -> anyhow::Result<()> {
    let code: usize = match e {
        GetError::NotInHeader => 1,
        GetError::WrongType => 2,
        GetError::Absent => 3,
        GetError::TooFewSamples => return Err(e.into()),
    };
    let nsmpl = header.sample_count();
    if code < nsmpl {
        rec.remove_format(header, b"GT")?;
    } else if e != GetError::NotInHeader {
        rec.update_genotypes(header, &vec![0; nsmpl])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use pasta::variant::vcf;

    use super::*;

    /// Runs [`set_genotypes`] on each line with every sample failing, returning the INFO
    /// column and the genotypes.
    fn set_all(mode: SetGts, lines: &[&str]) -> Vec<String> {
        set_all_with(mode, lines, |_, _| {})
    }

    /// [`set_all`], calling `edit` on each record before setting its genotypes.
    fn set_all_with(
        mode: SetGts,
        lines: &[&str],
        mut edit: impl FnMut(&Header, &mut Record),
    ) -> Vec<String> {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n##contig=<ID=1>\n\
              ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"c\">\n\
              ##INFO=<ID=AN,Number=1,Type=Integer,Description=\"n\">\n\
              ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        )
        .expect("the header parses");
        let mut rec = Record::new();
        let mut buf = GtBuffers::default();
        lines
            .iter()
            .map(|line| {
                let mut line = line.as_bytes().to_vec();
                vcf::parse(&mut header, &mut line, &mut rec).expect("the line parses");
                edit(&header, &mut rec);
                set_genotypes(&header, &mut rec, false, None, None, mode, &mut buf)
                    .expect("the genotypes are set");
                let mut out = Vec::new();
                vcf::format_record(&header, &rec, &mut out).expect("the record formats");
                let out = String::from_utf8(out).expect("UTF-8 output");
                out.trim_end()
                    .split('\t')
                    .skip(7)
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .collect()
    }

    #[test]
    fn missing_ac_and_an_stay_missing() {
        let lines = [
            "1\t1\t.\tA\tC\t.\t.\tAC=.;AN=.\tGT\t1/1\t./1",
            "1\t2\t.\tA\tC,G\t.\t.\tAC=.,3;AN=4\tGT\t1/2\t2/2",
            "1\t3\t.\tA\tC,G\t.\t.\tAC=2,.;AN=.\tGT\t1/0\t./1",
            "1\t4\t.\tA\tC\t.\t.\tAC=3;AN=4\tGT\t1/1\t0/1",
            "1\t5\t.\tA\tC\t.\t.\tAC=3;AN=4\tGT\t1/3\t1/.",
        ];
        assert_eq!(
            set_all(SetGts::Missing, &lines),
            [
                "AC=.;AN=. GT ./. ./.",
                "AC=.,0;AN=0 GT ./. ./.",
                "AC=0,.;AN=. GT ./. ./.",
                "AC=0;AN=0 GT ./. ./.",
                "AC=1;AN=1 GT ./. ./.",
            ]
        );
        assert_eq!(
            set_all(SetGts::Ref, &lines),
            [
                "AC=.;AN=. GT 0/0 0/0",
                "AC=.,0;AN=4 GT 0/0 0/0",
                "AC=0,.;AN=. GT 0/0 0/0",
                "AC=0;AN=4 GT 0/0 0/0",
                "AC=1;AN=5 GT 0/0 0/0",
            ]
        );
    }

    #[test]
    fn ac_ending_early_is_left_alone() {
        // AC with fewer values than ALT alleles (a vector end in BCF) is not updated.
        let line = ["1\t1\t.\tA\tC,G\t.\t.\tAN=4\tGT\t1/2\t2/2"];
        let out = set_all_with(SetGts::Missing, &line, |h, rec| {
            rec.update_info_i32(h, b"AC", &[1, INT32_VECTOR_END])
                .expect("AC is set");
        });
        assert_eq!(out, ["AN=0;AC=1 GT ./. ./."]);
    }
}
