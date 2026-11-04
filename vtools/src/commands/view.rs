//! `vtools view`: convert, subset, and filter VCF/BCF files (bcftools view).

mod subset;

use std::{
    io::{self, IsTerminal},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, HeaderRecordType, Record,
    io::{OutputType, WriteError, Writer},
    log,
    readlist::{self, ListSource},
    synced::{Overlap, Regions, SyncedOptions, SyncedReader},
    utils::CalcAcError,
};

use self::subset::{AlleleKind, Subsetter, VariantTypes};
use crate::options::later;
use crate::{
    commands::head::{apply_verbosity, strtol_whole},
    filter::{self, FilterLogic},
    output::{self, OutputSetting},
    pipeline,
    raw_input::strerror,
    version,
};

const USAGE: &str = "
About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.
Usage:   vtools view [options] <in.vcf.gz> [region1 [...]]

Output options:
    -G, --drop-genotypes              Drop individual genotype information (after subsetting if -s option set)
    -h, --header-only                 Print only the header in VCF output (equivalent to vtools head)
    -H, --no-header                   Suppress the header in VCF output
        --with-header                 Print both header and records in VCF output [default]
    -l, --compression-level [0-9]     Compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]
        --no-version                  Do not append version and command line to the header
    -o, --output FILE                 Output file name [stdout]
    -O, --output-type u|b|v|z[0-9]    u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -r, --regions REGION              Restrict to comma-separated list of regions
    -R, --regions-file FILE           Restrict to regions listed in FILE
        --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -t, --targets [^]REGION           Similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix
    -T, --targets-file [^]FILE        Similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix
        --targets-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT                 Use multithreading with INT worker threads [0]
        --verbosity INT               Verbosity level

Subset options:
    -A, --trim-unseen-allele          Remove '<*>' or '<NON_REF>' at variant (-A) or at all (-AA) sites
    -a, --trim-alt-alleles            Trim ALT alleles not seen in the genotype fields (or their subset with -s/-S)
    -I, --no-update                   Do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)
    -s, --samples [^]LIST             Comma separated list of samples to include (or exclude with \"^\" prefix). Be careful
                                        when combining filtering with sample subsetting as filtering comes (usually) first.
                                        If unsure, split sample subsetting and filtering in two commands, using -Ou when piping.
    -S, --samples-file [^]FILE        File of samples to include (or exclude with \"^\" prefix)
        --force-samples               Only warn about unknown subset samples

Filter options:
    -c/C, --min-ac/--max-ac INT[:TYPE]     Minimum/maximum count for non-reference (nref), 1st alternate (alt1), least frequent
                                               (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
    -f,   --apply-filters LIST             Require at least one of the listed FILTER strings (e.g. \"PASS,.\")
    -g,   --genotype [^]hom|het|miss       Require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude such sites
    -i/e, --include/--exclude EXPR         Select/exclude sites for which the expression is true (see man page for details)
    -k/n, --known/--novel                  Select known/novel sites only (ID is not/is '.')
    -m/M, --min-alleles/--max-alleles INT  Minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
    -p/P, --phased/--exclude-phased        Select/exclude sites where all samples are phased
    -q/Q, --min-af/--max-af FLOAT[:TYPE]   Minimum/maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent
                                               (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
    -u/U, --uncalled/--exclude-uncalled    Select/exclude sites without a called genotype
    -v/V, --types/--exclude-types LIST     Select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other,overlap [null]
    -x/X, --private/--exclude-private      Select/exclude sites where the non-reference alleles are exclusive (private) to the subset samples
    -W,   --write-index[=FMT]              Automatically index the output files [off]

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin), then regions to view (needs an index).
    inputs: Vec<String>,

    /// Drop individual genotype information (after subsetting if -s is set).
    #[arg(short = 'G', long = "drop-genotypes")]
    drop_genotypes: bool,

    /// Print only the header in VCF output.
    #[arg(short = 'h', long = "header-only", action = clap::ArgAction::Count)]
    header_only: u8,

    /// Suppress the header in VCF output.
    #[arg(short = 'H', long = "no-header", action = clap::ArgAction::Count)]
    no_header: u8,

    /// Print both header and records in VCF output [default].
    #[arg(long = "with-header", action = clap::ArgAction::Count)]
    with_header: u8,

    /// Compression level: 0 uncompressed, 1 best speed, 9 best compression.
    #[arg(
        short = 'l',
        long = "compression-level",
        value_name = "0-9",
        action = clap::ArgAction::Append,
        allow_hyphen_values = true
    )]
    compression_level: Vec<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Output file name [stdout].
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "output-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(
        short = 'O',
        long = "output-type",
        value_name = "u|b|v|z[0-9]",
        action = clap::ArgAction::Append,
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

    /// Restrict to regions listed in FILE.
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

    /// Similar to -r but streams rather than index-jumps. Exclude regions with "^" prefix.
    #[arg(
        short = 't',
        long = "targets",
        value_name = "[^]REGION",
        allow_hyphen_values = true
    )]
    targets: Option<String>,

    /// Similar to -R but streams rather than index-jumps. Exclude regions with "^" prefix.
    #[arg(
        short = 'T',
        long = "targets-file",
        value_name = "[^]FILE",
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
    #[arg(long = "verbosity", value_name = "INT", allow_hyphen_values = true)]
    verbosity: Option<String>,

    /// Remove '<*>' or '<NON_REF>' at variant (-A) or at all (-AA) sites.
    #[arg(short = 'A', long = "trim-unseen-allele", action = clap::ArgAction::Count)]
    trim_unseen_allele: u8,

    /// Trim ALT alleles not seen in the genotype fields (or their subset with -s/-S).
    #[arg(short = 'a', long = "trim-alt-alleles")]
    trim_alt_alleles: bool,

    /// Do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN).
    #[arg(short = 'I', long = "no-update")]
    no_update: bool,

    /// Comma separated list of samples to include (or exclude with "^" prefix).
    #[arg(
        short = 's',
        long = "samples",
        value_name = "[^]LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to include (or exclude with "^" prefix).
    #[arg(
        short = 'S',
        long = "samples-file",
        value_name = "[^]FILE",
        allow_hyphen_values = true
    )]
    samples_file: Option<String>,

    /// Only warn about unknown subset samples.
    #[arg(long = "force-samples")]
    force_samples: bool,

    /// Minimum count for non-reference (nref), 1st alternate (alt1), least frequent (minor),
    /// most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref].
    #[arg(
        short = 'c',
        long = "min-ac",
        value_name = "INT[:TYPE]",
        allow_hyphen_values = true
    )]
    min_ac: Option<String>,

    /// Maximum count, as for --min-ac.
    #[arg(
        short = 'C',
        long = "max-ac",
        value_name = "INT[:TYPE]",
        allow_hyphen_values = true
    )]
    max_ac: Option<String>,

    /// Require at least one of the listed FILTER strings (e.g. "PASS,.").
    #[arg(
        short = 'f',
        long = "apply-filters",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    apply_filters: Option<String>,

    /// Require one or more hom/het/missing genotype or, if prefixed with "^", exclude such
    /// sites.
    #[arg(
        short = 'g',
        long = "genotype",
        value_name = "[^]hom|het|miss",
        allow_hyphen_values = true
    )]
    genotype: Option<String>,

    /// Select sites for which the expression is true.
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Select known sites only (ID is not '.').
    #[arg(short = 'k', long = "known")]
    known: bool,

    /// Select novel sites only (ID is '.').
    #[arg(short = 'n', long = "novel")]
    novel: bool,

    /// Minimum number of alleles listed in REF and ALT.
    #[arg(
        short = 'm',
        long = "min-alleles",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    min_alleles: Option<String>,

    /// Maximum number of alleles listed in REF and ALT.
    #[arg(
        short = 'M',
        long = "max-alleles",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    max_alleles: Option<String>,

    /// Select sites where all samples are phased.
    #[arg(short = 'p', long = "phased")]
    phased: bool,

    /// Exclude sites where all samples are phased.
    #[arg(short = 'P', long = "exclude-phased")]
    exclude_phased: bool,

    /// Minimum frequency, as for --min-ac.
    #[arg(
        short = 'q',
        long = "min-af",
        value_name = "FLOAT[:TYPE]",
        allow_hyphen_values = true
    )]
    min_af: Option<String>,

    /// Maximum frequency, as for --min-ac.
    #[arg(
        short = 'Q',
        long = "max-af",
        value_name = "FLOAT[:TYPE]",
        allow_hyphen_values = true
    )]
    max_af: Option<String>,

    /// Select sites without a called genotype.
    #[arg(short = 'u', long = "uncalled")]
    uncalled: bool,

    /// Exclude sites without a called genotype.
    #[arg(short = 'U', long = "exclude-uncalled")]
    exclude_uncalled: bool,

    /// Select comma-separated list of variant types: snps,indels,mnps,ref,bnd,other,overlap.
    #[arg(
        short = 'v',
        long = "types",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    types: Option<String>,

    /// Exclude comma-separated list of variant types.
    #[arg(
        short = 'V',
        long = "exclude-types",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    exclude_types: Option<String>,

    /// Select sites where the non-reference alleles are exclusive (private) to the subset
    /// samples.
    #[arg(short = 'x', long = "private")]
    private: bool,

    /// Exclude sites where the non-reference alleles are exclusive to the subset samples.
    #[arg(short = 'X', long = "exclude-private")]
    exclude_private: bool,

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
    #[arg(long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// `-g` genotype requirements (bcftools' `GT_NEED_*`/`GT_NO_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GtRequirement {
    NeedHom,
    NeedHet,
    NoHom,
    NoHet,
    NeedMissing,
    NoMissing,
}

/// Include or exclude (bcftools' `FLT_INCLUDE`/`FLT_EXCLUDE`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Selection {
    Include,
    Exclude,
}

/// What `-h`, `-H`, and `--with-header` ask to write: `-h` keeps the records from being read,
/// `-H` the header from being written, and `--with-header` undoes both.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum HeaderOutput {
    /// The header and the records (the default).
    HeaderAndRecords,
    /// The records alone (`-H`).
    Records,
    /// The header alone (`-h`).
    Header,
    /// Nothing (`-h` and `-H`).
    Nothing,
}

impl HeaderOutput {
    /// Reads `-h`, `-H`, and `--with-header` in command-line order.
    fn from_matches(matches: &ArgMatches) -> Self {
        let mut print_header = true;
        let mut header_only = false;
        let mut flags: Vec<(usize, char)> = Vec::new();
        for (id, c) in [
            ("header_only", 'h'),
            ("no_header", 'H'),
            ("with_header", 'w'),
        ] {
            // Counted flags have a default with an index; only real occurrences count.
            if matches.get_count(id) > 0
                && let Some(idx) = matches.indices_of(id)
            {
                flags.extend(idx.map(|i| (i, c)));
            }
        }
        flags.sort_unstable();
        for (_, c) in flags {
            match c {
                'h' => header_only = true,
                'H' => print_header = false,
                _ => {
                    print_header = true;
                    header_only = false;
                }
            }
        }
        match (print_header, header_only) {
            (true, false) => Self::HeaderAndRecords,
            (false, false) => Self::Records,
            (true, true) => Self::Header,
            (false, true) => Self::Nothing,
        }
    }

    fn writes_header(self) -> bool {
        matches!(self, Self::HeaderAndRecords | Self::Header)
    }

    fn reads_records(self) -> bool {
        matches!(self, Self::HeaderAndRecords | Self::Records)
    }
}

/// Collects `-O` and `-l` in command-line order.
fn output_settings(args: &Args, matches: &ArgMatches) -> Vec<OutputSetting> {
    let mut settings: Vec<(usize, OutputSetting)> = Vec::new();
    if let Some(idx) = matches.indices_of("output_type") {
        for (i, v) in idx.zip(&args.output_type) {
            settings.push((i, OutputSetting::Type(v.clone())));
        }
    }
    if let Some(idx) = matches.indices_of("compression_level") {
        for (i, v) in idx.zip(&args.compression_level) {
            settings.push((i, OutputSetting::Level(v.clone())));
        }
    }
    settings.sort_by_key(|(i, _)| *i);
    settings.into_iter().map(|(_, s)| s).collect()
}

/// C's `strtol(s, &end, 10)` into an `int`, with nothing after the number (an empty text is 0).
fn parse_int_strict(s: &str, what: &str) -> anyhow::Result<i64> {
    match strtol_whole(s) {
        Some(v) => Ok(i64::from(v as i32)),
        None => bail!("Could not parse argument: --{what} {s}"),
    }
}

/// How a `-c/-C/-q/-Q` threshold's number is read: C's `%d` or `%f`.
#[derive(Clone, Copy)]
enum NumberKind {
    Int,
    Float,
}

/// Parses `INT[:TYPE]` or `FLOAT[:TYPE]` as bcftools does with `sscanf(s, "%d:%8s")` (or
/// `%f`), falling back to the number alone.
///
/// The number is returned as an `f64`: a `%d` value already truncated to a C `int`, or a `%f`
/// value (decimal, hex, `inf`, or `nan`, as `strtod` reads them).
fn parse_count_type(s: &str, what: &str, kind: NumberKind) -> anyhow::Result<(f64, AlleleKind)> {
    let bytes = s.as_bytes();
    let start = bytes.iter().take_while(|c| c.is_ascii_whitespace()).count();
    let (value, end) = match kind {
        NumberKind::Float => pasta::variant::number::strtod(bytes, start),
        NumberKind::Int => {
            let (v, end) = pasta::variant::number::strtol(bytes, start);
            (f64::from(v as i32), end)
        }
    };
    if end == start {
        bail!("Error: Could not parse --{what} {s}");
    }
    let allele = match s[end..].strip_prefix(':') {
        Some(rest) => {
            let word: String = rest
                .trim_start()
                .chars()
                .take_while(|c| !c.is_whitespace())
                .take(8)
                .collect();
            if word.is_empty() {
                AlleleKind::NonRef
            } else {
                AlleleKind::parse(&word)?
            }
        }
        None => AlleleKind::NonRef,
    };
    Ok((value, allele))
}

fn parse_overlap(s: Option<&str>, default: Overlap, what: &str) -> anyhow::Result<Overlap> {
    match s {
        None => Ok(default),
        Some(s) => match Overlap::parse(s) {
            Some(o) => Ok(o),
            None => bail!("Could not parse: --{what} {s}"),
        },
    }
}

/// The selection of an include flag and its exclude flag (`-x/-X`, `-u/-U`, `-p/-P`), which
/// cannot both be given.
fn selection(include: bool, exclude: bool, flags: &str) -> anyhow::Result<Option<Selection>> {
    if include && exclude {
        bail!("Only one of {flags} can be given.");
    }
    Ok(if include {
        Some(Selection::Include)
    } else if exclude {
        Some(Selection::Exclude)
    } else {
        None
    })
}

/// Parses `-g [^]hom|het|miss`.
fn parse_genotype(g: &str) -> anyhow::Result<GtRequirement> {
    Ok(match g.to_ascii_lowercase().as_str() {
        "hom" => GtRequirement::NeedHom,
        "het" => GtRequirement::NeedHet,
        "miss" => GtRequirement::NeedMissing,
        "^hom" => GtRequirement::NoHom,
        "^het" => GtRequirement::NoHet,
        "^miss" => GtRequirement::NoMissing,
        _ => bail!(
            "The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"{g}\"."
        ),
    })
}

/// Reads the filter and subset options into the record selection, failing at the first invalid
/// one in the order bcftools checks them. The samples and variant types are set later, once the
/// header is read.
fn subsetter(args: &Args) -> anyhow::Result<Subsetter> {
    let mut sub = Subsetter::default();
    if let Some(m) = &args.min_alleles {
        sub.min_alleles = parse_int_strict(m, "min-alleles")?;
    }
    if let Some(m) = &args.max_alleles {
        sub.max_alleles = parse_int_strict(m, "max-alleles")?;
    }
    sub.known = args.known;
    sub.novel = args.novel;
    if args.trim_alt_alleles {
        sub.trim_alts = true;
        sub.calc_ac = true;
    }
    sub.trim_star_allele = args.trim_unseen_allele;
    sub.update_info = !args.no_update;
    sub.sites_only = args.drop_genotypes;
    if let Some(s) = &args.min_ac {
        let (v, k) = parse_count_type(s, "min-ac", NumberKind::Int)?;
        sub.min_ac = Some((v as i32, k));
        sub.calc_ac = true;
    }
    if let Some(s) = &args.max_ac {
        let (v, k) = parse_count_type(s, "max-ac", NumberKind::Int)?;
        sub.max_ac = Some((v as i32, k));
        sub.calc_ac = true;
    }
    if let Some(s) = &args.min_af {
        let (v, k) = parse_count_type(s, "min-af", NumberKind::Float)?;
        sub.min_af = Some((v as f32, k));
        sub.calc_ac = true;
    }
    if let Some(s) = &args.max_af {
        let (v, k) = parse_count_type(s, "max-af", NumberKind::Float)?;
        sub.max_af = Some((v as f32, k));
        sub.calc_ac = true;
    }
    sub.private_vars = selection(args.private, args.exclude_private, "-x or -X")?;
    sub.uncalled = selection(args.uncalled, args.exclude_uncalled, "-u or -U")?;
    sub.phased = selection(args.phased, args.exclude_phased, "-p or -P")?;
    if sub.private_vars.is_some() || sub.uncalled.is_some() {
        sub.calc_ac = true;
    }
    if let Some(g) = &args.genotype {
        sub.gt_type = Some(parse_genotype(g)?);
    }
    Ok(sub)
}

/// The reader's regions and targets: `-r/-R` (else the regions after the input) and `-t/-T`,
/// the later of each pair winning.
fn synced_options(
    args: &Args,
    matches: &ArgMatches,
    positional_regions: &[String],
    threads: i64,
) -> anyhow::Result<SyncedOptions> {
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
    let region_spec = later(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    );
    let mut regions = None;
    if let Some((spec, source)) = region_spec {
        let mut r = Regions::new(spec.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = regions_overlap;
        regions = Some(r);
    } else if !positional_regions.is_empty() {
        let joined = positional_regions.join(",");
        let mut r = Regions::new(joined.as_bytes(), ListSource::Inline)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {joined}"))?;
        r.overlap = Overlap::Record;
        regions = Some(r);
    }
    let target_spec = later(
        matches,
        ("targets", args.targets.as_deref()),
        ("targets_file", args.targets_file.as_deref()),
    );
    let mut targets = None;
    let mut targets_exclude = false;
    if let Some((spec, source)) = target_spec {
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

    if threads < 0 || !output::threads_start(threads as usize) {
        bail!("Failed to create threads");
    }
    Ok(SyncedOptions {
        regions,
        targets,
        targets_exclude,
        apply_filters: args.apply_filters.as_ref().map(|f| f.as_bytes().to_vec()),
        threads: threads as usize,
    })
}

/// Opens the input, reporting a failure as bcftools does.
fn open_reader(fname: &str, options: SyncedOptions) -> anyhow::Result<SyncedReader> {
    match SyncedReader::open(fname, options) {
        Ok(r) => Ok(r),
        Err(e) => {
            let name = if fname == "-" {
                "standard input"
            } else {
                fname
            };
            bail!("Failed to read from {name}: {e}");
        }
    }
}

/// Opens the output (`-o`, `-O`, `-l`); returns it with its type.
fn create_writer(
    args: &Args,
    matches: &ArgMatches,
    threads: usize,
) -> anyhow::Result<(Writer, OutputType)> {
    let (out_type, level) =
        output::resolve(&output_settings(args, matches), args.output.as_deref())?;
    let out_path = args.output.as_deref().map(output::data_path);
    match Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, threads),
    ) {
        Ok(w) => Ok((w, out_type)),
        Err(e) => {
            let name = out_path
                .as_deref()
                .map_or_else(|| "-".into(), |p| p.display().to_string());
            let msg = strerror(&e);
            log::error(
                "hts_open_format",
                format_args!("Failed to open file \"{name}\" : {msg}"),
            );
            bail!(
                "{}: {msg}",
                args.output
                    .as_deref()
                    .map_or_else(|| "-".into(), |p| p.display().to_string())
            )
        }
    }
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    let header_output = HeaderOutput::from_matches(matches);

    let filter_expr = FilterLogic::single(&args.include, &args.exclude)?;

    // getopt rejects a bad -O, -l, or -W value before anything is read.
    output::requested_type(&output_settings(&args, matches))?;
    if let Some(wi) = &args.write_index {
        output::parse_write_index(wi)?;
    }

    // C's `(int) strtol(optarg, 0, 0)`; a negative count fails in `bcf_sr_set_threads`.
    let threads = args
        .threads
        .as_deref()
        .map_or(0, |t| i64::from(output::thread_count(t)));

    let mut sub = subsetter(&args)?;
    // -s and -S: the later one wins.
    let sample_spec = later(
        matches,
        ("samples", args.samples.as_deref()),
        ("samples_file", args.samples_file.as_deref()),
    );
    if sample_spec.is_some() && sub.update_info {
        sub.calc_ac = true;
    }

    let (fname, positional_regions) = match args.inputs.split_first() {
        Some((fname, regions)) => (fname.as_str(), regions),
        None if io::stdin().is_terminal() => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
        None => ("-", &[][..]),
    };
    let options = synced_options(&args, matches, positional_regions, threads)?;
    let threads = options.threads;
    let mut reader = open_reader(fname, options)?;

    // init_data
    if sub.calc_ac && sub.update_info {
        let h = reader.header_mut();
        h.append(b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">")
            .context("[init_data] Failed to add \"AC\" INFO header")?;
        h.append(b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">")
            .context("[init_data] Failed to add \"AN\" INFO header")?;
    }
    if !args.no_version {
        version::append_version(reader.header_mut(), "view")?;
    }

    let mut sample_names: Vec<Vec<u8>> = Vec::new();
    if let Some((spec, source)) = sample_spec {
        let unknown = if args.force_samples {
            UnknownSample::Skip
        } else {
            UnknownSample::Fail
        };
        sample_names = select_samples(reader.header(), spec, source, unknown)?;
        if sample_names.is_empty() {
            eprintln!("Warn: subsetting has removed all samples");
            sub.sites_only = true;
        }
    }

    if let Some(t) = &args.types {
        if args.exclude_types.is_some() {
            bail!("Error: only supply one of --include-types, --exclude-types options");
        }
        sub.include = VariantTypes::parse(t)?;
    } else if let Some(t) = &args.exclude_types {
        sub.exclude = VariantTypes::parse(t)?;
    }

    let (mut writer, out_type) = create_writer(&args, matches, threads)?;
    let result = write_records(
        &mut reader,
        &mut writer,
        &mut sub,
        filter_expr,
        &sample_names,
        header_output,
        out_type,
        threads,
        &args,
    );
    finish(&mut reader, writer, result, header_output, &args)
}

/// Ends the run once the records are written: reports a fatal error or the reader's pending
/// one, and closes the output.
fn finish(
    reader: &mut SyncedReader,
    writer: Writer,
    result: anyhow::Result<u8>,
    header_output: HeaderOutput,
    args: &Args,
) -> anyhow::Result<ExitCode> {
    let code = match result {
        Ok(code) => code,
        Err(e) => {
            // bcftools exits from error() without flushing its output buffers.
            writer.abandon();
            // HTSlib's `bcf_calc_ac` reports a bad genotype allele and calls exit(1).
            if e.downcast_ref::<CalcAcError>().is_some() {
                return Ok(ExitCode::from(1));
            }
            return Err(e);
        }
    };
    // With -h, bcftools never reads a record, nor reports the reader's pending error.
    let code = if code == 0 && header_output.reads_records() {
        match reader.take_pending_error() {
            Some(e) => {
                eprintln!("Error: {e}");
                e.code() as u8
            }
            None => 0,
        }
    } else {
        code
    };
    match writer.finish() {
        Ok(()) => {}
        Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
        // Records the compression threads failed to write before closing.
        Err(WriteError::Deferred(e)) => return Err(write_error(WriteError::Io(e), args)),
        // bcftools' message, without the I/O error (which `{:#}` would print twice).
        Err(WriteError::Io(e)) => {
            return Err(output::write_failure(e, || {
                anyhow::anyhow!("Error: close failed {}", close_name(args))
            }));
        }
        Err(e) => {
            return Err(
                anyhow::Error::new(e).context(format!("Error: close failed {}", close_name(args)))
            );
        }
    }
    Ok(ExitCode::from(code))
}

/// Writes the header and the selected records; returns the exit status of a reader error (or
/// 0), or a fatal error.
#[expect(
    clippy::too_many_arguments,
    reason = "the state of bcftools view's main"
)]
fn write_records(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    sub: &mut Subsetter,
    filter_expr: Option<(&str, FilterLogic)>,
    sample_names: &[Vec<u8>],
    header_output: HeaderOutput,
    out_type: OutputType,
    threads: usize,
    args: &Args,
) -> anyhow::Result<u8> {
    let header = reader.header();
    let hnull = if sub.sites_only {
        let (mut h, _) = header
            .subset(&[])
            .map_err(|_| anyhow::anyhow!("Error occurred while subsetting samples"))?;
        h.remove(HeaderRecordType::Format, None);
        Some(h)
    } else {
        None
    };
    let mut hsub = None;
    if !sample_names.is_empty() {
        let refs: Vec<&[u8]> = sample_names.iter().map(Vec::as_slice).collect();
        let (h, imap) = header
            .subset(&refs)
            .map_err(|_| anyhow::anyhow!("Error occurred while subsetting samples"))?;
        if sample_names.len() != h.sample_count() {
            for (i, &m) in imap.iter().enumerate() {
                if m < 0 {
                    bail!(
                        "Error: No such sample: \"{}\"",
                        String::from_utf8_lossy(&sample_names[i])
                    );
                }
            }
        }
        sub.imap = imap;
        hsub = Some(h);
    }
    sub.sample_subset = !sample_names.is_empty();

    let mut filter = match filter_expr {
        Some((expr, logic)) => Some((filter::init(header, expr)?, logic)),
        None => None,
    };

    let out_hdr: &Header = hnull.as_ref().or(hsub.as_ref()).unwrap_or(header);
    if let Some(wi) = &args.write_index {
        let fmt = output::parse_write_index(wi)?;
        let ok = match output::index_request(fmt, out_type, args.output.as_deref()) {
            Some(req) => writer.request_index(req).is_ok(),
            None => false,
        };
        if !ok {
            bail!("Error: failed to initialise index for {}", out_name(args));
        }
    }
    if header_output.writes_header() {
        if writer.write_header(out_hdr).is_err() {
            bail!("[main_vcfview] Error: cannot write to {}", out_name(args));
        }
    } else if matches!(out_type, OutputType::Bcf | OutputType::BcfUncompressed) {
        bail!("BCF output requires header, cannot proceed with -H");
    }
    writer.start_index(out_hdr)?;

    if !header_output.reads_records() {
        return Ok(0);
    }

    let hsub_ref = hsub.as_ref();
    let hnull_ref = hnull.as_ref();
    let subset_mode = hnull.is_some() || hsub.is_some();
    let plain = filter.is_none() && !sub.is_active() && reader.is_plain_stream() && !subset_mode;
    let streamed = if plain && threads > 0 {
        pipeline::stream(reader, writer, threads)
    } else {
        Ok(pipeline::Streamed::NoThreads)
    };
    // Without threads, or when none could be started.
    let result = match streamed {
        Ok(pipeline::Streamed::Done) => Ok(()),
        Ok(pipeline::Streamed::NoThreads) => {
            serial(reader, writer, sub, &mut filter, hsub_ref, hnull_ref)
        }
        Err(e) => Err(e),
    };
    let code = match result {
        Ok(()) => 0u8,
        Err(ViewError::Synced(e)) => {
            eprintln!("Error: {e}");
            e.code() as u8
        }
        Err(ViewError::Write(e)) => return Err(write_error(e, args)),
        Err(ViewError::Fatal(e)) => return Err(e),
    };
    Ok(code)
}

/// bcftools view's error for a record it cannot write (a closed pipe is not reported).
fn write_error(e: WriteError, args: &Args) -> anyhow::Error {
    output::write_failure(e, || {
        anyhow::anyhow!("[main_vcfview] Error: cannot write to {}", out_name(args))
    })
}

fn out_name(args: &Args) -> String {
    args.output
        .as_deref()
        .map_or_else(|| "(null)".into(), |p| p.display().to_string())
}

/// The output's name in bcftools view's error for a failed close, which names standard output
/// `stdout`.
fn close_name(args: &Args) -> String {
    args.output
        .as_deref()
        .map_or_else(|| "stdout".into(), |p| p.display().to_string())
}

/// An error while viewing records: a reader error (reported with its exit status), a record
/// that could not be written (reported with the output's name), or a fatal one.
pub enum ViewError {
    Synced(pasta::variant::synced::SyncedError),
    Write(WriteError),
    Fatal(anyhow::Error),
}

impl From<anyhow::Error> for ViewError {
    fn from(e: anyhow::Error) -> Self {
        Self::Fatal(e)
    }
}

/// Selects, subsets, and writes the records on this thread (the main loop of bcftools view).
/// The output header is `hnull` (sites only), else `hsub` (a sample subset), else the input's.
fn serial(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    sub: &mut Subsetter,
    filter: &mut Option<(Box<dyn filter::RecordFilter + Send>, FilterLogic)>,
    hsub: Option<&Header>,
    hnull: Option<&Header>,
) -> Result<(), ViewError> {
    let subset_mode = hnull.is_some() || hsub.is_some();
    let mut rec = Box::new(Record::new());
    loop {
        match reader.next_record_boxed(&mut rec) {
            Ok(true) => {}
            Ok(false) => return Ok(()),
            Err(e) => return Err(ViewError::Synced(e)),
        }
        if !rec.errors().is_empty() && subset_mode {
            return Err(ViewError::Fatal(anyhow::anyhow!(
                "Undefined tags in the header, cannot proceed in the sample subset mode."
            )));
        }
        let header = reader.header();
        if !sub.keep(header, hsub, &mut rec, filter)? {
            continue;
        }
        let out_hdr = hnull.or(hsub).unwrap_or(header);
        writer
            .write_record(out_hdr, &mut rec)
            .map_err(ViewError::Write)?;
    }
}

/// What `-s/-S` does with a sample the header lacks.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum UnknownSample {
    /// Fail (the default).
    Fail,
    /// Warn and skip it (`--force-samples`).
    Skip,
}

/// Builds the subset sample list as bcftools view's `init_data` does.
fn select_samples(
    header: &Header,
    spec: &str,
    source: ListSource,
    unknown: UnknownSample,
) -> anyhow::Result<Vec<Vec<u8>>> {
    let (exclude, body) = match spec.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, spec),
    };
    let Ok(list) = readlist::read_list(body.as_bytes(), source) else {
        bail!("Could not read the list: \"{body}\"");
    };
    let mut out = Vec::new();
    if exclude {
        let mut excluded = std::collections::HashSet::new();
        for s in list {
            if header.sample_id(&s).is_none() {
                match unknown {
                    UnknownSample::Skip => eprintln!(
                        "Warn: exclude called for sample that does not exist in header: \"{}\"... skipping",
                        String::from_utf8_lossy(&s)
                    ),
                    UnknownSample::Fail => bail!(
                        "Error: exclude called for sample that does not exist in header: \"{}\". Use \"--force-samples\" to ignore this error.",
                        String::from_utf8_lossy(&s)
                    ),
                }
            }
            excluded.insert(s);
        }
        for s in header.samples() {
            if !excluded.contains(&s[..]) {
                out.push(s.to_vec());
            }
        }
    } else {
        for s in list {
            if header.sample_id(&s).is_none() {
                match unknown {
                    UnknownSample::Skip => {
                        eprintln!(
                            "Warn: subset called for sample that does not exist in header: \"{}\"... skipping",
                            String::from_utf8_lossy(&s)
                        );
                        continue;
                    }
                    UnknownSample::Fail => bail!(
                        "Error: subset called for sample that does not exist in header: \"{}\". Use \"--force-samples\" to ignore this error.",
                        String::from_utf8_lossy(&s)
                    ),
                }
            }
            out.push(s);
        }
    }
    Ok(out)
}
