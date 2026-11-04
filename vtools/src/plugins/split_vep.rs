//! `+split-vep`: query structured annotations such as INFO/CSQ written by VEP, INFO/BCSQ by
//! `bcftools csq`, or INFO/ANN by SnpEff (bcftools' `plugins/split-vep.c`): select
//! transcripts by type and consequence severity, extract subfields into INFO tags (`-c`) or
//! print them with a format string (`-f`, through [`crate::convert`]).

mod init;
mod process;

use std::{
    collections::{HashMap, HashSet},
    io::{self, IsTerminal},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::bail;
use clap::{ArgMatches, CommandFactory, FromArgMatches, Parser, builder::TypedValueParser as _};
use pasta::variant::synced::SyncedReader;

use super::Exit;
use crate::{
    commands::head::apply_verbosity,
    convert::Convert,
    filter::{Filter, FilterLogic, bre::Bre},
    shortopt,
};

/// The description `plugin -lv` prints.
pub const ABOUT: &str = "Query structured annotations such as INFO/CSQ produced by VEP, BCSQ by bcftools/csq, or ANN by SnpEff.\n";

/// The default severity scale (`-S -` prints it).
const DEFAULT_SEVERITY: &str = "\
# Default consequence substrings ordered in ascending order by severity.\n\
# Consequences with the same severity can be put on the same line in arbitrary order.\n\
# See also https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html\n\
intergenic\n\
feature_truncation feature_elongation\n\
regulatory\n\
TF_binding_site TFBS\n\
downstream upstream\n\
non_coding_transcript non_coding\n\
intron NMD_transcript\n\
non_coding_transcript_exon\n\
5_prime_utr 3_prime_utr\n\
coding_sequence mature_miRNA\n\
stop_retained start_retained synonymous\n\
incomplete_terminal_codon\n\
splice_region\n\
missense inframe protein_altering\n\
transcript_amplification\n\
exon_loss\n\
disruptive\n\
start_lost stop_lost stop_gained frameshift\n\
splice_acceptor splice_donor\n\
transcript_ablation\n";

/// The default subfield types (`--columns-types -` prints them).
const DEFAULT_COLUMN_TYPES: &str = "\
# Default CSQ subfield types, unlisted fields are type String.\n\
# Note that the name search is done using regular expressions, with\n\
# \"^\" and \"$\" appended automatically\n\
DISTANCE                   Integer\n\
Distance                   Integer\n\
STRAND                     Integer\n\
TSL                        Integer\n\
GENE_PHENO                 Integer\n\
HGVS_OFFSET                Integer\n\
.*_POPS                    String\n\
AF                         Float\n\
.*_AF                      Float\n\
MAX_AF_.*                  Float\n\
MOTIF_POS                  Integer\n\
MOTIF_SCORE_CHANGE         Float\n\
existing_InFrame_oORFs     Integer\n\
existing_OutOfFrame_oORFs  Integer\n\
existing_uORFs             Integer\n\
SpliceAI_pred_DP_.*        Integer\n\
SpliceAI_pred_DS_.*        Float\n";

/// The plugin's usage text.
pub fn usage() -> String {
    "\n\
About: Query structured annotations such as INFO/BCSQ produced by bcftools/csq, INFO/CSQ produced by VEP, or\n\
\x20  INFO/ANN produced by SnpEff.\n\
\x20  For more information and pointers see http://samtools.github.io/bcftools/howtos/plugin.split-vep.html\n\
Usage: vtools +split-vep [Plugin Options]\n\
Plugin options:\n\
\x20  -a, --annotation STR            INFO annotation to parse; auto-detected in precedence order CSQ, BCSQ, ANN\n\
\x20  -A, --all-fields DELIM          Output all fields replacing the selected annotation tag, e.g. %CSQ, %BCSQ or %ANN,\n\
\x20                                    in the -f filtering expression using the field delimiter DELIM. This can be \"tab\",\n\
\x20                                    \"space\" or an arbitrary string.\n\
\x20  -c, --columns [LIST|-][:TYPE]   Extract the fields listed either as 0-based indexes or names, \"-\" to extract all\n\
\x20                                    fields. See --columns-types for the defaults. Supported types are String/Str,\n\
\x20                                    Integer/Int and Float/Real. Unlisted fields are set to String. Existing header\n\
\x20                                    definitions will not be overwritten, remove first with `vtools annotate -x`\n\
\x20      --columns-types -|FILE      Pass \"-\" to print the default -c types or FILE to override the presets\n\
\x20      --csq-field STR,            \n\
\x20      --consequence-field STR     Use STR as the consequence field name [Consequence|Annotation]\n\
\x20  -d, --duplicate                 Output per transcript/allele consequences on a new line rather rather than\n\
\x20                                    as comma-separated fields on a single line\n\
\x20  -f, --format STR                Create non-VCF output; similar to `vtools query -f` but drops lines w/o consequence\n\
\x20  -g, --gene-list [+]FILE         Consider only features listed in FILE, or prioritize if FILE is prefixed with \"+\"\n\
\x20      --gene-list-fields LIST     Fields to match against by the -g list, by default gene names [SYMBOL,Gene,gene,Gene_Name,Gene_ID]\n\
\x20  -H, --print-header              Print header, -HH to omit column indices\n\
\x20  -l, --list                      Parse the VCF header and list the annotation fields\n\
\x20  -p, --annot-prefix STR          Before doing anything else, prepend STR to all CSQ fields to avoid tag name conflicts\n\
\x20  -s, --select TR:CSQ[:PRN]       Select transcripts to extract by type and/or consequence severity, see also -S and -x\n\
\x20                                    TR, filter transcripts:   all,worst,primary,pick,mane,EXPRESSION [all]\n\
\x20                                    CSQ, filter consequences: any,missense,missense+,etc [any]\n\
\x20                                    PRN, print consequences:  all,worst [all]\n\
\x20                                  TR transcript selection\n\
\x20                                      all       .. list all transcripts\n\
\x20                                      worst     .. list only one transcript with the worst consequence (see -S)\n\
\x20                                      primary   .. list transcripts marked as CANONICAL=YES\n\
\x20                                      pick      .. as PICK=1\n\
\x20                                      mane      .. as MANE_SELECT!=\"\"\n\
\x20                                    or an EXPRESSION in the form of \"<FIELD><OPERATOR><VALUE>\", where\n\
\x20                                      FIELD     .. field name (e.g. \"CANONICAL\")\n\
\x20                                      OPERATOR  .. string comparison (=,!=), regex matching (~,!~)\n\
\x20                                      VALUE     .. required string value (e.g. \"YES\")\n\
\x20                                  CSQ consequence filtering, selects transcripts by CSQ severity (see \"-S -\")\n\
\x20                                      missense  .. selects only transcripts with a missense variant\n\
\x20                                      missense+ .. transcripts with a missense consequence or more severe\n\
\x20                                  PRN controls what consequences are actually printed\n\
\x20                                      all   .. print all consequences if multiple per transcript (e.g., start_lost&splice_region)\n\
\x20                                      worst .. print the worst consequence per transcript (e.g., start_lost above)\n\
\x20  -S, --severity -|FILE           Pass \"-\" to print the default severity scale or FILE to override\n\
\x20                                    the default scale\n\
\x20  -u, --allow-undef-tags          Print \".\" for undefined tags\n\
\x20  -x, --drop-sites                Drop sites without consequences after -s is applied (the default with -f)\n\
\x20  -X, --keep-sites                Do not drop sites without consequences (the default without -f)\n\
Common options:\n\
\x20  -e, --exclude EXPR              Exclude sites and samples for which the expression is true\n\
\x20  -i, --include EXPR              Include sites and samples for which the expression is true\n\
\x20      --no-version                Do not append version and command line to the header\n\
\x20  -o, --output FILE               Output file name [stdout]\n\
\x20  -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n\
\x20  -r, --regions REG               Restrict to comma-separated list of regions\n\
\x20  -R, --regions-file FILE         Restrict to regions listed in a file\n\
\x20      --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n\
\x20  -t, --targets REG               Similar to -r but streams rather than index-jumps\n\
\x20  -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n\
\x20      --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n\
\x20      --threads INT               Use multithreading with INT worker threads [0] (a vtools addition)\n\
\x20  -v, --verbosity INT             Verbosity level\n\
\x20  -W, --write-index[=FMT]         Automatically index the output files [off]\n\
\n\
Examples:\n\
\x20  # List available fields of the INFO/CSQ or INFO/BCSQ annotation\n\
\x20  vtools +split-vep -l file.vcf.gz\n\
\n\
\x20  # List available fields of INFO/BCSQ when both INFO/CSQ and INFO/BCSQ are present\n\
\x20  vtools +split-vep -l file.vcf.gz -a BCSQ\n\
\n\
\x20  # List the default severity scale\n\
\x20  vtools +split-vep -S -\n\
\n\
\x20  # Extract Consequence, IMPACT and gene SYMBOL of the most severe consequence into\n\
\x20  # INFO annotations starting with the prefix \"vep\". For brevity, the columns can\n\
\x20  # be given also as 0-based indexes\n\
\x20  vtools +split-vep -c Consequence,IMPACT,SYMBOL -s worst -p vep file.vcf.gz\n\
\x20  vtools +split-vep -c 1-3 -s worst -p vep file.vcf.gz\n\
\n\
\x20  # Same as above but use the text output of the \"vtools query\" format\n\
\x20  vtools +split-vep -s worst -f '%CHROM %POS %Consequence %IMPACT %SYMBOL\\n' file.vcf.gz\n\
\n\
\x20  # Print all subfields (tab-delimited) in place of %CSQ, each consequence on a new line\n\
\x20  vtools +split-vep -f '%CHROM %POS %CSQ\\n' -d -A tab file.vcf.gz\n\
\n\
\x20  # Extract gnomAD_AF subfield into a new INFO/gnomAD_AF annotation of Type=Float so that\n\
\x20  # numeric filtering can be used.\n\
\x20  vtools +split-vep -c gnomAD_AF:Float file.vcf.gz -i'gnomAD_AF<0.001'\n\
\n\
\x20  # Similar to above, but add the annotation only if the consequence severity is missense\n\
\x20  # or equivalent. In order to drop sites with different consequences completely, provide\n\
\x20  # the -x option. See the online documentation referenced above for more examples.\n\
\x20  vtools +split-vep -c gnomAD_AF:Float -s :missense    file.vcf.gz\n\
\x20  vtools +split-vep -c gnomAD_AF:Float -s :missense -x file.vcf.gz\n\
\n\
\x20  See also http://samtools.github.io/bcftools/howtos/plugin.split-vep.html\n\
\n"
    .to_owned()
}

/// The options (bcftools' getopt string `o:O:i:e:r:R:t:T:lS:s:c:p:a:f:dA:xXuHg:W::v:`).
#[derive(Debug, Parser)]
#[command(
    no_binary_name = true,
    disable_help_flag = true,
    disable_version_flag = true,
    args_override_self = true
)]
struct Args {
    /// The input file (`-` or none for standard input).
    inputs: Vec<String>,
    /// Drop sites without consequences after -s is applied (the default with -f).
    #[arg(short = 'x', long = "drop-sites")]
    drop_sites: bool,
    /// Do not drop sites without consequences (the default without -f).
    #[arg(short = 'X', long = "keep-sites")]
    keep_sites: bool,
    /// Output all fields in place of the annotation tag of -f, with this delimiter.
    #[arg(short = 'A', long = "all-fields", allow_hyphen_values = true)]
    all_fields: Option<String>,
    /// Output per transcript/allele consequences on a new line.
    #[arg(short = 'd', long = "duplicate")]
    duplicate: bool,
    /// Create non-VCF output, as `query -f`.
    #[arg(short = 'f', long = "format", allow_hyphen_values = true)]
    format: Option<String>,
    /// Consider only features listed in FILE, or prioritize them with a "+" prefix.
    #[arg(short = 'g', long = "gene-list", allow_hyphen_values = true)]
    gene_list: Option<String>,
    /// Fields to match against by the -g list.
    #[arg(long = "gene-list-fields", allow_hyphen_values = true)]
    gene_list_fields: Option<String>,
    /// INFO annotation to parse.
    #[arg(short = 'a', long = "annotation", allow_hyphen_values = true)]
    annotation: Option<String>,
    /// Prepend STR to all CSQ fields.
    #[arg(short = 'p', long = "annot-prefix", allow_hyphen_values = true)]
    annot_prefix: Option<String>,
    /// The consequence field name.
    #[arg(
        long = "csq-field",
        visible_alias = "consequence-field",
        allow_hyphen_values = true
    )]
    csq_field: Option<String>,
    /// Extract the fields listed as 0-based indexes or names, "-" for all.
    #[arg(short = 'c', long = "columns", allow_hyphen_values = true)]
    columns: Option<String>,
    /// "-" to print the default -c types, or a file to override them.
    #[arg(long = "columns-types", allow_hyphen_values = true)]
    columns_types: Option<String>,
    /// Select transcripts by type and/or consequence severity.
    #[arg(short = 's', long = "select", allow_hyphen_values = true)]
    select: Option<String>,
    /// "-" to print the default severity scale, or a file to override it.
    #[arg(short = 'S', long = "severity", allow_hyphen_values = true)]
    severity: Option<String>,
    /// List the annotation fields.
    #[arg(short = 'l', long = "list")]
    list: bool,
    /// Include sites and samples for which the expression is true.
    #[arg(short = 'i', long = "include", allow_hyphen_values = true)]
    include: Vec<String>,
    /// Exclude sites and samples for which the expression is true.
    #[arg(short = 'e', long = "exclude", allow_hyphen_values = true)]
    exclude: Vec<String>,
    /// Output file name [stdout].
    #[arg(short = 'o', long = "output", value_parser = clap::builder::OsStringValueParser::new().map(PathBuf::from))]
    output: Option<PathBuf>,
    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", action = clap::ArgAction::Append)]
    output_type: Vec<String>,
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
    /// Print "." for undefined tags.
    #[arg(short = 'u', long = "allow-undef-tags")]
    allow_undef_tags: bool,
    /// Print header, -HH to omit column indices.
    #[arg(short = 'H', long = "print-header", action = clap::ArgAction::Count)]
    print_header: u8,
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
    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,
    /// Not in bcftools' option string: getopt reports `-h` as an invalid option.
    #[arg(short = 'h', hide = true, action = clap::ArgAction::Count)]
    invalid_h: u8,
    /// Not in bcftools' option string: getopt reports `-?` as an invalid option.
    #[arg(short = '?', hide = true, action = clap::ArgAction::Count)]
    invalid_q: u8,
    /// Worker threads for reading, parsing, and writing (a vtools addition).
    #[arg(long = "threads", allow_hyphen_values = true)]
    threads: Option<String>,
}

/// The type of an extracted subfield (`BCF_HT_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ColType {
    Str,
    Int,
    Real,
    Flag,
}

impl ColType {
    /// Parses a `-c NAME:TYPE` type.
    fn parse(s: &str) -> Option<Self> {
        let is = |names: &[&str]| names.iter().any(|n| s.eq_ignore_ascii_case(n));
        if is(&["string", "str"]) {
            Some(Self::Str)
        } else if is(&["float", "real"]) {
            Some(Self::Real)
        } else if is(&["integer", "int"]) {
            Some(Self::Int)
        } else if is(&["flag"]) {
            Some(Self::Flag)
        } else {
            None
        }
    }

    fn header_name(self) -> &'static str {
        match self {
            Self::Str => "String",
            Self::Int => "Integer",
            Self::Real => "Float",
            Self::Flag => "Flag",
        }
    }
}

/// An extracted subfield (bcftools' `annot_t`).
#[derive(Debug)]
struct Annot {
    /// The INFO tag (the subfield's name).
    tag: String,
    /// The subfield index, `None` for the whole annotation (the raw tag).
    idx: Option<usize>,
    ty: ColType,
    /// The output header ID of the tag.
    id: Option<u32>,
    /// The value being built.
    value: Vec<u8>,
}

/// The `--select` transcript expression operator.
#[derive(Debug)]
enum TrOp {
    Eq(Vec<u8>),
    Ne(Vec<u8>),
    Re(Bre),
    NotRe(Bre),
}

/// `--select FIELD<OP>VALUE` (bcftools' `select_tr_t`).
#[derive(Debug)]
struct TrExpr {
    field: String,
    idx: usize,
    op: TrOp,
}

/// Which transcripts are selected (bcftools' `SELECT_TR_*`).
#[derive(Debug)]
enum SelectTr {
    All,
    Worst,
    Expr(TrExpr),
}

/// The consequence severity scale and the terms seen so far (bcftools' `scale` and
/// `csq2severity`).
#[derive(Debug, Default)]
struct Severity {
    scale: Vec<Vec<u8>>,
    map: HashMap<Vec<u8>, i32>,
    key: Vec<u8>,
}

/// `-g`: the genes and the fields to look them up in.
#[derive(Debug)]
struct Genes {
    names: HashSet<Vec<u8>>,
    fields: Vec<usize>,
    prioritize: bool,
}

/// Everything `init_data` sets up.
struct SplitVep {
    vep_tag: String,
    fields: Vec<String>,
    field2idx: HashMap<String, usize>,
    annot_prefix: Option<String>,
    csq_idx: usize,
    /// The consequence field's name as given or detected, for error messages.
    csq_field: String,
    severity: Severity,
    min_severity: i32,
    max_severity: i32,
    prn_csq_worst: bool,
    select_tr: SelectTr,
    annots: Vec<Annot>,
    column_str: Option<String>,
    column_types: Option<String>,
    raw_vep_request: bool,
    filter: Option<(Filter, FilterLogic)>,
    convert: Option<Convert>,
    genes: Option<Genes>,
    duplicate: bool,
    drop_sites: bool,
    has_format: bool,
    warned: Warned,
    csq: process::Csq,
    farr: Vec<f32>,
    iarr: Vec<i32>,
    list_tr: Vec<usize>,
    line: Vec<u8>,
    /// How much of each VCF line is parsed.
    max_unpack: pasta::variant::Unpack,
}

/// Warnings printed only once.
#[derive(Debug, Default)]
struct Warned {
    nfields: bool,
    too_few_fields: bool,
    real: bool,
    int: bool,
}

/// An option bcftools' getopt loop checks as it reads it.
#[derive(Clone, Copy, Debug)]
enum Checked<'a> {
    Shared(super::Checked<'a>),
    /// An option getopt does not know (`-h`, `-?`): bcftools prints the usage.
    Invalid(char),
}

/// Checks the options in command-line order, failing at the first bad one as bcftools'
/// getopt loop does.
fn check_options(args: &Args, matches: &ArgMatches) -> anyhow::Result<()> {
    use super::Checked as S;
    let events = super::in_order(
        matches,
        &[
            ("verbosity", &args.verbosity, |v| {
                Checked::Shared(S::Verbosity(v))
            }),
            ("output_type", &args.output_type, |v| {
                Checked::Shared(S::OutputType(v))
            }),
            ("include", &args.include, |_| Checked::Shared(S::Filter)),
            ("exclude", &args.exclude, |_| Checked::Shared(S::Filter)),
            ("regions_overlap", &args.regions_overlap, |v| {
                Checked::Shared(S::RegionsOverlap(v))
            }),
            ("targets_overlap", &args.targets_overlap, |v| {
                Checked::Shared(S::TargetsOverlap(v))
            }),
            ("write_index", &args.write_index, |v| {
                Checked::Shared(S::WriteIndex(v))
            }),
        ],
        &[
            ("invalid_h", args.invalid_h, Checked::Invalid('h')),
            ("invalid_q", args.invalid_q, Checked::Invalid('?')),
        ],
    );
    let mut filters = 0;
    for event in events {
        match event {
            Checked::Invalid(c) => {
                eprintln!("split-vep: invalid option -- '{c}'");
                return Err(Exit::error(usage()));
            }
            Checked::Shared(S::Verbosity(v)) => apply_verbosity(v)?,
            Checked::Shared(c) => super::check_shared(c, &mut filters)?,
        }
    }
    Ok(())
}

/// Runs `+split-vep` (bcftools' `run`).
///
/// # Errors
///
/// Returns bcftools' fatal errors.
pub fn run(argv: &[String]) -> anyhow::Result<ExitCode> {
    let argv = shortopt::rewrite_strings(&Args::command(), "split-vep", argv)
        .map_err(|msg| Exit::error(format!("{msg}{}", usage())))?;
    let matches = Args::command()
        .try_get_matches_from(argv)
        .map_err(|_| Exit::error(usage()))?;
    let args = Args::from_arg_matches(&matches).map_err(|_| Exit::error(usage()))?;
    check_options(&args, &matches)?;
    if args.print_header > 0 && args.format.is_none() {
        bail!("Error: the -H header printing is supported only with -f");
    }
    if args.all_fields.is_some() && args.format.is_none() {
        bail!("Error: the -A option must be used with -f");
    }
    if matches!(args.severity.as_deref(), Some("?" | "-")) {
        return Err(Exit::error(DEFAULT_SEVERITY));
    }
    if args.columns_types.as_deref() == Some("-") {
        return Err(Exit::error(DEFAULT_COLUMN_TYPES));
    }
    let fname = match args.inputs.as_slice() {
        [] => {
            if io::stdin().is_terminal() {
                return Err(Exit::error(usage()));
            }
            "-".to_owned()
        }
        [f] => f.clone(),
        _ => return Err(Exit::error(usage())),
    };
    let drop_sites = sites_option(&matches);
    let threads = crate::output::threads(args.threads.as_deref()).unwrap_or(0);
    let lists = super::Lists {
        regions: crate::options::later(
            &matches,
            ("regions", args.regions.as_deref()),
            ("regions_file", args.regions_file.as_deref()),
        )
        .map(|(spec, source)| (spec.to_owned(), source)),
        targets: crate::options::later(
            &matches,
            ("targets", args.targets.as_deref()),
            ("targets_file", args.targets_file.as_deref()),
        )
        .map(|(spec, source)| (spec.to_owned(), source)),
    };
    let options = super::synced_options(
        &lists,
        args.regions_overlap.last().map(String::as_str),
        args.targets_overlap.last().map(String::as_str),
        threads,
    )?;
    let mut reader =
        SyncedReader::open(&fname, options).map_err(|e| anyhow::anyhow!("Error: {e}"))?;
    let mut hdr_out = reader.header().dup()?;
    let mut sv = init::init_data(&args, &mut reader, &mut hdr_out)?;
    if args.list {
        let mut out = String::new();
        for (i, f) in sv.fields.iter().enumerate() {
            out.push_str(&format!("{i}\t{f}\n"));
        }
        print!("{out}");
        return Ok(ExitCode::SUCCESS);
    }
    let drop_sites = sv.resolve_drop_sites(&args, drop_sites)?;
    sv.drop_sites = drop_sites;
    process::run(sv, &args, reader, hdr_out, threads)?;
    Ok(ExitCode::SUCCESS)
}

/// The later of `-x` (drop) and `-X` (keep), if any.
fn sites_option(matches: &ArgMatches) -> Option<bool> {
    let x = matches
        .get_flag("drop_sites")
        .then(|| matches.index_of("drop_sites"))
        .flatten();
    let big_x = matches
        .get_flag("keep_sites")
        .then(|| matches.index_of("keep_sites"))
        .flatten();
    match (x, big_x) {
        (Some(a), Some(b)) => Some(a > b),
        (Some(_), None) => Some(true),
        (None, Some(_)) => Some(false),
        (None, None) => None,
    }
}

impl SplitVep {
    /// bcftools' checks of `-c/-f/-s` and the `-x/-X` default.
    fn resolve_drop_sites(&self, args: &Args, drop_sites: Option<bool>) -> anyhow::Result<bool> {
        let mut drop_sites = drop_sites;
        if args.format.is_none()
            && self.column_str.is_none()
            && self.min_severity == process::SELECT_CSQ_ANY
            && self.max_severity == process::SELECT_CSQ_ANY
        {
            if !matches!(self.select_tr, SelectTr::Expr(_)) {
                bail!(
                    "Error: none of the -c,-f,-s options was given, why not use \"vtools view\" instead?"
                );
            }
            match drop_sites {
                // `-s mane` and the like are expected to keep only sites with such transcripts.
                None => drop_sites = Some(true),
                Some(false) => bail!(
                    "Error: the option -X has no effect without -c,-f, why not use \"vtools view\" instead?"
                ),
                Some(true) => {}
            }
        }
        Ok(drop_sites.unwrap_or(args.format.is_some()))
    }
}
