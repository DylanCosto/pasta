//! The command line of `gtcheck`: its options and bcftools' option loop, which reads them
//! in command-line order.

use anyhow::bail;
use clap::{ArgMatches, parser::ValueSource};
use pasta::variant::{
    number::{strtod, strtol},
    readlist::ListSource,
    synced::Overlap,
};

use super::kernels::Tag;
use crate::commands::head::apply_verbosity;

pub(super) const USAGE: &str = "
About:   Check sample identity. With no -g BCF given, multi-sample cross-check is performed.
Usage:   vtools gtcheck [options] [-g <genotypes.vcf.gz>] <query.vcf.gz>

Options:
        --distinctive-sites            Find sites that can distinguish between at least NUM sample pairs.
                  NUM[,MEM[,TMP]]          If the number is smaller or equal to 1, it is interpreted as the fraction of pairs.
                                           The optional MEM string sets the maximum memory used for in-memory sorting [500M]
                                           and TMP is a prefix of temporary files used by external sorting [/tmp/bcftools-gtcheckXXXXXX]
        --dry-run                      Stop after first record to estimate required time
    -E, --error-probability INT        Phred-scaled probability of genotyping error, 0 for faster but less accurate results [40]
    -e, --exclude [qry|gt]:EXPR        Exclude sites for which the expression is true
    -g, --genotypes FILE               Genotypes to compare against
    -H, --homs-only                    Homozygous genotypes only, useful with low coverage data (requires -g)
    -i, --include [qry|gt]:EXPR        Include sites for which the expression is true
        --keep-refs                    Keep monoallelic sites with no alternate allele
        --n-matches INT                Print only top INT matches for each sample (sorted by average score), 0 for unlimited.
                                           Use negative value to sort by HWE probability rather than by discordance [0]
        --no-HWE-prob                  Disable calculation of HWE probability
    -o, --output FILE                  Write output to a file [standard output]
    -O, --output-type t|z              t: plain tab-delimited text output, z: compressed [t]
    -p, --pairs LIST                   Comma-separated sample pairs to compare (qry,gt[,qry,gt..] with -g or qry,qry[,qry,qry..] w/o)
    -P, --pairs-file FILE              File with tab-delimited sample pairs to compare (qry,gt with -g or qry,qry w/o)
    -r, --regions REGION               Restrict to comma-separated list of regions
    -R, --regions-file FILE            Restrict to regions listed in a file
        --regions-overlap 0|1|2        Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --samples [qry|gt]:LIST        List of query or -g samples, \"-\" to select all samples (by default all samples are compared)
    -S, --samples-file [qry|gt]:FILE   File with the query or -g samples to compare
    -t, --targets REGION               Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE            Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2        Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT                  Use multithreading with INT worker threads [0]
    -u, --use TAG1[,TAG2]              Which tag to use in the query file (TAG1) and the -g file (TAG2) [PL,GT]
    -v, --verbosity INT                Verbosity level
Examples:
   # Check discordance of all samples from B against all samples in A
   vtools gtcheck -g A.bcf B.bcf

   # Limit comparisons to the given list of samples
   vtools gtcheck -s gt:a1,a2,a3 -s qry:b1,b2 -g A.bcf B.bcf

   # Compare only two pairs a1,b1 and a1,b2
   vtools gtcheck -p a1,b1,a1,b2 -g A.bcf B.bcf

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// The query VCF/BCF file (standard input when absent and not a terminal).
    pub(super) inputs: Vec<String>,

    /// Output comparison for all sites (not implemented in bcftools 1.24: an error).
    #[arg(short = 'a', long = "all-sites", action = clap::ArgAction::Count)]
    all_sites: u8,

    /// Min inter- and max intra-sample error (not implemented in bcftools 1.24: an error).
    #[arg(short = 'c', long = "cluster", value_name = "MIN,MAX", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    cluster: Vec<String>,

    /// Find sites that can distinguish between at least NUM sample pairs.
    #[arg(long = "distinctive-sites", value_name = "NUM[,MEM[,TMP]]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    distinctive_sites: Vec<String>,

    /// Stop after first record to estimate required time.
    #[arg(long = "dry-run", action = clap::ArgAction::Count)]
    dry_run: u8,

    /// Phred-scaled probability of genotyping error, 0 for faster but less accurate results [40].
    #[arg(short = 'E', long = "error-probability", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    error_probability: Vec<String>,

    /// Exclude sites for which the expression is true ([qry|gt]:EXPR).
    #[arg(short = 'e', long = "exclude", value_name = "EXPR", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    exclude: Vec<String>,

    /// Genotypes to compare against.
    #[arg(short = 'g', long = "genotypes", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    genotypes: Vec<String>,

    /// Deprecated (an error).
    #[arg(short = 'G', long = "GTs-only", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    gts_only: Vec<String>,

    /// Homozygous genotypes only, useful with low coverage data (requires -g).
    #[arg(short = 'H', long = "homs-only", action = clap::ArgAction::Count)]
    homs_only: u8,

    /// Include sites for which the expression is true ([qry|gt]:EXPR).
    #[arg(short = 'i', long = "include", value_name = "EXPR", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    include: Vec<String>,

    /// Keep monoallelic sites with no alternate allele.
    #[arg(long = "keep-refs", action = clap::ArgAction::Count)]
    keep_refs: u8,

    /// Print only top INT matches for each sample, 0 for unlimited; negative to sort by HWE.
    #[arg(long = "n-matches", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    n_matches: Vec<String>,

    /// Disable calculation of HWE probability.
    #[arg(long = "no-HWE-prob", action = clap::ArgAction::Count)]
    no_hwe_prob: u8,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output: Vec<String>,

    /// t: plain tab-delimited text output, z: compressed [t].
    #[arg(short = 'O', long = "output-type", value_name = "t|z", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output_type: Vec<String>,

    /// Comma-separated sample pairs to compare (`--plot` is an old alias).
    #[arg(short = 'p', long = "pairs", alias = "plot", value_name = "LIST", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    pairs: Vec<String>,

    /// File with tab-delimited sample pairs to compare.
    #[arg(short = 'P', long = "pairs-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    pairs_file: Vec<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(short = 'r', long = "regions", value_name = "REGION", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions: Vec<String>,

    /// Restrict to regions listed in a file.
    #[arg(short = 'R', long = "regions-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions_file: Vec<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(long = "regions-overlap", value_name = "0|1|2", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions_overlap: Vec<String>,

    /// List of query or -g samples ([qry|gt]:LIST), "-" for all.
    #[arg(short = 's', long = "samples", value_name = "LIST", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    samples: Vec<String>,

    /// File with the query or -g samples to compare ([qry|gt]:FILE).
    #[arg(short = 'S', long = "samples-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    samples_file: Vec<String>,

    /// Deprecated (an error).
    #[arg(long = "target-sample", value_name = "STR", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    target_sample: Vec<String>,

    /// Similar to -r but streams rather than index-jumps.
    #[arg(short = 't', long = "targets", value_name = "REGION", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    targets: Vec<String>,

    /// Similar to -R but streams rather than index-jumps.
    #[arg(short = 'T', long = "targets-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    targets_file: Vec<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0].
    #[arg(long = "targets-overlap", value_name = "0|1|2", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    targets_overlap: Vec<String>,

    /// Worker threads for the pair sums and decompression (a vtools addition) [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Which tag to use in the query file (TAG1) and the -g file (TAG2) [PL,GT].
    #[arg(short = 'u', long = "use", value_name = "TAG1[,TAG2]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    use_tags: Vec<String>,

    /// Verbosity level.
    #[arg(short = 'v', long = "verbosity", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    verbosity: Vec<String>,

    /// Print help.
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Count)]
    help: u8,
}

/// Which sites a filter expression keeps, as bcftools combines `-i` and `-e` bits
/// (`FLT_INCLUDE` 1, `FLT_EXCLUDE` 2): only a pure `-e` excludes.
pub(super) const FLT_INCLUDE: u8 = 1;
pub(super) const FLT_EXCLUDE: u8 = 2;

/// The settings of bcftools' `args_t` that its option loop fills.
#[derive(Debug)]
pub(super) struct Settings {
    pub(super) output_fname: String,
    pub(super) gzip: bool,
    pub(super) qry_filter: Option<String>,
    pub(super) gt_filter: Option<String>,
    pub(super) qry_filter_logic: u8,
    pub(super) gt_filter_logic: u8,
    pub(super) gt_err: i32,
    pub(super) qry_use: Option<Tag>,
    pub(super) gt_use: Option<Tag>,
    pub(super) ntop: usize,
    pub(super) sort_by_hwe: bool,
    pub(super) calc_hwe: bool,
    pub(super) dry_run: bool,
    pub(super) distinctive: f64,
    pub(super) es_max_mem: String,
    pub(super) es_tmp_prefix: String,
    pub(super) hom_only: bool,
    pub(super) keep_ref: bool,
    pub(super) gt_fname: Option<String>,
    pub(super) gt_samples: Option<(String, ListSource)>,
    pub(super) qry_samples: Option<(String, ListSource)>,
    pub(super) pair_samples: Option<String>,
    /// A file once `-P` named one, even when a later `-p` gives the list (as in bcftools);
    /// likewise for the regions and targets.
    pub(super) pair_samples_source: ListSource,
    pub(super) regions: Option<String>,
    pub(super) regions_source: ListSource,
    pub(super) targets: Option<String>,
    pub(super) targets_source: ListSource,
    pub(super) regions_overlap: Overlap,
    pub(super) targets_overlap: Overlap,
    pub(super) threads: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            output_fname: "-".to_owned(),
            gzip: false,
            qry_filter: None,
            gt_filter: None,
            qry_filter_logic: 0,
            gt_filter_logic: 0,
            gt_err: 40,
            qry_use: None,
            gt_use: None,
            ntop: 0,
            sort_by_hwe: false,
            calc_hwe: true,
            dry_run: false,
            distinctive: 0.0,
            es_max_mem: "500M".to_owned(),
            es_tmp_prefix: "/tmp/bcftools-gtcheck".to_owned(),
            hom_only: false,
            keep_ref: false,
            gt_fname: None,
            gt_samples: None,
            qry_samples: None,
            pair_samples: None,
            pair_samples_source: ListSource::Inline,
            regions: None,
            regions_source: ListSource::Inline,
            targets: None,
            targets_source: ListSource::Inline,
            regions_overlap: Overlap::Record,
            targets_overlap: Overlap::Pos,
            threads: 0,
        }
    }
}

/// An option value with its command-line position.
type Indexed<'a> = (usize, &'a str, &'a str);

/// Collects the options with values (`id`, value) in command-line order.
fn ordered<'a>(matches: &ArgMatches, args: &'a Args) -> Vec<Indexed<'a>> {
    let lists: [(&'static str, &'a [String]); 25] = [
        ("cluster", &args.cluster),
        ("distinctive_sites", &args.distinctive_sites),
        ("error_probability", &args.error_probability),
        ("exclude", &args.exclude),
        ("genotypes", &args.genotypes),
        ("gts_only", &args.gts_only),
        ("include", &args.include),
        ("n_matches", &args.n_matches),
        ("output", &args.output),
        ("output_type", &args.output_type),
        ("pairs", &args.pairs),
        ("pairs_file", &args.pairs_file),
        ("regions", &args.regions),
        ("regions_file", &args.regions_file),
        ("regions_overlap", &args.regions_overlap),
        ("samples", &args.samples),
        ("samples_file", &args.samples_file),
        ("target_sample", &args.target_sample),
        ("targets", &args.targets),
        ("targets_file", &args.targets_file),
        ("targets_overlap", &args.targets_overlap),
        ("use_tags", &args.use_tags),
        ("verbosity", &args.verbosity),
        ("all_sites", &[]),
        ("dry_run", &[]),
    ];
    let given = |id: &str| matches.value_source(id) == Some(ValueSource::CommandLine);
    let mut out = Vec::new();
    for (id, values) in lists {
        let Some(idx) = matches.indices_of(id).filter(|_| given(id)) else {
            continue;
        };
        if values.is_empty() {
            out.extend(idx.map(|i| (i, id, "")));
        } else {
            out.extend(idx.zip(values).map(|(i, v)| (i, id, v.as_str())));
        }
    }
    for id in ["homs_only", "keep_refs", "no_hwe_prob", "help"] {
        if let Some(idx) = matches.indices_of(id).filter(|_| given(id)) {
            out.extend(idx.map(|i| (i, id, "")));
        }
    }
    out.sort_by_key(|(i, ..)| *i);
    out
}

/// Whether the value is a whole C `strtol` number: the value when it is.
pub(super) fn whole_strtol(s: &str) -> Option<i64> {
    let (v, end) = strtol(s.as_bytes(), 0);
    (end == s.len()).then_some(v)
}

/// Strips a case-insensitive prefix such as `gt:`.
fn strip_prefix_ci<'a>(s: &'a str, prefix: &str) -> Option<&'a str> {
    let n = prefix.len();
    (s.len() >= n && s.as_bytes()[..n].eq_ignore_ascii_case(prefix.as_bytes())).then(|| &s[n..])
}

/// Parses `-O` (bcftools' option loop; the compression level is checked but unused).
fn parse_output_type(s: &str, st: &mut Settings) -> anyhow::Result<()> {
    let b = s.as_bytes();
    match b.first() {
        // bcftools reads the level from past the end of an empty value (README).
        None => bail!("The output type \"\" not recognised"),
        Some(b't') => st.gzip = false,
        Some(b'z') => st.gzip = true,
        _ => match whole_strtol(s) {
            Some(v) if (0..=9).contains(&v) => {}
            _ => bail!("The output type \"{s}\" not recognised"),
        },
    }
    if b.len() > 1 {
        match whole_strtol(&s[1..]) {
            Some(v) if (0..=9).contains(&v) => {}
            _ => bail!("Could not parse argument: --output-type {}", &s[1..]),
        }
    }
    Ok(())
}

/// Parses `-u TAG1[,TAG2]`.
fn parse_use(s: &str, st: &mut Settings) -> anyhow::Result<()> {
    let list: Vec<&str> = s.split(',').collect();
    if list.len() > 2 {
        bail!("Failed to parse --use {s}");
    }
    let tag = |t: &str| {
        if t.eq_ignore_ascii_case("GT") {
            Ok(Tag::Gt)
        } else if t.eq_ignore_ascii_case("PL") {
            Ok(Tag::Pl)
        } else {
            Err(anyhow::anyhow!(
                "Failed to parse --use {s}; only GT and PL are supported"
            ))
        }
    };
    let qry = tag(list[0])?;
    st.qry_use = Some(qry);
    st.gt_use = Some(match list.get(1) {
        Some(t) => tag(t)?,
        None => qry,
    });
    Ok(())
}

/// Parses `--distinctive-sites NUM[,MEM[,TMP]]`.
fn parse_distinctive(s: &str, st: &mut Settings) -> anyhow::Result<()> {
    let (v, end) = strtod(s.as_bytes(), 0);
    st.distinctive = v;
    if end < s.len() {
        if s.as_bytes()[end] != b',' {
            bail!("Could not parse: --distinctive-sites {s}");
        }
        // bcftools keeps the rest, ",TMP" included, as the memory string, so it reads MEM
        // without its unit when TMP is given; the unit is read here.
        let rest = &s[end + 1..];
        match rest.split_once(',') {
            Some((mem, tmp)) => {
                st.es_max_mem = mem.to_owned();
                st.es_tmp_prefix = tmp.to_owned();
            }
            None => st.es_max_mem = rest.to_owned(),
        }
    }
    st.gt_err = 0;
    Ok(())
}

/// Sets the filter of `-i`/`-e` (`logic` is [`FLT_INCLUDE`] or [`FLT_EXCLUDE`]).
fn set_filter(s: &str, logic: u8, opt: char, st: &mut Settings) -> anyhow::Result<()> {
    if let Some(expr) = strip_prefix_ci(s, "gt:") {
        if st.gt_filter.is_some() {
            bail!(
                "Error: only one -i or -e expression can be given for gt:, and they cannot be combined"
            );
        }
        st.gt_filter = Some(expr.to_owned());
        st.gt_filter_logic |= logic;
    } else if let Some(expr) = strip_prefix_ci(s, "qry:") {
        if st.qry_filter.is_some() {
            bail!(
                "Error: only one -i or -e expression can be given for qry:, and they cannot be combined"
            );
        }
        st.qry_filter = Some(expr.to_owned());
        st.qry_filter_logic |= logic;
    } else if opt == 'e'
        && let Some(v) = whole_strtol(s)
    {
        // The old -e, --error-probability.
        st.gt_err = v as i32;
        eprintln!(
            "[warning] auto-detected the old format --error-probability option, please switch from -e to -E."
        );
    } else {
        st.gt_filter = Some(s.to_owned());
        st.qry_filter = Some(s.to_owned());
        st.gt_filter_logic |= logic;
        st.qry_filter_logic |= logic;
    }
    Ok(())
}

/// Sets `-s`/`-S` (a `gt:` or `qry:` list); a file stays a file once `-S` named one, as in
/// bcftools.
fn set_samples(s: &str, source: ListSource, st: &mut Settings) -> anyhow::Result<()> {
    let keep = |old: &Option<(String, ListSource)>| match (old, source) {
        (Some((_, ListSource::File)), _) | (_, ListSource::File) => ListSource::File,
        _ => ListSource::Inline,
    };
    if let Some(list) = strip_prefix_ci(s, "gt:") {
        st.gt_samples = Some((list.to_owned(), keep(&st.gt_samples)));
    } else if let Some(list) = strip_prefix_ci(s, "qry:") {
        st.qry_samples = Some((list.to_owned(), keep(&st.qry_samples)));
    } else {
        bail!("Which one? Query samples (qry:{s}) or genotype samples (gt:{s})?");
    }
    Ok(())
}

/// bcftools' option loop, in command-line order; `None` when the usage is to be printed.
pub(super) fn parse_options(matches: &ArgMatches, args: &Args) -> anyhow::Result<Option<Settings>> {
    let mut st = Settings::default();
    for (_, id, v) in ordered(matches, args) {
        match id {
            "verbosity" => apply_verbosity(v)?,
            "output" => st.output_fname = v.to_owned(),
            "output_type" => parse_output_type(v, &mut st)?,
            "exclude" => set_filter(v, FLT_EXCLUDE, 'e', &mut st)?,
            "include" => set_filter(v, FLT_INCLUDE, 'i', &mut st)?,
            "error_probability" => match whole_strtol(v) {
                Some(n) => st.gt_err = n as i32,
                None => bail!("Could not parse: --error-probability {v}"),
            },
            "use_tags" => parse_use(v, &mut st)?,
            "n_matches" => {
                let Some(n) = whole_strtol(v) else {
                    bail!("Could not parse: --n-matches {v}");
                };
                let n = n as i32;
                st.sort_by_hwe = n < 0 || st.sort_by_hwe;
                st.ntop = n.unsigned_abs() as usize;
            }
            "no_hwe_prob" => st.calc_hwe = false,
            "target_sample" => bail!("The option -S, --target-sample has been deprecated"),
            "dry_run" => st.dry_run = true,
            "distinctive_sites" => parse_distinctive(v, &mut st)?,
            "cluster" => {
                bail!("The -c option is to be implemented, please open an issue on github")
            }
            "gts_only" => bail!("The option -G, --GTs-only has been deprecated"),
            "all_sites" => {
                bail!("The -a option is to be implemented, please open an issue on github")
            }
            "homs_only" => st.hom_only = true,
            "genotypes" => st.gt_fname = Some(v.to_owned()),
            "samples" => set_samples(v, ListSource::Inline, &mut st)?,
            "samples_file" => set_samples(v, ListSource::File, &mut st)?,
            "pairs" => st.pair_samples = Some(v.to_owned()),
            "pairs_file" => {
                st.pair_samples = Some(v.to_owned());
                st.pair_samples_source = ListSource::File;
            }
            "regions" => st.regions = Some(v.to_owned()),
            "regions_file" => {
                st.regions = Some(v.to_owned());
                st.regions_source = ListSource::File;
            }
            "targets" => st.targets = Some(v.to_owned()),
            "targets_file" => {
                st.targets = Some(v.to_owned());
                st.targets_source = ListSource::File;
            }
            "regions_overlap" => {
                st.regions_overlap = Overlap::parse(v)
                    .ok_or_else(|| anyhow::anyhow!("Could not parse: --regions-overlap {v}"))?;
            }
            "keep_refs" => st.keep_ref = true,
            "targets_overlap" => {
                st.targets_overlap = Overlap::parse(v)
                    .ok_or_else(|| anyhow::anyhow!("Could not parse: --targets-overlap {v}"))?;
            }
            "help" => return Ok(None),
            _ => unreachable!("option {id} is listed in `ordered`"),
        }
    }
    if let Some(t) = &args.threads {
        let n = crate::output::cap_threads(whole_strtol(t).unwrap_or(0));
        st.threads = usize::try_from(n).unwrap_or(0);
    }
    Ok(Some(st))
}
