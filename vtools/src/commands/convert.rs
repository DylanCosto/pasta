//! `vtools convert`: convert VCF/BCF to other formats and back (bcftools convert): gVCF block
//! expansion, IMPUTE2 gen/sample and hap/legend/sample, SHAPEIT hap/sample, and TSV to VCF.

mod impute;
mod text;
mod tsv;

use std::{
    io::{self, IsTerminal},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::bail;
use clap::ArgMatches;
use pasta::variant::{
    Header, Record,
    io::{WriteError, Writer},
    readlist::{self, ListSource},
    synced::{Overlap, Regions, SyncedOptions, SyncedReader},
};

use self::text::Reference;
use crate::options::later;
use crate::{
    commands::head::apply_verbosity,
    convert::Convert,
    filter::{self, FilterLogic, RecordFilter},
    output::{self, OutputSetting},
    raw_input::strerror,
    version,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin); further arguments are ignored, as in
    /// bcftools.
    inputs: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Select sites for which the expression is true.
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

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

    /// List of samples to include (sample names with --tsv2vcf).
    #[arg(
        short = 's',
        long = "samples",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to include (sample names with --tsv2vcf).
    #[arg(
        short = 'S',
        long = "samples-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    samples_file: Option<String>,

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

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Output file name [stdout].
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Vec<String>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "u|b|v|z[0-9]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output_type: Vec<String>,

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

    /// Convert IMPUTE2 gen/sample to VCF: <PREFIX>|<GEN-FILE>,<SAMPLE-FILE>.
    #[arg(
        short = 'G',
        long = "gensample2vcf",
        value_name = "FILES",
        allow_hyphen_values = true
    )]
    gensample2vcf: Vec<String>,

    /// Convert to IMPUTE2 gen/sample: <PREFIX>|<GEN-FILE>,<SAMPLE-FILE>.
    #[arg(
        short = 'g',
        long = "gensample",
        value_name = "FILES",
        allow_hyphen_values = true
    )]
    gensample: Vec<String>,

    /// Use 3*N+6 column format instead of the old 3*N+5 column format.
    #[arg(long = "3N6")]
    gen_3n6: bool,

    /// Tag to take values for .gen file: GT,PL,GP [GT].
    #[arg(long = "tag", value_name = "STRING", allow_hyphen_values = true)]
    tag: Option<String>,

    /// Deprecated, use --3N6.
    #[arg(long = "chrom")]
    chrom: bool,

    /// Keep duplicate positions.
    #[arg(long = "keep-duplicates")]
    keep_duplicates: bool,

    /// Output sex column in the sample-file, input format is: Sample\t[MF].
    #[arg(long = "sex", value_name = "FILE", allow_hyphen_values = true)]
    sex: Option<String>,

    /// Output VCF IDs instead of CHROM:POS_REF_ALT.
    #[arg(long = "vcf-ids")]
    vcf_ids: bool,

    /// Expand gVCF reference blocks.
    #[arg(long = "gvcf2vcf")]
    gvcf2vcf: bool,

    /// Reference sequence in fasta format.
    #[arg(
        short = 'f',
        long = "fasta-ref",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    fasta_ref: Option<String>,

    /// Convert SHAPEIT hap/sample to VCF: <PREFIX>|<HAP-FILE>,<SAMPLE-FILE>.
    #[arg(
        long = "hapsample2vcf",
        value_name = "FILES",
        allow_hyphen_values = true
    )]
    hapsample2vcf: Vec<String>,

    /// Convert to SHAPEIT hap/sample: <PREFIX>|<HAP-FILE>,<SAMPLE-FILE>.
    #[arg(long = "hapsample", value_name = "FILES", allow_hyphen_values = true)]
    hapsample: Vec<String>,

    /// Convert haploid genotypes to diploid homozygotes.
    #[arg(long = "haploid2diploid")]
    haploid2diploid: bool,

    /// Convert IMPUTE2 hap/legend/sample to VCF: <PREFIX>|<HAP-FILE>,<LEGEND-FILE>,<SAMPLE-FILE>.
    #[arg(
        short = 'H',
        long = "haplegendsample2vcf",
        value_name = "FILES",
        allow_hyphen_values = true
    )]
    haplegendsample2vcf: Vec<String>,

    /// Convert to IMPUTE2 hap/legend/sample: <PREFIX>|<HAP-FILE>,<LEGEND-FILE>,<SAMPLE-FILE>.
    #[arg(
        short = 'h',
        long = "haplegendsample",
        value_name = "FILES",
        allow_hyphen_values = true
    )]
    haplegendsample: Vec<String>,

    /// Convert a TSV file to VCF.
    #[arg(long = "tsv2vcf", value_name = "FILE", allow_hyphen_values = true)]
    tsv2vcf: Vec<String>,

    /// Columns of the input tsv file [ID,CHROM,POS,AA].
    #[arg(
        short = 'c',
        long = "columns",
        value_name = "STRING",
        allow_hyphen_values = true
    )]
    columns: Option<String>,

    /// Print help.
    #[arg(short = '?')]
    help: bool,

    /// Not an option of bcftools convert: reported as unrecognized, with the usage.
    #[arg(long = "help", hide = true)]
    long_help: bool,
}

const USAGE: &str = "
About:   Converts VCF/BCF to other formats and back. See man page for file
         formats details. When specifying output files explicitly instead
         of with PREFIX, one can use '-' for stdout and '.' to suppress.
Usage:   vtools convert [OPTIONS] INPUT_FILE

VCF input options:
   -e, --exclude EXPR             Exclude sites for which the expression is true
   -i, --include EXPR             Select sites for which the expression is true
   -r, --regions REGION           Restrict to comma-separated list of regions
   -R, --regions-file FILE        Restrict to regions listed in a file
       --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -s, --samples LIST             List of samples to include
   -S, --samples-file FILE        File of samples to include
   -t, --targets REGION           Similar to -r but streams rather than index-jumps
   -T, --targets-file FILE        Similar to -R but streams rather than index-jumps
       --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]

General options:
       --no-version               Do not append version and command line to the header
   -o, --output FILE              Output file name [stdout]
   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
       --threads INT              Use multithreading with INT worker threads [0]
   -v, --verbosity INT            Verbosity level
   -W, --write-index[=FMT]        Automatically index the output files [off]

GEN/SAMPLE conversion (input/output from IMPUTE2):
   -G, --gensample2vcf ...        <PREFIX>|<GEN-FILE>,<SAMPLE-FILE>
   -g, --gensample ...            <PREFIX>|<GEN-FILE>,<SAMPLE-FILE>
       --3N6                      Use 3*N+6 column format instead of the old 3*N+5 column format
       --tag STRING               Tag to take values for .gen file: GT,PL,GL,GP [GT]
       --chrom                    Output chromosome in first column instead of CHROM:POS_REF_ALT
       --keep-duplicates          Keep duplicate positions
       --sex FILE                 Output sex column in the sample-file, input format is: Sample\\t[MF]
       --vcf-ids                  Output VCF IDs in second column instead of CHROM:POS_REF_ALT

gVCF conversion:
       --gvcf2vcf                 Expand gVCF reference blocks
   -f, --fasta-ref FILE           Reference sequence in fasta format

HAP/SAMPLE conversion (output from SHAPEIT):
       --hapsample2vcf ...        <PREFIX>|<HAP-FILE>,<SAMPLE-FILE>
       --hapsample ...            <PREFIX>|<HAP-FILE>,<SAMPLE-FILE>
       --haploid2diploid          Convert haploid genotypes to diploid homozygotes
       --sex FILE                 Output sex column in the sample-file, input format is: Sample\\t[MF]
       --vcf-ids                  Output VCF IDs instead of CHROM:POS_REF_ALT

HAP/LEGEND/SAMPLE conversion:
   -H, --haplegendsample2vcf ...  <PREFIX>|<HAP-FILE>,<LEGEND-FILE>,<SAMPLE-FILE>
   -h, --haplegendsample ...      <PREFIX>|<HAP-FILE>,<LEGEND-FILE>,<SAMPLE-FILE>
       --haploid2diploid          Convert haploid genotypes to diploid homozygotes
       --sex FILE                 Output sex column in the sample-file, input format is: Sample\\t[MF]
       --vcf-ids                  Output VCF IDs instead of CHROM:POS_REF_ALT

TSV conversion:
       --tsv2vcf FILE
   -c, --columns STRING           Columns of the input tsv file, see man page for details [ID,CHROM,POS,AA]
   -f, --fasta-ref FILE           Reference sequence in fasta format
   -s, --samples LIST             List of sample names
   -S, --samples-file FILE        File of sample names
";

/// The conversion to run (bcftools' `convert_func`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Mode {
    VcfToVcf,
    GvcfToVcf,
    VcfToGensample,
    GensampleToVcf,
    VcfToHapsample,
    HapsampleToVcf,
    VcfToHaplegendsample,
    HaplegendsampleToVcf,
    TsvToVcf,
}

/// The parsed command line of `convert`.
pub struct Settings {
    mode: Mode,
    infname: String,
    /// The output name: `-o`, or the files of `-g`, `-h`, `--hapsample`.
    outfname: String,
    filter_expr: Option<(String, FilterLogic)>,
    regions: Option<(String, ListSource)>,
    regions_overlap: Overlap,
    targets: Option<(String, ListSource)>,
    targets_overlap: Overlap,
    samples: Option<(String, ListSource)>,
    sex: Option<String>,
    tag: Option<String>,
    columns: Option<String>,
    fasta_ref: Option<String>,
    /// The `-O` occurrences, resolved with the output name when a VCF/BCF output is opened.
    output_settings: Vec<OutputSetting>,
    threads: usize,
    /// A negative `--threads` (C's `int`), for which bcftools fails to start the thread pool
    /// of a VCF/BCF input.
    negative_threads: Option<i32>,
    write_index: Option<output::WriteIndex>,
    record_cmd_line: bool,
    vcf_ids: bool,
    hap2dip: bool,
    gen_3n6: bool,
    keep_duplicates: bool,
}

/// Collects option values of several ids in command-line order.
fn ordered<'a>(
    matches: &ArgMatches,
    ids: &[(&'a str, &'a [String])],
) -> Vec<(usize, &'a str, &'a String)> {
    let mut out = Vec::new();
    for &(id, values) in ids {
        if let Some(idx) = matches.indices_of(id) {
            for (i, v) in idx.zip(values) {
                out.push((i, id, v));
            }
        }
    }
    out.sort_by_key(|(i, _, _)| *i);
    out
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

/// The conversion and the input and output names it names (`None` for the input of a
/// VCF/BCF conversion), in command-line order: the last one wins.
fn conversion(args: &Args, matches: &ArgMatches) -> (Mode, Option<String>, String) {
    let mut mode = Mode::VcfToVcf;
    let mut infname: Option<String> = None;
    let mut outfname = "-".to_owned();
    // A flag, given a value so that it takes part in the ordering.
    let gvcf = if args.gvcf2vcf {
        vec![String::new()]
    } else {
        Vec::new()
    };
    for (_, id, value) in ordered(
        matches,
        &[
            ("gensample", &args.gensample),
            ("gensample2vcf", &args.gensample2vcf),
            ("tsv2vcf", &args.tsv2vcf),
            ("hapsample2vcf", &args.hapsample2vcf),
            ("gvcf2vcf", &gvcf),
            ("hapsample", &args.hapsample),
            ("haplegendsample2vcf", &args.haplegendsample2vcf),
            ("haplegendsample", &args.haplegendsample),
            ("output", &args.output),
        ],
    ) {
        match id {
            "gensample" => (mode, outfname) = (Mode::VcfToGensample, value.clone()),
            "gensample2vcf" => (mode, infname) = (Mode::GensampleToVcf, Some(value.clone())),
            "tsv2vcf" => (mode, infname) = (Mode::TsvToVcf, Some(value.clone())),
            "hapsample2vcf" => (mode, infname) = (Mode::HapsampleToVcf, Some(value.clone())),
            "gvcf2vcf" => mode = Mode::GvcfToVcf,
            "hapsample" => (mode, outfname) = (Mode::VcfToHapsample, value.clone()),
            "haplegendsample2vcf" => {
                (mode, infname) = (Mode::HaplegendsampleToVcf, Some(value.clone()));
            }
            "haplegendsample" => (mode, outfname) = (Mode::VcfToHaplegendsample, value.clone()),
            _ => outfname = value.clone(),
        }
    }
    (mode, infname, outfname)
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    if args.long_help {
        eprintln!("convert: unrecognized option '--help'");
    }
    if args.help || args.long_help {
        eprintln!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    let filter_expr =
        FilterLogic::single(&args.include, &args.exclude)?.map(|(e, logic)| (e.to_owned(), logic));
    if args.chrom {
        bail!("The --chrom option has been deprecated, please use --3N6 instead");
    }

    let (mode, infname, outfname) = conversion(&args, matches);
    let infname = match infname {
        Some(f) => f,
        None => match args.inputs.first() {
            Some(f) => f.clone(),
            None if !io::stdin().is_terminal() => "-".to_owned(),
            None => {
                eprintln!("{USAGE}");
                return Ok(ExitCode::from(1));
            }
        },
    };

    let write_index = match &args.write_index {
        Some(w) => Some(output::parse_write_index(w)?),
        None => None,
    };
    let type_settings: Vec<OutputSetting> = args
        .output_type
        .iter()
        .map(|t| OutputSetting::Type(t.clone()))
        .collect();
    // bcftools' option loop rejects an unknown -O type; the name's extension and the level
    // apply only where a VCF/BCF output is opened (not for -g, -h, --hapsample).
    output::requested_type(&type_settings)?;

    // bcftools' `strtol` into an `int`.
    let threads = args
        .threads
        .as_deref()
        .map_or(0, crate::output::thread_count);
    let settings = Settings {
        mode,
        infname,
        outfname,
        filter_expr,
        regions: later(
            matches,
            ("regions", args.regions.as_deref()),
            ("regions_file", args.regions_file.as_deref()),
        )
        .map(|(spec, source)| (spec.to_owned(), source)),
        regions_overlap: parse_overlap(
            args.regions_overlap.as_deref(),
            Overlap::Record,
            "regions-overlap",
        )?,
        targets: later(
            matches,
            ("targets", args.targets.as_deref()),
            ("targets_file", args.targets_file.as_deref()),
        )
        .map(|(spec, source)| (spec.to_owned(), source)),
        targets_overlap: parse_overlap(
            args.targets_overlap.as_deref(),
            Overlap::Pos,
            "targets-overlap",
        )?,
        samples: later(
            matches,
            ("samples", args.samples.as_deref()),
            ("samples_file", args.samples_file.as_deref()),
        )
        .map(|(spec, source)| (spec.to_owned(), source)),
        sex: args.sex,
        tag: args.tag,
        columns: args.columns,
        fasta_ref: args.fasta_ref,
        output_settings: type_settings,
        threads: usize::try_from(threads).unwrap_or(0),
        negative_threads: (threads < 0).then_some(threads),
        write_index,
        record_cmd_line: !args.no_version,
        vcf_ids: args.vcf_ids,
        hap2dip: args.haploid2diploid,
        gen_3n6: args.gen_3n6,
        keep_duplicates: args.keep_duplicates,
    };
    match settings.mode {
        Mode::VcfToVcf => vcf_to_vcf(&settings)?,
        Mode::GvcfToVcf => gvcf_to_vcf(&settings)?,
        Mode::VcfToGensample => impute::vcf_to_gensample(&settings)?,
        Mode::VcfToHaplegendsample => impute::vcf_to_haplegendsample(&settings)?,
        Mode::VcfToHapsample => impute::vcf_to_hapsample(&settings)?,
        Mode::GensampleToVcf => tsv::gensample_to_vcf(&settings)?,
        Mode::HaplegendsampleToVcf => tsv::haplegendsample_to_vcf(&settings)?,
        Mode::HapsampleToVcf => tsv::hapsample_to_vcf(&settings)?,
        Mode::TsvToVcf => tsv::tsv_to_vcf(&settings)?,
    }
    Ok(ExitCode::SUCCESS)
}

/// An opened VCF/BCF input with its filter and format (bcftools' `open_vcf`).
struct Input {
    reader: SyncedReader,
    filter: Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    convert: Option<Convert>,
}

impl Input {
    /// Applies `-i/-e`; `true` when the record is kept.
    fn pass(&mut self, rec: &mut Record) -> anyhow::Result<bool> {
        let Some((flt, logic)) = &mut self.filter else {
            return Ok(true);
        };
        let pass = flt.pass(self.reader.header(), rec)?;
        Ok(match logic {
            FilterLogic::Include => pass,
            FilterLogic::Exclude => !pass,
        })
    }
}

fn open_vcf(s: &Settings, format_str: Option<&str>) -> anyhow::Result<Input> {
    if let Some(n) = s.negative_threads {
        bail!("Could not initialize --threads {n}");
    }
    if !crate::output::threads_start(s.threads) {
        bail!("Could not initialize --threads {}", s.threads);
    }
    let mut regions = None;
    if let Some((spec, source)) = &s.regions {
        let mut r = Regions::new(spec.as_bytes(), *source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = s.regions_overlap;
        regions = Some(r);
    }
    let mut targets = None;
    let mut targets_exclude = false;
    if let Some((spec, source)) = &s.targets {
        let body = match spec.strip_prefix('^') {
            Some(rest) => {
                targets_exclude = true;
                rest
            }
            None => spec.as_str(),
        };
        let mut t = Regions::new(body.as_bytes(), *source)
            .map_err(|_| anyhow::anyhow!("Failed to read the targets: {spec}"))?;
        t.overlap = s.targets_overlap;
        targets = Some(t);
    }
    let mut reader = SyncedReader::open(
        &s.infname,
        SyncedOptions {
            regions,
            targets,
            targets_exclude,
            apply_filters: None,
            threads: s.threads,
        },
    )
    .map_err(|e| anyhow::anyhow!("Failed to open {}: {e}", s.infname))?;

    let filter = match &s.filter_expr {
        Some((expr, logic)) => Some((filter::init(reader.header(), expr)?, *logic)),
        None => None,
    };

    let mut order: Option<Vec<usize>> = None;
    if let Some((list, source)) = &s.samples
        && list != "-"
    {
        let ret = reader
            .header_mut()
            .set_samples(Some(list.as_bytes()), *source)
            .map_err(|_| anyhow::anyhow!("Error parsing the sample list"))?;
        if ret > 0 {
            bail!("Sample name mismatch: sample #{ret} not found in the header");
        }
        if !list.starts_with('^') {
            // The sample order may differ from the header's when not negated.
            let Ok(names) = readlist::read_list(list.as_bytes(), *source) else {
                bail!("Could not parse {list}");
            };
            let header = reader.header();
            if names.len() != header.sample_count() {
                bail!(
                    "The number of samples does not match, perhaps some are present multiple times?"
                );
            }
            order = Some(
                names
                    .iter()
                    .map(|n| header.sample_id(n).map_or(0, |i| i as usize))
                    .collect(),
            );
        }
    }
    let convert = match format_str {
        Some(f) => Some(Convert::new(reader.header(), order.as_deref(), f)?),
        None => None,
    };
    Ok(Input {
        reader,
        filter,
        convert,
    })
}

/// The output name as bcftools prints it in messages (`-` for standard output).
fn out_display(s: &Settings) -> &str {
    &s.outfname
}

/// Opens the VCF/BCF output and writes the header (`set_wmode`, `hts_open`, `bcf_hdr_write`,
/// and `init_index2`).
fn open_output(s: &Settings, header: &Header, func: &str) -> anyhow::Result<Writer> {
    // bcftools' set_wmode, which names standard output `-` in its messages.
    let (output_type, level) = output::resolve(&s.output_settings, Some(Path::new(&s.outfname)))?;
    let path = (s.outfname != "-").then(|| PathBuf::from(&s.outfname));
    let data = path.as_deref().map(output::data_path);
    let mut writer = Writer::create(
        data.as_deref(),
        &output::write_options(output_type, level, s.threads),
    )
    .map_err(|e| anyhow::anyhow!("Can't write to \"{}\": {}", out_display(s), strerror(&e)))?;
    if let Some(fmt) = s.write_index {
        let ok = match output::index_request(fmt, output_type, path.as_deref()) {
            Some(req) => writer.request_index(req).is_ok(),
            None => false,
        };
        if !ok {
            bail!("Error: failed to initialise index for {}", out_display(s));
        }
    }
    if writer.write_header(header).is_err() {
        bail!("[{func}] Error: cannot write to {}", out_display(s));
    }
    Ok(writer)
}

/// Writes a record, mapping failures to bcftools' message.
fn write_record(
    writer: &mut Writer,
    header: &Header,
    rec: &mut Record,
    s: &Settings,
    func: &str,
) -> anyhow::Result<()> {
    match writer.write_record(header, rec) {
        Ok(()) => Ok(()),
        Err(WriteError::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => {
            Err(anyhow::Error::new(e))
        }
        Err(_) => bail!("[{func}] Error: cannot write to {}", out_display(s)),
    }
}

/// How a converter words a failed close of its output.
#[derive(Clone, Copy, Debug)]
enum CloseMessage {
    /// `Close failed: FILE` (the converters from IMPUTE2 files).
    Short,
    /// `[FUNC] Error: close failed .. FILE`.
    Named,
}

/// Closes the output, saving the index; `func` names the converter in messages.
fn finish_output(
    writer: Writer,
    s: &Settings,
    func: &str,
    message: CloseMessage,
) -> anyhow::Result<()> {
    let out = out_display(s);
    let close_failed = || match message {
        CloseMessage::Short => anyhow::anyhow!("Close failed: {out}"),
        CloseMessage::Named => anyhow::anyhow!("[{func}] Error: close failed .. {out}"),
    };
    match writer.finish() {
        Ok(()) => Ok(()),
        Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
        Err(WriteError::Deferred(e)) => Err(output::write_failure(e, || {
            anyhow::anyhow!("[{func}] Error: cannot write to {out}")
        })),
        Err(WriteError::Io(e)) => Err(output::write_failure(e, close_failed)),
        Err(e) => Err(anyhow::Error::new(e).context(close_failed())),
    }
}

/// Plain VCF/BCF conversion (bcftools' `vcf_to_vcf`); the header is written unchanged.
fn vcf_to_vcf(s: &Settings) -> anyhow::Result<()> {
    let mut input = open_vcf(s, None)?;
    let mut writer = open_output(s, input.reader.header(), "vcf_to_vcf")?;
    let result = (|| -> anyhow::Result<()> {
        let mut rec = Record::new();
        loop {
            match input.reader.next_record(&mut rec) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => bail!("Error: {e}"),
            }
            if !input.pass(&mut rec)? {
                continue;
            }
            write_record(
                &mut writer,
                input.reader.header(),
                &mut rec,
                s,
                "vcf_to_vcf",
            )?;
        }
        if let Some(e) = input.reader.take_pending_error() {
            bail!("Error: {e}");
        }
        Ok(())
    })();
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    finish_output(writer, s, "vcf_to_vcf", CloseMessage::Named)
}

/// Whether ALT marks a gVCF reference block: `<*>`, `<X>`, or `<NON_REF>`.
fn is_gvcf_record(rec: &Record) -> bool {
    match rec.allele_count() {
        1 => true,
        n => {
            rec.allele(1).first() == Some(&b'<')
                && (1..n).any(|i| matches!(rec.allele(i), b"<*>" | b"<X>" | b"<NON_REF>"))
        }
    }
}

/// Expands gVCF reference blocks into one record per position (bcftools' `gvcf_to_vcf`).
fn gvcf_to_vcf(s: &Settings) -> anyhow::Result<()> {
    let Some(ref_fname) = &s.fasta_ref else {
        bail!("--gvcf2vcf requires the --fasta-ref option");
    };
    let mut reference = Reference::open(Path::new(ref_fname))
        .map_err(|_| anyhow::anyhow!("Could not load the fai index for reference {ref_fname}"))?;
    let mut input = open_vcf(s, None)?;
    if s.record_cmd_line {
        version::append_version(input.reader.header_mut(), "convert")?;
    }
    let mut writer = open_output(s, input.reader.header(), "gvcf_to_vcf")?;
    let result = expand_gvcf(s, &mut input, &mut writer, &mut reference);
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    finish_output(writer, s, "gvcf_to_vcf", CloseMessage::Named)
}

fn expand_gvcf(
    s: &Settings,
    input: &mut Input,
    writer: &mut Writer,
    reference: &mut Reference,
) -> anyhow::Result<()> {
    const FUNC: &str = "gvcf_to_vcf";
    // A one-record lookahead (bcftools' vcfbuf of two sites).
    let mut cur = Record::new();
    let mut next = Record::new();
    let mut have_next = match input.reader.next_record(&mut next) {
        Ok(b) => b,
        Err(e) => bail!("Error: {e}"),
    };
    let mut warned = false;
    let mut ends: Vec<i32> = Vec::new();
    let mut alleles: Vec<Vec<u8>> = Vec::new();
    let mut base = Vec::new();
    while have_next {
        std::mem::swap(&mut cur, &mut next);
        have_next = match input.reader.next_record(&mut next) {
            Ok(b) => b,
            Err(e) => bail!("Error: {e}"),
        };
        let header = input.reader.header();
        cur.unpack(pasta::variant::Unpack::STR);

        // Is it a gVCF block: ALT one of ., <*>, <X>, <NON_REF>, and INFO/END present?
        let mut end1 = 0i64;
        if is_gvcf_record(&cur) && cur.get_info_i32(header, b"END", &mut ends).ok() == Some(1) {
            end1 = i64::from(ends[0]);
            if have_next && cur.rid() == next.rid() && end1 > next.pos() {
                if !warned {
                    pasta::variant::log::warn(
                        "gvcf_next_line",
                        format_args!(
                            "Malformed gVCF: INFO/END at {}:{} overlaps the next record",
                            String::from_utf8_lossy(cur.chrom_safe(header)),
                            cur.pos() + 1
                        ),
                    );
                    warned = true;
                }
                end1 = if cur.pos() < next.pos() {
                    next.pos()
                } else {
                    0
                };
            }
        }

        if !input.pass(&mut cur)? || end1 == 0 {
            write_record(writer, input.reader.header(), &mut cur, s, FUNC)?;
            continue;
        }

        let header = input.reader.header();
        cur.remove_info(header, b"END")
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        alleles.clear();
        alleles.extend(cur.alleles().map(<[u8]>::to_vec));
        let chrom = cur.chrom_safe(header).to_vec();
        for pos in cur.pos()..end1 {
            cur.set_pos(pos);
            if !reference.fetch(&chrom, pos, &mut base)? {
                bail!(
                    "faidx_fetch_seq failed at {}:{}",
                    String::from_utf8_lossy(&chrom),
                    pos + 1
                );
            }
            // strncpy of the fetched base over the start of REF.
            for (dst, &b) in alleles[0].iter_mut().zip(&base) {
                *dst = b;
            }
            cur.update_alleles(header, &alleles);
            write_record(writer, header, &mut cur, s, FUNC)?;
        }
    }
    if let Some(e) = input.reader.take_pending_error() {
        bail!("Error: {e}");
    }
    Ok(())
}

/// Splits `-g`/`-h`/`--hapsample` output names on commas (bcftools' `hts_readlist(..., 0)`).
fn split_names(names: &str) -> Vec<String> {
    names.split(',').map(str::to_owned).collect()
}
