//! `vtools csq`: haplotype-aware consequence calling (bcftools csq, `csq.c`).
//!
//! The structure follows bcftools. Each record is buffered and looked up in four position
//! indexes of the GFF3 gene models: coding sequences, UTRs, exons (for splice sites), and
//! transcripts. Alleles in coding sequences extend a haplotype tree per transcript (one leaf
//! per sample haplotype, or one haplotype without samples); when a transcript ends, each
//! haplotype is spliced, translated, and broken into independent consequences, which go to
//! the records they come from. Records are written once no active transcript needs them,
//! with INFO/BCSQ and a FORMAT/BCSQ bitmask of the consequences each haplotype carries.

mod annotate;
mod gencode;
mod hap;
mod splice;
mod translate;
mod types;
mod vbuf;

use std::{
    collections::{HashMap, VecDeque},
    fs::File,
    io::{self, BufWriter, IsTerminal, Write},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::{ArgMatches, builder::TypedValueParser as _};
use pasta::variant::{
    Header, Record, Unpack,
    io::{WriteError, Writer},
    log,
    readlist::{self, ListSource},
    synced::{Overlap, Regions, SyncedOptions, SyncedReader},
};

use self::{
    gencode::{Gencode, GencodeChoice},
    hap::{Hap, TrAux},
    types::RecRef,
    vbuf::{TrHeap, VBuf},
};
use crate::options::later_file_sticky;
use crate::{
    commands::norm::reference::Reference,
    filter::{self, FilterLogic, RecordFilter},
    gff::{Gff, GffOptions},
    output::{self, OutputSetting},
    raw_input::strerror,
    version,
};

/// The number of reference bases padding each transcript, to avoid boundary effects
/// (bcftools' `N_REF_PAD`).
const N_REF_PAD: usize = 10;

/// The consequence bits each FORMAT/BCSQ integer holds (bcftools packs 30 per `int32`).
const BCSQ_BITS_PER_INT: i32 = 30;

/// bcftools' `REGIDX_MAX`: flushing up to here flushes everything.
const REGIDX_MAX: u32 = 2_147_483_646;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file (`-` or none for stdin).
    inputs: Vec<String>,

    /// Reference file in fasta format.
    #[arg(short = 'f', long = "fasta-ref", value_name = "FILE", value_parser = clap::builder::OsStringValueParser::new().map(PathBuf::from))]
    fasta_ref: Option<PathBuf>,

    /// GFF3 annotation file.
    #[arg(short = 'g', long = "gff-annot", value_name = "FILE", value_parser = clap::builder::OsStringValueParser::new().map(PathBuf::from))]
    gff_annot: Option<PathBuf>,

    /// Abbreviate protein-changing predictions to max INT aminoacids.
    #[arg(
        short = 'B',
        long = "trim-protein-seq",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    trim_protein_seq: Vec<String>,

    /// Deprecated: the same as `-B 1`.
    #[arg(short = 'b', long = "brief-predictions", action = clap::ArgAction::Count)]
    brief_predictions: u8,

    /// Specify the genetic code table to use, 'l' to print a list [0].
    #[arg(short = 'C', long = "genetic-code", value_name = "INT|l")]
    genetic_code: Option<String>,

    /// Use this tag instead of the default BCSQ.
    #[arg(
        short = 'c',
        long = "custom-tag",
        value_name = "STRING",
        allow_hyphen_values = true
    )]
    custom_tag: Option<String>,

    /// Also check transcript-level consequences after feature hits [0].
    #[arg(
        short = 'G',
        long = "greedy",
        value_name = "0|1",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    greedy: Vec<String>,

    /// Localized predictions, consider only one VCF record at a time.
    #[arg(short = 'l', long = "local-csq")]
    local_csq: bool,

    /// Maximum number of per-haplotype consequences to consider for each site [15].
    #[arg(
        short = 'n',
        long = "ncsq",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    ncsq: Vec<String>,

    /// How to handle unphased heterozygous genotypes: a|m|r|R|s [r].
    #[arg(
        short = 'p',
        long = "phase",
        value_name = "a|m|r|R|s",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    phase: Vec<String>,

    /// Dump the parsed GFF file (for debugging purposes).
    #[arg(long = "dump-gff", value_name = "FILE.gz", value_parser = clap::builder::OsStringValueParser::new().map(PathBuf::from))]
    dump_gff: Option<PathBuf>,

    /// Run even if some sanity checks fail.
    #[arg(long = "force")]
    force: bool,

    /// Unify chromosome naming by stripping a prefix in VCF,GFF,fasta, respectively [0].
    #[arg(
        long = "unify-chr-names",
        value_name = "0|LIST",
        allow_hyphen_values = true
    )]
    unify_chr_names: Option<String>,

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

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_name = "FILE", value_parser = clap::builder::OsStringValueParser::new().map(PathBuf::from))]
    output: Option<PathBuf>,

    /// b|u|z|v: BCF or VCF, compressed or not; t: plain tab-delimited text; 0-9: level [v].
    #[arg(short = 'O', long = "output-type", value_name = "b|u|z|v|t[0-9]", action = clap::ArgAction::Append)]
    output_type: Vec<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(short = 'r', long = "regions", value_name = "REGION")]
    regions: Option<String>,

    /// Restrict to regions listed in a file.
    #[arg(short = 'R', long = "regions-file", value_name = "FILE")]
    regions_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(
        long = "regions-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    regions_overlap: Vec<String>,

    /// Samples to include or "-" to apply all variants and ignore samples.
    #[arg(
        short = 's',
        long = "samples",
        value_name = "-|LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// Samples to include.
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
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
    #[arg(short = 'T', long = "targets-file", value_name = "FILE")]
    targets_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0].
    #[arg(
        long = "targets-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    targets_overlap: Vec<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(
        long = "threads",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    threads: Vec<String>,

    /// Verbosity level 0-6 [1].
    #[arg(
        short = 'v',
        long = "verbosity",
        visible_alias = "verbose",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,

    /// Deprecated, use -v.
    #[arg(short = 'q', long = "quiet", hide = true, action = clap::ArgAction::Count)]
    quiet: u8,

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
    #[arg(short = 'h', long = "help", short_alias = '?', action = clap::ArgAction::Count)]
    help: u8,
}

/// How unphased heterozygous genotypes are handled (`-p`; bcftools' `PHASE_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Phase {
    /// `r`: require phased genotypes.
    Require,
    /// `m`: merge all genotypes into a single haplotype.
    Merge,
    /// `a`: take genotypes as they are.
    AsIs,
    /// `s`: skip unphased hets.
    Skip,
    /// `R`: non-reference haplotypes where possible.
    NonRef,
    /// `-s -` or no samples: apply all variants, ignore samples.
    DropGt,
}

/// The options that drive consequence calling.
#[derive(Debug)]
struct Opts {
    phase: Phase,
    verbosity: i32,
    local_csq: bool,
    greedy: bool,
    /// Twice `--ncsq`: the number of consequence bits per sample.
    ncsq2_max: i32,
    /// The number of FORMAT/BCSQ integers per sample.
    nfmt_bcsq: usize,
    /// `-B`, 0 for off.
    brief_predictions: i32,
    force: bool,
    bcsq_tag: Vec<u8>,
    /// Tab-delimited text output (`-Ot`).
    text: bool,
}

/// Warnings printed once (bcftools' `warned` and static flags).
#[derive(Debug, Default)]
struct Warned {
    ref_allele_mismatch: usize,
    faidx_fetch_failed: usize,
    overlaps: bool,
    multiploid: bool,
    small_ref_padding: bool,
    symbolic: bool,
    missing_chr_fai: bool,
    missing_chr_gff: bool,
}

/// Where the output goes.
enum Output {
    Vcf { writer: Box<Writer>, name: String },
    Text(Box<dyn Write + Send>),
}

/// The record being annotated: position, length, and alleles, copied so that the record
/// itself can sit in the buffer.
#[derive(Debug, Default)]
struct CurRec {
    serial: u64,
    pos: i64,
    rlen: i64,
    /// The alleles as bcftools stores them: `REF\0ALT1\0...`.
    als: Vec<u8>,
    starts: Vec<usize>,
}

impl CurRec {
    fn fill(&mut self, rec: &Record) {
        self.pos = rec.pos();
        self.rlen = rec.rlen();
        self.als.clear();
        self.starts.clear();
        for a in rec.alleles() {
            self.starts.push(self.als.len());
            self.als.extend_from_slice(a);
            self.als.push(0);
        }
    }

    fn n_allele(&self) -> usize {
        self.starts.len()
    }

    /// Allele `i` up to its end.
    fn allele_str(&self, i: usize) -> &[u8] {
        let s = &self.als[self.starts[i]..];
        &s[..memchr::memchr(0, s).unwrap_or(s.len())]
    }

    /// Allele `i` and the alleles after it (a C pointer into `d.als`).
    fn allele_raw(&self, i: usize) -> &[u8] {
        &self.als[self.starts[i]..]
    }

    fn rec_ref(&self) -> RecRef {
        RecRef {
            serial: self.serial,
            pos: self.pos,
        }
    }
}

/// The genotypes of the current record (bcftools' `gt_arr`, which `test_cds` changes and
/// `csq_stage` reads again).
#[derive(Debug, Default)]
struct Genotypes {
    /// As read.
    orig: Vec<i32>,
    /// The count of `get_genotypes` once read (negative without GT).
    ret: Option<i32>,
    /// As `test_cds` changes them.
    arr: Vec<i32>,
    /// Whether `arr` differs from `orig`.
    dirty: bool,
    /// The samples and haplotypes carrying each allele, for `csq_stage`.
    carriers: Vec<Vec<(u32, u32)>>,
    carriers_built: bool,
}

impl Genotypes {
    /// Undoes the changes to `arr` (bcftools reads the genotypes again).
    fn restore(&mut self) {
        if self.dirty {
            self.arr.clone_from(&self.orig);
            self.dirty = false;
        }
    }

    /// Changes a genotype of `arr`.
    fn set(&mut self, i: usize, g: i32) {
        self.arr[i] = g;
        self.dirty = true;
    }
}

/// The source of a chromosome name, for `--unify-chr-names` (bcftools' `CHR_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ChrSource {
    Vcf = 0,
    Gff = 1,
    Fai = 2,
}

/// The consequence caller (bcftools' `args_t`).
struct Csq {
    opts: Opts,
    gencode: &'static Gencode,
    gff: Gff,
    /// The end of each transcript, for the heap of active transcripts.
    tr_ends: Vec<u32>,
    fai: Reference,
    fa_fname: String,
    gff_fname: String,
    hdr: Header,
    nsamples: usize,
    /// The number of samples with FORMAT/BCSQ (0 without genotypes).
    hdr_nsmpl: usize,
    /// The header indexes of the samples to process.
    smpl: Vec<usize>,
    filter: Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    out: Output,
    chr_prefix: [Option<Vec<u8>>; 3],
    unify_err: String,
    // per-chromosome state
    prev_rid: i32,
    prev_pos: i64,
    rid: i32,
    chr_vcf: Vec<u8>,
    chr_gff: Vec<u8>,
    chr_fai: Vec<u8>,
    // annotation state
    active_tr: TrHeap,
    aux: Vec<Option<Box<TrAux>>>,
    hap: Hap,
    vbufs: VecDeque<VBuf>,
    vbuf_pool: Vec<VBuf>,
    next_vbuf_id: u64,
    pos2vbuf: HashMap<i32, u64>,
    rm_tr: Vec<u32>,
    serial: u64,
    cur: CurRec,
    gt: Genotypes,
    ncsq2_small_warned: i32,
    warned: Warned,
    str_buf: Vec<u8>,
    fmt_i32: Vec<i32>,
}

/// C's `atoi`.
fn atoi(s: &str) -> i64 {
    pasta::variant::number::strtol(s.as_bytes(), 0).0
}

/// C's `strtol` that must consume the whole string.
fn strtol_strict(s: &str) -> Option<i64> {
    let (v, end) = pasta::variant::number::strtol(s.as_bytes(), 0);
    (end == s.len()).then_some(v)
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

/// The usage text, printed for `-h` and argument errors.
const USAGE: &str = "
About: Haplotype-aware consequence caller.
Usage: vtools csq [OPTIONS] in.vcf

Required options:
   -f, --fasta-ref FILE              Reference file in fasta format
   -g, --gff-annot FILE              GFF3 annotation file

CSQ options:
   -B, --trim-protein-seq INT        Abbreviate protein-changing predictions to max INT aminoacids
   -C, --genetic-code INT|l          Specify the genetic code table to use, 'l' to print a list [0]
   -c, --custom-tag STRING           Use this tag instead of the default BCSQ
   -G, --greedy 0|1                  Also check transcript-level consequences after feature hits [0]
   -l, --local-csq                   Localized predictions, consider only one VCF record at a time
   -n, --ncsq INT                    Maximum number of per-haplotype consequences to consider for each site [15]
   -p, --phase a|m|r|R|s             How to handle unphased heterozygous genotypes: [r]
                                       a: take GTs as is, create haplotypes regardless of phase (0/1 -> 0|1)
                                       m: merge *all* GTs into a single haplotype (0/1 -> 1, 1/2 -> 1)
                                       r: require phased GTs, throw an error on unphased het GTs
                                       R: create non-reference haplotypes if possible (0/1 -> 1|1, 1/2 -> 1|2)
                                       s: skip unphased hets
GFF options:
       --dump-gff FILE.gz            Dump the parsed GFF file (for debugging purposes)
       --force                       Run even if some sanity checks fail
       --unify-chr-names 0|LIST      Unify chromosome naming by stripping a prefix in VCF,GFF,fasta, respectively [0]
                                     (e.g., \"chr,Chr,-\" trims \"chr\" in VCF and \"Chr\" in GFF, fasta is unchanged)
General options:
   -e, --exclude EXPR                Exclude sites for which the expression is true
   -i, --include EXPR                Select sites for which the expression is true
       --no-version                  Do not append version and command line to the header
   -o, --output FILE                 Write output to a file [standard output]
   -O, --output-type b|u|z|v|t[0-9]  b: compressed BCF, u: uncompressed BCF, z: compressed VCF
                                     v: uncompressed VCF, t: plain tab-delimited text output, 0-9: compression level [v]
   -r, --regions REGION              Restrict to comma-separated list of regions
   -R, --regions-file FILE           Restrict to regions listed in a file
       --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -s, --samples -|LIST              Samples to include or \"-\" to apply all variants and ignore samples
   -S, --samples-file FILE           Samples to include
   -t, --targets REGION              Similar to -r but streams rather than index-jumps
   -T, --targets-file FILE           Similar to -R but streams rather than index-jumps
       --targets-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
       --threads INT                 Use multithreading with <int> worker threads [0]
   -v, --verbosity INT               Verbosity level 0-6 [1]
   -W, --write-index[=FMT]           Automatically index the output files [off]

Example:
   vtools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.87.gff3.gz in.vcf

   # GFF3 annotation files can be downloaded from Ensembl. e.g. for human:
   http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
   http://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/

";

/// The output type of the `-O` occurrences (bcftools' `case 'O'` of `main_csq`).
#[derive(Debug, Default)]
struct OutputTypeOption {
    /// `t`: the tab-delimited text output.
    text: bool,
    ty: Option<char>,
    level: Option<i64>,
}

impl OutputTypeOption {
    /// Applies one `-O` value.
    fn apply(&mut self, v: &str) -> anyhow::Result<()> {
        let mut chars = v.chars();
        match chars.next() {
            Some('t') => {
                self.text = true;
                self.ty = None;
            }
            Some(c @ ('b' | 'u' | 'z' | 'v')) => {
                self.text = false;
                self.ty = Some(c);
            }
            _ => match strtol_strict(v) {
                Some(n) if (0..=9).contains(&n) => self.level = Some(n),
                _ => bail!("The output type \"{v}\" not recognised"),
            },
        }
        if v.len() > 1 {
            let rest = &v[1..];
            match strtol_strict(rest) {
                Some(n) if (0..=9).contains(&n) => self.level = Some(n),
                _ => bail!("Could not parse argument: --output-type {rest}"),
            }
        }
        Ok(())
    }

    /// The VCF/BCF output settings.
    fn settings(&self) -> Vec<OutputSetting> {
        let mut setting = self.ty.unwrap_or('v').to_string();
        if let Some(l) = self.level {
            setting.push_str(&l.to_string());
        }
        vec![OutputSetting::Type(setting)]
    }
}

/// One option occurrence that bcftools checks while parsing the command line.
#[derive(Clone, Copy, Debug)]
enum OptionValue<'a> {
    Help,
    Quiet,
    Threads(&'a str),
    Brief,
    TrimProteinSeq(&'a str),
    Greedy(&'a str),
    Verbosity(&'a str),
    Phase(&'a str),
    Ncsq(&'a str),
    OutputType(&'a str),
    Filter(&'a str, FilterLogic),
    RegionsOverlap(&'a str),
    TargetsOverlap(&'a str),
    WriteIndex(&'a str),
}

/// The option values bcftools checks while parsing the command line.
#[derive(Debug)]
struct ParsedOptions {
    verbosity: i32,
    /// `-b`/`-B`, 0 for off.
    brief: i32,
    greedy: bool,
    phase: Phase,
    ncsq2_max: i32,
    out_type: OutputTypeOption,
    filter: Option<(String, FilterLogic)>,
    threads: i32,
    regions_overlap: Overlap,
    targets_overlap: Overlap,
    write_index: Option<output::WriteIndex>,
}

/// The occurrences of the options bcftools checks as it reads them, in command-line order.
fn ordered_options<'a>(args: &'a Args, matches: &ArgMatches) -> Vec<(usize, OptionValue<'a>)> {
    let mut events: Vec<(usize, OptionValue<'a>)> = Vec::new();
    let mut add = |id: &str, values: &'a [String], make: fn(&'a str) -> OptionValue<'a>| {
        if let Some(idx) = matches.indices_of(id) {
            events.extend(idx.zip(values).map(|(i, v)| (i, make(v))));
        }
    };
    add("threads", &args.threads, OptionValue::Threads);
    add(
        "trim_protein_seq",
        &args.trim_protein_seq,
        OptionValue::TrimProteinSeq,
    );
    add("greedy", &args.greedy, OptionValue::Greedy);
    add("verbosity", &args.verbosity, OptionValue::Verbosity);
    add("phase", &args.phase, OptionValue::Phase);
    add("ncsq", &args.ncsq, OptionValue::Ncsq);
    add("output_type", &args.output_type, OptionValue::OutputType);
    add("include", &args.include, |v| {
        OptionValue::Filter(v, FilterLogic::Include)
    });
    add("exclude", &args.exclude, |v| {
        OptionValue::Filter(v, FilterLogic::Exclude)
    });
    add(
        "regions_overlap",
        &args.regions_overlap,
        OptionValue::RegionsOverlap,
    );
    add(
        "targets_overlap",
        &args.targets_overlap,
        OptionValue::TargetsOverlap,
    );
    add("write_index", &args.write_index, OptionValue::WriteIndex);
    // A count has a default value, with an index of its own: only real occurrences count.
    for (id, n, event) in [
        ("help", args.help, OptionValue::Help),
        ("quiet", args.quiet, OptionValue::Quiet),
        (
            "brief_predictions",
            args.brief_predictions,
            OptionValue::Brief,
        ),
    ] {
        if n > 0
            && let Some(idx) = matches.indices_of(id)
        {
            events.extend(idx.map(|i| (i, event)));
        }
    }
    events.sort_unstable_by_key(|e| e.0);
    events
}

/// Parses the options in command-line order, failing at the first bad one (or at `-h`) as
/// bcftools' `getopt` loop does.
fn parse_options(args: &Args, matches: &ArgMatches) -> anyhow::Result<ParsedOptions> {
    let mut p = ParsedOptions {
        verbosity: 1,
        brief: 0,
        greedy: false,
        phase: Phase::Require,
        ncsq2_max: 2 * (16 - 1),
        out_type: OutputTypeOption::default(),
        filter: None,
        threads: 0,
        regions_overlap: Overlap::Record,
        targets_overlap: Overlap::Pos,
        write_index: None,
    };
    for (_, event) in ordered_options(args, matches) {
        match event {
            OptionValue::Help => bail!("{USAGE}"),
            OptionValue::Quiet => {
                bail!("Error: the -q option has been deprecated, use -v, --verbose instead.")
            }
            OptionValue::Threads(t) => match strtol_strict(t) {
                // bcftools keeps the value in an int; a negative count fails later.
                Some(n) => p.threads = (n as i32).min(crate::output::max_threads()),
                None => bail!("Could not parse argument: --threads  {t}"),
            },
            OptionValue::Brief => {
                p.brief = 1;
                eprintln!(
                    "Warning: The -b option will be removed in future versions. Please use -B 1 instead."
                );
            }
            OptionValue::TrimProteinSeq(b) => match strtol_strict(b) {
                Some(n) if n >= 1 => p.brief = n as i32,
                _ => bail!("Could not parse argument: --trim-protein-seq {b}"),
            },
            OptionValue::Greedy(g) => match strtol_strict(g) {
                Some(n @ 0..=1) => p.greedy = n == 1,
                _ => bail!("Could not parse argument: --greedy {g}"),
            },
            OptionValue::Verbosity(v) => {
                let n = atoi(v) as i32;
                if n < 0 {
                    bail!("Error: expected integer with -v, --verbosity");
                }
                p.verbosity = n;
                if n > 3 {
                    log::set_verbosity(n);
                }
            }
            OptionValue::Phase(v) => {
                p.phase = match v.as_bytes().first() {
                    Some(b'a') => Phase::AsIs,
                    Some(b'm') => Phase::Merge,
                    Some(b'r') => Phase::Require,
                    Some(b'R') => Phase::NonRef,
                    Some(b's') => Phase::Skip,
                    _ => bail!("The -p code \"{v}\" not recognised"),
                }
            }
            OptionValue::Ncsq(n) => {
                p.ncsq2_max = (atoi(n) as i32).wrapping_mul(2);
                if p.ncsq2_max <= 0 {
                    bail!("Expected positive integer with -n, got {n}");
                }
            }
            OptionValue::OutputType(v) => p.out_type.apply(v)?,
            OptionValue::Filter(expr, logic) => {
                if p.filter.is_some() {
                    bail!(
                        "Error: only one -i or -e expression can be given, and they cannot be combined"
                    );
                }
                p.filter = Some((expr.to_owned(), logic));
            }
            OptionValue::RegionsOverlap(o) => {
                p.regions_overlap = parse_overlap(Some(o), Overlap::Record, "regions-overlap")?;
            }
            OptionValue::TargetsOverlap(o) => {
                p.targets_overlap = parse_overlap(Some(o), Overlap::Pos, "targets-overlap")?;
            }
            OptionValue::WriteIndex(w) => p.write_index = Some(output::parse_write_index(w)?),
        }
    }
    Ok(p)
}

/// The prefix of a chromosome name bcftools suggests stripping (bcftools' `test_prefix`).
fn test_prefix(seq: &[u8]) -> Vec<u8> {
    for p in [
        &b"chromosome_"[..],
        b"chromosome",
        b"chrom_",
        b"chrom",
        b"chr_",
        b"chr",
    ] {
        if seq.len() >= p.len() && seq[..p.len()].eq_ignore_ascii_case(p) {
            return seq[..p.len()].to_vec();
        }
    }
    b"-".to_vec()
}

/// The samples to process (bcftools' `smpl_ilist_init` with `SMPL_STRICT`): header indexes
/// in header order.
fn sample_list(header: &Header, list: Option<(&str, ListSource)>) -> anyhow::Result<Vec<usize>> {
    let n = header.sample_count();
    let Some((list, source)) = list else {
        return Ok((0..n).collect());
    };
    let (negate, body) = match list.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, list),
    };
    let Ok(names) = readlist::read_list(body.as_bytes(), source) else {
        bail!("Could not parse {list}");
    };
    let mut selected = vec![false; n];
    for name in &names {
        // the first column, up to an unescaped space
        let mut end = name.len();
        for (i, &c) in name.iter().enumerate() {
            if c.is_ascii_whitespace()
                && name[..i].iter().rev().take_while(|&&b| b == b'\\').count() % 2 == 0
            {
                end = i;
                break;
            }
        }
        let name = &name[..end];
        match header.sample_id(name) {
            Some(i) => selected[i as usize] = true,
            None => bail!("No such sample: \"{}\"", String::from_utf8_lossy(name)),
        }
    }
    Ok((0..n).filter(|&i| selected[i] != negate).collect())
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let ParsedOptions {
        verbosity,
        brief,
        greedy,
        phase,
        ncsq2_max,
        out_type,
        filter: filter_expr,
        threads,
        regions_overlap,
        targets_overlap,
        write_index,
    } = parse_options(&args, matches)?;
    let text = out_type.text;
    let out_settings = out_type.settings();

    let gencode = match gencode::choose(args.genetic_code.as_deref())? {
        GencodeChoice::Table(g) => g,
        GencodeChoice::List => {
            // bcftools prints the tables with printf, whose errors it ignores.
            let mut out = io::stdout().lock();
            let _ = gencode::print_tables(&mut out).and_then(|()| out.flush());
            return Ok(ExitCode::SUCCESS);
        }
    };

    let fname = match args.inputs.as_slice() {
        [] if !io::stdin().is_terminal() => "-".to_owned(),
        [f] => f.clone(),
        _ => bail!("{USAGE}"),
    };
    let Some(fa_path) = &args.fasta_ref else {
        bail!("Missing the --fasta-ref option");
    };
    let Some(gff_path) = &args.gff_annot else {
        bail!("Missing the --gff option");
    };

    let (mut reader, threads) = open_reader(
        &args,
        matches,
        &fname,
        (regions_overlap, targets_overlap),
        threads,
    )?;

    let opts = Opts {
        phase,
        verbosity,
        local_csq: args.local_csq,
        greedy,
        ncsq2_max,
        nfmt_bcsq: (1 + (ncsq2_max - 1) / BCSQ_BITS_PER_INT) as usize,
        brief_predictions: brief,
        force: args.force,
        bcsq_tag: args
            .custom_tag
            .as_deref()
            .unwrap_or("BCSQ")
            .as_bytes()
            .to_vec(),
        text,
    };
    let samples = later_file_sticky(
        matches,
        ("samples", args.samples.as_deref()),
        ("samples_file", args.samples_file.as_deref()),
    )
    .map(|(spec, source)| (spec.to_owned(), source));
    let mut csq = Csq::init(
        &args,
        opts,
        gencode,
        &mut reader,
        (fa_path, gff_path),
        filter_expr,
        samples,
        (&out_settings, threads, write_index),
    )?;

    if let Err(e) = csq.process_all(&mut reader) {
        // bcftools exits from error() without flushing its output.
        if let Output::Vcf { writer, .. } = csq.out {
            writer.abandon();
        }
        return Err(e);
    }
    csq.finish(&args, write_index)?;
    Ok(ExitCode::SUCCESS)
}

/// Opens the input with the `-r/-R` regions and `-t/-T` targets, and returns it with the
/// number of extra threads.
fn open_reader(
    args: &Args,
    matches: &ArgMatches,
    fname: &str,
    (regions_overlap, targets_overlap): (Overlap, Overlap),
    threads: i32,
) -> anyhow::Result<(SyncedReader, usize)> {
    let mut targets = None;
    let mut targets_exclude = false;
    if let Some((spec, source)) = later_file_sticky(
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
    let mut regions = None;
    if let Some((spec, source)) = later_file_sticky(
        matches,
        ("regions", args.regions.as_deref()),
        ("regions_file", args.regions_file.as_deref()),
    ) {
        let mut r = Regions::new(spec.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
        r.overlap = regions_overlap;
        regions = Some(r);
    }
    let Ok(threads) = usize::try_from(threads) else {
        bail!("Failed to create {threads} extra threads");
    };
    if !crate::output::threads_start(threads) {
        bail!("Failed to create {threads} extra threads");
    }
    let reader = match SyncedReader::open(
        fname,
        SyncedOptions {
            regions,
            targets,
            targets_exclude,
            apply_filters: None,
            threads,
        },
    ) {
        Ok(r) => r,
        Err(e) => {
            let name = if fname == "-" {
                "standard input"
            } else {
                fname
            };
            bail!("Failed to read from {name}: {e}");
        }
    };
    Ok((reader, threads))
}

impl Csq {
    /// Annotates the records of the input and flushes the buffers at the end.
    fn process_all(&mut self, reader: &mut SyncedReader) -> anyhow::Result<()> {
        let mut rec = Box::new(Record::new());
        loop {
            match reader.next_record_boxed(&mut rec) {
                Ok(true) => {
                    self.sync_header(reader.header())?;
                    self.process(&mut rec)?;
                }
                Ok(false) => break,
                Err(e) => bail!("Error: {e}"),
            }
        }
        if let Some(e) = reader.take_pending_error() {
            bail!("Error: {e}");
        }
        self.hap_flush(REGIDX_MAX)?;
        self.vbuf_flush(REGIDX_MAX)
    }

    /// Loads the reference, the GFF, the samples, and opens the output (bcftools'
    /// `init_data`).
    #[expect(
        clippy::too_many_arguments,
        reason = "the inputs of bcftools' init_data"
    )]
    fn init(
        args: &Args,
        mut opts: Opts,
        gencode: &'static Gencode,
        reader: &mut SyncedReader,
        (fa_path, gff_path): (&PathBuf, &PathBuf),
        filter_expr: Option<(String, FilterLogic)>,
        samples: Option<(String, ListSource)>,
        (out_settings, threads, write_index): (&[OutputSetting], usize, Option<output::WriteIndex>),
    ) -> anyhow::Result<Self> {
        let fai = Reference::open(fa_path)?;
        let gff = Gff::parse(
            gff_path,
            GffOptions {
                verbosity: opts.verbosity,
                force: opts.force,
                dump: args.dump_gff.clone(),
            },
        )?;
        let (chr_prefix, unify_err) =
            init_chr_names(args.unify_chr_names.as_deref(), reader.header(), &gff, &fai)?;
        let filter = match &filter_expr {
            Some((expr, logic)) => Some((filter::init(reader.header(), expr)?, *logic)),
            None => None,
        };
        let nsamples = reader.header().sample_count();
        if nsamples == 0 {
            opts.phase = Phase::DropGt;
        }
        let mut smpl = Vec::new();
        match &samples {
            Some((s, _)) if s == "-" => opts.phase = Phase::DropGt,
            _ => {
                smpl = sample_list(
                    reader.header(),
                    samples.as_ref().map(|(s, src)| (s.as_str(), *src)),
                )?;
            }
        }
        let hdr_nsmpl = if opts.phase == Phase::DropGt {
            0
        } else {
            nsamples
        };
        let output_name = args
            .output
            .as_deref()
            .map_or_else(|| "standard output".to_owned(), |p| p.display().to_string());
        // As in bcftools, the BCSQ lines go into the reader's header, so tags and contigs
        // that the VCF parser adds later (undeclared in the header) get the next IDs.
        let out = if opts.text {
            Output::Text(open_text_output(args)?)
        } else {
            open_vcf_output(
                args,
                &opts,
                reader.header_mut(),
                hdr_nsmpl,
                (out_settings, threads, write_index),
                output_name,
            )?
        };
        let hdr = dup_header(reader.header())?;
        if opts.verbosity > 0 {
            eprintln!("Calling...");
        }
        let tr_ends = gff.transcripts.iter().map(|t| t.end).collect();
        let ntr = gff.transcripts.len();
        Ok(Self {
            opts,
            gencode,
            gff,
            tr_ends,
            fai,
            fa_fname: fa_path.display().to_string(),
            gff_fname: gff_path.display().to_string(),
            hdr,
            nsamples,
            hdr_nsmpl,
            smpl,
            filter,
            out,
            chr_prefix,
            unify_err,
            prev_rid: -1,
            prev_pos: -1,
            rid: -1,
            chr_vcf: Vec::new(),
            chr_gff: Vec::new(),
            chr_fai: Vec::new(),
            active_tr: TrHeap::default(),
            aux: (0..ntr).map(|_| None).collect(),
            hap: Hap::default(),
            vbufs: VecDeque::new(),
            vbuf_pool: Vec::new(),
            next_vbuf_id: 0,
            pos2vbuf: HashMap::new(),
            rm_tr: Vec::new(),
            serial: 0,
            cur: CurRec::default(),
            gt: Genotypes::default(),
            ncsq2_small_warned: 0,
            warned: Warned::default(),
            str_buf: Vec::new(),
            fmt_i32: Vec::new(),
        })
    }

    /// Takes over the lines the VCF parser added to the reader's header for undeclared
    /// contigs and tags (bcftools formats records with the reader's header itself).
    fn sync_header(&mut self, reader_hdr: &Header) -> anyhow::Result<()> {
        if reader_hdr.record_count() != self.hdr.record_count() {
            self.hdr = dup_header(reader_hdr)?;
        }
        Ok(())
    }

    /// Converts a chromosome name between naming conventions (bcftools'
    /// `unify_chr_name`).
    fn unify_chr_name(
        &self,
        chr: &[u8],
        src: ChrSource,
        dst: ChrSource,
    ) -> anyhow::Result<Vec<u8>> {
        let (ps, pd) = (
            &self.chr_prefix[src as usize],
            &self.chr_prefix[dst as usize],
        );
        if ps.is_none() && pd.is_none() {
            return Ok(chr.to_vec());
        }
        let mut name = chr;
        if let Some(p) = ps {
            match chr.strip_prefix(p.as_slice()) {
                Some(rest) => name = rest,
                None => bail!(
                    "Error: failed to unify chr names, cannot strip \"{}\" from \"{}\"",
                    String::from_utf8_lossy(p),
                    String::from_utf8_lossy(chr)
                ),
            }
        }
        let mut out = pd.clone().unwrap_or_default();
        out.extend_from_slice(name);
        Ok(out)
    }

    /// Checks that a new chromosome is in the reference and the GFF.
    fn new_chromosome(&mut self, rec: &Record) -> anyhow::Result<()> {
        self.chr_vcf = rec.chrom_safe(&self.hdr).to_vec();
        // A common error is different naming conventions in the fasta and the VCF.
        let chr_fai = self.unify_chr_name(&self.chr_vcf, ChrSource::Vcf, ChrSource::Fai)?;
        if !self.fai.has_seq(&chr_fai) {
            missing_chr(
                &self.opts,
                &self.unify_err,
                &chr_fai,
                &self.fa_fname,
                &mut self.warned.missing_chr_fai,
            )?;
        }
        let chr_gff = self.unify_chr_name(&self.chr_vcf, ChrSource::Vcf, ChrSource::Gff)?;
        if !self.gff.has_seq(&chr_gff) {
            missing_chr(
                &self.opts,
                &self.unify_err,
                &chr_gff,
                &self.gff_fname,
                &mut self.warned.missing_chr_gff,
            )?;
        }
        self.chr_fai = chr_fai;
        self.chr_gff = chr_gff;
        Ok(())
    }

    /// Annotates and buffers one record (bcftools' `process`).
    fn process(&mut self, rec: &mut Record) -> anyhow::Result<()> {
        rec.unpack(Unpack::STR);
        let rid = rec.rid();
        let pos = rec.pos();
        if self.prev_rid != rid {
            self.prev_rid = rid;
            self.prev_pos = pos;
            self.new_chromosome(rec)?;
        }
        if self.prev_pos > pos {
            bail!(
                "Error: The file is not sorted, {}:{} comes before {}:{}",
                String::from_utf8_lossy(&self.chr_vcf),
                self.prev_pos + 1,
                String::from_utf8_lossy(rec.chrom_safe(&self.hdr)),
                pos + 1
            );
        }
        let n_allele = rec.allele_count();
        let mut call_csq = true;
        if n_allele < 2 {
            call_csq = false; // no alternate allele
        } else if n_allele == 2 {
            let a = rec.allele(1);
            if a.first() == Some(&b'*') || a.get(1) == Some(&b'*') {
                call_csq = false; // gVCF, not an alt allele
            }
        }
        if call_csq && let Some((f, logic)) = &mut self.filter {
            let pass = f.pass(&self.hdr, rec)?;
            call_csq = match logic {
                FilterLogic::Include => pass,
                FilterLogic::Exclude => !pass,
            };
        }
        // bcftools flushes the previous chromosome only before a record it annotates, so one
        // it does not annotate could be buffered with the previous chromosome's records, which
        // the buffers keyed by position mix up ("This should not happen"); flush it first.
        if self.rid != rid {
            self.hap_flush(REGIDX_MAX)?;
            self.vbuf_flush(REGIDX_MAX)?;
        }
        self.rid = rid;
        if !call_csq {
            if self.opts.text {
                return Ok(()); // not a VCF output
            }
            self.vbuf_push(rec);
            self.hap_flush((pos - 1) as u32)?;
            self.vbuf_flush((pos - 1) as u32)?;
            return Ok(());
        }
        let mut cur = std::mem::take(&mut self.cur);
        cur.fill(rec);
        cur.serial = self.vbuf_push(rec);
        self.gt.ret = None;
        let result = self.annotate(&cur);
        self.cur = cur;
        result?;
        if pos > 0 {
            self.hap_flush((pos - 1) as u32)?;
            self.vbuf_flush((pos - 1) as u32)?;
        }
        Ok(())
    }

    /// Looks up the current record in the GFF indexes.
    fn annotate(&mut self, cur: &CurRec) -> anyhow::Result<()> {
        if cur.allele_str(1).first() != Some(&b'<') {
            let mut hit = if self.opts.local_csq {
                self.test_cds_local(cur)?
            } else {
                self.test_cds(cur)?
            };
            hit |= self.test_utr(cur)?;
            hit |= self.test_splice(cur)?;
            if !hit || self.opts.greedy {
                self.test_tscript(cur)?;
            }
        } else {
            self.test_symbolic_alt(cur)?;
        }
        Ok(())
    }

    /// Closes the output (bcftools' `destroy_data`).
    fn finish(self, args: &Args, write_index: Option<output::WriteIndex>) -> anyhow::Result<()> {
        if self.ncsq2_small_warned != 0 {
            eprintln!(
                "Note: Some samples had too many consequences to be represented in {} bytes. If you need to record them all,\n      the limit can be increased by running with `--ncsq {}`.",
                self.opts.nfmt_bcsq / 8,
                1 + self.ncsq2_small_warned / 2
            );
        }
        let name = args
            .output
            .as_deref()
            .map_or_else(|| "stdout".to_owned(), |p| p.display().to_string());
        match self.out {
            Output::Vcf { writer, name: out } => match writer.finish() {
                Ok(()) => {}
                Err(WriteError::Index(_)) if write_index.is_some() => {
                    bail!("Error: cannot write to index {}", index_name(args))
                }
                Err(WriteError::Deferred(e)) => {
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[vbuf_flush] Error: cannot write to {out}")
                    }));
                }
                Err(WriteError::Io(e)) => {
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("Error: close failed .. {name}")
                    }));
                }
                Err(e) => {
                    return Err(
                        anyhow::Error::new(e).context(format!("Error: close failed .. {name}"))
                    );
                }
            },
            Output::Text(mut out) => out
                .flush()
                .with_context(|| format!("Error: close failed .. {name}"))?,
        }
        Ok(())
    }
}

/// Fails on a chromosome missing from the reference or the GFF `fname`, or with `--force`
/// warns about it once.
fn missing_chr(
    opts: &Opts,
    unify_err: &str,
    chr: &[u8],
    fname: &str,
    warned: &mut bool,
) -> anyhow::Result<()> {
    if !opts.force {
        bail!(
            "Error: the chromosome \"{}\" is not present in {fname}\n       {unify_err}\n",
            String::from_utf8_lossy(chr)
        );
    }
    if !*warned {
        *warned = true;
        eprintln!(
            "Warning: the chromosome \"{}\" is not present in {fname}. This warning is printed only once.",
            String::from_utf8_lossy(chr)
        );
    }
    Ok(())
}

/// The index file name, for messages.
fn index_name(args: &Args) -> String {
    let fmt = args
        .write_index
        .last()
        .and_then(|w| output::parse_write_index(w).ok())
        .unwrap_or(output::WriteIndex::Csi);
    let mut out_type = OutputTypeOption::default();
    for v in &args.output_type {
        if out_type.apply(v).is_err() {
            break;
        }
    }
    let ty = output::resolve(&out_type.settings(), args.output.as_deref())
        .map_or(pasta::variant::io::OutputType::Vcf, |(t, _)| t);
    output::index_request(fmt, ty, args.output.as_deref())
        .map_or_else(|| "(null)".to_owned(), |r| r.path.display().to_string())
}

/// Copies the reader's header for writing records. bcftools uses the reader's header itself,
/// so the warnings reparsing its lines would repeat (such as an invalid `-c` tag name) are
/// dropped.
fn dup_header(hdr: &Header) -> anyhow::Result<Header> {
    let (dup, _repeated) = log::capture(|| hdr.dup());
    dup.context("Failed to duplicate the header")
}

/// Opens the `-Ot` output and writes its header lines.
fn open_text_output(args: &Args) -> anyhow::Result<Box<dyn Write + Send>> {
    let mut out: Box<dyn Write + Send> = match &args.output {
        Some(p) => {
            let file = File::create(p).map_err(|e| {
                anyhow::anyhow!("Failed to write to {}: {}", p.display(), strerror(&e))
            })?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(output::stdout())),
    };
    writeln!(
        out,
        "# This file was produced by: vtools +csq({}+pasta-{})",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION")
    )?;
    let argv: Vec<String> = std::env::args().skip(1).collect();
    write!(
        out,
        "# The command line was:\tvtools +{}",
        argv.first().map_or("csq", |s| s)
    )?;
    for a in argv.iter().skip(1) {
        write!(out, " {a}")?;
    }
    writeln!(out)?;
    writeln!(out, "# LOG\t[2]Message")?;
    writeln!(
        out,
        "# CSQ\t[2]Sample\t[3]Haplotype\t[4]Chromosome\t[5]Position\t[6]Consequence"
    )?;
    Ok(out)
}

/// Opens the VCF/BCF output and writes the header with the BCSQ lines.
fn open_vcf_output(
    args: &Args,
    opts: &Opts,
    hdr: &mut Header,
    hdr_nsmpl: usize,
    (out_settings, threads, write_index): (&[OutputSetting], usize, Option<output::WriteIndex>),
    output_name: String,
) -> anyhow::Result<Output> {
    let (out_type, level) = output::resolve(out_settings, args.output.as_deref())?;
    let out_path = args.output.as_deref().map(output::data_path);
    let mut writer = match Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, threads),
    ) {
        Ok(w) => w,
        Err(e) => {
            log::error(
                "hts_open_format",
                format_args!("Failed to open file \"{output_name}\" : {}", strerror(&e)),
            );
            bail!(
                "[init_data] Error: cannot write to {output_name}: {}",
                strerror(&e)
            )
        }
    };
    if !args.no_version {
        version::append_version(hdr, "csq")?;
    }
    let tag = String::from_utf8_lossy(&opts.bcsq_tag).into_owned();
    let line = format!(
        "##INFO=<ID={tag},Number=.,Type=String,Description=\"{} consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change\">",
        if opts.local_csq {
            "Local"
        } else {
            "Haplotype-aware"
        }
    );
    hdr.append(line.as_bytes())
        .with_context(|| format!("Could not add the header line: {line}"))?;
    if hdr_nsmpl > 0 {
        let line = format!(
            "##FORMAT=<ID={tag},Number=.,Type=Integer,Description=\"Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \\\"bcftools query -f'[%CHROM\\t%POS\\t%SAMPLE\\t%TBCSQ\\n]'\\\" to translate.\">"
        );
        hdr.append(line.as_bytes())
            .with_context(|| format!("Could not add the header line: {line}"))?;
    }
    if let Some(fmt) = write_index {
        let ok = output::index_request(fmt, out_type, args.output.as_deref())
            .is_some_and(|req| writer.request_index(req).is_ok());
        if !ok {
            if writer.write_header(hdr).is_err() {
                bail!("[init_data] Error: cannot write the header to {output_name}");
            }
            writer.abandon();
            bail!(
                "Error: failed to initialise index for {}",
                args.output
                    .as_deref()
                    .map_or_else(|| "(null)".to_owned(), |p| p.display().to_string())
            );
        }
    }
    if writer.write_header(hdr).is_err() {
        bail!("[init_data] Error: cannot write the header to {output_name}");
    }
    Ok(Output::Vcf {
        writer: Box::new(writer),
        name: output_name,
    })
}

/// The chromosome prefixes of `--unify-chr-names`, or a suggestion for the error message
/// (bcftools' `init_chr_names`).
fn init_chr_names(
    unify: Option<&str>,
    header: &Header,
    gff: &Gff,
    fai: &Reference,
) -> anyhow::Result<([Option<Vec<u8>>; 3], String)> {
    let mut prefix: [Option<Vec<u8>>; 3] = [None, None, None];
    let mut err = "check if --unify-chr-names or --force could help".to_owned();
    if let Some(u) = unify
        && let Ok(list) = readlist::read_list(u.as_bytes(), ListSource::Inline)
    {
        if list.len() != 3 {
            bail!("Error: expected three strings, got --unify-chr-names {u}");
        }
        for (i, p) in list.into_iter().enumerate() {
            if p != b"-" {
                prefix[i] = Some(p);
            }
        }
        return Ok((prefix, err));
    }
    let names = header.contig_names();
    let Some(&seq_vcf) = names.first() else {
        return Ok((prefix, err));
    };
    let seq_gff = gff.seq_name(0).unwrap_or_default();
    let seq_fa = fai.seq_name(0).unwrap_or_default();
    if seq_vcf == seq_fa && seq_vcf == seq_gff {
        return Ok((prefix, err));
    }
    // Different prefixes or a different order: suggest --unify-chr-names if possible.
    let (pv, pg, pf) = (
        test_prefix(seq_vcf),
        test_prefix(seq_gff),
        test_prefix(seq_fa),
    );
    let strip = |s: &'_ [u8], p: &[u8]| -> Vec<u8> {
        if p == b"-" {
            s.to_vec()
        } else {
            s[p.len()..].to_vec()
        }
    };
    let (sv, sg, sf) = (strip(seq_vcf, &pv), strip(seq_gff, &pg), strip(seq_fa, &pf));
    let same_chr = sv == sg && sg == sf && sf == sv;
    let l = |s: &[u8]| String::from_utf8_lossy(s).into_owned();
    err = if same_chr {
        format!(
            "the first sequence name in VCF/GFF/fasta is {}/{}/{}, try to run with --unify-chr-names {},{},{}\n",
            l(seq_vcf),
            l(seq_gff),
            l(seq_fa),
            l(&pv),
            l(&pg),
            l(&pf)
        )
    } else {
        format!(
            "the first sequence name in VCF/GFF/fasta is {}/{}/{}, check if running with --unify-chr-names or --force coud help\n",
            l(seq_vcf),
            l(seq_gff),
            l(seq_fa)
        )
    };
    Ok((prefix, err))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use clap::{Args as _, Command, FromArgMatches as _};

    /// Runs `csq -p a` with `options` on `records` (CHROM, POS, ALT, QUAL), on two sequences
    /// with one protein-coding gene each; returns CHROM, POS, and INFO of each record.
    fn csq(
        name: &str,
        records: &[(&str, usize, &str, &str)],
        options: &[&str],
    ) -> anyhow::Result<Vec<String>> {
        let dir = std::env::temp_dir().join(format!("vtools-csq-{name}-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        // A CDS of 34 codons at 1-102, then intergenic sequence.
        let seq = format!("ATG{}TAA{}", "GCA".repeat(32), "C".repeat(198));
        let (fa, gff, vcf, out) = (
            dir.join("ref.fa"),
            dir.join("genes.gff3"),
            dir.join("in.vcf"),
            dir.join("out.vcf"),
        );
        let mut text = String::new();
        let mut genes = String::new();
        let mut lines = String::new();
        for chr in ["chr1", "chr2"] {
            text += &format!(">{chr}\n{seq}\n");
            genes += &format!(
                "{chr}\t.\tgene\t1\t102\t.\t+\t.\tID=gene:g{chr};Name=g{chr};biotype=protein_coding\n\
                 {chr}\t.\ttranscript\t1\t102\t.\t+\t.\tID=transcript:t{chr};Parent=gene:g{chr};biotype=protein_coding\n\
                 {chr}\t.\texon\t1\t102\t.\t+\t.\tParent=transcript:t{chr}\n\
                 {chr}\t.\tCDS\t1\t102\t.\t+\t0\tParent=transcript:t{chr}\n"
            );
        }
        for &(chr, pos, alt, qual) in records {
            let r = &seq[pos - 1..pos];
            lines += &format!("{chr}\t{pos}\t.\t{r}\t{alt}\t{qual}\t.\t.\tGT\t0|1\n");
        }
        fs::write(&fa, text)?;
        fs::write(&gff, genes)?;
        fs::write(
            &vcf,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=300>\n##contig=<ID=chr2,length=300>\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n{lines}"
            ),
        )?;
        let path = |p: &std::path::Path| p.to_str().expect("UTF-8 path").to_owned();
        let mut argv = vec!["csq".to_owned(), "--no-version".to_owned()];
        for (o, v) in [("-f", &fa), ("-g", &gff), ("-o", &out)] {
            argv.extend([o.to_owned(), path(v)]);
        }
        argv.extend(["-p", "a"].iter().chain(options).map(|s| (*s).to_owned()));
        argv.push(path(&vcf));
        let matches = super::Args::augment_args(Command::new("csq")).try_get_matches_from(&argv)?;
        let result = super::run(super::Args::from_arg_matches(&matches)?, &matches);
        let text = fs::read_to_string(&out).unwrap_or_default();
        fs::remove_dir_all(&dir)?;
        result?;
        Ok(text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<_> = l.split('\t').collect();
                format!("{}:{} {}", f[0], f[1], f[7])
            })
            .collect())
    }

    #[test]
    fn an_unannotated_record_of_a_new_chromosome_flushes_the_previous_one() -> anyhow::Result<()> {
        // chr2:10 is not annotated (no ALT, or excluded by -e) and shares its position with a
        // buffered chr1 record that is not the last: bcftools buffers it before flushing chr1
        // and stops with "This should not happen".
        for (alt, options) in [(".", &[][..]), ("T", &["-e", "QUAL=5"][..])] {
            let records = [
                ("chr1", 10, "A", "."),
                ("chr1", 20, "A", "."),
                ("chr2", 10, alt, "5"),
                ("chr2", 50, "A", "."),
            ];
            let lines = csq("chr", &records, options)?;
            assert_eq!(lines.len(), 4, "{lines:?}");
            assert_eq!(lines[2], "chr2:10 .");
            assert!(lines[0].starts_with("chr1:10 BCSQ="), "{lines:?}");
            assert!(lines[3].starts_with("chr2:50 BCSQ="), "{lines:?}");
        }
        Ok(())
    }
}
