//! `vtools consensus`: create a consensus sequence by applying VCF variants to a reference FASTA
//! (bcftools consensus, `consensus.c`).
//!
//! The structure follows bcftools: the FASTA is streamed line by line into a buffer of the
//! current sequence, the variants of that sequence are read through the index (each FASTA
//! entry sets the synced reader's regions to its sequence, or to `chr:from-to` for headers in
//! that form) and applied to the buffer once it covers them, and the buffer is written in lines
//! of 60 bases up to the last applied variant.

mod apply;
mod genotype;

use std::{
    fs::File,
    io::{self, BufRead, BufWriter, IsTerminal, Write},
    path::Path,
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, Record, Unpack,
    index::split_idx_name,
    io::InputFormat,
    log,
    readlist::{self, ListSource},
    regidx::RegIdx,
    synced::{IndexRequirement, MAX_CSI_COOR, Overlap, Regions, SyncedReaders},
};

use crate::{
    commands::head::apply_verbosity,
    filter::{self, FilterLogic, RecordFilter},
    raw_input::strerror,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file, BGZF-compressed and indexed; further arguments are ignored, as in
    /// bcftools.
    inputs: Vec<String>,

    /// Write a chain file for liftover.
    #[arg(short = 'c', long = "chain", value_name = "FILE")]
    chain: Option<String>,

    /// Replace positions absent from VCF with CHAR.
    #[arg(
        short = 'a',
        long = "absent",
        value_name = "CHAR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    absent: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Reference sequence in fasta format.
    #[arg(short = 'f', long = "fasta-ref", value_name = "FILE")]
    fasta_ref: Option<String>,

    /// Choose which allele to use from the FORMAT/GT field: N (1, 2, ..), R, A, I, LR, LA,
    /// SR, SA, or NpIu (e.g. 2pIu); case-insensitive.
    #[arg(
        short = 'H',
        long = "haplotype",
        value_name = "WHICH",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    haplotype: Vec<String>,

    /// Select sites for which the expression is true.
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

    /// Output IUPAC codes based on FORMAT/GT, use -s/-S to subset samples.
    #[arg(short = 'I', long = "iupac-codes")]
    iupac_codes: bool,

    /// Instead of removing sequence, insert character CHAR for deletions.
    #[arg(
        long = "mark-del",
        value_name = "CHAR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    mark_del: Vec<String>,

    /// Highlight insertions in uppercase (uc), lowercase (lc), or use CHAR.
    #[arg(
        long = "mark-ins",
        value_name = "uc|lc|CHAR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    mark_ins: Vec<String>,

    /// Highlight substitutions in uppercase (uc), lowercase (lc), or use CHAR.
    #[arg(
        long = "mark-snv",
        value_name = "uc|lc|CHAR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    mark_snv: Vec<String>,

    /// Replace regions according to the next --mask-with option [N].
    #[arg(short = 'm', long = "mask", value_name = "FILE", action = clap::ArgAction::Append)]
    mask: Vec<String>,

    /// Replace with CHAR (skips overlapping variants); change to uppercase (uc) or lowercase
    /// (lc).
    #[arg(
        long = "mask-with",
        value_name = "CHAR|uc|lc",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    mask_with: Vec<String>,

    /// Output CHAR instead of skipping a missing genotype "./.".
    #[arg(
        short = 'M',
        long = "missing",
        value_name = "CHAR",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    missing: Vec<String>,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<String>,

    /// Prefix to add to output sequence names.
    #[arg(
        short = 'p',
        long = "prefix",
        value_name = "STRING",
        allow_hyphen_values = true
    )]
    prefix: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(
        long = "regions-overlap",
        value_name = "0|1|2",
        action = clap::ArgAction::Append
    )]
    regions_overlap: Vec<String>,

    /// Comma-separated list of samples to include, "-" to ignore samples and use REF,ALT.
    #[arg(
        short = 's',
        long = "samples",
        visible_alias = "sample",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to include.
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
    samples_file: Option<String>,

    /// Use INT threads to decompress the input (a vtools addition) [0].
    #[arg(long = "threads", value_name = "INT")]
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

    /// Print help.
    #[arg(short = 'h', long = "help", short_alias = '?', action = clap::ArgAction::Count)]
    help: u8,

    /// In bcftools' option string without a meaning: prints the usage.
    #[arg(short = '1', hide = true, action = clap::ArgAction::Count)]
    one: u8,
}

const USAGE: &str = "
About: Create consensus sequence by applying VCF variants to a reference fasta
       file. By default, the program will apply all ALT variants. Using the
       --samples (and, optionally, --haplotype) option will apply genotype
       (or haplotype) calls from FORMAT/GT. The program ignores allelic depth
       information, such as INFO/AD or FORMAT/AD.
Usage: vtools consensus [OPTIONS] <file.vcf.gz>
Options:
    -c, --chain FILE               Write a chain file for liftover
    -a, --absent CHAR              Replace positions absent from VCF with CHAR
    -e, --exclude EXPR             Exclude sites for which the expression is true (see man page for details)
    -f, --fasta-ref FILE           Reference sequence in fasta format
    -H, --haplotype WHICH          Choose which allele to use from the FORMAT/GT field, note
                                   the codes are case-insensitive:
                                       N: N={1,2,3,..} is the index of the allele from GT, regardless of phasing (e.g. \"2\")
                                       R: REF allele in het genotypes
                                       A: ALT allele
                                       I: IUPAC code for all genotypes
                                       LR,LA: longer allele and REF/ALT if equal length
                                       SR,SA: shorter allele and REF/ALT if equal length
                                       NpIu: index of the allele for phased and IUPAC code for unphased GTs (e.g. \"2pIu\")
    -i, --include EXPR             Select sites for which the expression is true (see man page for details)
    -I, --iupac-codes              Output IUPAC codes based on FORMAT/GT, use -s/-S to subset samples
        --mark-del CHAR            Instead of removing sequence, insert character CHAR for deletions
        --mark-ins uc|lc|CHAR      Highlight insertions in uppercase (uc), lowercase (lc), or use CHAR, leaving the rest as is
        --mark-snv uc|lc|CHAR      Highlight substitutions in uppercase (uc), lowercase (lc), or use CHAR, leaving the rest as is
    -m, --mask FILE                Replace regions according to the next --mask-with option. The default is --mask-with N
        --mask-with CHAR|uc|lc     Replace with CHAR (skips overlapping variants); change to uppercase (uc) or lowercase (lc)
    -M, --missing CHAR             Output CHAR instead of skipping a missing genotype \"./.\"
    -o, --output FILE              Write output to a file [standard output]
    -p, --prefix STRING            Prefix to add to output sequence names
        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --samples LIST             Comma-separated list of samples to include, \"-\" to ignore samples and use REF,ALT
    -S, --samples-file FILE        File of samples to include
        --threads INT              Use INT threads to decompress the input (a vtools addition) [0]
    -v, --verbosity INT            Verbosity level
Examples:
   # Get the consensus for one region. The fasta header lines are then expected
   # in the form \">chr:from-to\".
   samtools faidx ref.fa 8:11870-11890 | vtools consensus in.vcf.gz > out.fa

   # See also http://samtools.github.io/bcftools/howtos/consensus-sequence.html
";

/// Prints the usage and returns bcftools' exit status for it.
fn usage() -> ExitCode {
    eprintln!("{USAGE}");
    ExitCode::from(1)
}

/// How `--mark-ins` and `--mark-snv` highlight bases (bcftools' `TO_UPPER`, `TO_LOWER`, or a
/// character).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Mark {
    Upper,
    Lower,
    Char(u8),
}

impl Mark {
    fn parse(value: &str, option: &str) -> anyhow::Result<Self> {
        let b = value.as_bytes();
        if value.eq_ignore_ascii_case("uc") {
            Ok(Self::Upper)
        } else if value.eq_ignore_ascii_case("lc") {
            Ok(Self::Lower)
        } else if b.len() == 1 && b[0] > 32 && b[0] < 127 {
            Ok(Self::Char(b[0]))
        } else {
            bail!("The argument is not recognised: --{option} {value}")
        }
    }
}

/// What a mask does to the bases it covers (bcftools' `mask_t.with`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MaskWith {
    Upper,
    Lower,
    /// Replace with a character; variants overlapping the mask are skipped (`MASK_SKIP`).
    Char(u8),
}

/// One `-m` file and its `--mask-with` (bcftools' `mask_t`).
struct Mask {
    idx: RegIdx<()>,
    with: MaskWith,
}

/// `-H` bits (bcftools' `PICK_*`).
const PICK_REF: u8 = 1;
const PICK_ALT: u8 = 2;
const PICK_LONG: u8 = 4;
const PICK_SHORT: u8 = 8;
const PICK_IUPAC: u8 = 16;

/// The options that drive variant application.
#[derive(Debug, Default)]
struct Options {
    /// `-I`.
    output_iupac: bool,
    /// The `-H` haplotype number (0 when none).
    haplotype: i64,
    /// The `-H` `PICK_*` bits.
    allele: u8,
    missing: Option<u8>,
    absent: Option<u8>,
    mark_del: Option<u8>,
    mark_ins: Option<Mark>,
    mark_snv: Option<Mark>,
    prefix: Vec<u8>,
    regions_overlap: Overlap,
}

/// Parses one `-H` value into `opts` (bcftools' `case 'H'`).
fn parse_haplotype(opts: &mut Options, value: &str) -> anyhow::Result<()> {
    let pick = match value.to_ascii_lowercase().as_str() {
        "r" => Some(PICK_REF),
        "a" => Some(PICK_ALT),
        "l" | "lr" => Some(PICK_LONG | PICK_REF),
        "s" | "sr" => Some(PICK_SHORT | PICK_REF),
        "la" => Some(PICK_LONG | PICK_ALT),
        "sa" => Some(PICK_SHORT | PICK_ALT),
        "i" => Some(PICK_IUPAC),
        "1piu" => {
            opts.haplotype = 1;
            Some(PICK_IUPAC)
        }
        "2piu" => {
            opts.haplotype = 2;
            Some(PICK_IUPAC)
        }
        _ => None,
    };
    if let Some(p) = pick {
        opts.allele |= p;
        return Ok(());
    }
    let b = value.as_bytes();
    let (n, end) = pasta::variant::number::strtol(b, 0);
    let rest = &value[end..];
    if end == 0 || (!rest.is_empty() && !rest.eq_ignore_ascii_case("pIu")) {
        bail!(
            "Error: Could not parse \"--haplotype {value}\", expected number of number followed with \"pIu\""
        );
    }
    // C stores the value in an int.
    opts.haplotype = i64::from(n as i32);
    if !rest.is_empty() {
        opts.allele |= PICK_IUPAC;
    }
    if opts.haplotype <= 0 {
        bail!("Error: Expected positive integer with --haplotype");
    }
    Ok(())
}

/// The first character of a single-character option (`-a`, `-M`).
fn single_char(value: &str, option: char) -> anyhow::Result<Option<u8>> {
    let b = value.as_bytes();
    if b.len() > 1 {
        bail!("Expected single character with -{option}, got \"{value}\"");
    }
    Ok(b.first().copied())
}

/// One option occurrence whose value bcftools checks while parsing the command line.
#[derive(Clone, Copy, Debug)]
enum OptionValue<'a> {
    Help,
    Verbosity(&'a str),
    MarkDel(&'a str),
    MarkIns(&'a str),
    MarkSnv(&'a str),
    RegionsOverlap(&'a str),
    Filter(&'a str, FilterLogic),
    Mask(&'a str),
    MaskWith(&'a str),
    Absent(&'a str),
    Missing(&'a str),
    Haplotype(&'a str),
}

/// The options parsed from the command line, besides the plain values of [`Args`].
#[derive(Debug)]
struct ParsedOptions {
    opts: Options,
    /// The `-i/-e` expression.
    filter: Option<(String, FilterLogic)>,
    /// The `-m` files with their `--mask-with`.
    masks: Vec<(String, MaskWith)>,
    /// `-h` was given (bcftools prints the usage when it reaches it).
    help: bool,
}

/// The occurrences of the options bcftools validates as it reads them, in command-line order.
fn ordered_options<'a>(args: &'a Args, matches: &ArgMatches) -> Vec<(usize, OptionValue<'a>)> {
    let mut events: Vec<(usize, OptionValue<'a>)> = Vec::new();
    let mut add = |id: &str, values: &'a [String], make: fn(&'a str) -> OptionValue<'a>| {
        if let Some(idx) = matches.indices_of(id) {
            events.extend(idx.zip(values).map(|(i, v)| (i, make(v))));
        }
    };
    add("verbosity", &args.verbosity, OptionValue::Verbosity);
    add("mark_del", &args.mark_del, OptionValue::MarkDel);
    add("mark_ins", &args.mark_ins, OptionValue::MarkIns);
    add("mark_snv", &args.mark_snv, OptionValue::MarkSnv);
    add(
        "regions_overlap",
        &args.regions_overlap,
        OptionValue::RegionsOverlap,
    );
    add("include", &args.include, |v| {
        OptionValue::Filter(v, FilterLogic::Include)
    });
    add("exclude", &args.exclude, |v| {
        OptionValue::Filter(v, FilterLogic::Exclude)
    });
    add("mask", &args.mask, OptionValue::Mask);
    add("mask_with", &args.mask_with, OptionValue::MaskWith);
    add("absent", &args.absent, OptionValue::Absent);
    add("missing", &args.missing, OptionValue::Missing);
    add("haplotype", &args.haplotype, OptionValue::Haplotype);
    // A count has a default value, with an index of its own: only real occurrences count.
    for (id, n) in [("help", args.help), ("one", args.one)] {
        if n > 0
            && let Some(idx) = matches.indices_of(id)
        {
            events.extend(idx.map(|i| (i, OptionValue::Help)));
        }
    }
    events.sort_unstable_by_key(|e| e.0);
    events
}

/// Parses the options in command-line order, failing at the first bad one as bcftools'
/// `getopt` loop does. Stops at `-h`.
fn parse_options(args: &Args, matches: &ArgMatches) -> anyhow::Result<ParsedOptions> {
    let mut parsed = ParsedOptions {
        opts: Options {
            regions_overlap: Overlap::Record,
            prefix: args.prefix.clone().unwrap_or_default().into_bytes(),
            output_iupac: args.iupac_codes,
            ..Options::default()
        },
        filter: None,
        masks: Vec::new(),
        help: false,
    };
    let opts = &mut parsed.opts;
    for (_, event) in ordered_options(args, matches) {
        match event {
            OptionValue::Help => {
                parsed.help = true;
                break;
            }
            OptionValue::Verbosity(v) => apply_verbosity(v)?,
            OptionValue::MarkDel(m) => opts.mark_del = m.bytes().next(),
            OptionValue::MarkIns(m) => opts.mark_ins = Some(Mark::parse(m, "mark-ins")?),
            OptionValue::MarkSnv(m) => opts.mark_snv = Some(Mark::parse(m, "mark-snv")?),
            OptionValue::RegionsOverlap(o) => {
                opts.regions_overlap = Overlap::parse(o)
                    .with_context(|| format!("Could not parse: --regions-overlap {o}"))?;
            }
            OptionValue::Filter(expr, logic) => {
                if parsed.filter.is_some() {
                    bail!(
                        "Error: only one -i or -e expression can be given, and they cannot be combined"
                    );
                }
                parsed.filter = Some((expr.to_owned(), logic));
            }
            OptionValue::Mask(file) => parsed.masks.push((file.to_owned(), MaskWith::Char(b'N'))),
            OptionValue::MaskWith(value) => {
                let Some(last) = parsed.masks.last_mut() else {
                    bail!("The --mask-with option must follow --mask");
                };
                last.1 = if value.eq_ignore_ascii_case("uc") {
                    MaskWith::Upper
                } else if value.eq_ignore_ascii_case("lc") {
                    MaskWith::Lower
                } else if value.len() == 1 {
                    MaskWith::Char(value.as_bytes()[0])
                } else {
                    bail!(
                        "Expected \"lc\", \"uc\", or a single character with the --mask-with option"
                    );
                };
            }
            OptionValue::Absent(a) => opts.absent = single_char(a, 'a')?,
            OptionValue::Missing(m) => opts.missing = single_char(m, 'M')?,
            OptionValue::Haplotype(h) => parse_haplotype(opts, h)?,
        }
    }
    Ok(parsed)
}

/// The samples a list selects (bcftools' `smpl_ilist_init` with `SMPL_NONE|SMPL_VERBOSE`):
/// their header indexes in VCF order, and the number of list entries that matched (bcftools'
/// `smpl->n`, which counts repeated names).
fn sample_list(
    header: &Header,
    spec: &str,
    source: ListSource,
) -> anyhow::Result<(Vec<usize>, usize)> {
    let nsmpl = header.sample_count();
    let (negate, body) = match spec.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, spec),
    };
    let list = readlist::read_list(body.as_bytes(), source)
        .map_err(|_| anyhow::anyhow!("Could not parse {spec}"))?;
    let mut selected = vec![false; nsmpl];
    let mut n = 0usize;
    for item in &list {
        // The name ends at the first unescaped whitespace.
        let mut end = item.len();
        for (i, &c) in item.iter().enumerate() {
            if c.is_ascii_whitespace() {
                let escapes = item[..i].iter().rev().take_while(|&&b| b == b'\\').count();
                if escapes % 2 == 0 {
                    end = i;
                    break;
                }
            }
        }
        let name = &item[..end];
        match header.sample_id(name) {
            Some(i) => {
                selected[i as usize] = true;
                n += 1;
            }
            None => eprintln!("No such sample: \"{}\"", String::from_utf8_lossy(name)),
        }
    }
    let idx: Vec<usize> = (0..nsmpl).filter(|&i| selected[i] != negate).collect();
    if negate {
        n = nsmpl.saturating_sub(n);
    }
    Ok((idx, n))
}

/// The VCF input: the synced reader, the `-i/-e` filter, and the record read ahead but not
/// yet applied (bcftools' `vcf_rbuf`, which never holds more than one: a record is unread
/// only to stop reading).
struct VcfInput {
    sr: SyncedReaders,
    filter: Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    pending: Option<Record>,
    pool: Vec<Record>,
}

impl VcfInput {
    /// Returns the next record without taking it, reading one that passes `-i/-e` when none is
    /// pending (bcftools' `next_vcf_line`; leaving the record is `unread_vcf_line`). `None` at
    /// the end of the current regions.
    fn peek(&mut self) -> anyhow::Result<Option<&Record>> {
        if self.pending.is_none() {
            let mut rec = self.pool.pop().unwrap_or_default();
            if self.read(&mut rec)? {
                self.pending = Some(rec);
            } else {
                self.pool.push(rec);
            }
        }
        Ok(self.pending.as_ref())
    }

    fn read(&mut self, rec: &mut Record) -> anyhow::Result<bool> {
        loop {
            let n = self
                .sr
                .next_line()
                .map_err(|e| anyhow::anyhow!("Error: {e}"))?;
            if n == 0 {
                if let Some(e) = self.sr.take_pending_error() {
                    bail!("Error: {e}");
                }
                return Ok(false);
            }
            self.sr.swap_line(0, rec);
            if let Some((f, logic)) = &mut self.filter {
                let pass = f.pass(self.sr.header(0), rec)?;
                if pass != (*logic == FilterLogic::Include) {
                    continue;
                }
            }
            return Ok(true);
        }
    }

    /// Takes the record [`Self::peek`] returned.
    fn take(&mut self) -> Record {
        self.pending.take().expect("take follows a successful peek")
    }

    /// Returns a record's buffers for reuse.
    fn recycle(&mut self, rec: Record) {
        self.pool.push(rec);
    }

    /// Drops the pending record (bcftools' `vcf_rbuf.n = 0`).
    fn clear(&mut self) {
        if let Some(rec) = self.pending.take() {
            self.pool.push(rec);
        }
    }
}

/// An ungapped block of a chain and the gaps that follow it.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct ChainBlock {
    len: i64,
    ref_gap: i64,
    alt_gap: i64,
}

/// The ungapped blocks between the reference and the consensus of one sequence (bcftools'
/// `chain_t`).
#[derive(Debug, Default)]
struct Chain {
    blocks: Vec<ChainBlock>,
    ori_pos: i64,
    ref_last_block_ori: i64,
    alt_last_block_ori: i64,
}

impl Chain {
    fn new(ori_pos: i64) -> Self {
        Self {
            blocks: Vec::new(),
            ori_pos,
            ref_last_block_ori: ori_pos,
            alt_last_block_ori: ori_pos,
        }
    }

    /// Records a gap (bcftools' `push_chain_gap`).
    fn push_gap(&mut self, ref_start: i64, ref_len: i64, alt_start: i64, alt_len: i64) {
        match self.blocks.last_mut() {
            Some(last) if ref_start <= self.ref_last_block_ori => {
                // Back-to-back with the previous variant.
                last.ref_gap += ref_len;
                last.alt_gap += alt_len;
            }
            _ => self.blocks.push(ChainBlock {
                len: ref_start - self.ref_last_block_ori,
                ref_gap: ref_len,
                alt_gap: alt_len,
            }),
        }
        self.ref_last_block_ori = ref_start + ref_len;
        self.alt_last_block_ori = alt_start + alt_len;
    }
}

/// The consensus state (bcftools' `args_t`, minus the VCF reader).
struct Consensus {
    opts: Options,
    header: Header,
    /// The samples whose genotypes give IUPAC codes (without `-H`).
    iupac_samples: Option<Vec<usize>>,
    /// The sample of `-H`.
    isample: Option<usize>,
    /// The header ID of FORMAT/GT.
    gt_id: Option<i32>,

    /// The buffered sequence (bcftools' `fa_buf`).
    fa_buf: Vec<u8>,
    /// The longest FASTA line read so far, which sets the buffer capacity kept between
    /// sequences.
    longest_line: usize,
    /// Start of the buffer in the original sequence.
    fa_ori_pos: i64,
    /// The last protected position, to avoid applying overlapping variants.
    fa_frz_pos: i64,
    /// Length difference between the modified and original sequence up to `fa_frz_pos`.
    fa_mod_off: i64,
    /// The buffer offset of `fa_frz_pos`.
    fa_frz_mod: i64,
    /// With `--absent`, the buffer from `fa_frz_mod` up to here already holds the absent
    /// character, so repeated fills skip it (bcftools refills it for every FASTA line).
    absent_filled: usize,
    /// The region's start in the original sequence.
    fa_beg_pos: i64,
    /// The region's end in the original sequence (0 when unknown).
    fa_end_pos: i64,
    /// The length of the sequence read so far.
    fa_length: i64,
    /// The last position read from the FASTA.
    fa_src_pos: i64,
    /// The base before an insertion or SNP that follows a deletion, to check its REF.
    prev_base: u8,
    prev_base_pos: i64,
    prev_is_insert: bool,

    rid: i32,
    chr: Vec<u8>,
    masks: Vec<Mask>,
    chain: Option<Chain>,
    chain_id: i64,
    chain_out: Option<BufWriter<File>>,
    chain_name: String,
    out: BufWriter<Box<dyn Write>>,
    out_name: String,
    /// Output lines being written (reused).
    lines_out: Vec<u8>,
    napplied: i64,
    warned_haplotype: bool,
    warned_multibase_del: bool,
    scratch: apply::Scratch,
}

/// `--absent` positions up to the end (bcftools' `HTS_POS_MAX`).
const POS_MAX: i64 = i64::MAX;

/// FASTA output line width.
const LINE_WIDTH: usize = 60;

/// Output lines are formatted in batches of this many bytes, which stay in cache.
const FLUSH_BATCH: usize = 1024 * LINE_WIDTH;

/// A buffer longer than this is written whole when no record left can change it. Shorter ones
/// follow bcftools' protected position, so this also bounds how much more than bcftools is
/// written before an error.
const FLUSH_LONG: usize = 1024 * LINE_WIDTH;

/// The buffer capacity kept for the next sequence is at least this; a larger one, left by a
/// long record, is released.
const FA_BUF_KEEP: usize = 1 << 20;

impl Consensus {
    fn write_out(&mut self, data: &[u8]) -> anyhow::Result<()> {
        self.out
            .write_all(data)
            .with_context(|| format!("Could not write: {}", self.out_name))
    }

    /// Writes the buffer in lines of 60: whole lines with [`FlushMode::Lines`], all of it with
    /// [`FlushMode::All`] (bcftools' `flush_fa_buffer`). [`FlushMode::Lines`] is only for when
    /// no record left can change the buffer.
    fn flush_fa_buffer(&mut self, mode: FlushMode) -> anyhow::Result<()> {
        if self.fa_buf.is_empty() {
            return Ok(());
        }
        let len = self.fa_buf.len();
        let nlines = match mode {
            FlushMode::Lines => {
                // bcftools writes the lines before the protected position, so a whole sequence
                // is held when it has no variants.
                let frz = self.fa_frz_pos - self.fa_ori_pos;
                let protected = if frz <= 0 {
                    0
                } else {
                    ((frz - 1) / LINE_WIDTH as i64) as usize
                };
                // A long buffer is written whole, keeping a base, so that a record past it
                // finds it non-empty and is skipped as in bcftools' longer buffer. Short ones
                // are not, so that what is written before an error matches bcftools.
                let long = if len > FLUSH_LONG {
                    (len - 1) / LINE_WIDTH
                } else {
                    0
                };
                (len / LINE_WIDTH).min(protected).max(long)
            }
            FlushMode::All => len / LINE_WIDTH,
        };
        let nwr = nlines * LINE_WIDTH;
        let mut lines = std::mem::take(&mut self.lines_out);
        let mut written = Ok(());
        for block in self.fa_buf[..nwr].chunks(FLUSH_BATCH) {
            lines.clear();
            lines.resize(block.len() / LINE_WIDTH * (LINE_WIDTH + 1), b'\n');
            let (dsts, _) = lines.as_chunks_mut::<{ LINE_WIDTH + 1 }>();
            let (srcs, _) = block.as_chunks::<LINE_WIDTH>();
            for (dst, src) in dsts.iter_mut().zip(srcs) {
                dst[..LINE_WIDTH].copy_from_slice(src);
            }
            written = self
                .out
                .write_all(&lines)
                .with_context(|| format!("Could not write: {}", self.out_name));
            if written.is_err() {
                break;
            }
        }
        self.lines_out = lines;
        written?;
        self.fa_ori_pos += nwr as i64;
        self.fa_frz_mod -= nwr as i64;
        self.absent_filled = self.absent_filled.saturating_sub(nwr);
        if mode == FlushMode::Lines {
            // Not finished with this sequence yet: keep the rest.
            self.fa_buf.drain(..nwr);
            return Ok(());
        }
        if nwr == len {
            self.fa_buf.clear();
            return Ok(());
        }
        let rest = std::mem::take(&mut self.fa_buf);
        self.write_out(&rest[nwr..])?;
        self.write_out(b"\n")?;
        self.fa_buf = rest;
        self.fa_ori_pos += (len - nwr) as i64 - self.fa_mod_off;
        self.fa_mod_off = 0;
        self.fa_buf.clear();
        self.absent_filled = 0;
        Ok(())
    }

    /// Fills the gap since the last applied variant up to `pos` with the `--absent` character
    /// (bcftools' `apply_absent`).
    fn apply_absent(&mut self, pos: i64) {
        let Some(absent) = self.opts.absent else {
            return;
        };
        // bcftools also returns for `pos <= fa_ori_pos`: that is the region start until lines
        // are written, and then lies before `fa_frz_pos`. With lines written past it, and the
        // offset of an insertion, `fa_ori_pos` can pass a later variant.
        if self.fa_buf.is_empty() || pos <= self.fa_frz_pos + 1 || pos <= self.fa_beg_pos {
            return;
        }
        let len = self.fa_buf.len() as i64;
        let off = (pos - self.fa_ori_pos).saturating_add(self.fa_mod_off);
        let ie = if pos != 0 && off < len { off } else { len };
        let ib = self.fa_frz_mod.max(0) as usize;
        let ie = ie as usize;
        let start = ib.max(self.absent_filled);
        if start < ie {
            self.fa_buf[start..ie].fill(absent);
            self.absent_filled = ie;
        }
    }

    /// Protects the reference allele of `pos..pos+rlen` (bcftools' `freeze_ref`).
    fn freeze_ref(&mut self, pos: i64, rlen: i64) {
        if self.fa_frz_pos >= pos + rlen - 1 {
            return;
        }
        self.fa_frz_pos = pos + rlen - 1;
        self.fa_frz_mod = pos - self.fa_ori_pos + self.fa_mod_off + rlen;
    }

    /// Writes the chain of the finished sequence (bcftools' `print_chain`).
    fn print_chain(&mut self) -> anyhow::Result<()> {
        let (Some(chain), Some(out)) = (&self.chain, &mut self.chain_out) else {
            return Ok(());
        };
        let ref_end_pos = self.fa_length + chain.ori_pos;
        let last_block_size = ref_end_pos - chain.ref_last_block_ori;
        let alt_end_pos = chain.alt_last_block_ori + last_block_size;
        let score: i64 = chain.blocks.iter().map(|b| b.len).sum::<i64>() + last_block_size;
        self.chain_id += 1;
        let chr = String::from_utf8_lossy(&self.chr);
        let ori = chain.ori_pos;
        let id = self.chain_id;
        let mut write = || -> io::Result<()> {
            writeln!(
                out,
                "chain {score} {chr} {ref_end_pos} + {ori} {ref_end_pos} {chr} {alt_end_pos} + {ori} {alt_end_pos} {id}"
            )?;
            for b in &chain.blocks {
                writeln!(out, "{} {} {}", b.len, b.ref_gap, b.alt_gap)?;
            }
            write!(out, "{last_block_size}\n\n")
        };
        write().with_context(|| format!("Could not write: {}", self.chain_name))
    }

    /// Starts a new FASTA entry (bcftools' `init_region`): `line` is the header without `>`.
    fn init_region(&mut self, input: &mut VcfInput, line: &[u8]) -> anyhow::Result<()> {
        let (name, from, to) = parse_fasta_header(line);
        self.chr.clear();
        self.chr.extend_from_slice(name);
        self.rid = self.header.contig_id(name).map_or(-1, |r| r as i32);
        if self.rid < 0 {
            eprintln!(
                "Warning: Sequence \"{}\" not in {}",
                String::from_utf8_lossy(name),
                input.sr.file_name(0)
            );
        }
        self.prev_base_pos = -1;
        self.fa_buf.clear();
        // Do not hold a long record's high-water mark for the rest of the run. Lines as long as
        // the longest so far are likely to come again (an unwrapped FASTA), and releasing room
        // for them only makes the next sequence fault it in again: a record across a line end
        // holds two, and growth doubles that, so room for four is kept.
        let keep = FA_BUF_KEEP.max(4 * self.longest_line);
        if self.fa_buf.capacity() > keep {
            self.fa_buf.shrink_to(keep);
        }
        self.fa_length = 0;
        self.fa_beg_pos = from;
        self.fa_end_pos = to;
        self.fa_ori_pos = from;
        self.fa_src_pos = from;
        self.fa_mod_off = 0;
        self.fa_frz_pos = -1;
        self.fa_frz_mod = -1;
        self.absent_filled = 0;
        input.clear();

        // The synced reader's regions are 1-based.
        let beg = if from == 0 { 1 } else { from + 1 };
        let end = if to == 0 { MAX_CSI_COOR - 1 } else { to + 1 };
        let mut regions = Regions::empty();
        regions.add(name, beg, end);
        regions.sort_and_merge();
        input.sr.set_regions_overlap(self.opts.regions_overlap);
        input.sr.set_regions(regions);

        let mut hdr_line = Vec::with_capacity(line.len() + self.opts.prefix.len() + 2);
        hdr_line.push(b'>');
        hdr_line.extend_from_slice(&self.opts.prefix);
        // bcftools prints the line as a C string.
        hdr_line.extend_from_slice(line.split(|&c| c == 0).next().unwrap_or_default());
        hdr_line.push(b'\n');
        self.write_out(&hdr_line)?;
        if self.chain_out.is_some() {
            self.chain = Some(Chain::new(self.fa_ori_pos));
        }
        Ok(())
    }

    /// Applies the masks to a FASTA line of `len` bases that ends at `fa_src_pos`, whose bases
    /// (up to a NUL) were appended to the buffer at `at` (bcftools' `mask_region`).
    fn mask_region(&mut self, at: usize, len: usize) {
        let seq = &mut self.fa_buf[at..];
        let start = self.fa_src_pos - len as i64;
        let end = self.fa_src_pos;
        let (beg32, end32) = (clamp_u32(start), clamp_u32(end));
        let last = seq.len() as i64 - 1;
        for mask in &mut self.masks {
            for r in mask.idx.overlap(&self.chr, beg32, end32) {
                let idx_start = (i64::from(r.beg) - start).max(0);
                let idx_end = (i64::from(r.end) - start).min(last);
                if idx_start > idx_end {
                    continue;
                }
                let span = &mut seq[idx_start as usize..=idx_end as usize];
                match mask.with {
                    MaskWith::Upper => span.make_ascii_uppercase(),
                    MaskWith::Lower => span.make_ascii_lowercase(),
                    MaskWith::Char(c) => span.fill(c),
                }
            }
        }
    }

    /// Whether the record is on the current sequence and within its region.
    fn in_region(&self, rec: &Record) -> bool {
        rec.rid() == self.rid && (self.fa_end_pos == 0 || rec.pos() <= self.fa_end_pos)
    }

    /// The original-sequence position just past the buffer.
    fn buffer_end(&self) -> i64 {
        self.fa_ori_pos + self.fa_buf.len() as i64 - self.fa_mod_off
    }

    /// Applies the records of the current sequence that pass `keep` (bcftools' loops over
    /// `next_vcf_line` at the end of a sequence and of the FASTA). The first record that
    /// does not is dropped, as bcftools drops it.
    fn apply_rest(
        &mut self,
        input: &mut VcfInput,
        keep: impl Fn(&Self, &Record) -> bool,
    ) -> anyhow::Result<()> {
        while self.rid >= 0 {
            let Some(rec) = input.peek()? else {
                break;
            };
            let apply = keep(self, rec);
            let mut rec = input.take();
            let result = if apply {
                self.apply_variant(&mut rec)
            } else {
                Ok(())
            };
            input.recycle(rec);
            result?;
            if !apply {
                break;
            }
        }
        Ok(())
    }

    /// Streams the FASTA, applying the variants (bcftools' `consensus`).
    fn run(&mut self, input: &mut VcfInput, mut fasta: Box<dyn BufRead>) -> anyhow::Result<()> {
        let mut line = Vec::new();
        loop {
            let at = self.fa_buf.len();
            let len = match read_fasta_line(&mut fasta, &mut self.fa_buf, &mut line)? {
                None => break,
                Some(FastaLine::Sequence(len)) => len,
                Some(FastaLine::Header) => {
                    // A new sequence: finish the previous one.
                    self.print_chain()?;
                    // Apply the cached variants and those missed because the FASTA was short.
                    self.apply_rest(input, Self::in_region)?;
                    // bcftools then calls apply_absent up to the next pending record of this
                    // sequence; none is left, so it passes 0, which does nothing.
                    self.flush_fa_buffer(FlushMode::All)?;
                    self.init_region(input, &line)?;
                    continue;
                }
            };
            self.longest_line = self.longest_line.max(len);
            self.fa_length += len as i64;
            self.fa_src_pos += len as i64;
            if !self.masks.is_empty() {
                self.mask_region(at, len);
            }

            let mut have_rec = false;
            while self.rid >= 0 {
                let Some(rec) = input.peek()? else {
                    have_rec = false;
                    break;
                };
                have_rec = true;
                if !self.in_region(rec) || self.buffer_end() <= rec.pos() {
                    // A record on another sequence, or well beyond the buffer: flush.
                    have_rec = false;
                    break;
                }
                if self.buffer_end() < rec.pos() + rec.rlen() {
                    // The buffer does not cover the record yet: read more.
                    break;
                }
                let mut rec = input.take();
                let result = self.apply_variant(&mut rec);
                input.recycle(rec);
                result?;
            }
            if !have_rec {
                let end = self.buffer_end();
                self.apply_absent(end);
                // No record left can change the buffer.
                self.flush_fa_buffer(FlushMode::Lines)?;
            }
        }
        self.apply_rest(input, |c, rec| {
            c.in_region(rec) && c.buffer_end() > rec.pos()
        })?;
        self.print_chain()?;
        self.apply_absent(POS_MAX);
        self.flush_fa_buffer(FlushMode::All)?;
        eprintln!("Applied {} variants", self.napplied);
        Ok(())
    }

    /// Flushes the outputs after [`Consensus::run`] gave `result`; like bcftools, what was
    /// produced before an error is written.
    fn finish(&mut self, result: anyhow::Result<()>) -> anyhow::Result<()> {
        let flushed = self.out.flush();
        if let Some(c) = &mut self.chain_out {
            c.flush()
                .with_context(|| format!("Close failed: {}", self.chain_name))?;
        }
        result?;
        flushed.with_context(|| format!("Close failed: {}", self.out_name))
    }
}

/// A line of the FASTA.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FastaLine {
    /// A `>` line.
    Header,
    /// A sequence line of this many bytes.
    Sequence(usize),
}

/// Reads the next FASTA line as HTSlib's `bgzf_getline` does (without the `\n` and a trailing
/// `\r`): a header, without its `>`, into `header`, or sequence appended to `seq` up to a NUL
/// (bcftools appends the line as a C string). Sequence lines are copied straight from the read
/// buffer. `None` at the end of the input.
fn read_fasta_line(
    input: &mut impl BufRead,
    seq: &mut Vec<u8>,
    header: &mut Vec<u8>,
) -> io::Result<Option<FastaLine>> {
    let start = seq.len();
    let mut kind = None;
    let mut len = 0usize;
    let mut nul = false;
    let mut last = 0u8;
    loop {
        let buf = match input.fill_buf() {
            Ok(b) => b,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        };
        let Some(&first) = buf.first() else {
            break;
        };
        let is_header = *kind.get_or_insert_with(|| {
            header.clear();
            first == b'>'
        });
        let (piece, eol) = match memchr::memchr(b'\n', buf) {
            Some(i) => (&buf[..i], Some(i)),
            None => (buf, None),
        };
        if let Some(&c) = piece.last() {
            last = c;
        }
        len += piece.len();
        if is_header {
            header.extend_from_slice(piece);
        } else if !nul {
            match memchr::memchr(0, piece) {
                Some(z) => {
                    seq.extend_from_slice(&piece[..z]);
                    nul = true;
                }
                None => seq.extend_from_slice(piece),
            }
        }
        let used = eol.map_or(buf.len(), |i| i + 1);
        input.consume(used);
        if eol.is_some() {
            break;
        }
    }
    let Some(is_header) = kind else {
        return Ok(None);
    };
    if len > 0 && last == b'\r' {
        len -= 1;
        if is_header {
            header.pop();
        } else if !nul {
            seq.pop();
        }
    }
    if is_header {
        header.remove(0);
        return Ok(Some(FastaLine::Header));
    }
    debug_assert!(seq.len() - start <= len);
    Ok(Some(FastaLine::Sequence(len)))
}

/// How much of the buffer [`Consensus::flush_fa_buffer`] writes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FlushMode {
    /// Whole lines before the protected position, or all whole lines of a long buffer but its
    /// last base (the sequence goes on).
    Lines,
    /// Everything (the end of a sequence).
    All,
}

fn clamp_u32(v: i64) -> u32 {
    v.clamp(0, i64::from(u32::MAX)) as u32
}

/// Splits a FASTA header (without `>`) into the sequence name and the 0-based `from`, `to` of a
/// `chr:from-to` header, 0 when absent (bcftools' `init_region`).
fn parse_fasta_header(line: &[u8]) -> (&[u8], i64, i64) {
    let line = line
        .iter()
        .position(|&c| c == 0)
        .map_or(line, |i| &line[..i]);
    let name_end = line
        .iter()
        .position(|&c| c.is_ascii_whitespace() || c == 0x0b || c == b':')
        .unwrap_or(line.len());
    let name = &line[..name_end];
    if name_end == line.len() {
        return (name, 0, 0);
    }
    let at = |i: usize| line.get(i).copied().unwrap_or(0);
    let ss = name_end + 1;
    let (from, se) = pasta::variant::number::strtol(line, ss);
    if se == ss || at(se) != b'-' {
        return (name, 0, 0);
    }
    let from = from - 1;
    let ss = se + 1;
    let (to, se) = pasta::variant::number::strtol(line, ss);
    if se == ss || (at(se) != 0 && !is_c_space(at(se))) {
        return (name, 0, 0);
    }
    // C keeps these in an `hts_pos_t`; negative values make no sense, so clamp them.
    (name, from.max(0), (to - 1).max(0))
}

/// C's `isspace`.
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Whether choosing an allele from the `-H` sample's GT depends on the number of values per
/// sample (bcftools' `fmt->n`, the largest ploidy of all samples at the record): `-H N` with
/// N > 1 reads past a shorter genotype's end, and `-H NpIu` tests the last value's phase.
fn uses_max_ploidy(opts: &Options) -> bool {
    if opts.allele == PICK_IUPAC {
        opts.haplotype != 0
    } else {
        !opts.output_iupac && opts.haplotype > 1
    }
}

/// Makes the reader parse only the samples consensus uses, and none in REF,ALT mode (a vtools
/// speed-up: bcftools parses every sample). Only without `-i/-e`, whose sample terms see all
/// samples, and not when the `-H` mode depends on the ploidy of the other samples (VCF
/// parsing sizes GT by the kept samples). Remaps `isample` and `iupac_samples` to the kept
/// samples.
fn read_only_used_samples(
    sr: &mut SyncedReaders,
    opts: &Options,
    isample: &mut Option<usize>,
    iupac_samples: &mut Option<Vec<usize>>,
) -> anyhow::Result<()> {
    let nsmpl = sr.header(0).sample_count();
    let chosen: Vec<usize> = match (&isample, &iupac_samples) {
        (Some(_), _) if uses_max_ploidy(opts) => return Ok(()),
        (Some(i), _) => vec![*i],
        (None, Some(v)) => v.clone(),
        (None, None) => {
            if nsmpl > 0 {
                sr.set_max_unpack(Unpack::STR | Unpack::FLT | Unpack::INFO);
            }
            return Ok(());
        }
    };
    // BCF records decode FORMAT lazily; subsetting them would only add work.
    let is_vcf = sr.file_mut(0).format() == InputFormat::Vcf;
    let samples = sr.header(0).samples();
    if !is_vcf || chosen.len() == nsmpl || chosen.iter().any(|&i| samples[i].contains(&b',')) {
        return Ok(());
    }
    let list = chosen
        .iter()
        .map(|&i| &samples[i][..])
        .collect::<Vec<_>>()
        .join(&b',');
    sr.header_mut(0)
        .set_samples(Some(&list), ListSource::Inline)
        .context("Error parsing the sample list")?;
    // `chosen` is in VCF order, as the kept samples are.
    *isample = isample.map(|_| 0);
    *iupac_samples = iupac_samples.as_ref().map(|v| (0..v.len()).collect());
    Ok(())
}

/// Opens a plain or compressed FASTA, `-` for standard input.
fn open_fasta(path: &str) -> anyhow::Result<Box<dyn BufRead>> {
    /// Read buffer size.
    const BUF: usize = 1 << 17;
    let open = || -> io::Result<Box<dyn BufRead>> {
        let file: Box<dyn io::Read> = if path == "-" {
            Box::new(io::stdin())
        } else {
            Box::new(File::open(path)?)
        };
        let mut reader = io::BufReader::with_capacity(BUF, file);
        Ok(if reader.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
            Box::new(io::BufReader::with_capacity(
                BUF,
                flate2::bufread::MultiGzDecoder::new(reader),
            ))
        } else {
            Box::new(reader)
        })
    };
    open().map_err(|_| anyhow::anyhow!("Error reading {path}"))
}

/// Chooses the samples whose genotypes are applied (bcftools' `smpl_ilist_init` calls in
/// `init_data`) and prints the note on what is applied: the `-H` sample, or the samples
/// whose genotypes give IUPAC codes; neither in REF,ALT mode.
fn select_samples(
    args: &Args,
    header: &Header,
    opts: &Options,
) -> anyhow::Result<(Option<usize>, Option<Vec<usize>>)> {
    let nsmpl = header.sample_count();
    let mut smpl: Option<(Vec<usize>, usize)> = None;
    if let Some(file) = &args.samples_file {
        let s = sample_list(header, file, ListSource::File)?;
        if s.1 == 0 {
            bail!("No matching sample found");
        }
        smpl = Some(s);
    } else if let Some(list) = args.samples.as_deref().filter(|s| *s != "-") {
        let s = sample_list(header, list, ListSource::Inline)?;
        if s.1 == 0 {
            bail!("No matching sample found");
        }
        smpl = Some(s);
    } else if args.samples.is_none() && nsmpl > 0 {
        smpl = Some(((0..nsmpl).collect(), nsmpl));
    }
    let mut isample = None;
    let mut iupac_samples = None;
    let ignoring = if nsmpl > 0 { ", ignoring samples" } else { "" };
    match smpl {
        Some((idx, n)) if opts.haplotype != 0 || opts.allele != 0 => {
            if n > 1 {
                bail!("Too many samples, only one can be used with -H; check the -s,-S options");
            }
            isample = idx.first().copied();
        }
        Some((idx, n)) => {
            if n == 1 {
                let name = idx.first().map_or(&b""[..], |&i| &header.samples()[i][..]);
                eprintln!(
                    "Note: applying IUPAC codes based on FORMAT/GT in sample {}",
                    String::from_utf8_lossy(name)
                );
            } else {
                eprintln!("Note: applying IUPAC codes based on FORMAT/GT in {n} samples");
            }
            iupac_samples = Some(idx);
        }
        None if opts.output_iupac => {
            eprintln!("Note: applying IUPAC codes based on REF,ALT{ignoring}");
        }
        None => eprintln!("Note: applying REF,ALT variants{ignoring}"),
    }
    Ok((isample, iupac_samples))
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let Some((mut cons, mut input, ref_fname)) = setup(&args, matches)? else {
        return Ok(usage());
    };
    let fasta = open_fasta(&ref_fname)?;
    let result = cons.run(&mut input, fasta);
    cons.finish(result)?;
    Ok(ExitCode::SUCCESS)
}

/// Parses the options and opens the inputs and outputs (bcftools' `init_data`); returns the
/// state, the VCF and the FASTA's name, or `None` to print the usage.
fn setup(
    args: &Args,
    matches: &ArgMatches,
) -> anyhow::Result<Option<(Consensus, VcfInput, String)>> {
    let ParsedOptions {
        opts,
        filter: filter_expr,
        masks,
        help,
    } = parse_options(args, matches)?;
    if help {
        return Ok(None);
    }
    let Some(fname) = args.inputs.first().cloned() else {
        return Ok(None);
    };
    let ref_fname = match &args.fasta_ref {
        Some(f) => f.clone(),
        None if !io::stdin().is_terminal() => "-".to_owned(),
        None => return Ok(None),
    };

    // init_data
    let mut sr = SyncedReaders::new();
    sr.set_index_requirement(IndexRequirement::Required);
    sr.set_threads(crate::output::threads(args.threads.as_deref()).unwrap_or(0));
    if let Err(e) = sr.add_reader(&fname) {
        let name = if fname == "-" {
            "standard input"
        } else {
            fname.as_str()
        };
        bail!("Failed to read from {name}: {e}");
    }
    let header = sr.header(0).clone();
    let (mut isample, mut iupac_samples) = select_samples(args, &header, &opts)?;
    let mut mask_list = Vec::with_capacity(masks.len());
    for (file, with) in masks {
        // `hts_open` reads `FILE` of `FILE##idx##INDEX`.
        let data = split_idx_name(&file).0;
        if let Err(e) = File::open(data) {
            log::error(
                "hts_open_format",
                format_args!("Failed to open file \"{data}\" : {}", strerror(&e)),
            );
            bail!("Failed to initialize mask regions");
        }
        let idx = RegIdx::from_file(Path::new(&file)).map_err(|e| {
            eprintln!("{e}");
            anyhow::anyhow!("Failed to initialize mask regions")
        })?;
        mask_list.push(Mask { idx, with });
    }
    let chain_out = match &args.chain {
        Some(c) => Some(BufWriter::new(File::create(c).map_err(|e| {
            anyhow::anyhow!("Failed to create {c}: {}", strerror(&e))
        })?)),
        None => None,
    };
    let (out, out_name): (Box<dyn Write>, String) = match &args.output {
        Some(o) => (
            Box::new(
                File::create(o)
                    .map_err(|e| anyhow::anyhow!("Failed to create {o}: {}", strerror(&e)))?,
            ),
            o.clone(),
        ),
        None => (Box::new(crate::output::stdout()), "(null)".to_owned()),
    };
    let filter = match filter_expr {
        Some((expr, logic)) => Some((filter::init(&header, &expr)?, logic)),
        None => None,
    };
    let gt_id = header
        .id(b"GT")
        .filter(|&id| header.id_exists(pasta::variant::HeaderRecordType::Format, id))
        .map(|id| id as i32);
    if filter.is_none() {
        read_only_used_samples(&mut sr, &opts, &mut isample, &mut iupac_samples)?;
    }
    let input = VcfInput {
        sr,
        filter,
        pending: None,
        pool: Vec::new(),
    };
    let cons = Consensus {
        opts,
        header,
        iupac_samples,
        isample,
        gt_id,
        fa_buf: Vec::new(),
        longest_line: 0,
        fa_ori_pos: 0,
        fa_frz_pos: 0,
        fa_mod_off: 0,
        fa_frz_mod: 0,
        absent_filled: 0,
        fa_beg_pos: 0,
        fa_end_pos: 0,
        fa_length: 0,
        fa_src_pos: 0,
        prev_base: 0,
        prev_base_pos: 0,
        prev_is_insert: false,
        rid: -1,
        chr: Vec::new(),
        masks: mask_list,
        chain: None,
        chain_id: 0,
        chain_out,
        chain_name: args.chain.clone().unwrap_or_default(),
        out: BufWriter::with_capacity(1 << 16, out),
        out_name,
        lines_out: Vec::new(),
        napplied: 0,
        warned_haplotype: false,
        warned_multibase_del: false,
        scratch: apply::Scratch::default(),
    };

    Ok(Some((cons, input, ref_fname)))
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use clap::{Args as _, Command, FromArgMatches as _};
    use pasta::{
        bgzf,
        variant::{
            index::{
                IndexFormat,
                tabix::{TabixBuilder, TabixConfig},
            },
            io::Reader,
        },
    };

    use super::*;

    /// A scratch directory, removed when dropped.
    struct Dir(PathBuf);

    impl Drop for Dir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// Writes `text` as a bgzipped VCF with a TBI index; returns its path.
    fn indexed_vcf(dir: &Path, text: &str) -> anyhow::Result<String> {
        let path = dir.join("in.vcf.gz");
        let mut w = bgzf::io::Writer::new(File::create(&path)?);
        w.write_all(text.as_bytes())?;
        w.finish()?;
        let fname = path.to_str().expect("UTF-8 path").to_owned();
        let mut reader = Reader::open(&fname, 0)?;
        let mut builder = TabixBuilder::new(TabixConfig::VCF, 0);
        let mut line = Vec::new();
        while reader.read_line(&mut line)? {
            builder.push_line(&line, reader.virtual_position())?;
        }
        builder
            .finish(reader.virtual_position())
            .index
            .save(&dir.join("in.vcf.gz.tbi"), IndexFormat::Tbi)?;
        Ok(fname)
    }

    /// Runs `consensus` with `options` on the FASTA `fasta`; returns the output FASTA.
    fn consensus(dir: &Path, vcf: &str, fasta: &str, options: &[&str]) -> anyhow::Result<String> {
        Ok(consensus_buffer(dir, vcf, fasta, options)?.0)
    }

    /// Runs `consensus` as [`consensus`] does; also returns the capacity the sequence buffer
    /// keeps at the end.
    fn consensus_buffer(
        dir: &Path,
        vcf: &str,
        fasta: &str,
        options: &[&str],
    ) -> anyhow::Result<(String, usize)> {
        let fa = dir.join("ref.fa");
        std::fs::write(&fa, fasta)?;
        let out = dir.join("out.fa");
        let (fa, out) = (
            fa.to_str().expect("UTF-8 path"),
            out.to_str().expect("UTF-8 path"),
        );
        let argv = ["consensus", "-f", fa, "-o", out];
        let matches = Args::augment_args(Command::new("consensus"))
            .try_get_matches_from(argv.iter().chain(options).chain([&vcf]))?;
        let args = Args::from_arg_matches(&matches)?;
        let (mut cons, mut input, ref_fname) =
            setup(&args, &matches)?.expect("the arguments name the inputs");
        let result = cons.run(&mut input, open_fasta(&ref_fname)?);
        cons.finish(result)?;
        Ok((std::fs::read_to_string(out)?, cons.fa_buf.capacity()))
    }

    /// A scratch directory named after `test`.
    fn scratch_dir(test: &str) -> anyhow::Result<Dir> {
        let name = format!("vtools-consensus-{test}-{}", std::process::id());
        let dir = Dir(std::env::temp_dir().join(name));
        std::fs::create_dir_all(&dir.0)?;
        Ok(dir)
    }

    /// `seq` as a FASTA entry with lines of `width`.
    fn fasta_entry(name: &str, seq: &str, width: usize) -> String {
        let lines: Vec<_> = seq
            .as_bytes()
            .chunks(width)
            .map(String::from_utf8_lossy)
            .collect();
        format!(">{name}\n{}\n", lines.join("\n"))
    }

    /// A VCF header with the contigs `c1` and `c2`, then `records`.
    fn vcf_text(records: &str) -> String {
        format!(
            "##fileformat=VCFv4.2\n##contig=<ID=c1>\n##contig=<ID=c2>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"
        )
    }

    #[test]
    fn sequences_are_written_while_no_record_is_pending() -> anyhow::Result<()> {
        // bcftools holds a sequence from its start, or from the last applied variant, to the
        // next variant or its end; vtools writes it as it reads it.
        let dir = scratch_dir("streaming")?;
        let seq = "ACGTTGCAGT".repeat(200_000);
        let fasta = fasta_entry("c1", &seq, 70);
        let last = seq.len();
        let cases = [
            (String::new(), seq.clone()),
            (
                "c1\t1\t.\tA\tG\t.\t.\t.\n".to_owned(),
                format!("G{}", &seq[1..]),
            ),
            (
                format!("c1\t{last}\t.\tT\tC\t.\t.\t.\n"),
                format!("{}C", &seq[..last - 1]),
            ),
            (
                "c1\t1000000\t.\tT\tTAAA\t.\t.\t.\n".to_owned(),
                format!("{}AAA{}", &seq[..1_000_000], &seq[1_000_000..]),
            ),
        ];
        for (records, expected) in cases {
            let vcf = indexed_vcf(&dir.0, &vcf_text(&records))?;
            let (out, capacity) = consensus_buffer(&dir.0, &vcf, &fasta, &[])?;
            assert!(
                out == fasta_entry("c1", &expected, LINE_WIDTH),
                "wrong output with {records:?}"
            );
            assert!(
                capacity <= FA_BUF_KEEP,
                "buffer of {capacity} bytes with {records:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn a_long_record_is_released_for_the_next_sequence() -> anyhow::Result<()> {
        let dir = scratch_dir("release")?;
        let seq = "ACGTTGCAGT".repeat(200_000);
        let fasta = fasta_entry("c1", &seq, 70) + &fasta_entry("c2", &seq[..100], 70);
        // The buffer holds the whole deletion before applying it.
        let del = &seq[9..1_500_000];
        let records = format!("c1\t10\t.\t{del}\t{}\t.\t.\t.\n", &del[..1]);
        let vcf = indexed_vcf(&dir.0, &vcf_text(&records))?;
        let (out, capacity) = consensus_buffer(&dir.0, &vcf, &fasta, &[])?;
        let expected = fasta_entry(
            "c1",
            &format!("{}{}", &seq[..10], &seq[1_500_000..]),
            LINE_WIDTH,
        ) + &fasta_entry("c2", &seq[..100], LINE_WIDTH);
        assert!(out == expected, "wrong output");
        assert!(capacity <= FA_BUF_KEEP, "buffer of {capacity} bytes");
        Ok(())
    }

    #[test]
    fn long_lines_keep_the_buffer_for_the_next_sequence() -> anyhow::Result<()> {
        // Releasing the room an unwrapped FASTA's lines need only makes the next sequence
        // fault it in again.
        let dir = scratch_dir("keep")?;
        let seq = "ACGTTGCAGT".repeat(200_000);
        let fasta = format!(">c1\n{seq}\n>c2\n{}\n", &seq[..100]);
        let vcf = indexed_vcf(&dir.0, &vcf_text(""))?;
        let (out, capacity) = consensus_buffer(&dir.0, &vcf, &fasta, &[])?;
        let expected =
            fasta_entry("c1", &seq, LINE_WIDTH) + &fasta_entry("c2", &seq[..100], LINE_WIDTH);
        assert!(out == expected, "wrong output");
        assert!(capacity >= seq.len(), "buffer of {capacity} bytes");
        Ok(())
    }

    #[test]
    fn absent_fills_after_lines_written_past_an_insertion() -> anyhow::Result<()> {
        // Once lines past a long insertion are written, the buffer start no longer bounds the
        // positions --absent may fill.
        let dir = scratch_dir("absent")?;
        let seq = "ACGTTGCAGT".repeat(20_000);
        let insert = "CAGT".repeat(20_000);
        let records = format!("c1\t100\t.\tT\tT{insert}\t.\t.\t.\nc1\t150000\t.\tT\tA\t.\t.\t.\n");
        let vcf = indexed_vcf(&dir.0, &vcf_text(&records))?;
        let fasta = fasta_entry("c1", &seq, 60);
        let out = consensus(&dir.0, &vcf, &fasta, &["-a", "N"])?;
        let n = |len: usize| "N".repeat(len);
        let expected = format!(
            "{}T{insert}{}A{}",
            n(99),
            n(150_000 - 101),
            n(seq.len() - 150_000)
        );
        assert!(
            out == fasta_entry("c1", &expected, LINE_WIDTH),
            "wrong output"
        );
        Ok(())
    }

    #[test]
    fn missing_block_is_clipped_to_the_region() -> anyhow::Result<()> {
        // A ./. gVCF block at 11-150 masks only the part of it in each `>chr:from-to` region.
        // bcftools copies the block's REF from before the buffer of a region that starts inside
        // it, and past the buffer of one that ends inside it, and crashes.
        let dir = Dir(std::env::temp_dir().join(format!(
            "vtools-consensus-missing-block-{}",
            std::process::id()
        )));
        std::fs::create_dir_all(&dir.0)?;
        let seq = "ACGTTGCAGT".repeat(20);
        let vcf = indexed_vcf(
            &dir.0,
            &format!(
                "##fileformat=VCFv4.2\n##contig=<ID=c1,length=200>\n\
                 ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                 c1\t11\t.\t{}\t<NON_REF>\t.\t.\tEND=150\tGT\t./.\n",
                &seq[10..11]
            ),
        )?;
        let region = |from: usize, to: usize| format!(">c1:{from}-{to}\n{}\n", &seq[from - 1..to]);
        let n = |len: usize| "N".repeat(len);
        let options = ["-s", "S1", "-H", "1", "-M", "N"];

        let whole = consensus(&dir.0, &vcf, &format!(">c1\n{seq}\n"), &options)?;
        let masked = format!("{}{}{}", &seq[..10], n(140), &seq[150..]);
        let lines: Vec<_> = masked
            .as_bytes()
            .chunks(LINE_WIDTH)
            .map(String::from_utf8_lossy)
            .collect();
        assert_eq!(whole, format!(">c1\n{}\n", lines.join("\n")));

        let fasta = region(1, 20) + &region(21, 40) + &region(141, 180);
        let regions = consensus(&dir.0, &vcf, &fasta, &options)?;
        let expected = format!(
            ">c1:1-20\n{}{}\n>c1:21-40\n{}\n>c1:141-180\n{}{}\n",
            &seq[..10],
            n(10),
            n(20),
            n(10),
            &seq[150..180]
        );
        assert_eq!(regions, expected);
        Ok(())
    }

    #[test]
    fn fasta_headers() {
        assert_eq!(parse_fasta_header(b"chr1"), (&b"chr1"[..], 0, 0));
        assert_eq!(parse_fasta_header(b"chr1 desc"), (&b"chr1"[..], 0, 0));
        assert_eq!(parse_fasta_header(b"22:11-20"), (&b"22"[..], 10, 19));
        assert_eq!(parse_fasta_header(b"22:11-20 x"), (&b"22"[..], 10, 19));
        assert_eq!(parse_fasta_header(b"22:11-20x"), (&b"22"[..], 0, 0));
        assert_eq!(parse_fasta_header(b"22:11"), (&b"22"[..], 0, 0));
        // A description of the form "from-to" counts as a range, as in bcftools.
        assert_eq!(parse_fasta_header(b"1 5-9"), (&b"1"[..], 4, 8));
    }

    #[test]
    fn chain_gaps() {
        let mut c = Chain::new(0);
        c.push_gap(10, 0, 10, 2);
        c.push_gap(10, 3, 12, 0);
        c.push_gap(20, 1, 19, 1);
        let block = |len, ref_gap, alt_gap| ChainBlock {
            len,
            ref_gap,
            alt_gap,
        };
        assert_eq!(c.blocks, [block(10, 3, 2), block(7, 1, 1)]);
        assert_eq!((c.ref_last_block_ori, c.alt_last_block_ori), (21, 20));
    }

    #[test]
    fn haplotype_values() {
        let mut o = Options::default();
        parse_haplotype(&mut o, "3pIu").unwrap();
        assert_eq!((o.haplotype, o.allele), (3, PICK_IUPAC));
        let mut o = Options::default();
        parse_haplotype(&mut o, "la").unwrap();
        assert_eq!(o.allele, PICK_LONG | PICK_ALT);
        assert!(parse_haplotype(&mut Options::default(), "0").is_err());
        assert!(parse_haplotype(&mut Options::default(), "2x").is_err());
    }
}
