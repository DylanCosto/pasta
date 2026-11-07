//! `ntools view`: convert, filter, count, and query alignments (samtools view).

use std::{
    collections::HashSet,
    fs,
    io::{self, Write},
    mem,
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use bstr::BString;
use pasta::bam;
use pasta::cram;
use pasta::fasta;
use pasta::sam::{
    self as sam,
    alignment::{io::Write as _, record::Flags},
    header::record::value::{Map, map::ReferenceSequence},
};

use crate::{
    alignment::{
        AlignmentReader, CramQueries, Queries, RawEncoder, ReadOptions, Rec, encode_cram_record,
    },
    bam_raw,
    bam_sam::{self, Names},
    depth::format_g,
    expr::Expression,
    filter::{self, Filter, ReadGroupFilter, ReadGroupPresence, Subsample, TagFilter},
    flags,
    global_args::{IndexRequest, InputFormatOptions, OutputFormatOptions, Verbosity},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    parallel_cram::{self, Source},
    program,
    regions::{RegionSet, record_interval},
    source::{self, Format},
};

/// The region name samtools uses for unplaced, unmapped records.
const UNMAPPED_REGION: &str = "*";

/// The notes of samtools' long help (`-?`), for the options ntools supports.
const NOTES: &str = "\
Notes:

1. The input format (SAM, BAM, or CRAM) is detected automatically.

2. The file given with -t lists a reference name and length on each line, separated by a tab
   (the first two columns of a .fai file). It is used for SAM input without @SQ lines.

3. SAM->BAM conversion:  ntools view -b -o out.bam in.sam

4. BAM->SAM conversion:  ntools view -h in.bam

5. A region is `chr1', `chr2:1,000', or `chr3:1000-2,000'. Regions need a sorted, indexed BAM
   or CRAM file.

6. -M/--use-index computes overlaps with -L BED regions and region arguments through the index.
   This is faster, omits duplicates, and outputs the reads in file order.

7. -e takes a filter expression, e.g. 'mapq >= 30 && [NM] < 5 && qname =~ \"^read\"'. Fields:
   qname, flag (and flag.paired, flag.proper_pair, flag.unmap, flag.munmap, flag.reverse,
   flag.mreverse, flag.read1, flag.read2, flag.secondary, flag.qcfail, flag.dup,
   flag.supplementary), rname, refid, pos, endpos, mapq, cigar, ncigar, rnext, mrname, mrefid,
   pnext, mpos, tlen, seq, qual, qlen, rlen, sclen, hclen, library, and tags as [XX].
   Functions: length, min, max, avg, exists, default, exp, log, pow, sqrt. Operators, by
   precedence: ! ~ (unary + -), * / %, + -, &, ^, |, < <= > >=, == != =~ !~, && ||. Regular
   expressions are POSIX extended.
";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, after_long_help = NOTES)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin) [default: stdin, unless it is a terminal].
    input: Option<String>,

    /// Regions to query (`chr1`, `chr1:100-200`, or `*` for unplaced unmapped reads); requires an
    /// index.
    regions: Vec<String>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Output format.
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Output BAM.
    #[arg(short = 'b', long = "bam")]
    bam: bool,

    /// Output CRAM (requires -T).
    #[arg(short = 'C', long = "cram")]
    cram: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    uncompressed: bool,

    /// Use fast BAM compression (level 1).
    #[arg(short = '1', long = "fast")]
    fast: bool,

    #[command(flatten)]
    output_format_options: OutputFormatOptions,

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    /// Include the header in SAM output.
    #[arg(short = 'h', long = "with-header", overrides_with = "no_header")]
    with_header: bool,

    /// Print the header only.
    #[arg(short = 'H', long = "header-only", overrides_with = "no_header")]
    header_only: bool,

    /// Print no header in SAM output (cancels an earlier -h or -H).
    #[arg(long = "no-header", overrides_with_all = ["with_header", "header_only"])]
    no_header: bool,

    /// Write the numbers of records processed, accepted, and rejected to this file, as JSON.
    #[arg(long = "save-counts", value_name = "FILE")]
    save_counts: Option<PathBuf>,

    /// Print only the count of matching records.
    #[arg(short = 'c', long = "count")]
    count: bool,

    /// Write records that fail the filters to this file.
    #[arg(short = 'U', long = "unoutput", visible_alias = "output-unselected")]
    unselected_output: Option<PathBuf>,

    /// Set the unmapped flag on reads that fail the filters (clearing their CIGAR, mapping
    /// quality, and template length) and write them to the output.
    #[arg(short = 'p', long = "unmap")]
    unmap: bool,

    /// Retrieve complete pairs: also output the mates, wherever they are, of paired reads in
    /// the regions or -L BED file. Requires an index.
    #[arg(short = 'P', long = "fetch-pairs")]
    fetch_pairs: bool,

    /// Reference names and lengths (the first two columns of a `.fai` file), for SAM input
    /// without `@SQ` lines.
    #[arg(short = 't', long = "fai-reference")]
    fai_reference: Option<PathBuf>,

    /// Only include reads overlapping this BED file.
    #[arg(short = 'L', long = "target-file", visible_alias = "targets-file")]
    bed: Option<PathBuf>,

    /// Only include reads overlapping this BED file, read through the index (`-L FILE -M`).
    #[arg(long = "region-file", visible_alias = "regions-file")]
    regions_file: Option<PathBuf>,

    /// Merge overlapping regions and output each record once, in file order.
    #[arg(short = 'M', long = "use-index")]
    merge_regions: bool,

    /// The input is followed by its index file, then the regions.
    #[arg(
        short = 'X',
        long = "customized-index",
        visible_alias = "customised-index"
    )]
    customized_index: bool,

    /// The index file given after the input with `-X`.
    #[arg(skip)]
    index: Option<PathBuf>,

    // The flag options may be repeated, and their values are ORed, as in samtools.
    /// Only include reads with all of these flags set.
    #[arg(short = 'f', long = "require-flags", value_parser = flags::parse)]
    require_flags: Vec<Flags>,

    /// Exclude reads with any of these flags set.
    #[arg(short = 'F', long = "excl-flags", visible_alias = "exclude-flags", value_parser = flags::parse)]
    exclude_flags: Vec<Flags>,

    /// Exclude reads with all of these flags set.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_all_flags: Vec<Flags>,

    /// Only include reads with any of these flags set.
    #[arg(long = "rf", visible_aliases = ["incl-flags", "include-flags"], value_parser = flags::parse)]
    require_any_flags: Vec<Flags>,

    /// Set these flags on the records written.
    #[arg(long = "add-flags", value_parser = flags::parse)]
    add_flags: Vec<Flags>,

    /// Clear these flags on the records written.
    #[arg(long = "remove-flags", value_parser = flags::parse)]
    remove_flags: Vec<Flags>,

    /// Only include reads in this read group, or in no read group (may be repeated).
    #[arg(short = 'r', long = "read-group", visible_alias = "readgroup")]
    read_groups: Vec<String>,

    /// Only include reads in a read group listed in this file (whitespace-separated IDs), or in
    /// no read group; with `^FILE`, exclude the listed read groups instead.
    #[arg(
        short = 'R',
        long = "read-group-file",
        visible_alias = "readgroup-file",
        value_name = "[^]FILE"
    )]
    read_group_files: Vec<String>,

    /// Exclude reads without a read group.
    #[arg(
        short = 'n',
        long = "exclude-no-read-group",
        visible_aliases = ["excl-no-read-group", "excl-no-readgroup", "exclude-no-readgroup"]
    )]
    exclude_no_read_group: bool,

    /// Only include reads whose read group belongs to this library (`LB`).
    #[arg(short = 'l', long = "library")]
    library: Option<String>,

    /// Exclude reads with mapping quality below this value.
    #[arg(
        short = 'q',
        long = "min-MQ",
        visible_aliases = ["min-mq", "min-mapq"],
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<u8>
    )]
    min_mapping_quality: u8,

    /// Only include reads covering at least this many query bases (by CIGAR, with hard clips).
    #[arg(
        short = 'm',
        long = "min-qlen",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    min_query_length: i64,

    /// Only include reads matching this filter expression (see the notes in the long help).
    #[arg(
        short = 'e',
        long = "expr",
        visible_alias = "expression",
        allow_hyphen_values = true
    )]
    expression: Option<String>,

    /// Subsample: keep a FRAC fraction of templates, seeded by INT or, with `auto`, by the
    /// header (`INT.FRAC` or `auto.FRAC`; --subsample FRAC --subsample-seed INT).
    #[arg(
        short = 's',
        value_name = "INT.FRAC",
        value_parser = SeededFraction::parse,
        allow_hyphen_values = true
    )]
    seeded_subsample: Option<SeededFraction>,

    /// Keep only this fraction of templates.
    #[arg(long = "subsample", value_name = "FLOAT", value_parser = parse_fraction)]
    subsample: Option<f64>,

    /// The --subsample seed, or `auto` to derive it from the header.
    #[arg(
        long = "subsample-seed",
        value_name = "INT|auto",
        value_parser = SubsampleSeed::parse,
        allow_hyphen_values = true
    )]
    subsample_seed: Option<SubsampleSeed>,

    /// The subsampling options given, in command-line order (see [`Args::record_order`]).
    #[arg(skip)]
    subsample_order: Vec<SubsampleOption>,

    /// Only include reads whose names are listed in this file (may be repeated).
    #[arg(short = 'N', long = "qname-file", visible_alias = "QNAME-file")]
    qname_files: Vec<PathBuf>,

    /// Only include reads with this tag, and optionally this value (`TAG[:VALUE]`).
    #[arg(short = 'd', long = "tag")]
    tags: Vec<String>,

    /// Only include reads with this tag whose value is listed in FILE (`TAG:FILE`).
    #[arg(short = 'D', long = "tag-file", value_name = "TAG:FILE")]
    tag_files: Vec<String>,

    /// Remove these tags from output records (comma-separated; may be repeated); `^TAGS`
    /// keeps only those tags instead.
    #[arg(short = 'x', long = "remove-tag", value_name = "[^]TAGS")]
    remove_tags: Vec<String>,

    /// Keep only these tags in output records (comma-separated; may be repeated); overrides -x.
    #[arg(long = "keep-tag", value_name = "TAGS")]
    keep_tags: Vec<String>,

    /// Collapse the backward (`B`) CIGAR operation.
    #[arg(short = 'B', long = "remove-B")]
    remove_backward_ops: bool,

    /// Index the output file (CSI for BAM, CRAI for CRAM).
    #[arg(long = "write-index")]
    write_index: bool,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decoding, formatting, and compression.
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Ignored (the input format is detected automatically).
    #[arg(short = 'S')]
    sam_input: bool,

    #[command(flatten)]
    verbosity: Verbosity,

    /// Print long help, including notes on regions and filter expressions.
    #[arg(short = '?', long = "help", action = clap::ArgAction::HelpLong)]
    help: Option<bool>,
}

/// How records are reached.
enum Traversal {
    /// Every record, in file order.
    All,
    /// samtools regions, one indexed query each.
    Regions(Vec<String>),
    /// Merged intervals through the index, each record once.
    Set(RegionSet),
    /// The records of these intervals and the mates of their paired reads (`-P`).
    Pairs(RegionSet),
}

/// A subsampling option, which samtools applies in command-line order.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SubsampleOption {
    /// `-s INT.FRAC`.
    Seeded,
    /// `--subsample FRAC`.
    Fraction,
    /// `--subsample-seed INT`.
    Seed,
}

impl Args {
    /// Records the order of options that samtools applies in command-line order, from the parsed
    /// arguments of this command.
    pub fn record_order(&mut self, matches: &clap::ArgMatches) {
        let mut options: Vec<(usize, SubsampleOption)> = [
            ("seeded_subsample", SubsampleOption::Seeded),
            ("subsample", SubsampleOption::Fraction),
            ("subsample_seed", SubsampleOption::Seed),
        ]
        .into_iter()
        .filter_map(|(id, option)| Some((matches.indices_of(id)?.max()?, option)))
        .collect();

        options.sort_unstable_by_key(|&(i, _)| i);
        self.subsample_order = options.into_iter().map(|(_, option)| option).collect();
    }
}

pub fn run(mut args: Args) -> anyhow::Result<()> {
    let input = source::input_or_stdin(args.input.take())?;

    // `--region-file FILE` is `-L FILE -M`.
    if let Some(regions_file) = args.regions_file.take() {
        args.bed = Some(regions_file);
        args.merge_regions = true;
    }

    if args.no_header {
        args.with_header = false;
        args.header_only = false;
    }

    if args.customized_index {
        if args.regions.is_empty() {
            bail!("Incorrect number of arguments for -X option. Aborting.");
        }

        args.index = Some(PathBuf::from(args.regions.remove(0)));
    }

    if args.unmap && args.unselected_output.is_some() {
        bail!("Options --unoutput and --unmap are mutually exclusive.");
    }

    if args.count && args.fetch_pairs {
        bail!("The options -P and -c cannot be combined");
    }

    let selected = output_kind(&args, args.output.as_deref())?;
    let format = selected.format;
    let reference = reference(&args);

    let mut read_options = ReadOptions {
        threads: args.threads,
        reference: reference.clone(),
    };
    args.input_format_options.configure(&mut read_options)?;

    let (mut reader, path) = AlignmentReader::open(&input, &read_options)?;
    args.input_format_options.apply(&mut reader)?;
    let mut header = reader.read_header()?;

    if reader.format() == Format::Sam
        && header.reference_sequences().is_empty()
        && let Some(fai) = fai_path(&args)
    {
        add_fai_reference_sequences(&mut header, &fai)?;
    }

    let (traversal, mut filter) = plan(&args, &header, &reader, path.as_deref())?;

    // Records are decoded with the input header (CRAM stores read groups as indexes into its
    // `@RG` lines); `header` becomes the output header, which `-r`/`-R` prune.
    let input_header = header.clone();
    configure_filter(&mut filter, &args, &mut header)?;
    let filter = filter;

    let output_edits = OutputEdits::new(&args)?;
    let edits = Edits::new(&args);

    if args.count {
        let (fields, tags) = count_fields(&filter, &traversal, edits.as_ref());
        reader.set_required_fields(fields, tags);
    }

    // With one thread, SAM text is written faster serially, formatted on the writer's thread
    // while the next records decode; everything else is faster decoded on a worker.
    let min_threads = if format.is_sam() && !args.count { 2 } else { 1 };

    let parallel = args.threads >= min_threads
        && output_edits.is_none()
        && args.unselected_output.is_none()
        && args.save_counts.is_none()
        && edits.is_none();

    // view chooses for itself when CRAM is decoded on workers.
    if !parallel {
        reader.set_cram_decoding_threads(0);
    }

    // Indexed CRAM queries are decoded on worker threads too, except `-P`'s and the unmapped
    // region's. `-H` queries nothing, so needs no index.
    let cram_queries = match (&traversal, path.as_deref(), &reader) {
        (Traversal::Regions(_) | Traversal::Set(_), Some(path), AlignmentReader::Cram(..))
            if parallel && !args.header_only =>
        {
            CramQueries::new(
                path,
                args.index.as_deref(),
                &input_header,
                &queries(&traversal),
            )?
        }
        _ => None,
    };

    let source = match (&traversal, &cram_queries) {
        _ if !parallel => None,
        (Traversal::All, _) => Some(Source::All),
        (_, Some(queries)) => Some(Source::Queries(queries)),
        (_, None) => None,
    };

    let decoding = source.map_or(Decoding::Serial, |source| {
        Decoding::Parallel(Workers {
            source,
            threads: args.threads,
        })
    });

    if args.count {
        let mut selector = Selector::new(&filter, edits);
        let counts = count(
            &mut reader,
            path.as_deref(),
            &args,
            &input_header,
            &traversal,
            &mut selector,
            decoding,
        )?;
        write_count(args.output.as_deref(), counts.accepted)?;

        if let Some(dst) = &args.save_counts {
            write_counts(dst, counts)?;
        }

        return Ok(());
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = write_options(&args, format, reference)?;

    // The `-U` output takes its own format and header policy, as samtools resolves them.
    let unselected = args
        .unselected_output
        .as_deref()
        .map(|dst| output_kind(&args, Some(dst)).map(|kind| (dst, kind)))
        .transpose()?;

    // As samtools does, both outputs are checked for indexing before anything is written.
    let mut index_requests = Vec::new();

    if args.write_index {
        index_requests.push(IndexRequest::for_output(args.output.as_deref(), format)?);

        if let Some((dst, kind)) = unselected {
            index_requests.push(IndexRequest::for_output(Some(dst), kind.format)?);
        }
    }

    let mut writer = AlignmentWriter::create(args.output.as_deref(), &write_options, &header)?;

    if selected.writes_header {
        writer.write_header(&header)?;
    }

    let mut rejected = match unselected {
        Some((dst, kind)) => {
            let options = WriteOptions {
                format: kind.format,
                ..write_options.clone()
            };

            let mut unselected_writer = AlignmentWriter::create(Some(dst), &options, &header)?;

            if kind.writes_header {
                unselected_writer.write_header(&header)?;
            }

            Rejected::Written(Box::new(unselected_writer))
        }
        None if args.unmap => Rejected::Unmapped,
        None => Rejected::Dropped,
    };

    let mut counts = Counts::default();

    if !args.header_only {
        let parallel_kind = writer.preformatted_kind();

        if let (
            AlignmentReader::Cram(cram_reader, repository, ..),
            Decoding::Parallel(workers),
            Some(kind),
        ) = (&mut reader, decoding, parallel_kind)
        {
            write_cram_in_parallel(
                cram_reader,
                repository,
                &input_header,
                &filter,
                workers,
                kind,
                &mut writer,
            )?;
        } else {
            let mut sink = Sink {
                selector: Selector::new(&filter, edits),
                output_header: &header,
                writer: &mut writer,
                rejected: &mut rejected,
                output_edits: output_edits.as_ref(),
                encoder: RawEncoder::new(),
                edited: Vec::new(),
                counts: Counts::default(),
            };

            if let Traversal::Pairs(regions) = &traversal {
                fetch_pairs(
                    &mut reader,
                    path.as_deref(),
                    args.index.as_deref(),
                    &input_header,
                    regions,
                    &mut sink,
                )?;
            } else {
                visit(
                    &mut reader,
                    path.as_deref(),
                    args.index.as_deref(),
                    &input_header,
                    &traversal,
                    |rec| sink.emit(rec, &input_header),
                )?;
            }

            counts = sink.counts;
        }
    }

    writer.finish(&header)?;

    if let Rejected::Written(unselected) = rejected {
        unselected.finish(&header)?;
    }

    for request in index_requests {
        request.write(args.threads)?;
    }

    if let Some(dst) = &args.save_counts {
        write_counts(dst, counts)?;
    }

    Ok(())
}

/// The numbers of records the filters saw and accepted (`--save-counts`).
#[derive(Clone, Copy, Debug, Default)]
struct Counts {
    processed: u64,
    accepted: u64,
}

/// Writes the count of `-c` to `dst` (as given, like samtools), or to stdout when `dst` is `None`
/// or `-`.
fn write_count(dst: Option<&Path>, count: u64) -> anyhow::Result<()> {
    match dst.filter(|path| *path != Path::new("-")) {
        Some(path) => fs::write(path, format!("{count}\n"))
            .with_context(|| format!("can't create \"{}\"", path.display())),
        None => Ok(writeln!(io::stdout().lock(), "{count}")?),
    }
}

/// Writes `--save-counts` JSON, as samtools formats it.
fn write_counts(dst: &Path, counts: Counts) -> anyhow::Result<()> {
    let text = format!(
        "{{\n    \"records_processed\" : {},\n    \"records_filter_accepted\" : {},\n    \"records_filter_rejected\" : {}\n}}\n",
        counts.processed,
        counts.accepted,
        counts.processed - counts.accepted
    );

    fs::write(dst, text).with_context(|| format!("failed to open \"{}\"", dst.display()))
}

/// Sets the record filter's options (`-f`, `-F`, `-G`, `-q`, `-m`, `-r`/`-R`, `-n`, `-l`, `-s`,
/// `-d`/`-D`, `-N`, `-e`). As in samtools, `-r`/`-R` also remove the `@RG` lines of read groups
/// whose records are dropped.
fn configure_filter(
    filter: &mut Filter,
    args: &Args,
    header: &mut sam::Header,
) -> anyhow::Result<()> {
    filter.read_group_ids = ReadGroupFilter::from_args(&args.read_groups, &args.read_group_files)?;

    if let Some(read_groups) = &filter.read_group_ids {
        read_groups.retain_header_lines(header);
    }

    filter.read_group_presence = if args.exclude_no_read_group {
        ReadGroupPresence::Required
    } else {
        ReadGroupPresence::Optional
    };

    filter.require = flags::union(&args.require_flags).unwrap_or_else(Flags::empty);
    filter.exclude = flags::union(&args.exclude_flags).unwrap_or_else(Flags::empty);
    filter.exclude_all = flags::union(&args.exclude_all_flags);
    filter.require_any = flags::union(&args.require_any_flags);
    filter.min_mapping_quality = args.min_mapping_quality;
    filter.min_query_length = args.min_query_length;
    filter.read_groups = args
        .library
        .as_deref()
        .map(|library| filter::library_read_groups(header, library));
    filter.subsample = subsample(args, header)?;
    filter.tag = TagFilter::from_args(&args.tags, &args.tag_files)?;
    filter.names = if args.qname_files.is_empty() {
        None
    } else {
        let mut names = HashSet::new();

        for path in &args.qname_files {
            names.extend(filter::read_names(path)?);
        }

        Some(names)
    };
    filter.expression = args
        .expression
        .as_deref()
        .map(Expression::parse)
        .transpose()?;

    Ok(())
}

/// The subsample (`-s`, `--subsample`, `--subsample-seed`). As samtools does, its fraction and
/// seed are recorded in a `@CO` header line; an `auto` seed is the FNV-1a hash of the header
/// text, before the line (and the `@PG` line) is added.
fn subsample(args: &Args, header: &mut sam::Header) -> anyhow::Result<Option<Subsample>> {
    let (mut seed, mut fraction) = (SubsampleSeed::Auto, 0.0);

    // Later options win, as in samtools.
    let order = if args.subsample_order.is_empty() {
        &[
            SubsampleOption::Seeded,
            SubsampleOption::Fraction,
            SubsampleOption::Seed,
        ][..]
    } else {
        &args.subsample_order[..]
    };

    for option in order {
        match option {
            SubsampleOption::Seeded => {
                if let Some(SeededFraction {
                    seed: s,
                    fraction: f,
                }) = args.seeded_subsample
                {
                    (seed, fraction) = (s, f);
                }
            }
            SubsampleOption::Fraction => {
                if let Some(f) = args.subsample {
                    fraction = f;
                }
            }
            SubsampleOption::Seed => {
                if let Some(s) = args.subsample_seed {
                    seed = s;
                }
            }
        }
    }

    if fraction <= 0.0 {
        return Ok(None);
    }

    let seed = match seed {
        SubsampleSeed::Value(n) => n,
        SubsampleSeed::Auto => {
            let mut writer = sam::io::Writer::new(Vec::new());
            writer.write_header(header)?;
            fnv1a(&writer.into_inner())
        }
    };

    header.add_comment(format!(
        "Sub-sampled fraction={} seed={seed}",
        format_g(fraction, 6)
    ));

    Ok(Some(Subsample::new(seed, fraction)))
}

/// khash's `__ac_FNV1a_hash_string`.
fn fnv1a(s: &[u8]) -> u32 {
    s.iter().fold(2_166_136_261, |h: u32, &b| {
        (h ^ u32::from(b)).wrapping_mul(16_777_619)
    })
}

/// The seed of a subsample.
#[derive(Clone, Copy, Debug)]
enum SubsampleSeed {
    Value(u32),
    /// Derived from the header.
    Auto,
}

impl SubsampleSeed {
    /// Parses `auto` or an integer (`strtol` with base 0, so `0x` prefixes work), which `srand`
    /// truncates to 32 bits.
    fn parse(s: &str) -> anyhow::Result<Self> {
        if s.eq_ignore_ascii_case("auto") {
            return Ok(Self::Auto);
        }

        let n = parse_c_long(s).with_context(|| format!("Incorrect seed argument \"{s}\""))?;

        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "samtools stores the seed in a uint32_t"
        )]
        Ok(Self::Value(n as u32))
    }
}

/// `strtol(s, &end, 0)` that must consume the whole string.
fn parse_c_long(s: &str) -> Option<i64> {
    let (negative, digits) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };

    let n = if let Some(hex) = digits
        .strip_prefix("0x")
        .or_else(|| digits.strip_prefix("0X"))
    {
        i64::from_str_radix(hex, 16).ok()?
    } else if digits.len() > 1 && digits.starts_with('0') {
        i64::from_str_radix(&digits[1..], 8).ok()?
    } else {
        digits.parse().ok()?
    };

    Some(if negative { -n } else { n })
}

/// `-s INT.FRAC` or `-s auto.FRAC`.
#[derive(Clone, Copy, Debug)]
struct SeededFraction {
    seed: SubsampleSeed,
    fraction: f64,
}

impl SeededFraction {
    fn parse(s: &str) -> anyhow::Result<Self> {
        let invalid = || format!("Incorrect sampling argument \"{s}\"");

        let (seed, fraction) = if s.len() > 4 && s[..4].eq_ignore_ascii_case("auto") {
            (SubsampleSeed::Auto, &s[4..])
        } else {
            let i = s.find('.').with_context(invalid)?;
            let raw_seed = &s[..i];

            let seed: i64 = if raw_seed.is_empty() || raw_seed == "-" || raw_seed == "+" {
                0
            } else {
                raw_seed.parse().with_context(invalid)?
            };

            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "samtools stores the seed in a uint32_t"
            )]
            (SubsampleSeed::Value(seed as u32), &s[i..])
        };

        if !fraction.starts_with('.') {
            bail!(invalid());
        }

        let fraction = fraction.parse().with_context(invalid)?;
        Ok(Self { seed, fraction })
    }
}

/// `--subsample FLOAT`, between 0 and 1.
fn parse_fraction(s: &str) -> anyhow::Result<f64> {
    s.parse()
        .ok()
        .filter(|f| (0.0..=1.0).contains(f))
        .with_context(|| format!("Incorrect sampling argument \"{s}\""))
}

/// The CRAM reference: `-T`, or else the FASTA of a `-t` index (`ref.fa` for `ref.fa.fai`), as
/// HTSlib loads it.
fn reference(args: &Args) -> Option<PathBuf> {
    args.reference.clone().or_else(|| {
        let fai = args.fai_reference.as_deref()?;
        let fasta = fai.to_str()?.strip_suffix(".fai")?;
        Path::new(fasta).is_file().then(|| PathBuf::from(fasta))
    })
}

/// The reference index giving `@SQ` lines to SAM input without them: `-t`, or the `.fai` of the
/// `-T` FASTA when it exists.
fn fai_path(args: &Args) -> Option<PathBuf> {
    args.fai_reference.clone().or_else(|| {
        let mut fai = args.reference.clone()?.into_os_string();
        fai.push(".fai");
        let fai = PathBuf::from(fai);
        fai.is_file().then_some(fai)
    })
}

/// Adds an `@SQ` line for each line of a reference index (`NAME<TAB>LENGTH...`), as HTSlib does
/// when reading SAM without `@SQ` lines. Lines without a tab are skipped.
fn add_fai_reference_sequences(header: &mut sam::Header, path: &Path) -> anyhow::Result<()> {
    let text = fs::read(path).with_context(|| format!("reading {}", path.display()))?;

    for line in text.split(|&b| b == b'\n') {
        let Some(tab) = memchr::memchr(b'\t', line) else {
            continue;
        };

        let (name, rest) = (&line[..tab], &line[tab + 1..]);
        let digits: Vec<u8> = rest
            .iter()
            .copied()
            .skip_while(u8::is_ascii_whitespace)
            .take_while(u8::is_ascii_digit)
            .collect();

        let length = std::str::from_utf8(&digits)
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .and_then(NonZeroUsize::new)
            .with_context(|| {
                format!(
                    "{}: invalid length for {}",
                    path.display(),
                    String::from_utf8_lossy(name)
                )
            })?;

        header
            .reference_sequences_mut()
            .insert(BString::from(name), Map::<ReferenceSequence>::new(length));
    }

    Ok(())
}

/// The output options, with the compression level from `-u`, `-1`, or `--output-fmt-option`.
fn write_options(
    args: &Args,
    format: OutputFormat,
    reference: Option<PathBuf>,
) -> anyhow::Result<WriteOptions> {
    // An explicit `--output-fmt-option level=N` wins over `-u` and `-1`, as HTSlib applies it
    // after opening the output.
    let level = match args.output_format_options.explicit_level()? {
        Some(level) => level,
        None if args.uncompressed => 0,
        None if args.fast => 1,
        None => output::default_level(),
    };

    Ok(WriteOptions {
        format,
        level: output::parse_level(level)?,
        threads: args.threads,
        reference,
        cram: args.output_format_options.cram()?,
    })
}

/// Chooses the traversal and any region filter: `-M` merges the regions given; `-L` restricts to
/// a BED file, through the index when one exists and no regions are given; `-P` reads the
/// regions (intersected with `-L`) through the index, then their mates.
fn plan(
    args: &Args,
    header: &sam::Header,
    reader: &AlignmentReader,
    path: Option<&Path>,
) -> anyhow::Result<(Traversal, Filter)> {
    let mut filter = Filter::default();

    let bed = match &args.bed {
        Some(bed) => {
            let mut set = RegionSet::new(header);
            set.add_bed(header, bed)?;
            set.merge();
            Some(set)
        }
        None => None,
    };

    let regions = known_regions(header, &args.regions);

    let region_set = || -> anyhow::Result<RegionSet> {
        if regions.iter().any(|r| r == UNMAPPED_REGION) {
            bail!("-M and -P cannot be combined with the unmapped region *");
        }

        let mut set = RegionSet::new(header);
        for region in &regions {
            set.add_region(header, region)?;
        }
        set.merge();
        Ok(set)
    };

    if args.fetch_pairs {
        let set = match (args.regions.is_empty(), bed) {
            (false, Some(bed)) => region_set()?.intersection(&bed),
            (false, None) => region_set()?,
            (true, Some(bed)) => bed,
            (true, None) => bail!("No regions or BED file have been provided. Aborting."),
        };

        return Ok((Traversal::Pairs(set), filter));
    }

    if !args.regions.is_empty() {
        let traversal = if args.merge_regions {
            Traversal::Set(region_set()?)
        } else {
            Traversal::Regions(regions)
        };

        filter.regions = bed;
        return Ok((traversal, filter));
    }

    let Some(bed) = bed else {
        return Ok((Traversal::All, filter));
    };

    let indexed = path.is_some()
        && matches!(reader.format(), Format::Bam | Format::Cram)
        && (args.index.is_some() || path.is_some_and(has_index));

    if indexed {
        Ok((Traversal::Set(bed), filter))
    } else {
        filter.regions = Some(bed);
        Ok((Traversal::All, filter))
    }
}

/// Returns the regions that name a reference sequence of the header (or are `*`). Like samtools,
/// warns about the others and skips them, so that a list of only unknown regions selects no
/// records.
fn known_regions(header: &sam::Header, regions: &[String]) -> Vec<String> {
    let getid = |name: &[u8]| crate::regions::reference_sequence_index(header, name);

    regions
        .iter()
        .filter(|region| {
            let known = region.as_str() == UNMAPPED_REGION
                || crate::hts_region::parse_region(region, getid).is_some();

            if !known {
                eprintln!(
                    "[main_samview] region \"{region}\" specifies an invalid region or unknown reference. Continue anyway."
                );
            }

            known
        })
        .cloned()
        .collect()
}

fn has_index(path: &Path) -> bool {
    ["bai", "csi", "crai"].iter().any(|ext| {
        let mut s = path.as_os_str().to_owned();
        s.push(".");
        s.push(ext);
        Path::new(&s).exists()
    }) || path.with_extension("bai").exists()
}

/// Whether CRAM records are decoded on worker threads (whole-file reads and indexed queries
/// with `-@`).
#[derive(Clone, Copy)]
enum Decoding<'s, 'q> {
    Serial,
    Parallel(Workers<'s, 'q>),
}

/// The records worker threads decode, and how many workers decode them.
#[derive(Clone, Copy)]
struct Workers<'s, 'q> {
    source: Source<'s, 'q>,
    threads: usize,
}

/// Record edits samtools makes around filtering (`-B` after the expression, `-p` for rejected
/// records), which work on an encoded copy of each record.
struct Edits {
    remove_backward_ops: bool,
}

impl Edits {
    /// Returns `None` when no record is edited, so records are filtered as they are read.
    fn new(args: &Args) -> Option<Self> {
        // `-p` is ignored when counting.
        let unmaps = args.unmap && !args.count;

        (args.remove_backward_ops || unmaps).then_some(Self {
            remove_backward_ops: args.remove_backward_ops,
        })
    }
}

/// Applies the filter, and any edits, to records.
struct Selector<'f> {
    filter: &'f Filter,
    edits: Option<Edits>,
    encoder: RawEncoder,
    /// The edited copy of the current record.
    record: bam::Record,
}

impl<'f> Selector<'f> {
    fn new(filter: &'f Filter, edits: Option<Edits>) -> Self {
        Self {
            filter,
            edits,
            encoder: RawEncoder::new(),
            record: bam::Record::default(),
        }
    }

    /// Returns whether a record is selected, and the record as edited (`process_aln`): tested
    /// against the expression, then stripped of `B` operations, then tested against the other
    /// filters.
    fn select<'a>(&'a mut self, rec: Rec<'a>, header: &sam::Header) -> io::Result<(bool, Rec<'a>)> {
        let Some(edits) = &self.edits else {
            let selected = self.filter.matches(&rec, header)?;
            return Ok((selected, rec));
        };

        let mut buf = Vec::from(mem::take(&mut self.record));
        buf.clear();
        buf.extend_from_slice(rec.encode_raw(header, &mut self.encoder)?);
        self.record = bam::Record::try_from(buf)?;

        if !self
            .filter
            .matches_expression(&Rec::Bam(&self.record), header)?
        {
            return Ok((false, Rec::Bam(&self.record)));
        }

        if edits.remove_backward_ops {
            let mut buf = Vec::from(mem::take(&mut self.record));
            bam_raw::remove_backward_ops(&mut buf)?;
            self.record = bam::Record::try_from(buf)?;
        }

        let selected = self
            .filter
            .matches_fields(&Rec::Bam(&self.record), header)?;
        Ok((selected, Rec::Bam(&self.record)))
    }

    /// Unmaps the record last selected (`-p`); only valid with edits.
    fn unmap(&mut self) -> io::Result<Rec<'_>> {
        let mut buf = Vec::from(mem::take(&mut self.record));
        bam_raw::unmap(&mut buf)?;
        self.record = bam::Record::try_from(buf)?;
        Ok(Rec::Bam(&self.record))
    }
}

/// Where records that fail the filters go.
enum Rejected {
    Dropped,
    /// Unmapped and written to the output (`-p`).
    Unmapped,
    /// Written as they are to another file (`-U`).
    Written(Box<AlignmentWriter>),
}

/// Edits samtools makes to the records it writes (`adjust_tags`, `change_flag`).
struct OutputEdits {
    tags: TagEdit,
    add_flags: u16,
    remove_flags: u16,
}

/// Which tags the records written keep.
enum TagEdit {
    All,
    /// Drop these (`-x`).
    Remove(Vec<[u8; 2]>),
    /// Keep only these (`--keep-tag`, `-x ^TAGS`).
    Keep(Vec<[u8; 2]>),
}

impl OutputEdits {
    /// Returns `None` when records are written as they are.
    fn new(args: &Args) -> anyhow::Result<Option<Self>> {
        let mut remove = Vec::new();
        // A keep list, once given, applies even if empty (keeping no tags), as in samtools.
        let mut keep: Option<Vec<[u8; 2]>> = None;

        for raw in &args.remove_tags {
            match raw.strip_prefix('^') {
                Some(tags) => keep.get_or_insert_default().extend(parse_tags(tags)?),
                None => remove.extend(parse_tags(raw)?),
            }
        }

        for raw in &args.keep_tags {
            keep.get_or_insert_default().extend(parse_tags(raw)?);
        }

        let tags = if let Some(keep) = keep {
            TagEdit::Keep(keep)
        } else if !remove.is_empty() {
            TagEdit::Remove(remove)
        } else {
            TagEdit::All
        };

        let bits = |values: &[Flags]| flags::union(values).map_or(0, |f| f.bits());
        let add_flags = bits(&args.add_flags);
        let remove_flags = bits(&args.remove_flags);

        Ok(
            (!matches!(tags, TagEdit::All) || add_flags != 0 || remove_flags != 0).then_some(
                Self {
                    tags,
                    add_flags,
                    remove_flags,
                },
            ),
        )
    }

    /// Writes the edited copy of the raw record `src` to `dst`.
    fn apply(&self, src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
        match &self.tags {
            TagEdit::All => {
                dst.clear();
                dst.extend_from_slice(src);
            }
            TagEdit::Remove(tags) => {
                let keep_cg = tags.contains(CG) && has_long_cigar(src)?;
                bam_raw::strip_fields(src, dst, |tag| {
                    tags.contains(&tag) && !(keep_cg && tag == *CG)
                })?;
            }
            TagEdit::Keep(tags) => {
                let keep_cg = !tags.contains(CG) && has_long_cigar(src)?;
                bam_raw::strip_fields(src, dst, |tag| {
                    !(tags.contains(&tag) || (keep_cg && tag == *CG))
                })?;
            }
        }

        if self.add_flags != 0 || self.remove_flags != 0 {
            bam_raw::edit_flags(dst, self.add_flags, self.remove_flags)?;
        }

        Ok(())
    }
}

/// The tag holding a long CIGAR.
const CG: &[u8; 2] = b"CG";

/// Whether a raw record's CIGAR is held in its `CG` field (more than 65,535 operations). That field
/// is part of the alignment, not an auxiliary tag to filter: HTSlib moves it into the CIGAR when it
/// reads the record, so samtools' tag filters never see it.
fn has_long_cigar(src: &[u8]) -> io::Result<bool> {
    Ok(bam_raw::Record::new(src)?.long_cigar_field()?.is_some())
}

/// Filters records and writes them.
struct Sink<'a, 'f> {
    selector: Selector<'f>,
    /// The header records are written with (the input header decodes them).
    output_header: &'a sam::Header,
    writer: &'a mut AlignmentWriter,
    rejected: &'a mut Rejected,
    output_edits: Option<&'a OutputEdits>,
    encoder: RawEncoder,
    edited: Vec<u8>,
    counts: Counts,
}

impl Sink<'_, '_> {
    /// Filters and writes a record decoded with `header` (the input header).
    fn emit(&mut self, rec: Rec<'_>, header: &sam::Header) -> anyhow::Result<()> {
        let (selected, rec) = self.selector.select(rec, header)?;
        let output_header = self.output_header;

        self.counts.processed += 1;
        self.counts.accepted += u64::from(selected);

        if selected {
            if let Some(edits) = self.output_edits {
                let raw = rec.encode_raw(header, &mut self.encoder)?;
                edits.apply(raw, &mut self.edited)?;
                self.writer.write_raw(output_header, &self.edited)?;
            } else {
                self.writer.write(output_header, &rec)?;
            }

            return Ok(());
        }

        match self.rejected {
            Rejected::Dropped => {}
            Rejected::Written(unselected) => unselected.write(output_header, &rec)?,
            Rejected::Unmapped => {
                let rec = self.selector.unmap()?;
                self.writer.write(output_header, &rec)?;
            }
        }

        Ok(())
    }
}

/// The CRAM record fields and tags counting reads: those the filters read, and the alignment
/// interval for region queries. Edits re-encode whole records.
fn count_fields(
    filter: &Filter,
    traversal: &Traversal,
    edits: Option<&Edits>,
) -> (cram::io::reader::Fields, cram::io::reader::Tags) {
    use cram::io::reader::{Fields, Tags};

    if edits.is_some() {
        return (Fields::all(), Tags::All);
    }

    let (mut fields, tags) = filter.required_fields();

    if !matches!(traversal, Traversal::All) {
        fields |= Fields::REFERENCE_SEQUENCE_ID | Fields::ALIGNMENT_START | Fields::CIGAR;
    }

    (fields, tags)
}

fn count(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    args: &Args,
    header: &sam::Header,
    traversal: &Traversal,
    selector: &mut Selector<'_>,
    decoding: Decoding<'_, '_>,
) -> anyhow::Result<Counts> {
    let mut n = 0u64;
    let mut processed = 0u64;

    if let (AlignmentReader::Cram(cram_reader, repository, ..), Decoding::Parallel(workers)) =
        (&mut *reader, decoding)
    {
        let filter = selector.filter;

        parallel_cram::decode(
            cram_reader,
            header,
            repository,
            workers.source,
            workers.threads,
            |count: &mut u64, record| {
                if filter.matches(&Rec::Dyn(record), header)? {
                    *count += 1;
                }
                Ok(())
            },
            |count| {
                n += count;
                Ok(())
            },
        )?;

        return Ok(Counts {
            processed: n,
            accepted: n,
        });
    }

    visit(
        reader,
        path,
        args.index.as_deref(),
        header,
        traversal,
        |rec| {
            processed += 1;
            if selector.select(rec, header)?.0 {
                n += 1;
            }
            Ok(())
        },
    )?;

    Ok(Counts {
        processed,
        accepted: n,
    })
}

/// `-P` (`fetch_pairs_collect_mates`): reads the regions once to collect the names and mate
/// positions of the paired reads that pass the filters and whose mates lie outside the regions,
/// then reads the regions and the mate positions (and the unplaced reads, if a mate is
/// unplaced), writing the records that overlap the regions or have a collected name.
fn fetch_pairs(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    index: Option<&Path>,
    header: &sam::Header,
    regions: &RegionSet,
    sink: &mut Sink<'_, '_>,
) -> anyhow::Result<()> {
    let Some(path) = path else {
        bail!("-P needs an indexed file, not stdin");
    };

    let mut names: HashSet<BString> = HashSet::new();
    let mut mates = regions.clone();
    let mut unplaced_mates = false;

    reader.for_each_in_queries(path, index, header, &Queries::Set(regions), |rec| {
        if !rec.flags()?.is_segmented() {
            return Ok(());
        }

        let mate_id = rec.mate_reference_sequence_id(header)?;
        let mate_start = rec.as_dyn().mate_alignment_start().transpose()?;

        if let (Some(id), Some(start)) = (mate_id, mate_start)
            && mate_in_regions(regions, id, usize::from(start))
        {
            return Ok(());
        }

        let name = record_name(&rec);

        if !sink.selector.select(rec, header)?.0 {
            return Ok(());
        }

        names.insert(name);

        match (mate_id, mate_start) {
            (None, _) => unplaced_mates = true,
            (Some(id), Some(start)) => {
                let start = usize::from(start);
                mates.add_interval(id, start, start);
            }
            (Some(_), None) => {}
        }

        Ok(())
    })?;

    mates.merge();

    let mut emit = |rec: Rec<'_>| -> anyhow::Result<()> {
        let in_regions = record_interval(rec.as_dyn(), header)?
            .is_some_and(|(id, start, end)| regions.first_overlap(id, start, end).is_some());

        if in_regions || names.contains(&record_name(&rec)) {
            sink.emit(rec, header)?;
        }

        Ok(())
    };

    reader.for_each_in_queries(path, index, header, &Queries::Set(&mates), &mut emit)?;

    if unplaced_mates {
        let unmapped = [UNMAPPED_REGION.to_owned()];
        reader.for_each_in_queries(path, index, header, &Queries::Regions(&unmapped), emit)?;
    }

    Ok(())
}

/// Whether a mate at 1-based `position` is in the regions as samtools tests it: as the empty
/// interval at the position (`bed_overlap(mpos, mpos)`), which overlaps an interval only if the
/// interval also covers the base before it.
fn mate_in_regions(regions: &RegionSet, reference_sequence_id: usize, position: usize) -> bool {
    regions
        .first_overlap(reference_sequence_id, position, position)
        .is_some_and(|i| regions.intervals(reference_sequence_id)[i].0 < position)
}

fn record_name(rec: &Rec<'_>) -> BString {
    rec.as_dyn()
        .name()
        .map_or_else(|| BString::from("*"), |name| name.to_owned())
}

/// Parses a comma-separated tag list (`parse_aux_list`).
/// Parses a tag list as samtools' `parse_aux_list` does: two-character tags, each followed by a
/// comma or the end (so a trailing comma is accepted, and an empty list is a list).
fn parse_tags(raw: &str) -> anyhow::Result<Vec<[u8; 2]>> {
    let mut tags = Vec::new();
    let mut rest = raw.as_bytes();

    while let [a, b, tail @ ..] = rest {
        tags.push([*a, *b]);
        rest = tail;

        match rest {
            [b',', tail @ ..] => rest = tail,
            [] => {}
            _ => break,
        }
    }

    if !rest.is_empty() {
        bail!("Error parsing option, auxiliary tags should be exactly two characters long");
    }

    Ok(tags)
}

/// SAM text for the records of one CRAM container that pass the filter, and a buffer for
/// encoding records.
#[derive(Default)]
struct SamChunk {
    text: Vec<u8>,
    record: bam::Record,
}

/// `block_size`-prefixed BAM records for one CRAM container that pass the filter, and a buffer
/// for encoding records.
#[derive(Default)]
struct BamChunk {
    records: Vec<u8>,
    record: bam::Record,
}

/// Decodes CRAM records on worker threads, formatting them there as SAM or BAM.
fn write_cram_in_parallel<R>(
    reader: &mut cram::io::Reader<R>,
    repository: &fasta::Repository,
    header: &sam::Header,
    filter: &Filter,
    workers: Workers<'_, '_>,
    kind: OutputFormat,
    writer: &mut AlignmentWriter,
) -> anyhow::Result<()>
where
    R: io::Read + io::Seek + Send,
{
    // Records are encoded as BAM on the workers, so SAM output uses the raw BAM formatter.
    match kind {
        OutputFormat::Sam | OutputFormat::SamGz => {
            let names = Names::new(header);

            parallel_cram::decode(
                reader,
                header,
                repository,
                workers.source,
                workers.threads,
                |chunk: &mut SamChunk, record| {
                    encode_cram_record(header, record, &mut chunk.record)?;

                    if filter.matches(&Rec::Bam(&chunk.record), header)?
                        && !bam_sam::format_record(&names, chunk.record.as_ref(), &mut chunk.text)
                    {
                        let mut writer = sam::io::Writer::new(&mut chunk.text);
                        writer.write_alignment_record(header, &chunk.record)?;
                    }

                    Ok(())
                },
                |chunk| Ok(writer.write_sam_text(chunk.text)?),
            )
        }
        OutputFormat::Bam => parallel_cram::decode(
            reader,
            header,
            repository,
            workers.source,
            workers.threads,
            |chunk: &mut BamChunk, record| {
                encode_cram_record(header, record, &mut chunk.record)?;

                if filter.matches(&Rec::Bam(&chunk.record), header)? {
                    let raw = chunk.record.as_ref();
                    let block_size = u32::try_from(raw.len())
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                    chunk.records.extend(block_size.to_le_bytes());
                    chunk.records.extend_from_slice(raw);
                }

                Ok(())
            },
            |chunk| Ok(writer.write_bam_records(&chunk.records)?),
        ),
        OutputFormat::Cram => bail!("CRAM-to-CRAM output is written serially"),
    }
}

fn visit<F>(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    index: Option<&Path>,
    header: &sam::Header,
    traversal: &Traversal,
    f: F,
) -> anyhow::Result<()>
where
    F: FnMut(Rec<'_>) -> anyhow::Result<()>,
{
    if let Traversal::All = traversal {
        return reader.for_each(header, f);
    }

    let Some(path) = path else {
        bail!("region queries need an indexed file, not stdin");
    };

    reader.for_each_in_queries(path, index, header, &queries(traversal), f)
}

/// The indexed queries of a traversal that is not of every record.
fn queries(traversal: &Traversal) -> Queries<'_> {
    match traversal {
        Traversal::All => unreachable!("every record is read without an index"),
        Traversal::Regions(regions) => Queries::Regions(regions),
        Traversal::Set(set) | Traversal::Pairs(set) => Queries::Set(set),
    }
}

/// An output's format and whether it has a header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct OutputKind {
    format: OutputFormat,
    /// Always true for BAM and CRAM.
    writes_header: bool,
}

/// Resolves the format and header of an output (`-o` or `-U`) as samtools does for each:
/// `-O`, `-b`, or `-C` set the format of both outputs; otherwise each takes its extension, else
/// BAM with `-u` or `-1`, else SAM. SAM output has a header with `-h` or `-H`, or when the format
/// without `-O` would be BAM or CRAM (so `-O sam -o out.bam` writes one).
fn output_kind(args: &Args, dst: Option<&Path>) -> anyhow::Result<OutputKind> {
    let explicit = [
        args.output_format,
        args.bam.then_some(OutputFormat::Bam),
        args.cram.then_some(OutputFormat::Cram),
    ];

    let mut chosen = explicit.into_iter().flatten();

    // As in samtools, `-u` and `-1` (a compression level) default to BAM.
    let level_default = if args.uncompressed || args.fast {
        OutputFormat::Bam
    } else {
        OutputFormat::Sam
    };

    // The format without `-O`, which also decides whether SAM output has a header.
    let implied = if args.bam {
        OutputFormat::Bam
    } else if args.cram {
        OutputFormat::Cram
    } else {
        dst.and_then(OutputFormat::from_path)
            .unwrap_or(level_default)
    };

    let format = match (chosen.next(), chosen.next()) {
        (Some(a), Some(b)) if a != b => bail!("conflicting output formats: {a:?} and {b:?}"),
        (Some(format), _) => format,
        (None, _) => implied,
    };

    let writes_header =
        !format.is_sam() || !implied.is_sam() || args.with_header || args.header_only;

    Ok(OutputKind {
        format,
        writes_header,
    })
}

#[cfg(test)]
mod tag_list_tests {
    use super::parse_tags;

    #[test]
    fn test_parse_tags() {
        assert_eq!(parse_tags("").unwrap(), Vec::<[u8; 2]>::new());
        assert_eq!(parse_tags("XA").unwrap(), [*b"XA"]);
        assert_eq!(parse_tags("XA,").unwrap(), [*b"XA"]);
        assert_eq!(parse_tags("XA,NM").unwrap(), [*b"XA", *b"NM"]);

        for raw in [",XA", ",", "XAB", "X", "XA,,NM"] {
            assert!(parse_tags(raw).is_err(), "{raw}");
        }
    }
}

#[cfg(test)]
mod output_kind_tests {
    use clap::Parser;

    use super::{Args, OutputFormat, OutputKind, Path, output_kind};

    #[derive(Parser)]
    #[command(disable_help_flag = true)]
    struct Cli {
        #[command(flatten)]
        args: Args,
    }

    /// The kind of output `dst` with these view options.
    fn kind(options: &[&str], dst: &str) -> OutputKind {
        let argv = ["view"].iter().chain(options).chain(&["in.bam"]);
        let cli = Cli::try_parse_from(argv).expect("valid options");
        output_kind(&cli.args, Some(Path::new(dst))).expect("one output format")
    }

    fn sam(writes_header: bool) -> OutputKind {
        OutputKind {
            format: OutputFormat::Sam,
            writes_header,
        }
    }

    fn binary(format: OutputFormat) -> OutputKind {
        OutputKind {
            format,
            writes_header: true,
        }
    }

    #[test]
    fn test_output_kind() {
        // Each output takes its own extension.
        assert_eq!(kind(&[], "out.sam"), sam(false));
        assert_eq!(kind(&[], "out.bam"), binary(OutputFormat::Bam));
        assert_eq!(kind(&[], "out.cram"), binary(OutputFormat::Cram));
        assert_eq!(
            kind(&[], "out.bam##idx##out.csi"),
            binary(OutputFormat::Bam)
        );
        assert_eq!(kind(&["-h"], "out.sam"), sam(true));
        assert_eq!(kind(&["-H"], "out.sam"), sam(true));

        // -b, -C, and -O set every output's format.
        assert_eq!(kind(&["-b"], "out.sam"), binary(OutputFormat::Bam));
        assert_eq!(kind(&["-O", "bam"], "out.sam"), binary(OutputFormat::Bam));
        assert_eq!(kind(&["-C"], "out.bam"), binary(OutputFormat::Cram));

        // SAM output has a header when the format without -O is BAM or CRAM.
        assert_eq!(kind(&["-O", "sam"], "out.bam"), sam(true));
        assert_eq!(kind(&["-O", "sam"], "out.sam"), sam(false));
        assert_eq!(kind(&["-O", "sam", "-u"], "out"), sam(true));

        // -u and -1 choose BAM only for outputs without a known extension.
        assert_eq!(kind(&["-u"], "out"), binary(OutputFormat::Bam));
        assert_eq!(kind(&["-1"], "out.sam"), sam(false));
    }
}
