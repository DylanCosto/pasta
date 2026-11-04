//! `ntools collate`: group records by read name (samtools collate).
//!
//! Records are grouped by sorting on a hash of the name, so the records of a template are
//! adjacent: READ1 before READ2, as in samtools, and otherwise in input order. The order of
//! templates differs from samtools', which depends on its internal hashing; tools that consume
//! collated input (fastq, fixmate) only need the grouping.

use std::path::{Path, PathBuf};

use anyhow::bail;
use pasta::sam::{
    self as sam,
    alignment::record::Flags,
    header::record::value::{
        Map,
        map::header::tag::{GROUP_ORDER, SORT_ORDER},
    },
};

use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, OutputFormatOptions, Verbosity};
use crate::{
    alignment::{AlignmentReader, ReadOptions},
    commands::sort::{SortedRecords, create_temp_dir, parse_memory},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    sort_order::Order,
};

/// Memory for buffering records before spilling to temporary files.
const DEFAULT_MEMORY: &str = "768M";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Write `PREFIX.bam` (when neither -o nor -O is given).
    prefix: Option<PathBuf>,

    /// Write output to stdout.
    #[arg(short = 'O', conflicts_with = "output")]
    stdout: bool,

    /// Write output to this file.
    #[arg(short = 'o')]
    output: Option<PathBuf>,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Only keep primary alignments (samtools' fast mode). Unlike samtools, records without
    /// exactly one of READ1 and READ2 (single-end reads) are kept.
    #[arg(short = 'f')]
    fast: bool,

    /// BGZF compression level.
    #[arg(short = 'l', default_value_t = 1)]
    level: u8,

    /// Memory for buffering records before spilling to temporary files (K/M/G suffixes; at
    /// least 1M).
    #[arg(short = 'm', default_value = DEFAULT_MEMORY, value_parser = parse_memory)]
    memory: usize,

    /// Prefix for temporary files (ending in `/`, the directory to hold them).
    #[arg(short = 'T')]
    temp_prefix: Option<PathBuf>,

    /// Number of temporary files. Accepted for compatibility and ignored: it only tunes
    /// samtools' temporary files (ntools buffers up to -m).
    #[arg(short = 'n', value_name = "INT", default_value_t = 64)]
    temp_files: i32,

    /// Working reads stored (with -f). Accepted for compatibility and ignored: it only tunes
    /// samtools' temporary files (ntools buffers up to -m).
    #[arg(short = 'r', value_name = "INT", default_value_t = 10000)]
    reads_store: i32,

    /// Output format [default: with -o, from its extension, else SAM; otherwise BAM].
    #[arg(long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for sorting, compression, and decompression.
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

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    output_format_options: OutputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    // As in samtools, `-o FILE` takes the format of its extension, and SAM for any other name;
    // the prefix and `-O` write BAM.
    let output_format = args.output_format.unwrap_or(match &args.output {
        Some(path) if !args.stdout => OutputFormat::from_path(path).unwrap_or(OutputFormat::Sam),
        _ => OutputFormat::Bam,
    });

    let output: Option<PathBuf> = match (args.stdout, &args.output, &args.prefix) {
        (true, _, _) => None,
        (false, Some(path), _) => Some(path.clone()),
        (false, None, Some(prefix)) => {
            let mut path = prefix.as_os_str().to_owned();
            path.push(".");
            path.push(output_format.extension());
            Some(PathBuf::from(path))
        }
        (false, None, None) => bail!("give an output prefix, -o FILE, or -O for stdout"),
    };

    let mut options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };
    args.input_format_options.configure(&mut options)?;

    let (mut reader, _) = AlignmentReader::open(&args.input, &options)?;
    args.input_format_options.apply(&mut reader)?;
    let mut header = reader.read_header()?;

    let temp_dir = create_temp_dir(
        args.temp_prefix.as_deref(),
        output.as_deref().filter(|p| *p != Path::new("-")),
        false,
    )?;

    let not_primary = Flags::SECONDARY | Flags::SUPPLEMENTARY;
    let fast = args.fast;

    // samtools' fast mode also drops primary records without exactly one of READ1 and READ2,
    // which its usage ("only primary alignments") does not say: single-end input comes out
    // empty. ntools keeps them (README).
    let keep = |flags: Flags| !fast || !flags.intersects(not_primary);

    let sorted = SortedRecords::collect(
        &mut reader,
        &header,
        Order::Collate,
        args.memory.saturating_mul(args.threads + 1),
        args.threads,
        temp_dir.path(),
        |record| Ok(keep(record.flags()?)),
    )?;

    set_group_order(&mut header);

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format: output_format,
        level: output::parse_level(if args.uncompressed {
            0
        } else {
            args.output_format_options.level(args.level)?
        })?,
        threads: args.threads,
        reference: args.reference,
        cram: args.output_format_options.cram()?,
    };

    let mut writer = AlignmentWriter::create(output.as_deref(), &write_options, &header)?;
    writer.write_header(&header)?;
    sorted.write_to(&mut writer, &header)?;
    writer.finish(&header)?;

    Ok(())
}

/// Sets `@HD SO:unsorted GO:query`, as samtools collate does (keeping any `SS`).
fn set_group_order(header: &mut sam::Header) {
    let hd = header.header_mut().get_or_insert_with(Map::default);
    let fields = hd.other_fields_mut();
    fields.insert(SORT_ORDER, "unsorted".into());
    fields.insert(GROUP_ORDER, "query".into());
}
