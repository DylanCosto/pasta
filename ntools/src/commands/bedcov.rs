//! `ntools bedcov`: summed read depth over BED regions (samtools bedcov).

use std::{
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::{Context, bail};
use pasta::sam::alignment::record::Flags;

use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, Verbosity};
use crate::{
    alignment::{AlignmentIndex, AlignmentReader, Queries, ReadOptions},
    depth::{CountOptions, DEFAULT_EXCLUDE_FLAGS, DepthWindow},
    filter::Filter,
    flags, output,
    regions::{self, record_interval},
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// BED file of regions.
    bed: PathBuf,

    /// Indexed BAM or CRAM files; with `-X`, followed by their index files.
    #[arg(required = true)]
    inputs: Vec<String>,

    /// The inputs are followed by their index files, in the same order (`-X in.bed a.bam b.bam
    /// a.bai b.bai`).
    #[arg(short = 'X')]
    customized_index: bool,

    /// Exclude reads with mapping quality below this value.
    #[arg(
        short = 'Q',
        long = "min-MQ",
        visible_alias = "min-mq",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<u8>
    )]
    min_mapping_quality: u8,

    /// Remove these flags from the default exclusion list (UNMAP,SECONDARY,QCFAIL,DUP); may be
    /// repeated.
    #[arg(short = 'g', value_parser = flags::parse)]
    unexclude_flags: Vec<Flags>,

    /// Add these flags to the exclusion list; may be repeated.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_flags: Vec<Flags>,

    /// Do not count deletions and reference skips.
    #[arg(short = 'j')]
    no_deletions: bool,

    /// Also report the number of reference bases with at least this depth.
    #[arg(short = 'd')]
    depth_threshold: Option<u32>,

    /// Leave out reads that start where the previous one did once this many are in the pileup
    /// (HTSlib's pileup depth limit; at least the `-d` threshold).
    #[arg(long = "max-depth", value_name = "INT", allow_negative_numbers = true)]
    max_depth: Option<i64>,

    /// Also report read counts.
    #[arg(short = 'c')]
    read_counts: bool,

    /// Print a header line.
    #[arg(short = 'H')]
    print_header: bool,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional decompression threads.
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

/// Per-input results for one region.
#[derive(Clone, Copy, Debug, Default)]
struct RegionCounts {
    depth_sum: u64,
    deep_bases: u64,
    reads: u64,
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    let mut options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };
    args.input_format_options.configure(&mut options)?;

    let (sources, indexes) = if args.customized_index {
        if !args.inputs.len().is_multiple_of(2) {
            bail!("odd number of filenames detected! Each BAM file should have an index file");
        }

        let (sources, indexes) = args.inputs.split_at(args.inputs.len() / 2);
        (sources, Some(indexes))
    } else {
        (args.inputs.as_slice(), None)
    };

    let mut inputs = Vec::with_capacity(sources.len());

    for (i, src) in sources.iter().enumerate() {
        let (mut reader, path) = AlignmentReader::open(src, &options)?;
        args.input_format_options.apply(&mut reader)?;
        let header = reader.read_header()?;
        let path = path.with_context(|| format!("{src}: bedcov needs indexed files, not stdin"))?;
        let explicit = indexes.map(|indexes| Path::new(&indexes[i]));
        let index = match AlignmentIndex::read(reader.format(), &path, explicit) {
            Ok(index) => index,
            // samtools bedcov exits with status 2 when it cannot load an index.
            Err(e) => {
                eprintln!("ntools: {e:#}");
                return Ok(ExitCode::from(2));
            }
        };
        inputs.push((reader, index, header));
    }

    let exclude = (DEFAULT_EXCLUDE_FLAGS
        - flags::union(&args.unexclude_flags).unwrap_or_else(Flags::empty))
        | flags::union(&args.exclude_flags).unwrap_or_else(Flags::empty);

    let filter = Filter {
        exclude,
        min_mapping_quality: args.min_mapping_quality,
        ..Filter::default()
    };

    // Plain or gzip/BGZF compressed, as samtools reads it with `gzopen`.
    let text = regions::read_text_file(&args.bed)?;

    let mut out = BufWriter::new(output::raw_stdout()?);

    if args.print_header {
        write_header(&mut out, &args, sources)?;
    }

    let scan = Scan {
        filter,
        count_options: CountOptions {
            include_deletions: !args.no_deletions,
            include_skips: !args.no_deletions,
            ..CountOptions::default()
        },
        depth_threshold: args.depth_threshold,
        // samtools raises the limit to the -d threshold; a negative one acts as 0.
        max_count: args.max_depth.map(|max_depth| {
            let threshold = args.depth_threshold.map_or(-1, i64::from);
            usize::try_from(max_depth.max(threshold)).unwrap_or(0)
        }),
    };

    let mut counts = vec![RegionCounts::default(); inputs.len()];

    let mut status = ExitCode::SUCCESS;

    for line in text.split(|&c| c == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);

        if line.is_empty()
            || line.starts_with(b"#")
            || line.starts_with(b"track ")
            || line.starts_with(b"browser ")
        {
            continue;
        }

        let Some((name, start, end)) = parse_bed_line(line, &inputs[0].2) else {
            eprintln!("Errors in BED line '{}'", String::from_utf8_lossy(line));
            status = ExitCode::from(2);
            continue;
        };

        for ((reader, index, header), counts) in inputs.iter_mut().zip(counts.iter_mut()) {
            *counts = count_region(reader, index, header, &scan, (name, start, end))?;
        }

        out.write_all(line)?;

        for c in &counts {
            write!(out, "\t{}", c.depth_sum)?;
        }

        if args.depth_threshold.is_some() {
            for c in &counts {
                write!(out, "\t{}", c.deep_bases)?;
            }
        }

        if args.read_counts {
            for c in &counts {
                write!(out, "\t{}", c.reads)?;
            }
        }

        writeln!(out)?;
    }

    out.flush()?;

    Ok(status)
}

/// A BED line's reference sequence name and 0-based start and end, read as samtools does: the
/// name up to the first whitespace, then two integers. `None` if they are missing, the name is
/// not a reference sequence of `header`, or the end comes before the start.
fn parse_bed_line<'l>(
    line: &'l [u8],
    header: &pasta::sam::Header,
) -> Option<(&'l str, usize, usize)> {
    let name_len = line.iter().position(u8::is_ascii_whitespace)?;
    let (name, rest) = (&line[..name_len], &line[name_len + 1..]);

    header.reference_sequences().get_index_of(name)?;
    let name = str::from_utf8(name).ok()?;

    let (start, rest) = scan_integer(rest)?;
    let (end, _) = scan_integer(rest)?;

    if end < start {
        return None;
    }

    // Like HTSlib's iterators, negative coordinates count from 0.
    let clamp = |n: i64| usize::try_from(n.max(0)).unwrap_or(usize::MAX);
    Some((name, clamp(start), clamp(end)))
}

/// `sscanf("%ld")`: skips leading whitespace and reads an optionally signed integer prefix;
/// returns it and the rest of the string.
fn scan_integer(s: &[u8]) -> Option<(i64, &[u8])> {
    let s = s.trim_ascii_start();
    let sign_len = usize::from(s.starts_with(b"+") || s.starts_with(b"-"));
    let digits = s[sign_len..]
        .iter()
        .take_while(|c| c.is_ascii_digit())
        .count();

    if digits == 0 {
        return None;
    }

    let (number, rest) = s.split_at(sign_len + digits);
    let n = str::from_utf8(number).ok()?.parse().ok()?;
    Some((n, rest))
}

/// How every region is counted.
struct Scan {
    filter: Filter,
    count_options: CountOptions,
    /// `-d`.
    depth_threshold: Option<u32>,
    /// HTSlib's pileup `maxcnt` (`--max-depth`).
    max_count: Option<usize>,
}

/// Sums depth over the 0-based half-open region `[start, end)` of one input.
fn count_region(
    reader: &mut AlignmentReader,
    index: &AlignmentIndex,
    header: &pasta::sam::Header,
    scan: &Scan,
    (name, start, end): (&str, usize, usize),
) -> anyhow::Result<RegionCounts> {
    let Scan {
        filter,
        count_options,
        depth_threshold,
        max_count,
    } = scan;
    let (count_options, depth_threshold) = (*count_options, *depth_threshold);

    let mut counts = RegionCounts::default();

    let Some(id) = header.reference_sequences().get_index_of(name.as_bytes()) else {
        return Ok(counts);
    };

    if end <= start {
        return Ok(counts);
    }

    let (first, last) = (start + 1, end);
    let region = format!("{name}:{first}-{last}");
    let regions = [region];

    let mut window = match *max_count {
        Some(max_count) => DepthWindow::with_max_count(count_options, max_count),
        None => DepthWindow::new(count_options),
    };

    let mut emit = |rid: usize, position: usize, cell: crate::depth::Cell| -> io::Result<()> {
        if rid == id && (first..=last).contains(&position) {
            counts.depth_sum += u64::from(cell.depth);

            // As in samtools, the threshold applies to aligned bases, without deletions and
            // reference skips.
            let aligned = if count_options.include_deletions {
                cell.depth - cell.gaps
            } else {
                cell.depth
            };

            if depth_threshold.is_some_and(|d| aligned >= d) {
                counts.deep_bases += 1;
            }
        }

        Ok(())
    };

    let mut reads = 0;

    reader.for_each_in_indexed_queries(index, header, &Queries::Regions(&regions), |record| {
        if !filter.matches(&record, header)? {
            return Ok(());
        }

        // samtools counts reads as its pileup buffers them.
        if window.add(&record, header, &mut emit)?
            && let Some((rid, s, e)) = record_interval(record.as_dyn(), header)?
            && rid == id
            && s <= last
            && e >= first
        {
            reads += 1;
        }

        Ok(())
    })?;

    window.finish(&mut emit)?;
    counts.reads = reads;

    Ok(counts)
}

fn write_header<W: Write>(out: &mut W, args: &Args, sources: &[String]) -> io::Result<()> {
    write!(out, "#chrom\tchromStart\tchromEnd")?;

    for src in sources {
        write!(out, "\t{src}_cov")?;
    }

    if args.depth_threshold.is_some() {
        for src in sources {
            write!(out, "\t{src}_depth")?;
        }
    }

    if args.read_counts {
        for src in sources {
            write!(out, "\t{src}_count")?;
        }
    }

    writeln!(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed_line() -> anyhow::Result<()> {
        let header: pasta::sam::Header = "@SQ\tSN:c1\tLN:100\n".parse()?;

        // Extra columns, which bedcov prints back, may hold any bytes.
        assert_eq!(
            parse_bed_line(b"c1\t10\t20\tname\t\xff", &header),
            Some(("c1", 10, 20))
        );
        assert_eq!(parse_bed_line(b"c1 -5   +20", &header), Some(("c1", 0, 20)));
        assert_eq!(parse_bed_line(b"c1\t20\t10", &header), None);
        assert_eq!(parse_bed_line(b"c1\t20", &header), None);
        assert_eq!(parse_bed_line(b"c2\t10\t20", &header), None);
        assert_eq!(parse_bed_line(b"c1", &header), None);

        Ok(())
    }
}
