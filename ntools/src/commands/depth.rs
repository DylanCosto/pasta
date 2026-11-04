//! `ntools depth`: per-position read depth (samtools depth).

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    thread,
};

use anyhow::{Context, bail};
use crossbeam_channel::{Receiver, Sender};
use pasta::cram::io::reader::{Fields, Tags};
use pasta::sam::{self as sam, alignment::record::Flags};

use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, Verbosity};
use crate::{
    alignment::{AlignmentReader, Queries, RawEncoder, ReadOptions, Rec},
    bam_raw::Record,
    depth::{CountOptions, DEFAULT_EXCLUDE_FLAGS, DepthWindow, used_query_length},
    file_list,
    filter::Filter,
    flags, output,
    regions::{self, RegionSet},
};

/// `(reference sequence ID, 1-based position, depth)`.
type Point = (usize, usize, u32);

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Coordinate-sorted SAM, BAM, or CRAM files; with `-X`, followed by their index files.
    #[arg(required_unless_present = "file_list")]
    inputs: Vec<String>,

    /// Read the input file names from this file, one per line, instead of the arguments.
    #[arg(short = 'f')]
    file_list: Option<PathBuf>,

    /// The inputs are followed by their index files, in the same order (`-X a.bam b.bam a.bai
    /// b.bai`), for `-r`.
    #[arg(short = 'X')]
    customized_index: bool,

    /// Output all positions (-a), including unused reference sequences (-aa).
    #[arg(short = 'a', action = clap::ArgAction::Count)]
    all_positions: u8,

    /// Only report this region.
    #[arg(short = 'r')]
    region: Option<String>,

    /// Only report positions in this BED file.
    #[arg(short = 'b')]
    bed: Option<PathBuf>,

    /// Print a header line.
    #[arg(short = 'H')]
    print_header: bool,

    /// Exclude bases with base quality below this value.
    #[arg(
        short = 'q',
        long = "min-BQ",
        alias = "min-bq",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<u8>
    )]
    min_base_quality: u8,

    /// Exclude reads with mapping quality below this value.
    #[arg(
        short = 'Q',
        long = "min-MQ",
        alias = "min-mq",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<u8>
    )]
    min_mapping_quality: u8,

    /// Do not count bases of a read that overlap its mate (the read seen first keeps them).
    #[arg(short = 's')]
    remove_overlaps: bool,

    /// Count deletions.
    #[arg(short = 'J')]
    include_deletions: bool,

    /// Remove these flags from the default exclusion list (UNMAP,SECONDARY,QCFAIL,DUP); may be
    /// repeated.
    #[arg(short = 'g', value_parser = flags::parse)]
    unexclude_flags: Vec<Flags>,

    /// Add these flags to the exclusion list; may be repeated.
    #[arg(short = 'G', long = "excl-flags", value_parser = flags::parse)]
    exclude_flags: Vec<Flags>,

    /// Only include reads with all of these flags set; may be repeated.
    #[arg(long = "require-flags", value_parser = flags::parse)]
    require_flags: Vec<Flags>,

    /// Only include reads with at least one of these flags set; may be repeated.
    #[arg(long = "incl-flags", value_parser = flags::parse)]
    include_flags: Vec<Flags>,

    /// Ignored (samtools' former depth limit).
    #[arg(
        short = 'd',
        short_alias = 'm',
        hide = true,
        action = clap::ArgAction::Append,
        allow_hyphen_values = true
    )]
    max_depth: Vec<String>,

    /// Exclude reads shorter than this.
    #[arg(
        short = 'l',
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    min_read_length: usize,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o')]
    output: Option<PathBuf>,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional decompression threads per input.
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };
    args.input_format_options.configure(&mut options)?;

    let (sources, indexes) = input_files(&args)?;
    let mut inputs = Vec::with_capacity(sources.len());

    for (i, src) in sources.iter().enumerate() {
        let (mut reader, path) = AlignmentReader::open(src, &options)?;
        args.input_format_options.apply(&mut reader)?;
        // As in samtools (CRAM_OPT_DECODE_MD set to 0); depth never reads MD or NM.
        reader.keep_stored_md_nm();
        let (fields, tags) = required_fields(&args);
        reader.set_required_fields(fields, tags);
        let header = reader.read_header()?;
        let index = indexes.as_ref().map(|indexes| PathBuf::from(&indexes[i]));
        inputs.push((reader, path, header, index));
    }

    let header = inputs[0].2.clone();

    for ((_, _, input_header, _), src) in inputs.iter().zip(&sources).skip(1) {
        super::check_reference_sequences("depth", &sources[0], &header, src, input_header)?;
    }

    let exclude = (DEFAULT_EXCLUDE_FLAGS
        - flags::union(&args.unexclude_flags).unwrap_or_else(Flags::empty))
        | flags::union(&args.exclude_flags).unwrap_or_else(Flags::empty);

    let selection = Selection {
        filter: Filter {
            require: flags::union(&args.require_flags).unwrap_or_else(Flags::empty),
            exclude,
            min_mapping_quality: args.min_mapping_quality,
            ..Filter::default()
        },
        include: flags::union(&args.include_flags).unwrap_or_else(Flags::empty),
        min_read_length: args.min_read_length,
    };

    let count_options = CountOptions {
        min_base_quality: args.min_base_quality,
        include_deletions: args.include_deletions,
        track_quality: false,
        include_skips: false,
    };

    let domain = Domain::new(&args, &header)?;

    let sink: output::Sink = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => output::raw_stdout()?,
    };

    let mut out = Output {
        writer: BufWriter::with_capacity(1 << 17, sink),
        names: header
            .reference_sequences()
            .keys()
            .map(|name| name.to_vec())
            .collect(),
        buf: Vec::with_capacity(OUTPUT_BUFFER_LEN + 4096),
    };

    if args.print_header {
        write!(out.writer, "#CHROM\tPOS")?;
        for src in &sources {
            write!(out.writer, "\t{src}")?;
        }
        writeln!(out.writer)?;
    }

    let region_queries = domain.queries();
    let overlaps = if args.remove_overlaps {
        Overlaps::Remove
    } else {
        Overlaps::Count
    };

    thread::scope(|scope| -> anyhow::Result<()> {
        let mut streams = Vec::with_capacity(inputs.len());
        let mut handles = Vec::with_capacity(inputs.len());
        let batch_len = crate::depth::batch_len(inputs.len());

        for (reader, path, input_header, index) in &mut inputs {
            let (tx, rx) = crossbeam_channel::bounded::<Vec<Point>>(4);
            let selection = &selection;
            let region_queries = &region_queries;

            handles.push(scope.spawn(move || {
                read_depths(
                    reader,
                    path.as_deref(),
                    index.as_deref(),
                    input_header,
                    region_queries.as_ref(),
                    selection,
                    count_options,
                    overlaps,
                    batch_len,
                    &tx,
                )
            }));

            streams.push(Stream::new(rx));
        }

        let written = domain.write(&header, &mut streams, &mut out);

        // Dropping the streams lets reader threads stop early if writing failed.
        drop(streams);

        for handle in handles {
            handle
                .join()
                .map_err(|_| anyhow::anyhow!("depth thread panicked"))??;
        }

        written
    })?;

    out.finish()?;

    Ok(())
}

/// Whether bases of a read that overlap its mate are counted (`-s` removes them).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Overlaps {
    Count,
    Remove,
}

/// Computes nonzero depths of one input and sends them in order.
#[expect(
    clippy::too_many_arguments,
    reason = "one call site, all distinct inputs"
)]
fn read_depths(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    index: Option<&Path>,
    header: &sam::Header,
    queries: Option<&RegionSet>,
    selection: &Selection,
    count_options: CountOptions,
    overlaps: Overlaps,
    batch_len: usize,
    tx: &Sender<Vec<Point>>,
) -> anyhow::Result<()> {
    let mut window = DepthWindow::new(count_options);
    let mut batch = Vec::with_capacity(batch_len);

    // The writer stops reading once it has every position it reports (say, the end of a `-r`
    // region), or when it fails, which it reports itself. Either way, reading stops quietly.
    let writer_stopped = std::cell::Cell::new(false);

    // Small enough to inline into the depth window's loops; sending is out of line.
    let mut emit = |id: usize, position: usize, cell: crate::depth::Cell| -> io::Result<()> {
        batch.push((id, position, cell.depth));

        if batch.len() >= batch_len {
            send_batch(&mut batch, batch_len, tx, &writer_stopped)?;
        }

        Ok(())
    };

    let mut mates = (overlaps == Overlaps::Remove).then(MateEnds::default);
    let mut encoder = RawEncoder::new();
    let mut cigar = Vec::new();

    let mut visit = |record: Rec<'_>| -> anyhow::Result<()> {
        if !selection.matches(&record, header)? {
            return Ok(());
        }

        let clip = match &mut mates {
            Some(mates) => {
                let raw = Record::new(record.encode_raw(header, &mut encoder)?)?;
                mates.clip(&raw, &mut cigar)?
            }
            None => 0,
        };

        window.add_clipped(&record, header, clip, &mut emit)?;
        Ok(())
    };

    let read = match (queries, path) {
        (Some(set), Some(path)) => {
            reader.for_each_in_queries(path, index, header, &Queries::Set(set), &mut visit)
        }
        _ => reader.for_each(header, &mut visit),
    }
    .and_then(|()| Ok(window.finish(&mut emit)?));

    if writer_stopped.get() {
        return Ok(());
    }

    read?;

    if !batch.is_empty() {
        // The writer may have stopped; it reports why.
        let _ = tx.send(batch);
    }

    Ok(())
}

/// Sends a full batch of depths, or records that the writer stopped.
#[cold]
#[inline(never)]
fn send_batch(
    batch: &mut Vec<Point>,
    batch_len: usize,
    tx: &Sender<Vec<Point>>,
    writer_stopped: &std::cell::Cell<bool>,
) -> io::Result<()> {
    let full = std::mem::replace(batch, Vec::with_capacity(batch_len));

    if tx.send(full).is_err() {
        writer_stopped.set(true);
        return Err(io::Error::new(
            io::ErrorKind::BrokenPipe,
            "depth writer stopped",
        ));
    }

    Ok(())
}

/// The input files and, with `-X`, their index files.
fn input_files(args: &Args) -> anyhow::Result<(Vec<String>, Option<Vec<String>>)> {
    if let Some(list) = &args.file_list {
        if args.customized_index {
            bail!("The -f option cannot be combined with -X");
        }

        return Ok((file_list::read(list)?, None));
    }

    if args.customized_index {
        // samtools does not catch an odd count (`nfiles%1`) and pairs the files up wrongly.
        if !args.inputs.len().is_multiple_of(2) {
            bail!("-X needs one index specified per bam file");
        }

        let n = args.inputs.len() / 2;
        return Ok((args.inputs[..n].to_vec(), Some(args.inputs[n..].to_vec())));
    }

    Ok((args.inputs.clone(), None))
}

/// The CRAM fields depth reads, as samtools selects them: the rest are not decoded.
fn required_fields(args: &Args) -> (Fields, Tags) {
    // Removing overlaps compares mates by name and position, on a BAM encoding of each record.
    if args.remove_overlaps {
        return (Fields::all(), Tags::All);
    }

    let mut fields = Fields::FLAGS
        | Fields::REFERENCE_SEQUENCE_ID
        | Fields::ALIGNMENT_START
        | Fields::MAPPING_QUALITY
        | Fields::CIGAR;

    if args.min_base_quality > 0 {
        fields |= Fields::QUALITY_SCORES;
    }

    // The read length comes from the sequence.
    if args.min_read_length > 0 {
        fields |= Fields::SEQUENCE;
    }

    (fields, Tags::All)
}

/// Which reads are counted.
struct Selection {
    filter: Filter,
    /// Reads need at least one of these flags (`--incl-flags`; none when empty).
    include: Flags,
    min_read_length: usize,
}

impl Selection {
    fn matches(&self, record: &Rec<'_>, header: &sam::Header) -> io::Result<bool> {
        if !self.filter.matches(record, header)? {
            return Ok(false);
        }

        if !self.include.is_empty() && !record.flags()?.intersects(self.include) {
            return Ok(false);
        }

        Ok(self.min_read_length == 0 || used_query_length(record)? >= self.min_read_length)
    }
}

/// The 0-based end of each paired read whose mate may still overlap it, by name
/// (samtools depth `-s`).
#[derive(Default)]
struct MateEnds(HashMap<Vec<u8>, i64>);

impl MateEnds {
    /// Returns the 1-based position before which the record's bases overlap its mate seen
    /// earlier (0 for none), and records its end for its mate otherwise.
    fn clip(&mut self, record: &Record<'_>, cigar: &mut Vec<u32>) -> io::Result<usize> {
        const PAIRED: u16 = 0x1;
        const MATE_UNMAPPED: u16 = 0x8;

        let flags = record.flags();

        if flags & PAIRED == 0 || flags & MATE_UNMAPPED != 0 {
            return Ok(0);
        }

        if let Some(end) = self.0.remove(record.name()) {
            // The mate's 0-based exclusive end is the 1-based position of its last base + 1.
            return Ok(usize::try_from(end).unwrap_or(0) + 1);
        }

        record.read_cigar(cigar)?;
        let span: i64 = cigar
            .iter()
            .filter(|&&op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8))
            .map(|&op| i64::from(op >> 4))
            .sum();
        let end = i64::from(record.position()) + span.max(1);

        // Mates known not to overlap are not recorded.
        let mate_position = i64::from(record.mate_position());

        if mate_position == -1
            || (record.reference_sequence_id() == record.mate_reference_sequence_id()
                && mate_position <= end)
        {
            self.0.insert(record.name().to_vec(), end);
        }

        Ok(0)
    }
}

/// One input's ordered depths.
struct Stream {
    rx: Receiver<Vec<Point>>,
    batch: Vec<Point>,
    i: usize,
}

impl Stream {
    fn new(rx: Receiver<Vec<Point>>) -> Self {
        Self {
            rx,
            batch: Vec::new(),
            i: 0,
        }
    }

    fn peek(&mut self) -> Option<Point> {
        while self.i >= self.batch.len() {
            self.batch = self.rx.recv().ok()?;
            self.i = 0;
        }

        Some(self.batch[self.i])
    }

    fn advance(&mut self) {
        self.i += 1;
    }
}

/// Rows are formatted into a buffer of about this many bytes before each write.
const OUTPUT_BUFFER_LEN: usize = 1 << 17;

struct Output {
    writer: BufWriter<output::Sink>,
    names: Vec<Vec<u8>>,
    buf: Vec<u8>,
}

impl Output {
    fn write(&mut self, id: usize, position: usize, depths: &[u32]) -> io::Result<()> {
        self.buf.extend_from_slice(&self.names[id]);
        self.buf.push(b'\t');
        push_number(&mut self.buf, position);

        for &depth in depths {
            self.buf.push(b'\t');
            push_number(&mut self.buf, depth);
        }

        self.buf.push(b'\n');

        if self.buf.len() >= OUTPUT_BUFFER_LEN {
            self.writer.write_all(&self.buf)?;
            self.buf.clear();
        }

        Ok(())
    }

    /// Writes the buffered rows and flushes the writer.
    fn finish(&mut self) -> io::Result<()> {
        self.writer.write_all(&self.buf)?;
        self.buf.clear();
        self.writer.flush()
    }
}

fn push_number<N: lexical_core::ToLexical>(dst: &mut Vec<u8>, n: N) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

/// Which positions are reported.
enum Domain {
    /// Every position with depth; with `-a`, all positions of used reference sequences, and with
    /// `-aa`, of all reference sequences.
    All { all_positions: u8 },
    /// Positions inside these intervals; with `-a`, including zero depths on the region's
    /// reference sequence (`-r`) or on reference sequences with depth, and with `-aa`, on all.
    Intervals {
        set: RegionSet,
        all_positions: u8,
        /// Whether the intervals came from a region (`-r`), so reads can be queried by index.
        indexed: bool,
    },
}

impl Domain {
    fn new(args: &Args, header: &sam::Header) -> anyhow::Result<Self> {
        let region = args
            .region
            .as_ref()
            .map(|region| -> anyhow::Result<RegionSet> {
                // A region without an end runs past `@SQ LN`, as HTSlib's does, to the reads
                // that overhang the end of the reference sequence.
                let Some((id, start, end)) = regions::parse_region(header, region) else {
                    bail!("invalid region or unknown reference sequence: {region}");
                };

                let mut set = RegionSet::new(header);
                set.add_interval(id, start, end.unwrap_or(OPEN_END));
                set.merge();
                Ok(set)
            })
            .transpose()?;

        let bed = args
            .bed
            .as_ref()
            .map(|bed| -> anyhow::Result<RegionSet> {
                let mut set = RegionSet::new(header);
                set.add_bed(header, bed)?;
                set.merge();
                Ok(set)
            })
            .transpose()?;

        // As in samtools, a region and a BED file report the positions in both.
        let set = match (region, bed) {
            (Some(region), Some(bed)) => region.intersection(&bed),
            (Some(set), None) | (None, Some(set)) => set,
            (None, None) => {
                return Ok(Self::All {
                    all_positions: args.all_positions,
                });
            }
        };

        Ok(Self::Intervals {
            set,
            all_positions: args.all_positions,
            indexed: args.region.is_some(),
        })
    }

    /// Intervals to read through the index, if any. A region (`-r`) is queried; a BED file is
    /// filtered while streaming, as samtools does without an index.
    fn queries(&self) -> Option<RegionSet> {
        match self {
            Self::Intervals {
                set, indexed: true, ..
            } => Some(set.clone()),
            _ => None,
        }
    }

    fn write(
        &self,
        header: &sam::Header,
        streams: &mut [Stream],
        out: &mut Output,
    ) -> anyhow::Result<()> {
        let mut depths = vec![0; streams.len()];
        let lengths: Vec<usize> = header
            .reference_sequences()
            .values()
            .map(|rs| rs.length().get())
            .collect();

        match self {
            Self::All { all_positions } => {
                let mut current: Option<(usize, usize)> = None;

                while let Some((id, position)) = next_position(streams, &mut depths) {
                    if *all_positions > 0 {
                        let from = match current {
                            Some((cid, cpos)) if cid == id => cpos + 1,
                            Some((cid, cpos)) => {
                                write_zeros(out, cid, cpos + 1, lengths[cid], streams.len())?;
                                if *all_positions > 1 {
                                    write_unused(out, &lengths, cid + 1..id, streams.len())?;
                                }
                                1
                            }
                            None => {
                                if *all_positions > 1 {
                                    write_unused(out, &lengths, 0..id, streams.len())?;
                                }
                                1
                            }
                        };

                        write_zeros(out, id, from, position - 1, streams.len())?;
                    }

                    out.write(id, position, &depths)?;
                    current = Some((id, position));
                }

                if *all_positions > 0 {
                    let next_unused = match current {
                        Some((cid, cpos)) => {
                            write_zeros(out, cid, cpos + 1, lengths[cid], streams.len())?;
                            cid + 1
                        }
                        None => 0,
                    };

                    if *all_positions > 1 {
                        write_unused(out, &lengths, next_unused..lengths.len(), streams.len())?;
                    }
                }
            }
            Self::Intervals {
                set,
                all_positions,
                indexed,
            } => {
                let mut pending = next_position(streams, &mut depths);

                for (id, &length) in lengths.iter().enumerate() {
                    if set.intervals(id).is_empty() {
                        continue;
                    }

                    // Skip depths of earlier reference sequences, to see whether this one has
                    // any: with a BED file alone, `-a` reports zeros only on those that do.
                    while let Some((pid, _)) = pending
                        && pid < id
                    {
                        pending = next_position(streams, &mut depths);
                    }

                    let zeros = match *all_positions {
                        0 => false,
                        1 => *indexed || pending.is_some_and(|(pid, _)| pid == id),
                        _ => true,
                    };

                    // As in samtools, the depths of reads past `@SQ LN` are reported, and so,
                    // with `-a`, are the zeros between them, but not zeros after the last one.
                    for &(start, end) in set.intervals(id) {
                        let mut position = start;

                        while position <= end {
                            // Skip depths before this position.
                            while let Some((pid, ppos)) = pending
                                && (pid, ppos) < (id, position)
                            {
                                pending = next_position(streams, &mut depths);
                            }

                            match pending {
                                Some((pid, ppos)) if pid == id && ppos <= end => {
                                    if zeros {
                                        write_zeros(out, id, position, ppos - 1, streams.len())?;
                                    }
                                    out.write(id, ppos, &depths)?;
                                    position = ppos + 1;
                                    pending = next_position(streams, &mut depths);
                                }
                                Some((pid, _)) if pid == id => {
                                    if zeros {
                                        write_zeros(out, id, position, end, streams.len())?;
                                    }
                                    position = end + 1;
                                }
                                _ => {
                                    if zeros {
                                        let end = end.min(length);
                                        write_zeros(out, id, position, end, streams.len())?;
                                    }
                                    position = end + 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

/// The end of a region given without one (`HTS_POS_MAX`).
const OPEN_END: usize = crate::hts_region::POS_MAX as usize;

/// Pops the smallest position across streams, filling `depths` with each input's depth there.
fn next_position(streams: &mut [Stream], depths: &mut [u32]) -> Option<(usize, usize)> {
    if let ([stream], [depth]) = (&mut *streams, &mut *depths) {
        let (id, position, d) = stream.peek()?;
        stream.advance();
        *depth = d;
        return Some((id, position));
    }

    let next = streams
        .iter_mut()
        .filter_map(|s| s.peek().map(|(id, pos, _)| (id, pos)))
        .min()?;

    for (stream, depth) in streams.iter_mut().zip(depths.iter_mut()) {
        *depth = match stream.peek() {
            Some((id, pos, d)) if (id, pos) == next => {
                stream.advance();
                d
            }
            _ => 0,
        };
    }

    Some(next)
}

/// Writes zero depths for every position of the reference sequences in `ids`.
fn write_unused(
    out: &mut Output,
    lengths: &[usize],
    ids: std::ops::Range<usize>,
    input_count: usize,
) -> io::Result<()> {
    let (start, end) = (ids.start, ids.end);

    for (id, &length) in lengths.iter().enumerate().take(end).skip(start) {
        write_zeros(out, id, 1, length, input_count)?;
    }

    Ok(())
}

fn write_zeros(
    out: &mut Output,
    id: usize,
    start: usize,
    end: usize,
    input_count: usize,
) -> io::Result<()> {
    let zeros = vec![0; input_count];

    for position in start..=end {
        out.write(id, position, &zeros)?;
    }

    Ok(())
}
