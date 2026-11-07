//! `ntools sort`: sort alignments by coordinate, read name, or tag (samtools sort).

use std::{
    fs::{self, File},
    io,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use pasta::{bgzf::io::writer::CompressionLevel, sam};

use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions, Rec},
    bgzf_stream::BgzfWriter,
    global_args::IndexRequest,
    output::{self, AlignmentWriter, OutputFormat, WriteOptions, split_index_path},
    program,
    runs::{self, Chunk, SpillFile},
    sort_order::{Order, TagTieBreak},
    source::STDIN,
    template_coordinate::{self, HexIds, Libraries},
};

/// samtools' default per-thread memory limit.
const DEFAULT_MEMORY_PER_THREAD: &str = "768M";

/// The smallest memory limit samtools sort accepts: smaller limits spill a file every few
/// records.
const MIN_MEMORY: usize = 1 << 20;

/// The most spill files merged at once, however many file descriptors are free. Each open spill
/// file holds a compressed and a decompressed BGZF block, and more files make a deeper merge
/// heap. At the default memory limit, 256 spill files hold about 200 GB of records, so larger
/// inputs merge in passes.
const MAX_MERGE_FAN_IN: usize = 256;

/// File descriptors left free for everything other than spill files while merging (the output's
/// helpers, reference files, an index written afterwards).
const SPARE_DESCRIPTORS: u64 = 16;

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    #[arg(default_value = STDIN)]
    input: String,

    /// Write output to this file instead of stdout (`FILE##idx##INDEX` also names the index for
    /// --write-index).
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Output format [default: from the -o extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// BGZF compression level for BAM output.
    #[arg(short = 'l', default_value_t = output::default_level())]
    level: u8,

    /// Sort unaligned reads by minimiser (samtools' -M), which ntools does not support.
    #[arg(short = 'M', hide = true)]
    minimiser: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Maximum memory per thread before spilling to temporary files (K/M/G suffixes; at least
    /// 1M).
    #[arg(short = 'm', default_value = DEFAULT_MEMORY_PER_THREAD, value_parser = parse_memory)]
    memory_per_thread: usize,

    /// Sort by read name, comparing numbers naturally.
    #[arg(short = 'n', overrides_with = "template_coordinate")]
    name_natural: bool,

    /// Sort by read name, bytewise (also with -n).
    #[arg(short = 'N', overrides_with = "template_coordinate")]
    name_ascii: bool,

    /// Sort by template-coordinate: the unclipped 5' positions of both ends of the template
    /// (the mate's from its MC tag, as added by fixmate -m), library, CB, MI, then read name.
    /// The last of -n, -N, and this option applies; -t overrides it.
    #[arg(
        long = "template-coordinate",
        overrides_with_all = ["name_natural", "name_ascii"]
    )]
    template_coordinate: bool,

    /// Sort by the value of TAG (records without it first), then by position, or by read name
    /// with -n/-N.
    #[arg(short = 't', value_name = "TAG", value_parser = parse_tag)]
    tag: Option<[u8; 2]>,

    /// Do not use the reverse strand. Accepted for compatibility: in samtools it only affects
    /// minimiser sorting (-M), which ntools does not support.
    #[arg(short = 'R')]
    forward_only: bool,

    /// Squash homopolymers. Accepted for compatibility: in samtools it only affects minimiser
    /// sorting (-M), which ntools does not support.
    #[arg(short = 'H')]
    squash_homopolymers: bool,

    /// Prefix for temporary files, or an existing directory to hold them [default: next to the
    /// output, or in the current directory].
    #[arg(short = 'T')]
    temp_prefix: Option<PathBuf>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    #[command(flatten)]
    output_format_options: crate::global_args::OutputFormatOptions,

    /// Additional threads for sorting, compression, and decompression.
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,

    /// Index the output (coordinate order only): CSI for BAM, or the index named with
    /// `-o FILE##idx##INDEX` (BAI if it ends in .bai); CRAI for CRAM.
    #[arg(long = "write-index")]
    write_index: bool,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.minimiser {
        bail!("sort -M (minimiser sort) is not supported");
    }

    let order = order(&args);

    let (output, index) = match args.output.as_deref() {
        Some(path) => {
            let (output, index) = split_index_path(path);
            (Some(output), index)
        }
        None => (None, None),
    };

    let format = args
        .output_format
        .or_else(|| output.as_deref().and_then(OutputFormat::from_path))
        .unwrap_or(OutputFormat::Bam);

    let index_request = if !args.write_index {
        None
    } else if order == Order::Coordinate {
        Some(IndexRequest::new(output.as_deref(), index, format)?)
    } else {
        eprintln!(
            "[W::bam_sort] Ignoring --write-index as it only works for position sorted files."
        );
        None
    };

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let memory_limit = args.memory_per_thread.saturating_mul(args.threads + 1);
    let temp_dir = create_temp_dir(args.temp_prefix.as_deref(), output.as_deref(), true)?;

    let sorted = SortedRecords::collect(
        &mut reader,
        &header,
        order,
        memory_limit,
        args.threads,
        temp_dir.path(),
        |_| Ok(true),
    )?;

    set_sort_order(&mut header, order);

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { args.level })?,
        threads: args.threads,
        reference: args.reference,
        cram: args.output_format_options.cram()?,
    };

    let mut writer = AlignmentWriter::create(output.as_deref(), &write_options, &header)?;
    writer.write_header(&header)?;
    sorted.write_to(&mut writer, &header)?;
    writer.finish(&header)?;

    if let Some(request) = index_request {
        request.write(args.threads)?;
    }

    Ok(())
}

/// The sort order, as samtools picks it: the last of `-n`/`-N` (bytewise if `-N` is given at
/// all) and `--template-coordinate`, except that `-t` sorts by tag, then by name if `-n`/`-N`
/// came last.
fn order(args: &Args) -> Order {
    let then = name_tie_break(args.name_natural, args.name_ascii);

    match (args.tag, then) {
        (Some(tag), then) => Order::Tag { tag, then },
        (None, TagTieBreak::NameNatural) => Order::NameNatural,
        (None, TagTieBreak::NameAscii) => Order::NameAscii,
        (None, TagTieBreak::Coordinate) if args.template_coordinate => Order::TemplateCoordinate,
        (None, TagTieBreak::Coordinate) => Order::Coordinate,
    }
}

/// The name order of `-n` and `-N` (bytewise if `-N` is given at all), or `Coordinate` for
/// neither.
pub fn name_tie_break(name_natural: bool, name_ascii: bool) -> TagTieBreak {
    if name_ascii {
        TagTieBreak::NameAscii
    } else if name_natural {
        TagTieBreak::NameNatural
    } else {
        TagTieBreak::Coordinate
    }
}

/// Parses a tag name; samtools uses its first two characters.
pub fn parse_tag(s: &str) -> anyhow::Result<[u8; 2]> {
    match *s.as_bytes() {
        [a, b, ..] => Ok([a, b]),
        _ => bail!("invalid tag: {s}"),
    }
}

/// Records read and sorted in memory, with any sorted runs spilled to temporary files.
pub struct SortedRecords<'d> {
    chunk: Chunk,
    /// Sorted runs on disk, in input order, which merging relies on to keep ties stable.
    spill_paths: Vec<PathBuf>,
    /// How many merged spill files have been written, to name the next one.
    merged_spills: usize,
    order: Order,
    threads: usize,
    temp_dir: &'d Path,
}

impl<'d> SortedRecords<'d> {
    /// Reads every record that `keep` accepts, spilling sorted runs to `temp_dir` whenever the
    /// buffered records reach `memory_limit` bytes.
    pub fn collect<K>(
        reader: &mut AlignmentReader,
        header: &sam::Header,
        order: Order,
        memory_limit: usize,
        threads: usize,
        temp_dir: &'d Path,
        mut keep: K,
    ) -> anyhow::Result<Self>
    where
        K: FnMut(&Rec<'_>) -> anyhow::Result<bool>,
    {
        let mut sorted = Self {
            chunk: Chunk::new(order),
            spill_paths: Vec::new(),
            merged_spills: 0,
            order,
            threads,
            temp_dir,
        };

        let mut encoder = RawEncoder::new();

        // Template-coordinate keys depend on the header, so records carry them.
        let libraries =
            (order == Order::TemplateCoordinate).then(|| Libraries::from_header(header));
        let mut keyed = Vec::new();

        reader.for_each(header, |record| {
            if !keep(&record)? {
                return Ok(());
            }

            // Non-BAM records are encoded once so sorting and spilling work on raw bytes.
            let raw = record.encode_raw(header, &mut encoder)?;

            match &libraries {
                None => sorted.chunk.push(raw),
                Some(libraries) => {
                    keyed.clear();
                    let library = libraries.get(raw, HexIds::Ignore)?;
                    template_coordinate::push_keyed(&mut keyed, raw, library)?;
                    sorted.chunk.push(&keyed);
                }
            }

            if sorted.chunk.size() >= memory_limit {
                sorted.spill()?;
            }

            Ok(())
        })?;

        Ok(sorted)
    }

    fn spill(&mut self) -> anyhow::Result<()> {
        let index = self.spill_paths.len();
        let path = runs::spill(
            &mut self.chunk,
            self.order,
            self.temp_dir,
            index,
            self.threads,
        )?;
        self.spill_paths.push(path);
        Ok(())
    }

    /// Writes the records in order. `header` must have the input's reference sequences.
    pub fn write_to(
        self,
        writer: &mut AlignmentWriter,
        header: &sam::Header,
    ) -> anyhow::Result<()> {
        let order = self.order;
        let fan_in = merge_fan_in(free_descriptors());

        self.for_each_sorted(fan_in, |record| {
            writer.write_raw(header, unkeyed(order, record))
        })
    }

    /// Calls `f` with each record in order, merging at most `fan_in` spill files at once.
    fn for_each_sorted<F>(mut self, fan_in: usize, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(&[u8]) -> io::Result<()>,
    {
        let order = self.order;

        if self.spill_paths.is_empty() {
            self.chunk.drain_sorted(order, self.threads, f)?;
            return Ok(());
        }

        if !self.chunk.is_empty() {
            self.spill()?;
        }

        self.reduce_spills(fan_in)?;

        let mut sources = open_spills(&self.spill_paths, self.threads)?;

        runs::merge(&mut sources, order, |record| Ok(f(record)?))
    }

    /// Merges consecutive spill files into new ones until at most `fan_in` remain, so that the
    /// final merge keeps no more than `fan_in` files open. Each merge also opens its output, so
    /// it reads at most `fan_in - 1` files.
    ///
    /// A merged file takes the place of the files it merges, and merges break ties by file
    /// order, so the result is the same as merging every file at once. Only as many files are
    /// merged as needed to get down to `fan_in`.
    fn reduce_spills(&mut self, fan_in: usize) -> anyhow::Result<()> {
        // A merge of fewer than two files into a third would not reduce the count.
        let fan_in = fan_in.max(3);

        while self.spill_paths.len() > fan_in {
            let paths = std::mem::take(&mut self.spill_paths);
            // Merging n files removes n - 1.
            let mut excess = paths.len() - fan_in;
            let mut rest = paths.as_slice();

            while excess > 0 && rest.len() > 1 {
                let n = (excess + 1).min(fan_in - 1).min(rest.len());
                let (group, tail) = rest.split_at(n);
                let path = self.merge_spills(group)?;
                self.spill_paths.push(path);
                excess -= n - 1;
                rest = tail;
            }

            self.spill_paths.extend_from_slice(rest);
        }

        Ok(())
    }

    /// Merges spill files into a new one, removes them, and returns the new file's path.
    fn merge_spills(&mut self, paths: &[PathBuf]) -> anyhow::Result<PathBuf> {
        let path = self
            .temp_dir
            .join(format!("merged-{:04}.bam", self.merged_spills));
        self.merged_spills += 1;

        let mut sources = open_spills(paths, self.threads)?;
        let file = File::create(&path)
            .with_context(|| format!("creating temporary file {}", path.display()))?;
        let mut writer = BgzfWriter::new(file, CompressionLevel::FAST, self.threads);

        runs::merge(&mut sources, self.order, |record| {
            Ok(runs::write_raw_record(&mut writer, record)?)
        })?;

        writer.finish()?;
        drop(sources);

        for path in paths {
            fs::remove_file(path)
                .with_context(|| format!("removing temporary file {}", path.display()))?;
        }

        Ok(path)
    }
}

/// Opens spill files to merge, sharing `threads` decompression workers between them. With more
/// files than threads, each is decompressed on the merging thread.
fn open_spills(paths: &[PathBuf], threads: usize) -> anyhow::Result<Vec<SpillFile>> {
    let workers = threads / paths.len().max(1);

    paths
        .iter()
        .map(|path| SpillFile::open(path, workers))
        .collect::<Result<Vec<_>, _>>()
        .context("opening temporary files")
}

/// How many spill files to merge at once given `free` file descriptors (`None` for no limit):
/// as many as are free, up to [`MAX_MERGE_FAN_IN`].
fn merge_fan_in(free: Option<u64>) -> usize {
    free.map_or(MAX_MERGE_FAN_IN, |free| {
        usize::try_from(free).map_or(MAX_MERGE_FAN_IN, |free| free.min(MAX_MERGE_FAN_IN))
    })
}

/// The file descriptors that may still be opened: the soft `RLIMIT_NOFILE` limit less those open
/// now and [`SPARE_DESCRIPTORS`]. `None` if there is no limit or it cannot be read.
fn free_descriptors() -> Option<u64> {
    let mut limit = libc::rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };

    // SAFETY: `getrlimit` only writes the limits into the struct it is given.
    if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &raw mut limit) } != 0
        || limit.rlim_cur == libc::RLIM_INFINITY
    {
        return None;
    }

    // Listing the open descriptors opens one more, which is counted as spare.
    let open = fs::read_dir("/dev/fd").map_or(0, Iterator::count) as u64;
    Some(limit.rlim_cur.saturating_sub(open + SPARE_DESCRIPTORS))
}

/// A sorted record without the key that [`Order::TemplateCoordinate`] appends.
pub fn unkeyed(order: Order, record: &[u8]) -> &[u8] {
    match order {
        Order::TemplateCoordinate => template_coordinate::record(record),
        _ => record,
    }
}

/// Sets `@HD SO` (and `SS`, and `GO` for template-coordinate) for `order` as samtools does,
/// creating `@HD` if needed, and removes `GO` otherwise.
pub fn set_sort_order(header: &mut sam::Header, order: Order) {
    let (sort_order, subsort_order) = order.header_values();
    let group_order = (order == Order::TemplateCoordinate).then_some("query");
    set_sort_fields(header, sort_order, subsort_order.as_deref(), group_order);
}

/// Sets `@HD SO`, `SS`, and `GO` (removing `SS` and `GO` when they are `None`), creating `@HD`
/// if needed (`sam_hdr_update_sort` and `set_sort_order` in samtools).
pub fn set_sort_fields(
    header: &mut sam::Header,
    sort_order: &str,
    subsort_order: Option<&str>,
    group_order: Option<&str>,
) {
    use pasta::sam::header::record::value::{
        Map,
        map::header::tag::{GROUP_ORDER, SORT_ORDER, SUBSORT_ORDER},
    };

    let hd = header.header_mut().get_or_insert_with(Map::default);
    let fields = hd.other_fields_mut();
    fields.insert(SORT_ORDER, sort_order.into());

    for (tag, value) in [(GROUP_ORDER, group_order), (SUBSORT_ORDER, subsort_order)] {
        match value {
            Some(value) => {
                fields.insert(tag, value.into());
            }
            None => {
                fields.shift_remove(&tag);
            }
        }
    }
}

/// Creates the directory for spill files: in the `-T` prefix's directory, named after its file
/// name, or else next to the output. A prefix ending in `/` names the directory itself, and so,
/// with `existing_dir`, does a prefix that is an existing directory, as samtools sort treats
/// them (samtools collate takes any other prefix as a file name prefix).
pub fn create_temp_dir(
    prefix: Option<&Path>,
    output: Option<&Path>,
    existing_dir: bool,
) -> anyhow::Result<tempfile::TempDir> {
    let (dir, name) = match prefix {
        Some(prefix)
            if prefix.as_os_str().as_encoded_bytes().ends_with(b"/")
                || (existing_dir && prefix.is_dir()) =>
        {
            (Some(prefix), "ntools-sort")
        }
        Some(prefix) => (
            prefix.parent().filter(|p| !p.as_os_str().is_empty()),
            prefix
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("ntools-sort"),
        ),
        None => (
            output
                .and_then(Path::parent)
                .filter(|p| !p.as_os_str().is_empty()),
            "ntools-sort",
        ),
    };

    let dir = dir.unwrap_or_else(|| Path::new("."));

    tempfile::Builder::new()
        .prefix(&format!("{name}."))
        .tempdir_in(dir)
        .with_context(|| format!("creating a temporary directory in {}", dir.display()))
}

/// Parses a memory size with an optional K, M, or G suffix (powers of 1024), as samtools does.
/// Like samtools sort, rejects sizes below [`MIN_MEMORY`].
pub fn parse_memory(s: &str) -> anyhow::Result<usize> {
    let s = s.trim();
    let (digits, multiplier) = match s.char_indices().last() {
        Some((i, 'k' | 'K')) => (&s[..i], 1 << 10),
        Some((i, 'm' | 'M')) => (&s[..i], 1 << 20),
        Some((i, 'g' | 'G')) => (&s[..i], 1 << 30),
        _ => (s, 1),
    };

    let n: usize = digits
        .parse()
        .with_context(|| format!("invalid memory size: {s}"))?;

    let n = n
        .checked_mul(multiplier)
        .with_context(|| format!("memory size overflows: {s}"))?;

    if n < MIN_MEMORY {
        // samtools' wording, which shows sizes over 9K in K.
        let (n, suffix) = if n > 9 << 10 { (n >> 10, "K") } else { (n, "") };
        bail!(
            "-m setting ({n}{suffix} bytes) is less than the minimum required (1M). It should be \
             an integer followed by K, M or G, and much higher when sorting a large file"
        );
    }

    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_memory() -> anyhow::Result<()> {
        assert_eq!(parse_memory("768M")?, 768 << 20);
        assert_eq!(parse_memory("2g")?, 2 << 30);
        assert_eq!(parse_memory("1M")?, 1 << 20);
        assert_eq!(parse_memory("1024K")?, 1 << 20);
        assert_eq!(parse_memory("1048576")?, 1 << 20);
        assert!(parse_memory("1X").is_err());

        // Below samtools' 1M minimum.
        for (s, shown) in [
            ("0", "0"),
            ("1", "1"),
            ("2", "2"),
            ("1023K", "1023K"),
            ("1048575", "1023K"),
        ] {
            let e = parse_memory(s).expect_err(s).to_string();
            assert!(
                e.starts_with(&format!(
                    "-m setting ({shown} bytes) is less than the minimum required (1M)"
                )),
                "{s}: {e}"
            );
        }

        Ok(())
    }

    #[test]
    fn test_create_temp_dir() -> anyhow::Result<()> {
        let base = tempfile::tempdir()?;
        let dir = base.path().join("tdir");
        fs::create_dir(&dir)?;

        // sort: an existing directory holds the spills.
        let temp = create_temp_dir(Some(&dir), None, true)?;
        assert_eq!(temp.path().parent(), Some(dir.as_path()));

        // collate: an existing directory is a file name prefix, as in samtools' bamshuf.
        let temp = create_temp_dir(Some(&dir), None, false)?;
        assert_eq!(temp.path().parent(), Some(base.path()));
        assert!(
            temp.path()
                .file_name()
                .is_some_and(|name| { name.to_string_lossy().starts_with("tdir.") })
        );

        // A trailing slash names the directory for both.
        let mut slash = dir.clone().into_os_string();
        slash.push("/");
        let temp = create_temp_dir(Some(Path::new(&slash)), None, false)?;
        assert_eq!(temp.path().parent(), Some(dir.as_path()));

        // Otherwise the prefix's file name names the directory, beside it.
        let temp = create_temp_dir(Some(&dir.join("pre")), None, true)?;
        assert_eq!(temp.path().parent(), Some(dir.as_path()));
        assert!(
            temp.path()
                .file_name()
                .is_some_and(|name| { name.to_string_lossy().starts_with("pre.") })
        );

        Ok(())
    }

    /// An unmapped-strand raw BAM record without sequence on reference sequence 0 at `pos`.
    fn raw_record(pos: i32, name: &str) -> Vec<u8> {
        let mut record = Vec::new();
        record.extend_from_slice(&0i32.to_le_bytes());
        record.extend_from_slice(&pos.to_le_bytes());
        record.push(u8::try_from(name.len() + 1).expect("short name"));
        record.extend_from_slice(&[0; 23]);
        record.extend_from_slice(name.as_bytes());
        record.push(0);
        record
    }

    #[test]
    fn test_for_each_sorted_merges_in_passes() -> anyhow::Result<()> {
        // Every key is shared by many records, so stability shows in the output.
        let inputs = [
            (
                Order::Coordinate,
                (0..500)
                    .map(|i| raw_record((i * 7) % 13, &format!("r{i}")))
                    .collect::<Vec<_>>(),
            ),
            (
                Order::NameNatural,
                (0..500)
                    .map(|i| raw_record(i, &format!("r{}", (i * 7) % 13)))
                    .collect(),
            ),
        ];

        for (order, records) in &inputs {
            let mut expected = records.clone();
            expected.sort_by(|a, b| order.compare(a, b));

            // 25 spill files.
            for fan_in in [0, 3, 4, 7, 24, 25, MAX_MERGE_FAN_IN] {
                let dir = tempfile::tempdir()?;
                let mut sorted = SortedRecords {
                    chunk: Chunk::new(*order),
                    spill_paths: Vec::new(),
                    merged_spills: 0,
                    order: *order,
                    threads: 0,
                    temp_dir: dir.path(),
                };

                for batch in records.chunks(20) {
                    for record in batch {
                        sorted.chunk.push(record);
                    }

                    sorted.spill()?;
                }

                let mut actual = Vec::new();
                sorted.for_each_sorted(fan_in, |record| {
                    actual.push(record.to_vec());
                    Ok(())
                })?;

                assert!(actual == expected, "{order:?}, fan-in {fan_in}");

                // Merged files replace the files they merge.
                let files = fs::read_dir(dir.path())?.count();
                assert!(
                    files <= fan_in.clamp(3, 25),
                    "{files} files, fan-in {fan_in}"
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_merge_fan_in() {
        assert_eq!(merge_fan_in(None), MAX_MERGE_FAN_IN);
        assert_eq!(merge_fan_in(Some(40)), 40);
        assert_eq!(merge_fan_in(Some(u64::MAX)), MAX_MERGE_FAN_IN);
    }
}
