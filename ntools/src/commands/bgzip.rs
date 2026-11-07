//! `ntools bgzip`: HTSlib's `bgzip`, which compresses files to BGZF, decompresses them, and
//! builds or uses their `.gzi` indexes.
//!
//! Blocks are cut where HTSlib cuts them (65280 bytes, ending at a line end for text, or as an
//! index says for `-g`), so `.gzi` files list the same uncompressed offsets. `-l` levels map onto
//! zflate's as HTSlib maps zlib's onto libdeflate's; without one, level 6, the default for
//! variant output. zflate's bytes differ from libdeflate's, so the output is not HTSlib's byte for
//! byte; it decompresses to the same data.

mod blocks;
mod detect;

use std::{
    fs::{self, File, FileTimes, OpenOptions},
    io::{self, BufRead, BufReader, BufWriter, IsTerminal, Read, Seek, SeekFrom, Write},
    path::Path,
    process::ExitCode,
    time::{Duration, SystemTime, UNIX_EPOCH},
};

use anyhow::{Context, bail};
use flate2::read::MultiGzDecoder;
use pasta::{
    bgzf::{gzi, io::raw, io::writer::CompressionLevel},
    util::io::Rewound,
};

use self::blocks::{
    BLOCK_SIZE, BlockPlan, BlockWriter, Frame, Inflated, Ordered, Sink, frame_isize, inflate,
};
use crate::{
    getopt::{self, LongOption, ShortOption},
    output,
};

/// Arguments, parsed as HTSlib's `bgzip` parses them (see [`parse`]).
#[derive(Debug, clap::Args)]
pub struct Args {
    /// Options and files, as for `bgzip`.
    #[arg(trailing_var_arg = true, allow_hyphen_values = true, num_args = 0..)]
    args: Vec<String>,
}

const USAGE: &str = "Usage:   ntools bgzip [OPTIONS] [FILE] ...
Options:
   -b, --offset INT           decompress at virtual file pointer (0-based uncompressed offset)
   -c, --stdout               write on standard output, keep original files unchanged
   -d, --decompress           decompress
   -f, --force                overwrite files without asking
   -g, --rebgzip              use an index file to bgzip a file
   -h, --help                 give this help
   -i, --index                compress and create BGZF index
   -I, --index-name FILE      name of BGZF index file [file.gz.gzi]
   -k, --keep                 don't delete input files during operation
   -l, --compress-level INT   Compression level to use when compressing; 0 to 9, or -1 for default [-1]
   -o, --output FILE          write to file, keep original files unchanged
   -r, --reindex              (re)index compressed file
   -s, --size INT             decompress INT bytes (uncompressed size)
   -t, --test                 test integrity of compressed file
       --binary               Don't align blocks with text lines
   -@, --threads INT          number of compression threads to use [1]
";

const SHORT: &[ShortOption<'static>] = &[
    ('b', "offset", true),
    ('c', "stdout", false),
    ('d', "decompress", false),
    ('f', "force", false),
    ('g', "rebgzip", false),
    ('h', "help", false),
    ('i', "index", false),
    ('I', "index-name", true),
    ('k', "keep", false),
    ('l', "compress-level", true),
    ('o', "output", true),
    ('r', "reindex", false),
    ('s', "size", true),
    ('t', "test", false),
    ('@', "threads", true),
];

const LONG: &[LongOption<'static>] = &[
    ("help", "help", false),
    ("offset", "offset", true),
    ("stdout", "stdout", false),
    ("decompress", "decompress", false),
    ("force", "force", false),
    ("index", "index", false),
    ("index-name", "index-name", true),
    ("compress-level", "compress-level", true),
    ("reindex", "reindex", false),
    ("rebgzip", "rebgzip", false),
    ("size", "size", true),
    ("threads", "threads", true),
    ("test", "test", false),
    ("version", "version", false),
    ("keep", "keep", false),
    ("binary", "binary", false),
    ("output", "output", true),
];

/// The options, with HTSlib's names and meanings (several options set more than one).
#[derive(Debug)]
struct Options {
    compress: bool,
    compress_level: i64,
    /// `-c`: write to stdout.
    pstdout: bool,
    /// `-f`, counted.
    is_forced: u32,
    test: bool,
    index: bool,
    rebgzip: bool,
    reindex: bool,
    keep: bool,
    binary: bool,
    /// `-b`.
    start: i64,
    /// `-b` plus `-s`, or -1.
    end: i64,
    index_fname: Option<String>,
    write_fname: Option<String>,
    threads: i64,
    files: Vec<String>,
}

enum Parsed {
    Run(Options),
    /// Print this to stdout and succeed (`-h`, `--version`).
    Print(String),
}

/// Applies the options in command-line order, as `bgzip`'s `getopt_long` loop does.
fn parse(args: &[String]) -> anyhow::Result<Parsed> {
    let (opts, files) = getopt::getopt(args, SHORT, LONG)?;

    let mut o = Options {
        compress: true,
        compress_level: -1,
        pstdout: false,
        is_forced: 0,
        test: false,
        index: false,
        rebgzip: false,
        reindex: false,
        keep: false,
        binary: false,
        start: 0,
        end: -1,
        index_fname: None,
        write_fname: None,
        threads: 1,
        files,
    };
    let mut size = -1;

    for (key, value) in opts {
        let v = value.unwrap_or_default();

        match key {
            "decompress" => o.compress = false,
            "stdout" => o.pstdout = true,
            "offset" => {
                o.start = getopt::atol(&v);
                o.compress = false;
                o.pstdout = true;
            }
            "size" => {
                size = getopt::atol(&v);
                o.pstdout = true;
            }
            "force" => o.is_forced += 1,
            "index" => o.index = true,
            "index-name" => o.index_fname = Some(v),
            "compress-level" => o.compress_level = getopt::atol(&v),
            "rebgzip" => o.rebgzip = true,
            "reindex" => {
                o.reindex = true;
                o.compress = false;
            }
            "threads" => o.threads = getopt::atol(&v),
            "test" => {
                o.test = true;
                o.compress = false;
                o.reindex = false;
            }
            "keep" => o.keep = true,
            "output" => o.write_fname = Some(v),
            "binary" => o.binary = true,
            "version" => {
                return Ok(Parsed::Print(format!(
                    "bgzip (ntools) {}\n",
                    env!("CARGO_PKG_VERSION")
                )));
            }
            "help" => return Ok(Parsed::Print(usage())),
            _ => unreachable!("every option is handled"),
        }
    }

    if size >= 0 {
        o.end = o.start + size;
    }

    Ok(Parsed::Run(o))
}

fn usage() -> String {
    format!("\nVersion: ntools {}\n{USAGE}", env!("CARGO_PKG_VERSION"))
}

/// The exit status for the usage printed after an error.
fn usage_error() -> ExitCode {
    eprint!("{}", usage());
    ExitCode::FAILURE
}

/// Runs `bgzip`.
pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    let o = match parse(&args.args) {
        Ok(Parsed::Run(o)) => o,
        Ok(Parsed::Print(text)) => {
            print!("{text}");
            return Ok(ExitCode::SUCCESS);
        }
        Err(e) => {
            eprintln!("bgzip: {e}");
            return Ok(usage_error());
        }
    };

    if o.end >= 0 && o.end < o.start {
        bail!("[bgzip] Illegal region: [{}, {}]", o.start, o.end);
    }

    if (o.index || o.reindex) && o.rebgzip {
        bail!("[bgzip] Can't produce a index and rebgzip simultaneously");
    }

    if o.rebgzip && o.index_fname.is_none() {
        bail!("[bgzip] Index file name expected with rebgzip.  See -I option.");
    }

    if (o.index || o.reindex)
        && o.write_fname.is_none()
        && o.index_fname.is_some()
        && o.files.len() > 1
    {
        bail!("[bgzip] Cannot specify index filename with multiple data file on index, reindex.");
    }

    Bgzip::new(o)?.run()
}

/// The state `bgzip`'s loop over its inputs keeps.
struct Bgzip {
    o: Options,
    /// The worker threads (`-@` above 1), or 0.
    workers: usize,
    /// The exit status: 0, or 2 after a declined overwrite or an unknown extension.
    ret: u8,
    /// The `-o` file, once opened: compressed output shared by every input...
    shared_writer: Option<BlockWriter<output::Sink>>,
    /// ... or decompressed output.
    shared_out: Option<output::Sink>,
}

/// The action for one input.
enum Step {
    Next,
    /// End the command with this status.
    Exit(ExitCode),
}

impl Bgzip {
    fn new(mut o: Options) -> anyhow::Result<Self> {
        if let Some(name) = &o.write_fname {
            if o.pstdout {
                bail!("[bgzip] Cannot write to {name} and stdout at the same time.");
            } else if name.is_empty() || name == "-" {
                // stdout has special handling, so treat it as `-c`.
                o.pstdout = true;
                o.write_fname = None;
            }
        }

        let workers = if o.threads > 1 {
            usize::try_from(o.threads).unwrap_or(1)
        } else {
            0
        };

        Ok(Self {
            o,
            workers,
            ret: 0,
            shared_writer: None,
            shared_out: None,
        })
    }

    fn run(mut self) -> anyhow::Result<ExitCode> {
        // No files means stdin; so does `-`.
        let inputs: Vec<Option<String>> = if self.o.files.is_empty() {
            vec![None]
        } else {
            let files = std::mem::take(&mut self.o.files);
            files.into_iter().map(|f| (f != "-").then_some(f)).collect()
        };

        for input in &inputs {
            let step = if self.o.compress {
                self.compress(input.as_deref())?
            } else if self.o.reindex {
                self.reindex(input.as_deref())?
            } else {
                self.decompress(input.as_deref())?
            };

            if let Step::Exit(code) = step {
                return Ok(code);
            }
        }

        if let Some(writer) = self.shared_writer.take() {
            let (_, index) = writer.finish().context("Output close failed")?;

            if self.o.index {
                let name = self.o.write_fname.as_deref().unwrap_or_default();
                let path = match &self.o.index_fname {
                    Some(path) => path.clone(),
                    None => format!("{name}.gzi"),
                };
                write_index(&path, index)?;
            }
        }

        if let Some(mut out) = self.shared_out.take() {
            out.flush()?;
        }

        Ok(ExitCode::from(self.ret))
    }

    /// Compresses one input (stdin for `None`).
    fn compress(&mut self, input: Option<&str>) -> anyhow::Result<Step> {
        let o = &self.o;

        if !(-1..=9).contains(&o.compress_level) {
            eprintln!("[bgzip] Invalid compress-level: {}", o.compress_level);
            return Ok(Step::Exit(ExitCode::FAILURE));
        }

        let src: Box<dyn Read> = match input {
            Some(path) => match File::open(path) {
                Ok(f) => Box::new(f),
                Err(e) => {
                    eprintln!("[bgzip] {}: {path}", strerror(&e));
                    return Ok(Step::Exit(ExitCode::FAILURE));
                }
            },
            None => Box::new(io::stdin().lock()),
        };

        let level = compression_level(o.compress_level);
        let mut statfilename = None;

        let mut writer = if o.write_fname.is_some() {
            match self.shared_writer.take() {
                Some(writer) => writer,
                None => {
                    let name = o.write_fname.as_deref().unwrap_or_default();
                    match File::create(name) {
                        Ok(f) => BlockWriter::new(file_sink(f), level, self.workers),
                        Err(e) => {
                            eprintln!("[bgzip] can't create {name}: {}", strerror(&e));
                            return Ok(Step::Exit(ExitCode::FAILURE));
                        }
                    }
                }
            }
        } else if let Some(path) = input {
            if o.pstdout {
                BlockWriter::new(output::raw_stdout()?, level, self.workers)
            } else {
                let name = format!("{path}.gz");

                let file = match create(&name, o.is_forced > 0) {
                    Created::File(f) => f,
                    Created::Declined => {
                        self.ret = 2;
                        return Ok(Step::Next);
                    }
                    Created::Failed(e) => {
                        eprintln!("[bgzip] can't create {name}: {}", strerror(&e));
                        return Ok(Step::Exit(ExitCode::FAILURE));
                    }
                };

                statfilename = Some(name);
                BlockWriter::new(file_sink(file), level, self.workers)
            }
        } else if !o.pstdout && io::stdout().is_terminal() {
            return Ok(Step::Exit(usage_error()));
        } else if o.index && o.index_fname.is_none() {
            eprintln!("[bgzip] Index file name expected when writing to stdout");
            return Ok(Step::Exit(ExitCode::FAILURE));
        } else {
            BlockWriter::new(output::raw_stdout()?, level, self.workers)
        };

        if o.rebgzip {
            let name = o.index_fname.as_deref().unwrap_or_default();
            let index =
                gzi::fs::read(name).with_context(|| format!("Could not load index: {name}"))?;
            let mut plan = BlockPlan::new(index.as_ref());
            feed(src, |data| writer.write_planned(data, &mut plan))?;
        } else {
            compress_stream(src, &mut writer, o.binary)?;
        }

        if o.write_fname.is_some() {
            self.shared_writer = Some(writer);
        } else {
            let (_, index) = writer.finish().context("Output close failed")?;

            if o.index {
                let path = match (&o.index_fname, input) {
                    (Some(path), _) => path.clone(),
                    (None, Some(path)) => format!("{path}.gz.gzi"),
                    (None, None) => bail!(
                        "Can not write index for stdin data without index filename, use -I \
                         option to set index file."
                    ),
                };
                write_index(&path, index)?;
            }
        }

        if let (Some(name), Some(path)) = (&statfilename, input) {
            copy_times(path, name);
        }

        if let Some(path) = input
            && !o.pstdout
            && !o.keep
            && o.write_fname.is_none()
        {
            fs::remove_file(path).ok();
        }

        Ok(Step::Next)
    }

    /// Decompresses (or tests, with `-t`) one input.
    fn decompress(&mut self, input: Option<&str>) -> anyhow::Result<Step> {
        let file = match input {
            Some(path) => match File::open(path) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("[bgzip] Could not open {path}: {}", strerror(&e));
                    return Ok(Step::Exit(ExitCode::FAILURE));
                }
            },
            None if !self.o.pstdout && io::stdin().is_terminal() => {
                return Ok(Step::Exit(usage_error()));
            }
            None => stdin_file()?,
        };

        let mut src = Source::open(file)?;

        if src.compression == Compression::None {
            match input {
                Some(path) => eprintln!("[bgzip] {path}: not a compressed file -- ignored"),
                None => eprintln!("[bgzip] stdin is not compressed -- ignored"),
            }
            return Ok(Step::Exit(ExitCode::FAILURE));
        }

        let mut statfilename = None;

        let out: output::Sink = if self.o.test {
            Box::new(io::sink())
        } else if let Some(out) = self.shared_out.take() {
            out
        } else if let Some(path) = input
            && !self.o.pstdout
        {
            let mut is_forced = self.o.is_forced;

            let name = match check_name_and_extension(path, &mut is_forced) {
                Ok(name) => name,
                Err(NameError::NoExtension) => return Ok(Step::Exit(ExitCode::FAILURE)),
                Err(NameError::Declined) => {
                    self.ret = 2;
                    return Ok(Step::Next);
                }
            };

            // An `-o` file is opened once and overwritten without asking.
            let (target, forced) = match &self.o.write_fname {
                Some(w) => (w.as_str(), true),
                None => (name.as_str(), is_forced > 0),
            };

            let out = match create(target, forced) {
                Created::File(f) => file_sink(f),
                Created::Declined => {
                    self.ret = 2;
                    return Ok(Step::Next);
                }
                Created::Failed(e) => {
                    eprintln!("[bgzip] can't create {name}: {}", strerror(&e));
                    return Ok(Step::Exit(ExitCode::FAILURE));
                }
            };

            statfilename = Some(name);
            out
        } else if let (None, Some(name)) = (input, &self.o.write_fname) {
            match File::create(name) {
                Ok(f) => file_sink(f),
                Err(e) => {
                    eprintln!("[bgzip] can't create {name}: {}", strerror(&e));
                    return Ok(Step::Exit(ExitCode::FAILURE));
                }
            }
        } else {
            output::raw_stdout()?
        };

        let mut out = self.read_into(&mut src, input, out)?;

        if self.o.write_fname.is_some() && !self.o.test {
            self.shared_out = Some(out);
        } else {
            out.flush()?;
        }

        if let (Some(name), Some(path)) = (&statfilename, input)
            && self.o.write_fname.is_none()
        {
            copy_times(path, name);
        }

        if let Some(path) = input
            && !self.o.pstdout
            && !self.o.test
            && !self.o.keep
            && self.o.write_fname.is_none()
        {
            fs::remove_file(path).ok();
        }

        Ok(Step::Next)
    }

    /// Decompresses `src` from `-b` for `-s` bytes (or to its end) into `out`, and returns it.
    fn read_into(
        &self,
        src: &mut Source,
        input: Option<&str>,
        out: output::Sink,
    ) -> anyhow::Result<output::Sink> {
        let o = &self.o;
        let start = u64::try_from(o.start).unwrap_or(0);
        let mut skip = 0;

        if start > 0 {
            let index = match (&o.index_fname, input) {
                (Some(name), _) => {
                    gzi::fs::read(name).with_context(|| format!("Could not load index: {name}"))?
                }
                (None, Some(path)) => gzi::fs::read(format!("{path}.gzi"))
                    .with_context(|| format!("Could not load index: {path}.gzi"))?,
                (None, None) => bail!(
                    "The -b option requires -I when reading from stdin (and stdin must be \
                     seekable)"
                ),
            };

            skip = src
                .seek_uncompressed(index.as_ref(), start)
                .with_context(|| format!("Could not seek to {start}-th (uncompressd) byte"))?;
        }

        let limit = u64::try_from(o.end)
            .ok()
            .map(|end| end.saturating_sub(start));
        let sink = DataSink {
            out,
            skip,
            remaining: limit,
        };

        let sink = match src.compression {
            Compression::Bgzf => read_bgzf(src, sink, self.workers)?,
            Compression::Gzip => read_gzip(src, sink)?,
            Compression::None => unreachable!("uncompressed input is ignored"),
        };

        Ok(sink.out)
    }

    /// Rebuilds the `.gzi` index of one input.
    fn reindex(&mut self, input: Option<&str>) -> anyhow::Result<Step> {
        let file = match input {
            Some(path) => {
                File::open(path).with_context(|| format!("[bgzip] Could not open file: {path}"))?
            }
            None => {
                if self.o.index_fname.is_none() {
                    bail!("[bgzip] Index file name expected when reading from stdin");
                }
                stdin_file()?
            }
        };

        let mut src = Source::open(file)?;

        let index = match src.compression {
            // HTSlib reads uncompressed data through, with no blocks to index.
            Compression::None => Vec::new(),
            Compression::Gzip => {
                bail!("Is the file gzipped or bgzipped? The latter is required for indexing.")
            }
            Compression::Bgzf => index_bgzf(&mut src, self.workers)
                .context("Is the file gzipped or bgzipped? The latter is required for indexing.")?,
        };

        let path = match (&self.o.index_fname, input) {
            (Some(path), _) => path.clone(),
            (None, Some(path)) => format!("{path}.gzi"),
            (None, None) => unreachable!("stdin needs an index name"),
        };

        write_index(&path, index)?;

        Ok(Step::Next)
    }
}

/// The compression level for `-l`: each zlib level 0 to 9 as HTSlib maps it onto libdeflate's
/// (`bgzf_compress`), which zflate keeps; -1, the default for the VCF and other text `bgzip`
/// mostly compresses, [`CompressionLevel::VARIANTS`] (level 6).
fn compression_level(level: i64) -> CompressionLevel {
    let Ok(level) = usize::try_from(level) else {
        return CompressionLevel::VARIANTS;
    };
    let level = [0, 1, 2, 3, 5, 6, 7, 8, 10, 12][level];

    CompressionLevel::new(level).expect("levels 0 to 12 are valid")
}

/// Compresses `src` as `bgzip` does without `-g`: text in blocks that end at line ends (and
/// with the header lines in blocks of their own), anything else in full blocks.
fn compress_stream<W>(
    src: Box<dyn Read>,
    writer: &mut BlockWriter<W>,
    binary: bool,
) -> io::Result<()>
where
    W: Write + Send + 'static,
{
    let mut src = src;
    let mut prefix = vec![0; detect::PEEK_LEN];
    let n = read_full(&mut src, &mut prefix)?;
    prefix.truncate(n);

    let textual = !binary && detect::is_textual(&prefix);
    let src = Rewound::new(prefix, src);

    if !textual {
        return feed(src, |data| writer.write(data));
    }

    let mut src = src;
    let mut buffer = vec![0; BLOCK_SIZE];
    let mut n = 0;
    let mut in_header = true;
    let mut long_line = false;

    loop {
        let c = read_full(&mut src, &mut buffer[n..])?;

        if c == 0 {
            break;
        }

        let c2 = c + n;
        let mut flush = false;

        if in_header && (long_line || buffer[0] == b'@' || buffer[0] == b'#') {
            // Scan forward to find the last header line.
            let mut last_start = 0;
            n = 0;

            while n < c2 {
                n += 1;

                if buffer[n - 1] != b'\n' {
                    continue;
                }

                last_start = n;

                if n < c2 && !(buffer[n] == b'@' || buffer[n] == b'#') {
                    in_header = false;
                    break;
                }
            }

            if last_start == 0 {
                n = c2;
                long_line = true;
            } else {
                n = last_start;
                flush = true;
                long_line = false;
            }
        } else {
            // Scan backwards to find the last newline.
            match memchr::memrchr(b'\n', &buffer[..c2]) {
                Some(i) => {
                    flush = true;
                    n = i + 1;
                }
                None => n = c2,
            }
        }

        // `n` is either the end of the buffer (with no flush) or the first byte after a newline
        // and a flush point.
        writer.write(&buffer[..n])?;

        if flush {
            writer.flush()?;
        }

        buffer.copy_within(n..c2, 0);
        n = c2 - n;
    }

    // Trailing data.
    writer.write(&buffer[..n])
}

/// Reads `src` in [`BLOCK_SIZE`] chunks into `write`.
fn feed<R, F>(mut src: R, mut write: F) -> io::Result<()>
where
    R: Read,
    F: FnMut(&[u8]) -> io::Result<()>,
{
    let mut buffer = vec![0; BLOCK_SIZE];

    loop {
        match read_full(&mut src, &mut buffer)? {
            0 => return Ok(()),
            n => write(&buffer[..n])?,
        }
    }
}

/// Reads until `buf` is full or the stream ends, as HTSlib's `hread` does, and returns the
/// number of bytes read.
fn read_full<R: Read + ?Sized>(src: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut len = 0;

    while len < buf.len() {
        match src.read(&mut buf[len..]) {
            Ok(0) => break,
            Ok(n) => len += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    Ok(len)
}

/// How an input is compressed, as `bgzf_compression` says.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Compression {
    None,
    Gzip,
    Bgzf,
}

/// An input opened for decompression.
struct Source {
    reader: Rewound<BufReader<File>>,
    compression: Compression,
    /// The file offset of the next frame.
    offset: u64,
    /// A frame already read (by a seek), to decompress first.
    pending: Option<Vec<u8>>,
}

impl Source {
    fn open(file: File) -> io::Result<Self> {
        let mut reader = BufReader::with_capacity(1 << 17, file);
        let mut magic = vec![0; 18];
        let n = read_full(&mut reader, &mut magic)?;
        magic.truncate(n);

        // `bgzf_read_init`: compressed if it has a full gzip header; BGZF if that has a `BC`
        // extra subfield.
        let compression = if n == 18 && magic.starts_with(&[0x1f, 0x8b]) {
            if magic[3] & 4 != 0 && magic[12..16] == *b"BC\x02\x00" {
                Compression::Bgzf
            } else {
                Compression::Gzip
            }
        } else {
            Compression::None
        };

        Ok(Self {
            reader: Rewound::new(magic, reader),
            compression,
            offset: 0,
            pending: None,
        })
    }

    /// `bgzf_useek`: moves to the block holding uncompressed offset `pos`, by a `.gzi` index, and
    /// returns how far into the block `pos` is.
    fn seek_uncompressed(&mut self, index: &[(u64, u64)], pos: u64) -> io::Result<usize> {
        if self.compression != Compression::Bgzf {
            return Err(io::Error::other("only BGZF files can seek"));
        }

        let i = index.partition_point(|&(_, u)| u <= pos);
        let (coffset, uoffset) = if i == 0 { (0, 0) } else { index[i - 1] };

        self.reader.seek(SeekFrom::Start(coffset))?;
        self.offset = coffset;

        let mut frame = Vec::new();
        let len = if raw::read_frame(&mut self.reader, &mut frame)? {
            frame_isize(&frame)
        } else {
            0
        };

        let skip = usize::try_from(pos - uoffset).unwrap_or(usize::MAX);

        if skip > len {
            return Err(io::Error::other("the offset is outside the data"));
        }

        if !frame.is_empty() {
            self.pending = Some(frame);
        }

        Ok(skip)
    }
}

/// Writes decompressed data, skipping the start of the first block (`-b`) and stopping after
/// `-s` bytes.
struct DataSink {
    out: output::Sink,
    skip: usize,
    remaining: Option<u64>,
}

impl DataSink {
    fn write(&mut self, data: &[u8]) -> io::Result<()> {
        let data = &data[self.skip.min(data.len())..];
        self.skip = 0;

        let data = match &mut self.remaining {
            Some(remaining) => {
                let n = data
                    .len()
                    .min(usize::try_from(*remaining).unwrap_or(usize::MAX));
                *remaining -= n as u64;
                &data[..n]
            }
            None => data,
        };

        self.out.write_all(data)
    }

    fn is_done(&self) -> bool {
        self.remaining == Some(0)
    }
}

impl Sink<Inflated> for DataSink {
    fn accept(&mut self, result: Inflated) -> io::Result<()> {
        match result {
            Ok(data) => self.write(&data),
            Err((offset, e)) => Err(block_error(offset, &e)),
        }
    }
}

fn block_error(offset: u64, e: &io::Error) -> io::Error {
    io::Error::new(
        e.kind(),
        format!("Error in block starting at offset {offset}({offset:X}): {e}"),
    )
}

/// Decompresses BGZF frames on the worker threads, warning (as HTSlib does) when the data ends
/// without the end-of-file marker block.
fn read_bgzf(src: &mut Source, sink: DataSink, workers: usize) -> anyhow::Result<DataSink> {
    // Until the data covers `-s`: what the frames read so far hold past `-b`.
    let mut wanted = sink.remaining;
    let mut skip = sink.skip as u64;
    let mut pipeline = Ordered::new(workers, inflate, sink);
    let mut last_empty = false;
    let mut hit_eof = false;

    let mut frame = src.pending.take();

    loop {
        if wanted == Some(0) {
            break;
        }

        let frame = match frame.take() {
            Some(frame) => frame,
            None => {
                let mut buf = Vec::new();

                match raw::read_frame(&mut src.reader, &mut buf) {
                    Ok(true) => buf,
                    Ok(false) => {
                        hit_eof = true;
                        break;
                    }
                    Err(e) => {
                        drop(pipeline);
                        return Err(block_error(src.offset, &e).into());
                    }
                }
            }
        };

        let len = frame_isize(&frame) as u64;
        last_empty = len == 0;

        if let Some(w) = &mut wanted {
            let n = len.saturating_sub(skip);
            *w = w.saturating_sub(n);
            skip = skip.saturating_sub(len);
        }

        let offset = src.offset;
        src.offset += frame.len() as u64;
        pipeline.push((offset, frame))?;
    }

    let sink = pipeline.finish()?;

    if hit_eof && !last_empty {
        eprintln!("[W::bgzf_read_block] EOF marker is absent. The input may be truncated");
    }

    Ok(sink)
}

/// Decompresses gzip that is not BGZF.
fn read_gzip(src: &mut Source, mut sink: DataSink) -> anyhow::Result<DataSink> {
    let mut decoder = MultiGzDecoder::new(&mut src.reader);
    let mut buffer = vec![0; BLOCK_SIZE];

    while !sink.is_done() {
        let n = read_full(&mut decoder, &mut buffer).context("Error in gzip data")?;

        if n == 0 {
            break;
        }

        sink.write(&buffer[..n])?;
    }

    Ok(sink)
}

/// A decompressed block's file offset and uncompressed size.
type Checked = Result<(u64, usize), (u64, io::Error)>;

fn check(frame: Frame) -> Checked {
    let offset = frame.0;
    inflate(frame).map(|data| (offset, data.len()))
}

/// Builds a `.gzi` index while reading: the offsets of each non-empty block after the first.
#[derive(Default)]
struct IndexSink {
    started: bool,
    uncompressed_offset: u64,
    index: Vec<(u64, u64)>,
}

impl Sink<Checked> for IndexSink {
    fn accept(&mut self, result: Checked) -> io::Result<()> {
        let (offset, len) = result.map_err(|(offset, e)| block_error(offset, &e))?;

        if len > 0 {
            if self.started {
                self.index.push((offset, self.uncompressed_offset));
            }

            self.started = true;
            self.uncompressed_offset += len as u64;
        }

        Ok(())
    }
}

fn index_bgzf(src: &mut Source, workers: usize) -> io::Result<Vec<(u64, u64)>> {
    let mut pipeline = Ordered::new(workers, check, IndexSink::default());
    let mut last_empty = false;

    loop {
        let mut frame = Vec::new();

        match raw::read_frame(&mut src.reader, &mut frame) {
            Ok(true) => {}
            Ok(false) => break,
            Err(e) => return Err(block_error(src.offset, &e)),
        }

        last_empty = frame_isize(&frame) == 0;
        let offset = src.offset;
        src.offset += frame.len() as u64;
        pipeline.push((offset, frame))?;
    }

    let sink = pipeline.finish()?;

    if !last_empty {
        eprintln!("[W::bgzf_read_block] EOF marker is absent. The input may be truncated");
    }

    Ok(sink.index)
}

fn write_index(path: &str, index: Vec<(u64, u64)>) -> anyhow::Result<()> {
    // Flushed explicitly, as dropping a `BufWriter` (as `gzi::fs::write` does) ignores errors.
    let write = || -> io::Result<()> {
        let mut out = BufWriter::new(File::create(path)?);
        gzi::io::Writer::new(&mut out).write_index(&index.into())?;
        out.flush()
    };

    write().with_context(|| format!("Could not write index to '{path}'"))
}

/// The result of creating an output file.
enum Created {
    File(File),
    /// It exists, and overwriting it was declined (or could not be asked).
    Declined,
    Failed(io::Error),
}

/// Creates `name`, asking before overwriting it unless `forced` (`bgzf_open(name, "wx")`, or
/// `O_EXCL`, then `confirm_overwrite`).
fn create(name: &str, forced: bool) -> Created {
    let result = if forced {
        File::create(name)
    } else {
        OpenOptions::new().write(true).create_new(true).open(name)
    };

    match result {
        Ok(f) => Created::File(f),
        Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {
            let question =
                format!("[bgzip] {name} already exists; do you wish to overwrite (y or n)? ");

            if confirm(&question) {
                match File::create(name) {
                    Ok(f) => Created::File(f),
                    Err(e) => Created::Failed(e),
                }
            } else {
                Created::Declined
            }
        }
        Err(e) => Created::Failed(e),
    }
}

/// Asks a yes-or-no question on a terminal; anywhere else the answer is no.
fn confirm(question: &str) -> bool {
    if !io::stdin().is_terminal() {
        return false;
    }

    eprint!("{question}");
    let mut line = String::new();

    io::stdin().lock().read_line(&mut line).is_ok() && line.starts_with(['y', 'Y'])
}

enum NameError {
    /// Exit with status 1.
    NoExtension,
    /// Skip the file; exit with status 2.
    Declined,
}

/// `check_name_and_extension`: the name to decompress `path` to (without its extension), if its
/// extension is `gz`, `bgz`, or `bgzf`, or decompressing it is forced (each `-f` forces one
/// file) or confirmed.
fn check_name_and_extension(path: &str, is_forced: &mut u32) -> Result<String, NameError> {
    let ext_start = path
        .rfind(['.', '/'])
        .filter(|&i| i > 0 && path.as_bytes()[i] == b'.');

    let Some(dot) = ext_start else {
        eprintln!("[bgzip] can't find an extension in {path} -- please rename");
        return Err(NameError::NoExtension);
    };

    let (name, ext) = (&path[..dot], &path[dot + 1..]);
    let known = ["gz", "bgz", "bgzf"]
        .iter()
        .any(|k| k.eq_ignore_ascii_case(ext));

    let mut confirmed = || {
        if *is_forced > 0 {
            *is_forced -= 1;
            return true;
        }

        confirm(&format!(
            "[bgzip] .{ext} is not a known extension; do you wish to decompress to {name} (y or n)? "
        ))
    };

    if known || confirmed() {
        Ok(name.to_owned())
    } else {
        eprintln!("[bgzip] unknown extension .{ext} -- declining to decompress to {name}");
        Err(NameError::Declined)
    }
}

/// Gives `dst` the access and modification times of `src`, in whole seconds, as HTSlib does.
fn copy_times(src: &str, dst: &str) {
    fn seconds(t: io::Result<SystemTime>) -> io::Result<SystemTime> {
        let d = t?.duration_since(UNIX_EPOCH).map_err(io::Error::other)?;
        Ok(UNIX_EPOCH + Duration::from_secs(d.as_secs()))
    }

    let Ok(metadata) = fs::metadata(src) else {
        eprintln!("[bgzip] Failed to get file specification.");
        return;
    };

    let result = (|| {
        let times = FileTimes::new()
            .set_accessed(seconds(metadata.accessed())?)
            .set_modified(seconds(metadata.modified())?);
        OpenOptions::new()
            .write(true)
            .open(Path::new(dst))?
            .set_times(times)
    })();

    if result.is_err() {
        eprintln!("[bgzip] Failed to set file specifications.");
        eprintln!("[bgzip] Failed to set file specification.");
    }
}

fn file_sink(file: File) -> output::Sink {
    Box::new(BufWriter::with_capacity(1 << 17, file))
}

/// stdin as a file, which can seek if it is one.
fn stdin_file() -> anyhow::Result<File> {
    use std::os::fd::AsFd;

    let fd = io::stdin()
        .as_fd()
        .try_clone_to_owned()
        .context("duplicating stdin")?;

    Ok(File::from(fd))
}

/// The C library's message for an I/O error (`strerror(errno)`), without Rust's suffix.
fn strerror(e: &io::Error) -> String {
    let message = e.to_string();

    match message.find(" (os error") {
        Some(i) => message[..i].to_owned(),
        None => message,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    fn options(s: &str) -> anyhow::Result<Options> {
        match parse(&args(s))? {
            Parsed::Run(o) => Ok(o),
            Parsed::Print(_) => bail!("expected options"),
        }
    }

    #[test]
    fn test_parse() -> anyhow::Result<()> {
        let o = options("-dc -@4 a.gz -l9 b.gz")?;
        assert!(!o.compress && o.pstdout);
        assert_eq!((o.threads, o.compress_level), (4, 9));
        assert_eq!(o.files, ["a.gz", "b.gz"]);

        // -b decompresses to stdout; -s alone only writes to stdout.
        let o = options("-b 10 -s 5 x.gz")?;
        assert!(!o.compress && o.pstdout);
        assert_eq!((o.start, o.end), (10, 15));
        let o = options("-s 5 x")?;
        assert!(o.compress && o.pstdout);

        // -t after -r tests; -r after -t reindexes.
        let o = options("-r -t x.gz")?;
        assert!(o.test && !o.reindex);
        let o = options("-t -r x.gz")?;
        assert!(o.test && o.reindex);

        assert!(matches!(parse(&args("-c -h x"))?, Parsed::Print(_)));
        assert!(parse(&args("-x")).is_err());

        Ok(())
    }

    #[test]
    fn test_check_name_and_extension() {
        let mut forced = 0;
        assert!(
            matches!(check_name_and_extension("a/b.vcf.gz", &mut forced), Ok(n) if n == "a/b.vcf")
        );
        assert!(matches!(check_name_and_extension("x.BGZF", &mut forced), Ok(n) if n == "x"));
        assert!(matches!(
            check_name_and_extension("noext", &mut forced),
            Err(NameError::NoExtension)
        ));
        assert!(matches!(
            check_name_and_extension("dir.d/f", &mut forced),
            Err(NameError::NoExtension)
        ));
        assert!(matches!(
            check_name_and_extension(".gz", &mut forced),
            Err(NameError::NoExtension)
        ));

        let mut forced = 1;
        assert!(matches!(check_name_and_extension("z.foo", &mut forced), Ok(n) if n == "z"));
        assert_eq!(forced, 0);
    }

    #[test]
    fn test_compression_level() {
        assert_eq!(compression_level(0), CompressionLevel::NONE);

        assert_eq!(compression_level(-1), CompressionLevel::VARIANTS);
        assert_eq!(compression_level(-1).get(), 6);
        assert_eq!(compression_level(4).get(), 5);
        assert_eq!(compression_level(6).get(), 7);
        assert_eq!(compression_level(9).get(), 12);
    }

    fn compress_text(text: &[u8]) -> io::Result<Vec<usize>> {
        let mut writer = BlockWriter::new(Vec::new(), CompressionLevel::FAST, 0);
        compress_stream(Box::new(io::Cursor::new(text.to_vec())), &mut writer, false)?;
        let (out, _) = writer.finish()?;

        let mut src = &out[..];
        let mut frame = Vec::new();
        let mut sizes = Vec::new();

        while raw::read_frame(&mut src, &mut frame)? {
            sizes.push(frame_isize(&frame));
        }

        Ok(sizes)
    }

    #[test]
    fn test_compress_stream() -> io::Result<()> {
        // Header lines go in a block of their own, then blocks end at line ends.
        let mut text = b"#h1\n#h2\n".to_vec();
        let line = b"chr1\t100\t200\n";
        while text.len() < 150_000 {
            text.extend_from_slice(line);
        }

        let sizes = compress_text(&text)?;
        assert_eq!(sizes[0], 8);
        assert!(
            sizes[1..sizes.len() - 2]
                .iter()
                .all(|n| n % line.len() == 0)
        );
        assert_eq!(sizes.iter().sum::<usize>(), text.len());

        // Binary data fills blocks.
        let binary: Vec<u8> = (0..150_000u32).map(|i| (i % 7) as u8).collect();
        assert_eq!(compress_text(&binary)?, [65280, 65280, 19440, 0]);

        Ok(())
    }
}
