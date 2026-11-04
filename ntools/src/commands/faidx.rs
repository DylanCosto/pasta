//! `ntools faidx` and `ntools fqidx`: index FASTA or FASTQ files and extract regions from them
//! (samtools faidx and fqidx).
//!
//! Indexes are built as HTSlib's `fai_build_core` builds them (the same `.fai`, and a `.gzi` for
//! BGZF-compressed input), and regions are parsed as `hts_parse_region` parses them.

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    process::ExitCode,
    sync::Mutex,
};

use anyhow::{Context, bail};
use pasta::bgzf;
use pasta::fasta::fai;

use crate::global_args::{OutputFormatOptions, Verbosity};
use crate::{hts_region, sequence::reverse_complement};

/// samtools faidx's BGZF compression level for compressed output.
const DEFAULT_BGZF_LEVEL: u8 = 4;

/// samtools' default line length (`DEFAULT_FASTA_LINE_LEN`).
const DEFAULT_LINE_LENGTH: i64 = 60;

/// The length of a BGZF block's header (with its `BC` subfield), as in `BLOCK_HEADER_LENGTH`.
const BGZF_HEADER_LEN: usize = 18;

/// The length of a BGZF block's footer (CRC32 and ISIZE), as in `BLOCK_FOOTER_LENGTH`.
const BGZF_FOOTER_LEN: usize = 8;

/// The format of the sequences in a file and its index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequenceFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// FASTA or FASTQ file (plain or BGZF-compressed).
    input: PathBuf,

    /// Regions to extract (`chr1`, `chr1:100-200`, `{name:with:colons}:1-10`). With none (and no
    /// -r), builds the index.
    regions: Vec<String>,

    /// Write output to this file (BGZF-compressed for `.gz`, `.bgz`, or `.bgzf`), or with no
    /// regions, the index.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Line length of extracted sequences (0 for one line) [default: the input's].
    #[arg(short = 'n', long = "length", allow_negative_numbers = true)]
    length: Option<i64>,

    /// Continue after a missing region.
    #[arg(short = 'c', long = "continue")]
    continue_on_missing: bool,

    /// Read regions from this file, one per line.
    #[arg(short = 'r', long = "region-file")]
    region_file: Option<PathBuf>,

    /// Reverse complement the extracted sequences.
    #[arg(short = 'i', long = "reverse-complement")]
    reverse_complement: bool,

    /// Strand marks added to names: rc (/rc on the reverse strand), no, sign ((+) and (-)), or
    /// custom,<pos>,<neg>.
    #[arg(long = "mark-strand", default_value = "rc")]
    mark_strand: String,

    /// The index file [default: <input>.fai].
    #[arg(long = "fai-idx")]
    fai_idx: Option<PathBuf>,

    /// The BGZF index file [default: <input>.gzi].
    #[arg(long = "gzi-idx")]
    gzi_idx: Option<PathBuf>,

    /// The file and its index are FASTQ (faidx only; fqidx always reads FASTQ).
    #[arg(short = 'f', long = "fastq")]
    fastq: bool,

    /// Also index the output file.
    #[arg(long = "write-index")]
    write_index: bool,

    /// Additional threads (accepted for compatibility).
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,

    /// Print help.
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,

    #[command(flatten)]
    output_format_options: OutputFormatOptions,

    #[command(flatten)]
    verbosity: Verbosity,
}

/// One indexed sequence (`faidx1_t`).
#[derive(Clone, Debug, Eq, PartialEq)]
struct Entry {
    name: Vec<u8>,
    len: u64,
    seq_offset: u64,
    line_blen: u64,
    line_len: u64,
    qual_offset: u64,
}

/// A `.fai` index.
#[derive(Debug, Default)]
struct Index {
    entries: Vec<Entry>,
    by_name: HashMap<Vec<u8>, usize>,
}

impl Index {
    /// `fai_insert_index`: ignores (with a warning) a name seen before.
    fn insert(&mut self, entry: Entry) -> anyhow::Result<()> {
        if entry.line_len < entry.line_blen {
            bail!(
                "Malformed line: width ({}) less than base count ({})",
                entry.line_len,
                entry.line_blen
            );
        }

        if self.by_name.contains_key(&entry.name) {
            eprintln!(
                "[W::fai_insert_index] Ignoring duplicate sequence \"{}\" at byte offset {}",
                String::from_utf8_lossy(&entry.name),
                entry.seq_offset
            );
            return Ok(());
        }

        self.by_name.insert(entry.name.clone(), self.entries.len());
        self.entries.push(entry);
        Ok(())
    }

    /// `fai_save`
    fn write<W: Write>(&self, format: SequenceFormat, out: &mut W) -> io::Result<()> {
        for e in &self.entries {
            out.write_all(&e.name)?;

            match format {
                SequenceFormat::Fasta => writeln!(
                    out,
                    "\t{}\t{}\t{}\t{}",
                    e.len, e.seq_offset, e.line_blen, e.line_len
                )?,
                SequenceFormat::Fastq => writeln!(
                    out,
                    "\t{}\t{}\t{}\t{}\t{}",
                    e.len, e.seq_offset, e.line_blen, e.line_len, e.qual_offset
                )?,
            }
        }

        Ok(())
    }

    /// `fai_read`: FASTA indexes have four numbers after the name, FASTQ ones five.
    fn read(path: &Path, format: SequenceFormat) -> anyhow::Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut index = Self::default();

        for (n, line) in BufReader::new(file).split(b'\n').enumerate() {
            let line = line?;
            let name_end = line
                .iter()
                .position(|b| b.is_ascii_whitespace())
                .unwrap_or(line.len());
            let name = line[..name_end].to_vec();

            let numbers: Vec<u64> = line[name_end..]
                .split(|b| b.is_ascii_whitespace())
                .filter(|t| !t.is_empty())
                .map_while(|t| std::str::from_utf8(t).ok()?.parse().ok())
                .collect();

            let needed = match format {
                SequenceFormat::Fasta => 4,
                SequenceFormat::Fastq => 5,
            };

            if numbers.len() < needed {
                if format == SequenceFormat::Fastq && numbers.len() == 4 {
                    bail!(
                        "Possibly this is a FASTA index, try using faidx.  Problem in {} line {}",
                        path.display(),
                        n + 1
                    );
                }

                bail!(
                    "Could not understand {format:?} index {} line {}",
                    path.display(),
                    n + 1
                );
            }

            index.insert(Entry {
                name,
                len: numbers[0],
                seq_offset: numbers[1],
                line_blen: numbers[2],
                line_len: numbers[3],
                qual_offset: numbers.get(4).copied().filter(|_| needed == 5).unwrap_or(0),
            })?;
        }

        Ok(index)
    }
}

/// A sequence file's bytes, plain or BGZF compressed (with block positions for `.gzi`).
pub(crate) enum Source {
    Plain(File),
    Bgzf {
        file: File,
        /// (compressed, uncompressed) offsets of blocks after the first.
        blocks: Vec<(u64, u64)>,
        /// The last block read (shared by threads reading one reference).
        last: Mutex<CachedBlock>,
    },
}

impl Source {
    /// Opens a file, telling plain, BGZF, and (unsupported) plain gzip files apart.
    pub(crate) fn open(path: &Path, gzi: Option<&Path>) -> anyhow::Result<Self> {
        let mut file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut magic = [0; 18];
        let n = read_full(&mut file, &mut magic)?;
        let file = File::open(path)?;

        if n < 2 || magic[..2] != [0x1f, 0x8b] {
            return Ok(Self::Plain(file));
        }

        if !is_bgzf_header(&magic[..n]) {
            bail!("Cannot index files compressed with gzip, please use bgzip");
        }

        let blocks = match gzi {
            Some(gzi) if gzi.exists() => bgzf::gzi::fs::read(gzi)
                .with_context(|| format!("Failed to load .gzi index: {}", gzi.display()))?
                .as_ref()
                .to_vec(),
            _ => Vec::new(),
        };

        Ok(Self::Bgzf {
            file,
            blocks,
            last: Mutex::default(),
        })
    }

    pub(crate) fn is_compressed(&self) -> bool {
        matches!(self, Self::Bgzf { .. })
    }

    /// Reads `len` uncompressed bytes from `offset`.
    fn read_at(&self, offset: u64, len: usize, dst: &mut Vec<u8>) -> anyhow::Result<()> {
        dst.clear();

        match self {
            Self::Plain(file) => {
                dst.resize(len, 0);
                read_exact_at(file, dst, offset).context("unexpected end of file")?;
            }
            Self::Bgzf { file, blocks, last } => {
                // The last block starting at or before the offset.
                let i = blocks.partition_point(|&(_, u)| u <= offset);
                let (mut coffset, mut uoffset) = if i == 0 { (0, 0) } else { blocks[i - 1] };
                let mut last = last.lock().expect("BGZF block cache lock poisoned");

                while dst.len() < len {
                    let Some((block, size)) = last.load(file, coffset)? else {
                        bail!("unexpected end of file");
                    };

                    let end = uoffset + block.len() as u64;

                    if end > offset {
                        let start = offset.saturating_sub(uoffset) as usize;
                        let take = (block.len() - start).min(len - dst.len());
                        dst.extend_from_slice(&block[start..start + take]);
                    }

                    coffset += size;
                    uoffset = end;
                }
            }
        }

        Ok(())
    }
}

/// Whether a gzip header has BGZF's `BC` extra subfield.
fn is_bgzf_header(h: &[u8]) -> bool {
    h.len() >= 18 && h[3] & 4 != 0 && h[12] == b'B' && h[13] == b'C'
}

fn read_full(file: &mut File, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;

    while n < buf.len() {
        match file.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }

    Ok(n)
}

/// Reads and inflates BGZF blocks, reusing one buffer for their compressed frames.
#[derive(Default)]
struct BlockReader {
    frame: Vec<u8>,
}

impl BlockReader {
    /// Reads and inflates the BGZF block at `offset` into `dst`; returns its compressed size, or
    /// `None` at the end of the file.
    fn read(&mut self, file: &File, offset: u64, dst: &mut Vec<u8>) -> anyhow::Result<Option<u64>> {
        let mut header = [0; BGZF_HEADER_LEN];

        // Only a file ending exactly at a block boundary ends cleanly; part of a header (as part
        // of a block) means the file is truncated, as in HTSlib's `bgzf_read_block`.
        match read_at_most(file, &mut header, offset)? {
            0 => return Ok(None),
            BGZF_HEADER_LEN => {}
            n => bail!(
                "Failed to read BGZF header at offset {offset}: read {n} of {BGZF_HEADER_LEN} bytes"
            ),
        }

        let size = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;

        if !is_bgzf_header(&header) || size < BGZF_HEADER_LEN + BGZF_FOOTER_LEN {
            bail!("invalid BGZF block at offset {offset}");
        }

        self.frame.resize(size, 0);
        let n = read_at_most(file, &mut self.frame, offset)?;

        if n < size {
            bail!("Failed to read BGZF block data at offset {offset}: read {n} of {size} bytes");
        }

        bgzf::io::raw::inflate_frame(&self.frame, dst)
            .with_context(|| format!("invalid BGZF block at offset {offset}"))?;

        Ok(Some(size as u64))
    }
}

/// The last BGZF block read, which the next request reuses when it starts in the same block, as
/// HTSlib keeps its current block.
#[derive(Default)]
pub(crate) struct CachedBlock {
    reader: BlockReader,
    /// The compressed offset and size of the block in `data`.
    position: Option<(u64, u64)>,
    data: Vec<u8>,
}

impl CachedBlock {
    /// The uncompressed data and compressed size of the block at `offset`, or `None` at the end
    /// of the file.
    fn load(&mut self, file: &File, offset: u64) -> anyhow::Result<Option<(&[u8], u64)>> {
        if let Some((start, size)) = self.position
            && start == offset
        {
            return Ok(Some((&self.data, size)));
        }

        self.position = None;

        let Some(size) = self.reader.read(file, offset, &mut self.data)? else {
            return Ok(None);
        };

        self.position = Some((offset, size));
        Ok(Some((&self.data, size)))
    }
}

fn read_at_most(file: &File, buf: &mut [u8], offset: u64) -> io::Result<usize> {
    let mut n = 0;

    while n < buf.len() {
        match read_at(file, &mut buf[n..], offset + n as u64)? {
            0 => break,
            k => n += k,
        }
    }

    Ok(n)
}

/// The uncompressed bytes of a source as a stream, recording BGZF block offsets.
struct Bytes<'a> {
    source: &'a Source,
    buf: Vec<u8>,
    pos: usize,
    /// The uncompressed offset of `buf[0]`.
    base: u64,
    next_coffset: u64,
    plain: Option<BufReader<&'a File>>,
    reader: BlockReader,
    blocks: Vec<(u64, u64)>,
}

impl<'a> Bytes<'a> {
    fn new(source: &'a Source) -> Self {
        let plain = match source {
            Source::Plain(file) => Some(BufReader::with_capacity(1 << 16, file)),
            Source::Bgzf { .. } => None,
        };

        Self {
            source,
            buf: Vec::new(),
            pos: 0,
            base: 0,
            next_coffset: 0,
            plain,
            reader: BlockReader::default(),
            blocks: Vec::new(),
        }
    }

    /// `bgzf_utell`
    fn tell(&self) -> u64 {
        self.base + self.pos as u64
    }

    fn next(&mut self) -> anyhow::Result<Option<u8>> {
        while self.pos >= self.buf.len() {
            self.base += self.buf.len() as u64;
            self.pos = 0;

            match (&mut self.plain, self.source) {
                (Some(reader), _) => {
                    self.buf.clear();
                    let chunk = reader.fill_buf()?;

                    if chunk.is_empty() {
                        return Ok(None);
                    }

                    self.buf.extend_from_slice(chunk);
                    let n = chunk.len();
                    reader.consume(n);
                }
                (None, Source::Bgzf { file, .. }) => {
                    let coffset = self.next_coffset;

                    let Some(size) = self.reader.read(file, coffset, &mut self.buf)? else {
                        // `buf` still holds the last block, which is empty for an EOF marker.
                        if !self.buf.is_empty() {
                            eprintln!(
                                "[W::bgzf_read_block] EOF marker is absent. The input may be truncated"
                            );
                        }

                        self.buf.clear();
                        return Ok(None);
                    };

                    // As HTSlib, the empty end-of-file block is not indexed.
                    if coffset > 0 && !self.buf.is_empty() {
                        self.blocks.push((coffset, self.base));
                    }

                    self.next_coffset += size;
                }
                (None, Source::Plain(_)) => unreachable!("plain sources have a reader"),
            }
        }

        let b = self.buf[self.pos];
        self.pos += 1;
        Ok(Some(b))
    }
}

/// HTSlib's `isgraph`.
fn is_graph(c: u8) -> bool {
    c > b' ' && c <= b'~'
}

/// A built index: the file's format, the index, and BGZF block offsets for a `.gzi`.
type Built = (SequenceFormat, Index, Vec<(u64, u64)>);

/// `fai_build_core`: indexes a FASTA or FASTQ file.
fn build(source: &Source) -> anyhow::Result<Built> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum State {
        OutRead,
        InName,
        InSeq,
        SeqEnd,
        InQual,
    }

    let mut bytes = Bytes::new(source);
    let mut index = Index::default();
    let mut format: Option<SequenceFormat> = None;
    let mut state = State::OutRead;
    let mut read_done = false;
    let mut line_num = 1;
    let mut name: Vec<u8> = Vec::new();
    let (mut seq_offset, mut qual_offset) = (0, 0);
    let (mut seq_len, mut qual_len) = (0u64, 0u64);
    let (mut char_len, mut line_len) = (0u64, 0u64);

    let entry = |name: &[u8], seq_len, line_len, char_len, seq_offset, qual_offset| Entry {
        name: name.to_vec(),
        len: seq_len,
        seq_offset,
        line_blen: char_len,
        line_len,
        qual_offset,
    };

    // As in HTSlib, every iteration reads one new character.
    while let Some(b) = bytes.next()? {
        match state {
            State::OutRead => match b {
                b'>' => {
                    if format == Some(SequenceFormat::Fastq) {
                        bail!("Found '>' in a FASTQ file, error at line {line_num}");
                    }

                    format = Some(SequenceFormat::Fasta);
                    state = State::InName;
                }
                b'@' => {
                    if format == Some(SequenceFormat::Fasta) {
                        bail!("Found '@' in a FASTA file, error at line {line_num}");
                    }

                    format = Some(SequenceFormat::Fastq);
                    state = State::InName;
                }
                b'\r' => {
                    if bytes.next()? == Some(b'\n') {
                        line_num += 1;
                    } else {
                        bail!(
                            "Format error, carriage return not followed by new line at line {line_num}"
                        );
                    }
                }
                b'\n' => line_num += 1,
                _ => bail!(
                    "Format error, unexpected \"{}\" at line {line_num}",
                    char::from(b)
                ),
            },
            State::InName => {
                if read_done {
                    index.insert(entry(
                        &name,
                        seq_len,
                        line_len,
                        char_len,
                        seq_offset,
                        qual_offset,
                    ))?;
                    read_done = false;
                }

                name.clear();
                let mut ch = Some(b);

                while let Some(x) = ch {
                    if !x.is_ascii_whitespace() {
                        name.push(x);
                    } else if !name.is_empty() || x == b'\n' {
                        break;
                    }

                    ch = bytes.next()?;
                }

                let Some(mut x) = ch else {
                    bail!(
                        "The last entry '{}' has no sequence at line {line_num}",
                        String::from_utf8_lossy(&name)
                    );
                };

                while x != b'\n' {
                    match bytes.next()? {
                        Some(y) => x = y,
                        None => break,
                    }
                }

                state = State::InSeq;
                seq_len = 0;
                qual_len = 0;
                char_len = 0;
                line_len = 0;
                seq_offset = bytes.tell();
                line_num += 1;
            }
            State::InSeq => {
                if format == Some(SequenceFormat::Fasta) {
                    if b == b'\n' {
                        state = State::OutRead;
                        line_num += 1;
                        continue;
                    } else if b == b'>' {
                        state = State::InName;
                        continue;
                    }
                } else if b == b'+' {
                    state = State::InQual;
                    skip_line(&mut bytes)?;
                    qual_offset = bytes.tell();
                    line_num += 1;
                    continue;
                } else if b == b'\n' {
                    bail!(
                        "Inlined empty line is not allowed in sequence '{}' at line {line_num}",
                        String::from_utf8_lossy(&name)
                    );
                }

                if format == Some(SequenceFormat::Fasta) {
                    read_done = true;
                }

                let (ll, cl) = scan_line(&mut bytes, b)?;
                seq_len += cl;

                if line_len == 0 {
                    line_len = ll;
                    char_len = cl;
                } else if line_len > ll {
                    state = if format == Some(SequenceFormat::Fasta) {
                        State::OutRead
                    } else {
                        State::SeqEnd
                    };
                } else if line_len < ll {
                    bail!(
                        "Different line length in sequence '{}' at line {line_num}",
                        String::from_utf8_lossy(&name)
                    );
                }

                line_num += 1;
            }
            State::SeqEnd => {
                if b == b'+' {
                    state = State::InQual;
                    skip_line(&mut bytes)?;
                    qual_offset = bytes.tell();
                    line_num += 1;
                } else {
                    bail!(
                        "Format error, expecting '+', got '{}' at line {line_num}",
                        char::from(b)
                    );
                }
            }
            State::InQual => {
                if b == b'\n' {
                    if !read_done {
                        bail!(
                            "Inlined empty line is not allowed in quality of sequence '{}' at line {line_num}",
                            String::from_utf8_lossy(&name)
                        );
                    }

                    state = State::OutRead;
                    line_num += 1;
                    continue;
                } else if b == b'@' && read_done {
                    state = State::InName;
                    continue;
                }

                let (ll, cl) = scan_line(&mut bytes, b)?;
                qual_len += cl;

                if line_len < ll {
                    bail!(
                        "Quality line length too long in '{}' at line {line_num}",
                        String::from_utf8_lossy(&name)
                    );
                } else if qual_len == seq_len {
                    read_done = true;
                } else if qual_len > seq_len {
                    bail!(
                        "Quality length longer than sequence in '{}' at line {line_num}",
                        String::from_utf8_lossy(&name)
                    );
                } else if line_len > ll {
                    bail!(
                        "Quality line length too short in '{}' at line {line_num}",
                        String::from_utf8_lossy(&name)
                    );
                }

                line_num += 1;
            }
        }
    }

    if !read_done {
        bail!("File truncated at line {line_num}");
    }

    index.insert(entry(
        &name,
        seq_len,
        line_len,
        char_len,
        seq_offset,
        qual_offset,
    ))?;

    Ok((
        format.unwrap_or(SequenceFormat::Fasta),
        index,
        std::mem::take(&mut bytes.blocks),
    ))
}

/// Skips to the end of the line.
fn skip_line(bytes: &mut Bytes<'_>) -> anyhow::Result<()> {
    while let Some(b) = bytes.next()? {
        if b == b'\n' {
            break;
        }
    }

    Ok(())
}

/// A line from its first byte: its length with the line ending, and its bases (visible
/// characters).
fn scan_line(bytes: &mut Bytes<'_>, first: u8) -> anyhow::Result<(u64, u64)> {
    let (mut ll, mut cl) = (1, u64::from(is_graph(first)));

    while let Some(b) = bytes.next()? {
        if b == b'\n' {
            break;
        }

        ll += 1;
        cl += u64::from(is_graph(b));
    }

    Ok((ll + 1, cl))
}

/// `fai_build3`: builds and writes `fai` (and `gzi` for BGZF input).
fn build_and_write(input: &Path, fai_path: &Path, gzi_path: &Path) -> anyhow::Result<()> {
    let source = Source::open(input, None)?;
    let (format, index, blocks) = build(&source)?;

    // Both indexes are complete before either replaces an existing one.
    let gzi = if source.is_compressed() {
        let index = blocks.into();
        let gzi = PendingFile::write(gzi_path, |out| {
            bgzf::gzi::io::Writer::new(out).write_index(&index)
        })
        .with_context(|| format!("Failed to make bgzf index {}", gzi_path.display()))?;
        Some(gzi)
    } else {
        None
    };

    let fai = PendingFile::write(fai_path, |out| index.write(format, out))
        .with_context(|| format!("Failed to write {format:?} index {}", fai_path.display()))?;

    if let Some(gzi) = gzi {
        gzi.commit()
            .with_context(|| format!("Failed to make bgzf index {}", gzi_path.display()))?;
    }

    fai.commit()
        .with_context(|| format!("Failed to write {format:?} index {}", fai_path.display()))
}

/// A file written in full to a temporary file next to its destination, which replaces the
/// destination on [`PendingFile::commit`], so a failed or abandoned write leaves an existing file
/// as it was. Destinations that are not regular files (`/dev/stdout`, FIFOs) are written in place.
enum PendingFile {
    Temporary {
        file: tempfile::NamedTempFile,
        dst: PathBuf,
    },
    InPlace,
}

impl PendingFile {
    fn write<F>(dst: &Path, write: F) -> io::Result<Self>
    where
        F: FnOnce(&mut BufWriter<&mut File>) -> io::Result<()>,
    {
        // Replace the target of a symbolic link, not the link.
        let dst = match fs::canonicalize(dst) {
            Ok(path) => path,
            Err(e) if e.kind() == io::ErrorKind::NotFound => dst.to_path_buf(),
            Err(e) => return Err(e),
        };

        if fs::metadata(&dst).is_ok_and(|m| !m.is_file()) {
            let mut file = File::create(&dst)?;
            let mut out = BufWriter::new(&mut file);
            write(&mut out)?;
            out.flush()?;
            return Ok(Self::InPlace);
        }

        let dir = match dst.parent() {
            Some(dir) if !dir.as_os_str().is_empty() => dir,
            _ => Path::new("."),
        };

        let name = dst.file_name().unwrap_or_default().to_string_lossy();
        let prefix = format!(".{name}.");
        let mut builder = tempfile::Builder::new();
        builder.prefix(&prefix).suffix(".tmp");

        // As a newly created file (`0666` less the umask), not `tempfile`'s `0600`.
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            builder.permissions(fs::Permissions::from_mode(0o666));
        }

        let mut file = builder.tempfile_in(dir)?;
        let mut out = BufWriter::new(file.as_file_mut());
        write(&mut out)?;
        out.flush()?;
        drop(out);

        Ok(Self::Temporary { file, dst })
    }

    fn commit(self) -> io::Result<()> {
        match self {
            Self::Temporary { file, dst } => file.persist(dst).map(drop).map_err(|e| e.error),
            Self::InPlace => Ok(()),
        }
    }
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// Builds a FASTA or FASTQ file's `.fai` (and `.gzi`) next to it.
pub(crate) fn build_index(src: &Path) -> anyhow::Result<()> {
    build_and_write(src, &with_suffix(src, ".fai"), &with_suffix(src, ".gzi"))
}

/// Strand marks for names (`--mark-strand`).
struct StrandMarks {
    forward: String,
    reverse: String,
}

impl StrandMarks {
    fn parse(s: &str) -> anyhow::Result<Self> {
        let (forward, reverse) = match s {
            "no" => ("", ""),
            "sign" => ("(+)", "(-)"),
            "rc" => ("", "/rc"),
            _ => match s.strip_prefix("custom,") {
                Some(custom) => custom.split_once(',').unwrap_or((custom, "")),
                None => bail!("[faidx] Unknown --mark-strand option \"{s}\""),
            },
        };

        Ok(Self {
            forward: forward.into(),
            reverse: reverse.into(),
        })
    }
}

/// The writer for extracted sequences, kept concrete so that finishing BGZF output (its last
/// block and EOF marker) reports errors rather than leaving them to `Drop`.
enum Sink {
    Plain(BufWriter<Box<dyn Write>>),
    Bgzf(bgzf::io::Writer<File>),
}

impl Sink {
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.write_all(buf),
            Self::Bgzf(writer) => writer.write_all(buf),
        }
    }

    /// Flushes plain output, or writes BGZF output's last block and EOF marker.
    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut writer) => writer.flush(),
            Self::Bgzf(writer) => writer.finish().map(drop),
        }
    }
}

/// Where extracted sequences go: plain or BGZF, tracking offsets for `--write-index`.
struct Output {
    writer: Sink,
    written: u64,
    index: Option<Index>,
}

impl Output {
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        self.written += buf.len() as u64;
        self.writer.write_all(buf)
    }
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    let format = if args.fastq {
        SequenceFormat::Fastq
    } else {
        SequenceFormat::Fasta
    };

    run_with_format(args, format)
}

pub fn run_fqidx(args: Args) -> anyhow::Result<ExitCode> {
    run_with_format(args, SequenceFormat::Fastq)
}

fn run_with_format(args: Args, format: SequenceFormat) -> anyhow::Result<ExitCode> {
    let strand = StrandMarks::parse(&args.mark_strand)?;

    let line_length = match args.length {
        Some(n) if n < 0 => {
            eprintln!("[faidx] bad line length '{n}', using default:{DEFAULT_LINE_LENGTH}");
            Some(DEFAULT_LINE_LENGTH)
        }
        n => n,
    };

    let gzi_path = args
        .gzi_idx
        .clone()
        .unwrap_or_else(|| with_suffix(&args.input, ".gzi"));

    // With no regions, build the index (to -o, if given).
    if args.regions.is_empty() && args.region_file.is_none() {
        let fai_path = args
            .fai_idx
            .clone()
            .or_else(|| args.output.clone())
            .unwrap_or_else(|| with_suffix(&args.input, ".fai"));

        build_and_write(&args.input, &fai_path, &gzi_path)
            .with_context(|| format!("[faidx] Could not build fai index {}", fai_path.display()))?;

        return Ok(ExitCode::SUCCESS);
    }

    let fai_path = args
        .fai_idx
        .clone()
        .unwrap_or_else(|| with_suffix(&args.input, ".fai"));

    let mut source = Source::open(&args.input, Some(&gzi_path))?;

    // `FAI_CREATE`: build a missing index (or a missing .gzi for BGZF input).
    if !fai_path.exists() || (source.is_compressed() && !gzi_path.exists()) {
        build_and_write(&args.input, &fai_path, &gzi_path)
            .with_context(|| format!("[faidx] Could not load fai index {}", fai_path.display()))?;
        source = Source::open(&args.input, Some(&gzi_path))?;
    }

    let index = Index::read(&fai_path, format)
        .with_context(|| format!("[faidx] Could not load fai index {}", fai_path.display()))?;

    let writer = match &args.output {
        Some(path) => {
            if path == &args.input {
                bail!("[faidx] Same input/output : {}", path.display());
            }

            let file = File::create(path).with_context(|| {
                format!("[faidx] Cannot open \"{}\" for writing", path.display())
            })?;

            let compressed = path
                .extension()
                .is_some_and(|ext| ext == "gz" || ext == "bgz" || ext == "bgzf");

            if compressed {
                if args.write_index {
                    bail!("--write-index with compressed output is not supported");
                }

                let level = args.output_format_options.level(DEFAULT_BGZF_LEVEL)?;
                let level = bgzf::io::writer::CompressionLevel::new(level)
                    .with_context(|| format!("invalid compression level: {level}"))?;

                Sink::Bgzf(
                    bgzf::io::writer::Builder::default()
                        .set_compression_level(level)
                        .build_from_writer(file),
                )
            } else {
                Sink::Plain(BufWriter::new(Box::new(file)))
            }
        }
        None => Sink::Plain(BufWriter::new(Box::new(io::stdout().lock()))),
    };

    let mut out = Output {
        writer,
        written: 0,
        index: (args.write_index && args.output.is_some()).then(Index::default),
    };

    let extract = Extract {
        index: &index,
        source: &source,
        format,
        line_length,
        reverse: args.reverse_complement,
        strand: &strand,
        continue_on_missing: args.continue_on_missing,
    };

    let mut regions = Vec::new();

    if let Some(path) = &args.region_file {
        let text = fs::read(path).with_context(|| {
            format!("[faidx] Failed to open \"{}\" for reading.", path.display())
        })?;

        for line in text.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            regions.push(String::from_utf8_lossy(line).into_owned());
        }

        // A final newline does not start another region.
        if text.ends_with(b"\n") {
            regions.pop();
        }
    }

    regions.extend(args.regions.iter().cloned());

    let mut status = ExitCode::SUCCESS;

    for region in &regions {
        if !extract.write(region, &mut out)? {
            status = ExitCode::FAILURE;
            break;
        }
    }

    let output_name = args
        .output
        .as_deref()
        .map_or_else(|| "-".into(), |path| path.display().to_string());

    out.writer
        .finish()
        .with_context(|| format!("[faidx] Failed to close output {output_name}"))?;

    if let (Some(index), Some(path)) = (&out.index, &args.output) {
        let fai_path = with_suffix(path, ".fai");
        PendingFile::write(&fai_path, |fai| index.write(format, fai))
            .and_then(PendingFile::commit)
            .with_context(|| format!("[faidx] Failed to write index {}", fai_path.display()))?;
    }

    Ok(status)
}

/// Region extraction settings.
struct Extract<'a> {
    index: &'a Index,
    source: &'a Source,
    format: SequenceFormat,
    line_length: Option<i64>,
    reverse: bool,
    strand: &'a StrandMarks,
    continue_on_missing: bool,
}

impl Extract<'_> {
    /// `write_output`: one region; returns `false` if it failed (and processing should stop).
    fn write(&self, region: &str, out: &mut Output) -> anyhow::Result<bool> {
        let parsed = hts_region::parse_region_unclamped(region, |name| {
            self.index.by_name.get(name).copied()
        });

        let entry = parsed.map(|(id, _, _)| &self.index.entries[id]);

        let wrap = match self.line_length {
            Some(n) => n,
            None => entry.map_or(-1, |e| e.line_blen as i64),
        };
        let wrap = if wrap <= 0 { usize::MAX } else { wrap as usize };

        let mark = if self.reverse {
            &self.strand.reverse
        } else {
            &self.strand.forward
        };

        let prefix = match self.format {
            SequenceFormat::Fasta => '>',
            SequenceFormat::Fastq => '@',
        };

        let header = format!("{prefix}{region}{mark}\n");
        let header_len = header.len() as u64;
        out.write_all(header.as_bytes())?;

        let Some((id, beg, end)) = parsed else {
            eprintln!(
                "[W::fai_get_val] Reference {region} not found in FASTA file, returning empty sequence"
            );
            eprintln!("[faidx] Failed to fetch sequence in {region}");
            return Ok(self.continue_on_missing);
        };

        let entry = &self.index.entries[id];
        let len = entry.len as i64;
        let (fbeg, fend) = (beg.min(len), end.min(len));
        let fbeg = fbeg.min(fend);

        // A start of 0 (`a:0-4`, a BED-style start) gives HTSlib a start of -1, which
        // `fai_retrieve` rejects: an empty sequence, with the warning below.
        let (fbeg, fend) = if beg < 0 { (0, 0) } else { (fbeg, fend) };

        let mut seq = Vec::new();
        self.fetch(entry, entry.seq_offset, fbeg, fend, &mut seq)?;

        if self.reverse {
            reverse_complement(&mut seq);
        }

        let seq_len = seq.len();

        // samtools warns for the sequence and again for the qualities of FASTQ.
        let warn = || {
            if seq_len == 0 {
                eprintln!("[faidx] Zero length sequence: {region}");
            } else if end < hts_region::POS_MAX && seq_len as i64 != end - beg {
                eprintln!("[faidx] Truncated sequence: {region}");
            }
        };

        warn();

        let wrapped = write_wrapped(out, &seq, wrap)?;
        let seq_offset = out.written - wrapped;
        let mut qual_offset = 0;

        if self.format == SequenceFormat::Fastq {
            let mut qual = Vec::new();
            self.fetch(entry, entry.qual_offset, fbeg, fend, &mut qual)?;
            warn();

            if self.reverse {
                qual.reverse();
            }

            out.write_all(b"+\n")?;
            qual_offset = out.written;
            write_wrapped(out, &qual, wrap)?;
        }

        if let Some(index) = &mut out.index {
            index.insert(Entry {
                name: format!("{region}{mark}").into_bytes(),
                len: seq_len as u64,
                seq_offset: seq_offset.max(header_len),
                line_blen: seq_len.min(wrap) as u64,
                line_len: seq_len.min(wrap) as u64 + 1,
                qual_offset,
            })?;
        }

        Ok(true)
    }

    /// `fai_retrieve`: the bases (or qualities) of `[beg, end)` from the block at `offset`.
    fn fetch(
        &self,
        entry: &Entry,
        offset: u64,
        beg: i64,
        end: i64,
        dst: &mut Vec<u8>,
    ) -> anyhow::Result<()> {
        dst.clear();

        if beg >= end {
            return Ok(());
        }

        if entry.line_blen == 0 {
            bail!("Invalid line length in index: 0");
        }

        let (beg, end) = (beg as u64, end as u64);
        let position = |p: u64| offset + p / entry.line_blen * entry.line_len + p % entry.line_blen;
        let first = position(beg);
        let last = position(end - 1);

        let mut raw = Vec::new();
        self.source
            .read_at(first, (last - first + 1) as usize, &mut raw)
            .context("Failed to retrieve block")?;

        // Drop the line endings between lines.
        let mut pos = beg;
        let mut i = 0;

        while pos < end {
            let in_line = (entry.line_blen - pos % entry.line_blen).min(end - pos) as usize;
            dst.extend_from_slice(&raw[i..i + in_line]);
            pos += in_line as u64;
            i += in_line + (entry.line_len - entry.line_blen) as usize;
        }

        Ok(())
    }
}

/// Writes lines of up to `width` characters; returns the bytes written.
fn write_wrapped(out: &mut Output, data: &[u8], width: usize) -> io::Result<u64> {
    let start = out.written;

    for line in data.chunks(width) {
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }

    Ok(out.written - start)
}

/// Reads the bases of the 1-based closed interval `[start, end]` into `sequence`, reading the
/// whole span at once and dropping line endings, as samtools does.
/// Reads the bases at the 1-based positions `start..=end` of a FASTA sequence (plain or BGZF),
/// skipping line endings.
pub(crate) fn read_bases(
    source: &Source,
    record: &fai::Record,
    start: usize,
    end: usize,
    raw: &mut Vec<u8>,
    sequence: &mut Vec<u8>,
) -> anyhow::Result<()> {
    sequence.clear();

    if start > end || start == 0 {
        return Ok(());
    }

    let line_bases = record.line_base_count().get();
    let line_width = record.line_width().get();

    // The file offset of the base at a 0-based position.
    let offset = |position: u64| {
        record.position() + position / line_bases * line_width + position % line_bases
    };

    let first = offset(u64::try_from(start - 1)?);
    let last = offset(u64::try_from(end - 1)?);
    let len = usize::try_from(last - first + 1)?;

    source.read_at(first, len, raw)?;

    // Copy whole lines between line endings rather than filtering byte by byte.
    for line in raw.split(|&b| b == b'\n') {
        sequence.extend_from_slice(line.strip_suffix(b"\r").unwrap_or(line));
    }

    Ok(())
}

#[cfg(unix)]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;
    file.read_exact_at(buf, offset)
}

#[cfg(unix)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<usize> {
    use std::os::unix::fs::FileExt;
    file.read_at(buf, offset)
}

#[cfg(not(unix))]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::io::{Seek, SeekFrom};

    let mut file = file;
    file.seek(SeekFrom::Start(offset))?;
    file.read_exact(buf)
}

#[cfg(not(unix))]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<usize> {
    use std::io::{Seek, SeekFrom};

    let mut file = file;
    file.seek(SeekFrom::Start(offset))?;
    file.read(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index_of(content: &[u8]) -> anyhow::Result<(SequenceFormat, Vec<u8>)> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x");
        fs::write(&path, content)?;
        let source = Source::open(&path, None)?;
        let (format, index, _) = build(&source)?;
        let mut out = Vec::new();
        index.write(format, &mut out)?;
        Ok((format, out))
    }

    #[test]
    fn test_build() -> anyhow::Result<()> {
        let (format, fai) = index_of(b">a desc\nACGT\nAC\n>b\nGG\n")?;
        assert_eq!(format, SequenceFormat::Fasta);
        assert_eq!(fai, b"a\t6\t8\t4\t5\nb\t2\t19\t2\t3\n");

        let (format, fai) = index_of(b"@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n")?;
        assert_eq!(format, SequenceFormat::Fastq);
        assert_eq!(fai, b"r1\t4\t4\t4\t5\t11\nr2\t2\t20\t2\t3\t25\n");

        assert!(index_of(b"@r1\nACGT\n+\nII\n").is_err());
        Ok(())
    }

    /// Writes `content` BGZF compressed to `path` and builds its `.gzi`; returns the source.
    fn bgzf_source(path: &Path, content: &[u8]) -> anyhow::Result<Source> {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(content)?;
        fs::write(path, writer.finish()?)?;

        let gzi = with_suffix(path, ".gzi");
        build_and_write(path, &with_suffix(path, ".fai"), &gzi)?;
        Source::open(path, Some(&gzi))
    }

    #[test]
    fn test_read_at_bgzf_spans_blocks() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x.fa.gz");
        let mut content = b">s\n".to_vec();
        content.extend((0..200_000u32).map(|i| b"ACGT"[(i % 7 % 4) as usize]));
        content.push(b'\n');
        let source = bgzf_source(&path, &content)?;

        let Source::Bgzf { blocks, .. } = &source else {
            panic!("expected a BGZF source");
        };
        assert!(blocks.len() >= 3, "{} blocks", blocks.len());

        // Nearby requests, requests within one block, and requests across block boundaries.
        let mut dst = Vec::new();

        for (offset, len) in [
            (0, 10),
            (5, 10),
            (65_000, 2_000),
            (65_500, 100),
            (199_994, 10),
        ] {
            source.read_at(offset, len, &mut dst)?;
            let offset = offset as usize;
            assert_eq!(dst, content[offset..offset + len], "{offset}+{len}");
        }

        assert!(source.read_at(199_995, 10, &mut dst).is_err());
        Ok(())
    }

    #[test]
    fn test_read_at_bgzf_reuses_block() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x.fa.gz");
        let source = bgzf_source(&path, b">s\nACGTACGTAC\n")?;

        let mut dst = Vec::new();
        source.read_at(3, 4, &mut dst)?;
        assert_eq!(dst, b"ACGT");

        // A request in the block just read does not read the file again, as HTSlib keeps its
        // current block.
        File::options().write(true).open(&path)?.set_len(0)?;
        source.read_at(9, 4, &mut dst)?;
        assert_eq!(dst, b"GTAC");
        Ok(())
    }

    #[test]
    fn test_read_at_bgzf_short_block_size() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x.fa.gz");
        let source = bgzf_source(&path, b">s\nACGT\n")?;

        // A block size too small for a BGZF block's header and footer.
        let mut bytes = fs::read(&path)?;
        bytes[16..18].copy_from_slice(&3u16.to_le_bytes());
        fs::write(&path, bytes)?;

        let mut dst = Vec::new();
        assert!(source.read_at(0, 4, &mut dst).is_err());
        Ok(())
    }

    #[test]
    fn test_read_at_bgzf_checksum_mismatch() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x.fa.gz");
        let source = bgzf_source(&path, b">s\nACGT\n")?;

        // Damage the first block's CRC32, as HTSlib's `bgzf_uncompress` checks it.
        let mut bytes = fs::read(&path)?;
        let size = usize::from(u16::from_le_bytes([bytes[16], bytes[17]])) + 1;
        bytes[size - BGZF_FOOTER_LEN] ^= 0xff;
        fs::write(&path, bytes)?;

        let mut dst = Vec::new();
        assert!(source.read_at(3, 4, &mut dst).is_err());
        Ok(())
    }

    /// The offsets of the BGZF blocks in `data`, and its length.
    fn block_offsets(data: &[u8]) -> Vec<usize> {
        let mut offsets = Vec::new();
        let mut offset = 0;

        while offset < data.len() {
            offsets.push(offset);
            offset += usize::from(u16::from_le_bytes([data[offset + 16], data[offset + 17]])) + 1;
        }

        offsets.push(offset);
        offsets
    }

    #[test]
    fn test_build_truncated_bgzf() -> anyhow::Result<()> {
        let mut content = b">s\n".to_vec();
        content.extend((0..150_000u32).map(|i| b"ACGT"[(i % 7 % 4) as usize]));
        content.push(b'\n');

        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(&content)?;
        let data = writer.finish()?;

        let offsets = block_offsets(&data);
        assert!(offsets.len() >= 4, "{offsets:?}");
        let (second, third) = (offsets[1], offsets[2]);

        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x.fa.gz");

        let build_cut = |len: usize| -> anyhow::Result<Built> {
            fs::write(&path, &data[..len])?;
            build(&Source::open(&path, None)?)
        };

        // Part of a header, of a block's data, or of its footer is an error, as in HTSlib.
        for len in [
            second + 1,
            second + 17,
            second + 18,
            second + 100,
            third - 4,
            third - 1,
        ] {
            assert!(
                build_cut(len).is_err(),
                "truncated at {len} of {}",
                data.len()
            );
        }

        // A file ending at a block boundary (here, without its EOF marker) ends cleanly.
        let (_, index, _) = build_cut(third)?;
        assert_eq!(index.entries.len(), 1);
        assert!(index.entries[0].len < 150_000);

        let (_, index, _) = build_cut(data.len())?;
        assert_eq!(index.entries[0].len, 150_000);

        Ok(())
    }

    #[test]
    fn test_build_and_write_keeps_index_on_failure() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("x.fa.gz");
        bgzf_source(
            &path,
            &[b">s\n".as_slice(), &[b'A'; 100_000], b"\n"].concat(),
        )?;

        let (fai, gzi) = (with_suffix(&path, ".fai"), with_suffix(&path, ".gzi"));
        let (fai_before, gzi_before) = (fs::read(&fai)?, fs::read(&gzi)?);

        // Truncate the input inside its second block's header.
        let data = fs::read(&path)?;
        let second = block_offsets(&data)[1];
        fs::write(&path, &data[..second + 10])?;

        assert!(build_and_write(&path, &fai, &gzi).is_err());
        assert_eq!(fs::read(&fai)?, fai_before);
        assert_eq!(fs::read(&gzi)?, gzi_before);

        // No temporary files are left behind.
        let mut names: Vec<_> = fs::read_dir(dir.path())?
            .map(|entry| entry.map(|e| e.file_name()))
            .collect::<io::Result<_>>()?;
        names.sort();
        assert_eq!(names, ["x.fa.gz", "x.fa.gz.fai", "x.fa.gz.gzi"]);

        Ok(())
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn test_sink_finish_reports_bgzf_errors() -> anyhow::Result<()> {
        // `/dev/full` fails every write, so the last block and EOF marker cannot be written.
        let file = File::options().write(true).open("/dev/full")?;
        let mut sink = Sink::Bgzf(bgzf::io::Writer::new(file));
        sink.write_all(b">s\nACGT\n")?;
        assert!(sink.finish().is_err());

        let file = File::options().write(true).open("/dev/full")?;
        let sink = Sink::Bgzf(bgzf::io::Writer::new(file));
        assert!(sink.finish().is_err(), "the EOF marker alone");

        Ok(())
    }

    #[test]
    fn test_strand_marks() -> anyhow::Result<()> {
        let marks = StrandMarks::parse("custom,+,-")?;
        assert_eq!((marks.forward.as_str(), marks.reverse.as_str()), ("+", "-"));
        assert!(StrandMarks::parse("bogus").is_err());
        Ok(())
    }
}
