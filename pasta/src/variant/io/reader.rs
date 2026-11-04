//! Opening VCF and BCF inputs with HTSlib's format detection, reading headers and records, and
//! iterating over indexed regions.

use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom},
    num::NonZero,
    path::{Path, PathBuf},
    sync::{Mutex, PoisonError},
};

use crate::bgzf::{self, VirtualPosition, io::writer::BGZF_EOF};
use crate::util::inputs::Opened;
use crate::variant::{
    Header, HeaderError, Record, bcf,
    index::{IndexFormat, QueryPlan, locate, split_idx_name, tabix::Tabix},
    log,
    vcf::{self, ParseError},
};

/// The length of the BGZF EOF marker block.
const EOF_LEN: usize = BGZF_EOF.len();

/// The data format of an input.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum InputFormat {
    /// VCF text.
    Vcf,
    /// BCF.
    Bcf,
}

/// The compression of an input.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Compression {
    /// Uncompressed.
    None,
    /// gzip that is not BGZF (not indexable).
    Gzip,
    /// BGZF.
    Bgzf,
}

/// Why an input could not be opened or read.
#[derive(Debug)]
#[non_exhaustive]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The input is neither VCF nor BCF (HTSlib's `file_type_error`).
    UnknownFormat,
    /// The input is binary data in no format HTSlib opens ([`is_unknown_binary`]): its
    /// `hts_open` fails with ENOEXEC, displayed as C's "Exec format error".
    UnknownBinary,
    /// The header could not be read (HTSlib's `header_error`).
    Header(HeaderError),
    /// A VCF line could not be parsed (HTSlib's `vcf_parse_error`).
    Parse(ParseError),
    /// A BCF record is malformed (HTSlib's `bcf_read_error`).
    Bcf(io::Error),
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::UnknownFormat => f.write_str("unknown file type"),
            Self::UnknownBinary => f.write_str("Exec format error"),
            Self::Header(_) => f.write_str("could not parse header"),
            Self::Parse(_) => f.write_str("VCF parse error"),
            Self::Bcf(_) => f.write_str("BCF read error"),
        }
    }
}

impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) | Self::Bcf(e) => Some(e),
            Self::Header(e) => Some(e),
            Self::Parse(e) => Some(e),
            Self::UnknownFormat | Self::UnknownBinary => None,
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

/// A file or standard input, replaying bytes read during format detection.
struct Input {
    /// `None` for standard input.
    file: Option<File>,
    prefix: Vec<u8>,
    prefix_pos: usize,
    /// How a BGZF input is checked for the EOF marker block when reading reaches its end, where
    /// HTSlib warns when it is missing.
    eof_check: Option<EofCheck>,
    /// The file recorded as an input while it is open.
    _opened: Option<Opened>,
}

/// How [`Input`] checks a BGZF input for the EOF marker block at its end.
enum EofCheck {
    /// The last bytes read from a stream, which cannot be checked up front.
    Stream(EofWatch),
    /// The end of a named regular file, which reading may reach after seeking.
    File(String),
}

impl EofCheck {
    /// Reports the input when it lacks the EOF marker block, now that reading has reached its
    /// end.
    fn finish(&self) {
        match self {
            Self::Stream(w) => w.finish(),
            Self::File(path) => {
                if let Ok(Some(false)) = has_bgzf_eof(Path::new(path)) {
                    report_missing_eof(path);
                }
            }
        }
    }
}

/// The tail of a BGZF stream, compared with the EOF marker block at the end.
struct EofWatch {
    /// The input's name, `-` for standard input.
    name: String,
    tail: [u8; EOF_LEN],
    len: usize,
}

impl EofWatch {
    fn new(name: String) -> Self {
        Self {
            name,
            tail: [0; EOF_LEN],
            len: 0,
        }
    }

    fn push(&mut self, data: &[u8]) {
        let n = data.len();
        if n >= EOF_LEN {
            self.tail.copy_from_slice(&data[n - EOF_LEN..]);
        } else {
            self.tail.copy_within(n.., 0);
            self.tail[EOF_LEN - n..].copy_from_slice(data);
        }
        self.len = (self.len + n).min(EOF_LEN);
    }

    /// Reports the input when the bytes read did not end with the EOF marker block.
    fn finish(&self) {
        if self.len < EOF_LEN || self.tail != BGZF_EOF {
            report_missing_eof(&self.name);
        }
    }
}

/// The first input reported by [`report_missing_eof`].
static MISSING_EOF: Mutex<Option<String>> = Mutex::new(None);

/// Prints HTSlib's warning for a BGZF input whose reading reached its end without the EOF
/// marker block (`bgzf_read_block`), and records `name` (`-` for standard input) for
/// [`take_missing_eof`]: every block was intact, but the input may have been cut at a block
/// boundary, as a writer killed between two blocks (or one that never writes the marker)
/// leaves it. HTSlib only warns; vtools then fails, after writing its output.
pub fn report_missing_eof(name: &str) {
    log::warn(
        "bgzf_read_block",
        format_args!("EOF marker is absent. The input may be truncated"),
    );
    let mut first = MISSING_EOF.lock().unwrap_or_else(PoisonError::into_inner);
    if first.is_none() {
        *first = Some(name.to_owned());
    }
}

/// Returns the first input reported by [`report_missing_eof`] since the last call, if any.
pub fn take_missing_eof() -> Option<String> {
    MISSING_EOF
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .take()
}

/// A reader of BGZF data that reports, when reading reaches its end, a stream that lacks the
/// EOF marker block ([`report_missing_eof`]), as HTSlib's `bgzf_read_block` warns. It wraps the
/// compressed bytes of an input read from its start to its end, without seeking.
pub struct EofWatched<R> {
    inner: R,
    watch: Option<EofWatch>,
}

impl<R: Read> EofWatched<R> {
    /// Watches `inner`, named `name` in the report (`-` for standard input).
    pub fn new(inner: R, name: &str) -> Self {
        Self {
            inner,
            watch: Some(EofWatch::new(name.to_owned())),
        }
    }
}

impl<R: Read> Read for EofWatched<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n > 0 {
            if let Some(w) = &mut self.watch {
                w.push(&buf[..n]);
            }
        } else if !buf.is_empty()
            && let Some(w) = self.watch.take()
        {
            w.finish();
        }
        Ok(n)
    }
}

impl Read for Input {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = if self.prefix_pos < self.prefix.len() {
            let n = (self.prefix.len() - self.prefix_pos).min(buf.len());
            buf[..n].copy_from_slice(&self.prefix[self.prefix_pos..self.prefix_pos + n]);
            self.prefix_pos += n;
            n
        } else {
            match &mut self.file {
                Some(f) => f.read(buf)?,
                None => io::stdin().read(buf)?,
            }
        };
        if n > 0 {
            if let Some(EofCheck::Stream(w)) = &mut self.eof_check {
                w.push(&buf[..n]);
            }
        } else if !buf.is_empty()
            && let Some(check) = self.eof_check.take()
        {
            // Once, as HTSlib's `no_eof_block` flag does.
            check.finish();
        }
        Ok(n)
    }
}

impl Seek for Input {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match &mut self.file {
            Some(f) => {
                self.prefix.clear();
                self.prefix_pos = 0;
                f.seek(pos)
            }
            None => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "cannot seek standard input",
            )),
        }
    }
}

/// A reader that keeps failing once it has failed, as HTSlib keeps a BGZF handle's error: a
/// gzip decoder that hit a corrupt deflate stream or a checksum mismatch reads on as if from a
/// new member, or reports the end.
struct StickyError<R> {
    inner: R,
    /// The kind and text of the first error.
    error: Option<(io::ErrorKind, String)>,
}

impl<R> StickyError<R> {
    fn new(inner: R) -> Self {
        Self { inner, error: None }
    }
}

impl<R: Read> Read for StickyError<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if let Some((kind, text)) = &self.error {
            return Err(io::Error::new(*kind, text.clone()));
        }
        self.inner.read(buf).inspect_err(|e| {
            if e.kind() != io::ErrorKind::Interrupted {
                self.error = Some((e.kind(), e.to_string()));
            }
        })
    }
}

enum Bgzf {
    Single(bgzf::io::Reader<Input>),
    Multi(bgzf::io::MultithreadedReader<Input>),
}

enum Stream {
    Plain(BufReader<Input>),
    Gzip(Box<BufReader<StickyError<flate2::read::MultiGzDecoder<Input>>>>),
    Bgzf(Box<Bgzf>),
}

impl Read for Stream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Plain(r) => r.read(buf),
            Self::Gzip(r) => r.read(buf),
            Self::Bgzf(b) => match &mut **b {
                Bgzf::Single(r) => r.read(buf),
                Bgzf::Multi(r) => r.read(buf),
            },
        }
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        match self {
            Self::Bgzf(b) => match &mut **b {
                Bgzf::Single(r) => r.read_exact(buf),
                Bgzf::Multi(r) => r.read_exact(buf),
            },
            Self::Plain(r) => r.read_exact(buf),
            Self::Gzip(r) => r.read_exact(buf),
        }
    }
}

impl BufRead for Stream {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Plain(r) => r.fill_buf(),
            Self::Gzip(r) => r.fill_buf(),
            Self::Bgzf(b) => match &mut **b {
                Bgzf::Single(r) => r.fill_buf(),
                Bgzf::Multi(r) => r.fill_buf(),
            },
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::Plain(r) => r.consume(amt),
            Self::Gzip(r) => r.consume(amt),
            Self::Bgzf(b) => match &mut **b {
                Bgzf::Single(r) => r.consume(amt),
                Bgzf::Multi(r) => r.consume(amt),
            },
        }
    }
}

/// A VCF or BCF input (HTSlib's `htsFile` opened for reading).
pub struct Reader {
    stream: Stream,
    format: InputFormat,
    compression: Compression,
    path: Option<String>,
}

impl std::fmt::Debug for Reader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reader")
            .field("format", &self.format)
            .field("compression", &self.compression)
            .field("path", &self.path)
            .finish_non_exhaustive()
    }
}

/// What [`Reader::open_impl`] accepts.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Content {
    /// VCF or BCF only.
    Variants,
    /// Any text, read as lines.
    AnyText,
}

/// The number of decompressed bytes inspected to detect the format.
const PEEK: usize = 1024;

/// The number of raw bytes read to find them.
const RAW_PEEK: usize = 1 << 16;

/// How an input is opened.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Opening {
    /// For the first time: a FIFO waits for its writer.
    First,
    /// Once more ([`open_again`]).
    Again,
}

/// Returns whether HTSlib's `hts_hopen` refuses an input whose first (decompressed) bytes are
/// `peek`, failing with ENOEXEC ("Exec format error"): they are neither empty, nor text, nor
/// the start of a binary format it opens (BAM, BCF, CRAM). Other inputs that are not VCF or BCF
/// open and fail when their header is read.
#[must_use]
pub fn is_unknown_binary(peek: &[u8]) -> bool {
    const MAGICS: [&[u8]; 4] = [b"CRAM", b"BAM\x01", b"BCF\x04", b"BCF\x02"];
    // HTSlib's `is_text_only`, over the bytes `hts_detect_format2` looks at.
    let peek = &peek[..peek.len().min(PEEK)];
    let is_text = peek
        .iter()
        .all(|&c| c >= b' ' || matches!(c, b'\t' | b'\r' | b'\n'));
    !is_text && !MAGICS.iter().any(|m| peek.starts_with(m))
}

/// Returns C's ENOEXEC error, "Exec format error", with which HTSlib's `hts_open` refuses
/// unknown binary data.
#[must_use]
pub fn exec_format_error() -> io::Error {
    #[cfg(unix)]
    {
        io::Error::from_raw_os_error(libc::ENOEXEC)
    }
    #[cfg(not(unix))]
    {
        io::Error::other("Exec format error")
    }
}

/// Opens a file that was read once already, as bcftools opens some inputs twice. A regular
/// file opens as usual. A stream ([`is_stream`]) is opened without waiting for a writer: the
/// rest of a pipe the first reading left is read, as HTSlib reads it, and a FIFO whose writer
/// is gone reads as empty, where HTSlib's `open` waits forever.
///
/// # Errors
///
/// Returns an I/O error when the file cannot be opened.
pub fn open_again(path: &Path) -> io::Result<File> {
    #[cfg(unix)]
    {
        use std::os::{fd::AsRawFd, unix::fs::OpenOptionsExt};

        if !std::fs::metadata(path)?.is_file() {
            let file = std::fs::OpenOptions::new()
                .read(true)
                .custom_flags(libc::O_NONBLOCK)
                .open(path)?;
            // Reads wait for data again, as long as a writer is left.
            let fd = file.as_raw_fd();
            // SAFETY: `fd` is the open descriptor of `file`, which outlives both calls.
            let flags = unsafe { libc::fcntl(fd, libc::F_GETFL) };
            // SAFETY: as above, with flags the first call returned.
            if flags < 0 || unsafe { libc::fcntl(fd, libc::F_SETFL, flags & !libc::O_NONBLOCK) } < 0
            {
                return Err(io::Error::last_os_error());
            }
            return Ok(file);
        }
    }
    File::open(path)
}

impl Reader {
    /// Opens a file (or standard input for `-`) and detects its format as HTSlib's
    /// `hts_detect_format` does. `threads` adds BGZF decompression threads. A `##idx##` suffix
    /// (an explicit index name) is ignored here.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be opened or is neither VCF nor BCF.
    pub fn open(path: &str, threads: usize) -> Result<Self, ReadError> {
        Self::open_impl(path, threads, Content::Variants, Opening::First)
    }

    /// [`Self::open`] for an input that was read once already ([`open_again`]): standard input
    /// goes on from where it is, and a stream is opened without waiting for a writer.
    ///
    /// # Errors
    ///
    /// As [`Self::open`].
    pub fn open_again(path: &str, threads: usize) -> Result<Self, ReadError> {
        Self::open_impl(path, threads, Content::Variants, Opening::Again)
    }

    /// Opens a text file (BGZF, gzip, or plain; `-` for standard input) whatever lines it
    /// holds, for [`Self::read_line`] and tabix queries ([`Self::query_next_line`]), as HTSlib
    /// opens a tabix-indexed annotation file. Its format is reported as VCF.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be opened.
    pub fn open_text(path: &str, threads: usize) -> Result<Self, ReadError> {
        Self::open_impl(path, threads, Content::AnyText, Opening::First)
    }

    fn open_impl(
        path: &str,
        threads: usize,
        content: Content,
        opening: Opening,
    ) -> Result<Self, ReadError> {
        let (data_path, _) = split_idx_name(path);
        let file = match opening {
            _ if data_path == "-" => None,
            Opening::First => Some(File::open(data_path)?),
            Opening::Again => Some(open_again(Path::new(data_path))?),
        };
        let opened = file.as_ref().and_then(Opened::new);
        let label = file.as_ref().map(|_| data_path.to_owned());
        let mut input = Input {
            file,
            prefix: Vec::new(),
            prefix_pos: 0,
            eof_check: None,
            _opened: opened,
        };
        // Read the first bytes, then replay them (files rewind instead).
        let mut head = vec![0u8; RAW_PEEK];
        let n = read_up_to(&mut input, &mut head)?;
        head.truncate(n);
        let is_gz = head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b;
        let compression = if is_gz {
            if head.len() >= 18 && head[3] & 4 != 0 && &head[12..16] == b"BC\x02\x00" {
                Compression::Bgzf
            } else {
                Compression::Gzip
            }
        } else {
            Compression::None
        };
        let peek = if is_gz {
            let mut d = flate2::read::MultiGzDecoder::new(&head[..]);
            let mut out = vec![0u8; PEEK];
            let mut got = 0;
            while got < PEEK {
                match d.read(&mut out[got..]) {
                    Ok(0) | Err(_) => break,
                    Ok(k) => got += k,
                }
            }
            out.truncate(got);
            out
        } else {
            head[..head.len().min(PEEK)].to_vec()
        };
        let format = if peek.len() >= 4 && peek[3] <= 4 && peek.starts_with(b"BCF\x02") {
            InputFormat::Bcf
        } else if peek.starts_with(b"##fileformat=VCF") || content == Content::AnyText {
            InputFormat::Vcf
        } else if is_unknown_binary(&peek) {
            return Err(ReadError::UnknownBinary);
        } else {
            return Err(ReadError::UnknownFormat);
        };
        // A named regular file rewinds. Standard input, and a pipe, FIFO, or device opened by
        // name (`<(...)`, `/dev/stdin`), cannot: their first bytes are replayed.
        let rewinds = match &input.file {
            Some(f) => f.metadata()?.is_file(),
            None => false,
        };
        if rewinds {
            input.seek(SeekFrom::Start(0))?;
            // VCF and BCF inputs are checked by their readers (the synced reader up front, as
            // HTSlib's `bgzf_check_EOF` does). A text file read as it goes (a tabix-indexed
            // targets or annotation file) is checked when reading reaches its end.
            if compression == Compression::Bgzf && content == Content::AnyText {
                input.eof_check = Some(EofCheck::File(data_path.to_owned()));
            }
        } else {
            input.prefix = head;
            // A regular file on standard input is checked up front, as HTSlib seeks it.
            let checked = input.file.is_none() && stdin_has_bgzf_eof().ok().flatten().is_some();
            if compression == Compression::Bgzf && !checked {
                let name = label.clone().unwrap_or_else(|| "-".to_owned());
                input.eof_check = Some(EofCheck::Stream(EofWatch::new(name)));
            }
        }
        let stream = match compression {
            Compression::None => Stream::Plain(BufReader::with_capacity(1 << 17, input)),
            Compression::Gzip => Stream::Gzip(Box::new(BufReader::with_capacity(
                1 << 17,
                StickyError::new(flate2::read::MultiGzDecoder::new(input)),
            ))),
            Compression::Bgzf => Stream::Bgzf(Box::new(match NonZero::new(threads) {
                None => Bgzf::Single(bgzf::io::Reader::new(input)),
                // Without a thread to start, the input is read on this thread.
                Some(n) => match bgzf::io::MultithreadedReader::try_with_worker_count(n, input) {
                    Ok(reader) => Bgzf::Multi(reader),
                    Err(e) => Bgzf::Single(bgzf::io::Reader::new(e.into_inner())),
                },
            })),
        };
        Ok(Self {
            stream,
            format,
            compression,
            path: label,
        })
    }

    /// Returns the data format.
    #[must_use]
    pub fn format(&self) -> InputFormat {
        self.format
    }

    /// Returns the compression.
    #[must_use]
    pub fn compression(&self) -> Compression {
        self.compression
    }

    /// Returns the file path (`None` for standard input).
    #[must_use]
    pub fn path(&self) -> Option<&str> {
        self.path.as_deref()
    }

    /// Reads the header (HTSlib's `bcf_hdr_read`). For VCF, contigs named in a tabix or CSI
    /// index next to the file but missing from the header are added, as HTSlib does.
    ///
    /// # Errors
    ///
    /// Returns an error when the header is missing or malformed.
    pub fn read_header(&mut self) -> Result<Header, ReadError> {
        match self.format {
            InputFormat::Bcf => bcf::read_header(&mut self.stream).map_err(|e| {
                if e.kind() == io::ErrorKind::InvalidData {
                    ReadError::Header(HeaderError::MissingSampleLine)
                } else {
                    ReadError::Io(e)
                }
            }),
            InputFormat::Vcf => self.read_vcf_header(),
        }
    }

    fn read_vcf_header(&mut self) -> Result<Header, ReadError> {
        let mut txt = Vec::new();
        let mut line = Vec::new();
        while self.read_line(&mut line)? {
            if line.is_empty() {
                continue;
            }
            if line[0] != b'#' {
                log::error("vcf_hdr_read", format_args!("No sample line"));
                return Err(ReadError::Header(HeaderError::MissingSampleLine));
            }
            txt.extend_from_slice(&line);
            txt.push(b'\n');
            if line.get(1) != Some(&b'#') {
                break;
            }
        }
        if txt.is_empty() {
            log::error("vcf_hdr_read", format_args!("Could not read the header"));
            return Err(ReadError::Header(HeaderError::MissingSampleLine));
        }
        let mut header = Header::parse(&txt).map_err(ReadError::Header)?;
        if let Some(path) = &self.path
            && let Some(idx) = locate(path, IndexFormat::Tbi)
            && let Ok(index) = crate::variant::index::Index::read(&idx)
        {
            match Tabix::from_index(index) {
                Some(tbx) => {
                    for name in tbx.names() {
                        if header.contig_id(name).is_none() {
                            let mut rec = crate::variant::HeaderRecord::structured("contig");
                            rec.push_field("ID", name, crate::variant::Quoting::Bare);
                            header.add_record(rec).map_err(ReadError::Header)?;
                        }
                    }
                }
                None => log::error(
                    "index_load",
                    format_args!("Invalid index header for {path}"),
                ),
            }
        }
        Ok(header)
    }

    /// Reads one line without its line ending (HTSlib's `hts_getline`), returning `false` at
    /// the end of the input.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn read_line(&mut self, out: &mut Vec<u8>) -> io::Result<bool> {
        out.clear();
        let mut found_newline = false;
        loop {
            let buf = self.stream.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            match memchr::memchr(b'\n', buf) {
                Some(i) => {
                    out.extend_from_slice(&buf[..i]);
                    self.stream.consume(i + 1);
                    found_newline = true;
                    break;
                }
                None => {
                    let n = buf.len();
                    out.extend_from_slice(buf);
                    self.stream.consume(n);
                }
            }
        }
        if !found_newline && out.is_empty() {
            return Ok(false);
        }
        let strip_cr = match self.compression {
            Compression::None => found_newline,
            _ => true,
        };
        if strip_cr && out.last() == Some(&b'\r') {
            out.pop();
        }
        Ok(true)
    }

    /// Reads and parses the next record (HTSlib's `bcf_read`), returning `false` at the end.
    /// VCF lines with undefined tags extend the header, as HTSlib does.
    ///
    /// # Errors
    ///
    /// Returns an error for an unparsable line or a malformed BCF record.
    pub fn read_record(
        &mut self,
        header: &mut Header,
        line: &mut Vec<u8>,
        rec: &mut Record,
    ) -> Result<bool, ReadError> {
        match self.format {
            InputFormat::Vcf => {
                if !self.read_line(line)? {
                    return Ok(false);
                }
                if let Err(e) = vcf::parse(header, line, rec) {
                    self.check_gzip_rest()?;
                    return Err(ReadError::Parse(e));
                }
                Ok(true)
            }
            InputFormat::Bcf => match bcf::read_record(&mut self.stream, header, rec) {
                Ok(bcf::ReadOutcome::Record) => Ok(true),
                Ok(bcf::ReadOutcome::Eof) => Ok(false),
                Err(e) => Err(ReadError::Bcf(e)),
            },
        }
    }

    /// Reads the rest of a plain gzip (not BGZF) input, discarding it, to check what is left of
    /// its deflate streams and each member's CRC32 and length; other inputs are not read. A
    /// corrupt deflate stream may inflate to garbage before the error shows, and a VCF line of
    /// that garbage that does not parse ends the reading: call this then, so that the damaged
    /// input is a read failure rather than an unparsable line. (A BGZF block is checked before
    /// its data is used.)
    ///
    /// # Errors
    ///
    /// Returns the error of decompressing the rest.
    pub fn check_gzip_rest(&mut self) -> io::Result<()> {
        if self.compression == Compression::Gzip {
            io::copy(&mut self.stream, &mut io::sink())?;
        }
        Ok(())
    }

    /// Returns buffered input bytes, reading more when the buffer is empty (empty at the end),
    /// for callers that split lines themselves.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.stream.fill_buf()
    }

    /// Marks `amt` bytes returned by [`Self::fill_buf`] as read.
    pub fn consume(&mut self, amt: usize) {
        self.stream.consume(amt);
    }

    /// Reads one BCF record without checking it (HTSlib's `bcf_read1_core`).
    ///
    /// # Errors
    ///
    /// Returns an error for a truncated record.
    pub fn read_bcf_raw(&mut self, rec: &mut Record) -> io::Result<bool> {
        Ok(bcf::read_record_raw(&mut self.stream, rec)? == bcf::ReadOutcome::Record)
    }

    /// Reads up to `max` bytes of whole BCF records into `out` (each 32-byte prefix and blocks),
    /// for decoding elsewhere; returns the number of records.
    ///
    /// # Errors
    ///
    /// Returns an error for a truncated record or a failed read; `out` then holds the whole
    /// records read before it.
    pub fn read_bcf_batch(
        &mut self,
        out: &mut Vec<u8>,
        max_records: usize,
        max_bytes: usize,
    ) -> io::Result<usize> {
        let mut n = 0;
        while n < max_records && out.len() < max_bytes {
            let mut head = [0u8; 8];
            let got = read_up_to(&mut self.stream, &mut head)?;
            if got == 0 {
                break;
            }
            if got != 8 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "truncated BCF record",
                ));
            }
            let l_shared = u32::from_le_bytes(head[0..4].try_into().expect("4 bytes")) as usize;
            let l_indiv = u32::from_le_bytes(head[4..8].try_into().expect("4 bytes")) as usize;
            let start = out.len();
            out.extend_from_slice(&head);
            out.resize(start + 8 + l_shared + l_indiv, 0);
            match read_up_to(&mut self.stream, &mut out[start + 8..]) {
                Ok(got) if got == l_shared + l_indiv => {}
                Ok(_) => {
                    out.truncate(start);
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        "truncated BCF record",
                    ));
                }
                Err(e) => {
                    out.truncate(start);
                    return Err(e);
                }
            }
            n += 1;
        }
        Ok(n)
    }

    /// Returns the BGZF virtual offset of the next byte (HTSlib's `bgzf_tell`), or 0 for
    /// other compressions.
    #[must_use]
    pub fn virtual_position(&self) -> u64 {
        match &self.stream {
            Stream::Bgzf(b) => u64::from(match &**b {
                Bgzf::Single(r) => r.virtual_position(),
                Bgzf::Multi(r) => r.virtual_position(),
            }),
            _ => 0,
        }
    }

    /// Seeks to a BGZF virtual offset (HTSlib's `bgzf_seek`).
    ///
    /// # Errors
    ///
    /// Returns an error for non-BGZF input or standard input.
    pub fn seek(&mut self, voffset: u64) -> io::Result<()> {
        use crate::bgzf::io::Seek as _;
        let vp = VirtualPosition::from(voffset);
        match &mut self.stream {
            Stream::Bgzf(b) => match &mut **b {
                Bgzf::Single(r) => r.seek(vp).map(|_| ()),
                Bgzf::Multi(r) => r.seek_to_virtual_position(vp).map(|_| ()),
            },
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "only BGZF input can seek",
            )),
        }
    }

    /// Returns the next line of a tabix query (HTSlib's `hts_itr_next` with `tbx_readrec`),
    /// with its target and interval.
    ///
    /// # Errors
    ///
    /// Returns an error for an unparsable line or one on a sequence not in the index.
    pub fn query_next_line(
        &mut self,
        q: &mut Query,
        tbx: &mut Tabix,
        line: &mut Vec<u8>,
    ) -> io::Result<bool> {
        let meta_char = tbx.conf.meta_char;
        self.query_next(q, |r| {
            loop {
                if !r.read_line(line)? {
                    return Ok(None);
                }
                if line.is_empty() || i32::from(line[0]) != meta_char {
                    break;
                }
            }
            let (tid, intv) = tbx.line_interval(line)?;
            Ok(Some((i64::from(tid), intv.beg, intv.end)))
        })
    }

    /// Returns the next BCF record of an index query (HTSlib's `hts_itr_next` with
    /// `bcf_readrec`), checked against `header`.
    ///
    /// # Errors
    ///
    /// Returns an error for a malformed record.
    pub fn query_next_bcf(
        &mut self,
        q: &mut Query,
        header: &Header,
        rec: &mut Record,
    ) -> io::Result<bool> {
        self.query_next(q, |r| {
            if bcf::read_record_raw(&mut r.stream, rec)? == bcf::ReadOutcome::Eof {
                return Ok(None);
            }
            if bcf::check_record(Some(header), rec).is_err() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "BCF read error"));
            }
            Ok(Some((
                i64::from(rec.rid()),
                rec.pos(),
                rec.pos() + rec.rlen(),
            )))
        })
    }

    fn query_next<F>(&mut self, q: &mut Query, mut read: F) -> io::Result<bool>
    where
        F: FnMut(&mut Self) -> io::Result<Option<(i64, i64, i64)>>,
    {
        if q.finished {
            return Ok(false);
        }
        let chunks = &q.plan.chunks;
        loop {
            if q.curr_off == 0 || (q.i >= 0 && q.curr_off >= chunks[q.i as usize].1) {
                if q.i == chunks.len() as isize - 1 {
                    break;
                }
                let next = chunks[(q.i + 1) as usize];
                if q.i < 0 || chunks[q.i as usize].1 != next.0 {
                    self.seek(next.0)?;
                    q.curr_off = self.virtual_position();
                }
                q.i += 1;
            }
            match read(self) {
                Ok(Some((tid, beg, end))) => {
                    q.curr_off = self.virtual_position();
                    if q.plan.is_past(tid, beg) {
                        break;
                    }
                    if q.plan.overlaps(beg, end) {
                        return Ok(true);
                    }
                }
                Ok(None) => break,
                Err(e) => {
                    q.finished = true;
                    return Err(e);
                }
            }
        }
        q.finished = true;
        Ok(false)
    }
}

/// The progress of an index query (HTSlib's `hts_itr_t` state).
#[derive(Clone, Debug, Default)]
pub struct Query {
    plan: QueryPlan,
    i: isize,
    curr_off: u64,
    finished: bool,
}

impl Query {
    /// Starts a query from a plan.
    #[must_use]
    pub fn new(plan: QueryPlan) -> Self {
        let finished = plan.chunks.is_empty();
        Self {
            plan,
            i: -1,
            curr_off: 0,
            finished,
        }
    }

    /// Returns whether the query is done.
    #[must_use]
    pub fn is_finished(&self) -> bool {
        self.finished
    }
}

fn read_up_to<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// Checks whether a BGZF file ends with the EOF marker block (HTSlib's `bgzf_check_EOF`).
/// Returns `None` when the check cannot be made: for a pipe, FIFO, or device, which is not
/// opened again (a FIFO would block with no writer left).
///
/// # Errors
///
/// Returns an I/O error.
pub fn has_bgzf_eof(path: &Path) -> io::Result<Option<bool>> {
    if !std::fs::metadata(path)?.is_file() {
        return Ok(None);
    }
    let mut f = File::open(path)?;
    let len = f.seek(SeekFrom::End(0))?;
    if len < EOF_LEN as u64 {
        return Ok(Some(false));
    }
    f.seek(SeekFrom::End(-(EOF_LEN as i64)))?;
    let mut buf = [0u8; EOF_LEN];
    f.read_exact(&mut buf)?;
    Ok(Some(buf == BGZF_EOF))
}

/// [`has_bgzf_eof`] for standard input: HTSlib seeks standard input too, so the check is made
/// when it is a regular file (`< file`), and `None` is returned for a pipe or a terminal. The
/// read does not move the shared file offset.
///
/// # Errors
///
/// Returns an I/O error.
pub fn stdin_has_bgzf_eof() -> io::Result<Option<bool>> {
    #[cfg(unix)]
    {
        use std::os::{fd::AsFd, unix::fs::FileExt};

        let f = File::from(io::stdin().as_fd().try_clone_to_owned()?);
        let meta = f.metadata()?;
        if !meta.is_file() {
            return Ok(None);
        }
        let len = meta.len();
        if len < EOF_LEN as u64 {
            return Ok(Some(false));
        }
        let mut buf = [0u8; EOF_LEN];
        f.read_exact_at(&mut buf, len - EOF_LEN as u64)?;
        Ok(Some(buf == BGZF_EOF))
    }
    #[cfg(not(unix))]
    {
        Ok(None)
    }
}

/// Returns whether `name` is read as a stream, which cannot be opened again to read it from
/// the start: standard input (`-`), a pipe (`<(...)`, `/dev/stdin`), a FIFO, or a device.
/// What was read is gone, and opening a FIFO again with no writer left waits forever.
#[must_use]
pub fn is_stream(name: &str) -> bool {
    name == "-" || std::fs::metadata(data_path(name)).is_ok_and(|m| !m.is_file() && !m.is_dir())
}

/// Returns the data path of a name that may carry an `##idx##` index name.
#[must_use]
pub fn data_path(name: &str) -> PathBuf {
    PathBuf::from(split_idx_name(name).0)
}

#[cfg(all(test, unix))]
mod tests {
    use std::{io::Write, os::fd::AsRawFd, thread};

    use super::*;
    use crate::variant::io::{OutputType, WriteOptions, Writer};

    #[test]
    fn eof_watched_reports_a_stream_without_the_eof_marker() {
        let mut w = crate::bgzf::io::Writer::new(Vec::new());
        w.write_all(b"chr1\t1\t2\n").expect("write");
        let data = w.finish().expect("finish");
        // Only this test reports missing markers in this crate's tests.
        let _ = take_missing_eof();
        let read = |data: &[u8], name: &str| {
            let mut text = Vec::new();
            flate2::read::MultiGzDecoder::new(EofWatched::new(data, name))
                .read_to_end(&mut text)
                .expect("inflate");
            text
        };
        assert_eq!(read(&data, "intact"), b"chr1\t1\t2\n");
        assert_eq!(take_missing_eof(), None);
        let cut = &data[..data.len() - EOF_LEN];
        assert_eq!(read(cut, "cut"), b"chr1\t1\t2\n");
        assert_eq!(take_missing_eof().as_deref(), Some("cut"));
        assert_eq!(take_missing_eof(), None);
    }

    /// A VCF longer than the bytes read to detect the format, so reading goes on past them.
    fn long_vcf() -> Vec<u8> {
        let mut vcf = b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            .to_vec();
        for pos in 1..=5000 {
            vcf.extend_from_slice(format!("chr1\t{pos}\t.\tA\tC\t.\t.\t.\n").as_bytes());
        }
        assert!(vcf.len() > RAW_PEEK);
        vcf
    }

    /// Converts a VCF to `output_type` through files in the temporary directory.
    fn convert(vcf: &[u8], output_type: OutputType) -> Vec<u8> {
        // A directory per call, as tests run at the same time.
        static CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let call = CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!("pasta-pipe-{}-{call}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let (src, dst) = (dir.join("in.vcf"), dir.join("out"));
        std::fs::write(&src, vcf).expect("write the test VCF");
        let mut reader = Reader::open(&src.to_string_lossy(), 0).expect("open the test VCF");
        let mut header = reader.read_header().expect("read the header");
        let options = WriteOptions {
            output_type,
            level: None,
            threads: 0,
        };
        let mut writer = Writer::create(Some(&dst), &options).expect("create the output");
        writer.write_header(&header).expect("write the header");
        let (mut line, mut rec) = (Vec::new(), Record::new());
        while reader
            .read_record(&mut header, &mut line, &mut rec)
            .expect("read a record")
        {
            writer
                .write_record(&header, &mut rec)
                .expect("write a record");
        }
        writer.finish().expect("finish the output");
        let out = std::fs::read(&dst).expect("read the output");
        let _ = std::fs::remove_dir_all(&dir);
        out
    }

    /// Passes `data` through a pipe named `/dev/fd/N`, as `<(...)` names one, and returns the
    /// positions read from it.
    fn positions_through_pipe(data: Vec<u8>) -> Result<Vec<i64>, ReadError> {
        let (rx, mut tx) = io::pipe()?;
        // The writer stops with a broken pipe if the reader fails early.
        let feeder = thread::spawn(move || tx.write_all(&data));
        let path = format!("/dev/fd/{}", rx.as_raw_fd());
        let read = || -> Result<Vec<i64>, ReadError> {
            let mut reader = Reader::open(&path, 0)?;
            let mut header = reader.read_header()?;
            let (mut line, mut rec) = (Vec::new(), Record::new());
            let mut positions = Vec::new();
            while reader.read_record(&mut header, &mut line, &mut rec)? {
                positions.push(rec.pos() + 1);
            }
            Ok(positions)
        };
        let positions = read();
        drop(rx);
        let _ = feeder.join().expect("the writer thread does not panic");
        positions
    }

    #[test]
    fn reads_vcf_bgzf_and_bcf_from_a_named_pipe() {
        let vcf = long_vcf();
        let want: Vec<i64> = (1..=5000).collect();
        for (output_type, data) in [
            (OutputType::Vcf, vcf.clone()),
            (OutputType::VcfBgzf, convert(&vcf, OutputType::VcfBgzf)),
            (OutputType::Bcf, convert(&vcf, OutputType::Bcf)),
        ] {
            let same = positions_through_pipe(data).map(|got| got == want);
            assert!(matches!(same, Ok(true)), "{output_type:?}: {same:?}");
        }
    }

    #[test]
    fn tells_streams_from_files() {
        let (rx, _tx) = io::pipe().expect("create a pipe");
        assert!(is_stream("-"));
        assert!(is_stream(&format!("/dev/fd/{}", rx.as_raw_fd())));
        let dir = std::env::temp_dir().join(format!("pasta-stream-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let file = dir.join("x.vcf");
        std::fs::write(&file, b"").expect("write the test file");
        let file = file.to_string_lossy();
        let streams = [
            is_stream(&file),
            is_stream(&format!("{file}##idx##x.csi")),
            is_stream(&dir.to_string_lossy()),
        ];
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(streams, [false; 3]);
    }

    /// Writes `text` as plain gzip members, one per part, the CRC32 of the last one wrong.
    fn gzip_with_bad_last_crc(parts: &[&[u8]]) -> Vec<u8> {
        let mut data = Vec::new();
        for part in parts {
            let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(0));
            gz.write_all(part).expect("compress");
            data.extend_from_slice(&gz.finish().expect("finish"));
        }
        let crc = data.len() - 8;
        data[crc] ^= 0xff;
        data
    }

    #[test]
    fn plain_gzip_checks_each_member_and_a_line_that_does_not_parse_reads_the_rest() {
        let dir = std::env::temp_dir().join(format!("pasta-gzip-crc-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let header = b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let read_all = |name: &str, parts: &[&[u8]]| {
            let path = dir.join(name);
            std::fs::write(&path, gzip_with_bad_last_crc(parts)).expect("write the test file");
            let mut reader = Reader::open(&path.to_string_lossy(), 0).expect("open");
            assert_eq!(reader.compression(), Compression::Gzip);
            let mut header = reader.read_header().expect("the header");
            let (mut line, mut rec, mut n) = (Vec::new(), Record::new(), 0);
            let end = loop {
                match reader.read_record(&mut header, &mut line, &mut rec) {
                    Ok(true) => n += 1,
                    other => break other,
                }
            };
            // A failed read stays failed, as HTSlib's error does.
            let again = reader.read_record(&mut header, &mut line, &mut rec);
            (n, end, again)
        };
        // The second member's CRC32 is wrong: the records before its end are read, then the
        // read fails.
        let (n, end, again) = read_all("crc.vcf.gz", &[header, b"chr1\t1\t.\tA\tC\t.\t.\t.\n"]);
        assert_eq!(n, 1);
        assert!(matches!(end, Err(ReadError::Io(_))), "{end:?}");
        assert!(matches!(again, Err(ReadError::Io(_))), "{again:?}");
        // A line that does not parse in a member whose CRC32 is wrong: a read failure, not a
        // parse error.
        let bad = b"chr1\tx\t.\tA\tC\t.\t.\t.\nchr1\t3\t.\tA\tC\t.\t.\t.\n";
        let (n, end, _) = read_all(
            "line.vcf.gz",
            &[header, b"chr1\t1\t.\tA\tC\t.\t.\t.\n", bad],
        );
        assert_eq!(n, 1);
        assert!(matches!(end, Err(ReadError::Io(_))), "{end:?}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_pipe_is_not_checked_for_the_eof_marker() {
        let (rx, tx) = io::pipe().expect("create a pipe");
        let path = format!("/dev/fd/{}", rx.as_raw_fd());
        assert_eq!(has_bgzf_eof(Path::new(&path)).expect("no I/O error"), None);
        drop(tx);
    }

    #[test]
    fn a_fifo_without_a_writer_is_read_again_as_empty() {
        let dir = std::env::temp_dir().join(format!("pasta-fifo-again-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let fifo = dir.join("fifo");
        let c = std::ffi::CString::new(fifo.to_string_lossy().as_bytes()).expect("no NUL");
        // SAFETY: `c` is a NUL-terminated path, as mkfifo requires.
        assert_eq!(unsafe { libc::mkfifo(c.as_ptr(), 0o600) }, 0, "mkfifo");
        let name = fifo.to_string_lossy().into_owned();
        // Opening the FIFO with `open` would wait for a writer forever.
        let (tx, rx) = std::sync::mpsc::channel();
        thread::spawn(move || tx.send(Reader::open_again(&name, 0).map(drop)));
        let again = rx.recv_timeout(std::time::Duration::from_secs(10));
        let _ = std::fs::remove_dir_all(&dir);
        let again = again.expect("returns instead of waiting for a writer");
        assert!(matches!(again, Err(ReadError::UnknownFormat)), "{again:?}");
    }

    #[test]
    fn a_pipe_read_again_goes_on_where_the_first_reading_stopped() {
        // A bgzipped VCF much longer than the first reading takes (its data compresses badly).
        let mut vcf = long_vcf();
        let mut x = 1u64;
        for pos in 5001..=100_000 {
            x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let line = format!("chr1\t{pos}\tr{x}\tA\tC\t.\t.\t.\n");
            vcf.extend_from_slice(line.as_bytes());
        }
        let data = convert(&vcf, OutputType::VcfBgzf);
        assert!(data.len() > 8 * RAW_PEEK, "{} bytes", data.len());
        let (rx, mut tx) = io::pipe().expect("create a pipe");
        // The writer stops with a broken pipe once both readings are done.
        let feeder = thread::spawn(move || tx.write_all(&data));
        let path = format!("/dev/fd/{}", rx.as_raw_fd());
        let mut first = Reader::open(&path, 0).expect("the pipe opens");
        first.read_header().expect("the header is read");
        drop(first);
        // HTSlib reads the rest too: BGZF data cut mid-way is binary it does not know.
        let again = Reader::open_again(&path, 0).map(drop);
        drop(rx);
        let _ = feeder.join().expect("the writer thread does not panic");
        assert!(matches!(again, Err(ReadError::UnknownBinary)), "{again:?}");
    }

    #[test]
    fn tells_unknown_binary_data_from_text_and_known_formats() {
        assert!(is_unknown_binary(b"\x00\x01\x02binary"));
        for known in [
            &b""[..],
            b"chr1\t1\ttext\r\n",
            b"BCF\x02\x02\x00",
            b"BAM\x01\x00",
            b"CRAM\x03",
        ] {
            assert!(!is_unknown_binary(known), "{known:?}");
        }
        // HTSlib looks at the first 1024 bytes only.
        let mut late = vec![b'a'; PEEK];
        late.push(0);
        assert!(!is_unknown_binary(&late));
    }
}
