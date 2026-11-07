//! Writing VCF and BCF (bcftools' `-O v|z|b|u`), with BGZF compression on worker threads and
//! indexing while writing (`--write-index`).
//!
//! Only [`Writer::finish`] completes an output. A writer dropped before that is abandoned, as
//! [`Writer::abandon`] abandons it, so that an early return on an error leaves the output as
//! bcftools' `error()` leaves it, with or without compression threads.

use std::{
    collections::VecDeque,
    fs::File,
    io::{self, Write},
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

use super::bgzf_out::BlockWriter;
use crate::bgzf::io::writer::CompressionLevel;
use crate::variant::{
    Header, HeaderTextFlavor, Record, bcf,
    index::{
        Index, IndexFormat, PushCheck, PushError, TBI_MIN_SHIFT, TBI_N_LVLS, adjust_csi_settings,
        tabix,
    },
    vcf::{self, FormatError},
};

/// An output type (bcftools' `-O`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum OutputType {
    /// Uncompressed VCF (`-Ov`).
    Vcf,
    /// BGZF-compressed VCF (`-Oz`).
    VcfBgzf,
    /// BGZF-compressed BCF (`-Ob`).
    Bcf,
    /// BCF in uncompressed BGZF blocks (`-Ou`).
    BcfUncompressed,
}

impl OutputType {
    /// Returns whether the output is VCF text.
    #[must_use]
    pub fn is_vcf(self) -> bool {
        matches!(self, Self::Vcf | Self::VcfBgzf)
    }

    /// Returns whether the output is BGZF.
    #[must_use]
    pub fn is_bgzf(self) -> bool {
        !matches!(self, Self::Vcf)
    }
}

/// Options for [`Writer::create`].
#[derive(Clone, Debug)]
pub struct WriteOptions {
    /// The output type.
    pub output_type: OutputType,
    /// The zlib-style compression level, 0-9, mapped to the compressor's level as HTSlib maps
    /// it; `None` for pasta's default, zflate's fast level 5 (see [`CompressionLevel`]).
    ///
    /// [`CompressionLevel`]: crate::bgzf::io::writer::CompressionLevel
    pub level: Option<u8>,
    /// BGZF compression threads (0 compresses on the calling thread).
    pub threads: usize,
}

/// An index to build while writing (bcftools' `--write-index`).
#[derive(Clone, Debug)]
pub struct IndexRequest {
    /// The index file to write.
    pub path: PathBuf,
    /// The minimum shift: 0 for TBI (VCF only), otherwise CSI.
    pub min_shift: u32,
}

/// An error writing a record.
#[derive(Debug)]
#[non_exhaustive]
pub enum WriteError {
    /// An I/O error.
    Io(io::Error),
    /// An I/O error on output handed to the compression threads before [`Writer::finish`],
    /// which they reported only then. It is a failed write of records, not of closing: a
    /// serial run (and bcftools, whose threads notice it sooner) reports it while writing.
    Deferred(io::Error),
    /// The record could not be encoded as BCF.
    Encode(bcf::EncodeError),
    /// The record could not be formatted as VCF.
    Format(FormatError),
    /// The record could not be indexed.
    Index(PushError),
    /// The output was written, but its index could not be saved (HTSlib's `bcf_idx_save`
    /// failing). Displayed as the I/O error.
    SaveIndex(io::Error),
}

impl std::fmt::Display for WriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) | Self::Deferred(e) => write!(f, "{e}"),
            Self::Encode(e) => write!(f, "{e}"),
            Self::Format(e) => write!(f, "{e}"),
            Self::Index(e) => write!(f, "{e}"),
            Self::SaveIndex(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for WriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) | Self::Deferred(e) => Some(e),
            Self::Encode(e) => Some(e),
            Self::Format(e) => Some(e),
            Self::Index(e) => Some(e),
            Self::SaveIndex(e) => Some(e),
        }
    }
}

impl From<io::Error> for WriteError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

/// The position and interval of one record in pre-formatted output, for indexing.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct RecordSpan {
    /// The record's length in bytes.
    pub len: u32,
    /// The contig ID.
    pub rid: i32,
    /// The 0-based position.
    pub pos: i64,
    /// The reference length.
    pub rlen: i64,
}

type Sink = Box<dyn Write + Send>;

/// A file buffer with HTSlib's `hwrite` semantics: writes that fit are buffered, others flush
/// the buffer and go out directly when at least half the capacity. The capacity is HTSlib's:
/// 128 KiB for pipes, the file system block size otherwise. Matching it makes output lost on
/// a fatal error (which bcftools exits on without flushing) the same as bcftools'.
struct HFile {
    inner: Sink,
    buf: Vec<u8>,
    cap: usize,
}

impl HFile {
    fn new(inner: Sink, cap: usize) -> Self {
        Self {
            inner,
            buf: Vec::with_capacity(cap),
            cap,
        }
    }
}

impl HFile {
    /// Appends what `fill` writes (one record) as [`Write::write_all`] would, without first
    /// formatting it into a separate buffer: the record goes straight into the buffer, and the
    /// buffer is then written out exactly as `hwrite` would have written it. On an error from
    /// `fill`, its partial output is dropped.
    fn append_with<E: From<io::Error>>(
        &mut self,
        fill: impl FnOnce(&mut Vec<u8>) -> Result<(), E>,
    ) -> Result<(), E> {
        let before = self.buf.len();
        if let Err(e) = fill(&mut self.buf) {
            self.buf.truncate(before);
            return Err(e);
        }
        let len = self.buf.len() - before;
        if before == 0 {
            // hwrite2 with an empty buffer: written directly when it fills the buffer.
            if len >= self.cap {
                self.inner.write_all(&self.buf)?;
                self.buf.clear();
            }
            return Ok(());
        }
        // hwrite keeps a buffer it fills exactly, and writes it on the next call.
        if len <= self.cap.saturating_sub(before) {
            return Ok(());
        }
        // The buffer was topped up to its capacity and written; then the rest.
        self.inner.write_all(&self.buf[..self.cap])?;
        let rest = self.buf.len() - self.cap;
        if rest * 2 >= self.cap {
            self.inner.write_all(&self.buf[self.cap..])?;
            self.buf.clear();
        } else {
            self.buf.copy_within(self.cap.., 0);
            self.buf.truncate(rest);
        }
        Ok(())
    }
}

impl Write for HFile {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.write_all(data)?;
        Ok(data.len())
    }

    fn write_all(&mut self, data: &[u8]) -> io::Result<()> {
        // `hwrite`: fill the buffer with what fits (unless it is empty and the data does
        // not fit), then `hwrite2` flushes it and writes or buffers the rest.
        let room = self.cap - self.buf.len();
        let copied = if data.len() >= room && self.buf.is_empty() {
            0
        } else {
            let n = room.min(data.len());
            self.buf.extend_from_slice(&data[..n]);
            if n == data.len() {
                return Ok(());
            }
            n
        };
        let rest = &data[copied..];
        self.inner.write_all(&self.buf)?;
        self.buf.clear();
        if rest.len() * 2 >= self.cap {
            self.inner.write_all(rest)?;
        } else {
            self.buf.extend_from_slice(rest);
        }
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.write_all(&self.buf)?;
        self.buf.clear();
        self.inner.flush()
    }
}

/// HTSlib's `hFILE` buffer size for a file descriptor (`blksize`).
fn hfile_capacity(file: &File) -> usize {
    const MAX: usize = 128 * 1024;
    #[cfg(unix)]
    {
        use std::os::unix::fs::{FileTypeExt, MetadataExt};
        match file.metadata() {
            Ok(m) if m.file_type().is_fifo() => MAX,
            Ok(m) if m.blksize() > 0 => m.blksize() as usize,
            _ => MAX,
        }
    }
    #[cfg(not(unix))]
    {
        let _ = file;
        MAX
    }
}

enum Output {
    Plain(HFile),
    Bgzf(BlockWriter<Sink>),
}

#[derive(Clone, Copy)]
struct Pending {
    tid: i32,
    beg: i64,
    end: i64,
    start: (u64, u16),
    stop: (u64, u16),
}

struct OtfIndex {
    index: Index,
    /// The checks of each record as it is written, which the index pushes later.
    check: PushCheck,
    format: IndexFormat,
    path: PathBuf,
    is_vcf: bool,
    pending: VecDeque<Pending>,
}

/// The number of writers dropped without [`Writer::finish`] or [`Writer::abandon`].
static DROPPED_UNFINISHED: AtomicUsize = AtomicUsize::new(0);

/// Returns how many writers were dropped without [`Writer::finish`] or [`Writer::abandon`],
/// which abandoned their output. A program that succeeded has dropped none; any other count
/// is a missing `finish`, which truncates the output.
///
/// The count is for the whole process and never reset: it is a diagnostic for a program that
/// runs one command, such as a debug assertion before it exits.
#[must_use]
pub fn dropped_unfinished() -> usize {
    DROPPED_UNFINISHED.load(Ordering::Relaxed)
}

/// Counts a writer dropped unfinished, unless [`Writer::finish`] or [`Writer::abandon`] forgets
/// it first.
struct Unfinished;

impl Drop for Unfinished {
    fn drop(&mut self) {
        DROPPED_UNFINISHED.fetch_add(1, Ordering::Relaxed);
    }
}

/// A VCF or BCF output (HTSlib's `htsFile` opened for writing).
///
/// Dropping it before [`Self::finish`] abandons the output, as [`Self::abandon`] does.
pub struct Writer {
    out: Output,
    output_type: OutputType,
    index_request: Option<IndexRequest>,
    index: Option<OtfIndex>,
    buf: Vec<u8>,
    unfinished: Unfinished,
}

impl std::fmt::Debug for Writer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Writer")
            .field("output_type", &self.output_type)
            .finish_non_exhaustive()
    }
}

/// The compression level for a zlib-style level, as HTSlib's `bgzf_compress` maps it onto
/// libdeflate's, which zflate keeps: 1-9 become `{1,2,3,5,6,7,8,10,12}`, and 0 writes stored
/// blocks. Without a level (`None`), [`CompressionLevel::VARIANTS`], level 6, where HTSlib's
/// default compresses at 7.
pub(crate) fn htslib_level(level: Option<u8>) -> CompressionLevel {
    match level {
        Some(level) => {
            let level = [0, 1, 2, 3, 5, 6, 7, 8, 10, 12][usize::from(level.min(9))];
            CompressionLevel::new(level).unwrap_or_default()
        }
        None => CompressionLevel::VARIANTS,
    }
}

/// Standard output, written to its file descriptor without Rust's line buffering. Like
/// HTSlib, it writes to descriptor 1 itself rather than to a copy, which would take one more
/// descriptor.
#[cfg(unix)]
struct Stdout(std::mem::ManuallyDrop<File>);

#[cfg(unix)]
impl Stdout {
    fn new() -> Self {
        use std::os::fd::FromRawFd;
        // SAFETY: descriptor 1 is open for the life of the process (Rust's runtime opens
        // /dev/null in its place if it is closed at startup), and `ManuallyDrop` never closes
        // it, so the `File` only borrows it.
        Self(std::mem::ManuallyDrop::new(unsafe {
            File::from_raw_fd(libc::STDOUT_FILENO)
        }))
    }
}

#[cfg(unix)]
impl Write for Stdout {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }

    fn write_vectored(&mut self, bufs: &[io::IoSlice<'_>]) -> io::Result<usize> {
        self.0.write_vectored(bufs)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

/// Standard output and the size of HTSlib's buffer for it.
fn stdout_sink() -> io::Result<(Sink, usize)> {
    #[cfg(unix)]
    {
        let out = Stdout::new();
        let cap = hfile_capacity(&out.0);
        Ok((Box::new(out), cap))
    }
    #[cfg(not(unix))]
    {
        Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "standard output",
        ))
    }
}

impl Writer {
    /// Creates an output file, or standard output for `None` or `-`.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be created, or when it is a regular file opened as
    /// an input ([`crate::util::inputs::check_output`]), which is left untouched.
    pub fn create(path: Option<&Path>, options: &WriteOptions) -> io::Result<Self> {
        let (sink, cap): (Sink, usize) = match path {
            Some(p) if p != Path::new("-") => {
                crate::util::inputs::check_output(p)?;
                let file = File::create(p)?;
                let cap = hfile_capacity(&file);
                (Box::new(file), cap)
            }
            _ => stdout_sink()?,
        };
        Ok(Self::with_capacity(sink, options, cap))
    }

    /// Creates an output over any writer.
    #[must_use]
    pub fn from_sink(sink: Box<dyn Write + Send>, options: &WriteOptions) -> Self {
        Self::with_capacity(sink, options, 1 << 17)
    }

    fn with_capacity(sink: Sink, options: &WriteOptions, cap: usize) -> Self {
        let out = match options.output_type {
            OutputType::Vcf => Output::Plain(HFile::new(sink, cap)),
            t => {
                let level = if t == OutputType::BcfUncompressed {
                    CompressionLevel::NONE
                } else {
                    htslib_level(options.level)
                };
                Output::Bgzf(BlockWriter::new(
                    Box::new(HFile::new(sink, cap)),
                    level,
                    options.threads,
                ))
            }
        };
        Self {
            out,
            output_type: options.output_type,
            index_request: None,
            index: None,
            buf: Vec::with_capacity(1 << 16),
            unfinished: Unfinished,
        }
    }

    /// Returns the output type.
    #[must_use]
    pub fn output_type(&self) -> OutputType {
        self.output_type
    }

    /// Requests an index built while writing; call before [`Self::write_header`].
    ///
    /// # Errors
    ///
    /// Returns an error for uncompressed VCF output (HTSlib: "Indexing is only supported on
    /// BGZF-compressed files").
    pub fn request_index(&mut self, request: IndexRequest) -> io::Result<()> {
        if !self.output_type.is_bgzf() {
            crate::variant::log::error(
                "bcf_idx_init",
                format_args!("Indexing is only supported on BGZF-compressed files"),
            );
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cannot index uncompressed output",
            ));
        }
        self.index_request = Some(request);
        Ok(())
    }

    /// Writes the header (HTSlib's `bcf_hdr_write`), ending its BGZF block, and starts the
    /// requested index.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        let data = if self.output_type.is_vcf() {
            let mut text = Vec::new();
            header.format_into(&mut text, HeaderTextFlavor::Vcf);
            text
        } else {
            bcf::encode_header(header)
        };
        match &mut self.out {
            Output::Plain(w) => w.write_all(&data)?,
            Output::Bgzf(w) => {
                w.write_all(&data)?;
                w.flush_block()?;
            }
        }
        self.start_index(header)
    }

    /// Starts the requested index at the current position, if it has not started yet. Called
    /// by [`Self::write_header`]; call it directly when no header is written.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn start_index(&mut self, header: &Header) -> io::Result<()> {
        if let Some(request) = self.index_request.take() {
            self.init_index(header, request)?;
        }
        Ok(())
    }

    fn init_index(&mut self, header: &Header, request: IndexRequest) -> io::Result<()> {
        let Output::Bgzf(w) = &mut self.out else {
            unreachable!("indexes need BGZF output");
        };
        let (block, _) = w.tell();
        let offset0 = w.wait_address(block)? << 16;
        let max_len = (0..header.contig_count() as u32)
            .filter_map(|i| header.contig_length(i))
            .max()
            .unwrap_or(0);
        let max_len = if max_len == 0 {
            (1i64 << 31) - 1
        } else {
            max_len
        };
        let is_vcf = self.output_type.is_vcf();
        let (index, format) = if is_vcf {
            if request.min_shift == 0 {
                let mut idx = Index::new(0, IndexFormat::Tbi, offset0, TBI_MIN_SHIFT, TBI_N_LVLS);
                idx.set_meta(tabix::vcf_meta());
                (idx, IndexFormat::Tbi)
            } else {
                let mut min_shift = request.min_shift;
                let mut n_lvls = (31 - min_shift).div_ceil(3);
                adjust_csi_settings(max_len, &mut min_shift, &mut n_lvls);
                let mut idx = Index::new(0, IndexFormat::Csi, offset0, min_shift, n_lvls);
                idx.set_meta(tabix::vcf_meta());
                (idx, IndexFormat::Csi)
            }
        } else {
            let mut min_shift = if request.min_shift == 0 {
                14
            } else {
                request.min_shift
            };
            let mut n_lvls = 0;
            adjust_csi_settings(max_len, &mut min_shift, &mut n_lvls);
            let nids = header.contigs.iter().filter(|c| c.is_some()).count();
            (
                Index::new(nids, IndexFormat::Csi, offset0, min_shift, n_lvls),
                IndexFormat::Csi,
            )
        };
        self.index = Some(OtfIndex {
            index,
            check: PushCheck::new(),
            format,
            path: request.path,
            is_vcf,
            pending: VecDeque::new(),
        });
        Ok(())
    }

    /// Writes a record (HTSlib's `bcf_write`): formatted for VCF output, synced and encoded for
    /// BCF output.
    ///
    /// # Errors
    ///
    /// Returns an error (already reported as HTSlib reports it) when the record does not match
    /// the header or cannot be written.
    pub fn write_record(&mut self, header: &Header, rec: &mut Record) -> Result<(), WriteError> {
        bcf::check_sample_count(header, rec).map_err(WriteError::Encode)?;
        if let (Output::Plain(out), None) = (&mut self.out, &self.index) {
            return out.append_with(|buf| {
                vcf::format_record(header, rec, buf).map_err(WriteError::Format)
            });
        }
        let mut buf = std::mem::take(&mut self.buf);
        buf.clear();
        let result = if self.output_type.is_vcf() {
            vcf::format_record(header, rec, &mut buf).map_err(WriteError::Format)
        } else {
            bcf::encode_record(header, rec, &mut buf).map_err(WriteError::Encode)
        };
        let result = result.and_then(|()| {
            let span = RecordSpan {
                len: buf.len() as u32,
                rid: rec.rid(),
                pos: rec.pos(),
                rlen: rec.rlen(),
            };
            self.write_spans(header, &buf, &[span])
        });
        self.buf = buf;
        result
    }

    /// Writes records formatted or encoded elsewhere (for example on worker threads), with
    /// their spans for indexing.
    ///
    /// # Errors
    ///
    /// Returns an I/O or indexing error.
    pub fn write_spans(
        &mut self,
        header: &Header,
        data: &[u8],
        spans: &[RecordSpan],
    ) -> Result<(), WriteError> {
        let Some(index) = &mut self.index else {
            match &mut self.out {
                Output::Plain(w) => w.write_all(data)?,
                // Compressed VCF blocks end at records (`vcf_write`), BCF blocks are full
                // (`bcf_write`).
                Output::Bgzf(w) if self.output_type.is_vcf() => {
                    w.write_records(data, spans.iter().map(|s| s.len as usize))?;
                }
                Output::Bgzf(w) => w.write_all(data)?,
            }
            return Ok(());
        };
        let Output::Bgzf(w) = &mut self.out else {
            unreachable!("indexes need BGZF output");
        };
        let mut off = 0usize;
        for span in spans {
            let rec = &data[off..off + span.len as usize];
            off += span.len as usize;
            if index.is_vcf {
                w.flush_try(rec.len())?;
            }
            let start = w.tell();
            w.write_all(rec)?;
            let stop = w.tell();
            let tid = if index.is_vcf {
                let name = header.contig_name(span.rid).unwrap_or(b"(unknown)");
                index.index.tbi_name(span.rid, name)
            } else {
                span.rid
            };
            let (beg, end) = (span.pos, span.pos + span.rlen);
            index
                .check
                .check(&index.index, tid, beg, end)
                .map_err(WriteError::Index)?;
            index.pending.push_back(Pending {
                tid,
                beg,
                end,
                start,
                stop,
            });
        }
        resolve(index, w, Resolve::Written)
    }

    /// Writes bytes as they are (into the BGZF stream for compressed output), for callers that
    /// copy header text or records without parsing them.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn write_raw(&mut self, data: &[u8]) -> io::Result<()> {
        match &mut self.out {
            Output::Plain(w) => w.write_all(data),
            Output::Bgzf(w) => w.write_all(data),
        }
    }

    /// Ends the current BGZF block and appends complete BGZF frames without recompressing them
    /// (HTSlib's `bgzf_flush` and `bgzf_raw_write`). Not for indexed or uncompressed output.
    ///
    /// # Errors
    ///
    /// Returns an error for uncompressed or indexed output, or an I/O error.
    pub fn write_bgzf_frames(&mut self, frames: &[u8]) -> io::Result<()> {
        match &mut self.out {
            Output::Bgzf(w) if self.index.is_none() => w.write_frames(frames),
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "raw BGZF frames need unindexed BGZF output",
            )),
        }
    }

    /// Drops the writer without flushing buffered output, as bcftools' `error()` exits. The
    /// full BGZF blocks handed to compression threads are written first, so that the output is
    /// the same as without threads, which write each block as it fills; the block being filled
    /// and the EOF marker are not, nor is the index. Dropping the writer does the same; this
    /// says that it is meant.
    pub fn abandon(self) {
        let Self { unfinished, .. } = self;
        std::mem::forget(unfinished);
    }

    /// Flushes buffered output, writes the BGZF EOF marker, and saves the index.
    ///
    /// # Errors
    ///
    /// Returns an I/O or indexing error, [`WriteError::Deferred`] when the compression threads
    /// failed on a block written before this call, or [`WriteError::SaveIndex`] when the output
    /// was written but its index could not be saved.
    pub fn finish(self) -> Result<(), WriteError> {
        let Self {
            out,
            index,
            unfinished,
            ..
        } = self;
        std::mem::forget(unfinished);
        match out {
            Output::Plain(mut w) => {
                w.flush()?;
                Ok(())
            }
            Output::Bgzf(mut w) => {
                // The blocks before this one were handed over before closing.
                let (queued, _) = w.tell();
                finish_bgzf(&mut w, index).map_err(|e| match e {
                    WriteError::Io(e) if w.failed_block().is_some_and(|b| b < queued) => {
                        WriteError::Deferred(e)
                    }
                    e => e,
                })
            }
        }
    }
}

/// Flushes the last block, writes the EOF marker, and saves the index.
fn finish_bgzf(w: &mut BlockWriter<Sink>, index: Option<OtfIndex>) -> Result<(), WriteError> {
    w.flush_block()?;
    let Some(mut index) = index else {
        w.finish()?;
        return Ok(());
    };
    let (block, _) = w.tell();
    let end = w.wait_address(block)? << 16;
    resolve(&mut index, w, Resolve::All)?;
    index.index.amend_last(end);
    index.index.finish(end);
    w.finish()?;
    index
        .index
        .save(&index.path, index.format)
        .map_err(WriteError::SaveIndex)
}

/// Which pending index entries [`resolve`] pushes.
#[derive(Clone, Copy, Debug)]
enum Resolve {
    /// The entries whose blocks are already written, so that their offsets are known.
    Written,
    /// All of them, waiting for their blocks to be compressed.
    All,
}

/// Pushes pending index entries to the index.
fn resolve(
    index: &mut OtfIndex,
    w: &mut BlockWriter<Sink>,
    which: Resolve,
) -> Result<(), WriteError> {
    while let Some(p) = index.pending.front().copied() {
        let need = p.start.0.max(p.stop.0);
        let addr = match which {
            Resolve::All => Some(w.wait_address(need)?),
            Resolve::Written => w.address(need)?,
        };
        if addr.is_none() {
            break;
        }
        let voff = |w: &mut BlockWriter<Sink>, pos: (u64, u16)| -> io::Result<u64> {
            Ok((w.wait_address(pos.0)? << 16) | u64::from(pos.1))
        };
        let start = voff(w, p.start)?;
        let stop = voff(w, p.stop)?;
        index.index.amend_last(start);
        index
            .index
            .push(p.tid, p.beg, p.end, stop)
            .map_err(WriteError::Index)?;
        index.pending.pop_front();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::bgzf_out::{
        BLOCK_SIZE,
        tests::{FailAfter, within_timeout},
    };
    use super::*;

    #[test]
    fn a_threaded_indexed_output_reports_a_write_error() {
        let header = Header::parse(b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .expect("valid header");
        let mut data = Vec::new();
        let mut spans = Vec::new();
        for pos in 0..50_000 {
            let before = data.len();
            data.extend_from_slice(format!("chr1\t{}\t.\tA\tC\t.\t.\t.\n", pos + 1).as_bytes());
            spans.push(RecordSpan {
                len: (data.len() - before) as u32,
                rid: 0,
                pos,
                rlen: 1,
            });
        }
        for threads in [1, 2, 4, 8] {
            let header = header.clone();
            let (data, spans) = (data.clone(), spans.clone());
            // Where the error was reported: by a write, or by `finish`.
            let result = within_timeout(move || -> Result<(), (&'static str, WriteError)> {
                let options = WriteOptions {
                    output_type: OutputType::VcfBgzf,
                    level: Some(0),
                    threads,
                };
                // The error comes on the second block, with many blocks after it.
                let sink = Box::new(FailAfter {
                    left: BLOCK_SIZE + 100,
                });
                let mut writer = Writer::from_sink(sink, &options);
                let writing = |e: WriteError| ("write", e);
                writer
                    .request_index(IndexRequest {
                        path: PathBuf::from("never-written.csi"),
                        min_shift: 14,
                    })
                    .map_err(|e| writing(e.into()))?;
                writer
                    .write_header(&header)
                    .map_err(|e| writing(e.into()))?;
                let mut off = 0;
                for chunk in spans.chunks(1000) {
                    let len: usize = chunk.iter().map(|s| s.len as usize).sum();
                    writer
                        .write_spans(&header, &data[off..off + len], chunk)
                        .map_err(writing)?;
                    off += len;
                }
                writer.finish().map_err(|e| ("finish", e))
            });
            // A block written before `finish` failed: `finish` reports it as deferred, not as
            // a failure to close.
            let full = |e: &io::Error| e.kind() == io::ErrorKind::StorageFull;
            assert!(
                matches!(&result, Err(("write", WriteError::Io(e)) | ("finish", WriteError::Deferred(e))) if full(e)),
                "threads {threads}: {result:?}"
            );
        }
    }

    /// A writer whose bytes stay readable after it is handed over.
    #[derive(Clone, Default)]
    struct Shared(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

    impl Write for Shared {
        fn write(&mut self, data: &[u8]) -> io::Result<usize> {
            self.0.lock().expect("not poisoned").extend_from_slice(data);
            Ok(data.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn abandoning_writes_the_same_blocks_with_threads() {
        // Blocks that do not compress, so that the file buffer is written several times; the
        // block being filled at the end is lost, with what the file buffer holds.
        let mut x = 12345u32;
        let data: Vec<u8> = (0..20 * BLOCK_SIZE + 100)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                (x >> 16) as u8
            })
            .collect();
        let abandoned = |threads| {
            let data = data.clone();
            within_timeout(move || {
                let sink = Shared::default();
                let options = WriteOptions {
                    output_type: OutputType::Bcf,
                    level: None,
                    threads,
                };
                let mut writer = Writer::from_sink(Box::new(sink.clone()), &options);
                writer.write_raw(&data).expect("write to memory");
                writer.abandon();
                sink.0.lock().expect("not poisoned").clone()
            })
        };
        let serial = abandoned(0);
        assert!(!serial.is_empty(), "the file buffer was written");
        assert!(
            serial.len() < data.len(),
            "the last block and buffer are lost"
        );
        for threads in [1, 2, 4] {
            let threaded = abandoned(threads);
            assert!(
                threaded == serial,
                "threads {threads}: {} bytes, {} without threads",
                threaded.len(),
                serial.len()
            );
        }
    }

    #[test]
    fn dropping_an_unfinished_writer_abandons_it() {
        // Full blocks, then one being filled when the writer is dropped; and a writer dropped
        // before anything was written to it (an error found after the output was opened).
        let data = vec![b'x'; 3 * BLOCK_SIZE + 100];
        for output_type in [OutputType::VcfBgzf, OutputType::Bcf, OutputType::Vcf] {
            for len in [0, data.len()] {
                let data = data[..len].to_vec();
                let written = move |threads, abandon: bool| {
                    let data = data.clone();
                    within_timeout(move || {
                        let sink = Shared::default();
                        let options = WriteOptions {
                            output_type,
                            level: Some(0),
                            threads,
                        };
                        let mut writer = Writer::from_sink(Box::new(sink.clone()), &options);
                        writer.write_raw(&data).expect("write to memory");
                        let before = dropped_unfinished();
                        if abandon {
                            writer.abandon();
                        } else {
                            drop(writer);
                            assert!(dropped_unfinished() > before, "the drop is counted");
                            // Whatever threads left running would write comes before they
                            // drop the sink.
                            while std::sync::Arc::strong_count(&sink.0) > 1 {
                                std::thread::yield_now();
                            }
                        }
                        sink.0.lock().expect("not poisoned").clone()
                    })
                };
                let serial = written(0, true);
                if len == 0 {
                    assert!(serial.is_empty(), "{output_type:?}: no EOF marker");
                }
                for threads in [0, 1, 2, 4] {
                    let dropped = written(threads, false);
                    assert!(
                        dropped == serial,
                        "{output_type:?}, {len} bytes, threads {threads}: {} bytes, {} abandoned",
                        dropped.len(),
                        serial.len()
                    );
                }
            }
        }
    }

    #[test]
    fn an_unsorted_record_stops_an_indexed_output_where_a_serial_run_does() {
        let header = Header::parse(b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .expect("valid header");
        const N: usize = 50_000;
        // Records with one out of order (unsorted, or back on a contig left before), in the
        // middle or last; each case is written one record at a time and in batches.
        let records = |bad: usize, back_to_chr1: bool| -> (Vec<u8>, Vec<RecordSpan>) {
            let mut data = Vec::new();
            let mut spans = Vec::new();
            for i in 0..N {
                let (rid, pos) = match i {
                    _ if i == bad && back_to_chr1 => (0, 5),
                    _ if i == bad => (1, 10),
                    _ if i < N / 4 => (0, i as i64),
                    _ => (1, i as i64),
                };
                let before = data.len();
                data.extend_from_slice(
                    format!("chr{}\t{}\t.\tA\tC\t.\t.\t.\n", rid + 1, pos + 1).as_bytes(),
                );
                spans.push(RecordSpan {
                    len: (data.len() - before) as u32,
                    rid,
                    pos,
                    rlen: 1,
                });
            }
            (data, spans)
        };
        // A VCF index names a contig again when it comes back (`hts_idx_tbi_name`), so only
        // BCF finds a contig that is not contiguous.
        let cases = [
            (OutputType::VcfBgzf, N / 2, false),
            (OutputType::VcfBgzf, N - 1, false),
            (OutputType::Bcf, N / 2, false),
            (OutputType::Bcf, N - 1, false),
            (OutputType::Bcf, N / 2, true),
        ];
        for (output_type, bad, back_to_chr1) in cases {
            let (data, spans) = records(bad, back_to_chr1);
            let written = |threads, batch: usize| {
                let (header, data, spans) = (header.clone(), data.clone(), spans.clone());
                within_timeout(move || {
                    let sink = Shared::default();
                    let options = WriteOptions {
                        output_type,
                        level: Some(0),
                        threads,
                    };
                    let mut writer = Writer::from_sink(Box::new(sink.clone()), &options);
                    writer
                        .request_index(IndexRequest {
                            path: PathBuf::from("never-written.csi"),
                            min_shift: 14,
                        })
                        .expect("BGZF output");
                    writer.write_header(&header).expect("write to memory");
                    let mut off = 0;
                    let failed = spans.chunks(batch).any(|chunk| {
                        let len: usize = chunk.iter().map(|s| s.len as usize).sum();
                        let r = writer.write_spans(&header, &data[off..off + len], chunk);
                        off += len;
                        r.is_err()
                    });
                    assert!(failed, "the bad record fails");
                    writer.abandon();
                    sink.0.lock().expect("not poisoned").clone()
                })
            };
            let serial = written(0, 1);
            assert!(
                !serial.is_empty(),
                "blocks before the bad record are written"
            );
            for threads in [0, 1, 2, 4] {
                for batch in [1, 1000] {
                    let got = written(threads, batch);
                    assert!(
                        got == serial,
                        "{output_type:?}, bad record {bad}, threads {threads}, batch {batch}: \
                         {} bytes, {} serially",
                        got.len(),
                        serial.len()
                    );
                }
            }
        }
    }

    /// Returns the uncompressed sizes of the BGZF blocks of `data`, without the EOF marker.
    fn block_sizes(mut data: &[u8]) -> Vec<usize> {
        let mut sizes = Vec::new();
        while !data.is_empty() {
            let bsize = usize::from(u16::from_le_bytes([data[16], data[17]])) + 1;
            let isize = u32::from_le_bytes(data[bsize - 4..bsize].try_into().expect("4 bytes"));
            sizes.push(isize as usize);
            data = &data[bsize..];
        }
        assert_eq!(sizes.pop(), Some(0), "the EOF marker ends the output");
        sizes
    }

    #[test]
    fn compressed_vcf_blocks_end_at_records() {
        let header = Header::parse(b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .expect("valid header");
        // Records of many lengths, two of them longer than a block.
        let mut data = Vec::new();
        let mut spans = Vec::new();
        for pos in 0..3000i64 {
            let fill = match pos {
                1000 => 2 * BLOCK_SIZE + 10,
                1001 => BLOCK_SIZE,
                _ => (pos * 7919 % 1500) as usize,
            };
            let before = data.len();
            let line = format!("chr1\t{}\t.\tA\tC\t.\t.\tX={}\n", pos + 1, "x".repeat(fill));
            data.extend_from_slice(line.as_bytes());
            spans.push(RecordSpan {
                len: (data.len() - before) as u32,
                rid: 0,
                pos,
                rlen: 1,
            });
        }
        // HTSlib's layout: `bgzf_flush_try` ends the block before a record that does not fit,
        // and `bgzf_write` fills blocks up to `BGZF_BLOCK_SIZE`.
        let mut expected = Vec::new();
        let mut filled = 0;
        for span in &spans {
            let len = span.len as usize;
            if filled > 0 && filled + len > BLOCK_SIZE {
                expected.push(filled);
                filled = 0;
            }
            filled += len;
            while filled >= BLOCK_SIZE {
                expected.push(BLOCK_SIZE);
                filled -= BLOCK_SIZE;
            }
        }
        expected.push(filled);
        for threads in [0, 2] {
            for indexed in [false, true] {
                let sink = Shared::default();
                let options = WriteOptions {
                    output_type: OutputType::VcfBgzf,
                    level: Some(1),
                    threads,
                };
                let mut writer = Writer::from_sink(Box::new(sink.clone()), &options);
                let dir = std::env::temp_dir();
                let index = dir.join(format!("pasta-blocks-{}-{threads}.csi", std::process::id()));
                if indexed {
                    let request = IndexRequest {
                        path: index.clone(),
                        min_shift: 14,
                    };
                    writer.request_index(request).expect("BGZF output");
                }
                writer.write_header(&header).expect("write the header");
                if indexed {
                    writer.start_index(&header).expect("start the index");
                }
                let mut off = 0;
                for chunk in spans.chunks(700) {
                    let len: usize = chunk.iter().map(|s| s.len as usize).sum();
                    writer
                        .write_spans(&header, &data[off..off + len], chunk)
                        .expect("write to memory");
                    off += len;
                }
                writer.finish().expect("write to memory");
                let _ = std::fs::remove_file(&index);
                let out = sink.0.lock().expect("not poisoned");
                let sizes = block_sizes(&out);
                let header_len = sizes[0];
                assert_eq!(
                    &sizes[1..],
                    expected,
                    "threads {threads}, indexed {indexed}"
                );
                assert_eq!(header_len + data.len(), sizes.iter().sum::<usize>());
            }
        }
    }

    #[test]
    fn an_index_that_cannot_be_saved_is_reported_as_such() {
        let dir = std::env::temp_dir().join(format!("pasta-save-index-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        let out = dir.join("out.vcf.gz");
        let header = Header::parse(b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            .expect("valid header");
        let options = WriteOptions {
            output_type: OutputType::VcfBgzf,
            level: None,
            threads: 0,
        };
        let mut writer = Writer::create(Some(&out), &options).expect("create the output");
        writer
            .request_index(IndexRequest {
                path: dir.join("missing").join("out.vcf.gz.csi"),
                min_shift: 14,
            })
            .expect("BGZF output can be indexed");
        writer.write_header(&header).expect("write the header");
        writer.start_index(&header).expect("start the index");
        let result = writer.finish();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(
            matches!(result, Err(WriteError::SaveIndex(ref e)) if e.kind() == io::ErrorKind::NotFound),
            "{result:?}"
        );
    }
}
