//! Raw access to VCF/BCF inputs, for commands that copy BGZF blocks without recompressing them
//! (bcftools `reheader` and `concat --naive`): the file or standard input, its detected
//! format, and readers of its BGZF frames or decompressed bytes.

use std::{
    collections::VecDeque,
    fs::File,
    io::{self, BufRead, BufReader, Read},
    num::NonZero,
    path::Path,
};

use pasta::{
    bgzf::{self, io::raw},
    util::threads,
    variant::{
        Header, Record, bcf,
        io::{
            Writer,
            reader::{is_unknown_binary, open_again, report_missing_eof},
        },
        log,
    },
};

use crate::output;

/// The text of an error as C's `strerror` prints it: without the ` (os error N)` that Rust's
/// I/O errors append.
pub fn strerror(e: &impl std::fmt::Display) -> String {
    let s = e.to_string();
    match s.find(" (os error") {
        Some(i) => s[..i].to_owned(),
        None => s,
    }
}

/// Prints HTSlib's messages for a file `hts_open` cannot open, `[E::hts_hopen]` first for a
/// directory (which opens but cannot be read).
pub fn log_open_error(path: &str, e: &io::Error) {
    if e.kind() == io::ErrorKind::IsADirectory {
        log::error("hts_hopen", format_args!("Failed to open file {path}"));
    }
    log::error(
        "hts_open_format",
        format_args!("Failed to open file \"{path}\" : {}", strerror(e)),
    );
}

/// The compression of an input (HTSlib's `htsFormat.compression`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    /// Plain text or uncompressed BCF.
    None,
    /// gzip that is not BGZF.
    Gzip,
    /// BGZF.
    Bgzf,
}

/// The data format of an input (HTSlib's `htsFormat.format`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// VCF text.
    Vcf,
    /// BCF.
    Bcf,
    /// Neither (HTSlib opens it as text or binary data and fails later).
    Other,
}

/// The bytes of a file or of standard input, with the bytes read for detection replayed.
pub struct Source {
    prefix: Vec<u8>,
    pos: usize,
    inner: Box<dyn Read + Send>,
}

impl Read for Source {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos < self.prefix.len() {
            let n = (self.prefix.len() - self.pos).min(buf.len());
            buf[..n].copy_from_slice(&self.prefix[self.pos..self.pos + n]);
            self.pos += n;
            return Ok(n);
        }
        self.inner.read(buf)
    }
}

/// An opened input and its detected format.
pub struct RawInput {
    pub source: BufReader<Source>,
    pub compression: Compression,
    pub format: Format,
    /// The first decompressed bytes, as HTSlib peeks them to detect the format.
    peek: Vec<u8>,
}

/// How many raw bytes are read to detect the format: one maximal BGZF frame.
const PEEK: usize = 1 << 16;

/// How many decompressed bytes HTSlib's `hts_detect_format2` looks at.
const HTS_PEEK: usize = 1024;

/// Reports why HTSlib cannot open a file of an unknown binary format, as `hts_hopen` and
/// `hts_open_format` do.
pub fn report_exec_format_error(path: &str) {
    pasta::variant::log::error("hts_hopen", format_args!("Failed to open file {path}"));
    pasta::variant::log::error(
        "hts_open_format",
        format_args!("Failed to open file \"{path}\" : {EXEC_FORMAT_ERROR}"),
    );
}

/// C's `strerror(ENOEXEC)`.
pub const EXEC_FORMAT_ERROR: &str = "Exec format error";

impl RawInput {
    /// Opens `path` (`-` for standard input) and detects its format as HTSlib's
    /// `hts_detect_format` does for VCF and BCF.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be opened or read.
    pub fn open(path: &str) -> io::Result<Self> {
        let inner: Box<dyn Read + Send> = if path == "-" {
            Box::new(io::stdin())
        } else {
            Box::new(File::open(path)?)
        };
        Self::detect(inner)
    }

    /// [`Self::open`] for an input that was read once already (pasta's [`open_again`]): a
    /// stream is read on from where it is, without waiting for the writer of a FIFO.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be opened or read.
    pub fn open_again(path: &str) -> io::Result<Self> {
        let inner: Box<dyn Read + Send> = if path == "-" {
            Box::new(io::stdin())
        } else {
            Box::new(open_again(Path::new(path))?)
        };
        Self::detect(inner)
    }

    /// Reads the first bytes of `inner` and detects its format.
    fn detect(mut inner: Box<dyn Read + Send>) -> io::Result<Self> {
        let mut prefix = vec![0u8; PEEK];
        let mut n = 0;
        while n < PEEK {
            match inner.read(&mut prefix[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        prefix.truncate(n);
        let is_gz = prefix.starts_with(&[0x1f, 0x8b]);
        let compression = if is_gz {
            if prefix.len() >= 18 && prefix[3] & 4 != 0 && &prefix[12..16] == b"BC\x02\x00" {
                Compression::Bgzf
            } else {
                Compression::Gzip
            }
        } else {
            Compression::None
        };
        let mut head = Vec::new();
        match compression {
            Compression::Bgzf => {
                let mut src = &prefix[..];
                let mut frame = Vec::new();
                if raw::read_frame(&mut src, &mut frame).unwrap_or(false) {
                    let _ = raw::inflate_frame(&frame, &mut head);
                }
            }
            Compression::None => head.extend_from_slice(&prefix[..prefix.len().min(1024)]),
            Compression::Gzip => {
                // As much of the decompressed start as the peeked bytes hold.
                let mut d = flate2::read::MultiGzDecoder::new(&prefix[..]);
                head.resize(1024, 0);
                let mut n = 0;
                while n < head.len() {
                    match d.read(&mut head[n..]) {
                        Ok(0) | Err(_) => break,
                        Ok(k) => n += k,
                    }
                }
                head.truncate(n);
            }
        }
        let format = if head.len() >= 4 && head.starts_with(b"BCF\x02") && head[3] <= 4 {
            Format::Bcf
        } else if head.starts_with(b"##fileformat=VCF") {
            Format::Vcf
        } else {
            Format::Other
        };
        head.truncate(HTS_PEEK);
        Ok(Self {
            source: BufReader::with_capacity(
                1 << 17,
                Source {
                    prefix,
                    pos: 0,
                    inner,
                },
            ),
            compression,
            format,
            peek: head,
        })
    }

    /// Returns whether HTSlib's `hts_hopen` refuses the input, failing with ENOEXEC ("Exec
    /// format error"): its first bytes are neither empty, nor text, nor a binary format that
    /// `hts_hopen` opens (BAM, BCF, CRAM; index files are refused too). Other inputs that are
    /// not VCF or BCF open and fail when their header is read.
    #[must_use]
    pub fn is_unknown_binary(&self) -> bool {
        is_unknown_binary(&self.peek)
    }

    /// Returns a reader of the decompressed bytes, with `threads` decompression threads for
    /// BGZF (fewer when the system cannot start that many, none when it cannot start any).
    #[must_use]
    pub fn into_decompressed(self, threads: usize) -> Box<dyn BufRead + Send> {
        match self.compression {
            Compression::Bgzf => match NonZero::new(threads) {
                Some(n) => {
                    match bgzf::io::MultithreadedReader::try_with_worker_count(n, self.source) {
                        Ok(reader) => Box::new(BufReader::with_capacity(1 << 17, reader)),
                        Err(e) => Box::new(bgzf::io::Reader::new(e.into_inner())),
                    }
                }
                None => Box::new(bgzf::io::Reader::new(self.source)),
            },
            Compression::Gzip => Box::new(BufReader::with_capacity(
                1 << 17,
                flate2::read::MultiGzDecoder::new(self.source),
            )),
            Compression::None => Box::new(self.source),
        }
    }
}

/// Whether a raw BGZF frame holds no data (its ISIZE, the last four bytes, is 0), as the
/// end-of-file marker does.
fn is_empty_frame(raw: &[u8]) -> bool {
    raw.ends_with(&[0, 0, 0, 0])
}

/// Reads BGZF frames one at a time, decompressing them on request (HTSlib's
/// `bgzf_read_block` and `bgzf_raw_read`).
pub struct FrameReader<R: Read> {
    inner: R,
    /// The input's name, for [`report_missing_eof`].
    name: String,
    /// The current frame, compressed.
    pub frame: Vec<u8>,
    /// The current frame's data, after [`Self::read_block`].
    pub data: Vec<u8>,
    /// Whether the last block [`Self::read_block`] read was empty (HTSlib's `last_block_eof`).
    last_block_empty: bool,
}

impl<R: Read> FrameReader<R> {
    /// Reads the frames of `inner`, the input named `name` (`-` for standard input).
    pub fn new(inner: R, name: &str) -> Self {
        Self {
            inner,
            name: name.to_owned(),
            frame: Vec::new(),
            data: Vec::new(),
            last_block_empty: false,
        }
    }

    /// Reads the next frame; `false` at the end of the input.
    ///
    /// # Errors
    ///
    /// Returns an error for a truncated or malformed frame.
    pub fn next_frame(&mut self) -> io::Result<bool> {
        raw::read_frame(&mut self.inner, &mut self.frame)
    }

    /// Reads and decompresses the next non-empty frame into [`Self::data`] (HTSlib's
    /// `bgzf_read_block`, which skips empty blocks such as the end-of-file markers of
    /// concatenated files); at the end of the input, `data` is empty, and the input is reported
    /// ([`report_missing_eof`]) when the last block was not empty (no end-of-file marker).
    ///
    /// # Errors
    ///
    /// Returns an error for a truncated, malformed, or corrupt frame.
    pub fn read_block(&mut self) -> io::Result<()> {
        loop {
            if !self.next_frame()? {
                if !self.last_block_empty {
                    report_missing_eof(&self.name);
                    // Once, as HTSlib's `no_eof_block` flag does.
                    self.last_block_empty = true;
                }
                self.data.clear();
                return Ok(());
            }
            raw::inflate_frame(&self.frame, &mut self.data)?;
            self.last_block_empty = self.data.is_empty();
            if !self.data.is_empty() {
                return Ok(());
            }
        }
    }

    /// Reads the next frame into [`Self::frame`] to copy as it is (HTSlib's `bgzf_raw_read`,
    /// but whole frames, so a truncated frame is an error rather than copied); `false` at the
    /// end of the input, which is reported ([`report_missing_eof`]) when the last frame was not
    /// empty.
    ///
    /// # Errors
    ///
    /// Returns an error for a truncated or malformed frame.
    pub fn next_raw_frame(&mut self) -> io::Result<bool> {
        if self.next_frame()? {
            self.last_block_empty = is_empty_frame(&self.frame);
            return Ok(true);
        }
        if !self.last_block_empty {
            report_missing_eof(&self.name);
            self.last_block_empty = true;
        }
        Ok(false)
    }

    /// Whether the last frame read was empty, as the end-of-file marker is, or the end of the
    /// input was reached and reported.
    pub fn last_block_empty(&self) -> bool {
        self.last_block_empty
    }

    /// Returns the underlying reader, positioned after the last frame read.
    pub fn inner_mut(&mut self) -> &mut R {
        &mut self.inner
    }
}

/// A BGZF frame and its decompressed data.
#[derive(Default)]
pub struct Frame {
    pub raw: Vec<u8>,
    pub data: Vec<u8>,
}

/// Reads BGZF frames in order and decompresses them, in batches on `threads` worker threads.
pub struct InflatingReader<R: Read> {
    inner: R,
    /// The input's name, for [`report_missing_eof`].
    name: String,
    threads: usize,
    ready: std::collections::VecDeque<Frame>,
    spare: Vec<Frame>,
    eof: bool,
    /// An error to report once the frames before it are taken.
    pending_error: Option<io::Error>,
    /// Whether the last frame read was empty, as the end-of-file marker is.
    last_frame_empty: bool,
    /// Whether the input ended without an end-of-file marker, to warn about once the frames
    /// before the end are taken.
    missing_eof: bool,
}

/// Frames per worker thread in one batch.
const FRAMES_PER_THREAD: usize = 16;

impl<R: Read> InflatingReader<R> {
    /// Reads the frames of `inner`, the input named `name` (`-` for standard input).
    pub fn new(inner: R, name: &str, threads: usize) -> Self {
        Self {
            inner,
            name: name.to_owned(),
            threads,
            ready: std::collections::VecDeque::new(),
            spare: Vec::new(),
            eof: false,
            pending_error: None,
            last_frame_empty: false,
            missing_eof: false,
        }
    }

    /// Returns the next frame, or `None` at the end of the input, which is reported
    /// ([`report_missing_eof`]) when it has no end-of-file marker.
    ///
    /// # Errors
    ///
    /// Returns an error for a truncated, malformed, or corrupt frame (after the frames before
    /// it).
    pub fn next_frame(&mut self) -> io::Result<Option<Frame>> {
        if self.ready.is_empty() && !self.eof {
            self.fill();
        }
        if let Some(f) = self.ready.pop_front() {
            return Ok(Some(f));
        }
        if let Some(e) = self.pending_error.take() {
            return Err(e);
        }
        if std::mem::take(&mut self.missing_eof) {
            report_missing_eof(&self.name);
        }
        Ok(None)
    }

    /// Returns a frame's buffers for reuse.
    pub fn recycle(&mut self, frame: Frame) {
        self.spare.push(frame);
    }

    fn fill(&mut self) {
        let n = (self.threads.max(1)) * FRAMES_PER_THREAD;
        let mut batch: Vec<Frame> = Vec::with_capacity(n);
        let mut error = None;
        while batch.len() < n {
            let mut f = self.spare.pop().unwrap_or_default();
            match raw::read_frame(&mut self.inner, &mut f.raw) {
                Ok(true) => {
                    self.last_frame_empty = is_empty_frame(&f.raw);
                    batch.push(f);
                }
                Ok(false) => {
                    self.eof = true;
                    self.missing_eof = !self.last_frame_empty;
                    break;
                }
                Err(e) => {
                    self.eof = true;
                    error = Some(e);
                    break;
                }
            }
        }
        let results: Vec<io::Result<()>> = if self.threads == 0 || batch.len() < 2 {
            batch
                .iter_mut()
                .map(|f| raw::inflate_frame(&f.raw, &mut f.data))
                .collect()
        } else {
            let chunk = batch.len().div_ceil(self.threads);
            std::thread::scope(|s| {
                // A chunk whose thread cannot be started is inflated here.
                let handles: Vec<_> = batch
                    .chunks_mut(chunk)
                    .map(|c| {
                        threads::spawn_scoped_or_run(s, move || {
                            c.iter_mut()
                                .map(|f| raw::inflate_frame(&f.raw, &mut f.data))
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().expect("inflating threads do not panic"))
                    .collect()
            })
        };
        for (f, r) in batch.into_iter().zip(results) {
            if let Err(e) = r {
                self.eof = true;
                error = Some(e);
                break;
            }
            self.ready.push_back(f);
        }
        self.pending_error = error;
    }
}

/// Bytes of emitted input dropped from the window at a time.
const TRIM_BYTES: usize = 1 << 20;

/// Copies a BGZF input to a [`Writer`], passing whole frames through compressed and
/// recompressing only around the parts the caller replaces.
///
/// The caller walks the decompressed input (records or lines) by offset, reports the parts it
/// changes with [`Self::replace`], and calls [`Self::advance`] as it goes; frames holding only
/// unchanged data are written as they are.
pub struct FrameCopier<R: Read> {
    frames: InflatingReader<R>,
    /// Frames not yet fully written, with the offset of their first decompressed byte.
    queue: VecDeque<(Frame, u64)>,
    /// Decompressed bytes from `win_start`.
    win: Vec<u8>,
    win_start: u64,
    /// Everything before this offset has been written (or skipped).
    emitted: u64,
    eof: bool,
    /// Whether reading stopped at a truncated, malformed, or corrupt frame (HTSlib's
    /// `errcode`, which makes closing the input fail), or [`Self::copy_bcf_records`] stopped at
    /// a truncated or malformed record.
    read_failed: bool,
}

impl<R: Read> FrameCopier<R> {
    pub fn new(frames: InflatingReader<R>) -> Self {
        Self {
            frames,
            queue: VecDeque::new(),
            win: Vec::new(),
            win_start: 0,
            emitted: 0,
            eof: false,
            read_failed: false,
        }
    }

    /// Whether reading stopped at a truncated, malformed, or corrupt frame, or at a truncated
    /// or malformed BCF record, rather than at the end of the input.
    pub fn read_failed(&self) -> bool {
        self.read_failed
    }

    /// The offset just past the data read so far.
    pub fn end(&self) -> u64 {
        self.win_start + self.win.len() as u64
    }

    /// Reads frames until the data reaches `upto`; `false` when the input ends (or fails)
    /// first.
    pub fn fill_to(&mut self, upto: u64) -> bool {
        while self.end() < upto {
            if self.eof {
                return false;
            }
            match self.frames.next_frame() {
                Ok(Some(f)) => {
                    let start = self.end();
                    self.win.extend_from_slice(&f.data);
                    self.queue.push_back((f, start));
                }
                Ok(None) => self.eof = true,
                Err(_) => {
                    self.eof = true;
                    self.read_failed = true;
                }
            }
        }
        true
    }

    /// Returns the decompressed bytes in `from..to`, which must have been read and not yet
    /// trimmed.
    pub fn slice(&self, from: u64, to: u64) -> &[u8] {
        &self.win[(from - self.win_start) as usize..(to - self.win_start) as usize]
    }

    /// Returns the end of the line starting at `p` (just past its newline, or the end of the
    /// input for a last line without one), or `None` at the end of the input.
    pub fn line_end(&mut self, p: u64) -> Option<u64> {
        let mut from = p;
        loop {
            if from < self.end()
                && let Some(i) =
                    memchr::memchr(b'\n', &self.win[(from - self.win_start) as usize..])
            {
                return Some(from + i as u64 + 1);
            }
            from = self.end();
            if !self.fill_to(from + 1) {
                return (self.end() > p).then(|| self.end());
            }
        }
    }

    /// Skips the data before `offset` (a header that is written separately).
    pub fn skip_to(&mut self, offset: u64) {
        self.emitted = offset;
        while let Some((f, fs)) = self.queue.front() {
            if fs + f.data.len() as u64 <= offset {
                self.pop_frame();
            } else {
                break;
            }
        }
    }

    /// Writes the input as it is up to `upto`: whole frames raw, partial ones recompressed.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn emit_original(&mut self, writer: &mut Writer, upto: u64) -> io::Result<()> {
        while self.emitted < upto {
            let Some((f, fs)) = self.queue.front() else {
                break;
            };
            let fs = *fs;
            let fe = fs + f.data.len() as u64;
            if fe <= self.emitted {
                self.pop_frame();
                continue;
            }
            if self.emitted == fs && fe <= upto {
                writer.write_bgzf_frames(&f.raw)?;
                self.emitted = fe;
                self.pop_frame();
            } else {
                let to = fe.min(upto);
                writer.write_raw(self.slice(self.emitted, to))?;
                self.emitted = to;
                if to == fe {
                    self.pop_frame();
                }
            }
        }
        Ok(())
    }

    /// Writes `bytes` in place of the input in `from..to`.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn replace(
        &mut self,
        writer: &mut Writer,
        from: u64,
        to: u64,
        bytes: &[u8],
    ) -> io::Result<()> {
        self.emit_original(writer, from)?;
        writer.write_raw(bytes)?;
        self.skip_to(to);
        Ok(())
    }

    /// Notes that everything before `p` has been checked: whole frames up to there are
    /// written, and data no longer needed is dropped.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn advance(&mut self, writer: &mut Writer, p: u64) -> io::Result<()> {
        let mut upto = self.emitted;
        for (f, fs) in &self.queue {
            let fe = fs + f.data.len() as u64;
            if fe > p {
                break;
            }
            upto = upto.max(fe);
        }
        self.emit_original(writer, upto)?;
        let keep = self.emitted.min(p);
        let drop = (keep - self.win_start) as usize;
        if drop >= TRIM_BYTES {
            self.win.drain(..drop);
            self.win_start = keep;
        }
        Ok(())
    }

    fn pop_frame(&mut self) {
        if let Some((f, _)) = self.queue.pop_front() {
            self.frames.recycle(f);
        }
    }

    /// Reads a BCF header at the start of the input; returns it and the offset of the first
    /// record.
    pub fn read_bcf_header(&mut self) -> Option<(Header, u64)> {
        let mut need = 9u64;
        if self.fill_to(need) {
            need += u64::from(u32::from_le_bytes(
                self.win[5..9].try_into().expect("4 bytes"),
            ));
            self.fill_to(need);
        }
        let hdr = bcf::read_header(&mut &self.win[..]).ok()?;
        Some((hdr, need))
    }

    /// Copies the BCF records from `start` as HTSlib reads and writes them: each record is
    /// read and checked against `hdr` (`bcf_read`), handed to `on_record`, and written as it
    /// was unless reading or `on_record` changed it, in which case it is encoded for
    /// `out_hdr`. Reading stops at a truncated or malformed record, after writing the records
    /// before it, as bcftools' loops over `bcf_read` do; unlike them (which take such a record
    /// for the end of the input), this marks the reading failed (see [`Self::read_failed`]).
    /// `write_err` makes the error for a record that cannot be encoded.
    ///
    /// # Errors
    ///
    /// Returns the errors of `on_record`, encoding, and writing.
    pub fn copy_bcf_records(
        &mut self,
        start: u64,
        hdr: &Header,
        out_hdr: &Header,
        writer: &mut Writer,
        write_err: &dyn Fn() -> anyhow::Error,
        mut on_record: impl FnMut(&mut Record) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let mut rec = Record::new();
        let mut scratch = Vec::new();
        let mut p = start;
        self.skip_to(start);
        loop {
            if !self.fill_to(p + 32) {
                // Bytes short of a record's fixed fields are a truncated record.
                self.read_failed |= self.end() > p;
                break;
            }
            let fixed = self.slice(p, p + 32);
            let le32 = |i: usize| u32::from_le_bytes(fixed[i..i + 4].try_into().expect("4 bytes"));
            let l_shared = u64::from(le32(0));
            let l_indiv = u64::from(le32(4));
            let rid = le32(8) as i32;
            let rlen = le32(16) as i32;
            let n_fmt = fixed[31];
            let len = 8 + l_shared + l_indiv;
            if !self.fill_to(p + len) {
                self.read_failed = true;
                break;
            }
            let mut src = self.slice(p, p + len);
            if !matches!(
                bcf::read_record_raw(&mut src, &mut rec),
                Ok(bcf::ReadOutcome::Record)
            ) || bcf::check_record(Some(hdr), &mut rec).is_err()
            {
                self.read_failed = true;
                break;
            }
            on_record(&mut rec)?;
            let bytes = self.slice(p, p + len);
            let shared_end = (8 + l_shared) as usize;
            let unchanged = rec.rid() == rid
                && rec.rlen() == i64::from(rlen)
                && rec.format_count() == usize::from(n_fmt)
                && rec.shared_bytes() == &bytes[32..shared_end]
                && rec.indiv_bytes() == &bytes[shared_end..];
            if !unchanged {
                scratch.clear();
                bcf::encode_record(out_hdr, &mut rec, &mut scratch).map_err(|_| write_err())?;
                self.replace(writer, p, p + len, &scratch)
                    .map_err(|e| output::write_failure(e, write_err))?;
            }
            p += len;
            self.advance(writer, p)
                .map_err(|e| output::write_failure(e, write_err))?;
        }
        self.emit_original(writer, p)
            .map_err(|e| output::write_failure(e, write_err))?;
        Ok(())
    }
}

/// Creates a FIFO for tests in a new directory, named after `name` and the process, under the
/// temporary directory. Returns the directory, which the caller removes, and the FIFO.
#[cfg(test)]
pub fn make_fifo(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
    use std::os::unix::ffi::OsStrExt;

    let dir = std::env::temp_dir().join(format!("vtools-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("create the test directory");
    let path = dir.join("fifo");
    let c = std::ffi::CString::new(path.as_os_str().as_bytes()).expect("no NUL in the path");
    // SAFETY: `c` is a NUL-terminated path, as mkfifo requires.
    let rc = unsafe { libc::mkfifo(c.as_ptr(), 0o600) };
    assert_eq!(rc, 0, "mkfifo: {}", io::Error::last_os_error());
    (dir, path)
}

/// Writes `text` as `file` for tests in a new directory, named after `name` and the process,
/// under the temporary directory. Returns the directory, which the caller removes, and the
/// file's name as `FILE##idx##INDEX`, with an index that does not exist.
#[cfg(test)]
pub fn write_idx_named(name: &str, file: &str, text: &str) -> (std::path::PathBuf, String) {
    let dir = std::env::temp_dir().join(format!("vtools-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("create the test directory");
    let path = dir.join(file);
    std::fs::write(&path, text).expect("write the test file");
    let named = format!("{}##idx##{}", path.display(), dir.join("none").display());
    (dir, named)
}

/// Runs `f` on another thread for tests, failing if it does not return within ten seconds (it
/// would wait forever, for the writer of a FIFO say). `what` names what `f` does, for the
/// failure message.
#[cfg(test)]
pub fn within_seconds<T: Send + 'static>(what: &str, f: impl FnOnce() -> T + Send + 'static) -> T {
    use std::sync::mpsc::RecvTimeoutError;
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || tx.send(f()));
    match rx.recv_timeout(std::time::Duration::from_secs(10)) {
        Ok(value) => value,
        Err(RecvTimeoutError::Timeout) => panic!("{what} hangs: no return within ten seconds"),
        Err(RecvTimeoutError::Disconnected) => panic!("{what} panicked"),
    }
}

/// A directory for test files under the temporary directory, named after `name` and the
/// process, and removed when dropped, whether the test passes or fails.
#[cfg(test)]
pub struct TestDir(std::path::PathBuf);

#[cfg(test)]
impl TestDir {
    /// Creates the directory.
    pub fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("vtools-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        Self(dir)
    }
}

#[cfg(test)]
impl std::ops::Deref for TestDir {
    type Target = std::path::Path;

    fn deref(&self) -> &std::path::Path {
        &self.0
    }
}

#[cfg(test)]
impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn bgzf(data: &[u8]) -> Vec<u8> {
        let mut w = bgzf::io::Writer::new(Vec::new());
        w.write_all(data).expect("writing to memory works");
        w.finish().expect("finishing in memory works")
    }

    const BCF_HDR: &[u8] = b"##fileformat=VCFv4.2\n##contig=<ID=a>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    /// Returns three records as uncompressed BCF (POS 1, 2, 3), and where each record starts.
    fn raw_bcf() -> (Vec<u8>, Vec<usize>) {
        let mut hdr = Header::parse(BCF_HDR).expect("the header parses");
        let mut out = bcf::encode_header(&hdr);
        let mut starts = Vec::new();
        let mut rec = Record::new();
        for pos in 1..=3 {
            starts.push(out.len());
            let mut line = format!("a\t{pos}\t.\tA\tG\t.\t.\tDP={pos}").into_bytes();
            pasta::variant::vcf::parse(&mut hdr, &mut line, &mut rec).expect("the line parses");
            bcf::encode_record(&hdr, &mut rec, &mut out).expect("the record encodes");
        }
        (out, starts)
    }

    /// Copies the records of bgzipped `data` with [`FrameCopier::copy_bcf_records`]; returns
    /// whether reading failed and the positions written.
    fn copy_records(data: &[u8]) -> (bool, Vec<i64>) {
        use pasta::variant::io::{OutputType, WriteOptions};

        let dir = TestDir::new(&format!("copy-records-{}", data.len()));
        let path = dir.join("out.bcf");
        let compressed = bgzf(data);
        let mut copier = FrameCopier::new(InflatingReader::new(&compressed[..], "test", 0));
        let (hdr, start) = copier.read_bcf_header().expect("the header reads");
        let mut writer = Writer::create(
            Some(&path),
            &WriteOptions {
                output_type: OutputType::Bcf,
                level: None,
                threads: 0,
            },
        )
        .expect("create the output");
        writer.write_header(&hdr).expect("write the header");
        copier
            .copy_bcf_records(
                start,
                &hdr,
                &hdr,
                &mut writer,
                &|| anyhow::anyhow!("write"),
                |_| Ok(()),
            )
            .expect("the records copy");
        writer.finish().expect("finish the output");
        let file = File::open(&path).expect("open the output");
        let mut src = BufReader::new(bgzf::io::Reader::new(file));
        let hdr = bcf::read_header(&mut src).expect("the output header reads");
        let mut rec = Record::new();
        let mut pos = Vec::new();
        while let bcf::ReadOutcome::Record =
            bcf::read_record(&mut src, &hdr, &mut rec).expect("the output records read")
        {
            pos.push(rec.pos() + 1);
        }
        (copier.read_failed(), pos)
    }

    #[test]
    fn copying_bcf_records_fails_at_a_truncated_or_malformed_record() {
        let (data, starts) = raw_bcf();
        assert_eq!(copy_records(&data), (false, vec![1, 2, 3]));
        // Cut inside the last record's fixed fields, and after them.
        assert_eq!(copy_records(&data[..starts[2] + 20]), (true, vec![1, 2]));
        assert_eq!(copy_records(&data[..data.len() - 3]), (true, vec![1, 2]));
        let with = |at: usize, bytes: [u8; 4]| {
            let mut d = data.clone();
            d[at..at + 4].copy_from_slice(&bytes);
            d
        };
        // The second record's shared block runs past the data, or is shorter than its fixed
        // fields; its contig is not in the header.
        let s = starts[1];
        assert_eq!(
            copy_records(&with(s, 1000u32.to_le_bytes())),
            (true, vec![1])
        );
        assert_eq!(copy_records(&with(s, 10u32.to_le_bytes())), (true, vec![1]));
        assert_eq!(
            copy_records(&with(s + 8, 9i32.to_le_bytes())),
            (true, vec![1])
        );
    }

    #[test]
    fn frame_reader_skips_empty_blocks() {
        // Two bgzipped files joined with cat: an end-of-file block sits between the data.
        let mut joined = bgzf(b"abc");
        joined.extend(bgzf(b"def"));
        let mut fr = FrameReader::new(&joined[..], "test");
        fr.read_block().expect("the first block inflates");
        assert_eq!(fr.data, b"abc");
        fr.read_block().expect("the second block inflates");
        assert_eq!(fr.data, b"def");
        fr.read_block().expect("the end is not an error");
        assert!(fr.data.is_empty());
    }
}
