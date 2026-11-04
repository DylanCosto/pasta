//! Text inputs and outputs of `convert`: line readers for plain or gzipped files (HTSlib's
//! `hts_getline`), gzipped or plain outputs (`bgzf_open` with `wg`/`wu`), and the reference
//! FASTA (`faidx`).

use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    num::NonZero,
    path::Path,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};

use pasta::{
    bgzf,
    core::{Position, Region},
    fasta::{self, fai},
};

use crate::{
    faidx::{Faidx, LoadError},
    raw_input::strerror,
};

/// A text file read line by line, gzip-compressed or not.
pub struct LineReader {
    inner: Box<dyn BufRead>,
}

impl LineReader {
    /// Opens a file, or standard input for `-`, detecting gzip compression.
    pub fn open(path: &str) -> io::Result<Self> {
        let raw: Box<dyn Read> = if path == "-" {
            Box::new(io::stdin())
        } else {
            Box::new(File::open(path)?)
        };
        let mut buffered = BufReader::with_capacity(1 << 17, raw);
        let head = buffered.fill_buf()?;
        let inner: Box<dyn BufRead> = if head.starts_with(&[0x1f, 0x8b]) {
            Box::new(BufReader::with_capacity(
                1 << 17,
                flate2::bufread::MultiGzDecoder::new(buffered),
            ))
        } else {
            Box::new(buffered)
        };
        Ok(Self { inner })
    }

    /// Reads the next line without its newline (and a carriage return before it). Returns
    /// `false` at the end and for an empty line, as `hts_getline() <= 0` does.
    pub fn next_line(&mut self, line: &mut Vec<u8>) -> io::Result<bool> {
        line.clear();
        if self.inner.read_until(b'\n', line)? == 0 {
            return Ok(false);
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.len() > 1 && line.last() == Some(&b'\r') {
            line.pop();
        }
        Ok(!line.is_empty())
    }
}

/// How a text output is compressed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TextCompression {
    /// `bgzf_open(..., "wg")`: gzip.
    Gzip,
    /// `bgzf_open(..., "wu")`: plain text.
    Plain,
}

/// A text output file (or standard output for `-`).
pub struct TextWriter {
    inner: Option<Sink>,
    name: String,
    open: Arc<AtomicBool>,
}

enum Sink {
    Gzip(Box<flate2::write::GzEncoder<BufWriter<Gate>>>),
    /// Gzip output compressed on worker threads: BGZF, as HTSlib writes a `wg` file given a
    /// thread pool.
    Bgzf(Box<bgzf::io::MultithreadedWriter<Gate>>),
    /// The same BGZF output compressed on the calling thread, when no worker thread starts.
    BgzfSingle(Box<bgzf::io::Writer<Gate>>),
    Plain(BufWriter<Gate>),
}

/// Passes writes through until closed, then drops them: buffered output is lost on a fatal
/// error, as with bcftools, which exits without closing its BGZF files.
struct Gate {
    inner: Box<dyn Write + Send>,
    open: Arc<AtomicBool>,
}

impl Write for Gate {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if self.open.load(Ordering::Acquire) {
            self.inner.write(buf)
        } else {
            Ok(buf.len())
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.open.load(Ordering::Acquire) {
            self.inner.flush()
        } else {
            Ok(())
        }
    }
}

impl TextWriter {
    pub fn create(name: &str, compression: TextCompression) -> io::Result<Self> {
        Self::create_threaded(name, compression, 0)
    }

    /// Creates an output whose gzip compression, if any, runs on `threads` workers (none for
    /// 0), as bcftools gives its hap files the thread pool; on the calling thread when no worker
    /// starts (the output is the same).
    pub fn create_threaded(
        name: &str,
        compression: TextCompression,
        threads: usize,
    ) -> io::Result<Self> {
        let raw: Box<dyn Write + Send> = if name == "-" {
            Box::new(crate::output::stdout())
        } else {
            Box::new(File::create(name)?)
        };
        let open = Arc::new(AtomicBool::new(true));
        let gate = Gate {
            inner: raw,
            open: Arc::clone(&open),
        };
        let inner = match (compression, NonZero::new(threads)) {
            (TextCompression::Gzip, Some(workers)) => {
                match bgzf::io::MultithreadedWriter::try_with_worker_count(workers, gate) {
                    Ok(w) => Sink::Bgzf(Box::new(w)),
                    Err(e) => Sink::BgzfSingle(Box::new(bgzf::io::Writer::new(e.into_inner()))),
                }
            }
            // BGZF's block size: what bcftools has written before a fatal error.
            (TextCompression::Gzip, None) => Sink::Gzip(Box::new(flate2::write::GzEncoder::new(
                BufWriter::with_capacity(0xff00, gate),
                flate2::Compression::default(),
            ))),
            (TextCompression::Plain, _) => Sink::Plain(BufWriter::with_capacity(0xff00, gate)),
        };
        Ok(Self {
            inner: Some(inner),
            name: name.to_owned(),
            open,
        })
    }

    pub fn write(&mut self, data: &[u8]) -> anyhow::Result<()> {
        let result = match self.inner.as_mut() {
            Some(Sink::Gzip(w)) => w.write_all(data),
            Some(Sink::Bgzf(w)) => w.write_all(data),
            Some(Sink::BgzfSingle(w)) => w.write_all(data),
            Some(Sink::Plain(w)) => w.write_all(data),
            None => Ok(()),
        };
        result.map_err(|e| {
            let msg = format!("Error writing {}: {}", self.name, strerror(&e));
            anyhow::Error::new(e).context(msg)
        })
    }

    pub fn close(mut self) -> anyhow::Result<()> {
        let result = match self.inner.take() {
            Some(Sink::Gzip(w)) => w.finish().and_then(|mut raw| raw.flush()),
            Some(Sink::Bgzf(mut w)) => w.finish().and_then(|mut raw| raw.flush()),
            Some(Sink::BgzfSingle(w)) => w.finish().and_then(|mut raw| raw.flush()),
            Some(Sink::Plain(mut w)) => w.flush(),
            None => Ok(()),
        };
        result.map_err(|e| {
            let msg = format!("Error closing {}: {}", self.name, strerror(&e));
            anyhow::Error::new(e).context(msg)
        })
    }
}

impl Drop for TextWriter {
    /// Unclosed (on an error), buffered output is dropped.
    fn drop(&mut self) {
        self.open.store(false, Ordering::Release);
    }
}

/// Whether an output name ends with `.gz` (case-insensitive).
pub fn has_gz_suffix(name: &str) -> bool {
    name.len() >= 3 && name[name.len() - 3..].eq_ignore_ascii_case(".gz")
}

/// Bases read from the reference at a time.
const REF_WINDOW: u64 = 1 << 20;

/// An indexed FASTA reference (HTSlib's `faidx_t`), read in windows.
pub struct Reference {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    /// The cached window: sequence name, 0-based start, and bases.
    cached: Option<(Vec<u8>, u64, Vec<u8>)>,
}

impl Reference {
    /// Opens a FASTA file, plain or bgzipped, building its index if missing (HTSlib's
    /// `fai_load`).
    ///
    /// # Errors
    ///
    /// Returns an error, after HTSlib's messages, if the file or its index cannot be read,
    /// built, or saved.
    pub fn open(path: &Path) -> Result<Self, LoadError> {
        Ok(Self {
            reader: Faidx::load(path)?.into_indexed_reader(),
            cached: None,
        })
    }

    /// The sequence names and lengths, in index order (`faidx_iseq`, `faidx_seq_len`).
    pub fn sequences(&self) -> Vec<(Vec<u8>, u64)> {
        self.reader
            .index()
            .as_ref()
            .iter()
            .map(|r| (r.name().to_vec(), r.length()))
            .collect()
    }

    /// Fetches the base at 0-based `pos` into `out` (HTSlib's `faidx_fetch_seq` of one base):
    /// empty beyond the end of the sequence. Returns `false` for an unknown sequence.
    pub fn fetch(&mut self, name: &[u8], pos: i64, out: &mut Vec<u8>) -> anyhow::Result<bool> {
        out.clear();
        let pos = pos.max(0) as u64;
        if let Some((n, start, bases)) = &self.cached
            && n == name
            && pos >= *start
            && pos < start + bases.len() as u64
        {
            out.push(bases[(pos - start) as usize]);
            return Ok(true);
        }
        let Some(length) = self
            .reader
            .index()
            .as_ref()
            .iter()
            .find(|r| r.name() == name)
            .map(fai::Record::length)
        else {
            return Ok(false);
        };
        if pos >= length {
            return Ok(true);
        }
        let end = (pos + REF_WINDOW).min(length);
        let start_pos = Position::try_from(pos as usize + 1)?;
        let end_pos = Position::try_from(end as usize)?;
        let record = self
            .reader
            .query(&Region::new(name.to_vec(), start_pos..=end_pos))?;
        let bases = record.sequence().as_ref().to_vec();
        if let Some(&b) = bases.first() {
            out.push(b);
        }
        self.cached = Some((name.to_vec(), pos, bases));
        Ok(true)
    }
}
