//! The text output of `roh`: a file or standard output, plain or gzip-compressed
//! (`bgzf_open(.., "wu")` or `"wg"`).

use std::{
    cell::Cell,
    fs::File,
    io::{self, BufWriter, Write},
    rc::Rc,
};

use flate2::write::GzEncoder;

use crate::raw_input::strerror;

/// Passes writes through until closed, then drops them: output still buffered when an error
/// ends the program is lost, as with bcftools, which exits without closing its output.
struct Gate {
    inner: Box<dyn Write>,
    open: Rc<Cell<bool>>,
}

impl Write for Gate {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if self.open.get() {
            self.inner.write(buf)
        } else {
            Ok(buf.len())
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if self.open.get() {
            self.inner.flush()
        } else {
            Ok(())
        }
    }
}

enum Sink {
    Plain(BufWriter<Gate>),
    Gzip(Box<GzEncoder<BufWriter<Gate>>>),
}

/// Whether the output is compressed (`-Oz`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    Plain,
    Gzip,
}

/// The output of `roh`.
pub struct OutFile {
    sink: Option<Sink>,
    open: Rc<Cell<bool>>,
}

impl OutFile {
    /// Opens the output: standard output for `stdout` (as bcftools names it) or `-`.
    ///
    /// # Errors
    ///
    /// Returns bcftools' error when the file cannot be created.
    pub fn open(name: &str, compression: Compression) -> anyhow::Result<Self> {
        let raw: Box<dyn Write> = if name == "stdout" || name == "-" {
            Box::new(crate::output::stdout())
        } else {
            Box::new(
                File::create(name)
                    .map_err(|e| anyhow::anyhow!("Failed to open {name}: {}", strerror(&e)))?,
            )
        };
        let open = Rc::new(Cell::new(true));
        let gate = Gate {
            inner: raw,
            open: Rc::clone(&open),
        };
        let buf = BufWriter::with_capacity(1 << 17, gate);
        let sink = if compression == Compression::Gzip {
            Sink::Gzip(Box::new(GzEncoder::new(
                buf,
                flate2::Compression::default(),
            )))
        } else {
            Sink::Plain(buf)
        };
        Ok(Self {
            sink: Some(sink),
            open,
        })
    }

    /// Writes everything out and closes the output.
    ///
    /// # Errors
    ///
    /// Returns write errors.
    pub fn close(mut self) -> io::Result<()> {
        match self.sink.take() {
            Some(Sink::Plain(mut w)) => w.flush(),
            Some(Sink::Gzip(w)) => w.finish()?.flush(),
            None => Ok(()),
        }
    }
}

impl Write for OutFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.sink {
            Some(Sink::Plain(w)) => w.write(buf),
            Some(Sink::Gzip(w)) => w.write(buf),
            None => Ok(buf.len()),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match &mut self.sink {
            Some(Sink::Plain(w)) => w.write_all(buf),
            Some(Sink::Gzip(w)) => w.write_all(buf),
            None => Ok(()),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.sink {
            Some(Sink::Plain(w)) => w.flush(),
            Some(Sink::Gzip(w)) => w.flush(),
            None => Ok(()),
        }
    }
}

impl Drop for OutFile {
    /// Unclosed (on an error), buffered output is dropped.
    fn drop(&mut self) {
        self.open.set(false);
    }
}
