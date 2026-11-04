#![doc(hidden)]

pub mod index;
mod num;

use std::io::{self, Write};

use crate::bgzf;

use self::index::write_index;
use crate::csi::Index;

/// A CSI writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::io::Writer<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a CSI writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::csi;
    /// let writer = csi::io::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::io::Writer::new(writer),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::csi;
    /// let writer = csi::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::io::Writer<W> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::csi;
    /// let mut writer = csi::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::io::Writer<W> {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::csi;
    /// let writer = csi::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::io::Writer<W> {
        self.inner
    }

    /// Attempts to finish the output stream.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::csi;
    /// let mut writer = csi::io::Writer::new(Vec::new());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.inner.try_finish()
    }

    /// Writes a coordinate-sorted index (CSI).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::csi;
    /// let index = csi::Index::default();
    /// let mut writer = csi::io::Writer::new(Vec::new());
    /// writer.write_index(&index)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index)
    }
}
