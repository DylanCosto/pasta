//! Indexed BGZF reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead, Read, Seek, SeekFrom};

use super::Reader;
use crate::bgzf::{VirtualPosition, gzi};

/// An indexed BGZF reader.
pub struct IndexedReader<R> {
    inner: Reader<R>,
    index: gzi::Index,
    /// The uncompressed position, unless a seek failed.
    position: Option<u64>,
}

impl<R> IndexedReader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let reader = bgzf::io::IndexedReader::new(io::empty(), gzi::Index::default());
    /// let inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let mut reader = bgzf::io::IndexedReader::new(io::empty(), gzi::Index::default());
    /// let inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    /// Unwraps and returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let reader = bgzf::io::IndexedReader::new(io::empty(), gzi::Index::default());
    /// let inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Creates an indexed BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let reader = bgzf::io::IndexedReader::new(io::empty(), gzi::Index::default());
    /// ```
    pub fn new(inner: R, index: gzi::Index) -> Self {
        Self {
            inner: Reader::new(inner),
            index,
            position: Some(0),
        }
    }

    /// Returns the current position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let reader = bgzf::io::IndexedReader::new(io::empty(), gzi::Index::default());
    /// assert_eq!(reader.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.inner.position()
    }

    /// Returns the current virtual position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let reader = bgzf::io::IndexedReader::new(io::empty(), gzi::Index::default());
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.inner.virtual_position()
    }

    /// Returns the associated index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf::{self as bgzf, gzi};
    /// let index = gzi::Index::default();
    /// let reader = bgzf::io::IndexedReader::new(io::empty(), index.clone());
    /// assert_eq!(reader.index(), &index);
    /// ```
    pub fn index(&self) -> &gzi::Index {
        &self.index
    }
}

impl<R> Read for IndexedReader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.advance(n);
        Ok(n)
    }
}

impl<R> BufRead for IndexedReader<R>
where
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt);
        self.advance(amt);
    }
}

impl<R> IndexedReader<R> {
    fn advance(&mut self, n: usize) {
        if let Some(position) = &mut self.position {
            *position += n as u64;
        }
    }
}

/// Seeks by uncompressed position.
///
/// [`SeekFrom::End`] is unsupported ([`io::ErrorKind::Unsupported`]), as the index does not give
/// the stream's uncompressed size.
impl<R> Seek for IndexedReader<R>
where
    R: Read + Seek,
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        use super::seek::{end_unsupported, offset_position, resolve_with_index};

        let pos = match (pos, self.position) {
            (SeekFrom::Start(pos), _) => pos,
            (SeekFrom::Current(offset), Some(position)) => offset_position(position, offset)?,
            (SeekFrom::End(_), _) => return Err(end_unsupported()),
            // After a failed seek, the position is found from the index.
            (SeekFrom::Current(_), None) => {
                resolve_with_index(&self.index, self.inner.virtual_position(), pos)?
            }
        };

        match self.inner.seek_by_uncompressed_position(&self.index, pos) {
            Ok(pos) => {
                self.position = Some(pos);
                Ok(pos)
            }
            Err(e) => {
                self.position = None;
                Err(e)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        io::{Cursor, Write},
        num::NonZero,
    };

    use super::*;
    use crate::bgzf::io::MultithreadedReader;

    /// Returns BGZF data of `data` in blocks of 1000 bytes, and its index.
    fn indexed_data(data: &[u8]) -> io::Result<(Vec<u8>, gzi::Index)> {
        let mut writer = crate::bgzf::io::Writer::new(Vec::new());
        let mut entries = Vec::new();

        for (i, chunk) in data.chunks(1000).enumerate() {
            if i > 0 {
                entries.push((writer.get_ref().len() as u64, (i * 1000) as u64));
            }
            writer.write_all(chunk)?;
            writer.flush()?;
        }

        Ok((writer.finish()?, gzi::Index::from(entries)))
    }

    #[test]
    fn test_seek_relative() -> io::Result<()> {
        let data: Vec<u8> = (0..10_000u32).map(|i| (i % 251) as u8).collect();
        let (compressed, index) = indexed_data(&data)?;
        let byte_at = |pos: u64| data[pos as usize];

        let mut reader = IndexedReader::new(Cursor::new(compressed.clone()), index.clone());
        assert_eq!(reader.stream_position()?, 0);
        reader.seek(SeekFrom::Start(7000))?;
        assert_eq!(reader.stream_position()?, 7000);
        let mut buf = [0; 1500];
        reader.read_exact(&mut buf)?;
        assert_eq!(reader.stream_position()?, 8500);
        assert_eq!(reader.seek(SeekFrom::Current(-3000))?, 5500);
        reader.read_exact(&mut buf[..1])?;
        assert_eq!(buf[0], byte_at(5500));
        reader.read_to_end(&mut Vec::new())?;
        assert_eq!(reader.stream_position()?, 10_000);
        assert_eq!(
            reader.seek(SeekFrom::Current(-10_001)).unwrap_err().kind(),
            io::ErrorKind::InvalidInput
        );
        assert_eq!(
            reader.seek(SeekFrom::End(0)).unwrap_err().kind(),
            io::ErrorKind::Unsupported
        );
        assert_eq!(reader.stream_position()?, 10_000);

        // `bgzf::io::Reader` and `MultithreadedReader`, which seek with an index they are given.
        let mut reader = crate::bgzf::io::Reader::new(Cursor::new(compressed.clone()));
        let mut mt_reader = MultithreadedReader::with_worker_count(
            NonZero::new(2).unwrap(),
            Cursor::new(compressed),
        );

        for reader in [&mut reader as &mut dyn SeekRead, &mut mt_reader] {
            assert_eq!(reader.seek_with_index(&index, SeekFrom::Start(7000))?, 7000);
            reader.read_exact(&mut buf)?;
            assert_eq!(reader.seek_with_index(&index, SeekFrom::Current(0))?, 8500);
            assert_eq!(
                reader.seek_with_index(&index, SeekFrom::Current(-3000))?,
                5500
            );
            reader.read_exact(&mut buf[..1])?;
            assert_eq!(buf[0], byte_at(5500));
            assert_eq!(
                reader
                    .seek_with_index(&index, SeekFrom::End(0))
                    .unwrap_err()
                    .kind(),
                io::ErrorKind::Unsupported
            );
        }

        Ok(())
    }

    trait SeekRead: crate::bgzf::io::Seek + Read {}

    impl<T> SeekRead for T where T: crate::bgzf::io::Seek + Read {}
}
