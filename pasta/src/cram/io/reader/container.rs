pub mod block;
pub mod compression_header;
pub mod header;
pub mod slice;

use std::io::{self, Read};

use self::{block::read_block_as, header::read_header};
pub use self::{
    block::{Block, read_block},
    compression_header::read_compression_header,
    slice::Slice,
    slice::read_slice,
};
use crate::cram::container::{CompressionHeader, Header};

/// A CRAM container.
#[derive(Default)]
pub struct Container {
    pub(crate) header: Header,
    pub(crate) src: Vec<u8>,
}

impl Container {
    /// Returns the header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Returns the raw container data that follows the header.
    ///
    /// This is the compression header block followed by the slices. The header landmarks are
    /// offsets of the slices into this data.
    pub fn src(&self) -> &[u8] {
        &self.src
    }

    /// Returns the compression header.
    pub fn compression_header(&self) -> io::Result<CompressionHeader> {
        let end = self
            .header
            .landmarks
            .first()
            .copied()
            .unwrap_or(self.src.len());

        let mut src = self.src.get(..end).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid compression header bounds",
            )
        })?;

        read_compression_header(&mut src)
    }

    pub(crate) fn slice_at_landmark(&self, landmark: u64) -> io::Result<Slice<'_>> {
        let landmark =
            usize::try_from(landmark).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let landmarks = self.header.landmarks();
        let i = landmarks
            .iter()
            .position(|&n| n == landmark)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid slice landmark"))?;
        self.slice(i)
    }

    /// Returns the number of slices.
    pub(crate) fn slice_count(&self) -> usize {
        self.header.landmarks.len()
    }

    /// Returns the `i`th slice.
    ///
    /// # Errors
    ///
    /// Returns an error if the landmarks put the slice outside the container.
    ///
    /// # Panics
    ///
    /// Panics if `i` is not less than [`Self::slice_count`].
    pub(crate) fn slice(&self, i: usize) -> io::Result<Slice<'_>> {
        let landmarks = &self.header.landmarks;
        let start = landmarks[i];
        let end = landmarks.get(i + 1).copied().unwrap_or(self.src.len());
        let mut src = self
            .src
            .get(start..end)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid slice bounds"))?;
        read_slice(&mut src)
    }

    /// Returns the iterator over slices.
    pub fn slices(&self) -> impl Iterator<Item = io::Result<Slice<'_>>> + '_ {
        (0..self.slice_count()).map(|i| self.slice(i))
    }
}

pub fn read_container<R>(reader: &mut R, container: &mut Container) -> io::Result<usize>
where
    R: Read,
{
    match read_header(reader, &mut container.header)? {
        0 => Ok(0),
        len => {
            container.src.resize(len, 0);
            reader.read_exact(&mut container.src)?;
            Ok(len)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_slices_with_a_landmark_past_the_end() {
        let mut container = Container {
            src: vec![0; 8],
            ..Default::default()
        };
        container.header.landmarks = vec![4, 16];

        let results: Vec<_> = container.slices().collect();
        assert_eq!(results.len(), 2);

        let e = results[1].as_ref().err().expect("an error");
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_compression_header_with_a_landmark_past_the_end() {
        let mut container = Container {
            src: vec![0; 8],
            ..Default::default()
        };
        container.header.landmarks = vec![16];

        let e = container.compression_header().expect_err("an error");
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
    }
}
