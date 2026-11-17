#![doc(hidden)]

pub mod bins;
mod metadata;

use std::{
    error, fmt,
    io::{self, Read},
    num,
};

use self::bins::read_bins;
pub use self::metadata::read_metadata;
use crate::csi::{
    binning_index::index::{ReferenceSequence, reference_sequence::index::BinnedIndex},
    io::reader::num::read_i32_le,
};

/// The most entries reserved up front for a count read from an index.
///
/// Counts in an index are untrusted, and a corrupt one would otherwise reserve gigabytes before
/// any entry is read. A longer list grows as its entries are read. Every bin of a BAI or tabix
/// index (37,450 at most) fits.
pub(crate) const MAX_INITIAL_CAPACITY: usize = 1 << 16;

/// Returns the number of entries to reserve up front for a count read from an index: the
/// count, up to [`MAX_INITIAL_CAPACITY`].
pub(crate) fn initial_capacity(count: usize) -> usize {
    count.min(MAX_INITIAL_CAPACITY)
}

/// Makes room in `entries` for one more entry read from an index.
///
/// The list grows as `push` grows it, but a failed allocation returns an `OutOfMemory` error
/// instead of aborting: each entry takes more memory than the bytes it is read from, so a
/// corrupt count over a long run of compressed zeros can hold more entries than fit in memory.
pub(crate) fn reserve_entry<T>(entries: &mut Vec<T>, what: &str) -> io::Result<()> {
    entries
        .try_reserve(1)
        .map_err(|_| out_of_memory_error(entries.len() + 1, what))
}

/// Makes room in `entries` for one more entry read from an index, as [`reserve_entry`] does for
/// a list.
pub(crate) fn reserve_map_entry<K, V, S>(
    entries: &mut indexmap::IndexMap<K, V, S>,
    what: &str,
) -> io::Result<()> {
    entries
        .try_reserve(1)
        .map_err(|_| out_of_memory_error(entries.len() + 1, what))
}

fn out_of_memory_error(count: usize, what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::OutOfMemory,
        format!("could not allocate {count} {what}: out of memory"),
    )
}

/// An error returned when CSI reference sequences fail to be read.
#[derive(Debug)]
pub enum ReadError {
    /// An I/O error.
    Io(io::Error),
    /// The reference sequence count is invalid.
    InvalidReferenceSequenceCount(num::TryFromIntError),
    /// The bins are invalid.
    InvalidBins(bins::ReadError),
}

impl error::Error for ReadError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidReferenceSequenceCount(e) => Some(e),
            Self::InvalidBins(e) => Some(e),
        }
    }
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidReferenceSequenceCount(_) => write!(f, "invalid reference sequence count"),
            Self::InvalidBins(_) => write!(f, "invalid bins"),
        }
    }
}

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
    depth: u8,
) -> Result<Vec<ReferenceSequence<BinnedIndex>>, ReadError>
where
    R: Read,
{
    let reference_sequence_count = read_reference_sequence_count(reader)?;

    let mut reference_sequences = Vec::with_capacity(initial_capacity(reference_sequence_count));

    for _ in 0..reference_sequence_count {
        let reference_sequence = read_reference_sequence(reader, depth)?;
        reserve_entry(&mut reference_sequences, "reference sequences").map_err(ReadError::Io)?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

fn read_reference_sequence_count<R>(reader: &mut R) -> Result<usize, ReadError>
where
    R: Read,
{
    read_i32_le(reader)
        .map_err(ReadError::Io)
        .and_then(|n| usize::try_from(n).map_err(ReadError::InvalidReferenceSequenceCount))
}

fn read_reference_sequence<R>(
    reader: &mut R,
    depth: u8,
) -> Result<ReferenceSequence<BinnedIndex>, ReadError>
where
    R: Read,
{
    let (bins, index, metadata) = read_bins(reader, depth).map_err(ReadError::InvalidBins)?;
    Ok(ReferenceSequence::new(bins, index, metadata))
}
