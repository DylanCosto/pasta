use std::io::{self, SeekFrom};

use crate::bgzf::{VirtualPosition, gzi};

/// A seekable BGZF reader.
pub trait Seek {
    /// Seeks the stream to the given virtual position.
    fn seek_to_virtual_position(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition>;

    /// Seeks the stream to the given position using an index.
    fn seek_with_index(&mut self, index: &gzi::Index, pos: SeekFrom) -> io::Result<u64>;
}

/// Returns the uncompressed position that `pos` seeks to, from the stream's virtual position and
/// its index.
///
/// # Errors
///
/// [`SeekFrom::End`] is unsupported, as the index does not give the stream's uncompressed size.
/// [`SeekFrom::Current`] fails with [`io::ErrorKind::InvalidInput`] when it seeks before the
/// start, or when the index does not list the current block (e.g., at the end of the stream).
pub(crate) fn resolve_with_index(
    index: &gzi::Index,
    virtual_position: VirtualPosition,
    pos: SeekFrom,
) -> io::Result<u64> {
    match pos {
        SeekFrom::Start(pos) => Ok(pos),
        SeekFrom::Current(offset) => {
            let current = uncompressed_position(index, virtual_position)?;
            offset_position(current, offset)
        }
        SeekFrom::End(_) => Err(end_unsupported()),
    }
}

/// Returns the uncompressed position at a virtual position, which is at the start of a block
/// the index lists.
fn uncompressed_position(index: &gzi::Index, virtual_position: VirtualPosition) -> io::Result<u64> {
    let (compressed_pos, block_data_pos) = virtual_position.into();

    let block_start = if compressed_pos == 0 {
        Some(0)
    } else {
        let entries = index.as_ref();
        entries
            .binary_search_by_key(&compressed_pos, |&(c, _)| c)
            .ok()
            .map(|i| entries[i].1)
    };

    block_start
        .map(|start| start + u64::from(block_data_pos))
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "the current BGZF block is not in the index",
            )
        })
}

/// Returns `current` moved by `offset`.
///
/// # Errors
///
/// Fails with [`io::ErrorKind::InvalidInput`] when the result is before the start (or
/// overflows), as [`std::io::Seek`] requires.
pub(crate) fn offset_position(current: u64, offset: i64) -> io::Result<u64> {
    current.checked_add_signed(offset).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid seek to a negative or overflowing position",
        )
    })
}

/// The error of seeking a BGZF stream relative to its end, whose uncompressed size is unknown.
pub(crate) fn end_unsupported() -> io::Error {
    io::Error::new(
        io::ErrorKind::Unsupported,
        "seeking relative to the end of a BGZF stream is unsupported",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_with_index() -> Result<(), Box<dyn std::error::Error>> {
        let index = gzi::Index::from(vec![(8, 21), (13, 55)]);
        let at = |c, u| VirtualPosition::try_from((c, u));

        assert_eq!(
            resolve_with_index(&index, at(13, 2)?, SeekFrom::Start(5))?,
            5
        );
        assert_eq!(
            resolve_with_index(&index, at(0, 3)?, SeekFrom::Current(0))?,
            3
        );
        assert_eq!(
            resolve_with_index(&index, at(8, 3)?, SeekFrom::Current(0))?,
            24
        );
        assert_eq!(
            resolve_with_index(&index, at(13, 2)?, SeekFrom::Current(-57))?,
            0
        );
        assert_eq!(
            resolve_with_index(&index, at(13, 2)?, SeekFrom::Current(10))?,
            67
        );

        let kind = |r: io::Result<u64>| r.map_err(|e| e.kind()).err();
        let before_start = resolve_with_index(&index, at(13, 2)?, SeekFrom::Current(-58));
        assert_eq!(kind(before_start), Some(io::ErrorKind::InvalidInput));
        let unindexed = resolve_with_index(&index, at(21, 0)?, SeekFrom::Current(0));
        assert_eq!(kind(unindexed), Some(io::ErrorKind::InvalidInput));
        let end = resolve_with_index(&index, at(0, 0)?, SeekFrom::End(0));
        assert_eq!(kind(end), Some(io::ErrorKind::Unsupported));

        Ok(())
    }
}
