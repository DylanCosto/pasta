use std::{
    io::{self, BufRead, Read},
    vec,
};

use crate::bgzf;

use super::IndexedRecords;
use crate::csi::binning_index::index::{Header, reference_sequence::bin::Chunk};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// How a format reader classifies a record read from a query.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RecordMatch {
    /// The record does not intersect the region.
    Skip,
    /// The record intersects the region.
    Include,
    /// The record starts past the region, as does every later record in its chunk.
    NextChunk,
}

/// A query reader.
///
/// This reader returns the uncompressed data between all the given chunks.
pub struct Query<'r, R> {
    reader: &'r mut R,
    chunks: vec::IntoIter<Chunk>,
    chunks_are_sorted: bool,
    state: State,
}

impl<'r, R> Query<'r, R> {
    /// Creates a query reader.
    pub fn new(reader: &'r mut R, chunks: Vec<Chunk>) -> Self {
        Self {
            reader,
            chunks_are_sorted: chunks_are_sorted(&chunks),
            chunks: chunks.into_iter(),
            state: State::Seek,
        }
    }

    /// Applies the match of a complete record just read, returning whether to yield it.
    ///
    /// A record past the region ends the query only when every later chunk follows the current
    /// one in the file. Otherwise, only the rest of the current chunk is skipped, preserving
    /// custom chunk order and overlapping chunks, including their duplicate records.
    #[doc(hidden)]
    pub fn accept(&mut self, record_match: RecordMatch) -> bool {
        match record_match {
            RecordMatch::Include => return true,
            RecordMatch::NextChunk if self.chunks_are_sorted => self.state = State::Done,
            RecordMatch::NextChunk => self.skip_current_chunk(),
            RecordMatch::Skip => {}
        }
        false
    }

    /// Discards unread bytes in the active chunk, retaining all later chunks.
    ///
    /// Call this only after consuming a complete record. If no chunk is active, this is a no-op.
    #[doc(hidden)]
    pub fn skip_current_chunk(&mut self) {
        if matches!(self.state, State::Read(_)) {
            self.state = State::Seek;
        }
    }

    /// Sets the query aside, releasing the reader: [`Self::resume`] carries it on with the same
    /// reader, which must not be read or moved in between.
    pub fn suspend(self) -> SuspendedQuery {
        SuspendedQuery {
            chunks: self.chunks,
            chunks_are_sorted: self.chunks_are_sorted,
            state: self.state,
        }
    }

    /// Carries on a query set aside by [`Self::suspend`] with the reader it was reading.
    pub fn resume(reader: &'r mut R, suspended: SuspendedQuery) -> Self {
        let SuspendedQuery {
            chunks,
            chunks_are_sorted,
            state,
        } = suspended;
        Self {
            reader,
            chunks,
            chunks_are_sorted,
            state,
        }
    }
}

/// A query set aside without its reader ([`Query::suspend`]).
pub struct SuspendedQuery {
    chunks: vec::IntoIter<Chunk>,
    chunks_are_sorted: bool,
    state: State,
}

pub(crate) fn chunks_are_sorted(chunks: &[Chunk]) -> bool {
    chunks
        .windows(2)
        .all(|pair| pair[0].end() <= pair[1].start())
}

impl<R> Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Creates an iterator that parses indexed records.
    pub fn indexed_records(self, header: &Header) -> IndexedRecords<Self> {
        IndexedRecords::new(self, header)
    }
}

impl<R> Read for Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = match self.chunks.next() {
                        Some(chunk) => {
                            self.reader.seek_to_virtual_position(chunk.start())?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    }
                }
                State::Read(chunk_end) => {
                    if self.reader.virtual_position() < chunk_end {
                        return self.reader.fill_buf();
                    } else {
                        self.state = State::Seek;
                    }
                }
                State::Done => return Ok(&[]),
            }
        }
    }

    fn consume(&mut self, amt: usize) {
        self.reader.consume(amt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chunk(start: u64, end: u64) -> Chunk {
        Chunk::new(
            bgzf::VirtualPosition::from(start),
            bgzf::VirtualPosition::from(end),
        )
    }

    #[test]
    fn test_accept_ends_query_for_chunks_in_file_order() {
        let mut reader = bgzf::io::Reader::new(io::Cursor::new(Vec::new()));
        let mut query = Query::new(&mut reader, vec![chunk(0, 1), chunk(2, 3)]);

        assert!(query.accept(RecordMatch::Include));
        assert!(!query.accept(RecordMatch::Skip));
        assert!(!query.accept(RecordMatch::NextChunk));
        assert!(matches!(query.state, State::Done));
    }

    #[test]
    fn test_suspend_and_resume_carry_on_the_query() -> io::Result<()> {
        // Two blocks of four bytes; the query reads the second half of the first and the
        // second.
        let mut data = Vec::new();
        let mut second = 0;
        for block in [b"abcd", b"efgh"] {
            second = data.len() as u64;
            let mut writer = bgzf::io::Writer::new(Vec::new());
            io::Write::write_all(&mut writer, block)?;
            let mut compressed = writer.finish()?;
            // Without the EOF block.
            compressed.truncate(compressed.len() - 28);
            data.extend(compressed);
        }
        let mut reader = bgzf::io::Reader::new(io::Cursor::new(data));
        let chunks = vec![chunk(2, 4), chunk(second << 16, (second << 16) | 4)];
        let mut suspended = Query::new(&mut reader, chunks).suspend();
        let mut buf = [0; 1];
        let mut read = Vec::new();
        loop {
            // A byte at a time, setting the query aside after each.
            let mut query = Query::resume(&mut reader, suspended);
            if query.read(&mut buf)? == 0 {
                break;
            }
            read.push(buf[0]);
            suspended = query.suspend();
        }
        assert_eq!(read, b"cdefgh");
        Ok(())
    }

    #[test]
    fn test_accept_skips_only_the_current_chunk_for_unordered_chunks() {
        let mut reader = bgzf::io::Reader::new(io::Cursor::new(Vec::new()));
        let mut query = Query::new(&mut reader, vec![chunk(2, 3), chunk(0, 1)]);
        query.state = State::Read(bgzf::VirtualPosition::from(3));

        assert!(!query.accept(RecordMatch::NextChunk));
        assert!(matches!(query.state, State::Seek));
    }
}
