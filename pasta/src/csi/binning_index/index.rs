//! Coordinate-sorted index and fields.

mod builder;
pub mod header;
pub mod reference_sequence;

use std::io;

use crate::bgzf;
use crate::core::{Position, region::Interval};

pub use self::{builder::Builder, header::Header, reference_sequence::ReferenceSequence};
use super::{BinningIndex, calculate_max_position, index::reference_sequence::bin::Chunk};

/// A binning index.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Index<I> {
    min_shift: u8,
    depth: u8,
    header: Option<Header>,
    reference_sequences: Vec<ReferenceSequence<I>>,
    unplaced_unmapped_record_count: Option<u64>,
}

impl<I> Index<I>
where
    I: reference_sequence::Index,
{
    /// Returns a builder to create an index from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::csi;
    /// let builder = csi::Index::builder();
    /// ```
    pub fn builder() -> Builder<I> {
        Builder::default()
    }

    /// Returns a list of indexed reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::csi;
    /// let index = csi::Index::default();
    /// assert!(index.reference_sequences().is_empty());
    /// ```
    pub fn reference_sequences(&self) -> &[ReferenceSequence<I>] {
        &self.reference_sequences
    }
}

impl<I> Default for Index<I>
where
    I: reference_sequence::Index,
{
    fn default() -> Self {
        Self::builder().build()
    }
}

impl<I> BinningIndex for Index<I>
where
    I: reference_sequence::Index,
{
    fn min_shift(&self) -> u8 {
        self.min_shift
    }

    fn depth(&self) -> u8 {
        self.depth
    }

    fn header(&self) -> Option<&Header> {
        self.header.as_ref()
    }

    fn reference_sequences(&self) -> Box<dyn Iterator<Item = &dyn super::ReferenceSequence> + '_> {
        Box::new(
            self.reference_sequences
                .iter()
                .map(|reference_sequence| reference_sequence as &dyn super::ReferenceSequence),
        )
    }

    fn unplaced_unmapped_record_count(&self) -> Option<u64> {
        self.unplaced_unmapped_record_count
    }

    fn query(&self, reference_sequence_id: usize, interval: Interval) -> io::Result<Vec<Chunk>> {
        use super::optimize_chunks;

        let reference_sequence = self
            .reference_sequences()
            .get(reference_sequence_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence ID: {reference_sequence_id}"),
                )
            })?;

        let Some((start, _)) = resolve_interval(self.min_shift(), self.depth(), interval)? else {
            return Ok(Vec::new());
        };

        let query_bins = reference_sequence
            .query(self.min_shift(), self.depth(), interval)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        let min_offset = reference_sequence.min_offset(self.min_shift(), self.depth(), start);
        let merged_chunks = optimize_chunks(&chunks, min_offset);

        Ok(merged_chunks)
    }

    fn last_first_record_start_position(&self) -> Option<bgzf::VirtualPosition> {
        self.reference_sequences
            .iter()
            .rev()
            .find_map(|rs| rs.first_record_in_last_linear_bin_start_position())
    }
}

/// Resolves `interval` to the positions an index with `min_shift` and `depth` can hold, as
/// htslib's `hts_itr_query` does: an end past the index's maximum position (or no end) is clamped
/// to it, e.g. `chr1:1-999999999` for a whole chromosome, and an interval that starts past it
/// holds nothing (`None`).
fn resolve_interval<I>(
    min_shift: u8,
    depth: u8,
    interval: I,
) -> io::Result<Option<(Position, Position)>>
where
    I: Into<Interval>,
{
    let interval = interval.into();

    let start = interval.start().unwrap_or(Position::MIN);
    let max_position = calculate_max_position(min_shift, depth)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "unsupported max position"))?;

    if start > max_position {
        return Ok(None);
    }

    let end = interval
        .end()
        .map_or(max_position, |end| end.min(max_position));

    Ok(Some((start, end)))
}

#[cfg(test)]
mod query_offset_tests {
    use super::*;
    use crate::csi::binning_index::{Indexer, index::reference_sequence::index::BinnedIndex};

    #[test]
    fn test_query_generated_csi_preserves_spanning_records()
    -> Result<(), Box<dyn std::error::Error>> {
        // Keep every chunk large enough that index optimization cannot fold away the leaf bins.
        let starts = [8, 100, 130, 200, 300, 500];
        let ends = [400, 110, 240, 210, 310, 510];
        let position = |n| bgzf::VirtualPosition::new(n * 100000, 0).unwrap();
        let chunk = |i| Chunk::new(position(i as u64 + 1), position(i as u64 + 2));
        let mut indexer = Indexer::<BinnedIndex>::new(4, 2).unwrap();
        for (i, (start, end)) in starts.iter().zip(ends).enumerate() {
            let (start, end) = (Position::new(*start).unwrap(), Position::new(end).unwrap());
            indexer.add_record(Some((0, start, end, true)), chunk(i))?;
        }

        let original = indexer.build(1);
        let mut writer = crate::csi::io::Writer::new(Vec::new());
        writer.write_index(&original)?;
        let data = writer.into_inner().finish()?;
        let decoded = crate::csi::io::Reader::new(&data[..]).read_index()?;

        for index in [&original, &decoded] {
            for query in 1..=1024 {
                let interval = Position::new(query).unwrap()..=Position::new(query).unwrap();
                let chunks = index.query(0, interval.into())?;

                for (i, (&start, end)) in starts.iter().zip(ends).enumerate() {
                    let pos = chunk(i).start();
                    let is_read = chunks.iter().any(|c| c.start() <= pos && pos < c.end());
                    // A query reads the records that span it, but not the long first record (in
                    // bin 0, ending at 400) in a later window that it does not reach.
                    assert!(
                        is_read || !(start..=end).contains(&query),
                        "lost {i} at {query}"
                    );
                    assert!(
                        !is_read || i > 0 || !(497..=512).contains(&query),
                        "{query}"
                    );
                }
            }
        }

        Ok(())
    }

    /// Like htslib, an end past the index's maximum position (2^29 for BAI) is clamped
    /// to it, and a start past it gives no chunks, instead of failing the query.
    #[test]
    fn test_query_clamps_to_max_position() -> Result<(), Box<dyn std::error::Error>> {
        use crate::csi::binning_index::index::reference_sequence::index::LinearIndex;

        const MIN_SHIFT: u8 = 14;
        const DEPTH: u8 = 5;

        let mut indexer = Indexer::<LinearIndex>::new(MIN_SHIFT, DEPTH).unwrap();
        let chunk = Chunk::new(
            bgzf::VirtualPosition::new(100, 0).unwrap(),
            bgzf::VirtualPosition::new(200, 0).unwrap(),
        );
        let (start, end) = (Position::new(1000).unwrap(), Position::new(2000).unwrap());
        indexer.add_record(Some((0, start, end, true)), chunk)?;
        let index = indexer.build(1);

        let max_position = Position::new(1 << 29).unwrap();
        let whole = index.query(0, (..).into())?;
        assert_eq!(whole, [chunk]);

        for end in [
            max_position,
            Position::new(999_999_999).unwrap(),
            Position::MAX,
        ] {
            assert_eq!(
                index.query(0, (Position::MIN..=end).into())?,
                whole,
                "{end}"
            );
        }

        let past_max = max_position.checked_add(1).unwrap();
        assert!(index.query(0, (past_max..).into())?.is_empty());
        assert!(
            index
                .query(0, (past_max..=Position::MAX).into())?
                .is_empty()
        );

        // Still an error: a reference sequence the index does not have.
        assert!(index.query(1, (..).into()).is_err());

        Ok(())
    }
}
