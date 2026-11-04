//! Binning index indexer.

mod add_record_error;

use crate::bgzf;
use crate::core::Position;
use indexmap::IndexMap;

pub use self::add_record_error::AddRecordError;
use super::{
    calculate_max_position,
    index::{
        Header, Index, ReferenceSequence,
        reference_sequence::{self, bin::Chunk},
    },
};

/// A binning index indexer.
#[derive(Debug)]
pub struct Indexer<I> {
    min_shift: u8,
    depth: u8,
    header: Option<Header>,
    reference_sequences: Vec<ReferenceSequence<I>>,
    // For each reference sequence, the start of the first record overlapping each linear window
    // (`1 << min_shift` positions), or `None`. A reference sequence's list is empty until its
    // first record.
    window_offsets: Vec<Vec<Option<bgzf::VirtualPosition>>>,
    // The reference sequence of the previous placed record.
    current_reference_sequence_id: Option<usize>,
    unplaced_unmapped_record_count: u64,
    prev_start_position: Position,
    max_position: Position,
}

impl<I> Indexer<I>
where
    I: reference_sequence::Index + Default,
{
    /// Creates a binning index indexer.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::csi::binning_index::{index::reference_sequence::index::BinnedIndex, Indexer};
    ///
    /// let indexer = Indexer::<BinnedIndex>::new(14, 5).expect("invalid min shift and/or depth");
    /// ```
    pub fn new(min_shift: u8, depth: u8) -> Option<Self> {
        let max_position = calculate_max_position(min_shift, depth)?;

        Some(Self {
            min_shift,
            depth,
            header: None,
            reference_sequences: Vec::new(),
            window_offsets: Vec::new(),
            current_reference_sequence_id: None,
            unplaced_unmapped_record_count: 0,
            prev_start_position: Position::MIN,
            max_position,
        })
    }

    /// Sets a tabix header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::csi::binning_index::{
    ///     index::{reference_sequence::index::BinnedIndex, Header},
    ///     Indexer,
    /// };
    ///
    /// let header = Header::default();
    /// let indexer = Indexer::<BinnedIndex>::default().set_header(header);
    /// ```
    pub fn set_header(mut self, header: Header) -> Self {
        self.header = Some(header);
        self
    }

    /// Adds a record.
    ///
    /// As in HTSlib, each reference sequence's records must be sorted by start position and form
    /// one continuous block, but the blocks may come in any order of reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::bgzf;
    /// use pasta::core::Position;
    /// use pasta::csi::binning_index::{
    ///     index::reference_sequence::{bin::Chunk, index::BinnedIndex},
    ///     Indexer,
    /// };
    ///
    /// let mut indexer = Indexer::<BinnedIndex>::default();
    ///
    /// let reference_sequence_id = 0;
    /// let start = Position::try_from(8)?;
    /// let end = Position::try_from(13)?;
    /// let is_mapped = true;
    /// let chunk = Chunk::new(
    ///     bgzf::VirtualPosition::from(144),
    ///     bgzf::VirtualPosition::from(233),
    /// );
    ///
    /// indexer.add_record(Some((reference_sequence_id, start, end, is_mapped)), chunk)?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_record(
        &mut self,
        alignment_context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> Result<(), AddRecordError> {
        let Some((reference_sequence_id, start, end, is_mapped)) = alignment_context else {
            self.unplaced_unmapped_record_count += 1;
            return Ok(());
        };

        if start > end {
            return Err(AddRecordError::InvalidInterval { start, end });
        }

        if end > self.max_position {
            return Err(AddRecordError::EndExceedsMaxPosition {
                end,
                max_position: self.max_position,
            });
        }

        if self.current_reference_sequence_id != Some(reference_sequence_id) {
            if let Some(current_reference_sequence_id) = self.current_reference_sequence_id
                && self
                    .window_offsets
                    .get(reference_sequence_id)
                    .is_some_and(|offsets| !offsets.is_empty())
            {
                return Err(AddRecordError::OutOfOrderReferenceSequenceId {
                    reference_sequence_id,
                    current_reference_sequence_id,
                });
            }

            if reference_sequence_id >= self.reference_sequences.len() {
                self.add_reference_sequences_until(reference_sequence_id);
            }

            self.current_reference_sequence_id = Some(reference_sequence_id);
            self.prev_start_position = Position::MIN;
        }

        if start < self.prev_start_position {
            return Err(AddRecordError::OutOfOrderStartPosition {
                start,
                prev_start_position: self.prev_start_position,
            });
        }

        self.prev_start_position = start;

        let reference_sequence = &mut self.reference_sequences[reference_sequence_id];
        reference_sequence.update(self.min_shift, self.depth, start, end, is_mapped, chunk);

        // Records arrive in file order, so the first record to reach a window has its lowest
        // offset.
        let window_offsets = &mut self.window_offsets[reference_sequence_id];
        let first_window = (usize::from(start) - 1) >> self.min_shift;
        let last_window = (usize::from(end) - 1) >> self.min_shift;

        if window_offsets.len() <= last_window {
            window_offsets.resize(last_window + 1, None);
        }

        for offset in &mut window_offsets[first_window..=last_window] {
            offset.get_or_insert(chunk.start());
        }

        Ok(())
    }

    /// Builds a binning index.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::csi::binning_index::{index::reference_sequence::index::BinnedIndex, Indexer};
    /// let indexer = Indexer::<BinnedIndex>::default();
    /// let index = indexer.build(0);
    /// ```
    pub fn build(mut self, reference_sequence_count: usize) -> Index<I> {
        if reference_sequence_count > self.reference_sequences.len() {
            // SAFETY: `reference_sequence_count` is nonzero.
            self.add_reference_sequences_until(reference_sequence_count - 1);
        }

        let reference_sequences = self
            .reference_sequences
            .into_iter()
            .zip(self.window_offsets)
            .map(|(reference_sequence, window_offsets)| {
                reference_sequence.optimize(self.depth, &window_offsets)
            })
            .collect();

        let mut builder = Index::builder()
            .set_min_shift(self.min_shift)
            .set_depth(self.depth)
            .set_reference_sequences(reference_sequences)
            .set_unplaced_unmapped_record_count(self.unplaced_unmapped_record_count);

        if let Some(header) = self.header {
            builder = builder.set_header(header);
        }

        builder.build()
    }

    fn add_reference_sequences_until(&mut self, reference_sequence_id: usize) {
        let new_len = reference_sequence_id + 1;

        self.reference_sequences.resize_with(new_len, || {
            ReferenceSequence::new(IndexMap::new(), Default::default(), None)
        });

        self.window_offsets.resize_with(new_len, Vec::new);
    }
}

impl<I> Default for Indexer<I>
where
    I: reference_sequence::Index + Default,
{
    fn default() -> Self {
        Self {
            min_shift: 14,
            depth: 5,
            header: None,
            reference_sequences: Vec::new(),
            window_offsets: Vec::new(),
            current_reference_sequence_id: None,
            unplaced_unmapped_record_count: 0,
            prev_start_position: Position::MIN,
            max_position: const { Position::new(1 << 29).unwrap() },
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bgzf;

    use super::*;
    use crate::csi::binning_index::index::reference_sequence::{Bin, Metadata, index::LinearIndex};

    #[test]
    fn test_default() {
        let indexer = Indexer::<LinearIndex>::default();

        assert_eq!(indexer.min_shift, 14);
        assert_eq!(indexer.depth, 5);
        assert!(indexer.header.is_none());
        assert!(indexer.reference_sequences.is_empty());
        assert_eq!(indexer.unplaced_unmapped_record_count, 0);
        assert_eq!(indexer.prev_start_position, Position::MIN);
        assert_eq!(
            indexer.max_position,
            const { Position::new(1 << 29).unwrap() }
        );
    }

    #[test]
    fn test_add_record_with_invalid_interval() {
        let mut indexer = Indexer::<LinearIndex>::default();

        let start = const { Position::new(13).unwrap() };
        let end = const { Position::new(8).unwrap() };
        let alignment_context = Some((0, start, end, true));
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
        );

        assert_eq!(
            indexer.add_record(alignment_context, chunk),
            Err(AddRecordError::InvalidInterval { start, end })
        );
    }

    #[test]
    fn test_add_record_with_end_exceeds_max_position() {
        let mut indexer = Indexer::<LinearIndex>::new(1, 1).unwrap(); // max position = 16

        let alignment_context = Some((0, Position::MIN, Position::MAX, true));
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
        );

        assert_eq!(
            indexer.add_record(alignment_context, chunk),
            Err(AddRecordError::EndExceedsMaxPosition {
                end: Position::MAX,
                max_position: const { Position::new(16).unwrap() }
            })
        );
    }

    #[test]
    fn test_add_record_with_out_of_order_reference_sequence_id() -> Result<(), AddRecordError> {
        let mut indexer = Indexer::<LinearIndex>::default();

        let start = const { Position::new(8).unwrap() };
        let end = const { Position::new(13).unwrap() };

        for (i, reference_sequence_id) in [1, 0].into_iter().enumerate() {
            let chunk = Chunk::new(
                bgzf::VirtualPosition::from(i as u64),
                bgzf::VirtualPosition::from(i as u64 + 1),
            );
            indexer.add_record(Some((reference_sequence_id, start, end, true)), chunk)?;
        }

        let alignment_context = Some((1, start, end, true));
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(2),
            bgzf::VirtualPosition::from(3),
        );
        assert_eq!(
            indexer.add_record(alignment_context, chunk),
            Err(AddRecordError::OutOfOrderReferenceSequenceId {
                reference_sequence_id: 1,
                current_reference_sequence_id: 0
            })
        );

        Ok(())
    }

    #[test]
    fn test_build_with_reference_sequence_blocks_out_of_order() -> Result<(), AddRecordError> {
        // Each reference sequence's records, as (start, end, length in the file).
        const RECORDS: [&[(usize, usize, u64)]; 3] = [
            &[(5, 8, 9), (70000, 70100, 12)],
            &[(1, 2, 5)],
            &[(3, 5, 7), (20000, 20004, 3), (20001, 150000, 70000)],
        ];

        // The records of each reference sequence in `order`, as (reference sequence ID, start,
        // end, chunk), laid out one block after the other in the file.
        let layout = |order: &[usize]| {
            let mut offset = 0;
            let mut records = Vec::new();

            for &id in order {
                for &(start, end, len) in RECORDS[id] {
                    let start = Position::new(start).expect("start is nonzero");
                    let end = Position::new(end).expect("end is nonzero");
                    let chunk = Chunk::new(
                        bgzf::VirtualPosition::from(offset << 16),
                        bgzf::VirtualPosition::from((offset + len) << 16),
                    );
                    records.push((id, start, end, chunk));
                    offset += len;
                }
            }

            records
        };

        let index = |records: &[(usize, Position, Position, Chunk)]| {
            let mut indexer = Indexer::<LinearIndex>::default();

            for &(id, start, end, chunk) in records {
                indexer.add_record(Some((id, start, end, true)), chunk)?;
            }

            Ok::<_, AddRecordError>(indexer.build(RECORDS.len() + 1))
        };

        // The blocks of reference sequences 2, 0 and 1, in that order, index as they would one
        // reference sequence at a time.
        let mut records = layout(&[2, 0, 1]);
        let actual = index(&records)?;
        records.sort_by_key(|&(id, ..)| id);
        let expected = index(&records)?;

        assert_eq!(actual, expected);
        assert_eq!(actual.reference_sequences().len(), RECORDS.len() + 1);

        Ok(())
    }

    #[test]
    fn test_add_record_with_out_of_order_start_position() -> Result<(), AddRecordError> {
        let mut indexer = Indexer::<LinearIndex>::default();

        let alignment_context = Some((
            0,
            const { Position::new(5).unwrap() },
            const { Position::new(8).unwrap() },
            true,
        ));
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(0),
            bgzf::VirtualPosition::from(1),
        );
        indexer.add_record(alignment_context, chunk)?;

        let alignment_context = Some((
            1,
            const { Position::new(3).unwrap() },
            const { Position::new(8).unwrap() },
            true,
        ));
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(1),
            bgzf::VirtualPosition::from(2),
        );
        indexer.add_record(alignment_context, chunk)?;

        let alignment_context = Some((1, Position::MIN, const { Position::new(5).unwrap() }, true));
        let chunk = Chunk::new(
            bgzf::VirtualPosition::from(2),
            bgzf::VirtualPosition::from(3),
        );
        assert_eq!(
            indexer.add_record(alignment_context, chunk),
            Err(AddRecordError::OutOfOrderStartPosition {
                start: Position::MIN,
                prev_start_position: const { Position::new(3).unwrap() }
            })
        );

        Ok(())
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let mut indexer = Indexer::<LinearIndex>::default();

        indexer.add_record(
            Some((0, Position::try_from(8)?, Position::try_from(13)?, true)),
            Chunk::new(
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
            ),
        )?;

        indexer.add_record(
            Some((
                0,
                Position::try_from(121393)?,
                Position::try_from(196418)?,
                false,
            )),
            Chunk::new(
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(3473408),
            ),
        )?;

        let actual = indexer.build(1);

        let expected = {
            let bins = [
                (
                    4681,
                    Bin::new(vec![Chunk::new(
                        bgzf::VirtualPosition::from(0),
                        bgzf::VirtualPosition::from(9),
                    )]),
                ),
                (
                    73,
                    Bin::new(vec![Chunk::new(
                        bgzf::VirtualPosition::from(9),
                        bgzf::VirtualPosition::from(3473408),
                    )]),
                ),
            ]
            .into_iter()
            .collect();

            let index = vec![
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
                bgzf::VirtualPosition::from(9),
            ];

            let metadata = Metadata::new(
                bgzf::VirtualPosition::from(0),
                bgzf::VirtualPosition::from(3473408),
                1,
                1,
            );

            Index::builder()
                .set_reference_sequences(vec![ReferenceSequence::new(bins, index, Some(metadata))])
                .set_unplaced_unmapped_record_count(0)
                .build()
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_build_with_reference_sequence_count() {
        let index = Indexer::<LinearIndex>::default().build(2);
        assert_eq!(index.reference_sequences().len(), 2);
    }

    #[test]
    fn test_build_with_no_reference_sequences() {
        let header = crate::csi::binning_index::index::header::Builder::vcf().build();

        let actual = Indexer::<LinearIndex>::default()
            .set_header(header.clone())
            .build(0);

        let expected = Index::builder()
            .set_header(header)
            .set_unplaced_unmapped_record_count(0)
            .build();

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_build_with_nonstandard_bin_factors() {
        const MIN_SHIFT: u8 = 4;
        const DEPTH: u8 = 2;

        let actual = Indexer::<LinearIndex>::new(MIN_SHIFT, DEPTH)
            .unwrap()
            .build(0);

        let expected = Index::builder()
            .set_min_shift(MIN_SHIFT)
            .set_depth(DEPTH)
            .set_unplaced_unmapped_record_count(0)
            .build();

        assert_eq!(actual, expected);
    }
}
