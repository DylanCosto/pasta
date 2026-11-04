use std::{io, num::NonZero, path::Path};

use crate::bgzf;
use crate::core::Position;
use crate::csi::binning_index::index::reference_sequence::bin::Chunk;
use crate::sam::{
    self as sam,
    alignment::Record as _,
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};
use bstr::ByteSlice;

use crate::bam::{Index, Record, index::Indexer, io::Reader};

/// Indexes a BAM file.
///
/// The input must be coordinate-sorted and marked as such in the SAM header, i.e., `SO:coordinate`.
///
/// This typically returns the index as a BAM index (BAI); however, if the length of a reference
/// sequence is longer than 2<sup>29</sup>, this returns a coordinate-sorted index (CSI) instead.
///
/// See also [`crate::bam::bai::fs::write`] and [`pasta::csi::fs::write`] to write the resulting
/// [`Index`] to a file.
///
/// # Examples
///
/// ```no_run
/// use pasta::bam;
/// let _index = bam::fs::index("sample.bam")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = super::open(src)?;
    index_inner(&mut reader)
}

fn index_inner<R>(reader: &mut Reader<R>) -> io::Result<Index>
where
    R: bgzf::io::Read,
{
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid sort order: expected {:?}, got {:?}",
                Some(COORDINATE.as_bstr()),
                header
                    .header()
                    .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
            ),
        ));
    }

    let mut record = Record::default();

    let max_reference_sequence_length = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length())
        .max()
        .unwrap_or(NonZero::<usize>::MIN);

    // SAFETY: `max_reference_sequence_length` is nonzero.
    let max_position = Position::new(max_reference_sequence_length.get()).unwrap();

    let mut indexer = Indexer::builder()
        .set_max_position_hint(max_position)
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match alignment_context(&record)? {
            Some((id, start, end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };

        indexer
            .add_record(alignment_context, chunk)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        start_position = end_position;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .map(|sort_order| sort_order == COORDINATE)
        .unwrap_or_default()
}

fn alignment_context(record: &Record) -> io::Result<Option<(usize, Position, Position)>> {
    let (Some(id), Some(start)) = (
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };

    // As HTSlib's `bam_endpos`: an unmapped record ends where it starts, whatever its CIGAR.
    let end = if record.flags().is_unmapped() {
        Some(start)
    } else {
        record.alignment_end().transpose()?
    };

    Ok(end.map(|end| (id, start, end)))
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use crate::csi::{BinningIndex, binning_index::ReferenceSequence as _};
    use crate::sam::{
        alignment::{RecordBuf, io::Write, record::Flags},
        header::record::value::{
            Map,
            map::{self, ReferenceSequence, builder::BuildError},
        },
    };
    use bstr::BString;

    use super::*;

    #[test]
    fn test_alignment_context_of_an_unmapped_record_with_a_cigar()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::sam::alignment::record::cigar::{Op, op::Kind};

        let header = sam::Header::builder()
            .set_reference_sequences(
                [(
                    BString::from("sq0"),
                    Map::<ReferenceSequence>::new(const { NonZero::new(100_000).unwrap() }),
                )]
                .into_iter()
                .collect(),
            )
            .build();

        let start = Position::try_from(16_380)?;
        let mut data = Vec::new();

        {
            let mut writer = crate::bam::io::Writer::from(&mut data);
            writer.write_header(&header)?;

            for flags in [Flags::empty(), Flags::UNMAPPED] {
                let record = RecordBuf::builder()
                    .set_flags(flags)
                    .set_reference_sequence_id(0)
                    .set_alignment_start(start)
                    .set_cigar([Op::new(Kind::Match, 50)].into_iter().collect())
                    .set_sequence(vec![b'A'; 50].into())
                    .build();

                writer.write_alignment_record(&header, &record)?;
            }
        }

        let mut reader = crate::bam::io::Reader::from(&data[..]);
        reader.read_header()?;
        let mut record = Record::default();

        // A mapped record spans its CIGAR; an unmapped one ends where it starts, as in HTSlib's
        // `bam_endpos`, whatever its CIGAR says.
        reader.read_record(&mut record)?;
        let end = Position::try_from(16_429)?;
        assert_eq!(alignment_context(&record)?, Some((0, start, end)));

        reader.read_record(&mut record)?;
        assert_eq!(alignment_context(&record)?, Some((0, start, start)));

        Ok(())
    }

    #[test]
    fn test_index() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = crate::bam::io::Writer::new(Vec::new());

        let header = sam::Header::builder()
            .set_header(
                Map::<map::Header>::builder()
                    .insert(SORT_ORDER, COORDINATE)
                    .build()?,
            )
            .set_reference_sequences(
                [(
                    BString::from("sq0"),
                    Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
                )]
                .into_iter()
                .collect(),
            )
            .build();

        writer.write_header(&header)?;

        let record = RecordBuf::builder()
            .set_flags(Flags::default())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .build();

        writer.write_alignment_record(&header, &record)?;

        let record = Record::default();
        writer.write_record(&header, &record)?;

        writer.try_finish()?;

        let data = writer.into_inner().into_inner();

        let mut reader = Reader::new(&data[..]);

        let Index::Bai(index) = index_inner(&mut reader)? else {
            panic!();
        };

        assert_eq!(index.min_shift(), 14);
        assert_eq!(index.depth(), 5);
        assert!(index.header().is_none());

        let reference_sequences = index.reference_sequences();
        assert_eq!(reference_sequences.len(), 1);

        let sq0 = &reference_sequences[0];
        let bins = sq0.bins();
        assert_eq!(bins.len(), 1);
        assert!(bins.get(&4681).is_some());

        assert_eq!(
            sq0.metadata().map(|metadata| (
                metadata.mapped_record_count(),
                metadata.unmapped_record_count()
            )),
            Some((1, 0))
        );

        assert_eq!(index.unplaced_unmapped_record_count(), Some(1));

        Ok(())
    }

    #[test]
    fn test_is_coordinate_sorted() -> Result<(), BuildError> {
        let header = sam::Header::default();
        assert!(!is_coordinate_sorted(&header));

        let header = sam::Header::builder()
            .set_header(
                Map::<map::Header>::builder()
                    .insert(SORT_ORDER, COORDINATE)
                    .build()?,
            )
            .build();

        assert!(is_coordinate_sorted(&header));

        Ok(())
    }
}
