use std::io::{self, Read};

use crate::csi::binning_index::index::reference_sequence::{Bin, Metadata};
use indexmap::IndexMap;

use crate::tabix::io::reader::num::{read_i32_le, read_u32_le};

pub(super) fn read_bins<R>(reader: &mut R) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: Read,
{
    use crate::csi::io::reader::index::reference_sequences::{
        bins::{invalid_bin_id_error, read_chunks},
        initial_capacity, read_metadata,
    };

    use crate::tabix::index::DEPTH;

    const MAX_ID: usize = Bin::max_id(DEPTH).unwrap();
    const METADATA_ID: usize = Bin::metadata_id(DEPTH).unwrap();

    let bin_count = read_bin_count(reader)?;

    let mut bins = IndexMap::with_capacity(initial_capacity(bin_count));
    let mut metadata = None;

    for _ in 0..bin_count {
        let id = read_u32_le(reader).and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let is_duplicate = if id == METADATA_ID {
            let m =
                read_metadata(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            metadata.replace(m).is_some()
        } else if id > MAX_ID {
            return Err(invalid_bin_id_error(id, MAX_ID));
        } else {
            let chunks =
                read_chunks(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let bin = Bin::new(chunks);

            bins.insert(id, bin).is_some()
        };

        if is_duplicate {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate bin ID: {id}"),
            ));
        }
    }

    Ok((bins, metadata))
}

fn read_bin_count<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    read_i32_le(reader)
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_bins() -> io::Result<()> {
        let src = [
            0x00, 0x00, 0x00, 0x00, // n_bin = 0
        ];
        let (actual_bins, actual_metadata) = read_bins(&mut &src[..])?;
        assert!(actual_bins.is_empty());
        assert!(actual_metadata.is_none());

        let src = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x00, 0x00, 0x00, 0x00, // bins[0].bin = 0
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
            0x4a, 0x92, 0x00, 0x00, // bins[1].bin = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_unmapped = 0
        ];
        let (actual_bins, actual_metadata) = read_bins(&mut &src[..])?;
        assert_eq!(actual_bins.len(), 1);
        assert!(actual_bins.get(&0).is_some());
        assert!(actual_metadata.is_some());

        let src = [
            0x01, 0x00, 0x00, 0x00, // n_bin = 1
        ];
        assert!(matches!(
            read_bins(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let src = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x00, 0x00, 0x00, 0x00, // bins[0].bin = 0
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
            0x00, 0x00, 0x00, 0x00, // bins[1].bin = 0
            0x00, 0x00, 0x00, 0x00, // bins[1].n_chunk = 0
        ];
        assert!(matches!(
            read_bins(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        let src = [
            0x02, 0x00, 0x00, 0x00, // n_bin = 2
            0x4a, 0x92, 0x00, 0x00, // bins[0].bin = 37450
            0x02, 0x00, 0x00, 0x00, // bins[0].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[0].metadata.n_unmapped = 0
            0x4a, 0x92, 0x00, 0x00, // bins[1].bin = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].metadata.n_unmapped = 0
        ];
        assert!(matches!(
            read_bins(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_read_bins_with_an_out_of_range_bin_id() {
        for id in [37449u32, 37451, 201364042] {
            let mut src = vec![0x01, 0x00, 0x00, 0x00]; // n_bin = 1
            src.extend(id.to_le_bytes()); // bins[0].bin
            src.extend([0x00, 0x00, 0x00, 0x00]); // bins[0].n_chunk = 0

            assert!(matches!(
                read_bins(&mut &src[..]),
                Err(e) if e.kind() == io::ErrorKind::InvalidData
                    && e.to_string() == format!("invalid bin ID: {id} (max 37448)")
            ));
        }
    }

    #[test]
    fn test_read_bins_with_a_bin_count_past_the_data() {
        let src = [
            0xff, 0xff, 0xff, 0x7f, // n_bin = 2147483647
        ];
        assert!(matches!(
            read_bins(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
    }
}
