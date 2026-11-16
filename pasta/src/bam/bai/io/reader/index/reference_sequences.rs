mod bins;
mod intervals;

use std::io::{self, Read};

use crate::csi::{
    binning_index::index::{ReferenceSequence, reference_sequence::index::LinearIndex},
    io::reader::index::reference_sequences::MAX_INITIAL_CAPACITY,
};

use self::{bins::read_bins, intervals::read_intervals};
use crate::bam::io::reader::num::read_u32_le;

pub(super) fn read_reference_sequences<R>(
    reader: &mut R,
) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: Read,
{
    let n_ref = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut references = Vec::with_capacity(n_ref.min(MAX_INITIAL_CAPACITY));

    for _ in 0..n_ref {
        let (bins, metadata) = read_bins(reader)?;
        let intervals = read_intervals(reader)?;
        references.push(ReferenceSequence::new(bins, intervals, metadata));
    }

    Ok(references)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_reference_sequences_with_a_count_past_the_data() {
        let data = [
            0xff, 0xff, 0xff, 0xff, // n_ref = 4294967295
        ];
        let mut reader = &data[..];
        assert!(matches!(
            read_reference_sequences(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
    }
}
