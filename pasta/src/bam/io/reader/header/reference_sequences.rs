mod reference_sequence;

use std::io::{self, Read};

use crate::sam::header::ReferenceSequences;

use self::reference_sequence::read_reference_sequence;
use crate::bam::io::reader::num::read_u32_le;

pub(super) fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: Read,
{
    let reference_sequence_count = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = ReferenceSequences::with_capacity(reference_sequence_count);

    for _ in 0..reference_sequence_count {
        let (name, reference_sequence) = read_reference_sequence(reader)?;

        // A duplicate name would shift the reference sequence IDs that follow it.
        if reference_sequences.contains_key(&name) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Duplicate entry \"{name}\" in target list"),
            ));
        }

        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use crate::sam::header::record::value::{Map, map::ReferenceSequence};
    use bstr::BString;

    use super::*;

    #[test]
    fn test_read_reference_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader)?;

        let expected: ReferenceSequences = [(
            BString::from("sq0"),
            Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
        )]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_reference_sequences_with_duplicate_names() {
        let data = [
            0x02, 0x00, 0x00, 0x00, // n_ref = 2
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
            0x04, 0x00, 0x00, 0x00, // ref[1].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[1].name = "sq0\x00"
            0x0d, 0x00, 0x00, 0x00, // ref[1].l_ref = 13
        ];

        let mut reader = &data[..];

        assert!(matches!(
            read_reference_sequences(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
