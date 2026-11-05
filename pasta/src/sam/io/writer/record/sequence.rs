use std::io::{self, Write};

use super::MISSING;
use crate::bam::record::codec::decoder::decode_bases;
use crate::sam::alignment::record::{Sequence, SequenceRef, sequence_ref::FourBitPacked};

pub(super) fn write_sequence<W>(
    writer: &mut W,
    read_length: usize,
    sequence: SequenceRef<'_>,
) -> io::Result<()>
where
    W: Write,
{
    // § 1.4.10 "`SEQ`" (2022-08-22): "This field can be a '*' when the sequence is not stored."
    if sequence.is_empty() {
        writer.write_all(&[MISSING])?;
        return Ok(());
    }

    // § 1.4.10 "`SEQ`" (2022-08-22): "If not a '*', the length of the sequence must equal the sum
    // of lengths of `M`/`I`/`S`/`=`/`X` operations in `CIGAR`."
    if read_length > 0 && sequence.len() != read_length {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "read length-sequence length mismatch",
        ));
    }

    match sequence {
        SequenceRef::FourBitPacked(sequence) => write_four_bit_packed_sequence(writer, &sequence)?,
        SequenceRef::Raw(sequence) => write_raw_sequence(writer, sequence)?,
        SequenceRef::Sequence(sequence) => write_generic_sequence(writer, sequence)?,
    }

    Ok(())
}

fn write_four_bit_packed_sequence<W>(writer: &mut W, sequence: &FourBitPacked) -> io::Result<()>
where
    W: Write,
{
    // An even number of bases, so that each chunk starts on a byte.
    const CHUNK_SIZE: usize = 256;

    let src = sequence.as_ref();
    let base_count = sequence.len().min(src.len() * 2);

    let mut buf = [0; CHUNK_SIZE];

    for start in (0..base_count).step_by(CHUNK_SIZE) {
        let len = (base_count - start).min(CHUNK_SIZE);
        decode_bases(&src[start / 2..], &mut buf[..len]);
        writer.write_all(&buf[..len])?;
    }

    Ok(())
}

fn write_raw_sequence<W>(writer: &mut W, sequence: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if all_valid_bases(sequence) {
        writer.write_all(sequence)
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidInput))
    }
}

fn write_generic_sequence<W, S>(writer: &mut W, sequence: S) -> io::Result<()>
where
    W: Write,
    S: Sequence,
{
    for base in sequence.iter() {
        if !is_valid_base(base) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        writer.write_all(&[base])?;
    }

    Ok(())
}

// § 1.4 "The alignment section: mandatory fields" (2024-11-06): `[A-Za-z=.]+`.
fn is_valid_base(b: u8) -> bool {
    // Branch free, so that `all_valid_bases` vectorizes: setting bit 5 lowercases letters.
    let is_letter = (b | 0x20).wrapping_sub(b'a') < 26;
    is_letter | (b == b'=') | (b == b'.')
}

/// Checks bases in fixed-size chunks without short-circuiting within a chunk, which lets the
/// check vectorize while still stopping early on invalid input.
fn all_valid_bases(sequence: &[u8]) -> bool {
    const CHUNK_LEN: usize = 64;

    sequence.chunks(CHUNK_LEN).all(|chunk| {
        chunk
            .iter()
            .fold(true, |valid, &b| valid & is_valid_base(b))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::alignment::record_buf::Sequence as SequenceBuf;

    #[test]
    fn test_write_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        buf.clear();
        let sequence = SequenceBuf::default();
        let s = SequenceRef::Sequence(Box::new(&sequence));
        write_sequence(&mut buf, 0, s)?;
        assert_eq!(buf, b"*");

        buf.clear();
        let sequence = SequenceBuf::from(b"ACGT");
        let s = SequenceRef::Sequence(Box::new(&sequence));
        write_sequence(&mut buf, 4, s)?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        let sequence = SequenceBuf::from(b"ACGT");
        let s = SequenceRef::Sequence(Box::new(&sequence));
        write_sequence(&mut buf, 0, s)?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        let sequence = SequenceBuf::from(b"ACGT");
        let s = SequenceRef::Sequence(Box::new(&sequence));
        assert!(matches!(
            write_sequence(&mut buf, 1, s),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        buf.clear();
        let sequence = SequenceBuf::from(vec![b'!']);
        let s = SequenceRef::Sequence(Box::new(&sequence));
        assert!(matches!(
            write_sequence(&mut buf, 1, s),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }

    #[test]
    fn test_write_four_bit_packed_sequence() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, sequence: &FourBitPacked, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_four_bit_packed_sequence(buf, sequence)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &FourBitPacked::new(&[0x12, 0x40], 3), b"ACG")?;
        t(&mut buf, &FourBitPacked::new(&[0x12, 0x48], 4), b"ACGT")?;
        t(
            &mut buf,
            &FourBitPacked::new(&[0x12, 0x48, 0x00], 4),
            b"ACGT",
        )?;

        let src = [0x12, 0x48].repeat(96);
        let expected = b"ACGT".repeat(96);
        t(&mut buf, &FourBitPacked::new(&src, 384), &expected)?;

        Ok(())
    }

    #[test]
    fn test_write_raw_sequence() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_raw_sequence(&mut buf, b"ACGT")?;
        assert_eq!(buf, b"ACGT");

        buf.clear();
        assert!(matches!(
            write_raw_sequence(&mut buf, &[0xf0, 0x9f, 0x8d, 0x9c]),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid_base() {
        for b in (b'A'..=b'Z').chain(b'a'..=b'z').chain(*b"=.") {
            assert!(is_valid_base(b));
        }

        assert!(!is_valid_base(b'!'));
    }
}
