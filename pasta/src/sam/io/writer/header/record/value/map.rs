mod header;
mod program;
mod read_group;
mod reference_sequence;
mod tag;

use std::io::{self, Write};

use self::tag::write_tag;
pub(crate) use self::{
    header::write_header, program::write_program, read_group::write_read_group,
    reference_sequence::write_reference_sequence,
};
use crate::sam::header::record::value::map::OtherFields;

const SEPARATOR: u8 = b':';

fn write_field<W, K>(writer: &mut W, tag: K, value: &[u8]) -> io::Result<()>
where
    W: Write,
    K: AsRef<[u8; 2]>,
{
    use crate::sam::io::writer::header::record::write_delimiter;

    write_delimiter(writer)?;
    write_tag(writer, *tag.as_ref())?;
    write_separator(writer)?;
    write_value(writer, value)?;

    Ok(())
}

fn write_other_fields<W, S>(writer: &mut W, other_fields: &OtherFields<S>) -> io::Result<()>
where
    W: Write,
{
    for (tag, value) in other_fields {
        write_field(writer, tag, value)?;
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[SEPARATOR])
}

fn write_value<W>(writer: &mut W, value: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if is_valid_value(value) {
        writer.write_all(value)
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid value"))
    }
}

pub(super) fn is_valid_value(buf: &[u8]) -> bool {
    // § 1.3 "The header section" (2023-05-24): `[ -~]+`. htslib reads and writes any value,
    // including an empty one (e.g., `SM:`) and other bytes (e.g., UTF-8), so that a header read
    // leniently can be written, this only rejects what would break the line.
    !buf.iter().any(|&b| matches!(b, b'\t' | b'\n' | 0x00))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_valid_value() {
        assert!(is_valid_value(b"sq0"));
        assert!(is_valid_value(b""));
        assert!(!is_valid_value(&[0x00]));
        assert!(!is_valid_value(b"sq\t0"));
    }
}
