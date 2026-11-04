use std::io::{self, Write};

use crate::sam::header::record::value::map::reference_sequence::tag;

pub(super) fn write_name_field<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: Write,
{
    use crate::sam::io::writer::header::record::{value::map::write_separator, write_delimiter};

    write_delimiter(writer)?;
    writer.write_all(tag::NAME.as_ref())?;
    write_separator(writer)?;
    write_value(writer, name)?;

    Ok(())
}

fn write_value<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if is_valid_name(name) {
        writer.write_all(name)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid reference sequence name",
        ))
    }
}

// The SAM specification restricts names (§ 1.2.1 "Character set restrictions": `[:rname:∧*=][:rname:]*`),
// but htslib reads and writes any `SN` value, including an empty one, so this only rejects what
// would break the line, as for other values.
fn is_valid_name(name: &[u8]) -> bool {
    super::super::is_valid_value(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_valid_name() {
        assert!(is_valid_name(b"sq0"));
        assert!(is_valid_name(b"sq0*"));
        assert!(is_valid_name(b"sq0="));
        assert!(is_valid_name(b""));
        assert!(is_valid_name(b"sq 0"));
        assert!(is_valid_name(b"*sq0"));

        assert!(!is_valid_name(b"sq\t0"));
        assert!(!is_valid_name(b"sq0\n"));
    }
}
