mod version;

use std::io::{self, Write};

use self::version::write_version_field;
use super::{write_field, write_other_fields};
use crate::sam::header::record::value::{
    Map,
    map::{
        Header,
        header::{RawVersion, tag},
    },
};

pub(crate) fn write_header<W>(writer: &mut W, header: &Map<Header>) -> io::Result<()>
where
    W: Write,
{
    // A header read leniently keeps a missing or unparsable `VN` as it was read, as htslib does.
    match header.raw_version() {
        None => write_version_field(writer, header.version())?,
        Some(RawVersion::Missing) => {}
        Some(RawVersion::Invalid(raw)) => write_field(writer, tag::VERSION, raw)?,
    }

    write_other_fields(writer, header.other_fields())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;
    use crate::sam::header::record::value::map::{
        OtherFields,
        header::{Version, tag::SORT_ORDER},
    };

    fn header(raw_version: Option<RawVersion>) -> Map<Header> {
        let mut other_fields = OtherFields::new();
        other_fields.insert(SORT_ORDER, BString::from("coordinate"));

        Map {
            inner: Header {
                version: Version::default(),
                raw_version,
            },
            other_fields,
        }
    }

    fn write(header: &Map<Header>) -> io::Result<Vec<u8>> {
        let mut buf = Vec::new();
        write_header(&mut buf, header)?;
        Ok(buf)
    }

    #[test]
    fn test_write_header_with_raw_version() -> io::Result<()> {
        assert_eq!(write(&header(None))?, b"\tVN:1.6\tSO:coordinate");
        assert_eq!(
            write(&header(Some(RawVersion::Missing)))?,
            b"\tSO:coordinate"
        );
        assert_eq!(
            write(&header(Some(RawVersion::Invalid(BString::from("abc")))))?,
            b"\tVN:abc\tSO:coordinate"
        );
        assert_eq!(
            write(&header(Some(RawVersion::Invalid(BString::default()))))?,
            b"\tVN:\tSO:coordinate"
        );
        Ok(())
    }

    #[test]
    fn test_write_header_after_version_mut() -> io::Result<()> {
        let mut header = header(Some(RawVersion::Missing));
        *header.version_mut() = Version::new(1, 5);
        assert_eq!(write(&header)?, b"\tVN:1.5\tSO:coordinate");
        Ok(())
    }
}
