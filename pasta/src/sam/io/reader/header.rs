//! SAM header reader.

use std::io::{self, BufRead, Read};

use bstr::ByteSlice;

use super::read_line;
use crate::sam::{Header, header};

/// A SAM header reader.
///
/// This is created by calling [`super::Reader::header_reader`].
pub struct Reader<R> {
    inner: R,
    is_eol: bool,
}

impl<R> Reader<R> {
    pub(super) fn new(inner: R) -> Self {
        Self {
            inner,
            is_eol: true,
        }
    }
}

impl<R> Read for Reader<R>
where
    R: BufRead,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;

        if !src.is_empty() {
            self.is_eol = false;
        }

        self.consume(amt);

        Ok(amt)
    }
}

impl<R> BufRead for Reader<R>
where
    R: BufRead,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        const PREFIX: u8 = b'@';
        const LINE_FEED: u8 = b'\n';

        let src = self.inner.fill_buf()?;

        let buf = if self.is_eol && src.first().map(|&b| b != PREFIX).unwrap_or(true) {
            &[]
        } else if let Some(i) = src.as_bstr().find_byte(LINE_FEED) {
            self.is_eol = true;
            &src[..=i]
        } else {
            self.is_eol = false;
            src
        };

        Ok(buf)
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt);
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: BufRead,
{
    let mut reader = Reader::new(reader);

    // Accept what htslib accepts.
    let mut parser = header::Parser::lenient();
    let mut buf = Vec::new();

    while read_line(&mut reader, &mut buf)? != 0 {
        // htslib's SAM reader (`sam_hdr_build_from_sam_file`) accepts an untabbed `@CO` line as
        // an empty comment.
        if buf == b"@CO" {
            buf.push(b'\t');
        }

        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        buf.clear();
    }

    Ok(parser.finish())
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use super::*;
    use crate::sam::header::record::value::{
        Map,
        map::{self, ReferenceSequence, header::Version},
    };

    #[test]
    fn test_read_header_with_no_header() -> io::Result<()> {
        let data = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];
        assert!(read_header(&mut reader)?.is_empty());
        Ok(())
    }

    #[test]
    fn test_read_header_with_no_records() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n";
        let mut reader = data.as_bytes();

        let actual = read_header(&mut reader)?;

        let expected = crate::sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        use std::io::BufReader;

        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());

        let actual = read_header(&mut reader)?;

        let expected = crate::sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_nonconforming_records() -> io::Result<()> {
        // Duplicate `@PG` and `@RG` IDs, a late `@HD` without `VN`, an empty value, and an
        // untabbed `@CO`, all of which htslib's SAM reader accepts.
        let data = b"@SQ\tSN:sq0\tLN:8\n@HD\tSO:coordinate\n@RG\tID:rg0\tSM:\n@RG\tID:rg0\tSM:s1\n@PG\tID:pg0\n@PG\tID:pg0\n@CO\n*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];

        let header = read_header(&mut reader)?;
        assert!(header.header().is_some());
        assert_eq!(header.reference_sequences().len(), 1);
        assert_eq!(header.read_groups().len(), 1);
        let ids: Vec<_> = header.programs().as_ref().keys().collect();
        assert_eq!(ids, ["pg0"]);
        assert_eq!(header.comments(), [""]);

        // Lines are written back as read, but with the `@HD` first and a tab after `@CO`, as
        // htslib does.
        let mut writer = crate::sam::io::Writer::new(Vec::new());
        writer.write_header(&header)?;
        assert_eq!(
            writer.into_inner(),
            b"@HD\tSO:coordinate\n@SQ\tSN:sq0\tLN:8\n@RG\tID:rg0\tSM:\n@RG\tID:rg0\tSM:s1\n@PG\tID:pg0\n@PG\tID:pg0\n@CO\t\n"
        );

        Ok(())
    }

    #[test]
    fn test_read_header_with_invalid_reference_sequence_length() {
        // htslib rejects a non-positive `LN` in a SAM file.
        let mut reader = &b"@SQ\tSN:sq0\tLN:0\n"[..];
        assert!(matches!(
            read_header(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
