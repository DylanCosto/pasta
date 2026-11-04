//! BAM header reader.

pub(crate) mod magic_number;
mod reference_sequences;
pub mod sam_header;

use std::io::{self, BufRead, Read};

use crate::sam::{self as sam, header::ReferenceSequences};

use self::{magic_number::read_magic_number, reference_sequences::read_reference_sequences};
use super::num::read_u32_le;
use crate::bam::io::MAGIC_NUMBER;

/// A BAM header reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the magic number.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use pasta::bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// let magic_number = header_reader.read_magic_number()?;
    /// assert_eq!(magic_number, *b"BAM\x01");
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        read_magic_number(&mut self.inner)
    }

    /// Returns a SAM header reader.
    ///
    /// The caller is responsible for discarding any extra padding in the header text, e.g., using
    /// [`sam_header::Reader::discard_to_end`].
    ///
    /// The position of the stream is expected to be directly after the magic number.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{io::Read, fs::File};
    /// use pasta::bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    ///
    /// let mut raw_sam_header_reader = header_reader.raw_sam_header_reader()?;
    ///
    /// let mut buf = Vec::new();
    /// raw_sam_header_reader.read_to_end(&mut buf)?;
    ///
    /// raw_sam_header_reader.discard_to_end()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn raw_sam_header_reader(&mut self) -> io::Result<sam_header::Reader<&mut R>> {
        let len = read_u32_le(&mut self.inner).map(u64::from)?;
        Ok(sam_header::Reader::new(&mut self.inner, len))
    }

    /// Reads the reference sequences.
    ///
    /// The position of the stream is expected to be directly after the raw SAM header.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{io, fs::File};
    /// use pasta::bam;
    ///
    /// let mut reader = File::open("sample.bam").map(bam::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    ///
    /// let mut raw_sam_header_reader = header_reader.raw_sam_header_reader()?;
    /// io::copy(&mut raw_sam_header_reader, &mut io::sink())?;
    /// raw_sam_header_reader.discard_to_end()?;
    ///
    /// let _reference_sequences = header_reader.read_reference_sequences()?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        read_reference_sequences(&mut self.inner)
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: Read,
{
    let mut header_reader = Reader::new(reader);
    read_header_inner(&mut header_reader)
}

fn read_header_inner<R>(reader: &mut Reader<R>) -> io::Result<sam::Header>
where
    R: Read,
{
    reader
        .read_magic_number()
        .and_then(magic_number::validate)?;

    let mut raw_sam_header_reader = reader.raw_sam_header_reader()?;
    let mut text = Vec::new();
    // This stops at a NUL, which pads the text in some files.
    raw_sam_header_reader.read_to_end(&mut text)?;
    raw_sam_header_reader.discard_to_end()?;

    let reference_sequences = reader.read_reference_sequences()?;

    // As in htslib, the binary reference sequence dictionary is authoritative: it decides the
    // reference sequences and their lengths, whatever the `@SQ` records say.
    let mut parser = sam::header::Parser::lenient();
    parser.set_binary_reference_sequences(reference_sequences);
    parser.allow_other_record_types();

    let mut lines = &text[..];
    let mut line = Vec::new();

    while read_line(&mut lines, &mut line)? != 0 {
        parser
            .parse_partial(&line)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    parser.set_text(text.into());

    Ok(parser.finish())
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use crate::sam::{
        self as sam,
        header::record::value::{
            Map,
            map::{self, header::Version},
        },
    };

    use super::*;

    fn put_u32_le(dst: &mut Vec<u8>, n: u32) {
        dst.extend(n.to_le_bytes());
    }

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut src = Vec::new();
        src.extend(MAGIC_NUMBER); // magic
        put_u32_le(&mut src, 27); // l_text
        src.extend(b"@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n"); // text
        put_u32_le(&mut src, 1); // n_ref
        put_u32_le(&mut src, 4); // ref[0].l_name
        src.extend(b"sq0\x00"); // ref[0].name
        put_u32_le(&mut src, 8); // ref[0].l_ref

        let mut reader = &src[..];
        let actual = read_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq0",
                Map::<map::ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_missing_sam_header_reference_sequence_dictionary() -> io::Result<()> {
        let mut src = Vec::new();
        src.extend(MAGIC_NUMBER); // magic
        put_u32_le(&mut src, 11); // l_text
        src.extend(b"@HD\tVN:1.6\n"); // text
        put_u32_le(&mut src, 1); // n_ref
        put_u32_le(&mut src, 4); // ref[0].l_name
        src.extend(b"sq0\x00"); // ref[0].name
        put_u32_le(&mut src, 8); // ref[0].l_ref

        let mut reader = &src[..];
        let actual = read_header(&mut reader)?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence(
                "sq0",
                Map::<map::ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    /// Builds a BAM header with the given text and binary reference sequence dictionary.
    fn bam_header(text: &[u8], reference_sequences: &[(&[u8], u32)]) -> Vec<u8> {
        let mut src = Vec::new();
        src.extend(MAGIC_NUMBER);
        put_u32_le(&mut src, u32::try_from(text.len()).unwrap());
        src.extend(text);
        put_u32_le(&mut src, u32::try_from(reference_sequences.len()).unwrap());

        for (name, len) in reference_sequences {
            put_u32_le(&mut src, u32::try_from(name.len() + 1).unwrap());
            src.extend(*name);
            src.push(0x00);
            put_u32_le(&mut src, *len);
        }

        src
    }

    /// The binary reference sequence dictionary of the header fixtures.
    const REFERENCE_SEQUENCES: [(&[u8], u32); 2] = [(b"chr1", 2000), (b"chr2", 1500)];

    fn names_and_lengths(header: &sam::Header) -> Vec<(&[u8], usize)> {
        header
            .reference_sequences()
            .iter()
            .map(|(name, map)| (name.as_slice(), usize::from(map.length())))
            .collect()
    }

    // Fixture `hdr_lnmismatch`: the text `LN` of chr2 is 1501, the binary length 1500.
    #[test]
    fn test_read_header_with_mismatched_reference_sequence_length() -> io::Result<()> {
        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1501\n@RG\tID:rg1\tSM:s1\n";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;

        // As in htslib, the binary dictionary is authoritative.
        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500)]
        );
        assert_eq!(header.read_groups().len(), 1);

        Ok(())
    }

    // Fixture `hdr_sqmissing`: chr2 has no `@SQ` record.
    #[test]
    fn test_read_header_with_missing_sam_header_reference_sequence() -> io::Result<()> {
        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@RG\tID:rg1\tSM:s1\n";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;

        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500)]
        );

        Ok(())
    }

    #[test]
    fn test_read_header_with_text_reference_sequences_htslib_overrides() -> io::Result<()> {
        use crate::sam::header::record::value::map::reference_sequence::tag;

        // A zero `LN` for a binary reference sequence (htslib's `view -c` does not check it), a
        // different order, and an `@SQ` record that is not in the binary dictionary, which htslib
        // appends.
        let text = b"@SQ\tSN:chr2\tLN:0\tM5:0123456789abcdef0123456789abcdef\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr3\tLN:10\n";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;

        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500), (b"chr3", 10)]
        );
        assert!(
            header.reference_sequences()[&b"chr2"[..]]
                .other_fields()
                .contains_key(&tag::MD5_CHECKSUM)
        );

        // No binary dictionary: the `@SQ` records are used.
        let src = bam_header(b"@SQ\tSN:chr1\tLN:2000\n", &[]);
        let header = read_header(&mut &src[..])?;
        assert_eq!(names_and_lengths(&header), [(&b"chr1"[..], 2000)]);

        Ok(())
    }

    // Fixture `hdr_duppg`.
    #[test]
    fn test_read_header_with_duplicate_program_ids() -> io::Result<()> {
        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n@RG\tID:rg1\tSM:s1\n@PG\tID:bwa\tPN:bwa\n@PG\tID:bwa\tPN:bwa\tVN:2\n";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;

        // htslib looks up the last record with an ID.
        let programs = header.programs().as_ref();
        assert_eq!(programs.len(), 1);
        assert!(
            programs[&b"bwa"[..]]
                .other_fields()
                .contains_key(&crate::sam::header::record::value::map::program::tag::VERSION)
        );

        // Both lines are written back, and a new program follows the last one.
        assert_eq!(write_header(&header)?, text);

        let mut header = header;
        header.add_program_record("ntools", Default::default())?;
        let mut expected = text.to_vec();
        expected.extend_from_slice(b"@PG\tID:ntools\tPN:ntools\tPP:bwa\n");
        assert_eq!(write_header(&header)?, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_writes_text_as_read_until_changed() -> io::Result<()> {
        use crate::sam::header::record::value::map::header::tag;

        // Fixtures `hdr_lnmismatch` and `hdr_sqmissing`, with an `@HD` that is not first.
        let text = b"@SQ\tSN:chr2\tLN:1501\tAS:x\n@HD\tVN:1.6\tSO:unsorted\tSO:x\n@CO\tc\n";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let mut header = read_header(&mut &src[..])?;

        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500)]
        );

        // As samtools does with `--no-PG`.
        assert_eq!(write_header(&header)?, text);

        // Once the header changes, it is written as htslib rebuilds it: `@HD` first, the binary
        // length in the first `LN`, and the missing `@SQ` after the last `@SQ`.
        header
            .header_mut()
            .as_mut()
            .expect("@HD")
            .other_fields_mut()
            .insert(tag::SORT_ORDER, "coordinate".into());

        assert_eq!(
            write_header(&header)?,
            b"@HD\tVN:1.6\tSO:coordinate\tSO:x\n@SQ\tSN:chr2\tLN:1500\tAS:x\n@SQ\tSN:chr1\tLN:2000\n@CO\tc\n"
        );

        Ok(())
    }

    #[test]
    fn test_read_header_with_duplicate_binary_names() {
        let src = bam_header(b"", &[(b"chr1", 2000), (b"chr1", 1500)]);
        assert!(matches!(
            read_header(&mut &src[..]),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    fn write_header(header: &sam::Header) -> io::Result<Vec<u8>> {
        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(header)?;
        Ok(writer.into_inner())
    }

    // Fixture `hdr_duprg`.
    #[test]
    fn test_read_header_with_duplicate_read_group_ids() -> io::Result<()> {
        use crate::sam::header::record::value::map::read_group::tag;

        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n@RG\tID:rg1\tSM:s1\n@RG\tID:rg1\tSM:s2\n";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;

        assert_eq!(header.read_groups().len(), 1);
        assert_eq!(
            header.read_groups()[&b"rg1"[..]]
                .other_fields()
                .get(&tag::SAMPLE)
                .map(|s| s.as_slice()),
            Some(&b"s1"[..])
        );

        Ok(())
    }

    // Fixture `hdr_novn`, and `VN:abc`, a late `@HD`, and an empty value.
    #[test]
    fn test_read_header_with_nonconforming_records() -> io::Result<()> {
        use crate::sam::header::record::value::map::header::tag;

        for text in [
            &b"@HD\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n@RG\tID:rg1\tSM:s1\n"[..],
            b"@HD\tVN:abc\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n",
            b"@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n@HD\tVN:1.6\tSO:coordinate\n@RG\tID:rg1\tSM:\n",
        ] {
            let src = bam_header(text, &REFERENCE_SEQUENCES);
            let header = read_header(&mut &src[..])?;

            let hd = header.header().expect("@HD is kept");
            assert_eq!(
                hd.other_fields()
                    .get(&tag::SORT_ORDER)
                    .map(|s| s.as_slice()),
                Some(&b"coordinate"[..])
            );
            assert_eq!(header.reference_sequences().len(), 2);
        }

        Ok(())
    }

    // Fixtures `hdr_empty` and `hdr_nul`.
    #[test]
    fn test_read_header_with_empty_or_nul_padded_text() -> io::Result<()> {
        let src = bam_header(b"", &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;
        assert!(header.header().is_none());
        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500)]
        );

        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:1500\n\x00\x00\x00\x00";
        let src = bam_header(text, &REFERENCE_SEQUENCES);
        let header = read_header(&mut &src[..])?;
        assert!(header.header().is_some());
        assert_eq!(header.reference_sequences().len(), 2);

        Ok(())
    }
}
