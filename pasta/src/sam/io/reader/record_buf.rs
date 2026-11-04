//! SAM record buf reader.

pub mod cigar;
pub mod data;
pub mod flags;
pub mod mapping_quality;
pub mod name;
pub mod position;
pub mod quality_scores;
pub mod reference_sequence_id;
pub mod sequence;
pub mod template_length;

use std::{
    error, fmt,
    io::{self, BufRead},
};

use self::{
    cigar::parse_cigar, data::parse_data, mapping_quality::parse_mapping_quality, name::parse_name,
    position::parse_alignment_start, quality_scores::parse_quality_scores,
    sequence::parse_sequence,
};
pub(crate) use self::{flags::parse_flags, template_length::parse_template_length};
use super::{normalize::Normalizer, read_line};
use crate::sam::{Header, alignment::RecordBuf};

pub(crate) fn read_record_buf<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    header: &Header,
    normalizer: &mut Normalizer,
    record: &mut RecordBuf,
) -> io::Result<usize>
where
    R: BufRead,
{
    buf.clear();

    match read_line(reader, buf)? {
        0 => Ok(0),
        n => {
            parse_record_buf(buf, header, normalizer, record)?;
            Ok(n)
        }
    }
}

/// An error when a raw SAM record fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The name is invalid.
    InvalidName(name::ParseError),
    /// The flags are invalid.
    InvalidFlags(flags::ParseError),
    /// The reference sequence ID is invalid.
    InvalidReferenceSequenceId(reference_sequence_id::ParseError),
    /// The position is invalid.
    InvalidPosition(position::ParseError),
    /// The mapping quality is invalid.
    InvalidMappingQuality(mapping_quality::ParseError),
    /// The CIGAR is invalid.
    InvalidCigar(cigar::ParseError),
    /// The mate reference sequence ID is invalid.
    InvalidMateReferenceSequenceId(reference_sequence_id::ParseError),
    /// The mate position is invalid.
    InvalidMatePosition(position::ParseError),
    /// The template length is invalid.
    InvalidTemplateLength(template_length::ParseError),
    /// The sequence is invalid.
    InvalidSequence(sequence::ParseError),
    /// The quality scores are invalid.
    InvalidQualityScores(quality_scores::ParseError),
    /// The data is invalid.
    InvalidData(data::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidName(e) => Some(e),
            Self::InvalidFlags(e) => Some(e),
            Self::InvalidReferenceSequenceId(e) => Some(e),
            Self::InvalidPosition(e) => Some(e),
            Self::InvalidMappingQuality(e) => Some(e),
            Self::InvalidCigar(e) => Some(e),
            Self::InvalidMateReferenceSequenceId(e) => Some(e),
            Self::InvalidMatePosition(e) => Some(e),
            Self::InvalidTemplateLength(e) => Some(e),
            Self::InvalidSequence(e) => Some(e),
            Self::InvalidQualityScores(e) => Some(e),
            Self::InvalidData(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidName(_) => write!(f, "invalid name"),
            Self::InvalidFlags(_) => write!(f, "invalid flags"),
            Self::InvalidReferenceSequenceId(_) => write!(f, "invalid reference sequence ID"),
            Self::InvalidPosition(_) => write!(f, "invalid position"),
            Self::InvalidMappingQuality(_) => write!(f, "invalid mapping quality"),
            Self::InvalidCigar(_) => write!(f, "invalid CIGAR"),
            Self::InvalidMateReferenceSequenceId(_) => {
                write!(f, "invalid mate reference sequence ID")
            }
            Self::InvalidMatePosition(_) => write!(f, "invalid mate position"),
            Self::InvalidTemplateLength(_) => write!(f, "invalid template length"),
            Self::InvalidSequence(_) => write!(f, "invalid sequence"),
            Self::InvalidQualityScores(_) => write!(f, "invalid quality scores"),
            Self::InvalidData(_) => write!(f, "invalid data"),
        }
    }
}

pub(crate) fn parse_record_buf(
    mut src: &[u8],
    header: &Header,
    normalizer: &mut Normalizer,
    record: &mut RecordBuf,
) -> io::Result<()> {
    fn invalid_data(e: ParseError) -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, e)
    }

    const MISSING: &[u8] = b"*";
    const EQ: &[u8] = b"=";

    match next_field(&mut src) {
        MISSING => *record.name_mut() = None,
        field => {
            parse_name(field, record.name_mut())
                .map(Some)
                .map_err(|e| invalid_data(ParseError::InvalidName(e)))?;
        }
    };

    let field = next_field(&mut src);
    *record.flags_mut() =
        parse_flags(field).map_err(|e| invalid_data(ParseError::InvalidFlags(e)))?;

    let field = next_field(&mut src);
    let (reference_sequence_id, is_unknown) =
        normalizer.reference_sequence_id(header, field, false)?;
    *record.reference_sequence_id_mut() = reference_sequence_id;

    let field = next_field(&mut src);
    *record.alignment_start_mut() =
        parse_alignment_start(field).map_err(|e| invalid_data(ParseError::InvalidPosition(e)))?;

    // Unlike HTSlib, the position of a reference sequence missing from the header is dropped
    // with it (see `normalize`).
    if is_unknown {
        *record.alignment_start_mut() = None;
    }

    let field = next_field(&mut src);
    *record.mapping_quality_mut() = parse_mapping_quality(field)
        .map_err(|e| invalid_data(ParseError::InvalidMappingQuality(e)))?;

    record.cigar_mut().as_mut().clear();
    let field = next_field(&mut src);
    // HTSlib reads any CIGAR starting with `*` as missing.
    let is_cigar_missing = field.starts_with(MISSING);
    if !is_cigar_missing {
        parse_cigar(field, record.cigar_mut())
            .map_err(|e| invalid_data(ParseError::InvalidCigar(e)))?;
    }

    normalizer.normalize_record_buf(record, is_cigar_missing);

    // RNEXT `=` is the RNAME as written, even if POS 0 dropped it (see `normalize`).
    let (mate_reference_sequence_id, is_mate_unknown) = match next_field(&mut src) {
        EQ => (reference_sequence_id, is_unknown),
        field => normalizer.reference_sequence_id(header, field, true)?,
    };
    *record.mate_reference_sequence_id_mut() = mate_reference_sequence_id;

    let field = next_field(&mut src);
    *record.mate_alignment_start_mut() = parse_alignment_start(field)
        .map_err(|e| invalid_data(ParseError::InvalidMatePosition(e)))?;

    if is_mate_unknown {
        *record.mate_alignment_start_mut() = None;
    }

    normalizer.normalize_record_buf_mate(record);

    let field = next_field(&mut src);
    *record.template_length_mut() = parse_template_length(field)
        .map_err(|e| invalid_data(ParseError::InvalidTemplateLength(e)))?;

    record.sequence_mut().as_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_sequence(field, record.sequence_mut())
            .map_err(|e| invalid_data(ParseError::InvalidSequence(e)))?;
    }

    record.quality_scores_mut().as_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_quality_scores(field, record.sequence().len(), record.quality_scores_mut())
            .map_err(|e| invalid_data(ParseError::InvalidQualityScores(e)))?;
    }

    record.data_mut().clear();
    parse_data(src, record.data_mut()).map_err(|e| invalid_data(ParseError::InvalidData(e)))?;

    Ok(())
}

fn next_field<'a>(src: &mut &'a [u8]) -> &'a [u8] {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let (field, rest) = if let Some(i) = memchr(DELIMITER, src) {
        let (field, rest) = src.split_at(i);
        (field, &rest[1..])
    } else {
        src.split_at(src.len())
    };

    *src = rest;

    field
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use bstr::ByteSlice;

    use super::*;
    use crate::sam::alignment::record::Flags;

    fn header() -> Header {
        use crate::sam::header::record::value::{Map, map::ReferenceSequence};

        let mut sq1 = Map::<ReferenceSequence>::new(const { NonZero::new(13).unwrap() });
        sq1.other_fields_mut().insert(
            crate::sam::header::record::value::map::reference_sequence::tag::ALTERNATIVE_NAMES,
            "1,one".into(),
        );

        Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence("sq1", sq1)
            .build()
    }

    fn parse(header: &Header, s: &[u8]) -> io::Result<RecordBuf> {
        let mut record = RecordBuf::default();
        parse_record_buf(s, header, &mut Normalizer::default(), &mut record)?;
        Ok(record)
    }

    #[test]
    fn test_parse_with_data() -> io::Result<()> {
        use crate::sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};

        let header = Header::default();
        let s = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tNH:i:1\tCO:Z:ndls";
        let record = parse(&header, s)?;

        let expected = RecordBuf::builder()
            .set_data(
                [
                    (Tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
                    (Tag::COMMENT, Value::from("ndls")),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        assert_eq!(record, expected);

        Ok(())
    }

    #[test]
    fn test_parse_mate_reference_sequence_id() -> io::Result<()> {
        let header = header();

        let mate = |rnext: &str| -> io::Result<Option<usize>> {
            let s = format!("r\t1\tsq0\t1\t60\t4M\t{rnext}\t5\t0\tACGT\tNDLS");
            parse(&header, s.as_bytes()).map(|record| record.mate_reference_sequence_id())
        };

        assert_eq!(mate("=")?, Some(0));
        assert_eq!(mate("sq0")?, Some(0));
        assert_eq!(mate("sq1")?, Some(1));
        assert_eq!(mate("one")?, Some(1));
        assert_eq!(mate("*")?, None);
        // htslib: "unrecognized mate reference name"; treated as unmapped.
        assert_eq!(mate("sq2")?, None);

        Ok(())
    }

    // The normalization of htslib's `sam_parse1`.
    #[test]
    fn test_parse_normalizes_as_sam_parse1() -> io::Result<()> {
        let header = header();

        // (line, flags, reference sequence ID, mate reference sequence ID)
        type Case = (&'static [u8], u16, Option<usize>, Option<usize>);

        let cases: [Case; 11] = [
            // RNAME `*` without 0x4
            (b"a\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", 0x4, None, None),
            (b"b\t65\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", 0x45, None, None),
            // placed, CIGAR `*`
            (
                b"c\t0\tsq0\t10\t30\t*\t*\t0\t0\tACGT\tIIII",
                0x4,
                Some(0),
                None,
            ),
            // RNAME with POS 0
            (b"d\t0\tsq0\t0\t0\t4M\t*\t0\t0\tACGT\tIIII", 0x4, None, None),
            // normal
            (
                b"e\t0\tsq0\t2\t30\t4M\t*\t0\t0\tACGT\tIIII",
                0x0,
                Some(0),
                None,
            ),
            // RNEXT with PNEXT 0
            (
                b"f\t1\tsq0\t3\t30\t4M\tsq0\t0\t0\tACGT\tIIII",
                0x1,
                Some(0),
                None,
            ),
            // RNEXT `=` of a record with POS 0: the mate stays on sq0 (HTSlib drops it)
            (
                b"g\t1\tsq0\t0\t30\t4M\t=\t4\t0\tACGT\tIIII",
                0x5,
                None,
                Some(0),
            ),
            // unknown RNAME, RNEXT `=`
            (
                b"h\t1\tsqX\t5\t30\t4M\t=\t4\t0\tACGT\tIIII",
                0x5,
                None,
                None,
            ),
            // unknown RNEXT
            (
                b"i\t1\tsq0\t5\t30\t4M\tsqY\t4\t0\tACGT\tIIII",
                0x1,
                Some(0),
                None,
            ),
            // alternative names
            (
                b"j\t1\t1\t5\t30\t4M\tone\t4\t0\tACGT\tIIII",
                0x1,
                Some(1),
                Some(1),
            ),
            // already unmapped
            (
                b"k\t4\tsq1\t5\t0\t*\t=\t5\t0\tACGT\tIIII",
                0x4,
                Some(1),
                Some(1),
            ),
        ];

        for (line, flags, id, mate_id) in cases {
            let record = parse(&header, line)?;

            // Unlike HTSlib, a position on a reference sequence missing from the header is
            // dropped with it.
            if line.contains_str("sqX") {
                assert!(record.alignment_start().is_none());
                assert!(record.mate_alignment_start().is_none());
            }

            if line.contains_str("sqY") {
                assert!(record.mate_alignment_start().is_none());
            }

            assert_eq!(record.flags(), Flags::from(flags), "{}", line.as_bstr());
            assert_eq!(record.reference_sequence_id(), id, "{}", line.as_bstr());
            assert_eq!(
                record.mate_reference_sequence_id(),
                mate_id,
                "{}",
                line.as_bstr()
            );
        }

        Ok(())
    }

    #[test]
    fn test_parse_with_a_reference_sequence_name_and_no_reference_sequences() {
        let header = Header::default();
        let s = b"r\t0\tsq0\t1\t60\t4M\t*\t0\t0\tACGT\tNDLS";
        assert!(parse(&header, s).is_err());

        // RNEXT alone is looked up (and not found) without an error.
        let s = b"r\t4\t*\t0\t0\t*\tsq0\t1\t0\tACGT\tNDLS";
        assert!(parse(&header, s).is_ok());
    }
}
