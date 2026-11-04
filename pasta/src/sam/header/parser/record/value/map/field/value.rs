use std::{error, fmt};

use bstr::{BStr, ByteSlice};

use crate::sam::header::parser::Context;

/// An error returned when a SAM header record field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is missing.
    Missing,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing => write!(f, "missing input"),
        }
    }
}

/// Parses a field value, which runs to the next tab or the end of the line.
///
/// The SAM specification requires a value (`[ -~]+`), but htslib accepts an empty one (e.g.,
/// `SM:`), so an empty value is only an error when the context is not lenient.
pub fn parse_value<'a>(src: &mut &'a [u8], ctx: &Context) -> Result<&'a BStr, ParseError> {
    const DELIMITER: u8 = b'\t';

    let i = src.as_bstr().find_byte(DELIMITER).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);

    *src = rest;

    if buf.is_empty() && !ctx.is_lenient() {
        Err(ParseError::Missing)
    } else {
        Ok(buf.as_bstr())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() {
        let ctx = Context::default();

        let mut src = &b"ndls"[..];
        assert_eq!(parse_value(&mut src, &ctx), Ok(b"ndls".as_bstr()));

        let mut src = &b""[..];
        assert_eq!(parse_value(&mut src, &ctx), Err(ParseError::Missing));

        let mut src = &b"\tSM:s0"[..];
        assert_eq!(parse_value(&mut src, &ctx), Err(ParseError::Missing));
    }

    #[test]
    fn test_parse_value_with_lenient_context() {
        let ctx = Context::lenient();

        let mut src = &b""[..];
        assert_eq!(parse_value(&mut src, &ctx), Ok(b"".as_bstr()));

        let mut src = &b"\tSM:s0"[..];
        assert_eq!(parse_value(&mut src, &ctx), Ok(b"".as_bstr()));
        assert_eq!(src, b"\tSM:s0");
    }
}
