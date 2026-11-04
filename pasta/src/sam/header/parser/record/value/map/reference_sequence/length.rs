use std::{error, fmt, num::NonZero};

/// An error returned when a SAM header reference sequence record length value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid(lexical_core::Error),
    /// The value is 0.
    Zero,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
            Self::Zero => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(_) => write!(f, "invalid input"),
            Self::Zero => write!(f, "invalid value: 0"),
        }
    }
}

pub(super) fn parse_length(src: &mut &[u8]) -> Result<NonZero<usize>, ParseError> {
    let (n, i) = lexical_core::parse_partial::<usize>(src).map_err(ParseError::Invalid)?;
    *src = &src[i..];
    NonZero::new(n).ok_or(ParseError::Zero)
}

/// Parses an `LN` value as htslib does (`strtoll`): leading whitespace and a sign are allowed, and
/// anything after the digits is ignored. The length must be positive.
pub(super) fn parse_length_lenient(src: &[u8]) -> Result<NonZero<usize>, ParseError> {
    let src = src.trim_ascii_start();

    let (is_negative, src) = match src.split_first() {
        Some((b'-', rest)) => (true, rest),
        Some((b'+', rest)) => (false, rest),
        _ => (false, src),
    };

    let end = src
        .iter()
        .position(|b| !b.is_ascii_digit())
        .unwrap_or(src.len());
    let digits = &src[..end];

    if digits.is_empty() {
        return Err(ParseError::Zero);
    }

    let n: usize = lexical_core::parse(digits).map_err(ParseError::Invalid)?;

    if is_negative {
        Err(ParseError::Zero)
    } else {
        NonZero::new(n).ok_or(ParseError::Zero)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_length() {
        let mut src = &b"8"[..];
        assert_eq!(
            parse_length(&mut src),
            Ok(const { NonZero::new(8).unwrap() })
        );

        let mut src = &b"538522340430300790495419781092981030533"[..];
        assert!(matches!(
            parse_length(&mut src),
            Err(ParseError::Invalid(_))
        ));

        let mut src = &b"0"[..];
        assert_eq!(parse_length(&mut src), Err(ParseError::Zero));
    }

    #[test]
    fn test_parse_length_lenient() {
        assert_eq!(
            parse_length_lenient(b"8"),
            Ok(const { NonZero::new(8).unwrap() })
        );
        assert_eq!(
            parse_length_lenient(b" +8kb"),
            Ok(const { NonZero::new(8).unwrap() })
        );

        assert_eq!(parse_length_lenient(b""), Err(ParseError::Zero));
        assert_eq!(parse_length_lenient(b"0"), Err(ParseError::Zero));
        assert_eq!(parse_length_lenient(b"-5"), Err(ParseError::Zero));
        assert_eq!(parse_length_lenient(b"NA"), Err(ParseError::Zero));
        assert!(matches!(
            parse_length_lenient(b"538522340430300790495419781092981030533"),
            Err(ParseError::Invalid(_))
        ));
    }
}
