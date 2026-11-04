mod version;

use std::{error, fmt};

use bstr::{BStr, BString};

pub(crate) use self::version::parse_version;
use super::field::{consume_delimiter, consume_separator, parse_tag, parse_value, value};
use crate::sam::header::{
    parser::Context,
    record::value::{
        Map,
        map::{
            self, Header, OtherFields,
            header::{RawVersion, Tag, Version, tag},
            tag::Other,
        },
    },
};

/// An error returned when a SAM header header record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidField(super::field::ParseError),
    InvalidTag(super::field::tag::ParseError),
    InvalidValue(value::ParseError),
    MissingVersion,
    InvalidVersion(version::ParseError),
    InvalidOther(Other<tag::Standard>, value::ParseError),
    DuplicateTag(Tag),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::InvalidTag(e) => Some(e),
            Self::InvalidOther(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(_) => write!(f, "invalid field"),
            Self::InvalidTag(_) => write!(f, "invalid tag"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
            Self::MissingVersion => write!(f, "missing version ({}) field", tag::VERSION),
            Self::InvalidVersion(_) => write!(f, "invalid version ({})", tag::VERSION),
            Self::InvalidOther(tag, _) => write!(f, "invalid other ({tag})"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {tag}"),
        }
    }
}

pub(crate) fn parse_header(src: &mut &[u8], ctx: &Context) -> Result<Map<Header>, ParseError> {
    let mut version = None;
    let mut raw_version = None;

    let mut other_fields = OtherFields::new();

    while !src.is_empty() {
        consume_delimiter(src).map_err(ParseError::InvalidField)?;
        let tag = parse_tag(src).map_err(ParseError::InvalidTag)?;
        consume_separator(src).map_err(ParseError::InvalidField)?;

        match tag {
            tag::VERSION if ctx.is_lenient() => {
                let buf = parse_value(src, ctx).map_err(ParseError::InvalidValue)?;

                // Like htslib, use the first `VN` and accept any value.
                if version.is_none() && raw_version.is_none() {
                    match parse_version(buf) {
                        Ok(v) => version = Some(v),
                        Err(_) => raw_version = Some(RawVersion::Invalid(buf.into())),
                    }
                }
            }
            tag::VERSION => {
                parse_value(src, ctx)
                    .map_err(ParseError::InvalidValue)
                    .and_then(|buf| parse_version(buf).map_err(ParseError::InvalidVersion))
                    .and_then(|v| try_replace(&mut version, ctx, tag::VERSION, v))?;
            }
            Tag::Other(t) => parse_other(src, ctx, t)
                .and_then(|value| try_insert(&mut other_fields, ctx, t, value))?,
        }
    }

    let version = match version {
        Some(version) => version,
        None if ctx.is_lenient() => {
            raw_version.get_or_insert(RawVersion::Missing);
            Version::default()
        }
        None => return Err(ParseError::MissingVersion),
    };

    Ok(Map {
        inner: Header {
            version,
            raw_version,
        },
        other_fields,
    })
}

fn parse_other<'a>(
    src: &mut &'a [u8],
    ctx: &Context,
    tag: Other<tag::Standard>,
) -> Result<&'a BStr, ParseError> {
    parse_value(src, ctx).map_err(|e| ParseError::InvalidOther(tag, e))
}

fn try_replace<T>(
    option: &mut Option<T>,
    ctx: &Context,
    tag: Tag,
    value: T,
) -> Result<(), ParseError> {
    if option.replace(value).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(tag))
    } else {
        Ok(())
    }
}

fn try_insert<V>(
    other_fields: &mut OtherFields<tag::Standard>,
    ctx: &Context,
    tag: map::tag::Other<tag::Standard>,
    value: V,
) -> Result<(), ParseError>
where
    V: Into<BString>,
{
    if ctx.is_lenient() {
        // htslib looks up the first of duplicate tags.
        other_fields.entry(tag).or_insert_with(|| value.into());
        Ok(())
    } else if other_fields.insert(tag, value.into()).is_some() && !ctx.allow_duplicate_tags() {
        Err(ParseError::DuplicateTag(Tag::Other(tag)))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_header() {
        let mut src = &b"\tVN:1.6"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_header(&mut src, &ctx),
            Ok(Map::<Header>::new(Version::new(1, 6)))
        );
    }

    #[test]
    fn test_parse_header_with_missing_version() {
        let mut src = &b"\tSO:coordinate"[..];
        let ctx = Context::default();
        assert_eq!(
            parse_header(&mut src, &ctx),
            Err(ParseError::MissingVersion)
        );
    }

    #[test]
    fn test_parse_header_with_lenient_context() {
        let ctx = Context::lenient();

        let header = parse_header(&mut &b"\tSO:coordinate"[..], &ctx).unwrap();
        assert_eq!(header.version(), Version::default());
        assert_eq!(header.inner.raw_version, Some(RawVersion::Missing));
        assert_eq!(
            header.other_fields().get(&tag::SORT_ORDER),
            Some(&BString::from("coordinate"))
        );

        let header = parse_header(&mut &b"\tVN:abc"[..], &ctx).unwrap();
        assert_eq!(header.version(), Version::default());
        assert_eq!(
            header.inner.raw_version,
            Some(RawVersion::Invalid(BString::from("abc")))
        );

        let header = parse_header(&mut &b"\tVN:"[..], &ctx).unwrap();
        assert_eq!(
            header.inner.raw_version,
            Some(RawVersion::Invalid(BString::default()))
        );

        // The first of duplicate tags is used.
        let header = parse_header(
            &mut &b"\tVN:1.5\tSO:unsorted\tVN:1.6\tSO:coordinate"[..],
            &ctx,
        )
        .unwrap();
        assert_eq!(header, {
            let mut expected = Map::<Header>::new(Version::new(1, 5));
            expected
                .other_fields_mut()
                .insert(tag::SORT_ORDER, BString::from("unsorted"));
            expected
        });

        let header = parse_header(&mut &b"\tVN:abc\tVN:1.6"[..], &ctx).unwrap();
        assert_eq!(
            header.inner.raw_version,
            Some(RawVersion::Invalid(BString::from("abc")))
        );
    }

    #[test]
    fn test_parse_header_with_strict_context() {
        let ctx = Context::default();
        assert!(matches!(
            parse_header(&mut &b"\tVN:abc"[..], &ctx),
            Err(ParseError::InvalidVersion(_))
        ));
        assert_eq!(
            parse_header(&mut &b"\tVN:1.6\tSO:"[..], &ctx),
            Err(ParseError::InvalidOther(
                tag::SORT_ORDER,
                value::ParseError::Missing
            ))
        );
    }
}
