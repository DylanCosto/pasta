//! Region strings resolved against a header or index (HTSlib's `hts_parse_region` and
//! `hts_parse_reg64`): `chr`, `chr:beg`, `chr:beg-end`, `chr:beg-`, `chr:-end`, `{name}:...` for
//! names containing colons, and positions with commas and `k`/`m`/`g` suffixes.

use crate::variant::{
    log,
    number::{ParseDecimalFlags, parse_decimal},
};

/// The exclusive end meaning "to the end of the sequence" (HTSlib's `HTS_POS_MAX`).
pub const POS_MAX: i64 = ((i32::MAX as i64) << 32) | i32::MAX as i64;

/// How [`parse_region`] reads a region (HTSlib's `HTS_PARSE_*` flags).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct RegionFlags {
    /// The string is a comma-separated list; positions take no comma separators
    /// (`HTS_PARSE_LIST`).
    pub list: bool,
    /// `chr:100` is the single base 100, `chr:-100` is `chr:1-100`, and `chr:100-` runs to the
    /// end (`HTS_PARSE_ONE_COORD`, bcftools' convention).
    pub one_coord: bool,
}

/// A parsed region: the target ID and a 0-based half-open interval.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct ParsedRegion {
    /// The target ID.
    pub tid: u32,
    /// The 0-based start.
    pub beg: i64,
    /// The exclusive end ([`POS_MAX`] for the end of the sequence).
    pub end: i64,
    /// The number of bytes consumed, including a trailing comma in list mode.
    pub consumed: usize,
}

/// Parses a region against a name lookup (HTSlib's `hts_parse_region`).
///
/// Returns `None` (after reporting the reason as HTSlib does) for an unknown name or a
/// malformed or empty interval.
pub fn parse_region<F>(s: &[u8], name_to_id: F, flags: RegionFlags) -> Option<ParsedRegion>
where
    F: Fn(&[u8]) -> Option<u32>,
{
    let dec_flags = ParseDecimalFlags {
        thousands_separator: !flags.list,
    };
    let mut s_len = s.len();
    let mut s_end = s.len();
    let mut colon: Option<usize> = None;
    let mut quoted = 0usize;
    let mut start = 0usize;
    if s.first() == Some(&b'{') {
        let Some(close) = s.iter().position(|&c| c == b'}') else {
            log::error(
                "hts_parse_region",
                format_args!("Mismatching braces in \"{}\"", String::from_utf8_lossy(s)),
            );
            return None;
        };
        start = 1;
        s_len -= 1;
        if s.get(close + 1) == Some(&b':') {
            colon = Some(close + 1);
        }
        quoted = 1;
        if flags.list
            && let Some(comma) = s[close..]
                .iter()
                .position(|&c| c == b',')
                .map(|p| p + close)
        {
            s_len = comma - start;
            s_end = comma + 1;
        }
    } else {
        if flags.list
            && let Some(comma) = s.iter().position(|&c| c == b',')
        {
            s_len = comma;
            s_end = comma + 1;
        }
        colon = s[..s_len].iter().rposition(|&c| c == b':');
    }
    let name_all = &s[start..start + s_len - quoted];

    let Some(colon) = colon else {
        let tid = name_to_id(name_all)?;
        return Some(ParsedRegion {
            tid,
            beg: 0,
            end: POS_MAX,
            consumed: s_end,
        });
    };

    if quoted == 0 {
        let whole = &s[start..start + s_len];
        if let Some(tid) = name_to_id(whole) {
            if name_to_id(&s[start..colon]).is_some() {
                log::error(
                    "hts_parse_region",
                    format_args!(
                        "Range is ambiguous. Use {{{}}} or {{{}}}{} instead",
                        String::from_utf8_lossy(whole),
                        String::from_utf8_lossy(&s[start..colon]),
                        String::from_utf8_lossy(&s[colon..])
                    ),
                );
                return None;
            }
            return Some(ParsedRegion {
                tid,
                beg: 0,
                end: POS_MAX,
                consumed: s_end,
            });
        }
    }
    let tid = name_to_id(&s[start..colon - quoted])?;

    let at = |i: usize| s.get(i).copied().unwrap_or(0);
    let (b, mut hyphen) = parse_decimal(s, colon + 1, dec_flags);
    let mut beg = b - 1;
    let mut end;
    if beg < 0 {
        // HTSlib's own test: `chr1:0-5` (beg == -1) is accepted, `chr1:-5-10` is not.
        if beg != -1 && at(hyphen) == b'-' && at(colon + 1) != 0 {
            log::error("hts_parse_region", format_args!("Coordinates must be > 0"));
            return None;
        }
        let c = at(hyphen);
        if c.is_ascii_digit() || c == 0 || c == b',' {
            end = if beg == -1 { POS_MAX } else { -(beg + 1) };
            return Some(ParsedRegion {
                tid,
                beg: 0,
                end,
                consumed: s_end,
            });
        } else if beg < -1 {
            log::error(
                "hts_parse_region",
                format_args!(
                    "Unexpected string \"{}\" after region",
                    String::from_utf8_lossy(&s[hyphen..])
                ),
            );
            return None;
        }
    }
    let c = at(hyphen);
    if c == 0 || (flags.list && c == b',') {
        end = if flags.one_coord { beg + 1 } else { POS_MAX };
    } else if c == b'-' {
        let (e, h) = parse_decimal(s, hyphen + 1, dec_flags);
        end = e;
        hyphen = h;
        if at(hyphen) != 0 && at(hyphen) != b',' {
            log::error(
                "hts_parse_region",
                format_args!(
                    "Unexpected string \"{}\" after region",
                    String::from_utf8_lossy(&s[hyphen..])
                ),
            );
            return None;
        }
    } else {
        log::error(
            "hts_parse_region",
            format_args!(
                "Unexpected string \"{}\" after region",
                String::from_utf8_lossy(&s[hyphen..])
            ),
        );
        return None;
    }
    if end == 0 {
        end = POS_MAX;
    }
    if beg < 0 {
        beg = 0;
    }
    if beg >= end {
        return None;
    }
    Some(ParsedRegion {
        tid,
        beg,
        end,
        consumed: s_end,
    })
}

/// Parses the interval of a region without resolving the name (HTSlib's `hts_parse_reg64`):
/// returns the 0-based start, exclusive end, and the position of the colon (the name's end),
/// or `None` for a malformed or empty interval.
#[must_use]
pub fn parse_reg(s: &[u8]) -> Option<(i64, i64, usize)> {
    let Some(colon) = s.iter().rposition(|&c| c == b':') else {
        return Some((0, POS_MAX, s.len()));
    };
    let flags = ParseDecimalFlags {
        thousands_separator: true,
    };
    let (b, hyphen) = parse_decimal(s, colon + 1, flags);
    let beg = (b - 1).max(0);
    let end = match s.get(hyphen) {
        None => POS_MAX,
        Some(b'-') => parse_decimal(s, hyphen + 1, flags).0,
        Some(_) => return None,
    };
    if beg >= end {
        return None;
    }
    Some((beg, end, colon))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lookup(name: &[u8]) -> Option<u32> {
        match name {
            b"chr1" => Some(0),
            b"chr1:100-200" => Some(1),
            b"HLA-A*01:01" => Some(2),
            _ => None,
        }
    }

    fn parse(s: &str, flags: RegionFlags) -> Option<(u32, i64, i64)> {
        parse_region(s.as_bytes(), lookup, flags).map(|r| (r.tid, r.beg, r.end))
    }

    #[test]
    fn samtools_and_bcftools_conventions() {
        let sam = RegionFlags::default();
        let bcf = RegionFlags {
            one_coord: true,
            ..RegionFlags::default()
        };
        assert_eq!(parse("chr1", sam), Some((0, 0, POS_MAX)));
        assert_eq!(parse("chr1:1,000-2,000", sam), Some((0, 999, 2000)));
        assert_eq!(parse("chr1:100", sam), Some((0, 99, POS_MAX)));
        assert_eq!(parse("chr1:100", bcf), Some((0, 99, 100)));
        assert_eq!(parse("chr1:-100", bcf), Some((0, 0, 100)));
        assert_eq!(parse("chr1:100-", bcf), Some((0, 99, POS_MAX)));
        assert_eq!(parse("{HLA-A*01:01}:5-6", sam), Some((2, 4, 6)));
        assert_eq!(parse("chr1:100-200", sam), None, "ambiguous");
        assert_eq!(parse("chrX:1-2", sam), None);
        assert_eq!(parse("chr1:0-5", sam), Some((0, 0, 5)));
        assert_eq!(parse("chr1:0-0", sam), Some((0, 0, POS_MAX)));
        assert_eq!(parse("chr1:-5-10", sam), None);
    }

    #[test]
    fn list_mode_stops_at_commas() {
        let flags = RegionFlags {
            list: true,
            ..RegionFlags::default()
        };
        let r = parse_region(b"chr1:5-6,chr1", lookup, flags).expect("the region parses");
        assert_eq!((r.beg, r.end, r.consumed), (4, 6, 9));
    }
}
