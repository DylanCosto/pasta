//! Parsing of samtools-style flag arguments.

use anyhow::bail;
use pasta::sam::alignment::record::Flags;

/// Flag names accepted by samtools, in bit order.
const NAMES: [(&str, Flags); 12] = [
    ("PAIRED", Flags::SEGMENTED),
    ("PROPER_PAIR", Flags::PROPERLY_SEGMENTED),
    ("UNMAP", Flags::UNMAPPED),
    ("MUNMAP", Flags::MATE_UNMAPPED),
    ("REVERSE", Flags::REVERSE_COMPLEMENTED),
    ("MREVERSE", Flags::MATE_REVERSE_COMPLEMENTED),
    ("READ1", Flags::FIRST_SEGMENT),
    ("READ2", Flags::LAST_SEGMENT),
    ("SECONDARY", Flags::SECONDARY),
    ("QCFAIL", Flags::QC_FAIL),
    ("DUP", Flags::DUPLICATE),
    ("SUPPLEMENTARY", Flags::SUPPLEMENTARY),
];

/// Parses a flag set as HTSlib's `bam_str2flag` does: an integer (decimal, `0x` hexadecimal or
/// `0` octal, with an optional sign) or else a comma-separated list of names (`UNMAP,SECONDARY`).
/// Like `strtol`, an integer is read from the start of the string and the rest is ignored, so
/// `0x80,0x10` is 0x80.
pub fn parse(s: &str) -> anyhow::Result<Flags> {
    let s = s.trim();

    if let Some(n) = parse_integer_prefix(s) {
        // samtools keeps flags in an `int` and tests them against 16-bit record flags.
        return Ok(Flags::from_bits_retain(n as u16));
    }

    let mut flags = Flags::empty();

    for name in s.split(',') {
        let name = name.trim();

        match NAMES.iter().find(|(n, _)| n.eq_ignore_ascii_case(name)) {
            Some((_, flag)) => flags |= *flag,
            None => bail!("unknown flag name: {name}"),
        }
    }

    Ok(flags)
}

/// The union of the values of a flag option given several times, as samtools ORs them, or `None`
/// if it was not given.
pub fn union(values: &[Flags]) -> Option<Flags> {
    values.iter().copied().reduce(|a, b| a | b)
}

/// Reads an integer from the start of `s` as `strtol(s, &end, 0)` does, or returns `None` if it
/// starts with no digits. Values beyond `i64` saturate, as `strtol` clamps to `LONG_MAX`.
fn parse_integer_prefix(s: &str) -> Option<i64> {
    let (negative, rest) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };

    let bytes = rest.as_bytes();

    let (radix, digits) = match bytes {
        [b'0', b'x' | b'X', d, ..] if d.is_ascii_hexdigit() => (16, &bytes[2..]),
        [b'0', ..] => (8, bytes),
        _ => (10, bytes),
    };

    let len = digits
        .iter()
        .take_while(|&&b| char::from(b).is_digit(radix))
        .count();

    if len == 0 {
        return None;
    }

    let magnitude = digits[..len].iter().fold(0i64, |n, &b| {
        let d = char::from(b).to_digit(radix).expect("a digit of the radix");
        n.saturating_mul(i64::from(radix))
            .saturating_add(i64::from(d))
    });

    Some(if negative {
        magnitude.saturating_neg()
    } else {
        magnitude
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> anyhow::Result<()> {
        assert_eq!(parse("4")?, Flags::UNMAPPED);
        assert_eq!(parse("0x904")?, Flags::from_bits_retain(0x904));
        assert_eq!(parse("010")?, Flags::MATE_UNMAPPED);
        assert_eq!(parse("0")?, Flags::empty());
        assert_eq!(
            parse("UNMAP,secondary")?,
            Flags::UNMAPPED | Flags::SECONDARY
        );
        assert!(parse("NOPE").is_err());
        assert!(parse("DUP,0x200").is_err());
        assert!(parse("UNMAP,").is_err());
        Ok(())
    }

    #[test]
    fn test_parse_reads_an_integer_prefix_as_strtol_does() -> anyhow::Result<()> {
        assert_eq!(parse("0x80,0x10")?, Flags::from_bits_retain(0x80));
        assert_eq!(parse("4,UNMAP")?, Flags::UNMAPPED);
        assert_eq!(parse("0xZZ")?, Flags::empty());
        assert_eq!(parse("09")?, Flags::empty());
        assert_eq!(parse("+12")?, Flags::from_bits_retain(12));
        assert_eq!(parse("-4")?, Flags::from_bits_retain(0xfffc));
        assert_eq!(parse("0x10000")?, Flags::empty());
        Ok(())
    }
}
