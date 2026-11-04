//! Whether `bgzip` treats its input as text, whose blocks then end at line ends: the part of
//! HTSlib's `hts_detect_format` that tells uncompressed text formats from everything else.

/// How much of the input `hts_detect_format` looks at.
pub const PEEK_LEN: usize = 1024;

/// Whether `bgzip` compresses an input starting with `s` (its first [`PEEK_LEN`] bytes, or all
/// of it if shorter) in text mode: the input is uncompressed and detected as plain text, SAM,
/// VCF, BED, FASTA, or FASTQ.
pub fn is_textual(s: &[u8]) -> bool {
    if is_compressed(s) || s.is_empty() {
        return false;
    }

    // `hts_detect_format` tests these in turn, and stops at the first format it recognizes. No
    // input matches both a binary and a text test, so they can be grouped.
    let is_binary_format =
        // CRAM.
        (s.len() >= 6 && s.starts_with(b"CRAM") && (1..=7).contains(&s[4]) && s[5] <= 7)
        // BAM, BAI, BCF, CSI, or TBI; anything else with this control byte is not text either.
        || (s.len() >= 4 && s[3] <= 4)
        || (s.len() >= 8 && s.starts_with(b"d4\xdd\xdd"))
        || starts_with_nonblank(s, b"{\"htsget\":")
        || (s.len() > 8 && s.starts_with(b"crypt4gh"));

    if is_binary_format {
        return false;
    }

    let is_text_format = s.starts_with(b"##fileformat=VCF")
        || [b"@HD\t", b"@SQ\t", b"@RG\t", b"@PG\t", b"@CO\t"]
            .iter()
            .any(|tag| s.starts_with(*tag))
        || ((s[0] == b'>' || s[0] == b'@') && is_fastaq(s))
        || is_tabbed_format(s);

    is_text_format || is_text_only(s)
}

/// Whether the first line looks like SAM or BED.
fn is_tabbed_format(s: &[u8]) -> bool {
    const SAM: &[u8] = b"ZiZiiCZiiZZOOOOOOOOOOOOOOOOOOOO+";

    parse_tabbed_text(s).is_some_and(|(columns, complete)| {
        !columns.is_empty()
            && (colmatch(&columns, SAM) >= 9 + 2 * usize::from(complete)
                || colmatch(&columns, b"Zii+") >= 3)
    })
}

/// Whether `s` starts with a gzip, bzip2, xz, or zstd magic number.
fn is_compressed(s: &[u8]) -> bool {
    s.starts_with(&[0x1f, 0x8b])
        || (s.len() >= 10
            && s.starts_with(b"BZh")
            && (s[4..10] == *b"\x31\x41\x59\x26\x53\x59"
                || s[4..10] == *b"\x17\x72\x45\x38\x50\x90"))
        || s.starts_with(b"\xfd7zXZ\0")
        || s.starts_with(b"\x28\xb5\x2f\xfd")
}

/// `cmp_nonblank(key, s) == 0`: `s` starts with `key`, ignoring whitespace in `s`.
fn starts_with_nonblank(s: &[u8], key: &[u8]) -> bool {
    let mut bytes = s.iter().filter(|&&b| !is_space(b));
    key.iter().all(|k| bytes.next() == Some(k))
}

/// C's `isspace`.
fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn is_text_only(s: &[u8]) -> bool {
    s.iter()
        .all(|&b| b >= b' ' || b == b'\t' || b == b'\r' || b == b'\n')
}

/// Whether `s` looks like FASTA or FASTQ: a textual first line, then bases up to the end of the
/// second line (or of `s`).
fn is_fastaq(s: &[u8]) -> bool {
    let Some(eol) = memchr::memchr(b'\n', s) else {
        // A very long first line.
        return is_text_only(s);
    };

    if !is_text_only(&s[..eol]) {
        return false;
    }

    let rest = &s[eol + 1..];
    let n = rest
        .iter()
        .take_while(|&&b| SEQ_NT16_TABLE[usize::from(b)] != 15 || b.eq_ignore_ascii_case(&b'N'))
        .count();

    if rest[..n].contains(&b'=') {
        return false;
    }

    matches!(rest.get(n), None | Some(b'\r' | b'\n'))
}

/// HTSlib's `seq_nt16_table` (15 for bytes that are not bases).
const SEQ_NT16_TABLE: [u8; 256] = {
    let mut table = [15; 256];
    let codes = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;

    while i < codes.len() {
        table[codes[i] as usize] = i as u8;
        table[codes[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }

    // N is 15 like non-bases (callers test for it), U is T, and '0'-'3' are color-space codes.
    table[b'U' as usize] = 8;
    table[b'u' as usize] = 8;
    table[b'0' as usize] = 1;
    table[b'1' as usize] = 2;
    table[b'2' as usize] = 4;
    table[b'3' as usize] = 8;
    table
};

/// `parse_tabbed_text`: the types of the first line's tab-separated columns (`i` integer, `s`
/// strand, `C` CIGAR, `O` SAM optional field, `Z` anything), and whether the line ended. `None`
/// for a control (or, as C `char` is signed, a non-ASCII) byte.
fn parse_tabbed_text(s: &[u8]) -> Option<(Vec<u8>, bool)> {
    const DIGIT: u8 = 1;
    const LEADING_SIGN: u8 = 2;
    const CIGAR_OPERATOR: u8 = 4;
    const OTHER: u8 = 8;
    const COLUMN_LEN: usize = 24;
    const CIGAR_OPS: &[u8] = b"MIDNSHP=XB";

    let mut columns = Vec::new();
    let mut start = 0;
    let mut seen = 0;

    for (i, &b) in s.iter().enumerate() {
        if (b' '..0x80).contains(&b) {
            seen |= if b.is_ascii_digit() {
                DIGIT
            } else if (b == b'+' || b == b'-') && i == start {
                LEADING_SIGN
            } else if CIGAR_OPS.contains(&b) && i > start && s[i - 1].is_ascii_digit() {
                CIGAR_OPERATOR
            } else {
                OTHER
            };
        } else if matches!(b, b'\t' | b'\r' | b'\n') {
            let field = &s[start..i];

            let column = if seen == DIGIT || seen == LEADING_SIGN | DIGIT {
                b'i'
            } else if seen == DIGIT | CIGAR_OPERATOR {
                b'C'
            } else if field.len() == 1 {
                match field[0] {
                    b'*' => b'C',
                    b'+' | b'-' | b'.' => b's',
                    _ => b'Z',
                }
            } else if field.len() >= 5 && field[2] == b':' && field[4] == b':' {
                b'O'
            } else {
                b'Z'
            };

            columns.push(column);

            if b != b'\t' || columns.len() >= COLUMN_LEN - 1 {
                return Some((columns, true));
            }

            start = i + 1;
            seen = 0;
        } else {
            return None;
        }
    }

    Some((columns, false))
}

/// `colmatch`: how many of `columns` match `pattern` (`Z` matches anything; `+` ends the match),
/// or 0 on a mismatch.
fn colmatch(columns: &[u8], pattern: &[u8]) -> usize {
    for (i, &c) in columns.iter().enumerate() {
        match pattern.get(i) {
            Some(b'+') => return i,
            Some(&p) if c == p || p == b'Z' => {}
            _ => return 0,
        }
    }

    columns.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_textual() {
        // Text formats.
        assert!(is_textual(b"1\n2\n3\n"));
        assert!(is_textual(b"##fileformat=VCFv4.2\n#CHROM\n"));
        assert!(is_textual(b"@HD\tVN:1.6\n"));
        assert!(is_textual(b">chr1\nACGTN\n"));
        assert!(is_textual(b"@r1\nACGT\n+\nIIII\n"));
        assert!(is_textual(b"chr1\t10\t20\tname\n"));
        assert!(is_textual(b"caf\xc3\xa9\n"));

        // Everything else.
        assert!(!is_textual(b""));
        assert!(!is_textual(b"\x1f\x8b\x08\x04"));
        assert!(!is_textual(b"BAM\x01\0\0\0\0"));
        assert!(!is_textual(b"CRAM\x03\x01xxxx"));
        assert!(!is_textual(b" { \"htsget\" : {}}"));
        assert!(!is_textual(b"\x28\xb5\x2f\xfd"));
        assert!(!is_textual(b"text\0with a NUL"));
    }

    #[test]
    fn test_is_fastaq() {
        assert!(is_fastaq(b">x\nACGT"));
        assert!(is_fastaq(b">x\nacgtn\r\n"));
        assert!(!is_fastaq(b">x\nACGT=\n"));
        assert!(!is_fastaq(b">x\nhello\n"));
        assert!(is_fastaq(b">a long first line without an end"));
    }

    #[test]
    fn test_parse_tabbed_text() {
        assert_eq!(
            parse_tabbed_text(b"r1\t0\tchr1\t+5\t60\t10M\t*\t-\tNM:i:0\n"),
            Some((b"ZiZiiCCsO".to_vec(), true))
        );
        assert_eq!(parse_tabbed_text(b"a\tb"), Some((b"Z".to_vec(), false)));
        assert_eq!(parse_tabbed_text(b"a\x01\n"), None);
        assert_eq!(parse_tabbed_text(b"\xc3\xa9\n"), None);
    }

    #[test]
    fn test_colmatch() {
        assert_eq!(colmatch(b"Zii", b"Zii+"), 3);
        assert_eq!(colmatch(b"ZiiZ", b"Zii+"), 3);
        assert_eq!(colmatch(b"Zi", b"Zii+"), 2);
        assert_eq!(colmatch(b"ii", b"Zii+"), 2);
        assert_eq!(colmatch(b"Zs", b"Zii+"), 0);
    }
}
