//! Splitting an expression into tokens (bcftools' `filters_next_token`).

use pasta::variant::number;

use super::token::{Fraction, Op};

/// What [`next_token`] found.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Lexed {
    /// An operator, bracket, or function with one argument (the position moved past it).
    Op(Op),
    /// A function with a variable number of comma-separated arguments (`binom`, `fisher`,
    /// `perl.NAME`); the position is on the opening bracket.
    VarArgs(Op),
    /// `N_PASS(` or `F_PASS(`; the position is on the opening bracket.
    Pass(Fraction),
    /// A value of `len` bytes at the position (0 at the end or on an unparseable character).
    Val(usize),
    /// An external value (`{}`, `{str}`, ...) of `len` bytes.
    Ext(usize),
    /// An unterminated quote.
    MissingQuotes,
}

#[inline]
fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

#[inline]
fn starts_with_ci(s: &[u8], i: usize, prefix: &[u8]) -> bool {
    s.len() >= i + prefix.len() && s[i..i + prefix.len()].eq_ignore_ascii_case(prefix)
}

/// Functions recognised by name, with how far to advance (onto the opening bracket).
const FUNCTIONS: &[(&[u8], usize, Op)] = &[
    (b"SMPL_COUNT(", 10, Op::SCount),
    (b"SMPL_MAX(", 8, Op::SMax),
    (b"SMPL_MIN(", 8, Op::SMin),
    (b"SMPL_MEAN(", 9, Op::SAvg),
    (b"SMPL_MEDIAN(", 11, Op::SMedian),
    (b"SMPL_AVG(", 8, Op::SAvg),
    (b"SMPL_STDEV(", 10, Op::SStdev),
    (b"SMPL_SUM(", 8, Op::SSum),
    (b"sCOUNT(", 6, Op::SCount),
    (b"sMAX(", 4, Op::SMax),
    (b"sMIN(", 4, Op::SMin),
    (b"sMEAN(", 5, Op::SAvg),
    (b"sMEDIAN(", 7, Op::SMedian),
    (b"sAVG(", 4, Op::SAvg),
    (b"sSTDEV(", 6, Op::SStdev),
    (b"sSUM(", 4, Op::SSum),
    (b"MAX(", 3, Op::Max),
    (b"MIN(", 3, Op::Min),
    (b"MEAN(", 4, Op::Avg),
    (b"MEDIAN(", 6, Op::Median),
    (b"AVG(", 3, Op::Avg),
    (b"STDEV(", 5, Op::Stdev),
    (b"SUM(", 3, Op::Sum),
    (b"ABS(", 3, Op::Abs),
    (b"COUNT(", 5, Op::Cnt),
    (b"STRLEN(", 6, Op::Len),
    (b"BINOM(", 5, Op::Binom),
    // bcftools compares only "FISHER" here, then steps onto the bracket.
    (b"FISHER", 6, Op::Fisher),
    (b"PHRED(", 5, Op::Phred),
    (b"%MAX(", 4, Op::Max),
    (b"%MIN(", 4, Op::Min),
    (b"%AVG(", 4, Op::Avg),
    (b"%SUM(", 4, Op::Sum),
];

/// Reads the next token of `s` at `*pos`, skipping leading whitespace. Operators and functions
/// advance `*pos`; values do not (the caller consumes their length).
pub(super) fn next_token(s: &[u8], pos: &mut usize) -> Lexed {
    let mut tmp = *pos;
    while at(s, tmp).is_ascii_whitespace() {
        tmp += 1;
    }
    *pos = tmp;
    let start = tmp;

    // A number: d.ddde[+-]dd (strtod would eat a sign).
    if at(s, start).is_ascii_digit() || at(s, start) == b'.' {
        let (_, end) = number::strtod(s, start);
        if end != start && !at(s, end).is_ascii_alphanumeric() {
            return Lexed::Val(end - start);
        }
    }

    for &(name, advance, op) in FUNCTIONS {
        if starts_with_ci(s, tmp, name) {
            *pos += advance;
            return match op {
                Op::Binom | Op::Fisher => Lexed::VarArgs(op),
                _ => Lexed::Op(op),
            };
        }
    }
    if starts_with_ci(s, tmp, b"INFO/") {
        tmp += 5;
    }
    if starts_with_ci(s, tmp, b"FORMAT/") {
        tmp += 7;
    }
    if starts_with_ci(s, tmp, b"FMT/") {
        tmp += 4;
    }
    if starts_with_ci(s, tmp, b"PERL.") {
        *pos += 5;
        return Lexed::VarArgs(Op::PerlSub);
    }
    if starts_with_ci(s, tmp, b"N_PASS(") {
        *pos += 6;
        return Lexed::Pass(Fraction::Count);
    }
    if starts_with_ci(s, tmp, b"F_PASS(") {
        *pos += 6;
        return Lexed::Pass(Fraction::Fraction);
    }
    if starts_with_ci(s, tmp, b"%ILEN") {
        return Lexed::Val(5);
    }
    for ext in [&b"{}"[..], b"{STR}", b"{INT}", b"{FLOAT}"] {
        if starts_with_ci(s, tmp, ext) {
            return Lexed::Ext(ext.len());
        }
    }

    if at(s, tmp) == b'@' {
        // A file name.
        while !matches!(at(s, tmp), 0 | b'=' | b'!') && !at(s, tmp).is_ascii_whitespace() {
            tmp += 1;
        }
        return Lexed::Val(tmp - start);
    }

    let mut square_brackets = 0;
    while tmp < s.len() {
        let c = s[tmp];
        if square_brackets == 0
            && (matches!(
                c,
                b'"' | b'\''
                    | b'<'
                    | b'>'
                    | b'='
                    | b'!'
                    | b'&'
                    | b'|'
                    | b'('
                    | b')'
                    | b'+'
                    | b'*'
                    | b'-'
                    | b'/'
                    | b'~'
                    | b'%'
            ) || c.is_ascii_whitespace())
        {
            break;
        }
        if c == b']' {
            if square_brackets > 0 {
                tmp += 1;
            }
            break;
        }
        if c == b'[' {
            square_brackets += 1;
        }
        tmp += 1;
    }
    if tmp > start {
        return Lexed::Val(tmp - start);
    }

    let c = at(s, tmp);
    if c == b'"' || c == b'\'' {
        tmp += 1;
        while tmp < s.len() && s[tmp] != c {
            tmp += 1;
        }
        if tmp >= s.len() {
            return Lexed::MissingQuotes;
        }
        return Lexed::Val(tmp - start + 1);
    }
    let next = at(s, tmp + 1);
    let (op, len) = match (c, next) {
        (b'!', b'=') => (Op::Ne, 2),
        (b'!', b'~') => (Op::NLike, 2),
        (b'<', b'=') => (Op::Le, 2),
        (b'<', _) => (Op::Lt, 1),
        (b'>', b'=') => (Op::Be, 2),
        (b'>', _) => (Op::Bt, 1),
        (b'=', b'=') => (Op::Eq, 2),
        (b'=', _) => (Op::Eq, 1),
        (b'(', _) => (Op::Lft, 1),
        (b')', _) => (Op::Rgt, 1),
        (b'&', b'&') => (Op::AndVec, 2),
        (b'|', b'|') => (Op::OrVec, 2),
        (b'&', _) => (Op::And, 1),
        (b'|', _) => (Op::Or, 1),
        (b'+', _) => (Op::Add, 1),
        (b'-', _) => (Op::Sub, 1),
        (b'*', _) => (Op::Mult, 1),
        (b'/', _) => (Op::Div, 1),
        (b'~', _) => (Op::Like, 1),
        (b'%', _) => (Op::Modulo, 1),
        _ => return Lexed::Val(0),
    };
    *pos += len;
    Lexed::Op(op)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokens(s: &str) -> Vec<String> {
        let s = s.as_bytes();
        let mut pos = 0;
        let mut out = Vec::new();
        loop {
            match next_token(s, &mut pos) {
                Lexed::Val(0) => break,
                Lexed::Val(n) | Lexed::Ext(n) => {
                    out.push(String::from_utf8_lossy(&s[pos..pos + n]).into_owned());
                    pos += n;
                }
                Lexed::MissingQuotes => {
                    out.push("MissingQuotes".into());
                    break;
                }
                other => out.push(format!("{other:?}")),
            }
        }
        out
    }

    #[test]
    fn splits_expressions() {
        assert_eq!(
            tokens("QUAL>=30 && FMT/AD[0:1]>2"),
            [
                "QUAL",
                "Op(Be)",
                "30",
                "Op(AndVec)",
                "FMT/AD[0:1]",
                "Op(Bt)",
                "2"
            ]
        );
        assert_eq!(
            tokens("MAX(DP)<1e3"),
            ["Op(Max)", "Op(Lft)", "DP", "Op(Rgt)", "Op(Lt)", "1e3"]
        );
        assert_eq!(
            tokens("ID=@a.txt | 'x y'"),
            ["ID", "Op(Eq)", "@a.txt", "Op(Or)", "'x y'"]
        );
        assert_eq!(tokens("1a"), ["1a"]);
        assert_eq!(tokens("\"abc"), ["MissingQuotes"]);
    }
}
