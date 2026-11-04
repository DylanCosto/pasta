use std::{
    collections::HashMap,
    io::{self, Write},
};

use super::Type;
use crate::cram::io::writer::num::{write_u8, write_u32_le, write_uint7};

const NUL: u8 = 0x00;

pub fn encode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    if let Some(buf) = src.strip_suffix(&[NUL]) {
        src = buf;
    }

    let names: Vec<_> = src.split(|&b| b == NUL).collect();

    write_header(&mut dst, src.len(), names.len())?;

    let mut names_indices = HashMap::new();
    let mut diffs = Vec::with_capacity(names.len());
    let mut max_token_count = 0;

    if let Some(name) = names.first() {
        let diff = build_first_diff(name);
        max_token_count = max_token_count.max(diff.tokens.len());
        diffs.push(diff);
    }

    for (i, name) in names.iter().enumerate().skip(1) {
        let diff = build_diff(&diffs, &names_indices, i, name);
        names_indices.entry(name).or_insert(i);
        max_token_count = max_token_count.max(diff.tokens.len());
        diffs.push(diff);
    }

    let mut token_writer = TokenWriter::default();

    for diff in &diffs {
        let token = match diff.mode {
            Mode::Dup(delta) => Token::Dup(delta),
            Mode::Diff(delta) => Token::Diff(delta),
        };

        token_writer.write_token(&token)?;
    }

    encode_token_byte_streams(&mut dst, &token_writer)?;

    for i in 0..max_token_count {
        let mut token_writer = TokenWriter::default();

        for diff in &diffs {
            if diff.is_dup() {
                continue;
            }

            if let Some(token) = diff.tokens.get(i) {
                token_writer.write_token(token)?;
            }
        }

        encode_token_byte_streams(&mut dst, &token_writer)?;
    }

    Ok(dst)
}

fn write_header<W>(writer: &mut W, src_len: usize, names_count: usize) -> io::Result<()>
where
    W: Write,
{
    let ulen =
        u32::try_from(src_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, ulen)?;

    let n_names =
        u32::try_from(names_count).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, n_names)?;

    // TODO: use_arith
    write_u8(writer, 0)?;

    Ok(())
}

#[derive(Debug)]
enum Token {
    String(Vec<u8>),
    Char(u8),
    PaddedDigits(u32, usize),
    // ...
    Dup(usize),
    Diff(usize),
    Digits(u32),
    Delta(u8),
    Delta0(u8),
    Match,
    // ...
    End,
}

impl Token {
    fn ty(&self) -> Type {
        match self {
            Self::String(_) => Type::String,
            Self::Char(_) => Type::Char,
            Self::PaddedDigits(..) => Type::Digits0,
            // ...
            Self::Dup(_) => Type::Dup,
            Self::Diff(_) => Type::Diff,
            Self::Digits(_) => Type::Digits,
            Self::Delta(..) => Type::Delta,
            Self::Delta0(..) => Type::Delta0,
            Self::Match => Type::Match,
            // ...
            Self::End => Type::End,
        }
    }
}

enum Mode {
    Diff(usize),
    Dup(usize),
}

struct Diff {
    mode: Mode,
    raw_tokens: Vec<Vec<u8>>,
    tokens: Vec<Token>,
    /// The tokens as the decoder resolves them (a match or a delta becomes the value it
    /// stands for), which later names are encoded against.
    resolved: Vec<Resolved>,
}

/// A token's value after decoding.
#[derive(Clone, Copy, Debug)]
enum Resolved {
    Digits(u32),
    PaddedDigits(u32, usize),
    Other,
}

impl Diff {
    fn new(mode: Mode) -> Self {
        Self {
            mode,
            raw_tokens: Vec::new(),
            tokens: Vec::new(),
            resolved: Vec::new(),
        }
    }

    fn is_dup(&self) -> bool {
        matches!(self.mode, Mode::Dup(_))
    }
}

/// Splits a name into runs of letters, runs of digits, and single other characters, as
/// htscodecs' tokenizer does, so that numbers can be stored as deltas from the previous name.
fn tokenize(b: &[u8]) -> impl Iterator<Item = &[u8]> {
    use std::iter;

    let mut start = 0;

    iter::from_fn(move || {
        let first = *b.get(start)?;

        let len = if first.is_ascii_alphabetic() {
            b[start..]
                .iter()
                .take_while(|c| c.is_ascii_alphabetic())
                .count()
        } else if first.is_ascii_digit() {
            b[start..].iter().take_while(|c| c.is_ascii_digit()).count()
        } else {
            1
        };

        let token = &b[start..start + len];
        start += len;
        Some(token)
    })
}

/// A token for `raw` on its own: zero-padded digits, digits, a character, or a string.
fn new_token(raw: &[u8]) -> (Token, Resolved) {
    if let Some(n) = parse_digits0(raw) {
        (
            Token::PaddedDigits(n, raw.len()),
            Resolved::PaddedDigits(n, raw.len()),
        )
    } else if let Some(n) = parse_digits(raw) {
        (Token::Digits(n), Resolved::Digits(n))
    } else if let [c] = raw {
        (Token::Char(*c), Resolved::Other)
    } else {
        (Token::String(raw.into()), Resolved::Other)
    }
}

/// A token for `raw` given the previous name's token in the same position: a match, a delta
/// from its number, or a new token.
fn token_after(raw: &[u8], prev_raw: &[u8], prev: Resolved) -> (Token, Resolved) {
    if raw == prev_raw {
        return (Token::Match, prev);
    }

    match prev {
        // Printed without padding, so the number must have no leading zero (a `0` is padded
        // digits).
        Resolved::Digits(n) if raw.first() != Some(&b'0') => {
            if let Some(m) = parse_digits(raw)
                && let Some(delta) = small_delta(n, m)
            {
                return (Token::Delta(delta), Resolved::Digits(m));
            }
        }
        // Printed in exactly the previous width (htscodecs keeps it fixed), which must give back
        // `raw`: a number that grows wider starts a new token.
        Resolved::PaddedDigits(n, width) => {
            if let Some(m) = parse_digits(raw)
                && let Some(delta) = small_delta(n, m)
                && raw.len() == width
                && digit_count(m) <= width
            {
                return (Token::Delta0(delta), Resolved::PaddedDigits(m, width));
            }
        }
        _ => {}
    }

    new_token(raw)
}

fn small_delta(n: u32, m: u32) -> Option<u8> {
    m.checked_sub(n).and_then(|delta| u8::try_from(delta).ok())
}

fn digit_count(n: u32) -> usize {
    n.checked_ilog10().map_or(1, |log| log as usize + 1)
}

fn build_first_diff(name: &[u8]) -> Diff {
    let mut diff = Diff::new(Mode::Diff(0));

    for raw_token in tokenize(name) {
        let (token, resolved) = new_token(raw_token);
        diff.raw_tokens.push(raw_token.into());
        diff.tokens.push(token);
        diff.resolved.push(resolved);
    }

    diff.tokens.push(Token::End);

    diff
}

fn build_diff(
    diffs: &[Diff],
    names_indices: &HashMap<&[u8], usize>,
    i: usize,
    name: &[u8],
) -> Diff {
    if let Some(j) = names_indices.get(name) {
        // A duplicate resolves to the name it copies.
        let mut diff = Diff::new(Mode::Dup(i - j));
        diff.raw_tokens = diffs[*j].raw_tokens.clone();
        diff.resolved = diffs[*j].resolved.clone();
        return diff;
    }

    let mut diff = Diff::new(Mode::Diff(1));
    let prev_diff = &diffs[i - 1];

    for (j, raw_token) in tokenize(name).enumerate() {
        let (token, resolved) = match (prev_diff.raw_tokens.get(j), prev_diff.resolved.get(j)) {
            (Some(prev_raw), Some(&prev)) => token_after(raw_token, prev_raw, prev),
            _ => new_token(raw_token),
        };

        diff.raw_tokens.push(raw_token.into());
        diff.tokens.push(token);
        diff.resolved.push(resolved);
    }

    diff.tokens.push(Token::End);

    diff
}

/// Zero-padded digits: digits with a leading zero, including a lone `0` (which htscodecs decodes
/// to nothing as plain digits).
fn parse_digits0(s: &[u8]) -> Option<u32> {
    if s.starts_with(b"0") {
        parse_u32(s).ok()
    } else {
        None
    }
}

fn parse_digits(s: &[u8]) -> Option<u32> {
    if s.iter().all(u8::is_ascii_digit) {
        parse_u32(s).ok()
    } else {
        None
    }
}

fn parse_u32(src: &[u8]) -> io::Result<u32> {
    lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[derive(Default)]
struct TokenWriter {
    type_writer: Vec<u8>,
    string_writer: Vec<u8>,
    char_writer: Vec<u8>,
    digits0_writer: Vec<u8>,
    dz_len_writer: Vec<u8>,
    dup_writer: Vec<u8>,
    diff_writer: Vec<u8>,
    digits_writer: Vec<u8>,
    delta_writer: Vec<u8>,
    delta0_writer: Vec<u8>,
}

impl TokenWriter {
    fn write_token(&mut self, token: &Token) -> io::Result<()> {
        write_type(&mut self.type_writer, token.ty())?;

        match token {
            Token::String(s) => {
                self.string_writer.write_all(s)?;
                self.string_writer.write_all(&[NUL])?;
            }
            Token::Char(b) => write_u8(&mut self.char_writer, *b)?,
            Token::PaddedDigits(n, width) => {
                // d
                write_u32_le(&mut self.digits0_writer, *n)?;

                let l = u8::try_from(*width)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u8(&mut self.dz_len_writer, l)?;
            }
            // ...
            Token::Dup(delta) => {
                let n = u32::try_from(*delta)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u32_le(&mut self.dup_writer, n)?;
            }
            Token::Diff(delta) => {
                let n = u32::try_from(*delta)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u32_le(&mut self.diff_writer, n)?;
            }
            Token::Digits(n) => {
                write_u32_le(&mut self.digits_writer, *n)?;
            }
            Token::Delta(delta) => {
                write_u8(&mut self.delta_writer, *delta)?;
            }
            Token::Delta0(delta) => {
                write_u8(&mut self.delta0_writer, *delta)?;
            }
            Token::Match => {}
            // ...
            Token::End => {}
        }

        Ok(())
    }
}

fn write_type<W>(writer: &mut W, ty: Type) -> io::Result<()>
where
    W: Write,
{
    let n = match ty {
        Type::Type => 0,
        Type::String => 1,
        Type::Char => 2,
        Type::Digits0 => 3,
        Type::DZLen => 4,
        Type::Dup => 5,
        Type::Diff => 6,
        Type::Digits => 7,
        Type::Delta => 8,
        Type::Delta0 => 9,
        Type::Match => 10,
        Type::Nop => 11,
        Type::End => 12,
    };

    write_u8(writer, n)
}

fn encode_token_byte_streams<W>(writer: &mut W, token_writer: &TokenWriter) -> io::Result<()>
where
    W: Write,
{
    encode_token_byte_stream(writer, Type::Type, &token_writer.type_writer)?;
    encode_token_byte_stream(writer, Type::String, &token_writer.string_writer)?;
    encode_token_byte_stream(writer, Type::Char, &token_writer.char_writer)?;
    encode_token_byte_stream(writer, Type::Digits0, &token_writer.digits0_writer)?;
    encode_token_byte_stream(writer, Type::DZLen, &token_writer.dz_len_writer)?;
    encode_token_byte_stream(writer, Type::Dup, &token_writer.dup_writer)?;
    encode_token_byte_stream(writer, Type::Diff, &token_writer.diff_writer)?;
    encode_token_byte_stream(writer, Type::Digits, &token_writer.digits_writer)?;
    encode_token_byte_stream(writer, Type::Delta, &token_writer.delta_writer)?;
    encode_token_byte_stream(writer, Type::Delta0, &token_writer.delta0_writer)?;

    Ok(())
}

fn encode_token_byte_stream<W>(writer: &mut W, ty: Type, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    use crate::cram::codecs::rans_nx16::{self, Flags};

    if buf.is_empty() {
        return Ok(());
    }

    if ty == Type::Type {
        write_u8(writer, 0x80)?;
    } else {
        write_type(writer, ty)?;
    }

    // As htscodecs does, try a few rANS variants and keep the smallest.
    const FLAGS: [u8; 8] = [0, 1, 64, 65, 128, 129, 192, 193];

    let mut cdata: Option<Vec<u8>> = None;

    for flags in FLAGS {
        let data = rans_nx16::encode(Flags::from_bits_truncate(flags), buf)?;

        if cdata.as_ref().is_none_or(|best| data.len() < best.len()) {
            cdata = Some(data);
        }
    }

    let cdata = cdata.expect("at least one variant is tried");

    let clen =
        u32::try_from(cdata.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(writer, clen)?;

    writer.write_all(&cdata)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_round_trip_and_size() -> std::io::Result<()> {
        let mut names: Vec<String> = Vec::new();

        for i in 0..5000 {
            names.push(format!("short{:09}", i / 2));
        }

        for i in 0..3000u32 {
            names.push(format!(
                "A00123:8:H5KJXDSXY:1:{}:{}:{}",
                1101 + i / 1000,
                (i * 37) % 32000,
                1000 + i % 7
            ));
        }

        // Numbers that change width, go down, lose padding, or are too large for 32 bits.
        for name in [
            "r099",
            "r100",
            "r101",
            "r0",
            "r1",
            "r00",
            "r01",
            "x9",
            "x10",
            "x08",
            "x09",
            "x10",
            "y5000",
            "y4999",
            "z4294967295",
            "z4294967296",
            "z18446744073709551616",
            "",
            "a",
            "-",
            "a-b.c_d",
            "q01q",
            "q02q",
            "q02q",
        ] {
            names.push(name.into());
        }

        let mut src = Vec::new();

        for name in &names {
            src.extend_from_slice(name.as_bytes());
            src.push(0);
        }

        let encoded = super::encode(&src)?;
        let decoded = super::super::decode(&encoded)?;
        assert_eq!(decoded, src);

        // Sequential names cost next to nothing.
        let mut src = Vec::new();

        for i in 0..10_000 {
            src.extend_from_slice(format!("short{:09}", i / 2).as_bytes());
            src.push(0);
        }

        assert!(super::encode(&src)?.len() < 500);

        Ok(())
    }
}
