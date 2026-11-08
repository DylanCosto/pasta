use std::{
    collections::hash_map::Entry,
    io::{self, Write},
    ops::Range,
};

use rustc_hash::FxHashMap;

use super::Type;
use crate::cram::io::writer::num::{write_u8, write_u32_le, write_uint7};

const NUL: u8 = 0x00;

/// Compresses NUL-terminated read names with the name tokenizer, coding its token streams with
/// rANS Nx16 at a CRAM compression `level` (1–9), as HTSlib's `tok3_encode_names` does.
pub fn encode(mut src: &[u8], level: u8) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    if let Some(buf) = src.strip_suffix(&[NUL]) {
        src = buf;
    }

    let names: Vec<_> = src.split(|&b| b == NUL).collect();

    write_header(&mut dst, src.len(), names.len())?;

    let diffs = Diffs::build(&names);
    let mut token_writer = TokenWriter::default();

    for diff in &diffs.diffs {
        let token = match diff.mode {
            Mode::Dup(delta) => Token::Dup(delta),
            Mode::Diff(delta) => Token::Diff(delta),
        };

        token_writer.write_token(&token)?;
    }

    encode_token_byte_streams(&mut dst, &token_writer, level)?;

    for i in 0..diffs.max_token_count {
        let mut token_writer = TokenWriter::default();

        // A duplicate has no tokens.
        for diff in &diffs.diffs {
            if let Some(token) = diffs.tokens[diff.tokens.clone()].get(i) {
                token_writer.write_token(token)?;
            }
        }

        encode_token_byte_streams(&mut dst, &token_writer, level)?;
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
enum Token<'a> {
    String(&'a [u8]),
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

impl Token<'_> {
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

/// How a name is encoded.
struct Diff {
    mode: Mode,
    /// The name's tokens in [`Diffs::tokens`]; none for a duplicate.
    tokens: Range<usize>,
    /// The name's raw and resolved tokens in [`Diffs::raw`] and [`Diffs::resolved`], which the next
    /// name is encoded against: its own, or for a duplicate those of the name it copies.
    view: Range<usize>,
}

/// The diffs of a list of names, with their tokens in flat lists that the diffs index.
struct Diffs<'a> {
    diffs: Vec<Diff>,
    tokens: Vec<Token<'a>>,
    raw: Vec<&'a [u8]>,
    /// The tokens as the decoder resolves them (a match or a delta becomes the value it stands
    /// for).
    resolved: Vec<Resolved>,
    /// The most tokens of any name, counting its end.
    max_token_count: usize,
}

/// A token's value after decoding.
#[derive(Clone, Copy, Debug)]
enum Resolved {
    Digits(u32),
    PaddedDigits(u32, usize),
    Other,
}

impl<'a> Diffs<'a> {
    fn build(names: &[&'a [u8]]) -> Self {
        let mut diffs = Self {
            diffs: Vec::with_capacity(names.len()),
            tokens: Vec::new(),
            raw: Vec::new(),
            resolved: Vec::new(),
            max_token_count: 0,
        };

        let mut indices: FxHashMap<&[u8], usize> = FxHashMap::default();

        for (i, &name) in names.iter().enumerate() {
            match indices.entry(name) {
                Entry::Occupied(entry) => {
                    // A duplicate resolves to the name it copies.
                    let j = *entry.get();
                    let view = diffs.diffs[j].view.clone();

                    diffs.diffs.push(Diff {
                        mode: Mode::Dup(i - j),
                        tokens: 0..0,
                        view,
                    });
                }
                Entry::Vacant(entry) => {
                    entry.insert(i);

                    // The first name stands alone; the others are encoded against the one before.
                    match i.checked_sub(1) {
                        Some(h) => {
                            let prev = diffs.diffs[h].view.clone();
                            diffs.push(Mode::Diff(1), name, prev);
                        }
                        None => diffs.push(Mode::Diff(0), name, 0..0),
                    }
                }
            }
        }

        diffs
    }

    /// Adds a name, encoding each token against the token in the same position of `prev`, a
    /// range of [`Self::raw`].
    fn push(&mut self, mode: Mode, name: &'a [u8], prev: Range<usize>) {
        let (tokens_start, view_start) = (self.tokens.len(), self.raw.len());

        for (j, raw) in tokenize(name).enumerate() {
            let k = prev.start + j;

            let (token, resolved) = if k < prev.end {
                token_after(raw, self.raw[k], self.resolved[k])
            } else {
                new_token(raw)
            };

            self.tokens.push(token);
            self.raw.push(raw);
            self.resolved.push(resolved);
        }

        self.tokens.push(Token::End);
        self.max_token_count = self.max_token_count.max(self.tokens.len() - tokens_start);

        self.diffs.push(Diff {
            mode,
            tokens: tokens_start..self.tokens.len(),
            view: view_start..self.raw.len(),
        });
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
fn new_token(raw: &[u8]) -> (Token<'_>, Resolved) {
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
        (Token::String(raw), Resolved::Other)
    }
}

/// A token for `raw` given the previous name's token in the same position: a match, a delta
/// from its number, or a new token.
fn token_after<'a>(raw: &'a [u8], prev_raw: &[u8], prev: Resolved) -> (Token<'a>, Resolved) {
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
    fn write_token(&mut self, token: &Token<'_>) -> io::Result<()> {
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

fn encode_token_byte_streams<W>(
    writer: &mut W,
    token_writer: &TokenWriter,
    level: u8,
) -> io::Result<()>
where
    W: Write,
{
    let streams = [
        (Type::Type, &token_writer.type_writer),
        (Type::String, &token_writer.string_writer),
        (Type::Char, &token_writer.char_writer),
        (Type::Digits0, &token_writer.digits0_writer),
        (Type::DZLen, &token_writer.dz_len_writer),
        (Type::Dup, &token_writer.dup_writer),
        (Type::Diff, &token_writer.diff_writer),
        (Type::Digits, &token_writer.digits_writer),
        (Type::Delta, &token_writer.delta_writer),
        (Type::Delta0, &token_writer.delta0_writer),
    ];

    for (ty, buf) in streams {
        encode_token_byte_stream(writer, ty, buf, level)?;
    }

    Ok(())
}

/// The rANS Nx16 flags htscodecs tries on a token stream of type `ty` (the table `R` of
/// `compress` in `tokenise_name3.c`). HTSlib caps the level at 3 when it codes names with rANS
/// (`cram_compress_by_method`), so only two of the table's rows apply: one for levels 1 and 2,
/// and one for 3 and above. The number streams hold 4-byte values, which striping (`0x08`)
/// splits into a stream per byte.
///
/// One difference remains: htscodecs codes each striped stream's chunks with whichever of order
/// 1, RLE and PACK the flags allow comes out smallest, where pasta's rANS encoder codes them with
/// order 0. On 2 million Illumina read names that search saved 0.1% of the CRAM file for 4.5% more
/// time writing it.
fn rans_flags(ty: Type, level: u8) -> &'static [u8] {
    if level < 3 {
        match ty {
            Type::Type => &[128],
            Type::String => &[129],
            Type::Digits0 | Type::Dup | Type::Diff | Type::Digits => &[8],
            Type::Delta0 => &[128],
            Type::Char | Type::DZLen | Type::Delta | Type::Match | Type::Nop | Type::End => &[0],
        }
    } else {
        match ty {
            Type::Type => &[192, 0],
            Type::String => &[129, 1],
            Type::Digits0 => &[136, 0],
            Type::Dup | Type::Digits => &[200],
            Type::Diff => &[136],
            Type::Delta0 => &[128],
            Type::Char | Type::DZLen | Type::Delta | Type::Match | Type::Nop | Type::End => &[0],
        }
    }
}

fn encode_token_byte_stream<W>(writer: &mut W, ty: Type, buf: &[u8], level: u8) -> io::Result<()>
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

    // Try htscodecs' flags for the stream and keep the smallest result. Striping needs a length
    // that divides by 4, which the number streams always have; plain order 0 stands in if no flags
    // are left.
    let flags = rans_flags(ty, level)
        .iter()
        .copied()
        .filter(|&flags| buf.len().is_multiple_of(4) || flags & Flags::STRIPE.bits() == 0);

    let mut cdata: Option<Vec<u8>> = None;

    for flags in flags {
        let data = rans_nx16::encode(Flags::from_bits_truncate(flags), buf)?;

        if cdata.as_ref().is_none_or(|best| data.len() < best.len()) {
            cdata = Some(data);
        }
    }

    let cdata = match cdata {
        Some(cdata) => cdata,
        None => rans_nx16::encode(Flags::empty(), buf)?,
    };

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

        for level in 1..=9 {
            let encoded = super::encode(&src, level)?;
            let decoded = super::super::decode(&encoded)?;
            assert_eq!(decoded, src, "level {level}");
        }

        // Sequential names cost next to nothing.
        let mut src = Vec::new();

        for i in 0..10_000 {
            src.extend_from_slice(format!("short{:09}", i / 2).as_bytes());
            src.push(0);
        }

        assert!(super::encode(&src, 5)?.len() < 500);

        Ok(())
    }

    /// Illumina names in an order where neighbours share little leave long number streams, which
    /// compress far better striped into a stream per byte, as htscodecs codes them: trying only
    /// unstriped flags, as this encoder once did, made these names over a fifth larger.
    /// A stream whose flags all stripe, but whose length does not divide by 4, falls back to
    /// plain order 0.
    #[test]
    fn test_unstripable_stream_falls_back_to_order_0() -> std::io::Result<()> {
        use crate::cram::codecs::rans_nx16;

        let buf = [1, 2, 3, 4, 5, 6];
        let mut dst = Vec::new();
        super::encode_token_byte_stream(&mut dst, super::Type::Digits, &buf, 1)?;

        // A type byte and a 1-byte length precede the rANS stream.
        let data = &dst[2..];
        assert_eq!(data[0] & rans_nx16::Flags::STRIPE.bits(), 0);
        assert_eq!(rans_nx16::decode(data, buf.len())?, buf);

        Ok(())
    }

    #[test]
    fn test_number_streams_are_striped() -> std::io::Result<()> {
        let mut src = Vec::new();
        let mut x = 0x9e37_79b9_u32;

        for _ in 0..10_000 {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;

            let name = format!(
                "A00744:46:HV3C3DSXX:{}:{}:{}:{}",
                1 + x % 4,
                1101 + (x >> 3) % 1500,
                1000 + (x >> 7) % 30000,
                1000 + (x >> 11) % 36000
            );

            src.extend_from_slice(name.as_bytes());
            src.push(0);
        }

        for level in [1, 5] {
            let encoded = super::encode(&src, level)?;
            assert_eq!(super::super::decode(&encoded)?, src, "level {level}");
            assert!(
                encoded.len() < 64_000,
                "level {level}: {} bytes",
                encoded.len()
            );
        }

        Ok(())
    }
}
