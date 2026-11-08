mod trie;

use std::{io, ops::Range};

use self::trie::Trie;
use super::Type;
use crate::cram::io::writer::num::{write_u8, write_u32_le, write_uint7};

const NUL: u8 = 0x00;

/// The token positions a name can use, counting position 0 (whether it is a duplicate, and of
/// which name) and its end. htscodecs' decoder accepts no more (`MAX_TOKENS`).
const MAX_TOKENS: usize = 128;

/// The most digits read as one number, as htscodecs reads them, so that any fits in 32 bits.
const MAX_DIGITS: usize = 9;

/// The types of the streams of a token position, in the order they are written: the types of its
/// tokens, then their values.
const STREAM_TYPES: [Type; 10] = [
    Type::Type,
    Type::String,
    Type::Char,
    Type::Digits0,
    Type::DZLen,
    Type::Dup,
    Type::Diff,
    Type::Digits,
    Type::Delta,
    Type::Delta0,
];

/// The flag on the first stream of a token position.
const NEW_POSITION: u8 = 0x80;

/// The flag on a stream that is a copy of an earlier one.
const DUPLICATE: u8 = 0x40;

/// Compresses NUL-terminated read names with the name tokenizer, coding its token streams with
/// rANS Nx16 at a CRAM compression `level` (1–9), as HTSlib's `tok3_encode_names` does.
pub fn encode(src: &[u8], level: u8) -> io::Result<Vec<u8>> {
    let names: Vec<_> = match src.strip_suffix(&[NUL]) {
        Some(names) => names.split(|&b| b == NUL).collect(),
        None if src.is_empty() => Vec::new(),
        None => src.split(|&b| b == NUL).collect(),
    };

    let mut encoder = Encoder::new();

    for &name in &names {
        encoder.push(name)?;
    }

    // The decoded length: each name and its NUL.
    let len = names.iter().map(|name| name.len() + 1).sum();

    let mut dst = Vec::new();
    write_header(&mut dst, len, names.len())?;
    encoder.write_streams(&mut dst, level)?;

    Ok(dst)
}

fn write_header(dst: &mut Vec<u8>, len: usize, names_count: usize) -> io::Result<()> {
    let ulen = u32::try_from(len).map_err(invalid_input)?;
    write_u32_le(dst, ulen)?;

    let n_names = u32::try_from(names_count).map_err(invalid_input)?;
    write_u32_le(dst, n_names)?;

    // TODO: use_arith
    write_u8(dst, 0)?;

    Ok(())
}

/// The token streams of names so far, and what htscodecs keeps to encode the next name.
struct Encoder<'a> {
    trie: Trie<'a>,
    names: Vec<Context<'a>>,
    tokens: Vec<Token>,
    positions: Vec<Position>,
}

/// A name as later names are encoded against it (htscodecs' `last_context`).
#[derive(Clone)]
struct Context<'a> {
    name: &'a [u8],
    /// Its tokens in [`Encoder::tokens`], from position 1.
    tokens: Range<usize>,
}

/// A token as later names compare against it (htscodecs' `last_context_tok`).
#[derive(Clone, Copy, Debug)]
enum Token {
    /// Letters and punctuation, at `start..start + len` of its name.
    Alpha {
        start: usize,
        len: usize,
    },
    Char(u8),
    /// A number padded with zeros to `width` digits.
    Digits0 {
        value: u32,
        width: usize,
    },
    Digits(u32),
}

/// The streams of a token position.
#[derive(Default)]
struct Position {
    /// The stream of each type of [`STREAM_TYPES`].
    streams: [Vec<u8>; STREAM_TYPES.len()],
    /// The numbers stored as deltas from the number in this position of the name they are
    /// encoded against, and those that could not be (htscodecs' `token_dcount` and
    /// `token_icount`). A number is stored as a delta only while the second count is less than
    /// the first plus 5.
    delta_count: u32,
    digits_count: u32,
}

impl<'a> Encoder<'a> {
    fn new() -> Self {
        Self {
            trie: Trie::new(),
            names: Vec::new(),
            tokens: Vec::new(),
            positions: vec![Position::default()],
        }
    }

    /// Encodes a name as htscodecs' `encode_name` does: against the last name equal to it or that
    /// it begins, if any; otherwise the last name sharing the prefix its format suggests, if any;
    /// otherwise the name before it.
    fn push(&mut self, name: &'a [u8]) -> io::Result<()> {
        let n = self.names.len();
        let id = u32::try_from(n).map_err(invalid_input)?;

        let format = Format::guess(name);
        let visit = self.trie.visit(name, id, format.prefix_len);
        let exact = visit.full.filter(|&m| m != id);

        let m = exact
            .or(visit.prefix)
            .map_or(n.saturating_sub(1), |m| m as usize);

        let distance = (n - m) as u32;

        if exact.is_some() && self.names[m].name.len() == name.len() {
            self.positions[0].push_u32(Type::Dup, distance);

            let tokens = self.names[m].tokens.clone();
            self.names.push(Context { name, tokens });

            return Ok(());
        }

        self.positions[0].push_u32(Type::Diff, distance);

        // The first name sharing the prefix is encoded against itself, which is nothing.
        let prev = match self.names.get(m) {
            Some(prev) => prev.clone(),
            None => Context {
                name: &[],
                tokens: 0..0,
            },
        };

        let tokens_start = self.tokens.len();
        let mut i = 0;
        let mut pos = 1;

        let prev_token = |tokens: &[Token], pos: usize| {
            let k = prev.tokens.start + pos - 1;
            (k < prev.tokens.end).then(|| tokens[k])
        };

        // The prefix of a known format is one token.
        if format.is_fixed {
            let prefix = &name[..format.prefix_len];
            let position = position(&mut self.positions, pos)?;

            match prev_token(&self.tokens, pos) {
                Some(Token::Alpha { len, .. })
                    if len == prefix.len() && prev.name.starts_with(prefix) =>
                {
                    position.push_type(Type::Match);
                }
                _ => position.push_string(prefix),
            }

            self.tokens.push(Token::Alpha {
                start: 0,
                len: prefix.len(),
            });

            i = prefix.len();
            pos += 1;
        }

        while i < name.len() {
            let prev_token = prev_token(&self.tokens, pos);
            let position = position(&mut self.positions, pos)?;
            let (token, end) = push_token(position, name, i, prev.name, prev_token);

            self.tokens.push(token);
            i = end;
            pos += 1;
        }

        position(&mut self.positions, pos)?.push_type(Type::End);

        self.names.push(Context {
            name,
            tokens: tokens_start..self.tokens.len(),
        });

        Ok(())
    }

    /// Writes each token position's streams, as htscodecs' `tok3_encode_names` does.
    fn write_streams(&self, dst: &mut Vec<u8>, level: u8) -> io::Result<()> {
        struct Block {
            pos: u8,
            ty: Type,
            data: Vec<u8>,
        }

        let mut blocks = Vec::new();

        for (pos, position) in self.positions.iter().enumerate() {
            let [types, values @ ..] = &position.streams;

            // A decoder takes a type stream that is missing for the first type of the position's
            // first stream, then matches.
            let is_implicit = types
                .split_first()
                .is_some_and(|(_, rest)| rest.iter().all(|&ty| ty == u8::from(Type::Match)))
                && values.iter().any(|stream| !stream.is_empty());

            for (&ty, stream) in STREAM_TYPES.iter().zip(&position.streams) {
                if stream.is_empty() || (ty == Type::Type && is_implicit) {
                    continue;
                }

                blocks.push(Block {
                    pos: pos as u8,
                    ty,
                    data: compress(ty, stream, level)?,
                });
            }
        }

        let mut prev_pos = None;

        for (i, block) in blocks.iter().enumerate() {
            let mut ty = u8::from(block.ty);

            if prev_pos != Some(block.pos) {
                ty |= NEW_POSITION;
                prev_pos = Some(block.pos);
            }

            // A stream coded the same as an earlier one refers to it.
            let original = blocks[..i]
                .iter()
                .find(|other| block.data.len() > 4 && other.data == block.data);

            match original {
                Some(original) => {
                    dst.extend([ty | DUPLICATE, original.pos, u8::from(original.ty)]);
                }
                None => {
                    dst.push(ty);
                    dst.extend_from_slice(&block.data);
                }
            }
        }

        Ok(())
    }
}

/// The streams of token position `pos`, which follows those so far.
fn position(positions: &mut Vec<Position>, pos: usize) -> io::Result<&mut Position> {
    if pos == positions.len() {
        if pos >= MAX_TOKENS {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "too many tokens in a read name",
            ));
        }

        positions.push(Position::default());
    }

    Ok(&mut positions[pos])
}

/// Encodes the token of `name` at `start` against the token in the same position of the previous
/// name, if any, returning the token and its end.
fn push_token(
    position: &mut Position,
    name: &[u8],
    start: usize,
    prev_name: &[u8],
    prev: Option<Token>,
) -> (Token, usize) {
    let rest = &name[start..];
    let b = rest[0];

    if b.is_ascii_alphabetic() {
        let len = 1 + rest[1..]
            .iter()
            .take_while(|c| c.is_ascii_alphabetic() || c.is_ascii_punctuation())
            .count();

        // A single letter is a character.
        if len > 1 {
            let s = &rest[..len];

            match prev {
                Some(Token::Alpha { start, len })
                    if prev_name.get(start..start + len) == Some(s) =>
                {
                    position.push_type(Type::Match);
                }
                _ => position.push_string(s),
            }

            return (Token::Alpha { start, len }, start + len);
        }
    } else if b.is_ascii_digit() {
        let (len, value) = parse_digits(rest);

        // A number with a leading zero keeps its width, as does one as wide as a padded number
        // before it.
        let token = match prev {
            Some(Token::Digits0 { value: prev, width }) if width == len => {
                match value.checked_sub(prev) {
                    Some(0) => position.push_type(Type::Match),
                    Some(delta @ 1..=255) => position.push_u8(Type::Delta0, delta as u8),
                    _ => position.push_digits0(value, len),
                }

                Token::Digits0 { value, width: len }
            }
            _ if b == b'0' => {
                position.push_digits0(value, len);
                Token::Digits0 { value, width: len }
            }
            Some(Token::Digits(prev)) => {
                match value.checked_sub(prev) {
                    Some(0) => position.push_type(Type::Match),
                    Some(delta @ 1..=255) if position.digits_count < position.delta_count + 5 => {
                        position.push_u8(Type::Delta, delta as u8);
                        position.delta_count += 1;
                    }
                    _ => {
                        position.push_u32(Type::Digits, value);
                        position.digits_count += 1;
                    }
                }

                Token::Digits(value)
            }
            _ => {
                position.push_u32(Type::Digits, value);
                Token::Digits(value)
            }
        };

        return (token, start + len);
    }

    match prev {
        Some(Token::Char(c)) if c == b => position.push_type(Type::Match),
        _ => position.push_char(b),
    }

    (Token::Char(b), start + 1)
}

/// Reads up to [`MAX_DIGITS`] digits at the start of `src`, returning their count and value.
fn parse_digits(src: &[u8]) -> (usize, u32) {
    let len = src
        .iter()
        .take(MAX_DIGITS)
        .take_while(|b| b.is_ascii_digit())
        .count();

    let value = src[..len]
        .iter()
        .fold(0, |n, &b| n * 10 + u32::from(b - b'0'));

    (len, value)
}

impl Position {
    fn stream(&mut self, ty: Type) -> &mut Vec<u8> {
        &mut self.streams[usize::from(ty)]
    }

    fn push_type(&mut self, ty: Type) {
        self.stream(Type::Type).push(u8::from(ty));
    }

    fn push_string(&mut self, s: &[u8]) {
        self.push_type(Type::String);

        let stream = self.stream(Type::String);
        stream.extend_from_slice(s);
        stream.push(NUL);
    }

    fn push_char(&mut self, c: u8) {
        self.push_u8(Type::Char, c);
    }

    fn push_digits0(&mut self, value: u32, width: usize) {
        self.push_u32(Type::Digits0, value);
        // At most `MAX_DIGITS`.
        self.stream(Type::DZLen).push(width as u8);
    }

    fn push_u8(&mut self, ty: Type, n: u8) {
        self.push_type(ty);
        self.stream(ty).push(n);
    }

    fn push_u32(&mut self, ty: Type, n: u32) {
        self.push_type(ty);
        self.stream(ty).extend(n.to_le_bytes());
    }
}

/// What htscodecs guesses of a name's format from its shape (in `search_trie`): how long a prefix
/// to look for an earlier name sharing, and whether it is a token of its own.
#[derive(Debug, Eq, PartialEq)]
struct Format {
    prefix_len: usize,
    is_fixed: bool,
}

impl Format {
    fn guess(name: &[u8]) -> Self {
        let d = name.strip_prefix(b"@").unwrap_or(name);
        let f = usize::from(name.first() == Some(&b'>'));
        let is_hex = |b: u8| b.is_ascii_digit() || (b'a'..=b'f').contains(&b);

        if d.len() > 70 && d[f] == b'm' && d[7] == b'_' && d[f + 14] == b'_' && d[f + 61] == b'/' {
            // PacBio
            Self {
                prefix_len: 60,
                is_fixed: false,
            }
        } else if d.len() == 17 && d[f + 5] == b':' && d[f + 11] == b':' {
            // Ion Torrent
            Self {
                prefix_len: 6,
                is_fixed: true,
            }
        } else if d.len() > 37
            && [8, 13, 18, 23].iter().all(|&i| d[f + i] == b'-')
            && is_hex(d[f])
            && is_hex(d[f + 35])
        {
            // Oxford Nanopore: a UUID.
            Self {
                prefix_len: 37,
                is_fixed: true,
            }
        } else {
            // Illumina: the prefix before the lane, tile and coordinates.
            let end = name.iter().position(|&b| b <= b' ').unwrap_or(name.len());

            let fourth_last_colon = name[..end]
                .iter()
                .rposition(|&b| b == b':')
                .and_then(|i| name[..i].iter().rposition(|&b| b == b':'))
                .and_then(|i| name[..i].iter().rposition(|&b| b == b':'))
                .and_then(|i| name[..i].iter().rposition(|&b| b == b':'));

            match fourth_last_colon {
                Some(i) => Self {
                    prefix_len: i + 1,
                    is_fixed: true,
                },
                None => Self {
                    prefix_len: usize::MAX,
                    is_fixed: false,
                },
            }
        }
    }
}

/// The rANS Nx16 flags htscodecs tries on a token stream of type `ty` (the table `R` of
/// `compress` in `tokenise_name3.c`). HTSlib caps the level at 3 when it codes names with rANS
/// (`cram_compress_by_method`), so only two of the table's rows apply: one for levels 1 and 2,
/// and one for 3 and above. The number streams hold 4-byte values, which striping (`0x08`)
/// splits into a stream per byte.
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

/// Codes a token stream with each of htscodecs' flags for its type and keeps the smallest, after
/// its length.
fn compress(ty: Type, src: &[u8], level: u8) -> io::Result<Vec<u8>> {
    use crate::cram::codecs::rans_nx16::{self, Flags};

    // Striping needs a length that divides by 4, which the number streams always have; plain
    // order 0 stands in if no flags are left.
    let flags = rans_flags(ty, level)
        .iter()
        .copied()
        .filter(|&flags| src.len().is_multiple_of(4) || flags & Flags::STRIPE.bits() == 0);

    let mut cdata: Option<Vec<u8>> = None;

    for flags in flags {
        let data = rans_nx16::encode(Flags::from_bits_truncate(flags), src)?;

        if cdata.as_ref().is_none_or(|best| data.len() < best.len()) {
            cdata = Some(data);
        }
    }

    let cdata = match cdata {
        Some(cdata) => cdata,
        None => rans_nx16::encode(Flags::empty(), src)?,
    };

    let clen = u32::try_from(cdata.len()).map_err(invalid_input)?;

    let mut dst = Vec::with_capacity(cdata.len() + 5);
    write_uint7(&mut dst, clen)?;
    dst.extend(cdata);

    Ok(dst)
}

fn invalid_input<E>(e: E) -> io::Error
where
    E: Into<Box<dyn std::error::Error + Send + Sync>>,
{
    io::Error::new(io::ErrorKind::InvalidInput, e)
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
        let dst = super::compress(super::Type::Digits, &buf, 1)?;

        // A 1-byte length precedes the rANS stream.
        let data = &dst[1..];
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

    /// The output is htscodecs' byte for byte, which picks the same names to encode against,
    /// the same tokens, and the same ways of coding their streams.
    #[test]
    fn test_encode_matches_htscodecs() -> std::io::Result<()> {
        let src = include_bytes!("../fixtures/names-formats.bin");

        for (level, expected) in [
            (
                1,
                &include_bytes!("../fixtures/tok3-rans-formats-1.bin")[..],
            ),
            (
                3,
                &include_bytes!("../fixtures/tok3-rans-formats-3.bin")[..],
            ),
        ] {
            assert_eq!(super::encode(src, level)?, expected, "level {level}");
        }

        Ok(())
    }

    #[test]
    fn test_format_guess() {
        use super::Format;

        let fixed = |prefix_len| Format {
            prefix_len,
            is_fixed: true,
        };

        assert_eq!(
            Format::guess(b"A00744:46:HV3C3DSXX:2:1131:20076:26021"),
            fixed(20)
        );
        assert_eq!(Format::guess(b"A:1:2:3:4"), fixed(2));
        assert_eq!(Format::guess(b"ZEN10:00027:00321"), fixed(6));
        assert_eq!(
            Format::guess(b"f33d30d5-6eb8-4115-8f46-154c2620a5da_Basecall_1D_template"),
            fixed(37)
        );
        assert_eq!(
            Format::guess(b"m64011_190830_220126/1/ccs"),
            Format {
                prefix_len: usize::MAX,
                is_fixed: false
            }
        );
        assert_eq!(
            Format::guess(
                b"m140415_143853_42175_c100635972550000001823121909121417_s1_p0/553/3100_11230"
            ),
            Format {
                prefix_len: 60,
                is_fixed: false
            }
        );
        assert_eq!(
            Format::guess(b""),
            Format {
                prefix_len: usize::MAX,
                is_fixed: false
            }
        );
    }

    #[test]
    fn test_encode_with_too_many_tokens() -> std::io::Result<()> {
        // Each "a1" is a character and a number: two tokens.
        let name = |n| "a1".repeat(n) + "\0";

        assert!(super::encode(name(63).as_bytes(), 3).is_ok());
        assert!(super::encode(name(64).as_bytes(), 3).is_err());

        let src = name(63) + &name(63);
        let encoded = super::encode(src.as_bytes(), 3)?;
        assert_eq!(super::super::decode(&encoded)?, src.as_bytes());

        Ok(())
    }

    #[test]
    fn test_encode_empty() -> std::io::Result<()> {
        assert_eq!(super::encode(b"", 3)?, [0, 0, 0, 0, 0, 0, 0, 0, 0]);
        Ok(())
    }
}
