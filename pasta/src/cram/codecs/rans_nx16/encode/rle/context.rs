use std::{error, fmt, io, num::NonZero};

use super::super::order_0;
use crate::cram::{codecs::rans_nx16::ALPHABET_SIZE, io::writer::num::write_uint7};

pub struct Context {
    pub alphabet: [bool; ALPHABET_SIZE],
    pub dst: Vec<u8>,
}

#[derive(Debug)]
pub enum BuildContextError {
    EmptyAlphabet,
}

impl error::Error for BuildContextError {}

impl fmt::Display for BuildContextError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyAlphabet => write!(f, "empty alphabet"),
        }
    }
}

pub fn build_context(src: &[u8]) -> Result<Context, BuildContextError> {
    let alphabet = count_savings(src).map(|n| n > 0);

    let n = alphabet.iter().filter(|&&a| a).count();
    let symbol_count = NonZero::new(n).ok_or(BuildContextError::EmptyAlphabet)?;

    let mut dst = Vec::new();
    write_alphabet(&mut dst, &alphabet, symbol_count);

    Ok(Context { alphabet, dst })
}

/// Counts, for each symbol, how many more times it repeats the one before than it does not
/// (htscodecs' `rle_find_syms`, which counts the first symbol as a run break, too). A block is
/// under 2^32 bytes, so this fits.
///
/// Counting into one table, each repeat in a run would wait for the last to be stored. As
/// htscodecs does, this counts alternate positions into four tables and adds them up.
fn count_savings(src: &[u8]) -> [i64; ALPHABET_SIZE] {
    let step = |prev: u8, sym: u8| if sym == prev { 1 } else { -1 };

    let Some((&first, rest)) = src.split_first() else {
        return [0; ALPHABET_SIZE];
    };

    let mut tables = [[0; ALPHABET_SIZE]; 4];
    tables[0][usize::from(first)] = -1;

    let mut prev = first;
    let (quads, remainder) = rest.as_chunks::<4>();

    for &[a, b, c, d] in quads {
        tables[0][usize::from(a)] += step(prev, a);
        tables[1][usize::from(b)] += step(a, b);
        tables[2][usize::from(c)] += step(b, c);
        tables[3][usize::from(d)] += step(c, d);
        prev = d;
    }

    for &sym in remainder {
        tables[0][usize::from(sym)] += step(prev, sym);
        prev = sym;
    }

    let [mut savings, t1, t2, t3] = tables;

    for (((n, a), b), c) in savings.iter_mut().zip(t1).zip(t2).zip(t3) {
        *n += a + b + c;
    }

    savings
}

fn write_alphabet(
    dst: &mut Vec<u8>,
    alphabet: &[bool; ALPHABET_SIZE],
    symbol_count: NonZero<usize>,
) {
    write_symbol_count(dst, symbol_count);

    for (sym, _) in alphabet.iter().enumerate().filter(|(_, a)| **a) {
        // SAFETY: `sym` < `ALPHABET_SIZE`.
        dst.push(sym as u8);
    }
}

fn write_symbol_count(dst: &mut Vec<u8>, symbol_count: NonZero<usize>) {
    let mut n = symbol_count.get();

    if n >= ALPHABET_SIZE {
        n = 0;
    }

    // SAFETY: `n` < `ALPHABET_SIZE`.
    dst.push(n as u8);
}

/// Writes the run-length metadata (the alphabet and run lengths), compressed with order-0 rANS
/// when that makes it smaller, as htscodecs does.
pub fn write_context(
    dst: &mut Vec<u8>,
    ctx: &Context,
    compressed_size: usize,
    state_count: usize,
) -> io::Result<()> {
    let len =
        u32::try_from(ctx.dst.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let compressed_size = u32::try_from(compressed_size)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut buf = Vec::new();
    let order_0_ctx = order_0::build_context(&ctx.dst)?;
    order_0::write_context(&mut buf, &order_0_ctx)?;
    order_0::encode(&ctx.dst, &order_0_ctx, state_count, &mut buf)?;

    if buf.len() < ctx.dst.len() {
        write_uint7(dst, len << 1)?;
        write_uint7(dst, compressed_size)?;

        let n =
            u32::try_from(buf.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(dst, n)?;

        dst.extend(&buf);
    } else {
        write_uint7(dst, (len << 1) | 1)?;
        write_uint7(dst, compressed_size)?;
        dst.extend(&ctx.dst);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_context() -> Result<(), BuildContextError> {
        let src = b"nnndlllllss";
        let ctx = build_context(src)?;

        let mut expected_alphabet = [false; ALPHABET_SIZE];
        expected_alphabet[usize::from(b'l')] = true;
        expected_alphabet[usize::from(b'n')] = true;
        assert_eq!(ctx.alphabet, expected_alphabet);

        assert_eq!(
            ctx.dst,
            [
                0x02, // symbol count
                b'l', b'n', // alphabet
            ]
        );

        Ok(())
    }

    #[test]
    fn test_count_savings() {
        for len in [0, 1, 2, 5, 8, 1001] {
            let src: Vec<u8> = (0..len).map(|i| b"aaab"[(i * 7 / 5) % 4]).collect();

            let mut expected = [0i64; ALPHABET_SIZE];
            let mut prev = None;

            for &sym in &src {
                expected[usize::from(sym)] += if prev == Some(sym) { 1 } else { -1 };
                prev = Some(sym);
            }

            assert_eq!(count_savings(&src), expected, "len {len}");
        }
    }

    #[test]
    fn test_write_alphabet() {
        let mut alphabet = [false; ALPHABET_SIZE];
        alphabet[usize::from(b'n')] = true;

        let symbol_count = const { NonZero::new(1).unwrap() };

        let mut dst = Vec::new();
        write_alphabet(&mut dst, &alphabet, symbol_count);

        assert_eq!(
            dst,
            [
                0x01, // symbol count = 1
                b'n', // alphabet
            ]
        );
    }

    #[test]
    fn test_write_symbol_count() {
        let mut dst = Vec::new();

        dst.clear();
        write_symbol_count(&mut dst, const { NonZero::new(1).unwrap() });
        assert_eq!(dst, [0x01]);

        dst.clear();
        write_symbol_count(&mut dst, const { NonZero::new(256).unwrap() });
        assert_eq!(dst, [0x00]);
    }

    #[test]
    fn test_write_context() -> Result<(), Box<dyn std::error::Error>> {
        let src = b"nnndlllllss";
        let ctx = build_context(src)?;

        let mut dst = Vec::new();
        write_context(&mut dst, &ctx, 5, 4)?;

        assert_eq!(
            dst,
            [
                0x07, // (context size = 3, is_compressed = false)
                0x05, // compressed size
                0x02, // symbol count
                b'l', b'n', // alphabet
            ]
        );

        Ok(())
    }
}
