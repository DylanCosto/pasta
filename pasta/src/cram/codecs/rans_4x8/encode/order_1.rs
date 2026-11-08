use std::io::{self, Write};

use super::{EncodingSymbol, ReverseWriter, header, order_0, write_header, write_states};
use crate::cram::{
    codecs::{
        histogram,
        rans_4x8::{ALPHABET_SIZE, LOWER_BOUND, Order, STATE_COUNT},
    },
    io::writer::num::write_u8,
};

const CONTEXT_SIZE: usize = 2;
const NUL: u8 = 0x00;

/// Frequencies by context, then symbol (`ALPHABET_SIZE` rows), on the heap: only the rows of
/// contexts that occur are filled.
type RawFrequencies = Vec<[u32; ALPHABET_SIZE]>;
type Frequencies = Vec<[u16; ALPHABET_SIZE]>; // F
type CumulativeFrequencies = Frequencies; // C

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    // § 2.2.1 "rANS entropy encoding: Interleaving" (2023-03-15): "We do not permit Order-1
    // encoding of data streams smaller than 4 bytes."
    if src.len() < STATE_COUNT {
        return Err(io::Error::from(io::ErrorKind::InvalidInput));
    }

    let raw_frequencies = build_raw_frequencies(src);
    let contexts = build_alphabet(&raw_frequencies);
    let last = *src.last().expect("the source is not empty");
    let frequencies = normalize_frequencies(&raw_frequencies, last);
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies, &contexts);

    // Steps for the (context, symbol) pairs present, in one row per context that occurs.
    let mut rows = [0u8; ALPHABET_SIZE];
    let mut symbols: Vec<[EncodingSymbol; ALPHABET_SIZE]> = Vec::new();

    for (i, (fs, gs)) in frequencies.iter().zip(&cumulative_frequencies).enumerate() {
        if !contexts[i] {
            continue;
        }

        let mut row = [EncodingSymbol::default(); ALPHABET_SIZE];

        for ((symbol, &f), &g) in row.iter_mut().zip(fs).zip(gs) {
            if f > 0 {
                *symbol = EncodingSymbol::new(f, g);
            }
        }

        // At most 256 rows, so an index fits a byte.
        rows[i] = symbols.len() as u8;
        symbols.push(row);
    }

    // A row for each context, absent ones an empty one, so a step takes one lookup.
    let empty_row = [EncodingSymbol::default(); ALPHABET_SIZE];
    let rows: [&[EncodingSymbol; ALPHABET_SIZE]; ALPHABET_SIZE] = std::array::from_fn(|i| {
        if contexts[i] {
            &symbols[usize::from(rows[i])]
        } else {
            &empty_row
        }
    });

    let put = |state: &mut u32, i: u8, j: u8, buf: &mut ReverseWriter| {
        *state = rows[usize::from(i)][usize::from(j)].put(*state, buf);
    };

    let mut buf = ReverseWriter::for_symbols(src.len());
    let mut states = [LOWER_BOUND; STATE_COUNT];

    let [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4] = split_chunks(src);

    // § 2.2.1 "rANS entropy encoding: Interleaving" (2023-03-15): "Any remainder, when the input
    // buffer is not divisible by 4, is processed ... by the 4th rANS state."
    for syms in chunk_4.windows(CONTEXT_SIZE).rev() {
        put(&mut states[3], syms[0], syms[1], &mut buf);
    }

    // The decoder interleaves the states at each position, reading from a single stream, so
    // encode positions in reverse with the states in reverse at each one.
    let chunk_size = chunk_0.len();

    for k in (0..chunk_size.saturating_sub(1)).rev() {
        put(&mut states[3], chunk_3[k], chunk_3[k + 1], &mut buf);
        put(&mut states[2], chunk_2[k], chunk_2[k + 1], &mut buf);
        put(&mut states[1], chunk_1[k], chunk_1[k + 1], &mut buf);
        put(&mut states[0], chunk_0[k], chunk_0[k + 1], &mut buf);
    }

    let chunks = [chunk_0, chunk_1, chunk_2, chunk_3];

    for (state, chunk) in states.iter_mut().rev().zip(chunks.iter().rev()) {
        put(state, NUL, chunk[0], &mut buf);
    }

    let mut dst = vec![0; header::SIZE];

    write_frequencies(&mut dst, &frequencies)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.as_slice());

    let compressed_size = dst[header::SIZE..].len();
    // SAFETY: `dst.len() >= header::SIZE`.
    let header_dst = dst.first_chunk_mut().unwrap();
    write_header(header_dst, Order::One, compressed_size, src.len())?;

    Ok(dst)
}

fn write_frequencies<W>(writer: &mut W, frequencies: &[[u16; ALPHABET_SIZE]]) -> io::Result<()>
where
    W: Write,
{
    let alphabet = build_alphabet(frequencies);

    let mut iter = alphabet.iter().zip(frequencies).enumerate();
    let mut prev_sym = None;

    while let Some((sym, (&a, f))) = iter.next() {
        if !a {
            continue;
        }

        // SAFETY: `sym < ALPHABET_SIZE`.
        write_u8(writer, sym as u8)?;

        if sym > 0 && prev_sym == Some(sym - 1) {
            let i = sym + 1;

            let len = alphabet[i..]
                .iter()
                .position(|&a| !a)
                .unwrap_or(ALPHABET_SIZE - i);

            // SAFETY: `len < ALPHABET_SIZE`.
            write_u8(writer, len as u8)?;

            order_0::write_frequencies(writer, f)?;

            for (sym, (_, g)) in iter.by_ref().take(len) {
                order_0::write_frequencies(writer, g)?;
                prev_sym = Some(sym);
            }

            continue;
        }

        order_0::write_frequencies(writer, f)?;

        prev_sym = Some(sym);
    }

    write_u8(writer, NUL)?;

    Ok(())
}

fn build_alphabet<T>(frequencies: &[[T; ALPHABET_SIZE]]) -> [bool; ALPHABET_SIZE]
where
    T: Copy + Default + PartialEq,
{
    let mut alphabet = [false; ALPHABET_SIZE];

    for (a, fs) in alphabet.iter_mut().zip(frequencies) {
        *a = fs.iter().any(|&g| g != T::default());
    }

    alphabet
}

fn build_raw_frequencies(src: &[u8]) -> RawFrequencies {
    let mut frequencies = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    // § 2.1.2 "Frequency table: Order-1 encoding": "We use the ASCII value (`\0`) as the starting
    // context for each interleaved rANS state..."
    let chunk_size = src.len() / STATE_COUNT;

    for chunk in src.chunks_exact(chunk_size).take(STATE_COUNT) {
        // SAFETY: `chunk.len() > 0`.
        let (i, j) = (usize::from(NUL), usize::from(chunk[0]));
        frequencies[i][j] += 1;
    }

    let rows = <&mut [_; ALPHABET_SIZE]>::try_from(&mut frequencies[..]).expect("a row per symbol");
    histogram::add_order_1(src, rows);

    frequencies
}

/// Normalizes each context's frequencies to add up to 4095, as htscodecs' `rans_compress_O1`
/// does, so that the output is the same.
///
/// htscodecs counts the last symbol, `last`, once more in its context's total (in `hist1_4`), so
/// that context gets a row even if no symbol follows it.
// Kept out of `encode`, whose loop otherwise compiles slower.
#[inline(never)]
fn normalize_frequencies(raw_frequencies: &RawFrequencies, last: u8) -> Frequencies {
    let mut frequencies = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for (i, (f, g)) in raw_frequencies.iter().zip(&mut frequencies).enumerate() {
        let total = f.iter().sum::<u32>() + u32::from(i == usize::from(last));

        if total > 0 {
            *g = normalize_row(f, total);
        }
    }

    frequencies
}

/// Normalizes a context's frequencies `raw_frequencies`, which add up to `total` (or one less),
/// as htscodecs does: each is scaled by [`order_0::TOTAL`] over `total` in floating point,
/// rounded down but kept above 0, and the most frequent symbol (the first, if several) takes up
/// the difference, unless that would take at least half of it, when all are scaled by 0.98
/// first.
fn normalize_row(raw_frequencies: &[u32; ALPHABET_SIZE], total: u32) -> [u16; ALPHABET_SIZE] {
    use order_0::TOTAL;

    let mut scale = f64::from(TOTAL) / f64::from(total);
    let mut frequencies = *raw_frequencies;

    loop {
        let (mut sum, mut max, mut max_index) = (1, 0, 0);

        for (i, f) in frequencies.iter_mut().enumerate().filter(|(_, f)| **f > 0) {
            if *f > max {
                (max, max_index) = (*f, i);
            }

            *f = ((f64::from(*f) * scale) as u32).max(1);
            sum += *f;
        }

        let f = &mut frequencies[max_index];

        if sum < TOTAL {
            *f += TOTAL - sum;
        } else if sum - TOTAL >= *f / 2 {
            scale = 0.98;
            continue;
        } else {
            *f -= sum - TOTAL;
        }

        // Each is at most `TOTAL` - 1.
        return frequencies.map(|f| f as u16);
    }
}

fn build_cumulative_frequencies(
    frequencies: &Frequencies,
    contexts: &[bool; ALPHABET_SIZE],
) -> CumulativeFrequencies {
    let mut cumulative_frequencies = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for ((f, g), _) in frequencies
        .iter()
        .zip(&mut cumulative_frequencies)
        .zip(contexts)
        .filter(|(_, c)| **c)
    {
        *g = order_0::build_cumulative_frequencies(f);
    }

    cumulative_frequencies
}

fn split_chunks(buf: &[u8]) -> [&[u8]; 5] {
    let chunk_size = buf.len() / STATE_COUNT;

    let (left_chunk, right_chunk) = buf.split_at(2 * chunk_size);
    let (chunk_0, chunk_1) = left_chunk.split_at(chunk_size);
    let (chunk_2, chunk_3_4) = right_chunk.split_at(chunk_size);
    let (chunk_3, _) = chunk_3_4.split_at(chunk_size);

    // The last chunk includes the last byte of chunk 3 for context.
    let chunk_4 = &chunk_3_4[chunk_size - 1..];

    [chunk_0, chunk_1, chunk_2, chunk_3, chunk_4]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07):
        // "abracadabraabracadabraabracadabraabracadabrad".
        let mut frequencies = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        frequencies[usize::from(NUL)][usize::from(b'a')] = 4;
        frequencies[usize::from(b'a')][usize::from(b'a')] = 3;
        frequencies[usize::from(b'a')][usize::from(b'b')] = 8;
        frequencies[usize::from(b'a')][usize::from(b'c')] = 4;
        frequencies[usize::from(b'a')][usize::from(b'd')] = 5;
        frequencies[usize::from(b'b')][usize::from(b'r')] = 8;
        frequencies[usize::from(b'c')][usize::from(b'a')] = 4;
        frequencies[usize::from(b'd')][usize::from(b'a')] = 4;
        frequencies[usize::from(b'r')][usize::from(b'a')] = 8;

        let mut dst = Vec::new();
        write_frequencies(&mut dst, &frequencies)?;

        let expected = [
            0x00, // symbols[0] = '\0' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['\0']['a'] = 4
            0x00, // }
            b'a', // symbols[0] = 'a' {
            b'a', //   symbols[1] = 'a'
            0x03, //     frequencies['a']['a'] = 3
            b'b', //   symbols[1] = 'b'
            0x02, //     run length = 2
            0x08, //     frequencies['a']['b'] = 8
            0x04, //     frequencies['a']['c'] = 4
            0x05, //     frequencies['a']['d'] = 5
            0x00, // }
            b'b', // symbols[0] = 'b' {
            0x02, //   run length = 2
            b'r', //   symbols[1] = 'r'
            0x08, //     frequencies['b']['r'] = 8
            0x00, // }
            //    // symbols[0] = 'c' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['c']['a'] = 4
            0x00, // }
            //    // symbols[0] = 'd' {
            b'a', //   symbols[1] = 'a'
            0x04, //     frequencies['d']['a'] = 4
            0x00, // }
            b'r', // symbols[0] = 'r' {
            b'a', //   symbols[1] = 'a'
            0x08, //     frequencies['r']['a'] = 8
            0x00, // }
            0x00, // EOF
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        let src = b"abracadabraabracadabraabracadabraabracadabrad";
        let frequencies = build_raw_frequencies(src);

        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(NUL)] = true;
        expected[usize::from(b'a')] = true;
        expected[usize::from(b'b')] = true;
        expected[usize::from(b'c')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'r')] = true;

        assert_eq!(build_alphabet(&frequencies), expected);

        let src = b"noodles";
        let frequencies = build_raw_frequencies(src);

        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(NUL)] = true;
        expected[usize::from(b'n')] = true;
        expected[usize::from(b'o')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'l')] = true;
        expected[usize::from(b'e')] = true;

        assert_eq!(build_alphabet(&frequencies), expected);
    }

    #[test]
    fn test_build_raw_frequencies() {
        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        let src = b"abracadabraabracadabraabracadabraabracadabrad";

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4;
        expected[usize::from(b'a')][usize::from(b'a')] = 3;
        expected[usize::from(b'a')][usize::from(b'b')] = 8;
        expected[usize::from(b'a')][usize::from(b'c')] = 4;
        expected[usize::from(b'a')][usize::from(b'd')] = 5;
        expected[usize::from(b'b')][usize::from(b'r')] = 8;
        expected[usize::from(b'c')][usize::from(b'a')] = 4;
        expected[usize::from(b'd')][usize::from(b'a')] = 4;
        expected[usize::from(b'r')][usize::from(b'a')] = 8;

        assert_eq!(build_raw_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() {
        // The example of § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07), which
        // rounds differently from htscodecs: it has 1639 and 1023 for `a` before `b` and `d`.
        let src = b"abracadabraabracadabraabracadabraabracadabrad";
        let raw_frequencies = build_raw_frequencies(src);

        let mut expected = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 4095;
        expected[usize::from(b'a')][usize::from(b'a')] = 614;
        expected[usize::from(b'a')][usize::from(b'b')] = 1638;
        expected[usize::from(b'a')][usize::from(b'c')] = 819;
        expected[usize::from(b'a')][usize::from(b'd')] = 1024;
        expected[usize::from(b'b')][usize::from(b'r')] = 4095;
        expected[usize::from(b'c')][usize::from(b'a')] = 4095;
        expected[usize::from(b'd')][usize::from(b'a')] = 4095;
        expected[usize::from(b'r')][usize::from(b'a')] = 4095;

        assert_eq!(normalize_frequencies(&raw_frequencies, b'd'), expected);

        // The last symbol's context has a row even if no symbol follows it there.
        let raw_frequencies = build_raw_frequencies(b"aaab");
        let frequencies = normalize_frequencies(&raw_frequencies, b'b');
        assert_eq!(frequencies[usize::from(b'b')][usize::from(NUL)], 4095);
    }

    #[test]
    fn test_split_chunks() {
        assert_eq!(
            split_chunks(&[1, 2, 3, 4]),
            [&[1][..], &[2][..], &[3][..], &[4][..], &[4][..]]
        );

        assert_eq!(
            split_chunks(&[1, 2, 3, 4, 5]),
            [&[1][..], &[2][..], &[3][..], &[4][..], &[4, 5][..]]
        );

        // § 2.1.2 "Frequency table: Order-1 encoding" (563e8ab 2025-04-07)
        assert_eq!(
            split_chunks(b"abracadabraabracadabraabracadabraabracadabrad"),
            [
                &b"abracadabra"[..],
                &b"abracadabra"[..],
                &b"abracadabra"[..],
                &b"abracadabra"[..],
                &b"ad"[..],
            ]
        );
    }
}
