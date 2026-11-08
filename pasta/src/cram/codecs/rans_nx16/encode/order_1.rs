use std::io;

use super::{EncodingSymbol, LOWER_BOUND, ReverseWriter, order_0, write_alphabet, write_states};
use crate::cram::{
    codecs::{histogram, rans_nx16::ALPHABET_SIZE},
    io::writer::num::{write_u8, write_uint7},
};

const CONTEXT_SIZE: usize = 2;
const NUL: u8 = 0x00;

/// The frequency table shifts htscodecs chooses between (`TF_SHIFT_O1_FAST` and `TF_SHIFT_O1`).
const FAST_BITS: u32 = 10;
const BITS: u32 = 12;

/// Frequency tables longer than this (with their header byte) are compressed with order-0 rANS
/// when that saves space, as htscodecs does.
const MAX_UNCOMPRESSED_TABLE_LEN: usize = 1000;

/// The state count of a compressed frequency table (htscodecs' `rans_compress_O0_4x16`).
const TABLE_STATE_COUNT: usize = 4;

/// Frequencies by context, then symbol (`ALPHABET_SIZE` rows), on the heap: rows of absent
/// contexts stay zero.
type Frequencies = Vec<[u32; ALPHABET_SIZE]>;

pub(super) struct Context {
    alphabet: [bool; ALPHABET_SIZE],
    /// The frequencies encoded with: each context's sum to `1 << bits`.
    frequencies: Frequencies,
    bits: u32,
    /// The frequency tables as written after the header byte, compressed or not.
    table: Vec<u8>,
    is_table_compressed: bool,
}

/// Builds the frequency tables as htscodecs' `encode_freq1` does: each context's is normalized
/// to the power of two its count rounds to (written so, it costs a byte or two per symbol less
/// than at the full total, which decoders scale back up), and to 10-bit totals instead of
/// 12-bit ones when that is estimated to cost no more than 1% in size.
pub(super) fn build_context(src: &[u8], state_count: usize) -> io::Result<Context> {
    let alphabet = build_alphabet(src);
    let mut frequencies = build_frequencies(src, state_count);

    // The last symbol counts once more, as in htscodecs' `hist1_4`.
    let mut totals = [0; ALPHABET_SIZE];

    for (total, fs) in totals.iter_mut().zip(&frequencies) {
        *total = fs.iter().sum();
    }

    if let Some(&last) = src.last() {
        totals[usize::from(last)] += 1;
    }

    let (bits, maxes) = choose_bits(&frequencies, &totals);

    for i in contexts(&alphabet) {
        let max = if bits == FAST_BITS {
            maxes[i].min(1 << FAST_BITS)
        } else {
            maxes[i]
        };

        order_0::normalize_frequencies(&mut frequencies[i], totals[i], max)?;
    }

    let mut table = Vec::new();
    write_alphabet(&mut table, &alphabet)?;
    write_frequencies(&mut table, &alphabet, &frequencies)?;

    for i in contexts(&alphabet) {
        let max = frequencies[i].iter().sum();
        order_0::scale_frequencies(&mut frequencies[i], max, 1 << bits);
    }

    let mut is_table_compressed = false;

    if 1 + table.len() > MAX_UNCOMPRESSED_TABLE_LEN {
        let ctx = order_0::build_context(&table)?;
        let mut compressed_table = Vec::new();
        order_0::write_context(&mut compressed_table, &ctx)?;
        order_0::encode(&table, &ctx, TABLE_STATE_COUNT, &mut compressed_table)?;

        if compressed_table.len() + 6 < 1 + table.len() {
            let len = |n: usize| {
                u32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            };

            let mut buf = Vec::with_capacity(compressed_table.len() + 10);
            write_uint7(&mut buf, len(table.len())?)?;
            write_uint7(&mut buf, len(compressed_table.len())?)?;
            buf.extend(compressed_table);

            table = buf;
            is_table_compressed = true;
        }
    }

    Ok(Context {
        alphabet,
        frequencies,
        bits,
        table,
        is_table_compressed,
    })
}

pub(super) fn write_context(dst: &mut Vec<u8>, ctx: &Context) -> io::Result<()> {
    // The shift is 10 or 12, so this fits a byte.
    write_u8(
        dst,
        (ctx.bits << 4) as u8 | u8::from(ctx.is_table_compressed),
    )?;
    dst.extend(&ctx.table);
    Ok(())
}

/// Chooses 10-bit or 12-bit frequency totals from an estimate of the encoded size with each, and
/// the power of two to write each context's frequencies at (htscodecs' `rans_compute_shift`).
fn choose_bits(
    frequencies: &Frequencies,
    totals: &[u32; ALPHABET_SIZE],
) -> (u32, [u32; ALPHABET_SIZE]) {
    let (fast_total, total) = (1u32 << FAST_BITS, 1u32 << BITS);

    let (mut e10, mut e12) = (0.0, 0.0);
    let mut max_total = 0;
    let mut maxes = [0; ALPHABET_SIZE];

    for ((fs, &t), max) in frequencies.iter().zip(totals).zip(&mut maxes) {
        if t == 0 {
            continue;
        }

        let mut max_val = t.next_power_of_two();

        // The symbols whose frequencies would be raised to 1.
        let raised = |limit| fs.iter().filter(|&&f| f > 0 && max_val / f > limit).count();
        let l10 = ((fast_total as usize + raised(fast_total)) as f64).ln();
        let l12 = ((total as usize + raised(total)) as f64).ln();

        let (t_fast, t_slow) = (
            f64::from(fast_total) / f64::from(t),
            f64::from(total) / f64::from(t),
        );
        let mut symbol_count = 0;

        for &f in fs.iter().filter(|&&f| f > 0) {
            symbol_count += 1;

            let f = f64::from(f);
            e10 -= f * (fast_log((f * t_fast).max(1.0)) - l10);
            e12 -= f * (fast_log((f * t_slow).max(1.0)) - l12);

            // An estimate of the cost of the frequency tables.
            e10 += 1.3;
            e12 += 4.7;
        }

        if symbol_count < 64 && max_val > 128 {
            max_val /= 2;
        }

        if max_val > 1024 {
            max_val /= 2;
        }

        *max = max_val.min(total);
        max_total = max_total.max(*max);
    }

    let bits = if e10 / e12 < 1.01 || max_total <= fast_total {
        FAST_BITS
    } else {
        BITS
    };

    (bits, maxes)
}

/// An approximation of the natural logarithm from the bits of a double (htscodecs' `fast_log`),
/// which this matches, so that the same shift is chosen.
fn fast_log(a: f64) -> f64 {
    (a.to_bits() as i64 - 4_606_921_278_410_026_770) as f64 * 1.539_095_918_623_324e-16
}

pub fn encode(src: &[u8], ctx: &Context, state_count: usize, dst: &mut Vec<u8>) -> io::Result<()> {
    let frequencies = &ctx.frequencies;
    let cumulative_frequencies = build_cumulative_frequencies(frequencies, &ctx.alphabet);

    // Steps for the (context, symbol) pairs present, in one row per context that occurs.
    let mut rows = [0u8; ALPHABET_SIZE];
    let mut symbols: Vec<[EncodingSymbol; ALPHABET_SIZE]> = Vec::new();

    for i in contexts(&ctx.alphabet) {
        let (fs, gs) = (&frequencies[i], &cumulative_frequencies[i]);
        let mut row = [EncodingSymbol::default(); ALPHABET_SIZE];

        for ((symbol, &f), &g) in row.iter_mut().zip(fs).zip(gs) {
            if f > 0 {
                *symbol = EncodingSymbol::new(f, g, ctx.bits);
            }
        }

        // At most 256 rows, so an index fits a byte.
        rows[i] = symbols.len() as u8;
        symbols.push(row);
    }

    // A row for each context, absent ones an empty one, so a step takes one lookup.
    let empty_row = [EncodingSymbol::default(); ALPHABET_SIZE];
    let rows: [&[EncodingSymbol; ALPHABET_SIZE]; ALPHABET_SIZE] = std::array::from_fn(|i| {
        if ctx.alphabet[i] {
            &symbols[usize::from(rows[i])]
        } else {
            &empty_row
        }
    });

    let mut buf = ReverseWriter::for_symbols(src.len());

    match state_count {
        4 => write_states(dst, &encode_states::<4>(src, &rows, &mut buf))?,
        32 => write_states(dst, &encode_states::<32>(src, &rows, &mut buf))?,
        _ => unreachable!("state count is 4 or 32"),
    }

    dst.extend(buf.as_slice());

    Ok(())
}

/// Encodes `src` on `N` states, each a chunk of it, the last also taking the rest, returning the
/// final states.
fn encode_states<const N: usize>(
    src: &[u8],
    rows: &[&[EncodingSymbol; ALPHABET_SIZE]; ALPHABET_SIZE],
    buf: &mut ReverseWriter,
) -> [u32; N] {
    let put = |state: &mut u32, i: u8, j: u8, buf: &mut ReverseWriter| {
        *state = rows[usize::from(i)][usize::from(j)].put(*state, buf);
    };

    let mut states = [LOWER_BOUND; N];

    let (chunks, remainder) = split_chunks(src, N);

    for syms in remainder.windows(CONTEXT_SIZE).rev() {
        put(&mut states[N - 1], syms[0], syms[1], buf);
    }

    let Some(chunk_size) = chunks.first().map(|chunk| chunk.len()) else {
        return states;
    };

    let chunks: [&[u8]; N] = std::array::from_fn(|i| chunks[i]);

    // The decoder interleaves the states at each position, reading from a single stream, so
    // encode positions in reverse with the states in reverse at each one.
    for k in (0..chunk_size - 1).rev() {
        for (state, chunk) in states.iter_mut().zip(&chunks).rev() {
            put(state, chunk[k], chunk[k + 1], buf);
        }
    }

    for (state, chunk) in states.iter_mut().zip(&chunks).rev() {
        put(state, NUL, chunk[0], buf);
    }

    states
}

fn write_frequencies(
    dst: &mut Vec<u8>,
    alphabet: &[bool; ALPHABET_SIZE],
    frequencies: &Frequencies,
) -> io::Result<()> {
    for (_, fs) in alphabet.iter().zip(frequencies).filter(|(a, _)| **a) {
        let mut iter = alphabet.iter().zip(fs).filter(|(b, _)| **b).peekable();

        while let Some((_, &f)) = iter.next() {
            write_uint7(dst, f)?;

            if f == 0 {
                let mut len = 0;

                while let Some((_, g)) = iter.peek()
                    && **g == 0
                {
                    len += 1;
                    iter.next();
                }

                let n = u8::try_from(len)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                write_u8(dst, n)?;
            }
        }
    }

    Ok(())
}

fn build_alphabet(src: &[u8]) -> [bool; ALPHABET_SIZE] {
    let mut alphabet = [false; ALPHABET_SIZE];
    alphabet[usize::from(NUL)] = true;

    for &b in src {
        let i = usize::from(b);
        alphabet[i] = true;
    }

    alphabet
}

fn build_frequencies(src: &[u8], state_count: usize) -> Frequencies {
    let mut frequencies = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    let chunk_size = src.len() / state_count;

    for chunk in src.chunks_exact(chunk_size).take(state_count) {
        let (i, j) = (usize::from(NUL), usize::from(chunk[0]));
        frequencies[i][j] += 1;
    }

    let rows = <&mut [_; ALPHABET_SIZE]>::try_from(&mut frequencies[..]).expect("a row per symbol");
    histogram::add_order_1(src, rows);

    frequencies
}

/// The contexts that can occur: the symbols in `alphabet` (which includes NUL).
fn contexts(alphabet: &[bool; ALPHABET_SIZE]) -> impl Iterator<Item = usize> + '_ {
    alphabet
        .iter()
        .enumerate()
        .filter(|(_, a)| **a)
        .map(|(i, _)| i)
}

// Only the rows of contexts in the alphabet are touched: a fresh table is zeroed lazily, so
// rows never read cost nothing.

fn build_cumulative_frequencies(
    frequencies: &Frequencies,
    alphabet: &[bool; ALPHABET_SIZE],
) -> Frequencies {
    let mut cumulative_frequencies = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];

    for i in contexts(alphabet) {
        cumulative_frequencies[i] = order_0::build_cumulative_frequencies(&frequencies[i]);
    }

    cumulative_frequencies
}

fn split_chunks(src: &[u8], state_count: usize) -> (Vec<&[u8]>, &[u8]) {
    let (q, r) = (src.len() / state_count, src.len() % state_count);
    let i = q * state_count;

    // SAFETY: `i` <= `src.len()`.
    let chunks = src[..i].chunks(q).collect();
    let remainder = if r > 0 { &src[i - 1..] } else { &[] };

    (chunks, remainder)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        let src = b"abra";
        let frequencies = build_frequencies(src, 4);
        let alphabet = build_alphabet(src);

        let mut dst = Vec::new();
        write_frequencies(&mut dst, &alphabet, &frequencies)?;

        let expected = [
            0x00, // frequencies[NUL][NUL] = 0
            0x00, // run length = 0
            0x02, // frequencies[NUL]['a'] = 2
            0x01, // frequencies[NUL]['b'] = 1
            0x01, // frequencies[NUL]['r'] = 1
            0x00, // frequencies['a'][NUL] = 0
            0x01, // run length = 1
            0x01, // frequencies['a']['b'] = 1
            0x00, // frequencies['a']['r'] = 0
            0x00, // run length = 0
            0x00, // frequencies['b'][NUL] = 0
            0x02, // run length = 2
            0x01, // frequencies['b']['r'] = 1
            0x00, // frequencies['r'][NUL] = 0
            0x00, // run length = 0
            0x01, // frequencies['r']['a'] = 1
            0x00, // frequencies['r']['b'] = 0
            0x01, // run length = 1
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        let src = b"abracadabra";

        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(NUL)] = true;
        expected[usize::from(b'a')] = true;
        expected[usize::from(b'b')] = true;
        expected[usize::from(b'c')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'r')] = true;

        assert_eq!(build_alphabet(src), expected);
    }

    #[test]
    fn test_build_frequencies() {
        let src = b"abracadabra";

        let mut expected = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];
        expected[usize::from(NUL)][usize::from(b'a')] = 1;
        expected[usize::from(NUL)][usize::from(b'c')] = 1;
        expected[usize::from(NUL)][usize::from(b'd')] = 1;
        expected[usize::from(NUL)][usize::from(b'r')] = 1;
        expected[usize::from(b'a')][usize::from(b'b')] = 2;
        expected[usize::from(b'a')][usize::from(b'c')] = 1;
        expected[usize::from(b'a')][usize::from(b'd')] = 1;
        expected[usize::from(b'b')][usize::from(b'r')] = 2;
        expected[usize::from(b'c')][usize::from(b'a')] = 1;
        expected[usize::from(b'd')][usize::from(b'a')] = 1;
        expected[usize::from(b'r')][usize::from(b'a')] = 2;

        assert_eq!(build_frequencies(src, 4), expected);
    }

    #[test]
    fn test_split_chunks() {
        assert_eq!(
            split_chunks(b"abracadabra", 4),
            (
                vec![&b"ab"[..], &b"ra"[..], &b"ca"[..], &b"da"[..]],
                &b"abra"[..]
            )
        );

        assert_eq!(
            split_chunks(b"ndls", 2),
            (vec![&b"nd"[..], &b"ls"[..]], &[][..])
        );
    }
}
