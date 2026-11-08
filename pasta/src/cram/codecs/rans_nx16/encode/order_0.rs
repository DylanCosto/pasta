use std::io;

use super::{EncodingSymbol, LOWER_BOUND, ReverseWriter, write_alphabet, write_states};
use crate::cram::{
    codecs::{histogram, rans_nx16::ALPHABET_SIZE},
    io::writer::num::write_uint7,
};

const NORMALIZATION_BITS: u32 = 12;

pub(super) type Frequencies = [u32; ALPHABET_SIZE];

pub(super) struct Context {
    /// The frequencies as written: normalized to the power of two at or above the input's length,
    /// up to `1 << NORMALIZATION_BITS`, so that a short input has a short table (as htscodecs
    /// writes them). Decoders scale them up.
    table: Frequencies,
    /// The frequencies encoded with, summing to `1 << NORMALIZATION_BITS`.
    frequencies: Frequencies,
}

pub(super) fn build_context(src: &[u8]) -> io::Result<Context> {
    let mut table = build_frequencies(src);

    let len =
        u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let total = len.next_power_of_two().min(1 << NORMALIZATION_BITS);
    normalize_frequencies(&mut table, len, total)?;

    let mut frequencies = table;
    scale_frequencies(&mut frequencies, total, 1 << NORMALIZATION_BITS);

    Ok(Context { table, frequencies })
}

pub(super) fn write_context(dst: &mut Vec<u8>, ctx: &Context) -> io::Result<()> {
    let alphabet = build_alphabet(&ctx.table);
    write_alphabet(dst, &alphabet)?;

    write_frequencies(dst, &ctx.table)?;

    Ok(())
}

pub fn encode(src: &[u8], ctx: &Context, state_count: usize, dst: &mut Vec<u8>) -> io::Result<()> {
    let frequencies = &ctx.frequencies;
    let cumulative_frequencies = build_cumulative_frequencies(frequencies);

    let mut symbols = [EncodingSymbol::default(); ALPHABET_SIZE];

    for ((symbol, &f), &g) in symbols
        .iter_mut()
        .zip(frequencies)
        .zip(&cumulative_frequencies)
    {
        if f > 0 {
            *symbol = EncodingSymbol::new(f, g, NORMALIZATION_BITS);
        }
    }

    let mut buf = ReverseWriter::for_symbols(src.len());

    match state_count {
        4 => write_states(dst, &encode_states::<4>(src, &symbols, &mut buf))?,
        32 => write_states(dst, &encode_states::<32>(src, &symbols, &mut buf))?,
        _ => unreachable!("state count is 4 or 32"),
    }

    dst.extend(buf.as_slice());

    Ok(())
}

/// Encodes `src` on `N` states, symbol `i` on state `i % N`, returning the final states.
fn encode_states<const N: usize>(
    src: &[u8],
    symbols: &[EncodingSymbol; ALPHABET_SIZE],
    buf: &mut ReverseWriter,
) -> [u32; N] {
    let mut states = [LOWER_BOUND; N];

    // Encoding goes backward: the last, partial round first, then whole rounds, each with its
    // states in reverse.
    let (rounds, last_round) = src.as_chunks::<N>();

    for (state, &sym) in states.iter_mut().zip(last_round).rev() {
        *state = symbols[usize::from(sym)].put(*state, buf);
    }

    for round in rounds.iter().rev() {
        for (state, &sym) in states.iter_mut().zip(round).rev() {
            *state = symbols[usize::from(sym)].put(*state, buf);
        }
    }

    states
}

fn write_frequencies(dst: &mut Vec<u8>, frequencies: &Frequencies) -> io::Result<()> {
    for &f in frequencies {
        if f > 0 {
            write_uint7(dst, f)?;
        }
    }

    Ok(())
}

fn build_alphabet(frequencies: &Frequencies) -> [bool; ALPHABET_SIZE] {
    frequencies.map(|f| f > 0)
}

fn build_frequencies(src: &[u8]) -> Frequencies {
    histogram::order_0(src)
}

/// Scales `frequencies`, which sum to `sum`, to sum to `total`, rounding each to the nearest and
/// keeping it at least 1, then moving what is over or under to the most frequent symbol; if that
/// would take more than half of it, the scaling is done once more from the rounded frequencies
/// (htscodecs' `normalise_freq`, which this matches exactly, so that tables are the same size).
pub(super) fn normalize_frequencies(
    frequencies: &mut Frequencies,
    sum: u32,
    total: u32,
) -> io::Result<()> {
    if sum == 0 {
        return Ok(());
    }

    let mut size = u64::from(sum);
    let mut retried = false;

    loop {
        let ratio = (u64::from(total) << 31) / size + (1 << 30) / size;
        let (mut max, mut max_index, mut new_size) = (0, 0, 0);

        for (i, f) in frequencies.iter_mut().enumerate() {
            if *f == 0 {
                continue;
            }

            if max < *f {
                (max, max_index) = (*f, i);
            }

            // At most `total` (a `u32`) once scaled.
            *f = (((u64::from(*f) * ratio) >> 31) as u32).max(1);
            new_size += u64::from(*f);
        }

        let adjust = i64::from(total) - new_size as i64;
        let f_max = i64::from(frequencies[max_index]);

        if adjust > 0 {
            frequencies[max_index] += adjust as u32;
        } else if adjust < 0 {
            if f_max > -adjust && (retried || f_max / 2 >= -adjust) {
                frequencies[max_index] = (f_max + adjust) as u32;
            } else if !retried {
                retried = true;
                size = new_size;
                continue;
            } else {
                // Take the excess from the other symbols, keeping each at least 1.
                let mut adjust = adjust + f_max - 1;
                frequencies[max_index] = 1;

                for f in frequencies.iter_mut() {
                    if adjust == 0 {
                        break;
                    }

                    if *f < 2 {
                        continue;
                    }

                    let step = if i64::from(*f) > -adjust {
                        adjust
                    } else {
                        1 - i64::from(*f)
                    };

                    *f = (i64::from(*f) + step) as u32;
                    adjust -= step;
                }
            }
        }

        return if frequencies[max_index] > 0 {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "cannot normalize frequencies",
            ))
        };
    }
}

/// Scales `frequencies` from summing to `total` to summing to `new_total`, both powers of two
/// (htscodecs' `normalise_freq_shift`).
pub(super) fn scale_frequencies(frequencies: &mut Frequencies, total: u32, new_total: u32) {
    if total == 0 || total >= new_total {
        return;
    }

    let shift = new_total.trailing_zeros() - total.trailing_zeros();

    for f in frequencies {
        *f <<= shift;
    }
}

pub(super) fn build_cumulative_frequencies(frequencies: &Frequencies) -> [u32; ALPHABET_SIZE] {
    let mut cumulative_frequencies = [0; ALPHABET_SIZE];
    let mut f = cumulative_frequencies[0];

    for (next_f, g) in cumulative_frequencies[1..].iter_mut().zip(frequencies) {
        *next_f = f + g;
        f = *next_f;
    }

    cumulative_frequencies
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        let src = b"abracadabra";
        let frequencies = build_frequencies(src);

        let mut dst = Vec::new();
        write_frequencies(&mut dst, &frequencies)?;

        let expected = [
            0x05, // frequencies[b'a'] = 5
            0x02, // frequencies[b'b'] = 2
            0x01, // frequencies[b'c'] = 1
            0x01, // frequencies[b'd'] = 1
            0x02, // frequencies[b'r'] = 2
        ];

        assert_eq!(dst, expected);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        let src = b"abracadabra";
        let frequencies = build_frequencies(src);

        let mut expected = [false; ALPHABET_SIZE];

        for &sym in src {
            expected[usize::from(sym)] = true;
        }

        assert_eq!(build_alphabet(&frequencies), expected);
    }

    #[test]
    fn test_build_frequencies() {
        let src = b"abracadabra";

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(build_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() -> io::Result<()> {
        // The values htscodecs' `normalise_freq` gives.
        fn t(total: u32, expected: [u32; 5]) -> io::Result<()> {
            let mut frequencies = build_frequencies(b"abracadabra");
            normalize_frequencies(&mut frequencies, 11, total)?;
            let actual = b"abcdr".map(|sym| frequencies[usize::from(sym)]);
            assert_eq!(actual, expected, "total {total}");
            Ok(())
        }

        t(16, [10, 2, 1, 1, 2])?;
        t(4096, [1864, 744, 372, 372, 744])?;

        // Too many rare symbols for the most frequent one to make up the difference: the
        // scaling is redone, then the rest is taken from the other symbols.
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[0] = 1000;
        frequencies[1..40].fill(1);
        normalize_frequencies(&mut frequencies, 1039, 64)?;
        assert_eq!(frequencies[0], 25);
        assert_eq!(frequencies.iter().sum::<u32>(), 64);

        Ok(())
    }

    #[test]
    fn test_build_context_writes_a_table_scaled_to_the_input() -> io::Result<()> {
        let ctx = build_context(b"abracadabra")?;
        assert_eq!(ctx.table.iter().sum::<u32>(), 16);
        assert_eq!(ctx.frequencies.iter().sum::<u32>(), 4096);
        assert_eq!(ctx.frequencies[usize::from(b'a')], 10 << 8);
        Ok(())
    }
}
