use std::io::{self, Write};

use super::{EncodingSymbol, ReverseWriter, header, write_header, write_states};
use crate::cram::{
    codecs::{
        histogram,
        rans_4x8::{ALPHABET_SIZE, LOWER_BOUND, Order, STATE_COUNT},
    },
    io::writer::num::{write_itf8, write_u8},
};

type RawFrequencies = [u32; ALPHABET_SIZE];
type Frequencies = [u16; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    let raw_frequencies = build_raw_frequencies(src);
    let frequencies = normalize_frequencies(&raw_frequencies);
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut symbols = [EncodingSymbol::default(); ALPHABET_SIZE];

    for ((symbol, &f), &g) in symbols
        .iter_mut()
        .zip(&frequencies)
        .zip(&cumulative_frequencies)
    {
        if f > 0 {
            *symbol = EncodingSymbol::new(f, g);
        }
    }

    let mut buf = ReverseWriter::for_symbols(src.len());
    let mut states = [LOWER_BOUND; STATE_COUNT];

    for (i, &sym) in src.iter().enumerate().rev() {
        let state = &mut states[i % STATE_COUNT];
        *state = symbols[usize::from(sym)].put(*state, &mut buf);
    }

    let mut dst = vec![0; header::SIZE];

    write_frequencies(&mut dst, &frequencies)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.as_slice());

    let compressed_size = dst[header::SIZE..].len();
    // SAFETY: `dst.len() >= header::SIZE`.
    let header_dst = dst.first_chunk_mut().unwrap();
    write_header(header_dst, Order::Zero, compressed_size, src.len())?;

    Ok(dst)
}

pub fn write_frequencies<W>(writer: &mut W, frequencies: &Frequencies) -> io::Result<()>
where
    W: Write,
{
    const NUL: u8 = 0x00;

    let mut iter = frequencies.iter().enumerate();
    let mut prev_sym = None;

    while let Some((sym, &f)) = iter.next() {
        if f == 0 {
            continue;
        }

        // SAFETY: `sym < ALPHABET_SIZE`.
        write_u8(writer, sym as u8)?;

        if sym > 0 && prev_sym == Some(sym - 1) {
            let i = sym + 1;

            let len = frequencies[i..]
                .iter()
                .position(|&g| g == 0)
                .unwrap_or(ALPHABET_SIZE - i);

            // SAFETY: `len < ALPHABET_SIZE`.
            write_u8(writer, len as u8)?;

            write_itf8(writer, i32::from(f))?;

            for (sym, &g) in iter.by_ref().take(len) {
                write_itf8(writer, i32::from(g))?;
                prev_sym = Some(sym);
            }

            continue;
        }

        write_itf8(writer, i32::from(f))?;

        prev_sym = Some(sym);
    }

    write_u8(writer, NUL)?;

    Ok(())
}

fn build_raw_frequencies(src: &[u8]) -> RawFrequencies {
    histogram::order_0(src)
}

/// The total that normalized frequencies are scaled to (htscodecs' `TOTFREQ`). They add up to one
/// less, 4095 (§ 2.1 "Frequency table").
pub(super) const TOTAL: u32 = 1 << 12;

/// 0.98 in 1.31 fixed point (htscodecs').
const SCALE_98: u64 = 2_104_533_975;

/// Normalizes `raw_frequencies` to add up to 4095 as htscodecs' `rans_compress_O0` does, so that
/// the output is the same: each is scaled by [`TOTAL`] over their sum in 1.31 fixed point, rounded
/// down but kept above 0, and the most frequent symbol (the first, if several) takes up the
/// difference, unless that would take more than half of it, when all are scaled by 0.98 first.
pub(super) fn normalize_frequencies(raw_frequencies: &RawFrequencies) -> Frequencies {
    let len: u64 = raw_frequencies.iter().map(|&f| u64::from(f)).sum();

    if len == 0 {
        return [0; ALPHABET_SIZE];
    }

    let mut scale = (u64::from(TOTAL) << 31) / len + (1 << 30) / len;
    let mut frequencies = *raw_frequencies;

    loop {
        let (mut sum, mut max, mut max_index) = (1, 0, 0);

        for (i, f) in frequencies.iter_mut().enumerate().filter(|(_, f)| **f > 0) {
            if *f > max {
                (max, max_index) = (*f, i);
            }

            // At most `TOTAL`.
            *f = ((u64::from(*f) * scale) >> 31).max(1) as u32;
            sum += *f;
        }

        let f = &mut frequencies[max_index];

        if sum < TOTAL {
            *f += TOTAL - sum;
        } else if sum - TOTAL > *f / 2 {
            scale = SCALE_98;
            continue;
        } else {
            *f -= sum - TOTAL;
        }

        // Each is at most `TOTAL` - 1.
        return frequencies.map(|f| f as u16);
    }
}

pub(super) fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
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
    fn test_encode() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-0 encoding (2020-07-22)
        let data = b"abracadabra";
        let actual = encode(data)?;

        let expected = [
            // header
            0x00, // order = 0
            0x1f, 0x00, 0x00, 0x00, // compressed size = 31
            0x0b, 0x00, 0x00, 0x00, // uncompressed size = 11
            // frequency table
            0x61, // sym = 'a'
            0x87, 0x47, // f['a'] = 1863
            0x62, // sym = 'b'
            0x02, // rle = 2
            0x82, 0xe8, // f['b'] = 744
            0x81, 0x74, // f['c'] = 372
            0x81, 0x74, // f['d'] = 372
            0x72, // 'r'
            0x82, 0xe8, // f['r'] = 744
            0x00, // end
            // compressed blob
            0xd2, 0x02, 0xa4, 0x42, // states[0]
            0x0d, 0x3a, 0x52, 0x21, // states[1]
            0xd0, 0xfe, 0xa1, 0x42, // states[2]
            0x40, 0xa6, 0x6a, 0x02, // states[3]
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        fn t(src: &[u8], expected: &[u8]) -> io::Result<()> {
            let raw_frequencies = build_raw_frequencies(src);
            let frequencies = normalize_frequencies(&raw_frequencies);

            let mut actual = Vec::new();
            write_frequencies(&mut actual, &frequencies)?;

            assert_eq!(actual, expected);

            Ok(())
        }

        // § 2.1.1 "Frequency table: Order-0 encoding" (2023-03-15)
        t(
            b"abracadabra",
            &[
                b'a', // sym = 'a'
                0x87, 0x47, // f['a'] = 1863
                b'b', // sym = 'b'
                0x02, // rle = 2
                0x82, 0xe8, // f['b'] = 744
                0x81, 0x74, // f['c'] = 372
                0x81, 0x74, // f['d'] = 372
                b'r', // 'r'
                0x82, 0xe8, // f['r'] = 744
                0x00, // end
            ],
        )?;

        t(
            &[0x01],
            &[
                0x01, // sym = 0x01
                0x8f, 0xff, // f[0x01] = 4095
                0x00, // end
            ],
        )?;

        t(
            &[0xfd, 0xfe, 0xff],
            &[
                0xfd, // sym = 0xfd
                0x85, 0x55, // f[0xfd] = 1365
                0xfe, // sym = 0xfe
                0x01, // rle = 1
                0x85, 0x55, // f[0xfe] = 1365
                0x85, 0x55, // f[0xff] = 1365
                0x00, // end
            ],
        )?;

        Ok(())
    }

    #[test]
    fn test_build_raw_frequencies() {
        // § 2.1.1 "Frequency table: Order-0 encoding" (2023-03-15)
        let src = b"abracadabra";

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(build_raw_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() {
        let mut raw_frequencies = [0; ALPHABET_SIZE];
        raw_frequencies[usize::from(b'a')] = 5;
        raw_frequencies[usize::from(b'b')] = 2;
        raw_frequencies[usize::from(b'c')] = 1;
        raw_frequencies[usize::from(b'd')] = 1;
        raw_frequencies[usize::from(b'r')] = 2;

        let actual = normalize_frequencies(&raw_frequencies);

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 1863;
        expected[usize::from(b'b')] = 744;
        expected[usize::from(b'c')] = 372;
        expected[usize::from(b'd')] = 372;
        expected[usize::from(b'r')] = 744;

        assert_eq!(actual, expected);

        let raw_frequencies = [0; ALPHABET_SIZE];
        assert_eq!(normalize_frequencies(&raw_frequencies), [0; ALPHABET_SIZE]);
    }
}
