#[cfg(target_arch = "x86_64")]
mod avx2;

use std::io;

use super::{SymbolTable, read_alphabet, read_states, renormalize_states, state_renormalize};
use crate::cram::{codecs::rans_nx16::ALPHABET_SIZE, io::reader::num::read_uint7};

const NORMALIZATION_BITS: u32 = 12;

/// The fewest whole rounds decoded with SIMD. Its table takes time to pack that a shorter
/// output does not repay: measured, 2 kB decode about as fast either way, 512 bytes slower with
/// SIMD and 8 kB faster.
#[cfg(target_arch = "x86_64")]
const SIMD_MIN_ROUNDS: usize = 128;

pub fn decode(src: &mut &[u8], dst: &mut [u8], state_count: usize) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    let table = SymbolTable::new(&frequencies, NORMALIZATION_BITS)?;

    match state_count {
        4 => decode_with_states::<4, false>(src, dst, &table),
        32 => decode_with_states::<32, true>(src, dst, &table),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid rANS Nx16 state count",
        )),
    }
}

/// Decodes with `N` states, using SIMD where the CPU supports it if `SIMD` is set.
fn decode_with_states<const N: usize, const SIMD: bool>(
    src: &mut &[u8],
    dst: &mut [u8],
    table: &SymbolTable,
) -> io::Result<()> {
    let mut states = read_states::<N>(src)?;

    let (chunks, tail) = dst.as_chunks_mut::<N>();

    let simd_decoded =
        SIMD && !chunks.is_empty() && try_decode_simd(table, &mut states, src, chunks)?;

    if !simd_decoded {
        for chunk in chunks {
            for (d, state) in chunk.iter_mut().zip(&mut states) {
                *d = table.decode_symbol(state, NORMALIZATION_BITS);
            }

            // The state steps are independent, so they run ahead of the ordered renormalization
            // reads.
            renormalize_states(&mut states, src)?;
        }
    }

    for (d, state) in tail.iter_mut().zip(&mut states) {
        *d = table.decode_symbol(state, NORMALIZATION_BITS);
        *state = state_renormalize(*state, src)?;
    }

    Ok(())
}

/// Decodes every whole round with SIMD, if the CPU supports it for `N` states and the table
/// allows, and returns whether it did.
#[cfg(target_arch = "x86_64")]
fn try_decode_simd<const N: usize>(
    table: &SymbolTable,
    states: &mut [u32; N],
    src: &mut &[u8],
    chunks: &mut [[u8; N]],
) -> io::Result<bool> {
    let Ok(states) = <&mut [u32; 32]>::try_from(states.as_mut_slice()) else {
        return Ok(false);
    };

    if chunks.len() < SIMD_MIN_ROUNDS || !std::arch::is_x86_feature_detected!("avx2") {
        return Ok(false);
    }

    let Some(entries) = avx2::Entries::new(table) else {
        return Ok(false);
    };

    // `N` is 32, so these are the same chunks.
    let (rounds, rest) = chunks.as_flattened_mut().as_chunks_mut::<32>();
    debug_assert!(rest.is_empty());

    // SAFETY: AVX2 is available.
    unsafe { avx2::decode(&entries, states, src, rounds)? };

    Ok(true)
}

#[cfg(not(target_arch = "x86_64"))]
fn try_decode_simd<const N: usize>(
    _table: &SymbolTable,
    _states: &mut [u32; N],
    _src: &mut &[u8],
    _chunks: &mut [[u8; N]],
) -> io::Result<bool> {
    Ok(false)
}

pub fn normalize_frequencies(frequencies: &mut [u32; ALPHABET_SIZE], bits: u32) {
    // Sum in u64: a malformed table must not overflow before the symbol table rejects it.
    let mut sum: u64 = frequencies.iter().copied().map(u64::from).sum();

    if sum == 0 || sum >= (1 << bits) {
        return;
    }

    let mut shift = 0;

    while sum < (1 << bits) {
        sum *= 2;
        shift += 1;
    }

    for f in frequencies {
        *f <<= shift;
    }
}

fn read_frequencies(src: &mut &[u8]) -> io::Result<[u32; ALPHABET_SIZE]> {
    let alphabet = read_alphabet(src)?;

    let mut frequencies = [0; ALPHABET_SIZE];

    for (i, frequency) in alphabet.iter().zip(&mut frequencies) {
        if *i {
            *frequency = read_uint7(src)?;
        }
    }

    normalize_frequencies(&mut frequencies, NORMALIZATION_BITS);

    Ok(frequencies)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cram::{
        codecs::rans_nx16::{Flags, encode::encode},
        io::reader::num::{read_u8, read_uint7_as},
    };

    #[test]
    fn test_decode_round_trips_32_states() -> io::Result<()> {
        // Quality-like runs, uniform bytes that renormalize nearly every state every round, and
        // one symbol, whose frequency is the whole table.
        let mut state = 1u32;
        let mut next = move || {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            state >> 8
        };

        let qualities: Vec<u8> = (0..200_003)
            .scan(30u8, |q, _| {
                let r = next();
                if r % 7 == 0 {
                    *q = 2 + (r % 40) as u8;
                }
                Some(*q + 33)
            })
            .collect();
        let uniform: Vec<u8> = (0..100_001).map(|_| next() as u8).collect();
        let constant = vec![b'#'; 5_000];

        for src in [&qualities[..], &uniform[..], &constant[..]] {
            for len in [31, 32, 33, 63, 64, 65, 1_000, 4_097, src.len()] {
                let src = &src[..len];
                let encoded = encode(Flags::N32, src)?;
                assert_eq!(super::super::decode(&encoded, len)?, src, "len = {len}");

                // Bytes that code no smaller are stored as they are.
                if Flags::from(encoded[0]).contains(Flags::CAT) {
                    continue;
                }

                // With and without SIMD, the output and what is left of the input are the same.
                let (scalar, scalar_rest) = decode_32::<false>(&encoded, len)?;
                let (simd, simd_rest) = decode_32::<true>(&encoded, len)?;
                assert_eq!(scalar, src, "len = {len}");
                assert_eq!(simd, src, "len = {len}");
                assert_eq!(simd_rest, scalar_rest, "len = {len}");

                // The decoder reads the stream to its end, so a truncated one runs out of input.
                for cut in [1, 2, 3, 70] {
                    let truncated = &encoded[..encoded.len().saturating_sub(cut)];

                    for result in [
                        decode_32::<false>(truncated, len),
                        decode_32::<true>(truncated, len),
                    ] {
                        assert!(
                            result.is_err_and(|e| e.kind() == io::ErrorKind::UnexpectedEof),
                            "len = {len}, cut = {cut}"
                        );
                    }
                }
            }
        }

        Ok(())
    }

    /// A table whose frequencies do not fill it leaves slots with a frequency of 0, which the
    /// SIMD table cannot hold: such a stream is decoded without SIMD, to the same bytes.
    #[test]
    fn test_decode_with_a_table_simd_cannot_pack() -> io::Result<()> {
        let frequencies = [0; ALPHABET_SIZE];
        let table = SymbolTable::new(&frequencies, NORMALIZATION_BITS)?;

        #[cfg(target_arch = "x86_64")]
        assert!(avx2::Entries::new(&table).is_none());

        let mut src = Vec::new();

        for i in 0..32u32 {
            src.extend((0x0001_0000 + i * 0x0101).to_le_bytes());
        }

        src.extend([0xa5; 512]);

        let decode = |simd: bool| {
            let mut dst = [0; 100];
            let mut input = &src[..];

            if simd {
                decode_with_states::<32, true>(&mut input, &mut dst, &table)?;
            } else {
                decode_with_states::<32, false>(&mut input, &mut dst, &table)?;
            }

            Ok::<_, io::Error>((dst, input.len()))
        };

        assert_eq!(decode(true)?, decode(false)?);

        Ok(())
    }

    /// Decodes a 32-state order-0 stream with or without SIMD, returning the bytes and the
    /// length of the input left.
    fn decode_32<const SIMD: bool>(mut src: &[u8], len: usize) -> io::Result<(Vec<u8>, usize)> {
        let flags = read_u8(&mut src)?;
        assert_eq!(Flags::from(flags), Flags::N32);
        let size: usize = read_uint7_as(&mut src)?;
        assert_eq!(size, len);

        let frequencies = read_frequencies(&mut src)?;
        let table = SymbolTable::new(&frequencies, NORMALIZATION_BITS)?;
        let mut dst = vec![0; len];
        decode_with_states::<32, SIMD>(&mut src, &mut dst, &table)?;

        Ok((dst, src.len()))
    }

    #[test]
    fn test_decode_rejects_oversized_frequency_tables() {
        // Alphabet {'a', 'c'}; 'a' = 4095 and 'c' = f, with no states following: a total of 4096
        // is valid (and reading the states fails), but 4097 is not.
        for (f, expected) in [
            (0x01, io::ErrorKind::UnexpectedEof),
            (0x02, io::ErrorKind::InvalidData),
        ] {
            let src = [b'a', b'c', 0x00, 0x9f, 0x7f, f];
            let mut dst = [0; 8];
            let actual = decode(&mut &src[..], &mut dst, 4).unwrap_err().kind();
            assert_eq!(actual, expected, "f = {f}");
        }
    }
}
