pub mod context;

pub use self::context::{Context, build_context, write_context};
use crate::cram::codecs::rans_nx16::ALPHABET_SIZE;

pub fn encode(src: &[u8], ctx: &Context) -> Vec<u8> {
    let mapping_table = &ctx.mapping_table;

    match ctx.symbol_count.get() {
        1 => Vec::new(),
        2 => pack::<8>(src, mapping_table),
        3..=4 => pack::<4>(src, mapping_table),
        5..=16 => pack::<2>(src, mapping_table),
        _ => unreachable!(),
    }
}

/// Packs `N` symbols, mapped to `8 / N` bits each, into each byte, the first in the low bits.
fn pack<const N: usize>(src: &[u8], mapping_table: &[u8; ALPHABET_SIZE]) -> Vec<u8> {
    const BITS: usize = u8::BITS as usize;

    let pack_chunk = |chunk: &[u8]| {
        chunk.iter().enumerate().fold(0, |d, (i, &sym)| {
            d | (mapping_table[usize::from(sym)] << ((BITS / N) * i))
        })
    };

    let (chunks, remainder) = src.as_chunks::<N>();
    let mut dst = Vec::with_capacity(src.len().div_ceil(N));
    dst.extend(chunks.iter().map(|chunk| pack_chunk(chunk)));

    if !remainder.is_empty() {
        dst.push(pack_chunk(remainder));
    }

    dst
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_8() {
        let src = *b"nndd";

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'n')] = 0;
        mapping_table[usize::from(b'd')] = 1;

        assert_eq!(pack::<8>(&src, &mapping_table), [0b00001100]);
    }

    #[test]
    fn test_pack_4() {
        let src = *b"ndls";

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'n')] = 0;
        mapping_table[usize::from(b'd')] = 1;
        mapping_table[usize::from(b'l')] = 2;
        mapping_table[usize::from(b's')] = 3;

        assert_eq!(pack::<4>(&src, &mapping_table), [0b11100100]);
    }

    #[test]
    fn test_pack_2() {
        let src = *b"nd";

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'n')] = 0;
        mapping_table[usize::from(b'o')] = 1;
        mapping_table[usize::from(b'd')] = 2;
        mapping_table[usize::from(b'l')] = 3;
        mapping_table[usize::from(b's')] = 4;

        assert_eq!(pack::<2>(&src, &mapping_table), [0b00100000]);
    }
}
