//! Decoding 4-bit BAM bases to their characters, 32 bases at a time with SSSE3 (x86-64, when
//! the CPU has it) or NEON (aarch64), and a byte at a time through a table for the rest.

/// SAM's bases by their 4-bit BAM codes (HTSlib's `seq_nt16_str`).
const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// The two characters of each byte of two 4-bit bases, the first from the high nibble.
static PAIRS: [[u8; 2]; 256] = {
    let mut table = [[0; 2]; 256];
    let mut i = 0;

    while i < 256 {
        table[i] = [BASES[i >> 4], BASES[i & 0x0f]];
        i += 1;
    }

    table
};

/// Decodes the first `dst.len()` 4-bit bases of `src` (two a byte, the first in the high nibble)
/// to their characters (`=ACMGRSVTWYHKDBN`) in `dst`.
///
/// # Panics
///
/// Panics if `src` holds fewer than `dst.len()` bases.
///
/// # Examples
///
/// ```
/// use pasta::bam::record::codec::decoder::decode_bases;
/// let mut dst = [0; 3];
/// decode_bases(&[0x12, 0x40], &mut dst);
/// assert_eq!(&dst, b"ACG");
/// ```
pub fn decode_bases(src: &[u8], dst: &mut [u8]) {
    assert!(
        src.len() >= dst.len().div_ceil(2),
        "{} bytes hold fewer than {} bases",
        src.len(),
        dst.len()
    );

    let done = decode_blocks(src, dst);
    decode_pairs(&src[done / 2..], &mut dst[done..]);
}

/// Decodes `dst.len()` bases a byte at a time. `src` holds at least that many.
fn decode_pairs(src: &[u8], dst: &mut [u8]) {
    let (pairs, rest) = dst.as_chunks_mut::<2>();

    for (pair, &byte) in pairs.iter_mut().zip(src) {
        *pair = PAIRS[usize::from(byte)];
    }

    if let [last] = rest {
        *last = PAIRS[usize::from(src[pairs.len()])][0];
    }
}

/// Decodes as many whole blocks of 32 bases as fit in `dst` with SIMD, returning how many bases
/// it decoded (0 without SIMD). `src` holds at least `dst.len()` bases.
#[cfg(target_arch = "x86_64")]
fn decode_blocks(src: &[u8], dst: &mut [u8]) -> usize {
    if dst.len() < BLOCK || !std::arch::is_x86_feature_detected!("ssse3") {
        return 0;
    }

    // SAFETY: SSSE3 is available.
    unsafe { decode_blocks_ssse3(src, dst) }
}

#[cfg(target_arch = "aarch64")]
fn decode_blocks(src: &[u8], dst: &mut [u8]) -> usize {
    // SAFETY: NEON is part of the aarch64 baseline.
    unsafe { decode_blocks_neon(src, dst) }
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
fn decode_blocks(_: &[u8], _: &mut [u8]) -> usize {
    0
}

/// The bases a SIMD block decodes: 16 bytes of input.
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
const BLOCK: usize = 32;

/// Splits each byte into its two nibbles, looks each up in [`BASES`] with a byte shuffle, and
/// interleaves the high (first) and low (second) characters.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn decode_blocks_ssse3(src: &[u8], dst: &mut [u8]) -> usize {
    use std::arch::x86_64::{
        _mm_and_si128, _mm_loadu_si128, _mm_set1_epi8, _mm_shuffle_epi8, _mm_srli_epi16,
        _mm_storeu_si128, _mm_unpackhi_epi8, _mm_unpacklo_epi8,
    };

    let blocks = dst.len() / BLOCK;
    // SAFETY: `BASES` has 16 bytes.
    let table = unsafe { _mm_loadu_si128(BASES.as_ptr().cast()) };
    let low_nibbles = _mm_set1_epi8(0x0f);

    // `src` has at least as many 16-byte chunks as `dst` has blocks.
    let (src, _) = src.as_chunks::<{ BLOCK / 2 }>();
    let (dst, _) = dst.as_chunks_mut::<BLOCK>();

    for (src, dst) in src.iter().zip(dst) {
        // SAFETY: `src` has 16 bytes.
        let packed = unsafe { _mm_loadu_si128(src.as_ptr().cast()) };
        let first = _mm_shuffle_epi8(table, _mm_and_si128(_mm_srli_epi16(packed, 4), low_nibbles));
        let second = _mm_shuffle_epi8(table, _mm_and_si128(packed, low_nibbles));

        // SAFETY: `dst` has 32 bytes.
        unsafe {
            _mm_storeu_si128(dst.as_mut_ptr().cast(), _mm_unpacklo_epi8(first, second));
            _mm_storeu_si128(
                dst[16..].as_mut_ptr().cast(),
                _mm_unpackhi_epi8(first, second),
            );
        }
    }

    blocks * BLOCK
}

/// As [`decode_blocks_ssse3`], with a table lookup per nibble and an interleaving store.
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn decode_blocks_neon(src: &[u8], dst: &mut [u8]) -> usize {
    use std::arch::aarch64::{
        uint8x16x2_t, vandq_u8, vdupq_n_u8, vld1q_u8, vqtbl1q_u8, vshrq_n_u8, vst2q_u8,
    };

    let blocks = dst.len() / BLOCK;
    // SAFETY: `BASES` has 16 bytes.
    let table = unsafe { vld1q_u8(BASES.as_ptr()) };
    let low_nibbles = vdupq_n_u8(0x0f);

    // `src` has at least as many 16-byte chunks as `dst` has blocks.
    let (src, _) = src.as_chunks::<{ BLOCK / 2 }>();
    let (dst, _) = dst.as_chunks_mut::<BLOCK>();

    for (src, dst) in src.iter().zip(dst) {
        // SAFETY: `src` has 16 bytes.
        let packed = unsafe { vld1q_u8(src.as_ptr()) };
        let first = vqtbl1q_u8(table, vshrq_n_u8::<4>(packed));
        let second = vqtbl1q_u8(table, vandq_u8(packed, low_nibbles));

        // SAFETY: `dst` has 32 bytes.
        unsafe { vst2q_u8(dst.as_mut_ptr(), uint8x16x2_t(first, second)) };
    }

    blocks * BLOCK
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decode_each(src: &[u8], len: usize) -> Vec<u8> {
        (0..len)
            .map(|i| BASES[usize::from((src[i / 2] >> ((!i & 1) << 2)) & 0x0f)])
            .collect()
    }

    #[test]
    fn test_decode_bases_matches_decoding_each_base() {
        // Every byte value, at every offset of a SIMD block, and lengths around block edges.
        let src: Vec<u8> = (0..=255).chain((0..=255).rev()).chain(0..=63).collect();

        for start in 0..17 {
            let src = &src[start..];

            for len in (0..=100).chain([255, 256, 257, 511, 512, 513, 2 * src.len() - 1]) {
                let mut dst = vec![0; len];
                decode_bases(src, &mut dst);
                assert_eq!(dst, decode_each(src, len), "start {start}, length {len}");
            }
        }
    }

    #[test]
    fn test_decode_pairs_matches_decoding_each_base() {
        let src: Vec<u8> = (0..=255).collect();

        for len in [0, 1, 2, 31, 32, 33, 511, 512] {
            let mut dst = vec![0; len];
            decode_pairs(&src, &mut dst);
            assert_eq!(dst, decode_each(&src, len), "length {len}");
        }
    }

    #[test]
    #[should_panic(expected = "1 bytes hold fewer than 3 bases")]
    fn test_decode_bases_with_too_few_bytes() {
        decode_bases(&[0x12], &mut [0; 3]);
    }
}
