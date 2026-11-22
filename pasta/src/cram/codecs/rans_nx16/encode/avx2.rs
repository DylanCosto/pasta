//! AVX2 encoding of 32-state rANS Nx16 streams.
//!
//! A round encodes one symbol per state, 8 states per vector: the symbols' steps are loaded and
//! turned from a step per state into a vector per field, the states that renormalize put their
//! low words out in state order, and all advance together. This is the same computation as the
//! scalar loop, so the two can hand off between rounds.

use std::{
    arch::x86_64::{
        __m128i, __m256i, _mm_loadu_si128, _mm_storeu_si128, _mm256_add_epi32, _mm256_and_si256,
        _mm256_blend_epi32, _mm256_blendv_epi8, _mm256_castsi256_ps, _mm256_castsi256_si128,
        _mm256_cmpgt_epi32, _mm256_loadu_si256, _mm256_movemask_ps, _mm256_mul_epu32,
        _mm256_mullo_epi32, _mm256_packus_epi32, _mm256_permute4x64_epi64,
        _mm256_permutevar8x32_epi32, _mm256_set_m128i, _mm256_set1_epi32, _mm256_srli_epi32,
        _mm256_srli_epi64, _mm256_srlv_epi32, _mm256_storeu_si256, _mm256_sub_epi32,
        _mm256_unpackhi_epi32, _mm256_unpackhi_epi64, _mm256_unpacklo_epi32, _mm256_unpacklo_epi64,
    },
    ptr,
};

use super::{EncodingSymbol, ReverseWriter};
use crate::cram::codecs::rans_nx16::ALPHABET_SIZE;

pub(super) const STATE_COUNT: usize = 32;
const LANES: usize = 8;
const VECTORS: usize = STATE_COUNT / LANES;

/// A step is loaded as one 128-bit value: `x_max`, `rcp_freq`, `bias`, then `cmpl_freq` and
/// `rcp_shift` in the halves of the last 32 bits.
const _: () = assert!(size_of::<EncodingSymbol>() == size_of::<__m128i>());

/// For each 8-bit mask of lanes that renormalize, the lanes in order, at the end of the vector:
/// with `n` of them, lane `LANES - n + i` takes the `i`th. The lanes before them take lane 0,
/// unused.
static COMPRESS: [[u32; LANES]; 256] = build_compress();

const fn build_compress() -> [[u32; LANES]; 256] {
    let mut table = [[0; LANES]; 256];
    let mut mask = 0usize;

    while mask < 256 {
        let mut i = LANES - mask.count_ones() as usize;
        let mut lane = 0;

        while lane < LANES {
            if mask & (1 << lane) != 0 {
                table[mask][i] = lane as u32;
                i += 1;
            }

            lane += 1;
        }

        mask += 1;
    }

    table
}

/// Encodes every whole round of `src` (symbol `i` on state `i % 32`), last round first, as the
/// scalar loop does.
///
/// # Safety
///
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(super) unsafe fn encode_order_0(
    rounds: &[[u8; STATE_COUNT]],
    symbols: &[EncodingSymbol; ALPHABET_SIZE],
    states: &mut [u32; STATE_COUNT],
    buf: &mut ReverseWriter,
) {
    let mut rounds = rounds.iter().rev();

    let next_round = || {
        let round = rounds.next()?;
        Some(std::array::from_fn(|j| &symbols[usize::from(round[j])]))
    };

    // SAFETY: AVX2 is available, as this function requires.
    unsafe { put_rounds(states, buf, next_round) };
}

/// Encodes every symbol of the 32 `chunks` but each one's first (state `j` encodes `chunks[j]`,
/// each symbol in the context of the one before it), last symbols first, as the scalar loop
/// does.
///
/// # Safety
///
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(super) unsafe fn encode_order_1(
    chunks: &[&[u8]; STATE_COUNT],
    rows: &[&[EncodingSymbol; ALPHABET_SIZE]; ALPHABET_SIZE],
    states: &mut [u32; STATE_COUNT],
    buf: &mut ReverseWriter,
) {
    let chunk_size = chunks[0].len();
    assert!(chunks.iter().all(|chunk| chunk.len() == chunk_size));

    let Some(last) = chunk_size.checked_sub(1) else {
        return;
    };

    // The symbol after the one each state encodes next: its step is found in that one's row.
    let mut next: [u8; STATE_COUNT] = std::array::from_fn(|j| chunks[j][last]);
    let mut positions = (0..last).rev();

    let next_round = || {
        let k = positions.next()?;

        Some(std::array::from_fn(|j| {
            // SAFETY: `k` is below `chunk_size`, the length of every chunk.
            let sym = unsafe { *chunks[j].get_unchecked(k) };
            let step = &rows[usize::from(sym)][usize::from(next[j])];
            next[j] = sym;
            step
        }))
    };

    // SAFETY: AVX2 is available, as this function requires.
    unsafe { put_rounds(states, buf, next_round) };
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn load_states(states: &[u32; STATE_COUNT]) -> [__m256i; VECTORS] {
    std::array::from_fn(|v| {
        // SAFETY: `states` has `LANES` values from `v * LANES`.
        unsafe { _mm256_loadu_si256(states.as_ptr().add(v * LANES).cast()) }
    })
}

#[inline]
#[target_feature(enable = "avx2")]
unsafe fn store_states(s: &[__m256i; VECTORS], states: &mut [u32; STATE_COUNT]) {
    for (v, &s) in s.iter().enumerate() {
        // SAFETY: `states` has `LANES` values from `v * LANES`.
        unsafe { _mm256_storeu_si256(states.as_mut_ptr().add(v * LANES).cast(), s) };
    }
}

/// The high 32 bits of each lane's unsigned 64-bit product.
#[inline]
#[target_feature(enable = "avx2")]
fn mulhi_epu32(a: __m256i, b: __m256i) -> __m256i {
    // A 64-bit product of the even lanes has its high half in an odd lane, and one of the odd
    // lanes, shifted down to be multiplied, has it in the odd lane it came from.
    let even = _mm256_mul_epu32(a, b);
    let odd = _mm256_mul_epu32(_mm256_srli_epi64::<32>(a), _mm256_srli_epi64::<32>(b));
    _mm256_blend_epi32::<0b1010_1010>(_mm256_srli_epi64::<32>(even), odd)
}

/// Encodes rounds until `next_round` gives no more: `steps[j]` on state `j`, as
/// [`EncodingSymbol::put`] does for states 31 down to 0.
#[inline]
#[target_feature(enable = "avx2")]
unsafe fn put_rounds<'a, F>(
    states: &mut [u32; STATE_COUNT],
    buf: &mut ReverseWriter,
    mut next_round: F,
) where
    F: FnMut() -> Option<[&'a EncodingSymbol; STATE_COUNT]>,
{
    let one = _mm256_set1_epi32(1);
    let low_word = _mm256_set1_epi32(0xffff);

    // SAFETY: AVX2 is available, as this function requires.
    let mut s = unsafe { load_states(states) };

    while let Some(steps) = next_round() {
        let mut words = [[0u8; 16]; VECTORS];
        let mut counts = [0; VECTORS];

        for (v, s) in s.iter_mut().enumerate() {
            // SAFETY: a step is 16 bytes (see the assertion above).
            let step =
                |k: usize| unsafe { _mm_loadu_si128(ptr::from_ref(steps[v * LANES + k]).cast()) };

            // Two steps to a vector, those of lanes `k` and `k + 4`, so that unpacking within
            // each half leaves the fields in lane order.
            let t0 = _mm256_set_m128i(step(4), step(0));
            let t1 = _mm256_set_m128i(step(5), step(1));
            let t2 = _mm256_set_m128i(step(6), step(2));
            let t3 = _mm256_set_m128i(step(7), step(3));

            let x_max_rcp_01 = _mm256_unpacklo_epi32(t0, t1);
            let bias_last_01 = _mm256_unpackhi_epi32(t0, t1);
            let x_max_rcp_23 = _mm256_unpacklo_epi32(t2, t3);
            let bias_last_23 = _mm256_unpackhi_epi32(t2, t3);

            let x_max = _mm256_unpacklo_epi64(x_max_rcp_01, x_max_rcp_23);
            let rcp_freq = _mm256_unpackhi_epi64(x_max_rcp_01, x_max_rcp_23);
            let bias = _mm256_unpacklo_epi64(bias_last_01, bias_last_23);
            let last = _mm256_unpackhi_epi64(bias_last_01, bias_last_23);

            let cmpl_freq = _mm256_and_si256(last, low_word);
            let rcp_shift = _mm256_srli_epi32::<16>(last);

            // `x >= x_max`: a state is below 2^31, so the signed comparison is the unsigned
            // one (an `x_max` of 2^31, of a symbol with every slot, wraps to no state reaching
            // it).
            let x = *s;
            let renormalizes = _mm256_cmpgt_epi32(x, _mm256_sub_epi32(x_max, one));
            let mask = _mm256_movemask_ps(_mm256_castsi256_ps(renormalizes)) as usize;

            // The low words of the lanes that renormalize, in lane order, as the last 16-bit
            // words of `words[v]`.
            // SAFETY: `COMPRESS[mask]` has `LANES` values.
            let compress = unsafe { _mm256_loadu_si256(COMPRESS[mask].as_ptr().cast()) };
            let out = _mm256_permutevar8x32_epi32(_mm256_and_si256(x, low_word), compress);
            let out = _mm256_permute4x64_epi64::<0b1000>(_mm256_packus_epi32(out, out));

            // SAFETY: `words[v]` is 16 bytes.
            unsafe {
                _mm_storeu_si128(words[v].as_mut_ptr().cast(), _mm256_castsi256_si128(out));
            }
            counts[v] = mask.count_ones() as usize;

            let x = _mm256_blendv_epi8(x, _mm256_srli_epi32::<16>(x), renormalizes);
            let q = _mm256_srlv_epi32(mulhi_epu32(x, rcp_freq), rcp_shift);
            *s = _mm256_add_epi32(_mm256_add_epi32(x, bias), _mm256_mullo_epi32(q, cmpl_freq));
        }

        // State 31's word is written first, so it comes last: so do the last vector's words.
        for (words, count) in words.iter().zip(counts).rev() {
            buf.push_last_words(words, count);
        }
    }

    // SAFETY: AVX2 is available, as this function requires.
    unsafe { store_states(&s, states) };
}
