//! What the AVX2 decoders of 32-state rANS Nx16 streams share: the layout of the states in
//! vectors, and their renormalization.

use std::{
    arch::x86_64::{
        __m256i, _mm_loadu_si128, _mm256_blendv_epi8, _mm256_castsi256_ps, _mm256_cmpeq_epi32,
        _mm256_cvtepu16_epi32, _mm256_loadu_si256, _mm256_min_epu32, _mm256_movemask_ps,
        _mm256_or_si256, _mm256_permutevar8x32_epi32, _mm256_set1_epi32, _mm256_setzero_si256,
        _mm256_slli_epi32, _mm256_storeu_si256,
    },
    io,
};

pub(super) const STATE_COUNT: usize = 32;
pub(super) const LANES: usize = 8;
pub(super) const VECTORS: usize = STATE_COUNT / LANES;

/// The input a round may read: two bytes per state.
const ROUND_WINDOW: usize = 2 * STATE_COUNT;

/// For each 8-bit mask of lanes that renormalize, the input word each lane takes: lane `k` takes
/// word `popcount(mask & ((1 << k) - 1))`. Lanes that do not renormalize take word 0, unused.
static EXPAND: [[u32; LANES]; 256] = build_expand();

const fn build_expand() -> [[u32; LANES]; 256] {
    let mut table = [[0; LANES]; 256];
    let mut mask = 0;

    while mask < 256 {
        let mut word = 0;
        let mut lane = 0;

        while lane < LANES {
            if mask & (1 << lane) != 0 {
                table[mask][lane] = word;
                word += 1;
            }

            lane += 1;
        }

        mask += 1;
    }

    table
}

/// Loads the states, 8 per vector in order.
///
/// # Safety
///
/// The CPU must support AVX2.
#[inline]
#[target_feature(enable = "avx2")]
pub(super) unsafe fn load_states(states: &[u32; STATE_COUNT]) -> [__m256i; VECTORS] {
    std::array::from_fn(|v| {
        // SAFETY: `states` has `LANES` values from `v * LANES`.
        unsafe { _mm256_loadu_si256(states.as_ptr().add(v * LANES).cast()) }
    })
}

/// Stores the states [`load_states`] loaded.
///
/// # Safety
///
/// The CPU must support AVX2.
#[inline]
#[target_feature(enable = "avx2")]
pub(super) unsafe fn store_states(s: &[__m256i; VECTORS], states: &mut [u32; STATE_COUNT]) {
    for (v, &s) in s.iter().enumerate() {
        // SAFETY: `states` has `LANES` values from `v * LANES`.
        unsafe { _mm256_storeu_si256(states.as_mut_ptr().add(v * LANES).cast(), s) };
    }
}

/// Renormalizes the states in order, each taking the next input word if it is below 2^15, and
/// advances `src` past the words taken, as the scalar `renormalize_states` does.
///
/// # Errors
///
/// Returns an [`io::ErrorKind::UnexpectedEof`] error if `src` has fewer words than the states
/// take.
///
/// # Safety
///
/// The CPU must support AVX2.
#[inline]
#[target_feature(enable = "avx2")]
pub(super) unsafe fn renormalize(s: &mut [__m256i; VECTORS], src: &mut &[u8]) -> io::Result<()> {
    let renormalize_max = _mm256_set1_epi32((1 << 15) - 1);

    let mut renormalize = [_mm256_setzero_si256(); VECTORS];
    let mut masks = [0; VECTORS];
    let mut needed = 0;

    for ((state, renormalize), mask) in s.iter().zip(&mut renormalize).zip(&mut masks) {
        // `min(s, 2^15 - 1) == s` is an unsigned `s < 2^15`.
        *renormalize = _mm256_cmpeq_epi32(_mm256_min_epu32(*state, renormalize_max), *state);
        *mask = _mm256_movemask_ps(_mm256_castsi256_ps(*renormalize)) as usize;
        needed += 2 * mask.count_ones() as usize;
    }

    // Near the end of the input, the loads read a padded copy. The padding only reaches lanes
    // that do not renormalize.
    let tail;

    let window: &[u8; ROUND_WINDOW] = match src.first_chunk() {
        Some(window) => window,
        None if needed <= src.len() => {
            let mut padded = [0; ROUND_WINDOW];
            padded[..src.len()].copy_from_slice(src);
            tail = padded;
            &tail
        }
        None => return Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
    };

    let mut n = 0;

    for ((state, renormalize), mask) in s.iter_mut().zip(renormalize).zip(masks) {
        // SAFETY: `n <= 2 * LANES * v` before vector `v`, so the 16-byte load ends within
        // `window` (`2 * STATE_COUNT` bytes), and `EXPAND[mask]` has `LANES` values.
        let words = unsafe {
            let words = _mm256_cvtepu16_epi32(_mm_loadu_si128(window.as_ptr().add(n).cast()));
            let expand = _mm256_loadu_si256(EXPAND[mask].as_ptr().cast());
            _mm256_permutevar8x32_epi32(words, expand)
        };

        let renormalized = _mm256_or_si256(_mm256_slli_epi32::<16>(*state), words);
        *state = _mm256_blendv_epi8(*state, renormalized, renormalize);

        n += 2 * mask.count_ones() as usize;
    }

    debug_assert_eq!(n, needed);
    *src = &src[n..];

    Ok(())
}
