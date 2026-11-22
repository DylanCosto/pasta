//! AVX2 decoding of 32-state order-1 rANS Nx16 streams.
//!
//! Each round decodes one symbol per state, 8 states per vector: the context offsets and table
//! entries are gathered, the states advance together, and renormalization expands the next input
//! words into the lanes that need them, in state order. This is the same computation as the
//! scalar loop, so the two can hand off between rounds.

use std::{
    arch::x86_64::{
        __m256i, _mm_cvtsi32_si128, _mm256_add_epi32, _mm256_and_si256, _mm256_i32gather_epi32,
        _mm256_loadu_si256, _mm256_mullo_epi32, _mm256_set1_epi32, _mm256_srl_epi32,
        _mm256_srli_epi32, _mm256_storeu_si256,
    },
    io,
};

use super::{
    super::avx2::{LANES, STATE_COUNT, VECTORS, load_states, renormalize, store_states},
    Tables,
};

/// Decodes rounds `0..chunk_size` of the 32 interleaved streams into `dst` (state `j` writes
/// `dst[j * chunk_size + i]` in round `i`), updating `states`, `prev_syms`, and `src` as the scalar
/// loop does.
///
/// # Errors
///
/// Returns an [`io::ErrorKind::UnexpectedEof`] error if `src` ends before the last round, as the
/// scalar loop does.
///
/// # Safety
///
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(super) unsafe fn decode(
    tables: &Tables,
    bits: u32,
    states: &mut [u32; STATE_COUNT],
    prev_syms: &mut [u8; STATE_COUNT],
    src: &mut &[u8],
    dst: &mut [u8],
    chunk_size: usize,
) -> io::Result<()> {
    assert_eq!(dst.len(), STATE_COUNT * chunk_size);

    let slot_mask = _mm256_set1_epi32((1 << bits) - 1);
    let shift = _mm_cvtsi32_si128(bits as i32);
    let sym_mask = _mm256_set1_epi32(0xff);
    let frequency_mask = _mm256_set1_epi32(0xfff);
    let one = _mm256_set1_epi32(1);

    let offsets = tables.offsets.as_ptr().cast::<i32>();
    let entries = tables.entries.as_ptr().cast::<i32>();

    // SAFETY: AVX2 is available, as this function requires.
    let mut s = unsafe { load_states(states) };

    let mut prev: [__m256i; VECTORS] = std::array::from_fn(|v| {
        let syms: [u32; LANES] = std::array::from_fn(|k| u32::from(prev_syms[v * LANES + k]));
        // SAFETY: `syms` has `LANES` values.
        unsafe { _mm256_loadu_si256(syms.as_ptr().cast()) }
    });

    let mut syms = [0u32; STATE_COUNT];

    for round in 0..chunk_size {
        for v in 0..VECTORS {
            // SAFETY: `prev` lanes are symbols (< 256 = `offsets.len()`), and every offset plus a
            // slot (< `1 << bits`) is within `entries` (see `Tables`). Both tables are far smaller
            // than `i32::MAX` elements.
            let entry = unsafe {
                let offset = _mm256_i32gather_epi32::<4>(offsets, prev[v]);
                let i = _mm256_add_epi32(offset, _mm256_and_si256(s[v], slot_mask));
                _mm256_i32gather_epi32::<4>(entries, i)
            };

            let frequency = _mm256_add_epi32(
                _mm256_and_si256(_mm256_srli_epi32::<8>(entry), frequency_mask),
                one,
            );
            let offset = _mm256_srli_epi32::<20>(entry);
            let x = _mm256_srl_epi32(s[v], shift);
            s[v] = _mm256_add_epi32(_mm256_mullo_epi32(frequency, x), offset);

            prev[v] = _mm256_and_si256(entry, sym_mask);

            // SAFETY: `syms` has `LANES` values from `v * LANES`.
            unsafe { _mm256_storeu_si256(syms.as_mut_ptr().add(v * LANES).cast(), prev[v]) };
        }

        for (j, &sym) in syms.iter().enumerate() {
            dst[j * chunk_size + round] = sym as u8;
        }

        // SAFETY: AVX2 is available, as this function requires.
        unsafe { renormalize(&mut s, src)? };
    }

    // SAFETY: AVX2 is available, as this function requires.
    unsafe { store_states(&s, states) };

    for (v, &prev) in prev.iter().enumerate() {
        // SAFETY: `syms` has `LANES` values from `v * LANES`.
        unsafe { _mm256_storeu_si256(syms.as_mut_ptr().add(v * LANES).cast(), prev) };
    }

    for (prev_sym, &sym) in prev_syms.iter_mut().zip(&syms) {
        *prev_sym = sym as u8;
    }

    Ok(())
}
