//! AVX2 decoding of 32-state order-0 rANS Nx16 streams.
//!
//! Each round decodes one symbol per state, 8 states per vector: the table entries are
//! gathered, the states advance together, the 32 symbols are written as they lie in the output,
//! and renormalization expands the next input words into the lanes that need them, in state
//! order. This is the same computation as the scalar loop, so the two can hand off between
//! rounds.

use std::{
    arch::x86_64::{
        _mm256_add_epi32, _mm256_and_si256, _mm256_i32gather_epi32, _mm256_mullo_epi32,
        _mm256_packus_epi16, _mm256_packus_epi32, _mm256_permutevar8x32_epi32, _mm256_set1_epi32,
        _mm256_setr_epi32, _mm256_setzero_si256, _mm256_srli_epi32, _mm256_storeu_si256,
    },
    io,
};

use super::{
    super::{
        SymbolTable,
        avx2::{STATE_COUNT, VECTORS, load_states, renormalize, store_states},
    },
    NORMALIZATION_BITS,
};

/// The most a frequency or an offset can be and still pack into an [`Entries`] value.
const MAX_FREQUENCY: u32 = 1 << NORMALIZATION_BITS;

/// The slots of a table.
const TABLE_SIZE: usize = 1 << NORMALIZATION_BITS;

/// A [`SymbolTable`] with a slot's symbol and step in one value, for a single gather:
/// `offset << 20 | (frequency - 1) << 8 | symbol`.
pub(super) struct Entries([u32; TABLE_SIZE]);

impl Entries {
    /// Packs `table`, or returns `None` if a step of it does not fit: a frequency of 0 or more
    /// than 4096, or an offset of 4096 or more, which no table of a valid stream has.
    pub(super) fn new(table: &SymbolTable) -> Option<Self> {
        let mut entries = [0; TABLE_SIZE];
        let mut fits = true;

        // Every slot is packed, whether it fits or not, so that the loop has no way out and the
        // compiler can vectorize it.
        for ((entry, &sym), &step) in entries.iter_mut().zip(&table.symbols).zip(&table.steps) {
            let (frequency, offset) = (step >> 16, step & 0xffff);
            let frequency = frequency.wrapping_sub(1);

            fits &= frequency < MAX_FREQUENCY && offset < MAX_FREQUENCY;
            *entry = (offset << 20) | ((frequency & 0xfff) << 8) | u32::from(sym);
        }

        (table.steps.len() == TABLE_SIZE && fits).then_some(Self(entries))
    }
}

/// Decodes a round into each of `dst` (state `j` writes byte `j` of its round), updating `states`
/// and `src` as the scalar loop does.
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
    entries: &Entries,
    states: &mut [u32; STATE_COUNT],
    src: &mut &[u8],
    dst: &mut [[u8; STATE_COUNT]],
) -> io::Result<()> {
    let slot_mask = _mm256_set1_epi32((1 << NORMALIZATION_BITS) - 1);
    let sym_mask = _mm256_set1_epi32(0xff);
    let one = _mm256_set1_epi32(1);

    // The two packs below work within each half of a vector, which leaves the symbols of the
    // four vectors' low halves before those of their high halves.
    let in_state_order = _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7);

    let table = entries.0.as_ptr().cast::<i32>();

    // SAFETY: AVX2 is available, as this function requires.
    let mut s = unsafe { load_states(states) };

    for round in dst {
        let mut syms = [_mm256_setzero_si256(); VECTORS];

        for (s, syms) in s.iter_mut().zip(&mut syms) {
            // SAFETY: a slot is below `1 << NORMALIZATION_BITS`, the length of the table.
            let entry =
                unsafe { _mm256_i32gather_epi32::<4>(table, _mm256_and_si256(*s, slot_mask)) };

            let frequency = _mm256_add_epi32(
                _mm256_and_si256(_mm256_srli_epi32::<8>(entry), slot_mask),
                one,
            );
            let offset = _mm256_srli_epi32::<20>(entry);
            let x = _mm256_srli_epi32::<{ NORMALIZATION_BITS as i32 }>(*s);
            *s = _mm256_add_epi32(_mm256_mullo_epi32(frequency, x), offset);

            *syms = _mm256_and_si256(entry, sym_mask);
        }

        let bytes = _mm256_packus_epi16(
            _mm256_packus_epi32(syms[0], syms[1]),
            _mm256_packus_epi32(syms[2], syms[3]),
        );

        // SAFETY: `round` is `STATE_COUNT` bytes, a vector's worth.
        unsafe {
            _mm256_storeu_si256(
                round.as_mut_ptr().cast(),
                _mm256_permutevar8x32_epi32(bytes, in_state_order),
            );
        }

        // SAFETY: AVX2 is available, as this function requires.
        unsafe { renormalize(&mut s, src)? };
    }

    // SAFETY: AVX2 is available, as this function requires.
    unsafe { store_states(&s, states) };

    Ok(())
}
