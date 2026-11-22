use std::io;

use super::{
    NORMALIZATION_BITS, TABLE_SIZE, build_table, decode_symbol, read_states, renormalize_states,
    state_renormalize,
};
use crate::cram::{
    codecs::rans_4x8::{ALPHABET_SIZE, STATE_COUNT},
    io::reader::num::{read_itf8_as, read_u8},
};

type Frequencies = [u16; ALPHABET_SIZE]; // F

/// The shortest output decoded with a [`SlotTable`]. A shorter one does not repay the time the
/// larger table takes to fill: measured, 128 bytes decode faster without it, 2 kB as fast or
/// faster with it, and the lengths between about the same either way.
const SLOT_TABLE_MIN_LEN: usize = 2048;

pub fn decode(src: &mut &[u8], dst: &mut [u8]) -> io::Result<()> {
    let frequencies = read_frequencies(src)?;
    validate_frequencies(&frequencies)?;

    if dst.len() < SLOT_TABLE_MIN_LEN {
        let mut symbols = [0; TABLE_SIZE];
        let mut packed_frequencies = [0; ALPHABET_SIZE];
        build_table(&frequencies, &mut symbols, &mut packed_frequencies);

        decode_with(src, dst, |state| {
            decode_symbol(&symbols, &packed_frequencies, state)
        })
    } else {
        let table = SlotTable::new(&frequencies);
        decode_with(src, dst, |state| table.decode_symbol(state))
    }
}

/// Decodes `dst` with the four interleaved states, each symbol by `decode_symbol`.
#[inline(always)]
fn decode_with<F>(src: &mut &[u8], dst: &mut [u8], decode_symbol: F) -> io::Result<()>
where
    F: Fn(&mut u32) -> u8,
{
    let mut states = read_states(src)?;
    let (chunks, remainder) = dst.as_chunks_mut::<STATE_COUNT>();

    for chunk in chunks {
        for (d, state) in chunk.iter_mut().zip(&mut states) {
            *d = decode_symbol(state);
        }

        // The state steps are independent, so they run ahead of the ordered renormalization reads.
        renormalize_states(&mut states, src)?;
    }

    for (d, state) in remainder.iter_mut().zip(&mut states) {
        *d = decode_symbol(state);
        *state = state_renormalize(*state, src)?;
    }

    Ok(())
}

/// A decoding table with a step for each cumulative frequency slot.
///
/// `symbols` maps a slot to its symbol, as [`build_table`] fills it, and `steps[slot]` packs
/// `frequency << 16 | (slot - cumulative_frequency)` for that symbol. [`build_table`]'s
/// frequencies are indexed by symbol instead, so a step waits for the symbol to load; here the
/// two loads are independent. One table of slots is 16 kB, which order 1, with a table for each
/// context, cannot afford.
struct SlotTable {
    symbols: [u8; TABLE_SIZE],
    steps: [u32; TABLE_SIZE],
}

impl SlotTable {
    /// Builds the table of `frequencies`, whose total is at most 4096 (see
    /// [`validate_frequencies`]).
    fn new(frequencies: &Frequencies) -> Self {
        let mut symbols = [0; TABLE_SIZE];
        let mut steps = [0; TABLE_SIZE];
        let mut cumulative_frequency = 0;

        for (sym, &f) in (0..=u8::MAX).zip(frequencies) {
            let next_cumulative_frequency = cumulative_frequency + usize::from(f);
            let range = cumulative_frequency..next_cumulative_frequency;

            symbols[range.clone()].fill(sym);

            for (offset, step) in (0..).zip(&mut steps[range]) {
                *step = (u32::from(f) << 16) | offset;
            }

            cumulative_frequency = next_cumulative_frequency;
        }

        // Slots at or past the total frequency resolve to the last symbol, as in
        // [`build_table`], and count on from its own slots.
        let last_frequency = u32::from(frequencies[usize::from(u8::MAX)]);
        symbols[cumulative_frequency..].fill(u8::MAX);

        for (offset, step) in (last_frequency..).zip(&mut steps[cumulative_frequency..]) {
            *step = (last_frequency << 16) | offset;
        }

        Self { symbols, steps }
    }

    /// Decodes a symbol and advances the state, before renormalization, as [`decode_symbol`]
    /// does.
    #[inline(always)]
    fn decode_symbol(&self, state: &mut u32) -> u8 {
        let slot = (*state & 0x0fff) as usize;
        let step = self.steps[slot];

        // Valid tables cannot overflow here; wrapping keeps malformed initial states panic-free.
        *state = (step >> 16)
            .wrapping_mul(*state >> NORMALIZATION_BITS)
            .wrapping_add(step & 0xffff);

        self.symbols[slot]
    }
}

pub(super) fn read_frequencies(src: &mut &[u8]) -> io::Result<Frequencies> {
    const NUL: u8 = 0x00;

    let mut frequencies = [0; ALPHABET_SIZE];

    let mut sym = read_u8(src)?;
    let mut prev_sym = sym;

    loop {
        let f = read_itf8_as(src)?;
        frequencies[usize::from(sym)] = f;

        sym = read_u8(src)?;

        if sym == NUL {
            break;
        }

        if sym - 1 == prev_sym {
            let len = read_u8(src)?;

            for _ in 0..len {
                let f = read_itf8_as(src)?;
                frequencies[usize::from(sym)] = f;
                // A malformed run can pass 255; wrap (as release builds did) instead of panicking.
                sym = sym.wrapping_add(1);
            }
        }

        prev_sym = sym;
    }

    Ok(frequencies)
}

/// Rejects tables whose total exceeds 4096, which cannot come from a valid encoder and would
/// overflow the cumulative frequencies.
fn validate_frequencies(frequencies: &Frequencies) -> io::Result<()> {
    const TOTAL_FREQUENCY: u32 = 1 << NORMALIZATION_BITS;

    let total: u32 = frequencies.iter().copied().map(u32::from).sum();

    if total > TOTAL_FREQUENCY {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid rANS order-0 frequency table",
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_rejects_oversized_frequency_tables() {
        // 'a' = 4095 and 'c' = f, with no states following: a total of 4096 is valid (and reading
        // the states fails), but 4097 is not.
        for (f, expected) in [
            (0x01, io::ErrorKind::UnexpectedEof),
            (0x02, io::ErrorKind::InvalidData),
        ] {
            let src = [b'a', 0x8f, 0xff, b'c', f, 0x00];
            let mut dst = [0; 8];
            let actual = decode(&mut &src[..], &mut dst).unwrap_err().kind();
            assert_eq!(actual, expected, "f = {f}");
        }
    }

    /// Outputs on either side of [`SLOT_TABLE_MIN_LEN`] are decoded with different tables, to
    /// the same bytes.
    #[test]
    fn test_decode_round_trips_with_either_table() -> io::Result<()> {
        use crate::cram::codecs::rans_4x8::{self, Order};

        let mut x = 0x9e37_79b9_u32;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            x
        };

        // Quality-like bytes, and bytes of every value.
        let qualities: Vec<u8> = (0..70_000).map(|_| 33 + (next() % 41) as u8).collect();
        let uniform: Vec<u8> = (0..70_000).map(|_| next() as u8).collect();

        for data in [&qualities[..], &uniform[..]] {
            for len in [
                1,
                3,
                4,
                5,
                SLOT_TABLE_MIN_LEN - 1,
                SLOT_TABLE_MIN_LEN,
                SLOT_TABLE_MIN_LEN + 1,
                4_099,
                data.len(),
            ] {
                let data = &data[..len];
                let encoded = rans_4x8::encode(Order::Zero, data)?;
                assert_eq!(rans_4x8::decode(&encoded)?, data, "len = {len}");
            }
        }

        Ok(())
    }

    /// The two tables give the same symbols and states, also for the slots past the total of
    /// frequencies that add up to less than 4096.
    #[test]
    fn test_slot_table_matches_build_table() {
        let mut frequencies = [0; ALPHABET_SIZE];
        frequencies[usize::from(b'a')] = 1000;
        frequencies[usize::from(b'c')] = 3;
        frequencies[usize::from(u8::MAX)] = 77;

        let mut symbols = [0; TABLE_SIZE];
        let mut packed_frequencies = [0; ALPHABET_SIZE];
        build_table(&frequencies, &mut symbols, &mut packed_frequencies);

        let table = SlotTable::new(&frequencies);

        for slot in 0..TABLE_SIZE as u32 {
            for high in [0, 1, 0x7ff, 0xf_ffff] {
                let state = (high << NORMALIZATION_BITS) | slot;

                let mut expected_state = state;
                let expected = decode_symbol(&symbols, &packed_frequencies, &mut expected_state);

                let mut actual_state = state;
                let actual = table.decode_symbol(&mut actual_state);

                assert_eq!(
                    (actual, actual_state),
                    (expected, expected_state),
                    "slot {slot}"
                );
            }
        }
    }

    #[test]
    fn test_read_frequencies() -> io::Result<()> {
        let src = [
            b'a', // symbol = 'a'
            0x05, // frequencies['a'] = 5
            b'b', // symbol = 'b'
            0x02, // run length = 2
            0x02, // frequencies['b'] = 2
            0x01, // frequencies['c'] = 1
            0x01, // frequencies['d'] = 1
            b'r', // symbol = 'r'
            0x02, // frequencies['r'] = 2
            0x00, // EOF
        ];

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(read_frequencies(&mut &src[..])?, expected);

        Ok(())
    }
}
