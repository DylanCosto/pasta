mod header;
mod order_0;
mod order_1;

use std::{
    cell::Cell,
    io::{self, Write},
};

use crate::cram::io::writer::num::write_u32_le;

use self::header::write_header;
use super::{LOWER_BOUND, Order, STATE_COUNT};

pub fn encode(order: Order, src: &[u8]) -> io::Result<Vec<u8>> {
    match order {
        Order::Zero => order_0::encode(src),
        Order::One => order_1::encode(src),
    }
}

fn write_states<W>(writer: &mut W, states: &[u32; STATE_COUNT]) -> io::Result<()>
where
    W: Write,
{
    for state in states {
        write_u32_le(writer, *state)?;
    }

    Ok(())
}

/// The bits of the frequency totals.
const BITS: u32 = 12;

/// A symbol's encoding step with a precomputed reciprocal, so encoding multiplies instead of
/// dividing (ryg_rans' byte-wise `RansEncSymbol`, which htscodecs also uses). It gives the same
/// states and bytes as [`state_renormalize`] and [`state_step`].
#[derive(Clone, Copy, Debug, Default)]
struct EncodingSymbol {
    x_max: u32,
    rcp_freq: u32,
    bias: u32,
    // 16 bits each, so that a symbol takes 16 bytes.
    cmpl_freq: u16,
    rcp_shift: u16,
}

impl EncodingSymbol {
    /// The step for a symbol with frequency `f` (nonzero) and cumulative frequency `g`, with
    /// frequencies summing to `1 << BITS`.
    fn new(f: u16, g: u16) -> Self {
        let (f, g) = (u32::from(f), u32::from(g));
        let x_max = (LOWER_BOUND >> (BITS - 8)) * f;
        // At most `1 << BITS` - 1, as `f` is at least 1.
        let cmpl_freq = ((1 << BITS) - f) as u16;

        if f < 2 {
            Self {
                x_max,
                rcp_freq: u32::MAX,
                rcp_shift: 0,
                bias: g + (1 << BITS) - 1,
                cmpl_freq,
            }
        } else {
            let shift = 32 - (f - 1).leading_zeros();
            let rcp_freq = (1u64 << (shift + 31)).div_ceil(u64::from(f));

            Self {
                x_max,
                rcp_freq: rcp_freq as u32,
                // At most `BITS`.
                rcp_shift: (shift - 1) as u16,
                bias: g,
                cmpl_freq,
            }
        }
    }

    /// Renormalizes `x` (writing bytes to `dst`) and encodes the symbol.
    #[inline]
    fn put(&self, mut x: u32, dst: &mut ReverseWriter) -> u32 {
        while x >= self.x_max {
            dst.push(x as u8);
            x >>= 8;
        }

        let q = ((u64::from(x) * u64::from(self.rcp_freq)) >> 32) as u32 >> self.rcp_shift;
        x + self.bias + q * u32::from(self.cmpl_freq)
    }
}

thread_local! {
    /// The buffer of the last [`ReverseWriter`] on this thread, so that the next does not
    /// allocate and zero one for each block.
    static SCRATCH: Cell<Vec<u8>> = const { Cell::new(Vec::new()) };
}

/// The encoded bytes, written from the end of a buffer toward its start: encoders emit them last
/// first, so they end up in order without being reversed.
struct ReverseWriter {
    buf: Vec<u8>,
    end: usize,
    pos: usize,
}

impl ReverseWriter {
    /// A writer for the bytes of encoding `len` symbols: at most 2 each, as a state below
    /// `LOWER_BOUND << 8` (31 bits) is at most 2 bytes over a symbol's lowest bound (19 bits).
    fn for_symbols(len: usize) -> Self {
        let end = 2 * len + 16;
        let mut buf = SCRATCH.take();

        if buf.len() < end {
            buf.resize(end, 0);
        }

        Self { buf, end, pos: end }
    }

    #[inline]
    fn push(&mut self, b: u8) {
        self.pos -= 1;
        self.buf[self.pos] = b;
    }

    fn as_slice(&self) -> &[u8] {
        &self.buf[self.pos..self.end]
    }
}

impl Drop for ReverseWriter {
    fn drop(&mut self) {
        SCRATCH.set(std::mem::take(&mut self.buf));
    }
}

#[cfg(test)]
fn state_renormalize(mut s: u32, f: u16, dst: &mut Vec<u8>) -> u32 {
    while s >= (LOWER_BOUND >> 4) * u32::from(f) {
        dst.push((s & 0xff) as u8);
        s >>= 8;
    }

    s
}

#[cfg(test)]
fn state_step(s: u32, f: u16, g: u16) -> u32 {
    let (q, r) = (s / u32::from(f), s % u32::from(f));
    (q << 12) + r + u32::from(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encoding_symbol_matches_division() {
        let total = 1u16 << BITS;

        for f in [1, 2, 3, 7, 100, 1000, total / 2, total - 1, total] {
            let g = (total - f) / 2;
            let symbol = EncodingSymbol::new(f, g);

            for x in [
                LOWER_BOUND,
                LOWER_BOUND + 1,
                0x12_3456,
                0x7fff_ffff,
                u32::MAX,
            ] {
                let (mut a, mut b) = (Vec::new(), ReverseWriter::for_symbols(1));
                let expected = state_step(state_renormalize(x, f, &mut a), f, g);
                assert_eq!(symbol.put(x, &mut b), expected, "f {f}, x {x:#x}");
                a.reverse();
                assert_eq!(a, b.as_slice(), "f {f}, x {x:#x}");
            }
        }
    }

    #[test]
    fn test_encode_with_order_0() -> io::Result<()> {
        let data = b"noodles";
        let actual = encode(Order::Zero, data)?;

        let expected = [
            0x00, 0x25, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x64, 0x82, 0x49, 0x65, 0x00,
            0x82, 0x49, 0x6c, 0x82, 0x49, 0x6e, 0x82, 0x49, 0x6f, 0x00, 0x84, 0x92, 0x73, 0x82,
            0x49, 0x00, 0xe2, 0x06, 0x83, 0x18, 0x74, 0x7b, 0x41, 0x0c, 0x2b, 0xa9, 0x41, 0x0c,
            0x25, 0x31, 0x80, 0x03,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_with_order_1() -> io::Result<()> {
        let data = b"noodles";
        let actual = encode(Order::One, data)?;

        // htscodecs' (`rans_compress`), with a row for the context of the last symbol, `s`.
        let expected = [
            0x01, 0x40, 0x00, 0x00, 0x00, 0x07, 0x00, 0x00, 0x00, 0x00, 0x64, 0x84, 0x00, 0x6e,
            0x84, 0x00, 0x6f, 0x00, 0x87, 0xff, 0x00, 0x64, 0x6c, 0x8f, 0xff, 0x00, 0x65, 0x00,
            0x73, 0x8f, 0xff, 0x00, 0x6c, 0x65, 0x8f, 0xff, 0x00, 0x6e, 0x6f, 0x8f, 0xff, 0x00,
            0x6f, 0x00, 0x64, 0x87, 0xff, 0x6f, 0x88, 0x00, 0x00, 0x73, 0x00, 0x8f, 0xff, 0x00,
            0x00, 0x00, 0x04, 0x00, 0x02, 0x02, 0x28, 0x00, 0x01, 0x02, 0x28, 0x00, 0x01, 0x02,
            0x60, 0x00, 0x02,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
