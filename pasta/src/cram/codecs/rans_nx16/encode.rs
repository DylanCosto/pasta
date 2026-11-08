pub mod bit_pack;
mod order_0;
mod order_1;
mod rle;
mod stripe;

use std::{
    borrow::Cow,
    cell::Cell,
    io::{self, Write},
};

use super::{ALPHABET_SIZE, Flags};
use crate::cram::io::writer::num::{write_u8, write_u32_le, write_uint7};

// § 3 "rANS Nx16" (2023-03-15): "The lower-bound and initial encoder state _L_ is [...] 0x8000."
const LOWER_BOUND: u32 = 0x8000;

/// The fewest bytes coded with order 1.
const MIN_ORDER_1_LEN: usize = 8;

pub fn encode(mut flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    // As htscodecs does, inputs of 20 bytes or fewer are not striped.
    if src.len() <= 20 {
        flags.remove(Flags::STRIPE);
    }

    let mut src = Cow::from(src);
    let mut dst = Vec::new();

    write_flags(&mut dst, flags)?;

    if flags.has_uncompressed_size() {
        write_uncompressed_size(&mut dst, src.len())?;
    }

    let state_count = flags.state_count();

    if flags.is_striped() {
        let buf = stripe::encode(flags, &src)?;
        dst.extend(&buf);
        return Ok(dst);
    }

    if flags.is_bit_packed() {
        match bit_pack::build_context(&src) {
            Ok(ctx) => {
                src = Cow::from(bit_pack::encode(&src, &ctx));
                bit_pack::write_context(&mut dst, &ctx, src.len())?;
            }
            Err(
                bit_pack::context::BuildContextError::EmptyAlphabet
                | bit_pack::context::BuildContextError::TooManySymbols(_),
            ) => {
                flags.remove(Flags::PACK);
                dst[0] = u8::from(flags);
            }
        }
    }

    if flags.is_rle() {
        let literals = match rle::build_context(&src) {
            Ok(mut ctx) => rle::encode(&src, &mut ctx).map(|literals| Some((ctx, literals)))?,
            Err(rle::context::BuildContextError::EmptyAlphabet) => None,
        };

        match literals {
            // Like htscodecs, run-length encode only when it saves at least 1%.
            Some((ctx, literals))
                if ((literals.len() + ctx.dst.len()) as f64) < 0.99 * src.len() as f64 =>
            {
                rle::write_context(&mut dst, &ctx, literals.len(), state_count)?;
                src = Cow::from(literals);
            }
            _ => {
                flags.remove(Flags::RLE);
                dst[0] = u8::from(flags);
            }
        }
    }

    // As htscodecs does, fewer than 8 bytes are coded with order 0.
    if src.len() < MIN_ORDER_1_LEN {
        flags.remove(Flags::ORDER);
        dst[0] = u8::from(flags);
    }

    // Too few bytes for the states, let alone a frequency table, are stored as they are.
    if src.len() < state_count {
        flags.remove(Flags::ORDER);
        flags.insert(Flags::CAT);
        dst[0] = u8::from(flags);
    }

    if flags.is_uncompressed() {
        dst.write_all(&src)?;
        return Ok(dst);
    }

    let start = dst.len();

    if flags.order() == 0 {
        let ctx = order_0::build_context(&src)?;
        order_0::write_context(&mut dst, &ctx)?;
        order_0::encode(&src, &ctx, state_count, &mut dst)?;
    } else {
        let ctx = order_1::build_context(&src, state_count)?;
        order_1::write_context(&mut dst, &ctx)?;
        order_1::encode(&src, &ctx, state_count, &mut dst)?;
    }

    // As htscodecs does, data that codes to no fewer bytes is stored as it is.
    if dst.len() - start >= src.len() {
        flags.remove(Flags::ORDER);
        flags.insert(Flags::CAT);
        dst[0] = u8::from(flags);
        dst.truncate(start);
        dst.extend_from_slice(&src);
    }

    Ok(dst)
}

fn write_flags(dst: &mut Vec<u8>, flags: Flags) -> io::Result<()> {
    write_u8(dst, u8::from(flags))
}

fn write_uncompressed_size(dst: &mut Vec<u8>, size: usize) -> io::Result<()> {
    let n = u32::try_from(size).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(dst, n)
}

fn write_states(dst: &mut Vec<u8>, states: &[u32]) -> io::Result<()> {
    for &state in states {
        write_u32_le(dst, state)?;
    }

    Ok(())
}

fn write_alphabet(dst: &mut Vec<u8>, alphabet: &[bool; ALPHABET_SIZE]) -> io::Result<()> {
    const NUL: u8 = 0x00;

    let mut iter = alphabet.iter().enumerate();
    let mut prev_sym = None;

    while let Some((sym, &a)) = iter.next() {
        if !a {
            continue;
        }

        // SAFETY: `sym < ALPHABET_SIZE`.
        write_u8(dst, sym as u8)?;

        if sym > 0 && prev_sym == Some(sym - 1) {
            let i = sym + 1;

            let len = alphabet[i..]
                .iter()
                .position(|&b| !b)
                .unwrap_or(ALPHABET_SIZE - i);

            // SAFETY: `len < ALPHABET_SIZE`.
            write_u8(dst, len as u8)?;
            for _ in iter.by_ref().take(len) {}
        }

        prev_sym = Some(sym);
    }

    write_u8(dst, NUL)?;

    Ok(())
}

/// A symbol's encoding step with a precomputed reciprocal, so encoding multiplies instead of
/// dividing (ryg_rans' `RansEncSymbol`, which htscodecs also uses). It gives the same states as
/// [`state_renormalize`] and [`state_step`].
///
/// It is 16 bytes, as htscodecs' is, so that an order-1 table's steps do not straddle cache lines.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct EncodingSymbol {
    x_max: u32,
    rcp_freq: u32,
    bias: u32,
    /// `(1 << bits) - f`, at most 4095.
    cmpl_freq: u16,
    rcp_shift: u16,
}

impl EncodingSymbol {
    /// The step for a symbol with frequency `f` (nonzero) and cumulative frequency `g`, with
    /// frequencies summing to `1 << bits` (at most 2^12).
    pub(super) fn new(f: u32, g: u32, bits: u32) -> Self {
        let x_max = (1 << (31 - bits)) * f;
        let cmpl_freq = ((1 << bits) - f) as u16;

        if f < 2 {
            Self {
                x_max,
                rcp_freq: u32::MAX,
                bias: g + (1 << bits) - 1,
                cmpl_freq,
                rcp_shift: 0,
            }
        } else {
            let shift = 32 - (f - 1).leading_zeros();
            let rcp_freq = (1u64 << (shift + 31)).div_ceil(u64::from(f));

            Self {
                x_max,
                rcp_freq: rcp_freq as u32,
                bias: g,
                cmpl_freq,
                rcp_shift: (shift - 1) as u16,
            }
        }
    }

    /// Renormalizes `x` (writing to `dst`) and encodes the symbol.
    ///
    /// A state is below 2^31 and a symbol's bound at least 2^19, so renormalizing shifts out one
    /// word at most. Whether it does is hard to predict, so this does not branch on it.
    #[inline]
    pub(super) fn put(&self, x: u32, dst: &mut ReverseWriter) -> u32 {
        let renormalizes = x >= self.x_max;
        dst.push_if(renormalizes, x as u16);
        let x = if renormalizes { x >> 16 } else { x };

        let q = ((u64::from(x) * u64::from(self.rcp_freq)) >> 32) as u32 >> self.rcp_shift;
        x + self.bias + q * u32::from(self.cmpl_freq)
    }
}

thread_local! {
    /// The buffer of the last [`ReverseWriter`] on this thread, so that the next does not
    /// allocate and zero one for each block.
    static SCRATCH: Cell<Vec<u8>> = const { Cell::new(Vec::new()) };
}

/// The encoded 16-bit words, written from the end of a buffer toward its start, each little
/// endian: encoders emit them last first, so they end up in order without being reversed.
pub(super) struct ReverseWriter {
    buf: Vec<u8>,
    end: usize,
    pos: usize,
}

impl ReverseWriter {
    /// A writer for the words of encoding `len` symbols: at most one each, as a state below
    /// `LOWER_BOUND << 16` (31 bits) is at most 16 bits over a symbol's lowest bound (19 bits).
    pub(super) fn for_symbols(len: usize) -> Self {
        let end = 2 * len + 16;
        let mut buf = SCRATCH.take();

        if buf.len() < end {
            buf.resize(end, 0);
        }

        Self { buf, end, pos: end }
    }

    /// Writes `n` if `is_written`. The word goes before the last either way, which the slack of
    /// [`Self::for_symbols`] has room for, and counts only if written.
    #[inline]
    fn push_if(&mut self, is_written: bool, n: u16) {
        self.buf[self.pos - 2..self.pos].copy_from_slice(&n.to_le_bytes());
        self.pos -= 2 * usize::from(is_written);
    }

    pub(super) fn as_slice(&self) -> &[u8] {
        &self.buf[self.pos..self.end]
    }
}

impl Drop for ReverseWriter {
    fn drop(&mut self) {
        SCRATCH.set(std::mem::take(&mut self.buf));
    }
}

#[cfg(test)]
fn state_step(s: u32, f: u32, g: u32, bits: u32) -> u32 {
    let (q, r) = (s / f, s % f);
    (q << bits) + r + g
}

#[cfg(test)]
fn state_renormalize(mut s: u32, f: u32, bits: u32, dst: &mut Vec<u8>) -> u32 {
    while s >= (1 << (31 - bits)) * f {
        write_u16_be(dst, (s & 0xffff) as u16);
        s >>= 16;
    }

    s
}

#[cfg(test)]
fn write_u16_be(dst: &mut Vec<u8>, n: u16) {
    let buf = n.to_be_bytes();
    dst.extend(buf);
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_encoding_symbol_matches_division() {
        use super::{EncodingSymbol, state_renormalize, state_step};

        for bits in [10, 12] {
            let total = 1u32 << bits;

            for f in [1, 2, 3, 7, 100, 1000, total / 2, total - 1] {
                let g = (total - f) / 2;
                let symbol = EncodingSymbol::new(f, g, bits);

                for x in [0x8000, 0x8001, 0x1_2345, 0x7fff_ffff, 0xffff_ffff, 1 << 20] {
                    let (mut a, mut b) = (Vec::new(), ReverseWriter::for_symbols(1));
                    let expected = state_step(state_renormalize(x, f, bits, &mut a), f, g, bits);
                    assert_eq!(symbol.put(x, &mut b), expected, "f={f} x={x:#x}");
                    // The bytes in the order they are written to the output.
                    a.reverse();
                    assert_eq!(a, b.as_slice());
                }
            }
        }
    }

    use super::*;

    // The expected bytes of the order-0, order-1, and PACK tests are htscodecs'
    // (`rans_compress_4x16`), which compresses this input in each mode.
    #[test]
    fn test_encode_order_0() -> io::Result<()> {
        let actual = encode(Flags::empty(), &b"noodles".repeat(30))?;

        let expected = [
            0x00, 0x81, 0x52, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x24, 0x24,
            0x24, 0x24, 0x4c, 0x24, 0xd8, 0x38, 0x1b, 0x00, 0x59, 0x1c, 0x0d, 0x00, 0xa0, 0xca,
            0x03, 0x00, 0x7b, 0xa1, 0x08, 0x00, 0xef, 0xa0, 0x10, 0xd6, 0xd0, 0xb3, 0x23, 0xcf,
            0xdf, 0x1f, 0xdc, 0x98, 0x2c, 0x49, 0x8f, 0x5b, 0x56, 0xca, 0xd6, 0xb5, 0x3d, 0x74,
            0x81, 0x97, 0x91, 0xf9, 0xdf, 0x4b, 0x06, 0x12, 0xff, 0xf5, 0x47, 0x72, 0x56, 0xb7,
            0xd3, 0xec, 0x44, 0x70, 0x23, 0x3a, 0xaf, 0x74, 0xea, 0x7f, 0xe0, 0x66, 0x82, 0x9d,
            0xb0, 0x1b, 0xf6, 0xa0, 0xbc, 0xe5, 0x00, 0x06, 0x40, 0xb4, 0x40, 0x3f, 0x80, 0x28,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER, &b"noodles".repeat(30))?;

        let expected = [
            0x01, 0x81, 0x52, 0xa0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00,
            0x00, 0x00, 0x01, 0x00, 0x01, 0x01, 0x01, 0x01, 0x00, 0x02, 0x20, 0x00, 0x02, 0x00,
            0x05, 0x20, 0x00, 0x01, 0x20, 0x00, 0x03, 0x00, 0x04, 0x20, 0x00, 0x00, 0x00, 0x00,
            0x20, 0x00, 0x02, 0x20, 0x00, 0x00, 0x00, 0x03, 0x20, 0x00, 0x01, 0xaa, 0x01, 0x01,
            0x00, 0x2a, 0x80, 0x00, 0x00, 0x2a, 0x83, 0x00, 0x00, 0x55, 0x02, 0x01, 0x00, 0x00,
            0xaa, 0x00, 0xaa, 0x00, 0xaa, 0x00, 0x54,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    /// As in htscodecs, inputs of 20 bytes or fewer are not striped: the stream is the one the
    /// flags give without STRIPE.
    #[test]
    fn test_encode_stripe_skips_short_inputs() -> io::Result<()> {
        assert_eq!(
            encode(Flags::STRIPE, b"noodles")?,
            encode(Flags::empty(), b"noodles")?
        );
        assert_eq!(
            encode(Flags::STRIPE | Flags::RLE, b"noodles")?,
            encode(Flags::RLE, b"noodles")?
        );

        Ok(())
    }

    #[test]
    fn test_encode_uncompressed() -> io::Result<()> {
        let actual = encode(Flags::CAT, b"noodles")?;

        let expected = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_rle() -> io::Result<()> {
        let actual = encode(Flags::CAT | Flags::RLE, b"noooooooodles")?;

        let expected = [
            0x60, // flags = CAT | RLE
            0x0d, // uncompressed length = 13
            0x07, 0x06, 0x01, 0x6f, 0x07, 0x6e, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_pack() -> io::Result<()> {
        let actual = encode(Flags::PACK, &b"noodles".repeat(30))?;

        let expected = [
            0x80, 0x81, 0x52, 0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x69, 0x04, 0x12, 0x20,
            0x35, 0x43, 0x44, 0x00, 0x51, 0x00, 0x14, 0x12, 0x12, 0x12, 0x12, 0x12, 0x12, 0xab,
            0xea, 0xd3, 0x06, 0x0f, 0x21, 0x06, 0x01, 0x6f, 0x13, 0x1b, 0x01, 0x79, 0x47, 0xfc,
            0x00, 0x72, 0xd8, 0x5a, 0x9d, 0x60, 0xe5, 0xf3, 0x5d, 0xde, 0x23, 0x04, 0x87, 0x9f,
            0x0d, 0x79, 0xa6, 0x55, 0x55, 0x46, 0x3f, 0xf3, 0xea, 0xc3, 0x51, 0x40, 0xb2, 0x00,
            0x94, 0x00, 0x68, 0xc0, 0x9c,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_alphabet() -> io::Result<()> {
        const NUL: u8 = 0x00;

        fn t(src: &[u8], expected: &[u8]) -> io::Result<()> {
            let mut alphabet = [false; ALPHABET_SIZE];

            for &sym in src {
                alphabet[usize::from(sym)] = true;
            }

            let mut actual = Vec::new();
            write_alphabet(&mut actual, &alphabet)?;

            assert_eq!(actual, expected);

            Ok(())
        }

        t(b"abracadabra", &[b'a', b'b', 0x02, b'r', NUL])?;
        t(&[0x01], &[0x01, NUL])?;
        t(&[0xfd, 0xfe, 0xff], &[0xfd, 0xfe, 0x01, NUL])?;

        Ok(())
    }
}
