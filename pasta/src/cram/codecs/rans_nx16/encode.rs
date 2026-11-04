pub mod bit_pack;
mod order_0;
mod order_1;
mod rle;
mod stripe;

use std::{
    borrow::Cow,
    io::{self, Write},
};

use super::{ALPHABET_SIZE, Flags};
use crate::cram::io::writer::num::{write_u8, write_u32_le, write_uint7};

// § 3 "rANS Nx16" (2023-03-15): "The lower-bound and initial encoder state _L_ is [...] 0x8000."
const LOWER_BOUND: u32 = 0x8000;

pub fn encode(mut flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut src = Cow::from(src);
    let mut dst = Vec::new();

    write_flags(&mut dst, flags)?;

    if flags.has_uncompressed_size() {
        write_uncompressed_size(&mut dst, src.len())?;
    }

    let state_count = flags.state_count();

    if flags.is_striped() {
        let buf = stripe::encode(&src)?;
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

    if src.len() < state_count {
        flags.remove(Flags::ORDER);
        flags.insert(Flags::CAT);
        dst[0] = u8::from(flags);
    }

    if flags.is_uncompressed() {
        dst.write_all(&src)?;
    } else if flags.order() == 0 {
        let ctx = order_0::build_context(&src);
        order_0::write_context(&mut dst, &ctx)?;
        order_0::encode(&src, &ctx, state_count, &mut dst)?;
    } else {
        let ctx = order_1::build_context(&src, state_count);
        order_1::write_context(&mut dst, &ctx)?;
        order_1::encode(&src, &ctx, state_count, &mut dst)?;
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
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct EncodingSymbol {
    x_max: u32,
    rcp_freq: u32,
    rcp_shift: u32,
    bias: u32,
    cmpl_freq: u32,
}

impl EncodingSymbol {
    /// The step for a symbol with frequency `f` (nonzero) and cumulative frequency `g`, with
    /// frequencies summing to `1 << bits`.
    pub(super) fn new(f: u32, g: u32, bits: u32) -> Self {
        let x_max = (1 << (31 - bits)) * f;
        let cmpl_freq = (1 << bits) - f;

        if f < 2 {
            Self {
                x_max,
                rcp_freq: u32::MAX,
                rcp_shift: 0,
                bias: g + (1 << bits) - 1,
                cmpl_freq,
            }
        } else {
            let shift = 32 - (f - 1).leading_zeros();
            let rcp_freq = (1u64 << (shift + 31)).div_ceil(u64::from(f));

            Self {
                x_max,
                rcp_freq: rcp_freq as u32,
                rcp_shift: shift - 1,
                bias: g,
                cmpl_freq,
            }
        }
    }

    /// Renormalizes `x` (writing to `dst`) and encodes the symbol.
    #[inline]
    pub(super) fn put(&self, mut x: u32, dst: &mut Vec<u8>) -> u32 {
        while x >= self.x_max {
            write_u16_be(dst, (x & 0xffff) as u16);
            x >>= 16;
        }

        let q = ((u64::from(x) * u64::from(self.rcp_freq)) >> 32) as u32 >> self.rcp_shift;
        x + self.bias + q * self.cmpl_freq
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
                    let (mut a, mut b) = (Vec::new(), Vec::new());
                    let expected = state_step(state_renormalize(x, f, bits, &mut a), f, g, bits);
                    assert_eq!(symbol.put(x, &mut b), expected, "f={f} x={x:#x}");
                    assert_eq!(a, b);
                }
            }
        }
    }

    use super::*;

    #[test]
    fn test_encode_order_0() -> io::Result<()> {
        let actual = encode(Flags::empty(), b"noodles")?;

        let expected = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x84, 0x49, 0x84, 0x49, 0x84,
            0x49, 0x84, 0x49, 0x89, 0x13, 0x84, 0x49, 0x1b, 0xa7, 0x18, 0x00, 0xe9, 0x4a, 0x0c,
            0x00, 0x31, 0x6d, 0x0c, 0x00, 0x08, 0x80, 0x03, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER, b"noodles")?;

        let expected = [
            0x01, 0x07, 0xc0, 0x00, 0x64, 0x65, 0x00, 0x6c, 0x6e, 0x6f, 0x00, 0x73, 0x00, 0x00,
            0x00, 0x88, 0x00, 0x00, 0x01, 0x88, 0x00, 0x90, 0x00, 0x00, 0x00, 0x00, 0x02, 0xa0,
            0x00, 0x00, 0x02, 0x00, 0x05, 0xa0, 0x00, 0x00, 0x01, 0xa0, 0x00, 0x00, 0x03, 0x00,
            0x04, 0xa0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x90, 0x00, 0x00, 0x02, 0x90, 0x00, 0x00,
            0x00, 0x00, 0x06, 0x00, 0x04, 0x02, 0x00, 0x00, 0x08, 0x01, 0x00, 0x00, 0x08, 0x01,
            0x00, 0x00, 0x00, 0x02, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_stripe() -> io::Result<()> {
        let actual = encode(Flags::STRIPE, b"noodles")?;

        let expected = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, // chunk count = 4
            0x03, 0x03, 0x03, 0x02, // compressed sizes = [3, 3, 3, 2]
            0x30, 0x6e, 0x6c, // chunks[0]
            0x30, 0x6f, 0x65, // chunks[1]
            0x30, 0x6f, 0x73, // chunks[2]
            0x30, 0x64, // chunks[3]
        ];

        assert_eq!(actual, expected);

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
        let actual = encode(Flags::PACK, b"noodles")?;

        let expected = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x04, 0x05, 0x00, 0x12, 0x43, 0x00,
            0x88, 0x00, 0x88, 0x00, 0x88, 0x00, 0x88, 0x00, 0x00, 0x0c, 0x02, 0x00, 0x00, 0x00,
            0x02, 0x00, 0x00, 0x08, 0x02, 0x00, 0x00, 0x04, 0x02, 0x00,
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
