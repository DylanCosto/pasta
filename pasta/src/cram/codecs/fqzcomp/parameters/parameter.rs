mod flags;

use std::{io, num::NonZero};

pub use self::flags::Flags;
use self::flags::read_flags;
use super::read_array;
use crate::cram::io::reader::num::{read_u8, read_u16_le};

// § 6.2 "FQZComp Data Stream" (2023-03-15).
const QUALITIES_TABLE_SIZE: usize = 256;
const POSITIONS_TABLE_SIZE: usize = 1024;
const DELTAS_TABLE_SIZE: usize = 256;

/// A parameter set, with its tables as the decoder's context takes them.
pub struct Parameter {
    pub context: u16,
    flags: Flags,
    symbol_count: NonZero<usize>,
    /// The quality bits the context keeps, as a mask.
    pub q_mask: u32,
    pub q_shift: u8,
    pub q_loc: u8,
    pub s_loc: u8,
    /// The score each symbol stands for: itself, without a quality map.
    pub quality_map: [u8; QUALITIES_TABLE_SIZE],
    pub qualities_table: [u32; QUALITIES_TABLE_SIZE],
    /// Shifted into place in the context, and zero without a positions table.
    pub positions_table: [u32; POSITIONS_TABLE_SIZE],
    /// Shifted into place in the context, and zero without a deltas table.
    pub deltas_table: [u32; DELTAS_TABLE_SIZE],
}

impl Parameter {
    pub fn flags(&self) -> Flags {
        self.flags
    }

    pub fn symbol_count(&self) -> NonZero<usize> {
        self.symbol_count
    }
}

pub fn fqz_decode_single_param(src: &mut &[u8]) -> io::Result<Parameter> {
    let context = read_u16_le(src)?;
    let flags = read_flags(src)?;
    let max_symbol = read_max_symbol(src)?;

    let (q_bits, q_shift) = read_u4x2(src)?;
    let (q_loc, s_loc) = read_u4x2(src)?;
    let (p_loc, d_loc) = read_u4x2(src)?;

    let mut quality_map = std::array::from_fn(|i| i as u8);

    if flags.has_quality_map() {
        let symbols = read_quality_map(src, max_symbol)?;
        quality_map[..symbols.len()].copy_from_slice(symbols);
    }

    let mut qualities_table = std::array::from_fn(|i| i as u32);

    if flags.has_qualities_table() {
        read_table(src, 0, &mut qualities_table)?;
    }

    let mut positions_table = [0; POSITIONS_TABLE_SIZE];

    if flags.has_positions_table() {
        read_table(src, p_loc, &mut positions_table)?;
    }

    let mut deltas_table = [0; DELTAS_TABLE_SIZE];

    if flags.has_deltas_table() {
        read_table(src, d_loc, &mut deltas_table)?;
    }

    Ok(Parameter {
        context,
        flags,
        symbol_count: max_to_count(max_symbol),
        q_mask: (1 << q_bits) - 1,
        q_shift,
        q_loc,
        s_loc,
        quality_map,
        qualities_table,
        positions_table,
        deltas_table,
    })
}

fn read_u4x2(src: &mut &[u8]) -> io::Result<(u8, u8)> {
    let n = read_u8(src)?;
    Ok((n >> 4, n & 0x0f))
}

fn read_max_symbol(src: &mut &[u8]) -> io::Result<u8> {
    read_u8(src)
}

fn read_quality_map<'a>(src: &mut &'a [u8], max_symbol: u8) -> io::Result<&'a [u8]> {
    let len = usize::from(max_symbol);

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    Ok(buf)
}

/// Reads a run-length encoded table into `table`, each value shifted left by `shift`.
fn read_table(src: &mut &[u8], shift: u8, table: &mut [u32]) -> io::Result<()> {
    let values = read_array(src, table.len())?;

    for (entry, value) in table.iter_mut().zip(values) {
        *entry = u32::from(value) << shift;
    }

    Ok(())
}

fn max_to_count(n: u8) -> NonZero<usize> {
    let m = usize::from(n) + 1;
    // SAFETY: `m > 0`.
    NonZero::new(m).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_u4x2() -> io::Result<()> {
        assert_eq!(read_u4x2(&mut &[0x0f][..])?, (0x00, 0x0f));
        Ok(())
    }
}
