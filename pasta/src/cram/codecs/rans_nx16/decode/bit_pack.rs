mod context;

use std::io;

use self::context::Context;
pub use self::context::read_context;

/// The most symbols a bit-packed stream has: 4 bits each.
const MAX_SYMBOL_COUNT: usize = 16;

pub fn decode(src: &[u8], ctx: &Context<'_>) -> io::Result<Vec<u8>> {
    let symbols_per_byte = match ctx.symbol_count.get() {
        1 => 0,
        2 => 8,
        3..=4 => 4,
        5..=MAX_SYMBOL_COUNT => 2,
        n => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected bit pack symbol count to be <= {MAX_SYMBOL_COUNT}, got {n}"),
            ));
        }
    };

    // As in htscodecs, a code past the symbols is 0.
    let mut mapping_table = [0; MAX_SYMBOL_COUNT];

    for (d, s) in mapping_table.iter_mut().zip(ctx.mapping_table) {
        *d = *s;
    }

    if symbols_per_byte == 0 {
        return Ok(vec![mapping_table[0]; ctx.uncompressed_size]);
    }

    if src.len() < ctx.uncompressed_size.div_ceil(symbols_per_byte) {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let mut dst = vec![0; ctx.uncompressed_size];
    unpack(src, &mapping_table, symbols_per_byte, &mut dst);

    Ok(dst)
}

fn unpack(src: &[u8], mapping_table: &[u8; MAX_SYMBOL_COUNT], chunk_size: usize, dst: &mut [u8]) {
    const BITS: usize = u8::BITS as usize;

    let shift = BITS / chunk_size;
    let mask = (1 << shift) - 1;

    for (mut s, chunk) in src.iter().copied().zip(dst.chunks_mut(chunk_size)) {
        for d in chunk {
            *d = mapping_table[usize::from(s & mask)];
            s >>= shift;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unpack() {
        fn mapping_table(symbols: &[u8]) -> [u8; MAX_SYMBOL_COUNT] {
            let mut table = [0; MAX_SYMBOL_COUNT];
            table[..symbols.len()].copy_from_slice(symbols);
            table
        }

        let mut dst = [0; 4];
        unpack(&[0b00001100], &mapping_table(b"nd"), 8, &mut dst);
        assert_eq!(&dst, b"nndd");

        let mut dst = [0; 4];
        unpack(&[0b11100100], &mapping_table(b"ndls"), 4, &mut dst);
        assert_eq!(&dst, b"ndls");

        let mut dst = [0; 2];
        unpack(&[0b00100000], &mapping_table(b"nodles"), 2, &mut dst);
        assert_eq!(&dst, b"nd");
    }

    fn context(mapping_table: &[u8], uncompressed_size: usize) -> Context<'_> {
        Context {
            symbol_count: std::num::NonZero::new(mapping_table.len()).unwrap(),
            mapping_table,
            uncompressed_size,
        }
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        assert_eq!(decode(&[], &context(b"n", 3))?, b"nnn");
        assert_eq!(decode(&[0b00001100], &context(b"nd", 4))?, b"nndd");
        assert_eq!(decode(&[0b11100100, 0b01], &context(b"ndls", 5))?, b"ndlsd");
        assert_eq!(decode(&[0x20, 0x45], &context(b"nodles", 4))?, b"ndse");

        Ok(())
    }

    #[test]
    fn test_decode_with_a_code_past_the_symbols() -> io::Result<()> {
        // As in htscodecs, such a code is 0.
        assert_eq!(decode(&[0b11100100], &context(b"ndl", 4))?, b"ndl\0");
        assert_eq!(decode(&[0xf1], &context(b"nodle", 2))?, b"o\0");

        Ok(())
    }

    #[test]
    fn test_decode_with_too_few_bytes() {
        for (src, ctx) in [
            (&[0x00][..], context(b"nd", 9)),
            (&[0x00, 0x00][..], context(b"ndls", 9)),
            (&[][..], context(b"nodles", 1)),
        ] {
            assert!(matches!(
                decode(src, &ctx),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
            ));
        }
    }

    #[test]
    fn test_decode_with_too_many_symbols() {
        let mapping_table = [0; 17];

        assert!(matches!(
            decode(&[0x00], &context(&mapping_table, 1)),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }
}
