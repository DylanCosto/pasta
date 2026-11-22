mod context;

use std::io;

use self::context::Context;
pub use self::context::read_context;

/// The most symbols a bit-packed stream has: 4 bits each.
const MAX_SYMBOL_COUNT: usize = 16;

pub fn decode(src: &[u8], ctx: &Context<'_>) -> io::Result<Vec<u8>> {
    // As in htscodecs, a code past the symbols is 0.
    let mut mapping_table = [0; MAX_SYMBOL_COUNT];

    for (d, s) in mapping_table.iter_mut().zip(ctx.mapping_table) {
        *d = *s;
    }

    match ctx.symbol_count.get() {
        1 => Ok(vec![mapping_table[0]; ctx.uncompressed_size]),
        2 => unpack::<8>(src, &mapping_table, ctx.uncompressed_size),
        3..=4 => unpack::<4>(src, &mapping_table, ctx.uncompressed_size),
        5..=MAX_SYMBOL_COUNT => unpack::<2>(src, &mapping_table, ctx.uncompressed_size),
        n => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("expected bit pack symbol count to be <= {MAX_SYMBOL_COUNT}, got {n}"),
        )),
    }
}

/// Unpacks `len` symbols, `N` to a byte from its low bits.
///
/// A byte's symbols are looked up together, in a table of what each of the 256 bytes unpacks
/// to.
fn unpack<const N: usize>(
    src: &[u8],
    mapping_table: &[u8; MAX_SYMBOL_COUNT],
    len: usize,
) -> io::Result<Vec<u8>> {
    let bits = u8::BITS as usize / N;
    let mask = (1 << bits) - 1;

    if src.len() < len.div_ceil(N) {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let table: [[u8; N]; 256] =
        std::array::from_fn(|b| std::array::from_fn(|i| mapping_table[(b >> (i * bits)) & mask]));

    let mut dst = vec![0; len];
    let (chunks, tail) = dst.as_chunks_mut::<N>();

    for (chunk, &b) in chunks.iter_mut().zip(src) {
        *chunk = table[usize::from(b)];
    }

    if let Some(&b) = src.get(chunks.len()) {
        tail.copy_from_slice(&table[usize::from(b)][..tail.len()]);
    }

    Ok(dst)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unpack() -> io::Result<()> {
        fn mapping_table(symbols: &[u8]) -> [u8; MAX_SYMBOL_COUNT] {
            let mut table = [0; MAX_SYMBOL_COUNT];
            table[..symbols.len()].copy_from_slice(symbols);
            table
        }

        assert_eq!(
            unpack::<8>(&[0b00001100], &mapping_table(b"nd"), 4)?,
            b"nndd"
        );
        assert_eq!(
            unpack::<8>(&[0b10001100, 0b101], &mapping_table(b"nd"), 11)?,
            b"nnddnnnddnd"
        );

        assert_eq!(
            unpack::<4>(&[0b11100100], &mapping_table(b"ndls"), 4)?,
            b"ndls"
        );
        assert_eq!(
            unpack::<4>(&[0b11100100, 0b0110], &mapping_table(b"ndls"), 6)?,
            b"ndlsld"
        );

        assert_eq!(
            unpack::<2>(&[0b00100000], &mapping_table(b"nodles"), 2)?,
            b"nd"
        );
        assert_eq!(
            unpack::<2>(&[0x20, 0x45, 0x03], &mapping_table(b"nodles"), 5)?,
            b"ndsel"
        );

        // Bytes past the symbols are not read.
        assert_eq!(
            unpack::<4>(&[0b11100100, 0xff], &mapping_table(b"ndls"), 3)?,
            b"ndl"
        );
        assert!(unpack::<2>(&[], &mapping_table(b"nodles"), 0)?.is_empty());

        Ok(())
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
