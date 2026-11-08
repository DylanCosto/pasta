pub mod context;

use std::io;

pub use self::context::{Context, build_context, write_context};
use crate::cram::io::writer::num::write_uint7;

pub fn encode(src: &[u8], ctx: &mut Context) -> io::Result<Vec<u8>> {
    let mut dst = Vec::with_capacity(src.len());
    let mut rest = src;

    while let Some((&sym, tail)) = rest.split_first() {
        dst.push(sym);
        rest = tail;

        if ctx.alphabet[usize::from(sym)] {
            let len = rest.iter().position(|&s| s != sym).unwrap_or(rest.len());
            rest = &rest[len..];

            let n =
                u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            write_uint7(&mut ctx.dst, n)?;
        }
    }

    Ok(dst)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> Result<(), Box<dyn std::error::Error>> {
        let src = b"nnndlllllss";
        let mut ctx = build_context(src)?;
        let buf = encode(src, &mut ctx)?;

        assert_eq!(buf, [b'n', b'd', b'l', b's', b's']);

        assert_eq!(
            ctx.dst,
            [
                0x02, // symbol count
                b'l', b'n', // alphabet
                0x02, 0x04, // run-lengths
            ]
        );

        Ok(())
    }
}
