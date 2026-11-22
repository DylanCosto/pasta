mod context;

use std::io;

pub use self::context::{Context, read_context};
use crate::cram::{
    codecs::rans_nx16::ALPHABET_SIZE,
    io::reader::num::{read_u8, read_uint7_as},
};

/// The longest run written as one block of this many bytes, wherever there is room for the
/// block: what it writes past the run is written over by what follows.
const BLOCK_LEN: usize = 32;

/// Expands the runs of `src`: a symbol of the run-length alphabet is followed, in the context,
/// by how many more times it repeats. A run past the end of the output is cut short.
///
/// A symbol without a run is written as a run of one, and a run length of one byte is taken
/// or left without a branch, so that the loop branches on little that the data decides.
pub fn decode(src: &[u8], ctx: &Context<'_>) -> io::Result<Vec<u8>> {
    let mut run_lens = &ctx.src[..];
    let rle_alphabet = read_rle_alphabet(&mut run_lens)?;

    let mut dst = vec![0; ctx.len];
    let mut i = 0;

    for &sym in src {
        let Some(rest) = dst.get_mut(i..).filter(|rest| !rest.is_empty()) else {
            break;
        };

        let has_run = rle_alphabet[usize::from(sym)];

        let run_len = match run_lens.first() {
            Some(&n) if n < 0x80 => {
                let k = usize::from(has_run);
                run_lens = &run_lens[k..];
                usize::from(n) * k
            }
            _ if has_run => read_uint7_as(&mut run_lens)?,
            _ => 0,
        };

        let len = run_len.saturating_add(1);

        match rest.first_chunk_mut::<BLOCK_LEN>() {
            Some(block) if len <= BLOCK_LEN => *block = [sym; BLOCK_LEN],
            _ => {
                let n = len.min(rest.len());
                rest[..n].fill(sym);
            }
        }

        i = i.saturating_add(len);
    }

    if i < dst.len() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Ok(dst)
}

fn read_rle_alphabet(src: &mut &[u8]) -> io::Result<[bool; ALPHABET_SIZE]> {
    let mut alphabet = [false; ALPHABET_SIZE];

    let symbol_count = {
        let n = read_u8(src).map(usize::from)?;
        if n == 0 { ALPHABET_SIZE } else { n }
    };

    for _ in 0..symbol_count {
        let sym = read_u8(src)?;
        alphabet[usize::from(sym)] = true;
    }

    Ok(alphabet)
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use super::*;

    fn context(src: &[u8], len: usize) -> Context<'_> {
        Context {
            src: Cow::from(src),
            len,
        }
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        // One symbol with runs, `o`, which repeats 7 more times.
        assert_eq!(
            decode(b"nodles", &context(&[0x01, b'o', 0x07], 13))?,
            b"noooooooodles"
        );

        // Two, with a run of 1 + 300 (0x82 0x2c) and runs of 1 + 0.
        let mut expected = vec![b'a'; 301];
        expected.extend(b"bcab");
        assert_eq!(
            decode(
                b"abcab",
                &context(&[0x02, b'a', b'b', 0x82, 0x2c, 0x00, 0x00, 0x00], 305)
            )?,
            expected
        );

        // No symbol count is every symbol.
        let mut src = vec![0x00];
        src.extend(0..=u8::MAX);
        src.extend([0x02, 0x01]);
        assert_eq!(decode(b"xy", &context(&src, 5))?, b"xxxyy");

        assert!(decode(b"", &context(&[0x01, b'o'], 0))?.is_empty());

        Ok(())
    }

    #[test]
    fn test_decode_with_a_run_past_the_end() -> io::Result<()> {
        assert_eq!(decode(b"no", &context(&[0x01, b'o', 0x07], 4))?, b"nooo");
        Ok(())
    }

    #[test]
    fn test_decode_with_too_few_symbols_or_run_lengths() {
        for (src, ctx) in [
            (&b"nodl"[..], context(&[0x01, b'o', 0x07], 13)),
            (&b"nodles"[..], context(&[0x01, b'o'], 13)),
            (&b"nodles"[..], context(&[0x01, b'o', 0x82], 13)),
            (&b"nodles"[..], context(&[0x02, b'o'], 13)),
        ] {
            assert!(matches!(
                decode(src, &ctx),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
            ));
        }
    }
}
