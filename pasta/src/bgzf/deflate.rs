use std::io;

use super::io::writer::CompressionLevel;
use crate::deflate::{with_compressor, with_decompressor};

/// Inflates `src` into `dst`, which it must fill exactly, as its block's ISIZE says.
pub(crate) fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    let result = with_decompressor(|d| d.deflate_decompress(src, dst))?;
    filled(result, dst.len())
}

/// A decompression's result as an I/O error, and an error unless it filled the `len` bytes
/// expected: a short block would leave the end of its buffer unwritten.
fn filled(result: Result<usize, zflate::DecompressionError>, len: usize) -> io::Result<()> {
    match result {
        Ok(n) if n == len => Ok(()),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data size mismatch",
        )),
        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

/// [`decode`] for several streams at once, each into its `dst`, with zflate's batch decoder,
/// which on 64-bit Arm decodes two at a time: each result is what `decode` gives that stream.
///
/// # Errors
///
/// Fails, for all of them, when the decompressor cannot be allocated.
pub(crate) fn decode_batch(
    src: &[&[u8]],
    dst: &mut [&mut [u8]],
) -> io::Result<Vec<io::Result<()>>> {
    let mut results = vec![Ok(0); src.len()];
    with_decompressor(|d| d.deflate_decompress_batch(src, dst, &mut results))?;
    Ok(results
        .into_iter()
        .zip(dst.iter())
        .map(|(r, d)| filled(r, d.len()))
        .collect())
}

pub(crate) fn encode(
    src: &[u8],
    compression_level: CompressionLevel,
    dst: &mut Vec<u8>,
) -> io::Result<u32> {
    dst.resize(zflate::deflate_compress_bound(src.len()), 0);
    let len = with_compressor(compression_level.to_zflate(), |c| {
        c.deflate_compress(src, dst)
    })?
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.truncate(len);

    Ok(crc32(src))
}

pub(crate) fn crc32(src: &[u8]) -> u32 {
    zflate::crc32(0, src)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc32() {
        assert_eq!(crc32(b""), 0x00000000);
        assert_eq!(crc32(b"noodles"), 0x802a58a1);
    }

    #[test]
    fn test_decode_must_fill_its_buffer() -> io::Result<()> {
        let data = b"noodles-bgzf ".repeat(100);
        let mut compressed = Vec::new();
        encode(&data, CompressionLevel::default(), &mut compressed)?;
        let mut out = vec![0; data.len() + 1];
        let err = decode(&compressed, &mut out).expect_err("one byte short of the buffer");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        let mut outs = [vec![0; data.len()], vec![0; data.len() + 1]];
        let mut dsts: Vec<&mut [u8]> = outs.iter_mut().map(|o| &mut o[..]).collect();
        let results = decode_batch(&[&compressed, &compressed], &mut dsts)?;
        assert!(results[0].is_ok());
        assert_eq!(
            results[1].as_ref().map_err(io::Error::kind),
            Err(io::ErrorKind::InvalidData)
        );
        Ok(())
    }

    #[test]
    fn test_encode_decode() -> io::Result<()> {
        let data = b"noodles-bgzf ".repeat(1000);
        for level in [0, 1, 5, 6, 12] {
            let level = CompressionLevel::new(level).expect("valid level");
            let mut compressed = Vec::new();
            assert_eq!(encode(&data, level, &mut compressed)?, crc32(&data));
            let mut out = vec![0; data.len()];
            decode(&compressed, &mut out)?;
            assert_eq!(out, data);
        }
        Ok(())
    }
}
