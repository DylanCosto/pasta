use std::io;

use flate2::Compression;

use crate::deflate::{with_compressor, with_decompressor};

/// Inflates `src` into `dst`, which it must fill exactly, as the block's size says (HTSlib
/// rejects a block that inflates to less).
pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    match with_decompressor(|d| d.gzip_decompress(src, dst))? {
        Ok(n) if n == dst.len() => Ok(()),
        Ok(_) => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "gzip block size mismatch",
        )),
        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

pub fn encode(compression_level: Compression, src: &[u8]) -> io::Result<Vec<u8>> {
    let level = i32::try_from(compression_level.level())
        .map(libdeflate_level)
        .ok()
        .and_then(|n| u8::try_from(n).ok())
        .and_then(zflate::CompressionLevel::new)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid compression level"))?;

    let mut dst = vec![0; zflate::gzip_compress_bound(src.len())];
    let len = with_compressor(level, |c| c.gzip_compress(src, &mut dst))?
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.truncate(len);

    Ok(dst)
}

/// Maps a zlib compression level to libdeflate's 1–12 scale as HTSlib does
/// (`libdeflate_deflate`): 1 to 4 stay, and 5, 6, 7, 8 and 9 become 6, 7, 9, 10 and 12.
fn libdeflate_level(level: i32) -> i32 {
    let level = if level > 0 { level } else { 6 };
    let level = level * 123 / 100;
    (level + i32::from(level >= 8)).min(12)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_libdeflate_level() {
        let actual: Vec<_> = (0..=9).map(libdeflate_level).collect();
        assert_eq!(actual, [7, 1, 2, 3, 4, 6, 7, 9, 10, 12]);
    }
}
