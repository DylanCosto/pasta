use std::io;

use flate2::Compression;

#[cfg(feature = "libdeflate")]
pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use crate::libdeflate::Decompressor;

    let mut decoder = Decompressor::try_new()?;

    decoder
        .gzip_decompress(src, dst)
        .map(|_| ())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(not(feature = "libdeflate"))]
pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use std::io::Read;

    use flate2::bufread::GzDecoder;

    let mut decoder = GzDecoder::new(src);
    decoder.read_exact(dst)
}

#[cfg(feature = "libdeflate")]
pub fn encode(compression_level: Compression, src: &[u8]) -> io::Result<Vec<u8>> {
    use libdeflater::CompressionLvl;

    use crate::libdeflate::Compressor;

    let lvl = i32::try_from(compression_level.level())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .map(libdeflate_level)
        .and_then(|n| {
            CompressionLvl::new(n).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid libdeflate compression level",
                )
            })
        })?;

    let mut encoder = Compressor::try_new(lvl)?;

    let max_len = encoder.gzip_compress_bound(src.len());
    let mut dst = vec![0; max_len];

    let len = encoder
        .gzip_compress(src, &mut dst)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    dst.resize(len, 0);

    Ok(dst)
}

/// Maps a zlib compression level to libdeflate's 1–12 scale as HTSlib does
/// (`libdeflate_deflate`): 1 to 4 stay, and 5, 6, 7, 8 and 9 become 6, 7, 9, 10 and 12.
#[cfg(feature = "libdeflate")]
fn libdeflate_level(level: i32) -> i32 {
    let level = if level > 0 { level } else { 6 };
    let level = level * 123 / 100;
    (level + i32::from(level >= 8)).min(12)
}

#[cfg(not(feature = "libdeflate"))]
pub fn encode(compression_level: Compression, src: &[u8]) -> io::Result<Vec<u8>> {
    use std::io::Write;

    use flate2::write::GzEncoder;

    let mut encoder = GzEncoder::new(Vec::new(), compression_level);
    encoder.write_all(src)?;
    encoder.finish()
}

#[cfg(all(test, feature = "libdeflate"))]
mod tests {
    use super::*;

    #[test]
    fn test_libdeflate_level() {
        let actual: Vec<_> = (0..=9).map(libdeflate_level).collect();
        assert_eq!(actual, [7, 1, 2, 3, 4, 6, 7, 9, 10, 12]);
    }
}
