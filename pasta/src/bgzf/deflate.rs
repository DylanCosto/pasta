use std::io;

#[cfg(feature = "libdeflate")]
thread_local! {
    // libdeflate's decompressor builds large tables when created; reuse one per thread, made on
    // first use so that a failed allocation is an error.
    static DECOMPRESSOR: std::cell::RefCell<Option<crate::libdeflate::Decompressor>> =
        const { std::cell::RefCell::new(None) };
}

#[cfg(feature = "libdeflate")]
pub(crate) fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    DECOMPRESSOR.with(|d| {
        let mut d = d.borrow_mut();
        let d = match &mut *d {
            Some(d) => d,
            empty => empty.insert(crate::libdeflate::Decompressor::try_new()?),
        };
        d.deflate_decompress(src, dst)
            .map(|_| ())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use zlib_rs::{Inflate, InflateFlush, Status};

    const HAS_ZLIB_HEADER: bool = false;
    const WINDOW_BITS: u8 = 15;

    let mut decoder = Inflate::new(HAS_ZLIB_HEADER, WINDOW_BITS);

    let status = decoder
        .decompress(src, dst, InflateFlush::Finish)
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidData))?;

    if status == Status::StreamEnd {
        Ok(())
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
}

#[cfg(feature = "libdeflate")]
pub(crate) fn encode(
    src: &[u8],
    compression_level: libdeflater::CompressionLvl,
    dst: &mut Vec<u8>,
) -> io::Result<u32> {
    use crate::libdeflate::Compressor;

    thread_local! {
        // Compressors allocate their match finders when created; reuse one per thread and level.
        static COMPRESSORS: std::cell::RefCell<Vec<(i32, Compressor)>> =
            const { std::cell::RefCell::new(Vec::new()) };
    }

    let level = i32::from(compression_level);
    let len = COMPRESSORS.with(|c| {
        let mut c = c.borrow_mut();
        let i = match c.iter().position(|(l, _)| *l == level) {
            Some(i) => i,
            None => {
                c.push((level, Compressor::try_new(compression_level)?));
                c.len() - 1
            }
        };
        let encoder = &mut c[i].1;
        let max_len = encoder.deflate_compress_bound(src.len());
        dst.resize(max_len, 0);
        encoder
            .deflate_compress(src, dst)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;

    dst.truncate(len);

    Ok(crc32(src))
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn encode(src: &[u8], compression_level: i32, dst: &mut Vec<u8>) -> io::Result<u32> {
    use zlib_rs::{Deflate, DeflateFlush, Status, compress_bound};

    use crate::bgzf::io::writer::{COMPRESSION_LEVEL_0_OVERHEAD, MAX_BUF_SIZE};

    const COMPRESSION_LEVEL_0: i32 = 0;
    const HAS_ZLIB_HEADER: bool = false;
    const WINDOW_BITS: u8 = 15;
    const MAX_COMPRESSED_SIZE: usize = MAX_BUF_SIZE + COMPRESSION_LEVEL_0_OVERHEAD;

    let max_len = compress_bound(src.len());
    dst.resize(max_len, 0);

    for level in [compression_level, COMPRESSION_LEVEL_0] {
        let mut encoder = Deflate::new(level, HAS_ZLIB_HEADER, WINDOW_BITS);

        let status = encoder
            .compress(src, dst, DeflateFlush::Finish)
            .map_err(|_| io::Error::from(io::ErrorKind::InvalidInput))?;

        if status == Status::StreamEnd {
            let compressed_size = encoder.total_out() as usize;

            if compressed_size > MAX_COMPRESSED_SIZE {
                continue;
            }

            dst.truncate(compressed_size);

            return Ok(crc32(src));
        } else {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }
    }

    // When `src.len() <= MAX_BUF_SIZE`, the output is guaranteed to be <= `MAX_COMPRESSED_SIZE`.
    unreachable!();
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn crc32(src: &[u8]) -> u32 {
    const START: u32 = 0;

    zlib_rs::crc32::crc32(START, src)
}

#[cfg(feature = "libdeflate")]
pub(crate) fn crc32(src: &[u8]) -> u32 {
    let mut crc = libdeflater::Crc::new();
    crc.update(src);
    crc.sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc32() {
        assert_eq!(crc32(b""), 0x00000000);
        assert_eq!(crc32(b"noodles"), 0x802a58a1);
    }
}
