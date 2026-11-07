//! zflate's compressors and decompressor, one of each per thread: creating one allocates its
//! tables, so each thread makes a decompressor and a compressor per level on first use and
//! reuses them. A failed allocation (an address-space limit, say) is an I/O error, not an
//! abort.

use std::{cell::RefCell, io};

use zflate::{AllocError, CompressionLevel, Compressor, Decompressor};

thread_local! {
    static DECOMPRESSOR: RefCell<Option<Decompressor>> = const { RefCell::new(None) };
    static COMPRESSORS: RefCell<Vec<Compressor>> = const { RefCell::new(Vec::new()) };
}

/// Runs `f` with this thread's decompressor.
///
/// # Errors
///
/// Fails when the decompressor cannot be allocated.
pub(crate) fn with_decompressor<R>(f: impl FnOnce(&mut Decompressor) -> R) -> io::Result<R> {
    DECOMPRESSOR.with_borrow_mut(|d| {
        let d = match d {
            Some(d) => d,
            empty => empty.insert(
                Decompressor::try_new().map_err(|e| out_of_memory("deflate decompressor", e))?,
            ),
        };
        Ok(f(d))
    })
}

/// Runs `f` with this thread's compressor at `level`.
///
/// # Errors
///
/// Fails when the compressor cannot be allocated.
pub(crate) fn with_compressor<R>(
    level: CompressionLevel,
    f: impl FnOnce(&mut Compressor) -> R,
) -> io::Result<R> {
    COMPRESSORS.with_borrow_mut(|c| {
        let i = match c.iter().position(|c| c.level() == level) {
            Some(i) => i,
            None => {
                c.push(
                    Compressor::try_new(level)
                        .map_err(|e| out_of_memory("deflate compressor", e))?,
                );
                c.len() - 1
            }
        };
        Ok(f(&mut c[i]))
    })
}

/// The error for tables zflate could not allocate.
fn out_of_memory(what: &str, _: AllocError) -> io::Error {
    io::Error::new(
        io::ErrorKind::OutOfMemory,
        format!("could not allocate the {what}: out of memory"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_round_trips() -> io::Result<()> {
        let data = b"ACGT".repeat(1000);
        let mut raw = vec![0; zflate::deflate_compress_bound(data.len())];
        let n = with_compressor(CompressionLevel::DEFAULT, |c| {
            c.deflate_compress(&data, &mut raw)
        })?
        .map_err(io::Error::other)?;
        let mut out = vec![0; data.len()];
        assert_eq!(
            with_decompressor(|d| d.deflate_decompress(&raw[..n], &mut out))?,
            Ok(data.len())
        );
        assert_eq!(out, data);
        assert_eq!(
            with_decompressor(|d| d.deflate_decompress(b"\xff\xff", &mut out))?,
            Err(zflate::DecompressionError::BadData)
        );
        Ok(())
    }
}
