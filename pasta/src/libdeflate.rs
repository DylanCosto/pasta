//! libdeflate's compressor and decompressor, created fallibly.
//!
//! `libdeflater`'s constructors panic when libdeflate cannot allocate its tables (an
//! address-space limit, say); these wrappers return an error instead, and otherwise behave as
//! `libdeflater`'s types do, with its error types.

use std::{ffi::c_void, io, ptr::NonNull};

use libdeflate_sys::{
    libdeflate_alloc_compressor, libdeflate_alloc_decompressor, libdeflate_compressor,
    libdeflate_decompressor, libdeflate_deflate_compress, libdeflate_deflate_compress_bound,
    libdeflate_deflate_decompress, libdeflate_free_compressor, libdeflate_free_decompressor,
    libdeflate_gzip_compress, libdeflate_gzip_compress_bound, libdeflate_gzip_decompress,
    libdeflate_result, libdeflate_result_LIBDEFLATE_INSUFFICIENT_SPACE,
    libdeflate_result_LIBDEFLATE_SUCCESS,
};
use libdeflater::{CompressionError, CompressionLvl, DecompressionError};

/// The error for tables libdeflate could not allocate.
fn out_of_memory(what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::OutOfMemory,
        format!("could not allocate the {what}: out of memory"),
    )
}

/// A libdeflate decompressor.
pub(crate) struct Decompressor(NonNull<libdeflate_decompressor>);

impl Decompressor {
    /// Allocates a decompressor.
    ///
    /// # Errors
    ///
    /// Fails when libdeflate cannot allocate it.
    pub(crate) fn try_new() -> io::Result<Self> {
        // SAFETY: the allocation takes no arguments; a null result is handled below.
        let p = unsafe { libdeflate_alloc_decompressor() };
        NonNull::new(p)
            .map(Self)
            .ok_or_else(|| out_of_memory("deflate decompressor"))
    }

    /// Inflates raw DEFLATE data into `out`, returning the length written.
    pub(crate) fn deflate_decompress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, DecompressionError> {
        let mut written = 0;
        // SAFETY: the decompressor is live, the input and output pointers come from slices of
        // the lengths given, and `written` is a local the call writes to.
        let ret = unsafe {
            libdeflate_deflate_decompress(
                self.0.as_ptr(),
                data.as_ptr().cast::<c_void>(),
                data.len(),
                out.as_mut_ptr().cast::<c_void>(),
                out.len(),
                &mut written,
            )
        };
        decompression_result(ret, written)
    }

    /// Inflates gzip data into `out`, returning the length written.
    pub(crate) fn gzip_decompress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, DecompressionError> {
        let mut written = 0;
        // SAFETY: as in `deflate_decompress`.
        let ret = unsafe {
            libdeflate_gzip_decompress(
                self.0.as_ptr(),
                data.as_ptr().cast::<c_void>(),
                data.len(),
                out.as_mut_ptr().cast::<c_void>(),
                out.len(),
                &mut written,
            )
        };
        decompression_result(ret, written)
    }
}

/// libdeflate's result as `libdeflater` reports it: data that is not DEFLATE, or not of the
/// expected length, is bad data.
fn decompression_result(
    ret: libdeflate_result,
    written: usize,
) -> Result<usize, DecompressionError> {
    #[expect(
        non_upper_case_globals,
        reason = "libdeflate-sys names its constants after the C enum"
    )]
    match ret {
        libdeflate_result_LIBDEFLATE_SUCCESS => Ok(written),
        libdeflate_result_LIBDEFLATE_INSUFFICIENT_SPACE => {
            Err(DecompressionError::InsufficientSpace)
        }
        _ => Err(DecompressionError::BadData),
    }
}

impl Drop for Decompressor {
    fn drop(&mut self) {
        // SAFETY: the pointer came from libdeflate_alloc_decompressor and is freed once.
        unsafe { libdeflate_free_decompressor(self.0.as_ptr()) };
    }
}

/// A libdeflate compressor.
pub(crate) struct Compressor(NonNull<libdeflate_compressor>);

impl Compressor {
    /// Allocates a compressor for `level`.
    ///
    /// # Errors
    ///
    /// Fails when libdeflate cannot allocate it.
    pub(crate) fn try_new(level: CompressionLvl) -> io::Result<Self> {
        // SAFETY: the level is one libdeflate accepts (CompressionLvl checks it); a null result
        // is handled below.
        let p = unsafe { libdeflate_alloc_compressor(i32::from(level)) };
        NonNull::new(p)
            .map(Self)
            .ok_or_else(|| out_of_memory("deflate compressor"))
    }

    /// The largest raw DEFLATE output of `n` bytes.
    pub(crate) fn deflate_compress_bound(&mut self, n: usize) -> usize {
        // SAFETY: the compressor is live.
        unsafe { libdeflate_deflate_compress_bound(self.0.as_ptr(), n) }
    }

    /// Compresses `data` as raw DEFLATE into `out`, returning the length written.
    pub(crate) fn deflate_compress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, CompressionError> {
        // SAFETY: the compressor is live, and the pointers come from slices of the lengths
        // given.
        let n = unsafe {
            libdeflate_deflate_compress(
                self.0.as_ptr(),
                data.as_ptr().cast::<c_void>(),
                data.len(),
                out.as_mut_ptr().cast::<c_void>(),
                out.len(),
            )
        };
        if n == 0 {
            Err(CompressionError::InsufficientSpace)
        } else {
            Ok(n)
        }
    }

    /// The largest gzip output of `n` bytes.
    pub(crate) fn gzip_compress_bound(&mut self, n: usize) -> usize {
        // SAFETY: the compressor is live.
        unsafe { libdeflate_gzip_compress_bound(self.0.as_ptr(), n) }
    }

    /// Compresses `data` as gzip into `out`, returning the length written.
    pub(crate) fn gzip_compress(
        &mut self,
        data: &[u8],
        out: &mut [u8],
    ) -> Result<usize, CompressionError> {
        // SAFETY: as in `deflate_compress`.
        let n = unsafe {
            libdeflate_gzip_compress(
                self.0.as_ptr(),
                data.as_ptr().cast::<c_void>(),
                data.len(),
                out.as_mut_ptr().cast::<c_void>(),
                out.len(),
            )
        };
        if n == 0 {
            Err(CompressionError::InsufficientSpace)
        } else {
            Ok(n)
        }
    }
}

impl Drop for Compressor {
    fn drop(&mut self) {
        // SAFETY: the pointer came from libdeflate_alloc_compressor and is freed once.
        unsafe { libdeflate_free_compressor(self.0.as_ptr()) };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_round_trips() -> io::Result<()> {
        let data = b"ACGT".repeat(1000);
        let mut c = Compressor::try_new(CompressionLvl::default())?;
        let mut d = Decompressor::try_new()?;

        let mut raw = vec![0; c.deflate_compress_bound(data.len())];
        let n = c
            .deflate_compress(&data, &mut raw)
            .map_err(io::Error::other)?;
        let mut out = vec![0; data.len()];
        assert_eq!(d.deflate_decompress(&raw[..n], &mut out), Ok(data.len()));
        assert_eq!(out, data);
        assert_eq!(
            d.deflate_decompress(b"\xff\xff", &mut out),
            Err(DecompressionError::BadData)
        );

        let mut gz = vec![0; c.gzip_compress_bound(data.len())];
        let n = c.gzip_compress(&data, &mut gz).map_err(io::Error::other)?;
        out.fill(0);
        assert_eq!(d.gzip_decompress(&gz[..n], &mut out), Ok(data.len()));
        assert_eq!(out, data);
        Ok(())
    }
}
