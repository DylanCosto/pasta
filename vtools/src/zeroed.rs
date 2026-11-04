//! Zeroed tables sized by command-line options (`stats -d`, `stats -u`), allocated so that a size
//! the system cannot give is reported as an error instead of aborting the program.

use std::{
    alloc::{self, Layout},
    fmt,
};

/// A table that the system could not allocate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AllocError {
    bytes: usize,
}

impl fmt::Display for AllocError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "could not allocate {} bytes of memory", self.bytes)
    }
}

impl std::error::Error for AllocError {}

/// Returns `n` zeroed counts, as `vec![0u64; n]` does: the pages are zeroed by the system as they
/// are first used, so a large table that is mostly unused costs little.
///
/// # Errors
///
/// Fails when the system cannot give the memory (`ulimit -v`, say), or `n` counts do not fit the
/// address space.
pub fn counts(n: usize) -> Result<Vec<u64>, AllocError> {
    let too_large = AllocError {
        bytes: n.saturating_mul(size_of::<u64>()),
    };
    let layout = Layout::array::<u64>(n).map_err(|_| too_large)?;
    if layout.size() == 0 {
        return Ok(Vec::new());
    }
    // SAFETY: the layout has a nonzero size, checked just above.
    let ptr = unsafe { alloc::alloc_zeroed(layout) }.cast::<u64>();
    if ptr.is_null() {
        return Err(too_large);
    }
    // SAFETY: `ptr` was allocated by the global allocator with the layout of `n` u64s, which is
    // what a `Vec<u64>` of capacity `n` uses, and all its `n` values are zero, a valid u64.
    Ok(unsafe { Vec::from_raw_parts(ptr, n, n) })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_are_zeroed() {
        assert_eq!(counts(0), Ok(Vec::new()));
        let mut v = counts(1000).unwrap();
        assert!(v.iter().all(|&c| c == 0));
        v.push(1);
        assert_eq!(v.len(), 1001);
    }

    #[test]
    fn an_impossible_size_is_an_error() {
        let e = counts(usize::MAX / 2).unwrap_err();
        assert!(e.to_string().starts_with("could not allocate"), "{e}");
    }
}
