//! FASTA filesystem operations.

mod index;

use std::{io, path::Path};

pub use self::index::index;
use super::fai;

/// Reads an associated FASTA index.
///
/// This attempts to read an associated index at `<src>.fai`.
///
/// # Examples
///
/// ```no_run
/// use pasta::fasta;
/// let index = fasta::fs::read_associated_index("src.fa")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read_associated_index<P>(src: P) -> io::Result<fai::Index>
where
    P: AsRef<Path>,
{
    const FAI_EXT: &str = "fai";

    let fai_src = src.as_ref().with_added_extension(FAI_EXT);
    fai::fs::read(fai_src)
}
