//! gzip index filesystem operations.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::Path,
};

use super::{
    Index,
    io::{Reader, Writer},
};

/// Reads the entire contents of a GZ index.
///
/// This is a convenience function and is equivalent to opening the given path and reading the
/// index.
///
/// # Examples
///
/// ```no_run
/// use pasta::bgzf::gzi;
/// let index = gzi::fs::read("in.gz.gzi")?;
/// # Ok::<_, std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}

/// Writes a GZ index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// use pasta::bgzf::gzi;
/// let index = gzi::Index::default();
/// gzi::fs::write("in.gz.gzi", &index)?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(BufWriter::new).map(Writer::new)?;
    writer.write_index(index)?;
    writer.get_mut().flush()
}

// The only test writes to `/dev/full`, which macOS and the BSDs lack.
#[cfg(all(test, target_os = "linux"))]
mod tests {
    use super::*;

    #[test]
    fn test_write_reports_a_failed_final_write() {
        // `/dev/full` fails every write, and the whole index is still buffered when `write`
        // returns.
        let result = write("/dev/full", &Index::default());
        assert_eq!(
            result.map_err(|e| e.kind()),
            Err(io::ErrorKind::StorageFull)
        );
    }
}
