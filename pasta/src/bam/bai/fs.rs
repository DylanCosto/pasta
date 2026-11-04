//! BAI filesystem operations.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::Path,
};

use super::{
    Index,
    io::{Reader, Writer},
};

/// Reads the entire contents of a BAM index.
///
/// This is a convenience function and is equivalent to opening the file at the given path, reading
/// the header, and reading the index.
///
/// # Examples
///
/// ```no_run
/// use pasta::bam::bai;
/// let index = bai::fs::read("sample.bam.bai")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}

/// Writes a BAM index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path, writing
/// the header, and writing the index.
///
/// # Examples
///
/// ```no_run
/// use pasta::bam::bai;
/// let index = bai::Index::default();
/// bai::fs::write("sample.bam.bai", &index)?;
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
