//! FAI filesystem operations.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Write},
    path::Path,
};

use super::{
    Index,
    io::{Reader, Writer},
};

/// Reads the entire contents of a FASTA index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// parsing each record.
///
/// # Examples
///
/// ```no_run
/// use pasta::fasta::fai;
/// let index = fai::fs::read("reference.fa.fai")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(BufReader::new).map(Reader::new)?;
    reader.read_index()
}

/// Writes a FASTA index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// use pasta::fasta::fai;
/// let index = fai::Index::default();
/// fai::fs::write("reference.fa.fai", &index)?;
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
        let line_bases = std::num::NonZero::new(80).expect("80 is not zero");
        let line_width = std::num::NonZero::new(81).expect("81 is not zero");
        let record = crate::fasta::fai::Record::new("sq0", 13, 5, line_bases, line_width);
        let index = Index::from(vec![record]);
        let result = write("/dev/full", &index);
        assert_eq!(
            result.map_err(|e| e.kind()),
            Err(io::ErrorKind::StorageFull)
        );
    }
}
