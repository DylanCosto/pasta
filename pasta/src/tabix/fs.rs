//! tabix filesystem operations.

use std::{fs::File, io, path::Path};

use super::{
    Index,
    io::{Reader, Writer},
};

/// Reads the entire contents of a tabix index.
///
/// This is a convenience function and is equivalent to opening the file at the given path and
/// reading the index.
///
/// # Examples
///
/// ```no_run
/// use pasta::tabix as tabix;
/// let index = tabix::fs::read("sample.vcf.gz.tbi")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn read<P>(src: P) -> io::Result<Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_index()
}

/// Writes a tabix index to a file.
///
/// This is a convenience function and is equivalent to creating a file at the given path and
/// writing the index.
///
/// # Examples
///
/// ```no_run
/// use pasta::csi::binning_index::index::Header;
/// use pasta::tabix as tabix;
/// let index = tabix::Index::builder().set_header(Header::default()).build();
/// tabix::fs::write("sample.vcf.gz.tbi", &index)?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn write<P>(dst: P, index: &Index) -> io::Result<()>
where
    P: AsRef<Path>,
{
    let mut writer = File::create(dst).map(Writer::new)?;
    writer.write_index(index)?;
    writer.try_finish()
}

// The only test writes to `/dev/full`, which macOS and the BSDs lack.
#[cfg(all(test, target_os = "linux"))]
mod tests {
    use super::*;

    #[test]
    fn test_write_reports_a_failed_final_write() {
        // `/dev/full` fails every write, and the whole index is still buffered when `write`
        // returns.
        let index = Index::builder()
            .set_header(crate::csi::binning_index::index::Header::default())
            .build();
        let result = write("/dev/full", &index);
        assert_eq!(
            result.map_err(|e| e.kind()),
            Err(io::ErrorKind::StorageFull)
        );
    }
}
