//! Indexed FASTA files opened as HTSlib's `fai_load` opens them: plain or BGZF-compressed,
//! with the `.fai` index (and a BGZF file's `.gzi` index) built and saved when missing.

use std::{
    fmt,
    fs::File,
    io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    os::unix::fs::FileExt,
    path::{Path, PathBuf},
};

use pasta::{
    bgzf::{self, gzi},
    fasta::{self, fai},
    variant::{io::reader::is_stream, log},
};

use crate::raw_input::strerror;

/// The first 16 bytes of every BGZF block, as HTSlib's `bgzf_compression` checks them: the
/// gzip header with no mtime, the `BC` extra field, and its length.
const BGZF_MAGIC: [u8; 16] = [
    0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
];

/// How a FASTA file is compressed (HTSlib's `bgzf_compression`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Compression {
    None,
    Gzip,
    Bgzf,
}

/// A FASTA file and its index could not be loaded. HTSlib's messages have been printed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LoadError;

impl fmt::Display for LoadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("failed to load the FASTA index")
    }
}

impl std::error::Error for LoadError {}

/// An indexed FASTA file open for random access.
pub struct Faidx {
    path: PathBuf,
    index: fai::Index,
    reader: fasta::io::BufReader<File>,
}

impl Faidx {
    /// Opens a FASTA file and its index (`fai_load`). A missing `.fai`, or the `.gzi` of a
    /// BGZF file, is built from the file and saved next to it.
    ///
    /// A stream (a pipe, FIFO, or standard input) is refused before it is opened: it cannot
    /// be read again after indexing, and a FIFO opened again with no writer waits forever.
    ///
    /// # Errors
    ///
    /// Returns [`LoadError`], after printing HTSlib's messages, if the file or its index
    /// cannot be read, built, or saved.
    pub fn load(path: &Path) -> Result<Self, LoadError> {
        let fai_path = with_suffix(path, ".fai");
        let gzi_path = with_suffix(path, ".gzi");
        if is_stream(&path.to_string_lossy()) {
            if fai_path.exists() {
                fail_open_fasta(path);
            } else {
                // What HTSlib reports when it cannot save the index next to a pipe.
                log::error(
                    "fai_build3_core",
                    format_args!(
                        "Failed to open FASTA index {} : No such file or directory",
                        fai_path.display()
                    ),
                );
            }
            return Err(LoadError);
        }
        let index = read_fai(&fai_path)?;
        let file = match File::open(path) {
            Ok(f) => f,
            Err(e) => {
                if index.is_some() {
                    fail_open_fasta(path);
                } else {
                    fail_open_file(path, &e);
                }
                return Err(LoadError);
            }
        };
        // HTSlib peeks at the file while opening it, after reading the index or, without
        // one, while building it.
        let compression = compression(&file).map_err(|e| {
            if index.is_some() {
                fail_open_fasta(path);
            } else {
                fail_open_file(path, &e);
            }
            LoadError
        })?;
        // A BGZF file whose `.gzi` is missing is indexed again, `.fai` and all.
        let saved = match index {
            Some(index) if compression == Compression::None => Some((index, None)),
            Some(index) => read_gzi(&gzi_path, compression)?.map(|gzi| (index, Some(gzi))),
            None => None,
        };
        let (index, gzi) = match saved {
            Some(saved) => saved,
            None => build(path, &file, compression, &fai_path, &gzi_path)?,
        };
        let reader = match gzi {
            Some(gzi) => fasta::io::BufReader::Bgzf(bgzf::io::IndexedReader::new(file, gzi)),
            None => fasta::io::BufReader::Uncompressed(BufReader::new(file)),
        };
        Ok(Self {
            path: path.to_owned(),
            index,
            reader,
        })
    }

    /// The index.
    pub fn index(&self) -> &fai::Index {
        &self.index
    }

    /// Fills `buf` with the (uncompressed) bytes of the file from `offset`.
    ///
    /// # Errors
    ///
    /// Returns an error if the bytes cannot be read, of kind
    /// [`io::ErrorKind::UnexpectedEof`] if the file ends first.
    pub fn read_at(&mut self, offset: u64, buf: &mut [u8]) -> io::Result<()> {
        match &mut self.reader {
            fasta::io::BufReader::Uncompressed(r) => r.get_ref().read_exact_at(buf, offset),
            fasta::io::BufReader::Bgzf(r) => {
                r.seek(SeekFrom::Start(offset))?;
                r.read_exact(buf)
            }
        }
    }

    /// Opens the file again, with the same index, for another thread.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened.
    pub fn reopen(&self) -> io::Result<Self> {
        let file = File::open(&self.path)?;
        let reader = match &self.reader {
            fasta::io::BufReader::Bgzf(r) => {
                fasta::io::BufReader::Bgzf(bgzf::io::IndexedReader::new(file, r.index().clone()))
            }
            fasta::io::BufReader::Uncompressed(_) => {
                fasta::io::BufReader::Uncompressed(BufReader::new(file))
            }
        };
        Ok(Self {
            path: self.path.clone(),
            index: self.index.clone(),
            reader,
        })
    }

    /// The file as an indexed reader of whole records (for CRAM reference repositories).
    pub fn into_indexed_reader(self) -> fasta::io::IndexedReader<fasta::io::BufReader<File>> {
        fasta::io::IndexedReader::new(self.reader, self.index)
    }
}

/// `path` with `suffix` appended to its file name.
fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

fn fail_open_fasta(path: &Path) {
    log::error(
        "fai_load3_core",
        format_args!("Failed to open FASTA file {}", path.display()),
    );
}

fn fail_open_file(path: &Path, e: &io::Error) {
    log::error(
        "fai_build3_core",
        format_args!(
            "Failed to open the file {} : {}",
            path.display(),
            strerror(e)
        ),
    );
}

/// Reads a `.fai` index, or `None` if there is none.
fn read_fai(fai_path: &Path) -> Result<Option<fai::Index>, LoadError> {
    let file = match File::open(fai_path) {
        Ok(f) => f,
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(e) => {
            log::error(
                "fai_load3_core",
                format_args!(
                    "Failed to open FASTA index {}: {}",
                    fai_path.display(),
                    strerror(&e)
                ),
            );
            return Err(LoadError);
        }
    };
    fai::io::Reader::new(BufReader::new(file))
        .read_index()
        .map(Some)
        .map_err(|_| {
            log::error(
                "fai_load3_core",
                format_args!("Failed to read FASTA index {}", fai_path.display()),
            );
            LoadError
        })
}

/// Reads the `.gzi` index of a compressed file with a `.fai` index, or `None` if a BGZF file
/// has none (it is then indexed again).
fn read_gzi(gzi_path: &Path, compression: Compression) -> Result<Option<gzi::Index>, LoadError> {
    match gzi::fs::read(gzi_path) {
        Ok(gzi) => Ok(Some(gzi)),
        Err(e) if e.kind() == io::ErrorKind::NotFound && compression == Compression::Bgzf => {
            Ok(None)
        }
        Err(e) => {
            if e.raw_os_error().is_some() {
                log::error(
                    "bgzf_index_load",
                    format_args!("Error opening {} : {}", gzi_path.display(), strerror(&e)),
                );
            }
            log::error(
                "fai_load3_core",
                format_args!("Failed to load .gzi index: {}", gzi_path.display()),
            );
            Err(LoadError)
        }
    }
}

/// How `file` is compressed, from its first bytes.
fn compression(file: &File) -> io::Result<Compression> {
    let mut head = [0u8; BGZF_MAGIC.len()];
    let mut n = 0;
    while n < head.len() {
        match file.read_at(&mut head[n..], n as u64)? {
            0 => break,
            k => n += k,
        }
    }
    Ok(if n < 2 || head[..2] != BGZF_MAGIC[..2] {
        Compression::None
    } else if n == head.len() && head == BGZF_MAGIC {
        Compression::Bgzf
    } else {
        Compression::Gzip
    })
}

/// Indexes a FASTA file (`fai_build3_core`), saving the `.fai` and, for BGZF, the `.gzi`
/// first. As in HTSlib, an index that cannot be saved is an error.
fn build(
    path: &Path,
    file: &File,
    compression: Compression,
    fai_path: &Path,
    gzi_path: &Path,
) -> Result<(fai::Index, Option<gzi::Index>), LoadError> {
    let index = match compression {
        Compression::Gzip => {
            // HTSlib's BGZF reader stops at the first gzip block it cannot index.
            log::error("fai_build_core", format_args!("File truncated at line 1"));
            log::error(
                "fai_build3_core",
                format_args!("Cannot index files compressed with gzip, please use bgzip"),
            );
            return Err(LoadError);
        }
        Compression::Bgzf => index_records(bgzf::io::Reader::new(BufReader::new(file))),
        Compression::None => index_records(BufReader::new(file)),
    };
    let index = index.map_err(|e| {
        log::error("fai_build_core", format_args!("{e}"));
        LoadError
    })?;
    let gzi = if compression == Compression::Bgzf {
        let gzi = scan_gzi(file).map_err(|e| {
            log::error(
                "fai_build3_core",
                format_args!(
                    "Failed to index the file {} : {}",
                    path.display(),
                    strerror(&e)
                ),
            );
            LoadError
        })?;
        save_gzi(gzi_path, &gzi)?;
        Some(gzi)
    } else {
        None
    };
    save_fai(fai_path, &index)?;
    Ok((index, gzi))
}

/// Saves a `.gzi` index (`bgzf_index_dump`).
fn save_gzi(gzi_path: &Path, gzi: &gzi::Index) -> Result<(), LoadError> {
    let saved = File::create(gzi_path)
        .map_err(|e| ("Error opening", e))
        .and_then(|f| {
            let mut w = gzi::io::Writer::new(BufWriter::new(f));
            w.write_index(gzi)
                .and_then(|()| w.get_mut().flush())
                .map_err(|e| ("Error writing", e))
        });
    saved.map_err(|(what, e)| {
        log::error(
            "bgzf_index_dump",
            format_args!("{what} {} : {}", gzi_path.display(), strerror(&e)),
        );
        log::error(
            "fai_build3_core",
            format_args!("Failed to make bgzf index {}", gzi_path.display()),
        );
        LoadError
    })
}

/// Saves a `.fai` index.
fn save_fai(fai_path: &Path, index: &fai::Index) -> Result<(), LoadError> {
    let saved = File::create(fai_path)
        .map_err(|e| ("Failed to open", e))
        .and_then(|f| {
            let mut w = fai::io::Writer::new(BufWriter::new(f));
            w.write_index(index)
                .and_then(|()| w.get_mut().flush())
                .map_err(|e| ("Failed to write", e))
        });
    saved.map_err(|(what, e)| {
        log::error(
            "fai_build3_core",
            format_args!(
                "{what} FASTA index {} : {}",
                fai_path.display(),
                strerror(&e)
            ),
        );
        LoadError
    })
}

fn index_records(reader: impl io::BufRead) -> io::Result<fai::Index> {
    let mut indexer = fasta::io::Indexer::new(reader);
    let mut records = Vec::new();
    while let Some(r) = indexer.index_record()? {
        records.push(r);
    }
    Ok(fai::Index::from(records))
}

/// Indexes the BGZF blocks of a file as `bgzip -r` and HTSlib's `fai_build` do: the
/// compressed and uncompressed offsets of every block that holds data, except the first.
fn scan_gzi(file: &File) -> io::Result<gzi::Index> {
    let mut entries = Vec::new();
    let (mut coff, mut uoff) = (0u64, 0u64);
    let mut header = [0u8; 18];
    let mut footer = [0u8; 4];
    loop {
        match file.read_exact_at(&mut header, coff) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e),
        }
        // BSIZE, the block size minus 1, is in the BC extra field of a BGZF header.
        if header[..4] != BGZF_MAGIC[..4] || header[12..14] != *b"BC" {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not BGZF"));
        }
        let bsize = u64::from(u16::from_le_bytes([header[16], header[17]])) + 1;
        file.read_exact_at(&mut footer, coff + bsize - 4)?;
        let isize = u64::from(u32::from_le_bytes(footer));
        if isize > 0 {
            entries.push((coff, uoff));
        }
        coff += bsize;
        uoff += isize;
    }
    // The first block is implicit.
    if !entries.is_empty() {
        entries.remove(0);
    }
    Ok(entries.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raw_input::make_fifo;

    const FASTA: &[u8] = b">a\nACGTACGTACG\n>b\nTTTT\n";

    fn bgzip(data: &[u8], block: usize) -> Vec<u8> {
        let mut w = bgzf::io::Writer::new(Vec::new());
        for chunk in data.chunks(block) {
            w.write_all(chunk).expect("writing to memory");
            w.flush().expect("writing to memory");
        }
        w.finish().expect("writing to memory")
    }

    /// A new, empty directory for one test, which the test removes.
    fn test_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("vtools-faidx-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("creating the test directory");
        dir
    }

    fn fetch_all(f: &mut Faidx) -> Vec<Vec<u8>> {
        let records: Vec<fai::Record> = f.index().as_ref().to_vec();
        records
            .iter()
            .map(|r| {
                let mut buf = vec![0; r.length() as usize];
                f.read_at(r.position(), &mut buf)
                    .expect("reading the bases");
                buf
            })
            .collect()
    }

    #[test]
    fn loads_plain_and_bgzf_files_building_their_indexes() {
        let dir = test_dir("load");
        let plain = dir.join("p.fa");
        std::fs::write(&plain, FASTA).expect("writing");
        let gz = dir.join("g.fa.gz");
        std::fs::write(&gz, bgzip(FASTA, 7)).expect("writing");
        for path in [&plain, &gz] {
            let mut f = Faidx::load(path).expect("loading");
            assert_eq!(
                fetch_all(&mut f),
                [&b"ACGTACGTACG"[..], b"TTTT"],
                "{path:?}"
            );
            assert!(with_suffix(path, ".fai").exists(), "{path:?}");
            // Loading again reads the saved indexes.
            let mut f = Faidx::load(path).expect("loading again");
            assert_eq!(fetch_all(&mut f)[1], b"TTTT", "{path:?}");
        }
        assert!(with_suffix(&gz, ".gzi").exists());
        assert!(!with_suffix(&plain, ".gzi").exists());
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }

    #[test]
    fn rebuilds_a_missing_gzi() {
        let dir = test_dir("gzi");
        let gz = dir.join("g.fa.gz");
        std::fs::write(&gz, bgzip(FASTA, 5)).expect("writing");
        drop(Faidx::load(&gz).expect("loading"));
        let gzi_path = with_suffix(&gz, ".gzi");
        let gzi = std::fs::read(&gzi_path).expect("reading .gzi");
        std::fs::remove_file(&gzi_path).expect("removing .gzi");
        let mut f = Faidx::load(&gz).expect("loading without .gzi");
        assert_eq!(fetch_all(&mut f)[0], b"ACGTACGTACG");
        assert_eq!(std::fs::read(&gzi_path).expect("reading .gzi"), gzi);
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }

    #[test]
    fn the_gzi_has_htslibs_entries() {
        let dir = test_dir("scan");
        let gz = dir.join("g.fa.gz");
        // Five blocks of data (the last three bytes alone) and the end-of-file block.
        std::fs::write(&gz, bgzip(FASTA, 5)).expect("writing");
        let file = File::open(&gz).expect("opening");
        let gzi = scan_gzi(&file).expect("scanning");
        let uoffs: Vec<u64> = gzi.as_ref().iter().map(|&(_, u)| u).collect();
        assert_eq!(uoffs, [5, 10, 15, 20]);
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }

    #[test]
    fn an_index_that_cannot_be_saved_is_an_error() {
        use std::os::unix::fs::PermissionsExt;

        let dir = test_dir("unsaved");
        let plain = dir.join("p.fa");
        std::fs::write(&plain, FASTA).expect("writing");
        let gz = dir.join("g.fa.gz");
        std::fs::write(&gz, bgzip(FASTA, 7)).expect("writing");
        let set_mode = |mode| {
            std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(mode))
                .expect("changing the permissions of the test directory");
        };
        set_mode(0o555);
        // A user who may write anyway (root) has nothing to test.
        if File::create(dir.join("probe")).is_err() {
            for path in [&plain, &gz] {
                assert!(Faidx::load(path).is_err(), "{path:?}");
            }
            set_mode(0o755);
            // With its indexes saved, a file in a read-only directory loads.
            for path in [&plain, &gz] {
                drop(Faidx::load(path).expect("loading"));
            }
            set_mode(0o555);
            for path in [&plain, &gz] {
                let mut f = Faidx::load(path).expect("loading with saved indexes");
                assert_eq!(fetch_all(&mut f)[1], b"TTTT", "{path:?}");
            }
        }
        set_mode(0o755);
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }

    #[test]
    fn refuses_gzip_and_streams() {
        let dir = test_dir("refuse");
        let gz = dir.join("g.fa.gz");
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(FASTA).expect("compressing");
        std::fs::write(&gz, enc.finish().expect("compressing")).expect("writing");
        assert!(Faidx::load(&gz).is_err());
        assert!(!with_suffix(&gz, ".fai").exists());

        assert!(Faidx::load(Path::new("-")).is_err());
        std::fs::remove_dir_all(&dir).expect("removing the test directory");

        let (fifo_dir, fifo) = make_fifo("faidx");
        // Nothing writes to the FIFO: opening it would wait forever.
        let fifo_path = fifo.clone();
        let loaded = crate::raw_input::within_seconds("loading a FIFO with no writer", move || {
            Faidx::load(&fifo_path).is_ok()
        });
        assert!(!loaded);
        assert!(!with_suffix(&fifo, ".fai").exists());
        std::fs::remove_dir_all(&fifo_dir).expect("removing the test directory");
    }
}
