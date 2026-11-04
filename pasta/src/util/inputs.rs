//! The regular files this process has open as inputs, so that an output naming one of them is
//! refused before it is truncated (as `-o in.vcf.gz in.vcf.gz` would otherwise destroy the input
//! and leave a file with a header and no records).

use std::{
    fs::{File, Metadata},
    io,
    path::Path,
    sync::{Mutex, PoisonError},
};

/// The device and inode of every regular file open as an input, once per [`Opened`].
static OPEN: Mutex<Vec<(u64, u64)>> = Mutex::new(Vec::new());

/// The device and inode of a regular file; `None` for anything else (a pipe, a device).
fn identity(meta: &Metadata) -> Option<(u64, u64)> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        meta.is_file().then(|| (meta.dev(), meta.ino()))
    }
    #[cfg(not(unix))]
    {
        let _ = meta;
        None
    }
}

/// The identity of the regular file at `path`, following symbolic links.
fn path_identity(path: &Path) -> Option<(u64, u64)> {
    std::fs::metadata(path).ok().as_ref().and_then(identity)
}

fn open_inputs() -> std::sync::MutexGuard<'static, Vec<(u64, u64)>> {
    OPEN.lock().unwrap_or_else(PoisonError::into_inner)
}

/// A regular file open as an input, recorded for [`check_output`] until this is dropped. Keep
/// it as long as the file is open.
#[derive(Debug)]
#[must_use = "the input is forgotten when this is dropped"]
pub struct Opened((u64, u64));

impl Opened {
    /// Records `file`, just opened as an input; `None` when it is not a regular file.
    pub fn new(file: &File) -> Option<Self> {
        let id = file.metadata().ok().as_ref().and_then(identity)?;
        open_inputs().push(id);
        Some(Self(id))
    }
}

impl Drop for Opened {
    fn drop(&mut self) {
        let mut open = open_inputs();
        if let Some(i) = open.iter().position(|&id| id == self.0) {
            open.swap_remove(i);
        }
    }
}

fn refused() -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        "the output file is also an input",
    )
}

/// Checks that the output `path` is not a regular file open as an input ([`Opened`]). A path
/// that does not exist yet, or is not a regular file, passes.
///
/// # Errors
///
/// Returns an [`io::ErrorKind::InvalidInput`] error when `path` is one of the inputs.
pub fn check_output(path: &Path) -> io::Result<()> {
    match path_identity(path) {
        Some(id) if open_inputs().contains(&id) => Err(refused()),
        _ => Ok(()),
    }
}

/// [`check_output`] against the named `inputs`, for a command that closes its inputs before it
/// opens the output and reads them again afterwards.
///
/// # Errors
///
/// As [`check_output`].
pub fn check_output_against<P: AsRef<Path>>(path: &Path, inputs: &[P]) -> io::Result<()> {
    let Some(id) = path_identity(path) else {
        return Ok(());
    };
    if inputs.iter().any(|p| path_identity(p.as_ref()) == Some(id)) {
        return Err(refused());
    }
    Ok(())
}
