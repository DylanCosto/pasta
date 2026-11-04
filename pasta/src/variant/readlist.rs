//! Lists given on the command line or in files (HTSlib's `hts_readlist` and `hts_readlines`).

use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

use crate::variant::io::reader::EofWatched;

/// Where a list comes from.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ListSource {
    /// The argument itself is a comma-separated list.
    Inline,
    /// The argument names a file (plain or gzip-compressed) with one item per line.
    File,
}

/// Reads a list as HTSlib's `hts_readlist` does: a comma-separated string (empty items kept),
/// or the non-empty lines of a file (with a trailing `\r` removed).
///
/// # Errors
///
/// Returns an error when the file cannot be read.
pub fn read_list(s: &[u8], source: ListSource) -> io::Result<Vec<Vec<u8>>> {
    match source {
        ListSource::Inline => Ok(s.split(|&c| c == b',').map(<[u8]>::to_vec).collect()),
        ListSource::File => read_lines(Path::new(&*String::from_utf8_lossy(s))),
    }
}

/// Reads the non-empty lines of a plain or gzip-compressed file, without line endings
/// (HTSlib's `hts_readlines` on a file; `-` reads standard input).
///
/// # Errors
///
/// Returns an error when the file cannot be read.
pub fn read_lines(path: &Path) -> io::Result<Vec<Vec<u8>>> {
    Ok(read_all_lines(path)?
        .into_iter()
        .filter(|l| !l.is_empty())
        .collect())
}

/// Reads a list as HTSlib's `hts_readlines` does: the non-empty lines of a file, or, when no
/// such file opens and the argument starts with `:`, the comma-separated items after the colon
/// (empty items kept). bcftools reads its `-l`/`--file-list` options this way.
///
/// # Errors
///
/// Returns the error of opening or reading the file when the argument is not a `:` list.
pub fn read_lines_or_inline(s: &[u8]) -> io::Result<Vec<Vec<u8>>> {
    read_lines(Path::new(&*String::from_utf8_lossy(s))).or_else(|e| match s.strip_prefix(b":") {
        Some(list) => read_list(list, ListSource::Inline),
        None => Err(e),
    })
}

/// Opens a plain or gzip-compressed text file for reading (`-` reads standard input), as
/// HTSlib's `hts_open(fname, "r")` does for text, to read large files line by line. A BGZF file
/// read to its end without the EOF marker block is reported ([`report_missing_eof`]).
///
/// [`report_missing_eof`]: crate::variant::io::reader::report_missing_eof
///
/// # Errors
///
/// Returns an error when the file cannot be opened or read.
pub fn open_text(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let reader: Box<dyn Read> = if path == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path)?)
    };
    let mut reader = BufReader::new(reader);
    let head = reader.fill_buf()?;
    let is_gzip = head.starts_with(&[0x1f, 0x8b]);
    let is_bgzf = is_gzip && head.len() >= 16 && head[3] & 4 != 0 && &head[12..16] == b"BC\x02\x00";
    Ok(if is_bgzf {
        let watched = EofWatched::new(reader, &path.to_string_lossy());
        Box::new(BufReader::new(flate2::bufread::MultiGzDecoder::new(
            BufReader::new(watched),
        )))
    } else if is_gzip {
        Box::new(BufReader::new(flate2::bufread::MultiGzDecoder::new(reader)))
    } else {
        Box::new(reader)
    })
}

/// Reads all lines of a plain or gzip-compressed file, empty ones included, without line
/// endings (a trailing `\r` is removed) or anything after a NUL byte.
///
/// # Errors
///
/// Returns an error when the file cannot be read.
pub fn read_all_lines(path: &Path) -> io::Result<Vec<Vec<u8>>> {
    read_all_lines_from(&mut open_text(path)?)
}

/// Reads all lines of an input opened with [`open_text`] as [`read_all_lines`] reads them, for
/// callers that tell a file that cannot be opened from one that cannot be read.
///
/// # Errors
///
/// Returns an error when the input cannot be read (a truncated or corrupt gzip member, say).
pub fn read_all_lines_from(reader: &mut impl BufRead) -> io::Result<Vec<Vec<u8>>> {
    let mut out = Vec::new();
    let mut line = Vec::new();
    while read_line(reader, &mut line)? {
        out.push(line.clone());
    }
    Ok(out)
}

/// Reads the next line into `line` as [`read_all_lines`] reads each: without its line ending
/// (a trailing `\r` is removed) or anything after a NUL byte. Returns `false` at the end.
///
/// # Errors
///
/// Returns an error when the input cannot be read.
pub fn read_line(reader: &mut impl BufRead, line: &mut Vec<u8>) -> io::Result<bool> {
    line.clear();
    if reader.read_until(b'\n', line)? == 0 {
        return Ok(false);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
    // HTSlib copies each line with strdup, so it ends at a NUL byte.
    if let Some(nul) = memchr::memchr(0, line) {
        line.truncate(nul);
    }
    Ok(true)
}
