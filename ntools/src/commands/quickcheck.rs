//! `ntools quickcheck`: check that alignment files look intact (samtools quickcheck).

use std::{
    fs::File,
    io::{self, BufRead},
    os::{
        fd::AsFd,
        unix::fs::{FileExt, FileTypeExt},
    },
    path::{Path, PathBuf},
    process::ExitCode,
};

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    bam_blocks::BGZF_EOF,
    cram_file::{EOF_2, EOF_3, FileStart},
    source::{Compression, Format, Input, Stream},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Files to check.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Print the names of files that fail to stdout; repeat for progress on stderr.
    #[arg(short = 'v', action = clap::ArgAction::Count)]
    verbose: u8,

    /// Do not print messages about failures to stderr (unless -vv).
    #[arg(short = 'q')]
    quiet: bool,

    /// Do not require @SQ lines (for unmapped files).
    #[arg(short = 'u')]
    unmapped: bool,
}

/// The ways a file can fail, as samtools reports them: the exit status is the bitwise OR of the
/// failures of all files.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Failure {
    Open,
    NotSequence,
    BadHeader,
    NoEofBlock,
    Close,
}

impl Failure {
    fn bit(self) -> u8 {
        match self {
            Self::Open => 2,
            Self::NotSequence => 4,
            Self::BadHeader => 8,
            Self::NoEofBlock => 16,
            Self::Close => 32,
        }
    }
}

/// Reports failures and progress for one file.
struct Report<'a> {
    args: &'a Args,
    name: String,
    state: u8,
}

impl Report<'_> {
    /// Records a failure, with samtools' message unless `-q` (and not `-vv`).
    fn fail(&mut self, failure: Failure, message: &str) {
        self.state |= failure.bit();

        if !self.args.quiet || self.args.verbose >= 2 {
            eprintln!("{} {message}", self.name);
        }
    }

    /// Prints a progress message with `-vvv`.
    fn progress(&self, message: &str) {
        if self.args.verbose >= 3 {
            eprintln!("{} {message}", self.name);
        }
    }
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    if args.verbose >= 2 {
        eprintln!("verbosity set to {}", args.verbose);
    }

    let mut status = 0;

    for path in &args.inputs {
        let mut report = Report {
            args: &args,
            name: path.display().to_string(),
            state: 0,
        };

        if args.verbose >= 3 {
            eprintln!("checking {}", report.name);
        }

        check(path, &mut report);

        if report.state != 0 && args.verbose >= 1 {
            println!("{}", report.name);
        }

        status |= report.state;
    }

    Ok(ExitCode::from(status))
}

/// `main_quickcheck`'s checks of one file: that it opens, is alignment data, has a header with
/// reference sequences (unless `-u`), and ends with an EOF block where its format has one.
fn check(path: &Path, report: &mut Report<'_>) {
    let Some(src) = path.to_str() else {
        report.fail(Failure::Open, "could not be opened for reading.");
        return;
    };

    let Ok(mut input) = Input::open(src) else {
        report.fail(Failure::Open, "could not be opened for reading.");
        return;
    };

    // The input is read once, in order, so it can be a pipe. The EOF block is checked through
    // another handle to the same file, where it can seek.
    let Ok(eof_check) = EofCheck::new(&mut input) else {
        report.fail(Failure::Open, "could not be opened for reading.");
        return;
    };

    // HTSlib reads a CRAM file's header when it opens the file.
    let cram_header = match input.format {
        Format::Cram => match FileStart::read(&mut input.reader).and_then(|start| start.header()) {
            Ok(header) => Some(header),
            Err(_) => {
                report.fail(Failure::Open, "could not be opened for reading.");
                return;
            }
        },
        _ => None,
    };

    if report.args.verbose >= 3 {
        eprintln!("opened {}", report.name);
    }

    if eof_check.is_empty {
        // HTSlib does not recognize an empty file as any format.
        report.fail(Failure::NotSequence, "was not identified as sequence data.");
        report.progress(NO_EOF_BLOCK_IN_FORMAT);
        return;
    }

    report.progress("is sequence data");

    let compression = input.compression;

    let header = match cram_header {
        Some(header) => Ok(header),
        None => AlignmentReader::from_input(input, &ReadOptions::default())
            .and_then(|mut reader| reader.read_header()),
    };

    match header.map(|header| header.reference_sequences().len()) {
        Ok(0) if !report.args.unmapped => {
            report.fail(Failure::BadHeader, "had no targets in header.");
        }
        Ok(n) => report.progress(&format!("has {n} targets in header.")),
        Err(e) => {
            report.fail(
                Failure::BadHeader,
                "caused an error whilst reading its header.",
            );

            // HTSlib's BGZF reader keeps the error of a truncated block, so closing the file
            // fails too.
            if compression == Compression::Bgzf && is_truncation(&e) {
                report.state |= Failure::Close.bit();
            }
        }
    }

    match eof_check.state() {
        Ok(EofState::Present) => report.progress("has good EOF block."),
        Ok(EofState::Missing) => report.fail(
            Failure::NoEofBlock,
            "was missing EOF block when one should be present.",
        ),
        Ok(EofState::NotApplicable) => report.progress(NO_EOF_BLOCK_IN_FORMAT),
        Ok(EofState::NotSeekable) => {
            report.progress("cannot be checked for EOF block as it is not seekable.");
        }
        Err(_) => report.fail(
            Failure::NoEofBlock,
            "caused an error whilst checking for EOF block.",
        ),
    }
}

const NO_EOF_BLOCK_IN_FORMAT: &str =
    "cannot be checked for EOF block because its filetype does not contain one.";

/// Whether an error comes from data that ends early.
fn is_truncation(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|cause| cause.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::UnexpectedEof)
}

#[derive(Debug, Eq, PartialEq)]
enum EofState {
    Present,
    Missing,
    /// The format has no EOF marker (uncompressed SAM, CRAM before 2.1).
    NotApplicable,
    /// The input is a pipe, which cannot be read from its end.
    NotSeekable,
}

/// What `hts_check_EOF` needs, taken from an input before its header is read.
struct EofCheck {
    format: Format,
    compression: Compression,
    /// The CRAM file definition's major and minor version.
    cram_version: Option<(u8, u8)>,
    /// Whether the input holds no data at all.
    is_empty: bool,
    /// Another handle to the input's file, if it can seek.
    file: Option<File>,
}

impl EofCheck {
    fn new(input: &mut Input) -> io::Result<Self> {
        let start = input.reader.fill_buf()?;
        let is_empty = start.is_empty();
        let cram_version = match start.get(4..6) {
            Some(&[major, minor]) if input.format == Format::Cram => Some((major, minor)),
            _ => None,
        };

        let file = match input.reader.get_ref() {
            Stream::File(file, _) => Some(file.try_clone()?),
            Stream::Pipe(_) => None,
            // stdin can be a file (`- < in.bam`), which HTSlib reads as one.
            Stream::Stdin(stdin) => {
                let file = File::from(stdin.get_ref().as_fd().try_clone_to_owned()?);
                let file_type = file.metadata()?.file_type();
                (!file_type.is_fifo() && !file_type.is_socket()).then_some(file)
            }
        };

        Ok(Self {
            format: input.format,
            compression: input.compression,
            cram_version,
            is_empty,
            file,
        })
    }

    /// `hts_check_EOF`
    fn state(&self) -> io::Result<EofState> {
        let eof: Option<&[u8]> = match (self.format, self.compression) {
            (Format::Cram, _) => match self.cram_version {
                Some((3.., _)) => Some(&EOF_3),
                Some((2, 1..)) => Some(&EOF_2),
                Some(_) => return Ok(EofState::NotApplicable),
                None => return Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
            },
            (_, Compression::Bgzf) => None,
            (_, Compression::None | Compression::Gzip) => return Ok(EofState::NotApplicable),
        };

        let Some(file) = &self.file else {
            return Ok(EofState::NotSeekable);
        };

        let present = match eof {
            Some(eof) => ends_with_eof_container(file, eof)?,
            None => ends_with(file, &BGZF_EOF)?.is_some_and(|buf| buf == BGZF_EOF),
        };

        Ok(if present {
            EofState::Present
        } else {
            EofState::Missing
        })
    }
}

/// Reads the last `tail.len()` bytes of a file without moving its position, or returns `None` if
/// the file is shorter.
fn ends_with(file: &File, tail: &[u8]) -> io::Result<Option<Vec<u8>>> {
    let len = file.metadata()?.len();

    let Some(offset) = len.checked_sub(tail.len() as u64) else {
        return Ok(None);
    };

    let mut buf = vec![0; tail.len()];
    file.read_exact_at(&mut buf, offset)?;

    Ok(Some(buf))
}

/// `cram_check_EOF`: whether the file ends with the EOF container `eof`.
fn ends_with_eof_container(file: &File, eof: &[u8]) -> io::Result<bool> {
    Ok(ends_with(file, eof)?.is_some_and(|buf| matches_eof_container(&buf, eof)))
}

/// Like HTSlib, compares only the low 4 bits of byte 8 (the last ITF-8 byte of the reference
/// sequence ID), which early Java and C implementations wrote differently.
fn matches_eof_container(buf: &[u8], eof: &[u8]) -> bool {
    const MASKED: usize = 8;

    buf.len() == eof.len()
        && buf.iter().zip(eof).enumerate().all(|(i, (&a, &b))| {
            if i == MASKED {
                a & 0x0f == b & 0x0f
            } else {
                a == b
            }
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_failure_bits() {
        let all = [
            Failure::Open,
            Failure::NotSequence,
            Failure::BadHeader,
            Failure::NoEofBlock,
            Failure::Close,
        ]
        .iter()
        .fold(0, |bits, failure| bits | failure.bit());

        assert_eq!(all, 2 | 4 | 8 | 16 | 32);
    }

    #[test]
    fn test_matches_eof_container() {
        // HTSlib writes the CRAM 2.1 EOF container's reference sequence ID as ff ff ff ff 0f.
        let mut htslib_eof = EOF_2;
        htslib_eof[8] = 0x0f;
        assert!(matches_eof_container(&htslib_eof, &EOF_2));
        assert!(matches_eof_container(&EOF_3, &EOF_3));

        let mut other = EOF_3;
        other[20] ^= 1;
        assert!(!matches_eof_container(&other, &EOF_3));
        assert!(!matches_eof_container(&EOF_3[1..], &EOF_3));
    }
}
