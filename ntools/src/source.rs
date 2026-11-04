//! Opening alignment inputs from paths or stdin and detecting their format.

use std::{
    fmt,
    fs::File,
    io::{self, BufReader, IsTerminal, Read, Seek, SeekFrom},
    path::PathBuf,
};

use anyhow::{Context, bail};
use flate2::bufread::MultiGzDecoder;
use pasta::util::{
    inputs::Opened,
    io::{Rewound, read_detection_prefix},
};

/// The input buffer size.
const BUFFER_CAPACITY: usize = 1 << 17;

/// The smallest input buffer, for a regular file shorter than [`BUFFER_CAPACITY`].
const MIN_BUFFER_CAPACITY: usize = 1 << 13;

/// The path that names stdin, as in samtools.
pub const STDIN: &str = "-";

/// The input path, or stdin when none is given, as samtools reads stdin when no input is named
/// and stdin is not a terminal.
pub fn input_or_stdin(input: Option<String>) -> anyhow::Result<String> {
    match input {
        Some(input) => Ok(input),
        None if io::stdin().is_terminal() => {
            bail!("no input file given, and stdin is a terminal (use `-` to read it)")
        }
        None => Ok(STDIN.to_owned()),
    }
}

/// An alignment file format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Sam,
    Bam,
    Cram,
}

/// How a stream is compressed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    None,
    Bgzf,
    /// gzip that is not BGZF (say, from `gzip`), which HTSlib also reads (for SAM).
    Gzip,
}

/// A raw input stream: a file, or a pipe (stdin, or a file that cannot seek, such as a FIFO)
/// with the bytes read for format detection put back in front.
pub enum Stream {
    /// A regular file, recorded as an input while it is open.
    File(
        File,
        #[expect(dead_code, reason = "held for its drop, which forgets the input")] Option<Opened>,
    ),
    Pipe(Rewound<File>),
    Stdin(Rewound<io::Stdin>),
}

impl Read for Stream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::File(f, _) => f.read(buf),
            Self::Pipe(p) => p.read(buf),
            Self::Stdin(s) => s.read(buf),
        }
    }
}

impl Seek for Stream {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match self {
            Self::File(f, _) => f.seek(pos),
            Self::Pipe(_) => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "the input is not seekable",
            )),
            Self::Stdin(_) => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "stdin is not seekable",
            )),
        }
    }
}

/// An opened input with its detected format.
pub struct Input {
    pub path: Option<PathBuf>,
    pub reader: BufReader<Stream>,
    pub format: Format,
    pub compression: Compression,
    /// The format the input looks like, when it is neither SAM, BAM nor CRAM (detected SAM that
    /// is really VCF, FASTA, or other text).
    pub other_format: Option<OtherFormat>,
}

/// A format that an input detected as SAM looks like instead, named in the error that rejects
/// it as an alignment input.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OtherFormat {
    Vcf,
    Bcf,
    Fasta,
    Fastq,
    /// Text whose first line is neither a SAM header line nor a SAM record.
    Text,
}

impl fmt::Display for OtherFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind = match self {
            Self::Vcf => "it looks like VCF",
            Self::Bcf => "it looks like BCF",
            Self::Fasta => "it looks like FASTA",
            Self::Fastq => "it looks like FASTQ",
            Self::Text => "its first line is neither a SAM header line nor a record",
        };

        write!(f, "not SAM, BAM or CRAM ({kind})")
    }
}

impl Input {
    /// Opens `src` (`-` for stdin) and detects its format from its leading bytes.
    pub fn open(src: &str) -> anyhow::Result<Self> {
        let detecting = || format!("detecting the format of {src}");

        let (path, prefix, stream, capacity) = if src == STDIN {
            let mut stdin = io::stdin();
            let mut prefix = read_detection_prefix(&mut stdin).with_context(detecting)?;
            extend_text_prefix(&mut stdin, &mut prefix).with_context(detecting)?;
            (
                None,
                prefix.clone(),
                Stream::Stdin(Rewound::new(prefix, stdin)),
                BUFFER_CAPACITY,
            )
        } else {
            let mut file = File::open(src).with_context(|| format!("opening {src}"))?;
            let mut prefix = read_detection_prefix(&mut file).with_context(detecting)?;
            extend_text_prefix(&mut file, &mut prefix).with_context(detecting)?;

            // A file that cannot seek back (a FIFO, `<(...)`) is read like a pipe.
            let (stream, capacity) = match file.seek(SeekFrom::Start(0)) {
                Ok(_) => {
                    // `BufReader` zeroes its whole buffer on the first read (`Stream` has no
                    // `read_buf`), so a buffer longer than the file costs memory for nothing:
                    // 128 KiB for each of the inputs of `depth` or `coverage`.
                    let len = file.metadata().map_or(u64::MAX, |m| m.len());
                    let capacity = usize::try_from(len)
                        .unwrap_or(usize::MAX)
                        .clamp(MIN_BUFFER_CAPACITY, BUFFER_CAPACITY);
                    let opened = Opened::new(&file);
                    (Stream::File(file, opened), capacity)
                }
                Err(_) => (
                    Stream::Pipe(Rewound::new(prefix.clone(), file)),
                    BUFFER_CAPACITY,
                ),
            };

            (Some(PathBuf::from(src)), prefix, stream, capacity)
        };

        let (format, compression) = detect(&prefix).with_context(detecting)?;
        let other_format = match (format, compression) {
            (Format::Sam, Compression::None) => sniff_other_format(&prefix),
            (Format::Sam, _) => sniff_other_format(&decompressed_prefix(&prefix)),
            _ => None,
        };
        let reader = BufReader::with_capacity(capacity, stream);

        Ok(Self {
            path,
            reader,
            format,
            compression,
            other_format,
        })
    }
}

/// Detects the format and compression from the start of a stream, as read by
/// [`read_detection_prefix`].
pub fn detect(src: &[u8]) -> io::Result<(Format, Compression)> {
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];
    const BAM_MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";
    const CRAM_MAGIC_NUMBER: [u8; 4] = *b"CRAM";

    if src.starts_with(&GZIP_MAGIC_NUMBER) {
        let mut decoder = MultiGzDecoder::new(src);
        let mut buf = [0; BAM_MAGIC_NUMBER.len()];

        let format = match decoder.read_exact(&mut buf) {
            Ok(()) if buf == BAM_MAGIC_NUMBER => Format::Bam,
            _ => Format::Sam,
        };

        let compression = if is_bgzf(src) {
            Compression::Bgzf
        } else {
            Compression::Gzip
        };

        Ok((format, compression))
    } else if src.starts_with(&CRAM_MAGIC_NUMBER) {
        Ok((Format::Cram, Compression::None))
    } else if src.starts_with(&BAM_MAGIC_NUMBER) {
        bail_io("uncompressed (raw) BAM streams are not supported")
    } else {
        Ok((Format::Sam, Compression::None))
    }
}

/// How much of an uncompressed text input is read to tell SAM from other text formats.
const TEXT_PREFIX_LEN: usize = 1 << 12;

/// The number of lines [`sniff_other_format`] looks at.
const TEXT_PREFIX_LINES: usize = 3;

/// Reads more of an uncompressed text input into `prefix`, up to [`TEXT_PREFIX_LINES`] lines or
/// [`TEXT_PREFIX_LEN`] bytes, so that [`sniff_other_format`] sees whole lines. Binary and
/// compressed inputs are left as they are.
fn extend_text_prefix<R: Read>(reader: &mut R, prefix: &mut Vec<u8>) -> io::Result<()> {
    if prefix.starts_with(&[0x1f, 0x8b])
        || prefix.starts_with(b"CRAM")
        || prefix.starts_with(b"BAM\x01")
    {
        return Ok(());
    }

    let mut buf = [0; TEXT_PREFIX_LEN];

    while prefix.len() < TEXT_PREFIX_LEN && line_count(prefix) < TEXT_PREFIX_LINES {
        let n = match reader.read(&mut buf[..TEXT_PREFIX_LEN - prefix.len()]) {
            Ok(0) => break,
            Ok(n) => n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        };
        prefix.extend_from_slice(&buf[..n]);
    }

    Ok(())
}

fn line_count(src: &[u8]) -> usize {
    src.iter().filter(|&&b| b == b'\n').count()
}

/// Decompresses as much of a gzip prefix as it holds, up to [`TEXT_PREFIX_LEN`] bytes.
fn decompressed_prefix(src: &[u8]) -> Vec<u8> {
    let mut buf = Vec::new();
    // The prefix may end inside a gzip member, so the error at its end is expected; the bytes
    // decompressed before it are kept.
    let _ = MultiGzDecoder::new(src)
        .take(TEXT_PREFIX_LEN as u64)
        .read_to_end(&mut buf);
    buf
}

/// Returns the format the (decompressed) start of an input detected as SAM looks like, if it is
/// clearly not SAM. Only complete lines are judged, and a header-less SAM file, whose records
/// have tab-separated fields, is never rejected.
pub fn sniff_other_format(src: &[u8]) -> Option<OtherFormat> {
    if src.starts_with(b"##fileformat=VCF") {
        return Some(OtherFormat::Vcf);
    } else if src.starts_with(b"BCF\x02") {
        return Some(OtherFormat::Bcf);
    }

    let mut lines = src.split(|&b| b == b'\n');
    let first = lines.next().filter(|_| src.contains(&b'\n'))?;
    let first = first.strip_suffix(b"\r").unwrap_or(first);

    if first.is_empty() || first.contains(&b'\t') {
        return None;
    }

    match first[0] {
        b'>' => Some(OtherFormat::Fasta),
        // A FASTQ record's third line starts with `+`, which no SAM line does.
        b'@' if lines.nth(1).is_some_and(|third| third.starts_with(b"+")) => {
            Some(OtherFormat::Fastq)
        }
        // A SAM header line, say `@CO` with no text.
        b'@' => None,
        _ => Some(OtherFormat::Text),
    }
}

/// Whether a gzip stream starts with a BGZF block header (the `BC` extra subfield).
fn is_bgzf(src: &[u8]) -> bool {
    const FEXTRA: u8 = 0x04;
    src.get(3).is_some_and(|flags| flags & FEXTRA != 0) && src.get(12..14) == Some(b"BC")
}

fn bail_io<T>(msg: &str) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, msg.to_owned()))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use pasta::bgzf;

    use super::*;

    fn bgzf_of(data: &[u8]) -> io::Result<Vec<u8>> {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(data)?;
        writer.finish()
    }

    #[test]
    fn test_detect() -> io::Result<()> {
        let cases: [(Vec<u8>, Format, Compression); 5] = [
            (b"@HD\tVN:1.6\n".to_vec(), Format::Sam, Compression::None),
            (Vec::new(), Format::Sam, Compression::None),
            (b"CRAM\x03\x01".to_vec(), Format::Cram, Compression::None),
            (bgzf_of(b"BAM\x01\0\0\0\0")?, Format::Bam, Compression::Bgzf),
            (bgzf_of(b"@HD\tVN:1.6\n")?, Format::Sam, Compression::Bgzf),
        ];

        for (src, format, compression) in cases {
            let got = detect(&src)?;
            assert_eq!(got, (format, compression), "input {src:?}");
        }

        Ok(())
    }
}
