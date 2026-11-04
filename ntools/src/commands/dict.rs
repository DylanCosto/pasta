//! `ntools dict`: write a SAM sequence dictionary for a FASTA file (samtools dict).

use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use flate2::bufread::MultiGzDecoder;
use md5::{Digest, Md5};

use crate::source::{self, Compression, Input};

/// The decompressed input's buffer size.
const BUFFER_CAPACITY: usize = 1 << 17;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// FASTA file (`-` for stdin) [default: stdin, unless it is a terminal].
    input: Option<String>,

    /// Assembly identifier (AS).
    #[arg(short = 'a', long = "assembly")]
    assembly: Option<String>,

    /// Species (SP).
    #[arg(short = 's', long = "species")]
    species: Option<String>,

    /// URI (UR) [default: file:// and the input's absolute path].
    #[arg(short = 'u', long = "uri")]
    uri: Option<String>,

    /// Add an AN tag naming the sequence with `chr` added or removed.
    #[arg(short = 'A', long = "alias", visible_alias = "alternative-name")]
    alias: bool,

    /// Add AH:* to the sequences named in this alternate locus (`.alt`) file.
    #[arg(short = 'l', long = "alt")]
    alt: Option<PathBuf>,

    /// Do not print the @HD line.
    #[arg(short = 'H', long = "no-header")]
    no_header: bool,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let alternates = args.alt.as_deref().map(read_alt_names).transpose()?;

    let input = source::input_or_stdin(args.input.clone())?;
    let src = input.as_str();

    let uri = match &args.uri {
        Some(uri) => Some(uri.clone()),
        // As samtools, stdin has no URI.
        None if src == source::STDIN => None,
        None => {
            // As samtools' realpath(), symbolic links are resolved.
            let path = std::fs::canonicalize(src).with_context(|| format!("resolving {src}"))?;
            Some(format!("file://{}", path.display()))
        }
    };

    let mut reader = open_fasta(src)?;

    let sink: Box<dyn Write> = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(sink);

    if !args.no_header {
        writeln!(out, "@HD\tVN:1.0\tSO:unsorted")?;
    }

    let mut current: Option<Sequence> = None;
    let mut seen_sequence = false;
    let mut line = Vec::new();

    loop {
        line.clear();

        let n = reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("reading {src}"))?;

        if n == 0 {
            break;
        }

        let text = line.strip_suffix(b"\n").unwrap_or(&line);
        let text = text.strip_suffix(b"\r").unwrap_or(text);

        if let Some(definition) = text.strip_prefix(b">") {
            if let Some(sequence) = current.take() {
                sequence.write(&mut out, uri.as_deref(), &args, alternates.as_ref())?;
            }

            seen_sequence = true;

            // As in samtools, the name ends at the first whitespace.
            let name_end = definition
                .iter()
                .position(u8::is_ascii_whitespace)
                .unwrap_or(definition.len());

            current = Some(Sequence::new(&definition[..name_end]));
        } else if let Some(sequence) = &mut current {
            sequence.add(text);
        }
    }

    match current {
        Some(sequence) => sequence.write(&mut out, uri.as_deref(), &args, alternates.as_ref())?,
        // samtools writes the empty dictionary silently.
        None if !seen_sequence => eprintln!("[W::dict] No sequences found in {src}"),
        None => {}
    }

    out.flush()?;

    Ok(())
}

/// Opens a FASTA file (`-` for stdin), plain or gzip or BGZF compressed, as samtools reads it
/// through zlib's `gzopen`.
fn open_fasta(src: &str) -> anyhow::Result<Box<dyn BufRead>> {
    let input = Input::open(src)?;

    Ok(match input.compression {
        Compression::None => Box::new(input.reader),
        Compression::Bgzf | Compression::Gzip => Box::new(BufReader::with_capacity(
            BUFFER_CAPACITY,
            MultiGzDecoder::new(input.reader),
        )),
    })
}

/// A sequence being read: its name, length, and the MD5 of its uppercased bases.
struct Sequence {
    name: Vec<u8>,
    length: usize,
    hasher: Md5,
    upper: Vec<u8>,
}

impl Sequence {
    fn new(name: &[u8]) -> Self {
        Self {
            name: name.to_vec(),
            length: 0,
            hasher: Md5::new(),
            upper: Vec::new(),
        }
    }

    fn add(&mut self, line: &[u8]) {
        self.upper.clear();
        self.upper.extend(
            line.iter()
                .filter(|b| !b.is_ascii_whitespace())
                .map(u8::to_ascii_uppercase),
        );

        self.length += self.upper.len();
        self.hasher.update(&self.upper);
    }

    fn write<W: Write>(
        self,
        out: &mut W,
        uri: Option<&str>,
        args: &Args,
        alternates: Option<&HashSet<Vec<u8>>>,
    ) -> io::Result<()> {
        let digest = self.hasher.finalize();

        out.write_all(b"@SQ\tSN:")?;
        out.write_all(&self.name)?;
        write!(out, "\tLN:{}\tM5:", self.length)?;

        for byte in digest {
            write!(out, "{byte:02x}")?;
        }

        if alternates.is_some_and(|names| names.contains(&self.name)) {
            out.write_all(b"\tAH:*")?;
        }

        if args.alias {
            write_alias(out, &self.name)?;
        }

        if let Some(uri) = uri {
            write!(out, "\tUR:{uri}")?;
        }

        if let Some(assembly) = &args.assembly {
            write!(out, "\tAS:{assembly}")?;
        }

        if let Some(species) = &args.species {
            write!(out, "\tSP:{species}")?;
        }

        writeln!(out)
    }
}

/// Writes the AN tag of `-A`: the name without its `chr` prefix, or with one added, and for
/// mitochondria both the `M` and `MT` forms.
fn write_alias<W: Write>(out: &mut W, name: &[u8]) -> io::Result<()> {
    let base = match name.strip_prefix(b"chr") {
        Some(base) => {
            out.write_all(b"\tAN:")?;
            out.write_all(base)?;
            base
        }
        None => {
            out.write_all(b"\tAN:chr")?;
            out.write_all(name)?;
            name
        }
    };

    match base {
        b"M" => out.write_all(b",chrMT,MT"),
        b"MT" => out.write_all(b",chrM,M"),
        _ => Ok(()),
    }
}

/// Reads the sequence names of an alternate locus file: the first field of each line that is
/// not a header (`@`) line. The file may be gzip-compressed.
fn read_alt_names(path: &Path) -> anyhow::Result<HashSet<Vec<u8>>> {
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let file = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut text = Vec::new();

    if reader.fill_buf()?.starts_with(&GZIP_MAGIC_NUMBER) {
        MultiGzDecoder::new(reader).read_to_end(&mut text)?;
    } else {
        reader.read_to_end(&mut text)?;
    }

    Ok(text
        .split(|&b| b == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
        .filter(|line| !line.is_empty() && !line.starts_with(b"@"))
        .map(|line| line.split(|&b| b == b'\t').next().unwrap_or(line).to_vec())
        .collect())
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::write::GzEncoder;
    use pasta::bgzf;

    use super::*;

    #[test]
    fn test_open_fasta_decompresses() -> anyhow::Result<()> {
        let fasta = b">a desc\nACGT\nac\n>b\nGG\n";

        let mut gzip = GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(fasta)?;
        let gzip = gzip.finish()?;

        let mut bgzf = bgzf::io::Writer::new(Vec::new());
        bgzf.write_all(fasta)?;
        let bgzf = bgzf.finish()?;

        let dir = tempfile::tempdir()?;

        for (name, data) in [
            ("plain.fa", fasta.to_vec()),
            ("gzip.fa.gz", gzip),
            ("bgzf.fa.gz", bgzf),
        ] {
            let path = dir.path().join(name);
            std::fs::write(&path, data)?;

            let mut text = Vec::new();
            open_fasta(path.to_str().context("path")?)?.read_to_end(&mut text)?;
            assert_eq!(text, fasta, "{name}");
        }

        Ok(())
    }

    #[test]
    fn test_open_fasta_truncated_gzip() -> anyhow::Result<()> {
        let mut gzip = GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(b">a\nACGT\n")?;
        let gzip = gzip.finish()?;

        let dir = tempfile::tempdir()?;
        let path = dir.path().join("t.fa.gz");
        std::fs::write(&path, &gzip[..gzip.len() - 4])?;

        let mut text = Vec::new();
        assert!(
            open_fasta(path.to_str().context("path")?)?
                .read_to_end(&mut text)
                .is_err()
        );

        Ok(())
    }
}
