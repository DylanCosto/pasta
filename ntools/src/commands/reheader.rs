//! `ntools reheader`: replace a BAM or CRAM file's header without recompressing records
//! (samtools reheader).

mod cram;

use std::{
    env,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use anyhow::{Context, bail};
use pasta::sam;

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    bam_blocks,
    cram_file::Containers,
    output, program,
    source::{Format, Input},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// SAM file whose header lines become the new header, then the input BAM or CRAM file
    /// (only the input with -c).
    #[arg(required = true, num_args = 1..=2)]
    files: Vec<PathBuf>,

    /// Do not add an @PG line to the header.
    #[arg(short = 'P', long = "no-PG")]
    no_pg: bool,

    /// Modify the CRAM file directly, if the new header fits.
    #[arg(short = 'i', long = "in-place")]
    in_place: bool,

    /// Pass the input's header, as SAM, to the standard input of the shell command CMD, and use
    /// its output as the new header.
    #[arg(short = 'c', long = "command", value_name = "CMD")]
    command: Option<String>,
}

/// Where the new header comes from.
enum HeaderSource<'a> {
    /// A SAM file's header lines.
    File(Box<sam::Header>),
    /// The output of a shell command given the input's header (`-c`).
    Command(&'a str),
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let (source, input) = match (&args.command, args.files.as_slice()) {
        (Some(command), [input]) => (HeaderSource::Command(command), input),
        (None, [header, input]) => (
            HeaderSource::File(Box::new(read_sam_header(header)?)),
            input,
        ),
        _ => bail!(
            "Usage: ntools reheader [-P] in.header.sam in.bam > out.bam\n   \
             or  ntools reheader [-P] -i in.header.sam file.cram\n   \
             or  ntools reheader -c CMD in.bam\n   \
             or  ntools reheader -c CMD in.cram"
        ),
    };

    // The input is read once, in order, so it can be a pipe.
    let src = input.to_str().context("input path is not valid UTF-8")?;
    let Input { reader, format, .. } = Input::open(src)?;

    let new_header = |input_header: &sam::Header| -> anyhow::Result<sam::Header> {
        let mut header = match &source {
            HeaderSource::File(header) => sam::Header::clone(header),
            HeaderSource::Command(command) => external_header(input_header, command)?,
        };

        if !args.no_pg {
            program::add(&mut header)?;
        }

        Ok(header)
    };

    let read_error = || format!("failed to read the header for '{src}'");

    match format {
        Format::Sam => bail!("input file '{src}' must be BAM or CRAM"),
        Format::Bam if args.in_place => bail!("cannot reheader BAM '{src}' in-place"),
        Format::Bam => {
            let (input_header, body) = bam_blocks::split(reader, input).with_context(read_error)?;
            let header = new_header(&input_header)?;

            let mut dst = BufWriter::new(output::raw_stdout()?);
            bam_blocks::write_header(&mut dst, &header)?;
            body.copy_to(&mut dst)?;
            bam_blocks::write_eof(&mut dst)?;
            dst.flush()?;

            Ok(())
        }
        Format::Cram => {
            let (start, mut containers) = Containers::open(reader).with_context(read_error)?;
            let header = new_header(&start.header().with_context(read_error)?)?;

            if args.in_place {
                return cram::reheader_in_place(input, &header);
            }

            let mut dst = BufWriter::new(output::raw_stdout()?);
            cram::reheader(&start, &mut containers, &header, &mut dst)
                .with_context(|| format!("copying '{src}'"))?;
            dst.flush()?;

            Ok(())
        }
    }
}

/// Runs `( command ) < header.sam` with the input's header in a temporary file, and reads the
/// new header from its output (`-c`).
fn external_header(header: &sam::Header, command: &str) -> anyhow::Result<sam::Header> {
    let temp_dir = env::var_os("TMPDIR").map_or_else(|| PathBuf::from("/tmp"), PathBuf::from);
    let mut temp = tempfile::Builder::new()
        .prefix("reheader")
        .tempfile_in(&temp_dir)
        .with_context(|| format!("fail to open temp file in '{}'", temp_dir.display()))?;

    sam::io::Writer::new(&mut temp)
        .write_header(header)
        .context("failed to write the header to the temp file")?;
    temp.flush()?;

    let output = Command::new("sh")
        .arg("-c")
        .arg(format!("( {command} ) < {}", shell_quote(temp.path())))
        .stdin(Stdio::null())
        .stderr(Stdio::inherit())
        .output()
        .with_context(|| format!("failed to run external command '{command}'"))?;

    if !output.status.success() {
        bail!("Non-zero exit code returned by command '{command}'");
    }

    sam::io::Reader::new(output.stdout.as_slice())
        .read_header()
        .with_context(|| format!("failed to read the header from '{command}'"))
}

/// Quotes a path for `sh`.
fn shell_quote(path: &Path) -> String {
    format!("'{}'", path.display().to_string().replace('\'', r"'\''"))
}

/// Reads the header of a SAM, BAM, or CRAM file.
pub fn read_sam_header(path: &Path) -> anyhow::Result<sam::Header> {
    let src = path.to_str().context("header path is not valid UTF-8")?;
    let (mut reader, _) = AlignmentReader::open(src, &ReadOptions::default())?;

    reader
        .read_header()
        .with_context(|| format!("reading header from {}", path.display()))
}
