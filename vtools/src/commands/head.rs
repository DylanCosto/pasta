//! `vtools head`: print a VCF/BCF header and the first records (bcftools head).

use std::{
    io::{self, IsTerminal, Write},
    process::ExitCode,
};

use anyhow::bail;
use pasta::variant::{
    Header, HeaderTextFlavor, Record,
    io::{Compression, ReadError, Reader},
    log, vcf,
};

use crate::options::last_index;
use crate::raw_input::{EXEC_FORMAT_ERROR, report_exec_format_error, strerror};

const USAGE: &str = "
About: Displays VCF/BCF headers and optionally the first few variant records
Usage: vtools head [OPTION]... [FILE]

Options:
  -h, --headers INT      Display INT header lines [all]
  -n, --records INT      Display INT variant record lines [none]
  -s, --samples INT      Display INT records starting with the #CHROM header line [none]
  -v, --verbosity INT    Verbosity level

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file [default: stdin].
    inputs: Vec<String>,

    /// Display INT header lines [all].
    #[arg(
        short = 'h',
        long = "headers",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    headers: Option<String>,

    /// Display INT variant record lines [none].
    #[arg(
        short = 'n',
        long = "records",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    records: Option<String>,

    /// Display INT records starting with the #CHROM header line [none].
    #[arg(
        short = 's',
        long = "samples",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    verbosity: Option<String>,

    /// Print help.
    #[arg(long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// C's `strtoull(s, NULL, 0)`: decimal, `0x` hex, or `0` octal; 0 when nothing parses.
pub fn strtoull0(s: &str) -> u64 {
    let s = s.trim_start();
    let (neg, s) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    let (radix, digits) = if (s.starts_with("0x") || s.starts_with("0X"))
        && s[2..].chars().next().is_some_and(|c| c.is_ascii_hexdigit())
    {
        (16, &s[2..])
    } else if s.starts_with('0') {
        (8, s)
    } else {
        (10, s)
    };
    let end = digits
        .char_indices()
        .find(|(_, c)| !c.is_digit(radix))
        .map_or(digits.len(), |(i, _)| i);
    let v =
        u64::from_str_radix(&digits[..end], radix).unwrap_or(if end == 0 { 0 } else { u64::MAX });
    if neg { v.wrapping_neg() } else { v }
}

/// C's `strtol(s, &end, 10)` checked with `*end == '\0'`: the value when the whole text is a
/// number, after optional white space and sign. An empty text is 0, as in C.
pub fn strtol_whole(s: &str) -> Option<i64> {
    let (v, end) = pasta::variant::number::strtol(s.as_bytes(), 0);
    (end == s.len()).then_some(v)
}

/// Applies bcftools' `--verbosity` (`apply_verbosity`): only levels above 3 change anything.
pub fn apply_verbosity(value: &str) -> anyhow::Result<()> {
    parse_verbosity(value).map(|_| ())
}

/// Applies a `--verbosity` value as [`apply_verbosity`] does and returns the level. The value
/// is read with C's `strtol` into an `int`, as bcftools reads it: leading spaces and a sign
/// are allowed, an empty value is 0, and the rest of the value must be digits.
pub fn parse_verbosity(value: &str) -> anyhow::Result<i32> {
    let b = value.as_bytes();
    let (v, end) = pasta::variant::number::strtol(b, 0);
    // C's conversion of the `long` to `int` keeps the low 32 bits.
    let v = v as i32;
    if end != b.len() || v < 0 {
        bail!("Could not parse argument: --verbosity {value}");
    }
    if v > 3 {
        log::set_verbosity(v);
    }
    Ok(v)
}

/// The header lines to print.
#[derive(Clone, Copy, Debug)]
struct HeaderLines {
    /// How many of the first lines to print; `None` prints them all.
    first: Option<u64>,
    /// Whether the `#CHROM` line is printed too, if it is not among the first lines (`-s`).
    chrom_line: bool,
}

pub fn run(args: Args, matches: &clap::ArgMatches) -> anyhow::Result<ExitCode> {
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    let mut nrecords = 0u64;
    // -n and -s both set the record count: the later one wins.
    let last = |id: &str| last_index(matches, id);
    if let Some(s) = &args.samples {
        nrecords = strtoull0(s);
    }
    if let Some(n) = &args.records
        && (args.samples.is_none() || last("records") > last("samples"))
    {
        nrecords = strtoull0(n);
    }
    let chrom_line = args.samples.is_some();
    let header_lines = HeaderLines {
        first: match &args.headers {
            Some(h) => Some(strtoull0(h)),
            None if chrom_line => Some(0),
            None => None,
        },
        chrom_line,
    };
    if args.inputs.is_empty() && io::stdin().is_terminal() {
        print!("{USAGE}");
        return Ok(ExitCode::SUCCESS);
    } else if args.inputs.len() > 1 {
        eprint!("{USAGE}");
        return Ok(ExitCode::FAILURE);
    }
    let fname = args.inputs.first().map_or("-", String::as_str);
    let (mut reader, mut header) = open(fname)?;

    let mut out = io::BufWriter::new(crate::output::stdout());
    let mut text = Vec::new();
    header.format_into(&mut text, HeaderTextFlavor::Vcf);
    write_header_lines(&mut out, &text, header_lines)?;
    if nrecords > 0 {
        write_records(
            &mut out,
            &mut reader,
            &mut header,
            fname,
            nrecords,
            &mut text,
        )?;
    }
    out.flush()?;
    Ok(ExitCode::SUCCESS)
}

/// Opens the input and reads its header, failing as bcftools head does.
fn open(fname: &str) -> anyhow::Result<(Reader, Header)> {
    let cant_read_headers = || {
        if fname == "-" {
            anyhow::anyhow!("[main_vcfhead] Can't read headers")
        } else {
            anyhow::anyhow!("[main_vcfhead] Can't read headers from \"{fname}\"")
        }
    };
    let mut reader = match Reader::open(fname, 0) {
        Ok(r) => r,
        // HTSlib opens any text (or a format it knows) and fails when reading the header; it
        // refuses binary data it does not recognize.
        Err(ReadError::UnknownBinary) => {
            report_exec_format_error(fname);
            bail!("[main_vcfhead] Can't open \"{fname}\": {EXEC_FORMAT_ERROR}");
        }
        Err(ReadError::UnknownFormat) => {
            log::error(
                "bcf_hdr_read",
                format_args!("Input is not detected as bcf or vcf format"),
            );
            return Err(cant_read_headers());
        }
        Err(e) => {
            let msg = strerror(&e);
            if fname == "-" {
                bail!("[main_vcfhead] Can't open standard input: {msg}");
            }
            if let ReadError::Io(_) = e {
                log::error(
                    "hts_open_format",
                    format_args!("Failed to open file \"{fname}\" : {msg}"),
                );
            }
            bail!("[main_vcfhead] Can't open \"{fname}\": {msg}");
        }
    };
    let Ok(header) = reader.read_header() else {
        return Err(cant_read_headers());
    };
    Ok((reader, header))
}

/// Writes the selected lines of the header text.
fn write_header_lines(out: &mut impl Write, text: &[u8], lines: HeaderLines) -> io::Result<()> {
    let Some(nheaders) = lines.first else {
        return out.write_all(cstr(text));
    };
    let samples = lines.chrom_line;
    if nheaders == 0 && !samples {
        return Ok(());
    }
    let mut lim = Some(0usize);
    let mut samples_printed = false;
    for _ in 0..nheaders {
        let Some(l) = lim else { break };
        if samples && text[l..].starts_with(b"#CHROM\t") {
            samples_printed = true;
        }
        lim = memchr::memchr(b'\n', &text[l..]).map(|i| l + i + 1);
    }
    if nheaders > 0 {
        let end = lim.unwrap_or(text.len());
        out.write_all(cstr(&text[..end]))?;
    }
    if let Some(mut l) = lim
        && samples
        && !samples_printed
    {
        while l < text.len() {
            if text[l..].starts_with(b"#CHROM\t") {
                out.write_all(cstr(&text[l..]))?;
                break;
            }
            match memchr::memchr(b'\n', &text[l..]) {
                Some(i) => l += i + 1,
                None => break,
            }
        }
    }
    Ok(())
}

/// Writes up to `nrecords` records as VCF lines, formatting each into `text`. An unparsable
/// VCF line ends the records without failing, as in bcftools.
///
/// # Errors
///
/// Fails after the records before it for a failed read (a truncated or corrupt BGZF block, a
/// malformed BCF record), which bcftools reports and takes for the end with status 0. Also
/// returns write errors.
fn write_records(
    out: &mut impl Write,
    reader: &mut Reader,
    header: &mut Header,
    fname: &str,
    nrecords: u64,
    text: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let mut rec = Record::new();
    let mut line = Vec::new();
    let mut n = 0u64;
    while n < nrecords {
        let (result, msgs) = log::capture(|| reader.read_record(header, &mut line, &mut rec));
        let reported = !msgs.is_empty();
        log::replay(msgs);
        match result {
            Ok(true) => {}
            Ok(false) => {
                if reader.compression() == Compression::Bgzf {
                    super::index::warn_if_no_eof(fname);
                }
                break;
            }
            Err(e) => {
                // HTSlib reports every failed read; a failed BGZF block (a truncated or
                // corrupt file) is reported by `bgzf_read_block`.
                if !reported {
                    report_read_error(&e, reader.compression());
                }
                if let ReadError::Io(_) | ReadError::Bcf(_) = e {
                    out.flush()?;
                    if fname == "-" {
                        bail!("[main_vcfhead] Error reading standard input");
                    }
                    bail!("[main_vcfhead] Error reading \"{fname}\"");
                }
                break;
            }
        }
        text.clear();
        if vcf::format_record(header, &rec, text).is_ok() {
            out.write_all(cstr(text))?;
        } else {
            eprintln!("[main_vcfhead] Record #{} is invalid", n + 1);
        }
        n += 1;
    }
    Ok(())
}

/// Reports a failed read of compressed input that nothing has reported yet (HTSlib's
/// `bgzf_read_block` errors). Other failures, such as a VCF line with too few columns, are
/// silent in HTSlib too.
fn report_read_error(e: &ReadError, compression: Compression) {
    if let ReadError::Io(e) | ReadError::Bcf(e) = e
        && compression != Compression::None
    {
        log::error(
            "bgzf_read_block",
            format_args!("Failed to read BGZF block data: {e}"),
        );
    }
}

/// The text up to its first NUL, as `fputs` writes it.
fn cstr(s: &[u8]) -> &[u8] {
    &s[..memchr::memchr(0, s).unwrap_or(s.len())]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verbosity_is_read_as_strtol_reads_it() {
        for (value, level) in [("", 0), ("0", 0), (" 2", 2), ("+3", 3), ("1", 1)] {
            assert_eq!(parse_verbosity(value).ok(), Some(level), "value {value:?}");
        }
        for value in ["x", "3x", " ", "-1", "9223372036854775807"] {
            assert!(parse_verbosity(value).is_err(), "value {value:?}");
        }
    }
}
