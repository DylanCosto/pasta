//! `vtools reheader`: replace the header, rename samples, or update contigs from a `.fai`
//! (bcftools reheader).
//!
//! Bgzipped VCF keeps its record blocks as they are: only the header blocks are rewritten and
//! the rest of the file is copied without recompressing, as bcftools does. Plain VCF is copied
//! line by line. BCF records are checked against the new header and written again.

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, IsTerminal, Write},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, HeaderTextFlavor, Quoting, Record, Unpack, bcf,
    index::split_idx_name,
    io::{OutputType, WriteOptions, Writer},
    log,
    readlist::{self, ListSource},
};

use super::head::apply_verbosity;
use crate::options::last_index;
use crate::{
    output,
    raw_input::{
        Compression, Format, FrameCopier, FrameReader, InflatingReader, RawInput,
        report_exec_format_error, strerror,
    },
};

const USAGE: &str = "
About:   Modify header of VCF/BCF files, change sample names.
Usage:   vtools reheader [OPTIONS] <in.vcf.gz>

Options:
    -f, --fai FILE             Update sequences and their lengths from the .fai file
    -h, --header FILE          New header
    -o, --output FILE          Write output to a file [standard output]
    -n, --samples-list LIST    New sample names given as a comma-separated list
    -N, --samples-file FILE    New sample names in a file, see the man page for details
    -T, --temp-prefix PATH     Ignored; was template for temporary file name
        --threads INT          Use multithreading with INT worker threads (BCF only) [0]
    -v, --verbosity INT        Verbosity level

Example:
   # Write out the header to be modified
   vtools view -h old.bcf > header.txt

   # Edit the header using your favorite text editor
   vi header.txt

   # Reheader the file
   vtools reheader -h header.txt -o new.bcf old.bcf

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file [default: stdin].
    inputs: Vec<String>,

    /// Update sequences and their lengths from the .fai file.
    #[arg(
        short = 'f',
        long = "fai",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    fai: Option<String>,

    /// New header.
    #[arg(
        short = 'h',
        long = "header",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    header: Option<String>,

    /// Write output to a file [standard output].
    #[arg(
        short = 'o',
        long = "output",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<PathBuf>,

    /// New sample names in a file (also `-s`).
    #[arg(
        short = 'N',
        long = "samples-file",
        short_alias = 's',
        alias = "samples",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    samples_file: Option<String>,

    /// New sample names given as a comma-separated list.
    #[arg(
        short = 'n',
        long = "samples-list",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    samples_list: Option<String>,

    /// Ignored; was template for temporary file name.
    #[arg(
        short = 'T',
        long = "temp-prefix",
        value_name = "PATH",
        allow_hyphen_values = true
    )]
    _temp_prefix: Option<String>,

    /// Use multithreading with INT worker threads (BCF only) [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    verbosity: Option<String>,

    /// Print help.
    #[arg(long = "help")]
    help: bool,
}

/// The new sample names: a file or a list (the last of `-s/-N/-n` wins).
struct SampleSpec {
    spec: String,
    source: ListSource,
}

/// The sequences of a `.fai` file, in file order (HTSlib's `faidx_t` for lengths).
struct Fai {
    names: Vec<Vec<u8>>,
    lengths: HashMap<Vec<u8>, i64>,
}

impl Fai {
    /// Reads a `.fai` file as HTSlib's `fai_load3` and `fai_read` do: NAME (up to the first
    /// whitespace), then LENGTH, OFFSET, LINEBASES, and LINEWIDTH as unsigned numbers. Logs
    /// HTSlib's errors and returns `None` when the file is unreadable or a line is malformed.
    fn load(path: &str) -> Option<Self> {
        let data = match std::fs::read(path) {
            Ok(d) => d,
            Err(e) => {
                // HTSlib tries to build the missing index from the file itself, which fails.
                let msg = e.to_string();
                let msg = msg.find(" (os error").map_or(msg.as_str(), |i| &msg[..i]);
                log::error(
                    "fai_build3_core",
                    format_args!("Failed to open the file {path} : {msg}"),
                );
                return None;
            }
        };
        let fai = Self::parse(path, &data);
        if fai.is_none() {
            log::error(
                "fai_load3_core",
                format_args!("Failed to read FASTA index {path}"),
            );
        }
        fai
    }

    /// HTSlib's `fai_read` and `fai_insert_index` on the bytes of a `.fai` file.
    fn parse(path: &str, data: &[u8]) -> Option<Self> {
        let mut names = Vec::new();
        let mut lengths = HashMap::new();
        let mut lnum = 1usize;
        for line in data.split_inclusive(|&c| c == b'\n') {
            let name_end = line
                .iter()
                .position(|&c| c == 0 || is_c_space(c))
                .unwrap_or(line.len());
            let name = &line[..name_end];
            let rest = line.get(name_end + 1..).unwrap_or_default();
            let mut nums = [0u64; 4];
            let mut p = 0;
            for n in &mut nums {
                let Some((v, end)) = scan_u64(rest, p) else {
                    log::error(
                        "fai_read",
                        format_args!("Could not understand FASTA index {path} line {lnum}"),
                    );
                    return None;
                };
                *n = v;
                p = end;
            }
            let [len, seq_offset, line_blen, line_len] = nums;
            const HTS_POS_MAX: u64 = (1 << 62) - 1;
            if line_len > HTS_POS_MAX || line_blen > HTS_POS_MAX {
                log::error(
                    "fai_insert_index",
                    format_args!("Malformed line: length out of bounds"),
                );
                return None;
            }
            if line_len < line_blen {
                log::error(
                    "fai_insert_index",
                    format_args!(
                        "Malformed line: width ({line_len}) less than base count ({line_blen})"
                    ),
                );
                return None;
            }
            if line_len - line_blen > u64::from(u32::MAX) {
                log::error(
                    "fai_insert_index",
                    format_args!(
                        "Malformed line: difference between width ({line_len}) and base count ({line_blen}) too large"
                    ),
                );
                return None;
            }
            if lengths.contains_key(name) {
                log::warn(
                    "fai_insert_index",
                    format_args!(
                        "Ignoring duplicate sequence \"{}\" at byte offset {seq_offset}",
                        String::from_utf8_lossy(name)
                    ),
                );
            } else {
                names.push(name.to_vec());
                // HTSlib keeps the length as an unsigned number and prints it as a signed one.
                lengths.insert(name.to_vec(), len as i64);
            }
            if line.last() == Some(&b'\n') {
                lnum += 1;
            }
        }
        Some(Self { names, lengths })
    }

    fn len(&self, name: &[u8]) -> Option<i64> {
        self.lengths.get(name).copied()
    }
}

pub fn run(args: Args, matches: &clap::ArgMatches) -> anyhow::Result<ExitCode> {
    if args.help {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    // A negative count means no threads.
    let threads = crate::output::threads(args.threads.as_deref()).unwrap_or(0);
    let fname = match args.inputs.first() {
        Some(f) => f.as_str(),
        None if !io::stdin().is_terminal() => "-",
        None => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
    };
    let fai = match &args.fai {
        Some(f) => match Fai::load(f) {
            Some(fai) => Some(fai),
            None => bail!("Could not parse {f}"),
        },
        None => None,
    };
    // -s/-N and -n set the same option: the last one given wins.
    let last = |id: &str| last_index(matches, id);
    let samples = match (&args.samples_list, &args.samples_file) {
        (Some(_), Some(file)) if last("samples_file") > last("samples_list") => Some(SampleSpec {
            spec: file.clone(),
            source: ListSource::File,
        }),
        (Some(list), _) => Some(SampleSpec {
            spec: list.clone(),
            source: ListSource::Inline,
        }),
        (None, Some(file)) => Some(SampleSpec {
            spec: file.clone(),
            source: ListSource::File,
        }),
        (None, None) => None,
    };
    if samples.is_none() && args.header.is_none() && fai.is_none() {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    let edits = Edits {
        samples,
        header: args.header,
        fai,
    };

    let read_failed = || {
        let name = if fname == "-" {
            "standard input"
        } else {
            fname
        };
        anyhow::anyhow!("Failed to read from {name}")
    };
    // `hts_open` reads `FILE` of `FILE##idx##INDEX`.
    let data = split_idx_name(fname).0;
    let input = match RawInput::open(data) {
        Ok(i) if i.is_unknown_binary() => {
            report_exec_format_error(data);
            return Err(read_failed());
        }
        Ok(i) => i,
        Err(e) => {
            output::log_open_failure(data, &e);
            return Err(read_failed());
        }
    };
    let out = args.output.as_deref();
    match (input.format, input.compression) {
        (Format::Vcf, Compression::Bgzf) => reheader_vcf_gz(input, fname, out, &edits)?,
        (Format::Vcf, Compression::None) => reheader_vcf(input, fname, out, &edits)?,
        (Format::Vcf, Compression::Gzip) => bail!(
            "Error: cannot reheader gzip-compressed files, first convert with `vtools view --output-type` to a supported format"
        ),
        _ => reheader_bcf(input, fname, out, &edits, threads)?,
    }
    Ok(ExitCode::SUCCESS)
}

/// The requested changes to the header text.
struct Edits {
    samples: Option<SampleSpec>,
    header: Option<String>,
    fai: Option<Fai>,
}

impl Edits {
    /// Applies the edits to the header text as bcftools does: the new header file replaces it,
    /// then the `.fai` updates contigs, then samples are renamed or replaced.
    fn apply(&self, hdr: &mut Vec<u8>) -> anyhow::Result<()> {
        let samples = match &self.samples {
            Some(s) => Some(read_samples(s)?),
            None => None,
        };
        if let Some(h) = &self.header {
            *hdr = read_header_file(h)?;
        }
        if let Some(fai) = &self.fai {
            *hdr = update_from_fai(fai, hdr)?;
        }
        if let Some(samples) = samples {
            set_samples(&samples, hdr)?;
        }
        Ok(())
    }
}

/// bcftools' `read_samples` (`hts_readlist`). An empty file is an empty list, which leaves
/// the samples as they are.
fn read_samples(s: &SampleSpec) -> anyhow::Result<Vec<Vec<u8>>> {
    match readlist::read_list(s.spec.as_bytes(), s.source) {
        Ok(list) => Ok(list),
        Err(_) => {
            let (opt, what) = match s.source {
                ListSource::File => ("--samples-file", "file"),
                ListSource::Inline => ("--samples-list", "list"),
            };
            bail!("Error parsing the {opt} {what} \"{}\"", s.spec)
        }
    }
}

/// bcftools' `read_header_file`: the lines up to the first empty one, trailing whitespace
/// removed, ending with a newline. `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
fn read_header_file(fname: &str) -> anyhow::Result<Vec<u8>> {
    let data = split_idx_name(fname).0;
    let lines = match readlist::read_all_lines(Path::new(data)) {
        Ok(lines) => lines,
        Err(e) => {
            output::log_open_failure(data, &e);
            bail!("Could not read: {fname}");
        }
    };
    let mut hdr = Vec::new();
    for line in lines {
        if line.is_empty() {
            break;
        }
        hdr.extend_from_slice(&line);
        hdr.push(b'\n');
    }
    while hdr.last().is_some_and(|c| c.is_ascii_whitespace()) {
        hdr.pop();
    }
    hdr.push(b'\n');
    Ok(hdr)
}

/// HTSlib's `is_escaped`: whether `s[i]` follows an odd number of backslashes after `min`.
fn is_escaped(s: &[u8], min: usize, i: usize) -> bool {
    let mut n = 0;
    let mut j = i;
    while j > min && s[j - 1] == b'\\' {
        n += 1;
        j -= 1;
    }
    n % 2 == 1
}

fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// bcftools' `copy_and_update_contig_line`: rewrites the `##contig` line at `line` with the
/// `.fai` length into `dst`, or drops it. Returns the offset where copying resumes.
fn copy_and_update_contig_line(
    fai: &Fai,
    s: &[u8],
    line: usize,
    seen: &mut std::collections::HashSet<Vec<u8>>,
    dst: &mut Vec<u8>,
) -> anyhow::Result<usize> {
    let mut q = line + 9;
    let mut end = q;
    while at(s, end) != 0 && at(s, end) != b'\n' {
        end += 1;
    }
    let mut nopen = 1;
    let mut chr: Option<(Vec<u8>, i64)> = None;
    let mut tmp: Vec<u8> = Vec::new();
    while at(s, q) != 0 && at(s, q) != b'\n' && nopen > 0 {
        q += 1;
        let mut p = q;
        while matches!(at(s, q), b' ' | b'\t') {
            p += 1;
            q += 1;
        }
        if p == q && (at(s, q).is_ascii_alphabetic() || at(s, q) == b'_') {
            q += 1;
            while at(s, q).is_ascii_alphanumeric() || matches!(at(s, q), b'_' | b'.') {
                q += 1;
            }
        }
        let n = q - p;
        let mut m = 0;
        while matches!(at(s, q), b' ' | b'\t') {
            q += 1;
            m += 1;
        }
        if at(s, q) != b'=' || n == 0 {
            let cut = |from: usize| String::from_utf8_lossy(&s[from.min(end)..end]).into_owned();
            bail!(
                "Could not parse the line: {} [{}][{}]",
                cut(line),
                cut(p),
                cut(q)
            );
        }
        let key = s[p..q - m].to_vec();
        q += 1;
        p = q;
        while matches!(at(s, q), b' ' | b'\t') {
            p += 1;
            q += 1;
        }
        let quoted = at(s, p) == b'"';
        if quoted {
            p += 1;
            q += 1;
        }
        while at(s, q) != 0 && at(s, q) != b'\n' {
            if quoted {
                if at(s, q) == b'"' && !is_escaped(s, p, q) {
                    break;
                }
            } else {
                if at(s, q) == b'<' {
                    nopen += 1;
                }
                if at(s, q) == b'>' {
                    nopen -= 1;
                }
                if nopen == 0 || (at(s, q) == b',' && nopen == 1) {
                    break;
                }
            }
            q += 1;
        }
        let mut r = q;
        while r > p && s[r - 1] == b' ' {
            r -= 1;
        }
        let val = &s[p..r];
        if quoted && at(s, q) == b'"' {
            q += 1;
        }
        if at(s, q) == b'>' {
            nopen -= 1;
            q += 1;
        }
        if key == b"length" {
            continue;
        }
        if key == b"ID" {
            if seen.contains(val) {
                continue;
            }
            let Some(len) = fai.len(val) else {
                return Ok(end);
            };
            seen.insert(val.to_vec());
            chr = Some((val.to_vec(), len));
            continue;
        }
        tmp.push(b',');
        tmp.extend_from_slice(&key);
        tmp.push(b'=');
        if quoted {
            tmp.push(b'"');
        }
        tmp.extend_from_slice(val);
        if quoted {
            tmp.push(b'"');
        }
    }
    let Some((name, len)) = chr else {
        return Ok(end);
    };
    dst.extend_from_slice(b"##contig=<ID=");
    dst.extend_from_slice(&name);
    dst.extend_from_slice(format!(",length={len}").as_bytes());
    dst.extend_from_slice(&tmp);
    dst.push(b'>');
    Ok(q)
}

/// bcftools' `update_from_fai`: updates the lengths of contigs in the `.fai`, removes the
/// others, and adds the missing ones before `#CHROM`.
fn update_from_fai(fai: &Fai, hdr: &[u8]) -> anyhow::Result<Vec<u8>> {
    let mut out = Vec::with_capacity(hdr.len());
    let mut seen = std::collections::HashSet::new();
    let mut beg = 0usize;
    while beg < hdr.len() {
        let Some(off) = memchr::memmem::find(&hdr[beg..], b"\n##contig=<") else {
            break;
        };
        let tmp = beg + off;
        out.extend_from_slice(&hdr[beg..=tmp]);
        let l_prev = out.len();
        beg = copy_and_update_contig_line(fai, hdr, tmp + 1, &mut seen, &mut out)?;
        if l_prev == out.len() {
            out.pop();
        }
    }
    let Some(off) = memchr::memmem::find(&hdr[beg.min(hdr.len())..], b"\n#CHROM") else {
        bail!("Failed to parse the header, #CHROM not found");
    };
    let tmp = beg + off;
    out.extend_from_slice(&hdr[beg..=tmp]);
    for name in &fai.names {
        if seen.contains(name) {
            continue;
        }
        out.extend_from_slice(b"##contig=<ID=");
        out.extend_from_slice(name);
        out.extend_from_slice(format!(",length={}>\n", fai.lengths[name]).as_bytes());
    }
    out.extend_from_slice(&hdr[tmp + 1..]);
    Ok(out)
}

/// bcftools' `set_sample_pairs`: when every item is an "old new" pair, renames the samples on
/// the `#CHROM` line starting at `idx`. Returns whether the items were pairs.
fn set_sample_pairs(samples: &[Vec<u8>], hdr: &mut Vec<u8>, idx: usize) -> bool {
    let mut map: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for s in samples {
        let mut key = Vec::new();
        let mut val = Vec::new();
        let mut escaped = false;
        let mut i = 0;
        while i < s.len() && s[i] != 0 {
            let c = s[i];
            if c == b'\\' && !escaped {
                escaped = true;
                i += 1;
                continue;
            }
            if is_c_space(c) && !escaped {
                break;
            }
            key.push(c);
            escaped = false;
            i += 1;
        }
        if i >= s.len() {
            return false;
        }
        while i < s.len() && is_c_space(s[i]) {
            i += 1;
        }
        while i < s.len() {
            let c = s[i];
            if c == b'\\' && !escaped {
                escaped = true;
                i += 1;
                continue;
            }
            if is_c_space(c) && !escaped {
                break;
            }
            val.push(c);
            escaped = false;
            i += 1;
        }
        map.insert(key, val);
    }
    while hdr.last().is_some_and(|&c| is_c_space(c)) {
        hdr.pop();
    }
    let line = hdr.split_off(idx);
    let ncols = line.split(|&c| c == b'\t').count();
    let mut out = Vec::with_capacity(line.len());
    for (n, col) in line.split(|&c| c == b'\t').enumerate() {
        // Columns 10 and on are samples; the last column is looked up whatever it is.
        let rename = n >= 9 || n + 1 == ncols;
        let new = if rename { map.get(col) } else { None };
        out.extend_from_slice(new.map_or(col, Vec::as_slice));
        if n + 1 < ncols {
            out.push(b'\t');
        }
    }
    hdr.extend_from_slice(&out);
    hdr.push(b'\n');
    true
}

/// C's `isspace`.
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
}

/// One `sscanf` `%llu` conversion from `s[i..]`: leading whitespace, an optional sign, then
/// decimal digits (negative values wrap and overflow saturates, as `strtoull` does). Returns
/// the value and the offset after it, or `None` when there are no digits.
fn scan_u64(s: &[u8], mut i: usize) -> Option<(u64, usize)> {
    while i < s.len() && is_c_space(s[i]) {
        i += 1;
    }
    let neg = match s.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let start = i;
    let mut v: u64 = 0;
    let mut overflow = false;
    while let Some(d) = s.get(i).filter(|c| c.is_ascii_digit()) {
        match v
            .checked_mul(10)
            .and_then(|v| v.checked_add(u64::from(d - b'0')))
        {
            Some(n) => v = n,
            None => overflow = true,
        }
        i += 1;
    }
    if i == start {
        return None;
    }
    let v = if overflow { u64::MAX } else { v };
    Some((if neg { v.wrapping_neg() } else { v }, i))
}

/// bcftools' `set_samples`: renames ("old new" pairs) or replaces the samples on the
/// `#CHROM` line.
fn set_samples(samples: &[Vec<u8>], hdr: &mut Vec<u8>) -> anyhow::Result<()> {
    const FULL: &[u8] = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    const SITES: &[u8] = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    // The start of the #CHROM line (the header ends with a newline).
    let mut i = hdr.len() as isize - 2;
    let mut ncols = 0usize;
    while i >= 0 && hdr[i as usize] != b'\n' {
        if hdr[i as usize] == b'\t' {
            ncols += 1;
        }
        i -= 1;
    }
    let start = (i + 1) as usize;
    if i < 0 || !hdr[start..].starts_with(FULL) {
        if i > 0 && hdr[start..].starts_with(SITES) {
            bail!("Error: missing FORMAT fields, cowardly refusing to add samples");
        }
        bail!(
            "Could not parse the header: {}",
            String::from_utf8_lossy(hdr)
        );
    }
    if set_sample_pairs(samples, hdr, start) {
        return Ok(());
    }
    if ncols != samples.len() + 8 {
        eprintln!(
            "Warning: different number of samples: {} vs {}",
            samples.len(),
            ncols as isize - 8
        );
    }
    // The samples start after the ninth tab; a line without samples ends after FORMAT
    // (bcftools reads past the end of the header there).
    let mut tabs = 0;
    let mut j = start;
    while j < hdr.len() && hdr[j] != b'\n' {
        if hdr[j] == b'\t' {
            tabs += 1;
            if tabs == 9 {
                break;
            }
        }
        j += 1;
    }
    hdr.truncate(j);
    for s in samples {
        hdr.push(b'\t');
        hdr.extend_from_slice(s);
    }
    hdr.push(b'\n');
    Ok(())
}

fn out_name(out: Option<&Path>) -> String {
    out.map_or_else(|| "-".to_owned(), |p| p.display().to_string())
}

/// bcftools' `reheader_vcf_gz`: rewrites the header blocks and copies the rest of the file
/// without recompressing it.
fn reheader_vcf_gz(
    input: RawInput,
    fname: &str,
    out: Option<&Path>,
    edits: &Edits,
) -> anyhow::Result<()> {
    let mut fr = FrameReader::new(input.source, split_idx_name(fname).0);
    let read_err = |e: io::Error| anyhow::anyhow!("Failed to read {fname}: {}", strerror(&e));
    fr.read_block().map_err(read_err)?;
    if fr.data.is_empty() {
        bail!("Failed to read {fname}: Success");
    }
    if fr.data[0] != b'#' {
        bail!(
            "Could not parse the header, expected '#', found '{}'",
            fr.data[0] as char
        );
    }
    let mut hdr = Vec::new();
    let mut skip = 1usize;
    let next = |fr: &mut FrameReader<_>| -> anyhow::Result<()> {
        fr.read_block()
            .map_err(|_| anyhow::anyhow!("Error reading {fname}"))
    };
    loop {
        if at(&fr.data, skip) == b'\n' {
            skip += 1;
            if skip >= fr.data.len() {
                hdr.extend_from_slice(&fr.data[..skip.min(fr.data.len())]);
                next(&mut fr)?;
                if fr.data.is_empty() {
                    break;
                }
                skip = 0;
            }
            if at(&fr.data, skip) != b'#' {
                hdr.extend_from_slice(&fr.data[..skip]);
                break;
            }
        }
        skip += 1;
        if skip >= fr.data.len() {
            hdr.extend_from_slice(&fr.data);
            next(&mut fr)?;
            if fr.data.is_empty() {
                break;
            }
            skip = 0;
        }
    }
    edits.apply(&mut hdr)?;

    let mut writer = Writer::create(
        out,
        &WriteOptions {
            output_type: OutputType::VcfBgzf,
            level: None,
            threads: 0,
        },
    )
    .map_err(|e| anyhow::anyhow!("{}: {}", out_name(out), strerror(&e)))?;
    writer.write_raw(&hdr)?;
    if fr.data.len() > skip {
        writer.write_raw(&fr.data[skip..])?;
    }
    // Copy the rest of the file as it is, in batches of whole frames. bcftools copies raw
    // bytes (`bgzf_raw_read`), so a truncated last block goes into the output, which is then
    // corrupt, with status 0. Reading frame by frame finds it without decompressing anything:
    // the complete frames before it are written and reheader fails.
    let mut buf = Vec::with_capacity(COPY_BATCH + (1 << 16));
    let mut failed = None;
    loop {
        match fr.next_raw_frame() {
            Ok(true) => buf.extend_from_slice(&fr.frame),
            Ok(false) => break,
            Err(e) => {
                failed = Some(e);
                break;
            }
        }
        if buf.len() >= COPY_BATCH {
            writer.write_bgzf_frames(&buf)?;
            buf.clear();
        }
    }
    writer.write_bgzf_frames(&buf)?;
    writer.finish().map_err(|e| {
        let msg = format!("Error closing {}: {}", out_name(out), strerror(&e));
        output::write_failure(e, || anyhow::anyhow!(msg))
    })?;
    if let Some(e) = failed {
        log::error(
            "bgzf_read_block",
            format_args!("Failed to read BGZF block data: {}", strerror(&e)),
        );
        bail!("Error reading {fname}");
    }
    Ok(())
}

/// How many bytes of whole BGZF frames [`reheader_vcf_gz`] copies at a time.
const COPY_BATCH: usize = 1 << 20;

/// bcftools' `reheader_vcf`: the header lines, edited, then the other lines as they are.
fn reheader_vcf(
    input: RawInput,
    fname: &str,
    out: Option<&Path>,
    edits: &Edits,
) -> anyhow::Result<()> {
    let mut src = input.source;
    let mut hdr = Vec::new();
    let mut line = Vec::new();
    let mut first: Option<Vec<u8>> = None;
    let read_err = |_| anyhow::anyhow!("[reheader_vcf] Error: close failed .. {fname}");
    while read_text_line(&mut src, &mut line).map_err(read_err)? {
        line.push(b'\n');
        if line[0] != b'#' {
            first = Some(std::mem::take(&mut line));
            break;
        }
        hdr.extend_from_slice(&line);
    }
    edits.apply(&mut hdr)?;
    let sink: Box<dyn Write> = match out {
        Some(p) => match File::create(p) {
            Ok(f) => Box::new(f),
            Err(e) => bail!("{}: {}", p.display(), strerror(&e)),
        },
        None => Box::new(output::stdout()),
    };
    let mut w = io::BufWriter::with_capacity(1 << 17, sink);
    let write_err = |n: usize| anyhow::anyhow!("Failed to write {n} bytes");
    w.write_all(&hdr)
        .map_err(|e| output::write_failure(e, || write_err(hdr.len())))?;
    if let Some(l) = first {
        w.write_all(&l)
            .map_err(|e| output::write_failure(e, || write_err(l.len())))?;
        while read_text_line(&mut src, &mut line).map_err(read_err)? {
            line.push(b'\n');
            w.write_all(&line)
                .map_err(|e| output::write_failure(e, || write_err(line.len())))?;
        }
    }
    w.flush().map_err(|e| {
        output::write_failure(e, || {
            anyhow::anyhow!("[reheader_vcf] Error: close failed .. {}", out_name(out))
        })
    })?;
    Ok(())
}

/// Reads a line without its newline (HTSlib's `hts_getline` on uncompressed text: a `\r`
/// before the newline is removed too). Returns `false` at the end.
fn read_text_line<R: BufRead>(r: &mut R, line: &mut Vec<u8>) -> io::Result<bool> {
    line.clear();
    let n = r.read_until(b'\n', line)?;
    if n == 0 {
        return Ok(false);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
        if line.last() == Some(&b'\r') {
            line.pop();
        }
    }
    Ok(true)
}

/// bcftools' `strip_header`: the new header with the old header's IDX numbering, so records
/// keep their IDs.
fn strip_header(src: &Header, dst: &Header) -> anyhow::Result<Header> {
    let is_id_kind = |k: HeaderRecordType| {
        matches!(
            k,
            HeaderRecordType::Filter
                | HeaderRecordType::Info
                | HeaderRecordType::Format
                | HeaderRecordType::Contig
        )
    };
    let mut out = Header::new();
    for r in dst.records() {
        if is_id_kind(r.kind()) {
            continue;
        }
        let _ = out.add_record(r.dup());
    }
    for r in src.records() {
        if !is_id_kind(r.kind()) {
            continue;
        }
        let Some(id) = r.get(b"ID") else {
            continue;
        };
        let Some(d) = dst.get_record(r.kind(), b"ID", Some(id), None) else {
            continue;
        };
        let mut tmp = d.dup();
        if let Some(idx) = r.get(b"IDX") {
            let j = pasta::variant::number::atoi(idx);
            tmp.push_field("IDX", j.to_string().as_bytes(), Quoting::Bare);
        }
        if out.add_record(tmp).is_err() {
            bail!("[strip_header] Failed to update header");
        }
    }
    for r in dst.records() {
        if !is_id_kind(r.kind()) {
            continue;
        }
        let Some(id) = r.get(b"ID") else {
            continue;
        };
        if out.get_record(r.kind(), b"ID", Some(id), None).is_none() {
            let _ = out.add_record(r.dup());
        }
    }
    for s in dst.samples() {
        let _ = out.add_sample(s);
    }
    Ok(out)
}

/// Whether an ID of the input header means the same in the output header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum IdCheck {
    Same,
    Undefined,
    Broken,
}

/// Per-ID checks of FILTER, INFO, and FORMAT IDs and contigs, computed once per header pair.
struct IdTable {
    contigs: Vec<bool>,
    kinds: [Vec<IdCheck>; 3],
}

impl IdTable {
    fn new(hdr: &Header, out: &Header) -> Self {
        let contigs = (0..hdr.contig_count() as i32)
            .map(|rid| {
                rid < out.contig_count() as i32 && hdr.contig_name(rid) == out.contig_name(rid)
            })
            .collect();
        let kind_table = |kind: HeaderRecordType, check_key: bool| -> Vec<IdCheck> {
            (0..hdr.id_count() as u32)
                .map(|id| {
                    if id as usize >= out.id_count() {
                        return IdCheck::Undefined;
                    }
                    if check_key && out.id_name(id).is_none() {
                        return IdCheck::Undefined;
                    }
                    if !out.id_exists(kind, id) {
                        return IdCheck::Undefined;
                    }
                    if hdr.id_name(id) != out.id_name(id) {
                        return IdCheck::Broken;
                    }
                    IdCheck::Same
                })
                .collect()
        };
        Self {
            contigs,
            kinds: [
                kind_table(HeaderRecordType::Filter, false),
                kind_table(HeaderRecordType::Info, true),
                kind_table(HeaderRecordType::Format, true),
            ],
        }
    }

    fn check(&self, slot: usize, id: i32) -> IdCheck {
        usize::try_from(id)
            .ok()
            .and_then(|i| self.kinds[slot].get(i).copied())
            .unwrap_or(IdCheck::Undefined)
    }
}

fn id_name(h: &Header, id: i32) -> String {
    u32::try_from(id)
        .ok()
        .and_then(|i| h.id_name(i))
        .map_or_else(
            || "(null)".to_owned(),
            |n| String::from_utf8_lossy(n).into_owned(),
        )
}

/// The output header of bcftools' `reheader_bcf`: the input header text, edited and parsed,
/// with the old IDX numbering when the header was replaced.
fn new_bcf_header(hdr: &Header, edits: &Edits) -> anyhow::Result<Header> {
    let mut htxt = Vec::new();
    hdr.format_into(&mut htxt, HeaderTextFlavor::Bcf);
    edits.apply(&mut htxt)?;
    let Ok(hdr_out) = Header::parse(&htxt) else {
        bail!("An error occurred while parsing the header");
    };
    if edits.header.is_some() {
        return strip_header(hdr, &hdr_out);
    }
    Ok(hdr_out)
}

/// bcftools' `reheader_bcf`: reads the records, checks their IDs against the new header, and
/// writes them unchanged. BGZF input keeps its compressed blocks where the records in them are
/// unchanged by reading.
fn reheader_bcf(
    input: RawInput,
    fname: &str,
    out: Option<&Path>,
    edits: &Edits,
    threads: usize,
) -> anyhow::Result<()> {
    if input.format != Format::Bcf {
        pasta::variant::log::error(
            "bcf_hdr_read",
            format_args!("Input is not detected as bcf or vcf format"),
        );
        bail!("Failed to read the header: {fname}");
    }
    if !output::threads_start(threads) {
        bail!("Could not initialize threading");
    }
    let name = out.map_or_else(|| "standard output".to_owned(), |p| p.display().to_string());
    let create = |is_compressed: bool| {
        Writer::create(
            out,
            &WriteOptions {
                output_type: if is_compressed {
                    OutputType::Bcf
                } else {
                    OutputType::BcfUncompressed
                },
                level: None,
                threads,
            },
        )
        .map_err(|e| {
            output::log_open_failure(&out_name(out), &e);
            anyhow::anyhow!("{}: {}", out_name(out), strerror(&e))
        })
    };
    let mut writer;
    let mut read_failed = false;
    let result = if input.compression == Compression::Bgzf {
        let mut copier = FrameCopier::new(InflatingReader::new(
            input.source,
            split_idx_name(fname).0,
            threads,
        ));
        let Some((hdr, header_end)) = copier.read_bcf_header() else {
            bail!("Failed to read the header: {fname}");
        };
        let hdr_out = new_bcf_header(&hdr, edits)?;
        writer = create(true)?;
        if writer.write_header(&hdr_out).is_err() {
            writer.abandon();
            bail!("[reheader_bcf] Error: cannot write the header to {name}");
        }
        let table = IdTable::new(&hdr, &hdr_out);
        let write_err = || anyhow::anyhow!("[reheader_bcf] Error: cannot write to {name}");
        let copied =
            copier.copy_bcf_records(header_end, &hdr, &hdr_out, &mut writer, &write_err, |rec| {
                check_record(&hdr, &hdr_out, &table, rec)?;
                if bcf::check_sample_count(&hdr_out, rec).is_err() {
                    return Err(write_err());
                }
                Ok(())
            });
        read_failed = copier.read_failed();
        copied
    } else {
        let is_compressed = input.compression != Compression::None;
        let mut src = input.into_decompressed(threads);
        let hdr = bcf::read_header(&mut src)
            .map_err(|_| anyhow::anyhow!("Failed to read the header: {fname}"))?;
        let hdr_out = new_bcf_header(&hdr, edits)?;
        writer = create(is_compressed)?;
        if writer.write_header(&hdr_out).is_err() {
            writer.abandon();
            bail!("[reheader_bcf] Error: cannot write the header to {name}");
        }
        let table = IdTable::new(&hdr, &hdr_out);
        let mut rec = Record::new();
        (|| -> anyhow::Result<()> {
            loop {
                match bcf::read_record(&mut src, &hdr, &mut rec) {
                    Ok(bcf::ReadOutcome::Record) => {}
                    Ok(bcf::ReadOutcome::Eof) => break,
                    // bcftools' `bcf_read` loop stops here as at the end, and closing the
                    // input reports nothing, so the records after a truncated or malformed one
                    // are dropped with status 0. vtools fails as it does for a truncated BGZF
                    // block, after writing the records before it (as `copy_bcf_records` does
                    // for BGZF input).
                    Err(_) => {
                        read_failed = true;
                        break;
                    }
                }
                check_record(&hdr, &hdr_out, &table, &mut rec)?;
                if let Err(e) = writer.write_record(&hdr_out, &mut rec) {
                    return Err(output::write_failure(e, || {
                        anyhow::anyhow!("[reheader_bcf] Error: cannot write to {name}")
                    }));
                }
            }
            Ok(())
        })()
    };
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    if let Err(e) = writer.finish() {
        return Err(output::write_failure(e, || {
            anyhow::anyhow!("[reheader_bcf] Error: failed to close the file {name}")
        }));
    }
    // The records read before a truncated or corrupt block are written; closing the input then
    // fails, as HTSlib keeps the read error (and, unlike bcftools, after a truncated or
    // malformed record).
    if read_failed {
        bail!("[reheader_bcf] Error: close failed .. {fname}");
    }
    Ok(())
}

/// The ID checks of bcftools' `reheader_bcf`.
fn check_record(
    hdr: &Header,
    hdr_out: &Header,
    table: &IdTable,
    rec: &mut Record,
) -> anyhow::Result<()> {
    rec.unpack(Unpack::ALL);
    let rid = rec.rid();
    let ok = usize::try_from(rid)
        .ok()
        .and_then(|r| table.contigs.get(r).copied())
        .unwrap_or(false);
    if !ok {
        bail!(
            "The CHROM is not defined: \"{}\"",
            hdr.contig_name(rid)
                .map_or_else(|| "(null)".into(), String::from_utf8_lossy)
        );
    }
    let check = |slot: usize, what: &str, id: i32| -> anyhow::Result<()> {
        match table.check(slot, id) {
            IdCheck::Same => Ok(()),
            IdCheck::Broken => bail!(
                "FIXME: Broken {what} ids: {} vs {}",
                id_name(hdr, id),
                id_name(hdr_out, id)
            ),
            IdCheck::Undefined => {
                let tag = if slot == 0 { "" } else { " tag" };
                bail!("The {what}{tag} is not defined: \"{}\"", id_name(hdr, id))
            }
        }
    };
    for &id in rec.filters() {
        check(0, "FILTER", id)?;
    }
    for f in rec.info_fields() {
        check(1, "INFO", f.key())?;
    }
    for f in rec.format_fields() {
        check(2, "FORMAT", f.id())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Read;

    use super::*;

    #[test]
    fn idx_names_read_the_data_file_of_a_header() {
        let (dir, name) =
            crate::raw_input::write_idx_named("reheader-idx", "h.vcf", "##fileformat=VCFv4.2\n\n");
        let hdr = read_header_file(&name);
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(
            hdr.expect("the data file is read"),
            b"##fileformat=VCFv4.2\n"
        );
    }

    const HDR: &[u8] = b"##fileformat=VCFv4.2\n##contig=<ID=a,length=5,assembly=\"x,y\">\n##contig=<ID=b,length=7>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

    fn fai(entries: &[(&str, i64)]) -> Fai {
        Fai {
            names: entries.iter().map(|(n, _)| n.as_bytes().to_vec()).collect(),
            lengths: entries
                .iter()
                .map(|&(n, l)| (n.as_bytes().to_vec(), l))
                .collect(),
        }
    }

    #[test]
    fn samples_are_replaced_or_renamed() {
        let mut h = HDR.to_vec();
        set_samples(&[b"X".to_vec(), b"Y".to_vec(), b"Z".to_vec()], &mut h).unwrap();
        assert!(h.ends_with(b"FORMAT\tX\tY\tZ\n"));

        let mut h = HDR.to_vec();
        set_samples(&[b"S2 new\\ name".to_vec()], &mut h).unwrap();
        assert!(h.ends_with(b"FORMAT\tS1\tnew name\n"));

        let mut h =
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_vec();
        let err = set_samples(&[b"X".to_vec()], &mut h).unwrap_err();
        assert!(err.to_string().contains("missing FORMAT"));

        let mut h =
            b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
                .to_vec();
        set_samples(&[b"X".to_vec(), b"Y".to_vec()], &mut h).unwrap();
        assert!(h.ends_with(b"\tINFO\tFORMAT\tX\tY\n"), "samples are added");
    }

    #[test]
    fn fai_lines_parse_as_htslib_reads_them() {
        let fai = Fai::parse("x.fai", b"chr1 100\t6\t60\t61\nchr2\t-5\t1\t2\t3\n")
            .expect("both lines parse");
        assert_eq!(fai.names, [b"chr1".to_vec(), b"chr2".to_vec()]);
        assert_eq!(fai.len(b"chr1"), Some(100), "the name ends at the space");
        assert_eq!(fai.len(b"chr2"), Some(-5), "unsigned values wrap");
        assert!(
            Fai::parse("x.fai", b"chr1 x\t1\t2\t3\t4\n").is_none(),
            "word after the name"
        );

        let fai = Fai::parse("x.fai", b"a\t1\t2\t3\t4\na\t9\t2\t3\t4").expect("parses");
        assert_eq!(fai.len(b"a"), Some(1), "the first of duplicates wins");

        assert!(
            Fai::parse("x.fai", b"a\t1\t2\t3\t4\n\n").is_none(),
            "empty line"
        );
        assert!(
            Fai::parse("x.fai", b"a\t1\t2\t3\n").is_none(),
            "four columns"
        );
        assert!(
            Fai::parse("x.fai", b"a\t1\t2\t61\t60\n").is_none(),
            "width < bases"
        );
        assert!(Fai::parse("x.fai", b"").is_some_and(|f| f.names.is_empty()));
    }

    #[test]
    fn scan_u64_follows_strtoull() {
        assert_eq!(scan_u64(b"  +12x", 0), Some((12, 5)));
        assert_eq!(scan_u64(b"-1", 0), Some((u64::MAX, 2)));
        assert_eq!(
            scan_u64(b"99999999999999999999999", 0),
            Some((u64::MAX, 23))
        );
        assert_eq!(scan_u64(b" x", 0), None);
    }

    #[test]
    fn contigs_follow_the_fai() {
        let out = update_from_fai(&fai(&[("b", 70), ("c", 9)]), HDR).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(!text.contains("ID=a"), "{text}");
        assert!(text.contains("##contig=<ID=b,length=70>\n"), "{text}");
        assert!(text.contains("##contig=<ID=c,length=9>\n#CHROM"), "{text}");

        let out = update_from_fai(&fai(&[("a", 50)]), HDR).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.contains("##contig=<ID=a,length=50,assembly=\"x,y\">\n"),
            "{text}"
        );
    }

    /// A scratch directory for the files of one test, removed when dropped.
    struct Scratch(PathBuf);

    impl Scratch {
        fn new(name: &str) -> Self {
            let dir =
                std::env::temp_dir().join(format!("vtools-reheader-{name}-{}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("create the test directory");
            Self(dir)
        }

        /// Writes `data` to `name` in the directory and returns its path.
        fn file(&self, name: &str, data: &[u8]) -> String {
            let path = self.0.join(name);
            std::fs::write(&path, data).expect("write the test file");
            path.to_str().expect("UTF-8 path").to_owned()
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    const RECORDS: [&str; 3] = [
        "a\t1\t.\tA\tG\t.\t.\tDP=1\tGT\t0/1\t1/1",
        "a\t2\t.\tC\tT\t.\t.\tDP=2\tGT\t0/0\t0/1",
        "b\t3\t.\tG\tA\t.\t.\tDP=3\tGT\t1/1\t0/0",
    ];

    /// Renames the samples to X and Y.
    fn rename() -> Edits {
        Edits {
            samples: Some(SampleSpec {
                spec: "X,Y".to_owned(),
                source: ListSource::Inline,
            }),
            header: None,
            fai: None,
        }
    }

    /// [`HDR`] with GT.
    const GT_HDR: &[u8] = b"##fileformat=VCFv4.2\n##contig=<ID=a,length=5>\n##contig=<ID=b,length=7>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

    /// Returns [`GT_HDR`] and [`RECORDS`] as uncompressed (not BGZF) BCF, and where each
    /// record starts.
    fn raw_bcf() -> (Vec<u8>, Vec<usize>) {
        let mut hdr = Header::parse(GT_HDR).expect("the header parses");
        let mut out = bcf::encode_header(&hdr);
        let mut starts = Vec::new();
        let mut rec = Record::new();
        for r in RECORDS {
            starts.push(out.len());
            let mut line = r.as_bytes().to_vec();
            pasta::variant::vcf::parse(&mut hdr, &mut line, &mut rec).expect("the line parses");
            bcf::encode_record(&hdr, &mut rec, &mut out).expect("the record encodes");
        }
        (out, starts)
    }

    /// The samples and record positions of a BGZF BCF file.
    fn read_bcf(path: &str) -> (Vec<Box<[u8]>>, Vec<i64>) {
        let file = File::open(path).expect("open the output");
        let mut src = io::BufReader::new(pasta::bgzf::io::Reader::new(file));
        let hdr = bcf::read_header(&mut src).expect("the output header reads");
        let mut rec = Record::new();
        let mut pos = Vec::new();
        while let bcf::ReadOutcome::Record =
            bcf::read_record(&mut src, &hdr, &mut rec).expect("the output records read")
        {
            pos.push(rec.pos() + 1);
        }
        (hdr.samples().to_vec(), pos)
    }

    fn reheader_bcf_file(dir: &Scratch, input: &str) -> (anyhow::Result<()>, String) {
        let out = dir.0.join("out.bcf");
        let input_raw = RawInput::open(input).expect("open the input");
        let r = reheader_bcf(input_raw, input, Some(&out), &rename(), 0);
        (r, out.to_str().expect("UTF-8 path").to_owned())
    }

    #[test]
    fn uncompressed_bcf_is_copied_whole() {
        let dir = Scratch::new("raw-ok");
        let (data, _) = raw_bcf();
        let input = dir.file("in.bcf", &data);
        let (r, out) = reheader_bcf_file(&dir, &input);
        r.expect("reheader succeeds");
        assert_eq!(
            read_bcf(&out),
            (
                vec![Box::from(&b"X"[..]), Box::from(&b"Y"[..])],
                vec![1, 2, 3]
            )
        );
    }

    #[test]
    fn truncated_uncompressed_bcf_fails_after_the_records_before() {
        let dir = Scratch::new("raw-truncated");
        let (data, _) = raw_bcf();
        let input = dir.file("in.bcf", &data[..data.len() - 5]);
        let (r, out) = reheader_bcf_file(&dir, &input);
        let e = r.expect_err("a truncated record fails");
        assert_eq!(
            e.to_string(),
            format!("[reheader_bcf] Error: close failed .. {input}")
        );
        assert_eq!(read_bcf(&out).1, [1, 2]);
    }

    #[test]
    fn malformed_uncompressed_bcf_record_fails() {
        let dir = Scratch::new("raw-malformed");
        let (mut data, starts) = raw_bcf();
        // The second record's CHROM becomes contig 9, which the header lacks.
        data[starts[1] + 8..starts[1] + 12].copy_from_slice(&9i32.to_le_bytes());
        let input = dir.file("in.bcf", &data);
        let (r, out) = reheader_bcf_file(&dir, &input);
        assert!(r.is_err(), "a malformed record fails");
        assert_eq!(read_bcf(&out).1, [1]);
    }

    #[test]
    fn truncated_or_malformed_bgzf_bcf_fails_after_the_records_before() {
        use std::io::Write;

        let dir = Scratch::new("bgzf-bad-record");
        let (data, starts) = raw_bcf();
        let mut malformed = data.clone();
        // The second record's shared block runs past the data.
        malformed[starts[1]..starts[1] + 4].copy_from_slice(&1000u32.to_le_bytes());
        for (name, data, want) in [
            ("truncated.bcf", &data[..data.len() - 5], &[1, 2][..]),
            ("malformed.bcf", &malformed[..], &[1][..]),
        ] {
            let mut w = pasta::bgzf::io::Writer::new(Vec::new());
            w.write_all(data).expect("writing to memory works");
            let input = dir.file(name, &w.finish().expect("finishing in memory works"));
            let (r, out) = reheader_bcf_file(&dir, &input);
            let e = r.expect_err("a truncated or malformed record fails");
            assert_eq!(
                e.to_string(),
                format!("[reheader_bcf] Error: close failed .. {input}")
            );
            assert_eq!(read_bcf(&out).1, want, "{name}");
        }
    }

    /// Returns [`GT_HDR`] and [`RECORDS`] as bgzipped VCF, the header and each record in a
    /// block of its own, with or without the end-of-file marker.
    fn vcf_gz(eof: bool) -> Vec<u8> {
        use pasta::bgzf::io::{raw, writer::CompressionLevel};

        let mut out = Vec::new();
        let mut frame = Vec::new();
        raw::deflate_frame(GT_HDR, CompressionLevel::default(), &mut frame).expect("deflate");
        out.extend_from_slice(&frame);
        for r in RECORDS {
            let line = format!("{r}\n");
            raw::deflate_frame(line.as_bytes(), CompressionLevel::default(), &mut frame)
                .expect("deflate");
            out.extend_from_slice(&frame);
        }
        if eof {
            out.extend_from_slice(&raw::EOF_MARKER);
        }
        out
    }

    /// Runs reheader on bgzipped VCF `data`, returning the result and the decompressed output.
    fn reheader_gz(name: &str, data: &[u8]) -> (anyhow::Result<()>, String, String) {
        let dir = Scratch::new(name);
        let input = dir.file("in.vcf.gz", data);
        let out = dir.0.join("out.vcf.gz");
        let r = reheader_vcf_gz(
            RawInput::open(&input).expect("open the input"),
            &input,
            Some(&out),
            &rename(),
        );
        let mut text = String::new();
        pasta::bgzf::io::Reader::new(File::open(&out).expect("open the output"))
            .read_to_string(&mut text)
            .expect("the output is complete BGZF");
        (r, text, input)
    }

    #[test]
    fn bgzipped_vcf_record_blocks_are_copied() {
        let (r, text, _) = reheader_gz("gz-ok", &vcf_gz(true));
        r.expect("reheader succeeds");
        assert!(text.contains("\tFORMAT\tX\tY\n"), "{text}");
        assert_eq!(text.lines().filter(|l| !l.starts_with('#')).count(), 3);
        // Only a missing end-of-file marker (a warning).
        let (r, text, _) = reheader_gz("gz-no-eof", &vcf_gz(false));
        r.expect("reheader succeeds");
        assert_eq!(text.lines().filter(|l| !l.starts_with('#')).count(), 3);
    }

    #[test]
    fn truncated_bgzipped_vcf_fails_without_copying_the_cut_block() {
        let data = vcf_gz(true);
        let cut = &data[..data.len() - pasta::bgzf::io::raw::EOF_MARKER.len() - 5];
        let (r, text, input) = reheader_gz("gz-truncated", cut);
        let e = r.expect_err("a truncated block fails");
        assert_eq!(e.to_string(), format!("Error reading {input}"));
        assert!(text.ends_with(&format!("{}\n", RECORDS[1])), "{text}");
    }
}
