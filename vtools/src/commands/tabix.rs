//! `vtools tabix`: bcftools' small tabix subcommand: index a bgzipped tab-delimited file
//! (GFF, BED, SAM, VCF, or custom columns), print it, or print the lines overlapping regions.

use std::{
    fs::File,
    io::{self, BufRead, Read, Write},
    path::Path,
    process::ExitCode,
};

use pasta::{
    bgzf::{self, VirtualPosition},
    variant::{
        index::{
            IndexFormat,
            tabix::{PRESET_GENERIC, PRESET_SAM, PRESET_UCSC, Tabix, TabixBuilder, TabixConfig},
        },
        log,
        region::{RegionFlags, parse_region},
    },
};

use crate::{
    commands::index::warn_if_no_eof,
    raw_input::{Compression, RawInput},
};

const USAGE: &str = "
Usage: vtools tabix [options] <in.gz> [reg1 [...]]

Options: -p STR    preset: gff, bed, sam or vcf [gff]
         -s INT    column number for sequence names (suppressed by -p) [1]
         -b INT    column number for region start [4]
         -e INT    column number for region end (if no end, set INT to -b) [5]
         -0        specify coordinates are zero-based
         -S INT    skip first INT lines [0]
         -c CHAR   skip lines starting with CHAR [null]
         -a        print all records
         -f        force to overwrite existing index
         -m INT    set the minimal interval size to 1<<INT; 0 for the old tabix index [0]

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// The bgzipped file, then regions to print.
    inputs: Vec<String>,

    /// Coordinates are zero-based.
    #[arg(short = '0')]
    zero_based: bool,

    /// Overwrite an existing index.
    #[arg(short = 'f')]
    force: bool,

    /// Print all records.
    #[arg(short = 'a')]
    all: bool,

    /// Preset: gff, bed, sam, or vcf [gff].
    #[arg(short = 'p', value_name = "STR", allow_hyphen_values = true)]
    preset: Option<String>,

    /// Column of sequence names [1].
    #[arg(short = 's', value_name = "INT", allow_hyphen_values = true)]
    seq_col: Option<String>,

    /// Column of region starts [4].
    #[arg(short = 'b', value_name = "INT", allow_hyphen_values = true)]
    begin_col: Option<String>,

    /// Column of region ends [5].
    #[arg(short = 'e', value_name = "INT", allow_hyphen_values = true)]
    end_col: Option<String>,

    /// Skip the first INT lines.
    #[arg(short = 'S', value_name = "INT", allow_hyphen_values = true)]
    skip_lines: Option<String>,

    /// Skip lines starting with CHAR.
    #[arg(short = 'c', value_name = "CHAR", allow_hyphen_values = true)]
    meta_char: Option<String>,

    /// Minimal interval size 1<<INT; 0 for a TBI index [0].
    #[arg(short = 'm', value_name = "INT", allow_hyphen_values = true)]
    min_shift: Option<String>,
}

/// HTSlib's `tbx_conf_gff`.
const GFF: TabixConfig = TabixConfig {
    preset: PRESET_GENERIC,
    sc: 1,
    bc: 4,
    ec: 5,
    meta_char: b'#' as i32,
    line_skip: 0,
};

/// HTSlib's `tbx_conf_sam`.
const SAM: TabixConfig = TabixConfig {
    preset: PRESET_SAM,
    sc: 3,
    bc: 4,
    ec: 0,
    meta_char: b'@' as i32,
    line_skip: 0,
};

/// The size of the output buffer of `-a` and of region queries.
const OUT_BUFFER: usize = 1 << 17;

fn atoi(s: &str) -> i32 {
    pasta::variant::number::atoi(s.as_bytes()) as i32
}

pub fn run(args: Args, matches: &clap::ArgMatches) -> anyhow::Result<ExitCode> {
    let mut conf = GFF;
    let mut detect = true;
    // getopt applies the options in command-line order: -p resets the columns.
    let mut opts: Vec<(usize, char, &str)> = Vec::new();
    for (id, c, value) in [
        ("preset", 'p', &args.preset),
        ("seq_col", 's', &args.seq_col),
        ("begin_col", 'b', &args.begin_col),
        ("end_col", 'e', &args.end_col),
        ("skip_lines", 'S', &args.skip_lines),
        ("meta_char", 'c', &args.meta_char),
    ] {
        if let (Some(v), Some(i)) = (value, matches.index_of(id)) {
            opts.push((i, c, v));
        }
    }
    if args.zero_based
        && let Some(i) = matches.index_of("zero_based")
    {
        opts.push((i, '0', ""));
    }
    opts.sort_by_key(|o| o.0);
    for (_, c, v) in opts {
        match c {
            '0' => conf.preset |= PRESET_UCSC,
            's' => conf.sc = atoi(v),
            'b' => conf.bc = atoi(v),
            'e' => conf.ec = atoi(v),
            'c' => conf.meta_char = v.as_bytes().first().map_or(0, |&b| i32::from(b)),
            'S' => conf.line_skip = atoi(v),
            _ => {
                conf = match v {
                    "gff" => GFF,
                    "bed" => TabixConfig::BED,
                    "sam" => SAM,
                    "vcf" => TabixConfig::VCF,
                    _ => {
                        eprintln!("The type '{v}' not recognised");
                        return Ok(ExitCode::from(1));
                    }
                };
                detect = false;
            }
        }
    }
    let min_shift = args.min_shift.as_deref().map_or(-1, atoi);
    let Some(fname) = args.inputs.first() else {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    };
    if args.all {
        return Ok(print_all(fname));
    }
    if args.inputs.len() < 2 {
        if detect {
            let lower = fname.to_ascii_lowercase();
            if lower.ends_with(".gff.gz") {
                conf = GFF;
            } else if lower.ends_with(".bed.gz") {
                conf = TabixConfig::BED;
            } else if lower.ends_with(".sam.gz") {
                conf = SAM;
            } else if lower.ends_with(".vcf.gz") {
                conf = TabixConfig::VCF;
            }
        }
        if !args.force {
            let idx = format!("{fname}{}", if min_shift <= 0 { ".tbi" } else { ".csi" });
            if File::open(&idx).is_ok() {
                eprintln!("[E::main_tabix] the index file exists; use option '-f' to overwrite");
                return Ok(ExitCode::from(1));
            }
        }
        if build_index(fname, min_shift, conf).is_err() {
            eprintln!(
                "tbx_index_build failed: Is the file bgzip-compressed? Was wrong -p [type] option used?"
            );
            return Ok(ExitCode::from(1));
        }
        return Ok(ExitCode::SUCCESS);
    }
    Ok(query(fname, &args.inputs[1..]))
}

/// Prints every line (`-a`, reading without the index).
fn print_all(fname: &str) -> ExitCode {
    let Ok(input) = RawInput::open(fname) else {
        eprintln!("[E::main_tabix] could not open '{fname}'");
        return ExitCode::from(1);
    };
    let is_compressed = input.compression != Compression::None;
    let is_bgzf = input.compression == Compression::Bgzf;
    let mut src = input.into_decompressed(0);
    let mut out = io::BufWriter::with_capacity(OUT_BUFFER, crate::output::stdout());
    let mut line = Vec::new();
    loop {
        line.clear();
        match src.read_until(b'\n', &mut line) {
            Ok(0) => {
                if is_bgzf {
                    warn_if_no_eof(fname);
                }
                break;
            }
            Err(e) => {
                if is_compressed {
                    report_block_error(&e);
                }
                break;
            }
            Ok(_) => {}
        }
        if line.last() == Some(&b'\n') {
            line.pop();
            if line.last() == Some(&b'\r') {
                line.pop();
            }
        } else if is_compressed && line.last() == Some(&b'\r') {
            line.pop();
        }
        // `puts` stops at a NUL.
        let end = line.iter().position(|&c| c == 0).unwrap_or(line.len());
        if out.write_all(&line[..end]).is_err() || out.write_all(b"\n").is_err() {
            return ExitCode::SUCCESS;
        }
    }
    let _ = out.flush();
    ExitCode::SUCCESS
}

/// Reports a failed read of a compressed block, as HTSlib's `bgzf_read_block` does.
fn report_block_error(e: &io::Error) {
    log::error(
        "bgzf_read_block",
        format_args!(
            "Failed to read BGZF block data: {}",
            crate::raw_input::strerror(e)
        ),
    );
}

/// Reads lines from a BGZF file with their virtual offsets (HTSlib's `bgzf_getline` and
/// `bgzf_tell`).
struct Lines<R = File> {
    reader: bgzf::io::Reader<R>,
    fname: String,
    /// The end of the file was reached (and reported when the EOF marker is absent).
    at_end: bool,
}

impl Lines {
    fn open(fname: &str) -> io::Result<Self> {
        Ok(Lines::new(File::open(fname)?, fname))
    }
}

impl<R: Read> Lines<R> {
    fn new(inner: R, fname: &str) -> Self {
        Self {
            reader: bgzf::io::Reader::new(inner),
            fname: fname.to_owned(),
            at_end: false,
        }
    }

    /// Reads one line without its newline (or a `\r` before it); `false` at the end.
    fn next(&mut self, line: &mut Vec<u8>) -> io::Result<bool> {
        line.clear();
        let mut found = false;
        loop {
            let buf = self.reader.fill_buf().inspect_err(report_block_error)?;
            if buf.is_empty() {
                break;
            }
            match memchr::memchr(b'\n', buf) {
                Some(i) => {
                    line.extend_from_slice(&buf[..i]);
                    self.reader.consume(i + 1);
                    found = true;
                    break;
                }
                None => {
                    let n = buf.len();
                    line.extend_from_slice(buf);
                    self.reader.consume(n);
                }
            }
        }
        if !found && line.is_empty() {
            if !self.at_end {
                self.at_end = true;
                warn_if_no_eof(&self.fname);
            }
            return Ok(false);
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        Ok(true)
    }

    fn tell(&self) -> u64 {
        u64::from(self.reader.virtual_position())
    }
}

/// HTSlib's `tbx_index_build`: indexes `fname` into `fname.tbi`, or `fname.csi` with a
/// positive minimum shift.
fn build_index(fname: &str, min_shift: i32, conf: TabixConfig) -> Result<(), ()> {
    let Ok(input) = RawInput::open(fname) else {
        return Err(());
    };
    // HTSlib 1.24's `tbx_index_build3` fails silently on input that is not BGZF.
    if input.compression != Compression::Bgzf {
        return Err(());
    }
    // The input is read on from what detection read: a pipe or FIFO cannot be opened again.
    let mut lines = Lines::new(input.source, fname);
    let mut builder = TabixBuilder::new(conf, min_shift.max(0) as u32);
    let mut line = Vec::new();
    while lines.next(&mut line).map_err(|_| ())? {
        builder.push_line(&line, lines.tell()).map_err(|_| ())?;
    }
    let fmt = builder.format();
    let tbx = builder.finish(lines.tell());
    let suffix = if fmt == IndexFormat::Csi {
        "csi"
    } else {
        "tbi"
    };
    tbx.index
        .save(Path::new(&format!("{fname}.{suffix}")), fmt)
        .map_err(|_| ())
}

/// Prints the lines overlapping each region (HTSlib's `tbx_itr_querys` and
/// `tbx_bgzf_itr_next`).
fn query(fname: &str, regions: &[String]) -> ExitCode {
    let Some(idx) = pasta::variant::index::locate(fname, IndexFormat::Tbi) else {
        log::error(
            "idx_find_and_load",
            format_args!("Could not retrieve index file for '{fname}'"),
        );
        return ExitCode::from(1);
    };
    pasta::variant::index::warn_if_older(Path::new(fname), &idx);
    let mut tbx = match Tabix::load(&idx) {
        Ok(t) => t,
        Err(e) => {
            pasta::variant::index::report_load_error(&idx, &e);
            return ExitCode::from(1);
        }
    };
    let Ok(mut lines) = Lines::open(fname) else {
        return ExitCode::from(1);
    };
    let mut out = io::BufWriter::with_capacity(OUT_BUFFER, crate::output::stdout());
    let mut line = Vec::new();
    for reg in regions {
        let plan = if reg == "." {
            let Some(off) = tbx.index.first_offset() else {
                continue;
            };
            Plan::All(off)
        } else if reg == "*" {
            continue;
        } else {
            let Some(r) = parse_region(
                reg.as_bytes(),
                |name| tbx.name_to_id(name),
                RegionFlags::default(),
            ) else {
                continue;
            };
            let Some(p) = tbx.index.query(r.tid as usize, r.beg, r.end) else {
                continue;
            };
            Plan::Query(p)
        };
        if print_query(&mut lines, &mut tbx, &plan, &mut line, &mut out).is_err() {
            break;
        }
    }
    let _ = out.flush();
    ExitCode::SUCCESS
}

enum Plan {
    /// Everything from this virtual offset (`HTS_IDX_START`).
    All(u64),
    Query(pasta::variant::index::QueryPlan),
}

fn print_query(
    lines: &mut Lines,
    tbx: &mut Tabix,
    plan: &Plan,
    line: &mut Vec<u8>,
    out: &mut impl Write,
) -> io::Result<()> {
    let mut emit = |line: &[u8]| -> io::Result<()> {
        let end = line.iter().position(|&c| c == 0).unwrap_or(line.len());
        out.write_all(&line[..end])?;
        out.write_all(b"\n")
    };
    match plan {
        Plan::All(off) => {
            lines.reader.seek(VirtualPosition::from(*off))?;
            while lines.next(line)? {
                if tbx.line_interval(line).is_err() {
                    break;
                }
                emit(line)?;
            }
        }
        Plan::Query(p) => {
            let mut prev_end: Option<u64> = None;
            'chunks: for &(beg, end) in &p.chunks {
                if prev_end != Some(beg) {
                    lines.reader.seek(VirtualPosition::from(beg))?;
                }
                while lines.tell() < end {
                    if !lines.next(line)? {
                        break 'chunks;
                    }
                    let Ok((tid, intv)) = tbx.line_interval(line) else {
                        break 'chunks;
                    };
                    if p.is_past(i64::from(tid), intv.beg) {
                        break 'chunks;
                    }
                    if p.overlaps(intv.beg, intv.end) {
                        emit(line)?;
                    }
                }
                prev_end = Some(end);
            }
        }
    }
    Ok(())
}
