//! `vtools index`: build CSI/TBI indexes for bgzipped VCF and BCF, or report the record counts
//! they hold (bcftools index).

use std::{
    io::{self, IsTerminal, Write},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Record,
    index::{
        Index, IndexFormat, adjust_csi_settings, locate, report_load_error,
        tabix::{Tabix, TabixBuilder, TabixConfig},
        warn_if_older,
    },
    io::{
        Compression, InputFormat, ReadError, Reader, has_bgzf_eof,
        reader::{report_missing_eof, take_missing_eof},
        stdin_has_bgzf_eof,
    },
};

use super::head::{apply_verbosity, strtol_whole};
use crate::options::last_index;
use crate::raw_input::report_exec_format_error;

const USAGE: &str = "
About:   Index bgzip compressed VCF/BCF files for random access.
Usage:   vtools index [options] <in.bcf>|<in.vcf.gz>

Indexing options:
    -c, --csi                generate CSI-format index for VCF/BCF files [default]
    -f, --force              overwrite index if it already exists
    -m, --min-shift INT      set minimal interval size for CSI indices to 2^INT [14]
    -o, --output FILE        optional output index file name
    -t, --tbi                generate TBI-format index for VCF files
        --threads INT        use multithreading with INT worker threads [0]
    -v, --verbosity INT      verbosity level

Stats options:
    -a, --all            with --stats, print stats for all contigs even when zero
    -n, --nrecords       print number of records based on existing index file
    -s, --stats          print per contig stats based on existing index file

";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input BCF or bgzipped VCF (`-` for stdin, which needs -o).
    inputs: Vec<String>,

    /// Generate a CSI index [default].
    #[arg(short = 'c', long = "csi")]
    csi: bool,

    /// Generate a TBI index (VCF only).
    #[arg(short = 't', long = "tbi")]
    tbi: bool,

    /// Overwrite an existing index.
    #[arg(short = 'f', long = "force")]
    force: bool,

    /// Minimal interval size for CSI indexes is 2^INT [14].
    #[arg(
        short = 'm',
        long = "min-shift",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    min_shift: Option<String>,

    /// Print per-contig record counts from the index.
    #[arg(short = 's', long = "stats")]
    stats: bool,

    /// Print the number of records from the index.
    #[arg(short = 'n', long = "nrecords")]
    nrecords: bool,

    /// With --stats, list contigs without records too.
    #[arg(short = 'a', long = "all")]
    all: bool,

    /// Worker threads for decompression.
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

    /// Index file name [FILE.csi or FILE.tbi].
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "output-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<String>,

    /// Print help.
    #[arg(long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,
}

const PER_CONTIG: u32 = 1;
const ALL_CONTIGS: u32 = 2;
const TOTAL: u32 = 4;

/// The default `--min-shift` of a CSI index: bins of 2^14 bases at the finest level.
const CSI_MIN_SHIFT: i64 = 14;

/// The length of a BGZF block header, which [`check_bgzf`] reads to recognize BGZF.
const BGZF_HEADER_LEN: usize = 18;

/// C's `strtol(s, &end, 10)` into an `int`, with nothing left after the number (an empty text
/// is 0).
fn strtol_strict(s: &str) -> Option<i64> {
    strtol_whole(s).map(|v| i64::from(v as i32))
}

pub fn run(args: Args, matches: &clap::ArgMatches) -> anyhow::Result<ExitCode> {
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    // In command-line order: -c clears the TBI choice, -t sets it with a min-shift of 0, and
    // -m sets the min-shift.
    let last = |id: &str| last_index(matches, id);
    let mut tbi = args.tbi && (!args.csi || last("tbi") > last("csi"));
    let mut min_shift = CSI_MIN_SHIFT;
    if let Some(m) = &args.min_shift {
        let Some(v) = strtol_strict(m) else {
            bail!("Could not parse argument: --min-shift {m}");
        };
        min_shift = v;
    }
    if args.tbi && (args.min_shift.is_none() || last("tbi") > last("min_shift")) {
        min_shift = 0;
    }
    // A negative thread count leaves the build unthreaded.
    let threads = match &args.threads {
        Some(t) => match strtol_strict(t) {
            Some(n) => usize::try_from(crate::output::cap_threads(n)).unwrap_or(0),
            None => bail!("Could not parse argument: --threads {t}"),
        },
        None => 0,
    };
    let mut stats = 0;
    if args.stats {
        stats |= PER_CONTIG;
    }
    if args.nrecords {
        stats |= TOTAL;
    }
    if args.all {
        stats |= ALL_CONTIGS;
    }
    if min_shift == 0 {
        tbi = true;
    }
    if stats > TOTAL {
        eprintln!("[E::main_vcfindex] expected only one of --stats or --nrecords options");
        return Ok(ExitCode::from(1));
    }
    if tbi && min_shift > 0 {
        eprintln!("[E::main_vcfindex] min-shift option only expected for CSI indices ");
        return Ok(ExitCode::from(1));
    }
    if !(0..=30).contains(&min_shift) {
        eprintln!("[E::main_vcfindex] expected min_shift in range [0,30] ({min_shift})");
        return Ok(ExitCode::from(1));
    }
    let fname = match args.inputs.first() {
        Some(f) => f.clone(),
        None if !io::stdin().is_terminal() => "-".to_owned(),
        None => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
    };
    if stats != 0 {
        return Ok(ExitCode::from(index_stats(&fname, stats)?));
    }

    let idx_fname = match &args.output {
        Some(o) => PathBuf::from(o),
        None => {
            if fname == "-" {
                eprintln!(
                    "[E::main_vcfindex] must specify an output path for index file when reading VCF/BCF from stdin"
                );
                return Ok(ExitCode::from(1));
            }
            PathBuf::from(format!("{fname}.{}", if tbi { "tbi" } else { "csi" }))
        }
    };
    if !args.force {
        if let Ok(idx_meta) = std::fs::metadata(&idx_fname) {
            let newer = std::fs::metadata(&fname)
                .and_then(|m| m.modified())
                .ok()
                .zip(idx_meta.modified().ok())
                .is_some_and(|(data, idx)| mtime_secs(data) > mtime_secs(idx));
            if !newer {
                eprintln!(
                    "[E::main_vcfindex] the index file exists. Please use '-f' to overwrite {}",
                    idx_fname.display()
                );
                return Ok(ExitCode::from(1));
            }
        }
        check_bgzf(&fname)?;
    }
    match build_index(&fname, &idx_fname, min_shift as u32, threads) {
        Ok(()) => {
            // `-f` indexes a file without the EOF marker block on purpose ("use -f to index
            // anyway"): HTSlib's warning is all, without the failure other commands report.
            let _ = take_missing_eof();
            Ok(ExitCode::SUCCESS)
        }
        Err(BuildError::Open) => bail!("index: failed to open \"{fname}\""),
        Err(BuildError::Format) => {
            bail!("index: \"{fname}\" is in a format that cannot be usefully indexed")
        }
        Err(BuildError::Failed) => bail!("index: failed to create index for \"{fname}\""),
    }
}

fn mtime_secs(t: std::time::SystemTime) -> u64 {
    t.duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
}

/// HTSlib's warning when reading reaches the end of a BGZF file whose last block is not the
/// EOF marker (`bgzf_read_block`), and the failure vtools reports for it when the command is
/// done ([`report_missing_eof`]). Standard input is checked when it is a regular file (a pipe
/// is checked as it is read).
pub fn warn_if_no_eof(fname: &str) {
    let has_eof = if fname == "-" {
        stdin_has_bgzf_eof()
    } else {
        has_bgzf_eof(&pasta::variant::io::reader::data_path(fname))
    };
    if let Ok(Some(false)) = has_eof {
        report_missing_eof(fname);
    }
}

/// The checks bcftools makes before indexing without `-f`.
fn check_bgzf(fname: &str) -> anyhow::Result<()> {
    let mut head = [0u8; BGZF_HEADER_LEN];
    let (n, has_eof) = if fname == "-" {
        stdin_head(&mut head)?
    } else {
        let mut f = match std::fs::File::open(fname) {
            Ok(f) => f,
            Err(_) => bail!("index: failed to open {fname}"),
        };
        let n = read_full(&mut f, &mut head);
        (n, has_bgzf_eof(Path::new(fname))?)
    };
    let is_bgzf = n >= BGZF_HEADER_LEN
        && head[0] == 0x1f
        && head[1] == 0x8b
        && head[3] & 4 != 0
        && &head[12..16] == b"BC\x02\x00";
    if !is_bgzf {
        bail!("index: the file is not BGZF compressed, cannot index: {fname}");
    }
    // HTSlib cannot check a pipe (`bgzf_check_EOF` returns 2), which counts as truncated.
    if has_eof != Some(true) {
        bail!("index: the input is probably truncated, use -f to index anyway: {fname}");
    }
    Ok(())
}

/// Reads up to `buf.len()` bytes; returns how many were read.
fn read_full(r: &mut impl std::io::Read, buf: &mut [u8]) -> usize {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) | Err(_) => break,
            Ok(k) => n += k,
        }
    }
    n
}

/// The first bytes of standard input and whether it ends with the BGZF EOF marker. A regular
/// file is read in place, leaving its offset for the indexing; a pipe is consumed, but then
/// the checks fail anyway, as they do in bcftools.
fn stdin_head(head: &mut [u8; BGZF_HEADER_LEN]) -> anyhow::Result<(usize, Option<bool>)> {
    #[cfg(unix)]
    {
        use std::os::{fd::AsFd, unix::fs::FileExt};

        let f = std::fs::File::from(io::stdin().as_fd().try_clone_to_owned()?);
        if f.metadata()?.is_file() {
            let pos = std::io::Seek::stream_position(&mut &f)?;
            let mut n = 0;
            while n < head.len() {
                match f.read_at(&mut head[n..], pos + n as u64) {
                    Ok(0) | Err(_) => break,
                    Ok(k) => n += k,
                }
            }
            return Ok((n, stdin_has_bgzf_eof()?));
        }
    }
    Ok((read_full(&mut io::stdin().lock(), head), None))
}

enum BuildError {
    Open,
    Format,
    Failed,
}

/// Builds and saves an index (HTSlib's `bcf_index_build3`).
fn build_index(
    fname: &str,
    idx_fname: &Path,
    min_shift: u32,
    threads: usize,
) -> Result<(), BuildError> {
    let mut reader = match Reader::open(fname, threads) {
        Ok(r) => r,
        Err(ReadError::Io(_)) => return Err(BuildError::Open),
        Err(ReadError::UnknownBinary) => {
            report_exec_format_error(fname);
            return Err(BuildError::Open);
        }
        Err(_) => return Err(BuildError::Format),
    };
    if reader.compression() != Compression::Bgzf {
        return Err(BuildError::Format);
    }
    match reader.format() {
        InputFormat::Bcf => {
            if min_shift == 0 {
                pasta::variant::log::error(
                    "bcf_index_build3",
                    format_args!("TBI indices for BCF files are not supported"),
                );
                return Err(BuildError::Failed);
            }
            let idx = bcf_index(&mut reader, min_shift).ok_or(BuildError::Failed)?;
            warn_if_no_eof(fname);
            idx.save(idx_fname, IndexFormat::Csi)
                .map_err(|_| BuildError::Failed)
        }
        InputFormat::Vcf => {
            let tbx = tabix_index(&mut reader, min_shift).ok_or(BuildError::Failed)?;
            warn_if_no_eof(fname);
            let fmt = if min_shift > 0 {
                IndexFormat::Csi
            } else {
                IndexFormat::Tbi
            };
            tbx.index
                .save(idx_fname, fmt)
                .map_err(|_| BuildError::Failed)
        }
    }
}

/// Indexes a BCF file (HTSlib's `bcf_index`).
pub fn bcf_index(reader: &mut Reader, min_shift: u32) -> Option<Index> {
    let header = reader.read_header().ok()?;
    let (n_lvls, min_shift, nids) = csi_levels(&header, min_shift, 0);
    let mut idx = Index::new(
        nids,
        IndexFormat::Csi,
        reader.virtual_position(),
        min_shift,
        n_lvls,
    );
    let mut rec = Record::new();
    let mut header = header;
    let mut line = Vec::new();
    loop {
        match reader.read_record(&mut header, &mut line, &mut rec) {
            Ok(true) => {}
            Ok(false) => break,
            Err(_) => return None,
        }
        idx.push(
            rec.rid(),
            rec.pos(),
            rec.pos() + rec.rlen(),
            reader.virtual_position(),
        )
        .ok()?;
    }
    idx.finish(reader.virtual_position());
    Some(idx)
}

/// The CSI levels for a header's longest contig (HTSlib's `idx_calc_n_lvls_ids`), with the
/// adjusted minimum shift and the number of contigs.
pub fn csi_levels(header: &Header, mut min_shift: u32, start_lvls: u32) -> (u32, u32, usize) {
    let mut max_len = 0i64;
    let mut nids = 0;
    for i in 0..header.contig_count() as u32 {
        if let Some(len) = header.contig_length(i) {
            max_len = max_len.max(len);
            nids += 1;
        }
    }
    if max_len == 0 {
        max_len = (1i64 << 31) - 1;
    }
    let mut n_lvls = start_lvls;
    adjust_csi_settings(max_len, &mut min_shift, &mut n_lvls);
    (n_lvls, min_shift, nids)
}

/// Indexes a bgzipped VCF (HTSlib's `tbx_index` with the VCF configuration).
fn tabix_index(reader: &mut Reader, min_shift: u32) -> Option<Tabix> {
    let mut builder = TabixBuilder::new(TabixConfig::VCF, min_shift);
    let mut line = Vec::new();
    loop {
        match reader.read_line(&mut line) {
            Ok(true) => {}
            Ok(false) => break,
            Err(_) => return None,
        }
        builder.push_line(&line, reader.virtual_position()).ok()?;
    }
    Some(builder.finish(reader.virtual_position()))
}

/// The data file and index that `index -s/-n` reads, from the name given.
struct StatsName {
    fname: String,
    /// The index named with `FILE##idx##INDEX`, or the name itself when it ends in `.csi` or
    /// `.tbi`.
    fnidx: Option<String>,
    /// Whether the name given is the index's.
    idx_only: bool,
}

impl StatsName {
    fn parse(fname_in: &str) -> Self {
        let mut fname = fname_in.to_owned();
        let mut fnidx: Option<String> = None;
        let mut idx_only = false;
        let lower = fname.to_ascii_lowercase();
        if let Some(i) = fname.find("##idx##") {
            fnidx = Some(fname[i + 7..].to_owned());
            fname.truncate(i);
        } else if fname.len() > 4 && (lower.ends_with(".csi") || lower.ends_with(".tbi")) {
            fnidx = Some(fname.clone());
            fname.truncate(fname.len() - 4);
            idx_only = true;
        }
        Self {
            fname,
            fnidx,
            idx_only,
        }
    }
}

/// What `index -s/-n` reads: a tabix or CSI index, and with `-s` the data file's header and
/// reader.
struct StatsSource {
    header: Option<Header>,
    reader: Option<Reader>,
    tbx: Option<Tabix>,
    idx: Option<Index>,
}

/// Prints index statistics (bcftools' `vcf_index_stats`); returns the exit status.
///
/// # Errors
///
/// Returns the error of a failed write, which bcftools ignores.
fn index_stats(fname_in: &str, stats: u32) -> io::Result<u8> {
    let name = StatsName::parse(fname_in);
    match open_stats_source(&name, stats) {
        Some(source) => print_index_stats(&name.fname, source, stats),
        None => Ok(1),
    }
}

/// Loads the index (and with `-s` the data file's header), reporting a failure.
fn open_stats_source(name: &StatsName, stats: u32) -> Option<StatsSource> {
    let fname = name.fname.as_str();
    let fnidx = name.fnidx.as_deref();
    let mut source = StatsSource {
        header: None,
        reader: None,
        tbx: None,
        idx: None,
    };
    if stats & PER_CONTIG != 0 {
        let idx_only = name.idx_only && !Path::new(fname).exists();
        let fmt = if idx_only {
            let ix = fnidx.unwrap_or("");
            if ix.to_ascii_lowercase().ends_with(".tbi") {
                InputFormat::Vcf
            } else {
                InputFormat::Bcf
            }
        } else {
            let (r, h) = open_with_header(fname)?;
            source.header = Some(h);
            let fmt = r.format();
            source.reader = Some(r);
            fmt
        };
        match fmt {
            InputFormat::Vcf => match load_tabix(fname, fnidx) {
                Some(t) => source.tbx = Some(t),
                None => {
                    eprintln!("Could not load index for VCF: {fname}");
                    return None;
                }
            },
            InputFormat::Bcf => match load_csi(fname, fnidx) {
                Some(i) => source.idx = Some(i),
                None => {
                    eprintln!("Could not load index for BCF file: {fname}");
                    return None;
                }
            },
        }
    } else if let Some(ix) = fnidx {
        // The text from the last '.', as C's `strrchr` finds it.
        match ix.rfind('.').map(|i| &ix[i..]) {
            Some(".tbi") => source.tbx = load_tabix(fname, Some(ix)),
            Some(".csi") => source.idx = load_csi(fname, Some(ix)),
            _ => {}
        }
        if source.tbx.is_none() && source.idx.is_none() {
            eprintln!("Could not load index file '{ix}'");
            return None;
        }
    } else if fname.ends_with(".bcf") {
        source.idx = load_csi(fname, None);
    } else if fname.len() > 7 && fname.ends_with(".vcf.gz") {
        source.tbx = load_tabix(fname, None);
    }
    if source.tbx.is_none() && source.idx.is_none() {
        eprintln!("No index file could be found for '{fname}'. Use 'vtools index' to create one");
        return None;
    }
    Some(source)
}

/// Opens the data file and reads its header, reporting a failure.
fn open_with_header(fname: &str) -> Option<(Reader, Header)> {
    let mut r = match Reader::open(fname, 0) {
        Ok(r) => r,
        Err(ReadError::UnknownBinary) => {
            report_exec_format_error(fname);
            eprintln!("Could not read {fname}");
            return None;
        }
        Err(ReadError::UnknownFormat) => {
            pasta::variant::log::error(
                "bcf_hdr_read",
                format_args!("Input is not detected as bcf or vcf format"),
            );
            eprintln!("Could not read the header: {fname}");
            return None;
        }
        Err(ReadError::Io(e)) => {
            pasta::variant::log::error(
                "hts_open_format",
                format_args!(
                    "Failed to open file \"{fname}\" : {}",
                    crate::raw_input::strerror(&e)
                ),
            );
            eprintln!("Could not read {fname}");
            return None;
        }
        Err(_) => {
            eprintln!("Could not read the header: {fname}");
            return None;
        }
    };
    match r.read_header() {
        Ok(h) => Some((r, h)),
        Err(_) => {
            eprintln!("Could not read the header: {fname}");
            None
        }
    }
}

/// Prints the record counts of the index: per contig (`-s`) or in total (`-n`).
fn print_index_stats(fname: &str, source: StatsSource, stats: u32) -> io::Result<u8> {
    let StatsSource {
        mut header,
        mut reader,
        tbx,
        idx,
    } = source;
    let nseq = match (&tbx, &idx) {
        (Some(t), _) => t.names().len(),
        (None, Some(i)) => i.nseq(),
        (None, None) => 0,
    };
    if tbx.is_none() && header.is_none() {
        eprintln!("Warning: cannot determine contig names given the .csi index alone");
    }
    let index = tbx
        .as_ref()
        .map_or_else(|| idx.as_ref().expect("an index"), |t| &t.index);
    let mut out = io::BufWriter::new(crate::output::stdout());
    let mut sum = 0u64;
    let mut ret = 0;
    for tid in 0..nseq {
        let stat = index.stat(tid);
        let records = stat.map_or(0, |s| s.0);
        sum += records;
        if stats & TOTAL != 0 || (records == 0 && stats & ALL_CONTIGS == 0) {
            continue;
        }
        let name: Vec<u8> = match (&tbx, &header) {
            (Some(t), _) => t.names()[tid].clone(),
            (None, Some(h)) => h.contig_name(tid as i32).unwrap_or(b"(null)").to_vec(),
            (None, None) => b"n/a".to_vec(),
        };
        let length = header
            .as_ref()
            .and_then(|h| h.definition(HeaderRecordType::Contig, &name))
            .and_then(|r| r.get(b"length"))
            .map_or_else(|| b".".to_vec(), <[u8]>::to_vec);
        out.write_all(&name)?;
        out.write_all(b"\t")?;
        out.write_all(&length)?;
        out.write_all(b"\t")?;
        match stat {
            Some(_) => writeln!(out, "{records}")?,
            None => writeln!(out, ".")?,
        }
    }
    if sum == 0
        && let (Some(r), Some(h)) = (&mut reader, &mut header)
    {
        let mut rec = Record::new();
        let mut line = Vec::new();
        if matches!(r.read_record(h, &mut line, &mut rec), Ok(true)) {
            eprintln!(
                "index of {fname} does not contain any count metadata. Please re-index with a newer version of vtools or tabix."
            );
            ret = 1;
        }
    }
    if stats & TOTAL != 0 && ret == 0 {
        writeln!(out, "{sum}")?;
    }
    out.flush()?;
    Ok(ret)
}

/// Loads a tabix-style index (`tbx_index_load2`): the given one, else `FILE.csi` or `FILE.tbi`.
fn load_tabix(fname: &str, fnidx: Option<&str>) -> Option<Tabix> {
    let path = match fnidx {
        Some(i) => PathBuf::from(i),
        None => locate_or_report(fname, IndexFormat::Tbi)?,
    };
    load_reporting(fname, &path, Tabix::load)
}

/// Loads a CSI index (`bcf_index_load2`): the given one, else `FILE.csi`.
fn load_csi(fname: &str, fnidx: Option<&str>) -> Option<Index> {
    let path = match fnidx {
        Some(i) => PathBuf::from(i),
        None => locate_or_report(fname, IndexFormat::Csi)?,
    };
    load_reporting(fname, &path, Index::read)
}

/// Loads an index as HTSlib's `hts_idx_load3` does: warning when it is older than the data
/// file, and reporting a failure.
fn load_reporting<T>(
    fname: &str,
    path: &Path,
    load: impl FnOnce(&Path) -> io::Result<T>,
) -> Option<T> {
    warn_if_older(Path::new(fname), path);
    load(path).inspect_err(|e| report_load_error(path, e)).ok()
}

/// Finds the index next to `fname`, reporting a missing one as HTSlib's `idx_find_and_load`.
fn locate_or_report(fname: &str, fmt: IndexFormat) -> Option<PathBuf> {
    let path = locate(fname, fmt);
    if path.is_none() {
        pasta::variant::log::error(
            "idx_find_and_load",
            format_args!("Could not retrieve index file for '{fname}'"),
        );
    }
    path
}
