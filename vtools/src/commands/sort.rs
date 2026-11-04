//! `vtools sort`: sort VCF/BCF files by contig, position, and alleles (bcftools sort).
//!
//! Records are packed into one memory arena with their sort keys. Memory is counted as bcftools
//! counts its packed records, so a given `--max-mem` spills the same blocks. When the input
//! fits, the arena is sorted and written directly; otherwise sorted chunks spill to BGZF
//! temporary files that are merged with bcftools' layered scheme (32 files per layer).
//! The order is bcftools': contig ID, position, then REF and ALT joined by tabs (compared as
//! bytes), with ties kept in input order.

use std::{
    cmp::Ordering,
    ffi::{CStr, CString},
    fs::File,
    io::{self, IsTerminal, Read, Write},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::{Context, bail};
use pasta::{
    bgzf::{self, io::writer::CompressionLevel},
    variant::{
        Header, Record, Unpack, bcf,
        io::{Compression, ReadError, Reader, WriteError, Writer},
        log,
    },
};

use super::head::apply_verbosity;
use crate::{
    output::{self, OutputSetting},
    raw_input::{report_exec_format_error, strerror},
};

const USAGE: &str = "
About:   Sort VCF/BCF file.
Usage:   vtools sort [OPTIONS] <FILE.vcf>

Options:
    -m, --max-mem FLOAT[kMG]       Maximum memory to use [768M]
    -o, --output FILE              Output file name [stdout]
    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -T, --temp-dir DIR             Temporary files [/tmp/bcftools.XXXXXX]
    -v, --verbosity INT            Verbosity level
    -W, --write-index[=FMT]        Automatically index the output files [off]

";

/// bcftools' `MAX_TMP_FILES_PER_LAYER`.
const MAX_TMP_FILES_PER_LAYER: u32 = 32;
/// bcftools' `MERGE_LAYERS`.
const MERGE_LAYERS: usize = 12;
/// The memory limit without `-m` (bcftools' 768M).
const DEFAULT_MAX_MEM: u64 = 768_000_000;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file [default: stdin].
    inputs: Vec<String>,

    /// Maximum memory to use [768M].
    #[arg(
        short = 'm',
        long = "max-mem",
        value_name = "FLOAT[kMG]",
        allow_hyphen_values = true
    )]
    max_mem: Option<String>,

    /// Output file name [stdout].
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "output-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "u|b|v|z[0-9]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output_type: Vec<String>,

    /// Directory prefix for temporary files [$TMPDIR/bcftools.XXXXXX].
    #[arg(
        short = 'T',
        long = "temp-dir",
        value_name = "DIR",
        allow_hyphen_values = true
    )]
    temp_dir: Option<String>,

    /// Verbosity level.
    #[arg(
        short = 'v',
        long = "verbosity",
        value_name = "INT",
        allow_hyphen_values = true
    )]
    verbosity: Option<String>,

    /// Automatically index the output files [off]: csi (default) or tbi.
    #[arg(
        short = 'W',
        long = "write-index",
        value_name = "FMT",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = ""
    )]
    write_index: Option<String>,

    /// Print help.
    #[arg(short = 'h', long = "help", short_alias = '?')]
    help: bool,
}

/// bcftools' `parse_mem_string`: a number with an optional k, M, or G suffix (powers of 1000).
fn parse_mem_string(s: &str) -> anyhow::Result<u64> {
    let bytes = s.as_bytes();
    let (mem, end) = pasta::variant::number::strtod(bytes, 0);
    if end == 0 {
        bail!("Could not parse the memory string: \"{s}\"");
    }
    let suffix = &s[end..];
    let mem = if suffix.eq_ignore_ascii_case("k") {
        mem * 1e3
    } else if suffix.eq_ignore_ascii_case("m") {
        mem * 1e6
    } else if suffix.eq_ignore_ascii_case("g") {
        mem * 1e9
    } else {
        mem
    };
    Ok(c_size_t(mem))
}

/// C's conversion of a double to `size_t` as GCC compiles it for x86-64: truncation, with
/// negative values wrapping (so `-m -5` asks for nearly 2^64 bytes), NaN giving 2^63, and
/// values of 2^64 or more (infinity) giving 0.
fn c_size_t(x: f64) -> u64 {
    const TWO_64: f64 = 18_446_744_073_709_551_616.0;
    if x.is_nan() {
        1 << 63
    } else if x < 0.0 {
        #[expect(
            clippy::cast_possible_truncation,
            reason = "cvttsd2si truncates, saturating to i64::MIN"
        )]
        let v = x as i64;
        v.cast_unsigned()
    } else if x >= TWO_64 {
        0
    } else {
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "x is in [0, 2^64), where the conversion truncates"
        )]
        let v = x as u64;
        v
    }
}

/// The memory bcftools sorts in: 90% of `-m` (the rest is for its bookkeeping), which it
/// allocates up front. More than `PTRDIFF_MAX` bytes (a negative `-m`) cannot be allocated.
fn sort_memory(max_mem: u64) -> anyhow::Result<usize> {
    #[expect(
        clippy::cast_precision_loss,
        reason = "bcftools multiplies the size_t by a double"
    )]
    let limit = c_size_t(max_mem as f64 * 0.9);
    match usize::try_from(limit) {
        Ok(l) if isize::try_from(l).is_ok() => Ok(l),
        _ => bail!("Error: could not allocate {limit} bytes of memory, try reducing --max-mem"),
    }
}

/// The temporary directory (bcftools' `init_tmp_prefix` and `mkdtemp`), removed with what it
/// holds when dropped.
struct TempDir {
    path: PathBuf,
    files: Vec<PathBuf>,
    count: usize,
    /// Whether a closed output pipe removes the directory too (see [`tmp_on_sigpipe`]).
    removed_on_sigpipe: bool,
}

impl TempDir {
    fn create(prefix: Option<&str>) -> anyhow::Result<Self> {
        let template = match prefix {
            Some(p) => format!("{p}XXXXXX"),
            None => {
                let dir = std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".to_owned());
                format!("{dir}/bcftools.XXXXXX")
            }
        };
        let c = CString::new(template.clone()).context("invalid temporary directory name")?;
        let mut buf = c.into_bytes_with_nul();
        // SAFETY: `buf` is a writable NUL-terminated template, as mkdtemp requires.
        let p = unsafe { libc::mkdtemp(buf.as_mut_ptr().cast()) };
        if p.is_null() {
            let e = io::Error::last_os_error();
            // mkdtemp leaves its last try in the template, which bcftools prints.
            let tried = CStr::from_bytes_until_nul(&buf)
                .map_or_else(|_| template.clone(), |s| s.to_string_lossy().into_owned());
            bail!("mkdtemp({tried}) failed: {}", strerror(&e));
        }
        // SAFETY: mkdtemp returned `buf`, still NUL-terminated.
        let path = unsafe { CStr::from_ptr(buf.as_ptr().cast()) };
        let path = PathBuf::from(path.to_string_lossy().into_owned());
        Ok(Self {
            path,
            files: Vec::new(),
            count: 0,
            removed_on_sigpipe: false,
        })
    }

    /// Has the directory and its files removed when a write to a closed pipe kills the
    /// program, which leaves them behind in bcftools.
    fn remove_on_sigpipe(&mut self) {
        tmp_on_sigpipe::remove_on_death(&self.path);
        self.removed_on_sigpipe = true;
    }

    /// Returns a new file name in the directory (bcftools' `open_tmp_file`).
    fn next_name(&mut self) -> PathBuf {
        let p = self.path.join(format!("{:05}", self.count));
        self.count += 1;
        if self.removed_on_sigpipe {
            tmp_on_sigpipe::set_files(self.count);
        }
        self.files.push(p.clone());
        p
    }

    fn remove(&mut self, path: &Path) {
        let _ = std::fs::remove_file(path);
        self.files.retain(|f| f != path);
    }

    /// Removes the files and the directory (bcftools' `clean_files`).
    fn clean(self) {
        eprintln!("Cleaning");
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        for f in &self.files {
            let _ = std::fs::remove_file(f);
        }
        let _ = std::fs::remove_dir(&self.path);
        if self.removed_on_sigpipe {
            tmp_on_sigpipe::keep_on_death();
        }
    }
}

/// Removing the temporary directory when a write to a closed pipe (`| head`) kills the program:
/// a SIGPIPE handler removes the directory and its files, named as [`TempDir::next_name`] names
/// them, and then lets SIGPIPE's default action end the program as it ends bcftools.
mod tmp_on_sigpipe {
    use std::{
        ffi::{CStr, CString},
        os::unix::ffi::OsStrExt,
        path::Path,
        ptr,
        sync::atomic::{AtomicPtr, AtomicUsize, Ordering},
    };

    /// The directory to remove, NUL-terminated, or null. A directory is never freed, as the
    /// handler may be reading it.
    static DIR: AtomicPtr<libc::c_char> = AtomicPtr::new(ptr::null_mut());
    /// How many file names the directory has given out.
    static FILES: AtomicUsize = AtomicUsize::new(0);

    /// Installs the handler that removes `dir` and the files it gives out.
    pub(super) fn remove_on_death(dir: &Path) {
        let Ok(dir) = CString::new(dir.as_os_str().as_bytes()) else {
            return;
        };
        FILES.store(0, Ordering::Release);
        DIR.store(dir.into_raw(), Ordering::Release);
        // SAFETY: `action` is zeroed and then filled in as sigaction(2) requires; the handler
        // only makes async-signal-safe calls.
        unsafe {
            let mut action: libc::sigaction = std::mem::zeroed();
            action.sa_sigaction =
                remove_and_die as extern "C" fn(libc::c_int) as libc::sighandler_t;
            libc::sigemptyset(&mut action.sa_mask);
            libc::sigaction(libc::SIGPIPE, &action, ptr::null_mut());
        }
    }

    /// Records that the directory has given out `n` file names.
    pub(super) fn set_files(n: usize) {
        FILES.store(n, Ordering::Release);
    }

    /// Leaves the directory alone: it has been removed.
    pub(super) fn keep_on_death() {
        DIR.store(ptr::null_mut(), Ordering::Release);
    }

    /// The SIGPIPE handler. SIGPIPE stays blocked while it runs, so the signal it raises once
    /// the default action is back ends the program as soon as it returns.
    extern "C" fn remove_and_die(_: libc::c_int) {
        let dir = DIR.load(Ordering::Acquire);
        if !dir.is_null() {
            // SAFETY: a non-null `DIR` is a NUL-terminated path that is never freed.
            let bytes = unsafe { CStr::from_ptr(dir) }.to_bytes();
            let mut buf = [0u8; libc::PATH_MAX as usize];
            for i in 0..FILES.load(Ordering::Acquire) {
                if let Some(name) = file_name(bytes, i, &mut buf) {
                    // SAFETY: `name` is NUL-terminated; unlink is async-signal-safe.
                    unsafe { libc::unlink(name.as_ptr()) };
                }
            }
            // SAFETY: `dir` is NUL-terminated; rmdir is async-signal-safe.
            unsafe { libc::rmdir(dir) };
        }
        // SAFETY: signal and raise are async-signal-safe, and SIG_DFL installs no handler.
        unsafe {
            libc::signal(libc::SIGPIPE, libc::SIG_DFL);
            libc::raise(libc::SIGPIPE);
        }
    }

    /// Writes the name of file `i` of `dir` (`DIR/00042`, as `format!("{:05}")` pads it) into
    /// `buf`, NUL-terminated, without allocating; `None` when it does not fit.
    pub(super) fn file_name<'b>(dir: &[u8], i: usize, buf: &'b mut [u8]) -> Option<&'b CStr> {
        let mut digits = [0u8; 20];
        let mut len = 0;
        let mut n = i;
        loop {
            digits[len] = b'0' + (n % 10) as u8;
            len += 1;
            n /= 10;
            if n == 0 {
                break;
            }
        }
        let width = len.max(5);
        let end = dir.len() + 1 + width;
        let name = buf.get_mut(..=end)?;
        name[..dir.len()].copy_from_slice(dir);
        name[dir.len()] = b'/';
        let (pad, number) = name[dir.len() + 1..end].split_at_mut(width - len);
        pad.fill(b'0');
        for (d, &digit) in number.iter_mut().zip(digits[..len].iter().rev()) {
            *d = digit;
        }
        name[end] = 0;
        CStr::from_bytes_with_nul(name).ok()
    }
}

/// One record in the arena: its sort key and where its alleles and packed data are.
#[derive(Clone, Copy, Debug)]
struct Entry {
    rid: i32,
    key_len: u32,
    pos: i64,
    /// Where the alleles start in the arena's data; the packed record follows them.
    data: usize,
}

/// What bcftools sort reserves in its memory block for any record besides the alleles and the
/// encoded blocks (the `delta` of `buf_push`): 16 bytes for the position, contig, and QUAL;
/// the longest varint of each field's C type (`varint_size`), 10 bytes for the 64-bit `rlen`,
/// 5 each for `n_info`, `n_allele`, `n_fmt`, and `n_sample` (counted as `int`), and 10 each
/// for the `size_t` lengths of the shared and individual blocks; and a flat 8 for alignment
/// (`_align_up`). 16 + 10 + 4 * 5 + 2 * 10 + 8 = 74.
const BLOCK_RESERVE: usize = 74;

/// Records held in memory (bcftools' `mem_block` and `buf`). It holds each record's alleles
/// (the tie-breaking sort key) followed by its packed form.
#[derive(Default)]
struct Arena {
    data: Vec<u8>,
    entries: Vec<Entry>,
    /// The bytes the records take in bcftools' memory block (its `mem`), which decides when the
    /// block is full, so that a given `-m` spills as many temporary files as bcftools does.
    block_used: usize,
}

impl Arena {
    /// Adds a record, returning `false` when it does not fit in the `limit` bytes of bcftools'
    /// memory block: the block is then written out with this record, as bcftools writes it.
    fn push(&mut self, rec: &mut Record, limit: usize) -> bool {
        let data = self.data.len();
        push_key(rec, &mut self.data);
        let key_len = self.data.len() - data;
        bcf::pack_record(rec, &mut self.data);
        self.entries.push(Entry {
            rid: rec.rid(),
            key_len: key_len as u32,
            pos: rec.pos(),
            data,
        });
        let need = BLOCK_RESERVE
            + rec.encoded_alleles_len()
            + rec.shared_bytes().len()
            + rec.indiv_bytes().len();
        if self.block_used + need > limit {
            return false;
        }
        self.block_used = self.block_used.next_multiple_of(8) + block_packed_len(rec, key_len);
        true
    }

    /// The record's alleles joined by tabs.
    fn key(&self, e: &Entry) -> &[u8] {
        &self.data[e.data..e.data + e.key_len as usize]
    }

    /// The data from the record's packed form on.
    fn packed(&self, e: &Entry) -> &[u8] {
        &self.data[e.data + e.key_len as usize..]
    }

    /// Sorts stably by contig, position, and alleles (bcftools'
    /// `cmp_packed_bcf_pos_ref_alt_stable`).
    fn sort(&mut self) {
        let data = &self.data;
        let key = |e: &Entry| &data[e.data..e.data + e.key_len as usize];
        self.entries.sort_by(|a, b| {
            a.rid
                .cmp(&b.rid)
                .then(a.pos.cmp(&b.pos))
                .then_with(|| key(a).cmp(key(b)))
        });
    }

    fn clear(&mut self) {
        self.data.clear();
        self.entries.clear();
        self.block_used = 0;
    }
}

/// The length of bcftools' packed form of a record (`pack_bcf_data`): its own length, the
/// position, contig, and QUAL in 24 bytes, the alleles joined by tabs and terminated, the other
/// fixed fields and the block lengths as varints, then the encoded blocks.
fn block_packed_len(rec: &Record, key_len: usize) -> usize {
    let rlen = rec.rlen();
    let rlen_zigzag = ((rlen << 1) ^ (rlen >> 63)).cast_unsigned();
    let shared = rec.shared_bytes().len();
    let indiv = rec.indiv_bytes().len();
    24 + key_len
        + 1
        + varint_len(rlen_zigzag)
        + varint_len(rec.info_count() as u64)
        + varint_len(rec.allele_count() as u64)
        + varint_len(rec.format_count() as u64)
        + varint_len(rec.sample_count() as u64)
        + varint_len(shared as u64)
        + varint_len(indiv as u64)
        + shared
        + indiv
}

/// The bytes of `v` as an unsigned LEB128 varint.
fn varint_len(v: u64) -> usize {
    (u64::BITS - (v | 1).leading_zeros()).div_ceil(7) as usize
}

/// Appends the record's alleles joined by tabs: the tie-breaking sort key.
fn push_key(rec: &mut Record, out: &mut Vec<u8>) {
    rec.unpack(Unpack::STR);
    for (i, a) in rec.alleles().enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(a);
    }
}

/// A BGZF temporary file, compressed at level 1 as bcftools' are.
type TempWriter = bgzf::io::Writer<io::BufWriter<File>>;

fn temp_writer(file: File) -> TempWriter {
    bgzf::io::writer::Builder::default()
        .set_compression_level(CompressionLevel::FAST)
        .build_from_writer(io::BufWriter::new(file))
}

/// A sorted temporary file (bcftools' `blk_t`), in creation order.
struct Block {
    path: PathBuf,
    /// The creation index, which breaks ties between blocks so that the sort is stable.
    idx: usize,
}

/// An error that bcftools reports and then cleans up after.
fn fail(msg: impl Into<String>) -> anyhow::Error {
    anyhow::anyhow!(msg.into())
}

/// bcftools' errors for an input that cannot be opened (`hts_open` fails, with HTSlib's
/// messages) or that is neither VCF nor BCF (`hts_open` succeeds and `bcf_hdr_read` fails).
fn open_error(fname: &str, e: &ReadError) -> anyhow::Error {
    match e {
        ReadError::UnknownBinary => {
            report_exec_format_error(fname);
            fail(format!("Could not read {fname}"))
        }
        ReadError::UnknownFormat => {
            log::error(
                "bcf_hdr_read",
                format_args!("Input is not detected as bcf or vcf format"),
            );
            fail(format!("Could not read VCF/BCF headers from {fname}"))
        }
        ReadError::Io(e) => {
            output::log_open_failure(fname, e);
            fail(format!("Could not read {fname}"))
        }
        _ => fail(format!("Could not read {fname}")),
    }
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    if args.help {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    }
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    let max_mem = match &args.max_mem {
        Some(m) => parse_mem_string(m)?,
        None => DEFAULT_MAX_MEM,
    };
    let settings: Vec<OutputSetting> = args
        .output_type
        .iter()
        .map(|t| OutputSetting::Type(t.clone()))
        .collect();
    // bcftools checks -O while parsing the options, before making the temporary directory.
    output::requested_type(&settings)?;
    let write_index = args
        .write_index
        .as_deref()
        .map(output::parse_write_index)
        .transpose()?;
    let fname = match args.inputs.first() {
        Some(f) => f.clone(),
        None if !io::stdin().is_terminal() => "-".to_owned(),
        None => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
    };
    let limit = sort_memory(max_mem)?;

    let mut tmp = TempDir::create(args.temp_dir.as_deref())?;
    tmp.remove_on_sigpipe();
    eprintln!("Writing to {}", tmp.path.display());
    let output = Output {
        path: args.output.as_deref(),
        settings: &settings,
        write_index,
    };
    match sort_to_output(&fname, limit, &output, &mut tmp) {
        Ok(()) => {
            tmp.clean();
            Ok(ExitCode::SUCCESS)
        }
        // A closed pipe while SIGPIPE was blocked: the directory goes, without the clean-up
        // message, and `main` dies of SIGPIPE (see `crate::sigpipe`).
        Err(e) if is_broken_pipe(&e) => Err(e),
        // bcftools' `error()` exits without its clean-up message (and leaves the temporary
        // files, which are removed here all the same).
        Err(e) if e.is::<NotCleaned>() => Err(e),
        Err(e) => {
            // bcftools' `clean_files_and_throw`: the message, then the clean-up.
            let msg = format!("{e:#}");
            if msg.ends_with('\n') {
                eprint!("{msg}");
            } else {
                eprintln!("{msg}");
            }
            tmp.clean();
            Ok(ExitCode::from(255))
        }
    }
}

fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|cause| cause.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::BrokenPipe)
}

/// Where and how the sorted records are written.
struct Output<'a> {
    path: Option<&'a Path>,
    settings: &'a [OutputSetting],
    write_index: Option<output::WriteIndex>,
}

/// Reads and sorts the input, then writes it (bcftools' `sort_blocks` and `merge_to_output`).
fn sort_to_output(
    fname: &str,
    limit: usize,
    out: &Output<'_>,
    tmp: &mut TempDir,
) -> anyhow::Result<()> {
    let settings = out.settings;
    let write_index = out.write_index;
    let mut sorter = Sorter {
        arena: Arena::default(),
        limit,
        blocks: Vec::new(),
        layers: [0; MERGE_LAYERS],
    };
    let mut header = sorter.read_input(fname, tmp)?;

    let (out_type, level) = output::resolve(settings, out.path)?;
    let out_name = out
        .path
        .map_or_else(|| "-".to_owned(), |p| p.display().to_string());
    let out_path = out.path.map(output::data_path);
    let mut writer = Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, 0),
    )
    .map_err(|e| {
        output::log_open_failure(&out_name, &e);
        fail(format!("[merge_to_output] Error: cannot open {out_name}"))
    })?;
    // bcftools writes the last block to a file too, so it counts it.
    eprintln!(
        "Merging {} temporary files",
        sorter.blocks.len() + usize::from(!sorter.arena.entries.is_empty())
    );
    if writer.write_header(&header).is_err() {
        writer.abandon();
        bail!("[merge_blocks] Error: cannot write to {out_name}");
    }
    if let Some(fmt) = write_index {
        let ok = match output::index_request(fmt, out_type, out.path) {
            Some(req) => writer.request_index(req).is_ok(),
            None => false,
        };
        if !ok {
            writer.abandon();
            return Err(
                NotCleaned(format!("Error: failed to initialise index for {out_name}")).into(),
            );
        }
        writer.start_index(&header)?;
    }
    let result = if sorter.blocks.is_empty() {
        sorter.write_arena(&mut header, &mut writer)
    } else {
        let blocks = std::mem::take(&mut sorter.blocks);
        merge_blocks(
            &blocks,
            Some(&sorter.arena),
            tmp,
            &mut header,
            &mut BlockOutput::Final(&mut writer),
        )
    };
    if let Err(e) = result {
        writer.abandon();
        // A record the output cannot take (HTSlib has logged why) is bcftools' write error;
        // other errors (reading temporary files, a closed pipe) are reported as they are.
        if e.is::<WriteError>() && !is_broken_pipe(&e) {
            bail!("[merge_blocks] Error: cannot write to {out_name}");
        }
        return Err(e);
    }
    let cannot_write_index = || -> anyhow::Error {
        let index = write_index
            .and_then(|fmt| output::index_request(fmt, out_type, out.path))
            .map(|req| req.path.display().to_string())
            .unwrap_or_default();
        NotCleaned(format!("Error: cannot write to index {index}")).into()
    };
    // bcftools prints "Done" before it closes the output, so also when the close fails or its
    // flush finds the reader of the output gone (SIGPIPE). With an index it saves the index
    // first, at the end of merge_blocks; pasta saves it when closing, so "Done" waits for that.
    if write_index.is_none() {
        eprintln!("Done");
    }
    let result = writer.finish();
    if write_index.is_some() && !matches!(result, Err(WriteError::SaveIndex(_))) {
        eprintln!("Done");
    }
    match result {
        Ok(()) => Ok(()),
        Err(WriteError::SaveIndex(_) | WriteError::Index(_)) => Err(cannot_write_index()),
        Err(e) => Err(output::write_failure(e, || {
            fail(format!("Close failed: {out_name}"))
        })),
    }
}

/// An error bcftools reports with `error()` rather than `clean_files_and_throw()`, so without
/// the "Cleaning" message.
#[derive(Debug)]
struct NotCleaned(String);

impl std::fmt::Display for NotCleaned {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for NotCleaned {}

struct Sorter {
    arena: Arena,
    limit: usize,
    blocks: Vec<Block>,
    /// The number of blocks in each merge layer (bcftools' `tmp_layers`).
    layers: [u32; MERGE_LAYERS],
}

impl Sorter {
    /// Reads every record into the arena, spilling sorted blocks when it is full (bcftools'
    /// `sort_blocks`). Returns the header, with any contigs and tags the records added.
    fn read_input(&mut self, fname: &str, tmp: &mut TempDir) -> anyhow::Result<Header> {
        let mut reader = Reader::open(fname, 0).map_err(|e| open_error(fname, &e))?;
        let mut header = reader
            .read_header()
            .map_err(|_| fail(format!("Could not read VCF/BCF headers from {fname}")))?;
        let mut rec = Record::new();
        let mut line = Vec::new();
        loop {
            match reader.read_record(&mut header, &mut line, &mut rec) {
                Ok(true) => {}
                Ok(false) => {
                    // HTSlib's warning for a file cut at a BGZF block boundary (a pipe is
                    // checked as it is read).
                    if reader.compression() == Compression::Bgzf {
                        super::index::warn_if_no_eof(fname);
                    }
                    break;
                }
                Err(_) => bail!("Error encountered while parsing the input"),
            }
            if !rec.errors().is_empty() {
                bail!(
                    "Error encountered while parsing the input at {}:{}",
                    String::from_utf8_lossy(rec.chrom_safe(&header)),
                    rec.pos() + 1
                );
            }
            if !self.arena.push(&mut rec, self.limit) {
                self.flush(tmp, &mut header)?;
            }
        }
        // The last block stays in memory, merged with the temporary files. bcftools writes it
        // to a file like the others, merging full layers first, so that is done here too.
        if !self.arena.entries.is_empty() && self.layers[0] >= MAX_TMP_FILES_PER_LAYER {
            self.partial_merge(tmp, &mut header)?;
        }
        self.arena.sort();
        Ok(header)
    }

    /// Sorts the arena and writes it to a new temporary file (bcftools' `buf_flush`).
    fn flush(&mut self, tmp: &mut TempDir, header: &mut Header) -> anyhow::Result<()> {
        if self.arena.entries.is_empty() {
            return Ok(());
        }
        self.arena.sort();
        if self.layers[0] >= MAX_TMP_FILES_PER_LAYER {
            self.partial_merge(tmp, header)?;
        }
        let path = tmp.next_name();
        let idx = tmp.count - 1;
        let file = File::create_new(&path).map_err(|e| {
            fail(format!(
                "Cannot write {}: {}\n",
                path.display(),
                strerror(&e)
            ))
        })?;
        let mut w = temp_writer(file);
        let write_err = || {
            fail(format!(
                "[buf_flush] Error: cannot write to {}",
                path.display()
            ))
        };
        for e in &self.arena.entries {
            let packed = self.arena.packed(e);
            let len = bcf::packed_len(packed).expect("packed in the arena");
            w.write_all(&packed[..len]).map_err(|_| write_err())?;
        }
        w.finish().map_err(|_| {
            fail(format!(
                "[buf_flush] Error: close failed .. {}",
                path.display()
            ))
        })?;
        self.blocks.push(Block { path, idx });
        self.layers[0] += 1;
        self.arena.clear();
        Ok(())
    }

    /// Merges full layers into one file of the next layer (bcftools' `do_partial_merge`).
    fn partial_merge(&mut self, tmp: &mut TempDir, header: &mut Header) -> anyhow::Result<()> {
        let mut to_layer = 0;
        let mut to_merge = 0usize;
        while to_layer < MERGE_LAYERS && self.layers[to_layer] >= MAX_TMP_FILES_PER_LAYER {
            to_merge += self.layers[to_layer] as usize;
            self.layers[to_layer] = 0;
            to_layer += 1;
        }
        if to_layer == MERGE_LAYERS {
            to_layer = MERGE_LAYERS - 1;
        }
        let path = tmp.next_name();
        let idx = tmp.count - 1;
        let file = File::create_new(&path).map_err(|e| {
            fail(format!(
                "Cannot write {}: {}\n",
                path.display(),
                strerror(&e)
            ))
        })?;
        let mut w = temp_writer(file);
        let from = self.blocks.len() - to_merge;
        let merged: Vec<Block> = self.blocks.drain(from..).collect();
        merge_blocks(&merged, None, tmp, header, &mut BlockOutput::Temp(&mut w))?;
        w.finish()
            .map_err(|_| fail(format!("Close failed: {}", path.display())))?;
        self.blocks.push(Block { path, idx });
        self.layers[to_layer] += 1;
        Ok(())
    }

    /// Writes the sorted arena (the in-memory path: no temporary files).
    fn write_arena(&mut self, header: &mut Header, writer: &mut Writer) -> anyhow::Result<()> {
        let mut rec = Record::new();
        for e in &self.arena.entries {
            bcf::unpack_record(self.arena.packed(e), &mut rec)?;
            writer.write_record(header, &mut rec)?;
        }
        Ok(())
    }
}

/// Where merged records go: the output, or a temporary file of a partial merge.
enum BlockOutput<'a> {
    Final(&'a mut Writer),
    Temp(&'a mut TempWriter),
}

/// A block being merged: its reader and current record.
struct Source {
    /// The temporary file, or `None` for the in-memory block.
    reader: Option<bgzf::io::Reader<File>>,
    /// The next entry of the in-memory block.
    entry: usize,
    idx: usize,
    rec: Record,
    data: Vec<u8>,
    key: Vec<u8>,
}

impl Source {
    /// Reads the next packed record; `false` at the end of the block.
    fn next(&mut self, arena: Option<&Arena>) -> io::Result<bool> {
        let Some(reader) = self.reader.as_mut() else {
            let arena = arena.expect("the in-memory block is merged with the arena");
            let Some(e) = arena.entries.get(self.entry) else {
                return Ok(false);
            };
            self.entry += 1;
            bcf::unpack_record(arena.packed(e), &mut self.rec)?;
            self.key.clear();
            self.key.extend_from_slice(arena.key(e));
            return Ok(true);
        };
        self.data.clear();
        self.data.resize(bcf::PACKED_HEADER_LEN, 0);
        let n = read_up_to(reader, &mut self.data)?;
        if n == 0 {
            return Ok(false);
        }
        if n != bcf::PACKED_HEADER_LEN {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated temporary file",
            ));
        }
        let len = bcf::packed_len(&self.data).expect("the fixed part was read");
        self.data.resize(len, 0);
        reader.read_exact(&mut self.data[bcf::PACKED_HEADER_LEN..])?;
        bcf::unpack_record(&self.data, &mut self.rec)?;
        self.key.clear();
        push_key(&mut self.rec, &mut self.key);
        Ok(true)
    }

    fn cmp(&self, other: &Self) -> Ordering {
        self.rec
            .rid()
            .cmp(&other.rec.rid())
            .then(self.rec.pos().cmp(&other.rec.pos()))
            .then_with(|| self.key.cmp(&other.key))
            .then(self.idx.cmp(&other.idx))
    }
}

fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// Merges sorted blocks, and the sorted in-memory block when given, into `out` and removes the
/// blocks' files (bcftools' `merge_blocks`).
fn merge_blocks(
    blocks: &[Block],
    arena: Option<&Arena>,
    tmp: &mut TempDir,
    header: &mut Header,
    out: &mut BlockOutput<'_>,
) -> anyhow::Result<()> {
    let mut sources = Vec::with_capacity(blocks.len());
    for b in blocks {
        let file = File::open(&b.path).map_err(|e| {
            fail(format!(
                "Could not read {}: {}",
                b.path.display(),
                strerror(&e)
            ))
        })?;
        let mut s = Source {
            reader: Some(bgzf::io::Reader::new(file)),
            entry: 0,
            idx: b.idx,
            rec: Record::new(),
            data: Vec::new(),
            key: Vec::new(),
        };
        let read_err = |b: &Block| fail(format!("Error reading {}", b.path.display()));
        if s.next(arena).map_err(|_| read_err(b))? {
            sources.push(s);
        }
    }
    if arena.is_some() {
        // The in-memory block holds the last records read, so it sorts after the files.
        let mut s = Source {
            reader: None,
            entry: 0,
            idx: usize::MAX,
            rec: Record::new(),
            data: Vec::new(),
            key: Vec::new(),
        };
        if s.next(arena)? {
            sources.push(s);
        }
    }
    // A binary min-heap of source indices (bcftools' `khp_blk_t`).
    let mut heap: Vec<usize> = Vec::with_capacity(sources.len());
    for i in 0..sources.len() {
        heap.push(i);
        let last = heap.len() - 1;
        sift_up(&mut heap, &sources, last);
    }
    while let Some(&top) = heap.first() {
        let s = &mut sources[top];
        match out {
            BlockOutput::Final(w) => w.write_record(header, &mut s.rec)?,
            BlockOutput::Temp(w) => w.write_all(&s.data)?,
        }
        let more = s
            .next(arena)
            .map_err(|_| fail(format!("Error reading {}", blocks_path(blocks, s.idx))))?;
        if more {
            sift_down(&mut heap, &sources, 0);
        } else {
            let last = heap.pop().expect("the heap is not empty");
            if !heap.is_empty() {
                heap[0] = last;
                sift_down(&mut heap, &sources, 0);
            }
        }
    }
    for b in blocks {
        tmp.remove(&b.path);
    }
    Ok(())
}

fn blocks_path(blocks: &[Block], idx: usize) -> String {
    blocks
        .iter()
        .find(|b| b.idx == idx)
        .map_or_else(String::new, |b| b.path.display().to_string())
}

fn sift_up(heap: &mut [usize], sources: &[Source], mut i: usize) {
    while i > 0 {
        let parent = (i - 1) / 2;
        if sources[heap[i]].cmp(&sources[heap[parent]]) == Ordering::Less {
            heap.swap(i, parent);
            i = parent;
        } else {
            break;
        }
    }
}

fn sift_down(heap: &mut [usize], sources: &[Source], mut i: usize) {
    loop {
        let l = 2 * i + 1;
        let r = l + 1;
        let mut min = i;
        if l < heap.len() && sources[heap[l]].cmp(&sources[heap[min]]) == Ordering::Less {
            min = l;
        }
        if r < heap.len() && sources[heap[r]].cmp(&sources[heap[min]]) == Ordering::Less {
            min = r;
        }
        if min == i {
            break;
        }
        heap.swap(i, min);
        i = min;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn varints_take_as_many_bytes_as_bcftools_packs() {
        // pack_unsigned writes 7 bits a byte. Within HTSlib's 8-bit n_fmt this is the one or
        // two bytes the FORMAT count used to be counted as.
        for n in 0..=255_u64 {
            assert_eq!(varint_len(n), if n < 0x80 { 1 } else { 2 }, "{n}");
        }
        for (v, len) in [
            (16_383, 2),
            (16_384, 3),
            (u64::from(u32::MAX), 5),
            (u64::MAX, 10),
        ] {
            assert_eq!(varint_len(v), len, "{v}");
        }
    }

    #[test]
    fn memory_strings_parse_as_bcftools() {
        assert_eq!(parse_mem_string("768M").unwrap(), 768_000_000);
        assert_eq!(parse_mem_string("1.5k").unwrap(), 1500);
        assert_eq!(parse_mem_string("2g").unwrap(), 2_000_000_000);
        assert_eq!(parse_mem_string("100").unwrap(), 100);
        assert!(parse_mem_string("x").is_err());
        // An unknown suffix is ignored, as strcasecmp finds no match.
        assert_eq!(parse_mem_string("1Mb").unwrap(), 1);
        assert_eq!(parse_mem_string("0x10k").unwrap(), 16_000);
        // Negative values wrap, infinity converts to 0.
        assert_eq!(parse_mem_string("-5").unwrap(), u64::MAX - 4);
        assert_eq!(parse_mem_string("inf").unwrap(), 0);
    }

    /// Sorts the VCF `body` (after a header with `samples` samples) in `max_mem` bytes as far as
    /// the files to merge, returning how many bcftools reports ("Merging N temporary files").
    fn merged_files(name: &str, samples: usize, body: &str, max_mem: u64) -> usize {
        let dir = std::env::temp_dir().join(format!("vtools-sort-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n\
             ##INFO=<ID=S,Number=1,Type=String,Description=\"s\">\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
             ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        );
        if samples > 0 {
            text.push_str("\tFORMAT");
            (0..samples).for_each(|i| text.push_str(&format!("\ts{i}")));
        }
        text.push('\n');
        text.push_str(body);
        let input = dir.join("in.vcf");
        std::fs::write(&input, text).unwrap();
        let mut tmp = TempDir::create(Some(&format!("{}/tmp.", dir.display()))).unwrap();
        let mut sorter = Sorter {
            arena: Arena::default(),
            limit: sort_memory(max_mem).unwrap(),
            blocks: Vec::new(),
            layers: [0; MERGE_LAYERS],
        };
        sorter
            .read_input(input.to_str().unwrap(), &mut tmp)
            .unwrap();
        let files = sorter.blocks.len() + usize::from(!sorter.arena.entries.is_empty());
        tmp.clean();
        std::fs::remove_dir_all(&dir).unwrap();
        files
    }

    #[test]
    fn memory_is_counted_as_bcftools_counts_it() {
        // The counts are bcftools 1.24's for the same files and -m.
        let small = |n: usize| -> String {
            (0..n)
                .map(|i| format!("chr1\t{}\t.\tA\tC\t50\tPASS\t.\n", (n - i) * 10))
                .collect()
        };
        assert_eq!(merged_files("small-300", 0, &small(130), 300), 26);
        // 33 blocks: bcftools merges the full first layer before it writes the last block.
        assert_eq!(merged_files("small-200", 0, &small(130), 200), 2);

        // Long and multiallelic alleles, INFO strings, and END (a longer varint rlen).
        let mixed: String = (0..60)
            .map(|i| {
                let pos = 1000 + (i * 7919) % 997 * 10;
                let reference = if i % 5 == 0 {
                    "A".repeat(16)
                } else {
                    "A".to_owned()
                };
                let alts = ["C", "C,G", "C,G,TT"][i % 3];
                let mut info = Vec::new();
                if i * 7 % 40 != 0 {
                    info.push(format!("S={}", "x".repeat(i * 7 % 40)));
                }
                if i % 4 == 0 {
                    info.push(format!("END={}", pos + 200));
                }
                let info = if info.is_empty() {
                    ".".to_owned()
                } else {
                    info.join(";")
                };
                format!("chr1\t{pos}\t.\t{reference}\t{alts}\t50\tPASS\t{info}\n")
            })
            .collect();
        for (max_mem, files) in [(700, 8), (1000, 6), (2000, 3)] {
            assert_eq!(
                merged_files("mixed", 0, &mixed, max_mem),
                files,
                "-m {max_mem}"
            );
        }

        // 130 samples: the sample count and the per-sample block take two-byte varints in
        // bcftools' packed form, so two records share a block from 395 bytes on (-m 439); with
        // 394 the second one does not fit and is written out with the first.
        let samples: String = (0..3)
            .map(|i| {
                format!(
                    "chr1\t{}\t.\tA\tC\t50\tPASS\t.\tDP{}\n",
                    30 - i * 10,
                    "\t7".repeat(130)
                )
            })
            .collect();
        assert_eq!(merged_files("samples-438", 130, &samples, 438), 2);
        assert_eq!(merged_files("samples-439", 130, &samples, 439), 1);
    }

    #[test]
    fn memory_beyond_ptrdiff_max_cannot_be_allocated() {
        assert_eq!(sort_memory(768_000_000).unwrap(), 691_200_000);
        assert_eq!(sort_memory(0).unwrap(), 0);
        let err = sort_memory(parse_mem_string("-5").unwrap()).unwrap_err();
        assert_eq!(
            err.to_string(),
            "Error: could not allocate 16602069666338596864 bytes of memory, try reducing --max-mem"
        );
    }

    #[test]
    fn sigpipe_handler_names_the_files_as_the_directory_does() {
        // Never dropped: the directory does not exist, and nothing is to be removed.
        let mut tmp = std::mem::ManuallyDrop::new(TempDir {
            path: PathBuf::from("/tmp/bcftools.AbC123"),
            files: Vec::new(),
            count: 0,
            removed_on_sigpipe: false,
        });
        let dir = tmp.path.as_os_str().as_encoded_bytes().to_vec();
        let mut buf = [0u8; 64];
        for i in 0..100_123 {
            let want = tmp.next_name();
            let got = tmp_on_sigpipe::file_name(&dir, i, &mut buf).expect("the name fits");
            assert_eq!(
                got.to_bytes(),
                want.as_os_str().as_encoded_bytes(),
                "file {i}"
            );
        }
        let long = [b'd'; 58];
        assert!(
            tmp_on_sigpipe::file_name(&long, 0, &mut buf).is_none(),
            "64 bytes and the NUL need more than 64"
        );
        assert!(
            tmp_on_sigpipe::file_name(&long[..57], 0, &mut buf).is_some(),
            "63 bytes and the NUL fit in 64"
        );
    }
}
