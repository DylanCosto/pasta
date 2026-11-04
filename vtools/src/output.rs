//! bcftools' shared output options: `-o`, `-O u|b|v|z[0-9]`, `-l`, `--threads`, and
//! `-W/--write-index` (bcftools' `set_wmode`, `write_index_parse`, and `init_index2`).

use std::path::{Path, PathBuf};

use anyhow::bail;
use pasta::variant::io::{IndexRequest, OutputType, WriteOptions};

/// bcftools' output type bits (`FT_*`).
const FT_GZ: u32 = 1;
const FT_VCF: u32 = 2;
const FT_BCF: u32 = 4;

/// Separates the data file name from the index file name in an output name (HTSlib's
/// `HTS_IDX_DELIM`).
const IDX_DELIM: &str = "##idx##";

/// The `min_shift` of a CSI index that bcftools writes (a TBI index has 0).
const CSI_MIN_SHIFT: u32 = 14;

/// One `-O` or `-l` occurrence, in command-line order.
#[derive(Clone, Debug)]
pub enum OutputSetting {
    /// `-O TYPE`.
    Type(String),
    /// `-l LEVEL`.
    Level(String),
}

/// Resolves the output type and level from the `-O`/`-l` occurrences in order and the output
/// name, as bcftools' option loop and `set_wmode` do.
pub fn resolve(
    settings: &[OutputSetting],
    output: Option<&Path>,
) -> anyhow::Result<(OutputType, Option<u8>)> {
    let (file_type, clevel) = parse_settings(settings)?;
    set_wmode(file_type, output, clevel)
}

/// Returns the output type `-O`/`-l` ask for (bcftools' `output_type`), before an output
/// name's extension applies and without the level check that opening an output makes.
pub fn requested_type(settings: &[OutputSetting]) -> anyhow::Result<OutputType> {
    parse_settings(settings).map(|(file_type, _)| type_of(file_type))
}

/// Parses the `-O`/`-l` occurrences in order into bcftools' `FT_*` bits and level.
fn parse_settings(settings: &[OutputSetting]) -> anyhow::Result<(u32, i64)> {
    let mut file_type = FT_VCF;
    let mut clevel: i64 = -1;
    for s in settings {
        match s {
            OutputSetting::Type(t) => {
                let b = t.as_bytes();
                match b.first() {
                    Some(b'b') => file_type = FT_BCF | FT_GZ,
                    Some(b'u') => file_type = FT_BCF,
                    Some(b'z') => file_type = FT_VCF | FT_GZ,
                    Some(b'v') => file_type = FT_VCF,
                    // An empty type is rejected here: bcftools would go on to read past the
                    // end of the argument.
                    _ => match strtol_whole_int(b) {
                        Some(n) if !t.is_empty() && (0..=9).contains(&n) => clevel = n,
                        _ => bail!("The output type \"{t}\" not recognised"),
                    },
                }
                if b.len() > 1 {
                    // Bytes, as the first character may be wider than one byte.
                    let rest = &b[1..];
                    match strtol_whole_int(rest) {
                        Some(n) if (0..=9).contains(&n) => clevel = n,
                        _ => bail!(
                            "Could not parse argument: --compression-level {}",
                            String::from_utf8_lossy(rest)
                        ),
                    }
                }
            }
            OutputSetting::Level(l) => {
                let Some(n) = strtol_whole_int(l.as_bytes()) else {
                    bail!("Could not parse argument: --compression-level {l}");
                };
                clevel = n;
                file_type |= FT_GZ;
            }
        }
    }
    Ok((file_type, clevel))
}

/// C's `strtol(s, &end, 10)` into an `int`, when the whole text is a number (an empty text is
/// 0); a larger number is truncated as C's conversion does.
fn strtol_whole_int(s: &[u8]) -> Option<i64> {
    let (v, end) = pasta::variant::number::strtol(s, 0);
    (end == s.len()).then_some(i64::from(v as i32))
}

/// The output type of bcftools' `FT_*` bits (`hts_bcf_wmode`).
fn type_of(file_type: u32) -> OutputType {
    if file_type == FT_BCF {
        OutputType::BcfUncompressed
    } else if file_type & FT_BCF != 0 {
        OutputType::Bcf
    } else if file_type & FT_GZ != 0 {
        OutputType::VcfBgzf
    } else {
        OutputType::Vcf
    }
}

/// bcftools' `set_wmode`: the output name's extension overrides `-O`.
fn set_wmode(
    file_type: u32,
    output: Option<&Path>,
    clevel: i64,
) -> anyhow::Result<(OutputType, Option<u8>)> {
    let name = output.map(|p| p.to_string_lossy().into_owned());
    let base = name.as_deref().map(|n| {
        n.find("##idx##")
            .map_or(n, |i| &n[..i])
            .to_ascii_lowercase()
    });
    let ty = match base.as_deref() {
        Some(b) if b.ends_with(".bcf") => type_of(if file_type & FT_BCF != 0 {
            file_type
        } else {
            FT_BCF | FT_GZ
        }),
        Some(b) if b.ends_with(".vcf") => OutputType::Vcf,
        Some(b) if b.ends_with(".vcf.gz") || b.ends_with(".vcf.bgz") => OutputType::VcfBgzf,
        _ => type_of(file_type),
    };
    if (0..=9).contains(&clevel) {
        if ty == OutputType::BcfUncompressed {
            bail!(
                "Error: compression level ({clevel}) cannot be set on uncompressed streams ({})",
                name.as_deref().unwrap_or("(null)")
            );
        }
        // HTSlib ignores a level on uncompressed VCF ("w5" opens plain text).
        return Ok((ty, Some(clevel as u8)));
    }
    Ok((ty, None))
}

/// The most threads a `--threads` count starts; a larger count is cut to it. HTSlib has no limit
/// (`hts_tpool_init` starts every thread asked for, so bcftools hangs or runs out of memory on a
/// huge count), and vtools' worker pools, batches in flight, and mpileup's thread budget grow
/// with the count. Output does not depend on the count.
pub const MAX_THREADS: i32 = 1024;

/// Parses `--threads` as bcftools does (`int n_threads = strtol(optarg, 0, 0)`; 0 without the
/// option), cut to [`max_threads`]. A negative count, which HTSlib cannot create, is returned as
/// the error for the command to report or ignore as bcftools does.
pub fn threads(value: Option<&str>) -> Result<usize, i64> {
    let n = value.map_or(0, thread_count);
    usize::try_from(n).map_err(|_| i64::from(n))
}

/// A `--threads` value read as bcftools reads it (`(int) strtol(s, NULL, 0)`), with a positive
/// count cut to [`max_threads`]; zero and negative counts are kept for the command to handle as
/// bcftools does.
pub fn thread_count(s: &str) -> i32 {
    strtol_int(s).min(max_threads())
}

/// Cuts a `--threads` count that a command parsed on its own (strictly, or without the `int`
/// truncation) to [`max_threads`]; zero and negative counts are kept.
pub fn cap_threads(n: i64) -> i64 {
    n.min(i64::from(max_threads()))
}

/// The cap of a `--threads` count: [`MAX_THREADS`], or fewer under an address-space limit
/// (`ulimit -v`), where each thread is given 64 MB of it (see [`crate::limits`]).
pub fn max_threads() -> i32 {
    crate::limits::max_threads().map_or(MAX_THREADS, |n| {
        i32::try_from(n).unwrap_or(i32::MAX).min(MAX_THREADS)
    })
}

/// Whether a thread can be started, when `threads` asks for any. Worker pools start fewer threads
/// when the system cannot start them all, and work without threads when it cannot start one; the
/// commands whose bcftools thread pool fails to start then (`view`, `annotate`, `norm`, ...)
/// check this first to fail as bcftools does.
pub fn threads_start(threads: usize) -> bool {
    threads == 0 || pasta::util::threads::spawn(|| {}).is_ok_and(|t| t.join().is_ok())
}

/// C's `(int) strtol(s, NULL, 0)`, as bcftools reads its `--threads`: an optional sign, then a
/// decimal, `0x` hex, or `0` octal number, clamped to `long` and then truncated to `int` (so
/// `4294967297` is 1); 0 when no number parses.
fn strtol_int(s: &str) -> i32 {
    let s = s.trim_start_matches([' ', '\t', '\n', '\x0b', '\x0c', '\r']);
    let (negative, digits) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };
    // A second sign or a space after the sign is not a number.
    if !digits.as_bytes().first().is_some_and(u8::is_ascii_digit) {
        return 0;
    }
    // strtoull saturates at u64::MAX; strtol then clamps to the long range.
    let magnitude = crate::commands::head::strtoull0(digits);
    let long = if negative {
        0i64.checked_sub_unsigned(magnitude).unwrap_or(i64::MIN)
    } else {
        i64::try_from(magnitude).unwrap_or(i64::MAX)
    };
    // The conversion to `int` keeps the low 32 bits, as C compilers do.
    long as i32
}

/// Builds the writer options.
pub fn write_options(ty: OutputType, level: Option<u8>, threads: usize) -> WriteOptions {
    WriteOptions {
        output_type: ty,
        level,
        threads,
    }
}

/// The index format of `--write-index[=FMT]`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WriteIndex {
    /// CSI (the default).
    Csi,
    /// TBI (VCF only; BCF gets CSI).
    Tbi,
}

/// Parses `--write-index`'s optional value (bcftools' `write_index_parse`).
pub fn parse_write_index(value: &str) -> anyhow::Result<WriteIndex> {
    match value {
        "" | "csi" | "=csi" => Ok(WriteIndex::Csi),
        "tbi" | "=tbi" => Ok(WriteIndex::Tbi),
        _ => bail!("Unsupported index format '{value}'"),
    }
}

/// Plans the index file for `--write-index` (bcftools' `init_index2`): `OUT.csi` or
/// `OUT.tbi`, or the name after `##idx##`.
pub fn index_request(
    fmt: WriteIndex,
    ty: OutputType,
    output: Option<&Path>,
) -> Option<IndexRequest> {
    let (mut min_shift, suffix) = if fmt == WriteIndex::Tbi && ty.is_vcf() {
        (0, "tbi")
    } else {
        (CSI_MIN_SHIFT, "csi")
    };
    let name = output?.to_string_lossy().into_owned();
    if name.is_empty() || name == "-" {
        return None;
    }
    let path = match name.find(IDX_DELIM) {
        Some(i) => {
            let idx = name[i + IDX_DELIM.len()..].to_owned();
            if idx.ends_with(".tbi") {
                min_shift = 0;
            }
            PathBuf::from(idx)
        }
        None => PathBuf::from(format!("{name}.{suffix}")),
    };
    Some(IndexRequest { path, min_shift })
}

/// Returns the data file name of an output name that may carry `##idx##INDEX`.
pub fn data_path(output: &Path) -> PathBuf {
    let name = output.to_string_lossy();
    PathBuf::from(name.find(IDX_DELIM).map_or(&*name, |i| &name[..i]))
}

/// The error to return for a failed write: the write's own error when the reader of the output
/// went away (`| head`), for `main` to die of SIGPIPE, or `message` (bcftools' error) otherwise.
pub fn write_failure<E>(e: E, message: impl FnOnce() -> anyhow::Error) -> anyhow::Error
where
    E: std::error::Error + Send + Sync + 'static,
{
    let mut cause: Option<&(dyn std::error::Error + 'static)> = Some(&e);
    while let Some(c) = cause {
        if c.downcast_ref::<std::io::Error>()
            .is_some_and(|io| io.kind() == std::io::ErrorKind::BrokenPipe)
        {
            return anyhow::Error::new(e);
        }
        cause = c.source();
    }
    message()
}

/// A failed read of the input (a truncated or corrupt BGZF block, a malformed BCF record) that
/// bcftools takes for the end of the input, writing out the records before it with status 0.
/// The command closes its output as bcftools does, so the records before the damaged block
/// are written in full, and then fails with this message.
#[derive(Debug)]
pub struct InputReadFailure(pub String);

impl std::fmt::Display for InputReadFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for InputReadFailure {}

/// Standard output for a command's output, line-buffered as [`std::io::stdout`] is, whose
/// writes report every error. Rust's standard output takes `EBADF` for success, so a closed
/// standard output (see `closed_stdout`) would swallow the output with status 0.
pub fn stdout() -> std::io::LineWriter<Stdout> {
    std::io::LineWriter::new(Stdout::new())
}

/// Descriptor 1, written to directly; see [`stdout`].
pub struct Stdout(std::mem::ManuallyDrop<std::fs::File>);

impl Stdout {
    fn new() -> Self {
        use std::os::fd::FromRawFd;
        // SAFETY: descriptor 1 is open for the life of the process (Rust's start-up code or
        // `closed_stdout` fills it when it is closed at start), and `ManuallyDrop` never closes
        // it, so the `File` only borrows it.
        Self(std::mem::ManuallyDrop::new(unsafe {
            std::fs::File::from_raw_fd(libc::STDOUT_FILENO)
        }))
    }
}

impl std::io::Write for Stdout {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.write(buf)
    }

    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        self.0.write_vectored(bufs)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

/// Logs HTSlib's errors for a file that cannot be opened (`hts_open`): a directory opens, but
/// detecting its format fails first (`hts_hopen`).
pub fn log_open_failure(name: &str, e: &std::io::Error) {
    if e.kind() == std::io::ErrorKind::IsADirectory {
        pasta::variant::log::error("hts_hopen", format_args!("Failed to open file {name}"));
    }
    pasta::variant::log::error(
        "hts_open_format",
        format_args!(
            "Failed to open file \"{name}\" : {}",
            crate::raw_input::strerror(e)
        ),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strtol_int_is_c() {
        let cases = [
            ("4", 4),
            (" +4x", 4),
            ("0x10", 16),
            ("010", 8),
            ("-1", -1),
            ("-0x1", -1),
            ("4294967297", 1),
            ("-4294967295", 1),
            ("99999999999999999999", -1),
            ("-99999999999999999999", 0),
            ("", 0),
            ("x", 0),
            ("--1", 0),
            ("- 1", 0),
        ];
        for (s, want) in cases {
            assert_eq!(strtol_int(s), want, "strtol_int({s:?})");
        }
        assert_eq!(threads(None), Ok(0));
        assert_eq!(threads(Some("4294967298")), Ok(2));
        assert_eq!(threads(Some("-2")), Err(-2));
    }

    #[test]
    fn huge_thread_counts_are_capped() {
        // 99999999999 is 1215752191 as a C int.
        for s in ["99999999999", "1025", "0x7fffffff"] {
            assert_eq!(thread_count(s), max_threads(), "thread_count({s:?})");
            assert_eq!(
                threads(Some(s)),
                Ok(max_threads() as usize),
                "threads({s:?})"
            );
        }
        assert!(max_threads() <= MAX_THREADS);
        assert_eq!(thread_count("1024"), max_threads());
        assert_eq!(thread_count("4"), 4);
        assert_eq!(thread_count("0"), 0);
        assert_eq!(thread_count("x"), 0);
        assert_eq!(thread_count("-1"), -1);
        // Past i32::MAX, C's truncation can give a negative count, which stays negative.
        assert_eq!(thread_count("4294967295"), -1);
        assert_eq!(cap_threads(99_999_999_999), i64::from(max_threads()));
        assert_eq!(cap_threads(3), 3);
        assert_eq!(cap_threads(-5), -5);
    }
}
