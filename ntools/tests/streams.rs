//! `cat`, `reheader`, and `quickcheck` on inputs that cannot seek (stdin and FIFOs), and on CRAM
//! files with corrupt container headers; `view` and `stats` reading stdin when given no input;
//! writes to a closed pipe.

use std::{
    ffi::CString,
    fs,
    io::{self, Read, Write},
    os::unix::{ffi::OsStrExt, process::ExitStatusExt},
    path::{Path, PathBuf},
    process::{Command, ExitStatus, Output, Stdio},
    thread,
    time::{Duration, Instant},
};

/// How long a command may run before the test fails as hung.
const TIMEOUT: Duration = Duration::from_secs(60);

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

/// How a test passes an input file to a command.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Mode {
    File,
    /// `-`, with the file's bytes written to a pipe on stdin.
    Stdin,
    /// `/dev/stdin`, on a pipe.
    DevStdin,
    /// A named pipe that a thread writes the file to.
    Fifo,
}

const MODES: [Mode; 4] = [Mode::File, Mode::Stdin, Mode::DevStdin, Mode::Fifo];

/// Runs ntools with `args`, then the input, given as `mode` says, and returns its output.
fn run(dir: &Path, args: &[&str], input: &Path, mode: Mode) -> io::Result<Output> {
    let data = fs::read(input)?;
    let fifo = dir.join(format!(
        "{}.fifo",
        input.file_stem().unwrap_or_default().display()
    ));

    let input_arg = match mode {
        Mode::File => input.to_path_buf(),
        Mode::Stdin => PathBuf::from("-"),
        Mode::DevStdin => PathBuf::from("/dev/stdin"),
        Mode::Fifo => {
            make_fifo(&fifo)?;
            fifo.clone()
        }
    };

    let mut command = Command::new(NTOOLS);
    command
        .current_dir(dir)
        .args(args)
        .arg(&input_arg)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let piped = matches!(mode, Mode::Stdin | Mode::DevStdin);
    command.stdin(if piped { Stdio::piped() } else { Stdio::null() });

    let mut child = command.spawn()?;

    let writer = match mode {
        Mode::File => None,
        Mode::Stdin | Mode::DevStdin => {
            let mut stdin = child.stdin.take().expect("stdin is piped");
            Some(thread::spawn(move || stdin.write_all(&data)))
        }
        Mode::Fifo => {
            let fifo = fifo.clone();
            Some(thread::spawn(move || fs::write(fifo, data)))
        }
    };

    let deadline = Instant::now() + TIMEOUT;

    while child.try_wait()?.is_none() {
        if Instant::now() > deadline {
            child.kill()?;
            panic!("ntools {args:?} {input_arg:?} did not finish within {TIMEOUT:?}");
        }

        thread::sleep(Duration::from_millis(10));
    }

    let output = child.wait_with_output()?;

    // A command that stops reading early closes the pipe, which is not an error here.
    if let Some(writer) = writer {
        let _ = writer.join().expect("writer thread panicked");
    }

    if mode == Mode::Fifo {
        fs::remove_file(&fifo)?;
    }

    Ok(output)
}

fn make_fifo(path: &Path) -> io::Result<()> {
    let path = CString::new(path.as_os_str().as_bytes()).map_err(io::Error::other)?;

    // SAFETY: `path` is a valid NUL-terminated string that outlives the call.
    if unsafe { libc::mkfifo(path.as_ptr(), 0o600) } == 0 {
        Ok(())
    } else {
        Err(io::Error::last_os_error())
    }
}

/// Runs ntools with `args` and nothing more, with the bytes of `input` on a pipe on stdin.
fn run_piped(dir: &Path, args: &[&str], input: &Path) -> io::Result<Output> {
    let data = fs::read(input)?;

    let mut child = Command::new(NTOOLS)
        .current_dir(dir)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;

    let mut stdin = child.stdin.take().expect("stdin is piped");
    let writer = thread::spawn(move || stdin.write_all(&data));
    let output = child.wait_with_output()?;
    let _ = writer.join().expect("writer thread panicked");

    Ok(output)
}

/// Runs ntools on files, failing the test if it fails.
fn ntools(dir: &Path, args: &[&str]) -> io::Result<()> {
    let output = Command::new(NTOOLS).current_dir(dir).args(args).output()?;
    assert!(
        output.status.success(),
        "ntools {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    Ok(())
}

/// Writes a reference, a SAM file, and BAM and CRAM 3.0 copies (with several containers) to
/// `dir`.
fn write_inputs(dir: &Path) -> io::Result<()> {
    let sequence: String = (0..1000).map(|i| ['A', 'C', 'G', 'T'][i * 7 % 4]).collect();
    fs::write(dir.join("ref.fa"), format!(">chr1\n{sequence}\n"))?;

    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    sam.push_str("@RG\tID:rg1\tSM:sample\n");

    for i in 1..=200 {
        let pos = i * 4;
        let bases = &sequence[pos - 1..pos + 19];
        sam.push_str(&format!(
            "r{i}\t0\tchr1\t{pos}\t60\t20M\t*\t0\t0\t{bases}\t{}\tRG:Z:rg1\n",
            "I".repeat(20)
        ));
    }

    fs::write(dir.join("in.sam"), sam)?;
    fs::write(
        dir.join("header.sam"),
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n",
    )?;

    ntools(dir, &["faidx", "ref.fa"])?;
    ntools(dir, &["view", "--no-PG", "-b", "-o", "in.bam", "in.sam"])?;
    ntools(
        dir,
        &[
            "view",
            "--no-PG",
            "-T",
            "ref.fa",
            "-O",
            "cram,version=3.0,seqs_per_slice=50",
            "-o",
            "in.cram",
            "in.sam",
        ],
    )
}

fn temp_dir() -> io::Result<tempfile::TempDir> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path())?;
    Ok(dir)
}

#[test]
fn test_cat_and_reheader_read_streams() -> io::Result<()> {
    let dir = temp_dir()?;
    let dir = dir.path();

    let commands: [&[&str]; 3] = [
        &["cat", "--no-PG"],
        &["reheader", "-P", "header.sam"],
        &["reheader", "-P", "-c", "sed s/sample/renamed/"],
    ];

    for args in commands {
        for input in ["in.bam", "in.cram"] {
            let expected = run(dir, args, &dir.join(input), Mode::File)?;
            assert!(expected.status.success(), "{args:?} {input}");
            assert!(!expected.stdout.is_empty(), "{args:?} {input}");

            for mode in MODES {
                let actual = run(dir, args, &dir.join(input), mode)?;
                assert!(
                    actual.status.success(),
                    "{args:?} {input} {mode:?}: {}",
                    String::from_utf8_lossy(&actual.stderr)
                );
                assert_eq!(actual.stdout, expected.stdout, "{args:?} {input} {mode:?}");
            }
        }
    }

    Ok(())
}

/// With no input argument, `view` and `stats` read stdin, as samtools does when stdin is
/// not a terminal (`bwa mem ... | ntools view -b -o out.bam`).
#[test]
fn test_view_and_stats_read_stdin_without_input() -> io::Result<()> {
    let dir = temp_dir()?;
    let dir = dir.path();

    for input in ["in.sam", "in.bam", "in.cram"] {
        let args: &[&str] = &["view", "-c", "-T", "ref.fa"];
        let output = run_piped(dir, args, &dir.join(input))?;
        assert!(
            output.status.success(),
            "{input}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, b"200\n", "{input}");

        let args: &[&str] = &["view", "--no-PG", "-b", "-T", "ref.fa", "-o", "piped.bam"];
        let output = run_piped(dir, args, &dir.join(input))?;
        assert!(output.status.success(), "{input}");
        let args: &[&str] = &["view", "--no-PG", "-T", "ref.fa"];
        let expected = run(dir, args, &dir.join(input), Mode::File)?;
        let actual = run(dir, args, &dir.join("piped.bam"), Mode::File)?;
        assert_eq!(actual.stdout, expected.stdout, "{input}");

        let output = run_piped(dir, &["stats", "--reference", "ref.fa"], &dir.join(input))?;
        assert!(output.status.success(), "{input}");
        let stats = String::from_utf8_lossy(&output.stdout);
        assert!(
            stats.contains("SN\traw total sequences:\t200\t"),
            "{input}: {stats}"
        );
    }

    Ok(())
}

/// A BAM cut inside a BGZF block is an error, not a corrupt output with a fresh EOF
/// marker; one cut between blocks is copied with samtools' truncation warning.
#[test]
fn test_cat_and_reheader_detect_truncated_bam() -> io::Result<()> {
    let dir = temp_dir()?;
    let dir = dir.path();

    let bam = fs::read(dir.join("in.bam"))?;
    let eof_start = bam.len() - 28;
    // Inside the last block of records, and inside the EOF marker's block header.
    fs::write(dir.join("cut-block.bam"), &bam[..eof_start - 20])?;
    fs::write(dir.join("cut-eof.bam"), &bam[..eof_start + 10])?;
    fs::write(dir.join("no-eof.bam"), &bam[..eof_start])?;

    let commands: [&[&str]; 2] = [&["cat", "--no-PG"], &["reheader", "-P", "header.sam"]];

    for args in commands {
        for mode in [Mode::File, Mode::Stdin] {
            for input in ["cut-block.bam", "cut-eof.bam"] {
                let output = run(dir, args, &dir.join(input), mode)?;
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{args:?} {input} {mode:?}");
                assert!(
                    stderr.contains("truncated"),
                    "{args:?} {input} {mode:?}: {stderr}"
                );
            }

            let expected = run(dir, args, &dir.join("in.bam"), mode)?;
            let output = run(dir, args, &dir.join("no-eof.bam"), mode)?;
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(output.status.success(), "{args:?} {mode:?}: {stderr}");
            assert!(
                stderr.contains("EOF marker is absent"),
                "{args:?} {mode:?}: {stderr}"
            );
            assert_eq!(output.stdout, expected.stdout, "{args:?} {mode:?}");
        }
    }

    Ok(())
}

#[test]
fn test_quickcheck_reads_streams() -> io::Result<()> {
    let dir = temp_dir()?;
    let dir = dir.path();

    for input in ["in.sam", "in.bam", "in.cram"] {
        for mode in MODES {
            let output = run(dir, &["quickcheck", "-vvv"], &dir.join(input), mode)?;
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(output.status.success(), "{input} {mode:?}: {stderr}");
            assert!(
                stderr.contains("has 1 targets in header."),
                "{input} {mode:?}: {stderr}"
            );

            // As in HTSlib, a pipe's EOF block cannot be checked, and is not a failure.
            let eof_message = match (input, mode) {
                ("in.sam", _) => "because its filetype does not contain one.",
                (_, Mode::File) => "has good EOF block.",
                _ => "as it is not seekable.",
            };
            assert!(stderr.contains(eof_message), "{input} {mode:?}: {stderr}");
        }
    }

    // An empty pipe is not sequence data (4), as with an empty file.
    fs::write(dir.join("empty.bam"), "")?;

    for mode in MODES {
        let output = run(dir, &["quickcheck"], &dir.join("empty.bam"), mode)?;
        assert_eq!(output.status.code(), Some(4), "{mode:?}");
    }

    Ok(())
}

/// Returns `src` with the low bit of the byte at `offset` flipped.
fn flip_bit(src: &[u8], offset: usize) -> Vec<u8> {
    let mut dst = src.to_vec();
    dst[offset] ^= 1;
    dst
}

/// The offset of the first data container: after the 26-byte file definition and the header
/// container.
fn first_data_container(cram: &[u8]) -> usize {
    /// Reads an ITF-8 integer (LTF-8 too, for the small values here).
    fn read_itf8(src: &[u8], offset: &mut usize) -> usize {
        let b0 = src[*offset];
        let extra = usize::try_from(b0.leading_ones()).expect("small").min(4);
        let mut n = usize::from(b0 & (0xff >> (extra + 1)));

        for &b in &src[*offset + 1..*offset + 1 + extra] {
            n = (n << 8) | usize::from(b);
        }

        *offset += 1 + extra;
        n
    }

    let len = u32::from_le_bytes(cram[26..30].try_into().expect("4 bytes"));
    let mut offset = 30;

    // Reference sequence ID, start, span, record count, record counter, base count, and block
    // count, then the landmarks and the CRC32.
    for _ in 0..7 {
        read_itf8(cram, &mut offset);
    }

    for _ in 0..read_itf8(cram, &mut offset) {
        read_itf8(cram, &mut offset);
    }

    offset + 4 + usize::try_from(len).expect("valid length")
}

#[test]
fn test_cat_and_reheader_check_container_header_crc32s() -> io::Result<()> {
    let dir = temp_dir()?;
    let dir = dir.path();

    let cram = fs::read(dir.join("in.cram"))?;
    let data_container = first_data_container(&cram);

    // The header container's span, and the first data container's reference sequence ID.
    fs::write(dir.join("bad-header.cram"), flip_bit(&cram, 32))?;
    fs::write(
        dir.join("bad-data.cram"),
        flip_bit(&cram, data_container + 4),
    )?;
    // Cut inside the first data container.
    fs::write(dir.join("truncated.cram"), &cram[..data_container + 40])?;

    let commands: [&[&str]; 2] = [&["cat", "--no-PG"], &["reheader", "-P", "header.sam"]];

    for args in commands {
        for (input, message) in [
            ("bad-header.cram", "container header CRC32 mismatch"),
            ("bad-data.cram", "container header CRC32 mismatch"),
            ("truncated.cram", "truncated container"),
        ] {
            for mode in [Mode::File, Mode::Stdin] {
                let output = run(dir, args, &dir.join(input), mode)?;
                let stderr = String::from_utf8_lossy(&output.stderr);
                assert_eq!(output.status.code(), Some(1), "{args:?} {input} {mode:?}");
                assert!(
                    stderr.contains(message),
                    "{args:?} {input} {mode:?}: {stderr}"
                );
            }
        }
    }

    let output = run(
        dir,
        &["quickcheck"],
        &dir.join("bad-header.cram"),
        Mode::File,
    )?;
    assert_eq!(
        output.status.code(),
        Some(2),
        "a CRAM file that cannot be opened"
    );

    Ok(())
}

/// Runs ntools with `args` in `dir`, with standard output (or, with `stderr_closed`, standard
/// error) a pipe whose reading end is closed. Returns the exit status and what it wrote to the
/// other stream.
fn run_into_closed_pipe(
    dir: &Path,
    args: &[&str],
    stderr_closed: bool,
) -> io::Result<(ExitStatus, String)> {
    let (reader, writer) = io::pipe()?;
    drop(reader);
    let (mut other, other_writer) = io::pipe()?;
    let mut command = Command::new(NTOOLS);
    command.args(args).current_dir(dir).stdin(Stdio::null());
    if stderr_closed {
        command.stderr(writer).stdout(other_writer);
    } else {
        command.stdout(writer).stderr(other_writer);
    }
    let mut child = command.spawn()?;
    // Close this process's copies of the write ends, so that reading stops when ntools ends.
    drop(command);
    let mut text = String::new();
    other.read_to_string(&mut text)?;
    Ok((child.wait()?, text))
}

#[test]
fn test_closed_output_pipe_kills_with_sigpipe() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();

    // Enough records for `sort -n -m 1M` and `collate -r 1000` to spill to temporary files.
    let mut sam = String::from("@SQ\tSN:c\tLN:1000000\n");
    for i in 0..30_000u32 {
        let pos = i * 30 + 1;
        sam.push_str(&format!(
            "r{i}\t0\tc\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
        ));
    }
    fs::write(dir.join("in.sam"), sam)?;

    let prefix = format!("{}/tmp", dir.display());
    let cases: &[&[&str]] = &[
        &["view", "in.sam"],
        &["view", "-b", "-@", "2", "in.sam"],
        &[
            "sort", "-n", "-m", "1M", "-T", &prefix, "-O", "sam", "in.sam",
        ],
        &["collate", "-r", "1000", "-T", &prefix, "-O", "in.sam"],
        &["depth", "in.sam"],
        &["flagstat", "in.sam"],
    ];

    for args in cases {
        let (status, err) = run_into_closed_pipe(dir, args, false)?;
        // samtools is killed by SIGPIPE, as ntools is.
        assert_eq!(status.signal(), Some(libc::SIGPIPE), "{args:?}: {status}");
        assert!(!err.contains("ntools:"), "{args:?} reports no error: {err}");
        let left: Vec<_> = fs::read_dir(dir)?
            .map(|e| e.map(|e| e.file_name()))
            .collect::<io::Result<_>>()?;
        assert_eq!(left, ["in.sam"], "{args:?} removes its temporary files");
    }

    let (status, _) = run_into_closed_pipe(dir, &["view", "nosuch.bam"], true)?;
    assert_eq!(
        status.signal(),
        Some(libc::SIGPIPE),
        "an error on a closed standard error: {status}"
    );

    Ok(())
}

#[test]
fn test_output_naming_an_input_is_refused() -> io::Result<()> {
    let tmp = temp_dir()?;
    let dir = tmp.path();
    let cases: &[&[&str]] = &[
        &["view", "-b", "-o", "in.bam", "in.bam"],
        &["view", "-o", "in.sam", "in.sam"],
        &["sort", "-o", "in.bam", "in.bam"],
    ];
    for args in cases {
        let input = args[args.len() - 1];
        let before = fs::read(dir.join(input))?;
        let out = Command::new(NTOOLS).current_dir(dir).args(*args).output()?;
        assert!(!out.status.success(), "{args:?}: {out:?}");
        assert!(
            String::from_utf8_lossy(&out.stderr).contains("also an input"),
            "{args:?}: {out:?}"
        );
        assert_eq!(
            fs::read(dir.join(input))?,
            before,
            "{args:?} leaves the input"
        );
    }
    Ok(())
}
