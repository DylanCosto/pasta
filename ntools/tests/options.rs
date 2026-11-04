//! Command lines that samtools accepts or rejects: where `sort -T` puts its spill files, the
//! smallest `-m`, repeated options.

use std::{
    ffi::CString,
    fs, io,
    os::unix::ffi::OsStrExt,
    path::Path,
    process::{Command, Output},
    time::{Duration, SystemTime},
};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

fn run(dir: &Path, args: &[&str]) -> io::Result<Output> {
    Command::new(NTOOLS).current_dir(dir).args(args).output()
}

/// Runs ntools in `dir`, failing the test if it fails, and returns its standard output.
fn ntools(dir: &Path, args: &[&str]) -> io::Result<String> {
    let output = run(dir, args)?;
    assert!(
        output.status.success(),
        "ntools {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).map_err(io::Error::other)
}

/// Writes `in.bam`: `n` unsorted, unpaired records on one reference sequence.
fn write_bam(dir: &Path, n: usize) -> io::Result<()> {
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n");

    for i in 0..n {
        let pos = (i * 7919) % 99000 + 1;
        sam.push_str(&format!(
            "r{i}\t0\tchr1\t{pos}\t60\t50M\t*\t0\t0\t{}\t{}\n",
            "ACGT".repeat(13).get(..50).unwrap_or_default(),
            "I".repeat(50)
        ));
    }

    fs::write(dir.join("in.sam"), sam)?;
    ntools(dir, &["view", "--no-PG", "-b", "-o", "in.bam", "in.sam"])?;
    Ok(())
}

fn set_mtime(path: &Path, time: SystemTime) -> io::Result<()> {
    fs::File::open(path)?.set_modified(time)
}

fn mtime(path: &Path) -> io::Result<SystemTime> {
    fs::metadata(path)?.modified()
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

/// `sort -T DIR` with an existing directory spills inside DIR, as samtools sort does,
/// not in a `DIR.XXXXXX` directory beside it. Creating and removing an entry changes a
/// directory's modification time, which shows where the spill directory was made.
#[test]
fn test_sort_spills_inside_an_existing_temp_dir() -> io::Result<()> {
    let base = tempfile::tempdir()?;
    let dir = base.path().join("work");
    let tdir = dir.join("tdir");
    fs::create_dir_all(&tdir)?;

    // About 2.4 MB of records: spills with -m 1M.
    write_bam(&dir, 20_000)?;
    let expected = ntools(&dir, &["sort", "--no-PG", "-O", "sam", "in.bam"])?;

    for prefix in ["tdir", "tdir/"] {
        let old = SystemTime::UNIX_EPOCH + Duration::from_secs(1_000_000_000);
        set_mtime(&dir, old)?;
        set_mtime(&tdir, old)?;

        let args = [
            "sort", "--no-PG", "-O", "sam", "-m", "1M", "-T", prefix, "in.bam",
        ];
        assert_eq!(ntools(&dir, &args)?, expected, "-T {prefix}");

        assert_ne!(mtime(&tdir)?, old, "-T {prefix}: nothing was made in tdir");
        assert_eq!(
            mtime(&dir)?,
            old,
            "-T {prefix}: something was made beside tdir"
        );
        assert_eq!(fs::read_dir(&tdir)?.count(), 0, "-T {prefix}");
    }

    // A prefix that is not a directory names the spill directory, in the prefix's directory.
    let old = SystemTime::UNIX_EPOCH + Duration::from_secs(1_000_000_000);
    set_mtime(&tdir, old)?;
    let args = [
        "sort", "--no-PG", "-O", "sam", "-m", "1M", "-T", "tdir/pre", "in.bam",
    ];
    assert_eq!(ntools(&dir, &args)?, expected);
    assert_ne!(mtime(&tdir)?, old);

    // A FIFO named by -T is not a directory.
    make_fifo(&dir.join("fifo"))?;
    let args = [
        "sort", "--no-PG", "-O", "sam", "-m", "1M", "-T", "fifo", "in.bam",
    ];
    assert_eq!(ntools(&dir, &args)?, expected);

    Ok(())
}

/// `sort` and `collate` reject `-m` below samtools' 1M minimum instead of spilling a file
/// every few records.
#[test]
fn test_sort_and_collate_reject_memory_below_1m() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();
    write_bam(dir, 100)?;

    for command in ["sort", "collate"] {
        for (memory, shown) in [("1", "1"), ("2", "2"), ("1023K", "1023K"), ("0", "0")] {
            let output = run(dir, &[command, "-m", memory, "-o", "out.bam", "in.bam"])?;
            assert_eq!(output.status.code(), Some(1), "{command} -m {memory}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(
                stderr.contains(&format!(
                    "-m setting ({shown} bytes) is less than the minimum required (1M)"
                )),
                "{command} -m {memory}: {stderr}"
            );
            assert!(!dir.join("out.bam").exists(), "{command} -m {memory}");
        }

        for memory in ["1M", "1024K", "1048576"] {
            ntools(dir, &[command, "-m", memory, "-o", "out.bam", "in.bam"])?;
            fs::remove_file(dir.join("out.bam"))?;
        }
    }

    Ok(())
}

/// Writes `flags.bam`: sorted paired records with every combination of a few flags and mapping
/// qualities, and `names.txt` and `names2.txt`, which list some of their names.
fn write_flags_bam(dir: &Path) -> io::Result<()> {
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n");
    let mut i = 0;

    for flags in [0x41, 0x51, 0x81, 0x91, 0x141, 0x841, 0x481, 0x291] {
        for mapq in [0, 10, 20, 30, 60] {
            i += 1;
            sam.push_str(&format!(
                "r{i}\t{flags}\tchr1\t{}\t{mapq}\t4M\t=\t1\t0\tACGT\tIIII\n",
                i * 10
            ));
        }
    }

    fs::write(dir.join("flags.sam"), sam)?;
    fs::write(dir.join("names.txt"), "r1\nr2\nr3\n")?;
    fs::write(dir.join("names2.txt"), "r3\nr10\n")?;
    ntools(
        dir,
        &["view", "--no-PG", "-b", "-o", "flags.bam", "flags.sam"],
    )?;
    ntools(dir, &["index", "flags.bam"])?;
    Ok(())
}

/// Repeated `-f`/`-F`/`-G` (and the other flag options) are ORed, as in samtools, and
/// other repeated options take their last value, rather than being clap errors.
#[test]
fn test_repeated_options() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();
    write_flags_bam(dir)?;

    let count = |args: &[&str]| -> io::Result<String> {
        let mut args = args.to_vec();
        args.push("flags.bam");
        ntools(dir, &args)
    };

    let cases: [(&[&str], &[&str]); 8] = [
        (
            &["view", "-c", "-F", "16", "-F", "64"],
            &["view", "-c", "-F", "80"],
        ),
        (
            &["view", "-c", "-f", "64", "-f", "16"],
            &["view", "-c", "-f", "80"],
        ),
        (
            &["view", "-c", "-G", "16", "-G", "64"],
            &["view", "-c", "-G", "80"],
        ),
        (
            &["view", "-c", "--rf", "16", "--rf", "0x400"],
            &["view", "-c", "--rf", "0x410"],
        ),
        (
            &["view", "-c", "-q", "10", "-q", "30"],
            &["view", "-c", "-q", "30"],
        ),
        (
            &["view", "-c", "-e", "mapq > 0", "-e", "mapq > 20"],
            &["view", "-c", "-q", "21"],
        ),
        (
            &["depth", "-a", "-q", "10", "-q", "30", "-Q", "0", "-Q", "20"],
            &["depth", "-a", "-Q", "20"],
        ),
        (
            &[
                "depth", "-G", "16", "-G", "0x40", "-g", "DUP", "-g", "QCFAIL",
            ],
            &["depth", "-G", "0x50", "-g", "0x600"],
        ),
    ];

    for (repeated, single) in cases {
        assert_eq!(count(repeated)?, count(single)?, "{repeated:?}");
    }

    // Read names from every -N file: r1, r2, r3 and r10.
    assert_eq!(
        count(&["view", "-c", "-N", "names.txt", "-N", "names2.txt"])?,
        "4\n"
    );

    // The values differ between the repeated and single options, so the comparisons above
    // check something.
    assert_ne!(
        count(&["view", "-c", "-F", "16"])?,
        count(&["view", "-c", "-F", "80"])?
    );
    assert_ne!(
        count(&["view", "-c", "-q", "10"])?,
        count(&["view", "-c", "-q", "30"])?
    );

    // A repeated flag without a value is no error either.
    assert_eq!(count(&["view", "-c", "-c"])?, "40\n");

    // View edits flags with the ORed values too.
    let edited = ntools(
        dir,
        &[
            "view",
            "--add-flags",
            "0x400",
            "--add-flags",
            "0x200",
            "flags.bam",
            "chr1:1-10",
        ],
    )?;
    assert_eq!(
        edited.split('\t').nth(1),
        Some(&*format!("{}", 0x41 | 0x600))
    );

    Ok(())
}

/// `bam2fq` is samtools' alias for `fastq`.
#[test]
fn test_bam2fq_is_an_alias_of_fastq() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    write_bam(dir.path(), 3)?;

    let expected = ntools(dir.path(), &["fastq", "in.bam"])?;
    assert!(!expected.is_empty());
    assert_eq!(ntools(dir.path(), &["bam2fq", "in.bam"])?, expected);

    // The command list shows the alias.
    let help = ntools(dir.path(), &["--help"])?;
    assert!(help.contains("[alias: bam2fq]"), "{help}");

    Ok(())
}

/// Help wraps to the terminal, at most 80 columns, however wide the terminal is.
#[test]
fn test_help_wraps_to_80_columns() -> io::Result<()> {
    let help = |args: &[&str]| -> io::Result<String> {
        let output = Command::new(NTOOLS)
            .args(args)
            .env("COLUMNS", "300")
            .output()?;
        String::from_utf8(output.stdout).map_err(io::Error::other)
    };

    let top = help(&["--help"])?;
    let commands = top
        .lines()
        .skip_while(|line| !line.starts_with("Commands:"))
        .skip(1)
        .take_while(|line| !line.is_empty())
        .filter_map(|line| line.split_whitespace().next())
        // bgzip prints HTSlib's bgzip usage as is.
        .filter(|&name| name != "bgzip" && name != "help");

    for name in commands {
        let text = help(&[name, "--help"])?;
        assert!(text.contains("Usage:"), "{name}: {text}");

        for line in text.lines() {
            assert!(line.chars().count() <= 80, "{name}: {line}");
        }
    }

    Ok(())
}

/// Messages that tell the user to run a command name ntools's, not samtools'.
#[test]
fn test_messages_name_ntools_commands() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    fs::write(
        dir.path().join("in.sam"),
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100\n\
         r\t99\tc\t1\t60\t4M\t=\t10\t13\tACGT\tIIII\n\
         r\t147\tc\t10\t60\t4M\t=\t1\t-13\tACGT\tIIII\n",
    )?;

    let output = run(dir.path(), &["markdup", "in.sam", "out.bam"])?;
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "ntools: [markdup] error, no MC tag. Please run ntools fixmate on file first.\n"
    );

    Ok(())
}

/// `--version` names the samtools version ntools is compatible with, after `ntools VERSION`.
#[test]
fn test_version_names_samtools() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let version = ntools(dir.path(), &["--version"])?;
    let lines: Vec<_> = version.lines().collect();

    assert_eq!(
        lines,
        [
            concat!("ntools ", env!("CARGO_PKG_VERSION")),
            concat!("Using pasta ", env!("CARGO_PKG_VERSION")),
            "Compatible with samtools 1.24",
        ]
    );

    Ok(())
}

/// Numeric options take negative values, as with samtools' getopt and atoi: `view -q -5` keeps
/// every read. `sort -M` (minimiser sort) is reported as unsupported.
#[test]
fn test_negative_numbers_and_sort_minimiser() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    write_bam(dir.path(), 3)?;

    assert_eq!(
        ntools(dir.path(), &["view", "-c", "-q", "-5", "in.bam"])?,
        "3\n"
    );
    assert_eq!(
        ntools(dir.path(), &["view", "-c", "-q", "61", "in.bam"])?,
        "0\n"
    );
    assert_eq!(
        ntools(dir.path(), &["view", "-c", "-@", "-1", "in.bam"])?,
        "3\n"
    );

    let output = run(dir.path(), &["sort", "-M", "in.bam"])?;
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "ntools: sort -M (minimiser sort) is not supported\n"
    );

    Ok(())
}
