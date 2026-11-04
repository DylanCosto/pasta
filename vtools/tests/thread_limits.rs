//! When the system cannot start all the threads asked for (an address-space limit, `ulimit -v`,
//! or a process limit, `ulimit -u`), vtools runs with the threads it can start and the output of
//! a serial run, or fails as bcftools does with "Failed to create threads" (status 255) where
//! bcftools' own thread pool fails; it never panics. Tables sized by options (`stats -d`) that
//! cannot be allocated fail with a message and status 255 instead of aborting.
//!
//! Linux only: macOS rejects `ulimit -v`, and its `ulimit -u` does not count threads.
#![cfg(target_os = "linux")]

use std::{
    fs,
    io::Read as _,
    os::unix::process::ExitStatusExt,
    path::{Path, PathBuf},
    process::{Command, Output, Stdio},
};

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
    ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
    chr1\t10\t.\tG\tT\t50\tPASS\tDP=5\tGT\t1/1\t0/1\n\
    chr1\t20\t.\tA\tC,G\t20\tPASS\tDP=7\tGT\t0/1\t1/2\n\
    chr1\t30\t.\tAT\tA\t40\tPASS\tDP=3\tGT\t0/0\t0/1\n";

/// Records after those of [`VCF`], for `concat`.
const LATER_VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
    ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
    chr1\t40\t.\tC\tA\t30\tPASS\tDP=4\tGT\t0/1\t0/0\n\
    chr1\t50\t.\tG\tA\t60\tPASS\tDP=9\tGT\t1/1\t1/1\n";

/// The address-space limit these tests run under, in KB.
const ADDRESS_SPACE_KB: &str = "8000000";

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "vtools-thread-limits-{name}-{}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).expect("create the test directory");
        fs::write(dir.join("in.vcf"), VCF).expect("write the input");
        fs::write(dir.join("later.vcf"), LATER_VCF).expect("write the second input");
        for name in ["in", "later"] {
            let made = Command::new(env!("CARGO_BIN_EXE_vtools"))
                .args(["view", "-Oz", "-o"])
                .arg(format!("{name}.vcf.gz"))
                .arg("--write-index")
                .arg(format!("{name}.vcf"))
                .current_dir(&dir)
                .status()
                .expect("run vtools");
            assert!(made.success(), "bgzip and index {name}.vcf: {made}");
        }
        Self(dir)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Runs vtools with `args` in `dir` after the bash command `limit` (such as `ulimit -v N`).
fn run_limited(dir: &Path, limit: &str, args: &[&str]) -> Output {
    run_wrapped(dir, limit, "", args)
}

/// Runs vtools with `args` in `dir` after the bash command `limit`, started by `wrapper` (such as
/// `taskset -c 0`) when it is not empty.
fn run_wrapped(dir: &Path, limit: &str, wrapper: &str, args: &[&str]) -> Output {
    Command::new("bash")
        .arg("-c")
        .arg(format!("{limit} && exec {wrapper} \"$0\" \"$@\""))
        .arg(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .current_dir(dir)
        .stdin(Stdio::null())
        .env_remove("MALLOC_ARENA_MAX")
        .env_remove("GLIBC_TUNABLES")
        .output()
        .expect("run vtools through bash")
}

/// Whether a process can run in a user namespace of its own (`unshare -U`), where `ulimit -u`
/// counts only its own threads, so that the limit leaves room for exactly so many, whatever
/// else the user runs.
fn user_namespaces() -> bool {
    Command::new("unshare")
        .args(["-U", "true"])
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Runs vtools with `args` in `dir`, in a user namespace of its own where it may start `room`
/// threads besides its main one (see [`user_namespaces`]).
fn run_with_room(dir: &Path, room: usize, args: &[&str]) -> Output {
    Command::new("unshare")
        .args(["-U", "bash", "-c"])
        .arg(format!("ulimit -Su {} && exec \"$0\" \"$@\"", room + 1))
        .arg(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .current_dir(dir)
        .stdin(Stdio::null())
        .output()
        .expect("run vtools through unshare")
}

/// Whether the test cannot limit the threads of vtools: the process limit does not hold for
/// root.
fn is_root() -> bool {
    // SAFETY: geteuid has no preconditions.
    unsafe { libc::geteuid() == 0 }
}

/// Writes a reference `ref.fa` and two indexed BAM files, `in1.bam` and `in2.bam`, of 2,000
/// reads of 50 bp every 25 bp: more reads than a batch of a reader holds.
fn write_alignments(dir: &Path) {
    use pasta::{bam, sam, sam::alignment::io::Write as _};

    const LENGTH: usize = 60_000;
    let pattern = b"ACGTTGCAAGCTTCGA";
    let reference: Vec<u8> = pattern.iter().copied().cycle().take(LENGTH).collect();
    let mut fasta = b">c1\n".to_vec();
    for line in reference.chunks(60) {
        fasta.extend_from_slice(line);
        fasta.push(b'\n');
    }
    fs::write(dir.join("ref.fa"), fasta).expect("write the reference");
    let mut text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:{LENGTH}\n");
    for (i, start) in (0..50_000).step_by(25).enumerate() {
        let mut seq = reference[start..start + 50].to_vec();
        // A mismatch in every fourth read, so that some positions have an alternate allele.
        if i % 4 == 0 {
            seq[10] = if seq[10] == b'A' { b'C' } else { b'A' };
        }
        text += &format!(
            "r{i}\t0\tc1\t{}\t60\t50M\t*\t0\t0\t{}\t{}\n",
            start + 1,
            String::from_utf8_lossy(&seq),
            "I".repeat(50)
        );
    }
    for name in ["in1.bam", "in2.bam"] {
        let path = dir.join(name);
        let mut reader = sam::io::Reader::new(text.as_bytes());
        let header = reader.read_header().expect("parse the SAM header");
        let file = fs::File::create(&path).expect("create the BAM file");
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).expect("write the BAM header");
        for record in reader.record_bufs(&header) {
            let record = record.expect("parse a SAM record");
            writer
                .write_alignment_record(&header, &record)
                .expect("write a BAM record");
        }
        writer.try_finish().expect("finish the BAM file");
        let bam::Index::Bai(index) = bam::fs::index(&path).expect("index the BAM file") else {
            panic!("a short contig gets a BAI index");
        };
        bam::bai::fs::write(dir.join(format!("{name}.bai")), &index).expect("write the index");
    }
}

/// `mpileup` runs on the inputs of [`write_alignments`]: all of them, and regions.
const MPILEUP: &[&[&str]] = &[
    &[
        "mpileup",
        "--no-version",
        "-a",
        "AD,DP",
        "--threads",
        "{T}",
        "-f",
        "ref.fa",
        "in1.bam",
        "in2.bam",
    ],
    &[
        "mpileup",
        "--no-version",
        "--threads",
        "{T}",
        "-r",
        "c1:100-20000,c1:30000-",
        "-f",
        "ref.fa",
        "in1.bam",
        "in2.bam",
    ],
];

fn with_threads<'a>(args: &[&'a str], threads: &'a str) -> Vec<&'a str> {
    args.iter()
        .map(|&a| if a == "{T}" { threads } else { a })
        .collect()
}

#[track_caller]
fn assert_no_panic(out: &Output, what: &str) {
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.signal().is_none() && out.status.code() != Some(101),
        "{what}: {}; stderr: {err}",
        out.status
    );
    assert!(!err.contains("panicked"), "{what} panicked: {err}");
}

const COMMANDS: &[&[&str]] = &[
    &[
        "view",
        "--no-version",
        "-Oz",
        "--threads",
        "{T}",
        "in.vcf.gz",
    ],
    &[
        "filter",
        "--no-version",
        "-e",
        "QUAL<30",
        "--threads",
        "{T}",
        "in.vcf.gz",
    ],
    &[
        "norm",
        "--no-version",
        "-m-",
        "-Ob",
        "--threads",
        "{T}",
        "in.vcf.gz",
    ],
    &[
        "annotate",
        "--no-version",
        "-x",
        "INFO/DP",
        "--threads",
        "{T}",
        "in.vcf.gz",
    ],
    &["stats", "--threads", "{T}", "in.vcf.gz"],
    &[
        "merge",
        "--no-version",
        "--force-samples",
        "--threads",
        "{T}",
        "in.vcf.gz",
        "in.vcf.gz",
    ],
    &[
        "+fill-tags",
        "--no-version",
        "--threads",
        "{T}",
        "in.vcf.gz",
        "--",
        "-t",
        "AN,AC",
    ],
];

/// The output of `out` without the line of `stats` that holds its command line.
fn output_text(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| !l.starts_with("# The command"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn many_threads_under_an_address_space_limit_run_with_the_serial_output() {
    let dir = TestDir::new("address-space");
    let limit = format!("ulimit -v {ADDRESS_SPACE_KB}");
    // Also on one CPU: glibc sizes its arenas by the machine's CPUs, not those the process may
    // run on.
    let pinned = Path::new("/usr/bin/taskset").exists() || Path::new("/bin/taskset").exists();
    // And with GLIBC_TUNABLES set for something other than the arenas, as some systems set it.
    let tuned = "env GLIBC_TUNABLES=glibc.malloc.tcache_count=7";
    let wrappers: &[&str] = if pinned {
        &["", "taskset -c 0", tuned]
    } else {
        &["", tuned]
    };
    for (args, wrapper) in COMMANDS
        .iter()
        .flat_map(|args| wrappers.iter().map(move |w| (args, w)))
    {
        let what = format!("{} {wrapper}", args.join(" "));
        let serial = run_wrapped(dir.path(), &limit, wrapper, &with_threads(args, "0"));
        assert!(
            serial.status.success(),
            "{what} with 0 threads: {}",
            serial.status
        );
        for threads in ["64", "1024"] {
            let out = run_wrapped(dir.path(), &limit, wrapper, &with_threads(args, threads));
            let what = format!("{what} with {threads} threads under ulimit -v");
            assert_no_panic(&out, &what);
            assert!(
                out.status.success(),
                "{what}: {}; stderr: {}",
                out.status,
                String::from_utf8_lossy(&out.stderr)
            );
            assert_eq!(output_text(&out), output_text(&serial), "{what}");
        }
    }
}

#[test]
fn no_thread_to_start_fails_as_bcftools_or_runs_serially() {
    if is_root() {
        eprintln!("skipped: the process limit does not hold for root");
        return;
    }
    let dir = TestDir::new("processes");
    // With one process allowed, and more running, no thread can start.
    let limit = "ulimit -Su 1";
    // mpileup reads each input on a thread of its own, or on the main thread when that thread
    // cannot start; bcftools reads them all on one.
    write_alignments(dir.path());
    for args in MPILEUP {
        let serial = run_limited(dir.path(), "true", &with_threads(args, "0"));
        assert!(
            serial.status.success(),
            "{}: {}",
            args.join(" "),
            serial.status
        );
        for threads in ["0", "4"] {
            let what = format!(
                "{} with {threads} threads and none to start",
                args.join(" ")
            );
            let out = run_limited(dir.path(), limit, &with_threads(args, threads));
            assert_no_panic(&out, &what);
            assert!(
                out.status.success(),
                "{what}: {}; stderr: {}",
                out.status,
                String::from_utf8_lossy(&out.stderr)
            );
            assert_eq!(output_text(&out), output_text(&serial), "{what}");
        }
    }
    for args in COMMANDS {
        let what = format!("{} with no thread to start", args.join(" "));
        let out = run_limited(dir.path(), limit, &with_threads(args, "4"));
        assert_no_panic(&out, &what);
        if args[0] == "filter" || args[0] == "+fill-tags" {
            // bcftools threads only their output, and goes on without threads.
            let serial = run_limited(dir.path(), "true", &with_threads(args, "0"));
            assert!(out.status.success(), "{what}: {}", out.status);
            assert_eq!(output_text(&out), output_text(&serial), "{what}");
        } else {
            // bcftools fails to create its reader's thread pool.
            assert_eq!(out.status.code(), Some(255), "{what}");
            let err = String::from_utf8_lossy(&out.stderr);
            assert!(err.contains("Failed to create threads"), "{what}: {err}");
        }
    }
}

/// Commands, besides [`COMMANDS`] and [`MPILEUP`], that start their workers after other
/// threads (the reader's, the writer's), with the file they write (`None` for standard output).
/// `pileup.vcf.gz` is the first run of [`MPILEUP`].
const LATE_WORKERS: &[(&[&str], Option<&str>)] = &[
    (&["stats", "--threads", "{T}", "in.vcf.gz"], None),
    (
        &[
            "stats",
            "--threads",
            "{T}",
            "-r",
            "chr1:15-100",
            "in.vcf.gz",
        ],
        None,
    ),
    (
        &[
            "concat",
            "--no-version",
            "-Oz",
            "--threads",
            "{T}",
            "in.vcf.gz",
            "later.vcf.gz",
        ],
        None,
    ),
    (
        &[
            "convert",
            "--hapsample",
            "hs.gz",
            "--threads",
            "{T}",
            "in.vcf.gz",
        ],
        Some("hs.gz.hap.gz"),
    ),
    (
        &[
            "call",
            "--no-version",
            "-m",
            "--threads",
            "{T}",
            "pileup.vcf.gz",
        ],
        None,
    ),
];

/// What a run wrote: standard output (without the line of `stats` that holds its command line)
/// or the file `written`, decompressed.
fn written_text(dir: &Path, out: &Output, written: Option<&str>) -> String {
    let Some(name) = written else {
        return String::from_utf8_lossy(&decompressed(&out.stdout)).into_owned();
    };
    let data = fs::read(dir.join(name)).unwrap_or_default();
    String::from_utf8_lossy(&decompressed(&data)).into_owned()
}

/// `data` decompressed when it is gzip (BGZF), with the line of `stats` that holds its command
/// line left out.
fn decompressed(data: &[u8]) -> Vec<u8> {
    let mut text = Vec::new();
    if data.starts_with(&[0x1f, 0x8b]) {
        flate2::read::MultiGzDecoder::new(data)
            .read_to_end(&mut text)
            .expect("decompress the output");
    } else {
        text = data.to_vec();
    }
    text.split_inclusive(|&b| b == b'\n')
        .filter(|l| !l.starts_with(b"# The command"))
        .flatten()
        .copied()
        .collect()
}

#[test]
fn workers_that_cannot_start_leave_the_work_to_a_serial_run() {
    // The command's check that a thread starts passes, but the threads of the reader and the
    // writer take the room before the workers start: the run goes on without them, with the
    // output of --threads 0, where it used to fail with "Failed to create threads" after
    // writing part of its output.
    if is_root() || !user_namespaces() {
        eprintln!("skipped: no user namespace of its own to limit the threads of vtools in");
        return;
    }
    let dir = TestDir::new("room");
    write_alignments(dir.path());
    let mut pileup = with_threads(MPILEUP[0], "0");
    pileup.extend(["-Oz", "-o", "pileup.vcf.gz"]);
    let made = run_limited(dir.path(), "true", &pileup);
    assert!(made.status.success(), "mpileup for call: {}", made.status);
    let commands = COMMANDS
        .iter()
        .chain(MPILEUP)
        .map(|&args| (args, None))
        .chain(LATE_WORKERS.iter().copied());
    for (args, written) in commands {
        let serial = run_limited(dir.path(), "true", &with_threads(args, "0"));
        assert!(
            serial.status.success(),
            "{}: {}",
            args.join(" "),
            serial.status
        );
        let expected = written_text(dir.path(), &serial, written);
        // mpileup runs without a thread at all, as bcftools does; the other commands fail
        // then as bcftools does (see above).
        let least = usize::from(args[0] != "mpileup");
        for room in least..=6 {
            let what = format!("{} with room for {room} threads", args.join(" "));
            if let Some(name) = written {
                let _ = fs::remove_file(dir.path().join(name));
            }
            let out = run_with_room(dir.path(), room, &with_threads(args, "4"));
            assert_no_panic(&out, &what);
            assert!(
                out.status.success(),
                "{what}: {}; stderr: {}",
                out.status,
                String::from_utf8_lossy(&out.stderr)
            );
            assert_eq!(written_text(dir.path(), &out, written), expected, "{what}");
        }
    }
}

#[test]
fn call_runs_serially_when_its_workers_cannot_start() {
    // call checks that a thread starts, joins it, and starts its workers: the thread it joined
    // may still count against the process limit then, so with room for one thread the workers
    // cannot always start. Reading plain VCF, nothing else takes the room.
    if is_root() || !user_namespaces() {
        eprintln!("skipped: no user namespace of its own to limit the threads of vtools in");
        return;
    }
    let dir = TestDir::new("call-room");
    write_alignments(dir.path());
    let mut pileup = with_threads(MPILEUP[0], "0");
    pileup.extend(["-o", "pileup.vcf"]);
    let made = run_limited(dir.path(), "true", &pileup);
    assert!(made.status.success(), "mpileup for call: {}", made.status);
    let args = [
        "call",
        "--no-version",
        "-mv",
        "--threads",
        "{T}",
        "pileup.vcf",
    ];
    let serial = run_limited(dir.path(), "true", &with_threads(&args, "0"));
    assert!(serial.status.success(), "call: {}", serial.status);
    for run in 0..20 {
        let what = format!("call with room for 1 thread (run {run})");
        let out = run_with_room(dir.path(), 1, &with_threads(&args, "4"));
        assert_no_panic(&out, &what);
        assert!(
            out.status.success(),
            "{what}: {}; stderr: {}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        );
        assert_eq!(output_text(&out), output_text(&serial), "{what}");
    }
}

#[test]
fn an_option_sized_table_too_large_to_allocate_is_an_error() {
    let dir = TestDir::new("tables");
    let limit = format!("ulimit -v {ADDRESS_SPACE_KB}");
    for args in [
        &["stats", "-d", "0,2000000000,1", "in.vcf.gz"][..],
        &["stats", "-u", "DP:0:100:2000000000", "in.vcf.gz"],
    ] {
        let what = args.join(" ");
        let out = run_limited(dir.path(), &limit, args);
        assert_no_panic(&out, &what);
        assert_eq!(out.status.code(), Some(255), "{what}");
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(err.contains("could not allocate"), "{what}: {err}");
    }
}

#[test]
fn workers_without_room_for_their_tables_leave_the_work_to_fewer() {
    let dir = TestDir::new("worker-tables");
    let limit = format!("ulimit -v {ADDRESS_SPACE_KB}");
    // Each set of tables takes 3.2 GB of address space: the main one fits, not one per worker.
    let args = [
        "stats",
        "-d",
        "0,200000000,1",
        "--threads",
        "{T}",
        "in.vcf.gz",
    ];
    let serial = run_limited(dir.path(), &limit, &with_threads(&args, "0"));
    assert!(
        serial.status.success(),
        "stats -d with 0 threads: {}",
        serial.status
    );
    for threads in ["1", "4"] {
        let out = run_limited(dir.path(), &limit, &with_threads(&args, threads));
        let what = format!("stats -d with {threads} threads under ulimit -v");
        assert_no_panic(&out, &what);
        assert!(
            out.status.success(),
            "{what}: {}; stderr: {}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        );
        assert_eq!(output_text(&out), output_text(&serial), "{what}");
    }
}
