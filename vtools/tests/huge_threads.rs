//! A huge `--threads` count (99999999999, 1215752191 as bcftools' C `int`) is cut to
//! 1024 threads, so commands finish with the output of a serial run instead of starting threads
//! until the system runs out of them.

use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output, Stdio},
    thread,
    time::{Duration, Instant},
};

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
    ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
    chr1\t10\t.\tG\tT\t50\tPASS\tDP=5\tGT\t1/1\t0/1\n\
    chr1\t20\t.\tA\tC,G\t20\tPASS\tDP=7\tGT\t0/1\t1/2\n\
    chr1\t30\t.\tAT\tA\t40\tPASS\tDP=3\tGT\t0/0\t0/1\n";

const HUGE: &str = "99999999999";

/// How long one command may take; a few seconds at most once the count is capped.
const TIME_LIMIT: Duration = Duration::from_secs(30);

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("vtools-huge-threads-{name}-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("create the test directory");
        fs::write(dir.join("in.vcf"), VCF).expect("write the input");
        let made = Command::new(env!("CARGO_BIN_EXE_vtools"))
            .args(["view", "-Oz", "-o", "in.vcf.gz", "--write-index", "in.vcf"])
            .current_dir(&dir)
            .status()
            .expect("run vtools");
        assert!(made.success(), "bgzip and index the input: {made}");
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

/// Runs vtools with `args`, `{T}` replaced by `threads`, killing it after [`TIME_LIMIT`].
fn run(dir: &Path, args: &[&str], threads: &str) -> Output {
    let args: Vec<String> = args.iter().map(|a| a.replace("{T}", threads)).collect();
    let mut child = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(&args)
        .current_dir(dir)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run vtools");
    // Read the output on threads, so that a full pipe cannot stall the child.
    let mut stdout = child.stdout.take().expect("piped stdout");
    let mut stderr = child.stderr.take().expect("piped stderr");
    let out = thread::spawn(move || read_all(&mut stdout));
    let err = thread::spawn(move || read_all(&mut stderr));
    let start = Instant::now();
    let status = loop {
        if let Some(status) = child.try_wait().expect("wait for vtools") {
            break status;
        }
        if start.elapsed() > TIME_LIMIT {
            let _ = child.kill();
            let _ = child.wait();
            panic!("vtools {} ran for over {TIME_LIMIT:?}", args.join(" "));
        }
        thread::sleep(Duration::from_millis(20));
    };
    Output {
        status,
        stdout: out.join().expect("read stdout"),
        stderr: err.join().expect("read stderr"),
    }
}

fn read_all(r: &mut impl std::io::Read) -> Vec<u8> {
    let mut buf = Vec::new();
    r.read_to_end(&mut buf).expect("read the child's output");
    buf
}

#[test]
fn huge_thread_counts_finish_with_the_serial_output() {
    let dir = TestDir::new("commands");
    let cases: &[&[&str]] = &[
        &["view", "--no-version", "--threads", "{T}", "in.vcf.gz"],
        &[
            "view",
            "--no-version",
            "-Oz",
            "--threads",
            "{T}",
            "in.vcf.gz",
        ],
        &["view", "--no-version", "-Ob", "--threads", "{T}", "in.vcf"],
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
            "query",
            "-f",
            "%POS[ %GT]\\n",
            "--threads",
            "{T}",
            "in.vcf.gz",
        ],
        &[
            "convert",
            "--no-version",
            "-Oz",
            "--threads",
            "{T}",
            "in.vcf.gz",
        ],
        &[
            "norm",
            "--no-version",
            "-m-",
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
        &[
            "concat",
            "--no-version",
            "-a",
            "--threads",
            "{T}",
            "in.vcf.gz",
            "in.vcf.gz",
        ],
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
            "isec",
            "--no-version",
            "-n=2",
            "-w1",
            "--threads",
            "{T}",
            "in.vcf.gz",
            "in.vcf.gz",
        ],
        &["stats", "--threads", "{T}", "in.vcf.gz"],
        &["gtcheck", "--threads", "{T}", "in.vcf.gz"],
        &[
            "roh",
            "-G30",
            "--AF-dflt",
            "0.4",
            "--threads",
            "{T}",
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
        &[
            "+setGT",
            "--no-version",
            "--threads",
            "{T}",
            "in.vcf.gz",
            "--",
            "-t",
            "a",
            "-n",
            "0",
        ],
        &["+counts", "--threads", "{T}", "in.vcf.gz"],
    ];
    for args in cases {
        let what = args.join(" ");
        let serial = run(dir.path(), args, "0");
        assert!(
            serial.status.success(),
            "{what} with 0 threads: {}",
            serial.status
        );
        let huge = run(dir.path(), args, HUGE);
        assert!(
            huge.status.success(),
            "{what} with {HUGE} threads: {}; stderr: {}",
            huge.status,
            String::from_utf8_lossy(&huge.stderr)
        );
        // stats and gtcheck print their command line, and gtcheck a timing; compare the rest.
        let keep = |out: &[u8]| -> String {
            String::from_utf8_lossy(out)
                .lines()
                .filter(|l| !l.starts_with("# The command") && !l.starts_with("INFO\tTime"))
                .map(|l| l.replace(HUGE, "0"))
                .collect::<Vec<_>>()
                .join("\n")
        };
        assert_eq!(
            keep(&huge.stdout),
            keep(&serial.stdout),
            "{what}: the output with {HUGE} threads"
        );
    }
}

#[test]
fn huge_thread_counts_index_as_a_serial_run_does() {
    let dir = TestDir::new("index");
    for (threads, index) in [("0", "t0.csi"), (HUGE, "huge.csi")] {
        let out = run(
            dir.path(),
            &["index", "--threads", "{T}", "-o", index, "in.vcf.gz"],
            threads,
        );
        assert!(
            out.status.success(),
            "index --threads {threads}: {}",
            out.status
        );
    }
    assert_eq!(
        fs::read(dir.path().join("huge.csi")).expect("read the index"),
        fs::read(dir.path().join("t0.csi")).expect("read the index"),
        "the index built with {HUGE} threads"
    );
}
