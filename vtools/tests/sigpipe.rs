//! A write to a pipe whose reader went away (`| head`) kills vtools with SIGPIPE and no message,
//! as it kills bcftools.

use std::{
    fs,
    io::{self, Read},
    os::unix::process::ExitStatusExt,
    path::{Path, PathBuf},
    process::{Command, ExitStatus, Stdio},
};

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
    chr1\t20\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n\
    chr1\t10\t.\tG\tT\t50\tPASS\t.\tGT\t1/1\n";

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("vtools-sigpipe-{name}-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("create the test directory");
        fs::write(dir.join("in.vcf"), VCF).expect("write the input");
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

/// Runs vtools with `args` in `dir`, with standard output (or, with `Closed::Stderr`, standard
/// error) a pipe whose reading end is closed. Returns the exit status and what it wrote to the
/// other stream.
fn run_into_closed_pipe(dir: &Path, args: &[&str], closed: Closed) -> (ExitStatus, String) {
    let (reader, writer) = io::pipe().expect("make a pipe");
    drop(reader);
    let (mut other, other_writer) = io::pipe().expect("make a pipe");
    let mut command = Command::new(env!("CARGO_BIN_EXE_vtools"));
    command.args(args).current_dir(dir).stdin(Stdio::null());
    match closed {
        Closed::Stdout => command.stdout(writer).stderr(other_writer),
        Closed::Stderr => command.stderr(writer).stdout(other_writer),
    };
    let mut child = command.spawn().expect("run vtools");
    // Close this process's copies of the write ends, so that reading stops when vtools ends.
    drop(command);
    let mut text = String::new();
    other
        .read_to_string(&mut text)
        .expect("read the other stream");
    (child.wait().expect("wait for vtools"), text)
}

#[derive(Clone, Copy)]
enum Closed {
    Stdout,
    Stderr,
}

#[track_caller]
fn assert_sigpipe(status: ExitStatus, what: &str) {
    assert_eq!(
        status.signal(),
        Some(libc::SIGPIPE),
        "{what}: {status} (bcftools is killed by SIGPIPE)"
    );
}

#[test]
fn closed_output_pipe_kills_with_sigpipe() {
    let dir = TestDir::new("stdout");
    let cases: &[&[&str]] = &[
        &["view", "in.vcf"],
        &["view", "--threads", "2", "-Oz", "in.vcf"],
        &["view", "-Ob", "in.vcf"],
        &["query", "-f", "%POS\\n", "in.vcf"],
        &["norm", "-m-", "in.vcf"],
        &["stats", "in.vcf"],
        &["+fill-tags", "in.vcf"],
    ];
    for args in cases {
        let (status, err) = run_into_closed_pipe(dir.path(), args, Closed::Stdout);
        assert_sigpipe(status, &args.join(" "));
        assert_eq!(err, "", "{args:?} writes nothing to standard error");
    }
}

#[test]
fn sort_removes_its_temporary_directory_before_sigpipe_kills_it() {
    let dir = TestDir::new("sort");
    let prefix = format!("{}/tmp.", dir.path().display());
    for max_mem in ["768M", "100"] {
        let args = ["sort", "-m", max_mem, "-T", &prefix, "in.vcf"];
        let (status, err) = run_into_closed_pipe(dir.path(), &args, Closed::Stdout);
        assert_sigpipe(status, &format!("sort -m {max_mem}"));
        assert!(err.starts_with(&format!("Writing to {prefix}")), "{err}");
        // bcftools prints "Done" before it closes the output, whose flush is the write that
        // finds the pipe closed, and no clean-up message.
        assert!(err.ends_with("\nDone\n"), "Done before the close: {err}");
        let left: Vec<_> = fs::read_dir(dir.path())
            .expect("list the test directory")
            .map(|e| e.expect("a directory entry").file_name())
            .filter(|name| name != "in.vcf")
            .collect();
        assert!(left.is_empty(), "sort -m {max_mem} left {left:?}");
    }
}

#[test]
fn closed_stderr_pipe_kills_with_sigpipe() {
    let dir = TestDir::new("stderr");
    let (status, out) = run_into_closed_pipe(dir.path(), &["view", "nosuch.vcf"], Closed::Stderr);
    assert_sigpipe(status, "view nosuch.vcf");
    assert_eq!(out, "");
}

#[test]
fn roh_prints_its_summary_before_sigpipe_kills_it() {
    let dir = TestDir::new("roh");
    let sorted = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
        chr1\t10\t.\tG\tT\t50\tPASS\t.\tGT\t1/1\n\
        chr1\t20\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n";
    fs::write(dir.path().join("sorted.vcf"), sorted).expect("write the input");
    let (status, err) =
        run_into_closed_pipe(dir.path(), &["roh", "-G30", "sorted.vcf"], Closed::Stdout);
    assert_sigpipe(status, "roh -G30");
    // bcftools prints the summary before it closes (and flushes) its output.
    assert_eq!(
        err,
        "Number of target samples: 1\n\
         Number of --estimate-AF samples: 0\n\
         Number of sites in the buffer/overlap: unlimited\n\
         Number of lines total/processed: 2/0\n\
         Number of lines filtered/no AF/no alt/multiallelic/dup: 0/2/0/0/0\n\
         No usable sites were found.\n \
         Consider using one of the AF options.\n"
    );
}
