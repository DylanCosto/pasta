//! With standard output closed (`>&-`), a command that writes to it fails with a write error,
//! as bcftools does, instead of exiting with status 0 having written nothing, and no output goes
//! into a file vtools opened for something else.

use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output, Stdio},
};

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
    chr1\t10\t.\tG\tT\t50\tPASS\t.\tGT\t1/1\n\
    chr1\t20\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n";

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "vtools-closed-stdout-{name}-{}",
            std::process::id()
        ));
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

/// Runs vtools with `args` in `dir` with standard output closed, and with standard input closed
/// too when `close_stdin` is set.
fn run_closed(dir: &Path, args: &[&str], close_stdin: bool) -> Output {
    let script = if close_stdin {
        "exec \"$0\" \"$@\" <&- >&-"
    } else {
        "exec \"$0\" \"$@\" >&-"
    };
    Command::new("sh")
        .arg("-c")
        .arg(script)
        .arg(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .current_dir(dir)
        .stdin(Stdio::null())
        .output()
        .expect("run vtools through sh")
}

#[track_caller]
fn assert_write_error(out: &Output, what: &str) {
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(255),
        "{what}: {} (bcftools exits with 255); stderr: {err}",
        out.status
    );
    assert!(!err.is_empty(), "{what} reports the write error");
}

#[test]
fn closed_stdout_fails_the_commands_that_write_to_it() {
    let dir = TestDir::new("commands");
    let cases: &[&[&str]] = &[
        &["view", "in.vcf"],
        &["view", "-Ob", "in.vcf"],
        &["view", "--threads", "2", "-Oz", "in.vcf"],
        &["view", "-o", "-", "in.vcf"],
        &["query", "-f", "%POS\\n", "in.vcf"],
        &["query", "-l", "in.vcf"],
        &["norm", "-m-", "in.vcf"],
        &["sort", "-T", "./tmp.", "in.vcf"],
        &["head", "in.vcf"],
        &["stats", "in.vcf"],
        &["index", "-s", "in.vcf.gz"],
        &["index", "-n", "in.vcf.gz"],
        &["+counts", "in.vcf"],
        &["+fill-tags", "in.vcf"],
    ];
    for args in cases {
        let out = run_closed(dir.path(), args, false);
        assert_write_error(&out, &args.join(" "));
    }
    // The output fits in the file buffer, so the write fails when the output is closed.
    let out = run_closed(dir.path(), &["view", "in.vcf"], false);
    assert_eq!(
        String::from_utf8_lossy(&out.stderr),
        "Error: close failed stdout\n"
    );
}

#[test]
fn closed_stdin_and_stdout_fail_the_write() {
    let dir = TestDir::new("stdin");
    let out = run_closed(dir.path(), &["view", "in.vcf"], true);
    assert_write_error(&out, "view with standard input closed too");
}

#[test]
fn output_never_goes_into_a_file_opened_for_something_else() {
    let dir = TestDir::new("reuse");
    let input = fs::read(dir.path().join("in.vcf.gz")).expect("read the input");
    // bcftools opens the input at the free descriptor 1, so /dev/stdout is the input, which it
    // truncates and overwrites with its output.
    for target in ["/dev/stdout", "/dev/fd/1"] {
        for ty in ["-Ov", "-Oz"] {
            let args = ["view", ty, "-o", target, "in.vcf.gz"];
            let out = run_closed(dir.path(), &args, false);
            assert_write_error(&out, &args.join(" "));
            assert_eq!(
                fs::read(dir.path().join("in.vcf.gz")).expect("read the input"),
                input,
                "{args:?} leaves its input as it was"
            );
        }
    }
}

#[test]
fn closed_stdout_does_not_fail_other_outputs() {
    let dir = TestDir::new("other");
    let out = run_closed(dir.path(), &["view", "-o", "out.vcf", "in.vcf.gz"], false);
    assert!(out.status.success(), "{}", out.status);
    let written = fs::read_to_string(dir.path().join("out.vcf")).expect("read the output");
    assert!(written.ends_with("chr1\t20\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n"));
    let out = run_closed(dir.path(), &["index", "-f", "in.vcf.gz"], false);
    assert!(out.status.success(), "{}", out.status);
}
