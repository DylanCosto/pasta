//! Commands that read several inputs in step reject inputs whose reference sequences differ
//! from the first's, rather than putting reads on the wrong sequence.

use std::{
    fs,
    path::Path,
    process::{Command, Output},
};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

/// Covers chr1.
const A: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
@SQ\tSN:chr2\tLN:100
a1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII
";

/// The same sequences in another order, covering chr2.
const B: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr2\tLN:100
@SQ\tSN:chr1\tLN:100
b1\t0\tchr2\t20\t60\t4M\t*\t0\t0\tACGT\tIIII
";

fn ntools(dir: &Path, args: &[&str]) -> Output {
    Command::new(NTOOLS)
        .args(args)
        .current_dir(dir)
        .output()
        .expect("run ntools")
}

/// Writes `A.bam` and `B.bam` to a new directory.
fn inputs() -> tempfile::TempDir {
    let dir = tempfile::tempdir().expect("create the test directory");

    for (name, sam) in [("A", A), ("B", B)] {
        fs::write(dir.path().join(format!("{name}.sam")), sam).expect("write the input");
        let output = ntools(
            dir.path(),
            &[
                "view",
                "-b",
                "-o",
                &format!("{name}.bam"),
                &format!("{name}.sam"),
            ],
        );
        assert!(output.status.success(), "{output:?}");
    }

    dir
}

fn assert_rejected(command: &str, args: &[&str]) {
    let dir = inputs();
    let mut all = vec![command];
    all.extend(args);
    all.extend(["A.bam", "B.bam"]);

    let output = ntools(dir.path(), &all);

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{command}: {stderr}");
    assert!(
        stderr.contains(&format!(
            "[{command}] the reference sequences of B.bam differ from those of A.bam"
        )),
        "{command}: {stderr}"
    );
    assert!(output.stdout.is_empty(), "{command}: {:?}", output.stdout);
}

#[test]
fn test_depth_rejects_differing_reference_sequences() {
    assert_rejected("depth", &[]);
}

#[test]
fn test_coverage_rejects_differing_reference_sequences() {
    assert_rejected("coverage", &[]);
}

#[test]
fn test_cat_rejects_differing_reference_sequences() {
    assert_rejected("cat", &[]);
}

#[test]
fn test_depth_accepts_matching_reference_sequences() {
    let dir = inputs();

    let output = ntools(dir.path(), &["depth", "A.bam", "A.bam"]);

    assert!(output.status.success(), "{output:?}");
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        "chr1\t10\t1\t1\nchr1\t11\t1\t1\nchr1\t12\t1\t1\nchr1\t13\t1\t1\n"
    );
}

/// Only chr1, A's first sequence, so its reads keep their sequence when read after A.
const C: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:100
c1\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII
";

#[test]
fn test_depth_accepts_an_input_with_the_first_inputs_first_sequences() {
    let dir = inputs();
    fs::write(dir.path().join("C.sam"), C).expect("write the input");
    let output = ntools(dir.path(), &["view", "-b", "-o", "C.bam", "C.sam"]);
    assert!(output.status.success(), "{output:?}");

    let output = ntools(dir.path(), &["depth", "A.bam", "C.bam"]);
    assert!(output.status.success(), "{output:?}");
    let depth = String::from_utf8_lossy(&output.stdout);
    assert!(depth.lines().any(|l| l == "chr1\t30\t0\t1"), "{depth}");

    // The other way round, C's header lacks A's chr2.
    let output = ntools(dir.path(), &["depth", "C.bam", "A.bam"]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
}
