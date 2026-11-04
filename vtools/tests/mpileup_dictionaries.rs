//! `mpileup` rejects files whose reference sequences differ from the first file's, rather
//! than piling their reads up on the wrong sequence.

use std::{fs, process::Command};

const REFERENCE: &str = ">chr1\nACGTACGTAC\n>chr2\nACGTACGTAC\n";
const INDEX: &str = "chr1\t10\t6\t10\t11\nchr2\t10\t23\t10\t11\n";

const A: &str = "@SQ\tSN:chr1\tLN:10
@SQ\tSN:chr2\tLN:10
a1\t0\tchr1\t2\t30\t4M\t*\t0\t0\tCGTA\tIIII
";

/// The same sequences in another order, covering chr2.
const B: &str = "@SQ\tSN:chr2\tLN:10
@SQ\tSN:chr1\tLN:10
b1\t0\tchr2\t2\t30\t4M\t*\t0\t0\tCTTA\tIIII
";

#[test]
fn test_mpileup_rejects_differing_reference_sequences() {
    let dir = std::env::temp_dir().join(format!("vtools-mpileup-dict-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("create the test directory");
    fs::write(dir.join("ref.fa"), REFERENCE).expect("write the reference");
    fs::write(dir.join("ref.fa.fai"), INDEX).expect("write the reference index");
    fs::write(dir.join("A.sam"), A).expect("write A");
    fs::write(dir.join("B.sam"), B).expect("write B");

    let output = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["mpileup", "--no-version", "-f", "ref.fa", "A.sam", "B.sam"])
        .current_dir(&dir)
        .output()
        .expect("run vtools");
    let _ = fs::remove_dir_all(&dir);

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains("[mpileup] the reference sequences of B.sam differ from those of A.sam"),
        "{stderr}"
    );
    assert!(output.stdout.is_empty(), "{:?}", output.stdout);
}
