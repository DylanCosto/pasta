//! `mpileup` on SAM text that htslib's `sam_parse1` normalizes as it reads it.

use std::{fs, path::PathBuf, process::Command};

const REFERENCE: &str = ">c1\nACGTACGTAC\n";
const INDEX: &str = "c1\t10\t4\t10\t11\n";

/// `x` is on a reference sequence missing from the header, which bcftools reads as unmapped
/// with a warning, and `z` has no CIGAR, which makes it unmapped too.
const SAM: &str = "@SQ\tSN:c1\tLN:10
x\t0\tchrX\t2\t30\t4M\t*\t0\t0\tACGT\tIIII
y\t0\tc1\t2\t30\t4M\t*\t0\t0\tCGTA\tIIII
z\t0\tc1\t3\t30\t*\t*\t0\t0\tGTAC\tIIII
";

/// bcftools 1.24 `mpileup --no-version -f ref.fa in.sam`, without the `##` lines.
const EXPECTED: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tin.sam
c1\t2\t.\tC\t<*>\t0\t.\tDP=1;I16=1,0,0,0,40,1600,0,0,30,900,0,0,0,0,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,3,30:1,0
c1\t3\t.\tG\t<*>\t0\t.\tDP=1;I16=1,0,0,0,40,1600,0,0,30,900,0,0,1,1,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,3,30:1,0
c1\t4\t.\tT\t<*>\t0\t.\tDP=1;I16=1,0,0,0,40,1600,0,0,30,900,0,0,1,1,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,3,30:1,0
c1\t5\t.\tA\t<*>\t0\t.\tDP=1;I16=1,0,0,0,40,1600,0,0,30,900,0,0,0,0,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,3,30:1,0
";

struct TestDir(PathBuf);

impl Drop for TestDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn test_mpileup_reads_unknown_reference_sequence_names_as_unmapped() {
    let dir =
        TestDir(std::env::temp_dir().join(format!("vtools-mpileup-sam-{}", std::process::id())));
    fs::create_dir_all(&dir.0).expect("create the test directory");
    fs::write(dir.0.join("ref.fa"), REFERENCE).expect("write the reference");
    fs::write(dir.0.join("ref.fa.fai"), INDEX).expect("write the reference index");
    fs::write(dir.0.join("in.sam"), SAM).expect("write the input");

    let output = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["mpileup", "--no-version", "-f", "ref.fa", "in.sam"])
        .current_dir(&dir.0)
        .output()
        .expect("run vtools");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{}: {stderr}", output.status);

    let stdout = String::from_utf8(output.stdout).expect("UTF-8");
    let records: String = stdout
        .lines()
        .filter(|line| !line.starts_with("##"))
        .map(|line| format!("{line}\n"))
        .collect();

    assert_eq!(records, EXPECTED);
    assert!(
        stderr
            .contains("[W::sam_parse1] unrecognized reference name \"chrX\"; treated as unmapped"),
        "{stderr}"
    );
}
