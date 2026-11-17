//! `mpileup` on SAM text that htslib's `sam_parse1` normalizes as it reads it, and on reads
//! whose CIGAR puts them in the pileup past their last base.

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

/// Reads with a zero-length match operation where the query is used up, so the pileup has
/// an entry at the query position just past the end of the sequence: the base there is
/// unknown (`N`) and its quality 0, as samtools mpileup prints it.
const ZERO_LENGTH_MATCHES: [&str; 6] = [
    "a\t0\tc1\t2\t60\t2M0M2D\t*\t0\t0\tAC\tII\n",
    "a\t0\tc1\t2\t60\t0M2D\t*\t0\t0\t*\t*\n",
    "a\t0\tc1\t2\t60\t1S0M3D\t*\t0\t0\tA\tI\n",
    "a\t0\tc1\t2\t60\t2M1I0M3D\t*\t0\t0\tACG\tIII\n",
    "a\t0\tc1\t2\t60\t0=3N\t*\t0\t0\t*\t*\n",
    "a\t0\tc1\t2\t60\t3D0X3D\t*\t0\t0\t*\t*\n",
];

/// `2M0M2D` and a reference read, with `-B -Q 0 -a AD`: at position 4 the first read is an
/// `N` of quality 0 (counted at the lowest quality, 4), so it adds nothing to `AD` or `QS`.
const EXPECTED_ZERO_LENGTH_MATCH: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tin.sam
c1\t2\t.\tC\tA,<*>\t0\t.\tDP=2;I16=0,1,1,0,40,1600,40,1600,50,2500,60,3600,0,0,0,0;QS=0.5,0.5,0;SGB=-0.379885;RPBZ=1;MQBZ=1;MQSBZ=-1;BQBZ=0;SCBZ=0;MQ0F=0\tPL:AD\t34,0,34,37,37,71:1,1,0
c1\t3\t.\tG\tC,<*>\t0\t.\tDP=2;I16=0,1,1,0,40,1600,40,1600,50,2500,60,3600,1,1,0,0;QS=0.5,0.5,0;SGB=-0.379885;RPBZ=1;MQBZ=1;MQSBZ=-1;BQBZ=0;SCBZ=0;MQ0F=0\tPL:AD\t34,0,34,37,37,71:1,1,0
c1\t4\t.\tT\t<*>\t0\t.\tDP=2;I16=0,1,1,0,40,1600,0,0,50,2500,60,3600,1,1,-1,1;QS=1,0;MQ0F=0\tPL:AD\t0,3,37:1,0
c1\t5\t.\tA\t<*>\t0\t.\tDP=1;I16=0,1,0,0,40,1600,0,0,50,2500,0,0,0,0,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,3,40:1,0
";

fn mpileup(dir: &TestDir, sam: &str, args: &[&str]) -> (std::process::ExitStatus, String, String) {
    fs::write(dir.0.join("in.sam"), format!("@SQ\tSN:c1\tLN:10\n{sam}")).expect("write the input");
    let output = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["mpileup", "--no-version", "-f", "ref.fa"])
        .args(args)
        .arg("in.sam")
        .current_dir(&dir.0)
        .output()
        .expect("run vtools");
    let stdout = String::from_utf8(output.stdout).expect("UTF-8");
    let records = stdout
        .lines()
        .filter(|line| !line.starts_with("##"))
        .map(|line| format!("{line}\n"))
        .collect();
    (
        output.status,
        records,
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

#[test]
fn test_mpileup_reads_a_zero_length_match_past_the_sequence_as_an_unknown_base() {
    let dir = TestDir(
        std::env::temp_dir().join(format!("vtools-mpileup-zero-length-{}", std::process::id())),
    );
    fs::create_dir_all(&dir.0).expect("create the test directory");
    fs::write(dir.0.join("ref.fa"), REFERENCE).expect("write the reference");
    fs::write(dir.0.join("ref.fa.fai"), INDEX).expect("write the reference index");

    for sam in ZERO_LENGTH_MATCHES {
        for args in [
            &[][..],
            &["-B", "-Q", "0"],
            &["-C", "50", "-A"],
            &["--indels-cns"],
        ] {
            let (status, _, stderr) = mpileup(&dir, sam, args);
            assert!(status.success(), "{sam:?} {args:?}: {status}: {stderr}");
        }
    }

    let (status, records, stderr) = mpileup(
        &dir,
        &format!(
            "{}r\t16\tc1\t2\t50\t4M\t*\t0\t0\tCGTA\tIIII\n",
            ZERO_LENGTH_MATCHES[0]
        ),
        &["-B", "-Q", "0", "-a", "AD"],
    );
    assert!(status.success(), "{status}: {stderr}");
    assert_eq!(records, EXPECTED_ZERO_LENGTH_MATCH);
}
