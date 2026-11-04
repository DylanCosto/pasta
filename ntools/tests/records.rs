//! SAM records that htslib's `sam_parse1` normalizes as it reads them.

use std::{
    io::Write as _,
    process::{Command, Output, Stdio},
};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

/// The example: of six records, samtools 1.24 counts two (`e`, `f`) as mapped.
const UNMAPPED: &str = "@SQ\tSN:c1\tLN:1000
a\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
b\t65\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
c\t0\tc1\t10\t30\t*\t*\t0\t0\tACGT\tIIII
d\t0\tc1\t0\t0\t4M\t*\t0\t0\tACGT\tIIII
e\t0\tc1\t20\t30\t4M\t*\t0\t0\tACGT\tIIII
f\t1\tc1\t30\t30\t4M\tc1\t0\t0\tACGT\tIIII
";

/// A record on a reference sequence missing from the header, e.g., a decoy filtered out of it.
/// samtools keeps its POS (`* 10`), placing it on no reference sequence; ntools drops the
/// position with the reference sequence (README).
const UNKNOWN_REFERENCE: &str = "@SQ\tSN:c1\tLN:1000
x\t0\tchrX\t10\t30\t4M\t*\t0\t0\tACGT\tIIII
y\t0\tc1\t10\t30\t4M\t*\t0\t0\tACGT\tIIII
";

fn ntools(args: &[&str], stdin: &[u8]) -> Output {
    let mut child = Command::new(NTOOLS)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("ntools runs");

    child
        .stdin
        .take()
        .expect("stdin is piped")
        .write_all(stdin)
        .expect("stdin is written");

    let output = child.wait_with_output().expect("ntools runs");

    assert!(
        output.status.success(),
        "ntools {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    output
}

fn stdout(output: &Output) -> &str {
    std::str::from_utf8(&output.stdout).expect("UTF-8")
}

#[test]
fn test_view_normalizes_placement_as_samtools() {
    let output = ntools(&["view", "-"], UNMAPPED.as_bytes());

    assert_eq!(
        stdout(&output),
        "a\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
b\t69\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
c\t4\tc1\t10\t30\t*\t*\t0\t0\tACGT\tIIII
d\t4\t*\t0\t0\t4M\t*\t0\t0\tACGT\tIIII
e\t0\tc1\t20\t30\t4M\t*\t0\t0\tACGT\tIIII
f\t1\tc1\t30\t30\t4M\t*\t0\t0\tACGT\tIIII
"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("[W::sam_parse1] mapped query must have a CIGAR; treated as unmapped"),
        "{stderr}"
    );
}

#[test]
fn test_flagstat_counts_normalized_records_as_samtools() {
    let output = ntools(&["flagstat", "-"], UNMAPPED.as_bytes());

    assert_eq!(
        stdout(&output),
        "6 + 0 in total (QC-passed reads + QC-failed reads)
6 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2 + 0 mapped (33.33% : N/A)
2 + 0 primary mapped (33.33% : N/A)
2 + 0 paired in sequencing
1 + 0 read1
0 + 0 read2
0 + 0 properly paired (0.00% : N/A)
1 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
1 + 0 with mate mapped to a different chr
1 + 0 with mate mapped to a different chr (mapQ>=5)
"
    );

    // BAM written from the SAM keeps the normalized flags.
    let bam = ntools(&["view", "-b", "-"], UNMAPPED.as_bytes());
    let output = ntools(&["view", "-c", "-F", "4", "-"], &bam.stdout);
    assert_eq!(stdout(&output), "2\n");
}

#[test]
fn test_unknown_reference_sequence_name_is_unmapped() {
    let output = ntools(&["view", "-"], UNKNOWN_REFERENCE.as_bytes());

    assert_eq!(
        stdout(&output),
        "x\t4\t*\t0\t30\t4M\t*\t0\t0\tACGT\tIIII
y\t0\tc1\t10\t30\t4M\t*\t0\t0\tACGT\tIIII
"
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "[W::sam_parse1] unrecognized reference name \"chrX\"; treated as unmapped\n"
    );

    let output = ntools(
        &["sort", "-O", "sam", "--no-PG", "-"],
        UNKNOWN_REFERENCE.as_bytes(),
    );

    assert_eq!(
        stdout(&output),
        "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:c1\tLN:1000
y\t0\tc1\t10\t30\t4M\t*\t0\t0\tACGT\tIIII
x\t4\t*\t0\t30\t4M\t*\t0\t0\tACGT\tIIII
"
    );
}

/// `fastq --no-sc` remaps base modifications to the remaining bases, in the original
/// orientation, dropping calls on clipped bases with their `ML` values and updating `MN`.
/// (samtools 1.24 keeps the tags of the untrimmed sequence.)
#[test]
fn test_fastq_no_sc_remaps_base_modifications() {
    let sam = "@SQ\tSN:c1\tLN:1000
f\t0\tc1\t100\t60\t2S4M\t*\t0\t0\tCACGCC\tABCDEF\tMM:Z:C+m,0,0;\tML:B:C,10,200\tMN:i:6
r\t16\tc1\t100\t60\t4M2S\t*\t0\t0\tGGCGTG\tFEDCBA\tMM:Z:C+m,0,0;\tML:B:C,10,200
";

    let output = ntools(&["fastq", "--no-sc", "-T", "MM,ML,MN", "-"], sam.as_bytes());
    assert!(output.status.success());

    let headers: Vec<String> = String::from_utf8_lossy(&output.stdout)
        .lines()
        .step_by(4)
        .map(String::from)
        .collect();

    assert_eq!(
        headers,
        [
            "@f\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:4",
            "@r\tMM:Z:C+m,0;\tML:B:C,200",
        ]
    );
}
