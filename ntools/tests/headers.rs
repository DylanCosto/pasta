//! Headers that break the SAM specification but that htslib reads.

use std::{path::PathBuf, process::Command};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// `hdr_duppg.bam` has two `@PG` records with `ID:bwa`, as merges and re-runs leave, and 62
/// records. samtools 1.24 reads it silently.
#[test]
fn test_view_count_with_duplicate_program_ids() {
    let output = Command::new(NTOOLS)
        .args(["view", "-c"])
        .arg(data("hdr_duppg.bam"))
        .output()
        .expect("ntools runs");

    assert!(
        output.status.success(),
        "ntools view -c: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(String::from_utf8_lossy(&output.stdout), "62\n");
    assert_eq!(String::from_utf8_lossy(&output.stderr), "");
}

/// Both `@PG` records are kept as written, as samtools keeps them, and a new `@PG` record links to
/// the last, which htslib resolves `bwa` to.
#[test]
fn test_view_header_with_duplicate_program_ids() {
    const TEXT: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:2000
@SQ\tSN:chr2\tLN:1500
@RG\tID:rg1\tSM:s1
@PG\tID:bwa\tPN:bwa
@PG\tID:bwa\tPN:bwa\tVN:2
";

    let view_header = |args: &[&str]| {
        let output = Command::new(NTOOLS)
            .args(["view", "-H"])
            .args(args)
            .arg(data("hdr_duppg.bam"))
            .output()
            .expect("ntools runs");

        assert!(
            output.status.success(),
            "ntools view -H: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        String::from_utf8(output.stdout).expect("UTF-8")
    };

    assert_eq!(view_header(&["--no-PG"]), TEXT);

    let header = view_header(&[]);
    let (text, program) = header.split_at(TEXT.len());
    assert_eq!(text, TEXT);
    assert!(
        program.starts_with("@PG\tID:ntools\tPN:ntools\tPP:bwa\tVN:"),
        "{program}"
    );
    assert_eq!(program.lines().count(), 1);
}

/// An input that is not SAM, BAM or CRAM is rejected with its name and what it looks like, with
/// status 1 as samtools fails to read its header, while a SAM file with a bare `@CO` line or no
/// header at all is read.
#[test]
fn test_other_formats_are_named_in_errors() -> std::io::Result<()> {
    use std::io::Write;

    let dir = tempfile::tempdir()?;
    let dir = dir.path();
    let record = "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(b"##fileformat=VCFv4.2\n#CHROM\tPOS\n")?;
    let inputs: [(&str, Vec<u8>, &str); 5] = [
        (
            "t.vcf",
            b"##fileformat=VCFv4.2\n".to_vec(),
            "it looks like VCF",
        ),
        ("t.vcf.gz", gz.finish()?, "it looks like VCF"),
        ("t.fa", b">chr1\nACGT\n".to_vec(), "it looks like FASTA"),
        (
            "t.fq",
            b"@r1\nACGT\n+\nIIII\n".to_vec(),
            "it looks like FASTQ",
        ),
        (
            "t.txt",
            b"hello world\n".to_vec(),
            "its first line is neither",
        ),
    ];

    for (name, contents, kind) in &inputs {
        std::fs::write(dir.join(name), contents)?;

        for command in ["view", "flagstat", "sort", "fastq", "depth"] {
            let output = Command::new(NTOOLS)
                .current_dir(dir)
                .args([command, name])
                .output()?;
            let stderr = String::from_utf8_lossy(&output.stderr);

            assert_eq!(output.status.code(), Some(1), "{command} {name}: {stderr}");
            assert!(
                stderr.starts_with(&format!("ntools: {name}: not SAM, BAM or CRAM ({kind}")),
                "{command} {name}: {stderr}"
            );
        }
    }

    std::fs::write(dir.join("bad.sam"), "@HD\tVN:1.6\nr1\t4\t*\n")?;
    let output = Command::new(NTOOLS)
        .current_dir(dir)
        .args(["view", "bad.sam"])
        .output()?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.starts_with("ntools: bad.sam: reading SAM record 1: unexpected end of line"),
        "{stderr}"
    );

    std::fs::write(dir.join("co.sam"), format!("@HD\tVN:1.6\n@CO\n{record}"))?;
    std::fs::write(dir.join("bare.sam"), record)?;

    for name in ["co.sam", "bare.sam"] {
        let output = Command::new(NTOOLS)
            .current_dir(dir)
            .args(["view", name])
            .output()?;
        assert!(
            output.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(String::from_utf8_lossy(&output.stdout), record, "{name}");
    }

    Ok(())
}
