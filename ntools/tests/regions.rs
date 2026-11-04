//! Regions and BED intervals that end past the largest position an index can hold.

use std::{fs, io, path::Path, process::Command};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

/// Runs ntools in `dir`, failing the test if it fails, and returns its standard output.
fn ntools(dir: &Path, args: &[&str]) -> io::Result<String> {
    let output = Command::new(NTOOLS).current_dir(dir).args(args).output()?;
    assert!(
        output.status.success(),
        "ntools {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).map_err(io::Error::other)
}

/// Writes a reference, and indexed BAM (BAI and CSI) and CRAM copies of records on two reference
/// sequences, to `dir`.
fn write_inputs(dir: &Path) -> io::Result<()> {
    let sequence: String = (0..2000).map(|i| ['A', 'C', 'G', 'T'][i * 7 % 4]).collect();
    fs::write(
        dir.join("ref.fa"),
        format!(">chr1\n{sequence}\n>chr2\n{sequence}\n"),
    )?;

    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    sam.push_str("@SQ\tSN:chr1\tLN:2000\n@SQ\tSN:chr2\tLN:2000\n");

    for (name, count) in [("chr1", 41), ("chr2", 17)] {
        for i in 1..=count {
            let pos = i * 40;
            let bases = &sequence[pos - 1..pos + 19];
            sam.push_str(&format!(
                "{name}r{i}\t0\t{name}\t{pos}\t60\t20M\t*\t0\t0\t{bases}\t{}\n",
                "I".repeat(20)
            ));
        }
    }

    fs::write(dir.join("in.sam"), sam)?;
    fs::write(dir.join("whole.bed"), "chr1\t0\t999999999\nchr2\t0\t2000\n")?;

    ntools(dir, &["faidx", "ref.fa"])?;
    ntools(dir, &["view", "--no-PG", "-b", "-o", "in.bam", "in.sam"])?;
    ntools(dir, &["index", "in.bam"])?;
    fs::copy(dir.join("in.bam"), dir.join("csi.bam"))?;
    ntools(dir, &["index", "-c", "csi.bam"])?;
    ntools(
        dir,
        &[
            "view", "--no-PG", "-T", "ref.fa", "-C", "-o", "in.cram", "in.sam",
        ],
    )?;
    ntools(dir, &["index", "in.cram"])?;

    Ok(())
}

/// As in samtools, a region end past the index's maximum position (2^29 for BAI), a
/// common way to say "the whole chromosome", is clamped instead of failing the command, and a
/// region that starts past it is empty.
#[test]
fn test_regions_past_the_index_maximum_are_clamped() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();
    write_inputs(dir)?;

    let count = |input: &str, region: &str| -> io::Result<String> {
        let mut args = vec!["view", "-c"];

        if input.ends_with(".cram") {
            args.extend(["-T", "ref.fa"]);
        }

        args.extend([input, region]);
        ntools(dir, &args).map(|out| out.trim().to_string())
    };

    for input in ["in.bam", "csi.bam", "in.cram"] {
        assert_eq!(count(input, "chr1")?, "41", "{input}");
        assert_eq!(count(input, "chr1:1-999999999")?, "41", "{input}");
        assert_eq!(count(input, "chr1:1000-999999999")?, "17", "{input}");
        assert_eq!(count(input, "chr1:600000000")?, "0", "{input}");
        assert_eq!(count(input, "chr1:600000000-700000000")?, "0", "{input}");
    }

    // A BED interval to 999999999.
    for input in ["in.bam", "csi.bam"] {
        let out = ntools(dir, &["view", "-c", "-L", "whole.bed", input])?;
        assert_eq!(out.trim(), "58", "{input}");

        let out = ntools(dir, &["bedcov", "whole.bed", input])?;
        assert_eq!(
            out, "chr1\t0\t999999999\t820\nchr2\t0\t2000\t340\n",
            "{input}"
        );

        let out = ntools(dir, &["depth", "-r", "chr1:1-999999999", input])?;
        assert_eq!(out.lines().count(), 820, "{input}");
    }

    Ok(())
}

/// A region query on an input without an index names the index files looked for and how to make
/// one, with samtools's exit statuses (2 for `bedcov`, 1 otherwise); an index named by replacing
/// the extension (`in.bai`) is found, as HTSlib finds it.
#[test]
fn test_missing_index_names_the_files_looked_for() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();
    write_inputs(dir)?;

    fs::copy(dir.join("in.bam"), dir.join("noidx.bam"))?;
    fs::copy(dir.join("in.cram"), dir.join("noidx.cram"))?;

    let bam = "ntools: no index for noidx.bam (looked for noidx.bam.bai, noidx.bai, noidx.bam.csi, \
               noidx.csi); run \"ntools index noidx.bam\"\n";
    let cram = "ntools: no index for noidx.cram (looked for noidx.cram.crai, noidx.crai); run \
                \"ntools index noidx.cram\"\n";

    let cases: [(&[&str], i32, &str); 5] = [
        (&["view", "noidx.bam", "chr1:1-100"], 1, bam),
        (&["depth", "-r", "chr1:1-100", "noidx.bam"], 1, bam),
        (&["coverage", "-r", "chr1:1-100", "noidx.bam"], 1, bam),
        (&["bedcov", "whole.bed", "noidx.bam"], 2, bam),
        (
            &["view", "-T", "ref.fa", "noidx.cram", "chr1:1-100"],
            1,
            cram,
        ),
    ];

    for (args, status, stderr) in cases {
        let output = Command::new(NTOOLS).current_dir(dir).args(args).output()?;
        assert_eq!(output.status.code(), Some(status), "ntools {args:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            stderr,
            "ntools {args:?}"
        );
    }

    fs::copy(dir.join("in.bam.bai"), dir.join("noidx.bai"))?;
    let count = ntools(dir, &["view", "-c", "noidx.bam", "chr1:1-100"])?;
    assert_eq!(count.trim(), "2");

    Ok(())
}
