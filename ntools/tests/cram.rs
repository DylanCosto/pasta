//! CRAM output and input against the reference: header checksums, the CRAM format options, and
//! many reference sequences.

use std::{
    fmt::Write as _,
    fs, io,
    path::Path,
    process::{Command, Output},
    time::{Duration, Instant},
};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

fn run(dir: &Path, args: &[&str]) -> io::Result<Output> {
    Command::new(NTOOLS).current_dir(dir).args(args).output()
}

/// Runs ntools in `dir`, failing the test if it fails, and returns its standard output.
fn ntools(dir: &Path, args: &[&str]) -> io::Result<String> {
    let output = run(dir, args)?;
    assert!(
        output.status.success(),
        "ntools {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).map_err(io::Error::other)
}

/// A pseudorandom sequence of `len` bases.
fn bases(seed: u64, len: usize) -> String {
    let mut state = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);

    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4] as char
        })
        .collect()
}

/// Writes `ref.fa` with sequences `c0`, `c1`, `c2` (120, 300, 120 bases) and returns them.
fn write_reference(dir: &Path) -> io::Result<Vec<String>> {
    let sequences: Vec<String> = [120, 300, 120]
        .into_iter()
        .enumerate()
        .map(|(i, len)| bases(i as u64, len))
        .collect();

    let mut fasta = String::new();

    for (i, sequence) in sequences.iter().enumerate() {
        writeln!(fasta, ">c{i}\n{sequence}").map_err(io::Error::other)?;
    }

    fs::write(dir.join("ref.fa"), fasta)?;

    Ok(sequences)
}

/// SAM text of three reads on each of `c0`, `c1`, and `c2` (sorted), with `c1`'s `@SQ` given
/// `c1_fields`, and the reads' `MD` and `NM`.
fn sam(sequences: &[String], c1_fields: &str) -> String {
    let mut sam = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c0\tLN:120\n@SQ\tSN:c1\tLN:300{c1_fields}\n\
         @SQ\tSN:c2\tLN:120\n"
    );

    for (i, sequence) in sequences.iter().enumerate() {
        for j in 1..=3 {
            let start = j * 10;
            let read = &sequence[start - 1..start - 1 + 30];
            sam.push_str(&format!(
                "r{i}_{j}\t0\tc{i}\t{start}\t60\t30M\t*\t0\t0\t{read}\t*\tMD:Z:30\tNM:i:0\n"
            ));
        }
    }

    sam
}

/// The records of SAM text, without the header.
fn records(sam: &str) -> Vec<&str> {
    sam.lines().filter(|line| !line.starts_with('@')).collect()
}

/// An `@SQ M5` that is not the checksum of the `-T` sequence the records are stored
/// against fails CRAM output, as in samtools, rather than writing a CRAM file that cannot be read
/// back with that reference.
#[test]
fn test_cram_output_fails_on_an_m5_that_does_not_match_the_reference() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let sequences = write_reference(dir.path())?;
    let input = sam(&sequences, "\tM5:00004142844c68362f611ca22888656a");
    fs::write(dir.path().join("in.sam"), &input)?;

    let output = run(
        dir.path(),
        &["view", "-C", "-T", "ref.fa", "-o", "out.cram", "in.sam"],
    )?;
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("SQ header M5 tag discrepancy for reference 'c1'"),
        "{stderr}"
    );

    // ignore_md5=1 writes it, with the reference's checksum, so that it reads back.
    ntools(
        dir.path(),
        &[
            "view",
            "-C",
            "-T",
            "ref.fa",
            "--output-fmt-option",
            "ignore_md5=1",
            "-o",
            "out.cram",
            "in.sam",
        ],
    )?;
    let header = ntools(dir.path(), &["view", "-H", "out.cram"])?;
    assert!(
        !header.contains("00004142844c68362f611ca22888656a"),
        "{header}"
    );
    let output = ntools(dir.path(), &["view", "-T", "ref.fa", "out.cram"])?;
    assert_eq!(records(&output), records(&input));

    // A matching M5 is kept, and an LN that differs from the reference gets HTSlib's warning.
    let md5 = ntools(dir.path(), &["view", "-H", "out.cram"])?
        .lines()
        .find(|line| line.starts_with("@SQ\tSN:c1\t"))
        .and_then(|line| line.split('\t').find(|field| field.starts_with("M5:")))
        .map(str::to_owned)
        .expect("an M5 for c1");
    let input = sam(&sequences, &format!("\t{md5}")).replace("SN:c2\tLN:120", "SN:c2\tLN:121");
    fs::write(dir.path().join("in.sam"), &input)?;
    let output = run(
        dir.path(),
        &["view", "-C", "-T", "ref.fa", "-o", "out.cram", "in.sam"],
    )?;
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("[W::sanitise_SQ_lines] Header @SQ length mismatch for ref c2, 121 vs 120"),
        "{stderr}"
    );

    Ok(())
}

/// `store_md=1` and `store_nm=1` keep the `MD` and `NM` tags that CRAM output otherwise
/// leaves out for readers to regenerate, as in samtools.
#[test]
fn test_cram_output_stores_md_and_nm() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let sequences = write_reference(dir.path())?;
    fs::write(dir.path().join("in.sam"), sam(&sequences, ""))?;

    let stored_tags = |options: &[&str]| -> io::Result<String> {
        let mut args = vec!["view", "-C", "-T", "ref.fa", "-o", "out.cram"];

        for option in options {
            args.extend(["--output-fmt-option", option]);
        }

        args.push("in.sam");
        ntools(dir.path(), &args)?;

        let output = ntools(
            dir.path(),
            &[
                "view",
                "--input-fmt-option",
                "decode_md=0",
                "-T",
                "ref.fa",
                "out.cram",
            ],
        )?;
        let first = records(&output).first().copied().unwrap_or_default();
        Ok(first.split('\t').skip(11).collect::<Vec<_>>().join("\t"))
    };

    assert_eq!(stored_tags(&[])?, "");
    assert_eq!(stored_tags(&["store_md=1"])?, "MD:Z:30");
    assert_eq!(stored_tags(&["store_nm=1"])?, "NM:i:0");
    assert_eq!(
        stored_tags(&["store_md=1", "store_nm=1"])?,
        "MD:Z:30\tNM:i:0"
    );
    assert_eq!(stored_tags(&["store_md=1", "store_md=0"])?, "");

    Ok(())
}

/// A slice whose reference checksum does not match `-T` names the region.
#[test]
fn test_cram_slice_checksum_mismatch_names_the_region() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let sequences = write_reference(dir.path())?;
    let read = &sequences[1][9..39];
    fs::write(
        dir.path().join("in.sam"),
        format!("@SQ\tSN:c1\tLN:300\nr1\t0\tc1\t10\t60\t30M\t*\t0\t0\t{read}\t*\n"),
    )?;
    ntools(
        dir.path(),
        &["view", "-C", "-T", "ref.fa", "-o", "out.cram", "in.sam"],
    )?;

    // c1 with other bases.
    let wrong: String = sequences[1]
        .chars()
        .map(|c| match c {
            'A' => 'C',
            'C' => 'A',
            c => c,
        })
        .collect();
    fs::write(dir.path().join("wrong.fa"), format!(">c1\n{wrong}\n"))?;

    let output = run(dir.path(), &["view", "-T", "wrong.fa", "out.cram"])?;
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("reference sequence checksum mismatch at c1:10-39"),
        "{stderr}"
    );
    assert!(!stderr.contains('['), "no byte arrays: {stderr}");

    Ok(())
}

/// CRAM output and input take time linear in the number of reference sequences
/// (formerly quadratic: 14 s each way for 50,000 contigs).
#[test]
fn test_cram_with_many_reference_sequences() -> io::Result<()> {
    const N: usize = 100_000;

    let dir = tempfile::tempdir()?;
    let mut fasta = String::new();
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n");

    for i in 0..N {
        writeln!(sam, "@SQ\tSN:ctg{i}\tLN:100").map_err(io::Error::other)?;
    }

    for i in 0..N {
        let sequence = bases(i as u64, 100);
        writeln!(fasta, ">ctg{i}\n{sequence}").map_err(io::Error::other)?;
        writeln!(
            sam,
            "r{i}\t0\tctg{i}\t11\t60\t50M\t*\t0\t0\t{}\t*",
            &sequence[10..60]
        )
        .map_err(io::Error::other)?;
    }

    fs::write(dir.path().join("ref.fa"), fasta)?;
    fs::write(dir.path().join("in.sam"), &sam)?;

    // Generous for unoptimized builds; the quadratic version took minutes.
    let limit = Duration::from_secs(60);

    let start = Instant::now();
    ntools(
        dir.path(),
        &["view", "-C", "-T", "ref.fa", "-o", "out.cram", "in.sam"],
    )?;
    let elapsed = start.elapsed();
    assert!(elapsed < limit, "writing took {elapsed:?}");

    let start = Instant::now();
    let output = ntools(dir.path(), &["view", "-T", "ref.fa", "out.cram"])?;
    let elapsed = start.elapsed();
    assert!(elapsed < limit, "reading took {elapsed:?}");

    assert_eq!(records(&output).len(), N);
    let last = records(&sam).last().copied().unwrap_or_default();
    // With the MD and NM tags that readers add.
    let expected = format!("{last}\tMD:Z:50\tNM:i:0");
    assert_eq!(records(&output).last().copied(), Some(expected.as_str()));

    Ok(())
}

/// Without `-T`, CRAM output of unaligned records needs no reference, and aligned
/// records take theirs from the local FASTA files `@SQ UR` names, as in samtools. Aligned records
/// with no reference available are stored without one, with a warning, as in samtools.
#[test]
fn test_cram_output_without_a_given_reference() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let sequences = write_reference(dir.path())?;
    ntools(dir.path(), &["faidx", "ref.fa"])?;

    fs::write(
        dir.path().join("unaligned.sam"),
        "@HD\tVN:1.6\tSO:unsorted\nq1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )?;
    ntools(
        dir.path(),
        &["view", "-C", "-o", "unaligned.cram", "unaligned.sam"],
    )?;
    let output = ntools(dir.path(), &["view", "unaligned.cram"])?;
    assert_eq!(records(&output), ["q1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"]);

    let reference = dir.path().join("ref.fa");
    let input = sam(&sequences, "").replace("\tLN:", &format!("\tUR:{}\tLN:", reference.display()));
    fs::write(dir.path().join("ur.sam"), &input)?;
    ntools(dir.path(), &["view", "-C", "-o", "ur.cram", "ur.sam"])?;
    let output = ntools(dir.path(), &["view", "ur.cram"])?;
    assert_eq!(records(&output), records(&input));

    let input = sam(&sequences, "");
    fs::write(dir.path().join("in.sam"), &input)?;
    let output = run(dir.path(), &["view", "-C", "-o", "out.cram", "in.sam"])?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert_eq!(stderr.lines().count(), 1, "{stderr}");
    assert!(stderr.contains("No reference for c0"), "{stderr}");
    let output = ntools(dir.path(), &["view", "out.cram"])?;
    assert_eq!(records(&output), records(&input));

    Ok(())
}

/// Reading a CRAM file whose `@SQ UR` reference is gone names what was tried and the fix, and
/// fails as samtools does.
#[test]
fn test_cram_input_without_its_reference() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let sequences = write_reference(dir.path())?;
    fs::write(dir.path().join("in.sam"), sam(&sequences, ""))?;
    ntools(
        dir.path(),
        &["view", "-C", "-T", "ref.fa", "-o", "in.cram", "in.sam"],
    )?;
    fs::remove_file(dir.path().join("ref.fa"))?;
    fs::remove_file(dir.path().join("ref.fa.fai"))?;

    let output = run(dir.path(), &["view", "-o", "out.bam", "in.cram"])?;
    assert_eq!(output.status.code(), Some(1));
    // As samtools does, the UR tag names the reference by its canonical path (on macOS the
    // temporary directory is under the symbolic link /var).
    let uri = dir.path().canonicalize()?.join("ref.fa");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        format!(
            "ntools: in.cram: no reference for c0 (tried UR:{}, not found); give one with \
             -T/--reference\n",
            uri.display()
        )
    );

    Ok(())
}

/// A command that fails before writing a record leaves no output file it created, as samtools
/// does, but never removes a file that was already there.
#[test]
fn test_failed_output_leaves_no_file() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let sequences = write_reference(dir.path())?;
    fs::write(dir.path().join("in.sam"), sam(&sequences, ""))?;
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         c0\t1\t.\tA\tC\t.\t.\t.\n",
    )?;

    for args in [
        &["view", "-b", "-o", "out.bam", "in.vcf"][..],
        &["view", "-C", "-T", "missing.fa", "-o", "out.bam", "in.sam"],
    ] {
        let output = run(dir.path(), args)?;
        assert!(!output.status.success(), "{args:?}");
        assert!(!dir.path().join("out.bam").exists(), "{args:?}");
    }

    fs::write(dir.path().join("out.bam"), "kept")?;
    let output = run(dir.path(), &["view", "-b", "-o", "out.bam", "in.vcf"])?;
    assert!(!output.status.success());
    assert!(dir.path().join("out.bam").exists());

    Ok(())
}
