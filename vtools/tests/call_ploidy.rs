//! `call` genotypes ploidy 0, 1 and 2 only: a sample with any other ploidy stops either caller
//! with an error and status 255 at the first record where it has it, instead of asserting, and
//! -c, which sizes its tables by the largest ploidy, needs that one to be 1 or 2. Ploidies are
//! bytes, as in bcftools.

use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

const VCF: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
    ##INFO=<ID=QS,Number=R,Type=Float,Description=\"Quality sums\">\n\
    ##INFO=<ID=I16,Number=16,Type=Float,Description=\"Counts\">\n\
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
    chr1\t5\t.\tA\tC\t.\t.\tQS=10,0;I16=5,5,0,0,300,9000,0,0,600,36000,0,0,100,2000,0,0\tPL\t0,90,100\t0,30,60\n\
    chr1\t10\t.\tA\tC\t.\t.\tQS=0,10;I16=0,0,5,5,0,0,300,9000,0,0,600,36000,0,0,100,2000\tPL\t100,100,0\t50,0,50\n";

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("vtools-call-ploidy-{name}-{}", std::process::id()));
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

/// Runs `vtools call CALLER --ploidy-file` with a file holding `ploidy`.
fn call_with_ploidy(dir: &TestDir, caller: &str, ploidy: &str) -> Output {
    call_with_samples(dir, caller, ploidy, None)
}

/// Runs `vtools call CALLER --ploidy-file`, with `-S` and a file holding `samples` if given.
fn call_with_samples(dir: &TestDir, caller: &str, ploidy: &str, samples: Option<&str>) -> Output {
    fs::write(dir.path().join("ploidy.txt"), ploidy).expect("write the ploidy file");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vtools"));
    cmd.args([
        "call",
        "--no-version",
        caller,
        "--ploidy-file",
        "ploidy.txt",
    ]);
    if let Some(samples) = samples {
        fs::write(dir.path().join("samples.txt"), samples).expect("write the samples file");
        cmd.args(["-S", "samples.txt"]);
    }
    cmd.arg("in.vcf")
        .current_dir(dir.path())
        .output()
        .expect("run vtools")
}

#[test]
fn unsupported_ploidies_fail() {
    let dir = TestDir::new("unsupported");
    let files = [
        "* * * * 3\n",
        "* * * * 4\n",
        // A byte of 255.
        "* * * * -1\n",
        // Only a region of the sequence, and not the default: the record at 5 is written first.
        "chr1 8 12 M 3\n* * * M 2\n",
    ];
    let consensus_only = [
        // -c computes every sample with the largest ploidy: bcftools reads past its tables
        // with 0, and computes nonsense with 3 although no sample is of sex X.
        "* * * * 0\n",
        "* * * X 3\n* * * M 1\n* * * F 2\n",
    ];
    let runs = files
        .iter()
        .flat_map(|&f| [("-m", f), ("-c", f)])
        .chain(consensus_only.iter().map(|&f| ("-c", f)));
    for (caller, file) in runs {
        let out = call_with_ploidy(&dir, caller, file);
        let err = String::from_utf8_lossy(&out.stderr);
        assert_eq!(
            out.status.code(),
            Some(255),
            "{caller} with {file:?}: {err}"
        );
        assert!(
            err.contains("Unsupported ploidy"),
            "{caller} with {file:?}: {err}"
        );
        if caller == "-m" {
            assert!(err.contains(" at chr1:"), "{caller} with {file:?}: {err}");
        }
        // The calls before the first record with the ploidy are kept, as after a failed read.
        let body = String::from_utf8_lossy(&out.stdout);
        let calls: Vec<&str> = body
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap_or_default())
            .collect();
        // -c stops before calling: it sizes its tables by the largest ploidy.
        let expected: &[&str] = if caller == "-m" && file.starts_with("chr1 8") {
            &["5"]
        } else {
            &[]
        };
        assert_eq!(calls, expected, "{caller} with {file:?}: {body}");
    }
}

/// The genotypes of the last record `call` wrote.
fn last_genotypes(out: &Output) -> Vec<String> {
    let body = String::from_utf8_lossy(&out.stdout);
    let last = body.lines().last().unwrap_or_default();
    last.split('\t')
        .skip(9)
        .map(|s| s.split(':').next().unwrap_or_default().to_owned())
        .collect()
}

#[test]
fn supported_ploidies_are_called() {
    let dir = TestDir::new("supported");
    // As bcftools calls them.
    for (file, caller, gts) in [
        ("* * * * 0\n", "-m", [".", "."]),
        ("* * * * 1\n", "-m", ["1", "1"]),
        ("* * * * 2\n", "-m", ["1/1", "1/1"]),
        // Every sample is of the last sex declared, so none is triploid.
        ("* * * X 3\n* * * M 1\n* * * F 2\n", "-m", ["1/1", "1/1"]),
        // bcftools takes a ploidy modulo 256, as a byte.
        ("* * * * 258\n", "-m", ["1/1", "1/1"]),
        ("* * * * 4294967298\n", "-m", ["1/1", "1/1"]),
        ("* * * * 1\n", "-c", ["1", "1"]),
        ("* * * * 2\n", "-c", ["1/1", "0/1"]),
        ("* * * * 258\n", "-c", ["1/1", "0/1"]),
    ] {
        let out = call_with_ploidy(&dir, caller, file);
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(out.status.success(), "{caller} with {file:?}: {err}");
        assert_eq!(last_genotypes(&out), gts, "{caller} with {file:?}");
    }
}

#[test]
fn unused_unsupported_ploidies_are_called() {
    let dir = TestDir::new("unused");
    // As bcftools calls them.
    for (file, samples, gts) in [
        // The samples of a set ploidy keep the largest one, 3, until the ploidy of a sex
        // changes, which F's does at the first record, from 3 to 2.
        ("* * * M 3\n* * * F 2\n", "S1 2\nS2 1\n", ["1/1", "1"]),
        // A region on a sequence the input does not have.
        (
            "chrZ 1 100 F 3\n* * * F 2\n* * * M 1\n",
            "S1 F\nS2 F\n",
            ["1/1", "1/1"],
        ),
    ] {
        let out = call_with_samples(&dir, "-m", file, Some(samples));
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(out.status.success(), "{file:?} and {samples:?}: {err}");
        assert_eq!(last_genotypes(&out), gts, "{file:?} and {samples:?}");
    }
}
