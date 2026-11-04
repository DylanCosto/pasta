//! `gtcheck --distinctive-sites NUM,MEM,TMP` holds at most MEM bytes of sites and writes the
//! others to a temporary file named from TMP, which gives the same report.

use std::{
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
};

/// Four samples and sites where different pairs differ.
const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
    chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/0\t1/1\t0/0\t0/1\n\
    chr1\t20\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t1/1\t0/0\n\
    chr1\t30\t.\tA\tG\t.\t.\t.\tGT\t1/1\t0/0\t0/0\t1/1\n\
    chr1\t40\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t1/1\t1/1\n\
    chr1\t50\t.\tA\tG\t.\t.\t.\tGT\t0/0\t1/1\t1/1\t0/0\n\
    chr1\t60\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\t0/0\t0/0\n";

const PAIRS: &str = "S1\tS2\nS1\tS3\nS1\tS4\nS2\tS3\nS2\tS4\nS3\tS4\n";

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        Self::with_input(name, VCF, PAIRS)
    }

    fn with_input(name: &str, vcf: &str, pairs: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "vtools-gtcheck-distinctive-{name}-{}",
            std::process::id()
        ));
        fs::create_dir_all(dir.join("tmp")).expect("create the test directory");
        fs::write(dir.join("in.vcf"), vcf).expect("write the input");
        fs::write(dir.join("pairs.txt"), pairs).expect("write the pairs");
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

fn gtcheck(dir: &TestDir, distinctive: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args([
            "gtcheck",
            "-P",
            "pairs.txt",
            "--distinctive-sites",
            distinctive,
            "in.vcf",
        ])
        .current_dir(dir.path())
        .output()
        .expect("run vtools")
}

/// The DS lines of a report.
fn ds_lines(out: &Output) -> Vec<String> {
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "gtcheck: {err}");
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| l.starts_with("DS\t"))
        .map(str::to_owned)
        .collect()
}

#[test]
fn sites_past_mem_go_to_a_temporary_file() {
    let dir = TestDir::new("spill");
    let in_memory = ds_lines(&gtcheck(&dir, "1"));
    assert_eq!(in_memory.len(), 6, "{in_memory:?}");
    // One byte holds one site: every other one is written out.
    let spilled = ds_lines(&gtcheck(&dir, "1,1,tmp/ds"));
    assert_eq!(spilled, in_memory);
    let left = fs::read_dir(dir.path().join("tmp"))
        .expect("read the temporary directory")
        .count();
    assert_eq!(left, 0, "the temporary file is removed");
}

#[test]
fn an_unusable_tmp_fails_only_when_the_sites_exceed_mem() {
    let dir = TestDir::new("no-tmp");
    let out = gtcheck(&dir, "1,1,nosuchdir/ds");
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(255), "{err}");
    assert!(err.contains("temporary file"), "{err}");
    // The sites fit in 500 MB, so no temporary file is needed (bcftools always makes one).
    assert_eq!(ds_lines(&gtcheck(&dir, "1,500M,nosuchdir/ds")).len(), 6);
}

/// `nsamples` samples at `nsites` sites of random genotypes, and every pair of them.
fn many_samples(nsamples: usize, nsites: usize) -> (String, String) {
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for i in 0..nsamples {
        vcf.push_str(&format!("\tS{i}"));
    }
    vcf.push('\n');
    let mut x = 1u64;
    for pos in 1..=nsites {
        vcf.push_str(&format!("chr1\t{pos}\t.\tA\tG\t.\t.\t.\tGT"));
        for _ in 0..nsamples {
            x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            vcf.push_str(["\t0/0", "\t0/1", "\t1/1"][(x >> 33) as usize % 3]);
        }
        vcf.push('\n');
    }
    let mut pairs = String::new();
    for i in 0..nsamples {
        for j in i + 1..nsamples {
            pairs.push_str(&format!("S{i}\tS{j}\n"));
        }
    }
    (vcf, pairs)
}

#[test]
fn many_runs_are_merged_in_passes() {
    // 19,900 pairs make records of 2,512 bytes, of which 417 runs are merged at a time; one
    // byte of MEM makes a run of each of the 600 sites, so they are merged in two passes.
    let (vcf, pairs) = many_samples(200, 600);
    let dir = TestDir::with_input("passes", &vcf, &pairs);
    let in_memory = ds_lines(&gtcheck(&dir, "1"));
    assert!(in_memory.len() > 100, "{}", in_memory.len());
    for mem in ["1", "64k"] {
        let spilled = ds_lines(&gtcheck(&dir, &format!("1,{mem},tmp/ds")));
        assert!(spilled == in_memory, "MEM {mem} gives the in-memory report");
    }
    let left = fs::read_dir(dir.path().join("tmp"))
        .expect("read the temporary directory")
        .count();
    assert_eq!(left, 0, "the temporary files are removed");
}
