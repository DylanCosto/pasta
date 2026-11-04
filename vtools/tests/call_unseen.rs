//! `call -m` at a site whose best genotypes include the unseen allele `<*>`, which is not kept:
//! bcftools writes genotypes that point at it (`1/1` with no ALT) and NaN GP. vtools writes the
//! dropped allele missing (`./.`, `0/.`), counts only the alleles it writes in AC and AN, and
//! writes GQ and GP as for a sample without data.

use std::process::{Command, Stdio};

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
    ##INFO=<ID=QS,Number=R,Type=Float,Description=\"Quality sums\">\n\
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
    ##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Quality sums\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
    chr1\t16\t.\tG\tT,<*>\t0\t.\tQS=0.256778,0.358242,0.38498\tPL:QS\t\
    204,145,59,39,97,0:146,309,32\t146,127,229,6,0,45:272,457,247\t109,4,133,0,13,128:370,34,476\n";

fn call(args: &[&str]) -> Vec<String> {
    let mut child = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["call", "--no-version", "-m", "-a", "GQ,GP"])
        .args(args)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("run vtools");
    std::io::Write::write_all(&mut child.stdin.take().expect("stdin"), VCF.as_bytes())
        .expect("write the input");
    let out = child.wait_with_output().expect("wait for vtools");
    assert!(out.status.success(), "{out:?}");
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

#[test]
fn a_genotype_of_the_dropped_unseen_allele_is_missing() {
    let lines = call(&[]);
    assert_eq!(lines.len(), 1, "{lines:?}");
    let f: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(f[4], ".", "{lines:?}");
    assert!(f[7].contains("AN=2"), "{lines:?}");
    let gt = |s: &str| s.split(':').next().unwrap_or_default().to_owned();
    assert_eq!(
        [gt(f[9]), gt(f[10]), gt(f[11])],
        ["./.", "0/.", "0/."],
        "{lines:?}"
    );
    assert!(!lines[0].contains("nan"), "{lines:?}");
}

#[test]
fn a_haploid_genotype_of_the_dropped_unseen_allele_is_missing() {
    let lines = call(&["--ploidy", "1"]);
    let f: Vec<&str> = lines[0].split('\t').collect();
    let gt = |s: &str| s.split(':').next().unwrap_or_default().to_owned();
    assert_eq!(
        [gt(f[9]), gt(f[10]), gt(f[11])],
        [".", ".", "0"],
        "{lines:?}"
    );
    assert!(f[7].contains("AN=1"), "{lines:?}");
    assert!(!lines[0].contains("nan"), "{lines:?}");
}
