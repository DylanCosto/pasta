//! `concat -l` compares the phase of the overlap by the alleles the genotype indexes name, so
//! that two chunks listing a site's ALT alleles in different orders (`C,G` with `1|2`, `G,C`
//! with `2|1`) agree. bcftools compares the indexes, takes such sites for a phase mismatch, and
//! flips the phase of everything after the overlap.

use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

const HEADER: &str = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "vtools-concat-ligate-{name}-{}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).expect("create the test directory");
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

/// Writes `records` as `NAME.vcf.gz` with a CSI index in `dir`.
fn chunk(dir: &Path, name: &str, records: &str) {
    fs::write(
        dir.join(format!("{name}.vcf")),
        format!("{HEADER}{records}"),
    )
    .expect("write the chunk");
    let made = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["view", "-Oz", "--write-index", "-o"])
        .arg(format!("{name}.vcf.gz"))
        .arg(format!("{name}.vcf"))
        .current_dir(dir)
        .status()
        .expect("run vtools");
    assert!(made.success(), "compress and index {name}: {made}");
}

#[test]
fn alt_alleles_in_another_order_keep_the_phase() {
    let t = TestDir::new("alt-order");
    let dir = t.path();
    chunk(
        dir,
        "a",
        "1\t10\t.\tA\tT\t.\t.\t.\tGT\t0|1\n\
         1\t20\t.\tA\tC,G\t.\t.\t.\tGT\t1|2\n\
         1\t25\t.\tA\tC,G\t.\t.\t.\tGT\t1|2\n\
         1\t30\t.\tA\tT\t.\t.\t.\tGT\t0|1\n",
    );
    chunk(
        dir,
        "b",
        "1\t20\t.\tA\tG,C\t.\t.\t.\tGT\t2|1\n\
         1\t25\t.\tA\tG,C\t.\t.\t.\tGT\t2|1\n\
         1\t30\t.\tA\tT\t.\t.\t.\tGT\t0|1\n\
         1\t40\t.\tA\tT\t.\t.\t.\tGT\t1|0\n",
    );
    let out = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["concat", "--no-version", "-l", "a.vcf.gz", "b.vcf.gz"])
        .current_dir(dir)
        .output()
        .expect("run vtools");
    assert!(out.status.success(), "{out:?}");
    let text = String::from_utf8_lossy(&out.stdout);
    let gts: Vec<(&str, &str)> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[1], f[9])
        })
        .collect();
    // The overlap agrees, so the second chunk's phase is kept, in the first chunk's phase set.
    assert_eq!(
        gts,
        [
            ("10", "0|1:10"),
            ("20", "1|2:10"),
            ("25", "1|2:10"),
            ("30", "0|1:99:10"),
            ("40", "1|0:10"),
        ],
        "{text}"
    );
}

/// The BGZF blocks of `data`, as (start, size).
fn blocks(data: &[u8]) -> Vec<(usize, usize)> {
    let mut blocks = Vec::new();
    let mut p = 0;
    while p < data.len() {
        let size = usize::from(u16::from_le_bytes([data[p + 16], data[p + 17]])) + 1;
        blocks.push((p, size));
        p += size;
    }
    blocks
}

#[test]
fn a_damaged_first_data_block_fails_instead_of_dropping_the_chunk() {
    let t = TestDir::new("damaged");
    let dir = t.path();
    chunk(
        dir,
        "a",
        "1\t10\t.\tA\tT\t.\t.\t.\tGT\t0|1\n1\t20\t.\tA\tT\t.\t.\t.\tGT\t1|0\n",
    );
    chunk(
        dir,
        "b",
        "1\t20\t.\tA\tT\t.\t.\t.\tGT\t1|0\n1\t30\t.\tA\tT\t.\t.\t.\tGT\t0|1\n",
    );
    chunk(dir, "empty", "");
    let data = fs::read(dir.join("b.vcf.gz")).expect("read the chunk");
    let blocks = blocks(&data);
    assert!(blocks.len() >= 3, "header, data and EOF blocks: {blocks:?}");
    // The first data block, with its CRC32 flipped, and the file cut inside it.
    let (start, size) = blocks[1];
    let mut bad = data.clone();
    bad[start + size - 8] ^= 0xff;
    fs::write(dir.join("bad.vcf.gz"), bad).expect("write the damaged chunk");
    fs::write(dir.join("cut.vcf.gz"), &data[..start + size / 2]).expect("write the cut chunk");
    for f in ["bad", "cut"] {
        fs::copy(
            dir.join("b.vcf.gz.csi"),
            dir.join(format!("{f}.vcf.gz.csi")),
        )
        .expect("copy the index");
        let name = format!("{f}.vcf.gz");
        let out = Command::new(env!("CARGO_BIN_EXE_vtools"))
            .args(["concat", "--no-version", "-l", "a.vcf.gz", &name])
            .current_dir(dir)
            .output()
            .expect("run vtools");
        assert_eq!(out.status.code(), Some(255), "{name}: {out:?}");
        assert!(
            String::from_utf8_lossy(&out.stderr).contains(&name),
            "{name}: {out:?}"
        );
    }
    // A chunk without records is still skipped.
    let out = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args([
            "concat",
            "--no-version",
            "-l",
            "a.vcf.gz",
            "empty.vcf.gz",
            "b.vcf.gz",
        ])
        .current_dir(dir)
        .output()
        .expect("run vtools");
    assert!(out.status.success(), "{out:?}");
    let text = String::from_utf8_lossy(&out.stdout);
    assert_eq!(
        text.lines().filter(|l| !l.starts_with('#')).count(),
        3,
        "{text}"
    );
}
