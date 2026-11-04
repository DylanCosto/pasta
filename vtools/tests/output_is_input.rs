//! `-o` naming an input is refused before the output is opened, leaving the input untouched.
//! bcftools truncates the input first, and the commands then read an empty file.

use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "vtools-output-is-input-{name}-{}",
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

#[test]
fn output_naming_an_input_is_refused() {
    let t = TestDir::new("refused");
    let dir = t.path();
    fs::write(
        dir.join("in.vcf"),
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t10\t.\tA\tT\t30\t.\t.\n1\t20\t.\tA\tC,G\t40\t.\t.\n",
    )
    .expect("write the VCF");
    for (flag, f) in [("-Oz", "in.vcf.gz"), ("-Ob", "in.bcf")] {
        let made = Command::new(env!("CARGO_BIN_EXE_vtools"))
            .args(["view", flag, "-o", f, "in.vcf"])
            .current_dir(dir)
            .status()
            .expect("run vtools");
        assert!(made.success(), "{f}: {made}");
    }
    let commands: &[&[&str]] = &[
        &["view", "-i", "QUAL>20"],
        &["norm", "-m-any"],
        &["filter", "-i", "QUAL>20"],
        &["annotate", "-x", "INFO"],
        &["concat"],
    ];
    for f in ["in.vcf.gz", "in.bcf"] {
        let before = fs::read(dir.join(f)).expect("read the input");
        for args in commands {
            let out = Command::new(env!("CARGO_BIN_EXE_vtools"))
                .args(*args)
                .args(["-o", f, f])
                .current_dir(dir)
                .output()
                .expect("run vtools");
            assert!(!out.status.success(), "{args:?} {f}: {out:?}");
            assert!(
                String::from_utf8_lossy(&out.stderr).contains("also an input"),
                "{args:?} {f}: {out:?}"
            );
            assert_eq!(
                fs::read(dir.join(f)).expect("read the input"),
                before,
                "{args:?} {f}"
            );
        }
    }
}
