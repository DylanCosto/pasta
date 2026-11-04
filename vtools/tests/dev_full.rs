//! With standard output or standard error on /dev/full, the usage, version, listings, and
//! messages bcftools prints with printf/fprintf (ignoring write errors) end with bcftools'
//! status instead of a panic.
//!
//! Linux only: other systems have no /dev/full.
#![cfg(target_os = "linux")]

use std::{
    fs::{self, File},
    path::PathBuf,
    process::{Command, Stdio},
};

const VCF: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
    ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    chr1\t10\t.\tG\tT\t50\tPASS\tCSQ=T|missense_variant|ABC\n";

/// Which stream goes to /dev/full.
#[derive(Clone, Copy, Debug)]
enum Full {
    Stdout,
    Stderr,
}

fn input(name: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("vtools-dev-full-{name}-{}.vcf", std::process::id()));
    fs::write(&path, VCF).expect("write the input");
    path
}

/// Runs vtools with one stream on /dev/full and the other discarded; returns the status.
fn status(full: Full, args: &[&str]) -> Option<i32> {
    let dev_full = || {
        Stdio::from(
            File::options()
                .write(true)
                .open("/dev/full")
                .expect("/dev/full"),
        )
    };
    let (stdout, stderr) = match full {
        Full::Stdout => (dev_full(), Stdio::null()),
        Full::Stderr => (Stdio::null(), dev_full()),
    };
    Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(stderr)
        .status()
        .expect("run vtools")
        .code()
}

#[test]
fn usage_version_and_listings_ignore_a_full_stdout() {
    let vcf = input("stdout");
    let vcf = vcf.to_str().expect("a UTF-8 path");
    // The statuses bcftools 1.24 exits with.
    let cases: &[(&[&str], i32)] = &[
        (&["--version"], 0),
        (&["version"], 0),
        (&["-v"], 0),
        (&["--version-only"], 0),
        (&["help"], 0),
        (&["-h"], 0),
        (&["--help"], 0),
        (&["plugin", "-l"], 0),
        (&["plugin", "-lv"], 0),
        (&["+fill-tags", "-V"], 0),
        (&["+split-vep", "-l", vcf], 0),
        (&["csq", "--genetic-code", "l"], 0),
        (&["mpileup", "-X", "list"], 1),
    ];
    for (args, want) in cases {
        assert_eq!(
            status(Full::Stdout, args),
            Some(*want),
            "vtools {args:?} >/dev/full"
        );
    }
    let _ = fs::remove_file(vcf);
}

#[test]
fn messages_ignore_a_full_stderr() {
    let vcf = input("stderr");
    let vcf = vcf.to_str().expect("a UTF-8 path");
    let cases: &[(&[&str], i32)] = &[
        (&[], 1),
        (&["view", "/nonexistent/x.vcf"], 255),
        (&["view", "--nosuch"], 1),
        (&["help", "view"], 255),
        (&["nosuch"], 1),
        (&["mpileup", "-X", "list"], 1),
        (&["filter", "-i", "BAD", vcf], 255),
        // A warning on the way to a successful run.
        (&["view", "-H", "-i", "INFO/CSQ=\"a,b\"", vcf], 0),
    ];
    for (args, want) in cases {
        assert_eq!(
            status(Full::Stderr, args),
            Some(*want),
            "vtools {args:?} 2>/dev/full"
        );
    }
    let _ = fs::remove_file(vcf);
}
