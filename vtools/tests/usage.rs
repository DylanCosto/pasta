//! Exit statuses of command lines that bcftools rejects with its usage.

use std::process::{Command, Output, Stdio};

fn vtools(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .stdin(Stdio::null())
        .output()
        .expect("run vtools")
}

/// `--` before any plugin name ends bcftools' options with no plugin named: the plugin
/// interface's usage, with status 1, not the word after `--` taken as the plugin's name.
#[test]
fn plugin_with_double_dash_before_the_name_prints_usage() {
    for args in [
        &["plugin", "--", "fill-tags"][..],
        &["plugin", "--", "no-such-plugin"],
        &["plugin", "--"],
    ] {
        let out = vtools(args);
        assert_eq!(out.status.code(), Some(1), "vtools {args:?}");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains("Run user defined plugin"), "{stderr}");
        assert!(out.stdout.is_empty(), "vtools {args:?}");
    }
}

/// An empty `-O` is not an output type. bcftools reads a compression level from past the end
/// of the empty value (the next argument) and fails with status 255 unless that argument is
/// empty or a digit; vtools always fails, before any other check (mpileup's `-f`).
#[test]
fn empty_output_type_is_rejected() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("usage-empty-output-type");
    std::fs::create_dir_all(&dir).expect("create the test directory");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
         1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n",
    )
    .expect("write the VCF");
    let vcf = vcf.to_str().expect("UTF-8 path");
    for args in [
        &["gtcheck", "-O", "", vcf][..],
        &["mpileup", "-O", "", vcf],
        &["view", "-O", "", vcf],
    ] {
        let out = vtools(args);
        assert_eq!(out.status.code(), Some(255), "vtools {args:?}");
        assert!(out.stdout.is_empty(), "vtools {args:?}");
    }
}

/// An unrecognized command fails as in bcftools, with a hint naming the closest command.
#[test]
fn unrecognized_command_suggests_a_similar_one() {
    let out = vtools(&["veiw"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&out.stderr),
        "[E::main] unrecognized command 'veiw'\n  tip: a similar command exists: 'view'\n"
    );

    let out = vtools(&["xyzzy"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&out.stderr),
        "[E::main] unrecognized command 'xyzzy'\n"
    );
}
