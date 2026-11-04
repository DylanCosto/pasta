//! CSI and TBI files written as bcftools 1.24 writes them.

use std::path::PathBuf;

use pasta::variant::index::{Index, IndexFormat};

/// The uncompressed sizes of a BGZF file's blocks, from their headers and footers.
fn block_sizes(data: &[u8]) -> Vec<u32> {
    let mut sizes = Vec::new();
    let mut rest = data;
    while !rest.is_empty() {
        let bsize = usize::from(u16::from_le_bytes([rest[16], rest[17]])) + 1;
        let isize = &rest[bsize - 4..bsize];
        sizes.push(u32::from_le_bytes(isize.try_into().expect("four bytes")));
        rest = &rest[bsize..];
    }
    sizes
}

/// Reads an index bcftools wrote, saves it again, and returns the bytes of both files.
fn save_again(name: &str, fmt: IndexFormat) -> (Vec<u8>, Vec<u8>) {
    let golden = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/index")
        .join(name);
    let idx = Index::read(&golden).expect("the index reads");
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    idx.save(&path, fmt).expect("the index is saved");
    let got = std::fs::read(&path).expect("the saved index is read");
    let want = std::fs::read(&golden).expect("the index is read");
    (got, want)
}

#[test]
fn saved_indexes_have_htslib_blocks() {
    // A record at 400 Mb gives a linear index of 191 KiB: three blocks of up to 0xff00 bytes.
    let (got, want) = save_again("far.vcf.gz.tbi", IndexFormat::Tbi);
    assert_eq!(block_sizes(&want), [0xff00, 0xff00, 64902, 0]);
    assert_eq!(block_sizes(&got), block_sizes(&want));
}

/// bcftools is built with libdeflate, whose output pasta matches only with libdeflate too.
#[test]
#[cfg(feature = "libdeflate")]
fn saved_indexes_have_htslib_bytes() {
    // 50 records over 400 Mb, enough data for the compression level to show.
    for (name, fmt) in [
        ("scattered.vcf.gz.csi", IndexFormat::Csi),
        ("scattered.vcf.gz.tbi", IndexFormat::Tbi),
        ("far.vcf.gz.tbi", IndexFormat::Tbi),
    ] {
        let (got, want) = save_again(name, fmt);
        assert!(
            got == want,
            "{name}: {} bytes, want {}",
            got.len(),
            want.len()
        );
    }
}
