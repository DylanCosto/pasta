//! Decoding CRAM data series stored in the core data block with bit codecs.
//!
//! samtools writes these series to external blocks, so the fixtures in `data/cram_codecs` are
//! uncompressed CRAM 3.0 files written by samtools 1.24 whose encodings were then rewritten to
//! core-block codecs. samtools 1.24 decodes each one to the records in the matching SAM file.
//!
//! - `subexp.cram`: read lengths with SUBEXP (offset 0, k 0).
//! - `multi-subexp-k2.cram`: read lengths with SUBEXP (offset -1, k 2).
//! - `multi-gamma.cram`: read lengths with GAMMA (offset 2).
//! - `multi-beta.cram`: read lengths with BETA (offset -1, 8 bits).
//! - `multi-subexp-qsbeta.cram`: read lengths with SUBEXP (offset 0, k 4) and quality scores
//!   with byte BETA (offset 5, 7 bits), interleaved in the core data block.

use std::{error::Error, fs, path::PathBuf};

use pasta::cram;
use pasta::sam::{self as sam, alignment::RecordBuf};

fn data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_codecs")
        .join(name)
}

fn read_cram(name: &str) -> Result<Vec<RecordBuf>, Box<dyn Error>> {
    let src = fs::read(data_path(name))?;
    let mut reader = cram::io::Reader::new(&src[..]);
    let header = reader.read_header()?;
    let records = reader.records(&header).collect::<Result<_, _>>()?;
    Ok(records)
}

fn read_sam(name: &str) -> Result<Vec<RecordBuf>, Box<dyn Error>> {
    let src = fs::read(data_path(name))?;
    let mut reader = sam::io::Reader::new(&src[..]);
    let header = reader.read_header()?;
    let records = reader.record_bufs(&header).collect::<Result<_, _>>()?;
    Ok(records)
}

#[test]
fn test_decode_core_codecs() -> Result<(), Box<dyn Error>> {
    let cases = [
        ("subexp.cram", "subexp.sam"),
        ("multi-subexp-k2.cram", "multi.sam"),
        ("multi-gamma.cram", "multi.sam"),
        ("multi-beta.cram", "multi.sam"),
        ("multi-subexp-qsbeta.cram", "multi.sam"),
    ];

    for (cram_name, sam_name) in cases {
        let actual = read_cram(cram_name)?;
        let expected = read_sam(sam_name)?;
        assert!(!expected.is_empty(), "{sam_name}");
        assert_eq!(actual, expected, "{cram_name}");
    }

    Ok(())
}
