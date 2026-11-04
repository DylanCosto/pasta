//! SAM records read as HTSlib's `sam_parse1` normalizes them, and SAM.gz unmapped
//! queries.

use std::{
    fs,
    io::{self, Write as _},
    path::PathBuf,
};

use pasta::bam;
use pasta::bgzf;
use pasta::sam::{
    self as sam,
    alignment::{
        io::{Read as _, Write as _},
        record::Flags,
    },
};
use pasta::util::alignment;

/// The example records with what samtools 1.24 prints for them, except where pasta
/// differs on purpose (see `sam::io::reader::normalize`).
const CASES: &[(&str, &str)] = &[
    // RNAME `*` without 0x4
    ("a\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", "4\t*\t*"),
    ("b\t65\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", "69\t*\t*"),
    // placed, CIGAR `*`
    ("c\t0\tc1\t10\t30\t*\t*\t0\t0\tACGT\tIIII", "4\tc1\t*"),
    // RNAME with POS 0
    ("d\t0\tc1\t0\t0\t4M\t*\t0\t0\tACGT\tIIII", "4\t*\t*"),
    // normal
    ("e\t0\tc1\t20\t30\t4M\t*\t0\t0\tACGT\tIIII", "0\tc1\t*"),
    // RNEXT with PNEXT 0
    ("f\t1\tc1\t30\t30\t4M\tc1\t0\t0\tACGT\tIIII", "1\tc1\t*"),
    // RNEXT `=` of a record with POS 0: the mate stays on c1 (HTSlib: `*`)
    ("g\t1\tc1\t0\t30\t4M\t=\t40\t0\tACGT\tIIII", "5\t*\tc1\t40"),
    // unknown RNAME
    // (POS is dropped with RNAME; HTSlib keeps 10)
    (
        "x\t0\tchrX\t10\t30\t4M\t*\t0\t0\tACGT\tIIII",
        "4\t*\t*\t0\t0",
    ),
    // unknown RNEXT
    // (PNEXT is dropped with RNEXT; HTSlib keeps 40)
    (
        "i\t1\tc1\t50\t30\t4M\tchrY\t40\t0\tACGT\tIIII",
        "1\tc1\t*\t0",
    ),
    // unknown RNAME with RNEXT `=`
    (
        "h\t1\tchrX\t50\t30\t4M\t=\t40\t0\tACGT\tIIII",
        "5\t*\t*\t0\t0",
    ),
    // `@SQ AN` alternative names
    (
        "j\t1\tone\t50\t30\t4M\t1\t40\t0\tACGT\tIIII",
        "1\tc1\tc1\t40",
    ),
];

const HEADER: &str = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\tAN:1,one\n";

fn input() -> String {
    let mut s = String::from(HEADER);

    for (line, _) in CASES {
        s.push_str(line);
        s.push('\n');
    }

    s
}

/// FLAG, RNAME and RNEXT; with 4 fields, also PNEXT; with 5, also POS (before PNEXT).
fn placement(record: &sam::Record, fields: usize) -> io::Result<String> {
    let name = |n: Option<&bstr::BStr>| n.map_or_else(|| "*".to_string(), |n| n.to_string());
    let position = |p: Option<io::Result<pasta::core::Position>>| {
        p.transpose().map(|p| p.map_or(0, usize::from).to_string())
    };

    let mut s = format!(
        "{}\t{}\t{}",
        u16::from(record.flags()?),
        name(record.reference_sequence_name()),
        name(record.mate_reference_sequence_name()),
    );

    if fields == 5 {
        s = format!("{s}\t{}", position(record.alignment_start())?);
    }

    if fields >= 4 {
        s = format!("{s}\t{}", position(record.mate_alignment_start())?);
    }

    Ok(s)
}

#[test]
fn test_read_normalized_record() -> io::Result<()> {
    let src = input();
    let mut reader = sam::io::Reader::new(src.as_bytes());
    let header = reader.read_header()?;

    let mut record = sam::Record::default();

    for (line, expected) in CASES {
        reader.read_normalized_record(&header, &mut record)?;
        let fields = expected.split('\t').count();
        assert_eq!(placement(&record, fields)?, *expected, "{line}");
    }

    assert_eq!(reader.read_normalized_record(&header, &mut record)?, 0);

    Ok(())
}

#[test]
fn test_read_normalized_record_writes_htslib_bam() -> io::Result<()> {
    let src = input();
    let mut reader = sam::io::Reader::new(src.as_bytes());
    let header = reader.read_header()?;

    let mut writer = bam::io::Writer::from(Vec::new());
    let mut record = sam::Record::default();
    let mut mapped = 0;

    while reader.read_normalized_record(&header, &mut record)? != 0 {
        writer.write_alignment_record(&header, &record)?;

        if !record.flags()?.is_unmapped() {
            mapped += 1;
        }
    }

    // samtools flagstat: 4 mapped (e, f, i, j).
    assert_eq!(mapped, 4);

    let src = writer.get_ref().clone();
    let mut reader = bam::io::Reader::from(&src[..]);
    let records: Vec<_> = reader.records().collect::<io::Result<_>>()?;

    let x = &records[7];
    assert!(x.reference_sequence_id().is_none());
    assert!(x.flags().is_unmapped());
    // The position is dropped with the reference sequence (HTSlib keeps 10).
    assert!(x.alignment_start().is_none());

    Ok(())
}

#[test]
fn test_read_record_buf_and_alignment_records_are_normalized() -> io::Result<()> {
    let src = input();

    let mut reader = sam::io::Reader::new(src.as_bytes());
    let header = reader.read_header()?;
    let flags: Vec<_> = reader
        .record_bufs(&header)
        .map(|result| result.map(|record| u16::from(record.flags())))
        .collect::<io::Result<_>>()?;

    let expected: Vec<u16> = CASES
        .iter()
        .map(|(_, e)| e.split('\t').next().unwrap().parse().unwrap())
        .collect();

    assert_eq!(flags, expected);

    let mut reader = sam::io::Reader::new(src.as_bytes());
    let header = reader.read_alignment_header()?;
    let flags: Vec<_> = reader
        .alignment_records(&header)
        .map(|result| result.and_then(|record| record.flags().map(u16::from)))
        .collect::<io::Result<_>>()?;

    assert_eq!(flags, expected);

    Ok(())
}

#[test]
fn test_read_normalized_record_with_no_reference_sequences() {
    let src = "r\t0\tc1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    let mut reader = sam::io::Reader::new(src.as_bytes());
    let header = reader.read_header().unwrap();
    let mut record = sam::Record::default();

    // samtools: "no SQ lines present in the header".
    assert!(reader.read_normalized_record(&header, &mut record).is_err());

    let src = "r\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
    let mut reader = sam::io::Reader::new(src.as_bytes());
    let header = reader.read_header().unwrap();
    assert!(reader.read_normalized_record(&header, &mut record).is_ok());
    assert_eq!(record.flags().unwrap(), Flags::UNMAPPED);
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(name: &str) -> io::Result<Self> {
        let path = std::env::temp_dir().join(format!("pasta-{name}-{}", std::process::id()));
        fs::create_dir_all(&path)?;
        Ok(Self(path))
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

// A placed unmapped mate on the last reference sequence is not unplaced.
#[test]
fn test_sam_gz_query_unmapped_returns_only_unplaced_records() -> io::Result<()> {
    let dir = TempDir::new("psam-002")?;
    let src = dir.0.join("x.sam.gz");

    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
                pair\t73\tchr1\t100\t60\t10M\t=\t100\t0\tACGTACGTAC\tIIIIIIIIII\n\
                pair\t133\tchr1\t100\t0\t*\t=\t100\t0\tACGTACGTAC\tIIIIIIIIII\n\
                unplaced\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";

    let mut writer = bgzf::io::Writer::new(fs::File::create(&src)?);
    writer.write_all(text.as_bytes())?;
    writer.finish()?;

    let index = sam::fs::index(&src)?;
    pasta::csi::fs::write(dir.0.join("x.sam.gz.csi"), &index)?;

    let mut reader = sam::io::indexed_reader::Builder::default().build_from_path(&src)?;
    reader.read_header()?;
    let names: Vec<_> = reader
        .query_unmapped()?
        .map(|result| result.map(|record| record.name().map(|n| n.to_string())))
        .collect::<io::Result<_>>()?;

    assert_eq!(names, [Some(String::from("unplaced"))]);

    let mut reader = alignment::io::indexed_reader::Builder::default().build_from_path(&src)?;
    let header = reader.read_header()?;
    let names: Vec<_> = reader
        .query_unmapped(&header)?
        .map(|result| result.map(|record| record.name().map(|n| n.to_string())))
        .collect::<io::Result<_>>()?;

    assert_eq!(names, [Some(String::from("unplaced"))]);

    Ok(())
}
