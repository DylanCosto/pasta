//! VCF parse/format and BCF round trips against outputs of bcftools 1.24.

use std::io::Cursor;

use pasta::variant::{
    Header, Record, Unpack, bcf,
    io::{OutputType, WriteOptions, Writer},
    log, vcf,
};

/// Edge cases: multi-sample, multi-allelic, missing values, Number=A/R/G, flags, long floats,
/// haploid and mixed ploidy, phased GT, symbolic alleles, breakends, undefined tags.
const INPUT: &str = include_str!("data/variant_edge.vcf");
/// `bcftools view --no-version` of [`INPUT`].
const EXPECTED: &str = include_str!("data/variant_edge.expected.vcf");

fn split_header(text: &str) -> (&str, Vec<&str>) {
    let body = text.find("\n#CHROM").expect("a #CHROM line");
    let end = text[body + 1..].find('\n').expect("a newline") + body + 2;
    (&text[..end], text[end..].lines().collect())
}

fn parse_all(header: &mut Header, lines: &[&str]) -> Vec<Record> {
    lines
        .iter()
        .map(|l| {
            let mut rec = Record::new();
            let mut buf = l.as_bytes().to_vec();
            let ((), _) = log::capture(|| {
                vcf::parse(header, &mut buf, &mut rec).expect("the line parses");
            });
            rec
        })
        .collect()
}

#[test]
fn vcf_parse_format_matches_bcftools() {
    let (header_text, lines) = split_header(INPUT);
    let mut header = Header::parse(header_text.as_bytes()).expect("the header parses");
    let records = parse_all(&mut header, &lines);
    let (_, expected) = split_header(EXPECTED);
    for (rec, want) in records.iter().zip(&expected) {
        let mut out = Vec::new();
        vcf::format_record(&header, rec, &mut out).expect("the record formats");
        assert_eq!(String::from_utf8_lossy(&out), format!("{want}\n"));
    }
}

#[test]
fn bcf_round_trip_preserves_text() {
    let (header_text, lines) = split_header(INPUT);
    let mut header = Header::parse(header_text.as_bytes()).expect("the header parses");
    let mut records = parse_all(&mut header, &lines);
    let mut data = Vec::new();
    for rec in &mut records {
        rec.clear_errors();
        bcf::encode_record(&header, rec, &mut data).expect("the record encodes");
    }
    let (_, expected) = split_header(EXPECTED);
    let mut cursor = Cursor::new(&data[..]);
    for want in expected {
        let mut rec = Record::new();
        let outcome = bcf::read_record(&mut cursor, &header, &mut rec).expect("the record reads");
        assert_eq!(outcome, bcf::ReadOutcome::Record);
        // Unpacked and packed formatting agree.
        for unpack in [Unpack::empty(), Unpack::ALL] {
            rec.unpack(unpack);
            let mut out = Vec::new();
            vcf::format_record(&header, &rec, &mut out).expect("the record formats");
            assert_eq!(String::from_utf8_lossy(&out), format!("{want}\n"));
        }
    }
}

#[test]
fn edits_sync_like_htslib() {
    let (header_text, lines) = split_header(INPUT);
    let mut header = Header::parse(header_text.as_bytes()).expect("the header parses");
    let mut rec = parse_all(&mut header, &lines[..1]).remove(0);
    rec.update_info_i32(&header, b"DP", &[12345])
        .expect("DP is defined");
    rec.update_info_flag(&header, b"DB", false)
        .expect("DB is defined");
    rec.update_id(Some(b"rs9"));
    rec.update_filter(&[header.id(b"q10").expect("q10 is defined") as i32]);
    rec.update_format_i32(&header, b"GQ", &[1, 2, 3])
        .expect("GQ is defined");
    let mut text = Vec::new();
    vcf::format_record(&header, &rec, &mut text).expect("the record formats");
    let want = "1\t100\trs9\tA\tC,G\t30.1235\tq10\tDP=12345;AF=0.123457,0.5;AD=1,2,3\t\
GT:GQ:PL:AD\t0/1:1:0,10,100,20,200,300:1,2,3\t1|2:2:.:.\t./.:3:1,2:.\n";
    assert_eq!(String::from_utf8_lossy(&text), want);
    let mut data = Vec::new();
    bcf::encode_record(&header, &mut rec, &mut data).expect("the record encodes");
    let mut back = Record::new();
    bcf::read_record(&mut &data[..], &header, &mut back).expect("the record reads");
    text.clear();
    vcf::format_record(&header, &back, &mut text).expect("the record formats");
    assert_eq!(String::from_utf8_lossy(&text), want);
}

#[test]
fn writer_round_trips_through_bgzf() {
    let (header_text, lines) = split_header(INPUT);
    let mut header = Header::parse(header_text.as_bytes()).expect("the header parses");
    let mut records = parse_all(&mut header, &lines);
    let path = std::env::temp_dir().join(format!("pasta-variant-{}.vcf.gz", std::process::id()));
    let options = WriteOptions {
        output_type: OutputType::VcfBgzf,
        level: None,
        threads: 2,
    };
    let mut writer = Writer::create(Some(&path), &options).expect("the file is created");
    writer.write_header(&header).expect("the header is written");
    for rec in &mut records {
        writer
            .write_record(&header, rec)
            .expect("the record is written");
    }
    writer.finish().expect("the file is finished");
    let mut reader = pasta::variant::io::Reader::open(path.to_str().expect("a UTF-8 path"), 0)
        .expect("the file opens");
    let mut header2 = reader.read_header().expect("the header reads");
    let mut rec = Record::new();
    let mut line = Vec::new();
    let mut out = header2.to_vcf_text();
    while reader
        .read_record(&mut header2, &mut line, &mut rec)
        .expect("the records read")
    {
        vcf::format_record(&header2, &rec, &mut out).expect("the record formats");
    }
    let _ = std::fs::remove_file(&path);
    let mut want = header.to_vcf_text();
    for l in split_header(EXPECTED).1 {
        want.extend_from_slice(l.as_bytes());
        want.push(b'\n');
    }
    assert_eq!(
        String::from_utf8_lossy(&out),
        String::from_utf8_lossy(&want)
    );
}

#[test]
fn diploid_gt_only_lines_encode_like_the_general_parser() {
    let header_text = "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";
    let mut header = Header::parse(header_text.as_bytes()).expect("the header parses");
    // The first line takes the one-pass path; the others fall back to the general steps.
    let lines = [
        "1\t1\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1/0\t9|9",
        "1\t2\t.\tA\tC\t.\t.\t.\tGT\t0|1\t./.\t9|9",
        "1\t3\t.\tA\tC\t.\t.\t.\tGT\t0|1\t10|1\t1/1",
        "1\t4\t.\tA\tC\t.\t.\t.\tGT\t0\t1|0\t1/1",
        "1\t5\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1|0",
    ];
    let records = parse_all(&mut header, &lines[..4]);
    // GT is ID 1 (after PASS): int8 key 1, two int8 values per sample.
    assert_eq!(
        records[0].indiv_bytes(),
        [0x11, 1, 0x21, 3, 5, 4, 2, 21, 21]
    );
    assert_eq!(
        records[1].indiv_bytes(),
        [0x11, 1, 0x21, 3, 5, 0, 0, 21, 21]
    );
    for (rec, line) in records.iter().zip(lines) {
        let mut out = Vec::new();
        vcf::format_record(&header, rec, &mut out).expect("the record formats");
        assert_eq!(String::from_utf8_lossy(&out), format!("{line}\n"));
    }
    // Too few samples is still the general parser's error.
    let mut rec = Record::new();
    let mut buf = lines[4].as_bytes().to_vec();
    let (result, _) = log::capture(|| vcf::parse(&mut header, &mut buf, &mut rec));
    assert!(result.is_err());
}
