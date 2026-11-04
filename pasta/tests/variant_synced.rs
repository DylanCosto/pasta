//! The multi-file synced reader's pairing against `bcftools isec -n+1` (bcftools 1.24) on two
//! files with several records per position.

use std::path::PathBuf;

use pasta::variant::{
    Unpack,
    readlist::ListSource,
    synced::{IndexRequirement, Pairing, SyncedError, SyncedReaders},
};

fn data(name: &str) -> String {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/data");
    p.push(name);
    p.to_string_lossy().into_owned()
}

/// Reads both files in step, printing bcftools isec's sites list: the first record's CHROM,
/// POS, REF, ALT, and which files have a record.
fn sites(pairing: Pairing) -> String {
    let mut sr = SyncedReaders::new();
    sr.set_pairing(pairing);
    sr.set_index_requirement(IndexRequirement::AllowMissing);
    sr.add_reader(&data("synced_a.vcf"))
        .expect("the first file opens");
    sr.add_reader(&data("synced_b.vcf"))
        .expect("the second file opens");
    let mut out = String::new();
    while sr.next_line().expect("the files read") > 0 {
        let first = (0..2).find(|&i| sr.has_line(i)).expect("a file has a line");
        let bits: String = (0..2)
            .map(|i| if sr.has_line(i) { '1' } else { '0' })
            .collect();
        let (rec, header) = sr.line_with_header(first).expect("the line is there");
        rec.unpack(Unpack::STR);
        let chrom = String::from_utf8_lossy(header.contig_name(rec.rid()).expect("a contig"));
        let alleles: Vec<String> = rec
            .alleles()
            .map(|a| String::from_utf8_lossy(a).into_owned())
            .collect();
        let alt = if alleles.len() > 1 {
            alleles[1..].join(",")
        } else {
            ".".to_owned()
        };
        out.push_str(&format!(
            "{chrom}\t{}\t{}\t{alt}\t{bits}\n",
            rec.pos() + 1,
            alleles[0]
        ));
    }
    out
}

#[test]
fn exact_pairing_matches_bcftools() {
    let want = "\
chr1\t100\tA\tG\t10
chr1\t100\tA\tT\t10
chr1\t100\tAC\tA\t10
chr1\t100\tA\tC\t01
chr1\t100\tA\tAT\t01
chr1\t200\tA\tG,T\t10
chr1\t200\tA\tT\t01
chr1\t300\tA\tAT\t10
chr1\t300\tA\tATT\t01
chr1\t400\tA\t.\t10
chr1\t400\tA\tG\t01
chr1\t500\tC\tG\t11
chr1\t600\tG\t<DEL>\t11
chr1\t600\tG\t<DEL>\t10
chr1\t700\tT\tC\t11
chr1\t700\tT\tC\t10
chr1\t800\tTA\tT\t10
chr1\t800\tT\tC\t01
chr1\t900\tGT\tGA\t10
chr1\t900\tG\tA\t01
chr1\t1000\tA\t*,G\t10
chr1\t1000\tA\tG\t01
chr1\t1100\tACG\tA,AC\t11
chr1\t1200\tC\tT\t11
chr1\t1300\tG\tA\t11
chr2\t50\tC\tT\t01
";
    assert_eq!(sites(Pairing::empty()), want);
}

#[test]
fn snp_and_indel_pairing_matches_bcftools() {
    let want = "\
chr1\t100\tA\tG\t11
chr1\t100\tA\tT\t10
chr1\t100\tAC\tA\t11
chr1\t200\tA\tG,T\t11
chr1\t300\tA\tAT\t11
chr1\t400\tA\t.\t10
chr1\t400\tA\tG\t01
chr1\t500\tC\tG\t11
chr1\t600\tG\t<DEL>\t11
chr1\t600\tG\t<DEL>\t10
chr1\t700\tT\tC\t11
chr1\t700\tT\tC\t10
chr1\t800\tTA\tT\t10
chr1\t800\tT\tC\t01
chr1\t900\tGT\tGA\t11
chr1\t1000\tA\t*,G\t11
chr1\t1100\tACG\tA,AC\t11
chr1\t1200\tC\tT\t11
chr1\t1300\tG\tA\t11
chr2\t50\tC\tT\t01
";
    assert_eq!(sites(Pairing::BOTH), want);
}

#[test]
fn any_pairing_matches_bcftools() {
    let want = "\
chr1\t100\tA\tG\t11
chr1\t100\tA\tT\t11
chr1\t100\tAC\tA\t10
chr1\t200\tA\tG,T\t11
chr1\t300\tA\tAT\t11
chr1\t400\tA\t.\t11
chr1\t500\tC\tG\t11
chr1\t600\tG\t<DEL>\t11
chr1\t600\tG\t<DEL>\t10
chr1\t700\tT\tC\t11
chr1\t700\tT\tC\t10
chr1\t800\tTA\tT\t11
chr1\t900\tGT\tGA\t11
chr1\t1000\tA\t*,G\t11
chr1\t1100\tACG\tA,AC\t11
chr1\t1200\tC\tT\t11
chr1\t1300\tG\tA\t11
chr2\t50\tC\tT\t01
";
    assert_eq!(sites(Pairing::ANY), want);
}

#[test]
fn streaming_several_files_needs_permission() {
    let mut sr = SyncedReaders::new();
    sr.add_reader(&data("synced_a.vcf"))
        .expect("one file streams");
    assert!(matches!(
        sr.add_reader(&data("synced_b.vcf")),
        Err(SyncedError::ApiUsage)
    ));
    let mut sr = SyncedReaders::new();
    sr.set_index_requirement(IndexRequirement::Required);
    assert!(matches!(
        sr.add_reader(&data("synced_a.vcf")),
        Err(SyncedError::NotBgzf)
    ));
}

#[test]
fn sample_lists_follow_htslib() {
    let mut sr = SyncedReaders::new();
    sr.add_reader(&data("synced_a.vcf"))
        .expect("one file streams");
    assert_eq!(sr.set_samples(b"-", ListSource::Inline).ok(), Some(true));
    assert_eq!(sr.samples(), [b"A".to_vec()]);
    assert_eq!(sr.set_samples(b"^A", ListSource::Inline).ok(), Some(false));
    assert_eq!(sr.set_samples(b"B,A", ListSource::Inline).ok(), Some(true));
    assert_eq!(sr.sample_map(0), [0]);
    // An empty file names no samples, unlike `-`.
    let empty = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("synced_no_samples.txt");
    std::fs::write(&empty, b"").expect("the temporary file is written");
    let path = empty.to_string_lossy().into_owned();
    assert_eq!(
        sr.set_samples(path.as_bytes(), ListSource::File).ok(),
        Some(false)
    );
}
