//! A BGZF block whose data is corrupt (here, a wrong CRC32) in the middle of an input fails the
//! command with status 255 and a message, after the records before it. bcftools prints HTSlib's
//! checksum message and, in these commands, takes the failed read for the end of the input (or
//! of the side file), exiting with status 0 with the records after it lost.

use std::{
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, ExitStatus, Output},
};

/// A directory for one test's files, removed when dropped.
struct TestDir(PathBuf);

impl TestDir {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "vtools-corrupt-input-{name}-{}",
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

/// The number of records: enough for several BGZF blocks of every file.
const RECORDS: u32 = 4000;

/// A VCF of [`RECORDS`] records with GTs, PLs, and an indel every fifth record.
fn vcf_text() -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n\
         ##INFO=<ID=QS,Number=R,Type=Float,Description=\"QS\">\n\
         ##INFO=<ID=PAD,Number=1,Type=String,Description=\"Padding\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
    );
    for i in 1..=RECORDS {
        let (r, a) = if i % 5 == 0 { ("AC", "A") } else { ("A", "G") };
        let (g1, p1) = if i % 3 == 0 {
            ("0/1", "30,0,40")
        } else {
            ("0/0", "0,30,60")
        };
        s.push_str(&format!(
            "chr1\t{}\t.\t{r}\t{a}\t50\tPASS\tQS=0.6,0.4;PAD={i:060}\tGT:PL\t{g1}:{p1}\t1/1:60,30,0\n",
            i * 100
        ));
    }
    s
}

/// Writes `text` as BGZF.
fn write_bgzf(path: &Path, text: &[u8]) {
    let mut w = pasta::bgzf::io::Writer::new(fs::File::create(path).expect("create"));
    w.write_all(text).expect("write");
    w.finish().expect("finish");
}

/// The start of the middle block of the BGZF data `data` (of at least four blocks), and its size.
fn middle_block(data: &[u8], name: &Path) -> (usize, usize) {
    let mut blocks = Vec::new();
    let mut p = 0;
    while p < data.len() {
        let size = usize::from(u16::from_le_bytes([data[p + 16], data[p + 17]])) + 1;
        blocks.push((p, size));
        p += size;
    }
    assert!(
        blocks.len() >= 4,
        "{} has {} blocks",
        name.display(),
        blocks.len()
    );
    blocks[blocks.len() / 2]
}

/// Copies the BGZF file `good` to `bad` with the CRC32 of its middle block changed.
fn corrupt(good: &Path, bad: &Path) {
    let mut data = fs::read(good).expect("read");
    let (start, size) = middle_block(&data, good);
    data[start + size - 8] ^= 0xff;
    fs::write(bad, data).expect("write");
}

/// Copies the BGZF file `good` to `cut` up to its middle block, as a writer killed between two
/// blocks leaves it: every block is intact, and the EOF marker block is missing.
fn cut_at_block(good: &Path, cut: &Path) {
    let data = fs::read(good).expect("read");
    let (start, _) = middle_block(&data, good);
    fs::write(cut, &data[..start]).expect("write");
}

/// The number of records of the bgzipped VCF `path` in the blocks before its middle block:
/// the records a reader gets before the damaged block of a [`corrupt`] copy.
fn records_before_middle_block(path: &Path) -> usize {
    let data = fs::read(path).expect("read");
    let (start, _) = middle_block(&data, path);
    let mut text = Vec::new();
    pasta::bgzf::io::Reader::new(&data[..start])
        .read_to_end(&mut text)
        .expect("inflate the blocks");
    text.split_inclusive(|&c| c == b'\n')
        .filter(|l| l.ends_with(b"\n") && !l.starts_with(b"#"))
        .count()
}

/// The records of VCF text: what follows the header lines (lines starting with `#`), which
/// name the input file.
fn records(text: &[u8]) -> &[u8] {
    let mut p = 0;
    while text[p..].starts_with(b"#") {
        p += text[p..]
            .iter()
            .position(|&c| c == b'\n')
            .map_or(text.len() - p, |i| i + 1);
    }
    &text[p..]
}

/// The number of records (lines not starting with `#`) of VCF text.
fn count_records(text: &[u8]) -> usize {
    text.split(|&c| c == b'\n')
        .filter(|l| !l.is_empty() && !l.starts_with(b"#"))
        .count()
}

fn vtools(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .current_dir(dir)
        .output()
        .expect("run vtools")
}

#[track_caller]
fn assert_success(dir: &Path, args: &[&str]) {
    let out = vtools(dir, args);
    assert!(
        out.status.success(),
        "{args:?}: {}; stderr: {}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
}

/// `args` with every `{}` replaced by `f`.
fn with_file(args: &[&str], f: &str) -> Vec<String> {
    args.iter().map(|a| a.replace("{}", f)).collect()
}

/// Runs `args` with `good` and then with `bad` for every `{}` in them: the first succeeds,
/// the second fails with status 255 and a message. Returns both outputs, or what went wrong.
fn run_on_corrupt(
    dir: &Path,
    args: &[&str],
    good: &str,
    bad: &str,
) -> Result<(Output, Output), String> {
    let good_args = with_file(args, good);
    let good_args: Vec<&str> = good_args.iter().map(String::as_str).collect();
    let good_out = vtools(dir, &good_args);
    if !good_out.status.success() {
        return Err(format!(
            "{good_args:?}: {}; stderr: {}",
            good_out.status,
            String::from_utf8_lossy(&good_out.stderr)
        ));
    }
    let bad_args = with_file(args, bad);
    let bad_args: Vec<&str> = bad_args.iter().map(String::as_str).collect();
    let bad_out = vtools(dir, &bad_args);
    let err = String::from_utf8_lossy(&bad_out.stderr);
    if bad_out.status.code() != Some(255) {
        return Err(format!("{bad_args:?}: {}; stderr: {err}", bad_out.status));
    }
    if err.trim().is_empty() {
        return Err(format!("{bad_args:?} does not report the failed read"));
    }
    Ok((good_out, bad_out))
}

/// Checks every case of a test before failing, so that each case that fails is listed.
#[track_caller]
fn assert_no_failures(failures: &[String]) {
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// [`run_on_corrupt`] for every argument list of `cases`; returns what went wrong.
fn check_fails_on_corrupt(dir: &Path, cases: &[&[&str]], good: &str, bad: &str) -> Vec<String> {
    cases
        .iter()
        .filter_map(|args| run_on_corrupt(dir, args, good, bad).err())
        .collect()
}

/// [`run_on_corrupt`] for a command that writes VCF text: the failed run has written the
/// records before the damaged block in full, as bcftools writes them with status 0, rather
/// than output cut at a buffer boundary. That is a prefix of the good run's output that ends at
/// a line, with at least `min_records` records.
fn check_keeps_records_before(
    dir: &Path,
    args: &[&str],
    (good, bad): (&str, &str),
    min_records: usize,
) -> Option<String> {
    let (good_out, bad_out) = match run_on_corrupt(dir, args, good, bad) {
        Ok(outs) => outs,
        Err(e) => return Some(e),
    };
    let text = records(&bad_out.stdout);
    let n = count_records(text);
    if !records(&good_out.stdout).starts_with(text)
        || text.last() != Some(&b'\n')
        || n < min_records
    {
        let tail = String::from_utf8_lossy(&text[text.len().saturating_sub(40)..]).into_owned();
        return Some(format!(
            "{args:?} on {bad}: {n} records (at least {min_records} expected), ending {tail:?}"
        ));
    }
    None
}

/// The good and corrupt inputs: `good.vcf.gz`, `good.bcf`, `bad.vcf.gz`, `bad.bcf` (indexed),
/// a tabix-indexed `good.tsv.gz` and `bad.tsv.gz` of CHROM, FROM, TO, and a value per record,
/// `good.map.gz` and `bad.map.gz` (genetic maps), and `xv.hdr`.
fn inputs(name: &str) -> TestDir {
    let t = TestDir::new(name);
    let dir = t.path();
    fs::write(dir.join("good.vcf"), vcf_text()).expect("write the VCF");
    assert_success(dir, &["view", "-Oz", "-o", "good.vcf.gz", "good.vcf"]);
    assert_success(dir, &["view", "-Ob", "-o", "good.bcf", "good.vcf"]);
    let mut tab = String::new();
    let mut map = String::from("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n");
    for i in 1..=RECORDS {
        let pos = i * 100;
        tab.push_str(&format!("chr1\t{pos}\t{}\tv{i:060}\n", pos + 2));
        map.push_str(&format!("{pos} 1.0 {}.{i:060}\n", i / 10));
    }
    write_bgzf(&dir.join("good.tsv.gz"), tab.as_bytes());
    write_bgzf(&dir.join("good.map.gz"), map.as_bytes());
    assert_success(dir, &["tabix", "-s1", "-b2", "-e3", "good.tsv.gz"]);
    for f in ["vcf.gz", "bcf"] {
        assert_success(dir, &["index", &format!("good.{f}")]);
        corrupt(
            &dir.join(format!("good.{f}")),
            &dir.join(format!("bad.{f}")),
        );
        fs::copy(
            dir.join(format!("good.{f}.csi")),
            dir.join(format!("bad.{f}.csi")),
        )
        .expect("copy the index");
    }
    corrupt(&dir.join("good.tsv.gz"), &dir.join("bad.tsv.gz"));
    fs::copy(dir.join("good.tsv.gz.tbi"), dir.join("bad.tsv.gz.tbi")).expect("copy the index");
    corrupt(&dir.join("good.map.gz"), &dir.join("bad.map.gz"));
    fs::write(
        dir.join("xv.hdr"),
        "##INFO=<ID=XV,Number=1,Type=String,Description=\"Value\">\n",
    )
    .expect("write the header");
    t
}

#[test]
fn concat_fails_on_a_corrupt_block() {
    let t = inputs("concat");
    let dir = t.path();
    let mut failures = Vec::new();
    for f in ["vcf.gz", "bcf"] {
        let (good, bad) = (format!("good.{f}"), format!("bad.{f}"));
        let cases: &[&[&str]] = &[
            &["concat", "-Oz", "{}"],
            &["concat", "-Ob", "{}"],
            &["concat", "--threads", "2", "-Oz", "{}"],
        ];
        failures.extend(check_fails_on_corrupt(dir, cases, &good, &bad));
    }
    assert_no_failures(&failures);
}

#[test]
fn concat_writes_the_records_before_a_corrupt_block() {
    let t = inputs("concat-records");
    let dir = t.path();
    let before = records_before_middle_block(&dir.join("good.vcf.gz"));
    let mut failures = Vec::new();
    for (f, min) in [("vcf.gz", before), ("bcf", 1)] {
        let files = (&*format!("good.{f}"), &*format!("bad.{f}"));
        for args in [
            &["concat", "{}"][..],
            &["concat", "--threads", "2", "-Ov", "{}"],
        ] {
            failures.extend(check_keeps_records_before(dir, args, files, min));
        }
    }
    assert_no_failures(&failures);
}

#[test]
fn call_fails_on_a_corrupt_block() {
    let t = inputs("call");
    let dir = t.path();
    // The synced reader has read the record after the last one called.
    let before = records_before_middle_block(&dir.join("good.vcf.gz")) - 1;
    let mut failures = Vec::new();
    for (f, min) in [("vcf.gz", before), ("bcf", 1)] {
        let files = (&*format!("good.{f}"), &*format!("bad.{f}"));
        for args in [
            &["call", "-m", "{}"][..],
            &["call", "-m", "--threads", "2", "{}"],
            &["call", "-m", "-r", "chr1", "{}"],
        ] {
            failures.extend(check_keeps_records_before(dir, args, files, min));
        }
        let cases: &[&[&str]] = &[&["call", "-m", "-Oz", "{}"]];
        failures.extend(check_fails_on_corrupt(dir, cases, files.0, files.1));
    }
    assert_no_failures(&failures);
}

#[test]
fn head_fails_on_a_corrupt_block() {
    let t = inputs("head");
    let dir = t.path();
    let before = records_before_middle_block(&dir.join("good.vcf.gz"));
    let mut failures = Vec::new();
    for (f, min) in [("vcf.gz", before), ("bcf", 1)] {
        let files = (&*format!("good.{f}"), &*format!("bad.{f}"));
        let args = &["head", "-n", "100000", "{}"];
        failures.extend(check_keeps_records_before(dir, args, files, min));
    }
    assert_no_failures(&failures);
}

#[test]
fn side_files_with_a_corrupt_block_fail() {
    let t = inputs("side");
    let dir = t.path();
    let cases: &[&[&str]] = &[
        &["view", "-T", "{}", "good.vcf.gz"],
        &["view", "-T", "{}", "--threads", "2", "good.bcf"],
        &["norm", "-m-", "-T", "{}", "good.vcf.gz"],
        &["call", "-m", "-T", "{}", "good.vcf.gz"],
        &["stats", "-E", "{}", "good.vcf.gz"],
        &[
            "annotate",
            "-a",
            "{}",
            "-h",
            "xv.hdr",
            "-c",
            "CHROM,POS,-,XV",
            "good.vcf.gz",
        ],
        // Read through regidx, which failed before too.
        &[
            "annotate",
            "-a",
            "{}",
            "-h",
            "xv.hdr",
            "-c",
            "CHROM,FROM,TO,XV",
            "good.vcf.gz",
        ],
    ];
    let mut failures = check_fails_on_corrupt(dir, cases, "good.tsv.gz", "bad.tsv.gz");
    let cases: &[&[&str]] = &[&["roh", "-G30", "--AF-dflt", "0.4", "-m", "{}", "good.vcf.gz"]];
    failures.extend(check_fails_on_corrupt(
        dir,
        cases,
        "good.map.gz",
        "bad.map.gz",
    ));
    assert_no_failures(&failures);
}

/// How the input reaches a command in [`check_fails_on_missing_eof`].
#[derive(Clone, Copy, Debug)]
enum Feed {
    /// The file named by `{}`, and also on standard input for `-`.
    File,
    /// Standard input through a pipe, which cannot be checked for the EOF marker up front.
    Pipe,
}

/// The VCF text of a command's output, decompressing BGZF or BCF output with `vtools view`.
fn output_text(out: Vec<u8>) -> Vec<u8> {
    if !out.starts_with(&[0x1f, 0x8b]) {
        return out;
    }
    let mut view = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(["view", "-"])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("run vtools view");
    let mut stdin = view.stdin.take().expect("piped");
    let feeder = std::thread::spawn(move || stdin.write_all(&out));
    let text = view.wait_with_output().expect("view the output").stdout;
    feeder
        .join()
        .expect("feed")
        .expect("write the output to view");
    text
}

/// Runs `args` with `file` for `{}`, fed as `feed` says, and returns its output and stderr.
fn run_fed(dir: &Path, args: &[&str], file: &str, feed: Feed) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vtools"));
    cmd.args(with_file(args, file)).current_dir(dir);
    match feed {
        Feed::File => cmd
            .stdin(fs::File::open(dir.join(file)).expect("open the input"))
            .output()
            .expect("run vtools"),
        Feed::Pipe => {
            let mut child = cmd
                .stdin(std::process::Stdio::piped())
                .stdout(std::process::Stdio::piped())
                .stderr(std::process::Stdio::piped())
                .spawn()
                .expect("run vtools");
            let mut stdin = child.stdin.take().expect("piped");
            let data = fs::read(dir.join(file)).expect("read the input");
            // The command may stop reading early; what it did not read does not matter.
            let feeder = std::thread::spawn(move || {
                let _ = stdin.write_all(&data);
            });
            let out = child.wait_with_output().expect("wait for vtools");
            feeder.join().expect("feed");
            out
        }
    }
}

/// Runs `args` on `file`, a bgzipped VCF or BCF without the EOF marker block: it writes the
/// `records` records of the blocks it has, prints HTSlib's warning for the missing marker, and
/// then fails with status 255 and a message naming what is wrong. Returns what went wrong.
fn check_fails_on_missing_eof(
    dir: &Path,
    args: &[&str],
    (file, feed): (&str, Feed),
    records: usize,
) -> Option<String> {
    let out = run_fed(dir, args, file, feed);
    let err = String::from_utf8_lossy(&out.stderr);
    let n = count_records(&output_text(out.stdout));
    let warned = err.contains("[W::bgzf_read_block] EOF marker is absent");
    let failed = err.contains("Error: no BGZF EOF marker; the input may be truncated");
    (out.status.code() != Some(255) || !warned || !failed || n != records).then(|| {
        format!(
            "{args:?} on {file} ({feed:?}): {}, {n} of {records} records; stderr: {err}",
            out.status
        )
    })
}

/// The inputs of the missing EOF marker tests: `cut.vcf.gz` (the intact blocks before the
/// middle block of `good.vcf.gz`), `noeof.bcf` (`good.bcf` without the marker alone), and
/// tabix-indexed `cut.tsv.gz` and `noeof.tsv.gz` made the same way (with the index of the
/// intact file).
fn missing_eof_inputs(name: &str) -> TestDir {
    let t = inputs(name);
    let dir = t.path();
    cut_at_block(&dir.join("good.vcf.gz"), &dir.join("cut.vcf.gz"));
    cut_at_block(&dir.join("good.tsv.gz"), &dir.join("cut.tsv.gz"));
    for f in ["bcf", "tsv.gz"] {
        let data = fs::read(dir.join(format!("good.{f}"))).expect("read");
        fs::write(dir.join(format!("noeof.{f}")), &data[..data.len() - 28]).expect("write");
    }
    for f in ["cut", "noeof"] {
        fs::copy(
            dir.join("good.tsv.gz.tbi"),
            dir.join(format!("{f}.tsv.gz.tbi")),
        )
        .expect("copy the index");
    }
    t
}

#[test]
fn a_file_cut_at_a_block_boundary_fails_after_its_records() {
    let t = missing_eof_inputs("cut");
    let dir = t.path();
    let records = records_before_middle_block(&dir.join("good.vcf.gz"));
    let all = RECORDS as usize;
    let mut failures = Vec::new();
    let cases: &[(&[&str], Feed)] = &[
        (&["concat", "{}"], Feed::File),
        (&["concat", "-Ob", "{}"], Feed::File),
        (&["concat", "--threads", "2", "-Ov", "{}"], Feed::File),
        (&["concat", "-Oz", "{}"], Feed::File),
        (&["concat", "--naive", "{}"], Feed::File),
        (&["sort", "{}"], Feed::File),
        (&["sort", "-"], Feed::File),
        (&["sort", "-"], Feed::Pipe),
        (&["head", "-n", "100000", "{}"], Feed::File),
        (&["head", "-n", "100000", "-"], Feed::Pipe),
        (&["reheader", "-s", "names.txt", "{}"], Feed::File),
        (&["tabix", "-a", "{}"], Feed::File),
        // The synced reader cannot check a pipe up front: it reads to the end.
        (&["view", "-"], Feed::Pipe),
        (&["view", "--threads", "2", "-"], Feed::Pipe),
    ];
    fs::write(dir.join("names.txt"), "S1 R1\n").expect("write the names");
    for &(args, feed) in cases {
        failures.extend(check_fails_on_missing_eof(
            dir,
            args,
            ("cut.vcf.gz", feed),
            records,
        ));
    }
    for (args, feed) in [
        (&["concat", "{}"][..], Feed::File),
        (&["concat", "--naive", "{}"], Feed::File),
        (&["head", "-n", "100000", "{}"], Feed::File),
        (&["reheader", "-s", "names.txt", "{}"], Feed::File),
        (&["view", "-"], Feed::Pipe),
    ] {
        failures.extend(check_fails_on_missing_eof(
            dir,
            args,
            ("noeof.bcf", feed),
            all,
        ));
    }
    assert_no_failures(&failures);
}

#[test]
fn a_side_file_read_to_its_end_without_the_eof_marker_fails() {
    let t = missing_eof_inputs("cut-side");
    let dir = t.path();
    let mut failures = Vec::new();
    // Read to the end: the blocks after the cut are missing, so the targets run out before the
    // records do (and a regions list without an index is read whole).
    let cases: &[&[&str]] = &[
        &["view", "-T", "{}", "good.vcf.gz"],
        &["view", "-T", "^{}", "good.vcf.gz"],
        &["stats", "-E", "{}", "good.vcf.gz"],
        &[
            "annotate",
            "-a",
            "{}",
            "-h",
            "xv.hdr",
            "-c",
            "CHROM,FROM,TO,XV",
            "good.vcf.gz",
        ],
    ];
    for args in cases {
        let out = run_fed(dir, args, "cut.tsv.gz", Feed::File);
        let err = String::from_utf8_lossy(&out.stderr);
        if out.status.code() != Some(255)
            || !err.contains("Error: no BGZF EOF marker; the input may be truncated: cut.tsv.gz")
        {
            failures.push(format!("{args:?}: {}; stderr: {err}", out.status));
        }
    }
    assert_no_failures(&failures);
}

/// Runs `args` on `file` and checks that it succeeds without a word about the EOF marker.
fn check_quiet_success(dir: &Path, args: &[&str], file: &str) -> Option<String> {
    let out = run_fed(dir, args, file, Feed::File);
    let err = String::from_utf8_lossy(&out.stderr);
    (!out.status.success() || err.contains("EOF marker"))
        .then(|| format!("{args:?} on {file}: {}; stderr: {err}", out.status))
}

#[test]
fn a_missing_eof_marker_is_not_reported_where_reading_stops_before_the_end() {
    let t = missing_eof_inputs("cut-early");
    let dir = t.path();
    for f in ["cut.vcf.gz", "noeof.bcf"] {
        let index = if f.ends_with(".bcf") {
            "good.bcf.csi"
        } else {
            "good.vcf.gz.csi"
        };
        fs::copy(dir.join(index), dir.join(format!("{f}.csi"))).expect("copy the index");
    }
    let mut failures = Vec::new();
    for f in ["cut.vcf.gz", "noeof.bcf"] {
        failures.extend(check_quiet_success(dir, &["head", "-n", "3", "{}"], f));
        failures.extend(check_quiet_success(dir, &["head", "-h", "1", "{}"], f));
    }
    // An index query ends at its last block, as HTSlib's does.
    failures.extend(check_quiet_success(
        dir,
        &["tabix", "{}", "chr1:100-5000"],
        "cut.vcf.gz",
    ));
    failures.extend(check_quiet_success(
        dir,
        &["tabix", "{}", "chr1:100-5000"],
        "noeof.tsv.gz",
    ));
    // A tabix-indexed targets file is read by sequence, as far as the records go.
    failures.extend(check_quiet_success(
        dir,
        &["view", "-T", "noeof.tsv.gz", "{}"],
        "good.vcf.gz",
    ));
    assert_no_failures(&failures);
}

#[test]
fn index_force_indexes_a_file_without_the_eof_marker() {
    let t = missing_eof_inputs("cut-index");
    let dir = t.path();
    let out = run_fed(
        dir,
        &["index", "-f", "-o", "cut.csi", "{}"],
        "cut.vcf.gz",
        Feed::File,
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.success()
            && err.contains("EOF marker is absent")
            && !err.contains("Error")
            && dir.join("cut.csi").exists(),
        "{}; stderr: {err}",
        out.status
    );
    // Without -f, as in bcftools.
    let out = run_fed(
        dir,
        &["index", "-o", "cut2.csi", "{}"],
        "cut.vcf.gz",
        Feed::File,
    );
    assert_eq!(out.status.code(), Some(255));
}

/// Copies the BGZF file `good` to `bad` as plain gzip (one deflate stream) with 16 bytes in the
/// middle of the stream overwritten, so that it fails or inflates to garbage there.
fn spoil_as_plain_gzip(good: &Path, bad: &Path) {
    let mut text = Vec::new();
    pasta::bgzf::io::Reader::new(fs::File::open(good).expect("open"))
        .read_to_end(&mut text)
        .expect("inflate");
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&text).expect("deflate");
    let mut data = gz.finish().expect("finish");
    let mid = data.len() / 2;
    data[mid..mid + 16].fill(0xff);
    fs::write(bad, data).expect("write");
}

/// Runs vtools with standard output written to the file `out` in `dir`, so that the output is
/// buffered as a file's (in blocks of the file system's size) rather than a pipe's, and returns
/// the status, the output, and standard error.
fn vtools_to_file(dir: &Path, args: &[&str], out: &str) -> (ExitStatus, Vec<u8>, Vec<u8>) {
    let path = dir.join(out);
    let run = Command::new(env!("CARGO_BIN_EXE_vtools"))
        .args(args)
        .current_dir(dir)
        .stdout(fs::File::create(&path).expect("create the output"))
        .output()
        .expect("run vtools");
    let mut stdout = fs::read(&path).expect("read the output");
    // The command line `stats` reports, which differs.
    let command_line = b"# The command line was:";
    if let Some(at) = stdout
        .windows(command_line.len())
        .position(|w| w == command_line)
    {
        let end = stdout[at..]
            .iter()
            .position(|&c| c == b'\n')
            .map_or(stdout.len(), |i| at + i + 1);
        stdout.drain(at..end);
    }
    (run.status, stdout, run.stderr)
}

#[test]
fn threads_give_the_serial_output_on_a_corrupt_block() {
    // The records read before the damaged block are written as a serial run writes them, and
    // the output is dropped as it drops it: without flushing the file buffer, the BGZF block
    // being filled, or the EOF marker, but with the full blocks written. The same holds for
    // plain gzip input, which fails at a malformed line.
    let t = inputs("threads");
    let dir = t.path();
    spoil_as_plain_gzip(&dir.join("good.vcf.gz"), &dir.join("spoilt.vcf.gz"));
    let commands: &[&[&str]] = &[
        &["norm", "--no-version", "-m-"],
        &["norm", "--no-version", "-m-", "-Oz"],
        &["norm", "--no-version", "-m+", "-Ob"],
        &["view", "--no-version"],
        &["filter", "--no-version", "-i", "QUAL>40", "-Oz"],
        &["annotate", "--no-version", "-x", "INFO/PAD"],
        &["query", "-f", "%POS\\n"],
        &["stats"],
        &["call", "--no-version", "-m"],
        &["concat", "--no-version", "-a", "-Oz"],
        &["+fill-tags", "--no-version"],
        &["+fill-tags", "--no-version", "-Ob"],
        &["convert", "--no-version", "-Oz"],
    ];
    let mut failures = Vec::new();
    for input in ["bad.vcf.gz", "bad.bcf", "spoilt.vcf.gz"] {
        for command in commands {
            let run = |threads: &str| {
                let mut args = command.to_vec();
                args.extend(["--threads", threads, input]);
                vtools_to_file(dir, &args, &format!("out.{threads}"))
            };
            let (serial, threaded) = (run("0"), run("4"));
            if threaded != serial {
                failures.push(format!(
                    "{command:?} on {input}: {} with {} bytes with --threads 4, {} with {} \
                     bytes without; stderr {:?} and {:?}",
                    threaded.0,
                    threaded.1.len(),
                    serial.0,
                    serial.1.len(),
                    String::from_utf8_lossy(&threaded.2),
                    String::from_utf8_lossy(&serial.2),
                ));
            }
        }
    }
    assert_no_failures(&failures);
}

#[test]
fn threads_leave_the_output_of_a_late_error_empty() {
    // An expression with an undefined tag fails after the output is opened, before anything
    // is written to it. bcftools leaves the output empty, with or without threads; the
    // compression threads of a dropped output used to close it with the BGZF EOF marker.
    let t = TestDir::new("late-error");
    let dir = t.path();
    fs::write(dir.join("good.vcf"), vcf_text()).expect("write the VCF");
    assert_success(dir, &["view", "-Oz", "-o", "good.vcf.gz", "good.vcf"]);
    assert_success(dir, &["view", "-Ob", "-o", "good.bcf", "good.vcf"]);
    let commands: &[&[&str]] = &[
        &["filter", "-i", "NOSUCH>0", "good.vcf.gz"],
        &[
            "merge",
            "--force-samples",
            "-i",
            "NOSUCH:sum",
            "good.vcf.gz",
            "good.bcf",
        ],
    ];
    let mut failures = Vec::new();
    for command in commands {
        for output_type in ["-Oz", "-Ob", "-Ou", "-Ov"] {
            for threads in ["0", "4"] {
                let mut args = command.to_vec();
                args.extend([output_type, "--threads", threads]);
                let (status, out, _) = vtools_to_file(dir, &args, &format!("out.{threads}"));
                if status.code() != Some(255) || !out.is_empty() {
                    failures.push(format!("{args:?}: {status} with {} bytes", out.len()));
                }
            }
        }
    }
    assert_no_failures(&failures);
}

/// Writes `text` as one plain gzip member (not BGZF) at compression `level`.
fn gzip_member(text: &[u8], level: u32) -> Vec<u8> {
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(level));
    gz.write_all(text).expect("compress");
    gz.finish().expect("finish")
}

/// Changes the POS of the record at `pos` to `x...` in a gzip member written uncompressed
/// (level 0), whose stored blocks hold the text as it is: the line no longer parses, and the
/// member's CRC32 no longer matches.
fn spoil_stored_record(member: &mut [u8], text: &[u8], pos: u32) {
    let line = format!("\nchr1\t{pos}\t");
    let at = text
        .windows(line.len())
        .position(|w| w == line.as_bytes())
        .expect("the record is there");
    // Enough context to be found in one stored block.
    let context = &text[at..at + 40];
    let found = member
        .windows(context.len())
        .position(|w| w == context)
        .expect("the stored text is there");
    member[found + 6] = b'x';
}

/// The plain gzip inputs: `gz_bad.vcf.gz` (one uncompressed member whose middle record's POS is
/// spoilt), `gz_multi.vcf.gz` (a compressed member with the first half, then such a member with
/// the rest), and `gz_deflate.vcf.gz` (a compressed member with a byte of its deflate stream
/// changed), and the intact `gz_good.vcf.gz`.
fn gzip_inputs(name: &str) -> TestDir {
    let t = TestDir::new(name);
    let dir = t.path();
    let text = vcf_text().into_bytes();
    fs::write(dir.join("gz_good.vcf.gz"), gzip_member(&text, 6)).expect("write");
    let middle = RECORDS / 2 * 100;
    let mut bad = gzip_member(&text, 0);
    spoil_stored_record(&mut bad, &text, middle);
    fs::write(dir.join("gz_bad.vcf.gz"), bad).expect("write");
    let quarter = format!("\nchr1\t{}\t", RECORDS / 4 * 100);
    let half = text
        .windows(quarter.len())
        .position(|w| w == quarter.as_bytes())
        .expect("the record is there")
        + 1;
    let mut multi = gzip_member(&text[..half], 6);
    let mut second = gzip_member(&text[half..], 0);
    spoil_stored_record(&mut second, &text[half..], RECORDS * 3 / 4 * 100);
    multi.extend_from_slice(&second);
    fs::write(dir.join("gz_multi.vcf.gz"), multi).expect("write");
    let mut deflate = gzip_member(&text, 6);
    let mid = deflate.len() / 2;
    deflate[mid] ^= 0x55;
    fs::write(dir.join("gz_deflate.vcf.gz"), deflate).expect("write");
    t
}

#[test]
fn a_corrupt_plain_gzip_input_is_a_read_failure() {
    let t = gzip_inputs("gzip");
    let dir = t.path();
    // The command, its status for a failed read (that of bcftools' synced reader for its
    // `BCF read error` in view), and what the message says.
    let cases: &[(&[&str], i32, &str)] = &[
        (&["head", "-n", "100000", "{}"], 255, "Error reading"),
        (&["call", "-m", "{}"], 255, "BCF read error"),
        (
            &["call", "-m", "--threads", "2", "{}"],
            255,
            "BCF read error",
        ),
        (&["view", "{}"], 9, "BCF read error"),
        (&["view", "--threads", "2", "{}"], 9, "BCF read error"),
        (&["query", "-f", "%POS\\n", "{}"], 255, "BCF read error"),
        (&["stats", "{}"], 255, "BCF read error"),
        (&["stats", "--threads", "2", "{}"], 255, "BCF read error"),
        (
            &["norm", "-m-", "--threads", "2", "{}"],
            255,
            "BCF read error",
        ),
        (&["norm", "-m-", "{}"], 255, "BCF read error"),
    ];
    let mut failures = Vec::new();
    for &(args, status, message) in cases {
        let a = with_file(args, "gz_good.vcf.gz");
        let good = vtools(dir, &a.iter().map(String::as_str).collect::<Vec<_>>());
        if !good.status.success() {
            failures.push(format!("{args:?} on gz_good.vcf.gz: {}", good.status));
        }
        for f in ["gz_bad.vcf.gz", "gz_multi.vcf.gz"] {
            let a = with_file(args, f);
            let out = vtools(dir, &a.iter().map(String::as_str).collect::<Vec<_>>());
            let err = String::from_utf8_lossy(&out.stderr);
            if out.status.code() != Some(status) || !err.contains(message) {
                failures.push(format!("{args:?} on {f}: {}; stderr: {err}", out.status));
            }
        }
        // What a corrupt deflate stream inflates to before the error shows may parse, and
        // fail the command on its content instead: it fails either way.
        let a = with_file(args, "gz_deflate.vcf.gz");
        let out = vtools(dir, &a.iter().map(String::as_str).collect::<Vec<_>>());
        if out.status.success() || out.stderr.is_empty() {
            failures.push(format!("{args:?} on gz_deflate.vcf.gz: {}", out.status));
        }
    }
    // concat and sort fail on the spoilt line (as in bcftools) or the read; the output of
    // head is the records before the spoilt line, in full.
    for f in ["gz_bad.vcf.gz", "gz_multi.vcf.gz", "gz_deflate.vcf.gz"] {
        for args in [&["concat", "{}"][..], &["sort", "{}"]] {
            let a = with_file(args, f);
            let out = vtools(dir, &a.iter().map(String::as_str).collect::<Vec<_>>());
            if out.status.code() != Some(255) {
                failures.push(format!("{args:?} on {f}: {}", out.status));
            }
        }
    }
    let out = vtools(dir, &["head", "-n", "100000", "gz_bad.vcf.gz"]);
    let n = count_records(&out.stdout);
    if n != RECORDS as usize / 2 - 1 {
        failures.push(format!("head on gz_bad.vcf.gz wrote {n} records"));
    }
    assert_no_failures(&failures);
}
