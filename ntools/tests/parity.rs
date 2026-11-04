//! Small differences from samtools 1.24 found in the second pre-launch bug hunt: `sort -n` and
//! `merge -n` on names that differ only in zero padding, BGZF SAM output, `depth` past
//! `@SQ LN`, `collate`'s order within a template, `faidx` regions starting at 0, usage without
//! arguments, and `dict` reading stdin.

use std::{
    fs,
    io::{self, Write},
    path::Path,
    process::{Command, Output, Stdio},
};

const NTOOLS: &str = env!("CARGO_BIN_EXE_ntools");

fn run(dir: &Path, args: &[&str]) -> io::Result<Output> {
    Command::new(NTOOLS)
        .current_dir(dir)
        .args(args)
        .stdin(Stdio::null())
        .output()
}

/// Runs ntools in `dir`, failing the test if it fails, and returns its standard output.
fn ntools(dir: &Path, args: &[&str]) -> io::Result<String> {
    let output = run(dir, args)?;
    assert!(
        output.status.success(),
        "ntools {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).map_err(io::Error::other)
}

/// The records of SAM text, without the header.
fn records(sam: &str) -> Vec<&str> {
    sam.lines().filter(|line| !line.starts_with('@')).collect()
}

/// The first `n` fields of each record of SAM text, joined with `/`.
fn fields(sam: &str, n: usize) -> Vec<String> {
    records(sam)
        .iter()
        .map(|line| line.split('\t').take(n).collect::<Vec<_>>().join("/"))
        .collect()
}

/// An unmapped record named `name` with flags `flags`.
fn unmapped(name: &str, flags: u16) -> String {
    format!("{name}\t{flags}\t*\t0\t0\t*\t*\t0\t0\tA\tI\n")
}

/// Samtools' `strnum_cmp` compares names that differ only in zero padding as equal, so
/// they sort by the READ1/READ2 flags and then in input order.
#[test]
fn test_sort_and_merge_n_compare_zero_padded_names_as_equal() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();

    let mut sam = String::from("@SQ\tSN:c1\tLN:1000\n");
    for name in [
        "00", "0", "r007", "r7", "r07", "r7", "x01y2", "x1y02", "x1y2",
    ] {
        sam.push_str(&unmapped(name, 4));
    }
    sam.push_str(&unmapped("q01", 133));
    sam.push_str(&unmapped("q1", 69));
    sam.push_str(&unmapped("q01", 69));
    fs::write(dir.join("in.sam"), sam)?;

    // As samtools 1.24 sorts them.
    let sorted = ntools(dir, &["sort", "-n", "-O", "sam", "in.sam"])?;
    assert_eq!(
        fields(&sorted, 2),
        [
            "00/4", "0/4", "q1/69", "q01/69", "q01/133", "r007/4", "r7/4", "r07/4", "r7/4",
            "x01y2/4", "x1y02/4", "x1y2/4",
        ]
    );

    let header = "@HD\tVN:1.6\tSO:queryname\n";
    let a: String = ["r007", "r7", "x1"].map(|n| unmapped(n, 4)).concat();
    let b: String = ["r07", "r0007", "x01", "x001"]
        .map(|n| unmapped(n, 4))
        .concat();
    fs::write(dir.join("a.sam"), format!("{header}{a}"))?;
    fs::write(dir.join("b.sam"), format!("{header}{b}"))?;

    // Equal names come from the files in the order they are given, as in samtools.
    for (inputs, expected) in [
        (
            ["a.sam", "b.sam"],
            ["r007", "r7", "r07", "r0007", "x1", "x01", "x001"],
        ),
        (
            ["b.sam", "a.sam"],
            ["r07", "r0007", "r007", "r7", "x01", "x001", "x1"],
        ),
    ] {
        let mut args = vec!["merge", "-n", "-O", "sam", "-o", "-"];
        args.extend(inputs);
        let merged = ntools(dir, &args)?;
        assert_eq!(fields(&merged, 1), expected, "{inputs:?}");
    }

    Ok(())
}

/// Whether `data` starts with a BGZF block header, and its text decompressed.
fn bgzf_text(data: &[u8]) -> io::Result<(bool, String)> {
    use io::Read;

    let is_bgzf = data.len() >= 16 && data[..4] == [0x1f, 0x8b, 8, 4] && data[12..14] == *b"BC";
    let mut text = Vec::new();
    flate2::read::MultiGzDecoder::new(data).read_to_end(&mut text)?;
    Ok((is_bgzf, String::from_utf8_lossy(&text).into_owned()))
}

/// `-O sam.gz`, and an output name ending in `.sam.gz` (also with `-O sam`), give BGZF
/// SAM, as in samtools; `collate -o FILE` takes its format from the extension, else SAM.
#[test]
fn test_sam_gz_output_is_bgzf_sam() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();

    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n");
    for i in 1..=3 {
        sam.push_str(&format!(
            "r{i}\t0\tc1\t{}\t30\t4M\t*\t0\t0\tACGT\tIIII\n",
            i * 10
        ));
    }
    fs::write(dir.join("in.sam"), &sam)?;
    ntools(dir, &["view", "-b", "-o", "in.bam", "in.sam"])?;

    // The command, the output it writes, and whether that has a header.
    let cases: &[(&[&str], &str, bool)] = &[
        (&["view", "-o", "a.sam.gz"], "a.sam.gz", false),
        (&["view", "-h", "-o", "b.sam.gz"], "b.sam.gz", true),
        (&["view", "-O", "sam.gz", "-o", "c.out"], "c.out", false),
        (&["view", "-O", "SAM.GZ", "-o", "d.out"], "d.out", false),
        (&["view", "-O", "sam", "-o", "e.SAM.gz"], "e.SAM.gz", false),
        (&["view", "-O", "sam.gz", "-o", "f.bam"], "f.bam", true),
        (
            &["view", "-@", "2", "-h", "-o", "g.sam.gz"],
            "g.sam.gz",
            true,
        ),
        (&["sort", "-o", "h.sam.gz"], "h.sam.gz", true),
        (
            &["sort", "-O", "sam", "-@", "2", "-o", "i.sam.gz"],
            "i.sam.gz",
            true,
        ),
        (&["merge", "-o", "j.sam.gz"], "j.sam.gz", true),
        (&["collate", "-o", "k.sam.gz"], "k.sam.gz", true),
        // samtools names it `prefix.sam`, hiding its compression (README).
        (
            &["collate", "--output-fmt", "sam.gz"],
            "prefix.sam.gz",
            true,
        ),
    ];

    for &(args, output, with_header) in cases {
        let mut args = args.to_vec();
        args.extend(["--no-PG", "in.bam"]);
        if args[0] == "collate" && !args.contains(&"-o") {
            args.push("prefix");
        }
        ntools(dir, &args)?;

        let (is_bgzf, text) = bgzf_text(&fs::read(dir.join(output))?)?;
        assert!(is_bgzf, "{args:?} writes BGZF");
        assert_eq!(text.starts_with("@HD"), with_header, "{args:?}: {text}");
        let mut names: Vec<_> = fields(&text, 1);
        names.sort();
        assert_eq!(names, ["r1", "r2", "r3"], "{args:?}");
    }

    // Names HTSlib does not take as BGZF SAM, and `collate -o` of other names, write plain SAM.
    for (args, output) in [
        (&["view", "-o", "l.sam.GZ"][..], "l.sam.GZ"),
        (&["view", "-o", "m.sam.bgz"], "m.sam.bgz"),
        (&["collate", "-o", "n.sam"], "n.sam"),
        (&["collate", "-o", "o.out"], "o.out"),
    ] {
        let mut args = args.to_vec();
        args.push("in.bam");
        ntools(dir, &args)?;
        let data = fs::read(dir.join(output))?;
        assert!(
            data.starts_with(b"@HD") || data.starts_with(b"r1"),
            "{args:?} writes plain SAM"
        );
    }

    // `out.BAM` is BAM, as in HTSlib.
    ntools(dir, &["view", "-o", "p.BAM", "in.bam"])?;
    let (is_bgzf, text) = bgzf_text(&fs::read(dir.join("p.BAM"))?)?;
    assert!(is_bgzf && text.starts_with("BAM\x01"), "p.BAM is BAM");

    Ok(())
}

/// `depth -r`, `-b`, and `-a` report the positions of reads that overhang `@SQ LN`, as
/// samtools does, with `-a` filling zeros between them but not after the last one.
#[test]
fn test_depth_reports_positions_past_the_reference_end() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();

    // r1 covers 98-102 and r2 99-100 and 105-107 (skipping 101-104), past LN 100; r3 starts past
    // LN, at 110.
    fs::write(
        dir.join("in.sam"),
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100\n@SQ\tSN:d\tLN:50\n\
         r1\t0\tc\t98\t30\t5M\t*\t0\t0\tAAAAA\t*\n\
         r2\t0\tc\t99\t30\t2M4N3M\t*\t0\t0\tAAAAA\t*\n\
         r3\t0\tc\t110\t30\t3M\t*\t0\t0\tAAA\t*\n\
         r4\t0\td\t5\t30\t3M\t*\t0\t0\tAAA\t*\n",
    )?;
    ntools(dir, &["view", "-b", "-o", "in.bam", "in.sam"])?;
    ntools(dir, &["index", "in.bam"])?;
    fs::write(dir.join("a.bed"), "c\t95\t200\n")?;
    fs::write(dir.join("b.bed"), "c\t107\t108\nc\t111\t130\n")?;

    let covered = "98:1 99:2 100:2 101:1 102:1 103:0 104:0 105:1 106:1 107:1 110:1 111:1 112:1";
    let filled = "98:1 99:2 100:2 101:1 102:1 103:0 104:0 105:1 106:1 107:1 108:0 109:0 \
                  110:1 111:1 112:1";

    // The positions samtools 1.24 prints, as `POS:DEPTH` on `c`.
    let cases: &[(&[&str], String)] = &[
        (&["-r", "c"], covered.into()),
        (&["-b", "a.bed"], covered.into()),
        (&["-a", "-b", "a.bed"], format!("96:0 97:0 {filled}")),
        (
            &["-r", "c:95-105"],
            "98:1 99:2 100:2 101:1 102:1 103:0 104:0 105:1".into(),
        ),
        (
            &["-a", "-r", "c:101-120"],
            "101:1 102:1 103:0 104:0 105:1 106:1 107:1 108:0 109:0 110:1 111:1 112:1".into(),
        ),
        (
            &["-a", "-r", "c:108-120"],
            "108:0 109:0 110:1 111:1 112:1".into(),
        ),
        // r3 is not read for this region, so nothing is past 107.
        (
            &["-a", "-r", "c:90-108"],
            "90:0 91:0 92:0 93:0 94:0 95:0 96:0 97:0 \
             98:1 99:2 100:2 101:1 102:1 103:0 104:0 105:1 106:1 107:1"
                .into(),
        ),
        (&["-a", "-r", "c:150-160"], String::new()),
        // 108 is filled, as a depth (r3's) follows it; nothing follows 130.
        (&["-a", "-b", "b.bed"], "108:0 112:1".into()),
        (&["-a", "-r", "c"], {
            let zeros: Vec<_> = (1..98).map(|p| format!("{p}:0")).collect();
            format!("{} {filled}", zeros.join(" "))
        }),
    ];

    for (args, expected) in cases {
        let mut command = vec!["depth"];
        command.extend_from_slice(args);
        command.push("in.bam");
        let output = ntools(dir, &command)?;
        let positions: Vec<_> = output
            .lines()
            .filter_map(|line| line.strip_prefix("c\t"))
            .map(|rest| rest.replace('\t', ":"))
            .collect();
        assert_eq!(positions.join(" "), *expected, "depth {args:?}");
    }

    Ok(())
}

/// `collate` writes READ1 before READ2 within a template, as samtools does. `collate -f`
/// keeps every primary record; samtools also drops those without exactly one of READ1 and READ2,
/// such as single-end reads, which its usage does not say (README).
#[test]
fn test_collate_orders_read1_first_and_fast_mode_keeps_primaries() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();

    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n");
    for (name, flags, pos) in [
        ("p1", 163, 10), // READ2 first in coordinate order
        ("se1", 0, 20),  // single-end
        ("both", 193, 27),
        ("p2", 131, 28),
        ("p2", 323, 29),  // READ1 secondary
        ("p2", 2115, 30), // READ1 supplementary
        ("p1", 83, 50),
        ("p2", 67, 60),
        ("r1only", 64, 61),
    ] {
        sam.push_str(&format!(
            "{name}\t{flags}\tc1\t{pos}\t30\t4M\t*\t0\t0\tACGT\tIIII\n"
        ));
    }
    fs::write(dir.join("in.sam"), sam)?;

    // Each template's records, in output order, as `NAME/FLAGS`, by template name.
    let templates = |args: &[&str]| -> io::Result<Vec<String>> {
        let mut command = vec!["collate", "-O", "--output-fmt", "sam"];
        command.extend_from_slice(args);
        command.extend(["in.sam", "prefix"]);
        let output = ntools(dir, &command)?;
        let mut groups: Vec<String> = Vec::new();
        let mut last = "";
        for line in records(&output) {
            let name = line.split('\t').next().unwrap_or_default();
            let record = fields(line, 2).join("");
            if name == last {
                if let Some(group) = groups.last_mut() {
                    group.push(' ');
                    group.push_str(&record);
                }
            } else {
                groups.push(record);
            }
            last = name;
        }
        groups.sort();
        Ok(groups)
    };

    assert_eq!(
        templates(&[])?,
        [
            "both/193",
            "p1/83 p1/163",
            "p2/323 p2/2115 p2/67 p2/131",
            "r1only/64",
            "se1/0",
        ]
    );
    assert_eq!(
        templates(&["-f"])?,
        [
            "both/193",
            "p1/83 p1/163",
            "p2/67 p2/131",
            "r1only/64",
            "se1/0"
        ]
    );

    Ok(())
}

/// A `faidx` region starting at 0 followed by an end (`a:0-4`, `a:0-0`, BED-style
/// coordinates) prints an empty sequence with samtools' warning; `a:0` is the whole sequence.
#[test]
fn test_faidx_region_starting_at_zero_is_empty() -> io::Result<()> {
    let dir = tempfile::tempdir()?;
    let dir = dir.path();
    fs::write(dir.join("x.fa"), ">a\nACGTACGTAC\nGGGGG\n")?;
    fs::write(dir.join("x.fq"), "@a\nACGTACGTAC\n+\nIIIIIIIIII\n")?;
    ntools(dir, &["faidx", "x.fa"])?;
    ntools(dir, &["fqidx", "x.fq"])?;

    for region in ["a:0-4", "a:0-0", "a:0-", "a:0-20", "{a}:0-4"] {
        let output = run(dir, &["faidx", "x.fa", region])?;
        assert!(output.status.success(), "{region}");
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            format!(">{region}\n")
        );
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!("[faidx] Zero length sequence: {region}\n"),
        );
    }

    // samtools warns for the sequence and for the qualities.
    let output = run(dir, &["fqidx", "x.fq", "a:0-4"])?;
    assert!(output.status.success());
    assert_eq!(String::from_utf8_lossy(&output.stdout), "@a:0-4\n+\n");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "[faidx] Zero length sequence: a:0-4\n".repeat(2)
    );

    assert_eq!(
        ntools(dir, &["faidx", "x.fa", "a:0", "a:-4", "a:1-4"])?,
        ">a:0\nACGTACGTAC\nGGGGG\n>a:-4\nACGT\n>a:1-4\nACGT\n"
    );

    Ok(())
}

/// Run without arguments, these commands print their usage to standard output and exit 0, as in
/// samtools.
#[test]
fn test_usage_without_arguments_exits_0() -> io::Result<()> {
    let dir = tempfile::tempdir()?;

    for command in [
        "coverage", "depth", "faidx", "fixmate", "flagstat", "fqidx", "idxstats", "import",
        "index", "merge",
    ] {
        let output = run(dir.path(), &[command])?;
        assert!(output.status.success(), "{command}: {}", output.status);
        let usage = String::from_utf8_lossy(&output.stdout);
        assert!(
            usage.contains(&format!("Usage: ntools {command}")),
            "{command}: {usage}"
        );
        assert!(output.stderr.is_empty(), "{command} writes no error");
    }

    // Other arguments, or none for a command that reads stdin, are not this case.
    let output = run(dir.path(), &["flagstat", "nosuch.bam"])?;
    assert_eq!(output.status.code(), Some(1));

    Ok(())
}

/// `dict` reads stdin when no input is given, as samtools does (`samtools dict < x.fa`).
#[test]
fn test_dict_reads_stdin_without_input() -> io::Result<()> {
    let dir = tempfile::tempdir()?;

    let mut child = Command::new(NTOOLS)
        .current_dir(dir.path())
        .arg("dict")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    child
        .stdin
        .take()
        .expect("stdin is piped")
        .write_all(b">a\nACGT\n")?;
    let output = child.wait_with_output()?;

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:a\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\n"
    );

    Ok(())
}

/// Opens a pseudo-terminal and returns its controller and terminal ends.
fn open_terminal() -> io::Result<(std::os::fd::OwnedFd, std::os::fd::OwnedFd)> {
    use std::os::fd::FromRawFd;

    let (mut controller, mut terminal) = (-1, -1);

    // SAFETY: the out-pointers are valid; the name, termios, and winsize arguments may be null.
    let result = unsafe {
        libc::openpty(
            &mut controller,
            &mut terminal,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        )
    };

    if result != 0 {
        return Err(io::Error::last_os_error());
    }

    // SAFETY: `openpty` succeeded, so both are open file descriptors that nothing else owns.
    unsafe {
        Ok((
            std::os::fd::OwnedFd::from_raw_fd(controller),
            std::os::fd::OwnedFd::from_raw_fd(terminal),
        ))
    }
}

/// With no input and stdin a terminal, these commands print their usage, on samtools' stream
/// and with its exit status, rather than waiting on the terminal.
#[test]
fn test_usage_without_input_on_a_terminal() -> io::Result<()> {
    let dir = tempfile::tempdir()?;

    for (args, stdout, status) in [
        (&["sort", "-n"][..], true, 0),
        (&["bam2fq"], true, 0),
        (&["reset"], true, 0),
        (&["view"], true, 0),
        (&["view", "-c"], false, 1),
        (&["stats"], true, 1),
        (&["samples"], false, 1),
    ] {
        let (_controller, terminal) = open_terminal()?;
        let output = Command::new(NTOOLS)
            .current_dir(dir.path())
            .args(args)
            .stdin(terminal)
            .output()?;

        assert_eq!(output.status.code(), Some(status), "{args:?}");
        let (usage, other) = if stdout {
            (&output.stdout, &output.stderr)
        } else {
            (&output.stderr, &output.stdout)
        };
        let usage = String::from_utf8_lossy(usage);
        assert!(usage.contains("Usage: ntools "), "{args:?}: {usage}");
        assert!(other.is_empty(), "{args:?}");
    }

    // Run without arguments, mpileup prints its usage whatever stdin is.
    let output = run(dir.path(), &["mpileup"])?;
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("Usage: ntools mpileup"));

    Ok(())
}
