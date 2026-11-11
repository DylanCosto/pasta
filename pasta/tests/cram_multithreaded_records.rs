//! Records decoded on several threads match the records decoded on one.
//!
//! The fixtures in `data/cram_multislice` were written by samtools 1.24 from `multislice.sam`
//! (mapped records on `sq0` and `sq1` of `ref.fa`, then unmapped records) with 7 records per
//! slice and 4 slices per container:
//!
//! - `multislice-3.0.cram`: CRAM 3.0, reference-based
//!   (`-O cram,version=3.0,seqs_per_slice=7,slices_per_container=4`).
//! - `multislice-3.1-no-ref.cram`: CRAM 3.1, reference-free (`...,version=3.1,...,no_ref=1`).
//! - `multislice-3.1-embed-ref.cram`: CRAM 3.1, with embedded references
//!   (`...,version=3.1,...,embed_ref=1`).

use std::{error::Error, fs, io, num::NonZero, path::PathBuf};

use pasta::core::Position;
use pasta::cram::{
    self as cram,
    file_definition::Version,
    io::writer::compression::{Options, Profile},
};
use pasta::fasta;
use pasta::sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::data::field::Value,
    },
};

const WORKER_COUNTS: [usize; 4] = [1, 2, 3, 8];

/// The length of `sq0`, long enough for a CIGAR of more than 65,535 operations.
const SQ0_LEN: usize = 50_000;

#[test]
fn test_written_records_match_single_threaded() -> Result<(), Box<dyn Error>> {
    let reference = reference_sequence(SQ0_LEN);
    let header = header();
    let records = records(&reference);

    for (version, profile) in [
        (Version::new(3, 0), Profile::Normal),
        (Version::new(3, 1), Profile::Normal),
        (Version::new(3, 1), Profile::Archive),
    ] {
        let src = write(&header, &reference, &records, version, profile)?;
        let repository = repository(&reference);

        let expected = read(&src, &repository, None)?;
        assert_eq!(expected.len(), records.len());

        for workers in WORKER_COUNTS {
            let actual = read(&src, &repository, Some(workers))?;
            assert_eq!(
                actual, expected,
                "{version:?} {profile:?} {workers} workers"
            );
        }
    }

    Ok(())
}

#[test]
fn test_samtools_records_match_single_threaded() -> Result<(), Box<dyn Error>> {
    let repository = fasta_repository(&data_path("cram_multislice/ref.fa"))?;

    for name in [
        "cram_multislice/multislice-3.0.cram",
        "cram_multislice/multislice-3.1-no-ref.cram",
        "cram_multislice/multislice-3.1-embed-ref.cram",
        "cram_codecs/multi-subexp-qsbeta.cram",
    ] {
        let src = fs::read(data_path(name))?;
        let expected = read(&src, &repository, None)?;
        assert!(!expected.is_empty(), "{name}");

        for workers in WORKER_COUNTS {
            let actual = read(&src, &repository, Some(workers))?;
            assert_eq!(actual, expected, "{name} {workers} workers");
        }
    }

    // The fixtures have containers of several slices.
    let src = fs::read(data_path("cram_multislice/multislice-3.0.cram"))?;
    let mut reader = cram::io::Reader::new(&src[..]);
    reader.read_header()?;
    let mut container = cram::io::reader::Container::default();
    let mut most_slices = 0;

    while reader.read_container(&mut container)? != 0 {
        most_slices = most_slices.max(container.header().landmarks().len());
    }

    assert_eq!(most_slices, 4);

    Ok(())
}

#[test]
fn test_truncated_input_ends_with_the_error_after_the_same_records() -> Result<(), Box<dyn Error>> {
    let reference = reference_sequence(SQ0_LEN);
    let header = header();
    let records = records(&reference);
    let src = write(
        &header,
        &reference,
        &records,
        Version::new(3, 1),
        Profile::Normal,
    )?;
    let repository = repository(&reference);

    for len in (0..src.len()).step_by(src.len() / 41 + 1) {
        let truncated = &src[..len];

        // The file definition or the header is cut.
        if cram_reader(truncated, &repository).read_header().is_err() {
            continue;
        }

        let expected = read_until_error(truncated, &repository, None)?;

        for workers in WORKER_COUNTS {
            let actual = read_until_error(truncated, &repository, Some(workers))?;
            assert_eq!(actual.0, expected.0, "{len} bytes, {workers} workers");
            assert_eq!(
                actual.1.as_ref().map(io::Error::kind),
                expected.1.as_ref().map(io::Error::kind),
                "{len} bytes, {workers} workers"
            );
        }
    }

    // Cut in the middle of the last container of records.
    let truncated = &src[..src.len() - 200];
    let mut reader = cram_reader(truncated, &repository);
    let header = reader.read_header()?;
    let workers = NonZero::new(4).unwrap();
    let mut results = reader.multithreaded_records(&header, workers);
    let mut n = 0;

    let e = loop {
        match results.next() {
            Some(Ok(_)) => n += 1,
            Some(Err(e)) => break e,
            None => panic!("no error after {n} records"),
        }
    };

    assert_eq!(e.kind(), io::ErrorKind::UnexpectedEof);
    assert!(n > 0 && n < records.len(), "{n} records");
    assert!(results.next().is_none(), "the records end after an error");

    Ok(())
}

#[test]
fn test_corrupt_slice_is_an_error_after_the_slices_before_it() -> Result<(), Box<dyn Error>> {
    let reference = reference_sequence(SQ0_LEN);
    let header = header();
    let records = records(&reference);
    let mut src = write(
        &header,
        &reference,
        &records,
        Version::new(3, 0),
        Profile::Normal,
    )?;
    let repository = repository(&reference);

    // A byte in the middle of the records, which fails its block's checksum.
    let i = src.len() / 2;
    src[i] ^= 0xff;

    let expected = read_until_error(&src, &repository, None)?;
    assert!(expected.1.is_some(), "the corruption is detected");
    assert!(!expected.0.is_empty());

    for workers in WORKER_COUNTS {
        let actual = read_until_error(&src, &repository, Some(workers))?;
        assert_eq!(actual.0, expected.0, "{workers} workers");
        assert_eq!(
            actual.1.as_ref().map(io::Error::kind),
            expected.1.as_ref().map(io::Error::kind),
            "{workers} workers"
        );
    }

    Ok(())
}

#[test]
fn test_dropping_the_records_early_leaves_the_reader_usable() -> Result<(), Box<dyn Error>> {
    let reference = reference_sequence(SQ0_LEN);
    let header = header();
    let records = records(&reference);
    let src = write(
        &header,
        &reference,
        &records,
        Version::new(3, 1),
        Profile::Normal,
    )?;
    let repository = repository(&reference);

    let mut reader = cram_reader(&src[..], &repository);
    let header = reader.read_header()?;
    let workers = NonZero::new(4).unwrap();

    let first: Vec<_> = reader
        .multithreaded_records(&header, workers)
        .take(3)
        .collect::<io::Result<_>>()?;
    assert_eq!(first.len(), 3);

    // The records read ahead are gone, and the rest are read from where reading stopped.
    let rest: Vec<_> = reader
        .multithreaded_records(&header, workers)
        .collect::<io::Result<_>>()?;
    assert!(rest.len() < records.len() - 3);

    Ok(())
}

fn data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn cram_reader<'a>(src: &'a [u8], repository: &fasta::Repository) -> cram::io::Reader<&'a [u8]> {
    cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_reader(src)
}

/// Reads every record, on one thread with [`cram::io::Reader::records`] without `workers`.
fn read(
    src: &[u8],
    repository: &fasta::Repository,
    workers: Option<usize>,
) -> io::Result<Vec<RecordBuf>> {
    let mut reader = cram_reader(src, repository);
    let header = reader.read_header()?;

    match workers.and_then(NonZero::new) {
        Some(workers) => reader.multithreaded_records(&header, workers).collect(),
        None => reader.records(&header).collect(),
    }
}

/// Reads the records up to the first error, and returns them with the error.
fn read_until_error(
    src: &[u8],
    repository: &fasta::Repository,
    workers: Option<usize>,
) -> io::Result<(Vec<RecordBuf>, Option<io::Error>)> {
    let mut reader = cram_reader(src, repository);
    let header = reader.read_header()?;

    let results: Box<dyn Iterator<Item = io::Result<RecordBuf>>> =
        match workers.and_then(NonZero::new) {
            Some(workers) => Box::new(reader.multithreaded_records(&header, workers)),
            None => Box::new(reader.records(&header)),
        };

    let mut records = Vec::new();

    for result in results {
        match result {
            Ok(record) => records.push(record),
            Err(e) => return Ok((records, Some(e))),
        }
    }

    Ok((records, None))
}

/// Mapped records on `sq0` (one with a CIGAR of more than 65,535 operations) and on `sq1`,
/// which the repository lacks, so its slices embed a reference, and unmapped records, in
/// reference-free slices.
fn records(reference: &[u8]) -> Vec<RecordBuf> {
    use Kind::{Deletion, Insertion, Match, SoftClip};

    let mut records = Vec::new();

    for i in 0..300 {
        let start = 1 + i * 97 % (SQ0_LEN - 200);
        let ops: &[(Kind, usize)] = match i % 3 {
            0 => &[(Match, 120)],
            1 => &[(SoftClip, 4), (Match, 40), (Insertion, 2), (Match, 30)],
            _ => &[(Match, 50), (Deletion, 3), (Match, 60), (SoftClip, 7)],
        };
        records.push(mapped(reference, 0, start, ops, &format!("m{i}")));
    }

    // A CIGAR too long for BAM's 16-bit operation count.
    let ops: Vec<_> = (0..35_000)
        .flat_map(|_| [(Match, 1), (Insertion, 1)])
        .collect();
    records.push(mapped(reference, 0, 1001, &ops, "long-cigar"));

    records.sort_by_key(|record| record.alignment_start());

    for i in 0..40 {
        let ops = [(Match, 80)];
        records.push(mapped(reference, 1, 1 + i * 13, &ops, &format!("sq1-{i}")));
    }

    for i in 0..60 {
        let sequence: Vec<u8> = (0..1 + i % 37).map(|j| b"ACGTN"[(i + j) % 5]).collect();
        records.push(
            RecordBuf::builder()
                .set_name(format!("u{i}"))
                .set_flags(Flags::UNMAPPED)
                .set_quality_scores(vec![20; sequence.len()].into())
                .set_sequence(sequence.into())
                .build(),
        );
    }

    records
}

/// A record on `reference_sequence_id` at `start` with `ops`, its bases from `sq0`'s reference
/// with a mismatch every 17 bases.
fn mapped(
    reference: &[u8],
    reference_sequence_id: usize,
    start: usize,
    ops: &[(Kind, usize)],
    name: &str,
) -> RecordBuf {
    let mut sequence = Vec::new();
    let mut i = start - 1;

    for &(kind, len) in ops {
        match kind {
            Kind::Match => {
                sequence.extend(&reference[i..i + len]);
                i += len;
            }
            Kind::Deletion => i += len,
            _ => sequence.extend(b"T".repeat(len)),
        }
    }

    for j in (2..sequence.len()).step_by(17) {
        sequence[j] = if sequence[j] == b'A' { b'C' } else { b'A' };
    }

    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(reference_sequence_id)
        .set_alignment_start(Position::new(start).unwrap())
        .set_mapping_quality(MappingQuality::new(60).unwrap())
        .set_cigar(ops.iter().map(|&(kind, len)| Op::new(kind, len)).collect())
        .set_quality_scores(vec![30; sequence.len()].into())
        .set_sequence(sequence.into())
        .set_data(
            [
                (Tag::from(*b"XI"), Value::from(7)),
                (Tag::from(*b"XZ"), Value::from(name)),
            ]
            .into_iter()
            .collect(),
        )
        .build()
}

/// A pseudorandom reference.
fn reference_sequence(len: usize) -> Vec<u8> {
    let mut state = 7u32;

    (0..len)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            b"ACGT"[(state >> 16) as usize % 4]
        })
        .collect()
}

fn header() -> sam::Header {
    format!("@SQ\tSN:sq0\tLN:{SQ0_LEN}\n@SQ\tSN:sq1\tLN:{SQ0_LEN}\n")
        .parse()
        .expect("valid header")
}

/// A repository with `sq0` only.
fn repository(reference: &[u8]) -> fasta::Repository {
    let definition = fasta::record::Definition::new("sq0", None);
    let record = fasta::Record::new(definition, reference.to_vec().into());
    fasta::Repository::new(vec![record])
}

fn fasta_repository(path: &PathBuf) -> io::Result<fasta::Repository> {
    let mut reader = fasta::io::Reader::new(io::BufReader::new(fs::File::open(path)?));
    let records = reader.records().collect::<io::Result<Vec<_>>>()?;
    Ok(fasta::Repository::new(records))
}

fn write(
    header: &sam::Header,
    reference: &[u8],
    records: &[RecordBuf],
    version: Version,
    profile: Profile,
) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut options = Options::new(version, profile);
    options.records_per_slice = 23;

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository(reference))
        .set_compression_options(options)
        .build_from_writer(Vec::new());

    writer.write_header(header)?;

    for record in records {
        writer.write_alignment_record(header, record)?;
    }

    writer.try_finish(header)?;

    Ok(writer.into_inner())
}
