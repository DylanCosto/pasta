//! Reconstructing records in bulk matches the generic, per-base conversion.

use std::{error::Error, io};

use pasta::bam;
use pasta::core::Position;
use pasta::cram;
use pasta::fasta;
use pasta::sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::data::field::{Value, value::Array},
    },
};

const REFERENCE_LEN: usize = 1000;

#[test]
fn test_records_match_generic_conversion() -> Result<(), Box<dyn Error>> {
    let reference = reference_sequence(REFERENCE_LEN);
    let records = records(&reference);
    let header = header();
    let src = write(&header, &reference, &records)?;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository(&reference))
        .build_from_reader(&src[..]);
    reader.read_header()?;
    let actual: Vec<RecordBuf> = reader.records(&header).collect::<io::Result<_>>()?;

    let mut expected = Vec::new();
    let mut sam_via_view = Vec::new();
    let mut sam_via_record = Vec::new();

    let mut reader = cram::io::Reader::new(&src[..]);
    reader.read_header()?;
    let mut container = cram::io::reader::Container::default();
    let mut sequence = Vec::new();
    let mut bam_record = Vec::new();

    while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            let records = slice.records(
                repository(&reference),
                &header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;

            for record in &records {
                expected.push(RecordBuf::try_from_alignment_record(&header, record)?);

                let mut writer = sam::io::Writer::new(&mut sam_via_record);
                writer.write_alignment_record(&header, record)?;

                let view = record.with_sequence(&mut sequence);
                let mut writer = sam::io::Writer::new(&mut sam_via_view);
                writer.write_alignment_record(&header, &view)?;

                let mut generic_bam = bam::io::Writer::from(Vec::new());
                generic_bam.write_alignment_record(&header, &view)?;
                record.encode_bam_record(&header, &mut bam_record)?;
                assert_eq!(bam_record, generic_bam.get_ref()[size_of::<u32>()..]);
            }
        }
    }

    assert_eq!(actual.len(), records.len());
    assert_eq!(actual, expected);
    assert_eq!(sam_via_view, sam_via_record);

    Ok(())
}

#[test]
fn test_records_with_required_fields_match_full_decoding() -> Result<(), Box<dyn Error>> {
    use cram::io::reader::Fields;

    let reference = reference_sequence(REFERENCE_LEN);
    let records = records(&reference);
    let header = header();
    let src = write(&header, &reference, &records)?;

    let read = |fields: Fields| -> io::Result<Vec<RecordBuf>> {
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository(&reference))
            .set_required_fields(fields)
            .build_from_reader(&src[..]);
        reader.read_header()?;
        reader.records(&header).collect()
    };

    let expected = read(Fields::all())?;

    for fields in [
        Fields::empty(),
        Fields::FLAGS,
        Fields::FLAGS | Fields::MAPPING_QUALITY,
        Fields::REFERENCE_SEQUENCE_ID | Fields::ALIGNMENT_START | Fields::CIGAR,
        Fields::NAME,
        Fields::DATA,
        Fields::READ_GROUP,
        Fields::DATA | Fields::READ_GROUP,
        Fields::SEQUENCE | Fields::DATA,
        Fields::SEQUENCE | Fields::DATA | Fields::READ_GROUP,
        Fields::SEQUENCE,
        Fields::QUALITY_SCORES,
    ] {
        let actual = read(fields)?;
        assert_eq!(actual.len(), expected.len(), "{fields:?}");

        for (actual, expected) in actual.iter().zip(&expected) {
            assert_eq!(actual.flags(), expected.flags(), "{fields:?}");

            if fields.contains(Fields::MAPPING_QUALITY) {
                assert_eq!(actual.mapping_quality(), expected.mapping_quality());
            }

            if fields.contains(Fields::CIGAR) {
                assert_eq!(
                    actual.reference_sequence_id(),
                    expected.reference_sequence_id()
                );
                assert_eq!(actual.alignment_start(), expected.alignment_start());
                assert_eq!(actual.cigar(), expected.cigar(), "{fields:?}");
            }

            if fields.contains(Fields::NAME) {
                assert_eq!(actual.name(), expected.name());
            }

            // The read group is a data series of its own, and missing MD and NM tags are only
            // rebuilt with the sequence and the data.
            let expected_data: Vec<_> = expected
                .data()
                .iter()
                .filter(|(tag, _)| {
                    if *tag == Tag::READ_GROUP {
                        fields.contains(Fields::READ_GROUP)
                    } else if [Tag::MISMATCHED_POSITIONS, Tag::EDIT_DISTANCE].contains(tag) {
                        fields.contains(Fields::SEQUENCE | Fields::DATA)
                    } else {
                        fields.contains(Fields::DATA)
                    }
                })
                .collect();
            assert_eq!(
                actual.data().iter().collect::<Vec<_>>(),
                expected_data,
                "{fields:?}"
            );

            if fields.contains(Fields::SEQUENCE) {
                assert_eq!(actual.sequence(), expected.sequence());
            } else {
                assert!(actual.sequence().is_empty(), "{fields:?}");
            }

            if fields.contains(Fields::QUALITY_SCORES) {
                assert_eq!(actual.quality_scores(), expected.quality_scores());
            }
        }
    }

    Ok(())
}

#[test]
fn test_decode_blocks_for_skips_unread_blocks() -> Result<(), Box<dyn Error>> {
    use cram::io::reader::{Fields, Tags};

    let reference = reference_sequence(REFERENCE_LEN);
    let header = header();
    let src = write(&header, &reference, &records(&reference))?;

    let mut reader = cram::io::Reader::new(&src[..]);
    reader.read_header()?;
    let mut container = cram::io::reader::Container::default();

    while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;
            let (_, all) = slice.decode_blocks()?;
            let (_, flags) =
                slice.decode_blocks_for(&compression_header, Fields::FLAGS, &Tags::All)?;
            let (_, every) =
                slice.decode_blocks_for(&compression_header, Fields::all(), &Tags::All)?;
            let (_, no_tags) = slice.decode_blocks_for(
                &compression_header,
                Fields::all(),
                &Tags::Only(Vec::new()),
            )?;

            assert!(flags.len() < all.len(), "{} < {}", flags.len(), all.len());
            assert_eq!(every.len(), all.len());
            assert!(
                no_tags.len() < all.len(),
                "{} < {}",
                no_tags.len(),
                all.len()
            );
        }
    }

    Ok(())
}

/// Records exercising every kind of read feature: substitutions (including against lowercase
/// reference bases), insertions, deletions, skips, soft clips, missing sequences, and unmapped
/// reads, with integer, string, array, and read group data fields.
fn records(reference: &[u8]) -> Vec<RecordBuf> {
    use Kind::{Deletion, Insertion, Match, Skip, SoftClip};

    let mut records = Vec::new();

    for (i, ops) in [
        &[(Match, 50)][..],
        &[
            (SoftClip, 3),
            (Match, 20),
            (Insertion, 2),
            (Match, 10),
            (Deletion, 4),
            (Match, 15),
        ],
        &[(Match, 10), (Skip, 30), (Match, 10), (SoftClip, 5)],
    ]
    .into_iter()
    .enumerate()
    {
        for start in [1, 100, 401, 900] {
            let mut sequence = sequence_for(reference, start, ops);

            // Mismatches against the reference become substitution features.
            for j in (2..sequence.len()).step_by(17) {
                sequence[j] = match sequence[j].to_ascii_uppercase() {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }

            records.push(
                RecordBuf::builder()
                    .set_name(format!("r{i}-{start}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::new(start).unwrap())
                    .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
                    .set_cigar(ops.iter().map(|&(kind, len)| Op::new(kind, len)).collect())
                    .set_quality_scores(vec![30; sequence.len()].into())
                    .set_sequence(sequence.into())
                    .set_data(
                        [
                            (Tag::from(*b"XI"), Value::from(7)),
                            (Tag::from(*b"XZ"), Value::from("a b~")),
                            (Tag::READ_GROUP, Value::from("rg0")),
                            (
                                Tag::from(*b"XB"),
                                Value::Array(Array::UInt8(vec![1, 2, 255])),
                            ),
                        ]
                        .into_iter()
                        .collect(),
                    )
                    .build(),
            );
        }
    }

    records.sort_by_key(|record| record.alignment_start());

    records.push(
        RecordBuf::builder()
            .set_name("unmapped")
            .set_flags(Flags::UNMAPPED)
            .set_sequence(b"ACGTN".to_vec().into())
            .set_quality_scores(vec![20; 5].into())
            .build(),
    );

    records
}

fn sequence_for(reference: &[u8], start: usize, ops: &[(Kind, usize)]) -> Vec<u8> {
    let mut sequence = Vec::new();
    let mut i = start - 1;

    for &(kind, len) in ops {
        match kind {
            Kind::Match => {
                sequence.extend(&reference[i..i + len]);
                i += len;
            }
            Kind::Deletion | Kind::Skip => i += len,
            _ => sequence.extend(b"N".repeat(len)),
        }
    }

    sequence
}

/// A pseudorandom reference with a lowercase (soft-masked) stretch.
fn reference_sequence(len: usize) -> Vec<u8> {
    let mut state = 7u32;

    (0..len)
        .map(|i| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            let base = b"ACGT"[(state >> 16) as usize % 4];
            if (380..520).contains(&i) {
                base.to_ascii_lowercase()
            } else {
                base
            }
        })
        .collect()
}

fn header() -> sam::Header {
    format!("@SQ\tSN:sq0\tLN:{REFERENCE_LEN}\n@RG\tID:rg0\n")
        .parse()
        .expect("valid header")
}

fn repository(reference: &[u8]) -> fasta::Repository {
    let definition = fasta::record::Definition::new("sq0", None);
    let record = fasta::Record::new(definition, reference.to_vec().into());
    fasta::Repository::new(vec![record])
}

fn write(
    header: &sam::Header,
    reference: &[u8],
    records: &[RecordBuf],
) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository(reference))
        .build_from_writer(Vec::new());

    writer.write_header(header)?;

    for record in records {
        writer.write_alignment_record(header, record)?;
    }

    writer.try_finish(header)?;

    Ok(writer.into_inner())
}
