mod container;
mod file_id;
mod format_version;
mod magic_number;

use std::{
    io::{self, Write},
    thread,
};

use crate::fasta;
use crate::sam::{self as sam, header::ReferenceSequences};

use self::{
    container::write_container, file_id::write_file_id, format_version::write_format_version,
    magic_number::write_magic_number,
};
use crate::cram::{FileDefinition, NormalizedSequenceDigest, calculate_normalized_sequence_digest};

pub fn write_header<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    file_definition: &FileDefinition,
    header: &sam::Header,
    threads: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_file_definition(writer, file_definition)?;
    write_file_header(writer, reference_sequence_repository, header, threads)?;
    Ok(())
}

/// Writes the file header container of `header` as it is.
pub(crate) fn write_header_container<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    write_container(writer, header)
}

pub fn write_file_definition<W>(writer: &mut W, file_definition: &FileDefinition) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;
    write_format_version(writer, file_definition.version())?;
    write_file_id(writer, file_definition.file_id())?;
    Ok(())
}

/// Writes the file header container, calculating the missing reference sequence checksums on up
/// to `threads` threads.
pub fn write_file_header<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    threads: usize,
) -> io::Result<()>
where
    W: Write,
{
    let mut header = header.clone();
    let releases_sequences = super::is_coordinate_sorted(&header);

    add_missing_reference_sequence_checksums(
        reference_sequence_repository,
        header.reference_sequences_mut(),
        releases_sequences,
        threads,
    )?;

    write_container(writer, &header)
}

/// Adds the missing MD5 checksums of the reference sequences, a sequence at a time on each of up
/// to `threads` threads.
///
/// A sequence the repository has not cached is read a part at a time where the repository reads
/// ranges (an indexed FASTA), so that a thread holds a part of it rather than all; otherwise it is
/// loaded, and then released from the cache if `releases_sequences` (see
/// `Writer::use_reference_sequence`).
pub(crate) fn add_missing_reference_sequence_checksums(
    reference_sequence_repository: &fasta::Repository,
    reference_sequences: &mut ReferenceSequences,
    releases_sequences: bool,
    threads: usize,
) -> io::Result<()> {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use crate::sam::header::record::value::map::reference_sequence::{Md5Checksum, tag};

    let missing: Vec<(usize, Vec<u8>)> = reference_sequences
        .iter()
        .enumerate()
        .filter(|(_, (_, map))| !map.other_fields().contains_key(&tag::MD5_CHECKSUM))
        .map(|(i, (name, _))| (i, name.to_vec()))
        .collect();

    let calculate = |name: &[u8]| {
        calculate_checksum(
            reference_sequence_repository,
            name,
            releases_sequences,
            CHECKSUM_PART_LEN,
        )
    };

    let checksums: Vec<io::Result<Option<[u8; 16]>>> = if threads <= 1 || missing.len() <= 1 {
        missing.iter().map(|(_, name)| calculate(name)).collect()
    } else {
        let next = AtomicUsize::new(0);

        let mut results: Vec<_> = thread::scope(|scope| {
            let handles: Vec<_> = (0..threads.min(missing.len()))
                .map(|_| {
                    scope.spawn(|| {
                        let mut results = Vec::new();

                        while let Some((_, name)) =
                            missing.get(next.fetch_add(1, Ordering::Relaxed))
                        {
                            results.push((name.as_slice(), calculate(name)));
                        }

                        results
                    })
                })
                .collect();

            handles
                .into_iter()
                .flat_map(|handle| handle.join().expect("a checksum thread does not panic"))
                .collect()
        });

        // In the order of the sequences.
        let order: std::collections::HashMap<&[u8], usize> = missing
            .iter()
            .enumerate()
            .map(|(k, (_, name))| (name.as_slice(), k))
            .collect();
        results.sort_by_key(|(name, _)| order[name]);
        results.into_iter().map(|(_, result)| result).collect()
    };

    for ((i, _), checksum) in missing.iter().zip(checksums) {
        // As in HTSlib, a sequence missing from the reference gets no checksum.
        if let Some(checksum) = checksum? {
            let (_, map) = reference_sequences
                .get_index_mut(*i)
                .expect("the index is of a reference sequence");

            let fields = map.other_fields_mut();
            let (j, _) = fields.insert_full(
                tag::MD5_CHECKSUM,
                Md5Checksum::from(checksum).to_string().into(),
            );

            // Before a `UR` field, where HTSlib puts it.
            if let Some(k) = fields.get_index_of(&tag::URI)
                && k < j
            {
                fields.move_index(j, k);
            }
        }
    }

    Ok(())
}

/// The bases read at a time to calculate a checksum.
const CHECKSUM_PART_LEN: usize = 1 << 22;

/// Calculates the checksum of the reference sequence `name`, or `None` if the repository lacks it.
///
/// Where the repository reads ranges, this reads `part_len` bases at a time and caches none.
/// Otherwise it reads the whole sequence: unsorted records keep every sequence loaded, so it is
/// cached; sorted records load each in turn and release it, and this may run while they do, so it
/// is not.
fn calculate_checksum(
    reference_sequence_repository: &fasta::Repository,
    name: &[u8],
    releases_sequences: bool,
    part_len: usize,
) -> io::Result<Option<[u8; 16]>> {
    use crate::core::Position;

    if !reference_sequence_repository.is_cached(name) {
        let position = |n: usize| {
            Position::new(n)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid position"))
        };

        let mut digest = NormalizedSequenceDigest::default();
        let mut start = 0;

        // A part ends early at the end of the sequence.
        while let Some(part) = reference_sequence_repository
            .get_range(
                name,
                (position(start + 1)?..=position(start + part_len)?).into(),
            )
            .transpose()?
        {
            digest.update(&part[..]);

            if part.len() < part_len {
                return Ok(Some(digest.finalize()));
            }

            start += part_len;
        }
    }

    if releases_sequences {
        reference_sequence_repository
            .normalized_md5(name)
            .transpose()
    } else {
        let sequence = reference_sequence_repository.get(name).transpose()?;
        Ok(sequence.map(|sequence| calculate_normalized_sequence_digest(&sequence[..])))
    }
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

    #[test]
    fn test_add_missing_reference_sequence_checksums() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZero;

        use fasta::record::{Definition, Sequence};
        use sam::header::record::value::{
            Map,
            map::{ReferenceSequence, reference_sequence::tag},
        };

        let reference_sequences = vec![
            fasta::Record::new(
                Definition::new("sq0", None),
                Sequence::from(b"TTCACCCA".to_vec()),
            ),
            fasta::Record::new(
                Definition::new("sq1", None),
                Sequence::from(b"GATCTTACTTTTT".to_vec()),
            ),
        ];

        let sq0_md5_checksum = BString::from("be19336b7e15968f7ac7dc82493d9cd8");
        let sq1_md5_checksum = BString::from("d80f22a19aeeb623b3e4f746c762f21d");

        let repository = fasta::Repository::new(reference_sequences);

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::builder()
                    .set_length(const { NonZero::new(13).unwrap() })
                    .insert(tag::MD5_CHECKSUM, sq1_md5_checksum.clone())
                    .build()?,
            )
            .build();

        for threads in [1, 4] {
            let mut header = header.clone();

            add_missing_reference_sequence_checksums(
                &repository,
                header.reference_sequences_mut(),
                false,
                threads,
            )?;

            let sq0 = header.reference_sequences().get(&b"sq0"[..]);
            assert_eq!(
                sq0.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
                Some(&sq0_md5_checksum)
            );

            let sq1 = header.reference_sequences().get(&b"sq1"[..]);
            assert_eq!(
                sq1.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
                Some(&sq1_md5_checksum)
            );
        }

        Ok(())
    }

    #[test]
    fn test_calculate_checksum_reads_ranges_of_an_indexed_reference() -> io::Result<()> {
        use std::io::Cursor;

        use fasta::{io::Indexer, repository::adapters::IndexedReader};

        let src = b">sq0\nTTCA\nCCCA\n>sq1\nGATcttact\nTTTT\n".to_vec();

        let mut indexer = Indexer::new(Cursor::new(&src));
        let mut records = Vec::new();

        while let Some(record) = indexer.index_record().map_err(io::Error::other)? {
            records.push(record);
        }

        let reader = fasta::io::IndexedReader::new(Cursor::new(src), records.into());
        let repository = fasta::Repository::new(IndexedReader::new(reader));

        for part_len in [1, 3, 100] {
            let checksum = calculate_checksum(&repository, b"sq1", false, part_len)?;
            assert_eq!(
                checksum,
                Some(calculate_normalized_sequence_digest(b"GATCTTACTTTTT"))
            );
            assert!(!repository.is_cached(b"sq1"));
        }

        assert_eq!(calculate_checksum(&repository, b"sq2", false, 3)?, None);

        Ok(())
    }
}
