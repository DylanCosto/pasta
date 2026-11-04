//! CRAM decoding across a worker pool, one container per task.
//!
//! A reader thread reads containers, workers decode them and fold each record into a
//! per-container value (SAM text, raw BAM, counts), and the calling thread receives the values in
//! container order. Containers are independent, so this is equivalent to decoding serially.
//! Indexed queries are decoded the same way, a slice per task, by pasta.

use std::{
    collections::BTreeMap,
    io::{self, Read, Seek},
    num::NonZeroUsize,
    thread,
};

use anyhow::Context;
use pasta::cram::{
    self as cram,
    io::reader::{Fields, MdNm, Tags},
};
use pasta::fasta;
use pasta::sam;

use crate::alignment::CramQueries;

/// Containers in flight per worker: read and not yet merged, whether queued, being decoded, or
/// decoded and waiting for an earlier container. This bounds memory, even behind a slow
/// container, while keeping the workers busy.
const IN_FLIGHT_PER_WORKER: usize = 4;

/// The records to decode.
#[derive(Clone, Copy)]
pub enum Source<'s, 'q> {
    /// Every record, in file order.
    All,
    /// The records each indexed query visits, query by query.
    Queries(&'s CramQueries<'q>),
}

/// Decodes the records of `source` on `workers` threads.
///
/// `fold` is called on a worker for each record with an accumulator per container (or, for
/// queries, per slice); `merge` is called on this thread with each accumulator, in order. An
/// error from `merge` stops decoding.
pub fn decode<R, T, F, M>(
    reader: &mut cram::io::Reader<R>,
    header: &sam::Header,
    repository: &fasta::Repository,
    source: Source<'_, '_>,
    workers: usize,
    fold: F,
    merge: M,
) -> anyhow::Result<()>
where
    R: Read + Seek + Send,
    T: Default + Send,
    F: Fn(&mut T, &cram::Record<'_>) -> io::Result<()> + Sync,
    M: FnMut(T) -> anyhow::Result<()>,
{
    let workers = NonZeroUsize::new(workers).unwrap_or(NonZeroUsize::MIN);

    match source {
        Source::All => decode_all(reader, header, repository, workers, fold, merge),
        Source::Queries(queries) => reader.fold_queries(
            header,
            queries.index(),
            queries.regions(),
            workers,
            |acc, i, record| {
                if queries.visits(i, record, header)? {
                    fold(acc, record)
                } else {
                    Ok(())
                }
            },
            merge,
        ),
    }
}

/// Decodes every container on `workers` threads.
fn decode_all<R, T, F, M>(
    reader: &mut cram::io::Reader<R>,
    header: &sam::Header,
    repository: &fasta::Repository,
    workers: NonZeroUsize,
    fold: F,
    mut merge: M,
) -> anyhow::Result<()>
where
    R: Read + Send,
    T: Default + Send,
    F: Fn(&mut T, &cram::Record<'_>) -> io::Result<()> + Sync,
    M: FnMut(T) -> anyhow::Result<()>,
{
    let workers = workers.get();
    let md_nm = reader.md_nm();
    let fields = reader.required_fields();
    let tags = reader.required_tags().clone();

    let in_flight = workers * IN_FLIGHT_PER_WORKER;

    thread::scope(|scope| {
        // One permit per container that may be in flight; this thread returns each once the
        // container is merged. Neither channel can then hold more than `in_flight` containers.
        let (permits_tx, permits_rx) = crossbeam_channel::bounded(in_flight);

        for _ in 0..in_flight {
            permits_tx
                .send(())
                .expect("the permit channel has room for every permit");
        }

        let (work_tx, work_rx) =
            crossbeam_channel::bounded::<(u64, cram::io::reader::Container)>(in_flight);
        let (done_tx, done_rx) = crossbeam_channel::bounded::<(u64, io::Result<T>)>(in_flight);

        let read = scope.spawn(move || -> io::Result<()> {
            for id in 0.. {
                // The permits run out (the channel closes) only when the consumer stopped; it
                // reports why.
                if permits_rx.recv().is_err() {
                    break;
                }

                let mut container = cram::io::reader::Container::default();

                if reader.read_container(&mut container)? == 0 {
                    break;
                }

                if work_tx.send((id, container)).is_err() {
                    // The consumer stopped; it reports why.
                    break;
                }
            }

            Ok(())
        });

        for _ in 0..workers {
            let work_rx = work_rx.clone();
            let done_tx = done_tx.clone();
            let fold = &fold;
            let tags = &tags;

            scope.spawn(move || {
                for (id, container) in work_rx {
                    let _guard = PanicGuard { id, done: &done_tx };
                    let result =
                        decode_container(&container, header, repository, md_nm, fields, tags, fold);

                    if done_tx.send((id, result)).is_err() {
                        break;
                    }
                }
            });
        }

        // Only workers hold these now, so the channels close when the workers finish.
        drop(work_rx);
        drop(done_tx);

        let mut pending = BTreeMap::new();
        let mut next = 0;

        for (id, result) in &done_rx {
            pending.insert(id, result);

            while let Some(result) = pending.remove(&next) {
                merge(result.with_context(|| format!("decoding CRAM container {next}"))?)?;
                next += 1;

                // The container is merged: let another one in. The reader may have finished.
                let _ = permits_tx.send(());
            }
        }

        read.join()
            .map_err(|_| anyhow::anyhow!("CRAM reader thread panicked"))?
            .context("reading CRAM containers")
    })
}

/// Sends an error in place of a container if decoding it panics: a container that never
/// arrives would leave the consumer waiting for it, and the reader for permits, forever.
struct PanicGuard<'a, T> {
    id: u64,
    done: &'a crossbeam_channel::Sender<(u64, io::Result<T>)>,
}

impl<T> Drop for PanicGuard<'_, T> {
    fn drop(&mut self) {
        if thread::panicking() {
            let error = io::Error::other("CRAM decoding thread panicked");
            let _ = self.done.send((self.id, Err(error)));
        }
    }
}

fn decode_container<T, F>(
    container: &cram::io::reader::Container,
    header: &sam::Header,
    repository: &fasta::Repository,
    md_nm: MdNm,
    fields: Fields,
    tags: &Tags,
    fold: &F,
) -> io::Result<T>
where
    T: Default,
    F: Fn(&mut T, &cram::Record<'_>) -> io::Result<()>,
{
    let compression_header = container.compression_header()?;
    let mut acc = T::default();

    for result in container.slices() {
        let slice = result?;
        let (core_data_src, external_data_srcs) =
            slice.decode_blocks_for(&compression_header, fields, tags)?;

        let records = slice.records_with_fields(
            repository.clone(),
            header,
            &compression_header,
            &core_data_src,
            &external_data_srcs,
            md_nm,
            fields,
            tags,
        )?;

        for record in &records {
            fold(&mut acc, record)?;
        }
    }

    Ok(acc)
}

#[cfg(test)]
mod tests {
    use std::{
        io::Cursor,
        sync::atomic::{AtomicUsize, Ordering},
        time::Duration,
    };

    use pasta::cram::{
        file_definition::Version,
        io::writer::compression::{Options, Profile},
    };
    use pasta::sam::alignment::{
        Record as _, RecordBuf, io::Write as _, record::Flags, record_buf::QualityScores,
    };

    use super::*;

    /// Containers whose decoding has started.
    static STARTED: AtomicUsize = AtomicUsize::new(0);

    /// The record names of a container. Creating one counts a started container.
    struct Names(Vec<Vec<u8>>);

    impl Default for Names {
        fn default() -> Self {
            STARTED.fetch_add(1, Ordering::SeqCst);
            Self(Vec::new())
        }
    }

    #[test]
    fn test_decode_all_bounds_containers_in_flight() -> anyhow::Result<()> {
        let header = sam::Header::default();
        let mut options = Options::new(Version::new(3, 0), Profile::Normal);
        options.records_per_slice = 2;

        let mut writer = cram::io::writer::Builder::default()
            .set_compression_options(options)
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;

        let expected: Vec<Vec<u8>> = (0..200).map(|i| format!("r{i:03}").into_bytes()).collect();

        for name in &expected {
            let record = RecordBuf::builder()
                .set_name(name.clone())
                .set_flags(Flags::UNMAPPED)
                .set_sequence(b"ACGT".to_vec().into())
                .set_quality_scores(QualityScores::from(vec![30; 4]))
                .build();
            writer.write_alignment_record(&header, &record)?;
        }

        writer.try_finish(&header)?;

        let mut reader = cram::io::Reader::new(Cursor::new(writer.into_inner()));
        let header = reader.read_header()?;
        let workers = NonZeroUsize::new(2).expect("nonzero");

        let mut names = Vec::new();
        let mut merged = 0;
        let mut most_retained = 0;

        // The first container is slow, so the others finish first and wait for it.
        decode_all(
            &mut reader,
            &header,
            &fasta::Repository::default(),
            workers,
            |acc: &mut Names, record| {
                let name = record.name().map(|name| name.to_vec()).unwrap_or_default();

                if name == expected[0] {
                    thread::sleep(Duration::from_millis(200));
                }

                acc.0.push(name);
                Ok(())
            },
            |acc| {
                // Containers started and not merged before this one, this one included.
                most_retained = most_retained.max(STARTED.load(Ordering::SeqCst) - merged);
                merged += 1;
                names.extend(acc.0);
                Ok(())
            },
        )?;

        assert_eq!(names, expected);
        assert_eq!(merged, 100);

        let bound = workers.get() * IN_FLIGHT_PER_WORKER;
        assert!(
            most_retained <= bound,
            "{most_retained} containers retained, bound {bound}"
        );

        Ok(())
    }
}
