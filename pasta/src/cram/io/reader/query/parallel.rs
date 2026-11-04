//! Region queries decoded on several threads.

use std::{
    collections::BTreeMap,
    io::{self, Read, Seek, SeekFrom},
    num::NonZeroUsize,
    panic::{self, AssertUnwindSafe},
    sync::Arc,
    thread,
};

use crossbeam_channel::{Receiver, Sender};

use super::{index_record_intersects, record_intersects};
use crate::core::{Region, region::Interval};
use crate::cram::{
    self as cram, crai,
    io::reader::{Container, Fields, MdNm, Reader, Tags},
};
use crate::fasta;
use crate::sam;
use crate::util::threads;

/// Slices in flight per worker: read and not yet merged, whether queued, being decoded, or
/// decoded and waiting for an earlier slice. This bounds memory, even behind a slow slice, while
/// keeping the workers busy.
pub(crate) const IN_FLIGHT_PER_WORKER: usize = 4;

/// A decoded slice's accumulator, or the panic of the worker that decoded it.
type Outcome<T> = thread::Result<io::Result<T>>;

/// A slice to decode for a region: the container holding it (shared by the region's other
/// slices in it) and its landmark.
struct Task {
    region: usize,
    container: Arc<Container>,
    landmark: u64,
}

/// What decoding a slice takes besides the slice.
struct Decoder<'a> {
    header: &'a sam::Header,
    repository: fasta::Repository,
    md_nm: MdNm,
    fields: Fields,
    tags: Tags,
}

impl<R> Reader<R>
where
    R: Read + Seek + Send,
{
    /// Decodes the records that intersect each of `regions` on `workers` threads, in the order
    /// [`Self::query`] gives them region by region.
    ///
    /// Each slice a region's query reads is a task. `fold` is called on a worker with the task's
    /// accumulator, the index of the task's region in `regions` and each record of the slice
    /// that intersects the region. `merge` is called on this thread with each task's
    /// accumulator, in query order, so folding and then merging gives what folding every record
    /// of every query in turn would. A record in more than one region is folded for each. At
    /// most a few slices per worker are read and not yet merged at any time, however slow one
    /// slice is to decode.
    ///
    /// Records are decoded as the reader decodes them ([`Self::set_required_fields`],
    /// [`Self::set_required_tags`], [`Self::set_md_nm`]). To read them as BAM records, encode
    /// them with [`cram::Record::encode_bam_record`] in `fold`.
    ///
    /// # Errors
    ///
    /// Returns an error if a region names a reference sequence the header does not have, a
    /// container cannot be read or decoded, `fold` fails, or `merge` fails. Decoding stops at the
    /// first error in query order.
    ///
    /// # Panics
    ///
    /// If `fold` panics, the panic is resumed on this thread once the slices before it are
    /// merged.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io, num::NonZeroUsize};
    /// use pasta::cram::{self as cram, crai};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let regions = ["sq0:8-13".parse()?];
    /// let workers = NonZeroUsize::new(4).unwrap();
    ///
    /// let mut n = 0;
    ///
    /// reader.fold_queries(
    ///     &header,
    ///     &index,
    ///     &regions,
    ///     workers,
    ///     |count: &mut u64, _, _| {
    ///         *count += 1;
    ///         Ok(())
    ///     },
    ///     |count| {
    ///         n += count;
    ///         Ok::<_, io::Error>(())
    ///     },
    /// )?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn fold_queries<T, F, M, E>(
        &mut self,
        header: &sam::Header,
        index: &crai::Index,
        regions: &[Region],
        workers: NonZeroUsize,
        fold: F,
        mut merge: M,
    ) -> Result<(), E>
    where
        T: Default + Send,
        F: Fn(&mut T, usize, &cram::Record<'_>) -> io::Result<()> + Sync,
        M: FnMut(T) -> Result<(), E>,
        E: From<io::Error>,
    {
        let queries = resolve(header, regions)?;

        let decoder = Decoder {
            header,
            repository: self.reference_sequence_repository.clone(),
            md_nm: self.md_nm(),
            fields: self.required_fields(),
            tags: self.required_tags().clone(),
        };

        thread::scope(|scope| {
            // The workers start first, fewer of them when the system cannot start them all.
            let idle = threads::start_idle_scoped(scope, workers)?;
            let in_flight = idle.len() * IN_FLIGHT_PER_WORKER;
            let (work_tx, work_rx) = crossbeam_channel::bounded::<(u64, Task)>(in_flight);
            let (done_tx, done_rx) = crossbeam_channel::bounded::<(u64, Outcome<T>)>(in_flight);

            for t in idle {
                let work_rx = work_rx.clone();
                let done_tx = done_tx.clone();
                let (decoder, queries, fold) = (&decoder, &queries, &fold);

                t.run(move || {
                    for (id, task) in work_rx {
                        let (reference_sequence_id, interval) = queries[task.region];

                        // A panic is sent on to be raised on the calling thread in order: a
                        // slice that never arrives would leave it waiting forever.
                        let outcome = panic::catch_unwind(AssertUnwindSafe(|| {
                            decoder.decode(&task, reference_sequence_id, interval, fold)
                        }));

                        if done_tx.send((id, outcome)).is_err() {
                            // The caller stopped; it reports why.
                            break;
                        }
                    }
                });
            }

            drop(work_rx);

            // One permit per slice that may be in flight; this thread returns each once the
            // slice is merged. Neither channel can then hold more than `in_flight` slices.
            let (permits_tx, permits_rx) = crossbeam_channel::bounded(in_flight);

            for _ in 0..in_flight {
                permits_tx
                    .send(())
                    .expect("the permit channel has room for every permit");
            }

            {
                let done_tx = done_tx.clone();
                let queries = &queries;

                // Returning the error drops the work channel, which stops the workers.
                thread::Builder::new()
                    .spawn_scoped(scope, move || {
                        // A read error takes the place of the next task, so that it is reported
                        // after every task before it.
                        if let Err((id, e)) = self.send_tasks(index, queries, &permits_rx, &work_tx)
                        {
                            let _ = done_tx.send((id, Ok(Err(e))));
                        }
                    })
                    .map_err(threads::failed)?;
            }

            // Only the threads hold this now, so the channel closes when they finish.
            drop(done_tx);

            let mut pending = BTreeMap::new();
            let mut next = 0;

            for (id, outcome) in &done_rx {
                pending.insert(id, outcome);

                while let Some(outcome) = pending.remove(&next) {
                    let result = outcome.unwrap_or_else(|payload| panic::resume_unwind(payload));
                    merge(result?)?;
                    next += 1;

                    // The slice is merged: let another one in. The reader may have finished.
                    let _ = permits_tx.send(());
                }
            }

            Ok(())
        })
    }

    /// Reads the containers of every query's slices in order, sending one task per slice. A
    /// container is read once for consecutive slices in it, even across queries. Each task, and
    /// the read before it, waits for a permit. On a read error, returns it with the ID the next
    /// task would have had, in that task's place.
    fn send_tasks(
        &mut self,
        index: &crai::Index,
        queries: &[(usize, Interval)],
        permits: &Receiver<()>,
        tx: &Sender<(u64, Task)>,
    ) -> Result<(), (u64, io::Error)> {
        let mut id = 0;
        let mut last: Option<(u64, Arc<Container>)> = None;
        // A permit taken for a task that was not sent (the region ended) is kept for the next.
        let mut has_permit = false;

        for (region, &(reference_sequence_id, interval)) in queries.iter().enumerate() {
            let index_records = index
                .iter()
                .filter(|record| index_record_intersects(record, reference_sequence_id, interval));

            for index_record in index_records {
                // The permits run out (the channel closes) only when the caller stopped; it
                // reports why.
                if !has_permit {
                    if permits.recv().is_err() {
                        return Ok(());
                    }

                    has_permit = true;
                }

                let offset = index_record.offset();

                let container = match &last {
                    Some((last_offset, container)) if *last_offset == offset => {
                        Arc::clone(container)
                    }
                    _ => {
                        last = None;

                        let Some(container) =
                            self.read_container_at(offset).map_err(|e| (id, e))?
                        else {
                            // As in a query, the end of the stream ends the region.
                            break;
                        };

                        let container = Arc::new(container);
                        last = Some((offset, Arc::clone(&container)));
                        container
                    }
                };

                let task = Task {
                    region,
                    container,
                    landmark: index_record.landmark(),
                };

                if tx.send((id, task)).is_err() {
                    // The workers stopped; the caller reports why.
                    return Ok(());
                }

                has_permit = false;
                id += 1;
            }
        }

        Ok(())
    }

    /// Reads the container at `offset`, or `None` at the end of the stream.
    fn read_container_at(&mut self, offset: u64) -> io::Result<Option<Container>> {
        self.seek(SeekFrom::Start(offset))?;

        let mut container = Container::default();

        match self.read_container(&mut container)? {
            0 => Ok(None),
            _ => Ok(Some(container)),
        }
    }
}

impl Decoder<'_> {
    /// Decodes a task's slice and folds the records that intersect the region into a new
    /// accumulator.
    fn decode<T, F>(
        &self,
        task: &Task,
        reference_sequence_id: usize,
        interval: Interval,
        fold: &F,
    ) -> io::Result<T>
    where
        T: Default,
        F: Fn(&mut T, usize, &cram::Record<'_>) -> io::Result<()>,
    {
        let compression_header = task.container.compression_header()?;
        let slice = task.container.slice_at_landmark(task.landmark)?;
        let (core_data_src, external_data_srcs) =
            slice.decode_blocks_for(&compression_header, self.fields, &self.tags)?;

        let records = slice.records_with_fields(
            self.repository.clone(),
            self.header,
            &compression_header,
            &core_data_src,
            &external_data_srcs,
            self.md_nm,
            self.fields,
            &self.tags,
        )?;

        let mut acc = T::default();

        for record in &records {
            if record_intersects(record, reference_sequence_id, interval) {
                fold(&mut acc, task.region, record)?;
            }
        }

        Ok(acc)
    }
}

/// Resolves each region's reference sequence name to its index in the header.
fn resolve(header: &sam::Header, regions: &[Region]) -> io::Result<Vec<(usize, Interval)>> {
    regions
        .iter()
        .map(|region| {
            let reference_sequence_id = header
                .reference_sequences()
                .get_index_of(region.name())
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "invalid reference sequence name",
                    )
                })?;

            Ok((reference_sequence_id, region.interval()))
        })
        .collect()
}
