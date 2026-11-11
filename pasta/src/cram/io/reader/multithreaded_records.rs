//! An iterator over records that decodes slices on several threads.

use std::{
    collections::VecDeque,
    io::{self, Read},
    num::NonZero,
    panic,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread, vec,
};

use crossbeam_channel::Receiver;

use super::{Container, Fields, MdNm, Reader, Tags, query::parallel::IN_FLIGHT_PER_WORKER};
use crate::bgzf::io::worker;
use crate::cram;
use crate::fasta;
use crate::sam::{self as sam, alignment::RecordBuf};
use crate::util::threads;

/// A slice's records or why they could not be read or decoded, or the panic of the job that
/// decoded them.
type Decoded = thread::Result<io::Result<Vec<RecordBuf>>>;

/// What decoding a slice takes besides the slice, shared by the jobs.
struct Decoder {
    header: sam::Header,
    repository: fasta::Repository,
    md_nm: MdNm,
    fields: Fields,
    tags: Tags,
}

/// An iterator over records of a CRAM reader that decodes slices on several threads.
///
/// This is created by calling [`Reader::multithreaded_records`]. It gives the records
/// [`Reader::records`] gives, in the same order, up to the first error.
///
/// Containers are read on the calling thread, and each slice is decoded by a job on a pool of
/// workers. At most a few slices per worker are read and not yet given out at any time, however
/// slow one slice is to decode.
///
/// A panic of a decoding job is re-raised on the calling thread when the iterator reaches that
/// slice. After an error, the iterator ends.
pub struct MultithreadedRecords<'r, R> {
    reader: &'r mut Reader<R>,
    decoder: Arc<Decoder>,
    /// The workers, or `None` when not even one could be started, so that slices are decoded on
    /// the calling thread.
    pool: Option<threads::Pool>,
    /// The most slices sent for decoding and not yet given out.
    max_in_flight: usize,
    /// The slices sent for decoding, in file order.
    in_flight: VecDeque<Receiver<Decoded>>,
    /// The container whose slices are being sent, and the index of the next one to send.
    container: Option<(Arc<Container>, usize)>,
    records: vec::IntoIter<RecordBuf>,
    /// Whether no more containers are read: the stream ended, or reading or decoding failed.
    is_done: bool,
    /// Tells the jobs left when the iterator is dropped to skip their slices.
    is_dropped: Arc<AtomicBool>,
}

impl<'r, R> MultithreadedRecords<'r, R>
where
    R: Read,
{
    pub(crate) fn new(
        reader: &'r mut Reader<R>,
        header: &sam::Header,
        worker_count: NonZero<usize>,
    ) -> Self {
        let decoder = Decoder {
            header: header.clone(),
            repository: reader.reference_sequence_repository.clone(),
            md_nm: reader.md_nm,
            fields: reader.required_fields,
            tags: reader.required_tags.clone(),
        };

        // Without workers, each slice is decoded on this thread once it is needed.
        let pool = threads::pool(worker_count).ok();
        let max_in_flight = pool
            .as_ref()
            .map_or(1, |pool| pool.current_num_threads() * IN_FLIGHT_PER_WORKER);

        Self {
            reader,
            decoder: Arc::new(decoder),
            pool,
            max_in_flight,
            in_flight: VecDeque::with_capacity(max_in_flight),
            container: None,
            records: Vec::new().into_iter(),
            is_done: false,
            is_dropped: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Reads containers and sends their slices for decoding until as many slices are in flight
    /// as may be. A read error takes the place of the next slice, so that it is given out after
    /// every slice before it.
    fn send_slices(&mut self) {
        while self.in_flight.len() < self.max_in_flight {
            let (container, i) = match self.container.take() {
                Some((container, i)) if i < container.slice_count() => (container, i),
                _ if self.is_done => return,
                _ => {
                    let mut container = Container::default();

                    match self.reader.read_container(&mut container) {
                        Ok(0) => {
                            self.is_done = true;
                            return;
                        }
                        Ok(_) => {
                            // A container without slices is passed over by the next pass.
                            self.container = Some((Arc::new(container), 0));
                            continue;
                        }
                        Err(e) => {
                            self.is_done = true;
                            self.in_flight.push_back(ready(Ok(Err(e))));
                            return;
                        }
                    }
                }
            };

            let rx = self.send_slice(Arc::clone(&container), i);
            self.in_flight.push_back(rx);
            self.container = Some((container, i + 1));
        }
    }

    /// Decodes the `i`th slice of `container` on a worker, or here without workers, and returns
    /// the channel of its records.
    fn send_slice(&self, container: Arc<Container>, i: usize) -> Receiver<Decoded> {
        let (tx, rx) = crossbeam_channel::bounded(1);
        let decoder = Arc::clone(&self.decoder);
        let is_dropped = Arc::clone(&self.is_dropped);

        let job = move || {
            if is_dropped.load(Ordering::Relaxed) {
                return;
            }

            // A panic is sent on to be raised on the calling thread: a pool thread aborts the
            // process on a panic, and a slice that never arrives would end the records early.
            let decoded = worker::catch(|| decoder.decode(&container, i));
            let _ = tx.send(decoded);
        };

        match &self.pool {
            Some(pool) => pool.spawn(job),
            None => job(),
        }

        rx
    }

    /// Stops reading, after an error or a panic.
    fn stop(&mut self) {
        self.is_done = true;
        self.container = None;
        self.in_flight.clear();
        self.is_dropped.store(true, Ordering::Relaxed);
    }
}

impl<R> Iterator for MultithreadedRecords<'_, R>
where
    R: Read,
{
    type Item = io::Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.records.next() {
                return Some(Ok(record));
            }

            // The workers decode the slices after this one while its records are given out.
            self.send_slices();
            let rx = self.in_flight.pop_front()?;

            // Every job sends its result, a panic included, unless the iterator was dropped.
            let result = match rx.recv() {
                Ok(Ok(result)) => result,
                Ok(Err(payload)) => {
                    self.stop();
                    panic::resume_unwind(payload);
                }
                Err(_) => Err(io::Error::other("CRAM decoding job ended without a result")),
            };

            match result {
                Ok(records) => self.records = records.into_iter(),
                Err(e) => {
                    self.stop();
                    return Some(Err(e));
                }
            }
        }
    }
}

impl<R> Drop for MultithreadedRecords<'_, R> {
    fn drop(&mut self) {
        self.is_dropped.store(true, Ordering::Relaxed);
    }
}

impl Decoder {
    /// Decodes the `i`th slice of `container` as the reader would.
    fn decode(&self, container: &Container, i: usize) -> io::Result<Vec<RecordBuf>> {
        let compression_header = container.compression_header()?;
        let slice = container.slice(i)?;
        let (core_data_src, external_data_srcs) =
            slice.decode_blocks_for(&compression_header, self.fields, &self.tags)?;

        let records = slice
            .records_with_fields(
                self.repository.clone(),
                &self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                self.md_nm,
                self.fields,
                &self.tags,
            )?
            .into_iter()
            .map(cram::Record::try_into_record_buf)
            .collect::<io::Result<Vec<_>>>()?;

        #[cfg(test)]
        for record in &records {
            worker::panic_if_marked(record.name().unwrap_or_default());
        }

        Ok(records)
    }
}

/// Returns a channel that already holds `decoded`.
fn ready(decoded: Decoded) -> Receiver<Decoded> {
    let (tx, rx) = crossbeam_channel::bounded(1);
    tx.send(decoded)
        .expect("the channel has room for one result");
    rx
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf::io::worker::{MARKER_PANIC, PANIC_MARKER, assert_panics_with, within_timeout};
    use crate::cram::io::writer::compression::Options;
    use crate::sam::alignment::{io::Write as _, record::Flags};

    /// Unmapped records, 2 per slice, one of which is named with the test panic marker.
    fn marked_data() -> io::Result<Vec<u8>> {
        let header = sam::Header::default();
        let options = Options {
            records_per_slice: 2,
            ..Options::default()
        };

        let mut writer = cram::io::writer::Builder::default()
            .set_compression_options(options)
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;

        for i in 0..100 {
            let name = if i == 37 {
                PANIC_MARKER.to_vec()
            } else {
                format!("r{i}").into_bytes()
            };
            let record = RecordBuf::builder()
                .set_name(name)
                .set_flags(Flags::UNMAPPED)
                .set_sequence(b"ACGT".to_vec().into())
                .set_quality_scores(vec![30; 4].into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }

        writer.try_finish(&header)?;
        Ok(writer.into_inner())
    }

    #[test]
    fn a_panic_on_a_worker_is_raised_again_in_order() -> io::Result<()> {
        let data = marked_data()?;

        for workers in [1, 2, 4] {
            let workers = NonZero::new(workers).unwrap();
            let data = data.clone();

            within_timeout(move || {
                let mut reader = Reader::new(&data[..]);
                let header = reader.read_header().unwrap();
                let mut records = reader.multithreaded_records(&header, workers);

                // The records of the slices before the marked one come first.
                for i in 0..36 {
                    let record = records.next().unwrap().unwrap();
                    assert_eq!(record.name().unwrap(), format!("r{i}").as_bytes());
                }

                assert_panics_with(MARKER_PANIC, || records.next());
            });
        }

        Ok(())
    }
}
