//! Writing output records, formatted (VCF) or encoded (BCF) on worker threads with
//! `--threads`, and written in order.

use std::{collections::BTreeMap, num::NonZero, sync::Arc, thread};

use anyhow::bail;
use crossbeam_channel::{Receiver, Sender};
use pasta::{
    util::threads,
    variant::{
        Header, Record, bcf,
        io::{RecordSpan, WriteError, Writer},
        log::{self, Message},
        vcf,
    },
};

use crate::{
    output,
    workers::{self, ReportPanic},
};

/// Records per batch handed to a worker, at most.
const BATCH_RECORDS: usize = 1024;

/// Bytes of records ([`Record::allocated_bytes`]) shared among the batches in flight: each
/// batch is cut at [`BATCH_RECORDS`] or at its share of this, whichever comes first, so that
/// records of many samples do not hold memory in proportion to the sample count.
const IN_FLIGHT_BYTES: usize = 8 << 20;

/// The least bytes of records a batch is cut at, so that at many threads the hand-off of
/// each batch to a worker still costs little against its work.
const MIN_BATCH_BYTES: usize = 512 << 10;

/// The most records the workers have returned that are kept for reuse.
const MAX_SPARE_RECORDS: usize = 4 * BATCH_RECORDS;

/// The most bytes of records the workers have returned that are kept for reuse; records
/// past it are dropped.
const MAX_SPARE_BYTES: usize = 4 << 20;

/// What the workers turn records into.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Encoding {
    /// VCF text.
    Vcf,
    /// BCF records.
    Bcf,
}

struct Job {
    seq: u64,
    recs: Vec<Record>,
}

struct Done {
    seq: u64,
    recs: Vec<Record>,
    data: Vec<u8>,
    spans: Vec<RecordSpan>,
    messages: Vec<Message>,
    /// Whether a record could not be formatted (the batch stops there).
    failed: bool,
}

/// What a worker sends back.
enum Reply {
    Done(Done),
    /// The worker panicked: its batch never comes back.
    Panicked,
}

/// Formats or encodes a record onto `data`, returning its span (the write part of HTSlib's
/// `bcf_write`).
fn emit(
    header: &Header,
    rec: &mut Record,
    data: &mut Vec<u8>,
    encoding: Encoding,
) -> Option<RecordSpan> {
    bcf::check_sample_count(header, rec).ok()?;
    let before = data.len();
    let ok = match encoding {
        Encoding::Vcf => vcf::format_record(header, rec, data).is_ok(),
        Encoding::Bcf => bcf::encode_record(header, rec, data).is_ok(),
    };
    if !ok {
        data.truncate(before);
        return None;
    }
    Some(RecordSpan {
        len: (data.len() - before) as u32,
        rid: rec.rid(),
        pos: rec.pos(),
        rlen: rec.rlen(),
    })
}

fn worker(header: &Header, jobs: &Receiver<Job>, done: &Sender<Reply>, encoding: Encoding) {
    let _guard = ReportPanic::new(done, || Reply::Panicked);
    while let Ok(Job { seq, mut recs }) = jobs.recv() {
        let mut data = Vec::with_capacity(recs.len() * 256);
        let mut spans = Vec::with_capacity(recs.len());
        let (failed, messages) = log::capture(|| {
            for rec in &mut recs {
                match emit(header, rec, &mut data, encoding) {
                    Some(span) => spans.push(span),
                    None => return true,
                }
            }
            false
        });
        workers::panic_if_marked(&data);
        let d = Done {
            seq,
            recs,
            data,
            spans,
            messages,
            failed,
        };
        if done.send(Reply::Done(d)).is_err() {
            break;
        }
    }
}

/// The error of a failed write: a closed pipe passes through, anything else is the command's
/// message.
fn write_error(e: WriteError, message: &str) -> anyhow::Error {
    output::write_failure(e, || anyhow::anyhow!("{message}"))
}

/// Output records, written directly or through worker threads.
pub(crate) struct Emitter {
    pool: Option<Pool>,
    /// The message of a failed write (bcftools' `error()` after `bcf_write`).
    failed: String,
}

struct Pool {
    jobs: Option<Sender<Job>>,
    done: Receiver<Reply>,
    workers: Vec<thread::JoinHandle<()>>,
    batch: Vec<Record>,
    /// The bytes of the records in `batch`.
    batch_bytes: usize,
    /// The bytes of records a batch is cut at.
    max_batch_bytes: usize,
    next_seq: u64,
    next_write: u64,
    in_flight: usize,
    max_in_flight: usize,
    pending: BTreeMap<u64, Done>,
    /// Whether a record could not be written, which ends the output where a serial run ends
    /// it: nothing queued after it is written.
    stopped: bool,
    /// Records returned by the workers, for reuse.
    spare: Vec<Record>,
    /// The bytes of the records in `spare`.
    spare_bytes: usize,
}

impl Pool {
    /// Keeps records the workers are done with for reuse, up to [`MAX_SPARE_RECORDS`] and
    /// [`MAX_SPARE_BYTES`]; the rest are dropped.
    fn keep_spare(&mut self, recs: Vec<Record>) {
        for rec in recs {
            let bytes = rec.allocated_bytes();
            if self.spare.len() >= MAX_SPARE_RECORDS || self.spare_bytes + bytes > MAX_SPARE_BYTES {
                continue;
            }
            self.spare_bytes += bytes;
            self.spare.push(rec);
        }
    }

    /// Takes a record kept for reuse.
    fn take_spare(&mut self) -> Option<Record> {
        let rec = self.spare.pop()?;
        self.spare_bytes -= rec.allocated_bytes();
        Some(rec)
    }

    /// Waits for the next finished batch (`None` if the workers stopped). When a worker
    /// panicked, stops the others and re-raises its panic.
    fn recv(&mut self) -> Option<Done> {
        match self.done.recv() {
            Ok(Reply::Done(done)) => Some(done),
            Ok(Reply::Panicked) => {
                self.jobs = None;
                workers::rethrow(self.workers.drain(..), "a norm output worker panicked")
            }
            Err(_) => None,
        }
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        self.jobs = None;
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
    }
}

impl Emitter {
    /// Writes records directly, or with `threads` workers when it is above zero. A failed
    /// write reports `norm`'s message for the output `name`.
    pub(crate) fn new(header: &Header, writer: &Writer, threads: usize, name: String) -> Self {
        let failed = format!("[flush_buffer] Error: cannot write to {name}");
        if threads == 0 {
            return Self { pool: None, failed };
        }
        let header = Arc::new(header.clone());
        let encoding = if writer.output_type().is_vcf() {
            Encoding::Vcf
        } else {
            Encoding::Bcf
        };
        // Fewer workers when the system cannot start `threads`, and none (records written
        // directly) when it cannot start any.
        let n = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
        let started = threads::start(n, |count| {
            let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(count.get() * 2);
            let (done_tx, done_rx) = crossbeam_channel::unbounded::<Reply>();
            let workers = threads::spawn_all(count.get(), |_| {
                let (jobs, done, header) = (job_rx.clone(), done_tx.clone(), Arc::clone(&header));
                move || worker(&header, &jobs, &done, encoding)
            })?;
            Ok((workers, job_tx, done_rx))
        });
        let Ok((workers, job_tx, done_rx)) = started else {
            return Self { pool: None, failed };
        };
        let threads = workers.len();
        Self {
            pool: Some(Pool {
                jobs: Some(job_tx),
                done: done_rx,
                workers,
                batch: Vec::with_capacity(BATCH_RECORDS),
                batch_bytes: 0,
                max_batch_bytes: (IN_FLIGHT_BYTES / (threads * 2)).max(MIN_BATCH_BYTES),
                next_seq: 0,
                next_write: 0,
                in_flight: 0,
                max_in_flight: threads * 2,
                pending: BTreeMap::new(),
                stopped: false,
                spare: Vec::new(),
                spare_bytes: 0,
            }),
            failed,
        }
    }

    /// Reports `message` when a record cannot be written.
    #[must_use]
    pub(crate) fn with_write_error(mut self, message: String) -> Self {
        self.failed = message;
        self
    }

    /// The error for records that could not be written, such as a
    /// [`WriteError::Deferred`] from closing the output.
    pub(crate) fn write_failure(&self, e: WriteError) -> anyhow::Error {
        write_error(e, &self.failed)
    }

    /// Takes a record the workers are done with, if any.
    pub(crate) fn spare(&mut self) -> Option<Record> {
        self.pool.as_mut()?.take_spare()
    }

    /// Writes a record, or queues it for a worker. Returns the record when it was written
    /// directly, for reuse.
    pub(crate) fn write(
        &mut self,
        header: &Header,
        writer: &mut Writer,
        mut rec: Record,
    ) -> anyhow::Result<Option<Record>> {
        let Some(pool) = &mut self.pool else {
            if let Err(e) = writer.write_record(header, &mut rec) {
                return Err(write_error(e, &self.failed));
            }
            return Ok(Some(rec));
        };
        pool.batch_bytes += rec.allocated_bytes();
        pool.batch.push(rec);
        if pool.batch.len() >= BATCH_RECORDS || pool.batch_bytes >= pool.max_batch_bytes {
            self.send(header, writer)?;
        }
        Ok(None)
    }

    /// Sends the current batch to the workers, writing finished batches while too many are
    /// in flight.
    fn send(&mut self, header: &Header, writer: &mut Writer) -> anyhow::Result<()> {
        let full = |p: &Option<Pool>| p.as_ref().is_some_and(|p| p.in_flight >= p.max_in_flight);
        if self.pool.as_ref().is_none_or(|p| p.batch.is_empty()) {
            return Ok(());
        }
        while full(&self.pool) {
            self.receive_one(header, writer)?;
        }
        let pool = self.pool.as_mut().expect("a pool sends batches");
        // The next batch is likely the size of this one.
        let capacity = pool.batch.len();
        let recs = std::mem::replace(&mut pool.batch, Vec::with_capacity(capacity));
        pool.batch_bytes = 0;
        let job = Job {
            seq: pool.next_seq,
            recs,
        };
        pool.next_seq += 1;
        pool.in_flight += 1;
        if pool
            .jobs
            .as_ref()
            .expect("jobs are open until the end")
            .send(job)
            .is_err()
        {
            bail!("worker threads stopped");
        }
        Ok(())
    }

    /// Waits for one batch and writes the batches that are ready in order.
    fn receive_one(&mut self, header: &Header, writer: &mut Writer) -> anyhow::Result<()> {
        let pool = self.pool.as_mut().expect("a pool receives batches");
        let Some(done) = pool.recv() else {
            bail!("worker threads stopped");
        };
        pool.pending.insert(done.seq, done);
        self.write_ready(header, writer)
    }

    /// Writes the batches that are ready, in order.
    fn write_ready(&mut self, header: &Header, writer: &mut Writer) -> anyhow::Result<()> {
        let pool = self.pool.as_mut().expect("a pool writes batches");
        while let Some(done) = pool.pending.remove(&pool.next_write) {
            pool.next_write += 1;
            pool.in_flight -= 1;
            let Done {
                recs,
                data,
                spans,
                messages,
                failed,
                ..
            } = done;
            // One record at a time, as a serial run writes them, so that the output held in
            // the file buffer (and lost on a fatal error) is the same.
            let mut r = Ok(());
            let mut off = 0;
            for span in &spans {
                let end = off + span.len as usize;
                r = writer.write_spans(header, &data[off..end], std::slice::from_ref(span));
                if r.is_err() {
                    break;
                }
                off = end;
            }
            log::replay(messages);
            if let Err(e) = r {
                pool.stopped = true;
                return Err(write_error(e, &self.failed));
            }
            if failed {
                pool.stopped = true;
                bail!("{}", self.failed);
            }
            pool.keep_spare(recs);
        }
        Ok(())
    }

    /// Writes what was queued before a fatal error, as a serial run would have written it
    /// before stopping, so that the caller can drop the writer with the same output lost.
    /// Stops at the first batch that cannot be written, and writes nothing when the error
    /// was a record that could not be written: a serial run stops at that record.
    pub(crate) fn abort(&mut self, header: &Header, writer: &mut Writer) {
        if self.pool.as_ref().is_none_or(|p| p.stopped) || self.send(header, writer).is_err() {
            return;
        }
        loop {
            if self.write_ready(header, writer).is_err() {
                return;
            }
            let pool = self.pool.as_mut().expect("a pool receives batches");
            if pool.in_flight == 0 {
                return;
            }
            let Some(done) = pool.recv() else {
                return;
            };
            pool.pending.insert(done.seq, done);
        }
    }

    /// Writes everything queued.
    pub(crate) fn finish(&mut self, header: &Header, writer: &mut Writer) -> anyhow::Result<()> {
        if self.pool.is_none() {
            return Ok(());
        }
        self.send(header, writer)?;
        while self.pool.as_ref().is_some_and(|p| p.in_flight > 0) {
            self.receive_one(header, writer)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::{io, sync::Mutex};

    use pasta::variant::io::{OutputType, WriteOptions};

    use super::*;

    /// An output whose bytes the test reads back.
    #[derive(Clone, Default)]
    struct SharedSink(Arc<Mutex<Vec<u8>>>);

    impl io::Write for SharedSink {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.lock().expect("lock poisoned").extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn header() -> Header {
        Header::parse(
            b"##fileformat=VCFv4.2\n##contig=<ID=1>\n\
              ##FORMAT=<ID=X,Number=.,Type=Integer,Description=\"x\">\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n",
        )
        .expect("the header parses")
    }

    /// A record at `pos` whose FORMAT/X has `values` values, held in an edited field as
    /// `merge` builds records.
    fn record(header: &mut Header, pos: usize, values: usize) -> Record {
        let mut line = format!("1\t{pos}\t.\tA\tC\t.\t.\t.\tX\t0").into_bytes();
        let mut rec = Record::new();
        vcf::parse(header, &mut line, &mut rec).expect("the line parses");
        rec.update_format_i32(header, b"X", &vec![100_000; values])
            .expect("X is a FORMAT tag");
        rec
    }

    /// Writes `n` records of `values` values with `threads` workers, checking the batch and
    /// spare limits after each record; returns the BCF output.
    fn write_all(threads: usize, n: usize, values: usize) -> Vec<u8> {
        let mut header = header();
        let recs: Vec<Record> = (1..=n)
            .map(|pos| record(&mut header, pos, values))
            .collect();
        let sink = SharedSink::default();
        let options = WriteOptions {
            output_type: OutputType::BcfUncompressed,
            level: None,
            threads: 0,
        };
        let mut writer = Writer::from_sink(Box::new(sink.clone()), &options);
        writer.write_header(&header).expect("the header is written");
        let mut emitter = Emitter::new(&header, &writer, threads, "-".into());
        for rec in recs {
            emitter
                .write(&header, &mut writer, rec)
                .expect("the record is written");
            if let Some(pool) = &emitter.pool {
                assert!(pool.batch.len() < BATCH_RECORDS, "{}", pool.batch.len());
                assert!(
                    pool.batch_bytes < pool.max_batch_bytes,
                    "{} of {} bytes in the batch",
                    pool.batch_bytes,
                    pool.max_batch_bytes
                );
                let spare: usize = pool.spare.iter().map(Record::allocated_bytes).sum();
                assert_eq!(spare, pool.spare_bytes);
                assert!(spare <= MAX_SPARE_BYTES, "{spare} spare bytes");
            }
        }
        emitter
            .finish(&header, &mut writer)
            .expect("the records are written");
        writer.finish().expect("the output is finished");
        sink.0.lock().expect("lock poisoned").clone()
    }

    #[test]
    fn a_panicking_worker_ends_the_wait() {
        // One worker panics on its batch while the others wait for jobs, holding the results
        // channel open: `finish` must stop waiting and re-raise the panic.
        let panicked =
            crate::raw_input::within_seconds("norm output with a panicking worker", || {
                std::panic::catch_unwind(|| {
                    let mut header = header();
                    let mut line = format!("1\t1\t{}\tA\tC\t.\t.\t.\tX\t0", workers::PANIC_MARKER)
                        .into_bytes();
                    let mut rec = Record::new();
                    vcf::parse(&mut header, &mut line, &mut rec).expect("the line parses");
                    let options = WriteOptions {
                        output_type: OutputType::Vcf,
                        level: None,
                        threads: 0,
                    };
                    let mut writer = Writer::from_sink(Box::new(SharedSink::default()), &options);
                    let mut emitter = Emitter::new(&header, &writer, 3, "-".into());
                    let _ = emitter.write(&header, &mut writer, rec);
                    let _ = emitter.finish(&header, &mut writer);
                })
                .is_err()
            });
        assert!(panicked, "the worker's panic is re-raised");
    }

    #[test]
    fn batches_of_large_records_are_cut_by_bytes() {
        // About 800 kB per record: a batch takes a handful at one thread.
        let serial = write_all(0, 24, 200_000);
        let threaded = write_all(1, 24, 200_000);
        assert!(serial == threaded, "the threaded output differs");
        assert_eq!(write_all(4, 24, 200_000), serial);
    }

    #[test]
    fn batches_of_small_records_are_cut_by_count() {
        let mut header = header();
        let options = WriteOptions {
            output_type: OutputType::BcfUncompressed,
            level: None,
            threads: 0,
        };
        let mut writer = Writer::from_sink(Box::new(SharedSink::default()), &options);
        let mut emitter = Emitter::new(&header, &writer, 1, "-".into());
        for pos in 1..BATCH_RECORDS {
            let rec = record(&mut header, pos, 1);
            emitter
                .write(&header, &mut writer, rec)
                .expect("the record is written");
        }
        let pool = emitter.pool.as_ref().expect("a threaded emitter");
        assert_eq!(pool.batch.len(), BATCH_RECORDS - 1);
        let rec = record(&mut header, BATCH_RECORDS, 1);
        emitter
            .write(&header, &mut writer, rec)
            .expect("the record is written");
        let pool = emitter.pool.as_ref().expect("a threaded emitter");
        assert!(pool.batch.is_empty(), "{} records left", pool.batch.len());
        emitter
            .finish(&header, &mut writer)
            .expect("the records are written");
    }

    #[test]
    fn spare_records_are_kept_within_the_byte_limit() {
        let mut header = header();
        let options = WriteOptions {
            output_type: OutputType::BcfUncompressed,
            level: None,
            threads: 0,
        };
        let writer = Writer::from_sink(Box::new(SharedSink::default()), &options);
        let mut emitter = Emitter::new(&header, &writer, 1, "-".into());
        let pool = emitter.pool.as_mut().expect("a threaded emitter");
        // About 1.6 MB each: two fit.
        let recs = (1..=4)
            .map(|pos| record(&mut header, pos, 400_000))
            .collect();
        pool.keep_spare(recs);
        assert_eq!(pool.spare.len(), 2);
        let small = (1..=8).map(|pos| record(&mut header, pos, 1)).collect();
        pool.keep_spare(small);
        assert_eq!(pool.spare.len(), 10);
        while emitter.spare().is_some() {}
        let pool = emitter.pool.as_ref().expect("a threaded emitter");
        assert_eq!(pool.spare_bytes, 0);
    }

    #[test]
    fn a_record_that_cannot_be_written_ends_the_output_as_serially() {
        // Record `bad` is out of order, which the index finds as it is written; a serial run
        // stops there, so aborting writes none of the batches queued after it.
        let output = |threads: usize, bad: usize| {
            let mut header = header();
            let sink = SharedSink::default();
            let options = WriteOptions {
                output_type: OutputType::VcfBgzf,
                level: Some(0),
                threads: 0,
            };
            let mut writer = Writer::from_sink(Box::new(sink.clone()), &options);
            writer
                .request_index(pasta::variant::io::IndexRequest {
                    path: "never-written.csi".into(),
                    min_shift: 14,
                })
                .expect("BGZF output");
            writer.write_header(&header).expect("the header is written");
            let mut emitter = Emitter::new(&header, &writer, threads, "-".into());
            let written = (1..=60_000).try_for_each(|i| {
                let pos = if i == bad { 1 } else { i };
                let rec = record(&mut header, pos, 1);
                emitter.write(&header, &mut writer, rec).map(drop)
            });
            let result = written.and_then(|()| emitter.finish(&header, &mut writer));
            assert!(result.is_err(), "the unsorted record fails");
            emitter.abort(&header, &mut writer);
            writer.abandon();
            sink.0.lock().expect("lock poisoned").clone()
        };
        for bad in [20_000, 60_000] {
            let serial = output(0, bad);
            assert!(!serial.is_empty(), "blocks before record {bad} are written");
            for threads in [1, 3] {
                let threaded = output(threads, bad);
                assert!(
                    threaded == serial,
                    "record {bad}, threads {threads}: {} bytes, {} serially",
                    threaded.len(),
                    serial.len()
                );
            }
        }
    }
}
