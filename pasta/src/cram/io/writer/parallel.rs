//! Containers encoded and compressed on worker threads.
//!
//! A container goes through three steps (see [`EncodedContainer`]). Workers convert its records
//! and encode it; this thread then chooses its blocks' compression methods, in container order,
//! since each choice learns from the blocks before it; workers compress it; and this thread writes
//! the containers in order. The output is what encoding every container on this thread gives.
//!
//! Choosing usually only decides on a method. Every so often a block is a trial, compressed with
//! every candidate method before the next block is chosen for; those compressions go to the
//! workers ahead of other work, and this thread goes on taking records meanwhile.

use std::{
    collections::VecDeque,
    fmt,
    io::{self, Write},
    mem,
    num::NonZero,
    panic::{self, AssertUnwindSafe},
    sync::Arc,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Select, Sender, TryRecvError};

use super::{
    Batch, Context,
    compression::Compressor,
    container::{Choosing, EncodedContainer, encode_container, finish_container},
};
use crate::sam;
use crate::util::threads;

/// Containers in flight per worker: submitted and not yet written, whether queued, being encoded
/// or compressed, or done and waiting for an earlier container. This bounds memory while keeping
/// the workers busy.
const IN_FLIGHT_PER_WORKER: usize = 2;

/// A job's result, or its panic.
type Outcome<T> = thread::Result<io::Result<T>>;

/// A container in flight, oldest first.
enum InFlight {
    /// Being encoded.
    Encoding(Receiver<Outcome<EncodedContainer>>),
    /// Its compression being chosen, waiting for trials ([`Workers::trying`]).
    Trying,
    /// Its compression chosen, being compressed.
    Finishing(Receiver<Outcome<Vec<u8>>>),
}

/// A container whose compression is being chosen, with the trials it waits for.
struct Trying {
    choosing: Choosing,
    jobs: Vec<Receiver<Outcome<Option<Vec<u8>>>>>,
    /// The results of the first jobs.
    results: Vec<Option<Vec<u8>>>,
}

impl Trying {
    /// Takes the results that are in, in order, and returns whether all are.
    fn poll(&mut self) -> bool {
        while let Some(rx) = self.jobs.get(self.results.len()) {
            match rx.try_recv() {
                Ok(outcome) => {
                    let result = unwrap(outcome).expect("a trial job does not fail");
                    self.results.push(result);
                }
                Err(TryRecvError::Empty) => return false,
                Err(TryRecvError::Disconnected) => unreachable!("every job sends its result"),
            }
        }

        true
    }

    /// The job whose result is waited for next, if any.
    fn next_job(&self) -> Option<&Receiver<Outcome<Option<Vec<u8>>>>> {
        self.jobs.get(self.results.len())
    }
}

/// The workers of a CRAM writer and the containers they work on.
pub(super) struct Workers {
    pool: Pool,
    /// Containers in flight. The first `chosen` are [`InFlight::Finishing`], the next may be
    /// [`InFlight::Trying`], and the rest are [`InFlight::Encoding`].
    in_flight: VecDeque<InFlight>,
    chosen: usize,
    trying: Option<Trying>,
    /// The writer's context and header as the workers see them, shared by the containers.
    context: Option<Arc<Context>>,
    header: Option<Arc<sam::Header>>,
}

impl fmt::Debug for Workers {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Workers")
            .field("workers", &self.pool.threads.len())
            .field("in_flight", &self.in_flight.len())
            .finish_non_exhaustive()
    }
}

impl Workers {
    /// Starts `n` workers, or fewer when the system cannot start that many.
    ///
    /// # Errors
    ///
    /// Fails when not even one worker can be started.
    pub(super) fn start(n: NonZero<usize>) -> io::Result<Self> {
        Ok(Self {
            pool: Pool::start(n)?,
            in_flight: VecDeque::new(),
            chosen: 0,
            trying: None,
            context: None,
            header: None,
        })
    }

    /// The number of workers.
    pub(super) fn thread_count(&self) -> usize {
        self.pool.threads.len()
    }

    /// The most containers in flight between calls.
    pub(super) fn capacity(&self) -> usize {
        self.pool.threads.len() * IN_FLIGHT_PER_WORKER
    }

    /// Hands a container of `records` (not empty), the first of which is record `record_counter`
    /// of the stream, to the workers.
    ///
    /// The workers keep the context as it is on the first call, and a copy of `header` until a
    /// call with a different one.
    pub(super) fn submit(
        &mut self,
        ctx: &Context,
        header: &sam::Header,
        record_counter: u64,
        batch: Batch,
    ) {
        let ctx = Arc::clone(self.context.get_or_insert_with(|| Arc::new(ctx.clone())));

        let header = match &self.header {
            Some(h) if **h == *header => Arc::clone(h),
            _ => Arc::clone(self.header.insert(Arc::new(header.clone()))),
        };

        let rx = self.pool.spawn(Priority::Normal, move || {
            let mut batch = batch;
            // Encoding uses the reference sequences, which may be released from the cache.
            let reference_sequences = mem::take(&mut batch.reference_sequences);
            let mut records = super::convert(&header, batch, &reference_sequences)?;
            encode_container(&ctx, &header, record_counter, &mut records)
        });

        self.in_flight.push_back(InFlight::Encoding(rx));
    }

    /// Chooses the compression of the containers encoded so far and writes the ones compressed
    /// so far, in order, waiting until at most `max` containers are in flight.
    ///
    /// # Errors
    ///
    /// The first error of a container, in order, or of writing to `writer`. The containers before
    /// the one that failed are written first.
    ///
    /// # Panics
    ///
    /// A panic on a worker is resumed here.
    pub(super) fn advance<W>(
        &mut self,
        writer: &mut W,
        mut compressor: Option<&mut Compressor>,
        max: usize,
    ) -> io::Result<()>
    where
        W: Write,
    {
        loop {
            if let Err(e) = self.choose(compressor.as_deref_mut()) {
                self.write_chosen(writer)?;
                return Err(e);
            }

            while let Some(InFlight::Finishing(rx)) = self.in_flight.front() {
                match rx.try_recv() {
                    Ok(outcome) => {
                        self.in_flight.pop_front();
                        self.chosen -= 1;
                        writer.write_all(&unwrap(outcome)?)?;
                    }
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => unreachable!("every job sends its result"),
                }
            }

            if self.in_flight.len() <= max {
                return Ok(());
            }

            // Wait for the oldest container to be compressed or for what choosing waits for,
            // whichever comes first.
            let mut select = Select::new();

            if let Some(InFlight::Finishing(rx)) = self.in_flight.front() {
                select.recv(rx);
            }

            match self.in_flight.get(self.chosen) {
                Some(InFlight::Encoding(rx)) => {
                    select.recv(rx);
                }
                Some(InFlight::Trying) => {
                    let trying = self.trying.as_ref().expect("a container is being tried");

                    if let Some(rx) = trying.next_job() {
                        select.recv(rx);
                    }
                }
                _ => {}
            }

            select.ready();
        }
    }

    /// Waits for the containers whose compression is chosen and writes them, in order.
    fn write_chosen<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: Write,
    {
        while let Some(InFlight::Finishing(rx)) = self.in_flight.front() {
            let outcome = rx.recv().expect("every job sends its result");
            self.in_flight.pop_front();
            self.chosen -= 1;
            writer.write_all(&unwrap(outcome)?)?;
        }

        Ok(())
    }

    /// Chooses the compression of the containers in order, as far as it can without waiting, and
    /// hands the ones chosen for to the workers to compress.
    fn choose(&mut self, mut compressor: Option<&mut Compressor>) -> io::Result<()> {
        loop {
            let mut trying = match self.in_flight.get(self.chosen) {
                Some(InFlight::Encoding(rx)) => match rx.try_recv() {
                    Ok(outcome) => Trying {
                        choosing: Choosing::new(unwrap(outcome)?),
                        jobs: Vec::new(),
                        results: Vec::new(),
                    },
                    Err(TryRecvError::Empty) => return Ok(()),
                    Err(TryRecvError::Disconnected) => unreachable!("every job sends its result"),
                },
                Some(InFlight::Trying) => {
                    let trying = self.trying.as_mut().expect("a container is being tried");

                    if !trying.poll() {
                        return Ok(());
                    }

                    self.trying.take().expect("a container is being tried")
                }
                _ => return Ok(()),
            };

            if let Some(compressor) = compressor.as_deref_mut() {
                // Trials end even when none of their blocks had a method to try, and so no jobs.
                let results = mem::take(&mut trying.results);
                trying.choosing.end_trials(compressor, results);
                trying.jobs.clear();

                if let Some(jobs) = trying.choosing.step(compressor) {
                    // Ahead of other work: the next containers wait for these.
                    trying.jobs = jobs
                        .into_iter()
                        .map(|job| self.pool.spawn(Priority::Urgent, move || Ok(job.run())))
                        .collect();

                    self.trying = Some(trying);
                    self.in_flight[self.chosen] = InFlight::Trying;
                    continue;
                }
            }

            let container = trying.choosing.into_container();

            let rx = self
                .pool
                .spawn(Priority::Normal, move || finish_container(container));

            self.in_flight[self.chosen] = InFlight::Finishing(rx);
            self.chosen += 1;
        }
    }
}

/// Returns a job's result, resuming its panic.
fn unwrap<T>(outcome: Outcome<T>) -> io::Result<T> {
    outcome.unwrap_or_else(|payload| panic::resume_unwind(payload))
}

/// A job for a worker.
type Job = Box<dyn FnOnce() + Send>;

#[derive(Clone, Copy)]
enum Priority {
    Normal,
    /// Run before every normal job not yet started.
    Urgent,
}

/// Worker threads that run jobs, urgent ones first. Dropping the pool waits for the job each
/// worker is running, if any, and drops the rest.
struct Pool {
    queues: Option<(Sender<Job>, Sender<Job>)>,
    threads: Vec<JoinHandle<()>>,
}

impl Pool {
    fn start(n: NonZero<usize>) -> io::Result<Self> {
        threads::start(n, |count| {
            let (urgent_tx, urgent_rx) = crossbeam_channel::unbounded::<Job>();
            let (normal_tx, normal_rx) = crossbeam_channel::unbounded::<Job>();

            // On a shortfall, the queues are dropped here, which stops the threads started.
            let threads = threads::spawn_all(count.get(), |_| {
                let (urgent, normal) = (urgent_rx.clone(), normal_rx.clone());
                move || work(&urgent, &normal)
            })?;

            Ok(Self {
                queues: Some((normal_tx, urgent_tx)),
                threads,
            })
        })
    }

    /// Runs `f` on a worker, returning a channel that receives its result or its panic.
    fn spawn<T, F>(&self, priority: Priority, f: F) -> Receiver<Outcome<T>>
    where
        T: Send + 'static,
        F: FnOnce() -> io::Result<T> + Send + 'static,
    {
        let (tx, rx) = crossbeam_channel::bounded(1);

        let job: Job = Box::new(move || {
            // The receiver may be gone: the writer was dropped.
            let _ = tx.send(panic::catch_unwind(AssertUnwindSafe(f)));
        });

        let (normal, urgent) = self.queues.as_ref().expect("queues exist until drop");

        let queue = match priority {
            Priority::Normal => normal,
            Priority::Urgent => urgent,
        };

        queue
            .send(job)
            .expect("workers run until the pool is dropped");

        rx
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        drop(self.queues.take());

        for handle in mem::take(&mut self.threads) {
            let _ = handle.join();
        }
    }
}

/// Runs jobs, urgent ones first, until the queues close.
fn work(urgent: &Receiver<Job>, normal: &Receiver<Job>) {
    loop {
        let job = crossbeam_channel::select_biased! {
            recv(urgent) -> job => job,
            recv(normal) -> job => job,
        };

        match job {
            Ok(job) => job(),
            Err(_) => return,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::bam;
    use crate::core::Position;
    use crate::cram::{
        container::BlockContentEncoderMap,
        file_definition::Version,
        io::writer::{
            Builder, Writer,
            compression::{Options, Profile},
        },
    };
    use crate::fasta;
    use crate::sam::alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::data::field::Value,
    };

    const SQ0: &[u8] = b"ACGTTGCAAGGCTTACCGATTGACCATGGTACCGTA";

    fn repository() -> fasta::Repository {
        // sq2 is missing, so its records are stored verbatim against an embedded reference.
        fasta::Repository::new(vec![
            fasta::Record::new(
                fasta::record::Definition::new("sq0", None),
                fasta::record::Sequence::from(SQ0.repeat(30)),
            ),
            fasta::Record::new(
                fasta::record::Definition::new("sq1", None),
                fasta::record::Sequence::from(b"TTGCAGGA".repeat(100)),
            ),
        ])
    }

    fn header() -> sam::Header {
        "@SQ\tSN:sq0\tLN:1080\n@SQ\tSN:sq1\tLN:800\n@SQ\tSN:sq2\tLN:500\n@RG\tID:rg0\n@RG\tID:rg1\n"
            .parse()
            .unwrap()
    }

    /// Records of every kind: mapped on each sequence (in pairs, with mismatches, clips, and
    /// indels), unmapped, without qualities or bases, and with tags.
    fn records(n: usize) -> Vec<RecordBuf> {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = move |m: u64| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            (state >> 33) % m
        };

        (0..n)
            .map(|i| {
                let len = 20 + next(40) as usize;
                let mut sequence: Vec<u8> = (0..len).map(|_| b"ACGT"[next(4) as usize]).collect();
                let quality_scores: Vec<u8> = (0..len).map(|_| 10 + next(30) as u8).collect();

                let id = (i / (n / 3 + 1)).min(2);
                let start = 1 + (i % 40) * 10 + next(5) as usize;

                let cigar = match next(4) {
                    0 => vec![Op::new(Kind::Match, len)],
                    1 => vec![Op::new(Kind::SoftClip, 3), Op::new(Kind::Match, len - 3)],
                    2 => vec![
                        Op::new(Kind::Match, 8),
                        Op::new(Kind::Insertion, 2),
                        Op::new(Kind::Match, len - 10),
                    ],
                    _ => vec![
                        Op::new(Kind::Match, 8),
                        Op::new(Kind::Deletion, 3),
                        Op::new(Kind::Match, len - 8),
                    ],
                };

                // Mostly matching the reference, so that substitutions are few.
                if id == 0 && next(2) == 0 {
                    let offset = (start - 1) % SQ0.len();
                    let reference = SQ0.repeat(3);
                    sequence.copy_from_slice(&reference[offset..offset + len]);
                    sequence[next(len as u64) as usize] = b'N';
                }

                let mut builder = RecordBuf::builder()
                    .set_name(format!("read{}", i / 2))
                    .set_sequence(sequence.into())
                    .set_data(
                        [
                            (Tag::READ_GROUP, Value::from(format!("rg{}", i % 2))),
                            (Tag::ALIGNMENT_HIT_COUNT, Value::from(next(3) as u8)),
                        ]
                        .into_iter()
                        .collect(),
                    );

                if i % 17 != 0 {
                    builder = builder.set_quality_scores(quality_scores.into());
                }

                if i % 23 == 5 {
                    builder.set_flags(Flags::UNMAPPED).build()
                } else {
                    let mut flags = Flags::SEGMENTED;

                    if i % 2 == 1 {
                        flags |= Flags::LAST_SEGMENT;
                    } else {
                        flags |= Flags::FIRST_SEGMENT;
                    }

                    builder
                        .set_flags(flags)
                        .set_reference_sequence_id(id)
                        .set_alignment_start(Position::new(start).unwrap())
                        .set_cigar(cigar.into_iter().collect())
                        .set_mate_reference_sequence_id(id)
                        .set_mate_alignment_start(Position::new(start).unwrap())
                        .build()
                }
            })
            .collect()
    }

    fn to_bam(header: &sam::Header, record: &RecordBuf) -> io::Result<bam::Record> {
        let mut encoder = bam::io::Writer::from(Vec::new());
        encoder.write_alignment_record(header, record)?;
        let buf = encoder.into_inner();
        bam::Record::try_from(buf[size_of::<u32>()..].to_vec())
    }

    /// A configuration of the writer: how blocks are compressed, and slices per container.
    #[derive(Clone, Copy, Debug)]
    enum Config {
        Choose(Options, usize),
        Fixed,
    }

    fn build<W: Write>(config: Config, workers: usize, inner: W) -> Writer<W> {
        let mut builder = Builder::default().set_reference_sequence_repository(repository());

        builder = match config {
            Config::Choose(options, _) => builder.set_compression_options(options),
            Config::Fixed => {
                builder.set_block_content_encoder_map(BlockContentEncoderMap::default())
            }
        };

        if let Some(n) = NonZero::new(workers) {
            builder = builder.set_worker_count(n);
        }

        let mut writer = builder.build_from_writer(inner);

        match config {
            Config::Choose(options, slices) => {
                writer.context.records_per_slice = options.records_per_slice / slices;
            }
            Config::Fixed => {
                writer.context.records_per_slice = 40;
                writer.records_per_container = 40;
            }
        }

        writer
    }

    /// Writes `records` with `workers`, as BAM records but every third as a record buffer.
    fn write(config: Config, workers: usize, records: &[RecordBuf]) -> io::Result<Vec<u8>> {
        let header = header();
        let mut writer = build(config, workers, Vec::new());
        writer.write_header(&header)?;

        for (i, record) in records.iter().enumerate() {
            if i % 3 == 0 {
                writer.write_alignment_record(&header, record)?;
            } else {
                writer.write_bam_record(&header, &to_bam(&header, record)?)?;
            }

            if let Some(workers) = &writer.workers {
                assert!(workers.in_flight.len() <= workers.capacity());
            }
        }

        writer.try_finish(&header)?;
        Ok(writer.into_inner())
    }

    fn options(version: Version, profile: Profile, records_per_slice: usize) -> Options {
        let mut options = Options::new(version, profile);
        options.records_per_slice = records_per_slice;
        options.bases_per_slice = records_per_slice * 35;
        options
    }

    #[test]
    fn test_workers_write_what_this_thread_writes() -> io::Result<()> {
        let records = records(3000);

        let configs = [
            // Over 100 containers, for several rounds of trials.
            Config::Choose(options(Version::new(3, 1), Profile::Normal, 25), 1),
            Config::Choose(options(Version::new(3, 0), Profile::Fast, 40), 1),
            // Several slices per container.
            Config::Choose(options(Version::new(3, 1), Profile::Archive, 120), 3),
            // Gzip levels 2 to 4 without rANS leave some trials with no method to try.
            Config::Choose(
                Options {
                    level: 3,
                    use_rans: false,
                    ..options(Version::new(3, 1), Profile::Normal, 25)
                },
                1,
            ),
            Config::Fixed,
        ];

        for config in configs {
            let expected = write(config, 0, &records)?;

            let mut reader = crate::cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository())
                .build_from_reader(Cursor::new(&expected));
            let header = reader.read_header()?;
            let names = reader
                .records(&header)
                .map(|record| record.map(|record| record.name().map(|name| name.to_vec())))
                .collect::<io::Result<Vec<_>>>()?;
            let expected_names: Vec<_> = records
                .iter()
                .map(|record| record.name().map(|name| name.to_vec()))
                .collect();
            assert_eq!(names, expected_names, "{config:?}");

            for workers in [1, 2, 3, 8] {
                let actual = write(config, workers, &records)?;
                assert!(actual == expected, "{config:?} with {workers} workers");
            }
        }

        Ok(())
    }

    #[test]
    fn test_workers_stop_at_the_first_error() -> io::Result<()> {
        let header = header();
        let records = records(400);
        let config = Config::Choose(options(Version::new(3, 1), Profile::Normal, 25), 1);

        // Record 230 names a reference sequence the header lacks.
        let mut bad = to_bam(&header, &records[230])?;
        bad.0[..4].copy_from_slice(&7i32.to_le_bytes());

        for workers in [0, 1, 4] {
            let mut writer = build(config, workers, Vec::new());
            writer.write_header(&header)?;

            let mut result = Ok(());

            for (i, record) in records.iter().enumerate() {
                let record = if i == 230 {
                    bad.clone()
                } else {
                    to_bam(&header, record)?
                };

                result = writer.write_bam_record(&header, &record);

                if result.is_err() {
                    break;
                }
            }

            let result = result.and_then(|()| writer.try_finish(&header));
            let e = result.expect_err("a record fails");
            assert_eq!(
                e.kind(),
                io::ErrorKind::InvalidData,
                "{workers} workers: {e}"
            );

            // On the workers, the writer fails from then on.
            if workers > 0 {
                let record = to_bam(&header, &records[0])?;
                assert!(writer.write_bam_record(&header, &record).is_err());
                assert!(writer.try_finish(&header).is_err());
            }

            // The containers before the one that failed are written, and none after.
            let expected = {
                let mut writer = build(config, 0, Vec::new());
                writer.write_header(&header)?;

                for record in &records[..230] {
                    writer.write_bam_record(&header, &to_bam(&header, record)?)?;
                }

                writer.into_inner()
            };

            assert!(writer.get_ref() == &expected, "{workers} workers");
        }

        Ok(())
    }

    /// Fails every write after `len` bytes.
    struct FailingWriter {
        len: usize,
    }

    impl Write for FailingWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if buf.len() > self.len {
                return Err(io::Error::other("disk full"));
            }

            self.len -= buf.len();
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn test_workers_report_a_write_error() -> io::Result<()> {
        let header = header();
        let records = records(2000);
        let config = Config::Choose(options(Version::new(3, 1), Profile::Normal, 25), 1);

        for workers in [0, 1, 4] {
            let mut writer = build(config, workers, FailingWriter { len: 20_000 });
            writer.write_header(&header)?;

            let mut result = Ok(());

            for record in &records {
                result = writer.write_bam_record(&header, &to_bam(&header, record)?);

                if result.is_err() {
                    break;
                }
            }

            let e = result
                .and_then(|()| writer.try_finish(&header))
                .expect_err("the stream fails");
            assert_eq!(e.to_string(), "disk full", "{workers} workers");
        }

        Ok(())
    }
}
