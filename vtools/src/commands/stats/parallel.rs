//! Collecting statistics on worker threads (`--threads`).
//!
//! Every worker has its own [`Collector`]; they are merged at the end, which gives the serial
//! result because every statistic is a sum. Records are handed out in batches:
//!
//! - When the input is read straight through (no regions, targets, or `-f`), this thread only
//!   reads raw data: VCF lines or BCF records. Workers parse and check them against a snapshot
//!   of the header. A VCF line that needs a header change (an undefined tag or contig) stops
//!   its batch there; this thread parses the rest of that batch itself, in order, updating the
//!   header as HTSlib does, and later batches get the new header.
//! - Otherwise this thread reads records through the [`SyncedReader`] and workers collect them.
//!
//! Messages are captured per record and printed in input order, so stderr matches a serial
//! run; an error stops at the same record as a serial run.

use std::{collections::BTreeMap, num::NonZero, panic, sync::Arc, thread};

use crossbeam_channel::{Receiver, Sender};
use pasta::{
    util::threads,
    variant::{
        Header, Record, bcf,
        io::InputFormat,
        log::{self, Message},
        synced::{SyncedError, SyncedReader},
        vcf::{self, ParseError},
    },
};

use super::{
    Filter, Setup, collect::Collector, collect::Fatal, compile_filter, indel_context::IndelContext,
    process, read_error,
};
use crate::{
    pipeline::{LastLineCr, next_line, read_vcf_batch},
    workers::{self, ReportPanic},
};

/// Records per batch when this thread reads records, at most.
const BATCH_RECORDS: usize = 256;
/// Bytes per batch of raw input, and of records ([`Record::allocated_bytes`]) when this
/// thread reads records: a batch of records is cut at [`BATCH_RECORDS`] or at this many
/// bytes, whichever comes first, so that records of many samples do not hold memory in
/// proportion to the sample count.
const BATCH_BYTES: usize = 1 << 20;

/// Whether a batch of `len` records holding `bytes` bytes is full.
fn records_batch_full(len: usize, bytes: usize) -> bool {
    len >= BATCH_RECORDS || bytes >= BATCH_BYTES
}

/// What a worker processes.
enum Input {
    /// Records read by this thread; the first `len` are used.
    Records(Vec<Record>, usize),
    /// Whole VCF lines.
    Vcf(Vec<u8>),
    /// Whole BCF records.
    Bcf(Vec<u8>),
}

struct Job {
    seq: u64,
    input: Input,
    header: Arc<Header>,
}

/// Why a raw batch stopped before its end.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Stop {
    /// The VCF line needs a header update: this thread continues from it.
    HeaderUpdate,
    /// The VCF line is malformed.
    Parse,
    /// The BCF record is malformed.
    Bcf,
}

/// Where and why a batch stopped.
enum Halt {
    /// Collecting the record failed.
    Fatal(Fatal),
    /// Reading the record failed (raw batches), at this byte offset.
    Stop(Stop, usize),
}

struct Done {
    seq: u64,
    input: Input,
    /// The messages of each record, by index in the batch.
    messages: Vec<(usize, Vec<Message>)>,
    /// The record where the batch stopped.
    halt: Option<(usize, Halt)>,
}

/// What a worker sends back.
enum Reply {
    Done(Done),
    /// The worker panicked: its batch never comes back.
    Panicked,
}

/// The batches handed out and not yet replayed. Results that arrive before their turn wait
/// here and still count, so the window bounds the results held as well as the batches being
/// collected, however long an earlier batch takes.
///
/// The tests check this type, not [`Driver::drive`], so `drive` must make every window check
/// (when to hand out, wait, or stop) through it and keep no counts of its own.
struct Window<T> {
    size: usize,
    next_seq: u64,
    next_replay: u64,
    waiting: BTreeMap<u64, T>,
}

impl<T> Window<T> {
    fn new(size: usize) -> Self {
        Self {
            size,
            next_seq: 0,
            next_replay: 0,
            waiting: BTreeMap::new(),
        }
    }

    /// Whether no more batches may be handed out until one is replayed.
    fn is_full(&self) -> bool {
        self.next_seq - self.next_replay >= self.size as u64
    }

    /// Whether every batch handed out was replayed. Otherwise, once [`Self::next_ready`]
    /// returns `None`, the result of the next batch in order has not arrived yet, so waiting
    /// for a result cannot block forever.
    fn is_empty(&self) -> bool {
        self.next_seq == self.next_replay
    }

    /// Takes the sequence number of the next batch handed out.
    fn send(&mut self) -> u64 {
        let seq = self.next_seq;
        self.next_seq += 1;
        seq
    }

    /// Keeps the result of batch `seq` until its turn.
    fn receive(&mut self, seq: u64, done: T) {
        self.waiting.insert(seq, done);
    }

    /// Takes the next batch in order, if its result has arrived.
    fn next_ready(&mut self) -> Option<(u64, T)> {
        let done = self.waiting.remove(&self.next_replay)?;
        let seq = self.next_replay;
        self.next_replay += 1;
        Some((seq, done))
    }
}

/// Prints the messages of reading and of collecting a batch in record order, up to record
/// `upto` (inclusive) when given.
fn replay(read: Vec<(usize, Vec<Message>)>, work: Vec<(usize, Vec<Message>)>, upto: Option<usize>) {
    let mut work = work.into_iter().peekable();
    for (i, m) in read {
        if upto.is_some_and(|u| i > u) {
            break;
        }
        while let Some((_, w)) = work.next_if(|(j, _)| *j < i) {
            log::replay(w);
        }
        log::replay(m);
    }
    for (j, w) in work {
        if upto.is_some_and(|u| j > u) {
            break;
        }
        log::replay(w);
    }
}

/// What identifies a header version: header lines are only ever added while reading.
fn header_version(h: &Header) -> (usize, usize, usize) {
    (h.id_count(), h.contig_count(), h.record_count())
}

/// A worker's state.
struct Worker<'a> {
    setup: &'a Setup,
    collector: Collector,
    filter: Option<Filter>,
    rec: Record,
    line: Vec<u8>,
    last_line_cr: LastLineCr,
}

impl Worker<'_> {
    /// Processes a batch, stopping at the first record that fails.
    fn run(&mut self, job: Job) -> Done {
        let header = &*job.header;
        let cfg = &self.setup.cfg;
        let mut messages = Vec::new();
        let mut halt = None;
        let mut input = job.input;
        match &mut input {
            Input::Records(records, len) => {
                for (i, rec) in records[..*len].iter_mut().enumerate() {
                    let (r, m) = log::capture(|| {
                        process(&mut self.collector, cfg, &mut self.filter, header, rec)
                    });
                    if !m.is_empty() {
                        messages.push((i, m));
                    }
                    if let Err(e) = r {
                        halt = Some((i, Halt::Fatal(e)));
                        break;
                    }
                }
            }
            Input::Vcf(buf) => {
                workers::panic_if_marked(buf);
                let mut pos = 0;
                let mut i = 0;
                loop {
                    let start = pos;
                    let Some(l) = next_line(buf, &mut pos, self.last_line_cr) else {
                        break;
                    };
                    self.line.clear();
                    self.line.extend_from_slice(l);
                    let (r, m) = log::capture(|| {
                        match vcf::parse_shared(header, &mut self.line, &mut self.rec) {
                            Ok(()) => {}
                            Err(ParseError::HeaderUpdateNeeded) => {
                                return Err(Halt::Stop(Stop::HeaderUpdate, start));
                            }
                            Err(_) => return Err(Halt::Stop(Stop::Parse, start)),
                        }
                        process(
                            &mut self.collector,
                            cfg,
                            &mut self.filter,
                            header,
                            &mut self.rec,
                        )
                        .map_err(Halt::Fatal)
                    });
                    // The messages of a line handed back to this thread are printed when it
                    // parses the line again.
                    let handed_back = matches!(r, Err(Halt::Stop(Stop::HeaderUpdate, _)));
                    if !m.is_empty() && !handed_back {
                        messages.push((i, m));
                    }
                    if let Err(h) = r {
                        halt = Some((i, h));
                        break;
                    }
                    i += 1;
                }
            }
            Input::Bcf(buf) => {
                workers::panic_if_marked(buf);
                let mut pos = 0;
                let mut i = 0;
                while pos < buf.len() {
                    let start = pos;
                    let mut slice = &buf[pos..];
                    let before = slice.len();
                    let (r, m) = log::capture(|| {
                        if !matches!(
                            bcf::read_record_raw(&mut slice, &mut self.rec),
                            Ok(bcf::ReadOutcome::Record)
                        ) || bcf::check_record(Some(header), &mut self.rec).is_err()
                        {
                            return Err(Halt::Stop(Stop::Bcf, start));
                        }
                        process(
                            &mut self.collector,
                            cfg,
                            &mut self.filter,
                            header,
                            &mut self.rec,
                        )
                        .map_err(Halt::Fatal)
                    });
                    pos += before - slice.len();
                    if !m.is_empty() {
                        messages.push((i, m));
                    }
                    if let Err(h) = r {
                        halt = Some((i, h));
                        break;
                    }
                    i += 1;
                }
            }
        }
        Done {
            seq: job.seq,
            input,
            messages,
            halt,
        }
    }
}

/// Collects the statistics with `threads` workers, or fewer when the threads or the workers'
/// tables (`-d`, `-u`) cannot all be had. `None` when no worker thread starts or no worker gets
/// its tables: the caller then collects on its own thread, as without threads.
pub(super) fn collect(
    reader: &mut SyncedReader,
    setup: &Setup,
    indel_ctx: Option<&IndelContext>,
    threads: usize,
) -> Result<Option<Collector>, Fatal> {
    let cfg = &setup.cfg;
    let max_unpack = setup.max_unpack();
    let last_line_cr = LastLineCr::of(reader.reader_mut().compression());
    let tables = || {
        Collector::try_new(
            cfg,
            &setup.user_types,
            indel_ctx.and_then(IndelContext::reopen),
            None,
        )
    };
    let new_worker = |collector: Collector, header: &Header| -> Result<Worker<'_>, Fatal> {
        let mut rec = Record::new();
        rec.set_max_unpack(max_unpack);
        Ok(Worker {
            setup,
            collector,
            filter: compile_filter(setup, header)?,
            rec,
            line: Vec::new(),
            last_line_cr,
        })
    };
    // Records this thread parses itself, after a header change. Its tables are the ones a run
    // without threads needs, so failing to allocate them is the error.
    let mut own = new_worker(tables()?, reader.header())?;
    let raw = reader.is_plain_stream();
    // The closure owns `job_tx` (it drops it below), so a panic on this thread drops it while
    // unwinding and the workers' `recv` ends; their sends cannot block either, since the window
    // never has more batches out than `done_tx` holds. So the scope's join cannot wait forever.
    thread::scope(|scope| {
        // The threads start first, fewer of them when the system cannot start `threads`, and
        // then get their workers. The command checked that a thread starts, but other threads
        // (the reader's) may have taken the room since.
        let Ok(idle) = threads::start_idle_scoped(
            scope,
            NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN),
        ) else {
            return Ok(None);
        };
        let mut workers = Vec::with_capacity(idle.len());
        for _ in 0..idle.len() {
            let Ok(collector) = tables() else {
                // Half of the workers that got their tables go on, which leaves room for the
                // batches in flight.
                let keep = workers.len() / 2;
                workers.truncate(keep);
                break;
            };
            workers.push(new_worker(collector, reader.header())?);
        }
        if workers.is_empty() {
            return Ok(None);
        }
        let threads = workers.len();
        let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(threads * 2);
        let (done_tx, done_rx) = crossbeam_channel::bounded::<Reply>(threads * 2);
        let handles: Vec<_> = workers
            .into_iter()
            .zip(idle)
            .map(|(mut w, t)| {
                let job_rx = job_rx.clone();
                let done_tx = done_tx.clone();
                t.run(move || {
                    let _guard = ReportPanic::new(&done_tx, || Reply::Panicked);
                    while let Ok(job) = job_rx.recv() {
                        if done_tx.send(Reply::Done(w.run(job))).is_err() {
                            break;
                        }
                    }
                    w.collector
                })
            })
            .collect();
        drop(done_tx);
        drop(job_rx);
        let mut driver = Driver {
            reader,
            job_tx: &job_tx,
            done_rx: &done_rx,
            threads,
            own: &mut own,
            raw,
        };
        let result = driver.drive();
        drop(job_tx);
        // Drain what workers still send so they can finish.
        while done_rx.recv().is_ok() {}
        // A worker's panic is raised here, before any result.
        let collectors: Vec<Collector> = handles
            .into_iter()
            .map(|h| h.join().unwrap_or_else(|p| panic::resume_unwind(p)))
            .collect();
        result?;
        let mut merged = own.collector;
        for c in &collectors {
            merged.merge(c);
        }
        Ok(Some(merged))
    })
}

/// This thread's part: reading batches, handing them out, and replaying results in order.
struct Driver<'a, 'w> {
    reader: &'a mut SyncedReader,
    job_tx: &'a Sender<Job>,
    done_rx: &'a Receiver<Reply>,
    threads: usize,
    own: &'a mut Worker<'w>,
    raw: bool,
}

impl Driver<'_, '_> {
    fn drive(&mut self) -> Result<(), Fatal> {
        let format = self.reader.reader_mut().format();
        let mut header = Arc::new(self.reader.header().clone());
        let mut version = header_version(self.reader.header());
        let mut spare_records: Vec<Vec<Record>> = Vec::new();
        let mut spare_bufs: Vec<Vec<u8>> = Vec::new();
        let mut carry = Vec::new();
        // The messages of reading each record of the batches not yet replayed.
        let mut read_messages: BTreeMap<u64, Vec<(usize, Vec<Message>)>> = BTreeMap::new();
        let mut window = Window::new(self.threads * 2);
        let mut eof = false;
        let mut read_err: Option<Fatal> = None;
        loop {
            while !eof && !window.is_full() {
                let mut messages = Vec::new();
                let input = if !self.raw {
                    let mut records = spare_records.pop().unwrap_or_default();
                    records.resize_with(BATCH_RECORDS, Record::new);
                    let (mut len, mut bytes) = (0, 0);
                    while !records_batch_full(len, bytes) {
                        let (r, m) = log::capture(|| self.reader.next_record(&mut records[len]));
                        if !m.is_empty() {
                            messages.push((len, m));
                        }
                        match r {
                            Ok(true) => {
                                bytes += records[len].allocated_bytes();
                                len += 1;
                            }
                            Ok(false) => {
                                eof = true;
                                break;
                            }
                            Err(e) => {
                                read_err = Some(read_error(&e));
                                eof = true;
                                break;
                            }
                        }
                    }
                    Input::Records(records, len)
                } else {
                    let mut buf = spare_bufs.pop().unwrap_or_default();
                    buf.clear();
                    match format {
                        InputFormat::Vcf => {
                            match read_vcf_batch(self.reader, &mut carry, &mut buf) {
                                Ok(more) => eof = !more,
                                Err(_) => {
                                    // Whole lines read before the error are still collected.
                                    let keep = memchr::memrchr(b'\n', &buf).map_or(0, |i| i + 1);
                                    buf.truncate(keep);
                                    read_err = Some(read_error(&SyncedError::BcfRead));
                                    eof = true;
                                }
                            }
                            Input::Vcf(buf)
                        }
                        InputFormat::Bcf => {
                            let r = self.reader.reader_mut().read_bcf_batch(
                                &mut buf,
                                usize::MAX,
                                BATCH_BYTES,
                            );
                            match r {
                                Ok(0) => eof = true,
                                Ok(_) => {}
                                Err(_) => {
                                    // The records of a batch with a truncated one are dropped.
                                    buf.clear();
                                    read_err = Some(read_error(&SyncedError::BcfRead));
                                    eof = true;
                                }
                            }
                            Input::Bcf(buf)
                        }
                    }
                };
                let v = header_version(self.reader.header());
                if v != version {
                    header = Arc::new(self.reader.header().clone());
                    version = v;
                }
                let seq = window.send();
                read_messages.insert(seq, messages);
                let job = Job {
                    seq,
                    input,
                    header: Arc::clone(&header),
                };
                if self.job_tx.send(job).is_err() {
                    return Err(Fatal::error("Error: stats workers stopped\n".into()));
                }
            }
            if window.is_empty() {
                break;
            }
            let done = match self.done_rx.recv() {
                Ok(Reply::Done(done)) => done,
                // `collect` re-raises the worker's panic once every worker has stopped.
                Ok(Reply::Panicked) => {
                    return Err(Fatal::error("Error: a stats worker panicked\n".into()));
                }
                Err(_) => return Err(Fatal::error("Error: stats workers stopped\n".into())),
            };
            window.receive(done.seq, done);
            while let Some((seq, done)) = window.next_ready() {
                let read = read_messages.remove(&seq).unwrap_or_default();
                let Done {
                    input,
                    messages,
                    halt,
                    ..
                } = done;
                match halt {
                    None => replay(read, messages, None),
                    Some((i, Halt::Fatal(e))) => {
                        replay(read, messages, Some(i));
                        return Err(e);
                    }
                    Some((i, Halt::Stop(stop, start))) => {
                        replay(read, messages, Some(i));
                        match stop {
                            Stop::Parse => {
                                return Err(read_error(&self.reader.vcf_parse_error()));
                            }
                            Stop::Bcf => return Err(read_error(&SyncedError::BcfRead)),
                            Stop::HeaderUpdate => {
                                let Input::Vcf(buf) = &input else {
                                    unreachable!("only VCF lines need header updates")
                                };
                                self.finish_vcf_batch(&buf[start..])?;
                                header = Arc::new(self.reader.header().clone());
                                version = header_version(self.reader.header());
                            }
                        }
                    }
                }
                match input {
                    Input::Records(records, _) => spare_records.push(records),
                    Input::Vcf(buf) | Input::Bcf(buf) => spare_bufs.push(buf),
                }
            }
        }
        match read_err {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    /// Parses and collects the rest of a VCF batch on this thread, updating the header.
    fn finish_vcf_batch(&mut self, rest: &[u8]) -> Result<(), Fatal> {
        let w = &mut *self.own;
        let cfg = &w.setup.cfg;
        let mut pos = 0;
        while let Some(l) = next_line(rest, &mut pos, w.last_line_cr) {
            w.line.clear();
            w.line.extend_from_slice(l);
            if vcf::parse(self.reader.header_mut(), &mut w.line, &mut w.rec).is_err() {
                return Err(read_error(&self.reader.vcf_parse_error()));
            }
            process(
                &mut w.collector,
                cfg,
                &mut w.filter,
                self.reader.header(),
                &mut w.rec,
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use clap::{Args as _, Command, FromArgMatches as _};

    use super::*;
    use crate::raw_input::{TestDir, within_seconds};

    #[test]
    fn a_panicking_worker_ends_the_run() -> anyhow::Result<()> {
        // One worker panics on its batch while the others wait for jobs, holding the results
        // channel open: the main thread must stop waiting and re-raise the panic.
        let dir = TestDir::new("stats-panic");
        let input = dir.join("panic.vcf");
        std::fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 1\t1\t{}\tA\tG\t.\t.\t.\n",
                workers::PANIC_MARKER
            ),
        )?;
        // Collecting only: `run` would hold standard output, which the test harness then
        // waits for should the run hang.
        let name = input.display().to_string();
        let panicked = within_seconds("stats with a panicking worker", move || {
            panic::catch_unwind(|| -> anyhow::Result<()> {
                let command = super::super::Args::augment_args(Command::new("stats"));
                let matches = command.try_get_matches_from(["stats", "--threads", "3", &name])?;
                let args = super::super::Args::from_arg_matches(&matches)?;
                let options = super::super::LoopOptions::parse(&args, &matches)?;
                let files = [name.clone()];
                let run = super::super::Run {
                    files: &files,
                    sample_spec: None,
                    filter: options.filter,
                    depth: options.depth,
                    user: options.user,
                    verbose_sites: options.verbose_sites,
                };
                let synced = super::super::synced_options(
                    &args,
                    &matches,
                    options.regions_overlap,
                    options.targets_overlap,
                    3,
                )?;
                let mut reader = SyncedReader::open(&name, synced)?;
                let setup = super::super::single_setup(&args, run, &name, reader.header())?;
                reader.set_max_unpack(setup.max_unpack());
                let _ = collect(&mut reader, &setup, None, 3);
                Ok(())
            })
            .is_err()
        });
        assert!(panicked, "the worker's panic is re-raised");
        Ok(())
    }

    #[test]
    fn record_batches_are_cut_by_count_or_bytes() {
        assert!(!records_batch_full(BATCH_RECORDS - 1, BATCH_BYTES - 1));
        assert!(records_batch_full(BATCH_RECORDS, 0));
        // A few records of many samples fill a batch.
        assert!(records_batch_full(3, BATCH_BYTES));
    }

    #[test]
    fn a_slow_first_batch_holds_back_reading() {
        let mut window = Window::new(4);
        let mut sent = Vec::new();
        while !window.is_full() {
            sent.push(window.send());
        }
        assert_eq!(sent, [0, 1, 2, 3]);
        // Every later batch finishes while batch 0 is still being collected.
        for seq in 1..4 {
            window.receive(seq, seq);
            assert_eq!(window.next_ready(), None);
            assert!(window.is_full(), "batch {seq} finished before batch 0");
        }
        window.receive(0, 0);
        let replayed: Vec<_> = std::iter::from_fn(|| window.next_ready()).collect();
        assert_eq!(replayed, [(0, 0), (1, 1), (2, 2), (3, 3)]);
        assert!(window.is_empty());
    }

    #[test]
    fn the_window_holds_at_most_its_size_while_batch_0_lags() {
        // Batch 0 finishes last; the others finish as soon as they are sent. A driver that
        // always waits for a result while the window is not empty must still finish, and
        // never hold more than the window.
        const BATCHES: u64 = 100;
        let mut window = Window::new(4);
        let (mut sent, mut running, mut replayed) = (0, Vec::new(), Vec::new());
        let mut most_held = 0;
        loop {
            while sent < BATCHES && !window.is_full() {
                running.push(window.send());
                sent += 1;
            }
            if window.is_empty() {
                break;
            }
            // The driver waits for a result here: one must still be running, or it would
            // wait forever (`remove` panics on an empty list).
            let i = running.iter().rposition(|&s| s != 0).unwrap_or(0);
            let seq = running.remove(i);
            window.receive(seq, seq);
            most_held = most_held.max(window.waiting.len());
            replayed.extend(std::iter::from_fn(|| window.next_ready()).map(|(s, _)| s));
        }
        assert!(most_held <= 4, "held {most_held} results");
        assert_eq!(replayed, (0..BATCHES).collect::<Vec<_>>());
    }
}
