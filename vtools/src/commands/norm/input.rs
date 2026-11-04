//! Reading records, with VCF lines parsed on worker threads with `--threads`.
//!
//! Workers parse batches of lines against a snapshot of the header. A line that needs a
//! header change (an undefined contig or tag) or does not parse stops its batch; the main
//! thread parses the rest of that batch itself, updating the header and reporting what HTSlib
//! reports. Records are returned with HTSlib's synced-reader lookahead: a record is returned
//! only once one at a different position has been read, so a read error drops the records at
//! the last position, as in a serial run. A failed read ends the input: the whole lines read
//! before it are still parsed and returned, and the error follows them.

use std::{collections::BTreeMap, collections::VecDeque, num::NonZero, sync::Arc, thread};

use crossbeam_channel::{Receiver, Sender};
use pasta::{
    util::threads,
    variant::{
        Header, Record, Unpack,
        io::InputFormat,
        log::{self, Message},
        synced::{SyncedError, SyncedReader},
        vcf,
    },
};

use crate::{
    pipeline::{LastLineCr, next_line, read_vcf_batch},
    workers::{self, ReportPanic},
};

struct Job {
    seq: u64,
    input: Vec<u8>,
    header: Arc<Header>,
    max_unpack: Unpack,
    spare: Vec<Record>,
}

struct Done {
    seq: u64,
    input: Vec<u8>,
    recs: Vec<Record>,
    /// The messages of each parsed record.
    messages: Vec<(usize, Vec<Message>)>,
    /// Where parsing stopped (byte offset into the input), for the main thread to go on.
    stop: Option<usize>,
}

/// What a worker sends back.
enum Reply {
    Done(Done),
    /// The worker panicked: its batch never comes back.
    Panicked,
}

fn worker(jobs: &Receiver<Job>, done: &Sender<Reply>, last_line_cr: LastLineCr) {
    let _guard = ReportPanic::new(done, || Reply::Panicked);
    let mut line = Vec::new();
    while let Ok(job) = jobs.recv() {
        let Job {
            seq,
            input,
            header,
            max_unpack,
            mut spare,
        } = job;
        workers::panic_if_marked(&input);
        let mut recs = Vec::new();
        let mut messages = Vec::new();
        let mut stop = None;
        let mut pos = 0;
        loop {
            let start = pos;
            let Some(l) = next_line(&input, &mut pos, last_line_cr) else {
                break;
            };
            line.clear();
            line.extend_from_slice(l);
            let mut rec = spare.pop().unwrap_or_default();
            rec.set_max_unpack(max_unpack);
            let (ok, msgs) = log::capture(|| vcf::parse_shared(&header, &mut line, &mut rec));
            match ok {
                Ok(()) => {
                    if !msgs.is_empty() {
                        messages.push((recs.len(), msgs));
                    }
                    recs.push(rec);
                }
                Err(_) => {
                    // The main thread parses from here, reporting what it reports.
                    stop = Some(start);
                    break;
                }
            }
        }
        let d = Done {
            seq,
            input,
            recs,
            messages,
            stop,
        };
        if done.send(Reply::Done(d)).is_err() {
            break;
        }
    }
}

/// Reads VCF lines in batches and parses them on worker threads.
pub(crate) struct ParallelVcf {
    jobs: Option<Sender<Job>>,
    done: Receiver<Reply>,
    workers: Vec<thread::JoinHandle<()>>,
    snapshot: Arc<Header>,
    last_line_cr: LastLineCr,
    carry: Vec<u8>,
    eof: bool,
    next_seq: u64,
    next_take: u64,
    in_flight: usize,
    max_in_flight: usize,
    pending: BTreeMap<u64, Done>,
    /// Parsed records not returned yet, in input order.
    ready: VecDeque<Record>,
    /// Input left to parse on this thread: a buffer and the offset to start at.
    serial: Option<(Vec<u8>, usize)>,
    line: Vec<u8>,
    error: Option<SyncedError>,
    /// The failed read that ended the input, reported once the lines before it are returned.
    read_error: Option<SyncedError>,
    recycled: Vec<Vec<u8>>,
    /// How much of each line is parsed (empty: everything).
    max_unpack: Unpack,
}

impl Drop for ParallelVcf {
    fn drop(&mut self) {
        self.jobs = None;
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
    }
}

impl ParallelVcf {
    /// Starts parsing with `threads` workers (fewer when the system cannot start that many),
    /// when the input is a plain VCF stream; `None` reads on the calling thread, as it does when
    /// the system cannot start any.
    pub(crate) fn new(reader: &mut SyncedReader, threads: usize) -> Option<Self> {
        let n = NonZero::new(threads)?;
        if !reader.is_plain_stream() || reader.reader_mut().format() != InputFormat::Vcf {
            return None;
        }
        let last_line_cr = LastLineCr::of(reader.reader_mut().compression());
        let (workers, job_tx, done_rx) = threads::start(n, |count| {
            let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(count.get() * 2);
            let (done_tx, done_rx) = crossbeam_channel::unbounded::<Reply>();
            let workers = threads::spawn_all(count.get(), |_| {
                let (jobs, done) = (job_rx.clone(), done_tx.clone());
                move || worker(&jobs, &done, last_line_cr)
            })?;
            Ok((workers, job_tx, done_rx))
        })
        .ok()?;
        let threads = workers.len();
        Some(Self {
            jobs: Some(job_tx),
            done: done_rx,
            workers,
            snapshot: Arc::new(reader.header().clone()),
            last_line_cr,
            carry: Vec::new(),
            eof: false,
            next_seq: 0,
            next_take: 0,
            in_flight: 0,
            max_in_flight: threads * 2,
            pending: BTreeMap::new(),
            ready: VecDeque::new(),
            serial: None,
            line: Vec::new(),
            error: None,
            read_error: None,
            recycled: Vec::new(),
            max_unpack: Unpack::empty(),
        })
    }

    /// Limits how much of each line is parsed, as [`SyncedReader::set_max_unpack`] does.
    pub(crate) fn set_max_unpack(&mut self, which: Unpack) {
        self.max_unpack = which;
    }

    /// Returns the next record into `rec` (swapping it out), `false` at the end (HTSlib's
    /// `bcf_sr_next_line` on a plain stream).
    pub(crate) fn next(
        &mut self,
        reader: &mut SyncedReader,
        rec: &mut Record,
        spare: &mut Vec<Record>,
    ) -> Result<bool, SyncedError> {
        loop {
            let complete = match (self.ready.front(), self.ready.back()) {
                (Some(f), Some(b)) => f.rid() != b.rid() || f.pos() != b.pos(),
                _ => false,
            };
            if !complete && self.error.is_none() && self.exhausted() {
                // The lines before a failed read are all returned: the error follows them.
                self.error = self.read_error.take();
            }
            if complete || (self.error.is_none() && self.exhausted()) {
                let Some(front) = self.ready.pop_front() else {
                    return Ok(false);
                };
                spare.push(std::mem::replace(rec, front));
                return Ok(true);
            }
            if let Some(e) = self.error.take() {
                // The records at the last position are lost, as in HTSlib.
                spare.extend(self.ready.drain(..));
                return Err(e);
            }
            self.advance(reader, spare)?;
        }
    }

    /// Whether everything has been read and parsed.
    fn exhausted(&self) -> bool {
        self.eof && self.in_flight == 0 && self.serial.is_none()
    }

    /// Parses more records into `ready`.
    fn advance(
        &mut self,
        reader: &mut SyncedReader,
        spare: &mut Vec<Record>,
    ) -> Result<(), SyncedError> {
        if let Some((input, start)) = self.serial.take() {
            self.parse_serially(reader, input, start, spare);
            return Ok(());
        }
        while !self.eof && self.in_flight < self.max_in_flight {
            let mut buf = self.recycled.pop().unwrap_or_default();
            buf.clear();
            match read_vcf_batch(reader, &mut self.carry, &mut buf) {
                Ok(more) => self.eof = !more,
                Err(_) => {
                    // The line being read when the block failed is lost.
                    buf.truncate(memchr::memrchr(b'\n', &buf).map_or(0, |i| i + 1));
                    self.carry.clear();
                    self.read_error = Some(SyncedError::BcfRead);
                    self.eof = true;
                }
            }
            if buf.is_empty() {
                break;
            }
            let n = spare.len().min(4096);
            let job = Job {
                seq: self.next_seq,
                input: buf,
                header: Arc::clone(&self.snapshot),
                max_unpack: self.max_unpack,
                spare: spare.split_off(spare.len() - n),
            };
            let Some(jobs) = &self.jobs else { break };
            if jobs.send(job).is_err() {
                return Err(SyncedError::BcfRead);
            }
            self.next_seq += 1;
            self.in_flight += 1;
        }
        if self.in_flight == 0 {
            return Ok(());
        }
        while !self.pending.contains_key(&self.next_take) {
            let done = match self.done.recv() {
                Ok(Reply::Done(done)) => done,
                Ok(Reply::Panicked) => {
                    self.jobs = None;
                    workers::rethrow(self.workers.drain(..), "a norm parsing worker panicked")
                }
                Err(_) => return Err(SyncedError::BcfRead),
            };
            self.pending.insert(done.seq, done);
        }
        let done = self
            .pending
            .remove(&self.next_take)
            .expect("the batch was received");
        self.next_take += 1;
        self.in_flight -= 1;
        let Done {
            input,
            recs,
            messages,
            stop,
            ..
        } = done;
        let mut msgs = messages.into_iter().peekable();
        for (i, r) in recs.into_iter().enumerate() {
            while let Some((_, m)) = msgs.next_if(|(j, _)| *j == i) {
                log::replay(m);
            }
            self.ready.push_back(r);
        }
        match stop {
            Some(start) => self.serial = Some((input, start)),
            None => self.recycled.push(input),
        }
        Ok(())
    }

    /// Parses the rest of a stopped batch on this thread, with the live header.
    fn parse_serially(
        &mut self,
        reader: &mut SyncedReader,
        input: Vec<u8>,
        mut pos: usize,
        spare: &mut Vec<Record>,
    ) {
        while let Some(l) = next_line(&input, &mut pos, self.last_line_cr) {
            self.line.clear();
            self.line.extend_from_slice(l);
            let mut rec = spare.pop().unwrap_or_default();
            rec.set_max_unpack(self.max_unpack);
            if vcf::parse(reader.header_mut(), &mut self.line, &mut rec).is_err() {
                spare.push(rec);
                self.error = Some(reader.vcf_parse_error());
                return;
            }
            self.ready.push_back(rec);
        }
        // Later batches parse against the updated header.
        self.snapshot = Arc::new(reader.header().clone());
        self.recycled.push(input);
    }
}
