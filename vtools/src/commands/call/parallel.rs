//! Calling on worker threads (`--threads`): the main thread reads batches of raw VCF lines
//! or BCF records; workers parse them, call the sites, and format or encode the output; the
//! main thread writes the batches in order and keeps the gVCF blocks.
//!
//! Everything a worker does to a record depends on that record alone, so the output is the
//! same as a serial run. The options that need records in order (`-C alleles`, `-i`,
//! streamed `-t` regions, `-r`) or a ploidy that depends on earlier records (an explicit
//! ploidy in `-S` for an undeclared sex) use the serial loop. The one exception is `-c` on
//! records without INFO/I16, which take the I16 values of the last record called with them
//! ([`CarriedI16`]).
//!
//! An unreadable record ends the input quietly, as in bcftools call, and HTSlib's synced
//! reader has by then read the records that share a position with the last one it returned
//! without returning them. The main thread replays that buffering ([`Lookahead`]), so it
//! writes, and reports the errors of, the records a serial run calls and no others.

use std::{
    collections::{BTreeMap, VecDeque},
    io,
    num::NonZero,
    thread,
};

use anyhow::bail;
use crossbeam_channel::{Receiver, Sender};
use pasta::{
    util::threads,
    variant::{
        Header, Record, Unpack, bcf,
        io::{InputFormat, RecordSpan, WriteError, Writer},
        log::{self, Message},
        synced::{SyncedError, SyncedReader},
        vcf::{self, ParseError},
    },
};

use crate::{
    output,
    pipeline::{LastLineCr, next_line},
    workers::ReportPanic,
};

use super::{
    Caller, Config, EofMarker, MainLoop, Prepared, Shared,
    ccall::I16Batch,
    check_ploidy,
    gvcf::{CalledSite, Gvcf, GvcfItem, GvcfRanges, GvcfSite},
    ploidy::Ploidy,
    prepare,
    samples::SampleSet,
    targets::Targets,
    write_one,
};

/// Bytes of input per batch, at least (VCF input comes in whole reader buffers); fewer in
/// tests, so that test inputs span batches.
const BATCH_BYTES: usize = if cfg!(test) { 1 << 10 } else { 1 << 20 };

/// Whether the parallel path handles these options and this input.
pub fn supported(state: &MainLoop<'_>, reader: &SyncedReader) -> bool {
    reader.is_plain_stream()
        && state.eof == EofMarker::Present
        && !state.shared.flags.constrain_alleles
        && !state.config.insert_missed
        && state.targets.as_deref().is_none_or(|t| !t.has_alleles())
        && state.set.sample2sex[..state.set.nsamples.min(state.set.sample2sex.len())]
            .iter()
            .all(|&s| s >= 0)
}

/// What workers read.
struct Ctx<'a> {
    rhdr: &'a Header,
    hdr: &'a Header,
    config: &'a Config,
    shared: &'a Shared,
    set: &'a SampleSet,
    ploidy: &'a Ploidy,
    targets: Option<&'a Targets>,
    gvcf: Option<&'a GvcfRanges>,
    out_vcf: bool,
    last_line_cr: LastLineCr,
}

enum Input {
    Vcf(Vec<u8>),
    Bcf(Vec<u8>),
}

struct Job {
    seq: u64,
    input: Input,
    /// Whether reading the input failed right after this batch.
    read_failed: bool,
    /// gVCF sites freed by the main thread, for the worker to fill again.
    #[expect(
        clippy::vec_box,
        reason = "sites travel boxed in gVCF items, and are reused box and all"
    )]
    sites: Vec<Box<GvcfSite>>,
    /// The I16 values that `-c` gives the batch's first records without INFO/I16.
    anno16: Vec<f32>,
}

/// The most freed gVCF sites handed to one job.
const SITES_PER_JOB: usize = 1 << 14;

/// One output record of a batch.
enum Entry {
    /// Formatted or encoded, next in the batch data.
    Data(RecordSpan),
    /// A gVCF item; `Write` items have their record, next in the batch data.
    Gvcf(GvcfItem, Option<RecordSpan>),
}

/// The position of a record, as the synced reader compares them.
type Key = (i32, i64);

/// A worker's result for one batch. Records are numbered from 0 in the batch.
struct Done {
    seq: u64,
    input: Input,
    /// The job's `read_failed`.
    read_failed: bool,
    /// The job's `anno16`.
    anno16: Vec<f32>,
    /// What the batch did with the I16 values of `-c`.
    i16: Option<I16Batch>,
    data: Vec<u8>,
    /// The output entries, with the number of the record each comes from.
    entries: VecDeque<(u32, Entry)>,
    /// The position of each record read.
    keys: Vec<Key>,
    /// Messages of reading each record, by record number (`keys.len()` for the failed read).
    read_messages: Vec<(u32, Vec<Message>)>,
    /// Messages of calling each record, by record number.
    call_messages: Vec<(u32, Vec<Message>)>,
    /// A fatal error calling a record; the records after it are only read.
    fatal: Option<(u32, anyhow::Error)>,
    /// Whether reading stopped at an unreadable record (after those in `keys`).
    failed: bool,
    /// Whether that was a failed read (a truncated or corrupt BGZF block, a malformed BCF
    /// record), which is fatal, rather than an unparsable VCF line, which ends the input.
    read_error: bool,
}

/// How [`run`] ended.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Called {
    /// Every record was called.
    Done,
    /// No worker thread could be started, and nothing was read: the caller runs the serial
    /// loop.
    NoThreads,
}

/// Runs the main loop with calling on worker threads, each with a copy of `caller`.
///
/// # Errors
///
/// Returns the first fatal error in input order, after writing the records before it.
pub fn run(
    state: &mut MainLoop<'_>,
    reader: &mut SyncedReader,
    writer: &mut Writer,
    caller: &Caller,
) -> anyhow::Result<Called> {
    let threads = state.config.threads;
    let rhdr = reader.header().clone();
    let format = reader.reader_mut().format();
    let last_line_cr = LastLineCr::of(reader.reader_mut().compression());
    let hdr = state.hdr;
    let out_name = state.out_name.to_owned();
    let mut gvcf = state.gvcf.take();
    let ranges = gvcf.as_deref().map(|g| g.ranges.clone());
    let ctx = Ctx {
        rhdr: &rhdr,
        hdr,
        config: state.config,
        shared: state.shared,
        set: state.set,
        ploidy: state.ploidy,
        targets: state.targets.as_deref(),
        gvcf: ranges.as_ref(),
        out_vcf: writer.output_type().is_vcf(),
        last_line_cr,
    };
    let result = thread::scope(|scope| -> anyhow::Result<Called> {
        // Fewer workers when the system cannot start `threads`. The command checked that a
        // thread starts, but other threads (the reader's, the writer's) may have taken the room
        // since.
        let Ok(idle) = threads::start_idle_scoped(
            scope,
            NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN),
        ) else {
            return Ok(Called::NoThreads);
        };
        let threads = idle.len();
        let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(threads * 2);
        let (done_tx, done_rx) = crossbeam_channel::bounded::<Reply>(threads * 2);
        for t in idle {
            let job_rx = job_rx.clone();
            let done_tx = done_tx.clone();
            let mut caller = caller.clone();
            let ctx = &ctx;
            t.run(move || worker(ctx, &mut caller, &job_rx, &done_tx));
        }
        drop(done_tx);
        drop(job_rx);
        let mut out = Output {
            writer: &mut *writer,
            gvcf: gvcf.as_deref_mut(),
            hdr,
            out_name: &out_name,
            site_pool: Vec::new(),
            recycled: Vec::new(),
            held: VecDeque::new(),
            lookahead: Lookahead::default(),
            released: 0,
            next_base: 0,
        };
        let r = drive(reader, &mut out, format, &job_tx, &done_rx, threads);
        drop(job_tx);
        // Let the workers finish what they have.
        while done_rx.recv().is_ok() {}
        r.map(|()| Called::Done)
    });
    state.gvcf = gvcf;
    if result? == Called::NoThreads {
        return Ok(Called::NoThreads);
    }
    state.finish(writer)?;
    Ok(Called::Done)
}

fn drive(
    reader: &mut SyncedReader,
    out: &mut Output<'_>,
    format: InputFormat,
    job_tx: &Sender<Job>,
    done_rx: &Receiver<Reply>,
    threads: usize,
) -> anyhow::Result<()> {
    let mut carry: Vec<u8> = Vec::new();
    let mut pending: BTreeMap<u64, Done> = BTreeMap::new();
    let mut carried = CarriedI16::default();
    let (mut next_seq, mut next_write, mut in_flight) = (0u64, 0u64, 0usize);
    let mut eof = false;
    loop {
        while !eof && in_flight < threads * 2 {
            let mut buf = out.recycled.pop().unwrap_or_default();
            buf.clear();
            let mut read_failed = false;
            let input = match format {
                InputFormat::Vcf => {
                    match read_vcf_batch(reader, &mut carry, &mut buf) {
                        Ok(more) => eof = !more,
                        Err(_) => {
                            eof = true;
                            read_failed = true;
                        }
                    }
                    if buf.is_empty() && !read_failed {
                        break;
                    }
                    Input::Vcf(buf)
                }
                InputFormat::Bcf => {
                    let n =
                        match reader
                            .reader_mut()
                            .read_bcf_batch(&mut buf, usize::MAX, BATCH_BYTES)
                        {
                            Ok(n) => n,
                            Err(_) => {
                                read_failed = true;
                                0
                            }
                        };
                    if n == 0 {
                        eof = true;
                        if buf.is_empty() && !read_failed {
                            break;
                        }
                    }
                    Input::Bcf(buf)
                }
            };
            let pool = &mut out.site_pool;
            let sites = pool.split_off(pool.len().saturating_sub(SITES_PER_JOB));
            send(
                job_tx,
                Job {
                    seq: next_seq,
                    input,
                    read_failed,
                    sites,
                    anno16: carried.values.clone(),
                },
            )?;
            next_seq += 1;
            in_flight += 1;
        }
        if in_flight == 0 {
            return out.finish();
        }
        let done = match done_rx.recv() {
            Ok(Reply::Done(done)) => *done,
            // The scope re-raises the worker's panic once every thread has stopped.
            Ok(Reply::Panicked) => bail!("a calling worker panicked"),
            Err(_) => bail!("worker threads stopped"),
        };
        pending.insert(done.seq, done);
        if let Some(job) = carried.check(&mut pending) {
            send(job_tx, job)?;
        }
        while next_write < carried.seq
            && let Some(done) = pending.remove(&next_write)
        {
            next_write += 1;
            in_flight -= 1;
            if out.add(done)? == Flow::Stop {
                // An unparsable line ends the input quietly, unless it is garbage from a
                // corrupt deflate stream of a plain gzip input.
                if reader.reader_mut().check_gzip_rest().is_err() {
                    let msg = format!("Error: {}", SyncedError::BcfRead);
                    return Err(output::InputReadFailure(msg).into());
                }
                return Ok(());
            }
        }
    }
}

/// Hands a job to the workers.
fn send(job_tx: &Sender<Job>, job: Job) -> anyhow::Result<()> {
    if job_tx.send(job).is_err() {
        bail!("worker threads stopped");
    }
    Ok(())
}

/// The I16 values that `-c` gives a record without INFO/I16: those of the last record called
/// with I16 (none before any).
///
/// A batch starts from the values known when it is read. The values it leaves depend on them
/// only when a record used them, so the values before each batch are known once the batches
/// before it are called; a batch whose records used other values is called again from those,
/// and the batches after it are checked once it is.
#[derive(Debug, Default)]
struct CarriedI16 {
    /// The values before batch `seq`.
    values: Vec<f32>,
    /// The first batch not checked (the batches before it can be written once called).
    seq: u64,
}

impl CarriedI16 {
    /// Checks the batches called in input order, up to one to call again, which it takes out
    /// of `pending`: returns its job.
    fn check(&mut self, pending: &mut BTreeMap<u64, Done>) -> Option<Job> {
        // Bit for bit, so that NaN values match.
        let same = |a: &[f32], b: &[f32]| {
            a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
        };
        while let Some(done) = pending.get_mut(&self.seq) {
            if let Some(I16Batch { took_start, last }) = done.i16.take() {
                if took_start && !same(&done.anno16, &self.values) {
                    // The values it leaves may be wrong too.
                    let done = pending.remove(&self.seq).expect("found above");
                    return Some(Job {
                        seq: done.seq,
                        input: done.input,
                        read_failed: done.read_failed,
                        sites: Vec::new(),
                        anno16: self.values.clone(),
                    });
                }
                if let Some(last) = last {
                    self.values = last;
                }
            }
            self.seq += 1;
        }
        None
    }
}

/// HTSlib's synced-reader buffering with one streamed reader (`_reader_fill_buffer` and
/// `next_line`): before returning a record, it reads on until a record at another position,
/// and while the last record read is at the same position as the first one buffered (on
/// another sequence), it reads one more before each record it returns.
#[derive(Debug)]
struct Lookahead {
    /// The positions of the records read but not returned.
    buffered: VecDeque<Key>,
    /// The number of records returned.
    returned: u64,
    /// Whether the next record is read before returning one.
    filling: bool,
}

impl Default for Lookahead {
    fn default() -> Self {
        Self {
            buffered: VecDeque::new(),
            returned: 0,
            filling: true,
        }
    }
}

impl Lookahead {
    /// Reads a record at `key`; returns the number of records returned before the next read.
    fn read(&mut self, key: Key) -> u64 {
        self.buffered.push_back(key);
        if self.buffered[0] != key {
            self.filling = false;
        }
        while !self.filling {
            self.buffered.pop_front();
            self.returned += 1;
            // The next `next_line`: `_reader_fill_buffer` returns at once when the last
            // record buffered is at another position than the first.
            self.filling = match (self.buffered.front(), self.buffered.back()) {
                (Some(first), Some(last)) => first.1 == last.1,
                _ => true,
            };
        }
        self.returned
    }

    /// At the end of the input: every record read is returned.
    fn finish(&mut self) -> u64 {
        self.returned += self.buffered.len() as u64;
        self.buffered.clear();
        self.returned
    }
}

/// Whether to go on reading.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Flow {
    Continue,
    Stop,
}

/// A batch whose output is not all written.
struct Held {
    done: Done,
    /// The number of the batch's first record in the input.
    base: u64,
    /// The number of records in the batch.
    len: u64,
    /// The data offset of the next entry.
    off: usize,
    /// The next calling messages to print.
    call_message: usize,
}

/// Writes the batches' output in input order, a record once the serial loop would call it.
struct Output<'a> {
    writer: &'a mut Writer,
    gvcf: Option<&'a mut Gvcf>,
    hdr: &'a Header,
    out_name: &'a str,
    /// gVCF sites freed by the blocks, for the workers.
    #[expect(
        clippy::vec_box,
        reason = "sites travel boxed in gVCF items, and are reused box and all"
    )]
    site_pool: Vec<Box<GvcfSite>>,
    /// Input buffers to read into again.
    recycled: Vec<Vec<u8>>,
    held: VecDeque<Held>,
    lookahead: Lookahead,
    /// The number of records called (whose output may be written).
    released: u64,
    /// The number of the next batch's first record.
    next_base: u64,
}

impl Output<'_> {
    /// Takes the next batch in input order.
    fn add(&mut self, mut done: Done) -> anyhow::Result<Flow> {
        let keys = std::mem::take(&mut done.keys);
        let read_messages = std::mem::take(&mut done.read_messages);
        let (failed, read_error) = (done.failed, done.read_error);
        let base = self.next_base;
        self.next_base += keys.len() as u64;
        self.held.push_back(Held {
            done,
            base,
            len: keys.len() as u64,
            off: 0,
            call_message: 0,
        });
        let mut messages = read_messages.into_iter().peekable();
        for (i, &key) in keys.iter().enumerate() {
            while let Some((_, m)) = messages.next_if(|(j, _)| *j as usize == i) {
                log::replay(m);
            }
            let returned = self.lookahead.read(key);
            self.release(returned)?;
        }
        // The failed read's messages.
        for (_, m) in messages {
            log::replay(m);
        }
        self.write()?;
        if read_error {
            let msg = format!("Error: {}", SyncedError::BcfRead);
            return Err(output::InputReadFailure(msg).into());
        }
        if failed {
            return Ok(Flow::Stop);
        }
        Ok(Flow::Continue)
    }

    /// At the end of the input: writes everything left.
    fn finish(&mut self) -> anyhow::Result<()> {
        let returned = self.lookahead.finish();
        self.release(returned)?;
        self.write()
    }

    /// Calls the records up to `upto` (excluded): prints their calling messages, and stops
    /// at a fatal error, after writing the records before it.
    fn release(&mut self, upto: u64) -> anyhow::Result<()> {
        while self.released < upto {
            let g = self.released;
            let h = self
                .held
                .iter_mut()
                .find(|h| g < h.base + h.len)
                .expect("released records are in held batches");
            let local = (g - h.base) as u32;
            let msgs = &mut h.done.call_messages;
            while h.call_message < msgs.len() && msgs[h.call_message].0 == local {
                log::replay(std::mem::take(&mut msgs[h.call_message].1));
                h.call_message += 1;
            }
            if h.done.fatal.as_ref().is_some_and(|(i, _)| *i == local) {
                let (_, e) = h.done.fatal.take().expect("checked above");
                self.write()?;
                return Err(e);
            }
            self.released += 1;
        }
        Ok(())
    }

    /// Writes the entries of the records called so far.
    fn write(&mut self) -> anyhow::Result<()> {
        while let Some(h) = self.held.front_mut() {
            let limit = self.released.saturating_sub(h.base).min(h.len) as u32;
            let mut run: Vec<RecordSpan> = Vec::new();
            let mut run_start = h.off;
            while h.done.entries.front().is_some_and(|(r, _)| *r < limit) {
                let (_, entry) = h.done.entries.pop_front().expect("checked above");
                match entry {
                    Entry::Data(span) => {
                        h.off += span.len as usize;
                        run.push(span);
                    }
                    Entry::Gvcf(item, span) => {
                        let data = &h.done.data;
                        write_spans(
                            self.writer,
                            self.hdr,
                            self.out_name,
                            &data[run_start..h.off],
                            &run,
                        )?;
                        run.clear();
                        let g = self.gvcf.as_deref_mut().expect("gVCF items come with -g");
                        let (writer, hdr, out_name) = (&mut *self.writer, self.hdr, self.out_name);
                        let freed = g.advance(hdr, Some(item), |line| {
                            write_one(writer, hdr, line, out_name)
                        })?;
                        self.site_pool.extend(freed);
                        if let Some(span) = span {
                            let len = span.len as usize;
                            write_spans(
                                self.writer,
                                self.hdr,
                                self.out_name,
                                &data[h.off..h.off + len],
                                &[span],
                            )?;
                            h.off += len;
                        }
                        run_start = h.off;
                    }
                }
            }
            write_spans(
                self.writer,
                self.hdr,
                self.out_name,
                &h.done.data[run_start..h.off],
                &run,
            )?;
            if u64::from(limit) < h.len {
                break;
            }
            let h = self.held.pop_front().expect("checked above");
            let (Input::Vcf(b) | Input::Bcf(b)) = h.done.input;
            self.recycled.push(b);
        }
        Ok(())
    }
}

/// Writes formatted records (`bcf_write1`).
fn write_spans(
    writer: &mut Writer,
    hdr: &Header,
    out_name: &str,
    data: &[u8],
    spans: &[RecordSpan],
) -> anyhow::Result<()> {
    if spans.is_empty() {
        return Ok(());
    }
    match writer.write_spans(hdr, data, spans) {
        Ok(()) => Ok(()),
        Err(WriteError::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => {
            Err(anyhow::Error::new(e))
        }
        Err(_) => bail!("[main_vcfcall] Error: failed to write to {out_name}"),
    }
}

/// The record could not be read: the input ends there.
struct Unreadable;

/// Parses a VCF line (`line` is a copy of `raw`, used as scratch space).
fn read_vcf(
    ctx: &Ctx<'_>,
    line: &mut Vec<u8>,
    raw: &[u8],
    rec: &mut Record,
) -> Result<(), Unreadable> {
    match vcf::parse_shared(ctx.rhdr, line, rec) {
        Ok(()) => Ok(()),
        Err(ParseError::HeaderUpdateNeeded) => {
            // The record names an undefined tag or contig, which bcftools call rejects
            // (after the parser's warning) when it calls the record.
            let mut h = ctx.rhdr.clone();
            line.clear();
            line.extend_from_slice(raw);
            vcf::parse(&mut h, line, rec).map_err(|_| Unreadable)
        }
        Err(_) => Err(Unreadable),
    }
}

/// Decodes the next BCF record of `input`.
fn read_bcf(ctx: &Ctx<'_>, input: &mut &[u8], rec: &mut Record) -> Result<(), Unreadable> {
    match bcf::read_record_raw(input, rec) {
        Ok(bcf::ReadOutcome::Record) => {}
        _ => return Err(Unreadable),
    }
    bcf::check_record(Some(ctx.rhdr), rec).map_err(|_| Unreadable)
}

/// What a worker sends back.
enum Reply {
    Done(Box<Done>),
    /// The worker panicked; its job never completes.
    Panicked,
}

fn worker(ctx: &Ctx<'_>, caller: &mut Caller, job_rx: &Receiver<Job>, done_tx: &Sender<Reply>) {
    // Tells the main thread when this worker panics, so that it stops instead of waiting for
    // a job that never completes.
    let _report = ReportPanic::new(done_tx, || Reply::Panicked);
    let mut rec = Record::new();
    let mut line = Vec::new();
    let mut ploidy = PloidyCache::new(ctx);
    let mut dp = Vec::new();
    while let Ok(job) = job_rx.recv() {
        let mut sites = job.sites;
        caller.start_batch(&job.anno16);
        let mut done = Done {
            seq: job.seq,
            input: job.input,
            read_failed: job.read_failed,
            anno16: job.anno16,
            i16: None,
            data: Vec::with_capacity(BATCH_BYTES * 2),
            entries: VecDeque::new(),
            keys: Vec::new(),
            read_messages: Vec::new(),
            call_messages: Vec::new(),
            fatal: None,
            failed: false,
            read_error: false,
        };
        let mut pos = 0usize;
        loop {
            let (read, msgs) = match &done.input {
                Input::Vcf(input) => {
                    let Some(l) = next_line(input, &mut pos, ctx.last_line_cr) else {
                        break;
                    };
                    line.clear();
                    line.extend_from_slice(l);
                    log::capture(|| read_vcf(ctx, &mut line, l, &mut rec))
                }
                Input::Bcf(input) => {
                    if pos >= input.len() {
                        break;
                    }
                    let mut slice = &input[pos..];
                    let before = slice.len();
                    let r = log::capture(|| read_bcf(ctx, &mut slice, &mut rec));
                    pos += before - slice.len();
                    r
                }
            };
            let index = done.keys.len() as u32;
            if !msgs.is_empty() {
                done.read_messages.push((index, msgs));
            }
            if read.is_err() {
                // The serial loop's synced reader fails on a malformed BCF record, and ends
                // the input at an unparsable VCF line.
                done.failed = true;
                done.read_error = matches!(done.input, Input::Bcf(_));
                break;
            }
            done.keys.push((rec.rid(), rec.pos()));
            if done.fatal.is_some() {
                // The serial loop reads ahead of the record that failed.
                continue;
            }
            let (result, msgs) = log::capture(|| {
                let bufs = Buffers {
                    dp: &mut dp,
                    sites: &mut sites,
                    data: &mut done.data,
                };
                process(ctx, caller, &mut rec, &mut ploidy, bufs)
            });
            if !msgs.is_empty() {
                done.call_messages.push((index, msgs));
            }
            match result {
                Ok(Some(entry)) => done.entries.push_back((index, entry)),
                Ok(None) => {}
                Err(e) => done.fatal = Some((index, e)),
            }
        }
        if !done.failed && job.read_failed {
            done.failed = true;
            done.read_error = true;
        }
        done.i16 = caller.end_batch();
        if done_tx.send(Reply::Done(Box::new(done))).is_err() {
            break;
        }
    }
}

/// Per-worker ploidy: sample ploidies for the sex ploidies at a position, recomputed when
/// those change.
struct PloidyCache {
    sex2ploidy: Vec<i32>,
    last: Vec<i32>,
    samples: Vec<u8>,
    /// Whether [`check_ploidy`] has passed `samples` as they stand.
    checked: bool,
}

impl PloidyCache {
    fn new(ctx: &Ctx<'_>) -> Self {
        let pmax = ctx.ploidy.max();
        let nsex = ctx.ploidy.nsex();
        Self {
            sex2ploidy: vec![0; nsex],
            last: vec![pmax; nsex],
            samples: vec![pmax as u8; ctx.set.nsamples],
            checked: false,
        }
    }

    /// Sets the sample ploidies at a position, failing on one the callers cannot genotype.
    fn update(&mut self, ctx: &Ctx<'_>, chr: &[u8], pos: i64) -> anyhow::Result<()> {
        if ctx.ploidy.nsex() > 0 {
            ctx.ploidy.query(chr, pos as u32, &mut self.sex2ploidy);
            if self.sex2ploidy != self.last {
                for (p, &sex) in self.samples.iter_mut().zip(&ctx.set.sample2sex) {
                    *p = self.sex2ploidy[sex as usize] as u8;
                }
                std::mem::swap(&mut self.sex2ploidy, &mut self.last);
                self.checked = false;
            }
        }
        if !self.checked {
            check_ploidy(ctx.hdr, chr, pos, &self.samples)?;
            self.checked = true;
        }
        Ok(())
    }
}

/// A worker's buffers for [`process`].
struct Buffers<'a> {
    /// FORMAT/DP of gVCF output.
    dp: &'a mut Vec<i32>,
    /// Freed gVCF sites to fill.
    #[expect(
        clippy::vec_box,
        reason = "sites travel boxed in gVCF items, and are reused box and all"
    )]
    sites: &'a mut Vec<Box<GvcfSite>>,
    /// The batch's formatted or encoded records.
    data: &'a mut Vec<u8>,
}

/// The main loop's work on one record, up to formatting it.
fn process(
    ctx: &Ctx<'_>,
    caller: &mut Caller,
    rec: &mut Record,
    ploidy: &mut PloidyCache,
    bufs: Buffers<'_>,
) -> anyhow::Result<Option<Entry>> {
    if !rec.errors().is_empty() {
        bail!("Error: could not parse the input VCF");
    }
    if let Some(t) = ctx.targets
        && !t.matches(ctx.rhdr, rec)
    {
        return Ok(None);
    }
    if let Some(map) = &ctx.set.samples_map {
        rec.subset(map);
    }
    let Prepared::Call { unseen } = prepare(ctx.config, ctx.shared, rec) else {
        return Ok(None);
    };
    rec.unpack(Unpack::ALL);
    ploidy.update(ctx, rec.chrom_safe(ctx.rhdr), rec.pos())?;
    let ret = caller.call(ctx.hdr, rec, ctx.shared, &ploidy.samples, unseen, None)?;
    if ctx.shared.flags.var_only && ret == 0 && ctx.gvcf.is_none() {
        return Ok(None);
    }
    let item = match ctx.gvcf {
        Some(g) => {
            let mut spare = bufs.sites.pop();
            let item = g.classify(
                ctx.hdr,
                rec,
                CalledSite::from_call(ret),
                bufs.dp,
                &mut spare,
            )?;
            bufs.sites.extend(spare);
            if matches!(item, GvcfItem::Collapse { .. }) {
                return Ok(Some(Entry::Gvcf(item, None)));
            }
            Some(item)
        }
        None => None,
    };
    let span = emit(ctx, rec, bufs.data)?;
    Ok(Some(match item {
        Some(item) => Entry::Gvcf(item, Some(span)),
        None => Entry::Data(span),
    }))
}

/// Formats or encodes a record onto `data`.
fn emit(ctx: &Ctx<'_>, rec: &mut Record, data: &mut Vec<u8>) -> anyhow::Result<RecordSpan> {
    let fail = || anyhow::anyhow!("[main_vcfcall] Error: failed to write to -");
    bcf::check_sample_count(ctx.hdr, rec).map_err(|_| fail())?;
    let before = data.len();
    let ok = if ctx.out_vcf {
        vcf::format_record(ctx.hdr, rec, data).is_ok()
    } else {
        bcf::encode_record(ctx.hdr, rec, data).is_ok()
    };
    if !ok {
        data.truncate(before);
        return Err(fail());
    }
    Ok(RecordSpan {
        len: (data.len() - before) as u32,
        rid: rec.rid(),
        pos: rec.pos(),
        rlen: rec.rlen(),
    })
}

/// Reads about [`BATCH_BYTES`] of whole lines into `out`; returns `false` at the end of input.
///
/// # Errors
///
/// Returns a read error, with `out` holding the whole lines read before it (the serial
/// reader cannot read the line the error cuts).
fn read_vcf_batch(
    reader: &mut SyncedReader,
    carry: &mut Vec<u8>,
    out: &mut Vec<u8>,
) -> io::Result<bool> {
    out.append(carry);
    let r = reader.reader_mut();
    loop {
        let buf = match r.fill_buf() {
            Ok(buf) => buf,
            Err(e) => {
                let whole = memchr::memrchr(b'\n', out).map_or(0, |i| i + 1);
                out.truncate(whole);
                return Err(e);
            }
        };
        if buf.is_empty() {
            return Ok(false);
        }
        let n = buf.len();
        out.extend_from_slice(buf);
        r.consume(n);
        if out.len() >= BATCH_BYTES
            && let Some(last) = memchr::memrchr(b'\n', out)
        {
            carry.extend_from_slice(&out[last + 1..]);
            out.truncate(last + 1);
            return Ok(true);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{fmt::Write as _, fs, path::Path};

    use clap::{Args as _, Command, FromArgMatches as _};

    use super::*;

    /// Runs `call` on `input`; returns the output records.
    fn call(input: &Path, options: &[&str]) -> anyhow::Result<String> {
        let out = input.with_extension(format!("{}.out", options.join("")));
        let out_name = out.to_str().expect("UTF-8 path");
        let argv = [
            "call",
            "--no-version",
            "-o",
            out_name,
            input.to_str().expect("UTF-8 path"),
        ];
        let matches = super::super::Args::augment_args(Command::new("call"))
            .try_get_matches_from(argv.iter().chain(options))?;
        super::super::run(super::super::Args::from_arg_matches(&matches)?, &matches)?;
        let text = fs::read_to_string(&out)?;
        fs::remove_file(&out)?;
        Ok(text
            .lines()
            .filter(|l| !l.starts_with("##"))
            .collect::<Vec<_>>()
            .join("\n"))
    }

    #[test]
    fn records_without_i16_take_the_values_of_earlier_batches() -> anyhow::Result<()> {
        // Records without I16 take those of the last record with them, which can be several
        // batches before (a batch of VCF input is one 128 KiB buffer here); some batches
        // start with records that read their own.
        let mut vcf = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
             ##INFO=<ID=I16,Number=16,Type=Float,Description=\"I16\">\n\
             ##INFO=<ID=QS,Number=R,Type=Float,Description=\"QS\">\n\
             ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
        );
        // Long IDs make few records fill a batch.
        let id = "x".repeat(200);
        for i in 0..6000 {
            let i16 = if i % 1201 == 0 {
                let d = i / 1201;
                format!("I16={d},2,3,4,100,1000,100,1000,240,7200,120,3600,10,100,10,100;")
            } else {
                String::new()
            };
            let pl = if i % 3 == 0 { "0,30,255" } else { "200,0,200" };
            writeln!(
                vcf,
                "1\t{}\t{id}\tA\tC\t0\t.\t{i16}QS=0.6,0.4\tPL\t{pl}\t0,30,255",
                i + 1
            )?;
        }
        let dir = std::env::temp_dir().join(format!("vtools-call-i16-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let input = dir.join("in.vcf");
        fs::write(&input, vcf)?;
        assert!(
            fs::metadata(&input)?.len() > 10 << 17,
            "input spans batches"
        );
        for variants_only in [&[][..], &["-v"]] {
            let serial = call(&input, &[&["-c"], variants_only].concat())?;
            assert!(
                serial.contains("DP4=4,2,3,4"),
                "the values of a record with I16"
            );
            for threads in ["--threads=2", "--threads=4"] {
                let options = [&["-c", threads], variants_only].concat();
                assert!(call(&input, &options)? == serial, "{options:?}");
            }
        }
        fs::remove_file(&input)?;
        fs::remove_dir(&dir)?;
        Ok(())
    }

    #[test]
    fn i16_stored_as_integers_or_short_is_carried_as_in_a_serial_run() -> anyhow::Result<()> {
        // A BCF whose INFO/I16 is stored as integers on some records (a header changed from
        // Integer to Float: HTSlib reads the bits of the integers as floats), and has fewer
        // than 16 values on others (which replace only the first values carried), with
        // records without I16 after them across batches. The short record soon after one
        // stored as integers is in a batch read before the values it starts from are known.
        let header = |ty: &str| {
            Header::parse(
                format!(
                    "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
                     ##INFO=<ID=I16,Number=16,Type={ty},Description=\"I16\">\n\
                     ##INFO=<ID=QS,Number=R,Type=Float,Description=\"QS\">\n\
                     ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                )
                .as_bytes(),
            )
        };
        let (mut float_hdr, mut int_hdr) = (header("Float")?, header("Integer")?);
        let dir = std::env::temp_dir().join(format!("vtools-call-i16-int-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let input = dir.join("in.bcf");
        let options = pasta::variant::io::WriteOptions {
            output_type: pasta::variant::io::OutputType::Bcf,
            level: None,
            threads: 0,
        };
        let mut writer = Writer::create(Some(&input), &options)?;
        writer.write_header(&float_hdr)?;
        let mut rec = Record::new();
        for i in 0..3000 {
            let full = "2,3,4,100,1000,100,1000,240,7200,120,3600,10,100,10,100;";
            let (i16, as_int) = match i % 701 {
                0 => (format!("I16={},{full}", i / 701 + 1), false),
                350 => (format!("I16={},{full}", i / 701 + 1), true),
                175 | 375 | 525 => ("I16=7,6,5,4;".to_owned(), false),
                _ => (String::new(), false),
            };
            let pl = if i % 3 == 0 { "0,30,255" } else { "200,0,200" };
            let line = format!(
                "1\t{}\t.\tA\tC\t0\t.\t{i16}QS=0.6,0.4\tPL\t{pl}\t0,30,255",
                i + 1
            );
            let hdr = if as_int { &mut int_hdr } else { &mut float_hdr };
            vcf::parse(hdr, &mut line.into_bytes(), &mut rec)?;
            writer.write_record(&float_hdr, &mut rec)?;
        }
        writer.finish()?;
        let serial = call(&input, &["-c"])?;
        for values in [
            "DP4=1,2,3,4;MQ=33",
            "DP4=7,6,5,4;MQ=22",
            "DP4=0,0,0,0;MQ=0",
            "DP4=7,6,5,4;MQ=0",
        ] {
            assert!(serial.contains(values), "{values}");
        }
        for threads in ["--threads=2", "--threads=4"] {
            assert!(call(&input, &["-c", threads])? == serial, "{threads}");
        }
        fs::remove_file(&input)?;
        fs::remove_dir(&dir)?;
        Ok(())
    }

    /// The records returned before each read, for records at these positions.
    fn returned(keys: &[Key]) -> Vec<u64> {
        let mut l = Lookahead::default();
        keys.iter().map(|&k| l.read(k)).collect()
    }

    #[test]
    fn lookahead_holds_back_the_records_at_the_last_position() {
        // Each record is returned once one at another position is read.
        assert_eq!(returned(&[(0, 1), (0, 2), (0, 3)]), [0, 1, 2]);
        // Records at one position wait for the next position.
        assert_eq!(returned(&[(0, 1), (0, 1), (0, 2), (0, 3)]), [0, 0, 2, 3]);
        assert_eq!(Lookahead::default().finish(), 0);
    }

    #[test]
    fn lookahead_reads_ahead_at_the_same_position_on_another_sequence() {
        // After (0,5),(0,5),(1,5): returning the first (0,5) leaves (0,5),(1,5) buffered, at
        // one position, so the reader reads (1,6) before returning the second.
        assert_eq!(returned(&[(0, 5), (0, 5), (1, 5), (1, 6)]), [0, 0, 1, 3]);
        let mut l = Lookahead::default();
        for k in [(0, 5), (1, 5)] {
            l.read(k);
        }
        assert_eq!(l.finish(), 2);
    }
}
