//! Parallel conversion: VCF lines are parsed and BCF records decoded on worker threads, then
//! formatted or encoded there, and written in order.
//!
//! Workers parse against a snapshot of the header. A line that needs a header change (an
//! undefined contig or tag) stops its batch there; the main thread parses the rest of that
//! batch itself, updating the header and reporting what HTSlib reports, and later batches get
//! the new header. Messages from workers are captured and printed in input order, so stderr
//! matches a serial run.

use std::{
    collections::{BTreeMap, VecDeque, btree_map::IntoValues},
    io,
    num::NonZero,
    panic,
    sync::Arc,
};

use crossbeam_channel::{Receiver, Sender};
use pasta::{
    util::threads,
    variant::{
        Header, Record, bcf,
        io::{Compression, InputFormat, RecordSpan, Writer},
        log::{self, Message},
        synced::{SyncedError, SyncedReader},
        vcf::{self, ParseError},
    },
};

use crate::{
    commands::view::ViewError,
    workers::{self, ReportPanic},
};

/// Bytes of input per batch.
const BATCH_BYTES: usize = 1 << 20;

enum Input {
    /// Whole lines, each ending with `\n` except perhaps the last of the file.
    Vcf(Vec<u8>),
    /// Whole BCF records.
    Bcf(Vec<u8>),
}

/// What workers turn records into.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Encoding {
    /// VCF text.
    Vcf,
    /// BCF records.
    Bcf,
}

impl Encoding {
    fn of(writer: &Writer) -> Self {
        if writer.output_type().is_vcf() {
            Self::Vcf
        } else {
            Self::Bcf
        }
    }
}

/// The input and output formats of a run.
#[derive(Clone, Copy, Debug)]
struct Formats {
    input: InputFormat,
    /// What happens to a `\r` ending the last line (see [`next_line`]).
    last_line_cr: LastLineCr,
    output: Encoding,
}

struct Job {
    seq: u64,
    input: Input,
    header: Arc<Header>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StopKind {
    /// The record needs a header update; the main thread parses from here.
    HeaderUpdate,
    /// The VCF line is malformed.
    Parse,
    /// The BCF record is malformed.
    Bcf,
    /// The record could not be written.
    Write,
}

/// What a worker sends back.
enum Reply {
    Done(Done),
    /// The worker panicked: its batch never comes back.
    Panicked,
}

struct Done {
    seq: u64,
    input: Input,
    data: Vec<u8>,
    spans: Vec<RecordSpan>,
    /// Messages of records, by record index within the batch.
    messages: Vec<(usize, Vec<Message>)>,
    /// Where the batch stopped, by record index (and byte offset into the input).
    stop: Option<(usize, usize, StopKind)>,
}

/// Records not yet written: the trailing ones at the same position. HTSlib's synced reader
/// returns a record only after reading the next one at a different position, so a read error
/// drops the records at the last position; holding them back here does the same.
#[derive(Default)]
struct Holdback {
    data: Vec<u8>,
    spans: Vec<RecordSpan>,
}

impl Holdback {
    /// Adds records and writes all but the trailing group at one position.
    fn push(
        &mut self,
        reader: &SyncedReader,
        writer: &mut Writer,
        data: &[u8],
        spans: &[RecordSpan],
    ) -> Result<(), ViewError> {
        if spans.is_empty() {
            return Ok(());
        }
        self.data.extend_from_slice(data);
        self.spans.extend_from_slice(spans);
        let last = *self.spans.last().expect("records were added");
        let mut k = self.spans.len();
        while k > 0 && self.spans[k - 1].rid == last.rid && self.spans[k - 1].pos == last.pos {
            k -= 1;
        }
        if k == 0 {
            return Ok(());
        }
        let bytes: usize = self.spans[..k].iter().map(|s| s.len as usize).sum();
        write_spans(reader, writer, &self.data[..bytes], &self.spans[..k])?;
        self.data.drain(..bytes);
        self.spans.drain(..k);
        Ok(())
    }

    /// Writes everything held (at the end of the input).
    fn flush(&mut self, reader: &SyncedReader, writer: &mut Writer) -> Result<(), ViewError> {
        write_spans(reader, writer, &self.data, &self.spans)?;
        self.data.clear();
        self.spans.clear();
        Ok(())
    }
}

/// How [`stream`] ended.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Streamed {
    /// Every record was streamed.
    Done,
    /// No worker thread could be started, and nothing was read: the caller streams the records
    /// itself.
    NoThreads,
}

/// Streams every record of `reader` to `writer` (bcftools view with no filtering), using
/// `threads` workers (fewer when the system cannot start that many).
pub fn stream(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    threads: usize,
) -> Result<Streamed, ViewError> {
    let formats = Formats {
        input: reader.reader_mut().format(),
        last_line_cr: LastLineCr::of(reader.reader_mut().compression()),
        output: Encoding::of(writer),
    };
    // Fewer workers when the system cannot start `threads`.
    let n = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
    let started = threads::start(n, |count| {
        let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(count.get() * 2);
        let (done_tx, done_rx) = crossbeam_channel::bounded::<Reply>(count.get() * 2);
        let workers = threads::spawn_all(count.get(), |_| {
            let job_rx = job_rx.clone();
            let done_tx = done_tx.clone();
            move || worker(&job_rx, &done_tx, formats)
        })?;
        Ok((workers, job_tx, done_rx))
    });
    let Ok((workers, job_tx, done_rx)) = started else {
        return Ok(Streamed::NoThreads);
    };
    let threads = workers.len();
    let result = drive(reader, writer, formats, &job_tx, &done_rx, threads);
    // After an error, a worker stops at its next result instead of working through the queue.
    drop(job_tx);
    drop(done_rx);
    for w in workers {
        if let Err(panic) = w.join() {
            panic::resume_unwind(panic);
        }
    }
    result.map(|()| Streamed::Done)
}

fn drive(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    formats: Formats,
    job_tx: &Sender<Job>,
    done_rx: &Receiver<Reply>,
    threads: usize,
) -> Result<(), ViewError> {
    let snapshot = Arc::new(reader.header().clone());
    let mut carry: Vec<u8> = Vec::new();
    let mut batches = InFlight::default();
    let mut eof = false;
    let mut recycled: Vec<Vec<u8>> = Vec::new();
    let mut held = Holdback::default();
    // A failed read ends the input: the whole records before it are still written, in order,
    // and the error is reported after them, as a serial run does.
    let mut read_failed = false;
    loop {
        while !eof && batches.len() < threads * 2 {
            let mut buf = recycled.pop().unwrap_or_default();
            buf.clear();
            let input = match formats.input {
                InputFormat::Vcf => {
                    match read_vcf_batch(reader, &mut carry, &mut buf) {
                        Ok(more) => eof = !more,
                        Err(_) => {
                            // The line being read when the block failed is lost.
                            buf.truncate(memchr::memrchr(b'\n', &buf).map_or(0, |i| i + 1));
                            carry.clear();
                            read_failed = true;
                            eof = true;
                        }
                    }
                    if buf.is_empty() {
                        break;
                    }
                    Input::Vcf(buf)
                }
                InputFormat::Bcf => {
                    if reader
                        .reader_mut()
                        .read_bcf_batch(&mut buf, usize::MAX, BATCH_BYTES)
                        .is_err()
                    {
                        read_failed = true;
                        eof = true;
                    }
                    if buf.is_empty() {
                        eof = true;
                        break;
                    }
                    Input::Bcf(buf)
                }
            };
            job_tx
                .send(Job {
                    seq: batches.next_seq(),
                    input,
                    header: Arc::clone(&snapshot),
                })
                .map_err(|_| workers_stopped())?;
        }
        if batches.is_empty() {
            if read_failed {
                // HTSlib's synced reader drops the records at the last position.
                return Err(ViewError::Synced(SyncedError::BcfRead));
            }
            return held.flush(reader, writer);
        }
        batches.receive(done_rx)?;
        while let Some(done) = batches.take() {
            let Some((input, start)) =
                finish_batch(reader, writer, done, &mut recycled, &mut held)?
            else {
                continue;
            };
            // Finish on this thread, from the inputs of the batches already sent.
            let mut chunks = VecDeque::from([(input, start)]);
            for done in batches.take_all(done_rx)? {
                let (Input::Vcf(b) | Input::Bcf(b)) = done.input;
                chunks.push_back((b, 0));
            }
            let source = TailSource {
                chunks,
                format: formats.input,
                last_line_cr: formats.last_line_cr,
                line: Vec::new(),
                partial: std::mem::take(&mut carry),
                read_failed,
            };
            return serial_tail(reader, writer, std::mem::take(&mut held), source);
        }
    }
}

/// The batches sent to the workers and not yet taken, numbered in input order. Results come
/// back in any order and are taken in input order.
///
/// Each count is derived from the sequence numbers and the results held, so a result already
/// received is never also waited for.
#[derive(Default)]
struct InFlight {
    /// The number of the next batch to send.
    next_send: u64,
    /// The number of the next batch to take; every batch before it was taken.
    next_take: u64,
    /// Results received and not taken yet, by sequence number.
    received: BTreeMap<u64, Done>,
}

impl InFlight {
    /// Numbers the next batch to send.
    fn next_seq(&mut self) -> u64 {
        let seq = self.next_send;
        self.next_send += 1;
        seq
    }

    /// The batches sent and not taken: those received and those still being worked on.
    fn len(&self) -> usize {
        (self.next_send - self.next_take) as usize
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The batches sent whose results have not been received.
    fn unreceived(&self) -> usize {
        self.len() - self.received.len()
    }

    /// Waits for the next result, in any order.
    fn receive(&mut self, done_rx: &Receiver<Reply>) -> Result<(), ViewError> {
        let done = match done_rx.recv() {
            Ok(Reply::Done(done)) => done,
            // `stream` re-raises the worker's panic once every thread has stopped.
            Ok(Reply::Panicked) => {
                return Err(ViewError::Fatal(anyhow::anyhow!("a view worker panicked")));
            }
            Err(_) => return Err(workers_stopped()),
        };
        debug_assert!(
            (self.next_take..self.next_send).contains(&done.seq),
            "batch {} was not sent or was taken (sent {}, taken {})",
            done.seq,
            self.next_send,
            self.next_take
        );
        let seq = done.seq;
        let earlier = self.received.insert(seq, done);
        debug_assert!(earlier.is_none(), "batch {seq} was received twice");
        Ok(())
    }

    /// Takes the next batch in input order, if it has been received.
    fn take(&mut self) -> Option<Done> {
        let done = self.received.remove(&self.next_take)?;
        self.next_take += 1;
        Some(done)
    }

    /// Waits for every result not yet received and takes all the batches, in input order.
    fn take_all(&mut self, done_rx: &Receiver<Reply>) -> Result<IntoValues<u64, Done>, ViewError> {
        while self.unreceived() > 0 {
            self.receive(done_rx)?;
        }
        self.next_take = self.next_send;
        Ok(std::mem::take(&mut self.received).into_values())
    }
}

/// The error for a channel to the workers that closed: a worker thread ended early.
fn workers_stopped() -> ViewError {
    ViewError::Fatal(anyhow::anyhow!("worker threads stopped"))
}

/// Reads about [`BATCH_BYTES`] of whole lines into `out`; returns `false` at the end of input.
pub(crate) fn read_vcf_batch(
    reader: &mut SyncedReader,
    carry: &mut Vec<u8>,
    out: &mut Vec<u8>,
) -> io::Result<bool> {
    out.append(carry);
    let r = reader.reader_mut();
    loop {
        let buf = r.fill_buf()?;
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

/// What `hts_getline` does with a `\r` ending a last line without a newline (a `\r` before
/// a newline is always stripped).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum LastLineCr {
    /// Kept, reading a plain file.
    Keep,
    /// Stripped, reading through BGZF or gzip.
    Strip,
}

impl LastLineCr {
    /// The behaviour for input with this compression.
    pub(crate) fn of(compression: Compression) -> Self {
        if compression == Compression::None {
            Self::Keep
        } else {
            Self::Strip
        }
    }
}

/// Splits the next line off `input` at `*pos`, as `hts_getline` returns it.
pub(crate) fn next_line<'a>(
    input: &'a [u8],
    pos: &mut usize,
    last_line_cr: LastLineCr,
) -> Option<&'a [u8]> {
    if *pos >= input.len() {
        return None;
    }
    let rest = &input[*pos..];
    let (mut line, newline) = match memchr::memchr(b'\n', rest) {
        Some(i) => {
            *pos += i + 1;
            (&rest[..i], true)
        }
        None => {
            *pos = input.len();
            (rest, false)
        }
    };
    if (newline || last_line_cr == LastLineCr::Strip) && line.last() == Some(&b'\r') {
        line = &line[..line.len() - 1];
    }
    Some(line)
}

fn worker(job_rx: &Receiver<Job>, done_tx: &Sender<Reply>, formats: Formats) {
    let _guard = ReportPanic::new(done_tx, || Reply::Panicked);
    let mut rec = Record::new();
    let mut line = Vec::new();
    while let Ok(job) = job_rx.recv() {
        let (Input::Vcf(input) | Input::Bcf(input)) = &job.input;
        workers::panic_if_marked(input);
        let header = &*job.header;
        let mut data = Vec::with_capacity(BATCH_BYTES * 2);
        let mut spans = Vec::new();
        let mut messages = Vec::new();
        let mut stop = None;
        let mut index = 0usize;
        let mut pos = 0usize;
        loop {
            let start = pos;
            let outcome = match &job.input {
                Input::Vcf(input) => {
                    let Some(l) = next_line(input, &mut pos, formats.last_line_cr) else {
                        break;
                    };
                    line.clear();
                    line.extend_from_slice(l);
                    let (result, msgs) =
                        log::capture(|| match vcf::parse_shared(header, &mut line, &mut rec) {
                            Ok(()) => emit(header, &mut rec, &mut data, formats.output),
                            Err(ParseError::HeaderUpdateNeeded) => Err(StopKind::HeaderUpdate),
                            Err(_) => Err(StopKind::Parse),
                        });
                    (result, msgs)
                }
                Input::Bcf(input) => {
                    if pos >= input.len() {
                        break;
                    }
                    let mut slice = &input[pos..];
                    let before = slice.len();
                    let (result, msgs) = log::capture(|| {
                        match bcf::read_record_raw(&mut slice, &mut rec) {
                            Ok(bcf::ReadOutcome::Record) => {}
                            _ => return Err(StopKind::Bcf),
                        }
                        if bcf::check_record(Some(header), &mut rec).is_err() {
                            return Err(StopKind::Bcf);
                        }
                        emit(header, &mut rec, &mut data, formats.output)
                    });
                    pos += before - slice.len();
                    (result, msgs)
                }
            };
            let (result, msgs) = outcome;
            match result {
                Ok(span) => {
                    if !msgs.is_empty() {
                        messages.push((index, msgs));
                    }
                    spans.push(span);
                }
                Err(kind) => {
                    if kind != StopKind::HeaderUpdate && !msgs.is_empty() {
                        messages.push((index, msgs));
                    }
                    stop = Some((index, start, kind));
                    break;
                }
            }
            index += 1;
        }
        let done = Done {
            seq: job.seq,
            input: job.input,
            data,
            spans,
            messages,
            stop,
        };
        if done_tx.send(Reply::Done(done)).is_err() {
            break;
        }
    }
}

/// Formats or encodes a parsed record onto `data`, returning its span.
fn emit(
    header: &Header,
    rec: &mut Record,
    data: &mut Vec<u8>,
    encoding: Encoding,
) -> Result<RecordSpan, StopKind> {
    if bcf::check_sample_count(header, rec).is_err() {
        return Err(StopKind::Write);
    }
    let before = data.len();
    let ok = match encoding {
        Encoding::Vcf => vcf::format_record(header, rec, data).is_ok(),
        Encoding::Bcf => bcf::encode_record(header, rec, data).is_ok(),
    };
    if !ok {
        data.truncate(before);
        return Err(StopKind::Write);
    }
    Ok(RecordSpan {
        len: (data.len() - before) as u32,
        rid: rec.rid(),
        pos: rec.pos(),
        rlen: rec.rlen(),
    })
}

/// Writes a finished batch up to where it stopped. Returns the batch input and the offset to
/// continue from serially when it stopped on a record the main thread must handle.
fn finish_batch(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    done: Done,
    recycled: &mut Vec<Vec<u8>>,
    held: &mut Holdback,
) -> Result<Option<(Vec<u8>, usize)>, ViewError> {
    let Done {
        input,
        data,
        spans,
        messages,
        stop,
        ..
    } = done;
    // Messages of the records before the stop (and of a malformed stopping record), in input
    // order; data is written up to the stop.
    let mut msgs = messages.into_iter().peekable();
    let mut offset = 0usize;
    for (i, span) in spans.iter().enumerate() {
        while let Some((_, m)) = msgs.next_if(|(j, _)| *j == i) {
            log::replay(m);
        }
        offset += span.len as usize;
    }
    held.push(reader, writer, &data[..offset], &spans)?;
    let input_buf = match input {
        Input::Vcf(b) | Input::Bcf(b) => b,
    };
    let Some((_, start, kind)) = stop else {
        recycled.push(input_buf);
        return Ok(None);
    };
    match kind {
        StopKind::Parse | StopKind::Bcf => {
            for (_, m) in msgs {
                log::replay(m);
            }
            Err(ViewError::Synced(if kind == StopKind::Parse {
                reader.vcf_parse_error()
            } else {
                SyncedError::BcfRead
            }))
        }
        // The stopping record is processed again serially, with its messages.
        StopKind::Write | StopKind::HeaderUpdate => Ok(Some((input_buf, start))),
    }
}

/// A parsed record waiting in the lookahead buffer, or one already encoded by a worker.
enum Entry {
    Encoded(Vec<u8>, RecordSpan),
    Parsed(Box<Record>),
}

impl Entry {
    fn position(&self) -> (i32, i64) {
        match self {
            Self::Encoded(_, s) => (s.rid, s.pos),
            Self::Parsed(r) => (r.rid(), r.pos()),
        }
    }
}

/// Where the serial tail reads records: the rest of the stopped batch, the inputs of batches
/// already dispatched, then the reader.
struct TailSource {
    chunks: VecDeque<(Vec<u8>, usize)>,
    format: InputFormat,
    last_line_cr: LastLineCr,
    line: Vec<u8>,
    /// The start of the next line, read before switching to serial processing.
    partial: Vec<u8>,
    /// The batch reader hit a read error: it is reported after the chunks.
    read_failed: bool,
}

impl TailSource {
    /// Reads and parses the next record into `rec` (`false` at the end).
    fn next(&mut self, reader: &mut SyncedReader, rec: &mut Record) -> Result<bool, ViewError> {
        loop {
            let Some((chunk, pos)) = self.chunks.front_mut() else {
                if self.read_failed {
                    return Err(ViewError::Synced(SyncedError::BcfRead));
                }
                return self.next_from_reader(reader, rec);
            };
            match self.format {
                InputFormat::Vcf => {
                    let Some(l) = next_line(chunk, pos, self.last_line_cr) else {
                        self.chunks.pop_front();
                        continue;
                    };
                    self.line.clear();
                    self.line.extend_from_slice(l);
                    if vcf::parse(reader.header_mut(), &mut self.line, rec).is_err() {
                        return Err(ViewError::Synced(reader.vcf_parse_error()));
                    }
                    return Ok(true);
                }
                InputFormat::Bcf => {
                    if *pos >= chunk.len() {
                        self.chunks.pop_front();
                        continue;
                    }
                    let mut slice = &chunk[*pos..];
                    let before = slice.len();
                    let ok = matches!(
                        bcf::read_record_raw(&mut slice, rec),
                        Ok(bcf::ReadOutcome::Record)
                    ) && bcf::check_record(Some(reader.header()), rec).is_ok();
                    *pos += before - slice.len();
                    if !ok {
                        return Err(ViewError::Synced(SyncedError::BcfRead));
                    }
                    return Ok(true);
                }
            }
        }
    }

    fn next_from_reader(
        &mut self,
        reader: &mut SyncedReader,
        rec: &mut Record,
    ) -> Result<bool, ViewError> {
        match self.format {
            InputFormat::Vcf => {
                let got = reader
                    .reader_mut()
                    .read_line(&mut self.line)
                    .map_err(|_| ViewError::Synced(SyncedError::BcfRead))?;
                if !self.partial.is_empty() {
                    let mut full = std::mem::take(&mut self.partial);
                    if !got && self.last_line_cr == LastLineCr::Strip && full.last() == Some(&b'\r')
                    {
                        full.pop();
                    }
                    full.extend_from_slice(&self.line);
                    self.line = full;
                } else if !got {
                    return Ok(false);
                }
                if vcf::parse(reader.header_mut(), &mut self.line, rec).is_err() {
                    return Err(ViewError::Synced(reader.vcf_parse_error()));
                }
                Ok(true)
            }
            InputFormat::Bcf => {
                let got = reader
                    .reader_mut()
                    .read_bcf_raw(rec)
                    .map_err(|_| ViewError::Synced(SyncedError::BcfRead))?;
                if got && bcf::check_record(Some(reader.header()), rec).is_err() {
                    return Err(ViewError::Synced(SyncedError::BcfRead));
                }
                Ok(got)
            }
        }
    }
}

/// Processes the rest of the input on this thread, with HTSlib's lookahead: a record is
/// written only after the next record at a different position has been read.
fn serial_tail(
    reader: &mut SyncedReader,
    writer: &mut Writer,
    held: Holdback,
    mut source: TailSource,
) -> Result<(), ViewError> {
    let mut buffer: VecDeque<Entry> = VecDeque::new();
    let mut off = 0usize;
    for span in held.spans {
        let len = span.len as usize;
        buffer.push_back(Entry::Encoded(held.data[off..off + len].to_vec(), span));
        off += len;
    }
    let mut eof = false;
    loop {
        let full = buffer.len() > 1
            && buffer.front().map(Entry::position) != buffer.back().map(Entry::position);
        if !full {
            while !eof {
                let mut rec = Box::new(Record::new());
                if !source.next(reader, &mut rec)? {
                    eof = true;
                    break;
                }
                buffer.push_back(Entry::Parsed(rec));
                if buffer.front().map(Entry::position) != buffer.back().map(Entry::position) {
                    break;
                }
            }
        }
        let Some(entry) = buffer.pop_front() else {
            return Ok(());
        };
        match entry {
            Entry::Encoded(data, span) => write_spans(reader, writer, &data, &[span])?,
            Entry::Parsed(mut rec) => writer
                .write_record(reader.header(), &mut rec)
                .map_err(ViewError::Write)?,
        }
    }
}

fn write_spans(
    reader: &SyncedReader,
    writer: &mut Writer,
    data: &[u8],
    spans: &[RecordSpan],
) -> Result<(), ViewError> {
    writer
        .write_spans(reader.header(), data, spans)
        .map_err(ViewError::Write)
}

#[cfg(test)]
mod tests {
    use std::{fmt::Write as _, fs, path::Path};

    use clap::{Args as _, Command, FromArgMatches as _};

    use super::*;
    use crate::{commands::view, raw_input, workers};

    /// A result of batch `seq` with no records.
    fn done(seq: u64) -> Done {
        Done {
            seq,
            input: Input::Vcf(Vec::new()),
            data: Vec::new(),
            spans: Vec::new(),
            messages: Vec::new(),
            stop: None,
        }
    }

    #[test]
    fn take_all_waits_only_for_results_not_received() {
        // Batch 2 comes back before batches 0 and 1; batch 0 is taken (and stops, say), so
        // only batch 1 is still to come.
        let (done_tx, done_rx) = crossbeam_channel::unbounded();
        let mut batches = InFlight::default();
        for want in 0..3 {
            assert_eq!(batches.next_seq(), want);
        }
        for seq in [2, 0] {
            done_tx
                .send(Reply::Done(done(seq)))
                .expect("the receiver is open");
            assert!(batches.receive(&done_rx).is_ok(), "batch {seq} is received");
        }
        assert_eq!(batches.take().map(|d| d.seq), Some(0));
        assert!(batches.take().is_none(), "batch 1 was not received");
        assert_eq!((batches.len(), batches.unreceived()), (2, 1));
        // Batch 1 is the last result sent: waiting for a second one would fail.
        done_tx
            .send(Reply::Done(done(1)))
            .expect("the receiver is open");
        drop(done_tx);
        let Ok(rest) = batches.take_all(&done_rx) else {
            panic!("take_all waits for a result already received");
        };
        assert_eq!(rest.map(|d| d.seq).collect::<Vec<_>>(), [1, 2]);
        assert!(batches.is_empty());
    }

    #[test]
    fn take_all_waits_for_every_result_when_none_was_received() {
        let (done_tx, done_rx) = crossbeam_channel::unbounded();
        let mut batches = InFlight::default();
        for _ in 0..3 {
            batches.next_seq();
        }
        for seq in [1, 2, 0] {
            done_tx
                .send(Reply::Done(done(seq)))
                .expect("the receiver is open");
        }
        drop(done_tx);
        let Ok(all) = batches.take_all(&done_rx) else {
            panic!("every result was sent");
        };
        assert_eq!(all.map(|d| d.seq).collect::<Vec<_>>(), [0, 1, 2]);
        assert!(batches.take().is_none());
    }

    /// Runs `view` with `options` on `input` into `out`.
    fn run_view(input: &Path, out: &Path, options: &[&str]) -> anyhow::Result<()> {
        let mut argv: Vec<String> = ["view", "--no-version", "-o"].map(String::from).to_vec();
        argv.push(out.display().to_string());
        argv.extend(options.iter().map(|o| o.to_string()));
        argv.push(input.display().to_string());
        let matches = view::Args::augment_args(Command::new("view")).try_get_matches_from(argv)?;
        view::run(view::Args::from_arg_matches(&matches)?, &matches).map(drop)
    }

    /// Runs `view` with `options` on `input` into `out`, failing if it takes too long.
    fn view_into(input: &Path, out: &Path, options: &[&str]) -> anyhow::Result<()> {
        let (input, out) = (input.to_owned(), out.to_owned());
        let options: Vec<String> = options.iter().map(|o| o.to_string()).collect();
        raw_input::within_seconds("view", move || {
            let options: Vec<&str> = options.iter().map(String::as_str).collect();
            run_view(&input, &out, &options)
        })
    }

    #[test]
    fn a_panicking_worker_ends_the_run() -> anyhow::Result<()> {
        // One worker panics on its batch while the others wait for jobs, holding the results
        // channel open: the main thread must stop waiting and re-raise the panic.
        let dir = raw_input::TestDir::new("pipeline-panic");
        let input = dir.join("panic.vcf");
        fs::write(
            &input,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 1\t1\t{}\tA\tG\t.\t.\t.\n",
                workers::PANIC_MARKER
            ),
        )?;
        for output in ["-Ov", "-Ou"] {
            let (input, out) = (input.clone(), dir.join("out"));
            let panicked = raw_input::within_seconds("view with a panicking worker", move || {
                panic::catch_unwind(panic::AssertUnwindSafe(|| {
                    run_view(&input, &out, &[output, "--threads", "3"])
                }))
                .is_err()
            });
            assert!(panicked, "{output}: the worker's panic is re-raised");
        }
        Ok(())
    }

    #[test]
    fn threads_finish_after_a_late_header_update() -> anyhow::Result<()> {
        // About 1.2 MB of lines, so two batches, with the only undefined INFO tag near the end
        // of the first: the short second batch is usually back first.
        let dir = raw_input::TestDir::new("pipeline-late");
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
             ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n",
        );
        for i in 0..24_000 {
            let extra = if i == 22_000 { ";XX=1" } else { "" };
            writeln!(
                text,
                "1\t{}\trs{i}\tA\tG\t50\tPASS\tDP={i}{extra}\tGT\t0/1\t1/1\t0/0",
                i + 1
            )?;
        }
        let input = dir.join("late.vcf");
        fs::write(&input, text)?;
        for threads in ["2", "4"] {
            // VCF output has every record, as without threads.
            let (serial, threaded) = (dir.join("t0.vcf"), dir.join(format!("t{threads}.vcf")));
            view_into(&input, &serial, &[])?;
            view_into(&input, &threaded, &["--threads", threads])?;
            assert!(
                fs::read(&serial)? == fs::read(&threaded)?,
                "--threads {threads}"
            );
            // BCF output fails at the record, as bcf_write does.
            let out = dir.join(format!("t{threads}.bcf"));
            let e = view_into(&input, &out, &["-Ou", "--threads", threads])
                .expect_err("BCF output fails at the undefined tag");
            assert_eq!(
                e.to_string(),
                format!("[main_vcfview] Error: cannot write to {}", out.display())
            );
        }
        Ok(())
    }
}
