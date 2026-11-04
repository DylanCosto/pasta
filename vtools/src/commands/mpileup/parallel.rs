//! Piling up on several threads with the output of a serial run.
//!
//! Each span (a region, or a chromosome of indexed inputs) is cut into chunks, leaving out
//! what the inputs' indexes show no reads in (see [`plan_chunks`]). A chunk other
//! than a span's first starts reading a margin before its first position, to rebuild the
//! pileup's state there (the buffered reads, BAQ already applied to them, and mates waiting
//! for each other), and writes records from its first position on. Chunks run on worker
//! threads; the main thread takes them in order and checks that each chunk's state at its
//! first position is the state the previous chunk ended with (by a fingerprint of the whole
//! pileup state). If it is, the chunk's output is the serial output; if not (a read longer
//! than the margin had BAQ applied before the margin, say), the main thread carries on the
//! previous chunk's pileup through the chunk instead. Every chunk's walker opens each input,
//! so the chunks piled up ahead of the one being written are as many as the file descriptor
//! limit allows.
//!
//! Deep pileups are subsampled with `drand48`: each chunk takes the generator state after the
//! draws of the chunks before it, waiting for them when it first needs a number.

use std::{collections::BTreeMap, num::NonZero, sync::Arc, thread};

use pasta::{
    util::threads,
    variant::{
        Header, Record, bcf,
        io::{RecordSpan, Writer},
        vcf,
    },
};

use super::{
    SpanRange,
    errmod::Rng48,
    input::{AlnIndex, Descriptors, Occupancy},
    pileup::{
        Carried, Engine, InputFile, RngGate, RngSource, Shared, Sink, Span, WalkRange, Walker,
    },
};
use crate::workers::ReportPanic;

/// The shortest and longest chunk.
const MIN_CHUNK: i64 = 20_000;
const MAX_CHUNK: i64 = 200_000;

/// File descriptors left free for what is opened after the pileup starts, besides the
/// inputs (the output's index, say).
const SPARE_DESCRIPTORS: u64 = 16;

/// One chunk of a span.
#[derive(Debug)]
pub(super) struct Chunk {
    /// What the chunk's readers query.
    query: Span,
    /// The span's region (0-based, inclusive).
    beg: u32,
    end: u32,
    /// Whether this is the span's first chunk, which starts with the span.
    first: bool,
    /// Where BAQ starts being applied (the margin's start).
    realn_from: i64,
    /// The chunk's first position and the next chunk's.
    start: i64,
    next: i64,
}

/// How the output stores records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// VCF text.
    Vcf,
    /// BCF records.
    Bcf,
}

/// Records formatted or encoded into memory.
struct BufferSink {
    encoding: Encoding,
    data: Vec<u8>,
    spans: Vec<RecordSpan>,
}

impl Sink for BufferSink {
    fn emit(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<()> {
        let start = self.data.len();
        match self.encoding {
            Encoding::Vcf => vcf::format_record(hdr, rec, &mut self.data)
                .map_err(|e| anyhow::anyhow!("failed to format a record: {e:?}"))?,
            Encoding::Bcf => bcf::encode_record(hdr, rec, &mut self.data)
                .map_err(|e| anyhow::anyhow!("failed to encode a record: {e:?}"))?,
        }
        self.spans.push(RecordSpan {
            len: (self.data.len() - start) as u32,
            rid: rec.rid(),
            pos: rec.pos(),
            rlen: rec.rlen(),
        });
        Ok(())
    }
}

/// Records written as they come.
struct DirectSink<'w> {
    writer: &'w mut Writer,
}

impl Sink for DirectSink<'_> {
    fn emit(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<()> {
        self.writer.write_record(hdr, rec)?;
        Ok(())
    }
}

/// A chunk piled up by a worker.
struct ChunkResult {
    walker: Walker,
    engine: Engine,
    out: BufferSink,
    fp_in: u128,
    /// The BQBZ and MQSBZ values the chunk started with.
    initial_carried: Carried,
    draws: u64,
}

/// The margin before a chunk: long enough to hold every read that BAQ may be applied to
/// (reads longer than `-M` are left alone) and mates that overlap.
fn margin(shared: &Shared<'_>) -> i64 {
    let conf = shared.conf;
    if conf.flag & super::pileup::REALN == 0 {
        return 2_000;
    }
    (4 * i64::from(conf.max_read_len)).clamp(2_000, 200_000)
}

/// Cuts the spans into chunks, leaving out what the inputs' indexes show that no query
/// reads from (see [`plan_chunks`]).
fn plan(shared: &Shared<'_>, spans: &[SpanRange], threads: usize) -> Vec<Chunk> {
    plan_inputs(
        shared.files,
        &shared.header.lengths,
        spans,
        threads,
        margin(shared),
    )
}

/// [`plan_chunks`] for `files`, from their indexes (an input without one may have records
/// anywhere).
pub(super) fn plan_inputs(
    files: &[InputFile],
    lengths: &[i64],
    spans: &[SpanRange],
    threads: usize,
    m: i64,
) -> Vec<Chunk> {
    let occupancy: Vec<Option<Occupancy<'_>>> = files
        .iter()
        .map(|file| file.index.as_ref().map(AlnIndex::occupancy))
        .collect();
    let may_read = |tid: usize, beg: u32, end: u32| {
        let len = lengths.get(tid).copied().unwrap_or(0);
        occupancy
            .iter()
            .any(|o| o.as_ref().is_none_or(|o| o.may_read(tid, beg, end, len)))
    };
    plan_chunks(spans, lengths, threads, m, may_read)
}

/// Cuts the spans into chunks of a sequence of `lengths[tid]` bases, with a margin of `m`
/// before each but a span's first. `may_read(tid, beg, end)` tells whether a query of the
/// inputs from `beg` to `end` (0-based, inclusive) may read a record, from their indexes.
///
/// A span that no query reads from is left out: its walker would read nothing, and write
/// nothing. A chunk that no query reads from is read by the chunk before it, which then
/// reads on to the next chunk's first position, as the serial run does; chunks at the start
/// of a span that no query reads from are read by the first that it reads from. Each chunk
/// left thus has records to read, and opens every input once (twice if the main thread
/// carries on its pileup), however few positions the inputs cover.
pub(super) fn plan_chunks(
    spans: &[SpanRange],
    lengths: &[i64],
    threads: usize,
    m: i64,
    may_read: impl Fn(usize, u32, u32) -> bool,
) -> Vec<Chunk> {
    let mut chunks = Vec::new();
    for (span, beg, end) in spans {
        let whole = |chunks: &mut Vec<Chunk>| {
            chunks.push(Chunk {
                query: span.clone(),
                beg: *beg,
                end: *end,
                first: true,
                realn_from: i64::from(*beg),
                start: i64::from(*beg),
                next: i64::MAX,
            });
        };
        let Span::Region {
            tid,
            beg: query_beg,
            end: query_end,
        } = *span
        else {
            whole(&mut chunks);
            continue;
        };
        if !may_read(tid, query_beg, query_end) {
            continue;
        }
        let len = lengths.get(tid).copied().unwrap_or(0);
        let hi = i64::from(*end).min(len - 1);
        let total = hi - i64::from(*beg) + 1;
        let size = (total / (threads as i64 * 8))
            .clamp(MIN_CHUNK, MAX_CHUNK)
            .max(4 * m);
        if total <= 2 * size {
            whole(&mut chunks);
            continue;
        }
        let n = (total + size - 1) / size;
        let span_first = chunks.len();
        // Whether the last chunk has records to read.
        let mut last_reads = false;
        for j in 0..n {
            let start = i64::from(*beg) + j * size;
            let next = if j + 1 < n { start + size } else { i64::MAX };
            // What the chunk writes records for.
            let last_pos = if next == i64::MAX {
                query_end
            } else {
                (next - 1) as u32
            };
            let reads = may_read(tid, start as u32, last_pos);
            if chunks.len() > span_first && !(reads && last_reads) {
                let last: &mut Chunk = chunks.last_mut().expect("a chunk of the span");
                last.next = next;
                last_reads |= reads;
                continue;
            }
            let from = if j == 0 {
                i64::from(*beg)
            } else {
                (start - m).max(i64::from(*beg))
            };
            chunks.push(Chunk {
                query: Span::Region {
                    tid,
                    beg: from as u32,
                    end: *end,
                },
                beg: *beg,
                end: *end,
                first: j == 0,
                realn_from: from,
                start,
                next,
            });
            last_reads = reads;
        }
    }
    chunks
}

/// How many chunks may be piled up ahead of the one being written: two per thread, to keep
/// the threads busy while a slow chunk holds up the others, as far as `free` file descriptors
/// allow (`None` for no limit). Each chunk's walker keeps every input open (`readers`) until
/// the chunk after it is written, and the main thread's walker one more set, so `window + 1`
/// walkers have their inputs open at once, and as many as one per thread and the main
/// thread's are opening them (each one input at a time, with the indexes already loaded).
/// Returns `None` if not even one chunk fits; the serial run then opens each input once.
pub fn window(threads: usize, readers: Descriptors, free: Option<u64>) -> Option<usize> {
    let most = (2 * threads).max(2);
    let Some(free) = free else {
        return Some(most);
    };
    let needed = |window: usize| {
        let walkers = window as u64 + 1;
        let opening = window.min(threads) as u64 + 1;
        walkers
            .saturating_mul(readers.held)
            .saturating_add(opening.saturating_mul(readers.opening))
    };
    (1..=most).rev().find(|&w| needed(w) <= free)
}

/// The file descriptors that may still be opened: the soft `RLIMIT_NOFILE` limit less those
/// open now and [`SPARE_DESCRIPTORS`]. `None` if there is no limit.
pub fn free_descriptors() -> Option<u64> {
    let mut limit = libc::rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };
    // SAFETY: `getrlimit` only writes the limits into the struct it is given.
    if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut limit) } != 0
        || limit.rlim_cur == libc::RLIM_INFINITY
    {
        return None;
    }
    // Listing the open descriptors opens one more, which is counted as spare.
    let open = std::fs::read_dir("/dev/fd").map_or(0, Iterator::count) as u64;
    Some(limit.rlim_cur.saturating_sub(open + SPARE_DESCRIPTORS))
}

/// Piles up one chunk.
fn run_chunk<'scope, 'env>(
    scope: &'scope thread::Scope<'scope, 'env>,
    shared: &'env Shared<'env>,
    chunk: &Chunk,
    index: usize,
    gate: &Arc<RngGate>,
    encoding: Encoding,
) -> anyhow::Result<ChunkResult> {
    let mut walker = Walker::start(scope, shared, chunk.query.clone())?;
    let mut engine = Engine::new(
        &shared.conf.call,
        shared.nsmpl,
        Rng48::new(0),
        RngSource::Gate(Arc::clone(gate), index),
    );
    let mut out = BufferSink {
        encoding,
        data: Vec::new(),
        spans: Vec::new(),
    };
    let mut fp_in = 0;
    if !chunk.first {
        let lead = WalkRange {
            beg: chunk.beg,
            end: chunk.end,
            realn_from: chunk.realn_from,
            out_from: chunk.start,
            stop_at: chunk.start,
            warm_from: chunk.realn_from,
        };
        walker.walk(shared, &mut engine, &mut out, &lead)?;
        fp_in = walker.fingerprint();
    }
    let initial_carried = engine.carried();
    engine.set_carried(initial_carried);
    let body = WalkRange {
        beg: chunk.beg,
        end: chunk.end,
        realn_from: chunk.realn_from,
        out_from: chunk.start,
        stop_at: chunk.next,
        warm_from: i64::MAX,
    };
    walker.walk(shared, &mut engine, &mut out, &body)?;
    let draws = engine.draws();
    Ok(ChunkResult {
        walker,
        engine,
        out,
        fp_in,
        initial_carried,
        draws,
    })
}

/// How [`run`] ended.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Piled {
    /// Every span was piled up.
    Done,
    /// No worker thread could be started, and nothing was read: the caller runs serially.
    NoThreads,
}

/// Piles up the spans on `threads` workers, with at most `window` chunks (see [`window`])
/// waiting to be written, writing the output of a serial run.
///
/// # Errors
///
/// Returns an error if an input cannot be read or the output written.
pub fn run(
    shared: &Shared<'_>,
    spans: &[(Span, u32, u32)],
    threads: usize,
    window: usize,
    seed: i64,
    writer: &mut Writer,
    encoding: Encoding,
) -> anyhow::Result<Piled> {
    let chunks = plan(shared, spans, threads);
    let n = chunks.len();
    let gate = Arc::new(RngGate::new(Rng48::new(seed), n));
    thread::scope(|scope| -> anyhow::Result<Piled> {
        let (task_tx, task_rx) = crossbeam_channel::unbounded::<usize>();
        let (res_tx, res_rx) =
            crossbeam_channel::unbounded::<(usize, anyhow::Result<ChunkResult>)>();
        // Fewer workers when the system cannot start them all, and none when it cannot start
        // one.
        let idle = match NonZero::new(threads.min(n)).map(|w| threads::start_idle_scoped(scope, w))
        {
            Some(Ok(idle)) => idle,
            Some(Err(_)) => return Ok(Piled::NoThreads),
            None => Vec::new(),
        };
        for t in idle {
            let task_rx = task_rx.clone();
            let res_tx = res_tx.clone();
            let (chunks, gate) = (&chunks, &gate);
            t.run(move || {
                for i in task_rx {
                    // A panic reports the chunk as failed, so that the main thread stops
                    // waiting for it (the scope then propagates the panic).
                    let _guard = ReportPanic::new(&res_tx, || {
                        (i, Err(anyhow::anyhow!("a pileup worker panicked")))
                    });
                    let busy = shared.conf.call.threads.hold();
                    let result = run_chunk(scope, shared, &chunks[i], i, gate, encoding);
                    drop(busy);
                    if res_tx.send((i, result)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(res_tx);
        for i in 0..n.min(window) {
            let _ = task_tx.send(i);
        }
        let result = drive(
            scope, shared, &chunks, &gate, writer, &res_rx, &task_tx, window,
        );
        drop(task_tx);
        if result.is_err() {
            gate.abort();
        }
        result.map(|()| Piled::Done)
    })
}

/// Takes the chunks' results in order, checking and writing them.
#[expect(
    clippy::too_many_arguments,
    reason = "the main thread's view of the workers"
)]
fn drive<'scope, 'env>(
    scope: &'scope thread::Scope<'scope, 'env>,
    shared: &'env Shared<'env>,
    chunks: &[Chunk],
    gate: &RngGate,
    writer: &mut Writer,
    res_rx: &crossbeam_channel::Receiver<(usize, anyhow::Result<ChunkResult>)>,
    task_tx: &crossbeam_channel::Sender<usize>,
    window: usize,
) -> anyhow::Result<()> {
    let mut waiting: BTreeMap<usize, anyhow::Result<ChunkResult>> = BTreeMap::new();
    // The pileup whose state is that of the serial run, paused before the next chunk.
    let mut exact: Option<(Walker, Engine)> = None;
    // The BQBZ and MQSBZ values of the serial run, which its indel records repeat from the
    // last SNP with an alternate allele (also across regions).
    let mut carried = Carried::default();
    for (j, chunk) in chunks.iter().enumerate() {
        let result = loop {
            if let Some(r) = waiting.remove(&j) {
                break r;
            }
            let (i, r) = res_rx
                .recv()
                .map_err(|_| anyhow::anyhow!("a pileup worker stopped"))?;
            waiting.insert(i, r);
        };
        // The next chunk starts once this chunk's walker or the previous one is closed, so
        // that at most `window + 1` walkers have their inputs open.
        let refill = || {
            if j + window < chunks.len() {
                let _ = task_tx.send(j + window);
            }
        };
        let res = result?;
        // The chunk's records are right if what it used of the state carried into it is the
        // serial run's.
        let carried_ok = carried.same(&res.initial_carried, res.engine.carried_used());
        let previous = exact.take();
        let state_ok = chunk.first
            || previous
                .as_ref()
                .is_some_and(|(w, _)| w.fingerprint() == res.fp_in);
        if previous.is_none() && !chunk.first {
            // The span ended before this chunk.
            res.walker.close();
            refill();
            gate.finish(j, 0);
            continue;
        }
        if state_ok && carried_ok {
            if let Some((walker, _)) = previous {
                walker.close();
            }
            refill();
            write(writer, shared.out_header, &res.out)?;
            gate.finish(j, res.draws);
            carried = carried.updated(&res.engine.carried(), res.engine.carried_fresh());
            let ChunkResult {
                walker, mut engine, ..
            } = res;
            engine.set_carried(carried);
            exact = unfinished(walker, engine);
            continue;
        }
        // Run this chunk as the serial run would: carry on the previous chunk's pileup, or
        // start the span afresh.
        res.walker.close();
        refill();
        let (mut walker, mut engine) = match previous {
            Some(p) if !chunk.first => p,
            _ => (
                Walker::start(scope, shared, chunk.query.clone())?,
                Engine::new(
                    &shared.conf.call,
                    shared.nsmpl,
                    Rng48::new(0),
                    RngSource::Known,
                ),
            ),
        };
        engine.set_rng(gate.wait(j));
        engine.set_carried(carried);
        let before = engine.draws();
        let range = WalkRange {
            beg: chunk.beg,
            end: chunk.end,
            realn_from: i64::from(chunk.beg),
            out_from: chunk.start,
            stop_at: chunk.next,
            warm_from: i64::MAX,
        };
        // A thread of the budget, as each worker holds: the workers hold at most N of the 2N
        // and the helpers never wait, so one comes back soon.
        let busy = shared.conf.call.threads.hold();
        walker.walk(shared, &mut engine, &mut DirectSink { writer }, &range)?;
        drop(busy);
        gate.finish(j, engine.draws() - before);
        carried = carried.updated(&engine.carried(), engine.carried_fresh());
        engine.set_carried(carried);
        exact = unfinished(walker, engine);
    }
    Ok(())
}

/// The pileup to carry on into the next chunk, unless its span has been read to the end.
fn unfinished(walker: Walker, engine: Engine) -> Option<(Walker, Engine)> {
    if walker.is_done() {
        walker.close();
        return None;
    }
    Some((walker, engine))
}

fn write(writer: &mut Writer, hdr: &Header, out: &BufferSink) -> anyhow::Result<()> {
    writer.write_spans(hdr, &out.data, &out.spans)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const BAMS_8: Descriptors = Descriptors {
        held: 8,
        opening: 0,
    };

    #[test]
    fn the_window_is_two_chunks_per_thread_without_a_tight_limit() {
        assert_eq!(window(4, BAMS_8, None), Some(8));
        assert_eq!(window(4, BAMS_8, Some(1_000)), Some(8));
        assert_eq!(window(1, BAMS_8, None), Some(2));
    }

    #[test]
    fn the_window_shrinks_to_the_free_file_descriptors() {
        // Window w keeps w + 1 sets of the inputs open; BAM inputs open nothing else.
        assert_eq!(window(4, BAMS_8, Some(39)), Some(3), "4 sets = 32");
        assert_eq!(window(4, BAMS_8, Some(40)), Some(4), "5 sets = 40");
        assert_eq!(window(4, BAMS_8, Some(16)), Some(1));
        // 60 inputs with 1024 descriptors (bcftools opens each once) at --threads 8.
        let bams_60 = Descriptors {
            held: 60,
            opening: 0,
        };
        assert_eq!(window(8, bams_60, Some(1024 - 16 - 8)), Some(15));
        // 4 BAM and 4 CRAM inputs: each walker opening its inputs (as many as one per thread
        // and the main thread's) has one CRAM reference's FASTA index open at a time.
        let mixed = Descriptors {
            held: 4 + 4 * 2,
            opening: 1,
        };
        assert_eq!(
            window(4, mixed, Some(4 * 12 + 4)),
            Some(3),
            "4 sets + 4 opening"
        );
        assert_eq!(
            window(4, mixed, Some(4 * 12 + 3)),
            Some(2),
            "3 sets + 3 opening"
        );
        let few_held = Descriptors {
            held: 1,
            opening: 1,
        };
        assert_eq!(window(2, few_held, Some(8)), Some(4), "5 sets + 3 opening");
    }

    /// chr1, chr2 and chrM of GRCh38.
    const LENGTHS: [i64; 3] = [248_956_422, 242_193_529, 16_569];

    /// Whole-sequence spans, as `plan_spans` makes without regions.
    fn whole_sequences() -> Vec<SpanRange> {
        (0..LENGTHS.len())
            .map(|tid| {
                let end = super::super::regidx::REGIDX_MAX;
                (Span::Region { tid, beg: 0, end }, 0, u32::MAX)
            })
            .collect()
    }

    /// Plans whole sequences at --threads 4 with a 2 kb margin, where the inputs have records
    /// on `occupied` (sequence, first and last position).
    fn plan_with(occupied: &[(usize, u32, u32)]) -> Vec<Chunk> {
        plan_chunks(&whole_sequences(), &LENGTHS, 4, 2_000, |tid, beg, end| {
            occupied
                .iter()
                .any(|&(t, b, e)| t == tid && b <= end && beg <= e)
        })
    }

    #[test]
    fn empty_sequences_are_not_planned() {
        assert!(plan_with(&[]).is_empty());
    }

    #[test]
    fn a_sequence_with_records_in_one_window_is_one_chunk() {
        // The example: reads in one 10 kb window of chr1, which would make 1,245 chunks
        // of 200 kb, each opening every input.
        let chunks = plan_with(&[(0, 1_000_000, 1_010_100)]);
        assert_eq!(chunks.len(), 1, "{chunks:?}");
        let c = &chunks[0];
        assert!(c.first);
        assert_eq!((c.start, c.next, c.realn_from), (0, i64::MAX, 0));
        assert!(matches!(c.query, Span::Region { tid: 0, beg: 0, .. }));
    }

    #[test]
    fn chunks_without_records_are_read_by_the_chunk_before() {
        let chunks = plan_with(&[(0, 1_000_000, 1_000_100), (0, 100_000_050, 100_000_200)]);
        assert_eq!(chunks.len(), 2, "{chunks:?}");
        let (a, b) = (&chunks[0], &chunks[1]);
        assert!(a.first && !b.first);
        // The second chunk starts where it would have, with its margin, and the first reads
        // on to it.
        assert_eq!(b.start, 100_000_000);
        assert_eq!(b.realn_from, 100_000_000 - 2_000);
        assert!(matches!(
            b.query,
            Span::Region {
                tid: 0,
                beg: 99_998_000,
                ..
            }
        ));
        assert_eq!((a.start, a.next), (0, b.start));
        assert_eq!(b.next, i64::MAX);
    }

    #[test]
    fn sequences_with_records_throughout_are_cut_as_before() {
        let all: Vec<_> = (0..LENGTHS.len()).map(|t| (t, 0, u32::MAX)).collect();
        let chunks = plan_with(&all);
        // 200 kb chunks of chr1 and chr2; chrM is one chunk.
        let n = |len: i64| (len + MAX_CHUNK - 1) / MAX_CHUNK;
        assert_eq!(chunks.len() as i64, n(LENGTHS[0]) + n(LENGTHS[1]) + 1);
        for pair in chunks.windows(2) {
            if !pair[1].first {
                assert_eq!(pair[0].next, pair[1].start);
                assert_eq!(pair[1].start - pair[0].start, MAX_CHUNK);
            }
        }
    }

    #[test]
    fn every_input_is_read_without_an_index() {
        let spans = vec![(Span::All, 0, u32::MAX)];
        let chunks = plan_chunks(&spans, &LENGTHS, 4, 2_000, |_, _, _| false);
        assert_eq!(chunks.len(), 1);
        assert!(matches!(chunks[0].query, Span::All));
    }

    #[test]
    fn no_window_fits_when_two_sets_of_inputs_do_not() {
        assert_eq!(window(4, BAMS_8, Some(15)), None);
        assert_eq!(window(8, BAMS_8, Some(0)), None);
    }
}
