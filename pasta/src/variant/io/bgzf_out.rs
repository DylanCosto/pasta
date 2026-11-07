//! A BGZF writer that reports block addresses, so records' virtual offsets can be indexed while
//! writing, with optional compression on worker threads.
//!
//! Blocks hold up to [`BLOCK_SIZE`] bytes, as HTSlib's do. A block's compressed address is known
//! once every earlier block is compressed; [`BlockWriter::address`] returns it when known.
//!
//! A panic on a compression thread is re-raised on the calling thread by the next call that
//! hands over a block or waits for one.
//!
//! Only [`BlockWriter::finish`] completes the output. Dropping the writer before that abandons
//! it as a serial run that exits does: the blocks handed over are written, but not the block
//! being filled nor the EOF marker, and the inner writer is not flushed.

use std::{
    collections::BTreeMap,
    io::{self, Write},
    num::NonZero,
    panic,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender, TryRecvError};

use crate::{
    bgzf::io::{
        worker,
        writer::{BGZF_EOF, CompressionLevel, write_frame},
    },
    util::threads,
};

/// The uncompressed size of a full block (HTSlib's `BGZF_BLOCK_SIZE`).
pub const BLOCK_SIZE: usize = 0xff00;

/// Compresses a block into a complete BGZF frame.
fn compress_frame(
    data: &[u8],
    level: CompressionLevel,
    scratch: &mut Vec<u8>,
    out: &mut Vec<u8>,
) -> io::Result<()> {
    let crc = crate::bgzf::deflate::encode(data, level, scratch)?;
    out.clear();
    write_frame(out, scratch, crc, data.len())?;
    Ok(())
}

struct Job {
    seq: u64,
    data: Vec<u8>,
    /// Whether `data` holds complete BGZF frames to pass through as they are.
    raw: bool,
}

struct Done {
    seq: u64,
    /// The frame, or the panic of the worker that compressed it.
    frame: thread::Result<io::Result<Vec<u8>>>,
    data: Vec<u8>,
}

enum Mode<W: Write + Send + 'static> {
    Single {
        inner: W,
        scratch: Vec<u8>,
        frame: Vec<u8>,
    },
    Multi {
        jobs: Option<Sender<Job>>,
        recycled: Receiver<Vec<u8>>,
        sizes: Receiver<(u64, u64)>,
        workers: Vec<JoinHandle<()>>,
        writer: Option<JoinHandle<Result<W, Failure>>>,
        /// Set by [`BlockWriter::finish`] before it closes the job queue, for the writer thread
        /// to end the output with the EOF marker and flush it; without it the output is
        /// abandoned.
        finishing: Arc<AtomicBool>,
    },
    /// The writer thread stopped on an error (or the writer finished): every later call
    /// returns a copy of the error.
    Failed(io::Error),
}

/// The error that stopped the writer thread.
struct Failure {
    error: io::Error,
    /// The block it failed to write; `None` for the EOF marker and the final flush.
    block: Option<u64>,
}

/// Returns an error with the kind and message of `e`, to report it again.
fn copy_error(e: &io::Error) -> io::Error {
    io::Error::new(e.kind(), e.to_string())
}

/// Logs HTSlib's error for a block that could not be written (`bgzf_flush`; its compression
/// threads do not log). A closed pipe is not logged: bcftools dies of SIGPIPE before it could.
fn log_block_failure(e: &io::Error) {
    if e.kind() != io::ErrorKind::BrokenPipe {
        crate::variant::log::error("bgzf_flush", format_args!("File write failed (wrong size)"));
    }
}

/// Writes the EOF marker and flushes `inner`, logging HTSlib's error on failure (`bgzf_close`)
/// unless the reader of a pipe went away.
fn write_eof(inner: &mut impl Write) -> io::Result<()> {
    let result = inner.write_all(&BGZF_EOF).and_then(|()| inner.flush());
    if let Err(e) = &result
        && e.kind() != io::ErrorKind::BrokenPipe
    {
        crate::variant::log::error("bgzf_close", format_args!("File write failed"));
    }
    result
}

/// Records the compressed size of block `seq`, which follows the blocks known so far.
fn push_size(addresses: &mut Vec<u64>, seq: u64, size: u64) {
    debug_assert_eq!(seq as usize + 1, addresses.len(), "sizes arrive in order");
    let last = *addresses.last().expect("addresses start with block 0");
    addresses.push(last + size);
}

/// A BGZF writer that tracks block numbers and compressed addresses.
pub struct BlockWriter<W: Write + Send + 'static> {
    block: Vec<u8>,
    /// The number of the block being filled.
    block_no: u64,
    level: CompressionLevel,
    mode: Mode<W>,
    /// `addresses[k]` is the compressed address of block `k`, for the blocks known so far.
    addresses: Vec<u64>,
    /// The block the writer thread failed to write, once its error is known.
    failed_block: Option<u64>,
}

impl<W: Write + Send + 'static> BlockWriter<W> {
    /// Creates a writer compressing at `level` on `threads` worker threads (0 compresses on the
    /// calling thread, as it does when the system cannot start any; it may start fewer).
    pub fn new(inner: W, level: CompressionLevel, threads: usize) -> Self {
        let mode = match NonZero::new(threads) {
            None => Mode::Single {
                inner,
                scratch: Vec::new(),
                frame: Vec::new(),
            },
            Some(n) => spawn(inner, level, n),
        };
        Self {
            block: Vec::with_capacity(BLOCK_SIZE),
            block_no: 0,
            level,
            mode,
            addresses: vec![0],
            failed_block: None,
        }
    }

    /// Returns the block the compression threads failed to write, once that error has been
    /// returned: `None` while they work, after a failure on the EOF marker or the final flush,
    /// and without threads (where every error is returned by the call that wrote the block).
    #[must_use]
    pub fn failed_block(&self) -> Option<u64> {
        self.failed_block
    }

    /// Returns the current position as (block number, offset in block).
    #[inline]
    pub fn tell(&self) -> (u64, u16) {
        (self.block_no, self.block.len() as u16)
    }

    /// Returns the compressed address of block `block_no`, if known yet.
    ///
    /// # Errors
    ///
    /// Returns the I/O error that stopped the writer thread.
    #[inline]
    pub fn address(&mut self, block_no: u64) -> io::Result<Option<u64>> {
        self.poll()?;
        Ok(self.addresses.get(block_no as usize).copied())
    }

    /// Waits until the address of block `block_no` is known and returns it.
    ///
    /// # Errors
    ///
    /// Returns the I/O error that stopped the writer thread before it wrote the block.
    pub fn wait_address(&mut self, block_no: u64) -> io::Result<u64> {
        loop {
            if let Some(a) = self.addresses.get(block_no as usize) {
                return Ok(*a);
            }
            match &self.mode {
                Mode::Single { .. } => {
                    return Err(io::Error::other("block address requested before writing"));
                }
                Mode::Multi { sizes, .. } => match sizes.recv() {
                    Ok((seq, size)) => push_size(&mut self.addresses, seq, size),
                    Err(_) => return Err(self.writer_stopped()),
                },
                Mode::Failed(e) => return Err(copy_error(e)),
            }
        }
    }

    /// Records the sizes the writer thread has reported.
    ///
    /// # Errors
    ///
    /// Returns the I/O error that stopped the writer thread.
    fn poll(&mut self) -> io::Result<()> {
        let sizes = match &self.mode {
            Mode::Single { .. } => return Ok(()),
            Mode::Multi { sizes, .. } => sizes,
            Mode::Failed(e) => return Err(copy_error(e)),
        };
        loop {
            match sizes.try_recv() {
                Ok((seq, size)) => push_size(&mut self.addresses, seq, size),
                Err(TryRecvError::Empty) => return Ok(()),
                // The writer thread only stops early, before `finish`, on an error.
                Err(TryRecvError::Disconnected) => return Err(self.writer_stopped()),
            }
        }
    }

    /// Joins the threads after the writer thread stopped early, and returns the error that
    /// stopped it (HTSlib's writer thread failing a `bgzf_write`). Later calls return copies.
    ///
    /// # Panics
    ///
    /// Re-raises the panic of a compression job or of the writer thread.
    fn writer_stopped(&mut self) -> io::Error {
        if let Mode::Failed(e) = &self.mode {
            return copy_error(e);
        }
        let stopped = io::Error::other("BGZF compression threads stopped");
        let Mode::Multi {
            jobs,
            workers,
            writer,
            ..
        } = std::mem::replace(&mut self.mode, Mode::Failed(copy_error(&stopped)))
        else {
            unreachable!("only the threaded writer stops");
        };
        // Dropping the job queue stops the workers once they reach its end; each compresses
        // at most one more block, whose result the stopped writer thread no longer takes.
        drop(jobs);
        let error = match join(workers, writer) {
            Some(Err(Failure { error, block })) => {
                self.failed_block = block;
                error
            }
            Some(Ok(_)) | None => stopped,
        };
        self.mode = Mode::Failed(copy_error(&error));
        error
    }

    /// Writes data, starting new blocks as they fill (HTSlib's `bgzf_write`).
    ///
    /// # Errors
    ///
    /// Returns an I/O error from compression or the underlying writer.
    pub fn write_all(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            let n = (BLOCK_SIZE - self.block.len()).min(data.len());
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];
            if self.block.len() == BLOCK_SIZE {
                self.flush_block()?;
            }
        }
        Ok(())
    }

    /// Writes records, starting a new block before each one that would not fit in the current
    /// block (HTSlib's `vcf_write`: `bgzf_flush_try` before `bgzf_write`), so that blocks end at
    /// records unless a record alone is larger than a block. `lens` are the lengths of the
    /// records `data` holds.
    ///
    /// # Errors
    ///
    /// Returns an I/O error from compression or the underlying writer.
    pub fn write_records(
        &mut self,
        mut data: &[u8],
        lens: impl IntoIterator<Item = usize>,
    ) -> io::Result<()> {
        for len in lens {
            let (rec, rest) = data.split_at(len);
            data = rest;
            self.flush_try(len)?;
            self.write_all(rec)?;
        }
        debug_assert!(data.is_empty(), "the lengths cover the data");
        Ok(())
    }

    /// Starts a new block if `size` more bytes would not fit in this one (HTSlib's
    /// `bgzf_flush_try`).
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn flush_try(&mut self, size: usize) -> io::Result<()> {
        if self.block.len() + size > BLOCK_SIZE {
            self.flush_block()?;
        }
        Ok(())
    }

    /// Ends the current block if it has data (HTSlib's `bgzf_flush`).
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn flush_block(&mut self) -> io::Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }
        let seq = self.block_no;
        self.block_no += 1;
        match &mut self.mode {
            Mode::Single {
                inner,
                scratch,
                frame,
            } => {
                compress_frame(&self.block, self.level, scratch, frame)?;
                inner.write_all(frame).inspect_err(log_block_failure)?;
                let size = frame.len() as u64;
                self.block.clear();
                push_size(&mut self.addresses, seq, size);
            }
            Mode::Multi { jobs, recycled, .. } => {
                let next = recycled
                    .try_recv()
                    .unwrap_or_else(|_| Vec::with_capacity(BLOCK_SIZE));
                let data = std::mem::replace(&mut self.block, next);
                self.block.clear();
                let sent = jobs
                    .as_ref()
                    .expect("jobs are open until finish")
                    .send(Job {
                        seq,
                        data,
                        raw: false,
                    });
                if sent.is_err() {
                    return Err(self.writer_stopped());
                }
                self.poll()?;
            }
            Mode::Failed(e) => return Err(copy_error(e)),
        }
        Ok(())
    }

    /// Ends the current block and appends complete BGZF frames as they are (HTSlib's
    /// `bgzf_flush` and `bgzf_raw_write`), for copying compressed data without recompressing.
    /// The frames count as one block for [`Self::tell`] and [`Self::address`].
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn write_frames(&mut self, frames: &[u8]) -> io::Result<()> {
        self.flush_block()?;
        if frames.is_empty() {
            return Ok(());
        }
        let seq = self.block_no;
        self.block_no += 1;
        match &mut self.mode {
            Mode::Single { inner, .. } => {
                inner.write_all(frames)?;
                push_size(&mut self.addresses, seq, frames.len() as u64);
            }
            Mode::Multi { jobs, .. } => {
                let sent = jobs
                    .as_ref()
                    .expect("jobs are open until finish")
                    .send(Job {
                        seq,
                        data: frames.to_vec(),
                        raw: true,
                    });
                if sent.is_err() {
                    return Err(self.writer_stopped());
                }
                self.poll()?;
            }
            Mode::Failed(e) => return Err(copy_error(e)),
        }
        Ok(())
    }

    /// Flushes the last block, writes the EOF marker, and returns the inner writer and the
    /// address of the EOF marker. The writer cannot be used afterwards, but
    /// [`Self::failed_block`] still tells where a threaded write failed.
    ///
    /// # Errors
    ///
    /// Returns an I/O error.
    pub fn finish(&mut self) -> io::Result<(W, u64)> {
        self.flush_block()?;
        let end_block = self.block_no;
        let mode = std::mem::replace(
            &mut self.mode,
            Mode::Failed(io::Error::other("BGZF writer finished")),
        );
        match mode {
            Mode::Single { mut inner, .. } => {
                write_eof(&mut inner)?;
                let end = self.addresses[end_block as usize];
                Ok((inner, end))
            }
            Mode::Multi {
                jobs,
                sizes,
                workers,
                writer,
                finishing,
                ..
            } => {
                finishing.store(true, Ordering::Release);
                drop(jobs);
                let joined = join(workers, writer).expect("the writer thread runs until finish");
                let inner = match joined {
                    Ok(inner) => inner,
                    Err(Failure { error, block }) => {
                        self.failed_block = block;
                        return Err(error);
                    }
                };
                for (seq, size) in sizes.try_iter() {
                    push_size(&mut self.addresses, seq, size);
                }
                let end = *self
                    .addresses
                    .get(end_block as usize)
                    .ok_or_else(|| io::Error::other("missing block sizes"))?;
                Ok((inner, end))
            }
            Mode::Failed(e) => Err(e),
        }
    }
}

impl<W: Write + Send + 'static> Drop for BlockWriter<W> {
    /// Abandons an unfinished output, as bcftools' `error()` exits: the blocks handed to the
    /// compression threads are written, as a serial run writes each block when it is handed
    /// over, but not the block being filled nor the EOF marker, and the inner writer is not
    /// flushed. The threads are joined, so that nothing is written after the drop.
    fn drop(&mut self) {
        let Mode::Multi {
            jobs,
            workers,
            writer,
            ..
        } = &mut self.mode
        else {
            return;
        };
        // The workers compress what is queued and stop; the writer thread writes it and stops
        // without the EOF marker, as `finishing` is not set.
        drop(jobs.take());
        let (workers, writer) = (std::mem::take(workers), writer.take());
        if thread::panicking() {
            for w in workers {
                let _ = w.join();
            }
            if let Some(w) = writer {
                let _ = w.join();
            }
        } else {
            // A write error goes unreported, as the output is abandoned on an error already.
            let _ = join(workers, writer);
        }
    }
}

/// Joins the workers, once the job queue is closed, and then the writer thread, and returns its
/// result (`None` if it was joined before).
///
/// # Panics
///
/// Re-raises the panic of any of the threads.
fn join<W>(
    workers: Vec<JoinHandle<()>>,
    writer: Option<JoinHandle<Result<W, Failure>>>,
) -> Option<Result<W, Failure>> {
    for w in workers {
        w.join().unwrap_or_else(|p| panic::resume_unwind(p));
    }
    writer.map(|w| w.join().unwrap_or_else(|p| panic::resume_unwind(p)))
}

/// Starts the compression workers and the writer thread, with fewer workers when the system
/// cannot start `threads`, or compresses on the calling thread when it cannot start any (the
/// output is the same).
fn spawn<W: Write + Send + 'static>(
    inner: W,
    level: CompressionLevel,
    threads: NonZero<usize>,
) -> Mode<W> {
    let single = |inner| Mode::Single {
        inner,
        scratch: Vec::new(),
        frame: Vec::new(),
    };
    let started = threads::start(threads, |count| {
        let threads = count.get();
        let (job_tx, job_rx) = crossbeam_channel::bounded::<Job>(threads * 4);
        let (done_tx, done_rx) = crossbeam_channel::bounded::<Done>(threads * 4);
        let workers = threads::spawn_all(threads, |_| {
            let job_rx = job_rx.clone();
            let done_tx = done_tx.clone();
            move || {
                let mut scratch = Vec::new();
                while let Ok(Job { seq, data, raw }) = job_rx.recv() {
                    let (frame, data) = if raw {
                        (Ok(Ok(data)), Vec::new())
                    } else {
                        let frame = worker::catch(|| {
                            worker::panic_if_marked(&data);
                            let mut frame = Vec::with_capacity(BLOCK_SIZE / 2);
                            compress_frame(&data, level, &mut scratch, &mut frame).map(|()| frame)
                        });
                        (frame, data)
                    };
                    // A worker that panicked stops: the writer thread re-raises the panic when
                    // it reaches that block. Jobs are taken in order, so every earlier block
                    // has been taken by a worker that sends its result.
                    let panicked = frame.is_err();
                    if done_tx.send(Done { seq, frame, data }).is_err() || panicked {
                        break;
                    }
                }
            }
        })?;
        Ok((threads, workers, job_tx, done_rx))
    });
    let Ok((threads, workers, job_tx, done_rx)) = started else {
        return single(inner);
    };
    let (recycle_tx, recycle_rx) = crossbeam_channel::bounded::<Vec<u8>>(threads * 8);
    let (size_tx, size_rx) = crossbeam_channel::unbounded::<(u64, u64)>();
    // The first compression or write error, or a panic of a compression job, ends the writer
    // thread. That closes the sizes channel, so the main thread stops waiting for sizes and joins
    // this thread for the error (or re-raises the panic), and closes the results channel, so
    // workers stop. The output goes to the thread once it runs, so that it stays here when the
    // thread cannot be started.
    let (inner_tx, inner_rx) = crossbeam_channel::bounded::<W>(1);
    let finishing = Arc::new(AtomicBool::new(false));
    let finished = Arc::clone(&finishing);
    let writer = threads::spawn(move || -> Result<W, Failure> {
        let mut inner = inner_rx
            .recv()
            .expect("the output is sent once the thread runs");
        let mut pending: BTreeMap<u64, Done> = BTreeMap::new();
        let mut next = 0u64;
        while let Ok(done) = done_rx.recv() {
            pending.insert(done.seq, done);
            while let Some(Done { frame, data, .. }) = pending.remove(&next) {
                let failure = |error| Failure {
                    error,
                    block: Some(next),
                };
                let frame = frame
                    .unwrap_or_else(|p| panic::resume_unwind(p))
                    .map_err(failure)?;
                inner.write_all(&frame).map_err(failure)?;
                let _ = size_tx.send((next, frame.len() as u64));
                if data.capacity() > 0 {
                    let _ = recycle_tx.try_send(data);
                }
                next += 1;
            }
        }
        // An abandoned output (the writer dropped before `finish`) is left as it is.
        if finished.load(Ordering::Acquire) {
            write_eof(&mut inner).map_err(|error| Failure { error, block: None })?;
        }
        Ok(inner)
    });
    let writer = match writer {
        Ok(writer) => writer,
        Err(_) => {
            // The workers stop once the job channel closes.
            drop(job_tx);
            join::<W>(workers, None);
            return single(inner);
        }
    };
    inner_tx
        .send(inner)
        .expect("the writer thread waits for the output");
    Mode::Multi {
        jobs: Some(job_tx),
        recycled: recycle_rx,
        sizes: size_rx,
        workers,
        writer: Some(writer),
        finishing,
    }
}

#[cfg(test)]
pub(super) mod tests {
    use std::io::Read;

    pub(in crate::variant::io) use crate::bgzf::io::worker::within_timeout;

    use super::*;
    use crate::bgzf::io::worker::{MARKER_PANIC, PANIC_MARKER, assert_panics_with};

    /// A writer that fails with `StorageFull` (ENOSPC) once it has taken `left` bytes.
    pub(in crate::variant::io) struct FailAfter {
        pub left: usize,
    }

    impl Write for FailAfter {
        fn write(&mut self, data: &[u8]) -> io::Result<usize> {
            if self.left == 0 {
                return Err(io::ErrorKind::StorageFull.into());
            }
            let n = data.len().min(self.left);
            self.left -= n;
            Ok(n)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A writer that panics once it has taken `left` bytes.
    struct PanicAfter {
        left: usize,
    }

    impl Write for PanicAfter {
        fn write(&mut self, data: &[u8]) -> io::Result<usize> {
            self.left = (self.left.checked_sub(data.len())).expect("the writer met its limit");
            Ok(data.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_panic_on_a_thread_is_raised_again() {
        // Whether a block's address is waited for (as an index does), and whether the panic is
        // a compression job's or the writer thread's.
        for threads in [1, 2, 4] {
            for wait in [false, true] {
                for in_writer in [false, true] {
                    within_timeout(move || {
                        let message = if in_writer {
                            "the writer met its limit"
                        } else {
                            MARKER_PANIC
                        };
                        assert_panics_with(message, || {
                            let inner = PanicAfter {
                                left: if in_writer {
                                    3 * BLOCK_SIZE
                                } else {
                                    usize::MAX
                                },
                            };
                            let mut w = BlockWriter::new(inner, CompressionLevel::NONE, threads);
                            w.write_all(&vec![b'x'; 2 * BLOCK_SIZE])?;
                            if !in_writer {
                                w.write_all(PANIC_MARKER)?;
                                w.flush_block()?;
                            }
                            w.write_all(&vec![b'x'; 40 * BLOCK_SIZE])?;
                            if wait {
                                let (block, _) = w.tell();
                                w.wait_address(block)?;
                            }
                            w.finish().map(drop)
                        });
                    });
                }
            }
        }
    }

    #[test]
    fn a_write_error_reaches_every_waiter() {
        for threads in [1, 2, 4, 8] {
            let (write, wait, finish) = within_timeout(move || {
                // The error comes on the second block, with many blocks after it.
                let inner = FailAfter {
                    left: BLOCK_SIZE + 100,
                };
                let mut w = BlockWriter::new(inner, CompressionLevel::NONE, threads);
                let write = w.write_all(&vec![b'x'; 40 * BLOCK_SIZE]);
                let (block, _) = w.tell();
                let wait = w.wait_address(block);
                (write, wait, w.finish().map(|(_, end)| end))
            });
            if let Err(e) = write {
                assert_eq!(e.kind(), io::ErrorKind::StorageFull, "threads {threads}");
            }
            let kinds = (wait.map_err(|e| e.kind()), finish.map_err(|e| e.kind()));
            let full = Err(io::ErrorKind::StorageFull);
            assert_eq!(kinds, (full, full), "threads {threads}");
        }
    }

    #[test]
    fn a_failure_is_traced_to_its_block() {
        // Three full blocks, then a partial one that only `finish` flushes.
        let data = vec![b'x'; 3 * BLOCK_SIZE + 100];
        let mut serial = BlockWriter::new(Vec::new(), CompressionLevel::NONE, 0);
        serial.write_all(&data).expect("writing to memory");
        let (out, eof) = serial.finish().expect("writing to memory");
        let block_starts = [serial.addresses[1], serial.addresses[3], eof];
        assert_eq!(eof as usize + BGZF_EOF.len(), out.len());
        // The bytes written before the failure, and the block that fails.
        let cases = [
            (block_starts[0], Some(1)),
            (block_starts[1], Some(3)),
            (block_starts[2], None),
        ];
        for threads in [1, 2, 4] {
            for (left, failed) in cases {
                let data = data.clone();
                let (queued, finish, block) = within_timeout(move || {
                    let left = left as usize;
                    let mut w =
                        BlockWriter::new(FailAfter { left }, CompressionLevel::NONE, threads);
                    // The failure may be seen here already, or only by `finish`.
                    let _ = w.write_all(&data);
                    let (queued, _) = w.tell();
                    let finish = w.finish().map(drop).map_err(|e| e.kind());
                    (queued, finish, w.failed_block())
                });
                assert_eq!(queued, 3, "threads {threads}, left {left}");
                assert_eq!(
                    finish,
                    Err(io::ErrorKind::StorageFull),
                    "threads {threads}, left {left}"
                );
                assert_eq!(block, failed, "threads {threads}, left {left}");
            }
        }
    }

    #[test]
    fn writes_stop_at_the_first_error() {
        for threads in [1, 4] {
            let write = within_timeout(move || {
                let mut w =
                    BlockWriter::new(FailAfter { left: 0 }, CompressionLevel::NONE, threads);
                let block = vec![b'x'; BLOCK_SIZE];
                // Writes must fail once the error is seen, or this loop never ends.
                loop {
                    if let Err(e) = w.write_all(&block) {
                        return e;
                    }
                }
            });
            assert_eq!(
                write.kind(),
                io::ErrorKind::StorageFull,
                "threads {threads}"
            );
        }
    }

    #[test]
    fn round_trip_and_addresses() {
        for threads in [0, 3] {
            let mut w = BlockWriter::new(Vec::new(), CompressionLevel::default(), threads);
            let data: Vec<u8> = (0..200_000u32).map(|i| (i % 251) as u8).collect();
            w.write_all(&data).unwrap();
            let pos = w.tell();
            assert_eq!(pos.0, 3);
            let (out, end) = w.finish().unwrap();
            assert_eq!(end as usize, out.len() - BGZF_EOF.len());
            let mut r = crate::bgzf::io::Reader::new(&out[..]);
            let mut back = Vec::new();
            r.read_to_end(&mut back).unwrap();
            assert_eq!(back, data);
        }
    }

    #[test]
    fn raw_frames_pass_through_in_order() {
        let mut frame_src = BlockWriter::new(Vec::new(), CompressionLevel::default(), 0);
        frame_src.write_all(b"middle ").unwrap();
        let (frames, end) = frame_src.finish().unwrap();
        let frames = &frames[..end as usize];
        for threads in [0, 2] {
            let mut w = BlockWriter::new(Vec::new(), CompressionLevel::default(), threads);
            w.write_all(b"first ").unwrap();
            w.write_frames(frames).unwrap();
            w.write_frames(&[]).unwrap();
            w.write_all(b"last").unwrap();
            let (out, _) = w.finish().unwrap();
            let mut r = crate::bgzf::io::Reader::new(&out[..]);
            let mut back = Vec::new();
            r.read_to_end(&mut back).unwrap();
            assert_eq!(back, b"first middle last", "threads {threads}");
        }
    }
}
