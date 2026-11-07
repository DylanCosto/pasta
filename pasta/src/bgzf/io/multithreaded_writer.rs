//! Multithreaded BGZF writer.

mod builder;

use std::{
    io::{self, Write},
    mem,
    num::NonZero,
    panic,
    thread::{self, JoinHandle},
};

use bytes::BytesMut;
use crossbeam_channel::{Receiver, Sender};

pub use self::builder::Builder;
use super::{
    worker,
    writer::{CompressionLevel, MAX_BUF_SIZE},
};

/// The result of compressing a block, or the panic of the job that compressed it.
type BufferedRx = Receiver<thread::Result<io::Result<FrameParts>>>;
type WriteTx = Sender<Message>;
type WriteRx = Receiver<Message>;

/// A request to the thread that writes the compressed blocks, handled in the order sent.
enum Message {
    /// Write the block this job compresses.
    Block(BufferedRx),
    /// Flush the inner writer, then acknowledge on the flushed channel.
    Flush,
}

struct FrameParts {
    compressed_data: Vec<u8>,
    crc32: u32,
    uncompressed_size: usize,
}

enum State<W> {
    Running {
        writer_handle: JoinHandle<io::Result<W>>,
        write_tx: WriteTx,
        /// Acknowledges each [`Message::Flush`]; disconnected once the writer thread ends.
        flushed_rx: Receiver<()>,
    },
    /// The writer was finished, or failed with the error of this kind and message.
    Done(Option<(io::ErrorKind, String)>),
}

/// A multithreaded BGZF writer.
///
/// This is much more basic than [`super::Writer`] but uses a thread pool to compress block data.
///
/// A panic on a compression thread is re-raised on the calling thread by a later write, flush,
/// or [`Self::finish`].
///
/// Once writing fails (e.g., the inner writer's disk is full), every later write, flush, and
/// [`Self::finish`] returns an error of the same kind, as does using the writer once it is
/// finished.
pub struct MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    pool: crate::util::threads::Pool,
    state: State<W>,
    buf: BytesMut,
    compression_level: CompressionLevel,
}

impl<W> MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    /// Creates a multithreaded BGZF writer with a default worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let writer = bgzf::io::MultithreadedWriter::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self::with_worker_count(NonZero::<usize>::MIN, inner)
    }

    /// Creates a multithreaded BGZF writer with a worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    /// use pasta::bgzf;
    /// let writer = bgzf::io::MultithreadedWriter::with_worker_count(
    ///     NonZero::<usize>::MIN,
    ///     io::sink(),
    /// );
    /// ```
    ///
    /// # Panics
    ///
    /// Panics when the threads cannot be started; [`Self::try_with_worker_count`] returns that
    /// error.
    pub fn with_worker_count(worker_count: NonZero<usize>, inner: W) -> Self {
        Builder::default()
            .set_worker_count(worker_count)
            .build_from_writer(inner)
    }

    /// Creates a multithreaded BGZF writer with a worker count, or with fewer workers when the
    /// system cannot start that many threads.
    ///
    /// # Errors
    ///
    /// Fails when not even one worker thread, or the thread that writes the compressed blocks,
    /// can be started, returning `inner` in the error.
    pub fn try_with_worker_count(
        worker_count: NonZero<usize>,
        inner: W,
    ) -> Result<Self, crate::util::threads::StartError<W>> {
        Builder::default()
            .set_worker_count(worker_count)
            .try_build_from_writer(inner)
    }

    /// Finishes the output stream by flushing any remaining buffers.
    ///
    /// This shuts down the writer and deflater workers and appends the final BGZF EOF block.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let mut writer = bgzf::io::MultithreadedWriter::new(io::sink());
    /// writer.finish()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(&mut self) -> io::Result<W> {
        self.send_buffered()?;
        self.finish_inner()
    }

    fn finish_inner(&mut self) -> io::Result<W> {
        let state = mem::replace(&mut self.state, State::Done(None));

        match state {
            State::Running {
                writer_handle,
                write_tx,
                ..
            } => {
                drop(write_tx);
                let result = writer_handle
                    .join()
                    .unwrap_or_else(|p| panic::resume_unwind(p));

                if let Err(e) = &result {
                    // What could not be written is dropped, and later calls fail as this one.
                    self.buf.clear();
                    self.state = State::Done(Some((e.kind(), e.to_string())));
                }

                result
            }
            State::Done(failure) => {
                self.state = State::Done(failure);
                Err(self.done_error())
            }
        }
    }

    /// Returns the error of using the writer once it is finished or has failed.
    fn done_error(&self) -> io::Error {
        match &self.state {
            State::Done(Some((kind, message))) => {
                io::Error::new(*kind, format!("BGZF writer failed earlier: {message}"))
            }
            _ => io::Error::other("BGZF writer is finished"),
        }
    }

    fn remaining(&self) -> usize {
        MAX_BUF_SIZE - self.buf.len()
    }

    fn has_remaining(&self) -> bool {
        self.buf.len() < MAX_BUF_SIZE
    }

    /// Hands the buffered data, if any, to a compression job without waiting for it.
    fn send_buffered(&mut self) -> io::Result<()> {
        if matches!(self.state, State::Done(_)) {
            Err(self.done_error())
        } else if self.buf.is_empty() {
            Ok(())
        } else {
            self.send()
        }
    }

    fn send(&mut self) -> io::Result<()> {
        let State::Running { write_tx, .. } = &self.state else {
            return Err(self.done_error());
        };

        let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

        if write_tx.send(Message::Block(buffered_rx)).is_err() {
            return self.finish_inner().map(|_| ());
        }

        let src = self.buf.split().freeze();
        let compression_level = self.compression_level;

        self.pool.spawn(move || {
            let result = worker::catch(|| {
                worker::panic_if_marked(&src);
                compress(&src, compression_level)
            });
            buffered_tx.send(result).ok();
        });

        Ok(())
    }
}

impl<W> Drop for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn drop(&mut self) {
        if !matches!(self.state, State::Done(_)) {
            // A panic of a worker is not raised again while the thread is already unwinding.
            if thread::panicking() {
                let _ = worker::catch(|| self.finish());
            } else {
                let _ = self.finish();
            }
        }
    }
}

impl<W> Write for MultithreadedWriter<W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &self.state {
            State::Done(_) => return Err(self.done_error()),
            // The writer thread ends early only when it fails; collect its error without waiting.
            State::Running { writer_handle, .. } if writer_handle.is_finished() => {
                return self.finish_inner().and_then(|_| Err(self.done_error()));
            }
            State::Running { .. } => {}
        }

        let amt = self.remaining().min(buf.len());
        self.buf.extend_from_slice(&buf[..amt]);

        if !self.has_remaining() {
            self.send_buffered()?;
        }

        Ok(amt)
    }

    /// Writes every block buffered or sent so far to the inner writer and flushes it.
    fn flush(&mut self) -> io::Result<()> {
        self.send_buffered()?;

        let State::Running {
            write_tx,
            flushed_rx,
            ..
        } = &self.state
        else {
            return Err(self.done_error());
        };

        // The writer thread drops its acknowledging sender when it ends, so this cannot hang.
        if write_tx.send(Message::Flush).is_err() || flushed_rx.recv().is_err() {
            return self.finish_inner().and_then(|_| Err(self.done_error()));
        }

        Ok(())
    }
}

/// Starts the thread that writes the compressed blocks, once it is sent the writer on `start_rx`.
fn spawn_writer<W>(
    start_rx: Receiver<W>,
    write_rx: WriteRx,
    flushed_tx: Sender<()>,
) -> io::Result<JoinHandle<io::Result<W>>>
where
    W: Write + Send + 'static,
{
    use super::writer::{BGZF_EOF, write_frame};

    crate::util::threads::spawn(move || {
        let mut writer = start_rx
            .recv()
            .expect("the writer is sent once the thread runs");
        while let Ok(message) = write_rx.recv() {
            let buffered_rx = match message {
                Message::Block(buffered_rx) => buffered_rx,
                Message::Flush => {
                    writer.flush()?;
                    // The caller waits for this, unless the writer is being dropped.
                    flushed_tx.send(()).ok();
                    continue;
                }
            };
            // A job sends its result unless it panicked, which `worker::catch` turns into a
            // result; this thread re-raises it, and `finish` raises it again.
            let result = buffered_rx
                .recv()
                .map_err(|_| io::Error::other("BGZF compression job ended without a result"))?;
            let FrameParts {
                compressed_data,
                crc32,
                uncompressed_size,
            } = result.unwrap_or_else(|p| panic::resume_unwind(p))?;

            write_frame(&mut writer, &compressed_data, crc32, uncompressed_size)?;
        }

        writer.write_all(&BGZF_EOF)?;

        Ok(writer)
    })
}

fn compress(src: &[u8], compression_level: CompressionLevel) -> io::Result<FrameParts> {
    use crate::bgzf::deflate;

    let mut dst = Vec::new();
    let crc32 = deflate::encode(src, compression_level, &mut dst)?;

    Ok(FrameParts {
        compressed_data: dst,
        crc32,
        uncompressed_size: src.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf::io::worker::{MARKER_PANIC, PANIC_MARKER, assert_panics_with, within_timeout};

    /// A writer that panics on its first write.
    struct Panicking;

    impl Write for Panicking {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            panic!("the writer panicked");
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A writer that fails every write, as a full disk or a closed pipe does.
    struct Failing;

    impl Write for Failing {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(io::ErrorKind::BrokenPipe, "the pipe closed"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn calls_after_a_write_error_return_errors() {
        let is_broken_pipe = |e: &io::Error| e.kind() == io::ErrorKind::BrokenPipe;

        for workers in [1, 2, 4] {
            let workers = NonZero::new(workers).unwrap();
            within_timeout(move || {
                let new_failed_writer = || {
                    let mut writer = MultithreadedWriter::with_worker_count(workers, Failing);
                    let block = vec![b'x'; MAX_BUF_SIZE];
                    let e = (0..256)
                        .find_map(|_| writer.write_all(&block).err())
                        .expect("a write fails");
                    assert!(is_broken_pipe(&e), "{e}");
                    writer
                };

                let mut writer = new_failed_writer();
                assert!(writer.finish().is_err_and(|e| is_broken_pipe(&e)));
                assert!(writer.finish().is_err_and(|e| is_broken_pipe(&e)));

                let mut writer = new_failed_writer();
                assert!(writer.flush().is_err_and(|e| is_broken_pipe(&e)));

                let mut writer = new_failed_writer();
                assert!(writer.write(b"x").is_err_and(|e| is_broken_pipe(&e)));
                drop(writer);
            });
        }

        // A finished writer is not written to again.
        let mut writer = MultithreadedWriter::new(io::sink());
        writer.finish().expect("an empty stream is finished");
        assert!(writer.write(b"x").is_err());
        assert!(writer.finish().is_err());
    }

    #[test]
    fn calls_after_the_writer_thread_fails_return_errors() {
        let is_broken_pipe = |e: &io::Error| e.kind() == io::ErrorKind::BrokenPipe;

        within_timeout(move || {
            let mut writer = MultithreadedWriter::new(Failing);
            writer.write_all(&[b'x'; MAX_BUF_SIZE]).unwrap();
            let State::Running { writer_handle, .. } = &writer.state else {
                panic!("the writer is running");
            };
            while !writer_handle.is_finished() {
                thread::yield_now();
            }

            // A small write notices the failure without buffering.
            assert!(writer.write(b"x").is_err_and(|e| is_broken_pipe(&e)));
            assert!(writer.flush().is_err_and(|e| is_broken_pipe(&e)));
        });

        within_timeout(move || {
            let mut writer = MultithreadedWriter::new(Failing);
            writer.write_all(b"x").unwrap();
            // With buffered data and with none.
            assert!(writer.flush().is_err_and(|e| is_broken_pipe(&e)));
            assert!(writer.flush().is_err_and(|e| is_broken_pipe(&e)));
            assert!(writer.write(b"x").is_err_and(|e| is_broken_pipe(&e)));
        });
    }

    #[test]
    fn flush_writes_the_data_to_the_inner_writer() {
        use std::sync::{Arc, Mutex};

        /// A writer whose bytes, and whether it was flushed, stay readable while it is owned.
        #[derive(Clone, Default)]
        struct Shared(Arc<Mutex<(Vec<u8>, bool)>>);

        impl Write for Shared {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().0.extend_from_slice(buf);
                Ok(buf.len())
            }

            fn flush(&mut self) -> io::Result<()> {
                self.0.lock().unwrap().1 = true;
                Ok(())
            }
        }

        let inner = Shared::default();
        let mut writer =
            MultithreadedWriter::with_worker_count(NonZero::new(4).unwrap(), inner.clone());
        let data = vec![b'n'; 3 * MAX_BUF_SIZE + 5];
        writer.write_all(&data).unwrap();
        writer.flush().unwrap();

        let (written, flushed) = inner.0.lock().unwrap().clone();
        assert!(flushed);
        let mut decoded = Vec::new();
        io::Read::read_to_end(
            &mut crate::bgzf::io::Reader::new(&written[..]),
            &mut decoded,
        )
        .unwrap();
        assert_eq!(decoded, data);

        writer.finish().unwrap();
    }

    #[test]
    fn a_panic_on_a_thread_is_raised_again() {
        for workers in [1, 2, 4] {
            let workers = NonZero::new(workers).unwrap();
            // A compression job's panic (which used to abort the process).
            within_timeout(move || {
                assert_panics_with(MARKER_PANIC, || {
                    let mut writer = MultithreadedWriter::with_worker_count(workers, io::sink());
                    writer.write_all(&vec![b'x'; 2 * MAX_BUF_SIZE])?;
                    writer.write_all(PANIC_MARKER)?;
                    writer.flush()?;
                    writer.write_all(&vec![b'x'; 40 * MAX_BUF_SIZE])?;
                    writer.finish().map(drop)
                });
            });
            // The writer thread's panic, re-raised as it was.
            within_timeout(move || {
                assert_panics_with("the writer panicked", || {
                    let mut writer = MultithreadedWriter::with_worker_count(workers, Panicking);
                    writer.write_all(&vec![b'x'; 40 * MAX_BUF_SIZE])?;
                    writer.finish().map(drop)
                });
            });
        }
    }
}
