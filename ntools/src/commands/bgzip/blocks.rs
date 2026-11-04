//! BGZF blocks for `bgzip`: compression with HTSlib's block boundaries, and frame-by-frame
//! decompression, each on a pool of worker threads that keeps the blocks in order.

use std::{
    io::{self, Write},
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};
use pasta::bgzf::io::{raw, writer::CompressionLevel};

/// HTSlib's `BGZF_BLOCK_SIZE`: the most data `bgzf_write` puts in a block.
pub const BLOCK_SIZE: usize = 0xff00;

/// HTSlib's `BGZF_MAX_BLOCK_SIZE`: the block size `bgzf_block_write` uses past the end of the
/// index it follows.
const MAX_BLOCK_SIZE: usize = 0x10000;

/// How many blocks may be queued or in flight per worker.
const QUEUE_DEPTH_PER_WORKER: usize = 4;

/// Runs `work` on each pushed item, on `threads` worker threads (or on the calling thread for 0),
/// and passes the results to `sink` in the order the items were pushed.
pub struct Ordered<T, R, S> {
    state: State<T, R, S>,
}

enum State<T, R, S> {
    Inline {
        work: fn(T) -> R,
        sink: S,
    },
    Threaded {
        jobs: Sender<(T, Sender<R>)>,
        order: Sender<Receiver<R>>,
        workers: Vec<JoinHandle<()>>,
        sink: JoinHandle<io::Result<S>>,
    },
    Done,
}

/// Consumes results in order.
pub trait Sink<R>: Send + 'static {
    /// Handles the next result.
    ///
    /// # Errors
    ///
    /// Returns an error that ends the pipeline.
    fn accept(&mut self, result: R) -> io::Result<()>;
}

impl<T, R, S> Ordered<T, R, S>
where
    T: Send + 'static,
    R: Send + 'static,
    S: Sink<R>,
{
    pub fn new(threads: usize, work: fn(T) -> R, sink: S) -> Self {
        if threads == 0 {
            return Self {
                state: State::Inline { work, sink },
            };
        }

        let (jobs, jobs_rx) = crossbeam_channel::unbounded::<(T, Sender<R>)>();
        let (order, order_rx) =
            crossbeam_channel::bounded::<Receiver<R>>(threads * QUEUE_DEPTH_PER_WORKER);

        let workers = (0..threads)
            .map(|_| {
                let jobs_rx = jobs_rx.clone();
                thread::spawn(move || {
                    for (item, result_tx) in jobs_rx {
                        result_tx.send(work(item)).ok();
                    }
                })
            })
            .collect();

        let sink = thread::spawn(move || {
            let mut sink = sink;

            for result_rx in order_rx {
                let Ok(result) = result_rx.recv() else {
                    break;
                };

                sink.accept(result)?;
            }

            Ok(sink)
        });

        Self {
            state: State::Threaded {
                jobs,
                order,
                workers,
                sink,
            },
        }
    }

    /// Queues `item`, blocking while the queue is full.
    ///
    /// # Errors
    ///
    /// Returns the sink's error, if it has stopped.
    pub fn push(&mut self, item: T) -> io::Result<()> {
        match &mut self.state {
            State::Inline { work, sink } => sink.accept(work(item)),
            State::Threaded { jobs, order, .. } => {
                let (result_tx, result_rx) = crossbeam_channel::bounded(1);

                // The receiver goes first, so that a full queue holds back new work.
                if order.send(result_rx).is_ok() && jobs.send((item, result_tx)).is_ok() {
                    Ok(())
                } else {
                    // The sink stopped: report its error.
                    self.stop().map(|_| ())
                }
            }
            State::Done => Err(io::Error::other("the bgzip pipeline already stopped")),
        }
    }

    /// Waits for every queued item and returns the sink.
    ///
    /// # Errors
    ///
    /// Returns the sink's error.
    pub fn finish(mut self) -> io::Result<S> {
        self.stop()
    }

    fn stop(&mut self) -> io::Result<S> {
        match std::mem::replace(&mut self.state, State::Done) {
            State::Inline { sink, .. } => Ok(sink),
            State::Threaded {
                jobs,
                order,
                workers,
                sink,
            } => {
                drop(order);
                drop(jobs);

                for worker in workers {
                    worker.join().expect("a bgzip worker panicked");
                }

                sink.join().expect("the bgzip writer panicked")
            }
            State::Done => Err(io::Error::other("the bgzip pipeline already stopped")),
        }
    }
}

impl<T, R, S> Drop for Ordered<T, R, S> {
    fn drop(&mut self) {
        if let State::Threaded {
            jobs,
            order,
            workers,
            sink,
        } = std::mem::replace(&mut self.state, State::Done)
        {
            drop(order);
            drop(jobs);

            for worker in workers {
                worker.join().ok();
            }

            sink.join().ok();
        }
    }
}

/// A block to compress.
struct Job {
    data: Vec<u8>,
    level: CompressionLevel,
}

/// A compressed block: the frame and its uncompressed size.
type Compressed = io::Result<(Vec<u8>, usize)>;

fn compress(job: Job) -> Compressed {
    let mut frame = Vec::with_capacity(job.data.len() / 2);
    raw::deflate_frame(&job.data, job.level, &mut frame)?;
    Ok((frame, job.data.len()))
}

/// Writes frames and records where each block after the first starts.
pub struct FrameSink<W> {
    writer: W,
    compressed_offset: u64,
    uncompressed_offset: u64,
    index: Vec<(u64, u64)>,
}

impl<W> Sink<Compressed> for FrameSink<W>
where
    W: Write + Send + 'static,
{
    fn accept(&mut self, result: Compressed) -> io::Result<()> {
        let (frame, len) = result?;

        if self.compressed_offset > 0 {
            self.index
                .push((self.compressed_offset, self.uncompressed_offset));
        }

        self.writer.write_all(&frame)?;
        self.compressed_offset += frame.len() as u64;
        self.uncompressed_offset += len as u64;

        Ok(())
    }
}

/// The block boundaries of an index for `bgzf_block_write` (`bgzip -g`): the uncompressed
/// offset of each block, starting at 0.
pub struct BlockPlan {
    starts: Vec<u64>,
    current: usize,
}

impl BlockPlan {
    pub fn new(index: &[(u64, u64)]) -> Self {
        let starts = std::iter::once(0).chain(index.iter().map(|&(_, u)| u));

        Self {
            starts: starts.collect(),
            current: 0,
        }
    }

    /// The size of the current block: up to the next indexed block, or the BGZF maximum.
    fn block_size(&self) -> io::Result<usize> {
        match self.starts.get(self.current..self.current + 2) {
            Some(&[start, end]) => usize::try_from(end.saturating_sub(start))
                .ok()
                .filter(|&n| (1..=MAX_BLOCK_SIZE).contains(&n))
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid index")),
            _ => Ok(MAX_BLOCK_SIZE),
        }
    }
}

/// HTSlib's BGZF writer, as `bgzip` drives it: data is cut into blocks of [`BLOCK_SIZE`] bytes
/// unless flushed earlier (or cut as an index says), and each block is compressed on the worker
/// threads.
pub struct BlockWriter<W> {
    block: Vec<u8>,
    level: CompressionLevel,
    pipeline: Ordered<Job, Compressed, FrameSink<W>>,
}

impl<W> BlockWriter<W>
where
    W: Write + Send + 'static,
{
    /// Creates a writer with `threads` compression workers (0 compresses on the calling thread).
    pub fn new(writer: W, level: CompressionLevel, threads: usize) -> Self {
        let sink = FrameSink {
            writer,
            compressed_offset: 0,
            uncompressed_offset: 0,
            index: Vec::new(),
        };

        Self {
            block: Vec::with_capacity(BLOCK_SIZE),
            level,
            pipeline: Ordered::new(threads, compress, sink),
        }
    }

    /// `bgzf_write`.
    pub fn write(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            let n = (BLOCK_SIZE - self.block.len()).min(data.len());
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];

            if self.block.len() == BLOCK_SIZE {
                self.flush()?;
            }
        }

        Ok(())
    }

    /// `bgzf_block_write`: cuts blocks where `plan` says.
    pub fn write_planned(&mut self, mut data: &[u8], plan: &mut BlockPlan) -> io::Result<()> {
        while !data.is_empty() {
            let size = plan.block_size()?;
            let n = size.saturating_sub(self.block.len()).min(data.len());
            self.block.extend_from_slice(&data[..n]);
            data = &data[n..];

            if self.block.len() >= size {
                self.flush()?;
                plan.current += 1;
            }
        }

        Ok(())
    }

    /// Ends the current block, if it has data (`bgzf_flush`).
    pub fn flush(&mut self) -> io::Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }

        let data = std::mem::replace(&mut self.block, Vec::with_capacity(BLOCK_SIZE));

        self.pipeline.push(Job {
            data,
            level: self.level,
        })
    }

    /// Writes the last block and the end-of-file marker, and returns the writer and the index
    /// (the offsets of each block after the first, as in a `.gzi` file).
    pub fn finish(mut self) -> io::Result<(W, Vec<(u64, u64)>)> {
        self.flush()?;

        let FrameSink {
            mut writer, index, ..
        } = self.pipeline.finish()?;

        writer.write_all(&raw::EOF_MARKER)?;
        writer.flush()?;

        Ok((writer, index))
    }
}

/// A frame to decompress: its file offset and bytes.
pub type Frame = (u64, Vec<u8>);

/// A decompressed block, or the offset of the block that failed and the error.
pub type Inflated = Result<Vec<u8>, (u64, io::Error)>;

pub fn inflate((offset, frame): Frame) -> Inflated {
    let mut data = Vec::new();

    match raw::inflate_frame(&frame, &mut data) {
        Ok(()) => Ok(data),
        Err(e) => Err((offset, e)),
    }
}

/// The uncompressed size a frame's trailer gives.
pub fn frame_isize(frame: &[u8]) -> usize {
    frame
        .last_chunk::<4>()
        .map_or(0, |b| u32::from_le_bytes(*b) as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frames(src: &[u8]) -> io::Result<Vec<Vec<u8>>> {
        let mut src = src;
        let mut frames = Vec::new();
        let mut frame = Vec::new();

        while raw::read_frame(&mut src, &mut frame)? {
            frames.push(frame.clone());
        }

        Ok(frames)
    }

    fn block_sizes(src: &[u8]) -> io::Result<Vec<usize>> {
        Ok(frames(src)?.iter().map(|f| frame_isize(f)).collect())
    }

    #[test]
    fn test_block_writer() -> io::Result<()> {
        let data: Vec<u8> = (0..200_000u32).map(|i| (i % 13) as u8).collect();

        for threads in [0, 3] {
            let mut writer = BlockWriter::new(Vec::new(), CompressionLevel::default(), threads);
            writer.write(&data[..10])?;
            writer.flush()?;
            writer.write(&data[10..])?;
            let (out, index) = writer.finish()?;

            assert_eq!(block_sizes(&out)?, [10, 65280, 65280, 65280, 4150, 0]);
            assert_eq!(index.len(), 4);
            assert_eq!(index[0].1, 10);
            assert_eq!(index[3].1, 10 + 3 * 65280);

            let mut inflated = Vec::new();
            for frame in frames(&out)? {
                inflated.extend(inflate((0, frame)).map_err(|(_, e)| e)?);
            }
            assert_eq!(inflated, data);
        }

        Ok(())
    }

    #[test]
    fn test_write_planned() -> io::Result<()> {
        let mut plan = BlockPlan::new(&[(0, 100), (0, 250)]);
        let mut writer = BlockWriter::new(Vec::new(), CompressionLevel::default(), 0);
        writer.write_planned(&[b'x'; 70_000], &mut plan)?;
        let (out, _) = writer.finish()?;
        assert_eq!(block_sizes(&out)?, [100, 150, 65536, 69_750 - 65536, 0]);
        Ok(())
    }

    #[test]
    fn test_empty() -> io::Result<()> {
        let (out, index) = BlockWriter::new(Vec::new(), CompressionLevel::default(), 2).finish()?;
        assert_eq!(out, raw::EOF_MARKER);
        assert!(index.is_empty());
        Ok(())
    }
}
