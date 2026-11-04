use std::{
    collections::VecDeque,
    io::{self, BufRead, Read, Seek, SeekFrom},
    mem,
    num::NonZero,
    panic,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};

use super::{Block, reader::frame::block_initialize, worker};
use crate::{
    bgzf::{VirtualPosition, gzi},
    util::threads,
};

/// A decompressed block or a read error, or the panic of the job or of the reader thread that
/// produced it.
type Buffered = thread::Result<io::Result<Buffer>>;
type BufferedRx = Receiver<Buffered>;
type ReadTx = Sender<BufferedRx>;
type ReadRx = Receiver<BufferedRx>;
type RecycleTx = Sender<Buffer>;
type RecycleRx = Receiver<Buffer>;

const NON_WORKER_COUNT: usize = 2;

/// Blocks in flight per worker.
const BLOCKS_PER_WORKER: usize = 4;

enum State<R> {
    Paused((R, threads::Pool)),
    Running {
        reader_handle: JoinHandle<(R, threads::Pool)>,
        read_rx: ReadRx,
        recycle_tx: RecycleTx,
    },
    Done,
}

#[derive(Debug, Default)]
struct Buffer {
    buf: Vec<u8>,
    block: Block,
}

/// A multithreaded BGZF reader.
///
/// This is a multithreaded BGZF reader that uses a thread pool to decompress block data. It places
/// the inner reader on its own thread to read raw frames asynchronously.
///
/// A panic of a decompression job or of the inner reader is re-raised on the calling thread when
/// it reaches that block, or by [`Self::finish`].
pub struct MultithreadedReader<R> {
    state: State<R>,
    worker_count: NonZero<usize>,
    position: u64,
    // The end of the last block with data (or the position sought to): where the stream's
    // virtual position stays at its end, before any empty blocks such as the EOF marker.
    data_end: u64,
    buffer: Buffer,
    /// The blocks the reader thread had read ahead when it was paused, which are read before the
    /// inner reader is read again.
    read_ahead: VecDeque<Buffered>,
    /// Whether the reader thread could not be started, so that blocks are read and inflated on
    /// the calling thread from then on.
    serial: bool,
}

impl<R> MultithreadedReader<R> {
    /// Returns the current position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// assert_eq!(reader.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the current virtual position of the stream.
    ///
    /// At the end of the stream, this is the end of the last block with data, before any empty
    /// blocks (such as the EOF marker) that follow it, as HTSlib's `bgzf_tell` gives.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::MIN);
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.buffer.block.virtual_position()
    }

    /// Shuts down the reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let mut reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// reader.finish()?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn finish(&mut self) -> io::Result<R> {
        let state = mem::replace(&mut self.state, State::Done);

        // Blocks read ahead, and the errors among them, are dropped as blocks still in flight are.
        let read_ahead = mem::take(&mut self.read_ahead);

        match state {
            State::Paused((inner, _)) => match read_ahead_error(read_ahead) {
                Some(e) => Err(e),
                None => Ok(inner),
            },
            State::Running {
                reader_handle,
                read_rx,
                recycle_tx,
            } => {
                drop(recycle_tx);

                let err = read_rx_drain(&read_rx);
                let (inner, _) = reader_handle
                    .join()
                    .unwrap_or_else(|p| panic::resume_unwind(p));

                match err.unwrap_or_else(|p| panic::resume_unwind(p)) {
                    Some(e) => Err(e),
                    None => Ok(inner),
                }
            }
            State::Done => panic!("invalid state"),
        }
    }
}

impl<R> MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    /// Creates a multithreaded BGZF reader with a worker count of 1.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self::with_worker_count(NonZero::<usize>::MIN, inner)
    }

    /// Creates a multithreaded BGZF reader with a worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    /// use pasta::bgzf;
    /// let reader = bgzf::io::MultithreadedReader::with_worker_count(
    ///     NonZero::<usize>::MIN,
    ///     io::empty(),
    /// );
    /// ```
    ///
    /// # Panics
    ///
    /// Panics when not even one worker thread can be started; [`Self::try_with_worker_count`]
    /// returns that error.
    pub fn with_worker_count(worker_count: NonZero<usize>, inner: R) -> Self {
        Self::try_with_worker_count(worker_count, inner)
            .map_err(threads::StartError::into_error)
            .expect("a worker thread starts")
    }

    /// Creates a multithreaded BGZF reader with a worker count, or with fewer workers when the
    /// system cannot start that many threads (see [`crate::util::threads::pool`]).
    ///
    /// # Errors
    ///
    /// Fails when not even one worker thread can be started, returning `inner` in the error.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    /// use pasta::bgzf;
    /// let reader = bgzf::io::MultithreadedReader::try_with_worker_count(
    ///     NonZero::<usize>::MIN,
    ///     io::empty(),
    /// )?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_with_worker_count(
        worker_count: NonZero<usize>,
        inner: R,
    ) -> Result<Self, threads::StartError<R>> {
        let pool = match threads::pool(worker_count) {
            Ok(pool) => pool,
            Err(e) => return Err(threads::StartError::new(e, inner)),
        };
        let worker_count = NonZero::new(pool.current_num_threads()).unwrap_or(worker_count);

        Ok(Self {
            state: State::Paused((inner, pool)),
            worker_count,
            position: 0,
            data_end: 0,
            buffer: Buffer::default(),
            read_ahead: VecDeque::new(),
            serial: false,
        })
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// This stops the reader thread. The blocks it had already read from the underlying reader
    /// are kept and are read first, so that reading on loses nothing; the underlying reader is
    /// then read from wherever it is. Its position is thus past the blocks that were read ahead,
    /// not at [`Self::position`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    /// let mut reader = bgzf::io::MultithreadedReader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        self.pause();

        match &mut self.state {
            State::Paused((inner, _)) => inner,
            _ => panic!("invalid state"),
        }
    }

    /// Starts the reader thread, unless it runs.
    ///
    /// # Errors
    ///
    /// Fails, leaving the reader paused, when the thread cannot be started.
    fn resume(&mut self) -> io::Result<()> {
        if matches!(self.state, State::Running { .. }) {
            return Ok(());
        }

        let state = mem::replace(&mut self.state, State::Done);

        let State::Paused((inner, pool)) = state else {
            panic!("invalid state");
        };

        let worker_count = self.worker_count.get();
        let buffer_count = worker_count * BLOCKS_PER_WORKER + NON_WORKER_COUNT;

        let (read_tx, read_rx) = crossbeam_channel::bounded(buffer_count);
        let (recycle_tx, recycle_rx) = crossbeam_channel::bounded(buffer_count);

        for _ in 0..buffer_count {
            recycle_tx.send(Buffer::default()).unwrap();
        }

        // The reader and the pool go to the thread once it runs, so that they stay here when it
        // cannot be started.
        let (start_tx, start_rx) = crossbeam_channel::bounded(1);
        let reader_handle = match spawn_reader(start_rx, read_tx, recycle_rx) {
            Ok(handle) => handle,
            Err(e) => {
                self.state = State::Paused((inner, pool));
                return Err(e);
            }
        };
        start_tx
            .send((inner, pool))
            .expect("the reader thread waits for its reader");

        self.state = State::Running {
            reader_handle,
            read_rx,
            recycle_tx,
        };

        Ok(())
    }

    fn pause(&mut self) {
        if matches!(self.state, State::Paused(_)) {
            return;
        }

        let state = mem::replace(&mut self.state, State::Done);

        let State::Running {
            reader_handle,
            read_rx,
            recycle_tx,
        } = state
        else {
            panic!("invalid state");
        };

        drop(recycle_tx);

        // The blocks read ahead (and their errors and panics) are kept, to be read in order.
        while let Ok(buffered_rx) = read_rx.recv() {
            let buffered = buffered_rx
                .recv()
                .unwrap_or_else(|_| Ok(Err(job_ended_error())));
            self.read_ahead.push_back(buffered);
        }

        let (inner, pool) = reader_handle
            .join()
            .unwrap_or_else(|p| panic::resume_unwind(p));

        self.state = State::Paused((inner, pool));
    }

    fn read_block(&mut self) -> io::Result<()> {
        loop {
            let next = if let Some(buffered) = self.read_ahead.pop_front() {
                unwrap_buffered(buffered).map(Some)?
            } else {
                // Without its reader thread, which the system could not start, the reader reads
                // and inflates each block here.
                if self.serial || self.resume().is_err() {
                    self.serial = true;
                    return self.read_block_here();
                }

                let State::Running { read_rx, .. } = &self.state else {
                    panic!("invalid state");
                };

                recv_buffer(read_rx)?
            };

            if let Some(mut buffer) = next {
                buffer.block.set_position(self.position);
                self.position += buffer.block.size();

                let prev_buffer = mem::replace(&mut self.buffer, buffer);
                if let State::Running { recycle_tx, .. } = &self.state {
                    recycle_tx.send(prev_buffer).ok();
                }

                if self.buffer.block.data().len() > 0 {
                    self.data_end = self.position;
                    break;
                }
            } else {
                block_initialize(&mut self.buffer.block, 0, 0);
                self.buffer.block.set_position(self.data_end);
                break;
            }
        }

        Ok(())
    }
}

impl<R> MultithreadedReader<R>
where
    R: Read,
{
    /// Reads and inflates the next block with data on this thread, as the reader thread and the
    /// pool would.
    fn read_block_here(&mut self) -> io::Result<()> {
        use super::reader::frame::{parse_block, read_frame_into};

        let State::Paused((inner, _)) = &mut self.state else {
            panic!("invalid state");
        };

        loop {
            // The current block is used up: its buffers take the next one.
            let mut buffer = mem::take(&mut self.buffer);
            let read = read_frame_into(inner, &mut buffer.buf).and_then(|frame| {
                frame
                    .map(|()| parse_block(&buffer.buf, &mut buffer.block))
                    .transpose()
            });
            match read {
                Ok(Some(())) => {
                    buffer.block.set_position(self.position);
                    self.position += buffer.block.size();
                    self.buffer = buffer;

                    if self.buffer.block.data().len() > 0 {
                        self.data_end = self.position;
                        return Ok(());
                    }
                }
                Ok(None) => {
                    self.buffer = buffer;
                    block_initialize(&mut self.buffer.block, 0, 0);
                    self.buffer.block.set_position(self.data_end);
                    return Ok(());
                }
                Err(e) => {
                    self.buffer = buffer;
                    return Err(e);
                }
            }
        }
    }
}

impl<R> Drop for MultithreadedReader<R> {
    fn drop(&mut self) {
        if !matches!(self.state, State::Done) {
            // A panic of a worker is not raised again while the thread is already unwinding.
            if thread::panicking() {
                let _ = worker::catch(|| self.finish());
            } else {
                let _ = self.finish();
            }
        }
    }
}

impl<R> Read for MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        use super::reader::default_read_exact;

        if let Some(src) = self.buffer.block.data().as_ref().get(..buf.len()) {
            buf.copy_from_slice(src);
            self.consume(src.len());
            Ok(())
        } else {
            default_read_exact(self, buf)
        }
    }
}

impl<R> BufRead for MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.buffer.block.data().has_remaining() {
            self.read_block()?;
        }

        Ok(self.buffer.block.data().as_ref())
    }

    fn consume(&mut self, amt: usize) {
        self.buffer.block.data_mut().consume(amt);
    }
}

impl<R> crate::bgzf::io::Read for MultithreadedReader<R>
where
    R: Read + Send + 'static,
{
    fn virtual_position(&self) -> VirtualPosition {
        self.buffer.block.virtual_position()
    }
}

impl<R> crate::bgzf::io::BufRead for MultithreadedReader<R> where R: Read + Send + 'static {}

impl<R> crate::bgzf::io::Seek for MultithreadedReader<R>
where
    R: Read + Send + Seek + 'static,
{
    fn seek_to_virtual_position(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        let (cpos, upos) = pos.into();

        self.get_mut().seek(SeekFrom::Start(cpos))?;
        self.read_ahead.clear();
        self.position = cpos;
        self.data_end = cpos;

        if let Err(e) = self.read_block() {
            block_initialize(&mut self.buffer.block, 0, 0);
            self.buffer.block.set_position(cpos);
            return Err(e);
        }

        let data = self.buffer.block.data_mut();

        if usize::from(upos) <= data.len() {
            data.set_position(usize::from(upos));
        } else {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        Ok(pos)
    }

    fn seek_with_index(&mut self, index: &gzi::Index, pos: SeekFrom) -> io::Result<u64> {
        use crate::bgzf::io::seek::resolve_with_index;

        let pos = resolve_with_index(index, self.virtual_position(), pos)?;
        let virtual_position = index.query(pos)?;
        self.seek_to_virtual_position(virtual_position)?;
        Ok(pos)
    }
}

/// Returns the next block, or `None` once the reader thread has read the last one.
///
/// # Panics
///
/// Re-raises the panic of the job that decompressed the block or of the reader thread.
fn recv_buffer(read_rx: &ReadRx) -> io::Result<Option<Buffer>> {
    let Ok(buffered_rx) = read_rx.recv() else {
        return Ok(None);
    };

    // Every job sends its result, a panic included.
    let buffered = buffered_rx.recv().map_err(|_| job_ended_error())?;

    unwrap_buffered(buffered).map(Some)
}

/// Returns the block, or its read error.
///
/// # Panics
///
/// Re-raises the panic of the job that decompressed the block or of the reader thread.
fn unwrap_buffered(buffered: Buffered) -> io::Result<Buffer> {
    buffered.unwrap_or_else(|p| panic::resume_unwind(p))
}

fn job_ended_error() -> io::Error {
    io::Error::other("BGZF decompression job ended without a result")
}

/// Returns the first read error of blocks read ahead.
///
/// # Panics
///
/// Re-raises the first panic among them.
fn read_ahead_error(read_ahead: VecDeque<Buffered>) -> Option<io::Error> {
    let mut err = None;

    for buffered in read_ahead {
        match buffered {
            Ok(Err(e)) => {
                err.get_or_insert(e);
            }
            Err(p) => panic::resume_unwind(p),
            Ok(Ok(_)) => {}
        }
    }

    err
}

/// Starts the thread that reads frames and hands them to the pool, once it is sent the reader and
/// the pool on `start_rx`.
fn spawn_reader<R>(
    start_rx: Receiver<(R, threads::Pool)>,
    read_tx: ReadTx,
    recycle_rx: RecycleRx,
) -> io::Result<JoinHandle<(R, threads::Pool)>>
where
    R: Read + Send + 'static,
{
    use super::reader::frame::{parse_block, read_frame_into};

    threads::spawn(move || {
        let (mut reader, pool) = start_rx
            .recv()
            .expect("the reader and the pool are sent once the thread runs");
        while let Ok(mut buffer) = recycle_rx.recv() {
            // A panic of the inner reader is sent like a read error, and ends reading.
            let error = match worker::catch(|| read_frame_into(&mut reader, &mut buffer.buf)) {
                Ok(Ok(Some(()))) => None,
                Ok(Ok(None)) => break,
                Ok(Err(e)) => Some(Ok(Err(e))),
                Err(p) => Some(Err(p)),
            };

            if let Some(error) = error {
                let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);
                let _ = buffered_tx.send(error);
                let _ = read_tx.send(buffered_rx);
                break;
            }

            let (buffered_tx, buffered_rx) = crossbeam_channel::bounded(1);

            pool.spawn(move || {
                let result = worker::catch(|| {
                    parse_block(&buffer.buf, &mut buffer.block)?;
                    worker::panic_if_marked(buffer.block.data().as_ref());
                    Ok(buffer)
                });
                let _ = buffered_tx.send(result);
            });

            if read_tx.send(buffered_rx).is_err() {
                break;
            }
        }

        (reader, pool)
    })
}

/// Waits for the blocks still being read and returns the first read error, or the first panic of
/// a job or of the reader thread.
fn read_rx_drain(read_rx: &ReadRx) -> thread::Result<Option<io::Error>> {
    let mut err = None;
    let mut panic = None;

    while let Ok(buffered_rx) = read_rx.recv() {
        match buffered_rx.recv() {
            Ok(Ok(Err(e))) => {
                err.get_or_insert(e);
            }
            Ok(Err(p)) => {
                panic.get_or_insert(p);
            }
            Ok(Ok(Ok(_))) | Err(_) => {}
        }
    }

    match panic {
        Some(p) => Err(p),
        None => Ok(err),
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Write};

    use super::*;
    use crate::bgzf::io::{
        Seek,
        worker::{MARKER_PANIC, PANIC_MARKER, assert_panics_with, within_timeout},
        writer::BGZF_EOF,
    };

    /// A reader of `data` that panics once it has given `left` bytes.
    struct PanicAfter {
        data: Cursor<Vec<u8>>,
        left: usize,
    }

    impl Read for PanicAfter {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.data.read(buf)?;
            self.left = (self.left.checked_sub(n)).expect("the reader met its limit");
            Ok(n)
        }
    }

    /// Returns BGZF data of many blocks, one of which holds the test panic marker.
    fn marked_data() -> io::Result<Vec<u8>> {
        let mut writer = crate::bgzf::io::Writer::new(Vec::new());
        writer.write_all(&vec![b'x'; 1 << 17])?;
        writer.flush()?;
        writer.write_all(PANIC_MARKER)?;
        writer.flush()?;
        writer.write_all(&vec![b'x'; 1 << 22])?;
        writer.finish()
    }

    #[test]
    fn a_panic_on_a_thread_is_raised_again() -> io::Result<()> {
        let data = marked_data()?;
        for workers in [1, 2, 4] {
            let workers = NonZero::new(workers).unwrap();
            // A decompression job's panic (which used to abort the process).
            let marked = data.clone();
            within_timeout(move || {
                assert_panics_with(MARKER_PANIC, || {
                    let inner = Cursor::new(marked);
                    let mut reader = MultithreadedReader::with_worker_count(workers, inner);
                    reader.read_to_end(&mut Vec::new())
                });
            });
            // The inner reader's panic (which used to end the stream early, as if it were the
            // end of the input).
            let mut plain = crate::bgzf::io::Writer::new(Vec::new());
            plain.write_all(&vec![b'x'; 1 << 22])?;
            let plain = plain.finish()?;
            within_timeout(move || {
                assert_panics_with("the reader met its limit", || {
                    let inner = PanicAfter {
                        left: plain.len() / 2,
                        data: Cursor::new(plain),
                    };
                    let mut reader = MultithreadedReader::with_worker_count(workers, inner);
                    reader.read_to_end(&mut Vec::new())
                });
            });
        }
        Ok(())
    }

    #[test]
    fn get_mut_loses_no_blocks_read_ahead() -> io::Result<()> {
        // Many blocks, each of distinct data.
        let data: Vec<u8> = (0..(1u32 << 20)).flat_map(u32::to_le_bytes).collect();
        let mut writer = crate::bgzf::io::Writer::new(Vec::new());
        writer.write_all(&data)?;
        let compressed = writer.finish()?;

        for workers in [1, 2, 4] {
            let workers = NonZero::new(workers).unwrap();
            let mut reader =
                MultithreadedReader::with_worker_count(workers, Cursor::new(compressed.clone()));

            let mut buf = vec![0; 100_000];
            reader.read_exact(&mut buf)?;

            // As a caller that reports its progress does, between reads.
            let inner_position = reader.get_mut().position();
            assert!(inner_position > reader.position());
            reader.get_mut();

            let mut rest = Vec::new();
            reader.read_to_end(&mut rest)?;
            buf.extend(rest);
            assert!(
                buf == data,
                "with {workers} workers, the data is read whole"
            );
            assert_eq!(reader.finish()?.position(), compressed.len() as u64);
        }

        Ok(())
    }

    #[test]
    fn test_fill_buf_with_frame_reader_error() {
        static DATA: &[u8] = &[
            // EOF block with invalid BSIZE (255)
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0xff, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = MultithreadedReader::new(Cursor::new(DATA));

        assert!(matches!(
            reader.fill_buf(),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
    }

    #[test]
    fn test_seek_to_virtual_position_with_uncompressed_position_gt_data_len()
    -> Result<(), crate::bgzf::virtual_position::TryFromU64U16TupleError> {
        #[rustfmt::skip]
        static DATA: &[u8] = &[
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = MultithreadedReader::new(Cursor::new(DATA));

        assert!(matches!(
            reader.seek_to_virtual_position(VirtualPosition::try_from((0, 8))?),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_seek_to_virtual_position() -> Result<(), Box<dyn std::error::Error>> {
        use crate::bgzf::io::Seek;

        #[rustfmt::skip]
        static DATA: &[u8] = &[
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        // HTSlib's `bgzf_tell` at the end: before the EOF block.
        const EOF_VIRTUAL_POSITION: VirtualPosition = VirtualPosition::new(35, 0).unwrap();

        let mut reader = MultithreadedReader::new(Cursor::new(DATA));

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(reader.virtual_position(), EOF_VIRTUAL_POSITION);

        reader.seek_to_virtual_position(const { VirtualPosition::new(0, 3).unwrap() })?;

        buf.clear();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"dles");
        assert_eq!(reader.virtual_position(), EOF_VIRTUAL_POSITION);

        Ok(())
    }

    #[test]
    fn test_virtual_position_at_end_is_before_empty_blocks()
    -> Result<(), Box<dyn std::error::Error>> {
        #[rustfmt::skip]
        const BLOCK: [u8; 35] = [
            // b"noodles"
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
        ];
        // A block, an empty block (at 35), a block (at 63), and two empty blocks (at 98).
        let data = [&BLOCK[..], &BGZF_EOF, &BLOCK, &BGZF_EOF, &BGZF_EOF].concat();
        let end = VirtualPosition::new(98, 0).unwrap();

        let mut reader = MultithreadedReader::new(Cursor::new(data.clone()));
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;
        assert_eq!(buf, b"noodlesnoodles");
        assert_eq!(reader.virtual_position(), end);

        reader.seek_to_virtual_position(VirtualPosition::new(35, 0).unwrap())?;
        buf.clear();
        reader.read_to_end(&mut buf)?;
        assert_eq!(buf, b"noodles");
        assert_eq!(reader.virtual_position(), end);

        // The same without the reader thread (as when it cannot be started).
        let mut reader = MultithreadedReader::new(Cursor::new(data));
        reader.serial = true;
        buf.clear();
        reader.read_to_end(&mut buf)?;
        assert_eq!(buf, b"noodlesnoodles");
        assert_eq!(reader.virtual_position(), end);

        reader.seek_to_virtual_position(VirtualPosition::new(35, 0).unwrap())?;
        buf.clear();
        reader.read_to_end(&mut buf)?;
        assert_eq!(buf, b"noodles");
        assert_eq!(reader.virtual_position(), end);
        assert!(
            matches!(reader.state, State::Paused(_)),
            "no thread was started"
        );

        Ok(())
    }
}
