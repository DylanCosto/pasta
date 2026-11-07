use std::{io::Write, num::NonZero};

use bytes::BytesMut;

use super::MultithreadedWriter;
use crate::{bgzf::io::writer::CompressionLevel, util::threads};

/// A multithreaded BGZF writer builder.
pub struct Builder {
    compression_level: CompressionLevel,
    worker_count: NonZero<usize>,
}

impl Builder {
    /// Sets the compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::bgzf::{self as bgzf, io::writer::CompressionLevel};
    ///
    /// let builder = bgzf::io::multithreaded_writer::Builder::default()
    ///     .set_compression_level(CompressionLevel::default());
    /// ```
    pub fn set_compression_level(mut self, compression_level: CompressionLevel) -> Self {
        self.compression_level = compression_level;
        self
    }

    /// Sets the worker count.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    /// use pasta::bgzf;;
    ///
    /// let builder = bgzf::io::multithreaded_writer::Builder::default()
    ///     .set_worker_count(NonZero::<usize>::MIN);
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZero<usize>) -> Self {
        self.worker_count = worker_count;
        self
    }

    /// Builds a multithreaded BGZF writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    ///
    /// let builder = bgzf::io::multithreaded_writer::Builder::default()
    ///     .build_from_writer(io::sink());
    /// ```
    ///
    /// # Panics
    ///
    /// Panics when the threads cannot be started; [`Self::try_build_from_writer`] returns that
    /// error.
    pub fn build_from_writer<W>(self, writer: W) -> MultithreadedWriter<W>
    where
        W: Write + Send + 'static,
    {
        self.try_build_from_writer(writer)
            .map_err(threads::StartError::into_error)
            .expect("the writer threads start")
    }

    /// Builds a multithreaded BGZF writer from a writer, with fewer workers than set when the
    /// system cannot start that many threads (see [`crate::util::threads::pool`]).
    ///
    /// # Errors
    ///
    /// Fails when not even one worker thread, or the thread that writes the compressed blocks,
    /// can be started, returning `writer` in the error.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::bgzf;
    ///
    /// let writer = bgzf::io::multithreaded_writer::Builder::default()
    ///     .try_build_from_writer(io::sink())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_build_from_writer<W>(
        self,
        writer: W,
    ) -> Result<MultithreadedWriter<W>, threads::StartError<W>>
    where
        W: Write + Send + 'static,
    {
        use super::{State, spawn_writer};

        let pool = match threads::pool(self.worker_count) {
            Ok(pool) => pool,
            Err(e) => return Err(threads::StartError::new(e, writer)),
        };
        let worker_count = pool.current_num_threads();

        let (write_tx, write_rx) = crossbeam_channel::bounded(worker_count);
        // The writer goes to its thread once it runs, so that it is returned when the thread
        // cannot be started.
        let (start_tx, start_rx) = crossbeam_channel::bounded(1);
        let (flushed_tx, flushed_rx) = crossbeam_channel::bounded(1);
        let writer_handle = match spawn_writer(start_rx, write_rx, flushed_tx) {
            Ok(handle) => handle,
            Err(e) => return Err(threads::StartError::new(e, writer)),
        };
        start_tx
            .send(writer)
            .expect("the writer thread waits for its writer");

        Ok(MultithreadedWriter {
            pool,
            state: State::Running {
                writer_handle,
                write_tx,
                flushed_rx,
            },
            buf: BytesMut::new(),
            compression_level: self.compression_level,
        })
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            compression_level: CompressionLevel::default(),
            worker_count: NonZero::<usize>::MIN,
        }
    }
}
