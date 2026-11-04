//! Starting worker threads when the system may not start as many as asked for.
//!
//! A thread cannot start when the process reaches its thread or address-space limit (`ulimit -u`,
//! `ulimit -v`; glibc reserves 64 MB of address space for each thread's allocation arena). The
//! helpers here report that as an [`io::Error`] instead of panicking, and [`start`] retries a
//! group of workers with fewer threads, so that a large thread count runs with the threads the
//! system has room for.

use std::{
    fmt, io,
    num::NonZero,
    panic::{self, AssertUnwindSafe},
    process,
    thread::{self, JoinHandle, Scope, ScopedJoinHandle},
};

/// A group of threads that could not all be started: how many were, and why the next one was
/// not.
pub struct Shortfall {
    started: usize,
    error: io::Error,
    /// Waits for the threads that were started, once they have been told to stop.
    joins: Vec<Box<dyn FnOnce() + Send>>,
}

impl fmt::Debug for Shortfall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Shortfall")
            .field("started", &self.started)
            .field("error", &self.error)
            .finish_non_exhaustive()
    }
}

impl Shortfall {
    /// A group that stopped after `started` threads, the next one failing with `error`.
    pub fn new(started: usize, error: io::Error) -> Self {
        Self {
            started,
            error,
            joins: Vec::new(),
        }
    }

    /// A group that stopped after the threads of `handles`, the next one failing with `error`.
    /// [`start`] joins them before it tries again, so that their resources are free by then.
    pub fn with_threads<T: Send + 'static>(handles: Vec<JoinHandle<T>>, error: io::Error) -> Self {
        Self {
            started: handles.len(),
            error,
            joins: handles
                .into_iter()
                .map(|h| -> Box<dyn FnOnce() + Send> {
                    Box::new(move || {
                        let _ = h.join();
                    })
                })
                .collect(),
        }
    }

    /// The number of threads that were started.
    pub fn started(&self) -> usize {
        self.started
    }

    /// The error that stopped the next thread, worded as [`failed`] words it.
    pub fn into_error(self) -> io::Error {
        failed(self.error)
    }
}

/// Threads that could not be started, with the value that was handed in to be used by them (an
/// input or output stream), so that the caller can go on without threads.
pub struct StartError<T> {
    error: io::Error,
    value: T,
}

impl<T> StartError<T> {
    /// The failure `error` (worded by [`failed`]), returning `value`.
    pub fn new(error: io::Error, value: T) -> Self {
        Self { error, value }
    }

    /// Returns the value that was handed in.
    pub fn into_inner(self) -> T {
        self.value
    }

    /// Returns the error, dropping the value.
    pub fn into_error(self) -> io::Error {
        self.error
    }
}

impl<T> fmt::Debug for StartError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("StartError")
            .field("error", &self.error)
            .finish_non_exhaustive()
    }
}

impl<T> fmt::Display for StartError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.error.fmt(f)
    }
}

impl<T> std::error::Error for StartError<T> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.error)
    }
}

impl<T> From<StartError<T>> for io::Error {
    fn from(e: StartError<T>) -> Self {
        e.error
    }
}

/// The error for a thread that could not be started, worded as bcftools words it
/// ("Failed to create threads"), with the system's reason.
pub fn failed(error: io::Error) -> io::Error {
    io::Error::new(error.kind(), format!("Failed to create threads: {error}"))
}

/// Runs `start` with `n` threads; when it falls short, runs it again with half the threads it
/// managed to start (at least one), and so on. Halving leaves room for the threads started after
/// this group. `start` must make the threads of an attempt that falls short stop (by dropping the
/// channels they wait on, say) when it returns; the threads of its [`Shortfall`] are joined
/// before the next attempt.
///
/// # Errors
///
/// The failure of the last attempt (see [`failed`]) when not even one thread starts.
pub fn start<T>(
    n: NonZero<usize>,
    mut start: impl FnMut(NonZero<usize>) -> Result<T, Shortfall>,
) -> io::Result<T> {
    let mut count = n;
    loop {
        match start(count) {
            Ok(started) => return Ok(started),
            Err(mut shortfall) => {
                for join in shortfall.joins.drain(..) {
                    join();
                }
                let fewer = NonZero::new(shortfall.started / 2).unwrap_or(NonZero::<usize>::MIN);
                if fewer >= count {
                    return Err(shortfall.into_error());
                }
                count = fewer;
            }
        }
    }
}

/// Spawns `count` threads, the `i`th running `make(i)`.
///
/// # Errors
///
/// A [`Shortfall`] when a thread cannot be started, holding the threads started before it, which
/// must end once the caller's attempt is dropped (see [`start`]).
pub fn spawn_all<T, F>(
    count: usize,
    mut make: impl FnMut(usize) -> F,
) -> Result<Vec<JoinHandle<T>>, Shortfall>
where
    T: Send + 'static,
    F: FnOnce() -> T + Send + 'static,
{
    let mut handles = Vec::with_capacity(count);
    for i in 0..count {
        match thread::Builder::new().spawn(make(i)) {
            Ok(handle) => handles.push(handle),
            Err(e) => return Err(Shortfall::with_threads(handles, e)),
        }
    }
    Ok(handles)
}

/// A job to run on a thread of a scope.
type ScopedJob<'scope, T> = Box<dyn FnOnce() -> T + Send + 'scope>;

/// A thread started in a scope that waits for one job (see [`start_idle_scoped`]). Dropping it
/// unused ends the thread.
pub struct ScopedIdle<'scope, T> {
    job: crossbeam_channel::Sender<ScopedJob<'scope, T>>,
    handle: ScopedJoinHandle<'scope, Option<T>>,
}

impl<'scope, T> ScopedIdle<'scope, T> {
    /// Runs `f` on the thread.
    pub fn run(self, f: impl FnOnce() -> T + Send + 'scope) -> ScopedRunning<'scope, T> {
        self.job
            .send(Box::new(f))
            .expect("an idle thread waits for its job");
        ScopedRunning(self.handle)
    }
}

/// A job running on a thread of a scope (see [`ScopedIdle::run`]).
pub struct ScopedRunning<'scope, T>(ScopedJoinHandle<'scope, Option<T>>);

impl<T> ScopedRunning<'_, T> {
    /// Waits for the job and returns its result, or the panic it raised.
    ///
    /// # Errors
    ///
    /// The payload of the job's panic.
    pub fn join(self) -> thread::Result<T> {
        self.0
            .join()
            .map(|result| result.expect("a running thread was given its job"))
    }
}

/// Starts `n` threads in `scope` that each wait for one job, or fewer as [`start`] chooses them
/// when the system cannot start `n`. Starting the threads before their jobs exist lets a caller
/// size its channels and hand out its per-thread state for the threads it actually has.
///
/// # Errors
///
/// [`failed`] when not even one thread starts.
pub fn start_idle_scoped<'scope, T: Send + 'scope>(
    scope: &'scope Scope<'scope, '_>,
    n: NonZero<usize>,
) -> io::Result<Vec<ScopedIdle<'scope, T>>> {
    start(n, |count| {
        let mut idle: Vec<ScopedIdle<'scope, T>> = Vec::with_capacity(count.get());
        for started in 0..count.get() {
            let (job, jobs) = crossbeam_channel::bounded::<ScopedJob<'scope, T>>(1);
            match thread::Builder::new().spawn_scoped(scope, move || jobs.recv().ok().map(|f| f()))
            {
                Ok(handle) => idle.push(ScopedIdle { job, handle }),
                Err(e) => {
                    // The threads started so far end without a job; they are waited for here,
                    // so that their resources are free for the next attempt.
                    for t in idle {
                        drop(t.job);
                        let _ = t.handle.join();
                    }
                    return Err(Shortfall::new(started, e));
                }
            }
        }
        Ok(idle)
    })
}

/// A job run on a thread of a scope, or on the calling thread when no thread could be started
/// (see [`spawn_scoped_or_run`]).
pub enum ScopedTask<'scope, T> {
    /// The job runs on a thread.
    Spawned(ScopedJoinHandle<'scope, T>),
    /// The job ran on the calling thread, with this result.
    Done(T),
}

impl<T> ScopedTask<'_, T> {
    /// Waits for the job and returns its result, or the panic it raised on its thread.
    ///
    /// # Errors
    ///
    /// The payload of the job's panic.
    pub fn join(self) -> thread::Result<T> {
        match self {
            Self::Spawned(handle) => handle.join(),
            Self::Done(result) => Ok(result),
        }
    }
}

/// Runs `f` on a new thread in `scope`, or right away on the calling thread when the thread
/// cannot be started: for work that a thread only speeds up.
pub fn spawn_scoped_or_run<'scope, T, F>(
    scope: &'scope Scope<'scope, '_>,
    f: F,
) -> ScopedTask<'scope, T>
where
    T: Send + 'scope,
    F: FnOnce() -> T + Send + 'scope,
{
    // The job goes to the thread once it runs, so that it stays here when it cannot start.
    let (job, jobs) = crossbeam_channel::bounded::<F>(1);
    let started = thread::Builder::new().spawn_scoped(scope, move || {
        let f = jobs.recv().expect("the job is sent once the thread runs");
        f()
    });
    match started {
        Ok(handle) => {
            job.send(f).expect("the thread waits for its job");
            ScopedTask::Spawned(handle)
        }
        Err(_) => ScopedTask::Done(f()),
    }
}

/// Spawns one thread running `f`.
///
/// # Errors
///
/// [`failed`] when the thread cannot be started.
pub fn spawn<T, F>(f: F) -> io::Result<JoinHandle<T>>
where
    T: Send + 'static,
    F: FnOnce() -> T + Send + 'static,
{
    thread::Builder::new().spawn(f).map_err(failed)
}

/// A job for a [`Pool`].
type Job = Box<dyn FnOnce() + Send + 'static>;

/// A fixed group of threads that run the jobs sent to them, each job on whichever thread is free
/// first (see [`pool`]).
///
/// An idle thread sleeps until a job arrives. (A rayon pool's idle threads spin between jobs,
/// which costs as much CPU time as the work itself when the jobs are as short as inflating one
/// BGZF block.) Dropping the pool lets its threads finish the jobs already sent and end; it does
/// not wait for them.
#[derive(Debug)]
pub struct Pool {
    jobs: crossbeam_channel::Sender<Job>,
    thread_count: usize,
}

impl Pool {
    /// The number of threads.
    pub fn current_num_threads(&self) -> usize {
        self.thread_count
    }

    /// Runs `job` on a thread of the pool. A job that panics aborts the process, as on a rayon
    /// pool, so jobs catch their own panics.
    pub fn spawn(&self, job: impl FnOnce() + Send + 'static) {
        self.jobs
            .send(Box::new(job))
            .expect("the pool's threads wait for jobs while it exists");
    }
}

/// Runs the jobs sent on `jobs` until the pool is dropped.
fn run_jobs(jobs: crossbeam_channel::Receiver<Job>) {
    while let Ok(job) = jobs.recv() {
        if panic::catch_unwind(AssertUnwindSafe(job)).is_err() {
            process::abort();
        }
    }
}

/// Builds a [`Pool`] of `n` threads, or of fewer as [`start`] chooses them when the system cannot
/// start `n`.
///
/// # Errors
///
/// [`failed`] when not even one thread starts.
pub fn pool(n: NonZero<usize>) -> io::Result<Pool> {
    start(n, |count| {
        // The threads of an attempt that falls short end when `jobs` is dropped on return.
        let (jobs, rx) = crossbeam_channel::unbounded();
        let handles = spawn_all(count.get(), |_| {
            let rx = rx.clone();
            move || run_jobs(rx)
        })?;
        Ok(Pool {
            jobs,
            thread_count: handles.len(),
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn start_halves_what_could_start() {
        let mut tried = Vec::new();
        let result = start(NonZero::new(100).unwrap(), |count| {
            tried.push(count.get());
            if count.get() > 10 {
                // Up to 30 threads start.
                Err(Shortfall::new(
                    count.get().min(31) - 1,
                    io::Error::from(io::ErrorKind::WouldBlock),
                ))
            } else {
                Ok(count.get())
            }
        });
        assert_eq!(result.unwrap(), 7);
        assert_eq!(tried, [100, 15, 7]);
    }

    #[test]
    fn start_fails_when_no_thread_starts() {
        let mut attempts = 0;
        let result: io::Result<()> = start(NonZero::new(8).unwrap(), |_| {
            attempts += 1;
            Err(Shortfall::new(
                0,
                io::Error::from(io::ErrorKind::WouldBlock),
            ))
        });
        let e = result.unwrap_err();
        assert!(e.to_string().starts_with("Failed to create threads"), "{e}");
        // 8 threads, then one.
        assert_eq!(attempts, 2);
    }

    #[test]
    fn idle_threads_run_their_jobs() {
        thread::scope(|scope| {
            let mut idle = start_idle_scoped(scope, NonZero::new(3).unwrap()).unwrap();
            assert_eq!(idle.len(), 3);
            // An idle thread left without a job ends.
            drop(idle.pop());
            let running: Vec<_> = idle
                .into_iter()
                .enumerate()
                .map(|(i, t)| t.run(move || i * 10))
                .collect();
            let results: Vec<usize> = running.into_iter().map(|r| r.join().unwrap()).collect();
            assert_eq!(results, [0, 10]);
        });
    }

    #[test]
    fn a_task_runs_on_a_thread_or_here() {
        let mut data = [1, 2, 3];
        thread::scope(|scope| {
            let (a, b) = data.split_at_mut(1);
            let first = spawn_scoped_or_run(scope, move || {
                a[0] += 1;
                a[0]
            });
            let second = spawn_scoped_or_run(scope, move || b.iter().sum::<i32>());
            assert_eq!(first.join().unwrap(), 2);
            assert_eq!(second.join().unwrap(), 5);
        });
        assert_eq!(data, [2, 2, 3]);
    }

    #[test]
    fn pool_runs_jobs() {
        let pool = pool(NonZero::new(3).unwrap()).unwrap();
        assert_eq!(pool.current_num_threads(), 3);
        let (tx, rx) = crossbeam_channel::unbounded();
        for i in 0..10 {
            let tx = tx.clone();
            pool.spawn(move || tx.send(i).unwrap());
        }
        drop(tx);
        let mut results: Vec<i32> = rx.iter().collect();
        results.sort_unstable();
        assert_eq!(results, (0..10).collect::<Vec<_>>());
    }
}
