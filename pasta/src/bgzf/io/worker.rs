//! Hands the panic of a compression or decompression job to the thread waiting for its result.
//!
//! A job that panics never sends its result. Its peers keep the result channel open, so the
//! waiting thread would wait forever (and a panic on a pool thread aborts the process). Jobs run
//! under [`catch`] and send what it returns; the waiting thread re-raises a panic with
//! [`std::panic::resume_unwind`].

use std::{
    panic::{self, AssertUnwindSafe},
    thread,
};

/// Runs a job, returning its panic instead of unwinding. The job's state must not be used after
/// a panic: the panic is re-raised on the thread that waits for the result.
pub(crate) fn catch<T>(job: impl FnOnce() -> T) -> thread::Result<T> {
    panic::catch_unwind(AssertUnwindSafe(job))
}

/// The bytes that make a job panic in test builds (see [`panic_if_marked`]).
#[cfg(test)]
pub(crate) const PANIC_MARKER: &[u8] = b"pasta-test-worker-panic";

/// Panics in test builds when `data` holds [`PANIC_MARKER`], so that tests can make a real job
/// panic; does nothing in other builds.
#[inline(always)]
pub(crate) fn panic_if_marked(data: &[u8]) {
    #[cfg(test)]
    assert!(
        memchr::memmem::find(data, PANIC_MARKER).is_none(),
        "{MARKER_PANIC}"
    );
    #[cfg(not(test))]
    let _ = data;
}

/// Runs `f` on its own thread and returns its result, failing (rather than hanging) when it
/// takes more than 30 seconds.
#[cfg(test)]
#[track_caller]
pub(crate) fn within_timeout<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
    use std::{sync::mpsc, time::Duration};

    let (tx, rx) = mpsc::channel();
    thread::spawn(move || tx.send(f()));
    match rx.recv_timeout(Duration::from_secs(30)) {
        Ok(value) => value,
        Err(mpsc::RecvTimeoutError::Timeout) => panic!("hung for 30 seconds"),
        Err(mpsc::RecvTimeoutError::Disconnected) => panic!("panicked"),
    }
}

/// The message of the panic [`panic_if_marked`] raises.
#[cfg(test)]
pub(crate) const MARKER_PANIC: &str = "a job met the test panic marker";

/// Runs `f`, which must panic with `message`, and fails otherwise.
#[cfg(test)]
#[track_caller]
pub(crate) fn assert_panics_with<T: std::fmt::Debug>(message: &str, f: impl FnOnce() -> T) {
    match catch(f) {
        Ok(value) => panic!("returned {value:?} instead of panicking"),
        Err(panic) => {
            let raised = (panic.downcast_ref::<&str>().copied())
                .or_else(|| panic.downcast_ref::<String>().map(String::as_str));
            assert_eq!(raised, Some(message));
        }
    }
}
