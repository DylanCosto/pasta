//! Helpers for worker threads that send their results to a main thread over a channel.

use std::thread;

use crossbeam_channel::Sender;

/// Sends a message on a channel when the thread holding it panics, so that a main thread
/// waiting for that worker's result stops instead of waiting forever: the other workers keep
/// the channel open, so it never disconnects. The main thread then re-raises the panic once it
/// has joined the workers.
pub(crate) struct ReportPanic<'a, T, F: FnMut() -> T> {
    tx: &'a Sender<T>,
    message: F,
}

impl<'a, T, F: FnMut() -> T> ReportPanic<'a, T, F> {
    /// Sends `message()` on `tx` if the thread panics while this is held.
    pub(crate) fn new(tx: &'a Sender<T>, message: F) -> Self {
        Self { tx, message }
    }
}

impl<T, F: FnMut() -> T> Drop for ReportPanic<'_, T, F> {
    fn drop(&mut self) {
        if thread::panicking() {
            let _ = self.tx.send((self.message)());
        }
    }
}

/// The text that makes a worker panic in test builds (see [`panic_if_marked`]).
#[cfg(test)]
pub(crate) const PANIC_MARKER: &str = "vtools-test-worker-panic";

/// Panics in test builds when `input` holds [`PANIC_MARKER`], so that tests can make a real
/// worker panic; does nothing in other builds.
#[inline(always)]
pub(crate) fn panic_if_marked(input: &[u8]) {
    #[cfg(test)]
    assert!(
        memchr::memmem::find(input, PANIC_MARKER.as_bytes()).is_none(),
        "a worker met the test panic marker"
    );
    #[cfg(not(test))]
    let _ = input;
}

/// Joins `workers`, then re-raises the first panic among them, or panics with `what` if none
/// panicked. For a main thread told by [`ReportPanic`] that a worker panicked; the channel of
/// jobs must be closed first, so that the other workers stop.
pub(crate) fn rethrow(workers: impl IntoIterator<Item = thread::JoinHandle<()>>, what: &str) -> ! {
    let mut first = None;
    for w in workers {
        if let Err(panic) = w.join() {
            first.get_or_insert(panic);
        }
    }
    match first {
        Some(panic) => std::panic::resume_unwind(panic),
        None => panic!("{what}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_panic_is_reported_once_and_a_normal_end_is_not() {
        let (tx, rx) = crossbeam_channel::unbounded();
        let fine = {
            let tx = tx.clone();
            thread::spawn(move || {
                let _guard = ReportPanic::new(&tx, || "panicked");
            })
        };
        assert!(fine.join().is_ok());
        let failed = {
            let tx = tx.clone();
            thread::spawn(move || {
                let _guard = ReportPanic::new(&tx, || "panicked");
                panic_if_marked(PANIC_MARKER.as_bytes());
            })
        };
        assert!(failed.join().is_err());
        drop(tx);
        assert_eq!(rx.iter().collect::<Vec<_>>(), ["panicked"]);
    }

    #[test]
    fn rethrow_raises_the_workers_panic() {
        let workers = vec![
            thread::spawn(|| {}),
            thread::spawn(|| std::panic::panic_any(7_u32)),
        ];
        let raised = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            rethrow(workers, "none panicked")
        }))
        .expect_err("rethrow always panics");
        assert_eq!(raised.downcast_ref::<u32>(), Some(&7));
    }
}
