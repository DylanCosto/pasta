//! One budget of CPU threads for the threads that pile up and the threads that realign the
//! reads of heavy indel positions, so that together they keep a bounded number of cores busy
//! (twice `--threads` with worker threads, `--threads` in a serial run).

use std::sync::{Arc, Condvar, Mutex};

/// CPU threads shared between the threads that pile up and the threads they borrow for heavy
/// indel positions. Clones share the budget.
#[derive(Clone, Debug)]
pub struct ThreadBudget(Arc<Pool>);

#[derive(Debug)]
struct Pool {
    free: Mutex<usize>,
    freed: Condvar,
}

impl ThreadBudget {
    /// A budget of `threads` threads (at least one, the thread that piles up).
    pub fn new(threads: usize) -> Self {
        Self(Arc::new(Pool {
            free: Mutex::new(threads.max(1)),
            freed: Condvar::new(),
        }))
    }

    /// Waits for a thread of the budget, which is held until the guard is dropped.
    pub fn hold(&self) -> Held<'_> {
        let pool = &*self.0;
        let mut free = pool.free.lock().expect("thread budget lock poisoned");
        while *free == 0 {
            free = pool.freed.wait(free).expect("thread budget lock poisoned");
        }
        *free -= 1;
        Held { pool, count: 1 }
    }

    /// Takes up to `most` threads of the budget, as many as are free now, without waiting.
    pub fn take_up_to(&self, most: usize) -> Held<'_> {
        let pool = &*self.0;
        let mut free = pool.free.lock().expect("thread budget lock poisoned");
        let count = most.min(*free);
        *free -= count;
        Held { pool, count }
    }
}

/// Threads taken from a [`ThreadBudget`], given back when dropped.
#[must_use = "the threads are given back when the guard is dropped"]
#[derive(Debug)]
pub struct Held<'a> {
    pool: &'a Pool,
    count: usize,
}

impl Held<'_> {
    /// How many threads are held.
    pub fn count(&self) -> usize {
        self.count
    }
}

impl Drop for Held<'_> {
    fn drop(&mut self) {
        if self.count == 0 {
            return;
        }
        // Lock poisoning cannot happen here: no code panics while holding the lock.
        if let Ok(mut free) = self.pool.free.lock() {
            *free += self.count;
            self.pool.freed.notify_all();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn borrowed_threads_never_exceed_the_budget() {
        let budget = ThreadBudget::new(4);
        let worker = budget.hold();
        let extra = budget.take_up_to(8);
        assert_eq!(extra.count(), 3, "the worker holds one of four");
        assert_eq!(budget.take_up_to(1).count(), 0);
        drop(extra);
        assert_eq!(budget.take_up_to(2).count(), 2);
        drop(worker);
        assert_eq!(budget.take_up_to(8).count(), 4);
    }

    #[test]
    fn no_threads_leaves_one_for_the_pileup_and_none_to_borrow() {
        let budget = ThreadBudget::new(0);
        let _main = budget.hold();
        assert_eq!(budget.take_up_to(4).count(), 0);
    }

    #[test]
    fn a_worker_waits_for_borrowed_threads_to_come_back() {
        let budget = ThreadBudget::new(1);
        let borrowed = budget.take_up_to(1);
        assert_eq!(borrowed.count(), 1);
        std::thread::scope(|scope| {
            let waiter = scope.spawn(|| budget.hold().count());
            std::thread::sleep(std::time::Duration::from_millis(20));
            assert!(!waiter.is_finished(), "the budget is spent");
            drop(borrowed);
            assert_eq!(waiter.join().expect("the waiter panicked"), 1);
        });
    }
}
