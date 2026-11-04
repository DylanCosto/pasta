//! A write to a pipe or FIFO whose reader went away (`| head`) ends ntools as it ends samtools:
//! killed by SIGPIPE, with no message (status 141 in the shell), so that `set -o pipefail`
//! sees that the output was cut short.
//!
//! Rust ignores SIGPIPE before `main`, so such a write returns `EPIPE`, which the commands pass
//! up as an error. ntools keeps it so and raises SIGPIPE in [`die`] only once the command has
//! returned, so its temporary files (`sort` and `collate` spills, `faidx` index files) are
//! removed first; samtools, killed in the write itself, leaves them behind.

/// Ends the program as a write to a closed pipe ends samtools: killed by SIGPIPE. The signal is
/// given its default action and unblocked first, so this ends the program however SIGPIPE was
/// set up by the parent.
pub fn die() -> ! {
    // SAFETY: SIG_DFL installs no handler, and `set` is a local signal set that sigemptyset
    // initialises before it is used.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
        let mut set: libc::sigset_t = std::mem::zeroed();
        libc::sigemptyset(&mut set);
        libc::sigaddset(&mut set, libc::SIGPIPE);
        libc::pthread_sigmask(libc::SIG_UNBLOCK, &set, std::ptr::null_mut());
        libc::raise(libc::SIGPIPE);
    }
    // Not reached: the signal's default action ends the process.
    std::process::exit(128 + libc::SIGPIPE)
}
