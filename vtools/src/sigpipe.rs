//! SIGPIPE as bcftools has it: a write to a pipe or FIFO whose reader went away (`| head`)
//! kills the program, with no message and the status of death by SIGPIPE (141 in the shell).
//!
//! Rust ignores SIGPIPE before `main`, so that such a write returns `EPIPE` instead; [`restore`]
//! puts back the default action, which C programs start with. Then a closed pipe stops the
//! program in whichever thread writes to it (a compression thread, a writer thread, or a write
//! to standard error), as in bcftools, and the paths that pass a broken pipe up as an error are
//! reached only when SIGPIPE is blocked; [`die`] ends those the same way. vtools writes to no
//! pipe or socket of its own (no child processes, no network), so nothing needs to survive a
//! closed pipe; `sort` removes its temporary files before it dies (see its `tmp_on_sigpipe`).
//!
//! bcftools started with SIGPIPE ignored or blocked (by a parent that left it so) gets `EPIPE`,
//! reports a write error of its own and exits with 255; vtools still dies of SIGPIPE then
//! (README), as its writers pass a closed pipe up instead of bcftools' messages.

/// Restores the default action of SIGPIPE (death), which Rust's start-up code replaced with
/// `SIG_IGN`.
pub fn restore() {
    // SAFETY: setting a standard signal's disposition to SIG_DFL installs no handler.
    unsafe { libc::signal(libc::SIGPIPE, libc::SIG_DFL) };
}

/// Ends the program as a write to a closed pipe ends bcftools when SIGPIPE has its default
/// action: killed by SIGPIPE. It is unblocked and given its default action first, so this ends
/// the program also when SIGPIPE was blocked, where bcftools reports a write error instead.
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
