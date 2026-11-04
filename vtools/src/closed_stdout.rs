//! A closed standard output (`vtools view x.bcf >&-`) fails the commands that write to it, as it
//! fails bcftools, instead of swallowing their output.
//!
//! Before `main`, Rust's start-up code opens `/dev/null` in place of a closed standard
//! descriptor, so writes to it would succeed and the command would exit with status 0 having
//! written nothing. [`occupy`] runs first, as a constructor, and puts a descriptor of the root
//! directory, open only for reading, in its place: every write to descriptor 1 then fails with
//! `EBADF`, as in bcftools, and the command reports its write error. Opening `/dev/stdout` for
//! output fails too (it is a directory), and no file vtools opens later (an input, an index, a
//! temporary file) can be given descriptor 1, so output never goes into it. bcftools does not
//! guard against that: its first input gets descriptor 1, and `-o /dev/stdout` truncates that
//! input and writes the output over it (README).
//!
//! Writes through `io::stdout` still succeed, as Rust's standard output takes `EBADF` for
//! success; the commands write their output through `output::stdout`, which does not.

/// Runs [`occupy`] before `main`, ahead of Rust's start-up code, which would see the closed
/// descriptor and open `/dev/null` there: from `.init_array` on Linux, and from
/// `__mod_init_func`, which dyld runs before `main`, on macOS.
#[cfg(any(target_os = "linux", target_os = "macos"))]
#[used]
#[cfg_attr(target_os = "linux", unsafe(link_section = ".init_array"))]
#[cfg_attr(target_os = "macos", unsafe(link_section = "__DATA,__mod_init_func"))]
static OCCUPY: extern "C" fn() = occupy;

/// Puts the root directory, open for reading, at descriptor 1 when it is closed.
#[cfg(any(target_os = "linux", target_os = "macos"))]
extern "C" fn occupy() {
    // SAFETY: plain system calls on descriptors this function opens or checks; no Rust state is
    // touched, so running before `main` is sound.
    unsafe {
        if libc::fcntl(libc::STDOUT_FILENO, libc::F_GETFD) != -1 {
            return;
        }
        let fd = libc::open(c"/".as_ptr(), libc::O_RDONLY | libc::O_DIRECTORY);
        // When standard input is closed too, the directory gets descriptor 0: move it to 1,
        // leaving 0 closed for Rust's start-up code to fill with /dev/null.
        if fd >= 0 && fd != libc::STDOUT_FILENO {
            libc::dup2(fd, libc::STDOUT_FILENO);
            libc::close(fd);
        }
    }
}
