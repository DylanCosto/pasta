//! Fitting threads and malloc arenas into an address-space limit (`ulimit -v`).
//!
//! glibc gives threads their own allocation arenas, up to eight per core, and reserves 64 MB of
//! address space for each; every thread also reserves its stack. vtools runs several pools of
//! `--threads` threads (decompression, parsing, formatting, compression), so under a limit of a
//! few GB the arenas and stacks alone would take the whole address space: then threads cannot
//! start and allocations fail (Rust aborts on a failed allocation). Under a limit,
//! [`fit_arenas`] caps the arenas at a quarter of it, and [`max_threads`] caps `--threads` so
//! that each thread has 64 MB of it. Without a limit nothing changes.

use std::sync::OnceLock;

/// The address space glibc reserves for each arena on 64-bit systems (`HEAP_MAX_SIZE`).
#[cfg(all(target_os = "linux", target_env = "gnu"))]
const ARENA_RESERVE: u64 = 64 << 20;

/// The address space a `--threads` thread is given under a limit: its stack and its share of
/// the batches in flight, in each of the pools it is counted in.
const THREAD_RESERVE: u64 = 64 << 20;

/// The process' address-space limit (`RLIMIT_AS`), `None` when there is none.
pub fn address_space() -> Option<u64> {
    static LIMIT: OnceLock<Option<u64>> = OnceLock::new();
    *LIMIT.get_or_init(|| {
        let mut limit = libc::rlimit {
            rlim_cur: 0,
            rlim_max: 0,
        };
        // SAFETY: getrlimit writes the limit into the local struct it is given.
        let read = unsafe { libc::getrlimit(libc::RLIMIT_AS, &mut limit) } == 0;
        (read && limit.rlim_cur != libc::RLIM_INFINITY).then_some(limit.rlim_cur)
    })
}

/// The most threads one `--threads` pool may start under the address-space limit (`None` without
/// one): one per 64 MB of it, at least one.
pub fn max_threads() -> Option<u64> {
    address_space().map(|limit| (limit / THREAD_RESERVE).max(1))
}

/// Caps glibc's malloc arenas so that they take at most a quarter of the address-space limit,
/// when the process has one and neither `MALLOC_ARENA_MAX` nor `glibc.malloc.arena_max` in
/// `GLIBC_TUNABLES` sets the cap (they are left to rule).
/// Call it before any thread starts.
pub fn fit_arenas() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    {
        let tunables = std::env::var_os("GLIBC_TUNABLES").unwrap_or_default();
        if std::env::var_os("MALLOC_ARENA_MAX").is_some()
            || tunables
                .to_string_lossy()
                .contains("glibc.malloc.arena_max")
        {
            return;
        }
        let Some(limit) = address_space() else {
            return;
        };
        // glibc's own cap is eight arenas per online CPU of the machine (not only those the
        // process may run on, which `taskset` or a cgroup limits), so a limit that leaves room
        // for those keeps it.
        // SAFETY: sysconf only reads a system setting.
        let online = unsafe { libc::sysconf(libc::_SC_NPROCESSORS_ONLN) };
        let cores = u64::try_from(online).unwrap_or(1).max(1);
        let arenas = (limit / 4 / ARENA_RESERVE).clamp(1, 8 * cores);
        let arenas = libc::c_int::try_from(arenas).unwrap_or(libc::c_int::MAX);
        // SAFETY: M_ARENA_MAX only changes how many arenas malloc creates from now on; no other
        // thread runs this early.
        unsafe { libc::mallopt(libc::M_ARENA_MAX, arenas) };
    }
}
