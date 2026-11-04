//! The `##vtools_<command>Version` and `##vtools_<command>Command` header lines, written where
//! bcftools writes `##bcftools_<command>Version/Command` (bcftools' `bcf_hdr_append_version`).

use std::env;

use anyhow::Context;
use pasta::variant::Header;

/// Appends this program's version and command line to the header for command `cmd` (such as
/// `view`).
pub fn append_version(header: &mut Header, cmd: &str) -> anyhow::Result<()> {
    let version = format!(
        "##vtools_{cmd}Version={}+pasta-{}",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION")
    );
    header
        .append(version.as_bytes())
        .context("[bcf_hdr_append_version] Failed to add program information to header")?;

    let mut line = format!("##vtools_{cmd}Command=");
    for (i, arg) in env::args().skip(1).enumerate() {
        if i > 0 {
            line.push(' ');
        }
        if i > 0 && arg.contains(' ') {
            line.push('\'');
            line.push_str(&arg);
            line.push('\'');
        } else {
            line.push_str(&arg);
        }
    }
    line.push_str("; Date=");
    line.push_str(&ctime());
    header
        .append(line.as_bytes())
        .context("[bcf_hdr_append_version] Failed to add program information to header")?;
    Ok(())
}

/// Returns the local time as C's `ctime` formats it, without the newline.
fn ctime() -> String {
    // SAFETY: `time` with a null pointer only returns the current time.
    let now = unsafe { libc::time(std::ptr::null_mut()) };
    let mut buf = [0 as libc::c_char; 32];
    // SAFETY: `ctime_r` writes at most 26 bytes (including the NUL) into the 32-byte buffer.
    let p = unsafe { libc::ctime_r(&now, buf.as_mut_ptr()) };
    if p.is_null() {
        return String::new();
    }
    // SAFETY: `ctime_r` returned a pointer to the NUL-terminated string in `buf`.
    let s = unsafe { std::ffi::CStr::from_ptr(buf.as_ptr()) };
    s.to_string_lossy().trim_end().to_owned()
}
