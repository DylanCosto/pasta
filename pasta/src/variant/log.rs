//! HTSlib-style diagnostics: `[W::function] message` lines on stderr.
//!
//! Messages are worded and attributed (the function name in brackets) as HTSlib words them, so
//! stderr matches bcftools. The verbosity follows HTSlib's `hts_verbose` (3, warnings, by
//! default).
//!
//! Work done on worker threads can [`capture`] its messages and [`replay`] them later in input
//! order, so parallel parsing reports exactly what serial parsing would. Messages HTSlib prints
//! only once per process (`static int warned`) go through [`warn_once`] with a flag that is
//! checked when the message is finally emitted, not when it is captured.

use std::{
    cell::RefCell,
    fmt,
    io::Write,
    sync::atomic::{AtomicBool, AtomicI32, Ordering},
};

/// HTSlib's default `hts_verbose` (`HTS_LOG_WARNING`).
const DEFAULT_VERBOSITY: i32 = 3;

static VERBOSITY: AtomicI32 = AtomicI32::new(DEFAULT_VERBOSITY);

/// A message severity, with HTSlib's numeric log levels.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Level {
    /// `[E::...]`, shown at verbosity 1 and above.
    Error,
    /// `[W::...]`, shown at verbosity 3 and above.
    Warning,
    /// `[I::...]`, shown at verbosity 4 and above.
    Info,
    /// `[D::...]`, shown at verbosity 5 and above.
    Debug,
}

impl Level {
    fn threshold(self) -> i32 {
        match self {
            Self::Error => 1,
            Self::Warning => 3,
            Self::Info => 4,
            Self::Debug => 5,
        }
    }

    fn tag(self) -> char {
        match self {
            Self::Error => 'E',
            Self::Warning => 'W',
            Self::Info => 'I',
            Self::Debug => 'D',
        }
    }
}

/// Sets HTSlib's verbosity (`hts_verbose`).
pub fn set_verbosity(level: i32) {
    VERBOSITY.store(level, Ordering::Relaxed);
}

/// Returns HTSlib's verbosity (`hts_verbose`).
#[must_use]
pub fn verbosity() -> i32 {
    VERBOSITY.load(Ordering::Relaxed)
}

/// Returns whether messages of `level` are shown.
#[must_use]
pub fn enabled(level: Level) -> bool {
    level.threshold() <= verbosity()
}

/// A diagnostic captured on a worker thread, to be replayed in order.
#[derive(Debug)]
pub struct Message {
    once: Option<&'static AtomicBool>,
    text: String,
}

thread_local! {
    static CAPTURED: RefCell<Option<Vec<Message>>> = const { RefCell::new(None) };
}

/// Runs `f`, capturing the diagnostics it emits on this thread instead of printing them.
pub fn capture<R>(f: impl FnOnce() -> R) -> (R, Vec<Message>) {
    let previous = CAPTURED.with(|c| c.borrow_mut().replace(Vec::new()));
    let result = f();
    let messages = CAPTURED.with(|c| {
        let mut slot = c.borrow_mut();
        let messages = slot.take().unwrap_or_default();
        *slot = previous;
        messages
    });
    (result, messages)
}

/// Prints captured diagnostics, applying once-only flags now.
pub fn replay(messages: Vec<Message>) {
    for message in messages {
        if let Some(flag) = message.once
            && flag.swap(true, Ordering::Relaxed)
        {
            continue;
        }
        emit_text(message.text);
    }
}

fn emit_text(text: String) {
    let captured = CAPTURED.with(|c| {
        if let Some(buf) = c.borrow_mut().as_mut() {
            buf.push(Message { once: None, text });
            None
        } else {
            Some(text)
        }
    });
    if let Some(text) = captured {
        let mut stderr = std::io::stderr().lock();
        let _ = stderr.write_all(text.as_bytes());
    }
}

fn format(level: Level, func: &str, args: fmt::Arguments<'_>) -> String {
    format!("[{}::{}] {}\n", level.tag(), func, args)
}

/// Emits a message at `level` attributed to the HTSlib function `func`.
pub fn log(level: Level, func: &str, args: fmt::Arguments<'_>) {
    if enabled(level) {
        emit_text(format(level, func, args));
    }
}

/// Emits a line as it stands, whatever the verbosity (a program's own `fprintf(stderr, ...)`),
/// so that work on worker threads can [`capture`] it like HTSlib's messages.
pub fn line(args: fmt::Arguments<'_>) {
    emit_text(format!("{args}\n"));
}

/// Emits a warning (`[W::func] ...`).
pub fn warn(func: &str, args: fmt::Arguments<'_>) {
    log(Level::Warning, func, args);
}

/// Emits an error message (`[E::func] ...`).
pub fn error(func: &str, args: fmt::Arguments<'_>) {
    log(Level::Error, func, args);
}

/// Emits a warning only the first time `flag` is seen unset, like HTSlib's `static int warned`.
pub fn warn_once(flag: &'static AtomicBool, func: &str, args: fmt::Arguments<'_>) {
    log_once(flag, Level::Warning, func, args);
}

/// Emits a message only the first time `flag` is seen unset.
pub fn log_once(flag: &'static AtomicBool, level: Level, func: &str, args: fmt::Arguments<'_>) {
    if !enabled(level) {
        return;
    }
    let text = format(level, func, args);
    let pending = CAPTURED.with(|c| {
        if let Some(buf) = c.borrow_mut().as_mut() {
            buf.push(Message {
                once: Some(flag),
                text,
            });
            None
        } else {
            Some(text)
        }
    });
    if let Some(text) = pending
        && !flag.swap(true, Ordering::Relaxed)
    {
        let mut stderr = std::io::stderr().lock();
        let _ = stderr.write_all(text.as_bytes());
    }
}

/// Returns whether a once-only flag has fired. Workers use this to skip formatting messages that
/// can no longer print.
#[must_use]
pub fn fired(flag: &AtomicBool) -> bool {
    flag.load(Ordering::Relaxed)
}

/// Quotes and escapes `s` as HTSlib's `hts_strprint` does into a buffer of `buflen` bytes,
/// truncating with `...` when it does not fit.
#[must_use]
pub fn strprint(s: &[u8], quote: u8, buflen: usize) -> String {
    let qlen = usize::from(quote != 0);
    let mut out: Vec<u8> = Vec::with_capacity(s.len().min(buflen) + 4);
    if quote != 0 {
        out.push(quote);
    }
    for &c in s {
        let (esc, clen) = match c {
            b'\n' => (b'n', 2),
            b'\r' => (b'r', 2),
            b'\t' => (b't', 2),
            0 => (b'0', 2),
            b'\\' => (b'\\', 2),
            _ if c == quote && quote != 0 => (c, 2),
            _ if (0x20..0x7f).contains(&c) => (c, 1),
            _ => (c, 4),
        };
        if out.len() + clen + qlen >= buflen {
            while out.len() + 3 + qlen >= buflen {
                out.pop();
            }
            if quote != 0 {
                out.push(quote);
            }
            out.extend_from_slice(b"...");
            return String::from_utf8_lossy(&out).into_owned();
        }
        match clen {
            4 => out.extend_from_slice(format!("\\x{c:02X}").as_bytes()),
            2 => {
                out.push(b'\\');
                out.push(esc);
            }
            _ => out.push(c),
        }
    }
    if quote != 0 {
        out.push(quote);
    }
    String::from_utf8_lossy(&out).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strprint_escapes_and_truncates() {
        assert_eq!(strprint(b"a\tb", b'"', 320), "\"a\\tb\"");
        assert_eq!(strprint(b"\x01", b'"', 320), "\"\\x01\"");
        let long = vec![b'x'; 400];
        let s = strprint(&long, b'"', 320);
        assert!(s.ends_with("\"..."));
        assert_eq!(s.len(), 319);
    }

    #[test]
    fn capture_and_replay_respect_once() {
        static FLAG: AtomicBool = AtomicBool::new(false);
        let ((), messages) = capture(|| {
            warn_once(&FLAG, "f", format_args!("x"));
            warn_once(&FLAG, "f", format_args!("x"));
        });
        assert_eq!(messages.len(), 2);
        assert!(!fired(&FLAG));
    }
}
