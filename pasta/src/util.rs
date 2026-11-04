//! Utilities: unified interfaces for reading and writing alignment (BAM/CRAM/SAM) data, and
//! starting worker threads within the system's limits. Variant data has its own module,
//! [`crate::variant`].

pub mod alignment;
pub mod inputs;
pub mod io;
pub mod threads;
