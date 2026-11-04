//! CRAM index.

pub mod fs;
pub mod io;
pub mod record;

pub use self::record::Record;

/// A CRAM index.
pub type Index = Vec<Record>;
