//! FASTA index (FAI) and fields.

pub mod fs;
mod index;
pub mod io;
mod record;

pub use self::{index::Index, record::Record};
