//! Alignment format utilities.

pub mod fs;
mod index;
pub mod io;
pub mod iter;
mod record;

pub use self::{index::Index, record::Record};
