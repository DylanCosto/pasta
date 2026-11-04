use std::io;

use crate::bgzf::VirtualPosition;

/// A BGZF reader.
pub trait Read: io::Read {
    /// Returns the current virtual position.
    fn virtual_position(&self) -> VirtualPosition;
}
