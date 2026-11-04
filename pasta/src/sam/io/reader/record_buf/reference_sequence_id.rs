//! SAM record buf reference sequence ID reader.

use std::{error, fmt};

use bstr::BString;

/// An error when a raw SAM record reference sequence ID fails to parse.
///
/// Records are normalized as HTSlib's `sam_parse1` normalizes them, so a name missing from the
/// header is read as `*` rather than failing with this error.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The reference sequence is not in the reference sequence dictionary.
    MissingReferenceSequenceDictionaryEntry(BString),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceDictionaryEntry(name) => {
                write!(f, "missing reference sequence dictionary entry: {name}")
            }
        }
    }
}
