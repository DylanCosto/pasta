//! VCF text: parsing lines into records and formatting records as lines, byte-compatible with
//! HTSlib.

pub(crate) mod format;
mod format_fields;
mod parse;

pub use self::{
    format::{FormatError, format_array, format_genotype, format_record},
    parse::{ParseError, parse, parse_shared},
};
