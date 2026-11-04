//! Reading and writing VCF and BCF files.

pub(crate) mod bgzf_out;
pub mod reader;
pub mod writer;

pub use self::{
    reader::{
        Compression, InputFormat, Query, ReadError, Reader, has_bgzf_eof, stdin_has_bgzf_eof,
    },
    writer::{IndexRequest, OutputType, RecordSpan, WriteError, WriteOptions, Writer},
};
