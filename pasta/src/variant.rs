//! VCF and BCF variant data, ported from HTSlib 1.24's `vcf.c` with the same semantics.
//!
//! This is the engine under `vtools`: a [`Header`] (HTSlib's `bcf_hdr_t`) and a [`Record`]
//! (`bcf1_t`) that holds BCF-encoded data with lazy unpacking, plus VCF parsing and formatting
//! that reproduce HTSlib's output byte for byte, including its number formatting, missing and
//! vector-end sentinels, and warnings.
//!
//! Names follow HTSlib where a porter would look for them: each item's documentation names the
//! HTSlib function it corresponds to.

pub mod bcf;
pub mod encoding;
pub mod header;
pub mod index;
pub mod io;
pub mod log;
pub mod number;
pub mod readlist;
pub mod record;
pub mod regidx;
pub mod region;
pub mod synced;
pub mod utils;
pub mod vcf;

pub use self::{
    header::{
        Header, HeaderError, HeaderRecord, HeaderRecordType, HeaderTextFlavor, Number, Quoting,
        VCF_DEFAULT_VERSION, VCF44, VCF45, ValueType,
    },
    record::{
        FormatField, Genotype, GetError, InfoField, Record, RecordError, TypedValues, Unpack,
        UpdateError, Variant, VariantMatch, VariantType,
    },
};
