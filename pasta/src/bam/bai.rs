//! BAM index (BAI) and fields.
//!
//! A BAM index (BAI) is used with an associated coordinate-sorted BAM file that allows random
//! access to records, e.g., [querying].
//!
//! The index contains a list of reference sequences parallel to the one defined in the BAM file.
//! Each indexed reference sequence has a calculated set of hierarchical bins at different
//! granularities. The bins then define a list of physical file positions in the BAM to search for
//! overlapping records.
//!
//! When reading entire BAM files sequentially, a BAM index is not necessary.
//!
//! [querying]: crate::bam::io::Reader::query
//!
//! # Examples
//!
//! ## Reading a BAM index
//!
//! ```no_run
//! # use std::io;
//! use pasta::bam::bai;
//! let index = bai::fs::read("sample.bam.bai")?;
//! # Ok::<(), io::Error>(())
//! ```

pub mod fs;
pub mod io;

use crate::csi::binning_index::{self, index::reference_sequence::index::LinearIndex};

const DEPTH: u8 = 5;

const MAGIC_NUMBER: [u8; 4] = *b"BAI\x01";

/// A BAM index.
pub type Index = binning_index::Index<LinearIndex>;
