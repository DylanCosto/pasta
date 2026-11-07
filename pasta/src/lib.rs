//! **pasta** reads, writes, and processes high-throughput sequencing data: SAM, BAM, and CRAM
//! alignments, VCF and BCF variant calls, BGZF compression, CSI/BAI/CRAI/tabix indexes, and FASTA
//! references, plus an HTSlib-compatible pileup engine.
//!
//! It began as a copy of the [noodles](https://github.com/zaeleus/noodles) crates that `ntools`
//! used, with that project's changes, and is developed independently from there.

pub mod bam;
pub mod bgzf;
pub mod core;
pub mod cram;
pub mod csi;
mod deflate;
pub mod fasta;
pub mod pileup;
pub mod sam;
pub mod tabix;
pub mod util;
pub mod variant;
