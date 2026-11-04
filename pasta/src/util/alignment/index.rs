use crate::bam;
use crate::cram::crai;
use crate::csi;

/// An alignment index.
pub enum Index {
    /// A SAM index.
    Sam(csi::Index),
    /// A BAM index.
    Bam(bam::Index),
    /// A CRAM index.
    Cram(crai::Index),
}
