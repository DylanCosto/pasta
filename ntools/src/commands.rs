//! Subcommands.

pub mod addreplacerg;
pub mod ampliconclip;
pub mod bedcov;
pub mod bgzip;
pub mod calmd;
pub mod cat;
pub mod checksum;
pub mod collate;
pub mod consensus;
pub mod coverage;
pub mod cram_size;
pub mod depth;
pub mod dict;
pub mod faidx;
pub mod fastq;
pub mod fixmate;
pub mod flags;
pub mod flagstat;
pub mod head;
pub mod idxstats;
pub mod import;
pub mod index;
pub mod markdup;
pub mod merge;
pub mod mpileup;
pub mod quickcheck;
pub mod reference;
pub mod reheader;
pub mod reset;
pub mod rmdup;
pub mod samples;
pub mod sort;
pub mod split;
pub mod stats;
pub mod view;

use std::fmt::Display;

use anyhow::bail;
use pasta::sam;

/// Fails unless `header`'s reference sequences are `first`'s, or the first of them: the same
/// names and lengths, in the same order.
///
/// Commands that read several inputs in step label every input's records with the first
/// input's reference sequence names, by index, so a differing dictionary would put reads on
/// the wrong sequence. (samtools does so; ntools rejects such inputs.) An input with fewer
/// sequences than the first, all at the same indexes, is labelled correctly.
pub(crate) fn check_reference_sequences(
    command: &str,
    first_input: impl Display,
    first: &sam::Header,
    input: impl Display,
    header: &sam::Header,
) -> anyhow::Result<()> {
    let (a, b) = (first.reference_sequences(), header.reference_sequences());

    let same = b.len() <= a.len()
        && a.iter()
            .zip(b)
            .all(|((a_name, a), (b_name, b))| a_name == b_name && a.length() == b.length());

    if !same {
        bail!(
            "[{command}] the reference sequences of {input} differ from those of {first_input} \
             (they must be the first input's, or the first of them, in the same order)"
        );
    }

    Ok(())
}
