//! Consequence bits, the per-record consequence (bcftools' `vcsq_t`, `csq_t`), and their
//! text form (`kput_vcsq`).

use crate::gff::{Gff, gf, type2gff_string};

pub(super) const CSQ_PRINTED_UPSTREAM: u32 = 1 << 0;
pub(super) const CSQ_SYNONYMOUS_VARIANT: u32 = 1 << 1;
pub(super) const CSQ_MISSENSE_VARIANT: u32 = 1 << 2;
pub(super) const CSQ_STOP_LOST: u32 = 1 << 3;
pub(super) const CSQ_STOP_GAINED: u32 = 1 << 4;
pub(super) const CSQ_INFRAME_DELETION: u32 = 1 << 5;
pub(super) const CSQ_INFRAME_INSERTION: u32 = 1 << 6;
pub(super) const CSQ_FRAMESHIFT_VARIANT: u32 = 1 << 7;
pub(super) const CSQ_SPLICE_ACCEPTOR: u32 = 1 << 8;
pub(super) const CSQ_SPLICE_DONOR: u32 = 1 << 9;
pub(super) const CSQ_START_LOST: u32 = 1 << 10;
pub(super) const CSQ_SPLICE_REGION: u32 = 1 << 11;
pub(super) const CSQ_STOP_RETAINED: u32 = 1 << 12;
pub(super) const CSQ_UTR5: u32 = 1 << 13;
pub(super) const CSQ_UTR3: u32 = 1 << 14;
pub(super) const CSQ_NON_CODING: u32 = 1 << 15;
pub(super) const CSQ_INTRON: u32 = 1 << 16;
pub(super) const CSQ_INFRAME_ALTERING: u32 = 1 << 18;
/// Adds `*` in front of the consequence string.
pub(super) const CSQ_UPSTREAM_STOP: u32 = 1 << 19;
/// Removes start/stop consequences of incomplete CDSs.
pub(super) const CSQ_INCOMPLETE_CDS: u32 = 1 << 20;
pub(super) const CSQ_CODING_SEQUENCE: u32 = 1 << 21;
pub(super) const CSQ_ELONGATION: u32 = 1 << 22;
pub(super) const CSQ_TRUNCATION: u32 = 1 << 23;
pub(super) const CSQ_START_RETAINED: u32 = 1 << 24;

/// Haplotype-aware consequences, printed in one record only (the others get `@POS`).
pub(super) const CSQ_COMPOUND: u32 = CSQ_SYNONYMOUS_VARIANT
    | CSQ_MISSENSE_VARIANT
    | CSQ_STOP_LOST
    | CSQ_STOP_GAINED
    | CSQ_INFRAME_DELETION
    | CSQ_INFRAME_INSERTION
    | CSQ_FRAMESHIFT_VARIANT
    | CSQ_START_LOST
    | CSQ_STOP_RETAINED
    | CSQ_INFRAME_ALTERING
    | CSQ_INCOMPLETE_CDS
    | CSQ_UPSTREAM_STOP
    | CSQ_START_RETAINED
    | CSQ_ELONGATION
    | CSQ_TRUNCATION;
pub(super) const CSQ_START_STOP: u32 =
    CSQ_STOP_LOST | CSQ_STOP_GAINED | CSQ_STOP_RETAINED | CSQ_START_LOST | CSQ_START_RETAINED;

/// Whether the strand is printed (bcftools' `CSQ_PRN_STRAND`).
fn prn_strand(csq: u32) -> bool {
    csq & CSQ_COMPOUND != 0
        && csq
            & (CSQ_SPLICE_ACCEPTOR
                | CSQ_SPLICE_DONOR
                | CSQ_SPLICE_REGION
                | CSQ_ELONGATION
                | CSQ_TRUNCATION)
            == 0
}
/// Consequences that print the transcript (bcftools' `CSQ_PRN_TSCRIPT`).
pub(super) const CSQ_PRN_TSCRIPT: u32 = !(CSQ_INTRON | CSQ_NON_CODING);
/// Consequences that get `NMD_transcript` (bcftools' `CSQ_PRN_NMD`).
const CSQ_PRN_NMD: u32 = !(CSQ_INTRON | CSQ_NON_CODING);

/// The names of the consequence bits (bcftools' `csq_strings`).
const CSQ_STRINGS: [Option<&str>; 25] = [
    None,
    Some("synonymous"),
    Some("missense"),
    Some("stop_lost"),
    Some("stop_gained"),
    Some("inframe_deletion"),
    Some("inframe_insertion"),
    Some("frameshift"),
    Some("splice_acceptor"),
    Some("splice_donor"),
    Some("start_lost"),
    Some("splice_region"),
    Some("stop_retained"),
    Some("5_prime_utr"),
    Some("3_prime_utr"),
    Some("non_coding"),
    Some("intron"),
    None,
    Some("inframe_altering"),
    None,
    None,
    Some("coding_sequence"),
    Some("feature_elongation"),
    Some("feature_truncation"),
    Some("start_retained"),
];

/// A buffered record: its serial number (bcftools compares `bcf1_t` pointers) and position.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct RecRef {
    pub serial: u64,
    pub pos: i64,
}

/// A buffered record's slot: the buffer of its position and its index there.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct VrecHandle {
    pub vbuf: u64,
    pub slot: u32,
}

/// The variant string of a consequence, as a `kstring_t`: bcftools tests the pointer (`s`)
/// when comparing consequences but the length (`l`) when printing, and the two can disagree.
#[derive(Clone, Debug, Default)]
pub(super) struct VStr {
    pub s: Option<Vec<u8>>,
    /// Whether `l` is nonzero.
    pub shown: bool,
}

impl VStr {
    /// A string that is set and printed.
    pub(super) fn new(s: Vec<u8>) -> Self {
        Self {
            s: Some(s),
            shown: true,
        }
    }
}

/// A consequence of one VCF record (bcftools' `vcsq_t`).
#[derive(Clone, Debug, Default)]
pub(super) struct VCsq {
    /// bcftools' one-bit `strand` field: forward is 1, reverse and unknown 0.
    pub strand: u32,
    /// `CSQ_*` bits.
    pub ty: u32,
    pub trid: u32,
    pub vcf_ial: u32,
    /// The biotype ([`gf`]).
    pub biotype: u32,
    /// The gene; `None` for a gene without a name, as bcftools compares name pointers.
    pub gene: Option<u32>,
    /// With [`CSQ_PRINTED_UPSTREAM`], the record that carries the consequence (`@POS`).
    pub ref_: Option<RecRef>,
    pub vstr: VStr,
}

/// A consequence tied to a haplotype (bcftools' `csq_t`).
#[derive(Clone, Debug, Default)]
pub(super) struct CsqT {
    pub vrec: Option<VrecHandle>,
    /// The index of the consequence in its record, for FORMAT/BCSQ.
    pub idx: i32,
    pub ty: VCsq,
}

impl CsqT {
    /// A consequence of a transcript at a position (the fields bcftools sets after `memset`).
    pub(super) fn new(ty: u32, tr: &crate::gff::Transcript, gff: &Gff, vcf_ial: u32) -> Self {
        Self {
            vrec: None,
            idx: 0,
            ty: VCsq {
                strand: tr.strand.code() & 1,
                ty,
                trid: tr.id,
                vcf_ial,
                biotype: tr.ty,
                gene: gene_ptr(gff, tr.gene),
                ref_: None,
                vstr: VStr::default(),
            },
        }
    }
}

/// The gene of a consequence: bcftools stores the gene's name pointer, NULL for genes without
/// a gene line.
pub(super) fn gene_ptr(gff: &Gff, gene: u32) -> Option<u32> {
    gff.genes[gene as usize].name.as_ref().map(|_| gene)
}

/// Formats a consequence like "inframe_deletion|XYZ|ENST01|+|5TY>5I|121ACG>A+124TA>T"
/// (bcftools' `kput_vcsq`).
pub(super) fn kput_vcsq(gff: &Gff, csq: &VCsq, out: &mut Vec<u8>) {
    let mut csq_type = csq.ty;

    // Remove start/stop from incomplete CDS, but only if there is another consequence.
    if csq_type & CSQ_INCOMPLETE_CDS != 0
        && csq_type & !(CSQ_START_STOP | CSQ_INCOMPLETE_CDS | CSQ_UPSTREAM_STOP) != 0
    {
        csq_type &= !(CSQ_START_STOP | CSQ_INCOMPLETE_CDS);
    }
    // Remove missense from start/stops.
    if csq_type & CSQ_START_STOP != 0 && csq_type & CSQ_MISSENSE_VARIANT != 0 {
        csq_type &= !CSQ_MISSENSE_VARIANT;
    }
    if csq_type & CSQ_PRINTED_UPSTREAM != 0
        && let Some(r) = csq.ref_
    {
        out.push(b'@');
        out.extend_from_slice(((r.pos + 1) as i32).to_string().as_bytes());
        return;
    }
    if csq_type & CSQ_UPSTREAM_STOP != 0 {
        out.push(b'*');
    }
    let mut has_csq = false;
    for (i, s) in CSQ_STRINGS.iter().enumerate().skip(1) {
        let Some(s) = s else { continue };
        if csq_type & (1 << i) != 0 {
            if has_csq {
                out.push(b'&');
            }
            has_csq = true;
            out.extend_from_slice(s.as_bytes());
        }
    }
    if csq.biotype == gf::NMD && csq_type & CSQ_PRN_NMD != 0 {
        if has_csq {
            out.push(b'&');
        }
        out.extend_from_slice(b"NMD_transcript");
    }
    out.push(b'|');
    if let Some(g) = csq.gene
        && let Some(name) = &gff.genes[g as usize].name
    {
        out.extend_from_slice(name);
    }
    out.push(b'|');
    if csq_type & CSQ_PRN_TSCRIPT != 0 {
        out.extend_from_slice(gff.transcript_name(csq.trid));
    }
    out.push(b'|');
    out.extend_from_slice(type2gff_string(csq.biotype).as_bytes());
    if prn_strand(csq_type) || csq.vstr.shown {
        out.extend_from_slice(match csq.strand {
            1 => b"|+",
            0 => b"|-",
            _ => b"|.",
        });
    }
    if csq.vstr.shown
        && let Some(s) = &csq.vstr.s
    {
        out.extend_from_slice(s);
    }
}
