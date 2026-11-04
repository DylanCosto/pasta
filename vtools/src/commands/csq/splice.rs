//! Consequences around exon boundaries: splice sites, starts, stops, and UTRs (bcftools'
//! `splice_t` and the `splice_*` functions).
//!
//! The arithmetic mixes bcftools' unsigned 32-bit coordinates (`ref_beg`, `ref_end`, exon
//! bounds) with signed ones (`pos`, `tbeg`, `tend`); the comparisons follow C's conversions.

use super::{
    Csq, CurRec, N_REF_PAD,
    translate::{KStr, strncmp_eq},
    types::{
        CSQ_CODING_SEQUENCE, CSQ_ELONGATION, CSQ_FRAMESHIFT_VARIANT, CSQ_INFRAME_DELETION,
        CSQ_SPLICE_ACCEPTOR, CSQ_SPLICE_DONOR, CSQ_SPLICE_REGION, CSQ_START_LOST,
        CSQ_START_RETAINED, CSQ_STOP_LOST, CSQ_SYNONYMOUS_VARIANT, CSQ_TRUNCATION, CSQ_UTR3,
        CSQ_UTR5, CsqT,
    },
};
use crate::gff::{Gff, N_SPLICE_DONOR, N_SPLICE_REGION_INTRON, Strand, UtrKind};

/// REF equals ALT: not applicable, skip.
pub(super) const SPLICE_VAR_REF: i32 = 0;
/// Splice acceptor or similar: set and done, does not overlap the region.
pub(super) const SPLICE_OUTSIDE: i32 = 1;
/// Overlaps the coding region: set, but a coding prediction is needed.
pub(super) const SPLICE_INSIDE: i32 = 2;
/// An indel overlaps the region boundary: set, but the consequence is unknown.
pub(super) const SPLICE_OVERLAP: i32 = 3;

const DONOR: u32 = N_SPLICE_DONOR;
const INTRON: u32 = N_SPLICE_REGION_INTRON;

/// A variant allele near a feature boundary (bcftools' `splice_t`).
pub(super) struct Splice<'a> {
    pub tr: u32,
    tr_beg: u32,
    tr_end: u32,
    tr_strand: Strand,
    tr_id: u32,
    /// The transcript's padded reference, when it is loaded.
    tr_ref: Option<&'a [u8]>,
    pub vcf_pos: i32,
    pub vcf_rlen: i32,
    pub vcf_alen: i32,
    pub vcf_ial: i32,
    /// The REF allele and the alleles after it (bcftools reads `rec->d.allele[0]` up to
    /// `rlen`, past its end when END makes `rlen` longer).
    vcf_ref: &'a [u8],
    vcf_alt: &'a [u8],
    /// Check the distance from the exon start (forward) or end (reverse).
    pub check_acceptor: bool,
    /// The first coding exon: check the first (forward) or last (reverse) codon.
    pub check_start: bool,
    /// The last coding exon: check the last (forward) or first (reverse) codon.
    pub check_stop: bool,
    pub check_donor: bool,
    /// Check for splices at this end (not in the first or last exon).
    pub check_region_beg: bool,
    pub check_region_end: bool,
    /// Check splice sites only outside UTRs.
    pub check_utr: bool,
    /// Set `kref`, `kalt` (and check for synonymous events).
    pub set_refalt: bool,
    pub csq: u32,
    /// The numbers of bases trimmed from the beginning and end of REF and ALT.
    pub tbeg: i32,
    pub tend: i32,
    /// The REF range that matters for consequences.
    pub ref_beg: u32,
    pub ref_end: u32,
    pub kref: KStr,
    pub kalt: KStr,
}

/// The byte at `i` of a C string pointer's buffer, 0 outside it.
#[inline]
fn at(s: &[u8], i: i64) -> u8 {
    usize::try_from(i)
        .ok()
        .and_then(|i| s.get(i))
        .copied()
        .unwrap_or(0)
}

/// C's `strlen`.
fn strlen(s: &[u8]) -> usize {
    memchr::memchr(0, s).unwrap_or(s.len())
}

/// `n` bytes of `s` from `from` (0 outside it), for `kputsn`.
fn bytes(s: &[u8], from: i64, n: i64) -> impl Iterator<Item = u8> + '_ {
    (0..n.max(0)).map(move |i| at(s, from + i))
}

impl<'a> Splice<'a> {
    /// A splice of the current record's REF in transcript `tr` (bcftools' `splice_init`).
    pub(super) fn new(cur: &'a CurRec, gff: &Gff, tr: u32, tr_ref: Option<&'a [u8]>) -> Self {
        let mut s = Self {
            tr,
            tr_beg: 0,
            tr_end: 0,
            tr_strand: Strand::Unknown,
            tr_id: 0,
            tr_ref,
            vcf_pos: cur.pos as i32,
            vcf_rlen: cur.rlen as i32,
            vcf_alen: 0,
            vcf_ial: 0,
            vcf_ref: cur.allele_raw(0),
            vcf_alt: &[],
            check_acceptor: false,
            check_start: false,
            check_stop: false,
            check_donor: false,
            check_region_beg: false,
            check_region_end: false,
            check_utr: false,
            set_refalt: false,
            csq: 0,
            tbeg: 0,
            tend: 0,
            ref_beg: 0,
            ref_end: 0,
            kref: KStr::default(),
            kalt: KStr::default(),
        };
        s.set_tr(gff, tr);
        s
    }

    /// Switches to another transcript (`splice.tr = ...`).
    pub(super) fn set_tr(&mut self, gff: &Gff, tr: u32) {
        let t = &gff.transcripts[tr as usize];
        self.tr = tr;
        self.tr_beg = t.beg;
        self.tr_end = t.end;
        self.tr_strand = t.strand;
        self.tr_id = t.id;
    }

    /// Selects the alternate allele (`splice.vcf.alt = rec->d.allele[ial]`).
    pub(super) fn set_alt(&mut self, cur: &'a CurRec, ial: usize) {
        self.vcf_alt = cur.allele_raw(ial);
        self.vcf_ial = ial as i32;
    }

    fn tref(&self) -> &[u8] {
        self.tr_ref.unwrap_or(&[])
    }

    /// Builds `kref` and `kalt` over the region starting at `beg` (`len > 0`, filled from
    /// the left) or ending at `beg` (`len < 0`) (bcftools' `splice_build_hap`).
    fn build_hap(&mut self, beg: u32, len: i32) {
        let (rlen, alen, rbeg, abeg): (i64, i64, i64, i64);
        let mut rlen_v: i64;
        let mut alen_v: i64;
        if len < 0 {
            rlen = -i64::from(len);
            alen = rlen;
            rbeg = i64::from(beg) - rlen + 1;
            let mut dlen = i64::from(self.vcf_alen) - i64::from(self.vcf_rlen);
            if dlen < 0 && beg < self.ref_end {
                // incomplete del, beg is in the middle
                dlen += i64::from(self.ref_end - beg);
            }
            abeg = rbeg + dlen;
        } else if beg < self.tr_beg {
            // Very short exons and introns; not real biology (test/csq/ENSCAFT00000047742).
            rbeg = i64::from(self.tr_beg);
            abeg = rbeg;
            rlen = 0;
            alen = 0;
        } else {
            rbeg = i64::from(beg);
            abeg = rbeg;
            rlen = i64::from(len);
            alen = rlen;
        }
        rlen_v = rlen;
        alen_v = alen;
        let pad = N_REF_PAD as i64;
        let tr_beg = i64::from(self.tr_beg);
        let tr_end = i64::from(self.tr_end);
        let pos = i64::from(self.vcf_pos);
        let vrlen = i64::from(self.vcf_rlen);
        let valen = i64::from(self.vcf_alen);
        let tref = self.tr_ref.unwrap_or(&[]);
        self.kref.clear();
        self.kalt.clear();

        // the part before vcf.ref, in vcf.ref, and after vcf.ref
        let roff = if rbeg < pos {
            self.kref.put(bytes(tref, pad + rbeg - tr_beg, pos - rbeg));
            0
        } else {
            rbeg - pos
        };
        if roff < vrlen && (self.kref.len() as i64) < rlen_v {
            let mut n = vrlen - roff;
            if n > rlen_v - self.kref.len() as i64 {
                n = rlen_v - self.kref.len() as i64;
            }
            self.kref.put(bytes(self.vcf_ref, roff, n));
        }
        let end = pos + vrlen; // just after the ref allele
        if (self.kref.len() as i64) < rlen_v {
            let kl = self.kref.len() as i64;
            if end + rlen_v - kl - 1 > tr_end {
                // trim, the requested sequence is too long
                rlen_v -= end + rlen_v - kl - 1 - tr_end;
            }
            if kl < rlen_v {
                self.kref.put(bytes(tref, pad + end - tr_beg, rlen_v - kl));
            }
        }

        let mut aoff = if abeg < pos {
            self.kalt.put(bytes(tref, pad + abeg - tr_beg, pos - abeg));
            0
        } else {
            abeg - pos
        };
        if aoff < valen && (self.kalt.len() as i64) < alen_v {
            let mut n = valen - aoff;
            if n > alen_v - self.kalt.len() as i64 {
                n = alen_v - self.kalt.len() as i64;
            }
            self.kalt.put(bytes(self.vcf_alt, aoff, n));
            aoff -= n;
        }
        if aoff < 0 {
            aoff = 0;
        } else {
            aoff -= 1;
        }
        if (self.kalt.len() as i64) < alen_v {
            let kl = self.kalt.len() as i64;
            if end + alen_v + aoff - kl - 1 > tr_end {
                alen_v -= end + alen_v + aoff - kl - 1 - tr_end;
            }
            if alen_v > 0 && alen_v > kl {
                self.kalt
                    .put(bytes(tref, aoff + pad + end - tr_beg, alen_v - kl));
            }
        }
    }

    /// `kref`/`kalt` as the pair compared by bcftools when `set_refalt` is on and both were
    /// written (`ref && alt`: a string never written is a NULL pointer in C, as for an ALT
    /// that `build_hap` fills with no bases).
    fn refalt(&self) -> Option<(&KStr, &KStr)> {
        (self.set_refalt && !self.kref.is_null() && !self.kalt.is_null())
            .then_some((&self.kref, &self.kalt))
    }
}

impl Csq {
    /// The UTRs overlapping `beg..=end` of the current sequence.
    fn utr_overlaps(&mut self, beg: u32, end: u32) -> Vec<u32> {
        self.gff
            .idx_utr
            .overlap(&self.chr_gff, beg, end)
            .map(|r| *r.payload)
            .collect()
    }

    /// Stages a UTR consequence for the first of `utrs` in transcript `trid`; returns its
    /// type, or 0 when none is (bcftools' `csq_stage_utr`).
    fn csq_stage_utr(
        &mut self,
        cur: &CurRec,
        utrs: &[u32],
        trid: u32,
        ty: u32,
        ial: i32,
    ) -> anyhow::Result<u32> {
        for &u in utrs {
            let utr = self.gff.utrs[u as usize];
            let tr = &self.gff.transcripts[utr.tr as usize];
            if tr.id != trid {
                continue;
            }
            let t = if utr.which == UtrKind::Prime5 {
                CSQ_UTR5
            } else {
                CSQ_UTR3
            } | ty;
            let mut csq = CsqT::new(t, tr, &self.gff, ial as u32);
            let id = tr.id;
            self.claim_tscript_allele(cur, id, ial as u32)?;
            self.csq_stage(cur, &mut csq)?;
            return Ok(csq.ty.ty);
        }
        Ok(0)
    }

    /// Stages a splice consequence (bcftools' `csq_stage_splice`).
    fn csq_stage_splice(&mut self, cur: &CurRec, splice: &Splice<'_>) -> anyhow::Result<()> {
        if splice.csq == 0 {
            return Ok(());
        }
        let tr = &self.gff.transcripts[splice.tr as usize];
        let mut csq = CsqT::new(splice.csq, tr, &self.gff, splice.vcf_ial as u32);
        let id = tr.id;
        self.claim_tscript_allele(cur, id, splice.vcf_ial as u32)?;
        self.csq_stage(cur, &mut csq)
    }

    /// Checks for an adjacent UTR of the splice's transcript over `beg..=end` and stages it
    /// (the `check_utr` blocks); returns the staged type.
    fn splice_utr(
        &mut self,
        cur: &CurRec,
        splice: &Splice<'_>,
        beg: u32,
        end: u32,
    ) -> anyhow::Result<u32> {
        let utrs = self.utr_overlaps(beg, end);
        if utrs.is_empty() {
            return Ok(0);
        }
        self.csq_stage_utr(cur, &utrs, splice.tr_id, splice.csq, splice.vcf_ial)
    }

    /// An insertion (bcftools' `splice_csq_ins`).
    fn splice_csq_ins(
        &mut self,
        cur: &CurRec,
        s: &mut Splice<'_>,
        ex_beg: u32,
        ex_end: u32,
    ) -> anyhow::Result<i32> {
        // coordinates that matter for consequences, eg AC>ACG trimmed to C>CG, 1bp before
        // and after the inserted bases
        if s.tbeg != 0 || at(s.vcf_ref, 0) != at(s.vcf_alt, 0) {
            s.ref_beg = (s.vcf_pos + s.tbeg - 1) as u32;
            s.ref_end = (s.vcf_pos + s.vcf_rlen - s.tend) as u32;
        } else {
            if s.tend != 0 {
                s.tend -= 1;
            }
            s.ref_beg = s.vcf_pos as u32;
            s.ref_end = (s.vcf_pos + s.vcf_rlen - s.tend) as u32;
        }
        if s.ref_beg >= ex_end {
            // fully outside, beyond the exon
            if s.check_utr {
                let b = s.ref_beg.wrapping_add(1);
                if self.splice_utr(cur, s, b, b)? != 0 {
                    return Ok(SPLICE_OUTSIDE); // overlaps utr
                }
            }
            if !s.check_region_end {
                return Ok(SPLICE_OUTSIDE);
            }
            if s.set_refalt {
                s.build_hap(ex_end.wrapping_add(1), INTRON as i32);
            }
            if s.ref_beg < ex_end.wrapping_add(INTRON) && s.ref_end > ex_end.wrapping_add(DONOR) {
                s.csq |= CSQ_SPLICE_REGION;
                if s.refalt()
                    .is_some_and(|(r, a)| strncmp_eq(r, 0, a, 0, INTRON as usize))
                {
                    s.csq |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
            if s.ref_beg < ex_end.wrapping_add(DONOR) {
                if s.check_donor && s.tr_strand == Strand::Fwd {
                    s.csq |= CSQ_SPLICE_DONOR;
                }
                if s.check_acceptor && s.tr_strand == Strand::Rev {
                    s.csq |= CSQ_SPLICE_ACCEPTOR;
                }
                if s.refalt()
                    .is_some_and(|(r, a)| strncmp_eq(r, 0, a, 0, DONOR as usize))
                {
                    s.csq |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
            self.csq_stage_splice(cur, s)?;
            return Ok(SPLICE_OUTSIDE);
        }
        if s.ref_end < ex_beg || (s.ref_end == ex_beg && !s.check_region_beg) {
            // fully outside, before the exon
            if s.check_utr {
                let e = s.ref_end.wrapping_sub(1);
                if self.splice_utr(cur, s, e, e)? != 0 {
                    return Ok(SPLICE_OUTSIDE); // overlaps utr
                }
            }
            if !s.check_region_beg {
                return Ok(SPLICE_OUTSIDE);
            }
            if s.set_refalt {
                s.build_hap(ex_beg.wrapping_sub(INTRON), INTRON as i32);
            }
            if s.ref_end > ex_beg.wrapping_sub(INTRON) && s.ref_beg < ex_beg.wrapping_sub(DONOR) {
                s.csq |= CSQ_SPLICE_REGION;
                if s.refalt()
                    .is_some_and(|(r, a)| strncmp_eq(r, 0, a, 0, INTRON as usize))
                {
                    s.csq |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
            if s.ref_end > ex_beg.wrapping_sub(DONOR) {
                if s.check_donor && s.tr_strand == Strand::Rev {
                    s.csq |= CSQ_SPLICE_DONOR;
                }
                if s.check_acceptor && s.tr_strand == Strand::Fwd {
                    s.csq |= CSQ_SPLICE_ACCEPTOR;
                }
                let off = (INTRON - DONOR) as usize;
                if s.refalt()
                    .is_some_and(|(r, a)| strncmp_eq(r, off, a, off, DONOR as usize))
                {
                    s.csq |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
            self.csq_stage_splice(cur, s)?;
            return Ok(SPLICE_OUTSIDE);
        }
        // overlaps the exon or inside the exon
        if s.ref_beg <= ex_beg.wrapping_add(2) {
            // in the first 3bp
            if s.check_region_beg {
                s.csq |= CSQ_SPLICE_REGION;
            }
            if s.tr_strand == Strand::Fwd {
                if s.check_start {
                    s.csq |= CSQ_START_LOST;
                }
            } else if s.check_stop {
                s.csq |= CSQ_STOP_LOST;
            }
        }
        if s.ref_end > ex_end.wrapping_sub(2) {
            if s.check_region_end {
                s.csq |= CSQ_SPLICE_REGION;
            }
            if s.tr_strand == Strand::Rev {
                if s.check_start {
                    s.csq |= CSQ_START_LOST;
                }
            } else if s.check_stop {
                s.csq |= CSQ_STOP_LOST;
            }
        }
        if s.set_refalt {
            // Make sure the variant will not end up left aligned to avoid overlapping records.
            if s.ref_beg < s.vcf_pos as u32 {
                // caused by too much trimming from right
                let dlen = (s.vcf_pos as u32).wrapping_sub(s.ref_beg) as i32;
                s.tbeg += dlen;
                if s.tbeg + s.tend == s.vcf_rlen {
                    s.tend -= dlen;
                }
                s.ref_beg = s.vcf_pos as u32;
            }
            if s.ref_end == ex_beg {
                s.tend -= 1; // prevent zero-length ref allele
            }
            let (b, l) = (s.ref_beg, s.vcf_alen - s.tend - s.tbeg + 1);
            s.build_hap(b, l);
            s.vcf_rlen -= s.tbeg + s.tend - 1;
            if s.vcf_rlen >= 0 && s.kref.len() > s.vcf_rlen as usize {
                s.kref.truncate(s.vcf_rlen as usize);
            }
        }
        self.csq_stage_splice(cur, s)?;
        Ok(SPLICE_INSIDE)
    }

    /// Whether a deletion at the start codon can be shifted to keep it (bcftools'
    /// `shifted_del_synonymous`).
    fn shifted_del_synonymous(&mut self, s: &Splice<'_>, ex_beg: u32, ex_end: u32) -> bool {
        // We know the VCF record overlaps the exon, but does it overlap the start codon?
        if s.tr_strand == Strand::Rev && (s.vcf_pos + s.vcf_rlen + 2) as u32 <= ex_end {
            return false;
        }
        if s.tr_strand == Strand::Fwd && s.vcf_pos as u32 >= ex_beg.wrapping_add(3) {
            return false;
        }
        let ref_len = strlen(s.vcf_ref) as i32;
        let alt_len = strlen(s.vcf_alt) as i32;
        if ref_len <= alt_len {
            return false;
        }
        let ndel = ref_len - alt_len;
        let pad = N_REF_PAD as i32;
        let tref = s.tref();
        let ref_allele = &s.vcf_ref[..ref_len as usize];
        let ptr_vcf = &ref_allele[alt_len as usize..]; // the first deleted base
        let cmp_from = |from: i64| -> bool {
            // the deleted sequence can be replaced
            ptr_vcf
                .iter()
                .enumerate()
                .all(|(i, &c)| at(tref, from + i as i64) == c)
        };
        match s.tr_strand {
            Strand::Rev => {
                let vcf_ref_end = s.vcf_pos + ref_len - 1;
                let tr_ref_end = s.tr_end as i32 + pad;
                if vcf_ref_end + ndel > tr_ref_end {
                    self.warn_small_ref_padding(s);
                    return false;
                }
                let from = i64::from(pad) + i64::from(vcf_ref_end) + 1 - i64::from(s.tr_beg);
                cmp_from(from)
            }
            Strand::Fwd => {
                let vcf_block_beg = s.vcf_pos + ref_len - 2 * ndel;
                if vcf_block_beg < 0 {
                    return false;
                }
                if ((pad + vcf_block_beg) as u32) < ex_beg {
                    self.warn_small_ref_padding(s);
                    return false;
                }
                let from = i64::from(pad) + i64::from(vcf_block_beg) - i64::from(s.tr_beg);
                cmp_from(from)
            }
            Strand::Unknown => true,
        }
    }

    fn warn_small_ref_padding(&mut self, s: &Splice<'_>) {
        if !self.warned.small_ref_padding {
            eprintln!(
                "Warning: Could not verify synonymous start/stop at {}:{} due to small N_REF_PAD. (Improve me?)",
                String::from_utf8_lossy(&self.chr_vcf),
                s.vcf_pos + 1
            );
            self.warned.small_ref_padding = true;
        }
    }

    /// A deletion (bcftools' `splice_csq_del`).
    fn splice_csq_del(
        &mut self,
        cur: &CurRec,
        s: &mut Splice<'_>,
        ex_beg: u32,
        ex_end: u32,
    ) -> anyhow::Result<i32> {
        if s.check_start && self.shifted_del_synonymous(s, ex_beg, ex_end) {
            // synonymous start
            s.csq |= CSQ_START_RETAINED;
            return Ok(SPLICE_OVERLAP);
        }
        // coordinates that matter for consequences, eg AC>ACG trimmed to C>CG
        s.ref_beg = (s.vcf_pos + s.tbeg - 1) as u32; // 1b before the deleted base
        s.ref_end = (s.vcf_pos + s.vcf_rlen - s.tend - 1) as u32; // the last deleted base

        if s.ref_beg.wrapping_add(1) < ex_beg {
            // the part before the exon; ref_beg is off by -1
            if s.check_region_beg {
                let mut csq = 0;
                if s.check_utr {
                    csq = self.splice_utr(cur, s, s.ref_beg, ex_beg.wrapping_sub(1))?;
                }
                if csq == 0 {
                    if s.set_refalt {
                        // filling from the right (from the left fails ENST00000341065/frame3)
                        s.build_hap(ex_beg.wrapping_sub(INTRON), INTRON as i32);
                    }
                    if s.ref_end >= ex_beg.wrapping_sub(INTRON)
                        && s.ref_beg < ex_beg.wrapping_sub(DONOR)
                    {
                        s.csq |= CSQ_SPLICE_REGION;
                        if s.refalt()
                            .is_some_and(|(r, a)| strncmp_eq(r, 0, a, 0, INTRON as usize))
                        {
                            s.csq |= CSQ_SYNONYMOUS_VARIANT;
                        }
                    }
                    if s.ref_end >= ex_beg.wrapping_sub(DONOR) {
                        if s.check_donor && s.tr_strand == Strand::Rev {
                            s.csq |= CSQ_SPLICE_DONOR;
                        }
                        if s.check_acceptor && s.tr_strand == Strand::Fwd {
                            s.csq |= CSQ_SPLICE_ACCEPTOR;
                        }
                        let noff = (INTRON - DONOR) as usize;
                        if s.refalt().is_some_and(|(r, a)| {
                            noff < r.len()
                                && noff < a.len()
                                && strncmp_eq(r, noff, a, noff, DONOR as usize)
                        }) {
                            s.csq |= CSQ_SYNONYMOUS_VARIANT;
                        }
                    }
                }
            }
            if s.ref_end >= ex_beg {
                s.tbeg = s.ref_beg.wrapping_sub(s.vcf_pos as u32).wrapping_add(1) as i32;
                s.ref_beg = ex_beg.wrapping_sub(1);
                if s.tbeg + s.tend == s.vcf_alen {
                    // the deletion overlaps ex_beg and cannot be easily realigned to the right
                    if s.tend == 0 {
                        s.csq |= CSQ_CODING_SEQUENCE;
                        return Ok(SPLICE_OVERLAP);
                    }
                    s.tend -= 1;
                }
            }
        }
        if ex_end < s.ref_end {
            // the part after the exon
            if s.check_region_end {
                let mut csq = 0;
                if s.check_utr {
                    csq = self.splice_utr(cur, s, ex_end.wrapping_add(1), s.ref_end)?;
                }
                if csq == 0 {
                    if s.set_refalt {
                        // ref,alt positioned at the first intron base
                        s.build_hap(ex_end.wrapping_add(1), INTRON as i32);
                    }
                    if s.ref_beg < ex_end.wrapping_add(INTRON)
                        && s.ref_end > ex_end.wrapping_add(DONOR)
                    {
                        s.csq |= CSQ_SPLICE_REGION;
                        if s.refalt()
                            .is_some_and(|(r, a)| strncmp_eq(r, 0, a, 0, INTRON as usize))
                        {
                            s.csq |= CSQ_SYNONYMOUS_VARIANT;
                        }
                    }
                    if s.ref_beg < ex_end.wrapping_add(DONOR) {
                        if s.check_donor && s.tr_strand == Strand::Fwd {
                            s.csq |= CSQ_SPLICE_DONOR;
                        }
                        if s.check_acceptor && s.tr_strand == Strand::Rev {
                            s.csq |= CSQ_SPLICE_ACCEPTOR;
                        }
                        let off = (INTRON - DONOR) as usize;
                        if s.refalt()
                            .is_some_and(|(r, a)| strncmp_eq(r, off, a, off, DONOR as usize))
                        {
                            s.csq |= CSQ_SYNONYMOUS_VARIANT;
                        }
                    }
                }
            }
            if s.ref_beg < ex_end {
                s.tend =
                    s.vcf_rlen - (s.ref_end.wrapping_sub(s.vcf_pos as u32).wrapping_add(1) as i32);
                s.ref_end = ex_end;
            }
        }
        if s.ref_end < ex_beg || s.ref_beg >= ex_end {
            self.csq_stage_splice(cur, s)?;
            return Ok(SPLICE_OUTSIDE);
        }
        if s.ref_beg < ex_beg.wrapping_add(2) {
            // ref_beg is off by -1
            if s.check_region_beg {
                s.csq |= CSQ_SPLICE_REGION;
            }
            if s.tr_strand == Strand::Fwd {
                if s.check_start {
                    s.csq |= CSQ_START_LOST;
                }
            } else if s.check_stop {
                s.csq |= CSQ_STOP_LOST;
            }
        }
        if s.ref_end > ex_end.wrapping_sub(3) {
            if s.check_region_end {
                s.csq |= CSQ_SPLICE_REGION;
            }
            if s.tr_strand == Strand::Rev {
                if s.check_start {
                    s.csq |= CSQ_START_LOST;
                }
            } else if s.check_stop {
                s.csq |= CSQ_STOP_LOST;
            }
        }
        if s.set_refalt {
            if s.tbeg > 0 {
                s.tbeg -= 1;
            }
            if s.vcf_rlen > s.tbeg + s.tend && s.vcf_alen > s.tbeg + s.tend {
                s.vcf_rlen -= s.tbeg + s.tend;
                s.vcf_alen -= s.tbeg + s.tend;
            }
            s.kref.clear();
            s.kref
                .put(bytes(s.vcf_ref, i64::from(s.tbeg), i64::from(s.vcf_rlen)));
            s.kalt.clear();
            s.kalt
                .put(bytes(s.vcf_alt, i64::from(s.tbeg), i64::from(s.vcf_alen)));
            if (s.ref_beg.wrapping_add(1) < ex_beg && s.ref_end >= ex_beg)
                || (s.ref_beg.wrapping_add(1) < ex_end && s.ref_end >= ex_end)
            {
                // ENST00000409523/long-overlapping-del.vcf
                let ref_beg = (s.ref_beg as i32).wrapping_add(s.kalt.len() as i32 - 1);
                if (ref_beg as u32) < s.ref_end {
                    s.csq |= if !s.ref_end.wrapping_sub(ref_beg as u32).is_multiple_of(3) {
                        CSQ_FRAMESHIFT_VARIANT
                    } else {
                        CSQ_INFRAME_DELETION
                    };
                }
                return Ok(SPLICE_OVERLAP);
            }
        }
        self.csq_stage_splice(cur, s)?;
        Ok(SPLICE_INSIDE)
    }

    /// A substitution of one or more bases (bcftools' `splice_csq_mnp`).
    fn splice_csq_mnp(
        &mut self,
        cur: &CurRec,
        s: &mut Splice<'_>,
        ex_beg: u32,
        ex_end: u32,
    ) -> anyhow::Result<i32> {
        // not a real variant, can be ignored: eg ACGT>ACGT
        if s.tbeg + s.tend == s.vcf_rlen {
            return Ok(SPLICE_VAR_REF);
        }
        s.ref_beg = (s.vcf_pos + s.tbeg) as u32;
        s.ref_end = (s.vcf_pos + s.vcf_rlen - s.tend - 1) as u32;

        let mut ret = SPLICE_INSIDE;
        if s.ref_beg < ex_beg {
            // the part before the exon
            if s.check_region_beg {
                let mut csq = 0;
                if s.check_utr {
                    csq = self.splice_utr(cur, s, s.ref_beg, ex_beg.wrapping_sub(1))?;
                }
                if csq == 0 {
                    if s.ref_end >= ex_beg.wrapping_sub(INTRON)
                        && s.ref_beg < ex_beg.wrapping_sub(DONOR)
                    {
                        s.csq |= CSQ_SPLICE_REGION;
                    }
                    if s.ref_end >= ex_beg.wrapping_sub(DONOR) {
                        if s.check_donor && s.tr_strand == Strand::Rev {
                            s.csq |= CSQ_SPLICE_DONOR;
                        }
                        if s.check_acceptor && s.tr_strand == Strand::Fwd {
                            s.csq |= CSQ_SPLICE_ACCEPTOR;
                        }
                    }
                }
            }
            if s.ref_end >= ex_beg {
                s.tbeg = s.ref_beg.wrapping_sub(s.vcf_pos as u32) as i32;
                s.ref_beg = ex_beg;
                ret = SPLICE_OVERLAP;
            }
        }
        if ex_end < s.ref_end {
            // the part after the exon
            if s.check_region_end {
                let mut csq = 0;
                if s.check_utr {
                    csq = self.splice_utr(cur, s, ex_end.wrapping_add(1), s.ref_end)?;
                }
                if csq == 0 {
                    if s.ref_beg <= ex_end.wrapping_add(INTRON)
                        && s.ref_end > ex_end.wrapping_add(DONOR)
                    {
                        s.csq |= CSQ_SPLICE_REGION;
                    }
                    if s.ref_beg <= ex_end.wrapping_add(DONOR) {
                        if s.check_donor && s.tr_strand == Strand::Fwd {
                            s.csq |= CSQ_SPLICE_DONOR;
                        }
                        if s.check_acceptor && s.tr_strand == Strand::Rev {
                            s.csq |= CSQ_SPLICE_ACCEPTOR;
                        }
                    }
                }
            }
            if s.ref_beg <= ex_end {
                s.tend =
                    s.vcf_rlen - (s.ref_end.wrapping_sub(s.vcf_pos as u32).wrapping_add(1) as i32);
                s.ref_end = ex_end;
                ret = SPLICE_OVERLAP;
            }
        }
        if s.ref_end < ex_beg || s.ref_beg > ex_end {
            self.csq_stage_splice(cur, s)?;
            return Ok(SPLICE_OUTSIDE);
        }
        if s.ref_beg < ex_beg.wrapping_add(3) {
            if s.check_region_beg {
                s.csq |= CSQ_SPLICE_REGION;
            }
            if s.tr_strand == Strand::Fwd {
                if s.check_start {
                    s.csq |= CSQ_START_LOST;
                }
            } else if s.tr_strand == Strand::Rev && s.check_stop {
                s.csq |= CSQ_STOP_LOST;
            }
        }
        if s.ref_end > ex_end.wrapping_sub(3) {
            if s.check_region_end {
                s.csq |= CSQ_SPLICE_REGION;
            }
            if s.tr_strand == Strand::Rev {
                if s.check_start {
                    s.csq |= CSQ_START_LOST;
                }
            } else if s.tr_strand == Strand::Fwd && s.check_stop {
                s.csq |= CSQ_STOP_LOST;
            }
        }
        if s.set_refalt {
            s.vcf_rlen -= s.tbeg + s.tend;
            s.vcf_alen -= s.tbeg + s.tend;
            s.kref.clear();
            s.kref
                .put(bytes(s.vcf_ref, i64::from(s.tbeg), i64::from(s.vcf_rlen)));
            s.kalt.clear();
            s.kalt
                .put(bytes(s.vcf_alt, i64::from(s.tbeg), i64::from(s.vcf_alen)));
        }
        self.csq_stage_splice(cur, s)?;
        Ok(ret)
    }

    /// The consequences of a variant near the exon `ex_beg..=ex_end` (bcftools'
    /// `splice_csq`): one of the `SPLICE_*` codes.
    pub(super) fn splice_csq(
        &mut self,
        cur: &CurRec,
        s: &mut Splice<'_>,
        ex_beg: u32,
        ex_end: u32,
    ) -> anyhow::Result<i32> {
        s.vcf_alen = strlen(s.vcf_alt) as i32;

        let mut rlen1 = s.vcf_rlen - 1;
        let mut alen1 = s.vcf_alen - 1;
        let mut i = 0;
        // trim from the right, then from the left
        while i <= rlen1 && i <= alen1 {
            if at(s.vcf_ref, i64::from(rlen1 - i)) != at(s.vcf_alt, i64::from(alen1 - i)) {
                break;
            }
            i += 1;
        }
        s.tend = i;
        rlen1 -= i;
        alen1 -= i;
        i = 0;
        while i <= rlen1 && i <= alen1 {
            if at(s.vcf_ref, i64::from(i)) != at(s.vcf_alt, i64::from(i)) {
                break;
            }
            i += 1;
        }
        s.tbeg = i;

        let mut rtrim = s.vcf_rlen - s.tbeg - s.tend;
        let mut atrim = s.vcf_alen - s.tbeg - s.tend;
        if at(s.vcf_alt, 0) == b'<' {
            rtrim = 0;
            atrim = 0;
        }
        if s.vcf_rlen == s.vcf_alen {
            return self.splice_csq_mnp(cur, s, ex_beg, ex_end);
        }
        if rtrim > 1 && atrim > 1 {
            // complex (bcftools' `splice_csq_complex`)
            s.csq |= if s.vcf_rlen > s.vcf_alen {
                CSQ_TRUNCATION
            } else {
                CSQ_ELONGATION
            };
            return self.splice_csq_mnp(cur, s, ex_beg, ex_end);
        }
        if s.vcf_rlen < s.vcf_alen {
            return self.splice_csq_ins(cur, s, ex_beg, ex_end);
        }
        if s.vcf_rlen > s.vcf_alen {
            return self.splice_csq_del(cur, s, ex_beg, ex_end);
        }
        Ok(SPLICE_VAR_REF)
    }
}
