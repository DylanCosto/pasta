//! Haplotype trees of coding transcripts (bcftools' `hap_node_t`, `tscript_t`, `hap_init`,
//! `hap_finalize`, and `hap_add_csq`): each transcript collects the variants of every
//! haplotype in a tree whose leaves are the haplotypes; when the transcript is done, each
//! haplotype's spliced sequence is translated and broken into independent consequences.

use anyhow::bail;

use super::{
    Csq, CurRec, N_REF_PAD,
    splice::{SPLICE_OUTSIDE, SPLICE_OVERLAP, SPLICE_VAR_REF, Splice},
    translate::{Fill, KStr, SplicedSeq, cds_translate},
    types::{
        CSQ_CODING_SEQUENCE, CSQ_COMPOUND, CSQ_FRAMESHIFT_VARIANT, CSQ_INCOMPLETE_CDS,
        CSQ_INFRAME_ALTERING, CSQ_INFRAME_DELETION, CSQ_INFRAME_INSERTION, CSQ_MISSENSE_VARIANT,
        CSQ_PRINTED_UPSTREAM, CSQ_START_LOST, CSQ_START_RETAINED, CSQ_STOP_GAINED, CSQ_STOP_LOST,
        CSQ_STOP_RETAINED, CSQ_SYNONYMOUS_VARIANT, CSQ_UPSTREAM_STOP, CsqT, RecRef, VCsq, VStr,
        gene_ptr,
    },
};
use crate::gff::{Strand, TRIM_3PRIME, TRIM_5PRIME};

/// The kind of a haplotype tree node (bcftools' `HAP_*`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum HapType {
    /// A coding variant.
    #[default]
    Cds,
    /// The root.
    Root,
    /// A start, stop, or splice variant.
    Sss,
}

/// A node of a haplotype tree (bcftools' `hap_node_t`).
#[derive(Clone, Debug, Default)]
pub(super) struct HapNode {
    /// The CDS segment from the parent node to this one, ending with the alternate allele.
    pub seq: Option<Vec<u8>>,
    /// The variant as `REF>ALT`.
    pub var: Vec<u8>,
    pub ty: HapType,
    pub csq: u32,
    /// The alternate minus the reference length.
    pub dlen: i32,
    /// The variant's position (0-based).
    pub rbeg: u32,
    pub rlen: i32,
    /// The variant's position on the spliced reference transcript.
    pub sbeg: u32,
    /// The CDS the variant overlaps.
    pub icds: u32,
    /// The children (indexes into [`TrAux::nodes`]).
    pub child: Vec<u32>,
    /// The record whose alleles `cur_child` maps.
    pub cur_rec: Option<u64>,
    pub rec: Option<RecRef>,
    pub vcf_ial: i32,
    /// The number of haplotypes ending in this node.
    pub nend: u32,
    /// The child of each allele of `cur_rec`, or -1.
    pub cur_child: Vec<i32>,
    /// The haplotype's consequences, one per record.
    pub csq_list: Vec<CsqT>,
}

/// The state of a coding transcript being annotated (bcftools' `tscript_t`).
#[derive(Debug, Default)]
pub(super) struct TrAux {
    /// The reference sequence padded with [`N_REF_PAD`] bases on both ends.
    pub ref_: Vec<u8>,
    /// The spliced reference sequence, padded likewise.
    pub sref: Option<Vec<u8>>,
    /// The haplotype tree; the root is node 0 when present.
    pub nodes: Vec<HapNode>,
    pub has_root: bool,
    /// The leaf of each haplotype, two per sample.
    pub hap: Vec<Option<u32>>,
}

impl TrAux {
    /// The root node (index 0).
    pub(super) const ROOT: u32 = 0;
}

/// A node on the traversal stack (bcftools' `hstack_t`).
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct HStack {
    node: u32,
    ichild: i32,
    dlen: i32,
    slen: usize,
}

/// The traversal of a transcript's haplotype tree (bcftools' `hap_t`).
#[derive(Debug, Default)]
pub(super) struct Hap {
    stack: Vec<HStack>,
    tr: u32,
    /// The spliced haplotype sequence on the reference strand.
    sseq: Vec<u8>,
    pub tseq: KStr,
    pub tref: KStr,
    pub tseq_stop: KStr,
    pub tref_stop: KStr,
    sbeg: u32,
    upstream_stop: bool,
    /// The length of the spliced reference, padding included (`tscript_t.nsref`).
    nsref: i64,
}

/// Builds the spliced reference of a transcript (bcftools' `tscript_splice_ref`).
///
/// # Errors
///
/// A CDS whose length wrapped around (it ends before it starts in the GFF, or is shorter
/// than its phase), where bcftools copies gigabytes and crashes.
pub(super) fn splice_ref(aux: &TrAux, tr_beg: u32, cds: &[(u32, u32)]) -> anyhow::Result<Vec<u8>> {
    if let Some(&(beg, l)) = cds.iter().find(|&&(_, l)| i32::try_from(l).is_err()) {
        bail!(
            "Error: the CDS at {} has a negative length ({}), check the GFF",
            beg + 1,
            l as i32
        );
    }
    let len: usize = cds.iter().map(|&(_, l)| l as usize).sum();
    let mut sref = Vec::with_capacity(len + 2 * N_REF_PAD);
    let get =
        |from: usize, n: usize| (from..from + n).map(|i| aux.ref_.get(i).copied().unwrap_or(0));
    sref.extend(get((cds[0].0 - tr_beg) as usize, N_REF_PAD));
    for &(beg, l) in cds {
        sref.extend(get(N_REF_PAD + (beg - tr_beg) as usize, l as usize));
    }
    let (lbeg, llen) = cds[cds.len() - 1];
    sref.extend(get(
        N_REF_PAD + (lbeg - tr_beg) as usize + llen as usize,
        N_REF_PAD,
    ));
    Ok(sref)
}

/// What [`Csq::hap_init`] did with a variant.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum HapInit {
    /// The node was filled.
    Added,
    /// The variant overlaps the previous one of the haplotype.
    Overlapping,
    /// Not a coding variant (intronic, or REF equals ALT).
    Silent,
}

/// Appends `n` bases of `s` from `from` (0 outside it), as `kputsn_` does.
fn put_bases(out: &mut Vec<u8>, s: &[u8], from: i64, n: i64) {
    for i in 0..n.max(0) {
        let k = from + i;
        out.push(
            usize::try_from(k)
                .ok()
                .and_then(|k| s.get(k))
                .copied()
                .unwrap_or(0),
        );
    }
}

impl Csq {
    /// Fills a haplotype node for allele `ial` of the current record in CDS `cds_id`
    /// (bcftools' `hap_init`).
    pub(super) fn hap_init(
        &mut self,
        cur: &CurRec,
        aux: &TrAux,
        parent: &HapNode,
        child: &mut HapNode,
        cds_id: u32,
        ial: usize,
    ) -> anyhow::Result<HapInit> {
        let cds = self.gff.cds[cds_id as usize];
        let (tr_beg, strand, trim, ncds) = {
            let tr = &self.gff.transcripts[cds.tr as usize];
            (tr.beg, tr.strand, tr.trim, tr.cds.len() as u32)
        };
        child.icds = cds.icds;
        child.vcf_ial = ial as i32;

        let mut splice = Splice::new(cur, &self.gff, cds.tr, Some(&aux.ref_));
        splice.set_alt(cur, ial);
        splice.check_acceptor = true;
        splice.check_donor = true;
        splice.set_refalt = true;
        splice.check_utr = true;
        if trim & TRIM_5PRIME == 0 {
            match strand {
                Strand::Fwd => splice.check_start = child.icds == 0,
                Strand::Rev => splice.check_start = child.icds == ncds - 1,
                Strand::Unknown => {}
            }
        }
        if trim & TRIM_3PRIME == 0 {
            match strand {
                Strand::Fwd => splice.check_stop = child.icds == ncds - 1,
                Strand::Rev => splice.check_stop = child.icds == 0,
                Strand::Unknown => {}
            }
        }
        if splice.check_start {
            // Do not check starts in incomplete CDS, defined as not starting with M.
            let at = |i: i64| {
                usize::try_from(i)
                    .ok()
                    .and_then(|i| aux.ref_.get(i))
                    .copied()
                    .unwrap_or(0)
            };
            let off = N_REF_PAD as i64 + i64::from(cds.beg) - i64::from(tr_beg);
            match strand {
                Strand::Fwd => {
                    if self.gencode.dna2stop([at(off), at(off + 1), at(off + 2)]) != b'M' {
                        splice.check_start = false;
                    }
                }
                Strand::Rev => {
                    let o = off + i64::from(cds.len) - 3;
                    if self.gencode.cdna2stop([at(o), at(o + 1), at(o + 2)]) != b'M' {
                        splice.check_start = false;
                    }
                }
                Strand::Unknown => {}
            }
        }
        if child.icds != 0 {
            splice.check_region_beg = true;
        }
        if child.icds != ncds - 1 {
            splice.check_region_end = true;
        }
        let ret = self.splice_csq(
            cur,
            &mut splice,
            cds.beg,
            cds.beg.wrapping_add(cds.len).wrapping_sub(1),
        )?;
        if ret == SPLICE_VAR_REF {
            return Ok(HapInit::Silent);
        }
        if ret == SPLICE_OUTSIDE || ret == SPLICE_OVERLAP {
            if splice.csq == 0 {
                return Ok(HapInit::Silent); // fully intronic, no csq
            }
            // splice_region/acceptor/donor
            child.seq = None;
            child.sbeg = 0;
            child.rbeg = cur.pos as u32;
            child.rlen = 0;
            child.dlen = 0;
            child.var = var_string(cur, ial);
            child.ty = HapType::Sss;
            child.csq = splice.csq;
            child.rec = Some(cur.rec_ref());
            return Ok(HapInit::Added);
        }
        // synonymous&splice,frame could become synonymous&frame,splice
        splice.csq &= !CSQ_SYNONYMOUS_VARIANT;

        let mut dbeg = 0usize;
        if splice.ref_beg < cds.beg {
            // The record overlaps the exon boundary, but the variant should fit inside.
            dbeg = (cds.beg - splice.ref_beg) as usize;
            let l = splice.kref.len();
            splice.kref.set_len(l.wrapping_sub(dbeg));
            splice.ref_beg = cds.beg;
        }

        let tr_cds = |i: usize| {
            let c = &self.gff.cds[self.gff.transcripts[cds.tr as usize].cds[i] as usize];
            (c.beg, c.len)
        };
        let tref = &aux.ref_;
        let pad = N_REF_PAD as i64;
        let mut str = Vec::new();
        // bcftools asserts that the parent is not HAP_SSS (see hap_finalize); such a parent,
        // like the root, contributes no sequence.
        if parent.ty == HapType::Cds {
            let mut i = parent.icds as usize;
            if i != cds.icds as usize {
                // the variant is on a new exon, finish up the previous
                let len = (i64::from(tr_cds(i).1) - i64::from(parent.rbeg) - i64::from(parent.rlen)
                    + i64::from(tr_cds(i).0)) as i32;
                if len > 0 {
                    let from =
                        pad + i64::from(parent.rbeg) + i64::from(parent.rlen) - i64::from(tr_beg);
                    put_bases(&mut str, tref, from, i64::from(len));
                }
            }
            // append any skipped non-variant exons
            i += 1;
            while i < cds.icds as usize {
                let from = pad + i64::from(tr_cds(i).0) - i64::from(tr_beg);
                put_bases(&mut str, tref, from, i64::from(tr_cds(i).1));
                i += 1;
            }
            if parent.icds == child.icds {
                let len = (i64::from(splice.ref_beg)
                    - i64::from(parent.rbeg)
                    - i64::from(parent.rlen)) as i32;
                if len < 0 {
                    return Ok(HapInit::Overlapping); // overlapping variants
                }
                let from =
                    pad + i64::from(parent.rbeg) + i64::from(parent.rlen) - i64::from(tr_beg);
                put_bases(&mut str, tref, from, i64::from(len));
            } else {
                let from = pad + i64::from(cds.beg) - i64::from(tr_beg);
                put_bases(
                    &mut str,
                    tref,
                    from,
                    i64::from(splice.ref_beg) - i64::from(cds.beg),
                );
            }
        }
        // kputs(kalt.s + dbeg): up to the NUL
        let mut k = dbeg;
        loop {
            let c = splice.kalt.at(k);
            if c == 0 {
                break;
            }
            str.push(c);
            k += 1;
        }

        child.seq = Some(str);
        child.sbeg = cds.pos.wrapping_add(splice.ref_beg - cds.beg);
        child.rbeg = splice.ref_beg;
        child.rlen = splice.kref.len() as i32;
        child.ty = HapType::Cds;
        child.rec = Some(cur.rec_ref());
        child.csq = splice.csq;

        let rlen = cur.allele_str(0).len() as i32;
        let alen = cur.allele_str(ial).len() as i32;
        child.dlen = alen - rlen;
        child.var = var_string(cur, ial);

        // The whole CDS is modified or deleted.
        if i64::from(child.rbeg) + i64::from(child.rlen) > i64::from(cds.beg) + i64::from(cds.len) {
            child.ty = HapType::Sss;
            if child.csq == 0 {
                child.csq |= CSQ_CODING_SEQUENCE;
            }
        }
        Ok(HapInit::Added)
    }

    /// Translates the haplotypes of a finished transcript and records their consequences
    /// (bcftools' `hap_finalize`).
    pub(super) fn hap_finalize(&mut self, tr: u32, aux: &mut TrAux) -> anyhow::Result<()> {
        let mut hap = std::mem::take(&mut self.hap);
        hap.tr = tr;
        let result = self.hap_finalize_inner(&mut hap, aux);
        self.hap = hap;
        result
    }

    fn hap_finalize_inner(&mut self, hap: &mut Hap, aux: &mut TrAux) -> anyhow::Result<()> {
        let (strand, tr_beg) = {
            let t = &self.gff.transcripts[hap.tr as usize];
            (t.strand, t.beg)
        };
        if aux.sref.is_none() {
            let cds: Vec<(u32, u32)> = self.gff.transcripts[hap.tr as usize]
                .cds
                .iter()
                .map(|&c| (self.gff.cds[c as usize].beg, self.gff.cds[c as usize].len))
                .collect();
            aux.sref = Some(splice_ref(aux, tr_beg, &cds)?);
        }
        let sref = aux.sref.take().unwrap_or_default();
        hap.nsref = sref.len() as i64;

        hap.stack.clear();
        hap.stack.push(HStack {
            node: TrAux::ROOT,
            ichild: -1,
            slen: 0,
            dlen: 0,
        });
        hap.sseq.clear();
        hap.tseq.clear();
        hap.tseq_stop.clear();

        let result = self.hap_walk(hap, aux, &sref, strand);
        aux.sref = Some(sref);
        result
    }

    /// Walks the tree depth-first, building each haplotype's spliced sequence.
    fn hap_walk(
        &mut self,
        hap: &mut Hap,
        aux: &mut TrAux,
        sref: &[u8],
        strand: Strand,
    ) -> anyhow::Result<()> {
        let mut istack: isize = 0;
        while istack >= 0 {
            let is = istack as usize;
            let node_id = hap.stack[is].node;
            let nchild = aux.nodes[node_id as usize].child.len() as i32;
            hap.stack[is].ichild += 1;
            if hap.stack[is].ichild >= nchild {
                istack -= 1;
                continue;
            }
            let child = aux.nodes[node_id as usize].child[hap.stack[is].ichild as usize];
            istack += 1;
            let is = istack as usize;
            if hap.stack.len() <= is {
                hap.stack.push(HStack::default());
            }
            let prev = hap.stack[is - 1];
            hap.sseq.truncate(prev.slen);
            let cnode = &aux.nodes[child as usize];
            if cnode.ty == HapType::Cds
                && let Some(seq) = &cnode.seq
            {
                hap.sseq.extend_from_slice(seq);
            }
            hap.stack[is] = HStack {
                node: child,
                ichild: -1,
                slen: hap.sseq.len(),
                dlen: prev.dlen + cnode.dlen,
            };
            if cnode.nend == 0 {
                continue; // not a leaf node
            }
            self.hap_break(hap, aux, sref, strand, is)?;
        }
        Ok(())
    }

    /// Breaks the spliced sequence of the haplotype ending at stack position `istack` into
    /// independent consequences (the leaf part of bcftools' `hap_finalize`).
    fn hap_break(
        &mut self,
        hap: &mut Hap,
        aux: &mut TrAux,
        sref: &[u8],
        strand: Strand,
        istack: usize,
    ) -> anyhow::Result<()> {
        let pad = N_REF_PAD as i64;
        let node = hap.stack[istack].node;
        let total_dlen = i64::from(hap.stack[istack].dlen);
        let sref_m = sref.len() as i64;
        let sseq_m = sref_m - 2 * pad + total_dlen;
        hap.upstream_stop = false;

        let mut dlen: i32 = 0;
        let mut indel = false;
        hap.sbeg = node_at(aux, hap, 1).sbeg;
        // bcftools asserts that no node on the stack is HAP_SSS, but hap_add_csq turns a
        // start-lost frameshift node into one after its own leaf was done. Its descendants
        // still carry its sequence and lengths, so they are broken up as in release builds
        // (where the assert is compiled out): a part starting at that node reports its
        // splice/start/stop consequence.
        match strand {
            Strand::Fwd => {
                let mut ibeg: isize = -1;
                for i in 1..=istack {
                    let n = node_at(aux, hap, i);
                    let ndl = n.dlen;
                    dlen += ndl;
                    if ndl != 0 {
                        indel = true;
                    }
                    // This condition extends compound variants.
                    if i < istack {
                        if dlen % 3 != 0 {
                            // frameshift
                            if ibeg == -1 {
                                ibeg = i as isize;
                            }
                            continue;
                        }
                        // the last base of the current variant vs the first base of the next
                        // variant: are they in the same codon? (forward strand)
                        let mut icur = node2sbeg(aux, hap, i) as i32;
                        let inext = node2sbeg(aux, hap, i + 1) as i32;
                        if ndl > 0 {
                            icur += ndl;
                        } else if ndl < 0 {
                            icur += 1;
                        }
                        if icur / 3 == inext / 3 {
                            if ibeg == -1 {
                                ibeg = i as isize;
                            }
                            continue;
                        }
                    }
                    if ibeg < 0 {
                        ibeg = i as isize;
                    }
                    let ib = ibeg as usize;
                    let ioff = node2soff(aux, hap, ib) as i32;
                    let icur = node2sbeg(aux, hap, ib) as i32;
                    let rbeg = node2rbeg(aux, hap, ib) as i32;
                    let rend = node2rend(aux, hap, i) as i32;
                    let to = hap.stack[i].slen;
                    self.translate_pair(
                        hap,
                        sref,
                        (ioff, to),
                        sseq_m,
                        (icur, rbeg, rend),
                        strand,
                        dlen,
                    )?;
                    self.hap_add_csq(hap, aux, node, 0, ib, i, dlen, indel)?;
                    ibeg = -1;
                    dlen = 0;
                    indel = false;
                }
            }
            Strand::Rev => {
                let mut ibeg: isize = -1;
                for i in (1..=istack).rev() {
                    let n = node_at(aux, hap, i);
                    let ndl = n.dlen;
                    dlen += ndl;
                    if ndl != 0 {
                        indel = true;
                    }
                    if i > 1 {
                        if dlen % 3 != 0 {
                            if ibeg == -1 {
                                ibeg = i as isize;
                            }
                            continue;
                        }
                        // the last base of the current variant vs the first base of the next
                        // variant: are they in the same codon? (reverse strand)
                        let mut icur = (sseq_m - 1 - node2sbeg(aux, hap, i)) as i32;
                        let mut inext = (sseq_m - 1 - node2sbeg(aux, hap, i - 1)) as i32;
                        if ndl > 0 {
                            icur += ndl - 1;
                        } else if ndl < 0 {
                            icur -= ndl;
                        }
                        let pdl = node_at(aux, hap, i - 1).dlen;
                        if pdl > 0 {
                            inext -= pdl;
                        }
                        if icur / 3 == inext / 3 {
                            if ibeg == -1 {
                                ibeg = i as isize;
                            }
                            continue;
                        }
                    }
                    if ibeg < 0 {
                        ibeg = i as isize;
                    }
                    let ib = ibeg as usize;
                    let ioff = node2soff(aux, hap, i) as i32;
                    let icur = node2sbeg(aux, hap, i) as i32;
                    let rbeg = node2rbeg(aux, hap, i) as i32;
                    let rend = node2rend(aux, hap, ib) as i32;
                    let to = hap.stack[ib].slen;
                    self.translate_pair(
                        hap,
                        sref,
                        (ioff, to),
                        sseq_m,
                        (icur, rbeg, rend),
                        strand,
                        dlen,
                    )?;
                    self.hap_add_csq(hap, aux, node, sseq_m as i32, i, ib, dlen, indel)?;
                    ibeg = -1;
                    dlen = 0;
                    indel = false;
                }
            }
            Strand::Unknown => {}
        }
        Ok(())
    }

    /// Translates the haplotype part `hap.sseq[from..to]` and the reference part it replaces
    /// (the two `cds_translate` calls of bcftools' `hap_finalize`).
    #[expect(
        clippy::too_many_arguments,
        reason = "the arguments of two cds_translate calls"
    )]
    fn translate_pair(
        &self,
        hap: &mut Hap,
        sref: &[u8],
        (ioff, to): (i32, usize),
        sseq_m: i64,
        (icur, rbeg, rend): (i32, i32, i32),
        strand: Strand,
        dlen: i32,
    ) -> anyhow::Result<()> {
        let pad = N_REF_PAD as i64;
        let mut fill = if dlen % 3 != 0 { Fill::ToEnd } else { Fill::No };
        // alt
        let alt: &[u8] = if hap.sseq.is_empty() {
            // splice site overlap, see #1475227917
            fill = Fill::No;
            &[]
        } else {
            usize::try_from(ioff)
                .ok()
                .and_then(|from| hap.sseq.get(from..to))
                .unwrap_or(&[])
        };
        cds_translate(
            self.gencode,
            sref,
            &SplicedSeq {
                bases: alt,
                total: sseq_m,
            },
            icur as u32,
            rbeg as u32,
            rend as u32,
            strand,
            &mut hap.tseq,
            &mut hap.tseq_stop,
            fill,
        )?;
        // ref
        let from = (pad + i64::from(rbeg)) as usize;
        let to = (pad + i64::from(rend)) as usize;
        let rseq = sref.get(from..to.max(from)).unwrap_or(&[]);
        cds_translate(
            self.gencode,
            sref,
            &SplicedSeq {
                bases: rseq,
                total: sref.len() as i64 - 2 * pad,
            },
            rbeg as u32,
            rbeg as u32,
            rend as u32,
            strand,
            &mut hap.tref,
            &mut hap.tref_stop,
            fill,
        )
    }

    /// Checks the translated sequences for stops and amino acid changes and sets the
    /// consequence bits (the first half of bcftools' `hap_add_csq`). Returns the bits to
    /// remove and whether to leave out the amino acid change.
    fn hap_csq_type(
        hap: &mut Hap,
        aux: &mut TrAux,
        csq: &mut CsqT,
        (ibeg, iend): (usize, usize),
        dlen: i32,
        indel: bool,
    ) -> (u32, bool) {
        let mut rm_csq = 0u32;
        let mut suppress_aa = false;
        for i in ibeg..=iend {
            csq.ty.ty |= node_at(aux, hap, i).csq & CSQ_COMPOUND;
        }
        if dlen == 0 && indel {
            csq.ty.ty |= CSQ_INFRAME_ALTERING;
        }
        let has_upstream_stop = hap.upstream_stop;
        if node_at(aux, hap, ibeg).ty != HapType::Sss {
            // check for truncating stops
            if let Some(i) = (0..hap.tref_stop.len()).find(|&i| hap.tref_stop.at(i) == b'*') {
                hap.tref.truncate(i + 1);
                hap.tref_stop.truncate(i + 1);
            }
            if let Some(i) = (0..hap.tseq_stop.len()).find(|&i| hap.tseq_stop.at(i) == b'*') {
                hap.tseq.truncate(i + 1);
                hap.tseq_stop.truncate(i + 1);
                hap.upstream_stop = true;
            }
            let ref_last = hap.tref_stop.at(hap.tref_stop.len().wrapping_sub(1));
            let seq_last = hap.tseq_stop.at(hap.tseq_stop.len().wrapping_sub(1));
            if csq.ty.ty & CSQ_STOP_LOST != 0 {
                if ref_last == b'*' && ref_last == seq_last {
                    rm_csq |= CSQ_STOP_LOST;
                    csq.ty.ty |= CSQ_STOP_RETAINED;
                } else if ref_last != b'*' {
                    // CDS 3' incomplete, can also be missense or a change to a stop
                    if seq_last == b'*' {
                        rm_csq |= CSQ_STOP_GAINED;
                        csq.ty.ty |= CSQ_STOP_RETAINED;
                    } else {
                        csq.ty.ty |= CSQ_INCOMPLETE_CDS;
                    }
                }
            }
            if csq.ty.ty & CSQ_START_LOST != 0 && ref_last == b'M' && ref_last == seq_last {
                rm_csq |= CSQ_START_LOST;
                csq.ty.ty |= CSQ_START_RETAINED;
            }
            if dlen != 0 {
                if dlen % 3 != 0 {
                    csq.ty.ty |= CSQ_FRAMESHIFT_VARIANT;
                } else if dlen < 0 {
                    csq.ty.ty |= CSQ_INFRAME_DELETION;
                } else {
                    csq.ty.ty |= CSQ_INFRAME_INSERTION;
                }
                if ref_last != b'*' && seq_last == b'*' {
                    csq.ty.ty |= CSQ_STOP_GAINED;
                }
            } else {
                let mut aa_change = false;
                for i in 0..hap.tref.len() {
                    if hap.tref.at(i) == hap.tseq.at(i) {
                        continue;
                    }
                    aa_change = true;
                    if hap.tref_stop.at(i) == b'*' {
                        csq.ty.ty |= CSQ_STOP_LOST;
                    } else if hap.tseq_stop.at(i) == b'*' {
                        csq.ty.ty |= CSQ_STOP_GAINED;
                    } else {
                        csq.ty.ty |= CSQ_MISSENSE_VARIANT;
                    }
                }
                if !aa_change {
                    csq.ty.ty |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
        }
        // Compound inframe variants: is the stop codon before the frame is restored?
        if ibeg != iend
            && csq.ty.ty & (CSQ_INFRAME_DELETION | CSQ_INFRAME_INSERTION | CSQ_INFRAME_ALTERING)
                != 0
            && hap.tseq_stop.at(hap.tseq_stop.len().wrapping_sub(1)) == b'*'
        {
            rm_csq |= CSQ_INFRAME_DELETION | CSQ_INFRAME_INSERTION | CSQ_INFRAME_ALTERING;
            csq.ty.ty |= CSQ_FRAMESHIFT_VARIANT | CSQ_STOP_GAINED;
        }
        if csq.ty.ty & CSQ_FRAMESHIFT_VARIANT != 0 && csq.ty.ty & CSQ_START_LOST != 0 {
            rm_csq |= CSQ_FRAMESHIFT_VARIANT;
            // The amino-acid prediction from the original start is misleading.
            suppress_aa = true;
            if ibeg == iend {
                let n = hap.stack[ibeg].node;
                aux.nodes[n as usize].ty = HapType::Sss;
            }
        }
        if has_upstream_stop {
            csq.ty.ty |= CSQ_UPSTREAM_STOP;
        }
        csq.ty.ty &= !rm_csq;
        (rm_csq, suppress_aa)
    }

    /// Appends a consequence of the haplotype part `ibeg..=iend` of the stack to the leaf
    /// `node` and pushes it to the records (bcftools' `hap_add_csq`).
    #[expect(
        clippy::too_many_arguments,
        reason = "line-for-line port of hap_add_csq"
    )]
    fn hap_add_csq(
        &mut self,
        hap: &mut Hap,
        aux: &mut TrAux,
        node: u32,
        tlen: i32,
        ibeg: usize,
        iend: usize,
        dlen: i32,
        indel: bool,
    ) -> anyhow::Result<()> {
        let (strand, trid, gene, biotype) = {
            let t = &self.gff.transcripts[hap.tr as usize];
            (t.strand, t.id, gene_ptr(&self.gff, t.gene), t.ty)
        };
        if !matches!(strand, Strand::Fwd | Strand::Rev) {
            bail!("csq: assertion failed in hap_add_csq: unknown strand");
        }
        let ref_node = if strand == Strand::Fwd { ibeg } else { iend };
        let ref_rec = node_rec(aux, hap, ref_node)?;

        let mut csq = CsqT {
            vrec: None,
            idx: 0,
            ty: VCsq {
                strand: strand.code() & 1,
                ty: 0,
                trid,
                vcf_ial: aux.nodes[node as usize].vcf_ial as u32,
                biotype,
                gene,
                ref_: None,
                vstr: VStr::default(),
            },
        };
        let (rm_csq, suppress_aa) =
            Self::hap_csq_type(hap, aux, &mut csq, (ibeg, iend), dlen, indel);

        if node_at(aux, hap, ibeg).ty == HapType::Sss {
            let first = node_at(aux, hap, ibeg);
            csq.ty.ty |= first.csq & !rm_csq;
            csq.ty.ref_ = first.rec;
            csq.ty.biotype = biotype;
            let rec = node_rec(aux, hap, ibeg)?;
            self.csq_push(&mut csq, rec)?;
            aux.nodes[node as usize].csq_list.push(csq);
            return Ok(());
        }

        let mut str = Vec::new();
        if suppress_aa {
            // Leave amino_acid_change empty and keep only the DNA haplotype:
            // "...|strand||dna_change".
            str.extend_from_slice(b"||");
        } else {
            // C does this arithmetic unsigned (uint32_t for the reference positions, size_t
            // for the haplotype ones), so a negative difference wraps before the division.
            let pad = N_REF_PAD as u32;
            let aa_rbeg = if strand == Strand::Fwd {
                (node2rbeg(aux, hap, ibeg) as u32 / 3 + 1) as i32
            } else {
                let rend = node2rend(aux, hap, iend) as u32;
                ((hap.nsref as u32).wrapping_sub(2 * pad).wrapping_sub(rend) / 3 + 1) as i32
            };
            let aa_sbeg = if strand == Strand::Fwd {
                (node2sbeg(aux, hap, ibeg) as u64 / 3 + 1) as i32
            } else {
                ((i64::from(tlen) as u64).wrapping_sub(node2send(hap, iend) as u64) / 3 + 1) as i32
            };
            str.push(b'|');
            str.extend_from_slice(aa_rbeg.to_string().as_bytes());
            self.kprint_aa_prediction(aa_rbeg, &hap.tref, &hap.tref_stop, &mut str);
            if csq.ty.ty & CSQ_SYNONYMOUS_VARIANT == 0 {
                str.push(b'>');
                str.extend_from_slice(aa_sbeg.to_string().as_bytes());
                self.kprint_aa_prediction(aa_sbeg, &hap.tseq, &hap.tseq_stop, &mut str);
            }
            str.push(b'|');
        }
        // the dna variant string
        for i in ibeg..=iend {
            if i > ibeg {
                str.push(b'+');
            }
            let n = node_at(aux, hap, i);
            let pos = n.rec.map_or(0, |r| r.pos);
            str.extend_from_slice(((pos + 1) as i32).to_string().as_bytes());
            str.extend_from_slice(&n.var);
        }
        csq.ty.vstr = VStr::new(str.clone());
        self.csq_push(&mut csq, ref_rec)?;
        let main_type = csq.ty.ty;
        aux.nodes[node as usize].csq_list.push(csq);

        // csq are printed at one position only for combined variants, the rest is silent
        // and references the first
        for i in ibeg..=iend {
            let ncsq = node_at(aux, hap, i).csq;
            let nrec = node_rec(aux, hap, i)?;
            if ncsq & !CSQ_COMPOUND != 0 {
                let mut tmp = CsqT {
                    vrec: None,
                    idx: 0,
                    ty: VCsq {
                        strand: strand.code() & 1,
                        ty: ncsq & !CSQ_COMPOUND & !rm_csq,
                        trid,
                        vcf_ial: 0,
                        biotype,
                        gene,
                        ref_: None,
                        vstr: VStr::new(str.clone()),
                    },
                };
                self.csq_push(&mut tmp, nrec)?;
                aux.nodes[node as usize].csq_list.push(tmp);
            }
            if i != ref_node && (main_type & CSQ_COMPOUND != 0 || ncsq & !CSQ_COMPOUND == 0) {
                let mut tmp = CsqT {
                    vrec: None,
                    idx: 0,
                    ty: VCsq {
                        strand: strand.code() & 1,
                        ty: CSQ_PRINTED_UPSTREAM | ncsq,
                        trid,
                        vcf_ial: 0,
                        biotype,
                        gene,
                        ref_: Some(ref_rec),
                        vstr: VStr::default(),
                    },
                };
                self.csq_push(&mut tmp, nrec)?;
                aux.nodes[node as usize].csq_list.push(tmp);
            }
        }
        Ok(())
    }

    /// Appends an amino acid sequence, shortened with `-B` (bcftools'
    /// `kprint_aa_prediction`).
    pub(super) fn kprint_aa_prediction(&self, beg: i32, aa: &KStr, stop: &KStr, out: &mut Vec<u8>) {
        let brief = self.opts.brief_predictions;
        if brief == 0 || aa.len() as i32 - brief < 3 {
            out.extend_from_slice(aa.as_bytes());
            return;
        }
        let mut len = aa.len();
        if stop.at(len.wrapping_sub(1)) == b'*' {
            len -= 1;
        }
        let n = len.min(brief as usize);
        out.extend_from_slice(&aa.as_bytes()[..n]);
        out.extend_from_slice(b"..");
        out.extend_from_slice((beg + len as i32).to_string().as_bytes());
    }
}

/// The node at stack position `i`.
fn node_at<'a>(aux: &'a TrAux, hap: &Hap, i: usize) -> &'a HapNode {
    &aux.nodes[hap.stack[i].node as usize]
}

/// The record of the node at stack position `i`.
fn node_rec(aux: &TrAux, hap: &Hap, i: usize) -> anyhow::Result<RecRef> {
    node_at(aux, hap, i)
        .rec
        .ok_or_else(|| anyhow::anyhow!("csq: haplotype node without a record"))
}

/// The position of a variant within the trimmed query transcript (bcftools' `node2soff`).
fn node2soff(aux: &TrAux, hap: &Hap, i: usize) -> i64 {
    let n = node_at(aux, hap, i);
    hap.stack[i].slen as i64 - (i64::from(n.rlen) + i64::from(n.dlen))
}

/// The position of a variant within the query transcript (bcftools' `node2sbeg`).
fn node2sbeg(aux: &TrAux, hap: &Hap, i: usize) -> i64 {
    i64::from(hap.sbeg) + node2soff(aux, hap, i)
}

/// bcftools' `node2send`.
fn node2send(hap: &Hap, i: usize) -> i64 {
    i64::from(hap.sbeg) + hap.stack[i].slen as i64
}

/// The position of a variant on the reference transcript (bcftools' `node2rbeg`).
fn node2rbeg(aux: &TrAux, hap: &Hap, i: usize) -> i64 {
    i64::from(node_at(aux, hap, i).sbeg)
}

/// bcftools' `node2rend`.
fn node2rend(aux: &TrAux, hap: &Hap, i: usize) -> i64 {
    let n = node_at(aux, hap, i);
    i64::from(n.sbeg) + i64::from(n.rlen)
}

/// `REF>ALT` of allele `ial`.
pub(super) fn var_string(cur: &CurRec, ial: usize) -> Vec<u8> {
    let r = cur.allele_str(0);
    let a = cur.allele_str(ial);
    let mut v = Vec::with_capacity(r.len() + a.len() + 1);
    v.extend_from_slice(r);
    v.push(b'>');
    v.extend_from_slice(a);
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splice_ref_joins_the_cds_with_padding() {
        let aux = TrAux {
            ref_: b"NNNNNNNNNNAAACCCGGGTTTNNNNNNNNNN".to_vec(),
            ..TrAux::default()
        };
        let sref = splice_ref(&aux, 0, &[(0, 3), (6, 3)]).expect("valid CDS lengths");
        assert_eq!(sref, b"NNNNNNNNNNAAAGGGTTTNNNNNNN");
    }

    #[test]
    fn splice_ref_rejects_a_wrapped_cds_length() {
        let aux = TrAux {
            ref_: vec![b'A'; 40],
            ..TrAux::default()
        };
        let err = splice_ref(&aux, 0, &[(5, 2u32.wrapping_sub(5))]).expect_err("negative length");
        assert_eq!(
            err.to_string(),
            "Error: the CDS at 6 has a negative length (-3), check the GFF"
        );
    }
}
