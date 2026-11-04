//! The per-record lookups: coding sequences (`test_cds`, `test_cds_local`), UTRs
//! (`test_utr`), splice sites (`test_splice`), whole transcripts (`test_tscript`), and
//! symbolic alleles (`test_symbolic_alt`).

use std::io::Write;

use anyhow::bail;

use super::{
    Csq, CurRec, N_REF_PAD, Output, Phase,
    hap::{HapInit, HapNode, HapType, TrAux, splice_ref},
    splice::{SPLICE_INSIDE, SPLICE_OVERLAP, Splice},
    translate::{Fill, SplicedSeq, cds_translate},
    types::{
        CSQ_CODING_SEQUENCE, CSQ_COMPOUND, CSQ_ELONGATION, CSQ_FRAMESHIFT_VARIANT,
        CSQ_INCOMPLETE_CDS, CSQ_INFRAME_DELETION, CSQ_INFRAME_INSERTION, CSQ_INTRON,
        CSQ_MISSENSE_VARIANT, CSQ_NON_CODING, CSQ_START_LOST, CSQ_STOP_GAINED, CSQ_STOP_LOST,
        CSQ_STOP_RETAINED, CSQ_SYNONYMOUS_VARIANT, CSQ_TRUNCATION, CSQ_UTR3, CSQ_UTR5, CsqT, VStr,
    },
    vbuf::{GT_VECTOR_END, gt_allele, gt_is_missing, gt_is_phased},
};
use crate::{
    commands::norm::reference::MissingSeq,
    gff::{Strand, UtrKind, is_coding},
};

/// Whether an allele is symbolic or a spanning deletion.
fn is_symbolic(a: &[u8]) -> bool {
    matches!(a.first(), Some(b'<' | b'*'))
}

impl Csq {
    /// The overlaps of `idx` (one of the GFF indexes) with `beg..=end` on the current
    /// sequence.
    fn overlaps(&mut self, which: Index, beg: u32, end: u32) -> Vec<u32> {
        let idx = match which {
            Index::Cds => &mut self.gff.idx_cds,
            Index::Utr => &mut self.gff.idx_utr,
            Index::Exon => &mut self.gff.idx_exon,
            Index::Tscript => &mut self.gff.idx_tscript,
        };
        idx.overlap(&self.chr_gff, beg, end)
            .map(|r| *r.payload)
            .collect()
    }

    /// The query range of a record: the deliberate off-by-one extension of `rlen` accounts
    /// for insertions.
    fn rec_range(cur: &CurRec) -> (u32, u32) {
        (cur.pos as u32, (cur.pos + cur.rlen) as u32)
    }

    /// Loads a transcript's reference (bcftools' `tscript_init_ref`); `false` when the
    /// sequence is missing and `--force` is on.
    fn tscript_init_ref(&mut self, tr: u32) -> anyhow::Result<Option<Vec<u8>>> {
        let (tbeg, tend) = {
            let t = &self.gff.transcripts[tr as usize];
            (t.beg, t.end)
        };
        let pad = N_REF_PAD as i64;
        let pad_beg = if i64::from(tbeg) >= pad {
            pad
        } else {
            i64::from(tbeg)
        };
        let mut seq = Vec::new();
        let missing = if self.warned.faidx_fetch_failed > 0 && self.opts.verbosity < 2 {
            MissingSeq::Silent
        } else {
            MissingSeq::Report
        };
        let found = self.fai.fetch_with(
            &self.chr_fai,
            i64::from(tbeg) - pad_beg,
            i64::from(tend) + pad,
            &mut seq,
            missing,
        )?;
        if !found {
            let chr = String::from_utf8_lossy(&self.chr_fai).into_owned();
            if !self.opts.force {
                bail!(
                    "Error: unable to fetch the region of the fasta reference {chr}:{}-{}",
                    tbeg as i32 + 1,
                    tend as i32 + 1
                );
            } else if self.opts.verbosity != 0
                && (self.warned.faidx_fetch_failed == 0 || self.opts.verbosity > 1)
            {
                eprintln!(
                    "Warning: unable to fetch the region of the fasta reference {chr}:{}-{}",
                    tbeg as i32 + 1,
                    tend as i32 + 1
                );
                if self.opts.verbosity < 2 {
                    eprintln!(
                        "         This message is printed only once, the verbosity can be increased with `--verbosity 2`"
                    );
                }
            }
            self.warned.faidx_fetch_failed += 1;
            return Ok(None);
        }
        let len = seq.len() as i64;
        let pad_end = len - (i64::from(tend) - i64::from(tbeg) + 1 + pad_beg);
        if pad_beg + pad_end != 2 * pad {
            let mut r = Vec::with_capacity((i64::from(tend) - i64::from(tbeg)) as usize + 21);
            r.resize((pad - pad_beg).max(0) as usize, b'N');
            r.extend_from_slice(&seq);
            r.resize(r.len() + (pad - pad_end).max(0) as usize, b'N');
            seq = r;
        }
        Ok(Some(seq))
    }

    /// Checks the REF allele against the transcript's reference (bcftools'
    /// `sanity_check_ref`); `false` on a mismatch with `--force`.
    fn sanity_check_ref(&mut self, cur: &CurRec, tr: u32, aux: &TrAux) -> anyhow::Result<bool> {
        let tbeg = self.gff.transcripts[tr as usize].beg;
        let mut vbeg = 0i64;
        let mut rbeg = cur.pos - i64::from(tbeg) + N_REF_PAD as i64;
        if rbeg < 0 {
            vbeg += rbeg.abs();
            rbeg = 0;
        }
        let r = aux.ref_.get(rbeg as usize..).unwrap_or(&[]);
        let v = cur.allele_str(0).get(vbeg as usize..).unwrap_or(&[]);
        for (&rc, &vc) in r.iter().zip(v) {
            if rc == 0 {
                break;
            }
            if rc != vc && !rc.eq_ignore_ascii_case(&vc) {
                let msg = format!(
                    "the fasta reference does not match the VCF REF allele at {}:{} .. fasta={} vcf={}",
                    String::from_utf8_lossy(&self.chr_vcf),
                    cur.pos + vbeg + 1,
                    char::from(rc),
                    char::from(vc)
                );
                if !self.opts.force {
                    bail!("Error: {msg}");
                } else if self.opts.verbosity != 0
                    && (self.warned.ref_allele_mismatch == 0 || self.opts.verbosity > 1)
                {
                    eprintln!("Warning: {msg}");
                    if self.opts.verbosity < 2 {
                        eprintln!(
                            "         This message is printed only once, the verbosity can be increased with `--verbosity 2`"
                        );
                    }
                }
                self.warned.ref_allele_mismatch += 1;
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Writes a `LOG` line to the text output.
    fn text_log(&mut self, msg: &str) -> anyhow::Result<()> {
        if let Output::Text(out) = &mut self.out {
            writeln!(out, "LOG\t{msg}")?;
        }
        Ok(())
    }

    /// Localized predictions, one record at a time (bcftools' `test_cds_local`).
    pub(super) fn test_cds_local(&mut self, cur: &CurRec) -> anyhow::Result<bool> {
        let (beg, end) = Self::rec_range(cur);
        let cds_ids = self.overlaps(Index::Cds, beg, end);
        if cds_ids.is_empty() {
            return Ok(false);
        }
        let mut ret = false;
        let root = HapNode {
            ty: HapType::Root,
            ..HapNode::default()
        };
        for cds_id in cds_ids {
            let tr = self.gff.cds[cds_id as usize].tr;
            let (tr_ty, tr_id) = {
                let t = &self.gff.transcripts[tr as usize];
                (t.ty, t.id)
            };
            if !is_coding(tr_ty) {
                continue;
            }
            ret = true;
            if self.aux[tr as usize].is_none() {
                let Some(ref_) = self.tscript_init_ref(tr)? else {
                    self.claim_tscript_alleles(cur, tr_id)?;
                    continue;
                };
                let mut aux = TrAux {
                    ref_,
                    ..TrAux::default()
                };
                let cds: Vec<(u32, u32)> = self.gff.transcripts[tr as usize]
                    .cds
                    .iter()
                    .map(|&c| (self.gff.cds[c as usize].beg, self.gff.cds[c as usize].len))
                    .collect();
                let tbeg = self.gff.transcripts[tr as usize].beg;
                aux.sref = Some(splice_ref(&aux, tbeg, &cds)?);
                self.aux[tr as usize] = Some(Box::new(aux));
                // only to clean the reference afterwards
                self.active_tr.insert(tr, &self.tr_ends);
            }
            let aux = self.aux[tr as usize]
                .take()
                .expect("the transcript was just initialized");
            let result = self.test_cds_local_tr(cur, cds_id, &root, &aux);
            self.aux[tr as usize] = Some(aux);
            result?;
        }
        Ok(ret)
    }

    /// [`Self::test_cds_local`] for one CDS of an initialized transcript.
    fn test_cds_local_tr(
        &mut self,
        cur: &CurRec,
        cds_id: u32,
        root: &HapNode,
        aux: &TrAux,
    ) -> anyhow::Result<()> {
        let tr = self.gff.cds[cds_id as usize].tr;
        let tr_id = self.gff.transcripts[tr as usize].id;
        if !self.sanity_check_ref(cur, tr, aux)? {
            self.claim_tscript_alleles(cur, tr_id)?;
            return Ok(());
        }
        let sref = aux.sref.as_deref().unwrap_or(&[]);
        let pad = N_REF_PAD as i64;
        for i in 1..cur.n_allele() {
            if is_symbolic(cur.allele_str(i)) {
                continue;
            }
            let mut node = HapNode::default();
            let hap_ret = self.hap_init(cur, aux, root, &mut node, cds_id, i)?;
            if hap_ret != HapInit::Silent {
                self.claim_tscript_allele(cur, tr_id, i as u32)?; // new or overlapping variant
            }
            if hap_ret != HapInit::Added {
                continue;
            }
            let t = &self.gff.transcripts[tr as usize];
            let strand = t.strand;
            let mut csq = CsqT::new(0, t, &self.gff, i as u32);
            let mut csq_type = node.csq;

            if node.ty == HapType::Sss {
                csq.ty.ty = csq_type;
                self.csq_stage(cur, &mut csq)?;
                continue;
            }
            let seq = node.seq.as_deref().unwrap_or(&[]);
            let alen = seq.len() as i32;
            let fill = if node.dlen % 3 != 0 && alen != 0 {
                Fill::ToEnd
            } else {
                Fill::No
            };
            let mut hap = std::mem::take(&mut self.hap);
            let nsref = sref.len() as i64;
            cds_translate(
                self.gencode,
                sref,
                &SplicedSeq {
                    bases: seq,
                    total: nsref - 2 * pad + i64::from(node.dlen),
                },
                node.sbeg,
                node.sbeg,
                node.sbeg.wrapping_add(node.rlen as u32),
                strand,
                &mut hap.tseq,
                &mut hap.tseq_stop,
                fill,
            )?;
            let rfrom = (pad + i64::from(node.sbeg)) as usize;
            let rseq = sref
                .get(rfrom..rfrom + node.rlen.max(0) as usize)
                .unwrap_or(&[]);
            cds_translate(
                self.gencode,
                sref,
                &SplicedSeq {
                    bases: rseq,
                    total: nsref - 2 * pad,
                },
                node.sbeg,
                node.sbeg,
                node.sbeg.wrapping_add(node.rlen as u32),
                strand,
                &mut hap.tref,
                &mut hap.tref_stop,
                fill,
            )?;
            let result = self.local_csq(
                cur,
                &mut hap,
                &node,
                csq_type,
                &mut csq,
                (nsref, alen, strand),
            );
            self.hap = hap;
            match result? {
                Some(t) => csq_type = t,
                None => continue,
            }
            if csq_type & !CSQ_COMPOUND != 0 {
                csq.ty.ty = csq_type & !CSQ_COMPOUND;
                csq.ty.vstr.shown = false;
                self.csq_stage(cur, &mut csq)?;
            }
        }
        Ok(())
    }

    /// The consequences of one localized variant from its translations; `None` when it was
    /// staged as a start loss (the second half of bcftools' `test_cds_local` loop).
    fn local_csq(
        &mut self,
        cur: &CurRec,
        hap: &mut super::hap::Hap,
        node: &HapNode,
        mut csq_type: u32,
        csq: &mut CsqT,
        (nsref, alen, strand): (i64, i32, Strand),
    ) -> anyhow::Result<Option<u32>> {
        let pad = N_REF_PAD as i64;
        let (tref, tseq) = (&mut hap.tref, &mut hap.tseq);
        let (tref_stop, tseq_stop) = (&mut hap.tref_stop, &mut hap.tseq_stop);
        // check for truncating stops
        if let Some(j) = (0..tref_stop.len()).find(|&j| tref_stop.at(j) == b'*') {
            tref.truncate(j + 1);
            tref_stop.truncate(j + 1);
        }
        let j = (0..tseq_stop.len())
            .find(|&j| tseq_stop.at(j) == b'*')
            .unwrap_or(tseq_stop.len());
        if j != tseq.len() {
            tseq.truncate(j + 1);
            tseq_stop.truncate(j + 1);
        }
        let ref_last = || tref_stop.at(tref_stop.len().wrapping_sub(1));
        let seq_last = || tseq_stop.at(tseq_stop.len().wrapping_sub(1));
        if csq_type & CSQ_STOP_LOST != 0 {
            if ref_last() == b'*' && ref_last() == seq_last() {
                csq_type &= !CSQ_STOP_LOST;
                csq_type |= CSQ_STOP_RETAINED;
            } else if ref_last() != b'*' {
                // CDS 3' incomplete
                if seq_last() == b'*' {
                    csq_type &= !CSQ_STOP_GAINED;
                    csq_type |= CSQ_STOP_RETAINED;
                } else {
                    csq_type |= CSQ_INCOMPLETE_CDS;
                }
            }
        }
        if csq_type & CSQ_START_LOST != 0 && tref_stop.at(0) != b'M' {
            csq_type &= !CSQ_START_LOST;
        }
        if node.dlen != 0 {
            if node.dlen % 3 != 0 {
                csq_type |= CSQ_FRAMESHIFT_VARIANT;
            } else if node.dlen < 0 {
                csq_type |= CSQ_INFRAME_DELETION;
            } else {
                csq_type |= CSQ_INFRAME_INSERTION;
            }
            if ref_last() != b'*' && seq_last() == b'*' {
                csq_type |= CSQ_STOP_GAINED;
            }
            if csq_type & CSQ_START_LOST != 0 && csq_type & CSQ_FRAMESHIFT_VARIANT != 0 {
                // prevents "start_lost|..|+|1M>1?|4959GA>G" next to "start_lost|..|+"
                // (bcftools stages the consequence before clearing its type here)
                self.csq_stage(cur, csq)?;
                return Ok(None);
            }
        } else {
            let mut aa_change = false;
            for j in 0..tref.len() {
                if tref.at(j) == tseq.at(j) {
                    continue;
                }
                aa_change = true;
                if tref_stop.at(j) == b'*' {
                    csq_type |= CSQ_STOP_LOST;
                } else if tseq_stop.at(j) == b'*' {
                    csq_type |= CSQ_STOP_GAINED;
                } else {
                    csq_type |= CSQ_MISSENSE_VARIANT;
                }
            }
            if !aa_change {
                csq_type |= CSQ_SYNONYMOUS_VARIANT;
            }
        }
        if csq_type & CSQ_COMPOUND != 0 {
            // the aa variant string
            // C mixes int and uint32_t here, so a negative difference wraps before the
            // division.
            let base = (nsref - 2 * pad) as u32;
            let aa_rbeg = if strand == Strand::Fwd {
                (node.sbeg / 3 + 1) as i32
            } else {
                (base.wrapping_sub(node.sbeg).wrapping_sub(node.rlen as u32) / 3 + 1) as i32
            };
            let aa_sbeg = if strand == Strand::Fwd {
                (node.sbeg / 3 + 1) as i32
            } else {
                (base
                    .wrapping_add(node.dlen as u32)
                    .wrapping_sub(node.sbeg)
                    .wrapping_sub(alen as u32)
                    / 3
                    + 1) as i32
            };
            let mut str = vec![b'|'];
            str.extend_from_slice(aa_rbeg.to_string().as_bytes());
            self.kprint_aa_prediction(aa_rbeg, &hap.tref, &hap.tref_stop, &mut str);
            if csq_type & CSQ_SYNONYMOUS_VARIANT == 0 {
                str.push(b'>');
                str.extend_from_slice(aa_sbeg.to_string().as_bytes());
                self.kprint_aa_prediction(aa_sbeg, &hap.tseq, &hap.tseq_stop, &mut str);
            }
            str.push(b'|');
            str.extend_from_slice(((cur.pos + 1) as i32).to_string().as_bytes());
            str.extend_from_slice(&node.var);
            csq.ty.vstr = VStr::new(str);
            csq.ty.ty = csq_type & CSQ_COMPOUND;
            self.csq_stage(cur, csq)?;
        }
        Ok(Some(csq_type))
    }

    /// Haplotype-aware predictions: adds the record's alleles to the haplotype trees of the
    /// coding transcripts it overlaps (bcftools' `test_cds`).
    pub(super) fn test_cds(&mut self, cur: &CurRec) -> anyhow::Result<bool> {
        let (beg, end) = Self::rec_range(cur);
        let cds_ids = self.overlaps(Index::Cds, beg, end);
        if cds_ids.is_empty() {
            return Ok(false);
        }
        let mut ret = false;
        for cds_id in cds_ids {
            let tr = self.gff.cds[cds_id as usize].tr;
            let (tr_ty, tr_id, tr_end) = {
                let t = &self.gff.transcripts[tr as usize];
                (t.ty, t.id, t.end)
            };
            if !is_coding(tr_ty) {
                continue;
            }
            if let Some(vbuf) = self.vbufs.back_mut()
                && vbuf.keep_until < tr_end
            {
                vbuf.keep_until = tr_end;
            }
            ret = true;
            if self.aux[tr as usize].is_none() {
                // initialize the transcript and its haplotype tree, fetch the reference
                let Some(ref_) = self.tscript_init_ref(tr)? else {
                    self.claim_tscript_alleles(cur, tr_id)?;
                    continue;
                };
                let nhap = if self.opts.phase == Phase::DropGt {
                    1
                } else {
                    2 * self.smpl.len()
                };
                let aux = TrAux {
                    ref_,
                    sref: None,
                    nodes: vec![HapNode {
                        ty: HapType::Root,
                        nend: nhap as u32,
                        ..HapNode::default()
                    }],
                    has_root: true,
                    hap: vec![None; nhap],
                };
                self.aux[tr as usize] = Some(Box::new(aux));
                self.active_tr.insert(tr, &self.tr_ends);
            }
            let mut aux = self.aux[tr as usize]
                .take()
                .expect("the transcript was just initialized");
            let result = self.test_cds_tr(cur, cds_id, &mut aux);
            self.aux[tr as usize] = Some(aux);
            if result? {
                ret = true;
            }
        }
        Ok(ret)
    }

    /// [`Self::test_cds`] for one CDS of an initialized transcript.
    fn test_cds_tr(&mut self, cur: &CurRec, cds_id: u32, aux: &mut TrAux) -> anyhow::Result<bool> {
        let tr = self.gff.cds[cds_id as usize].tr;
        let tr_id = self.gff.transcripts[tr as usize].id;
        if !self.sanity_check_ref(cur, tr, aux)? {
            self.claim_tscript_alleles(cur, tr_id)?;
            return Ok(false);
        }
        if self.opts.phase == Phase::DropGt {
            return self.test_cds_drop_gt(cur, cds_id, aux);
        }
        // apply the VCF variants and extend the haplotype tree
        let ngts_total = self.load_genotypes(cur)?;
        self.gt.restore();
        let ngts = ngts_total / self.nsamples.max(1) as i32;
        if ngts != 1 && ngts != 2 {
            if self.opts.verbosity != 0 && (!self.warned.multiploid || self.opts.verbosity > 1) {
                eprintln!(
                    "Warning: Skipping site with non-diploid/non-haploid genotypes at {}:{}\t{}>{}.",
                    String::from_utf8_lossy(&self.chr_vcf),
                    cur.pos + 1,
                    String::from_utf8_lossy(cur.allele_str(0)),
                    String::from_utf8_lossy(cur.allele_str(1))
                );
                if !self.warned.multiploid {
                    eprintln!(
                        "         This message is printed only once, the verbosity can be increased with `--verbosity 2`"
                    );
                }
                self.warned.multiploid = true;
            }
            let msg = format!(
                "Warning: Skipping site with non-diploid/non-haploid genotypes at {}:{}\t{}>{}",
                String::from_utf8_lossy(&self.chr_vcf),
                cur.pos + 1,
                String::from_utf8_lossy(cur.allele_str(0)),
                String::from_utf8_lossy(cur.allele_str(1))
            );
            self.text_log(&msg)?;
            self.claim_tscript_alleles(cur, tr_id)?;
            return Ok(false);
        }
        let ngts = ngts as usize;
        for ismpl in 0..self.smpl.len() {
            let base = self.smpl[ismpl] * ngts;
            if gt_is_missing(self.gt.arr[base]) {
                continue;
            }
            if ngts > 1 {
                let (g0, g1) = (self.gt.arr[base], self.gt.arr[base + 1]);
                if !gt_is_missing(g1) && g1 != GT_VECTOR_END && gt_allele(g0) != gt_allele(g1) {
                    if self.opts.phase == Phase::Merge && gt_allele(g0) == 0 {
                        self.gt.set(base, g1);
                    }
                    let (g0, g1) = (self.gt.arr[base], self.gt.arr[base + 1]);
                    if !gt_is_phased(g0) && !gt_is_phased(g1) {
                        match self.opts.phase {
                            Phase::Require => bail!(
                                "Unphased heterozygous genotype at {}:{}, sample {}. See the --phase option.",
                                String::from_utf8_lossy(&self.chr_vcf),
                                cur.pos + 1,
                                String::from_utf8_lossy(&self.hdr.samples()[self.smpl[ismpl]])
                            ),
                            Phase::Skip => {
                                if !self.opts.greedy {
                                    continue;
                                }
                                for a in [gt_allele(g0), gt_allele(g1)] {
                                    if a > 0
                                        && (a as usize) < cur.n_allele()
                                        && !is_symbolic(cur.allele_str(a as usize))
                                    {
                                        self.claim_tscript_allele(cur, tr_id, a as u32)?;
                                    }
                                }
                                continue;
                            }
                            Phase::NonRef => {
                                if gt_allele(g0) == 0 {
                                    self.gt.set(base, g1);
                                } else if gt_allele(g1) == 0 {
                                    self.gt.set(base + 1, g0);
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
            for ihap in 0..ngts {
                let g = self.gt.arr[base + ihap];
                if gt_is_missing(g) || g == GT_VECTOR_END {
                    continue;
                }
                let i = 2 * ismpl + ihap;
                let ial = gt_allele(g);
                if ial == 0 {
                    continue;
                }
                let ial = ial as usize;
                if ial >= cur.n_allele() {
                    bail!("csq: assertion failed in test_cds: ial < rec->n_allele");
                }
                if is_symbolic(cur.allele_str(ial)) {
                    continue;
                }
                self.apply_allele(cur, cds_id, aux, i, ial, Some(ismpl))?;
            }
        }
        Ok(false)
    }

    /// Adds allele `ial` to haplotype `i` of a transcript's tree (the inner loop of
    /// bcftools' `test_cds`); `ismpl` is `None` without samples (`--samples -`).
    fn apply_allele(
        &mut self,
        cur: &CurRec,
        cds_id: u32,
        aux: &mut TrAux,
        i: usize,
        ial: usize,
        ismpl: Option<usize>,
    ) -> anyhow::Result<()> {
        let tr = self.gff.cds[cds_id as usize].tr;
        let tr_id = self.gff.transcripts[tr as usize].id;
        let parent = aux.hap[i].unwrap_or(TrAux::ROOT);
        let serial = cur.serial;
        {
            let p = &aux.nodes[parent as usize];
            if p.cur_rec == Some(serial) && p.cur_child.get(ial).copied().unwrap_or(-1) >= 0 {
                // this haplotype has been seen in another sample
                self.claim_tscript_allele(cur, tr_id, ial as u32)?;
                let c = p.child[p.cur_child[ial] as usize];
                aux.hap[i] = Some(c);
                aux.nodes[c as usize].nend += 1;
                aux.nodes[parent as usize].nend -= 1;
                return Ok(());
            }
        }
        let mut child = HapNode::default();
        let hap_ret = self.hap_init(
            cur,
            aux,
            &aux.nodes[parent as usize],
            &mut child,
            cds_id,
            ial,
        )?;
        if hap_ret != HapInit::Silent {
            self.claim_tscript_allele(cur, tr_id, ial as u32)?;
        }
        if hap_ret != HapInit::Added {
            // overlapping or intron variant, cannot apply
            if hap_ret == HapInit::Overlapping {
                self.warn_overlapping(cur, ial, ismpl)?;
            }
            return Ok(());
        }
        if child.ty == HapType::Sss {
            let t = &self.gff.transcripts[tr as usize];
            let mut csq = CsqT::new(child.csq, t, &self.gff, ial as u32);
            self.csq_stage(cur, &mut csq)?;
            return Ok(());
        }
        let p = &mut aux.nodes[parent as usize];
        if p.cur_rec != Some(serial) {
            p.cur_child.clear();
            p.cur_child.resize(cur.n_allele(), -1);
            p.cur_rec = Some(serial);
        }
        let j = p.child.len();
        let c = aux.nodes.len() as u32;
        let p = &mut aux.nodes[parent as usize];
        p.child.push(c);
        p.cur_child[ial] = j as i32;
        p.nend -= 1;
        child.nend += 1;
        aux.nodes.push(child);
        aux.hap[i] = Some(c);
        Ok(())
    }

    /// The overlapping-variants warning of `test_cds`.
    fn warn_overlapping(
        &mut self,
        cur: &CurRec,
        ial: usize,
        ismpl: Option<usize>,
    ) -> anyhow::Result<()> {
        let chr = String::from_utf8_lossy(&self.chr_vcf).into_owned();
        let r = String::from_utf8_lossy(cur.allele_str(0)).into_owned();
        let a = String::from_utf8_lossy(cur.allele_str(ial)).into_owned();
        let at = match ismpl {
            Some(s) => format!(
                "{chr}:{}, sample {}",
                cur.pos + 1,
                String::from_utf8_lossy(&self.hdr.samples()[self.smpl[s]])
            ),
            None => format!("{chr}:{}", cur.pos + 1),
        };
        if self.opts.verbosity != 0 && (!self.warned.overlaps || self.opts.verbosity > 1) {
            eprintln!("Warning: Skipping overlapping variants at {at}\t{r}>{a}.");
            if !self.warned.overlaps {
                eprintln!(
                    "         This message is printed only once, the verbosity can be increased with `--verbosity 2`"
                );
            }
            self.warned.overlaps = true;
        }
        self.text_log(&format!(
            "Warning: Skipping overlapping variants at {at}\t{r}>{a}"
        ))
    }

    /// [`Self::test_cds_tr`] without genotypes: one haplotype with the first alternate
    /// allele.
    fn test_cds_drop_gt(
        &mut self,
        cur: &CurRec,
        cds_id: u32,
        aux: &mut TrAux,
    ) -> anyhow::Result<bool> {
        let tr = self.gff.cds[cds_id as usize].tr;
        let tr_id = self.gff.transcripts[tr as usize].id;
        if is_symbolic(cur.allele_str(1)) {
            return Ok(false);
        }
        let parent = aux.hap[0].unwrap_or(TrAux::ROOT);
        let mut child = HapNode::default();
        let hap_ret =
            self.hap_init(cur, aux, &aux.nodes[parent as usize], &mut child, cds_id, 1)?;
        if hap_ret != HapInit::Silent {
            self.claim_tscript_allele(cur, tr_id, 1)?;
        }
        if hap_ret != HapInit::Added {
            // overlapping or intron variant, cannot apply
            if hap_ret == HapInit::Overlapping {
                self.warn_overlapping(cur, 1, None)?;
                return Ok(false);
            }
            // no csq in the CDS path; --greedy may still add intron via test_tscript
            return Ok(true);
        }
        if child.ty == HapType::Sss {
            let t = &self.gff.transcripts[tr as usize];
            let mut csq = CsqT::new(child.csq, t, &self.gff, 1);
            self.csq_stage(cur, &mut csq)?;
            return Ok(true);
        }
        let c = aux.nodes.len() as u32;
        let p = &mut aux.nodes[parent as usize];
        p.nend -= 1;
        p.child.clear();
        p.child.push(c);
        child.nend = 1;
        aux.nodes.push(child);
        aux.hap[0] = Some(c);
        Ok(false)
    }

    /// UTR consequences (bcftools' `test_utr`).
    pub(super) fn test_utr(&mut self, cur: &CurRec) -> anyhow::Result<bool> {
        let (beg, end) = Self::rec_range(cur);
        let utrs = self.overlaps(Index::Utr, beg, end);
        if utrs.is_empty() {
            return Ok(false);
        }
        let mut ret = false;
        let first_tr = self.gff.utrs[utrs[0] as usize].tr;
        let mut splice = Splice::new(cur, &self.gff, first_tr, None);
        for u in utrs {
            let utr = self.gff.utrs[u as usize];
            splice.set_tr(&self.gff, utr.tr);
            for i in 1..cur.n_allele() {
                if is_symbolic(cur.allele_str(i)) {
                    continue;
                }
                splice.set_alt(cur, i);
                splice.csq = 0;
                let r = self.splice_csq(cur, &mut splice, utr.beg, utr.end)?;
                if r != SPLICE_INSIDE && r != SPLICE_OVERLAP {
                    continue;
                }
                let t = &self.gff.transcripts[utr.tr as usize];
                let ty = if utr.which == UtrKind::Prime5 {
                    CSQ_UTR5
                } else {
                    CSQ_UTR3
                };
                let mut csq = CsqT::new(ty, t, &self.gff, i as u32);
                let id = t.id;
                self.claim_tscript_allele(cur, id, i as u32)?;
                self.csq_stage(cur, &mut csq)?;
                ret = true;
            }
        }
        Ok(ret)
    }

    /// Splice-site consequences of exons (bcftools' `test_splice`).
    pub(super) fn test_splice(&mut self, cur: &CurRec) -> anyhow::Result<bool> {
        let (beg, end) = Self::rec_range(cur);
        let exons = self.overlaps(Index::Exon, beg, end);
        if exons.is_empty() {
            return Ok(false);
        }
        let mut ret = false;
        let first_tr = self.gff.exons[exons[0] as usize].tr;
        let mut splice = Splice::new(cur, &self.gff, first_tr, None);
        splice.check_acceptor = true;
        splice.check_donor = true;
        for e in exons {
            let exon = self.gff.exons[e as usize];
            splice.set_tr(&self.gff, exon.tr);
            let t = &self.gff.transcripts[exon.tr as usize];
            if t.cds.is_empty() {
                continue; // not a coding transcript, no interest in splice sites
            }
            splice.check_region_beg = t.beg != exon.beg;
            splice.check_region_end = t.end != exon.end;
            let id = t.id;
            for i in 1..cur.n_allele() {
                if is_symbolic(cur.allele_str(i)) {
                    continue;
                }
                splice.set_alt(cur, i);
                splice.csq = 0;
                self.splice_csq(cur, &mut splice, exon.beg, exon.end)?;
                if splice.csq != 0 {
                    self.claim_tscript_allele(cur, id, i as u32)?;
                    ret = true;
                }
            }
        }
        Ok(ret)
    }

    /// Intronic and non-coding consequences (bcftools' `test_tscript`).
    pub(super) fn test_tscript(&mut self, cur: &CurRec) -> anyhow::Result<bool> {
        let (beg, end) = Self::rec_range(cur);
        let trs = self.overlaps(Index::Tscript, beg, end);
        if trs.is_empty() {
            return Ok(false);
        }
        let mut ret = false;
        let mut splice = Splice::new(cur, &self.gff, trs[0], None);
        for tr in trs {
            splice.set_tr(&self.gff, tr);
            let (ty, id, tbeg, tend) = {
                let t = &self.gff.transcripts[tr as usize];
                (t.ty, t.id, t.beg, t.end)
            };
            for i in 1..cur.n_allele() {
                if is_symbolic(cur.allele_str(i)) {
                    continue;
                }
                splice.set_alt(cur, i);
                splice.csq = 0;
                if is_coding(ty) && self.tscript_allele_is_claimed(cur, id, i as u32)? {
                    continue;
                }
                let r = self.splice_csq(cur, &mut splice, tbeg, tend)?;
                if r != SPLICE_INSIDE && r != SPLICE_OVERLAP {
                    continue; // SPLICE_OUTSIDE or SPLICE_REF
                }
                let t = &self.gff.transcripts[tr as usize];
                let cty = if is_coding(ty) {
                    CSQ_INTRON
                } else {
                    CSQ_NON_CODING
                };
                let mut csq = CsqT::new(cty, t, &self.gff, i as u32);
                self.csq_stage(cur, &mut csq)?;
                ret = true;
            }
        }
        Ok(ret)
    }

    /// Symbolic `<INS...>` and `<DEL...>` alleles (bcftools' `test_symbolic_alt`).
    pub(super) fn test_symbolic_alt(&mut self, cur: &CurRec) -> anyhow::Result<()> {
        if self.opts.verbosity != 0 && !self.warned.symbolic && self.opts.verbosity > 0 {
            eprintln!("Warning: The support for symbolic ALT insertions is experimental.");
            self.warned.symbolic = true;
        }
        // only insertions atm
        let beg = (cur.pos + 1) as u32;
        let end = beg;
        let alt = cur.allele_str(1);
        let prefix = |p: &[u8]| alt.len() >= p.len() && alt[..p.len()].eq_ignore_ascii_case(p);
        let csq_class = if prefix(b"<INS") {
            CSQ_ELONGATION
        } else if prefix(b"<DEL") {
            CSQ_TRUNCATION
        } else {
            return Ok(());
        };
        let mut hit = false;
        for c in self.overlaps(Index::Cds, beg, end) {
            let tr = self.gff.cds[c as usize].tr;
            let t = &self.gff.transcripts[tr as usize];
            let ty = if is_coding(t.ty) {
                CSQ_CODING_SEQUENCE
            } else {
                CSQ_NON_CODING
            } | csq_class;
            let mut csq = CsqT::new(ty, t, &self.gff, 1);
            let id = t.id;
            self.claim_tscript_allele(cur, id, 1)?;
            self.csq_stage(cur, &mut csq)?;
            hit = true;
        }
        for u in self.overlaps(Index::Utr, beg, end) {
            let utr = self.gff.utrs[u as usize];
            let t = &self.gff.transcripts[utr.tr as usize];
            let ty = if utr.which == UtrKind::Prime5 {
                CSQ_UTR5
            } else {
                CSQ_UTR3
            } | csq_class;
            let mut csq = CsqT::new(ty, t, &self.gff, 1);
            let id = t.id;
            self.claim_tscript_allele(cur, id, 1)?;
            self.csq_stage(cur, &mut csq)?;
            hit = true;
        }
        let exons = self.overlaps(Index::Exon, beg, end);
        if !exons.is_empty() {
            let mut splice =
                Splice::new(cur, &self.gff, self.gff.exons[exons[0] as usize].tr, None);
            splice.check_acceptor = true;
            splice.check_donor = true;
            for e in exons {
                let exon = self.gff.exons[e as usize];
                splice.set_tr(&self.gff, exon.tr);
                let t = &self.gff.transcripts[exon.tr as usize];
                if t.cds.is_empty() {
                    continue; // not a coding transcript, no interest in splice sites
                }
                splice.check_region_beg = t.beg != exon.beg;
                splice.check_region_end = t.end != exon.end;
                let id = t.id;
                splice.set_alt(cur, 1);
                splice.csq = csq_class;
                self.splice_csq(cur, &mut splice, exon.beg, exon.end)?;
                if splice.csq != 0 {
                    self.claim_tscript_allele(cur, id, 1)?;
                    hit = true;
                }
            }
        }
        if !hit || self.opts.greedy {
            let trs = self.overlaps(Index::Tscript, beg, end);
            if !trs.is_empty() {
                let mut splice = Splice::new(cur, &self.gff, trs[0], None);
                for tr in trs {
                    splice.set_tr(&self.gff, tr);
                    splice.set_alt(cur, 1);
                    splice.csq = csq_class;
                    let (ty, id, tbeg, tend) = {
                        let t = &self.gff.transcripts[tr as usize];
                        (t.ty, t.id, t.beg, t.end)
                    };
                    if is_coding(ty) && self.tscript_allele_is_claimed(cur, id, 1)? {
                        continue;
                    }
                    let r = self.splice_csq(cur, &mut splice, tbeg, tend)?;
                    if r != SPLICE_INSIDE && r != SPLICE_OVERLAP {
                        continue;
                    }
                    let t = &self.gff.transcripts[tr as usize];
                    let cty = if is_coding(ty) {
                        CSQ_INTRON
                    } else {
                        CSQ_NON_CODING
                    } | csq_class;
                    let mut csq = CsqT::new(cty, t, &self.gff, 1);
                    self.csq_stage(cur, &mut csq)?;
                }
            }
        }
        Ok(())
    }
}

/// The GFF indexes.
#[derive(Clone, Copy, Debug)]
enum Index {
    Cds,
    Utr,
    Exon,
    Tscript,
}
