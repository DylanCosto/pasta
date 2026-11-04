//! Buffered records and their consequences (bcftools' `vbuf_t`, `vrec_t`, `vbuf_push`,
//! `vbuf_flush`, `csq_push`, `csq_stage`, and `hap_flush`): records wait until every
//! transcript they overlap is finished, then are written with INFO/BCSQ and FORMAT/BCSQ (or
//! as text lines).

use std::{
    collections::HashSet,
    io::{self, Write},
};

use anyhow::{Context, bail};
use pasta::variant::{Record, io::WriteError};

use super::{
    BCSQ_BITS_PER_INT, Csq, CurRec, Output, Phase,
    hap::HapNode,
    types::{
        CSQ_COMPOUND, CSQ_ELONGATION, CSQ_INFRAME_DELETION, CSQ_INFRAME_INSERTION,
        CSQ_PRINTED_UPSTREAM, CSQ_PRN_TSCRIPT, CSQ_SPLICE_ACCEPTOR, CSQ_SPLICE_DONOR,
        CSQ_SPLICE_REGION, CSQ_START_LOST, CSQ_START_RETAINED, CSQ_START_STOP, CSQ_STOP_LOST,
        CSQ_STOP_RETAINED, CSQ_SYNONYMOUS_VARIANT, CSQ_TRUNCATION, CSQ_UPSTREAM_STOP, CsqT, RecRef,
        VCsq, VrecHandle, kput_vcsq,
    },
};

/// A delayed text-output consequence (bcftools' `txt_csq_t`).
#[derive(Clone, Copy, Debug)]
pub(super) struct TxtCsq {
    idx: i32,
    /// The header sample index, or -1.
    ismpl: i32,
    /// 1 or 2 for a haplotype, 0 for none.
    ihap: i32,
}

/// A buffered record with its consequences (bcftools' `vrec_t`).
pub(super) struct VRec {
    pub line: Record,
    pub serial: u64,
    pub pos: i64,
    /// Sample bitmasks, `nfmt_bcsq` values per sample.
    fmt_bm: Vec<u32>,
    nfmt: usize,
    pub vcsq: Vec<VCsq>,
    txt: Vec<TxtCsq>,
    /// The entries of `txt`, for finding duplicates.
    txt_seen: HashSet<(i32, i32, i32)>,
    /// Transcript/allele pairs handled by the CDS/UTR/splice logic (`--greedy`).
    claimed: Vec<(u32, u32)>,
}

/// The records at one position (bcftools' `vbuf_t`).
pub(super) struct VBuf {
    id: u64,
    vrec: Vec<VRec>,
    n: usize,
    /// The largest end of the coding transcripts the records overlap.
    pub keep_until: u32,
}

/// The heap of active transcripts by end (bcftools' `kheap` with `cmp_tscript`).
#[derive(Debug, Default)]
pub(super) struct TrHeap {
    pub dat: Vec<u32>,
}

impl TrHeap {
    fn less(ends: &[u32], a: u32, b: u32) -> bool {
        ends[a as usize] < ends[b as usize]
    }

    /// `khp_insert`.
    pub(super) fn insert(&mut self, tr: u32, ends: &[u32]) {
        self.dat.push(tr);
        let mut i = self.dat.len() - 1;
        while i > 0 && Self::less(ends, tr, self.dat[(i - 1) / 2]) {
            self.dat[i] = self.dat[(i - 1) / 2];
            i = (i - 1) / 2;
        }
        self.dat[i] = tr;
    }

    /// `khp_delete`.
    pub(super) fn delete(&mut self, ends: &[u32]) {
        let Some(last) = self.dat.pop() else { return };
        if self.dat.is_empty() {
            return;
        }
        self.dat[0] = last;
        let mut i = 0;
        loop {
            let (l, r) = (2 * i + 1, 2 * i + 2);
            let mut extreme = if l < self.dat.len() && Self::less(ends, self.dat[l], self.dat[i]) {
                l
            } else {
                i
            };
            if r < self.dat.len() && Self::less(ends, self.dat[r], self.dat[extreme]) {
                extreme = r;
            }
            if extreme == i {
                break;
            }
            self.dat.swap(i, extreme);
            i = extreme;
        }
    }
}

impl Csq {
    /// The buffer with ID `id`.
    fn vbuf_index(&self, id: u64) -> Option<usize> {
        let front = self.vbufs.front()?.id;
        let i = id.checked_sub(front)? as usize;
        (i < self.vbufs.len()).then_some(i)
    }

    /// The chromosome of the buffered records (bcftools names the record's; `chr_vcf` is
    /// already the next chromosome while the previous one is flushed).
    fn buffered_chr(&self) -> std::borrow::Cow<'_, str> {
        String::from_utf8_lossy(self.hdr.contig_name(self.rid).unwrap_or(b"(unknown)"))
    }

    /// The record `rec` (bcftools' `rec2vrec`).
    pub(super) fn rec2vrec(&self, rec: RecRef) -> anyhow::Result<VrecHandle> {
        let Some(&id) = self.pos2vbuf.get(&(rec.pos as i32)) else {
            bail!(
                "This should not happen. {}:{}",
                self.buffered_chr(),
                rec.pos + 1
            );
        };
        if let Some(bi) = self.vbuf_index(id) {
            let vbuf = &self.vbufs[bi];
            if let Some(slot) = vbuf.vrec[..vbuf.n]
                .iter()
                .position(|v| v.serial == rec.serial)
            {
                return Ok(VrecHandle {
                    vbuf: id,
                    slot: slot as u32,
                });
            }
        }
        bail!(
            "This should not happen.. {}:{}",
            self.buffered_chr(),
            rec.pos + 1
        );
    }

    fn vrec_mut(&mut self, h: VrecHandle) -> Option<&mut VRec> {
        let bi = self.vbuf_index(h.vbuf)?;
        self.vbufs[bi].vrec.get_mut(h.slot as usize)
    }

    /// The record that was just pushed.
    pub(super) fn cur_vrec(&mut self, cur: &CurRec) -> anyhow::Result<&mut VRec> {
        let h = self.rec2vrec(cur.rec_ref())?;
        self.vrec_mut(h)
            .ok_or_else(|| anyhow::anyhow!("csq: the current record is not buffered"))
    }

    /// Records that a transcript/allele pair was handled (bcftools'
    /// `claim_tscript_allele`).
    pub(super) fn claim_tscript_allele(
        &mut self,
        cur: &CurRec,
        trid: u32,
        ial: u32,
    ) -> anyhow::Result<()> {
        if !self.opts.greedy {
            return Ok(());
        }
        let vrec = self.cur_vrec(cur)?;
        if !vrec.claimed.contains(&(trid, ial)) {
            vrec.claimed.push((trid, ial));
        }
        Ok(())
    }

    /// Claims all real alternate alleles (bcftools' `claim_tscript_alleles`).
    pub(super) fn claim_tscript_alleles(&mut self, cur: &CurRec, trid: u32) -> anyhow::Result<()> {
        if !self.opts.greedy {
            return Ok(());
        }
        for ial in 1..cur.n_allele() {
            let a = cur.allele_str(ial);
            if a.first() == Some(&b'<') || a.first() == Some(&b'*') {
                continue;
            }
            self.claim_tscript_allele(cur, trid, ial as u32)?;
        }
        Ok(())
    }

    /// Whether a transcript/allele pair was handled (bcftools' `tscript_allele_is_claimed`).
    pub(super) fn tscript_allele_is_claimed(
        &mut self,
        cur: &CurRec,
        trid: u32,
        ial: u32,
    ) -> anyhow::Result<bool> {
        if !self.opts.greedy {
            return Ok(false);
        }
        Ok(self.cur_vrec(cur)?.claimed.contains(&(trid, ial)))
    }

    /// Adds a consequence to the record `rec` unless it is there already; returns whether it
    /// was a duplicate (bcftools' `csq_push`).
    pub(super) fn csq_push(&mut self, csq: &mut CsqT, rec: RecRef) -> anyhow::Result<bool> {
        let h = self.rec2vrec(rec)?;
        let vrec = self
            .vrec_mut(h)
            .ok_or_else(|| anyhow::anyhow!("csq: the record is not buffered"))?;
        csq.vrec = Some(h);
        let t = &mut csq.ty;
        if t.ty & CSQ_INFRAME_INSERTION != 0 && t.ty & CSQ_ELONGATION != 0 {
            t.ty &= !CSQ_INFRAME_INSERTION;
        }
        if t.ty & CSQ_INFRAME_DELETION != 0 && t.ty & CSQ_TRUNCATION != 0 {
            t.ty &= !CSQ_INFRAME_DELETION;
        }
        // if the variant overlaps donor/acceptor and also splice region, report only
        // donor/acceptor
        if t.ty & CSQ_SPLICE_REGION != 0 && t.ty & (CSQ_SPLICE_DONOR | CSQ_SPLICE_ACCEPTOR) != 0 {
            t.ty &= !CSQ_SPLICE_REGION;
        }
        let dup = if t.ty & CSQ_PRINTED_UPSTREAM != 0 {
            let mut found = None;
            for (i, v) in vrec.vcsq.iter_mut().enumerate() {
                // avoid both "@3632,stop_lost|..." and "stop_lost|...|3630T>A+3632A>C"
                if t.ty & CSQ_START_STOP != 0 && v.ty & CSQ_START_STOP != 0 {
                    *v = t.clone();
                    found = Some(i);
                    break;
                }
                if v.ty & CSQ_PRINTED_UPSTREAM == 0 {
                    continue;
                }
                if t.ref_ != v.ref_ {
                    continue;
                }
                found = Some(i);
                break;
            }
            found
        } else if t.ty & CSQ_COMPOUND != 0 {
            push_compound(t, &mut vrec.vcsq)
        } else {
            let mut found = None;
            for (i, v) in vrec.vcsq.iter_mut().enumerate() {
                if t.trid != v.trid && (t.ty | v.ty) & CSQ_PRN_TSCRIPT != 0 {
                    continue;
                }
                if t.biotype != v.biotype {
                    continue;
                }
                if v.ty & CSQ_COMPOUND == 0 {
                    v.ty |= t.ty;
                    found = Some(i);
                    break;
                }
                if v.ty == (v.ty | t.ty) {
                    found = Some(i);
                    break;
                }
            }
            found
        };
        match dup {
            Some(i) => {
                csq.idx = i as i32;
                Ok(true)
            }
            None => {
                csq.idx = vrec.vcsq.len() as i32;
                vrec.vcsq.push(csq.ty.clone());
                Ok(false)
            }
        }
    }

    /// The genotypes of the current record, decoded once (bcftools calls
    /// `bcf_get_genotypes`); returns its count, negative without GT.
    pub(super) fn load_genotypes(&mut self, cur: &CurRec) -> anyhow::Result<i32> {
        if let Some(n) = self.gt.ret {
            return Ok(n);
        }
        let h = self.rec2vrec(cur.rec_ref())?;
        let bi = self
            .vbuf_index(h.vbuf)
            .ok_or_else(|| anyhow::anyhow!("csq: the current record is not buffered"))?;
        let line = &mut self.vbufs[bi].vrec[h.slot as usize].line;
        let n = match line.get_genotypes(&self.hdr, &mut self.gt.orig) {
            Ok(n) => n as i32,
            Err(_) => -1,
        };
        self.gt.ret = Some(n);
        self.gt.arr.clone_from(&self.gt.orig);
        self.gt.dirty = false;
        self.gt.carriers_built = false;
        Ok(n)
    }

    /// Lists, per allele, the samples and haplotypes that carry it, in bcftools' loop order
    /// (samples in list order, then haplotypes).
    fn build_carriers(&mut self, ngt: usize) {
        let gt = &mut self.gt;
        for c in &mut gt.carriers {
            c.clear();
        }
        for &ismpl in &self.smpl {
            for j in 0..ngt {
                let g = gt.orig[ismpl * ngt + j];
                if gt_is_missing(g) || g == GT_VECTOR_END {
                    continue;
                }
                let ial = gt_allele(g);
                if ial <= 0 {
                    continue;
                }
                let ial = ial as usize;
                if gt.carriers.len() <= ial {
                    gt.carriers.resize_with(ial + 1, Vec::new);
                }
                gt.carriers[ial].push((ismpl as u32, j as u32));
            }
        }
        gt.carriers_built = true;
    }

    /// Adds a consequence of the current record and sets the sample bits of the samples
    /// that carry its allele (bcftools' `csq_stage`).
    pub(super) fn csq_stage(&mut self, cur: &CurRec, csq: &mut CsqT) -> anyhow::Result<()> {
        // Text output is delayed so it is printed from the same consequences as VCF output.
        if self.csq_push(csq, cur.rec_ref())? && self.opts.phase == Phase::DropGt {
            return Ok(()); // the consequence already exists
        }
        let mut ngt = 0;
        if self.opts.phase != Phase::DropGt {
            ngt = self.load_genotypes(cur)?;
            // bcftools re-reads the genotypes here, undoing test_cds' changes to them.
            self.gt.restore();
            if ngt > 0 {
                ngt /= self.nsamples as i32;
            }
        }
        if ngt <= 0 {
            self.text_stage(csq, -1, 0);
            return Ok(());
        }
        // bcftools asserts ngt<=2 here, but release builds go on with all ploidies.
        let ngt = ngt as usize;
        let Some(h) = csq.vrec else { return Ok(()) };
        if !self.gt.carriers_built {
            self.build_carriers(ngt);
        }
        let ial = csq.ty.vcf_ial as usize;
        let carriers = match self.gt.carriers.get_mut(ial) {
            Some(c) => std::mem::take(c),
            None => return Ok(()),
        };
        let mut skip = None;
        for &(ismpl, j) in &carriers {
            if skip == Some(ismpl) {
                continue;
            }
            let ismpl_u = ismpl as usize;
            if self.opts.text {
                self.text_stage(csq, ismpl as i32, j as i32 + 1);
                continue;
            }
            let icsq2 = 2 * csq.idx + j as i32;
            if icsq2 >= self.opts.ncsq2_max {
                // more than ncsq2_max consequences, so can't fit it in FMT
                self.warn_too_many(ismpl_u, cur.pos, icsq2 + 1, icsq2);
                skip = Some(ismpl);
                continue;
            }
            let nfmt_bcsq = self.opts.nfmt_bcsq;
            let (ival, ibit) = (
                (icsq2 / BCSQ_BITS_PER_INT) as usize,
                icsq2 % BCSQ_BITS_PER_INT,
            );
            if let Some(vrec) = self.vrec_mut(h) {
                if vrec.nfmt < 1 + ival {
                    vrec.nfmt = 1 + ival;
                }
                vrec.fmt_bm[ismpl_u * nfmt_bcsq + ival] |= 1u32 << ibit;
            }
        }
        self.gt.carriers[ial] = carriers;
        Ok(())
    }

    /// The warning of `csq_stage` for a sample with more consequences than `--ncsq`.
    fn warn_too_many(&mut self, ismpl: usize, pos: i64, keep: i32, icsq2: i32) {
        if self.opts.verbosity != 0 && (self.ncsq2_small_warned == 0 || self.opts.verbosity > 1) {
            eprintln!(
                "Warning: Too many consequences for sample {} at {}:{}, keeping the first {keep} and skipping the rest.",
                String::from_utf8_lossy(&self.hdr.samples()[ismpl]),
                String::from_utf8_lossy(&self.chr_vcf),
                pos + 1
            );
            if self.ncsq2_small_warned == 0 {
                eprintln!(
                    "         The limit can be increased by setting the --ncsq parameter. This warning is printed only once."
                );
            }
            self.ncsq2_small_warned = 1;
        }
        if self.ncsq2_small_warned < icsq2 {
            self.ncsq2_small_warned = icsq2;
        }
    }

    /// Queues a consequence for text output (bcftools' `text_stage`).
    fn text_stage(&mut self, csq: &CsqT, ismpl: i32, ihap: i32) {
        if !self.opts.text {
            return;
        }
        let Some(h) = csq.vrec else { return };
        let Some(vrec) = self.vrec_mut(h) else { return };
        if csq.idx < 0 || csq.idx as usize >= vrec.vcsq.len() {
            return;
        }
        if !vrec.txt_seen.insert((csq.idx, ismpl, ihap)) {
            return;
        }
        vrec.txt.push(TxtCsq {
            idx: csq.idx,
            ismpl,
            ihap,
        });
    }

    /// Queues a haplotype's consequences for text output (bcftools' `hap_stage_text`).
    fn hap_stage_text(&mut self, node: &HapNode, ismpl: i32, ihap: i32) {
        for csq in &node.csq_list {
            self.text_stage(csq, ismpl, ihap);
        }
    }

    /// Sets the sample bits of a haplotype's consequences (bcftools' `hap_stage_vcf`).
    fn hap_stage_vcf(&mut self, ismpl: usize, ihap: usize, node: &HapNode) {
        for csq in &node.csq_list {
            let Some(h) = csq.vrec else { continue };
            let icsq2 = 2 * csq.idx + ihap as i32;
            if icsq2 >= self.opts.ncsq2_max {
                // more than ncsq2_max consequences, so can't fit it in FMT
                if self.opts.verbosity != 0
                    && (self.ncsq2_small_warned == 0 || self.opts.verbosity > 1)
                {
                    let pos = self
                        .vbuf_index(h.vbuf)
                        .map_or(0, |bi| self.vbufs[bi].vrec[h.slot as usize].pos);
                    eprintln!(
                        "Warning: Too many consequences for sample {} at {}:{}, keeping the first {} and skipping the rest.",
                        String::from_utf8_lossy(&self.hdr.samples()[ismpl]),
                        String::from_utf8_lossy(self.hdr.contig_name(self.rid).unwrap_or(b"")),
                        pos + 1,
                        csq.idx
                    );
                    if self.ncsq2_small_warned == 0 {
                        eprintln!(
                            "         The limit can be increased by setting the --ncsq parameter. This warning is printed only once."
                        );
                    }
                }
                if self.ncsq2_small_warned < icsq2 {
                    self.ncsq2_small_warned = icsq2;
                }
                break;
            }
            let nfmt_bcsq = self.opts.nfmt_bcsq;
            let (ival, ibit) = (
                (icsq2 / BCSQ_BITS_PER_INT) as usize,
                icsq2 % BCSQ_BITS_PER_INT,
            );
            if let Some(vrec) = self.vrec_mut(h) {
                if vrec.nfmt < 1 + ival {
                    vrec.nfmt = 1 + ival;
                }
                vrec.fmt_bm[ismpl * nfmt_bcsq + ival] |= 1u32 << ibit;
            }
        }
    }

    /// Finishes the transcripts that end at or before `pos` (bcftools' `hap_flush`).
    pub(super) fn hap_flush(&mut self, pos: u32) -> anyhow::Result<()> {
        while let Some(&tr) = self.active_tr.dat.first() {
            if self.gff.transcripts[tr as usize].end > pos {
                break;
            }
            self.active_tr.delete(&self.tr_ends);
            if let Some(mut aux) = self.aux[tr as usize].take() {
                if aux.has_root && !aux.nodes[0].child.is_empty() {
                    // normal, non-localized calling
                    let result = self.hap_finalize(tr, &mut aux);
                    if result.is_ok() {
                        self.stage_haplotypes(&aux);
                    }
                    self.aux[tr as usize] = Some(aux);
                    result?;
                } else {
                    self.aux[tr as usize] = Some(aux);
                }
            }
            // Delete later: the records written by position still need it.
            self.rm_tr.push(tr);
        }
        Ok(())
    }

    /// Stages the consequences of each sample's haplotypes of a finished transcript.
    fn stage_haplotypes(&mut self, aux: &super::hap::TrAux) {
        let leaf = |i: usize| {
            aux.hap
                .get(i)
                .copied()
                .flatten()
                .map(|n| &aux.nodes[n as usize])
        };
        if self.opts.text {
            if self.opts.phase == Phase::DropGt {
                if let Some(n) = leaf(0) {
                    self.hap_stage_text(n, -1, 0);
                }
            } else {
                for i in 0..self.smpl.len() {
                    for j in 0..2 {
                        if let Some(n) = leaf(i * 2 + j) {
                            self.hap_stage_text(n, self.smpl[i] as i32, j as i32 + 1);
                        }
                    }
                }
            }
        } else if self.opts.phase != Phase::DropGt {
            for i in 0..self.smpl.len() {
                for j in 0..2 {
                    if let Some(n) = leaf(i * 2 + j) {
                        self.hap_stage_vcf(self.smpl[i], j, n);
                    }
                }
            }
        }
    }

    /// Buffers a record (bcftools' `vbuf_push`): swaps `rec` with a recycled record and
    /// returns the record's serial number.
    pub(super) fn vbuf_push(&mut self, rec: &mut Record) -> u64 {
        let pos = rec.pos();
        let new_pos = match self.vbufs.back() {
            None => true,
            Some(last) => last.n == 0 || last.vrec[0].pos != pos,
        };
        if new_pos {
            let mut vbuf = self.vbuf_pool.pop().unwrap_or_else(|| VBuf {
                id: 0,
                vrec: Vec::new(),
                n: 0,
                keep_until: 0,
            });
            vbuf.id = self.next_vbuf_id;
            self.next_vbuf_id += 1;
            vbuf.n = 0;
            vbuf.keep_until = 0;
            self.vbufs.push_back(vbuf);
        }
        let with_fmt = !self.opts.text && self.opts.phase != Phase::DropGt && !self.smpl.is_empty();
        let nfmt_len = self.hdr_nsmpl * self.opts.nfmt_bcsq;
        self.serial += 1;
        let serial = self.serial;
        let vbuf = self.vbufs.back_mut().expect("a buffer was just added");
        vbuf.n += 1;
        if vbuf.vrec.len() < vbuf.n {
            vbuf.vrec.push(VRec {
                line: Record::new(),
                serial: 0,
                pos: 0,
                fmt_bm: Vec::new(),
                nfmt: 0,
                vcsq: Vec::new(),
                txt: Vec::new(),
                txt_seen: HashSet::new(),
                claimed: Vec::new(),
            });
        }
        let vrec = &mut vbuf.vrec[vbuf.n - 1];
        vrec.vcsq.clear();
        vrec.txt.clear();
        vrec.txt_seen.clear();
        vrec.nfmt = 0;
        vrec.claimed.clear();
        if with_fmt {
            if vrec.fmt_bm.len() == nfmt_len {
                vrec.fmt_bm.fill(0);
            } else {
                vrec.fmt_bm.clear();
                vrec.fmt_bm.resize(nfmt_len, 0);
            }
        }
        vrec.serial = serial;
        vrec.pos = pos;
        std::mem::swap(rec, &mut vrec.line);
        let id = vbuf.id;
        self.pos2vbuf.insert(pos as i32, id);
        serial
    }

    /// Writes the buffered records before `pos` (bcftools' `vbuf_flush`).
    pub(super) fn vbuf_flush(&mut self, pos: u32) -> anyhow::Result<()> {
        while let Some(front) = self.vbufs.front() {
            if !self.opts.local_csq && !self.active_tr.dat.is_empty() && front.keep_until > pos {
                // cannot output buffered VCF lines until the active transcripts are gone
                break;
            }
            let mut vbuf = self.vbufs.pop_front().expect("the front buffer exists");
            let vbuf_pos = if vbuf.n > 0 { vbuf.vrec[0].pos } else { -1 };
            for i in 0..vbuf.n {
                self.write_vrec(&mut vbuf.vrec[i])?;
            }
            if vbuf_pos != -1 {
                self.pos2vbuf.remove(&(vbuf_pos as i32));
            }
            vbuf.n = 0;
            self.vbuf_pool.push(vbuf);
        }
        if !self.active_tr.dat.is_empty() {
            return Ok(());
        }
        for tr in std::mem::take(&mut self.rm_tr) {
            self.aux[tr as usize] = None;
        }
        Ok(())
    }

    /// Writes one buffered record with its consequences.
    fn write_vrec(&mut self, vrec: &mut VRec) -> anyhow::Result<()> {
        self.str_buf.clear();
        match &mut self.out {
            Output::Text(out) => {
                // tab-delimited text output
                for t in &vrec.txt {
                    if t.idx < 0 || t.idx as usize >= vrec.vcsq.len() {
                        continue;
                    }
                    let vcsq = &vrec.vcsq[t.idx as usize];
                    // CSQ_PRINTED_UPSTREAM is VCF-only bookkeeping.
                    if vcsq.ty & CSQ_PRINTED_UPSTREAM != 0 {
                        continue;
                    }
                    let b = &mut self.str_buf;
                    b.clear();
                    b.extend_from_slice(b"CSQ\t");
                    if t.ismpl >= 0 {
                        b.extend_from_slice(&self.hdr.samples()[t.ismpl as usize]);
                    } else {
                        b.push(b'-');
                    }
                    b.push(b'\t');
                    if t.ihap > 0 {
                        b.extend_from_slice(t.ihap.to_string().as_bytes());
                    } else {
                        b.push(b'-');
                    }
                    b.push(b'\t');
                    b.extend_from_slice(vrec.line.chrom_safe(&self.hdr));
                    b.push(b'\t');
                    b.extend_from_slice(((vrec.pos + 1) as i32).to_string().as_bytes());
                    b.push(b'\t');
                    kput_vcsq(&self.gff, vcsq, b);
                    b.push(b'\n');
                    out.write_all(b)?;
                }
                vrec.txt.clear();
                vrec.txt_seen.clear();
                vrec.vcsq.clear();
                vrec.nfmt = 0;
            }
            Output::Vcf { writer, name } => {
                if !vrec.vcsq.is_empty() {
                    for (j, v) in vrec.vcsq.iter().enumerate() {
                        if j > 0 {
                            self.str_buf.push(b',');
                        }
                        kput_vcsq(&self.gff, v, &mut self.str_buf);
                    }
                    vrec.line
                        .update_info_string(&self.hdr, &self.opts.bcsq_tag, Some(&self.str_buf))
                        .context("Failed to update INFO/BCSQ")?;
                    if self.hdr_nsmpl > 0 {
                        let nfmt_bcsq = self.opts.nfmt_bcsq;
                        let nfmt = vrec.nfmt;
                        if nfmt > 0 && nfmt < nfmt_bcsq {
                            for j in 1..self.hdr_nsmpl {
                                vrec.fmt_bm
                                    .copy_within(j * nfmt_bcsq..j * nfmt_bcsq + nfmt, j * nfmt);
                            }
                        }
                        let n = self.hdr_nsmpl * nfmt;
                        self.fmt_i32.clear();
                        self.fmt_i32
                            .extend(vrec.fmt_bm[..n].iter().map(|&v| v as i32));
                        vrec.line
                            .update_format_i32(&self.hdr, &self.opts.bcsq_tag, &self.fmt_i32)
                            .context("Failed to update FORMAT/BCSQ")?;
                    }
                    vrec.vcsq.clear();
                }
                if let Err(e) = writer.write_record(&self.hdr, &mut vrec.line) {
                    // A closed downstream pipe ends the run with SIGPIPE (see main).
                    if let WriteError::Io(io) = &e
                        && io.kind() == io::ErrorKind::BrokenPipe
                    {
                        return Err(anyhow::Error::new(e));
                    }
                    bail!("[vbuf_flush] Error: cannot write to {name}");
                }
            }
        }
        Ok(())
    }
}

/// Merges a haplotype-aware consequence into an equivalent one of the record (the
/// `CSQ_COMPOUND` branch of bcftools' `csq_push`); returns its index when it was merged.
fn push_compound(t: &mut VCsq, vcsq: &mut [VCsq]) -> Option<usize> {
    for (i, v) in vcsq.iter_mut().enumerate() {
        if t.trid != v.trid && (t.ty | v.ty) & CSQ_PRN_TSCRIPT != 0 {
            continue;
        }
        if t.biotype != v.biotype || t.gene != v.gene || t.vcf_ial != v.vcf_ial {
            continue;
        }
        if (t.ty & CSQ_UPSTREAM_STOP) ^ (v.ty & CSQ_UPSTREAM_STOP) != 0 {
            continue; // both must or mustn't have upstream_stop
        }
        if t.vstr.s.is_some() || v.vstr.s.is_some() {
            // splice_csq() can trigger start/stop events from an indel overlap and add_csq()
            // or test_cds_local() another from the sequence: keep one.
            let (Some(ts), Some(vs)) = (&t.vstr.s, &v.vstr.s) else {
                if t.ty & CSQ_START_STOP != 0 && v.ty & CSQ_START_STOP != 0 {
                    v.ty |= t.ty;
                    // remove stop_lost&synonymous if stop_retained set
                    if v.ty & CSQ_STOP_RETAINED != 0 {
                        v.ty &= !(CSQ_STOP_LOST | CSQ_SYNONYMOUS_VARIANT);
                    }
                    if v.ty & CSQ_START_RETAINED != 0 {
                        v.ty &= !(CSQ_START_LOST | CSQ_SYNONYMOUS_VARIANT);
                    }
                    if v.vstr.s.is_none() {
                        v.vstr = t.vstr.clone();
                    }
                    return Some(i);
                }
                continue;
            };
            if c_str(ts) != c_str(vs) {
                continue;
            }
        }
        v.ty |= t.ty;
        return Some(i);
    }
    None
}

/// A string up to its first NUL (`strcmp` compares that much).
fn c_str(s: &[u8]) -> &[u8] {
    &s[..memchr::memchr(0, s).unwrap_or(s.len())]
}

/// HTSlib's `bcf_int32_vector_end`.
pub(super) const GT_VECTOR_END: i32 = i32::MIN + 1;

/// HTSlib's `bcf_gt_is_missing`.
#[inline]
pub(super) fn gt_is_missing(g: i32) -> bool {
    (g >> 1) == 0
}

/// HTSlib's `bcf_gt_allele`.
#[inline]
pub(super) fn gt_allele(g: i32) -> i32 {
    (g >> 1) - 1
}

/// HTSlib's `bcf_gt_is_phased`.
#[inline]
pub(super) fn gt_is_phased(g: i32) -> bool {
    g & 1 != 0
}
