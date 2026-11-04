//! The `--indels-2.0` indel model: a port of bcftools' `bam2bcf_iaux.c`, which realigns
//! reads against up to two consensus haplotypes per sample with each indel type applied.

use std::sync::atomic::{AtomicBool, Ordering};

use super::{
    call::{CallAux, Plp, SEQ_NT16_INT, c_int},
    consensus::{CigarState, ReadCns, Trim},
    indel::{self, MAX_TYPES, MINUS_CONST, ProbalnBuf, RefSeq, Site},
};

/// A position whose realignment window would read outside a read (bcftools reads out of
/// bounds there): the position is skipped with this message.
#[derive(Debug)]
struct SkipSite(String);

impl std::fmt::Display for SkipSite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for SkipSite {}

static WARNED: AtomicBool = AtomicBool::new(false);

/// `indel_aux_t`
pub struct Iaux {
    iref_type: usize,
    types: Vec<i32>,
    max_ins_len: usize,
    left: i64,
    right: i64,
    inscns: Vec<u8>,
    rcns: ReadCns,
    ref_seq: Vec<u8>,
    qry_seq: Vec<u8>,
    qual: Vec<u8>,
    read_scores: Vec<i32>,
    ref_qual: [i32; MAX_TYPES],
    sum_qual: [i32; MAX_TYPES],
    probaln: ProbalnBuf,
    /// How many times consensus sequences were made.
    calls: u64,
}

impl Default for Iaux {
    fn default() -> Self {
        Self {
            iref_type: 0,
            types: Vec::new(),
            max_ins_len: 0,
            left: 0,
            right: 0,
            inscns: Vec::new(),
            rcns: ReadCns::default(),
            ref_seq: Vec::new(),
            qry_seq: Vec::new(),
            qual: Vec::new(),
            read_scores: Vec::new(),
            ref_qual: [0; MAX_TYPES],
            sum_qual: [0; MAX_TYPES],
            probaln: ProbalnBuf::default(),
            calls: 0,
        }
    }
}

impl Iaux {
    /// How many times consensus sequences were made, which reads and sets the candidate
    /// count bcftools carries between calls.
    pub fn calls(&self) -> u64 {
        self.calls
    }

    /// The candidate count the next consensus starts from.
    pub fn ncvar(&self) -> usize {
        self.rcns.ncvar()
    }

    /// Sets the candidate count the next consensus starts from.
    pub fn set_ncvar(&mut self, ncvar: usize) {
        self.rcns.set_ncvar(ncvar);
    }

    /// `iaux_init_types`: the indel types to test; 0 or less to skip the position.
    fn init_types(&mut self, site: &Site<'_, '_>, bca: &mut CallAux, chr: &[u8]) -> i32 {
        if !indel::has_indel(site.plp) {
            return 0;
        }
        bca.max_support = 0;
        self.sum_qual = [0; MAX_TYPES];
        let mut aux: Vec<u32> = vec![MINUS_CONST];
        let (mut indel_support_ok, mut n_alt, mut n_tot) = (false, 0i32, 0i32);
        let mut max_rd_len = 0i32;
        for s in site.plp.iter() {
            let nalt0 = aux.len();
            let mut ntot = 0;
            for p in s {
                ntot += 1;
                if p.indel != 0 {
                    aux.push(MINUS_CONST.wrapping_add(p.indel as u32));
                }
                let qlen = indel::cigar2qlen(p.b.cigar());
                if qlen > max_rd_len {
                    max_rd_len = qlen;
                }
            }
            let nalt = (aux.len() - nalt0) as i32;
            if bca.per_sample_flt {
                let frac = f64::from(nalt) / f64::from(ntot);
                if nalt as u32 >= bca.min_support && frac >= bca.min_frac {
                    indel_support_ok = true;
                }
                if nalt as u32 > bca.max_support && frac > 0.0 {
                    bca.max_support = nalt as u32;
                    bca.max_frac = frac as f32;
                }
            } else {
                n_alt += nalt;
                n_tot += ntot;
            }
        }
        if !bca.per_sample_flt
            && n_alt as u32 >= bca.min_support
            && f64::from(n_alt) / f64::from(n_tot) >= bca.min_frac
        {
            indel_support_ok = true;
        }
        if aux.len() == 1 || !indel_support_ok {
            return 0;
        }
        let pos = site.pos;
        let r = site.reference;
        let win = bca.indel_win_size;
        let i_end = pos + i64::from(if win < max_rd_len { win } else { max_rd_len });
        let mut n_n = 0i64;
        let mut i = pos;
        while i < i_end && r.at(i) != 0 {
            if r.at(i) == b'N' {
                n_n += 1;
            }
            i += 1;
        }
        if 2 * n_n > i - pos {
            return -1;
        }
        aux.sort_unstable();
        let mut n_types = 1;
        for i in 1..aux.len() {
            if aux[i] != aux[i - 1] {
                n_types += 1;
            }
        }
        if n_types >= MAX_TYPES {
            if !WARNED.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "Warning: excessive number of INDEL alleles at {}:{}, skipping. (This warning is printed only once)",
                    String::from_utf8_lossy(chr),
                    pos + 1
                );
            }
            return -1;
        }
        self.types.clear();
        self.max_ins_len = 0;
        let mut i = 0;
        while i < aux.len() {
            let isize = aux[i].wrapping_sub(MINUS_CONST) as i32;
            let mut j = i + 1;
            while j < aux.len() && aux[j] == aux[i] {
                j += 1;
            }
            let is_ok = if isize == 0 {
                self.iref_type = self.types.len();
                true
            } else {
                (j - i) as u32 >= bca.min_support
            };
            if is_ok {
                self.types.push(isize);
                if isize > 0 && isize as usize > self.max_ins_len {
                    self.max_ins_len = isize as usize;
                }
            }
            i = j;
        }
        if self.types.len() <= 1 {
            return 0;
        }
        self.init_ins_types(site);
        self.init_sequence_context(site, bca);
        self.types.len() as i32
    }

    /// `iaux_init_ins_types`: the majority sequence of each insertion type.
    fn init_ins_types(&mut self, site: &Site<'_, '_>) {
        let mi = self.max_ins_len;
        let nt = self.types.len();
        self.inscns.clear();
        self.inscns.resize(nt * mi, 0);
        if mi == 0 {
            return;
        }
        let mut aux = vec![0u32; 5 * nt * mi];
        for t in 0..nt {
            if self.types[t] <= 0 {
                continue;
            }
            for p in site.plp.iter().flatten() {
                if p.indel != self.types[t] {
                    continue;
                }
                for j in 0..p.indel {
                    let c = SEQ_NT16_INT[p.seqi(p.qpos + j + 1) as usize] as usize;
                    aux[5 * (t * mi + j as usize) + c] += 1;
                }
            }
        }
        for t in 0..nt {
            let mut i = 0;
            while (i as i32) < self.types[t] {
                let tmp = &aux[5 * (t * mi + i)..5 * (t * mi + i) + 5];
                let (mut max, mut max_j) = (tmp[0], 0usize);
                for (j, &v) in tmp.iter().enumerate().skip(1) {
                    if max < v {
                        max = v;
                        max_j = j;
                    }
                }
                self.inscns[t * mi + i] = if max != 0 { max_j as u8 } else { 4 };
                if max_j == 4 {
                    self.types[t] = 0;
                    break;
                }
                i += 1;
            }
        }
    }

    /// `iaux_init_sequence_context`
    fn init_sequence_context(&mut self, site: &Site<'_, '_>, bca: &mut CallAux) {
        let pos = site.pos;
        let r = site.reference;
        let win = i64::from(bca.indel_win_size);
        self.left = if pos > win { pos - win } else { 0 };
        self.right = pos + win;
        if self.types[0] < 0 {
            self.right -= i64::from(self.types[0]);
        }
        let mut i = pos;
        while i < self.right {
            if r.at(i) == 0 {
                break;
            }
            i += 1;
        }
        self.right = i;
        let l_run = indel::l_run(r, pos);
        for (i, &l) in self.types.iter().enumerate() {
            let mut q = bca.openq + bca.extq * (l.abs() - 1);
            let qh = if l_run >= 3 {
                (f64::from(bca.tandemq) * f64::from(l.abs()) / f64::from(l_run) + 0.499) as i32
            } else {
                1000
            };
            if q > qh {
                q = qh;
            }
            self.ref_qual[i] = q.min(255);
        }
        bca.indelreg = 0;
        let mi = self.max_ins_len;
        for (i, &t) in self.types.iter().enumerate() {
            if t == 0 {
                continue;
            }
            let ireg = if t > 0 {
                indel::est_indelreg(pos, r, t, Some(&self.inscns[i * mi..]))
            } else {
                indel::est_indelreg(pos, r, -t, None)
            };
            if ireg > bca.indelreg {
                bca.indelreg = ireg;
            }
        }
    }

    /// `iaux_align_read`
    fn align_read(
        &mut self,
        b: &pasta::pileup::Read,
        zq: Option<&[u8]>,
        nref: usize,
        bca: &CallAux,
    ) -> anyhow::Result<i32> {
        if b.flags() & 0x4 != 0 {
            return Ok(1);
        }
        let mut ref_beg = self.left;
        let ref_end = if self.right < ref_beg + nref as i64 {
            self.right
        } else {
            ref_beg + nref as i64 - 1
        };
        let l_qseq = b.sequence_len() as i64;
        let (qry_off1, qry_off2, ref_off1, ref_off2);
        if ref_beg > b.position() {
            let mut cig = CigarState::new(b);
            qry_off1 = cig.seek_fwd(&mut ref_beg, Trim::Left)?;
            ref_off1 = ref_beg - self.left;
            if ref_beg + (l_qseq - qry_off1) > ref_end {
                qry_off2 = ref_end - ref_beg + qry_off1;
                ref_off2 = ref_end - self.left;
            } else {
                qry_off2 = l_qseq - 1;
                ref_off2 = ref_off1 + qry_off2 - qry_off1;
            }
        } else {
            qry_off1 = 0;
            ref_off1 = b.position() - ref_beg;
            if b.position() + l_qseq - 1 > ref_end {
                ref_off2 = ref_end - self.left;
                qry_off2 = ref_off2 - ref_off1;
            } else {
                qry_off2 = l_qseq - 1;
                ref_off2 = ref_off1 + qry_off2 - qry_off1;
            }
        }
        // bcftools asserts these offsets, but release builds go on: windows that leave nothing
        // to align give lengths of zero or less, which probaln_glocal scores as 0 or INT_MIN.
        let qlen = qry_off2 - qry_off1 + 1;
        let rlen = ref_off2 - ref_off1 + 1;
        let score = if qlen < 0 || rlen < 0 {
            i32::MIN
        } else if qlen == 0 || rlen == 0 {
            0
        } else {
            let in_bounds =
                qry_off1 >= 0 && qry_off2 < l_qseq && ref_off1 >= 0 && (ref_off2 as usize) < nref;
            if !in_bounds {
                return Err(anyhow::Error::new(SkipSite(format!(
                    "bam2bcf_iaux: cannot align {} at {}: query {qry_off1}-{qry_off2} of {l_qseq}, reference {ref_off1}-{ref_off2} of {nref}",
                    String::from_utf8_lossy(b.name()),
                    self.left
                ))));
            }
            let (q1, q2) = (qry_off1 as usize, qry_off2 as usize);
            self.qry_seq.clear();
            self.qry_seq
                .extend((q1..=q2).map(|i| SEQ_NT16_INT[b.base(i) as usize]));
            let qq = b.quality_scores();
            self.qual.clear();
            for i in q1..=q2 {
                let v = match zq {
                    Some(bq) => i32::from(qq[i]) + (i32::from(bq[i]) - 64),
                    None => i32::from(qq[i]),
                };
                // A uint8_t clamp of the sum.
                let v = (v as u8).clamp(7, 30);
                self.qual.push(v);
            }
            let r1 = ref_off1 as usize;
            indel::probaln_glocal(
                &self.ref_seq[r1..r1 + rlen as usize],
                &self.qry_seq,
                &self.qual,
                1e-4,
                1e-2,
                10,
                &mut self.probaln,
            )
        };
        let adj = c_int(100.0 * f64::from(score) / qlen as f64 + 0.499) as f32 * bca.indel_bias;
        let adj = c_int(f64::from(adj));
        Ok(score.wrapping_shl(8) | adj.min(255))
    }

    /// `iaux_score_reads`
    fn score_reads(
        &mut self,
        plp: &[Plp<'_>],
        zqs: &[Option<&[u8]>],
        itype: usize,
        bca: &CallAux,
    ) -> anyhow::Result<()> {
        let nt = self.types.len();
        let t = self.types[itype];
        let mi = self.max_ins_len;
        for which in 0..2 {
            let cns = std::mem::take(&mut self.rcns.cns[which]);
            if cns.seq.is_empty() {
                self.rcns.cns[which] = cns;
                break;
            }
            let nseq = cns.seq.len() as i64;
            let ipos = cns.ipos as i64;
            let ref_len = (nseq + i64::from(t)).max(0) as usize;
            self.ref_seq.clear();
            self.ref_seq
                .extend_from_slice(&cns.seq[..(ipos + 1).min(nseq) as usize]);
            if t < 0 {
                let from = (ipos + 1 - i64::from(t)) as usize;
                let n = (nseq - ipos - 1 + i64::from(t)).max(0) as usize;
                self.ref_seq.extend_from_slice(
                    &cns.seq[from.min(cns.seq.len())..(from + n).min(cns.seq.len())],
                );
            } else {
                self.ref_seq
                    .extend_from_slice(&self.inscns[itype * mi..itype * mi + t as usize]);
                let from = (ipos + 1) as usize;
                self.ref_seq
                    .extend_from_slice(&cns.seq[from.min(cns.seq.len())..]);
            }
            self.ref_seq.resize(ref_len, 0);
            for (i, p) in plp.iter().enumerate() {
                let aln = self.align_read(p.b, zqs[i], ref_len, bca)?;
                let score = &mut self.read_scores[i * nt + itype];
                if which == 0 || *score > aln {
                    *score = aln;
                }
            }
            self.rcns.cns[which] = cns;
        }
        Ok(())
    }

    /// `iaux_eval_scored_reads`
    fn eval_scored_reads(&mut self, plp: &mut [Plp<'_>]) {
        let nt = self.types.len();
        for (i, p) in plp.iter_mut().enumerate() {
            let score = &self.read_scores[i * nt..(i + 1) * nt];
            let (mut alt_score, mut alt_j) = (i32::MAX, 0usize);
            for (j, &s) in score.iter().enumerate() {
                if j == self.iref_type {
                    continue;
                }
                if alt_score > s {
                    alt_score = s;
                    alt_j = j;
                }
            }
            let ref_score = score[self.iref_type];
            let (sc0, sc1, mut j0) = if alt_score < ref_score {
                (alt_score, ref_score, alt_j)
            } else {
                (ref_score, alt_score, self.iref_type)
            };
            let indel_q = (sc1 >> 8) - (sc0 >> 8);
            let seq_q = self.ref_qual[alt_j];
            let len_norm_q = sc0 & 0xff;
            let mut adj = if len_norm_q > 111 {
                j0 = self.iref_type;
                12
            } else {
                ((1.0 - f64::from(len_norm_q) / 111.0) * f64::from(indel_q) + 0.499) as i32
            };
            if adj > seq_q {
                adj = seq_q;
            }
            p.aux = (j0 as i32) << 16 | seq_q << 8 | adj;
            self.sum_qual[j0] += adj;
        }
    }

    /// `iaux_eval_best_indels`
    #[expect(
        clippy::needless_range_loop,
        reason = "a line-for-line port of the C loops"
    )]
    fn eval_best_indels(&mut self, site: &mut Site<'_, '_>, bca: &mut CallAux) -> i32 {
        let mi = self.max_ins_len;
        bca.maxins = mi as i32;
        bca.inscns.clear();
        bca.inscns.resize(mi * 4, 0);
        let nt = self.types.len();
        let sumq = &mut self.sum_qual;
        for t in 0..nt {
            sumq[t] = sumq[t] << 6 | t as i32;
        }
        for t in 1..nt {
            let mut j = t;
            while j > 0 && sumq[j] > sumq[j - 1] {
                sumq.swap(j, j - 1);
                j -= 1;
            }
        }
        let mut t = 0;
        while t < nt {
            if (sumq[t] & 0x3f) as usize == self.iref_type {
                break;
            }
            t += 1;
        }
        if t > 0 && t < nt {
            let tmp = sumq[t];
            while t > 0 {
                sumq[t] = sumq[t - 1];
                t -= 1;
            }
            sumq[0] = tmp;
        }
        let mut old2new = [0usize; MAX_TYPES];
        for t in 0..nt {
            let itype = (sumq[t] & 0x3f) as usize;
            old2new[itype] = t;
            if t >= 4 {
                continue;
            }
            bca.indel_types[t] = self.types[itype];
            if bca.indel_types[t] <= 0 {
                continue;
            }
            bca.inscns[t * mi..(t + 1) * mi]
                .copy_from_slice(&self.inscns[itype * mi..(itype + 1) * mi]);
        }
        let mut n_alt = 0;
        for p in site.plp.iter_mut().flatten() {
            let itype_old = (p.aux >> 16 & 0x3f) as usize;
            let itype_new = old2new[itype_old];
            p.aux = (itype_new as i32) << 16 | if itype_new >= 4 { 0 } else { p.aux & 0xffff };
            if itype_new > 0 {
                n_alt += 1;
            }
        }
        n_alt
    }

    /// `bcf_iaux_gap_prep`: whether to call an indel at the position.
    ///
    /// # Errors
    ///
    /// Returns bcftools' fatal errors (an unsupported CIGAR operation); a window that would
    /// read outside a read skips the position with a message instead.
    pub fn gap_prep(
        &mut self,
        site: &mut Site<'_, '_>,
        bca: &mut CallAux,
        chr: &[u8],
    ) -> anyhow::Result<bool> {
        match self.try_gap_prep(site, bca, chr) {
            Ok(v) => Ok(v),
            Err(e) => match e.downcast_ref::<SkipSite>() {
                Some(skip) => {
                    eprintln!("{skip}");
                    Ok(false)
                }
                None => Err(e),
            },
        }
    }

    /// [`Self::gap_prep`], returning every error.
    pub fn try_gap_prep(
        &mut self,
        site: &mut Site<'_, '_>,
        bca: &mut CallAux,
        chr: &[u8],
    ) -> anyhow::Result<bool> {
        let ntypes = self.init_types(site, bca, chr);
        if ntypes <= 0 {
            return Ok(false);
        }
        let r: RefSeq<'_> = site.reference;
        let window: &[u8] = r.0.get(self.left as usize..).unwrap_or(&[]);
        for s in 0..site.plp.len() {
            self.calls += 1;
            self.rcns.reset(site.pos, self.left, self.right);
            self.rcns.set_reads(&site.plp[s])?;
            self.rcns.get_consensus(window, &site.plp[s])?;
            let n = site.plp[s].len() * self.types.len();
            self.read_scores.clear();
            self.read_scores.resize(n, 0);
            // The BAQ adjustments of each read (bam_aux_get(b, "ZQ")), looked up once.
            let zqs: Vec<Option<&[u8]>> = site.plp[s].iter().map(super::call::cd::zq).collect();
            for j in 0..self.types.len() {
                self.score_reads(&site.plp[s], &zqs, j, bca)?;
            }
            self.eval_scored_reads(&mut site.plp[s]);
        }
        Ok(self.eval_best_indels(site, bca) > 0)
    }
}
