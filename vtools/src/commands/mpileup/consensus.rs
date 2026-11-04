//! Read consensus for `--indels-2.0`: a port of bcftools' `read_consensus.c` (up to two
//! haplotype consensus sequences around a position) and `cigar_state.h`.

use anyhow::bail;
use pasta::pileup::Read;

use super::{
    call::{Plp, SEQ_NT16_INT, SEQ_NT16_TABLE},
    indel::{CDEL, CHARD_CLIP, CINS, CPAD, CREF_SKIP, CSOFT_CLIP, is_match},
};

/// `NI`: alternative insertions or deletions kept per position.
const NI: usize = 10;
/// `BF_DEL`: the deletion slot of the base frequencies.
const BF_DEL: usize = 5;
/// `MAX_NCVAR`
const MAX_NCVAR: usize = 8;
/// `NHAP`
const NHAP: usize = 1 << MAX_NCVAR;

/// Where [`CigarState::seek_fwd`] moves a position that falls in a deletion (`trim_left`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trim {
    /// To the base after the deletion, trimming the read's left end.
    Left,
    /// To the base before the deletion, trimming the read's right end.
    #[expect(
        dead_code,
        reason = "bcftools 1.24 only trims left; kept with the port of `cstate_seek_fwd`"
    )]
    Right,
}

/// `cigar_state_t`
pub struct CigarState<'a> {
    b: &'a Read,
    icig: usize,
    pub iseq: i64,
    ref_pos: i64,
}

impl<'a> CigarState<'a> {
    pub fn new(b: &'a Read) -> Self {
        Self {
            b,
            icig: 0,
            iseq: 0,
            ref_pos: b.position(),
        }
    }

    /// `cstate_seek_fwd`: the query index of reference position `*pos`.
    ///
    /// # Errors
    ///
    /// Returns an error for an unsupported CIGAR operation.
    pub fn seek_fwd(&mut self, pos: &mut i64, trim: Trim) -> anyhow::Result<i64> {
        let trim_left = trim == Trim::Left;
        let cigar = self.b.cigar();
        while self.ref_pos <= *pos {
            if self.icig >= cigar.len() {
                if trim_left {
                    return Ok(-2);
                }
                *pos = self.ref_pos - 1;
                return Ok(self.iseq - 1);
            }
            let op = cigar[self.icig] & 0xf;
            let len = i64::from(cigar[self.icig] >> 4);
            if is_match(op) {
                if self.ref_pos + len > *pos {
                    return Ok(*pos - self.ref_pos + self.iseq);
                }
                self.ref_pos += len;
                self.iseq += len;
                self.icig += 1;
                continue;
            }
            if op == CINS || op == CSOFT_CLIP {
                self.iseq += len;
                self.icig += 1;
                continue;
            }
            if op == CDEL || op == CREF_SKIP {
                if self.ref_pos + len > *pos {
                    *pos = if trim_left {
                        self.ref_pos + len
                    } else {
                        self.ref_pos - 1
                    };
                    return Ok(if trim_left { self.iseq } else { self.iseq - 1 });
                }
                self.ref_pos += len;
                self.icig += 1;
                continue;
            }
            if op == CHARD_CLIP || op == CPAD {
                self.icig += 1;
                continue;
            }
            bail!("FIXME: not ready for CIGAR operator {op}");
        }
        if trim_left {
            *pos = self.b.position();
            return Ok(0);
        }
        Ok(-1)
    }

    /// `cstate_seek_op_fwd`
    ///
    /// # Errors
    ///
    /// Returns an error for an unsupported CIGAR operation.
    pub fn seek_op_fwd(&mut self, pos: i64, seek_op: u32, oplen: &mut i64) -> anyhow::Result<i64> {
        let cigar = self.b.cigar();
        while self.ref_pos <= pos {
            if self.icig >= cigar.len() {
                return Ok(-2);
            }
            let op = cigar[self.icig] & 0xf;
            let len = i64::from(cigar[self.icig] >> 4);
            if is_match(op) {
                if self.ref_pos + len <= pos {
                    self.ref_pos += len;
                    self.iseq += len;
                    self.icig += 1;
                    continue;
                }
                if seek_op == 0 {
                    return Ok(pos - self.ref_pos + self.iseq);
                }
                return Ok(-1);
            }
            if op == CINS || op == CSOFT_CLIP {
                if self.ref_pos == pos && seek_op == op {
                    *oplen = len;
                    return Ok(self.iseq);
                }
                if self.ref_pos >= pos {
                    return Ok(-1);
                }
                self.iseq += len;
                self.icig += 1;
                continue;
            }
            if op == CDEL || op == CREF_SKIP {
                if self.ref_pos == pos && seek_op == op {
                    *oplen = len;
                    return Ok(self.iseq);
                }
                if self.ref_pos >= pos {
                    return Ok(-1);
                }
                self.ref_pos += len;
                self.icig += 1;
                continue;
            }
            if op == CHARD_CLIP || op == CPAD {
                self.icig += 1;
                continue;
            }
            bail!("FIXME: not ready for CIGAR operator {op}");
        }
        Ok(if self.icig < cigar.len() { -1 } else { -2 })
    }
}

#[derive(Clone, Default)]
struct InsFreq {
    seq: [Vec<u8>; NI],
    len: [i32; NI],
    freq: [i32; NI],
    used: usize,
}

#[derive(Clone, Copy, Default)]
struct DelFreq {
    len: [i32; NI],
    freq: [i32; NI],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord)]
enum VType {
    Snv = 0,
    Ins = 1,
    Del = 2,
}

#[derive(Clone, Copy)]
struct CandidateVar {
    vtype: VType,
    pos: i64,
    idx: usize,
    which: usize,
    depth: i32,
    af_dev: f32,
}

/// `cvar_not_preferred`
fn not_preferred(a: &CandidateVar, b: &CandidateVar) -> bool {
    if a.af_dev == b.af_dev {
        return a.depth < b.depth;
    }
    a.af_dev > b.af_dev
}

/// klib's `kheap` with `cvar_not_preferred`: the root is the least preferred variant. As in
/// klib, deleting leaves the array's elements in place: the root of an empty heap is the last
/// element deleted, which bcftools may compare against.
#[derive(Default)]
struct Heap {
    dat: Vec<CandidateVar>,
    ndat: usize,
}

impl Heap {
    /// `khp_insert`
    fn insert(&mut self, v: CandidateVar) {
        self.ndat += 1;
        if self.ndat > self.dat.len() {
            self.dat.push(v);
        }
        let mut i = self.ndat - 1;
        while i > 0 && not_preferred(&v, &self.dat[(i - 1) / 2]) {
            self.dat[i] = self.dat[(i - 1) / 2];
            i = (i - 1) / 2;
        }
        self.dat[i] = v;
    }

    /// `khp_heapify`
    fn heapify(&mut self, i: usize) {
        let n = self.ndat;
        let l = 2 * i + 1;
        let r = 2 * i + 2;
        let mut extreme = if l < n && not_preferred(&self.dat[l], &self.dat[i]) {
            l
        } else {
            i
        };
        if r < n && not_preferred(&self.dat[r], &self.dat[extreme]) {
            extreme = r;
        }
        if extreme != i {
            self.dat.swap(i, extreme);
            self.heapify(extreme);
        }
    }

    /// `khp_delete`
    fn delete(&mut self) {
        if self.ndat == 0 {
            return;
        }
        self.ndat -= 1;
        self.dat[0] = self.dat[self.ndat];
        self.heapify(0);
    }
}

/// `cns_seq_t`: a consensus sequence (0-4) and the index of the position in it.
#[derive(Clone, Default)]
pub struct CnsSeq {
    pub seq: Vec<u8>,
    pub ipos: usize,
}

/// `read_cns_t`
#[derive(Default)]
pub struct ReadCns {
    pos: i64,
    beg: i64,
    end: i64,
    band: i32,
    max_del: i32,
    base_freq: Vec<[i32; 6]>,
    ins_freq: Vec<InsFreq>,
    del_freq: Vec<DelFreq>,
    heap: Heap,
    /// `ncvar`, which bcftools does not reset between calls: a call starts counting its
    /// candidates from the previous call's number, so fewer of them fit.
    ncvar: usize,
    cvar: Vec<CandidateVar>,
    hap_freq: Vec<i32>,
    cns_hap: [u32; 2],
    ncns: usize,
    pub cns: [CnsSeq; 2],
}

impl ReadCns {
    /// `rcns_reset`
    pub fn reset(&mut self, pos: i64, beg: i64, end: i64) {
        self.band = 0;
        self.pos = pos;
        self.beg = beg;
        self.end = end;
        for c in &mut self.cns {
            c.seq.clear();
            c.ipos = 0;
        }
        while self.heap.ndat > 0 {
            self.heap.delete();
        }
        let n = (end - beg + 1).max(0) as usize;
        self.base_freq.clear();
        self.base_freq.resize(n, [0; 6]);
        self.del_freq.clear();
        self.del_freq.resize(n, DelFreq::default());
        self.ins_freq.clear();
        self.ins_freq.resize(n, InsFreq::default());
    }

    fn add_base(&mut self, ref_pos: i64, nt16: u8) {
        let i = (ref_pos - self.beg) as usize;
        self.base_freq[i][SEQ_NT16_INT[nt16 as usize] as usize] += 1;
    }

    fn add_ins(&mut self, ref_pos: i64, seq_pos: i64, b: &Read, len: i64) {
        let i = (ref_pos - self.beg) as usize;
        let s: Vec<u8> = (0..len).map(|k| b.base((k + seq_pos) as usize)).collect();
        let ifrq = &mut self.ins_freq[i];
        let mut k = 0;
        while k < NI && k < ifrq.used {
            if ifrq.len[k] == len as i32 && ifrq.seq[k] == s {
                break;
            }
            k += 1;
        }
        if k >= NI {
            return;
        }
        if k == ifrq.used {
            ifrq.seq[k] = s;
            ifrq.len[k] = len as i32;
            ifrq.used += 1;
        }
        ifrq.freq[k] += 1;
    }

    fn add_del(&mut self, ref_pos: i64, len: i64) {
        let i = (ref_pos - self.beg) as usize;
        let mut n = (self.end - self.beg + 1) as usize;
        if i + (len as usize) + 1 < n {
            n = i + len as usize + 1;
        }
        for j in i + 1..n {
            self.base_freq[j][BF_DEL] += 1;
        }
        let dfrq = &mut self.del_freq[i];
        let mut k = 0;
        while k < NI && dfrq.len[k] != 0 {
            if i64::from(dfrq.len[k]) == len {
                break;
            }
            k += 1;
        }
        if k >= NI {
            return;
        }
        if dfrq.len[k] == 0 {
            dfrq.len[k] = len as i32;
        }
        dfrq.freq[k] += 1;
    }

    /// The candidate count the next call starts from.
    pub fn ncvar(&self) -> usize {
        self.ncvar
    }

    /// Sets the candidate count the next call starts from.
    pub fn set_ncvar(&mut self, ncvar: usize) {
        self.ncvar = ncvar;
    }

    /// `rcns_set_reads`
    ///
    /// # Errors
    ///
    /// Returns an error for an unsupported CIGAR operation.
    pub fn set_reads(&mut self, plp: &[Plp<'_>]) -> anyhow::Result<()> {
        for p in plp {
            let b = p.b;
            let (mut x, mut y) = (b.position(), 0i64);
            let mut local_band = 0i32;
            let mut local_band_max = 0i32;
            for &c in b.cigar() {
                let op = c & 0xf;
                let len = i64::from(c >> 4);
                if op == CSOFT_CLIP {
                    y += len;
                } else if is_match(op) {
                    if x < self.end && x + len > self.beg {
                        let j_beg = if self.beg > x { self.beg - x } else { 0 };
                        let j_end = if self.end < x + len - 1 {
                            self.end - x
                        } else {
                            len - 1
                        };
                        x += j_beg;
                        y += j_beg;
                        for _ in j_beg..=j_end {
                            self.add_base(x, b.base(y as usize));
                            x += 1;
                            y += 1;
                        }
                    } else {
                        x += len;
                        y += len;
                    }
                } else if op == CINS {
                    if x > self.beg && x < self.end {
                        local_band += p.indel;
                        self.add_ins(x - 1, y, b, len);
                    }
                    y += len;
                } else if op == CDEL || op == CREF_SKIP {
                    if x > self.beg && x + len - 1 <= self.end {
                        local_band += -p.indel;
                        self.add_del(x - 1, len);
                    }
                    x += len;
                } else if op == CHARD_CLIP {
                    continue;
                } else {
                    bail!("rcns_set_reads todo: unknown cigar operator {op}");
                }
                if local_band_max < local_band {
                    local_band_max = local_band;
                }
            }
            if self.band < local_band_max {
                self.band = local_band_max;
            }
        }
        Ok(())
    }

    fn register_variant(&mut self, var: CandidateVar) {
        let mut var = var;
        let free_slot;
        if self.ncvar == MAX_NCVAR {
            if let Some(root) = self.heap.dat.first()
                && not_preferred(&var, root)
            {
                return;
            }
            free_slot = self.heap.dat.first().map_or(0, |root| root.idx);
            self.heap.delete();
        } else {
            free_slot = self.ncvar;
            self.ncvar += 1;
        }
        var.idx = free_slot;
        self.heap.insert(var);
    }

    /// `select_candidate_variants`
    fn select_candidate_variants(&mut self, r: &[u8]) {
        const AF_TH: f32 = 0.1;
        let n = (self.end - self.beg + 1) as usize;
        let at = |i: usize| r.get(i).copied().unwrap_or(0);
        let mut i = 0;
        while i < n && at(i) != 0 {
            if i as i64 == self.pos - self.beg {
                i += 1;
                continue;
            }
            let bf = self.base_freq[i];
            let df = self.del_freq[i];
            let ifr = &self.ins_freq[i];
            let mut dp = bf[..4].iter().sum::<i32>();
            let mut j = 0;
            while j < NI && df.len[j] != 0 {
                dp += df.freq[j];
                j += 1;
            }
            for j in 0..ifr.used.min(NI) {
                dp += ifr.freq[j];
            }
            let mut vars = Vec::new();
            for (j, &cnt) in bf[..4].iter().enumerate() {
                if cnt == 0 || at(i) == b"ACGTN"[j] {
                    continue;
                }
                let af = cnt as f32 / dp as f32;
                if af > AF_TH && af < (1.0 - AF_TH) {
                    vars.push((VType::Snv, j, af));
                }
            }
            let mut j = 0;
            while j < NI && df.len[j] != 0 {
                let af = df.freq[j] as f32 / dp as f32;
                if af > AF_TH && af < (1.0 - AF_TH) {
                    vars.push((VType::Del, j, af));
                }
                j += 1;
            }
            for j in 0..ifr.used.min(NI) {
                let af = ifr.freq[j] as f32 / dp as f32;
                if af > AF_TH && af < (1.0 - AF_TH) {
                    vars.push((VType::Ins, j, af));
                }
            }
            for (vtype, which, af) in vars {
                self.register_variant(CandidateVar {
                    vtype,
                    pos: i as i64 + self.beg,
                    idx: 0,
                    which,
                    depth: dp,
                    af_dev: (0.5 - f64::from(af)).abs() as f32,
                });
            }
            i += 1;
        }
        // done: take the heap's variants, sorted by position, type, and which.
        self.cvar.clear();
        while self.heap.ndat > 0 {
            self.cvar.push(self.heap.dat[0]);
            self.heap.delete();
        }
        self.ncvar = self.cvar.len();
        self.cvar.sort_by(|a, b| {
            a.pos
                .cmp(&b.pos)
                .then(a.vtype.cmp(&b.vtype))
                .then(a.which.cmp(&b.which))
        });
        let i = (self.pos - self.beg) as usize;
        self.max_del = 0;
        if let Some(df) = self.del_freq.get(i) {
            let mut j = 0;
            while j < NI && df.len[j] != 0 {
                self.max_del = self.max_del.max(df.len[j]);
                j += 1;
            }
        }
    }

    /// `create_haplotype_frequency_spectrum`
    fn haplotype_frequencies(&mut self, plp: &[Plp<'_>]) -> anyhow::Result<()> {
        self.hap_freq.clear();
        self.hap_freq.resize(NHAP, 0);
        for p in plp {
            let mut cig = CigarState::new(p.b);
            let mut hap = 0usize;
            for (j, cvar) in self.cvar.iter().enumerate() {
                let mut len = 0i64;
                match cvar.vtype {
                    VType::Snv => {
                        let iseq = cig.seek_op_fwd(cvar.pos, 0, &mut len)?;
                        if iseq == -2 {
                            break;
                        }
                        if iseq == -1 {
                            continue;
                        }
                        let nt16 = p.b.base(iseq as usize);
                        if SEQ_NT16_INT[nt16 as usize] as usize == cvar.which {
                            hap |= 1 << j;
                        }
                    }
                    VType::Ins => {
                        let ifrq = &self.ins_freq[(cvar.pos - self.beg) as usize];
                        let iseq = cig.seek_op_fwd(cvar.pos + 1, CINS, &mut len)?;
                        if iseq == -2 {
                            break;
                        }
                        if iseq == -1 {
                            continue;
                        }
                        if len != i64::from(ifrq.len[cvar.which]) {
                            continue;
                        }
                        let s = &ifrq.seq[cvar.which];
                        if (0..s.len()).all(|k| p.b.base(iseq as usize + k) == s[k]) {
                            hap |= 1 << j;
                        }
                    }
                    VType::Del => {
                        let dfrq = &self.del_freq[(cvar.pos - self.beg) as usize];
                        let ret = cig.seek_op_fwd(cvar.pos + 1, CDEL, &mut len)?;
                        if ret == -2 {
                            break;
                        }
                        if ret == -1 {
                            continue;
                        }
                        if len != i64::from(dfrq.len[cvar.which]) {
                            continue;
                        }
                        hap |= 1 << j;
                    }
                }
            }
            self.hap_freq[hap] += 1;
        }
        Ok(())
    }

    /// `correct_haplotype_errors`
    fn correct_haplotype_errors(&mut self) {
        let mut freq: Vec<(usize, i32)> = (0..NHAP).map(|i| (i, self.hap_freq[i])).collect();
        let tot: i32 = self.hap_freq.iter().sum();
        // glibc's qsort is a stable merge sort here.
        freq.sort_by_key(|f| std::cmp::Reverse(f.1));
        let mut i = NHAP - 1;
        while i > 1 {
            if freq[i].1 == 0 {
                i -= 1;
                continue;
            }
            if freq[1].1 > tot - freq[0].1 - freq[1].1 {
                break;
            }
            let mut count = freq[i].1;
            let mut max_hap = 0usize;
            for j in 0..MAX_NCVAR {
                if freq[i].0 & (1 << j) == 0 {
                    continue;
                }
                let hap = freq[i].0 ^ (1 << j);
                if count < self.hap_freq[hap] {
                    count = self.hap_freq[hap];
                    max_hap = hap;
                }
            }
            if count == freq[i].1 {
                i -= 1;
                continue;
            }
            let count = freq[i].1;
            freq[i].1 = 0;
            self.hap_freq[freq[i].0] = 0;
            self.hap_freq[max_hap] += count;
            for j in i + 1..NHAP {
                if freq[j].1 == 0 {
                    break;
                }
                freq.swap(j - 1, j);
            }
            let mut j = i as isize - 1;
            while j >= 0 {
                let ju = j as usize;
                if freq[ju].0 == max_hap {
                    freq[ju].1 += count;
                }
                if freq[ju].1 < freq[ju + 1].1 {
                    freq.swap(ju, ju + 1);
                }
                j -= 1;
            }
            i -= 1;
        }
        let f0 = freq[0].1 as f32;
        let f1 = freq[1].1 as f32;
        self.ncns = if f64::from(f1 / (f0 + f1)) < 0.1 {
            1
        } else {
            2
        };
        let hap0 = freq[0].0 as u32;
        let hap1 = if self.ncns == 2 { freq[1].0 as u32 } else { 0 };
        self.cns_hap = [0, 0];
        let mut kept = Vec::with_capacity(self.cvar.len());
        for i in 0..MAX_NCVAR {
            if (hap0 | hap1) & (1 << i) == 0 {
                continue;
            }
            let j = kept.len();
            kept.push(self.cvar[i]);
            if hap0 & (1 << i) != 0 {
                self.cns_hap[0] |= 1 << j;
            }
            if hap1 & (1 << i) != 0 {
                self.cns_hap[1] |= 1 << j;
            }
        }
        self.cvar = kept;
        self.ncvar = self.cvar.len();
    }

    /// `apply_consensus_insertion`
    fn apply_consensus_insertion(&self, cns: &mut CnsSeq, j: usize, ivar: usize) {
        let ref_pos = self.beg + j as i64;
        if self.pos == ref_pos {
            return;
        }
        let mut ivar = ivar;
        while ivar < self.cvar.len() && self.cvar[ivar].pos == ref_pos {
            if self.cvar[ivar].vtype == VType::Ins {
                return;
            }
            ivar += 1;
        }
        let bf = &self.base_freq[j];
        let ifr = &self.ins_freq[j];
        let nreads: i32 = bf[..BF_DEL].iter().sum();
        let (mut max_freq, mut kmax) = (0, 0);
        for k in 0..ifr.used.min(NI) {
            if max_freq < ifr.freq[k] {
                max_freq = ifr.freq[k];
                kmax = k;
            }
        }
        if nreads > max_freq * 2 {
            return;
        }
        for &c in &ifr.seq[kmax] {
            cns.seq.push(SEQ_NT16_INT[c as usize]);
        }
    }

    /// `create_consensus`
    fn create_consensus(&self, r: &[u8], ith: usize) -> CnsSeq {
        let n = (self.end - self.beg + 1) as usize;
        let mut cns = CnsSeq::default();
        let mut prev_pos = 0i64;
        let mut ivar = 0usize;
        let at = |i: usize| r.get(i).copied().unwrap_or(0);
        let base_of_ref = |j: usize| SEQ_NT16_INT[SEQ_NT16_TABLE[at(j) as usize] as usize] as usize;
        let mut j = 0usize;
        while j < n {
            let ref_pos = self.beg + j as i64;
            if self.pos == ref_pos {
                cns.ipos = cns.seq.len();
            }
            while ivar < self.cvar.len() && self.cvar[ivar].pos < ref_pos {
                ivar += 1;
            }
            let bf = &self.base_freq[j];
            if ivar >= self.cvar.len() || self.cvar[ivar].pos != ref_pos {
                let (mut max_freq, mut kmax) = (0, base_of_ref(j));
                let nk = if ref_pos < self.pos || ref_pos > self.pos + i64::from(self.max_del) {
                    BF_DEL + 1
                } else {
                    BF_DEL
                };
                for (k, &f) in bf[..nk].iter().enumerate() {
                    if max_freq < f {
                        max_freq = f;
                        kmax = k;
                    }
                }
                if kmax != BF_DEL {
                    prev_pos = ref_pos;
                    cns.seq.push(kmax as u8);
                }
                self.apply_consensus_insertion(&mut cns, j, ivar);
                j += 1;
                continue;
            }
            let cv = self.cvar[ivar];
            let which = cv.which;
            if self.cns_hap[ith] & (1 << ivar) == 0 {
                let (mut max_freq, mut kmax) = (0, base_of_ref(j));
                for (k, &f) in bf.iter().enumerate() {
                    if cv.vtype == VType::Snv && cv.which == k {
                        continue;
                    }
                    if max_freq < f {
                        max_freq = f;
                        kmax = k;
                    }
                }
                if kmax != BF_DEL && (cns.seq.is_empty() || prev_pos != ref_pos) {
                    prev_pos = ref_pos;
                    cns.seq.push(kmax as u8);
                }
                self.apply_consensus_insertion(&mut cns, j, ivar);
                j += 1;
                continue;
            }
            if cv.vtype == VType::Snv {
                prev_pos = ref_pos;
                cns.seq.push(which as u8);
                self.apply_consensus_insertion(&mut cns, j, ivar);
                j += 1;
                continue;
            }
            if cns.seq.is_empty() || prev_pos != ref_pos {
                let (mut max_freq, mut kmax) = (0, base_of_ref(j));
                for (k, &f) in bf.iter().enumerate() {
                    if cv.vtype == VType::Snv && cv.which == k {
                        continue;
                    }
                    if max_freq < f {
                        max_freq = f;
                        kmax = k;
                    }
                }
                if kmax != BF_DEL {
                    prev_pos = ref_pos;
                    cns.seq.push(kmax as u8);
                }
            }
            if cv.vtype == VType::Ins {
                for &c in &self.ins_freq[j].seq[which] {
                    prev_pos = ref_pos;
                    cns.seq.push(SEQ_NT16_INT[c as usize]);
                }
            } else if cv.vtype == VType::Del {
                j += self.del_freq[j].len[which] as usize;
            }
            j += 1;
        }
        cns
    }

    /// `rcns_get_consensus`: up to two consensus sequences (the second empty when one
    /// suffices). `r` is the reference from the window start.
    ///
    /// # Errors
    ///
    /// Returns an error for an unsupported CIGAR operation.
    pub fn get_consensus(&mut self, r: &[u8], plp: &[Plp<'_>]) -> anyhow::Result<()> {
        self.select_candidate_variants(r);
        if self.cvar.is_empty() {
            self.cns_hap[0] = 0;
            self.ncns = 1;
        } else {
            self.haplotype_frequencies(plp)?;
            self.correct_haplotype_errors();
        }
        for i in 0..self.ncns {
            self.cns[i] = self.create_consensus(r, i);
        }
        if self.ncns < 2 {
            self.cns[1] = CnsSeq::default();
        }
        Ok(())
    }
}
