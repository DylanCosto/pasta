//! Indel likelihoods: a port of bcftools' `bam2bcf_indel.c` (the default model), with the
//! helpers it shares with the `--indels-2.0` and `--indels-cns` models, `str_finder.c`, and
//! the forward pass of HTSlib's `probaln_glocal`.

use std::sync::atomic::{AtomicUsize, Ordering};

use super::call::{CallAux, INDEL_NULL, NQUAL, Plp, SEQ_NT16_INT, SEQ_NT16_TABLE, c_int};

/// `MINUS_CONST`
pub const MINUS_CONST: u32 = 0x1000_0000;
/// `MAX_TYPES`
pub const MAX_TYPES: usize = 64;
/// Reads longer than this are aligned in a narrower window (and, in the default model, with
/// higher gap probabilities).
pub const LONG_READ_LEN: i32 = 1000;
/// The alignment score of a read that cannot be aligned to a candidate.
pub const NO_SCORE: i32 = 0xff_ffff;

/// A reference sequence as the C code sees it: NUL past the end.
#[derive(Clone, Copy)]
pub struct RefSeq<'a>(pub &'a [u8]);

impl RefSeq<'_> {
    #[inline]
    pub fn at(&self, i: i64) -> u8 {
        if i < 0 {
            return 0;
        }
        self.0.get(i as usize).copied().unwrap_or(0)
    }

    #[inline]
    pub fn len(&self) -> i64 {
        self.0.len() as i64
    }
}

// CIGAR operations.
pub const CMATCH: u32 = 0;
pub const CINS: u32 = 1;
pub const CDEL: u32 = 2;
pub const CREF_SKIP: u32 = 3;
pub const CSOFT_CLIP: u32 = 4;
pub const CHARD_CLIP: u32 = 5;
pub const CPAD: u32 = 6;
pub const CEQUAL: u32 = 7;
pub const CDIFF: u32 = 8;

#[inline]
pub fn is_match(op: u32) -> bool {
    op == CMATCH || op == CEQUAL || op == CDIFF
}

/// `bam_cigar2qlen`
pub fn cigar2qlen(cigar: &[u32]) -> i32 {
    cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 1 | 4 | 7 | 8))
        .map(|&c| (c >> 4) as i32)
        .sum()
}

/// `bam_cigar2rlen`
pub fn cigar2rlen(cigar: &[u32]) -> i32 {
    cigar
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|&c| (c >> 4) as i32)
        .sum()
}

/// The end of a deletion that [`tpos2qpos`] gives for a position inside it (`is_left`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DelEnd {
    /// The deletion's first position.
    Left,
    /// The position after the deletion.
    Right,
}

/// `tpos2qpos`: the query position of reference position `tpos`; returns it and the
/// reference position used (`*_tpos`).
pub fn tpos2qpos(pos: i32, cigar: &[u32], tpos: i32, del_end: DelEnd) -> (i32, i32) {
    let (mut x, mut y, mut last_y) = (pos, 0i32, 0i32);
    for &c in cigar {
        let op = c & 0xf;
        let l = (c >> 4) as i32;
        if is_match(op) {
            if pos > tpos {
                return (y, pos);
            }
            if x + l > tpos {
                return (y + (tpos - x), tpos);
            }
            x += l;
            y += l;
            last_y = y;
        } else if op == CINS || op == CSOFT_CLIP {
            y += l;
        } else if op == CDEL || op == CREF_SKIP {
            if x + l > tpos {
                let end = match del_end {
                    DelEnd::Left => x,
                    DelEnd::Right => x + l,
                };
                return (y, end);
            }
            x += l;
        }
    }
    (last_y, x)
}

/// `est_indelreg`: the length of reference sequence after `pos` the indel is ambiguous over.
pub fn est_indelreg(pos: i64, r: RefSeq<'_>, l: i32, ins4: Option<&[u8]>) -> i32 {
    let l = l.unsigned_abs() as i64;
    let (mut max, mut max_i, mut score) = (0i32, pos, 0i32);
    let mut i = pos + 1;
    let mut j = 0i64;
    while r.at(i) != 0 {
        let c = r.at(i).to_ascii_uppercase();
        let other = match ins4 {
            Some(ins) => b"ACGTN"[ins[(j % l) as usize] as usize],
            None => r.at(pos + 1 + j % l).to_ascii_uppercase(),
        };
        score += if c != other { -10 } else { 1 };
        if score < 0 {
            break;
        }
        if max < score {
            max = score;
            max_i = i;
        }
        i += 1;
        j += 1;
    }
    (max_i - pos) as i32
}

/// `get_pos`: the soft clip length weight, sequence length without clips, and scaled
/// position of the indel in the read.
pub fn get_pos(p: &Plp<'_>) -> (i32, i32, i32) {
    let cigar = p.b.cigar();
    let (mut sc_len, mut sc_dist, mut at_left) = (0i32, -1i32, true);
    let mut epos = p.qpos;
    let mut slen = p.l_qseq();
    for &c in cigar {
        let op = c & 0xf;
        let len = (c >> 4) as i32;
        if op == CSOFT_CLIP {
            slen -= len;
            if at_left {
                sc_len += len;
                epos -= sc_len;
                sc_dist = epos;
            } else {
                let rd = p.l_qseq() - len - p.qpos;
                if sc_dist < 0 || sc_dist > rd {
                    sc_dist = rd;
                    sc_len = len;
                }
            }
        } else if op != CHARD_CLIP {
            at_left = false;
        }
    }
    if p.indel > 0 && slen - (epos + p.indel) < epos {
        epos += p.indel - 1;
    }
    let epos_r = (f64::from(epos) / f64::from(slen + 1) * super::call::NPOS as f64) as i32;
    let sc_len_r = if sc_len != 0 {
        ((15.0 * f64::from(sc_len) / f64::from(sc_dist + 1)) as i32).min(99)
    } else {
        0
    };
    (sc_len_r, slen, epos_r)
}

/// `bcf_cgp_l_run`: the length of the homopolymer run around `pos`.
pub fn l_run(r: RefSeq<'_>, pos: i64) -> i32 {
    let c = SEQ_NT16_TABLE[r.at(pos + 1) as usize];
    if c == 15 {
        return 1;
    }
    let mut i = pos + 2;
    while r.at(i) != 0 {
        if SEQ_NT16_TABLE[r.at(i) as usize] != c {
            break;
        }
        i += 1;
    }
    let mut l_run = i;
    let mut i = pos;
    while i >= 0 {
        if SEQ_NT16_TABLE[r.at(i) as usize] != c {
            break;
        }
        i -= 1;
    }
    l_run -= i + 1;
    l_run as i32
}

/// `bcf_cgp_calc_cons`: the majority consensus of each insertion type; types whose consensus
/// has an N are set to 0.
pub fn calc_cons(plp: &[Vec<Plp<'_>>], types: &mut [i32], max_ins: i32) -> Vec<u8> {
    let n_types = types.len();
    let max_ins = max_ins as usize;
    let mut aux = vec![0i32; 5 * n_types * max_ins];
    for t in 0..n_types {
        if types[t] <= 0 {
            continue;
        }
        for p in plp.iter().flatten() {
            if p.indel == types[t] {
                for k in 1..=p.indel {
                    let c = SEQ_NT16_INT[p.seqi(p.qpos + k) as usize] as usize;
                    aux[(t * max_ins + (k as usize - 1)) * 5 + c] += 1;
                }
            }
        }
    }
    let mut inscns = vec![0u8; n_types * max_ins];
    for t in 0..n_types {
        let mut j = 0;
        while (j as i32) < types[t] {
            let ia = &aux[(t * max_ins + j) * 5..(t * max_ins + j) * 5 + 5];
            let (mut max, mut max_k) = (0, -1i32);
            for (k, &v) in ia.iter().enumerate() {
                if v > max {
                    max = v;
                    max_k = k as i32;
                }
            }
            inscns[t * max_ins + j] = if max != 0 { max_k as u8 } else { 4 };
            if max_k == 4 {
                types[t] = 0;
                break;
            }
            j += 1;
        }
    }
    inscns
}

/// `find_STR`: the short tandem repeats (homopolymers up to 8-mers) of `cons`, as
/// (start, end, repeat length).
pub fn find_str(cons: &[u8]) -> Vec<(i32, i32, i32)> {
    let len = cons.len();
    let mut reps: Vec<(i32, i32, i32)> = Vec::new();
    let mut w: u32 = 0;
    let mut i = 0;
    let mut j = 0;
    while i < len && j < 15 {
        w = (w << 2) | u32::from(cons[i]);
        for (rlen, jmin) in (1..=7).zip([1, 3, 5, 7, 9, 11, 13]) {
            let bits = 2 * rlen;
            let mask = (1u32 << bits) - 1;
            if j >= jmin && (w & mask) == ((w >> bits) & mask) {
                add_rep(&mut reps, cons, i, rlen);
            }
        }
        j += 1;
        i += 1;
    }
    while i < len {
        w = (w << 2) | u32::from(cons[i]);
        for rlen in (1..=8).rev() {
            let bits = 2 * rlen;
            let mask = if bits == 32 {
                u32::MAX
            } else {
                (1u32 << bits) - 1
            };
            let hi = if bits >= 32 { 0 } else { (w >> bits) & mask };
            if (w & mask) == hi {
                add_rep(&mut reps, cons, i, rlen);
                break;
            }
        }
        i += 1;
    }
    reps
}

/// `add_rep` (without pads or lower-case filtering, which bcftools does not use here).
fn add_rep(reps: &mut Vec<(i32, i32, i32)>, cons: &[u8], pos: usize, rlen: i32) {
    let pos = pos as i32;
    if let Some(&(start, end, _)) = reps.last()
        && start <= pos - rlen * 2 + 1
        && end >= pos
    {
        return;
    }
    let clen = cons.len() as i32;
    let mut cp1 = pos - rlen + 1;
    let mut cp2 = pos + 1;
    while cp2 < clen {
        if cons[cp1 as usize] != cons[cp2 as usize] {
            break;
        }
        cp1 += 1;
        cp2 += 1;
    }
    let end = pos + cp2 - (pos + 1);
    let start = pos + 1 - 2 * rlen;
    // Remove older entries entirely contained in the new one, from the tail.
    let mut k = reps.len();
    while k > 0 {
        let (s, e, _) = reps[k - 1];
        if e < start {
            break;
        }
        if s >= start {
            reps.remove(k - 1);
        }
        k -= 1;
    }
    reps.push((start, end, rlen));
}

/// `g_qual2prob`: `pow(10, -i/10.)` in single precision.
static QUAL2PROB: std::sync::LazyLock<[f32; 256]> =
    std::sync::LazyLock::new(|| std::array::from_fn(|i| 10f64.powf(-(i as f64) / 10.0) as f32));

/// Scratch rows for [`probaln_glocal`].
#[derive(Default)]
pub struct ProbalnBuf {
    rows: [Vec<f64>; 2],
    s: Vec<f64>,
    qual: Vec<f32>,
}

/// `probaln_glocal` without the backward pass (`state` and `q` NULL): the phred-scaled
/// likelihood of the query given the reference, or `i32::MIN` on invalid input.
pub fn probaln_glocal(
    reference: &[u8],
    query: &[u8],
    iqual: &[u8],
    d: f32,
    e: f32,
    bw_par: i32,
    buf: &mut ProbalnBuf,
) -> i32 {
    const EI: f64 = 0.25;
    const EM: f64 = 0.33333333333;
    let l_ref = reference.len() as i32;
    let l_query = query.len() as i32;
    if l_query >= i32::MAX - 2 {
        return i32::MIN;
    }
    if l_ref == 0 || l_query == 0 {
        return 0;
    }
    let mut bw = l_ref.max(l_query);
    if bw > bw_par {
        bw = bw_par;
    }
    if bw < (l_ref - l_query).abs() {
        bw = (l_ref - l_query).abs();
    }
    let bw2 = bw * 2 + 1;
    let i_dim = if bw2 < l_ref {
        bw2 as usize * 3 + 6
    } else {
        l_ref as usize * 3 + 6
    };
    let set_u = |i: i32, k: i32| -> usize {
        let x = (i - bw).max(0);
        ((k - x + 1) * 3) as usize
    };

    let table = &*QUAL2PROB;
    buf.qual.clear();
    buf.qual.extend(iqual.iter().map(|&q| table[q as usize]));
    let qual = &buf.qual;
    buf.s.clear();
    buf.s.resize(l_query as usize + 2, 0.0);
    let s = &mut buf.s;
    for row in &mut buf.rows {
        row.clear();
        row.resize(i_dim, 0.0);
    }

    let sm = 1.0 / (2.0 * f64::from(l_query) + 2.0);
    let si = sm;
    let mut m = [0f64; 9];
    m[0] = f64::from(1.0 - d - d) * (1.0 - sm);
    m[1] = f64::from(d) * (1.0 - sm);
    m[2] = m[1];
    m[3] = f64::from(1.0 - e) * (1.0 - si);
    m[4] = f64::from(e) * (1.0 - si);
    m[5] = 0.0;
    m[6] = f64::from(1.0 - e);
    m[7] = 0.0;
    m[8] = f64::from(e);
    let bm = f64::from((1.0 - d) / l_ref as f32);
    let bi = f64::from(d / l_ref as f32);

    // f[0] lives in rows[0]; f[1] in rows[1].
    s[0] = 1.0;
    {
        let fi = &mut buf.rows[1];
        let end = if l_ref < bw + 1 { l_ref } else { bw + 1 };
        let mut sum = 0.0;
        for k in 1..=end {
            let rk = reference[k as usize - 1];
            let ev = if rk > 3 || query[0] > 3 {
                1.0
            } else if rk == query[0] {
                1.0 - f64::from(qual[0])
            } else {
                f64::from(qual[0]) * EM
            };
            let u = set_u(1, k);
            fi[u] = ev * bm;
            fi[u + 1] = EI * bi;
            sum += fi[u] + fi[u + 1];
        }
        s[1] = sum;
    }
    // Rows f[i-1] and f[i] alternate between two buffers. A buffer holds stale values from
    // two rows before, so the cells a row reads but the previous row did not write (the
    // triples just outside its band, zero in HTSlib's calloc'ed matrix) are cleared.
    let mut cur = 1usize;
    let (mut beg, mut end) = (1, if l_ref < bw + 1 { l_ref } else { bw + 1 });
    for i in 2..=l_query {
        cur ^= 1;
        let (lo, hi) = buf.rows.split_at_mut(1);
        let (fi, fi1) = if cur == 0 {
            (&mut lo[0], &hi[0])
        } else {
            (&mut hi[0], &lo[0])
        };
        let qli = f64::from(qual[i as usize - 1]);
        let qyi = query[i as usize - 1];
        beg = 1.max(i - bw);
        end = l_ref.min(i + bw);
        let ev = [qli * EM, 1.0 - qli, 1.0, 1.0];
        let mm = 1.0 / s[i as usize - 1];
        let xm = [
            mm * m[0],
            mm * m[3],
            mm * m[6],
            EI * mm * m[1],
            EI * mm * m[4],
        ];
        let n = (end - beg + 1) as usize;
        let u0 = set_u(i, beg);
        let v0 = set_u(i - 1, beg - 1);
        let refs = &reference[beg as usize - 1..end as usize];
        let out = &mut fi[u0..u0 + 3 * n];
        let inp = &fi1[v0..v0 + 3 * n + 2];
        let (mut l_x0, mut l_x2) = (0.0f64, 0.0f64);
        let mut sum = 0.0;
        for ((o, y), &rk) in out
            .as_chunks_mut::<3>()
            .0
            .iter_mut()
            .zip(inp.windows(5).step_by(3))
            .zip(refs)
        {
            let cond = usize::from(rk > 3 || qyi > 3) * 2 + usize::from(rk == qyi);
            let z0 = xm[0] * y[0];
            let z1 = xm[1] * y[1];
            let z2 = xm[2] * y[2];
            let z3 = xm[3] * y[3];
            let z4 = xm[4] * y[4];
            let x0 = ev[cond] * (z0 + z1 + z2);
            let x1 = z3 + z4;
            let x2 = l_x0 + l_x2;
            o[0] = x0;
            o[1] = x1;
            o[2] = x2;
            sum += x0 + x1 + x2;
            l_x0 = m[2] * x0;
            l_x2 = m[8] * x2;
        }
        // Clear the triples next to this row's band, which the next row may read.
        if u0 >= 3 {
            fi[u0 - 3..u0].fill(0.0);
        }
        let after = u0 + 3 * n;
        let stop = (after + 3).min(fi.len());
        if after < stop {
            fi[after..stop].fill(0.0);
        }
        s[i as usize] = sum;
    }
    {
        let last = &buf.rows[cur];
        let mm = 1.0 / s[l_query as usize];
        let mut sum = 0.0;
        // Cells outside the last row's band are zero in HTSlib and add nothing.
        for k in beg..=end {
            let u = set_u(l_query, k);
            if u < 3 || u >= i_dim {
                continue;
            }
            sum += mm * last[u] * sm + mm * last[u + 1] * si;
        }
        s[l_query as usize + 1] = sum;
    }
    let mut p = 1.0f64;
    let mut pr1 = 0.0f64;
    for &si in s.iter().take(l_query as usize + 2) {
        p *= si;
        if p < 1e-100 {
            pr1 += -4.343 * p.ln();
            p = 1.0;
        }
    }
    pr1 += -4.343 * (p * f64::from(l_ref) * f64::from(l_query)).ln();
    // A NaN or infinite likelihood converts to INT_MIN, as on x86-64.
    c_int(pr1 + 0.499)
}

/// The data the indel callers share for one position: the reads of each sample.
pub struct Site<'p, 'a> {
    pub plp: &'p mut [Vec<Plp<'a>>],
    pub pos: i64,
    pub reference: RefSeq<'p>,
}

/// `bcf_cgp_find_types` of `bam2bcf_indel.c`: the indel sizes present, sorted; `None` to
/// skip the position.
fn find_types(site: &Site<'_, '_>, bca: &mut CallAux) -> Option<FoundTypes> {
    let n_reads: usize = site.plp.iter().map(Vec::len).sum();
    bca.max_support = 0;
    bca.max_frac = 0.0;
    let mut aux: Vec<u32> = Vec::with_capacity(n_reads + 1);
    aux.push(MINUS_CONST);
    let (mut max_rd_len, mut n_alt, mut n_tot) = (0i32, 0i32, 0i32);
    let mut indel_support_ok = false;
    for s in site.plp.iter() {
        let (mut na, mut nt) = (0i32, 0i32);
        for p in s {
            nt += 1;
            if p.indel != 0 {
                na += 1;
                aux.push(MINUS_CONST.wrapping_add(p.indel as u32));
            }
            let j = cigar2qlen(p.b.cigar());
            if j > max_rd_len {
                max_rd_len = j;
            }
        }
        let frac = f64::from(na) / f64::from(nt);
        if !indel_support_ok && na as u32 >= bca.min_support && frac >= bca.min_frac {
            indel_support_ok = true;
        }
        if na as u32 > bca.max_support && frac > 0.0 {
            bca.max_support = na as u32;
            bca.max_frac = frac as f32;
        }
        n_alt += na;
        n_tot += nt;
    }
    aux.sort_unstable();
    let mut n_types = 1;
    for i in 1..aux.len() {
        if aux[i] != aux[i - 1] {
            n_types += 1;
        }
    }
    if !bca.per_sample_flt {
        indel_support_ok = !(f64::from(n_alt) / f64::from(n_tot) < bca.min_frac
            || (n_alt as u32) < bca.min_support);
    }
    if n_types == 1 || !indel_support_ok {
        return None;
    }
    if n_types >= MAX_TYPES {
        eprintln!(
            "[bcf_cgp_find_types] excessive INDEL alleles at position {}. Skip the position.",
            site.pos + 1
        );
        return None;
    }
    let win = 2 * bca.indel_win_size;
    let i_end = site.pos + i64::from(if win < max_rd_len { win } else { max_rd_len });
    let mut n_n = 0i64;
    let mut i = site.pos;
    while i < i_end && site.reference.at(i) != 0 {
        if site.reference.at(i) == b'N' {
            n_n += 1;
        }
        i += 1;
    }
    if n_n * 2 > i - site.pos {
        return None;
    }
    let mut types = Vec::with_capacity(n_types);
    types.push(aux[0].wrapping_sub(MINUS_CONST) as i32);
    for i in 1..aux.len() {
        if aux[i] != aux[i - 1] {
            types.push(aux[i].wrapping_sub(MINUS_CONST) as i32);
        }
    }
    let ref_type = types.iter().position(|&t| t == 0).unwrap_or(types.len());
    Some(FoundTypes {
        types,
        ref_type,
        max_rd_len,
        n_reads,
    })
}

/// The indel types of a position.
pub struct FoundTypes {
    pub types: Vec<i32>,
    pub ref_type: usize,
    pub max_rd_len: i32,
    pub n_reads: usize,
}

/// `bcf_cgp_ref_sample`: a per-sample reference with mismatches present in at least 30% of
/// the reads at up to two loci masked as N (15).
#[expect(
    clippy::needless_range_loop,
    reason = "a line-for-line port of the C loops"
)]
fn ref_sample(site: &Site<'_, '_>, left: i64, right: i64) -> Vec<Vec<u8>> {
    let l = (right - left + 1) as usize;
    let mut ref0 = vec![0u8; l];
    for i in 0..(right - left) as usize {
        ref0[i] = SEQ_NT16_TABLE[site.reference.at(i as i64 + left) as usize];
    }
    let mut out = Vec::with_capacity(site.plp.len());
    let mut cns = vec![0u32; l];
    for s in site.plp.iter() {
        let mut r = vec![0u8; l];
        cns.fill(0);
        for p in s {
            let (mut x, mut y) = (i64::from(p.b.position() as i32), 0i64);
            for &c in p.b.cigar() {
                let op = c & 0xf;
                let len = i64::from(c >> 4);
                if is_match(op) {
                    if x + len >= left {
                        let mut j = (left - x).max(0);
                        let j_end = (right - x).min(len);
                        while j < j_end {
                            let base = p.b.base((y + j) as usize);
                            cns[(x + j - left) as usize] += if base == ref0[(x + j - left) as usize]
                            {
                                1
                            } else {
                                1 << 16
                            };
                            j += 1;
                        }
                    }
                    x += len;
                    y += len;
                } else if op == CDEL || op == CREF_SKIP {
                    x += len;
                } else if op == CINS || op == CSOFT_CLIP {
                    y += len;
                }
                if x > right {
                    break;
                }
            }
        }
        r[..(right - left) as usize].copy_from_slice(&ref0[..(right - left) as usize]);
        let (mut max, mut max2) = (0u32, 0u32);
        let (mut max_i, mut max2_i) = (-1i64, -1i64);
        for i in 0..(right - left) as usize {
            if cns[i] >> 16 >= max >> 16 {
                max2 = max;
                max2_i = max_i;
                max = cns[i];
                max_i = i as i64;
            } else if cns[i] >> 16 >= max2 >> 16 {
                max2 = cns[i];
                max2_i = i as i64;
            }
        }
        let frac = |v: u32| f64::from(v & 0xffff) / f64::from((v & 0xffff) + (v >> 16));
        if frac(max) >= 0.7 {
            max_i = -1;
        }
        if frac(max2) >= 0.7 {
            max2_i = -1;
        }
        if max_i >= 0 {
            r[max_i as usize] = 15;
        }
        if max2_i >= 0 {
            r[max2_i as usize] = 15;
        }
        out.push(r);
    }
    out
}

/// `est_seqQ` of `bam2bcf_indel.c`
fn est_seq_q(bca: &CallAux, l: i32, l_run: i32) -> i32 {
    let q = bca.openq + bca.extq * (l.abs() - 1);
    let qh = if l_run >= 3 {
        (f64::from(bca.tandemq) * f64::from(l.abs()) / f64::from(l_run) + 0.499) as i32
    } else {
        1000
    };
    q.min(qh)
}

/// Scratch buffers for the default model.
#[derive(Default)]
pub struct IndelBuf {
    pub probaln: ProbalnBuf,
    query: Vec<u8>,
    ref2: Vec<u8>,
    qq: Vec<u8>,
    score: Vec<i32>,
    jobs: Vec<AlignJob>,
}

/// `bcf_cgp_align_score` of `bam2bcf_indel.c`, with the read's coordinates in `job`.
#[expect(clippy::too_many_arguments, reason = "a line-for-line port")]
fn align_score(
    p: &Plp<'_>,
    bca: &CallAux,
    ty: i32,
    seg: &[u8],
    query: &[u8],
    job: &AlignJob,
    zq: Option<&[u8]>,
    qq: &mut Vec<u8>,
    pbuf: &mut ProbalnBuf,
) -> i32 {
    let AlignJob {
        r_start,
        r_end,
        long_read,
        tbeg,
        qbeg,
        qend,
        qpos,
        ..
    } = *job;
    let (d, e) = if long_read {
        (1e-3f32, 1e-1f32)
    } else {
        (1e-4f32, 1e-2f32)
    };
    let ty = ty.abs();
    let bw = ty + 3;
    let quals = p.b.quality_scores();
    qq.clear();
    for l in qbeg..qend {
        let l = l as usize;
        let qval = match zq {
            Some(bq) => i32::from(quals[l]) + (i32::from(bq[l]) - 64),
            None => i32::from(quals[l]),
        };
        qq.push(qval.clamp(7, 30) as u8);
    }
    let sc = probaln_glocal(seg, query, qq, d, e, bw, pbuf);
    if sc < 0 {
        return NO_SCORE;
    }
    let l = (c_int(100.0 * f64::from(sc) / f64::from(qend - qbeg) + 0.499)) as f32 * bca.indel_bias;
    let mut score = sc << 8 | c_int(f64::from(l)).min(255);

    let mut iscore = 0;
    for (start, end, rep_len) in find_str(seg) {
        if start <= qpos && end >= qpos {
            iscore += (end - start) / rep_len;
            if start + tbeg <= r_start || end + tbeg >= r_end {
                iscore += 2 * (end - start);
            }
        }
    }
    let l = (f64::from(score & 0xff) * 0.8 + f64::from(iscore * 2)) as i32;
    score = (score & !0xff) | l.min(255);
    score
}

/// One read to realign against a candidate haplotype (the arguments of
/// `bcf_cgp_align_score`).
#[derive(Debug, Clone, Copy)]
struct AlignJob {
    slot: usize,
    read: usize,
    zq: usize,
    r_start: i32,
    r_end: i32,
    long_read: bool,
    tbeg: i32,
    qbeg: i32,
    qend: i32,
    qpos: i32,
    seg_off: usize,
    seg_len: usize,
}

/// Alignment cells above which a position's realignments are spread over threads.
const PARALLEL_CELLS: usize = 4 << 20;

/// Scores the reads of one sample against one candidate haplotype in `buf.ref2`, on several
/// threads when there is much to do (a few positions, such as long repeats, take most of the
/// time) and threads of the budget are free.
fn run_align_jobs(
    jobs: &[AlignJob],
    plp: &[Plp<'_>],
    ty: i32,
    zqs: &[Option<&[u8]>],
    bca: &CallAux,
    buf: &mut IndelBuf,
) {
    let cells: usize = jobs
        .iter()
        .map(|j| (j.qend - j.qbeg).max(0) as usize * (ty.unsigned_abs() as usize * 2 + 7))
        .sum();
    let ref2 = &buf.ref2;
    let score_one = |j: &AlignJob, query: &mut Vec<u8>, qq: &mut Vec<u8>, pb: &mut ProbalnBuf| {
        let p = &plp[j.read];
        query.clear();
        query.extend((j.qbeg..j.qend).map(|l| SEQ_NT16_INT[p.seqi(l) as usize]));
        let seg = &ref2[j.seg_off..(j.seg_off + j.seg_len).min(ref2.len())];
        align_score(p, bca, ty, seg, query, j, zqs[j.zq], qq, pb)
    };
    let (mut query, mut qq, mut pb) = (
        std::mem::take(&mut buf.query),
        std::mem::take(&mut buf.qq),
        std::mem::take(&mut buf.probaln),
    );
    if cells < PARALLEL_CELLS || jobs.len() < 2 {
        for j in jobs {
            let score = score_one(j, &mut query, &mut qq, &mut pb);
            buf.score[j.slot] = score;
        }
        (buf.query, buf.qq, buf.probaln) = (query, qq, pb);
        return;
    }
    // The jobs are taken in turn by this thread and by helpers on the threads of the budget
    // it shares with the other pileup workers, which join as they become free.
    let next = AtomicUsize::new(0);
    let take = |query: &mut Vec<u8>, qq: &mut Vec<u8>, pb: &mut ProbalnBuf| {
        let k = next.fetch_add(1, Ordering::Relaxed);
        jobs.get(k).map(|j| (j.slot, score_one(j, query, qq, pb)))
    };
    let results: Vec<Vec<(usize, i32)>> = std::thread::scope(|scope| {
        let mut helpers = Vec::new();
        let mut mine = Vec::new();
        let mut can_start = true;
        loop {
            let left = jobs.len().saturating_sub(next.load(Ordering::Relaxed));
            while can_start && helpers.len() + 1 < left {
                let held = bca.threads.take_up_to(1);
                if held.count() == 0 {
                    break;
                }
                // Without a thread to start (the budget comes back with `held`), this thread
                // takes the jobs on its own.
                let helper = std::thread::Builder::new().spawn_scoped(scope, move || {
                    let _held = held;
                    let (mut query, mut qq, mut pb) =
                        (Vec::new(), Vec::new(), ProbalnBuf::default());
                    std::iter::from_fn(|| take(&mut query, &mut qq, &mut pb)).collect()
                });
                match helper {
                    Ok(helper) => helpers.push(helper),
                    Err(_) => can_start = false,
                }
            }
            match take(&mut query, &mut qq, &mut pb) {
                Some(score) => mine.push(score),
                None => break,
            }
        }
        std::iter::once(mine)
            .chain(
                helpers
                    .into_iter()
                    .map(|h| h.join().expect("an alignment thread panicked")),
            )
            .collect()
    });
    (buf.query, buf.qq, buf.probaln) = (query, qq, pb);
    for (slot, score) in results.into_iter().flatten() {
        buf.score[slot] = score;
    }
}

/// `bcf_cgp_compute_indelQ` of `bam2bcf_indel.c`: sets each read's `aux` and the site's
/// indel types. Returns the number of reads supporting an alternate allele.
#[expect(clippy::too_many_arguments, reason = "a line-for-line port")]
fn compute_indel_q(
    site: &mut Site<'_, '_>,
    bca: &mut CallAux,
    inscns: &[u8],
    l_run: i32,
    max_ins: i32,
    ref_type: usize,
    types: &[i32],
    score: &[i32],
) -> i32 {
    let n_types = types.len();
    let mut sumq = [0i32; MAX_TYPES];
    let mut sc = [0i32; MAX_TYPES];
    let mut k = 0usize;
    for s in site.plp.iter_mut() {
        for p in s.iter_mut() {
            let sct = &score[k * n_types..(k + 1) * n_types];
            k += 1;
            for t in 0..n_types {
                sc[t] = sct[t] << 6 | t as i32;
            }
            for t in 1..n_types {
                let mut j = t;
                while j > 0 && sc[j] < sc[j - 1] {
                    sc.swap(j, j - 1);
                    j -= 1;
                }
            }
            let (mut indel_q, seq_q);
            if (sc[0] & 0x3f) as usize == ref_type {
                indel_q = (sc[1] >> 14) - (sc[0] >> 14);
                seq_q = est_seq_q(bca, types[(sc[1] & 0x3f) as usize], l_run);
            } else {
                let mut t = 0;
                while t < n_types {
                    if (sc[t] & 0x3f) as usize == ref_type {
                        break;
                    }
                    t += 1;
                }
                indel_q = (sc[t.min(MAX_TYPES - 1)] >> 14) - (sc[0] >> 14);
                seq_q = est_seq_q(bca, types[(sc[0] & 0x3f) as usize], l_run);
            }
            let tmp = sc[0] >> 6 & 0xff;
            indel_q = if tmp > 111 {
                0
            } else {
                ((1.0 - f64::from(tmp) / 111.0) * f64::from(indel_q) + 0.499) as i32
            };
            if indel_q > seq_q {
                indel_q = seq_q;
            }
            let indel_q = indel_q.min(255);
            let seq_q = seq_q.min(255);
            p.aux = (sc[0] & 0x3f) << 16 | seq_q << 8 | indel_q;
            sumq[(sc[0] & 0x3f) as usize] += indel_q.min(seq_q);
        }
    }
    finish_types(site, bca, inscns, max_ins, ref_type, types, &mut sumq)
}

/// The end of `bcf_cgp_compute_indelQ` shared with the edlib model: picks the reference and
/// up to three alternate types by summed quality and re-labels the reads.
#[expect(
    clippy::needless_range_loop,
    reason = "a line-for-line port of the C loops"
)]
pub fn finish_types(
    site: &mut Site<'_, '_>,
    bca: &mut CallAux,
    inscns: &[u8],
    max_ins: i32,
    ref_type: usize,
    types: &[i32],
    sumq: &mut [i32; MAX_TYPES],
) -> i32 {
    let n_types = types.len();
    bca.maxins = max_ins;
    bca.inscns.clear();
    bca.inscns.resize(max_ins as usize * 4, 0);
    for t in 0..n_types {
        sumq[t] = sumq[t] << 6 | t as i32;
    }
    for t in 1..n_types {
        let mut j = t;
        while j > 0 && sumq[j] > sumq[j - 1] {
            sumq.swap(j, j - 1);
            j -= 1;
        }
    }
    let mut t = 0;
    while t < n_types {
        if (sumq[t] & 0x3f) as usize == ref_type {
            break;
        }
        t += 1;
    }
    if t > 0 && t < MAX_TYPES {
        let tmp = sumq[t];
        while t > 0 {
            sumq[t] = sumq[t - 1];
            t -= 1;
        }
        sumq[0] = tmp;
    }
    bca.indel_types = [INDEL_NULL; 4];
    let mi = max_ins as usize;
    for t in 0..4.min(n_types) {
        let it = (sumq[t] & 0x3f) as usize;
        bca.indel_types[t] = types[it];
        if mi > 0 {
            bca.inscns[t * mi..(t + 1) * mi].copy_from_slice(&inscns[it * mi..(it + 1) * mi]);
        }
    }
    let mut n_alt = 0;
    for p in site.plp.iter_mut().flatten() {
        let x = types[(p.aux >> 16 & 0x3f) as usize];
        let mut j = 0;
        while j < 4 {
            if x == bca.indel_types[j] {
                break;
            }
            j += 1;
        }
        p.aux = (j as i32) << 16 | if j == 4 { 0 } else { p.aux & 0xffff };
        if (p.aux >> 16 & 0x3f) > 0 {
            n_alt += 1;
        }
    }
    n_alt
}

/// Records the indel statistics of a read (the `t == 0` block of the gap preparations).
pub fn indel_stats(p: &Plp<'_>, bca: &mut CallAux) {
    let nqual_over_60 = NQUAL as f64 / 60.0;
    let imq = (f64::from(p.mapq().min(59)) * nqual_over_60) as usize;
    let (sc_len, _slen, epos) = get_pos(p);
    let (sc_len, epos) = (sc_len as usize, epos as usize);
    // bcftools asserts (in debug builds only) that these are in range: an insertion at the end
    // of a read can put `epos` past the end, where release builds write over the next array.
    // Such counts are dropped here.
    let bump = |counts: &mut [i32], i: usize| {
        if let Some(c) = counts.get_mut(i) {
            *c += 1;
        }
    };
    if p.indel != 0 {
        bca.ialt_mq[imq] += 1;
        bump(&mut bca.ialt_scl, sc_len);
        bump(&mut bca.ialt_pos, epos);
    } else {
        bca.iref_mq[imq] += 1;
        bump(&mut bca.iref_scl, sc_len);
        bump(&mut bca.iref_pos, epos);
    }
}

/// Whether the pileup has a read with an indel after the position.
pub fn has_indel(plp: &[Vec<Plp<'_>>]) -> bool {
    plp.iter().flatten().any(|p| p.indel != 0)
}

/// `bcf_call_gap_prep`: the default indel model. Returns whether to call an indel here.
pub fn gap_prep(site: &mut Site<'_, '_>, bca: &mut CallAux, buf: &mut IndelBuf) -> bool {
    if !has_indel(site.plp) {
        return false;
    }
    let Some(found) = find_types(site, bca) else {
        return false;
    };
    let FoundTypes {
        mut types,
        ref_type,
        max_rd_len: _,
        n_reads,
    } = found;
    let n_types = types.len();
    let pos = site.pos;
    let r = site.reference;
    let win = i64::from(bca.indel_win_size);
    let left = if pos > win { pos - win } else { 0 };
    let mut right = pos + win;
    if types[0] < 0 {
        right -= i64::from(types[0]);
    }
    let mut i = pos;
    while i < right {
        if r.at(i) == 0 {
            break;
        }
        i += 1;
    }
    right = i;

    let ref_sample = ref_sample(site, left, right);
    let l_run = l_run(r, pos);
    let max_ins = types[n_types - 1];
    let inscns = if max_ins > 0 {
        calc_cons(site.plp, &mut types, max_ins)
    } else {
        Vec::new()
    };
    let mi = max_ins.max(0) as usize;

    let max_ref2 = (right - left + 2 + 2 * i64::from(max_ins.max(-types[0]))) as usize;
    buf.ref2.clear();
    buf.ref2.resize(max_ref2, 0);
    buf.score.clear();
    buf.score.resize(n_reads * n_types, 0);
    bca.indelreg = 0;
    // The BAQ adjustments of each read (bam_aux_get(b, "ZQ")), looked up once.
    let zqs: Vec<Option<&[u8]>> = site.plp.iter().flatten().map(super::call::cd::zq).collect();

    for t in 0..n_types {
        let ir = if types[t] == 0 {
            0
        } else if types[t] > 0 {
            est_indelreg(pos, r, types[t], Some(&inscns[t * mi..]))
        } else {
            est_indelreg(pos, r, -types[t], None)
        };
        if ir > bca.indelreg {
            bca.indelreg = ir;
        }

        let mut k_read = 0usize;
        for s in 0..site.plp.len() {
            let mut k = 0usize;
            let mut j = left;
            while j <= pos {
                buf.ref2[k] = SEQ_NT16_INT[ref_sample[s][(j - left) as usize] as usize];
                k += 1;
                j += 1;
            }
            if types[t] <= 0 {
                j += i64::from(-types[t]);
            } else {
                for l in 0..types[t] as usize {
                    buf.ref2[k] = inscns[t * mi + l];
                    k += 1;
                }
            }
            while j < right && r.at(j) != 0 {
                buf.ref2[k] = SEQ_NT16_INT[ref_sample[s][(j - left) as usize] as usize];
                k += 1;
                j += 1;
            }
            while k < max_ref2 {
                buf.ref2[k] = 4;
                k += 1;
            }
            if right > j {
                right = j;
            }

            let mut jobs = std::mem::take(&mut buf.jobs);
            jobs.clear();
            for i in 0..site.plp[s].len() {
                let p = site.plp[s][i];
                let kk = k_read;
                k_read += 1;
                if t == 0 {
                    indel_stats(&p, bca);
                }
                if p.b.flags() & 0x4 != 0 {
                    continue;
                }
                let cigar = p.b.cigar();
                if cigar.iter().any(|&c| c & 0xf == CREF_SKIP) {
                    continue;
                }
                let long_read = p.l_qseq() > LONG_READ_LEN;
                let (mut left2, mut right2) = (left, right);
                if long_read {
                    if pos - left >= win {
                        left2 += win / 2;
                    }
                    if right - pos >= win {
                        right2 -= win / 2;
                    }
                }
                let bpos = p.b.position() as i32;
                let r_start = bpos;
                let r_end = cigar2rlen(cigar) - 1 + r_start;
                let (qbeg, mut tbeg) = tpos2qpos(bpos, cigar, left2 as i32, DelEnd::Right);
                let (qp, _) = tpos2qpos(bpos, cigar, pos as i32, DelEnd::Right);
                let qpos = qp - qbeg;
                let (qend, tend) = tpos2qpos(bpos, cigar, right2 as i32, DelEnd::Left);
                if types[t] < 0 {
                    let l = -types[t];
                    tbeg = if i64::from(tbeg - l) > left {
                        tbeg - l
                    } else {
                        left as i32
                    };
                }
                let slot = kk * n_types + t;
                if tend > tbeg {
                    jobs.push(AlignJob {
                        slot,
                        read: i,
                        zq: kk,
                        r_start,
                        r_end,
                        long_read,
                        tbeg,
                        qbeg,
                        qend,
                        qpos,
                        seg_off: (i64::from(tbeg) - left) as usize,
                        seg_len: (tend - tbeg + types[t].abs()) as usize,
                    });
                } else {
                    buf.score[slot] = NO_SCORE;
                }
            }
            run_align_jobs(&jobs, &site.plp[s], types[t], &zqs, bca, buf);
            buf.jobs = jobs;
        }
    }
    let score = std::mem::take(&mut buf.score);
    let n_alt = compute_indel_q(site, bca, &inscns, l_run, max_ins, ref_type, &types, &score);
    buf.score = score;
    n_alt > 0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// HTSlib's `probaln_glocal` forward pass over a whole calloc'ed matrix, as written.
    fn probaln_full(reference: &[u8], query: &[u8], iqual: &[u8], d: f32, e: f32, bw0: i32) -> i32 {
        const EI: f64 = 0.25;
        const EM: f64 = 0.33333333333;
        let (l_ref, l_query) = (reference.len() as i32, query.len() as i32);
        if l_ref == 0 || l_query == 0 {
            return 0;
        }
        let mut bw = l_ref.max(l_query).min(bw0);
        if bw < (l_ref - l_query).abs() {
            bw = (l_ref - l_query).abs();
        }
        let bw2 = bw * 2 + 1;
        let i_dim = if bw2 < l_ref {
            bw2 * 3 + 6
        } else {
            l_ref * 3 + 6
        } as usize;
        let mut f = vec![0f64; (l_query as usize + 1) * i_dim];
        let mut s = vec![0f64; l_query as usize + 2];
        let qual: Vec<f32> = iqual.iter().map(|&q| QUAL2PROB[q as usize]).collect();
        let set_u = |i: i32, k: i32| (((k - (i - bw).max(0)) + 1) * 3) as usize;
        let sm = 1.0 / (2.0 * f64::from(l_query) + 2.0);
        let si = sm;
        let m = [
            f64::from(1.0 - d - d) * (1.0 - sm),
            f64::from(d) * (1.0 - sm),
            f64::from(d) * (1.0 - sm),
            f64::from(1.0 - e) * (1.0 - si),
            f64::from(e) * (1.0 - si),
            0.0,
            f64::from(1.0 - e),
            0.0,
            f64::from(e),
        ];
        let bm = f64::from((1.0 - d) / l_ref as f32);
        let bi = f64::from(d / l_ref as f32);
        s[0] = 1.0;
        let mut sum = 0.0;
        for k in 1..=l_ref.min(bw + 1) {
            let rk = reference[k as usize - 1];
            let ev = if rk > 3 || query[0] > 3 {
                1.0
            } else if rk == query[0] {
                1.0 - f64::from(qual[0])
            } else {
                f64::from(qual[0]) * EM
            };
            let u = i_dim + set_u(1, k);
            f[u] = ev * bm;
            f[u + 1] = EI * bi;
            sum += f[u] + f[u + 1];
        }
        s[1] = sum;
        for i in 2..=l_query {
            let qli = f64::from(qual[i as usize - 1]);
            let qyi = query[i as usize - 1];
            let ev = [qli * EM, 1.0 - qli, 1.0, 1.0];
            let mm = 1.0 / s[i as usize - 1];
            let mut sum = 0.0;
            for k in 1.max(i - bw)..=l_ref.min(i + bw) {
                let rk = reference[k as usize - 1];
                let cond = usize::from(rk > 3 || qyi > 3) * 2 + usize::from(rk == qyi);
                let u = i as usize * i_dim + set_u(i, k);
                let v11 = (i as usize - 1) * i_dim + set_u(i - 1, k - 1);
                let v10 = (i as usize - 1) * i_dim + set_u(i - 1, k);
                let v01 = i as usize * i_dim + set_u(i, k - 1);
                f[u] = ev[cond]
                    * (mm * m[0] * f[v11] + mm * m[3] * f[v11 + 1] + mm * m[6] * f[v11 + 2]);
                f[u + 1] = EI * mm * m[1] * f[v10] + EI * mm * m[4] * f[v10 + 1];
                f[u + 2] = m[2] * f[v01] + m[8] * f[v01 + 2];
                sum += f[u] + f[u + 1] + f[u + 2];
            }
            s[i as usize] = sum;
        }
        let mm = 1.0 / s[l_query as usize];
        let mut sum = 0.0;
        for k in 1..=l_ref {
            let u = set_u(l_query, k);
            if u < 3 || u >= i_dim {
                continue;
            }
            let u = l_query as usize * i_dim + u;
            sum += mm * f[u] * sm + mm * f[u + 1] * si;
        }
        s[l_query as usize + 1] = sum;
        let (mut p, mut pr1) = (1.0f64, 0.0f64);
        for &x in &s {
            p *= x;
            if p < 1e-100 {
                pr1 += -4.343 * p.ln();
                p = 1.0;
            }
        }
        pr1 += -4.343 * (p * f64::from(l_ref) * f64::from(l_query)).ln();
        (pr1 + 0.499) as i32
    }

    #[test]
    fn banded_probaln_matches_the_full_matrix() {
        let mut rng = crate::commands::mpileup::errmod::Rng48::new(11);
        let mut next = |n: u32| (rng.next_f64() * f64::from(n)) as u32;
        let mut buf = ProbalnBuf::default();
        for _ in 0..2000 {
            let l_ref = 1 + next(80) as usize;
            let l_query = 1 + next(80) as usize;
            let reference: Vec<u8> = (0..l_ref).map(|_| next(5) as u8).collect();
            let query: Vec<u8> = (0..l_query).map(|_| next(5) as u8).collect();
            let qual: Vec<u8> = (0..l_query).map(|_| 7 + next(24) as u8).collect();
            let bw = 1 + next(12) as i32;
            let want = probaln_full(&reference, &query, &qual, 1e-4, 1e-2, bw);
            let got = probaln_glocal(&reference, &query, &qual, 1e-4, 1e-2, bw, &mut buf);
            assert_eq!(got, want, "ref {reference:?} query {query:?} bw {bw}");
        }
    }

    #[test]
    fn finds_homopolymers() {
        // AAAAA: a homopolymer from 0 to 4.
        let reps = find_str(&[0, 0, 0, 0, 0, 1, 2]);
        assert!(
            reps.iter().any(|&(s, e, l)| s == 0 && e == 4 && l == 1),
            "{reps:?}"
        );
    }

    #[test]
    fn tpos2qpos_through_a_deletion() {
        // 5M2D5M at 100
        let cigar = [5 << 4, 2 << 4 | CDEL, 5 << 4];
        assert_eq!(tpos2qpos(100, &cigar, 102, DelEnd::Right), (2, 102));
        assert_eq!(tpos2qpos(100, &cigar, 105, DelEnd::Left), (5, 105));
        assert_eq!(tpos2qpos(100, &cigar, 105, DelEnd::Right), (5, 107));
        assert_eq!(tpos2qpos(100, &cigar, 200, DelEnd::Right), (10, 112));
    }
}
