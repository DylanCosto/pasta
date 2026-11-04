//! The `--indels-cns` indel model: a port of bcftools' `bam2bcf_edlib.c` (consensus per indel
//! type, scored with edit distances) and its cut-down C translation of edlib (`edlib.c`).

use super::{
    call::{CallAux, Plp, SEQ_NT16_INT, SEQ_NT16_TABLE, c_int},
    indel::{
        self, CDEL, CINS, CMATCH, CREF_SKIP, CSOFT_CLIP, DelEnd, FoundTypes, LONG_READ_LEN,
        MAX_TYPES, MINUS_CONST, NO_SCORE, RefSeq, Site, calc_cons, cigar2qlen, cigar2rlen,
        est_indelreg, find_str, finish_types, indel_stats, l_run, tpos2qpos,
    },
};

/// `NI`: alternative insertion sequences kept per position.
const NI: usize = 100;
/// `MAX_INS`
const MAX_INS: usize = 8192;

/// Scratch buffers.
#[derive(Default)]
pub struct EdlibBuf {
    query: Vec<u8>,
    score: Vec<i32>,
}

/// `bcf_cgp_find_types` of `bam2bcf_edlib.c`: like the default model's, but only keeping
/// types with enough support.
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
    let mut i = 0;
    while i < aux.len() {
        let sz = aux[i].wrapping_sub(MINUS_CONST) as i32;
        let mut j = i + 1;
        while j < aux.len() && aux[j] == aux[i] {
            j += 1;
        }
        if sz == 0
            || ((j - i) as u32 >= bca.min_support
                && (bca.per_sample_flt || (j - i) as f64 / f64::from(n_tot) >= bca.min_frac))
        {
            types.push(sz);
        }
        i = j;
    }
    if types.len() <= 1 {
        return None;
    }
    let ref_type = types.iter().position(|&t| t == 0).unwrap_or(types.len());
    Some(FoundTypes {
        types,
        ref_type,
        max_rd_len,
        n_reads,
    })
}

/// `str_freq`: insertion sequences (0-4) and their frequencies.
#[derive(Clone, Default)]
struct StrFreq {
    str: Vec<Vec<u8>>,
    freq: Vec<i32>,
}

impl StrFreq {
    /// `bcf_cgp_append_cons`
    fn append(&mut self, s: &[u8], freq: i32) {
        let mut j = 0;
        while j < NI && j < self.str.len() {
            if self.str[j] == s {
                break;
            }
            j += 1;
        }
        if j >= NI {
            return;
        }
        if j == self.str.len() {
            self.str.push(s.to_vec());
            self.freq.push(0);
        }
        self.freq[j] += freq;
    }
}

/// ASCII `ACGTN*` to 0-5.
fn base6(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        b'*' => 5,
        _ => 4,
    }
}

/// The consensus sequences of one sample and indel type.
struct Consensus {
    cons: [Vec<u8>; 2],
    left_shift: i32,
    right_shift: i32,
    cpos_pos: i32,
}

/// `bcf_cgp_consensus`
#[expect(clippy::too_many_arguments, reason = "a line-for-line port")]
#[expect(
    clippy::needless_range_loop,
    reason = "a line-for-line port of the C loops"
)]
fn consensus(
    plp: &[Plp<'_>],
    pos: i64,
    bca: &CallAux,
    r: RefSeq<'_>,
    left: i64,
    right: i64,
    ty: i32,
    biggest_del: i32,
    band: &mut i32,
    pos_l: i64,
    pos_r: i64,
) -> Consensus {
    let width = (right - left + 1).max(1) as usize;
    let mut cons_base = vec![[0i32; 6]; width];
    let mut cons_ins = vec![StrFreq::default(); width];
    let mut ref_base = vec![[0i32; 6]; width];
    let mut ref_ins = vec![StrFreq::default(); width];
    let mut local_band_max = 0i32;
    let mut total_span_str = 0i32;
    let mut type_depth = 0i32;
    let mut ins = vec![0u8; MAX_INS];
    for p in plp {
        let b = p.b;
        let (mut x, mut y) = (b.position(), 0i64);
        let mut local_band = 0i32;
        for &c in b.cigar() {
            let op = c & 0xf;
            let len = i64::from(c >> 4);
            let mut skip_to = 0i64;
            match op {
                CSOFT_CLIP => y += len,
                0 | 7 | 8 => {
                    let mut j = 0;
                    while j < len {
                        if x >= left {
                            if x >= right {
                                break;
                            }
                            let base = b.base(y as usize);
                            if p.indel == ty {
                                cons_base[(x - left) as usize]
                                    [SEQ_NT16_INT[base as usize] as usize] += 1;
                            } else if x != pos + 1 {
                                ref_base[(x - left) as usize]
                                    [SEQ_NT16_INT[base as usize] as usize] += 1;
                            }
                        }
                        j += 1;
                        x += 1;
                        y += 1;
                    }
                }
                CINS => {
                    if x >= left && x < right {
                        local_band += p.indel;
                        if local_band_max < local_band {
                            local_band_max = local_band;
                        }
                    }
                    let mut j = 0;
                    while j < len {
                        if x >= left {
                            if x >= right {
                                break;
                            }
                            let base = b.base(y as usize);
                            if (j as usize) < MAX_INS {
                                ins[j as usize] = SEQ_NT16_INT[base as usize];
                            }
                        }
                        j += 1;
                        y += 1;
                    }
                    if x >= left && x < right {
                        let ilen = (j as usize).min(MAX_INS);
                        if p.indel == ty {
                            cons_ins[(x - left) as usize].append(&ins[..ilen], 1);
                            type_depth += i32::from(x == pos + 1);
                        } else if x != pos + 1 {
                            ref_ins[(x - left) as usize].append(&ins[..ilen], 1);
                        }
                    }
                }
                CDEL => {
                    if x >= left && x < right {
                        local_band += p.indel;
                        if local_band_max < -local_band {
                            local_band_max = -local_band;
                        }
                    }
                    let mut j = 0;
                    while j < len {
                        if x >= left {
                            if x >= right {
                                break;
                            }
                            if (p.indel == ty && !p.is_del)
                                || (p.indel == 0 && p.is_del && len == i64::from(-ty))
                            {
                                cons_base[(x - left) as usize][5] += 1;
                                type_depth += i32::from(x == pos + 1);
                            } else if x + len <= pos + 1 || (skip_to != 0 && x > skip_to) {
                                ref_base[(x - left) as usize][5] += 1;
                            } else if x <= pos && x + len > pos + 1 && x > skip_to {
                                skip_to = x + len;
                            }
                        }
                        j += 1;
                        x += 1;
                    }
                }
                _ => {}
            }
        }
        if b.position() <= pos_l && x >= pos_r {
            total_span_str += 1;
        }
        if *band < local_band_max {
            *band = local_band_max;
        }
    }

    let span = (right - left).max(0) as usize;
    for i in 0..span {
        let mut t: i32 = cons_base[i].iter().sum();
        t += cons_ins[i].freq.iter().sum::<i32>();
        let mut rr: i32 = ref_base[i].iter().sum();
        rr += ref_ins[i].freq.iter().sum::<i32>();
        let mut rfract = f64::from(rr - t * 2) * 0.75 / f64::from(rr + 1);
        if rfract < 1.01 / (f64::from(rr) + 1e-10) {
            rfract = 1.01 / (f64::from(rr) + 1e-10);
        }
        let ipos = i as i64 + left;
        if ipos > pos && ipos < pos + 1 - i64::from(biggest_del) {
            continue;
        }
        for k in 0..6 {
            cons_base[i][k] =
                (f64::from(cons_base[i][k]) + rfract * f64::from(ref_base[i][k])) as i32;
        }
        let ri = std::mem::take(&mut ref_ins[i]);
        for (s, &f) in ri.str.iter().zip(&ri.freq) {
            cons_ins[i].append(s, (rfract * f64::from(f)) as i32);
        }
    }

    // Merge insertions of the same length into their majority sequence.
    let mut insc = vec![[0i32; 5]; MAX_INS];
    for ci in cons_ins.iter_mut().take(span) {
        for j in 0..ci.str.len() {
            if ci.freq[j] == 0 {
                continue;
            }
            let lj = ci.str[j].len();
            for l in 0..lj {
                insc[l] = [0; 5];
                insc[l][ci.str[j][l] as usize] = ci.freq[j];
            }
            for k in j + 1..ci.str.len() {
                if ci.str[k].len() != lj || ci.freq[k] == 0 {
                    continue;
                }
                for l in 0..lj {
                    insc[l][ci.str[k][l] as usize] += ci.freq[k];
                }
                ci.freq[j] += ci.freq[k];
                ci.freq[k] = 0;
            }
            for l in 0..lj {
                let v = insc[l];
                let tot: i32 = v.iter().sum();
                let (mut max_v, mut base) = (0, 0u8);
                for (bi, &c) in v.iter().enumerate() {
                    if max_v < c {
                        max_v = c;
                        base = bi as u8;
                    }
                }
                ci.str[j][l] = if f64::from(max_v) > 0.6 * f64::from(tot) {
                    base
                } else {
                    4
                };
            }
        }
    }

    const CONS_CUTOFF: f64 = 0.40;
    const CONS_CUTOFF_DEL: f64 = 0.35;
    const CONS_CUTOFF2: f64 = 0.80;
    const CONS_CUTOFF_INC: f64 = 0.35;
    const CONS_CUTOFF_INC2: f64 = 0.80;
    const CONS_CUTOFF_INS: f64 = 0.60;

    let mut out = Consensus {
        cons: [Vec::new(), Vec::new()],
        left_shift: 0,
        right_shift: 0,
        cpos_pos: -1,
    };
    let mut heti = vec![0i8; MAX_INS];
    let mut hetd = vec![0i8; MAX_INS];
    let rel_pos = pos - left;
    for cnum in 0..2 {
        let mut cons = Vec::with_capacity(span + 16);
        for i in 0..span {
            let ii = i as i64;
            if ii > rel_pos && out.cpos_pos == -1 {
                out.cpos_pos = cons.len() as i32;
            }
            let cb = cons_base[i];
            let (mut max_v, mut max_v2, mut max_j, mut max_j2, mut tot) = (0, 0, 4usize, 4usize, 0);
            for (j, &c) in cb.iter().enumerate() {
                if max_v < c {
                    max_v2 = max_v;
                    max_j2 = max_j;
                    max_v = c;
                    max_j = j;
                } else if max_v2 < c {
                    max_v2 = c;
                    max_j2 = j;
                }
                tot += c;
            }
            let ci = &mut cons_ins[i];
            if ty > 0
                && ii + left == pos + 1
                && ci.str.first().is_none_or(|s| (s.len() as i32) < ty)
            {
                if ci.str.is_empty() {
                    ci.str.push(Vec::new());
                    ci.freq.push(0);
                }
                ci.str[0].resize(ty as usize, 4);
            }
            let (mut max_v_ins, mut max_j_ins, mut tot_ins) = (0, 0usize, 0);
            for j in 0..ci.str.len() {
                if ci.freq[j] == 0 {
                    continue;
                }
                if max_v_ins < ci.freq[j] {
                    max_v_ins = ci.freq[j];
                    max_j_ins = j;
                }
                tot_ins += ci.freq[j];
            }
            let tot_sum = tot;
            let always_ins = (ii == rel_pos + 1 && ty > 0)
                || f64::from(max_v_ins) > CONS_CUTOFF_INC2 * f64::from(tot_sum);
            let mut het_ins = false;
            if !always_ins && max_v_ins as u32 >= bca.min_support {
                if cnum == 0 {
                    het_ins = f64::from(max_v_ins) > CONS_CUTOFF_INC * f64::from(tot_sum);
                    if i < MAX_INS {
                        heti[i] = if het_ins {
                            1
                        } else if f64::from(max_v_ins) > 0.3 * f64::from(tot_sum) {
                            -1
                        } else {
                            0
                        };
                    }
                } else {
                    het_ins = i < MAX_INS && heti[i] == -1;
                }
            }
            if always_ins || het_ins {
                let ins_len = ci.str.get(max_j_ins).map_or(0, Vec::len);
                if f64::from(max_v_ins) > CONS_CUTOFF_INS * f64::from(tot_ins) {
                    for j in 0..ins_len {
                        if cnum == 0 {
                            if (cons.len() as i64) < rel_pos + i64::from(out.left_shift) {
                                out.left_shift += 1;
                            } else {
                                out.right_shift += 1;
                            }
                        }
                        cons.push(ci.str[max_j_ins][j]);
                    }
                } else {
                    cons.extend(std::iter::repeat_n(4u8, ins_len));
                }
            }
            let always_del = (ty < 0 && ii > rel_pos && ii <= rel_pos - i64::from(ty))
                || f64::from(cb[5]) > CONS_CUTOFF2 * f64::from(tot);
            let mut het_del = false;
            if !always_del && cb[5] as u32 >= bca.min_support {
                if cnum == 0 {
                    let mut tot2 = tot;
                    if ii > rel_pos && ii <= rel_pos - i64::from(biggest_del) {
                        tot2 = total_span_str - type_depth;
                    }
                    het_del = f64::from(cb[5]) >= CONS_CUTOFF_DEL * f64::from(tot2);
                    if i < MAX_INS {
                        hetd[i] = if ii > rel_pos && ii <= rel_pos - i64::from(biggest_del) {
                            0
                        } else if het_del {
                            1
                        } else if f64::from(cb[5]) >= 0.3 * f64::from(tot2) {
                            -1
                        } else {
                            0
                        };
                    }
                } else {
                    het_del = i < MAX_INS && hetd[i] == -1;
                    if max_j == 5 && !het_del {
                        max_v = max_v2;
                        max_j = max_j2;
                    }
                }
            }
            if always_del || het_del {
                if (cons.len() as i64) < rel_pos + i64::from(out.left_shift) {
                    out.left_shift -= 1;
                } else {
                    out.right_shift += 1;
                }
            } else if f64::from(max_v) > CONS_CUTOFF * f64::from(tot) {
                cons.push(max_j as u8);
            } else if max_v > 0 {
                cons.push(4);
            } else {
                let k = cons.len() as i64;
                cons.push(if left + k < r.len() {
                    base6(r.at(left + k))
                } else {
                    4
                });
            }
        }
        out.cons[cnum] = cons;
    }
    out
}

/// edlib's `transformSequences` and the HW-mode `edlibAlign` with `EDLIB_TASK_LOC`:
/// the edit distance and the first end location's span in the target, or `None`.
fn edlib_align(query: &[u8], target: &[u8]) -> Option<(i32, i32, i32)> {
    let mut letter_idx = [0u8; 8];
    let mut in_alpha = [false; 8];
    let mut alpha_len = 0u8;
    let mut map = |c: u8| -> u8 {
        let c = c as usize & 7;
        if !in_alpha[c] {
            in_alpha[c] = true;
            letter_idx[c] = alpha_len;
            alpha_len += 1;
        }
        letter_idx[c]
    };
    let q: Vec<u8> = query.iter().map(|&c| map(c)).collect();
    let t: Vec<u8> = target.iter().map(|&c| map(c)).collect();
    let alpha = alpha_len as usize;
    if q.is_empty() || t.is_empty() {
        // edlib leaves the start locations unset, which edlib_glocal treats as a failure.
        return None;
    }
    let max_blocks = q.len().div_ceil(64);
    let w = (max_blocks * 64 - q.len()) as i32;
    let peq = build_peq(alpha, &q);
    let mut k: i32 = 64;
    let (best, positions) = loop {
        let (best, positions) =
            myers_semi_global(&peq, w, max_blocks, q.len() as i32, &t, k, Mode::Hw);
        if best != -1 {
            break (best, positions);
        }
        k *= 2;
    };
    let end = positions[0];
    let start = if end == -1 {
        0
    } else {
        let rt: Vec<u8> = t[..=end as usize].iter().rev().copied().collect();
        let rq: Vec<u8> = q.iter().rev().copied().collect();
        let rpeq = build_peq(alpha, &rq);
        let (_, pos_shw) =
            myers_semi_global(&rpeq, w, max_blocks, q.len() as i32, &rt, best, Mode::Shw);
        end - pos_shw[pos_shw.len() - 1]
    };
    Some((best, start, end))
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Hw,
    Shw,
}

/// `buildPeq`: the query profile, with a wildcard row last.
fn build_peq(alpha: usize, q: &[u8]) -> Vec<u64> {
    let max_blocks = q.len().div_ceil(64);
    let mut peq = vec![0u64; (alpha + 1) * max_blocks];
    for symbol in 0..alpha {
        for b in 0..max_blocks {
            let mut w = 0u64;
            for r in (b * 64..(b + 1) * 64).rev() {
                w = (w << 1) + u64::from(r >= q.len() || q[r] as usize == symbol);
            }
            peq[symbol * max_blocks + b] = w;
        }
    }
    for b in 0..max_blocks {
        peq[alpha * max_blocks + b] = u64::MAX;
    }
    peq
}

#[derive(Clone, Copy, Default)]
struct Block {
    p: u64,
    m: u64,
    score: i32,
}

/// `calculateBlock`
#[inline]
fn calculate_block(pv: u64, mv: u64, eq: u64, hin: i32) -> (u64, u64, i32) {
    let hin_is_neg = ((hin >> 2) as u64) & 1;
    let xv = eq | mv;
    let eq = eq | hin_is_neg;
    let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
    let mut ph = mv | !(xh | pv);
    let mut mh = pv & xh;
    let mut hout = ((ph & (1 << 63)) >> 63) as i32;
    hout -= ((mh & (1 << 63)) >> 63) as i32;
    ph <<= 1;
    mh <<= 1;
    mh |= hin_is_neg;
    ph |= ((hin + 1) >> 1) as u64;
    (mh | !(xv | ph), ph & xv, hout)
}

/// `getBlockCellValues`
fn block_cell_values(block: Block) -> [i32; 64] {
    let mut scores = [0i32; 64];
    let mut score = block.score;
    let mut mask = 1u64 << 63;
    for s in scores.iter_mut().take(63) {
        *s = score;
        if block.p & mask != 0 {
            score -= 1;
        }
        if block.m & mask != 0 {
            score += 1;
        }
        mask >>= 1;
    }
    scores[63] = score;
    scores
}

fn all_block_cells_larger(block: Block, k: i32) -> bool {
    block_cell_values(block).iter().all(|&s| s > k)
}

/// `myersCalcEditDistanceSemiGlobal`: the best score and up to 100 end positions.
fn myers_semi_global(
    peq: &[u64],
    w: i32,
    max_blocks: usize,
    query_len: i32,
    target: &[u8],
    k: i32,
    mode: Mode,
) -> (i32, Vec<i32>) {
    const STRONG_REDUCE_NUM: usize = 2048;
    const MAX_POS: usize = 100;
    let mut k = k;
    let mut first_block = 0i32;
    let mut last_block = (((k + 1 + 63) / 64) as usize).min(max_blocks) as i32 - 1;
    let mut blocks = vec![Block::default(); max_blocks];
    if mode == Mode::Hw {
        k = k.min(query_len);
    }
    for b in 0..=last_block.max(-1) {
        if b < 0 {
            break;
        }
        blocks[b as usize] = Block {
            p: u64::MAX,
            m: 0,
            score: (b + 1) * 64,
        };
    }
    let mut best = -1i32;
    let mut positions: Vec<i32> = Vec::new();
    let start_hout = if mode == Mode::Hw { 0 } else { 1 };
    for (c, &tc) in target.iter().enumerate() {
        let peq_c = &peq[tc as usize * max_blocks..];
        let mut hout = start_hout;
        for b in first_block..=last_block {
            let bl = &mut blocks[b as usize];
            let (p, m, h) = calculate_block(bl.p, bl.m, peq_c[b as usize], hout);
            bl.p = p;
            bl.m = m;
            hout = h;
            bl.score += hout;
        }
        // bl and Peq_c now point at last_block.
        let lb = last_block;
        if (lb < max_blocks as i32 - 1)
            && lb >= 0
            && (blocks[lb as usize].score - hout <= k)
            && ((peq_c[(lb + 1) as usize] & 1) != 0 || hout < 0)
        {
            last_block += 1;
            let nb = last_block as usize;
            let (p, m, h) = calculate_block(u64::MAX, 0, peq_c[nb], hout);
            blocks[nb] = Block {
                p,
                m,
                score: blocks[nb - 1].score - hout + 64 + h,
            };
        } else {
            while last_block >= first_block && blocks[last_block as usize].score >= k + 64 {
                last_block -= 1;
            }
        }
        if c % STRONG_REDUCE_NUM == 0 {
            while last_block >= 0
                && last_block >= first_block
                && all_block_cells_larger(blocks[last_block as usize], k)
            {
                last_block -= 1;
            }
        }
        if mode == Mode::Hw && last_block == -1 {
            last_block += 1;
        }
        if mode != Mode::Hw {
            while first_block <= last_block && blocks[first_block as usize].score >= k + 64 {
                first_block += 1;
            }
            if c % STRONG_REDUCE_NUM == 0 {
                while first_block <= last_block
                    && all_block_cells_larger(blocks[first_block as usize], k)
                {
                    first_block += 1;
                }
            }
        }
        if last_block < first_block {
            return (best, positions);
        }
        if last_block == max_blocks as i32 - 1 {
            let col_score = blocks[last_block as usize].score;
            if col_score <= k && (best == -1 || col_score <= best) {
                if col_score != best {
                    positions.clear();
                    best = col_score;
                    k = best;
                }
                if positions.len() < MAX_POS {
                    positions.push(c as i32 - w);
                }
            }
        }
    }
    if last_block == max_blocks as i32 - 1 {
        let scores = block_cell_values(blocks[last_block as usize]);
        for i in 0..w as usize {
            let col_score = scores[i + 1];
            if col_score <= k && (best == -1 || col_score <= best) {
                if col_score != best {
                    positions.clear();
                    best = col_score;
                    k = best;
                }
                if positions.len() < MAX_POS {
                    positions.push(target.len() as i32 - w + i as i32);
                }
            }
        }
    }
    (best, positions)
}

/// `edlib_glocal`: the edit distance skewed by the deletion bias, scaled by `m`.
fn edlib_glocal(reference: &[u8], query: &[u8], m: f64, del_bias: f64) -> i32 {
    match edlib_align(query, reference) {
        Some((dist, start, end)) => {
            let t_len = end - start + 1;
            c_int(m * (f64::from(dist) - del_bias * f64::from(t_len - query.len() as i32)))
        }
        None => i32::MAX,
    }
}

/// `bcf_cgp_align_score` of `bam2bcf_edlib.c`.
#[expect(clippy::too_many_arguments, reason = "a line-for-line port")]
fn align_score(
    bca: &CallAux,
    ty: i32,
    ref1: &[u8],
    ref2: &[u8],
    query: &[u8],
    mut tbeg: i32,
    mut tend1: i32,
    mut tend2: i32,
    left: i32,
    del_bias: f64,
) -> i32 {
    let atype = ty.abs();
    let at1 = |i: i32| ref1.get(i as usize).copied().unwrap_or(4);
    let at2 = |i: i32| ref2.get(i as usize).copied().unwrap_or(4);
    let mut l = 0;
    while l < tend1 - tbeg && l < tend2 - tbeg {
        if at1(l + tbeg - left) != 4 || at2(l + tbeg - left) != 4 {
            break;
        }
        l += 1;
    }
    if l > atype {
        tbeg += l - atype;
    }
    let mut l = tend1 - tbeg - 1;
    while l >= 0 {
        if at1(l + tbeg - left) != 4 {
            break;
        }
        l -= 1;
    }
    let l = tend1 - tbeg - 1 - l;
    if l > atype {
        tend1 -= l - atype;
    }
    let mut l = tend2 - tbeg - 1;
    while l >= 0 {
        if at2(l + tbeg - left) != 4 {
            break;
        }
        l -= 1;
    }
    let l = tend2 - tbeg - 1 - l;
    if l > atype {
        tend2 -= l - atype;
    }
    let slice = |r: &[u8], a: i32, n: i32| -> Vec<u8> {
        (0..n.max(0))
            .map(|i| r.get((a + i) as usize).copied().unwrap_or(4))
            .collect()
    };
    let mm = 30.0;
    let off = tbeg - left;
    let s2 = slice(ref2, off, tend2 - tbeg);
    let mut sc2 = edlib_glocal(&s2, query, mm, del_bias);
    let s1 = slice(ref1, off, tend1 - tbeg);
    let sc1 = if tend1 != tend2 || s1 != s2[..s1.len().min(s2.len())] {
        edlib_glocal(&s1, query, mm, del_bias)
    } else {
        i32::MAX
    };
    if sc1 < 0 && sc2 < 0 {
        return NO_SCORE;
    }
    // The best of the two, ignoring failed alignments.
    if sc1 >= 0 && (sc2 < 0 || sc2 > sc1) {
        sc2 = sc1;
    }
    let qlen = query.len() as i32;
    let l = c_int(0.5 * (100.0 * f64::from(sc2) / f64::from(qlen) + 0.499));
    (sc2 << 8) | c_int((f64::from(l as f32 * bca.indel_bias) * 0.5).min(255.0))
}

/// `est_seqQ` of `bam2bcf_edlib.c`
fn est_seq_q(bca: &CallAux, l: i32, l_run: i32) -> i32 {
    let q = bca.openq + bca.extq * (l.abs() - 1);
    let qh = (f64::from(bca.tandemq) * f64::from(l.abs()) / f64::from(l_run) + 0.499) as i32;
    q.min(qh)
}

/// `bcf_cgp_compute_indelQ` of `bam2bcf_edlib.c`.
#[expect(clippy::too_many_arguments, reason = "a line-for-line port")]
fn compute_indel_q(
    site: &mut Site<'_, '_>,
    bca: &mut CallAux,
    inscns: &[u8],
    l_run: i32,
    max_ins: i32,
    ref_type: usize,
    types: &[i32],
    qavg: f64,
    score: &[i32],
) -> i32 {
    let n_types = types.len();
    let mut sumq = [0i32; MAX_TYPES];
    let mut sc = [0i32; MAX_TYPES];
    let mut kk = 0usize;
    for s in site.plp.iter_mut() {
        for p in s.iter_mut() {
            let sct = &score[kk * n_types..(kk + 1) * n_types];
            kk += 1;
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
            let (mut indel_q, mut seq_q);
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
                let indel_q1 = (sc[t.min(MAX_TYPES - 1)] >> 14) - (sc[0] >> 14);
                let mut t = 1;
                while t < n_types {
                    if (sc[t] & 0x3f) as usize == ref_type {
                        t += 1;
                        continue;
                    }
                    break;
                }
                if t == n_types {
                    t -= 1;
                }
                let indel_q2 = (sc[t] >> 14) - (sc[0] >> 14);
                seq_q = est_seq_q(bca, types[(sc[0] & 0x3f) as usize], l_run);
                indel_q =
                    (bca.vs_ref * indel_q1 as f32 + (1.0 - bca.vs_ref) * indel_q2 as f32) as i32;
            }
            let div = f64::from(bca.indel_bias) * 0.5;
            indel_q = (f64::from(indel_q) / div) as i32;
            if bca.poly_mqual {
                let qpos = p.qpos;
                let l_qseq = p.l_qseq();
                let mut min_q = p.qual(qpos);
                let base_l = p.seqi(if qpos + 1 < l_qseq { qpos + 1 } else { qpos });
                let mut l = qpos;
                while l >= 0 {
                    if p.seqi(l) != base_l {
                        break;
                    }
                    min_q = min_q.min(p.qual(l));
                    l -= 1;
                }
                // bam_seqi past the end reads the next byte of the record, as in C; treat it
                // as no base.
                let base = if qpos + 1 < l_qseq {
                    p.seqi(qpos + 1)
                } else {
                    0xff
                };
                let mut l = qpos + 1;
                while l < l_qseq {
                    min_q = min_q.min(p.qual(l));
                    if p.seqi(l) != base {
                        break;
                    }
                    l += 1;
                }
                let mq = f64::from(min_q);
                seq_q = (f64::from(seq_q) + (qavg / 20.0).min(mq - qavg / 10.0)) as i32;
                indel_q = (f64::from(indel_q) + (qavg / 20.0).min(mq - qavg / 5.0)) as i32;
                seq_q = seq_q.max(0);
                indel_q = indel_q.max(0);
            }
            let tmp = sc[0] >> 6 & 0xff;
            const TMP_MAGIC: f64 = 255.0;
            indel_q = if f64::from(tmp) > TMP_MAGIC {
                0
            } else {
                ((1.0 - f64::from(tmp) / TMP_MAGIC) * f64::from(indel_q) + 0.499) as i32
            };
            indel_q = indel_q.min(255);
            if indel_q > seq_q {
                indel_q = seq_q;
            }
            let indel_q = indel_q.min(255);
            let seq_q = seq_q.min(255);
            p.aux = (sc[0] & 0x3f) << 16 | seq_q << 8 | indel_q;
            sumq[(sc[0] & 0x3f) as usize] += indel_q;
        }
    }
    finish_types(site, bca, inscns, max_ins, ref_type, types, &mut sumq)
}

/// `bcf_edlib_gap_prep`: whether to call an indel at the position.
pub fn gap_prep(site: &mut Site<'_, '_>, bca: &mut CallAux, buf: &mut EdlibBuf) -> bool {
    if !indel::has_indel(site.plp) {
        return false;
    }
    let (mut qsum, mut qcount) = (0f64, 0f64);
    for p in site.plp.iter().flatten() {
        let kstart = (p.qpos - 50).max(0);
        let kend = (p.qpos + 50).min(p.l_qseq());
        for k in kstart..kend {
            qsum += f64::from(p.qual(k));
            qcount += 1.0;
        }
    }
    let qavg = (qsum + 1.0) / (qcount + 1.0);

    let Some(found) = find_types(site, bca) else {
        return false;
    };
    let FoundTypes {
        mut types,
        ref_type,
        max_rd_len,
        n_reads,
    } = found;
    let n_types = types.len();
    let pos = site.pos;
    let r = site.reference;
    let win = bca.indel_win_size;
    let mut max_indel = 20 * types[0].abs().max(types[n_types - 1].abs()) + win / 4;
    if max_indel > win {
        max_indel = win;
    }
    let max_indel = i64::from(max_indel);
    let left = if pos > max_indel { pos - max_indel } else { 0 };
    let mut right = pos + max_indel;
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

    let max_ins = types[n_types - 1];
    let mut l_run = l_run(r, pos);
    let l_run_base = i32::from(SEQ_NT16_TABLE[r.at(pos + 1) as usize]);
    let mut l_run_ins = 0i32;
    let inscns = if max_ins > 0 {
        calc_cons(site.plp, &mut types, max_ins)
    } else {
        Vec::new()
    };
    let mi = max_ins.max(0) as usize;

    buf.query.clear();
    buf.query.resize(
        (right - left + i64::from(max_rd_len) + i64::from(max_ins) + 2) as usize,
        0,
    );
    buf.score.clear();
    buf.score.resize(n_reads * n_types, 0);
    bca.indelreg = 0;
    let biggest_del = types.iter().copied().min().unwrap_or(0).min(0);
    let biggest_ins = types.iter().copied().max().unwrap_or(0).max(0);
    let mut band = biggest_ins - biggest_del;

    let (mut pos_l, mut pos_r) = (pos, pos);
    {
        let pstart = (pos - 30).max(0);
        let pmid = pos - pstart;
        let pend = r.len().min(pos + 30);
        let slice: Vec<u8> = (pstart..pend).map(|i| r.at(i)).collect();
        for (start, end, _) in find_str(&slice) {
            let (start, end) = (i64::from(start), i64::from(end));
            if end >= pmid && start <= pmid {
                if pos_l > pstart + start {
                    pos_l = pstart + start;
                }
                if pos_r < pstart + end {
                    pos_r = pstart + end;
                }
            }
        }
    }

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
        let mut kread = 0usize;
        for s in 0..site.plp.len() {
            let tc = consensus(
                &site.plp[s],
                pos,
                bca,
                r,
                left,
                right,
                types[t],
                biggest_del,
                &mut band,
                pos_l,
                pos_r,
            );
            let cpos = tc.cpos_pos.max(0) as usize;
            let c0 = &tc.cons[0];
            let k = c0.get(cpos).copied().unwrap_or(0);
            let mut j = 0;
            while j < types[t] {
                if c0.get(cpos + j as usize).copied().unwrap_or(0xff) != k {
                    break;
                }
                j += 1;
            }
            if j != 0 && j == types[t] {
                l_run_ins |= i32::from([1u8, 2, 4, 8, 0xf, 0][(k as usize).min(5)]);
            }
            if types[t] < 0 {
                l_run_ins |= 0xff;
            }
            let tot_str: i32 = find_str(c0).iter().map(|&(s, e, _)| e - s).sum();
            for i in 0..site.plp[s].len() {
                let p = site.plp[s][i];
                let slot = kread * n_types + t;
                kread += 1;
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
                let (mut left2, mut right2) = (left, right);
                let mut min_win_size = (-biggest_del).max(biggest_ins);
                min_win_size += tc.left_shift.abs() + tc.right_shift.abs();
                min_win_size += tot_str;
                min_win_size += 10;
                let mws = i64::from(min_win_size);
                if p.l_qseq() > LONG_READ_LEN {
                    if pos - left >= mws {
                        left2 = left2.max(pos - mws);
                    }
                    if right - pos >= mws {
                        right2 = right2.min(pos + mws);
                    }
                }
                let bpos = p.b.position() as i32;
                let _r_end = cigar2rlen(cigar) - 1 + bpos;
                let (qbeg, tbeg0) = tpos2qpos(bpos, cigar, left2 as i32, DelEnd::Right);
                let (qend, tend0) = tpos2qpos(bpos, cigar, right2 as i32, DelEnd::Left);
                for l in qbeg..qend {
                    buf.query[(l - qbeg) as usize] = SEQ_NT16_INT[p.seqi(l) as usize];
                }
                let wband = band + (-biggest_del).max(biggest_ins) * 2 + 20;
                let off = (left2 - left) as i32;
                let mut tend1 = left as i32 + tc.cons[0].len() as i32 - off;
                let mut tend2 = left as i32 + tc.cons[1].len() as i32 - off;
                tend1 = tend1.min(tend0 + wband);
                tend2 = tend2.min(tend0 + wband);
                let tbeg = (left2 as i32).max(tbeg0 - wband);
                if tend1 > tbeg && tend2 > tbeg {
                    let r1 = tc.cons[0].get(off as usize..).unwrap_or(&[]);
                    let r2 = tc.cons[1].get(off as usize..).unwrap_or(&[]);
                    buf.score[slot] = align_score(
                        bca,
                        types[t],
                        r1,
                        r2,
                        &buf.query[..(qend - qbeg).max(0) as usize],
                        tbeg,
                        tend1,
                        tend2,
                        left2 as i32,
                        f64::from(bca.del_bias),
                    );
                } else {
                    buf.score[slot] = NO_SCORE;
                }
            }
        }
    }
    if l_run_base & l_run_ins == 0 {
        l_run = 1;
    }
    let score = std::mem::take(&mut buf.score);
    let n_alt = compute_indel_q(
        site, bca, &inscns, l_run, max_ins, ref_type, &types, qavg, &score,
    );
    buf.score = score;
    let _ = CMATCH;
    n_alt > 0
}
