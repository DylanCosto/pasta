//! Genotype likelihoods and annotations from a pileup column: a port of bcftools'
//! `bam2bcf.c` (`bcf_call_glfgen`, `bcf_call_combine`, `bcf_call2bcf`) with the same
//! floating-point operations, so PL, I16, QS and the bias tests match bcftools.

use pasta::{
    pileup::Read,
    variant::{Header, Record},
};

use super::{
    budget::ThreadBudget,
    errmod::{ERROR_MODEL, Rng48, kf_erfc, kt_fisher_exact},
};

// The annotation flags (`B2B_*`).
pub const FMT_DP: u32 = 1 << 0;
pub const FMT_SP: u32 = 1 << 1;
pub const FMT_DV: u32 = 1 << 2;
pub const FMT_DP4: u32 = 1 << 3;
pub const FMT_DPR: u32 = 1 << 4;
pub const INFO_DPR: u32 = 1 << 5;
pub const FMT_AD: u32 = 1 << 6;
pub const FMT_ADF: u32 = 1 << 7;
pub const FMT_ADR: u32 = 1 << 8;
pub const INFO_AD: u32 = 1 << 9;
pub const INFO_ADF: u32 = 1 << 10;
pub const INFO_ADR: u32 = 1 << 11;
pub const INFO_SCR: u32 = 1 << 12;
pub const FMT_SCR: u32 = 1 << 13;
pub const INFO_VDB: u32 = 1 << 14;
pub const FMT_QS: u32 = 1 << 15;
pub const FMT_NMBZ: u32 = 1 << 16;
pub const INFO_NMBZ: u32 = 1 << 17;
pub const INFO_BQBZ: u32 = 1 << 18;
pub const INFO_MQBZ: u32 = 1 << 19;
pub const INFO_MQSBZ: u32 = 1 << 20;
pub const INFO_RPBZ: u32 = 1 << 21;
pub const INFO_SCBZ: u32 = 1 << 22;
pub const INFO_SGB: u32 = 1 << 23;
pub const FMT_QM: u32 = 1 << 24;
pub const INFO_NM: u32 = 1 << 25;
pub const INFO_MQ0F: u32 = 1 << 26;
pub const INFO_IDV: u32 = 1 << 27;
pub const INFO_IMF: u32 = 1 << 28;
pub const INFO_FS: u32 = 1 << 29;

/// The flags that need the allelic depth arrays.
pub const NEEDS_AD: u32 =
    INFO_DPR | FMT_DPR | INFO_AD | INFO_ADF | INFO_ADR | FMT_AD | FMT_ADF | FMT_ADR;

/// `B2B_INDEL_NULL`
pub const INDEL_NULL: i32 = 10000;
/// `B2B_MAX_ALLELES`
pub const MAX_ALLELES: usize = 5;
/// `B2B_N_NM`: the number of NM bins.
pub const N_NM: usize = 32;
/// `CAP_DIST`
const CAP_DIST: i32 = 25;
/// `DEF_MAPQ`: the mapping quality used for 255.
const DEF_MAPQ: i32 = 20;
/// `bca->npos`
pub const NPOS: usize = 100;
/// `bca->nqual`
pub const NQUAL: usize = 60;

/// What to do with ambiguous indel reads (`--ambig-reads`, `B2B_DROP` and so on).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AmbigReads {
    Drop,
    IncAd,
    IncAd0,
}

/// The `seq_nt16_int` table: 4-bit codes to 0-3 for A, C, G, T and 4 otherwise.
pub const SEQ_NT16_INT: [u8; 16] = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];

/// `seq_nt16_table`: ASCII to 4-bit codes.
pub const SEQ_NT16_TABLE: [u8; 256] = pasta::pileup::baq::SEQ_NT16_TABLE;

/// Per-read data cached while the read is piled up (bcftools' `plp_cd_t`), packed in the
/// read's user data: bit 0 soft clip, bit 1 indel, bit 2 realigned, bit 3 whether the read may
/// have a ZQ tag, bits 4-31 the sample, bits 32-63 the NM value (`PLP_NM`).
pub mod cd {
    pub const SOFT_CLIP: u64 = 1;
    pub const INDEL: u64 = 2;
    pub const REALN: u64 = 4;
    /// The read had a ZQ tag when it was read, or BAQ was applied to it.
    pub const HAS_ZQ: u64 = 8;

    /// `PLP_NM_UNSET` is never stored: NM is computed when the read is read.
    #[inline]
    pub fn sample(cd: u64) -> usize {
        ((cd & 0xffff_ffff) >> 4) as usize
    }

    /// The read's ZQ tag value (`bam_aux_get(b, "ZQ")` after the type), if it may have one.
    /// It runs to the end of the record, as the C pointer does: BAQ offsets of qualities above
    /// 191 wrap, and one may be a NUL.
    #[inline]
    pub fn zq<'a>(p: &super::Plp<'a>) -> Option<&'a [u8]> {
        if p.cd & HAS_ZQ == 0 {
            return None;
        }
        p.b.tag_value_tail(*b"ZQ")
    }

    #[inline]
    pub fn nm(cd: u64) -> i32 {
        (cd >> 32) as u32 as i32
    }

    /// Packs the sample, the flag bits (`SOFT_CLIP`, `INDEL`, `HAS_ZQ`), and NM.
    pub fn pack(sample: usize, flags: u64, nm: i32) -> u64 {
        debug_assert_eq!(flags & !(SOFT_CLIP | INDEL | HAS_ZQ), 0, "flags {flags:#x}");
        (sample as u64) << 4 | flags | (u64::from(nm as u32) << 32)
    }
}

/// `bam_pileup1_t` as bcftools sees it: one read at the current position, with the indel
/// caller's `aux`.
#[derive(Clone, Copy, Debug)]
pub struct Plp<'a> {
    pub b: &'a Read,
    pub qpos: i32,
    pub indel: i32,
    pub is_del: bool,
    pub is_refskip: bool,
    pub cd: u64,
    pub aux: i32,
}

impl Plp<'_> {
    #[inline]
    pub fn is_rev(&self) -> bool {
        self.b.flags() & 0x10 != 0
    }

    #[inline]
    pub fn l_qseq(&self) -> i32 {
        self.b.sequence_len() as i32
    }

    #[inline]
    pub fn qual(&self, i: i32) -> i32 {
        i32::from(self.b.quality_scores()[i as usize])
    }

    #[inline]
    pub fn seqi(&self, i: i32) -> u8 {
        self.b.base(i as usize)
    }

    #[inline]
    pub fn mapq(&self) -> i32 {
        i32::from(self.b.mapping_quality())
    }
}

/// `bcf_callaux_t`
pub struct CallAux {
    pub fmt_flag: u32,
    pub ambig_reads: AmbigReads,
    pub capq: i32,
    pub min_baseq: i32,
    pub max_baseq: i32,
    pub delta_baseq: i32,
    pub openq: i32,
    pub extq: i32,
    pub tandemq: i32,
    pub min_support: u32,
    pub max_support: u32,
    pub min_frac: f64,
    pub max_frac: f32,
    pub per_sample_flt: bool,
    pub ref_pos: [i32; NPOS],
    pub alt_pos: [i32; NPOS],
    pub ref_mq: [i32; NQUAL],
    pub alt_mq: [i32; NQUAL],
    pub ref_bq: [i32; NQUAL],
    pub alt_bq: [i32; NQUAL],
    pub fwd_mqs: [i32; NQUAL],
    pub rev_mqs: [i32; NQUAL],
    pub iref_pos: [i32; NPOS],
    pub ialt_pos: [i32; NPOS],
    pub iref_mq: [i32; NQUAL],
    pub ialt_mq: [i32; NQUAL],
    pub ref_scl: [i32; 100],
    pub alt_scl: [i32; 100],
    pub iref_scl: [i32; 100],
    pub ialt_scl: [i32; 100],
    pub indel_types: [i32; 4],
    pub indel_win_size: i32,
    pub indels_v20: bool,
    pub edlib: bool,
    pub seqq_offset: i32,
    pub maxins: i32,
    pub indelreg: i32,
    pub poly_mqual: bool,
    /// `inscns`: the inserted sequences of the indel types, `maxins` each.
    pub inscns: Vec<u8>,
    pub bases: Vec<u16>,
    pub indel_bias: f32,
    pub del_bias: f32,
    pub vs_ref: f32,
    /// The threads the realignments of a heavy position may borrow (`--threads`).
    pub threads: ThreadBudget,
    /// The site-wide NM histograms (`bca->ref_nm`, `bca->alt_nm`).
    pub ref_nm: [i32; N_NM],
    pub alt_nm: [i32; N_NM],
    pub nnm: [u32; 2],
    pub nm: [f32; 2],
    /// The `drand48` state errmod shuffles deep pileups with.
    pub rng: Rng48,
    /// How many numbers have been drawn from `rng`.
    pub draws: u64,
}

/// Annotation and model settings, as `mpileup` sets them on `bcf_callaux_t`.
#[derive(Clone, Debug)]
pub struct CallSettings {
    pub fmt_flag: u32,
    pub ambig_reads: AmbigReads,
    pub min_baseq: i32,
    pub max_baseq: i32,
    pub delta_baseq: i32,
    pub openq: i32,
    pub extq: i32,
    pub tandemq: i32,
    pub min_support: u32,
    pub min_frac: f64,
    pub per_sample_flt: bool,
    pub indel_win_size: i32,
    pub indels_v20: bool,
    pub edlib: bool,
    pub seqq_offset: i32,
    pub poly_mqual: bool,
    pub indel_bias: f32,
    pub del_bias: f32,
    pub vs_ref: f32,
    pub threads: ThreadBudget,
}

impl CallAux {
    /// `bcf_call_init` plus the settings mpileup copies in.
    pub fn new(s: &CallSettings, rng: Rng48) -> Self {
        Self {
            fmt_flag: s.fmt_flag,
            ambig_reads: s.ambig_reads,
            capq: 60,
            min_baseq: s.min_baseq,
            max_baseq: s.max_baseq,
            delta_baseq: s.delta_baseq,
            openq: s.openq,
            extq: s.extq,
            tandemq: s.tandemq,
            min_support: s.min_support,
            max_support: 0,
            min_frac: s.min_frac,
            max_frac: 0.0,
            per_sample_flt: s.per_sample_flt,
            ref_pos: [0; NPOS],
            alt_pos: [0; NPOS],
            ref_mq: [0; NQUAL],
            alt_mq: [0; NQUAL],
            ref_bq: [0; NQUAL],
            alt_bq: [0; NQUAL],
            fwd_mqs: [0; NQUAL],
            rev_mqs: [0; NQUAL],
            iref_pos: [0; NPOS],
            ialt_pos: [0; NPOS],
            iref_mq: [0; NQUAL],
            ialt_mq: [0; NQUAL],
            ref_scl: [0; 100],
            alt_scl: [0; 100],
            iref_scl: [0; 100],
            ialt_scl: [0; 100],
            indel_types: [0; 4],
            indel_win_size: s.indel_win_size,
            indels_v20: s.indels_v20,
            edlib: s.edlib,
            seqq_offset: s.seqq_offset,
            maxins: 0,
            indelreg: 0,
            poly_mqual: s.poly_mqual,
            inscns: Vec::new(),
            bases: Vec::new(),
            indel_bias: s.indel_bias,
            del_bias: s.del_bias,
            vs_ref: s.vs_ref,
            threads: s.threads.clone(),
            ref_nm: [0; N_NM],
            alt_nm: [0; N_NM],
            nnm: [0; 2],
            nm: [0.0; 2],
            rng,
            draws: 0,
        }
    }
}

/// `bcf_callret1_t`: one sample's results.
#[derive(Clone, Debug)]
pub struct CallRet {
    pub ori_depth: u32,
    pub mq0: u32,
    pub qm: [f64; MAX_ALLELES],
    pub adf: [i32; MAX_ALLELES],
    pub adr: [i32; MAX_ALLELES],
    pub scr: i32,
    pub qs: [i32; MAX_ALLELES],
    pub ref_nm: [i32; N_NM],
    pub alt_nm: [i32; N_NM],
    pub anno: [f64; 16],
    pub p: [f32; 25],
}

impl Default for CallRet {
    fn default() -> Self {
        Self {
            ori_depth: 0,
            mq0: 0,
            qm: [0.0; MAX_ALLELES],
            adf: [0; MAX_ALLELES],
            adr: [0; MAX_ALLELES],
            scr: 0,
            qs: [0; MAX_ALLELES],
            ref_nm: [0; N_NM],
            alt_nm: [0; N_NM],
            anno: [0.0; 16],
            p: [0.0; 25],
        }
    }
}

/// `bcf_call_t`: the site's results.
#[derive(Debug)]
pub struct Call {
    pub tid: i32,
    pub pos: i64,
    pub a: [i32; 5],
    pub qsum: [f32; MAX_ALLELES],
    pub n: usize,
    pub n_alleles: usize,
    pub ori_ref: i32,
    pub unseen: i32,
    pub anno: [f64; 16],
    pub depth: u32,
    pub ori_depth: u32,
    pub mq0: u32,
    pub pl: Vec<i32>,
    pub dp4: Vec<i32>,
    /// Totals in the first `MAX_ALLELES`, then per sample (compacted to `n_alleles` each by
    /// `combine`).
    pub adr: Vec<i32>,
    pub adf: Vec<i32>,
    pub scr: Vec<i32>,
    pub qs: Vec<i32>,
    pub qm: Vec<i32>,
    pub mwu_nm: Vec<f32>,
    pub vdb: f32,
    pub mwu_pos: f32,
    pub mwu_mq: f32,
    pub mwu_bq: f32,
    pub mwu_mqs: f32,
    pub mwu_sc: f32,
    pub nm: [f32; 2],
    pub seg_bias: f32,
    pub strand_bias: f32,
    /// Whether the last `combine` found an alternate allele (and computed the tests).
    pub has_alt: bool,
    /// Scratch for the allele string.
    tmp: Vec<u8>,
    itmp: Vec<i32>,
}

impl Call {
    pub fn new(n: usize) -> Self {
        Self {
            tid: 0,
            pos: 0,
            a: [0; 5],
            qsum: [0.0; MAX_ALLELES],
            n,
            n_alleles: 0,
            ori_ref: 0,
            unseen: 0,
            anno: [0.0; 16],
            depth: 0,
            ori_depth: 0,
            mq0: 0,
            pl: vec![0; 15 * n],
            dp4: vec![0; 4 * n],
            adr: vec![0; (n + 1) * MAX_ALLELES],
            adf: vec![0; (n + 1) * MAX_ALLELES],
            scr: vec![0; n + 1],
            qs: vec![0; n * MAX_ALLELES],
            qm: vec![0; n * MAX_ALLELES],
            mwu_nm: vec![0.0; n + 1],
            vdb: 0.0,
            mwu_pos: 0.0,
            mwu_mq: 0.0,
            mwu_bq: 0.0,
            mwu_mqs: 0.0,
            mwu_sc: 0.0,
            nm: [0.0; 2],
            seg_bias: 0.0,
            strand_bias: 0.0,
            has_alt: false,
            tmp: Vec::new(),
            itmp: Vec::new(),
        }
    }
}

/// `bcf_callaux_clean`
pub fn callaux_clean(bca: &mut CallAux, call: &mut Call, bcr: &mut [CallRet]) {
    bca.ref_pos.fill(0);
    bca.alt_pos.fill(0);
    bca.iref_pos.fill(0);
    bca.ialt_pos.fill(0);
    bca.ref_mq.fill(0);
    bca.alt_mq.fill(0);
    bca.iref_mq.fill(0);
    bca.ialt_mq.fill(0);
    bca.ref_bq.fill(0);
    bca.alt_bq.fill(0);
    bca.fwd_mqs.fill(0);
    bca.rev_mqs.fill(0);
    call.adf.fill(0);
    call.adr.fill(0);
    call.scr.fill(0);
    bca.ref_nm.fill(0);
    bca.alt_nm.fill(0);
    call.qm.fill(0);
    call.qs.fill(0);
    for r in bcr.iter_mut() {
        r.adf = [0; MAX_ALLELES];
        r.adr = [0; MAX_ALLELES];
        r.qs = [0; MAX_ALLELES];
        r.ref_nm = [0; N_NM];
        r.alt_nm = [0; N_NM];
    }
    bca.ref_scl.fill(0);
    bca.alt_scl.fill(0);
    bca.iref_scl.fill(0);
    bca.ialt_scl.fill(0);
    bca.nnm = [0; 2];
    bca.nm = [0.0; 2];
}

/// `get_aux_nm`: the read's NM (cached when it was read), less the base itself, binned.
fn get_aux_nm(p: &Plp<'_>, is_ref: bool) -> i32 {
    let nm = cd::nm(p.cd);
    if nm == -1 {
        return -1;
    }
    let nm = nm - if is_ref { 1 } else { 2 };
    nm.clamp(0, N_NM as i32 - 1)
}

/// The NM value bcftools' `get_aux_nm` derives from a read: the NM tag with indels counted
/// once and soft clips as mismatches; -1 without an NM tag.
pub fn read_nm(b: &Read) -> i32 {
    let Some(nm) = b.integer_tag(*b"NM") else {
        return -1;
    };
    let mut nm = nm;
    for &c in b.cigar() {
        let op = c & 0xf;
        let len = i64::from(c >> 4);
        if op == 4 {
            nm += len;
        } else if (op == 1 || op == 2) && len > 1 {
            nm -= len - 1;
        }
    }
    nm as i32
}

/// `get_position`: the position in the aligned part of the read, and the length and distance
/// of the nearest soft clip.
fn get_position(p: &Plp<'_>) -> (i32, i32, i32, i32) {
    let cigar = p.b.cigar();
    let mut edist = p.qpos + 1;
    let (mut sc_left, mut sc_right) = (0i32, 0i32);
    let (mut sc_left_dist, mut sc_right_dist) = (-1i32, -1i32);
    let mut i = 0;
    while i < cigar.len() {
        let op = cigar[i] & 0xf;
        if op == 5 {
            i += 1;
            continue;
        } else if op == 4 {
            sc_left += (cigar[i] >> 4) as i32;
        } else {
            break;
        }
        i += 1;
    }
    if sc_left != 0 {
        sc_left_dist = p.qpos + 1 - sc_left;
    }
    edist -= sc_left;
    let mut j = cigar.len() as isize - 1;
    while j >= i as isize {
        let op = cigar[j as usize] & 0xf;
        if op == 5 {
            j -= 1;
            continue;
        } else if op == 4 {
            sc_right += (cigar[j as usize] >> 4) as i32;
        } else {
            break;
        }
        j -= 1;
    }
    if sc_right != 0 {
        sc_right_dist = p.l_qseq() - sc_right - p.qpos;
    }
    let (mut sc_len, mut sc_dist) = (0, 0);
    if sc_left_dist >= 0 {
        if sc_right_dist < 0 || sc_left_dist < sc_right_dist {
            sc_len = sc_left;
            sc_dist = sc_left_dist;
        }
    } else if sc_right_dist >= 0 {
        sc_len = sc_right;
        sc_dist = sc_right_dist;
    }
    let len = p.l_qseq() - sc_left - sc_right;
    (edist, len, sc_len, sc_dist)
}

/// C's conversion of a double to `int` as x86 performs it (`cvttsd2si`): truncation, with
/// NaN and values out of range giving `INT_MIN`. Where bcftools converts values that can be
/// out of range (scores of failed alignments), its output depends on this.
#[inline]
pub fn c_int(x: f64) -> i32 {
    if x.is_nan() || !(-2_147_483_648.0..2_147_483_648.0).contains(&x) {
        i32::MIN
    } else {
        x as i32
    }
}

/// `qual2err` in bcftools.h: `pow(10, -0.1*qual)`.
fn qual2err(q: i32) -> f64 {
    10f64.powf(-0.1 * f64::from(q))
}

/// `bcf_call_glfgen`: one sample's likelihoods and annotations. `ref_base` is the 4-bit
/// reference base, or `None` when testing an indel. Returns the number of bases used, or -1
/// for an empty pileup.
pub fn glfgen(plps: &[Plp<'_>], ref_base: Option<u8>, bca: &mut CallAux, r: &mut CallRet) -> i32 {
    r.ori_depth = 0;
    r.mq0 = 0;
    r.anno = [0.0; 16];
    r.p = [0.0; 25];
    r.qm = [0.0; MAX_ALLELES];
    r.scr = 0;

    let (ref4, is_indel) = match ref_base {
        Some(b) => (i32::from(SEQ_NT16_INT[b as usize]), false),
        None => (4, true),
    };
    let n_reads = plps.len() as i32;
    if plps.is_empty() {
        return -1;
    }
    let want_ad = bca.fmt_flag & NEEDS_AD != 0;
    let want_nm = bca.fmt_flag & (FMT_NMBZ | INFO_NMBZ | INFO_NM) != 0;
    let want_sample_nm = bca.fmt_flag & FMT_NMBZ != 0;

    let indel_in_sample = bca.edlib && plps.iter().any(|p| p.indel != 0);

    let nqual_over_60 = NQUAL as f64 / 60.0;
    let mut adr_ref_missed = [0i32; 4];
    let mut adf_ref_missed = [0i32; 4];
    let mut ori_depth = 0u32;
    bca.bases.clear();
    for p in plps {
        if bca.fmt_flag & (INFO_SCR | FMT_SCR) != 0 && p.cd & cd::SOFT_CLIP != 0 {
            r.scr += 1;
        }
        if p.is_refskip || p.b.flags() & 0x4 != 0 {
            continue;
        }
        if p.is_del && !is_indel {
            continue;
        }
        let mut inm = -1;
        ori_depth += 1;
        let b: i32;
        let mut q: i32;
        let mut seq_q: i32;
        let base_q: i32;
        let is_diff: bool;
        if is_indel {
            let mut bb = p.aux >> 16 & 0x3f;
            q = p.aux & 0xff;
            seq_q = q;
            if bca.edlib && !indel_in_sample && p.aux & 0xff != 0 {
                q = p.qual(p.qpos).min(bca.max_baseq);
                seq_q = 99;
            }
            if !bca.indels_v20 && !bca.edlib {
                if p.indel == 0 && (q < n_reads / 2 || n_reads > 20) {
                    bb = 0;
                    q = p.qual(p.qpos);
                    seq_q = (3 * seq_q + 2 * q) / 8;
                }
                if n_reads > 20 && seq_q > 40 {
                    seq_q = 40;
                }
            }
            b = bb;
            is_diff = b != 0;
            if want_nm {
                inm = get_aux_nm(p, !is_diff);
                if inm >= 0 {
                    bca.nnm[usize::from(is_diff)] += 1;
                    bca.nm[usize::from(is_diff)] += inm as f32;
                }
            }
            if q < bca.min_baseq {
                if p.indel == 0 && b < 4 {
                    if p.is_rev() {
                        adr_ref_missed[b as usize] += 1;
                    } else {
                        adf_ref_missed[b as usize] += 1;
                    }
                }
                continue;
            }
            if bca.edlib {
                seq_q = seq_q.min(bca.seqq_offset - n_reads.min(20) * 5);
                if indel_in_sample && p.indel == 0 && b != 0 {
                    seq_q = seq_q.min(seq_q / 2 + 5);
                    q = (p.qual(p.qpos) / 4 + 10).min(q / 4 + 1);
                }
            }
            base_q = p.aux >> 8 & 0xff;
        } else {
            let mut bb = i32::from(p.seqi(p.qpos));
            bb = i32::from(
                SEQ_NT16_INT[if bb != 0 {
                    bb
                } else {
                    i32::from(ref_base.unwrap_or(0))
                } as usize],
            );
            b = bb;
            q = p.qual(p.qpos);
            if p.qpos > 0 && q > p.qual(p.qpos - 1) + bca.delta_baseq {
                q = p.qual(p.qpos - 1) + bca.delta_baseq;
            }
            if p.qpos + 1 < p.l_qseq() && q > p.qual(p.qpos + 1) + bca.delta_baseq {
                q = p.qual(p.qpos + 1) + bca.delta_baseq;
            }
            if q < bca.min_baseq {
                continue;
            }
            if q > bca.max_baseq {
                q = bca.max_baseq;
            }
            base_q = q;
            seq_q = 99;
            is_diff = !(ref4 < 4 && b == ref4);
            if want_nm {
                inm = get_aux_nm(p, !is_diff);
                if inm >= 0 {
                    bca.nnm[usize::from(is_diff)] += 1;
                    bca.nm[usize::from(is_diff)] += inm as f32;
                }
            }
        }
        let mut map_q = if p.mapq() < 255 { p.mapq() } else { DEF_MAPQ };
        if map_q == 0 {
            r.mq0 += 1;
        }
        if q > seq_q {
            q = seq_q;
        }
        map_q = map_q.min(bca.capq);
        if q > map_q {
            q = map_q;
        }
        q = q.clamp(4, 63);
        let is_rev = p.is_rev();
        bca.bases
            .push(((q << 5) as u16) | (u16::from(is_rev) << 4) | (b as u16));

        if b < 4 {
            let bi = b as usize;
            if bca.fmt_flag & FMT_QM != 0 {
                r.qm[bi] += qual2err(q);
            }
            r.qs[bi] += q;
            if want_ad {
                if is_rev {
                    r.adr[bi] += 1;
                } else {
                    r.adf[bi] += 1;
                }
            }
        }
        let d = usize::from(is_diff) << 1;
        r.anno[d | usize::from(is_rev)] += 1.0;
        let mut min_dist = p.l_qseq() - 1 - p.qpos;
        if min_dist > p.qpos {
            min_dist = p.qpos;
        }
        if min_dist > CAP_DIST {
            min_dist = CAP_DIST;
        }
        r.anno[4 | d] += f64::from(base_q);
        r.anno[4 | d | 1] += f64::from(base_q * base_q);
        r.anno[8 | d] += f64::from(map_q);
        r.anno[8 | d | 1] += f64::from(map_q * map_q);
        r.anno[12 | d] += f64::from(min_dist);
        r.anno[12 | d | 1] += f64::from(min_dist * min_dist);

        let base_q = base_q.min(59);
        let map_q = map_q.min(59);
        let mut epos = 0usize;
        let mut sc_len = 0i32;
        if bca.fmt_flag & (INFO_RPBZ | INFO_VDB | INFO_SCBZ) != 0 {
            let (pos, len, scl, sc_dist) = get_position(p);
            epos = (f64::from(pos) / f64::from(len + 1) * (NPOS - 1) as f64) as usize;
            if scl != 0 {
                sc_len = (15.0 * f64::from(scl) / f64::from(sc_dist + 1)) as i32;
                if sc_len > 99 {
                    sc_len = 99;
                }
            }
        }
        let imq = (f64::from(map_q) * nqual_over_60) as usize;
        let ibq = (f64::from(base_q) * nqual_over_60) as usize;
        if is_rev {
            bca.rev_mqs[imq] += 1;
        } else {
            bca.fwd_mqs[imq] += 1;
        }
        if !is_diff {
            bca.ref_pos[epos] += 1;
            bca.ref_bq[ibq] += 1;
            bca.ref_mq[imq] += 1;
            bca.ref_scl[sc_len as usize] += 1;
            if inm >= 0 {
                bca.ref_nm[inm as usize] += 1;
                if want_sample_nm {
                    r.ref_nm[inm as usize] += 1;
                }
            }
        } else {
            bca.alt_pos[epos] += 1;
            bca.alt_bq[ibq] += 1;
            bca.alt_mq[imq] += 1;
            bca.alt_scl[sc_len as usize] += 1;
            if inm >= 0 {
                bca.alt_nm[inm as usize] += 1;
                if want_sample_nm {
                    r.alt_nm[inm as usize] += 1;
                }
            }
        }
    }

    if want_ad && bca.ambig_reads == AmbigReads::IncAd0 {
        for i in 0..4 {
            r.adr[0] += adr_ref_missed[i];
            r.adf[0] += adf_ref_missed[i];
        }
    } else if want_ad && bca.ambig_reads == AmbigReads::IncAd {
        let dp: i32 = r.adr[..4].iter().sum();
        let dp_ambig: i32 = adr_ref_missed.iter().sum();
        if dp != 0 {
            for i in 0..4 {
                r.adr[i] += (dp_ambig as f32 * r.adr[i] as f32 / dp as f32).round() as i32;
            }
        }
        let dp: i32 = r.adf[..4].iter().sum();
        let dp_ambig: i32 = adf_ref_missed.iter().sum();
        if dp != 0 {
            for i in 0..4 {
                r.adf[i] += (dp_ambig as f32 * r.adf[i] as f32 / dp as f32).round() as i32;
            }
        }
    }

    r.ori_depth = ori_depth;
    let n = bca.bases.len();
    if n > 255 {
        bca.draws += n as u64 - 1;
    }
    ERROR_MODEL.cal(&mut bca.bases, 5, &mut r.p, &mut bca.rng);
    n as i32
}

/// `calc_vdb`: variant distance bias, `INFINITY` with fewer than two reads.
fn calc_vdb(pos: &[i32; NPOS]) -> f64 {
    const READLEN: i32 = 100;
    const PARAM: [[f32; 3]; 15] = [
        [3.0, 0.079, 18.0],
        [4.0, 0.09, 19.8],
        [5.0, 0.1, 20.5],
        [6.0, 0.11, 21.5],
        [7.0, 0.125, 21.6],
        [8.0, 0.135, 22.0],
        [9.0, 0.14, 22.2],
        [10.0, 0.153, 22.3],
        [15.0, 0.19, 22.8],
        [20.0, 0.22, 23.2],
        [30.0, 0.26, 23.4],
        [40.0, 0.29, 23.5],
        [50.0, 0.35, 23.65],
        [100.0, 0.5, 23.7],
        [200.0, 0.7, 23.7],
    ];
    let mut dp = 0i32;
    let mut mean_pos: f32 = 0.0;
    for (i, &p) in pos.iter().enumerate() {
        if p == 0 {
            continue;
        }
        dp += p;
        mean_pos += (p * i as i32) as f32;
    }
    if dp < 2 {
        return f64::INFINITY;
    }
    mean_pos /= dp as f32;
    let mut mean_diff: f32 = 0.0;
    for (i, &p) in pos.iter().enumerate() {
        if p == 0 {
            continue;
        }
        mean_diff =
            (f64::from(mean_diff) + f64::from(p) * f64::from(i as f32 - mean_pos).abs()) as f32;
    }
    mean_diff /= dp as f32;
    let ipos = mean_diff as i32;
    if dp == 2 {
        return f64::from((2 * READLEN - 2 * (ipos + 1) - 1) * (ipos + 1) / (READLEN - 1))
            / (f64::from(READLEN) * 0.5);
    }
    let i = if dp >= 200 {
        PARAM.len()
    } else {
        PARAM
            .iter()
            .position(|p| p[0] >= dp as f32)
            .unwrap_or(PARAM.len())
    };
    let (pscale, pshift): (f32, f32) = if i == PARAM.len() {
        (PARAM[i - 1][1], PARAM[i - 1][2])
    } else if i > 0 && PARAM[i][0] != dp as f32 {
        (
            ((f64::from(PARAM[i - 1][1]) + f64::from(PARAM[i][1])) * 0.5) as f32,
            ((f64::from(PARAM[i - 1][2]) + f64::from(PARAM[i][2])) * 0.5) as f32,
        )
    } else {
        (PARAM[i][1], PARAM[i][2])
    };
    0.5 * kf_erfc(-f64::from((mean_diff - pshift) * pscale))
}

/// `calc_mwu_biasZ` with `do_Z` set: the Mann-Whitney U z-score of `a` vs `b`.
fn calc_mwu_bias_z(a: &[i32], b: &[i32]) -> f64 {
    let n = a.len();
    let b_empty = b.iter().all(|&x| x == 0);
    let (mut e, mut l, mut na, mut nb) = (0i32, 0i32, 0i32, 0i32);
    let mut t: i64 = 0;
    if b_empty {
        for i in (0..n).rev() {
            na = na.wrapping_add(a[i]);
            t += i64::from(a[i].wrapping_mul(a[i]).wrapping_sub(1).wrapping_mul(a[i]));
        }
    } else {
        for i in (0..n).rev() {
            e = e.wrapping_add(a[i].wrapping_mul(b[i]));
            l = l.wrapping_add(a[i].wrapping_mul(nb));
            na = na.wrapping_add(a[i]);
            nb = nb.wrapping_add(b[i]);
            let p = a[i] + b[i];
            t += i64::from(p.wrapping_mul(p).wrapping_sub(1).wrapping_mul(p));
        }
    }
    if na == 0 || nb == 0 {
        return f64::INFINITY;
    }
    let u = f64::from(l) + f64::from(e) * 0.5;
    let m = f64::from(na.wrapping_mul(nb)) / 2.0;
    let nanb = f64::from(na.wrapping_mul(nb));
    let nn = (na + nb).wrapping_mul(na + nb - 1);
    let var2 = nanb / 12.0 * (f64::from(na + nb + 1) - t as f64 / f64::from(nn));
    if var2 <= 0.0 {
        return 0.0;
    }
    (u - m) / var2.sqrt()
}

/// `logsumexp2`
fn logsumexp2(a: f64, b: f64) -> f64 {
    if a > b {
        (1.0 + (b - a).exp()).ln() + a
    } else {
        (1.0 + (a - b).exp()).ln() + b
    }
}

/// `calc_SegBias`
fn calc_seg_bias(bcr: &[CallRet], call: &mut Call) {
    call.seg_bias = f32::INFINITY;
    let nr = (call.anno[2] + call.anno[3]) as i32;
    if nr == 0 {
        return;
    }
    let n = call.n as i32;
    let avg_dp = ((call.anno[0] + call.anno[1] + f64::from(nr)) / f64::from(n)) as i32;
    let mut m = (f64::from(nr) / f64::from(avg_dp) + 0.5).floor();
    if m > f64::from(n) {
        m = f64::from(n);
    } else if m == 0.0 {
        m = 1.0;
    }
    let f = m / 2.0 / f64::from(n);
    let p = f64::from(nr) / f64::from(n);
    let q = f64::from(nr) / m;
    let log2 = 2f64.ln();
    let mut sum = 0.0;
    for r in bcr.iter().take(call.n) {
        let oi = (r.anno[2] + r.anno[3]) as i32;
        let tmp = if oi != 0 {
            let mut tmp = logsumexp2((2.0 * (1.0 - f)).ln(), f.ln() + f64::from(oi) * log2 - q);
            tmp += f.ln() + f64::from(oi) * (q / p).ln() - q + p;
            tmp
        } else {
            (2.0 * f * (1.0 - f) * (-q).exp() + f * f * (-2.0 * q).exp() + (1.0 - f) * (1.0 - f))
                .ln()
                + p
        };
        sum += tmp;
    }
    call.seg_bias = sum as f32;
}

/// `bcf_call_combine`: orders the alleles, sets PL and the annotations. Returns -1 for an
/// indel without a supported alternate allele.
pub fn combine(calls: &[CallRet], bca: &CallAux, ref_base: Option<u8>, call: &mut Call) -> i32 {
    let n = call.n;
    let ref4: i32;
    if let Some(rb) = ref_base {
        ref4 = i32::from(SEQ_NT16_INT[rb as usize]).min(4);
        call.ori_ref = i32::from(SEQ_NT16_INT[rb as usize]);
    } else {
        call.ori_ref = -1;
        ref4 = 0;
    }
    let mut qsum = [0f32; MAX_ALLELES];
    for c in calls.iter().take(n) {
        let mut sum: f32 = 0.0;
        for &q in &c.qs[..4] {
            sum += q as f32;
        }
        if sum != 0.0 {
            for (qs, &q) in qsum.iter_mut().zip(&c.qs[..4]) {
                *qs += q as f32 / sum;
            }
        }
    }
    // Insertion sort of indices by qsum, ascending (stable).
    let mut ptr = [0usize, 1, 2, 3, 4];
    for i in 1..4 {
        let mut j = i;
        while j > 0 && qsum[ptr[j]] < qsum[ptr[j - 1]] {
            ptr.swap(j, j - 1);
            j -= 1;
        }
    }
    call.a = [-1; 5];
    call.qsum = [0.0; MAX_ALLELES];
    call.unseen = -1;
    call.a[0] = ref4;
    let mut j = 1usize;
    let mut i: i32 = 3;
    while i >= 0 {
        let ipos = ptr[i as usize];
        if ipos as i32 == ref4 {
            call.qsum[0] = qsum[ipos];
        } else {
            if qsum[ipos] == 0.0 {
                break;
            }
            call.qsum[j] = qsum[ipos];
            call.a[j] = ipos as i32;
            j += 1;
        }
        i -= 1;
    }
    if ref_base.is_some() {
        if ((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0 {
            call.unseen = j as i32;
            call.a[j] = ptr[i as usize] as i32;
            j += 1;
        }
        call.n_alleles = j;
    } else {
        call.n_alleles = j;
        if j == 1 {
            return -1;
        }
    }
    let has_alt = !(call.n_alleles == 2 && call.unseen != -1);
    call.has_alt = has_alt;

    let na = call.n_alleles;
    let x = na * (na + 1) / 2;
    let mut g = [0usize; 15];
    let mut z = 0;
    for i in 0..na {
        for j in 0..=i {
            g[z] = (call.a[j] * 5 + call.a[i]) as usize;
            z += 1;
        }
    }
    for (i, r) in calls.iter().take(n).enumerate() {
        let pl = &mut call.pl[x * i..x * (i + 1)];
        let mut min = f32::MAX;
        for &gj in &g[..x] {
            if min > r.p[gj] {
                min = r.p[gj];
            }
        }
        for (k, &gj) in g[..x].iter().enumerate() {
            let y = c_int(f64::from(r.p[gj] - min) + 0.499);
            pl[k] = y.min(255);
        }
    }
    for (i, r) in calls.iter().take(n).enumerate() {
        call.dp4[4 * i] = r.anno[0] as i32;
        call.dp4[4 * i + 1] = r.anno[1] as i32;
        call.dp4[4 * i + 2] = r.anno[2] as i32;
        call.dp4[4 * i + 3] = r.anno[3] as i32;
    }
    if bca.fmt_flag & (INFO_SCR | FMT_SCR) != 0 {
        for (i, r) in calls.iter().take(n).enumerate() {
            call.scr[0] += r.scr;
            call.scr[1 + i] = r.scr;
        }
    }
    if bca.fmt_flag & NEEDS_AD != 0 {
        for (i, r) in calls.iter().take(n).enumerate() {
            for jj in 0..na {
                let v = r.adr[call.a[jj] as usize];
                call.adr[jj] += v;
                call.adr[MAX_ALLELES + i * na + jj] = v;
                let v = r.adf[call.a[jj] as usize];
                call.adf[jj] += v;
                call.adf[MAX_ALLELES + i * na + jj] = v;
            }
        }
    }
    if bca.fmt_flag & FMT_QS != 0 {
        for (i, r) in calls.iter().take(n).enumerate() {
            for jj in 0..na {
                call.qs[i * na + jj] = r.qs[call.a[jj] as usize];
            }
        }
    }
    if bca.fmt_flag & FMT_QM != 0 {
        for (i, r) in calls.iter().take(n).enumerate() {
            for jj in 0..na {
                let jold = call.a[jj] as usize;
                let dp = call.adr[MAX_ALLELES + i * na + jj] + call.adf[MAX_ALLELES + i * na + jj];
                let val = if dp != 0 {
                    -4.3429 * (r.qm[jold] / f64::from(dp)).ln()
                } else {
                    0.0
                };
                // A depth from reads that add no quality (`--ambig-reads incAD0`) makes the
                // mean infinite, which converts to INT_MIN (missing).
                call.qm[i * na + jj] = c_int(val + 0.499);
            }
        }
    }

    call.anno = [0.0; 16];
    call.ori_depth = 0;
    call.depth = 0;
    call.mq0 = 0;
    for r in calls.iter().take(n) {
        call.depth = call
            .depth
            .wrapping_add((r.anno[0] + r.anno[1] + r.anno[2] + r.anno[3]) as u32);
        call.ori_depth += r.ori_depth;
        call.mq0 += r.mq0;
        for j in 0..16 {
            call.anno[j] += r.anno[j];
        }
    }

    if !has_alt {
        return 0;
    }
    if bca.fmt_flag & INFO_FS != 0 {
        call.strand_bias = kt_fisher_exact(
            call.anno[0] as i32,
            call.anno[1] as i32,
            call.anno[2] as i32,
            call.anno[3] as i32,
        )
        .0 as f32;
    }
    if bca.fmt_flag & INFO_SGB != 0 {
        calc_seg_bias(calls, call);
    }
    if call.ori_ref < 0 {
        if bca.fmt_flag & INFO_RPBZ != 0 {
            call.mwu_pos = calc_mwu_bias_z(&bca.iref_pos, &bca.ialt_pos) as f32;
        }
        if bca.fmt_flag & INFO_MQBZ != 0 {
            call.mwu_mq = calc_mwu_bias_z(&bca.iref_mq, &bca.ialt_mq) as f32;
        }
        if bca.fmt_flag & INFO_SCBZ != 0 {
            call.mwu_sc = calc_mwu_bias_z(&bca.iref_scl, &bca.ialt_scl) as f32;
        }
    } else {
        if bca.fmt_flag & INFO_RPBZ != 0 {
            call.mwu_pos = calc_mwu_bias_z(&bca.ref_pos, &bca.alt_pos) as f32;
        }
        if bca.fmt_flag & INFO_MQBZ != 0 {
            call.mwu_mq = calc_mwu_bias_z(&bca.ref_mq, &bca.alt_mq) as f32;
        }
        if bca.fmt_flag & INFO_BQBZ != 0 {
            call.mwu_bq = calc_mwu_bias_z(&bca.ref_bq, &bca.alt_bq) as f32;
        }
        if bca.fmt_flag & INFO_MQSBZ != 0 {
            call.mwu_mqs = calc_mwu_bias_z(&bca.fwd_mqs, &bca.rev_mqs) as f32;
        }
        if bca.fmt_flag & INFO_SCBZ != 0 {
            call.mwu_sc = calc_mwu_bias_z(&bca.ref_scl, &bca.alt_scl) as f32;
        }
    }
    if bca.fmt_flag & INFO_NMBZ != 0 {
        call.mwu_nm[0] = calc_mwu_bias_z(&bca.ref_nm, &bca.alt_nm) as f32;
    }
    if bca.fmt_flag & FMT_NMBZ != 0 {
        for (i, r) in calls.iter().take(n).enumerate() {
            let val = calc_mwu_bias_z(&r.ref_nm, &r.alt_nm) as f32;
            call.mwu_nm[i + 1] = if val == f32::INFINITY { 0.0 } else { val };
        }
    }
    if bca.fmt_flag & INFO_VDB != 0 {
        call.vdb = calc_vdb(&bca.alt_pos) as f32;
    }
    0
}

/// `bcf_call2bcf`: fills `rec` from the call. `reference` is the sequence (for indels).
///
/// # Errors
///
/// Returns an error if a field cannot be set (an annotation missing from the header).
pub fn call2bcf(
    bc: &mut Call,
    rec: &mut Record,
    hdr: &Header,
    fmt_flag: u32,
    bca: &CallAux,
    reference: &[u8],
) -> anyhow::Result<()> {
    rec.clear();
    rec.set_rid(bc.tid);
    rec.set_pos(bc.pos);
    rec.set_qual(Some(0.0));
    let ref_at = |i: i64| -> u8 { reference.get(i as usize).copied().unwrap_or(0) };

    let mut nals = 1usize;
    let mut has_alt = false;
    let mut tmp = std::mem::take(&mut bc.tmp);
    tmp.clear();
    if bc.ori_ref < 0 {
        tmp.push(ref_at(bc.pos));
        for j in 0..bca.indelreg {
            tmp.push(ref_at(bc.pos + 1 + i64::from(j)));
        }
        for i in 1..4 {
            if bc.a[i] < 0 {
                break;
            }
            tmp.push(b',');
            tmp.push(ref_at(bc.pos));
            let t = bca.indel_types[bc.a[i] as usize];
            if t < 0 {
                for j in -t..bca.indelreg {
                    tmp.push(ref_at(bc.pos + 1 + i64::from(j)));
                }
            } else {
                let off = bc.a[i] as usize * bca.maxins as usize;
                for j in 0..t as usize {
                    tmp.push(b"ACGTN"[bca.inscns[off + j] as usize]);
                }
                for j in 0..bca.indelreg {
                    tmp.push(ref_at(bc.pos + 1 + i64::from(j)));
                }
            }
            nals += 1;
            has_alt = true;
        }
    } else {
        tmp.push(b"ACGTN"[bc.ori_ref as usize]);
        for i in 1..5 {
            if bc.a[i] < 0 {
                break;
            }
            tmp.push(b',');
            if bc.unseen == i as i32 {
                tmp.extend_from_slice(b"<*>");
            } else {
                tmp.push(b"ACGT"[bc.a[i] as usize]);
                has_alt = true;
            }
            nals += 1;
        }
    }
    rec.update_alleles_str(hdr, &tmp);
    bc.tmp = tmp;
    let n_allele = rec.allele_count();

    if bc.ori_ref < 0 {
        rec.update_info_flag(hdr, b"INDEL", true)?;
        let mut idv = bca.max_support as i32;
        if fmt_flag & INFO_IMF != 0 {
            let max_frac: f32;
            if fmt_flag & NEEDS_AD != 0 && bca.edlib {
                let mut max_ad = 0;
                for k in 1..n_allele {
                    max_ad = max_ad.max(bc.adf[k] + bc.adr[k]);
                }
                max_frac = (f64::from(max_ad) / f64::from(bc.ori_depth)) as f32;
                idv = max_ad;
            } else {
                max_frac = bca.max_frac;
            }
            if fmt_flag & INFO_IDV != 0 {
                rec.update_info_i32(hdr, b"IDV", &[idv])?;
            }
            rec.update_info_f32(hdr, b"IMF", &[max_frac])?;
        } else if fmt_flag & INFO_IDV != 0 {
            rec.update_info_i32(hdr, b"IDV", &[idv])?;
        }
    }
    rec.update_info_i32(hdr, b"DP", &[bc.ori_depth as i32])?;
    if fmt_flag & INFO_ADF != 0 {
        rec.update_info_i32(hdr, b"ADF", &bc.adf[..n_allele])?;
    }
    if fmt_flag & INFO_ADR != 0 {
        rec.update_info_i32(hdr, b"ADR", &bc.adr[..n_allele])?;
    }
    if fmt_flag & (INFO_AD | INFO_DPR) != 0 {
        for i in 0..n_allele {
            bc.adf[i] += bc.adr[i];
        }
        if fmt_flag & INFO_AD != 0 {
            rec.update_info_i32(hdr, b"AD", &bc.adf[..n_allele])?;
        }
        if fmt_flag & INFO_DPR != 0 {
            rec.update_info_i32(hdr, b"DPR", &bc.adf[..n_allele])?;
        }
    }
    if fmt_flag & INFO_SCR != 0 {
        rec.update_info_i32(hdr, b"SCR", &bc.scr[..1])?;
    }
    let tmpf: [f32; 16] = std::array::from_fn(|i| bc.anno[i] as f32);
    rec.update_info_f32(hdr, b"I16", &tmpf)?;
    rec.update_info_f32(hdr, b"QS", &bc.qsum[..nals])?;

    if has_alt {
        if fmt_flag & INFO_VDB != 0 && bc.vdb != f32::INFINITY {
            rec.update_info_f32(hdr, b"VDB", &[bc.vdb])?;
        }
        if fmt_flag & INFO_SGB != 0 && bc.seg_bias != f32::INFINITY {
            rec.update_info_f32(hdr, b"SGB", &[bc.seg_bias])?;
        }
        if fmt_flag & INFO_NM != 0 && (bca.nnm[0] != 0 || bca.nnm[1] != 0) {
            for i in 0..2 {
                bc.nm[i] = if bca.nnm[i] != 0 {
                    bca.nm[i] / bca.nnm[i] as f32
                } else {
                    0.0
                };
            }
            rec.update_info_f32(hdr, b"NM", &bc.nm)?;
        }
        let fields: [(u32, &[u8], f32); 6] = [
            (INFO_RPBZ, b"RPBZ", bc.mwu_pos),
            (INFO_MQBZ, b"MQBZ", bc.mwu_mq),
            (INFO_MQSBZ, b"MQSBZ", bc.mwu_mqs),
            (INFO_BQBZ, b"BQBZ", bc.mwu_bq),
            (INFO_NMBZ, b"NMBZ", bc.mwu_nm[0]),
            (INFO_SCBZ, b"SCBZ", bc.mwu_sc),
        ];
        for (flag, key, v) in fields {
            if fmt_flag & flag != 0 && v != f32::INFINITY {
                rec.update_info_f32(hdr, key, &[v])?;
            }
        }
        if fmt_flag & INFO_FS != 0 && bc.strand_bias != f32::INFINITY {
            rec.update_info_f32(hdr, b"FS", &[bc.strand_bias])?;
        }
    }
    if fmt_flag & INFO_MQ0F != 0 {
        let v = if bc.ori_depth != 0 {
            bc.mq0 as f32 / bc.ori_depth as f32
        } else {
            0.0
        };
        rec.update_info_f32(hdr, b"MQ0F", &[v])?;
    }

    let n = bc.n;
    rec.update_format_i32(hdr, b"PL", &bc.pl[..nals * (nals + 1) / 2 * n])?;
    let mut ibuf = std::mem::take(&mut bc.itmp);
    if fmt_flag & FMT_DP != 0 {
        ibuf.clear();
        ibuf.extend(
            (0..n)
                .map(|i| bc.dp4[4 * i] + bc.dp4[4 * i + 1] + bc.dp4[4 * i + 2] + bc.dp4[4 * i + 3]),
        );
        rec.update_format_i32(hdr, b"DP", &ibuf)?;
    }
    if fmt_flag & FMT_DV != 0 {
        ibuf.clear();
        ibuf.extend((0..n).map(|i| bc.dp4[4 * i + 2] + bc.dp4[4 * i + 3]));
        rec.update_format_i32(hdr, b"DV", &ibuf)?;
    }
    if fmt_flag & FMT_SP != 0 {
        ibuf.clear();
        for i in 0..n {
            let (fwd_ref, rev_ref, fwd_alt, rev_alt) = (
                bc.dp4[4 * i],
                bc.dp4[4 * i + 1],
                bc.dp4[4 * i + 2],
                bc.dp4[4 * i + 3],
            );
            if fwd_ref + rev_ref < 2
                || fwd_alt + rev_alt < 2
                || fwd_ref + fwd_alt < 2
                || rev_ref + rev_alt < 2
            {
                ibuf.push(0);
            } else {
                let (_, two) = kt_fisher_exact(fwd_ref, rev_ref, fwd_alt, rev_alt);
                let x = c_int(-4.343 * two.ln() + 0.499);
                ibuf.push(x.min(255));
            }
        }
        rec.update_format_i32(hdr, b"SP", &ibuf)?;
    }
    if fmt_flag & FMT_DP4 != 0 {
        rec.update_format_i32(hdr, b"DP4", &bc.dp4[..n * 4])?;
    }
    let nsa = n * n_allele;
    if fmt_flag & FMT_ADF != 0 {
        rec.update_format_i32(hdr, b"ADF", &bc.adf[MAX_ALLELES..MAX_ALLELES + nsa])?;
    }
    if fmt_flag & FMT_ADR != 0 {
        rec.update_format_i32(hdr, b"ADR", &bc.adr[MAX_ALLELES..MAX_ALLELES + nsa])?;
    }
    if fmt_flag & (FMT_AD | FMT_DPR) != 0 {
        for i in 0..nsa {
            bc.adf[MAX_ALLELES + i] += bc.adr[MAX_ALLELES + i];
        }
        if fmt_flag & FMT_AD != 0 {
            rec.update_format_i32(hdr, b"AD", &bc.adf[MAX_ALLELES..MAX_ALLELES + nsa])?;
        }
        if fmt_flag & FMT_DPR != 0 {
            rec.update_format_i32(hdr, b"DPR", &bc.adf[MAX_ALLELES..MAX_ALLELES + nsa])?;
        }
    }
    if fmt_flag & FMT_SCR != 0 {
        rec.update_format_i32(hdr, b"SCR", &bc.scr[1..=n])?;
    }
    if fmt_flag & FMT_QS != 0 {
        rec.update_format_i32(hdr, b"QS", &bc.qs[..nsa])?;
    }
    if fmt_flag & FMT_QM != 0 {
        rec.update_format_i32(hdr, b"QM", &bc.qm[..nsa])?;
    }
    if has_alt && fmt_flag & FMT_NMBZ != 0 {
        rec.update_format_f32(hdr, b"NMBZ", &bc.mwu_nm[1..=n])?;
    }
    bc.itmp = ibuf;
    Ok(())
}
