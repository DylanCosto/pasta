//! The pairing of records at one position across readers (HTSlib's `bcf_sr_sort.c`).
//!
//! When several readers have records at the current position, the records are grouped into
//! variant sets of compatible records, one record per reader at most, and the sets are
//! returned one by one. Which records are compatible is set by [`Pairing`]; the order follows
//! HTSlib's greedy merging: the set covering the most readers first takes the best-scoring
//! compatible set, until no set can be merged, and is then output.

use std::{collections::VecDeque, ops::Range};

use bitflags::bitflags;

use super::ReaderState;
use crate::variant::{Record, Unpack, VariantType};

bitflags! {
    /// Which records at one position are returned together (HTSlib's `BCF_SR_PAIR_*`).
    ///
    /// The empty set means [`Pairing::EXACT`], which is what the synced reader uses when
    /// nothing is set (HTSlib's `COLLAPSE_NONE`).
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
    pub struct Pairing: u32 {
        /// SNPs pair with SNPs even when their alleles differ (`BCF_SR_PAIR_SNPS`,
        /// `COLLAPSE_SNPS`).
        const SNPS = 1;
        /// Indels pair with indels even when their alleles differ (`BCF_SR_PAIR_INDELS`,
        /// `COLLAPSE_INDELS`).
        const INDELS = 1 << 1;
        /// Any records pair (`BCF_SR_PAIR_ANY`, `COLLAPSE_ANY`).
        const ANY = 1 << 2;
        /// Records sharing at least one ALT allele pair (`BCF_SR_PAIR_SOME`,
        /// `COLLAPSE_SOME`).
        const SOME = 1 << 3;
        /// REF-only records pair with SNPs (`BCF_SR_PAIR_SNP_REF`).
        const SNP_REF = 1 << 4;
        /// REF-only records pair with indels (`BCF_SR_PAIR_INDEL_REF`).
        const INDEL_REF = 1 << 5;
        /// Only records with the same set of alleles pair (`BCF_SR_PAIR_EXACT`).
        const EXACT = 1 << 6;
        /// Records must also have the same ID (`BCF_SR_PAIR_ID`).
        const ID = 1 << 7;
        /// SNPs with SNPs and indels with indels (`BCF_SR_PAIR_BOTH`, `COLLAPSE_BOTH`).
        const BOTH = Self::SNPS.bits() | Self::INDELS.bits();
        /// As [`Self::BOTH`], and REF-only records with either (`BCF_SR_PAIR_BOTH_REF`).
        const BOTH_REF = Self::BOTH.bits() | Self::SNP_REF.bits() | Self::INDEL_REF.bits();
    }
}

/// Variant type bits of the pairing scores (`SR_*`).
const SR_REF: u32 = 1;
const SR_SNP: u32 = 2;
const SR_INDEL: u32 = 4;
const SR_OTHER: u32 = 8;

/// Identifies a buffered record by the address of its box, as HTSlib identifies records by
/// pointer. The address is only compared, never dereferenced.
pub(super) type RecId = usize;

/// Returns a buffered record's identity.
#[inline]
pub(super) fn rec_id(rec: &Record) -> RecId {
    std::ptr::from_ref(rec) as usize
}

/// A variant: identical records (by their allele string) across readers (`var_t`).
#[derive(Debug, Default)]
struct Var {
    /// The lookup key: the allele string, with a counter appended for duplicates within one
    /// reader.
    key: Vec<u8>,
    /// The allele string, such as `A>C,A>G`.
    text: Vec<u8>,
    /// The variant types: HTSlib's `VCF_*` bits, then the `SR_*` bits.
    ty: u32,
    nalt: i64,
    /// The readers with this variant, and their records.
    vcf: Vec<usize>,
    rec: Vec<RecId>,
    /// The groups containing the variant.
    mask: Vec<u64>,
}

/// A group of readers with the same combination of records at the position (`grp_t`).
#[derive(Debug, Default)]
struct Grp {
    key: Vec<u8>,
    var: Vec<usize>,
    nvcf: usize,
}

/// A set of compatible variants that can be output together (`varset_t`).
#[derive(Debug, Default, Clone)]
struct VarSet {
    var: Vec<usize>,
    /// The number of readers covered (HTSlib's `srt->cnt`).
    cnt: i64,
    mask: Vec<u64>,
}

/// HTSlib's `sr_sort_t`.
#[derive(Debug)]
pub(super) struct Sorter {
    pub(super) pair: Pairing,
    score: [u8; 256],
    initialized: bool,
    vars: Vec<Var>,
    nvar: usize,
    grps: Vec<Grp>,
    ngrp: usize,
    vsets: Vec<VarSet>,
    nvset: usize,
    /// Records in output order for each reader; `None` where a reader has no record.
    vcf_buf: Vec<VecDeque<Option<RecId>>>,
    /// The allele strings of one reader's records, `;`-separated (HTSlib's `str`).
    text: Vec<u8>,
    offs: Vec<(usize, usize)>,
    candidate: Vec<u8>,
    grp_key: Vec<u8>,
    grp_vars: Vec<usize>,
    chr: Option<Vec<u8>>,
    pos: i64,
    nsr: usize,
    /// The readers with records at the current position.
    pub(super) active: Vec<usize>,
}

impl Default for Sorter {
    fn default() -> Self {
        Self {
            pair: Pairing::empty(),
            score: [0; 256],
            initialized: false,
            vars: Vec::new(),
            nvar: 0,
            grps: Vec::new(),
            ngrp: 0,
            vsets: Vec::new(),
            nvset: 0,
            vcf_buf: Vec::new(),
            text: Vec::new(),
            offs: Vec::new(),
            candidate: Vec::new(),
            grp_key: Vec::new(),
            grp_vars: Vec::new(),
            chr: None,
            pos: 0,
            nsr: 0,
            active: Vec::new(),
        }
    }
}

#[inline]
fn score_index(a: u32, b: u32) -> usize {
    ((a << 4) | b) as usize
}

fn masks_intersect(a: &[u64], b: &[u64]) -> bool {
    a.iter().zip(b).any(|(x, y)| x & y != 0)
}

/// Compares the comma-separated parts of two allele strings, ignoring case (`strncasecmp`).
fn parts(s: &[u8]) -> impl Iterator<Item = &[u8]> {
    // An empty string has no parts, as in HTSlib's loops.
    s.split(|&c| c == b',')
        .take(if s.is_empty() { 0 } else { usize::MAX })
}

/// Whether two variants list the same alleles in any order (`multi_is_exact`).
fn multi_is_exact(a: &Var, b: &Var) -> bool {
    if a.nalt != b.nalt || a.text.len() != b.text.len() {
        return false;
    }
    parts(&a.text).all(|x| parts(&b.text).any(|y| x.eq_ignore_ascii_case(y)))
}

/// Whether two variants share an allele (`multi_is_subset`).
fn multi_is_subset(a: &Var, b: &Var) -> bool {
    parts(&a.text).any(|x| parts(&b.text).any(|y| x.eq_ignore_ascii_case(y)))
}

impl Sorter {
    /// Forgets the current position (`bcf_sr_sort_reset`).
    pub(super) fn reset(&mut self) {
        self.chr = None;
    }

    /// Starts a new position with reader `i` (`bcf_sr_sort_set_active`).
    pub(super) fn set_active(&mut self, i: usize) {
        self.active.clear();
        self.active.push(i);
    }

    /// Adds reader `i` at the current position (`bcf_sr_sort_add_active`).
    pub(super) fn add_active(&mut self, i: usize) {
        self.active.push(i);
    }

    /// Drops reader `i`'s queue (`bcf_sr_sort_remove_reader`).
    pub(super) fn remove_reader(&mut self, i: usize) {
        if i < self.vcf_buf.len() {
            self.vcf_buf.remove(i);
        }
    }

    /// HTSlib's `bcf_sr_init_scores`.
    fn init_scores(&mut self) {
        if self.pair.contains(Pairing::ANY) {
            self.pair |= Pairing::SNPS | Pairing::INDELS | Pairing::SNP_REF | Pairing::INDEL_REF;
        }
        if self.pair.contains(Pairing::SNPS) {
            self.score[score_index(SR_SNP, SR_SNP)] = 3;
        }
        if self.pair.contains(Pairing::INDELS) {
            self.score[score_index(SR_INDEL, SR_INDEL)] = 3;
        }
        if self.pair.contains(Pairing::SNP_REF) {
            self.score[score_index(SR_SNP, SR_REF)] = 2;
            self.score[score_index(SR_REF, SR_SNP)] = 2;
        }
        if self.pair.contains(Pairing::INDEL_REF) {
            self.score[score_index(SR_INDEL, SR_REF)] = 2;
            self.score[score_index(SR_REF, SR_INDEL)] = 2;
        }
        if self.pair.contains(Pairing::ANY) {
            for s in &mut self.score {
                if *s == 0 {
                    *s = 1;
                }
            }
        }
        // Combinations of types score as their best-scoring pair of single types.
        for i in 0..256u32 {
            if self.score[i as usize] != 0 {
                continue;
            }
            let mut max = 0;
            for jbit in 0..4 {
                let j = 1u32 << jbit;
                if i & (j << 4) == 0 {
                    continue;
                }
                for kbit in 0..4 {
                    let k = 1u32 << kbit;
                    if i & k == 0 {
                        continue;
                    }
                    max = max.max(self.score[score_index(j, k)]);
                }
            }
            self.score[i as usize] = max;
        }
    }

    /// HTSlib's `pairing_score`.
    fn pairing_score(&self, ivset: usize, jvset: usize) -> u32 {
        let iv = &self.vsets[ivset];
        let jv = &self.vsets[jvset];
        let mut min = u32::MAX;
        for &i in &iv.var {
            let ivar = &self.vars[i];
            for &j in &jv.var {
                let jvar = &self.vars[j];
                if self.pair.contains(Pairing::EXACT) {
                    if ivar.ty != jvar.ty {
                        continue;
                    }
                    if ivar.text == jvar.text || multi_is_exact(ivar, jvar) {
                        return u32::MAX;
                    }
                    continue;
                }
                if ivar.ty == jvar.ty && ivar.text == jvar.text {
                    return u32::MAX;
                }
                if ivar.ty & jvar.ty != 0 && multi_is_subset(ivar, jvar) {
                    return u32::MAX;
                }
                let score = u32::from(self.score[score_index(ivar.ty, jvar.ty)]);
                if score == 0 {
                    return 0;
                }
                min = min.min(score);
            }
        }
        if self.pair.contains(Pairing::EXACT) {
            return 0;
        }
        debug_assert!(min != u32::MAX, "variant sets are never empty");
        let cnt: usize = iv
            .var
            .iter()
            .chain(&jv.var)
            .map(|&v| self.vars[v].vcf.len())
            .sum();
        (1u32 << (28 + min)).wrapping_add(cnt as u32)
    }

    /// Removes a variant set, keeping its storage at the end (`remove_vset`).
    fn remove_vset(&mut self, j: usize) {
        self.vsets[j..self.nvset].rotate_left(1);
        self.nvset -= 1;
    }

    /// Merges two variant sets into the earlier one (`merge_vsets`); returns its index.
    fn merge_vsets(&mut self, a: usize, b: usize) -> usize {
        let (i, j) = if a > b { (b, a) } else { (a, b) };
        let (head, tail) = self.vsets.split_at_mut(j);
        let (iv, jv) = (&mut head[i], &tail[0]);
        for (x, y) in iv.mask.iter_mut().zip(&jv.mask) {
            *x |= y;
        }
        iv.var.extend_from_slice(&jv.var);
        iv.cnt += jv.cnt;
        self.remove_vset(j);
        i
    }

    /// Queues a variant set's records for output (`push_vset`).
    fn push_vset(&mut self, ivset: usize, nreaders: usize) {
        for buf in &mut self.vcf_buf[..nreaders] {
            buf.push_back(None);
        }
        for &v in &self.vsets[ivset].var {
            let var = &self.vars[v];
            for (&jvcf, &rec) in var.vcf.iter().zip(&var.rec) {
                if let Some(slot) = self.vcf_buf[jvcf].back_mut() {
                    *slot = Some(rec);
                }
            }
        }
        self.remove_vset(ivset);
    }

    /// Finds the variant with this key among those of the current position.
    fn find_var(&self, key: &[u8]) -> Option<usize> {
        self.vars[..self.nvar].iter().position(|v| v.key == key)
    }

    /// Builds the allele string of one record into `self.text` (as `bcf_sr_sort_set` does).
    fn append_record_text(&mut self, line: &mut Record, header: &crate::variant::Header) {
        line.unpack(Unpack::STR);
        if self.pair.contains(Pairing::ID) {
            self.text.extend_from_slice(line.id());
            self.text.push(b':');
        }
        let n_allele = line.allele_count();
        let mut end_pos: Option<i32> = None;
        for ivar in 1..n_allele {
            if ivar > 1 {
                self.text.push(b',');
            }
            self.text.extend_from_slice(line.allele(0));
            self.text.push(b'>');
            self.text.extend_from_slice(line.allele(ivar));
            // Symbolic alleles also use END, for several events (<DEL>s, say) at one POS.
            if line.allele(ivar).first() == Some(&b'<') {
                let end = *end_pos.get_or_insert_with(|| {
                    line.info_field(header, b"END")
                        .map_or(0, |f| f.int_value() as i32)
                });
                if end != 0 {
                    self.text.push(b'/');
                    self.text.extend_from_slice(end.to_string().as_bytes());
                }
            }
        }
        if n_allele == 1 {
            self.text.extend_from_slice(line.allele(0));
            self.text.extend_from_slice(b">.");
        }
    }

    /// Groups the records at a position into variant sets in output order
    /// (`bcf_sr_sort_set`).
    fn set(&mut self, readers: &mut [ReaderState], chr: &[u8], min_pos: i64) {
        if !self.initialized {
            if self.pair.is_empty() {
                self.pair = Pairing::EXACT;
            }
            self.init_scores();
            self.initialized = true;
        }
        self.nvar = 0;
        self.ngrp = 0;
        self.nvset = 0;
        for buf in &mut self.vcf_buf {
            buf.clear();
        }

        // Group the readers by their combination of records at the position.
        let active = std::mem::take(&mut self.active);
        for &ireader in &active {
            let reader = &mut readers[ireader];
            let rid = reader.header.contig_id(chr).map_or(-1, |r| r as i32);
            self.add_reader_vars(reader, ireader, rid, min_pos);
            self.make_group_key();
            let igrp = self.find_or_add_group();
            self.grps[igrp].nvcf += 1;
        }
        self.active = active;

        self.init_vsets();
        self.pair_vsets(readers.len());

        self.chr = Some(chr.to_vec());
        self.pos = min_pos;
    }

    /// Adds reader `ireader`'s records at the position to the variants, leaving their
    /// allele strings in `text` and `offs` and their variants in `grp_vars`.
    #[inline]
    fn add_reader_vars(
        &mut self,
        reader: &mut ReaderState,
        ireader: usize,
        rid: i32,
        min_pos: i64,
    ) {
        self.grp_vars.clear();
        self.offs.clear();
        self.text.clear();
        for irec in 1..=reader.nbuffer {
            let line = &mut *reader.buffer[irec];
            if line.rid() != rid || line.pos() != min_pos {
                break;
            }
            if !self.text.is_empty() {
                self.text.push(b';');
            }
            let beg = self.text.len();
            self.append_record_text(line, &reader.header);
            let var_end = self.text.len();
            self.offs.push((beg, var_end));
            let ivar = self.find_or_add_var(ireader, beg..var_end);
            let ty = line.variant_types().bits();
            let id = rec_id(line);
            let var = &mut self.vars[ivar];
            var.nalt = line.allele_count() as i64 - 1;
            var.ty = ty;
            var.vcf.push(ireader);
            var.rec.push(id);
            self.grp_vars.push(ivar);
        }
    }

    /// Returns the variant of reader `ireader`'s record whose allele string is `text[span]`,
    /// adding it when it is new.
    #[inline]
    fn find_or_add_var(&mut self, ireader: usize, span: Range<usize>) -> usize {
        // Duplicate records within one reader (dbSNP has them) are separate variants that can
        // still match between readers: a counter tells them apart.
        self.candidate.clear();
        self.candidate.extend_from_slice(&self.text[span.clone()]);
        let mut var_idx = 0;
        let found = loop {
            match self.find_var(&self.candidate) {
                None => break None,
                Some(i) if self.vars[i].vcf.last() != Some(&ireader) => break Some(i),
                Some(_) => {
                    self.candidate.truncate(span.len());
                    self.candidate
                        .extend_from_slice(var_idx.to_string().as_bytes());
                    var_idx += 1;
                }
            }
        };
        if let Some(i) = found {
            return i;
        }
        if self.nvar == self.vars.len() {
            self.vars.push(Var::default());
        }
        let var = &mut self.vars[self.nvar];
        var.key.clear();
        var.key.extend_from_slice(&self.candidate);
        var.text.clear();
        var.text.extend_from_slice(&self.text[span]);
        var.vcf.clear();
        var.rec.clear();
        self.nvar += 1;
        self.nvar - 1
    }

    /// Returns the group whose key is `grp_key`, adding it (with the variants in `grp_vars`)
    /// when it is new.
    #[inline]
    fn find_or_add_group(&mut self) -> usize {
        if let Some(i) = self.grps[..self.ngrp]
            .iter()
            .position(|g| g.key == self.grp_key)
        {
            return i;
        }
        if self.ngrp == self.grps.len() {
            self.grps.push(Grp::default());
        }
        let g = &mut self.grps[self.ngrp];
        g.key.clear();
        g.key.extend_from_slice(&self.grp_key);
        g.var.clear();
        g.var.extend_from_slice(&self.grp_vars);
        g.nvcf = 0;
        self.ngrp += 1;
        self.ngrp - 1
    }

    /// Marks which groups contain each variant, and starts with one variant set per variant.
    #[inline]
    fn init_vsets(&mut self) {
        let words = self.ngrp.div_ceil(64);
        for var in &mut self.vars[..self.nvar] {
            var.mask.clear();
            var.mask.resize(words, 0);
        }
        for (igrp, g) in self.grps[..self.ngrp].iter().enumerate() {
            for &v in &g.var {
                self.vars[v].mask[igrp / 64] |= 1 << (igrp % 64);
            }
        }
        for ivar in 0..self.nvar {
            if self.nvset == self.vsets.len() {
                self.vsets.push(VarSet::default());
            }
            let var = &mut self.vars[ivar];
            let vset = &mut self.vsets[self.nvset];
            vset.var.clear();
            vset.var.push(ivar);
            vset.cnt = var.vcf.len() as i64;
            vset.mask.clear();
            vset.mask.extend_from_slice(&var.mask);
            self.nvset += 1;
            var.ty = sr_type(var.ty);
        }
    }

    /// Pairs the variant sets, largest first, and moves them to the output in that order.
    #[inline]
    fn pair_vsets(&mut self, nreaders: usize) {
        while self.nvset > 0 {
            let mut imax = 0;
            for ivset in 1..self.nvset {
                if self.vsets[imax].cnt < self.vsets[ivset].cnt {
                    imax = ivset;
                }
            }
            let mut ipair = None;
            let mut max_score = 0;
            for ivset in 0..self.nvset {
                if masks_intersect(&self.vsets[imax].mask, &self.vsets[ivset].mask) {
                    continue;
                }
                let score = self.pairing_score(imax, ivset);
                if max_score < score {
                    max_score = score;
                    ipair = Some(ivset);
                }
            }
            if let Some(ipair) = ipair.filter(|&p| p != imax) {
                self.merge_vsets(imax, ipair);
                continue;
            }
            self.push_vset(imax, nreaders);
        }
    }

    /// The group key: the records' allele strings sorted and `;`-joined (`grp_create_key`).
    fn make_group_key(&mut self) {
        self.grp_key.clear();
        if self.text.is_empty() || self.offs.is_empty() {
            return;
        }
        let text = &self.text;
        self.offs
            .sort_unstable_by(|a, b| text[a.0..a.1].cmp(&text[b.0..b.1]));
        for (i, &(b, e)) in self.offs.iter().enumerate() {
            if i > 0 {
                self.grp_key.push(b';');
            }
            self.grp_key.extend_from_slice(&text[b..e]);
        }
    }

    /// Moves the next set of paired records into the readers' current slots
    /// (`bcf_sr_sort_next`); returns the number of readers with a record.
    pub(super) fn next(
        &mut self,
        readers: &mut [ReaderState],
        has_line: &mut [bool],
        chr: &[u8],
        min_pos: i64,
    ) -> usize {
        debug_assert!(!self.active.is_empty(), "a reader is active");
        if self.nsr != readers.len() {
            self.vcf_buf.resize_with(readers.len(), VecDeque::new);
            self.nsr = readers.len();
            self.chr = None;
        }
        if self.active.len() == 1 {
            if readers.len() > 1 {
                has_line.fill(false);
            }
            let a = self.active[0];
            let reader = &mut readers[a];
            super::rotate_left1(&mut reader.buffer[..=reader.nbuffer]);
            reader.nbuffer -= 1;
            has_line[a] = true;
            return 1;
        }
        if self.chr.as_deref() != Some(chr) || self.pos != min_pos {
            self.set(readers, chr, min_pos);
        }
        if self.vcf_buf.first().is_none_or(VecDeque::is_empty) {
            return 0;
        }
        let mut nret = 0;
        for (i, reader) in readers.iter_mut().enumerate() {
            let row = self.vcf_buf[i].pop_front().flatten();
            let found =
                row.and_then(|id| (1..=reader.nbuffer).find(|&j| rec_id(&reader.buffer[j]) == id));
            match found {
                Some(j) => {
                    let n = reader.nbuffer;
                    reader.buffer[..=j].rotate_right(1);
                    super::rotate_left1(&mut reader.buffer[1..=n]);
                    reader.nbuffer -= 1;
                    nret += 1;
                    has_line[i] = true;
                }
                // HTSlib asserts that the record is still buffered.
                None => has_line[i] = false,
            }
        }
        nret
    }
}

/// Converts HTSlib's variant types to the pairing types.
fn sr_type(ty: u32) -> u32 {
    let ty = VariantType::from_bits_truncate(ty);
    if ty.is_empty() {
        return SR_REF;
    }
    let mut out = 0;
    if ty.intersects(VariantType::SNP | VariantType::MNP) {
        out |= SR_SNP;
    }
    if ty.contains(VariantType::INDEL) {
        out |= SR_INDEL;
    }
    if ty.contains(VariantType::OTHER) {
        out |= SR_OTHER;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scores_follow_the_pairing_logic() {
        let mut s = Sorter {
            pair: Pairing::BOTH,
            ..Sorter::default()
        };
        s.init_scores();
        assert_eq!(s.score[score_index(SR_SNP, SR_SNP)], 3);
        assert_eq!(s.score[score_index(SR_SNP, SR_INDEL)], 0);
        assert_eq!(s.score[score_index(SR_SNP | SR_INDEL, SR_INDEL)], 3);
        let mut s = Sorter {
            pair: Pairing::ANY,
            ..Sorter::default()
        };
        s.init_scores();
        assert_eq!(s.score[score_index(SR_SNP, SR_INDEL)], 1);
        assert_eq!(s.score[score_index(SR_REF, SR_SNP)], 2);
    }

    #[test]
    fn allele_lists_compare_as_sets() {
        let var = |t: &str| Var {
            text: t.as_bytes().to_vec(),
            nalt: t.split(',').count() as i64,
            ..Var::default()
        };
        assert!(multi_is_exact(&var("A>C,A>G"), &var("A>g,A>C")));
        assert!(!multi_is_exact(&var("A>C,A>G"), &var("A>C,A>T")));
        assert!(multi_is_subset(&var("A>C,A>G"), &var("A>T,A>G")));
        assert!(!multi_is_subset(&var("A>C"), &var("A>T")));
    }
}
