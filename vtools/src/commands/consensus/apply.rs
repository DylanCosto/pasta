//! Applying one variant to the buffered sequence (bcftools' `apply_variant`).
//!
//! bcftools edits the record in place (its alleles, POS, and rlen); here the alleles are copied
//! into reusable buffers and POS and rlen kept in locals, which is equivalent because an applied
//! record is never read again.

use anyhow::bail;
use pasta::variant::{Record, Unpack, VariantMatch, VariantType, record::classify};

use super::{
    Consensus, Mark, MaskWith, PICK_ALT, PICK_IUPAC, PICK_LONG, PICK_REF, PICK_SHORT,
    genotype::{
        BT_INT8, INT8_VECTOR_END, gt_allele, gt_is_missing, gt_is_phased, set_iupac_allele,
    },
};

/// Reusable buffers for [`Consensus::apply_variant`].
#[derive(Debug, Default)]
pub struct Scratch {
    /// The record's alleles, REF first.
    als: Vec<Vec<u8>>,
    /// The alleles seen in the genotypes (bcftools' `iupac_als`).
    used: Vec<bool>,
    bitmask: Vec<u8>,
    /// One sample's GT bytes.
    gt: Vec<u8>,
    alleles_str: Vec<u8>,
    /// The REF and ALT being applied.
    ref_a: Vec<u8>,
    alt_a: Vec<u8>,
}

impl Scratch {
    fn load_alleles(&mut self, rec: &Record) {
        let n = rec.allele_count();
        self.als.resize_with(n.max(self.als.len()), Vec::new);
        self.als.truncate(n);
        for (dst, src) in self.als.iter_mut().zip(rec.alleles()) {
            dst.clear();
            dst.extend_from_slice(src);
        }
    }
}

/// The allele chosen from the genotypes (bcftools' `ialt`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Choice {
    /// Apply this allele.
    Allele(usize),
    /// A missing genotype with `--missing`: apply the missing character (`ialt == -1`).
    Missing,
    /// Leave the record out.
    Skip,
}

impl Choice {
    /// The choice for a missing genotype: [`Choice::Missing`] with `--missing`, else skip.
    fn missing(missing: Option<u8>) -> Self {
        if missing.is_some() {
            Self::Missing
        } else {
            Self::Skip
        }
    }

    /// The choice for bcftools' signed allele index (`-1` for missing).
    fn from_index(ialt: i64) -> Self {
        usize::try_from(ialt).map_or(Self::Missing, Self::Allele)
    }
}

/// Which of bcftools' three ways of choosing from one sample's GT applies.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Action {
    UseHap,
    UseIupac,
    PickOne,
}

/// `strncasecmp(a, b, n) == 0` on C strings held without their terminating NUL.
fn strncasecmp_eq(a: &[u8], b: &[u8], n: i64) -> bool {
    let n = usize::try_from(n).unwrap_or(usize::MAX);
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0).to_ascii_lowercase();
        let y = b.get(i).copied().unwrap_or(0).to_ascii_lowercase();
        if x != y {
            return false;
        }
        if x == 0 {
            break;
        }
    }
    true
}

/// bcftools' `mark_del`: `out` (ALT, or REF for a symbolic deletion) cut or padded with `mark`
/// to `rlen`.
fn mark_del(out: &mut Vec<u8>, rlen: i64, mark: u8) {
    let rlen = usize::try_from(rlen).unwrap_or(0);
    out.truncate(rlen);
    out.resize(rlen, mark);
}

/// bcftools' `mark_ins`: highlights the bases of `alt` beyond the length of `ref_a`.
fn mark_ins(ref_a: &[u8], alt: &mut [u8], mark: Mark) {
    let Some(tail) = alt.get_mut(ref_a.len()..) else {
        return;
    };
    match mark {
        Mark::Lower => tail.make_ascii_lowercase(),
        Mark::Upper => tail.make_ascii_uppercase(),
        Mark::Char(c) => tail.fill(c),
    }
}

/// bcftools' `mark_snv`: highlights the bases of `alt` that differ from `ref_a`.
fn mark_snv(ref_a: &[u8], alt: &mut [u8], mark: Mark) {
    for (a, &r) in alt.iter_mut().zip(ref_a) {
        if r.eq_ignore_ascii_case(a) {
            continue;
        }
        *a = match mark {
            Mark::Lower => a.to_ascii_lowercase(),
            Mark::Upper => a.to_ascii_uppercase(),
            Mark::Char(c) => c,
        };
    }
}

impl Consensus {
    fn seqname(&self, rec: &Record) -> String {
        String::from_utf8_lossy(rec.chrom_safe(&self.header)).into_owned()
    }

    /// Finds FORMAT/GT and checks its type, `None` when the record has none.
    /// FORMAT must be unpacked.
    fn gt_field<'r>(
        &self,
        rec: &'r Record,
    ) -> anyhow::Result<Option<&'r pasta::variant::FormatField>> {
        let Some(gt_id) = self.gt_id else {
            return Ok(None);
        };
        let Some(fmt) = rec
            .format_fields()
            .iter()
            .find(|f| f.id() == gt_id && !f.is_removed())
        else {
            return Ok(None);
        };
        if fmt.type_code() != BT_INT8 {
            bail!(
                "Todo: GT field represented with BCF_BT_INT8, too many alleles at {}:{}?",
                self.seqname(rec),
                rec.pos() + 1
            );
        }
        Ok(Some(fmt))
    }

    /// Marks the alleles of one sample's GT in `used` (bcftools' `iupac_add_gt`): whether any
    /// allele was set, or `None` for an allele beyond the record's.
    fn iupac_add_gt(used: &mut [bool], gt: &[u8]) -> Option<bool> {
        let mut is_set = false;
        for &v in gt {
            if gt_is_missing(v) {
                continue;
            }
            if v == INT8_VECTOR_END {
                break;
            }
            let ial = gt_allele(v);
            *usize::try_from(ial).ok().and_then(|i| used.get_mut(i))? = true;
            is_set = true;
        }
        Some(is_set)
    }

    /// IUPAC codes from the genotypes of several samples (bcftools' `args->iupac_GTs`).
    fn choose_iupac_samples(&mut self, rec: &mut Record) -> anyhow::Result<Choice> {
        let samples = self.iupac_samples.take().unwrap_or_default();
        let result = self.iupac_from_samples(rec, &samples);
        self.iupac_samples = Some(samples);
        result
    }

    fn iupac_from_samples(
        &mut self,
        rec: &mut Record,
        samples: &[usize],
    ) -> anyhow::Result<Choice> {
        rec.unpack(Unpack::FMT);
        // The GT bytes seen in the samples; the alleles follow from them.
        let mut seen = [false; 256];
        {
            let Some(fmt) = self.gt_field(rec)? else {
                return Ok(Choice::Skip);
            };
            let mark = |seen: &mut [bool; 256], gt: &[u8]| {
                for &v in gt {
                    if v == INT8_VECTOR_END {
                        break;
                    }
                    seen[usize::from(v)] = true;
                }
            };
            let data = rec.format_data(fmt);
            let size = fmt.size();
            if size > 0 && samples.len() * size == data.len() {
                // All samples.
                for gt in data.chunks_exact(size) {
                    mark(&mut seen, gt);
                }
            } else {
                for &i in samples {
                    mark(&mut seen, &data[i * size..(i + 1) * size]);
                }
            }
        }
        let mut used = std::mem::take(&mut self.scratch.used);
        used.clear();
        used.resize(rec.allele_count(), false);
        let mut is_set = false;
        let mut bad = false;
        for (v, _) in seen
            .iter()
            .enumerate()
            .filter(|(v, s)| **s && !gt_is_missing(*v as u8))
        {
            match usize::try_from(gt_allele(v as u8))
                .ok()
                .and_then(|i| used.get_mut(i))
            {
                Some(slot) => {
                    *slot = true;
                    is_set = true;
                }
                None => bad = true,
            }
        }
        if bad {
            self.scratch.used = used;
            bail!(
                "Invalid VCF, too few ALT alleles at {}:{}",
                self.seqname(rec),
                rec.pos() + 1
            );
        }
        let choice = if is_set {
            set_iupac_allele(
                &mut self.scratch.als,
                Some(&used),
                &mut self.scratch.bitmask,
            )
            .map_or(Choice::Missing, Choice::Allele)
        } else {
            Choice::missing(self.opts.missing)
        };
        self.scratch.used = used;
        Ok(choice)
    }

    /// The allele of the `-H` sample (bcftools' `args->isample >= 0` branch).
    fn choose_from_sample(&mut self, rec: &mut Record, isample: usize) -> anyhow::Result<Choice> {
        let n_allele = rec.allele_count();
        let mut gt = std::mem::take(&mut self.scratch.gt);
        gt.clear();
        rec.unpack(Unpack::FMT);
        match self.gt_field(rec)? {
            Some(fmt) => gt.extend_from_slice(rec.format_sample_data(fmt, isample)),
            None => {
                self.scratch.gt = gt;
                return Ok(Choice::Skip);
            }
        }
        let result = self.choose_from_gt(rec, &gt);
        self.scratch.gt = gt;
        match result? {
            Choice::Allele(0) => {
                // The reference allele.
                if self.opts.absent.is_some() {
                    self.freeze_ref(rec.pos(), rec.rlen());
                }
                Ok(Choice::Skip)
            }
            Choice::Allele(ialt) if n_allele <= ialt => bail!(
                "Broken VCF, too few alts at {}:{}",
                self.seqname(rec),
                rec.pos() + 1
            ),
            choice => Ok(choice),
        }
    }

    fn choose_from_gt(&mut self, rec: &Record, gt: &[u8]) -> anyhow::Result<Choice> {
        let n = gt.len();
        if n == 0 {
            return Ok(Choice::Skip);
        }
        let haplotype = self.opts.haplotype;
        let mut action = Action::UseHap;
        if self.opts.allele == PICK_IUPAC {
            if haplotype == 0 {
                action = Action::UseIupac;
            }
            if !gt_is_phased(gt[0]) && !gt_is_phased(gt[n - 1]) {
                action = Action::UseIupac;
            }
        } else if self.opts.output_iupac {
            action = Action::UseIupac;
        } else if haplotype == 0 {
            action = Action::PickOne;
        }
        match action {
            Action::UseHap => {
                if haplotype > n as i64 {
                    if gt_is_missing(gt[n - 1]) || gt_is_missing(gt[0]) {
                        return Ok(Choice::missing(self.opts.missing));
                    }
                    if !self.warned_haplotype {
                        eprintln!(
                            "Can't apply {haplotype}-th haplotype at {}:{}. (This warning is printed only once.)",
                            self.seqname(rec),
                            rec.pos() + 1
                        );
                        self.warned_haplotype = true;
                    }
                    return Ok(Choice::Skip);
                }
                let v = gt[(haplotype - 1) as usize];
                if gt_is_missing(v) || v == INT8_VECTOR_END {
                    return Ok(Choice::missing(self.opts.missing));
                }
                Ok(Choice::from_index(gt_allele(v)))
            }
            Action::UseIupac => {
                let mut used = std::mem::take(&mut self.scratch.used);
                used.clear();
                used.resize(rec.allele_count(), false);
                let set = Self::iupac_add_gt(&mut used, gt);
                let choice = match set {
                    None => {
                        self.scratch.used = used;
                        bail!(
                            "Invalid VCF, too few ALT alleles at {}:{}",
                            self.seqname(rec),
                            rec.pos() + 1
                        );
                    }
                    Some(false) => Choice::missing(self.opts.missing),
                    Some(true) => set_iupac_allele(
                        &mut self.scratch.als,
                        Some(&used),
                        &mut self.scratch.bitmask,
                    )
                    .map_or(Choice::Missing, Choice::Allele),
                };
                self.scratch.used = used;
                Ok(choice)
            }
            Action::PickOne => self.pick_one(rec, gt),
        }
    }

    /// `-H R/A/LR/LA/SR/SA` (bcftools' `pick_one`).
    fn pick_one(&self, rec: &Record, gt: &[u8]) -> anyhow::Result<Choice> {
        let mut ialt: i64 = 1;
        let mut is_hom = true;
        for (i, &v) in gt.iter().enumerate() {
            if gt_is_missing(v) {
                if self.opts.missing.is_none() {
                    return Ok(Choice::Skip);
                }
                ialt = -1;
                break;
            }
            if v == INT8_VECTOR_END {
                break;
            }
            ialt = gt_allele(v);
            if i > 0 && ialt != gt_allele(gt[i - 1]) {
                is_hom = false;
                break;
            }
        }
        if is_hom {
            return Ok(Choice::from_index(ialt));
        }
        let allele = self.opts.allele;
        // Whether allele `jalt` replaces `ialt` as -H R or -H A prefer.
        let prefer = |ialt: i64, jalt: i64| {
            (allele & PICK_REF != 0 && jalt == 0) || (allele & PICK_ALT != 0 && ialt == 0)
        };
        let mut prev_len = 0i64;
        for (i, &v) in gt.iter().enumerate() {
            if v == INT8_VECTOR_END {
                break;
            }
            let jalt = gt_allele(v);
            if rec.allele_count() as i64 <= jalt {
                bail!(
                    "Broken VCF, too few alts at {}:{}",
                    self.seqname(rec),
                    rec.pos() + 1
                );
            }
            if allele & (PICK_LONG | PICK_SHORT) != 0 {
                let len = if jalt == 0 {
                    rec.rlen()
                } else {
                    usize::try_from(jalt).map_or(0, |j| rec.allele(j).len() as i64)
                };
                let take = if i == 0 {
                    true
                } else if len == prev_len {
                    prefer(ialt, jalt)
                } else {
                    (allele & PICK_LONG != 0 && len > prev_len)
                        || (allele & PICK_SHORT != 0 && len < prev_len)
                };
                if take {
                    (ialt, prev_len) = (jalt, len);
                }
            } else if prefer(ialt, jalt) {
                ialt = jalt;
            }
        }
        Ok(Choice::from_index(ialt))
    }

    /// Whether a variant overlaps a mask that replaces bases (`MASK_SKIP`).
    fn masked(&mut self, rec: &Record) -> bool {
        let start = super::clamp_u32(rec.pos());
        let end = super::clamp_u32(rec.pos() + rec.rlen() - 1);
        let chr = &self.chr;
        self.masks
            .iter_mut()
            .any(|m| matches!(m.with, MaskWith::Char(_)) && m.idx.overlaps(chr, start, end))
    }

    /// Applies a record to the buffer (bcftools' `apply_variant`).
    pub(super) fn apply_variant(&mut self, rec: &mut Record) -> anyhow::Result<()> {
        if self.opts.absent.is_some() {
            self.apply_absent(rec.pos());
        }
        rec.unpack(Unpack::STR);
        if rec.allele_count() == 1 && self.opts.missing.is_none() && self.opts.absent.is_none() {
            return Ok(());
        }
        if !self.masks.is_empty() && self.masked(rec) {
            return Ok(());
        }
        self.scratch.load_alleles(rec);

        let mut choice = Choice::Allele(1);
        if self.iupac_samples.is_some() {
            choice = self.choose_iupac_samples(rec)?;
        } else if let Some(isample) = self.isample {
            choice = self.choose_from_sample(rec, isample)?;
        } else if self.opts.output_iupac && rec.allele_count() > 1 {
            // bcftools writes the codes into the longest ALT but applies the first one (its
            // index is a shadowed variable). The two differ only for multiallelic indels, where
            // the coded allele would merge misaligned repeats, so this keeps bcftools' choice.
            set_iupac_allele(&mut self.scratch.als, None, &mut self.scratch.bitmask);
        }

        let ialt = match choice {
            Choice::Skip => return Ok(()),
            Choice::Allele(_) if rec.allele_count() == 1 => {
                // A non-missing reference site.
                if self.opts.absent.is_some() {
                    self.freeze_ref(rec.pos(), rec.rlen());
                }
                return Ok(());
            }
            Choice::Allele(i) => i,
            Choice::Missing => return self.apply_missing(rec),
        };
        self.apply_allele(rec, ialt)
    }

    /// Applies a record with a missing genotype as REF replaced by the `--missing` character
    /// (bcftools' `ialt == -1` branch).
    fn apply_missing(&mut self, rec: &mut Record) -> anyhow::Result<()> {
        let mut ref_a = std::mem::take(&mut self.scratch.ref_a);
        let mut alt_a = std::mem::take(&mut self.scratch.alt_a);
        let (pos, rlen) = self.missing_alleles(rec, &mut ref_a, &mut alt_a);
        let result = self.apply_alleles(rec, 1, pos, rlen, &mut ref_a, &mut alt_a);
        self.scratch.ref_a = ref_a;
        self.scratch.alt_a = alt_a;
        result
    }

    /// Sets `ref_a` and `alt_a` for a record with a missing genotype: REF and the `--missing`
    /// character for every base of a reference block, else the first REF base and the
    /// character. Returns the POS and rlen they replace.
    ///
    /// bcftools copies the whole block from its buffer, reading before it for a block that
    /// starts before a `>chr:from-to` region and past it for one that ends after the region
    /// or sequence, and crashes. A block is clipped here to the buffer instead, as
    /// [`Self::apply_alleles`] trims other variants, so only the bases it replaces are copied.
    fn missing_alleles(
        &mut self,
        rec: &mut Record,
        ref_a: &mut Vec<u8>,
        alt_a: &mut Vec<u8>,
    ) -> (i64, i64) {
        let missing = self.opts.missing.unwrap_or(b'N');
        ref_a.clear();
        alt_a.clear();
        let ref_only = rec
            .has_variant_types(VariantType::empty(), VariantMatch::Exact)
            .0;
        if rec.rlen() <= 1 || !ref_only {
            let mut s = std::mem::take(&mut self.scratch.alleles_str);
            s.clear();
            s.push(rec.allele(0).first().copied().unwrap_or(0));
            s.push(b',');
            s.push(missing);
            // This also sets rlen from the new REF, as bcftools' update does.
            rec.update_alleles_str(&self.header, &s);
            self.scratch.alleles_str = s;
            ref_a.extend_from_slice(rec.allele(0));
            alt_a.extend_from_slice(rec.allele(1));
            return (rec.pos(), rec.rlen());
        }
        let (mut pos, mut rlen) = (rec.pos(), rec.rlen());
        // The position of the block within the modified sequence.
        let idx = pos - self.fa_ori_pos + self.fa_mod_off;
        if idx < 0 && pos > self.fa_frz_pos && rlen > -idx {
            // The block starts before the buffer: keep the part in it, as apply_alleles would.
            // (A block at or before the protected position is skipped as overlapping.)
            pos -= idx;
            rlen += idx;
        }
        // The block ends past the buffer; at least one base, for one that starts past it.
        rlen = rlen.min((self.buffer_end() - pos).max(1));
        let n = usize::try_from(rlen).unwrap_or(0);
        let start = usize::try_from(pos - self.fa_ori_pos + self.fa_mod_off)
            .unwrap_or(0)
            .min(self.fa_buf.len());
        let end = (start + n).min(self.fa_buf.len());
        ref_a.extend_from_slice(&self.fa_buf[start..end]);
        ref_a.resize(n, b'N');
        alt_a.resize(n, missing);
        (pos, rlen)
    }

    /// Applies allele `ialt` of the record (the rest of bcftools' `apply_variant`).
    fn apply_allele(&mut self, rec: &Record, ialt: usize) -> anyhow::Result<()> {
        let mut ref_a = std::mem::take(&mut self.scratch.ref_a);
        let mut alt_a = std::mem::take(&mut self.scratch.alt_a);
        ref_a.clone_from(&self.scratch.als[0]);
        alt_a.clone_from(&self.scratch.als[ialt]);
        let result = self.apply_alleles(rec, ialt, rec.pos(), rec.rlen(), &mut ref_a, &mut alt_a);
        self.scratch.ref_a = ref_a;
        self.scratch.alt_a = alt_a;
        result
    }

    /// [`Self::apply_allele`] with REF and ALT in `ref_a` and `alt_a` replacing `rlen` bases from
    /// `pos`, which it edits as bcftools edits the record's alleles, POS, and rlen.
    #[expect(
        clippy::too_many_lines,
        reason = "a line-for-line port of bcftools' apply_variant"
    )]
    fn apply_alleles(
        &mut self,
        rec: &Record,
        ialt: usize,
        mut pos: i64,
        mut rlen: i64,
        ref_a: &mut Vec<u8>,
        alt_a: &mut Vec<u8>,
    ) -> anyhow::Result<()> {
        let name = || String::from_utf8_lossy(rec.chrom_safe(&self.header)).into_owned();

        // For some variant types POS+REF refer to the base before the event.
        let var = if ialt == 0 {
            pasta::variant::Variant::default()
        } else {
            classify(ref_a, alt_a)
        };
        let mut trim_beg = false;
        let mut var_len = i64::from(var.n);
        if var.ty.contains(VariantType::INDEL) {
            // Normally an indel starts a base after POS, but not when the first REF base is
            // deleted.
            trim_beg = ref_a.first() == alt_a.first();
        } else if var.ty.contains(VariantType::OTHER) && alt_a.eq_ignore_ascii_case(b"<DEL>") {
            trim_beg = true;
            var_len = 1 - rlen;
        } else if var.ty.contains(VariantType::OTHER)
            && alt_a.len() >= 4
            && alt_a[..4].eq_ignore_ascii_case(b"<INS")
        {
            trim_beg = true;
        }

        // Overlapping variant?
        if ialt == 0 && pos <= self.fa_frz_pos && pos + rlen - 1 > self.fa_frz_pos {
            // The reference allele overlaps a previous deletion and extends beyond it: trim it.
            let ntrim = self.fa_frz_pos - pos + 1;
            let nref = ref_a.len() as i64;
            if ntrim >= nref {
                let name = name();
                bail!(
                    "Error: failed to trim overlapping variant at {name}:{}, ntrim={ntrim} >= nref={nref}. Is the VCF normalized?",
                    pos + 1
                );
            }
            pos += ntrim;
            rlen -= ntrim;
            ref_a.drain(..ntrim as usize);
            alt_a.clone_from(ref_a);
        }
        if pos <= self.fa_frz_pos {
            // Still fine for an insertion that does not follow another insertion.
            let overlap = pos < self.fa_frz_pos || !trim_beg || var_len == 0 || self.prev_is_insert;
            if overlap {
                let name = name();
                eprintln!(
                    "The site {name}:{} overlaps with another variant, skipping...",
                    pos + 1
                );
                return Ok(());
            }
        }

        let mut len_diff: i64;
        let mut alen: i64;
        // The position of the variant within the modified sequence.
        let mut idx = pos - self.fa_ori_pos + self.fa_mod_off;
        if idx < 0 {
            let back = -idx;
            if alt_a.first() == Some(&b'<') {
                pos -= idx + 1;
                rlen += idx + 1;
                idx = -1;
            } else if (ref_a.len() as i64) < back {
                let name = name();
                eprintln!(
                    "Warning: ignoring overlapping variant starting at {name}:{}",
                    pos + 1
                );
                return Ok(());
            } else if alt_a.len() as i64 > back {
                // Both alleles overlap the buffer.
                pos -= idx;
                rlen += idx;
                ref_a.drain(..back as usize);
                alt_a.drain(..back as usize);
                idx = 0;
            } else {
                // REF overlaps the buffer but ALT does not: leave one base before.
                pos -= idx + 1;
                rlen += idx + 1;
                ref_a.drain(..(back - 1) as usize);
                let last = alt_a.last().copied();
                alt_a.clear();
                alt_a.extend(last);
                idx = -1;
            }
        }
        let buf_len = self.fa_buf.len() as i64;
        if rlen > buf_len - idx {
            rlen = buf_len - idx;
            if alt_a.first() != Some(&b'<') && alt_a.len() as i64 > rlen {
                let name = name();
                eprintln!(
                    "Warning: trimming variant \"{}\" starting at {name}:{}",
                    String::from_utf8_lossy(alt_a),
                    pos + 1
                );
                let keep = usize::try_from(rlen).unwrap_or(0);
                alt_a.truncate(keep);
                if ialt == 0 {
                    // bcftools' REF and ALT are then the same string.
                    ref_a.truncate(keep);
                }
            }
        }
        // The variant is beyond the available sequence.
        if idx > 0 && idx >= buf_len {
            return Ok(());
        }

        if alt_a.first() == Some(&b'<') {
            let is_del = alt_a.eq_ignore_ascii_case(b"<DEL>");
            if !is_del
                && !alt_a.eq_ignore_ascii_case(b"<*>")
                && !alt_a.eq_ignore_ascii_case(b"<NON_REF>")
            {
                let name = name();
                bail!(
                    "Symbolic alleles other than <DEL>, <*> or <NON_REF> are currently not supported, e.g. \"{}\" at {name}:{}.\nPlease use filtering expressions to exclude such sites, for example by running with: -e 'ALT~\"<.*>\"'",
                    String::from_utf8_lossy(alt_a),
                    pos + 1
                );
            }
            if !is_del {
                // <*> or <NON_REF>: a gVCF block, evidence for the reference.
                self.freeze_ref(pos, rlen);
                return Ok(());
            }
            if ref_a.len() > 1 && !self.warned_multibase_del {
                let name = name();
                eprintln!(
                    "Warning: one REF base is expected with <DEL>, assuming the actual deletion starts at POS+1 at {name}:{}.\n         (This warning is printed only once.)",
                    pos + 1
                );
                self.warned_multibase_del = true;
            }
            if let Some(m) = self.opts.mark_del {
                alt_a.clone_from(ref_a);
                mark_del(alt_a, rlen, m);
                alen = rlen;
                len_diff = 0;
            } else {
                len_diff = 1 - rlen;
                // The first REF base precedes the event.
                alt_a.clone_from(ref_a);
                alen = 1;
            }
        } else {
            if idx >= 0 && !strncasecmp_eq(ref_a, &self.fa_buf[idx as usize..], rlen) {
                self.check_ref_after_deletion(ref_a, alt_a, idx, pos, rlen, &name())?;
            }
            alen = alt_a.len() as i64;
            len_diff = alen - rlen;
            if let Some(m) = self.opts.mark_del
                && len_diff < 0
            {
                mark_del(alt_a, rlen, m);
                alen = rlen;
                len_diff = 0;
            }
        }

        let alen_u = usize::try_from(alen).unwrap_or(0).min(alt_a.len());
        let case_base = self.fa_buf.get(idx.max(0) as usize).copied().unwrap_or(0);
        if case_base.to_ascii_uppercase() == case_base {
            alt_a[..alen_u].make_ascii_uppercase();
        } else {
            alt_a[..alen_u].make_ascii_lowercase();
        }
        if let Some(m) = self.opts.mark_ins
            && len_diff > 0
        {
            mark_ins(ref_a, alt_a, m);
        }
        if let Some(m) = self.opts.mark_snv {
            mark_snv(ref_a, alt_a, m);
        }

        // The buffer changes: forget what --absent filled.
        self.absent_filled = 0;
        if len_diff <= 0 {
            // A deletion or a same-length substitution.
            if idx + rlen < 0 || buf_len < idx + rlen {
                let name = name();
                bail!(
                    "Error: fasta buffer too short at {name}:{} (buf length {buf_len}, need {}). Is the reference fasta correct?",
                    pos + 1,
                    idx + rlen
                );
            }
            self.prev_base = self.byte_at(idx + rlen - 1);
            self.prev_base_pos = pos + rlen - 1;
            self.prev_is_insert = false;
            self.fa_frz_mod = idx + alen;
            for i in i64::from(trim_beg)..alen {
                self.set_byte(idx + i, alt_a[i as usize]);
            }
            if len_diff != 0 {
                let src = (idx + rlen).max(0) as usize;
                let dst = (idx + alen).max(0) as usize;
                self.fa_buf.copy_within(src.., dst);
                self.fa_buf.truncate((buf_len + len_diff) as usize);
            }
        } else {
            // An insertion.
            self.prev_is_insert = true;
            self.prev_base_pos = pos;
            let src = (idx + rlen).max(0) as usize;
            self.fa_buf.resize((buf_len + len_diff) as usize, 0);
            self.fa_buf
                .copy_within(src..buf_len as usize, src + len_diff as usize);
            // Bases the insertion leaves unchanged must not overwrite a preceding variant (as
            // in `C>T` followed by `C>CAA` at one position, which gives `TAA`).
            let mut ibeg = 0i64;
            while ibeg < alen
                && ref_a.get(ibeg as usize) == Some(&alt_a[ibeg as usize])
                && pos + ibeg <= self.prev_base_pos
            {
                ibeg += 1;
            }
            for i in ibeg..alen {
                self.set_byte(idx + i, alt_a[i as usize]);
            }
            self.fa_frz_mod = idx + alen - ibeg + 1;
        }
        if let Some(chain) = &mut self.chain
            && len_diff != 0
        {
            if strncasecmp_eq(ref_a, alt_a, 1) {
                // Indels include the base before the event: the block extends by one.
                chain.push_gap(pos + 1, rlen - 1, pos + 1 + self.fa_mod_off, alen - 1);
            } else {
                chain.push_gap(pos, rlen, pos + self.fa_mod_off, alen);
            }
        }
        self.fa_mod_off += len_diff;
        self.fa_frz_pos = pos + rlen - 1;
        self.napplied += 1;
        Ok(())
    }

    /// Accepts a REF that disagrees with the buffer only because a preceding deletion removed
    /// its first base (`AAC>A` then `C>CAA`), and fails otherwise.
    fn check_ref_after_deletion(
        &self,
        ref_a: &[u8],
        alt_a: &[u8],
        idx: i64,
        pos: i64,
        rlen: i64,
        name: &str,
    ) -> anyhow::Result<()> {
        let i = idx as usize;
        if self.prev_base_pos == pos
            && ref_a.first().map(u8::to_ascii_uppercase)
                == Some(self.prev_base.to_ascii_uppercase())
            && (rlen == 1
                || strncasecmp_eq(
                    ref_a.get(1..).unwrap_or_default(),
                    self.fa_buf.get(i + 1..).unwrap_or_default(),
                    rlen - 1,
                ))
        {
            return Ok(());
        }
        let buf = &self.fa_buf;
        let n = usize::try_from(rlen).unwrap_or(0);
        let (shown, next, rest) = if buf.len() - i > n {
            (&buf[i..i + n], buf[i + n], &buf[i + n + 1..])
        } else {
            (&buf[i..], b' ', &buf[buf.len()..])
        };
        bail!(
            "The fasta sequence does not match the REF allele at {name}:{}:\n   REF .vcf: [{}]\n   ALT .vcf: [{}]\n   REF .fa : [{}]{}{}",
            pos + 1,
            String::from_utf8_lossy(ref_a),
            String::from_utf8_lossy(alt_a),
            String::from_utf8_lossy(shown),
            next as char,
            String::from_utf8_lossy(rest)
        );
    }

    fn byte_at(&self, i: i64) -> u8 {
        usize::try_from(i)
            .ok()
            .and_then(|i| self.fa_buf.get(i))
            .copied()
            .unwrap_or(0)
    }

    /// Writes a buffer byte; bcftools writes before the buffer for some overlapping deletions,
    /// which is skipped here.
    fn set_byte(&mut self, i: i64, c: u8) {
        if let Some(b) = usize::try_from(i).ok().and_then(|i| self.fa_buf.get_mut(i)) {
            *b = c;
        }
    }
}
