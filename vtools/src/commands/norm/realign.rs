//! REF checks and fixes (`-c`), and left- or right-alignment of indels against the reference
//! (bcftools' `fix_ref`, `fix_dup_alt`, `realign`, `realign_left`, `realign_right`).

use anyhow::bail;
use pasta::variant::{Genotype, Header, Record, Unpack, VariantType};

use super::{CHECK_REF_EXIT, CHECK_REF_WARN, Normalizer, Realigned, eq_nocase, is_sequence_allele};

/// Buffers for realignment.
#[derive(Debug, Default)]
pub(super) struct Seqs {
    /// The sequences being trimmed and padded (bcftools' `tmp_als`).
    als: Vec<Vec<u8>>,
    /// The original symbolic alleles (bcftools' `tmp_sym`).
    sym: Vec<Vec<u8>>,
    /// Editable copies of the alleles.
    alleles: Vec<Vec<u8>>,
    reference: Vec<u8>,
    pad: Vec<u8>,
    joined: Vec<u8>,
}

/// bcftools' `iupac2first`: the first base an IUPAC code stands for (case kept), the base
/// itself for ACGTN, and 0 for anything else.
fn iupac2first(c: u8) -> u8 {
    match c {
        b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' => c,
        b'R' | b'W' | b'M' | b'D' | b'H' | b'V' => b'A',
        b'r' | b'w' | b'm' | b'd' | b'h' | b'v' => b'a',
        b'Y' | b'S' | b'B' => b'C',
        b'y' | b's' | b'b' => b'c',
        b'K' => b'G',
        b'k' => b'g',
        _ => 0,
    }
}

/// Replaces IUPAC ambiguity codes with their first base; returns how many changed (bcftools'
/// `replace_iupac_codes`).
fn replace_iupac_codes(seq: &mut [u8]) -> usize {
    let mut n = 0;
    for c in seq {
        let f = iupac2first(*c);
        if f != 0 && f != *c {
            *c = f;
            n += 1;
        }
    }
    n
}

/// Whether a sequence has characters other than ACGTN (bcftools' `has_non_acgtn`).
fn has_non_acgtn(seq: &[u8]) -> bool {
    seq.iter()
        .take_while(|&&c| c != 0)
        .any(|&c| iupac2first(c) != c)
}

/// bcftools' `nt_to_upper`.
fn nt_to_upper(seq: &mut [u8]) {
    for c in seq {
        if *c >= 97 {
            *c -= 32;
        }
    }
}

impl Normalizer {
    /// Fetches reference bases `beg..=end` of the record's contig into `seqs.reference`.
    /// Returns `false` when the contig is not in the reference.
    fn fetch(&mut self, hdr: &Header, rid: i32, beg: i64, end: i64) -> anyhow::Result<bool> {
        let fai = self.fai.as_mut().expect("realignment needs a reference");
        let name = hdr.contig_name(rid).unwrap_or(b"(unknown)");
        fai.fetch(name, beg, end, &mut self.seqs.reference)
    }

    /// Copies the record's alleles into `seqs.alleles`.
    fn load_alleles(&mut self, line: &mut Record) {
        line.unpack(Unpack::STR);
        let n = line.allele_count();
        let alleles = &mut self.seqs.alleles;
        if alleles.len() < n {
            alleles.resize_with(n, Vec::new);
        }
        for (i, dst) in alleles.iter_mut().take(n).enumerate() {
            dst.clear();
            dst.extend_from_slice(line.allele(i));
        }
    }

    /// Checks REF against the reference and fixes it (`-c s`): IUPAC codes, missing or `N`
    /// bases, and swapped REF/ALT. Returns whether the record changed (bcftools' `fix_ref`).
    pub(super) fn fix_ref(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        line: &mut Record,
    ) -> anyhow::Result<bool> {
        self.load_alleles(line);
        let n = line.allele_count();
        let reflen = self.seqs.alleles[0].len();
        let seqname = String::from_utf8_lossy(line.chrom_safe(hdr)).into_owned();
        let pos = line.pos();
        if !self.fetch(hdr, line.rid(), pos, pos + reflen as i64 - 1)? {
            bail!(
                "Error: unable to fetch {reflen} bp from the reference sequence at {seqname}:{}",
                pos + 1
            );
        }
        let got = self.seqs.reference.len();
        if got != reflen {
            bail!(
                "Error: unable to fetch {reflen} bp from the reference sequence at {seqname}:{} .. requested {reflen} bp, got {got}",
                pos + 1
            );
        }
        self.stats.ref_total += 1;
        let seqs = &mut self.seqs;
        replace_iupac_codes(&mut seqs.reference);
        let reference = &seqs.reference;
        let alleles = &mut seqs.alleles;

        let mut ret = false;
        for a in alleles.iter_mut().take(n) {
            if is_sequence_allele(a) && replace_iupac_codes(a) > 0 {
                ret = true;
            }
        }

        // Is the REF allele missing?
        if reflen == 1 && alleles[0][0] == b'.' {
            alleles[0][0] = reference[0];
            self.stats.ref_set += 1;
            line.update_alleles(out_hdr, &alleles[..n]);
            return Ok(true);
        }

        // If REF matches, only the IUPAC cleanup may need writing.
        if eq_nocase(&alleles[0], reference) {
            if ret {
                self.stats.ref_set += 1;
                line.update_alleles(out_hdr, &alleles[..n]);
            }
            return Ok(ret);
        }

        // Resolve N placeholders in REF before deciding whether REF and ALT are swapped.
        let mut ref_n_changed = false;
        for j in 0..reflen {
            if !matches!(alleles[0][j], b'N' | b'n') || matches!(reference[j], b'N' | b'n') {
                continue;
            }
            alleles[0][j] = reference[j];
            ref_n_changed = true;
            ret = true;
            for a in alleles.iter_mut().take(n).skip(1) {
                if !is_sequence_allele(a) || a.len() <= j || !matches!(a[j], b'N' | b'n') {
                    continue;
                }
                a[j] = reference[j];
            }
        }
        if ref_n_changed && eq_nocase(&alleles[0], reference) {
            self.stats.ref_set += 1;
            line.update_alleles(out_hdr, &alleles[..n]);
            return Ok(ret);
        }

        // Is REF/ALT simply swapped? Only same-length sequence alleles count.
        let ialt = (1..n)
            .find(|&i| {
                let a = &alleles[i];
                is_sequence_allele(a) && a.len() == reflen && eq_nocase(a, reference)
            })
            .unwrap_or(n);
        if ialt == n {
            // No ALT is evidence of a swap: set REF from the reference and make REF bases
            // copied into the ALTs follow.
            let old_ref = std::mem::take(&mut alleles[0]);
            for a in alleles.iter_mut().take(n).skip(1) {
                if !is_sequence_allele(a) {
                    continue;
                }
                let m = a.len().min(reflen);
                for j in 0..m {
                    if a[j].eq_ignore_ascii_case(&old_ref[j]) {
                        a[j] = reference[j];
                    }
                }
            }
            alleles[0] = old_ref;
            alleles[0].copy_from_slice(reference);
            self.stats.ref_set += 1;
            line.update_alleles(out_hdr, &alleles[..n]);
            return Ok(true);
        }

        // One of the ALTs matches the reference: assume a simple swap.
        let joined = &mut seqs.joined;
        joined.clear();
        joined.extend_from_slice(&alleles[ialt]);
        for (j, a) in alleles.iter().enumerate().take(n).skip(1) {
            joined.push(b',');
            joined.extend_from_slice(if j == ialt { &alleles[0] } else { a });
        }
        line.update_alleles_str(out_hdr, joined);
        self.stats.ref_swap += 1;

        // Swap genotypes.
        let gts = &mut self.scratch.i.a;
        let mut ni = 0;
        if line.get_genotypes(hdr, gts).is_ok_and(|n| n > 0) {
            let ialt = ialt as i32;
            for g in gts.iter_mut() {
                if *g == Genotype::unphased(0).0 {
                    *g = Genotype::unphased(ialt).0;
                    ni += 1;
                } else if *g == Genotype::phased(0).0 {
                    *g = Genotype::phased(ialt).0;
                    ni += 1;
                } else if *g == Genotype::unphased(ialt).0 {
                    *g = Genotype::unphased(0).0;
                } else if *g == Genotype::phased(ialt).0 {
                    *g = Genotype::phased(0).0;
                }
            }
            let _ = line.update_genotypes(out_hdr, gts);
        }

        // Update AC.
        let ac = &mut self.scratch.i.b;
        if let Ok(nac) = line.get_info_i32(hdr, b"AC", ac)
            && nac > 0
            && ialt <= nac
        {
            ac[ialt - 1] = ni;
            let _ = line.update_info_i32(out_hdr, b"AC", ac);
        }
        Ok(true)
    }

    /// Removes ALT alleles equal to REF, remapping genotypes (bcftools' `fix_dup_alt`).
    pub(super) fn fix_dup_alt(&mut self, hdr: &Header, out_hdr: &Header, line: &mut Record) {
        self.load_alleles(line);
        let nals_ori = line.allele_count();
        let mut map = vec![0usize; nals_ori];
        let mut j = 1;
        for (i, m) in map.iter_mut().enumerate().skip(1) {
            if eq_nocase(&self.seqs.alleles[0], &self.seqs.alleles[i]) {
                *m = 0;
            } else {
                *m = j;
                j += 1;
            }
        }
        let nals = j;
        let mut k = 1;
        for (i, &m) in map.iter().enumerate().skip(1) {
            if m != 0 {
                self.seqs.alleles.swap(k, i);
                k += 1;
            }
        }
        line.update_alleles(out_hdr, &self.seqs.alleles[..nals]);

        let gts = &mut self.scratch.i.a;
        if line.get_genotypes(hdr, gts).is_err() {
            return;
        }
        let mut changed = false;
        for g in gts.iter_mut() {
            let gt = Genotype(*g);
            if gt.is_missing() || gt.is_vector_end() {
                continue;
            }
            let ial = gt.allele();
            if (ial as usize) < nals_ori && ial as usize == map[ial as usize] {
                continue;
            }
            let new = if (0..nals_ori as i32).contains(&ial) {
                map[ial as usize] as i32
            } else {
                0
            };
            *g = if gt.is_phased() {
                Genotype::phased(new).0
            } else {
                Genotype::unphased(new).0
            };
            changed = true;
        }
        if changed {
            let _ = line.update_genotypes(out_hdr, gts);
        }
    }

    /// Checks REF and left- (or right-) aligns the record (bcftools' `realign`).
    pub(super) fn realign(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        line: &mut Record,
    ) -> anyhow::Result<Realigned> {
        line.unpack(Unpack::STR);
        self.old_rec_tag_init(hdr, line);

        // Sanity check REF.
        let mut reflen = line.allele(0).len();
        let pos = line.pos();
        let rid = line.rid();
        let chrom: &[u8] = line.chrom_safe(hdr);
        let ctg = || String::from_utf8_lossy(hdr.contig_name(rid).unwrap_or(b"(null)"));
        let seqname = || String::from_utf8_lossy(chrom);
        if !self.fetch(hdr, rid, pos, pos + reflen as i64 - 1)? {
            bail!("faidx_fetch_seq failed at {}:{}", ctg(), pos + 1);
        }
        let mut nref = self.seqs.reference.len();
        nt_to_upper(&mut self.seqs.reference);
        replace_iupac_codes(&mut self.seqs.reference);
        let check_ref = self.opts.check_ref;

        // Does VCF REF contain non-standard bases?
        if has_non_acgtn(line.allele(0)) {
            let (r, v) = (
                String::from_utf8_lossy(&self.seqs.reference),
                String::from_utf8_lossy(line.allele(0)),
            );
            if check_ref == CHECK_REF_EXIT {
                bail!(
                    "Non-ACGTN reference allele at {}:{} .. REF_SEQ:'{r}' vs VCF:'{v}'",
                    seqname(),
                    pos + 1
                );
            }
            if check_ref & CHECK_REF_WARN != 0 {
                eprintln!("NON_ACGTN_REF\t{}\t{}\t{v}", seqname(), pos + 1);
            }
            return Ok(Realigned::RefMismatch);
        }
        if !eq_nocase(&self.seqs.reference, line.allele(0)) {
            let (r, v) = (
                String::from_utf8_lossy(&self.seqs.reference),
                String::from_utf8_lossy(line.allele(0)),
            );
            if check_ref == CHECK_REF_EXIT {
                bail!(
                    "Reference allele mismatch at {}:{} .. REF_SEQ:'{r}' vs VCF:'{v}'",
                    seqname(),
                    pos + 1
                );
            }
            if check_ref & CHECK_REF_WARN != 0 {
                eprintln!("REF_MISMATCH\t{}\t{}\t{v}\t{r}", seqname(), pos + 1);
            }
            return Ok(Realigned::RefMismatch);
        }

        let n = line.allele_count();
        if n == 1 {
            // A REF-only record.
            if line.rlen() > 1 {
                let first = [line.allele(0)[0]];
                line.update_alleles(out_hdr, &[first]);
            }
            return Ok(Realigned::Ok);
        }
        if line.variant_types() == VariantType::BND {
            return Ok(Realigned::Symbolic);
        }

        // Copy each allele for trimming.
        let seqs = &mut self.seqs;
        if seqs.als.len() < n {
            seqs.als.resize_with(n, Vec::new);
            seqs.sym.resize_with(n, Vec::new);
        }
        let mut symbolic_alts = true;
        for i in 0..n {
            seqs.sym[i].clear();
            let allele = line.allele(i);
            if allele.first() == Some(&b'<') {
                // Only <DEL.*> and <DUP.*> are realigned.
                let mut sv_len: i64 = 0;
                if allele.starts_with(b"<DEL") {
                    sv_len = -line.rlen();
                } else if allele.starts_with(b"<DUP") {
                    let v = &mut self.scratch.i.a;
                    if line.get_info_i32(hdr, b"SVLEN", v) == Ok(1) {
                        sv_len = i64::from(v[0]);
                    }
                }
                if sv_len == 0 {
                    return Ok(Realigned::Symbolic);
                }
                let rlen = line.rlen();
                let fai = self.fai.as_mut().expect("realignment needs a reference");
                let name = hdr.contig_name(rid).unwrap_or(b"(unknown)");
                seqs.als[i].clear();
                if sv_len < 0 {
                    // A deletion: expand REF and replace ALT ("C <DEL>" becomes "CAT C").
                    if (nref as i64) < rlen {
                        reflen = rlen as usize;
                        if !fai.fetch(name, pos, pos + reflen as i64 - 1, &mut seqs.reference)? {
                            bail!("faidx_fetch_seq failed at {}:{}", ctg(), pos + 1);
                        }
                        nref = seqs.reference.len();
                        nt_to_upper(&mut seqs.reference);
                        replace_iupac_codes(&mut seqs.reference);
                        seqs.als[0].clear();
                        seqs.als[0].extend_from_slice(&seqs.reference);
                    }
                    let first = seqs.als[0].first().copied();
                    seqs.als[i].extend(first);
                } else {
                    // A duplication: "C <DUP>" becomes "C CAT".
                    if !fai.fetch(name, pos, pos + sv_len, &mut seqs.reference)? {
                        bail!("faidx_fetch_seq failed at {}:{}", ctg(), pos + 1);
                    }
                    nref = seqs.reference.len();
                    nt_to_upper(&mut seqs.reference);
                    replace_iupac_codes(&mut seqs.reference);
                    seqs.als[i].extend_from_slice(&seqs.reference);
                }
                seqs.sym[i].extend_from_slice(line.allele(i));
                continue;
            }
            if i > 0 {
                symbolic_alts = false;
            }
            if allele.first() == Some(&b'*') {
                return Ok(Realigned::SpanningDeletion);
            }
            if has_non_acgtn(allele) {
                let v = String::from_utf8_lossy(allele);
                if check_ref == CHECK_REF_EXIT {
                    bail!(
                        "Non-ACGTN alternate allele at {}:{} .. VCF:'{v}'",
                        seqname(),
                        pos + 1
                    );
                }
                if check_ref & CHECK_REF_WARN != 0 {
                    eprintln!("NON_ACGTN_ALT\t{}\t{}\t{v}", seqname(), pos + 1);
                }
                return Ok(Realigned::RefMismatch);
            }
            seqs.als[i].clear();
            seqs.als[i].extend_from_slice(allele);
            nt_to_upper(&mut seqs.als[i]);
            if i > 0
                && seqs.als[i].len() == seqs.als[0].len()
                && eq_nocase(&seqs.als[0], &seqs.als[i])
            {
                return Ok(Realigned::DupAllele);
            }
        }

        let left = self.is_left_align(hdr, line);
        let new_pos = if left {
            self.realign_left(hdr, line, n)?
        } else {
            self.realign_right(hdr, line, n)?
        };

        // Have the alleles changed? A <DEL> may have expanded REF: then it is trimmed back,
        // but the REF length keeps the whole deletion.
        let seqs = &mut self.seqs;
        let mut new_reflen = seqs.als[0].len();
        if symbolic_alts {
            seqs.als[0].truncate(1);
            new_reflen = line.rlen() as usize;
        }
        if new_pos == pos && eq_nocase(line.allele(0), &seqs.als[0]) {
            return Ok(Realigned::Ok);
        }

        self.old_rec_tag_set(out_hdr, line, 0)?;

        let seqs = &mut self.seqs;
        seqs.joined.clear();
        for i in 0..n {
            if i > 0 {
                seqs.joined.push(b',');
            }
            if seqs.sym[i].is_empty() {
                seqs.joined.extend_from_slice(&seqs.als[i]);
            } else {
                seqs.joined.extend_from_slice(&seqs.sym[i]);
            }
        }
        line.update_alleles_str(out_hdr, &seqs.joined);
        self.stats.changed += 1;

        // Update INFO/END if necessary.
        let end = &mut self.scratch.i.a;
        if (new_pos != pos || reflen != new_reflen) && line.get_info_i32(hdr, b"END", end) == Ok(1)
        {
            line.set_pos(new_pos);
            let v = (new_pos + new_reflen as i64) as i32;
            let _ = line.update_info_i32(out_hdr, b"END", &[v]);
        }
        line.set_pos(new_pos);
        Ok(Realigned::Ok)
    }

    /// Whether to left-align: always, except with `--right-align`, or within transcripts all
    /// on the forward strand with `-g` (bcftools' `is_left_align`).
    fn is_left_align(&self, hdr: &Header, line: &Record) -> bool {
        if self.opts.right_align {
            return false;
        }
        let Some(gff) = &self.gff else {
            return true;
        };
        let chr = line.chrom_safe(hdr);
        gff.left_align(chr, line.pos(), line.pos() + line.rlen())
    }

    /// Trims identical bases on the right, padding from the reference on the left, then on the
    /// left (bcftools' `realign_left`). Returns the new position.
    fn realign_left(&mut self, hdr: &Header, line: &Record, n: usize) -> anyhow::Result<i64> {
        let mut new_pos = line.pos();
        let aln_win = self.opts.aln_win;
        loop {
            // Is the rightmost base identical in all alleles?
            let als = &self.seqs.als;
            let last0 = als[0].last().map(u8::to_ascii_uppercase);
            let mut min_len = als[0].len();
            let mut same = true;
            for a in als.iter().take(n).skip(1) {
                if a.last().map(u8::to_ascii_uppercase) != last0 {
                    same = false;
                    break;
                }
                min_len = min_len.min(a.len());
            }
            if !same || last0.is_none() {
                break;
            }
            if min_len <= 1 && new_pos == 0 {
                break;
            }
            let mut pad_from_left = false;
            for a in self.seqs.als.iter_mut().take(n) {
                a.pop();
                if a.is_empty() {
                    pad_from_left = true;
                }
            }
            if pad_from_left {
                // Extend all alleles to the left by aln_win bases (unless close to the start);
                // the extra bases are trimmed below.
                let npad = new_pos.min(aln_win);
                let fai = self.fai.as_mut().expect("realignment needs a reference");
                let name = hdr.contig_name(line.rid()).unwrap_or(b"(unknown)");
                if !fai.fetch(name, new_pos - npad, new_pos - 1, &mut self.seqs.pad)? {
                    bail!(
                        "faidx_fetch_seq64 failed at {}:{}",
                        String::from_utf8_lossy(line.chrom_safe(hdr)),
                        new_pos - npad + 1
                    );
                }
                replace_iupac_codes(&mut self.seqs.pad);
                let pad = &self.seqs.pad;
                for a in self.seqs.als.iter_mut().take(n) {
                    a.splice(0..0, pad.iter().copied());
                }
                new_pos -= npad;
            }
        }

        // Trim from the left.
        let als = &mut self.seqs.als;
        let mut ntrim = 0;
        loop {
            let mut min_len = als[0].len() - ntrim;
            let mut same = true;
            for a in als.iter().take(n).skip(1) {
                if !als[0][ntrim].eq_ignore_ascii_case(&a[ntrim]) {
                    same = false;
                    break;
                }
                min_len = min_len.min(a.len() - ntrim);
            }
            if !same || min_len <= 1 {
                break;
            }
            ntrim += 1;
        }
        if ntrim > 0 {
            for a in als.iter_mut().take(n) {
                a.drain(..ntrim);
            }
            new_pos += ntrim as i64;
        }
        Ok(new_pos)
    }

    /// Trims identical bases on the left, padding from the reference on the right, then on the
    /// right (bcftools' `realign_right`). Returns the new position.
    fn realign_right(&mut self, hdr: &Header, line: &Record, n: usize) -> anyhow::Result<i64> {
        let mut new_pos = line.pos();
        let aln_win = self.opts.aln_win;
        let mut ntrim = 0usize;
        let mut npad_right = line.rlen();
        let mut has_indel = 0usize;
        loop {
            // Is the leftmost untrimmed base identical in all alleles?
            let als = &self.seqs.als;
            let Some(&b0) = als[0].get(ntrim) else {
                break;
            };
            let mut min_len = als[0].len() - ntrim;
            let mut same = true;
            for a in als.iter().take(n).skip(1) {
                if als[0].len() != a.len() {
                    has_indel = 1;
                }
                match a.get(ntrim) {
                    Some(b) if b.eq_ignore_ascii_case(&b0) => {}
                    _ => {
                        same = false;
                        break;
                    }
                }
                min_len = min_len.min(a.len() - ntrim);
            }
            if !same {
                break;
            }
            ntrim += 1;
            if min_len <= 1 {
                // Pad from the right. bcftools fetches aln_win+1 bases here and so repeats one
                // base at each later padding; this fetches the aln_win bases that follow.
                let fai = self.fai.as_mut().expect("realignment needs a reference");
                let name = hdr.contig_name(line.rid()).unwrap_or(b"(unknown)");
                let beg = line.pos() + npad_right;
                if !fai.fetch(name, beg, beg + aln_win - 1, &mut self.seqs.pad)? {
                    bail!(
                        "faidx_fetch_seq64 failed at {}:{}",
                        String::from_utf8_lossy(line.chrom_safe(hdr)),
                        new_pos + ntrim as i64
                    );
                }
                if self.seqs.pad.is_empty() {
                    // The end of the sequence: nothing to pad with. Alleles of one length
                    // that are the same to the end (a one-base <DEL>) keep the last base.
                    if has_indel == 0 {
                        ntrim -= 1;
                    }
                    break;
                }
                npad_right += aln_win;
                replace_iupac_codes(&mut self.seqs.pad);
                let pad = &self.seqs.pad;
                for a in self.seqs.als.iter_mut().take(n) {
                    a.extend_from_slice(pad);
                }
            }
        }
        let ntrim = ntrim.saturating_sub(has_indel);
        let als = &mut self.seqs.als;
        if ntrim > 0 {
            for a in als.iter_mut().take(n) {
                a.drain(..ntrim.min(a.len()));
            }
            new_pos += ntrim as i64;
        }

        // Trim from the right.
        loop {
            let last0 = als[0].last().map(u8::to_ascii_uppercase);
            let mut min_len = als[0].len();
            let mut same = last0.is_some();
            for a in als.iter().take(n).skip(1) {
                if a.last().map(u8::to_ascii_uppercase) != last0 {
                    same = false;
                    break;
                }
                min_len = min_len.min(a.len());
            }
            if !same || min_len <= 1 {
                break;
            }
            for a in als.iter_mut().take(n) {
                a.pop();
            }
        }
        Ok(new_pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_codes_become_their_first_base() {
        let mut s = b"ACGTNRYSWKMBDHVrykx*".to_vec();
        assert_eq!(replace_iupac_codes(&mut s), 13);
        assert_eq!(s, b"ACGTNACCAGACAAAacgx*");
    }

    #[test]
    fn non_acgtn_bases_are_found() {
        assert!(!has_non_acgtn(b"ACGTNacgtn"));
        assert!(has_non_acgtn(b"ACR"));
        assert!(has_non_acgtn(b"."));
        assert!(!has_non_acgtn(b"AC\0R"));
    }
}
