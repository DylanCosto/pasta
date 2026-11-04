//! Decomposing complex variants into their parts (`-a`; bcftools' `abuf.c` in SPLIT mode).
//!
//! Each ALT allele is split into atoms by trimming and greedy matching from the left, without
//! alignment. Each distinct atom becomes a record with one ALT (plus `*` where other alleles
//! overlap it), and a table maps the original alleles to the alleles of each new record: 1
//! when the original allele has the atom, 0 when it does not touch it, and 2 when it
//! overlaps it with something else.

use anyhow::bail;
use pasta::variant::{
    Genotype, Header, HeaderRecordType, Number, Record, TypedValues, Unpack, ValueType, encoding,
    number,
};

use super::{Normalizer, alleles2gt, strcasecmp, values::copy_string_field};

/// One part of an ALT allele (bcftools' `atom_t`).
#[derive(Clone, Debug, Default)]
struct Atom {
    ref_: Vec<u8>,
    alt: Vec<u8>,
    /// The index of the original ALT allele, 1-based.
    ial: usize,
    /// 0-based inclusive offsets into the original REF.
    beg: usize,
    end: usize,
    /// The length of the REF/ALT prefix shared, such as 1 for C>CA.
    plen: usize,
}

/// The atomization state (bcftools' `abuf_t` and `split_t`).
#[derive(Default)]
pub(super) struct Atomizer {
    atoms: Vec<Atom>,
    natoms: usize,
    /// The output atoms, as indexes into `atoms`.
    out_atoms: Vec<usize>,
    /// For each output row, how each original ALT maps to it (`nori` columns).
    tbl: Vec<u8>,
    nori: usize,
    /// Whether each output row needs the star allele.
    overlaps: Vec<bool>,
    /// The records made from the last pushed one.
    pub(super) out: Vec<Record>,
    history: Vec<u8>,
}

/// bcftools' `_atoms_inconsistent`: orders atoms by position, REF, and ALT.
fn atoms_inconsistent(a: &Atom, b: &Atom) -> std::cmp::Ordering {
    a.beg
        .cmp(&b.beg)
        .then_with(|| strcasecmp(&a.ref_, &b.ref_))
        .then_with(|| strcasecmp(&a.alt, &b.alt))
}

/// bcftools' `_atoms_overlap`: 0 for identical atoms, 1 for atoms that do not conflict
/// (an SNV next to an indel), 2 for conflicting ones.
fn atoms_overlap(a: &Atom, b: &Atom) -> u8 {
    if a.beg != b.beg {
        return 2;
    }
    if a.plen > 0 && a.plen >= b.ref_.len() {
        return 1;
    }
    if b.plen > 0 && b.plen >= a.ref_.len() {
        return 1;
    }
    if !a.ref_.eq_ignore_ascii_case(&b.ref_) {
        return 2;
    }
    if a.plen > 0 && a.plen >= b.alt.len() {
        return 1;
    }
    if b.plen > 0 && b.plen >= a.alt.len() {
        return 1;
    }
    if !a.alt.eq_ignore_ascii_case(&b.alt) {
        return 2;
    }
    0
}

/// How the values of original alleles that become one output allele combine (bcftools'
/// `merge_rule_t`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MergeRule {
    /// The value of the first allele.
    First,
    /// The sum of the values (`QS` and `AD`).
    Sum,
}

impl MergeRule {
    /// bcftools' hard-coded rule for a tag.
    fn of(tag: &[u8]) -> Self {
        if tag == b"QS" || tag == b"AD" {
            Self::Sum
        } else {
            Self::First
        }
    }
}

/// What [`Normalizer::atomize`] made of a record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Atomized {
    /// The record is output as it is: bcftools copies it (symbolic alleles, no ALT), or
    /// rebuilds it from a single atom that is the whole record and gets the same record.
    Unchanged,
    /// The records to output are in `out`.
    Split,
}

/// Whether rebuilding a biallelic record from its single, whole-record atom gives the record
/// back: every INFO and FORMAT field survives `atomize_info`/`atomize_format` with the same
/// values. Fields that would be dropped (Number=G INFO, wrong counts, mistyped values) or
/// reformatted (Number=A/R strings, Number=G FORMAT) send the record through the full path.
fn rebuild_is_identity(hdr: &Header, rec: &Record) -> bool {
    let infos = rec.info_fields();
    for (i, f) in infos.iter().enumerate() {
        if f.is_removed() || infos[..i].iter().any(|g| g.key() == f.key()) {
            return false;
        }
        let id = f.key() as u32;
        let number = hdr.id_number(HeaderRecordType::Info, id);
        let values = TypedValues::new(f.type_code(), rec.info_data(f));
        let ok = match hdr.id_type(HeaderRecordType::Info, id) {
            None => false,
            Some(_) if number == Some(Number::G) || is_unsplit_number(number) => false,
            Some(ValueType::Flag) => f.type_code() == 0,
            Some(ValueType::Integer | ValueType::Float) => {
                let is_int = hdr.id_type(HeaderRecordType::Info, id) == Some(ValueType::Integer);
                let nval = if is_int {
                    matches!(f.type_code(), 1..=3).then(|| {
                        (0..f.len())
                            .take_while(|&j| values.int(j) != Some(encoding::INT32_VECTOR_END))
                            .count()
                    })
                } else {
                    (f.type_code() == 5).then(|| {
                        (0..f.len())
                            .take_while(|&j| {
                                values.float(j).map(f32::to_bits)
                                    != Some(encoding::FLOAT_VECTOR_END)
                            })
                            .count()
                    })
                };
                match (nval, number) {
                    (None | Some(0), _) => false,
                    (Some(k), Some(Number::A)) => k == 1,
                    (Some(k), Some(Number::R)) => k == 2,
                    _ => true,
                }
            }
            Some(ValueType::String) => {
                let commas = values.bytes()[..f.len().min(values.bytes().len())]
                    .iter()
                    .take_while(|&&c| c != 0)
                    .filter(|&&c| c == b',')
                    .count();
                f.type_code() == 7
                    && !f.is_empty()
                    && match number {
                        Some(Number::A) => commas == 0,
                        Some(Number::R) => commas == 1,
                        _ => true,
                    }
            }
        };
        if !ok {
            return false;
        }
    }
    let fmts = rec.format_fields();
    if fmts.is_empty() {
        return true;
    }
    if hdr.sample_count() == 0 || rec.sample_count() != hdr.sample_count() {
        return false;
    }
    for (i, f) in fmts.iter().enumerate() {
        if f.is_removed() || fmts[..i].iter().any(|g| g.id() == f.id()) {
            return false;
        }
        let id = f.id() as u32;
        let ty = hdr.id_type(HeaderRecordType::Format, id);
        let number = hdr.id_number(HeaderRecordType::Format, id);
        let ok = if is_unsplit_number(number) {
            false
        } else if hdr.id_name(id) == Some(b"GT") {
            // Alleles 0 and 1 (raw 2..=5), missing alleles (0, 1), and vector ends map to
            // themselves; anything else is out of bounds or missing and handled in full. A
            // rebuilt record has GT first, as HTSlib's `bcf_update_format` puts it there.
            i == 0
                && ty == Some(ValueType::String)
                && f.type_code() == 1
                && rec
                    .format_data(f)
                    .iter()
                    .all(|&b| b <= 5 || b as i8 == encoding::INT8_VECTOR_END)
        } else {
            match ty {
                Some(ValueType::Integer) => {
                    matches!(f.type_code(), 1..=3) && fmt_number_is_identity(number, f.n())
                }
                Some(ValueType::Float) => {
                    f.type_code() == 5 && fmt_number_is_identity(number, f.n())
                }
                // A string field no sample has (zero bytes each) is dropped: HTSlib's
                // `bcf_update_format` removes a field given no values.
                Some(ValueType::String) => {
                    f.type_code() == 7
                        && f.n() > 0
                        && !matches!(number, Some(Number::A | Number::R | Number::G))
                }
                _ => false,
            }
        };
        if !ok {
            return false;
        }
    }
    true
}

/// Whether fields of this Number (VCF 4.4 and 4.5's P, LA, LG, LR, M) are left out of the
/// records `abuf.c` builds: its splitters handle only fixed, `.`, A, R and G counts.
fn is_unsplit_number(number: Option<Number>) -> bool {
    matches!(
        number,
        Some(Number::P | Number::LocalA | Number::LocalG | Number::LocalR | Number::M)
    )
}

/// Whether a numeric FORMAT field of a biallelic record with `n` values per sample is
/// rebuilt as it is.
fn fmt_number_is_identity(number: Option<Number>, n: usize) -> bool {
    match number {
        Some(Number::A) => n == 1,
        Some(Number::R) => n == 2,
        Some(Number::G) => false,
        _ => true,
    }
}

/// bcftools' `_is_acgtn`.
fn is_acgtn(seq: &[u8]) -> bool {
    seq.iter()
        .all(|c| matches!(c.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
}

impl Atomizer {
    fn new_atom(&mut self, refb: u8, altb: u8, i: usize, ial: usize, plen: usize) -> usize {
        if self.natoms == self.atoms.len() {
            self.atoms.push(Atom::default());
        }
        let a = &mut self.atoms[self.natoms];
        a.ref_.clear();
        a.ref_.push(refb);
        a.alt.clear();
        a.alt.push(altb);
        a.beg = i;
        a.end = i;
        a.ial = ial;
        a.plen = plen;
        self.natoms += 1;
        self.natoms - 1
    }

    /// Splits ALT allele `ial` into atoms (bcftools' `_atomize_allele`).
    fn atomize_allele(&mut self, rec: &Record, ial: usize) {
        let r = rec.allele(0);
        let a = rec.allele(ial);
        let (mut rlen, mut alen) = (r.len(), a.len());
        while rlen > 1 && alen > 1 && r[rlen - 1] == a[alen - 1] {
            rlen -= 1;
            alen -= 1;
        }
        let mlen = rlen.max(alen);
        let mut atom: Option<usize> = None;
        for i in 0..mlen {
            let refb = if i < rlen { r[i] } else { b'-' };
            let altb = if i < alen { a[i] } else { b'-' };
            if refb != altb {
                if refb == b'-' || altb == b'-' {
                    let at = &mut self.atoms[atom.expect("an atom precedes an indel")];
                    if altb != b'-' {
                        at.alt.push(altb);
                    }
                    if refb != b'-' {
                        at.ref_.push(refb);
                        at.end += 1;
                    }
                    continue;
                }
                atom = Some(self.new_atom(refb, altb, i, ial, 0));
                if rlen != alen && (i + 1 >= rlen || i + 1 >= alen) {
                    // The next base is an indel combined with an SNV, as in C>GGG.
                    atom = Some(self.new_atom(refb, refb, i, ial, 1));
                }
                continue;
            }
            if i + 1 >= rlen || i + 1 >= alen {
                // The next base is an indel.
                atom = Some(self.new_atom(refb, altb, i, ial, 0));
            }
        }
    }

    #[inline]
    fn tbl(&self, iout: usize, col: usize) -> u8 {
        self.tbl[iout * self.nori + col]
    }

    /// The allele of output row `iout` that original allele `ial` becomes (bcftools'
    /// `_split_table_get_ial`).
    fn get_ial(&self, iout: usize, ial: usize) -> usize {
        if ial == 0 {
            0
        } else {
            usize::from(self.tbl(iout, ial - 1))
        }
    }
}

impl Normalizer {
    /// Atomizes a record into `out` (bcftools' `abuf_push` and `_abuf_split`), or returns
    /// [`Atomized::Unchanged`] when the output would be the record itself.
    ///
    /// The new records are built with the reader's header, as bcftools hands `abuf` the
    /// header before duplicating it for output: tags that reading added to the header (ones
    /// the file did not declare) can be set, and fail later, on writing.
    pub(super) fn atomize(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        out: &mut Vec<Record>,
    ) -> anyhow::Result<Atomized> {
        rec.unpack(Unpack::ALL);
        out.clear();
        let n = rec.allele_count();
        if n < 2 || rec.alleles().any(|a| !is_acgtn(a)) {
            return Ok(Atomized::Unchanged);
        }

        let ab = &mut self.atomizer;
        ab.natoms = 0;
        for i in 1..n {
            ab.atomize_allele(rec, i);
        }
        let natoms = ab.natoms;
        if n == 2 && natoms == 1 {
            let atom = &ab.atoms[0];
            if atom.beg == 0
                && atom.ref_.len() == rec.allele(0).len()
                && atom.alt.len() == rec.allele(1).len()
                && self.opts.old_rec_tag.is_none()
                && rebuild_is_identity(hdr, rec)
            {
                return Ok(Atomized::Unchanged);
            }
        }
        ab.atoms[..natoms].sort_by(|a, b| atoms_inconsistent(a, b).then(a.ial.cmp(&b.ial)));

        // The table: one row per distinct atom.
        ab.nori = n - 1;
        ab.out_atoms.clear();
        ab.tbl.clear();
        for i in 0..natoms {
            if i > 0
                && atoms_inconsistent(&ab.atoms[i - 1], &ab.atoms[i]) == std::cmp::Ordering::Equal
            {
                continue;
            }
            ab.out_atoms.push(i);
            let row = ab.tbl.len();
            ab.tbl.resize(row + ab.nori, 0);
            ab.tbl[row + ab.atoms[i].ial - 1] = 1;
        }
        let nout = ab.out_atoms.len();
        ab.overlaps.clear();
        ab.overlaps.resize(nout, false);
        for i in 0..natoms {
            for j in 0..nout {
                let o = ab.out_atoms[j];
                if o == i {
                    continue;
                }
                let (atom, outa) = (&ab.atoms[i], &ab.atoms[o]);
                if atom.beg > outa.end {
                    continue;
                }
                if atom.end < outa.beg {
                    break;
                }
                let olap = atoms_overlap(atom, outa);
                let col = atom.ial - 1;
                ab.tbl[j * ab.nori + col] = if olap > 1 { 2 } else { 1 };
                ab.overlaps[j] = true;
            }
        }

        // CHROM-QUAL.
        let use_star = self.opts.use_star_allele;
        for iout in 0..nout {
            let mut o = self.take_record();
            let atom = &self.atomizer.atoms[self.atomizer.out_atoms[iout]];
            o.set_rid(rec.rid());
            o.set_pos(rec.pos() + atom.beg as i64);
            o.update_id(Some(rec.id()));
            let star = use_star && self.atomizer.overlaps[iout];
            if star {
                o.update_alleles(hdr, &[&atom.ref_[..], &atom.alt[..], &b"*"[..]]);
            } else {
                o.update_alleles(hdr, &[&atom.ref_[..], &atom.alt[..]]);
            }
            o.set_qual(rec.qual());
            o.update_filter(rec.filters());
            out.push(o);
        }

        // INFO.
        let mut keys = std::mem::take(&mut self.scratch.info_keys);
        keys.clear();
        keys.extend(rec.info_fields().iter().map(|f| f.key()));
        for &key in &keys {
            let tag = hdr.id_name(key as u32).unwrap_or_default();
            self.atomize_info(hdr, rec, key, tag, MergeRule::of(tag), out)?;
        }
        self.scratch.info_keys = keys;
        // The INFO tag showing the original record.
        if let Some(tag) = self.opts.old_rec_tag.take() {
            let r = self.atomize_history(hdr, rec, &tag, out);
            self.opts.old_rec_tag = Some(tag);
            r?;
        }
        // FORMAT.
        let mut ids = std::mem::take(&mut self.scratch.fmt_ids);
        ids.clear();
        ids.extend(rec.format_fields().iter().map(|f| f.id()));
        for &id in &ids {
            let tag = hdr.id_name(id as u32).unwrap_or_default();
            self.atomize_format(hdr, rec, id, tag, MergeRule::of(tag), out)?;
        }
        self.scratch.fmt_ids = ids;
        for o in out.iter_mut() {
            if o.sample_count() == 0 {
                o.set_sample_count(rec.sample_count());
            }
            o.sync();
        }
        Ok(Atomized::Split)
    }

    fn star(&self, iout: usize) -> bool {
        self.opts.use_star_allele && self.atomizer.overlaps[iout]
    }

    /// bcftools' `_split_table_set_info`.
    fn atomize_info(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        key: i32,
        tag: &[u8],
        rule: MergeRule,
        out: &mut [Record],
    ) -> anyhow::Result<()> {
        let ty = hdr.id_type(HeaderRecordType::Info, key as u32);
        let number = hdr.id_number(HeaderRecordType::Info, key as u32);
        if number == Some(Number::G) || is_unsplit_number(number) {
            return Ok(());
        }
        let n = rec.allele_count();
        let ab = &self.atomizer;
        let fail = |e: pasta::variant::UpdateError| {
            anyhow::anyhow!(
                "An error occurred while updating INFO/{} (errcode={e})",
                String::from_utf8_lossy(tag)
            )
        };
        match ty {
            Some(ValueType::Integer | ValueType::Float) => {
                let is_int = ty == Some(ValueType::Integer);
                let mut iv = std::mem::take(&mut self.scratch.i.a);
                let mut fv = std::mem::take(&mut self.scratch.f.a);
                let mut w = std::mem::take(&mut self.scratch.i.b);
                let mut fw = std::mem::take(&mut self.scratch.f.b);
                let nval = if is_int {
                    rec.get_info_i32(hdr, tag, &mut iv)
                } else {
                    rec.get_info_f32(hdr, tag, &mut fv)
                };
                let r = (|| -> anyhow::Result<()> {
                    let Ok(nval) = nval else {
                        return Ok(());
                    };
                    if (number == Some(Number::A) && nval != n - 1)
                        || (number == Some(Number::R) && nval != n)
                    {
                        return Ok(());
                    }
                    for (iout, o) in out.iter_mut().enumerate() {
                        let star = self.opts.use_star_allele && ab.overlaps[iout];
                        let atom = &ab.atoms[ab.out_atoms[iout]];
                        // Build the values as raw 32-bit words, as bcftools does.
                        w.clear();
                        let word = |k: usize| -> i32 {
                            if is_int {
                                iv[k]
                            } else {
                                fv[k].to_bits() as i32
                            }
                        };
                        let missing = if is_int {
                            encoding::INT32_MISSING
                        } else {
                            encoding::FLOAT_MISSING as i32
                        };
                        let vend = if is_int {
                            encoding::INT32_VECTOR_END
                        } else {
                            encoding::FLOAT_VECTOR_END as i32
                        };
                        match number {
                            Some(Number::A) => {
                                let v = word(atom.ial - 1);
                                w.push(if v == vend { missing } else { v });
                                if star {
                                    w.push(missing);
                                }
                            }
                            Some(Number::R) => {
                                w.push(word(0));
                                let iori = atom.ial;
                                let v = word(iori);
                                let mut v1 = if v == vend { missing } else { v };
                                if is_int && rule == MergeRule::Sum {
                                    for i in iori..ab.nori {
                                        if ab.tbl(iout, i) == 1 {
                                            v1 = v1.wrapping_add(iv[i + 1]);
                                        }
                                    }
                                }
                                w.push(v1);
                                if star {
                                    w.push(missing);
                                }
                            }
                            _ => w.extend((0..nval).map(word)),
                        }
                        let r = if is_int {
                            o.update_info_i32(hdr, tag, &w)
                        } else {
                            fw.clear();
                            fw.extend(w.iter().map(|&x| f32::from_bits(x as u32)));
                            o.update_info_f32(hdr, tag, &fw)
                        };
                        r.map_err(fail)?;
                    }
                    Ok(())
                })();
                self.scratch.i.a = iv;
                self.scratch.f.a = fv;
                self.scratch.i.b = w;
                self.scratch.f.b = fw;
                r
            }
            Some(ValueType::Flag) => {
                let set = rec.get_info_flag(hdr, tag).unwrap_or(false);
                for o in out.iter_mut() {
                    o.update_info_flag(hdr, tag, set).map_err(fail)?;
                }
                Ok(())
            }
            _ => {
                let Ok(s) = rec.get_info_string(hdr, tag) else {
                    return Ok(());
                };
                let dst = &mut self.scratch.bytes2;
                for (iout, o) in out.iter_mut().enumerate() {
                    let star = self.opts.use_star_allele && ab.overlaps[iout];
                    let atom = &ab.atoms[ab.out_atoms[iout]];
                    match number {
                        Some(Number::A) => {
                            dst.clear();
                            dst.push(b'.');
                            if star {
                                dst.extend_from_slice(b",.");
                            }
                            copy_string_field(s, atom.ial - 1, dst, 0);
                            o.update_info_string(hdr, tag, Some(dst)).map_err(fail)?;
                        }
                        Some(Number::R) => {
                            dst.clear();
                            dst.extend_from_slice(b".,.");
                            if star {
                                dst.extend_from_slice(b",.");
                            }
                            copy_string_field(s, 0, dst, 0);
                            copy_string_field(s, atom.ial, dst, 1);
                            o.update_info_string(hdr, tag, Some(dst)).map_err(fail)?;
                        }
                        _ => {
                            // An empty value gives no values, which removes the tag.
                            let s = Some(s).filter(|s| !s.is_empty());
                            o.update_info_string(hdr, tag, s).map_err(fail)?;
                        }
                    }
                }
                Ok(())
            }
        }
    }

    /// Sets the `--old-rec-tag` INFO on the atoms (bcftools' `_split_table_set_history`).
    fn atomize_history(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        tag: &[u8],
        out: &mut [Record],
    ) -> anyhow::Result<()> {
        // Don't overwrite the tag (as set by -m).
        if rec.get_info_string(hdr, tag).is_ok_and(|s| !s.is_empty()) {
            return Ok(());
        }
        let h = &mut self.atomizer.history;
        h.clear();
        h.extend_from_slice(rec.chrom_safe(hdr));
        h.push(b'|');
        number::push_int(h, rec.pos() + 1);
        h.push(b'|');
        h.extend_from_slice(rec.allele(0));
        h.push(b'|');
        for i in 1..rec.allele_count() {
            h.extend_from_slice(rec.allele(i));
            h.push(b',');
        }
        let len = h.len();
        h[len - 1] = b'|';
        let ab = &mut self.atomizer;
        for (iout, o) in out.iter_mut().enumerate() {
            ab.history.truncate(len);
            for j in 0..ab.nori {
                if ab.tbl[iout * ab.nori + j] != 1 {
                    continue;
                }
                number::push_uint(&mut ab.history, j as u64 + 1);
                ab.history.push(b',');
            }
            ab.history.pop();
            if o.update_info_string(hdr, tag, Some(&ab.history)).is_err() {
                bail!(
                    "An error occurred while updating INFO/{}",
                    String::from_utf8_lossy(tag)
                );
            }
        }
        Ok(())
    }

    /// bcftools' `_split_table_set_gt`.
    fn atomize_gt(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        out: &mut [Record],
    ) -> anyhow::Result<()> {
        let nsmpl = hdr.sample_count();
        if nsmpl == 0 {
            return Ok(());
        }
        let mut gt = std::mem::take(&mut self.scratch.i.a);
        let mut dst = std::mem::take(&mut self.scratch.i.out);
        let r = (|| -> anyhow::Result<()> {
            let ngt = match rec.get_genotypes(hdr, &mut gt) {
                Ok(n) if n > 0 => n,
                _ => return Ok(()),
            };
            let n = rec.allele_count();
            for (iout, o) in out.iter_mut().enumerate() {
                let star = self.star(iout);
                dst.clear();
                for &src in &gt[..ngt] {
                    let g = Genotype(src);
                    if g.is_vector_end() || g.is_missing() {
                        dst.push(src);
                        continue;
                    }
                    let iori = g.allele();
                    if iori < 0 || iori as usize >= n {
                        bail!(
                            "Out-of-bounds genotypes at {}:{}",
                            String::from_utf8_lossy(rec.chrom_safe(hdr)),
                            rec.pos() + 1
                        );
                    }
                    let ial = self.atomizer.get_ial(iout, iori as usize);
                    dst.push(if ial == 2 && !star {
                        Genotype::MISSING.0 | (src & 1)
                    } else if g.is_phased() {
                        Genotype::phased(ial as i32).0
                    } else {
                        Genotype::unphased(ial as i32).0
                    });
                }
                let _ = o.update_genotypes(hdr, &dst);
            }
            Ok(())
        })();
        self.scratch.i.a = gt;
        self.scratch.i.out = dst;
        r
    }

    /// bcftools' `_split_table_set_format`.
    fn atomize_format(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        id: i32,
        tag: &[u8],
        rule: MergeRule,
        out: &mut [Record],
    ) -> anyhow::Result<()> {
        let nsmpl = hdr.sample_count();
        if nsmpl == 0 {
            return Ok(());
        }
        if tag == b"GT" {
            return self.atomize_gt(hdr, rec, out);
        }
        let ty = hdr.id_type(HeaderRecordType::Format, id as u32);
        let number = hdr.id_number(HeaderRecordType::Format, id as u32);
        if is_unsplit_number(number) {
            return Ok(());
        }
        let n = rec.allele_count();
        let fail = |e: pasta::variant::UpdateError| {
            anyhow::anyhow!(
                "An error occurred while updating FORMAT/{} (errcode={e})",
                String::from_utf8_lossy(tag)
            )
        };
        let ab = &self.atomizer;
        if !matches!(ty, Some(ValueType::Integer | ValueType::Float)) {
            if number == Some(Number::G) {
                return Ok(());
            }
            let Ok((data, nval1)) = rec.get_format_chars(hdr, tag) else {
                return Ok(());
            };
            let buf = &mut self.scratch.bytes2;
            let s = &mut self.scratch.bytes;
            for (iout, o) in out.iter_mut().enumerate() {
                let star = self.opts.use_star_allele && ab.overlaps[iout];
                match number {
                    Some(Number::A | Number::R) => {
                        let is_r = number == Some(Number::R);
                        let ioff = usize::from(is_r);
                        let iori = ab.atoms[ab.out_atoms[iout]].ial - 1;
                        let w = if star { nval1 + 2 } else { nval1 };
                        buf.clear();
                        buf.resize(w * nsmpl, 0);
                        for i in 0..nsmpl {
                            let src = &data[i * nval1..(i + 1) * nval1];
                            s.clear();
                            s.push(b'.');
                            if star {
                                s.extend_from_slice(b",.");
                            }
                            if is_r {
                                s.extend_from_slice(b",.");
                                copy_string_field(src, 0, s, 0);
                            }
                            copy_string_field(src, iori + ioff, s, ioff);
                            let k = s.len().min(w);
                            buf[i * w..i * w + k].copy_from_slice(&s[..k]);
                        }
                        o.update_format_chars(hdr, tag, buf).map_err(fail)?;
                    }
                    _ => o.update_format_chars(hdr, tag, data).map_err(fail)?,
                }
            }
            return Ok(());
        }

        let is_int = ty == Some(ValueType::Integer);
        let mut iv = std::mem::take(&mut self.scratch.i.a);
        let mut fv = std::mem::take(&mut self.scratch.f.a);
        let mut w = std::mem::take(&mut self.scratch.i.b);
        let mut fw = std::mem::take(&mut self.scratch.f.b);
        let nval = if is_int {
            rec.get_format_i32(hdr, tag, &mut iv)
        } else {
            rec.get_format_f32(hdr, tag, &mut fv)
        };
        let r = (|| -> anyhow::Result<()> {
            let Ok(nval) = nval else {
                return Ok(());
            };
            if number == Some(Number::G) && nval != nsmpl * n && nval != nsmpl * n * (n + 1) / 2 {
                return Ok(());
            }
            if (number == Some(Number::A) && nval != nsmpl * (n - 1))
                || (number == Some(Number::R) && nval != nsmpl * n)
            {
                return Ok(());
            }
            let nval1 = nval / nsmpl;
            let word = |k: usize| -> i32 {
                if is_int {
                    iv[k]
                } else {
                    fv[k].to_bits() as i32
                }
            };
            let (missing, vend) = if is_int {
                (encoding::INT32_MISSING, encoding::INT32_VECTOR_END)
            } else {
                (
                    encoding::FLOAT_MISSING as i32,
                    encoding::FLOAT_VECTOR_END as i32,
                )
            };
            for (iout, o) in out.iter_mut().enumerate() {
                let star = self.opts.use_star_allele && ab.overlaps[iout];
                let atom = &ab.atoms[ab.out_atoms[iout]];
                w.clear();
                match number {
                    Some(Number::A) => {
                        let iori = atom.ial - 1;
                        for i in 0..nsmpl {
                            let v = word(i * nval1 + iori);
                            w.push(if v == vend { missing } else { v });
                            if star {
                                w.push(missing);
                            }
                        }
                    }
                    Some(Number::R) => {
                        let iori = atom.ial;
                        for i in 0..nsmpl {
                            let base = i * nval1;
                            w.push(word(base));
                            let mut v1 = word(base + iori);
                            if is_int && rule == MergeRule::Sum {
                                for j in iori..ab.nori {
                                    if ab.tbl(iout, j) == 1 {
                                        v1 = v1.wrapping_add(iv[base + j + 1]);
                                    }
                                }
                            }
                            w.push(v1);
                            if star {
                                w.push(missing);
                            }
                        }
                    }
                    Some(Number::G) => {
                        let iori = atom.ial;
                        let (i01, i11) = (alleles2gt(0, iori), alleles2gt(iori, iori));
                        for i in 0..nsmpl {
                            let base = i * nval1;
                            let src = |k: usize| word(base + k);
                            let ploidy = (0..nval1).take_while(|&k| src(k) != vend).count();
                            w.push(src(0));
                            let (mut nmiss, mut nend) = (0, 0);
                            if ploidy == n {
                                // Haploid.
                                w.push(src(iori));
                                if star {
                                    nmiss = 1;
                                    nend = 3;
                                } else {
                                    nend = 1;
                                }
                            } else if ploidy == nval1 {
                                w.push(src(i01));
                                w.push(src(i11));
                                if star {
                                    nmiss = 3;
                                }
                            } else if ploidy == 1 && src(0) == missing {
                                nend = if star { 5 } else { 2 };
                            } else {
                                bail!(
                                    "Incorrect number of values at {}:{} .. tag=FORMAT/{} Number=G nAlleles={n} nValues={ploidy}, {}-th sample",
                                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                                    rec.pos() + 1,
                                    String::from_utf8_lossy(tag),
                                    i + 1
                                );
                            }
                            w.extend(std::iter::repeat_n(missing, nmiss));
                            w.extend(std::iter::repeat_n(vend, nend));
                        }
                    }
                    _ => w.extend((0..nval).map(word)),
                }
                let r = if is_int {
                    o.update_format_i32(hdr, tag, &w)
                } else {
                    fw.clear();
                    fw.extend(w.iter().map(|&x| f32::from_bits(x as u32)));
                    o.update_format_f32(hdr, tag, &fw)
                };
                r.map_err(fail)?;
            }
            Ok(())
        })();
        self.scratch.i.a = iv;
        self.scratch.f.a = fv;
        self.scratch.i.b = w;
        self.scratch.f.b = fw;
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms(r: &[u8], a: &[u8]) -> Vec<(usize, String, String, usize)> {
        let mut rec = Record::new();
        let header = Header::new();
        rec.update_alleles(&header, &[r, a]);
        let mut ab = Atomizer::default();
        ab.atomize_allele(&rec, 1);
        ab.atoms[..ab.natoms]
            .iter()
            .map(|x| {
                (
                    x.beg,
                    String::from_utf8_lossy(&x.ref_).into_owned(),
                    String::from_utf8_lossy(&x.alt).into_owned(),
                    x.plen,
                )
            })
            .collect()
    }

    #[test]
    fn alleles_split_into_atoms() {
        assert_eq!(
            atoms(b"CC", b"TT"),
            [
                (0, "C".into(), "T".into(), 0),
                (1, "C".into(), "T".into(), 0)
            ]
        );
        assert_eq!(
            atoms(b"GCGT", b"GTGA"),
            [
                (1, "C".into(), "T".into(), 0),
                (3, "T".into(), "A".into(), 0)
            ]
        );
        assert_eq!(atoms(b"CA", b"C"), [(0, "CA".into(), "C".into(), 0)]);
        assert_eq!(
            atoms(b"C", b"GGG"),
            [
                (0, "C".into(), "G".into(), 0),
                (0, "C".into(), "CGG".into(), 1)
            ]
        );
    }

    #[test]
    fn overlapping_atoms_conflict_unless_snv_next_to_indel() {
        let atom = |r: &[u8], a: &[u8], plen| Atom {
            ref_: r.to_vec(),
            alt: a.to_vec(),
            plen,
            ..Atom::default()
        };
        assert_eq!(atoms_overlap(&atom(b"C", b"T", 0), &atom(b"C", b"T", 0)), 0);
        assert_eq!(atoms_overlap(&atom(b"C", b"T", 0), &atom(b"C", b"G", 0)), 2);
        assert_eq!(
            atoms_overlap(&atom(b"CC", b"C", 1), &atom(b"C", b"T", 0)),
            1
        );
        assert_eq!(
            atoms_overlap(&atom(b"T", b"TA", 0), &atom(b"T", b"TAAAA", 0)),
            2
        );
    }
}
