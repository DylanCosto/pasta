//! Joining biallelic records at one position into multiallelic ones (`-m+`; bcftools'
//! `mrows_*` and `merge_biallelics_to_multiallelic` with the `merge_*` functions).

use anyhow::bail;
use pasta::variant::{Genotype, Header, HeaderRecordType, Number, Record, Unpack, ValueType};

use super::{
    Collapse, Normalizer, alleles2gt,
    ploidy::{genotype_count, genotype_map, ploidy_of},
    values::{Value, copy_string_field},
};
use pasta::variant::VariantType;

/// The most bytes the joined values of a FORMAT field with a sample of ploidy above two may take
/// across the samples (1 GiB): `-m+` stops past it, as bcftools stops on any such field, rather
/// than build gigabytes. Half of what HTSlib's VCF reader allows all the FORMAT fields of a
/// record ("Excessive memory required by FORMAT fields"); tetraploid PL of 12 alleles in 100,000
/// samples takes about half of it.
const MAX_JOINED_FORMAT_BYTES: usize = 1 << 30;

/// [`MAX_JOINED_FORMAT_BYTES`] for strings (512 MiB). Their joined values are built in a row per
/// sample, a scratch row and the padded block of the record, so they take about four times
/// their bytes at the peak (3.2 GB for 772 MB of them).
const MAX_JOINED_STRING_BYTES: usize = MAX_JOINED_FORMAT_BYTES / 2;

/// The records waiting to be joined (bcftools' `mrows`).
#[derive(Default)]
pub(super) struct Joiner {
    pub(super) rows: Vec<Record>,
    first: usize,
    n: usize,
    /// For each joined row, its allele indexes in the joined record (bcftools' `maps`).
    maps: Vec<Vec<usize>>,
}

/// What [`Joiner::flush`] returns next.
pub(super) enum Flushed {
    /// A row written as it is.
    Single(usize),
    /// Rows `beg..beg+n` to join.
    Merge(usize, usize),
}

/// The variant type of a record as an integer, for ordering (bcftools compares the
/// `bcf_get_variant_types` values).
fn type_code(rec: &mut Record) -> u32 {
    rec.variant_types().bits()
}

impl Joiner {
    /// Whether a record is at a different position than the waiting ones (bcftools'
    /// `mrows_can_flush`).
    pub(super) fn can_flush(&self, rec: &Record) -> bool {
        if self.n == 0 {
            return false;
        }
        let first = &self.rows[self.first];
        first.rid() != rec.rid() || first.pos() != rec.pos()
    }

    /// Adds a record, keeping records of the same type together unless joining any (bcftools'
    /// `mrows_push`). Returns the record it replaces in the row storage, for reuse.
    pub(super) fn push(&mut self, rec: Record, collapse: Collapse) -> Option<Record> {
        if self.n == 0 {
            self.first = 0;
        }
        self.n += 1;
        let slot = self.n - 1;
        let old = if slot < self.rows.len() {
            Some(std::mem::replace(&mut self.rows[slot], rec))
        } else {
            self.rows.push(rec);
            None
        };
        if collapse == Collapse::Any {
            return old;
        }
        // Move the line up so that the same variant types end up together.
        let cur = type_code(&mut self.rows[slot]);
        let mut i = self.first + self.n - 1;
        while i > 0 {
            let prev = type_code(&mut self.rows[i - 1]);
            if prev <= cur {
                break;
            }
            self.rows.swap(i - 1, i);
            i -= 1;
        }
        old
    }

    /// Returns the next record or group of records to write (bcftools' `mrows_flush`).
    pub(super) fn flush(&mut self, collapse: Collapse) -> Option<Flushed> {
        if self.n == 0 {
            return None;
        }
        let beg = self.first;
        if self.n == 1 {
            self.n = 0;
            return Some(Flushed::Single(beg));
        }
        if collapse == Collapse::Any {
            let n = self.n - beg;
            self.n = 0;
            return Some(Flushed::Merge(beg, n));
        }
        // Join everything within the same category.
        let types: &[VariantType] = match collapse {
            Collapse::Snps => &[VariantType::SNP],
            Collapse::Indels => &[VariantType::INDEL],
            _ => &[
                VariantType::SNP,
                VariantType::MNP,
                VariantType::INDEL,
                VariantType::OTHER,
            ],
        };
        for &want in types {
            let mut i = beg;
            let mut ty = want;
            while i < beg + self.n {
                ty = self.rows[i].variant_types();
                if ty != want && !ty.is_empty() {
                    break;
                }
                i += 1;
            }
            if i == beg + 1 && !ty.is_empty() {
                // Just one line of this type, but lines of other types follow.
                self.n -= 1;
                self.first += 1;
                return Some(Flushed::Single(beg));
            }
            if i > beg {
                let n = i - beg;
                self.n -= n;
                self.first += n;
                return Some(Flushed::Merge(beg, n));
            }
        }
        self.n -= 1;
        self.first += 1;
        Some(Flushed::Single(beg))
    }
}

/// Merges the alleles of a record into the joined list `b`, filling `map` with each allele's
/// index in `b` (bcftools' `merge_alleles` from vcfmerge.c). Returns `false` when the REF
/// alleles disagree.
fn merge_alleles(a: &mut [Vec<u8>], map: &mut Vec<usize>, b: &mut Vec<Vec<u8>>) -> bool {
    map.clear();
    map.resize(a.len(), 0);
    let rla = if a[0].len() <= 1 { 1 } else { a[0].len() };
    let rlb = if b[0].len() <= 1 { 1 } else { b[0].len() };

    // The most common case: the same SNP.
    if a.len() == 2
        && b.len() == 2
        && rla == 1
        && rlb == 1
        && a[1].len() == 1
        && b[1].len() == 1
        && a[1][0] == b[1][0]
    {
        map[1] = 1;
        return true;
    }

    // REF prefixes must be identical.
    let m = rla.min(rlb);
    let prefix = |x: &[u8]| -> Vec<u8> {
        let mut v = x[..m.min(x.len())].to_vec();
        v.resize(m, 0);
        v
    };
    let (pa, pb) = (prefix(&a[0]), prefix(&b[0]));
    if pa != pb {
        if !pa.eq_ignore_ascii_case(&pb) {
            eprintln!(
                "The REF prefixes differ: {} vs {} ({rla},{rlb})",
                String::from_utf8_lossy(&a[0]),
                String::from_utf8_lossy(&b[0])
            );
            return false;
        }
        // Different case: change to uppercase.
        for x in a.iter_mut().chain(b.iter_mut()) {
            x.make_ascii_uppercase();
        }
    }

    // The b alleles need expanding.
    if rla > rlb {
        let tail = a[0][rlb..].to_vec();
        for x in b.iter_mut() {
            if matches!(x.first(), Some(b'<' | b'*')) {
                continue;
            }
            x.extend_from_slice(&tail);
        }
    }

    // Add the a alleles not present yet.
    for i in 1..a.len() {
        let ai: Vec<u8> = if rlb > rla && !matches!(a[i].first(), Some(b'<' | b'*')) {
            let mut v = a[i].clone();
            v.extend_from_slice(&b[0][rla..]);
            v
        } else {
            a[i].clone()
        };
        match (1..b.len()).find(|&j| ai.eq_ignore_ascii_case(&b[j])) {
            Some(j) => map[i] = j,
            None => {
                map[i] = b.len();
                b.push(ai);
            }
        }
    }
    true
}

/// Copies the per-allele fields of `src` to the joined alleles' places in `dst` (Number=A
/// from `jfrom` 1, Number=R from 0).
fn copy_allele_fields(src: &[u8], map: &[usize], jfrom: usize, dst: &mut Vec<u8>) {
    for (j, &m) in map.iter().enumerate().skip(jfrom) {
        copy_string_field(src, j - jfrom, dst, m - jfrom);
    }
}

/// Copies the per-genotype fields of `src` to the joined genotypes' places in `dst`.
fn copy_genotype_fields(src: &[u8], map: &[usize], dst: &mut Vec<u8>) {
    let mut kori = 0;
    for (iori, &inew) in map.iter().enumerate() {
        for &jnew in &map[..=iori] {
            copy_string_field(src, kori, dst, alleles2gt(inew, jnew));
            kori += 1;
        }
    }
}

/// Widens each of `nsmpl` runs of `n_ori` values to `n` values, padding with missing values
/// (bcftools' `ENLARGE_ARRAY`).
fn enlarge<T: Value>(vals: &mut Vec<T>, nsmpl: usize, n_ori: usize, n: usize) {
    let mut out = Vec::with_capacity(nsmpl * n);
    for s in 0..nsmpl {
        let src = &vals[(s * n_ori).min(vals.len())..((s + 1) * n_ori).min(vals.len())];
        out.extend_from_slice(src);
        out.extend(std::iter::repeat_n(T::missing(), n - src.len().min(n)));
    }
    *vals = out;
}

impl Normalizer {
    /// Joins `lines` into `dst` (bcftools' `merge_biallelics_to_multiallelic`).
    pub(super) fn merge_biallelics(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        dst: &mut Record,
        lines: &mut [Record],
    ) -> anyhow::Result<()> {
        for l in lines.iter_mut() {
            l.unpack(Unpack::ALL);
        }
        dst.clear();
        dst.set_rid(lines[0].rid());
        dst.set_pos(lines[0].pos());

        // QUAL: the maximum.
        let mut qual: Option<f32> = None;
        for l in lines.iter() {
            if let Some(q) = l.qual()
                && qual.is_none_or(|m| m < q)
            {
                qual = Some(q);
            }
        }
        dst.set_qual(qual);
        dst.update_id(Some(lines[0].id()));

        // Merge the alleles, mapping each line's alleles to the joined ones.
        let maps = &mut self.joiner.maps;
        if maps.len() < lines.len() {
            maps.resize_with(lines.len(), Vec::new);
        }
        let mut als: Vec<Vec<u8>> = lines[0].alleles().map(<[u8]>::to_vec).collect();
        maps[0].clear();
        maps[0].extend(0..als.len());
        for i in 1..lines.len() {
            let id = lines[i].id();
            if id != b"." {
                let id = id.to_vec();
                dst.add_id(&id);
            }
            let mut a: Vec<Vec<u8>> = lines[i].alleles().map(<[u8]>::to_vec).collect();
            if !merge_alleles(&mut a, &mut maps[i], &mut als) {
                bail!(
                    "Failed to merge alleles at {}:{}",
                    String::from_utf8_lossy(dst.chrom_safe(hdr)),
                    dst.pos() + 1
                );
            }
        }
        dst.update_alleles(out_hdr, &als);
        let nals = als.len();

        if !lines[0].filters().is_empty() {
            dst.update_filter(lines[0].filters());
        }
        for l in lines.iter().skip(1) {
            for &f in l.filters() {
                if f == 0 {
                    // PASS
                    if self.opts.strict_filter {
                        dst.update_filter(l.filters());
                        break;
                    }
                    continue;
                }
                dst.add_filter(f);
            }
        }

        // INFO
        let info_keys: Vec<i32> = lines[0].info_fields().iter().map(|f| f.key()).collect();
        for key in info_keys {
            let tag = hdr.id_name(key as u32).unwrap_or_default();
            match hdr.id_type(HeaderRecordType::Info, key as u32) {
                Some(ValueType::Integer) => {
                    self.merge_info_numeric::<i32>(hdr, out_hdr, lines, key, tag, dst, nals)?;
                }
                Some(ValueType::Float) => {
                    self.merge_info_numeric::<f32>(hdr, out_hdr, lines, key, tag, dst, nals)?;
                }
                Some(ValueType::Flag) => {
                    let set = lines[0].get_info_flag(hdr, tag).unwrap_or(false);
                    let _ = dst.update_info_flag(out_hdr, tag, set);
                }
                _ => self.merge_info_string(hdr, out_hdr, lines, key, tag, dst, nals),
            }
        }

        // FORMAT
        let gt_id = hdr.id(b"GT").map(|i| i as i32);
        dst.set_sample_count(lines[0].sample_count());
        let fmt_ids: Vec<i32> = lines[0].format_fields().iter().map(|f| f.id()).collect();
        for id in fmt_ids {
            let tag = hdr.id_name(id as u32).unwrap_or_default();
            if Some(id) == gt_id {
                self.merge_format_genotype(hdr, out_hdr, lines, dst)?;
                continue;
            }
            match hdr.id_type(HeaderRecordType::Format, id as u32) {
                Some(ValueType::Integer) => {
                    self.merge_format_numeric::<i32>(hdr, out_hdr, lines, id, tag, dst, nals)?;
                }
                Some(ValueType::Float) => {
                    self.merge_format_numeric::<f32>(hdr, out_hdr, lines, id, tag, dst, nals)?;
                }
                _ => self.merge_format_string(hdr, out_hdr, lines, id, tag, dst, nals)?,
            }
        }
        dst.sync();
        self.stats.joined += 1;
        Ok(())
    }

    /// bcftools' `merge_info_numeric`.
    #[expect(clippy::too_many_arguments, reason = "bcftools' merge_info_numeric")]
    fn merge_info_numeric<T: Value>(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        lines: &mut [Record],
        key: i32,
        tag: &[u8],
        dst: &mut Record,
        nals: usize,
    ) -> anyhow::Result<()> {
        let mut vals = std::mem::take(&mut T::buffers(&mut self.scratch).a);
        let mut vals2 = std::mem::take(&mut T::buffers(&mut self.scratch).b);
        let r = (|| -> anyhow::Result<()> {
            let Ok(nvals_ori) = T::get_info(&mut lines[0], hdr, tag, &mut vals) else {
                return Ok(());
            };
            let maps = &self.joiner.maps;
            let chrom = |l: &Record| String::from_utf8_lossy(l.chrom_safe(hdr)).into_owned();
            let n0 = lines[0].allele_count();
            let number = hdr.id_number(HeaderRecordType::Info, key as u32);
            let first_error = |expected: usize| {
                format!(
                    "vcfnorm: number of fields in first record at position {}:{} for INFO tag {} not as expected [found: {nvals_ori} vs expected:{expected}]",
                    chrom(&lines[0]),
                    lines[0].pos() + 1,
                    String::from_utf8_lossy(tag)
                )
            };
            let merge_error = |l: &Record| {
                format!(
                    "vcfnorm: could not merge INFO tag {} at position {}:{}",
                    String::from_utf8_lossy(tag),
                    chrom(l),
                    l.pos() + 1
                )
            };
            match number {
                Some(Number::A | Number::R) => {
                    let is_a = number == Some(Number::A);
                    let off = usize::from(is_a);
                    if nvals_ori != n0 - off {
                        bail!("{}", first_error(n0 - off));
                    }
                    let nvals = nals - off;
                    enlarge(&mut vals, 1, nvals_ori, nvals);
                    for i in 1..lines.len() {
                        let Ok(nvals2) = T::get_info(&mut lines[i], hdr, tag, &mut vals2) else {
                            continue;
                        };
                        if nvals2 != lines[i].allele_count() - off {
                            bail!("{}", merge_error(&lines[i]));
                        }
                        for (k, &v) in vals2.iter().enumerate().take(nvals2) {
                            if v.is_vector_end() {
                                break;
                            }
                            vals[maps[i][k + off] - off] = v;
                        }
                    }
                    let _ = T::update_info(dst, out_hdr, tag, &vals[..nvals]);
                }
                Some(Number::G) => {
                    if nvals_ori != n0 * (n0 + 1) / 2 {
                        eprintln!("todo: merge Number=G INFO fields for haploid sites");
                        bail!("{}", first_error(n0 * (n0 + 1) / 2));
                    }
                    let nvals = nals * (nals + 1) / 2;
                    enlarge(&mut vals, 1, nvals_ori, nvals);
                    for i in 1..lines.len() {
                        let Ok(nvals2) = T::get_info(&mut lines[i], hdr, tag, &mut vals2) else {
                            continue;
                        };
                        let ni = lines[i].allele_count();
                        if nvals2 != ni * (ni + 1) / 2 {
                            bail!("{}", merge_error(&lines[i]));
                        }
                        let mut k = 0;
                        for ia in 0..ni {
                            for ib in 0..=ia {
                                if vals2[k].is_vector_end() {
                                    break;
                                }
                                vals[alleles2gt(maps[i][ia], maps[i][ib])] = vals2[k];
                                k += 1;
                            }
                        }
                    }
                    let _ = T::update_info(dst, out_hdr, tag, &vals[..nvals]);
                }
                _ => {
                    let _ = T::update_info(dst, out_hdr, tag, &vals[..nvals_ori]);
                }
            }
            Ok(())
        })();
        T::buffers(&mut self.scratch).a = vals;
        T::buffers(&mut self.scratch).b = vals2;
        r
    }

    /// bcftools' `merge_info_string`.
    #[expect(clippy::too_many_arguments, reason = "bcftools' merge_info_string")]
    fn merge_info_string(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        lines: &mut [Record],
        key: i32,
        tag: &[u8],
        dst: &mut Record,
        nals: usize,
    ) {
        let maps = &self.joiner.maps;
        let s = &mut self.scratch.bytes;
        s.clear();
        let number = hdr.id_number(HeaderRecordType::Info, key as u32);
        match number {
            Some(Number::A | Number::R) => {
                let jfrom = usize::from(number == Some(Number::A));
                s.push(b'.');
                for _ in jfrom + 1..nals {
                    s.extend_from_slice(b",.");
                }
                for (i, l) in lines.iter_mut().enumerate() {
                    let n = l.allele_count();
                    let Some(field) = l.info_field_by_id(key).cloned() else {
                        continue;
                    };
                    let src = &l.info_data(&field)[..field.len().min(l.info_data(&field).len())];
                    copy_allele_fields(src, &maps[i][..n], jfrom, s);
                }
            }
            Some(Number::G) => {
                let ngts = nals * (nals + 1) / 2;
                s.push(b'.');
                for _ in 1..ngts {
                    s.extend_from_slice(b",.");
                }
                for (i, l) in lines.iter_mut().enumerate() {
                    let n = l.allele_count();
                    let Some(field) = l.info_field_by_id(key).cloned() else {
                        continue;
                    };
                    let src = &l.info_data(&field)[..field.len().min(l.info_data(&field).len())];
                    copy_genotype_fields(src, &maps[i][..n], s);
                }
            }
            _ => {
                if let Ok(v) = lines[0].get_info_string(hdr, tag) {
                    s.extend_from_slice(v);
                }
            }
        }
        let _ = dst.update_info_string(out_hdr, tag, Some(s));
    }

    /// bcftools' `merge_format_genotype`.
    fn merge_format_genotype(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        lines: &mut [Record],
        dst: &mut Record,
    ) -> anyhow::Result<()> {
        let nsmpl = hdr.sample_count();
        let mut gt = std::mem::take(&mut self.scratch.i.a);
        let mut gt2 = std::mem::take(&mut self.scratch.i.b);
        let r = (|| -> anyhow::Result<()> {
            let Ok(n) = lines[0].get_genotypes(hdr, &mut gt) else {
                return Ok(());
            };
            if n == 0 || nsmpl == 0 {
                return Ok(());
            }
            let mut ngts = n / nsmpl;
            let maps = &self.joiner.maps;
            for i in 1..lines.len() {
                let ngts2 = lines[i]
                    .get_genotypes(hdr, &mut gt2)
                    .map_or(0, |n| n / nsmpl);
                let ploidy_changed = ngts != ngts2;
                if ngts < ngts2 {
                    grow_ploidy(&mut gt, ngts, ngts2, nsmpl);
                    ngts = ngts2;
                }
                if ngts > ngts2 {
                    grow_ploidy(&mut gt2, ngts2, ngts, nsmpl);
                }
                let nals_i = lines[i].allele_count();
                for j in 0..nsmpl {
                    let g = &mut gt[j * ngts..(j + 1) * ngts];
                    let g2 = &gt2[j * ngts..(j + 1) * ngts];
                    for k2 in 0..ngts {
                        let v2 = Genotype(g2[k2]);
                        if v2.is_vector_end() {
                            if ploidy_changed && Genotype(g[k2]).is_missing() {
                                g[k2] = pasta::variant::encoding::INT32_VECTOR_END;
                            }
                            break;
                        }
                        if v2.is_missing() {
                            continue;
                        }
                        // Never overwrite with REF, unless the destination is missing.
                        let ial2 = v2.allele();
                        let v = Genotype(g[k2]);
                        if ial2 == 0 && !v.is_missing() && !v.is_vector_end() {
                            continue;
                        }
                        if ial2 < 0 || ial2 as usize >= nals_i {
                            bail!(
                                "Error at {}:{}: incorrect allele index {ial2}",
                                String::from_utf8_lossy(lines[i].chrom_safe(hdr)),
                                lines[i].pos() + 1
                            );
                        }
                        let ial = maps[i][ial2 as usize] as i32;
                        if v.is_vector_end() || v.is_missing() || v.allele() == 0 {
                            g[k2] = if !v.is_vector_end() && v.is_phased() {
                                Genotype::phased(ial).0
                            } else {
                                Genotype::unphased(ial).0
                            };
                        } else if let Some(k) = g.iter().position(|&x| {
                            let x = Genotype(x);
                            x.is_vector_end() || x.is_missing() || x.allele() == 0
                        }) {
                            // A conflict: the first line has an ALT here.
                            g[k] = Genotype::unphased(ial).0;
                        }
                    }
                }
            }
            let _ = dst.update_genotypes(out_hdr, &gt[..ngts * nsmpl]);
            Ok(())
        })();
        self.scratch.i.a = gt;
        self.scratch.i.b = gt2;
        r
    }

    /// bcftools' `merge_format_numeric`.
    #[expect(clippy::too_many_arguments, reason = "bcftools' merge_format_numeric")]
    fn merge_format_numeric<T: Value>(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        lines: &mut [Record],
        id: i32,
        tag: &[u8],
        dst: &mut Record,
        nals: usize,
    ) -> anyhow::Result<()> {
        let nsmpl = hdr.sample_count();
        let mut vals = std::mem::take(&mut T::buffers(&mut self.scratch).a);
        let mut vals2 = std::mem::take(&mut T::buffers(&mut self.scratch).b);
        let r = (|| -> anyhow::Result<()> {
            let Ok(n) = T::get_format(&mut lines[0], hdr, tag, &mut vals) else {
                return Ok(());
            };
            if n == 0 || nsmpl == 0 {
                return Ok(());
            }
            let nvals_ori = n / nsmpl;
            let maps = &self.joiner.maps;
            let merge_error = |l: &Record| {
                format!(
                    "vcfnorm: could not merge FORMAT tag {} at position {}:{}",
                    String::from_utf8_lossy(tag),
                    String::from_utf8_lossy(l.chrom_safe(hdr)),
                    l.pos() + 1
                )
            };
            let number = hdr.id_number(HeaderRecordType::Format, id as u32);
            match number {
                Some(Number::A | Number::R) => {
                    let off = usize::from(number == Some(Number::A));
                    let nvals = nals - off;
                    enlarge(&mut vals, nsmpl, nvals_ori, nvals);
                    for i in 1..lines.len() {
                        let Ok(n2) = T::get_format(&mut lines[i], hdr, tag, &mut vals2) else {
                            continue;
                        };
                        let nvals2 = n2 / nsmpl;
                        if nvals2 != lines[i].allele_count() - off {
                            bail!("{}", merge_error(&lines[i]));
                        }
                        for j in 0..nsmpl {
                            let d = &mut vals[j * nvals..(j + 1) * nvals];
                            for (k, &v) in vals2[j * nvals2..(j + 1) * nvals2].iter().enumerate() {
                                if v.is_vector_end() {
                                    break;
                                }
                                d[maps[i][k + off] - off] = v;
                            }
                        }
                    }
                    let _ = T::update_format(dst, out_hdr, tag, &vals[..nvals * nsmpl]);
                }
                Some(Number::G) => {
                    let n0 = lines[0].allele_count();
                    if nvals_ori > n0 * (n0 + 1) / 2 && ploidy_of(n0, nvals_ori).is_some() {
                        return join_any_ploidy(
                            hdr, out_hdr, lines, tag, dst, nals, maps, &mut vals2,
                        );
                    }
                    let mut diploid = vec![false; nsmpl];
                    let mut all_haploid = true;
                    if nvals_ori > n0 {
                        // The line may be diploid.
                        let ndiploid = n0 * (n0 + 1) / 2;
                        for (s, d) in diploid.iter_mut().enumerate() {
                            let sv = &vals[s * nvals_ori..(s + 1) * nvals_ori];
                            let k = sv
                                .iter()
                                .position(|v| v.is_vector_end())
                                .unwrap_or(nvals_ori);
                            if k == ndiploid {
                                *d = true;
                                all_haploid = false;
                            }
                        }
                    }
                    let mut nvals = nals * (nals + 1) / 2;
                    enlarge(&mut vals, nsmpl, nvals_ori, nvals);
                    for i in 1..lines.len() {
                        let Ok(n2) = T::get_format(&mut lines[i], hdr, tag, &mut vals2) else {
                            continue;
                        };
                        let nvals2 = n2 / nsmpl;
                        let ni = lines[i].allele_count();
                        let ndiploid = ni * (ni + 1) / 2;
                        let line_diploid = nvals2 == ndiploid;
                        if !(nvals2 == 1 || nvals2 == ni || nvals2 == ndiploid) {
                            if nvals2 > ndiploid && ploidy_of(ni, nvals2).is_some() {
                                return join_any_ploidy(
                                    hdr, out_hdr, lines, tag, dst, nals, maps, &mut vals2,
                                );
                            }
                            bail!("{}", merge_error(&lines[i]));
                        }
                        for j in 0..nsmpl {
                            let sv = &vals2[j * nvals2..(j + 1) * nvals2];
                            let mut smpl_diploid = line_diploid;
                            if smpl_diploid {
                                let k = sv.iter().position(|v| v.is_vector_end()).unwrap_or(nvals2);
                                if k != ndiploid {
                                    smpl_diploid = false;
                                }
                            }
                            if smpl_diploid && !diploid[j] {
                                diploid[j] = true;
                                all_haploid = false;
                            }
                            let d = &mut vals[j * nvals..(j + 1) * nvals];
                            if smpl_diploid {
                                let mut k = 0;
                                for ia in 0..ni {
                                    for ib in 0..=ia {
                                        d[alleles2gt(maps[i][ia], maps[i][ib])] = sv[k];
                                        k += 1;
                                    }
                                }
                            } else {
                                // bcftools reads n_allele values here, past a shorter row.
                                let flat = &vals2[j * nvals2..];
                                for k in 0..ni {
                                    d[maps[i][k]] = flat.get(k).copied().unwrap_or(T::missing());
                                }
                            }
                        }
                    }
                    if all_haploid {
                        // diploid_to_haploid
                        for s in 0..nsmpl {
                            vals.copy_within(s * nvals..s * nvals + nals, s * nals);
                        }
                        nvals = nals;
                    } else {
                        for (s, &d) in diploid.iter().enumerate() {
                            if !d {
                                vals[s * nvals + nals] = T::vector_end();
                            }
                        }
                    }
                    let _ = T::update_format(dst, out_hdr, tag, &vals[..nvals * nsmpl]);
                }
                _ => {
                    let _ = T::update_format(dst, out_hdr, tag, &vals[..nvals_ori * nsmpl]);
                }
            }
            Ok(())
        })();
        T::buffers(&mut self.scratch).a = vals;
        T::buffers(&mut self.scratch).b = vals2;
        r
    }

    /// bcftools' `merge_format_string`.
    #[expect(clippy::too_many_arguments, reason = "bcftools' merge_format_string")]
    fn merge_format_string(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        lines: &mut [Record],
        id: i32,
        tag: &[u8],
        dst: &mut Record,
        nals: usize,
    ) -> anyhow::Result<()> {
        let nsmpl = hdr.sample_count();
        let number = hdr.id_number(HeaderRecordType::Format, id as u32);
        if !matches!(number, Some(Number::A | Number::R | Number::G)) {
            if let Ok((data, _)) = lines[0].get_format_chars(hdr, tag) {
                let data = data.to_vec();
                let _ = dst.update_format_chars(out_hdr, tag, &data);
            }
            return Ok(());
        }
        if nsmpl == 0 {
            return Ok(());
        }
        let maps = &self.joiner.maps;
        let strs = &mut self.scratch.strings;
        if strs.len() < nsmpl {
            strs.resize_with(nsmpl, Vec::new);
        }
        for s in strs.iter_mut().take(nsmpl) {
            s.clear();
        }
        if number == Some(Number::G) {
            let n0 = lines[0].allele_count();
            let Ok((data, nret)) = lines[0].get_format_chars(hdr, tag) else {
                return Ok(());
            };
            let data = data.to_vec();
            // Each sample's ploidy, from the first line as in bcftools (which knows 1 and 2).
            let mut ploidy = vec![2; nsmpl];
            for i in 0..nsmpl {
                let src = &data[i * nret..(i + 1) * nret];
                let end = src.iter().position(|&c| c == 0).unwrap_or(nret);
                let mut nfields = 1 + src[..end].iter().filter(|&&c| c == b',').count();
                if nfields == n0 {
                    ploidy[i] = 1;
                    nfields = nals;
                } else if nfields == n0 * (n0 + 1) / 2 {
                    nfields = nals * (nals + 1) / 2;
                } else if let Some(p) = ploidy_of(n0, nfields)
                    // Each value takes at least `,.`, in a row as wide as the widest sample's.
                    && let Some(n) = genotype_count(nals, p).filter(|&n| {
                        n.saturating_mul(2).saturating_mul(nsmpl) <= MAX_JOINED_STRING_BYTES
                    })
                {
                    ploidy[i] = p;
                    nfields = n;
                } else {
                    bail!(
                        "The field {} at {}:{} neither diploid nor haploid?",
                        String::from_utf8_lossy(tag),
                        String::from_utf8_lossy(dst.chrom_safe(hdr)),
                        dst.pos() + 1
                    );
                }
                strs[i].push(b'.');
                for _ in 1..nfields {
                    strs[i].extend_from_slice(b",.");
                }
            }
            let merge_error = |l: &Record| {
                anyhow::anyhow!(
                    "vcfnorm: could not merge FORMAT tag {} at position {}:{}",
                    String::from_utf8_lossy(tag),
                    String::from_utf8_lossy(l.chrom_safe(hdr)),
                    l.pos() + 1
                )
            };
            // As for numbers, when a sample has values of a ploidy above two in some line, each
            // sample's values in every line must be of its ploidy in the first (or missing):
            // others cannot be placed. Without such a sample the joins are bcftools'.
            let polyploid = ploidy.iter().any(|&p| p > 2)
                || lines[1..].iter_mut().any(|l| {
                    let n = l.allele_count();
                    l.get_format_chars(hdr, tag).is_ok_and(|(data, nret)| {
                        (0..nsmpl).any(|k| {
                            let count = string_value_count(&data[k * nret..(k + 1) * nret]);
                            count > n * (n + 1) / 2 && ploidy_of(n, count).is_some()
                        })
                    })
                });
            if polyploid {
                for l in lines.iter_mut() {
                    let n = l.allele_count();
                    let Ok((data, nret)) = l.get_format_chars(hdr, tag) else {
                        continue;
                    };
                    let fits = ploidy.iter().enumerate().all(|(k, &p)| {
                        let count = string_value_count(&data[k * nret..(k + 1) * nret]);
                        count == 0 || genotype_count(n, p) == Some(count)
                    });
                    if !fits {
                        return Err(merge_error(l));
                    }
                }
            }
            let (mut gmap, mut picks, mut joined) = (Vec::new(), Vec::new(), Vec::new());
            for (i, l) in lines.iter_mut().enumerate() {
                let n = l.allele_count();
                let Ok((data, nret)) = l.get_format_chars(hdr, tag) else {
                    continue;
                };
                for (k, s) in strs.iter_mut().enumerate().take(nsmpl) {
                    let src = &data[k * nret..(k + 1) * nret];
                    match ploidy[k] {
                        1 => copy_allele_fields(src, &maps[i][..n], 0, s),
                        2 => copy_genotype_fields(src, &maps[i][..n], s),
                        p => {
                            genotype_map(&maps[i][..n], p, &mut gmap);
                            copy_genotype_strings(src, &gmap, s, &mut picks, &mut joined);
                            // Values past the limit stop too.
                            if s.len().saturating_mul(nsmpl) > MAX_JOINED_STRING_BYTES {
                                return Err(merge_error(l));
                            }
                        }
                    }
                }
            }
        } else {
            let jfrom = usize::from(number == Some(Number::A));
            for s in strs.iter_mut().take(nsmpl) {
                s.push(b'.');
                for _ in jfrom + 1..nals {
                    s.extend_from_slice(b",.");
                }
            }
            for (i, l) in lines.iter_mut().enumerate() {
                let n = l.allele_count();
                let Ok((data, nret)) = l.get_format_chars(hdr, tag) else {
                    continue;
                };
                for (k, s) in strs.iter_mut().enumerate().take(nsmpl) {
                    let src = &data[k * nret..(k + 1) * nret];
                    copy_allele_fields(src, &maps[i][..n], jfrom, s);
                }
            }
        }
        let max_len = strs.iter().take(nsmpl).map(Vec::len).max().unwrap_or(0);
        let out = &mut self.scratch.bytes;
        out.clear();
        for s in strs.iter().take(nsmpl) {
            out.extend_from_slice(s);
            out.extend(std::iter::repeat_n(0u8, max_len - s.len()));
        }
        let _ = dst.update_format_chars(out_hdr, tag, out);
        Ok(())
    }
}

/// The number of comma-separated values in one sample's string `src`, up to its NUL padding:
/// none for an empty string or a single `.`.
fn string_value_count(src: &[u8]) -> usize {
    let end = src.iter().position(|&c| c == 0).unwrap_or(src.len());
    match &src[..end] {
        b"" | b"." => 0,
        values => 1 + values.iter().filter(|&&c| c == b',').count(),
    }
}

/// Joins a Number=G FORMAT field whose values are of a ploidy above two in some sample, which
/// bcftools cannot do: it sizes the values for haploid or diploid samples only, and stops or
/// writes them at the wrong places. Each sample's ploidy comes from the number of its values
/// in the lines that have any (a single `.` has none). The joined values are missing where no
/// line has the genotype, and a later line's value of a genotype replaces an earlier one's, as
/// for diploid samples; a sample without values gets one missing value.
#[expect(
    clippy::too_many_arguments,
    reason = "the state of merge_format_numeric"
)]
fn join_any_ploidy<T: Value>(
    hdr: &Header,
    out_hdr: &Header,
    lines: &mut [Record],
    tag: &[u8],
    dst: &mut Record,
    nals: usize,
    maps: &[Vec<usize>],
    vals: &mut Vec<T>,
) -> anyhow::Result<()> {
    let nsmpl = hdr.sample_count();
    let merge_error = |l: &Record| {
        anyhow::anyhow!(
            "vcfnorm: could not merge FORMAT tag {} at position {}:{}",
            String::from_utf8_lossy(tag),
            String::from_utf8_lossy(l.chrom_safe(hdr)),
            l.pos() + 1
        )
    };
    // The number of values in a sample's row: up to the vector end, none for '.'.
    let row_len = |row: &[T]| {
        let k = row
            .iter()
            .position(|v| v.is_vector_end())
            .unwrap_or(row.len());
        if k == 1 && row[0].is_missing() { 0 } else { k }
    };
    // Each sample's ploidy (0 while unknown), the same in every line.
    let mut ploidy = vec![0; nsmpl];
    for l in lines.iter_mut() {
        let Ok(n) = T::get_format(l, hdr, tag, vals) else {
            continue;
        };
        let (stride, ni) = (n / nsmpl, l.allele_count());
        for (s, p) in ploidy.iter_mut().enumerate() {
            let k = row_len(&vals[s * stride..(s + 1) * stride]);
            if k == 0 || (ni < 2 && k == 1) {
                continue;
            }
            match ploidy_of(ni, k) {
                Some(q) if *p == 0 || *p == q => *p = q,
                _ => return Err(merge_error(l)),
            }
        }
    }
    let mut counts = Vec::with_capacity(nsmpl);
    for &p in &ploidy {
        let c = if p == 0 {
            Some(1)
        } else {
            genotype_count(nals, p)
        };
        counts.push(c.ok_or_else(|| merge_error(&lines[0]))?);
    }
    let width = counts.iter().copied().max().unwrap_or(1);
    if width
        .checked_mul(nsmpl)
        .and_then(|n| n.checked_mul(size_of::<T>()))
        .is_none_or(|bytes| bytes > MAX_JOINED_FORMAT_BYTES)
    {
        return Err(merge_error(&lines[0]));
    }
    let mut out = Vec::with_capacity(width * nsmpl);
    for &c in &counts {
        out.extend(std::iter::repeat_n(T::missing(), c));
        out.extend(std::iter::repeat_n(T::vector_end(), width - c));
    }
    // The joined index of each genotype of each line, by ploidy.
    let mut gmaps: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, l) in lines.iter_mut().enumerate() {
        let Ok(n) = T::get_format(l, hdr, tag, vals) else {
            continue;
        };
        let (stride, ni) = (n / nsmpl, l.allele_count());
        gmaps.clear();
        for (s, &p) in ploidy.iter().enumerate() {
            let row = &vals[s * stride..(s + 1) * stride];
            let k = row_len(row);
            if k == 0 || p == 0 {
                continue;
            }
            let d = &mut out[s * width..(s + 1) * width];
            if ni < 2 {
                // Only REF: its genotype is the first of every ploidy.
                d[0] = row[0];
                continue;
            }
            let g = match gmaps.iter().position(|(q, _)| *q == p) {
                Some(g) => g,
                None => {
                    let mut m = Vec::new();
                    genotype_map(&maps[i][..ni], p, &mut m);
                    gmaps.push((p, m));
                    gmaps.len() - 1
                }
            };
            for (&j, &v) in gmaps[g].1.iter().zip(&row[..k]) {
                d[j] = v;
            }
        }
    }
    let _ = T::update_format(dst, out_hdr, tag, &out);
    Ok(())
}

/// Copies field `j` of the comma-separated `src` (one sample's values in a line) into field
/// `gmap[j]` of `dst` where that is still missing (`.`), as [`copy_string_field`] does field by
/// field, in one pass over `dst`: the rows of high ploidies are long. `picks` and `joined` are
/// scratch buffers.
fn copy_genotype_strings(
    src: &[u8],
    gmap: &[usize],
    dst: &mut Vec<u8>,
    picks: &mut Vec<(usize, usize, usize)>,
    joined: &mut Vec<u8>,
) {
    let end = src.iter().position(|&c| c == 0).unwrap_or(src.len());
    picks.clear();
    let mut start = 0;
    for (field, &g) in src[..end].split(|&c| c == b',').zip(gmap) {
        if field != b"." {
            picks.push((g, start, start + field.len()));
        }
        start += field.len() + 1;
    }
    if picks.is_empty() {
        return;
    }
    picks.sort_unstable();
    joined.clear();
    let mut next = picks.iter().peekable();
    for (g, field) in dst.split(|&c| c == b',').enumerate() {
        if g > 0 {
            joined.push(b',');
        }
        match next.next_if(|p| p.0 == g) {
            Some(&(_, s, e)) if field == b"." => joined.extend_from_slice(&src[s..e]),
            _ => joined.extend_from_slice(field),
        }
    }
    std::mem::swap(dst, joined);
}

/// Widens genotype rows from `ngt_ori` to `ngt_new` alleles, padding with vector ends
/// (bcftools' `gt_array_grow_ploidy`).
fn grow_ploidy(gt: &mut Vec<i32>, ngt_ori: usize, ngt_new: usize, nsmpl: usize) {
    let mut out = Vec::with_capacity(ngt_new * nsmpl);
    for s in 0..nsmpl {
        let src = &gt[(s * ngt_ori).min(gt.len())..((s + 1) * ngt_ori).min(gt.len())];
        out.extend_from_slice(src);
        out.extend(std::iter::repeat_n(
            pasta::variant::encoding::INT32_VECTOR_END,
            ngt_new - src.len(),
        ));
    }
    *gt = out;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn merge(a: &[&str], b: &[&str]) -> (Vec<usize>, Vec<String>) {
        let mut a: Vec<Vec<u8>> = a.iter().map(|s| s.as_bytes().to_vec()).collect();
        let mut b: Vec<Vec<u8>> = b.iter().map(|s| s.as_bytes().to_vec()).collect();
        let mut map = Vec::new();
        assert!(merge_alleles(&mut a, &mut map, &mut b));
        let b = b
            .iter()
            .map(|x| String::from_utf8_lossy(x).into_owned())
            .collect();
        (map, b)
    }

    #[test]
    fn alleles_merge_as_vcfmerge_does() {
        assert_eq!(
            merge(&["A", "C"], &["A", "C"]),
            (vec![0, 1], vec!["A".into(), "C".into()])
        );
        assert_eq!(
            merge(&["A", "G"], &["A", "C"]),
            (vec![0, 2], vec!["A".into(), "C".into(), "G".into()])
        );
        // The shorter REF is extended.
        assert_eq!(
            merge(&["AT", "A"], &["A", "C"]),
            (vec![0, 2], vec!["AT".into(), "CT".into(), "A".into()])
        );
        assert_eq!(
            merge(&["A", "G", "<DEL>"], &["AT", "A"]),
            (
                vec![0, 2, 3],
                vec!["AT".into(), "A".into(), "GT".into(), "<DEL>".into()]
            )
        );
        let mut a = vec![b"C".to_vec(), b"T".to_vec()];
        let mut b = vec![b"A".to_vec(), b"G".to_vec(), b"T".to_vec()];
        assert!(!merge_alleles(&mut a, &mut Vec::new(), &mut b));
    }

    #[test]
    fn genotype_strings_are_copied_as_field_by_field() {
        let cases: [(&[u8], &[usize], &str); 5] = [
            (b"a,b,c,d\0\0", &[0, 4, 7, 9], ".,x,.,.,.,.,.,.,.,."),
            (b"a,.,c,d", &[9, 4, 7, 0], ".,.,.,.,.,.,.,.,.,."),
            (b".", &[0, 1, 2, 3], ".,.,.,."),
            (b"a,bb,,d", &[3, 1, 0, 5], "x,.,.,.,."),
            (b"a,b", &[1, 0, 2, 3], "."),
        ];
        let (mut picks, mut joined) = (Vec::new(), Vec::new());
        for (src, gmap, dst) in cases {
            let mut want = dst.as_bytes().to_vec();
            for (j, &g) in gmap.iter().enumerate() {
                copy_string_field(src, j, &mut want, g);
            }
            let mut got = dst.as_bytes().to_vec();
            copy_genotype_strings(src, gmap, &mut got, &mut picks, &mut joined);
            assert_eq!(
                String::from_utf8_lossy(&got),
                String::from_utf8_lossy(&want),
                "{gmap:?}"
            );
        }
    }
}
