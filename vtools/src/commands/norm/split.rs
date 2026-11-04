//! Splitting multiallelic records into biallelic ones (`-m-`; bcftools'
//! `split_multiallelic_to_biallelics` and the `split_*` functions).

use anyhow::bail;
use pasta::variant::{Genotype, Header, HeaderRecordType, Number, Record, Unpack, ValueType};

use super::{
    Normalizer, alleles2gt,
    ploidy::{ploidy_of, ref_alt_index},
    values::{Value, nth_field},
};

/// How the INFO/FORMAT value count is wrong, for the `--force` messages.
struct Mismatch<'a> {
    kind: &'a str,
    tag: &'a [u8],
    chrom: &'a [u8],
    pos: i64,
    expected: usize,
    found: usize,
}

impl Normalizer {
    /// Handles a field with the wrong number of values: with `--force`, warns once and
    /// returns `Ok` (the caller drops the field); otherwise fails with `error`.
    fn wrong_number(
        &mut self,
        m: &Mismatch<'_>,
        suffix: &str,
        error: String,
    ) -> anyhow::Result<()> {
        if !self.opts.force {
            bail!("{error}");
        }
        if !self.force_warned {
            eprintln!(
                "Warning: wrong number of fields in {}/{} at {}:{}, expected {}, found {}{suffix}\n         (This warning is printed only once.)",
                m.kind,
                String::from_utf8_lossy(m.tag),
                String::from_utf8_lossy(m.chrom),
                m.pos + 1,
                m.expected,
                m.found
            );
            self.force_warned = true;
        }
        Ok(())
    }

    /// Splits a multiallelic record into `self.split_lines`, one record per ALT.
    ///
    /// Each INFO/FORMAT tag is decoded once and written to every ALT's record, so that the
    /// work grows with the size of the record rather than with the ALT count times it (bcftools
    /// decodes every tag again for each ALT). Each record still receives its fields in the
    /// source's order.
    pub(super) fn split_multiallelic(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        line: &mut Record,
    ) -> anyhow::Result<()> {
        line.unpack(Unpack::ALL);
        let nalt = line.allele_count() - 1;
        let gt_id = hdr.id(b"GT").map(|i| i as i32);
        let mut out = std::mem::take(&mut self.split_lines);
        out.clear();
        let mut alleles = std::mem::take(&mut self.scratch.bytes2);
        alleles.clear();
        alleles.extend_from_slice(line.allele(0));
        alleles.push(b',');
        let rlen = alleles.len();
        for i in 0..nalt {
            let mut dst = self.take_record();
            dst.set_rid(line.rid());
            dst.set_pos(line.pos());
            dst.set_qual(line.qual());
            // IDs could belong to one ALT; they are kept for all.
            dst.update_id(Some(line.id()));
            alleles.truncate(rlen);
            alleles.extend_from_slice(line.allele(i + 1));
            dst.update_alleles_str(out_hdr, &alleles);
            if !line.filters().is_empty() {
                dst.update_filter(line.filters());
            }
            out.push(dst);
        }
        self.scratch.bytes2 = alleles;
        let mut info_keys = std::mem::take(&mut self.scratch.info_keys);
        info_keys.clear();
        info_keys.extend(line.info_fields().iter().map(|f| f.key()));
        for &key in &info_keys {
            let tag = hdr.id_name(key as u32).unwrap_or_default();
            match hdr.id_type(HeaderRecordType::Info, key as u32) {
                Some(ValueType::Integer) => {
                    self.split_info_numeric::<i32>(hdr, out_hdr, line, key, tag, &mut out)?;
                }
                Some(ValueType::Float) => {
                    self.split_info_numeric::<f32>(hdr, out_hdr, line, key, tag, &mut out)?;
                }
                Some(ValueType::Flag) => {
                    let set = line.get_info_flag(hdr, tag).unwrap_or(false);
                    for dst in &mut out {
                        let _ = dst.update_info_flag(out_hdr, tag, set);
                    }
                }
                _ => self.split_info_string(hdr, out_hdr, line, key, tag, &mut out),
            }
        }
        self.scratch.info_keys = info_keys;
        for (i, dst) in out.iter_mut().enumerate() {
            self.old_rec_tag_init(hdr, line);
            self.old_rec_tag_set(out_hdr, dst, i + 1)?;
            dst.set_sample_count(line.sample_count());
        }
        let mut fmt_ids = std::mem::take(&mut self.scratch.fmt_ids);
        fmt_ids.clear();
        fmt_ids.extend(line.format_fields().iter().map(|f| f.id()));
        for &id in &fmt_ids {
            let tag = hdr.id_name(id as u32).unwrap_or_default();
            if Some(id) == gt_id {
                self.split_format_genotype(hdr, out_hdr, line, &mut out);
                continue;
            }
            match hdr.id_type(HeaderRecordType::Format, id as u32) {
                Some(ValueType::Integer) => {
                    self.split_format_numeric::<i32>(hdr, out_hdr, line, id, tag, &mut out)?;
                }
                Some(ValueType::Float) => {
                    self.split_format_numeric::<f32>(hdr, out_hdr, line, id, tag, &mut out)?;
                }
                _ => self.split_format_string(hdr, out_hdr, line, id, tag, &mut out)?,
            }
        }
        self.scratch.fmt_ids = fmt_ids;
        for dst in &mut out {
            dst.sync();
        }
        self.split_lines = out;
        Ok(())
    }

    /// bcftools' `split_info_numeric`, for the records of all ALTs (`dsts[i]` is ALT `i + 1`).
    fn split_info_numeric<T: Value>(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        key: i32,
        tag: &[u8],
        dsts: &mut [Record],
    ) -> anyhow::Result<()> {
        let mut vals = std::mem::take(&mut T::buffers(&mut self.scratch).a);
        let r = self.split_info_numeric_into(hdr, out_hdr, src, key, tag, dsts, &mut vals);
        T::buffers(&mut self.scratch).a = vals;
        r
    }

    #[expect(clippy::too_many_arguments, reason = "bcftools' split_info_numeric")]
    fn split_info_numeric_into<T: Value>(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        key: i32,
        tag: &[u8],
        dsts: &mut [Record],
        vals: &mut Vec<T>,
    ) -> anyhow::Result<()> {
        let Ok(ret) = T::get_info(src, hdr, tag, vals) else {
            return Ok(());
        };
        let n = src.allele_count();
        let chrom = src.chrom_safe(hdr);
        let pos = src.pos();
        let mismatch = |expected| Mismatch {
            kind: "INFO",
            tag,
            chrom,
            pos,
            expected,
            found: ret,
        };
        let error = |expected: usize| {
            format!(
                "Error: wrong number of fields in INFO/{} at {}:{}, expected {expected}, found {ret}. Use --force to proceed anyway.",
                String::from_utf8_lossy(tag),
                String::from_utf8_lossy(chrom),
                pos + 1
            )
        };
        match hdr.id_number(HeaderRecordType::Info, key as u32) {
            Some(Number::A) => {
                if ret != n - 1 {
                    self.wrong_number(&mismatch(n - 1), "", error(n - 1))?;
                    return Ok(());
                }
                for (ialt, dst) in dsts.iter_mut().enumerate() {
                    let _ = T::update_info(dst, out_hdr, tag, &vals[ialt..=ialt]);
                }
            }
            Some(Number::R) => {
                if ret != n {
                    self.wrong_number(&mismatch(n), "", error(n))?;
                    return Ok(());
                }
                // With --keep-sum, the values summed into REF: all but REF and this ALT's.
                let summed = if self.opts.keep_sum_ad == Some(key) {
                    let len = src.info_field_by_id(key).map_or(0, |f| f.len());
                    &vals[..len.min(vals.len())]
                } else {
                    &[]
                };
                for (ialt, dst) in dsts.iter_mut().enumerate() {
                    let mut v0 = vals[0];
                    for (j, &v) in summed.iter().enumerate().skip(1) {
                        if j != ialt + 1 && !v.is_missing() && !v.is_vector_end() {
                            v0 = v0.add(v);
                        }
                    }
                    let _ = T::update_info(dst, out_hdr, tag, &[v0, vals[ialt + 1]]);
                }
            }
            Some(Number::G) => {
                let ng = n * (n + 1) / 2;
                if ret != ng {
                    self.wrong_number(&mismatch(ng), "", error(ng))?;
                    return Ok(());
                }
                for (ialt, dst) in dsts.iter_mut().enumerate() {
                    let picked = [
                        vals[0],
                        vals[alleles2gt(0, ialt + 1)],
                        vals[alleles2gt(ialt + 1, ialt + 1)],
                    ];
                    let _ = T::update_info(dst, out_hdr, tag, &picked);
                }
            }
            _ => {
                for dst in dsts {
                    let _ = T::update_info(dst, out_hdr, tag, &vals[..ret]);
                }
            }
        }
        Ok(())
    }

    /// bcftools' `split_info_string`, for the records of all ALTs.
    fn split_info_string(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        key: i32,
        tag: &[u8],
        dsts: &mut [Record],
    ) {
        let Ok(s) = src.get_info_string(hdr, tag) else {
            return;
        };
        let number = hdr.id_number(HeaderRecordType::Info, key as u32);
        let out = &mut self.scratch.bytes;
        for (ialt, dst) in dsts.iter_mut().enumerate() {
            out.clear();
            match number {
                Some(Number::A) => match nth_field(s, 0, ialt) {
                    Some((b, e)) => out.extend_from_slice(&s[b..e]),
                    None => out.push(b'.'),
                },
                Some(Number::R) => pick_fields(s, &[0, ialt], out),
                Some(Number::G) => {
                    let i0a = alleles2gt(0, ialt + 1);
                    let iaa = alleles2gt(ialt + 1, ialt + 1);
                    pick_fields(s, &[0, i0a - 1, iaa - i0a - 1], out);
                }
                _ => out.extend_from_slice(s),
            }
            let _ = dst.update_info_string(out_hdr, tag, Some(out));
        }
    }

    /// bcftools' `split_format_genotype`, for the records of all ALTs.
    fn split_format_genotype(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        dsts: &mut [Record],
    ) {
        let Ok(ngts) = src.get_genotypes(hdr, &mut self.scratch.i.a) else {
            return;
        };
        let nsmpl = hdr.sample_count();
        if ngts == 0 || nsmpl == 0 {
            return;
        }
        let ploidy = ngts / nsmpl;
        let gt = &self.scratch.i.a;
        let out = &mut self.scratch.i.out;
        for (ialt, dst) in dsts.iter_mut().enumerate() {
            let target = ialt as i32 + 1;
            out.clear();
            out.extend_from_slice(gt);
            for sample in out.chunks_mut(ploidy) {
                for g in sample.iter_mut() {
                    let v = Genotype(*g);
                    if v.is_vector_end() {
                        break;
                    }
                    if v.is_missing() || v.allele() == 0 {
                        continue;
                    }
                    let phased = *g & 1;
                    *g = if v.allele() == target {
                        Genotype::unphased(1).0 | phased
                    } else if self.opts.ma_use_ref_allele {
                        Genotype::unphased(0).0 | phased
                    } else {
                        Genotype::MISSING.0 | phased
                    };
                }
            }
            let _ = dst.update_genotypes(out_hdr, out);
        }
    }

    /// bcftools' `split_format_numeric`, for the records of all ALTs.
    fn split_format_numeric<T: Value>(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        id: i32,
        tag: &[u8],
        dsts: &mut [Record],
    ) -> anyhow::Result<()> {
        let mut vals = std::mem::take(&mut T::buffers(&mut self.scratch).a);
        let mut out = std::mem::take(&mut T::buffers(&mut self.scratch).out);
        let mut haploid = std::mem::take(&mut self.scratch.haploid);
        let r = self.split_format_numeric_into(
            hdr,
            out_hdr,
            src,
            id,
            tag,
            dsts,
            &mut vals,
            &mut out,
            &mut haploid,
        );
        T::buffers(&mut self.scratch).a = vals;
        T::buffers(&mut self.scratch).out = out;
        self.scratch.haploid = haploid;
        r
    }

    #[expect(clippy::too_many_arguments, reason = "bcftools' split_format_numeric")]
    fn split_format_numeric_into<T: Value>(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        id: i32,
        tag: &[u8],
        dsts: &mut [Record],
        vals: &mut Vec<T>,
        out: &mut Vec<T>,
        haploid: &mut Vec<bool>,
    ) -> anyhow::Result<()> {
        let Ok(nvals) = T::get_format(src, hdr, tag, vals) else {
            return Ok(());
        };
        let nsmpl = hdr.sample_count();
        if nvals == 0 || nsmpl == 0 {
            return Ok(());
        }
        if nvals == nsmpl {
            // One value per sample.
            for dst in dsts {
                let _ = T::update_format(dst, out_hdr, tag, &vals[..nsmpl]);
            }
            return Ok(());
        }
        let n = src.allele_count();
        let chrom = src.chrom_safe(hdr);
        let pos = src.pos();
        let mismatch = |expected| Mismatch {
            kind: "FMT",
            tag,
            chrom,
            pos,
            expected,
            found: nvals,
        };
        let suffix = ". Removing the field.";
        let error = |expected: usize| {
            format!(
                "Error: wrong number of fields in FMT/{} at {}:{}, expected {expected}, found {nvals}. Use --force to proceed anyway.",
                String::from_utf8_lossy(tag),
                String::from_utf8_lossy(chrom),
                pos + 1
            )
        };
        match hdr.id_number(HeaderRecordType::Format, id as u32) {
            Some(Number::A) => {
                if nvals != (n - 1) * nsmpl {
                    self.wrong_number(&mismatch((n - 1) * nsmpl), suffix, error((n - 1) * nsmpl))?;
                    return Ok(());
                }
                let per = nvals / nsmpl;
                for (ialt, dst) in dsts.iter_mut().enumerate() {
                    out.clear();
                    for s in vals.chunks(per) {
                        let v = s[ialt];
                        out.push(if v.is_missing() || v.is_vector_end() {
                            T::missing()
                        } else {
                            v
                        });
                    }
                    let _ = T::update_format(dst, out_hdr, tag, out);
                }
            }
            Some(Number::R) => {
                if nvals != n * nsmpl {
                    // bcftools' warning gives the Number=A count.
                    self.wrong_number(&mismatch((n - 1) * nsmpl), suffix, error(n * nsmpl))?;
                    return Ok(());
                }
                let per = nvals / nsmpl;
                let keep_sum = self.opts.keep_sum_ad == Some(id);
                for (ialt, dst) in dsts.iter_mut().enumerate() {
                    out.clear();
                    for s in vals.chunks(per) {
                        let mut v0 = s[0];
                        if keep_sum {
                            for (j, &v) in s.iter().enumerate().skip(1) {
                                if j != ialt + 1 && !v.is_missing() && !v.is_vector_end() {
                                    v0 = v0.add(v);
                                }
                            }
                        }
                        let v1 = s[ialt + 1];
                        out.push(v0);
                        out.push(if v1.is_vector_end() { T::missing() } else { v1 });
                    }
                    let _ = T::update_format(dst, out_hdr, tag, out);
                }
            }
            Some(Number::G) => {
                let ng = n * (n + 1) / 2;
                let per = nvals / nsmpl;
                // A ploidy above two in some sample, which bcftools takes for the wrong number
                // of values: each ALT gets the genotypes of REF and that ALT.
                let mut ploidy = Vec::new();
                if per > ng
                    && ploidy_of(n, per).is_some()
                    && sample_ploidies(vals, per, n, &mut ploidy)
                {
                    let width = ploidy.iter().max().map_or(1, |p| p + 1);
                    for (ialt, dst) in dsts.iter_mut().enumerate() {
                        out.clear();
                        for (s, &p) in vals.chunks(per).zip(&ploidy) {
                            // p + 1 values: a single '.' (p is 0), or each count of the ALT.
                            if p == 0 {
                                out.push(s[0]);
                            } else {
                                out.extend((0..=p).map(|k| s[ref_alt_index(ialt + 1, p, k)]));
                            }
                            out.extend(std::iter::repeat_n(T::vector_end(), width - p - 1));
                        }
                        let _ = T::update_format(dst, out_hdr, tag, out);
                    }
                    return Ok(());
                }
                if nvals != ng * nsmpl && nvals != n * nsmpl {
                    self.wrong_number(
                        &mismatch((n - 1) * nsmpl),
                        suffix,
                        format!(
                            "Error at {}:{}, the tag {} has wrong number of fields. Use --force to proceed anyway.",
                            String::from_utf8_lossy(chrom),
                            pos + 1,
                            String::from_utf8_lossy(tag)
                        ),
                    )?;
                    return Ok(());
                }
                let per = nvals / nsmpl;
                let all_haploid = per == n;
                // A sample is haploid when its diploid-sized vector ends early.
                haploid.clear();
                if !all_haploid {
                    haploid.extend(
                        vals.chunks(per)
                            .map(|s| s.iter().any(|v| v.is_vector_end())),
                    );
                }
                for (ialt, dst) in dsts.iter_mut().enumerate() {
                    let (i0a, iaa) = (alleles2gt(0, ialt + 1), alleles2gt(ialt + 1, ialt + 1));
                    out.clear();
                    for (k, s) in vals.chunks(per).enumerate() {
                        out.push(s[0]);
                        if all_haploid {
                            out.push(s[ialt + 1]);
                        } else if haploid[k] {
                            out.push(s[ialt + 1]);
                            out.push(T::vector_end());
                        } else {
                            out.push(s[i0a]);
                            out.push(s[iaa]);
                        }
                    }
                    let _ = T::update_format(dst, out_hdr, tag, out);
                }
            }
            _ => {
                for dst in dsts {
                    let _ = T::update_format(dst, out_hdr, tag, &vals[..nvals]);
                }
            }
        }
        Ok(())
    }

    /// bcftools' `split_format_string`, for the records of all ALTs.
    fn split_format_string(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        id: i32,
        tag: &[u8],
        dsts: &mut [Record],
    ) -> anyhow::Result<()> {
        let mut bufs = std::mem::take(&mut self.scratch.format_strings);
        let r = self.split_format_string_into(hdr, out_hdr, src, id, tag, dsts, &mut bufs);
        self.scratch.format_strings = bufs;
        r
    }

    #[expect(clippy::too_many_arguments, reason = "bcftools' split_format_string")]
    fn split_format_string_into(
        &mut self,
        hdr: &Header,
        out_hdr: &Header,
        src: &mut Record,
        id: i32,
        tag: &[u8],
        dsts: &mut [Record],
        bufs: &mut FormatStrings,
    ) -> anyhow::Result<()> {
        let n = src.allele_count();
        let chrom = src.chrom_safe(hdr);
        let pos = src.pos();
        let nsmpl = hdr.sample_count();
        let number = hdr.id_number(HeaderRecordType::Format, id as u32);
        let Ok((data, blen)) = src.get_format_chars(hdr, tag) else {
            return Ok(());
        };
        // A record with fewer samples than the header fails as the numeric reads do.
        if data.is_empty() || nsmpl == 0 || blen == 0 || data.len() < nsmpl * blen {
            return Ok(());
        }
        if !matches!(number, Some(Number::A | Number::R | Number::G)) {
            for dst in dsts {
                let _ = dst.update_format_chars(out_hdr, tag, data);
            }
            return Ok(());
        }
        let samples = || data.chunks_exact(blen).take(nsmpl);
        let FormatStrings {
            starts,
            samples: bounds,
            picked,
            ends,
            buf,
        } = bufs;
        starts.clear();
        bounds.clear();
        if number == Some(Number::G) {
            // Index the fields of each sample once: whether it is haploid or diploid does not
            // depend on the ALT.
            let ng = n * (n + 1) / 2;
            bounds.push(0);
            for sample in samples() {
                let end = sample.iter().position(|&c| c == 0).unwrap_or(blen);
                if end == 1 && sample[0] == b'.' {
                    // A missing value: left as it is, with no fields indexed.
                    bounds.push(starts.len());
                    continue;
                }
                let first = starts.len();
                starts.push(0);
                starts.extend(
                    sample[..end]
                        .iter()
                        .enumerate()
                        .filter(|&(_, &c)| c == b',')
                        .map(|(k, _)| offset(k + 1)),
                );
                starts.push(offset(end + 1));
                let nfields = starts.len() - first - 1;
                if nfields != ng && nfields != n && ploidy_of(n, nfields).is_none() {
                    if !self.opts.force {
                        bail!(
                            "Error: wrong number of fields in FMT/{} at {}:{}, expected {ng} or {n}, found {nfields}. Use --force to proceed anyway.",
                            String::from_utf8_lossy(tag),
                            String::from_utf8_lossy(chrom),
                            pos + 1,
                        );
                    }
                    if !self.force_warned {
                        eprintln!(
                            "Warning: wrong number of fields in FMT/{} at {}:{}, expected {ng} or {n}, found {nfields}. Removing the field.\n         (This warning is printed only once.)",
                            String::from_utf8_lossy(tag),
                            String::from_utf8_lossy(chrom),
                            pos + 1,
                        );
                        self.force_warned = true;
                    }
                    return Ok(());
                }
                bounds.push(starts.len());
            }
        }
        let is_missing = |s: usize| number == Some(Number::G) && bounds[s] == bounds[s + 1];
        let mut other_picks = Vec::new();
        for (ialt, dst) in dsts.iter_mut().enumerate() {
            let haploid_picks = [0, ialt + 1];
            let diploid_picks = [0, alleles2gt(0, ialt + 1), alleles2gt(ialt + 1, ialt + 1)];
            picked.clear();
            ends.clear();
            let mut maxlen = 0;
            for (s, sample) in samples().enumerate() {
                let start = picked.len();
                if is_missing(s) {
                    ends.push(start);
                    continue;
                }
                match number {
                    Some(Number::G) => {
                        let fields = &starts[bounds[s]..bounds[s + 1]];
                        let nfields = fields.len() - 1;
                        let picks: &[usize] = if nfields == n {
                            &haploid_picks
                        } else if nfields == n * (n + 1) / 2 {
                            &diploid_picks
                        } else {
                            // A ploidy above two: the genotypes of REF and this ALT.
                            let p = ploidy_of(n, nfields).expect("checked with the fields");
                            other_picks.clear();
                            other_picks.extend((0..=p).map(|k| ref_alt_index(ialt + 1, p, k)));
                            &other_picks
                        };
                        push_fields(sample, fields, picks, picked);
                    }
                    Some(Number::A) => match nth_field(sample, 0, ialt) {
                        Some((b, e)) => picked.extend_from_slice(&sample[b..e]),
                        None => picked.push(b'.'),
                    },
                    _ => pick_fields(sample, &[0, ialt], picked),
                }
                maxlen = maxlen.max(picked.len() - start);
                ends.push(picked.len());
            }
            // Each value NUL-padded to the longest one, as bcftools writes the fields back in
            // place and then shrinks them; a missing value keeps its original bytes.
            let width = maxlen.min(blen);
            buf.clear();
            let mut start = 0;
            for (s, sample) in samples().enumerate() {
                let value = if is_missing(s) {
                    &sample[..width]
                } else {
                    &picked[start..ends[s]]
                };
                start = ends[s];
                let k = value.len().min(width);
                buf.extend_from_slice(&value[..k]);
                buf.resize(buf.len() + width - k, 0);
            }
            let _ = dst.update_format_chars(out_hdr, tag, buf);
        }
        Ok(())
    }
}

/// The ploidy of each sample's `per` Number=G values of `n` alleles, 0 for none or a single
/// `.`; `false` when the number of values of a sample fits no ploidy.
fn sample_ploidies<T: Value>(vals: &[T], per: usize, n: usize, out: &mut Vec<usize>) -> bool {
    out.clear();
    for s in vals.chunks(per) {
        let k = s.iter().position(|v| v.is_vector_end()).unwrap_or(per);
        if k == 0 || (k == 1 && s[0].is_missing()) {
            out.push(0);
            continue;
        }
        let Some(p) = ploidy_of(n, k) else {
            return false;
        };
        out.push(p);
    }
    true
}

/// An offset within one sample's string, at most one past its end. The bytes per sample
/// come from the field's `i32` count, so the offset always fits.
fn offset(k: usize) -> u32 {
    u32::try_from(k).expect("a sample's FORMAT string is shorter than 2^31 bytes")
}

/// Reusable buffers for splitting string FORMAT fields.
#[derive(Debug, Default)]
pub(super) struct FormatStrings {
    /// Number=G: where each field of each sample starts, then one past the sample's end
    /// (32 bits each, as there are as many as fields of all samples).
    starts: Vec<u32>,
    /// Number=G: where each sample's entries in `starts` begin, then where the last ends. A
    /// sample with no entries is missing.
    samples: Vec<usize>,
    /// The values picked for one ALT, concatenated.
    picked: Vec<u8>,
    /// Where each sample's value in `picked` ends.
    ends: Vec<usize>,
    /// The NUL-padded values of one ALT.
    buf: Vec<u8>,
}

/// Appends fields `picks` (ascending) of `sample`, joined by commas, given where its fields
/// start (`starts`, then one past its end). Gives `.` when a picked field is the last one and
/// empty, as [`pick_fields`] does.
fn push_fields(sample: &[u8], starts: &[u32], picks: &[usize], out: &mut Vec<u8>) {
    let base = out.len();
    let last = starts.len() - 2;
    for (k, &f) in picks.iter().enumerate() {
        let (b, e) = (starts[f] as usize, starts[f + 1] as usize - 1);
        if f == last && b == e {
            out.truncate(base);
            out.push(b'.');
            return;
        }
        if k > 0 {
            out.push(b',');
        }
        out.extend_from_slice(&sample[b..e]);
    }
}

/// Appends comma-separated fields from `s` to `out` as bcftools' chained `string_move_nth`
/// calls pick them: each index counts from the field after the previous one picked. Any
/// missing field makes the result `.`.
fn pick_fields(s: &[u8], rel: &[usize], out: &mut Vec<u8>) {
    let base = out.len();
    let mut start = 0;
    for (k, &nth) in rel.iter().enumerate() {
        if k > 0 {
            out.push(b',');
        }
        match nth_field(s, start, nth) {
            Some((b, e)) => {
                out.extend_from_slice(&s[b..e]);
                start = e + 1;
            }
            None => {
                out.truncate(base);
                out.push(b'.');
                return;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pick_fields_counts_from_previous() {
        let mut out = b"x;".to_vec();
        pick_fields(b"a,b,c,d,e,f", &[0, 1], &mut out);
        assert_eq!(out, b"x;a,c");
        out.clear();
        pick_fields(b"a,b,c,d,e,f", &[0, 0, 1], &mut out);
        assert_eq!(out, b"a,b,d");
        out.truncate(1);
        pick_fields(b"a,b", &[0, 3], &mut out);
        assert_eq!(out, b"a.");
    }

    #[test]
    fn push_fields_matches_pick_fields() {
        // Absolute field indexes against pick_fields' relative ones, including empty fields
        // in the middle and an empty last field.
        for sample in [&b"a,bb,,c,d,e\0\0"[..], b"a,b,c,,e,", b"0,1,2,3,4,5"] {
            let end = sample.iter().position(|&c| c == 0).unwrap_or(sample.len());
            let mut starts = vec![0];
            starts.extend(
                (0..end)
                    .filter(|&k| sample[k] == b',')
                    .map(|k| offset(k + 1)),
            );
            starts.push(offset(end + 1));
            for picks in [[0, 1, 2], [0, 2, 5], [0, 3, 5], [0, 1, 4]] {
                let rel = [picks[0], picks[1] - picks[0] - 1, picks[2] - picks[1] - 1];
                let (mut want, mut got) = (b"p".to_vec(), b"p".to_vec());
                pick_fields(sample, &rel, &mut want);
                push_fields(sample, &starts, &picks, &mut got);
                assert_eq!(got, want, "sample {sample:?} picks {picks:?}");
            }
        }
    }
}
