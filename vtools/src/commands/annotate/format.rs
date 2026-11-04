//! The FORMAT column setters (bcftools' `setter_format_*`, `vcf_setter_format_*`, and
//! `core_setter_format_*`).

use anyhow::bail;
use pasta::variant::{
    Header, Number, Record, Unpack,
    encoding::{
        INT32_MISSING, INT32_VECTOR_END, float_is_missing, float_is_vector_end, float_missing,
        float_vector_end,
    },
    number::{push_int, strtod, strtol},
    record::GetError,
};

use super::{
    columns::{Column, Replace},
    line::AnnotLine,
    setters::{CopiedAlleles, Ctx, Set, Src, at, escape, lossy, site},
    vcmp::{AlleleList, Vcmp},
};

/// HTSlib's GT encodings.
const GT_MISSING: i32 = 0;

fn gt_phased(a: i32) -> i32 {
    ((a + 1) << 1) | 1
}

fn gt_unphased(a: i32) -> i32 {
    (a + 1) << 1
}

fn gt_is_missing(v: i32) -> bool {
    (v >> 1) == 0
}

/// A value of a FORMAT field, integer or float, with its missing and vector-end markers.
trait FmtValue: Copy {
    const MISSING: Self;
    const VECTOR_END: Self;
    fn is_missing(self) -> bool;
    fn is_vector_end(self) -> bool;
}

impl FmtValue for i32 {
    const MISSING: Self = INT32_MISSING;
    const VECTOR_END: Self = INT32_VECTOR_END;
    fn is_missing(self) -> bool {
        self == INT32_MISSING
    }
    fn is_vector_end(self) -> bool {
        self == INT32_VECTOR_END
    }
}

impl FmtValue for f32 {
    const MISSING: Self = f32::from_bits(pasta::variant::encoding::FLOAT_MISSING);
    const VECTOR_END: Self = f32::from_bits(pasta::variant::encoding::FLOAT_VECTOR_END);
    fn is_missing(self) -> bool {
        float_is_missing(self)
    }
    fn is_vector_end(self) -> bool {
        float_is_vector_end(self)
    }
}

/// Merges `nvals` values per annotated sample into the record's values `ori` (`ndst` per
/// sample, `ndst <= 0` when absent), following the replace mode and the sample map (the body
/// shared by bcftools' `core_setter_format_int` and `core_setter_format_real`). Returns the
/// values to set and their number per sample, or `None` to leave the record alone.
fn merge_samples<T: FmtValue>(
    replace: Replace,
    map: &[i32],
    vals: &[T],
    nvals: usize,
    ori: &mut Vec<T>,
    ndst: i64,
    out: &mut Vec<T>,
) -> Option<usize> {
    let nsmpl = map.len();
    if ndst <= 0 {
        if replace.has(Replace::NON_MISSING) {
            return None;
        }
        out.clear();
        for &isrc in map {
            if isrc == -1 {
                out.push(T::MISSING);
                out.extend(std::iter::repeat_n(T::VECTOR_END, nvals.saturating_sub(1)));
            } else {
                let s = isrc as usize * nvals;
                out.extend_from_slice(&vals[s..s + nvals]);
            }
        }
        return Some(nvals);
    }
    let ndst = ndst as usize;
    let use_new = |dst0: T, src0: T| {
        if replace.has(Replace::NON_MISSING) {
            !dst0.is_missing()
        } else if replace.has(Replace::MISSING) {
            dst0.is_missing()
        } else if replace.has(Replace::ALL) {
            !src0.is_missing()
        } else {
            true
        }
    };
    if ndst >= nvals {
        for (i, &isrc) in map.iter().enumerate() {
            if isrc == -1 {
                continue;
            }
            let src = &vals[isrc as usize * nvals..][..nvals];
            let dst = &mut ori[ndst * i..][..ndst];
            if !use_new(dst[0], src[0]) {
                continue;
            }
            dst[..nvals].copy_from_slice(src);
            dst[nvals..].fill(T::VECTOR_END);
        }
        std::mem::swap(out, ori);
        return Some(ndst);
    }
    out.clear();
    for (i, &isrc) in map.iter().enumerate().take(nsmpl) {
        let o = &ori[ndst * i..][..ndst];
        let take = isrc != -1 && use_new(o[0], vals[isrc as usize * nvals]);
        if take {
            out.extend_from_slice(&vals[isrc as usize * nvals..][..nvals]);
        } else {
            out.extend_from_slice(o);
            out.extend(std::iter::repeat_n(T::VECTOR_END, nvals - ndst));
        }
    }
    Some(nvals)
}

/// Copies one sample's per-allele or per-genotype values `src` into `dst` along `map` (from
/// the source index to the destination index, negative when dropped). A source shorter than
/// the map, as a single `.` is, gives missing values: bcftools reads past the sample's values.
fn remap_values<T: FmtValue>(map: &[i32], src: &[T], dst: &mut [T]) {
    for (j, &k) in map.iter().enumerate() {
        if let Ok(k) = usize::try_from(k)
            && let Some(d) = dst.get_mut(k)
        {
            *d = src.get(j).copied().unwrap_or(T::MISSING);
        }
    }
}

/// Rewrites the allele indexes of genotypes along `map` (from the source allele to the
/// record's, negative when the record lacks it), so that they name the same alleles. An allele
/// the record lacks, or an index past the map, becomes missing; phasing and vector ends stay.
fn remap_genotypes(map: &[i32], gts: &mut [i32]) {
    for g in gts {
        if *g == INT32_VECTOR_END || gt_is_missing(*g) {
            continue;
        }
        let ial = usize::try_from((*g >> 1) - 1).ok();
        let to = ial.and_then(|i| map.get(i)).copied().unwrap_or(-1);
        *g = (*g & 1) | if to >= 0 { (to + 1) << 1 } else { GT_MISSING };
    }
}

/// The alleles of a record once `REF` is copied from an annotation record: that REF, then the
/// record's ALTs.
struct WithRef<'a> {
    ref_: &'a [u8],
    line: &'a Record,
}

impl AlleleList for WithRef<'_> {
    fn count(&self) -> usize {
        self.line.allele_count()
    }

    fn allele_at(&self, i: usize) -> &[u8] {
        if i == 0 {
            self.ref_
        } else {
            self.line.allele(i)
        }
    }
}

/// Whether the alleles of `line` start with those of `rec` (compared undecoded unless both are
/// decoded), when the genotypes of `rec` name the same alleles in `line`.
#[inline]
fn starts_with_alleles(line: &Record, rec: &Record) -> bool {
    match (line.decoded_alleles(), rec.decoded_alleles()) {
        (Some(a), Some(b)) => a.starts_with(b),
        _ => {
            rec.allele_count() <= line.allele_count()
                && rec
                    .peek_alleles()
                    .zip(line.peek_alleles())
                    .all(|(a, b)| a == b)
        }
    }
}

/// Rewrites the genotypes `gts` copied from the annotation record `rec` to the alleles the
/// record `line` has once every column is set, through the allele map of Number=A/R values.
/// They stay when `line` takes the ALTs of `rec` or already starts with its alleles; REF alleles
/// that are not compatible map no allele.
// Out of line: inlined, it made the GT setter slower on every record (2-3% of annotate).
#[inline(never)]
fn remap_to_record_alleles(
    vcmp: &mut Vcmp,
    copied: CopiedAlleles,
    line: &mut Record,
    rec: &mut Record,
    gts: &mut [i32],
) {
    match copied {
        // With the REF of `rec` too, the alleles are then the same.
        CopiedAlleles::None | CopiedAlleles::Ref if starts_with_alleles(line, rec) => return,
        CopiedAlleles::Alts => return,
        _ => {}
    }
    let n = rec.allele_count();
    rec.unpack(Unpack::STR);
    line.unpack(Unpack::STR);
    let map = if copied == CopiedAlleles::Ref {
        let als = WithRef {
            ref_: rec.allele(0),
            line,
        };
        vcmp.map_ar_values(n, &als, &*rec)
    } else {
        vcmp.map_ar_values(n, &*line, &*rec)
    };
    match map {
        Some(map) if map.iter().enumerate().all(|(i, &m)| m == i as i32) => {}
        Some(map) => remap_genotypes(map, gts),
        None => remap_genotypes(&[], gts),
    }
}

/// The number of values per sample in a tab line's sample columns (bcftools' `count_vals`).
fn count_vals(tab: &AnnotLine, beg: i32, end: i32) -> usize {
    let mut nmax = 1;
    for i in beg..end {
        let s = tab.col(i);
        if s == b"." {
            continue;
        }
        let n = 1 + s.iter().filter(|&&c| c == b',').count();
        nmax = nmax.max(n);
    }
    nmax
}

/// The ploidy of each sample from the number of Number=G values (bcftools'
/// `determine_ploidy_*`): 1, 2, or 0 when all missing. Returns the largest, or `Err` with the
/// unexpected count.
fn determine_ploidy<T: FmtValue>(
    nals: usize,
    vals: &[T],
    nvals1: usize,
    out: &mut Vec<u8>,
    nsmpl: usize,
) -> Result<u8, usize> {
    let ndip = nals * (nals + 1) / 2;
    let mut max = 0;
    out.clear();
    for i in 0..nsmpl {
        let ptr = &vals[i * nvals1..][..nvals1];
        let mut has_value = false;
        let mut j = 0;
        while j < nvals1 {
            if ptr[j].is_vector_end() {
                break;
            }
            if !ptr[j].is_missing() {
                has_value = true;
            }
            j += 1;
        }
        if has_value {
            if j == ndip {
                out.push(2);
                max = 2;
            } else if j == nals {
                out.push(1);
                if max == 0 {
                    max = 1;
                }
            } else {
                return Err(j);
            }
        } else {
            out.push(0);
        }
    }
    Ok(max)
}

/// Whether a getter's status means the tag is absent (`-3`), which leaves the record alone,
/// rather than another error, which bcftools reports as "not done".
fn absent_or_error(e: &GetError) -> Set {
    if matches!(e, GetError::Absent) {
        Set::Done
    } else {
        Set::More
    }
}

impl Ctx {
    /// bcftools' `core_setter_format_int`.
    fn core_format_int(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        nvals: usize,
    ) -> Set {
        let Some(map) = self.smpl.map.as_deref() else {
            let n = nvals * self.smpl.nsmpl_annot;
            let r = line.update_format_i32(&self.hdr_out, &col.key_dst, &self.s.tmpi[..n]);
            return if r.is_ok() { Set::Done } else { Set::Failed };
        };
        let nsmpl = self.hdr_out.sample_count();
        let ndst = match line.get_format_i32(hdr, &col.key_dst, &mut self.s.tmpi2) {
            Ok(n) if n > 0 => (n / nsmpl) as i64,
            _ => 0,
        };
        let s = &mut self.s;
        let Some(n1) = merge_samples(
            col.replace,
            map,
            &s.tmpi,
            nvals,
            &mut s.tmpi2,
            ndst,
            &mut s.tmpi3,
        ) else {
            return Set::Done;
        };
        let r = line.update_format_i32(&self.hdr_out, &col.key_dst, &s.tmpi3[..n1 * nsmpl]);
        if r.is_ok() { Set::Done } else { Set::Failed }
    }

    /// bcftools' `core_setter_format_real`.
    fn core_format_real(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        nvals: usize,
    ) -> Set {
        let Some(map) = self.smpl.map.as_deref() else {
            let n = nvals * self.smpl.nsmpl_annot;
            let r = line.update_format_f32(&self.hdr_out, &col.key_dst, &self.s.tmpf[..n]);
            return if r.is_ok() { Set::Done } else { Set::Failed };
        };
        let nsmpl = self.hdr_out.sample_count();
        let ndst = match line.get_format_f32(hdr, &col.key_dst, &mut self.s.tmpf2) {
            Ok(n) if n > 0 => (n / nsmpl) as i64,
            _ => 0,
        };
        let s = &mut self.s;
        let Some(n1) = merge_samples(
            col.replace,
            map,
            &s.tmpf,
            nvals,
            &mut s.tmpf2,
            ndst,
            &mut s.tmpf3,
        ) else {
            return Set::Done;
        };
        let r = line.update_format_f32(&self.hdr_out, &col.key_dst, &s.tmpf3[..n1 * nsmpl]);
        if r.is_ok() { Set::Done } else { Set::Failed }
    }

    /// bcftools' `core_setter_format_str`, with the annotated samples' strings in `strs`.
    fn core_format_str(&mut self, col: &Column, hdr: &Header, line: &mut Record) -> Set {
        let s = &mut self.s;
        let Some(map) = self.smpl.map.as_deref() else {
            s.strs.pack(&mut s.tmpks);
            let r = line.update_format_chars(&self.hdr_out, &col.key_dst, &s.tmpks);
            return if r.is_ok() { Set::Done } else { Set::Failed };
        };
        let nsmpl = self.hdr_out.sample_count();
        s.strs2.clear();
        match line.get_format_chars(hdr, &col.key_dst) {
            Ok((data, n)) if n > 0 && !data.is_empty() => {
                for c in data.chunks(n) {
                    s.strs2
                        .push(&c[..c.iter().position(|&b| b == 0).unwrap_or(c.len())]);
                }
            }
            _ => {
                for _ in 0..nsmpl {
                    s.strs2.push(b".");
                }
            }
        }
        for (i, &isrc) in map.iter().enumerate() {
            if isrc == -1 {
                continue;
            }
            let src = s.strs.get(isrc as usize);
            let dst = s.strs2.get(i);
            let skip = if col.replace.has(Replace::NON_MISSING) {
                dst == b"."
            } else if col.replace.has(Replace::MISSING) {
                dst != b"."
            } else if col.replace.has(Replace::ALL) {
                src == b"."
            } else {
                false
            };
            if !skip {
                s.strs2.set(i, src);
            }
        }
        s.strs2.pack(&mut s.tmpks);
        let r = line.update_format_chars(&self.hdr_out, &col.key_dst, &s.tmpks);
        if r.is_ok() { Set::Done } else { Set::Failed }
    }

    /// Checks that a tab line has a column for each annotated sample.
    fn check_sample_columns(
        &self,
        col: &Column,
        hdr: &Header,
        line: &Record,
        tab: &AnnotLine,
    ) -> anyhow::Result<()> {
        if col.icol as usize + self.smpl.nsmpl_annot > tab.ncols() {
            bail!(
                "Incorrect number of values for {} at {}",
                lossy(&col.key_src),
                site(hdr, line)
            );
        }
        Ok(())
    }

    /// bcftools' `setter_format_int`.
    pub(super) fn set_format_int(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with FORMAT tags (yet?)");
        };
        self.check_sample_columns(col, hdr, line, tab)?;
        let n = self.smpl.nsmpl_annot;
        let nvals = count_vals(tab, col.icol, col.icol + n as i32);
        let tmpi = &mut self.s.tmpi;
        tmpi.clear();
        for ismpl in 0..n {
            let s = tab.col(col.icol + ismpl as i32);
            let start = tmpi.len();
            let mut p = 0;
            while at(s, p) != 0 {
                if at(s, p) == b'.' && (at(s, p + 1) == 0 || at(s, p + 1) == b',') {
                    tmpi.push(INT32_MISSING);
                    p += if at(s, p + 1) != 0 { 2 } else { 1 };
                    continue;
                }
                let (v, end) = strtol(s, p);
                if end == p || (at(s, end) != 0 && at(s, end) != b',') {
                    bail!(
                        "Could not parse {} (Type=Integer) at {} .. [{}]",
                        lossy(&col.key_src),
                        site(hdr, line),
                        lossy(tab.col(col.icol))
                    );
                }
                tmpi.push(v as i32);
                p = if at(s, end) != 0 { end + 1 } else { end };
            }
            tmpi.resize(start + nvals, INT32_VECTOR_END);
        }
        Ok(self.core_format_int(col, hdr, line, nvals))
    }

    /// bcftools' `setter_format_real`.
    pub(super) fn set_format_real(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with FORMAT tags (yet?)");
        };
        self.check_sample_columns(col, hdr, line, tab)?;
        let n = self.smpl.nsmpl_annot;
        let nvals = count_vals(tab, col.icol, col.icol + n as i32);
        let tmpf = &mut self.s.tmpf;
        tmpf.clear();
        for ismpl in 0..n {
            let s = tab.col(col.icol + ismpl as i32);
            let start = tmpf.len();
            let mut p = 0;
            while at(s, p) != 0 {
                if at(s, p) == b'.' && (at(s, p + 1) == 0 || at(s, p + 1) == b',') {
                    tmpf.push(float_missing());
                    p += if at(s, p + 1) != 0 { 2 } else { 1 };
                    continue;
                }
                let (v, end) = strtod(s, p);
                if end == p {
                    bail!(
                        "Could not parse {} at {} .. [{}]",
                        lossy(&col.key_src),
                        site(hdr, line),
                        lossy(tab.col(col.icol))
                    );
                }
                tmpf.push(v as f32);
                p = if at(s, end) != 0 { end + 1 } else { end };
            }
            tmpf.resize(start + nvals, float_vector_end());
        }
        Ok(self.core_format_real(col, hdr, line, nvals))
    }

    /// bcftools' `setter_format_str`.
    pub(super) fn set_format_str(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with FORMAT tags (yet?)");
        };
        self.check_sample_columns(col, hdr, line, tab)?;
        let mut esc = std::mem::take(&mut self.s.tmps2);
        self.s.strs.clear();
        for ismpl in 0..self.smpl.nsmpl_annot {
            escape(tab.col(col.icol + ismpl as i32), b":", &mut esc);
            self.s.strs.push(&esc);
        }
        self.s.tmps2 = esc;
        Ok(self.core_format_str(col, hdr, line))
    }

    /// Converts genotypes to strings of fixed width per sample (bcftools'
    /// `genotypes_to_string`), into `tmpks`.
    fn genotypes_to_string(&mut self, nsrc1: usize, nsmpl_dst: usize) {
        let mut blen = if nsrc1 > 1 { nsrc1 + 1 } else { 1 };
        let s = &mut self.s;
        'retry: loop {
            s.tmpks.clear();
            for idst in 0..nsmpl_dst {
                let isrc = self.smpl.map.as_deref().map_or(idst as i32, |m| m[idst]);
                if isrc == -1 {
                    s.tmpks.push(b'.');
                    s.tmpks.extend(std::iter::repeat_n(0, blen - 1));
                    continue;
                }
                let plen = s.tmpks.len();
                let ptr = &s.tmpi[isrc as usize * nsrc1..][..nsrc1];
                let mut i = 0;
                while i < nsrc1 && ptr[i] != INT32_VECTOR_END {
                    if i > 0 {
                        s.tmpks.push(if ptr[i] & 1 != 0 { b'|' } else { b'/' });
                    }
                    if gt_is_missing(ptr[i]) {
                        s.tmpks.push(b'.');
                    } else {
                        push_int(&mut s.tmpks, i64::from((ptr[i] >> 1) - 1));
                    }
                    i += 1;
                }
                if i == 0 {
                    s.tmpks.push(b'.');
                }
                if s.tmpks.len() - plen > blen {
                    blen *= 2;
                    continue 'retry;
                }
                s.tmpks.resize(plen + blen, 0);
            }
            break;
        }
    }

    /// bcftools' `vcf_setter_format_gt`, with the allele indexes rewritten to the record's
    /// alleles.
    pub(super) fn vcf_set_format_gt(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> anyhow::Result<Set> {
        let nsrc = match rec.get_genotypes(ahdr, &mut self.s.tmpi) {
            Ok(0) => return Ok(Set::More),
            Ok(n) => n,
            Err(e) => return Ok(absent_or_error(&e)),
        };
        // bcftools copies the allele indexes as they are, which name other alleles when the
        // record lists them in another order or lacks some (a deliberate difference).
        remap_to_record_alleles(
            &mut self.vcmp,
            self.copied_alleles,
            line,
            rec,
            &mut self.s.tmpi[..nsrc],
        );
        let nsmpl_src = ahdr.sample_count();
        let nsmpl_out = self.hdr_out.sample_count();
        if col.key_src != col.key_dst {
            self.genotypes_to_string(nsrc / nsmpl_src, nsmpl_out);
            let r = line.update_format_chars(&self.hdr_out, &col.key_dst, &self.s.tmpks);
            return Ok(if r.is_ok() { Set::Done } else { Set::Failed });
        }
        let Some(map) = self.smpl.map.as_deref() else {
            let r = line.update_genotypes(&self.hdr_out, &self.s.tmpi[..nsrc]);
            return Ok(if r.is_ok() { Set::Done } else { Set::Failed });
        };
        let s = &mut self.s;
        let ndst = match line.get_genotypes(hdr, &mut s.tmpi2) {
            Ok(n) if n > 0 => n / nsmpl_out,
            _ => 0,
        };
        let nsrc1 = nsrc / nsmpl_src;
        let replace = col.replace;
        let out: &[i32] = if ndst == 0 {
            if replace.has(Replace::NON_MISSING) {
                return Ok(Set::Done);
            }
            s.tmpi2.clear();
            for &isrc in map {
                if isrc == -1 {
                    s.tmpi2.push(GT_MISSING);
                    s.tmpi2.extend(std::iter::repeat_n(
                        INT32_VECTOR_END,
                        nsrc1.saturating_sub(1),
                    ));
                } else {
                    s.tmpi2
                        .extend_from_slice(&s.tmpi[isrc as usize * nsrc1..][..nsrc1]);
                }
            }
            &s.tmpi2
        } else if ndst >= nsrc1 {
            for (i, &isrc) in map.iter().enumerate() {
                if isrc == -1 {
                    continue;
                }
                let dst = &mut s.tmpi2[ndst * i..][..ndst];
                if replace.has(Replace::NON_MISSING) && gt_is_missing(dst[0]) {
                    continue;
                }
                if replace.has(Replace::MISSING) && !gt_is_missing(dst[0]) {
                    continue;
                }
                dst[..nsrc1].copy_from_slice(&s.tmpi[isrc as usize * nsrc1..][..nsrc1]);
                dst[nsrc1..].fill(INT32_VECTOR_END);
            }
            &s.tmpi2
        } else {
            s.tmpi3.clear();
            for (i, &isrc) in map.iter().enumerate() {
                let ori = &s.tmpi2[ndst * i..][..ndst];
                let keep_ori = isrc == -1
                    || (replace.has(Replace::NON_MISSING) && gt_is_missing(ori[0]))
                    || (replace.has(Replace::MISSING) && !gt_is_missing(ori[0]));
                if keep_ori {
                    s.tmpi3.extend_from_slice(ori);
                    s.tmpi3
                        .extend(std::iter::repeat_n(INT32_VECTOR_END, nsrc1 - ndst));
                } else {
                    s.tmpi3
                        .extend_from_slice(&s.tmpi[isrc as usize * nsrc1..][..nsrc1]);
                }
            }
            &s.tmpi3
        };
        let r = line.update_genotypes(&self.hdr_out, out);
        Ok(if r.is_ok() { Set::Done } else { Set::Failed })
    }

    /// The allele-aware part of bcftools' `vcf_setter_format_int` and `_real`: maps the
    /// source values (`src`, `nsrc1` per source sample) of a Number=A/R/G tag onto the
    /// record's alleles along the allele map the caller made with
    /// [`map_ar_values`](super::vcmp::Vcmp::map_ar_values). Returns the values and their
    /// number per sample.
    #[expect(clippy::too_many_arguments, reason = "mirrors the C setter's state")]
    fn map_format_alleles<T: FmtValue>(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &Record,
        rec: &Record,
        ahdr: &Header,
        src: &[T],
        nsrc1: usize,
        dst: &mut Vec<T>,
        ndst_total: Option<usize>,
    ) -> anyhow::Result<Option<usize>> {
        let nals = line.allele_count();
        let nsmpl_dst = line.sample_count();
        let nsmpl_src = ahdr.sample_count();
        let is_g = col.number == Some(Number::G);
        let mut ndst1 = ndst_total.map_or(0, |n| n / nsmpl_dst.max(1));
        if ndst_total.is_none() {
            if col.replace.has(Replace::NON_MISSING) {
                return Ok(None);
            }
            ndst1 = if is_g {
                nals * (nals + 1) / 2
            } else if col.number == Some(Number::A) {
                nals - 1
            } else {
                nals
            };
            dst.clear();
            dst.resize(ndst1 * nsmpl_dst, T::MISSING);
        }
        let mut map_dip_len = 0;
        if is_g {
            map_dip_len = self.vcmp.map_dip_g_values().len();
            let pld_src = determine_ploidy(
                rec.allele_count(),
                src,
                nsrc1,
                &mut self.s.src_pld,
                nsmpl_src,
            );
            if let Err(j) = pld_src {
                bail!(
                    "Unexpected number of {} values ({j}) for {} alleles at {}",
                    lossy(&col.key_src),
                    rec.allele_count(),
                    site(ahdr, rec)
                );
            }
            let pld_dst = match determine_ploidy(nals, dst, ndst1, &mut self.s.dst_pld, nsmpl_dst) {
                Ok(p) => p,
                Err(j) => bail!(
                    "Unexpected number of {} values ({j}) for {} alleles at {}",
                    lossy(&col.key_src),
                    nals,
                    site(hdr, line)
                ),
            };
            let new = if pld_dst == 1 {
                nals
            } else {
                nals * (nals + 1) / 2
            };
            if new != ndst1 {
                if ndst1 != 0 {
                    bail!(
                        "todo: {} ndst1!=ndst .. {new} {ndst1}  at {}",
                        lossy(&col.key_src),
                        site(hdr, line)
                    );
                }
                ndst1 = new;
                dst.resize(ndst1 * nsmpl_dst, T::MISSING);
            }
        } else if ndst1 == 0 {
            ndst1 = if col.number == Some(Number::A) {
                nals - 1
            } else {
                nals
            };
            dst.resize(ndst1 * nsmpl_dst, T::MISSING);
        }
        let map_hap = self.vcmp.map();
        let map_dip = &self.vcmp.map_dip()[..map_dip_len];
        let smap = self.smpl.map.as_deref();
        for i in 0..nsmpl_dst {
            // bcftools indexes the source by the destination sample and vice versa; the
            // sample map is from destination to source.
            let isrc = smap.map_or(i as i32, |m| m.get(i).copied().unwrap_or(-1));
            if isrc < 0 || isrc as usize >= nsmpl_src {
                continue;
            }
            let isrc = isrc as usize;
            let ptr_src = &src[isrc * nsrc1..][..nsrc1];
            let ptr_dst = &mut dst[i * ndst1..][..ndst1];
            if is_g {
                let (sp, dp) = (self.s.src_pld[isrc], self.s.dst_pld[i]);
                if sp > 0 && dp > 0 && sp != dp {
                    bail!("Sample ploidy differs at {}", site(hdr, line));
                }
                if dp == 0 {
                    ptr_dst.fill(T::MISSING);
                }
            }
            if !is_g || self.s.src_pld[isrc] == 1 {
                remap_values(map_hap, ptr_src, ptr_dst);
                if is_g {
                    for v in ptr_dst.iter_mut().skip(nals) {
                        *v = T::VECTOR_END;
                    }
                }
            } else {
                remap_values(map_dip, ptr_src, ptr_dst);
            }
        }
        Ok(Some(ndst1))
    }

    /// Whether the source alleles map onto the record's unchanged.
    fn alleles_unchanged(&self, rec: &Record, line: &Record) -> bool {
        rec.allele_count() == line.allele_count()
            && self
                .vcmp
                .map()
                .iter()
                .enumerate()
                .all(|(i, &m)| m == i as i32)
    }

    /// bcftools' `vcf_setter_format_int`.
    pub(super) fn vcf_set_format_int(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> anyhow::Result<Set> {
        let nsrc = match rec.get_format_i32(ahdr, &col.key_src, &mut self.s.tmpi) {
            Ok(0) => return Ok(Set::More),
            Ok(n) => n,
            Err(e) => return Ok(absent_or_error(&e)),
        };
        let nsrc1 = nsrc / ahdr.sample_count();
        if !matches!(col.number, Some(Number::G | Number::R | Number::A)) {
            return Ok(self.core_format_int(col, hdr, line, nsrc1));
        }
        rec.unpack(Unpack::STR);
        line.unpack(Unpack::STR);
        let nmap_hap = if matches!(col.number, Some(Number::G | Number::R)) {
            rec.allele_count()
        } else {
            rec.allele_count() - 1
        };
        if self.vcmp.map_ar_values(nmap_hap, &*line, &*rec).is_none() {
            bail!("REF alleles not compatible at {}", site(hdr, line));
        }
        if self.alleles_unchanged(rec, line) {
            return Ok(self.core_format_int(col, hdr, line, nsrc1));
        }
        let ndst = match line.get_format_i32(hdr, &col.key_dst, &mut self.s.tmpi2) {
            Ok(n) if n > 0 => Some(n),
            _ => None,
        };
        let src = std::mem::take(&mut self.s.tmpi);
        let mut dst = std::mem::take(&mut self.s.tmpi2);
        let r = self.map_format_alleles(col, hdr, line, rec, ahdr, &src, nsrc1, &mut dst, ndst);
        self.s.tmpi = src;
        let res = match r {
            Ok(Some(ndst1)) => {
                let n = ndst1 * line.sample_count();
                let r = line.update_format_i32(&self.hdr_out, &col.key_dst, &dst[..n]);
                Ok(if r.is_ok() { Set::Done } else { Set::Failed })
            }
            Ok(None) => Ok(Set::Done),
            Err(e) => Err(e),
        };
        self.s.tmpi2 = dst;
        res
    }

    /// bcftools' `vcf_setter_format_real`.
    pub(super) fn vcf_set_format_real(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> anyhow::Result<Set> {
        let nsrc = match rec.get_format_f32(ahdr, &col.key_src, &mut self.s.tmpf) {
            Ok(0) => return Ok(Set::More),
            Ok(n) => n,
            Err(e) => return Ok(absent_or_error(&e)),
        };
        let nsrc1 = nsrc / ahdr.sample_count();
        if !matches!(col.number, Some(Number::G | Number::R | Number::A)) {
            return Ok(self.core_format_real(col, hdr, line, nsrc1));
        }
        rec.unpack(Unpack::STR);
        line.unpack(Unpack::STR);
        let nmap_hap = if matches!(col.number, Some(Number::G | Number::R)) {
            rec.allele_count()
        } else {
            rec.allele_count() - 1
        };
        if self.vcmp.map_ar_values(nmap_hap, &*line, &*rec).is_none() {
            bail!("REF alleles not compatible at {}", site(hdr, line));
        }
        if self.alleles_unchanged(rec, line) {
            return Ok(self.core_format_real(col, hdr, line, nsrc1));
        }
        let ndst = match line.get_format_f32(hdr, &col.key_dst, &mut self.s.tmpf2) {
            Ok(n) if n > 0 => Some(n),
            _ => None,
        };
        let src = std::mem::take(&mut self.s.tmpf);
        let mut dst = std::mem::take(&mut self.s.tmpf2);
        let r = self.map_format_alleles(col, hdr, line, rec, ahdr, &src, nsrc1, &mut dst, ndst);
        self.s.tmpf = src;
        let res = match r {
            Ok(Some(ndst1)) => {
                let n = ndst1 * line.sample_count();
                let r = line.update_format_f32(&self.hdr_out, &col.key_dst, &dst[..n]);
                Ok(if r.is_ok() { Set::Done } else { Set::Failed })
            }
            Ok(None) => Ok(Set::Done),
            Err(e) => Err(e),
        };
        self.s.tmpf2 = dst;
        res
    }

    /// bcftools' `vcf_setter_format_str`, including `FMT/GT:=STR` (strings to genotypes).
    pub(super) fn vcf_set_format_str(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> anyhow::Result<Set> {
        let nsmpl_src = ahdr.sample_count();
        let nsrc1 = {
            let s = &mut self.s;
            let (data, n) = match rec.get_format_chars(ahdr, &col.key_src) {
                Ok((_, 0)) => return Ok(Set::More),
                Ok(v) => v,
                Err(e) => return Ok(absent_or_error(&e)),
            };
            s.strs.clear();
            s.tmps.clear();
            s.tmps.extend_from_slice(data);
            for c in data.chunks(n) {
                s.strs
                    .push(&c[..c.iter().position(|&b| b == 0).unwrap_or(c.len())]);
            }
            n
        };
        if col.key_dst != b"GT" {
            return Ok(self.core_format_str(col, hdr, line));
        }
        // FMT/GT:=STR: parse the strings as genotypes.
        let nsmpl_dst = self.hdr_out.sample_count();
        let mut ndst1 = 0;
        for c in self.s.tmps.chunks(nsrc1).take(nsmpl_src) {
            let n = 1 + c.iter().filter(|&&b| b == b'/' || b == b'|').count();
            ndst1 = ndst1.max(n);
        }
        line.unpack(Unpack::STR);
        let nals = line.allele_count() as i64;
        let s = &mut self.s;
        s.tmpi.clear();
        for idst in 0..nsmpl_dst {
            let isrc = self.smpl.map.as_deref().map_or(idst as i32, |m| m[idst]);
            let start = s.tmpi.len();
            if isrc == -1 {
                s.tmpi.push(GT_MISSING);
                s.tmpi.resize(start + ndst1, INT32_VECTOR_END);
                continue;
            }
            let gt = s.strs.get(isrc as usize);
            let mut is_phased = false;
            let mut beg = 0;
            while beg < gt.len() {
                let mut end = beg;
                while end < gt.len() && gt[end] != b'/' && gt[end] != b'|' {
                    end += 1;
                }
                if gt[beg] == b'.' && end - beg == 1 {
                    s.tmpi.push(GT_MISSING);
                } else {
                    if at(gt, end) == b'|' {
                        is_phased = true;
                    }
                    let (v, e) = strtol(gt, beg);
                    if e != end {
                        bail!(
                            "Could not parse the {} field at {} in {}",
                            lossy(&col.key_src),
                            site(ahdr, rec),
                            self.targets_fname
                        );
                    }
                    if v >= nals {
                        bail!(
                            "The source allele index is bigger than the number of destination alleles at {}",
                            site(ahdr, rec)
                        );
                    }
                    let v = v as i32;
                    s.tmpi.push(if is_phased {
                        gt_phased(v)
                    } else {
                        gt_unphased(v)
                    });
                }
                beg = if end < gt.len() { end + 1 } else { end };
            }
            s.tmpi.resize(start + ndst1, INT32_VECTOR_END);
        }
        let r = line.update_genotypes(&self.hdr_out, &self.s.tmpi);
        Ok(if r.is_ok() { Set::Done } else { Set::Failed })
    }
}

#[cfg(test)]
mod tests {
    use std::{fs, path::Path, process::ExitCode};

    use clap::{ArgMatches, Command, FromArgMatches};

    use super::*;
    use crate::commands::{annotate, view};

    #[test]
    fn remap_follows_the_allele_map() {
        // Source A,G onto destination A,C,G: REF to REF, G to the second ALT.
        let mut dst = [INT32_MISSING; 3];
        remap_values(&[0, 2], &[5, 7], &mut dst);
        assert_eq!(dst, [5, INT32_MISSING, 7]);
        let mut dst = [1, 2];
        remap_values(&[-1, 0, 9], &[5, 6, 7], &mut dst);
        assert_eq!(dst, [6, 2]);
    }

    #[test]
    fn remap_of_a_single_missing_value_writes_missing() {
        // A sample with AD=. (one value) at a site whose alleles were remapped.
        let mut dst = [1, 2, 3];
        remap_values(&[0, 2], &[INT32_MISSING], &mut dst);
        assert_eq!(dst, [INT32_MISSING, 2, INT32_MISSING]);
        // PL=. from source A,G onto destination A,C,G: AA, AG and GG are 0, 3 and 5.
        let mut dst = [1, 2, 3, 4, 5, 6];
        remap_values(&[0, 3, 5], &[INT32_MISSING], &mut dst);
        assert_eq!(dst, [INT32_MISSING, 2, 3, INT32_MISSING, 5, INT32_MISSING]);
        let mut dst = [1.5_f32; 2];
        remap_values(&[1, 0], &[f32::MISSING], &mut dst);
        assert!(dst.iter().all(|v| v.is_missing()), "{dst:?}");
    }

    #[test]
    fn genotypes_follow_the_allele_map() {
        let (u, p) = (gt_unphased, gt_phased);
        // Source A,C,G onto destination A,G: C is missing, G becomes 1.
        let mut gts = [
            u(1),
            p(2),
            u(0),
            GT_MISSING,
            u(2),
            INT32_VECTOR_END,
            u(5),
            1,
        ];
        remap_genotypes(&[0, -1, 1], &mut gts);
        assert_eq!(
            gts,
            [
                GT_MISSING,
                p(1),
                u(0),
                GT_MISSING,
                u(1),
                INT32_VECTOR_END,
                GT_MISSING,
                1
            ]
        );
        // A phased allele the destination lacks stays phased.
        let mut gts = [u(0), p(1)];
        remap_genotypes(&[], &mut gts);
        assert_eq!(gts, [GT_MISSING, 1]);
    }

    /// Runs the command `name` (its `Args` and `run`) with `argv`.
    fn run_command<A: clap::Args + FromArgMatches>(
        name: &'static str,
        argv: &[&str],
        run: fn(A, &ArgMatches) -> anyhow::Result<ExitCode>,
    ) -> anyhow::Result<()> {
        let matches = A::augment_args(Command::new(name))
            .try_get_matches_from(std::iter::once(name).chain(argv.iter().copied()))?;
        let code = run(A::from_arg_matches(&matches)?, &matches)?;
        anyhow::ensure!(code == ExitCode::SUCCESS, "{name} {argv:?} failed");
        Ok(())
    }

    /// Writes the header of the FORMAT tests and `line` as `dir/name.vcf.gz`, with an index.
    fn indexed_vcf(dir: &Path, name: &str, line: &str) -> anyhow::Result<String> {
        let text = dir.join(format!("{name}.vcf"));
        fs::write(
            &text,
            format!(
                "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                 ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n\
                 ##FORMAT=<ID=FR,Number=R,Type=Float,Description=\"Fractions\">\n\
                 ##FORMAT=<ID=XA,Number=A,Type=Integer,Description=\"Counts\">\n\
                 ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n{line}\n"
            ),
        )?;
        let gz = dir.join(format!("{name}.vcf.gz"));
        let gz = gz.to_str().expect("UTF-8 path").to_owned();
        let argv = ["--no-version", "-Oz", "--write-index=tbi", "-o", &gz];
        let argv = [&argv[..], &[text.to_str().expect("UTF-8 path")]].concat();
        run_command("view", &argv, view::run)?;
        Ok(gz)
    }

    /// The columns of the record that `annotate -a src -c column dst` writes.
    fn annotated_record(src: &str, column: &str, dst: &str) -> anyhow::Result<Vec<String>> {
        let out = format!("{dst}.{}.out", column.replace(['/', '+', ','], "_"));
        let argv = ["--no-version", "-a", src, "-c", column, "-o", &out, dst];
        run_command("annotate", &argv, annotate::run)?;
        let text = fs::read_to_string(&out)?;
        fs::remove_file(&out)?;
        let record = text.lines().last().unwrap_or_default();
        Ok(record.split('\t').map(str::to_owned).collect())
    }

    /// The sample columns of the record that `annotate -a src -c column dst` writes.
    fn annotated(src: &str, column: &str, dst: &str) -> anyhow::Result<Vec<String>> {
        let mut fields = annotated_record(src, column, dst)?;
        Ok(fields.split_off(9.min(fields.len())))
    }

    #[test]
    fn format_values_of_a_single_missing_value_are_missing_after_the_remap() -> anyhow::Result<()> {
        let dir = std::env::temp_dir().join(format!("vtools-fmt-dots-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        // Every source sample has a single '.' at a site whose alleles are in another order.
        let dots = indexed_vcf(
            &dir,
            "dots",
            "chr1\t10\t.\tA\tG,T\t.\t.\t.\tGT:AD:FR:XA:PL\t0/1:.:.:.:.\t0/1:.:.:.:.",
        )?;
        let dst = indexed_vcf(
            &dir,
            "dst",
            "chr1\t10\t.\tA\tT,G\t.\t.\t.\tGT:AD:PL\t0/1:1,2,3:1,2,3,4,5,6\t0/2:.:.",
        )?;
        let cases = [
            ("FMT/AD", ["0/1:.,.,.:1,2,3,4,5,6", "0/2:.,.,.:."]),
            ("+FMT/AD", ["0/1:.,.,.:1,2,3,4,5,6", "0/2:.,.,.:."]),
            ("FMT/FR", ["0/1:1,2,3:1,2,3,4,5,6:.,.,.", "0/2:.:.:.,.,."]),
            ("FMT/XA", ["0/1:1,2,3:1,2,3,4,5,6:.,.", "0/2:.:.:.,."]),
            ("FMT/PL", ["0/1:1,2,3:.,.,.,.,.,.", "0/2:.:.,.,.,.,.,."]),
        ];
        for (column, want) in cases {
            assert_eq!(annotated(&dots, column, &dst)?, want, "{column}");
        }
        // PL=. next to a haploid sample: the '.' sample is read as diploid from a row of three
        // values ('.' and two vector ends). The genotypes past the row are missing, and the
        // copied vector ends cut the values short, as in bcftools' output.
        let hap = indexed_vcf(
            &dir,
            "hap",
            "chr1\t10\t.\tA\tG,T\t.\t.\t.\tGT:PL\t1:1,2,3\t0/1:.",
        )?;
        let hap_dst = indexed_vcf(
            &dir,
            "hap_dst",
            "chr1\t10\t.\tA\tT,G\t.\t.\t.\tGT:PL\t1:7,8,9\t0/2:1,2,3,4,5,6",
        )?;
        assert_eq!(
            annotated(&hap, "FMT/PL", &hap_dst)?,
            ["1:1,3,2", "0/2:.,.,."]
        );
        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    #[test]
    fn copied_genotypes_name_the_same_alleles() -> anyhow::Result<()> {
        let dir = std::env::temp_dir().join(format!("vtools-fmt-gt-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let dst = indexed_vcf(&dir, "dst", "chr1\t10\t.\tA\tG,C\t.\t.\t.\tGT\t0/0\t0/0")?;
        // Source genotypes and what they become on ALT G,C: C is 2, G is 1, T is missing.
        let cases = [
            ("C,G", ["1/1", "0|2"], ["2/2", "0|1"]),
            ("C", ["0/1", "1|1"], ["0/2", "2|2"]),
            ("C,G,T", ["1/3", "3|2"], ["2/.", ".|1"]),
            ("C,G", ["1/2/.", "2"], ["2/1/.", "1"]),
            ("C,G", ["./.", ".|1"], ["./.", ".|2"]),
            ("T,C", ["0/1", "1|2"], ["0/.", ".|2"]),
            // The same alleles, and a prefix of them: the indexes are copied as they are.
            ("G,C", ["1/2", "2|0|1"], ["1/2", "2|0|1"]),
            ("G", ["1/1", "0|1"], ["1/1", "0|1"]),
        ];
        for (i, (alt, gts, want)) in cases.into_iter().enumerate() {
            let line = format!("chr1\t10\t.\tA\t{alt}\t.\t.\t.\tGT\t{}\t{}", gts[0], gts[1]);
            let src = indexed_vcf(&dir, &format!("src{i}"), &line)?;
            assert_eq!(
                annotated(&src, "FMT/GT", &dst)?,
                want,
                "ALT {alt}, GT {gts:?}"
            );
        }
        fs::remove_dir_all(&dir)?;
        Ok(())
    }

    #[test]
    fn copied_genotypes_follow_copied_alleles() -> anyhow::Result<()> {
        let dir = std::env::temp_dir().join(format!("vtools-fmt-gt-als-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let src = indexed_vcf(&dir, "src", "chr1\t10\t.\tA\tC,G\t.\t.\t.\tGT\t1/1\t0|2")?;
        let dst = indexed_vcf(&dir, "dst", "chr1\t10\t.\tA\tG,C\t.\t.\t.\tGT\t0/0\t0/0")?;
        let none = indexed_vcf(&dir, "none", "chr1\t10\t.\tA\t.\t.\t.\t.\tGT\t0/0\t0/0")?;
        // REF, ALT, and the genotypes, whichever comes first in the list: with the ALTs of the
        // annotation the indexes stay, with the record's own (+ALT finds some) they are mapped.
        let cases = [
            ("FMT/GT,REF,ALT", &dst, ["A", "C,G", "1/1", "0|2"]),
            ("FMT/GT,ALT", &dst, ["A", "C,G", "1/1", "0|2"]),
            ("ALT,FMT/GT", &dst, ["A", "C,G", "1/1", "0|2"]),
            ("FMT/GT,REF", &dst, ["A", "G,C", "2/2", "0|1"]),
            ("FMT/GT,+ALT", &dst, ["A", "G,C", "2/2", "0|1"]),
            ("+ALT,FMT/GT", &dst, ["A", "G,C", "2/2", "0|1"]),
            ("FMT/GT,+ALT", &none, ["A", "C,G", "1/1", "0|2"]),
            ("+ALT,FMT/GT", &none, ["A", "C,G", "1/1", "0|2"]),
            // ALT replaces the ALTs whether +ALT comes before or after it.
            ("FMT/GT,ALT,+ALT", &dst, ["A", "C,G", "1/1", "0|2"]),
            ("FMT/GT,+ALT,ALT", &dst, ["A", "C,G", "1/1", "0|2"]),
        ];
        for (column, dst, want) in cases {
            let got = annotated_record(&src, column, dst)?;
            let got = [3, 4, 9, 10].map(|i| got.get(i).map_or("", String::as_str));
            assert_eq!(got, want, "{column} onto {dst}");
        }
        // Paired at the position only: the REF of the annotation, or its ALTs, come along.
        let src = indexed_vcf(&dir, "src2", "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1|1")?;
        let dst = indexed_vcf(&dir, "dst2", "chr1\t10\t.\tC\tG\t.\t.\t.\tGT\t0/0\t0/0")?;
        for (column, want) in [
            ("FMT/GT,REF,ALT", ["A", "G", "0/1", "1|1"]),
            ("FMT/GT,ALT", ["C", "G", "0/1", "1|1"]),
        ] {
            let got = annotated_record(&src, column, &dst)?;
            let got = [3, 4, 9, 10].map(|i| got.get(i).map_or("", String::as_str));
            assert_eq!(got, want, "{column}");
        }
        fs::remove_dir_all(&dir)?;
        Ok(())
    }
}
