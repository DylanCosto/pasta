//! The FORMAT fields of an output line (bcftools merge's `merge_format`, `merge_GT`,
//! `merge_format_field`, `merge_format_string`, and the `-L` local-allele functions).

use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::bail;
use pasta::variant::{
    HeaderRecordType, Number, Record,
    encoding::{
        FLOAT_MISSING, FLOAT_VECTOR_END, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING,
        INT16_VECTOR_END, INT32_MISSING, INT32_VECTOR_END,
    },
    record::FormatField,
};

use super::{Merger, alleles::alleles2gt, maux::Inputs, maux::Maux, rules::MissingRule};

/// The BCF type codes (`BCF_BT_*`).
const BT_INT8: u8 = 1;
const BT_INT16: u8 = 2;
const BT_INT32: u8 = 3;
const BT_FLOAT: u8 = 5;
const BT_CHAR: u8 = 7;

/// The size of the PL-to-probability table (`PL2PROB_MAX`).
const PL2PROB_MAX: usize = 1024;

/// The largest row bcftools builds, in bytes.
const MAX_ROW_BYTES: usize = 2_147_483_647;

static WARNED_GT: AtomicBool = AtomicBool::new(false);
static WARNED_ROW: AtomicBool = AtomicBool::new(false);
static WARNED_STR: AtomicBool = AtomicBool::new(false);
static WARNED_LOCAL_G: AtomicBool = AtomicBool::new(false);
static WARNED_LOCAL_AR: AtomicBool = AtomicBool::new(false);
static WARNED_NO_PL: AtomicBool = AtomicBool::new(false);

/// The Number of a FORMAT tag as merging treats it (`BCF_VL_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Length {
    /// Fixed or variable: values are copied.
    Fixed,
    /// Number=., for strings.
    Var,
    A,
    G,
    R,
}

/// One value of a FORMAT field as HTSlib's typed macros read it.
#[derive(Clone, Copy, Debug)]
struct Val {
    /// The value sign-extended to 32 bits, or the float's bits.
    raw: i32,
    missing: bool,
    end: bool,
}

/// The values of a FORMAT field, all samples.
#[derive(Clone, Copy, Debug)]
struct Vals<'a> {
    ty: u8,
    data: &'a [u8],
    n: usize,
}

impl<'a> Vals<'a> {
    fn new(line: &'a Record, f: &'a FormatField) -> Self {
        Self {
            ty: f.type_code(),
            data: line.format_data(f),
            n: f.n(),
        }
    }

    /// Value `k` counted from the start of the field; `None` past its end.
    #[inline]
    fn at(&self, k: usize) -> Option<Val> {
        let d = self.data;
        Some(match self.ty {
            BT_INT8 => {
                let v = *d.get(k)? as i8;
                Val {
                    raw: i32::from(v),
                    missing: v == INT8_MISSING,
                    end: v == INT8_VECTOR_END,
                }
            }
            BT_INT16 => {
                let v = i16::from_le_bytes(d.get(k * 2..k * 2 + 2)?.try_into().ok()?);
                Val {
                    raw: i32::from(v),
                    missing: v == INT16_MISSING,
                    end: v == INT16_VECTOR_END,
                }
            }
            BT_INT32 => {
                let v = i32::from_le_bytes(d.get(k * 4..k * 4 + 4)?.try_into().ok()?);
                Val {
                    raw: v,
                    missing: v == INT32_MISSING,
                    end: v == INT32_VECTOR_END,
                }
            }
            BT_FLOAT => {
                let b = u32::from_le_bytes(d.get(k * 4..k * 4 + 4)?.try_into().ok()?);
                Val {
                    raw: b as i32,
                    missing: b == FLOAT_MISSING,
                    end: b == FLOAT_VECTOR_END,
                }
            }
            _ => return None,
        })
    }

    /// Value `k`, treating a value past the end of the field as a vector end.
    #[inline]
    fn get(&self, k: usize) -> Val {
        self.at(k).unwrap_or(Val {
            raw: 0,
            missing: false,
            end: true,
        })
    }

    #[inline]
    fn is_float(&self) -> bool {
        self.ty == BT_FLOAT
    }

    /// Whether `a < b` as values of this field's type.
    #[inline]
    fn less(&self, a: i32, b: i32) -> bool {
        if self.is_float() {
            f32::from_bits(a as u32) < f32::from_bits(b as u32)
        } else {
            a < b
        }
    }
}

/// The missing and vector-end words of an output field of a type.
#[inline]
fn sentinels(ty: u8) -> (u32, u32) {
    if ty == BT_FLOAT {
        (FLOAT_MISSING, FLOAT_VECTOR_END)
    } else {
        (INT32_MISSING as u32, INT32_VECTOR_END as u32)
    }
}

/// A constant of `-M` as a value of the output type (C's float to int conversion).
#[inline]
fn const_word(v: f32, ty: u8) -> u32 {
    if ty == BT_FLOAT {
        v.to_bits()
    } else {
        (v as i32) as u32
    }
}

/// What a reader's FORMAT tag is in the output, worked out once per dictionary ID.
#[derive(Clone, Copy, Debug)]
struct FmtKey {
    /// The output dictionary ID, `None` when the output header lacks the tag.
    out_id: Option<u32>,
    is_gt: bool,
    is_pl: bool,
    /// The tag's Number in the reader's header.
    number: Option<Number>,
}

/// Scratch space for the FORMAT fields.
#[derive(Debug, Default)]
pub(super) struct FormatState {
    /// Per reader, by dictionary ID.
    key_cache: Vec<Vec<Option<FmtKey>>>,
    /// For each output FORMAT row (0 is GT) and reader, the index of the field in the
    /// reader's record (`fmt_map`).
    map: Vec<Option<usize>>,
    /// The output dictionary ID of each row (`fmt_key`).
    keys: Vec<u32>,
    words: Vec<u32>,
    ints: Vec<i32>,
    floats: Vec<f32>,
    bytes: Vec<u8>,
    strs: Vec<Vec<u8>>,
    ac: Vec<i32>,
    /// Whether GT was set on the output line, with `ac` counted from it.
    gt_written: bool,
    /// One field of one record decoded to 32 bits.
    raw: Vec<i32>,
    /// One field of one record decoded to output words.
    decoded: Vec<u32>,
    /// `-L`: per sample, the local alleles (indexes of the reader's alleles), `L+1` each.
    laa: Vec<i32>,
    nlaa: usize,
    laa_dirty: bool,
    tmpd: Vec<f64>,
    tmpi: Vec<usize>,
    k2k: Vec<usize>,
    pl2prob: Vec<f64>,
}

impl Merger {
    /// Merges the FORMAT fields of the staged records (bcftools' `merge_format`).
    pub(super) fn merge_format(&mut self) -> anyhow::Result<()> {
        let n = self.ma.n;
        self.seen.clear();
        let fmt = &mut self.fmt;
        fmt.map.clear();
        fmt.map.resize(n, None);
        fmt.keys.clear();
        fmt.keys.push(0);
        fmt.gt_written = false;
        let mut has_gt = false;
        let mut has_pl: Option<usize> = None;
        let mut max_ifmt = 0;
        for i in 0..n {
            let Some(cur) = self.ma.buf[i].cur else {
                continue;
            };
            let src = self.ma.buf[i].src;
            let line = self.inp.line(src, i, cur);
            let hdr = self.inp.header(i);
            if fmt.key_cache.len() < n {
                fmt.key_cache.resize_with(n, Vec::new);
            }
            for (j, f) in line.format_fields().iter().enumerate() {
                if f.is_removed() {
                    continue;
                }
                let fk = fmt_key(
                    &mut fmt.key_cache[i],
                    &mut self.ids,
                    i,
                    f.id(),
                    hdr,
                    &self.out_hdr,
                );
                let Some(out_id) = fk.out_id else {
                    continue;
                };
                let ifmt = match self.seen.get(out_id) {
                    Some(v) => v as usize,
                    None => {
                        let ifmt = if fk.is_gt {
                            has_gt = true;
                            0
                        } else {
                            max_ifmt += 1;
                            fmt.map.resize((max_ifmt + 1) * n, None);
                            fmt.keys.push(out_id);
                            if fk.is_pl {
                                has_pl = Some(max_ifmt);
                            }
                            max_ifmt
                        };
                        self.seen.set(out_id, ifmt as i64);
                        ifmt
                    }
                };
                fmt.map[ifmt * n + i] = Some(j);
            }
            // Does the allele numbering change?
            let rec = &mut self.ma.buf[i].rec[cur];
            rec.als_differ = (1..line.allele_count()).any(|j| rec.map[j] != j);
        }

        if self.opts.local_alleles > 0 {
            fmt.laa_dirty = false;
            fmt.nlaa = 0;
            if self.out.allele_count() > self.opts.local_alleles + 1 {
                self.init_local_alleles(has_pl)?;
            }
        }

        self.out.set_sample_count(self.out_hdr.sample_count());
        if has_gt {
            self.merge_gt()?;
        }
        if !self.opts.keep_ac_an {
            self.update_an_ac()?;
        }
        for ifmt in 1..=max_ifmt {
            let key = self
                .out_hdr
                .id_name(self.fmt.keys[ifmt])
                .unwrap_or_default();
            let rule = self
                .missing_rules
                .binary_search_by(|(t, _)| t.as_slice().cmp(key))
                .ok()
                .map(|k| self.missing_rules[k].1);
            self.merge_format_field(ifmt, rule)?;
        }
        if self.fmt.laa_dirty {
            self.update_local_alleles();
        }
        Ok(())
    }

    /// The field of row `ifmt` (of the FORMAT map `map`) in reader `i`'s current record.
    fn fmt_field<'a>(
        map: &[Option<usize>],
        ma: &Maux,
        inp: &'a Inputs,
        ifmt: usize,
        i: usize,
    ) -> Option<(&'a Record, &'a FormatField)> {
        let j = map[ifmt * ma.n + i]?;
        let line = ma.line(inp, i)?;
        Some((line, &line.format_fields()[j]))
    }

    /// Sets INFO/AN and INFO/AC from the merged genotypes (bcftools' `update_AN_AC`).
    fn update_an_ac(&mut self) -> anyhow::Result<()> {
        let ac = &self.fmt.ac;
        if self.fmt.gt_written {
            let n_allele = self.out.allele_count();
            let an: i32 = ac[..n_allele.max(1)].iter().sum();
            if let Some(id) = self.ids.an {
                let _ = self.out.update_info_i32_id(&self.out_hdr, id, &[an]);
            }
            if let Some(id) = self.ids.ac {
                let _ = self
                    .out
                    .update_info_i32_id(&self.out_hdr, id, &ac[1..n_allele.max(1)]);
            }
        }
        Ok(())
    }

    /// Merges the genotypes, renumbering alleles (bcftools' `merge_GT`).
    fn merge_gt(&mut self) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            out,
            out_hdr,
            fmt,
            opts,
            ..
        } = self;
        let n = ma.n;
        let nsamples = out_hdr.sample_count();
        let mut nsize = 0;
        for i in 0..n {
            let Some((line, f)) = Self::fmt_field(&fmt.map, ma, inp, 0, i) else {
                continue;
            };
            let v = Vals::new(line, f);
            if !matches!(v.ty, BT_INT8 | BT_INT16 | BT_INT32) {
                bail!("Unexpected case: {}", v.ty);
            }
            nsize = nsize.max(max_used_gt_ploidy(&v, inp.header(i).sample_count()));
        }
        if nsize == 0 {
            nsize = 1;
        }
        let msize = 4 * nsize * nsamples;
        if msize > MAX_ROW_BYTES {
            if !WARNED_GT.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "Warning: Too many genotypes at {}:{}, requires {msize} bytes, skipping.",
                    String::from_utf8_lossy(out.chrom_safe(out_hdr)),
                    out.pos() + 1
                );
            }
            return Ok(());
        }
        let tmp = &mut fmt.ints;
        tmp.clear();
        tmp.resize(nsize * nsamples, 0);
        ma.smpl_ploidy.clear();
        ma.smpl_ploidy.resize(nsamples, 0);
        let default_gt = if opts.missing_to_ref { 2 } else { 0 };
        let mut ismpl = 0;
        for i in 0..n {
            let nsmpl = inp.header(i).sample_count();
            let Some((line, f)) = Self::fmt_field(&fmt.map, ma, inp, 0, i) else {
                // Missing values: assume the maximum ploidy.
                for j in 0..nsmpl {
                    let t = &mut tmp[(ismpl + j) * nsize..(ismpl + j + 1) * nsize];
                    t.fill(default_gt);
                    ma.smpl_ploidy[ismpl + j] += nsize as i32;
                }
                ismpl += nsmpl;
                continue;
            };
            let v = Vals::new(line, f);
            decode_raw(&v, nsmpl, &mut fmt.raw);
            let cur = ma.buf[i].cur.expect("the reader has a record");
            let rec = &ma.buf[i].rec[cur];
            let vn = v.n;
            let remap = |raw: i32| -> i32 {
                if raw >> 1 == 0 {
                    raw
                } else {
                    let al = (raw >> 1) - 1;
                    let al = if al <= 0 {
                        al + 1
                    } else {
                        rec.map.get(al as usize).map_or(al, |&m| m as i32) + 1
                    };
                    (al << 1) | (raw & 1)
                }
            };
            if vn == nsize && !fmt.raw.contains(&INT32_VECTOR_END) {
                // Every sample has the full ploidy: no padding to handle.
                let t = &mut tmp[ismpl * nsize..(ismpl + nsmpl) * nsize];
                if rec.als_differ {
                    for (t, &raw) in t.iter_mut().zip(&fmt.raw) {
                        *t = remap(raw);
                    }
                } else {
                    t.copy_from_slice(&fmt.raw);
                }
                for p in &mut ma.smpl_ploidy[ismpl..ismpl + nsmpl] {
                    *p += vn as i32;
                }
                ismpl += nsmpl;
                continue;
            }
            for j in 0..nsmpl {
                let t = &mut tmp[(ismpl + j) * nsize..(ismpl + j + 1) * nsize];
                let src = &fmt.raw[j * vn..(j + 1) * vn];
                let k = src
                    .iter()
                    .position(|&r| r == INT32_VECTOR_END)
                    .unwrap_or(vn);
                ma.smpl_ploidy[ismpl + j] += k as i32;
                if rec.als_differ {
                    for (t, &raw) in t.iter_mut().zip(&src[..k]) {
                        *t = remap(raw);
                    }
                } else {
                    t[..k].copy_from_slice(&src[..k]);
                }
                t[k..].fill(INT32_VECTOR_END);
            }
            ismpl += nsmpl;
        }
        if out.update_genotypes(out_hdr, tmp).is_err() {
            return Ok(());
        }
        // AC and AN as bcf_calc_ac counts them from the GT just set.
        let n_allele = out.allele_count();
        let ac = &mut fmt.ac;
        ac.clear();
        ac.resize(n_allele.max(1), 0);
        // Each sample's alleles are followed only by vector ends, so one pass over all
        // values counts them; most sites have only alleles 0 and 1.
        let (mut cnt, mut sum, mut max) = (0i32, 0i32, 0i32);
        for &gt in tmp.iter() {
            let a = (gt >> 1) - 1;
            let valid = gt != INT32_VECTOR_END && a >= 0;
            cnt += i32::from(valid);
            sum += if valid { a } else { 0 };
            max = max.max(if valid { a } else { 0 });
        }
        if max <= 1 && n_allele >= 2 {
            ac[0] = cnt - sum;
            ac[1] = sum;
        } else {
            for (s, gts) in tmp.chunks_exact(nsize).enumerate() {
                for &gt in gts {
                    let a = (gt >> 1) - 1;
                    if gt == INT32_VECTOR_END || a < 0 {
                        continue;
                    }
                    if a as usize >= n_allele {
                        eprintln!(
                            "[E::bcf_calc_ac] Incorrect allele (\"{a}\") in {} at {}:{}",
                            String::from_utf8_lossy(&out_hdr.samples()[s]),
                            String::from_utf8_lossy(out.chrom_safe(out_hdr)),
                            out.pos() + 1
                        );
                        return Err(super::ExitOne.into());
                    }
                    ac[a as usize] += 1;
                }
            }
        }
        fmt.gt_written = true;
        Ok(())
    }

    /// Merges one FORMAT field, remapping Number=A/G/R values to the output alleles
    /// (bcftools' `merge_format_field`).
    #[expect(
        clippy::too_many_lines,
        clippy::needless_range_loop,
        reason = "a line-for-line port of merge_format_field"
    )]
    fn merge_format_field(
        &mut self,
        ifmt: usize,
        mrule: Option<MissingRule>,
    ) -> anyhow::Result<()> {
        let n = self.ma.n;
        let out_id = self.fmt.keys[ifmt];
        let mut length = Length::Fixed;
        let mut nsize = 0usize;
        let mut ty: Option<u8> = None;
        let mut irdr = n;
        let out_nals = self.out.allele_count();
        for i in 0..n {
            let Some((_, f)) = Self::fmt_field(&self.fmt.map, &self.ma, &self.inp, ifmt, i) else {
                continue;
            };
            ty = Some(f.type_code());
            irdr = i;
            let number = self.fmt.key_cache[i]
                .get(f.id() as usize)
                .copied()
                .flatten()
                .and_then(|k| k.number);
            match number {
                Some(Number::G) => {
                    length = Length::G;
                    nsize = out_nals * (out_nals + 1) / 2;
                    break;
                }
                Some(Number::A) => {
                    length = Length::A;
                    nsize = out_nals.saturating_sub(1);
                    break;
                }
                Some(Number::R) => {
                    length = Length::R;
                    nsize = out_nals;
                    break;
                }
                Some(Number::Variable) if f.type_code() == BT_CHAR => {
                    nsize = nsize.max(f.n());
                    length = Length::Var;
                }
                _ => nsize = nsize.max(f.n()),
            }
        }
        let Some(mut ty) = ty else {
            return Ok(());
        };
        if self.fmt.nlaa > 0 && matches!(length, Length::A | Length::G | Length::R) {
            return self.merge_localized(irdr, ifmt, length);
        }
        if ty == BT_CHAR {
            return self.merge_format_string(ifmt, length, nsize);
        }

        let Self {
            inp,
            ma,
            out,
            out_hdr,
            fmt,
            ..
        } = self;
        let key = out_hdr.id_name(out_id).unwrap_or_default();
        let chrom = || String::from_utf8_lossy(out.chrom_safe(out_hdr)).into_owned();
        let nsamples = out_hdr.sample_count();
        let msize = 4 * nsamples * nsize;
        if msize > MAX_ROW_BYTES {
            if !WARNED_ROW.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "Warning: The row size is too big for FORMAT/{} at {}:{}, requires {msize} bytes, skipping.",
                    String::from_utf8_lossy(key),
                    chrom(),
                    out.pos() + 1
                );
            }
            return Ok(());
        }
        let words = &mut fmt.words;
        words.clear();
        words.resize(nsamples * nsize, 0);

        let mut ismpl = 0;
        for i in 0..n {
            let nsmpl = inp.header(i).sample_count();
            let field = Self::fmt_field(&fmt.map, ma, inp, ifmt, i);
            if let Some((line, f)) = field {
                ty = f.type_code();
                let nals_ori = line.allele_count();
                let n = f.n();
                let bad = match length {
                    Length::G => n != 1 && n != nals_ori * (nals_ori + 1) / 2 && n != nals_ori,
                    Length::A => n != 1 && n + 1 != nals_ori,
                    Length::R => n != 1 && n != nals_ori,
                    _ => false,
                };
                if bad {
                    let what = match length {
                        Length::G => "G",
                        Length::A => "A",
                        _ => "R",
                    };
                    bail!(
                        "Incorrect number of FORMAT/{} values at {}:{}, cannot merge. The tag is defined as Number={what}, but found\n{n} values and {nals_ori} alleles. See also http://samtools.github.io/bcftools/howtos/FAQ.html#incorrect-nfields",
                        String::from_utf8_lossy(key),
                        chrom(),
                        out.pos() + 1
                    );
                }
            }
            if !matches!(ty, BT_INT8 | BT_INT16 | BT_INT32 | BT_FLOAT) {
                bail!("Unexpected case: {ty}, {}", String::from_utf8_lossy(key));
            }
            let (missing, vector_end) = sentinels(ty);
            let Some((line, f)) = field else {
                // The field is not in this file.
                for j in 0..nsmpl {
                    let t = &mut words[(ismpl + j) * nsize..(ismpl + j + 1) * nsize];
                    if let Some((first, rest)) = t.split_first_mut() {
                        *first = missing;
                        rest.fill(vector_end);
                    }
                }
                ismpl += nsmpl;
                continue;
            };
            let v = Vals::new(line, f);
            let cur = ma.buf[i].cur.expect("the reader has a record");
            let rec = &ma.buf[i].rec[cur];
            let unchanged = matches!(length, Length::Fixed | Length::Var)
                || (line.allele_count() == out_nals && !rec.als_differ);
            if unchanged {
                decode_words(&v, nsmpl, &mut fmt.decoded);
                let vn = v.n;
                for j in 0..nsmpl {
                    let t = &mut words[(ismpl + j) * nsize..(ismpl + j + 1) * nsize];
                    let src = &fmt.decoded[j * vn..(j + 1) * vn];
                    let l = src
                        .iter()
                        .position(|&w| w == vector_end)
                        .unwrap_or(vn)
                        .min(nsize);
                    t[..l].copy_from_slice(&src[..l]);
                    t[l..].fill(vector_end);
                }
                ismpl += nsmpl;
                continue;
            }
            let unkn = ma.buf[i].unkn_allele;
            let map = &rec.map;
            for j in 0..nsmpl {
                let t = &mut words[(ismpl + j) * nsize..(ismpl + j + 1) * nsize];
                let base = j * v.n;
                let first = v.get(base);
                let tag_missing = if v.n == 1 {
                    first.missing
                } else {
                    first.missing && v.get(base + 1).end
                };
                if tag_missing {
                    t[0] = missing;
                    t[1..].fill(vector_end);
                    continue;
                }
                let sample_max = || {
                    let mut max = v.get(base).raw;
                    for l in 1..v.n {
                        let val = v.get(base + l).raw;
                        if v.less(max, val) {
                            max = val;
                        }
                    }
                    max as u32
                };
                if length == Length::G {
                    let haploid = ma.smpl_ploidy[ismpl + j] == 1;
                    let ngsize = if haploid {
                        out_nals
                    } else {
                        out_nals * (out_nals + 1) / 2
                    };
                    let fill = if unkn > 0 {
                        // Use the value of the unknown allele when there is one.
                        let iunkn = if haploid {
                            unkn
                        } else {
                            (unkn + 1) * (unkn + 2) / 2 - 1
                        };
                        v.get(base + iunkn).raw as u32
                    } else {
                        match mrule {
                            Some(MissingRule::Const(c)) => const_word(c, ty),
                            Some(MissingRule::Max) => sample_max(),
                            _ => missing,
                        }
                    };
                    let ng = ngsize.min(nsize);
                    t[..ng].fill(fill);
                    t[ng..].fill(vector_end);
                    if haploid {
                        for iori in 0..line.allele_count() {
                            let val = v.get(base + iori);
                            if val.end {
                                break;
                            }
                            if let Some(slot) = t.get_mut(map[iori]) {
                                *slot = if val.missing { missing } else { val.raw as u32 };
                            }
                        }
                    } else {
                        'outer: for iori in 0..line.allele_count() {
                            let inew = map[iori];
                            for jori in 0..=iori {
                                let jnew = map[jori];
                                let kori = iori * (iori + 1) / 2 + jori;
                                let knew = alleles2gt(inew, jnew);
                                let val = v.get(base + kori);
                                if val.end {
                                    break 'outer;
                                }
                                if let Some(slot) = t.get_mut(knew) {
                                    *slot = if val.missing { missing } else { val.raw as u32 };
                                }
                            }
                        }
                    }
                } else {
                    let ifrom = usize::from(length == Length::A);
                    let fill = if unkn > 0 {
                        // bcftools reads index `unkn` for Number=A too, one past the unknown
                        // allele's value (past the sample's values when it is the last).
                        v.get(base + unkn - ifrom).raw as u32
                    } else {
                        match mrule {
                            Some(MissingRule::Const(c)) => const_word(c, ty),
                            Some(MissingRule::Max) => sample_max(),
                            _ => missing,
                        }
                    };
                    t.fill(fill);
                    for (k, iori) in (ifrom..line.allele_count()).enumerate() {
                        let val = v.get(base + k);
                        if val.end {
                            break;
                        }
                        if let Some(slot) = t.get_mut(map[iori] - ifrom) {
                            *slot = if val.missing { missing } else { val.raw as u32 };
                        }
                    }
                }
            }
            ismpl += nsmpl;
        }
        if ty == BT_FLOAT {
            fmt.floats.clear();
            fmt.floats.extend(words.iter().map(|&w| f32::from_bits(w)));
            let _ = out.update_format_f32_id(out_hdr, out_id, &fmt.floats);
        } else {
            fmt.ints.clear();
            fmt.ints.extend(words.iter().map(|&w| w as i32));
            let _ = out.update_format_i32_id(out_hdr, out_id, &fmt.ints);
        }
        Ok(())
    }

    /// Merges a string FORMAT field (bcftools' `merge_format_string`).
    fn merge_format_string(
        &mut self,
        ifmt: usize,
        length: Length,
        nsize: usize,
    ) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            out,
            out_hdr,
            fmt,
            ..
        } = self;
        let out_id = fmt.keys[ifmt];
        let key = out_hdr.id_name(out_id).unwrap_or_default();
        let nsamples = out_hdr.sample_count();
        if fmt.strs.len() < nsamples {
            fmt.strs.resize_with(nsamples, Vec::new);
        }
        let mut nmax = 0;
        for s in &mut fmt.strs[..nsamples] {
            s.clear();
            s.push(b'.');
            if !matches!(length, Length::Fixed | Length::Var) {
                for _ in 1..nsize {
                    s.extend_from_slice(b",.");
                }
            }
            nmax = nmax.max(s.len());
        }
        let out_nals = out.allele_count();
        let mut ismpl = 0;
        for i in 0..ma.n {
            let nsmpl = inp.header(i).sample_count();
            let Some((line, f)) = Self::fmt_field(&fmt.map, ma, inp, ifmt, i) else {
                ismpl += nsmpl;
                continue;
            };
            let data = line.format_data(f);
            let size = f.size();
            let cur = ma.buf[i].cur.expect("the reader has a record");
            let rec = &ma.buf[i].rec[cur];
            if matches!(length, Length::Fixed | Length::Var)
                || (line.allele_count() == out_nals && !rec.als_differ)
            {
                for j in 0..nsmpl {
                    let s = &mut fmt.strs[ismpl];
                    ismpl += 1;
                    s.clear();
                    s.extend_from_slice(&data[j * size..(j + 1) * size]);
                    nmax = nmax.max(s.len());
                }
                continue;
            }
            if length == Length::G {
                bail!(
                    "[E::merge_format_string] Merging of Number=G FORMAT strings (in your case FORMAT/{}) is not supported yet, sorry!\nPlease open an issue on github if this feature is essential for you. However, note that using FORMAT strings is not\na good idea in general - it is slow to parse and does not compress well, it is better to use integer codes instead.\nIf you don't really need it, use `bcftools annotate -x` to remove the annotation before merging.",
                    String::from_utf8_lossy(key)
                );
            }
            let ifrom = usize::from(length == Length::A);
            for j in 0..nsmpl {
                let s = &mut fmt.strs[ismpl];
                ismpl += 1;
                let src = &data[j * size..(j + 1) * size];
                for iori in ifrom..line.allele_count() {
                    let inew = rec.map[iori] - ifrom;
                    let ret = super::site::copy_string_field(src, iori - ifrom, s, inew);
                    if ret < -1 {
                        bail!(
                            "[E::merge_format_string] fixme: internal error at {}:{} .. {ret}",
                            String::from_utf8_lossy(line.chrom_safe(inp.header(i))),
                            line.pos() + 1
                        );
                    }
                }
                nmax = nmax.max(s.len());
            }
        }
        let msize = nsamples * nmax;
        if msize > MAX_ROW_BYTES {
            if !WARNED_STR.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "Warning: The row size is too big for FORMAT/{} at {}:{}, requires {msize} bytes, skipping.",
                    String::from_utf8_lossy(key),
                    String::from_utf8_lossy(out.chrom_safe(out_hdr)),
                    out.pos() + 1
                );
            }
            return Ok(());
        }
        let bytes = &mut fmt.bytes;
        bytes.clear();
        bytes.resize(msize, 0);
        for (k, s) in fmt.strs[..nsamples].iter().enumerate() {
            bytes[k * nmax..k * nmax + s.len()].copy_from_slice(s);
        }
        let _ = out.update_format_chars(out_hdr, key, bytes);
        Ok(())
    }

    /// Chooses each sample's local alleles: all of a file's alleles when few enough, else
    /// the most likely by PL (bcftools' `init_local_alleles`).
    fn init_local_alleles(&mut self, ifmt_pl: Option<usize>) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            fmt,
            opts,
            out,
            ..
        } = self;
        let la = opts.local_alleles;
        let n_out = out.allele_count();
        if fmt.pl2prob.is_empty() {
            fmt.pl2prob = (0..PL2PROB_MAX)
                .map(|i| 10f64.powf(-0.1 * i as f64))
                .collect();
        }
        fmt.laa.clear();
        fmt.laa.resize(ma.nout_smpl * (1 + la), 0);
        fmt.tmpd.resize(n_out.max(1), 0.0);
        fmt.tmpi.resize(n_out.max(1), 0);
        fmt.k2k.resize(n_out.max(1), 0);
        let mut nlaa = 0;
        let mut ismpl = 0;
        for i in 0..ma.n {
            let nsmpl = inp.header(i).sample_count();
            let line = ma.line(inp, i);
            if let Some(line) = line {
                let nal = line.allele_count();
                if nlaa < nal - 1 {
                    nlaa = (nal - 1).min(la);
                }
                if nal <= la + 1 {
                    // Few enough alleles: all of them, in the output order.
                    let map = ma.cur_map(i);
                    let k2k = &mut fmt.k2k;
                    k2k.clear();
                    k2k.extend(0..nal);
                    for k in 1..nal {
                        let mut l = k;
                        while l > 0 && map[k2k[l]] < map[k2k[l - 1]] {
                            k2k.swap(l, l - 1);
                            l -= 1;
                        }
                    }
                    for _ in 0..nsmpl {
                        let ptr = &mut fmt.laa[(1 + la) * ismpl..(1 + la) * (ismpl + 1)];
                        for k in 0..nal {
                            ptr[k] = k2k[k] as i32;
                        }
                        ptr[nal..].fill(INT32_VECTOR_END);
                        ismpl += 1;
                    }
                    continue;
                }
            }
            let fmt_ori = ifmt_pl.and_then(|r| Self::fmt_field(&fmt.map, ma, inp, r, i));
            // Without PL, bcftools reads a stray pointer here and gives up (below).
            let absent = match ifmt_pl {
                Some(_) => fmt_ori.is_none(),
                None => false,
            };
            if line.is_none() || absent {
                for _ in 0..nsmpl {
                    let ptr = &mut fmt.laa[(1 + la) * ismpl..(1 + la) * (ismpl + 1)];
                    ptr[0] = INT32_MISSING;
                    ptr[1..].fill(INT32_VECTOR_END);
                    ismpl += 1;
                }
                continue;
            }
            let Some((line, f)) = fmt_ori else {
                if !WARNED_NO_PL.swap(true, Ordering::Relaxed) {
                    eprintln!(
                        "Warning: local alleles are determined from FORMAT/PL but the tag is missing, cannot apply --local-alleles"
                    );
                }
                fmt.nlaa = 0;
                return Ok(());
            };
            let hdr = inp.header(i);
            if hdr.id_number(HeaderRecordType::Format, f.id() as u32) != Some(Number::G) {
                bail!("FORMAT/PL must be defined as Number=G");
            }
            let nal = line.allele_count();
            if 2 * f.n() != nal * (nal + 1) {
                bail!("Todo: haploid PL to LPL");
            }
            let v = Vals::new(line, f);
            if !matches!(v.ty, BT_INT8 | BT_INT16 | BT_INT32) {
                bail!("Unexpected case: {}, PL", v.ty);
            }
            // An input line with an ALT given twice has more alleles than the output line.
            if fmt.tmpd.len() < nal {
                fmt.tmpd.resize(nal, 0.0);
                fmt.tmpi.resize(nal, 0);
            }
            let map = ma.cur_map(i);
            let mut src = 0;
            for _ in 0..nsmpl {
                let prob = &mut fmt.tmpd;
                prob[..nal].fill(0.0);
                for k in 0..nal {
                    for l in 0..=k {
                        let val = v.get(src);
                        if !val.missing && !val.end {
                            let idx = if val.raw >= 0 && (val.raw as usize) < PL2PROB_MAX {
                                val.raw as usize
                            } else {
                                PL2PROB_MAX - 1
                            };
                            let p = fmt.pl2prob[idx];
                            prob[k] += p;
                            prob[l] += p;
                        }
                        src += 1;
                    }
                }
                // Sort the ALTs by probability, descending; REF stays first.
                let idx = &mut fmt.tmpi;
                for (si, x) in idx[..nal - 1].iter_mut().enumerate() {
                    *x = si;
                }
                for si in 1..nal - 1 {
                    let mut sj = si;
                    while sj > 0 && prob[1 + idx[sj]] > prob[1 + idx[sj - 1]] {
                        idx.swap(sj, sj - 1);
                        sj -= 1;
                    }
                }
                let ptr = &mut fmt.laa[(1 + la) * ismpl..(1 + la) * (ismpl + 1)];
                ptr[0] = 0;
                let mut k = 1;
                while k <= la && k < nal {
                    ptr[k] = idx[k - 1] as i32 + 1;
                    k += 1;
                }
                let kmax = k;
                ptr[k..].fill(INT32_VECTOR_END);
                // Sort by output allele order, ascending.
                for k in 1..kmax {
                    let mut l = k;
                    while l > 0 && map[ptr[l] as usize] < map[ptr[l - 1] as usize] {
                        ptr.swap(l, l - 1);
                        l -= 1;
                    }
                }
                ismpl += 1;
            }
        }
        fmt.nlaa = nlaa;
        Ok(())
    }

    /// Writes FORMAT/LAA: each sample's local alleles as output allele indexes (bcftools'
    /// `update_local_alleles`).
    fn update_local_alleles(&mut self) {
        let Self {
            inp,
            ma,
            fmt,
            opts,
            out,
            out_hdr,
            ..
        } = self;
        let la = opts.local_alleles;
        let nlaa = fmt.nlaa;
        let mut ismpl = 0;
        for i in 0..ma.n {
            let nsmpl = inp.header(i).sample_count();
            let map = ma.buf[i].cur.map(|c| &ma.buf[i].rec[c].map);
            for _ in 0..nsmpl {
                let src = ismpl * (1 + la);
                let dst = ismpl * nlaa;
                let mut j = 0;
                if let Some(map) = map {
                    while j < nlaa {
                        let s = fmt.laa[src + j + 1];
                        if s == INT32_MISSING {
                            fmt.laa[dst + j] = INT32_MISSING;
                        } else if s == INT32_VECTOR_END {
                            break;
                        } else {
                            fmt.laa[dst + j] = map[s as usize] as i32;
                        }
                        j += 1;
                    }
                }
                if j == 0 {
                    fmt.laa[dst] = INT32_MISSING;
                    j = 1;
                }
                fmt.laa[dst + j..dst + nlaa.max(j)].fill(INT32_VECTOR_END);
                ismpl += 1;
            }
        }
        let nsamples = out_hdr.sample_count();
        let _ = out.update_format_i32(out_hdr, b"LAA", &fmt.laa[..nsamples * nlaa]);
    }

    /// Writes the localized `L<tag>` of a Number=G, A, or R field in row `row`, whose Number
    /// reader `irdr` gives (bcftools' `merge_localized_numberG_format_field` and
    /// `merge_localized_numberAR_format_field`).
    #[expect(
        clippy::needless_range_loop,
        reason = "a port of the localizing functions, indexing local alleles"
    )]
    fn merge_localized(&mut self, irdr: usize, row: usize, length: Length) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            fmt,
            opts,
            out,
            out_hdr,
            ..
        } = self;
        let la = opts.local_alleles;
        let nlaa = fmt.nlaa;
        let (_, f0) =
            Self::fmt_field(&fmt.map, ma, inp, row, irdr).expect("the reader has the field");
        let key = inp.header(irdr).id_name(f0.id() as u32).unwrap_or_default();
        let ty0 = f0.type_code();
        let nsize = match length {
            Length::G => (nlaa + 1) * (nlaa + 2) / 2,
            Length::R => nlaa + 1,
            _ => nlaa,
        };
        let nsamples = out_hdr.sample_count();
        let msize = 4 * nsamples * nsize;
        if msize > MAX_ROW_BYTES {
            let warned = if length == Length::G {
                &WARNED_LOCAL_G
            } else {
                &WARNED_LOCAL_AR
            };
            if !warned.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "Warning: The row size is too big for FORMAT/{} at {}:{}, requires {msize} bytes, skipping.",
                    String::from_utf8_lossy(key),
                    String::from_utf8_lossy(out.chrom_safe(out_hdr)),
                    out.pos() + 1
                );
            }
            return Ok(());
        }
        let words = &mut fmt.words;
        words.clear();
        words.resize(nsamples * nsize, 0);
        let ibeg = usize::from(length == Length::A);
        let mut ismpl = 0;
        for i in 0..ma.n {
            let nsmpl = inp.header(i).sample_count();
            let Some((line, f)) = Self::fmt_field(&fmt.map, ma, inp, row, i) else {
                if !matches!(ty0, BT_INT8 | BT_INT16 | BT_INT32 | BT_FLOAT) {
                    bail!("Unexpected case: {ty0}, {}", String::from_utf8_lossy(key));
                }
                let (missing, vector_end) = sentinels(ty0);
                for _ in 0..nsmpl {
                    let t = &mut words[ismpl * nsize..(ismpl + 1) * nsize];
                    if let Some((first, rest)) = t.split_first_mut() {
                        *first = missing;
                        rest.fill(vector_end);
                    }
                    ismpl += 1;
                }
                continue;
            };
            let nal = line.allele_count();
            if length == Length::G && 2 * f.n() != nal * (nal + 1) {
                bail!("Todo: localization of missing or haploid Number=G tags");
            }
            let v = Vals::new(line, f);
            if !matches!(v.ty, BT_INT8 | BT_INT16 | BT_INT32 | BT_FLOAT) {
                bail!(
                    "Unexpected case: {}, {}",
                    v.ty,
                    String::from_utf8_lossy(key)
                );
            }
            let (missing, vector_end) = sentinels(v.ty);
            for j in 0..nsmpl {
                let base = j * v.n;
                let t = &mut words[ismpl * nsize..(ismpl + 1) * nsize];
                let laa = &fmt.laa[(1 + la) * ismpl..(1 + la) * (ismpl + 1)];
                let mut tgt_idx = 0;
                let is_local = |a: i32| a != INT32_MISSING && a != INT32_VECTOR_END;
                if length == Length::G {
                    for ii in 0..=nlaa {
                        if !is_local(laa[ii]) {
                            break;
                        }
                        for ij in 0..=ii {
                            let src_idx = alleles2gt(laa[ii] as usize, laa[ij] as usize);
                            let val = v.get(base + src_idx);
                            if val.missing {
                                t[tgt_idx] = missing;
                            } else if val.end {
                                // Only this allele's genotypes end (a haploid sample's
                                // values), as in bcftools.
                                break;
                            } else {
                                t[tgt_idx] = val.raw as u32;
                            }
                            tgt_idx += 1;
                        }
                    }
                } else {
                    for ii in ibeg..=nlaa {
                        if !is_local(laa[ii]) {
                            break;
                        }
                        let src_idx = laa[ii] as usize - ibeg;
                        let val = v.get(base + src_idx);
                        if val.missing {
                            t[tgt_idx] = missing;
                        } else if val.end {
                            break;
                        } else {
                            t[tgt_idx] = val.raw as u32;
                        }
                        tgt_idx += 1;
                    }
                }
                if tgt_idx == 0 {
                    t[0] = missing;
                    tgt_idx = 1;
                }
                t[tgt_idx..].fill(vector_end);
                ismpl += 1;
            }
        }
        let mut name = Vec::with_capacity(key.len() + 1);
        name.push(b'L');
        name.extend_from_slice(key);
        if ty0 == BT_FLOAT {
            fmt.floats.clear();
            fmt.floats.extend(words.iter().map(|&w| f32::from_bits(w)));
            let _ = out.update_format_f32(out_hdr, &name, &fmt.floats);
        } else {
            fmt.ints.clear();
            fmt.ints.extend(words.iter().map(|&w| w as i32));
            let _ = out.update_format_i32(out_hdr, &name, &fmt.ints);
        }
        fmt.laa_dirty = true;
        Ok(())
    }
}

/// Returns the cached [`FmtKey`] of reader `i`'s FORMAT dictionary ID `id`.
#[inline]
fn fmt_key(
    cache: &mut Vec<Option<FmtKey>>,
    ids: &mut super::IdMaps,
    i: usize,
    id: i32,
    hdr: &pasta::variant::Header,
    out_hdr: &pasta::variant::Header,
) -> FmtKey {
    let k = id as usize;
    if k >= cache.len() {
        cache.resize(k + 1, None);
    }
    if let Some(fk) = cache[k] {
        return fk;
    }
    let name = hdr.id_name(id as u32).unwrap_or_default();
    let fk = FmtKey {
        out_id: ids.out_id(i, id, hdr, out_hdr),
        is_gt: name == b"GT",
        is_pl: name == b"PL",
        number: hdr.id_number(HeaderRecordType::Format, id as u32),
    };
    cache[k] = Some(fk);
    fk
}

/// Decodes the first `nsmpl` samples of an integer field to 32 bits as bcftools' typed macros
/// read them: values sign-extended (missing values included), vector ends as
/// [`INT32_VECTOR_END`].
fn decode_raw(v: &Vals<'_>, nsmpl: usize, out: &mut Vec<i32>) {
    let n = v.n * nsmpl;
    out.clear();
    match v.ty {
        BT_INT8 => out.extend(v.data[..n].iter().map(|&b| {
            let x = b as i8;
            if x == INT8_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(x)
            }
        })),
        BT_INT16 => out.extend(v.data[..n * 2].as_chunks::<2>().0.iter().map(|c| {
            let x = i16::from_le_bytes(*c);
            if x == INT16_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(x)
            }
        })),
        BT_INT32 => out.extend(
            v.data[..n * 4]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|c| i32::from_le_bytes(*c)),
        ),
        _ => out.resize(n, INT32_VECTOR_END),
    }
}

/// Decodes the first `nsmpl` samples of a numeric field to output words: integers widened
/// to 32 bits with their missing values and vector ends, floats as they are.
fn decode_words(v: &Vals<'_>, nsmpl: usize, out: &mut Vec<u32>) {
    let n = v.n * nsmpl;
    out.clear();
    match v.ty {
        BT_INT8 => out.extend(v.data[..n].iter().map(|&b| {
            let x = b as i8;
            (match x {
                INT8_MISSING => INT32_MISSING,
                INT8_VECTOR_END => INT32_VECTOR_END,
                _ => i32::from(x),
            }) as u32
        })),
        BT_INT16 => out.extend(v.data[..n * 2].as_chunks::<2>().0.iter().map(|c| {
            let x = i16::from_le_bytes(*c);
            (match x {
                INT16_MISSING => INT32_MISSING,
                INT16_VECTOR_END => INT32_VECTOR_END,
                _ => i32::from(x),
            }) as u32
        })),
        _ => out.extend(
            v.data[..n * 4]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|c| u32::from_le_bytes(*c)),
        ),
    }
}

/// The largest ploidy of a GT field's samples (bcftools' `max_used_gt_ploidy`).
fn max_used_gt_ploidy(v: &Vals<'_>, nsmpl: usize) -> usize {
    let mut max = 0;
    for i in 0..nsmpl {
        let base = i * v.n;
        let j = (0..v.n).find(|&j| v.get(base + j).end).unwrap_or(v.n);
        max = max.max(j);
        if j == v.n {
            break;
        }
    }
    max
}
