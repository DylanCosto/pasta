//! Reading values from a record into value tokens (bcftools' `filters_set_*`).

use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::bail;
use pasta::variant::{
    Header, Record, Unpack,
    encoding::{
        self, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING, INT16_VECTOR_END, INT32_MISSING,
        INT32_VECTOR_END,
    },
    utils::{self, AcSource},
    vcf,
};

use super::{
    eval::{Cache, Context},
    token::{Fraction, GtClass, Setter, Token, expand, missing, vector_end},
};

static GT_ALLELES_WARNED: AtomicBool = AtomicBool::new(false);
static AC_ALLELE_WARNED: AtomicBool = AtomicBool::new(false);

fn chrom<'h>(header: &'h Header, rec: &Record) -> &'h [u8] {
    header.contig_name(rec.rid()).unwrap_or(b"(null)")
}

fn lossy(s: &[u8]) -> String {
    String::from_utf8_lossy(s).into_owned()
}

/// Prints a warning to stderr the first time `flag` is used, as bcftools' static flags do.
fn warn_once(flag: &AtomicBool, args: std::fmt::Arguments<'_>) {
    if !flag.swap(true, Ordering::Relaxed) {
        eprint!("{args}");
    }
}

/// Sets a value token from the record.
pub(super) fn apply(
    setter: Setter,
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    match setter {
        Setter::Qual => {
            let q = rec.qual_raw();
            tok.values[0] = if encoding::float_is_missing(q) {
                missing()
            } else {
                f64::from(q)
            };
            tok.nvalues = 1;
        }
        Setter::Type => {
            let t = rec.variant_types().bits();
            // HTSlib's VCF_REF is 0; bcftools shifts the bits so that "ref" is 1.
            tok.values[0] = if t == 0 { 1.0 } else { f64::from(t << 1) };
            tok.nvalues = 1;
        }
        Setter::Info => set_info(rec, tok),
        Setter::InfoInt => set_info_int(ctx, header, rec, tok),
        Setter::InfoFloat => set_info_float(ctx, header, rec, tok),
        Setter::InfoString => set_info_string(ctx, header, rec, tok),
        Setter::InfoFlag => {
            let present = rec.info_fields().iter().any(|f| f.key() == tok.hdr_id);
            tok.values[0] = if present { 1.0 } else { 0.0 };
            tok.nvalues = 1;
        }
        Setter::FormatInt => set_format_int(ctx, header, rec, tok)?,
        Setter::FormatFloat => set_format_float(ctx, header, rec, tok)?,
        Setter::FormatString => set_format_string(ctx, header, rec, tok)?,
        Setter::Genotype(class) => set_genotype(header, rec, tok, class)?,
        Setter::GenotypeString => set_genotype_string(header, rec, tok)?,
        Setter::Ilen => set_ilen(rec, tok),
        Setter::RefString => {
            tok.str_value.clear();
            if rec.allele_count() > 0 {
                tok.str_value.extend_from_slice(rec.allele(0));
            }
            tok.nvalues = tok.str_value.len();
        }
        Setter::AltString => set_alt_string(rec, tok),
        Setter::Missing(fraction) => set_missing_count(header, rec, tok, fraction),
        Setter::NAlt => {
            tok.values[0] = rec.allele_count() as f64 - 1.0;
            tok.nvalues = 1;
        }
        Setter::Chrom => {
            tok.str_value.clear();
            tok.str_value.extend_from_slice(chrom(header, rec));
            tok.nvalues = tok.str_value.len();
            tok.is_str = true;
        }
        Setter::Pos => {
            tok.values[0] = (rec.pos() + 1) as f64;
            tok.nvalues = 1;
        }
        Setter::Ac | Setter::Af => set_ac_af(setter, ctx, header, rec, tok)?,
        Setter::An => {
            tok.nvalues = 1;
            tok.values[0] = if cache_ac_stats(ctx, header, rec)? {
                f64::from(ctx.site.an)
            } else {
                0.0
            };
        }
        Setter::Mac => set_mac(ctx, header, rec, tok)?,
        Setter::Maf => set_maf(ctx, header, rec, tok)?,
    }
    Ok(())
}

/// `ILEN`: the length change of each ALT, missing for symbolic alleles.
#[inline]
fn set_ilen(rec: &Record, tok: &mut Token) {
    if rec.allele_count() == 0 {
        tok.nvalues = 0;
        return;
    }
    let n = rec.allele_count() - 1;
    tok.nvalues = n;
    expand(&mut tok.values, n);
    let rlen = rec.allele(0).len() as f64;
    for i in 1..rec.allele_count() {
        let a = rec.allele(i);
        tok.values[i - 1] = if a.first() == Some(&b'<') {
            missing()
        } else {
            a.len() as f64 - rlen
        };
    }
}

/// The computed `AC` or `AF` of each ALT, or of the subscripted one.
#[inline]
fn set_ac_af(
    setter: Setter,
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    let ok = cache_ac_stats(ctx, header, rec)? && ctx.site.an != 0;
    if !ok {
        tok.nvalues = 0;
        return Ok(());
    }
    let n_allele = rec.allele_count();
    let get = |ctx: &Context, i: usize| {
        if setter == Setter::Ac {
            f64::from(ctx.site.ac[i])
        } else {
            ctx.site.af[i]
        }
    };
    // Here the index is relative to ALT, as for Number=A.
    if let Ok(idx) = usize::try_from(tok.idx) {
        tok.nvalues = 1;
        tok.values[0] = if idx + 1 < n_allele {
            get(ctx, idx + 1)
        } else {
            0.0
        };
    } else if n_allele == 1 {
        tok.nvalues = 1;
        tok.values[0] = 0.0;
    } else {
        expand(&mut tok.values, n_allele);
        for i in 1..n_allele {
            tok.values[i - 1] = get(ctx, i);
        }
        tok.nvalues = n_allele - 1;
    }
    Ok(())
}

/// `MAC`: the smallest non-zero allele count, 0 when one allele has them all.
#[inline]
fn set_mac(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    if !cache_ac_stats(ctx, header, rec)? || ctx.site.an == 0 {
        tok.nvalues = 0;
        return Ok(());
    }
    tok.nvalues = 1;
    let min = ctx.site.ac[..rec.allele_count()]
        .iter()
        .copied()
        .filter(|&c| c != 0)
        .min()
        .unwrap_or(i32::MAX);
    if min == ctx.site.an {
        tok.values[0] = 0.0;
    } else {
        tok.values[0] = f64::from(min);
    }
    Ok(())
}

/// `MAF`: the smallest non-zero allele frequency, 0 when one allele has them all.
#[inline]
fn set_maf(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    if !cache_ac_stats(ctx, header, rec)? || ctx.site.an == 0 {
        tok.nvalues = 0;
        return Ok(());
    }
    tok.nvalues = 1;
    let min = ctx.site.af[..rec.allele_count()]
        .iter()
        .copied()
        .filter(|&f| f != 0.0)
        .fold(f64::INFINITY, f64::min);
    tok.values[0] = if min == 1.0 { 0.0 } else { min };
    Ok(())
}

/// A `Number=1` INFO tag, read from the decoded field (`filters_set_info`).
fn set_info(rec: &Record, tok: &mut Token) {
    let Some(field) = rec.info_fields().iter().find(|f| f.key() == tok.hdr_id) else {
        tok.nvalues = 0;
        tok.str_value.clear();
        return;
    };
    match field.type_code() {
        // BCF_BT_CHAR
        7 => {
            let data = rec.info_data(field);
            let n = field.len().min(data.len());
            tok.str_value.clear();
            tok.str_value.extend_from_slice(&data[..n]);
            tok.nvalues = n;
        }
        // BCF_BT_FLOAT
        5 => {
            let v = field.float_value();
            if encoding::float_is_missing(v) {
                tok.nvalues = 0;
            } else {
                tok.values[0] = f64::from(v);
                tok.nvalues = 1;
            }
            tok.str_value.clear();
        }
        ty => {
            tok.str_value.clear();
            let v = field.int_value();
            if int_is_missing(ty, v) {
                tok.nvalues = 0;
            } else {
                tok.values[0] = v as f64;
                tok.nvalues = 1;
            }
        }
    }
}

/// Whether `v`, the value of a single-value INFO field of BCF type `ty`, is missing.
#[inline]
fn int_is_missing(ty: u8, v: i64) -> bool {
    match ty {
        1 => v == i64::from(INT8_MISSING),
        2 => v == i64::from(INT16_MISSING),
        3 => v == i64::from(INT32_MISSING),
        _ => false,
    }
}

/// One value of an INFO vector (bcftools' `bcf_get_info_value`): `None` when out of range,
/// missing, or past a vector end.
fn info_value(rec: &Record, id: i32, ivec: i32) -> Option<f64> {
    let field = rec.info_fields().iter().find(|f| f.key() == id)?;
    if field.len() == 1 {
        // bcftools gives a single value at any index and without the missing check, so
        // `AC[1]` of `AC=.` is -128; vtools reads it at index 0 only, and `.` as nothing.
        if ivec != 0 {
            return None;
        }
        return match field.type_code() {
            5 => {
                let v = field.float_value();
                (!encoding::float_is_missing(v)).then(|| f64::from(v))
            }
            ty @ 1..=3 => {
                let v = field.int_value();
                (!int_is_missing(ty, v)).then_some(v as f64)
            }
            _ => Some(0.0),
        };
    }
    let ivec = usize::try_from(ivec).unwrap_or(0);
    if ivec >= field.len() {
        return None;
    }
    let values = rec.info_values(field);
    if field.type_code() == 5 {
        for j in 0..ivec {
            if encoding::float_is_vector_end(values.float(j)?) {
                return None;
            }
        }
        let v = values.float(ivec)?;
        return (!encoding::float_is_missing(v)).then(|| f64::from(v));
    }
    let raw = |j: usize| -> Option<(bool, bool, i64)> {
        let d = values.bytes();
        Some(match field.type_code() {
            1 => {
                let v = d[j] as i8;
                (v == INT8_MISSING, v == INT8_VECTOR_END, i64::from(v))
            }
            2 => {
                let v = encoding::le_i16(&d[j * 2..]);
                (v == INT16_MISSING, v == INT16_VECTOR_END, i64::from(v))
            }
            3 => {
                let v = encoding::le_i32(&d[j * 4..]);
                (v == INT32_MISSING, v == INT32_VECTOR_END, i64::from(v))
            }
            _ => return None,
        })
    };
    for j in 0..ivec {
        if raw(j)?.1 {
            return None;
        }
    }
    let (is_missing, _, v) = raw(ivec)?;
    (!is_missing).then_some(v as f64)
}

/// The last index of a vector to visit for a list subscript (`end` in `filters_set_info_*`).
fn list_end(tok: &Token, nvalues: usize) -> usize {
    let end = if tok.idxs.last().is_some_and(|&l| l < 0) {
        nvalues - 1
    } else {
        tok.idxs.len() - 1
    };
    end.min(nvalues - 1)
}

/// Whether a list subscript selects index `i`.
#[inline]
fn selects(tok: &Token, i: usize) -> bool {
    i >= tok.idxs.len() || tok.idxs[i] != 0
}

/// The dictionary ID of the token's tag in `header`: the one found when parsing while it still
/// names the tag (unless the caller passes another header), saving a lookup by name.
fn tag_id(header: &Header, tok: &Token) -> Option<u32> {
    let tag = tok.tag.as_deref()?;
    match u32::try_from(tok.hdr_id) {
        Ok(id) if header.id_name(id) == Some(tag) => Some(id),
        _ => header.id(tag),
    }
}

fn set_info_int(ctx: &mut Context, header: &Header, rec: &mut Record, tok: &mut Token) {
    if tok.idx == -2 {
        let n = tag_id(header, tok)
            .and_then(|id| rec.get_info_i32_id(header, id, &mut ctx.tmpi).ok())
            .unwrap_or(0);
        if n == 0 {
            tok.nvalues = 0;
            return;
        }
        expand(&mut tok.values, n);
        let mut j = 0;
        for i in 0..=list_end(tok, n) {
            if selects(tok, i) {
                // bcftools copies a missing value as the number -2147483648, so `AC<0` holds
                // for `AC=.`; vtools reads it as missing, as for floats.
                let v = ctx.tmpi[i];
                tok.values[j] = if v == INT32_MISSING {
                    missing()
                } else {
                    f64::from(v)
                };
                j += 1;
            }
        }
        tok.nvalues = j;
    } else {
        match info_value(rec, tok.hdr_id, tok.idx) {
            Some(v) => {
                tok.values[0] = v;
                tok.nvalues = 1;
            }
            None => tok.nvalues = 0,
        }
    }
}

fn set_info_float(ctx: &mut Context, header: &Header, rec: &mut Record, tok: &mut Token) {
    if tok.idx == -2 {
        let n = tag_id(header, tok)
            .and_then(|id| rec.get_info_f32_id(header, id, &mut ctx.tmpf).ok())
            .unwrap_or(0);
        if n == 0 {
            tok.nvalues = 0;
            return;
        }
        expand(&mut tok.values, n);
        let mut j = 0;
        for i in 0..=list_end(tok, n) {
            if selects(tok, i) {
                let v = ctx.tmpf[i];
                tok.values[j] = if encoding::float_is_missing(v) {
                    missing()
                } else {
                    f64::from(v)
                };
                j += 1;
            }
        }
        tok.nvalues = j;
    } else {
        match info_value(rec, tok.hdr_id, tok.idx) {
            Some(v) => {
                tok.values[0] = v;
                tok.nvalues = 1;
            }
            None => tok.nvalues = 0,
        }
    }
}

fn set_info_string(ctx: &mut Context, header: &Header, rec: &mut Record, tok: &mut Token) {
    let tag = tok.tag.as_deref().unwrap_or_default();
    let Ok(s) = rec.get_info_string(header, tag) else {
        tok.nvalues = 0;
        tok.str_value.clear();
        return;
    };
    let n = s.len();
    if let Ok(idx) = usize::try_from(tok.idx) {
        // The idx-th comma-separated field.
        let mut i = 0;
        let mut ss = 0;
        while ss < n && i < idx {
            if s[ss] == b',' {
                i += 1;
            }
            ss += 1;
        }
        tok.str_value.clear();
        if ss == n || i != idx {
            tok.nvalues = 0;
            return;
        }
        let se = ss + s[ss..].iter().position(|&c| c == b',').unwrap_or(n - ss);
        tok.str_value.extend_from_slice(&s[ss..se]);
    } else if tok.idx == -2 && tok.idxs.first() == Some(&-1) {
        tok.str_value.clear();
        tok.str_value.extend_from_slice(s);
    } else if tok.idx == -2 {
        // Selected fields, each followed by a comma.
        ctx.strbuf.clear();
        let s = &s[..memchr::memchr(0, s).unwrap_or(n)];
        if n > 0 {
            let iend = list_end(tok, n);
            let mut beg = 0;
            for i in 0..=iend {
                if beg > s.len() {
                    break;
                }
                let end = beg
                    + s[beg..]
                        .iter()
                        .position(|&c| c == b',')
                        .unwrap_or(s.len() - beg);
                if selects(tok, i) {
                    ctx.strbuf.extend_from_slice(&s[beg..end]);
                    ctx.strbuf.push(b',');
                }
                beg = end + 1;
            }
        }
        std::mem::swap(&mut ctx.strbuf, &mut tok.str_value);
    }
    tok.nvalues = tok.str_value.len();
}

/// Fails unless the record has as many samples as the token (bcftools checks this for every
/// FORMAT value).
fn check_samples(header: &Header, rec: &Record, tok: &Token) -> anyhow::Result<()> {
    if rec.sample_count() != tok.nsamples {
        bail!(
            "Incorrect number of FORMAT fields at {}:{} .. {}, {} vs {}",
            lossy(chrom(header, rec)),
            rec.pos() + 1,
            lossy(tok.tag.as_deref().unwrap_or_default()),
            rec.sample_count(),
            tok.nsamples
        );
    }
    Ok(())
}

/// Fills a FORMAT token from per-sample values `src` (`nsrc1` per sample), converting each
/// with `conv` (missing and vector-end aware), as `filters_set_format_int/float` do.
fn fill_format<T: Copy>(
    ctx_mask: &[u64],
    src: &[T],
    nsrc1: usize,
    tok: &mut Token,
    conv: impl Fn(T) -> f64,
) {
    let nsamples = tok.nsamples;
    tok.nval1 = if tok.idx >= 0 {
        1
    } else if tok.nuidxs > 0 {
        tok.nuidxs
    } else {
        nsrc1
    };
    tok.nvalues = tok.nval1 * nsamples;
    expand(&mut tok.values, tok.nvalues);
    let nval1 = tok.nval1;
    if let Ok(idx) = usize::try_from(tok.idx) {
        for i in 0..nsamples {
            if !tok.usmpl[i] {
                continue;
            }
            tok.values[i] = if idx >= nsrc1 {
                missing()
            } else {
                conv(src[i * nsrc1 + idx])
            };
        }
    } else if tok.idx == -3 {
        for i in 0..nsamples {
            if !tok.usmpl[i] {
                continue;
            }
            let s = &src[i * nsrc1..(i + 1) * nsrc1];
            let dst = &mut tok.values[i * nval1..(i + 1) * nval1];
            let mut j = 0;
            for (k, &v) in s.iter().enumerate() {
                if k < 64 && ctx_mask[i] & (1u64 << k) != 0 {
                    dst[j] = conv(v);
                    j += 1;
                }
            }
            if j == 0 {
                dst[0] = missing();
                j = 1;
            }
            for d in &mut dst[j..] {
                *d = vector_end();
            }
        }
    } else {
        let kend = if tok.idxs.last().is_some_and(|&l| l < 0) {
            nval1
        } else {
            tok.idxs.len()
        };
        for i in 0..nsamples {
            if !tok.usmpl[i] {
                continue;
            }
            let s = &src[i * nsrc1..(i + 1) * nsrc1];
            let base = i * nval1;
            let mut j = 0;
            for k in 0..kend {
                if k < tok.idxs.len() && tok.idxs[k] == 0 {
                    continue;
                }
                if j >= nval1 {
                    break;
                }
                // bcftools reads past the sample's values here; they are missing.
                tok.values[base + j] = s.get(k).map_or_else(missing, |&v| conv(v));
                j += 1;
            }
            if j == 0 && nval1 > 0 {
                tok.values[base] = missing();
                j = 1;
            }
            for d in &mut tok.values[base + j..base + nval1] {
                *d = vector_end();
            }
        }
    }
}

fn set_format_int(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    check_samples(header, rec, tok)?;
    let got = tag_id(header, tok).map(|id| rec.get_format_i32_id(header, id, &mut ctx.tmpi));
    let nvals = match got {
        Some(Ok(n)) if tok.nsamples > 0 => n,
        _ => {
            tok.nvalues = 0;
            return Ok(());
        }
    };
    if tok.idx == -3 && !cache_genotypes(ctx, header, rec) {
        tok.nvalues = 0;
        return Ok(());
    }
    let nsrc1 = nvals / tok.nsamples;
    let tmpi = std::mem::take(&mut ctx.tmpi);
    // With `[GT]` subscripts bcftools copies the values without the missing check (a missing
    // value is the number -2147483648, so `AD[*:GT]<0` holds for `AD=.`); vtools reads them
    // as missing there too, as for floats and the other subscripts.
    fill_format(&ctx.gt.mask, &tmpi, nsrc1, tok, |v| match v {
        INT32_MISSING => missing(),
        INT32_VECTOR_END => vector_end(),
        v => f64::from(v),
    });
    ctx.tmpi = tmpi;
    Ok(())
}

fn set_format_float(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    check_samples(header, rec, tok)?;
    let got = tag_id(header, tok).map(|id| rec.get_format_f32_id(header, id, &mut ctx.tmpf));
    let nvals = match got {
        Some(Ok(n)) if tok.nsamples > 0 => n,
        _ => {
            tok.nvalues = 0;
            return Ok(());
        }
    };
    if tok.idx == -3 && !cache_genotypes(ctx, header, rec) {
        tok.nvalues = 0;
        return Ok(());
    }
    let nsrc1 = nvals / tok.nsamples;
    let tmpf = std::mem::take(&mut ctx.tmpf);
    fill_format(&ctx.gt.mask, &tmpf, nsrc1, tok, |v| {
        if encoding::float_is_missing(v) {
            missing()
        } else if encoding::float_is_vector_end(v) {
            vector_end()
        } else {
            f64::from(v)
        }
    });
    ctx.tmpf = tmpf;
    Ok(())
}

fn set_format_string(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
) -> anyhow::Result<()> {
    check_samples(header, rec, tok)?;
    tok.nvalues = 0;
    tok.str_value.clear();
    let tag = tok.tag.as_deref().unwrap_or_default();
    let Ok((data, _)) = rec.get_format_chars(header, tag) else {
        return Ok(());
    };
    tok.str_value.extend_from_slice(data);
    if tok.idx == -3 && !cache_genotypes(ctx, header, rec) {
        tok.str_value.clear();
        return Ok(());
    }
    let nstr = tok.str_value.len();
    if tok.nsamples == 0 {
        return Ok(());
    }
    tok.nvalues = nstr;
    tok.nval1 = nstr / tok.nsamples;
    let nval1 = tok.nval1;
    if nval1 == 0 {
        return Ok(());
    }
    for i in 0..tok.nsamples {
        if !tok.usmpl[i] {
            continue;
        }
        // Keep the selected comma-separated fields of the sample's block.
        ctx.strbuf.clear();
        let mut kept = false;
        let block = &tok.str_value[i * nval1..(i + 1) * nval1];
        let text = &block[..memchr::memchr(0, block).unwrap_or(nval1)];
        for (idx, field) in text.split(|&c| c == b',').enumerate() {
            let keep = if tok.idx >= 0 {
                tok.idx as usize == idx
            } else if tok.idx == -3 {
                idx < 64 && ctx.gt.mask[i] & (1u64 << idx) != 0
            } else if idx < tok.idxs.len() {
                tok.idxs[idx] != 0
            } else {
                tok.idxs.last().is_some_and(|&l| l < 0)
            };
            if keep {
                if kept {
                    ctx.strbuf.push(b',');
                }
                kept = true;
                ctx.strbuf.extend_from_slice(field);
                if tok.idx >= 0 {
                    break;
                }
            }
        }
        if !kept {
            ctx.strbuf.push(b'.');
        }
        let block = &mut tok.str_value[i * nval1..(i + 1) * nval1];
        let n = ctx.strbuf.len().min(nval1);
        block[..n].copy_from_slice(&ctx.strbuf[..n]);
        block[n..].fill(0);
    }
    Ok(())
}

/// The GT field's index among the record's FORMAT fields.
fn gt_index(header: &Header, rec: &Record) -> Option<usize> {
    let id = header.id(b"GT")? as i32;
    rec.format_fields().iter().position(|f| f.id() == id)
}

/// The shape of a sample's genotype: how many alleles were read (up to a missing one or the
/// vector end), and whether one was missing, one was REF, or two differed.
#[derive(Clone, Copy, Debug, Default)]
struct GtShape {
    n: usize,
    missing: bool,
    has_ref: bool,
    is_het: bool,
}

/// Reads the shape of one sample's GT values.
#[inline]
fn gt_shape(values: impl Iterator<Item = i32>) -> GtShape {
    let mut shape = GtShape::default();
    let mut jal = 0;
    for ial in values {
        if ial == INT32_VECTOR_END {
            break;
        }
        if ial >> 1 == 0 {
            shape.missing = true;
            break;
        }
        if ial >> 1 == 1 {
            shape.has_ref = true;
        }
        if shape.n > 0 && (ial >> 1) != (jal >> 1) {
            shape.is_het = true;
        }
        jal = ial;
        shape.n += 1;
    }
    shape
}

/// The class string of a genotype shape, NUL-padded.
#[inline]
fn gt_class_value(class: GtClass, g: GtShape) -> [u8; 4] {
    match class {
        GtClass::Presence => {
            if g.n == 0 || g.missing {
                *b"mis\0"
            } else if !g.has_ref || g.is_het {
                *b"alt\0"
            } else {
                *b"ref\0"
            }
        }
        _ if g.n == 0 || g.missing => *b".\0\0\0",
        GtClass::Zygosity => {
            if g.n == 1 {
                *b"hap\0"
            } else if !g.is_het {
                *b"hom\0"
            } else {
                *b"het\0"
            }
        }
        GtClass::Alleles => match (g.n == 1, g.is_het, g.has_ref) {
            (true, _, true) => *b"r\0\0\0",
            (true, _, false) => *b"a\0\0\0",
            (false, false, true) => *b"rr\0\0",
            (false, false, false) => *b"aa\0\0",
            (false, true, true) => *b"ra\0\0",
            (false, true, false) => *b"aA\0\0",
        },
    }
}

/// Writes each sample's class into `N`-byte blocks of `out`.
fn fill_gt_classes<const N: usize>(
    out: &mut [u8],
    data: &[u8],
    size: usize,
    n: usize,
    ty: u8,
    nsample: usize,
    class: GtClass,
) {
    for (i, dst) in out
        .as_chunks_mut::<N>()
        .0
        .iter_mut()
        .take(nsample)
        .enumerate()
    {
        let p = data.get(i * size..(i + 1) * size).unwrap_or_default();
        let shape = match ty {
            1 => gt_shape(p.get(..n).unwrap_or_default().iter().map(|&b| {
                let v = b as i8;
                if v == INT8_VECTOR_END {
                    INT32_VECTOR_END
                } else {
                    i32::from(v)
                }
            })),
            2 => gt_shape(
                p.get(..2 * n)
                    .unwrap_or_default()
                    .as_chunks::<2>()
                    .0
                    .iter()
                    .map(|&c| {
                        let v = i16::from_le_bytes(c);
                        if v == INT16_VECTOR_END {
                            INT32_VECTOR_END
                        } else {
                            i32::from(v)
                        }
                    }),
            ),
            _ => gt_shape(
                p.get(..4 * n)
                    .unwrap_or_default()
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|&c| i32::from_le_bytes(c)),
            ),
        };
        for (d, v) in dst.iter_mut().zip(gt_class_value(class, shape)) {
            *d = v;
        }
    }
}

/// Classifies each sample's genotype as a short string (`_filters_set_genotype`).
fn set_genotype(
    header: &Header,
    rec: &mut Record,
    tok: &mut Token,
    class: GtClass,
) -> anyhow::Result<()> {
    tok.nvalues = 0;
    tok.str_value.clear();
    rec.unpack(Unpack::FMT);
    let Some(k) = gt_index(header, rec) else {
        return Ok(());
    };
    let fmt = &rec.format_fields()[k];
    let nsmpl = header.sample_count();
    let nvals1 = if class == GtClass::Alleles { 3 } else { 4 };
    tok.str_value.resize(nvals1 * nsmpl, 0);
    let ty = fmt.type_code();
    if !(1..=3).contains(&ty) {
        bail!(
            "The GT type is not lineognised: {ty} at {}:{}",
            lossy(chrom(header, rec)),
            rec.pos() + 1
        );
    }
    let nsample = rec.sample_count().min(nsmpl);
    let data = rec.format_data(fmt);
    let (size, n) = (fmt.size(), fmt.n());
    if nvals1 == 3 {
        fill_gt_classes::<3>(&mut tok.str_value, data, size, n, ty, nsample, class);
    } else {
        fill_gt_classes::<4>(&mut tok.str_value, data, size, n, ty, nsample, class);
    }
    tok.nvalues = nvals1 * nsmpl;
    tok.nval1 = nvals1;
    Ok(())
}

/// Formats each sample's genotype into fixed-size NUL-padded blocks
/// (`filters_set_genotype_string`).
fn set_genotype_string(header: &Header, rec: &mut Record, tok: &mut Token) -> anyhow::Result<()> {
    tok.nvalues = 0;
    tok.str_value.clear();
    rec.unpack(Unpack::FMT);
    let Some(k) = gt_index(header, rec) else {
        return Ok(());
    };
    let fmt = &rec.format_fields()[k];
    let nsmpl = rec.sample_count();
    let mut blen = 4;
    'retry: loop {
        tok.str_value.clear();
        for i in 0..nsmpl {
            let start = tok.str_value.len();
            let data = rec.format_sample_data(fmt, i);
            if vcf::format_genotype(&mut tok.str_value, 0, fmt.type_code(), fmt.n(), data).is_err()
            {
                bail!(
                    "The GT type is not lineognised: {} at {}:{}",
                    fmt.type_code(),
                    lossy(chrom(header, rec)),
                    rec.pos() + 1
                );
            }
            tok.str_value.push(0);
            let len = tok.str_value.len() - start;
            if len > blen {
                // Too many alleles or too large a ploidy for the block.
                blen *= 2;
                continue 'retry;
            }
            tok.str_value.resize(start + blen, 0);
        }
        break;
    }
    tok.nvalues = tok.str_value.len();
    tok.nval1 = blen;
    Ok(())
}

fn set_alt_string(rec: &Record, tok: &mut Token) {
    tok.str_value.clear();
    let n_allele = rec.allele_count();
    if let Ok(idx) = usize::try_from(tok.idx) {
        if n_allele > idx + 1 {
            tok.str_value.extend_from_slice(rec.allele(idx + 1));
        } else {
            tok.str_value.push(b'.');
        }
    } else if tok.idx == -2 {
        if n_allele >= 2 {
            let end = if tok.idxs.last().is_some_and(|&l| l < 0) {
                n_allele - 2
            } else {
                (tok.idxs.len() - 1).min(n_allele - 2)
            };
            for i in 0..=end {
                if selects(tok, i) {
                    if !tok.str_value.is_empty() {
                        tok.str_value.push(b',');
                    }
                    tok.str_value.extend_from_slice(rec.allele(i + 1));
                }
            }
        }
    } else if n_allele > 1 {
        for i in 1..n_allele {
            if i > 1 {
                tok.str_value.push(b',');
            }
            tok.str_value.extend_from_slice(rec.allele(i));
        }
    } else if n_allele == 1 {
        tok.str_value.push(b'.');
    }
    tok.nvalues = tok.str_value.len();
}

/// `N_MISSING` and `F_MISSING` (`filters_set_nmissing`).
fn set_missing_count(header: &Header, rec: &mut Record, tok: &mut Token, fraction: Fraction) {
    rec.unpack(Unpack::FMT);
    if rec.sample_count() == 0 {
        tok.nvalues = 1;
        tok.values[0] = 0.0;
        return;
    }
    let Some(k) = gt_index(header, rec) else {
        tok.nvalues = 0;
        return;
    };
    let fmt = &rec.format_fields()[k];
    let data = rec.format_data(fmt);
    let (size, n) = (fmt.size(), fmt.n());
    let nsample = rec.sample_count();
    let all_used = tok
        .usmpl
        .get(..nsample)
        .is_some_and(|u| u.iter().all(|&u| u));
    let (nsmpl, nmissing) = match fmt.type_code() {
        // Diploid 8-bit genotypes of all samples, the common case: count without branches.
        1 if size == 2 && n == 2 && all_used => (
            nsample as u32,
            count_missing_diploid(data.get(..2 * nsample).unwrap_or(data)),
        ),
        1 => count_missing::<1>(data, size, n, nsample, tok, |b| {
            let v = i32::from(b[0] as i8);
            (v == i32::from(INT8_VECTOR_END), v)
        }),
        2 => count_missing::<2>(data, size, n, nsample, tok, |b| {
            let v = i16::from_le_bytes(b);
            (v == INT16_VECTOR_END, i32::from(v))
        }),
        3 => count_missing::<4>(data, size, n, nsample, tok, |b| {
            let v = i32::from_le_bytes(b);
            (v == INT32_VECTOR_END, v)
        }),
        _ => (0, 0),
    };
    tok.nsamples = 0;
    tok.nvalues = 1;
    tok.values[0] = match fraction {
        Fraction::Count => f64::from(nmissing),
        Fraction::Fraction if nsmpl > 0 => f64::from(nmissing) / f64::from(nsmpl),
        Fraction::Fraction => 0.0,
    };
}

/// Counts the samples of diploid 8-bit genotypes whose first allele is missing, or whose second
/// is missing and follows a first that is not the vector end.
#[inline]
fn count_missing_diploid(data: &[u8]) -> u32 {
    // An allele is missing when `value >> 1 == 0`, i.e. the byte is 0 or 1.
    const VECTOR_END: u8 = INT8_VECTOR_END as u8;
    data.as_chunks::<2>()
        .0
        .iter()
        .map(|&[a, b]| u32::from(a <= 1 || (a != VECTOR_END && b <= 1)))
        .sum()
}

/// Counts the used samples and those with a missing allele before the vector end, reading
/// `W`-byte GT values with `read` (which returns whether a value is the vector end, and the
/// value).
#[inline]
fn count_missing<const W: usize>(
    data: &[u8],
    size: usize,
    n: usize,
    nsample: usize,
    tok: &Token,
    read: impl Fn([u8; W]) -> (bool, i32),
) -> (u32, u32) {
    let (mut nsmpl, mut nmissing) = (0u32, 0u32);
    if size == 0 {
        return (0, 0);
    }
    for (i, p) in data.chunks_exact(size).take(nsample).enumerate() {
        if !tok.uses(i) {
            continue;
        }
        nsmpl += 1;
        let values = p.get(..n * W).unwrap_or_default().as_chunks::<W>().0;
        for &v in values {
            let (end, v) = read(v);
            if end {
                break;
            }
            if v >> 1 == 0 {
                nmissing += 1;
                break;
            }
        }
    }
    (nsmpl, nmissing)
}

/// Caches each sample's GT alleles as a bit mask (`filters_cache_genotypes`); `false` when
/// there are no genotypes or an allele is above 63.
pub(super) fn cache_genotypes(ctx: &mut Context, header: &Header, rec: &mut Record) -> bool {
    match ctx.gt.state {
        Cache::Valid => return true,
        Cache::Failed => return false,
        Cache::Stale => {}
    }
    ctx.gt.state = Cache::Failed;
    let nbuf = rec.get_genotypes(header, &mut ctx.gt.buf).unwrap_or(0);
    let nsample = rec.sample_count();
    if nbuf == 0 || nsample == 0 {
        return false;
    }
    let ngt1 = nbuf / nsample;
    ctx.gt.mask.clear();
    ctx.gt.mask.resize(nsample, 0);
    for i in 0..nsample {
        for &v in &ctx.gt.buf[i * ngt1..(i + 1) * ngt1] {
            if v >> 1 == 0 {
                continue;
            }
            if v == INT32_VECTOR_END {
                break;
            }
            let allele = (v >> 1) - 1;
            if !(0..=63).contains(&allele) {
                warn_once(
                    &GT_ALLELES_WARNED,
                    format_args!(
                        "Too many alleles, skipping GT filtering at this site {}:{}. (This warning is printed only once.)\n",
                        lossy(chrom(header, rec)),
                        rec.pos() + 1
                    ),
                );
                return false;
            }
            ctx.gt.mask[i] |= 1u64 << allele;
        }
    }
    ctx.gt.state = Cache::Valid;
    true
}

/// Caches allele counts and frequencies, from INFO/AC and AN when all samples are used, else
/// from the genotypes of the used samples (`filters_cache_AC_stats`). `false` on failure.
fn cache_ac_stats(ctx: &mut Context, header: &Header, rec: &mut Record) -> anyhow::Result<bool> {
    match ctx.site.state {
        Cache::Valid => return Ok(true),
        Cache::Failed => return Ok(false),
        Cache::Stale => {}
    }
    ctx.site.state = Cache::Failed;
    let n_allele = rec.allele_count();
    let site = &mut ctx.site;
    if ctx.sample_mask_set == 0 && utils::calc_ac(header, rec, &mut site.ac, AcSource::Info)? {
        site.n_samples = rec.sample_count();
        site.n_missing = -1;
        site.an = site.ac[..n_allele].iter().sum();
        site.af.resize(site.ac.len(), 0.0);
        for i in 0..n_allele {
            site.af[i] = if site.an != 0 {
                f64::from(site.ac[i]) / f64::from(site.an)
            } else {
                0.0
            };
        }
        site.state = Cache::Valid;
        return Ok(true);
    }
    site.n_samples = if ctx.sample_mask_set > 0 {
        ctx.sample_mask_set
    } else {
        rec.sample_count()
    };
    site.n_missing = 0;
    site.an = 0;
    site.ac.clear();
    site.ac.resize(n_allele.max(1), 0);
    site.af.resize(n_allele.max(1), 0.0);
    let ngt = rec.get_genotypes(header, &mut ctx.tmpi).unwrap_or(0);
    let nsample = rec.sample_count();
    let ngt1 = ngt.checked_div(nsample).unwrap_or(0);
    for i in 0..nsample {
        if ctx.sample_mask_set > 0 && !ctx.sample_mask[i] {
            continue;
        }
        let mut is_missing = true;
        for &v in &ctx.tmpi[i * ngt1..(i + 1) * ngt1] {
            if v >> 1 == 0 {
                continue;
            }
            if v == INT32_VECTOR_END {
                break;
            }
            let ial = (v >> 1) - 1;
            let Some(slot) = usize::try_from(ial).ok().filter(|&a| a < n_allele) else {
                warn_once(
                    &AC_ALLELE_WARNED,
                    format_args!(
                        "The allele index too large ({ial}), skipping GT parsing at this site {}:{}. (This warning is printed only once.)\n",
                        lossy(chrom(header, rec)),
                        rec.pos() + 1
                    ),
                );
                site.n_samples = 0;
                return Ok(false);
            };
            site.ac[slot] += 1;
            site.an += 1;
            is_missing = false;
        }
        if is_missing {
            site.n_missing += 1;
        }
    }
    for i in 0..n_allele {
        site.af[i] = if site.an != 0 {
            f64::from(site.ac[i]) / f64::from(site.an)
        } else {
            0.0
        };
    }
    site.state = Cache::Valid;
    Ok(true)
}
