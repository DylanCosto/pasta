//! The site fields of an output line: CHROM to QUAL, FILTER, and INFO (bcftools merge's
//! `merge_chrom2qual`, `merge_filter`, `merge_info`, and `merge_AGR_info_tag`).

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Number, Record, ValueType, VariantType,
    encoding::{FLOAT_MISSING, FLOAT_VECTOR_END, INT32_MISSING},
    record::{InfoField, TypedValues},
};

use super::{FilterLogic, Merger, alleles};

/// The buffer of one merged Number=A/G/R INFO tag (`AGR_info_t`).
#[derive(Debug, Default)]
struct AgrInfo {
    /// The tag's output dictionary ID.
    out_id: u32,
    /// The tag's type in the header of the file it was first seen in.
    ty: Option<ValueType>,
    nvals: usize,
    /// Whether the buffer holds this site's initial values (`nbuf != 0`).
    started: bool,
    /// Integer or float bits, one per value.
    words: Vec<u32>,
    /// Comma-separated strings.
    text: Vec<u8>,
}

impl AgrInfo {
    /// Sets the merged values on the output line.
    #[inline]
    fn write(
        &self,
        out_hdr: &Header,
        out: &mut Record,
        scratch_i: &mut Vec<i32>,
        scratch_f: &mut Vec<f32>,
    ) {
        let words = self.words.get(..self.nvals).unwrap_or_default();
        match self.ty {
            Some(ValueType::Integer) => {
                scratch_i.clear();
                scratch_i.extend(words.iter().map(|&w| w as i32));
                let _ = out.update_info_i32_id(out_hdr, self.out_id, scratch_i);
            }
            Some(ValueType::Float) => {
                scratch_f.clear();
                scratch_f.extend(words.iter().map(|&w| f32::from_bits(w)));
                let _ = out.update_info_f32_id(out_hdr, self.out_id, scratch_f);
            }
            Some(ValueType::String) => {
                // No values (Number=A without ALTs) removes the tag.
                let name = out_hdr.id_name(self.out_id).unwrap_or_default();
                let value = (self.nvals > 0).then_some(&self.text[..]);
                let _ = out.update_info_string(out_hdr, name, value);
            }
            _ => {}
        }
    }
}

/// How `merge_info` treats a reader's INFO tag, worked out once per dictionary ID.
#[derive(Clone, Copy, Debug)]
struct InfoKey {
    /// AC or AN, recomputed from the genotypes.
    skip: bool,
    out_id: Option<u32>,
    number: Option<Number>,
    /// The index of the tag's `-i` rule.
    rule: Option<usize>,
}

/// Scratch space for the site fields.
#[derive(Debug, Default)]
pub(super) struct SiteState {
    /// Per reader, by dictionary ID.
    keys: Vec<Vec<Option<InfoKey>>>,
    agr: Vec<AgrInfo>,
    nagr: usize,
    /// The merged IDs, `;`-separated, and the IDs seen.
    ids: Vec<u8>,
    id_seen: Vec<std::ops::Range<usize>>,
    /// The merged alleles that some staged record uses.
    al_idxs: Vec<bool>,
    filters: Vec<i32>,
    agr_map: Vec<usize>,
    /// The last contig looked up: reader, its contig ID, and the output's.
    rid_cache: Option<(usize, i32, i32)>,
}

impl Merger {
    /// Sets CHROM, POS, ID, QUAL, and the alleles of the output line from the staged records
    /// (bcftools' `merge_chrom2qual`).
    pub(super) fn merge_chrom2qual(&mut self) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            out,
            out_hdr,
            site,
            ..
        } = self;
        site.ids.clear();
        site.id_seen.clear();
        site.al_idxs.clear();
        site.al_idxs.resize(ma.als.len(), false);
        let mut qual: Option<f32> = None;
        let mut have_pos = false;
        for i in 0..ma.n {
            let Some(line) = ma.line(inp, i) else {
                continue;
            };
            let map = ma.cur_map(i);
            for &m in &map[1..line.allele_count()] {
                site.al_idxs[m] = true;
            }
            if !have_pos {
                let rid = match site.rid_cache {
                    Some((ci, crid, orid)) if ci == i && crid == line.rid() => orid,
                    _ => {
                        let chr = inp.header(i).contig_name(line.rid()).unwrap_or_default();
                        let Some(rid) = out_hdr.contig_id(chr) else {
                            bail!("Uh");
                        };
                        site.rid_cache = Some((i, line.rid(), rid as i32));
                        rid as i32
                    }
                };
                out.set_rid(rid);
                out.set_pos(line.pos());
                have_pos = true;
            }
            let id = line.id();
            if id != b"." {
                let seen = site.id_seen.iter().any(|r| &site.ids[r.clone()] == id);
                if !seen {
                    if !site.ids.is_empty() {
                        site.ids.push(b';');
                    }
                    let start = site.ids.len();
                    site.ids.extend_from_slice(id);
                    site.id_seen.push(start..site.ids.len());
                }
            }
            if let Some(q) = line.qual()
                && qual.is_none_or(|m| m < q)
            {
                qual = Some(q);
            }
        }
        out.set_qual(qual);
        out.update_id(Some(if site.ids.is_empty() {
            b"."
        } else {
            &site.ids[..]
        }));

        // Renumber the alleles: some merged ones may be unused by the staged records.
        let mut nout_als = 0;
        for i in 1..ma.als.len() {
            if !site.al_idxs[i] {
                continue;
            }
            nout_als += 1;
            for ir in 0..ma.n {
                ma.buf[ir].unkn_allele = 0;
                let Some(cur) = ma.buf[ir].cur else {
                    continue;
                };
                let src = ma.buf[ir].src;
                let line = inp.line_mut(src, ir, cur);
                for j in 1..line.allele_count() {
                    let rec = &mut ma.buf[ir].rec[cur];
                    if rec.map[j] == i {
                        rec.map[j] = nout_als;
                    }
                    if line.has_variant_type(j, VariantType::empty()) == Some(true)
                        && line.allele(j).first() == Some(&b'<')
                    {
                        ma.buf[ir].unkn_allele = j;
                    }
                }
            }
        }
        ma.out_als.clear();
        for i in 0..ma.als.len() {
            if i == 0 || site.al_idxs[i] {
                ma.out_als.push(ma.als.get(i));
            }
        }
        if ma.out_als.len() != nout_als + 1 {
            bail!(
                "Error: could not merge alleles at {}:{}, sanity check failed: {}!={}",
                String::from_utf8_lossy(out.chrom_safe(out_hdr)),
                out.pos() + 1,
                ma.out_als.len(),
                nout_als + 1
            );
        }
        alleles::normalize_alleles(&mut ma.out_als);
        out.update_alleles(out_hdr, ma.out_als.as_slice());
        Ok(())
    }

    /// Merges the FILTERs of the staged records (bcftools' `merge_filter`).
    pub(super) fn merge_filter(&mut self) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            out,
            out_hdr,
            site,
            opts,
            ids,
            ..
        } = self;
        if opts.filter_logic == FilterLogic::Remove {
            let any_pass =
                (0..ma.n).any(|i| ma.line(inp, i).is_some_and(|l| has_pass(l, inp.header(i))));
            if any_pass {
                let pass = out_hdr.id(b"PASS").map_or(0, |id| id as i32);
                out.add_filter(pass);
                return Ok(());
            }
        }
        site.filters.clear();
        for i in 0..ma.n {
            let Some(line) = ma.line(inp, i) else {
                continue;
            };
            let hdr = inp.header(i);
            for &f in line.filters() {
                let Some(id) = ids.out_id(i, f, hdr, out_hdr) else {
                    bail!(
                        "Error: The filter is not defined in the header: {}",
                        String::from_utf8_lossy(hdr.id_name(f as u32).unwrap_or_default())
                    );
                };
                let id = id as i32;
                if !site.filters.contains(&id) {
                    site.filters.push(id);
                }
            }
        }
        // PASS is not mixed with other filters.
        if site.filters.len() > 1
            && let Some(pass) = out_hdr.id(b"PASS")
            && let Some(k) = site.filters.iter().position(|&f| f == pass as i32)
        {
            site.filters.remove(k);
        }
        out.update_filter(&site.filters);
        Ok(())
    }

    /// Merges the INFO fields of the staged records: the first value of each tag, values
    /// combined by the `-i` rules, and Number=A/G/R values remapped to the output alleles
    /// (bcftools' `merge_info`).
    pub(super) fn merge_info(&mut self) -> anyhow::Result<()> {
        let Self {
            inp,
            ma,
            out,
            out_hdr,
            site,
            opts,
            ids,
            seen,
            rules,
            scratch_i,
            scratch_f,
            ..
        } = self;
        seen.clear();
        site.nagr = 0;
        for r in rules.iter_mut() {
            r.reset();
        }
        let nout_als = out.allele_count();
        for i in 0..ma.n {
            let Some(cur) = ma.buf[i].cur else {
                continue;
            };
            let src = ma.buf[i].src;
            let nfields = inp.line(src, i, cur).info_fields().len();
            if site.keys.len() < ma.n {
                site.keys.resize_with(ma.n, Vec::new);
            }
            for j in 0..nfields {
                let line = inp.line(src, i, cur);
                let inf = &line.info_fields()[j];
                if inf.is_removed() {
                    continue;
                }
                let hdr = inp.sr.header(i);
                let keys = &mut site.keys[i];
                let k = inf.key() as usize;
                if k >= keys.len() {
                    keys.resize(k + 1, None);
                }
                let info_key = match keys[k] {
                    Some(c) => c,
                    None => {
                        let name = hdr.id_name(inf.key() as u32).unwrap_or_default();
                        let c = InfoKey {
                            // AC and AN are recomputed from the genotypes.
                            skip: !opts.keep_ac_an && (name == b"AC" || name == b"AN"),
                            out_id: ids.out_id(i, inf.key(), hdr, out_hdr),
                            number: hdr.id_number(HeaderRecordType::Info, inf.key() as u32),
                            rule: rules.binary_search_by(|r| r.tag.as_slice().cmp(name)).ok(),
                        };
                        keys[k] = Some(c);
                        c
                    }
                };
                if info_key.skip {
                    continue;
                }
                let key = || hdr.id_name(inf.key() as u32).unwrap_or_default();
                let Some(out_id) = info_key.out_id else {
                    bail!(
                        "Error: The INFO field is not defined in the header: {}",
                        String::from_utf8_lossy(key())
                    );
                };
                let number = info_key.number;
                let is_agr = matches!(number, Some(Number::A | Number::G | Number::R));
                if let Some(k) = info_key.rule {
                    let map = &ma.buf[i].rec[cur].map;
                    // `bcf_get_info_values` reads the first field with the key, so a key
                    // given twice (`DP=1;DP=2`) adds the first value twice.
                    let first = line.info_fields()[..j]
                        .iter()
                        .find(|f| f.key() == inf.key())
                        .unwrap_or(inf);
                    rules[k].add_values(
                        hdr,
                        line,
                        first,
                        if is_agr { number } else { None },
                        map,
                        nout_als,
                        &mut site.agr_map,
                        scratch_i,
                        scratch_f,
                    )?;
                    continue;
                }
                let state = seen.get(out_id);
                if is_agr {
                    let idx = match state {
                        None => {
                            let nvals = match number {
                                Some(Number::A) => nout_als - 1,
                                Some(Number::G) => nout_als * (nout_als + 1) / 2,
                                _ => nout_als,
                            };
                            if site.nagr == site.agr.len() {
                                site.agr.push(AgrInfo::default());
                            }
                            let a = &mut site.agr[site.nagr];
                            a.out_id = out_id;
                            a.ty = hdr.id_type(HeaderRecordType::Info, inf.key() as u32);
                            a.nvals = nvals;
                            a.started = false;
                            seen.set(out_id, site.nagr as i64);
                            site.nagr += 1;
                            site.nagr - 1
                        }
                        // A tag copied first (Number=. in that file's header) is marked 0 when
                        // it was the first INFO copied; bcftools then merges into the first
                        // Number=A/G/R buffer, which is not written when there is none.
                        Some(0) if site.nagr == 0 => continue,
                        Some(v) => match usize::try_from(v) {
                            Ok(v) if v < site.nagr => v,
                            _ => bail!(
                                "Error occurred while processing INFO tag \"{}\" at {}:{}",
                                String::from_utf8_lossy(key()),
                                String::from_utf8_lossy(line.chrom_safe(hdr)),
                                line.pos() + 1
                            ),
                        },
                    };
                    let map = &ma.buf[i].rec[cur].map;
                    merge_agr_info_tag(line, hdr, inf, number, map, &mut site.agr[idx], key())?;
                    continue;
                }
                if state.is_none() {
                    out.push_info_copy(out_hdr, line, inf, out_id);
                    // An arbitrary negative value marks a copied tag.
                    seen.set(out_id, -(out.info_count() as i64 - 1));
                }
            }
        }
        for r in rules.iter_mut() {
            r.apply(out_hdr, out);
        }
        for a in &site.agr[..site.nagr] {
            a.write(out_hdr, out, scratch_i, scratch_f);
        }
        Ok(())
    }
}

/// Whether a record passes, as `bcf_has_filter(hdr, line, "PASS")` is true: no FILTER or
/// PASS among them, or PASS undefined in the header (HTSlib's -1).
fn has_pass(line: &Record, hdr: &Header) -> bool {
    let Some(pass) = hdr
        .id(b"PASS")
        .filter(|&id| hdr.id_exists(HeaderRecordType::Filter, id))
    else {
        return true;
    };
    (pass == 0 && line.filters().is_empty()) || line.filters().contains(&(pass as i32))
}

/// Copies a record's Number=A/G/R INFO values into the merged buffer at the output allele
/// indexes (bcftools' `merge_AGR_info_tag`).
#[expect(
    clippy::needless_range_loop,
    reason = "a port of merge_AGR_info_tag, indexing values by allele pairs"
)]
fn merge_agr_info_tag(
    line: &Record,
    hdr: &Header,
    inf: &InfoField,
    number: Option<Number>,
    map: &[usize],
    agr: &mut AgrInfo,
    key: &[u8],
) -> anyhow::Result<()> {
    let values = line.info_values(inf);
    let is_numeric = matches!(
        values,
        TypedValues::Int8(_)
            | TypedValues::Int16(_)
            | TypedValues::Int32(_)
            | TypedValues::Float(_)
    );
    let is_float = matches!(values, TypedValues::Float(_));
    let chrom = || String::from_utf8_lossy(line.chrom_safe(hdr)).into_owned();
    if !agr.started {
        if is_numeric {
            agr.words.clear();
            agr.words.resize(
                agr.nvals,
                if is_float {
                    FLOAT_MISSING
                } else {
                    INT32_MISSING as u32
                },
            );
            agr.started = agr.nvals > 0;
        } else if matches!(values, TypedValues::Char(_)) {
            agr.text.clear();
            agr.text.push(b'.');
            for _ in 1..agr.nvals {
                agr.text.extend_from_slice(b",.");
            }
            agr.started = true;
        } else {
            bail!(
                "Not ready for type [{}]: {} at {}",
                inf.type_code(),
                String::from_utf8_lossy(key),
                line.pos() + 1
            );
        }
    }
    let n_allele = line.allele_count();
    if is_numeric {
        // The raw value k and whether it is missing or the vector end.
        let get = |k: usize| -> Option<(u32, bool, bool)> {
            Some(match values {
                TypedValues::Float(d) => {
                    let b = u32::from_le_bytes(d.get(k * 4..k * 4 + 4)?.try_into().ok()?);
                    (b, b == FLOAT_MISSING, b == FLOAT_VECTOR_END)
                }
                _ => {
                    let raw = raw_int(&values, k)?;
                    (raw.0 as u32, raw.1, raw.2)
                }
            })
        };
        match number {
            Some(Number::A | Number::R) => {
                let ifrom = usize::from(number == Some(Number::A));
                for (k, iori) in (ifrom..n_allele).enumerate() {
                    let Some((v, missing, end)) = get(k) else {
                        break;
                    };
                    if end {
                        break;
                    }
                    if missing {
                        continue;
                    }
                    // A buffer shared with a tag of another Number may be shorter
                    // (bcftools writes past it).
                    if let Some(w) = agr.words.get_mut(map[iori] - ifrom) {
                        *w = v;
                    }
                }
            }
            _ => {
                'outer: for iori in 0..n_allele {
                    let inew = map[iori];
                    for jori in 0..=iori {
                        let jnew = map[jori];
                        let kori = iori * (iori + 1) / 2 + jori;
                        let Some((v, missing, end)) = get(kori) else {
                            break 'outer;
                        };
                        if end {
                            break 'outer;
                        }
                        if missing {
                            continue;
                        }
                        if let Some(w) = agr.words.get_mut(alleles::alleles2gt(inew, jnew)) {
                            *w = v;
                        }
                    }
                }
            }
        }
    } else {
        let src = line.info_data(inf);
        let err = || {
            anyhow::anyhow!(
                "Error at {}:{}: wrong number of fields in {}?",
                chrom(),
                line.pos() + 1,
                String::from_utf8_lossy(key)
            )
        };
        match number {
            Some(Number::A | Number::R) => {
                let ifrom = usize::from(number == Some(Number::A));
                for iori in ifrom..n_allele {
                    if copy_string_field(src, iori - ifrom, &mut agr.text, map[iori] - ifrom) != 0 {
                        return Err(err());
                    }
                }
            }
            _ => {
                for iori in 0..n_allele {
                    let inew = map[iori];
                    for jori in 0..=iori {
                        let jnew = map[jori];
                        let kori = iori * (iori + 1) / 2 + jori;
                        let knew = alleles::alleles2gt(inew, jnew);
                        if copy_string_field(src, kori, &mut agr.text, knew) != 0 {
                            return Err(err());
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Reads integer `k` of typed values as HTSlib's typed macros do: the raw value sign-extended,
/// and whether it is the type's missing value or vector end.
#[inline]
pub(super) fn raw_int(values: &TypedValues<'_>, k: usize) -> Option<(i32, bool, bool)> {
    use pasta::variant::encoding::{
        INT8_MISSING, INT8_VECTOR_END, INT16_MISSING, INT16_VECTOR_END, INT32_VECTOR_END,
    };
    Some(match *values {
        TypedValues::Int8(d) => {
            let v = *d.get(k)? as i8;
            (i32::from(v), v == INT8_MISSING, v == INT8_VECTOR_END)
        }
        TypedValues::Int16(d) => {
            let v = i16::from_le_bytes(d.get(k * 2..k * 2 + 2)?.try_into().ok()?);
            (i32::from(v), v == INT16_MISSING, v == INT16_VECTOR_END)
        }
        TypedValues::Int32(d) => {
            let v = i32::from_le_bytes(d.get(k * 4..k * 4 + 4)?.try_into().ok()?);
            (v, v == INT32_MISSING, v == INT32_VECTOR_END)
        }
        _ => return None,
    })
}

/// Copies comma-separated field `isrc` of `src` into field `idst` of `dst`, unless it is `.` or
/// `dst` already has a value there (bcftools' `copy_string_field`). Returns -1 when `src` lacks
/// the field, -2 when `dst` does.
pub(super) fn copy_string_field(src: &[u8], isrc: usize, dst: &mut Vec<u8>, idst: usize) -> i32 {
    let src_len = src.len();
    let (mut ith_src, mut start_src) = (0, 0);
    while ith_src < isrc && start_src < src_len {
        if src[start_src] == b',' {
            ith_src += 1;
        }
        start_src += 1;
    }
    if ith_src != isrc {
        return -1;
    }
    let mut end_src = start_src;
    while end_src < src_len && src[end_src] != 0 && src[end_src] != b',' {
        end_src += 1;
    }
    let nsrc_cpy = end_src - start_src;
    if nsrc_cpy == 1 && src[start_src] == b'.' {
        return 0;
    }
    let (mut ith_dst, mut start_dst) = (0, 0);
    while ith_dst < idst && start_dst < dst.len() {
        if dst[start_dst] == b',' {
            ith_dst += 1;
        }
        start_dst += 1;
    }
    if ith_dst != idst {
        return -2;
    }
    let mut end_dst = start_dst;
    while end_dst < dst.len() && dst[end_dst] != b',' {
        end_dst += 1;
    }
    if end_dst - start_dst > 1 || dst.get(start_dst) != Some(&b'.') {
        return 0;
    }
    dst.splice(start_dst..end_dst, src[start_src..end_src].iter().copied());
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copies_string_fields() {
        let mut dst = b".,.,.".to_vec();
        assert_eq!(copy_string_field(b"a,bb", 1, &mut dst, 2), 0);
        assert_eq!(dst, b".,.,bb");
        assert_eq!(copy_string_field(b"a,bb", 0, &mut dst, 0), 0);
        assert_eq!(dst, b"a,.,bb");
        // Values are not overwritten, and missing values are not copied.
        assert_eq!(copy_string_field(b"x", 0, &mut dst, 0), 0);
        assert_eq!(copy_string_field(b".", 0, &mut dst, 1), 0);
        assert_eq!(dst, b"a,.,bb");
        assert_eq!(copy_string_field(b"x", 3, &mut dst, 0), -1);
        assert_eq!(copy_string_field(b"x", 0, &mut dst, 5), -2);
    }
}
