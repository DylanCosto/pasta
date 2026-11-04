//! `-i/--info-rules` and `-M/--missing-rules` (bcftools merge's `info_rules_*` and
//! `missing_rules_*`).

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Number, Record, ValueType,
    encoding::{self, INT32_MISSING, INT32_VECTOR_END},
    header::{HeaderRecord, Quoting},
    record::{InfoField, TypedValues},
};

/// How an INFO rule combines the values of the merged records.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Method {
    Sum,
    Avg,
    Min,
    Max,
    Join,
}

/// The values an INFO rule collects (`rule->vals`), by the tag's type.
#[derive(Debug)]
pub(super) enum RuleValues {
    Int(Vec<i32>),
    Float(Vec<f32>),
    Str(Vec<u8>),
}

/// One `TAG:METHOD` of `-i` (`info_rule_t`).
#[derive(Debug)]
pub(super) struct InfoRule {
    pub(super) tag: Vec<u8>,
    /// The tag's dictionary ID in the output header.
    pub(super) out_id: u32,
    pub(super) method: Method,
    /// The number of values in a block (one block per merged record).
    pub(super) block_size: usize,
    pub(super) nblocks: usize,
    /// The number of values (or string bytes) collected.
    pub(super) nvals: usize,
    pub(super) vals: RuleValues,
}

/// What `-M` fills in for alleles a file lacks when it has no `<*>` (`MERGE_MISSING_*`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum MissingRule {
    /// A missing value (`.`).
    Dot,
    /// A constant.
    Const(f32),
    /// The largest of the sample's values.
    Max,
}

/// Splits `TAG:METHOD,TAG:METHOD` as bcftools does, checking the separators; `None` when they
/// do not alternate.
fn split_rules(s: &str) -> Option<Vec<(&str, &str)>> {
    let mut n = 0;
    for c in s.bytes() {
        match c {
            b':' => {
                n += 1;
                if n % 2 == 0 {
                    return None;
                }
            }
            b',' => {
                n += 1;
                if n % 2 == 1 {
                    return None;
                }
            }
            _ => {}
        }
    }
    if n % 2 == 0 {
        return None;
    }
    Some(
        s.split(',')
            .map(|r| r.split_once(':').unwrap_or((r, "")))
            .collect(),
    )
}

/// Whether an INFO tag is defined in the header.
fn info_exists(hdr: &Header, tag: &[u8]) -> bool {
    hdr.id(tag)
        .is_some_and(|id| hdr.id_exists(HeaderRecordType::Info, id))
}

/// Which default rules apply without `-i` and `-M`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Defaults {
    /// Those of plain VCF merging.
    Vcf,
    /// Those of gVCF merging (`-g`), which also merge the block statistics.
    Gvcf,
}

/// The rules used without `-i` (bcftools' defaults).
fn default_info_rules(hdr: &Header, defaults: Defaults) -> String {
    let gvcf = defaults == Defaults::Gvcf;
    let mut rules: Vec<&str> = Vec::new();
    let mut add = |tag: &str, rule: &'static str, when: bool| {
        if when && info_exists(hdr, tag.as_bytes()) {
            rules.push(rule);
        }
    };
    add("DP", "DP:sum", true);
    add("DP4", "DP4:sum", true);
    add("QS", "QS:sum", gvcf);
    add("MIN_DP", "MIN_DP:min", gvcf);
    add("MinDP", "MinDP:min", gvcf);
    add("I16", "I16:sum", gvcf);
    add("IDV", "IDV:max", gvcf);
    add("IMF", "IMF:max", gvcf);
    let no_samples = hdr.sample_count() == 0;
    add("AN", "AN:sum", no_samples);
    add("AC", "AC:sum", no_samples);
    rules.join(",")
}

/// Parses `-i` against the output header (bcftools' `info_rules_init`), redefining joined
/// Number=A/G/R tags as Number=. in the header. Returns the rules sorted by tag and whether AC
/// or AN has a rule (which keeps them from being recomputed).
pub(super) fn init_info_rules(
    spec: Option<&str>,
    hdr: &mut Header,
    defaults: Defaults,
) -> anyhow::Result<(Vec<InfoRule>, bool)> {
    if spec == Some("-") {
        return Ok((Vec::new(), false));
    }
    let owned;
    let spec = match spec {
        Some(s) => s,
        None => {
            owned = default_info_rules(hdr, defaults);
            if owned.is_empty() {
                return Ok((Vec::new(), false));
            }
            &owned
        }
    };
    let Some(pairs) = split_rules(spec) else {
        bail!("Could not parse INFO rules: \"{spec}\"");
    };
    let mut rules = Vec::with_capacity(pairs.len());
    let mut keep_ac_an = false;
    for (tag, method) in pairs {
        let id = hdr
            .id(tag.as_bytes())
            .filter(|&id| hdr.id_exists(HeaderRecordType::Info, id));
        let Some(id) = id else {
            bail!("The INFO tag is not defined in the header: \"{tag}\"");
        };
        let vals = match hdr.id_type(HeaderRecordType::Info, id) {
            Some(ValueType::Integer) => RuleValues::Int(Vec::new()),
            Some(ValueType::Float) => RuleValues::Float(Vec::new()),
            Some(ValueType::String) => RuleValues::Str(Vec::new()),
            other => bail!(
                "The INFO rule \"{tag}\" is not supported; the tag \"{tag}\" type is {}",
                other.map_or(15, |t| t as u8)
            ),
        };
        if tag == "AC" || tag == "AN" {
            keep_ac_an = true;
        }
        if method.is_empty() {
            bail!("Could not parse INFO rules, missing logic of \"{tag}\"");
        }
        let method = match method.to_ascii_lowercase().as_str() {
            "sum" => Method::Sum,
            "avg" => Method::Avg,
            "min" => Method::Min,
            "max" => Method::Max,
            "join" => Method::Join,
            _ => bail!("The rule logic \"{method}\" not recognised"),
        };
        if method != Method::Join && matches!(vals, RuleValues::Str(_)) {
            bail!(
                "Numeric operation \"{}\" requested on non-numeric field: {tag}",
                spec_method(spec, tag)
            );
        }
        let number = hdr.id_number(HeaderRecordType::Info, id);
        if !matches!(number, Some(Number::Fixed(_))) {
            let is_agr = matches!(number, Some(Number::A | Number::G | Number::R));
            if method == Method::Join && number != Some(Number::Variable) {
                set_number_variable(hdr, tag.as_bytes())?;
            }
            if method != Method::Join && !is_agr {
                bail!(
                    "Only fixed-length vectors are supported with -i {}:{tag}",
                    spec_method(spec, tag)
                );
            }
        }
        rules.push(InfoRule {
            tag: tag.as_bytes().to_vec(),
            out_id: id,
            method,
            block_size: 0,
            nblocks: 0,
            nvals: 0,
            vals,
        });
    }
    // A redefined (joined) tag keeps its ID; look it up again all the same.
    for r in &mut rules {
        r.out_id = hdr.id(&r.tag).unwrap_or(r.out_id);
    }
    rules.sort_by(|a, b| a.tag.cmp(&b.tag));
    Ok((rules, keep_ac_an))
}

/// The method text of a tag's rule, as the user wrote it.
fn spec_method<'a>(spec: &'a str, tag: &str) -> &'a str {
    spec.split(',')
        .filter_map(|r| r.split_once(':'))
        .find(|(t, _)| *t == tag)
        .map_or("", |(_, m)| m)
}

/// Redefines an INFO tag as Number=. (moving its line to the end, as bcftools does).
fn set_number_variable(hdr: &mut Header, tag: &[u8]) -> anyhow::Result<()> {
    let Some(rec) = hdr.definition(HeaderRecordType::Info, tag) else {
        bail!(
            "Uh, could not find the entry Number in the header record of {}",
            String::from_utf8_lossy(tag)
        );
    };
    let mut rec: HeaderRecord = rec.dup();
    let Some(i) = rec.find_key(b"Number") else {
        bail!(
            "Uh, could not find the entry Number in the header record of {}",
            String::from_utf8_lossy(tag)
        );
    };
    rec.set_field_value(i, b".", Quoting::Bare);
    hdr.remove(HeaderRecordType::Info, Some(tag));
    hdr.add_record(rec)?;
    Ok(())
}

/// Parses `-M` against the output header (bcftools' `missing_rules_init`). Without `-M`, gVCF
/// merging uses `PL:max,AD:0` for the tags the header defines.
pub(super) fn init_missing_rules(
    spec: Option<&str>,
    hdr: &Header,
    defaults: Defaults,
) -> anyhow::Result<Vec<(Vec<u8>, MissingRule)>> {
    let text = match spec {
        Some("-") => "PL:.,AD:.",
        Some(s) => s,
        None if defaults == Defaults::Gvcf => "PL:max,AD:0",
        None => return Ok(Vec::new()),
    };
    let shown = spec.unwrap_or("(null)");
    let Some(pairs) = split_rules(text) else {
        bail!("Could not parse --missing-rules: \"{shown}\"");
    };
    let mut rules = Vec::new();
    for (tag, method) in pairs {
        let defined = hdr
            .id(tag.as_bytes())
            .is_some_and(|id| hdr.id_exists(HeaderRecordType::Format, id));
        if !defined {
            if spec.is_some() {
                bail!("The FORMAT tag is not defined in the header: \"{tag}\"");
            }
            if method.is_empty() {
                bail!("Could not parse --missing-rules, missing logic of \"{tag}\"");
            }
            continue;
        }
        if method.is_empty() {
            bail!("Could not parse --missing-rules, missing logic of \"{tag}\"");
        }
        let rule = if method.eq_ignore_ascii_case(".") {
            MissingRule::Dot
        } else if method.eq_ignore_ascii_case("max") {
            MissingRule::Max
        } else {
            match strtod(method) {
                Some(v) => MissingRule::Const(v as f32),
                None => bail!("Could not parse --missing-rules: \"{shown}\""),
            }
        };
        rules.push((tag.as_bytes().to_vec(), rule));
    }
    rules.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(rules)
}

/// C's `strtod` over the whole string: `None` unless all of it is a number.
fn strtod(s: &str) -> Option<f64> {
    let t = s.trim_start();
    if t.is_empty() {
        return None;
    }
    if let Ok(v) = t.parse::<f64>() {
        return Some(v);
    }
    let lower = t.to_ascii_lowercase();
    let (neg, body) = match lower.as_bytes().first() {
        Some(b'-') => (true, &lower[1..]),
        Some(b'+') => (false, &lower[1..]),
        _ => (false, &lower[..]),
    };
    let v = match body {
        "inf" | "infinity" => f64::INFINITY,
        "nan" => f64::NAN,
        _ => {
            let hex = body.strip_prefix("0x")?;
            i64::from_str_radix(hex, 16).ok()? as f64
        }
    };
    Some(if neg { -v } else { v })
}

impl InfoRule {
    /// Forgets the values of the previous site (`info_rules_reset`).
    pub(super) fn reset(&mut self) {
        self.nblocks = 0;
        self.nvals = 0;
        self.block_size = 0;
    }

    /// Adds a record's values (bcftools' `info_rules_add_values`). `map` maps the record's
    /// alleles to the output alleles for Number=A/G/R tags (`number`, from the record's
    /// header); `nout_als` is the number of output alleles.
    #[expect(clippy::too_many_arguments, reason = "bcftools' info_rules_add_values")]
    pub(super) fn add_values(
        &mut self,
        hdr: &Header,
        line: &Record,
        inf: &InfoField,
        number: Option<Number>,
        map: &[usize],
        nout_als: usize,
        agr_map: &mut Vec<usize>,
        tmp_i: &mut Vec<i32>,
        tmp_f: &mut Vec<f32>,
    ) -> anyhow::Result<()> {
        let chrom = |line: &Record| String::from_utf8_lossy(line.chrom_safe(hdr)).into_owned();
        let tag = TagName(&self.tag);
        let ret = match &self.vals {
            RuleValues::Int(_) => info_i32(line, hdr, inf, tmp_i),
            RuleValues::Float(_) => info_f32(line, hdr, inf, tmp_f),
            RuleValues::Str(_) => info_str(line, hdr, inf).map(<[u8]>::len),
        };
        let ret = match ret {
            Ok(n) if n > 0 => n,
            Ok(n) => bail!(
                "FIXME: error parsing {tag} at {}:{} .. {n}",
                chrom(line),
                line.pos() + 1
            ),
            Err(e) => bail!(
                "FIXME: error parsing {tag} at {}:{} .. {e}",
                chrom(line),
                line.pos() + 1
            ),
        };
        self.nblocks += 1;

        if let RuleValues::Str(vals) = &mut self.vals {
            let s = info_str(line, hdr, inf).unwrap_or_default();
            vals.truncate(self.nvals);
            if self.nvals > 0 {
                vals.push(b',');
            }
            // strncpy: up to a NUL, padded with NULs.
            let end = s.iter().position(|&c| c == 0).unwrap_or(s.len());
            vals.extend_from_slice(&s[..end]);
            vals.resize(vals.len() + (s.len() - end), 0);
            self.nvals += ret + usize::from(self.nblocks > 1);
            return Ok(());
        }

        let n_allele = line.allele_count();
        let mut use_map = true;
        match number {
            Some(Number::A) => {
                if ret != n_allele - 1 {
                    bail!(
                        "Wrong number of {tag} fields at {}:{}",
                        chrom(line),
                        line.pos() + 1
                    );
                }
                agr_map.clear();
                agr_map.extend((0..ret).map(|i| map[i + 1] - 1));
                self.block_size = nout_als - 1;
            }
            Some(Number::R) => {
                if ret != n_allele {
                    bail!(
                        "Wrong number of {tag} fields at {}:{}",
                        chrom(line),
                        line.pos() + 1
                    );
                }
                agr_map.clear();
                agr_map.extend_from_slice(&map[..ret]);
                self.block_size = nout_als;
            }
            Some(Number::G) => {
                let ngt = n_allele * (n_allele + 1) / 2;
                if ret != n_allele && ret != ngt {
                    bail!(
                        "Wrong number of {tag} fields at {}:{}",
                        chrom(line),
                        line.pos() + 1
                    );
                }
                agr_map.clear();
                if ret == n_allele {
                    agr_map.extend_from_slice(&map[..ret]);
                    self.block_size = nout_als;
                } else {
                    for i in 0..n_allele {
                        for j in 0..=i {
                            agr_map.push(super::alleles::alleles2gt(map[i], map[j]));
                        }
                    }
                    self.block_size = nout_als * (nout_als + 1) / 2;
                }
            }
            _ => {
                if self.nblocks > 1 && ret != self.block_size {
                    bail!(
                        "Mismatch in number of values for INFO/{tag} at {}:{}",
                        chrom(line),
                        line.pos() + 1
                    );
                }
                self.block_size = ret;
                use_map = false;
            }
        }
        let start = self.nvals;
        self.nvals += self.block_size;
        match &mut self.vals {
            RuleValues::Int(vals) => {
                fill_block(
                    vals,
                    start,
                    self.block_size,
                    &tmp_i[..ret],
                    use_map.then_some(&agr_map[..]),
                    INT32_MISSING,
                );
            }
            RuleValues::Float(vals) => {
                fill_block(
                    vals,
                    start,
                    self.block_size,
                    &tmp_f[..ret],
                    use_map.then_some(&agr_map[..]),
                    encoding::float_missing(),
                );
            }
            RuleValues::Str(_) => unreachable!("strings are handled above"),
        }
        Ok(())
    }

    /// Sets the merged values on the output record (`rule->merger`).
    pub(super) fn apply(&mut self, hdr: &Header, out: &mut Record) {
        if self.nvals == 0 {
            return;
        }
        let ndim = self.block_size;
        let nblocks = self.nblocks;
        // Blocks of different sizes (a Number=G tag given per allele at a REF-only site and
        // per genotype elsewhere) make bcftools read past the values, into its zeroed buffer
        // (join only reads the values).
        let n = if self.method == Method::Join {
            self.nvals
        } else {
            self.nvals.max(ndim * nblocks)
        };
        match &mut self.vals {
            RuleValues::Int(v) => {
                v.truncate(self.nvals);
                v.resize(n, 0);
                let len = merge_ints(self.method, v, ndim, nblocks);
                let _ = out.update_info_i32_id(hdr, self.out_id, &v[..len]);
            }
            RuleValues::Float(v) => {
                v.truncate(self.nvals);
                v.resize(n, 0.0);
                let len = merge_floats(self.method, v, ndim, nblocks);
                let _ = out.update_info_f32_id(hdr, self.out_id, &v[..len]);
            }
            RuleValues::Str(v) => {
                if self.method == Method::Join {
                    v.truncate(self.nvals);
                    let _ = out.update_info_string(hdr, &self.tag, Some(v));
                }
            }
        }
    }
}

/// A tag name for messages, formatted only when used.
struct TagName<'a>(&'a [u8]);

impl std::fmt::Display for TagName<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(self.0))
    }
}

/// Checks that an INFO field's header type is `ty`, as `bcf_get_info_values` does (-2 when
/// not).
fn info_field<'a>(field: &'a InfoField, hdr: &Header, ty: ValueType) -> Result<&'a InfoField, i32> {
    if hdr.id_type(HeaderRecordType::Info, field.key() as u32) != Some(ty) {
        return Err(-2);
    }
    Ok(field)
}

/// `bcf_get_info_int32` on a record's field: missing values become [`INT32_MISSING`], a
/// vector end stops.
fn info_i32(
    line: &Record,
    hdr: &Header,
    field: &InfoField,
    out: &mut Vec<i32>,
) -> Result<usize, i32> {
    let field = info_field(field, hdr, ValueType::Integer)?;
    let values = line.info_values(field);
    out.clear();
    for j in 0..field.len() {
        match super::site::raw_int(&values, j) {
            Some((_, _, true)) => break,
            Some((_, true, _)) => out.push(INT32_MISSING),
            Some((v, ..)) => out.push(v),
            None => return Err(-2),
        }
    }
    Ok(out.len())
}

/// `bcf_get_info_float` on a record's field; a vector end stops.
fn info_f32(
    line: &Record,
    hdr: &Header,
    field: &InfoField,
    out: &mut Vec<f32>,
) -> Result<usize, i32> {
    let field = info_field(field, hdr, ValueType::Float)?;
    let TypedValues::Float(data) = line.info_values(field) else {
        return Err(-2);
    };
    out.clear();
    for c in data.as_chunks::<4>().0.iter().take(field.len()) {
        let bits = u32::from_le_bytes(*c);
        if bits == encoding::FLOAT_VECTOR_END {
            break;
        }
        out.push(f32::from_bits(bits));
    }
    Ok(out.len())
}

/// `bcf_get_info_string` on a record's field.
fn info_str<'a>(line: &'a Record, hdr: &Header, field: &'a InfoField) -> Result<&'a [u8], i32> {
    let field = info_field(field, hdr, ValueType::String)?;
    let data = line.info_data(field);
    Ok(&data[..field.len().min(data.len())])
}

/// Appends a block of values, mapped to output allele indexes when `map` is given (and the
/// rest missing), else copied.
fn fill_block<T: Copy>(
    vals: &mut Vec<T>,
    start: usize,
    block_size: usize,
    src: &[T],
    map: Option<&[usize]>,
    missing: T,
) {
    vals.truncate(start);
    vals.resize(start + block_size.max(src.len()), missing);
    match map {
        None => vals[start..start + src.len()].copy_from_slice(src),
        Some(map) => {
            for (i, &v) in src.iter().enumerate() {
                vals[start + map[i]] = v;
            }
        }
    }
    vals.truncate(start + block_size);
}

/// Merges integer blocks in place; returns the number of output values.
fn merge_ints(method: Method, v: &mut [i32], ndim: usize, nblocks: usize) -> usize {
    match method {
        Method::Sum => {
            for x in v.iter_mut() {
                if *x == INT32_MISSING {
                    *x = 0;
                }
            }
            for i in 1..nblocks {
                for j in 0..ndim {
                    v[j] = v[j].wrapping_add(v[j + i * ndim]);
                }
            }
            ndim
        }
        Method::Avg => {
            for x in v.iter_mut() {
                if *x == INT32_MISSING {
                    *x = 0;
                }
            }
            for j in 0..ndim {
                let sum: f64 = (0..nblocks).map(|i| f64::from(v[j + i * ndim])).sum();
                v[j] = (sum / nblocks as f64) as i32;
            }
            ndim
        }
        Method::Min | Method::Max => {
            let huge = if method == Method::Min {
                i32::MAX
            } else {
                i32::MIN
            };
            for x in v.iter_mut() {
                if *x == INT32_MISSING {
                    *x = huge;
                }
            }
            for i in 1..nblocks {
                for j in 0..ndim {
                    let o = v[j + i * ndim];
                    if (method == Method::Min && v[j] > o) || (method == Method::Max && v[j] < o) {
                        v[j] = o;
                    }
                }
            }
            for x in v.iter_mut() {
                if *x == huge {
                    *x = INT32_MISSING;
                }
            }
            ndim
        }
        Method::Join => {
            let mut idst = 0;
            for isrc in 0..v.len() {
                if v[isrc] == INT32_VECTOR_END {
                    break;
                }
                if v[isrc] == INT32_MISSING {
                    continue;
                }
                v[idst] = v[isrc];
                idst += 1;
            }
            idst
        }
    }
}

/// Merges float blocks in place; returns the number of output values.
fn merge_floats(method: Method, v: &mut [f32], ndim: usize, nblocks: usize) -> usize {
    let missing = |x: f32| encoding::float_is_missing(x);
    match method {
        Method::Sum => {
            for x in v.iter_mut() {
                if missing(*x) {
                    *x = 0.0;
                }
            }
            for i in 1..nblocks {
                for j in 0..ndim {
                    v[j] += v[j + i * ndim];
                }
            }
            ndim
        }
        Method::Avg => {
            for x in v.iter_mut() {
                if missing(*x) {
                    *x = 0.0;
                }
            }
            for j in 0..ndim {
                let sum: f64 = (0..nblocks).map(|i| f64::from(v[j + i * ndim])).sum();
                v[j] = (sum / nblocks as f64) as f32;
            }
            ndim
        }
        Method::Min | Method::Max => {
            let huge = if method == Method::Min {
                f32::INFINITY
            } else {
                f32::NEG_INFINITY
            };
            for x in v.iter_mut() {
                if missing(*x) {
                    *x = huge;
                }
            }
            for i in 1..nblocks {
                for j in 0..ndim {
                    let o = v[j + i * ndim];
                    if (method == Method::Min && v[j] > o) || (method == Method::Max && v[j] < o) {
                        v[j] = o;
                    }
                }
            }
            for x in v.iter_mut() {
                if *x == huge {
                    *x = encoding::float_missing();
                }
            }
            ndim
        }
        Method::Join => {
            let mut idst = 0;
            for isrc in 0..v.len() {
                if encoding::float_is_vector_end(v[isrc]) {
                    break;
                }
                if missing(v[isrc]) {
                    continue;
                }
                v[idst] = v[isrc];
                idst += 1;
            }
            idst
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_rules_like_bcftools() {
        assert_eq!(split_rules("DP:sum"), Some(vec![("DP", "sum")]));
        assert_eq!(
            split_rules("DP:sum,AC:max"),
            Some(vec![("DP", "sum"), ("AC", "max")])
        );
        assert_eq!(split_rules("DP"), None);
        assert_eq!(split_rules("DP:sum,"), None);
        assert_eq!(split_rules("DP:sum:max"), None);
        assert_eq!(split_rules(""), None);
        assert_eq!(split_rules("DP:"), Some(vec![("DP", "")]));
    }

    #[test]
    fn merges_blocks() {
        let mut v = vec![1, INT32_MISSING, 3, 4];
        assert_eq!(merge_ints(Method::Sum, &mut v, 2, 2), 2);
        assert_eq!(&v[..2], [4, 4]);
        let mut v = vec![1, INT32_MISSING, 3, INT32_MISSING];
        assert_eq!(merge_ints(Method::Max, &mut v, 2, 2), 2);
        assert_eq!(&v[..2], [3, INT32_MISSING]);
        let mut v = vec![1, 2, 4, 8];
        assert_eq!(merge_ints(Method::Avg, &mut v, 2, 2), 2);
        assert_eq!(&v[..2], [2, 5]);
        let mut v = vec![1, INT32_MISSING, 3];
        assert_eq!(merge_ints(Method::Join, &mut v, 1, 3), 2);
        assert_eq!(&v[..2], [1, 3]);
    }

    #[test]
    fn parses_strtod_values() {
        assert_eq!(strtod("0"), Some(0.0));
        assert_eq!(strtod(" 1.5"), Some(1.5));
        assert_eq!(strtod("1x"), None);
        assert_eq!(strtod("-inf"), Some(f64::NEG_INFINITY));
    }
}
