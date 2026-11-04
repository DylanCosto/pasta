//! Arithmetic, comparisons, and logic on tokens (bcftools' `VECTOR_ARITHMETICS`,
//! `CMP_VECTORS`, `cmp_vector_strings`, the `filters_cmp_*` comparators, and
//! `vector_logic_or/and`).

use anyhow::bail;
use pasta::variant::Record;

use super::{
    bre::Bre,
    token::{Op, Token, expand, is_absent, is_vector_end, missing},
};

/// Marks the samples a binary operation uses (the union of its operands') once, and clears the
/// per-sample results (`tok_init_samples`).
pub(super) fn init_samples(a: &Token, b: &Token, r: &mut Token) {
    if (a.nsamples > 0 || b.nsamples > 0) && r.nsamples == 0 {
        r.nsamples = if a.nsamples > 0 {
            a.nsamples
        } else {
            b.nsamples
        };
        r.usmpl = vec![false; r.nsamples];
        for i in 0..a.nsamples.min(r.nsamples) {
            r.usmpl[i] |= a.uses(i);
        }
        for i in 0..b.nsamples.min(r.nsamples) {
            r.usmpl[i] |= b.uses(i);
        }
    }
    if r.nsamples > 0 {
        if r.pass_samples.len() < r.nsamples {
            r.pass_samples.resize(r.nsamples, false);
        }
        r.pass_samples[..r.nsamples].fill(false);
    }
}

/// Sizes the result of an arithmetic operation like its larger operand (`tok_init_values`).
fn init_values(a: &Token, b: &Token, r: &mut Token) {
    let tok = if (a.nsamples > 0 || b.nsamples > 0) && (a.nsamples == 0 || b.nsamples == 0) {
        if a.nsamples > 0 { a } else { b }
    } else if a.nvalues > b.nvalues {
        a
    } else {
        b
    };
    r.nvalues = tok.nvalues;
    r.nval1 = tok.nval1;
    expand(&mut r.values, r.nvalues);
}

/// An arithmetic operator.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Arith {
    Add,
    Sub,
    Mult,
    Div,
    Modulo,
}

impl Arith {
    #[inline]
    fn apply(self, x: f64, y: f64) -> f64 {
        match self {
            Self::Add => x + y,
            Self::Sub => x - y,
            Self::Mult => x * y,
            Self::Div => x / y,
            // Integer modulo, as bcftools computes it; bcftools crashes on a zero divisor.
            Self::Modulo => (x as i32)
                .checked_rem(y as i32)
                .map_or_else(missing, f64::from),
        }
    }
}

/// Applies an arithmetic operator element-wise, broadcasting a single value
/// (`VECTOR_ARITHMETICS`).
pub(super) fn arithmetic(a: &Token, b: &Token, r: &mut Token, op: Arith) -> anyhow::Result<()> {
    let mut has_values = false;
    if a.nvalues > 0 && b.nvalues > 0 {
        init_values(a, b, r);
        init_samples(a, b, r);
        if a.nsamples == 0 && b.nsamples == 0 {
            if a.nvalues != b.nvalues && a.nvalues != 1 && b.nvalues != 1 {
                bail!(
                    "Cannot run numeric operator in -i/-e filtering on vectors of different lengths: {} vs {}",
                    a.nvalues,
                    b.nvalues
                );
            }
            for ir in 0..r.nvalues {
                let ia = if a.nvalues > 1 { ir } else { 0 };
                let ib = if b.nvalues > 1 { ir } else { 0 };
                let (x, y) = (a.value(ia), b.value(ib));
                if is_absent(x) || is_absent(y) {
                    r.values[ir] = missing();
                    continue;
                }
                has_values = true;
                r.values[ir] = op.apply(x, y);
            }
        } else if a.nsamples > 0 && b.nsamples > 0 {
            if a.nsamples != b.nsamples {
                bail!(
                    "Incompatible number samples in comparison: {} vs {}",
                    a.nsamples,
                    b.nsamples
                );
            }
            if a.nval1 != b.nval1 && a.nval1 != 1 && b.nval1 != 1 {
                bail!(
                    "Cannot run numeric operator in -i/-e filtering on vectors of different lengths: {} vs {}",
                    a.nval1,
                    b.nval1
                );
            }
            let n1 = r.nval1;
            expand(&mut r.values, r.nsamples * n1);
            for i in 0..r.nsamples {
                for ir in 0..n1 {
                    let ia = if a.nval1 > 1 { ir } else { 0 };
                    let ib = if b.nval1 > 1 { ir } else { 0 };
                    let (x, y) = (a.sample_value(i, ia), b.sample_value(i, ib));
                    let out = &mut r.values[i * n1 + ir];
                    if is_absent(x) || is_absent(y) {
                        *out = missing();
                        continue;
                    }
                    has_values = true;
                    *out = op.apply(x, y);
                }
            }
        } else {
            let (s, v, sample_first) = if a.nsamples > 0 {
                (a, b, true)
            } else {
                (b, a, false)
            };
            // bcftools asserts that the other side has one value, but is built without
            // assertions: its first value is used.
            let y = v.value(0);
            if !is_absent(y) {
                for i in 0..s.nvalues {
                    let x = s.value(i);
                    if is_absent(x) {
                        r.values[i] = missing();
                        continue;
                    }
                    has_values = true;
                    r.values[i] = if sample_first {
                        op.apply(x, y)
                    } else {
                        op.apply(y, x)
                    };
                }
            }
        }
    }
    if !has_values {
        r.nvalues = 0;
    }
    Ok(())
}

/// 2^24, above which not every integer is a float: larger values are compared as doubles.
const FLOAT_EXACT_MAX: f64 = 16_777_216.0;

/// A numeric comparison operator.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Cmp {
    Eq,
    Ne,
    Le,
    Lt,
    Bt,
    Be,
}

impl Cmp {
    #[inline]
    fn test<T: PartialOrd>(self, a: T, b: T) -> bool {
        match self {
            Self::Eq => a == b,
            Self::Ne => a != b,
            Self::Le => a <= b,
            Self::Lt => a < b,
            Self::Bt => a > b,
            Self::Be => a >= b,
        }
    }

    /// Compares as floats, as values read from VCF are, unless a value is too large for a
    /// float to hold exactly (bcftools issue #871).
    #[inline]
    fn numbers(self, a: f64, b: f64) -> bool {
        if a > FLOAT_EXACT_MAX || b > FLOAT_EXACT_MAX {
            self.test(a, b)
        } else {
            self.test(a as f32, b as f32)
        }
    }

    /// What to do with one or both sides missing: `[either, one, both]`.
    fn missing_logic(self) -> [bool; 3] {
        match self {
            Self::Eq => [true, false, true],
            Self::Ne => [true, true, false],
            _ => [false; 3],
        }
    }
}

#[inline]
fn miss_count(x: f64) -> usize {
    usize::from(is_absent(x))
}

/// Compares two numeric tokens, site-wise or per sample (`CMP_VECTORS`). A value matching
/// any value of the other side passes.
pub(super) fn compare_numbers(a: &Token, b: &Token, r: &mut Token, op: Cmp) -> anyhow::Result<()> {
    let ml = op.missing_logic();
    init_samples(a, b, r);
    if a.nsamples == 0 && b.nsamples == 0 {
        compare_site_values(a, b, r, op, ml);
        return Ok(());
    }

    if a.nvalues == 0 && b.nvalues == 0 {
        if ml[2] {
            pass_all_used(r);
        }
        return Ok(());
    }
    if a.nvalues == 0 || b.nvalues == 0 {
        compare_with_empty_side(a, b, r, ml);
        return Ok(());
    }
    if a.nsamples > 0 && b.nsamples > 0 {
        compare_sample_values(a, b, r, op, ml)?;
    } else {
        compare_samples_to_site(a, b, r, op, ml);
    }
    Ok(())
}

/// [`compare_numbers`] of two site-wise tokens.
#[inline]
fn compare_site_values(a: &Token, b: &Token, r: &mut Token, op: Cmp, ml: [bool; 3]) {
    if a.nvalues == 0 && b.nvalues == 0 {
        r.pass_site = ml[2];
    } else if a.nvalues == 0 || b.nvalues == 0 {
        let tok = if a.nvalues > 0 { a } else { b };
        for j in 0..tok.nvalues {
            if if is_absent(tok.value(j)) {
                ml[2]
            } else {
                ml[1]
            } {
                r.pass_site = true;
                break;
            }
        }
    } else {
        'outer: for i in 0..a.nvalues {
            let x = a.value(i);
            let amiss = miss_count(x);
            for j in 0..b.nvalues {
                let y = b.value(j);
                let nmiss = amiss + miss_count(y);
                let pass = if nmiss > 0 {
                    ml[nmiss]
                } else {
                    op.numbers(x, y)
                };
                if pass {
                    r.pass_site = true;
                    break 'outer;
                }
            }
        }
    }
}

/// [`compare_numbers`] with a per-sample operand when one side has no values.
#[inline]
fn compare_with_empty_side(a: &Token, b: &Token, r: &mut Token, ml: [bool; 3]) {
    let tok = if a.nvalues > 0 { a } else { b };
    if tok.nsamples == 0 {
        let miss = (0..tok.nvalues).any(|j| is_absent(tok.value(j)));
        if ml[usize::from(miss) + 1] {
            pass_all_used(r);
        }
    } else {
        for i in 0..tok.nsamples {
            if !r.uses(i) {
                continue;
            }
            let miss = (0..tok.nval1).any(|j| is_absent(tok.sample_value(i, j)));
            if ml[usize::from(miss) + 1] {
                r.pass_samples[i] = true;
                r.pass_site = true;
            }
        }
    }
}

/// [`compare_numbers`] of two per-sample tokens, sample by sample.
#[inline]
fn compare_sample_values(
    a: &Token,
    b: &Token,
    r: &mut Token,
    op: Cmp,
    ml: [bool; 3],
) -> anyhow::Result<()> {
    if a.nval1 != b.nval1 {
        bail!(
            "Incompatible number of per-sample values in comparison: {} vs {}",
            a.nval1,
            b.nval1
        );
    }
    if a.nsamples != b.nsamples {
        bail!(
            "Incompatible number samples in comparison: {} vs {}",
            a.nsamples,
            b.nsamples
        );
    }
    for i in 0..a.nsamples {
        if !a.uses(i) || !b.uses(i) {
            if let Some(u) = r.usmpl.get_mut(i) {
                *u = false;
            }
            continue;
        }
        for j in 0..a.nval1 {
            let (x, y) = (a.sample_value(i, j), b.sample_value(i, j));
            // With [GT] subscripts, the vectors are truncated (see bcftools' note).
            if (a.idx == -3 && is_vector_end(x)) || (b.idx == -3 && is_vector_end(y)) {
                break;
            }
            let mut nmiss = miss_count(x);
            if nmiss > 0 && !ml[0] {
                continue;
            }
            nmiss += miss_count(y);
            let pass = if nmiss > 0 {
                ml[nmiss]
            } else {
                op.numbers(x, y)
            };
            if pass {
                r.pass_samples[i] = true;
                r.pass_site = true;
                break;
            }
        }
    }
    Ok(())
}

/// [`compare_numbers`] of a per-sample token with a site-wise one, on either side.
#[inline]
fn compare_samples_to_site(a: &Token, b: &Token, r: &mut Token, op: Cmp, ml: [bool; 3]) {
    let (s, v) = if a.nsamples > 0 { (a, b) } else { (b, a) };
    let sample_first = a.nsamples > 0;
    // [GT]-subscripted vectors end early. bcftools checks the left side's subscript even
    // when the per-sample values are on the right (5!=AD[*:GT]), and so compares the
    // padding as missing; vtools checks the per-sample side.
    for i in 0..s.nsamples {
        if !r.uses(i) {
            continue;
        }
        'sample: for j in 0..s.nval1 {
            let x = s.sample_value(i, j);
            if s.idx == -3 && is_vector_end(x) {
                break;
            }
            let miss = miss_count(x);
            if miss > 0 && !ml[0] {
                continue;
            }
            for k in 0..v.nvalues {
                let y = v.value(k);
                let nmiss = miss + miss_count(y);
                let pass = if nmiss > 0 {
                    ml[nmiss]
                } else if sample_first {
                    op.numbers(x, y)
                } else {
                    op.numbers(y, x)
                };
                if pass {
                    r.pass_samples[i] = true;
                    r.pass_site = true;
                    break 'sample;
                }
            }
        }
    }
}

/// Passes every used sample, and the site.
fn pass_all_used(r: &mut Token) {
    for i in 0..r.nsamples {
        if r.uses(i) {
            r.pass_samples[i] = true;
            r.pass_site = true;
        }
    }
}

/// The comma-separated fields of `s` up to its first NUL, as bcftools' loops visit them: a
/// field starts wherever the previous one ended, so a trailing comma adds no empty field.
#[inline]
fn fields(s: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut pos = 0;
    std::iter::from_fn(move || {
        let rest = s.get(pos..)?;
        if rest.first().is_none_or(|&c| c == 0) {
            return None;
        }
        let len = rest
            .iter()
            .position(|&c| c == b',' || c == 0)
            .unwrap_or(rest.len());
        let field = &rest[..len];
        pos += if rest.get(len) == Some(&b',') {
            len + 1
        } else {
            // The end of the string: no more fields.
            s.len() + 1
        };
        Some(field)
    })
}

/// Whether any field is `.` (`_has_missing_string`).
fn has_missing_string(s: &[u8]) -> bool {
    fields(s).any(|f| f == b".")
}

/// Whether any field of `a` matches any field of `b`, or a missing-value rule applies
/// (`_match_vector_strings`).
fn match_vector_strings(a: &[u8], b: &[u8], op: Op, ml: [bool; 3]) -> bool {
    for x in fields(a) {
        let miss = usize::from(x == b".");
        for y in fields(b) {
            let nmiss = miss + usize::from(y == b".");
            if nmiss > 0 {
                if ml[nmiss] {
                    return true;
                }
            } else if (x == y) != (op == Op::Ne) {
                return true;
            }
        }
    }
    false
}

/// Whether any field of `s` matches the regex (or, with `!~`, does not)
/// (`_regex_vector_strings`).
fn regex_vector_strings(re: &Bre, s: &[u8], op: Op, ml: [bool; 3]) -> bool {
    let negate = op == Op::NLike;
    for f in fields(s) {
        let miss = f == b".";
        if (miss && ml[1]) != negate {
            return true;
        }
        if re.is_match(f) != negate {
            return true;
        }
    }
    false
}

/// Compares two string tokens with `==`, `!=`, `~`, or `!~` (`cmp_vector_strings`).
pub(super) fn compare_strings(a: &Token, b: &Token, r: &mut Token) -> anyhow::Result<()> {
    init_samples(a, b, r);
    let op = r.op;
    let regex = a.regex.as_ref().or(b.regex.as_ref());
    // bcftools asserts that the operator suits the strings, but is built without assertions:
    // other operators compare for equality, without the missing-value rules.
    let ml = match op {
        Op::Eq | Op::Like => [true, false, true],
        Op::Ne | Op::NLike => [true, true, false],
        _ => [false; 3],
    };

    if a.nsamples == 0 && b.nsamples == 0 {
        if a.nvalues == 0 && b.nvalues == 0 {
            r.pass_site = ml[2];
            return Ok(());
        }
        if a.nvalues == 0 || b.nvalues == 0 {
            let tok = if a.nvalues > 0 { a } else { b };
            if ml[usize::from(has_missing_string(&tok.str_value)) + 1] {
                r.pass_site = true;
            }
            return Ok(());
        }
        r.pass_site = match regex {
            None => match_vector_strings(&a.str_value, &b.str_value, op, ml),
            Some(re) => {
                let tok = if a.regex.is_some() { b } else { a };
                regex_vector_strings(re, &tok.str_value, op, ml)
            }
        };
        return Ok(());
    }

    if a.nvalues == 0 && b.nvalues == 0 {
        if ml[2] {
            pass_all_used(r);
        }
        return Ok(());
    }
    if a.nvalues == 0 || b.nvalues == 0 {
        let tok = if a.nvalues > 0 { a } else { b };
        if tok.nsamples == 0 {
            if ml[usize::from(has_missing_string(&tok.str_value)) + 1] {
                pass_all_used(r);
            }
        } else {
            for i in 0..tok.nsamples {
                if r.uses(i) && ml[usize::from(has_missing_string(tok.sample_str(i))) + 1] {
                    r.pass_samples[i] = true;
                    r.pass_site = true;
                }
            }
        }
        return Ok(());
    }

    // With both per sample (another unchecked assertion), y's first sample is compared.
    let (x, y) = if a.nsamples > 0 { (a, b) } else { (b, a) };
    // The usual case, one string such as "het": compare each sample's fields with it.
    let ys = &y.str_value[..memchr::memchr(0, &y.str_value).unwrap_or(y.str_value.len())];
    let single = regex.is_none() && !ys.is_empty() && !ys.contains(&b',');
    let y_miss = usize::from(ys == b".");
    for i in 0..x.nsamples {
        if !r.uses(i) {
            continue;
        }
        let block = x.sample_str(i);
        let matched = if single {
            fields(block).any(|f| {
                let nmiss = usize::from(f == b".") + y_miss;
                if nmiss > 0 {
                    ml[nmiss]
                } else {
                    (f == ys) != (op == Op::Ne)
                }
            })
        } else {
            match regex {
                Some(re) => regex_vector_strings(re, block, op, ml),
                None => match_vector_strings(block, &y.str_value, op, ml),
            }
        };
        if matched {
            r.pass_samples[i] = true;
            r.pass_site = true;
        }
    }
    Ok(())
}

/// `TYPE~"snp"`: a bitwise test of the variant types (`filters_cmp_bit_and`).
pub(super) fn bit_and(a: &Token, b: &Token, r: &mut Token) {
    let x = if a.nvalues > 0 {
        a.value(0)
    } else {
        a.threshold
    } as i32;
    let y = if b.nvalues > 0 {
        b.value(0)
    } else {
        b.threshold
    } as i32;
    r.pass_site = if r.op == Op::Like {
        x & y != 0
    } else {
        x & y == 0
    };
}

/// FILTER `==` (exact), `!=`, `~` (all listed present), and `!~` (`filters_cmp_filter`). The
/// values of `b` are the filter IDs of the expression; no values means `.`.
pub(super) fn compare_filter(b: &Token, r: &mut Token, rec: &Record) -> anyhow::Result<()> {
    let flt = rec.filters();
    let query = &b.values[..b.nvalues.min(b.values.len())];
    let all_present = || {
        query
            .iter()
            .all(|&q| flt.iter().any(|&f| f64::from(f) == q))
    };
    r.pass_site = match r.op {
        Op::NotIn => {
            if query.is_empty() {
                !flt.is_empty()
            } else {
                flt.is_empty() || !all_present()
            }
        }
        Op::In => {
            if query.is_empty() {
                flt.is_empty()
            } else {
                !flt.is_empty() && all_present()
            }
        }
        Op::Ne => query.len() != flt.len() || (!query.is_empty() && !all_present()),
        Op::Eq => query.len() == flt.len() && (query.is_empty() || all_present()),
        _ => bail!("Only ==, !=, ~, and !~ operators are supported for FILTER"),
    };
    Ok(())
}

/// The comma-separated items of a string as bcftools splits them for `@file` lookups: each
/// item runs to the next comma after its first character.
fn hash_items(s: &[u8]) -> impl Iterator<Item = &[u8]> {
    let s = &s[..memchr::memchr(0, s).unwrap_or(s.len())];
    let mut pos = 0;
    std::iter::from_fn(move || {
        if pos >= s.len() {
            return None;
        }
        let start = pos;
        let end = s[start + 1..]
            .iter()
            .position(|&c| c == b',')
            .map_or(s.len(), |k| start + 1 + k);
        pos = end + 1;
        Some(&s[start..end])
    })
}

/// A tag compared with the strings of a file (`filters_cmp_string_hash`).
pub(super) fn compare_hash(a: &Token, b: &Token, r: &mut Token) -> anyhow::Result<()> {
    let (a, b) = if b.hash.is_some() { (b, a) } else { (a, b) };
    if !matches!(r.op, Op::Eq | Op::Ne) {
        bail!("Only == and != operators are supported for strings read from a file");
    }
    if !b.is_str {
        bail!("Error: Lookup of values from a file is currenlty supported only for strings");
    }
    let Some(hash) = &a.hash else {
        bail!("Error: could not parse the expression with \"@file_name\" syntax (possible todo)");
    };
    let negate = r.op == Op::Ne;
    let value = &b.str_value[..b.nvalues.min(b.str_value.len())];
    if b.nsamples == 0 {
        let found = if b.idx >= 0 {
            hash.contains(&value[..memchr::memchr(0, value).unwrap_or(value.len())])
        } else {
            hash_items(value).any(|item| hash.contains(item))
        };
        r.pass_site = found != negate;
        return Ok(());
    }

    init_samples(a, b, r);
    r.pass_site = false;
    for i in 0..b.nsamples {
        if !r.uses(i) {
            continue;
        }
        let block = b.sample_str(i);
        let found = if b.idx >= 0 {
            hash.contains(&block[..memchr::memchr(0, block).unwrap_or(block.len())])
        } else {
            hash_items(block).any(|item| hash.contains(item))
        };
        let pass = found != negate;
        r.pass_samples[i] = pass;
        r.pass_site |= pass;
    }
    Ok(())
}

/// ID `==`, `!=`, `~`, and `!~`, or a lookup in `@file`, true when any of the record's
/// semicolon-separated IDs matches (`filters_cmp_id`).
pub(super) fn compare_id(a: &Token, b: &Token, r: &mut Token, rec: &Record) -> anyhow::Result<()> {
    let (a, b) = if b.hash.is_some() { (b, a) } else { (a, b) };
    let op = r.op;
    let mut pass = false;
    for id in rec.id().split(|&c| c == b';') {
        if let Some(hash) = &a.hash {
            if !matches!(op, Op::Eq | Op::Ne) {
                bail!("Only == and != operators are supported for strings read from a file");
            }
            pass = hash.contains(id);
        } else {
            if b.str_value.is_empty() {
                bail!("Error occurred while evaluating the expression");
            }
            if matches!(op, Op::Eq | Op::Ne) {
                pass = b.c_str() == id;
            } else {
                if !matches!(op, Op::Like | Op::NLike) {
                    bail!(
                        "Only the following operators are supported for querying ID: ==, !=, ~, !~; the operator type {} is not supported ((nil) (nil))",
                        op.code()
                    );
                }
                let Some(re) = a.regex.as_ref().or(b.regex.as_ref()) else {
                    bail!("fixme: regex initialization failed");
                };
                pass = re.is_match(id);
            }
        }
        if pass {
            break;
        }
    }
    if matches!(op, Op::Ne | Op::NLike) {
        pass = !pass;
    }
    r.pass_site = pass;
    Ok(())
}

/// `||` and `|` (`vector_logic_or`): with `||` a site-wise true side selects all samples; with
/// `|` only the samples that are true.
pub(super) fn logic_or(a: &Token, b: &Token, r: &mut Token) {
    init_samples(a, b, r);
    if !a.pass_site && !b.pass_site {
        return;
    }
    r.pass_site = true;
    if a.nsamples == 0 && b.nsamples == 0 {
        return;
    }
    let one_sided = a.nsamples == 0 || b.nsamples == 0;
    let sample_tok = if a.nsamples > 0 { a } else { b };
    if r.op == Op::OrVec {
        // These conditions keep QUAL>30 || FMT/GQ>30 from selecting all samples at sites with
        // QUAL<=30.
        if (a.nsamples == 0 && !a.pass_site) || (b.nsamples == 0 && !b.pass_site) {
            copy_used(r, |i| {
                sample_tok.pass_samples.get(i).copied().unwrap_or(false)
            });
        } else {
            copy_used(r, |_| true);
        }
        return;
    }
    if one_sided {
        copy_used(r, |i| {
            sample_tok.pass_samples.get(i).copied().unwrap_or(false)
        });
    } else {
        copy_used(r, |i| a.pass_samples[i] | b.pass_samples[i]);
    }
}

/// `&&` and `&` (`vector_logic_and`): with `&&` the sides may be true in different samples;
/// with `&` both must hold in one sample.
pub(super) fn logic_and(a: &Token, b: &Token, r: &mut Token) {
    init_samples(a, b, r);
    if !a.pass_site || !b.pass_site {
        return;
    }
    if a.nsamples == 0 && b.nsamples == 0 {
        r.pass_site = true;
        return;
    }
    if a.nsamples == 0 || b.nsamples == 0 {
        let tok = if a.nsamples > 0 { a } else { b };
        copy_used(r, |i| tok.pass_samples.get(i).copied().unwrap_or(false));
        r.pass_site = true;
        return;
    }
    if r.op == Op::AndVec {
        copy_used(r, |i| a.pass_samples[i] | b.pass_samples[i]);
        r.pass_site = true;
    } else {
        copy_used(r, |i| a.pass_samples[i] & b.pass_samples[i]);
        if (0..r.nsamples).any(|i| r.uses(i) && r.pass_samples[i]) {
            r.pass_site = true;
        }
    }
}

/// Sets the result of each used sample.
fn copy_used(r: &mut Token, f: impl Fn(usize) -> bool) {
    for i in 0..r.nsamples {
        if r.uses(i) {
            r.pass_samples[i] = f(i);
        }
    }
}
