//! Functions of the expression language (bcftools' `func_*`).

use anyhow::bail;
use pasta::variant::{Header, Record, encoding::INT32_VECTOR_END};

use super::{
    compare,
    eval::Context,
    kfunc,
    token::{Fraction, Func, Token, expand, is_absent, is_missing, missing},
};

/// `BCF_VL_R`.
const VL_R: u32 = 4;

fn lossy(s: &[u8]) -> String {
    String::from_utf8_lossy(s).into_owned()
}

/// The arguments of a function: the tokens on the stack.
pub(super) struct Args<'a> {
    pub(super) tokens: &'a [Token],
    pub(super) stack: &'a [usize],
}

impl Args<'_> {
    /// The token `k` places from the top of the stack (0 is the top).
    fn top(&self, k: usize) -> Option<&Token> {
        let n = self.stack.len();
        (k < n).then(|| &self.tokens[self.stack[n - 1 - k]])
    }

    /// The first of the last `nargs` tokens, and the second if there are two.
    fn first_of(&self, nargs: usize) -> Option<(&Token, Option<&Token>)> {
        let n = self.stack.len();
        let istack = n.checked_sub(nargs)?;
        let first = &self.tokens[*self.stack.get(istack)?];
        let second = (istack + 2 == n).then(|| &self.tokens[self.stack[istack + 1]]);
        Some((first, second))
    }
}

/// Applies a function token, returning how many stack entries it consumed.
pub(super) fn call(
    func: Func,
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    rtok: &mut Token,
    args: &Args<'_>,
    expr: &str,
) -> anyhow::Result<usize> {
    if func == Func::LogicOr || func == Func::LogicAnd {
        let (Some(b), Some(a)) = (args.top(0), args.top(1)) else {
            if func == Func::LogicAnd {
                bail!(
                    "Error occurred while processing the filter \"{expr}\". (nstack={})",
                    args.stack.len()
                );
            }
            bail!("Error occurred while processing the filter \"{expr}\"");
        };
        if func == Func::LogicOr {
            compare::logic_or(a, b, rtok);
        } else {
            compare::logic_and(a, b, rtok);
        }
        return Ok(2);
    }
    if matches!(func, Func::Binom | Func::Fisher) {
        return if func == Func::Binom {
            binom(ctx, header, rec, rtok, args)
        } else {
            fisher(ctx, header, rec, rtok, args)
        };
    }
    if func == Func::Perl {
        bail!(
            "\nPerl filtering requires running `configure --enable-perl-filters` at compile time.\n\n"
        );
    }
    let Some(tok) = args.top(0) else {
        bail!("Error parsing the expression");
    };
    match func {
        Func::Pass(fraction) => npass(rtok, tok, fraction)?,
        Func::Max => reduce(rtok, tok, Reduce::Max),
        Func::Min => reduce(rtok, tok, Reduce::Min),
        Func::Avg => reduce(rtok, tok, Reduce::Avg),
        Func::Sum => reduce(rtok, tok, Reduce::Sum),
        Func::Median => sorted_stat(ctx, rtok, tok, Reduce::Median),
        Func::Stdev => sorted_stat(ctx, rtok, tok, Reduce::Stdev),
        Func::SMax => per_sample(ctx, rtok, tok, Reduce::Max),
        Func::SMin => per_sample(ctx, rtok, tok, Reduce::Min),
        Func::SAvg => per_sample(ctx, rtok, tok, Reduce::Avg),
        Func::SSum => per_sample(ctx, rtok, tok, Reduce::Sum),
        Func::SMedian => per_sample(ctx, rtok, tok, Reduce::Median),
        Func::SStdev => per_sample(ctx, rtok, tok, Reduce::Stdev),
        Func::Abs => map_values(rtok, tok, f64::abs, "ABS")?,
        Func::Phred => {
            rtok.nsamples = tok.nsamples;
            let n = tok.nsamples.min(tok.pass_samples.len());
            if rtok.pass_samples.len() < n {
                rtok.pass_samples.resize(n, false);
            }
            rtok.pass_samples[..n].copy_from_slice(&tok.pass_samples[..n]);
            map_values(rtok, tok, |v| -4.342_944_819_03 * v.ln(), "PHRED")?;
        }
        Func::Count => count(rtok, tok)?,
        Func::SCount => {
            if tok.nsamples == 0 {
                count(rtok, tok)?;
            } else {
                if tok.is_str {
                    bail!("todo: Type=String for COUNT on FORMAT fields?");
                }
                start_per_sample(rtok, tok);
                for i in 0..tok.nsamples {
                    if rtok.uses(i) {
                        let cnt = (0..tok.nval1)
                            .filter(|&j| !is_absent(tok.sample_value(i, j)))
                            .count();
                        rtok.values[i] = cnt as f64;
                    }
                }
            }
        }
        Func::Strlen => strlen(rtok, tok),
        Func::LogicOr | Func::LogicAnd | Func::Binom | Func::Fisher | Func::Perl => {
            unreachable!("handled above")
        }
    }
    Ok(1)
}

/// `N_PASS()` and `F_PASS()`: the number or fraction of used samples that pass.
fn npass(rtok: &mut Token, tok: &Token, fraction: Fraction) -> anyhow::Result<()> {
    if tok.nsamples == 0 {
        bail!(
            "The function {} works with FORMAT fields",
            lossy(rtok.tag.as_deref().unwrap_or_default())
        );
    }
    let mut nsmpl = 0u32;
    let mut npass = 0u32;
    for i in 0..tok.nsamples {
        if !tok.uses(i) {
            continue;
        }
        nsmpl += 1;
        if tok.pass_samples.get(i).copied().unwrap_or(false) {
            npass += 1;
        }
    }
    expand(&mut rtok.values, 1);
    rtok.nsamples = 0;
    rtok.nvalues = 1;
    rtok.values[0] = match fraction {
        Fraction::Count => f64::from(npass),
        Fraction::Fraction if nsmpl > 0 => f64::from(npass) / f64::from(nsmpl),
        Fraction::Fraction => 0.0,
    };
    Ok(())
}

/// A reduction of a set of values.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Reduce {
    Max,
    Min,
    Avg,
    Sum,
    Median,
    Stdev,
}

impl Reduce {
    /// Reduces `vals` (non-missing values; reordered for the median), or `None` if empty.
    fn apply(self, vals: &mut [f64]) -> Option<f64> {
        let n = vals.len();
        if n == 0 {
            return None;
        }
        Some(match self {
            Self::Max => vals
                .iter()
                .fold(f64::NEG_INFINITY, |m, &v| if m < v { v } else { m }),
            Self::Min => vals
                .iter()
                .fold(f64::INFINITY, |m, &v| if m > v { v } else { m }),
            Self::Sum => vals.iter().sum(),
            Self::Avg => vals.iter().sum::<f64>() / n as f64,
            Self::Median => {
                if n == 1 {
                    vals[0]
                } else {
                    vals.sort_unstable_by(|a, b| {
                        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
                    });
                    if n % 2 == 1 {
                        vals[n / 2]
                    } else {
                        (vals[n / 2 - 1] + vals[n / 2]) * 0.5
                    }
                }
            }
            Self::Stdev => {
                if n == 1 {
                    0.0
                } else {
                    let avg = vals.iter().sum::<f64>() / n as f64;
                    let sdev: f64 = vals.iter().map(|v| (v - avg) * (v - avg)).sum();
                    (sdev / n as f64).sqrt()
                }
            }
        })
    }
}

/// Collects the non-missing values of a token (of its used samples) into `out`.
fn collect_values(tok: &Token, out: &mut Vec<f64>) {
    out.clear();
    if tok.nsamples > 0 {
        for i in 0..tok.nsamples {
            if !tok.uses(i) {
                continue;
            }
            out.extend(
                (0..tok.nval1)
                    .map(|j| tok.sample_value(i, j))
                    .filter(|&v| !is_absent(v)),
            );
        }
    } else {
        out.extend(
            (0..tok.nvalues)
                .map(|i| tok.value(i))
                .filter(|&v| !is_absent(v)),
        );
    }
}

/// MAX, MIN, AVG/MEAN, and SUM over all values (`func_max` and friends).
fn reduce(rtok: &mut Token, tok: &Token, how: Reduce) {
    rtok.nvalues = 0;
    if tok.nvalues == 0 {
        return;
    }
    let (mut max, mut min, mut sum, mut n) = (f64::NEG_INFINITY, f64::INFINITY, 0.0, 0usize);
    let mut visit = |v: f64| {
        if is_absent(v) {
            return;
        }
        n += 1;
        sum += v;
        if max < v {
            max = v;
        }
        if min > v {
            min = v;
        }
    };
    if tok.nsamples > 0 {
        for i in 0..tok.nsamples {
            if tok.uses(i) {
                for j in 0..tok.nval1 {
                    visit(tok.sample_value(i, j));
                }
            }
        }
    } else {
        for i in 0..tok.nvalues {
            visit(tok.value(i));
        }
    }
    if n == 0 {
        return;
    }
    rtok.values[0] = match how {
        Reduce::Max => max,
        Reduce::Min => min,
        Reduce::Sum => sum,
        _ => sum / n as f64,
    };
    rtok.nvalues = 1;
}

/// MEDIAN and STDEV over all values (`func_median`, `func_stddev`).
fn sorted_stat(ctx: &mut Context, rtok: &mut Token, tok: &Token, how: Reduce) {
    rtok.nvalues = 0;
    if tok.nvalues == 0 {
        return;
    }
    collect_values(tok, &mut ctx.scratch);
    if let Some(v) = how.apply(&mut ctx.scratch) {
        rtok.values[0] = v;
        rtok.nvalues = 1;
    }
}

/// Makes `rtok` one value per sample, using the samples of `tok`.
fn start_per_sample(rtok: &mut Token, tok: &Token) {
    rtok.nsamples = tok.nsamples;
    rtok.nvalues = tok.nsamples;
    rtok.nval1 = 1;
    expand(&mut rtok.values, rtok.nvalues);
    rtok.usmpl.clear();
    rtok.usmpl.extend_from_slice(&tok.usmpl);
    rtok.usmpl.resize(tok.nsamples, false);
}

/// The SMPL_* functions: one result per sample (`func_smpl_max` and friends).
fn per_sample(ctx: &mut Context, rtok: &mut Token, tok: &Token, how: Reduce) {
    if tok.nsamples == 0 {
        match how {
            Reduce::Median | Reduce::Stdev => sorted_stat(ctx, rtok, tok, how),
            _ => reduce(rtok, tok, how),
        }
        return;
    }
    start_per_sample(rtok, tok);
    for i in 0..tok.nsamples {
        if !rtok.uses(i) {
            continue;
        }
        ctx.scratch.clear();
        if tok.nvalues > 0 {
            ctx.scratch.extend(
                (0..tok.nval1)
                    .map(|j| tok.sample_value(i, j))
                    .filter(|&v| !is_absent(v)),
            );
        }
        rtok.values[i] = how.apply(&mut ctx.scratch).unwrap_or_else(missing);
    }
}

/// ABS and PHRED: a function of each value, keeping the shape (`func_abs`, `func_phred`).
fn map_values(
    rtok: &mut Token,
    tok: &Token,
    f: impl Fn(f64) -> f64,
    name: &str,
) -> anyhow::Result<()> {
    if tok.is_str {
        bail!("{name}() can be applied only on numeric values");
    }
    rtok.nsamples = tok.nsamples;
    rtok.nvalues = tok.nvalues;
    rtok.nval1 = tok.nval1;
    expand(&mut rtok.values, rtok.nvalues);
    if tok.nsamples > 0 {
        rtok.usmpl.clear();
        rtok.usmpl.extend_from_slice(&tok.usmpl);
        rtok.usmpl.resize(tok.nsamples, false);
    }
    if tok.nvalues == 0 {
        return Ok(());
    }
    if tok.nsamples > 0 {
        expand(&mut rtok.values, tok.nsamples * tok.nval1);
        for i in 0..tok.nsamples {
            if !tok.uses(i) {
                continue;
            }
            for j in 0..tok.nval1 {
                let v = tok.sample_value(i, j);
                rtok.values[i * tok.nval1 + j] = if is_absent(v) { missing() } else { f(v) };
            }
        }
    } else {
        for i in 0..tok.nvalues {
            let v = tok.value(i);
            rtok.values[i] = if is_missing(v) || (is_absent(v) && name == "PHRED") {
                missing()
            } else if is_absent(v) {
                v
            } else {
                f(v)
            };
        }
    }
    Ok(())
}

/// COUNT: the number of values, of passing samples, or of comma-separated strings
/// (`func_count`).
fn count(rtok: &mut Token, tok: &Token) -> anyhow::Result<()> {
    let cnt = if tok.tag.is_some() && tok.nsamples > 0 {
        if tok.is_str {
            bail!("todo: Type=String for COUNT on FORMAT fields?");
        }
        (0..tok.nsamples)
            .filter(|&i| tok.uses(i))
            .map(|i| {
                (0..tok.nval1)
                    .filter(|&j| !is_absent(tok.sample_value(i, j)))
                    .count()
            })
            .sum()
    } else if tok.nsamples > 0 {
        tok.pass_samples[..tok.nsamples.min(tok.pass_samples.len())]
            .iter()
            .filter(|&&p| p)
            .count()
    } else if tok.is_str {
        let s = &tok.str_value;
        usize::from(!s.is_empty()) + s.iter().filter(|&&c| c == b',').count()
    } else {
        tok.nvalues
    };
    rtok.nvalues = 1;
    rtok.values[0] = cnt as f64;
    Ok(())
}

/// STRLEN: the length of a string, or of each comma-separated string of a list
/// (`func_strlen`).
fn strlen(rtok: &mut Token, tok: &Token) {
    rtok.nvalues = 0;
    rtok.str_value.clear();
    if tok.str_value.is_empty() {
        return;
    }
    let s = tok.c_str();
    if tok.idx == -2 {
        let mut i = 0;
        let mut pos = 0;
        while pos < s.len() {
            let end = s[pos..]
                .iter()
                .position(|&c| c == b',')
                .map_or(s.len(), |k| pos + k);
            expand(&mut rtok.values, i + 1);
            rtok.values[i] = (end - pos) as f64;
            pos = if end < s.len() { end + 1 } else { end };
            i += 1;
        }
        rtok.nvalues = i;
    } else {
        rtok.values[0] = if s == b"." { 0.0 } else { s.len() as f64 };
        rtok.nvalues = 1;
    }
}

/// The two leading alleles of each sample's diploid genotype in `ctx.tmpi`, or `None` per
/// sample when missing or haploid. Fails on an allele index past the alleles.
fn diploid_alleles(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
) -> anyhow::Result<Option<usize>> {
    let ngt = rec.get_genotypes(header, &mut ctx.tmpi).unwrap_or(0);
    let nsample = rec.sample_count();
    let max_ploidy = ngt.checked_div(nsample).unwrap_or(0);
    if ngt == 0 || max_ploidy < 2 {
        return Ok(None);
    }
    Ok(Some(max_ploidy))
}

/// The allele indices of sample `i`'s first two GT values, or `None` when either is missing
/// or the sample is haploid.
fn sample_alleles(
    ctx: &Context,
    header: &Header,
    rec: &Record,
    max_ploidy: usize,
    i: usize,
) -> anyhow::Result<Option<(usize, usize)>> {
    let p = &ctx.tmpi[i * max_ploidy..];
    if p[0] >> 1 == 0 || p[1] >> 1 == 0 || p[1] == INT32_VECTOR_END {
        return Ok(None);
    }
    let n_allele = rec.allele_count() as i32;
    let (a1, a2) = ((p[0] >> 1) - 1, (p[1] >> 1) - 1);
    if a1 >= n_allele || a2 >= n_allele || a1 < 0 || a2 < 0 {
        bail!(
            "Incorrect allele index at {}:{}, sample {}",
            lossy(header.contig_name(rec.rid()).unwrap_or(b"(null)")),
            rec.pos() + 1,
            lossy(&header.samples()[i])
        );
    }
    Ok(Some((a1 as usize, a2 as usize)))
}

/// Value `k` of sample `i`, or missing past the sample's values.
fn sample_at(tok: &Token, i: usize, k: usize) -> f64 {
    if k < tok.nval1 && tok.nvalues > 0 {
        tok.sample_value(i, k)
    } else {
        missing()
    }
}

/// The two-sided binomial test of a pair of counts, missing when either is missing.
fn binom_value(x: f64, y: f64) -> f64 {
    if is_absent(x) || is_absent(y) {
        return missing();
    }
    let p = kfunc::binom_two_sided(x as i32, y as i32, 0.5);
    if p < 0.0 { missing() } else { p }
}

/// `binom(AD)` (the alleles taken from GT) or `binom(AD[:0],AD[:1])` (`func_binom`).
fn binom(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    rtok: &mut Token,
    args: &Args<'_>,
) -> anyhow::Result<usize> {
    let nargs = rtok.nargs;
    if nargs != 2 && nargs != 1 {
        bail!("Error: binom() takes one or two arguments");
    }
    let Some((tok, tok2)) = args.first_of(nargs) else {
        bail!("Error parsing the expression");
    };
    if tok.nsamples > 0 {
        rtok.nval1 = 1;
        rtok.nvalues = tok.nsamples;
        rtok.nsamples = tok.nsamples;
        expand(&mut rtok.values, rtok.nvalues);
        rtok.usmpl.clear();
        rtok.usmpl.extend_from_slice(&tok.usmpl);
        rtok.usmpl.resize(tok.nsamples, false);
        if nargs == 1 {
            let Some(max_ploidy) = diploid_alleles(ctx, header, rec)? else {
                for i in 0..rtok.nsamples {
                    if rtok.uses(i) {
                        rtok.values[i] = missing();
                    }
                }
                return Ok(nargs);
            };
            for i in 0..rtok.nsamples {
                if !rtok.uses(i) {
                    continue;
                }
                rtok.values[i] = match sample_alleles(ctx, header, rec, max_ploidy, i)? {
                    Some((a1, a2)) => binom_value(sample_at(tok, i, a1), sample_at(tok, i, a2)),
                    None => missing(),
                };
            }
        } else {
            let tok2 = tok2.expect("two arguments");
            if tok.nval1 != 1 || tok2.nval1 != 1 {
                bail!(
                    "Expected one value per binom() argument, found {} and {} at {}:{}",
                    tok.nval1,
                    tok2.nval1,
                    lossy(header.contig_name(rec.rid()).unwrap_or(b"(null)")),
                    rec.pos() + 1
                );
            }
            for i in 0..rtok.nsamples {
                if rtok.uses(i) {
                    rtok.values[i] = binom_value(sample_at(tok, i, 0), sample_at(tok2, i, 0));
                }
            }
        }
    } else {
        rtok.nvalues = 1;
        expand(&mut rtok.values, 1);
        let pair = match tok2 {
            None if tok.nvalues == 2 => Some((tok.value(0), tok.value(1))),
            Some(t2) if tok.nvalues == 1 && t2.nvalues == 1 => Some((tok.value(0), t2.value(0))),
            _ => None,
        };
        rtok.values[0] = pair.map_or_else(missing, |(x, y)| binom_value(x, y));
    }
    Ok(nargs)
}

/// Fisher's exact test of a 2x2 table, missing when any count is missing.
fn fisher_value(n: [f64; 4]) -> f64 {
    if n.iter().any(|&v| is_absent(v)) {
        return missing();
    }
    kfunc::fisher_exact_two(n[0] as i32, n[1] as i32, n[2] as i32, n[3] as i32)
}

/// `fisher(INFO/DP4)`, `fisher(ADF[...],ADR[...])`, and per-sample forms (`func_fisher`).
fn fisher(
    ctx: &mut Context,
    header: &Header,
    rec: &mut Record,
    rtok: &mut Token,
    args: &Args<'_>,
) -> anyhow::Result<usize> {
    let nargs = rtok.nargs;
    let Some((tok, tok2)) = args.first_of(nargs) else {
        bail!("Error parsing the expression");
    };
    if tok.nsamples == 0 {
        rtok.nvalues = 1;
        expand(&mut rtok.values, 1);
        let table = match (nargs, tok2) {
            (1, _) if tok.nvalues == 4 => {
                Some([tok.value(0), tok.value(1), tok.value(2), tok.value(3)])
            }
            (2, Some(t2)) if tok.nvalues >= 2 && t2.nvalues >= 2 => {
                Some([tok.value(0), t2.value(0), tok.value(1), t2.value(1)])
            }
            _ => None,
        };
        rtok.values[0] = table.map_or_else(missing, fisher_value);
        return Ok(nargs);
    }
    rtok.nval1 = 1;
    rtok.nvalues = tok.nsamples;
    rtok.nsamples = tok.nsamples;
    expand(&mut rtok.values, rtok.nvalues);
    rtok.usmpl.clear();
    rtok.usmpl.extend_from_slice(&tok.usmpl);
    rtok.usmpl.resize(tok.nsamples, false);
    let both_r = tok2.is_some_and(|t2| tok.vl_len == VL_R && t2.vl_len == VL_R);
    if nargs == 1 && tok.nval1 == 4 {
        for i in 0..rtok.nsamples {
            if rtok.uses(i) {
                rtok.values[i] = fisher_value([0, 1, 2, 3].map(|k| sample_at(tok, i, k)));
            }
        }
    } else if nargs == 2 && both_r && tok.nuidxs + tok2.map_or(0, |t| t.nuidxs) == 4 {
        let t2 = tok2.expect("two arguments");
        for i in 0..rtok.nsamples {
            if rtok.uses(i) {
                rtok.values[i] = fisher_value([
                    sample_at(tok, i, 0),
                    sample_at(tok, i, 1),
                    sample_at(t2, i, 0),
                    sample_at(t2, i, 1),
                ]);
            }
        }
    } else if nargs == 2 && both_r {
        let t2 = tok2.expect("two arguments");
        let Some(max_ploidy) = diploid_alleles(ctx, header, rec)? else {
            for i in 0..rtok.nsamples {
                if rtok.uses(i) {
                    rtok.values[i] = missing();
                }
            }
            return Ok(nargs);
        };
        for i in 0..rtok.nsamples {
            if !rtok.uses(i) {
                continue;
            }
            rtok.values[i] = match sample_alleles(ctx, header, rec, max_ploidy, i)? {
                Some((a1, a2)) => fisher_value([
                    sample_at(tok, i, a1),
                    sample_at(tok, i, a2),
                    sample_at(t2, i, a1),
                    sample_at(t2, i, a2),
                ]),
                None => missing(),
            };
        }
    } else {
        for i in 0..rtok.nsamples {
            if rtok.uses(i) {
                rtok.values[i] = missing();
            }
        }
    }
    Ok(nargs)
}
