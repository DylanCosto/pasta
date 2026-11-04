//! Evaluating a compiled expression on a record (bcftools' `filter_test`).

use anyhow::bail;
use pasta::variant::{Header, Record};

use super::{
    ExtType, ExtValue, Filter, FilterStatus, compare,
    compare::{Arith, Cmp},
    funcs::{self, Args},
    setters,
    token::{Comparator, Op, Token, expand},
};

/// Whether a per-record cache holds the current record's data.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum Cache {
    /// Not computed for this record.
    #[default]
    Stale,
    Valid,
    Failed,
}

/// Each sample's GT alleles as a bit mask (bcftools' `cached_GT`).
#[derive(Debug, Default)]
pub(super) struct GtCache {
    pub(super) state: Cache,
    pub(super) buf: Vec<i32>,
    pub(super) mask: Vec<u64>,
}

/// Allele counts and frequencies of a site (bcftools' `cached_site`).
#[derive(Debug, Default)]
pub(super) struct SiteCache {
    pub(super) state: Cache,
    pub(super) n_missing: i64,
    pub(super) n_samples: usize,
    pub(super) an: i32,
    /// Counts by allele, REF first.
    pub(super) ac: Vec<i32>,
    pub(super) af: Vec<f64>,
}

/// Scratch buffers and caches reused across records.
#[derive(Debug, Default)]
pub(super) struct Context {
    pub(super) tmpi: Vec<i32>,
    pub(super) tmpf: Vec<f32>,
    pub(super) strbuf: Vec<u8>,
    pub(super) scratch: Vec<f64>,
    pub(super) gt: GtCache,
    pub(super) site: SiteCache,
    /// The samples to use for values computed over samples (`filter_set_samples`).
    pub(super) sample_mask: Vec<bool>,
    /// The number of samples in the mask, 0 when no mask is set.
    pub(super) sample_mask_set: usize,
}

impl Context {
    pub(super) fn new(nsamples: usize) -> Self {
        Self {
            sample_mask: vec![false; nsamples],
            ..Self::default()
        }
    }
}

impl Filter {
    /// Sets the external values of `{}` tokens (`filter_test_ext`).
    pub(super) fn set_ext(&mut self, ext: &[ExtValue<'_>]) {
        for tok in &mut self.tokens {
            let Some(j) = tok.iext.checked_sub(1) else {
                continue;
            };
            let ty = self.ext.get(j).copied().flatten();
            match ext.get(j).copied().unwrap_or(ExtValue::Missing) {
                ExtValue::Missing => {
                    tok.is_missing = true;
                    tok.nvalues = 0;
                    if ty == Some(ExtType::Str) {
                        tok.str_value.clear();
                    }
                }
                value => {
                    tok.is_missing = false;
                    tok.nvalues = 1;
                    expand(&mut tok.values, 1);
                    match (ty, value) {
                        (Some(ExtType::Str), ExtValue::Str(s)) => {
                            tok.str_value.clear();
                            tok.str_value.extend_from_slice(s);
                            tok.nvalues = s.len();
                        }
                        (_, ExtValue::Int(v)) => tok.values[0] = f64::from(v),
                        (_, ExtValue::Float(v)) => tok.values[0] = f64::from(v),
                        _ => {}
                    }
                }
            }
        }
    }

    /// Evaluates the expression on a record (`filter_test`).
    pub(super) fn evaluate(&mut self, header: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        if self.status != FilterStatus::OK {
            bail!("Error: the caller did not check the filter status");
        }
        if !self.eval_unpack.is_empty() {
            rec.unpack(self.eval_unpack);
        }
        self.ctx.gt.state = Cache::Stale;
        self.ctx.site.state = Cache::Stale;
        let Self {
            expr,
            tokens,
            stack,
            ctx,
            nsamples,
            ..
        } = self;
        let nsamples = *nsamples;
        stack.clear();
        for i in 0..tokens.len() {
            let (before, rest) = tokens.split_at_mut(i);
            let tok = &mut rest[0];
            tok.pass_site = false;
            if tok.op == Op::Val {
                if let Some(setter) = tok.setter {
                    setters::apply(setter, ctx, header, rec, tok)?;
                } else if let Some(key) = &tok.key {
                    tok.str_value.clear();
                    tok.str_value.extend_from_slice(key);
                    tok.nvalues = key.len();
                } else if tok.is_constant {
                    tok.values[0] = tok.threshold;
                    tok.nvalues = 1;
                }
                // String tokens count characters as values; keep reads of them in bounds.
                expand(&mut tok.values, tok.nvalues);
                stack.push(i);
                continue;
            }
            if let Some(func) = tok.func {
                let args = Args {
                    tokens: before,
                    stack,
                };
                let nargs = funcs::call(func, ctx, header, rec, tok, &args, expr)?;
                let n = stack.len();
                if nargs == 0 || nargs > n {
                    bail!("Error occurred while processing the filter \"{expr}\"");
                }
                stack[n - nargs] = i;
                stack.truncate(n - nargs + 1);
                continue;
            }
            let n = stack.len();
            if n < 2 {
                bail!("Error occurred while processing the filter \"{expr}\" (1:{n})");
            }
            let (ia, ib) = (stack[n - 2], stack[n - 1]);
            binary(tok, before, ia, ib, rec)?;
            stack[n - 2] = i;
            stack.pop();
        }
        if stack.len() > 1 {
            bail!(
                "Error occurred while processing the filter \"{expr}\" (2:{})",
                stack.len()
            );
        }
        let Some(&top) = stack.first() else {
            bail!("Error occurred while processing the filter \"{expr}\"");
        };
        let top = &mut tokens[top];
        if nsamples > 0 && top.nsamples == 0 {
            let pass = top.pass_site;
            top.pass_samples.clear();
            top.pass_samples.resize(nsamples, pass);
        }
        Ok(top.pass_site)
    }

    /// The per-sample results of the last evaluation, when the expression reads FORMAT fields.
    pub(super) fn result_samples(&self) -> Option<&[bool]> {
        if self.nsamples == 0 {
            return None;
        }
        let top = self.tokens.get(*self.stack.first()?)?;
        Some(&top.pass_samples[..self.nsamples.min(top.pass_samples.len())])
    }

    /// The values of the last evaluation's result (`filter_get_doubles`).
    pub(super) fn result_values(&self) -> Option<(&[f64], usize)> {
        let top = self.tokens.get(*self.stack.first()?)?;
        (top.nvalues > 0).then(|| (&top.values[..top.nvalues], top.nval1))
    }

    /// Restricts sample-wise values to the samples set in `samples` (`filter_set_samples`).
    pub(super) fn restrict_samples(&mut self, samples: &[bool]) {
        self.ctx.sample_mask_set = 0;
        if self.nsamples == 0 {
            return;
        }
        for tok in &mut self.tokens {
            if tok.nsamples == 0 {
                continue;
            }
            for (u, &s) in tok.usmpl.iter_mut().zip(samples) {
                *u = s;
            }
        }
        for (m, &s) in self.ctx.sample_mask.iter_mut().zip(samples) {
            *m = s;
            if s {
                self.ctx.sample_mask_set += 1;
            }
        }
    }
}

/// Applies a binary operator token to the two tokens on top of the stack.
fn binary(
    tok: &mut Token,
    before: &[Token],
    ia: usize,
    ib: usize,
    rec: &Record,
) -> anyhow::Result<()> {
    let (a, b) = (&before[ia], &before[ib]);
    let arith = match tok.op {
        Op::Add => Some(Arith::Add),
        Op::Sub => Some(Arith::Sub),
        Op::Mult => Some(Arith::Mult),
        Op::Div => Some(Arith::Div),
        Op::Modulo => Some(Arith::Modulo),
        _ => None,
    };
    if let Some(op) = arith {
        return compare::arithmetic(a, b, tok, op);
    }

    // The comparator of the operator, else of the top operand, else of the one below.
    let comparator = if let Some(c) = tok.comparator {
        Some((c, b, a))
    } else if let Some(c) = b.comparator {
        Some((c, b, a))
    } else {
        a.comparator.map(|c| (c, a, b))
    };
    if let Some((c, x, y)) = comparator {
        return match c {
            Comparator::BitAnd => {
                compare::bit_and(x, y, tok);
                Ok(())
            }
            Comparator::Filter => compare::compare_filter(y, tok, rec),
            Comparator::StringHash => compare::compare_hash(x, y, tok),
            Comparator::Id => compare::compare_id(x, y, tok, rec),
        };
    }
    match u8::from(a.is_str) + u8::from(b.is_str) {
        2 => compare::compare_strings(a, b, tok),
        1 => bail!("Error: cannot use arithmetic operators to compare strings and numbers"),
        _ => {
            let op = match tok.op {
                Op::Eq => Cmp::Eq,
                Op::Ne => Cmp::Ne,
                Op::Le => Cmp::Le,
                Op::Lt => Cmp::Lt,
                Op::Bt => Cmp::Bt,
                Op::Be => Cmp::Be,
                other => bail!("todo: filter.c:4469 .. type={}", other.code()),
            };
            compare::compare_numbers(a, b, tok, op)
        }
    }
}
