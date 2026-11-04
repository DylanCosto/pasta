//! bcftools filter expressions (`-i/--include`, `-e/--exclude`), ported from bcftools 1.24's
//! `filter.c` without its Perl support.
//!
//! An expression is compiled into reverse Polish notation against a header ([`Filter::init`],
//! bcftools' `filter_init`) and evaluated per record ([`Filter::test`], `filter_test`). An
//! expression on FORMAT fields also decides per sample ([`Filter::sample_pass`], `filter_test`'s
//! `samples` output), which `filter -S`, `query`, and `+setGT` use.

pub(crate) mod bre;
mod compare;
mod eval;
mod funcs;
pub(crate) mod kfunc;
mod lexer;
mod parse;
mod setters;
#[cfg(test)]
mod tests;
mod token;

use std::ops::BitOrAssign;

use pasta::variant::{Header, Record, Unpack, ValueType};

/// A compiled filter expression (bcftools' `filter_t`).
pub trait RecordFilter {
    /// Evaluates the expression on a record (bcftools' `filter_test`); `true` when it holds.
    fn pass(&mut self, header: &Header, rec: &mut Record) -> anyhow::Result<bool>;

    /// The per-sample results of the last [`Self::pass`], when the expression has sample
    /// terms (bcftools' `smpl_pass`).
    fn sample_mask(&self) -> Option<&[bool]> {
        None
    }

    /// How much of a record the expression reads (bcftools' `filter_max_unpack`).
    fn max_unpack(&self) -> Unpack {
        Unpack::ALL
    }

    /// The values of the expression for the last [`Self::pass`], one per sample for sample
    /// expressions (bcftools' `filter_get_doubles`), which `query -f '%MAX(FMT/DP)'` prints.
    fn doubles(&self) -> &[f64] {
        &[]
    }
}

/// Whether records matching the expression are kept or dropped.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FilterLogic {
    /// `-i`: keep records for which the expression holds.
    Include,
    /// `-e`: drop records for which the expression holds.
    Exclude,
}

impl FilterLogic {
    /// The single `-i`/`-e` expression of a command that takes at most one, with its logic.
    ///
    /// # Errors
    ///
    /// Fails, with bcftools' message, when more than one expression is given.
    pub fn single<'a>(
        include: &'a [String],
        exclude: &'a [String],
    ) -> anyhow::Result<Option<(&'a str, Self)>> {
        if include.len() + exclude.len() > 1 {
            anyhow::bail!(
                "Error: only one -i or -e expression can be given, and they cannot be combined"
            );
        }
        Ok(include
            .first()
            .map(|e| (e.as_str(), Self::Include))
            .or_else(|| exclude.first().map(|e| (e.as_str(), Self::Exclude))))
    }
}

/// Compiles a filter expression against a header (bcftools' `filter_init`).
pub fn init(header: &Header, expression: &str) -> anyhow::Result<Box<dyn RecordFilter + Send>> {
    Ok(Box::new(Filter::init(header, expression)?))
}

/// Whether compilation stops at the first problem or records unknown tags (bcftools'
/// `filter_init` and `filter_parse`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ErrorMode {
    Exit,
    Collect,
}

/// Problems found by [`Filter::parse`] (bcftools' `FILTER_ERR_*`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct FilterStatus(u8);

impl FilterStatus {
    /// No problem (`FILTER_OK`).
    pub const OK: Self = Self(0);
    /// A name is neither a tag of the header nor a number (`FILTER_ERR_UNKN_TAGS`).
    pub const UNKNOWN_TAGS: Self = Self(1);
    /// The expression could not be parsed (`FILTER_ERR_OTHER`).
    pub const OTHER: Self = Self(2);
}

impl BitOrAssign for FilterStatus {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

/// The type of an external value, `{}`, `{str}`, `{int}`, or `{float}` (a `BCF_HT_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ExtType {
    /// A string.
    Str,
    /// An integer.
    Int,
    /// A float.
    Float,
}

impl ExtType {
    fn from_value_type(ty: ValueType) -> Self {
        match ty {
            ValueType::String => Self::Str,
            ValueType::Float => Self::Float,
            ValueType::Integer | ValueType::Flag => Self::Int,
        }
    }
}

/// The value of an external-value token for one evaluation ([`Filter::test_ext`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ExtValue<'a> {
    /// Missing (`.`).
    Missing,
    /// A string.
    Str(&'a [u8]),
    /// An integer.
    Int(i32),
    /// A float.
    Float(f32),
}

/// A compiled filter expression (bcftools' `filter_t`).
#[derive(Debug)]
pub struct Filter {
    expr: String,
    tokens: Vec<token::Token>,
    /// The evaluation stack: indices of tokens.
    stack: Vec<usize>,
    max_unpack: Unpack,
    /// What evaluating decodes: [`Self::max_unpack`] without the alleles and FORMAT when the
    /// tokens read only fixed fields and numeric FORMAT fields.
    eval_unpack: Unpack,
    /// The header's samples when the expression reads FORMAT fields, else 0.
    nsamples: usize,
    status: FilterStatus,
    undef_tags: Vec<String>,
    used_tags: Vec<String>,
    ext: Vec<Option<ExtType>>,
    ctx: eval::Context,
}

impl Filter {
    /// Compiles an expression, failing on any problem (bcftools' `filter_init`).
    ///
    /// # Errors
    ///
    /// Returns bcftools' message for a malformed expression or a name that is not a tag of
    /// the header.
    pub fn init(header: &Header, expr: &str) -> anyhow::Result<Self> {
        parse::compile(header, expr, ErrorMode::Exit)
    }

    /// Compiles an expression, recording unknown tags and unparseable expressions in
    /// [`Self::status`] instead of failing (bcftools' `filter_parse`).
    ///
    /// # Errors
    ///
    /// Returns bcftools' message for the other problems.
    pub fn parse(header: &Header, expr: &str) -> anyhow::Result<Self> {
        parse::compile(header, expr, ErrorMode::Collect)
    }

    /// The problems [`Self::parse`] found (bcftools' `filter_status`).
    pub fn status(&self) -> FilterStatus {
        self.status
    }

    /// The names that are neither tags nor numbers (bcftools' `filter_list_undef_tags`).
    pub fn undef_tags(&self) -> &[String] {
        &self.undef_tags
    }

    /// The tags and columns the expression reads, as `INFO/X`, `FORMAT/X`, or `QUAL`
    /// (bcftools' `filter_list_used_tags`).
    pub fn used_tags(&self) -> &[String] {
        &self.used_tags
    }

    /// The parts of a record the expression reads (bcftools' `filter_max_unpack`).
    pub fn max_unpack(&self) -> Unpack {
        self.max_unpack
    }

    /// The types of the external values, `None` where undetermined (bcftools'
    /// `filter_ext_types`).
    pub fn ext_types(&self) -> &[Option<ExtType>] {
        &self.ext
    }

    /// Evaluates the expression on a record (bcftools' `filter_test`).
    ///
    /// # Errors
    ///
    /// Returns bcftools' message for an expression that cannot be evaluated on the record.
    pub fn test(&mut self, header: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        self.evaluate(header, rec)
    }

    /// Evaluates the expression with values for its `{}` tokens (bcftools'
    /// `filter_test_ext`).
    ///
    /// # Errors
    ///
    /// As [`Self::test`].
    pub fn test_ext(
        &mut self,
        header: &Header,
        rec: &mut Record,
        ext: &[ExtValue<'_>],
    ) -> anyhow::Result<bool> {
        if !self.ext.is_empty() {
            self.set_ext(ext);
        }
        self.evaluate(header, rec)
    }

    /// The per-sample results of the last evaluation, `None` when the expression reads no
    /// FORMAT field (bcftools' `filter_test` samples).
    pub fn sample_pass(&self) -> Option<&[bool]> {
        self.result_samples()
    }

    /// The values of the last evaluation's result and their number per sample (bcftools'
    /// `filter_get_doubles`). Missing values are bcftools' NaN sentinels.
    pub fn doubles(&self) -> Option<(&[f64], usize)> {
        self.result_values()
    }

    /// Restricts the samples that sample-wise values (and computed AC/AF) use to those set in
    /// `samples` (bcftools' `filter_set_samples`).
    pub fn set_samples(&mut self, samples: &[bool]) {
        self.restrict_samples(samples);
    }

    /// Whether a value is bcftools' missing sentinel (`bcf_double_is_missing`).
    #[expect(
        dead_code,
        reason = "bcftools' filter API for reading `doubles`, which no command needs yet"
    )]
    pub fn is_missing(value: f64) -> bool {
        token::is_missing(value)
    }
}

impl RecordFilter for Filter {
    fn pass(&mut self, header: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        self.evaluate(header, rec)
    }

    fn sample_mask(&self) -> Option<&[bool]> {
        self.result_samples()
    }

    fn max_unpack(&self) -> Unpack {
        self.max_unpack
    }

    fn doubles(&self) -> &[f64] {
        self.result_values().map_or(&[], |(values, _)| values)
    }
}
