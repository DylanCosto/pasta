//! Filter tokens (bcftools' `token_t`) and the double sentinels for missing values.

use std::collections::HashSet;

use pasta::variant::ValueType;

use super::bre::Bre;

/// bcftools' `bcf_double_missing`.
const DOUBLE_MISSING: u64 = 0x7ff0_0000_0000_0001;
/// bcftools' `bcf_double_vector_end`.
const DOUBLE_VECTOR_END: u64 = 0x7ff0_0000_0000_0002;

/// A missing value (`bcf_double_set_missing`).
#[inline]
pub(super) fn missing() -> f64 {
    f64::from_bits(DOUBLE_MISSING)
}

/// The end of a shorter vector (`bcf_double_set_vector_end`).
#[inline]
pub(super) fn vector_end() -> f64 {
    f64::from_bits(DOUBLE_VECTOR_END)
}

/// `bcf_double_is_missing`.
#[inline]
pub(super) fn is_missing(x: f64) -> bool {
    x.to_bits() == DOUBLE_MISSING
}

/// `bcf_double_is_vector_end`.
#[inline]
pub(super) fn is_vector_end(x: f64) -> bool {
    x.to_bits() == DOUBLE_VECTOR_END
}

/// `bcf_double_is_missing_or_vector_end`.
#[inline]
pub(super) fn is_absent(x: f64) -> bool {
    let b = x.to_bits();
    b == DOUBLE_MISSING || b == DOUBLE_VECTOR_END
}

/// Grows `v` to hold at least `n` values (`hts_expand`).
#[inline]
pub(super) fn expand(v: &mut Vec<f64>, n: usize) {
    if v.len() < n {
        v.resize(n, 0.0);
    }
}

/// The kind of a token (bcftools' `TOK_*`), in bcftools' numbering, which indexes
/// [`OP_PREC`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Op {
    Val,
    Lft,
    Rgt,
    Le,
    Lt,
    Eq,
    Bt,
    Be,
    Ne,
    Or,
    And,
    Add,
    Sub,
    Mult,
    Div,
    Max,
    Min,
    Avg,
    AndVec,
    OrVec,
    Like,
    NLike,
    Sum,
    Abs,
    Len,
    Func,
    Cnt,
    PerlSub,
    Binom,
    Phred,
    Median,
    Stdev,
    SMax,
    SMin,
    SAvg,
    SMedian,
    SStdev,
    SSum,
    In,
    NotIn,
    Modulo,
    Ext,
    Fisher,
    SCount,
}

/// Operator precedence by [`Op`] (bcftools' `op_prec`).
const OP_PREC: [u8; 44] = [
    0, 1, 1, 5, 5, 5, 5, 5, 5, 2, 3, 6, 6, 7, 7, 8, 8, 8, 3, 2, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 7, 0, 8, 8,
];

impl Op {
    /// The operator's precedence.
    #[inline]
    pub(super) fn prec(self) -> u8 {
        OP_PREC[self as usize]
    }

    /// bcftools' numeric `TOK_*` code, for messages.
    pub(super) fn code(self) -> usize {
        self as usize
    }
}

/// Which genotype classes a GT token yields (`filters_set_genotype2/3/4`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum GtClass {
    /// `rr`, `ra`, `aa`, `aA`, `r`, `a`.
    Alleles,
    /// `hom`, `het`, `hap`.
    Zygosity,
    /// `mis`, `alt`, `ref`.
    Presence,
}

/// How a value token reads the record (bcftools' `setter`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Setter {
    Qual,
    Type,
    /// A `Number=1` INFO tag (`filters_set_info`).
    Info,
    InfoInt,
    InfoFloat,
    InfoString,
    InfoFlag,
    FormatInt,
    FormatFloat,
    FormatString,
    Genotype(GtClass),
    GenotypeString,
    Ilen,
    RefString,
    AltString,
    /// `N_MISSING` or `F_MISSING`.
    Missing(Fraction),
    NAlt,
    Chrom,
    Pos,
    Ac,
    An,
    Mac,
    Af,
    Maf,
}

/// A count or a fraction of samples (`N_*` or `F_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Fraction {
    Count,
    Fraction,
}

/// A function or logic operator applied to the stack (bcftools' `func`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Func {
    Pass(Fraction),
    Max,
    Min,
    Avg,
    Median,
    Stdev,
    Sum,
    Abs,
    Count,
    Strlen,
    Phred,
    Binom,
    Fisher,
    Perl,
    SCount,
    SMax,
    SMin,
    SAvg,
    SMedian,
    SStdev,
    SSum,
    LogicOr,
    LogicAnd,
}

/// A special comparison (bcftools' `comparator`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Comparator {
    /// `TYPE~` and `TYPE!~`.
    BitAnd,
    Filter,
    /// A tag compared with `@file`.
    StringHash,
    Id,
}

/// The header line kind a tag token refers to (`BCF_HL_*`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum LineKind {
    #[default]
    None,
    Filter,
    Info,
    Format,
}

/// One token of a compiled expression (bcftools' `token_t`): the fields fixed at parse time,
/// then the values of the last evaluation.
#[derive(Debug)]
pub(super) struct Token {
    pub(super) op: Op,
    pub(super) nargs: usize,
    /// A string constant, without quotes.
    pub(super) key: Option<Vec<u8>>,
    /// The tag name (or other name, for messages).
    pub(super) tag: Option<Vec<u8>>,
    pub(super) threshold: f64,
    pub(super) is_constant: bool,
    pub(super) hdr_id: i32,
    pub(super) hl_type: LineKind,
    pub(super) ht_type: Option<ValueType>,
    /// `BCF_VL_*` of the tag.
    pub(super) vl_len: u32,
    /// The vector index: >= 0 one index, -2 a list in `idxs`, -3 the GT alleles.
    pub(super) idx: i32,
    /// 1 to include an index, 0 to exclude it; a last -1 means "and all after".
    pub(super) idxs: Vec<i32>,
    pub(super) nuidxs: usize,
    /// The samples to use (FORMAT tokens and per-sample results).
    pub(super) usmpl: Vec<bool>,
    pub(super) nsamples: usize,
    pub(super) setter: Option<Setter>,
    pub(super) func: Option<Func>,
    pub(super) comparator: Option<Comparator>,
    pub(super) hash: Option<HashSet<Vec<u8>>>,
    pub(super) regex: Option<Bre>,
    /// 1-based index of an external value (`{}`), 0 for none.
    pub(super) iext: usize,

    pub(super) values: Vec<f64>,
    pub(super) str_value: Vec<u8>,
    pub(super) is_str: bool,
    pub(super) is_missing: bool,
    pub(super) pass_site: bool,
    pub(super) pass_samples: Vec<bool>,
    pub(super) nvalues: usize,
    pub(super) nval1: usize,
}

impl Token {
    /// A token of kind `op` with nothing set (`filters_init1`'s defaults).
    pub(super) fn new(op: Op) -> Self {
        Self {
            op,
            nargs: 0,
            key: None,
            tag: None,
            threshold: 0.0,
            is_constant: false,
            hdr_id: -1,
            hl_type: LineKind::None,
            ht_type: None,
            vl_len: 0,
            idx: 0,
            idxs: Vec::new(),
            nuidxs: 0,
            usmpl: Vec::new(),
            nsamples: 0,
            setter: None,
            func: None,
            comparator: None,
            hash: None,
            regex: None,
            iext: 0,
            values: Vec::new(),
            str_value: Vec::new(),
            is_str: false,
            is_missing: false,
            pass_site: false,
            pass_samples: Vec::new(),
            nvalues: 0,
            nval1: 0,
        }
    }

    /// A numeric constant.
    pub(super) fn constant(value: f64) -> Self {
        Self {
            threshold: value,
            is_constant: true,
            ..Self::new(Op::Val)
        }
    }

    /// Value `k`, or missing past the end.
    #[inline]
    pub(super) fn value(&self, k: usize) -> f64 {
        self.values.get(k).copied().unwrap_or_else(missing)
    }

    /// The values of sample `i` (`nval1` of them, missing past the end of the buffer).
    #[inline]
    pub(super) fn sample_value(&self, i: usize, j: usize) -> f64 {
        self.value(i * self.nval1 + j)
    }

    /// Whether sample `i` is used.
    #[inline]
    pub(super) fn uses(&self, i: usize) -> bool {
        self.usmpl.get(i).copied().unwrap_or(false)
    }

    /// The string value up to its first NUL (a C string).
    #[inline]
    pub(super) fn c_str(&self) -> &[u8] {
        let s = &self.str_value;
        &s[..memchr::memchr(0, s).unwrap_or(s.len())]
    }

    /// The block of sample `i` of a per-sample string value.
    #[inline]
    pub(super) fn sample_str(&self, i: usize) -> &[u8] {
        let start = (i * self.nval1).min(self.str_value.len());
        let end = (start + self.nval1).min(self.str_value.len());
        &self.str_value[start..end]
    }

    /// Selects all vector indices (`TAG[*]`).
    pub(super) fn set_all_indices(&mut self) {
        self.idxs = vec![-1];
        self.idx = -2;
    }
}
