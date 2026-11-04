//! bcftools' format-string engine (`convert.c`), shared by `query` and `convert`: a format
//! such as `%CHROM\t%POS[\t%SAMPLE=%GT]\n` is parsed once into a list of fields, then each
//! record is printed through it.
//!
//! Fields outside square brackets are printed once per record; fields inside are printed once
//! per sample. `%TAG` resolves to a fixed column (`%CHROM`, `%POS`, ...), an INFO tag outside
//! brackets, or a FORMAT tag inside them; `{N}` picks one value of a vector.

mod hap;
mod kfunc;
mod variantkey;

use std::io::Write;

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Record, Unpack, VariantType,
    encoding::{
        FLOAT_MISSING, FLOAT_VECTOR_END, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING,
        INT16_VECTOR_END, INT32_MISSING, INT32_VECTOR_END, le_i16, le_i32, le_u32,
    },
    number::{push_double, push_int},
    vcf::{format_array, format_genotype, format_record},
};

use crate::filter::{self, RecordFilter};

/// BCF type codes (HTSlib's `BCF_BT_*`).
const BT_INT8: u8 = 1;
const BT_INT16: u8 = 2;
const BT_INT32: u8 = 3;
const BT_FLOAT: u8 = 5;
const BT_CHAR: u8 = 7;

/// What a field prints (bcftools' `T_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Kind {
    Chrom,
    Pos,
    Id,
    Ref,
    Alt,
    Qual,
    Filter,
    Info,
    Format,
    Sample,
    Sep,
    IsTs,
    Type,
    Mask,
    Gt,
    Tgt,
    Line,
    ChromPosId,
    GtToProb3,
    PlToProb3,
    GpToProb3,
    FirstAlt,
    IupacGt,
    GtToHap,
    GtToHap2,
    Tbcsq,
    End,
    Pos0,
    End0,
    Rsx,
    Vkx,
    Pbinom,
    Npass,
    FilterExpr,
}

/// Where a field sits in the format string.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Scope {
    /// Outside square brackets: printed once per record.
    Site,
    /// Inside square brackets: printed once per sample (bcftools' `is_gtf`).
    Sample,
}

/// Per-field state that only some kinds have (bcftools' `fmt_t.usr`).
enum Extra {
    None,
    /// `%N_PASS(EXPR)`.
    Npass(Box<dyn RecordFilter + Send>),
    /// `%MAX(...)` and the other filter functions: the filter and, inside brackets, the values
    /// of the current record.
    FilterExpr {
        filter: Box<dyn RecordFilter + Send>,
        values: Vec<f64>,
    },
    /// `%TBCSQ`: the INFO/BCSQ consequences of the current record.
    Tbcsq(Tbcsq),
    /// `%PBINOM(TAG)`: the GT field of the current record, if usable.
    Pbinom {
        gt: Option<usize>,
    },
}

#[derive(Default)]
struct Tbcsq {
    /// The consequences, as ranges into `text`.
    items: Vec<(usize, usize)>,
    text: Vec<u8>,
    hap1: Vec<u8>,
    hap2: Vec<u8>,
}

/// One parsed field (bcftools' `fmt_t`).
struct Field {
    kind: Kind,
    key: Option<Vec<u8>>,
    /// The header ID of the tag, or -1.
    id: i32,
    /// Whether the field is inside square brackets ([`Scope::Sample`]).
    is_gt_field: bool,
    /// Whether the per-record lookup (`init_format`) was done.
    ready: bool,
    subscript: i32,
    /// The record's FORMAT field for this tag (an index into `Record::format_fields`).
    fmt: Option<usize>,
    extra: Extra,
}

/// A compiled format string (bcftools' `convert_t`).
pub struct Convert {
    fields: Vec<Field>,
    state: State,
}

/// Everything but the fields, so a field can be processed while the rest is borrowed.
struct State {
    /// Header sample indices, in output order.
    samples: Vec<usize>,
    max_unpack: Unpack,
    /// Whether a field reads ID, REF or ALT (or might: anything not known not to).
    reads_alleles: bool,
    format_str: String,
    undef_info_tag: Option<Vec<u8>>,
    print_filtered: Option<Vec<u8>>,
    allow_undef_tags: bool,
    header_samples: bool,
    no_hdr_indices: bool,
    /// `%IUPACGT`'s allele-to-base map of the current record.
    iupac: Vec<i8>,
    ibuf: Vec<i32>,
    fbuf: Vec<f32>,
}

/// The tag names that are VCF columns (bcftools' `_SET_NON_FORMAT_TAGS`).
fn non_format_tag(key: &[u8]) -> Option<Kind> {
    Some(match key {
        b"CHROM" => Kind::Chrom,
        b"POS" => Kind::Pos,
        b"POS0" => Kind::Pos0,
        b"END" => Kind::End,
        b"END0" => Kind::End0,
        b"ID" => Kind::Id,
        b"REF" => Kind::Ref,
        b"FIRST_ALT" => Kind::FirstAlt,
        b"QUAL" => Kind::Qual,
        b"TYPE" => Kind::Type,
        b"FILTER" => Kind::Filter,
        b"IS_TS" => Kind::IsTs,
        b"MASK" => Kind::Mask,
        b"LINE" => Kind::Line,
        _ => return None,
    })
}

/// The filter functions that can be printed (bcftools' `_SET_FILTER_EXPR`), matched
/// case-insensitively as prefixes.
const FILTER_FUNCS: [&[u8]; 28] = [
    b"MAX(",
    b"MIN(",
    b"MEAN(",
    b"MEDIAN(",
    b"AVG(",
    b"SUM(",
    b"ABS(",
    b"COUNT(",
    b"STDEV(",
    b"STRLEN(",
    b"BINOM(",
    b"PHRED(",
    b"SMPL_MAX(",
    b"SMPL_MIN(",
    b"SMPL_MEAN(",
    b"SMPL_MEDIAN(",
    b"SMPL_AVG(",
    b"SMPL_STDEV(",
    b"SMPL_SUM(",
    b"SMPL_COUNT(",
    b"sMAX(",
    b"sMIN(",
    b"sMEAN(",
    b"sMEDIAN(",
    b"sAVG(",
    b"sSTDEV(",
    b"sSUM(",
    b"sCOUNT(",
];

fn is_filter_func(s: &[u8]) -> bool {
    FILTER_FUNCS
        .iter()
        .any(|f| s.len() >= f.len() && s[..f.len()].eq_ignore_ascii_case(f))
}

fn is_tag_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// Whether `id` is defined for `kind` (HTSlib's `bcf_hdr_idinfo_exists`, false for -1).
fn id_exists(header: &Header, kind: HeaderRecordType, id: i32) -> bool {
    u32::try_from(id).is_ok_and(|id| header.id_exists(kind, id))
}

fn header_id(header: &Header, key: &[u8]) -> i32 {
    header.id(key).map_or(-1, |i| i as i32)
}

/// Parses `{N}` at `*q` (bcftools' `parse_subscript`); -1 when absent.
fn parse_subscript(s: &[u8], q: &mut usize) -> i32 {
    if s.get(*q) != Some(&b'{') {
        return -1;
    }
    let mut r = *q + 1;
    while r < s.len() && s[r] != b'}' && s[r].is_ascii_digit() {
        r += 1;
    }
    if s.get(r) != Some(&b'}') {
        return -1;
    }
    let idx = std::str::from_utf8(&s[*q + 1..r])
        .ok()
        .and_then(|d| d.parse::<i64>().ok())
        .unwrap_or(0) as i32;
    *q = r + 1;
    idx
}

impl Convert {
    /// Parses a format string (bcftools' `convert_init`). `samples` lists the header samples
    /// to print, in order; `None` prints all of them.
    ///
    /// # Errors
    ///
    /// Fails on a malformed format string, or when a filter expression in it does not compile.
    pub fn new(
        header: &Header,
        samples: Option<&[usize]>,
        format_str: &str,
    ) -> anyhow::Result<Self> {
        let mut conv = Self {
            fields: Vec::new(),
            state: State {
                samples: Vec::new(),
                max_unpack: Unpack::STR,
                reads_alleles: false,
                format_str: format_str.to_owned(),
                undef_info_tag: None,
                print_filtered: None,
                allow_undef_tags: false,
                header_samples: false,
                no_hdr_indices: false,
                iupac: Vec::new(),
                ibuf: Vec::new(),
                fbuf: Vec::new(),
            },
        };
        let s = format_str.as_bytes();
        let mut scope = Scope::Site;
        let mut p = 0;
        while p < s.len() {
            match s[p] {
                b'[' => {
                    scope = Scope::Sample;
                    p += 1;
                }
                b']' => {
                    scope = Scope::Site;
                    conv.register(header, None, Scope::Site, Kind::Sep)?;
                    p += 1;
                }
                b'%' => p = conv.parse_tag(header, s, p, scope)?,
                _ => p = conv.parse_sep(header, s, p, scope)?,
            }
        }
        if scope == Scope::Sample {
            bail!(
                "Could not parse the format string, missing the square bracket \"]\": {format_str}"
            );
        }
        conv.state.samples = match samples {
            Some(list) if !list.is_empty() => list.to_vec(),
            _ => (0..header.sample_count()).collect(),
        };
        Ok(conv)
    }

    fn parse_error(&self) -> anyhow::Error {
        // The trailing newline keeps a newline that ends the format string (`-a`) printed, as
        // bcftools' `error("...%s\n")` does.
        anyhow::anyhow!("Could not parse format string: {}\n", self.state.format_str)
    }

    /// Parses one `%TAG` at `p` (bcftools' `parse_tag`); returns the position after it.
    fn parse_tag(
        &mut self,
        header: &Header,
        s: &[u8],
        mut p: usize,
        scope: Scope,
    ) -> anyhow::Result<usize> {
        let is_vcf_column = s.get(p + 1) == Some(&b'/');
        if is_vcf_column {
            p += 1;
        }
        p += 1;
        let mut q = p;
        while q < s.len() && is_tag_char(s[q]) {
            q += 1;
        }
        if q == p {
            return Err(self.parse_error());
        }
        let name = &s[p..q];
        match scope {
            Scope::Sample if is_vcf_column => {
                if let Some(kind) = non_format_tag(name) {
                    self.register(header, Some(name), scope, kind)?;
                } else if name == b"ALT" {
                    let i = self.register(header, Some(name), scope, Kind::Alt)?;
                    self.fields[i].subscript = parse_subscript(s, &mut q);
                } else {
                    bail!(
                        "Could not parse tag: {} .. {}",
                        String::from_utf8_lossy(name),
                        self.state.format_str
                    );
                }
                Ok(q)
            }
            Scope::Sample => self.parse_sample_tag(header, s, p, q),
            Scope::Site => self.parse_site_tag(header, s, p, q),
        }
    }

    /// Parses the `%TAG` named `s[p..q]` inside square brackets (part of bcftools'
    /// `parse_tag`); returns the position after it.
    fn parse_sample_tag(
        &mut self,
        header: &Header,
        s: &[u8],
        p: usize,
        mut q: usize,
    ) -> anyhow::Result<usize> {
        const SCOPE: Scope = Scope::Sample;
        if is_filter_func(&s[p..]) {
            return self.set_filter_expr(header, s, p, SCOPE);
        }
        let name = &s[p..q];
        match name {
            b"SAMPLE" => {
                self.register(header, Some(b"SAMPLE"), SCOPE, Kind::Sample)?;
            }
            b"GT" => {
                self.register(header, Some(b"GT"), SCOPE, Kind::Gt)?;
            }
            b"TGT" => {
                self.register(header, Some(b"GT"), SCOPE, Kind::Tgt)?;
            }
            b"TBCSQ" => {
                let i = self.register(header, Some(b"BCSQ"), SCOPE, Kind::Tbcsq)?;
                let mut sub = parse_subscript(s, &mut q);
                if sub == -1 {
                    if s[q..].starts_with(b"{*}") {
                        sub = 0;
                        q += 3;
                    }
                } else {
                    sub += 1;
                }
                self.fields[i].subscript = sub;
            }
            b"IUPACGT" => {
                self.register(header, Some(b"GT"), SCOPE, Kind::IupacGt)?;
            }
            b"INFO" => {
                if s.get(q) != Some(&b'/') {
                    let id = header_id(header, name);
                    if id_exists(header, HeaderRecordType::Info, id) {
                        bail!(
                            "Could not parse format string \"{0}\". Did you mean %INFO/{1}?",
                            self.state.format_str,
                            String::from_utf8_lossy(name)
                        );
                    }
                    return Err(self.parse_error());
                }
                q = self.parse_info_slash(header, s, q, SCOPE)?;
            }
            b"PBINOM" => {
                if s.get(q) != Some(&b'(') {
                    bail!("Could not parse the expression: {}", self.state.format_str);
                }
                q += 1;
                let start = q;
                while q < s.len() && s[q] != b')' {
                    q += 1;
                }
                if q == start {
                    return Err(self.parse_error());
                }
                let tag = s[start..q].to_vec();
                self.register(header, Some(&tag), SCOPE, Kind::Pbinom)?;
                q += 1;
            }
            b"N_PASS" => bail!("N_PASS() must be placed outside the square brackets"),
            _ => {
                let i = self.register(header, Some(name), SCOPE, Kind::Format)?;
                self.fields[i].subscript = parse_subscript(s, &mut q);
            }
        }
        Ok(q)
    }

    /// Parses the `%TAG` named `s[p..q]` outside square brackets (part of bcftools'
    /// `parse_tag`); returns the position after it.
    fn parse_site_tag(
        &mut self,
        header: &Header,
        s: &[u8],
        p: usize,
        mut q: usize,
    ) -> anyhow::Result<usize> {
        const SCOPE: Scope = Scope::Site;
        let name = &s[p..q];
        if let Some(kind) = non_format_tag(name) {
            self.register(header, Some(name), SCOPE, kind)?;
            return Ok(q);
        }
        if is_filter_func(&s[p..]) {
            return self.set_filter_expr(header, s, p, SCOPE);
        }
        match name {
            b"ALT" => {
                let i = self.register(header, Some(name), SCOPE, Kind::Alt)?;
                self.fields[i].subscript = parse_subscript(s, &mut q);
            }
            b"_CHROM_POS_ID" => {
                self.register(header, Some(name), SCOPE, Kind::ChromPosId)?;
            }
            b"_GT_TO_PROB3" => {
                self.register(header, Some(name), SCOPE, Kind::GtToProb3)?;
            }
            b"_PL_TO_PROB3" => {
                self.register(header, Some(name), SCOPE, Kind::PlToProb3)?;
            }
            b"_GP_TO_PROB3" => {
                self.register(header, Some(name), SCOPE, Kind::GpToProb3)?;
            }
            b"_GT_TO_HAP" => {
                self.register(header, Some(name), SCOPE, Kind::GtToHap)?;
            }
            b"_GT_TO_HAP2" => {
                self.register(header, Some(name), SCOPE, Kind::GtToHap2)?;
            }
            b"RSX" => {
                self.register(header, Some(name), SCOPE, Kind::Rsx)?;
            }
            b"VKX" => {
                self.register(header, Some(name), SCOPE, Kind::Vkx)?;
            }
            b"PBINOM" => {
                bail!("Error: PBINOM() is currently supported only with FORMAT tags. (todo)")
            }
            b"INFO" => {
                if s.get(q) == Some(&b'/') {
                    q = self.parse_info_slash(header, s, q, SCOPE)?;
                } else {
                    self.register(header, None, SCOPE, Kind::Info)?;
                }
            }
            b"FORMAT" => {
                self.register(header, None, SCOPE, Kind::Format)?;
            }
            b"N_PASS" => {
                if s.get(q) != Some(&b'(') {
                    bail!("Could not parse the expression: {}", self.state.format_str);
                }
                q += 1;
                let start = q;
                let mut nopen = 1;
                while q < s.len() && nopen > 0 {
                    match s[q] {
                        b'(' => nopen += 1,
                        b')' => nopen -= 1,
                        _ => {}
                    }
                    q += 1;
                }
                if q == start || nopen > 0 {
                    return Err(self.parse_error());
                }
                let expr = s[start..q - 1].to_vec();
                self.register(header, Some(&expr), SCOPE, Kind::Npass)?;
            }
            _ => {
                let i = self.register(header, Some(name), SCOPE, Kind::Info)?;
                self.fields[i].subscript = parse_subscript(s, &mut q);
            }
        }
        Ok(q)
    }

    /// Parses the `TAG{N}` after `%INFO/`, with `q` at the slash.
    fn parse_info_slash(
        &mut self,
        header: &Header,
        s: &[u8],
        q: usize,
        scope: Scope,
    ) -> anyhow::Result<usize> {
        let p = q + 1;
        let mut q = p;
        while q < s.len() && is_tag_char(s[q]) {
            q += 1;
        }
        if q == p {
            return Err(self.parse_error());
        }
        let i = self.register(header, Some(&s[p..q]), scope, Kind::Info)?;
        self.fields[i].subscript = parse_subscript(s, &mut q);
        Ok(q)
    }

    /// Registers a printed filter function, up to the first `)` (bcftools' `set_filter_expr`).
    fn set_filter_expr(
        &mut self,
        header: &Header,
        s: &[u8],
        p: usize,
        scope: Scope,
    ) -> anyhow::Result<usize> {
        let Some(end) = s[p..].iter().position(|&c| c == b')') else {
            return Err(self.parse_error());
        };
        let key = s[p..=p + end].to_vec();
        self.register(header, Some(&key), scope, Kind::FilterExpr)?;
        Ok(p + key.len())
    }

    /// Parses literal text up to the next `[`, `]`, or `%`, with `\n`, `\t`, and `\X` escapes
    /// (bcftools' `parse_sep`).
    fn parse_sep(
        &mut self,
        header: &Header,
        s: &[u8],
        p: usize,
        scope: Scope,
    ) -> anyhow::Result<usize> {
        let mut q = p;
        let mut text = Vec::new();
        while q < s.len() && !matches!(s[q], b'[' | b']' | b'%') {
            if s[q] == b'\\' {
                q += 1;
                match s.get(q) {
                    Some(b'n') => text.push(b'\n'),
                    Some(b't') => text.push(b'\t'),
                    Some(&c) => text.push(c),
                    // bcftools copies the string's NUL: the field ends there, and is not empty.
                    None => {
                        text.push(0);
                        break;
                    }
                }
            } else {
                text.push(s[q]);
            }
            q += 1;
        }
        if text.is_empty() {
            return Err(self.parse_error());
        }
        self.register(header, Some(&text), scope, Kind::Sep)?;
        Ok(q)
    }

    /// Adds a field (bcftools' `register_tag`); returns its index.
    fn register(
        &mut self,
        header: &Header,
        key: Option<&[u8]>,
        scope: Scope,
        kind: Kind,
    ) -> anyhow::Result<usize> {
        // Keys are C strings: they end at a NUL.
        let key = key.map(|k| &k[..k.iter().position(|&c| c == 0).unwrap_or(k.len())]);
        let mut field = Field {
            kind,
            key: key.map(<[u8]>::to_vec),
            id: -1,
            is_gt_field: scope == Scope::Sample,
            ready: false,
            subscript: -1,
            fmt: None,
            extra: Extra::None,
        };
        // Non-FORMAT tags (CHROM, INFO tags, ...) may appear among the FORMAT tags.
        if let Some(key) = key {
            let id = header_id(header, key);
            if field.kind == Kind::Format && !id_exists(header, HeaderRecordType::Format, id) {
                if let Some(k) = non_format_tag(key) {
                    field.kind = k;
                } else {
                    match key {
                        b"ALT" => field.kind = Kind::Alt,
                        b"_CHROM_POS_ID" => field.kind = Kind::ChromPosId,
                        b"RSX" => field.kind = Kind::Rsx,
                        b"VKX" => field.kind = Kind::Vkx,
                        _ => {
                            if id_exists(header, HeaderRecordType::Info, id) {
                                field.kind = Kind::Info;
                            }
                        }
                    }
                }
            } else if field.kind == Kind::Pbinom {
                field.id = id;
                if !id_exists(header, HeaderRecordType::Format, id) {
                    bail!(
                        "No such FORMAT tag defined in the header: {}",
                        String::from_utf8_lossy(key)
                    );
                }
                field.extra = Extra::Pbinom { gt: None };
            } else if field.kind == Kind::Npass {
                let flt = filter::init(header, &String::from_utf8_lossy(key))?;
                self.state.max_unpack |= flt.max_unpack();
                field.extra = Extra::Npass(flt);
            } else if field.kind == Kind::FilterExpr {
                let flt = filter::init(header, &String::from_utf8_lossy(key))?;
                self.state.max_unpack |= flt.max_unpack();
                field.extra = Extra::FilterExpr {
                    filter: flt,
                    values: Vec::new(),
                };
            }
        }
        let unpack = match field.kind {
            Kind::End | Kind::End0 | Kind::Info => Unpack::INFO,
            Kind::Filter => Unpack::FLT,
            Kind::Format
            | Kind::Gt
            | Kind::Tgt
            | Kind::IupacGt
            | Kind::GtToHap
            | Kind::GtToHap2
            | Kind::Tbcsq
            | Kind::Line
            | Kind::Pbinom => Unpack::FMT,
            _ => Unpack::empty(),
        };
        self.state.max_unpack |= unpack;
        self.state.reads_alleles |= !matches!(
            field.kind,
            Kind::Chrom
                | Kind::Pos
                | Kind::Pos0
                | Kind::Qual
                | Kind::Sep
                | Kind::Mask
                | Kind::Sample
                | Kind::Format
                | Kind::Gt
                | Kind::Line
        );
        if field.kind == Kind::Tbcsq {
            field.extra = Extra::Tbcsq(Tbcsq::default());
        }
        if let Some(key) = key
            && field.kind == Kind::Info
        {
            field.id = header_id(header, key);
            if !id_exists(header, HeaderRecordType::Info, field.id) {
                field.id = -1;
                self.state.undef_info_tag = Some(key.to_vec());
            }
        }
        self.fields.push(field);
        Ok(self.fields.len() - 1)
    }

    /// Prints `.` for undefined tags instead of failing (`-u`).
    pub fn set_allow_undef_tags(&mut self) {
        self.state.allow_undef_tags = true;
    }

    /// Prefixes per-sample header columns with the sample name.
    pub fn set_header_samples(&mut self) {
        self.state.header_samples = true;
    }

    /// Prints `text` for each field of samples that fail the filter (`query -F`).
    pub fn set_print_filtered(&mut self, text: &[u8]) {
        self.state.print_filtered = Some(text.to_vec());
    }

    /// Omits the `[N]` column numbers from the header (`-HH`).
    pub fn set_no_hdr_indices(&mut self) {
        self.state.no_hdr_indices = true;
    }

    /// Adds a newline at the end unless the format has one (bcftools' `force_newline`).
    pub fn force_newline(&mut self, header: &Header) -> anyhow::Result<()> {
        let has_newline = self
            .fields
            .iter()
            .any(|f| f.key.as_ref().is_some_and(|k| k.contains(&b'\n')));
        if !has_newline {
            self.register(header, Some(b"\n"), Scope::Site, Kind::Sep)?;
        }
        Ok(())
    }

    /// Returns how much of each record the fields need (bcftools' `convert_max_unpack`).
    #[must_use]
    pub fn max_unpack(&self) -> Unpack {
        self.state.max_unpack
    }

    /// What formatting a line decodes: [`Self::max_unpack`] without ID, REF and ALT when no
    /// field reads them (bcftools always decodes them).
    #[must_use]
    pub fn line_unpack(&self) -> Unpack {
        if self.state.reads_alleles {
            self.state.max_unpack
        } else {
            self.state.max_unpack.difference(Unpack::STR)
        }
    }

    /// Whether the format reads an INFO tag (an `INFO/...` entry of bcftools'
    /// `convert_list_used_tags`).
    #[must_use]
    pub fn uses_info_tags(&self) -> bool {
        self.fields
            .iter()
            .any(|f| f.kind == Kind::Info && f.key.is_some())
    }

    /// Writes the header line (bcftools' `convert_header`): `#` and the column names, numbered
    /// `[1]`, `[2]`, ...; nothing when the format has `%LINE`.
    pub fn header_line(&mut self, header: &Header, out: &mut Vec<u8>) {
        if self.fields.iter().any(|f| f.kind == Kind::Line) {
            return;
        }
        let l_ori = out.len();
        loop {
            let mut icol = 0;
            let mut has_fmt_newline = false;
            let mut restart = false;
            out.push(b'#');
            let mut i = 0;
            while i < self.fields.len() {
                if self.fields[i].is_gt_field {
                    let mut j = i;
                    while j < self.fields.len() && self.fields[j].is_gt_field {
                        j += 1;
                    }
                    for &ks in &self.state.samples {
                        for f in &self.fields[i..j] {
                            if f.kind == Kind::Sep {
                                if let Some(k) = &f.key {
                                    for &c in k {
                                        if c == b'\n' {
                                            has_fmt_newline = true;
                                        } else {
                                            out.push(c);
                                        }
                                    }
                                }
                            } else {
                                icol += 1;
                                if !self.state.no_hdr_indices {
                                    out.extend_from_slice(format!("[{icol}]").as_bytes());
                                }
                                if self.state.header_samples {
                                    out.extend_from_slice(&header.samples()[ks]);
                                    out.push(b':');
                                }
                                push_key(out, f.key.as_deref());
                            }
                        }
                        if has_fmt_newline {
                            if !self.state.header_samples {
                                break;
                            }
                            // The per-sample output is split into lines, so sample names
                            // make no sense in the header.
                            self.state.header_samples = false;
                            out.truncate(l_ori);
                            restart = true;
                            break;
                        }
                    }
                    if restart {
                        break;
                    }
                    i = j;
                    continue;
                }
                let f = &self.fields[i];
                i += 1;
                if f.kind == Kind::Sep {
                    if let Some(k) = &f.key {
                        out.extend_from_slice(k);
                    }
                    continue;
                }
                icol += 1;
                if !self.state.no_hdr_indices {
                    out.extend_from_slice(format!("[{icol}]").as_bytes());
                }
                push_key(out, f.key.as_deref());
            }
            if restart {
                continue;
            }
            if has_fmt_newline {
                out.push(b'\n');
            }
            return;
        }
    }

    /// Formats a record into `out`, replacing its contents (bcftools' `convert_line`).
    /// `smpl_pass` marks the samples that passed a filter; the others are skipped, or printed
    /// with the `-F` text.
    ///
    /// # Errors
    ///
    /// Fails on an undefined tag (unless allowed) and on malformed values that bcftools treats
    /// as fatal.
    pub fn line(
        &mut self,
        header: &Header,
        rec: &mut Record,
        smpl_pass: Option<&[bool]>,
        out: &mut Vec<u8>,
    ) -> anyhow::Result<()> {
        if !self.state.allow_undef_tags
            && let Some(tag) = &self.state.undef_info_tag
        {
            let tag = String::from_utf8_lossy(tag);
            let mut msg = format!("Error: no such tag defined in the VCF header: INFO/{tag}");
            let id = header_id(header, tag.as_bytes());
            if id_exists(header, HeaderRecordType::Format, id) {
                msg.push_str(&format!(
                    ". FORMAT fields must be enclosed in square brackets, e.g. \"[ %{tag}]\""
                ));
            }
            bail!("{msg}");
        }
        rec.unpack(self.line_unpack());
        out.clear();
        let n = self.fields.len();
        let mut i = 0;
        while i < n {
            if !self.fields[i].is_gt_field {
                let f = &mut self.fields[i];
                match f.kind {
                    // Separators are most fields; print them without the dispatch.
                    Kind::Sep => match f.key.as_deref() {
                        Some(&[c]) => out.push(c),
                        Some(key) => out.extend_from_slice(key),
                        None => {}
                    },
                    Kind::Mask => {}
                    _ => process(f, &mut self.state, header, rec, None, out)?,
                }
                i += 1;
                continue;
            }
            let mut j = i;
            while j < n && self.fields[j].is_gt_field {
                self.fields[j].ready = false;
                j += 1;
            }
            if j - i <= FAST_MAX
                && self.fields[i..j].iter().all(|f| {
                    matches!(f.kind, Kind::Sep | Kind::Gt | Kind::Sample)
                        || (f.kind == Kind::Format && f.key.is_some())
                })
            {
                self.simple_block(i..j, header, rec, smpl_pass, out)?;
                i = j;
                continue;
            }
            for js in 0..self.state.samples.len() {
                let ks = self.state.samples[js];
                if let Some(pass) = smpl_pass
                    && !pass.get(ks).copied().unwrap_or(false)
                {
                    print_filtered(&self.fields[i..j], &self.state, out);
                    continue;
                }
                // A field that prints nothing (only %TBCSQ does) drops the sample's fields.
                let l_start = out.len();
                for k in i..j {
                    let f = &mut self.fields[k];
                    if f.kind == Kind::Mask {
                        continue;
                    }
                    let l = out.len();
                    process(f, &mut self.state, header, rec, Some(ks), out)?;
                    if l == out.len() {
                        out.truncate(l_start);
                        break;
                    }
                }
            }
            i = j;
        }
        Ok(())
    }
}

/// The most fields in a per-sample block printed by [`Convert::simple_block`].
const FAST_MAX: usize = 16;

/// What a field of a simple per-sample block prints.
#[derive(Clone, Copy)]
enum Simple {
    Sep,
    Sample,
    Missing,
    /// A FORMAT field (GT or another tag): its type, values per sample, and bytes per sample.
    Value {
        is_gt: bool,
        ty: u8,
        n: usize,
        size: usize,
    },
}

impl Convert {
    /// Prints a per-sample block of separators, `%SAMPLE`, `%GT`, and FORMAT tags, the common
    /// case, resolving the fields once per record rather than once per sample.
    fn simple_block(
        &mut self,
        range: std::ops::Range<usize>,
        header: &Header,
        rec: &Record,
        smpl_pass: Option<&[bool]>,
        out: &mut Vec<u8>,
    ) -> anyhow::Result<()> {
        let mut kinds = [Simple::Sep; FAST_MAX];
        let mut data: [&[u8]; FAST_MAX] = [&[]; FAST_MAX];
        let mut prepared = false;
        let fields = range.len();
        for &ks in &self.state.samples {
            if let Some(pass) = smpl_pass
                && !pass.get(ks).copied().unwrap_or(false)
            {
                print_filtered(&self.fields[range.clone()], &self.state, out);
                continue;
            }
            if !prepared {
                // bcftools looks the tags up at the first printed sample.
                for (k, f) in self.fields[range.clone()].iter_mut().enumerate() {
                    kinds[k] = match f.kind {
                        Kind::Sep => Simple::Sep,
                        Kind::Sample => Simple::Sample,
                        _ => {
                            init_format(f, &self.state, header, rec)?;
                            match f.fmt {
                                None => Simple::Missing,
                                Some(fi) => {
                                    let field = &rec.format_fields()[fi];
                                    data[k] = rec.format_data(field);
                                    Simple::Value {
                                        is_gt: f.kind == Kind::Gt,
                                        ty: field.type_code(),
                                        n: field.n(),
                                        size: field.size(),
                                    }
                                }
                            }
                        }
                    };
                }
                prepared = true;
            }
            let l_start = out.len();
            for k in 0..fields {
                let l = out.len();
                let f = &self.fields[range.start + k];
                match kinds[k] {
                    Simple::Sep => match f.key.as_deref() {
                        Some(&[c]) => out.push(c),
                        Some(key) => out.extend_from_slice(key),
                        None => {}
                    },
                    Simple::Sample => {
                        if let Some(name) = header.samples().get(ks) {
                            out.extend_from_slice(name);
                        }
                    }
                    Simple::Missing => out.push(b'.'),
                    Simple::Value { is_gt, ty, n, size } => {
                        let sample = &data[k][ks * size..(ks + 1) * size];
                        if is_gt {
                            let _ = format_genotype(out, 0, ty, n, sample);
                        } else {
                            format_value(ty, n, sample, f.subscript, out);
                        }
                    }
                }
                // A field that prints nothing drops the sample's fields.
                if l == out.len() {
                    out.truncate(l_start);
                    break;
                }
            }
        }
        Ok(())
    }
}

/// Prints the `-F` text for each field of a sample that failed the filter.
fn print_filtered(fields: &[Field], st: &State, out: &mut Vec<u8>) {
    let Some(text) = &st.print_filtered else {
        return;
    };
    for f in fields {
        if f.kind == Kind::Sep {
            push_key(out, f.key.as_deref());
        } else {
            out.extend_from_slice(text);
        }
    }
}

/// `kputs(key)`, printing `(null)` for no key as glibc's `%s` does.
fn push_key(out: &mut Vec<u8>, key: Option<&[u8]>) {
    match key {
        Some(k) => out.extend_from_slice(k),
        None => out.extend_from_slice(b"(null)"),
    }
}

fn location(header: &Header, rec: &Record) -> String {
    format!(
        "{}:{}",
        String::from_utf8_lossy(rec.chrom_safe(header)),
        rec.pos() + 1
    )
}

/// Appends a value as C's `%f`.
pub fn push_fixed6(out: &mut Vec<u8>, v: f64) {
    if v.is_nan() {
        out.extend_from_slice(if v.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
    } else if v.is_infinite() {
        out.extend_from_slice(if v < 0.0 { b"-inf" } else { b"inf" });
    } else {
        // Writing to a `Vec` cannot fail.
        let _ = write!(out, "{v:.6}");
    }
}

/// Appends field `idx` of a comma-separated string (bcftools' `_copy_field`), `.` when absent
/// or empty.
fn copy_field(src: &[u8], idx: i32, out: &mut Vec<u8>) {
    let len = src.len();
    let mut n = 0;
    let mut ibeg = 0;
    while ibeg < len && src[ibeg] != 0 && n < idx {
        if src[ibeg] == b',' {
            n += 1;
        }
        ibeg += 1;
    }
    if ibeg == len {
        out.push(b'.');
        return;
    }
    let mut iend = ibeg;
    while iend < len && src[iend] != 0 && src[iend] != b',' {
        iend += 1;
    }
    if iend > ibeg {
        out.extend_from_slice(&src[ibeg..iend]);
    } else {
        out.push(b'.');
    }
}

/// Reads value `idx` of an integer array as `i32`, mapping the narrow sentinels (bcftools'
/// `bcf_array_ivalue`).
fn array_ivalue(data: &[u8], ty: u8, idx: usize) -> i32 {
    match ty {
        BT_INT8 => {
            let v = data[idx] as i8;
            if v == INT8_MISSING {
                INT32_MISSING
            } else if v == INT8_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(v)
            }
        }
        BT_INT16 => {
            let v = le_i16(&data[idx * 2..]);
            if v == INT16_MISSING {
                INT32_MISSING
            } else if v == INT16_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(v)
            }
        }
        _ => le_i32(&data[idx * 4..]),
    }
}

/// Prints one field (the bcftools `process_*` handlers). `isample` is the header sample for
/// per-sample fields.
fn process(
    f: &mut Field,
    st: &mut State,
    header: &Header,
    rec: &mut Record,
    isample: Option<usize>,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    match f.kind {
        Kind::Chrom => out.extend_from_slice(rec.chrom_safe(header)),
        Kind::Pos => push_int(out, rec.pos() + 1),
        Kind::Pos0 => push_int(out, rec.pos()),
        Kind::End => push_int(out, rec.pos() + rec.rlen()),
        Kind::End0 => push_int(out, rec.pos() + rec.rlen() - 1),
        Kind::Id => out.extend_from_slice(rec.id()),
        Kind::Ref => out.extend_from_slice(allele0(rec)),
        Kind::Alt => process_alt(f, rec, out),
        Kind::FirstAlt => {
            if rec.allele_count() <= 1 {
                out.push(b'.');
            } else {
                out.extend_from_slice(rec.allele(1));
            }
        }
        Kind::Qual => match rec.qual() {
            None => out.push(b'.'),
            Some(q) => push_double(out, f64::from(q)),
        },
        Kind::Filter => process_filter(header, rec, out),
        Kind::Info => process_info(f, header, rec, out)?,
        Kind::Format => {
            if f.key.is_some() {
                process_format(f, st, header, rec, isample.unwrap_or(0), out)?;
            } else {
                process_complete_format(st, header, rec, out);
            }
        }
        Kind::Sample => {
            if let Some(name) = isample.and_then(|i| header.samples().get(i)) {
                out.extend_from_slice(name);
            }
        }
        Kind::Sep => match f.key.as_deref() {
            Some(&[c]) => out.push(c),
            Some(k) => out.extend_from_slice(k),
            None => {}
        },
        Kind::IsTs => process_is_ts(rec, out),
        Kind::Type => process_type(rec, out),
        Kind::Mask => {}
        Kind::Gt => process_gt(f, st, header, rec, isample.unwrap_or(0), out)?,
        Kind::Tgt => process_tgt(f, st, header, rec, isample.unwrap_or(0), out)?,
        Kind::IupacGt => process_iupac_gt(f, st, header, rec, isample.unwrap_or(0), out)?,
        Kind::Line => {
            let start = out.len();
            let _ = format_record(header, rec, out);
            if out.len() > start && out.last() == Some(&b'\n') {
                out.pop();
            }
        }
        Kind::ChromPosId => {
            let id = rec.id();
            if id != b"." {
                out.extend_from_slice(id);
            } else {
                out.extend_from_slice(rec.chrom_safe(header));
                out.push(b':');
                push_int(out, rec.pos() + 1);
            }
        }
        Kind::GtToProb3 => hap::gt_to_prob3(st, header, rec, out)?,
        Kind::PlToProb3 => hap::pl_to_prob3(st, header, rec, out)?,
        Kind::GpToProb3 => hap::gp_to_prob3(st, header, rec, out)?,
        Kind::GtToHap => hap::gt_to_hap(st, header, rec, out, hap::Haploid::Keep)?,
        Kind::GtToHap2 => hap::gt_to_hap(st, header, rec, out, hap::Haploid::ToDiploid)?,
        Kind::Tbcsq => process_tbcsq(f, st, header, rec, isample.unwrap_or(0), out)?,
        Kind::Rsx => {
            let id = rec.id();
            let digits = id.get(2..).unwrap_or(b"");
            let (v, _) = pasta::variant::number::strtol(digits, 0);
            let _ = write!(out, "{:08x}", v as u32);
        }
        Kind::Vkx => process_vkx(header, rec, out),
        Kind::Pbinom => process_pbinom(f, header, rec, isample.unwrap_or(0), out)?,
        Kind::Npass => process_npass(f, st, header, rec, out)?,
        Kind::FilterExpr => process_filter_expr(f, header, rec, isample.unwrap_or(0), out)?,
    }
    Ok(())
}

/// `%FILTER`: the filter names, separated by `;`, or `.`.
#[inline]
fn process_filter(header: &Header, rec: &Record, out: &mut Vec<u8>) {
    if rec.filters().is_empty() {
        out.push(b'.');
        return;
    }
    for (i, &id) in rec.filters().iter().enumerate() {
        if i > 0 {
            out.push(b';');
        }
        if let Ok(id) = u32::try_from(id)
            && let Some(name) = header.id_name(id)
        {
            out.extend_from_slice(name);
        }
    }
}

/// `%IS_TS`: 1 when the first ALT of a SNP or MNP is a transition, else 0.
#[inline]
fn process_is_ts(rec: &mut Record, out: &mut Vec<u8>) {
    let mut is_ts = false;
    if rec
        .variant_types()
        .intersects(VariantType::SNP | VariantType::MNP)
        && rec.allele_count() > 1
    {
        let a = acgt2int(allele0(rec).first().copied().unwrap_or(0));
        let b = acgt2int(rec.allele(1).first().copied().unwrap_or(0));
        is_ts = (a - b).abs() == 2;
    }
    out.push(if is_ts { b'1' } else { b'0' });
}

/// `%GT` of header sample `isample`.
#[inline]
fn process_gt(
    f: &mut Field,
    st: &State,
    header: &Header,
    rec: &Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    if !f.ready {
        init_format(f, st, header, rec)?;
    }
    if let Some(i) = f.fmt {
        let field = &rec.format_fields()[i];
        let data = rec.format_sample_data(field, isample);
        let _ = format_genotype(out, 0, field.type_code(), field.n(), data);
    } else {
        out.push(b'.');
    }
    Ok(())
}

/// `%VKX`: the VariantKey of CHROM, POS, REF, and the first ALT, in hexadecimal.
fn process_vkx(header: &Header, rec: &Record, out: &mut Vec<u8>) {
    let alt: &[u8] = if rec.allele_count() > 1 {
        rec.allele(1)
    } else {
        b""
    };
    let vk = variantkey::variantkey(rec.chrom_safe(header), rec.pos() as u32, allele0(rec), alt);
    let _ = write!(out, "{vk:016x}");
}

/// `%N_PASS(EXPR)`: how many of the printed samples pass the expression.
fn process_npass(
    f: &mut Field,
    st: &State,
    header: &Header,
    rec: &mut Record,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let Extra::Npass(flt) = &mut f.extra else {
        unreachable!("N_PASS fields carry a filter");
    };
    flt.pass(header, rec)?;
    let n = flt.sample_mask().map_or(0, |mask| {
        (0..st.samples.len())
            .filter(|&i| mask.get(i).copied().unwrap_or(false))
            .count()
    });
    push_double(out, n as f64);
    Ok(())
}

/// `%MAX(...)` and the other filter functions: the value of the expression, inside square
/// brackets the one of header sample `isample`, evaluated once per record.
#[inline]
fn process_filter_expr(
    f: &mut Field,
    header: &Header,
    rec: &mut Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let Extra::FilterExpr { filter, values } = &mut f.extra else {
        unreachable!("filter-function fields carry a filter");
    };
    let vals: &[f64] = if f.is_gt_field {
        if !f.ready {
            filter.pass(header, rec)?;
            values.clear();
            values.extend_from_slice(filter.doubles());
            f.ready = true;
        }
        values
    } else {
        filter.pass(header, rec)?;
        filter.doubles()
    };
    let mut i = isample;
    if i >= vals.len() {
        i = 0;
    }
    match vals.get(i) {
        Some(&v) => push_double(out, v),
        None => out.push(b'.'),
    }
    Ok(())
}

fn allele0(rec: &Record) -> &[u8] {
    if rec.allele_count() == 0 {
        b""
    } else {
        rec.allele(0)
    }
}

/// HTSlib's `bcf_acgt2int`.
fn acgt2int(c: u8) -> i32 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => -1,
    }
}

fn process_alt(f: &Field, rec: &Record, out: &mut Vec<u8>) {
    let n = rec.allele_count();
    if n <= 1 {
        out.push(b'.');
        return;
    }
    if f.subscript >= 0 {
        let i = f.subscript as usize + 1;
        if n > i {
            out.extend_from_slice(rec.allele(i));
        } else {
            out.push(b'.');
        }
        return;
    }
    for i in 1..n {
        if i > 1 {
            out.push(b',');
        }
        out.extend_from_slice(rec.allele(i));
    }
}

fn process_type(rec: &mut Record, out: &mut Vec<u8>) {
    let ty = rec.variant_types();
    let mut first = true;
    let mut put = |out: &mut Vec<u8>, s: &[u8]| {
        if !first {
            out.push(b',');
        }
        first = false;
        out.extend_from_slice(s);
    };
    if ty.is_empty() {
        put(out, b"REF");
    }
    for (flag, name) in [
        (VariantType::SNP, &b"SNP"[..]),
        (VariantType::MNP, b"MNP"),
        (VariantType::INDEL, b"INDEL"),
        (VariantType::OTHER, b"OTHER"),
        (VariantType::BND, b"BND"),
        (VariantType::OVERLAP, b"OVERLAP"),
    ] {
        if ty.contains(flag) {
            put(out, name);
        }
    }
}

/// Appends a single INFO value of type `ty` (`v1`), `.` when missing.
fn push_info_scalar(info_ty: u8, v1: i64, fval: f32, out: &mut Vec<u8>) -> anyhow::Result<()> {
    match info_ty {
        BT_INT8 if v1 == i64::from(INT8_MISSING) => out.push(b'.'),
        BT_INT16 if v1 == i64::from(INT16_MISSING) => out.push(b'.'),
        BT_INT32 if v1 == i64::from(INT32_MISSING) => out.push(b'.'),
        BT_INT8 | BT_INT16 | BT_INT32 => push_int(out, v1),
        BT_FLOAT => {
            if fval.to_bits() == FLOAT_MISSING {
                out.push(b'.');
            } else {
                push_double(out, f64::from(fval));
            }
        }
        BT_CHAR => out.push(v1 as u8),
        _ => bail!("Unexpected type {info_ty}"),
    }
    Ok(())
}

fn process_info(f: &Field, header: &Header, rec: &Record, out: &mut Vec<u8>) -> anyhow::Result<()> {
    if f.key.is_none() {
        // The whole INFO column.
        let mut first = true;
        for inf in rec.info_fields() {
            if inf.is_removed() {
                continue;
            }
            if !first {
                out.push(b';');
            }
            first = false;
            let Some(name) = u32::try_from(inf.key())
                .ok()
                .and_then(|k| header.id_name(k))
            else {
                continue;
            };
            out.extend_from_slice(name);
            if inf.is_empty() {
                continue;
            }
            out.push(b'=');
            if inf.len() == 1 {
                push_info_scalar(inf.type_code(), inf.int_value(), inf.float_value(), out)?;
            } else {
                format_array(out, inf.len(), inf.type_code(), rec.info_data(inf));
            }
        }
        if first {
            out.push(b'.');
        }
        return Ok(());
    }
    if f.id < 0 {
        out.push(b'.');
        return Ok(());
    }
    let Some(info) = rec.info_fields().iter().find(|inf| inf.key() == f.id) else {
        out.push(b'.');
        return Ok(());
    };
    if info.is_empty() {
        // A flag.
        out.push(b'1');
        return Ok(());
    }
    if info.len() == 1 {
        return push_info_scalar(info.type_code(), info.int_value(), info.float_value(), out);
    }
    let data = rec.info_data(info);
    if f.subscript >= 0 {
        let sub = f.subscript as usize;
        if info.len() <= sub {
            out.push(b'.');
            return Ok(());
        }
        match info.type_code() {
            BT_INT8 | BT_INT16 | BT_INT32 => {
                let v = array_ivalue(data, info.type_code(), sub);
                if v == INT32_MISSING || v == INT32_VECTOR_END {
                    out.push(b'.');
                } else {
                    push_int(out, i64::from(v));
                }
            }
            BT_FLOAT => {
                let v = le_u32(&data[sub * 4..]);
                if v == FLOAT_MISSING || v == FLOAT_VECTOR_END {
                    out.push(b'.');
                } else {
                    push_double(out, f64::from(f32::from_bits(v)));
                }
            }
            BT_CHAR => copy_field(data, f.subscript, out),
            ty => bail!("todo: type {ty}"),
        }
        return Ok(());
    }
    format_array(out, info.len(), info.type_code(), data);
    Ok(())
}

/// Looks up the FORMAT ID of a field's tag, keeping it once found: a definition never goes
/// away, but one may be added while reading VCF.
fn resolve_format_id(id: &mut i32, key: Option<&[u8]>, header: &Header) {
    if *id < 0 {
        let found = header_id(header, key.unwrap_or(b""));
        *id = if id_exists(header, HeaderRecordType::Format, found) {
            found
        } else {
            -1
        };
    }
}

/// Finds the record's FORMAT field of the tag (bcftools' `init_format`).
fn init_format(f: &mut Field, st: &State, header: &Header, rec: &Record) -> anyhow::Result<()> {
    resolve_format_id(&mut f.id, f.key.as_deref(), header);
    let key = f.key.as_deref().unwrap_or(b"");
    f.fmt = None;
    if f.id >= 0 {
        f.fmt = rec.format_fields().iter().position(|x| x.id() == f.id);
    } else if !st.allow_undef_tags {
        bail!(
            "Error: no such tag defined in the VCF header: FORMAT/{}",
            String::from_utf8_lossy(key)
        );
    }
    f.ready = true;
    Ok(())
}

fn process_format(
    f: &mut Field,
    st: &State,
    header: &Header,
    rec: &Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    if !f.ready {
        init_format(f, st, header, rec)?;
    }
    let Some(i) = f.fmt else {
        out.push(b'.');
        return Ok(());
    };
    let field = &rec.format_fields()[i];
    let data = rec.format_sample_data(field, isample);
    format_value(field.type_code(), field.n(), data, f.subscript, out);
    Ok(())
}

/// Prints one sample's values of a FORMAT field, or value `subscript` of them (the body of
/// bcftools' `process_format`).
fn format_value(ty: u8, n: usize, data: &[u8], subscript: i32, out: &mut Vec<u8>) {
    if subscript < 0 {
        format_array(out, n, ty, data);
        return;
    }
    let sub = subscript as usize;
    if n <= sub {
        out.push(b'.');
        return;
    }
    match ty {
        BT_FLOAT => {
            let v = le_u32(&data[sub * 4..]);
            if v == FLOAT_MISSING || v == FLOAT_VECTOR_END {
                out.push(b'.');
            } else {
                push_double(out, f64::from(f32::from_bits(v)));
            }
        }
        BT_CHAR => copy_field(data, subscript, out),
        ty => {
            let v = array_ivalue(data, ty, sub);
            if v == INT32_MISSING || v == INT32_VECTOR_END {
                out.push(b'.');
            } else {
                push_int(out, i64::from(v));
            }
        }
    }
}

/// `%FORMAT` outside brackets: the FORMAT column and all samples (bcftools'
/// `process_complete_format`).
fn process_complete_format(st: &State, header: &Header, rec: &Record, out: &mut Vec<u8>) {
    if st.samples.is_empty() {
        out.push(b'.');
        return;
    }
    if rec.format_count() == 0 {
        for _ in 0..=rec.sample_count() {
            out.extend_from_slice(b"\t.");
        }
        return;
    }
    let fmts = rec.format_fields();
    let mut gt_i = None;
    let mut first = true;
    for (i, f) in fmts.iter().enumerate() {
        if f.is_removed() || f.id() < 0 {
            continue;
        }
        if !first {
            out.push(b':');
        }
        first = false;
        let name = header.id_name(f.id() as u32).unwrap_or(b"");
        out.extend_from_slice(name);
        if name == b"GT" {
            gt_i = Some(i);
        }
    }
    if first {
        out.push(b'.');
    }
    for &s in &st.samples {
        out.push(b'\t');
        let mut first = true;
        for (i, f) in fmts.iter().enumerate() {
            if f.is_removed() {
                continue;
            }
            if !first {
                out.push(b':');
            }
            first = false;
            let data = rec.format_sample_data(f, s);
            if gt_i == Some(i) {
                let _ = format_genotype(out, 0, f.type_code(), f.n(), data);
            } else {
                format_array(out, f.n(), f.type_code(), data);
            }
        }
        if first {
            out.push(b'.');
        }
    }
}

/// Appends the alleles of one sample's 8-bit genotype (bcftools' `%TGT` loop).
fn push_gt_alleles(header: &Header, rec: &Record, x: &[u8], out: &mut Vec<u8>) {
    let n_allele = rec.allele_count();
    let mut l = 0;
    while l < x.len() && x[l] as i8 != INT8_VECTOR_END {
        let v = x[l] as i8;
        if l > 0 {
            out.push(if v & 1 != 0 { b'|' } else { b'/' });
        }
        if v >> 1 != 0 {
            let ial = (v >> 1) - 1;
            if ial < 0 || ial as usize >= n_allele {
                eprintln!(
                    "Warning: allele index {ial} >= n_allele {n_allele} at {}",
                    location(header, rec)
                );
                out.push(b'.');
                l += 1;
                continue;
            }
            out.extend_from_slice(rec.allele(ial as usize));
        } else {
            out.push(b'.');
        }
        l += 1;
    }
    if l == 0 {
        out.push(b'.');
    }
}

fn process_tgt(
    f: &mut Field,
    st: &State,
    header: &Header,
    rec: &Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    if !f.ready {
        init_format(f, st, header, rec)?;
    }
    let Some(i) = f.fmt else {
        out.push(b'.');
        return Ok(());
    };
    let field = &rec.format_fields()[i];
    if field.type_code() != BT_INT8 {
        eprintln!(
            "Warning: skipping GT field with unexpected type {} at {}",
            field.type_code(),
            location(header, rec)
        );
        out.push(b'.');
        return Ok(());
    }
    let x = &rec.format_sample_data(field, isample)[..field.n()];
    push_gt_alleles(header, rec, x, out);
    Ok(())
}

fn process_iupac_gt(
    f: &mut Field,
    st: &mut State,
    header: &Header,
    rec: &Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    const IUPAC: [[u8; 4]; 4] = [*b"AMRW", *b"MCSY", *b"RSGK", *b"WYKT"];
    if !f.ready {
        init_format(f, st, header, rec)?;
        if f.fmt.is_some() {
            st.iupac.clear();
            st.iupac.extend(rec.alleles().map(|a| {
                if a.len() != 1 {
                    -1
                } else {
                    i8::try_from(acgt2int(a[0])).unwrap_or(-1)
                }
            }));
        }
    }
    let Some(i) = f.fmt else {
        out.push(b'.');
        return Ok(());
    };
    let field = &rec.format_fields()[i];
    if field.type_code() != BT_INT8 {
        eprintln!(
            "Warning: skipping GT field with unexpected type {} at {}",
            field.type_code(),
            location(header, rec)
        );
        out.push(b'.');
        return Ok(());
    }
    let x = &rec.format_sample_data(field, isample)[..field.n()];
    let l = x
        .iter()
        .take_while(|&&v| v as i8 != INT8_VECTOR_END && v as i8 != INT8_MISSING)
        .count();
    if l == 2 {
        let n_allele = rec.allele_count() as i32;
        let ia = i32::from((x[0] as i8) >> 1) - 1;
        let ib = i32::from((x[1] as i8) >> 1) - 1;
        if (0..n_allele).contains(&ia) && (0..n_allele).contains(&ib) {
            let (da, db) = (st.iupac[ia as usize], st.iupac[ib as usize]);
            if da >= 0 && db >= 0 {
                out.push(IUPAC[da as usize][db as usize]);
                return Ok(());
            }
        }
    }
    push_gt_alleles(header, rec, x, out);
    Ok(())
}

fn process_tbcsq(
    f: &mut Field,
    st: &State,
    header: &Header,
    rec: &mut Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let Extra::Tbcsq(csq) = &mut f.extra else {
        unreachable!("TBCSQ fields carry their state");
    };
    if !f.ready {
        // init_format, inlined for the borrow of the field's state.
        resolve_format_id(&mut f.id, f.key.as_deref(), header);
        let key = f.key.as_deref().unwrap_or(b"");
        f.fmt = None;
        if f.id >= 0 {
            f.fmt = rec.format_fields().iter().position(|x| x.id() == f.id);
        } else if !st.allow_undef_tags {
            bail!(
                "Error: no such tag defined in the VCF header: FORMAT/{}",
                String::from_utf8_lossy(key)
            );
        }
        f.ready = true;
        csq.items.clear();
        csq.text.clear();
        if let Ok(s) = rec.get_info_string(header, key) {
            csq.text.extend_from_slice(s);
            let mut i = 0;
            let len = csq.text.len();
            loop {
                let start = i;
                while i < len && csq.text[i] != b',' {
                    i += 1;
                }
                csq.items.push((start, i));
                if i < len {
                    i += 1;
                }
                if i >= len {
                    break;
                }
            }
        }
    }
    let Some(fi) = f.fmt else {
        return Ok(());
    };
    if csq.items.is_empty() {
        return Ok(());
    }
    csq.hap1.clear();
    csq.hap2.clear();
    let mask: u32 = if f.subscript == 0 { 3 } else { 1 };
    let field = &rec.format_fields()[fi];
    let x = rec.format_sample_data(field, isample);
    let (nbits, width): (u32, usize) = match field.type_code() {
        BT_INT8 => (8, 1),
        BT_INT16 => (16, 2),
        BT_INT32 => (30, 4),
        ty => bail!("Unexpected type: {ty}"),
    };
    let value = |j: usize| -> u32 {
        match width {
            1 => u32::from(x[j]),
            2 => u32::from(u16::from_le_bytes([x[2 * j], x[2 * j + 1]])),
            _ => le_u32(&x[4 * j..]),
        }
    };
    let push_item = |dst: &mut Vec<u8>, idx: usize| {
        if let Some(&(a, b)) = csq.items.get(idx) {
            dst.extend_from_slice(&csq.text[a..b]);
            dst.push(b',');
        }
    };
    if f.subscript <= 0 || f.subscript == 1 {
        for j in 0..field.n() {
            let val = value(j);
            if val == 0 {
                continue;
            }
            for i in (0..nbits).step_by(2) {
                if val & (mask << i) != 0 {
                    push_item(&mut csq.hap1, (j * 30 + i as usize) / 2);
                }
            }
        }
    }
    if f.subscript < 0 || f.subscript == 2 {
        for j in 0..field.n() {
            let val = value(j);
            if val == 0 {
                continue;
            }
            for i in (1..nbits).step_by(2) {
                if val & (1 << i) != 0 {
                    push_item(&mut csq.hap2, (j * 30 + i as usize) / 2);
                }
            }
        }
    }
    if csq.hap1.is_empty() && csq.hap2.is_empty() {
        return Ok(());
    }
    csq.hap1.pop();
    csq.hap2.pop();
    let put = |out: &mut Vec<u8>, s: &[u8]| {
        if s.is_empty() {
            out.push(b'.');
        } else {
            out.extend_from_slice(s);
        }
    };
    if f.subscript < 0 {
        put(out, &csq.hap1);
        out.push(b'\t');
        put(out, &csq.hap2);
    } else if f.subscript < 2 {
        put(out, &csq.hap1);
    } else {
        put(out, &csq.hap2);
    }
    Ok(())
}

fn process_pbinom(
    f: &mut Field,
    header: &Header,
    rec: &Record,
    isample: usize,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let Extra::Pbinom { gt } = &mut f.extra else {
        unreachable!("PBINOM fields carry their state");
    };
    if !f.ready {
        let fmts = rec.format_fields();
        f.fmt = fmts.iter().position(|x| x.id() == f.id);
        let gt_id = header_id(header, b"GT");
        if !id_exists(header, HeaderRecordType::Format, gt_id) {
            bail!("Error: FORMAT/GT is not defined in the header");
        }
        *gt = fmts
            .iter()
            .position(|x| x.id() == gt_id)
            .filter(|&i| fmts[i].type_code() == BT_INT8);
        f.ready = true;
    }
    let (Some(ad_i), Some(gt_i)) = (f.fmt, *gt) else {
        out.push(b'.');
        return Ok(());
    };
    let fmts = rec.format_fields();
    let (ad, gtf) = (&fmts[ad_i], &fmts[gt_i]);
    if gtf.n() != 2 {
        out.push(b'.');
        return Ok(());
    }
    let gts = rec.format_sample_data(gtf, isample);
    let ad_data = rec.format_data(ad);
    let mut n = [0i32; 2];
    for i in 0..2 {
        let g = gts[i] as i8;
        if g >> 1 == 0 || g == INT8_VECTOR_END {
            out.push(b'.');
            return Ok(());
        }
        let al = i32::from(g >> 1) - 1;
        if al < 0 || al as usize >= rec.allele_count() || al as usize >= ad.n() {
            out.push(b'.');
            return Ok(());
        }
        let idx = al as usize + isample * ad.n();
        let val = match ad.type_code() {
            BT_INT8 | BT_INT16 | BT_INT32 => array_ivalue(ad_data, ad.type_code(), idx),
            _ => INT32_MISSING,
        };
        if val == INT32_MISSING || val == INT32_VECTOR_END {
            out.push(b'.');
            return Ok(());
        }
        n[i] = val;
    }
    if n[0] == n[1] {
        out.push(if n[0] == 0 { b'.' } else { b'0' });
    } else {
        let mut pval = kfunc::binom_two_sided(n[0], n[1], 0.5);
        // As a phred score.
        pval = if pval >= 1.0 {
            0.0
        } else if pval.is_nan() {
            // Negative counts give NaN, which bcftools' negation prints as `nan`.
            f64::NAN
        } else {
            -4.34294481903 * pval.ln()
        };
        push_double(out, pval);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> Header {
        Header::parse(
            b"##fileformat=VCFv4.2\n\
##contig=<ID=1>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        )
        .expect("the header parses")
    }

    fn record(h: &mut Header, line: &str) -> Record {
        let mut rec = Record::new();
        let mut l = line.as_bytes().to_vec();
        pasta::variant::vcf::parse(h, &mut l, &mut rec).expect("the line parses");
        rec
    }

    fn run(fmt: &str, line: &str) -> String {
        let mut h = header();
        let mut rec = record(&mut h, line);
        let mut conv = Convert::new(&h, None, fmt).expect("the format parses");
        let mut out = Vec::new();
        conv.line(&h, &mut rec, None, &mut out)
            .expect("the line formats");
        String::from_utf8(out).expect("the output is text")
    }

    const LINE: &str =
        "1\t100\trs5\tA\tG,T\t30\tPASS\tDP=7;AF=0.25,0.5\tGT:AD\t0/1:3,4,0\t1|2:.,1,9";

    #[test]
    fn fixed_and_per_sample_fields() {
        assert_eq!(
            run(
                "%CHROM:%POS %REF>%ALT{1} %DP %AF{1}[ %SAMPLE=%GT/%TGT/%AD{1}]",
                LINE
            ),
            "1:100 A>T 7 0.5 A=0/1/A/G/4 B=1|2/G|T/1"
        );
        assert_eq!(
            run("%INFO %FILTER %QUAL %END %POS0 %TYPE", LINE),
            "DP=7;AF=0.25,0.5 PASS 30 100 99 SNP"
        );
        assert_eq!(run("[%IUPACGT ]", LINE), "R K ");
    }

    #[test]
    fn header_line_numbers_columns() {
        let h = header();
        let mut conv = Convert::new(&h, None, "%CHROM[\t%GT]\n").expect("the format parses");
        conv.set_header_samples();
        let mut out = Vec::new();
        conv.header_line(&h, &mut out);
        assert_eq!(out, b"#[1]CHROM\t[2]A:GT\t[3]B:GT\n");
    }

    #[test]
    fn copy_field_picks_items() {
        let mut out = Vec::new();
        copy_field(b"a,bb,,c", 1, &mut out);
        copy_field(b"a,bb,,c", 2, &mut out);
        copy_field(b"a,bb,,c", 7, &mut out);
        assert_eq!(out, b"bb..");
    }
}
