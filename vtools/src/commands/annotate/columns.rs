//! Parsing `-c/--columns` (and `-C/--columns-file`, `-l/--merge-logic`) into the columns to
//! transfer (bcftools' `init_columns` and `init_merge_method`).

use std::collections::HashSet;

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecord, HeaderRecordType, Number, Quoting, ValueType,
    readlist::{self, ListSource},
    synced::Pairing,
};

use super::samples::{self, SampleMap};

/// How a column replaces values (bcftools' `REPLACE_*`, `SET_OR_APPEND`, `MATCH_VALUE`, and
/// `CARRY_OVER_MISSING` bits).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub(super) struct Replace(u8);

impl Replace {
    /// `+TAG`: replace only missing values.
    pub(super) const MISSING: Self = Self(1);
    /// `TAG`: replace missing and existing values.
    pub(super) const ALL: Self = Self(2);
    /// `-TAG`: replace only existing values, possibly with missing ones.
    pub(super) const NON_MISSING: Self = Self(4);
    /// `=TAG`: set if missing, append otherwise.
    pub(super) const SET_OR_APPEND: Self = Self(8);
    /// `~TAG`: only match the value, do not set it.
    pub(super) const MATCH_VALUE: Self = Self(16);
    /// `.TAG`: transfer missing source values too.
    pub(super) const CARRY_OVER_MISSING: Self = Self(32);

    /// Whether any of the bits of `other` is set.
    pub(super) fn has(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }

    /// Whether `.` (carry over missing) is combined with replacing (`TAG` or `-TAG`): a missing
    /// source value then clears the destination.
    pub(super) fn clears_with_missing(self) -> bool {
        self.has(Self::CARRY_OVER_MISSING) && self.has(Self::ALL | Self::NON_MISSING)
    }
}

impl std::ops::BitOr for Replace {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

/// How values of several overlapping annotation lines combine (`-l TAG:TYPE`, bcftools'
/// `MM_*`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub(super) enum MergeMethod {
    /// Use the first line only.
    #[default]
    First,
    /// Append all values.
    Append,
    /// Append unique values (strings).
    Unique,
    Sum,
    Avg,
    Min,
    Max,
    /// Append all values, missing ones too.
    AppendMissing,
}

/// The kind of annotation file (`-a`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(super) enum AnnotKind {
    /// A VCF or BCF file.
    Vcf,
    /// A tab-delimited file.
    Tab,
    /// No file: annotations from the record itself, or none.
    NoFile,
}

/// What a column sets in the record (the `setter` function of bcftools' `annot_col_t`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(super) enum Setter {
    // From tab-delimited annotation lines.
    Filter,
    Pos,
    Id,
    Qual,
    InfoFlag,
    InfoInt,
    InfoReal,
    InfoStr,
    FormatInt,
    FormatReal,
    FormatStr,
    // From VCF annotation records.
    VcfFilter,
    VcfId,
    VcfRef,
    VcfAlt,
    VcfQual,
    VcfInfoFlag,
    VcfInfoInt,
    VcfInfoReal,
    VcfInfoStr,
    VcfFormatGt,
    VcfFormatInt,
    VcfFormatReal,
    VcfFormatStr,
}

/// Where a VCF string column takes its value when it is not the INFO tag itself (bcftools'
/// `getter`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(super) enum Getter {
    /// A string INFO tag of the annotation record (`ID:=TAG`).
    InfoStr,
    /// The annotation record's ID, with `;` as `,`.
    Id,
    /// The annotation record's FILTER.
    Filter,
    /// The annotated record's own FILTER (`INFO/TAG:=FILTER` without `-a`).
    FilterLocal,
}

/// One column to transfer (bcftools' `annot_col_t`).
#[derive(Debug, Default)]
pub(super) struct Column {
    /// The column of a tab-delimited line (-1 for VCF annotations).
    pub(super) icol: i32,
    pub(super) replace: Replace,
    /// The `Number=` of the tag (`BCF_VL_FIXED` 0 when not a tag).
    pub(super) number: Option<Number>,
    pub(super) key_src: Vec<u8>,
    pub(super) key_dst: Vec<u8>,
    pub(super) setter: Option<Setter>,
    pub(super) getter: Option<Getter>,
    pub(super) merge: MergeMethod,
    /// Strings appended so far with `-l TAG:unique` (bcftools' `mm_str_hash`).
    pub(super) mm_seen: HashSet<Vec<u8>>,
    /// The merged string so far (bcftools' `mm_kstr`).
    pub(super) mm_str: Vec<u8>,
    /// The merged numbers so far (bcftools' `mm_dbl`, `mm_dbl_nused`).
    pub(super) mm_dbl: Vec<f64>,
    /// The number of merged lines, for averages (bcftools' `mm_dbl_ndat`).
    pub(super) mm_ndat: usize,
    /// Whether the column is done for the current record.
    pub(super) done: bool,
    /// The getter's output.
    pub(super) buf: Vec<u8>,
}

impl Column {
    fn new(icol: i32, replace: Replace, src: &[u8], dst: &[u8], setter: Option<Setter>) -> Self {
        Self {
            icol,
            replace,
            number: Some(Number::Fixed(0)),
            key_src: src.to_vec(),
            key_dst: dst.to_vec(),
            setter,
            ..Self::default()
        }
    }

    /// Whether the tag has one value per ALT (or per allele).
    pub(super) fn is_ar(&self) -> bool {
        matches!(self.number, Some(Number::A | Number::R))
    }
}

/// The tab-delimited columns with positions and alleles; -1 when absent.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub(super) struct TabColumns {
    pub(super) chr: i32,
    pub(super) beg: i32,
    pub(super) end: i32,
    pub(super) ref_: i32,
    pub(super) alt: i32,
}

impl Default for TabColumns {
    fn default() -> Self {
        Self {
            chr: -1,
            beg: -1,
            end: -1,
            ref_: -1,
            alt: -1,
        }
    }
}

/// The result of parsing `-c`.
#[derive(Debug, Default)]
pub(super) struct Columns {
    pub(super) cols: Vec<Column>,
    pub(super) tab: TabColumns,
    /// The column matched against the record's ID (`~ID`), -1 if none.
    pub(super) match_id: i32,
    /// The column matched against INFO/END (`~INFO/END`), -1 if none.
    pub(super) match_end: i32,
    /// `NEW:=OLD` renames, as `(old, new)`.
    pub(super) renames: Vec<(Vec<u8>, Vec<u8>)>,
    pub(super) samples: SampleMap,
}

/// What `init_columns` needs to know.
pub(super) struct ColumnsSpec<'a> {
    /// The `-c` list, or the `-C` file name.
    pub(super) columns: &'a str,
    pub(super) columns_is_file: bool,
    /// `-a`, if given.
    pub(super) targets_fname: Option<&'a str>,
    /// The annotation VCF's header, for VCF annotations.
    pub(super) tgts_hdr: Option<&'a Header>,
    /// `-s`/`-S`.
    pub(super) samples: Option<(&'a str, ListSource)>,
}

/// Parses the leading `+-=.` of a column name (bcftools' `set_replace_mode`).
fn replace_mode(s: &[u8]) -> (Replace, &[u8]) {
    let mut mode = 0u8;
    let mut i = 0;
    while i < s.len() {
        match s[i] {
            b'+' => mode |= Replace::MISSING.0,
            b'-' => mode |= Replace::NON_MISSING.0,
            b'=' => mode |= Replace::SET_OR_APPEND.0,
            b'.' => mode |= Replace::CARRY_OVER_MISSING.0,
            _ => break,
        }
        i += 1;
    }
    if mode == 0 {
        mode = Replace::ALL.0;
    }
    (Replace(mode), &s[i..])
}

fn starts_with_ci(s: &[u8], prefix: &[u8]) -> bool {
    s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}

fn lossy(s: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(s)
}

/// Rewrites `^INFO/TAG`, `^FMT/TAG`, and `^TAG` items as `INFO` or `FORMAT` with the tags to
/// skip (bcftools' `columns_complement`).
fn columns_complement(
    columns: &[u8],
    skip_info: &mut Option<HashSet<Vec<u8>>>,
    skip_fmt: &mut Option<HashSet<Vec<u8>>>,
) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    let push = |out: &mut Vec<u8>, s: &[u8]| {
        if !out.is_empty() {
            out.push(b',');
        }
        out.extend_from_slice(s);
    };
    for item in columns.split(|&c| c == b',') {
        if item.first() != Some(&b'^') {
            push(&mut out, item);
            continue;
        }
        if starts_with_ci(item, b"^INFO/") {
            if skip_info.is_none() {
                *skip_info = Some(HashSet::new());
                push(&mut out, b"INFO");
            }
            skip_info
                .as_mut()
                .expect("set above")
                .insert(item[6..].to_vec());
        } else if starts_with_ci(item, b"^FORMAT/") || starts_with_ci(item, b"^FMT/") {
            let n = if starts_with_ci(item, b"^FMT/") { 5 } else { 8 };
            if skip_fmt.is_none() {
                *skip_fmt = Some(HashSet::new());
                push(&mut out, b"FORMAT");
            }
            skip_fmt
                .as_mut()
                .expect("set above")
                .insert(item[n..].to_vec());
        } else {
            if skip_info.is_none() {
                *skip_info = Some(HashSet::new());
                push(&mut out, b"INFO");
            }
            skip_info
                .as_mut()
                .expect("set above")
                .insert(item[1..].to_vec());
        }
    }
    // bcftools stops at the end of the string: a trailing comma adds no empty item.
    if columns.last() == Some(&b',') && out.last() == Some(&b',') {
        out.pop();
    }
    out
}

/// Formats a definition line with a new ID, without `IDX` (bcftools' `bcf_hrec_format_rename`).
fn format_renamed(hrec: &HeaderRecord, tag: &[u8]) -> Vec<u8> {
    let mut out = b"##".to_vec();
    out.extend_from_slice(hrec.key());
    out.extend_from_slice(b"=<");
    let mut n = 0;
    for (k, v) in hrec.fields() {
        if k == b"IDX" {
            continue;
        }
        if n > 0 {
            out.push(b',');
        }
        out.extend_from_slice(k);
        out.push(b'=');
        out.extend_from_slice(if k == b"ID" { tag } else { v });
        n += 1;
    }
    out.push(b'>');
    out
}

/// Formats a header line without `IDX` and without the newline (bcftools' `bcf_hrec_format`).
pub(super) fn format_line(hrec: &HeaderRecord) -> Vec<u8> {
    let mut line = hrec.to_vcf_line();
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    line
}

/// Appends a line to the output header, failing as bcftools does.
pub(super) fn append_line(hdr_out: &mut Header, line: &[u8]) -> anyhow::Result<()> {
    if hdr_out.append(line).is_err() {
        bail!("[vcfannotate.c] failed to update the header");
    }
    Ok(())
}

/// The two families of bcftools' setters: `vcf_setter_*` copy values from a record of the
/// annotation VCF, `setter_*` parse them from a tab-delimited line (or the record itself).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SetterFamily {
    Vcf,
    Tab,
}

impl SetterFamily {
    /// The family for annotations with the header `tgts` (a VCF) or without one.
    fn of(tgts: Option<&Header>) -> Self {
        if tgts.is_some() { Self::Vcf } else { Self::Tab }
    }
}

/// The setter of an INFO tag of a type.
fn info_setter(ty: Option<ValueType>, family: SetterFamily, name: &[u8]) -> anyhow::Result<Setter> {
    use SetterFamily::{Tab, Vcf};
    Ok(match (ty, family) {
        (Some(ValueType::Flag), Vcf) => Setter::VcfInfoFlag,
        (Some(ValueType::Flag), Tab) => Setter::InfoFlag,
        (Some(ValueType::Integer), Vcf) => Setter::VcfInfoInt,
        (Some(ValueType::Integer), Tab) => Setter::InfoInt,
        (Some(ValueType::Float), Vcf) => Setter::VcfInfoReal,
        (Some(ValueType::Float), Tab) => Setter::InfoReal,
        (Some(ValueType::String), Vcf) => Setter::VcfInfoStr,
        (Some(ValueType::String), Tab) => Setter::InfoStr,
        (None, _) => bail!("The type of {} not recognised (15)", lossy(name)),
    })
}

/// The setter of a FORMAT tag of a type (not GT).
fn format_setter(
    ty: Option<ValueType>,
    family: SetterFamily,
    name: &[u8],
) -> anyhow::Result<Setter> {
    use SetterFamily::{Tab, Vcf};
    Ok(match (ty, family) {
        (Some(ValueType::Integer), Vcf) => Setter::VcfFormatInt,
        (Some(ValueType::Integer), Tab) => Setter::FormatInt,
        (Some(ValueType::Float), Vcf) => Setter::VcfFormatReal,
        (Some(ValueType::Float), Tab) => Setter::FormatReal,
        (Some(ValueType::String), Vcf) => Setter::VcfFormatStr,
        (Some(ValueType::String), Tab) => Setter::FormatStr,
        (ty, _) => bail!(
            "The type of {} not recognised ({})",
            lossy(name),
            ty.map_or(15, |t| t as u8)
        ),
    })
}

fn info_type(h: &Header, name: &[u8]) -> (Option<ValueType>, Option<Number>) {
    match h.id(name) {
        Some(id) => (
            h.id_type(HeaderRecordType::Info, id),
            h.id_number(HeaderRecordType::Info, id),
        ),
        None => (None, None),
    }
}

fn format_type(h: &Header, name: &[u8]) -> (Option<ValueType>, Option<Number>) {
    match h.id(name) {
        Some(id) => (
            h.id_type(HeaderRecordType::Format, id),
            h.id_number(HeaderRecordType::Format, id),
        ),
        None => (None, None),
    }
}

fn defined(h: &Header, kind: HeaderRecordType, name: &[u8]) -> bool {
    h.id(name).is_some_and(|id| h.id_exists(kind, id))
}

/// Reads `-C FILE`: one column per line, with an optional merge logic after whitespace,
/// which is added to `merge_logic`.
fn read_columns_file(path: &str, merge_logic: &mut Vec<u8>) -> anyhow::Result<Vec<u8>> {
    let Ok(lines) = readlist::read_list(path.as_bytes(), ListSource::File) else {
        bail!("Could not parse {path}");
    };
    let mut out = Vec::new();
    for line in &lines {
        let name_end = line
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(line.len());
        let name = &line[..name_end];
        let rest = line[name_end..].trim_ascii_start();
        if !rest.is_empty() {
            if !merge_logic.is_empty() {
                merge_logic.push(b',');
            }
            merge_logic.extend_from_slice(name);
            merge_logic.push(b':');
            merge_logic.extend_from_slice(rest);
        }
        if !out.is_empty() {
            out.push(b',');
        }
        out.extend_from_slice(name);
    }
    Ok(out)
}

/// bcftools' `init_columns`: parses the column list, adds the definitions of transferred
/// tags to `hdr_out`, and sets the pairing of VCF annotations.
#[expect(
    clippy::too_many_lines,
    reason = "a line-for-line port of init_columns"
)]
pub(super) fn init_columns(
    spec: &ColumnsSpec<'_>,
    hdr: &Header,
    hdr_out: &mut Header,
    pair_logic: &mut Option<Pairing>,
    merge_logic: &mut Vec<u8>,
) -> anyhow::Result<Columns> {
    let tgts = spec.tgts_hdr;
    let vcf = tgts.is_some();
    let targets_fname = spec.targets_fname.unwrap_or("(null)");
    let mut need_sample_map = false;
    let mut out = Columns {
        match_id: -1,
        match_end: -1,
        samples: samples::init_sample_map(tgts, hdr, spec.samples)?,
        ..Columns::default()
    };

    let mut columns = if spec.columns_is_file {
        read_columns_file(spec.columns, merge_logic)?
    } else {
        spec.columns.as_bytes().to_vec()
    };
    let mut skip_info = None;
    let mut skip_fmt = None;
    if vcf {
        columns = columns_complement(&columns, &mut skip_info, &mut skip_fmt);
    }

    let mut icol: i32 = -1;
    let items: Vec<&[u8]> = if columns.is_empty() {
        Vec::new()
    } else {
        let mut v: Vec<&[u8]> = columns.split(|&c| c == b',').collect();
        // A trailing comma ends the list without an empty item.
        if columns.last() == Some(&b',') {
            v.pop();
        }
        v
    };
    for item in items {
        let (mut replace, s) = replace_mode(item);
        icol += 1;
        let cols = &mut out.cols;
        if s.is_empty() || s == b"-" {
        } else if s.eq_ignore_ascii_case(b"CHROM") {
            out.tab.chr = icol;
        } else if s.eq_ignore_ascii_case(b"POS") {
            if replace == Replace::NON_MISSING && !vcf {
                cols.push(Column::new(icol, replace, s, s, Some(Setter::Pos)));
                out.match_end = icol;
            } else {
                out.tab.beg = icol;
            }
        } else if s.eq_ignore_ascii_case(b"FROM") || s.eq_ignore_ascii_case(b"BEG") {
            out.tab.beg = icol;
        } else if s.eq_ignore_ascii_case(b"TO") || s.eq_ignore_ascii_case(b"END") {
            out.tab.end = icol;
        } else if s.eq_ignore_ascii_case(b"REF") {
            if vcf {
                cols.push(Column::new(
                    0,
                    Replace::default(),
                    s,
                    s,
                    Some(Setter::VcfRef),
                ));
            } else {
                out.tab.ref_ = icol;
            }
        } else if s.eq_ignore_ascii_case(b"ALT") {
            if vcf {
                cols.push(Column::new(0, replace, s, s, Some(Setter::VcfAlt)));
                if pair_logic.is_none() {
                    *pair_logic = Some(Pairing::ANY);
                }
            } else {
                out.tab.alt = icol;
            }
        } else if s.eq_ignore_ascii_case(b"ID") || s.eq_ignore_ascii_case(b"~ID") {
            if replace.has(Replace::NON_MISSING) {
                bail!("Apologies, the -ID feature has not been implemented yet.");
            }
            if s[0] == b'~' {
                replace = Replace::MATCH_VALUE;
            }
            let setter = if vcf { Setter::VcfId } else { Setter::Id };
            cols.push(Column::new(icol, replace, s, s, Some(setter)));
            if replace.has(Replace::MATCH_VALUE) {
                out.match_id = icol;
                if vcf {
                    *pair_logic = Some(pair_logic.map_or(Pairing::ID, |p| p | Pairing::ID));
                }
            }
        } else if s.eq_ignore_ascii_case(b"~INFO/END") && !vcf {
            cols.push(Column::new(icol, Replace::MATCH_VALUE, s, s, None));
            out.match_end = icol;
        } else if s.eq_ignore_ascii_case(b"~POS") {
            bail!(
                "Error: the use of ~POS has been deprecated, use -POS to transfer the column POS."
            );
        } else if s[0] == b'~' {
            cols.push(Column::new(icol, Replace::MATCH_VALUE, &s[1..], b"", None));
        } else if s.eq_ignore_ascii_case(b"-POS") && !vcf {
            cols.push(Column::new(icol, replace, s, s, Some(Setter::Pos)));
            out.match_end = icol;
        } else if starts_with_ci(s, b"ID:=") {
            if !vcf {
                bail!("The annotation source must be a VCF for \"{}\"", lossy(s));
            }
            if replace.has(Replace::NON_MISSING) {
                bail!("Apologies, the -ID feature has not been implemented yet.");
            }
            let src = &s[4..];
            let src = if starts_with_ci(src, b"INFO/") {
                &src[5..]
            } else {
                src
            };
            let t = tgts.expect("a VCF annotation has a header");
            let (ty, _) = info_type(t, src);
            if !defined(t, HeaderRecordType::Info, src) {
                bail!(
                    "The INFO tag \"{}\" is not defined in {targets_fname}",
                    lossy(src)
                );
            }
            if ty != Some(ValueType::String) {
                bail!("Only Type=String tags can be used to annotate the ID column");
            }
            let mut col = Column::new(icol, replace, src, &s[..2], Some(Setter::VcfId));
            col.getter = Some(Getter::InfoStr);
            cols.push(col);
        } else if let Some(at) = s.windows(2).position(|w| w == b":=")
            && (spec.targets_fname.is_none() || starts_with_ci(&s[at + 2..], b"./"))
        {
            let (lhs, rhs) = (&s[..at], &s[at + 2..]);
            if starts_with_ci(lhs, b"INFO/")
                && (rhs.eq_ignore_ascii_case(b"FILTER") || rhs.eq_ignore_ascii_case(b"./FILTER"))
            {
                // Without -a, transferring the record's own FILTER.
                let dst = &lhs[5..];
                let mut col = Column::new(icol, replace, rhs, dst, Some(Setter::VcfInfoStr));
                col.getter = Some(Getter::FilterLocal);
                let line = format!(
                    "##INFO=<ID={},Number=.,Type=String,Description=\"Transferred FILTER column\">",
                    lossy(dst)
                );
                append_line(hdr_out, line.as_bytes())?;
                col.number = info_type(hdr_out, dst).1;
                cols.push(col);
            } else {
                out.renames.push((rhs.to_vec(), lhs.to_vec()));
            }
        } else if s.eq_ignore_ascii_case(b"FILTER") {
            if replace.has(Replace::NON_MISSING) {
                bail!("Apologies, the -FILTER feature has not been implemented yet.");
            }
            let setter = if vcf {
                Setter::VcfFilter
            } else {
                Setter::Filter
            };
            cols.push(Column::new(icol, replace, s, s, Some(setter)));
            if let Some(t) = tgts {
                let lines: Vec<Vec<u8>> = t
                    .records()
                    .filter(|r| r.kind() == HeaderRecordType::Filter)
                    .map(format_line)
                    .collect();
                for line in lines {
                    append_line(hdr_out, &line)?;
                }
            }
        } else if s.eq_ignore_ascii_case(b"QUAL") {
            if replace.has(Replace::NON_MISSING) {
                bail!("Apologies, the -QUAL feature has not been implemented yet.");
            }
            if replace.has(Replace::SET_OR_APPEND) {
                bail!("Apologies, the =QUAL feature has not been implemented yet.");
            }
            let setter = if vcf { Setter::VcfQual } else { Setter::Qual };
            cols.push(Column::new(icol, replace, s, s, Some(setter)));
        } else if vcf && s.eq_ignore_ascii_case(b"INFO") {
            if replace.has(Replace::NON_MISSING) {
                bail!("Apologies, the -INFO/TAG feature has not been implemented yet.");
            }
            if replace.has(Replace::SET_OR_APPEND) {
                bail!("Apologies, the =INFO feature has not been implemented yet.");
            }
            let t = tgts.expect("a VCF annotation has a header");
            let defs: Vec<(Vec<u8>, Vec<u8>)> = t
                .records()
                .filter(|r| r.kind() == HeaderRecordType::Info)
                .filter_map(|r| Some((r.id()?.to_vec(), format_line(r))))
                .collect();
            for (id, line) in defs {
                if skip_info.as_ref().is_some_and(|k| k.contains(&id)) {
                    continue;
                }
                append_line(hdr_out, &line)?;
                let (ty, number) = info_type(hdr_out, &id);
                let mut col = Column::new(-1, replace, &id, &id, None);
                col.number = number;
                col.setter = Some(info_setter(ty, SetterFamily::Vcf, s)?);
                cols.push(col);
            }
        } else if vcf && (s.eq_ignore_ascii_case(b"FORMAT") || s.eq_ignore_ascii_case(b"FMT")) {
            let t = tgts.expect("a VCF annotation has a header");
            need_sample_map = true;
            let defs: Vec<(Vec<u8>, Vec<u8>)> = t
                .records()
                .filter(|r| r.kind() == HeaderRecordType::Format)
                .filter_map(|r| Some((r.id()?.to_vec(), format_line(r))))
                .collect();
            for (id, line) in defs {
                if skip_fmt.as_ref().is_some_and(|k| k.contains(&id)) {
                    continue;
                }
                append_line(hdr_out, &line)?;
                let mut col = Column::new(-1, replace, &id, &id, None);
                col.setter = Some(if id.eq_ignore_ascii_case(b"GT") {
                    Setter::VcfFormatGt
                } else {
                    format_setter(format_type(hdr_out, &id).0, SetterFamily::Vcf, s)?
                });
                col.number = format_type(t, &id).1;
                cols.push(col);
            }
        } else if starts_with_ci(s, b"FORMAT/") || starts_with_ci(s, b"FMT/") {
            let body = &s[if starts_with_ci(s, b"FMT/") { 4 } else { 7 }..];
            let (key_dst, key_src) = match body.windows(2).position(|w| w == b":=") {
                Some(at) => {
                    let src = &body[at + 2..];
                    let src = if starts_with_ci(src, b"FORMAT/") {
                        &src[7..]
                    } else if starts_with_ci(src, b"FMT/") {
                        &src[4..]
                    } else {
                        src
                    };
                    (&body[..at], src)
                }
                None => (body, body),
            };
            let shown = &s[..s.len() - (body.len() - key_dst.len())];
            need_sample_map = true;
            if let Some(t) = tgts {
                let Some(hrec) = t.definition(HeaderRecordType::Format, key_src) else {
                    bail!(
                        "No such annotation \"{}\" in {targets_fname}",
                        lossy(key_src)
                    );
                };
                let line = format_renamed(hrec, key_dst);
                append_line(hdr_out, &line)?;
            }
            if !defined(hdr_out, HeaderRecordType::Format, key_dst) {
                bail!(
                    "The FORMAT tag \"{}\" is not defined in {targets_fname}, was the -h option provided?",
                    lossy(shown)
                );
            }
            let mut col = Column::new(if vcf { -1 } else { icol }, replace, key_src, key_dst, None);
            if !vcf {
                icol += out.samples.nsmpl_annot as i32 - 1;
            }
            col.setter = Some(if key_src.eq_ignore_ascii_case(b"GT") {
                if !vcf {
                    bail!("The FORMAT/GT field can be currently populated only from a VCF");
                }
                Setter::VcfFormatGt
            } else {
                format_setter(
                    format_type(hdr_out, key_dst).0,
                    SetterFamily::of(tgts),
                    shown,
                )?
            });
            if let Some(t) = tgts {
                col.number = format_type(t, key_src).1;
            }
            out.cols.push(col);
        } else {
            let col = info_column(s, icol, replace, tgts, targets_fname, hdr_out)?;
            out.cols.push(col);
        }
    }
    if !need_sample_map {
        out.samples.map = None;
    } else if !out.samples.any_match {
        bail!("No matching samples in source and destination file?");
    }
    Ok(out)
}

/// The INFO column `[INFO/]DST[:=[INFO/]SRC]` of `init_columns`.
fn info_column(
    s: &[u8],
    icol: i32,
    replace: Replace,
    tgts: Option<&Header>,
    targets_fname: &str,
    hdr_out: &mut Header,
) -> anyhow::Result<Column> {
    let vcf = tgts.is_some();
    if replace.has(Replace::NON_MISSING) {
        bail!("Apologies, the -INFO/TAG feature has not been implemented yet.");
    }
    if replace.has(Replace::SET_OR_APPEND) {
        if vcf {
            bail!(
                "Error: the =INFO/TAG feature is currently supported only with TAB annotation files and has limitations\n       (the annotation type is modified to \"Number=.\" and allele ordering is disregarded)"
            );
        }
        eprintln!(
            "Warning: the =INFO/TAG feature modifies the annotation to \"Number=.\" and disregards allele ordering"
        );
    }
    let mut col = Column::new(icol, replace, b"", b"", None);
    let mut explicit_dst = false;
    let mut explicit_src = false;
    let key_dst_full: &[u8] = if starts_with_ci(s, b"INFO/") {
        explicit_dst = true;
        &s[5..]
    } else if s.eq_ignore_ascii_case(b"~INFO/END") {
        explicit_dst = true;
        &s[6..]
    } else {
        s
    };
    let (key_dst, key_src) = match key_dst_full.windows(2).position(|w| w == b":=") {
        Some(at) => {
            let mut src = &key_dst_full[at + 2..];
            if starts_with_ci(src, b"INFO/") {
                src = &src[5..];
                explicit_src = true;
            } else if starts_with_ci(src, b"FMT/") || starts_with_ci(src, b"FORMA") {
                bail!(
                    "Did you mean \"FMT/{}\" rather than \"{}\"?",
                    lossy(s),
                    lossy(s)
                );
            } else if src.eq_ignore_ascii_case(b"FILTER") && vcf {
                col.getter = Some(Getter::Filter);
            }
            (&key_dst_full[..at], src)
        }
        None => (key_dst_full, key_dst_full),
    };
    // The name shown in messages: the column up to `:=`.
    let shown = &s[..s.len() - (key_dst_full.len() - key_dst.len())];
    col.key_src = key_src.to_vec();
    col.key_dst = key_dst.to_vec();

    if !defined(hdr_out, HeaderRecordType::Info, key_dst) {
        let Some(t) = tgts else {
            bail!(
                "The INFO tag \"{}\" is not defined in {targets_fname}, was the -h option provided?",
                lossy(key_dst)
            );
        };
        let line: Vec<u8> = if key_src.eq_ignore_ascii_case(b"ID") && !explicit_src {
            format!(
                "##INFO=<ID={},Number=.,Type=String,Description=\"Transferred ID column\">",
                lossy(key_dst)
            )
            .into_bytes()
        } else if key_src.eq_ignore_ascii_case(b"FILTER") && !explicit_src {
            format!(
                "##INFO=<ID={},Number=.,Type=String,Description=\"Transferred FILTER column\">",
                lossy(key_dst)
            )
            .into_bytes()
        } else {
            let Some(hrec) = t.definition(HeaderRecordType::Info, key_src) else {
                if !explicit_dst
                    && !explicit_src
                    && t.definition(HeaderRecordType::Format, key_src).is_some()
                {
                    bail!(
                        "Did you mean \"FMT/{}\" rather than \"{}\"?",
                        lossy(shown),
                        lossy(shown)
                    );
                }
                if let Some(eq) = key_src.iter().position(|&c| c == b'=') {
                    bail!(
                        "The INFO tag \"{}\" is not defined, is this what you want \"{}:={}\" ?",
                        lossy(&key_src[..eq]),
                        lossy(&key_src[..eq]),
                        lossy(&key_src[eq + 1..])
                    );
                }
                bail!(
                    "The INFO tag \"{}\" is not defined in {targets_fname}, was the -h option provided?",
                    lossy(key_src)
                );
            };
            format_renamed(hrec, key_dst)
        };
        append_line(hdr_out, &line)?;
    }
    if vcf {
        if key_src.eq_ignore_ascii_case(b"ID") && !explicit_src {
            col.getter = Some(Getter::Id);
        } else if key_src.eq_ignore_ascii_case(b"FILTER") && !explicit_src {
            col.getter = Some(Getter::Filter);
        }
    }
    let (ty, number) = info_type(hdr_out, key_dst);
    col.number = number;
    col.setter = Some(info_setter(ty, SetterFamily::of(tgts), s)?);
    if replace.has(Replace::SET_OR_APPEND) {
        // Change the definition to Number=.
        let Some(hrec) = hdr_out.definition(HeaderRecordType::Info, key_dst) else {
            bail!(
                "Uh, could not find the new tag \"{}\" in the header",
                lossy(key_dst)
            );
        };
        let mut hrec = hrec.dup();
        let Some(j) = hrec.find_key(b"Number") else {
            bail!(
                "Uh, could not find the entry Number in the header record of {}",
                lossy(key_dst)
            );
        };
        hrec.set_field_value(j, b".", Quoting::Bare);
        hdr_out.remove(HeaderRecordType::Info, Some(key_dst));
        if hdr_out.add_record(hrec).is_err() {
            bail!("[vcfannotate.c] failed to update the header");
        }
    }
    Ok(col)
}

/// bcftools' `init_merge_method`: applies `-l TAG:TYPE,...`. `has_tab_source` is whether
/// the annotations are tab-delimited with coordinates (`tgt_idx` or `tgts`). Returns whether
/// some column appends missing values for unmatched ALTs (`has_append_mode`).
pub(super) fn init_merge_method(
    cols: &mut [Column],
    merge_logic: &[u8],
    kind: AnnotKind,
    ref_idx: i32,
) -> anyhow::Result<bool> {
    let mut has_append_mode = false;
    if merge_logic.is_empty() {
        return Ok(false);
    }
    if kind == AnnotKind::Vcf {
        bail!("Error: the --merge-logic is intended for use with BED or TAB-delimited files only.");
    }
    if kind != AnnotKind::Tab {
        bail!(
            "Error: BEG,END (or FROM,TO) columns or REF,ALT columns are expected with the --merge-logic option."
        );
    }
    let all = lossy(merge_logic).into_owned();
    // bcftools stops at the end of the string: one trailing comma adds no empty item.
    let list = merge_logic.strip_suffix(b",").unwrap_or(merge_logic);
    for item in list.split(|&c| c == b',') {
        let Some(colon) = item.iter().rposition(|&c| c == b':') else {
            bail!("Error: could not parse the argument to --merge-logic: {all}");
        };
        let (tag, ty) = (&item[..colon], &item[colon + 1..]);
        let mm = match ty.to_ascii_lowercase().as_slice() {
            b"unique" => MergeMethod::Unique,
            b"first" => MergeMethod::First,
            b"append" => MergeMethod::Append,
            b"append-missing" => {
                if ref_idx != -1 {
                    has_append_mode = true;
                }
                MergeMethod::AppendMissing
            }
            b"sum" => MergeMethod::Sum,
            b"avg" => MergeMethod::Avg,
            b"min" => MergeMethod::Min,
            b"max" => MergeMethod::Max,
            _ => bail!(
                "Error: could not parse --merge-logic {all}, the logic \"{}\" is not recognised",
                lossy(ty)
            ),
        };
        let Some(col) = cols.iter_mut().find(|c| c.key_dst == tag) else {
            bail!("No such tag in the destination file: {}", lossy(tag));
        };
        if matches!(mm, MergeMethod::Append | MergeMethod::AppendMissing)
            && col.number != Some(Number::Variable)
        {
            bail!(
                "Error: --merge-logic append can be requested only for tags of variable length (Number=.)"
            );
        }
        col.merge = mm;
    }
    Ok(has_append_mode)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn replace_modes_parse_as_bcftools() {
        assert_eq!(replace_mode(b"+TAG"), (Replace::MISSING, b"TAG".as_slice()));
        assert_eq!(replace_mode(b"TAG"), (Replace::ALL, b"TAG".as_slice()));
        let (r, s) = replace_mode(b".-TAG");
        assert!(r.has(Replace::NON_MISSING) && r.clears_with_missing());
        assert_eq!(s, b"TAG");
    }

    #[test]
    fn complements_become_classes() {
        let (mut i, mut f) = (None, None);
        let out = columns_complement(b"ID,^INFO/AC,^FMT/DP,^AN", &mut i, &mut f);
        assert_eq!(out, b"ID,INFO,FORMAT");
        assert!(i.expect("skip set").contains(b"AN".as_slice()));
    }
}
