//! Setting up `+split-vep` from the header and the options (bcftools' `init_data` and the
//! functions it calls).

use std::{
    collections::{HashMap, HashSet},
    path::Path,
};

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType,
    index::split_idx_name,
    readlist::{self, ListSource},
    synced::SyncedReader,
};

use super::{
    Annot, Args, ColType, DEFAULT_COLUMN_TYPES, DEFAULT_SEVERITY, Genes, SelectTr, Severity,
    SplitVep, TrExpr, TrOp, Warned, process,
};
use crate::raw_input::log_open_error;
use crate::{
    convert::Convert,
    filter::{Filter, FilterLogic, FilterStatus, bre},
};

/// C's `isspace`.
pub(super) fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Whether a byte can be part of a tag name (bcftools' `valid_tag`).
fn is_tag_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// Makes a subfield name a valid tag name (bcftools' `sanitize_field_name`).
fn sanitize_field_name(s: &str) -> String {
    if s == "1000G" {
        return s.to_owned();
    }
    let mut out = String::with_capacity(s.len() + 1);
    if s.starts_with('.') || s.starts_with(|c: char| c.is_ascii_digit()) {
        out.push('_');
    }
    for &c in s.as_bytes() {
        out.push(if is_tag_char(c) { char::from(c) } else { '_' });
    }
    out
}

/// Removes surrounding white space and quotes (bcftools' `trim_field_name`).
fn trim_field_name(s: &str) -> &str {
    s.trim_matches(|c: char| c.is_ascii() && (is_space(c as u8) || c == '\'' || c == '"'))
}

impl SplitVep {
    /// Prefixes and sanitizes a subfield name (bcftools' `strdup_annot_prefix`).
    fn prefixed(&self, s: &str) -> String {
        prefixed(self.annot_prefix.as_deref(), s)
    }
}

fn prefixed(prefix: Option<&str>, s: &str) -> String {
    match prefix {
        Some(p) => sanitize_field_name(&format!("{p}{s}")),
        None => sanitize_field_name(s),
    }
}

/// Whether `%FIELD` occurs in the format string as a whole name (bcftools' `query_has_field`).
fn query_has_field(fmt: &str, field: &str) -> bool {
    let needle = format!("%{field}");
    let f = fmt.as_bytes();
    let mut from = 0;
    while let Some(i) = memchr::memmem::find(&f[from..], needle.as_bytes()) {
        let at = from + i;
        match f.get(at + needle.len()) {
            Some(&c) if is_tag_char(c) => from = at + 1,
            _ => return true,
        }
    }
    false
}

/// Opens the annotation of the header and sets everything up (bcftools' `init_data`).
pub(super) fn init_data(
    args: &Args,
    reader: &mut SyncedReader,
    hdr_out: &mut Header,
) -> anyhow::Result<SplitVep> {
    let hdr = reader.header();
    let vep_tag = match &args.annotation {
        Some(a) => a.clone(),
        None => default_vep_tag(hdr)?,
    };
    let (desc, vep_format) = read_vep_format(hdr, &vep_tag)?;
    let annot_prefix = args.annot_prefix.clone();
    let (fields, csq_field) = parse_vep_fields(
        &vep_format,
        annot_prefix.as_deref(),
        args.csq_field.clone(),
        &vep_tag,
        &desc,
    )?;
    let mut field2idx = HashMap::new();
    for (i, f) in fields.iter().enumerate() {
        if field2idx.contains_key(f) {
            eprintln!("Warning: duplicate INFO/{vep_tag} key \"{f}\"");
            continue;
        }
        field2idx.insert(f.clone(), i);
    }

    let severity = read_severity(args.severity.as_deref())?;
    let mut sv = SplitVep {
        vep_tag,
        fields,
        field2idx,
        annot_prefix,
        csq_idx: 0,
        csq_field: String::new(),
        severity,
        min_severity: process::SELECT_CSQ_ANY,
        max_severity: process::SELECT_CSQ_ANY,
        prn_csq_worst: false,
        select_tr: SelectTr::All,
        annots: Vec::new(),
        column_str: args.columns.clone(),
        column_types: args.columns_types.clone(),
        raw_vep_request: false,
        filter: None,
        convert: None,
        genes: None,
        duplicate: args.duplicate,
        drop_sites: false,
        has_format: args.format.is_some(),
        warned: Warned::default(),
        csq: process::Csq::default(),
        farr: Vec::new(),
        iarr: Vec::new(),
        list_tr: Vec::new(),
        line: Vec::new(),
        max_unpack: pasta::variant::Unpack::empty(),
    };
    sv.parse_select(args.select.as_deref().unwrap_or("all:any"), &vep_format)?;

    let tmp = sv.prefixed(&csq_field);
    let Some(&csq_idx) = sv.field2idx.get(&tmp) else {
        bail!(
            "The field \"{csq_field}\" is not present in INFO/{}: {desc}",
            sv.vep_tag
        );
    };
    sv.csq_idx = csq_idx;
    sv.csq_field = csq_field;

    let mut format_str = args.format.clone();
    if let Some(f) = &mut format_str {
        sv.parse_format_str(reader.header(), f, args.all_fields.as_deref());
    }
    let filter_expr = args
        .include
        .first()
        .map(|e| (e.as_str(), FilterLogic::Include))
        .or_else(|| {
            args.exclude
                .first()
                .map(|e| (e.as_str(), FilterLogic::Exclude))
        });
    if let Some((expr, logic)) = filter_expr {
        let filter = sv.parse_filter_str(hdr_out, expr)?;
        sv.filter = Some((filter, logic));
    }
    if sv.column_str.is_some() {
        sv.parse_column_str(hdr_out)?;
    }
    if let Some(f) = &format_str {
        let mut convert = Convert::new(hdr_out, None, f)?;
        if args.allow_undef_tags {
            convert.set_allow_undef_tags();
        }
        if args.print_header > 1 {
            convert.set_no_hdr_indices();
        }
        convert.force_newline(hdr_out)?;
        sv.convert = Some(convert);
    }
    if let Some(g) = &args.gene_list {
        sv.init_gene_list(g, args.gene_list_fields.as_deref())?;
    }
    for a in &mut sv.annots {
        a.id = hdr_out.id(a.tag.as_bytes());
    }

    let mut max_unpack = pasta::variant::Unpack::SHR;
    if let Some(c) = &sv.convert {
        max_unpack |= c.max_unpack();
    }
    if let Some((f, _)) = &sv.filter {
        max_unpack |= f.max_unpack();
    }
    if format_str.is_none() {
        max_unpack |= pasta::variant::Unpack::FMT;
    }
    reader.set_max_unpack(max_unpack);
    sv.max_unpack = max_unpack;
    Ok(sv)
}

/// The annotation tag when `-a` is not given: the first of INFO/CSQ, BCSQ, and ANN in the
/// header.
fn default_vep_tag(hdr: &Header) -> anyhow::Result<String> {
    let info_exists = |name: &str| {
        hdr.id(name.as_bytes())
            .is_some_and(|id| hdr.id_exists(HeaderRecordType::Info, id))
    };
    let has = [info_exists("CSQ"), info_exists("BCSQ"), info_exists("ANN")];
    if has.iter().filter(|&&h| h).count() > 1 {
        eprintln!(
            "Warning: multiple default annotation fields exist; the order of precedence is INFO/CSQ, BCSQ, ANN"
        );
    }
    Ok(match has {
        [true, _, _] => "CSQ".to_owned(),
        [false, true, _] => "BCSQ".to_owned(),
        [false, false, true] => "ANN".to_owned(),
        _ => bail!("Error: None of the INFO/CSQ, BCSQ, or ANN annotations was found in the header"),
    })
}

/// The Description of the annotation tag and its subfield list: the header must declare the
/// subfields after "Format: " or "Functional annotations: ".
fn read_vep_format(hdr: &Header, vep_tag: &str) -> anyhow::Result<(String, String)> {
    let Some(hrec) = hdr.definition(HeaderRecordType::Info, vep_tag.as_bytes()) else {
        bail!("The tag INFO/{vep_tag} not found in the header");
    };
    let Some(desc) = hrec.get(b"Description") else {
        bail!("No \"Description\" field was found for the tag INFO/{vep_tag} in the header");
    };
    let desc = String::from_utf8_lossy(desc).into_owned();
    let format = if let Some(i) = desc.find("Format:") {
        &desc[i + "Format:".len()..]
    } else if let Some(i) = desc.find("Functional annotations:") {
        &desc[i + "Functional annotations:".len()..]
    } else {
        bail!(
            "Expected \"Format: \" or \"Functional annotations: \" substring in the header INFO/{vep_tag}/Description, found: {desc}"
        );
    };
    let format = format.to_owned();
    Ok((desc, format))
}

/// Splits the subfield list into the prefixed field names, and finds the consequence field
/// (`--csq-field`, or the first Consequence or Annotation subfield).
fn parse_vep_fields(
    format: &str,
    annot_prefix: Option<&str>,
    mut csq_field: Option<String>,
    vep_tag: &str,
    desc: &str,
) -> anyhow::Result<(Vec<String>, String)> {
    let mut fields = Vec::new();
    for part in format.split('|') {
        let mut name = trim_field_name(part);
        if let Some(paren) = name.find('(') {
            name = trim_field_name(&name[..paren]);
        }
        if name.is_empty() {
            continue;
        }
        fields.push(prefixed(annot_prefix, name));
        if csq_field.is_none() && (name == "Consequence" || name == "Annotation") {
            csq_field = Some(name.to_owned());
        }
    }
    if fields.is_empty() {
        bail!("Could not parse Description of INFO/{vep_tag}: {desc}");
    }
    let Some(csq_field) = csq_field else {
        return Err(super::Exit::error(format!(
            "Error: Use the --csq-field option, unable to auto-detect consequence field (Consequence or Annotation) in: {format}"
        )));
    };
    Ok((fields, csq_field))
}

/// Reads the severity scale: the `-S` file up to its first empty line, or the default.
/// `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
fn read_severity(path: Option<&str>) -> anyhow::Result<Severity> {
    let mut text = Vec::new();
    match path {
        Some(p) => {
            let data = split_idx_name(p).0;
            let lines = readlist::read_all_lines(Path::new(data)).map_err(|e| {
                log_open_error(data, &e);
                anyhow::anyhow!("Cannot read {p}")
            })?;
            for l in lines {
                // hts_getline returns 0 for an empty line, which ends bcftools' loop.
                if l.is_empty() {
                    break;
                }
                text.extend_from_slice(&l);
                text.push(b'\n');
            }
        }
        None => text.extend_from_slice(DEFAULT_SEVERITY.as_bytes()),
    }
    let mut sev = Severity::default();
    let mut severity = 0;
    let mut ep = 0;
    while ep < text.len() {
        if text[ep] == b'#' {
            while ep < text.len() && text[ep] != b'\n' {
                ep += 1;
            }
            if ep == text.len() {
                break;
            }
            ep += 1;
            continue;
        }
        let bp = ep;
        while ep < text.len() && !is_space(text[ep]) {
            text[ep] = text[ep].to_ascii_lowercase();
            ep += 1;
        }
        let token = text[bp..ep].to_vec();
        sev.map.entry(token.clone()).or_insert(severity);
        sev.scale.push(token);
        if ep == text.len() {
            break;
        }
        if text[ep] == b'\n' {
            severity += 1;
        }
        ep += 1;
        while ep < text.len() && is_space(text[ep]) {
            ep += 1;
        }
    }
    Ok(sev)
}

impl SplitVep {
    /// Parses `-s TR:CSQ:PRN`.
    fn parse_select(&mut self, select: &str, vep_format: &str) -> anyhow::Result<()> {
        let cols: Vec<&str> = select.split(':').collect();
        let sel_tr = Some(cols[0]).filter(|s| !s.is_empty()).unwrap_or("all");
        let sel_csq = cols
            .get(1)
            .copied()
            .filter(|s| !s.is_empty())
            .unwrap_or("any");
        let prn_csq = cols
            .get(2)
            .copied()
            .filter(|s| !s.is_empty())
            .unwrap_or("all");

        self.select_tr = if sel_tr.eq_ignore_ascii_case("all") {
            SelectTr::All
        } else if sel_tr.eq_ignore_ascii_case("worst") {
            SelectTr::Worst
        } else if sel_tr.eq_ignore_ascii_case("primary") {
            SelectTr::Expr(self.tr_expr("CANONICAL=YES", vep_format)?)
        } else if sel_tr.eq_ignore_ascii_case("pick") {
            SelectTr::Expr(self.tr_expr("PICK=1", vep_format)?)
        } else if sel_tr.eq_ignore_ascii_case("mane") {
            SelectTr::Expr(self.tr_expr("MANE_SELECT!=\"\"", vep_format)?)
        } else {
            SelectTr::Expr(self.tr_expr(sel_tr, vep_format)?)
        };

        if !sel_csq.eq_ignore_ascii_case("any") {
            let (name, modifier) = if let Some(n) = sel_csq.strip_suffix('+') {
                (n, '+')
            } else if let Some(n) = sel_csq.strip_suffix('-') {
                (n, '-')
            } else {
                (sel_csq, '=')
            };
            let key = name.to_ascii_lowercase();
            let Some(&severity) = self.severity.map.get(key.as_bytes()) else {
                bail!(
                    "Error: the consequence \"{name}\" is not recognised. Run \"vtools +split-vep -S ?\" to see the default list."
                );
            };
            (self.min_severity, self.max_severity) = match modifier {
                '+' => (severity, i32::MAX),
                '-' => (0, severity),
                _ => (severity, severity),
            };
        }

        if prn_csq.eq_ignore_ascii_case("all") {
            self.prn_csq_worst = false;
        } else if prn_csq.eq_ignore_ascii_case("worst") {
            self.prn_csq_worst = true;
        } else {
            bail!("Error: could not parse \"{prn_csq}\" in the expression \"{select}\"");
        }
        Ok(())
    }

    /// Parses a `-s FIELD<OP>VALUE` transcript expression (bcftools' `init_select_tr_expr`).
    fn tr_expr(&self, expr: &str, vep_format: &str) -> anyhow::Result<TrExpr> {
        let b = expr.as_bytes();
        let mut found = None;
        for p in 0..b.len() {
            let op_len = match (b[p], b.get(p + 1)) {
                (b'=', _) | (b'~', _) => 1,
                (b'!', Some(b'=' | b'~')) => 2,
                _ => continue,
            };
            found = Some((p, op_len));
            break;
        }
        let Some((p, op_len)) = found else {
            bail!("Could not parse the expression: -s {expr}");
        };
        let field = self.prefixed(&expr[..p]);
        let rest = &b[p..];
        let quoted = rest.get(op_len) == Some(&b'"') && rest.last() == Some(&b'"');
        let start = (p + op_len + usize::from(quoted)).min(b.len());
        let mut value = b[start..].to_vec();
        if quoted {
            value.pop();
        }
        let regex = |v: &[u8]| {
            bre::compile(v, bre::Case::Sensitive).map_err(|_| {
                anyhow::anyhow!(
                    "Error: fail to compile the regular expression \"{}\"",
                    String::from_utf8_lossy(v)
                )
            })
        };
        let op = match (b[p], op_len) {
            (b'=', _) => TrOp::Eq(value),
            (b'~', _) => TrOp::Re(regex(&value)?),
            (_, _) if b[p + 1] == b'=' => TrOp::Ne(value),
            _ => TrOp::NotRe(regex(&value)?),
        };
        let Some(&idx) = self.field2idx.get(&field) else {
            bail!(
                "The field \"{field}\" was requested via \"{expr}\" but it is not present in INFO/{}: {vep_format}",
                self.vep_tag
            );
        };
        Ok(TrExpr { field, idx, op })
    }

    /// Finds the subfields named in the `-f` format string, which become columns (bcftools'
    /// `parse_format_str`); with `-A`, first replaces `%TAG` by all subfields.
    fn parse_format_str(&mut self, hdr: &Header, format: &mut String, all_fields: Option<&str>) {
        if let Some(delim) = all_fields {
            let delim = if delim.eq_ignore_ascii_case("tab") {
                "\t".to_owned()
            } else if delim.eq_ignore_ascii_case("space") {
                " ".to_owned()
            } else {
                format!("\\{delim}")
            };
            self.expand_csq_expression(format, &delim);
        }
        for f in &self.fields {
            if !query_has_field(format, f) {
                continue;
            }
            let is_info = hdr
                .id(f.as_bytes())
                .is_some_and(|id| hdr.id_exists(HeaderRecordType::Info, id));
            if is_info {
                eprintln!(
                    "Note: ambiguous key %{f}; using the {f} subfield of {}, not the INFO/{f} tag",
                    self.vep_tag
                );
            }
            let c = self.column_str.get_or_insert_with(String::new);
            if !c.is_empty() {
                c.push(',');
            }
            c.push_str(f);
        }
        if query_has_field(format, &self.vep_tag) {
            self.raw_vep_request = true;
        }
    }

    /// Replaces the first `%TAG` of the format with all subfields (bcftools'
    /// `expand_csq_expression`).
    fn expand_csq_expression(&self, format: &mut String, delim: &str) {
        let needle = format!("%{}", self.vep_tag);
        let Some(at) = format.find(&needle) else {
            return;
        };
        let end = at + needle.len();
        if format.as_bytes().get(end).is_some_and(|&c| is_tag_char(c)) {
            return;
        }
        let mut s = String::with_capacity(format.len() + self.fields.len() * 16);
        s.push_str(&format[..at]);
        for (i, f) in self.fields.iter().enumerate() {
            if i > 0 {
                s.push_str(delim);
            }
            s.push('%');
            s.push_str(f);
        }
        s.push_str(&format[end..]);
        *format = s;
    }

    /// Compiles the `-i/-e` expression; subfields it names that are not tags become columns
    /// (bcftools' `parse_filter_str`).
    fn parse_filter_str(&mut self, hdr_out: &mut Header, expr: &str) -> anyhow::Result<Filter> {
        let mut filter = Filter::parse(hdr_out, expr)?;
        let status = filter.status();
        if status != FilterStatus::OK {
            if status != FilterStatus::UNKNOWN_TAGS {
                bail!("Could not parse the expression: {expr}");
            }
            let mut c = self.column_str.clone().unwrap_or_default();
            for tag in filter.undef_tags() {
                if !self.field2idx.contains_key(tag) {
                    bail!(
                        "Error: the tag \"{tag}\" is not defined in the VCF header or in INFO/{}",
                        self.vep_tag
                    );
                }
                if !c.is_empty() {
                    c.push(',');
                }
                c.push_str(tag);
            }
            self.column_str = Some(c);
            self.parse_column_str(hdr_out)?;
            filter = Filter::init(hdr_out, expr)?;
        }
        let raw = format!("INFO/{}", self.vep_tag);
        if filter.used_tags().contains(&raw) {
            self.raw_vep_request = true;
        }
        Ok(filter)
    }

    /// Resolves the `-c` columns and adds their header lines (bcftools' `parse_column_str`).
    fn parse_column_str(&mut self, hdr_out: &mut Header) -> anyhow::Result<()> {
        if !self.annots.is_empty() {
            return Ok(());
        }
        let nfield = self.fields.len();
        let all = format!("0-{}", nfield as i64 - 1);
        let mut spec = self.column_str.clone().unwrap_or_default();
        if spec == "-" {
            spec.clone_from(&all);
        }
        let mut columns: Vec<(usize, Option<ColType>)> = Vec::new();
        let mut tokens: Vec<String> = split_tokens(&spec);
        let mut k = 0;
        while k < tokens.len() {
            let bp = tokens[k].clone();
            k += 1;
            if bp == "-" {
                let mut rest = vec![all.clone()];
                rest.extend(tokens.drain(k..));
                tokens = rest;
                k = 0;
                continue;
            }
            let tmp = self.prefixed(&bp);
            let (beg, end, ty) =
                if let Some(&i) = self.field2idx.get(&bp).or_else(|| self.field2idx.get(&tmp)) {
                    (i as i64, i as i64, None)
                } else if let Some(colon) = bp.rfind(':') {
                    let name = &bp[..colon];
                    let Some(&i) = self.field2idx.get(name) else {
                        bail!("No such column: \"{name}\"");
                    };
                    let Some(ty) = ColType::parse(&bp[colon + 1..]) else {
                        bail!(
                            "The type \"{}\" nor the column \"{bp}\" recognised",
                            &bp[colon + 1..]
                        );
                    };
                    (i as i64, i as i64, Some(ty))
                } else {
                    let b = bp.as_bytes();
                    let (beg, mut mp) = pasta::variant::number::strtol(b, 0);
                    let mut end = -1;
                    let mut ty = None;
                    if mp == b.len() {
                        end = beg;
                    } else if b[mp] == b'-' {
                        let (e, m) = pasta::variant::number::strtol(b, mp + 1);
                        end = e;
                        mp = m;
                    }
                    if mp < b.len() {
                        if b[mp] == b':' {
                            end = beg;
                            let t = &bp[mp + 1..];
                            let Some(t2) = ColType::parse(t) else {
                                bail!("The type \"{t}\" nor the column \"{bp}\" recognised");
                            };
                            ty = Some(t2);
                        } else if bp == self.vep_tag {
                            self.raw_vep_request = true;
                            continue;
                        } else {
                            bail!("No such column: \"{bp}\"");
                        }
                    }
                    if beg < 0 || end < beg {
                        bail!("No such column: \"{bp}\"");
                    }
                    (beg, end, ty)
                };
            for j in beg..=end {
                if j >= nfield as i64 {
                    bail!("The index is too big: {j}");
                }
                columns.push((j as usize, ty));
            }
        }
        // Prune duplicates, keeping the first.
        let mut seen = HashSet::new();
        columns.retain(|(j, _)| seen.insert(self.fields[*j].clone()));

        let mut types: Option<Vec<(crate::filter::bre::Bre, ColType, String)>> = None;
        let mut names = Vec::new();
        for (j, ty) in columns {
            let field = self.fields[j].clone();
            let (ty, ty_name) = match ty {
                Some(t) => (t, t.header_name().to_owned()),
                None => {
                    if types.is_none() {
                        types = Some(self.column_types()?);
                    }
                    let types = types.as_ref().expect("column types were read");
                    match types
                        .iter()
                        .find(|(re, _, _)| re.is_match(field.as_bytes()))
                    {
                        Some((_, t, name)) => (*t, name.clone()),
                        None => (ColType::Str, "String".to_owned()),
                    }
                }
            };
            let line = format!(
                "##INFO=<ID={field},Number=.,Type={ty_name},Description=\"The {field} field from INFO/{}\">",
                self.vep_tag
            );
            let _ = hdr_out.append(line.as_bytes());
            names.push(field.clone());
            self.annots.push(Annot {
                tag: field,
                idx: Some(j),
                ty,
                id: None,
                value: Vec::new(),
            });
        }
        self.column_str = (!names.is_empty()).then(|| names.join(","));
        if self.raw_vep_request && matches!(self.select_tr, SelectTr::All) {
            self.raw_vep_request = false;
        }
        if self.raw_vep_request {
            self.annots.push(Annot {
                tag: self.vep_tag.clone(),
                idx: None,
                ty: ColType::Str,
                id: None,
                value: Vec::new(),
            });
        }
        Ok(())
    }

    /// Reads the subfield types of `--columns-types` or the defaults (bcftools'
    /// `init_column2type`).
    fn column_types(&self) -> anyhow::Result<Vec<(crate::filter::bre::Bre, ColType, String)>> {
        let lines: Vec<Vec<u8>> = match self.column_types.as_deref() {
            Some(f) if f != "-" => readlist::read_lines_or_inline(f.as_bytes())
                .map_err(|_| anyhow::anyhow!("Failed to parse the column types"))?,
            _ => DEFAULT_COLUMN_TYPES
                .lines()
                .map(|l| l.as_bytes().to_vec())
                .collect(),
        };
        if lines.is_empty() {
            bail!("Failed to parse the column types");
        }
        let mut out = Vec::new();
        for line in &lines {
            if line.first() == Some(&b'#') {
                continue;
            }
            let text = String::from_utf8_lossy(line);
            let bad = || anyhow::anyhow!("Error: failed to parse the column type \"{text}\"");
            let name_end = line.iter().position(|&c| is_space(c)).ok_or_else(bad)?;
            let mut p = name_end;
            while p < line.len() && is_space(line[p]) {
                p += 1;
            }
            if p == line.len() {
                return Err(bad());
            }
            let mut pattern = b"^".to_vec();
            pattern.extend_from_slice(&line[..name_end]);
            pattern.push(b'$');
            let re = bre::compile(&pattern, bre::Case::Sensitive).map_err(|_| {
                anyhow::anyhow!(
                    "Error: fail to compile the column type regular expression \"{}\": {text}",
                    String::from_utf8_lossy(&pattern)
                )
            })?;
            let ty_name = String::from_utf8_lossy(&line[p..]).into_owned();
            let ty = match ty_name.as_str() {
                "Float" => ColType::Real,
                "Integer" => ColType::Int,
                "Flag" => ColType::Flag,
                "String" => ColType::Str,
                _ => bail!("Error: the column type \"{ty_name}\" is not supported: {text}"),
            };
            out.push((re, ty, ty_name));
        }
        if out.is_empty() {
            bail!("Failed to parse the column types");
        }
        Ok(out)
    }

    /// Reads `-g [+]FILE` (bcftools' `init_gene_list`).
    fn init_gene_list(&mut self, spec: &str, fields: Option<&str>) -> anyhow::Result<()> {
        let (prioritize, fname) = match spec.strip_prefix('+') {
            Some(f) => (true, f),
            None => (false, spec),
        };
        let lines = readlist::read_lines_or_inline(fname.as_bytes()).unwrap_or_default();
        if lines.is_empty() {
            bail!("Could not read the file provided with --gene-list {fname}");
        }
        let names: HashSet<Vec<u8>> = lines.into_iter().collect();
        let fields_str = fields.unwrap_or("SYMBOL,Gene,gene,Gene_Name,Gene_ID");
        let mut idx = Vec::new();
        for f in readlist::read_list(fields_str.as_bytes(), ListSource::Inline)? {
            let tmp = self.prefixed(&String::from_utf8_lossy(&f));
            if let Some(&i) = self.field2idx.get(&tmp) {
                idx.push(i);
            }
        }
        if idx.is_empty() {
            bail!(
                "None of the \"{fields_str}\" fields is present in INFO/{}",
                self.vep_tag
            );
        }
        self.genes = Some(Genes {
            names,
            fields: idx,
            prioritize,
        });
        Ok(())
    }
}

/// Splits a `-c` list at commas, without the empty item after a trailing comma.
fn split_tokens(s: &str) -> Vec<String> {
    let mut tokens: Vec<String> = s.split(',').map(str::to_owned).collect();
    if tokens.last().is_some_and(String::is_empty) {
        tokens.pop();
    }
    tokens
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_names_are_sanitized() {
        assert_eq!(sanitize_field_name("1000G"), "1000G");
        assert_eq!(sanitize_field_name("5prime"), "_5prime");
        assert_eq!(sanitize_field_name("GERP++_RS"), "GERP___RS");
        assert_eq!(trim_field_name(" 'Allele\" "), "Allele");
    }

    #[test]
    fn whole_names_in_format_strings() {
        assert!(query_has_field("%CHROM %AF\\n", "AF"));
        assert!(!query_has_field("%CHROM %AF_EUR\\n", "AF"));
        assert!(query_has_field("%AF_EUR %AF", "AF"));
    }

    #[test]
    fn severity_lines_and_comments() {
        let sev = read_severity(None).expect("default scale");
        assert_eq!(sev.map.get(&b"intergenic"[..]), Some(&0));
        assert_eq!(sev.map.get(&b"tfbs"[..]), Some(&3));
        assert_eq!(sev.map.get(&b"transcript_ablation"[..]), Some(&19));
    }

    #[test]
    fn column_lists() {
        assert_eq!(split_tokens("a,,b,"), ["a", "", "b"]);
        assert!(split_tokens("").is_empty());
    }
}
