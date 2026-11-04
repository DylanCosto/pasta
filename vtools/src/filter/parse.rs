//! Compiling an expression into reverse Polish notation (bcftools' `filter_init_`,
//! `filters_init1`, `parse_tag_idx`, and `parse_idxs`).

use std::{
    collections::HashSet,
    sync::atomic::{AtomicBool, Ordering},
};

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Number, Unpack, ValueType, VariantType, number,
    readlist::{self, ListSource},
};

use super::{
    ErrorMode, ExtType, Filter, FilterStatus, bre,
    eval::Context,
    lexer::{Lexed, next_token},
    token::{Comparator, Fraction, Func, GtClass, LineKind, Op, Setter, Token, expand, missing},
};

/// `BCF_VL_R`.
const VL_R: u32 = 4;

/// Printed once per process, as bcftools does with a static flag.
static COMMA_SEPARATOR_WARNED: AtomicBool = AtomicBool::new(false);

/// Whether `s` (case-insensitively) is a prefix of `name`, as `strncasecmp(s, name, len(s))`
/// tests it.
fn is_prefix_ci(s: &[u8], name: &[u8]) -> bool {
    s.len() <= name.len() && name[..s.len()].eq_ignore_ascii_case(s)
}

fn starts_with_ci(s: &[u8], prefix: &[u8]) -> bool {
    s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}

fn lossy(s: &[u8]) -> String {
    String::from_utf8_lossy(s).into_owned()
}

/// Expands `~/`, `~user`, and `$VAR` in a path (bcftools' `expand_path`).
pub(super) fn expand_path(path: &[u8]) -> Vec<u8> {
    if path.first() == Some(&b'~') {
        if path.len() == 1 || path[1] == b'/' {
            let mut out =
                std::env::var_os("HOME").map_or_else(|| b".".to_vec(), |h| h.into_encoded_bytes());
            out.extend_from_slice(&path[1..]);
            return out;
        }
        let end = memchr::memchr(b'/', path).unwrap_or(path.len());
        let mut out = home_of(&path[1..end]).unwrap_or_else(|| path[..end].to_vec());
        out.extend_from_slice(&path[end..]);
        return out;
    }
    if path.first() == Some(&b'$')
        && let Ok(name) = std::str::from_utf8(&path[1..])
        && let Some(v) = std::env::var_os(name)
    {
        return v.into_encoded_bytes();
    }
    path.to_vec()
}

/// The home directory of a user (`getpwnam`).
fn home_of(user: &[u8]) -> Option<Vec<u8>> {
    let name = std::ffi::CString::new(user).ok()?;
    // SAFETY: `name` is a valid NUL-terminated string; getpwnam returns NULL or a pointer to a
    // static record whose `pw_dir` is a valid C string, which is copied before any other
    // getpw* call can overwrite it.
    unsafe {
        let pw = libc::getpwnam(name.as_ptr());
        if pw.is_null() || (*pw).pw_dir.is_null() {
            return None;
        }
        Some(std::ffi::CStr::from_ptr((*pw).pw_dir).to_bytes().to_vec())
    }
}

/// Reads a list file (`hts_readlist(fname, 1, ...)`).
fn read_file_list(fname: &[u8]) -> Option<Vec<Vec<u8>>> {
    readlist::read_list(fname, ListSource::File).ok()
}

/// Splits a function's argument list at commas outside brackets (bcftools' `parse_tag_list`).
fn parse_tag_list(s: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut start = 0;
    for (i, &c) in s.iter().enumerate() {
        if c == b',' && depth == 0 {
            out.push(&s[start..i]);
            start = i + 1;
        }
        if c == b'[' || c == b'(' {
            depth += 1;
        }
        if (c == b']' || c == b')') && depth > 0 {
            depth -= 1;
        }
    }
    out.push(&s[start..]);
    out
}

/// Parses a vector or sample subscript (bcftools' `parse_idxs`): returns `idx` (>= 0 one
/// index, -2 a list in `idxs`, -3 the GT alleles) and the list.
fn parse_idxs(s: &[u8]) -> Option<(i32, Vec<i32>)> {
    if s.is_empty() || s == b"*" {
        return Some((-2, vec![-1]));
    }
    if s == b"GT" {
        return Some((-3, vec![-1]));
    }
    let (v, end) = number::strtol(s, 0);
    if v >= 0 && end == s.len() {
        return Some((i32::try_from(v).unwrap_or(i32::MAX), Vec::new()));
    }
    let mut idxs: Vec<i32> = Vec::new();
    let mut ibeg: i64 = -1;
    let mut beg = 0;
    while beg < s.len() {
        let (num, end) = number::strtol(s, beg);
        match s.get(end) {
            Some(b',') => beg = end + 1,
            None => beg = end,
            Some(b'-') => {
                beg = end + 1;
                ibeg = num;
                continue;
            }
            Some(_) => return None,
        }
        let num = usize::try_from(num).ok()?;
        if num >= idxs.len() {
            idxs.resize(num + 1, 0);
        }
        if let Ok(b) = usize::try_from(ibeg) {
            for v in idxs.iter_mut().take(num + 1).skip(b) {
                *v = 1;
            }
            ibeg = -1;
        }
        idxs[num] = 1;
    }
    if let Ok(b) = usize::try_from(ibeg) {
        if b >= idxs.len() {
            idxs.resize(b + 1, 0);
        }
        idxs[b] = -1;
    }
    Some((-2, idxs))
}

/// The state of compilation (the parts of bcftools' `filter_t` that parsing fills).
struct Builder<'h> {
    header: &'h Header,
    expr: &'h str,
    mode: ErrorMode,
    max_unpack: Unpack,
    status: FilterStatus,
    undef_tags: Vec<String>,
    used_tags: Vec<String>,
    ext: Vec<Option<ExtType>>,
}

impl Builder<'_> {
    fn add_used_tag(&mut self, prefix: &str, tag: &[u8]) {
        let name = format!("{prefix}{}", lossy(tag));
        if !self.used_tags.contains(&name) {
            self.used_tags.push(name);
        }
    }

    fn add_undef_tag(&mut self, tag: &[u8]) {
        let name = lossy(tag);
        if !self.undef_tags.contains(&name) {
            self.undef_tags.push(name);
        }
    }

    fn id_info(&self, kind: HeaderRecordType, id: i32) -> Option<u32> {
        u32::try_from(id)
            .ok()
            .and_then(|id| self.header.id_info(kind, id))
    }

    fn exists(&self, kind: HeaderRecordType, id: i32) -> bool {
        u32::try_from(id).is_ok_and(|id| self.header.id_exists(kind, id))
    }

    fn value_type(&self, kind: HeaderRecordType, id: i32) -> Option<ValueType> {
        u32::try_from(id)
            .ok()
            .and_then(|id| self.header.id_type(kind, id))
    }

    fn is_number_one(&self, kind: HeaderRecordType, id: i32) -> bool {
        u32::try_from(id)
            .ok()
            .and_then(|id| self.header.id_number(kind, id))
            == Some(Number::Fixed(1))
    }

    fn init_usmpl(&self, tok: &mut Token) {
        if tok.nsamples > 0 {
            return;
        }
        tok.nsamples = self.header.sample_count();
        tok.usmpl = vec![true; tok.nsamples];
    }

    /// `BCF_UN_INFO` when INFO/AC and INFO/AN are defined, else `BCF_UN_FMT`
    /// (`max_ac_an_unpack`).
    fn ac_an_unpack(&self) -> Unpack {
        let defined = |name: &[u8]| {
            self.header
                .id(name)
                .is_some_and(|id| self.header.id_exists(HeaderRecordType::Info, id))
        };
        if defined(b"AC") && defined(b"AN") {
            Unpack::INFO
        } else {
            Unpack::FMT
        }
    }

    /// Reads an external-value token (`filters_init1_ext`).
    fn init_ext(&mut self, s: &[u8]) -> Token {
        let mut tok = Token::new(Op::Val);
        self.ext.push(None);
        tok.iext = self.ext.len();
        // bcftools compares only the token's length, so "{}" matches any type name.
        if is_prefix_ci(s, b"{str}") {
            tok.ht_type = Some(ValueType::String);
            tok.is_str = true;
        } else if is_prefix_ci(s, b"{int}") {
            tok.ht_type = Some(ValueType::Integer);
        } else if is_prefix_ci(s, b"{float}") {
            tok.ht_type = Some(ValueType::Float);
        }
        *self.ext.last_mut().expect("an external value was added") =
            tok.ht_type.map(ExtType::from_value_type);
        tok
    }

    /// Reads one value token: a string constant, a file, a tag, a computed value, or a number
    /// (`filters_init1`).
    fn init1(&mut self, s: &[u8]) -> anyhow::Result<Token> {
        let mut tok = Token::new(Op::Val);

        if let Some(&quote @ (b'"' | b'\'')) = s.first() {
            if s.len() < 2 || s[s.len() - 1] != quote {
                bail!("TODO: [{}]", self.expr);
            }
            let key = s[1..s.len() - 1].to_vec();
            tok.values = vec![key.len() as f64];
            tok.nvalues = key.len();
            tok.is_missing = key == b".";
            tok.key = Some(key);
            tok.is_str = true;
            tok.ht_type = Some(ValueType::String);
            return Ok(tok);
        }

        if s.first() == Some(&b'@') {
            tok.tag = Some(s.to_vec());
            let fname = expand_path(&s[1..]);
            let Some(list) = read_file_list(&fname) else {
                bail!("Could not read: {}", lossy(&fname));
            };
            let mut hash = HashSet::new();
            for mut item in list {
                if let Some(end) = item.iter().position(u8::is_ascii_whitespace) {
                    item.truncate(end);
                }
                hash.insert(item);
            }
            tok.hash = Some(hash);
            return Ok(tok);
        }

        let mut name = s;
        let mut is_fmt = None;
        if starts_with_ci(name, b"FMT/") {
            name = &name[4..];
            is_fmt = Some(true);
        } else if starts_with_ci(name, b"FORMAT/") {
            name = &name[7..];
            is_fmt = Some(true);
        } else if starts_with_ci(name, b"INFO/") {
            name = &name[5..];
            is_fmt = Some(false);
        } else if self.init_fixed_field(name, &mut tok) {
            return Ok(tok);
        }

        // Does it have an array subscript?
        let (tag, subscript) = match name.last() {
            Some(b']') => match memchr::memchr(b'[', name) {
                Some(i) => (&name[..i], Some(&name[i + 1..])),
                None => (name, None),
            },
            _ => (name, None),
        };
        tok.hdr_id = self.header.id(tag).map_or(-1, |id| id as i32);
        let is_fmt = match is_fmt {
            Some(f) => f,
            None if tok.hdr_id >= 0 => {
                let is_info = self.exists(HeaderRecordType::Info, tok.hdr_id);
                let is_fmt = self.exists(HeaderRecordType::Format, tok.hdr_id);
                if is_info && is_fmt {
                    bail!(
                        "Error: ambiguous filtering expression, both INFO/{0} and FORMAT/{0} are defined in the VCF header.",
                        lossy(tag)
                    );
                }
                is_fmt
            }
            None => false,
        };
        let kind = if is_fmt {
            HeaderRecordType::Format
        } else {
            HeaderRecordType::Info
        };
        if let Some(sub) = subscript {
            self.parse_tag_idx(kind, tag, sub, &mut tok)?;
            if tok.idx == -3 && self.vl_len(HeaderRecordType::Format, tok.hdr_id) != VL_R {
                bail!("Error: GT subscripts can be used only with Number=R tags");
            }
        } else if is_fmt {
            self.init_usmpl(&mut tok);
        }

        tok.hl_type = if is_fmt {
            LineKind::Format
        } else {
            LineKind::Info
        };
        if tok.hdr_id >= 0 {
            tok.vl_len = self.vl_len(kind, tok.hdr_id);
        }
        if is_fmt {
            self.max_unpack |= Unpack::FMT;
        }
        if tok.hdr_id >= 0 {
            self.init_tag(&mut tok, tag, kind, subscript)?;
            tok.tag = Some(tag.to_vec());
            self.add_used_tag(if is_fmt { "FORMAT/" } else { "INFO/" }, tag);
            return Ok(tok);
        }
        if self.init_computed(tag, &mut tok) {
            return Ok(tok);
        }

        // A number: integers and floats are parsed separately; floats are compared as floats.
        let (v, end) = number::strtol(tag, 0);
        if end == tag.len() {
            tok.threshold = v as f64;
            tok.ht_type = Some(ValueType::Integer);
        } else {
            let (d, end) = number::strtod(tag, 0);
            tok.threshold = d;
            if strtod_range_error(tag, d, end) || end != name.len() {
                if self.mode == ErrorMode::Exit {
                    bail!(
                        "[filter.c:3615 filters_init1] Error: the tag \"{}\" is not defined in the VCF header",
                        lossy(tag)
                    );
                }
                self.status |= FilterStatus::UNKNOWN_TAGS;
                self.add_undef_tag(tag);
            }
            tok.ht_type = Some(ValueType::Float);
        }
        tok.is_constant = true;
        Ok(tok)
    }

    /// `BCF_VL_*` of an ID (`bcf_hdr_id2length`).
    fn vl_len(&self, kind: HeaderRecordType, id: i32) -> u32 {
        self.id_info(kind, id).map_or(0, |i| (i >> 8) & 0xf)
    }

    /// The fixed columns and computed values recognised by a prefix of their name, as bcftools'
    /// `strncasecmp(str, NAME, len)` does (so `Q` is QUAL). Returns whether `name` was one.
    fn init_fixed_field(&mut self, name: &[u8], tok: &mut Token) -> bool {
        if is_prefix_ci(name, b"QUAL") {
            tok.setter = Some(Setter::Qual);
            tok.tag = Some(b"QUAL".to_vec());
            tok.ht_type = Some(ValueType::Float);
            self.add_used_tag("", b"QUAL");
        } else if is_prefix_ci(name, b"TYPE") {
            tok.setter = Some(Setter::Type);
            tok.tag = Some(b"TYPE".to_vec());
            tok.ht_type = Some(ValueType::String);
        } else if is_prefix_ci(name, b"FILTER") {
            tok.comparator = Some(Comparator::Filter);
            tok.tag = Some(b"FILTER".to_vec());
            self.max_unpack |= Unpack::FLT;
            tok.hl_type = LineKind::Filter;
            tok.ht_type = Some(ValueType::String);
            self.add_used_tag("", b"FILTER");
        } else if is_prefix_ci(name, b"ID") || is_prefix_ci(name, b"%ID") {
            tok.comparator = Some(Comparator::Id);
            tok.tag = Some(b"ID".to_vec());
            tok.ht_type = Some(ValueType::String);
            self.add_used_tag("", b"ID");
        } else if is_prefix_ci(name, b"CHROM") {
            tok.setter = Some(Setter::Chrom);
            tok.tag = Some(b"CHROM".to_vec());
            tok.ht_type = Some(ValueType::String);
            self.add_used_tag("", b"CHROM");
        } else if is_prefix_ci(name, b"POS") {
            tok.setter = Some(Setter::Pos);
            tok.tag = Some(b"POS".to_vec());
            tok.ht_type = Some(ValueType::Integer);
            self.add_used_tag("", b"POS");
        } else if is_prefix_ci(name, b"REF") {
            tok.setter = Some(Setter::RefString);
            tok.is_str = true;
            tok.tag = Some(b"REF".to_vec());
            tok.ht_type = Some(ValueType::String);
            self.add_used_tag("", b"REF");
        } else if is_prefix_ci(name, b"ALT") {
            tok.setter = Some(Setter::AltString);
            tok.is_str = true;
            tok.tag = Some(b"ALT".to_vec());
            tok.ht_type = Some(ValueType::String);
            tok.set_all_indices();
            self.add_used_tag("", b"ALT");
        } else if is_prefix_ci(name, b"N_ALT") {
            tok.setter = Some(Setter::NAlt);
            tok.tag = Some(b"N_ALT".to_vec());
            tok.ht_type = Some(ValueType::Integer);
        } else if is_prefix_ci(name, b"N_SAMPLES") {
            tok.threshold = self.header.sample_count() as f64;
            tok.is_constant = true;
            tok.ht_type = Some(ValueType::Integer);
        } else if is_prefix_ci(name, b"N_MISSING") || is_prefix_ci(name, b"F_MISSING") {
            let (fraction, tag, ht) = if is_prefix_ci(name, b"N_MISSING") {
                (Fraction::Count, &b"N_MISSING"[..], ValueType::Integer)
            } else {
                (Fraction::Fraction, &b"F_MISSING"[..], ValueType::Float)
            };
            self.max_unpack |= Unpack::FMT;
            tok.setter = Some(Setter::Missing(fraction));
            tok.tag = Some(tag.to_vec());
            tok.ht_type = Some(ht);
            self.init_usmpl(tok);
        } else {
            return false;
        }
        true
    }

    /// Sets up a token for a tag defined in the header, an INFO or FORMAT tag as `kind` says,
    /// with the subscript in its name if any.
    fn init_tag(
        &mut self,
        tok: &mut Token,
        tag: &[u8],
        kind: HeaderRecordType,
        subscript: Option<&[u8]>,
    ) -> anyhow::Result<()> {
        let id = tok.hdr_id;
        let is_fmt = kind == HeaderRecordType::Format;
        let is_array = subscript.is_some();
        if is_fmt && tag == b"GT" {
            tok.setter = Some(Setter::GenotypeString);
            tok.is_str = true;
            tok.ht_type = Some(ValueType::String);
        } else if is_fmt {
            if !self.exists(HeaderRecordType::Format, id) {
                bail!("No such FORMAT field: {}", lossy(tag));
            }
            if !self.is_number_one(HeaderRecordType::Format, id) && !is_array {
                tok.set_all_indices();
            }
            match self.value_type(HeaderRecordType::Format, id) {
                Some(ValueType::Integer) => {
                    tok.setter = Some(Setter::FormatInt);
                    tok.ht_type = Some(ValueType::Integer);
                }
                Some(ValueType::Float) => {
                    tok.setter = Some(Setter::FormatFloat);
                    tok.ht_type = Some(ValueType::Float);
                }
                Some(ValueType::String) => {
                    tok.setter = Some(Setter::FormatString);
                    tok.ht_type = Some(ValueType::String);
                    tok.is_str = true;
                }
                _ => bail!("[filter.c:3491 filters_init1] FIXME"),
            }
        } else if !self.exists(HeaderRecordType::Info, id) {
            bail!("No such INFO field: {}", lossy(tag));
        } else {
            let ty = self.value_type(HeaderRecordType::Info, id);
            if ty == Some(ValueType::Flag) {
                tok.setter = Some(Setter::InfoFlag);
                tok.ht_type = Some(ValueType::Integer);
            } else {
                tok.ht_type = ty;
                if ty == Some(ValueType::String) {
                    tok.is_str = true;
                }
                if self.is_number_one(HeaderRecordType::Info, id) {
                    tok.setter = Some(Setter::Info);
                } else {
                    tok.setter = Some(match ty {
                        Some(ValueType::Integer) => Setter::InfoInt,
                        Some(ValueType::Float) => Setter::InfoFloat,
                        Some(ValueType::String) => Setter::InfoString,
                        _ => bail!("[filter.c:3516 filters_init1] FIXME"),
                    });
                    if !is_array {
                        tok.set_all_indices();
                    }
                }
            }
            self.max_unpack |= Unpack::INFO;
        }
        Ok(())
    }

    /// The values computed from the record when no such tag is defined: ALT, AN, AC, MAC, AF,
    /// MAF, and ILEN. Returns whether `tag` was one of them.
    fn init_computed(&mut self, tag: &[u8], tok: &mut Token) -> bool {
        let upper = tag.to_ascii_uppercase();
        match upper.as_slice() {
            b"ALT" => {
                tok.setter = Some(Setter::AltString);
                tok.is_str = true;
                tok.ht_type = Some(ValueType::String);
                tok.tag = Some(tag.to_vec());
                self.add_used_tag("", tag);
            }
            b"AN" => {
                self.max_unpack |= Unpack::FMT;
                tok.setter = Some(Setter::An);
                tok.tag = Some(b"AN".to_vec());
                tok.ht_type = Some(ValueType::Integer);
            }
            b"AC" => {
                self.max_unpack |= Unpack::FMT;
                tok.setter = Some(Setter::Ac);
                tok.tag = Some(b"AC".to_vec());
                tok.idx = -2;
                tok.ht_type = Some(ValueType::Integer);
            }
            b"MAC" => {
                self.max_unpack |= self.ac_an_unpack();
                tok.setter = Some(Setter::Mac);
                tok.tag = Some(b"MAC".to_vec());
                tok.ht_type = Some(ValueType::Integer);
            }
            b"AF" => {
                self.max_unpack |= self.ac_an_unpack();
                tok.setter = Some(Setter::Af);
                tok.tag = Some(b"AF".to_vec());
                tok.idx = -2;
                tok.ht_type = Some(ValueType::Float);
            }
            b"MAF" => {
                self.max_unpack |= self.ac_an_unpack();
                tok.setter = Some(Setter::Maf);
                tok.tag = Some(b"MAF".to_vec());
                tok.ht_type = Some(ValueType::Float);
            }
            b"ILEN" | b"%ILEN" => {
                self.max_unpack |= Unpack::STR;
                tok.setter = Some(Setter::Ilen);
                tok.tag = Some(b"ILEN".to_vec());
                tok.ht_type = Some(ValueType::Integer);
            }
            _ => return false,
        }
        true
    }

    /// Parses `TAG[...]`'s subscript: vector indices, sample indices, or a file of sample
    /// names (bcftools' `parse_tag_idx`). `kind` says whether `tag` is an INFO or FORMAT tag;
    /// `sub` is what follows `[`.
    fn parse_tag_idx(
        &self,
        kind: HeaderRecordType,
        tag: &[u8],
        sub: &[u8],
        tok: &mut Token,
    ) -> anyhow::Result<()> {
        let is_fmt = kind == HeaderRecordType::Format;
        let sub = sub.strip_suffix(b"]").unwrap_or(sub);
        let ori = lossy(sub);
        let nsmpl = self.header.sample_count();
        let colon = memchr::memrchr(b':', sub);
        let bad_index = || anyhow::anyhow!("Could not parse the index: {ori}");
        let mut sample_idxs: Option<(i32, Vec<i32>)> = None;

        if sub.first() == Some(&b'@') {
            if !is_fmt {
                bail!("Could not parse \"{ori}\". (Not a FORMAT tag yet a sample list provided.)");
            }
            let mut fname = expand_path(&sub[1..]);
            let mut list = read_file_list(&fname);
            if list.is_none()
                && let Some(c) = colon
            {
                let (idx, idxs) = parse_idxs(&sub[c + 1..]).ok_or_else(bad_index)?;
                tok.idx = idx;
                tok.idxs = idxs;
                if let Some(fc) = memchr::memrchr(b':', &fname) {
                    fname.truncate(fc);
                }
                list = read_file_list(&fname);
            }
            let Some(list) = list else {
                bail!("Could not read: {}", lossy(&fname));
            };
            tok.nsamples = nsmpl;
            tok.usmpl = vec![false; nsmpl];
            for name in &list {
                let Some(i) = self.header.sample_id(name) else {
                    bail!("No such sample in the VCF: \"{}\"", lossy(name));
                };
                tok.usmpl[i as usize] = true;
            }
            if colon.is_none() {
                tok.set_all_indices();
            }
        } else if let Some(c) = colon {
            if !is_fmt {
                bail!(
                    "Could not parse the index \"{ori}\". (Not a FORMAT tag yet sample index implied.)"
                );
            }
            sample_idxs = Some(parse_idxs(&sub[..c]).ok_or_else(bad_index)?);
            let (idx, idxs) = parse_idxs(&sub[c + 1..]).ok_or_else(bad_index)?;
            tok.idx = idx;
            tok.idxs = idxs;
        } else {
            let (idx1, idxs1) = parse_idxs(sub).ok_or_else(bad_index)?;
            if is_fmt {
                if idxs1 == [-1] {
                    tok.idxs = vec![-1];
                    tok.idx = idx1;
                } else if tok.hdr_id >= 0
                    && !self.is_number_one(HeaderRecordType::Format, tok.hdr_id)
                {
                    bail!(
                        "The FORMAT tag {0} can have multiple subfields, run as {0}[sample:subfield]",
                        lossy(tag)
                    );
                } else {
                    tok.idx = 0;
                }
                sample_idxs = Some((idx1, idxs1));
            } else {
                tok.idxs = idxs1;
                tok.idx = idx1;
            }
        }

        if let Some((idx1, idxs1)) = sample_idxs {
            tok.nsamples = nsmpl;
            tok.usmpl = vec![false; nsmpl];
            if let Ok(i) = usize::try_from(idx1) {
                if i >= nsmpl {
                    bail!("The sample index is too large: {ori}");
                }
                tok.usmpl[i] = true;
            } else {
                let mut i = 0;
                while i < idxs1.len() {
                    if idxs1[i] == 0 {
                        i += 1;
                        continue;
                    }
                    if idxs1[i] == -1 {
                        break;
                    }
                    if i >= nsmpl {
                        bail!("The sample index is too large: {ori}");
                    }
                    tok.usmpl[i] = true;
                    i += 1;
                }
                if idxs1.last() == Some(&-1) {
                    for u in tok.usmpl.iter_mut().skip(i) {
                        *u = true;
                    }
                }
            }
        }

        if tok.idxs.last().is_some_and(|&l| l != -1) {
            tok.nuidxs = tok.idxs.iter().filter(|&&v| v != 0).count();
        }
        Ok(())
    }
}

/// Whether glibc's `strtod` would set `ERANGE` for `s[..end]` parsing as `d`.
fn strtod_range_error(s: &[u8], d: f64, end: usize) -> bool {
    if end == 0 {
        return false;
    }
    let text = &s[..end];
    let mantissa_end = text
        .iter()
        .position(|&c| c == b'e' || c == b'E')
        .unwrap_or(text.len());
    let has_digits = text[..mantissa_end]
        .iter()
        .any(|&c| (b'1'..=b'9').contains(&c));
    let is_inf_text = text
        .iter()
        .find(|c| !matches!(c, b'+' | b'-' | b' ' | b'\t'))
        .is_some_and(|c| c.eq_ignore_ascii_case(&b'i'));
    let is_hex = text.windows(2).any(|w| w.eq_ignore_ascii_case(b"0x"));
    (d.is_infinite() && !is_inf_text) || (!is_hex && has_digits && (d == 0.0 || d.is_subnormal()))
}

/// Compiles an expression (bcftools' `filter_init_`).
pub(super) fn compile(header: &Header, expr: &str, mode: ErrorMode) -> anyhow::Result<Filter> {
    let s = expr.as_bytes();
    let mut b = Builder {
        header,
        expr,
        mode,
        max_unpack: Unpack::STR,
        status: FilterStatus::OK,
        undef_tags: Vec::new(),
        used_tags: Vec::new(),
        ext: Vec::new(),
    };

    let lead = s.iter().take_while(|c| c.is_ascii_whitespace()).count();
    if starts_with_ci(&s[lead..], b"perl:") {
        bail!(
            "\nPerl filtering requires running `configure --enable-perl-filters` at compile time.\n\n"
        );
    }

    let mut ops: Vec<Token> = Vec::new();
    let mut out: Vec<Token> = Vec::new();
    let mut pos = 0;
    let mut last_op: Option<Op> = None;
    while pos < s.len() {
        let lexed = next_token(s, &mut pos);
        let ret = match lexed {
            Lexed::MissingQuotes => bail!("Missing quotes in: {expr}"),
            Lexed::Op(Op::Lft) => {
                ops.push(Token::new(Op::Lft));
                Op::Lft
            }
            Lexed::Op(Op::Rgt) => {
                while ops.last().is_some_and(|t| t.op != Op::Lft) {
                    out.push(ops.pop().expect("an operator is on the stack"));
                }
                if ops.pop().is_none() {
                    bail!("Could not parse: {expr}");
                }
                Op::Rgt
            }
            Lexed::Ext(len) => {
                out.push(b.init_ext(&s[pos..pos + len]));
                pos += len;
                Op::Ext
            }
            Lexed::Pass(fraction) => {
                // The expression in the brackets is evaluated as a normal expression.
                let mut tok = Token::new(Op::Func);
                tok.threshold = -1.0;
                b.max_unpack |= Unpack::FMT;
                tok.func = Some(Func::Pass(fraction));
                tok.tag = Some(match fraction {
                    Fraction::Count => b"N_PASS".to_vec(),
                    Fraction::Fraction => b"F_PASS".to_vec(),
                });
                ops.push(tok);
                continue;
            }
            Lexed::VarArgs(op) => {
                pos = b.var_args(s, pos, op, &mut out)?;
                continue;
            }
            Lexed::Op(mut op) => {
                if op == Op::Sub && last_op != Some(Op::Val) && last_op != Some(Op::Rgt) {
                    // A unary minus: -value becomes -1*(value).
                    out.push(Token::constant(-1.0));
                    op = Op::Mult;
                } else {
                    while ops.last().is_some_and(|t| op.prec() < t.op.prec()) {
                        out.push(ops.pop().expect("an operator is on the stack"));
                    }
                }
                ops.push(Token::new(op));
                op
            }
            Lexed::Val(0) => {
                if pos < s.len() && !s[pos].is_ascii_whitespace() {
                    bail!("Could not parse the expression: [{expr}]");
                }
                break;
            }
            Lexed::Val(len) => {
                let text = &s[pos..pos + len];
                let text = text.strip_suffix(b",").unwrap_or(text);
                out.push(b.init1(text)?);
                pos += len;
                Op::Val
            }
        };
        last_op = Some(ret);
    }
    while let Some(tok) = ops.pop() {
        if tok.op == Op::Lft || tok.op == Op::Rgt {
            bail!("Could not parse the expression: [{expr}]");
        }
        out.push(tok);
    }

    if out.is_empty() {
        if mode == ErrorMode::Exit {
            bail!("[filter.c:4040 filter_init_] Error: could not parse the expression \"{expr}\"");
        }
        b.status |= FilterStatus::OTHER;
    }
    if b.status != FilterStatus::OK {
        return Ok(b.finish(out, 0));
    }

    determine_ext_types(&mut b, &mut out)?;
    postprocess(&b, &mut out)?;

    let nsamples = if b.max_unpack.contains(Unpack::FMT) {
        header.sample_count()
    } else {
        0
    };
    for tok in &mut out {
        if let Some(func) = op_func(tok.op) {
            tok.func = Some(func);
            tok.op = Op::Func;
        }
        expand(&mut tok.values, 1);
        if nsamples > 0 {
            tok.pass_samples = vec![false; nsamples];
        }
    }
    Ok(b.finish(out, nsamples))
}

/// The function an operator token stands for, evaluated as [`Op::Func`].
fn op_func(op: Op) -> Option<Func> {
    Some(match op {
        Op::Max => Func::Max,
        Op::Min => Func::Min,
        Op::Avg => Func::Avg,
        Op::Median => Func::Median,
        Op::Stdev => Func::Stdev,
        Op::Sum => Func::Sum,
        Op::Abs => Func::Abs,
        Op::Cnt => Func::Count,
        Op::Len => Func::Strlen,
        Op::Phred => Func::Phred,
        Op::Binom => Func::Binom,
        Op::Fisher => Func::Fisher,
        Op::PerlSub => Func::Perl,
        Op::SCount => Func::SCount,
        Op::SMax => Func::SMax,
        Op::SMin => Func::SMin,
        Op::SAvg => Func::SAvg,
        Op::SMedian => Func::SMedian,
        Op::SStdev => Func::SStdev,
        Op::SSum => Func::SSum,
        _ => return None,
    })
}

impl Builder<'_> {
    /// Reads the arguments of `binom(...)`, `fisher(...)`, or `perl.NAME(...)`, each a single
    /// value, and appends them and the function to `out`. Returns the position after `)`.
    fn var_args(
        &mut self,
        s: &[u8],
        pos: usize,
        op: Op,
        out: &mut Vec<Token>,
    ) -> anyhow::Result<usize> {
        let mut beg = pos;
        let mut nargs = 0;
        if op == Op::PerlSub {
            while beg < s.len() && (s[beg].is_ascii_alphanumeric() || s[beg] == b'_') {
                beg += 1;
            }
            if s.get(beg) != Some(&b'(') {
                bail!(
                    "[filter.c:3953] Could not parse the expression: {}",
                    self.expr
                );
            }
            let mut name = b"\"".to_vec();
            name.extend_from_slice(&s[pos..beg]);
            name.push(b'"');
            out.push(self.init1(&name)?);
            nargs += 1;
        }
        let Some(len) = memchr::memchr(b')', &s[beg..]) else {
            bail!(
                "[filter.c:3966] Could not parse the expression: {}",
                self.expr
            );
        };
        let end = beg + len;
        let inner = s.get(beg + 1..end).unwrap_or_default();
        for arg in parse_tag_list(inner) {
            nargs += 1;
            out.push(self.init1(arg)?);
        }
        let mut tok = Token::new(op);
        tok.nargs = nargs;
        tok.threshold = -1.0;
        out.push(tok);
        Ok(end + 1)
    }

    fn finish(self, tokens: Vec<Token>, nsamples: usize) -> Filter {
        let ntok = tokens.len();
        // bcftools always decodes ID, REF and ALT, and all of FORMAT; skip that when no token
        // reads the alleles and FORMAT values are looked up one field at a time.
        let eval_unpack = if tokens.iter().all(reads_fixed_or_numeric) {
            self.max_unpack.difference(Unpack::STR | Unpack::FMT)
        } else {
            self.max_unpack
        };
        Filter {
            expr: self.expr.to_owned(),
            tokens,
            stack: Vec::with_capacity(ntok),
            max_unpack: self.max_unpack,
            eval_unpack,
            nsamples,
            status: self.status,
            undef_tags: self.undef_tags,
            used_tags: self.used_tags,
            ext: self.ext,
            ctx: Context::new(nsamples),
        }
    }
}

/// Whether evaluating a token reads only what needs no decoding up front: no function, no ID or
/// FILTER comparison, and a value read from the fixed fields or a numeric FORMAT field (which
/// the record's accessors find without decoding the others).
fn reads_fixed_or_numeric(tok: &Token) -> bool {
    tok.func.is_none()
        && !matches!(tok.comparator, Some(Comparator::Id | Comparator::Filter))
        && tok.setter.is_none_or(|s| match s {
            Setter::Qual | Setter::Pos | Setter::Chrom | Setter::NAlt => true,
            // `[GT]` subscripts read the genotypes.
            Setter::FormatInt | Setter::FormatFloat => tok.idx != -3,
            _ => false,
        })
}

/// Infers the types of external values from the token before them, e.g. `STR_TAG={}` implies
/// a string (bcftools' `determine_ext_types`).
fn determine_ext_types(b: &mut Builder<'_>, out: &mut [Token]) -> anyhow::Result<()> {
    let n = out.len();
    let mut i = 0;
    while i < n {
        if out[i].iext == 0 || out[i].ht_type.is_some() {
            i += 1;
            continue;
        }
        if i == 0 || i + 1 == n || out[i - 1].ht_type.is_none() {
            break;
        }
        out[i].ht_type = out[i - 1].ht_type;
        i += 1;
    }
    if i != n {
        bail!(
            "[filter.c:3791 determine_ext_types] Error: unable to determine the type, use explicit notation: {}",
            b.expr
        );
    }
    for tok in out.iter_mut() {
        let Some(j) = tok.iext.checked_sub(1) else {
            continue;
        };
        let ty = tok.ht_type.map(ExtType::from_value_type);
        b.ext[j] = ty;
        if ty == Some(ExtType::Str) {
            tok.is_str = true;
        }
    }
    Ok(())
}

/// The kind of the token at `i`, if any.
fn op_at(out: &[Token], i: usize) -> Option<Op> {
    out.get(i).map(|t| t.op)
}

fn is_eq_ne(op: Option<Op>) -> bool {
    matches!(op, Some(Op::Eq | Op::Ne))
}

fn is_like(op: Option<Op>) -> bool {
    matches!(op, Some(Op::Like | Op::NLike))
}

/// Resolves what depends on neighbouring tokens: FILTER and TYPE names to IDs and bits, GT
/// classes, `"."` as a missing number, regexes, `@file` lookups, and the logic operators.
fn postprocess(b: &Builder<'_>, out: &mut [Token]) -> anyhow::Result<()> {
    let expr = b.expr;
    let header = b.header;
    let n = out.len();
    let parse_err = || anyhow::anyhow!("Could not parse the expression: {expr}");
    let mut i = 0;
    while i < n {
        if i + 1 < n && matches!(out[i].op, Op::Lt | Op::Bt) && out[i + 1].op == Op::Eq {
            bail!("Error parsing the expression: \"{expr}\"");
        }
        if out[i].hash.is_some() {
            let j = if op_at(out, i + 1) == Some(Op::Val) {
                i + 1
            } else {
                i.checked_sub(1).ok_or_else(parse_err)?
            };
            if out[j].comparator != Some(Comparator::Id) {
                if out[j].comparator.is_some() {
                    bail!(
                        "Error: could not parse the expression with \"@file_name\" syntax (possible todo)"
                    );
                }
                out[j].comparator = Some(Comparator::StringHash);
            }
        }
        match out[i].op {
            Op::Or | Op::OrVec => out[i].func = Some(Func::LogicOr),
            Op::And | Op::AndVec => out[i].func = Some(Func::LogicAnd),
            _ => {}
        }
        if matches!(out[i].op, Op::Eq | Op::Ne) {
            // Look for j="." and k numeric.
            let (mut j, mut k) = (i.checked_sub(1), i.checked_sub(2));
            if !j.is_some_and(|j| out[j].is_str) {
                (j, k) = (i.checked_sub(2), i.checked_sub(1));
            }
            if let Some(j) = j
                && out[j].is_str
                && out[j].key.as_deref() == Some(b".")
            {
                let k =
                    k.ok_or_else(|| anyhow::anyhow!("Error: could not parse the expression"))?;
                let set_missing = if out[k].hdr_id > 0 {
                    let kind = match out[k].hl_type {
                        LineKind::Info => Some(HeaderRecordType::Info),
                        LineKind::Format => Some(HeaderRecordType::Format),
                        LineKind::Filter => Some(HeaderRecordType::Filter),
                        LineKind::None => None,
                    };
                    let ty = kind.and_then(|kind| b.value_type(kind, out[k].hdr_id));
                    matches!(ty, Some(ValueType::Integer | ValueType::Float))
                } else if let Some(tag) = &out[k].tag {
                    tag == b"QUAL"
                } else {
                    bail!("Error: could not parse the expression");
                };
                if set_missing {
                    out[j].is_str = false;
                    out[j].is_missing = true;
                    out[j].values[0] = missing();
                }
            }
        }
        if matches!(out[i].op, Op::Like | Op::NLike) {
            let j = i.checked_sub(1).ok_or_else(|| {
                anyhow::anyhow!(
                    "Could not parse the expression, wrong value for regex operator: {expr}"
                )
            })?;
            let Some(key) = out[j].key.as_mut() else {
                bail!("Could not parse the expression, wrong value for regex operator: {expr}");
            };
            let len = key.len();
            let case = if len > 2 && key.ends_with(b"/i") && key[len - 3] != b'\\' {
                key.truncate(len - 2);
                bre::Case::Insensitive
            } else {
                bre::Case::Sensitive
            };
            match bre::compile(key, case) {
                Ok(re) => out[j].regex = Some(re),
                Err(_) => bail!(
                    "Could not compile the regex expression \"{}\": {expr}",
                    lossy(key)
                ),
            }
        }
        if out[i].is_str
            && out[i].op == Op::Val
            && out[i].key.as_ref().is_some_and(|k| k.contains(&b','))
        {
            let next_is_str = out.get(i + 1).is_some_and(|t| t.is_str);
            let print_note =
                is_eq_ne(op_at(out, i + 1)) || (next_is_str && is_eq_ne(op_at(out, i + 2)));
            if print_note && !COMMA_SEPARATOR_WARNED.swap(true, Ordering::Relaxed) {
                eprint!(
                    "Warning: comma is interpreted as a separator and OR logic is used in string comparisons.\n         (Search the manual for \"Comma in strings\" to learn more.)\n"
                );
            }
        }
        if out[i].op != Op::Val || out[i].tag.is_none() {
            i += 1;
            continue;
        }
        if out[i].setter == Some(Setter::Type) {
            i = resolve_type(out, i, expr)? + 1;
            continue;
        }
        if out[i].tag.as_deref() == Some(b"GT") {
            resolve_gt(out, i, expr)?;
            i += 1;
            continue;
        }
        if out[i].hl_type == LineKind::Filter {
            i = resolve_filter(header, out, i, expr)? + 1;
            continue;
        }
        i += 1;
    }
    Ok(())
}

/// Turns the string compared with TYPE into its variant-type bits. Returns the index of the
/// comparison operator.
fn resolve_type(out: &mut [Token], i: usize, expr: &str) -> anyhow::Result<usize> {
    let parse_err = || anyhow::anyhow!("Could not parse the expression: {expr}");
    if out.get(i + 1).is_none_or(|t| t.key.is_none()) {
        return Err(parse_err());
    }
    let (ival, itok) = if is_eq_ne(op_at(out, i + 1)) || is_like(op_at(out, i + 1)) {
        (i.checked_sub(1).ok_or_else(parse_err)?, i + 1)
    } else if is_eq_ne(op_at(out, i + 2)) || is_like(op_at(out, i + 2)) {
        (i + 1, i + 2)
    } else {
        bail!("[filter.c:4140 filter_init_] Could not parse the expression: {expr}");
    };
    let key = out[ival].key.take().ok_or_else(parse_err)?;
    let bits = match key.to_ascii_lowercase().as_slice() {
        b"snp" | b"snps" => VariantType::SNP.bits() << 1,
        b"indel" | b"indels" => VariantType::INDEL.bits() << 1,
        b"mnp" | b"mnps" => VariantType::MNP.bits() << 1,
        b"other" => VariantType::OTHER.bits() << 1,
        b"bnd" => VariantType::BND.bits() << 1,
        b"overlap" => VariantType::OVERLAP.bits() << 1,
        b"ref" => 1,
        _ => bail!("The type \"{}\" not recognised: {expr}", lossy(&key)),
    };
    let val = &mut out[ival];
    val.threshold = f64::from(bits);
    val.is_str = false;
    val.is_constant = true;
    val.tag = Some(key);
    if is_like(op_at(out, itok)) {
        out[itok].comparator = Some(Comparator::BitAnd);
    }
    Ok(itok)
}

/// Picks the genotype classes a GT token yields from the string it is compared with, and
/// normalises the string (`ar` to `ra`, `HOM` to `hom`).
fn resolve_gt(out: &mut [Token], i: usize, expr: &str) -> anyhow::Result<()> {
    let parse_err = || anyhow::anyhow!("Could not parse the expression: {expr}");
    if i + 1 == out.len() {
        return Err(parse_err());
    }
    let ival = if is_eq_ne(op_at(out, i + 1)) || is_like(op_at(out, i + 1)) {
        // With GT first (`GT=`), bcftools reads the token before the first one, which
        // has no string.
        i.checked_sub(1)
            .ok_or_else(|| anyhow::anyhow!("Comparison between samples is not supported, sorry!"))?
    } else if is_eq_ne(op_at(out, i + 2)) || is_like(op_at(out, i + 2)) {
        i + 1
    } else {
        bail!("[filter.c:4163 filter_init_] Could not parse the expression: {expr}");
    };
    let Some(key) = out[ival].key.as_mut() else {
        bail!("Comparison between samples is not supported, sorry!");
    };
    let lower = key.to_ascii_lowercase();
    let (class, norm): (GtClass, &[u8]) = match lower.as_slice() {
        b"hom" | b"het" | b"hap" => (GtClass::Zygosity, &lower),
        b"mis" | b"ref" | b"alt" => (GtClass::Presence, &lower),
        b"rr" => (GtClass::Alleles, b"rr"),
        b"ra" | b"ar" => (GtClass::Alleles, b"ra"),
        _ if key.as_slice() == b"aA" || key.as_slice() == b"Aa" => (GtClass::Alleles, b"aA"),
        b"aa" => (GtClass::Alleles, b"aa"),
        b"a" => (GtClass::Alleles, b"a"),
        b"r" => (GtClass::Alleles, b"r"),
        _ => return Ok(()),
    };
    *key = norm.to_vec();
    out[i].setter = Some(Setter::Genotype(class));
    Ok(())
}

/// Turns the string compared with FILTER into filter IDs, and `~`/`!~` into subset tests.
/// Returns the index to continue after.
fn resolve_filter(
    header: &Header,
    out: &mut [Token],
    mut i: usize,
    expr: &str,
) -> anyhow::Result<usize> {
    let parse_err = || anyhow::anyhow!("Could not parse the expression: {expr}");
    // With FILTER first (`FILTER=`), bcftools reads the token before the first one, which is
    // never a value.
    let unquoted = || {
        anyhow::anyhow!(
            "[filter.c:4194 filter_init_] Could not parse the expression, an unquoted string value perhaps? {expr}"
        )
    };
    if i + 1 == out.len() {
        return Err(parse_err());
    }
    let itok = i;
    let ival = match (op_at(out, i + 1), op_at(out, i + 2)) {
        (Some(Op::Eq | Op::Ne), _) => i.checked_sub(1).ok_or_else(unquoted)?,
        (Some(Op::Like), _) => {
            out[i + 1].op = Op::In;
            i.checked_sub(1).ok_or_else(unquoted)?
        }
        (Some(Op::NLike), _) => {
            out[i + 1].op = Op::NotIn;
            i.checked_sub(1).ok_or_else(unquoted)?
        }
        (_, Some(Op::Eq | Op::Ne)) => {
            i += 1;
            i
        }
        (_, Some(Op::Like)) => {
            out[i + 2].op = Op::In;
            i += 1;
            i
        }
        (_, Some(Op::NLike)) => {
            out[i + 2].op = Op::NotIn;
            i += 1;
            i
        }
        _ => bail!("[filter.c:4192 filter_init_] Could not parse the expression: {expr}"),
    };
    if out[ival].op != Op::Val || out[ival].key.is_none() {
        return Err(unquoted());
    }
    let tok = &mut out[ival];
    let key = tok.key.take().expect("the key was checked");
    tok.nvalues = 0;
    tok.values.clear();
    let mut has_missing = false;
    for name in key.split(|&c| c == b';') {
        if name == b"." {
            has_missing = true;
            continue;
        }
        let id = header
            .id(name)
            .filter(|&id| header.id_exists(HeaderRecordType::Filter, id));
        let Some(id) = id else {
            bail!(
                "The filter \"{}\" not present in the VCF header",
                lossy(name)
            );
        };
        tok.values.push(f64::from(id));
        tok.nvalues += 1;
    }
    if has_missing && tok.nvalues > 0 {
        bail!(
            "The FILTER expression cannot contain missing value AND filters: \"{}\" ({})",
            lossy(&key),
            tok.nvalues
        );
    }
    tok.tag = Some(key);
    let id = tok.hdr_id;
    out[itok].hdr_id = id;
    Ok(i)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subscripts() {
        assert_eq!(parse_idxs(b""), Some((-2, vec![-1])));
        assert_eq!(parse_idxs(b"*"), Some((-2, vec![-1])));
        assert_eq!(parse_idxs(b"GT"), Some((-3, vec![-1])));
        assert_eq!(parse_idxs(b"3"), Some((3, vec![])));
        assert_eq!(parse_idxs(b"0,2"), Some((-2, vec![1, 0, 1])));
        assert_eq!(parse_idxs(b"1-3"), Some((-2, vec![0, 1, 1, 1])));
        assert_eq!(parse_idxs(b"2-"), Some((-2, vec![0, 0, -1])));
        assert_eq!(parse_idxs(b"x"), None);
    }

    #[test]
    fn tag_lists() {
        assert_eq!(parse_tag_list(b"AD[0,1],DP"), [&b"AD[0,1]"[..], &b"DP"[..]]);
        assert_eq!(parse_tag_list(b""), [&b""[..]]);
    }

    #[test]
    fn paths() {
        assert_eq!(expand_path(b"a/b"), b"a/b");
        assert_eq!(
            expand_path(b"$VTOOLS_NO_SUCH_VARIABLE"),
            b"$VTOOLS_NO_SUCH_VARIABLE"
        );
    }
}
