//! Header lines (HTSlib's `bcf_hrec_t`) and their parser (`bcf_hdr_parse_line`).

use crate::variant::log;

/// The kind of a header line (HTSlib's `BCF_HL_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[repr(u8)]
pub enum HeaderRecordType {
    /// `##FILTER=<...>` (`BCF_HL_FLT`).
    Filter = 0,
    /// `##INFO=<...>` (`BCF_HL_INFO`).
    Info = 1,
    /// `##FORMAT=<...>` (`BCF_HL_FMT`).
    Format = 2,
    /// `##contig=<...>` (`BCF_HL_CTG`).
    Contig = 3,
    /// Any other structured line, `##KEY=<A=..,B=..>` (`BCF_HL_STR`).
    Structured = 4,
    /// A generic line, `##KEY=VALUE` (`BCF_HL_GEN`).
    Generic = 5,
}

impl HeaderRecordType {
    /// The index of FILTER/INFO/FORMAT in the shared ID dictionary's per-type slots, or `None`
    /// for other kinds.
    #[inline]
    #[must_use]
    pub fn id_slot(self) -> Option<usize> {
        match self {
            Self::Filter => Some(0),
            Self::Info => Some(1),
            Self::Format => Some(2),
            _ => None,
        }
    }
}

/// One header line: a generic `##key=value` line or a structured `##key=<k1=v1,k2=v2,...>`
/// line (HTSlib's `bcf_hrec_t`).
///
/// Structured values keep their surrounding quotes, as HTSlib stores them, so formatting the
/// record reproduces the input line.
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct HeaderRecord {
    pub(crate) kind: HeaderRecordType,
    pub(crate) key: Vec<u8>,
    pub(crate) value: Option<Vec<u8>>,
    pub(crate) keys: Vec<Vec<u8>>,
    pub(crate) vals: Vec<Vec<u8>>,
}

impl HeaderRecord {
    /// Creates a generic `##key=value` line.
    pub fn generic(key: impl Into<Vec<u8>>, value: impl Into<Vec<u8>>) -> Self {
        let mut rec = Self {
            kind: HeaderRecordType::Generic,
            key: key.into(),
            value: Some(value.into()),
            keys: Vec::new(),
            vals: Vec::new(),
        };
        rec.set_type();
        rec
    }

    /// Creates a structured line `##key=<>` with no fields yet; add them with
    /// [`Self::push_field`].
    pub fn structured(key: impl Into<Vec<u8>>) -> Self {
        Self {
            kind: HeaderRecordType::Structured,
            key: key.into(),
            value: None,
            keys: Vec::new(),
            vals: Vec::new(),
        }
    }

    /// Returns the kind of line.
    #[inline]
    #[must_use]
    pub fn kind(&self) -> HeaderRecordType {
        self.kind
    }

    /// Returns the key: `INFO`, `contig`, `fileformat`, ...
    #[inline]
    #[must_use]
    pub fn key(&self) -> &[u8] {
        &self.key
    }

    /// Returns the value of a generic line.
    #[inline]
    #[must_use]
    pub fn value(&self) -> Option<&[u8]> {
        self.value.as_deref()
    }

    /// Returns the number of structured fields.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Returns whether the line has no structured fields.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// Iterates over the structured `(key, value)` pairs, values with their quotes.
    pub fn fields(&self) -> impl Iterator<Item = (&[u8], &[u8])> {
        self.keys
            .iter()
            .zip(&self.vals)
            .map(|(k, v)| (k.as_slice(), v.as_slice()))
    }

    /// Returns the index of a structured field, compared case-insensitively
    /// (HTSlib's `bcf_hrec_find_key`).
    #[must_use]
    pub fn find_key(&self, key: &[u8]) -> Option<usize> {
        self.keys.iter().position(|k| k.eq_ignore_ascii_case(key))
    }

    /// Returns the value of a structured field (compared case-insensitively), with quotes.
    #[must_use]
    pub fn get(&self, key: &[u8]) -> Option<&[u8]> {
        self.find_key(key).map(|i| self.vals[i].as_slice())
    }

    /// Returns the value of the field at `i`.
    #[must_use]
    pub fn field_value(&self, i: usize) -> Option<&[u8]> {
        self.vals.get(i).map(Vec::as_slice)
    }

    /// Appends a structured field (HTSlib's `bcf_hrec_add_key` + `bcf_hrec_set_val`). A
    /// quoted value gets double quotes around it.
    pub fn push_field(&mut self, key: impl Into<Vec<u8>>, value: &[u8], quoting: Quoting) {
        self.keys.push(key.into());
        self.vals.push(quoted(value, quoting));
        if self.value.is_none() {
            self.set_type();
        }
    }

    /// Replaces the value of field `i` (HTSlib's `bcf_hrec_set_val`).
    ///
    /// # Panics
    ///
    /// Panics if `i` is not a field index.
    pub fn set_field_value(&mut self, i: usize, value: &[u8], quoting: Quoting) {
        self.vals[i] = quoted(value, quoting);
    }

    /// Removes the field at `i`.
    pub fn remove_field(&mut self, i: usize) {
        if i < self.keys.len() {
            self.keys.remove(i);
            self.vals.remove(i);
        }
    }

    /// Sets the value of a generic line.
    pub fn set_value(&mut self, value: impl Into<Vec<u8>>) {
        self.value = Some(value.into());
    }

    /// Returns a copy without the `IDX` field (HTSlib's `bcf_hrec_dup`).
    #[must_use]
    pub fn dup(&self) -> Self {
        let mut out = self.clone();
        let mut i = 0;
        while i < out.keys.len() {
            if out.keys[i] == b"IDX" {
                out.keys.remove(i);
                out.vals.remove(i);
            } else {
                i += 1;
            }
        }
        out
    }

    /// Appends the line as text, with its newline (HTSlib's `_bcf_hrec_format`). VCF text omits
    /// `IDX` fields; BCF header text keeps them.
    pub fn format_into(&self, out: &mut Vec<u8>, flavor: HeaderTextFlavor) {
        out.extend_from_slice(b"##");
        out.extend_from_slice(&self.key);
        out.push(b'=');
        if let Some(value) = &self.value {
            out.extend_from_slice(value);
        } else {
            out.push(b'<');
            let mut nout = 0;
            for (k, v) in self.keys.iter().zip(&self.vals) {
                if flavor == HeaderTextFlavor::Vcf && k == b"IDX" {
                    continue;
                }
                if nout > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(k);
                out.push(b'=');
                out.extend_from_slice(v);
                nout += 1;
            }
            out.push(b'>');
        }
        out.push(b'\n');
    }

    /// Returns the line as text (VCF flavor, without `IDX`), with its newline.
    #[must_use]
    pub fn to_vcf_line(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.format_into(&mut out, HeaderTextFlavor::Vcf);
        out
    }

    /// Sets the kind from the key and fields (HTSlib's `bcf_hrec_set_type`).
    pub(crate) fn set_type(&mut self) {
        self.kind = match self.key.as_slice() {
            b"contig" => HeaderRecordType::Contig,
            b"INFO" => HeaderRecordType::Info,
            b"FILTER" => HeaderRecordType::Filter,
            b"FORMAT" => HeaderRecordType::Format,
            _ if !self.keys.is_empty() => HeaderRecordType::Structured,
            _ => HeaderRecordType::Generic,
        };
    }

    /// Returns the value of the `ID` field, compared case-insensitively.
    #[must_use]
    pub fn id(&self) -> Option<&[u8]> {
        self.get(b"ID")
    }
}

/// Whether a structured value is written with surrounding double quotes.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Quoting {
    /// Written as given.
    Bare,
    /// Wrapped in double quotes.
    Quoted,
}

/// Which header text to produce: VCF (no `IDX` fields) or BCF (with them).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum HeaderTextFlavor {
    /// VCF header text, without `IDX` fields.
    Vcf,
    /// BCF header text, keeping `IDX` fields.
    Bcf,
}

fn quoted(value: &[u8], quoting: Quoting) -> Vec<u8> {
    match quoting {
        Quoting::Bare => value.to_vec(),
        Quoting::Quoted => {
            let mut v = Vec::with_capacity(value.len() + 2);
            v.push(b'"');
            v.extend_from_slice(value);
            v.push(b'"');
            v
        }
    }
}

/// The result of parsing one header line.
#[derive(Debug)]
pub(crate) enum ParsedLine {
    /// A parsed line, and the number of bytes consumed (including the newline).
    Record(HeaderRecord, usize),
    /// Not a `##` line; nothing consumed.
    NotMeta,
    /// A malformed `##` line (already reported), and the bytes to skip.
    Malformed(usize),
}

fn is_escaped(s: &[u8], min: usize, pos: usize) -> bool {
    let mut n = 0;
    let mut p = pos;
    while p > min && s[p - 1] == b'\\' {
        n += 1;
        p -= 1;
    }
    n % 2 == 1
}

#[inline]
fn c_at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// Parses one header line from the start of `text` (HTSlib's `bcf_hdr_parse_line`). The line
/// ends at a newline or at the end of `text` (or at a NUL byte, as in C).
pub(crate) fn parse_line(text: &[u8]) -> ParsedLine {
    let line_end = |from: usize| {
        let mut q = from;
        while c_at(text, q) != 0 && c_at(text, q) != b'\n' {
            q += 1;
        }
        q
    };
    let consumed = |q: usize| q + usize::from(c_at(text, q) != 0);

    if c_at(text, 0) != b'#' || c_at(text, 1) != b'#' {
        return ParsedLine::NotMeta;
    }
    let p = 2;
    let mut q = p;
    while c_at(text, q) != 0 && c_at(text, q) != b'=' && c_at(text, q) != b'\n' {
        q += 1;
    }
    let n = q - p;
    if c_at(text, q) != b'=' || n == 0 {
        return malformed(text, line_end(q));
    }
    let mut hrec = HeaderRecord {
        kind: HeaderRecordType::Generic,
        key: text[p..q].to_vec(),
        value: None,
        keys: Vec::new(),
        vals: Vec::new(),
    };
    q += 1;
    let mut p = q;
    if c_at(text, p) != b'<' {
        // Generic field, e.g. ##samtoolsVersion=0.1.18-r579
        let q = line_end(q);
        hrec.value = Some(text[p..q].to_vec());
        hrec.set_type();
        return ParsedLine::Record(hrec, consumed(q));
    }

    // Structured line, e.g. ##INFO=<ID=PV1,Number=1,Type=Float,Description="P-value">
    let mut nopen = 1;
    while c_at(text, q) != 0 && c_at(text, q) != b'\n' && nopen > 0 {
        q += 1;
        p = q;
        while c_at(text, q) == b' ' {
            p += 1;
            q += 1;
        }
        // ^[A-Za-z_][0-9A-Za-z_.]*$
        let c = c_at(text, q);
        if p == q && c != 0 && (c.is_ascii_alphabetic() || c == b'_') {
            q += 1;
            while {
                let c = c_at(text, q);
                c != 0 && (c.is_ascii_alphanumeric() || c == b'_' || c == b'.')
            } {
                q += 1;
            }
        }
        let n = q - p;
        let mut m = 0;
        while c_at(text, q) == b' ' {
            q += 1;
            m += 1;
        }
        if c_at(text, q) != b'=' || n == 0 {
            return malformed(text, line_end(q));
        }
        hrec.keys.push(text[p..q - m].to_vec());
        q += 1;
        p = q;
        while c_at(text, q) == b' ' {
            p += 1;
            q += 1;
        }
        let mut quoted = false;
        let mut ending = 0u8;
        match c_at(text, p) {
            b'"' => {
                quoted = true;
                ending = b'"';
                p += 1;
            }
            b'[' => {
                quoted = true;
                ending = b']';
            }
            _ => {}
        }
        if quoted {
            q += 1;
        }
        while c_at(text, q) != 0 && c_at(text, q) != b'\n' {
            let c = c_at(text, q);
            if quoted {
                if c == ending && !is_escaped(text, p, q) {
                    break;
                }
            } else {
                if c == b'<' {
                    nopen += 1;
                }
                if c == b'>' {
                    nopen -= 1;
                }
                if nopen == 0 {
                    break;
                }
                if c == b',' && nopen == 1 {
                    break;
                }
            }
            q += 1;
        }
        let mut r = q;
        if quoted && ending == b']' {
            if c_at(text, q) == ending {
                r += 1;
                q += 1;
                quoted = false;
            } else {
                log::error(
                    "bcf_hdr_parse_line",
                    format_args!(
                        "Missing ']' in header line {}",
                        log::strprint(&text[..q.min(text.len())], b'"', 320)
                    ),
                );
                return ParsedLine::Malformed(consumed(line_end(q)));
            }
        }
        while r > p && text[r - 1] == b' ' {
            r -= 1;
        }
        let value = text[p..r.max(p)].to_vec();
        hrec.vals.push(if quoted {
            let mut v = Vec::with_capacity(value.len() + 2);
            v.push(b'"');
            v.extend_from_slice(&value);
            v.push(b'"');
            v
        } else {
            value
        });
        if quoted && c_at(text, q) == ending {
            q += 1;
        }
        if c_at(text, q) == b'>' {
            if nopen > 0 {
                nopen -= 1;
            }
            q += 1;
        }
    }
    if nopen > 0 {
        // HTSlib prints the rest of the header text from this line on, up to its NUL.
        let end = text.iter().position(|&c| c == 0).unwrap_or(text.len());
        log::warn(
            "bcf_hdr_parse_line",
            format_args!(
                "Incomplete header line, trying to proceed anyway:\n\t[{}]\n\t[{}]",
                String::from_utf8_lossy(&text[..end]),
                i32::from(c_at(text, q) as i8)
            ),
        );
    }

    // Skip to end of line
    let mut nonspace = false;
    while c_at(text, q) != 0 && c_at(text, q) != b'\n' {
        nonspace |= !crate::variant::number::is_space(c_at(text, q));
        q += 1;
    }
    if nonspace {
        log::warn(
            "bcf_hdr_parse_line",
            format_args!(
                "Dropped trailing junk from header line '{}'",
                log::strprint(&text[..q], b'"', 320)
            ),
        );
    }
    hrec.set_type();
    ParsedLine::Record(hrec, consumed(q))
}

fn malformed(text: &[u8], q: usize) -> ParsedLine {
    log::error(
        "bcf_hdr_parse_line",
        format_args!(
            "Could not parse the header line: {}",
            log::strprint(&text[..q], b'"', 320)
        ),
    );
    ParsedLine::Malformed(q + usize::from(c_at(text, q) != 0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> HeaderRecord {
        match parse_line(s.as_bytes()) {
            ParsedLine::Record(r, _) => r,
            other => panic!("{s}: {other:?}"),
        }
    }

    #[test]
    fn structured_round_trip() {
        let line =
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth, \\\"raw\\\"\">\n";
        let r = parse(line);
        assert_eq!(r.kind(), HeaderRecordType::Info);
        assert_eq!(r.get(b"id"), Some(&b"DP"[..]));
        assert_eq!(r.to_vcf_line(), line.as_bytes());
    }

    #[test]
    fn generic_line() {
        let r = parse("##fileformat=VCFv4.2\n##x=1");
        assert_eq!(r.kind(), HeaderRecordType::Generic);
        assert_eq!(r.value(), Some(&b"VCFv4.2"[..]));
    }

    #[test]
    fn spaces_are_trimmed() {
        let r = parse("##FILTER=< ID = q10 , Description=\"x\">");
        assert_eq!(r.to_vcf_line(), b"##FILTER=<ID=q10,Description=\"x\">\n");
        // HTSlib rejects a space before the closing bracket.
        assert!(matches!(
            parse_line(b"##FILTER=<ID=q10,Description=\"x\" >"),
            ParsedLine::Malformed(_)
        ));
    }
}
