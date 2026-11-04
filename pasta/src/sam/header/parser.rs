mod context;
mod record;

use std::{collections::HashMap, error, fmt, hash::Hash, mem, str};

use bstr::{BString, ByteSlice};
use indexmap::IndexMap;

pub(crate) use self::context::Context;
use self::record::parse_record;
use super::{
    Header, Programs, ReadGroups, Record, ReferenceSequences,
    layout::{self, Kind},
    record::value::{
        Map,
        map::{self, ReadGroup, header::Version},
    },
};
use crate::variant::log;

/// An error returned when a raw SAM header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A header record is not on the first line.
    UnexpectedHeader,
    /// The record is invalid.
    InvalidRecord(record::ParseError),
    /// A reference sequence name is duplicated.
    DuplicateReferenceSequenceName(BString),
    /// A read group ID is duplicated.
    DuplicateReadGroupId(BString),
    /// A program ID is duplicated.
    DuplicateProgramId(BString),
    /// A comment record is invalid.
    InvalidComment,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidRecord(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedHeader => write!(f, "unexpected header (HD) record"),
            Self::InvalidRecord(_) => f.write_str("invalid record"),
            Self::DuplicateReferenceSequenceName(name) => {
                write!(f, "duplicate reference sequence name: {name}")
            }
            Self::DuplicateReadGroupId(id) => write!(f, "duplicate read group ID: {id}"),
            Self::DuplicateProgramId(id) => write!(f, "duplicate program ID: {id}"),
            Self::InvalidComment => f.write_str("invalid comment record"),
        }
    }
}

/// A SAM header parser.
///
/// By default, the parser follows the SAM specification. [`Parser::lenient`] returns a parser that
/// instead accepts what htslib 1.24 accepts, which is what the SAM, BAM, and CRAM readers use.
#[derive(Default)]
pub struct Parser {
    ctx: Context,
    header: Option<Map<map::Header>>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<BString>,
    binary_reference_sequences: Option<ReferenceSequences>,
    allow_other_record_types: bool,
    layout: layout::Builder,
    program_lines: HashMap<BString, usize>,
}

impl Parser {
    /// Returns a parser that accepts what htslib 1.24 accepts.
    ///
    /// Real files often break the SAM specification in ways that htslib (and so samtools)
    /// accepts. This parser accepts them too, and resolves them as htslib does:
    ///
    ///   * An `@HD` record need not be first. If there is more than one, the first is the
    ///     header's, as htslib looks it up. Writing the header puts it first, as htslib does.
    ///   * The `@HD` `VN` field may be missing or not `<major>.<minor>`. The header's version is
    ///     then the default version.
    ///   * Field values may be empty (e.g., `SM:`), and tags may repeat. The first of repeated
    ///     tags is used, except for `@SQ` `SN`, where the last is used, and `@SQ` `LN`, which may
    ///     only repeat with the same value. A repeated `@PG` `ID` warns.
    ///   * `@SQ` `LN` is read as `strtoll` reads it, but must still be positive.
    ///   * A duplicate `@RG` ID warns, and the first record is the header's, as htslib looks it
    ///     up.
    ///   * A duplicate `@PG` ID is accepted silently, and the last record is the header's, as
    ///     htslib looks it up, including for `PP` links.
    ///
    /// Every line is kept as read (see [`Header`]), including the records that are not the
    /// header's, so writing the header writes them back.
    ///
    /// A duplicate `@SQ` name, an `@SQ` record without a positive `LN`, and a record without a
    /// required name or ID are still errors, as they are in htslib.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{self, header::record::value::map::program::tag};
    ///
    /// let mut parser = sam::header::Parser::lenient();
    /// parser.parse_partial(b"@SQ\tSN:sq0\tLN:8")?;
    /// parser.parse_partial(b"@HD\tSO:coordinate")?;
    /// parser.parse_partial(b"@PG\tID:bwa\tPN:bwa")?;
    /// parser.parse_partial(b"@PG\tID:bwa\tPN:bwa\tVN:2")?;
    /// let header = parser.finish();
    ///
    /// assert!(header.header().is_some());
    /// let programs = header.programs().as_ref();
    /// assert_eq!(programs.len(), 1);
    /// let version = programs[&b"bwa"[..]].other_fields().get(&tag::VERSION);
    /// assert_eq!(version.map(|v| v.as_slice()), Some(&b"2"[..]));
    /// # Ok::<_, sam::header::ParseError>(())
    /// ```
    pub fn lenient() -> Self {
        Self {
            ctx: Context::lenient(),
            ..Default::default()
        }
    }

    /// Makes the given reference sequences, e.g., a BAM file's binary reference sequence
    /// dictionary, authoritative, as htslib does for BAM.
    ///
    /// The header's reference sequences are then these, in this order, and their lengths win over
    /// the lengths in `@SQ` records. An `@SQ` record with the same name adds its other fields,
    /// and its `LN` is not checked. Reference sequences without an `@SQ` record are kept as they
    /// are, and `@SQ` records for other names follow, as htslib appends them.
    pub(crate) fn set_binary_reference_sequences(
        &mut self,
        reference_sequences: ReferenceSequences,
    ) {
        self.binary_reference_sequences = Some(reference_sequences);
    }

    /// Accepts records of types other than `HD`, `SQ`, `RG`, `PG`, and `CO`, which are kept as
    /// lines, as htslib does when it parses a BAM or CRAM header (but not a SAM file's).
    pub(crate) fn allow_other_record_types(&mut self) {
        self.allow_other_record_types = true;
    }

    /// Sets the header text as read, which writing the header writes as is while the header is
    /// unchanged, as htslib does for BAM.
    pub(crate) fn set_text(&mut self, text: BString) {
        self.layout.set_text(text);
    }

    fn is_empty(&self) -> bool {
        self.header.is_none()
            && self.reference_sequences.is_empty()
            && self.read_groups.is_empty()
            && self.programs.as_ref().is_empty()
            && self.comments.is_empty()
    }

    /// Parses and adds a raw record to the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    /// let mut parser = sam::header::Parser::default();
    /// parser.parse_partial(b"@HD\tVN:1.6")?;
    /// # Ok::<_, sam::header::ParseError>(())
    /// ```
    pub fn parse_partial(&mut self, src: &[u8]) -> Result<(), ParseError> {
        if !self.ctx.is_lenient()
            && self.is_empty()
            && let Some(version) = extract_version(src)
        {
            self.ctx = Context::from(version);
        }

        let rewritten = self
            .binary_reference_sequences
            .as_ref()
            .and_then(|reference_sequences| {
                override_reference_sequence_length(src, reference_sequences)
            });
        let src = rewritten.as_deref().unwrap_or(src);

        if self.allow_other_record_types && is_other_record(src) {
            self.layout.push(Kind::Other, src);
            return Ok(());
        }

        let record = parse_record(src, &self.ctx).map_err(ParseError::InvalidRecord)?;

        match record {
            Record::Header(header) => {
                if self.ctx.is_lenient() {
                    // htslib uses the first `@HD` record, wherever it is, and moves it first.
                    if self.header.is_none() {
                        self.header = Some(header);
                        self.layout.push_front(Kind::Header, src);
                    } else {
                        self.layout.push(Kind::Other, src);
                    }
                } else if self.is_empty() {
                    self.header = Some(header);
                    self.layout.push(Kind::Header, src);
                } else {
                    return Err(ParseError::UnexpectedHeader);
                }
            }
            Record::ReferenceSequence(name, reference_sequence) => {
                try_insert(
                    &mut self.reference_sequences,
                    name.clone(),
                    reference_sequence,
                    ParseError::DuplicateReferenceSequenceName,
                )?;

                self.layout.push(Kind::ReferenceSequence(name), src);
            }
            Record::ReadGroup(id, read_group) if self.ctx.is_lenient() => {
                if insert_read_group(&mut self.read_groups, id.clone(), read_group) {
                    self.layout.push(Kind::ReadGroup(id), src);
                } else {
                    self.layout.push(Kind::Other, src);
                }
            }
            Record::ReadGroup(id, read_group) => {
                try_insert(
                    &mut self.read_groups,
                    id.clone(),
                    read_group,
                    ParseError::DuplicateReadGroupId,
                )?;

                self.layout.push(Kind::ReadGroup(id), src);
            }
            Record::Program(id, program) if self.ctx.is_lenient() => {
                // htslib looks up the last record with an ID, so an earlier one is only a line.
                if let Some(i) = self.program_lines.insert(id.clone(), self.layout.len()) {
                    self.layout.set_other(i);
                }

                self.programs.as_mut().insert(id.clone(), program);
                self.layout.push(Kind::Program(id), src);
            }
            Record::Program(id, program) => {
                try_insert(
                    self.programs.as_mut(),
                    id.clone(),
                    program,
                    ParseError::DuplicateProgramId,
                )?;

                self.layout.push(Kind::Program(id), src);
            }
            Record::Comment(comment) => {
                self.layout.push(Kind::Comment(self.comments.len()), src);
                self.comments.push(comment);
            }
        }

        Ok(())
    }

    /// Builds the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    /// let parser = sam::header::Parser::default();
    /// let header = parser.finish();
    /// assert!(header.is_empty());
    /// # Ok::<_, sam::header::ParseError>(())
    /// ```
    pub fn finish(self) -> Header {
        let reference_sequences = match self.binary_reference_sequences {
            Some(binary_reference_sequences) => {
                merge_reference_sequences(binary_reference_sequences, self.reference_sequences)
            }
            None => self.reference_sequences,
        };

        let mut header = Header {
            header: self.header,
            reference_sequences,
            read_groups: self.read_groups,
            programs: self.programs,
            comments: self.comments,
            layout: None,
        };

        let layout = self.layout.build(header.clone());
        header.set_layout(layout);

        header
    }
}

fn extract_version(src: &[u8]) -> Option<Version> {
    use self::record::value::map::header::parse_version;

    const RECORD_PREFIX: &[u8] = b"@HD\t";
    const DELIMITER: u8 = b'\t';
    const FIELD_PREFIX: &[u8] = b"VN:";

    if let Some(raw_value) = src.strip_prefix(RECORD_PREFIX) {
        for raw_field in raw_value.split(|&b| b == DELIMITER) {
            if let Some(s) = raw_field.strip_prefix(FIELD_PREFIX) {
                return parse_version(s).ok();
            }
        }
    }

    None
}

fn try_insert<K, V, F, E>(map: &mut IndexMap<K, V>, key: K, value: V, f: F) -> Result<(), E>
where
    K: Hash + Eq + Clone,
    F: FnOnce(K) -> E,
{
    use indexmap::map::Entry;

    match map.entry(key) {
        Entry::Vacant(e) => {
            e.insert(value);
            Ok(())
        }
        Entry::Occupied(e) => Err(f(e.key().clone())),
    }
}

/// Adds a read group as htslib does: a duplicate ID warns, and the first read group is kept.
///
/// This returns whether the read group was added.
fn insert_read_group(
    read_groups: &mut ReadGroups,
    id: BString,
    read_group: Map<ReadGroup>,
) -> bool {
    use indexmap::map::Entry;

    match read_groups.entry(id) {
        Entry::Vacant(e) => {
            e.insert(read_group);
            true
        }
        Entry::Occupied(e) => {
            log::warn(
                "sam_hrecs_update_hashes",
                format_args!("Duplicate entry \"{}\" in sam header", e.key()),
            );

            false
        }
    }
}

/// Returns whether a line is a record of a type other than `HD`, `SQ`, `RG`, `PG`, and `CO` that
/// htslib's `sam_hrecs_parse_lines` accepts: a two-letter type and `TG:value` fields.
fn is_other_record(src: &[u8]) -> bool {
    let [b'@', a, b, fields @ ..] = src else {
        return false;
    };

    if !a.is_ascii_alphabetic()
        || !b.is_ascii_alphabetic()
        || matches!(&[*a, *b], b"HD" | b"SQ" | b"RG" | b"PG" | b"CO")
    {
        return false;
    }

    let Some(fields) = fields.strip_prefix(b"\t") else {
        return false;
    };

    fields
        .split(|&b| b == b'\t')
        .all(|field| field.len() >= 3 && field[2] == b':')
}

/// Returns an `@SQ` record with the binary length of its reference sequence, if its `LN` (the
/// last, read as `strtoll` reads it) differs.
///
/// htslib's `sam_hrecs_update_hashes` does this for an `@SQ` record of a BAM file's binary
/// reference sequence (with the record's last `SN`): it sets the first `LN` to the binary
/// length, or adds one at the end.
fn override_reference_sequence_length(
    src: &[u8],
    reference_sequences: &ReferenceSequences,
) -> Option<Vec<u8>> {
    const RECORD_PREFIX: &[u8] = b"@SQ";
    const DELIMITER: u8 = b'\t';
    const NAME_PREFIX: &[u8] = b"SN:";
    const LENGTH_PREFIX: &[u8] = b"LN:";

    let fields = src
        .strip_prefix(RECORD_PREFIX)?
        .strip_prefix(&[DELIMITER])?;

    let name = fields
        .split(|&b| b == DELIMITER)
        .filter_map(|field| field.strip_prefix(NAME_PREFIX))
        .next_back()?;

    let length = reference_sequences.get(name.as_bstr())?.length().get();

    let text_length = fields
        .split(|&b| b == DELIMITER)
        .filter_map(|field| field.strip_prefix(LENGTH_PREFIX))
        .next_back()
        .map(strtoll);

    if text_length == Some(length as i128) {
        return None;
    }

    let mut dst = Vec::with_capacity(src.len() + 8);
    dst.extend_from_slice(RECORD_PREFIX);

    let mut is_set = false;

    for field in fields.split(|&b| b == DELIMITER) {
        dst.push(DELIMITER);

        if !is_set && field.starts_with(LENGTH_PREFIX) {
            dst.extend_from_slice(LENGTH_PREFIX);
            dst.extend_from_slice(length.to_string().as_bytes());
            is_set = true;
        } else {
            dst.extend_from_slice(field);
        }
    }

    if !is_set {
        dst.push(DELIMITER);
        dst.extend_from_slice(LENGTH_PREFIX);
        dst.extend_from_slice(length.to_string().as_bytes());
    }

    Some(dst)
}

/// Reads an integer as `strtoll` does, saturating, and 0 if there are no digits.
fn strtoll(src: &[u8]) -> i128 {
    let src = src.trim_ascii_start();

    let (sign, src) = match src.split_first() {
        Some((b'-', rest)) => (-1, rest),
        Some((b'+', rest)) => (1, rest),
        _ => (1, src),
    };

    let n = src
        .iter()
        .take_while(|b| b.is_ascii_digit())
        .fold(0i128, |n, &b| {
            (n * 10 + i128::from(b - b'0')).min(i128::from(i64::MAX) + 1)
        });

    sign * n
}

/// Merges the `@SQ` records into the authoritative (binary) reference sequences.
fn merge_reference_sequences(
    binary_reference_sequences: ReferenceSequences,
    mut text_reference_sequences: ReferenceSequences,
) -> ReferenceSequences {
    let mut reference_sequences = ReferenceSequences::with_capacity(
        binary_reference_sequences
            .len()
            .max(text_reference_sequences.len()),
    );

    for (name, mut reference_sequence) in binary_reference_sequences {
        if let Some(text_reference_sequence) = text_reference_sequences.get_mut(&name) {
            // Its length is already the binary length.
            mem::swap(&mut reference_sequence, text_reference_sequence);
        }

        reference_sequences.insert(name, reference_sequence);
    }

    for (name, reference_sequence) in text_reference_sequences {
        reference_sequences
            .entry(name)
            .or_insert(reference_sequence);
    }

    reference_sequences
}

/// Parses a raw SAM header.
///
/// # Examples
///
/// ```
/// use pasta::sam;
///
/// let s = "\
/// @HD\tVN:1.6\tSO:coordinate
/// @SQ\tSN:sq0\tLN:8
/// @SQ\tSN:sq1\tLN:13
/// ";
///
/// let header: sam::Header = s.parse()?;
///
/// assert!(header.header().is_some());
/// assert_eq!(header.reference_sequences().len(), 2);
/// assert!(header.read_groups().is_empty());
/// assert!(header.programs().as_ref().is_empty());
/// assert!(header.comments().is_empty());
/// # Ok::<(), sam::header::ParseError>(())
/// ```
pub(super) fn parse(s: &str) -> Result<Header, ParseError> {
    let mut parser = Parser::default();

    for line in s.lines() {
        parser.parse_partial(line.as_bytes())?;
    }

    Ok(parser.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZero;

        use crate::sam::header::record::value::map::{
            self, Map, Program, ReadGroup, ReferenceSequence,
            header::{self, Version},
            program,
        };

        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tndls
";

        let actual = parse(s)?;

        let expected = Header::builder()
            .set_header(
                Map::<map::Header>::builder()
                    .set_version(Version::new(1, 6))
                    .insert(header::tag::SORT_ORDER, "coordinate")
                    .build()?,
            )
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::new(const { NonZero::new(13).unwrap() }),
            )
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .add_program(
                "pg0",
                Map::<Program>::builder()
                    .insert(program::tag::NAME, "noodles")
                    .build()?,
            )
            .add_comment("ndls")
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_with_empty_input() -> Result<(), ParseError> {
        let header = parse("")?;

        assert!(header.header().is_none());
        assert!(header.reference_sequences().is_empty());
        assert!(header.read_groups().is_empty());
        assert!(header.programs().as_ref().is_empty());
        assert!(header.comments().is_empty());

        Ok(())
    }

    #[test]
    fn test_parse_without_hd() -> Result<(), ParseError> {
        let header = parse("@SQ\tSN:sq0\tLN:8\n")?;
        assert!(header.header().is_none());
        assert_eq!(header.reference_sequences().len(), 1);
        Ok(())
    }

    #[test]
    fn test_parse_with_multiple_hd() {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@HD\tVN:1.6\tSO:coordinate
";

        assert_eq!(parse(s), Err(ParseError::UnexpectedHeader));
    }

    #[test]
    fn test_parse_with_duplicate_reference_sequence_names() {
        let s = "\
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq0\tLN:8
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateReferenceSequenceName(BString::from(
                "sq0"
            )))
        );
    }

    #[test]
    fn test_parse_with_duplicate_read_group_ids() {
        let s = "\
@RG\tID:rg0
@RG\tID:rg0
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateReadGroupId(BString::from("rg0")))
        );
    }

    #[test]
    fn test_parse_with_duplicate_program_ids() {
        let s = "\
@PG\tID:pg0
@PG\tID:pg0
";
        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateProgramId(BString::from("pg0")))
        );
    }

    #[test]
    fn test_extract_version() {
        assert_eq!(extract_version(b"@HD\tVN:1.6"), Some(Version::new(1, 6)));
        assert_eq!(
            extract_version(b"@HD\tSO:coordinate\tVN:1.6"),
            Some(Version::new(1, 6))
        );
        assert!(extract_version(b"@HD\tVN:NA").is_none());
        assert!(extract_version(b"@SQ\tSN:sq0\tLN:8\tVN:1.6").is_none());
        assert!(extract_version(b"@CO\tVN:1.6").is_none());
    }

    fn parse_lenient(s: &str) -> Result<Header, ParseError> {
        let mut parser = Parser::lenient();

        for line in s.lines() {
            parser.parse_partial(line.as_bytes())?;
        }

        Ok(parser.finish())
    }

    fn program_ids(header: &Header) -> Vec<&[u8]> {
        header
            .programs()
            .as_ref()
            .keys()
            .map(|id| id.as_slice())
            .collect()
    }

    // Fixture `hdr_duppg`.
    #[test]
    fn test_parse_lenient_with_duplicate_program_ids() -> Result<(), Box<dyn std::error::Error>> {
        use crate::sam::header::record::value::map::program::tag;

        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:2000
@PG\tID:bwa\tPN:bwa
@PG\tID:bwa\tPN:bwa\tVN:2
";
        let mut header = parse_lenient(s)?;

        // The header's record is the last with the ID, as htslib looks it up, and both lines are
        // kept.
        assert_eq!(program_ids(&header), [&b"bwa"[..]]);
        assert_eq!(
            header.programs().as_ref()[&b"bwa"[..]]
                .other_fields()
                .get(&tag::VERSION),
            Some(&BString::from("2"))
        );
        assert_eq!(write(&header)?, s);

        // A new program follows the last line, as there are no chains.
        header.add_program_record("ntools", Map::default())?;
        assert_eq!(
            write(&header)?,
            format!("{s}@PG\tID:ntools\tPN:ntools\tPP:bwa\n")
        );

        // A `PP` link names the last record with the ID, as in htslib's `sam_hdr_link_pg`, so
        // this is a cycle through the second `bwa`, and the first `bwa` is the only line no `PP`
        // links to.
        let s = "\
@PG\tID:bwa\tPN:bwa
@PG\tID:ntools\tPN:ntools\tPP:bwa
@PG\tID:bwa\tPN:bwa\tPP:ntools
";
        let mut header = parse_lenient(s)?;
        header.add_program_record("ntools", Map::default())?;
        assert_eq!(
            write(&header)?,
            format!("{s}@PG\tID:ntools.1\tPN:ntools\tPP:bwa\n")
        );

        // Many duplicates add one program.
        let s: String = (0..20_000)
            .map(|i| format!("@PG\tID:bwa\tPN:bwa\tVN:{i}\n"))
            .collect();
        let mut header = parse_lenient(&s)?;
        header.add_program_record("ntools", Map::default())?;
        assert_eq!(program_ids(&header), [&b"bwa"[..], b"ntools"]);
        assert_eq!(write(&header)?.lines().count(), 20_001);

        // The strict parser still rejects duplicates.
        assert_eq!(
            parse(s.as_str()),
            Err(ParseError::DuplicateProgramId(BString::from("bwa")))
        );

        Ok(())
    }

    fn write(header: &Header) -> std::io::Result<String> {
        let mut writer = crate::sam::io::Writer::new(Vec::new());
        writer.write_header(header)?;
        Ok(String::from_utf8(writer.into_inner()).expect("UTF-8"))
    }

    // Fixture `hdr_duprg`.
    #[test]
    fn test_parse_lenient_with_duplicate_read_group_ids() -> Result<(), ParseError> {
        use crate::sam::header::record::value::map::read_group::tag;

        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:2000
@RG\tID:rg1\tSM:s1
@RG\tID:rg1\tSM:s2
";
        let header = parse_lenient(s)?;

        // htslib warns and looks up the first read group.
        assert_eq!(header.read_groups().len(), 1);
        assert_eq!(
            header.read_groups()[&b"rg1"[..]]
                .other_fields()
                .get(&tag::SAMPLE),
            Some(&BString::from("s1"))
        );

        Ok(())
    }

    // Fixture `hdr_novn`, and `VN:abc`.
    #[test]
    fn test_parse_lenient_with_missing_or_invalid_version() -> Result<(), ParseError> {
        use crate::sam::header::record::value::map::header::tag;

        for s in [
            "@HD\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n",
            "@HD\tVN:abc\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n",
        ] {
            let header = parse_lenient(s)?;
            let hd = header.header().expect("@HD is kept");
            assert_eq!(hd.version(), Version::default());
            assert_eq!(
                hd.other_fields().get(&tag::SORT_ORDER),
                Some(&BString::from("coordinate"))
            );
            assert_eq!(header.reference_sequences().len(), 1);

            assert!(matches!(parse(s), Err(ParseError::InvalidRecord(_))));
        }

        Ok(())
    }

    #[test]
    fn test_parse_lenient_with_late_or_repeated_hd() -> Result<(), ParseError> {
        use crate::sam::header::record::value::map::header::tag;

        let s = "@SQ\tSN:chr1\tLN:2000\n@HD\tVN:1.6\tSO:coordinate\n";
        let header = parse_lenient(s)?;
        assert_eq!(
            header.header().map(|hd| hd.version()),
            Some(Version::new(1, 6))
        );
        assert_eq!(parse(s), Err(ParseError::UnexpectedHeader));

        // htslib uses the first `@HD` record.
        let s = "@HD\tVN:1.6\tSO:coordinate\n@HD\tVN:1.5\tSO:unsorted\n";
        let header = parse_lenient(s)?;
        let hd = header.header().expect("@HD is kept");
        assert_eq!(hd.version(), Version::new(1, 6));
        assert_eq!(
            hd.other_fields().get(&tag::SORT_ORDER),
            Some(&BString::from("coordinate"))
        );

        Ok(())
    }

    #[test]
    fn test_parse_lenient_with_empty_values() -> Result<(), ParseError> {
        use crate::sam::header::record::value::map::read_group::tag;

        let s = "@HD\tVN:1.6\n@RG\tID:rg1\tSM:\n";
        let header = parse_lenient(s)?;
        assert_eq!(
            header.read_groups()[&b"rg1"[..]]
                .other_fields()
                .get(&tag::SAMPLE),
            Some(&BString::default())
        );
        assert!(matches!(parse(s), Err(ParseError::InvalidRecord(_))));

        Ok(())
    }

    #[test]
    fn test_parse_lenient_keeps_lines_as_read() -> Result<(), Box<dyn std::error::Error>> {
        use crate::sam::header::record::value::map::read_group::tag;

        // Duplicate `@RG` lines, a second `@HD`, repeated tags, empty values, and values that
        // are not written as they are read (`LN:2000kb`) are written back as read.
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:2000kb\tAS:a\tAS:b
@CO\tbetween
@HD\tVN:1.4
@RG\tID:rg1\tSM:s1\tSM:s2\tLB:
@RG\tID:rg1\tSM:s3
@PG\tID:pg0
";
        let mut header = parse_lenient(s)?;
        assert_eq!(write(&header)?, s);

        // A changed record keeps the fields that did not change, and new records follow the last
        // line of their type.
        header
            .read_groups_mut()
            .get_mut(&b"rg1"[..])
            .expect("rg1")
            .other_fields_mut()
            .insert(tag::SAMPLE, "x".into());
        header
            .read_groups_mut()
            .insert("rg2".into(), Map::default());
        header.add_comment("new");
        header.add_program_record("ntools", Map::default())?;

        assert_eq!(
            write(&header)?,
            "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:2000kb\tAS:a\tAS:b
@CO\tbetween
@CO\tnew
@HD\tVN:1.4
@RG\tID:rg1\tSM:x\tSM:s2\tLB:
@RG\tID:rg1\tSM:s3
@RG\tID:rg2
@PG\tID:pg0
@PG\tID:ntools\tPN:ntools\tPP:pg0
"
        );

        // htslib moves the first `@HD` first.
        let header = parse_lenient("@SQ\tSN:chr1\tLN:8\n@HD\tVN:1.6\n@HD\tVN:1.5\n")?;
        assert_eq!(
            write(&header)?,
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n@HD\tVN:1.5\n"
        );

        // Other record types are only accepted from BAM and CRAM, as in htslib.
        assert!(matches!(
            parse_lenient("@XY\tAB:1\n"),
            Err(ParseError::InvalidRecord(_))
        ));

        let mut parser = Parser::lenient();
        parser.allow_other_record_types();
        parser.parse_partial(b"@HD\tVN:1.6")?;
        parser.parse_partial(b"@XY\tAB:1")?;
        assert!(parser.parse_partial(b"@XY\tAB").is_err());
        let header = parser.finish();
        assert_eq!(write(&header)?, "@HD\tVN:1.6\n@XY\tAB:1\n");

        Ok(())
    }

    #[test]
    fn test_parse_lenient_with_errors_htslib_keeps() {
        // htslib rejects these too, when it parses the text (e.g., for a SAM file).
        assert!(matches!(
            parse_lenient("@SQ\tSN:chr1\tLN:0\n"),
            Err(ParseError::InvalidRecord(_))
        ));
        assert!(matches!(
            parse_lenient("@SQ\tSN:chr1\n"),
            Err(ParseError::InvalidRecord(_))
        ));
        assert_eq!(
            parse_lenient("@SQ\tSN:chr1\tLN:8\n@SQ\tSN:chr1\tLN:8\n"),
            Err(ParseError::DuplicateReferenceSequenceName(BString::from(
                "chr1"
            )))
        );
        assert!(matches!(
            parse_lenient("@RG\tSM:s1\n"),
            Err(ParseError::InvalidRecord(_))
        ));
        assert!(matches!(
            parse_lenient("@PG\tPN:x\n"),
            Err(ParseError::InvalidRecord(_))
        ));
        assert!(matches!(
            parse_lenient("@XY\tA:b\n"),
            Err(ParseError::InvalidRecord(_))
        ));
    }

    fn reference_sequences(entries: &[(&str, usize)]) -> ReferenceSequences {
        use std::num::NonZero;

        use crate::sam::header::record::value::map::ReferenceSequence;

        entries
            .iter()
            .map(|&(name, len)| {
                (
                    BString::from(name),
                    Map::<ReferenceSequence>::new(NonZero::new(len).unwrap()),
                )
            })
            .collect()
    }

    fn parse_with_binary_reference_sequences(
        s: &str,
        binary: &[(&str, usize)],
    ) -> Result<Header, ParseError> {
        let mut parser = Parser::lenient();
        parser.set_binary_reference_sequences(reference_sequences(binary));

        for line in s.lines() {
            parser.parse_partial(line.as_bytes())?;
        }

        Ok(parser.finish())
    }

    fn names_and_lengths(header: &Header) -> Vec<(&[u8], usize)> {
        header
            .reference_sequences()
            .iter()
            .map(|(name, map)| (name.as_slice(), usize::from(map.length())))
            .collect()
    }

    // Fixtures `hdr_lnmismatch` and `hdr_sqmissing`, and `LN:0`.
    #[test]
    fn test_parse_with_binary_reference_sequences() -> Result<(), ParseError> {
        use crate::sam::header::record::value::map::reference_sequence::tag;

        let binary = [("chr1", 2000), ("chr2", 1500)];

        // The binary length wins over a different, zero, missing, or invalid text `LN`.
        for sq2 in [
            "@SQ\tSN:chr2\tLN:1501\tAS:x",
            "@SQ\tSN:chr2\tLN:0\tAS:x",
            "@SQ\tSN:chr2\tAS:x",
        ] {
            let s = format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\n{sq2}\n");
            let header = parse_with_binary_reference_sequences(&s, &binary)?;
            assert_eq!(
                names_and_lengths(&header),
                [(&b"chr1"[..], 2000), (b"chr2", 1500)]
            );

            // The record's other fields are kept.
            assert_eq!(
                header.reference_sequences()[&b"chr2"[..]]
                    .other_fields()
                    .get(&tag::ASSEMBLY_ID),
                Some(&BString::from("x"))
            );
        }

        // A reference sequence missing from the text is kept.
        let header =
            parse_with_binary_reference_sequences("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\n", &binary)?;
        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500)]
        );

        // The binary order wins, and `@SQ` records for other names are appended, as htslib does.
        let s = "@SQ\tSN:chr3\tLN:10\n@SQ\tSN:chr2\tLN:1500\n@SQ\tSN:chr1\tLN:2000\n";
        let header = parse_with_binary_reference_sequences(s, &binary)?;
        assert_eq!(
            names_and_lengths(&header),
            [(&b"chr1"[..], 2000), (b"chr2", 1500), (b"chr3", 10)]
        );

        // As in htslib, `LN` may not repeat with different values.
        assert!(
            parse_with_binary_reference_sequences("@SQ\tSN:chr2\tLN:abc\tAS:x\tLN:7\n", &binary)
                .is_err()
        );

        // An `@SQ` record for another name still needs a positive length.
        assert!(parse_with_binary_reference_sequences("@SQ\tSN:chr3\tLN:0\n", &binary).is_err());

        Ok(())
    }

    #[test]
    fn test_override_reference_sequence_length() {
        let reference_sequences = reference_sequences(&[("sq0", 8)]);

        assert_eq!(
            override_reference_sequence_length(b"@SQ\tSN:sq0\tLN:0\tAS:a", &reference_sequences),
            Some(b"@SQ\tSN:sq0\tLN:8\tAS:a".to_vec())
        );
        assert!(
            override_reference_sequence_length(b"@SQ\tSN:sq0\tLN:8kb", &reference_sequences)
                .is_none()
        );
        assert_eq!(
            override_reference_sequence_length(b"@SQ\tSN:sq0", &reference_sequences),
            Some(b"@SQ\tSN:sq0\tLN:8".to_vec())
        );
        assert!(
            override_reference_sequence_length(b"@SQ\tSN:sq1\tLN:0", &reference_sequences)
                .is_none()
        );
        assert!(
            override_reference_sequence_length(b"@RG\tID:sq0\tLN:0", &reference_sequences)
                .is_none()
        );
    }
}
