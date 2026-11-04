//! The lines of a SAM header as read.
//!
//! htslib keeps every header line, with its tags as written and in its place, and only resolves
//! lookups by name or ID to one of them (`sam_hrecs_t`). A [`Header`] groups its records by
//! type and keeps one record per name or ID, so a header that was read also keeps a [`Layout`]:
//! its lines, in order, as read. Writing the header writes these lines, rewriting only the
//! records that have changed since, and adds new records after the last line of their type, as
//! htslib does.

use std::ops::Range;

use bstr::{BStr, BString, ByteSlice};

use super::Header;

/// The lines of a header as read.
#[derive(Clone, Debug, Default)]
pub(crate) struct Layout {
    lines: Vec<Line>,
    buf: Vec<u8>,
    snapshot: Header,
    text: Option<BString>,
}

/// A header line as read.
#[derive(Clone, Debug)]
pub(crate) struct Line {
    kind: Kind,
    src: Range<usize>,
}

/// What a header line is.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum Kind {
    /// The `@HD` record.
    Header,
    /// The `@SQ` record with this name.
    ReferenceSequence(BString),
    /// The `@RG` record with this ID.
    ReadGroup(BString),
    /// The `@PG` record with this ID.
    Program(BString),
    /// The comment at this index.
    Comment(usize),
    /// A line that is not one of the header's records, written as read: a later `@HD`, a later
    /// `@RG` or an earlier `@PG` with the same ID as another, or a line of another type.
    Other,
}

impl Layout {
    /// Returns the lines.
    pub(crate) fn lines(&self) -> &[Line] {
        &self.lines
    }

    /// Returns a line as read, without its line ending.
    pub(crate) fn line(&self, line: &Line) -> &BStr {
        self.buf[line.src.clone()].as_bstr()
    }

    /// Returns the header as read, which tells which records have changed since.
    pub(crate) fn snapshot(&self) -> &Header {
        &self.snapshot
    }

    /// Returns the header text as read, if it is to be written as is while the header is
    /// unchanged.
    ///
    /// htslib writes a BAM file's header text as read until something parses it, which adding
    /// an `@PG` line or changing the header does. This is that text, up to the first NUL.
    pub(crate) fn text(&self) -> Option<&BStr> {
        self.text.as_ref().map(|text| text.as_bstr())
    }
}

impl Layout {
    /// Adds a line for a new record at the end, and removes its record's earlier line, if any.
    pub(crate) fn push_record(&mut self, kind: Kind, src: &[u8]) {
        self.lines.retain(|line| line.kind != kind);
        self.text = None;

        let start = self.buf.len();
        self.buf.extend_from_slice(src);

        self.lines.push(Line {
            kind,
            src: start..self.buf.len(),
        });
    }
}

impl Layout {
    /// Removes the [`Kind::Other`] lines of a record type.
    pub(crate) fn remove_other_lines(&mut self, record_type: [u8; 2]) {
        let buf = &self.buf;
        let len = self.lines.len();

        self.lines.retain(|line| {
            line.kind != Kind::Other || self::record_type(&buf[line.src.clone()]) != record_type
        });

        if self.lines.len() != len {
            // The header has changed, so it is no longer written as read.
            self.text = None;
        }
    }
}

impl Line {
    pub(crate) fn kind(&self) -> &Kind {
        &self.kind
    }
}

/// Builds a [`Layout`] as a header is parsed.
#[derive(Debug, Default)]
pub(crate) struct Builder {
    lines: Vec<Line>,
    buf: Vec<u8>,
    text: Option<BString>,
}

impl Builder {
    /// Adds a line.
    pub(crate) fn push(&mut self, kind: Kind, src: &[u8]) {
        let line = self.line(kind, src);
        self.lines.push(line);
    }

    /// Adds a line before every other line, where htslib puts the first `@HD` line.
    pub(crate) fn push_front(&mut self, kind: Kind, src: &[u8]) {
        let line = self.line(kind, src);
        self.lines.insert(0, line);
    }

    /// Makes the line at index `i`, which a later line replaces as the header's record, an
    /// [`Kind::Other`] line.
    pub(crate) fn set_other(&mut self, i: usize) {
        self.lines[i].kind = Kind::Other;
    }

    pub(crate) fn len(&self) -> usize {
        self.lines.len()
    }

    /// Sets the header text to write as is while the header is unchanged.
    pub(crate) fn set_text(&mut self, text: BString) {
        self.text = Some(text);
    }

    /// Builds the layout of the given header, as read.
    pub(crate) fn build(self, snapshot: Header) -> Option<Layout> {
        if self.lines.is_empty() && self.text.is_none() {
            return None;
        }

        Some(Layout {
            lines: self.lines,
            buf: self.buf,
            snapshot,
            text: self.text,
        })
    }

    fn line(&mut self, kind: Kind, src: &[u8]) -> Line {
        let start = self.buf.len();
        self.buf.extend_from_slice(src);

        Line {
            kind,
            src: start..self.buf.len(),
        }
    }
}

/// Returns whether an [`Kind::Other`] line is a duplicate of a record that the header no longer
/// has: a later `@HD` without an `@HD`, or an `@RG` or `@PG` line whose ID is not a read group
/// or program. Such a line is left out, as htslib's `sam_hdr_remove_except` removes every line
/// with an ID.
pub(crate) fn is_removed_duplicate(header: &Header, line: &[u8]) -> bool {
    match record_type(line) {
        HEADER => header.header().is_none(),
        READ_GROUP => {
            field(line, ID).is_some_and(|id| !header.read_groups().contains_key(id.as_bstr()))
        }
        PROGRAM => {
            field(line, ID).is_some_and(|id| !header.programs().as_ref().contains_key(id.as_bstr()))
        }
        _ => false,
    }
}

const HEADER: [u8; 2] = *b"HD";
const READ_GROUP: [u8; 2] = *b"RG";
const PROGRAM: [u8; 2] = *b"PG";
const ID: [u8; 2] = *b"ID";

/// Returns a header line's record type, e.g., `SQ`.
pub(crate) fn record_type(line: &[u8]) -> [u8; 2] {
    match line {
        [b'@', a, b, ..] => [*a, *b],
        _ => [0; 2],
    }
}

/// Returns a header line's fields, as `(tag, value)` pairs.
pub(crate) fn fields(line: &[u8]) -> impl Iterator<Item = ([u8; 2], &[u8])> {
    line.get(3..)
        .unwrap_or_default()
        .split(|&b| b == b'\t')
        .skip(1)
        .filter_map(|field| match field {
            [a, b, b':', value @ ..] => Some(([*a, *b], value)),
            _ => None,
        })
}

/// Returns the value of the first field with the given tag in a header line.
pub(crate) fn field(line: &[u8], tag: [u8; 2]) -> Option<&[u8]> {
    fields(line)
        .find(|(t, _)| *t == tag)
        .map(|(_, value)| value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fields() {
        let line = b"@RG\tID:rg0\tSM:\tSM:s1";
        let actual: Vec<_> = fields(line).collect();
        assert_eq!(
            actual,
            [
                (*b"ID", &b"rg0"[..]),
                (*b"SM", &b""[..]),
                (*b"SM", &b"s1"[..])
            ]
        );

        assert_eq!(field(line, *b"SM"), Some(&b""[..]));
        assert_eq!(field(line, *b"PP"), None);
        assert_eq!(record_type(line), *b"RG");
    }
}
