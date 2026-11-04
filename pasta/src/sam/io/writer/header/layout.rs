//! Writing a header that was read, as htslib does (see [`Layout`]).

use std::{
    collections::{HashMap, HashSet},
    io::{self, Write},
};

use bstr::{BStr, ByteSlice};

use super::record::{
    write_comment, write_header, write_program, write_read_group, write_reference_sequence,
};
use crate::sam::{
    Header,
    header::layout::{self, Kind, Layout},
};

const REFERENCE_SEQUENCE: [u8; 2] = *b"SQ";
const READ_GROUP: [u8; 2] = *b"RG";
const PROGRAM: [u8; 2] = *b"PG";
const COMMENT: [u8; 2] = *b"CO";

/// Writes a header that was read.
///
/// While the header is unchanged, this writes its text as read, if it has one (BAM). Otherwise,
/// it writes the lines as read, rewriting the records that have changed, and adds new records
/// after the last line of their type, or at the end. A new `@HD` record goes first.
pub(super) fn write_layout<W>(writer: &mut W, header: &Header, layout: &Layout) -> io::Result<()>
where
    W: Write,
{
    if let Some(text) = layout.text()
        && header == layout.snapshot()
        && (header.reference_sequences().is_empty() || has_reference_sequence_lines(layout))
    {
        writer.write_all(text)?;

        if !text.is_empty() && !text.ends_with(b"\n") {
            writer.write_all(b"\n")?;
        }

        return Ok(());
    }

    let snapshot = layout.snapshot();
    let lines = layout.lines();

    let mut last_lines = HashMap::new();
    let mut has_header_line = false;
    let mut reference_sequence_names = HashSet::new();
    let mut read_group_ids = HashSet::new();
    let mut program_ids = HashSet::new();
    let mut comment_count = 0;

    for (i, line) in lines.iter().enumerate() {
        last_lines.insert(layout::record_type(layout.line(line)), i);

        match line.kind() {
            Kind::Header => has_header_line = true,
            Kind::ReferenceSequence(name) => {
                reference_sequence_names.insert(name.as_ref());
            }
            Kind::ReadGroup(id) => {
                read_group_ids.insert(id.as_ref());
            }
            Kind::Program(id) => {
                program_ids.insert(id.as_ref());
            }
            Kind::Comment(_) => comment_count += 1,
            Kind::Other => {}
        }
    }

    let known = Known {
        reference_sequence_names,
        read_group_ids,
        program_ids,
        comment_count,
    };

    if !has_header_line && let Some(map) = header.header() {
        write_header(writer, map)?;
    }

    let mut buf = Vec::new();
    let mut old_buf = Vec::new();

    for (i, line) in lines.iter().enumerate() {
        let src = layout.line(line);

        match line.kind() {
            Kind::Header => {
                if let Some(map) = header.header() {
                    write_header(&mut buf, map)?;

                    if let Some(old) = snapshot.header() {
                        write_header(&mut old_buf, old)?;
                    }
                }
            }
            Kind::ReferenceSequence(name) => {
                if let Some(map) = header.reference_sequences().get(name) {
                    write_reference_sequence(&mut buf, name, map)?;

                    if let Some(old) = snapshot.reference_sequences().get(name) {
                        write_reference_sequence(&mut old_buf, name, old)?;
                    }
                }
            }
            Kind::ReadGroup(id) => {
                if let Some(map) = header.read_groups().get(id) {
                    write_read_group(&mut buf, id, map)?;

                    if let Some(old) = snapshot.read_groups().get(id) {
                        write_read_group(&mut old_buf, id, old)?;
                    }
                }
            }
            Kind::Program(id) => {
                if let Some(map) = header.programs().as_ref().get(id) {
                    write_program(&mut buf, id, map)?;

                    if let Some(old) = snapshot.programs().as_ref().get(id) {
                        write_program(&mut old_buf, id, old)?;
                    }
                }
            }
            Kind::Comment(j) => {
                if let Some(comment) = header.comments().get(*j) {
                    write_comment(writer, comment)?;
                }
            }
            Kind::Other => {
                if !layout::is_removed_duplicate(header, src) {
                    writer.write_all(src)?;
                    writer.write_all(b"\n")?;
                }
            }
        }

        if !buf.is_empty() {
            write_line(writer, src, &buf, &old_buf)?;
            buf.clear();
            old_buf.clear();
        }

        let record_type = layout::record_type(src);

        if last_lines.get(&record_type) == Some(&i) {
            write_new_records(writer, header, &known, record_type)?;
        }
    }

    for record_type in [REFERENCE_SEQUENCE, READ_GROUP, PROGRAM, COMMENT] {
        if !last_lines.contains_key(&record_type) {
            write_new_records(writer, header, &known, record_type)?;
        }
    }

    Ok(())
}

/// The records that have a line.
struct Known<'a> {
    reference_sequence_names: HashSet<&'a BStr>,
    read_group_ids: HashSet<&'a BStr>,
    program_ids: HashSet<&'a BStr>,
    comment_count: usize,
}

fn has_reference_sequence_lines(layout: &Layout) -> bool {
    layout
        .lines()
        .iter()
        .any(|line| matches!(line.kind(), Kind::ReferenceSequence(_)))
}

/// Writes the records of a type that have no line.
fn write_new_records<W>(
    writer: &mut W,
    header: &Header,
    known: &Known<'_>,
    record_type: [u8; 2],
) -> io::Result<()>
where
    W: Write,
{
    match record_type {
        REFERENCE_SEQUENCE => {
            for (name, map) in header.reference_sequences() {
                if !known.reference_sequence_names.contains(name.as_bstr()) {
                    write_reference_sequence(writer, name, map)?;
                }
            }
        }
        READ_GROUP => {
            for (id, map) in header.read_groups() {
                if !known.read_group_ids.contains(id.as_bstr()) {
                    write_read_group(writer, id, map)?;
                }
            }
        }
        PROGRAM => {
            for (id, map) in header.programs().as_ref() {
                if !known.program_ids.contains(id.as_bstr()) {
                    write_program(writer, id, map)?;
                }
            }
        }
        COMMENT => {
            for comment in header.comments().iter().skip(known.comment_count) {
                write_comment(writer, comment)?;
            }
        }
        _ => {}
    }

    Ok(())
}

/// Writes a record's line.
///
/// `record` is the record as it is now, and `old` as it was read (or empty if it was not), both
/// written. An unchanged record is written as read. Otherwise, as htslib updates a line, a
/// changed field replaces the value of the first field with its tag, a removed field is dropped,
/// and a new field is added at the end; other fields stay as read.
fn write_line<W>(writer: &mut W, src: &[u8], record: &[u8], old: &[u8]) -> io::Result<()>
where
    W: Write,
{
    if record == old {
        writer.write_all(src)?;
        return writer.write_all(b"\n");
    }

    if old.is_empty() {
        return writer.write_all(record);
    }

    let record = record.strip_suffix(b"\n").unwrap_or(record);
    let old = old.strip_suffix(b"\n").unwrap_or(old);

    let fields: Vec<_> = layout::fields(record).collect();
    let old_fields: Vec<_> = layout::fields(old).collect();

    writer.write_all(&record[..3])?;

    let mut replaced = Vec::new();

    for (tag, value) in layout::fields(src) {
        let new_value = get(&fields, tag);

        let value = if new_value == get(&old_fields, tag) {
            value
        } else if let Some(new_value) = new_value {
            if replaced.contains(&tag) {
                value
            } else {
                replaced.push(tag);
                new_value
            }
        } else {
            continue;
        };

        write_field(writer, tag, value)?;
    }

    for (tag, value) in &fields {
        if layout::field(src, *tag).is_none() && get(&old_fields, *tag) != Some(value) {
            write_field(writer, *tag, value)?;
        }
    }

    writer.write_all(b"\n")
}

fn get<'a>(fields: &[([u8; 2], &'a [u8])], tag: [u8; 2]) -> Option<&'a [u8]> {
    fields
        .iter()
        .find(|(t, _)| *t == tag)
        .map(|(_, value)| *value)
}

fn write_field<W>(writer: &mut W, tag: [u8; 2], value: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(b"\t")?;
    writer.write_all(&tag)?;
    writer.write_all(b":")?;
    writer.write_all(value)
}
