use std::io::{self, Write};

use crate::cram::{container::compression_header::preservation_map::TagSets, io::writer::Record};

use super::write_array;

pub(super) fn write_tag_sets<W>(writer: &mut W, tag_sets: &TagSets) -> io::Result<()>
where
    W: Write,
{
    let buf = encode(tag_sets);
    write_array(writer, &buf)
}

fn encode(tag_sets: &TagSets) -> Vec<u8> {
    use crate::bam::record::codec::encoder::data::field::ty;

    const NUL: u8 = 0x00;

    let mut buf = Vec::new();

    for keys in tag_sets.iter() {
        for key in keys {
            let tag = key.tag();
            buf.extend(tag.as_ref());
            buf.push(ty::encode(key.ty()));
        }

        buf.push(NUL);
    }

    buf
}

pub(super) fn build_tag_sets(records: &[Record]) -> TagSets {
    use crate::cram::container::compression_header::preservation_map::tag_sets::Key;

    let mut tag_sets = Vec::new();

    for record in records {
        if find_tag_set(&tag_sets, record).is_none() {
            let set = record
                .data
                .iter()
                .map(|field| Key::new(field.tag, field.ty))
                .collect();

            tag_sets.push(set);
        }
    }

    tag_sets
}

/// Returns the index of the record's tag set (its tags and their types, in order) in
/// `tag_sets`, without collecting it.
pub(crate) fn find_tag_set(tag_sets: &TagSets, record: &Record) -> Option<usize> {
    use crate::cram::container::compression_header::preservation_map::tag_sets::Key;

    tag_sets.iter().position(|set| {
        let mut fields = record.data.iter();

        set.iter().all(|key| {
            fields
                .next()
                .is_some_and(|field| *key == Key::new(field.tag, field.ty))
        }) && fields.next().is_none()
    })
}

#[cfg(test)]
mod tests {
    use crate::sam::alignment::record::data::field::{Tag, Type};

    use super::*;
    use crate::cram::container::compression_header::preservation_map::tag_sets::Key;

    #[test]
    fn test_write_tag_sets() -> io::Result<()> {
        let mut buf = Vec::new();

        let tag_sets = vec![
            vec![Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8)],
            vec![
                Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8),
                Key::new(Tag::COMMENT, Type::String),
            ],
        ];

        write_tag_sets(&mut buf, &tag_sets)?;

        let expected = [
            0x0b, // data_len
            0x4e, 0x48, 0x63, // NH:c
            0x00, // nul
            0x4e, 0x48, 0x63, // NH:c
            0x43, 0x4f, 0x5a, // CO:Z
            0x00, // nul
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
