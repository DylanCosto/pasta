use std::io::{self, Write};

use rustc_hash::FxHashMap;

use super::write_byte_array_encoding;
use crate::{
    cram::{
        container::{
            block,
            compression_header::{
                Encoding, TagEncodings,
                encoding::codec::{Byte, ByteArray, Integer},
                preservation_map::tag_sets::Key,
            },
        },
        io::writer::{Record, collections::write_array, num::write_itf8},
    },
    sam::alignment::record::data::field::Type,
};

pub fn write_tag_encodings<W>(writer: &mut W, tag_encodings: &TagEncodings) -> io::Result<()>
where
    W: Write,
{
    let buf = encode(tag_encodings)?;
    write_array(writer, &buf)
}

fn encode(tag_encodings: &TagEncodings) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    encode_inner(&mut buf, tag_encodings)?;
    Ok(buf)
}

fn encode_inner<W>(writer: &mut W, tag_encodings: &TagEncodings) -> io::Result<()>
where
    W: Write,
{
    let len = i32::try_from(tag_encodings.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, len)?;

    // In block content ID order, so that the output does not depend on the map's order.
    let mut tag_encodings: Vec<_> = tag_encodings.iter().collect();
    tag_encodings.sort_unstable_by_key(|(block_content_id, _)| **block_content_id);

    for (&block_content_id, encoding) in tag_encodings {
        write_itf8(writer, block_content_id)?;
        write_byte_array_encoding(writer, encoding)?;
    }

    Ok(())
}

/// The byte that ends string values, as in HTSlib.
const STOP_BYTE: u8 = b'\t';

/// Builds the encodings of the records' tags, as HTSlib's `cram_encode_aux` does: each tag's
/// values go to its own external block, a value of a fixed size with no length (a Huffman code
/// of one symbol, which takes no bits), a string ended by a tab, and an array after its length.
/// A length before each value would be a varying byte between them, which compresses worse.
pub(super) fn build_tag_encodings(records: &[Record]) -> TagEncodings {
    // The type of each tag's values, and whether a string value holds the stop byte.
    let mut keys: FxHashMap<block::ContentId, (Type, bool)> = FxHashMap::default();

    for record in records {
        for field in record.data.iter() {
            let ty = field.ty;
            let block_content_id = block::ContentId::from(Key::new(field.tag, ty));

            let has_stop_byte = field.as_str().is_some_and(|s| s.contains(&STOP_BYTE));

            keys.entry(block_content_id).or_insert((ty, false)).1 |= has_stop_byte;
        }
    }

    keys.into_iter()
        .map(|(block_content_id, (ty, has_stop_byte))| {
            (
                block_content_id,
                build_tag_encoding(block_content_id, ty, has_stop_byte),
            )
        })
        .collect()
}

fn build_tag_encoding(
    block_content_id: block::ContentId,
    ty: Type,
    has_stop_byte: bool,
) -> Encoding<ByteArray> {
    let len = match ty {
        Type::Character | Type::Int8 | Type::UInt8 => Some(1),
        Type::Int16 | Type::UInt16 => Some(2),
        Type::Int32 | Type::UInt32 | Type::Float => Some(4),
        Type::String | Type::Hex | Type::Array => None,
    };

    if let Some(len) = len {
        return Encoding::new(ByteArray::ByteArrayLength {
            len_encoding: Encoding::new(Integer::Huffman {
                alphabet: vec![len],
                bit_lens: vec![0],
            }),
            value_encoding: Encoding::new(Byte::External { block_content_id }),
        });
    }

    // A string value is stored with its NUL, then the stop byte, unless a value (against the
    // SAM specification) holds the stop byte.
    if matches!(ty, Type::String | Type::Hex) && !has_stop_byte {
        return Encoding::new(ByteArray::ByteArrayStop {
            stop_byte: STOP_BYTE,
            block_content_id,
        });
    }

    Encoding::new(ByteArray::ByteArrayLength {
        len_encoding: Encoding::new(Integer::External { block_content_id }),
        value_encoding: Encoding::new(Byte::External { block_content_id }),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cram::io::writer::record::Data;
    use crate::sam::alignment::{record::data::field::Tag, record_buf::data::field::Value};

    fn data_of(fields: &[(Tag, Value)]) -> Data {
        let mut data = Data::default();

        for (tag, value) in fields {
            data.push(*tag, value).unwrap();
        }

        data
    }

    #[test]
    fn test_build_tag_encodings() {
        let nm = Tag::new(b'N', b'M');
        let mc = Tag::new(b'M', b'C');
        let xs = Tag::new(b'x', b's');
        let ml = Tag::new(b'M', b'L');

        let records = [
            Record {
                // A tab in a string is not valid SAM, so it is given as BAM input might hold it.
                data: Data::from_bam(
                    b"NMs\x03\x00MCZ151M\0xsZa\tb\0MLBC\x02\x00\x00\x00\x01\x02",
                    |_| true,
                )
                .unwrap(),
                ..Default::default()
            },
            Record {
                data: data_of(&[(xs, Value::from("c"))]),
                ..Default::default()
            },
        ];

        let encodings = build_tag_encodings(&records);
        let get = |tag, ty| &encodings[&block::ContentId::from(Key::new(tag, ty))];
        let length_in_block = |block_content_id| {
            Encoding::new(ByteArray::ByteArrayLength {
                len_encoding: Encoding::new(Integer::External { block_content_id }),
                value_encoding: Encoding::new(Byte::External { block_content_id }),
            })
        };

        let id = block::ContentId::from(Key::new(nm, Type::Int16));
        assert_eq!(
            get(nm, Type::Int16),
            &Encoding::new(ByteArray::ByteArrayLength {
                len_encoding: Encoding::new(Integer::Huffman {
                    alphabet: vec![2],
                    bit_lens: vec![0],
                }),
                value_encoding: Encoding::new(Byte::External {
                    block_content_id: id
                }),
            })
        );

        let id = block::ContentId::from(Key::new(mc, Type::String));
        assert_eq!(
            get(mc, Type::String),
            &Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: b'\t',
                block_content_id: id,
            })
        );

        // A value holds the stop byte.
        let id = block::ContentId::from(Key::new(xs, Type::String));
        assert_eq!(get(xs, Type::String), &length_in_block(id));

        let id = block::ContentId::from(Key::new(ml, Type::Array));
        assert_eq!(get(ml, Type::Array), &length_in_block(id));
    }
}
