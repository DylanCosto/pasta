//! A record's tags as BAM encodes them, which is how CRAM stores their values: for BAM input, the
//! record's own bytes, kept without decoding a value or allocating for each.

use std::io;

use crate::bam::record::codec::encoder::data::field::{ty, write_value};
use crate::sam::alignment::{
    record::data::field::{Tag, Type},
    record_buf::data::field::Value as ValueBuf,
};

/// The tags of a record, one after another, each a tag, a type and a value as BAM encodes them.
/// They are checked when added, so reading them back cannot fail.
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Data(Vec<u8>);

/// A tag of [`Data`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct Field<'a> {
    pub(crate) tag: Tag,
    pub(crate) ty: Type,
    /// The value as BAM encodes it: a string with its NUL, an array with its subtype and length.
    pub(crate) value: &'a [u8],
}

impl Data {
    /// Builds the tags of encoded BAM tags `src`, but for those `keep` rejects.
    pub(crate) fn from_bam<F>(src: &[u8], mut keep: F) -> io::Result<Self>
    where
        F: FnMut(&Field<'_>) -> bool,
    {
        let mut dst = Vec::new();
        let mut rest = src;
        // The start of the tags not yet copied, so that runs of them are copied at once.
        let mut copy_start = 0;

        while !rest.is_empty() {
            let start = src.len() - rest.len();
            let field = split_field(&mut rest)?;

            if !keep(&field) {
                if dst.capacity() == 0 {
                    dst.reserve_exact(src.len());
                }

                dst.extend_from_slice(&src[copy_start..start]);
                copy_start = src.len() - rest.len();
            }
        }

        if copy_start == 0 {
            // Nothing was left out.
            return Ok(Self(src.to_vec()));
        }

        dst.extend_from_slice(&src[copy_start..]);
        Ok(Self(dst))
    }

    /// Adds a tag.
    pub(crate) fn push(&mut self, tag: Tag, value: &ValueBuf) -> io::Result<()> {
        self.0.extend(tag.as_ref());
        self.0.push(ty::encode(value.ty()));
        write_value(&mut self.0, &value.into())
    }

    #[inline]
    pub(crate) fn iter(&self) -> impl Iterator<Item = Field<'_>> {
        let mut rest = &self.0[..];

        std::iter::from_fn(move || {
            (!rest.is_empty()).then(|| split_field(&mut rest).expect("tags are checked when added"))
        })
    }

    /// Keeps only the tags `f` accepts, moving the others' bytes in place, a run of them at a
    /// time.
    pub(crate) fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(&Field<'_>) -> bool,
    {
        let mut read = 0;
        let mut write = 0;
        // The start of the tags kept but not yet moved.
        let mut run_start = 0;

        while read < self.0.len() {
            let mut rest = &self.0[read..];
            let field = split_field(&mut rest).expect("tags are checked when added");
            let end = self.0.len() - rest.len();

            if !f(&field) {
                self.0.copy_within(run_start..read, write);
                write += read - run_start;
                run_start = end;
            }

            read = end;
        }

        if run_start > 0 {
            self.0.copy_within(run_start.., write);
            write += self.0.len() - run_start;
            self.0.truncate(write);
        }
    }
}

impl Field<'_> {
    /// The value of an integer tag.
    pub(crate) fn as_int(&self) -> Option<i64> {
        let v = self.value;

        match self.ty {
            Type::Int8 => Some(i64::from(v[0] as i8)),
            Type::UInt8 => Some(i64::from(v[0])),
            Type::Int16 => Some(i64::from(i16::from_le_bytes([v[0], v[1]]))),
            Type::UInt16 => Some(i64::from(u16::from_le_bytes([v[0], v[1]]))),
            Type::Int32 => Some(i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]]))),
            Type::UInt32 => Some(i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]]))),
            _ => None,
        }
    }

    /// The value of a string tag, without its NUL.
    pub(crate) fn as_str(&self) -> Option<&[u8]> {
        matches!(self.ty, Type::String | Type::Hex).then(|| &self.value[..self.value.len() - 1])
    }
}

/// Splits the next tag off `src`, checking that it is whole.
#[inline]
fn split_field<'a>(src: &mut &'a [u8]) -> io::Result<Field<'a>> {
    let eof = || io::Error::from(io::ErrorKind::UnexpectedEof);

    let (&[a, b, t], rest) = src.split_first_chunk::<3>().ok_or_else(eof)?;
    let tag = Tag::new(a, b);

    let ty = match t {
        b'A' => Type::Character,
        b'c' => Type::Int8,
        b'C' => Type::UInt8,
        b's' => Type::Int16,
        b'S' => Type::UInt16,
        b'i' => Type::Int32,
        b'I' => Type::UInt32,
        b'f' => Type::Float,
        b'Z' => Type::String,
        b'H' => Type::Hex,
        b'B' => Type::Array,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
    };

    let len = match ty {
        Type::Character | Type::Int8 | Type::UInt8 => 1,
        Type::Int16 | Type::UInt16 => 2,
        Type::Int32 | Type::UInt32 | Type::Float => 4,
        Type::String | Type::Hex => memchr::memchr(0, rest).ok_or_else(eof)? + 1,
        Type::Array => {
            let (&[subtype], tail) = rest.split_first_chunk::<1>().ok_or_else(eof)?;
            let (count, _) = tail.split_first_chunk::<4>().ok_or_else(eof)?;
            let count = u32::from_le_bytes(*count) as usize;

            let size = match subtype {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid array subtype",
                    ));
                }
            };

            count
                .checked_mul(size)
                .and_then(|n| n.checked_add(5))
                .ok_or_else(eof)?
        }
    };

    if rest.len() < len {
        return Err(eof());
    }

    let (value, tail) = rest.split_at(len);
    *src = tail;

    Ok(Field { tag, ty, value })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bam() -> io::Result<()> {
        let mut src = Vec::new();
        src.extend(b"NMC\x02");
        src.extend(b"RGZrg0\0");
        src.extend(b"MDZ10A5\0");
        src.extend(b"MLBC\x02\x00\x00\x00\x07\x08");

        let data = Data::from_bam(&src, |_| true)?;
        assert_eq!(data.0, src);
        assert_eq!(data.iter().count(), 4);

        let fields: Vec<_> = data.iter().map(|f| (f.tag, f.ty)).collect();
        assert_eq!(
            fields,
            [
                (Tag::new(b'N', b'M'), Type::UInt8),
                (Tag::new(b'R', b'G'), Type::String),
                (Tag::new(b'M', b'D'), Type::String),
                (Tag::new(b'M', b'L'), Type::Array),
            ]
        );

        let nm = data.iter().next().unwrap();
        assert_eq!(nm.as_int(), Some(2));
        assert_eq!(data.iter().nth(2).unwrap().as_str(), Some(&b"10A5"[..]));

        let data = Data::from_bam(&src, |f| f.tag != Tag::new(b'R', b'G'))?;
        assert_eq!(data.iter().count(), 3);
        assert!(data.iter().all(|f| f.tag != Tag::new(b'R', b'G')));

        // Truncated values.
        for bad in [
            &b"NMC"[..],
            b"MDZ10A5",
            b"MLBC\x02\x00\x00\x00\x07",
            b"NMx\x02",
        ] {
            assert!(Data::from_bam(bad, |_| true).is_err(), "{bad:?}");
        }

        Ok(())
    }

    #[test]
    fn test_push_and_retain() -> io::Result<()> {
        let mut data = Data::default();
        data.push(Tag::new(b'N', b'M'), &ValueBuf::UInt8(1))?;
        data.push(Tag::new(b'C', b'O'), &ValueBuf::from("note"))?;
        assert_eq!(data.0, b"NMC\x01COZnote\0");

        data.retain(|f| f.tag != Tag::new(b'N', b'M'));
        assert_eq!(data.0, b"COZnote\0");

        data.push(Tag::new(b'N', b'M'), &ValueBuf::UInt8(2))?;
        data.push(Tag::new(b'A', b'S'), &ValueBuf::UInt8(3))?;
        data.push(Tag::new(b'X', b'S'), &ValueBuf::UInt8(4))?;
        data.retain(|f| f.tag != Tag::new(b'N', b'M') && f.tag != Tag::new(b'X', b'S'));
        assert_eq!(data.0, b"COZnote\0ASC\x03");

        data.retain(|_| true);
        assert_eq!(data.0, b"COZnote\0ASC\x03");

        data.retain(|_| false);
        assert!(data.0.is_empty());

        Ok(())
    }
}
