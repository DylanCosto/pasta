use std::io;

use crate::cram::{
    container::{
        block,
        compression_header::encoding::{Decode, Encode},
    },
    huffman::CanonicalHuffmanDecoder,
    io::{
        BitReader, BitWriter,
        reader::{container::slice::records::ExternalDataReaders, num::read_itf8},
        writer::{container::slice::records::ExternalDataWriters, num::write_itf8},
    },
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Integer {
    External {
        block_content_id: block::ContentId,
    },
    Golomb {
        offset: i32,
        m: i32,
    },
    Huffman {
        alphabet: Vec<i32>,
        bit_lens: Vec<u32>,
    },
    Beta {
        offset: i32,
        len: u32,
    },
    Subexp {
        offset: i32,
        k: i32,
    },
    GolombRice {
        offset: i32,
        log2_m: i32,
    },
    Gamma {
        offset: i32,
    },
}

impl<'de> Decode<'de> for Integer {
    type Value = i32;

    fn decode(
        &self,
        core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
    ) -> io::Result<Self::Value> {
        match self {
            Self::External { block_content_id } => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                read_itf8(src)
            }
            Self::Huffman { alphabet, bit_lens } => {
                if let [symbol] = alphabet[..] {
                    Ok(symbol)
                } else {
                    let decoder = CanonicalHuffmanDecoder::new(alphabet, bit_lens)?;
                    decoder.decode(core_data_reader)
                }
            }
            Self::Beta { offset, len } => decode_beta(core_data_reader, *offset, *len),
            Self::Subexp { offset, k } => decode_subexp(core_data_reader, *offset, *k),
            Self::Gamma { offset } => decode_gamma(core_data_reader, *offset),
            // HTSlib does not decode these either (`cram_decoder_init` has no decoder for them).
            Self::Golomb { .. } | Self::GolombRice { .. } => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!("unsupported integer codec: {self:?}"),
            )),
        }
    }
}

/// Decodes a beta-coded value: `len` bits, less `offset`.
///
/// Like HTSlib, this reads up to 32 bits and wraps on overflow.
pub(super) fn decode_beta(reader: &mut BitReader<'_>, offset: i32, len: u32) -> io::Result<i32> {
    reader.read_i32(len).map(|n| n.wrapping_sub(offset))
}

/// Decodes a subexponential-coded value, as HTSlib's `cram_subexp_decode`.
///
/// A unary prefix of `u` 1 bits ends with a 0 bit. If `u` is 0, the value is the next `k` bits;
/// otherwise, it is 2^(u + k - 1) plus the next u + k - 1 bits. The result is that value less
/// `offset`.
fn decode_subexp(reader: &mut BitReader<'_>, offset: i32, k: i32) -> io::Result<i32> {
    let k = u32::try_from(k)
        .ok()
        .filter(|&k| k <= MAX_SUBEXP_K)
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("invalid subexp k: {k}"))
        })?;

    let mut u = 0;

    while reader.read_bit()? == 1 {
        u += 1;

        // The value must fit in 31 bits.
        if u > MAX_SUBEXP_K - k {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "subexp value is too large",
            ));
        }
    }

    // With u <= 31 - k, the value is less than 2^31.
    let n = if u == 0 {
        reader.read_i32(k)?
    } else {
        let len = u + k - 1;
        (1 << len) + reader.read_i32(len)?
    };

    Ok(n.wrapping_sub(offset))
}

/// Decodes an Elias gamma-coded value, as HTSlib's `cram_gamma_decode`.
///
/// A prefix of `n` 0 bits ends with a 1 bit, which is the most significant bit of the
/// (n + 1)-bit value. The result is that value less `offset`.
fn decode_gamma(reader: &mut BitReader<'_>, offset: i32) -> io::Result<i32> {
    let mut n = 0;

    while reader.read_bit()? == 0 {
        n += 1;

        // HTSlib reads the value into a 32-bit int.
        if n > MAX_GAMMA_LEN {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "gamma value is too large",
            ));
        }
    }

    let m = reader.read_i32(n)?;
    // For n = 31, this wraps, as in HTSlib.
    let x = (1_i32 << n) | m;

    Ok(x.wrapping_sub(offset))
}

/// The largest subexp `k` HTSlib accepts (`cram_subexp_decode_init`), which also bounds the bits
/// in a value.
const MAX_SUBEXP_K: u32 = 31;

/// The most bits after the gamma prefix that fit in an `i32`.
const MAX_GAMMA_LEN: u32 = 31;

impl Encode<'_> for Integer {
    type Value = i32;

    fn encode(
        &self,
        _core_data_writer: &mut BitWriter,
        external_data_writers: &mut ExternalDataWriters,
        value: Self::Value,
    ) -> io::Result<()> {
        match self {
            Self::External { block_content_id } => {
                let dst = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                write_itf8(dst, value)
            }
            // A code of one symbol takes no bits, e.g., the length of a tag value of a fixed size.
            Self::Huffman { alphabet, bit_lens } if bit_lens[..] == [0] => {
                if alphabet[..] == [value] {
                    Ok(())
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{value} is not the Huffman code's symbol, {alphabet:?}"),
                    ))
                }
            }
            // Otherwise, pasta writes integers to external blocks only.
            _ => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                format!("encoding integers is not supported for {self:?}"),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cram::container::compression_header::Encoding;

    #[test]
    fn test_decode() -> io::Result<()> {
        fn t(
            core_data: Option<&[u8]>,
            encoding: &Encoding<Integer>,
            expected: i32,
        ) -> io::Result<()> {
            let core_data = core_data.unwrap_or(&[0b10000000]);
            let mut core_data_reader = BitReader::new(core_data);

            let external_data = [0x0d];
            let mut external_data_readers = ExternalDataReaders::new();
            external_data_readers.insert(1, &external_data[..]);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, actual);

            Ok(())
        }

        t(
            None,
            &Encoding::new(Integer::External {
                block_content_id: 1,
            }),
            13,
        )?;
        t(
            None,
            &Encoding::new(Integer::Huffman {
                alphabet: vec![0x4e],
                bit_lens: vec![0],
            }),
            0x4e,
        )?;
        t(None, &Encoding::new(Integer::Beta { offset: 1, len: 3 }), 3)?;
        t(
            Some(&[0b00011010]),
            &Encoding::new(Integer::Gamma { offset: 5 }),
            8,
        )?;
        t(
            Some(&[0xff, 0xff, 0xff, 0xff]),
            &Encoding::new(Integer::Beta { offset: 0, len: 32 }),
            -1,
        )?;

        Ok(())
    }

    #[test]
    fn test_decode_subexp() -> io::Result<()> {
        fn t(core_data: &[u8], offset: i32, k: i32, expected: &[i32]) -> io::Result<()> {
            let mut core_data_reader = BitReader::new(core_data);
            let mut external_data_readers = ExternalDataReaders::new();
            let encoding = Encoding::new(Integer::Subexp { offset, k });

            for &value in expected {
                let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
                assert_eq!(actual, value, "offset {offset}, k {k}");
            }

            Ok(())
        }

        // k = 0: 0 => 0; 10 => 2^0 = 1; 110 1 => 2^1 + 1 = 3; 1110 10 => 2^2 + 2 = 6.
        t(&[0b01011011, 0b11010000], 0, 0, &[0, 1, 3, 6])?;

        // k = 2, offset = 1: 0 11 => 3; 10 01 => 2^2 + 1 = 5; 110 101 => 2^3 + 5 = 13.
        t(&[0b01110011, 0b10101000], 1, 2, &[2, 4, 12])?;

        // k = 31, u = 0: the next 31 bits.
        t(&[0x7f, 0xff, 0xff, 0xff], 0, 31, &[i32::MAX])?;

        // k = 0, u = 31: 2^30 + 30 bits, the largest value. An offset of -1 wraps, as in HTSlib.
        let max = [0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xfc];
        t(&max, 0, 0, &[i32::MAX])?;
        t(&max, -1, 0, &[i32::MIN])?;

        Ok(())
    }

    #[test]
    fn test_decode_with_invalid_input() {
        fn t(core_data: &[u8], codec: Integer, kind: io::ErrorKind) {
            let mut core_data_reader = BitReader::new(core_data);
            let mut external_data_readers = ExternalDataReaders::new();
            let encoding = Encoding::new(codec.clone());

            assert!(
                matches!(
                    encoding.decode(&mut core_data_reader, &mut external_data_readers),
                    Err(e) if e.kind() == kind
                ),
                "{codec:?}"
            );
        }

        // A unary prefix longer than the value allows.
        t(
            &[0xff; 5],
            Integer::Subexp { offset: 0, k: 0 },
            io::ErrorKind::InvalidData,
        );
        t(
            &[0x80],
            Integer::Subexp { offset: 0, k: 31 },
            io::ErrorKind::InvalidData,
        );
        t(
            &[0x00],
            Integer::Subexp { offset: 0, k: -1 },
            io::ErrorKind::InvalidData,
        );
        t(
            &[0x00],
            Integer::Subexp { offset: 0, k: 32 },
            io::ErrorKind::InvalidData,
        );
        t(
            &[0xf0],
            Integer::Subexp { offset: 0, k: 2 },
            io::ErrorKind::UnexpectedEof,
        );

        t(
            &[0x00; 5],
            Integer::Gamma { offset: 0 },
            io::ErrorKind::InvalidData,
        );
        t(
            &[0x00],
            Integer::Gamma { offset: 0 },
            io::ErrorKind::UnexpectedEof,
        );

        t(
            &[0x01],
            Integer::Huffman {
                alphabet: Vec::new(),
                bit_lens: Vec::new(),
            },
            io::ErrorKind::InvalidData,
        );
        t(
            &[0x01],
            Integer::Huffman {
                alphabet: vec![1, 2],
                bit_lens: vec![1, 40],
            },
            io::ErrorKind::InvalidData,
        );

        t(
            &[0x00],
            Integer::Golomb { offset: 0, m: 2 },
            io::ErrorKind::Unsupported,
        );
        t(
            &[0x00],
            Integer::GolombRice {
                offset: 0,
                log2_m: 1,
            },
            io::ErrorKind::Unsupported,
        );
    }

    #[test]
    fn test_decode_core_codecs_never_panic() {
        let offsets = [i32::MIN, -1, 0, 1, i32::MAX];

        let mut codecs = Vec::new();

        for offset in offsets {
            codecs.push(Integer::Gamma { offset });

            for k in -1..=33 {
                codecs.push(Integer::Subexp { offset, k });
            }

            for len in 0..=33 {
                codecs.push(Integer::Beta { offset, len });
            }
        }

        for codec in codecs {
            let encoding = Encoding::new(codec);

            for b in 0..=u8::MAX {
                let core_data = [b; 9];
                let mut core_data_reader = BitReader::new(&core_data);
                let mut external_data_readers = ExternalDataReaders::new();

                // Each decode consumes at least one bit or fails, except beta with no bits.
                for _ in 0..=core_data.len() * 8 {
                    if encoding
                        .decode(&mut core_data_reader, &mut external_data_readers)
                        .is_err()
                    {
                        break;
                    }
                }
            }
        }
    }

    #[test]
    fn test_encode() -> io::Result<()> {
        fn t(
            encoding: &Encoding<Integer>,
            value: i32,
            expected_core_data: &[u8],
            expected_external_data: &[u8],
        ) -> io::Result<()> {
            let mut core_data_writer = BitWriter::default();

            let block_content_id = 1;
            let mut external_data_writers = [(block_content_id, Vec::new())].into_iter().collect();

            encoding.encode(&mut core_data_writer, &mut external_data_writers, value)?;

            let actual_core_data = core_data_writer.finish()?;
            assert_eq!(actual_core_data, expected_core_data);

            let actual_external_data = &external_data_writers[&block_content_id];
            assert_eq!(actual_external_data, expected_external_data);

            Ok(())
        }

        t(
            &Encoding::new(Integer::External {
                block_content_id: 1,
            }),
            0x0d,
            &[],
            &[0x0d],
        )?;

        let constant = Encoding::new(Integer::Huffman {
            alphabet: vec![4],
            bit_lens: vec![0],
        });
        t(&constant, 4, &[], &[])?;
        assert!(matches!(
            t(&constant, 2, &[], &[]),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        assert!(matches!(
            t(&Encoding::new(Integer::Subexp { offset: 0, k: 0 }), 1, &[], &[]),
            Err(e) if e.kind() == io::ErrorKind::Unsupported
        ));

        Ok(())
    }
}
