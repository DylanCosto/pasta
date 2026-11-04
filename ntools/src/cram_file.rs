//! CRAM file structure at the container level: the file definition, container headers, and
//! header containers, for copying containers without decoding them (`reheader`, `cat`).
//!
//! Containers are read in order, without seeking, so inputs can be pipes. As in HTSlib, every
//! container header's CRC32 is checked; the blocks are copied without decoding them, so their
//! CRC32s are not (HTSlib checks a block's only when it decompresses it).
//!
//! CRAM 2.1 and 3.x are supported.

use std::io::{self, Read, Seek, SeekFrom, Write};

use anyhow::{Context, bail};
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use pasta::sam;

use crate::cram_raw::{invalid_data, read_itf8, read_len, read_ltf8, take};

/// The length of the file definition: magic number, version, and file ID.
pub const FILE_DEFINITION_LEN: usize = 26;

/// The block content type of the SAM header.
const FILE_HEADER_CONTENT_TYPE: u8 = 0;

/// Block compression methods.
pub const RAW: u8 = 0;
pub const GZIP: u8 = 1;

/// The reference sequence ID of multi-reference containers.
pub const MULTI_REFERENCE: i32 = -2;

/// The start position of EOF containers ("EOF" in ASCII).
const EOF_START: i32 = 0x0045_4f46;

/// The CRAM 3 EOF container.
pub const EOF_3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// The CRAM 2.1 EOF container.
pub const EOF_2: [u8; 30] = [
    0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
];

/// A CRAM file definition.
#[derive(Clone, Copy)]
pub struct FileDefinition(pub [u8; FILE_DEFINITION_LEN]);

impl FileDefinition {
    pub fn read<R: Read>(reader: &mut R) -> anyhow::Result<Self> {
        let mut buf = [0; FILE_DEFINITION_LEN];
        reader
            .read_exact(&mut buf)
            .context("reading the CRAM file definition")?;

        if &buf[..4] != b"CRAM" {
            bail!("not a CRAM file");
        }

        if !(2..=3).contains(&buf[4]) {
            bail!("unsupported CRAM version {}", buf[4]);
        }

        Ok(Self(buf))
    }

    pub fn major_version(&self) -> u8 {
        self.0[4]
    }

    /// Whether containers and blocks carry CRC32 checksums (CRAM 3).
    pub fn checksums(&self) -> Checksums {
        if self.major_version() >= 3 {
            Checksums::Crc32
        } else {
            Checksums::None
        }
    }
}

/// Whether CRAM containers and blocks end with a checksum.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Checksums {
    /// CRAM 3: a CRC32 of the preceding bytes.
    Crc32,
    /// CRAM 2.
    None,
}

/// A container header.
pub struct ContainerHeader {
    pub reference_sequence_id: i32,
    pub start: i32,
    pub span: i32,
    pub record_count: i32,
    /// The bytes from the reference sequence ID through the base count, as read.
    pub fields: Vec<u8>,
    pub block_count: usize,
    pub landmark_count: usize,
    /// The encoded header.
    pub raw: Vec<u8>,
    /// The length of the blocks that follow the header.
    pub data_len: usize,
}

impl ContainerHeader {
    /// Reads a container header, or returns `None` at the end of the stream.
    pub fn read<R: Read>(reader: &mut R, definition: &FileDefinition) -> io::Result<Option<Self>> {
        let mut raw = Vec::new();

        let mut len = [0; 4];

        match reader.read(&mut len[..1])? {
            0 => return Ok(None),
            _ => reader.read_exact(&mut len[1..])?,
        }

        raw.extend_from_slice(&len);

        // Reference sequence ID, start, span, record count.
        for _ in 0..4 {
            read_varint(reader, &mut raw, Varint::Itf8)?;
        }

        // Record counter and base count.
        let counter = if definition.major_version() >= 3 {
            Varint::Ltf8
        } else {
            Varint::Itf8
        };

        read_varint(reader, &mut raw, counter)?;
        read_varint(reader, &mut raw, Varint::Ltf8)?;

        let fields_end = raw.len();

        // Block count and landmark count.
        read_varint(reader, &mut raw, Varint::Itf8)?;
        read_varint(reader, &mut raw, Varint::Itf8)?;

        let mut src = &raw[fields_end..];
        let block_count = to_count(read_itf8(&mut src)?)?;
        let landmark_count = to_count(read_itf8(&mut src)?)?;

        for _ in 0..landmark_count {
            read_varint(reader, &mut raw, Varint::Itf8)?;
        }

        if definition.checksums() == Checksums::Crc32 {
            let mut crc = [0; 4];
            reader.read_exact(&mut crc)?;

            let expected = u32::from_le_bytes(crc);
            let actual = crc32fast::hash(&raw);

            if actual != expected {
                return Err(invalid_data(format!(
                    "container header CRC32 mismatch: stored {expected:08x}, computed {actual:08x}"
                )));
            }

            raw.extend_from_slice(&crc);
        }

        let data_len = to_count(i32::from_le_bytes(len))?;

        let mut src = &raw[4..fields_end];
        let reference_sequence_id = read_itf8(&mut src)?;
        let start = read_itf8(&mut src)?;
        let span = read_itf8(&mut src)?;
        let record_count = read_itf8(&mut src)?;

        Ok(Some(Self {
            reference_sequence_id,
            start,
            span,
            record_count,
            fields: raw[4..fields_end].to_vec(),
            block_count,
            landmark_count,
            raw,
            data_len,
        }))
    }

    /// The encoded length of the header and its blocks.
    pub fn total_len(&self) -> usize {
        self.raw.len() + self.data_len
    }

    /// Whether this is an EOF container (HTSlib's `empty_container`).
    pub fn is_eof(&self) -> bool {
        self.record_count == 0 && self.reference_sequence_id == -1 && self.start == EOF_START
    }
}

/// The start of a CRAM file: its file definition and header container.
pub struct FileStart {
    pub definition: FileDefinition,
    pub container: ContainerHeader,
    /// The header container's blocks.
    pub data: Vec<u8>,
}

impl FileStart {
    /// Reads the file definition and the header container, leaving `reader` at the first data
    /// container.
    pub fn read<R: Read>(reader: &mut R) -> anyhow::Result<Self> {
        let definition = FileDefinition::read(reader)?;
        let container = ContainerHeader::read(reader, &definition)
            .context("reading the CRAM header container")?
            .context("missing CRAM header container")?;

        // HTSlib can write a CRAM 2.1 header container whose length is 2 bytes short of its
        // padded block (it sizes the block's ITF-8 lengths before padding), and reads it by
        // reading the blocks, then any rest of the length. So does this.
        let mut data = Vec::new();

        for _ in 0..container.block_count {
            read_block(reader, definition.checksums(), &mut data)
                .context("reading the CRAM header container")?;
        }

        if let Some(padding) = container.data_len.checked_sub(data.len()) {
            let mut padding = vec![0; padding];
            reader
                .read_exact(&mut padding)
                .context("reading the CRAM header container")?;
            data.extend_from_slice(&padding);
        }

        Ok(Self {
            definition,
            container,
            data,
        })
    }

    /// The length of the file definition and the header container, which is the offset of the
    /// first data container.
    pub fn len(&self) -> usize {
        FILE_DEFINITION_LEN + self.container.raw.len() + self.data.len()
    }

    /// Decodes the SAM header.
    pub fn header(&self) -> anyhow::Result<sam::Header> {
        if self.definition.major_version() >= 3 {
            let mut src = Vec::with_capacity(self.len());
            src.extend_from_slice(&self.definition.0);
            src.extend_from_slice(&self.container.raw);
            src.extend_from_slice(&self.data);

            return Ok(pasta::cram::io::Reader::new(src.as_slice()).read_header()?);
        }

        // CRAM 2 has no checksums, and pasta does not read it: decode the header block here.
        let mut src = self.data.as_slice();
        let [method, content_type] = take(&mut src, 2)?.try_into().map_err(invalid_data)?;
        read_itf8(&mut src)?;
        let compressed_len = read_len(&mut src)?;
        read_len(&mut src)?;
        let compressed = take(&mut src, compressed_len)?;

        if content_type != FILE_HEADER_CONTENT_TYPE {
            bail!("the CRAM header container does not start with a file header block");
        }

        let data = match method {
            RAW => compressed.to_vec(),
            GZIP => {
                let mut data = Vec::new();
                MultiGzDecoder::new(compressed).read_to_end(&mut data)?;
                data
            }
            _ => bail!("unsupported CRAM 2 header block compression method {method}"),
        };

        let mut src = data.as_slice();
        let text_len = to_count(i32::from_le_bytes(
            take(&mut src, 4)?.try_into().map_err(invalid_data)?,
        ))?;
        let text = take(&mut src, text_len)?;

        Ok(sam::io::Reader::new(text).read_header()?)
    }
}

/// Reads one block, as it is, into `buf`, after what `buf` holds.
fn read_block<R: Read>(reader: &mut R, checksums: Checksums, buf: &mut Vec<u8>) -> io::Result<()> {
    // The compression method and content type.
    let mut prefix = [0; 2];
    reader.read_exact(&mut prefix)?;
    buf.extend_from_slice(&prefix);

    // The content ID, then the compressed and uncompressed sizes.
    read_varint(reader, buf, Varint::Itf8)?;
    let compressed_len_start = buf.len();
    read_varint(reader, buf, Varint::Itf8)?;
    let compressed_len = read_len(&mut &buf[compressed_len_start..])?;
    read_varint(reader, buf, Varint::Itf8)?;

    let crc32_len = if checksums == Checksums::Crc32 { 4 } else { 0 };
    let len = compressed_len + crc32_len;

    let start = buf.len();
    buf.resize(start + len, 0);
    reader.read_exact(&mut buf[start..])
}

/// Reads the blocks of the container whose header is `container` into `buf`, after what `buf`
/// holds.
pub fn read_data<R: Read>(
    reader: &mut R,
    container: &ContainerHeader,
    buf: &mut Vec<u8>,
) -> io::Result<()> {
    let len = container.data_len as u64;
    let n = reader.take(len).read_to_end(buf)?;

    if n as u64 == len {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("truncated container: {n} of {len} bytes"),
        ))
    }
}

/// Reads the containers after a CRAM file's header container, keeping track of their offsets
/// for messages and indexes.
pub struct Containers<R> {
    reader: R,
    definition: FileDefinition,
    /// The offset of the next byte to read.
    offset: u64,
    /// The offset of the container being read.
    container_offset: u64,
}

impl<R: Read> Containers<R> {
    /// Reads the start of a CRAM file, returning it with a reader of the containers after it.
    pub fn open(mut reader: R) -> anyhow::Result<(FileStart, Self)> {
        let start = FileStart::read(&mut reader)?;

        let offset = start.len() as u64;

        let containers = Self {
            reader,
            definition: start.definition,
            offset,
            container_offset: offset,
        };

        Ok((start, containers))
    }

    /// The offset of the next container from the start of the file.
    pub fn offset(&self) -> u64 {
        self.offset
    }

    /// Reads the next container header, or returns `None` at the end of the file.
    pub fn read_header(&mut self) -> io::Result<Option<ContainerHeader>> {
        self.container_offset = self.offset;

        let header = ContainerHeader::read(&mut self.reader, &self.definition).map_err(|e| {
            let e = if e.kind() == io::ErrorKind::UnexpectedEof {
                io::Error::new(e.kind(), "truncated container header")
            } else {
                e
            };

            self.error_at(e)
        })?;

        if let Some(header) = &header {
            self.offset += header.raw.len() as u64;
        }

        Ok(header)
    }

    /// Reads the blocks that follow `container`'s header into `buf`, after what `buf` holds.
    pub fn read_data(&mut self, container: &ContainerHeader, buf: &mut Vec<u8>) -> io::Result<()> {
        read_data(&mut self.reader, container, buf).map_err(|e| self.error_at(e))?;
        self.offset += container.data_len as u64;
        Ok(())
    }

    /// Skips the blocks that follow `container`'s header.
    pub fn skip_data(&mut self, container: &ContainerHeader) -> io::Result<()> {
        let len = container.data_len as u64;
        let n = io::copy(&mut (&mut self.reader).take(len), &mut io::sink())?;

        if n != len {
            return Err(self.error_at(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("truncated container: {n} of {len} bytes"),
            )));
        }

        self.offset += len;
        Ok(())
    }

    /// Copies every remaining container as it is, checking each header's CRC32.
    pub fn copy_to<W: Write>(&mut self, dst: &mut W) -> io::Result<()> {
        let mut buf = Vec::new();

        while let Some(container) = self.read_header()? {
            buf.clear();
            buf.extend_from_slice(&container.raw);
            self.read_data(&container, &mut buf)?;
            dst.write_all(&buf)?;
        }

        Ok(())
    }

    /// Adds the offset of the container being read to an error.
    fn error_at(&self, e: io::Error) -> io::Error {
        io::Error::new(
            e.kind(),
            format!(
                "reading the CRAM container at offset {}: {e}",
                self.container_offset
            ),
        )
    }
}

impl<R: Seek> Containers<R> {
    /// Moves to the container at `offset` from the start of the file.
    pub fn seek_to(&mut self, offset: u64) -> io::Result<()> {
        self.reader.seek(SeekFrom::Start(offset))?;
        self.offset = offset;
        self.container_offset = offset;
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum Varint {
    Itf8,
    Ltf8,
}

/// Reads one ITF-8 or LTF-8 integer's bytes into `raw`.
fn read_varint<R: Read>(reader: &mut R, raw: &mut Vec<u8>, kind: Varint) -> io::Result<()> {
    let mut b0 = [0; 1];
    reader.read_exact(&mut b0)?;
    raw.push(b0[0]);

    let ones = usize::try_from(b0[0].leading_ones()).unwrap_or(0);

    let extra = match kind {
        Varint::Itf8 => ones.min(4),
        Varint::Ltf8 => ones,
    };

    let start = raw.len();
    raw.resize(start + extra, 0);
    reader.read_exact(&mut raw[start..])?;

    // LTF-8 values are only checked to be well-formed.
    if let Varint::Ltf8 = kind {
        read_ltf8(&mut &raw[start - 1..])?;
    }

    Ok(())
}

fn to_count(n: i32) -> io::Result<usize> {
    usize::try_from(n).map_err(invalid_data)
}

fn to_i32(n: usize) -> io::Result<i32> {
    i32::try_from(n).map_err(invalid_data)
}

/// Encodes a header container (one gzip block) holding `header`.
pub fn encode_header_container(
    header: &sam::Header,
    definition: &FileDefinition,
) -> io::Result<Vec<u8>> {
    // Reference sequence ID, start, span, record count, record counter, and base count: all 0.
    const FIELDS: [u8; 6] = [0; 6];

    let block = encode_block(&header_data(header)?, GZIP, definition.checksums())?;
    let mut container =
        encode_container_header(&FIELDS, 1, &[], block.len(), definition.checksums())?;
    container.extend_from_slice(&block);

    Ok(container)
}

/// The header block's data: the length of the SAM header text, then the text.
pub fn header_data(header: &sam::Header) -> io::Result<Vec<u8>> {
    let mut text = Vec::new();
    sam::io::Writer::new(&mut text).write_header(header)?;

    let mut data = Vec::with_capacity(4 + text.len());
    data.extend_from_slice(&to_i32(text.len())?.to_le_bytes());
    data.extend_from_slice(&text);

    Ok(data)
}

/// Encodes a file header block, with a CRC32 for CRAM 3.
pub fn encode_block(data: &[u8], method: u8, checksums: Checksums) -> io::Result<Vec<u8>> {
    let compressed = match method {
        GZIP => {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::best());
            encoder.write_all(data)?;
            encoder.finish()?
        }
        _ => data.to_vec(),
    };

    let mut block = vec![method, FILE_HEADER_CONTENT_TYPE];
    write_itf8(&mut block, 0);
    write_itf8(&mut block, to_i32(compressed.len())?);
    write_itf8(&mut block, to_i32(data.len())?);
    block.extend_from_slice(&compressed);

    if checksums == Checksums::Crc32 {
        let crc = crc32fast::hash(&block);
        block.extend_from_slice(&crc.to_le_bytes());
    }

    Ok(block)
}

/// Encodes a container header from its raw fields (reference sequence ID through base count).
pub fn encode_container_header(
    fields: &[u8],
    block_count: i32,
    landmarks: &[i32],
    data_len: usize,
    checksums: Checksums,
) -> io::Result<Vec<u8>> {
    let mut buf = to_i32(data_len)?.to_le_bytes().to_vec();
    buf.extend_from_slice(fields);
    write_itf8(&mut buf, block_count);
    write_itf8(&mut buf, to_i32(landmarks.len())?);

    for &landmark in landmarks {
        write_itf8(&mut buf, landmark);
    }

    if checksums == Checksums::Crc32 {
        let crc = crc32fast::hash(&buf);
        buf.extend_from_slice(&crc.to_le_bytes());
    }

    Ok(buf)
}

/// Writes an ITF-8 integer.
fn write_itf8(dst: &mut Vec<u8>, n: i32) {
    let n = n.cast_unsigned();
    let [b0, b1, b2, b3] = n.to_be_bytes();

    if n < 0x80 {
        dst.push(b3);
    } else if n < 0x4000 {
        dst.extend_from_slice(&[0x80 | b2, b3]);
    } else if n < 0x20_0000 {
        dst.extend_from_slice(&[0xc0 | b1, b2, b3]);
    } else if n < 0x1000_0000 {
        dst.extend_from_slice(&[0xe0 | b0, b1, b2, b3]);
    } else {
        dst.extend_from_slice(&[
            0xf0 | (b0 >> 4),
            (b0 << 4) | (b1 >> 4),
            (b1 << 4) | (b2 >> 4),
            (b2 << 4) | (b3 >> 4),
            b3 & 0x0f,
        ]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_itf8() -> io::Result<()> {
        for n in [
            0,
            1,
            0x7f,
            0x80,
            0x3fff,
            0x4000,
            0x1f_ffff,
            0x20_0000,
            0x0fff_ffff,
            0x1000_0000,
            -1,
            i32::MAX,
        ] {
            let mut buf = Vec::new();
            write_itf8(&mut buf, n);
            let mut src = buf.as_slice();
            assert_eq!(read_itf8(&mut src)?, n);
            assert!(src.is_empty());
        }

        Ok(())
    }

    #[test]
    fn test_read_container_header() -> io::Result<()> {
        let definition = FileDefinition(*b"CRAM\x03\x00xxxxxxxxxxxxxxxxxxxx");
        let raw = encode_container_header(
            &[0x81, 0x00, 5, 7, 3, 0, 0x80, 0x80],
            2,
            &[0, 300],
            9,
            Checksums::Crc32,
        )?;

        let mut src = raw.as_slice();
        let header = ContainerHeader::read(&mut src, &definition)?.expect("header");

        assert!(src.is_empty());
        assert_eq!(header.reference_sequence_id, 0x100);
        assert_eq!(header.start, 5);
        assert_eq!(header.span, 7);
        assert_eq!(header.record_count, 3);
        assert_eq!(header.landmark_count, 2);
        assert_eq!(header.data_len, 9);
        assert_eq!(header.raw, raw);

        Ok(())
    }

    #[test]
    fn test_read_container_header_checks_crc32() -> io::Result<()> {
        let definition = FileDefinition(*b"CRAM\x03\x00xxxxxxxxxxxxxxxxxxxx");
        let raw = encode_container_header(&[0, 0, 0, 0, 0, 0], 1, &[0], 9, Checksums::Crc32)?;

        for i in 0..raw.len() {
            let mut corrupt = raw.clone();
            corrupt[i] ^= 0x01;

            let result = ContainerHeader::read(&mut corrupt.as_slice(), &definition);
            assert!(result.is_err(), "byte {i} flipped");
        }

        // CRAM 2 has no checksums.
        let definition = FileDefinition(*b"CRAM\x02\x01xxxxxxxxxxxxxxxxxxxx");
        let mut raw = encode_container_header(&[0; 6], 1, &[], 9, Checksums::None)?;
        raw[5] ^= 0x01;
        assert!(ContainerHeader::read(&mut raw.as_slice(), &definition)?.is_some());

        Ok(())
    }

    /// A CRAM file: the file definition, a header container, and `containers`.
    fn cram_file(
        definition: &FileDefinition,
        header: &sam::Header,
        containers: &[u8],
    ) -> io::Result<Vec<u8>> {
        let mut file = definition.0.to_vec();
        file.extend_from_slice(&encode_header_container(header, definition)?);
        file.extend_from_slice(containers);
        Ok(file)
    }

    #[test]
    fn test_file_start_header() -> anyhow::Result<()> {
        let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:8\n".parse()?;

        for version in [*b"\x02\x01", *b"\x03\x00", *b"\x03\x01"] {
            let mut definition = FileDefinition([0; FILE_DEFINITION_LEN]);
            definition.0[..4].copy_from_slice(b"CRAM");
            definition.0[4..6].copy_from_slice(&version);

            let file = cram_file(&definition, &header, &[])?;
            let start = FileStart::read(&mut file.as_slice())?;

            assert_eq!(start.len(), file.len(), "version {version:?}");
            assert_eq!(start.header()?, header, "version {version:?}");
        }

        Ok(())
    }

    #[test]
    fn test_file_start_reads_htslib_cram_2_1_header_containers() -> anyhow::Result<()> {
        let mut definition = FileDefinition([0; FILE_DEFINITION_LEN]);
        definition.0[..6].copy_from_slice(b"CRAM\x02\x01");

        // HTSlib pads a small header's block to 10,000 bytes, but sizes the container from the
        // block's 1-byte ITF-8 lengths before padding, so the length is 2 bytes short.
        let header: sam::Header = "@HD\tVN:1.6\n".parse()?;
        let mut data = header_data(&header)?;
        data.resize(10_000, 0);
        let block = encode_block(&data, RAW, Checksums::None)?;
        let container =
            encode_container_header(&[0; 6], 1, &[0], block.len() - 2, Checksums::None)?;

        let mut file = definition.0.to_vec();
        file.extend_from_slice(&container);
        file.extend_from_slice(&block);
        let data_start = file.len();
        file.extend_from_slice(&EOF_2);

        let (start, mut containers) = Containers::open(file.as_slice())?;
        assert_eq!(start.header()?, header);
        assert_eq!(containers.offset(), data_start as u64);

        let eof = containers.read_header()?.expect("EOF container");
        assert!(eof.is_eof());

        Ok(())
    }

    #[test]
    fn test_containers_copy_to() -> anyhow::Result<()> {
        let mut definition = FileDefinition([0; FILE_DEFINITION_LEN]);
        definition.0[..6].copy_from_slice(b"CRAM\x03\x01");
        let header: sam::Header = "@HD\tVN:1.6\n".parse()?;

        let mut containers =
            encode_container_header(&[0, 1, 2, 3, 4, 5], 1, &[0], 3, Checksums::Crc32)?;
        containers.extend_from_slice(b"abc");
        containers.extend_from_slice(&EOF_3);

        let file = cram_file(&definition, &header, &containers)?;
        let (_, mut reader) = Containers::open(file.as_slice())?;
        let mut dst = Vec::new();
        reader.copy_to(&mut dst)?;
        assert_eq!(dst, containers);
        assert_eq!(reader.offset(), file.len() as u64);

        // A container cut short, and one whose header is corrupt, are errors.
        let truncated = &file[..file.len() - EOF_3.len() - 1];
        let (_, mut reader) = Containers::open(truncated)?;
        let e = reader.copy_to(&mut Vec::new()).expect_err("truncated");
        assert_eq!(e.kind(), io::ErrorKind::UnexpectedEof);

        let mut corrupt = file.clone();
        let data_start = file.len() - containers.len();
        corrupt[data_start + 5] ^= 0x01;
        let (_, mut reader) = Containers::open(corrupt.as_slice())?;
        let e = reader.copy_to(&mut Vec::new()).expect_err("corrupt");
        assert!(
            e.to_string().contains(&format!("offset {data_start}")),
            "{e}"
        );
        assert!(e.to_string().contains("CRC32 mismatch"), "{e}");

        Ok(())
    }
}
