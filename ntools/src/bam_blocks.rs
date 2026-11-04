//! Rewriting BAM headers and concatenating BAM files by copying compressed BGZF blocks.
//!
//! Only the block holding the end of the header is decompressed and recompressed; every later
//! block is copied as is, as in `samtools reheader` and `samtools cat`.

use std::{
    io::{self, Read, Write},
    path::Path,
};

use anyhow::{Context, bail};
use flate2::read::DeflateDecoder;
use pasta::bam;
use pasta::bgzf;
use pasta::sam;

/// The BGZF end-of-file marker block.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// The BAM magic number.
const MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";

/// Length of a BGZF block header, up to and including `BSIZE`.
const BLOCK_HEADER_LEN: usize = 18;

/// Length of a BGZF block footer (`CRC32` and `ISIZE`).
const BLOCK_FOOTER_LEN: usize = 8;

/// A BAM file split at the end of its header.
pub struct Body<R> {
    reader: R,
    /// The input's path, for messages.
    path: std::path::PathBuf,
    /// Uncompressed record bytes that follow the header in its last block.
    leftover: Vec<u8>,
}

/// Reads a BAM file's header from `reader`, leaving it at the first block after the header.
/// The file is read in order, so it can be a pipe.
pub fn split<R: Read>(mut reader: R, path: &Path) -> anyhow::Result<(sam::Header, Body<R>)> {
    let mut data = Vec::new();
    let mut header_end = HeaderEnd::default();

    let len = loop {
        if let Some(len) = header_end.find(&data) {
            break len;
        }

        if data
            .get(..MAGIC_NUMBER.len())
            .is_some_and(|magic| magic != MAGIC_NUMBER)
        {
            bail!("reading header of {}: invalid BAM header", path.display());
        }

        let more = read_block(&mut reader, &mut data)
            .with_context(|| format!("reading header of {}", path.display()))?;

        if !more {
            bail!("reading header of {}: truncated BAM header", path.display());
        }
    };

    let header = bam::io::Reader::from(&data[..len])
        .read_header()
        .with_context(|| format!("reading header of {}", path.display()))?;

    data.drain(..len);

    Ok((
        header,
        Body {
            reader,
            path: path.to_path_buf(),
            leftover: data,
        },
    ))
}

impl<R: Read> Body<R> {
    /// Writes the records: the leftover bytes recompressed, then the remaining blocks copied.
    /// The input's EOF marker is not copied.
    ///
    /// Each block's framing is checked as it is copied, so an input cut inside a block is an
    /// error rather than a corrupt output with a fresh EOF marker; an input that ends between
    /// blocks without an EOF marker gets samtools' truncation warning.
    pub fn copy_to<W: Write>(mut self, dst: &mut W) -> anyhow::Result<()> {
        if !self.leftover.is_empty() {
            let mut writer = bgzf::io::Writer::new(&mut *dst);
            writer.write_all(&self.leftover)?;
            // Flushing writes the block; taking the inner writer skips the EOF marker that
            // finishing (or dropping) would write.
            writer.flush()?;
            writer.into_inner();
        }

        let has_eof = copy_blocks_without_eof(&mut self.reader, dst)
            .with_context(|| format!("copying {}", self.path.display()))?;

        if !has_eof {
            eprintln!(
                "[W::bam_hdr_read] EOF marker is absent. The input is probably truncated: {}",
                self.path.display()
            );
        }

        Ok(())
    }
}

/// Copies the BGZF blocks of `src` to `dst`, leaving out the EOF marker if `src` ends with one,
/// and returns whether it did.
///
/// Fails if `src` holds anything but whole BGZF blocks, e.g. if it is cut inside a block.
fn copy_blocks_without_eof<R: Read, W: Write>(src: &mut R, dst: &mut W) -> anyhow::Result<bool> {
    // A block is written only once another follows, so that a final EOF marker is left out.
    let mut block = Vec::new();
    let mut next = Vec::new();

    if !read_raw_block(src, &mut block)? {
        return Ok(false);
    }

    while read_raw_block(src, &mut next)? {
        dst.write_all(&block)?;
        std::mem::swap(&mut block, &mut next);
    }

    if block == BGZF_EOF {
        Ok(true)
    } else {
        dst.write_all(&block)?;
        Ok(false)
    }
}

/// Finds the end of a BAM header (magic number, text, and reference sequences) in its
/// uncompressed bytes as they are read, without parsing what was already scanned again.
#[derive(Default)]
struct HeaderEnd {
    /// The end of the fields scanned so far.
    position: usize,
    stage: Stage,
}

#[derive(Clone, Copy, Default)]
enum Stage {
    /// Before the text length.
    #[default]
    Text,
    /// Before the reference sequence count.
    ReferenceSequenceCount,
    /// Before the next of this many reference sequences.
    ReferenceSequences(u32),
}

impl HeaderEnd {
    /// Returns the length of the header if `src`, the start of a BAM stream, holds all of it.
    fn find(&mut self, src: &[u8]) -> Option<usize> {
        loop {
            match self.stage {
                Stage::Text => {
                    let text_len = read_len(src, self.position + MAGIC_NUMBER.len())?;
                    self.position += MAGIC_NUMBER.len() + 4 + text_len;
                    self.stage = Stage::ReferenceSequenceCount;
                }
                Stage::ReferenceSequenceCount => {
                    let n = read_u32(src, self.position)?;
                    self.position += 4;
                    self.stage = Stage::ReferenceSequences(n);
                }
                Stage::ReferenceSequences(0) => {
                    return (self.position <= src.len()).then_some(self.position);
                }
                Stage::ReferenceSequences(n) => {
                    let name_len = read_len(src, self.position)?;
                    // The name length, name, and sequence length.
                    self.position += 4 + name_len + 4;
                    self.stage = Stage::ReferenceSequences(n - 1);
                }
            }
        }
    }
}

fn read_u32(src: &[u8], position: usize) -> Option<u32> {
    let buf = src.get(position..position.checked_add(4)?)?;
    Some(u32::from_le_bytes(buf.try_into().ok()?))
}

fn read_len(src: &[u8], position: usize) -> Option<usize> {
    read_u32(src, position).and_then(|n| usize::try_from(n).ok())
}

/// Writes a BAM header in its own BGZF blocks, without an EOF marker.
pub fn write_header<W: Write>(dst: &mut W, header: &sam::Header) -> io::Result<()> {
    let mut bgzf_writer = bgzf::io::Writer::new(&mut *dst);
    bam::io::Writer::from(&mut bgzf_writer).write_header(header)?;
    bgzf_writer.flush()?;
    bgzf_writer.into_inner();
    Ok(())
}

pub fn write_eof<W: Write>(dst: &mut W) -> io::Result<()> {
    dst.write_all(&BGZF_EOF)
}

/// Reads a whole BGZF block, as stored, into `block` (replacing its contents), returning `false`
/// at the end of the stream.
fn read_raw_block<R: Read>(reader: &mut R, block: &mut Vec<u8>) -> anyhow::Result<bool> {
    let mut header = [0; BLOCK_HEADER_LEN];
    let mut len = 0;

    while len < header.len() {
        match reader.read(&mut header[len..]) {
            Ok(0) => break,
            Ok(n) => len += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e.into()),
        }
    }

    match len {
        0 => return Ok(false),
        BLOCK_HEADER_LEN => {}
        _ => bail!("truncated BGZF block header"),
    }

    if header[..4] != [0x1f, 0x8b, 0x08, 0x04] || header[12..14] != *b"BC" {
        bail!("not a BGZF block");
    }

    let block_len = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;

    if block_len < BLOCK_HEADER_LEN + BLOCK_FOOTER_LEN {
        bail!("invalid BGZF block size");
    }

    block.clear();
    block.extend_from_slice(&header);
    block.resize(block_len, 0);
    reader
        .read_exact(&mut block[BLOCK_HEADER_LEN..])
        .map_err(|e| match e.kind() {
            io::ErrorKind::UnexpectedEof => anyhow::anyhow!(
                "truncated BGZF block: expected {block_len} bytes; the input is probably truncated"
            ),
            _ => e.into(),
        })?;

    Ok(true)
}

/// Reads a BGZF block and appends its inflated data to `data`, returning `false` at the end of
/// the stream.
fn read_block<R: Read>(reader: &mut R, data: &mut Vec<u8>) -> anyhow::Result<bool> {
    let mut block = Vec::new();

    if !read_raw_block(reader, &mut block)? {
        return Ok(false);
    }

    let (deflated, footer) =
        block[BLOCK_HEADER_LEN..].split_at(block.len() - BLOCK_HEADER_LEN - BLOCK_FOOTER_LEN);
    let start = data.len();
    DeflateDecoder::new(deflated).read_to_end(data)?;

    let inflated = &data[start..];
    let crc32 = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
    let isize = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]);

    if usize::try_from(isize).ok() != Some(inflated.len()) || crc32fast::hash(inflated) != crc32 {
        bail!("BGZF block checksum or size mismatch");
    }

    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A reader that returns at most one byte per read, like a slow pipe.
    struct Trickle<'a>(&'a [u8]);

    impl Read for Trickle<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.0.len().min(buf.len()).min(1);
            buf[..n].copy_from_slice(&self.0[..n]);
            self.0 = &self.0[n..];
            Ok(n)
        }
    }

    /// A BAM file whose header spans several BGZF blocks and whose last header block also
    /// holds `records`.
    fn bam_file(header: &sam::Header, records: &[u8]) -> io::Result<Vec<u8>> {
        let mut data = Vec::new();
        bam::io::Writer::from(&mut data).write_header(header)?;
        data.extend_from_slice(records);

        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(&data)?;
        writer.finish()
    }

    #[test]
    fn test_split_reads_streams() -> anyhow::Result<()> {
        let mut text = String::from("@HD\tVN:1.6\n");

        for i in 0..4000 {
            text.push_str(&format!("@SQ\tSN:sequence{i}\tLN:1000\n"));
        }

        let header: sam::Header = text.parse()?;
        let records = b"not really records";
        let file = bam_file(&header, records)?;

        let (actual, body) = split(Trickle(&file), Path::new("in.bam"))?;
        assert_eq!(actual, header);
        assert_eq!(body.leftover, records);

        let mut dst = Vec::new();
        body.copy_to(&mut dst)?;

        let mut data = Vec::new();
        bgzf::io::Reader::new(dst.as_slice()).read_to_end(&mut data)?;
        assert_eq!(data, records);
        assert!(!dst.ends_with(&BGZF_EOF), "the EOF marker is left out");

        // A file cut inside the header.
        let e = split(&file[..file.len() / 2], Path::new("in.bam")).err();
        assert!(e.is_some_and(|e| format!("{e:#}").contains("truncated")));

        Ok(())
    }

    /// BGZF-compressed `data` in blocks of at most `block_len` uncompressed bytes, with an EOF
    /// marker.
    fn bgzf_blocks(data: &[u8], block_len: usize) -> io::Result<Vec<u8>> {
        let mut writer = bgzf::io::Writer::new(Vec::new());

        for chunk in data.chunks(block_len) {
            writer.write_all(chunk)?;
            writer.flush()?;
        }

        writer.finish()
    }

    #[test]
    fn test_copy_blocks_without_eof() -> anyhow::Result<()> {
        let data: Vec<u8> = (0..200_000u32).map(|i| (i % 251) as u8).collect();
        let src = bgzf_blocks(&data, 10_000)?;
        let body = &src[..src.len() - BGZF_EOF.len()];

        for (input, expected_eof) in [(&src[..], true), (body, false)] {
            let mut dst = Vec::new();
            let has_eof = copy_blocks_without_eof(
                &mut Trickle(&input[..100]).chain(&input[100..]),
                &mut dst,
            )?;
            assert_eq!(has_eof, expected_eof);
            assert_eq!(dst, body);
        }

        // Only a final EOF marker is left out.
        let mut doubled = src.clone();
        doubled.extend_from_slice(&src);
        let mut dst = Vec::new();
        assert!(copy_blocks_without_eof(&mut doubled.as_slice(), &mut dst)?);
        assert_eq!(dst, doubled[..doubled.len() - BGZF_EOF.len()]);

        Ok(())
    }

    #[test]
    fn test_copy_to_fails_on_truncated_block() -> anyhow::Result<()> {
        let header: sam::Header = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:1000\n".parse()?;
        let mut file = bam_file(&header, b"records")?;
        // Pseudorandom, so that the blocks hardly compress.
        let mut state = 1u32;
        let records: Vec<u8> = (0..100_000)
            .map(|_| {
                state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                (state >> 16) as u8
            })
            .collect();
        file.truncate(file.len() - BGZF_EOF.len());
        file.extend_from_slice(&bgzf_blocks(&records, 30_000)?);

        // Cut inside the last block of records, and inside the EOF marker's block header.
        for cut in [file.len() - 2000, file.len() - BGZF_EOF.len() + 10] {
            let (_, body) = split(&file[..cut], Path::new("trunc.bam"))?;
            let e = body.copy_to(&mut Vec::new()).unwrap_err();
            let message = format!("{e:#}");
            assert!(message.contains("trunc.bam"), "{message}");
            assert!(message.contains("truncated BGZF block"), "{message}");
        }

        // Cut between blocks: copied, with a warning.
        let (_, body) = split(&file[..file.len() - BGZF_EOF.len()], Path::new("in.bam"))?;
        body.copy_to(&mut Vec::new())?;

        // Not BGZF after the header.
        let mut garbage = file.clone();
        garbage.truncate(garbage.len() - BGZF_EOF.len());
        garbage.extend_from_slice(&[0; 40]);
        let (_, body) = split(garbage.as_slice(), Path::new("in.bam"))?;
        assert!(body.copy_to(&mut Vec::new()).is_err());

        Ok(())
    }
}
