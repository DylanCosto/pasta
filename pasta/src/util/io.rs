//! Reading the start of a stream to detect its format, then reading the stream from its start.

use std::io::{self, Read, Seek, SeekFrom};

use crate::bgzf::BGZF_HEADER_SIZE;

const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

/// How much of a gzip stream that is not BGZF to read: enough to decompress its first bytes.
const GZIP_PREFIX_LEN: usize = 1 << 16;

/// Reads the start of a stream for format detection: the whole first BGZF block, the start of
/// other gzip data, or a few bytes of anything else, stopping early at the end of the stream.
///
/// A pipe can return fewer bytes per read than detection needs (say, a BGZF header without its
/// compressed data), so this reads until it has enough.
pub fn read_detection_prefix<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let mut buf = Vec::with_capacity(BGZF_HEADER_SIZE);
    fill(reader, &mut buf, BGZF_HEADER_SIZE)?;

    if buf.starts_with(&GZIP_MAGIC_NUMBER) {
        let len = bgzf_block_size(&buf).unwrap_or(GZIP_PREFIX_LEN);
        fill(reader, &mut buf, len)?;
    }

    Ok(buf)
}

/// Returns the size of the BGZF block whose header starts `src`, if it is one.
fn bgzf_block_size(src: &[u8]) -> Option<usize> {
    const BGZF_SUBFIELD: [u8; 4] = [b'B', b'C', 0x02, 0x00];

    let header = src.get(..BGZF_HEADER_SIZE)?;

    if header[12..16] != BGZF_SUBFIELD {
        return None;
    }

    let bsize = u16::from_le_bytes([header[16], header[17]]);
    Some(usize::from(bsize) + 1)
}

/// Reads from `reader` until `buf` holds `len` bytes or the stream ends.
fn fill<R>(reader: &mut R, buf: &mut Vec<u8>, len: usize) -> io::Result<()>
where
    R: Read,
{
    if buf.len() >= len {
        return Ok(());
    }

    let missing = u64::try_from(len - buf.len()).map_err(io::Error::other)?;
    reader.by_ref().take(missing).read_to_end(buf)?;

    Ok(())
}

/// A reader whose first bytes were already read (for format detection) and are read again
/// first.
pub struct Rewound<R> {
    prefix: Vec<u8>,
    position: usize,
    inner: R,
}

impl<R> Rewound<R> {
    /// Creates a reader that reads `prefix`, then `inner`, which is positioned right after it.
    pub fn new(prefix: Vec<u8>, inner: R) -> Self {
        Self {
            prefix,
            position: 0,
            inner,
        }
    }

    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    fn remaining_prefix(&self) -> &[u8] {
        &self.prefix[self.position..]
    }
}

impl<R> Read for Rewound<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let prefix = self.remaining_prefix();

        if prefix.is_empty() {
            return self.inner.read(buf);
        }

        let n = prefix.len().min(buf.len());
        buf[..n].copy_from_slice(&prefix[..n]);
        self.position += n;

        Ok(n)
    }
}

impl<R> Seek for Rewound<R>
where
    R: Seek,
{
    /// Seeks the underlying reader; the prefix is not read again after a seek. Positions are
    /// those of the underlying reader, which starts at the start of the prefix.
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let pos = match pos {
            // The underlying reader is ahead of this one by the unread part of the prefix.
            SeekFrom::Current(n) => {
                let unread =
                    i64::try_from(self.remaining_prefix().len()).map_err(io::Error::other)?;
                SeekFrom::Current(n - unread)
            }
            pos => pos,
        };

        let position = self.inner.seek(pos)?;
        self.position = self.prefix.len();

        Ok(position)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    /// A reader that returns at most `n` bytes per read, like a pipe.
    struct Trickle<'a>(&'a [u8], usize);

    impl Read for Trickle<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.0.len().min(buf.len()).min(self.1);
            buf[..n].copy_from_slice(&self.0[..n]);
            self.0 = &self.0[n..];
            Ok(n)
        }
    }

    #[test]
    fn test_read_detection_prefix() -> io::Result<()> {
        use std::io::Write;

        let mut writer = crate::bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"BAM\x01")?;
        let src = writer.finish()?;
        let first_block_len = bgzf_block_size(&src).expect("a BGZF block");

        // The whole first block, even from a reader that returns one byte at a time.
        let prefix = read_detection_prefix(&mut Trickle(&src, 1))?;
        assert_eq!(prefix, &src[..first_block_len]);

        // Short streams end early.
        assert_eq!(read_detection_prefix(&mut Trickle(b"CRAM", 1))?, b"CRAM");
        assert!(read_detection_prefix(&mut Trickle(b"", 1))?.is_empty());

        Ok(())
    }

    #[test]
    fn test_rewound() -> io::Result<()> {
        let mut src = Cursor::new(b"abcdef".to_vec());
        let prefix = read_detection_prefix(&mut src)?;
        let mut reader = Rewound::new(prefix, src);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;
        assert_eq!(buf, b"abcdef");

        let mut src = Cursor::new(b"abcdef".to_vec());
        let mut prefix = vec![0; 3];
        src.read_exact(&mut prefix)?;
        let mut reader = Rewound::new(prefix, src);
        let mut buf = [0; 1];
        reader.read_exact(&mut buf)?;
        assert_eq!(reader.stream_position()?, 1);
        assert_eq!(reader.seek(SeekFrom::Current(2))?, 3);
        buf = [0];
        reader.read_exact(&mut buf)?;
        assert_eq!(&buf, b"d");

        Ok(())
    }
}
