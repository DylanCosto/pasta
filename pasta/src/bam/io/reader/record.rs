use std::{
    io::{self, Read},
    mem,
    ops::Range,
};

/// The length of the fixed-length fields that start every record, the smallest valid
/// `block_size`.
const MIN_BLOCK_SIZE: usize = 32;

/// The most a record buffer grows before the bytes to fill it have been read.
const MAX_STEP: usize = 1 << 20;

pub(super) fn read_record<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: Read,
{
    // Only a stream that ends at a record boundary is EOF. A `block_size` of 0 is not.
    let Some(block_size) = read_block_size(reader)? else {
        return Ok(0);
    };

    if block_size < MIN_BLOCK_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid record block size: {block_size} (expected at least {MIN_BLOCK_SIZE})"),
        ));
    }

    read_block(reader, buf, block_size).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("truncated record: expected {block_size} bytes (block size)"),
            )
        } else {
            e
        }
    })?;

    validate(buf)?;

    Ok(block_size)
}

/// Reads `len` bytes into `buf`.
///
/// `len` comes from the stream, and a corrupt one can be up to 4 GiB, so `buf` grows with the
/// bytes read, `MAX_STEP` at a time.
fn read_block<R>(reader: &mut R, buf: &mut Vec<u8>, len: usize) -> io::Result<()>
where
    R: Read,
{
    if len <= MAX_STEP {
        buf.resize(len, 0);
        return reader.read_exact(buf);
    }

    buf.clear();

    while buf.len() < len {
        let start = buf.len();
        buf.resize(start + (len - start).min(MAX_STEP), 0);
        reader.read_exact(&mut buf[start..])?;
    }

    Ok(())
}

/// Reads `block_size`, or returns `None` if the stream is at EOF.
fn read_block_size<R>(reader: &mut R) -> io::Result<Option<usize>>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];

    if !read_exact_or_eof(reader, &mut buf)? {
        return Ok(None);
    }

    let n = u32::from_le_bytes(buf);
    usize::try_from(n)
        .map(Some)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

/// Fills `buf`, or returns `false` if the stream is at EOF.
fn read_exact_or_eof<R>(reader: &mut R, mut buf: &mut [u8]) -> io::Result<bool>
where
    R: Read,
{
    let mut bytes_read = 0;

    while !buf.is_empty() {
        match reader.read(buf) {
            Ok(0) => break,
            Ok(n) => {
                buf = &mut buf[n..];
                bytes_read += n;
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    if bytes_read == 0 {
        Ok(false)
    } else if buf.is_empty() {
        Ok(true)
    } else {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated record block size",
        ))
    }
}

pub(crate) fn validate(src: &[u8]) -> io::Result<()> {
    const MIN_BUF_LENGTH: usize = 32;
    const NAME_LENGTH_INDEX: usize = 8;
    const CIGAR_OP_COUNT_RANGE: Range<usize> = 12..14;
    const READ_LENGTH_RANGE: Range<usize> = 16..20;

    if src.len() < MIN_BUF_LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let name_len = usize::from(src[NAME_LENGTH_INDEX]);

    let buf = &src[CIGAR_OP_COUNT_RANGE];
    // SAFETY: `buf.len() == mem::size_of::<u16>()`.
    let cigar_op_count = usize::from(u16::from_le_bytes(buf.try_into().unwrap()));

    let buf = &src[READ_LENGTH_RANGE];
    // SAFETY: `buf.len() == mem::size_of::<u32>()`.
    let base_count = usize::try_from(u32::from_le_bytes(buf.try_into().unwrap()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let quality_scores_end = MIN_BUF_LENGTH
        + name_len
        + (cigar_op_count * mem::size_of::<u32>())
        + base_count.div_ceil(2)
        + base_count;

    if src.len() < quality_scores_end {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_block_size() -> io::Result<()> {
        let data = [0x08, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader)?, Some(8));

        let data = [0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader)?, Some(0));

        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_block_size(&mut reader)?, None);

        let data = [0x08];
        let mut reader = &data[..];
        assert!(matches!(
            read_block_size(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }

    /// An unmapped record named `*`, `block_size`-prefixed.
    fn unmapped_record() -> Vec<u8> {
        let mut src = vec![0x22, 0x00, 0x00, 0x00]; // block_size = 34
        src.extend([
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ]);
        src
    }

    #[test]
    fn test_read_record() -> io::Result<()> {
        let record = unmapped_record();
        let mut src = record.repeat(2);
        let mut reader = &src[..];
        let mut buf = Vec::new();

        assert_eq!(read_record(&mut reader, &mut buf)?, 34);
        assert_eq!(buf, record[4..]);
        assert_eq!(read_record(&mut reader, &mut buf)?, 34);
        assert_eq!(buf, record[4..]);

        // Only a stream that ends at a record boundary is EOF.
        assert_eq!(read_record(&mut reader, &mut buf)?, 0);

        // A record larger than one read step.
        src.clear();
        let len = 3 * MAX_STEP + 5;
        src.extend(u32::try_from(len).unwrap().to_le_bytes());
        src.extend(&record[4..36]);
        src[4 + 8] = 1; // l_read_name = 1
        src.resize(4 + len, 0x01);
        src[4 + 32] = 0x00; // read_name = "\x00"
        let mut reader = &src[..];
        assert_eq!(read_record(&mut reader, &mut buf)?, len);
        assert_eq!(buf, src[4..]);
        assert_eq!(read_record(&mut reader, &mut buf)?, 0);

        Ok(())
    }

    #[test]
    fn test_read_record_with_block_size_below_fixed_fields() {
        // As in htslib, a `block_size` below the 32 bytes of fixed-length fields is an error, not
        // EOF, 0 included.
        for block_size in [0u32, 1, 31] {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&unmapped_record()[4..]);
            let mut reader = &src[..];
            let mut buf = Vec::new();

            assert!(matches!(
                read_record(&mut reader, &mut buf),
                Err(e) if e.kind() == io::ErrorKind::InvalidData
            ));
        }
    }

    #[test]
    fn test_read_record_with_truncated_record() {
        // A `block_size` larger than the rest of the stream.
        for block_size in [35, 1 << 24, u32::MAX] {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&unmapped_record()[4..]);
            let mut reader = &src[..];
            let mut buf = Vec::new();

            assert!(matches!(
                read_record(&mut reader, &mut buf),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
            ));

            // Memory grows with the bytes read, not with `block_size`.
            assert!(buf.capacity() <= MAX_STEP, "{}", buf.capacity());
        }
    }
}
