use std::{
    io::{self, Read},
    mem,
    ops::Range,
};

/// The length of the fixed-length fields that start every record, the smallest valid
/// `block_size`.
const MIN_BLOCK_SIZE: usize = 32;

/// The largest valid `block_size`. HTSlib reads it as an `int32_t`.
const MAX_BLOCK_SIZE: usize = i32::MAX as usize;

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

    // A BGZF stream of a few megabytes can decompress to gigabytes of zeros, so the bytes present
    // do not bound a corrupt `block_size`.
    if block_size > MAX_BLOCK_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid record block size: {block_size} (expected at most {MAX_BLOCK_SIZE})"),
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

    Ok(block_size)
}

/// Reads a record of `len` bytes into `buf` and checks its fixed-length fields.
///
/// `len` comes from the stream, so `buf` grows with the bytes read, `MAX_STEP` at a time, and the
/// fixed-length fields of a record longer than that are checked before the rest is read, as in
/// HTSlib.
fn read_block<R>(reader: &mut R, buf: &mut Vec<u8>, len: usize) -> io::Result<()>
where
    R: Read,
{
    if len <= MAX_STEP {
        buf.resize(len, 0);
        reader.read_exact(buf)?;
        return check_fixed_fields(buf, len);
    }

    buf.clear();
    buf.resize(MIN_BLOCK_SIZE, 0);
    reader.read_exact(buf)?;
    check_fixed_fields(buf, len)?;

    while buf.len() < len {
        // Up to the next multiple of `MAX_STEP`.
        let start = buf.len();
        let end = len.min((start / MAX_STEP + 1) * MAX_STEP);
        buf.resize(end, 0);
        reader.read_exact(&mut buf[start..])?;
    }

    Ok(())
}

/// Checks that a record's fixed-length fields, at the start of `src`, describe a read name with
/// its NUL and variable-length fields that fit in `block_size` bytes, as HTSlib does.
fn check_fixed_fields(src: &[u8], block_size: usize) -> io::Result<()> {
    const NAME_LENGTH_INDEX: usize = 8;
    const CIGAR_OP_COUNT_RANGE: Range<usize> = 12..14;
    const READ_LENGTH_RANGE: Range<usize> = 16..20;

    let name_len = src[NAME_LENGTH_INDEX];

    if name_len == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid record: l_read_name is 0 (expected at least 1)",
        ));
    }

    // SAFETY: The ranges have the lengths of the integers.
    let cigar_op_count = u16::from_le_bytes(src[CIGAR_OP_COUNT_RANGE].try_into().unwrap());
    let base_count = u32::from_le_bytes(src[READ_LENGTH_RANGE].try_into().unwrap());

    // At most 32 + 255 + 4 * (2^16 - 1) + 3 * 2^31, which does not overflow.
    let len = MIN_BLOCK_SIZE as u64
        + u64::from(name_len)
        + 4 * u64::from(cigar_op_count)
        + u64::from(base_count).div_ceil(2)
        + u64::from(base_count);

    if len > block_size as u64 {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid record: l_read_name = {name_len}, n_cigar_op = {cigar_op_count} and \
                 l_seq = {base_count} need {len} bytes, more than the block size ({block_size})"
            ),
        ))
    } else {
        Ok(())
    }
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
        for block_size in [35, 1 << 24, i32::MAX.unsigned_abs()] {
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

    /// A reader of `src`, then `pad` zero bytes, that counts the bytes read.
    struct CountingReader<'a> {
        inner: io::Chain<&'a [u8], io::Take<io::Repeat>>,
        bytes_read: u64,
    }

    impl<'a> CountingReader<'a> {
        fn new(src: &'a [u8], pad: u64) -> Self {
            Self {
                inner: src.chain(io::repeat(0).take(pad)),
                bytes_read: 0,
            }
        }
    }

    impl Read for CountingReader<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = self.inner.read(buf)?;
            self.bytes_read += n as u64;
            Ok(n)
        }
    }

    /// Enough zero bytes after a record's fixed fields that a reader that reads the whole block
    /// before it checks them reads far more than the fixed fields.
    const PAD: u64 = 64 << 20;

    #[test]
    fn test_read_record_with_block_size_above_i32_max() {
        // As in htslib, which reads `block_size` as an `int32_t`, a `block_size` above
        // `i32::MAX` is an error, found before the record is read.
        for block_size in [0x8000_0000, u32::MAX] {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&unmapped_record()[4..36]);
            let mut reader = CountingReader::new(&src, PAD);
            let mut buf = Vec::new();

            assert!(matches!(
                read_record(&mut reader, &mut buf),
                Err(e) if e.kind() == io::ErrorKind::InvalidData
            ));
            assert_eq!(reader.bytes_read, 4);
        }
    }

    #[test]
    fn test_read_record_with_fixed_fields_beyond_block_size() {
        fn fixed_fields(l_read_name: u8, n_cigar_op: u16, l_seq: u32) -> Vec<u8> {
            let mut src = unmapped_record()[4..36].to_vec();
            src[8] = l_read_name;
            src[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
            src[16..20].copy_from_slice(&l_seq.to_le_bytes());
            src
        }

        // As in htslib, the lengths in the fixed fields are checked against `block_size` before
        // the rest of the record is read, and a read name must have at least its NUL.
        let cases = [
            // l_read_name = 0
            (34, fixed_fields(0, 0, 0)),
            (i32::MAX.unsigned_abs(), fixed_fields(0, 0, 0)),
            // 32 + 3 > 34
            (34, fixed_fields(3, 0, 0)),
            // 32 + 2 + 4 * 65535 > 100
            (100, fixed_fields(2, u16::MAX, 0)),
            // l_seq = u32::MAX; htslib reads it as -1
            (i32::MAX.unsigned_abs(), fixed_fields(2, 0, u32::MAX)),
            // 32 + 2 + 2^30 + 2^31 > i32::MAX
            (i32::MAX.unsigned_abs(), fixed_fields(2, 0, 1 << 31)),
        ];

        for (block_size, fields) in cases {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&fields);
            let mut reader = CountingReader::new(&src, PAD);
            let mut buf = Vec::new();

            let result = read_record(&mut reader, &mut buf);
            assert!(
                matches!(&result, Err(e) if e.kind() == io::ErrorKind::InvalidData),
                "{block_size}: {result:?}"
            );

            // Only a block that fits in one read step is read before it is checked.
            if usize::try_from(block_size).unwrap() > MAX_STEP {
                assert_eq!(reader.bytes_read, 4 + 32, "{block_size}");
            }
        }

        // Fixed fields that fill the block exactly are valid.
        let mut src = 35u32.to_le_bytes().to_vec();
        src.extend(fixed_fields(1, 0, 1));
        src.extend([0x00, 0x10, 0xff]); // read_name = "\x00", seq = "A", qual = 0xff
        let mut reader = &src[..];
        let mut buf = Vec::new();
        assert!(matches!(read_record(&mut reader, &mut buf), Ok(35)));
    }
}
