//! Sorted runs of raw BAM records: in-memory chunks, spill files, and k-way merging.
//!
//! Every merge breaks ties by source order, and sources are numbered in input order, so sorting
//! is stable, as in samtools.

use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    fs::File,
    io::{self, BufReader, Read, Write},
    path::{Path, PathBuf},
    thread,
};

use pasta::bgzf::io::writer::CompressionLevel;

use crate::{
    bgzf_stream::{BgzfReader, BgzfWriter},
    sort_order::{Order, compare_same_hash, coordinate_key, name_hash, tag_offset},
};

/// The fixed-size prefix of every encoded BAM record, up to the read name.
const MIN_RECORD_LEN: usize = 32;

/// The largest valid `block_size`. HTSlib reads it as an `int32_t`.
const MAX_RECORD_LEN: usize = i32::MAX as usize;

/// The most a record buffer grows before the bytes to fill it have been read.
const MAX_RECORD_STEP: usize = 1 << 20;

/// Reads one `block_size`-prefixed raw record into `buf`. Returns `false` at end of stream.
pub fn read_raw_record<R: Read>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<bool> {
    let mut len = [0; 4];

    match reader.read(&mut len[..1])? {
        0 => return Ok(false),
        _ => reader.read_exact(&mut len[1..])?,
    }

    let len = usize::try_from(u32::from_le_bytes(len))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // A BGZF stream of a few megabytes can decompress to gigabytes of zeros, so the bytes present
    // do not bound a corrupt `block_size`.
    if !(MIN_RECORD_LEN..=MAX_RECORD_LEN).contains(&len) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid BAM record block size: {len} (expected {MIN_RECORD_LEN} to \
                 {MAX_RECORD_LEN})"
            ),
        ));
    }

    // The buffer grows with the bytes read, and the fixed fields of a record larger than one step
    // are checked before the rest is read.
    if len <= MAX_RECORD_STEP {
        buf.resize(len, 0);
        reader.read_exact(buf)?;
        check_fixed_fields(buf, len)?;
    } else {
        buf.clear();
        buf.resize(MIN_RECORD_LEN, 0);
        reader.read_exact(buf)?;
        check_fixed_fields(buf, len)?;

        while buf.len() < len {
            // Up to the next multiple of `MAX_RECORD_STEP`.
            let start = buf.len();
            let end = len.min((start / MAX_RECORD_STEP + 1) * MAX_RECORD_STEP);
            buf.resize(end, 0);
            reader.read_exact(&mut buf[start..])?;
        }
    }

    Ok(true)
}

/// Checks that a record's fixed fields, at the start of `record`, describe a read name with its
/// NUL and variable-length fields that fit in `len` bytes, as HTSlib does. Keys are read from
/// fixed offsets and the name, so those must be in bounds.
fn check_fixed_fields(record: &[u8], len: usize) -> io::Result<()> {
    let l_read_name = record[8];
    let n_cigar_op = u16::from_le_bytes([record[12], record[13]]);
    let l_seq = u32::from_le_bytes([record[16], record[17], record[18], record[19]]);

    // At most 32 + 255 + 4 * (2^16 - 1) + 3 * 2^31, which does not overflow.
    let need = MIN_RECORD_LEN as u64
        + u64::from(l_read_name)
        + 4 * u64::from(n_cigar_op)
        + u64::from(l_seq).div_ceil(2)
        + u64::from(l_seq);

    if l_read_name == 0 {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM record: l_read_name is 0 (expected at least 1)",
        ))
    } else if need > len as u64 {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid BAM record: l_read_name = {l_read_name}, n_cigar_op = {n_cigar_op} and \
                 l_seq = {l_seq} need {need} bytes, more than the block size ({len})"
            ),
        ))
    } else {
        Ok(())
    }
}

pub fn write_raw_record<W: Write>(writer: &mut W, record: &[u8]) -> io::Result<()> {
    let len =
        u32::try_from(record.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_all(&len.to_le_bytes())?;
    writer.write_all(record)
}

/// [`Entry::key`] for records without the sort tag of [`Order::Tag`].
const NO_TAG: u64 = u64::MAX;

/// The precomputed integer key of an order, if it has one: the coordinate key, the name hash, or
/// the offset of the sort tag's value.
fn integer_key(order: Order, record: &[u8]) -> u64 {
    match order {
        Order::Coordinate => coordinate_key(record),
        Order::Collate => name_hash(record),
        Order::Tag { tag, .. } => tag_offset(record, tag)
            .and_then(|offset| u64::try_from(offset).ok())
            .unwrap_or(NO_TAG),
        Order::NameNatural | Order::NameAscii | Order::TemplateCoordinate => 0,
    }
}

/// Compares records `a` and `b` by `order`, given their [`integer_key`]s, which decide most
/// comparisons without reading the records.
fn compare_keyed(order: Order, a_key: u64, a: &[u8], b_key: u64, b: &[u8]) -> Ordering {
    match order {
        Order::Coordinate => a_key.cmp(&b_key),
        // The key is a name hash; equal hashes fall back to comparing names.
        Order::Collate => a_key.cmp(&b_key).then_with(|| compare_same_hash(a, b)),
        Order::Tag { .. } => order.compare_tagged(a, tag_at(a, a_key), b, tag_at(b, b_key)),
        _ => order.compare(a, b),
    }
}

/// The sort tag value at the offset that `key` holds.
fn tag_at(record: &[u8], key: u64) -> Option<&[u8]> {
    usize::try_from(key)
        .ok()
        .filter(|_| key != NO_TAG)
        .map(|offset| &record[offset..])
}

#[derive(Clone, Copy, Debug)]
struct Entry {
    /// The order's integer key (see [`integer_key`]), so most comparisons compare integers.
    key: u64,
    start: usize,
    end: usize,
}

/// An in-memory batch of raw records in one contiguous buffer.
#[derive(Debug)]
pub struct Chunk {
    data: Vec<u8>,
    entries: Vec<Entry>,
    order: Order,
}

impl Chunk {
    pub fn new(order: Order) -> Self {
        Self {
            data: Vec::new(),
            entries: Vec::new(),
            order,
        }
    }

    pub fn push(&mut self, record: &[u8]) {
        let start = self.data.len();
        self.data.extend_from_slice(record);
        self.entries.push(Entry {
            key: integer_key(self.order, record),
            start,
            end: self.data.len(),
        });
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Approximate memory use, for comparison against the sort memory limit.
    pub fn size(&self) -> usize {
        self.data.len() + self.entries.len() * size_of::<Entry>()
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.entries.clear();
    }

    /// Sorts the chunk with `threads` additional threads and calls `f` with each record in order.
    pub fn drain_sorted<F>(&mut self, order: Order, threads: usize, mut f: F) -> io::Result<()>
    where
        F: FnMut(&[u8]) -> io::Result<()>,
    {
        let data = &self.data;
        let parts = threads + 1;
        let part_len = self.entries.len().div_ceil(parts).max(1);

        let compare = |a: &Entry, b: &Entry| {
            compare_keyed(
                order,
                a.key,
                &data[a.start..a.end],
                b.key,
                &data[b.start..b.end],
            )
        };

        // Stable sorts of contiguous parts, merged with ties going to the earlier part, give a
        // stable sort of the whole chunk.
        thread::scope(|scope| {
            for part in self.entries.chunks_mut(part_len) {
                scope.spawn(move || part.sort_by(compare));
            }
        });

        let parts: Vec<&[Entry]> = self.entries.chunks(part_len).collect();
        let mut heap = BinaryHeap::with_capacity(parts.len());

        for (source, part) in parts.iter().enumerate() {
            if let Some(&entry) = part.first() {
                heap.push(HeapEntry {
                    entry,
                    source,
                    index: 0,
                    compare: &compare,
                });
            }
        }

        while let Some(top) = heap.pop() {
            f(&data[top.entry.start..top.entry.end])?;

            let index = top.index + 1;

            if let Some(&entry) = parts[top.source].get(index) {
                heap.push(HeapEntry {
                    entry,
                    index,
                    ..top
                });
            }
        }

        self.clear();

        Ok(())
    }
}

struct HeapEntry<'c, C> {
    entry: Entry,
    source: usize,
    index: usize,
    compare: &'c C,
}

impl<C: Fn(&Entry, &Entry) -> Ordering> Ord for HeapEntry<'_, C> {
    fn cmp(&self, other: &Self) -> Ordering {
        // `BinaryHeap` is a max-heap, so the smallest record must compare greatest.
        (self.compare)(&self.entry, &other.entry)
            .then(self.source.cmp(&other.source))
            .reverse()
    }
}

impl<C: Fn(&Entry, &Entry) -> Ordering> PartialOrd for HeapEntry<'_, C> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<C: Fn(&Entry, &Entry) -> Ordering> PartialEq for HeapEntry<'_, C> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<C: Fn(&Entry, &Entry) -> Ordering> Eq for HeapEntry<'_, C> {}

/// Writes a chunk as a sorted spill file of raw records, and returns its path.
pub fn spill(
    chunk: &mut Chunk,
    order: Order,
    dir: &Path,
    index: usize,
    threads: usize,
) -> io::Result<PathBuf> {
    let path = dir.join(format!("{index:04}.bam"));
    let file = File::create(&path)?;
    let mut writer = BgzfWriter::new(file, CompressionLevel::FAST, threads);

    chunk.drain_sorted(order, threads, |record| {
        write_raw_record(&mut writer, record)
    })?;

    writer.finish()?;

    Ok(path)
}

/// A source of sorted raw records for a k-way merge.
pub trait RecordSource {
    /// Reads the next record into `buf`, returning `false` when exhausted.
    fn next_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool>;
}

/// A spill file written by [`spill`].
pub struct SpillFile(BgzfReader<BufReader<ReleasingFile>>);

impl SpillFile {
    /// Opens a spill file to read once, from start to end. The disk space of what has been read
    /// is released as reading goes on (see [`ReleasingFile`]), so the file cannot be read again.
    pub fn open(path: &Path, threads: usize) -> io::Result<Self> {
        let file = ReleasingFile::open(path, RELEASE_STEP)?;
        Ok(Self(BgzfReader::new(BufReader::new(file), threads)))
    }
}

/// How many bytes a [`SpillFile`] reads between releasing their disk space.
const RELEASE_STEP: u64 = 16 << 20;

/// A file read once, from start to end, that releases the disk space of what it has read.
///
/// A merge's temporary files then shrink as the merge goes, so a pipeline such as
/// `collate | sort` needs temporary space for about one copy of the reads, not one per command.
/// samtools gets the same by removing each `collate` bucket once it is read. Where releasing is
/// unsupported, this is a plain file.
struct ReleasingFile {
    file: File,
    /// Bytes read.
    read: u64,
    /// Bytes released, from the start of the file.
    released: u64,
    /// Bytes to read between releases.
    step: u64,
    /// Cleared if releasing fails, after which the file is only read.
    releasing: bool,
}

impl ReleasingFile {
    fn open(path: &Path, step: u64) -> io::Result<Self> {
        // Punching holes needs write access.
        let file = File::options().read(true).write(true).open(path)?;

        Ok(Self {
            file,
            read: 0,
            released: 0,
            step,
            releasing: cfg!(target_os = "linux"),
        })
    }
}

impl Read for ReleasingFile {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.file.read(buf)?;
        self.read += n as u64;

        if self.releasing && self.read - self.released >= self.step {
            self.releasing = punch_hole(&self.file, self.released, self.read - self.released);
            self.released = self.read;
        }

        Ok(n)
    }
}

/// Frees the disk space of `len` bytes at `offset`, keeping the file's size, and returns whether
/// it could.
#[cfg(target_os = "linux")]
fn punch_hole(file: &File, offset: u64, len: u64) -> bool {
    use std::os::fd::AsRawFd;

    let (Ok(offset), Ok(len)) = (libc::off_t::try_from(offset), libc::off_t::try_from(len)) else {
        return false;
    };

    // SAFETY: `fallocate` only changes the file behind the descriptor, which `file` keeps open.
    unsafe {
        libc::fallocate(
            file.as_raw_fd(),
            libc::FALLOC_FL_PUNCH_HOLE | libc::FALLOC_FL_KEEP_SIZE,
            offset,
            len,
        ) == 0
    }
}

#[cfg(not(target_os = "linux"))]
fn punch_hole(_: &File, _: u64, _: u64) -> bool {
    false
}

impl RecordSource for SpillFile {
    fn next_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool> {
        read_raw_record(&mut self.0, buf)
    }
}

struct Head {
    /// The record's [`integer_key`], computed once rather than in every comparison.
    key: u64,
    record: Vec<u8>,
    source: usize,
    order: Order,
}

impl Ord for Head {
    fn cmp(&self, other: &Self) -> Ordering {
        compare_keyed(self.order, self.key, &self.record, other.key, &other.record)
            .then(self.source.cmp(&other.source))
            .reverse()
    }
}

impl PartialOrd for Head {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Head {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Head {}

/// Merges sorted sources, calling `f` with each record in order. Ties go to the earlier source.
pub fn merge<S, F>(sources: &mut [S], order: Order, mut f: F) -> anyhow::Result<()>
where
    S: RecordSource,
    F: FnMut(&[u8]) -> anyhow::Result<()>,
{
    let mut heap = BinaryHeap::with_capacity(sources.len());

    for (source, s) in sources.iter_mut().enumerate() {
        let mut record = Vec::new();

        if s.next_record(&mut record)? {
            heap.push(Head {
                key: integer_key(order, &record),
                record,
                source,
                order,
            });
        }
    }

    while let Some(mut head) = heap.pop() {
        f(&head.record)?;

        if sources[head.source].next_record(&mut head.record)? {
            head.key = integer_key(order, &head.record);
            heap.push(head);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_releasing_file_releases_what_it_has_read() -> io::Result<()> {
        use std::{fs, os::unix::fs::MetadataExt};

        let dir = tempfile::tempdir()?;
        let path = dir.path().join("0000.bam");
        let data: Vec<u8> = (0..1u32 << 20).map(|i| (i % 251) as u8).collect();
        fs::write(&path, &data)?;
        let blocks = fs::metadata(&path)?.blocks();

        let mut file = ReleasingFile::open(&path, 64 << 10)?;
        let mut read = Vec::new();
        let mut buf = [0; 8192];

        loop {
            match file.read(&mut buf)? {
                0 => break,
                n => read.extend_from_slice(&buf[..n]),
            }
        }

        assert_eq!(read, data);

        // Only Linux punches holes, and not on every filesystem; elsewhere the file is only read.
        if file.releasing {
            let metadata = fs::metadata(&path)?;
            assert_eq!(metadata.len(), data.len() as u64);
            assert!(
                metadata.blocks() < blocks / 8,
                "{} of {blocks} blocks left",
                metadata.blocks()
            );
        }

        Ok(())
    }

    #[test]
    fn test_read_raw_record() -> io::Result<()> {
        let mut record = vec![0; MIN_RECORD_LEN];
        record[8] = 2; // l_read_name
        record.extend(b"*\x00");

        let mut src = Vec::new();
        write_raw_record(&mut src, &record)?;
        write_raw_record(&mut src, &record)?;

        // A record larger than one read step.
        let mut large = record.clone();
        large.resize(3 * MAX_RECORD_STEP + 5, 0x01);
        write_raw_record(&mut src, &large)?;

        let mut reader = &src[..];
        let mut buf = Vec::new();
        assert!(read_raw_record(&mut reader, &mut buf)?);
        assert_eq!(buf, record);
        assert!(read_raw_record(&mut reader, &mut buf)?);
        assert_eq!(buf, record);
        assert!(read_raw_record(&mut reader, &mut buf)?);
        assert_eq!(buf, large);
        assert!(!read_raw_record(&mut reader, &mut buf)?);

        Ok(())
    }

    #[test]
    fn test_read_raw_record_with_invalid_block_size() {
        let mut record = vec![0; MIN_RECORD_LEN];
        record[8] = 2; // l_read_name
        record.extend(b"*\x00");

        // Below the fixed-length fields, 0 included.
        for block_size in [0u32, 31] {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&record);
            let mut buf = Vec::new();

            assert!(matches!(
                read_raw_record(&mut &src[..], &mut buf),
                Err(e) if e.kind() == io::ErrorKind::InvalidData
            ));
        }

        // Larger than the rest of the stream: memory grows with the bytes read, not with
        // `block_size`.
        for block_size in [35, 1 << 24, i32::MAX.unsigned_abs()] {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&record);
            let mut buf = Vec::new();

            assert!(matches!(
                read_raw_record(&mut &src[..], &mut buf),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
            ));
            assert!(buf.capacity() <= MAX_RECORD_STEP, "{}", buf.capacity());
        }
    }

    #[test]
    fn test_read_raw_record_with_invalid_fixed_fields() {
        fn fixed_fields(l_read_name: u8, n_cigar_op: u16, l_seq: u32) -> Vec<u8> {
            let mut src = vec![0; MIN_RECORD_LEN];
            src[8] = l_read_name;
            src[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
            src[16..20].copy_from_slice(&l_seq.to_le_bytes());
            src
        }

        // As in htslib: a `block_size` above `i32::MAX`, a read name without its NUL, or lengths
        // that need more than `block_size` are errors, found from the block size and the fixed
        // fields alone, before the rest of the record is read.
        let cases = [
            (0x8000_0000, fixed_fields(2, 0, 0)),
            (u32::MAX, fixed_fields(2, 0, 0)),
            (34, fixed_fields(0, 0, 0)),
            (i32::MAX.unsigned_abs(), fixed_fields(0, 0, 0)),
            (34, fixed_fields(3, 0, 0)),
            (100, fixed_fields(2, u16::MAX, 0)),
            (i32::MAX.unsigned_abs(), fixed_fields(2, 0, u32::MAX)),
            (i32::MAX.unsigned_abs(), fixed_fields(2, 0, 1 << 31)),
        ];

        for (block_size, fields) in cases {
            let mut src = block_size.to_le_bytes().to_vec();
            src.extend(&fields);
            // Zeros after the fixed fields, which a reader that reads the whole block before it
            // checks them would read.
            let mut reader = (&src[..]).chain(io::repeat(0).take(64 << 20));
            let mut buf = Vec::new();

            let result = read_raw_record(&mut reader, &mut buf);
            assert!(
                matches!(&result, Err(e) if e.kind() == io::ErrorKind::InvalidData),
                "{block_size}: {result:?}"
            );

            if usize::try_from(block_size).unwrap() > MAX_RECORD_STEP {
                assert!(buf.capacity() <= MAX_RECORD_STEP, "{}", buf.capacity());
                assert_eq!(reader.into_inner().1.limit(), 64 << 20, "{block_size}");
            }
        }
    }
}
