//! Alignment inputs (SAM, BAM, CRAM): headers, streaming, and index queries, yielding each
//! record as its BAM encoding as HTSlib's `sam_read1` and `sam_itr_next` do.

use std::{
    collections::{HashMap, VecDeque},
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use pasta::{
    bam, bgzf,
    core::{Position, Region, region::Interval},
    cram::{self, io::reader::MdNm},
    csi::{self, BinningIndex as _},
    fasta, sam,
    variant::index::split_idx_name,
};

use crate::{faidx::Faidx, raw_input::strerror};

/// The input buffer size.
const BUFFER_CAPACITY: usize = 1 << 17;

/// An alignment file format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Sam,
    Bam,
    Cram,
    /// No data, or compressed data that decompresses to nothing (HTSlib's `empty_format`,
    /// which has no header to read).
    Empty,
}

/// File descriptors that readers of the inputs use.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Descriptors {
    /// Kept open while reading: the files, and the references CRAM files are decoded with.
    pub held: u64,
    /// Opened for a moment while a reader opens (a CRAM reference's FASTA index). The inputs
    /// are opened one after another, so this is the most any one input opens.
    pub opening: u64,
}

impl std::ops::AddAssign for Descriptors {
    fn add_assign(&mut self, other: Self) {
        self.held += other.held;
        self.opening = self.opening.max(other.opening);
    }
}

/// An input's index, loaded once and shared by the readers of the input.
#[derive(Debug)]
pub enum AlnIndex {
    Bam(bam::Index),
    Cram(cram::crai::Index),
}

/// The interval that [`AlnReader::read_region`] queries for `beg` to `end` (0-based,
/// inclusive) of a sequence of `len` bases, or `None` if it reads nothing.
fn query_interval(beg: u32, end: u32, len: i64) -> anyhow::Result<Option<Interval>> {
    if i64::from(beg) >= len {
        // Nothing to read past the end of the sequence (HTSlib's query is empty).
        return Ok(None);
    }
    // An end past the sequence (as in `chr:1-` or a whole-sequence span) is open, as HTSlib
    // clamps it; indexes cannot represent such large positions.
    let start = Position::try_from(beg as usize + 1)?;
    Ok(Some(if i64::from(end) + 1 >= len {
        if beg == 0 {
            Interval::from(..)
        } else {
            Interval::from(start..)
        }
    } else {
        Interval::from(start..=Position::try_from(end as usize + 1)?)
    }))
}

/// Where an input's index places records, to tell the parts of a sequence that a query of
/// the input reads nothing from (before opening it).
pub enum Occupancy<'a> {
    Bam(&'a bam::Index),
    /// The CRAI's slices on each sequence: whether there is one, as a query of the whole
    /// sequence reads every slice on it, and the merged 1-based extents of those that align
    /// bases, which a query of part of the sequence reads.
    Cram(HashMap<usize, (bool, Vec<(usize, usize)>)>),
}

impl AlnIndex {
    /// Where the index places records.
    pub fn occupancy(&self) -> Occupancy<'_> {
        let index = match self {
            Self::Bam(index) => return Occupancy::Bam(index),
            Self::Cram(index) => index,
        };
        let mut by_tid: HashMap<usize, (bool, Vec<(usize, usize)>)> = HashMap::new();
        for record in index {
            let Some(tid) = record.reference_sequence_id() else {
                continue;
            };
            let (any, extents) = by_tid.entry(tid).or_default();
            *any = true;
            // As `index_record_intersects`, which picks the slices a query reads.
            let span = record.alignment_span();
            if let Some(start) = record.alignment_start()
                && span > 0
            {
                let start = usize::from(start);
                extents.push((start, start.saturating_add(span - 1)));
            }
        }
        for (_, extents) in by_tid.values_mut() {
            extents.sort_unstable();
            let mut merged: Vec<(usize, usize)> = Vec::with_capacity(extents.len());
            for &(s, e) in extents.iter() {
                match merged.last_mut() {
                    Some(last) if s <= last.1.saturating_add(1) => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            *extents = merged;
        }
        Occupancy::Cram(by_tid)
    }
}

impl Occupancy<'_> {
    /// Whether [`AlnReader::read_region`] of `beg` to `end` (0-based, inclusive) of sequence
    /// `tid`, `len` bases long, may read a record: `false` only where it reads none.
    pub fn may_read(&self, tid: usize, beg: u32, end: u32, len: i64) -> bool {
        let Ok(interval) = query_interval(beg, end, len) else {
            return true;
        };
        let Some(interval) = interval else {
            return false;
        };
        match self {
            // A query of the whole sequence reads every bin's chunks (listing the bins it
            // reads takes a while).
            Self::Bam(index) if interval.start().is_none() && interval.end().is_none() => {
                fn any_chunk<I: csi::binning_index::index::reference_sequence::Index>(
                    index: &csi::binning_index::Index<I>,
                    tid: usize,
                ) -> bool {
                    index.reference_sequences().get(tid).is_none_or(|sequence| {
                        sequence.bins().values().any(|bin| !bin.chunks().is_empty())
                    })
                }
                match index {
                    bam::Index::Bai(index) => any_chunk(index, tid),
                    bam::Index::Csi(index) => any_chunk(index, tid),
                }
            }
            // The chunks of the bins that the query reads (an error is left to the query).
            Self::Bam(index) => !index.query(tid, interval).is_ok_and(|c| c.is_empty()),
            Self::Cram(by_tid) => {
                let Some((any, extents)) = by_tid.get(&tid) else {
                    return false;
                };
                let (Some(start), end) = (interval.start(), interval.end()) else {
                    // The whole sequence.
                    return *any;
                };
                let (start, end) = (usize::from(start), end.map_or(usize::MAX, usize::from));
                let i = extents.partition_point(|&(_, e)| e < start);
                extents.get(i).is_some_and(|&(s, _)| s <= end)
            }
        }
    }
}

/// The parts of an alignment header mpileup uses.
#[derive(Clone, Debug)]
pub struct AlnHeader {
    pub sam: sam::Header,
    /// The header text, as `sam_hdr_str` gives it (for the read groups).
    pub text: Vec<u8>,
    pub names: Vec<Vec<u8>>,
    pub lengths: Vec<i64>,
}

impl AlnHeader {
    fn new(sam: sam::Header) -> anyhow::Result<Self> {
        let mut writer = sam::io::Writer::new(Vec::new());
        writer.write_header(&sam)?;
        let text = writer.into_inner();
        let names = sam
            .reference_sequences()
            .keys()
            .map(|n| n.to_vec())
            .collect();
        let lengths = sam
            .reference_sequences()
            .values()
            .map(|rs| rs.length().get() as i64)
            .collect();
        Ok(Self {
            sam,
            text,
            names,
            lengths,
        })
    }

    /// The reference sequence ID of a name.
    pub fn tid(&self, name: &[u8]) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
}

/// A BGZF reader that decompresses on the calling thread or on a worker pool.
enum BgzfReader<R> {
    Single(bgzf::io::Reader<R>),
    Multi(bgzf::io::MultithreadedReader<R>),
}

impl<R: Read + Send + 'static> BgzfReader<R> {
    /// Decompresses on `threads` workers (fewer when the system cannot start that many), or on
    /// the calling thread with none or when it cannot start any.
    fn new(inner: R, threads: usize) -> Self {
        match NonZero::new(threads) {
            Some(n) => match bgzf::io::MultithreadedReader::try_with_worker_count(n, inner) {
                Ok(reader) => Self::Multi(reader),
                Err(e) => Self::Single(bgzf::io::Reader::new(e.into_inner())),
            },
            None => Self::Single(bgzf::io::Reader::new(inner)),
        }
    }
}

impl<R: Read + Send + 'static> Read for BgzfReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Single(r) => r.read(buf),
            Self::Multi(r) => r.read(buf),
        }
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        match self {
            Self::Single(r) => r.read_exact(buf),
            Self::Multi(r) => r.read_exact(buf),
        }
    }
}

impl<R: Read + Send + 'static> BufRead for BgzfReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Single(r) => r.fill_buf(),
            Self::Multi(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::Single(r) => r.consume(amt),
            Self::Multi(r) => r.consume(amt),
        }
    }
}

impl<R: Read + Send + 'static> bgzf::io::Read for BgzfReader<R> {
    fn virtual_position(&self) -> bgzf::VirtualPosition {
        match self {
            Self::Single(r) => r.virtual_position(),
            Self::Multi(r) => r.virtual_position(),
        }
    }
}

impl<R: Read + Send + 'static> bgzf::io::BufRead for BgzfReader<R> {}

impl<R: Read + Seek + Send + 'static> bgzf::io::Seek for BgzfReader<R> {
    fn seek_to_virtual_position(
        &mut self,
        pos: bgzf::VirtualPosition,
    ) -> io::Result<bgzf::VirtualPosition> {
        match self {
            Self::Single(r) => r.seek_to_virtual_position(pos),
            Self::Multi(r) => r.seek_to_virtual_position(pos),
        }
    }

    fn seek_with_index(&mut self, index: &bgzf::gzi::Index, pos: SeekFrom) -> io::Result<u64> {
        match self {
            Self::Single(r) => r.seek_with_index(index, pos),
            Self::Multi(r) => r.seek_with_index(index, pos),
        }
    }
}

/// An input stream: a file, or stdin.
enum Stream {
    File(File),
    Stdin(io::Stdin),
}

impl Read for Stream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::File(f) => f.read(buf),
            Self::Stdin(s) => s.read(buf),
        }
    }
}

impl Seek for Stream {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match self {
            Self::File(f) => f.seek(pos),
            Self::Stdin(_) => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "stdin is not seekable",
            )),
        }
    }
}

/// A prefix of bytes read for format detection, put back in front of the stream.
struct Prefixed {
    prefix: Vec<u8>,
    at: usize,
    inner: Stream,
}

impl Read for Prefixed {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.at < self.prefix.len() {
            let n = buf.len().min(self.prefix.len() - self.at);
            buf[..n].copy_from_slice(&self.prefix[self.at..self.at + n]);
            self.at += n;
            return Ok(n);
        }
        self.inner.read(buf)
    }
}

impl Seek for Prefixed {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.at = self.prefix.len();
        self.inner.seek(pos)
    }
}

type Raw = BufReader<Prefixed>;

/// An open alignment reader.
enum Inner {
    Bam(bam::io::Reader<BgzfReader<Raw>>),
    Cram(cram::io::Reader<Raw>, fasta::Repository),
    Sam(sam::io::Reader<Box<dyn BufRead + Send>>),
}

/// An alignment input opened for reading.
pub struct AlnReader {
    inner: Inner,
    path: Option<PathBuf>,
    /// The index named with `FILE##idx##INDEX`, used instead of the one next to the file.
    index: Option<PathBuf>,
    name: String,
}

/// Opens a FASTA reference repository for CRAM decoding.
pub fn reference_repository(src: Option<&Path>) -> anyhow::Result<fasta::Repository> {
    let Some(src) = src else {
        return Ok(fasta::Repository::default());
    };
    let reader = Faidx::load(src)
        .with_context(|| format!("opening reference {}", src.display()))?
        .into_indexed_reader();
    let repository =
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader));
    // The sequence dictionary beside it (`ref.dict`), whose checksums check the reference for
    // CRAM slices that carry none.
    repository
        .load_dictionary(src)
        .with_context(|| format!("reading the sequence dictionary of {}", src.display()))?;
    Ok(repository)
}

/// The local FASTA file a CRAM header's `@SQ UR` fields name, as HTSlib uses it when no
/// reference is given.
fn header_reference(header: &sam::Header) -> Option<PathBuf> {
    use sam::header::record::value::map::reference_sequence::tag::URI;
    header
        .reference_sequences()
        .values()
        .filter_map(|rs| rs.other_fields().get(&URI))
        .map(|uri| {
            let uri = uri.to_string();
            PathBuf::from(uri.strip_prefix("file://").unwrap_or(&uri))
        })
        .find(|path| {
            let mut index = path.clone().into_os_string();
            index.push(".fai");
            path.is_file() && Path::new(&index).is_file()
        })
}

/// mpileup's error for an input `src` (`data`, perhaps with `##idx##INDEX`) that cannot be
/// opened, after what HTSlib's `sam_open` prints. The messages are part of the error, so that
/// when the reader threads of several inputs fail, only the one reported is printed.
fn open_error(src: &str, data: &str, e: &io::Error) -> anyhow::Error {
    let why = strerror(e);
    let hopen = if e.kind() == io::ErrorKind::IsADirectory {
        format!("[E::hts_hopen] Failed to open file {data}\n")
    } else {
        String::new()
    };
    super::exit1(format!(
        "{hopen}[E::hts_open_format] Failed to open file \"{data}\" : {why}\n\
         [mpileup] failed to open {src}: {why}"
    ))
}

impl AlnReader {
    /// Opens `src` (`-` for stdin), detects its format, and reads its header. `reference` is
    /// the FASTA file for CRAM decoding; `threads` adds BGZF decompression workers.
    /// `FILE##idx##INDEX` opens `FILE` and queries it through `INDEX`, as HTSlib does.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or its header read.
    pub fn open(
        src: &str,
        reference: Option<&Path>,
        threads: usize,
    ) -> anyhow::Result<(Self, AlnHeader)> {
        let (data, index) = split_idx_name(src);
        let (stream, path) = if data == "-" {
            (Stream::Stdin(io::stdin()), None)
        } else {
            let file = File::open(data).map_err(|e| open_error(src, data, &e))?;
            (Stream::File(file), Some(PathBuf::from(data)))
        };
        let mut stream = stream;
        // Enough for one whole BGZF block, to look for the BAM magic inside it.
        let mut prefix = vec![0u8; 1 << 16];
        let mut n = 0;
        while n < prefix.len() {
            match stream.read(&mut prefix[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                // HTSlib reads the start of a file to open it: a directory fails here.
                Err(e) if path.is_some() => return Err(open_error(src, data, &e)),
                Err(e) => return Err(e.into()),
            }
        }
        prefix.truncate(n);
        let format = detect(&prefix);
        let raw = BufReader::with_capacity(
            BUFFER_CAPACITY,
            Prefixed {
                prefix,
                at: 0,
                inner: stream,
            },
        );
        let read_header = || format!("[mpileup] fail to read the header of {src}");
        let (inner, header) = match format {
            Format::Empty => bail!("{}", read_header()),
            Format::Bam => {
                let mut r = bam::io::Reader::from(BgzfReader::new(raw, threads));
                let h = r.read_header().with_context(read_header)?;
                (Inner::Bam(r), h)
            }
            Format::Cram => {
                let repo = reference_repository(reference)?;
                let mut r = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo.clone())
                    .build_from_reader(raw);
                r.set_md_nm(MdNm::Regenerate);
                let h = r.read_header().with_context(read_header)?;
                let repo = if reference.is_none()
                    && let Some(path) = header_reference(&h)
                {
                    let repo = reference_repository(Some(&path))?;
                    r.set_reference_sequence_repository(repo.clone());
                    repo
                } else {
                    repo
                };
                (Inner::Cram(r, repo), h)
            }
            Format::Sam => {
                let mut buf = raw;
                let is_gz = buf.fill_buf()?.starts_with(&[0x1f, 0x8b]);
                let mut inner: Box<dyn BufRead + Send> = if is_gz {
                    Box::new(BgzfReader::new(buf, threads))
                } else {
                    Box::new(buf)
                };
                // HTSlib reads an empty input as `empty_format`, which has no header.
                if inner.fill_buf()?.is_empty() {
                    bail!(read_header());
                }
                let mut r = sam::io::Reader::new(inner);
                let h = r.read_header().with_context(read_header)?;
                (Inner::Sam(r), h)
            }
        };
        let header = AlnHeader::new(header)?;
        Ok((
            Self {
                inner,
                path,
                index: index.map(PathBuf::from),
                name: src.to_owned(),
            },
            header,
        ))
    }

    /// Calls `f` with every record's BAM encoding in file order, from where `cursor` is, until
    /// `f` returns `false`: the next call carries on after that record. Returns whether the
    /// input was read to its end.
    ///
    /// # Errors
    ///
    /// Returns an error if a record cannot be read, or `f`'s error.
    pub fn read_all<F>(
        &mut self,
        header: &AlnHeader,
        cursor: &mut Cursor,
        mut f: F,
    ) -> anyhow::Result<bool>
    where
        F: FnMut(Vec<u8>) -> anyhow::Result<bool>,
    {
        match &mut self.inner {
            Inner::Bam(reader) => {
                let mut record = bam::Record::default();
                while reader.read_record(&mut record)? != 0 {
                    let raw = Vec::from(std::mem::take(&mut record));
                    if !f(raw)? {
                        return Ok(false);
                    }
                }
            }
            Inner::Cram(reader, repository) => {
                let CramAll {
                    container,
                    pending,
                    error,
                } = &mut **cursor.cram.get_or_insert_with(Box::default);
                // The records of the container read last that `f` has not had, and what
                // failed after them.
                while let Some(raw) = pending.pop_front() {
                    if !f(raw)? {
                        return Ok(false);
                    }
                }
                if let Some(e) = error.take() {
                    return Err(e);
                }
                while reader.read_container(container)? != 0 {
                    // Once `f` pauses, the rest of the container is kept for the next call.
                    let mut paused = false;
                    let result = read_container_records(container, repository, header, |raw| {
                        if paused {
                            pending.push_back(raw);
                        } else {
                            paused = !f(raw)?;
                        }
                        Ok(())
                    });
                    match result {
                        Ok(()) => {}
                        // A record after the pause cannot be read: it fails after the records
                        // before it, at the next call.
                        Err(e) if paused => *error = Some(e),
                        Err(e) => return Err(e),
                    }
                    if paused {
                        return Ok(false);
                    }
                }
            }
            Inner::Sam(reader) => {
                let mut record = sam::alignment::RecordBuf::default();
                let mut encoder = bam::io::Writer::from(Vec::new());
                while reader.read_record_buf(&header.sam, &mut record)? != 0 {
                    if !f(encode(&mut encoder, &header.sam, &record)?)? {
                        return Ok(false);
                    }
                }
            }
        }
        Ok(true)
    }

    /// Calls `f` with the records overlapping `[beg, end]` (0-based, inclusive) of reference
    /// sequence `tid`, as `sam_itr_querys` selects them through `index` (the input's index,
    /// from [`Self::load_index`]), from where `cursor` is, until `f` returns `false`: the next
    /// call carries on after that record. Returns whether the region was read to its end.
    ///
    /// # Errors
    ///
    /// Returns an error if there is no index for the input or a record cannot be read.
    #[expect(
        clippy::too_many_arguments,
        reason = "the region and where reading it stopped"
    )]
    pub fn read_region<F>(
        &mut self,
        header: &AlnHeader,
        index: Option<&AlnIndex>,
        tid: usize,
        beg: u32,
        end: u32,
        cursor: &mut Cursor,
        mut f: F,
    ) -> anyhow::Result<bool>
    where
        F: FnMut(Vec<u8>) -> anyhow::Result<bool>,
    {
        let fail = || format!("[mpileup] fail to load index for {}", self.name);
        let name = std::str::from_utf8(&header.names[tid]).context("invalid sequence name")?;
        let len = header.lengths.get(tid).copied().unwrap_or(0);
        let Some(interval) = query_interval(beg, end, len)? else {
            return Ok(true);
        };
        let region = Region::new(name, interval);
        match (&mut self.inner, index) {
            (Inner::Bam(reader), Some(AlnIndex::Bam(index))) => {
                let mut query = match cursor.query.take() {
                    Some(SuspendedQuery::Bam(suspended)) => {
                        bam::io::reader::Query::resume(reader, suspended)
                    }
                    _ => reader.query(&header.sam, index, &region)?,
                };
                let mut record = bam::Record::default();
                while query.read_record(&mut record)? != 0 {
                    let raw = Vec::from(std::mem::take(&mut record));
                    if !f(raw)? {
                        cursor.query = Some(SuspendedQuery::Bam(query.suspend()));
                        return Ok(false);
                    }
                }
            }
            (Inner::Cram(reader, _), Some(AlnIndex::Cram(index))) => {
                let mut query = match cursor.query.take() {
                    Some(SuspendedQuery::Cram(suspended)) => {
                        cram::io::reader::Query::resume(reader, &header.sam, index, *suspended)
                    }
                    _ => reader.query(&header.sam, index, &region)?,
                };
                let mut record = sam::alignment::RecordBuf::default();
                let mut encoder = bam::io::Writer::from(Vec::new());
                while query.read_record_buf(&mut record)? != 0 {
                    if !f(encode(&mut encoder, &header.sam, &record)?)? {
                        cursor.query = Some(SuspendedQuery::Cram(Box::new(query.suspend())));
                        return Ok(false);
                    }
                }
            }
            _ => bail!("{}", fail()),
        }
        Ok(true)
    }

    /// The file descriptors another reader of this input would use.
    pub fn descriptors(&self) -> Descriptors {
        match self.inner {
            Inner::Cram(..) => Descriptors {
                held: 2,
                opening: 1,
            },
            Inner::Bam(_) | Inner::Sam(_) => Descriptors {
                held: 1,
                opening: 0,
            },
        }
    }

    /// Loads the input's index for region queries (`None` if it has none that loads). The
    /// index file is closed once read.
    pub fn load_index(&self) -> Option<AlnIndex> {
        let path = self.path.as_deref()?;
        let index = self.index.as_deref();
        match &self.inner {
            Inner::Bam(_) => read_bam_index(path, index).ok().map(AlnIndex::Bam),
            Inner::Cram(..) => read_crai(path, index).ok().map(AlnIndex::Cram),
            Inner::Sam(_) => None,
        }
    }
}

/// Reads the BAI or CSI index `index`, told apart by content as HTSlib does (a BAI starts
/// with its magic number, a CSI is BGZF-compressed), or else the one next to `path`.
fn read_bam_index(path: &Path, index: Option<&Path>) -> io::Result<bam::Index> {
    const BAI_MAGIC: &[u8; 4] = b"BAI\x01";
    let Some(index) = index else {
        return bam::fs::read_associated_index(path);
    };
    let mut magic = [0; 4];
    File::open(index)?.read_exact(&mut magic)?;
    if &magic == BAI_MAGIC {
        bam::bai::fs::read(index).map(bam::Index::Bai)
    } else {
        pasta::csi::fs::read(index).map(bam::Index::Csi)
    }
}

/// Reads the CRAI index `index`, or else the one next to `path`.
fn read_crai(path: &Path, index: Option<&Path>) -> io::Result<cram::crai::Index> {
    match index {
        Some(index) => cram::crai::fs::read(index),
        None => cram::fs::read_associated_index(path),
    }
}

/// Where [`AlnReader::read_all`] or [`AlnReader::read_region`] stopped, to carry on from: a
/// new cursor starts at the beginning (of the region).
#[derive(Default)]
pub struct Cursor {
    /// The region query, set aside.
    query: Option<SuspendedQuery>,
    /// Reading a whole CRAM file.
    cram: Option<Box<CramAll>>,
}

enum SuspendedQuery {
    Bam(bam::io::reader::SuspendedQuery),
    Cram(Box<cram::io::reader::SuspendedQuery>),
}

/// Where reading a whole CRAM file stopped: the records of the container read last that were
/// not yet passed on, and the error that followed them.
#[derive(Default)]
struct CramAll {
    container: cram::io::reader::Container,
    pending: VecDeque<Vec<u8>>,
    error: Option<anyhow::Error>,
}

/// Calls `f` with the BAM encoding of every record of `container`.
fn read_container_records(
    container: &cram::io::reader::Container,
    repository: &fasta::Repository,
    header: &AlnHeader,
    mut f: impl FnMut(Vec<u8>) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let compression_header = container.compression_header()?;
    for result in container.slices() {
        let slice = result?;
        let (core, external) = slice.decode_blocks()?;
        let records = slice.records_with_md_nm(
            repository.clone(),
            &header.sam,
            &compression_header,
            &core,
            &external,
            MdNm::Regenerate,
        )?;
        for record in &records {
            let mut buf = Vec::new();
            record.encode_bam_record(&header.sam, &mut buf)?;
            f(buf)?;
        }
    }
    Ok(())
}

/// Encodes a record as BAM, without its `block_size`.
fn encode(
    encoder: &mut bam::io::Writer<Vec<u8>>,
    header: &sam::Header,
    record: &sam::alignment::RecordBuf,
) -> io::Result<Vec<u8>> {
    use sam::alignment::io::Write as _;
    encoder.get_mut().clear();
    encoder.write_alignment_record(header, record)?;
    Ok(encoder.get_ref()[4..].to_vec())
}

/// Detects the format from the first bytes, as `hts_detect_format` does for the formats
/// mpileup reads.
fn detect(prefix: &[u8]) -> Format {
    if prefix.is_empty() {
        return Format::Empty;
    }
    if prefix.starts_with(b"CRAM") {
        return Format::Cram;
    }
    if prefix.starts_with(&[0x1f, 0x8b]) {
        let mut decoder = bgzf::io::Reader::new(prefix);
        let mut magic = [0u8; 4];
        let (mut n, mut failed) = (0, false);
        while n < magic.len() {
            match decoder.read(&mut magic[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(_) => {
                    failed = true;
                    break;
                }
            }
        }
        return match &magic[..n] {
            b"BAM\x01" => Format::Bam,
            [] if !failed => Format::Empty,
            _ => Format::Sam,
        };
    }
    Format::Sam
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The BGZF end-of-file block: a gzip member holding no data.
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    #[test]
    fn detects_empty_inputs_as_htslib_does() {
        assert_eq!(detect(b""), Format::Empty);
        assert_eq!(detect(&BGZF_EOF), Format::Empty);
        assert_eq!(detect(b"@HD\tVN:1.6\n"), Format::Sam);
        assert_eq!(detect(b"CRAM\x03\x01"), Format::Cram);
        // A cut gzip member is not empty: it fails later, as a SAM file.
        assert_eq!(detect(&BGZF_EOF[..20]), Format::Sam);
    }

    #[test]
    fn tells_where_a_crai_query_reads_nothing() {
        let slice = |tid, start: usize, span| {
            cram::crai::Record::new(tid, Position::new(start), span, 0, 0, 0)
        };
        let index = AlnIndex::Cram(vec![
            // 1-based 1001 to 2000 and 1500 to 3000 of sequence 0.
            slice(Some(0), 1001, 1000),
            slice(Some(0), 1500, 1501),
            // Sequence 1 has a slice aligning nothing; a query of all of it reads the slice.
            slice(Some(1), 0, 0),
            slice(None, 0, 0),
        ]);
        let occupancy = index.occupancy();
        let len = 1_000_000;
        let may_read = |tid, beg, end| occupancy.may_read(tid, beg, end, len);
        assert!(may_read(0, 0, u32::MAX), "all of sequence 0");
        assert!(may_read(0, 0, 1000), "0-based 1000 is 1-based 1001");
        assert!(!may_read(0, 0, 999));
        assert!(may_read(0, 2999, 5000));
        assert!(!may_read(0, 3000, 5000));
        assert!(!may_read(0, 3000, u32::MAX), "the rest of sequence 0");
        assert!(may_read(1, 0, u32::MAX), "all of sequence 1");
        assert!(!may_read(1, 0, 999_998), "part of sequence 1");
        assert!(!may_read(2, 0, u32::MAX), "no slice on sequence 2");
        assert!(
            !may_read(0, 1_000_000, u32::MAX),
            "past the end of the sequence"
        );
    }
}
