//! CRAM writer.

pub(crate) mod builder;
mod collections;
pub mod compression;
pub(crate) mod container;
mod context;
pub(crate) mod header;
pub(crate) mod num;
mod parallel;
pub(crate) mod record;

use std::{
    collections::HashSet,
    io::{self, Write},
    mem,
    ops::Range,
    sync::Arc,
    thread,
};

use rustc_hash::FxHashMap;

use crate::bam;
use crate::fasta;
use crate::sam::{self as sam, alignment::io::Write as _};

pub use self::builder::Builder;
use self::{
    compression::Compressor,
    container::write_container,
    header::{write_file_definition, write_file_header, write_header},
    parallel::Workers,
};
pub(crate) use self::{
    context::Context,
    record::{ReadGroups, Record, ReferenceBases},
};
use crate::cram::FileDefinition;

pub(crate) const DEFAULT_SLICES_PER_CONTAINER: usize = 1;

/// Returns the reference sequence named `name`, or `None` if the repository lacks it. Missing names
/// are remembered, since the repository only caches sequences it has.
pub(crate) fn get_reference_sequence(
    repository: &fasta::Repository,
    missing_names: &mut HashSet<Vec<u8>>,
    name: &[u8],
) -> io::Result<Option<Arc<fasta::record::Sequence>>> {
    if missing_names.contains(name) {
        return Ok(None);
    }

    let sequence = repository.get(name).transpose()?;

    if sequence.is_none() {
        missing_names.insert(name.to_vec());
    }

    Ok(sequence)
}

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    use crate::sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};

    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .map(|sort_order| sort_order == COORDINATE)
        .unwrap_or_default()
}

/// A file header whose missing reference sequence checksums are calculated on another thread
/// (see [`Writer::write_header`]), and the containers written meanwhile, in order.
#[derive(Debug)]
struct PendingHeader {
    header: thread::JoinHandle<io::Result<sam::Header>>,
    buf: Vec<u8>,
}

/// The most bytes of containers held while the file header is pending: past this, writing waits
/// for the header.
const MAX_PENDING_LEN: usize = 256 << 20;

/// The reference sequences a writer used (see [`Writer::use_reference_sequence`]).
#[derive(Debug, Default)]
struct ReferenceSequences {
    /// The reference sequence of the last record on one.
    last: Option<usize>,
    /// The reference sequences released from the repository's cache.
    released: HashSet<usize>,
    /// Whether to keep every reference sequence in the cache.
    keeps_all: bool,
    /// The reference sequences checked against their `@SQ` `M5` (see
    /// [`validate_reference_sequence_checksum`]).
    validated: HashSet<usize>,
}

/// Checks that the reference sequence `name`, which records are stored against, matches its
/// `@SQ` `M5` field, if any, as HTSlib's `validate_md5` does: a CRAM file whose header checksum
/// is not that of the reference its records are stored against cannot be read back (and a reader
/// that finds the reference by its checksum gets another sequence).
fn validate_reference_sequence_checksum(
    repository: &fasta::Repository,
    name: &[u8],
    reference_sequence: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> io::Result<()> {
    use crate::sam::header::record::value::map::reference_sequence::{Md5Checksum, tag};
    use bstr::ByteSlice;

    let Some(value) = reference_sequence.other_fields().get(&tag::MD5_CHECKSUM) else {
        return Ok(());
    };

    let Some(actual) = repository.normalized_md5(name).transpose()? else {
        return Ok(());
    };

    // HTSlib compares the hexadecimal strings, so a malformed checksum never matches.
    let matches = value
        .to_str()
        .ok()
        .and_then(|s| s.parse::<Md5Checksum>().ok())
        .is_some_and(|expected| *expected == actual);

    if matches {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "SQ header M5 tag discrepancy for reference '{}': the header (@SQ M5) has {}, \
                 but the reference has {}; please use the correct reference, or consider using \
                 ignore_md5=1",
                name.as_bstr(),
                value.as_bstr(),
                Md5Checksum::from(actual),
            ),
        ))
    }
}

/// Records waiting for their container.
#[derive(Debug, Default)]
struct Batch {
    /// The records, in order.
    records: Vec<Pending>,
    /// The records converted already.
    converted: Vec<Record>,
    /// The data of the BAM records, one after another. Copying them here, rather than each into
    /// its own allocation, keeps the thread writing records from contending for the allocator
    /// with the workers that free them.
    bam_records: Vec<u8>,
    /// The reference sequences of the records, by ID, kept loaded until the container is
    /// encoded.
    reference_sequences: FxHashMap<usize, Option<Arc<fasta::record::Sequence>>>,
    /// The reference sequence ID of the last record, once there is one.
    last_reference_sequence_id: Option<Option<usize>>,
    /// Whether the records are on more than one reference sequence (or unmapped and not).
    has_many_reference_sequences: bool,
}

/// A record waiting for its container.
#[derive(Clone, Debug)]
enum Pending {
    /// The next record of [`Batch::converted`].
    Converted,
    /// A BAM record, converted on a worker: its data in [`Batch::bam_records`].
    Bam(Range<usize>),
}

impl Batch {
    /// A batch of converted records.
    fn from_converted(records: Vec<Record>) -> Self {
        let mut batch = Self {
            records: vec![Pending::Converted; records.len()],
            bam_records: Vec::new(),
            reference_sequences: FxHashMap::default(),
            ..Default::default()
        };

        for record in &records {
            batch.add_reference_sequence_id(record.reference_sequence_id);
        }

        batch.converted = records;
        batch
    }

    /// An empty batch with room for `records` records, likely as many of them converted and with
    /// as much BAM data as this batch.
    fn next(&self, records: usize) -> Self {
        Self {
            records: Vec::with_capacity(records),
            converted: Vec::with_capacity(self.converted.len()),
            bam_records: Vec::with_capacity(self.bam_records.len()),
            reference_sequences: FxHashMap::default(),
            last_reference_sequence_id: None,
            has_many_reference_sequences: false,
        }
    }

    fn len(&self) -> usize {
        self.records.len()
    }

    fn push(&mut self, record: Record) {
        self.add_reference_sequence_id(record.reference_sequence_id);
        self.records.push(Pending::Converted);
        self.converted.push(record);
    }

    /// Adds a BAM record on `reference_sequence_id`.
    fn push_bam_record(&mut self, record: &bam::Record, reference_sequence_id: Option<usize>) {
        self.add_reference_sequence_id(reference_sequence_id);
        let start = self.bam_records.len();
        self.bam_records.extend_from_slice(record.as_ref());
        let end = self.bam_records.len();
        self.records.push(Pending::Bam(start..end));
    }

    fn add_reference_sequence_id(&mut self, reference_sequence_id: Option<usize>) {
        if self
            .last_reference_sequence_id
            .is_some_and(|last| last != reference_sequence_id)
        {
            self.has_many_reference_sequences = true;
        }

        self.last_reference_sequence_id = Some(reference_sequence_id);
    }
}

/// Converts the records of a container.
/// Converts the records of `batch` not converted yet, on the batch's reference sequences, which
/// the caller holds (`reference_sequences`) while the container is encoded.
fn convert(
    header: &sam::Header,
    batch: Batch,
    reference_sequences: &FxHashMap<usize, Option<Arc<fasta::record::Sequence>>>,
) -> io::Result<Vec<Record>> {
    if batch.converted.len() == batch.records.len() {
        return Ok(batch.converted);
    }

    // The records' reference sequences are borrowed rather than each record cloning (and dropping)
    // one: workers counting references to the same sequence would contend for its cache line.
    let get_reference_sequence = |name: &[u8]| {
        let sequence = header
            .reference_sequences()
            .get_index_of(name)
            .and_then(|id| reference_sequences.get(&id))
            .and_then(|sequence| sequence.as_deref());

        Ok(sequence)
    };

    let mut converted = batch.converted.into_iter();
    let mut bam_record = bam::Record::default();
    let mut read_groups = ReadGroups::default();

    batch
        .records
        .into_iter()
        .map(|record| match record {
            Pending::Converted => Ok(converted.next().expect("a converted record is pending")),
            Pending::Bam(range) => {
                bam_record.0.clear();
                bam_record.0.extend_from_slice(&batch.bam_records[range]);
                Record::try_from_bam_record(
                    header,
                    &bam_record,
                    &mut read_groups,
                    get_reference_sequence,
                )
            }
        })
        .collect()
}

/// A CRAM writer.
///
/// [`Self::try_finish`] writes the pending container and the EOF container and reports any
/// error. A writer whose header was written is finished when it is dropped, as BAM and BGZF
/// writers are, but errors are then lost, so call [`Self::try_finish`] (or
/// [`sam::alignment::io::Write::finish`]). A writer dropped while its thread panics is left
/// unfinished, without the EOF container that marks a complete file.
///
/// As in HTSlib, writing a record on a reference sequence whose `@SQ` `M5` field (in the header
/// passed with the record) is not the checksum of the repository's sequence fails, since the file
/// would not decode with the reference its records are stored against, unless
/// [`compression::Options::ignore_md5`] is set. Each reference sequence is checked once.
///
/// With workers ([`Builder::set_worker_count`]), containers are encoded and compressed on them,
/// and the records of a container are written with the header passed with the call that
/// completes it. A failure is then reported by a later call than the one that wrote the
/// container's records, and every call after a failure fails.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use pasta::cram;
/// use pasta::sam::{self as sam, alignment::io::Write};
///
/// let mut writer = cram::io::Writer::new(io::sink());
///
/// let header = sam::Header::default();
/// writer.write_header(&header)?;
///
/// let record = sam::Record::default();
/// writer.write_alignment_record(&header, &record)?;
///
/// writer.try_finish(&header)?;
/// # Ok::<(), io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<W>
where
    W: Write,
{
    /// The output, until [`Self::into_inner`] takes it.
    inner: Option<W>,
    /// The header of [`Self::write_header`] or [`Self::write_file_header`], while the writer is
    /// to be finished on drop.
    unfinished_header: Option<sam::Header>,
    context: Context,
    /// Chooses each block's compression method, with [`Context::compression`].
    compressor: Option<Compressor>,
    batch: Batch,
    records_per_container: usize,
    /// The bases of the pending records.
    bases: usize,
    record_counter: u64,
    /// Reference sequences the repository lacks.
    missing_reference_sequence_names: HashSet<Vec<u8>>,
    /// The file header while its reference sequence checksums are calculated (see
    /// [`Self::write_header`]).
    pending_header: Option<PendingHeader>,
    /// The reference sequence of the last record converted on this thread, by name, as
    /// consecutive records are usually on the same one.
    last_reference_sequence: Option<(Vec<u8>, Option<Arc<fasta::record::Sequence>>)>,
    /// The read group of the last record converted on this thread.
    read_groups: ReadGroups,
    /// The reference sequences used so far.
    reference_sequences: ReferenceSequences,
    /// The workers that encode containers, if any; without them, this thread does.
    workers: Option<Workers>,
    /// Whether a container failed on the workers.
    failed: bool,
    /// Whether records on another reference sequence join the pending container (making
    /// slices with many reference sequences) rather than start a new one. Tests set it to make
    /// such slices.
    pub(crate) many_references: bool,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let writer = cram::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner
            .as_ref()
            .expect("the writer has its inner writer")
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        self.inner
            .as_mut()
            .expect("the writer has its inner writer")
    }

    /// Returns the underlying writer, without finishing the stream (see [`Self::try_finish`]).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(mut self) -> W {
        self.unfinished_header = None;
        // A file header still being completed is written, as far as it can be.
        let _ = self.write_pending_header(true);
        self.inner.take().expect("the writer has its inner writer")
    }
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a new CRAM writer with default options.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let writer = cram::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_from_writer(inner)
    }

    /// Attempts to finish the output stream by writing any pending containers and a final EOF
    /// container.
    ///
    /// A writer is also finished when it is dropped, but errors are then ignored. After this
    /// call, dropping the writer writes nothing more.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// use pasta::sam;
    ///
    /// let header = sam::Header::default();
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self, header: &sam::Header) -> io::Result<()> {
        use self::container::write_eof_container;

        // Whether or not this succeeds, a drop does not finish again.
        self.unfinished_header = None;

        self.flush(header)?;

        if self.workers.is_some() {
            self.advance_workers(0)?;
        }

        self.write_pending_header(true)?;

        write_eof_container(self.get_mut())
    }

    /// Writes a CRAM file definition.
    ///
    /// The file ID is set as a blank value (`[0x00; 20]`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// writer.write_file_definition()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_definition(&mut self) -> io::Result<()> {
        let file_definition = FileDefinition::new(self.context.version, Default::default());
        write_file_definition(self.get_mut(), &file_definition)
    }

    /// Writes a CRAM file header container.
    ///
    /// The position of the stream is expected to be directly after the file definition.
    ///
    /// Entries in the reference sequence dictionary that are missing MD5 checksums (`M5`) will
    /// automatically be calculated and added to the written record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// use pasta::sam;
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    /// writer.write_file_definition()?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_file_header(&header)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_header(&mut self, header: &sam::Header) -> io::Result<()> {
        let threads = self.workers.as_ref().map_or(1, Workers::thread_count);

        write_file_header(
            self.inner
                .as_mut()
                .expect("the writer has its inner writer"),
            &self.context.reference_sequence_repository,
            header,
            threads,
        )?;

        self.unfinished_header = Some(header.clone());

        Ok(())
    }

    /// Writes a SAM header.
    ///
    /// This writes the CRAM magic number, the file definition, and file header using the given SAM
    /// header.
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// use pasta::sam;
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        use crate::sam::header::record::value::map::reference_sequence::tag;

        let file_definition = FileDefinition::new(self.context.version, Default::default());
        let inner = self
            .inner
            .as_mut()
            .expect("the writer has its inner writer");

        let lacks_checksums = header
            .reference_sequences()
            .values()
            .any(|map| !map.other_fields().contains_key(&tag::MD5_CHECKSUM));

        match &self.workers {
            // The checksums, which for a large reference take long to calculate, are calculated on
            // another thread while the workers encode the first containers.
            Some(workers) if lacks_checksums => {
                write_file_definition(inner, &file_definition)?;

                let repository = self.context.reference_sequence_repository.clone();
                let threads = workers.thread_count();
                let mut pending = header.clone();
                let releases_sequences = is_coordinate_sorted(header);

                let handle =
                    thread::Builder::new()
                        .name("cram-checksums".into())
                        .spawn(move || {
                            header::add_missing_reference_sequence_checksums(
                                &repository,
                                pending.reference_sequences_mut(),
                                releases_sequences,
                                threads,
                            )?;

                            Ok(pending)
                        })?;

                self.pending_header = Some(PendingHeader {
                    header: handle,
                    buf: Vec::new(),
                });
            }
            _ => {
                let threads = self.workers.as_ref().map_or(1, Workers::thread_count);
                write_header(
                    inner,
                    &self.context.reference_sequence_repository,
                    &file_definition,
                    header,
                    threads,
                )?;
            }
        }

        self.unfinished_header = Some(header.clone());

        Ok(())
    }

    /// Writes the pending file header (see [`Self::write_header`]) and the containers written
    /// since, once its checksums are calculated, or waiting for them if `wait`.
    fn write_pending_header(&mut self, wait: bool) -> io::Result<()> {
        let Some(pending) = &self.pending_header else {
            return Ok(());
        };

        if !wait && !pending.header.is_finished() && pending.buf.len() < MAX_PENDING_LEN {
            return Ok(());
        }

        let pending = self.pending_header.take().expect("a header is pending");

        let header = pending
            .header
            .join()
            .map_err(|_| io::Error::other("the reference checksum thread panicked"))??;

        let inner = self
            .inner
            .as_mut()
            .expect("the writer has its inner writer");

        header::write_header_container(inner, &header)?;
        inner.write_all(&pending.buf)
    }

    /// Writes the containers the workers have finished, keeping at most `max` in flight.
    fn advance_workers(&mut self, max: usize) -> io::Result<()> {
        let workers = self.workers.as_mut().expect("the writer has workers");

        let result = match &mut self.pending_header {
            Some(pending) => workers.advance(&mut pending.buf, self.compressor.as_mut(), max),
            None => workers.advance(
                self.inner
                    .as_mut()
                    .expect("the writer has its inner writer"),
                self.compressor.as_mut(),
                max,
            ),
        };

        let result = result.and_then(|()| self.write_pending_header(false));
        self.check(result)
    }

    /// Writes a CRAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// use pasta::sam::{self as sam, alignment::io::Write};
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = cram::Record::default();
    /// writer.write_record(&header, &record)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(
        &mut self,
        header: &sam::Header,
        record: &crate::cram::Record<'_>,
    ) -> io::Result<()> {
        self.write_alignment_record(header, record)
    }

    /// Writes a BAM record.
    ///
    /// This writes what [`sam::alignment::io::Write::write_alignment_record`] writes, but with
    /// workers ([`Builder::set_worker_count`]), the record is converted on them.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::{bam, cram, sam};
    ///
    /// let mut writer = cram::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header)?;
    ///
    /// let record = bam::Record::default();
    /// writer.write_bam_record(&header, &record)?;
    ///
    /// writer.try_finish(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_bam_record(
        &mut self,
        header: &sam::Header,
        record: &bam::Record,
    ) -> io::Result<()> {
        if self.workers.is_none() {
            if self.failed {
                return Err(failed());
            }

            let repository = &self.context.reference_sequence_repository;
            let missing_names = &mut self.missing_reference_sequence_names;
            let last = &mut self.last_reference_sequence;

            // The last sequence is borrowed, not cloned, for each record.
            let record =
                Record::try_from_bam_record(header, record, &mut self.read_groups, |name| {
                    if last
                        .as_ref()
                        .is_none_or(|(last_name, _)| last_name[..] != name[..])
                    {
                        let sequence = get_reference_sequence(repository, missing_names, name)?;
                        *last = Some((name.to_vec(), sequence));
                    }

                    let last: &_ = last;
                    Ok(last.as_ref().and_then(|(_, sequence)| sequence.as_deref()))
                })?;

            let read_length = record.read_length;
            self.start_record(header, record.reference_sequence_id)?;
            self.use_reference_sequence(header, record.reference_sequence_id)?;
            self.batch.push(record);
            return self.record_added(header, read_length);
        }

        if self.failed {
            return Err(failed());
        }

        let read_length = record::read_length(record)?;
        let reference_sequence_id = record.reference_sequence_id().transpose()?;
        self.start_record(header, reference_sequence_id)?;
        self.use_reference_sequence(header, reference_sequence_id)?;
        self.batch.push_bam_record(record, reference_sequence_id);
        self.record_added(header, read_length)
    }

    /// Keeps the reference sequence of a record on `reference_sequence_id` loaded until its
    /// container is encoded.
    ///
    /// When coordinate-sorted records (`SO:coordinate`) move on to another reference sequence, the
    /// last one is released from the repository's cache, so they do not keep every sequence
    /// loaded: like HTSlib, which frees a reference sequence when no container uses it. A
    /// released sequence that comes back means the records are not sorted after all, and then
    /// none are released.
    fn use_reference_sequence(
        &mut self,
        header: &sam::Header,
        reference_sequence_id: Option<usize>,
    ) -> io::Result<()> {
        let Some(id) = reference_sequence_id else {
            return Ok(());
        };

        if self.batch.reference_sequences.contains_key(&id) {
            return Ok(());
        }

        let Some((name, reference_sequence)) = header.reference_sequences().get_index(id) else {
            // Converting the record fails.
            return Ok(());
        };

        let repository = &self.context.reference_sequence_repository;
        let sequence =
            get_reference_sequence(repository, &mut self.missing_reference_sequence_names, name)?;

        let used = &mut self.reference_sequences;

        if sequence.is_some() && !self.context.ignores_md5() && !used.validated.contains(&id) {
            validate_reference_sequence_checksum(repository, name, reference_sequence)?;
            used.validated.insert(id);
        }

        self.batch.reference_sequences.insert(id, sequence);

        if used.last == Some(id) {
            return Ok(());
        }

        if used.released.contains(&id) || !is_coordinate_sorted(header) {
            used.keeps_all = true;
        }

        if let Some(last) = used.last
            && !used.keeps_all
            && let Some((last_name, _)) = header.reference_sequences().get_index(last)
        {
            repository.release(last_name);
            used.released.insert(last);
        }

        used.last = Some(id);

        Ok(())
    }

    /// Ends the pending container before a record on `reference_sequence_id` if the record is on
    /// another reference sequence.
    ///
    /// As in HTSlib, a record on another reference sequence starts a new container, so that the
    /// slices of sorted data each have one reference sequence (and a reference checksum). Records
    /// join the pending container instead, making slices with many reference sequences, while the
    /// reference sequence changes after few records (e.g., unsorted data, or many small contigs).
    ///
    /// A record on a reference sequence the header lacks fails its container's conversion, so it
    /// ends nothing: on workers as without them, the containers before it are the ones written.
    fn start_record(
        &mut self,
        header: &sam::Header,
        reference_sequence_id: Option<usize>,
    ) -> io::Result<()> {
        let in_header =
            reference_sequence_id.is_none_or(|id| id < header.reference_sequences().len());

        if in_header
            && self
                .batch
                .last_reference_sequence_id
                .is_some_and(|last| last != reference_sequence_id)
            && !(self.many_references && self.is_small_container())
        {
            self.flush(header)?;
        }

        Ok(())
    }

    /// Whether the pending container is small enough to take records on another reference
    /// sequence, as HTSlib's `cram_put_bam_seq` decides for its `multi_seq` mode.
    fn is_small_container(&self) -> bool {
        self.batch.len() < self.context.records_per_slice / 4 + 10
    }

    /// Ends the container after a record of `read_length` bases was added if it is full.
    fn record_added(&mut self, header: &sam::Header, read_length: usize) -> io::Result<()> {
        self.bases += read_length;

        // A slice ends at its record or base limit, as in HTSlib, so long reads do not make
        // huge slices.
        if self.batch.len() >= self.records_per_container
            || self.bases >= self.context.bases_per_slice
        {
            self.flush(header)?;
        }

        Ok(())
    }

    fn flush(&mut self, header: &sam::Header) -> io::Result<()> {
        if self.failed {
            return Err(failed());
        }

        let record_count = u64::try_from(self.batch.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        // The next container takes records on other reference sequences if this one was small or
        // already had many reference sequences.
        self.many_references = self.batch.len() > 0
            && (self.is_small_container() || self.batch.has_many_reference_sequences);

        if let Some(workers) = &mut self.workers {
            if self.batch.len() == 0 {
                return Ok(());
            }

            let next = self.batch.next(self.records_per_container);
            let batch = mem::replace(&mut self.batch, next);

            workers.submit(&self.context, header, self.record_counter, batch);
            self.record_counter += record_count;
            self.bases = 0;

            let max = workers.capacity();
            return self.advance_workers(max);
        }

        let next = self.batch.next(self.records_per_container);
        let mut batch = mem::replace(&mut self.batch, next);
        let reference_sequences = mem::take(&mut batch.reference_sequences);
        let mut records = convert(header, batch, &reference_sequences)?;

        let result = write_container(
            self.inner
                .as_mut()
                .expect("the writer has its inner writer"),
            &self.context,
            self.compressor.as_mut(),
            header,
            self.record_counter,
            &mut records,
        );

        if let Err(e) = result {
            // The records stay pending.
            self.batch = Batch::from_converted(records);
            self.batch.reference_sequences = reference_sequences;
            return Err(e);
        }

        self.record_counter += record_count;
        self.bases = 0;

        Ok(())
    }

    /// Notes a failure of the workers.
    fn check(&mut self, result: io::Result<()>) -> io::Result<()> {
        if result.is_err() {
            self.failed = true;
            // The containers after the failed one are not written.
            self.workers = None;
        }

        result
    }
}

impl<W> Drop for Writer<W>
where
    W: Write,
{
    /// Finishes the stream (see [`Self::try_finish`]) if its header was written and it was not
    /// finished, ignoring errors. A panic of a worker is caught rather than resumed, since a
    /// destructor must not panic.
    fn drop(&mut self) {
        use std::panic::{self, AssertUnwindSafe};

        if self.inner.is_none() || thread::panicking() {
            return;
        }

        if let Some(header) = self.unfinished_header.take() {
            let _ = panic::catch_unwind(AssertUnwindSafe(|| self.try_finish(&header)));
        }
    }
}

/// The error of a call after a container failed.
fn failed() -> io::Error {
    io::Error::other("an earlier CRAM container failed")
}

impl<W> sam::alignment::io::Write for Writer<W>
where
    W: Write,
{
    fn write_alignment_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_alignment_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        if self.failed {
            return Err(failed());
        }

        let record = Record::try_from_alignment_record(header, record, |name| {
            get_reference_sequence(
                &self.context.reference_sequence_repository,
                &mut self.missing_reference_sequence_names,
                name,
            )
        })?;

        let read_length = record.read_length;
        self.start_record(header, record.reference_sequence_id)?;
        self.use_reference_sequence(header, record.reference_sequence_id)?;
        self.batch.push(record);
        self.record_added(header, read_length)
    }

    fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        self.try_finish(header)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use crate::core::Position;
    use crate::fasta;
    use crate::sam::alignment::record::{
        Flags,
        cigar::{Op, op::Kind},
        data::field::Tag,
    };

    use super::*;

    fn multi_slice_fixture() -> io::Result<(Vec<u8>, fasta::Repository, crate::cram::crai::Index)> {
        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(vec![b'A'; 64]),
        )]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:64\n".parse().unwrap();
        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.context.records_per_slice = 2;
        writer.write_header(&header)?;
        for i in 0..4 {
            let record = sam::alignment::RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::new(i + 1).unwrap())
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(b"AAAA".to_vec().into())
                .set_quality_scores(vec![30; 4].into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }
        writer.try_finish(&header)?;
        let src = writer.into_inner();
        let mut reader = crate::cram::io::Reader::new(Cursor::new(&src));
        reader.read_header()?;
        let offset = reader.get_ref().position();
        let mut container = crate::cram::io::reader::Container::default();
        reader.read_container(&mut container)?;
        assert_eq!(container.header().landmarks().len(), 2);
        let mut index = Vec::new();
        for (i, (&landmark, slice)) in container
            .header()
            .landmarks()
            .iter()
            .zip(container.slices())
            .enumerate()
        {
            let slice = slice?;
            let crate::cram::container::ReferenceSequenceContext::Some(context) =
                slice.header().reference_sequence_context()
            else {
                unreachable!()
            };
            index.push(crate::cram::crai::Record::new(
                Some(context.reference_sequence_id()),
                Some(context.alignment_start()),
                context.alignment_span(),
                offset,
                landmark as u64,
                (container
                    .header()
                    .landmarks()
                    .get(i + 1)
                    .copied()
                    .unwrap_or(container.src.len())
                    - landmark) as u64,
            ));
        }
        Ok((src, repository, index))
    }

    #[test]
    fn test_query_does_not_repeat_records_from_other_slices() -> io::Result<()> {
        let (src, repository, index) = multi_slice_fixture()?;
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;
        let region = "sq0:1-64".parse().unwrap();
        let records = reader
            .query(&header, &index, &region)?
            .records()
            .collect::<io::Result<Vec<_>>>()?;
        let names: Vec<_> = records
            .iter()
            .map(|r| r.name().unwrap().to_string())
            .collect();
        assert_eq!(names, ["r0", "r1", "r2", "r3"]);
        Ok(())
    }

    #[test]
    fn test_query_reads_bam_records_that_intersect_the_region() -> io::Result<()> {
        let (src, repository, index) = multi_slice_fixture()?;
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;

        // r0..r3 span [1, 4], [2, 5], [3, 6] and [4, 7], two to a slice.
        for (region, expected) in [
            ("sq0:1-2", &["r0", "r1"][..]),
            ("sq0:5-5", &["r1", "r2", "r3"][..]),
            ("sq0:7-64", &["r3"][..]),
            ("sq0:8-64", &[][..]),
        ] {
            let region = region.parse().unwrap();
            let mut query = reader.query(&header, &index, &region)?;
            let mut record = crate::bam::Record::default();
            let mut names = Vec::new();

            while query.read_record(&mut record)? != 0 {
                names.push(record.name().unwrap().to_string());
            }

            assert_eq!(names, expected, "{region:?}");
        }

        Ok(())
    }

    /// A CRAM with records on two reference sequences and unplaced unmapped records, in many
    /// containers of several slices each, with its index.
    fn many_containers_fixture()
    -> io::Result<(Vec<u8>, fasta::Repository, crate::cram::crai::Index)> {
        use crate::cram::container::ReferenceSequenceContext;

        let repository = fasta::Repository::new(vec![
            fasta::Record::new(
                fasta::record::Definition::new("sq0", None),
                fasta::record::Sequence::from(b"ACGT".repeat(32)),
            ),
            fasta::Record::new(
                fasta::record::Definition::new("sq1", None),
                fasta::record::Sequence::from(b"TTGCA".repeat(20)),
            ),
        ]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:128\n@SQ\tSN:sq1\tLN:100\n"
            .parse()
            .unwrap();

        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.context.records_per_slice = 3;
        writer.context.bases_per_slice = 60;
        writer.write_header(&header)?;

        let mut write = |name: String, id: Option<usize>, start: usize, len: usize| {
            let mut builder = sam::alignment::RecordBuf::builder()
                .set_name(name)
                .set_sequence(vec![b'A'; len].into())
                .set_quality_scores(vec![30; len].into());

            builder = match id {
                Some(id) => builder
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(id)
                    .set_alignment_start(Position::new(start).unwrap())
                    .set_cigar([Op::new(Kind::Match, len)].into_iter().collect()),
                None => builder.set_flags(Flags::UNMAPPED),
            };

            writer.write_alignment_record(&header, &builder.build())
        };

        for (id, length) in [(0, 128), (1, 100)] {
            for i in 0..40 {
                let len = 4 + (i * 7) % 13;
                let start = (1 + i * 3).min(length - len + 1);
                write(format!("sq{id}r{i}"), Some(id), start, len)?;
            }
        }

        for i in 0..5 {
            write(format!("u{i}"), None, 0, 6)?;
        }

        writer.try_finish(&header)?;
        let src = writer.into_inner();

        let mut reader = crate::cram::io::Reader::new(Cursor::new(&src));
        reader.read_header()?;
        let mut index = Vec::new();

        loop {
            let offset = reader.get_ref().position();
            let mut container = crate::cram::io::reader::Container::default();

            if reader.read_container(&mut container)? == 0 {
                break;
            }

            let landmarks = container.header().landmarks().to_vec();

            for (i, (&landmark, slice)) in landmarks.iter().zip(container.slices()).enumerate() {
                let slice = slice?;
                let end = landmarks.get(i + 1).copied().unwrap_or(container.src.len());
                let (id, start, span) = match slice.header().reference_sequence_context() {
                    ReferenceSequenceContext::Some(context) => (
                        Some(context.reference_sequence_id()),
                        Some(context.alignment_start()),
                        context.alignment_span(),
                    ),
                    ReferenceSequenceContext::None => (None, None, 0),
                    ReferenceSequenceContext::Many => {
                        unreachable!("a slice per reference sequence")
                    }
                };
                index.push(crate::cram::crai::Record::new(
                    id,
                    start,
                    span,
                    offset,
                    landmark as u64,
                    (end - landmark) as u64,
                ));
            }
        }

        Ok((src, repository, index))
    }

    #[test]
    fn test_suspended_query_reads_what_the_query_reads() -> io::Result<()> {
        use crate::cram::io::reader::Query;

        let (src, repository, index) = many_containers_fixture()?;
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;
        let mut record = crate::bam::Record::default();

        for region in ["sq0", "sq0:10-50", "sq1"] {
            let region = region.parse().unwrap();

            let mut expected = Vec::new();
            let mut query = reader.query(&header, &index, &region)?;
            while query.read_record(&mut record)? != 0 {
                expected.push(Vec::from(record.clone()));
            }
            assert!(!expected.is_empty());

            // A record at a time, setting the query aside after each.
            let mut actual = Vec::new();
            let mut suspended = reader.query(&header, &index, &region)?.suspend();
            loop {
                let mut query = Query::resume(&mut reader, &header, &index, suspended);
                if query.read_record(&mut record)? == 0 {
                    break;
                }
                actual.push(Vec::from(record.clone()));
                suspended = query.suspend();
            }

            assert_eq!(actual, expected);
        }

        Ok(())
    }

    #[test]
    fn test_fold_queries_matches_serial_queries() -> io::Result<()> {
        use std::num::NonZeroUsize;

        use crate::core::Region;
        use crate::cram::io::reader::Fields;

        let (src, repository, index) = many_containers_fixture()?;
        let containers: std::collections::BTreeSet<_> = index.iter().map(|r| r.offset()).collect();
        assert!(containers.len() > 4, "{} containers", containers.len());
        assert!(index.len() > containers.len(), "{} slices", index.len());

        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;

        let region_sets: &[&[&str]] = &[
            &["sq0"],
            &["sq0:10-50"],
            &["sq1:5-20", "sq0:100-128"],
            &["sq0:1-30", "sq0:20-60", "sq0:1-30"],
            &["sq0:200-300"],
            &["sq0", "sq1"],
            &[],
        ];

        for fields in [
            Fields::all(),
            Fields::REFERENCE_SEQUENCE_ID | Fields::ALIGNMENT_START | Fields::CIGAR,
        ] {
            reader.set_required_fields(fields);

            for regions in region_sets {
                let regions: Vec<Region> = regions.iter().map(|r| r.parse().unwrap()).collect();

                let mut expected = Vec::new();

                for (i, region) in regions.iter().enumerate() {
                    let mut query = reader.query(&header, &index, region)?;
                    query.for_each_record(|record| {
                        expected.push((i, record.name.as_ref().map(|n| n.to_vec())));
                        Ok::<_, io::Error>(())
                    })?;
                }

                if regions == ["sq0".parse().unwrap()] {
                    assert_eq!(expected.len(), 40);
                }

                for workers in [1, 2, 3, 8] {
                    let mut actual = Vec::new();

                    reader.fold_queries(
                        &header,
                        &index,
                        &regions,
                        NonZeroUsize::new(workers).unwrap(),
                        |names: &mut Vec<_>, i, record| {
                            names.push((i, record.name.as_ref().map(|n| n.to_vec())));
                            Ok(())
                        },
                        |names| {
                            actual.extend(names);
                            Ok::<_, io::Error>(())
                        },
                    )?;

                    assert_eq!(
                        actual, expected,
                        "{regions:?}, {workers} workers, {fields:?}"
                    );
                }
            }
        }

        Ok(())
    }

    #[test]
    fn test_fold_queries_stops_at_the_first_error() -> io::Result<()> {
        use std::num::NonZeroUsize;

        use crate::core::Region;

        let (src, repository, index) = many_containers_fixture()?;
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;
        let workers = NonZeroUsize::new(4).unwrap();
        let regions: Vec<Region> = vec!["sq0".parse().unwrap()];

        // `merge` fails on the third slice; nothing after it is merged.
        let mut merged = 0;
        let result = reader.fold_queries(
            &header,
            &index,
            &regions,
            workers,
            |_: &mut (), _, _| Ok(()),
            |()| {
                merged += 1;
                if merged == 3 {
                    Err(io::Error::other("stop"))
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(result.unwrap_err().to_string(), "stop");
        assert_eq!(merged, 3);

        // `fold` fails on a record; the slices before its slice are merged.
        let mut names = Vec::new();
        let result = reader.fold_queries(
            &header,
            &index,
            &regions,
            workers,
            |acc: &mut Vec<Vec<u8>>, _, record| {
                let name = record.name.as_ref().map(|n| n.to_vec()).unwrap_or_default();
                if name == b"sq0r20" {
                    return Err(io::Error::other("bad record"));
                }
                acc.push(name);
                Ok(())
            },
            |acc| {
                names.extend(acc);
                Ok::<_, io::Error>(())
            },
        );
        assert_eq!(result.unwrap_err().to_string(), "bad record");
        let prefix: Vec<_> = (0..names.len())
            .map(|i| format!("sq0r{i}").into_bytes())
            .collect();
        assert_eq!(names, prefix);
        assert!(
            !names.is_empty() && names.len() <= 20,
            "{} merged",
            names.len()
        );

        // An unknown reference sequence fails before anything is read.
        let regions: Vec<Region> = vec!["sq2".parse().unwrap()];
        let result = reader.fold_queries(
            &header,
            &index,
            &regions,
            workers,
            |_: &mut (), _, _| Ok(()),
            |()| Ok::<_, io::Error>(()),
        );
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::InvalidInput);

        Ok(())
    }

    #[test]
    fn test_fold_queries_bounds_slices_in_flight() -> io::Result<()> {
        use std::{
            num::NonZeroUsize,
            sync::atomic::{AtomicUsize, Ordering},
            thread,
            time::Duration,
        };

        use crate::core::Region;
        use crate::cram::io::reader::query::parallel::IN_FLIGHT_PER_WORKER;

        /// Slices whose decoding has started.
        static STARTED: AtomicUsize = AtomicUsize::new(0);

        /// The record names of a slice. Creating one counts a started slice.
        struct Names(Vec<Vec<u8>>);

        impl Default for Names {
            fn default() -> Self {
                STARTED.fetch_add(1, Ordering::SeqCst);
                Self(Vec::new())
            }
        }

        let (src, repository, index) = many_containers_fixture()?;
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;
        let regions: Vec<Region> = vec!["sq0".parse().unwrap(), "sq1".parse().unwrap()];
        // While one worker sleeps, the other decodes the rest.
        let workers = NonZeroUsize::new(2).unwrap();

        let mut names = Vec::new();
        let mut merged = 0;
        let mut most_retained = 0;

        // The first slice is slow, so the others finish first and wait for it.
        reader.fold_queries(
            &header,
            &index,
            &regions,
            workers,
            |acc: &mut Names, _, record| {
                let name = record.name.as_ref().map(|n| n.to_vec()).unwrap_or_default();

                if name == b"sq0r0" {
                    thread::sleep(Duration::from_millis(200));
                }

                acc.0.push(name);
                Ok(())
            },
            |acc| {
                // Slices started and not merged before this one, this one included.
                most_retained = most_retained.max(STARTED.load(Ordering::SeqCst) - merged);
                merged += 1;
                names.extend(acc.0);
                Ok::<_, io::Error>(())
            },
        )?;

        assert_eq!(names.len(), 80);
        assert!(merged > 3 * IN_FLIGHT_PER_WORKER, "{merged} slices");

        let bound = workers.get() * IN_FLIGHT_PER_WORKER;
        assert!(
            most_retained <= bound,
            "{most_retained} slices retained, bound {bound}"
        );

        Ok(())
    }

    #[test]
    fn test_fold_queries_resumes_a_fold_panic() -> io::Result<()> {
        use std::{
            num::NonZeroUsize,
            panic::{self, AssertUnwindSafe},
        };

        use crate::core::Region;

        let (src, repository, index) = many_containers_fixture()?;
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(src));
        let header = reader.read_header()?;
        let regions: Vec<Region> = vec!["sq0".parse().unwrap()];
        let workers = NonZeroUsize::new(2).unwrap();

        let mut names = Vec::new();

        let result = panic::catch_unwind(AssertUnwindSafe(|| {
            reader.fold_queries(
                &header,
                &index,
                &regions,
                workers,
                |acc: &mut Vec<Vec<u8>>, _, record| {
                    let name = record.name.as_ref().map(|n| n.to_vec()).unwrap_or_default();
                    assert_ne!(name, b"sq0r20", "bad record");
                    acc.push(name);
                    Ok(())
                },
                |acc| {
                    names.extend(acc);
                    Ok::<_, io::Error>(())
                },
            )
        }));

        let payload = result.expect_err("fold panicked");
        let message = payload
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or_default();
        assert!(message.contains("bad record"), "{message}");

        // The slices before the panicking one were merged.
        let prefix: Vec<_> = (0..names.len())
            .map(|i| format!("sq0r{i}").into_bytes())
            .collect();
        assert_eq!(names, prefix);
        assert!(!names.is_empty() && names.len() <= 20, "{}", names.len());

        Ok(())
    }

    #[test]
    fn test_write_records_on_a_reference_sequence_missing_from_the_repository() -> io::Result<()> {
        use crate::sam::alignment::RecordBuf;

        // The repository has sq0 but not sq1.
        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGTACGTAC".to_vec()),
        )]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:10\n@SQ\tSN:sq1\tLN:10\n".parse().unwrap();

        let record =
            |name: &str, id: Option<usize>, start: usize, cigar: &[(Kind, usize)], seq: &[u8]| {
                let mut builder = RecordBuf::builder()
                    .set_name(name)
                    .set_flags(if id.is_some() {
                        Flags::empty()
                    } else {
                        Flags::UNMAPPED
                    })
                    .set_sequence(seq.to_vec().into())
                    .set_quality_scores(vec![30; seq.len()].into());

                if let Some(id) = id {
                    builder = builder
                        .set_reference_sequence_id(id)
                        .set_alignment_start(Position::new(start).unwrap())
                        .set_cigar(
                            cigar
                                .iter()
                                .map(|&(kind, len)| Op::new(kind, len))
                                .collect(),
                        );
                }

                builder.build()
            };

        let records = [
            // A mismatch against sq0, which is stored verbatim with the sq1 records.
            record("r0", Some(0), 1, &[(Kind::Match, 4)], b"ACGA"),
            record(
                "r1",
                Some(1),
                2,
                &[
                    (Kind::HardClip, 2),
                    (Kind::SoftClip, 1),
                    (Kind::Match, 2),
                    (Kind::Insertion, 1),
                    (Kind::Deletion, 1),
                    (Kind::Match, 2),
                    (Kind::SoftClip, 1),
                ],
                b"TTGCAAG",
            ),
            record("r2", Some(1), 5, &[(Kind::Match, 5)], b""),
            record("r3", None, 0, &[], b"AC"),
        ];

        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.many_references = true;
        writer.write_header(&header)?;

        for record in &records {
            writer.write_alignment_record(&header, record)?;
        }

        writer.try_finish(&header)?;

        // The reader needs no sq1 (or sq0): the bases are stored verbatim.
        let mut reader = crate::cram::io::Reader::new(Cursor::new(writer.into_inner()));
        let actual_header = reader.read_header()?;
        let actual = reader
            .records(&actual_header)
            .collect::<io::Result<Vec<_>>>()?;

        assert_eq!(actual.len(), records.len());

        for (actual, expected) in actual.iter().zip(&records) {
            assert_eq!(actual.name(), expected.name());
            assert_eq!(
                actual.reference_sequence_id(),
                expected.reference_sequence_id()
            );
            assert_eq!(actual.alignment_start(), expected.alignment_start());
            assert_eq!(actual.cigar(), expected.cigar());
            assert_eq!(actual.sequence(), expected.sequence());
            assert!(actual.data().get(&Tag::MISMATCHED_POSITIONS).is_none());
        }

        Ok(())
    }

    #[test]
    fn test_write_records_with_an_embedded_consensus_reference() -> io::Result<()> {
        use crate::sam::alignment::RecordBuf;

        use crate::cram::io::reader::Container;

        // The repository lacks sq0, so the slice embeds a consensus of its records.
        let header: sam::Header = "@SQ\tSN:sq0\tLN:16\n".parse().unwrap();

        let records: Vec<_> = [(1, "ACGTAC"), (3, "GTTCGG"), (3, "GTACGG"), (5, "ACGGNA")]
            .into_iter()
            .enumerate()
            .map(|(i, (start, bases))| {
                RecordBuf::builder()
                    .set_name(format!("r{i}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::new(start).unwrap())
                    .set_cigar([Op::new(Kind::Match, bases.len())].into_iter().collect())
                    .set_sequence(bases.as_bytes().to_vec().into())
                    .set_quality_scores(vec![30; bases.len()].into())
                    .build()
            })
            .collect();

        let mut writer = Builder::default().build_from_writer(Vec::new());
        writer.write_header(&header)?;

        for record in &records {
            writer.write_alignment_record(&header, record)?;
        }

        writer.try_finish(&header)?;
        let src = writer.into_inner();

        let mut reader = crate::cram::io::Reader::new(Cursor::new(&src));
        let actual_header = reader.read_header()?;
        let mut container = Container::default();
        reader.read_container(&mut container)?;
        let slice = container.slices().next().expect("a slice")?;
        assert!(
            slice
                .header()
                .embedded_reference_bases_block_content_id()
                .is_some()
        );

        let mut reader = crate::cram::io::Reader::new(Cursor::new(&src));
        reader.read_header()?;
        let actual = reader
            .records(&actual_header)
            .collect::<io::Result<Vec<_>>>()?;

        assert_eq!(actual.len(), records.len());

        for (actual, expected) in actual.iter().zip(&records) {
            assert_eq!(actual.alignment_start(), expected.alignment_start());
            assert_eq!(actual.sequence(), expected.sequence());
            // No MD or NM from the consensus, and no `cF` tag.
            assert!(actual.data().is_empty());
        }

        Ok(())
    }

    #[test]
    fn test_read_with_required_tags() -> io::Result<()> {
        use crate::cram::io::reader::{Fields, Tags};
        use crate::sam::alignment::{
            RecordBuf,
            record::data::field::Tag,
            record_buf::{
                Data,
                data::field::{Value, value::Array},
            },
        };

        const HP: Tag = Tag::new(b'H', b'P');
        const XS: Tag = Tag::new(b'X', b'S');

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGTACGTACGTACGT".to_vec()),
        )]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:16\n@RG\tID:rg0\n".parse().unwrap();

        let data: Data = [
            (Tag::BASE_MODIFICATIONS, Value::from("C+m?,0;")),
            (
                Tag::BASE_MODIFICATION_PROBABILITIES,
                Value::Array(Array::UInt8(vec![200])),
            ),
            (HP, Value::from(1u8)),
            (XS, Value::from(7u8)),
            (Tag::READ_GROUP, Value::from("rg0")),
        ]
        .into_iter()
        .collect();

        // A mismatch at the third base, so that MD and NM (which the writer leaves out) are
        // regenerated.
        let record = RecordBuf::builder()
            .set_name("r0")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(1).unwrap())
            .set_cigar([Op::new(Kind::Match, 8)].into_iter().collect())
            .set_sequence(b"ACCTACGT".to_vec().into())
            .set_quality_scores(vec![30; 8].into())
            .set_data(data)
            .build();

        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;
        writer.write_alignment_record(&header, &record)?;
        writer.try_finish(&header)?;
        let src = writer.into_inner();

        let read = |fields: Fields, tags: Tags| -> io::Result<Vec<Tag>> {
            let mut reader = crate::cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_required_fields(fields)
                .set_required_tags(tags)
                .build_from_reader(&src[..]);
            reader.read_header()?;
            let record = reader.records(&header).next().unwrap()?;
            Ok(record.data().iter().map(|(tag, _)| tag).collect())
        };

        let all = read(Fields::all(), Tags::All)?;
        assert_eq!(
            all,
            [
                Tag::BASE_MODIFICATIONS,
                Tag::BASE_MODIFICATION_PROBABILITIES,
                HP,
                XS,
                Tag::MISMATCHED_POSITIONS,
                Tag::EDIT_DISTANCE,
                Tag::READ_GROUP,
            ]
        );

        // The rest keep their order.
        let except = Tags::Except(vec![
            Tag::BASE_MODIFICATIONS,
            Tag::BASE_MODIFICATION_PROBABILITIES,
        ]);
        assert_eq!(read(Fields::all(), except)?, &all[2..]);

        // NM is regenerated without MD.
        let only = Tags::Only(vec![HP, Tag::EDIT_DISTANCE]);
        assert_eq!(
            read(Fields::all(), only)?,
            [HP, Tag::EDIT_DISTANCE, Tag::READ_GROUP]
        );

        // The read group alone is no tag.
        assert_eq!(
            read(Fields::FLAGS | Fields::READ_GROUP, Tags::All)?,
            [Tag::READ_GROUP]
        );
        assert_eq!(
            read(Fields::FLAGS | Fields::DATA, Tags::Only(vec![XS]))?,
            [XS]
        );
        assert!(read(Fields::FLAGS, Tags::All)?.is_empty());

        Ok(())
    }

    #[test]
    fn test_write_records_without_the_md_and_nm_that_readers_regenerate() -> io::Result<()> {
        use crate::sam::alignment::{RecordBuf, record_buf::data::field::Value};

        const XS: Tag = Tag::new(b'X', b'S');

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGTACGTACGTACGT".to_vec()),
        )]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:16\n".parse().unwrap();

        // A mismatch at the third base. The second record's MD is wrong, so it is kept.
        let record = |md: &str| {
            RecordBuf::builder()
                .set_name("r0")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::MIN)
                .set_cigar([Op::new(Kind::Match, 8)].into_iter().collect())
                .set_sequence(b"ACCTACGT".to_vec().into())
                .set_quality_scores(vec![30; 8].into())
                .set_data(
                    [
                        (Tag::MISMATCHED_POSITIONS, Value::from(md)),
                        (Tag::EDIT_DISTANCE, Value::from(1u8)),
                        (XS, Value::from(7u8)),
                    ]
                    .into_iter()
                    .collect(),
                )
                .build()
        };

        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;
        writer.write_alignment_record(&header, &record("2G5"))?;
        writer.write_alignment_record(&header, &record("2T5"))?;
        writer.try_finish(&header)?;
        let src = writer.into_inner();

        assert!(!src.windows(3).any(|w| w == b"2G5"));
        assert!(src.windows(3).any(|w| w == b"2T5"));

        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(&src[..]);
        reader.read_header()?;

        let data: Vec<Vec<(Tag, Value)>> = reader
            .records(&header)
            .map(|result| {
                result.map(|record| record.data().iter().map(|(t, v)| (t, v.clone())).collect())
            })
            .collect::<io::Result<_>>()?;

        // Readers append the regenerated tags.
        assert_eq!(
            data,
            [
                vec![
                    (XS, Value::from(7u8)),
                    (Tag::MISMATCHED_POSITIONS, Value::from("2G5")),
                    (Tag::EDIT_DISTANCE, Value::from(1u8)),
                ],
                vec![
                    (Tag::MISMATCHED_POSITIONS, Value::from("2T5")),
                    (XS, Value::from(7u8)),
                    (Tag::EDIT_DISTANCE, Value::from(1u8)),
                ],
            ]
        );

        Ok(())
    }

    #[test]
    fn test_release_reference_sequences_of_sorted_records() -> io::Result<()> {
        use std::num::NonZero;

        use crate::sam::alignment::RecordBuf;

        let names = [b"sq0", b"sq1", b"sq2"];

        let write = |sort_order: &str, ids: &[usize], workers: usize| -> io::Result<_> {
            let repository = fasta::Repository::new(
                names
                    .iter()
                    .map(|name| {
                        fasta::Record::new(
                            fasta::record::Definition::new(&name[..], None),
                            fasta::record::Sequence::from(b"ACGTACGT".to_vec()),
                        )
                    })
                    .collect::<Vec<_>>(),
            );

            let header: sam::Header = format!(
                "@HD\tVN:1.6\tSO:{sort_order}\n@SQ\tSN:sq0\tLN:8\n@SQ\tSN:sq1\tLN:8\n@SQ\tSN:sq2\tLN:8\n"
            )
            .parse()
            .unwrap();

            // A container per record.
            let options = compression::Options {
                records_per_slice: 1,
                ..Default::default()
            };

            let mut builder = Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_compression_options(options);

            if let Some(n) = NonZero::new(workers) {
                builder = builder.set_worker_count(n);
            }

            let mut writer = builder.build_from_writer(Vec::new());
            writer.write_header(&header)?;

            for &id in ids {
                let record = RecordBuf::builder()
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(id)
                    .set_alignment_start(Position::MIN)
                    .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                    .set_sequence(b"ACGT".to_vec().into())
                    .build();

                writer.write_alignment_record(&header, &record)?;
            }

            writer.try_finish(&header)?;

            let cached: Vec<_> = names
                .iter()
                .map(|name| repository.is_cached(&name[..]))
                .collect();

            Ok((writer.into_inner(), cached))
        };

        let mut outputs = Vec::new();

        for workers in [0, 2] {
            let (src, cached) = write("coordinate", &[0, 0, 1, 2], workers)?;
            assert_eq!(cached, [false, false, true]);
            outputs.push(src);

            let (_, cached) = write("unsorted", &[0, 0, 1, 2], workers)?;
            assert_eq!(cached, [true, true, true]);

            // Records that are not sorted after all keep the rest loaded. (The first, released, may
            // still be in use on a worker, and then it is not loaded again.)
            let (_, cached) = write("coordinate", &[0, 1, 0, 2], workers)?;
            assert_eq!(cached[1..], [true, true]);
        }

        assert_eq!(outputs[0], outputs[1]);

        Ok(())
    }

    /// As in HTSlib, sorted records start a new container on each reference sequence,
    /// so that their slices have one reference sequence and a reference checksum; only short runs
    /// of records on a reference sequence (e.g., unsorted data, small contigs) share a container.
    #[test]
    fn test_write_starts_a_container_on_a_reference_sequence_change() -> io::Result<()> {
        use crate::cram::container::ReferenceSequenceContext;
        use crate::sam::alignment::RecordBuf;

        let names = ["sq0", "sq1", "sq2", "sq3"];
        let repository = fasta::Repository::new(
            names
                .iter()
                .map(|name| {
                    fasta::Record::new(
                        fasta::record::Definition::new(*name, None),
                        fasta::record::Sequence::from(b"ACGT".repeat(64)),
                    )
                })
                .collect::<Vec<_>>(),
        );
        let header: sam::Header = names
            .iter()
            .map(|name| format!("@SQ\tSN:{name}\tLN:256\n"))
            .collect::<String>()
            .parse()
            .unwrap();

        // Containers (or slices) of 100 records: sq0 and sq1 have long runs, then sq2 and sq3 short
        // ones, then sq0 and sq1 alternate (unsorted).
        let mut ids: Vec<usize> = [0; 150].into_iter().chain([1; 150]).collect();
        ids.extend([2; 5].into_iter().chain([3; 5]));
        ids.extend((0..20).map(|i| i % 2));

        let mut writer = Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());
        writer.context.records_per_slice = 100;
        writer.write_header(&header)?;

        for (i, &id) in ids.iter().enumerate() {
            let record = RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(Flags::empty())
                .set_reference_sequence_id(id)
                .set_alignment_start(Position::new(1 + i % 200).unwrap())
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(b"ACGT".to_vec().into())
                .set_quality_scores(vec![30; 4].into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }

        writer.try_finish(&header)?;
        let src = writer.into_inner();

        let mut reader = crate::cram::io::Reader::new(Cursor::new(&src));
        reader.read_header()?;
        let mut container = crate::cram::io::reader::Container::default();
        let mut slices = Vec::new();

        while reader.read_container(&mut container)? != 0 {
            for slice in container.slices() {
                let slice = slice?;
                let header = slice.header();

                slices.push(match header.reference_sequence_context() {
                    ReferenceSequenceContext::Some(context) => {
                        assert!(header.reference_md5().is_some());
                        (Some(context.reference_sequence_id()), header.record_count())
                    }
                    ReferenceSequenceContext::Many => (None, header.record_count()),
                    ReferenceSequenceContext::None => unreachable!(),
                });
            }
        }

        assert_eq!(
            slices,
            [
                (Some(0), 100),
                (Some(0), 50),
                (Some(1), 100),
                (Some(1), 50),
                // A new container after the long run, then many reference sequences.
                (Some(2), 5),
                (None, 25),
            ]
        );

        // The records read back in order.
        let mut reader = crate::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(Cursor::new(&src));
        let header = reader.read_header()?;
        let actual: Vec<_> = reader
            .records(&header)
            .map(|record| record.map(|record| record.reference_sequence_id()))
            .collect::<io::Result<_>>()?;
        let expected: Vec<_> = ids.into_iter().map(Some).collect();
        assert_eq!(actual, expected);

        Ok(())
    }

    /// An unmapped record named `r{i}`.
    fn unmapped_record(i: usize) -> crate::sam::alignment::RecordBuf {
        crate::sam::alignment::RecordBuf::builder()
            .set_name(format!("r{i}"))
            .set_flags(Flags::UNMAPPED)
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(vec![30; 4].into())
            .build()
    }

    /// The names of the records in `src`, which must end with an EOF container.
    fn read_names(src: &[u8]) -> io::Result<Vec<String>> {
        use self::container::EOF;

        assert!(src.ends_with(&EOF), "the stream ends with an EOF container");
        assert!(
            !src[..src.len() - EOF.len()].ends_with(&EOF),
            "the stream has one EOF container"
        );

        let mut reader = crate::cram::io::Reader::new(src);
        let header = reader.read_header()?;

        reader
            .records(&header)
            .map(|record| {
                let record = record?;
                Ok(record
                    .name()
                    .map(|name| name.to_string())
                    .unwrap_or_default())
            })
            .collect()
    }

    #[test]
    fn test_drop_finishes_the_stream() -> io::Result<()> {
        use std::num::NonZero;

        let header = sam::Header::default();
        let expected: Vec<_> = (0..25).map(|i| format!("r{i}")).collect();

        for workers in [0, 2] {
            // Several containers, so that some are in flight on the workers at the drop.
            let options = compression::Options {
                records_per_slice: 4,
                ..Default::default()
            };

            fn build_with(
                options: compression::Options,
                workers: usize,
                dst: &mut Vec<u8>,
            ) -> Writer<&mut Vec<u8>> {
                let mut builder = Builder::default().set_compression_options(options);

                if let Some(n) = NonZero::new(workers) {
                    builder = builder.set_worker_count(n);
                }

                builder.build_from_writer(dst)
            }

            let mut src = Vec::new();

            {
                let mut writer = build_with(options, workers, &mut src);
                writer.write_header(&header)?;

                for i in 0..25 {
                    writer.write_alignment_record(&header, &unmapped_record(i))?;
                }
            }

            assert_eq!(read_names(&src)?, expected, "{workers} workers");

            // A finished writer writes nothing more when dropped.
            let mut finished = Vec::new();

            {
                let mut writer = build_with(options, workers, &mut finished);
                writer.write_header(&header)?;

                for i in 0..25 {
                    writer.write_alignment_record(&header, &unmapped_record(i))?;
                }

                writer.try_finish(&header)?;
            }

            assert_eq!(finished, src);

            // Nothing is written without a header, and the inner writer is taken unfinished.
            let mut empty = Vec::new();
            drop(build_with(options, workers, &mut empty));
            assert!(empty.is_empty());

            let mut unfinished = Vec::new();
            let mut writer = build_with(options, workers, &mut unfinished);
            writer.write_header(&header)?;
            writer.into_inner();
            assert!(!unfinished.ends_with(&self::container::EOF));
        }

        Ok(())
    }

    #[test]
    fn test_drop_ignores_errors() {
        struct Failing;

        impl Write for Failing {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::other("disk full"))
            }

            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let header = sam::Header::default();
        let mut writer = Builder::default()
            .set_worker_count(std::num::NonZero::new(2).unwrap())
            .build_from_writer(Failing);
        writer.unfinished_header = Some(header.clone());
        writer.batch.push(
            Record::try_from_alignment_record(&header, &unmapped_record(0), |_| {
                Ok(None::<&fasta::record::Sequence>)
            })
            .unwrap(),
        );

        // Neither the pending container nor the EOF container can be written; dropping does not
        // panic.
        drop(writer);
    }

    #[test]
    fn test_write_fails_on_an_m5_that_does_not_match_the_reference() -> io::Result<()> {
        use std::num::NonZero;

        use crate::sam::alignment::RecordBuf;

        let repository = fasta::Repository::new(vec![
            fasta::Record::new(
                fasta::record::Definition::new("sq0", None),
                fasta::record::Sequence::from(b"ACGTACGT".to_vec()),
            ),
            fasta::Record::new(
                fasta::record::Definition::new("sq1", None),
                fasta::record::Sequence::from(b"TTGCAAGG".to_vec()),
            ),
        ]);

        let sq0_md5 = crate::sam::header::record::value::map::reference_sequence::Md5Checksum::from(
            crate::cram::calculate_normalized_sequence_digest(b"ACGTACGT"),
        );

        // sq0's M5 matches; sq1's is of another sequence.
        let header: sam::Header = format!(
            "@SQ\tSN:sq0\tLN:8\tM5:{sq0_md5}\n@SQ\tSN:sq1\tLN:8\tM5:00004142844c68362f611ca22888656a\n\
             @SQ\tSN:sq2\tLN:8\n"
        )
        .parse()
        .unwrap();

        let record = |reference_sequence_id| {
            RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(Position::MIN)
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(b"ACGT".to_vec().into())
                .set_quality_scores(vec![30; 4].into())
                .build()
        };

        for (workers, ignore_md5) in [(0, false), (2, false), (0, true), (2, true)] {
            let options = compression::Options {
                ignore_md5,
                ..Default::default()
            };

            let mut builder = Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_compression_options(options);

            if let Some(n) = NonZero::new(workers) {
                builder = builder.set_worker_count(n);
            }

            let mut writer = builder.build_from_writer(Vec::new());
            writer.write_header(&header)?;

            // A matching M5, and a sequence the repository lacks.
            writer.write_alignment_record(&header, &record(0))?;

            let mut bam_record = bam::Record::default();
            let mut encoder = bam::io::Writer::from(Vec::new());
            encoder.write_alignment_record(&header, &record(1))?;
            bam_record.0 = encoder.get_ref()[4..].to_vec();

            let result = writer.write_bam_record(&header, &bam_record);

            if ignore_md5 {
                result?;
                writer.try_finish(&header)?;
            } else {
                let e = result.unwrap_err();
                assert_eq!(e.kind(), io::ErrorKind::InvalidData);
                let message = e.to_string();
                assert!(
                    message.starts_with("SQ header M5 tag discrepancy for reference 'sq1'"),
                    "{message}"
                );
                assert!(
                    message.contains("has 00004142844c68362f611ca22888656a"),
                    "{message}"
                );
                writer.into_inner();
            }
        }

        Ok(())
    }

    #[test]
    fn test_write_records_with_stored_md_and_nm() -> io::Result<()> {
        use crate::sam::alignment::{RecordBuf, record_buf::data::field::Value};

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGTACGTACGTACGT".to_vec()),
        )]);
        let header: sam::Header = "@SQ\tSN:sq0\tLN:16\n".parse().unwrap();

        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_cigar([Op::new(Kind::Match, 8)].into_iter().collect())
            .set_sequence(b"ACCTACGT".to_vec().into())
            .set_quality_scores(vec![30; 8].into())
            .set_data(
                [
                    (Tag::MISMATCHED_POSITIONS, Value::from("2G5")),
                    (Tag::EDIT_DISTANCE, Value::from(1u8)),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        // The tags stored, read as stored.
        let write = |store_md, store_nm| -> io::Result<Vec<Tag>> {
            let options = compression::Options {
                store_md,
                store_nm,
                ..Default::default()
            };

            let mut writer = Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_compression_options(options)
                .build_from_writer(Vec::new());
            writer.write_header(&header)?;
            writer.write_alignment_record(&header, &record)?;
            writer.try_finish(&header)?;
            let src = writer.into_inner();

            let mut reader = crate::cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_md_nm(crate::cram::io::reader::MdNm::AsStored)
                .build_from_reader(&src[..]);
            reader.read_header()?;

            let record = reader.records(&header).next().expect("a record")?;
            Ok(record.data().iter().map(|(tag, _)| tag).collect())
        };

        assert!(write(false, false)?.is_empty());
        assert_eq!(write(true, false)?, [Tag::MISMATCHED_POSITIONS]);
        assert_eq!(write(false, true)?, [Tag::EDIT_DISTANCE]);
        assert_eq!(
            write(true, true)?,
            [Tag::MISMATCHED_POSITIONS, Tag::EDIT_DISTANCE]
        );

        Ok(())
    }
}
