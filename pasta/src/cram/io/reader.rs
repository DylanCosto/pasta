//! CRAM reader and record iterator.

mod builder;
pub(crate) mod collections;
pub(crate) mod container;
mod fields;
pub mod header;
mod md_nm;
mod multithreaded_records;
pub(crate) mod num;
pub(crate) mod query;
mod records;

use std::{
    io::{self, Read, Seek, SeekFrom},
    num::NonZero,
};

use crate::core::Region;
use crate::fasta;
use crate::sam;

pub use self::{
    builder::Builder,
    container::{Block, Container, read_block},
    fields::{Fields, Tags},
    md_nm::MdNm,
    multithreaded_records::MultithreadedRecords,
    query::{Query, SuspendedQuery, UnmappedQuery},
    records::Records,
};
use self::{container::read_container, header::read_header};
use crate::cram::{FileDefinition, crai};

/// A CRAM reader.
///
/// The CRAM format is comprised of four main parts: 1) a file definition, 2) a file header, 3) a
/// list of containers, and 4) an end-of-file (EOF) container.
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io};
/// use pasta::cram;
/// use pasta::fasta;
///
/// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
/// let header = reader.read_header()?;
///
/// for result in reader.records(&header) {
///     let record = result?;
///     // ...
/// }
///
/// # Ok::<_, io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
    reference_sequence_repository: fasta::Repository,
    md_nm: MdNm,
    required_fields: Fields,
    required_tags: Tags,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let reader = cram::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CRAM reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    pub(crate) fn reference_sequence_repository(&self) -> &fasta::Repository {
        &self.reference_sequence_repository
    }

    /// Replaces the reference sequence repository, for example with one found from the header's
    /// `@SQ` `UR` fields after reading it.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram;
    /// use pasta::fasta;
    /// let mut reader = cram::io::Reader::new(&[][..]);
    /// reader.set_reference_sequence_repository(fasta::Repository::default());
    /// ```
    pub fn set_reference_sequence_repository(&mut self, repository: fasta::Repository) {
        self.reference_sequence_repository = repository;
    }

    /// Returns whether decoded records get the MD and NM tags the file left out.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{self as cram, io::reader::MdNm};
    /// let reader = cram::io::Reader::new(&[][..]);
    /// assert_eq!(reader.md_nm(), MdNm::Regenerate);
    /// ```
    pub fn md_nm(&self) -> MdNm {
        self.md_nm
    }

    /// Sets whether decoded records get the MD and NM tags the file left out.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{self as cram, io::reader::MdNm};
    /// let mut reader = cram::io::Reader::new(&[][..]);
    /// reader.set_md_nm(MdNm::AsStored);
    /// assert_eq!(reader.md_nm(), MdNm::AsStored);
    /// ```
    pub fn set_md_nm(&mut self, md_nm: MdNm) {
        self.md_nm = md_nm;
    }

    /// Returns the record fields this reader decodes.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{self as cram, io::reader::Fields};
    /// let reader = cram::io::Reader::new(&[][..]);
    /// assert_eq!(reader.required_fields(), Fields::all());
    /// ```
    pub fn required_fields(&self) -> Fields {
        self.required_fields
    }

    /// Sets the record fields this reader decodes (`CRAM_OPT_REQUIRED_FIELDS`); see [`Fields`].
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{self as cram, io::reader::Fields};
    /// let mut reader = cram::io::Reader::new(&[][..]);
    /// reader.set_required_fields(Fields::FLAGS);
    /// assert_eq!(reader.required_fields(), Fields::FLAGS);
    /// ```
    pub fn set_required_fields(&mut self, fields: Fields) {
        self.required_fields = fields;
    }

    /// Returns the tags this reader decodes with [`Fields::DATA`].
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{self as cram, io::reader::Tags};
    /// let reader = cram::io::Reader::new(&[][..]);
    /// assert_eq!(reader.required_tags(), &Tags::All);
    /// ```
    pub fn required_tags(&self) -> &Tags {
        &self.required_tags
    }

    /// Sets the tags this reader decodes with [`Fields::DATA`]; see [`Tags`].
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{self as cram, io::reader::Tags};
    /// use pasta::sam::alignment::record::data::field::Tag;
    ///
    /// let mut reader = cram::io::Reader::new(&[][..]);
    /// reader.set_required_tags(Tags::Only(vec![Tag::EDIT_DISTANCE]));
    /// assert_eq!(reader.required_tags(), &Tags::Only(vec![Tag::EDIT_DISTANCE]));
    /// ```
    pub fn set_required_tags(&mut self, tags: Tags) {
        self.required_tags = tags;
    }

    /// Returns a CRAM header reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io::Read};
    /// use pasta::cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let mut header_reader = reader.header_reader();
    /// header_reader.read_magic_number()?;
    /// header_reader.read_format_version()?;
    /// header_reader.read_file_id()?;
    ///
    /// let mut container_reader = header_reader.container_reader()?;
    ///
    /// let _raw_header = {
    ///     let mut raw_sam_header_reader = container_reader.raw_sam_header_reader()?;
    ///     let mut raw_header = String::new();
    ///     raw_sam_header_reader.read_to_string(&mut raw_header)?;
    ///     raw_sam_header_reader.discard_to_end()?;
    ///     raw_header
    /// };
    ///
    /// container_reader.discard_to_end()?;
    /// Ok::<_, std::io::Error>(())
    /// ```
    pub fn header_reader(&mut self) -> header::Reader<&mut R> {
        header::Reader::new(&mut self.inner)
    }

    /// Reads the CRAM file definition.
    ///
    /// The CRAM magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let file_definition = reader.read_file_definition()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_definition(&mut self) -> io::Result<FileDefinition> {
        header::read_file_definition(&mut self.inner)
    }

    /// Reads the SAM header.
    ///
    /// The position of the stream is expected to be at the CRAM header container, i.e., directly
    /// after the file definition.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// reader.read_file_definition()?;
    ///
    /// let header = reader.read_file_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_header(&mut self) -> io::Result<sam::Header> {
        header::read_file_header(&mut self.inner)
    }

    /// Reads the SAM header.
    ///
    /// This verifies the CRAM magic number, discards the file definition, and reads and parses the
    /// file header as a SAM header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram;
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        read_header(&mut self.inner)
    }

    /// Reads a container.
    ///
    /// If successful, the container size is returned. If a container size of 0 is returned, the
    /// stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram::{self as cram, io::reader::Container};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_container(&mut self, container: &mut Container) -> io::Result<usize> {
        read_container(&mut self.inner, container)
    }

    /// Returns a iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a container.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram;
    /// use pasta::fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// for result in reader.records(&header) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'r, 'h: 'r>(&'r mut self, header: &'h sam::Header) -> Records<'r, 'h, R> {
        Records::new(self, header)
    }

    /// Returns an iterator over records starting from the current stream position that decodes
    /// slices on `worker_count` threads.
    ///
    /// It gives the records [`Self::records`] gives, in the same order; when a slice cannot be
    /// decoded, it gives the records of the slices before it in its container first, where
    /// [`Self::records`] gives none of that container's records. Containers are read on
    /// the calling thread, and each slice is decoded by a job on a pool of workers (fewer of them
    /// when the system cannot start that many threads, or none, decoding on the calling thread,
    /// when it cannot start any). A few slices per worker are read ahead of the records given
    /// out, so the stream is past them when the iterator is dropped.
    ///
    /// The stream is expected to be at the start of a container.
    ///
    /// A panic of a decoding job is re-raised on the calling thread when the iterator reaches
    /// that slice. After an error, the iterator ends.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io, num::NonZero};
    /// use pasta::cram;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    /// let header = reader.read_header()?;
    /// let worker_count = NonZero::new(4).unwrap();
    ///
    /// for result in reader.multithreaded_records(&header, worker_count) {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn multithreaded_records(
        &mut self,
        header: &sam::Header,
        worker_count: NonZero<usize>,
    ) -> MultithreadedRecords<'_, R> {
        MultithreadedRecords::new(self, header, worker_count)
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Seeks the underlying reader to the given position.
    ///
    /// Positions typically come from the associated CRAM index file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io::{self, SeekFrom};
    /// use pasta::cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// reader.seek(SeekFrom::Start(0))?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos)
    }

    /// Returns the current position of the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram;
    /// let mut reader = cram::io::Reader::new(io::empty());
    /// let position = reader.position()?;
    /// assert_eq!(position, 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn position(&mut self) -> io::Result<u64> {
        self.inner.stream_position()
    }

    /// Returns a reader over records that intersects the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram::{self as cram, crai};
    /// use pasta::fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i crai::Index,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, 'i, R>> {
        let reference_sequence_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "invalid reference sequence name",
                )
            })?;

        Ok(Query::new(
            self,
            header,
            index,
            reference_sequence_id,
            region.interval(),
        ))
    }

    /// Returns an iterator of unplaced unmapped records (`*`) after querying for the unmapped
    /// region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram::{self as cram, crai};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let query = reader.query_unmapped(&header, &index)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn query_unmapped<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &crai::Index,
    ) -> io::Result<UnmappedQuery<'r, 'h, R>> {
        let offset = index
            .iter()
            .find(|record| record.reference_sequence_id().is_none())
            .map(|record| record.offset());

        // Without unplaced records, there is nothing to read: the end of the stream has no
        // container (not even the EOF container) to read.
        let is_done = match offset {
            Some(offset) => {
                self.get_mut().seek(SeekFrom::Start(offset))?;
                false
            }
            None => true,
        };

        Ok(UnmappedQuery::new(self, header, is_done))
    }
}

impl<R> sam::alignment::io::Read<R> for Reader<R>
where
    R: Read,
{
    fn read_alignment_header(&mut self) -> io::Result<sam::Header> {
        self.read_header()
    }

    fn alignment_records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'a> {
        Box::new(
            self.records(header).map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            }),
        )
    }
}
