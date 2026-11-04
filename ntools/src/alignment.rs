//! Format-specific alignment readers behind one streaming interface.
//!
//! Each format keeps its fastest native path: BAM records stay lazy and can be forwarded as raw
//! bytes, SAM records are parsed lazily, and CRAM records are decoded into record buffers.

use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{self, BufRead, BufReader, Read},
    mem,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use pasta::bam;
use pasta::core::{Position, Region, region::Interval};
use pasta::cram::{
    self as cram,
    io::reader::{Fields, MdNm, Tags},
};
use pasta::fasta;
use pasta::sam::{
    self as sam,
    alignment::{RecordBuf, record::Flags},
};

use crate::{
    bgzf_stream::BgzfReader,
    parallel_cram,
    regions::{RegionSet, parse_query_region, record_interval, region_string},
    source::{Compression, Format, Input, STDIN, Stream},
};

/// The region name samtools uses for unplaced, unmapped records.
const UNMAPPED_REGION: &str = "*";

/// Options shared by every command that reads alignments.
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    /// Additional decompression threads.
    pub threads: usize,
    /// A FASTA reference for CRAM decoding.
    pub reference: Option<PathBuf>,
}

/// A borrowed alignment record from any supported format.
pub enum Rec<'a> {
    Bam(&'a bam::Record),
    Sam(&'a sam::Record),
    Buf(&'a RecordBuf),
    /// Any other record, e.g., a CRAM record decoded on a worker.
    Dyn(&'a dyn sam::alignment::Record),
}

impl Rec<'_> {
    /// Returns the record as a format-independent trait object.
    pub fn as_dyn(&self) -> &dyn sam::alignment::Record {
        match self {
            Self::Bam(r) => *r,
            Self::Sam(r) => *r,
            Self::Buf(r) => *r,
            Self::Dyn(r) => *r,
        }
    }

    pub fn flags(&self) -> io::Result<Flags> {
        match self {
            Self::Bam(r) => Ok(r.flags()),
            Self::Sam(r) => r.flags(),
            Self::Buf(r) => Ok(r.flags()),
            Self::Dyn(r) => r.flags(),
        }
    }

    /// Returns the mapping quality, with `255` (missing) as `None`.
    pub fn mapping_quality(&self) -> io::Result<Option<u8>> {
        match self {
            Self::Bam(r) => Ok(r.mapping_quality().map(u8::from)),
            Self::Sam(r) => r.mapping_quality().transpose().map(|q| q.map(u8::from)),
            Self::Buf(r) => Ok(r.mapping_quality().map(u8::from)),
            Self::Dyn(r) => r.mapping_quality().transpose().map(|q| q.map(u8::from)),
        }
    }

    pub fn reference_sequence_id(&self, header: &sam::Header) -> io::Result<Option<usize>> {
        self.as_dyn().reference_sequence_id(header).transpose()
    }

    pub fn mate_reference_sequence_id(&self, header: &sam::Header) -> io::Result<Option<usize>> {
        self.as_dyn().mate_reference_sequence_id(header).transpose()
    }

    /// Returns the record's BAM encoding (without `block_size`): a BAM record's own bytes, or
    /// any other record encoded into `encoder`.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be encoded.
    pub fn encode_raw<'r>(
        &'r self,
        header: &sam::Header,
        encoder: &'r mut RawEncoder,
    ) -> io::Result<&'r [u8]> {
        use sam::alignment::io::Write as _;

        match self {
            Self::Bam(r) => Ok(r.as_ref()),
            _ => {
                encoder.writer.get_mut().clear();
                encoder
                    .writer
                    .write_alignment_record(header, self.as_dyn())?;
                Ok(&encoder.writer.get_ref()[mem::size_of::<u32>()..])
            }
        }
    }
}

/// A reusable buffer for [`Rec::encode_raw`].
pub struct RawEncoder {
    writer: bam::io::Writer<Vec<u8>>,
}

impl RawEncoder {
    pub fn new() -> Self {
        Self {
            writer: bam::io::Writer::from(Vec::new()),
        }
    }
}

impl Default for RawEncoder {
    fn default() -> Self {
        Self::new()
    }
}

type BamReader = bam::io::Reader<BgzfReader<BufReader<Stream>>>;

/// An open alignment reader for SAM, BAM, or CRAM.
pub enum AlignmentReader {
    /// A SAM reader and the input's name (`-` for stdin), which its errors name.
    Sam(sam::io::Reader<Box<dyn BufRead + Send>>, String),
    Bam(BamReader),
    /// A CRAM reader, the reference repository it decodes with, where the reference comes
    /// from, and how many worker threads decode its containers (0 for none).
    Cram(
        cram::io::Reader<BufReader<Stream>>,
        fasta::Repository,
        ReferenceSource,
        usize,
    ),
}

/// Where a CRAM reader's reference sequences come from.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceSource {
    /// A FASTA file given with `--reference`.
    Given,
    /// The FASTA files named by the header's `@SQ` `UR` fields, as HTSlib finds them, for the
    /// input of the given name.
    Header(String),
}

/// Warns, as HTSlib's `bam_hdr_read` does, when a BAM file does not end with the BGZF EOF
/// marker, which usually means it was truncated. Streams that cannot seek are not checked.
fn warn_if_eof_marker_is_absent(stream: &Stream) {
    use std::os::unix::fs::FileExt;

    use crate::bam_blocks::BGZF_EOF;

    let Stream::File(file, _) = stream else {
        return;
    };

    let Ok(len) = file.metadata().map(|metadata| metadata.len()) else {
        return;
    };

    let mut buf = [0; BGZF_EOF.len()];
    let has_eof = len
        .checked_sub(BGZF_EOF.len() as u64)
        .is_some_and(|offset| file.read_exact_at(&mut buf, offset).is_ok() && buf == BGZF_EOF);

    if !has_eof {
        eprintln!("[W::bam_hdr_read] EOF marker is absent. The input is probably truncated");
    }
}

impl AlignmentReader {
    /// Opens `src` (`-` for stdin), detecting its format.
    pub fn open(src: &str, options: &ReadOptions) -> anyhow::Result<(Self, Option<PathBuf>)> {
        let input = Input::open(src)?;

        if let Some(other_format) = input.other_format {
            bail!("{src}: {other_format}");
        }

        let path = input.path.clone();
        let reader = Self::from_input(input, options)?;
        Ok((reader, path))
    }

    pub fn from_input(input: Input, options: &ReadOptions) -> anyhow::Result<Self> {
        let Input {
            path,
            reader,
            format,
            compression,
            ..
        } = input;
        let name = path
            .as_deref()
            .map_or_else(|| STDIN.to_owned(), |path| path.display().to_string());

        Ok(match (format, compression) {
            (Format::Bam, _) => {
                warn_if_eof_marker_is_absent(reader.get_ref());

                Self::Bam(bam::io::Reader::from(BgzfReader::new(
                    reader,
                    options.threads,
                )))
            }
            (Format::Sam, Compression::Bgzf) => {
                let inner: Box<dyn BufRead + Send> =
                    Box::new(BgzfReader::new(reader, options.threads));
                Self::Sam(sam::io::Reader::new(inner), name)
            }
            (Format::Sam, Compression::None) => {
                let inner: Box<dyn BufRead + Send> = Box::new(reader);
                Self::Sam(sam::io::Reader::new(inner), name)
            }
            (Format::Sam, Compression::Gzip) => {
                let decoder = flate2::bufread::MultiGzDecoder::new(reader);
                let inner: Box<dyn BufRead + Send> = Box::new(BufReader::new(decoder));
                Self::Sam(sam::io::Reader::new(inner), name)
            }
            (Format::Cram, _) => {
                let name = input_name(path.as_deref());

                // Without `--reference`, the header's is set once it is read.
                let (repository, source) = match &options.reference {
                    Some(src) => {
                        let adapter = indexed_reference(src)?;
                        let repository = fasta::Repository::new(ExplainedReference {
                            adapter: Some(Box::new(adapter)),
                            input: name,
                            tried: Tried::Given(src.clone()),
                        });
                        load_dictionary(&repository, src)?;
                        (repository, ReferenceSource::Given)
                    }
                    None => (fasta::Repository::default(), ReferenceSource::Header(name)),
                };

                Self::Cram(
                    cram::io::reader::Builder::default()
                        .set_reference_sequence_repository(repository.clone())
                        .build_from_reader(reader),
                    repository,
                    source,
                    options.threads,
                )
            }
        })
    }

    /// Stops CRAM decoding from adding the MD and NM tags the file left out, as samtools
    /// commands that set `CRAM_OPT_DECODE_MD` to 0 do. Other formats are unaffected.
    pub fn keep_stored_md_nm(&mut self) {
        self.set_md_nm(MdNm::AsStored);
    }

    /// Sets whether CRAM decoding adds the MD and NM tags the file left out (`decode_md`).
    /// Other formats are unaffected.
    pub fn set_md_nm(&mut self, md_nm: MdNm) {
        if let Self::Cram(reader, ..) = self {
            reader.set_md_nm(md_nm);
        }
    }

    /// Sets the record fields and tags CRAM decoding needs (`CRAM_OPT_REQUIRED_FIELDS`), so that
    /// it skips the data only other fields and tags use. The other fields of CRAM records are left unset; other
    /// formats are unaffected.
    pub fn set_required_fields(&mut self, fields: Fields, tags: Tags) {
        if let Self::Cram(reader, ..) = self {
            reader.set_required_fields(fields);
            reader.set_required_tags(tags);
        }
    }

    /// Sets how many worker threads decode CRAM containers for [`Self::for_each_while`] (0 for
    /// none); [`Self::open`] uses the thread count it is given. Other formats are unaffected.
    pub fn set_cram_decoding_threads(&mut self, threads: usize) {
        if let Self::Cram(.., n) = self {
            *n = threads;
        }
    }

    pub fn format(&self) -> Format {
        match self {
            Self::Sam(..) => Format::Sam,
            Self::Bam(_) => Format::Bam,
            Self::Cram(..) => Format::Cram,
        }
    }

    pub fn read_header(&mut self) -> anyhow::Result<sam::Header> {
        let header = match self {
            Self::Sam(r, name) => r
                .read_header()
                .with_context(|| format!("reading the header of {name}")),
            Self::Bam(r) => r.read_header().context("reading header"),
            Self::Cram(r, ..) => r.read_header().context("reading header"),
        }?;

        if let Self::Cram(r, repository, ReferenceSource::Header(input), _) = self {
            let adapter = header_references(&header)?;
            *repository = fasta::Repository::new(ExplainedReference {
                adapter: adapter.map(|a| Box::new(a) as Box<dyn fasta::repository::Adapter>),
                input: input.clone(),
                tried: Tried::Header(header_uris(&header)),
            });
            r.set_reference_sequence_repository(repository.clone());
        }

        Ok(header)
    }

    /// Calls `f` with every record in file order.
    pub fn for_each<F>(&mut self, header: &sam::Header, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        self.for_each_while(header, |rec| f(rec).map(|()| true))
    }

    /// Calls `f` with records in file order until it returns `false`.
    pub fn for_each_while<F>(&mut self, header: &sam::Header, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<bool>,
    {
        match self {
            Self::Bam(reader) => {
                let mut record = bam::Record::default();
                while reader.read_record(&mut record)? != 0 {
                    if !f(Rec::Bam(&record))? {
                        return Ok(());
                    }
                }
            }
            Self::Sam(reader, name) => {
                // Records are normalized as HTSlib's `sam_parse1` normalizes them.
                let mut record = sam::Record::default();
                for n in 1_u64.. {
                    let len = reader
                        .read_normalized_record(header, &mut record)
                        .with_context(|| format!("{name}: reading SAM record {n}"))?;

                    if len == 0 || !f(Rec::Sam(&record))? {
                        return Ok(());
                    }
                }
            }
            Self::Cram(reader, repository, _, threads)
                if *threads > 0 && reader.required_fields().is_all() =>
            {
                return for_each_in_parallel(reader, repository, header, *threads, f);
            }
            Self::Cram(reader, repository, ..) => {
                // Decoding containers directly avoids converting each record to a `RecordBuf`, and
                // encoding records as BAM lets callers use their raw BAM paths. Records with only
                // some fields decoded are passed as they are, since their encoding would be
                // incomplete anyway.
                let mut container = cram::io::reader::Container::default();
                let mut bam_record = bam::Record::default();
                let fields = reader.required_fields();
                let tags = reader.required_tags().clone();

                while reader.read_container(&mut container)? != 0 {
                    let compression_header = container.compression_header()?;

                    for result in container.slices() {
                        let slice = result?;
                        let (core_data_src, external_data_srcs) =
                            slice.decode_blocks_for(&compression_header, fields, &tags)?;

                        let records = slice.records_with_fields(
                            repository.clone(),
                            header,
                            &compression_header,
                            &core_data_src,
                            &external_data_srcs,
                            reader.md_nm(),
                            fields,
                            &tags,
                        )?;

                        for record in &records {
                            let rec = if fields.is_all() {
                                encode_cram_record(header, record, &mut bam_record)?;
                                Rec::Bam(&bam_record)
                            } else {
                                Rec::Dyn(record)
                            };

                            if !f(rec)? {
                                return Ok(());
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Calls `f` with the records of each indexed query, in query order.
    ///
    /// `index` names an index file explicitly (`-X`); otherwise the index next to `path` is used.
    pub fn for_each_in_queries<F>(
        &mut self,
        path: &Path,
        index: Option<&Path>,
        header: &sam::Header,
        queries: &Queries<'_>,
        f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let index = AlignmentIndex::read(self.format(), path, index)?;
        self.for_each_in_indexed_queries(&index, header, queries, f)
    }

    /// Like [`Self::for_each_in_queries`], with an index already read, for callers that query
    /// the same file many times.
    pub fn for_each_in_indexed_queries<F>(
        &mut self,
        index: &AlignmentIndex,
        header: &sam::Header,
        queries: &Queries<'_>,
        f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let queries = queries.to_list(header)?;
        self.run_queries(index, header, &queries, f)
    }

    /// Calls `f` with the records overlapping the 0-based, half-open interval `[start, end)` of
    /// one reference sequence, as `sam_itr_queryi` selects them. An `end` of `None` means the
    /// end of the reference sequence.
    pub fn for_each_in_interval<F>(
        &mut self,
        index: &AlignmentIndex,
        header: &sam::Header,
        reference_sequence_id: usize,
        start: usize,
        end: Option<usize>,
        f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let Some(name) = header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .map(|(name, _)| name)
        else {
            return Ok(());
        };

        let name = std::str::from_utf8(name).context("reference sequence name is not UTF-8")?;
        let start = Position::try_from(start + 1)?;

        let region = match end {
            Some(end) if end < start.get() => return Ok(()),
            Some(end) => Region::new(name, start..=Position::try_from(end)?),
            None if start == Position::MIN => Region::new(name, ..),
            None => Region::new(name, start..),
        };

        self.run_queries(index, header, &[Query::Region(region)], f)
    }

    fn run_queries<F>(
        &mut self,
        index: &AlignmentIndex,
        header: &sam::Header,
        queries: &[Query<'_>],
        mut f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let mut visit = |record: Rec<'_>, query: &Query| -> anyhow::Result<()> {
            if query.visits(record.as_dyn(), header)? {
                f(record)
            } else {
                Ok(())
            }
        };

        match (self, index) {
            (Self::Bam(reader), AlignmentIndex::Bam(index)) => {
                for query in queries {
                    match query {
                        Query::Unmapped => {
                            for result in reader.query_unmapped(index)? {
                                let record = result?;
                                visit(Rec::Bam(&record), query)?;
                            }
                        }
                        Query::Region(region) | Query::Interval { region, .. } => {
                            let mut query_reader = reader.query(header, index, region)?;
                            let mut record = bam::Record::default();

                            while query_reader.read_record(&mut record)? != 0 {
                                visit(Rec::Bam(&record), query)?;
                            }
                        }
                    }
                }
            }
            (Self::Cram(reader, ..), AlignmentIndex::Cram(index)) => {
                // Like BAM queries, CRAM queries yield raw BAM records, encoded from the decoded
                // CRAM fields directly.
                let mut record = bam::Record::default();

                for query in queries {
                    match query {
                        Query::Unmapped => {
                            let mut query_reader = reader.query_unmapped(header, index)?;

                            while query_reader.read_record(&mut record)? != 0 {
                                visit(Rec::Bam(&record), query)?;
                            }
                        }
                        // Records with only some fields decoded are visited as they are, since
                        // their BAM encoding would be incomplete (and is not needed).
                        Query::Region(region) | Query::Interval { region, .. }
                            if !reader.required_fields().is_all() =>
                        {
                            let mut query_reader = reader.query(header, index, region)?;
                            query_reader
                                .for_each_record(|record| visit(Rec::Dyn(record), query))?;
                        }
                        Query::Region(region) | Query::Interval { region, .. } => {
                            let mut query_reader = reader.query(header, index, region)?;

                            while query_reader.read_record(&mut record)? != 0 {
                                visit(Rec::Bam(&record), query)?;
                            }
                        }
                    }
                }
            }
            (Self::Sam(..), _) => {
                bail!("region queries on SAM input are not supported; use BAM or CRAM")
            }
            _ => bail!("the index does not match the input format"),
        }

        Ok(())
    }
}

/// Decodes the containers of a CRAM reader on `threads` workers, which encode the records as
/// BAM, and calls `f` with the records in file order until it returns `false`.
fn for_each_in_parallel<F>(
    reader: &mut cram::io::Reader<BufReader<Stream>>,
    repository: &fasta::Repository,
    header: &sam::Header,
    threads: usize,
    mut f: F,
) -> anyhow::Result<()>
where
    F: FnMut(Rec<'_>) -> anyhow::Result<bool>,
{
    /// A container's records, BAM-encoded one after another, and where each ends.
    #[derive(Default)]
    struct Encoded {
        buf: Vec<u8>,
        ends: Vec<usize>,
        record: Vec<u8>,
    }

    let mut stopped = false;
    let mut bam_record = bam::Record::default();

    let result = parallel_cram::decode(
        reader,
        header,
        repository,
        parallel_cram::Source::All,
        threads,
        |encoded: &mut Encoded, record| {
            record.encode_bam_record(header, &mut encoded.record)?;
            encoded.buf.extend_from_slice(&encoded.record);
            encoded.ends.push(encoded.buf.len());
            Ok(())
        },
        |encoded| {
            let mut start = 0;

            for end in encoded.ends {
                let mut buf = Vec::from(mem::take(&mut bam_record));
                buf.clear();
                buf.extend_from_slice(&encoded.buf[start..end]);
                bam_record = bam::Record::try_from(buf)?;
                start = end;

                if !f(Rec::Bam(&bam_record))? {
                    // An error stops decoding; this one is not reported.
                    stopped = true;
                    anyhow::bail!("stopped");
                }
            }

            Ok(())
        },
    );

    match result {
        Err(_) if stopped => Ok(()),
        result => result,
    }
}

/// Encodes a CRAM record as a BAM record into `dst`, reusing its buffer.
pub fn encode_cram_record(
    header: &sam::Header,
    record: &cram::Record<'_>,
    dst: &mut bam::Record,
) -> io::Result<()> {
    let mut buf = Vec::from(mem::take(dst));
    record.encode_bam_record(header, &mut buf)?;
    *dst = bam::Record::try_from(buf)?;
    Ok(())
}

/// A BAM (BAI or CSI) or CRAM (CRAI) index.
pub enum AlignmentIndex {
    Bam(bam::Index),
    Cram(cram::crai::Index),
}

impl AlignmentIndex {
    /// Reads the index named by `explicit`, or the one next to `path`, for an input of `format`.
    pub fn read(format: Format, path: &Path, explicit: Option<&Path>) -> anyhow::Result<Self> {
        match format {
            Format::Bam => read_bam_index(path, explicit).map(Self::Bam),
            Format::Cram => match explicit {
                Some(index) => cram::crai::fs::read(index)
                    .with_context(|| format!("reading index for {}", path.display())),
                None => read_associated_index(path, &["crai"], |index| cram::crai::fs::read(index)),
            }
            .map(Self::Cram),
            Format::Sam => bail!("region queries on SAM input are not supported; use BAM or CRAM"),
        }
    }
}

/// What to query through an index.
pub enum Queries<'a> {
    /// samtools regions, queried in the order given. Without `-M`, a record overlapping several
    /// regions is visited once per region.
    Regions(&'a [String]),
    /// Merged intervals (`-M`, `-L`), visiting each overlapping record once, in file order.
    Set(&'a RegionSet),
}

enum Query<'a> {
    Unmapped,
    Region(Region),
    Interval {
        region: Region,
        set: &'a RegionSet,
        reference_sequence_id: usize,
        index: usize,
    },
}

impl Query<'_> {
    /// Whether this query visits a record it reads. In a region set, a record is visited only
    /// while querying the first merged interval it overlaps, so each record is visited once and
    /// in file order.
    fn visits(
        &self,
        record: &dyn sam::alignment::Record,
        header: &sam::Header,
    ) -> io::Result<bool> {
        let Query::Interval {
            set,
            reference_sequence_id,
            index,
            ..
        } = self
        else {
            return Ok(true);
        };

        Ok(match record_interval(record, header)? {
            Some((id, start, end)) => {
                id == *reference_sequence_id && set.first_overlap(id, start, end) == Some(*index)
            }
            None => false,
        })
    }
}

/// Indexed CRAM queries, with their index, in the form
/// [`cram::io::Reader::fold_queries`] decodes on worker threads.
pub struct CramQueries<'a> {
    index: cram::crai::Index,
    queries: Vec<Query<'a>>,
    regions: Vec<Region>,
}

impl<'a> CramQueries<'a> {
    /// Lists `queries` and reads the index of the CRAM at `path` (`index` names one explicitly).
    /// Returns `None` if a query is for unplaced unmapped records (`*`), which only a serial
    /// query reads.
    pub fn new(
        path: &Path,
        index: Option<&Path>,
        header: &sam::Header,
        queries: &Queries<'a>,
    ) -> anyhow::Result<Option<Self>> {
        let queries = queries.to_list(header)?;

        let regions = queries
            .iter()
            .map(|query| match query {
                Query::Unmapped => None,
                Query::Region(region) | Query::Interval { region, .. } => Some(region.clone()),
            })
            .collect::<Option<Vec<_>>>();

        let Some(regions) = regions else {
            return Ok(None);
        };

        let AlignmentIndex::Cram(index) = AlignmentIndex::read(Format::Cram, path, index)? else {
            unreachable!("a CRAM index is read for CRAM");
        };

        Ok(Some(Self {
            index,
            queries,
            regions,
        }))
    }

    pub fn index(&self) -> &cram::crai::Index {
        &self.index
    }

    /// The region of each query.
    pub fn regions(&self) -> &[Region] {
        &self.regions
    }

    /// Whether query `i` (an index into [`Self::regions`]) visits a record it reads.
    pub fn visits(
        &self,
        i: usize,
        record: &cram::Record<'_>,
        header: &sam::Header,
    ) -> io::Result<bool> {
        self.queries[i].visits(record, header)
    }
}

impl<'a> Queries<'a> {
    fn to_list(&self, header: &sam::Header) -> anyhow::Result<Vec<Query<'a>>> {
        match *self {
            Self::Regions(raw_regions) => raw_regions
                .iter()
                .map(|raw| {
                    if raw == UNMAPPED_REGION {
                        Ok(Query::Unmapped)
                    } else {
                        parse_query_region(header, raw).map(Query::Region)
                    }
                })
                .collect(),
            Self::Set(set) => {
                let mut queries = Vec::new();

                for (id, name) in header.reference_sequences().keys().enumerate() {
                    for (i, &(start, end)) in set.intervals(id).iter().enumerate() {
                        queries.push(Query::Interval {
                            region: region_string(name, start, end)?,
                            set,
                            reference_sequence_id: id,
                            index: i,
                        });
                    }
                }

                Ok(queries)
            }
        }
    }
}

/// Reads an explicit BAI or CSI index (by extension), or the one associated with `path`.
fn read_bam_index(path: &Path, index: Option<&Path>) -> anyhow::Result<bam::Index> {
    match index {
        Some(index) => {
            read_bai_or_csi(index).with_context(|| format!("reading index for {}", path.display()))
        }
        None => read_associated_index(path, &["bai", "csi"], read_bai_or_csi),
    }
}

/// Reads the first index next to `path` that exists, looking, as HTSlib does, for `path` with
/// each of `extensions` added and then with its own extension replaced (`in.bam.bai`, then
/// `in.bai`). When there is none, the error names the files looked for and how to make one.
fn read_associated_index<T>(
    path: &Path,
    extensions: &[&str],
    read: impl Fn(&Path) -> io::Result<T>,
) -> anyhow::Result<T> {
    let candidates: Vec<PathBuf> = extensions
        .iter()
        .flat_map(|extension| {
            let replaced = path
                .extension()
                .is_some()
                .then(|| path.with_extension(extension));
            [Some(path.with_added_extension(extension)), replaced]
        })
        .flatten()
        .collect();

    for candidate in &candidates {
        match read(candidate) {
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            result => {
                return result.with_context(|| {
                    format!(
                        "reading index {} for {}",
                        candidate.display(),
                        path.display()
                    )
                });
            }
        }
    }

    let looked_for = candidates
        .iter()
        .map(|candidate| candidate.display().to_string())
        .collect::<Vec<_>>()
        .join(", ");

    bail!(
        "no index for {path} (looked for {looked_for}); run \"ntools index {path}\"",
        path = path.display()
    )
}

/// Reads a BAI or CSI index, told apart by content as HTSlib does: a BAI starts with its magic
/// number, and a CSI is BGZF-compressed.
fn read_bai_or_csi(path: &Path) -> io::Result<bam::Index> {
    const BAI_MAGIC: &[u8; 4] = b"BAI\x01";

    let mut magic = [0; 4];
    File::open(path)?.read_exact(&mut magic)?;

    if &magic == BAI_MAGIC {
        bam::bai::fs::read(path).map(bam::Index::Bai)
    } else {
        pasta::csi::fs::read(path).map(bam::Index::Csi)
    }
}

/// The local file an `@SQ` `UR` field names, as HTSlib's `cram_populate_ref` reads it: a `file:`
/// prefix is dropped, and other URLs are not local files.
fn local_uri_path(uri: &str) -> Option<PathBuf> {
    match uri.strip_prefix("file:") {
        Some(path) => Some(PathBuf::from(path)),
        None if uri.contains("://") => None,
        None => Some(PathBuf::from(uri)),
    }
}

/// Whether `path` is a file with a `.fai` index next to it, as HTSlib needs to read a `UR` file.
fn is_indexed_fasta(path: &Path) -> bool {
    let mut index = path.as_os_str().to_owned();
    index.push(".fai");
    path.is_file() && Path::new(&index).is_file()
}

/// Reads the reference sequences of a CRAM header from the local FASTA files its `@SQ` lines
/// name in their `UR` fields (those that exist, with an index), as HTSlib does when no reference
/// is given. Returns `None` if no line names one.
///
/// Each sequence is read from the file its own line names, or else from the first named file that
/// has it, since HTSlib finds a sequence in any file it has loaded.
pub fn header_repository(header: &sam::Header) -> anyhow::Result<Option<fasta::Repository>> {
    Ok(header_references(header)?.map(fasta::Repository::new))
}

/// See [`header_repository`].
fn header_references(header: &sam::Header) -> anyhow::Result<Option<HeaderReferences>> {
    use sam::header::record::value::map::reference_sequence::tag::URI;

    let mut paths: Vec<PathBuf> = Vec::new();
    let mut own_files = Vec::new();

    for (name, reference_sequence) in header.reference_sequences() {
        let path = reference_sequence
            .other_fields()
            .get(&URI)
            .and_then(|uri| local_uri_path(&uri.to_string()))
            .filter(|path| is_indexed_fasta(path));

        let file = path.map(|path| match paths.iter().position(|p| *p == path) {
            Some(i) => i,
            None => {
                paths.push(path);
                paths.len() - 1
            }
        });

        own_files.push((name.to_vec(), file));
    }

    if paths.is_empty() {
        return Ok(None);
    }

    let readers = paths
        .iter()
        .map(|path| {
            fasta::io::indexed_reader::Builder::default()
                .build_from_path(path)
                .with_context(|| format!("opening reference {}", path.display()))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    let files = {
        let names: Vec<HashSet<&[u8]>> = readers
            .iter()
            .map(|reader| {
                reader
                    .index()
                    .as_ref()
                    .iter()
                    .map(|record| record.name().as_ref())
                    .collect()
            })
            .collect();

        own_files
            .into_iter()
            .filter_map(|(name, own_file)| {
                let has = |&i: &usize| names[i].contains(name.as_slice());
                let file = own_file
                    .filter(has)
                    .or_else(|| (0..names.len()).find(has))?;
                Some((name, file))
            })
            .collect()
    };

    Ok(Some(HeaderReferences {
        readers: readers
            .into_iter()
            .map(fasta::repository::adapters::IndexedReader::new)
            .collect(),
        files,
    }))
}

/// The `@SQ` `UR` fields of a header, by reference sequence name.
fn header_uris(header: &sam::Header) -> HashMap<Vec<u8>, String> {
    use sam::header::record::value::map::reference_sequence::tag::URI;

    header
        .reference_sequences()
        .iter()
        .filter_map(|(name, reference_sequence)| {
            let uri = reference_sequence.other_fields().get(&URI)?;
            Some((name.to_vec(), uri.to_string()))
        })
        .collect()
}

/// A CRAM reader's reference. A sequence it lacks fails with what was tried and how to fix it,
/// as HTSlib names the file it tried and points to its documentation, where the CRAM reader would
/// only say that the sequence is missing.
struct ExplainedReference {
    adapter: Option<Box<dyn fasta::repository::Adapter>>,
    /// The CRAM file's name, for the error.
    input: String,
    tried: Tried,
}

/// Where [`ExplainedReference`] looked for reference sequences.
enum Tried {
    /// The FASTA file given with `--reference`.
    Given(PathBuf),
    /// The files the header's `@SQ` `UR` fields name, by reference sequence name.
    Header(HashMap<Vec<u8>, String>),
}

impl ExplainedReference {
    fn missing(&self, name: &[u8]) -> io::Error {
        use bstr::ByteSlice;

        let input = &self.input;
        let name = name.as_bstr();

        let message = match &self.tried {
            Tried::Given(src) => format!(
                "{input}: no reference for {name} (tried -T {}, no {name} in it); give the \
                 file's reference with -T/--reference",
                src.display()
            ),
            Tried::Header(uris) => {
                let tried = match uris.get(name.as_bytes()) {
                    Some(uri) => format!("tried UR:{uri}, {}", uri_problem(uri, name)),
                    None => String::from("no @SQ UR"),
                };

                format!("{input}: no reference for {name} ({tried}); give one with -T/--reference")
            }
        };

        io::Error::new(io::ErrorKind::NotFound, message)
    }
}

/// Why the `UR` file `uri` gave no reference sequence `name`.
fn uri_problem(uri: &str, name: &bstr::BStr) -> String {
    match local_uri_path(uri) {
        None => String::from("not a local file"),
        Some(path) if !path.exists() => String::from("not found"),
        Some(path) if !is_indexed_fasta(&path) => String::from("not an indexed (.fai) FASTA file"),
        Some(_) => format!("no {name} in it"),
    }
}

impl fasta::repository::Adapter for ExplainedReference {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let result = self.adapter.as_mut().and_then(|adapter| adapter.get(name));
        Some(result.unwrap_or_else(|| Err(self.missing(name))))
    }

    fn get_range(
        &mut self,
        name: &[u8],
        interval: Interval,
    ) -> Option<io::Result<fasta::record::Sequence>> {
        self.adapter.as_mut()?.get_range(name, interval)
    }

    fn len(&mut self, name: &[u8]) -> Option<io::Result<usize>> {
        self.adapter.as_mut()?.len(name)
    }
}

/// The name of an input for messages: its path, or `-` for stdin.
fn input_name(path: Option<&Path>) -> String {
    path.map_or_else(|| String::from("-"), |path| path.display().to_string())
}

/// A repository adapter that reads each reference sequence from one of several indexed FASTA
/// files.
struct HeaderReferences {
    readers: Vec<fasta::repository::adapters::IndexedReader<fasta::io::BufReader<File>>>,
    /// The index in `readers` of the file each sequence is read from, by name.
    files: HashMap<Vec<u8>, usize>,
}

impl fasta::repository::Adapter for HeaderReferences {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        let &i = self.files.get(name)?;
        self.readers[i].get(name)
    }

    fn get_range(
        &mut self,
        name: &[u8],
        interval: Interval,
    ) -> Option<io::Result<fasta::record::Sequence>> {
        let &i = self.files.get(name)?;
        self.readers[i].get_range(name, interval)
    }
}

/// Builds a reference sequence repository from an indexed FASTA, or an empty one.
pub fn reference_repository(src: Option<&Path>) -> anyhow::Result<fasta::Repository> {
    let Some(src) = src else {
        return Ok(fasta::Repository::default());
    };

    let repository = fasta::Repository::new(indexed_reference(src)?);
    load_dictionary(&repository, src)?;
    Ok(repository)
}

/// Reads the sequence dictionary beside a reference FASTA (`ref.dict`, from `ntools dict`, Picard
/// or GATK), whose checksums check the reference for CRAM slices that carry none.
fn load_dictionary(repository: &fasta::Repository, src: &Path) -> anyhow::Result<()> {
    repository
        .load_dictionary(src)
        .with_context(|| format!("reading the sequence dictionary of {}", src.display()))?;
    Ok(())
}

/// Opens an indexed FASTA file as a repository adapter, indexing it first if it has no index.
fn indexed_reference(
    src: &Path,
) -> anyhow::Result<fasta::repository::adapters::IndexedReader<fasta::io::BufReader<File>>> {
    crate::reference::ensure_index(src)?;

    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(src)
        .with_context(|| format!("opening reference {}", src.display()))?;

    Ok(fasta::repository::adapters::IndexedReader::new(reader))
}

#[cfg(test)]
mod tests {
    use std::fs;

    use pasta::sam::alignment::io::Write as _;

    use super::*;

    /// Writes a FASTA file with its `.fai` index.
    fn write_indexed_fasta(path: &Path, records: &[(&str, &str)]) -> io::Result<()> {
        let mut fasta = String::new();
        let mut index = String::new();

        for (name, sequence) in records {
            let header = format!(">{name}\n");
            let offset = fasta.len() + header.len();
            let length = sequence.len();
            index.push_str(&format!(
                "{name}\t{length}\t{offset}\t{length}\t{}\n",
                length + 1
            ));
            fasta.push_str(&header);
            fasta.push_str(sequence);
            fasta.push('\n');
        }

        fs::write(path, fasta)?;
        let mut index_path = path.as_os_str().to_owned();
        index_path.push(".fai");
        fs::write(index_path, index)
    }

    #[test]
    fn test_header_repository_reads_each_sequence_from_its_ur_file() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");
        let unindexed = dir.path().join("c.fa");
        write_indexed_fasta(
            &a,
            &[("chr1", "ACGTACGT"), ("chr3", "GGGG"), ("chr2", "TTTT")],
        )?;
        write_indexed_fasta(&b, &[("chr2", "CCCCAAAA")])?;
        fs::write(&unindexed, ">chr4\nAAAA\n")?;

        let header: sam::Header = format!(
            "@SQ\tSN:chr1\tLN:8\tUR:file://{a}\n\
             @SQ\tSN:chr2\tLN:8\tUR:{b}\n\
             @SQ\tSN:chr3\tLN:4\n\
             @SQ\tSN:chr4\tLN:4\tUR:{unindexed}\n\
             @SQ\tSN:chr5\tLN:4\tUR:https://example.invalid/a.fa\n",
            a = a.display(),
            b = b.display(),
            unindexed = unindexed.display(),
        )
        .parse()?;

        let repository = header_repository(&header)?.expect("UR files");
        let sequence = |name: &[u8]| -> anyhow::Result<Option<Vec<u8>>> {
            Ok(repository
                .get(name)
                .transpose()?
                .map(|sequence| AsRef::<[u8]>::as_ref(&*sequence).to_vec()))
        };

        // chr2 comes from its own file, though the first file has a chr2 too, and chr3, which
        // names no file, from the file that has it.
        assert_eq!(sequence(b"chr1")?.as_deref(), Some(&b"ACGTACGT"[..]));
        assert_eq!(sequence(b"chr2")?.as_deref(), Some(&b"CCCCAAAA"[..]));
        assert_eq!(sequence(b"chr3")?.as_deref(), Some(&b"GGGG"[..]));
        assert_eq!(sequence(b"chr4")?, None);
        assert_eq!(sequence(b"chr5")?, None);

        let header: sam::Header = "@SQ\tSN:chr1\tLN:8\n".parse()?;
        assert!(header_repository(&header)?.is_none());

        Ok(())
    }

    #[test]
    fn test_cram_reader_decodes_sequences_from_several_ur_files() -> anyhow::Result<()> {
        const CHR1: &str = "ACGTTGCAAGGCTTAACCGT";
        const CHR2: &str = "TTGGCCAATTGGCCAAGGTT";

        let dir = tempfile::tempdir()?;
        let a = dir.path().join("a.fa");
        let b = dir.path().join("b.fa");
        let both = dir.path().join("both.fa");
        write_indexed_fasta(&a, &[("chr1", CHR1)])?;
        write_indexed_fasta(&b, &[("chr2", CHR2)])?;
        write_indexed_fasta(&both, &[("chr1", CHR1), ("chr2", CHR2)])?;

        let header: sam::Header = format!(
            "@SQ\tSN:chr1\tLN:20\tUR:file://{}\n@SQ\tSN:chr2\tLN:20\tUR:{}\n",
            a.display(),
            b.display(),
        )
        .parse()?;

        // Reads that match the reference, so that CRAM stores none of their bases.
        let lines = [
            format!("r1\t0\tchr1\t3\t60\t10M\t*\t0\t0\t{}\t*", &CHR1[2..12]),
            format!("r2\t0\tchr2\t5\t60\t10M\t*\t0\t0\t{}\t*", &CHR2[4..14]),
        ];

        let cram = dir.path().join("in.cram");
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(reference_repository(Some(&both))?)
            .build_from_writer(File::create(&cram)?);
        writer.write_header(&header)?;

        for line in &lines {
            let record = sam::Record::try_from(line.as_bytes())?;
            writer.write_alignment_record(&header, &record)?;
        }

        writer.try_finish(&header)?;

        let (mut reader, _) =
            AlignmentReader::open(cram.to_str().expect("UTF-8 path"), &ReadOptions::default())?;
        let header = reader.read_header()?;
        let mut sequences = Vec::new();

        reader.for_each(&header, |record| {
            sequences.push(record.as_dyn().sequence().iter().collect::<Vec<u8>>());
            Ok(())
        })?;

        assert_eq!(
            sequences,
            [&CHR1.as_bytes()[2..12], &CHR2.as_bytes()[4..14]]
        );

        Ok(())
    }
}
