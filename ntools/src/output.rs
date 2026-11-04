//! Alignment writers for SAM, BAM, and CRAM output to a file or stdout.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::Context;
use pasta::bam;
use pasta::bgzf::io::writer::CompressionLevel;
use pasta::cram;
use pasta::sam::{self as sam, alignment::io::Write as _};

use crate::{
    alignment::{Rec, header_repository, reference_repository},
    bam_sam::{self, Names},
    bgzf_stream::BgzfWriter,
    parallel_sam::ParallelSamWriter,
    runs::write_raw_record,
};

const SAM_BUFFER_CAPACITY: usize = 1 << 17;

/// An output alignment format.
#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
pub enum OutputFormat {
    Sam,
    /// SAM compressed with BGZF (`-O sam.gz`, or an output path ending in `.sam.gz`).
    #[value(name = "sam.gz")]
    SamGz,
    Bam,
    Cram,
}

impl OutputFormat {
    /// Infers a format from an output path's extension, as HTSlib does (`find_file_extension`
    /// and `sam_open_mode`), ignoring any `##idx##INDEX` suffix: `sam`, `bam`, `cram`, or
    /// `sam.gz`, in any case, except that the `gz` of `sam.gz` must be lowercase.
    pub fn from_path(path: &Path) -> Option<Self> {
        let (path, _) = split_index_path(path);
        let path = path.as_os_str().as_encoded_bytes();

        // The last `.` of the file name, as HTSlib scans back for `.` and stops at `/`.
        let dot = |end: usize| {
            let i = path[..end].iter().rposition(|&c| c == b'.' || c == b'/')?;
            (path[i] == b'.').then_some(i)
        };

        let mut start = dot(path.len())?;

        // HTSlib takes `.gz` and `.bgz` with the extension before them.
        if matches!(&path[start + 1..], b"gz" | b"bgz") {
            start = dot(start)?;
        }

        let extension = &path[start + 1..];

        [
            (&b"sam"[..], Self::Sam),
            (b"bam", Self::Bam),
            (b"cram", Self::Cram),
            (b"sam.gz", Self::SamGz),
        ]
        .into_iter()
        .find_map(|(name, format)| extension.eq_ignore_ascii_case(name).then_some(format))
    }

    /// Whether the output is SAM text, compressed or not.
    pub fn is_sam(self) -> bool {
        matches!(self, Self::Sam | Self::SamGz)
    }

    /// The file name extension of this format, for output names made from a prefix (`collate`,
    /// `split`). samtools (`hts_format_file_extension`) names BGZF SAM `.sam` too, a name that
    /// hides its compression, so ntools uses `.sam.gz` (README).
    pub fn extension(self) -> &'static str {
        match self {
            Self::Sam => "sam",
            Self::SamGz => "sam.gz",
            Self::Bam => "bam",
            Self::Cram => "cram",
        }
    }
}

/// The delimiter between an output file name and its index file name (`HTS_IDX_DELIM`).
const INDEX_DELIMITER: &str = "##idx##";

/// Splits `FILE##idx##INDEX` into the file and the index; other paths have no index.
///
/// As in HTSlib, which strips the suffix whenever it opens a file, the writers here write to the
/// file part, and `--write-index` writes the index part.
pub fn split_index_path(path: &Path) -> (PathBuf, Option<PathBuf>) {
    match path.to_str().and_then(|s| s.split_once(INDEX_DELIMITER)) {
        Some((file, index)) => (PathBuf::from(file), Some(PathBuf::from(index))),
        None => (path.to_path_buf(), None),
    }
}

/// Options for creating an alignment writer.
#[derive(Clone, Debug)]
pub struct WriteOptions {
    pub format: OutputFormat,
    /// The BGZF compression level for BAM output.
    pub level: CompressionLevel,
    /// Additional threads for BGZF compression, for formatting SAM, or for encoding CRAM.
    pub threads: usize,
    /// A FASTA reference for CRAM encoding.
    pub reference: Option<PathBuf>,
    /// CRAM version, compression level, and methods.
    pub cram: cram::io::writer::compression::Options,
}

pub type Sink = Box<dyn Write + Send>;

/// Where SAM text goes: the output itself, or a BGZF writer to it.
pub enum SamSink {
    Plain(Sink),
    Bgzf(BgzfWriter<Sink>),
}

impl SamSink {
    /// Flushes the output and, for BGZF, writes the last block and the EOF marker.
    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Bgzf(w) => w.finish()?.flush(),
        }
    }
}

impl Write for SamSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Bgzf(w) => w.write(buf),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.write_all(buf),
            Self::Bgzf(w) => w.write_all(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Bgzf(w) => w.flush(),
        }
    }
}

/// An open alignment writer.
///
/// An output file it creates is removed if the writer is dropped before writing a record or
/// finishing, e.g., when the input or the reference fails: samtools, which checks them before
/// opening its output, then leaves no file either. A file that already existed (including a FIFO
/// or a device) is never removed, nor is stdout.
pub struct AlignmentWriter {
    writer: FormatWriter,
    /// Dropped after `writer`, so that a CRAM writer finishing on drop writes before removal.
    /// Boxed, as the path would make every writer bigger.
    created: Option<Box<CreatedFile>>,
}

/// An output file this process created, removed on drop unless kept.
struct CreatedFile {
    path: PathBuf,
    keep: bool,
}

impl CreatedFile {
    /// Opens `path` for writing, truncating it, and returns it with a guard if this created it.
    fn open(path: &Path) -> io::Result<(File, Option<Box<Self>>)> {
        match File::options().write(true).create_new(true).open(path) {
            Ok(file) => {
                let created = Box::new(Self {
                    path: path.to_path_buf(),
                    keep: false,
                });
                Ok((file, Some(created)))
            }
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => Ok((File::create(path)?, None)),
            Err(e) => Err(e),
        }
    }
}

impl Drop for CreatedFile {
    fn drop(&mut self) {
        if !self.keep {
            // Best effort: the command is already failing.
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

impl AlignmentWriter {
    /// Creates a writer to `dst` (without any `##idx##INDEX` suffix), or to stdout when `dst` is
    /// `None` or `-`, for records with reference sequence IDs relative to `header`.
    ///
    /// As in HTSlib, SAM output to a path ending in `.sam.gz` is compressed with BGZF even when
    /// the format is given as `sam`.
    pub fn create(
        dst: Option<&Path>,
        options: &WriteOptions,
        header: &sam::Header,
    ) -> anyhow::Result<Self> {
        let (writer, created) = FormatWriter::create(dst, options, header)?;
        Ok(Self { writer, created })
    }

    /// Keeps the output file from now on.
    fn keep(&mut self) {
        if let Some(created) = &mut self.created {
            created.keep = true;
        }
    }

    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.writer.write_header(header)
    }

    /// See [`FormatWriter::write_header_text`].
    pub fn write_header_text(&mut self, header: &sam::Header, text: &[u8]) -> io::Result<()> {
        self.writer.write_header_text(header, text)
    }

    /// See [`FormatWriter::write`].
    pub fn write(&mut self, header: &sam::Header, record: &Rec<'_>) -> io::Result<()> {
        self.writer.write(header, record)?;
        self.keep();
        Ok(())
    }

    /// See [`FormatWriter::write_raw`].
    pub fn write_raw(&mut self, header: &sam::Header, record: &[u8]) -> io::Result<()> {
        self.writer.write_raw(header, record)?;
        self.keep();
        Ok(())
    }

    /// See [`FormatWriter::preformatted_kind`].
    pub fn preformatted_kind(&self) -> Option<OutputFormat> {
        self.writer.preformatted_kind()
    }

    /// See [`FormatWriter::write_sam_text`].
    pub fn write_sam_text(&mut self, text: Vec<u8>) -> io::Result<()> {
        self.writer.write_sam_text(text)?;
        self.keep();
        Ok(())
    }

    /// See [`FormatWriter::write_bam_records`].
    pub fn write_bam_records(&mut self, records: &[u8]) -> io::Result<()> {
        self.writer.write_bam_records(records)?;
        self.keep();
        Ok(())
    }

    /// Flushes buffered records and writes any end-of-file marker.
    pub fn finish(self, header: &sam::Header) -> io::Result<()> {
        let Self { writer, created } = self;
        writer.finish(header)?;

        if let Some(mut created) = created {
            created.keep = true;
        }

        Ok(())
    }
}

/// An open alignment writer of one format.
enum FormatWriter {
    Sam {
        writer: sam::io::Writer<BufWriter<SamSink>>,
        /// For formatting BAM records directly; see [`bam_sam`].
        names: Names,
        line: Vec<u8>,
    },
    /// SAM formatted on worker threads. Records that are not raw BAM are encoded first.
    ParallelSam {
        writer: ParallelSamWriter,
        encoder: bam::io::Writer<Vec<u8>>,
    },
    Bam(bam::io::Writer<BgzfWriter<Sink>>),
    Cram(Box<cram::io::Writer<Sink>>, CramReference),
}

/// The reference of CRAM output, for the `@SQ` fields HTSlib's CRAM writer fills in.
pub enum CramReference {
    /// A FASTA file given with `--reference`.
    Given(GivenReference),
    /// Without `--reference`, the local FASTA files the header's `@SQ` `UR` fields name, as
    /// HTSlib finds them (see [`header_repository`]); possibly none, in which case records are
    /// stored without an external reference, as in samtools.
    Header(HeaderReference),
}

impl CramReference {
    /// Returns `header` with the `@SQ` fields HTSlib fills in (see [`GivenReference::annotate`]
    /// and [`HeaderReference::annotate`]).
    fn annotate(&self, header: &sam::Header) -> io::Result<sam::Header> {
        match self {
            Self::Given(reference) => reference.annotate(header),
            Self::Header(reference) => reference.annotate(header),
        }
    }

    /// Warns once, as samtools does, for the first aligned record whose reference sequence is
    /// not available without `--reference`: the CRAM writer stores the bases of such records
    /// without an external reference (see the README), which makes the file bigger. (With
    /// `--reference`, [`GivenReference::annotate`] warns about the sequences it lacks.)
    fn warn_if_unavailable(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        match self {
            Self::Given(_) => Ok(()),
            Self::Header(reference) => reference.warn_if_unavailable(header, record),
        }
    }
}

/// See [`CramReference::Header`].
pub struct HeaderReference {
    repository: pasta::fasta::Repository,
    /// Whether each reference sequence (by header index) is in the repository, once looked up.
    available: Vec<Option<bool>>,
    /// Whether a record without a reference has been warned about.
    warned: bool,
}

impl HeaderReference {
    fn new(header: &sam::Header) -> anyhow::Result<Self> {
        Ok(Self {
            repository: header_repository(header)?.unwrap_or_default(),
            available: vec![None; header.reference_sequences().len()],
            warned: false,
        })
    }

    /// Returns `header` with `M5` fields, where missing, from the `UR` files that have the
    /// sequence, as HTSlib does.
    fn annotate(&self, header: &sam::Header) -> io::Result<sam::Header> {
        use sam::header::record::value::map::reference_sequence::tag;

        let mut header = header.clone();

        for (name, reference_sequence) in header.reference_sequences_mut() {
            let fields = reference_sequence.other_fields_mut();

            if !fields.contains_key(&tag::MD5_CHECKSUM)
                && let Some(checksum) = self.repository.normalized_md5(name).transpose()?
            {
                let hex: String = checksum.iter().map(|b| format!("{b:02x}")).collect();
                fields.insert(tag::MD5_CHECKSUM, hex.into());
            }
        }

        Ok(header)
    }

    fn warn_if_unavailable(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        if self.warned || record.flags()?.is_unmapped() || record.alignment_start().is_none() {
            return Ok(());
        }

        let Some(id) = record.reference_sequence_id(header).transpose()? else {
            return Ok(());
        };

        // An invalid ID fails in the writer.
        let (Some(available), Some((name, _))) = (
            self.available.get_mut(id),
            header.reference_sequences().get_index(id),
        ) else {
            return Ok(());
        };

        let available = match *available {
            Some(available) => available,
            None => *available.insert(self.repository.get(name).transpose()?.is_some()),
        };

        if !available {
            eprintln!(
                "[W::cram_get_ref] No reference for {name} (no -T/--reference, and no @SQ UR file \
                 has it); the CRAM file will be bigger than using an external reference"
            );
            self.warned = true;
        }

        Ok(())
    }
}

/// See [`CramReference::Given`].
pub struct GivenReference {
    /// The reference FASTA's absolute path (`UR`), as HTSlib makes it (without resolving links).
    uri: String,
    /// The reference FASTA, read for checksums.
    path: PathBuf,
    /// Whether `M5` fields are replaced with the reference's checksums (`ignore_md5`).
    replaces_md5: bool,
}

impl GivenReference {
    /// Returns `header` with HTSlib's `@SQ` fields: `M5`, when missing, from the reference (if
    /// it has the sequence), and `UR`, the reference's absolute path, replacing any.
    ///
    /// As HTSlib does, this warns about sequences without `M5` that the reference lacks (the
    /// writer stores their records' bases verbatim) and about `LN` fields that differ from the
    /// reference's lengths. An `M5` that does not match the reference fails the records on its
    /// sequence (see `cram::io::Writer`), unless `ignore_md5` is set: the `M5` fields of the
    /// sequences the reference has are then replaced with its checksums, so that the file can be
    /// read back with this reference. (HTSlib keeps them, writing a header that does not
    /// describe the reference the records are stored against.)
    ///
    /// Each sequence is read only for its checksum and not kept (the writer's repository loads
    /// the ones records need), so one sequence at a time is in memory.
    fn annotate(&self, header: &sam::Header) -> io::Result<sam::Header> {
        use pasta::{core::Region, fasta};
        use sam::header::record::value::map::reference_sequence::tag;

        let mut header = header.clone();
        let mut missing = false;
        let mut reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(&self.path)?;

        for (name, reference_sequence) in header.reference_sequences_mut() {
            let length = reader.index().get(name).map(|record| record.length());

            if let Some(length) = length
                && u64::try_from(reference_sequence.length().get()).ok() != Some(length)
            {
                eprintln!(
                    "[W::sanitise_SQ_lines] Header @SQ length mismatch for ref {name}, {} vs {length}",
                    reference_sequence.length()
                );
            }

            let fields = reference_sequence.other_fields_mut();

            // A replaced checksum keeps its place.
            if !fields.contains_key(&tag::MD5_CHECKSUM) || (self.replaces_md5 && length.is_some()) {
                if length.is_some() {
                    let record = reader.query(&Region::new(name.to_vec(), ..))?;
                    let sequence = record.sequence().as_ref();
                    let checksum = cram::calculate_normalized_sequence_digest(sequence);
                    let hex: String = checksum.iter().map(|b| format!("{b:02x}")).collect();
                    fields.insert(tag::MD5_CHECKSUM, hex.into());
                } else {
                    eprintln!(
                        "[W::cram_get_ref] Reference file given, but ref '{name}' not present"
                    );
                    missing = true;
                }
            }

            fields.insert(tag::URI, self.uri.clone().into());
        }

        if missing {
            eprintln!(
                "[W::cram_write_SAM_hdr] NOTE: the CRAM file will be bigger than using an external reference"
            );
        }

        Ok(header)
    }
}

impl FormatWriter {
    /// See [`AlignmentWriter::create`]. This also returns the guard of the output file, if it
    /// created one, which removes the file if this fails.
    fn create(
        dst: Option<&Path>,
        options: &WriteOptions,
        header: &sam::Header,
    ) -> anyhow::Result<(Self, Option<Box<CreatedFile>>)> {
        let format = match options.format {
            OutputFormat::Sam
                if dst.and_then(OutputFormat::from_path) == Some(OutputFormat::SamGz) =>
            {
                OutputFormat::SamGz
            }
            format => format,
        };

        let dst = dst.map(|path| split_index_path(path).0);

        let (sink, created): (Sink, _) = match dst.as_deref() {
            Some(path) if path != Path::new("-") => {
                pasta::util::inputs::check_output(path)
                    .with_context(|| format!("creating {}", path.display()))?;
                let (file, created) = CreatedFile::open(path)
                    .with_context(|| format!("creating {}", path.display()))?;
                (Box::new(file), created)
            }
            _ => (raw_stdout()?, None),
        };

        let writer = match format {
            OutputFormat::Sam if options.threads > 0 => Self::ParallelSam {
                writer: ParallelSamWriter::new(sink, header, options.threads),
                encoder: bam::io::Writer::from(Vec::new()),
            },
            OutputFormat::Sam | OutputFormat::SamGz => {
                // BGZF SAM is formatted on this thread and compressed on the `threads` workers.
                let sink = match format {
                    OutputFormat::SamGz => {
                        SamSink::Bgzf(BgzfWriter::new(sink, options.level, options.threads))
                    }
                    _ => SamSink::Plain(sink),
                };

                Self::Sam {
                    writer: sam::io::Writer::new(BufWriter::with_capacity(
                        SAM_BUFFER_CAPACITY,
                        sink,
                    )),
                    names: Names::new(header),
                    line: Vec::new(),
                }
            }
            OutputFormat::Bam => Self::Bam(bam::io::Writer::from(BgzfWriter::new(
                sink,
                options.level,
                options.threads,
            ))),
            OutputFormat::Cram => {
                let (repository, reference) = match &options.reference {
                    Some(path) => (
                        reference_repository(Some(path))?,
                        CramReference::Given(GivenReference {
                            uri: std::path::absolute(path)
                                .unwrap_or_default()
                                .to_string_lossy()
                                .into_owned(),
                            path: path.clone(),
                            replaces_md5: options.cram.ignore_md5,
                        }),
                    ),
                    None => {
                        let reference = HeaderReference::new(header)?;
                        (
                            reference.repository.clone(),
                            CramReference::Header(reference),
                        )
                    }
                };

                let mut builder = cram::io::writer::Builder::default()
                    .set_reference_sequence_repository(repository)
                    .set_compression_options(options.cram);

                if let Some(worker_count) = NonZero::new(options.threads) {
                    builder = builder.set_worker_count(worker_count);
                }

                Self::Cram(Box::new(builder.build_from_writer(sink)), reference)
            }
        };

        Ok((writer, created))
    }

    fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam { writer, .. } => writer.write_header(header),
            Self::ParallelSam { writer, .. } => {
                let mut text = sam::io::Writer::new(Vec::new());
                text.write_header(header)?;
                writer.write_text(text.into_inner())
            }
            Self::Bam(w) => w.write_header(header),
            Self::Cram(w, reference) => w.write_header(&reference.annotate(header)?),
        }
    }

    /// Writes a header as the given SAM text, keeping its line order, as HTSlib writes a header it
    /// built or read. `header` must be `text` parsed; CRAM output writes that instead.
    fn write_header_text(&mut self, header: &sam::Header, text: &[u8]) -> io::Result<()> {
        match self {
            Self::Sam { writer, .. } => writer.get_mut().write_all(text),
            Self::ParallelSam { writer, .. } => writer.write_text(text.to_vec()),
            Self::Bam(w) => {
                let dst = w.get_mut();
                let len = |n: usize| {
                    u32::try_from(n)
                        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "header too long"))
                };

                dst.write_all(b"BAM\x01")?;
                dst.write_all(&len(text.len())?.to_le_bytes())?;
                dst.write_all(text)?;

                let reference_sequences = header.reference_sequences();
                dst.write_all(&len(reference_sequences.len())?.to_le_bytes())?;

                for (name, reference_sequence) in reference_sequences {
                    dst.write_all(&len(name.len() + 1)?.to_le_bytes())?;
                    dst.write_all(name)?;
                    dst.write_all(&[0])?;
                    dst.write_all(&len(reference_sequence.length().get())?.to_le_bytes())?;
                }

                Ok(())
            }
            Self::Cram(w, reference) => w.write_header(&reference.annotate(header)?),
        }
    }

    /// Writes a record. BAM records going to BAM output are forwarded without re-encoding, so
    /// `header` must have the same reference sequences as the header the record was read with.
    fn write(&mut self, header: &sam::Header, record: &Rec<'_>) -> io::Result<()> {
        match (self, record) {
            (Self::Bam(w), Rec::Bam(r)) => write_raw_record(w.get_mut(), r.as_ref()),
            (Self::ParallelSam { writer, .. }, Rec::Bam(r)) => writer.push(r.as_ref()),
            (Self::ParallelSam { writer, encoder }, _) => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(header, record.as_dyn())?;
                writer.push(&encoder.get_ref()[size_of::<u32>()..])
            }
            (
                Self::Sam {
                    writer,
                    names,
                    line,
                },
                Rec::Bam(r),
            ) => {
                line.clear();

                if bam_sam::format_record(names, r.as_ref(), line) {
                    writer.get_mut().write_all(line)
                } else {
                    writer.write_alignment_record(header, record.as_dyn())
                }
            }
            (Self::Sam { writer, .. }, _) => writer.write_alignment_record(header, record.as_dyn()),
            (Self::Bam(w), _) => w.write_alignment_record(header, record.as_dyn()),
            (Self::Cram(w, reference), Rec::Bam(r)) => {
                reference.warn_if_unavailable(header, record.as_dyn())?;
                w.write_bam_record(header, r)
            }
            (Self::Cram(w, reference), _) => {
                reference.warn_if_unavailable(header, record.as_dyn())?;
                w.write_alignment_record(header, record.as_dyn())
            }
        }
    }

    /// Writes an encoded BAM record (without `block_size`) whose reference sequence IDs are
    /// relative to `header`. BAM output copies it; other formats decode it first.
    fn write_raw(&mut self, header: &sam::Header, record: &[u8]) -> io::Result<()> {
        match self {
            Self::Bam(w) => write_raw_record(w.get_mut(), record),
            Self::ParallelSam { writer, .. } => writer.push(record),
            _ => {
                let record = bam::Record::try_from(record.to_vec())?;
                self.write(header, &Rec::Bam(&record))
            }
        }
    }

    /// Returns whether pre-formatted SAM text ([`Self::write_sam_text`]) or pre-encoded BAM
    /// records ([`Self::write_bam_records`]) can be written directly, and which.
    fn preformatted_kind(&self) -> Option<OutputFormat> {
        match self {
            Self::Sam { .. } | Self::ParallelSam { .. } => Some(OutputFormat::Sam),
            Self::Bam(_) => Some(OutputFormat::Bam),
            Self::Cram(..) => None,
        }
    }

    /// Writes SAM record lines formatted elsewhere. Only valid for SAM output.
    fn write_sam_text(&mut self, text: Vec<u8>) -> io::Result<()> {
        match self {
            Self::Sam { writer, .. } => writer.get_mut().write_all(&text),
            Self::ParallelSam { writer, .. } => writer.write_text(text),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "SAM text written to non-SAM output",
            )),
        }
    }

    /// Writes `block_size`-prefixed BAM records encoded elsewhere. Only valid for BAM output.
    fn write_bam_records(&mut self, records: &[u8]) -> io::Result<()> {
        match self {
            Self::Bam(w) => w.get_mut().write_all(records),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM records written to non-BAM output",
            )),
        }
    }

    /// Flushes buffered records and writes any end-of-file marker.
    fn finish(self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam { writer, .. } => writer
                .into_inner()
                .into_inner()
                .map_err(io::IntoInnerError::into_error)?
                .finish(),
            Self::ParallelSam { writer, .. } => writer.finish(),
            Self::Bam(w) => w.into_inner().finish()?.flush(),
            Self::Cram(mut w, _) => {
                w.finish(header)?;
                w.get_mut().flush()
            }
        }
    }
}

/// Returns stdout without Rust's line buffering, which scans every write for newlines; callers
/// buffer output themselves.
#[cfg(unix)]
pub fn raw_stdout() -> anyhow::Result<Sink> {
    use std::os::fd::AsFd;

    let fd = io::stdout()
        .as_fd()
        .try_clone_to_owned()
        .context("duplicating stdout")?;

    Ok(Box::new(File::from(fd)))
}

#[cfg(not(unix))]
pub fn raw_stdout() -> anyhow::Result<Sink> {
    Ok(Box::new(io::stdout()))
}

/// Parses a samtools compression level (`0`–`9`).
pub fn parse_level(level: u8) -> anyhow::Result<CompressionLevel> {
    CompressionLevel::new(level).with_context(|| format!("invalid compression level: {level}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_index_path() {
        assert_eq!(
            split_index_path(Path::new("out.bam##idx##out.custom.bai")),
            (
                PathBuf::from("out.bam"),
                Some(PathBuf::from("out.custom.bai"))
            )
        );
        assert_eq!(
            split_index_path(Path::new("out.bam")),
            (PathBuf::from("out.bam"), None)
        );
    }

    #[test]
    fn test_from_path() {
        assert_eq!(
            OutputFormat::from_path(Path::new("out.cram")),
            Some(OutputFormat::Cram)
        );
        assert_eq!(
            OutputFormat::from_path(Path::new("out.bam##idx##out.csi")),
            Some(OutputFormat::Bam)
        );
        assert_eq!(
            OutputFormat::from_path(Path::new("out##idx##out.sam")),
            None
        );

        // HTSlib's extensions: `.gz` and `.bgz` join the one before them, names compare in any
        // case but the `gz`, and the search stops at `/`.
        for (path, expected) in [
            ("out.sam.gz", Some(OutputFormat::SamGz)),
            ("out.SAM.gz", Some(OutputFormat::SamGz)),
            ("out.sam.gz##idx##out.csi", Some(OutputFormat::SamGz)),
            ("out.sam.GZ", None),
            ("out.sam.bgz", None),
            ("out.bam.gz", None),
            ("out.gz", None),
            ("out.BAM", Some(OutputFormat::Bam)),
            ("out.Sam", Some(OutputFormat::Sam)),
            ("dir.sam/out", None),
            ("dir.sam/out.gz", None),
            (".sam", Some(OutputFormat::Sam)),
        ] {
            assert_eq!(OutputFormat::from_path(Path::new(path)), expected, "{path}");
        }
    }
}
