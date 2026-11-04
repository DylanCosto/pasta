//! Reading several VCF/BCF files in step, position by position: a port of HTSlib's synced
//! reader (`synced_bcf_reader.c` and `bcf_sr_sort.c`).
//!
//! [`SyncedReaders`] is HTSlib's `bcf_srs_t`. It holds any number of readers and returns, on
//! each [`SyncedReaders::next_line`], the next position's records: one record per reader at
//! most, paired by their alleles according to [`Pairing`]. Records at one position that do not
//! pair are returned on successive calls, in HTSlib's order. Commands that combine files
//! (`isec`, `merge`, `annotate`, `stats`, `gtcheck`, `concat`) are built on it, and so are the
//! single-file commands, through [`SyncedReader`].
//!
//! # Using several readers
//!
//! Set the options first, then add the readers, then iterate:
//!
//! ```no_run
//! use pasta::variant::synced::{IndexRequirement, Pairing, SyncedReaders};
//!
//! # fn main() -> Result<(), pasta::variant::synced::SyncedError> {
//! let mut sr = SyncedReaders::new();
//! sr.set_pairing(Pairing::BOTH_REF); // bcf_sr_set_opt(sr, BCF_SR_PAIR_LOGIC, ...)
//! sr.set_index_requirement(IndexRequirement::Required); // BCF_SR_REQUIRE_IDX
//! for path in ["a.vcf.gz", "b.bcf"] {
//!     sr.add_reader(path)?; // bcf_sr_add_reader
//! }
//! while sr.next_line()? > 0 {
//!     for i in 0..sr.reader_count() {
//!         if let Some(rec) = sr.line(i) {
//!             // `rec` is reader i's record at this position (bcf_sr_get_line).
//!             let _ = (rec.pos(), sr.header(i));
//!         }
//!     }
//! }
//! # Ok(())
//! # }
//! ```
//!
//! The HTSlib calls map as follows:
//!
//! | HTSlib | here |
//! |---|---|
//! | `bcf_sr_init` | [`SyncedReaders::new`] |
//! | `bcf_sr_set_opt(BCF_SR_PAIR_LOGIC)`, `sr->collapse` | [`SyncedReaders::set_pairing`] |
//! | `bcf_sr_set_opt(BCF_SR_REQUIRE_IDX / ALLOW_NO_IDX)`, `sr->require_index` | [`SyncedReaders::set_index_requirement`] |
//! | `bcf_sr_set_opt(BCF_SR_REGIONS_OVERLAP / TARGETS_OVERLAP)` | [`SyncedReaders::set_regions_overlap`], [`SyncedReaders::set_targets_overlap`] |
//! | `bcf_sr_set_regions` | [`SyncedReaders::set_regions`], [`SyncedReaders::set_regions_spec`] |
//! | `bcf_sr_set_targets` | [`SyncedReaders::set_targets`], [`SyncedReaders::set_targets_spec`] |
//! | `sr->apply_filters` | [`SyncedReaders::set_apply_filters`] |
//! | `sr->max_unpack` | [`SyncedReaders::set_max_unpack`], or per reader [`SyncedReaders::set_reader_max_unpack`] |
//! | `bcf_sr_set_threads` | [`SyncedReaders::set_threads`] |
//! | `bcf_sr_add_reader` | [`SyncedReaders::add_reader`] |
//! | `bcf_sr_remove_reader` | [`SyncedReaders::remove_reader`] |
//! | `bcf_sr_set_samples`, `sr->samples`, `reader->samples` | [`SyncedReaders::set_samples`], [`SyncedReaders::samples`], [`SyncedReaders::sample_map`] |
//! | `bcf_sr_next_line` | [`SyncedReaders::next_line`] |
//! | `bcf_sr_has_line`, `bcf_sr_get_line` | [`SyncedReaders::has_line`], [`SyncedReaders::line`], [`SyncedReaders::line_mut`] |
//! | `sr->has_line[i] = 0` | [`SyncedReaders::clear_line`] |
//! | `bcf_sr_swap_line` | [`SyncedReaders::swap_line`] |
//! | `bcf_sr_region_done` | [`SyncedReaders::region_done`] |
//! | `bcf_sr_get_header` | [`SyncedReaders::header`], [`SyncedReaders::header_mut`] |
//! | `reader->fname` | [`SyncedReaders::file_name`] |
//! | `reader->buffer[1..=nbuffer]` | [`SyncedReaders::buffered`] |
//! | `reader->buffer[j]`, `reader->nbuffer`, dropping used slots | [`SyncedReaders::buffer_record`], [`SyncedReaders::buffer_record_mut`], [`SyncedReaders::buffered_count`], [`SyncedReaders::discard_buffered`] |
//! | `bcf_sr_seek` | [`SyncedReaders::seek`], [`SyncedReaders::seek_start`] |
//! | `sr->regions`, `sr->targets` | [`SyncedReaders::regions`], [`SyncedReaders::targets_mut`] |
//! | `sr->errnum` after the loop | [`SyncedReaders::take_pending_error`] |
//!
//! Errors that end `bcf_sr_next_line` (`errnum` set and 0 returned) are returned as
//! [`SyncedError`]s; [`SyncedError::code`] gives HTSlib's `bcf_sr_error` number, which
//! bcftools uses as an exit status, and its [`Display`](std::fmt::Display) is
//! `bcf_sr_strerror`.
//!
//! # Index use
//!
//! With [`IndexRequirement::Required`] (set implicitly by regions), every input must be
//! BGZF-compressed and indexed; the readers then visit the regions, or when none are given,
//! every sequence in the indexes (in the order the readers list them), through index jumps.
//! Without it, the readers stream, which HTSlib allows for one reader, or for several with
//! [`IndexRequirement::AllowMissing`] when they have the same sequences in the same order.
//! Targets (`-t`/`-T`) never need an index: they filter records as they stream past.
//!
//! # Pairing
//!
//! With one reader, records come back one by one in file order. With several, all records at
//! the current position are grouped as HTSlib's `bcf_sr_sort_set` does: each reader's records
//! form a group keyed by their allele strings (`REF>ALT`, with `/END` for symbolic alleles and
//! `ID:` with [`Pairing::ID`]); identical records pair; with [`Pairing::EXACT`] (the default)
//! only records with the same alleles (in any order) pair; otherwise records sharing an allele
//! pair first, then records whose types the pairing allows (SNPs with SNPs, and so on).

mod regions;
mod sort;

use std::{
    io,
    path::{Path, PathBuf},
    sync::atomic::AtomicBool,
};

use self::sort::Sorter;
#[cfg(test)]
use self::sort::rec_id;
pub use self::{
    regions::{MAX_CSI_COOR, MissedRegion, Overlap, Regions, RegionsError},
    sort::Pairing,
};
use crate::variant::{
    Header, Record, Unpack,
    index::{
        Index, IndexFormat, locate, report_load_error, split_idx_name, tabix::Tabix, warn_if_older,
    },
    io::{
        Compression, InputFormat, Query, ReadError, Reader, has_bgzf_eof,
        reader::exec_format_error, stdin_has_bgzf_eof,
    },
    log,
    readlist::{self, ListSource},
};

/// An error from the synced reader, with HTSlib's `bcf_sr_error` codes as exit statuses.
#[derive(Debug)]
#[non_exhaustive]
pub enum SyncedError {
    /// The input could not be opened (`bcf_sr_open_failed`, 11).
    Open(io::Error),
    /// Index jumps need a BGZF file (`not_bgzf`, 1).
    NotBgzf,
    /// Index jumps need an index (`idx_load_failed`, 2).
    IndexLoad,
    /// The input is not VCF/BCF (`file_type_error`, 3).
    FileType,
    /// The reader was set up in a way HTSlib refuses (`api_usage_error`, 4).
    ApiUsage,
    /// The header could not be read (`header_error`, 5).
    Header,
    /// The BGZF EOF marker is missing (`no_eof`, 6).
    NoEof,
    /// A VCF line could not be parsed (`vcf_parse_error`, 8).
    VcfParse,
    /// A BCF record could not be read (`bcf_read_error`, 9).
    BcfRead,
    /// Several unindexed files cannot be read together (`noidx_error`, 10).
    NoIndex,
    /// A region could not be sought (`bcf_sr_seek_error`, 12).
    Seek,
    /// The regions or targets could not be read (`bcf_sr_regions_error`, 13).
    Regions,
    /// The sample list could not be read (`bcf_sr_samples_error`, 14).
    Samples,
    /// A line of a tabix-indexed targets file could not be read (a truncated or corrupt BGZF
    /// block, say). HTSlib reports it, counts the record as outside the targets, and goes on,
    /// so the records after it are dropped; this error stops instead, with status 255 (that of
    /// bcftools' `error`), as HTSlib has no code for it.
    TargetsRead(RegionsError),
}

impl SyncedError {
    /// Returns HTSlib's `bcf_sr_error` value, which bcftools uses as its exit status (255 for
    /// [`Self::TargetsRead`], which has none).
    #[must_use]
    pub fn code(&self) -> i32 {
        match self {
            Self::NotBgzf => 1,
            Self::IndexLoad => 2,
            Self::FileType => 3,
            Self::ApiUsage => 4,
            Self::Header => 5,
            Self::NoEof => 6,
            Self::VcfParse => 8,
            Self::BcfRead => 9,
            Self::NoIndex => 10,
            Self::Open(_) => 11,
            Self::Seek => 12,
            Self::Regions => 13,
            Self::Samples => 14,
            Self::TargetsRead(_) => 255,
        }
    }
}

impl std::fmt::Display for SyncedError {
    /// HTSlib's `bcf_sr_strerror`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // C's strerror, without Rust's " (os error N)".
            Self::Open(e) => {
                let s = e.to_string();
                f.write_str(s.split(" (os error").next().unwrap_or(&s))
            }
            Self::NotBgzf => f.write_str("not compressed with bgzip"),
            Self::IndexLoad => f.write_str("could not load index"),
            Self::FileType => f.write_str("unknown file type"),
            Self::ApiUsage => f.write_str("API usage error"),
            Self::Header => f.write_str("could not parse header"),
            Self::NoEof => f.write_str("no BGZF EOF marker; file may be truncated"),
            Self::VcfParse => f.write_str("VCF parse error"),
            Self::BcfRead => f.write_str("BCF read error"),
            Self::NoIndex => f.write_str("merge of unindexed files failed"),
            Self::Seek => f.write_str("failed to seek to next location"),
            Self::Regions => f.write_str("failed to read regions/targets list"),
            Self::Samples => f.write_str("failed to read samples list"),
            Self::TargetsRead(e) => write!(f, "failed to read the targets: {e}"),
        }
    }
}

impl std::error::Error for SyncedError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Open(e) => Some(e),
            Self::TargetsRead(e) => Some(e),
            _ => None,
        }
    }
}

/// Whether the inputs must be indexed (HTSlib's `require_index`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub enum IndexRequirement {
    /// Stream when there are no regions (the default, 0); only one reader can stream.
    #[default]
    Optional,
    /// Every input must be BGZF-compressed and indexed (`BCF_SR_REQUIRE_IDX`); the readers
    /// jump through the regions, or through every indexed sequence.
    Required,
    /// Stream several unindexed inputs, which must list the same sequences in the same order
    /// (`BCF_SR_ALLOW_NO_IDX`, at the user's risk).
    AllowMissing,
}

/// Whether records overlapping the targets are kept or dropped.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub enum TargetsLogic {
    /// Keep records overlapping the targets.
    #[default]
    Include,
    /// Drop records overlapping the targets (a `^` prefix).
    Exclude,
}

/// How a reader reads records.
enum Access {
    Stream,
    Tabix(Box<Tabix>),
    Bcf(Box<Index>),
}

/// One input (HTSlib's `bcf_sr_t`).
struct ReaderState {
    file: Reader,
    access: Access,
    header: Header,
    /// The index query of the current region (HTSlib's `itr`).
    query: Option<Query>,
    /// Whether a streamed input has reached its end.
    eof: bool,
    fname: String,
    /// Slot 0 holds the current record; slots 1..=nbuffer the records read ahead, those at
    /// one position and then one past it.
    #[expect(
        clippy::vec_box,
        reason = "boxes keep records in place, for cheap rotation and identity by address"
    )]
    buffer: Vec<Box<Record>>,
    nbuffer: usize,
    /// FILTER IDs to require, -1 for `.` (HTSlib's `filter_ids`).
    filter_ids: Option<Vec<i32>>,
    /// The header index of each shared sample (HTSlib's `samples`).
    samples: Vec<i32>,
    /// How much of each VCF line to parse, when set for this reader alone.
    max_unpack: Option<Unpack>,
}

impl ReaderState {
    /// Whether the region (or the stream) has nothing left to return.
    fn is_done(&self) -> bool {
        self.query.is_none() && self.nbuffer == 0
    }

    /// Drops the first buffered record (HTSlib's `_reader_shift_buffer`).
    fn shift_buffer(&mut self) {
        if self.nbuffer == 0 {
            return;
        }
        let n = self.nbuffer;
        rotate_left1(&mut self.buffer[1..=n]);
        self.nbuffer -= 1;
    }
}

/// `buffer.rotate_left(1)` for the few records of a reader's buffer: swaps, which inline,
/// where the general rotation is a call.
#[inline]
fn rotate_left1(buffer: &mut [Box<Record>]) {
    for k in 1..buffer.len() {
        buffer.swap(k - 1, k);
    }
}

/// Whether a BGZF input is known to lack the EOF marker block (HTSlib's `bgzf_check_EOF`).
fn lacks_bgzf_eof(file: &Reader) -> bool {
    file.compression() == Compression::Bgzf
        && match file.path() {
            Some(p) => has_bgzf_eof(Path::new(p)),
            None => stdin_has_bgzf_eof(),
        }
        .ok()
        .flatten()
            == Some(false)
}

/// Checks that an unindexed input's header lists the same sequences, in the same order, as
/// the first input's, as streaming several inputs requires.
fn check_same_sequences(hdr0: &Header, header: &Header) -> Result<(), SyncedError> {
    if hdr0.contig_count() != header.contig_count() {
        log::error(
            "bcf_sr_add_hreader",
            format_args!(
                "Different number of sequences in the header, refusing to stream multiple unindexed files"
            ),
        );
        return Err(SyncedError::NoIndex);
    }
    for i in 0..hdr0.contig_count() as i32 {
        if hdr0.contig_name(i) != header.contig_name(i) {
            log::error(
                "bcf_sr_add_hreader",
                format_args!(
                    "Sequences in the header appear in different order, refusing to stream multiple unindexed files"
                ),
            );
            return Err(SyncedError::NoIndex);
        }
    }
    Ok(())
}

/// Warned once per process, as HTSlib does.
static NO_INDEX_WARNED: AtomicBool = AtomicBool::new(false);

/// Several VCF/BCF inputs read in step (HTSlib's `bcf_srs_t`). See the [module
/// documentation](self) for an overview.
pub struct SyncedReaders {
    readers: Vec<ReaderState>,
    has_line: Vec<bool>,
    sort: Sorter,
    regions: Option<Regions>,
    targets: Option<Regions>,
    targets_logic: TargetsLogic,
    targets_alleles: Option<usize>,
    require_index: IndexRequirement,
    streaming: bool,
    explicit_regs: bool,
    apply_filters: Option<Vec<u8>>,
    regions_overlap: Overlap,
    targets_overlap: Overlap,
    max_unpack: Unpack,
    threads: usize,
    /// How many readers read in step share `threads` ([`Self::set_shared_threads`]).
    sharing: usize,
    /// The readers opened so far, which take the remainder of `threads / sharing` in turn.
    opened: usize,
    samples: Vec<Vec<u8>>,
    line: Vec<u8>,
    chr: Vec<u8>,
    pending_error: Option<SyncedError>,
}

impl std::fmt::Debug for SyncedReaders {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SyncedReaders")
            .field(
                "readers",
                &self.readers.iter().map(|r| &r.fname).collect::<Vec<_>>(),
            )
            .field("streaming", &self.streaming)
            .finish_non_exhaustive()
    }
}

impl Default for SyncedReaders {
    fn default() -> Self {
        Self::new()
    }
}

impl SyncedReaders {
    /// Creates a reader set with HTSlib's defaults (`bcf_sr_init`): exact pairing, regions
    /// matched by record overlap, targets by position.
    #[must_use]
    pub fn new() -> Self {
        Self {
            readers: Vec::new(),
            has_line: Vec::new(),
            sort: Sorter::default(),
            regions: None,
            targets: None,
            targets_logic: TargetsLogic::Include,
            targets_alleles: None,
            require_index: IndexRequirement::Optional,
            streaming: false,
            explicit_regs: false,
            apply_filters: None,
            regions_overlap: Overlap::Record,
            targets_overlap: Overlap::Pos,
            max_unpack: Unpack::empty(),
            threads: 0,
            sharing: 1,
            opened: 0,
            samples: Vec::new(),
            line: Vec::new(),
            chr: Vec::new(),
            pending_error: None,
        }
    }

    /// Sets which records at one position pair (`BCF_SR_PAIR_LOGIC`, or bcftools'
    /// `sr->collapse`). The empty set means [`Pairing::EXACT`]. Takes effect at the first
    /// position with records in several readers.
    pub fn set_pairing(&mut self, pairing: Pairing) {
        self.sort.pair = pairing;
    }

    /// Sets whether inputs must be indexed (`BCF_SR_REQUIRE_IDX`, `BCF_SR_ALLOW_NO_IDX`).
    /// Call before adding readers.
    pub fn set_index_requirement(&mut self, requirement: IndexRequirement) {
        self.require_index = requirement;
    }

    /// Sets when a record overlaps a region (`BCF_SR_REGIONS_OVERLAP`; record overlap by
    /// default).
    pub fn set_regions_overlap(&mut self, overlap: Overlap) {
        self.regions_overlap = overlap;
        if let Some(r) = &mut self.regions {
            r.overlap = overlap;
        }
    }

    /// Sets when a record overlaps a target (`BCF_SR_TARGETS_OVERLAP`; position by default).
    pub fn set_targets_overlap(&mut self, overlap: Overlap) {
        self.targets_overlap = overlap;
        if let Some(t) = &mut self.targets {
            t.overlap = overlap;
        }
    }

    /// Requires records to have one of these comma-separated FILTER names, `.` standing for
    /// none (`sr->apply_filters`). Applies to readers added afterwards.
    pub fn set_apply_filters(&mut self, filters: Option<&[u8]>) {
        self.apply_filters = filters.map(<[u8]>::to_vec);
    }

    /// Limits how much of each VCF line is parsed (`sr->max_unpack`). Applies to records
    /// read afterwards.
    pub fn set_max_unpack(&mut self, which: Unpack) {
        self.max_unpack = which;
    }

    /// Limits how much of reader `i`'s VCF lines is parsed, overriding
    /// [`Self::set_max_unpack`] for that reader: for a reader whose records are only paired
    /// and never written, say.
    pub fn set_reader_max_unpack(&mut self, i: usize, which: Unpack) {
        let reader = &mut self.readers[i];
        reader.max_unpack = Some(which);
        for rec in &mut reader.buffer {
            rec.set_max_unpack(which);
        }
    }

    /// Sets the number of BGZF decompression threads of each reader added afterwards
    /// (`bcf_sr_set_threads`).
    pub fn set_threads(&mut self, threads: usize) {
        self.set_shared_threads(threads, 1);
    }

    /// Shares `threads` BGZF decompression threads among the `readers` readers read in step
    /// that are added afterwards, so that many inputs do not start `threads` each: each reader
    /// gets `threads / readers`, and the remainder goes one each to the readers added first.
    /// A reader without a thread decompresses on the calling thread.
    pub fn set_shared_threads(&mut self, threads: usize, readers: usize) {
        self.threads = threads;
        self.sharing = readers.max(1);
        self.opened = 0;
    }

    /// The decompression threads of the next reader added ([`Self::set_shared_threads`]).
    fn next_reader_threads(&mut self) -> usize {
        let extra = self.opened % self.sharing < self.threads % self.sharing;
        self.opened += 1;
        self.threads / self.sharing + usize::from(extra)
    }

    /// Sets the regions to visit through the indexes (`bcf_sr_set_regions`). Before any
    /// reader is added, this also requires indexes; afterwards it replaces the regions and
    /// starts over.
    pub fn set_regions(&mut self, mut regions: Regions) {
        regions.overlap = self.regions_overlap;
        if !self.readers.is_empty() || self.regions.is_some() {
            self.regions = Some(regions);
            self.seek_start();
            return;
        }
        self.regions = Some(regions);
        self.explicit_regs = true;
        self.require_index = IndexRequirement::Required;
    }

    /// Reads regions from a comma-separated list or a file and sets them
    /// ([`Self::set_regions`]).
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::Regions`] when the list cannot be read.
    pub fn set_regions_spec(&mut self, spec: &[u8], source: ListSource) -> Result<(), SyncedError> {
        let regions = Regions::new(spec, source).map_err(|_| SyncedError::Regions)?;
        self.set_regions(regions);
        Ok(())
    }

    /// Refuses targets once readers are added or targets set, as HTSlib does.
    fn check_targets_unset(&self) -> Result<(), SyncedError> {
        if !self.readers.is_empty() || self.targets.is_some() {
            log::error(
                "bcf_sr_set_targets",
                format_args!("Must call bcf_sr_set_targets() before bcf_sr_add_reader()"),
            );
            return Err(SyncedError::ApiUsage);
        }
        Ok(())
    }

    /// Sets the targets: records are streamed and kept (or with [`TargetsLogic::Exclude`]
    /// dropped) when they overlap a target (`bcf_sr_set_targets`). With `alleles_column`
    /// (0-based), records at a position are preferred when their type (indel or not) matches
    /// the alleles listed in that column of a tabix-indexed targets file.
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::ApiUsage`] when readers were already added or targets set.
    pub fn set_targets(
        &mut self,
        mut targets: Regions,
        logic: TargetsLogic,
        alleles_column: Option<usize>,
    ) -> Result<(), SyncedError> {
        self.check_targets_unset()?;
        targets.overlap = self.targets_overlap;
        self.targets = Some(targets);
        self.targets_logic = logic;
        self.targets_alleles = alleles_column;
        Ok(())
    }

    /// Reads targets from a comma-separated list or a file, either prefixed with `^` to
    /// exclude them, and sets them ([`Self::set_targets`]).
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::Regions`] when the list cannot be read, or the errors of
    /// [`Self::set_targets`].
    pub fn set_targets_spec(
        &mut self,
        spec: &[u8],
        source: ListSource,
        alleles_column: Option<usize>,
    ) -> Result<(), SyncedError> {
        self.check_targets_unset()?;
        let (logic, body) = match spec.strip_prefix(b"^") {
            Some(rest) => (TargetsLogic::Exclude, rest),
            None => (TargetsLogic::Include, spec),
        };
        let targets = Regions::new(body, source).map_err(|_| SyncedError::Regions)?;
        self.set_targets(targets, logic, alleles_column)
    }

    /// [`Self::add_reader`] for an input that was read once already, as bcftools reads some
    /// inputs twice: a stream is read on from where the first reading left it, without waiting
    /// for the writer of a FIFO ([`Reader::open_again`]).
    ///
    /// # Errors
    ///
    /// As [`Self::add_reader`].
    pub fn add_reader_again(&mut self, path: &str) -> Result<(), SyncedError> {
        let threads = self.next_reader_threads();
        let file = Reader::open_again(path, threads);
        self.add_opened_reader(path, file)
    }

    /// Opens an input and reads its header (`bcf_sr_add_reader`). `FILE##idx##INDEX` names
    /// the index explicitly. A missing BGZF EOF marker is warned about and kept for
    /// [`Self::take_pending_error`], as HTSlib leaves it in `errnum`.
    ///
    /// # Errors
    ///
    /// Returns HTSlib's errors for unreadable, unindexed, or malformed inputs, and for
    /// streaming several inputs without [`IndexRequirement::AllowMissing`].
    pub fn add_reader(&mut self, path: &str) -> Result<(), SyncedError> {
        let threads = self.next_reader_threads();
        let file = Reader::open(path, threads);
        self.add_opened_reader(path, file)
    }

    /// [`Self::add_reader`] with the input opened.
    fn add_opened_reader(
        &mut self,
        path: &str,
        file: Result<Reader, ReadError>,
    ) -> Result<(), SyncedError> {
        let (data_path, explicit_index) = split_idx_name(path);
        let open_failed = |e| {
            let err = SyncedError::Open(e);
            log::error(
                "hts_open_format",
                format_args!("Failed to open file \"{data_path}\" : {err}"),
            );
            err
        };
        let mut file = file.map_err(|e| match e {
            ReadError::Io(e) => open_failed(e),
            ReadError::UnknownBinary => {
                log::error("hts_hopen", format_args!("Failed to open file {data_path}"));
                open_failed(exec_format_error())
            }
            _ => SyncedError::FileType,
        })?;
        self.pending_error = None;
        if lacks_bgzf_eof(&file) {
            let p = file.path().unwrap_or("-");
            log::warn(
                "bcf_sr_add_hreader",
                format_args!("No BGZF EOF marker; file '{p}' may be truncated"),
            );
            self.pending_error = Some(SyncedError::NoEof);
        }
        let (access, header) =
            self.open_access(&mut file, data_path, explicit_index.map(PathBuf::from))?;
        let nreaders = self.readers.len() + 1;
        self.check_streaming(nreaders)?;
        let filter_ids = self.filter_ids(&header);
        // Without explicit regions, visit every sequence the readers know.
        if !self.explicit_regs && !self.streaming {
            let regions = self.regions.get_or_insert_with(Regions::empty);
            match &access {
                Access::Tabix(tbx) => {
                    for name in tbx.names() {
                        regions.add(name, -1, -1);
                    }
                }
                _ => {
                    for name in header.contig_names() {
                        regions.add(name, -1, -1);
                    }
                }
            }
            regions.sort_and_merge();
        }
        if self.require_index == IndexRequirement::AllowMissing && nreaders > 1 {
            check_same_sequences(&self.readers[0].header, &header)?;
        }
        self.readers.push(ReaderState {
            file,
            access,
            header,
            query: None,
            eof: false,
            fname: data_path.to_owned(),
            buffer: Vec::new(),
            nbuffer: 0,
            filter_ids,
            samples: Vec::new(),
            max_unpack: None,
        });
        self.has_line.push(false);
        Ok(())
    }

    /// Loads the index a new input needs (the explicit one, or the one found next to
    /// `data_path`) and reads its header; without a required index, streams it instead.
    fn open_access(
        &mut self,
        file: &mut Reader,
        data_path: &str,
        explicit_index: Option<PathBuf>,
    ) -> Result<(Access, Header), SyncedError> {
        if self.require_index != IndexRequirement::Required {
            let header = file.read_header().map_err(|_| SyncedError::Header)?;
            self.streaming = true;
            return Ok((Access::Stream, header));
        }
        if file.compression() != Compression::Bgzf {
            return Err(SyncedError::NotBgzf);
        }
        // HTSlib reports an index that cannot be read (`hts_idx_load3`).
        let load_failed = |idx_path: &Path, e: &io::Error| {
            report_load_error(idx_path, e);
            SyncedError::IndexLoad
        };
        let no_index = || {
            log::error(
                "idx_find_and_load",
                format_args!("Could not retrieve index file for '{data_path}'"),
            );
            SyncedError::IndexLoad
        };
        match file.format() {
            InputFormat::Vcf => {
                let idx_path = explicit_index
                    .or_else(|| locate(data_path, IndexFormat::Tbi))
                    .ok_or_else(no_index)?;
                warn_if_older(Path::new(data_path), &idx_path);
                let tbx = Tabix::load(&idx_path).map_err(|e| load_failed(&idx_path, &e))?;
                let header = file.read_header().map_err(|_| SyncedError::Header)?;
                Ok((Access::Tabix(Box::new(tbx)), header))
            }
            InputFormat::Bcf => {
                let header = file.read_header().map_err(|_| SyncedError::Header)?;
                let idx_path = explicit_index
                    .or_else(|| locate(data_path, IndexFormat::Csi))
                    .ok_or_else(no_index)?;
                warn_if_older(Path::new(data_path), &idx_path);
                let idx = Index::read(&idx_path).map_err(|e| load_failed(&idx_path, &e))?;
                Ok((Access::Bcf(Box::new(idx)), header))
            }
        }
    }

    /// Checks that streaming, once `nreaders` inputs are open, is allowed.
    fn check_streaming(&self, nreaders: usize) -> Result<(), SyncedError> {
        if self.streaming && nreaders > 1 {
            if self.require_index == IndexRequirement::AllowMissing {
                log::warn_once(
                    &NO_INDEX_WARNED,
                    "bcf_sr_add_hreader",
                    format_args!(
                        "Using multiple unindexed files may produce errors, make sure chromosomes are in the same order!"
                    ),
                );
            } else {
                log::error(
                    "bcf_sr_add_hreader",
                    format_args!(
                        "Must set require_index when the number of readers is greater than one"
                    ),
                );
                return Err(SyncedError::ApiUsage);
            }
        }
        if self.streaming && self.regions.is_some() {
            log::error(
                "bcf_sr_add_hreader",
                format_args!("Cannot tabix-jump in streaming mode"),
            );
            return Err(SyncedError::ApiUsage);
        }
        Ok(())
    }

    /// The header IDs of the FILTER names records must carry, `-1` for `.`. Names the header
    /// does not define are dropped; when none is left, HTSlib's `nfilter_ids` is 0 and no
    /// record is filtered.
    fn filter_ids(&self, header: &Header) -> Option<Vec<i32>> {
        self.apply_filters
            .as_ref()
            .map(|list| -> Vec<i32> {
                list.split(|&c| c == b',')
                    .filter_map(|f| {
                        if f == b"." {
                            Some(-1)
                        } else {
                            header.id(f).map(|i| i as i32)
                        }
                    })
                    .collect()
            })
            .filter(|ids| !ids.is_empty())
    }

    /// Closes reader `i` (`bcf_sr_remove_reader`); later readers move down by one.
    ///
    /// # Panics
    ///
    /// Panics when `i` is out of range.
    pub fn remove_reader(&mut self, i: usize) {
        assert!(
            i < self.readers.len(),
            "reader {i} of {}",
            self.readers.len()
        );
        self.sort.remove_reader(i);
        self.sort.reset();
        self.readers.remove(i);
        self.has_line.remove(i);
    }

    /// Restricts the samples to those every reader has (`bcf_sr_set_samples`): a list, or a
    /// file of names, or `-` for the first reader's samples; a `^` prefix excludes the listed
    /// samples instead. Samples missing from some reader are skipped with a warning. Returns
    /// `false` when no samples are left.
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::Samples`] when the list cannot be read.
    pub fn set_samples(&mut self, spec: &[u8], source: ListSource) -> Result<bool, SyncedError> {
        use std::collections::HashSet;

        let (exclude, body) = match spec.strip_prefix(b"^") {
            Some(rest) => (true, rest),
            None => (false, spec),
        };
        let list: Vec<Vec<u8>> = if exclude || spec != b"-" {
            readlist::read_list(body, source).map_err(|_| {
                log::error(
                    "bcf_sr_set_samples",
                    format_args!(
                        "Could not read the file: \"{}\"",
                        String::from_utf8_lossy(spec)
                    ),
                );
                SyncedError::Samples
            })?
        } else {
            Vec::new()
        };
        let first_samples = || -> Vec<Vec<u8>> {
            self.readers
                .first()
                .map(|r| r.header.samples().iter().map(|s| s.to_vec()).collect())
                .unwrap_or_default()
        };
        // An empty list (an empty file) leaves no samples, as in HTSlib.
        let (candidates, excluded): (Vec<Vec<u8>>, HashSet<Vec<u8>>) = if exclude {
            (first_samples(), list.into_iter().collect())
        } else if spec == b"-" {
            (first_samples(), HashSet::new())
        } else {
            (list, HashSet::new())
        };
        self.samples.clear();
        for name in candidates {
            if excluded.contains(&name) {
                continue;
            }
            match self
                .readers
                .iter()
                .find(|r| r.header.sample_id(&name).is_none())
            {
                Some(missing) => log::warn(
                    "bcf_sr_set_samples",
                    format_args!(
                        "The sample \"{}\" was not found in {}, skipping",
                        String::from_utf8_lossy(&name),
                        missing.fname
                    ),
                ),
                None => self.samples.push(name),
            }
        }
        if self.samples.is_empty() {
            if self.readers.len() > 1 {
                log::warn("bcf_sr_set_samples", format_args!("No samples in common"));
            }
            return Ok(false);
        }
        for r in &mut self.readers {
            r.samples = self
                .samples
                .iter()
                .map(|s| r.header.sample_id(s).map_or(-1, |i| i as i32))
                .collect();
        }
        Ok(true)
    }

    /// Returns the number of readers.
    #[inline]
    #[must_use]
    pub fn reader_count(&self) -> usize {
        self.readers.len()
    }

    /// Returns whether reader `i` has a record at the current position
    /// (`bcf_sr_has_line`).
    #[inline]
    #[must_use]
    pub fn has_line(&self, i: usize) -> bool {
        self.has_line[i]
    }

    /// Marks reader `i` as having no record at the current position, as bcftools does by
    /// clearing `has_line[i]` for filtered records.
    #[inline]
    pub fn clear_line(&mut self, i: usize) {
        self.has_line[i] = false;
    }

    /// Returns reader `i`'s record at the current position (`bcf_sr_get_line`).
    #[inline]
    #[must_use]
    pub fn line(&self, i: usize) -> Option<&Record> {
        self.has_line[i].then(|| &*self.readers[i].buffer[0])
    }

    /// Returns reader `i`'s record at the current position, for changes or lazy unpacking.
    #[inline]
    pub fn line_mut(&mut self, i: usize) -> Option<&mut Record> {
        self.has_line[i].then(|| &mut *self.readers[i].buffer[0])
    }

    /// Returns reader `i`'s record at the current position with its header.
    #[inline]
    pub fn line_with_header(&mut self, i: usize) -> Option<(&mut Record, &Header)> {
        if !self.has_line[i] {
            return None;
        }
        let r = &mut self.readers[i];
        Some((&mut *r.buffer[0], &r.header))
    }

    /// Returns reader `i`'s record at the current position with its header, for changes to
    /// both: [`Record::translate`] caches its translation tables in the source header.
    #[inline]
    pub fn line_with_header_mut(&mut self, i: usize) -> Option<(&mut Record, &mut Header)> {
        if !self.has_line[i] {
            return None;
        }
        let r = &mut self.readers[i];
        Some((&mut *r.buffer[0], &mut r.header))
    }

    /// Exchanges reader `i`'s current record with `rec` (`bcf_sr_swap_line`), to keep a
    /// record without copying it.
    pub fn swap_line(&mut self, i: usize, rec: &mut Record) {
        std::mem::swap(&mut *self.readers[i].buffer[0], rec);
    }

    /// [`Self::swap_line`] for a boxed record: swaps the boxes, not the records (a pointer
    /// swap, as HTSlib's reader hands out its buffer).
    pub fn swap_line_boxed(&mut self, i: usize, rec: &mut Box<Record>) {
        std::mem::swap(&mut self.readers[i].buffer[0], rec);
    }

    /// Returns the records reader `i` has read ahead (`buffer[1..=nbuffer]`): those at the
    /// current position not yet returned, then possibly one past it.
    pub fn buffered(&self, i: usize) -> impl ExactSizeIterator<Item = &Record> + '_ {
        let r = &self.readers[i];
        r.buffer[1..=r.nbuffer].iter().map(|b| &**b)
    }

    /// Returns the number of records reader `i` has read ahead (`reader->nbuffer`).
    #[inline]
    #[must_use]
    pub fn buffered_count(&self, i: usize) -> usize {
        self.readers[i].nbuffer
    }

    /// Returns slot `j` of reader `i`'s buffer (`reader->buffer[j]`): 0 holds the current
    /// record, valid while [`Self::has_line`], and `1..=`[`Self::buffered_count`] the records
    /// read ahead. bcftools merge works on these slots directly.
    ///
    /// # Panics
    ///
    /// Panics if `j` is past the records read ahead.
    #[inline]
    #[must_use]
    pub fn buffer_record(&self, i: usize, j: usize) -> &Record {
        let r = &self.readers[i];
        assert!(j <= r.nbuffer, "buffer slot {j} past {}", r.nbuffer);
        &r.buffer[j]
    }

    /// Returns slot `j` of reader `i`'s buffer for changes; see [`Self::buffer_record`]. A
    /// record swapped in must keep the slot's position (bcftools merge swaps in an empty
    /// record with the same `rid` and `pos`).
    ///
    /// # Panics
    ///
    /// Panics if `j` is past the records read ahead.
    #[inline]
    pub fn buffer_record_mut(&mut self, i: usize, j: usize) -> &mut Record {
        let r = &mut self.readers[i];
        assert!(j <= r.nbuffer, "buffer slot {j} past {}", r.nbuffer);
        &mut r.buffer[j]
    }

    /// Drops the first `n` records reader `i` has read ahead, once the caller has used them
    /// (bcftools merge's `clean_buffer`).
    ///
    /// # Panics
    ///
    /// Panics if `n` is more than [`Self::buffered_count`].
    pub fn discard_buffered(&mut self, i: usize, n: usize) {
        let r = &mut self.readers[i];
        assert!(n <= r.nbuffer, "discarding {n} of {} records", r.nbuffer);
        let nb = r.nbuffer;
        r.buffer[1..=nb].rotate_left(n);
        r.nbuffer -= n;
    }

    /// Returns whether reader `i` has nothing left in the current region
    /// (`bcf_sr_region_done`).
    #[must_use]
    pub fn region_done(&self, i: usize) -> bool {
        !self.has_line[i] && self.readers[i].nbuffer == 0
    }

    /// Returns reader `i`'s header (`bcf_sr_get_header`).
    #[inline]
    #[must_use]
    pub fn header(&self, i: usize) -> &Header {
        &self.readers[i].header
    }

    /// Returns reader `i`'s header for changes (adding lines before writing, say).
    #[inline]
    pub fn header_mut(&mut self, i: usize) -> &mut Header {
        &mut self.readers[i].header
    }

    /// Returns reader `i`'s file name, without an `##idx##` part (`reader->fname`).
    #[inline]
    #[must_use]
    pub fn file_name(&self, i: usize) -> &str {
        &self.readers[i].fname
    }

    /// Returns reader `i`'s underlying file, for fast paths that read raw lines or records
    /// while streaming one file without regions, targets, or filters.
    pub fn file_mut(&mut self, i: usize) -> &mut Reader {
        &mut self.readers[i].file
    }

    /// Returns the samples shared by all readers after [`Self::set_samples`].
    #[must_use]
    pub fn samples(&self) -> &[Vec<u8>] {
        &self.samples
    }

    /// Returns reader `i`'s header index of each of [`Self::samples`] (`reader->samples`).
    #[must_use]
    pub fn sample_map(&self, i: usize) -> &[i32] {
        &self.readers[i].samples
    }

    /// Returns the regions, set or derived from the indexes.
    #[must_use]
    pub fn regions(&self) -> Option<&Regions> {
        self.regions.as_ref()
    }

    /// Returns the targets, for the line and alleles of the current one.
    #[must_use]
    #[inline]
    pub fn targets(&self) -> Option<&Regions> {
        self.targets.as_ref()
    }

    /// Returns the targets, for their alleles or missed regions.
    pub fn targets_mut(&mut self) -> Option<&mut Regions> {
        self.targets.as_mut()
    }

    /// Returns whether the readers stream rather than jump through an index.
    #[must_use]
    pub fn is_streaming(&self) -> bool {
        self.streaming
    }

    /// Returns whether reader `i` is read straight through: streaming, without targets or
    /// FILTER requirements.
    #[must_use]
    pub fn is_plain_stream(&self, i: usize) -> bool {
        matches!(self.readers[i].access, Access::Stream)
            && self.targets.is_none()
            && self.readers[i].filter_ids.is_none()
    }

    /// Takes the error HTSlib leaves in `errnum` without stopping: a missing BGZF EOF
    /// marker of the last reader added, to be reported after reading everything.
    pub fn take_pending_error(&mut self) -> Option<SyncedError> {
        self.pending_error.take()
    }

    /// Returns the error [`Self::take_pending_error`] would take, leaving it in place (HTSlib's
    /// `errnum` as tools check it after each line).
    #[must_use]
    pub fn pending_error(&self) -> Option<&SyncedError> {
        self.pending_error.as_ref()
    }

    /// Restarts at the first region (`bcf_sr_seek(readers, NULL, 0)`).
    pub fn seek_start(&mut self) {
        self.sort.reset();
        if let Some(r) = &mut self.regions {
            r.rewind();
        }
    }

    /// Positions all readers at `seq:pos` (0-based) and the regions accordingly
    /// (`bcf_sr_seek`). Returns the number of readers with records on the sequence. Without
    /// regions (streaming), nothing happens.
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::Seek`] when an index cannot be queried or a line of a
    /// tabix-indexed regions file cannot be read.
    pub fn seek(&mut self, seq: &[u8], pos: i64) -> Result<usize, SyncedError> {
        let Some(regions) = &mut self.regions else {
            return Ok(0);
        };
        self.sort.reset();
        regions.rewind();
        regions.seek_to(seq, pos).map_err(|_| SyncedError::Seek)?;
        let mut n = 0;
        for i in 0..self.readers.len() {
            if self.reader_seek(i, seq, pos, MAX_CSI_COOR - 1)? {
                n += 1;
            }
        }
        Ok(n)
    }

    /// Resets reader `i` to iterate over `seq:start-end` (HTSlib's `_reader_seek`); `false`
    /// when the sequence is not in its index.
    fn reader_seek(
        &mut self,
        i: usize,
        seq: &[u8],
        start: i64,
        end: i64,
    ) -> Result<bool, SyncedError> {
        if end >= MAX_CSI_COOR {
            log::error(
                "_reader_seek",
                format_args!("The coordinate is out of csi index limit: {}", end + 1),
            );
            return Err(SyncedError::Seek);
        }
        let reader = &mut self.readers[i];
        reader.query = None;
        reader.nbuffer = 0;
        let plan = match &reader.access {
            Access::Tabix(tbx) => tbx
                .name_to_id(seq)
                .map(|tid| tbx.index.query(tid as usize, start, end + 1)),
            Access::Bcf(idx) => reader
                .header
                .contig_id(seq)
                .map(|tid| idx.query(tid as usize, start, end + 1)),
            Access::Stream => None,
        };
        match plan {
            None => Ok(false),
            Some(None) => {
                log::error(
                    "_reader_seek",
                    format_args!(
                        "Could not seek: {}:{}-{}",
                        String::from_utf8_lossy(seq),
                        start + 1,
                        end + 1
                    ),
                );
                Err(SyncedError::Seek)
            }
            Some(Some(plan)) => {
                reader.query = Some(Query::new(plan));
                Ok(true)
            }
        }
    }

    /// Moves to the next region when every reader is done with the current one (HTSlib's
    /// `_readers_next_region`); `false` when no regions are left.
    fn readers_next_region(&mut self) -> Result<bool, SyncedError> {
        if !self.readers.iter().all(ReaderState::is_done) {
            return Ok(true);
        }
        let regions = self.regions.as_mut().expect("called with regions");
        let prev_iseq = regions.iseq;
        let prev_end = regions.end;
        // A line of a tabix-indexed regions file that cannot be read ends the iteration with
        // HTSlib's `bcf_sr_seek_error`.
        if !regions.next_region().map_err(|_| SyncedError::Seek)? {
            return Ok(false);
        }
        regions.prev_end = if prev_iseq == regions.iseq {
            prev_end
        } else {
            -1
        };
        let (iseq, start, end) = regions.current().expect("a region is current");
        let seq = regions.seq_names()[iseq].clone();
        for i in 0..self.readers.len() {
            self.reader_seek(i, &seq, start, end)?;
        }
        Ok(true)
    }

    /// Reads one record into the reader's next free slot; `false` at the end of the region
    /// or stream.
    fn read_one(
        reader: &mut ReaderState,
        line: &mut Vec<u8>,
        slot: usize,
    ) -> Result<bool, SyncedError> {
        let rec = &mut *reader.buffer[slot];
        match &mut reader.access {
            Access::Stream => {
                if reader.eof {
                    return Ok(false);
                }
                match reader.file.read_record(&mut reader.header, line, rec) {
                    Ok(got) => {
                        reader.eof = !got;
                        Ok(got)
                    }
                    Err(ReadError::Parse(_)) => Err(SyncedError::VcfParse),
                    Err(_) => Err(SyncedError::BcfRead),
                }
            }
            Access::Tabix(tbx) => {
                let Some(q) = reader.query.as_mut() else {
                    return Ok(false);
                };
                let got = reader
                    .file
                    .query_next_line(q, tbx, line)
                    .map_err(|_| SyncedError::BcfRead)?;
                if got {
                    crate::variant::vcf::parse(&mut reader.header, line, rec)
                        .map_err(|_| SyncedError::VcfParse)?;
                }
                Ok(got)
            }
            Access::Bcf(_) => {
                let Some(q) = reader.query.as_mut() else {
                    return Ok(false);
                };
                let got = reader
                    .file
                    .query_next_bcf(q, &reader.header, rec)
                    .map_err(|_| SyncedError::BcfRead)?;
                if got {
                    rec.subset_format(&reader.header);
                }
                Ok(got)
            }
        }
    }

    /// Reads records until one at a different position than the first buffered one (HTSlib's
    /// `_reader_fill_buffer`).
    fn fill_buffer(&mut self, i: usize) -> Result<(), SyncedError> {
        let reader = &mut self.readers[i];
        if reader.nbuffer > 0 && reader.buffer[reader.nbuffer].pos() != reader.buffer[1].pos() {
            return Ok(());
        }
        if reader.query.is_none() && !self.streaming {
            return Ok(());
        }
        let mut at_end = false;
        loop {
            let slot = reader.nbuffer + 1;
            while reader.buffer.len() <= slot {
                let mut rec = Box::new(Record::new());
                rec.set_max_unpack(reader.max_unpack.unwrap_or(self.max_unpack));
                rec.set_pos(-1);
                reader.buffer.push(rec);
            }
            if !Self::read_one(reader, &mut self.line, slot)? {
                at_end = true;
                break;
            }
            // Skip records seen in the previous region, and recognize true variant overlaps
            // as opposed to record overlaps (TA>T versus A>-).
            if let Some(regions) = &self.regions {
                let (beg, end) = self.regions_overlap.interval(&mut reader.buffer[slot]);
                if beg <= regions.prev_end || end < regions.start || beg > regions.end {
                    continue;
                }
            }
            if let Some(ids) = &reader.filter_ids
                && !reader.buffer[slot].has_filter_in(ids)
            {
                continue;
            }
            reader.nbuffer += 1;
            let (last, first) = (&reader.buffer[reader.nbuffer], &reader.buffer[1]);
            if last.rid() != first.rid() || last.pos() != first.pos() {
                break;
            }
        }
        if at_end {
            reader.query = None;
        }
        if self.require_index == IndexRequirement::AllowMissing
            && reader.buffer.len() > 1
            && reader.buffer[reader.nbuffer].rid() < reader.buffer[1].rid()
        {
            log::error(
                "_reader_fill_buffer",
                format_args!(
                    "Sequences out of order, cannot stream multiple unindexed files: {}",
                    reader.fname
                ),
            );
            return Err(SyncedError::NoIndex);
        }
        Ok(())
    }

    /// Advances to the next position (or the next records at the same position) and makes
    /// each reader's record there current (`bcf_sr_next_line`). Returns the number of readers
    /// with a record ([`Self::has_line`]), or 0 when all are done.
    ///
    /// # Errors
    ///
    /// Returns read, parse, and seek errors, as HTSlib sets `errnum`.
    pub fn next_line(&mut self) -> Result<usize, SyncedError> {
        let Some(col) = self.targets_alleles else {
            return self.next_line_inner();
        };
        loop {
            let ret = self.next_line_inner()?;
            if ret == 0 {
                return Ok(0);
            }
            let Some(i) = self.has_line.iter().position(|&h| h) else {
                return Ok(ret);
            };
            let targets = self.targets.as_mut().expect("alleles need targets");
            if targets.match_alleles(col, &mut self.readers[i].buffer[0])? {
                return Ok(ret);
            }
            // Return this line as if it matched when no duplicate lines are left.
            let more = self
                .readers
                .iter()
                .zip(&self.has_line)
                .any(|(r, &h)| h && r.nbuffer > 0 && r.buffer[1].pos() == r.buffer[0].pos());
            if !more {
                return Ok(ret);
            }
        }
    }

    /// HTSlib's `next_line`.
    fn next_line_inner(&mut self) -> Result<usize, SyncedError> {
        let mut have_chr = false;
        let mut min_pos = i64::MAX;
        loop {
            if self.regions.is_some() && !self.readers_next_region()? {
                break;
            }
            let allow_no_idx = self.require_index == IndexRequirement::AllowMissing;
            let mut min_rid = i32::MAX;
            for i in 0..self.readers.len() {
                self.fill_buffer(i)?;
                let r = &self.readers[i];
                if allow_no_idx && r.nbuffer > 0 {
                    min_rid = min_rid.min(r.buffer[1].rid());
                }
            }
            for i in 0..self.readers.len() {
                let r = &self.readers[i];
                if r.nbuffer == 0 || (allow_no_idx && min_rid != r.buffer[1].rid()) {
                    continue;
                }
                let pos = r.buffer[1].pos();
                if min_pos > pos {
                    min_pos = pos;
                    let Some(name) = r.header.contig_name(r.buffer[1].rid()) else {
                        return Err(SyncedError::BcfRead);
                    };
                    self.chr.clear();
                    self.chr.extend_from_slice(name);
                    have_chr = true;
                    self.sort.set_active(i);
                } else if min_pos == pos {
                    self.sort.add_active(i);
                }
            }
            if min_pos == i64::MAX {
                if self.regions.is_none() {
                    break;
                }
                continue;
            }
            if !self.keep_targets(min_pos)? {
                min_pos = i64::MAX;
                have_chr = false;
                continue;
            }
            break;
        }
        if !have_chr {
            return Ok(0);
        }
        Ok(self
            .sort
            .next(&mut self.readers, &mut self.has_line, &self.chr, min_pos))
    }

    /// Drops the records at `min_pos` outside the targets; `false` when none is left.
    ///
    /// # Errors
    ///
    /// Returns [`SyncedError::TargetsRead`] when a line of a tabix-indexed targets file cannot
    /// be read.
    fn keep_targets(&mut self, min_pos: i64) -> Result<bool, SyncedError> {
        let Some(targets) = &mut self.targets else {
            return Ok(true);
        };
        let mut matched = false;
        for i in 0..self.readers.len() {
            let r = &mut self.readers[i];
            if r.nbuffer == 0 || r.buffer[1].pos() != min_pos {
                continue;
            }
            let (beg, end) = match self.targets_overlap {
                Overlap::Pos => (min_pos, min_pos),
                Overlap::Record => (min_pos, min_pos + r.buffer[1].rlen() - 1),
                Overlap::Variant => self.targets_overlap.interval(&mut r.buffer[1]),
            };
            // HTSlib reports a malformed targets line and counts the record as outside the
            // targets; the next record goes on from the next line. It does the same for a line
            // that cannot be read, which loses the records after it: that is an error here.
            let overlap = match targets.overlaps(&self.chr, beg, end) {
                Ok(overlap) => overlap,
                Err(e) if e.is_read_failure() => return Err(SyncedError::TargetsRead(e)),
                Err(_) => false,
            };
            if overlap == (self.targets_logic == TargetsLogic::Exclude) {
                r.shift_buffer();
            } else {
                matched = true;
            }
        }
        Ok(matched)
    }
}

/// Options for [`SyncedReader::open`].
#[derive(Debug, Default)]
pub struct SyncedOptions {
    /// Regions to visit through the index (`-r`/`-R`); `None` streams the whole file. Their
    /// [`Regions::overlap`] sets the overlap rule.
    pub regions: Option<Regions>,
    /// Targets to keep while streaming (`-t`/`-T`), with their overlap rule.
    pub targets: Option<Regions>,
    /// Whether the targets are excluded rather than kept (`^` prefix).
    pub targets_exclude: bool,
    /// FILTER names to require, `.` for none (`-f`).
    pub apply_filters: Option<Vec<u8>>,
    /// BGZF decompression threads.
    pub threads: usize,
}

/// A single VCF/BCF input read in order, restricted to regions and targets: a
/// [`SyncedReaders`] with one reader, as bcftools' single-file commands use it.
pub struct SyncedReader {
    inner: SyncedReaders,
}

impl std::fmt::Debug for SyncedReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SyncedReader").finish_non_exhaustive()
    }
}

impl SyncedReader {
    /// Opens an input and reads its header (HTSlib's `bcf_sr_add_reader`). With regions, the
    /// input must be BGZF-compressed and indexed (`FILE##idx##INDEX` names the index).
    ///
    /// # Errors
    ///
    /// Returns HTSlib's errors for unreadable, unindexed, or malformed inputs.
    pub fn open(path: &str, options: SyncedOptions) -> Result<Self, SyncedError> {
        let mut inner = SyncedReaders::new();
        inner.set_threads(options.threads);
        if let Some(t) = options.targets {
            inner.set_targets_overlap(t.overlap);
            let logic = if options.targets_exclude {
                TargetsLogic::Exclude
            } else {
                TargetsLogic::Include
            };
            inner.set_targets(t, logic, None)?;
        }
        if let Some(r) = options.regions {
            inner.set_regions_overlap(r.overlap);
            inner.set_regions(r);
        }
        inner.set_apply_filters(options.apply_filters.as_deref());
        inner.add_reader(path)?;
        Ok(Self { inner })
    }

    /// Returns the header.
    #[inline]
    #[must_use]
    pub fn header(&self) -> &Header {
        self.inner.header(0)
    }

    /// Returns the header for changes (adding lines before writing, for example).
    #[inline]
    pub fn header_mut(&mut self) -> &mut Header {
        self.inner.header_mut(0)
    }

    /// Returns the underlying reader, for fast paths that read raw lines or records while
    /// streaming without regions or targets.
    pub fn reader_mut(&mut self) -> &mut Reader {
        self.inner.file_mut(0)
    }

    /// The error for a VCF line that does not parse, read by a fast path from
    /// [`Self::reader_mut`]: [`SyncedError::VcfParse`], or [`SyncedError::BcfRead`] when the
    /// input is plain gzip whose rest does not decompress ([`Reader::check_gzip_rest`]), as the
    /// line is then likely garbage from a corrupt deflate stream.
    pub fn vcf_parse_error(&mut self) -> SyncedError {
        match self.reader_mut().check_gzip_rest() {
            Ok(()) => SyncedError::VcfParse,
            Err(_) => SyncedError::BcfRead,
        }
    }

    /// Returns whether records are read straight through (no regions, targets, or filters).
    #[must_use]
    pub fn is_plain_stream(&self) -> bool {
        self.inner.is_plain_stream(0)
    }

    /// Takes the error to report after reading everything (a missing BGZF EOF marker, which
    /// HTSlib reports but reads past).
    pub fn take_pending_error(&mut self) -> Option<SyncedError> {
        self.inner.take_pending_error()
    }

    /// Returns the error [`Self::take_pending_error`] would take, leaving it in place.
    #[must_use]
    pub fn pending_error(&self) -> Option<&SyncedError> {
        self.inner.pending_error()
    }

    /// Returns the sequences that HTSlib visits when several inputs have no explicit regions:
    /// the index's sequences in index order for indexed VCF (`tbx_seqnames`), otherwise the
    /// header's contigs (`bcf_hdr_seqnames`).
    #[must_use]
    pub fn seq_names(&self) -> Vec<&[u8]> {
        let reader = &self.inner.readers[0];
        match &reader.access {
            Access::Tabix(tbx) => tbx.names().iter().map(Vec::as_slice).collect(),
            Access::Bcf(_) | Access::Stream => reader.header.contig_names(),
        }
    }

    /// Replaces the regions of an input opened with regions, and starts again from the first
    /// one.
    ///
    /// # Panics
    ///
    /// Panics when the input was opened without regions (it has no index to jump with).
    pub fn set_regions(&mut self, regions: Regions) {
        assert!(
            !self.inner.is_streaming(),
            "set_regions needs an input opened with regions"
        );
        self.inner.set_regions(regions);
    }

    /// Limits how much of each VCF line is parsed (HTSlib's `bcf_srs_t.max_unpack`, see
    /// [`Record::set_max_unpack`]); empty, the default, parses everything.
    pub fn set_max_unpack(&mut self, which: Unpack) {
        self.inner.set_max_unpack(which);
    }

    /// Returns the underlying reader set.
    pub fn readers_mut(&mut self) -> &mut SyncedReaders {
        &mut self.inner
    }

    /// Reads the next record that passes the regions, targets, and filters (HTSlib's
    /// `bcf_sr_next_line` with one reader) into `rec`. Returns `false` at the end.
    ///
    /// # Errors
    ///
    /// Returns an error for unparsable or unreadable records.
    pub fn next_record(&mut self, rec: &mut Record) -> Result<bool, SyncedError> {
        if self.inner.next_line()? == 0 {
            return Ok(false);
        }
        self.inner.swap_line(0, rec);
        Ok(true)
    }

    /// [`Self::next_record`] into a boxed record, which swaps the boxes instead of the
    /// records: cheaper for every record read.
    ///
    /// # Errors
    ///
    /// Returns an error for unparsable or unreadable records.
    pub fn next_record_boxed(&mut self, rec: &mut Box<Record>) -> Result<bool, SyncedError> {
        if self.inner.next_line()? == 0 {
            return Ok(false);
        }
        self.inner.swap_line_boxed(0, rec);
        Ok(true)
    }
}

impl From<SyncedReaders> for SyncedReader {
    /// Wraps a reader set with one reader, set up by hand (an index requirement, say).
    fn from(inner: SyncedReaders) -> Self {
        Self { inner }
    }
}

// The reader set moves between threads with its owner.
const _: () = {
    const fn assert_send<T: Send>() {}
    assert_send::<SyncedReaders>();
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_codes_match_htslib() {
        assert_eq!(SyncedError::NotBgzf.code(), 1);
        assert_eq!(SyncedError::NoIndex.code(), 10);
        assert_eq!(SyncedError::Samples.code(), 14);
        let e = SyncedError::Open(io::Error::from_raw_os_error(2));
        assert_eq!(e.to_string(), "No such file or directory");
    }

    #[test]
    fn identity_survives_rotation() {
        let mut v: Vec<Box<Record>> = (0..3).map(|_| Box::new(Record::new())).collect();
        let id = rec_id(&v[2]);
        v.rotate_left(1);
        assert_eq!(rec_id(&v[1]), id);
    }

    /// Readers read in step share the threads instead of starting `threads` each.
    #[test]
    fn readers_share_the_threads() {
        let shares = |threads, readers, opened| {
            let mut sr = SyncedReaders::new();
            sr.set_shared_threads(threads, readers);
            (0..opened)
                .map(|_| sr.next_reader_threads())
                .collect::<Vec<_>>()
        };
        let many = shares(8, 2000, 2000);
        assert_eq!(many.iter().sum::<usize>(), 8);
        assert!(many.iter().all(|&t| t <= 1));
        assert_eq!(shares(8, 2, 2), [4, 4]);
        assert_eq!(shares(1, 2, 2), [1, 0]);
        assert_eq!(shares(5, 2, 4), [3, 2, 3, 2]);
        assert_eq!(shares(0, 3, 3), [0, 0, 0]);
        let mut sr = SyncedReaders::new();
        sr.set_threads(4);
        assert_eq!([sr.next_reader_threads(), sr.next_reader_threads()], [4, 4]);
    }
}
