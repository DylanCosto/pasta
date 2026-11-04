//! Choosing a compression method for each block, as HTSlib 1.24 does (`cram_compress_slice` and
//! `cram_compress_block3`).
//!
//! Every block content ID keeps metrics. For the first few blocks, and again every so often, each
//! candidate method compresses the block and the smallest result is kept; the method that does
//! best on average (weighted slightly in favor of faster methods) is then used alone until the
//! next trial. Methods that keep losing by a wide margin stop being tried.

use std::{collections::HashMap, io};

use crate::cram::{
    codecs::{aac, bzip2, fqzcomp, gzip, lzma, name_tokenizer, rans_4x8, rans_nx16},
    container::{
        block::{self, CompressionMethod, ContentType},
        compression_header::data_series_encodings::DataSeries,
    },
    file_definition::Version,
    io::writer::container::Block,
};

/// HTSlib's compression profiles (`fast`, `normal`, `small`, `archive`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Profile {
    /// Level 1, no name tokenizer, 10,000 records per slice.
    Fast,
    /// Level 5, rANS and (for CRAM 3.1) the name tokenizer.
    #[default]
    Normal,
    /// Level 6, adding bzip2 and fqzcomp, 25,000 records per slice.
    Small,
    /// Level 7, adding bzip2, lzma, fqzcomp, and arithmetic coding, 100,000 records per slice.
    Archive,
}

/// Which compression methods the writer may use, and how hard it tries, with the slice sizes and
/// HTSlib's other CRAM output options (`store_md`, `store_nm`, `ignore_md5`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Options {
    /// The CRAM version written. Version 3.0 cannot use rANS Nx16, the arithmetic coder, the
    /// name tokenizer, or fqzcomp.
    pub version: Version,
    /// The compression level (1–9; 0 stores blocks uncompressed).
    pub level: u8,
    /// Records per slice.
    pub records_per_slice: usize,
    /// Bases per slice: a slice also ends once its records hold this many bases (HTSlib's
    /// default is 500 per record of `records_per_slice`).
    pub bases_per_slice: usize,
    /// bzip2.
    pub use_bzip2: bool,
    /// lzma (xz).
    pub use_lzma: bool,
    /// rANS (4x8 for CRAM 3.0, Nx16 for 3.1).
    pub use_rans: bool,
    /// The name tokenizer for read names (CRAM 3.1).
    pub use_tok: bool,
    /// fqzcomp for quality scores (CRAM 3.1).
    pub use_fqz: bool,
    /// The adaptive arithmetic coder (CRAM 3.1).
    pub use_arith: bool,
    /// Keeps the `MD` tags that readers regenerate, which are otherwise left out (HTSlib's
    /// `store_md`).
    pub store_md: bool,
    /// Keeps the `NM` tags that readers regenerate, which are otherwise left out (HTSlib's
    /// `store_nm`).
    pub store_nm: bool,
    /// Stores records against reference sequences whose `@SQ` `M5` field does not match them,
    /// rather than failing (HTSlib's `ignore_md5`). The header keeps the mismatched checksums.
    pub ignore_md5: bool,
}

/// HTSlib's default compression level (`CRAM_DEFAULT_LEVEL`).
const DEFAULT_LEVEL: u8 = 5;

/// HTSlib's default records per slice (`SEQS_PER_SLICE`).
const DEFAULT_RECORDS_PER_SLICE: usize = 10_000;

/// HTSlib's default bases per slice for each record of a slice (`BASES_PER_SLICE`).
pub const BASES_PER_RECORD: usize = 500;

impl Options {
    /// The options of a profile for a CRAM version, as HTSlib sets them.
    pub fn new(version: Version, profile: Profile) -> Self {
        let is_3_1 = is_3_1_or_later(version);

        let mut options = Self {
            version,
            level: DEFAULT_LEVEL,
            records_per_slice: DEFAULT_RECORDS_PER_SLICE,
            bases_per_slice: DEFAULT_RECORDS_PER_SLICE * BASES_PER_RECORD,
            use_bzip2: false,
            use_lzma: false,
            use_rans: version.major() >= 3,
            use_tok: is_3_1,
            use_fqz: false,
            use_arith: false,
            store_md: false,
            store_nm: false,
            ignore_md5: false,
        };

        match profile {
            Profile::Fast => {
                options.level = 1;
                options.use_tok = false;
            }
            Profile::Normal => {}
            Profile::Small => {
                options.level = 6;
                options.use_bzip2 = true;
                options.use_fqz = true;
                options.records_per_slice = 25_000;
            }
            Profile::Archive => {
                options.level = 7;
                options.use_bzip2 = true;
                options.use_fqz = true;
                options.use_arith = true;
                options.use_lzma = true;
                options.records_per_slice = 100_000;
            }
        }

        options.bases_per_slice = options.records_per_slice * BASES_PER_RECORD;
        options
    }

    /// Whether CRAM 3.1 codecs may be used.
    fn is_3_1(&self) -> bool {
        is_3_1_or_later(self.version)
    }
}

fn is_3_1_or_later(version: Version) -> bool {
    (version.major(), version.minor()) >= (3, 1)
}

impl Default for Options {
    /// CRAM 3.1 with the normal profile, samtools' default output.
    fn default() -> Self {
        Self::new(Version::new(3, 1), Profile::Normal)
    }
}

// HTSlib's internal method numbers (`enum cram_block_method_int`). Bit `m` of a method mask
// selects method `m`.
const RAW: usize = 0;
const GZIP: usize = 1;
const BZIP2: usize = 2;
const LZMA: usize = 3;
const RANS0: usize = 4;
const RANS_PR0: usize = 5;
const ARITH_PR0: usize = 6;
const FQZ: usize = 7;
const TOK3: usize = 8;
const GZIP_1: usize = 12;
const RANS1: usize = 16;
const RANS_PR1: usize = 17;
const RANS_PR64: usize = 18;
const RANS_PR9: usize = 19;
const RANS_PR128: usize = 20;
const RANS_PR129: usize = 21;
const RANS_PR192: usize = 22;
const RANS_PR193: usize = 23;
const ARITH_PR1: usize = 25;
const ARITH_PR193: usize = 31;

const METHOD_COUNT: usize = 32;

/// The rANS Nx16 (and arithmetic coder) flags of `RANS_PR*` and `ARITH_PR*` methods.
const PR_FLAGS: [u8; 8] = [0, 1, 64, 9, 128, 129, 192, 193];

/// Relative costs (speed) of the methods, by number (`meth_cost`).
const METHOD_COSTS: [f64; METHOD_COUNT] = [
    1.00, 1.04, 1.07, 1.08, 1.00, 1.00, 1.04, 1.05, 1.05, 1.00, 1.00, 1.01, 1.01, 1.05, 1.05, 1.05,
    1.01, 1.01, 1.00, 1.03, 1.00, 1.01, 1.00, 1.01, 1.07, 1.04, 1.04, 1.04, 1.03, 1.04, 1.04, 1.04,
];

/// Blocks between trials.
const TRIAL_SPAN: i32 = 70;
/// Blocks per trial.
const TRIALS: i32 = 3;
/// A method is dropped after this many trials in a row losing by this much in total.
const MAX_FAILS: i32 = 4;
const MAX_DELTA: f64 = 0.20;

/// rANS Nx16 interleaves 32 states for blocks at least this large (HTSlib's `SIMD_AUTO`).
const RANS_32_STATE_MIN_LEN: usize = 50_000;

const fn bit(method: usize) -> u32 {
    1 << method
}

/// The candidate methods for the kinds of blocks, from the options (`cram_compress_slice`).
#[derive(Clone, Copy, Debug)]
struct MethodSets {
    /// Most data series and tag values.
    method: u32,
    /// Data series that need only entropy coding (the same as `method` from level 5).
    method_fast: u32,
    /// Quality scores.
    quality_scores: u32,
    /// Read names.
    names: u32,
}

impl MethodSets {
    fn new(options: &Options) -> Self {
        let level = options.level;
        let is_3_1 = options.is_3_1();

        // `GZIP_RLE` is left out: libdeflate builds of HTSlib drop it in favor of `GZIP_1`.
        let mut method = bit(GZIP);

        if options.use_bzip2 {
            method |= bit(BZIP2);
        }

        let method_rans = bit(RANS0) | bit(RANS1);
        let mut method_ranspr = method_rans;

        if options.use_rans {
            method_ranspr = bit(RANS_PR0) | bit(RANS_PR1);

            if level > 1 {
                method_ranspr |= bit(RANS_PR64) | bit(RANS_PR9) | bit(RANS_PR128) | bit(RANS_PR193);
            }

            if level > 5 {
                method_ranspr |= bit(RANS_PR129) | bit(RANS_PR192);
            }

            method |= if is_3_1 { method_ranspr } else { method_rans };
        }

        if options.use_arith && is_3_1 {
            method |= bit(ARITH_PR0) | bit(ARITH_PR1);

            if level > 1 {
                method |= (ARITH_PR1 + 1..=ARITH_PR193).fold(0, |m, i| m | bit(i));
            }
        }

        if options.use_lzma {
            method |= bit(LZMA);
        }

        let mut method_fast = method & !(bit(GZIP) | bit(BZIP2) | bit(LZMA));

        if level >= 5 {
            method |= bit(GZIP_1);
            method_fast = method;
        }

        if level == 1 {
            method &= !bit(GZIP);
            method |= bit(GZIP_1);
            method_fast = method;
        }

        let mut quality_scores = method;

        // HTSlib also tries other fqzcomp parameter sets, which fqzcomp here does not have.
        if is_3_1 && options.use_fqz {
            quality_scores |= bit(FQZ);
        }

        // Read names: best is generally xz, bzip2, zlib, then the name tokenizer.
        let mut names = method & !(method_rans | method_ranspr);

        if is_3_1 && options.use_tok {
            names |= bit(TOK3);
        }

        Self {
            method,
            method_fast,
            quality_scores,
            names,
        }
    }
}

/// What a block holds, as far as choosing its compression goes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Kind {
    QualityScores,
    Bases,
    Names,
    /// Inserted bases, stretches of bases, soft clips, and mate reference IDs.
    Sequence,
    TagValues,
    Other,
}

impl Kind {
    fn of(content_id: block::ContentId) -> Self {
        let id = |series: DataSeries| block::ContentId::from(series);

        match content_id {
            c if c == id(DataSeries::QualityScores) => Self::QualityScores,
            c if c == id(DataSeries::Bases) => Self::Bases,
            c if c == id(DataSeries::Names) => Self::Names,
            c if c == id(DataSeries::InsertionBases)
                || c == id(DataSeries::StretchesOfBases)
                || c == id(DataSeries::MateReferenceSequenceIds) =>
            {
                Self::Sequence
            }
            // Data series use small content IDs; tags use larger ones (from their key).
            c if c > 255 => Self::TagValues,
            _ => Self::Other,
        }
    }
}

/// Learned compression results for one block content ID (`cram_metrics`).
#[derive(Clone, Debug)]
struct Metrics {
    trial: i32,
    next_trial: i32,
    /// Accumulated compressed sizes by method.
    sizes: [f64; METHOD_COUNT],
    /// The method used between trials.
    method: usize,
    /// The candidate methods, narrowed as methods keep losing; 0 until the first trial.
    revised_methods: u32,
    consistency: i32,
    fails: [i32; METHOD_COUNT],
    extra: [f64; METHOD_COUNT],
    input_avg_size: f64,
    input_avg_delta: f64,
}

impl Default for Metrics {
    fn default() -> Self {
        Self {
            // Learn quicker at the start.
            trial: TRIALS - 1,
            next_trial: TRIAL_SPAN / 2,
            sizes: [0.0; METHOD_COUNT],
            method: RAW,
            revised_methods: 0,
            consistency: 0,
            fails: [0; METHOD_COUNT],
            extra: [0.0; METHOD_COUNT],
            input_avg_size: 0.0,
            input_avg_delta: 0.0,
        }
    }
}

/// The record data that fqzcomp needs: each record's quality score count.
pub(crate) struct QualityLayout<'a> {
    /// The read lengths, when every record stores its quality scores as an array.
    pub read_lengths: Option<&'a [usize]>,
}

impl Options {
    /// Compresses the core data block: with gzip at level 1 only at high levels.
    pub(crate) fn compress_core(
        &self,
        content_id: block::ContentId,
        src: &[u8],
    ) -> io::Result<Block> {
        if self.level > 5 && src.len() > 500 {
            let data = gzip::encode(flate2::Compression::new(1), src)?;

            if data.len() < src.len() {
                return Ok(block(
                    ContentType::CoreData,
                    content_id,
                    CompressionMethod::Gzip,
                    src,
                    data,
                ));
            }
        }

        Ok(raw_block(ContentType::CoreData, content_id, src))
    }
}

/// Learns the best method for each block content ID.
///
/// Methods are chosen block by block in the order the blocks are written ([`Self::choose`]), since
/// each choice depends on the blocks before it. Compressing with a chosen [`Method`] can then
/// happen on any thread; a [`Trial`] must report what its methods gave
/// ([`Self::end_trial`]) before the next block with its content ID is chosen for.
#[derive(Debug)]
pub(crate) struct Compressor {
    options: Options,
    methods: MethodSets,
    metrics: HashMap<block::ContentId, Metrics>,
}

impl Compressor {
    pub(crate) fn new(options: Options) -> Self {
        Self {
            methods: MethodSets::new(&options),
            options,
            metrics: HashMap::new(),
        }
    }

    /// Chooses how to compress the next external block with `content_id`, of `len` bytes.
    pub(crate) fn choose(&mut self, content_id: block::ContentId, len: usize) -> Choice {
        let options = self.options;
        let level = options.level;

        let method = |method, level| {
            Choice::Method(Method {
                options,
                content_id,
                method,
                level,
            })
        };

        if level == 0 || len == 0 {
            return method(RAW, level);
        }

        let kind = Kind::of(content_id);

        let methods = match kind {
            Kind::QualityScores => self.methods.quality_scores,
            Kind::Names => self.methods.names,
            Kind::Bases | Kind::Sequence | Kind::TagValues => self.methods.method,
            Kind::Other => self.methods.method_fast,
        };

        // At levels 1 and 2, HTSlib compresses quality scores and bases with level 1.
        let level = match kind {
            Kind::QualityScores if level < 3 => 1,
            Kind::Bases if level < 3 => 1,
            _ => level,
        };

        let metrics = self.metrics.entry(content_id).or_default();

        match begin(metrics, methods, len) {
            Step::Use(m) => method(m, level),
            Step::Trial(methods) => Choice::Trial(Trial {
                options,
                content_id,
                methods,
                level,
            }),
        }
    }

    /// Learns from a trial of the last block chosen for with the trial's content ID.
    pub(crate) fn end_trial(&mut self, trial: &Trial, sizes: TrialSizes) {
        let metrics = self
            .metrics
            .get_mut(&trial.content_id)
            .expect("a trial follows a choice for its content ID");

        end(metrics, trial.methods, trial.level, sizes.0);
    }
}

/// How to compress a block (see [`Compressor::choose`]).
#[derive(Clone, Copy, Debug)]
pub(crate) enum Choice {
    /// Compress the block with this method.
    Method(Method),
    /// Compress the block with each of the trial's methods, keep the smallest result, and report
    /// them all with [`Compressor::end_trial`].
    Trial(Trial),
}

/// A compression method for an external block.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Method {
    options: Options,
    content_id: block::ContentId,
    method: usize,
    level: u8,
}

impl Method {
    pub(crate) fn content_id(&self) -> block::ContentId {
        self.content_id
    }

    /// Compresses `src`, or returns `None` if the method cannot.
    pub(crate) fn encode(&self, src: &[u8], quality_layout: &QualityLayout<'_>) -> Option<Vec<u8>> {
        encode(&self.options, self.method, self.level, src, quality_layout)
    }

    /// Returns the block of `src`, compressed unless the method is raw, fails, or does not make it
    /// smaller.
    pub(crate) fn compress(&self, src: &[u8], quality_layout: &QualityLayout<'_>) -> Block {
        if self.method != RAW
            && let Some(data) = self.encode(src, quality_layout)
            && data.len() < src.len()
        {
            return block(
                ContentType::ExternalData,
                self.content_id,
                compression_method(self.method),
                src,
                data,
            );
        }

        raw_block(ContentType::ExternalData, self.content_id, src)
    }
}

/// The candidate methods to try on a block (see [`Choice::Trial`]).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Trial {
    options: Options,
    content_id: block::ContentId,
    methods: u32,
    level: u8,
}

/// The compressed sizes of a trial's block, for [`Compressor::end_trial`].
#[derive(Debug)]
pub(crate) struct TrialSizes([f64; METHOD_COUNT]);

impl Trial {
    /// The methods to try, in the order [`Self::pick`] takes their results.
    pub(crate) fn methods(&self) -> impl Iterator<Item = Method> + use<> {
        let Self {
            options,
            content_id,
            methods,
            level,
        } = *self;

        (0..METHOD_COUNT)
            .filter(move |&method| methods & bit(method) != 0)
            .map(move |method| Method {
                options,
                content_id,
                method,
                level,
            })
    }

    /// Returns the block of `src` with the smallest of `results` (what each of [`Self::methods`]
    /// gave, in order), or uncompressed if none is smaller, and the sizes for
    /// [`Compressor::end_trial`].
    pub(crate) fn pick(&self, src: &[u8], results: Vec<Option<Vec<u8>>>) -> (Block, TrialSizes) {
        let mut sizes = [f64::from(u32::MAX); METHOD_COUNT];
        let mut best: (usize, usize, Option<Vec<u8>>) = (RAW, src.len(), None);

        for (method, data) in self.methods().zip(results) {
            if let Some(data) = data {
                sizes[method.method] = data.len() as f64;

                if data.len() < best.1 {
                    best = (method.method, data.len(), Some(data));
                }
            }
        }

        let block = match best {
            (method, _, Some(data)) => block(
                ContentType::ExternalData,
                self.content_id,
                compression_method(method),
                src,
                data,
            ),
            (_, _, None) => raw_block(ContentType::ExternalData, self.content_id, src),
        };

        (block, TrialSizes(sizes))
    }
}

/// What [`begin`] decided for a block.
#[derive(Debug, Eq, PartialEq)]
enum Step {
    /// Use this method, without a trial.
    Use(usize),
    /// Try these methods; [`end`] takes their sizes.
    Trial(u32),
}

/// Decides whether a block of `len` bytes is a trial of `methods` and, if not, which method it
/// uses (the first half of `cram_compress_block3`).
fn begin(metrics: &mut Metrics, mut methods: u32, len: usize) -> Step {
    let size = len as f64;

    // Sudden changes in size (say, from sorted to unsorted data) trigger a new trial.
    if metrics.input_avg_size != 0.0
        && (size / 4.0 - 750.0 > metrics.input_avg_size
            || size < metrics.input_avg_size / 4.0 - 750.0)
        && (size - metrics.input_avg_size).abs() / 10.0 > metrics.input_avg_delta
    {
        metrics.next_trial = 0;
    }

    let trying = metrics.trial > 0 || {
        metrics.next_trial -= 1;
        metrics.next_trial <= 0
    };

    metrics.input_avg_delta =
        0.9 * (metrics.input_avg_delta + (size - metrics.input_avg_size).abs());
    metrics.input_avg_size += size * 0.2;
    metrics.input_avg_size *= 0.8;

    if !trying {
        return Step::Use(metrics.method);
    }

    if metrics.revised_methods != 0 {
        methods = metrics.revised_methods;
    } else {
        metrics.revised_methods = methods;
    }

    if metrics.next_trial <= 0 {
        metrics.next_trial = TRIAL_SPAN;
        metrics.trial = TRIALS;

        for size in &mut metrics.sizes {
            *size /= 2.0;
        }
    }

    Step::Trial(methods)
}

/// Learns from the compressed `sizes` of a trial of `methods` (the second half of
/// `cram_compress_block3`).
fn end(metrics: &mut Metrics, methods: u32, level: u8, sizes: [f64; METHOD_COUNT]) {
    // Small blocks are not to be trusted much: a fixed cost per block evens out small
    // differences.
    for (total, size) in metrics.sizes.iter_mut().zip(sizes) {
        *total += size + 2000.0;
    }

    metrics.trial -= 1;

    if metrics.trial == 0 {
        finish_trial(metrics, methods, level);
    }
}

/// Picks and applies a method for a block of `len` bytes, returning it and the compressed data,
/// or `None` to store the block uncompressed: what [`Compressor`] and [`Trial`] do, in one go.
#[cfg(test)]
fn select<F>(
    metrics: &mut Metrics,
    methods: u32,
    level: u8,
    len: usize,
    mut encode: F,
) -> io::Result<(usize, Option<Vec<u8>>)>
where
    F: FnMut(usize) -> Option<Vec<u8>>,
{
    let methods = match begin(metrics, methods, len) {
        Step::Use(RAW) => return Ok((RAW, None)),
        // A method that fails on this block leaves it uncompressed.
        Step::Use(method) => {
            return Ok(match encode(method) {
                Some(data) if data.len() < len => (method, Some(data)),
                _ => (RAW, None),
            });
        }
        Step::Trial(methods) => methods,
    };

    // Compress this block with every candidate, keeping the smallest.
    let mut sizes = [f64::from(u32::MAX); METHOD_COUNT];
    let mut best: (usize, usize, Option<Vec<u8>>) = (RAW, len, None);

    for (method, size) in sizes.iter_mut().enumerate() {
        if methods & bit(method) == 0 {
            continue;
        }

        if let Some(data) = encode(method) {
            *size = data.len() as f64;

            if data.len() < best.1 {
                best = (method, data.len(), Some(data));
            }
        }
    }

    end(metrics, methods, level, sizes);

    Ok((best.0, best.2))
}

/// Chooses the method to use until the next trial, and drops methods that keep losing.
fn finish_trial(metrics: &mut Metrics, mut methods: u32, level: u8) {
    // Scale the sizes by the methods' costs, depending on the level.
    let weight = match level {
        0..=1 => 4.0,
        2..=3 => 1.0,
        4..=6 => 0.5,
        7 => 1.0 / 3.0,
        _ => 0.0,
    };

    for (size, cost) in metrics.sizes.iter_mut().zip(METHOD_COSTS) {
        *size *= 1.0 + (cost - 1.0) * weight;
    }

    let mut best_method = RAW;
    let mut best_size = f64::from(i32::MAX);

    for (method, &size) in metrics.sizes.iter().enumerate() {
        if size == 0.0 || methods & bit(method) == 0 {
            continue;
        }

        if best_size > size {
            best_size = size;
            best_method = method;
        }
    }

    if best_method == metrics.method {
        metrics.next_trial = (f64::from(metrics.next_trial)
            * (1.0 + f64::from(metrics.consistency) / 4.0).min(2.0))
            as i32;
        metrics.consistency += 1;
    } else {
        metrics.consistency = 0;
    }

    metrics.method = best_method;

    let factor = if level >= 7 { 2 } else { 1 };

    for method in 0..METHOD_COUNT {
        if method == best_method {
            metrics.fails[method] = 0;
            metrics.extra[method] = 0.0;
        } else if best_size < metrics.sizes[method] {
            let ratio = metrics.sizes[method] / best_size - 1.0;

            metrics.fails[method] += 1;

            if metrics.fails[method] >= MAX_FAILS * factor {
                metrics.extra[method] += ratio;

                if metrics.extra[method] >= MAX_DELTA * f64::from(factor) {
                    methods &= !bit(method);
                }
            }

            // fqzcomp rarely changes.
            if method == FQZ {
                methods &= !bit(method);
            }
        }
    }

    metrics.revised_methods = methods;
}

/// Compresses `src` with an HTSlib method number, or returns `None` if the method cannot.
fn encode(
    options: &Options,
    method: usize,
    level: u8,
    src: &[u8],
    quality_layout: &QualityLayout<'_>,
) -> Option<Vec<u8>> {
    let level = u32::from(level);

    let result = match method {
        GZIP => gzip::encode(flate2::Compression::new(level), src),
        GZIP_1 => gzip::encode(flate2::Compression::new(1), src),
        BZIP2 => bzip2::encode(::bzip2::Compression::new(level.clamp(1, 9)), src),
        LZMA => lzma::encode(level, src),
        RANS0 => rans_4x8::encode(rans_4x8::Order::Zero, src),
        RANS1 => rans_4x8::encode(rans_4x8::Order::One, src),
        RANS_PR0 | RANS_PR1..=RANS_PR193 => {
            let mut flags =
                rans_nx16::Flags::from_bits_truncate(pr_flags(method, RANS_PR0, RANS_PR1));

            if src.len() >= RANS_32_STATE_MIN_LEN {
                flags.insert(rans_nx16::Flags::N32);
            }

            rans_nx16::encode(flags, src)
        }
        ARITH_PR0 | ARITH_PR1..=ARITH_PR193 => aac::encode(
            aac::Flags::from_bits_truncate(pr_flags(method, ARITH_PR0, ARITH_PR1)),
            src,
        ),
        TOK3 => name_tokenizer::encode(src),
        FQZ => match quality_layout.read_lengths {
            Some(lens) if options.is_3_1() => fqzcomp::encode(lens, src),
            _ => return None,
        },
        _ => return None,
    };

    result.ok()
}

/// The flags of a `*_PR*` method: order 0 for `pr0`, then [`PR_FLAGS`] from `pr1`.
fn pr_flags(method: usize, pr0: usize, pr1: usize) -> u8 {
    if method == pr0 {
        0
    } else {
        PR_FLAGS[method - pr1 + 1]
    }
}

fn compression_method(method: usize) -> CompressionMethod {
    match method {
        GZIP | GZIP_1 => CompressionMethod::Gzip,
        BZIP2 => CompressionMethod::Bzip2,
        LZMA => CompressionMethod::Lzma,
        RANS0 | RANS1 => CompressionMethod::Rans4x8,
        RANS_PR0 | RANS_PR1..=RANS_PR193 => CompressionMethod::RansNx16,
        ARITH_PR0 | ARITH_PR1..=ARITH_PR193 => CompressionMethod::AdaptiveArithmeticCoding,
        TOK3 => CompressionMethod::NameTokenizer,
        FQZ => CompressionMethod::Fqzcomp,
        _ => CompressionMethod::None,
    }
}

fn block(
    content_type: ContentType,
    content_id: block::ContentId,
    compression_method: CompressionMethod,
    src: &[u8],
    data: Vec<u8>,
) -> Block {
    Block {
        compression_method,
        content_type,
        content_id,
        uncompressed_size: src.len(),
        src: data,
    }
}

fn raw_block(content_type: ContentType, content_id: block::ContentId, src: &[u8]) -> Block {
    block(
        content_type,
        content_id,
        CompressionMethod::None,
        src,
        src.to_vec(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pr_flags() {
        assert_eq!(pr_flags(RANS_PR0, RANS_PR0, RANS_PR1), 0);
        assert_eq!(pr_flags(RANS_PR1, RANS_PR0, RANS_PR1), 1);
        assert_eq!(pr_flags(RANS_PR9, RANS_PR0, RANS_PR1), 9);
        assert_eq!(pr_flags(RANS_PR193, RANS_PR0, RANS_PR1), 193);
        assert_eq!(pr_flags(ARITH_PR193, ARITH_PR0, ARITH_PR1), 193);
    }

    #[test]
    fn test_method_sets() {
        let sets = MethodSets::new(&Options::default());
        assert!(sets.names & bit(TOK3) != 0);
        assert!(sets.names & bit(RANS_PR1) == 0);
        assert!(sets.method & bit(RANS_PR193) != 0 && sets.method & bit(RANS_PR129) == 0);
        assert!(sets.quality_scores & bit(FQZ) == 0);

        let sets = MethodSets::new(&Options::new(Version::new(3, 0), Profile::Normal));
        assert!(sets.method & bit(RANS1) != 0 && sets.method & bit(RANS_PR1) == 0);
        assert!(sets.names & bit(TOK3) == 0);

        let sets = MethodSets::new(&Options::new(Version::new(3, 1), Profile::Archive));
        assert!(sets.quality_scores & bit(FQZ) != 0);
        assert!(sets.method & bit(LZMA) != 0 && sets.method & bit(ARITH_PR193) != 0);
    }

    #[test]
    fn test_select_learns_the_best_method() -> io::Result<()> {
        let mut metrics = Metrics::default();
        let methods = bit(GZIP) | bit(RANS_PR0);

        // RANS_PR0 is always smaller.
        let encode = |method: usize| Some(vec![0; if method == RANS_PR0 { 10 } else { 20 }]);

        for _ in 0..2 {
            let (method, data) = select(&mut metrics, methods, 5, 100, encode)?;
            assert_eq!(method, RANS_PR0);
            assert_eq!(data.map(|d| d.len()), Some(10));
        }

        // After the first trial, only RANS_PR0 is used.
        assert_eq!(metrics.method, RANS_PR0);
        let mut calls = Vec::new();
        select(&mut metrics, methods, 5, 100, |method| {
            calls.push(method);
            encode(method)
        })?;
        assert_eq!(calls, [RANS_PR0]);

        Ok(())
    }

    /// Every method encodes short and awkward inputs to data that decodes back.
    #[test]
    fn test_every_method_round_trips() -> io::Result<()> {
        use crate::cram::codecs::{aac, bzip2, gzip, lzma, rans_4x8, rans_nx16};

        let options = Options::new(Version::new(3, 1), Profile::Archive);
        let methods = (0..METHOD_COUNT).filter(|&m| m != FQZ && m != TOK3);

        let mut inputs: Vec<Vec<u8>> = (1..=40)
            .map(|n| (0..n).map(|i| (i * 7 % 5) as u8).collect())
            .collect();
        inputs.push(vec![b'A'; 1000]);
        inputs.push((0..=255).collect());
        inputs.push((0..70_000u32).map(|i| (i % 17) as u8 + b'!').collect());

        let layout = QualityLayout { read_lengths: None };

        for method in methods {
            for src in &inputs {
                let Some(data) = encode(&options, method, 5, src, &layout) else {
                    continue;
                };

                let decoded = match compression_method(method) {
                    CompressionMethod::Gzip => {
                        let mut dst = vec![0; src.len()];
                        gzip::decode(&data, &mut dst)?;
                        dst
                    }
                    CompressionMethod::Bzip2 => {
                        let mut dst = vec![0; src.len()];
                        bzip2::decode(&data, &mut dst)?;
                        dst
                    }
                    CompressionMethod::Lzma => {
                        let mut dst = vec![0; src.len()];
                        lzma::decode(&data, &mut dst)?;
                        dst
                    }
                    CompressionMethod::Rans4x8 => rans_4x8::decode(&data)?,
                    CompressionMethod::RansNx16 => rans_nx16::decode(&data, src.len())?,
                    CompressionMethod::AdaptiveArithmeticCoding => aac::decode(&data, src.len())?,
                    _ => continue,
                };

                assert_eq!(&decoded, src, "method {method}, {} bytes", src.len());
            }
        }

        Ok(())
    }

    /// A block with a symbol frequent enough to overflow 32-bit frequency scaling (more than
    /// 2^20 occurrences) still round-trips.
    #[test]
    fn test_rans_with_a_dominant_symbol() -> io::Result<()> {
        use crate::cram::codecs::rans_nx16;

        let mut src = vec![b'F'; 2_000_000];
        src.extend(b"#,:".repeat(1000));

        for flags in [0, 1, 4, 5, 193] {
            let flags = rans_nx16::Flags::from_bits_truncate(flags);
            let data = rans_nx16::encode(flags, &src)?;
            assert_eq!(rans_nx16::decode(&data, src.len())?, src, "flags {flags:?}");
        }

        Ok(())
    }

    #[test]
    fn test_select_keeps_raw_when_nothing_is_smaller() -> io::Result<()> {
        let mut metrics = Metrics::default();
        let (method, data) = select(&mut metrics, bit(GZIP), 5, 10, |_| Some(vec![0; 20]))?;
        assert_eq!(method, RAW);
        assert!(data.is_none());
        Ok(())
    }
}
