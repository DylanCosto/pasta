//! `vtools gtcheck`: checks sample identity by the discordance of their genotypes (bcftools
//! gtcheck).
//!
//! Without `-g`, every pair of samples of the query file is compared (the cross-check);
//! with `-g`, every query sample against every sample of the second file, read in step
//! through the synced reader. `-s/-S` narrow either side, `-p/-P` name the pairs instead.

mod kernels;
mod options;
mod report;

use std::{
    fs::File,
    io::{self, IsTerminal, Write},
    process::ExitCode,
    time::Instant,
};

use anyhow::bail;
use clap::ArgMatches;
use flate2::{Compression, write::GzEncoder};
use pasta::variant::{
    Header, HeaderRecordType, Record, TypedValues, Unpack, ValueType,
    encoding::INT32_VECTOR_END,
    number::strtod,
    readlist::{self, ListSource},
    synced::{IndexRequirement, SyncedReaders},
    utils::{self, AcSource},
};

pub use self::options::Args;
use self::{
    kernels::{BATCH_SITES, Batch, Discordance, Genotypes, HweSums, Probs, Sums, Tag},
    options::{FLT_EXCLUDE, Settings, USAGE, parse_options},
    report::Distinctive,
};
use crate::{
    filter::{self, FilterLogic, RecordFilter},
    raw_input::strerror,
};

/// A requested exit status after a message already printed (HTSlib's `exit(1)`).
#[derive(Debug)]
struct ExitStatus(u8);

impl std::fmt::Display for ExitStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "exit {}", self.0)
    }
}

impl std::error::Error for ExitStatus {}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let Some(st) = parse_options(matches, &args)? else {
        eprint!("{USAGE}");
        return Ok(ExitCode::from(1));
    };
    let qry_fname = match args.inputs.first() {
        Some(f) => f.clone(),
        None if !io::stdin().is_terminal() => "-".to_owned(),
        None => {
            eprint!("{USAGE}");
            return Ok(ExitCode::from(1));
        }
    };
    if args.inputs.len() > 1 {
        bail!("Error: too many files given, run with -h for help");
    }
    if st.pair_samples.is_some() {
        if st.gt_samples.is_some() || st.qry_samples.is_some() {
            bail!("The -p/-P option cannot be combined with -s/-S");
        }
        if st.ntop != 0 {
            bail!("The --n-matches option cannot be combined with -p/-P");
        }
    }
    if st.distinctive != 0.0 && st.pair_samples.is_none() {
        bail!("The experimental option --distinctive-sites requires -p/-P");
    }
    if st.hom_only && st.gt_fname.is_none() {
        bail!("The option --homs-only requires --genotypes");
    }
    if st.distinctive != 0.0 && st.gt_err != 0 {
        bail!("The option --distinctive-sites cannot be combined with --error-probability");
    }

    let (mut gc, output) = match Gtcheck::init(&st, &qry_fname) {
        Ok(x) => x,
        Err(e) => return exit_status(e),
    };
    // The report is written as it is made: with every pair of a large cohort it would not
    // fit in memory.
    let mut out = io::BufWriter::with_capacity(OUT_BUFFER, output);
    write_header(&mut out)?;
    if let Err(e) = gc.run(&mut out) {
        // Nothing reaches the output when reading fails, as with bcftools (whose output
        // is still in its buffer when it exits): drop the header and the timing lines.
        drop(out.into_parts());
        return exit_status(e);
    }
    if !st.dry_run {
        gc.report(&mut out)?;
        gc.report_distinctive_sites(&mut out)?;
    }
    out.into_inner()
        .map_err(io::IntoInnerError::into_error)?
        .finish()?;
    Ok(ExitCode::SUCCESS)
}

/// Turns an [`ExitStatus`] into its exit code, passing other errors on.
fn exit_status(e: anyhow::Error) -> anyhow::Result<ExitCode> {
    match e.downcast_ref::<ExitStatus>() {
        Some(s) => Ok(ExitCode::from(s.0)),
        None => Err(e),
    }
}

/// Writes the comment lines naming the program, its command line, and the working
/// directory (bcftools' `print_header`).
fn write_header(out: &mut impl Write) -> io::Result<()> {
    let version = env!("CARGO_PKG_VERSION");
    write!(
        out,
        "# This file was produced by vtools ({version}+pasta-{version}), the command line was:\n# \t vtools gtcheck "
    )?;
    for a in std::env::args_os().skip(2) {
        out.write_all(b" ")?;
        out.write_all(a.as_encoded_bytes())?;
    }
    let cwd = std::env::current_dir().unwrap_or_default();
    write!(
        out,
        "\n# and the working directory was:\n# \t {}\n#\n",
        cwd.display()
    )
}

/// The size of the report's output buffer.
const OUT_BUFFER: usize = 1 << 17;

/// Where the report goes: plain text, or gzip with `-Oz`.
enum Output {
    Plain(Box<dyn Write>),
    Gzip(Box<GzEncoder<Box<dyn Write>>>),
}

impl Output {
    /// Opens `fname` (`-` for standard output) for plain text.
    fn plain(fname: &str) -> anyhow::Result<Self> {
        Ok(Self::Plain(Self::create(fname)?))
    }

    /// Opens `fname` (`-` for standard output) for gzip-compressed text (`-Oz`).
    fn gzip(fname: &str) -> anyhow::Result<Self> {
        let w = Self::create(fname)?;
        Ok(Self::Gzip(Box::new(GzEncoder::new(
            w,
            Compression::default(),
        ))))
    }

    /// Opens the file unbuffered (the report is written through a [`io::BufWriter`]),
    /// with its write errors naming it.
    fn create(fname: &str) -> anyhow::Result<Box<dyn Write>> {
        let inner: Box<dyn Write> = if fname == "-" {
            Box::new(crate::output::stdout())
        } else {
            match File::create(fname) {
                Ok(f) => Box::new(f),
                Err(e) => bail!(
                    "[init_data] Error: cannot write to {fname}: {}",
                    strerror(&e)
                ),
            }
        };
        Ok(Box::new(NamedSink {
            name: fname.to_owned(),
            inner,
        }))
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(w) => w.finish()?.flush(),
        }
    }
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gzip(w) => w.write(buf),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.write_all(buf),
            Self::Gzip(w) => w.write_all(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
        }
    }
}

/// The report's file, whose write errors read "Failed to write to NAME: REASON" wherever
/// they surface (a closed pipe keeps its kind, for `main` to die of SIGPIPE).
struct NamedSink {
    name: String,
    inner: Box<dyn Write>,
}

/// A failed write of [`NamedSink`].
#[derive(Debug)]
struct WriteFailed {
    name: String,
    source: io::Error,
}

impl std::fmt::Display for WriteFailed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Failed to write to {}", self.name)
    }
}

impl std::error::Error for WriteFailed {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.source)
    }
}

impl NamedSink {
    fn failed(&self, source: io::Error) -> io::Error {
        io::Error::new(
            source.kind(),
            WriteFailed {
                name: self.name.clone(),
                source,
            },
        )
    }
}

impl Write for NamedSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.inner.write(buf).map_err(|e| self.failed(e))
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        self.inner.write_all(buf).map_err(|e| self.failed(e))
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush().map_err(|e| self.failed(e))
    }
}

/// A pair of samples given with `-p/-P`, as header indices.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct Pair {
    qry: usize,
    gt: usize,
}

/// Where a kept record lies, for the overlap check of the query and `-g` records.
#[derive(Clone, Copy, Debug)]
struct Locus<'a> {
    /// The contig by name: the two headers can number their contigs differently.
    chrom: Option<&'a [u8]>,
    /// The first and last base, 0-based.
    beg: i64,
    end: i64,
}

impl<'a> Locus<'a> {
    fn new(header: &'a Header, rec: &Record) -> Self {
        Self {
            chrom: header.contig_name(rec.rid()),
            beg: rec.pos(),
            end: rec.pos() + rec.rlen() - 1,
        }
    }

    fn on(&self, chrom: Option<&[u8]>) -> bool {
        self.chrom.is_some() && self.chrom == chrom
    }
}

/// The overlap check of bcftools' `is_input_okay`: which of the kept query and `-g` records to
/// drop because it lies off `chrom` (the contig of the newest kept record) or ends before the
/// other begins, or `None` when the two overlap.
///
/// bcftools compares the header-local contig numbers of the two files, which misses identical
/// sites and pairs a stale record with a site of another contig when the headers list their
/// contigs in different orders; the names are compared here.
fn stale_record(chrom: Option<&[u8]>, recs: [Locus<'_>; 2]) -> Option<usize> {
    let [r0, r1] = recs;
    if !r0.on(chrom) {
        Some(0)
    } else if !r1.on(chrom) {
        Some(1)
    } else if r0.end < r1.beg {
        Some(0)
    } else if r1.end < r0.beg {
        Some(1)
    } else {
        None
    }
}

/// The counters of the INFO lines.
#[derive(Debug, Default)]
struct Counts {
    ncmp: u32,
    no_match: u32,
    not_biallelic: u32,
    mono: u32,
    no_data: u32,
    dip_gt: u32,
    dip_pl: u32,
    filter: u32,
    /// Sites by the tags used, `[query][genotypes]`, PL 0 and GT 1.
    nused: [[u32; 2]; 2],
}

type BoxedFilter = (Box<dyn RecordFilter + Send>, FilterLogic);

/// The state of a run (bcftools' `args_t` after `init_data`).
struct Gtcheck {
    sr: SyncedReaders,
    /// Whether `-g` was given (the second reader).
    two_files: bool,
    qry_filter: Option<BoxedFilter>,
    gt_filter: Option<BoxedFilter>,
    qry_tag: Tag,
    gt_tag: Tag,
    /// The header IDs of PL and GT ([`Tag::index`]) in the query and `-g` headers.
    tag_ids: [[Option<u32>; 2]; 2],
    /// The samples whose values each reader's sites need, when that is few enough of them to
    /// decode them alone.
    need: [Option<Vec<usize>>; 2],
    qry_names: Vec<Box<[u8]>>,
    nqry: usize,
    ngt: usize,
    qry_smpl: Option<Vec<usize>>,
    gt_smpl: Option<Vec<usize>>,
    cross_check: bool,
    pairs: Option<Vec<Pair>>,
    npairs: usize,
    probs: Option<Probs>,
    /// The `-g` genotypes compared (`--homs-only`).
    genotypes: Genotypes,
    keep_ref: bool,
    ntop: usize,
    sort_by_hwe: bool,
    /// `--dry-run`: stop after the first site compared.
    dry_run: bool,
    sums: Sums,
    batch: Option<Batch>,
    threads: usize,
    counts: Counts,
    /// The last records kept of each reader (bcftools' `rec[]`), which can be stale.
    recs: [Record; 2],
    has_rec: [bool; 2],
    arr: [Vec<i32>; 2],
    ac: Vec<i32>,
    warned_dip_gt: bool,
    warned_dip_pl: bool,
    distinctive: Option<Distinctive>,
}

/// Resolves `-s/-S` names to sorted header indices (bcftools' `init_samples`).
fn init_samples(
    list: &str,
    source: ListSource,
    header: &Header,
    fname: &str,
) -> anyhow::Result<Vec<usize>> {
    if list == "-" {
        return Ok((0..header.sample_count()).collect());
    }
    let names = match readlist::read_list(list.as_bytes(), source) {
        Ok(n) if !n.is_empty() => n,
        _ => bail!("Failed to parse {list}"),
    };
    let mut smpl = Vec::with_capacity(names.len());
    for name in &names {
        match header.sample_id(name) {
            Some(i) => smpl.push(i as usize),
            None => bail!(
                "No such sample in {fname}: [{}]",
                String::from_utf8_lossy(name)
            ),
        }
    }
    smpl.sort_unstable();
    if let Some(w) = smpl.windows(2).find(|w| w[0] == w[1]) {
        bail!(
            "Error: the sample \"{}\" is listed twice in {list}",
            String::from_utf8_lossy(&header.samples()[w[0]])
        );
    }
    Ok(smpl)
}

/// The tag a file's genotypes are read from (bcftools' `init_data`): `requested` if given
/// (and present), else `first` or else `second` of those in the header.
fn choose_tag(
    header: &Header,
    fname: &str,
    requested: Option<Tag>,
    first: Tag,
) -> anyhow::Result<Tag> {
    let has = |t: Tag| {
        header
            .id(match t {
                Tag::Gt => b"GT".as_slice(),
                Tag::Pl => b"PL",
            })
            .is_some()
    };
    match requested {
        None => {
            let second = match first {
                Tag::Gt => Tag::Pl,
                Tag::Pl => Tag::Gt,
            };
            if has(first) {
                Ok(first)
            } else if has(second) {
                Ok(second)
            } else {
                bail!("[E::init_data] Neither PL nor GT tag is present in the header of {fname}")
            }
        }
        Some(Tag::Gt) if !has(Tag::Gt) => {
            bail!("[E::init_data] The GT tag is not present in the header of {fname}")
        }
        Some(Tag::Pl) if !has(Tag::Pl) => {
            bail!("[E::init_data] The PL tag is not present in the header of {fname}")
        }
        Some(t) => Ok(t),
    }
}

/// Reads `-p/-P` pairs as header indices, sorted (bcftools' `init_data`).
fn read_pairs(
    spec: &str,
    source: ListSource,
    qry: (&Header, &str),
    gt: (&Header, &str),
) -> anyhow::Result<Vec<Pair>> {
    let items = match readlist::read_list(spec.as_bytes(), source) {
        Ok(n) if !n.is_empty() => n,
        _ => bail!("Failed to parse {spec}"),
    };
    if source == ListSource::Inline && items.len() % 2 != 0 {
        bail!("Expected even number of comma-delimited samples with -p");
    }
    let lookup = |(h, fname): (&Header, &str), name: &[u8]| -> anyhow::Result<usize> {
        match h.sample_id(name) {
            Some(i) => Ok(i as usize),
            None => bail!(
                "No such sample in {fname}: [{}]",
                String::from_utf8_lossy(name)
            ),
        }
    };
    let mut pairs = Vec::new();
    if source == ListSource::File {
        let is_space = |c: &u8| matches!(c, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r');
        for line in &items {
            let Some(cut) = line.iter().position(is_space) else {
                bail!("Could not parse {spec}: {}", String::from_utf8_lossy(line));
            };
            let q = lookup(qry, &line[..cut])?;
            let rest = &line[cut + 1..];
            let skip = rest.iter().take_while(|c| is_space(c)).count();
            let g = lookup(gt, &rest[skip..])?;
            pairs.push(Pair { qry: q, gt: g });
        }
    } else {
        for [q, g] in items.as_chunks::<2>().0 {
            let (q, g) = (lookup(qry, q)?, lookup(gt, g)?);

            pairs.push(Pair { qry: q, gt: g });
        }
    }
    pairs.sort_unstable();
    Ok(pairs)
}

/// The samples whose values are read from each reader, when they are at most half of the
/// file's samples: those of the pairs, or of `-s/-S`. With one file, both sides read it.
fn needed_samples(
    sr: &SyncedReaders,
    pairs: Option<&[Pair]>,
    qry_smpl: Option<&[usize]>,
    gt_smpl: Option<&[usize]>,
) -> [Option<Vec<usize>>; 2] {
    let (qry, gt): (Option<Vec<usize>>, Option<Vec<usize>>) = match pairs {
        Some(p) => (
            Some(p.iter().map(|p| p.qry).collect()),
            Some(p.iter().map(|p| p.gt).collect()),
        ),
        None => (
            qry_smpl.map(<[usize]>::to_vec),
            gt_smpl.map(<[usize]>::to_vec),
        ),
    };
    let few = |list: Option<Vec<usize>>, i: usize| {
        let mut list = list?;
        list.sort_unstable();
        list.dedup();
        (list.len() * 2 <= sr.header(i).sample_count()).then_some(list)
    };
    if sr.reader_count() > 1 {
        [few(qry, 0), few(gt, 1)]
    } else {
        let both = match (qry, gt) {
            (Some(mut q), Some(g)) => {
                q.extend(g);
                Some(q)
            }
            _ => None,
        };
        [few(both, 0), None]
    }
}

/// Reads an integer FORMAT field (GT or PL) as HTSlib's `bcf_get_format_int32` does, but
/// decodes only the values of `samples`: the others are left as they were, and never read.
/// Returns the number of values of all samples, `None` where HTSlib fails.
fn read_samples(
    rec: &mut Record,
    header: &Header,
    id: u32,
    tag: Tag,
    samples: &[usize],
    out: &mut Vec<i32>,
) -> Option<usize> {
    let want = match tag {
        Tag::Gt => ValueType::String,
        Tag::Pl => ValueType::Integer,
    };
    if header.id_type(HeaderRecordType::Format, id) != Some(want) {
        return None;
    }
    rec.unpack(Unpack::FMT);
    let rec = &*rec;
    let field = rec.format_fields().iter().find(|f| f.id() == id as i32)?;
    if field.is_removed() {
        return None;
    }
    let (n, nsmpl) = (field.n(), header.sample_count());
    if n == 0 || nsmpl == 0 {
        return Some(0);
    }
    out.resize(n * nsmpl, INT32_VECTOR_END);
    for &s in samples {
        let values = TypedValues::new(field.type_code(), rec.format_sample_data(field, s));
        let mut end = false;
        for (j, slot) in out[s * n..(s + 1) * n].iter_mut().enumerate() {
            if !end {
                *slot = values.int(j)?;
                end = *slot == INT32_VECTOR_END;
            } else {
                *slot = INT32_VECTOR_END;
            }
        }
    }
    Some(n * nsmpl)
}

/// bcftools' `parse_mem_string`, checked as `extsort` checks `--distinctive-sites`' MEM: the
/// number of bytes.
fn parse_mem_string(s: &str) -> anyhow::Result<usize> {
    let (mem, end) = strtod(s.as_bytes(), 0);
    if end == 0 {
        bail!("Could not parse the memory string: \"{s}\"");
    }
    let suffix = &s[end..];
    let mult = if suffix.eq_ignore_ascii_case("k") {
        1e3
    } else if suffix.eq_ignore_ascii_case("m") {
        1e6
    } else if suffix.eq_ignore_ascii_case("g") {
        1e9
    } else {
        1.0
    };
    // The value is stored in a size_t: anything that truncates to 0 is rejected.
    let bytes = mem * mult;
    if bytes.is_finite() && bytes.trunc() == 0.0 {
        bail!("Could not parse the memory string, expected positive number: {s}");
    }
    // A negative value wraps to a huge size_t, and an infinite one saturates: no limit.
    Ok(if bytes >= 1.0 {
        bytes as usize
    } else {
        usize::MAX
    })
}

/// The samples compared (bcftools' `init_data`): `-s/-S` lists or `-p/-P` pairs as header
/// indices, and the number of pairs.
struct Samples {
    nqry: usize,
    ngt: usize,
    qry_smpl: Option<Vec<usize>>,
    gt_smpl: Option<Vec<usize>>,
    /// Every query sample against every other (no `-g` samples given).
    cross_check: bool,
    pairs: Option<Vec<Pair>>,
    npairs: usize,
}

impl Samples {
    /// Resolves `-s/-S` and `-p/-P` against the query header and the `-g` header (the query
    /// header again without `-g`).
    fn new(st: &Settings, qry: (&Header, &str), gt: (&Header, &str)) -> anyhow::Result<Self> {
        let (qry_hdr, qry_fname) = qry;
        let (gt_hdr, gt_fname) = gt;
        let mut nqry = qry_hdr.sample_count();
        let mut qry_smpl = None;
        if let Some((list, source)) = &st.qry_samples {
            let s = init_samples(list, *source, qry_hdr, qry_fname)?;
            nqry = s.len();
            qry_smpl = Some(s);
        }
        let mut ngt = 0;
        let mut gt_smpl = None;
        let mut pairs = None;
        if let Some((list, source)) = &st.gt_samples {
            let s = init_samples(list, *source, gt_hdr, gt_fname)?;
            ngt = s.len();
            gt_smpl = Some(s);
        } else if let Some(spec) = &st.pair_samples {
            pairs = Some(read_pairs(spec, st.pair_samples_source, qry, gt)?);
        } else if st.gt_fname.is_some() {
            ngt = gt_hdr.sample_count();
        }
        let mut cross_check = false;
        if ngt == 0 {
            ngt = nqry;
            gt_smpl.clone_from(&qry_smpl);
            cross_check = true;
        }
        let npairs = match &pairs {
            Some(p) => p.len(),
            None if cross_check => nqry * (nqry + 1) / 2,
            None => ngt * nqry,
        };
        Ok(Self {
            nqry,
            ngt,
            qry_smpl,
            gt_smpl,
            cross_check,
            pairs,
            npairs,
        })
    }
}

/// Opens the query file and the `-g` file with the regions and targets (bcftools'
/// `init_data`).
fn open_readers(st: &Settings, qry_fname: &str) -> anyhow::Result<SyncedReaders> {
    let mut sr = SyncedReaders::new();
    sr.set_threads(st.threads);
    if let Some(r) = &st.regions {
        sr.set_regions_overlap(st.regions_overlap);
        if sr
            .set_regions_spec(r.as_bytes(), st.regions_source)
            .is_err()
        {
            bail!("Failed to read the regions: {r}");
        }
    }
    if let Some(t) = &st.targets {
        sr.set_targets_overlap(st.targets_overlap);
        if sr
            .set_targets_spec(t.as_bytes(), st.targets_source, None)
            .is_err()
        {
            bail!("Failed to read the targets: {t}");
        }
    }
    if st.gt_fname.is_some() {
        sr.set_index_requirement(IndexRequirement::Required);
    }
    if let Err(e) = sr.add_reader(qry_fname) {
        bail!("Failed to open {qry_fname}: {e}");
    }
    if let Some(g) = &st.gt_fname
        && let Err(e) = sr.add_reader(g)
    {
        let name = if g == "-" { "standard input" } else { g };
        bail!("Failed to read from {name}: {e}");
    }
    Ok(sr)
}

/// The `-i/-e` filters of the query and `-g` files.
fn init_filters(
    st: &Settings,
    sr: &SyncedReaders,
) -> anyhow::Result<(Option<BoxedFilter>, Option<BoxedFilter>)> {
    let logic = |bits: u8| {
        if bits == FLT_EXCLUDE {
            FilterLogic::Exclude
        } else {
            FilterLogic::Include
        }
    };
    let gt_filter = match (&st.gt_filter, st.gt_fname.is_some()) {
        (Some(e), true) => Some((filter::init(sr.header(1), e)?, logic(st.gt_filter_logic))),
        _ => None,
    };
    let qry_filter = match &st.qry_filter {
        Some(e) => Some((filter::init(sr.header(0), e)?, logic(st.qry_filter_logic))),
        None => None,
    };
    Ok((qry_filter, gt_filter))
}

impl Gtcheck {
    /// bcftools' `init_data`: opens the inputs and the output, and sets up the samples and
    /// the sums.
    fn init(st: &Settings, qry_fname: &str) -> anyhow::Result<(Self, Output)> {
        let sr = open_readers(st, qry_fname)?;
        let two_files = st.gt_fname.is_some();
        let gt_fname = st.gt_fname.as_deref().unwrap_or(qry_fname);
        if sr.header(0).sample_count() == 0 {
            bail!("No samples in {qry_fname}?");
        }
        if two_files && sr.header(1).sample_count() == 0 {
            bail!("No samples in {gt_fname}?");
        }
        let (qry_filter, gt_filter) = init_filters(st, &sr)?;

        let qry_tag = choose_tag(sr.header(0), qry_fname, st.qry_use, Tag::Pl)?;
        let gt_tag = if two_files {
            choose_tag(sr.header(1), gt_fname, st.gt_use, Tag::Gt)?
        } else {
            qry_tag
        };

        let qry_hdr = sr.header(0);
        let gt_hdr = if two_files { sr.header(1) } else { qry_hdr };
        let Samples {
            nqry,
            ngt,
            qry_smpl,
            gt_smpl,
            cross_check,
            pairs,
            npairs,
        } = Samples::new(st, (qry_hdr, qry_fname), (gt_hdr, gt_fname))?;

        let probs = (st.gt_err != 0).then(|| Probs::new(st.gt_err));
        let hwe = if st.calc_hwe {
            HweSums::On
        } else {
            HweSums::Off
        };
        let sums = Sums::new(npairs, probs.as_ref(), hwe);
        let distinctive = if st.distinctive == 0.0 {
            None
        } else {
            let mut d = Distinctive::new(st.distinctive, npairs)?;
            d.spill_to(parse_mem_string(&st.es_max_mem)?, &st.es_tmp_prefix);
            Some(d)
        };
        let output = if st.gzip {
            Output::gzip(&st.output_fname)?
        } else {
            Output::plain(&st.output_fname)?
        };
        let batch = pairs.is_none().then(|| {
            if cross_check {
                Batch::triangle(nqry)
            } else {
                Batch::rect(nqry, ngt)
            }
        });
        let ids = |h: &Header| [h.id(b"PL"), h.id(b"GT")];
        let tag_ids = [
            ids(sr.header(0)),
            if two_files {
                ids(sr.header(1))
            } else {
                [None, None]
            },
        ];
        let need = needed_samples(
            &sr,
            pairs.as_deref(),
            qry_smpl.as_deref(),
            gt_smpl.as_deref(),
        );
        let qry_names = sr.header(0).samples().to_vec();
        let gc = Self {
            sr,
            two_files,
            qry_filter,
            gt_filter,
            qry_tag,
            gt_tag,
            tag_ids,
            need,
            qry_names,
            nqry,
            ngt,
            qry_smpl,
            gt_smpl,
            cross_check,
            pairs,
            npairs,
            probs,
            genotypes: if st.hom_only {
                Genotypes::HomsOnly
            } else {
                Genotypes::All
            },
            keep_ref: st.keep_ref,
            ntop: st.ntop,
            sort_by_hwe: st.sort_by_hwe,
            dry_run: st.dry_run,
            sums,
            batch,
            threads: st.threads.max(1),
            counts: Counts::default(),
            recs: [Record::new(), Record::new()],
            has_rec: [false; 2],
            arr: [Vec::new(), Vec::new()],
            ac: Vec::new(),
            warned_dip_gt: false,
            warned_dip_pl: false,
            distinctive,
        };
        Ok((gc, output))
    }

    /// The header index of query sample `i`.
    fn qry_sample(&self, i: usize) -> usize {
        self.qry_smpl.as_ref().map_or(i, |s| s[i])
    }

    /// The index in the sums of query sample `i` against genotype sample `j` (without the
    /// cross-check or pairs).
    fn rect_index(&self, i: usize, j: usize) -> usize {
        self.batch
            .as_ref()
            .expect("the sample matrix has a batch")
            .pair_index(i, j)
    }

    /// The header index of genotype sample `j`.
    fn gt_sample(&self, j: usize) -> usize {
        self.gt_smpl.as_ref().map_or(j, |s| s[j])
    }

    /// Reads every site (bcftools' main loop), writing the timing line to `out`.
    fn run(&mut self, out: &mut impl Write) -> anyhow::Result<()> {
        loop {
            match self.sr.next_line() {
                Ok(0) => break,
                Ok(_) => {}
                Err(e) => bail!("Error: {e}"),
            }
            if !self.input_okay()? {
                continue;
            }
            let t0 = Instant::now();
            self.process_line()?;
            // bcftools prints this again for every site after the first compared one until
            // a second one is compared.
            if self.counts.ncmp == 1 {
                let secs = t0.elapsed().as_secs_f64();
                eprintln!("INFO:\tTime required to process one record .. {secs:.6} seconds");
                writeln!(
                    out,
                    "INFO\tTime required to process one record .. {secs:.6} seconds"
                )?;
                if self.dry_run {
                    break;
                }
            }
        }
        if let Some(e) = self.sr.take_pending_error() {
            bail!("Error: {e}");
        }
        if let Some(b) = &mut self.batch {
            b.flush(&mut self.sums, self.threads);
        }
        Ok(())
    }

    /// Prints bcftools' one-time note about a skipped site of reader `i`.
    fn skip_note(&self, at: (usize, i32, i64), msg: &str) {
        let (i, rid, pos) = at;
        let chrom = self.sr.header(i).contig_name(rid).unwrap_or(b"(null)");
        eprintln!(
            "INFO: skipping {}:{}, {msg}. (This is printed only once.)",
            String::from_utf8_lossy(chrom),
            pos + 1
        );
    }

    /// bcftools' `is_input_okay`: applies the filters, keeps the records, and checks that
    /// they overlap and are biallelic.
    fn input_okay(&mut self) -> anyhow::Result<bool> {
        let nreaders = self.sr.reader_count();
        // The reader of the newest kept record, whose contig is the current site's.
        let mut newest = None;
        let mut skip_filter = false;
        let mut last = (0, -1, 0);
        for i in 0..nreaders {
            let Some((rec, header)) = self.sr.line_with_header(i) else {
                continue;
            };
            last = (i, rec.rid(), rec.pos());
            let flt = if i == 0 {
                &mut self.qry_filter
            } else {
                &mut self.gt_filter
            };
            if let Some((f, logic)) = flt {
                let pass = f.pass(header, rec)? != (*logic == FilterLogic::Exclude);
                if !pass {
                    skip_filter = true;
                    continue;
                }
            }
            newest = Some(i);
            self.sr.swap_line(i, &mut self.recs[i]);
            self.has_rec[i] = true;
        }
        if skip_filter {
            self.counts.filter += 1;
            if self.counts.filter == 1 {
                self.skip_note(last, "filtered site");
            }
            return Ok(false);
        }
        if nreaders == 2 {
            let mut olap = self.has_rec[0] && self.has_rec[1];
            if olap {
                let locus = |i: usize| Locus::new(self.sr.header(i), &self.recs[i]);
                let chrom = newest.and_then(|i| locus(i).chrom);
                if let Some(c) = stale_record(chrom, [locus(0), locus(1)]) {
                    self.has_rec[c] = false;
                    olap = false;
                }
            }
            if !olap {
                self.counts.no_match += 1;
                if self.counts.no_match == 1 {
                    self.skip_note(last, "POS+ALT mismatch");
                }
                return Ok(false);
            }
        }
        let mut nmono = 0;
        for i in 0..nreaders {
            if self.recs[i].allele_count() > 2 {
                self.counts.not_biallelic += 1;
                if self.counts.not_biallelic == 1 {
                    self.skip_note(last, "not a biallelic site, run `vtools norm -m -` first");
                }
                return Ok(false);
            }
            if self.recs[i].variant_types().is_empty() {
                nmono += 1;
            }
        }
        if nmono == nreaders && !self.keep_ref {
            self.counts.mono += 1;
            if self.counts.mono == 1 {
                self.skip_note(last, "monoallelic site");
            }
            return Ok(false);
        }
        Ok(true)
    }

    /// bcftools' `set_data`: reads reader `i`'s genotypes or PLs, falling back to the other
    /// tag when the record lacks the one asked for; the values per sample, or `None` when
    /// the site is skipped.
    fn set_data(&mut self, i: usize, tag: &mut Tag) -> Option<usize> {
        let header = self.sr.header(i);
        let rec = &mut self.recs[i];
        let arr = &mut self.arr[i];
        let nsmpl = header.sample_count();
        for attempt in 0..2 {
            let id = self.tag_ids[i][tag.index()];
            let ret = match (id, &self.need[i]) {
                (None, _) => None,
                (Some(id), None) => rec.get_format_i32_id(header, id, arr).ok(),
                (Some(id), Some(need)) => read_samples(rec, header, id, *tag, need, arr),
            };
            let Some(n) = ret else {
                if attempt == 0 {
                    *tag = match *tag {
                        Tag::Gt => Tag::Pl,
                        Tag::Pl => Tag::Gt,
                    };
                    continue;
                }
                self.counts.no_data += 1;
                return None;
            };
            let per = tag.values();
            if n != per * nsmpl {
                let (warned, name) = match *tag {
                    Tag::Gt => (&mut self.warned_dip_gt, "GT"),
                    Tag::Pl => (&mut self.warned_dip_pl, "PL"),
                };
                if !*warned {
                    *warned = true;
                    let chrom = header.contig_name(rec.rid()).unwrap_or(b"(null)");
                    eprintln!(
                        "INFO: skipping {}:{}, only diploid FORMAT/{name} fields supported. (This is printed only once.)",
                        String::from_utf8_lossy(chrom),
                        rec.pos() + 1
                    );
                }
                match *tag {
                    Tag::Gt => self.counts.dip_gt += 1,
                    Tag::Pl => self.counts.dip_pl += 1,
                }
                return None;
            }
            return Some(per);
        }
        None
    }

    /// The HWE table of the site, from the allele counts of the `-g` record (or the query's).
    fn site_hwe(&mut self) -> anyhow::Result<[f64; 8]> {
        let i = usize::from(self.two_files);
        let header = self.sr.header(i);
        let rec = &mut self.recs[i];
        match utils::calc_ac(header, rec, &mut self.ac, AcSource::InfoThenGenotypes) {
            Ok(true) => {}
            Ok(false) => bail!("todo: bcf_calc_ac() failed"),
            Err(_) => return Err(ExitStatus(1).into()),
        }
        let ac1 = if rec.allele_count() < 2 {
            0
        } else {
            self.ac[1]
        };
        Ok(kernels::hwe_dsg([self.ac[0], ac1]))
    }

    /// bcftools' `process_line`: compares the samples at the current site.
    fn process_line(&mut self) -> anyhow::Result<()> {
        let mut qry_tag = self.qry_tag;
        let mut gt_tag = self.gt_tag;
        let Some(nq1) = self.set_data(0, &mut qry_tag) else {
            return Ok(());
        };
        let ng1 = if self.two_files {
            if let Some((f, logic)) = &mut self.gt_filter {
                let pass = f.pass(self.sr.header(1), &mut self.recs[1])?
                    != (*logic == FilterLogic::Exclude);
                if !pass {
                    self.counts.filter += 1;
                    return Ok(());
                }
            }
            let Some(n) = self.set_data(1, &mut gt_tag) else {
                return Ok(());
            };
            n
        } else {
            nq1
        };
        self.counts.ncmp += 1;
        self.counts.nused[qry_tag.index()][gt_tag.index()] += 1;
        let hwe = if self.sums.hwe.is_some() {
            self.site_hwe()?
        } else {
            [0.0; 8]
        };
        // With one file, the second sample of a pair is read as the first is: bcftools
        // keeps the configured tag there even when the site fell back to the other one,
        // reading GT values as PLs.
        let gt_tag = if self.two_files { gt_tag } else { qry_tag };
        if self.pairs.is_some() {
            return self.add_pairs(qry_tag, nq1, gt_tag, ng1, &hwe);
        }
        let qa = &self.arr[0];
        let ga = if self.two_files { &self.arr[1] } else { qa };
        let qs = self.qry_smpl.as_deref();
        let gs = self.gt_smpl.as_deref();
        let batch = self.batch.as_mut().expect("the sample matrix has a batch");
        batch.push(
            self.probs.as_ref(),
            qry_tag,
            |k| {
                let i = qs.map_or(k, |s| s[k]);
                &qa[i * nq1..(i + 1) * nq1]
            },
            gt_tag,
            |k| {
                let i = gs.map_or(k, |s| s[k]);
                &ga[i * ng1..(i + 1) * ng1]
            },
            self.genotypes,
            hwe,
        );
        if batch.len() == BATCH_SITES {
            batch.flush(&mut self.sums, self.threads);
        }
        Ok(())
    }

    /// The pairs of `-p/-P` at the current site.
    fn add_pairs(
        &mut self,
        qry_tag: Tag,
        nq1: usize,
        gt_tag: Tag,
        ng1: usize,
        hwe: &[f64; 8],
    ) -> anyhow::Result<()> {
        let pairs = self.pairs.as_deref().expect("pairs are set");
        let qa = &self.arr[0];
        let ga = if self.two_files { &self.arr[1] } else { qa };
        let sums = &mut self.sums;
        let genotypes = self.genotypes;
        match (&self.probs, &mut sums.disc) {
            (None, Discordance::Count(ndiff)) => {
                let mut nd = 0;
                let mut ds = self.distinctive.as_mut();
                if let Some(ds) = &mut ds {
                    ds.reset();
                }
                for (i, p) in pairs.iter().enumerate() {
                    let gd = kernels::dsg(gt_tag, &ga[p.gt * ng1..]);
                    if gd == 0 || genotypes.skips(gd) {
                        continue;
                    }
                    let qd = kernels::dsg(qry_tag, &qa[p.qry * nq1..]);
                    if qd == 0 {
                        continue;
                    }
                    let m = qd & gd;
                    if m == 0 {
                        ndiff[i] += 1;
                        if let Some(ds) = &mut ds {
                            nd += 1;
                            ds.insert(i);
                        }
                    } else if let Some((h, nm)) = &mut sums.hwe {
                        h[i] += hwe[usize::from(m)];
                        nm[i] += 1;
                    }
                    sums.ncnt[i] += 1;
                }
                if nd != 0
                    && let Some(ds) = ds
                {
                    let rec = &self.recs[0];
                    ds.push(nd, rec.rid(), rec.pos())?;
                }
            }
            (Some(probs), Discordance::Score(pdiff)) => {
                for (i, p) in pairs.iter().enumerate() {
                    let mut gp = [0.0; 3];
                    let gd = probs.dsg_prob(gt_tag, &ga[p.gt * ng1..], &mut gp);
                    if gd == 0 || genotypes.skips(gd) {
                        continue;
                    }
                    let mut qp = [0.0; 3];
                    let qd = probs.dsg_prob(qry_tag, &qa[p.qry * nq1..], &mut qp);
                    if qd == 0 {
                        continue;
                    }
                    pdiff[i] += kernels::min_sum(qp, gp[0], gp[1], gp[2]);
                    if let Some((h, nm)) = &mut sums.hwe {
                        let m = qd & gd;
                        h[i] += hwe[usize::from(m)];
                        nm[i] += u32::from(m != 0);
                    }
                    sums.ncnt[i] += 1;
                }
            }
            _ => unreachable!("the sums follow --error-probability"),
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use pasta::{
        bgzf,
        variant::{
            index::{
                IndexFormat,
                tabix::{TabixBuilder, TabixConfig},
            },
            io::Reader,
        },
    };

    use super::*;

    fn locus(chrom: &[u8], beg: i64, end: i64) -> Locus<'_> {
        Locus {
            chrom: Some(chrom),
            beg,
            end,
        }
    }

    #[test]
    fn stale_records_are_found_by_contig_name_and_position() {
        let (chr1, chr2) = (&b"chr1"[..], &b"chr2"[..]);
        assert_eq!(
            stale_record(Some(chr1), [locus(chr1, 9, 9), locus(chr1, 9, 9)]),
            None
        );
        assert_eq!(
            stale_record(Some(chr2), [locus(chr2, 9, 9), locus(chr1, 9, 9)]),
            Some(1)
        );
        assert_eq!(
            stale_record(Some(chr2), [locus(chr1, 9, 9), locus(chr2, 9, 9)]),
            Some(0)
        );
        assert_eq!(
            stale_record(Some(chr1), [locus(chr1, 4, 8), locus(chr1, 9, 9)]),
            Some(0)
        );
        assert_eq!(
            stale_record(Some(chr1), [locus(chr1, 9, 9), locus(chr1, 4, 8)]),
            Some(1)
        );
        // A deletion kept from before still overlaps a later SNP.
        assert_eq!(
            stale_record(Some(chr1), [locus(chr1, 4, 12), locus(chr1, 9, 9)]),
            None
        );
        assert_eq!(
            stale_record(None, [locus(chr1, 9, 9), locus(chr1, 9, 9)]),
            Some(0)
        );
    }

    /// A scratch directory of bgzipped, indexed VCFs, removed when dropped.
    struct Files(PathBuf);

    impl Files {
        fn new(name: &str) -> anyhow::Result<Self> {
            let dir =
                std::env::temp_dir().join(format!("vtools-gtcheck-{name}-{}", std::process::id()));
            std::fs::create_dir_all(&dir)?;
            Ok(Self(dir))
        }

        /// Writes `NAME.vcf.gz` and its TBI index with the contigs in the order given, one
        /// sample, and the records as `CHROM POS GT`.
        fn vcf(
            &self,
            name: &str,
            contigs: &[&str],
            sample: &str,
            recs: &[&str],
        ) -> anyhow::Result<String> {
            let mut text = String::from("##fileformat=VCFv4.2\n");
            for c in contigs {
                text += &format!("##contig=<ID={c},length=1000>\n");
            }
            text += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
            text += &format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n");
            for r in recs {
                let [chrom, pos, gt] = r.split(' ').collect::<Vec<_>>()[..] else {
                    anyhow::bail!("bad record {r:?}");
                };
                text += &format!("{chrom}\t{pos}\t.\tA\tG\t.\t.\t.\tGT\t{gt}\n");
            }
            let path = self.0.join(format!("{name}.vcf.gz"));
            let mut w = bgzf::io::Writer::new(File::create(&path)?);
            w.write_all(text.as_bytes())?;
            w.finish()?;
            let fname = path.to_str().expect("UTF-8 path").to_owned();
            index_vcf(&fname, &self.0.join(format!("{name}.vcf.gz.tbi")))?;
            Ok(fname)
        }
    }

    impl Drop for Files {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn index_vcf(fname: &str, idx: &Path) -> anyhow::Result<()> {
        let mut reader = Reader::open(fname, 0)?;
        let mut builder = TabixBuilder::new(TabixConfig::VCF, 0);
        let mut line = Vec::new();
        while reader.read_line(&mut line)? {
            builder.push_line(&line, reader.virtual_position())?;
        }
        builder
            .finish(reader.virtual_position())
            .index
            .save(idx, IndexFormat::Tbi)?;
        Ok(())
    }

    /// Runs `gtcheck -g GT QRY` and returns the counts and the DCv2 line.
    fn gtcheck(qry: &str, gt: &str) -> anyhow::Result<(u32, u32, String)> {
        let st = Settings {
            gt_fname: Some(gt.to_owned()),
            ..Settings::default()
        };
        let (mut gc, _) = Gtcheck::init(&st, qry)?;
        let mut out = Vec::new();
        gc.run(&mut out)?;
        out.clear();
        gc.report(&mut out)?;
        let dc = String::from_utf8(out)?
            .lines()
            .find(|l| l.starts_with("DCv2\t"))
            .unwrap_or_default()
            .to_owned();
        Ok((gc.counts.ncmp, gc.counts.no_match, dc))
    }

    #[test]
    fn identical_sites_match_across_contig_orders() -> anyhow::Result<()> {
        let files = Files::new("order")?;
        let qry = files.vcf("a", &["chr1", "chr2"], "S1", &["chr1 10 0/1"])?;
        let gt = files.vcf("b", &["chr2", "chr1"], "S2", &["chr1 10 0/1"])?;
        let (ncmp, no_match, dc) = gtcheck(&qry, &gt)?;
        assert_eq!((ncmp, no_match), (1, 0), "{dc}");
        assert_eq!(dc, "DCv2\tS1\tS2\t0.000000e+00\t6.931472e-01\t1\t1");
        Ok(())
    }

    #[test]
    fn stale_records_of_another_contig_are_not_paired() -> anyhow::Result<()> {
        let files = Files::new("stale")?;
        let qry = files.vcf(
            "a",
            &["chr1", "chr2"],
            "S1",
            &["chr1 10 0/1", "chr2 10 1/1"],
        )?;
        let gt = files.vcf("b", &["chr2", "chr1"], "S2", &["chr1 10 0/1"])?;
        let (ncmp, no_match, dc) = gtcheck(&qry, &gt)?;
        // chr1:10 is compared; chr2:10 has no -g record, not the stale one of chr1.
        assert_eq!((ncmp, no_match), (1, 1), "{dc}");
        assert_eq!(dc, "DCv2\tS1\tS2\t0.000000e+00\t6.931472e-01\t1\t1");
        Ok(())
    }

    /// A writer that takes `room` bytes and then fails with `kind`.
    struct Cramped {
        room: usize,
        kind: io::ErrorKind,
    }

    impl Write for Cramped {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.room == 0 {
                return Err(io::Error::new(self.kind, "no room"));
            }
            let n = buf.len().min(self.room);
            self.room -= n;
            Ok(n)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn cramped_sink(room: usize, kind: io::ErrorKind) -> io::BufWriter<Output> {
        let sink = NamedSink {
            name: "out.txt".to_owned(),
            inner: Box::new(Cramped { room, kind }),
        };
        io::BufWriter::with_capacity(16, Output::Plain(Box::new(sink)))
    }

    #[test]
    fn write_errors_name_the_output() {
        let mut out = cramped_sink(20, io::ErrorKind::Other);
        let e = out.write_all(&[b'x'; 100]).expect_err("only 20 bytes fit");
        assert_eq!(
            format!("{:#}", anyhow::Error::new(e)),
            "Failed to write to out.txt: no room"
        );
        // A closed pipe keeps its kind, for `main` to die of SIGPIPE.
        let mut out = cramped_sink(0, io::ErrorKind::BrokenPipe);
        let e = out
            .write_all(b"x")
            .and_then(|()| out.flush())
            .expect_err("nothing fits");
        assert_eq!(e.kind(), io::ErrorKind::BrokenPipe);
    }

    #[test]
    fn the_report_is_streamed_and_stops_at_a_write_error() -> anyhow::Result<()> {
        let files = Files::new("stream")?;
        let qry = files.vcf("a", &["chr1"], "S1", &["chr1 10 0/1", "chr1 20 1/1"])?;
        let gt = files.vcf("b", &["chr1"], "S2", &["chr1 10 0/1", "chr1 20 0/1"])?;
        let st = Settings {
            gt_fname: Some(gt),
            ..Settings::default()
        };
        let (mut gc, _) = Gtcheck::init(&st, &qry)?;
        let mut out = Vec::new();
        write_header(&mut out)?;
        gc.run(&mut out)?;
        let mut report = Vec::new();
        gc.report(&mut report)?;
        out.extend_from_slice(&report);
        let text = String::from_utf8(out)?;
        assert!(
            text.starts_with("# This file was produced by vtools ("),
            "{text}"
        );
        let lines: Vec<&str> = text.lines().collect();
        let timing = lines
            .iter()
            .position(|l| l.starts_with("INFO\tTime required to process one record .. "));
        let compared = lines.iter().position(|l| *l == "INFO\tsites-compared\t2");
        assert!(timing.is_some() && timing < compared, "{text}");
        assert_eq!(
            lines.last(),
            Some(&"DCv2\tS1\tS2\t9.210340e+00\t6.931472e-01\t2\t1")
        );

        // The DCv2 line does not fit: the report fails rather than ending quietly.
        let mut out = cramped_sink(report.len() - 10, io::ErrorKind::Other);
        let e = gc
            .report(&mut out)
            .and_then(|()| out.flush())
            .expect_err("the report does not fit");
        assert_eq!(e.to_string(), "Failed to write to out.txt");
        Ok(())
    }
}
