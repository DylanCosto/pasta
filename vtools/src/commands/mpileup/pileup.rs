//! The pileup loop: bcftools' `mplp_func` (read filtering), `pileup_constructor`,
//! `mplp_realn` (partial BAQ), and `mpileup_reg` (calling each position), over HTSlib's
//! pileup engine from `pasta::pileup`.

use std::{
    num::NonZero,
    sync::{Arc, Mutex},
    thread,
};

use crossbeam_channel::{Receiver, Sender};
use pasta::{
    pileup::{Overlaps, Pileup, Read, baq},
    util::threads,
    variant::{Header, Record},
};

use super::{
    call::{self, Call, CallAux, CallRet, CallSettings, Plp, SEQ_NT16_TABLE, cd, glfgen},
    edlib, iaux,
    indel::{self, IndelBuf, RefSeq, Site},
    input::{AlnHeader, AlnIndex, AlnReader, Cursor},
    reference::References,
    regidx::RegIdx,
    sample::Samples,
};

// bcftools' mpileup flags (`MPLP_*`).
pub const NO_ORPHAN: u32 = 1 << 3;
pub const REALN: u32 = 1 << 4;
pub const NO_INDEL: u32 = 1 << 5;
pub const REDO_BAQ: u32 = 1 << 6;
pub const ILLUMINA13: u32 = 1 << 7;
pub const PER_SAMPLE: u32 = 1 << 11;
pub const SMART_OVERLAPS: u32 = 1 << 12;
pub const REALN_PARTIAL: u32 = 1 << 13;

/// Reads per batch sent from a reader thread, at most.
const BATCH_SIZE: usize = 1024;

/// Bytes of encoded reads per batch sent from a reader thread: a batch is cut at its number
/// of reads or at this many bytes, whichever comes first, so that long reads do not fill the
/// queue of each input with gigabytes.
const BATCH_BYTES: usize = 1 << 20;

/// The settings the pileup needs (`mplp_conf_t`).
pub struct Conf {
    pub flag: u32,
    pub min_mq: i32,
    pub capq_thres: i32,
    pub max_depth: usize,
    pub max_indel_depth: i32,
    pub max_read_len: i32,
    pub skip_any_unset: u16,
    pub skip_all_unset: u16,
    pub skip_any_set: u16,
    pub skip_all_set: u16,
    /// `-t/-T`, and whether reads in them are included (not `^`).
    pub bed: Option<RegIdx>,
    pub bed_logic: bool,
    pub call: CallSettings,
}

/// One input file.
pub struct InputFile {
    pub name: String,
    /// The file's ID in [`Samples`].
    pub bam_id: usize,
    /// The reader opened to read the header, kept for a serial run as bcftools keeps it,
    /// and for an input that cannot be opened twice (standard input, a pipe, or a FIFO).
    pub opened: Mutex<Option<(AlnReader, AlnHeader)>>,
    /// The index, loaded when the inputs are opened if regions may be queried, and shared by
    /// every reader of the file (as bcftools keeps each input's index loaded).
    pub index: Option<AlnIndex>,
}

/// What the pileup of a region needs, shared by all threads.
pub struct Shared<'a> {
    pub conf: &'a Conf,
    pub files: &'a [InputFile],
    pub samples: &'a Samples,
    /// The first file's header: contig names for every file.
    pub header: &'a AlnHeader,
    pub refs: &'a References,
    pub out_header: &'a Header,
    pub cram_reference: Option<&'a std::path::Path>,
    /// BGZF decompression threads per input.
    pub io_threads: usize,
    pub nsmpl: usize,
}

/// What a job reads: everything, or a region (0-based, inclusive).
#[derive(Clone, Debug)]
pub enum Span {
    All,
    Region { tid: usize, beg: u32, end: u32 },
}

/// Receives the records of a job.
pub trait Sink {
    /// Takes a record (`flush_bcf_records`).
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be written.
    fn emit(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<()>;
}

/// How the drand48 state is found for a job: known from the start, or waited for at the
/// first draw (see [`RngGate`]).
pub enum RngSource {
    Known,
    Gate(Arc<RngGate>, usize),
}

/// Hands each job the drand48 state after the draws of the jobs before it. A job asks only
/// when it first needs a random number; until then its output does not depend on the state.
pub struct RngGate {
    start: call_rng::Rng48,
    state: std::sync::Mutex<(usize, u64, Vec<Option<u64>>)>,
    cond: std::sync::Condvar,
}

mod call_rng {
    pub use super::super::errmod::Rng48;
}

impl RngGate {
    pub fn new(start: call_rng::Rng48, njobs: usize) -> Self {
        Self {
            start,
            state: std::sync::Mutex::new((0, 0, vec![None; njobs])),
            cond: std::sync::Condvar::new(),
        }
    }

    /// Records that job `i` drew `draws` numbers.
    pub fn finish(&self, i: usize, draws: u64) {
        let mut st = self.state.lock().expect("lock poisoned");
        st.2[i] = Some(draws);
        while st.0 < st.2.len()
            && let Some(d) = st.2[st.0]
        {
            st.1 += d;
            st.0 += 1;
        }
        self.cond.notify_all();
    }

    /// Lets every waiting job go (after an error, when their output is not used).
    pub fn abort(&self) {
        let mut st = self.state.lock().expect("lock poisoned");
        st.0 = st.2.len();
        self.cond.notify_all();
    }

    /// The state job `i` starts with, waiting for the jobs before it.
    pub fn wait(&self, i: usize) -> call_rng::Rng48 {
        let mut st = self.state.lock().expect("lock poisoned");
        while st.0 < i {
            st = self.cond.wait(st).expect("lock poisoned");
        }
        self.start.advanced(st.1)
    }
}

/// A job's calling state (`bcf_callaux_t`, `bcf_callret1_t`, `bcf_call_t`).
pub struct Engine {
    bca: CallAux,
    bcr: Vec<CallRet>,
    bc: Call,
    rec: Record,
    indel_buf: IndelBuf,
    iaux: iaux::Iaux,
    edlib_buf: edlib::EdlibBuf,
    store: Vec<Vec<Plp<'static>>>,
    rng_source: RngSource,
    rng_ready: bool,
    /// Which parts of the carried state were recomputed since it was last set.
    fresh: CarriedUse,
    /// Which parts of the carried state from before that were used.
    used: CarriedUse,
}

/// State bcftools carries from one position to the next besides the pileup: the BQBZ and
/// MQSBZ values indel records repeat from the last SNP with an alternate allele, and the
/// candidate count `--indels-2.0` consensus calls start from.
#[derive(Clone, Copy, Debug, Default)]
pub struct Carried {
    pub bq: f32,
    pub mqs: f32,
    pub ncvar: usize,
}

/// Parts of the [`Carried`] state.
#[derive(Clone, Copy, Debug, Default)]
pub struct CarriedUse {
    pub bias: bool,
    pub ncvar: bool,
}

impl Carried {
    /// Whether the parts in `parts` are equal (the floats as they would be written).
    pub fn same(&self, other: &Self, parts: CarriedUse) -> bool {
        (!parts.bias
            || (self.bq.to_bits() == other.bq.to_bits()
                && self.mqs.to_bits() == other.mqs.to_bits()))
            && (!parts.ncvar || self.ncvar == other.ncvar)
    }

    /// This state with the parts in `parts` taken from `newer`.
    pub fn updated(&self, newer: &Self, parts: CarriedUse) -> Self {
        Self {
            bq: if parts.bias { newer.bq } else { self.bq },
            mqs: if parts.bias { newer.mqs } else { self.mqs },
            ncvar: if parts.ncvar { newer.ncvar } else { self.ncvar },
        }
    }
}

/// Reuses a vector's allocation for references of another lifetime.
fn recycle<'b>(mut v: Vec<Plp<'_>>) -> Vec<Plp<'b>> {
    v.clear();
    v.into_iter()
        .map(|_| unreachable!("the vector is empty"))
        .collect()
}

impl Engine {
    pub fn new(
        settings: &CallSettings,
        nsmpl: usize,
        rng: call_rng::Rng48,
        source: RngSource,
    ) -> Self {
        let rng_ready = matches!(source, RngSource::Known);
        Self {
            bca: CallAux::new(settings, rng),
            bcr: vec![CallRet::default(); nsmpl],
            bc: Call::new(nsmpl),
            rec: Record::new(),
            indel_buf: IndelBuf::default(),
            iaux: iaux::Iaux::default(),
            edlib_buf: edlib::EdlibBuf::default(),
            store: (0..nsmpl).map(|_| Vec::new()).collect(),
            rng_source: source,
            rng_ready,
            fresh: CarriedUse::default(),
            used: CarriedUse::default(),
        }
    }

    /// The state carried from one position to the next.
    pub fn carried(&self) -> Carried {
        Carried {
            bq: self.bc.mwu_bq,
            mqs: self.bc.mwu_mqs,
            ncvar: self.iaux.ncvar(),
        }
    }

    /// Sets the carried state, and starts tracking its use afresh.
    pub fn set_carried(&mut self, carried: Carried) {
        (self.bc.mwu_bq, self.bc.mwu_mqs) = (carried.bq, carried.mqs);
        self.iaux.set_ncvar(carried.ncvar);
        self.fresh = CarriedUse::default();
        self.used = CarriedUse::default();
    }

    /// Which parts of the carried state were recomputed since it was last set.
    pub fn carried_fresh(&self) -> CarriedUse {
        self.fresh
    }

    /// Which parts of the carried state from when it was last set were used.
    pub fn carried_used(&self) -> CarriedUse {
        self.used
    }

    /// Sets the drand48 state to continue with.
    pub fn set_rng(&mut self, rng: call_rng::Rng48) {
        self.bca.rng = rng;
        self.rng_ready = true;
    }

    /// How many random numbers the job drew.
    pub fn draws(&self) -> u64 {
        self.bca.draws
    }

    /// Before a deep pileup is shuffled: makes sure the random state is the right one.
    fn ensure_rng(&mut self) {
        if self.rng_ready {
            return;
        }
        if let RngSource::Gate(gate, i) = &self.rng_source {
            self.bca.rng = gate.wait(*i);
        }
        self.rng_ready = true;
    }
}

/// A prepared read, or the end of a span.
enum Msg {
    Reads(Vec<Read>),
    EndOfSpan,
}

/// `mplp_func` and `pileup_constructor`: filters a read and prepares it for the pileup, or
/// returns `None` to skip it.
fn prepare(raw: Vec<u8>, shared: &Shared<'_>, file: &InputFile) -> anyhow::Result<Option<Read>> {
    let conf = shared.conf;
    let mut b = Read::new(raw)?;
    let flag = b.flags();
    if b.reference_sequence_id() < 0 || flag & 0x4 != 0 {
        return Ok(None);
    }
    if conf.skip_any_unset != 0 && conf.skip_any_unset & flag != conf.skip_any_unset {
        return Ok(None);
    }
    if conf.skip_all_set != 0 && conf.skip_all_set & flag == conf.skip_all_set {
        return Ok(None);
    }
    if conf.skip_all_unset != 0 && conf.skip_all_unset & flag == 0 {
        return Ok(None);
    }
    if conf.skip_any_set != 0 && conf.skip_any_set & flag != 0 {
        return Ok(None);
    }
    // A mapped read without its sequence (SEQ `*`) whose CIGAR consumes query bases: bcftools
    // reads the bases and qualities from past the end of the record, so skip it instead.
    if b.sequence_len() == 0
        && b.cigar()
            .iter()
            .any(|&op| matches!(op & 0xf, 0 | 1 | 4 | 7 | 8) && op >> 4 != 0)
    {
        return Ok(None);
    }
    let tid = b.reference_sequence_id();
    // bcftools takes the contig names of every file from the first file's header and reads
    // past its end for other contigs; stop as it then does.
    if tid as usize >= shared.header.names.len() {
        eprintln!(
            "[E::mplp_func] {}: contig ID {tid} is not in the header of the first input",
            file.name
        );
        return Err(super::exit1("[mpileup] failed to read from input file"));
    }
    if let Some(bed) = &conf.bed {
        let name = &shared.header.names[tid as usize];
        let beg = b.position() as u32;
        let rlen = indel::cigar2rlen(b.cigar()).max(1);
        let end = (b.position() + i64::from(rlen) - 1) as u32;
        let mut overlap = bed.overlap(name, beg, end);
        if !conf.bed_logic {
            // Exclude only reads which are fully contained in an excluded region.
            overlap = !overlap
                || bed
                    .overlaps(name, beg, end)
                    .any(|(rb, re)| beg < rb || end > re);
        }
        if !overlap {
            return Ok(None);
        }
    }
    let sample = shared.samples.sample_id(file.bam_id, &b);
    if sample < 0 {
        return Ok(None);
    }
    if conf.flag & ILLUMINA13 != 0 {
        for q in b.quality_scores_mut() {
            *q = q.saturating_sub(31);
        }
    }
    let reference = if shared.refs.has_fasta() {
        shared.refs.get(tid)
    } else {
        None
    };
    if let Some(r) = &reference {
        if r.len() as i64 <= b.position() {
            eprintln!(
                "[mplp_func] Skipping because {} is outside of {} [ref:{}]",
                b.position(),
                r.len(),
                tid
            );
            return Ok(None);
        }
        if conf.capq_thres > 10 {
            match baq::cap_mapq(&b, r, conf.capq_thres) {
                None => return Ok(None),
                Some(q) if q < 0 => return Ok(None),
                Some(q) => {
                    if i32::from(b.mapping_quality()) > q {
                        b.set_mapping_quality(q as u8);
                    }
                }
            }
        }
    }
    if i32::from(b.mapping_quality()) < conf.min_mq {
        return Ok(None);
    }
    if conf.flag & NO_ORPHAN != 0 && flag & 0x1 != 0 && flag & 0x2 == 0 {
        return Ok(None);
    }
    // pileup_constructor
    let cigar = b.cigar();
    let soft_clip = cigar.iter().any(|&c| c & 0xf == indel::CSOFT_CLIP);
    let has_indel = conf.flag & REALN != 0
        && cigar
            .iter()
            .any(|&c| matches!(c & 0xf, indel::CINS | indel::CDEL | indel::CREF_SKIP));
    let nm = if conf.call.fmt_flag & (call::FMT_NMBZ | call::INFO_NMBZ | call::INFO_NM) != 0 {
        call::read_nm(&b)
    } else {
        -1
    };
    let has_zq = b.tag(*b"ZQ").is_some();
    let flags = [
        (soft_clip, cd::SOFT_CLIP),
        (has_indel, cd::INDEL),
        (has_zq, cd::HAS_ZQ),
    ]
    .iter()
    .filter(|(set, _)| *set)
    .fold(0, |acc, (_, bit)| acc | bit);
    b.set_user_data(cd::pack(sample as usize, flags, nm));
    Ok(Some(b))
}

/// Takes the reader of every input kept open, or opens it again, one input after another
/// before any reader thread starts: besides the descriptors the readers hold, at most one (a
/// CRAM reference's FASTA index) is open at a time, whatever the threads do.
///
/// # Errors
///
/// Returns the first input's error that cannot be opened, as reading it would report it.
fn open_readers(shared: &Shared<'_>) -> anyhow::Result<Vec<(AlnReader, AlnHeader)>> {
    shared
        .files
        .iter()
        .map(|file| {
            let opened = file.opened.lock().expect("input lock poisoned").take();
            match opened {
                Some(opened) => Ok(opened),
                None => AlnReader::open(&file.name, shared.cram_reference, shared.io_threads)
                    .map_err(read_error),
            }
        })
        .collect()
}

/// The error of an input that fails while it is read: bcftools reports HTSlib's error, then
/// exits with status 1.
fn read_error(e: anyhow::Error) -> anyhow::Error {
    if e.downcast_ref::<super::Exit1>().is_some() {
        return e;
    }
    eprintln!("{e:#}");
    super::exit1("[mpileup] failed to read from input file")
}

/// Reads one input's spans with its opened reader, a batch of prepared reads at a time: on a
/// thread of its own ([`Self::send_all`]), or on the thread that piles up when that thread
/// cannot start.
struct InputReader {
    reader: AlnReader,
    header: AlnHeader,
    spans: Vec<Span>,
    /// The span being read, and where in it.
    span: usize,
    cursor: Cursor,
    batch_size: usize,
    /// What follows the batch returned last: the end of its span, or the error that stopped
    /// it.
    after: Option<anyhow::Result<Msg>>,
    /// Whether reading has ended (after an error).
    ended: bool,
}

impl InputReader {
    fn new((reader, header): (AlnReader, AlnHeader), spans: Vec<Span>, batch_size: usize) -> Self {
        Self {
            reader,
            header,
            spans,
            span: 0,
            cursor: Cursor::default(),
            batch_size,
            after: None,
            ended: false,
        }
    }

    /// Sends every message to `tx`, until the receiver is gone.
    fn send_all(mut self, shared: &Shared<'_>, file: &InputFile, tx: &Sender<anyhow::Result<Msg>>) {
        while let Some(msg) = self.next_msg(shared, file) {
            if tx.send(msg).is_err() {
                return;
            }
        }
    }

    /// The next batch of reads, or the end of a span, or the error that ends reading (`None`
    /// after the last span or an error). The reads before a read that cannot be prepared are
    /// piled up first; those of the batch being filled when the input cannot be read are not.
    fn next_msg(&mut self, shared: &Shared<'_>, file: &InputFile) -> Option<anyhow::Result<Msg>> {
        if let Some(msg) = self.after.take() {
            return Some(msg);
        }
        if self.ended {
            return None;
        }
        let span = self.spans.get(self.span)?;
        let batch_size = self.batch_size;
        let mut batch = Vec::with_capacity(batch_size);
        let mut bytes = 0;
        let mut unprepared = false;
        let mut visit = |raw: Vec<u8>| -> anyhow::Result<bool> {
            let prepared = prepare(raw, shared, file).inspect_err(|_| unprepared = true)?;
            if let Some(read) = prepared {
                bytes += read.as_bytes().len();
                batch.push(read);
                // A full batch pauses reading.
                return Ok(batch.len() < batch_size && bytes < BATCH_BYTES);
            }
            Ok(true)
        };
        let read = match span {
            Span::All => self
                .reader
                .read_all(&self.header, &mut self.cursor, &mut visit),
            // As in bcftools, the contig ID of the first file's header is used for every file;
            // an ID past this file's contigs has no reads (HTSlib's index query).
            Span::Region { tid, .. } if *tid >= self.header.names.len() => Ok(true),
            Span::Region { tid, beg, end } => self.reader.read_region(
                &self.header,
                file.index.as_ref(),
                *tid,
                *beg,
                *end,
                &mut self.cursor,
                &mut visit,
            ),
        };
        match read {
            Ok(false) => Some(Ok(Msg::Reads(batch))),
            Ok(true) => {
                self.span += 1;
                self.cursor = Cursor::default();
                if batch.is_empty() {
                    return Some(Ok(Msg::EndOfSpan));
                }
                self.after = Some(Ok(Msg::EndOfSpan));
                Some(Ok(Msg::Reads(batch)))
            }
            Err(e) => {
                self.ended = true;
                if !unprepared || batch.is_empty() {
                    return Some(Err(e));
                }
                self.after = Some(Err(e));
                Some(Ok(Msg::Reads(batch)))
            }
        }
    }
}

/// Where an input's reads come from.
enum Feed {
    /// The input's reader thread.
    Thread(Receiver<anyhow::Result<Msg>>),
    /// The input, read on this thread when the reader thread cannot start.
    Here(Box<InputReader>),
}

impl Feed {
    /// Reads `input` on a thread started in `scope`, with `queue` batches read ahead, or on
    /// the calling thread when the thread cannot start. `closed` is dropped once the thread
    /// has closed the input.
    fn start<'scope, 'env>(
        scope: &'scope thread::Scope<'scope, 'env>,
        shared: &'env Shared<'env>,
        file: &'env InputFile,
        input: InputReader,
        queue: usize,
        closed: Option<Sender<()>>,
    ) -> Self {
        let Ok(mut idle) = threads::start_idle_scoped(scope, NonZero::<usize>::MIN) else {
            return Self::Here(Box::new(input));
        };
        let thread = idle.pop().expect("one thread was started");
        let (tx, rx) = crossbeam_channel::bounded(queue);
        thread.run(move || {
            // Dropped last, after the input is closed.
            let _closed = closed;
            input.send_all(shared, file, &tx);
        });
        Self::Thread(rx)
    }
}

/// A pileup error, reported as HTSlib does before bcftools gives up with status 1.
fn plp_error(e: pasta::pileup::Error) -> anyhow::Error {
    // HTSlib's messages start with a capital letter ("The input is not sorted ...").
    let msg = e.to_string();
    let mut chars = msg.chars();
    let msg: String = chars
        .next()
        .map(|c| c.to_ascii_uppercase())
        .into_iter()
        .chain(chars)
        .collect();
    eprintln!("[E::bam_plp_push] {msg}");
    super::exit1("[mpileup] failed to read from input file")
}

/// A pileup over one input, fed by its reader thread (or read on this one).
struct Source {
    pileup: Pileup,
    /// The input's index in [`Shared::files`].
    input: usize,
    feed: Feed,
    pending: std::vec::IntoIter<Read>,
    span_done: bool,
    /// A key of the last read pushed (its name, position, and flags), which pins where in
    /// the input the pileup is.
    last_key: u64,
}

/// The key of a read for [`Source::last_key`].
fn read_key(b: &Read) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::hash::DefaultHasher::new();
    (b.name(), b.position(), b.flags(), b.reference_sequence_id()).hash(&mut h);
    h.finish()
}

impl Source {
    fn new(shared: &Shared<'_>, input: usize, feed: Feed) -> Self {
        Self {
            pileup: new_pileup(shared.conf),
            input,
            feed,
            pending: Vec::new().into_iter(),
            span_done: false,
            last_key: 0,
        }
    }

    fn next_read(&mut self, shared: &Shared<'_>) -> anyhow::Result<Option<Read>> {
        loop {
            if let Some(r) = self.pending.next() {
                return Ok(Some(r));
            }
            if self.span_done {
                return Ok(None);
            }
            let msg = match &mut self.feed {
                Feed::Thread(rx) => rx.recv().ok(),
                Feed::Here(input) => input.next_msg(shared, &shared.files[self.input]),
            };
            match msg {
                Some(Ok(Msg::Reads(v))) => self.pending = v.into_iter(),
                Some(Ok(Msg::EndOfSpan)) => {
                    self.span_done = true;
                    return Ok(None);
                }
                Some(Err(e)) => return Err(read_error(e)),
                None => return Err(super::exit1("[mpileup] failed to read from input file")),
            }
        }
    }

    /// `bam_plp64_auto`: whether a position is ready.
    fn auto(&mut self, shared: &Shared<'_>) -> anyhow::Result<bool> {
        if self.pileup.next_position().map_err(plp_error)? {
            return Ok(true);
        }
        if self.pileup.is_eof() {
            return Ok(false);
        }
        loop {
            match self.next_read(shared)? {
                Some(read) => {
                    self.last_key = read_key(&read);
                    self.pileup.push(Some(read)).map_err(plp_error)?;
                    if self.pileup.next_position().map_err(plp_error)? {
                        return Ok(true);
                    }
                }
                None => {
                    self.pileup.push(None).map_err(plp_error)?;
                    return self.pileup.next_position().map_err(plp_error);
                }
            }
        }
    }
}

/// `bam_mplp_t`: the inputs' pileups at one position.
struct Multi {
    min_tid: i32,
    min_pos: i64,
    tids: Vec<i32>,
    positions: Vec<i64>,
    ready: Vec<bool>,
    active: Vec<bool>,
}

impl Multi {
    fn new(n: usize) -> Self {
        Self {
            min_tid: -1,
            min_pos: i64::MAX,
            tids: vec![-1; n],
            positions: vec![i64::MAX; n],
            ready: vec![false; n],
            active: vec![false; n],
        }
    }

    /// `bam_mplp64_auto`
    fn next(
        &mut self,
        shared: &Shared<'_>,
        sources: &mut [Source],
    ) -> anyhow::Result<Option<(i32, i64)>> {
        let mut new_min_tid = u32::MAX;
        let mut new_min_pos = i64::MAX;
        for (i, s) in sources.iter_mut().enumerate() {
            if self.positions[i] == self.min_pos && self.tids[i] == self.min_tid {
                self.ready[i] = s.auto(shared)?;
                if self.ready[i] {
                    let (tid, pos) = s.pileup.position();
                    self.tids[i] = tid;
                    self.positions[i] = pos;
                } else {
                    self.tids[i] = 0;
                    self.positions[i] = 0;
                }
            }
            if self.ready[i] {
                let tid = self.tids[i] as u32;
                if tid < new_min_tid {
                    new_min_tid = tid;
                    new_min_pos = self.positions[i];
                } else if tid == new_min_tid && self.positions[i] < new_min_pos {
                    new_min_pos = self.positions[i];
                }
            }
        }
        self.min_pos = new_min_pos;
        self.min_tid = new_min_tid as i32;
        if new_min_pos == i64::MAX {
            return Ok(None);
        }
        for i in 0..self.active.len() {
            self.active[i] =
                self.ready[i] && self.positions[i] == self.min_pos && self.tids[i] == self.min_tid;
        }
        Ok(Some((self.min_tid, self.min_pos)))
    }
}

/// `mplp_realn`: applies BAQ to the reads at a position with an indel nearby, unless the
/// position looks well aligned (partial BAQ).
fn realn(sources: &mut [Source], active: &[bool], conf: &Conf, reference: &[u8]) {
    let flag = conf.flag;
    let (mut has_indel, mut has_clip, mut nt) = (0i32, 0i32, 0i32);
    let (mut min_indel, mut max_indel) = (i32::MAX, i32::MIN);
    for (s, &a) in sources.iter().zip(active) {
        if !a {
            continue;
        }
        nt += s.pileup.entries().len() as i32;
        for e in s.pileup.entries() {
            let d = s.pileup.read(e).user_data();
            has_indel += i32::from(d & cd::INDEL != 0 || e.indel != 0);
            has_clip += i32::from(d & cd::SOFT_CLIP != 0);
            max_indel = max_indel.max(e.indel);
            min_indel = min_indel.min(e.indel);
        }
    }
    let ntf = f64::from(nt);
    if flag & REALN_PARTIAL != 0
        && (has_indel == 0
            || (f64::from(has_clip) < 0.2 * ntf
                && max_indel == min_indel
                && (f64::from(has_indel) < 0.1 * ntf || has_indel == 1)))
    {
        return;
    }
    let realn_dist = 40 + 10 * i32::from(nt < 40) + 10 * i32::from(nt < 20);
    let mode = baq::BaqMode {
        apply: true,
        extend: true,
        redo: flag & REDO_BAQ != 0,
    };
    for (s, &a) in sources.iter_mut().zip(active) {
        if !a {
            continue;
        }
        for j in 0..s.pileup.entries().len() {
            let e = s.pileup.entries()[j];
            let b = s.pileup.read_mut(&e);
            let d = b.user_data();
            if d & cd::REALN != 0 {
                continue;
            }
            b.set_user_data(d | cd::REALN);
            let l_qseq = b.sequence_len() as i32;
            if l_qseq > conf.max_read_len {
                continue;
            }
            let cig = b.cigar();
            let ncig = cig.len();
            if flag & REALN_PARTIAL != 0 && nt > 15 && ncig > 1 {
                let lr = l_qseq > 500;
                let (mut lm, mut rm, mut nm) = (0i32, 0i32, 0usize);
                for &c in cig {
                    let op = c & 0xf;
                    if lr && (op == indel::CHARD_CLIP || op == indel::CSOFT_CLIP) {
                        continue;
                    }
                    if indel::is_match(op) {
                        lm += (c >> 4) as i32;
                        nm += 1;
                    } else {
                        break;
                    }
                }
                if nm != ncig {
                    for &c in cig.iter().rev() {
                        let op = c & 0xf;
                        if lr && (op == indel::CHARD_CLIP || op == indel::CSOFT_CLIP) {
                            continue;
                        }
                        if indel::is_match(op) {
                            rm += (c >> 4) as i32;
                        } else {
                            break;
                        }
                    }
                    if lm >= realn_dist * 4 && rm >= realn_dist * 4 {
                        continue;
                    }
                    if lm >= realn_dist
                        && rm >= realn_dist
                        && f64::from(has_clip) < (0.15 + 0.05 * f64::from(i32::from(nt > 20))) * ntf
                    {
                        continue;
                    }
                }
            }
            if l_qseq > 500 {
                let rl = indel::cigar2rlen(cig);
                if i64::from((rl - l_qseq).abs()) * i64::from(l_qseq) >= 500_000 {
                    continue;
                }
            }
            let _ = baq::prob_realn(b, reference, mode);
            if b.tag(*b"ZQ").is_some() {
                b.set_user_data(b.user_data() | cd::HAS_ZQ);
            }
        }
    }
}

/// Runs the pileup over `spans`, emitting records to `sink` (`mpileup_reg` for each span).
///
/// # Errors
///
/// Returns an error if an input cannot be read or a record written.
pub fn run(
    shared: &Shared<'_>,
    spans: &[(Span, u32, u32)],
    engine: &mut Engine,
    sink: &mut dyn Sink,
) -> anyhow::Result<()> {
    let only: Vec<Span> = spans.iter().map(|(s, _, _)| s.clone()).collect();
    let readers = open_readers(shared)?;
    thread::scope(|scope| -> anyhow::Result<()> {
        let mut sources = Vec::with_capacity(shared.files.len());
        for (i, (file, reader)) in shared.files.iter().zip(readers).enumerate() {
            let input = InputReader::new(reader, only.clone(), BATCH_SIZE);
            let feed = Feed::start(scope, shared, file, input, 8, None);
            sources.push(Source::new(shared, i, feed));
        }
        let result = (|| {
            for (i, (_, beg, end)) in spans.iter().enumerate() {
                if i > 0 {
                    // bam_mplp_reset
                    for s in &mut sources {
                        s.pileup = new_pileup(shared.conf);
                        s.span_done = false;
                        s.pending = Vec::new().into_iter();
                    }
                }
                region(shared, &mut sources, *beg, *end, engine, sink)?;
                // Drain the rest of the span.
                for s in &mut sources {
                    while s.next_read(shared)?.is_some() {}
                }
            }
            Ok(())
        })();
        drop(sources);
        result
    })
}

/// Where a [`Walker`] works: the region's bounds (0-based, inclusive), where BAQ starts being
/// applied, where records start being written, and the position to stop before.
#[derive(Clone, Copy, Debug)]
pub struct WalkRange {
    pub beg: u32,
    pub end: u32,
    pub realn_from: i64,
    pub out_from: i64,
    pub stop_at: i64,
    /// Where positions before `out_from` start being called without output, to find the
    /// BQBZ and MQSBZ values indel records repeat (`i64::MAX` for none).
    pub warm_from: i64,
}

/// A resumable pileup over one span (`mpileup_reg`), which can stop before a position and
/// carry on later, possibly on another thread.
pub struct Walker {
    sources: Vec<Source>,
    multi: Multi,
    /// A position prepared but not yet processed.
    pending: Option<(i32, i64)>,
    done: bool,
    /// Disconnected by each reader thread once it has closed its input.
    closed: Vec<Receiver<()>>,
}

/// Reads per batch, and batches queued, for the readers of a [`Walker`] (which may be paused
/// while holding them).
const WALKER_BATCH: usize = 64;
const WALKER_QUEUE: usize = 4;

impl Walker {
    /// Opens the inputs and starts reading `span` on reader threads spawned in `scope`.
    ///
    /// # Errors
    ///
    /// Returns an error if an input cannot be opened.
    pub fn start<'scope, 'env>(
        scope: &'scope thread::Scope<'scope, 'env>,
        shared: &'env Shared<'env>,
        span: Span,
    ) -> anyhow::Result<Self> {
        let readers = open_readers(shared)?;
        let mut sources = Vec::with_capacity(shared.files.len());
        let mut closed = Vec::with_capacity(shared.files.len());
        for (i, (file, reader)) in shared.files.iter().zip(readers).enumerate() {
            let (closed_tx, closed_rx) = crossbeam_channel::bounded::<()>(0);
            let input = InputReader::new(reader, vec![span.clone()], WALKER_BATCH);
            let feed = Feed::start(scope, shared, file, input, WALKER_QUEUE, Some(closed_tx));
            if matches!(feed, Feed::Thread(_)) {
                closed.push(closed_rx);
            }
            sources.push(Source::new(shared, i, feed));
        }
        Ok(Self {
            multi: Multi::new(sources.len()),
            sources,
            pending: None,
            done: false,
            closed,
        })
    }

    /// Stops the reader threads and waits for them to close the inputs, so that their file
    /// descriptors are free for the next walker. A reader stops at its next batch.
    pub fn close(self) {
        let Self {
            sources, closed, ..
        } = self;
        drop(sources);
        for rx in closed {
            // Nothing is sent: this returns once the reader thread has finished.
            let _ = rx.recv();
        }
    }

    /// Whether the span has been read to its end.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Processes positions until the first at or after `range.stop_at`, which is kept for
    /// the next call, or the end of the span.
    ///
    /// # Errors
    ///
    /// Returns an error if an input cannot be read or a record written.
    pub fn walk(
        &mut self,
        shared: &Shared<'_>,
        engine: &mut Engine,
        sink: &mut dyn Sink,
        range: &WalkRange,
    ) -> anyhow::Result<()> {
        let conf = shared.conf;
        loop {
            let (tid, pos) = match self.pending.take() {
                Some(c) => c,
                None => match self.multi.next(shared, &mut self.sources)? {
                    Some(c) => c,
                    None => {
                        self.done = true;
                        return Ok(());
                    }
                },
            };
            if pos >= range.stop_at {
                self.pending = Some((tid, pos));
                return Ok(());
            }
            if pos < i64::from(range.beg) || pos > i64::from(range.end) {
                continue;
            }
            if let Some(bed) = &conf.bed
                && tid >= 0
            {
                let mut overlap =
                    bed.overlap(&shared.header.names[tid as usize], pos as u32, pos as u32);
                if !conf.bed_logic {
                    overlap = !overlap;
                }
                if !overlap {
                    continue;
                }
            }
            if pos < range.realn_from {
                continue;
            }
            let reference = shared.refs.get(tid);
            if let Some(r) = &reference
                && conf.flag & REALN != 0
            {
                realn(&mut self.sources, &self.multi.active, conf, r);
            }
            if pos < range.out_from {
                if pos >= range.warm_from {
                    warm(
                        shared,
                        &self.sources,
                        &self.multi.active,
                        pos,
                        reference.as_deref().map(Vec::as_slice),
                        engine,
                    );
                }
                continue;
            }
            column(
                shared,
                &self.sources,
                &self.multi.active,
                tid,
                pos,
                reference.as_deref().map(Vec::as_slice),
                engine,
                sink,
            )?;
        }
    }

    /// A hash of everything that decides the output from here on: each input's pileup and
    /// the last read it took, and the merge state. Two walkers over the same span with equal
    /// fingerprints write the same records from here on.
    pub fn fingerprint(&self) -> u128 {
        use std::hash::{Hash, Hasher};
        let mut halves = [0u64; 2];
        for (i, half) in halves.iter_mut().enumerate() {
            let mut h = std::hash::DefaultHasher::new();
            i.hash(&mut h);
            self.pending.hash(&mut h);
            self.done.hash(&mut h);
            let m = &self.multi;
            (
                m.min_tid,
                m.min_pos,
                &m.tids,
                &m.positions,
                &m.ready,
                &m.active,
            )
                .hash(&mut h);
            for s in &self.sources {
                s.pileup.hash_state(&mut h);
                s.last_key.hash(&mut h);
            }
            *half = h.finish();
        }
        u128::from(halves[0]) << 64 | u128::from(halves[1])
    }
}

fn new_pileup(conf: &Conf) -> Pileup {
    let overlaps = if conf.flag & SMART_OVERLAPS != 0 {
        Overlaps::Adjust
    } else {
        Overlaps::Ignore
    };
    Pileup::new(conf.max_depth, overlaps)
}

/// `mpileup_reg`: calls every position of one span within `[beg, end]`.
fn region(
    shared: &Shared<'_>,
    sources: &mut [Source],
    beg: u32,
    end: u32,
    engine: &mut Engine,
    sink: &mut dyn Sink,
) -> anyhow::Result<()> {
    let conf = shared.conf;
    let mut multi = Multi::new(sources.len());
    while let Some((tid, pos)) = multi.next(shared, sources)? {
        if pos < i64::from(beg) || pos > i64::from(end) {
            continue;
        }
        if let Some(bed) = &conf.bed
            && tid >= 0
        {
            let mut overlap =
                bed.overlap(&shared.header.names[tid as usize], pos as u32, pos as u32);
            if !conf.bed_logic {
                overlap = !overlap;
            }
            if !overlap {
                continue;
            }
        }
        let reference = shared.refs.get(tid);
        if let Some(r) = &reference
            && conf.flag & REALN != 0
        {
            realn(sources, &multi.active, conf, r);
        }
        column(
            shared,
            sources,
            &multi.active,
            tid,
            pos,
            reference.as_deref().map(Vec::as_slice),
            engine,
            sink,
        )?;
    }
    Ok(())
}

/// The active inputs' reads at a position, by sample (in vectors from the engine's store),
/// and their number.
#[inline]
fn gather<'s>(
    sources: &'s [Source],
    active: &[bool],
    engine: &mut Engine,
) -> (Vec<Vec<Plp<'s>>>, i32) {
    let mut gplp: Vec<Vec<Plp<'s>>> = std::mem::take(&mut engine.store)
        .into_iter()
        .map(recycle)
        .collect();
    let mut total_depth = 0i32;
    for (s, &a) in sources.iter().zip(active) {
        if !a {
            continue;
        }
        total_depth += s.pileup.entries().len() as i32;
        for e in s.pileup.entries() {
            let b = s.pileup.read(e);
            let d = b.user_data();
            gplp[cd::sample(d)].push(Plp {
                b,
                qpos: e.qpos as i32,
                indel: e.indel,
                is_del: e.is_del,
                is_refskip: e.is_refskip,
                cd: d,
                aux: 0,
            });
        }
    }
    (gplp, total_depth)
}

/// The reference base at `pos` in HTSlib's 4-bit encoding: `N` without a reference or past
/// its end.
#[inline]
fn ref_base16(reference: Option<&[u8]>, pos: i64) -> u8 {
    let ref0 = match reference {
        Some(r) if (pos as usize) < r.len() => r[pos as usize],
        _ => b'N',
    };
    SEQ_NT16_TABLE[ref0 as usize]
}

/// Calls the SNP of a position without writing it, for the BQBZ and MQSBZ values that indel
/// records repeat. The random state is left alone: subsampling does not change them.
fn warm(
    shared: &Shared<'_>,
    sources: &[Source],
    active: &[bool],
    pos: i64,
    reference: Option<&[u8]>,
    engine: &mut Engine,
) {
    let (mut gplp, total_depth) = gather(sources, active, engine);
    let ref16 = ref_base16(reference, pos);
    let Engine { bca, bcr, bc, .. } = engine;
    let (rng, draws) = (bca.rng, bca.draws);
    call::callaux_clean(bca, bc, bcr);
    for (g, r) in gplp.iter().zip(bcr.iter_mut()) {
        glfgen(g, Some(ref16), bca, r);
    }
    call::combine(bcr, bca, Some(ref16), bc);
    (bca.rng, bca.draws) = (rng, draws);
    if bc.has_alt {
        engine.fresh.bias = true;
    }
    // `--indels-2.0` consensus calls carry their candidate count to the next call.
    let conf = shared.conf;
    if bca_is_iaux(engine)
        && conf.flag & NO_INDEL == 0
        && total_depth < conf.max_indel_depth
        && let Some(r) = reference
    {
        let Engine {
            bca, bc, bcr, iaux, ..
        } = engine;
        call::callaux_clean(bca, bc, bcr);
        let mut site = Site {
            plp: &mut gplp,
            pos,
            reference: RefSeq(r),
        };
        let tid_name: &[u8] = b"";
        // Only the candidate count matters here; the chunk that owns the position reports
        // its errors and messages.
        let _ = iaux.try_gap_prep(&mut site, bca, tid_name);
    }
    engine.store = gplp.into_iter().map(recycle).collect();
}

/// Whether the engine uses the `--indels-2.0` model.
fn bca_is_iaux(engine: &Engine) -> bool {
    engine.bca.indels_v20 && !engine.bca.edlib
}

/// Calls one position: the SNP record, then the indel record if there is one.
#[expect(clippy::too_many_arguments, reason = "the position's context")]
fn column(
    shared: &Shared<'_>,
    sources: &[Source],
    active: &[bool],
    tid: i32,
    pos: i64,
    reference: Option<&[u8]>,
    engine: &mut Engine,
    sink: &mut dyn Sink,
) -> anyhow::Result<()> {
    let conf = shared.conf;
    let hdr = shared.out_header;
    let (mut gplp, total_depth) = gather(sources, active, engine);
    let ref16 = ref_base16(reference, pos);
    let need_rng = |engine: &Engine, gplp: &[Vec<Plp<'_>>]| -> bool {
        !engine.rng_ready && gplp.iter().any(|g| g.len() > 255)
    };
    if need_rng(engine, &gplp) {
        engine.ensure_rng();
    }
    let Engine { bca, bcr, bc, .. } = engine;
    call::callaux_clean(bca, bc, bcr);
    for (g, r) in gplp.iter().zip(bcr.iter_mut()) {
        glfgen(g, Some(ref16), bca, r);
    }
    bc.tid = tid;
    bc.pos = pos;
    call::combine(bcr, bca, Some(ref16), bc);
    if bc.has_alt {
        engine.fresh.bias = true;
    }
    let Engine { bca, bc, rec, .. } = engine;
    call::call2bcf(bc, rec, hdr, conf.call.fmt_flag, bca, &[])?;
    sink.emit(hdr, rec)?;

    if conf.flag & NO_INDEL == 0 && total_depth < conf.max_indel_depth {
        let Engine {
            bca,
            bcr,
            bc,
            rec,
            indel_buf,
            iaux,
            edlib_buf,
            ..
        } = engine;
        call::callaux_clean(bca, bc, bcr);
        let ok = match reference {
            None => false,
            Some(r) => {
                let mut site = Site {
                    plp: &mut gplp,
                    pos,
                    reference: RefSeq(r),
                };
                if bca.edlib {
                    edlib::gap_prep(&mut site, bca, edlib_buf)
                } else if bca.indels_v20 {
                    let calls = iaux.calls();
                    let ok = iaux.gap_prep(&mut site, bca, &shared.header.names[tid as usize])?;
                    if iaux.calls() != calls {
                        engine.used.ncvar |= !engine.fresh.ncvar;
                        engine.fresh.ncvar = true;
                    }
                    ok
                } else {
                    indel::gap_prep(&mut site, bca, indel_buf)
                }
            }
        };
        if ok {
            for (g, r) in gplp.iter().zip(bcr.iter_mut()) {
                glfgen(g, None, bca, r);
            }
            if call::combine(bcr, bca, None, bc) >= 0 {
                call::call2bcf(
                    bc,
                    rec,
                    hdr,
                    conf.call.fmt_flag,
                    bca,
                    reference.unwrap_or(&[]),
                )?;
                sink.emit(hdr, rec)?;
                engine.used.bias |= !engine.fresh.bias;
            }
        }
    }
    engine.store = gplp.into_iter().map(recycle).collect();
    Ok(())
}
