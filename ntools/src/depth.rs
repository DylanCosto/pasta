//! Per-position read depth over coordinate-sorted alignments, as samtools depth counts it.

use std::{
    cmp::Reverse,
    collections::{BinaryHeap, VecDeque},
    io,
};

use pasta::sam::{
    self as sam,
    alignment::record::{Flags, cigar::op::Kind},
};

use crate::{alignment::Rec, bam_raw::Record};

// BAM CIGAR operation codes.
const CIGAR_MATCH: u32 = 0;
const CIGAR_INSERTION: u32 = 1;
const CIGAR_DELETION: u32 = 2;
const CIGAR_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_SEQUENCE_MATCH: u32 = 7;
const CIGAR_SEQUENCE_MISMATCH: u32 = 8;

/// samtools' default filter for depth and coverage.
pub const DEFAULT_EXCLUDE_FLAGS: Flags = Flags::UNMAPPED
    .union(Flags::SECONDARY)
    .union(Flags::QC_FAIL)
    .union(Flags::DUPLICATE);

/// Which bases count toward depth.
#[derive(Clone, Copy, Debug, Default)]
pub struct CountOptions {
    /// Bases below this quality are not counted (when positive).
    pub min_base_quality: u8,
    /// Count deletions (`-J`).
    pub include_deletions: bool,
    /// Sum base qualities into [`Cell::quality_sum`] (for coverage).
    pub track_quality: bool,
    /// Count reference skips (`N`), as samtools bedcov does by default.
    pub include_skips: bool,
}

/// Depth at one position, and the summed base qualities of the bases counted there.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Cell {
    pub depth: u32,
    /// Deletions and reference skips at this position, whether or not `depth` counts them.
    pub gaps: u32,
    pub quality_sum: u64,
    /// Counted bases without a quality (no sequence, or `0xff`), so not in `quality_sum`.
    pub missing_quality: u32,
    /// Whether any read aligns here, counted or not (e.g., a filtered base or a deletion).
    /// samtools reports touched positions even at depth 0.
    pub touched: bool,
}

/// Gaps (deletions and reference skips) at least this long are kept as ranges rather than
/// per-position cells, and the rest of the read is counted when the window reaches it, so a
/// long intron or deletion does not grow the window.
const LONG_GAP: usize = 1 << 16;

/// The reference span a read (or the part of one after a gap kept as a range) may have in
/// per-position cells before its next gap, however short, is kept as a range too. The cells of
/// a read are then about as many as its bases, whatever its CIGAR says of the reference.
const MAX_DENSE_SPAN: usize = 1 << 22;

/// The most positions a reader thread of `depth` or `coverage` sends per message.
const MAX_BATCH_LEN: usize = 1 << 16;

/// The fewest positions a reader thread sends per message, however many inputs there are.
const MIN_BATCH_LEN: usize = 1 << 10;

/// The positions per message for each of `inputs` reader threads. Each input can have its
/// channel's batches and two more in memory (24 bytes or more a position), so many inputs
/// share about four full batches rather than each holding up to 9 MiB: 1,000 inputs of 100 kb
/// of depths held 2.5 GB at `-@0`, and hold 0.4 GB with batches of 1,024 positions.
pub fn batch_len(inputs: usize) -> usize {
    (4 * MAX_BATCH_LEN / inputs.max(1)).clamp(MIN_BATCH_LEN, MAX_BATCH_LEN)
}

/// A gap kept as a range, which adds to the cells it covers.
#[derive(Clone, Copy, Debug)]
struct LongGap {
    /// 1-based, and half-open with `end`.
    start: usize,
    end: usize,
    /// Whether it counts in [`Cell::gaps`]: not where a counted deletion overlaps the mate.
    is_gap: bool,
    /// Whether it counts in [`Cell::depth`].
    is_counted: bool,
}

impl LongGap {
    fn contains(&self, position: usize) -> bool {
        (self.start..self.end).contains(&position)
    }
}

/// The part of a read after a gap kept as a range, counted once the window reaches its start.
#[derive(Debug)]
struct Segment {
    /// 1-based reference position of its first operation.
    start: usize,
    ops: Vec<u32>,
    quality_scores: Vec<u8>,
    needs_quality: bool,
    clip: usize,
}

/// A sliding window of depths for one reference sequence at a time.
///
/// Records must arrive sorted by reference sequence and start, so positions before a record's
/// start are final and can be emitted.
#[derive(Debug)]
pub struct DepthWindow {
    reference_sequence_id: Option<usize>,
    /// 1-based position of `cells[0]`.
    start: usize,
    cells: VecDeque<Cell>,
    /// The gaps kept as ranges.
    long_gaps: Vec<LongGap>,
    /// Read parts after those gaps, not yet counted.
    segments: Vec<Segment>,
    options: CountOptions,
    /// [`MAX_DENSE_SPAN`], or less in tests.
    max_dense_span: usize,
    /// Reused per record: the CIGAR, BAM-encoded (length << 4 | operation).
    ops: Vec<u32>,
    quality_scores: Vec<u8>,
    cap: Option<PileupCap>,
    /// HTSlib's pileup reference sequence and 0-based position, as far as adding reads goes:
    /// the last added read's.
    pileup: (usize, usize),
}

impl Default for DepthWindow {
    fn default() -> Self {
        Self {
            reference_sequence_id: None,
            start: 0,
            cells: VecDeque::new(),
            long_gaps: Vec::new(),
            segments: Vec::new(),
            options: CountOptions::default(),
            max_dense_span: MAX_DENSE_SPAN,
            ops: Vec::new(),
            quality_scores: Vec::new(),
            cap: None,
            pileup: (0, 0),
        }
    }
}

impl DepthWindow {
    pub fn new(options: CountOptions) -> Self {
        Self {
            options,
            ..Self::default()
        }
    }

    /// A window that leaves out reads as HTSlib's pileup does with a maximum count (`maxcnt`,
    /// samtools coverage's `-d`).
    pub fn with_max_count(options: CountOptions, max_count: usize) -> Self {
        Self {
            options,
            cap: Some(PileupCap::new(max_count)),
            ..Self::default()
        }
    }

    /// Adds a record's aligned bases, first emitting every position before its start. Returns
    /// whether HTSlib's pileup would buffer the record (and call its constructor): not if the
    /// maximum count leaves it out, or if it is empty and starts at the pileup's position.
    ///
    /// `emit` receives `(reference sequence ID, position, cell)` for every touched position.
    pub fn add<F>(
        &mut self,
        record: &Rec<'_>,
        header: &sam::Header,
        emit: &mut F,
    ) -> io::Result<bool>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        self.add_clipped(record, header, 0, emit)
    }

    /// [`Self::add`], leaving out the aligned bases (and counted deletions) before the 1-based
    /// position `clip`, where they overlap the read's mate (samtools depth `-s`).
    pub fn add_clipped<F>(
        &mut self,
        record: &Rec<'_>,
        header: &sam::Header,
        clip: usize,
        emit: &mut F,
    ) -> io::Result<bool>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        // BAM records are read raw; others through the record trait.
        let raw = match record {
            Rec::Bam(b) => Some(Record::new(b.as_ref())?),
            _ => None,
        };

        let Some((id, start)) = record_position(record, raw.as_ref(), header)? else {
            return Ok(true);
        };

        if self.reference_sequence_id != Some(id) {
            self.finish(emit)?;
            self.reference_sequence_id = Some(id);
            self.start = start;
        } else if start < self.start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "input is not coordinate-sorted",
            ));
        }

        self.flush_before(start, emit)?;

        let needs_quality = self.load_record(record, raw.as_ref())?;

        let reach = Reach::of(&self.ops, self.max_dense_span);

        if !self.ops.is_empty() {
            check_query_length(reach.query_len, record, raw.as_ref())?;
        }
        let span = reach.span;

        if let Some(cap) = &mut self.cap
            && !cap.admit(id, start - 1, start - 1 + span)
        {
            return Ok(false);
        }

        let buffered = span > 0 || id != self.pileup.0 || start - 1 > self.pileup.1;
        self.pileup = (id, start - 1);

        self.count_bases(start, reach, needs_quality, clip);

        Ok(buffered)
    }

    /// Reads a record's CIGAR, BAM-encoded, and its qualities if counting needs them; returns
    /// whether it does.
    fn load_record(&mut self, record: &Rec<'_>, raw: Option<&Record<'_>>) -> io::Result<bool> {
        let r = record.as_dyn();

        match raw {
            Some(raw) => raw.read_cigar(&mut self.ops)?,
            None => {
                self.ops.clear();

                for result in r.cigar().iter() {
                    let op = result?;
                    let len = u32::try_from(op.len())
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                    self.ops.push(len << 4 | kind_code(op.kind()));
                }
            }
        }

        let needs_quality = self.options.min_base_quality > 0 || self.options.track_quality;

        self.quality_scores.clear();

        if needs_quality {
            match (&raw, record) {
                (Some(raw), _) => self.quality_scores.extend_from_slice(raw.quality_scores()),
                (None, Rec::Buf(r)) => self
                    .quality_scores
                    .extend_from_slice(r.quality_scores().as_ref()),
                (None, _) => {
                    for result in r.quality_scores().iter() {
                        self.quality_scores.push(result?);
                    }
                }
            }
        }

        Ok(needs_quality)
    }

    /// Counts the bases of the record in [`Self::ops`] and [`Self::quality_scores`], which
    /// starts at the 1-based `start`.
    fn count_bases(&mut self, start: usize, reach: Reach, needs_quality: bool, clip: usize) {
        // The reusable buffers are moved out while cells are updated, then put back.
        let ops = std::mem::take(&mut self.ops);
        let quality_scores = std::mem::take(&mut self.quality_scores);

        self.count(start, &ops, reach, &quality_scores, needs_quality, clip);

        // HTSlib's pileup still visits a read that covers no reference bases (say `4I`) at its
        // position, so that position is reported, with no depth from it.
        if reach.span == 0
            && let Some(cell) = self.cells.get_mut(start - self.start)
        {
            cell.touched = true;
        }

        self.ops = ops;
        self.quality_scores = quality_scores;
    }

    /// Counts the bases of `ops` (BAM-encoded CIGAR operations), which start at the 1-based
    /// `start`, up to the first gap kept as a range; the rest of the read becomes a [`Segment`].
    fn count(
        &mut self,
        start: usize,
        ops: &[u32],
        reach: Reach,
        quality_scores: &[u8],
        needs_quality: bool,
        clip: usize,
    ) {
        let CountOptions {
            min_base_quality,
            include_deletions,
            include_skips,
            ..
        } = self.options;

        // Operations from the first gap kept as a range on are set aside after the others are
        // counted.
        let (ops, deferred) = ops.split_at(reach.long_gap);
        let span = reach.near_span;

        // Size the window once for the reference span up to that gap (at least one
        // position, as `bam_endpos` gives a read without reference bases).
        let needed = start + span.max(1) - self.start;
        if needed > self.cells.len() {
            self.cells.resize(needed, Cell::default());
        }

        // One slice for the whole record (rotating the ring only when it has wrapped), so the
        // loops below run over plain slices.
        let window = self.cells.make_contiguous();

        let mut reference_position = start;
        let mut read_position = 0;

        for &op in ops {
            let len = (op >> 4) as usize;
            let offset = reference_position - self.start;
            // Bases before `clip` overlap the mate and are not counted.
            let clipped = clip.saturating_sub(reference_position).min(len);

            match op & 0xf {
                CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => {
                    // Clipped positions are still reported, at the depth of the other reads.
                    for cell in &mut window[offset..offset + clipped] {
                        cell.touched = true;
                    }

                    let cells = &mut window[offset + clipped..offset + len];
                    let first_base = read_position + clipped;

                    if needs_quality {
                        for (i, cell) in cells.iter_mut().enumerate() {
                            cell.touched = true;

                            // Missing qualities (0xff) count, as in samtools.
                            let q = match quality_scores.get(first_base + i) {
                                Some(&q) if q != 0xff => q,
                                _ => {
                                    cell.depth += 1;
                                    cell.missing_quality += 1;
                                    continue;
                                }
                            };

                            if q < min_base_quality {
                                continue;
                            }

                            cell.depth += 1;
                            cell.quality_sum += u64::from(q);
                        }
                    } else {
                        for cell in cells {
                            cell.touched = true;
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                    read_position += len;
                }
                CIGAR_DELETION if include_deletions => {
                    for cell in &mut window[offset..offset + clipped] {
                        cell.touched = true;
                    }

                    // As in samtools, a deletion takes the quality of the base after it.
                    let passes = !needs_quality
                        || quality_scores
                            .get(read_position)
                            .is_none_or(|&q| q == 0xff || q >= min_base_quality);

                    for cell in &mut window[offset + clipped..offset + len] {
                        cell.touched = true;
                        cell.gaps += 1;

                        if passes {
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                }
                CIGAR_DELETION => {
                    for cell in &mut window[offset..offset + len] {
                        cell.touched = true;
                        cell.gaps += 1;
                    }

                    reference_position += len;
                }
                CIGAR_SKIP => {
                    // Reference skips are within the read's span, so they are reported too.
                    for cell in &mut window[offset..offset + len] {
                        cell.touched = true;
                        cell.gaps += 1;

                        if include_skips {
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                }
                CIGAR_INSERTION | CIGAR_SOFT_CLIP => read_position += len,
                _ => {}
            }
        }

        if let Some((&gap, rest)) = deferred.split_first() {
            self.defer(
                reference_position,
                gap,
                rest,
                quality_scores,
                read_position,
                needs_quality,
                clip,
            );
        }
    }

    /// Keeps the gap `gap` (a BAM-encoded deletion or reference skip) at `start` as a range, and
    /// the operations after it as a [`Segment`] (with the qualities from `read_position` on).
    #[expect(
        clippy::too_many_arguments,
        reason = "the state of a record being counted"
    )]
    #[cold]
    fn defer(
        &mut self,
        start: usize,
        gap: u32,
        rest: &[u32],
        quality_scores: &[u8],
        read_position: usize,
        needs_quality: bool,
        clip: usize,
    ) {
        let len = (gap >> 4) as usize;
        let end = start + len;

        if gap & 0xf == CIGAR_SKIP {
            self.long_gaps.push(LongGap {
                start,
                end,
                is_gap: true,
                is_counted: self.options.include_skips,
            });
        } else if self.options.include_deletions {
            // As in `count`: the deletion takes the quality of the base after it, and where it
            // overlaps the mate it is reported but not counted.
            let passes = !needs_quality
                || quality_scores
                    .get(read_position)
                    .is_none_or(|&q| q == 0xff || q >= self.options.min_base_quality);

            let clipped = clip.saturating_sub(start).min(len);

            let parts = [
                (start, start + clipped, false, false),
                (start + clipped, end, true, passes),
            ];

            for (start, end, is_gap, is_counted) in parts {
                if start < end {
                    self.long_gaps.push(LongGap {
                        start,
                        end,
                        is_gap,
                        is_counted,
                    });
                }
            }
        } else {
            self.long_gaps.push(LongGap {
                start,
                end,
                is_gap: true,
                is_counted: false,
            });
        }

        if !rest.is_empty() {
            self.segments.push(Segment {
                start: end,
                ops: rest.to_vec(),
                quality_scores: quality_scores
                    .get(read_position..)
                    .unwrap_or_default()
                    .to_vec(),
                needs_quality,
                clip,
            });
        }
    }

    /// Counts the read parts that start at the window's first position.
    fn count_segments(&mut self) {
        while let Some(i) = self
            .segments
            .iter()
            .position(|segment| segment.start == self.start)
        {
            let segment = self.segments.swap_remove(i);
            self.count(
                segment.start,
                &segment.ops,
                Reach::of(&segment.ops, self.max_dense_span),
                &segment.quality_scores,
                segment.needs_quality,
                segment.clip,
            );
        }
    }

    /// Adds the gaps kept as ranges that cover `position` to its cell.
    fn add_long_gaps(&self, position: usize, cell: &mut Cell) {
        for gap in &self.long_gaps {
            if gap.contains(position) {
                cell.touched = true;
                cell.gaps += u32::from(gap.is_gap);
                cell.depth += u32::from(gap.is_counted);
            }
        }
    }

    /// The next position after the window's first where a gap kept as a range or a read part
    /// starts.
    fn next_event(&self) -> Option<usize> {
        let gaps = self.long_gaps.iter().map(|gap| gap.start);
        let segments = self.segments.iter().map(|segment| segment.start);
        gaps.chain(segments).filter(|&p| p > self.start).min()
    }

    fn flush_before<F>(&mut self, position: usize, emit: &mut F) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        let Some(id) = self.reference_sequence_id else {
            return Ok(());
        };

        // Without gaps kept as ranges (almost always), cells are emitted directly.
        if self.long_gaps.is_empty() && self.segments.is_empty() {
            while self.start < position {
                match self.cells.pop_front() {
                    Some(cell) => {
                        if cell.touched {
                            emit(id, self.start, cell)?;
                        }
                        self.start += 1;
                    }
                    None => {
                        self.start = position;
                        break;
                    }
                }
            }

            return Ok(());
        }

        self.flush_with_long_gaps(id, position, emit)
    }

    /// [`Self::flush_before`] while gaps kept as ranges or read parts after them are pending:
    /// kept apart so that the usual loop stays small enough to inline `emit`.
    #[cold]
    #[inline(never)]
    fn flush_with_long_gaps<F>(
        &mut self,
        id: usize,
        position: usize,
        emit: &mut F,
    ) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        while self.start < position {
            if !self.segments.is_empty() {
                self.count_segments();
            }

            let cell = match self.cells.pop_front() {
                Some(cell) => Some(cell),
                // Only gaps kept as ranges cover this position, if anything.
                None if self.long_gaps.iter().any(|gap| gap.contains(self.start)) => {
                    Some(Cell::default())
                }
                None => None,
            };

            let Some(mut cell) = cell else {
                // Nothing until the next gap or read part (or ever).
                self.start = self.next_event().map_or(position, |p| p.min(position));
                continue;
            };

            if !self.long_gaps.is_empty() {
                self.add_long_gaps(self.start, &mut cell);
            }

            if cell.touched {
                emit(id, self.start, cell)?;
            }

            self.start += 1;

            if !self.long_gaps.is_empty() {
                let start = self.start;
                self.long_gaps.retain(|gap| gap.end > start);
            }
        }

        Ok(())
    }

    /// Emits the remaining positions of the current reference sequence.
    pub fn finish<F>(&mut self, emit: &mut F) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        self.flush_before(usize::MAX, emit)?;
        self.reference_sequence_id = None;
        Ok(())
    }
}

/// HTSlib's pileup depth limit (`maxcnt`) for one input: a read that starts where the previous
/// one did is left out once `maxcnt` reads are buffered. The buffer holds the reads that have
/// not ended before the position just before the new read's start.
#[derive(Clone, Debug)]
struct PileupCap {
    max_count: usize,
    /// The pileup's reference sequence and 0-based position: the last added read's start.
    id: usize,
    position: usize,
    /// The 0-based exclusive ends of the buffered reads.
    ends: BinaryHeap<Reverse<usize>>,
}

impl PileupCap {
    fn new(max_count: usize) -> Self {
        Self {
            max_count,
            id: 0,
            position: 0,
            ends: BinaryHeap::new(),
        }
    }

    /// Returns whether a read at the 0-based `start` with the exclusive `end` enters the pileup
    /// (`bam_plp_push`).
    fn admit(&mut self, id: usize, start: usize, end: usize) -> bool {
        let new_reference = id != self.id;

        if new_reference {
            self.ends.clear();
        } else {
            // The pileup has dropped the reads that end before `start - 1`.
            while self.ends.peek().is_some_and(|&Reverse(e)| e < start) {
                self.ends.pop();
            }

            if start == self.position && self.ends.len() >= self.max_count {
                return false;
            }
        }

        // An empty read starting at the pileup's position is not buffered.
        if new_reference || end > self.position {
            self.ends.push(Reverse(end));
        }

        self.id = id;
        self.position = start;

        true
    }
}

/// A record's reference sequence ID and 1-based start, or `None` if it is unplaced; `raw` is
/// the record's raw view if it is BAM.
fn record_position(
    record: &Rec<'_>,
    raw: Option<&Record<'_>>,
    header: &sam::Header,
) -> io::Result<Option<(usize, usize)>> {
    let r = record.as_dyn();

    let position = match raw {
        Some(raw) => {
            let (Ok(id), Ok(position)) = (
                usize::try_from(raw.reference_sequence_id()),
                usize::try_from(raw.position()),
            ) else {
                return Ok(None);
            };

            (id, position + 1)
        }
        None => {
            let (Some(id), Some(start)) = (
                r.reference_sequence_id(header).transpose()?,
                r.alignment_start().transpose()?,
            ) else {
                return Ok(None);
            };

            (id, usize::from(start))
        }
    };

    Ok(Some(position))
}

/// samtools depth's `qlen_used`: the sequence length without leading and trailing soft clips,
/// or, without a sequence, the query length of the CIGAR.
pub fn used_query_length(record: &Rec<'_>) -> io::Result<usize> {
    let r = record.as_dyn();
    let ops = r.cigar().iter().collect::<io::Result<Vec<_>>>()?;
    let sequence_len = r.sequence().len();

    if sequence_len == 0 {
        return Ok(query_length(&ops));
    }

    let leading = ops
        .iter()
        .take_while(|op| op.kind() == Kind::SoftClip)
        .count();
    let soft_clipped: usize = ops[..leading]
        .iter()
        .chain(
            ops[leading..]
                .iter()
                .rev()
                .take_while(|op| op.kind() == Kind::SoftClip),
        )
        .map(|op| op.len())
        .sum();

    Ok(sequence_len.saturating_sub(soft_clipped))
}

/// HTSlib's `bam_cigar2qlen`: the query length of a record's CIGAR, soft clips included.
pub fn cigar_query_length(record: &Rec<'_>) -> io::Result<usize> {
    let ops = record
        .as_dyn()
        .cigar()
        .iter()
        .collect::<io::Result<Vec<_>>>()?;
    Ok(query_length(&ops))
}

fn query_length(ops: &[pasta::sam::alignment::record::cigar::Op]) -> usize {
    ops.iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len())
        .sum()
}

/// The BAM code of a CIGAR operation.
fn kind_code(kind: Kind) -> u32 {
    match kind {
        Kind::Match => CIGAR_MATCH,
        Kind::Insertion => CIGAR_INSERTION,
        Kind::Deletion => CIGAR_DELETION,
        Kind::Skip => CIGAR_SKIP,
        Kind::SoftClip => CIGAR_SOFT_CLIP,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => CIGAR_SEQUENCE_MATCH,
        Kind::SequenceMismatch => CIGAR_SEQUENCE_MISMATCH,
    }
}

/// How far a read's operations reach on the reference.
#[derive(Clone, Copy, Debug)]
struct Reach {
    /// The reference span of all operations.
    span: usize,
    /// The index of the first gap kept as a range, or the number of operations.
    long_gap: usize,
    /// The reference span of the operations before `long_gap`.
    near_span: usize,
    /// The bases of the read that all operations take (HTSlib's `bam_cigar2qlen`).
    query_len: u64,
}

impl Reach {
    /// A gap is kept as a range if it is long, or if the span would pass `max_dense_span` with
    /// it.
    fn of(ops: &[u32], max_dense_span: usize) -> Self {
        let mut span = 0;
        let mut long_gap = None;
        let mut query_len = 0;

        for (i, &op) in ops.iter().enumerate() {
            let (code, len) = (op & 0xf, (op >> 4) as usize);

            if consumes_query(code) {
                query_len += u64::from(op >> 4);
            }

            if consumes_reference(code) {
                if matches!(code, CIGAR_DELETION | CIGAR_SKIP)
                    && long_gap.is_none()
                    && (len >= LONG_GAP || span + len > max_dense_span)
                {
                    long_gap = Some((i, span));
                }

                span += len;
            }
        }

        let (long_gap, near_span) = long_gap.unwrap_or((ops.len(), span));

        Self {
            span,
            long_gap,
            near_span,
            query_len,
        }
    }
}

/// Fails if a record's CIGAR, which takes `query_len` bases, and its sequence differ in length,
/// as HTSlib's `bam_read1` does for a mapped record that has both. `raw` is the record's raw
/// view if it is BAM.
fn check_query_length(
    query_len: u64,
    record: &Rec<'_>,
    raw: Option<&Record<'_>>,
) -> io::Result<()> {
    let r = record.as_dyn();

    let sequence_len = match raw {
        Some(raw) => raw.sequence_len(),
        None => r.sequence().len(),
    };

    // A read without its sequence has none to compare.
    if query_len == sequence_len as u64 || sequence_len == 0 {
        return Ok(());
    }

    let is_unmapped = match raw {
        Some(raw) => Flags::from(raw.flags()).is_unmapped(),
        None => r.flags()?.is_unmapped(),
    };

    if is_unmapped {
        return Ok(());
    }

    let name = match raw {
        Some(raw) => String::from_utf8_lossy(raw.name()).into_owned(),
        None => r.name().map(|name| name.to_string()).unwrap_or_default(),
    };

    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("CIGAR and query sequence lengths differ for {name}"),
    ))
}

/// Whether a BAM CIGAR operation code consumes the read.
fn consumes_query(code: u32) -> bool {
    matches!(
        code,
        CIGAR_MATCH
            | CIGAR_INSERTION
            | CIGAR_SOFT_CLIP
            | CIGAR_SEQUENCE_MATCH
            | CIGAR_SEQUENCE_MISMATCH
    )
}

/// Whether a BAM CIGAR operation code consumes the reference.
fn consumes_reference(code: u32) -> bool {
    matches!(
        code,
        CIGAR_MATCH | CIGAR_DELETION | CIGAR_SKIP | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
    )
}

/// Formats like C's `printf("%.*g", precision, value)`.
pub fn format_g(value: f64, precision: usize) -> String {
    let precision = precision.max(1);

    if value == 0.0 {
        return "0".to_owned();
    }

    if !value.is_finite() {
        return if value.is_nan() {
            "nan".to_owned()
        } else if value > 0.0 {
            "inf".to_owned()
        } else {
            "-inf".to_owned()
        };
    }

    // The exponent after rounding to `precision` significant digits decides the style.
    let scientific = format!("{:.*e}", precision - 1, value);
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("scientific notation has an exponent");
    let exponent: i32 = exponent.parse().expect("valid exponent");

    let precision_i32 = i32::try_from(precision).unwrap_or(i32::MAX);

    if exponent < -4 || exponent >= precision_i32 {
        let mantissa = trim_fraction(mantissa);
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    } else {
        let decimals = usize::try_from(precision_i32 - 1 - exponent).unwrap_or(0);
        trim_fraction(&format!("{value:.decimals$}")).to_owned()
    }
}

fn trim_fraction(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use pasta::{
        core::Position,
        sam::{
            alignment::{RecordBuf, record::cigar::Op},
            header::record::value::{Map, map::ReferenceSequence},
        },
    };

    use super::*;

    fn header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "c",
                Map::<ReferenceSequence>::new(NonZero::new(1 << 28).unwrap()),
            )
            .build()
    }

    /// A read at `start` whose bases, as many as `cigar` says, have quality `quality`.
    fn record(start: usize, cigar: &[(Kind, usize)], quality: u8) -> RecordBuf {
        let ops: Vec<_> = cigar
            .iter()
            .map(|&(kind, len)| Op::new(kind, len))
            .collect();
        let len = query_length(&ops);

        RecordBuf::builder()
            .set_name("r")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::new(start).unwrap())
            .set_cigar(ops.into_iter().collect())
            .set_sequence(vec![b'A'; len].into())
            .set_quality_scores(vec![quality; len].into())
            .build()
    }

    /// The cells a window emits for `records`, each with the position its mate overlap ends at.
    fn cells(
        options: CountOptions,
        max_dense_span: usize,
        records: &[(RecordBuf, usize)],
    ) -> io::Result<Vec<(usize, Cell)>> {
        let header = header();
        let mut window = DepthWindow::new(options);
        window.max_dense_span = max_dense_span;

        let mut cells = Vec::new();
        let mut emit = |_, position, cell| {
            cells.push((position, cell));
            Ok(())
        };

        for (record, clip) in records {
            window.add_clipped(&Rec::Buf(record), &header, *clip, &mut emit)?;
        }

        window.finish(&mut emit)?;

        Ok(cells)
    }

    /// Every way the options count a deletion.
    fn all_options() -> Vec<CountOptions> {
        let mut options = Vec::new();

        for include_deletions in [false, true] {
            for min_base_quality in [0, 20, 40] {
                for track_quality in [false, true] {
                    options.push(CountOptions {
                        min_base_quality,
                        include_deletions,
                        track_quality,
                        include_skips: include_deletions,
                    });
                }
            }
        }

        options
    }

    #[test]
    fn test_add_counts_a_long_deletion_as_shorter_ones() -> io::Result<()> {
        use Kind::{Deletion, Insertion, Match, Skip, SoftClip};

        // The same alignments with a deletion of 70,000 whole, kept as a range, and in two
        // halves, counted in cells; with reads before, within, across and after it, one of
        // quality below the thresholds, and a mate overlap that ends within it.
        let reads = |deletion: &[(Kind, usize)]| {
            let mut long = vec![(SoftClip, 3), (Match, 10)];
            long.extend_from_slice(deletion);
            long.extend([
                (Match, 5),
                (Insertion, 2),
                (Match, 5),
                (Skip, 4),
                (Match, 3),
            ]);

            let mut poor = vec![(Match, 8)];
            poor.extend_from_slice(deletion);
            poor.push((Match, 8));

            [
                (record(90, &[(Match, 20)], 30), 0),
                (record(100, &long, 30), 0),
                (record(102, &poor, 10), 0),
                (record(105, &long, 30), 40_000),
                (record(30_000, &[(Match, 50)], 30), 0),
                (record(70_100, &[(Match, 30)], 30), 0),
                (record(140_500, &[(Match, 30)], 30), 0),
            ]
        };

        let whole = reads(&[(Deletion, 70_000)]);
        let halves = reads(&[(Deletion, 35_000), (Deletion, 35_000)]);

        for options in all_options() {
            let actual = cells(options, MAX_DENSE_SPAN, &whole)?;
            let expected = cells(options, MAX_DENSE_SPAN, &halves)?;

            // Positions 90 to 70,131, and 140,500 to 140,529.
            assert_eq!(actual.len(), 70_042 + 30, "{options:?}");
            assert!(actual == expected, "{options:?}");
        }

        Ok(())
    }

    #[test]
    fn test_add_with_a_dense_span_limit() -> io::Result<()> {
        use Kind::{Deletion, Insertion, Match, Skip};

        // Reads of many gaps: with a limit below their spans, gaps become ranges wherever the
        // span would pass it, and the cells are those of a window without a limit.
        let gapped: Vec<_> = (0..40)
            .flat_map(|i| {
                [
                    (Match, 3 + i % 4),
                    (if i % 3 == 0 { Skip } else { Deletion }, 5 + i % 7),
                    (Insertion, i % 2),
                ]
            })
            .filter(|&(_, len)| len > 0)
            .chain([(Match, 6)])
            .collect();

        let reads = [
            (record(10, &[(Match, 40)], 30), 0),
            (record(20, &gapped, 30), 0),
            (record(25, &gapped, 10), 0),
            (record(31, &gapped, 30), 200),
            (record(300, &[(Match, 100)], 30), 0),
            (record(700, &gapped, 30), 0),
        ];

        for options in all_options() {
            let expected = cells(options, usize::MAX, &reads)?;

            for max_dense_span in [1, 7, 50, 300] {
                let actual = cells(options, max_dense_span, &reads)?;
                assert!(actual == expected, "{options:?}, {max_dense_span}");
            }
        }

        Ok(())
    }

    #[test]
    fn test_add_rejects_cigar_and_sequence_lengths_that_differ() -> io::Result<()> {
        let options = CountOptions::default();

        // A CIGAR of 300 bases and a sequence of 4, as a corrupt record has.
        let mut corrupt = record(10, &[(Kind::Match, 4)], 30);
        *corrupt.cigar_mut() = [Op::new(Kind::Match, 300)].into_iter().collect();

        let Err(error) = cells(options, MAX_DENSE_SPAN, &[(corrupt.clone(), 0)]) else {
            panic!("the record was counted");
        };
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            error.to_string(),
            "CIGAR and query sequence lengths differ for r"
        );

        // A CIGAR that takes no bases differs from a sequence too.
        let mut deletion = record(10, &[(Kind::Match, 4)], 30);
        *deletion.cigar_mut() = [Op::new(Kind::Deletion, 5)].into_iter().collect();
        assert!(cells(options, MAX_DENSE_SPAN, &[(deletion, 0)]).is_err());

        // A read without a CIGAR has none to compare.
        let mut bare = record(10, &[(Kind::Match, 4)], 30);
        *bare.cigar_mut() = Default::default();
        assert!(cells(options, MAX_DENSE_SPAN, &[(bare, 0)]).is_ok());

        // HTSlib does not compare the lengths of an unmapped read.
        *corrupt.flags_mut() = Flags::UNMAPPED;
        assert_eq!(cells(options, MAX_DENSE_SPAN, &[(corrupt, 0)])?.len(), 300);

        // A read without its sequence has none to compare.
        let mut bare = record(10, &[(Kind::Match, 4)], 30);
        *bare.sequence_mut() = Vec::new().into();
        *bare.quality_scores_mut() = Vec::new().into();

        assert_eq!(cells(options, MAX_DENSE_SPAN, &[(bare, 0)])?.len(), 4);

        Ok(())
    }

    #[test]
    fn test_pileup_cap() {
        let mut cap = PileupCap::new(2);

        assert!(cap.admit(0, 10, 20));
        assert!(cap.admit(0, 10, 20));
        // Two reads are buffered at the pileup's position.
        assert!(!cap.admit(0, 10, 20));
        // The first read at a new position is always added.
        assert!(cap.admit(0, 11, 12));
        assert!(!cap.admit(0, 11, 30));
        // The reads before have ended; an empty read is buffered only when it comes first.
        assert!(cap.admit(0, 40, 40));
        assert!(cap.admit(0, 40, 40));
        assert!(cap.admit(0, 40, 50));
        assert!(!cap.admit(0, 40, 50));
        assert!(cap.admit(1, 40, 50));
    }

    #[test]
    fn test_format_g() {
        let cases = [
            (50.001_775, 6, "50.0018"),
            (17.718_2, 6, "17.7182"),
            (25.0, 3, "25"),
            (35.2, 3, "35.2"),
            (0.652_32, 6, "0.65232"),
            (100.0, 6, "100"),
            (0.000_012_34, 6, "1.234e-05"),
            (1_234_567.0, 6, "1.23457e+06"),
            (99.999_99, 3, "100"),
            (0.0, 6, "0"),
        ];

        for (value, precision, expected) in cases {
            assert_eq!(
                format_g(value, precision),
                expected,
                "{value} with {precision}"
            );
        }
    }
}
