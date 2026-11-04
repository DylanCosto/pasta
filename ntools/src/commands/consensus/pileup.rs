//! The consensus pileup engine, a port of samtools' `consensus_pileup.c`.
//!
//! This is not htslib's `bam_plp`: it visits one column per reference position and one more per
//! inserted base (`nth` > 0), and keeps per-read state such as the running quality of a read
//! that is padded across another read's insertion.
//!
//! Unlike samtools, a read's state does not depend on earlier reads: a recycled read struct
//! starts with the quality of a new one (samtools keeps the last read's, which seeds a read whose
//! first column is a deletion), and the quality of a padded or deleted position comes from the
//! read's own qualities (samtools reads a byte past them, from an earlier record).

use anyhow::bail;

use crate::bam_raw::{self, le_u16, le_u32};

/// BAM CIGAR operation codes.
pub const CMATCH: i32 = 0;
pub const CINS: i32 = 1;
pub const CDEL: i32 = 2;
pub const CREF_SKIP: i32 = 3;
pub const CSOFT_CLIP: i32 = 4;
pub const CHARD_CLIP: i32 = 5;
pub const CPAD: i32 = 6;
pub const CEQUAL: i32 = 7;
pub const CDIFF: i32 = 8;

const FUNMAP: u16 = 0x4;
const FREVERSE: u16 = 0x10;

/// `HTS_POS_MAX`
pub const POS_MAX: i64 = (1 << 62) - 1;

/// `seq_nt16_str`
const NT16: &[u8; 16] = b"NACMGRSVTWYHKDBN";

/// One read in the pileup (`pileup_t`), with its BAM record.
#[derive(Debug, Default)]
pub struct Entry {
    /// The record, without `block_size`, in `buf[..len]`; see the module docs for the rest.
    buf: Vec<u8>,
    len: usize,

    pub tid: i32,
    /// The 0-based alignment start.
    pub core_pos: i64,
    pub flag: u16,
    pub mapq: u8,
    pub l_qseq: i32,
    n_cigar: u32,
    cigar_off: usize,
    /// The CIGAR of a record that stores it in a `CG` field (more than 65,535 operations), as
    /// HTSlib reads it; empty otherwise.
    long_cigar: Vec<u32>,
    seq_off: usize,
    pub qual_off: usize,

    /// 1 if the read has finished, 2 or 3 at a reference skip (`N`).
    pub eof: i32,
    /// The quality of the current base.
    pub qual: i32,
    start: i32,
    /// The current base in ASCII (`*` for a deletion or pad).
    pub base: u8,
    /// Whether the current column is at or next to a reference skip.
    pub ref_skip: bool,
    padding: bool,
    /// The current base as a 4-bit code, or 16 for a deletion or pad.
    pub base4: i32,
    pos: i64,
    nth: i32,
    pub is_rev: bool,
    /// The index of the current base in the read.
    pub seq_offset: i32,
    cigar_ind: u32,
    cigar_op: i32,
    cigar_len: i32,
    first_del: bool,

    /// Client data (`p->cd`): the consensus caller's local NM array, if `has_nm`. The buffer
    /// is kept when the entry is recycled.
    nm: Vec<i32>,
    has_nm: bool,
}

impl Entry {
    /// Copies `raw` into this entry's buffer and parses its fixed fields.
    fn load(&mut self, raw: &[u8]) -> anyhow::Result<()> {
        if raw.len() < 32 {
            bail!("truncated BAM record");
        }

        if self.buf.len() < raw.len() {
            self.buf.resize(raw.len(), 0);
        }

        self.buf[..raw.len()].copy_from_slice(raw);
        self.len = raw.len();

        self.tid = le_u32(&raw[0..4]) as i32;
        self.core_pos = i64::from(le_u32(&raw[4..8]) as i32);
        self.mapq = raw[9];
        self.n_cigar = u32::from(le_u16(&raw[12..14]));
        self.flag = le_u16(&raw[14..16]);
        let l_qseq = le_u32(&raw[16..20]);

        let name_len = usize::from(raw[8]);
        self.cigar_off = 32 + name_len;
        self.seq_off = self.cigar_off + 4 * self.n_cigar as usize;
        self.qual_off = self.seq_off + (l_qseq as usize).div_ceil(2);
        self.l_qseq = i32::try_from(l_qseq)?;

        if self.qual_off + l_qseq as usize > raw.len() {
            bail!("truncated BAM record");
        }

        self.long_cigar.clear();

        // A `kSmN` placeholder stands for the CIGAR in a `CG` field.
        if self.n_cigar == 2
            && let Some(range) = bam_raw::Record::new(raw)?.long_cigar_field()?
        {
            self.long_cigar.extend(
                raw[range.start + 8..range.end]
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .map(|b| u32::from_le_bytes(*b)),
            );
            self.n_cigar = u32::try_from(self.long_cigar.len())?;
        }

        Ok(())
    }

    /// The local NM array set by [`Self::set_local_nm`], if any.
    pub fn local_nm(&self) -> Option<&[i32]> {
        self.has_nm.then_some(self.nm.as_slice())
    }

    /// Takes the local NM buffer (to be filled and given back with [`Self::set_local_nm`]),
    /// zeroed to `len` entries.
    pub fn take_nm_buffer(&mut self, len: usize) -> Vec<i32> {
        let mut nm = std::mem::take(&mut self.nm);
        self.has_nm = false;
        nm.clear();
        nm.resize(len, 0);
        nm
    }

    /// Sets the local NM array (`p->cd`).
    pub fn set_local_nm(&mut self, nm: Vec<i32>) {
        self.nm = nm;
        self.has_nm = true;
    }

    /// The record's bytes.
    pub fn record(&self) -> &[u8] {
        &self.buf[..self.len]
    }

    /// The `i`th raw CIGAR operation (from the `CG` field for a long CIGAR).
    pub fn cigar(&self, i: u32) -> u32 {
        if self.long_cigar.is_empty() {
            let off = self.cigar_off + 4 * i as usize;
            le_u32(&self.buf[off..off + 4])
        } else {
            self.long_cigar[i as usize]
        }
    }

    /// The number of CIGAR operations.
    pub fn n_cigar(&self) -> u32 {
        self.n_cigar
    }

    /// `bam_seqi`: the 4-bit code of base `i`.
    pub fn seqi(&self, i: i32) -> i32 {
        let byte = self.buf[self.seq_off + (i >> 1) as usize];
        i32::from((byte >> ((!i & 1) << 2)) & 0xf)
    }

    /// The quality of base `i`, clamped to the read's bases (samtools reads just outside its
    /// qualities there); 0 for a read without bases.
    pub fn qual_at(&self, i: i32) -> i32 {
        if self.l_qseq <= 0 {
            return 0;
        }

        let i = i.clamp(0, self.l_qseq - 1) as usize;
        i32::from(self.buf[self.qual_off + i])
    }

    /// The qualities, for in-place editing.
    pub fn quals_mut(&mut self) -> &mut [u8] {
        let start = self.qual_off;
        let end = start + self.l_qseq as usize;
        &mut self.buf[start..end]
    }

    /// The qualities.
    pub fn quals(&self) -> &[u8] {
        &self.buf[self.qual_off..self.qual_off + self.l_qseq as usize]
    }

    fn load_op(&mut self) -> i32 {
        let c = self.cigar(self.cigar_ind);
        self.cigar_op = (c & 0xf) as i32;
        self.cigar_len = (c >> 4) as i32;
        self.cigar_ind += 1;
        self.cigar_op
    }

    /// `get_next_base`: moves to the `nth` base at reference position `pos` and fills out
    /// `base`, `base4`, and `qual`. Returns 1 if a base was fetched, 0 if the read ran out, and
    /// -1 on an unknown CIGAR operation; `is_insert` is set to the number of inserted bases
    /// still to come in this column.
    fn next_base(&mut self, pos: i64, nth: i32, is_insert: &mut i32) -> i32 {
        let mut op = self.cigar_op;

        self.start -= i32::from(self.start > 0);

        if self.first_del && op != CPAD {
            self.first_del = false;
        }

        *is_insert = 0;

        // Find pos first.
        while self.pos < pos {
            self.nth = 0;

            if self.cigar_len == 0 {
                if self.cigar_ind >= self.n_cigar {
                    self.eof = 1;
                    return 0;
                }

                op = self.load_op();
            }

            if matches!(op, CMATCH | CEQUAL | CDIFF) && i64::from(self.cigar_len) <= pos - self.pos
            {
                self.seq_offset += self.cigar_len;
                self.pos += i64::from(self.cigar_len);
                self.cigar_len = 0;
            } else {
                match op {
                    CMATCH | CEQUAL | CDIFF => {
                        self.seq_offset += 1;
                        self.pos += 1;
                        self.cigar_len -= 1;
                    }
                    CDEL | CREF_SKIP => {
                        self.pos += 1;
                        self.cigar_len -= 1;
                    }
                    CINS | CSOFT_CLIP => {
                        self.seq_offset += self.cigar_len;
                        self.cigar_len = 0;
                    }
                    CPAD | CHARD_CLIP => self.cigar_len = 0,
                    _ => {
                        eprintln!("Unhandled cigar_op {op}");
                        return -1;
                    }
                }
            }
        }

        // Now at pos, find the nth base.
        while self.nth < nth {
            if self.cigar_len == 0 {
                if self.cigar_ind >= self.n_cigar {
                    self.eof = 1;
                    return 0;
                }

                op = self.load_op();
            }

            match op {
                CMATCH | CEQUAL | CDIFF | CSOFT_CLIP | CDEL | CREF_SKIP => break,
                CINS => {
                    self.seq_offset += 1;
                    self.cigar_len -= 1;
                    self.nth += 1;
                }
                CPAD => {
                    self.cigar_len -= 1;
                    self.nth += 1;
                }
                CHARD_CLIP => self.cigar_len = 0,
                _ => {
                    eprintln!("Unhandled cigar_op {op}");
                    return -1;
                }
            }
        }

        // Fill out the base and quality.
        self.ref_skip = false;

        if self.nth < nth && op != CINS {
            self.base = b'*';
            self.base4 = 16;
            self.padding = true;
            self.qual = if self.seq_offset < self.l_qseq {
                self.qual.min(self.qual_at(self.seq_offset + 1))
            } else {
                0
            };
        } else {
            self.padding = false;

            match op {
                CDEL | CPAD => {
                    self.base = b'*';
                    self.base4 = 16;
                    let q = if self.seq_offset + 1 < self.l_qseq {
                        self.qual_at(self.seq_offset + 1)
                    } else {
                        self.qual_at(self.seq_offset)
                    };
                    self.qual = self.qual.min(q);
                }
                CREF_SKIP => {
                    self.base = b'.';
                    self.base4 = 0;
                    self.qual = 0;
                    // The end of a fragment, but not of the sequence.
                    self.eof = if self.eof != 0 { 2 } else { 3 };
                    self.ref_skip = true;
                }
                _ => {
                    if (0..self.l_qseq).contains(&self.seq_offset) {
                        self.qual = i32::from(self.buf[self.qual_off + self.seq_offset as usize]);
                        self.base4 = self.seqi(self.seq_offset);
                        self.base = NT16[self.base4 as usize];
                    } else {
                        self.base = b'N';
                        self.base4 = 15;
                        self.qual = 0xff;
                    }
                }
            }
        }

        // Moving out of a reference skip into sequence again.
        if self.eof != 0 && self.base != b'.' {
            self.start = 1;
            self.ref_skip = true;
            self.eof = 0;
        }

        if self.start != 0 && self.cigar_op == CDEL {
            self.first_del = true;
        }

        // Check whether the next operation is an insertion of some sort.
        if self.cigar_len == 0 {
            if self.cigar_ind < self.n_cigar {
                op = self.load_op();

                if op == CREF_SKIP {
                    self.eof = 3;
                    self.ref_skip = true;
                }
            } else {
                self.eof = 1;
            }
        }

        match op {
            CPAD | CINS => *is_insert = self.cigar_len,
            CSOFT_CLIP => {
                // A last S operation is the end.
                let last = self.cigar_ind == self.n_cigar
                    || (self.cigar_ind + 1 == self.n_cigar
                        && (self.cigar(self.cigar_ind) & 0xf) as i32 == CHARD_CLIP);
                self.eof = i32::from(last);
            }
            CHARD_CLIP => self.eof = 1,
            _ => {}
        }

        1
    }
}

/// One pileup column, as passed to `seq_column`.
pub struct Column<'a> {
    entries: &'a [Entry],
    active: &'a [usize],
    /// The reference sequence ID of the first read.
    pub tid: i32,
    /// The 1-based reference position.
    pub pos: i64,
    /// 0 at the reference base, then 1, 2, ... for inserted bases.
    pub nth: i32,
    /// The number of reads in the column.
    pub depth: i32,
}

impl<'a> Column<'a> {
    /// The reads in pileup order.
    pub fn reads(&self) -> impl Iterator<Item = &'a Entry> + use<'a> {
        let entries = self.entries;
        self.active.iter().map(move |&i| &entries[i])
    }
}

/// The callbacks of `pileup_loop`.
pub trait Handler {
    /// `seq_init`: returns whether to keep a read entering the pileup.
    fn init(&mut self, entry: &mut Entry) -> anyhow::Result<bool>;

    /// `seq_column`: called for every column.
    fn column(&mut self, column: &Column<'_>) -> anyhow::Result<()>;
}

/// Filters applied as records are fetched (`readaln2`).
#[derive(Clone, Copy, Debug, Default)]
pub struct Filters {
    pub incl_flags: u16,
    pub excl_flags: u16,
    pub min_mqual: i32,
}

/// The state of `pileup_loop` between fetched records.
#[derive(Debug)]
pub struct Pileup {
    entries: Vec<Entry>,
    /// The active list, in order (`phead`).
    active: Vec<usize>,
    /// Recycled entries; the last is the head of `pfree`.
    free: Vec<usize>,
    /// The entry the next record is read into.
    pnew: usize,
    col: i64,
    nth: i32,
    last_ref: i32,
    filters: Filters,
    removed: Vec<usize>,
    /// The first and last (1-based) columns the handler looks at; see [`Self::set_window`].
    first_col: i64,
    last_col: i64,
}

impl Pileup {
    pub fn new(filters: Filters) -> Self {
        Self {
            entries: vec![Entry::default()],
            active: Vec::new(),
            free: Vec::new(),
            pnew: 0,
            col: 0,
            nth: 0,
            last_ref: -1,
            filters,
            removed: Vec::new(),
            first_col: 0,
            last_col: i64::MAX,
        }
    }

    /// Limits the traversal to the columns at 1-based positions `first..=last` (with the
    /// insertion columns after `last`), for a handler that ignores every other column.
    ///
    /// A read that crosses many threaded blocks (a long reference skip, say) is otherwise walked
    /// from its start to its end by every block. Here the columns after `last` are not visited,
    /// and runs of columns before `first` in which no read changes CIGAR operation are stepped
    /// over at once, leaving each read as visiting them would have. The handler may still see
    /// some columns before `first`.
    pub fn set_window(&mut self, first: i64, last: i64) {
        self.first_col = first;
        self.last_col = last;
    }

    /// Adds the next record (raw BAM, without `block_size`) in file order.
    pub fn push<H: Handler>(&mut self, raw: &[u8], handler: &mut H) -> anyhow::Result<()> {
        let pnew = self.pnew;
        self.entries[pnew].load(raw)?;

        let (flag, mapq, tid, core_pos) = {
            let b = &self.entries[pnew];
            (b.flag, b.mapq, b.tid, b.core_pos)
        };

        // readaln2's hard filters.
        let f = self.filters;
        if (f.incl_flags != 0 && flag & f.incl_flags == 0)
            || (f.excl_flags != 0 && flag & f.excl_flags != 0)
            || i32::from(mapq) < f.min_mqual
        {
            return Ok(());
        }

        if flag & FUNMAP != 0 || tid == -1 {
            return Ok(());
        }

        let pos = if tid == self.last_ref {
            core_pos + 1
        } else {
            POS_MAX
        };

        if self.col > pos {
            bail!("BAM/SAM file is not sorted by position. Aborting");
        }

        self.columns(pos, handler)?;
        self.col = pos;

        if tid != self.last_ref {
            self.last_ref = tid;
            self.nth = 0;
            self.col = core_pos + 1;
        }

        let p = &mut self.entries[pnew];
        p.eof = 0;
        p.start = 2;
        p.pos = core_pos;
        p.cigar_ind = 0;
        p.cigar_len = 0;
        p.cigar_op = -1;
        p.seq_offset = -1;
        p.first_del = false;
        p.qual = 0;
        p.is_rev = p.flag & FREVERSE != 0;

        if handler.init(p)? {
            self.active.push(pnew);
        } else {
            self.free.push(pnew);
        }

        self.pnew = match self.free.pop() {
            Some(i) => i,
            None => {
                self.entries.push(Entry::default());
                self.entries.len() - 1
            }
        };

        Ok(())
    }

    /// Processes the remaining columns at the end of the input.
    pub fn finish<H: Handler>(&mut self, handler: &mut H) -> anyhow::Result<()> {
        self.columns(POS_MAX, handler)
    }

    /// Processes the columns before `pos`.
    fn columns<H: Handler>(&mut self, pos: i64, handler: &mut H) -> anyhow::Result<()> {
        while self.col < pos && self.col <= self.last_col && !self.active.is_empty() {
            let mut is_insert = 0;
            let mut ins = 0;
            self.removed.clear();

            for &i in &self.active {
                let p = &mut self.entries[i];

                if p.next_base(self.col, self.nth, &mut ins) == 0 {
                    p.eof = 1;
                }

                if p.eof == 1 {
                    self.removed.push(i);
                }

                is_insert = is_insert.max(ins);
            }

            let column = Column {
                entries: &self.entries,
                active: &self.active,
                tid: self.entries[self.active[0]].tid,
                pos: self.col,
                nth: self.nth,
                depth: i32::try_from(self.active.len())?,
            };

            let result = handler.column(&column);

            // Remove the finished reads.
            if !self.removed.is_empty() {
                let entries = &mut self.entries;
                self.active.retain(|&i| entries[i].eof != 1);

                for &i in &self.removed {
                    self.entries[i].has_nm = false;
                    self.free.push(i);
                }
            }

            result?;

            if is_insert != 0 {
                self.nth += 1;
            } else {
                let had_insertions = self.nth != 0;
                self.nth = 0;
                self.col += 1;

                if !had_insertions && self.col < self.first_col.min(pos) {
                    self.skip_columns(self.first_col.min(pos));
                }
            }
        }

        Ok(())
    }

    /// Steps over unobserved columns from `col` towards `limit` (exclusive) while every read
    /// stays inside its current CIGAR operation, just after a column without insertions.
    ///
    /// Inside a match, deletion or reference skip, a column changes a read only by moving it one
    /// base on (with `start` counting down), as long as `eof` has already settled: 0 in a match
    /// or deletion, 2 or 3 in a skip. So the reads are moved `n` bases on here, and the next
    /// column, still inside each operation, is visited as usual and sets the rest (the base and
    /// quality, which in a deletion is the same minimum, and `eof` 2 in a skip).
    fn skip_columns(&mut self, limit: i64) {
        let mut n = limit - self.col;

        for &i in &self.active {
            let p = &self.entries[i];

            let settled = match p.cigar_op {
                CMATCH | CEQUAL | CDIFF | CDEL => p.eof == 0,
                CREF_SKIP => p.eof != 0,
                _ => false,
            };

            if !settled {
                return;
            }

            // Keep at least one base of the operation for the next column.
            n = n.min(i64::from(p.cigar_len) - 1);

            if n <= 0 {
                return;
            }
        }

        let Ok(n32) = i32::try_from(n) else {
            return;
        };

        for &i in &self.active {
            let p = &mut self.entries[i];
            p.pos += n;
            p.cigar_len -= n32;
            p.start = (p.start - n32).max(0);

            if matches!(p.cigar_op, CMATCH | CEQUAL | CDIFF) {
                p.seq_offset += n32;
            }
        }

        self.col += n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encodes a minimal BAM record.
    fn record(pos: i32, cigar: &[(u32, i32)], seq: &[u8], qual: &[u8]) -> Vec<u8> {
        let mut raw = Vec::new();
        raw.extend_from_slice(&0i32.to_le_bytes());
        raw.extend_from_slice(&pos.to_le_bytes());
        raw.push(2); // l_read_name
        raw.push(60);
        raw.extend_from_slice(&0u16.to_le_bytes());
        raw.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
        raw.extend_from_slice(&0u16.to_le_bytes());
        raw.extend_from_slice(&(seq.len() as u32).to_le_bytes());
        raw.extend_from_slice(&(-1i32).to_le_bytes());
        raw.extend_from_slice(&(-1i32).to_le_bytes());
        raw.extend_from_slice(&0i32.to_le_bytes());
        raw.extend_from_slice(b"r\0");

        for &(len, op) in cigar {
            raw.extend_from_slice(&((len << 4) | op as u32).to_le_bytes());
        }

        for pair in seq.chunks(2) {
            let code = |b: u8| NT16.iter().position(|&c| c == b).unwrap_or(15) as u8;
            let hi = code(pair[0]);
            let lo = pair.get(1).map_or(0, |&b| code(b));
            raw.push((hi << 4) | lo);
        }

        raw.extend_from_slice(qual);
        raw
    }

    #[derive(Default)]
    struct Collect(Vec<(i64, i32, String)>);

    impl Handler for Collect {
        fn init(&mut self, _: &mut Entry) -> anyhow::Result<bool> {
            Ok(true)
        }

        fn column(&mut self, column: &Column<'_>) -> anyhow::Result<()> {
            let bases = column.reads().map(|p| char::from(p.base)).collect();
            self.0.push((column.pos, column.nth, bases));
            Ok(())
        }
    }

    #[test]
    fn test_insertion_and_deletion_columns() -> anyhow::Result<()> {
        let mut pileup = Pileup::new(Filters::default());
        let mut collect = Collect::default();

        // 2M1I2M and 1M2D2M, both at position 1.
        let a = record(
            0,
            &[(2, CMATCH), (1, CINS), (2, CMATCH)],
            b"ACGTA",
            &[30; 5],
        );
        let b = record(0, &[(1, CMATCH), (2, CDEL), (2, CMATCH)], b"ATA", &[20; 3]);
        pileup.push(&a, &mut collect)?;
        pileup.push(&b, &mut collect)?;
        pileup.finish(&mut collect)?;

        let expected = [
            (1, 0, "AA"),
            (2, 0, "C*"),
            (2, 1, "G*"),
            (3, 0, "T*"),
            (4, 0, "AT"),
            (5, 0, "A"),
        ];

        let actual: Vec<_> = collect
            .0
            .iter()
            .map(|(pos, nth, s)| (*pos, *nth, s.as_str()))
            .collect();
        assert_eq!(actual, expected);

        Ok(())
    }

    /// A column's position, `nth`, and its reads' bases, qualities and reference skip flags.
    type ColumnState = (i64, i32, Vec<(u8, i32, bool)>);

    /// Every column.
    #[derive(Default)]
    struct CollectAll(Vec<ColumnState>);

    impl Handler for CollectAll {
        fn init(&mut self, _: &mut Entry) -> anyhow::Result<bool> {
            Ok(true)
        }

        fn column(&mut self, column: &Column<'_>) -> anyhow::Result<()> {
            let reads = column
                .reads()
                .map(|p| (p.base, p.qual, p.ref_skip))
                .collect();
            self.0.push((column.pos, column.nth, reads));
            Ok(())
        }
    }

    #[test]
    fn test_window_gives_the_same_columns() -> anyhow::Result<()> {
        let quals: Vec<u8> = (0..40u8).map(|i| 10 + (i % 6) * 5).collect();
        let records = [
            record(
                0,
                &[
                    (3, CMATCH),
                    (40, CREF_SKIP),
                    (2, CMATCH),
                    (1, CINS),
                    (3, CMATCH),
                ],
                b"ACGTTACGA",
                &quals[..9],
            ),
            record(
                1,
                &[(2, CSOFT_CLIP), (4, CMATCH), (30, CDEL), (5, CMATCH)],
                b"GGACGTACGTA",
                &quals[..11],
            ),
            record(
                3,
                &[
                    (2, CMATCH),
                    (2, CINS),
                    (6, CMATCH),
                    (1, CPAD),
                    (1, CINS),
                    (9, CMATCH),
                ],
                b"ACGTACGTACGTACGTACGT",
                &quals[..20],
            ),
            record(
                5,
                &[
                    (1, CMATCH),
                    (20, CREF_SKIP),
                    (1, CINS),
                    (2, CDEL),
                    (4, CMATCH),
                ],
                b"TACGT",
                &quals[..5],
            ),
            record(
                20,
                &[(10, CMATCH), (3, CREF_SKIP), (10, CMATCH)],
                b"ACGTACGTACGTACGTACGT",
                &quals[..20],
            ),
        ];

        let run_on =
            |records: &[Vec<u8>], window: Option<(i64, i64)>| -> anyhow::Result<CollectAll> {
                let mut pileup = Pileup::new(Filters::default());
                let mut collect = CollectAll::default();

                if let Some((first, last)) = window {
                    pileup.set_window(first, last);
                }

                for raw in records {
                    pileup.push(raw, &mut collect)?;
                }

                pileup.finish(&mut collect)?;
                Ok(collect)
            };
        let run = |window| run_on(&records, window);

        let all = run(None)?.0;

        for first in 1..60 {
            for last in first..60 {
                let expected: Vec<_> = all
                    .iter()
                    .filter(|(pos, _, _)| (first..=last).contains(pos))
                    .collect();
                let windowed = run(Some((first, last)))?.0;
                let actual: Vec<_> = windowed
                    .iter()
                    .filter(|(pos, _, _)| (first..=last).contains(pos))
                    .collect();
                assert_eq!(actual, expected, "window {first}..={last}");
                assert!(windowed.iter().all(|(pos, _, _)| *pos <= last));
            }
        }

        // The reference skip of a read on its own is stepped over.
        let skipped = run_on(&records[..1], Some((40, 50)))?.0;
        assert!(skipped.iter().filter(|c| c.0 < 40).count() < 5);

        Ok(())
    }

    #[test]
    fn test_unsorted_input_is_an_error() -> anyhow::Result<()> {
        let mut pileup = Pileup::new(Filters::default());
        let mut collect = Collect::default();

        pileup.push(
            &record(10, &[(5, CMATCH)], b"ACGTA", &[30; 5]),
            &mut collect,
        )?;
        pileup.push(
            &record(20, &[(5, CMATCH)], b"ACGTA", &[30; 5]),
            &mut collect,
        )?;
        assert!(
            pileup
                .push(&record(5, &[(5, CMATCH)], b"ACGTA", &[30; 5]), &mut collect)
                .is_err()
        );

        Ok(())
    }
}
