//! The repeat context of indels from a reference (bcftools stats `-F`; `indel_ctx_t` in
//! `vcfstats.c`): the most repeated k-mer next to an indel, to tell repeat-consistent indels
//! from inconsistent ones.

use std::{collections::HashMap, path::Path};

use pasta::variant::log;

use crate::faidx::Faidx;

/// The window after the indel's position that is searched for repeats.
const WIN_SIZE: usize = 50;
/// The longest repeat element counted (bcftools' `IRC_RLEN`).
pub const REP_LEN: usize = 10;

/// One sequence of a FASTA index.
#[derive(Clone, Copy, Debug)]
struct FaiEntry {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// A k-mer seen in the window (`_idc1_t`).
#[derive(Clone, Debug)]
struct Kmer {
    seq: Vec<u8>,
    pos: usize,
    cnt: u32,
}

/// The repeat-context finder: a FASTA file with its index, and k-mer scratch space.
pub struct IndelContext {
    path: String,
    faidx: Faidx,
    index: HashMap<Vec<u8>, FaiEntry>,
    kmers: Vec<Kmer>,
    buf: Vec<u8>,
    raw: Vec<u8>,
}

/// The repeat context of an indel (`indel_ctx_type`'s outputs).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Context {
    /// The indel length: negative for deletions.
    pub ndel: i32,
    /// How many times the repeat element occurs.
    pub nrep: u32,
    /// The length of the repeat element.
    pub nlen: usize,
}

/// A reference base that disagrees with REF (bcftools exits with this message).
#[derive(Debug)]
pub struct SanityError(pub String);

impl std::fmt::Debug for IndelContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IndelContext")
            .field("path", &self.path)
            .finish_non_exhaustive()
    }
}

impl IndelContext {
    /// Opens a FASTA file and its index, building the index when it is missing
    /// (`indel_ctx_init`). Returns `None`, after HTSlib's messages, when that fails.
    pub fn open(path: &str) -> Option<Self> {
        let faidx = Faidx::load(Path::new(path)).ok()?;
        let index = faidx
            .index()
            .as_ref()
            .iter()
            .map(|r| {
                (
                    r.name().to_vec(),
                    FaiEntry {
                        len: r.length(),
                        offset: r.position(),
                        line_bases: r.line_base_count().get(),
                        line_width: r.line_width().get(),
                    },
                )
            })
            .collect();
        Some(Self {
            path: path.to_owned(),
            faidx,
            index,
            kmers: Vec::new(),
            buf: Vec::new(),
            raw: Vec::new(),
        })
    }

    /// Opens the reference again, for another thread (`None` when it cannot be opened).
    pub fn reopen(&self) -> Option<Self> {
        Some(Self {
            path: self.path.clone(),
            faidx: self.faidx.reopen().ok()?,
            index: self.index.clone(),
            kmers: Vec::new(),
            buf: Vec::new(),
            raw: Vec::new(),
        })
    }

    /// Reads bases `beg..=end` (0-based) of a sequence into `self.buf`, clamped to the
    /// sequence as `faidx_fetch_seq` clamps them. Returns `false` when the sequence is absent.
    fn fetch(&mut self, chr: &[u8], mut beg: i64, mut end: i64) -> bool {
        let Some(e) = self.index.get(chr).copied() else {
            log::error(
                "faidx_adjust_position",
                format_args!(
                    "The sequence \"{}\" was not found",
                    String::from_utf8_lossy(chr)
                ),
            );
            return false;
        };
        let len = e.len as i64;
        if end < beg {
            beg = end;
        }
        beg = beg.clamp(0, len);
        if end < 0 {
            end = 0;
        } else if len <= end {
            end = len - 1;
        }
        self.buf.clear();
        let (beg, end) = (beg as u64, (end + 1).max(0) as u64);
        if end <= beg {
            return true;
        }
        let at = |p: u64| e.offset + p / e.line_bases * e.line_width + p % e.line_bases;
        let first = at(beg);
        let last = at(end - 1);
        self.raw.resize((last - first + 1) as usize, 0);
        if self.faidx.read_at(first, &mut self.raw).is_err() {
            return true;
        }
        for p in beg..end {
            self.buf.push(self.raw[(at(p) - first) as usize]);
        }
        true
    }

    /// Finds the repeat context of the indel `alt` at 1-based `pos` (`indel_ctx_type`).
    ///
    /// # Errors
    ///
    /// Fails when the reference disagrees with REF.
    pub fn classify(
        &mut self,
        chr: &[u8],
        pos: i64,
        ref_: &[u8],
        alt: &[u8],
    ) -> Result<Context, SanityError> {
        let ref_len = ref_.len();
        let alt_len = alt.iter().position(|&c| c == b',').unwrap_or(alt.len());
        let ndel = alt_len as i32 - ref_len as i32;
        if !self.fetch(chr, pos - 1, pos + WIN_SIZE as i64) {
            return Ok(Context {
                ndel,
                nrep: 0,
                nlen: 0,
            });
        }
        for c in &mut self.buf {
            if *c > 96 {
                *c -= 32;
            }
        }
        for (i, (&r, &f)) in ref_.iter().zip(&self.buf).enumerate() {
            if r != f && r.wrapping_sub(32) != f && !iupac_consistent(f, r) {
                return Err(SanityError(format!(
                    "\nSanity check failed, the reference sequence differs: {}:{}+{} .. {} vs {}\n",
                    String::from_utf8_lossy(chr),
                    pos,
                    i,
                    char::from(r),
                    char::from(f)
                )));
            }
        }
        // Beyond the fetched bases, bcftools reads the string terminator (and past it).
        self.buf.resize(WIN_SIZE + 2, 0);
        self.kmers.clear();
        for i in 0..WIN_SIZE {
            let kmax = if REP_LEN <= i { REP_LEN } else { i + 1 };
            for k in 0..kmax {
                let start = i + 1 - k;
                self.insert(start, k + 1, i - k);
            }
        }
        let (mut max_cnt, mut max_len) = (0u32, 0usize);
        for k in &self.kmers {
            if max_cnt < k.cnt || (max_cnt == k.cnt && max_len < k.seq.len()) {
                max_cnt = k.cnt;
                max_len = k.seq.len();
            }
        }
        Ok(Context {
            ndel,
            nrep: max_cnt,
            nlen: max_len,
        })
    }

    /// Binary search for a k-mer (`_indel_ctx_lookup`): its index when found, else the index
    /// of the last smaller one (`-1` when none).
    fn lookup(&self, seq: &[u8]) -> (isize, bool) {
        let (mut min, mut max) = (0isize, self.kmers.len() as isize - 1);
        while min <= max {
            let i = (min + max) / 2;
            let dat = &self.kmers[i as usize].seq;
            match strncmp(seq, dat, seq.len()) {
                std::cmp::Ordering::Less => max = i - 1,
                std::cmp::Ordering::Greater => min = i + 1,
                std::cmp::Ordering::Equal => {
                    if seq.len() == dat.len() {
                        return (i, true);
                    } else if seq.len() < dat.len() {
                        max = i - 1;
                    } else {
                        min = i + 1;
                    }
                }
            }
        }
        (max, false)
    }

    /// Counts the k-mer `buf[start..start + len]` ending at window offset `pos`
    /// (`_indel_ctx_insert`): consecutive tandem copies increase the count.
    fn insert(&mut self, start: usize, len: usize, pos: usize) {
        let seq = &self.buf[start..start + len];
        let (idat, hit) = self.lookup(seq);
        if !hit {
            if pos > 0 {
                return;
            }
            let at = (idat + 1) as usize;
            let kmer = Kmer {
                seq: seq.to_vec(),
                pos,
                cnt: 1,
            };
            self.kmers.insert(at, kmer);
            return;
        }
        let k = &mut self.kmers[idat as usize];
        if k.pos + len == pos {
            k.cnt += 1;
            k.pos = pos;
        }
    }
}

/// C's `strncmp` on byte strings that end at their first NUL (or their end).
fn strncmp(a: &[u8], b: &[u8], n: usize) -> std::cmp::Ordering {
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0);
        let y = b.get(i).copied().unwrap_or(0);
        if x != y {
            return x.cmp(&y);
        }
        if x == 0 {
            break;
        }
    }
    std::cmp::Ordering::Equal
}

/// Whether base `nt` is one of the bases of IUPAC code `iupac` (bcftools' `iupac_consistent`).
fn iupac_consistent(iupac: u8, nt: u8) -> bool {
    const IUPAC_MASK: [u8; 90] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 0, 7, 9, 0,
        10,
    ];
    if iupac > 89 {
        return false;
    }
    let mut nt = nt;
    if nt > 90 {
        nt = nt.wrapping_sub(32);
    }
    let nt = match nt {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        other => other,
    };
    IUPAC_MASK[usize::from(iupac)] & nt != 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac() {
        assert!(iupac_consistent(b'R', b'A'));
        assert!(iupac_consistent(b'R', b'g'));
        assert!(!iupac_consistent(b'R', b'C'));
        assert!(iupac_consistent(b'N', b'T'));
    }

    #[test]
    fn strncmp_stops_at_nul() {
        use std::cmp::Ordering;
        assert_eq!(strncmp(b"AC", b"AC", 2), Ordering::Equal);
        assert_eq!(strncmp(b"AC", b"A", 2), Ordering::Greater);
        assert_eq!(strncmp(b"A\0X", b"A\0Y", 3), Ordering::Equal);
    }
}
