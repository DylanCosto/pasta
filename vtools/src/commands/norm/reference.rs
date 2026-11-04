//! Random access to an indexed FASTA reference (HTSlib's `faidx_t`), with a window cache so
//! that the many short, nearby fetches of `norm` read the file only once per region.

use std::{collections::HashMap, path::Path};

use anyhow::{Context, bail};
use pasta::variant::log;

use crate::faidx::Faidx;

/// How many bases one cache fill reads at least.
const CACHE_BASES: i64 = 64 * 1024;

/// How many bases before the requested start a cache fill also reads, for the leftward fetches
/// of left-alignment.
const CACHE_BACK: i64 = 1024;

/// A sequence's index entry.
struct Entry {
    length: i64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// What [`Reference::fetch_with`] does about a sequence missing from the index.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MissingSeq {
    /// Prints HTSlib's message.
    Report,
    /// Prints nothing, as with `hts_verbose = 0`.
    Silent,
}

/// An indexed FASTA file (HTSlib's `fai_load`/`faidx_fetch_seq64`).
pub struct Reference {
    faidx: Faidx,
    entries: Vec<Entry>,
    names: HashMap<Vec<u8>, usize>,
    /// The sequence names in index order.
    order: Vec<Vec<u8>>,
    /// The cached window: sequence index, first base, and the bases.
    cache_seq: Option<usize>,
    cache_beg: i64,
    cache: Vec<u8>,
    raw: Vec<u8>,
}

impl Reference {
    /// Opens a FASTA file and its index, building the index when it is missing, as
    /// HTSlib's `fai_load` does.
    ///
    /// # Errors
    ///
    /// Returns bcftools' error, after HTSlib's messages, if the file or its index cannot be
    /// read or built.
    pub fn open(path: &Path) -> anyhow::Result<Self> {
        let Ok(faidx) = Faidx::load(path) else {
            bail!("Failed to load the fai index: {}", path.display());
        };
        let mut entries = Vec::new();
        let mut names = HashMap::new();
        let mut order = Vec::new();
        for rec in faidx.index().as_ref() {
            names.insert(rec.name().to_vec(), entries.len());
            order.push(rec.name().to_vec());
            entries.push(Entry {
                length: rec.length() as i64,
                offset: rec.position(),
                line_bases: rec.line_base_count().get(),
                line_width: rec.line_width().get(),
            });
        }
        Ok(Self {
            faidx,
            entries,
            names,
            order,
            cache_seq: None,
            cache_beg: 0,
            cache: Vec::new(),
            raw: Vec::new(),
        })
    }

    /// Fetches the bases `beg..=end` (0-based) of a sequence into `out`, clamped to the
    /// sequence as HTSlib's `faidx_fetch_seq64` clamps them. Returns `false` (after printing
    /// HTSlib's error) when the sequence is not in the index.
    pub fn fetch(
        &mut self,
        name: &[u8],
        beg: i64,
        end: i64,
        out: &mut Vec<u8>,
    ) -> anyhow::Result<bool> {
        self.fetch_with(name, beg, end, out, MissingSeq::Report)
    }

    /// Whether the index has a sequence (HTSlib's `faidx_has_seq`).
    pub fn has_seq(&self, name: &[u8]) -> bool {
        self.names.contains_key(name)
    }

    /// The name of the `i`-th sequence in index order (HTSlib's `faidx_iseq`).
    pub fn seq_name(&self, i: usize) -> Option<&[u8]> {
        self.order.get(i).map(Vec::as_slice)
    }

    /// [`Self::fetch`], printing HTSlib's message for a missing sequence only with
    /// [`MissingSeq::Report`].
    pub fn fetch_with(
        &mut self,
        name: &[u8],
        beg: i64,
        end: i64,
        out: &mut Vec<u8>,
        missing: MissingSeq,
    ) -> anyhow::Result<bool> {
        out.clear();
        let Some(&iseq) = self.names.get(name) else {
            if missing == MissingSeq::Silent {
                return Ok(false);
            }
            log::error(
                "faidx_adjust_position",
                format_args!(
                    "The sequence \"{}\" was not found",
                    String::from_utf8_lossy(name)
                ),
            );
            return Ok(false);
        };
        let len = self.entries[iseq].length;
        let (mut beg, mut end) = (beg, end);
        if end < beg {
            beg = end;
        }
        if beg < 0 {
            beg = 0;
        } else if len <= beg {
            beg = len;
        }
        if end < 0 {
            end = 0;
        } else if len <= end {
            end = len - 1;
        }
        let end = end + 1; // exclusive
        if end <= beg {
            return Ok(true);
        }
        let cached = self.cache_seq == Some(iseq)
            && beg >= self.cache_beg
            && end <= self.cache_beg + self.cache.len() as i64;
        if !cached {
            let fill_beg = (beg - CACHE_BACK).max(0);
            let fill_end = end.max(beg + CACHE_BASES).min(len);
            self.fill(iseq, fill_beg, fill_end)?;
        }
        let off = (beg - self.cache_beg) as usize;
        out.extend_from_slice(&self.cache[off..off + (end - beg) as usize]);
        Ok(true)
    }

    /// Reads the bases `beg..end` of a sequence into the cache.
    fn fill(&mut self, iseq: usize, beg: i64, end: i64) -> anyhow::Result<()> {
        let e = &self.entries[iseq];
        if e.line_bases == 0 {
            bail!("[E::fai_retrieve] Invalid line length in index: 0");
        }
        let file_pos = |b: i64| -> u64 {
            let b = b as u64;
            e.offset + b / e.line_bases * e.line_width + b % e.line_bases
        };
        let start = file_pos(beg);
        let stop = file_pos(end - 1) + 1;
        self.raw.resize((stop - start) as usize, 0);
        self.faidx
            .read_at(start, &mut self.raw)
            .context("[E::fai_retrieve] Failed to retrieve block: unexpected end of file")?;
        self.cache.clear();
        let mut start = 0;
        for i in memchr::memchr2_iter(b'\n', b'\r', &self.raw) {
            self.cache.extend_from_slice(&self.raw[start..i]);
            start = i + 1;
        }
        self.cache.extend_from_slice(&self.raw[start..]);
        if self.cache.len() < (end - beg) as usize {
            bail!("[E::fai_retrieve] Failed to retrieve block: unexpected end of file");
        }
        self.cache.truncate((end - beg) as usize);
        self.cache_seq = Some(iseq);
        self.cache_beg = beg;
        Ok(())
    }
}
