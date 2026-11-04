//! gVCF reference blocks, `bcftools call -g` (bcftools' `gvcf.c`).
//!
//! `gvcf_write` is split in two: [`GvcfRanges::classify`] looks at one called record (and
//! adds INFO/MIN_DP to reference sites written as they are), which needs no state and runs on
//! worker threads; [`Gvcf::advance`] then extends, starts, or flushes the open block in
//! input order.

use anyhow::{Context, bail};
use pasta::variant::{Header, Record};

use super::{Tag, get_code};

/// The minimum depths of `-g` (`dp_range`), and the tags gVCF output reads and writes.
#[derive(Clone, Debug)]
pub struct GvcfRanges {
    dp_range: Vec<i32>,
    tags: Option<GTags>,
}

/// The tags of gVCF blocks, resolved in the output header.
#[derive(Clone, Copy, Debug)]
struct GTags {
    dp: Tag,
    pl: Tag,
    qs: Tag,
    end: Tag,
    min_dp: Tag,
}

/// What the caller made of a site, as far as gVCF blocks are concerned.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CalledSite {
    /// A reference site, which may join a block.
    Reference,
    /// Any other site, written as it is.
    Other,
}

impl CalledSite {
    /// From a caller's return value, which is 1 for a reference site.
    pub fn from_call(ret: i32) -> Self {
        if ret == 1 {
            Self::Reference
        } else {
            Self::Other
        }
    }
}

/// What a reference site brings to a block.
#[derive(Clone, Debug, Default)]
pub struct GvcfSite {
    min_dp: i32,
    dp_range: usize,
    /// FORMAT/DP of each sample (whatever was read last, when absent).
    dp: Vec<i32>,
    /// FORMAT/PL, as `bcf_get_format_int32` returns it (count or error code).
    npl: i32,
    pl: Vec<i32>,
    nqsum: i32,
    qsum: Vec<f32>,
    ngts: i32,
    gts: Vec<i32>,
    /// The alleles, comma-separated.
    als: Vec<u8>,
    /// INFO/END (0-based), or the position.
    end: i64,
    /// Scratch space for reading INFO/END.
    end_buf: Vec<i32>,
}

/// A called record as gVCF output sees it.
#[derive(Debug)]
pub enum GvcfItem {
    /// A reference site that joins a block (the record itself is not written).
    Collapse {
        rid: i32,
        pos: i64,
        site: Box<GvcfSite>,
    },
    /// A record written as it is, after flushing the block if needed.
    Write { rid: i32, pos: i64 },
}

impl GvcfRanges {
    /// Parses the comma-separated minimum depths of `-g` (`gvcf_init`); `None` when they
    /// cannot be parsed.
    pub fn new(dp_ranges: &str) -> Option<Self> {
        let s = dp_ranges.as_bytes();
        if s.is_empty() {
            // gvcf_init sizes the list by commas, one entry here, and leaves it unset
            // (uninitialized memory, 0 in practice).
            return Some(Self {
                dp_range: vec![0],
                tags: None,
            });
        }
        let mut dp_range = Vec::new();
        let mut i = 0;
        while i < s.len() {
            let (v, end) = pasta::variant::number::strtol(s, i);
            if end == i {
                return None;
            }
            dp_range.push(v as i32);
            if s.get(end) == Some(&b',') && end + 1 < s.len() {
                i = end + 1;
                continue;
            } else if end == s.len() {
                break;
            }
            return None;
        }
        Some(Self {
            dp_range,
            tags: None,
        })
    }

    /// Resolves the tags in the output header; needed before any record is classified.
    pub fn resolve(&mut self, hdr: &Header) {
        self.tags = Some(GTags {
            dp: Tag::new(hdr, b"DP"),
            pl: Tag::new(hdr, b"PL"),
            qs: Tag::new(hdr, b"QS"),
            end: Tag::new(hdr, b"END"),
            min_dp: Tag::new(hdr, b"MIN_DP"),
        });
    }

    fn tags(&self) -> GTags {
        self.tags
            .expect("GvcfRanges::resolve runs before records are classified")
    }

    /// The per-record part of `gvcf_write`: whether a called record joins a block, and for
    /// a reference site written as is, INFO/MIN_DP. A site that joins a block is read into
    /// `spare` when one is given (its buffers are reused), or into a new site.
    ///
    /// # Errors
    ///
    /// Returns an error when INFO/MIN_DP cannot be set.
    pub fn classify(
        &self,
        hdr: &Header,
        rec: &mut Record,
        called: CalledSite,
        dp: &mut Vec<i32>,
        spare: &mut Option<Box<GvcfSite>>,
    ) -> anyhow::Result<GvcfItem> {
        let tags = self.tags();
        let nsmpl = hdr.sample_count();
        let (rid, pos) = (rec.rid(), rec.pos());
        if called == CalledSite::Other {
            return Ok(GvcfItem::Write { rid, pos });
        }
        let mut can_collapse = true;
        let mut dp_range = 0usize;
        let mut min_dp = 0i32;
        let ret = tags.dp.get_format_i32(hdr, rec, dp);
        if ret == nsmpl as i32 && nsmpl > 0 {
            min_dp = dp[..nsmpl].iter().copied().min().unwrap_or(0);
            dp_range = self
                .dp_range
                .iter()
                .position(|&r| min_dp < r)
                .unwrap_or(self.dp_range.len());
            if dp_range == 0 {
                can_collapse = false;
            }
        }
        if !can_collapse {
            if min_dp != 0 {
                tags.min_dp
                    .update_info_i32(hdr, rec, &[min_dp])
                    .context("Could not update INFO/MIN_DP")?;
            }
            return Ok(GvcfItem::Write { rid, pos });
        }
        let mut site = spare.take().unwrap_or_default();
        site.min_dp = min_dp;
        site.dp_range = dp_range;
        site.dp.clear();
        site.dp.extend_from_slice(dp);
        site.dp.resize(nsmpl, 0);
        site.npl = tags.pl.get_format_i32(hdr, rec, &mut site.pl);
        site.nqsum = tags.qs.get_info_f32(hdr, rec, &mut site.qsum);
        site.ngts = get_code(rec.get_genotypes(hdr, &mut site.gts));
        site.als.clear();
        for (i, a) in rec.alleles().enumerate() {
            if i > 0 {
                site.als.push(b',');
            }
            site.als.extend_from_slice(a);
        }
        site.end = if tags.end.get_info_i32(hdr, rec, &mut site.end_buf) == 1 {
            i64::from(site.end_buf[0]) - 1
        } else {
            pos
        };
        Ok(GvcfItem::Collapse { rid, pos, site })
    }
}

/// The open block (bcftools' `gvcf_t`).
#[derive(Debug)]
pub struct Gvcf {
    pub ranges: GvcfRanges,
    /// The DP range of the open block, 0 when no block is open.
    prev_range: usize,
    block: Box<GvcfSite>,
    rid: i32,
    start: i64,
    line: Record,
}

impl Gvcf {
    /// No open block.
    pub fn new(ranges: GvcfRanges) -> Self {
        Self {
            ranges,
            prev_range: 0,
            block: Box::default(),
            rid: -1,
            start: 0,
            line: Record::new(),
        }
    }

    /// Adds the INFO/END and INFO/MIN_DP definitions (`gvcf_update_header`).
    ///
    /// # Errors
    ///
    /// Returns an error when the header lines conflict with existing ones.
    pub fn update_header(hdr: &mut Header) -> anyhow::Result<()> {
        hdr.append(b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">")?;
        hdr.append(b"##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum per-sample depth in this gVCF block\">")?;
        Ok(())
    }

    /// The stateful part of `gvcf_write`: flushes the open block (through `write`) when the
    /// item does not extend it, then opens or extends a block with a reference site. With no
    /// item, flushes what is left. Returns a site whose buffers can be reused by
    /// [`GvcfRanges::classify`].
    ///
    /// # Errors
    ///
    /// Returns the errors of `write`, and an error for unexpected PL counts.
    pub fn advance<W>(
        &mut self,
        hdr: &Header,
        item: Option<GvcfItem>,
        mut write: W,
    ) -> anyhow::Result<Option<Box<GvcfSite>>>
    where
        W: FnMut(&mut Record) -> anyhow::Result<()>,
    {
        let nsmpl = hdr.sample_count();
        if item.is_none() && self.prev_range == 0 {
            return Ok(None);
        }
        let has_rec = item.is_some();
        let (rid, pos, site) = match item {
            None => (0, 0, None),
            Some(GvcfItem::Write { rid, pos }) => (rid, pos, None),
            Some(GvcfItem::Collapse { rid, pos, site }) => (rid, pos, Some(site)),
        };
        let mut needs_flush = site.is_none();
        let dp_range = site.as_ref().map_or(0, |s| s.dp_range);
        if self.prev_range != 0 && self.prev_range != dp_range {
            needs_flush = true;
        }
        if !has_rec || self.rid != rid || pos > self.block.end + 1 {
            needs_flush = true;
        }
        if self.prev_range != 0 && needs_flush {
            if has_rec && rid == self.rid && pos == self.block.end {
                self.block.end -= 1;
            }
            self.block.end += 1;
            self.flush(hdr, nsmpl, &mut write)?;
        }
        let Some(site) = site else {
            return Ok(None);
        };
        let dp_range = site.dp_range;
        let spare = if self.prev_range == 0 {
            self.rid = rid;
            self.start = pos;
            std::mem::replace(&mut self.block, site)
        } else {
            let b = &mut self.block;
            if b.min_dp > site.min_dp {
                b.min_dp = site.min_dp;
            }
            for (d, &t) in b.dp.iter_mut().zip(&site.dp) {
                if *d > t {
                    *d = t;
                }
            }
            if site.npl >= 0 {
                if site.npl != (nsmpl * 3) as i32 {
                    bail!("Unexpected number of PL fields");
                }
                // bcftools reads a PL buffer that was never filled when the block's first
                // site had no PL (and crashes); such blocks keep no PL here.
                // When the block's first site has another number of PLs (one per sample at
                // a site with only REF, say), bcftools reads and writes past them with
                // three per sample; only the values the block has are written out.
                if b.npl > 0 {
                    for i in 0..nsmpl {
                        let (k1, k2) = (3 * i + 1, 3 * i + 2);
                        let Some(&b1) = b.pl.get(k1) else {
                            break;
                        };
                        let (p1, p2) = (site.pl[k1], site.pl[k2]);
                        if b1 > p1 {
                            b.pl[k1] = p1;
                            if let Some(b2) = b.pl.get_mut(k2) {
                                *b2 = p2;
                            }
                        } else if b1 == p1
                            && let Some(b2) = b.pl.get_mut(k2)
                            && *b2 > p2
                        {
                            *b2 = p2;
                        }
                    }
                }
            } else {
                b.npl = 0;
            }
            b.end = site.end;
            site
        };
        self.prev_range = dp_range;
        Ok(Some(spare))
    }

    /// Writes the open block and closes it.
    fn flush<W>(&mut self, hdr: &Header, nsmpl: usize, write: &mut W) -> anyhow::Result<()>
    where
        W: FnMut(&mut Record) -> anyhow::Result<()>,
    {
        let b = &self.block;
        let line = &mut self.line;
        line.clear();
        line.set_rid(self.rid);
        line.set_pos(self.start);
        line.set_rlen(b.end - self.start);
        let tags = self.ranges.tags();
        line.update_alleles_str(hdr, &b.als);
        if self.start + 1 < b.end {
            tags.end
                .update_info_i32(hdr, line, &[b.end as i32])
                .context("Could not update INFO/END")?;
        }
        tags.min_dp
            .update_info_i32(hdr, line, &[b.min_dp])
            .context("Could not update INFO/MIN_DP")?;
        if b.nqsum > 0 {
            tags.qs
                .update_info_f32(hdr, line, &b.qsum[..b.nqsum as usize])
                .context("Could not update INFO/QS")?;
        }
        if b.ngts > 0 {
            line.update_genotypes(hdr, &b.gts[..b.ngts as usize])
                .context("Could not update FORMAT/GT")?;
        }
        if b.npl > 0 {
            tags.pl
                .update_format_i32(hdr, line, &b.pl[..b.npl as usize])
                .context("Could not update FORMAT/PL")?;
        }
        tags.dp
            .update_format_i32(hdr, line, &b.dp[..nsmpl])
            .context("Could not update FORMAT/DP")?;
        write(line)?;
        self.prev_range = 0;
        self.rid = -1;
        self.block.npl = 0;
        self.block.nqsum = 0;
        self.block.ngts = 0;
        Ok(())
    }
}
