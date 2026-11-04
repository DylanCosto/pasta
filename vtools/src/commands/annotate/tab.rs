//! Annotating from tab-delimited files: lines streamed from a tabix-indexed file
//! (`annotate_from_tab`, with `buffer_annot_lines`) or held in a region index when the file
//! has END columns (`annotate_from_regidx`), and the `{}` values of `-i/-e` expressions.

use anyhow::bail;
use pasta::variant::{
    Header, Record, Unpack,
    number::{strtod, strtol},
    regidx::RegIdx,
    synced::{Regions, RegionsError},
};

use super::{
    columns::{Column, MergeMethod, TabColumns},
    line::AnnotLine,
    setters::{Ctx, Set, Src, lossy, site},
};
use crate::filter::{ExtType, ExtValue, Filter, FilterLogic};

/// Where tab-delimited annotations come from.
pub(super) enum TabSource {
    /// A tabix-indexed file streamed along with the records (bcftools' `tgts`).
    Tabix(Box<Regions>),
    /// All lines in memory, for BEG,END columns (bcftools' `tgt_idx`).
    Idx(RegIdx<Vec<u8>>),
}

/// The tab-delimited annotation state.
pub(super) struct TabAnnot {
    pub(super) source: TabSource,
    pub(super) cols: TabColumns,
    pub(super) match_id: i32,
    pub(super) match_end: i32,
    pub(super) has_append_mode: bool,
    pub(super) aline_missing: AnnotLine,
    pub(super) min_overlap_ann: f32,
    pub(super) min_overlap_vcf: f32,
    /// Buffered lines (bcftools' `alines`): the first `nalines` are in use, though those
    /// ending before the current record are dead until the next compaction.
    alines: Vec<AnnotLine>,
    nalines: usize,
    /// The buffer size at which dead lines are next dropped.
    compact_at: usize,
    /// Matching lines as `ialt << 16 | iline` (bcftools' `srt_alines`).
    srt: Vec<u32>,
    match_end_str: Vec<u8>,
    ibuf: Vec<i32>,
    /// The record's sequence name, for region index queries.
    chrom: Vec<u8>,
}

/// An external value of a `-i/-e` expression: a column of the annotation line (bcftools'
/// `ext_t`).
#[derive(Debug)]
pub(super) struct Ext {
    pub(super) ty: Option<ExtType>,
    pub(super) icol: i32,
}

/// A `-i/-e` expression with `{}` values taken from annotation lines (bcftools'
/// `filter_ext`).
pub(super) struct ExtFilter {
    pub(super) filter: Filter,
    pub(super) ext: Vec<Ext>,
    pub(super) logic: FilterLogic,
}

impl ExtFilter {
    /// bcftools' `pass_filter_test_ext`.
    pub(super) fn pass(
        &mut self,
        hdr: &Header,
        line: &mut Record,
        ann: &AnnotLine,
    ) -> anyhow::Result<bool> {
        let mut vals: Vec<ExtValue<'_>> = Vec::with_capacity(self.ext.len());
        for e in &self.ext {
            let col = ann.col(e.icol);
            vals.push(match e.ty {
                Some(ExtType::Str) => ExtValue::Str(col),
                Some(ExtType::Int) => {
                    let (v, end) = strtol(col, 0);
                    if end != col.len() {
                        if col != b"." {
                            bail!(
                                "Error: could not parse the annotation file, expected an integer, found \"{}\"",
                                lossy(col)
                            );
                        }
                        ExtValue::Missing
                    } else {
                        ExtValue::Int(v as i32)
                    }
                }
                Some(ExtType::Float) => {
                    let (v, end) = strtod(col, 0);
                    if end != col.len() {
                        if col != b"." {
                            bail!(
                                "Error: could not parse the annotation file, expected a float, found \"{}\"",
                                lossy(col)
                            );
                        }
                        ExtValue::Missing
                    } else {
                        ExtValue::Float(v as f32)
                    }
                }
                None => ExtValue::Missing,
            });
        }
        let pass = self.filter.test_ext(hdr, line, &vals)?;
        Ok(if self.logic == FilterLogic::Exclude {
            !pass
        } else {
            pass
        })
    }
}

/// Whether one of the `;`-separated IDs of `a` is one of those of `b` (bcftools'
/// `strstr_match`).
fn ids_match(a: &[u8], b: &[u8]) -> bool {
    a.split(|&c| c == b';').any(|needle| {
        if needle.is_empty() {
            // strstr finds an empty needle only between two separators.
            b.windows(2).any(|w| w == b";;")
        } else {
            b.split(|&c| c == b';').any(|h| h == needle)
        }
    })
}

/// The arguments the annotation loops share.
pub(super) struct Annotate<'a> {
    pub(super) ctx: &'a mut Ctx,
    pub(super) cols: &'a mut [Column],
    pub(super) hdr: &'a Header,
    pub(super) ext: Option<&'a mut ExtFilter>,
    pub(super) keep_sites: bool,
    pub(super) mark_sites: bool,
}

/// Moves the lines of `alines[..*nalines]` that end at or after `pos` to the front, in order,
/// keeping the buffers of the others for reuse.
fn compact(alines: &mut [AnnotLine], nalines: &mut usize, pos: i64) {
    let mut live = 0;
    for i in 0..*nalines {
        if pos <= alines[i].end {
            if live != i {
                alines.swap(live, i);
            }
            live += 1;
        }
    }
    *nalines = live;
}

/// The buffer size below which dead annotation lines are not dropped.
const MIN_COMPACT: usize = 64;

/// Runs a setter, failing as bcftools' `fixme` message does.
fn run_setter(
    a: &mut Annotate<'_>,
    j: usize,
    line: &mut Record,
    src: Src<'_>,
    what: &str,
) -> anyhow::Result<Set> {
    let r = a.ctx.set(&mut a.cols[j], a.hdr, line, src)?;
    if r == Set::Failed {
        bail!(
            "fixme: Could not set {what}{} at {}",
            lossy(&a.cols[j].key_src),
            site(a.hdr, line)
        );
    }
    Ok(r)
}

impl TabAnnot {
    pub(super) fn new(source: TabSource, cols: TabColumns) -> Self {
        Self {
            source,
            cols,
            match_id: -1,
            match_end: -1,
            has_append_mode: false,
            aline_missing: AnnotLine::default(),
            min_overlap_ann: 0.0,
            min_overlap_vcf: 0.0,
            alines: vec![AnnotLine::default()],
            nalines: 0,
            compact_at: MIN_COMPACT,
            srt: Vec::new(),
            match_end_str: Vec::new(),
            ibuf: Vec::new(),
            chrom: Vec::new(),
        }
    }

    /// Annotates a record; returns whether an annotation line overlapped it.
    pub(super) fn annotate(
        &mut self,
        a: &mut Annotate<'_>,
        line: &mut Record,
    ) -> anyhow::Result<bool> {
        match self.source {
            TabSource::Idx(_) => self.annotate_from_regidx(a, line),
            TabSource::Tabix(_) => self.annotate_from_tab(a, line),
        }
    }

    /// bcftools' `annotate_from_regidx`.
    fn annotate_from_regidx(
        &mut self,
        a: &mut Annotate<'_>,
        line: &mut Record,
    ) -> anyhow::Result<bool> {
        let mut has_overlap = false;
        for c in a.cols.iter_mut() {
            c.done = false;
        }
        let TabSource::Idx(idx) = &mut self.source else {
            unreachable!("called for a region index");
        };
        let pos = line.pos();
        let vcf_end = pos + line.rlen() - 1;
        let clamp = |v: i64| v.clamp(0, i64::from(u32::MAX)) as u32;
        self.chrom.clear();
        self.chrom.extend_from_slice(line.chrom_safe(a.hdr));
        let tmp = &mut self.alines[0];
        for reg in idx.overlap(&self.chrom, clamp(pos), clamp(vcf_end)) {
            tmp.rid = line.rid();
            tmp.start = i64::from(reg.beg);
            tmp.end = i64::from(reg.end);
            let len_ann = (tmp.end - tmp.start + 1) as f32;
            let len_vcf = line.rlen() as f32;
            let isec = (tmp.end.min(vcf_end) - tmp.start.max(pos) + 1) as f32;
            if self.min_overlap_ann != 0.0 && self.min_overlap_ann > isec / len_ann {
                continue;
            }
            if self.min_overlap_vcf != 0.0 && self.min_overlap_vcf > isec / len_vcf {
                continue;
            }
            tmp.parse(reg.payload, &self.cols)?;
            if let Some(ext) = a.ext.as_deref_mut() {
                if !ext.pass(a.hdr, line, tmp)? {
                    continue;
                }
                has_overlap = true;
            }
            if a.cols.is_empty() && a.mark_sites {
                has_overlap = true;
            }
            for j in 0..a.cols.len() {
                if a.cols[j].done || a.cols[j].setter.is_none() {
                    continue;
                }
                if run_setter(a, j, line, Src::Tab(tmp), "")? == Set::Done {
                    a.cols[j].done = true;
                }
                has_overlap = true;
            }
        }
        for j in 0..a.cols.len() {
            let c = &a.cols[j];
            if c.done || c.merge == MergeMethod::First || c.setter.is_none() {
                continue;
            }
            run_setter(a, j, line, Src::Flush, "")?;
        }
        Ok(has_overlap)
    }

    /// bcftools' `buffer_annot_lines`.
    fn buffer_annot_lines(
        &mut self,
        a: &mut Annotate<'_>,
        line: &Record,
        start: i64,
        end: i64,
    ) -> anyhow::Result<()> {
        if self.nalines > 0 && self.alines[0].rid != line.rid() {
            self.nalines = 0;
        }
        // Lines that end before the record are dead. The annotation loop skips them, so they
        // are dropped only when the buffer has doubled since the last time (amortised
        // constant time per line; bcftools moves the buffer's tail for every dropped line,
        // which is quadratic after a long record).
        let pos = line.pos();
        if self.nalines >= self.compact_at {
            compact(&mut self.alines, &mut self.nalines, pos);
            self.compact_at = (2 * self.nalines).max(MIN_COMPACT);
        }
        if a.ext.is_none()
            && self.cols.ref_ == -1
            && self.alines[..self.nalines].iter().any(|l| pos <= l.end)
        {
            return Ok(());
        }
        let TabSource::Tabix(tgts) = &mut self.source else {
            unreachable!("called for a tabix source");
        };
        let chrom = line.chrom_safe(a.hdr);
        while found_or_fail(tgts.overlaps_line(chrom, start, end))? {
            if self.nalines + 1 == 0xffff {
                // Line indexes must fit 16 bits: drop the dead lines, then cap the live ones.
                compact(&mut self.alines, &mut self.nalines, pos);
                if self.nalines + 1 == 0xffff {
                    break;
                }
            }
            self.nalines += 1;
            if self.alines.len() < self.nalines {
                self.alines.push(AnnotLine::default());
            }
            let tmp = &mut self.alines[self.nalines - 1];
            let (iseq, s, e) = tgts.current().expect("an overlapping region is current");
            tmp.rid = line.rid();
            tmp.start = s;
            tmp.end = e;
            tmp.parse(tgts.current_line().unwrap_or_default(), &self.cols)?;
            if a.ext.is_some() || self.cols.ref_ != -1 {
                if !found_or_fail(tgts.next_region())? || tgts.current().map(|c| c.0) != Some(iseq)
                {
                    break;
                }
            } else {
                break;
            }
        }
        Ok(())
    }

    /// bcftools' `annotate_from_tab`.
    #[expect(
        clippy::too_many_lines,
        reason = "a line-for-line port of annotate_from_tab"
    )]
    fn annotate_from_tab(
        &mut self,
        a: &mut Annotate<'_>,
        line: &mut Record,
    ) -> anyhow::Result<bool> {
        let mut has_overlap = false;
        line.unpack(Unpack::STR);
        let end_pos = line.pos() + line.rlen() - 1;
        self.buffer_annot_lines(a, line, line.pos(), end_pos)?;
        self.srt.clear();
        if self.nalines >= 0xffff || line.allele_count() >= 0xffff {
            bail!(
                "Error: too many alleles or annotation lines in the buffer at {} (todo:skip?)",
                site(a.hdr, line)
            );
        }
        self.match_end_str.clear();
        if self.match_end >= 0 && line.get_info_i32(a.hdr, b"END", &mut self.ibuf).ok() == Some(1) {
            pasta::variant::number::push_int(&mut self.match_end_str, i64::from(self.ibuf[0]));
        }
        // The buffered lines are in file order, sorted by start (tabix requires it), so only
        // a slice of them can match: lines starting at the record with REF,ALT, else lines
        // starting at or before its end. Lines ending before it are skipped.
        let buffered = &self.alines[..self.nalines];
        let (lo, hi) = if self.cols.ref_ != -1 {
            let pos = line.pos();
            (
                buffered.partition_point(|l| l.start < pos),
                buffered.partition_point(|l| l.start <= pos),
            )
        } else {
            (0, buffered.partition_point(|l| l.start <= end_pos))
        };
        for i in lo..hi {
            let al = &self.alines[i];
            if line.pos() > al.end || end_pos < al.start {
                continue;
            }
            if self.cols.ref_ != -1 {
                if line.pos() != al.start || al.nals() == 0 {
                    continue;
                }
                let (ann_ref, line_ref) =
                    (super::vcmp::AlleleList::allele_at(al, 0), line.allele(0));
                if !a.ctx.vcmp.set_ref(line_ref, ann_ref) {
                    continue;
                }
                for j in 1..al.nals() {
                    let alt = super::vcmp::AlleleList::allele_at(al, j);
                    let ialt = if line.allele_count() == 1 && alt == b"." {
                        0
                    } else {
                        match a.ctx.vcmp.find_allele(&*line, 1, alt) {
                            Some(k) => k + 1,
                            None => continue,
                        }
                    };
                    if self.match_id >= 0 && !ids_match(line.id(), al.col(self.match_id)) {
                        continue;
                    }
                    if self.match_end >= 0
                        && !self.match_end_str.is_empty()
                        && self.match_end_str != al.col(self.match_end)
                    {
                        continue;
                    }
                    if let Some(ext) = a.ext.as_deref_mut()
                        && !ext.pass(a.hdr, line, al)?
                    {
                        continue;
                    }
                    self.srt.push(((ialt as u32) << 16) | i as u32);
                    has_overlap = true;
                    break;
                }
            } else if let Some(ext) = a.ext.as_deref_mut() {
                if ext.pass(a.hdr, line, al)? {
                    self.srt.push((0xffff << 16) | i as u32);
                    has_overlap = true;
                }
            } else {
                self.srt.push((0xffff << 16) | i as u32);
                has_overlap = true;
            }
        }
        if !has_overlap && a.ext.is_some() && !a.keep_sites {
            return Ok(has_overlap);
        }
        if self.has_append_mode {
            // Stable: equal keys keep their order, as bcftools' insertion sort does.
            self.srt.sort();
        }
        for c in a.cols.iter_mut() {
            c.done = false;
        }
        let mut ialt_exp: usize = 1;
        for k in 0..self.srt.len() {
            let ialt = (self.srt[k] >> 16) as usize;
            let ilin = (self.srt[k] & 0xffff) as usize;
            if self.has_append_mode {
                if ialt_exp > ialt {
                    continue;
                }
                while ialt_exp < ialt {
                    ialt_exp += 1;
                    self.append_missing(a, line)?;
                }
            }
            for j in 0..a.cols.len() {
                if a.cols[j].done || a.cols[j].setter.is_none() {
                    continue;
                }
                if run_setter(a, j, line, Src::Tab(&self.alines[ilin]), "")? == Set::Done {
                    a.cols[j].done = true;
                }
            }
            ialt_exp = ialt + 1;
        }
        if !self.srt.is_empty() {
            if self.has_append_mode && ialt_exp < line.allele_count() {
                while ialt_exp < line.allele_count() {
                    ialt_exp += 1;
                    self.append_missing(a, line)?;
                }
            }
            for j in 0..a.cols.len() {
                let c = &a.cols[j];
                if c.done || c.merge == MergeMethod::First || c.setter.is_none() {
                    continue;
                }
                run_setter(a, j, line, Src::Flush, "")?;
            }
        }
        Ok(has_overlap)
    }

    /// Appends missing values for an ALT without an annotation line (append-missing mode).
    fn append_missing(&self, a: &mut Annotate<'_>, line: &mut Record) -> anyhow::Result<()> {
        for j in 0..a.cols.len() {
            let c = &a.cols[j];
            if c.merge != MergeMethod::AppendMissing || c.done || c.setter.is_none() {
                continue;
            }
            if run_setter(a, j, line, Src::Tab(&self.aline_missing), "missing ")? == Set::Done {
                a.cols[j].done = true;
            }
        }
        Ok(())
    }
}

/// Whether a line of the annotation file overlaps (or was read): a malformed line was reported,
/// and bcftools goes on without it, as with no overlap. It does the same for a line that cannot
/// be read (a truncated or corrupt BGZF block), which loses the annotations after it: that is
/// an error here.
fn found_or_fail(r: Result<bool, RegionsError>) -> anyhow::Result<bool> {
    match r {
        Ok(found) => Ok(found),
        Err(e) if e.is_read_failure() => bail!("Error: failed to read the annotations: {e}"),
        Err(_) => Ok(false),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn id_lists_match_as_bcftools() {
        assert!(ids_match(b"rs1;rs2", b"rs3;rs2"));
        assert!(!ids_match(b"rs1", b"rs10;rs11"));
        assert!(ids_match(b"rs1", b"rs1"));
    }
}
