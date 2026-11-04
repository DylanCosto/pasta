//! gVCF blocks: records with INFO/END are held until other files' records or the block ends
//! cut them, then written as merged blocks (bcftools merge's `gvcf_set_alleles`,
//! `gvcf_write_block`, `gvcf_flush`, and `gvcf_stage`).

use anyhow::bail;

use super::{
    ExitOne, Merger,
    maux::{Source, info_ints, is_gvcf_block},
    to_u32,
};

/// C's `INT_MAX`, which bcftools uses as "no limit" for positions.
const INT_MAX: i64 = i32::MAX as i64;

/// Which of the cached gVCF blocks [`Merger::gvcf_flush`] writes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Flush {
    /// Those that end before the next record, or all of them on a new chromosome.
    BeforeNextRecord,
    /// All of them, at the end of the input.
    All,
}

impl Merger {
    /// Merges the alleles of the open gVCF blocks (bcftools' `gvcf_set_alleles`).
    fn gvcf_set_alleles(&mut self) -> anyhow::Result<()> {
        self.ma.als.clear();
        for i in 0..self.ma.n {
            if !self.ma.gvcf_active(i) {
                continue;
            }
            let Some(cur) = self.ma.buf[i].cur else {
                continue;
            };
            self.update_alleles(i, cur).map_err(|_| {
                let line = self.inp.line(self.ma.buf[i].src, i, cur);
                anyhow::anyhow!(
                    "Failed to merge alleles at {}:{}",
                    String::from_utf8_lossy(line.chrom_safe(self.inp.header(i))),
                    line.pos() + 1
                )
            })?;
        }
        Ok(())
    }

    /// Writes the open gVCF blocks over `start..=end` (0-based) as one line (bcftools'
    /// `gvcf_write_block`).
    fn gvcf_write_block(&mut self, chr: &[u8], start: i64, end: i64) -> anyhow::Result<()> {
        let n = self.ma.n;
        let mut ref_base = b'N';
        for i in 0..n {
            if !self.ma.gvcf_active(i) {
                continue;
            }
            if ref_base == b'N' && self.inp.gvcf_lines[i].pos() == start {
                ref_base = self.inp.ref_base(Source::Gvcf, i, 0);
            }
            self.inp.gvcf_lines[i].set_pos(start);
        }
        let mut min = INT_MAX;
        for i in 0..n {
            let g = &mut self.ma.gvcf.as_mut().expect("gVCF merging is on")[i];
            if !g.active {
                continue;
            }
            if g.end < start {
                g.active = false;
                self.ma.buf[i].cur = None;
                continue;
            }
            self.inp.gvcf_ref[i] = Some(ref_base);
            min = min.min(g.end);
        }
        if min == INT_MAX {
            self.ma.gvcf_min = 0;
            return Ok(());
        }

        self.gvcf_set_alleles()?;
        self.merge_chrom2qual()?;
        self.merge_filter()?;
        self.merge_info()?;
        self.merge_format()?;

        if self.out.allele(0).first() == Some(&b'N')
            && let Some(fai) = &mut self.fai
        {
            let pos = self.out.pos();
            let mut seq = Vec::new();
            if !fai.fetch(chr, pos, pos, &mut seq)? {
                return Err(ExitOne.into());
            }
            if let Some(&b) = seq.first() {
                let mut als: Vec<Vec<u8>> = self.out.alleles().map(<[u8]>::to_vec).collect();
                als[0][0] = b;
                self.out.update_alleles(&self.out_hdr, &als);
            }
        }

        // The block's END.
        let mut end = end;
        if end > start {
            end += 1;
            let e = i32::try_from(end).unwrap_or(i32::MAX);
            if let Some(id) = self.ids.end {
                let _ = self.out.update_info_i32_id(&self.out_hdr, id, &[e]);
            }
        } else if let Some(id) = self.ids.end {
            let _ = self.out.remove_info_id(&self.out_hdr, id);
        }
        self.trim_star_allele()?;
        self.write_out()?;

        // Close the blocks that end here and find the next block end.
        let mut min = INT_MAX;
        let gvcf_min = self.ma.gvcf_min;
        for i in 0..n {
            let g = &mut self.ma.gvcf.as_mut().expect("gVCF merging is on")[i];
            if !g.active {
                continue;
            }
            if g.end < end {
                g.active = false;
                self.ma.buf[i].cur = None;
                continue;
            }
            if gvcf_min < g.end + 1 && min > g.end + 1 {
                min = g.end + 1;
            }
        }
        self.ma.gvcf_min = if min == INT_MAX { 0 } else { min };
        Ok(())
    }

    /// Writes the gVCF blocks that end before the next record, or all of them with
    /// [`Flush::All`] or on a new chromosome (bcftools' `gvcf_flush`). Hanging blocks are cut
    /// at the next record:
    ///
    /// ```text
    /// A: 1 END=10, B: 3 END=7, C: 3 END=5
    /// out: 1 END=2 (A), 3 END=5 (A B C), 6 END=7 (A B), 8 END=10 (A)
    /// ```
    pub(super) fn gvcf_flush(&mut self, flush: Flush) -> anyhow::Result<()> {
        let Some(chr) = self.ma.chr.take() else {
            return Ok(());
        };
        let result = self.gvcf_flush_chr(&chr, flush);
        self.ma.chr = Some(chr);
        result
    }

    /// [`Self::gvcf_flush`] on the current chromosome.
    fn gvcf_flush_chr(&mut self, chr: &[u8], flush: Flush) -> anyhow::Result<()> {
        let mut flush_until = INT_MAX;
        if flush == Flush::BeforeNextRecord {
            let sr = &self.inp.sr;
            if let Some(i) = (0..self.ma.n).find(|&i| sr.has_line(i)) {
                let line = sr.line(i).expect("the reader has a line");
                if sr.header(i).contig_name(line.rid()) == Some(chr) {
                    flush_until = line.pos();
                }
            }
        }
        let mut start = if self.ma.gvcf_break >= 0 {
            self.ma.gvcf_break + 1
        } else {
            self.ma.pos
        };
        if let Some(regs) = &mut self.regs {
            let (mut rstart, mut rend) = (-1i64, -1i64);
            for r in regs.overlap(chr, to_u32(start), to_u32(flush_until)) {
                if rstart < 0 {
                    rstart = i64::from(r.beg);
                }
                rend = i64::from(r.end);
            }
            if rstart > start {
                start = rstart;
            }
            if rend < flush_until {
                flush_until = rend + 1;
            }
        }
        while self.ma.gvcf_min != 0 && start < flush_until {
            let tmp = self.ma.gvcf_min.min(flush_until);
            if start > tmp - 1 {
                break;
            }
            self.gvcf_write_block(chr, start, tmp - 1)?;
            start = tmp;
        }
        Ok(())
    }

    /// Opens gVCF blocks for the new records with INFO/END and points each reader at its
    /// records or its open block (bcftools' `gvcf_stage`).
    pub(super) fn gvcf_stage(&mut self, pos: i64) -> anyhow::Result<()> {
        self.ma.gvcf_break = -1;
        self.ma.gvcf_min = INT_MAX;
        for i in 0..self.ma.n {
            let g = &mut self.ma.gvcf.as_mut().expect("gVCF merging is on")[i];
            if g.active && g.end < pos {
                g.active = false;
            }
            if g.active {
                // A gVCF block cannot overlap another record.
                let e = g.end + 1;
                self.ma.gvcf_min = self.ma.gvcf_min.min(e);
                let b = &mut self.ma.buf[i];
                b.beg = 0;
                b.end = 1;
                b.cur = Some(0);
                continue;
            }
            let b = &self.ma.buf[i];
            if b.beg == b.end {
                continue;
            }
            let irec = b.beg;
            let line = self.inp.sr.buffer_record(i, irec);
            let hdr = self.inp.sr.header(i);
            let end = if is_gvcf_block(line)
                && info_ints(line, hdr, b"END", &mut self.scratch_i) == Some(1)
            {
                Some(i64::from(self.scratch_i[0]))
            } else {
                None
            };
            let Some(end) = end else {
                self.ma.gvcf_break = line.pos();
                continue;
            };
            if end == line.pos() + 1 {
                // POS and END are the same: a record without END.
                self.ma.gvcf_break = line.pos();
                continue;
            }
            if end <= line.pos() {
                bail!(
                    "Error: Incorrect END at {}:{} .. END={end}",
                    String::from_utf8_lossy(line.chrom_safe(hdr)),
                    line.pos() + 1
                );
            }
            // A new gVCF block: keep the record until its END, leaving an empty record at
            // this position in the reader's buffer.
            let g = &mut self.ma.gvcf.as_mut().expect("gVCF merging is on")[i];
            g.active = true;
            g.end = end - 1;
            let (rid, mpos) = (self.ma.buf[i].rid, self.ma.pos);
            let slot = self.inp.sr.buffer_record_mut(i, irec);
            std::mem::swap(slot, &mut self.inp.gvcf_lines[i]);
            slot.set_rid(rid);
            slot.set_pos(mpos);
            self.inp.gvcf_lines[i].set_pos(pos);
            self.inp.gvcf_ref[i] = None;
            let b = &mut self.ma.buf[i];
            b.src = Source::Gvcf;
            b.beg = 0;
            b.end = 1;
            b.cur = Some(0);
            self.ma.gvcf_min = self.ma.gvcf_min.min(end);
        }
        if self.ma.gvcf_min == INT_MAX {
            self.ma.gvcf_min = 0;
        }
        Ok(())
    }
}
