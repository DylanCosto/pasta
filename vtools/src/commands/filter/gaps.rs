//! `--SnpGap` and `--IndelGap`: filters on variants close to indels, which hold records back
//! until the next ones are known (bcftools' `buffered_filters`).

use std::collections::VecDeque;

use pasta::variant::{
    Header, Record, Unpack, VariantType,
    utils::{self, AcSource},
};

/// Set on a SNP once it is SnpGap-filtered (bcftools abuses bits of `d.var_type`).
const SNP_GAP_SET: u32 = 1 << 16;
/// Set on an indel in the current cluster.
const INDEL_GAP_SET: u32 = 1 << 15;
/// Set on an indel whose cluster was resolved (bcftools sets it but never reads it).
const INDEL_GAP_FLUSH: u32 = 1 << 14;

/// A record held back, with its variant types and flags.
struct Held {
    rec: Record,
    var_type: u32,
    /// The length of REF on the reference (`d.var[0].n`, set to `rlen`).
    ref_len: i64,
}

/// The gap filters' settings and FILTER IDs.
#[derive(Clone, Copy, Debug)]
pub struct GapFilters {
    /// SNPs within this many bases of a variant of `snp_gap_type` are filtered (0 for off;
    /// bcftools also applies a negative gap).
    pub snp_gap: i64,
    pub snp_gap_type: u32,
    pub snp_gap_id: i32,
    /// Indels within this many bases of each other are filtered, all but one (0 for off).
    pub indel_gap: i64,
    pub indel_gap_id: i32,
    /// With a soft filter, gap-filtered records are written; otherwise dropped.
    pub soft: bool,
}

/// The records held back for the gap filters (bcftools' `rbuf`).
pub struct GapBuffer {
    opts: GapFilters,
    held: VecDeque<Held>,
    spare: Vec<Record>,
    ac: Vec<i32>,
}

impl GapBuffer {
    pub fn new(opts: GapFilters) -> Self {
        Self {
            opts,
            held: VecDeque::with_capacity(64),
            spare: Vec::new(),
            ac: Vec::new(),
        }
    }

    /// An empty record to read the next record into, taking `rec`'s place.
    pub fn spare(&mut self) -> Record {
        self.spare.pop().unwrap_or_default()
    }

    /// Writes (or, without a soft filter, drops if gap-filtered) the first `n` records.
    fn flush(
        &mut self,
        n: usize,
        write: &mut dyn FnMut(&mut Record) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        for _ in 0..n {
            let Some(mut h) = self.held.pop_front() else {
                break;
            };
            // A record whose FILTER was never unpacked (one failing the expression, kept for
            // -S without a soft filter) has no gap filter: bcftools reads no FILTER IDs there.
            let dropped = !self.opts.soft
                && h.rec.unpacked().contains(Unpack::FLT)
                && h.rec.filters().iter().any(|&f| {
                    (self.opts.indel_gap != 0 && f == self.opts.indel_gap_id)
                        || (self.opts.snp_gap != 0 && f == self.opts.snp_gap_id)
                });
            if !dropped {
                write(&mut h.rec)?;
            }
            h.rec.clear();
            self.spare.push(h.rec);
        }
        Ok(())
    }

    /// Adds a record (or, with `None`, finishes): applies the gap filters to what is held and
    /// writes the records that can no longer change.
    pub fn push(
        &mut self,
        header: &Header,
        line: Option<Record>,
        write: &mut dyn FnMut(&mut Record) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let mut var_type = 0;
        let has_line = line.is_some();
        if let Some(mut rec) = line {
            // A new chromosome: flush everything.
            if self.held.back().is_some_and(|h| h.rec.rid() != rec.rid()) {
                self.flush(self.held.len(), write)?;
            }
            var_type = rec.variant_types().bits();
            let all = rec.variant_types_all().bits();
            let ref_len = rec.rlen();
            self.held.push_back(Held {
                rec,
                var_type: all,
                ref_len,
            });
        }

        let mut k_flush = 1;
        if self.opts.indel_gap != 0 {
            k_flush = self.find_indel_cluster(has_line);
            if k_flush > 0 || !has_line {
                self.resolve_indel_cluster(header, k_flush)?;
            }
        }

        if !has_line {
            return self.flush(self.held.len(), write);
        }

        let mut j_flush = 1;
        if self.opts.snp_gap != 0 {
            j_flush = 0;
            let last = self.held.len() - 1;
            let last_from = self.held[last].rec.pos();
            for i in 0..self.held.len() {
                let h = &mut self.held[i];
                let rec_to = h.rec.pos() + h.ref_len - 1;
                if rec_to + self.opts.snp_gap < last_from {
                    j_flush += 1;
                } else if var_type & self.opts.snp_gap_type != 0
                    && h.var_type & VariantType::SNP.bits() != 0
                    && h.var_type & SNP_GAP_SET == 0
                {
                    // This SNP has not been SnpGap-filtered yet.
                    h.var_type |= SNP_GAP_SET;
                    h.rec.add_filter(self.opts.snp_gap_id);
                } else if var_type & VariantType::SNP.bits() != 0
                    && h.var_type & self.opts.snp_gap_type != 0
                {
                    // The line being added is a SNP and needs to be filtered.
                    let line = &mut self.held[last];
                    line.var_type |= SNP_GAP_SET;
                    line.rec.add_filter(self.opts.snp_gap_id);
                    break;
                }
            }
        }
        self.flush(j_flush.min(k_flush), write)
    }

    /// Marks the indels of the leading cluster (indels within `indel_gap` of the previous
    /// one) and returns how many records it spans, or 0 to wait for more while the cluster
    /// may still grow.
    fn find_indel_cluster(&mut self, has_line: bool) -> usize {
        let mut k_flush = 0;
        let mut last_to: Option<i64> = None;
        let mut ended = false;
        for h in &mut self.held {
            let rec_from = h.rec.pos();
            if last_to.is_some_and(|t| t < rec_from) {
                ended = true;
                break;
            }
            k_flush += 1;
            if h.var_type & VariantType::INDEL.bits() == 0 {
                continue;
            }
            h.var_type |= INDEL_GAP_SET;
            last_to = Some(self.opts.indel_gap + h.rec.pos() + h.ref_len - 1);
        }
        if !ended && has_line && last_to.is_some() {
            k_flush = 0;
        }
        k_flush
    }

    /// Filters all but the best indel of the first `k_flush` records: the best QUAL, else the
    /// largest AC, else the first.
    fn resolve_indel_cluster(&mut self, header: &Header, k_flush: usize) -> anyhow::Result<()> {
        let mut max_ac = -1;
        let mut imax_ac: Option<usize> = None;
        let mut imax_qual: Option<usize> = None;
        let mut max_qual = -1.0f32;
        for i in 0..k_flush.min(self.held.len()) {
            let h = &mut self.held[i];
            if h.var_type & INDEL_GAP_SET == 0 {
                continue;
            }
            // A CalcAcError ends the program with HTSlib's status 1 (main).
            let ret = utils::calc_ac(
                header,
                &mut h.rec,
                &mut self.ac,
                AcSource::InfoThenGenotypes,
            )?;
            let ac1 = self.ac.get(1).copied().unwrap_or(0);
            if imax_ac.is_none() || (ret && max_ac < ac1) {
                max_ac = ac1;
                imax_ac = Some(i);
            }
            let qual = h.rec.qual_raw();
            if imax_qual.is_none() || max_qual < qual {
                max_qual = qual;
                imax_qual = Some(i);
            }
        }
        for i in 0..k_flush.min(self.held.len()) {
            let h = &mut self.held[i];
            if h.var_type & INDEL_GAP_SET == 0 {
                continue;
            }
            h.var_type |= INDEL_GAP_FLUSH;
            let best = if max_qual > 0.0 { imax_qual } else { imax_ac };
            if best != Some(i) {
                h.rec.add_filter(self.opts.indel_gap_id);
            }
        }
        Ok(())
    }
}
