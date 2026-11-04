//! The records buffered at the current position and the choice of which to merge into each
//! output line (bcftools merge's `maux_t`, `maux_reset`, `can_merge`, `stage_line`, and
//! `clean_buffer`).

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Record, ValueType, VariantType, encoding::INT32_VECTOR_END,
    synced::SyncedReaders,
};

use super::{Collapse, Merger, alleles::Alleles};

/// The record is merged (`SKIP_DONE`).
pub(super) const SKIP_DONE: u8 = 1;
/// The record is not compatible with this round, it is merged later (`SKIP_DIFF`).
pub(super) const SKIP_DIFF: u8 = 2;

/// Variant types shifted left by one so that REF gets its own bit (bcftools merge's masks).
pub(super) const SNP_MASK: u32 = (VariantType::SNP.bits() | VariantType::MNP.bits()) << 1;
pub(super) const INDEL_MASK: u32 = VariantType::INDEL.bits() << 1;
pub(super) const INS_MASK: u32 = VariantType::INS.bits() << 1;
pub(super) const DEL_MASK: u32 = VariantType::DEL.bits() << 1;
pub(super) const REF_MASK: u32 = 1;
pub(super) const OTHER_MASK: u32 = VariantType::OTHER.bits() << 1;

/// Where a reader's records at this position live (`buffer_t.lines`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(super) enum Source {
    /// The synced reader's buffer.
    #[default]
    Reader,
    /// The reader's cached gVCF block record.
    Gvcf,
}

/// The state of one buffered record (`maux1_t`).
#[derive(Debug, Default)]
pub(super) struct RecState {
    pub(super) skip: u8,
    /// Maps the record's alleles to the merged alleles.
    pub(super) map: Vec<usize>,
    pub(super) als_differ: bool,
    /// The record's variant types, shifted by one for REF.
    pub(super) var_types: u32,
}

/// One reader's records at this position (`buffer_t`).
#[derive(Debug, Default)]
pub(super) struct Buffer {
    pub(super) rid: i32,
    /// The records at this position are `beg..end`.
    pub(super) beg: usize,
    pub(super) end: usize,
    /// The index of the unknown allele (`<*>`, `<NON_REF>`) of the current record, or 0.
    pub(super) unkn_allele: usize,
    /// The record being merged.
    pub(super) cur: Option<usize>,
    pub(super) rec: Vec<RecState>,
    pub(super) src: Source,
    pub(super) var_types: u32,
}

/// A reader's open gVCF block (`gvcf_aux_t`).
#[derive(Debug, Default)]
pub(super) struct GvcfAux {
    /// The block's INFO/END, 0-based.
    pub(super) end: i64,
    pub(super) active: bool,
}

/// The records read but not yet merged, and the gVCF block records (the parts of bcftools'
/// `bcf_srs_t` and `gvcf_aux_t` merge reads records from).
pub(super) struct Inputs {
    pub(super) sr: SyncedReaders,
    /// Each reader's gVCF block record (`gaux[i].line`).
    pub(super) gvcf_lines: Vec<Record>,
    /// A REF first base set on a gVCF block record (bcftools writes it into the record).
    pub(super) gvcf_ref: Vec<Option<u8>>,
}

impl Inputs {
    /// Returns record `j` of reader `i` from `src`.
    #[inline]
    pub(super) fn line(&self, src: Source, i: usize, j: usize) -> &Record {
        match src {
            Source::Reader => self.sr.buffer_record(i, j),
            Source::Gvcf => &self.gvcf_lines[i],
        }
    }

    /// Returns record `j` of reader `i` from `src` for changes.
    #[inline]
    pub(super) fn line_mut(&mut self, src: Source, i: usize, j: usize) -> &mut Record {
        match src {
            Source::Reader => self.sr.buffer_record_mut(i, j),
            Source::Gvcf => &mut self.gvcf_lines[i],
        }
    }

    /// Returns reader `i`'s header.
    #[inline]
    pub(super) fn header(&self, i: usize) -> &Header {
        self.sr.header(i)
    }

    /// Copies a record's alleles into `out`, with the REF base set on a gVCF block record.
    pub(super) fn alleles_into(&self, src: Source, i: usize, j: usize, out: &mut Alleles) {
        let line = self.line(src, i, j);
        out.assign(line.alleles());
        if src == Source::Gvcf
            && let Some(b) = self.gvcf_ref[i]
            && let Some(first) = out.get_mut(0).first_mut()
        {
            *first = b;
        }
    }

    /// Returns the first REF base of a record, with the base set on a gVCF block record.
    pub(super) fn ref_base(&self, src: Source, i: usize, j: usize) -> u8 {
        if src == Source::Gvcf
            && let Some(b) = self.gvcf_ref[i]
        {
            return b;
        }
        self.line(src, i, j).allele(0).first().copied().unwrap_or(0)
    }
}

/// The merge state at the current position (`maux_t`).
#[derive(Debug, Default)]
pub(super) struct Maux {
    pub(super) n: usize,
    pub(super) pos: i64,
    pub(super) var_types: u32,
    /// The type of each merged allele.
    pub(super) als_types: Vec<u32>,
    pub(super) chr: Option<Vec<u8>>,
    /// The merged alleles (some may end up unused).
    pub(super) als: Alleles,
    /// The alleles of the output line.
    pub(super) out_als: Alleles,
    /// How many records have each merged allele.
    pub(super) cnt: Vec<i32>,
    /// The ploidy of each output sample at the current line.
    pub(super) smpl_ploidy: Vec<i32>,
    pub(super) buf: Vec<Buffer>,
    pub(super) gvcf: Option<Vec<GvcfAux>>,
    /// The smallest END+1 of the open gVCF blocks (1-based), or 0 when none is open.
    pub(super) gvcf_min: i64,
    /// The position of the record that breaks the gVCF blocks, or -1.
    pub(super) gvcf_break: i64,
    pub(super) nout_smpl: usize,
    /// A reader's alleles, copied for merging.
    pub(super) line_als: Alleles,
    pub(super) tmp: Vec<u8>,
    pub(super) map_tmp: Vec<usize>,
}

impl Maux {
    /// bcftools' `maux_init`, without the gVCF state that [`Merger::enable_gvcf`] adds.
    pub(super) fn new(sr: &SyncedReaders) -> Self {
        let n = sr.reader_count();
        let nout_smpl = (0..n).map(|i| sr.header(i).sample_count()).sum();
        Self {
            n,
            pos: -1,
            smpl_ploidy: vec![0; nout_smpl],
            buf: (0..n)
                .map(|_| Buffer {
                    rid: -1,
                    ..Buffer::default()
                })
                .collect(),
            nout_smpl,
            ..Self::default()
        }
    }

    /// Returns the current record of reader `i`, if any (`maux_get_line`).
    #[inline]
    pub(super) fn line<'a>(&self, inp: &'a Inputs, i: usize) -> Option<&'a Record> {
        let b = &self.buf[i];
        b.cur.map(|c| inp.line(b.src, i, c))
    }

    /// Returns the map of reader `i`'s current record.
    #[inline]
    pub(super) fn cur_map(&self, i: usize) -> &[usize] {
        let b = &self.buf[i];
        &b.rec[b.cur.expect("the reader has a current record")].map
    }

    /// Whether reader `i` has an open gVCF block.
    #[inline]
    pub(super) fn gvcf_active(&self, i: usize) -> bool {
        self.gvcf.as_ref().is_some_and(|g| g[i].active)
    }

    /// Makes sure reader `i` has room for `size` record states.
    fn expand(&mut self, i: usize, size: usize) {
        let rec = &mut self.buf[i].rec;
        if rec.len() < size {
            rec.resize_with(size, RecState::default);
        }
    }

    /// Finds each reader's records at the next position (bcftools' `maux_reset`).
    pub(super) fn reset(&mut self, sr: &SyncedReaders, rid_tab: &mut [i32]) {
        for i in 0..self.n {
            self.expand(i, sr.buffered_count(i) + 1);
        }
        self.als.clear();
        self.pos = -1;
        let mut chr: Option<&[u8]> = None;
        for i in 0..self.n {
            if let Some(line) = sr.line(i) {
                chr = sr.header(i).contig_name(line.rid());
                self.pos = line.pos();
                break;
            }
        }
        let mut new_chr = false;
        if let Some(c) = chr
            && self.chr.as_deref() != Some(c)
        {
            self.chr = Some(c.to_vec());
            new_chr = true;
        }
        for i in 0..self.n {
            if new_chr {
                let chr = self.chr.as_deref().expect("set above");
                rid_tab[i] = sr.header(i).contig_id(chr).map_or(-1, |r| r as i32);
            }
            let pos = self.pos;
            let b = &mut self.buf[i];
            b.rid = rid_tab[i];
            b.beg = usize::from(!sr.has_line(i));
            let mut j = b.beg;
            while j <= sr.buffered_count(i) {
                b.rec[j].skip = 0;
                b.rec[j].var_types = 0;
                let line = sr.buffer_record(i, j);
                if line.rid() != b.rid || line.pos() != pos {
                    break;
                }
                j += 1;
            }
            b.end = j;
            b.cur = None;
            if b.beg < b.end {
                b.src = Source::Reader;
                if let Some(g) = &mut self.gvcf {
                    // A gVCF block cannot overlap the next record.
                    g[i].active = false;
                }
            }
            if new_chr && let Some(g) = &mut self.gvcf {
                g[i].active = false;
            }
        }
    }
}

/// Whether a record is a gVCF block: it has an END past a one-base REF and only unobserved
/// ALTs (bcftools' `is_gvcf_block`, without its side effect; see [`Merger::enable_gvcf`]).
pub(super) fn is_gvcf_block(line: &Record) -> bool {
    if line.rlen() <= 1 || line.allele(0).len() as i64 == line.rlen() {
        return false;
    }
    line.allele_count() == 1 || line.alleles().skip(1).any(super::alleles::is_unobserved)
}

/// Reads an integer INFO field of an unpacked record (`bcf_get_info_int32`); `None` when the
/// tag is not an integer INFO tag or is absent.
pub(super) fn info_ints(
    line: &Record,
    hdr: &Header,
    tag: &[u8],
    out: &mut Vec<i32>,
) -> Option<usize> {
    let id = hdr.id(tag)?;
    if hdr.id_type(HeaderRecordType::Info, id) != Some(ValueType::Integer) {
        return None;
    }
    let field = line
        .info_fields()
        .iter()
        .find(|f| f.key() == id as i32 && !f.is_removed())?;
    let values = line.info_values(field);
    out.clear();
    for j in 0..field.len() {
        match values.int(j) {
            Some(INT32_VECTOR_END) | None => break,
            Some(v) => out.push(v),
        }
    }
    Some(out.len())
}

impl Merger {
    /// Turns gVCF merging on, as bcftools' `is_gvcf_block` does when it first sees a block.
    pub(super) fn enable_gvcf(&mut self) {
        if self.ma.gvcf.is_none() {
            self.opts.do_gvcf = true;
            self.ma.gvcf = Some((0..self.ma.n).map(|_| GvcfAux::default()).collect());
        }
    }

    /// [`is_gvcf_block`] on record `j` of reader `i`, turning gVCF merging on when it is one.
    pub(super) fn check_gvcf_block(&mut self, src: Source, i: usize, j: usize) -> bool {
        let is = is_gvcf_block(self.inp.line(src, i, j));
        if is {
            self.enable_gvcf();
        }
        is
    }

    /// Whether record `irec` of reader `i` can join the records selected so far given `-m`
    /// (bcftools' `types_compatible`).
    fn types_compatible(&mut self, selected: &mut Selected, i: usize, irec: usize) -> bool {
        let Self {
            inp,
            ma,
            vcmp,
            opts,
            scratch_i,
            ..
        } = self;
        let src = ma.buf[i].src;
        let rec_types = ma.buf[i].rec[irec].var_types;
        let mut end = -1;
        if rec_types & OTHER_MASK != 0 {
            let line = inp.line(src, i, irec);
            if info_ints(line, inp.header(i), b"END", scratch_i) == Some(1) {
                end = i64::from(scratch_i[0]);
            }
        }
        if selected.types == 0 {
            selected.end = end;
            selected.types = rec_types;
            return true;
        }
        let collapse = opts.collapse;
        let line = inp.line_mut(src, i, irec);
        // Open gVCF blocks join every line at their positions. bcftools also compares the
        // allele counts with theirs, so a record with a different count never gets merged:
        // it writes the block line again and again (an endless loop).
        if collapse == Collapse::EXACT
            && !selected.blocks_only
            && ma.als.len() != line.allele_count()
        {
            return false;
        }
        if collapse.contains(Collapse::ANY) {
            return true;
        }
        // REF and gVCF REF with no other alleles can be merged with anything.
        if selected.types & REF_MASK != 0 && selected.types & !REF_MASK == 0 {
            return true;
        }
        if rec_types & REF_MASK != 0 && rec_types & !REF_MASK == 0 {
            return true;
        }
        if collapse != Collapse::empty() && collapse != Collapse::EXACT {
            if collapse.intersects(Collapse::SNPS | Collapse::SNP_INS_DEL)
                && rec_types & SNP_MASK != 0
                && selected.types & SNP_MASK != 0
            {
                return true;
            }
            if collapse.contains(Collapse::INDELS)
                && rec_types & INDEL_MASK != 0
                && selected.types & INDEL_MASK != 0
            {
                return true;
            }
            if collapse.contains(Collapse::SNP_INS_DEL) {
                if rec_types & INS_MASK != 0 && selected.types & INS_MASK != 0 {
                    return true;
                }
                if rec_types & DEL_MASK != 0 && selected.types & DEL_MASK != 0 {
                    return true;
                }
            }
        }
        // -m none or exact matching: are the variants of the same type?
        let (x, y) = (selected.types, rec_types);
        if x & y == 0 {
            return false;
        }
        if x & y != x && x & y != y {
            return false;
        }
        if collapse == Collapse::EXACT && (x & y != x || x & y != y) {
            return false;
        }
        if !vcmp.set_ref(ma.als.get(0), line.allele(0)) {
            return false;
        }
        let maux_alts = &ma.als.as_slice()[1..];
        let n_allele = line.allele_count();
        if collapse == Collapse::EXACT {
            if !(1..n_allele).all(|k| vcmp.find_allele(maux_alts, line.allele(k)).is_some()) {
                return false;
            }
        } else {
            let mut any_present = false;
            for k in 1..n_allele {
                // A gVCF REF (<*> or <NON_REF>) does not count.
                if line.has_variant_type(k, VariantType::empty()) == Some(true) {
                    continue;
                }
                if vcmp.find_allele(maux_alts, line.allele(k)).is_some() {
                    any_present = true;
                    break;
                }
            }
            if !any_present {
                return false;
            }
        }
        if selected.types & OTHER_MASK != 0 && rec_types & OTHER_MASK != 0 && selected.end != end {
            return false;
        }
        true
    }

    /// Merges the alleles of record `irec` of reader `i` into the merged alleles (bcftools'
    /// `maux_update_alleles`).
    #[expect(
        clippy::needless_range_loop,
        reason = "a port of maux_update_alleles, indexing alleles and map together"
    )]
    pub(super) fn update_alleles(&mut self, i: usize, irec: usize) -> anyhow::Result<()> {
        let src = self.ma.buf[i].src;
        let collapse = self.opts.collapse;
        let ma = &mut self.ma;
        self.inp.alleles_into(src, i, irec, &mut ma.line_als);
        let n_allele = ma.line_als.len();
        if ma.als.len() == 0 {
            ma.als
                .assign(ma.line_als.as_slice().iter().map(Vec::as_slice));
            let nals = ma.als.len();
            ma.cnt.clear();
            ma.cnt.resize(nals, 1);
            ma.als_types.resize(nals, 0);
            let rec = &mut ma.buf[i].rec[irec];
            rec.map.clear();
            rec.map.extend(0..nals);
            let line = self.inp.line_mut(src, i, irec);
            for k in 0..nals {
                ma.als_types[k] = allele_type(line, k, collapse);
            }
            return Ok(());
        }
        let mut map = std::mem::take(&mut ma.map_tmp);
        let merged =
            super::alleles::merge_alleles(&mut ma.line_als, &mut map, &mut ma.als, &mut ma.tmp);
        let uppercased = match merged {
            Ok(u) => u,
            Err(_) => {
                ma.map_tmp = map;
                bail!(
                    "Failed to merge alleles at {}:{} in {}",
                    String::from_utf8_lossy(ma.chr.as_deref().unwrap_or_default()),
                    self.inp.line(src, i, irec).pos() + 1,
                    self.inp.sr.file_name(i)
                );
            }
        };
        if uppercased {
            self.uppercase_alleles(src, i, irec);
        }
        let ma = &mut self.ma;
        let nals = ma.als.len();
        ma.cnt.resize(nals, 0);
        ma.als_types.resize(nals, 0);
        let line = self.inp.line_mut(src, i, irec);
        for k in 1..n_allele {
            let ik = map[k];
            ma.als_types[ik] = allele_type(line, k, collapse);
            ma.cnt[ik] += 1;
        }
        ma.cnt[0] += 1;
        let rec = &mut ma.buf[i].rec[irec];
        std::mem::swap(&mut rec.map, &mut map);
        ma.map_tmp = map;
        Ok(())
    }

    /// Writes uppercased alleles back to the record, as bcftools' `merge_alleles` does.
    pub(super) fn uppercase_alleles(&mut self, src: Source, i: usize, irec: usize) {
        let als = self.ma.line_als.as_slice();
        let (sr, gvcf_lines) = (&mut self.inp.sr, &mut self.inp.gvcf_lines);
        // Take the record out of the buffer to update it next to its reader's header.
        let slot = match src {
            Source::Reader => sr.buffer_record_mut(i, irec),
            Source::Gvcf => &mut gvcf_lines[i],
        };
        let mut line = std::mem::replace(slot, Record::new());
        line.update_alleles(sr.header(i), als);
        let slot = match src {
            Source::Reader => sr.buffer_record_mut(i, irec),
            Source::Gvcf => &mut gvcf_lines[i],
        };
        *slot = line;
    }

    /// Selects the records that can be merged at this position given `-m`, counting how often
    /// each allele appears; `false` when all are merged (bcftools' `can_merge`).
    pub(super) fn can_merge(&mut self) -> anyhow::Result<bool> {
        let has_gaux = self.ma.gvcf.is_some();
        self.ma.als.clear();
        self.ma.cnt.clear();
        self.ma.var_types = 0;
        // The record whose ID is matched with -m id.
        let mut id: Option<(Source, usize, usize)> = None;
        let mut ref_base = b'N';
        let mut ntodo = 0;

        for i in 0..self.ma.n {
            self.ma.buf[i].var_types = 0;
            if has_gaux && self.ma.gvcf_active(i) {
                let b = &mut self.ma.buf[i];
                let beg = b.beg;
                b.rec[beg].skip = SKIP_DIFF;
                b.var_types |= REF_MASK;
                b.rec[beg].var_types = REF_MASK;
                self.ma.var_types |= REF_MASK;
                continue;
            }
            let (beg, end, src) = {
                let b = &self.ma.buf[i];
                (b.beg, b.end, b.src)
            };
            if beg < end && ref_base == b'N' {
                ref_base = self.inp.ref_base(src, i, beg);
            }
            for j in beg..end {
                if self.ma.buf[i].rec[j].skip & SKIP_DONE != 0 {
                    continue;
                }
                self.ma.buf[i].rec[j].skip = SKIP_DIFF;
                ntodo += 1;
                if self.opts.merge_by_id && id.is_none() {
                    id = Some((src, i, j));
                    continue;
                }
                if self.ma.buf[i].rec[j].var_types == 0 {
                    let collapse = self.opts.collapse;
                    let line = self.inp.line_mut(src, i, j);
                    let mut var_type = line.variant_types_all().bits();
                    if collapse == Collapse::SNP_INS_DEL {
                        var_type &= !VariantType::INDEL.bits();
                    }
                    var_type = if var_type != 0 {
                        var_type << 1
                    } else {
                        REF_MASK
                    };
                    if self.opts.do_gvcf && self.check_gvcf_block(src, i, j) {
                        var_type |= REF_MASK;
                    }
                    self.ma.buf[i].rec[j].var_types = var_type;
                }
                let t = self.ma.buf[i].rec[j].var_types;
                self.ma.var_types |= t;
                self.ma.buf[i].var_types |= t;
            }
        }
        if ntodo == 0 {
            return Ok(false);
        }

        let mut selected = Selected::default();
        for i in 0..self.ma.n {
            if has_gaux && self.ma.gvcf_active(i) {
                // gVCF records inherited from an upstream block have a stale REF and position.
                self.inp.gvcf_ref[i] = Some(ref_base);
                self.inp.gvcf_lines[i].set_pos(self.ma.pos);
                let beg = self.ma.buf[i].beg;
                self.update_alleles(i, beg)?;
                if selected.types == 0 {
                    selected.blocks_only = true;
                }
                selected.types |= REF_MASK;
                continue;
            }
            let (beg, end, src) = {
                let b = &self.ma.buf[i];
                (b.beg, b.end, b.src)
            };
            for j in beg..end {
                if self.ma.buf[i].rec[j].skip & SKIP_DONE != 0 {
                    continue;
                }
                let line_types = self.ma.buf[i].rec[j].var_types;
                if self.opts.merge_by_id {
                    let same_id = id.is_some_and(|(s0, i0, j0)| {
                        self.inp.line(s0, i0, j0).id() == self.inp.line(src, i, j).id()
                    });
                    if !same_id {
                        continue;
                    }
                } else if !self.types_compatible(&mut selected, i, j) {
                    continue;
                }
                selected.types |= line_types;
                selected.blocks_only = false;
                self.ma.buf[i].rec[j].skip = 0;
                self.update_alleles(i, j)?;
            }
        }
        Ok(true)
    }

    /// Picks at most one record per reader for the next output line, preferring those with the
    /// most frequent allele (bcftools' `stage_line`).
    pub(super) fn stage_line(&mut self) {
        let ma = &self.ma;
        let mut icnt: Option<usize> = None;
        for i in 1..ma.als.len() {
            if ma.als_types[i] & REF_MASK != 0 {
                continue;
            }
            if icnt.is_none_or(|c| ma.cnt[c] < ma.cnt[i]) {
                icnt = Some(i);
            }
        }
        let selected_type = icnt.map_or(REF_MASK, |c| ma.als_types[c]);
        let collapse = self.opts.collapse;

        for i in 0..self.ma.n {
            self.ma.buf[i].cur = None;
            let (beg, end, src) = {
                let b = &self.ma.buf[i];
                (b.beg, b.end, b.src)
            };
            if beg >= end {
                continue;
            }
            let mut found = None;
            for j in beg..end {
                if self.ma.buf[i].rec[j].skip != 0 {
                    let mut is_gvcf = self.ma.gvcf_active(i);
                    if !is_gvcf && self.check_gvcf_block(src, i, j) {
                        is_gvcf = true;
                    }
                    if is_gvcf && self.ma.buf[i].rec[j].skip != 0 && !self.ma.gvcf_active(i) {
                        continue;
                    }
                    if !is_gvcf {
                        continue;
                    }
                }
                if self.opts.merge_by_id {
                    found = Some(j);
                    break;
                }
                let b = &self.ma.buf[i];
                let rt = b.rec[j].var_types;
                if selected_type & b.var_types != 0 && selected_type & rt == 0 {
                    continue;
                }
                if selected_type & b.var_types == 0 && rt & REF_MASK != 0 {
                    found = Some(j);
                    break;
                }
                if selected_type == REF_MASK {
                    found = Some(j);
                    break;
                }
                let n_allele = self.inp.line(src, i, j).allele_count();
                if icnt.is_some_and(|c| b.rec[j].map[..n_allele].contains(&c)) {
                    found = Some(j);
                    break;
                }
            }
            if found.is_none() {
                if collapse == Collapse::empty() || collapse == Collapse::EXACT {
                    continue;
                }
                let b = &self.ma.buf[i];
                let vt = self.ma.var_types;
                found = (beg..end).find(|&j| {
                    if b.rec[j].skip != 0 {
                        return false;
                    }
                    if collapse.contains(Collapse::ANY) {
                        return true;
                    }
                    let lt = b.rec[j].var_types;
                    let snps = collapse.contains(Collapse::SNPS);
                    let indels = collapse.contains(Collapse::INDELS);
                    let sid = collapse.contains(Collapse::SNP_INS_DEL);
                    if (vt & SNP_MASK != 0 && lt & SNP_MASK != 0 && snps)
                        || (vt & INDEL_MASK != 0 && lt & INDEL_MASK != 0 && indels)
                        || (vt & INS_MASK != 0 && lt & INS_MASK != 0 && sid)
                        || (vt & DEL_MASK != 0 && lt & DEL_MASK != 0 && sid)
                    {
                        return true;
                    }
                    if lt & REF_MASK != 0 {
                        (vt & SNP_MASK != 0 && snps)
                            || (vt & INDEL_MASK != 0 && indels)
                            || (vt & INS_MASK != 0 && sid)
                            || (vt & DEL_MASK != 0 && sid)
                            || vt & REF_MASK != 0
                    } else if vt & REF_MASK != 0 {
                        (lt & SNP_MASK != 0 && snps)
                            || (lt & INDEL_MASK != 0 && indels)
                            || (lt & INS_MASK != 0 && sid)
                            || (lt & DEL_MASK != 0 && sid)
                    } else {
                        false
                    }
                });
            }
            self.ma.buf[i].cur = found;
        }

        // When only open gVCF blocks are staged, nothing gets done and bcftools writes the
        // same block line forever; stage the first selected record with them instead.
        let ma = &mut self.ma;
        let blocks_only = (0..ma.n).all(|i| ma.buf[i].cur.is_none() || ma.gvcf_active(i));
        if blocks_only
            && let Some((i, j)) = (0..ma.n).filter(|&i| !ma.gvcf_active(i)).find_map(|i| {
                let b = &ma.buf[i];
                (b.beg..b.end).find(|&j| b.rec[j].skip == 0).map(|j| (i, j))
            })
        {
            ma.buf[i].cur = Some(j);
        }

        // Mark the staged records as done.
        for b in &mut self.ma.buf {
            if let Some(c) = b.cur {
                b.rec[c].skip = SKIP_DONE;
            }
        }
    }

    /// Drops the merged records at this position from the readers' buffers (bcftools'
    /// `clean_buffer`).
    pub(super) fn clean_buffer(&mut self) {
        for ir in 0..self.ma.n {
            let pos = self.ma.pos;
            if let Some(g) = &mut self.ma.gvcf {
                let b = &mut self.ma.buf[ir];
                if g[ir].active {
                    if pos > g[ir].end {
                        g[ir].active = false;
                    } else if b.cur.is_none() {
                        // Re-activate an interrupted gVCF block.
                        b.cur = Some(b.beg);
                    }
                }
                if !g[ir].active {
                    b.cur = None;
                }
            }
            let sr = &mut self.inp.sr;
            let nbuffer = sr.buffered_count(ir);
            if nbuffer == 0 {
                continue;
            }
            let rid = self.ma.buf[ir].rid;
            let at = |j: usize| {
                let l = sr.buffer_record(ir, j);
                l.rid() == rid && l.pos() == pos
            };
            if !at(1) {
                continue;
            }
            let mut b = 2;
            while b <= nbuffer && at(b) {
                b += 1;
            }
            sr.discard_buffered(ir, b - 1);
        }
    }
}

/// The records chosen so far by `can_merge` (`selected_t`).
#[derive(Debug, Default)]
pub(super) struct Selected {
    types: u32,
    /// The END of the first record with a symbolic allele.
    end: i64,
    /// Only open gVCF blocks are selected.
    blocks_only: bool,
}

/// The shifted type of allele `k` of a record (bcftools' `maux_update_alleles`): only
/// `-m snp-ins-del` alone tells insertions from deletions.
fn allele_type(line: &mut Record, k: usize, collapse: Collapse) -> u32 {
    let mut t = line.variant(k).map_or(0, |v| v.ty.bits());
    if collapse == Collapse::SNP_INS_DEL {
        t &= !VariantType::INDEL.bits();
    }
    if t != 0 { t << 1 } else { REF_MASK }
}
