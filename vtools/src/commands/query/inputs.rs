//! The inputs of a query: one file, or several read in step as HTSlib's synced reader does.
//!
//! With several files, bcftools prints only the first file's records, but the other files
//! still decide two things: the chromosome order (every file is read through its index, one
//! chromosome at a time, in the order of the first file's index and then the others') and the
//! order of records that share a position (HTSlib's `bcf_sr_sort` emits matching records of
//! all files together, the ones present in most files first).

use std::collections::{HashMap, VecDeque, hash_map::Entry};

use pasta::variant::{
    Header, Record, Unpack, VariantType,
    io::InputFormat,
    synced::{Regions, SyncedError, SyncedReader},
};

/// HTSlib's `SR_*` variant classes used for pairing (`bcf_sr_sort.c`).
const SR_REF: u8 = 1;
const SR_SNP: u8 = 2;
const SR_INDEL: u8 = 4;
const SR_OTHER: u8 = 8;

/// A record waiting in an input's buffer, with an ID to find it again after sorting.
struct Buffered {
    id: u64,
    seq: usize,
    rec: Record,
}

/// One input file and the records read ahead from it.
struct Lane {
    reader: SyncedReader,
    /// Whether the input is VCF, whose INFO is not parsed under a `max_unpack` without it.
    vcf: bool,
    /// The records at the lane's current position, in file order.
    buf: VecDeque<Buffered>,
    /// The first record past them.
    peek: Option<Buffered>,
    done: bool,
}

/// The records of one position, grouped into rows as `bcf_sr_sort` emits them: each row has,
/// per lane, the ID of the record it holds.
struct Sorted {
    seq: usize,
    pos: i64,
    rows: VecDeque<Vec<Option<u64>>>,
}

/// The input files of a query.
pub struct Inputs {
    lanes: Vec<Lane>,
    /// The index of each region sequence, by name: the files are read one sequence at a time
    /// in this order.
    seq_index: HashMap<Vec<u8>, usize>,
    sorted: Option<Sorted>,
    next_id: u64,
    max_unpack: Unpack,
    /// Records to reuse.
    spare: Vec<Record>,
}

impl Inputs {
    /// Wraps a single input.
    pub fn single(reader: SyncedReader) -> Self {
        Self::new(vec![reader], None)
    }

    /// Reads several inputs in step. All of them were opened with regions; without explicit
    /// ones (`regions` is `None`), each visits the whole sequences of every input's index
    /// (HTSlib's `bcf_sr_add_reader` adding each reader's `tbx_seqnames`/`bcf_hdr_seqnames`).
    pub fn several(mut readers: Vec<SyncedReader>, regions: Option<&Regions>) -> Self {
        if let Some(r) = regions {
            return Self::new(readers, Some(r.seq_names()));
        }
        let names: Vec<Vec<Vec<u8>>> = readers
            .iter()
            .map(|r| r.seq_names().into_iter().map(<[u8]>::to_vec).collect())
            .collect();
        // Each input gets its own list (a list can stream from a file, so it is not cloned).
        let all = || {
            let mut all = Regions::default();
            for seqs in &names {
                for name in seqs {
                    all.add(name, -1, -1);
                }
                all.sort_and_merge();
            }
            all
        };
        for r in &mut readers {
            r.set_regions(all());
        }
        Self::new(readers, Some(all().seq_names()))
    }

    fn new(readers: Vec<SyncedReader>, seq_names: Option<&[Vec<u8>]>) -> Self {
        let seq_index = seq_names.map_or_else(HashMap::new, |names| {
            names
                .iter()
                .enumerate()
                .map(|(i, n)| (n.clone(), i))
                .collect()
        });
        let lanes = readers
            .into_iter()
            .map(|mut reader| Lane {
                vcf: reader.reader_mut().format() == InputFormat::Vcf,
                reader,
                buf: VecDeque::new(),
                peek: None,
                done: false,
            })
            .collect();
        Self {
            lanes,
            seq_index,
            sorted: None,
            next_id: 0,
            max_unpack: Unpack::empty(),
            spare: Vec::new(),
        }
    }

    /// Returns the first input's header.
    pub fn header(&self) -> &Header {
        self.lanes[0].reader.header()
    }

    /// Returns every input's header for changes, the first one first.
    pub fn headers_mut(&mut self) -> impl Iterator<Item = &mut Header> {
        self.lanes.iter_mut().map(|l| l.reader.header_mut())
    }

    /// Limits how much of each VCF line is parsed, in every input (HTSlib's
    /// `bcf_srs_t.max_unpack`).
    pub fn set_max_unpack(&mut self, which: Unpack) {
        self.max_unpack = which;
        for lane in &mut self.lanes {
            lane.reader.set_max_unpack(which);
        }
    }

    /// Takes the error to report after reading everything (the last input's, as HTSlib
    /// keeps one error number).
    pub fn take_pending_error(&mut self) -> Option<SyncedError> {
        self.lanes
            .iter_mut()
            .filter_map(|l| l.reader.take_pending_error())
            .last()
    }

    /// [`Self::next_record`] into a boxed record: with one input, the reader's box is swapped
    /// in instead of the record.
    ///
    /// # Errors
    ///
    /// Returns an error for unparsable or unreadable records of any input.
    pub fn next_record_boxed(&mut self, rec: &mut Box<Record>) -> Result<bool, SyncedError> {
        if self.lanes.len() == 1 {
            rec.set_max_unpack(self.max_unpack);
            return self.lanes[0].reader.next_record_boxed(rec);
        }
        self.next_record(rec)
    }

    /// Reads the next record of the first input (HTSlib's `bcf_sr_next_line` followed by
    /// `bcf_sr_has_line(files, 0)`). Returns `false` at the end.
    ///
    /// # Errors
    ///
    /// Returns an error for unparsable or unreadable records of any input.
    pub fn next_record(&mut self, rec: &mut Record) -> Result<bool, SyncedError> {
        if self.lanes.len() == 1 {
            // The reader keeps the record given in exchange, which must parse no more of
            // the next line than its own records do (an INFO/END outside `max_unpack` does
            // not extend the record into a region).
            rec.set_max_unpack(self.max_unpack);
            return self.lanes[0].reader.next_record(rec);
        }
        loop {
            for i in 0..self.lanes.len() {
                self.fill(i)?;
            }
            let Some(min) = self
                .lanes
                .iter()
                .filter_map(|l| l.buf.front().map(|b| (b.seq, b.rec.pos())))
                .min()
            else {
                return Ok(false);
            };
            let at_min = |l: &Lane| l.buf.front().is_some_and(|b| (b.seq, b.rec.pos()) == min);
            let active: Vec<usize> = (0..self.lanes.len())
                .filter(|&i| at_min(&self.lanes[i]))
                .collect();
            if let [only] = active[..] {
                let mut b = self.lanes[only]
                    .buf
                    .pop_front()
                    .expect("an active input has a record");
                if only == 0 {
                    std::mem::swap(rec, &mut b.rec);
                    self.spare.push(b.rec);
                    return Ok(true);
                }
                self.spare.push(b.rec);
                continue;
            }
            if self.sorted.as_ref().is_none_or(|s| (s.seq, s.pos) != min) {
                let rows = sort_position(&mut self.lanes, &active, self.max_unpack);
                self.sorted = Some(Sorted {
                    seq: min.0,
                    pos: min.1,
                    rows,
                });
            }
            let sorted = self.sorted.as_mut().expect("the position is sorted");
            // HTSlib ends the iteration here (`bcf_sr_sort_next` returns 0).
            let Some(row) = sorted.rows.pop_front() else {
                return Ok(false);
            };
            let mut first = None;
            for (i, id) in row.into_iter().enumerate() {
                let Some(id) = id else { continue };
                let buf = &mut self.lanes[i].buf;
                if let Some(k) = buf.iter().position(|b| b.id == id) {
                    let b = buf.remove(k).expect("the position was just found");
                    if i == 0 {
                        first = Some(b.rec);
                    } else {
                        self.spare.push(b.rec);
                    }
                }
            }
            if let Some(mut r) = first {
                std::mem::swap(rec, &mut r);
                self.spare.push(r);
                return Ok(true);
            }
        }
    }

    /// Reads lane `i`'s records at its next position, when it has none buffered (HTSlib's
    /// `_reader_fill_buffer`).
    fn fill(&mut self, i: usize) -> Result<(), SyncedError> {
        if !self.lanes[i].buf.is_empty() {
            return Ok(());
        }
        if let Some(b) = self.lanes[i].peek.take() {
            self.lanes[i].buf.push_back(b);
        }
        while !self.lanes[i].done {
            let mut rec = self.spare.pop().unwrap_or_default();
            rec.set_max_unpack(self.max_unpack);
            if !self.lanes[i].reader.next_record(&mut rec)? {
                self.spare.push(rec);
                self.lanes[i].done = true;
                break;
            }
            let lane = &self.lanes[i];
            let name = rec.chrom_safe(lane.reader.header());
            let seq = self.seq_index.get(name).copied().unwrap_or(usize::MAX);
            let b = Buffered {
                id: self.next_id,
                seq,
                rec,
            };
            self.next_id += 1;
            let lane = &mut self.lanes[i];
            match lane.buf.front() {
                Some(f) if (f.seq, f.rec.pos()) != (b.seq, b.rec.pos()) => {
                    lane.peek = Some(b);
                    break;
                }
                _ => lane.buf.push_back(b),
            }
        }
        Ok(())
    }
}

/// A distinct variant at the position (HTSlib's `var_t`).
struct Var {
    /// `REF>ALT,REF>ALT2` (with `/END` after symbolic alleles).
    text: Vec<u8>,
    nalt: usize,
    /// The `SR_*` class.
    class: u8,
    /// The lanes and records that have it.
    members: Vec<(usize, u64)>,
    /// The groups of lanes it is in.
    mask: Vec<bool>,
}

/// A set of variants paired across lanes (HTSlib's `varset_t`).
struct VarSet {
    vars: Vec<usize>,
    cnt: usize,
    mask: Vec<bool>,
}

/// Builds the variant text `bcf_sr_sort_set` matches records by.
fn variant_text(rec: &mut Record, end_id: Option<u32>, text: &mut Vec<u8>) {
    text.clear();
    rec.unpack(Unpack::STR);
    let n = rec.allele_count();
    let mut end_pos: Option<i32> = None;
    for i in 1..n {
        if i > 1 {
            text.push(b',');
        }
        text.extend_from_slice(rec.allele(0));
        text.push(b'>');
        text.extend_from_slice(rec.allele(i));
        if rec.allele(i).first() == Some(&b'<') {
            let end = *end_pos.get_or_insert_with(|| {
                let Some(id) = end_id else { return 0 };
                rec.unpack(Unpack::INFO);
                rec.info_fields()
                    .iter()
                    .find(|f| u32::try_from(f.key()) == Ok(id))
                    // Only a unique ID: the value may be truncated.
                    .map_or(0, |f| f.int_value() as i32)
            });
            if end != 0 {
                text.push(b'/');
                text.extend_from_slice(end.to_string().as_bytes());
            }
        }
    }
    if n == 1 {
        text.extend_from_slice(rec.allele(0));
        text.extend_from_slice(b">.");
    }
}

/// The `SR_*` class of a record's variant types.
fn sr_class(types: VariantType) -> u8 {
    if types.is_empty() {
        return SR_REF;
    }
    let mut class = 0;
    if types.intersects(VariantType::SNP | VariantType::MNP) {
        class |= SR_SNP;
    }
    if types.contains(VariantType::INDEL) {
        class |= SR_INDEL;
    }
    if types.contains(VariantType::OTHER) {
        class |= SR_OTHER;
    }
    class
}

/// The comma-separated alleles of `a` are all among `b`'s, ignoring case (HTSlib's
/// `multi_is_exact`).
fn multi_is_exact(a: &Var, b: &Var) -> bool {
    if a.nalt != b.nalt || a.text.len() != b.text.len() {
        return false;
    }
    a.text.split(|&c| c == b',').all(|x| {
        b.text
            .split(|&c| c == b',')
            .any(|y| x.eq_ignore_ascii_case(y))
    })
}

/// Groups the records of the active lanes at one position into output rows (HTSlib's
/// `bcf_sr_sort_set` with `BCF_SR_PAIR_EXACT`, the pairing for `COLLAPSE_NONE`).
fn sort_position(
    lanes: &mut [Lane],
    active: &[usize],
    max_unpack: Unpack,
) -> VecDeque<Vec<Option<u64>>> {
    let mut vars: Vec<Var> = Vec::new();
    let mut var_index: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut groups: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut group_vars: Vec<Vec<usize>> = Vec::new();
    let mut text = Vec::new();
    for &ilane in active {
        let lane = &mut lanes[ilane];
        // VCF lines parsed only up to FILTER have no INFO, so no END (HTSlib's `vcf_parse`
        // stops at `max_unpack`).
        let info_parsed =
            !lane.vcf || max_unpack.is_empty() || max_unpack.intersects(Unpack::INFO | Unpack::FMT);
        let end_id = info_parsed
            .then(|| lane.reader.header().id(b"END"))
            .flatten();
        let mut texts: Vec<Vec<u8>> = Vec::new();
        let mut grp: Vec<usize> = Vec::new();
        for b in &mut lane.buf {
            variant_text(&mut b.rec, end_id, &mut text);
            // Duplicate records of one lane are told apart by a counter, so that they can
            // still pair with another lane's duplicates.
            let mut key = text.clone();
            let mut idx = 0;
            let found = loop {
                match var_index.get(&key) {
                    Some(&v) if vars[v].members.last().map(|m| m.0) == Some(ilane) => {
                        key.truncate(text.len());
                        key.extend_from_slice(idx.to_string().as_bytes());
                        idx += 1;
                    }
                    found => break found.copied(),
                }
            };
            let v = found.unwrap_or_else(|| {
                vars.push(Var {
                    text: text.clone(),
                    nalt: 0,
                    class: 0,
                    members: Vec::new(),
                    mask: Vec::new(),
                });
                var_index.insert(key, vars.len() - 1);
                vars.len() - 1
            });
            let var = &mut vars[v];
            var.nalt = b.rec.allele_count().saturating_sub(1);
            var.class = sr_class(b.rec.variant_types());
            var.members.push((ilane, b.id));
            texts.push(text.clone());
            grp.push(v);
        }
        texts.sort_unstable();
        let key = texts.join(&b';');
        if let Entry::Vacant(e) = groups.entry(key) {
            e.insert(group_vars.len());
            group_vars.push(grp);
        }
    }

    let ngrp = group_vars.len();
    for var in &mut vars {
        var.mask = vec![false; ngrp];
    }
    for (g, gv) in group_vars.iter().enumerate() {
        for &v in gv {
            vars[v].mask[g] = true;
        }
    }
    let mut sets: Vec<VarSet> = vars
        .iter()
        .enumerate()
        .map(|(i, v)| VarSet {
            vars: vec![i],
            cnt: v.members.len(),
            mask: v.mask.clone(),
        })
        .collect();

    let mut rows = VecDeque::new();
    while !sets.is_empty() {
        let mut imax = 0;
        for (i, s) in sets.iter().enumerate().skip(1) {
            if sets[imax].cnt < s.cnt {
                imax = i;
            }
        }
        // With exact pairing, two sets pair only when they are in different groups and
        // share a variant of the same class with the same alleles.
        let ipair = (0..sets.len()).find(|&j| {
            let (a, b) = (&sets[imax], &sets[j]);
            !a.mask.iter().zip(&b.mask).any(|(x, y)| *x && *y)
                && a.vars.iter().any(|&x| {
                    b.vars.iter().any(|&y| {
                        let (vx, vy) = (&vars[x], &vars[y]);
                        vx.class == vy.class && (vx.text == vy.text || multi_is_exact(vx, vy))
                    })
                })
        });
        if let Some(j) = ipair.filter(|&j| j != imax) {
            let (lo, hi) = (imax.min(j), imax.max(j));
            let other = sets.remove(hi);
            let set = &mut sets[lo];
            for (m, o) in set.mask.iter_mut().zip(&other.mask) {
                *m |= *o;
            }
            set.vars.extend(other.vars);
            set.cnt += other.cnt;
            continue;
        }
        let set = sets.remove(imax);
        let mut row = vec![None; lanes.len()];
        for v in set.vars {
            for &(lane, id) in &vars[v].members {
                row[lane] = Some(id);
            }
        }
        rows.push_back(row);
    }
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    fn var(text: &str, nalt: usize) -> Var {
        Var {
            text: text.as_bytes().to_vec(),
            nalt,
            class: SR_SNP,
            members: Vec::new(),
            mask: Vec::new(),
        }
    }

    #[test]
    fn alleles_match_in_any_order_and_case() {
        assert!(multi_is_exact(&var("A>C,A>G", 2), &var("A>g,A>C", 2)));
        assert!(!multi_is_exact(&var("A>C,A>G", 2), &var("A>C,A>T", 2)));
        assert!(!multi_is_exact(&var("A>C", 1), &var("A>C,A>T", 2)));
    }

    #[test]
    fn classes_follow_htslib() {
        assert_eq!(sr_class(VariantType::empty()), SR_REF);
        assert_eq!(sr_class(VariantType::MNP), SR_SNP);
        assert_eq!(
            sr_class(VariantType::SNP | VariantType::INDEL),
            SR_SNP | SR_INDEL
        );
        assert_eq!(sr_class(VariantType::BND), 0);
    }
}
