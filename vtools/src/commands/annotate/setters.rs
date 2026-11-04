//! The column setters: copying values from an annotation line or record into a record
//! (bcftools' `setter_*` and `vcf_setter_*` functions for the site columns and INFO).
//!
//! FORMAT setters are in [`super::format`].

use anyhow::bail;
use pasta::variant::{
    Header, Number, Record, Unpack,
    encoding::{
        INT32_MISSING, INT32_VECTOR_END, float_is_missing, float_is_vector_end, float_missing,
    },
    number::{strtod, strtol},
    record::UpdateError,
};

use super::{
    columns::{Column, Getter, MergeMethod, Replace, Setter},
    line::AnnotLine,
    samples::SampleMap,
    vcmp::{AlleleList, Vcmp},
};

/// bcftools' missing double (`bcf_double_set_missing`).
pub(super) const DOUBLE_MISSING: u64 = 0x7FF0_0000_0000_0001;

/// Where a setter takes its values from (the `data` argument of bcftools' setters).
pub(super) enum Src<'a> {
    /// A tab-delimited annotation line.
    Tab(&'a AnnotLine),
    /// No line: flush values merged with `-l`, or annotate from the record itself.
    Flush,
    /// An annotation VCF record, with its header.
    Vcf(&'a mut Record, &'a Header),
}

/// What a setter did (bcftools' return values 0, 1, and negative).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Set {
    /// The column is done for this record.
    Done,
    /// More annotation lines may add to the column.
    More,
    /// The record could not be updated.
    Failed,
}

fn updated(r: Result<(), UpdateError>) -> Set {
    if r.is_ok() { Set::Done } else { Set::Failed }
}

/// A list of strings reused between records.
#[derive(Debug, Default)]
pub(super) struct StrList {
    items: Vec<Vec<u8>>,
    n: usize,
}

impl StrList {
    pub(super) fn clear(&mut self) {
        self.n = 0;
    }

    pub(super) fn push(&mut self, s: &[u8]) {
        if self.n == self.items.len() {
            self.items.push(Vec::new());
        }
        let item = &mut self.items[self.n];
        item.clear();
        item.extend_from_slice(s);
        self.n += 1;
    }

    pub(super) fn get(&self, i: usize) -> &[u8] {
        &self.items[i]
    }

    pub(super) fn set(&mut self, i: usize, s: &[u8]) {
        self.items[i].clear();
        self.items[i].extend_from_slice(s);
    }

    /// Packs the strings into fixed-width, NUL-padded blocks (the layout of HTSlib's
    /// `bcf_update_format_string`).
    pub(super) fn pack(&self, out: &mut Vec<u8>) {
        let width = self.items[..self.n].iter().map(Vec::len).max().unwrap_or(0);
        out.clear();
        for s in &self.items[..self.n] {
            out.extend_from_slice(s);
            out.resize(out.len() + width - s.len(), 0);
        }
    }
}

/// Buffers reused between records (bcftools' `tmpi`, `tmpf`, `tmps`, ...).
#[derive(Debug, Default)]
pub(super) struct Scratch {
    pub(super) tmpi: Vec<i32>,
    pub(super) tmpi2: Vec<i32>,
    pub(super) tmpi3: Vec<i32>,
    pub(super) tmpf: Vec<f32>,
    pub(super) tmpf2: Vec<f32>,
    pub(super) tmpf3: Vec<f32>,
    /// Values merged with `-l` (bcftools converts them to doubles in place).
    pub(super) tmpd: Vec<f64>,
    pub(super) tmps: Vec<u8>,
    pub(super) tmps2: Vec<u8>,
    pub(super) tmpks: Vec<u8>,
    pub(super) strs: StrList,
    pub(super) strs2: StrList,
    pub(super) src_pld: Vec<u8>,
    pub(super) dst_pld: Vec<u8>,
}

/// The allele columns of `-c` (`REF`, `ALT`) of an annotation VCF.
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct AlleleColumns {
    /// `REF` is listed.
    ref_: bool,
    /// How `ALT` is listed, if it is: without `+` if any of its columns is.
    alt: Option<Replace>,
}

/// The alleles a record takes from the annotation record paired with it, which decide the
/// alleles copied genotypes are mapped to.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum CopiedAlleles {
    /// The record keeps its alleles.
    None,
    /// The record takes the annotation's REF and keeps its ALTs.
    Ref,
    /// The record takes the annotation's ALTs, and maybe its REF.
    Alts,
}

impl AlleleColumns {
    pub(super) fn new(cols: &[Column]) -> Self {
        let mut out = Self::default();
        for col in cols {
            match col.setter {
                Some(Setter::VcfRef) => out.ref_ = true,
                // An `ALT` column replaces the ALTs before or after a `+ALT` one.
                Some(Setter::VcfAlt) if out.alt.is_none_or(|r| r.has(Replace::MISSING)) => {
                    out.alt = Some(col.replace);
                }
                _ => {}
            }
        }
        out
    }

    /// The alleles the record `line`, before any column is set, takes from its annotation
    /// record: [`Ctx::vcf_set_alt`] replaces the ALTs unless `+ALT` finds some.
    pub(super) fn copied(self, line: &Record) -> CopiedAlleles {
        if self
            .alt
            .is_some_and(|r| !r.has(Replace::MISSING) || line.allele_count() <= 1)
        {
            CopiedAlleles::Alts
        } else if self.ref_ {
            CopiedAlleles::Ref
        } else {
            CopiedAlleles::None
        }
    }
}

/// The state the setters share (the rest of bcftools' `args_t`).
pub(super) struct Ctx {
    pub(super) hdr_out: Header,
    pub(super) vcmp: Vcmp,
    pub(super) allele_cols: AlleleColumns,
    /// What the record being annotated takes from its annotation record, for FMT/GT.
    pub(super) copied_alleles: CopiedAlleles,
    pub(super) smpl: SampleMap,
    /// The `-a` file name, for messages.
    pub(super) targets_fname: String,
    pub(super) s: Scratch,
}

/// The `CHROM:POS` of a record, for messages.
pub(super) fn site(hdr: &Header, rec: &Record) -> String {
    format!(
        "{}:{}",
        String::from_utf8_lossy(rec.chrom_safe(hdr)),
        rec.pos() + 1
    )
}

pub(super) fn lossy(s: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(s)
}

/// The byte at `i`, or 0 past the end (C string semantics).
#[inline]
pub(super) fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

fn is_dot(s: &[u8]) -> bool {
    s == b"."
}

/// Percent-encodes the characters of `needle` (bcftools' `escape_string`).
pub(super) fn escape(s: &[u8], needle: &[u8], out: &mut Vec<u8>) {
    out.clear();
    for &c in s {
        if needle.contains(&c) {
            out.extend_from_slice(format!("%{c:X}").as_bytes());
        } else {
            out.push(c);
        }
    }
}

/// Copies field `isrc` of comma-separated `src` into field `idst` of `dst`, unless the
/// source is `.` or the destination is not `.` (bcftools' `copy_string_field`). Errors are
/// -1 (no such source field) and -2 (no such destination field).
fn copy_string_field(src: &[u8], isrc: usize, dst: &mut Vec<u8>, idst: usize) -> i32 {
    let mut ith = 0;
    let mut start_src = 0;
    while ith < isrc && start_src < src.len() {
        if src[start_src] == b',' {
            ith += 1;
        }
        start_src += 1;
    }
    if ith != isrc {
        return -1;
    }
    let mut end_src = start_src;
    while end_src < src.len() && src[end_src] != 0 && src[end_src] != b',' {
        end_src += 1;
    }
    let ncpy = end_src - start_src;
    if ncpy == 1 && src[start_src] == b'.' {
        return 0;
    }
    let mut ith = 0;
    let mut start_dst = 0;
    while ith < idst && start_dst < dst.len() {
        if dst[start_dst] == b',' {
            ith += 1;
        }
        start_dst += 1;
    }
    if ith != idst {
        return -2;
    }
    let mut end_dst = start_dst;
    while end_dst < dst.len() && dst[end_dst] != b',' {
        end_dst += 1;
    }
    if end_dst - start_dst > 1 || at(dst, start_dst) != b'.' {
        return 0;
    }
    dst.splice(start_dst..end_dst, src[start_src..end_src].iter().copied());
    0
}

/// What [`parse_info_values`] does with a missing value (`.`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Missing {
    /// Pass it on as `None` (`--merge-logic TAG:append-missing`, or a `.TAG` column).
    Keep,
    /// Leave it out.
    Skip,
}

impl Missing {
    /// Missing values are kept when the column appends them or carries them over.
    fn of(col: &Column) -> Self {
        if col.merge == MergeMethod::AppendMissing || col.replace.has(Replace::CARRY_OVER_MISSING) {
            Self::Keep
        } else {
            Self::Skip
        }
    }
}

/// Parses comma-separated values as bcftools' INFO setters do (the loops of
/// `setter_info_int` and `setter_info_real`), calling `push` with `None` for a kept missing
/// value. Fails where `parse` finds no value.
fn parse_info_values(
    s: &[u8],
    missing: Missing,
    mut parse: impl FnMut(&[u8], usize) -> (Option<f64>, usize),
    mut push: impl FnMut(Option<f64>),
) -> Result<(), ()> {
    let mut str_ = 0usize;
    let mut end = 0usize;
    while at(s, end) != 0 {
        if at(s, str_) == b'.' && (at(s, str_ + 1) == 0 || at(s, str_ + 1) == b',') {
            if missing == Missing::Keep {
                push(None);
            }
            if at(s, str_ + 1) == 0 {
                end = str_ + 1;
            }
            str_ += 2;
        } else {
            let (v, e) = parse(s, str_);
            end = e;
            if v.is_none() {
                return Err(());
            }
            push(v);
            str_ = end + 1;
        }
    }
    Ok(())
}

/// Parses an integer with `strtol`; it must end at a comma or the end of the values, as
/// `setter_info_int` checks.
fn parse_int(s: &[u8], start: usize) -> (Option<f64>, usize) {
    let (v, e) = strtol(s, start);
    if e == start || (at(s, e) != 0 && at(s, e) != b',') {
        (None, e)
    } else {
        (Some(f64::from(v as i32)), e)
    }
}

/// Parses a float with `strtod`, whatever follows it (`setter_info_real`).
fn parse_float(s: &[u8], start: usize) -> (Option<f64>, usize) {
    let (v, e) = strtod(s, start);
    if e == start {
        (None, e)
    } else {
        (Some(f64::from(v as f32)), e)
    }
}

/// Whether a source count fits a Number=A or Number=R tag (a single missing value always
/// does).
fn ar_count_ok(number: Option<Number>, n: usize, nals: usize, single_missing: bool) -> bool {
    match number {
        Some(Number::A) => n + 1 == nals || (n == 1 && single_missing),
        Some(Number::R) => n == nals || (n == 1 && single_missing),
        _ => true,
    }
}

impl Ctx {
    /// Runs a column's setter (bcftools' `col->setter(args, line, col, data)`).
    pub(super) fn set(
        &mut self,
        col: &mut Column,
        hdr: &Header,
        line: &mut Record,
        src: Src<'_>,
    ) -> anyhow::Result<Set> {
        let Some(setter) = col.setter else {
            return Ok(Set::Done);
        };
        match (setter, src) {
            (Setter::Filter, src) => self.set_filter(col, line, &src),
            (Setter::Pos, Src::Tab(tab)) => Self::set_pos(col, hdr, line, tab),
            (Setter::Id, src) => Self::set_id(col, line, &src),
            (Setter::Qual, src) => Self::set_qual(col, hdr, line, &src),
            (Setter::InfoFlag, src) => self.set_info_flag(col, hdr, line, &src),
            (Setter::InfoInt, src) => self.set_info_int(col, hdr, line, &src),
            (Setter::InfoReal, src) => self.set_info_real(col, hdr, line, &src),
            (Setter::InfoStr, src) => self.set_info_str(col, hdr, line, &src),
            (Setter::FormatInt, src) => self.set_format_int(col, hdr, line, &src),
            (Setter::FormatReal, src) => self.set_format_real(col, hdr, line, &src),
            (Setter::FormatStr, src) => self.set_format_str(col, hdr, line, &src),
            (Setter::VcfInfoStr, src) => self.vcf_set_info_str(col, hdr, line, src),
            (Setter::VcfId, Src::Vcf(rec, ahdr)) => Ok(self.vcf_set_id(col, line, rec, ahdr)),
            (Setter::VcfFilter, Src::Vcf(rec, ahdr)) => {
                Ok(self.vcf_set_filter(col, line, rec, ahdr))
            }
            (Setter::VcfRef, Src::Vcf(rec, _)) => Ok(self.vcf_set_ref(line, rec)),
            (Setter::VcfAlt, Src::Vcf(rec, _)) => Ok(self.vcf_set_alt(col, line, rec)),
            (Setter::VcfQual, Src::Vcf(rec, _)) => Ok(Self::vcf_set_qual(col, line, rec)),
            (Setter::VcfInfoFlag, Src::Vcf(rec, ahdr)) => {
                Ok(self.vcf_set_info_flag(col, line, rec, ahdr))
            }
            (Setter::VcfInfoInt, Src::Vcf(rec, ahdr)) => {
                self.vcf_set_info_int(col, hdr, line, rec, ahdr)
            }
            (Setter::VcfInfoReal, Src::Vcf(rec, ahdr)) => {
                self.vcf_set_info_real(col, hdr, line, rec, ahdr)
            }
            (Setter::VcfFormatGt, Src::Vcf(rec, ahdr)) => {
                self.vcf_set_format_gt(col, hdr, line, rec, ahdr)
            }
            (Setter::VcfFormatInt, Src::Vcf(rec, ahdr)) => {
                self.vcf_set_format_int(col, hdr, line, rec, ahdr)
            }
            (Setter::VcfFormatReal, Src::Vcf(rec, ahdr)) => {
                self.vcf_set_format_real(col, hdr, line, rec, ahdr)
            }
            (Setter::VcfFormatStr, Src::Vcf(rec, ahdr)) => {
                self.vcf_set_format_str(col, hdr, line, rec, ahdr)
            }
            (setter, _) => bail!(
                "Error: the --merge-logic option cannot be used with {} (yet?)",
                match setter {
                    Setter::Pos => "POS",
                    _ => "this column",
                }
            ),
        }
    }

    fn set_filter(
        &mut self,
        col: &Column,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with FILTER (yet?)");
        };
        let v = tab.col(col.icol);
        if is_dot(v) {
            if col.replace.clears_with_missing() {
                line.update_filter(&[]);
            }
            return Ok(Set::Done);
        }
        let Some(id) = self.hdr_out.id(v) else {
            bail!(
                "The FILTER \"{}\" is not defined in the header, was the -h option provided?",
                lossy(v)
            );
        };
        let id = id as i32;
        if col.replace.has(Replace::SET_OR_APPEND) {
            return Ok(if line.add_filter(id) {
                Set::More
            } else {
                Set::Done
            });
        }
        if !col.replace.has(Replace::MISSING) {
            line.update_filter(&[id]);
            return Ok(Set::Done);
        }
        line.unpack(Unpack::FLT);
        if line.filters().is_empty() {
            line.update_filter(&[id]);
        }
        Ok(Set::Done)
    }

    fn vcf_set_filter(
        &mut self,
        col: &Column,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> Set {
        rec.unpack(Unpack::FLT);
        line.unpack(Unpack::FLT);
        if rec.filters().is_empty() {
            if col.replace.clears_with_missing() {
                line.update_filter(&[]);
            }
            return Set::Done;
        }
        let ids = &mut self.s.tmpi;
        ids.clear();
        ids.extend(rec.filters().iter().map(|&f| {
            ahdr.id_name(f as u32)
                .and_then(|name| self.hdr_out.id(name))
                .map_or(-1, |i| i as i32)
        }));
        if col.replace.has(Replace::SET_OR_APPEND | Replace::MISSING) {
            if col.replace.has(Replace::MISSING) && !line.filters().is_empty() {
                return Set::Done;
            }
            for &id in ids.iter() {
                line.add_filter(id);
            }
            return Set::Done;
        }
        line.update_filter(ids);
        Set::Done
    }

    fn set_pos(
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        tab: &AnnotLine,
    ) -> anyhow::Result<Set> {
        let v = tab.col(col.icol);
        if is_dot(v) {
            return Ok(Set::Done);
        }
        let (pos, end) = strtol(v, 0);
        if end == 0 {
            bail!(
                "Could not parse -POS at {} .. [{}]",
                site(hdr, line),
                lossy(v)
            );
        }
        line.set_pos(i64::from(pos as i32) - 1);
        Ok(Set::Done)
    }

    fn set_id(col: &Column, line: &mut Record, src: &Src<'_>) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with ID (yet?)");
        };
        if col.replace.has(Replace::MATCH_VALUE) {
            return Ok(Set::Done);
        }
        let v = tab.col(col.icol);
        if is_dot(v) {
            return Ok(Set::Done);
        }
        Ok(Self::apply_id(col, line, v))
    }

    fn apply_id(col: &Column, line: &mut Record, id: &[u8]) -> Set {
        if col.replace.has(Replace::SET_OR_APPEND) {
            line.add_id(id);
            return Set::Done;
        }
        line.unpack(Unpack::STR);
        if !col.replace.has(Replace::MISSING) || line.id().is_empty() || line.id() == b"." {
            line.update_id(Some(id));
        }
        Set::Done
    }

    fn vcf_set_id(
        &mut self,
        col: &mut Column,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> Set {
        if col.replace.has(Replace::MATCH_VALUE) {
            return Set::Done;
        }
        let mut buf = std::mem::take(&mut col.buf);
        let id: &[u8] = if col.getter.is_some() {
            let n = self.get(col, ahdr, rec, None, &mut buf);
            if n <= 0 || (n == 1 && buf.first() == Some(&b'.')) {
                col.buf = buf;
                return Set::Done;
            }
            &buf
        } else {
            rec.unpack(Unpack::STR);
            if rec.id() == b"." {
                col.buf = buf;
                return Set::Done;
            }
            buf.clear();
            buf.extend_from_slice(rec.id());
            &buf
        };
        let r = Self::apply_id(col, line, id);
        col.buf = buf;
        r
    }

    /// Runs a column's getter into `out` (bcftools' `col->getter`), returning its length or a
    /// negative value.
    fn get(
        &self,
        col: &Column,
        ahdr: &Header,
        rec: &mut Record,
        local: Option<&mut Record>,
        out: &mut Vec<u8>,
    ) -> i64 {
        match col.getter {
            Some(Getter::InfoStr) => match rec.get_info_string(ahdr, &col.key_src) {
                Ok(v) => {
                    out.clear();
                    out.extend_from_slice(v);
                    v.len() as i64
                }
                Err(_) => -1,
            },
            Some(Getter::Id) => {
                rec.unpack(Unpack::STR);
                out.clear();
                out.extend(rec.id().iter().map(|&c| if c == b';' { b',' } else { c }));
                out.len() as i64
            }
            Some(Getter::Filter) => filter_to_str(ahdr, rec, out),
            Some(Getter::FilterLocal) => match local {
                Some(line) => filter_to_str(&self.hdr_out, line, out),
                None => filter_to_str(&self.hdr_out, rec, out),
            },
            None => -1,
        }
    }

    fn vcf_set_ref(&mut self, line: &mut Record, rec: &mut Record) -> Set {
        rec.unpack(Unpack::STR);
        line.unpack(Unpack::STR);
        if rec.allele(0) == line.allele(0) {
            return Set::Done;
        }
        let mut als: Vec<Vec<u8>> = Vec::with_capacity(line.allele_count());
        als.push(rec.allele(0).to_vec());
        als.extend(line.alleles().skip(1).map(<[u8]>::to_vec));
        line.update_alleles(&self.hdr_out, &als);
        Set::Done
    }

    fn vcf_set_alt(&mut self, col: &Column, line: &mut Record, rec: &mut Record) -> Set {
        rec.unpack(Unpack::STR);
        line.unpack(Unpack::STR);
        if line.allele_count() > 1 && col.replace.has(Replace::MISSING) {
            return Set::Done;
        }
        if rec.allele_count() == line.allele_count()
            && (1..rec.allele_count()).all(|i| rec.allele(i) == line.allele(i))
        {
            return Set::Done;
        }
        let mut als: Vec<Vec<u8>> = Vec::with_capacity(rec.allele_count());
        als.push(line.allele(0).to_vec());
        als.extend(rec.alleles().skip(1).map(<[u8]>::to_vec));
        line.update_alleles(&self.hdr_out, &als);
        Set::Done
    }

    fn set_qual(
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with QUAL (yet?)");
        };
        let v = tab.col(col.icol);
        if is_dot(v) {
            if col.replace.clears_with_missing() {
                line.set_qual(None);
            }
            return Ok(Set::Done);
        }
        if col.replace.has(Replace::MISSING) && line.qual().is_some() {
            return Ok(Set::Done);
        }
        let (q, end) = strtod(v, 0);
        if end == 0 {
            bail!(
                "Could not parse {} at {} .. [{}]",
                lossy(&col.key_src),
                site(hdr, line),
                lossy(v)
            );
        }
        line.set_qual(Some(q as f32));
        Ok(Set::Done)
    }

    fn vcf_set_qual(col: &Column, line: &mut Record, rec: &Record) -> Set {
        let Some(q) = rec.qual() else {
            if col.replace.clears_with_missing() {
                line.set_qual(None);
            }
            return Set::Done;
        };
        if col.replace.has(Replace::MISSING) && line.qual().is_some() {
            return Set::Done;
        }
        line.set_qual(Some(q));
        Set::Done
    }

    fn set_info_flag(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        let Src::Tab(tab) = src else {
            bail!("Error: the --merge-logic option cannot be used with INFO type=Flag (yet?)");
        };
        let v = tab.col(col.icol);
        let set = match v {
            b"." => {
                if col.replace.clears_with_missing() {
                    return Ok(updated(line.update_info_flag(
                        &self.hdr_out,
                        &col.key_dst,
                        false,
                    )));
                }
                return Ok(Set::Done);
            }
            b"1" => true,
            b"0" => false,
            _ => bail!(
                "Could not parse {} at {} .. [{}]",
                lossy(&col.key_src),
                site(hdr, line),
                lossy(v)
            ),
        };
        Ok(updated(line.update_info_flag(
            &self.hdr_out,
            &col.key_dst,
            set,
        )))
    }

    fn vcf_set_info_flag(
        &mut self,
        col: &Column,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> Set {
        // bcftools passes the getter's status as the flag, so an error sets it too.
        let flag = rec.get_info_flag(ahdr, &col.key_src).unwrap_or(true);
        let _ = line.update_info_flag(&self.hdr_out, &col.key_dst, flag);
        Set::Done
    }

    /// Fills a Number=A/R integer tag from `tmpi` (`ntmpi` values for the alleles `als`),
    /// matching alleles (bcftools' `setter_ARinfo_int32`).
    fn set_ar_info_int<L: AlleleList + ?Sized>(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        als: &L,
        ntmpi: usize,
    ) -> anyhow::Result<Set> {
        let nals = als.count();
        if nals == 0 {
            bail!("Cannot fill Number=R,A tags without --columns ..,REF,ALT,..");
        }
        let single_missing = self.s.tmpi.first() == Some(&INT32_MISSING);
        if !ar_count_ok(col.number, ntmpi, nals, single_missing) {
            bail!(
                "Incorrect number of values ({ntmpi}) for the {} tag at {}",
                lossy(&col.key_src),
                site(hdr, line)
            );
        }
        line.unpack(Unpack::STR);
        let ndst = if col.number == Some(Number::A) {
            line.allele_count() - 1
        } else {
            line.allele_count()
        };
        if self.vcmp.map_ar_values(ndst, als, &*line).is_none() {
            bail!("REF alleles not compatible at {}", site(hdr, line));
        }
        let map = self.vcmp.map();
        let ntmpi2 = line
            .get_info_i32(hdr, &col.key_dst, &mut self.s.tmpi2)
            .unwrap_or(0);
        let dst = &mut self.s.tmpi2;
        if ntmpi2 < ndst {
            dst.resize(ndst, INT32_MISSING);
        }
        for (i, &m) in map.iter().enumerate() {
            if m < 0 {
                if ntmpi2 < ndst {
                    dst[i] = INT32_MISSING;
                }
                continue;
            }
            if ntmpi2 == ndst
                && col.replace.has(Replace::MISSING)
                && dst[i] != INT32_MISSING
                && dst[i] != INT32_VECTOR_END
            {
                continue;
            }
            dst[i] = if (m as usize) < ntmpi {
                self.s.tmpi[m as usize]
            } else {
                INT32_MISSING
            };
        }
        let r = line.update_info_i32(&self.hdr_out, &col.key_dst, &dst[..ndst]);
        Ok(if r.is_ok() { Set::More } else { Set::Failed })
    }

    /// Fills a Number=A/R float tag from `tmpf` (bcftools' `setter_ARinfo_real`).
    fn set_ar_info_real<L: AlleleList + ?Sized>(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        als: &L,
        ntmpf: usize,
    ) -> anyhow::Result<Set> {
        let nals = als.count();
        if nals == 0 {
            bail!("Cannot fill Number=R,A tags without --columns ..,REF,ALT,..");
        }
        let single_missing = self.s.tmpf.first().is_some_and(|&f| float_is_missing(f));
        if !ar_count_ok(col.number, ntmpf, nals, single_missing) {
            bail!(
                "Incorrect number of values ({ntmpf}) for the {} tag at {}",
                lossy(&col.key_src),
                site(hdr, line)
            );
        }
        line.unpack(Unpack::STR);
        let ndst = if col.number == Some(Number::A) {
            line.allele_count() - 1
        } else {
            line.allele_count()
        };
        if self.vcmp.map_ar_values(ndst, als, &*line).is_none() {
            bail!("REF alleles not compatible at {}", site(hdr, line));
        }
        let map = self.vcmp.map();
        let ntmpf2 = line
            .get_info_f32(hdr, &col.key_dst, &mut self.s.tmpf2)
            .unwrap_or(0);
        let dst = &mut self.s.tmpf2;
        if ntmpf2 < ndst {
            dst.resize(ndst, float_missing());
        }
        for (i, &m) in map.iter().enumerate() {
            if m < 0 {
                if ntmpf2 < ndst {
                    dst[i] = float_missing();
                }
                continue;
            }
            if ntmpf2 == ndst
                && col.replace.has(Replace::MISSING)
                && !float_is_missing(dst[i])
                && !float_is_vector_end(dst[i])
            {
                continue;
            }
            dst[i] = if (m as usize) < ntmpf {
                self.s.tmpf[m as usize]
            } else {
                float_missing()
            };
        }
        let r = line.update_info_f32(&self.hdr_out, &col.key_dst, &dst[..ndst]);
        Ok(if r.is_ok() { Set::More } else { Set::Failed })
    }

    /// bcftools' `setter_info_int`, for tab-delimited lines and `-l` flushes.
    fn set_info_int(
        &mut self,
        col: &mut Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        if col.replace.has(Replace::SET_OR_APPEND) {
            col.merge = MergeMethod::Append;
        }
        let tab = match src {
            Src::Tab(t) => Some(*t),
            _ => None,
        };
        if tab.is_none()
            && !matches!(
                col.merge,
                MergeMethod::Sum
                    | MergeMethod::Avg
                    | MergeMethod::Min
                    | MergeMethod::Max
                    | MergeMethod::Append
                    | MergeMethod::AppendMissing
            )
        {
            bail!(
                "Error: at the moment only the sum,avg,min,max,append,append-missing options are supported with --merge-logic for INFO type=Integer"
            );
        }
        let carry = col.replace.has(Replace::CARRY_OVER_MISSING);
        if col.replace.has(Replace::SET_OR_APPEND) && col.mm_dbl.is_empty() {
            let n = line
                .get_info_i32(hdr, &col.key_dst, &mut self.s.tmpi)
                .unwrap_or(0);
            if n > 0 && (self.s.tmpi[0] != INT32_MISSING || carry) {
                col.mm_dbl.clear();
                col.mm_dbl
                    .extend(self.s.tmpi[..n].iter().map(|&v| f64::from(v)));
                col.mm_ndat = 1;
            }
        }
        let ntmpi;
        if let Some(tab) = tab {
            let v = tab.col(col.icol);
            let missing = Missing::of(col);
            if is_dot(v) && missing == Missing::Skip {
                return Ok(Set::More);
            }
            let tmpi = &mut self.s.tmpi;
            tmpi.clear();
            let ok = parse_info_values(v, missing, parse_int, |x| {
                tmpi.push(x.map_or(INT32_MISSING, |x| x as i32));
            });
            if ok.is_err() {
                bail!(
                    "Could not parse {} (Type=Integer) at {} .. [{}]",
                    lossy(&col.key_src),
                    site(hdr, line),
                    lossy(v)
                );
            }
            ntmpi = tmpi.len();
            if col.merge != MergeMethod::First {
                let tmpd = &mut self.s.tmpd;
                tmpd.clear();
                tmpd.extend(self.s.tmpi.iter().map(|&v| f64::from(v)));
                merge_values(col, tmpd)?;
                return Ok(Set::More);
            }
        } else {
            let tmpi = &mut self.s.tmpi;
            tmpi.clear();
            if col.merge == MergeMethod::Avg {
                let n = col.mm_ndat as f64;
                tmpi.extend(col.mm_dbl.iter().map(|&v| (v / n) as i32));
            } else {
                tmpi.extend(col.mm_dbl.iter().map(|&v| v as i32));
            }
            ntmpi = tmpi.len();
            col.mm_dbl.clear();
            col.mm_ndat = 0;
        }
        if col.is_ar() {
            let Some(tab) = tab else {
                bail!("Error: the --merge-logic option cannot be used with INFO tags Number=A,R");
            };
            return self.set_ar_info_int(col, hdr, line, tab, ntmpi);
        }
        if col.replace.has(Replace::MISSING) {
            let n = line
                .get_info_i32(hdr, &col.key_dst, &mut self.s.tmpi2)
                .unwrap_or(0);
            if n > 0 && self.s.tmpi2[0] != INT32_MISSING {
                return Ok(Set::Done);
            }
        }
        Ok(updated(line.update_info_i32(
            &self.hdr_out,
            &col.key_dst,
            &self.s.tmpi[..ntmpi],
        )))
    }

    fn vcf_set_info_int(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> anyhow::Result<Set> {
        let Ok(ntmpi) = rec.get_info_i32(ahdr, &col.key_src, &mut self.s.tmpi) else {
            return Ok(Set::Done);
        };
        if col.is_ar() {
            rec.unpack(Unpack::STR);
            return self.set_ar_info_int(col, hdr, line, &*rec, ntmpi);
        }
        if col.replace.has(Replace::MISSING) {
            let n = line
                .get_info_i32(hdr, &col.key_dst, &mut self.s.tmpi2)
                .unwrap_or(0);
            if n > 0 && self.s.tmpi2[0] != INT32_MISSING {
                return Ok(Set::Done);
            }
        }
        Ok(updated(line.update_info_i32(
            &self.hdr_out,
            &col.key_dst,
            &self.s.tmpi[..ntmpi],
        )))
    }

    /// bcftools' `setter_info_real`.
    fn set_info_real(
        &mut self,
        col: &mut Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        if col.replace.has(Replace::SET_OR_APPEND) {
            col.merge = MergeMethod::Append;
        }
        let tab = match src {
            Src::Tab(t) => Some(*t),
            _ => None,
        };
        if tab.is_none()
            && !matches!(
                col.merge,
                MergeMethod::Sum
                    | MergeMethod::Avg
                    | MergeMethod::Min
                    | MergeMethod::Max
                    | MergeMethod::Append
                    | MergeMethod::AppendMissing
            )
        {
            bail!(
                "Error: at the moment only the sum,avg,min,max,append,append-missing options are supported with --merge-logic for INFO type=Float"
            );
        }
        let carry = col.replace.has(Replace::CARRY_OVER_MISSING);
        let to_dbl = |f: f32| {
            if float_is_missing(f) {
                f64::from_bits(DOUBLE_MISSING)
            } else {
                f64::from(f)
            }
        };
        if col.replace.has(Replace::SET_OR_APPEND) && col.mm_dbl.is_empty() {
            let n = line
                .get_info_f32(hdr, &col.key_dst, &mut self.s.tmpf)
                .unwrap_or(0);
            if n > 0 && (!float_is_missing(self.s.tmpf[0]) || carry) {
                col.mm_dbl.clear();
                col.mm_dbl
                    .extend(self.s.tmpf[..n].iter().map(|&v| to_dbl(v)));
                col.mm_ndat = 1;
            }
        }
        let ntmpf;
        if let Some(tab) = tab {
            let v = tab.col(col.icol);
            let missing = Missing::of(col);
            if is_dot(v) && missing == Missing::Skip {
                return Ok(Set::More);
            }
            let tmpf = &mut self.s.tmpf;
            tmpf.clear();
            let ok = parse_info_values(v, missing, parse_float, |x| {
                tmpf.push(x.map_or_else(float_missing, |x| x as f32));
            });
            if ok.is_err() {
                bail!(
                    "Could not parse {} at {} .. [{}]",
                    lossy(&col.key_src),
                    site(hdr, line),
                    lossy(v)
                );
            }
            ntmpf = tmpf.len();
            if col.merge != MergeMethod::First {
                if col.mm_dbl.is_empty()
                    || matches!(col.merge, MergeMethod::Append | MergeMethod::AppendMissing)
                {
                    col.mm_dbl.extend(self.s.tmpf.iter().map(|&v| to_dbl(v)));
                } else {
                    // Sums and extremes of floats use the raw values, missing ones included.
                    let tmpd = &mut self.s.tmpd;
                    tmpd.clear();
                    tmpd.extend(self.s.tmpf.iter().map(|&v| f64::from(v)));
                    merge_values(col, tmpd)?;
                    return Ok(Set::More);
                }
                col.mm_ndat += 1;
                return Ok(Set::More);
            }
        } else {
            let tmpf = &mut self.s.tmpf;
            tmpf.clear();
            if col.merge == MergeMethod::Avg {
                let n = col.mm_ndat as f64;
                tmpf.extend(col.mm_dbl.iter().map(|&v| (v / n) as f32));
            } else {
                tmpf.extend(col.mm_dbl.iter().map(|&v| {
                    if v.to_bits() == DOUBLE_MISSING {
                        float_missing()
                    } else {
                        v as f32
                    }
                }));
            }
            ntmpf = tmpf.len();
            col.mm_dbl.clear();
            col.mm_ndat = 0;
        }
        if col.is_ar() {
            let Some(tab) = tab else {
                bail!("Error: the --merge-logic option cannot be used with INFO tags Number=A,R");
            };
            return self.set_ar_info_real(col, hdr, line, tab, ntmpf);
        }
        if col.replace.has(Replace::MISSING) {
            let n = line
                .get_info_f32(hdr, &col.key_dst, &mut self.s.tmpf2)
                .unwrap_or(0);
            if n > 0 && !float_is_missing(self.s.tmpf2[0]) {
                return Ok(Set::Done);
            }
        }
        Ok(updated(line.update_info_f32(
            &self.hdr_out,
            &col.key_dst,
            &self.s.tmpf[..ntmpf],
        )))
    }

    fn vcf_set_info_real(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        rec: &mut Record,
        ahdr: &Header,
    ) -> anyhow::Result<Set> {
        let Ok(ntmpf) = rec.get_info_f32(ahdr, &col.key_src, &mut self.s.tmpf) else {
            return Ok(Set::Done);
        };
        if col.is_ar() {
            rec.unpack(Unpack::STR);
            return self.set_ar_info_real(col, hdr, line, &*rec, ntmpf);
        }
        if col.replace.has(Replace::MISSING) {
            let n = line
                .get_info_f32(hdr, &col.key_dst, &mut self.s.tmpf2)
                .unwrap_or(0);
            if n > 0 && !float_is_missing(self.s.tmpf2[0]) {
                return Ok(Set::Done);
            }
        }
        Ok(updated(line.update_info_f32(
            &self.hdr_out,
            &col.key_dst,
            &self.s.tmpf[..ntmpf],
        )))
    }

    /// Fills a Number=A/R string tag from `tmps` (bcftools' `setter_ARinfo_string`).
    fn set_ar_info_str<L: AlleleList + ?Sized>(
        &mut self,
        col: &Column,
        hdr: &Header,
        line: &mut Record,
        als: &L,
    ) -> anyhow::Result<Set> {
        let nals = als.count();
        if nals == 0 {
            bail!("Cannot fill Number=R,A tags without --columns ..,REF,ALT,..");
        }
        let src = std::mem::take(&mut self.s.tmps);
        let nsrc = 1 + src.iter().filter(|&&c| c == b',').count();
        if !ar_count_ok(col.number, nsrc, nals, src == b".") {
            self.s.tmps = src;
            bail!(
                "Incorrect number of values ({nsrc}) for the {} tag at {}",
                lossy(&col.key_src),
                site(hdr, line)
            );
        }
        line.unpack(Unpack::STR);
        let ndst = if col.number == Some(Number::A) {
            line.allele_count() - 1
        } else {
            line.allele_count()
        };
        if self.vcmp.map_ar_values(ndst, als, &*line).is_none() {
            self.s.tmps = src;
            bail!("REF alleles not compatible at {}", site(hdr, line));
        }
        let map = self.vcmp.map();
        let dst = &mut self.s.tmpks;
        dst.clear();
        match line.get_info_string(hdr, &col.key_dst) {
            Ok(v) if v != b"." => dst.extend_from_slice(v),
            _ => {
                dst.push(b'.');
                for _ in 1..ndst {
                    dst.extend_from_slice(b",.");
                }
            }
        }
        // A lone `.`, accepted above for any number of alleles, is missing for each of them,
        // and missing values are not copied. (bcftools looks for the lone `.`'s second value
        // and fails with "Failed to copy a string field".)
        let all_missing = src == b".";
        for (i, &m) in map.iter().enumerate() {
            if m < 0 || all_missing {
                continue;
            }
            if col.replace.has(Replace::MISSING) {
                // Do not replace filled values.
                let mut n = 0;
                let mut p = 0;
                while p < dst.len() && n < i {
                    if dst[p] == b',' {
                        n += 1;
                    }
                    p += 1;
                }
                if at(dst, p) != b'.' || (at(dst, p + 1) != b',' && at(dst, p + 1) != 0) {
                    continue;
                }
            }
            if copy_string_field(&src, m as usize, dst, i) != 0 {
                self.s.tmps = src;
                bail!("[vcfannotate.c:setter_ARinfo_string] Failed to copy a string field");
            }
        }
        self.s.tmps = src;
        Ok(updated(line.update_info_string(
            &self.hdr_out,
            &col.key_dst,
            Some(&self.s.tmpks),
        )))
    }

    /// bcftools' `setter_info_str`.
    fn set_info_str(
        &mut self,
        col: &mut Column,
        hdr: &Header,
        line: &mut Record,
        src: &Src<'_>,
    ) -> anyhow::Result<Set> {
        if col.replace.has(Replace::MISSING)
            && !col.is_ar()
            && let Ok(v) = line.get_info_string(hdr, &col.key_dst)
            && !v.is_empty()
            && v != b"."
        {
            return Ok(Set::Done);
        }
        if col.replace.has(Replace::SET_OR_APPEND) {
            col.merge = MergeMethod::Unique;
        }
        let tab = match src {
            Src::Tab(t) => Some(*t),
            _ => None,
        };
        let mut escaped = std::mem::take(&mut self.s.tmps2);
        if let Some(tab) = tab {
            let v = tab.col(col.icol);
            if v.is_empty() {
                self.s.tmps2 = escaped;
                return Ok(Set::Done);
            }
            if is_dot(v)
                && col.merge != MergeMethod::AppendMissing
                && !col.replace.has(Replace::CARRY_OVER_MISSING)
            {
                self.s.tmps2 = escaped;
                return Ok(Set::More);
            }
            escape(v, b";=", &mut escaped);
        }
        let result = self.set_info_str_escaped(col, hdr, line, tab, &escaped);
        self.s.tmps2 = escaped;
        result
    }

    fn set_info_str_escaped(
        &mut self,
        col: &mut Column,
        hdr: &Header,
        line: &mut Record,
        tab: Option<&AnnotLine>,
        escaped: &[u8],
    ) -> anyhow::Result<Set> {
        if col.merge == MergeMethod::First {
            let Some(tab) = tab else {
                bail!(
                    "Error: the --merge-logic option cannot be used with INFO type=String (yet?)"
                );
            };
            self.s.tmps.clear();
            self.s.tmps.extend_from_slice(escaped);
            if col.is_ar() {
                return self.set_ar_info_str(col, hdr, line, tab);
            }
        } else {
            if col.is_ar() {
                bail!(
                    "Error: the --merge-logic option cannot be used with INFO tags Type=String,Number={{A,R,G}}"
                );
            }
            if tab.is_some() {
                if col.merge == MergeMethod::Unique {
                    if col.mm_seen.contains(escaped) {
                        return Ok(Set::More);
                    }
                    col.mm_seen.insert(escaped.to_vec());
                }
                if col.replace.has(Replace::SET_OR_APPEND)
                    && col.mm_str.is_empty()
                    && let Ok(v) = line.get_info_string(hdr, &col.key_dst)
                    && !v.is_empty()
                    && (col.replace.has(Replace::CARRY_OVER_MISSING) || v != b".")
                {
                    col.mm_str.extend_from_slice(v);
                }
                if !col.mm_str.is_empty() {
                    col.mm_str.push(b',');
                }
                col.mm_str.extend_from_slice(escaped);
                return Ok(Set::More);
            }
            if col.mm_str.is_empty() {
                return Ok(Set::Done);
            }
            self.s.tmps.clear();
            self.s.tmps.extend_from_slice(&col.mm_str);
            if col.merge == MergeMethod::Unique {
                col.mm_seen.clear();
            }
            col.mm_str.clear();
        }
        Ok(updated(line.update_info_string(
            &self.hdr_out,
            &col.key_dst,
            Some(&self.s.tmps),
        )))
    }

    /// bcftools' `vcf_setter_info_str` (also used without `-a` for `INFO/TAG:=FILTER`).
    fn vcf_set_info_str(
        &mut self,
        col: &mut Column,
        hdr: &Header,
        line: &mut Record,
        src: Src<'_>,
    ) -> anyhow::Result<Set> {
        let mut buf = std::mem::take(&mut self.s.tmps);
        let rec = match src {
            Src::Vcf(rec, ahdr) => {
                if col.getter.is_some() {
                    // `INFO/TAG:=./FILTER` reads the annotated record (bcftools'
                    // `args->current_rec`), the other getters the annotation record.
                    self.get(col, ahdr, rec, Some(&mut *line), &mut buf);
                } else {
                    match rec.get_info_string(ahdr, &col.key_src) {
                        Ok(v) => {
                            buf.clear();
                            buf.extend_from_slice(v);
                        }
                        Err(_) => {
                            self.s.tmps = buf;
                            return Ok(Set::Done);
                        }
                    }
                }
                Some(rec)
            }
            _ => {
                if col.getter == Some(Getter::FilterLocal) {
                    filter_to_str(&self.hdr_out, line, &mut buf);
                }
                None
            }
        };
        self.s.tmps = buf;
        if col.is_ar() {
            let Some(rec) = rec else {
                bail!("Error: the --merge-logic option cannot be used with INFO tags Number=A,R");
            };
            rec.unpack(Unpack::STR);
            return self.set_ar_info_str(col, hdr, line, &*rec);
        }
        if col.replace.has(Replace::MISSING)
            && let Ok(v) = line.get_info_string(hdr, &col.key_dst)
            && !v.is_empty()
            && v != b"."
        {
            return Ok(Set::Done);
        }
        Ok(updated(line.update_info_string(
            &self.hdr_out,
            &col.key_dst,
            Some(&self.s.tmps),
        )))
    }
}

/// Merges one annotation line's numbers into a column (the `-l sum/avg/min/max/append` part
/// of `setter_info_int` and `setter_info_real`).
fn merge_values(col: &mut Column, vals: &[f64]) -> anyhow::Result<()> {
    if col.mm_dbl.is_empty()
        || matches!(col.merge, MergeMethod::Append | MergeMethod::AppendMissing)
    {
        col.mm_dbl.extend_from_slice(vals);
    } else {
        if vals.len() != col.mm_dbl.len() {
            bail!("Error: cannot merge fields of unequal length");
        }
        for (m, &v) in col.mm_dbl.iter_mut().zip(vals) {
            match col.merge {
                MergeMethod::Sum | MergeMethod::Avg => *m += v,
                MergeMethod::Min if *m > v => *m = v,
                MergeMethod::Max if *m < v => *m = v,
                _ => {}
            }
        }
    }
    col.mm_ndat += 1;
    Ok(())
}

/// Formats a record's FILTER column as a comma-separated list, `.` when empty (bcftools'
/// `vcf_getter_filter2str_core`).
pub(super) fn filter_to_str(hdr: &Header, rec: &mut Record, out: &mut Vec<u8>) -> i64 {
    rec.unpack(Unpack::FLT);
    out.clear();
    if rec.filters().is_empty() {
        out.push(b'.');
    }
    for (i, &f) in rec.filters().iter().enumerate() {
        if i > 0 {
            out.push(b',');
        }
        out.extend_from_slice(hdr.id_name(f as u32).unwrap_or(b""));
    }
    out.len() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn string_fields_copy_into_missing_slots() {
        let mut dst = b".,.,x".to_vec();
        assert_eq!(copy_string_field(b"a,bb,c", 1, &mut dst, 1), 0);
        assert_eq!(dst, b".,bb,x");
        assert_eq!(copy_string_field(b"a,bb,c", 0, &mut dst, 2), 0);
        assert_eq!(dst, b".,bb,x");
        assert_eq!(copy_string_field(b"a", 3, &mut dst, 0), -1);
    }

    #[test]
    fn info_values_parse_as_bcftools() {
        let mut v = Vec::new();
        let ok = parse_info_values(b"1,.,3", Missing::Skip, parse_int, |x| v.push(x));
        assert!(ok.is_ok());
        assert_eq!(v, [Some(1.0), Some(3.0)]);
        v.clear();
        let ok = parse_info_values(b"1,.,3", Missing::Keep, parse_int, |x| v.push(x));
        assert!(ok.is_ok());
        assert_eq!(v, [Some(1.0), None, Some(3.0)]);
        assert!(parse_info_values(b"1x", Missing::Skip, parse_int, |_| {}).is_err());
        assert!(parse_info_values(b"5,", Missing::Skip, parse_int, |_| {}).is_err());
    }

    #[test]
    fn escapes_use_upper_case_hex() {
        let mut out = Vec::new();
        escape(b"a;b=c", b";=", &mut out);
        assert_eq!(out, b"a%3Bb%3Dc");
    }
}
