//! A VCF/BCF record (HTSlib's `bcf1_t`): the BCF-encoded shared (site) and individual
//! (sample) blocks, with lazy unpacking into decoded fields.
//!
//! Records are reused: readers clear and refill them, so steady-state reading allocates
//! nothing. Fields are decoded on demand with [`Record::unpack`] (HTSlib's `bcf_unpack`); edits
//! through the `update_*` methods mark parts dirty, and [`Record::sync`] (HTSlib's `bcf1_sync`)
//! re-encodes them for BCF output exactly as HTSlib does.

mod access;
pub(crate) mod rlen;
mod translate;
mod update;
mod variant_type;

use bitflags::bitflags;

pub use self::{
    access::{Genotype, GetError, TypedValues},
    update::{PassOnEmpty, UpdateError},
    variant_type::{Variant, VariantMatch, VariantType, classify},
};
use crate::variant::{
    Header,
    encoding::{self, BcfType, enc_int1, enc_size, enc_vchar, enc_vint, type_shift},
};

bitflags! {
    /// Which parts of a record to decode (HTSlib's `BCF_UN_*`).
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
    pub struct Unpack: u32 {
        /// ID, REF, and ALT (`BCF_UN_STR`).
        const STR = 1;
        /// FILTER, and everything before it (`BCF_UN_FLT`).
        const FLT = 2;
        /// INFO, and everything before it (`BCF_UN_INFO`).
        const INFO = 4;
        /// All site fields (`BCF_UN_SHR`).
        const SHR = 1 | 2 | 4;
        /// FORMAT and sample data (`BCF_UN_FMT`).
        const FMT = 8;
        /// Everything (`BCF_UN_ALL`).
        const ALL = 1 | 2 | 4 | 8;
    }
}

bitflags! {
    /// Problems found while parsing or reading a record (HTSlib's `BCF_ERR_*`).
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
    pub struct RecordError: u32 {
        /// The contig is not defined in the header (`BCF_ERR_CTG_UNDEF`).
        const CTG_UNDEF = 1;
        /// A tag is not defined in the header (`BCF_ERR_TAG_UNDEF`).
        const TAG_UNDEF = 2;
        /// Wrong number of columns (`BCF_ERR_NCOLS`).
        const NCOLS = 4;
        /// A size limit was exceeded (`BCF_ERR_LIMITS`).
        const LIMITS = 8;
        /// An invalid character (`BCF_ERR_CHAR`).
        const CHAR = 16;
        /// An invalid contig (`BCF_ERR_CTG_INVALID`).
        const CTG_INVALID = 32;
        /// An invalid tag (`BCF_ERR_TAG_INVALID`).
        const TAG_INVALID = 64;
    }
}

impl RecordError {
    /// Describes the errors as HTSlib's `bcf_strerror` does into a buffer of `maxbuffer` bytes.
    #[must_use]
    pub fn describe(self, maxbuffer: usize) -> String {
        const DESCRIPTIONS: [(u32, &str); 7] = [
            (1, "Contig not defined in header"),
            (2, "Tag not defined in header"),
            (4, "Incorrect number of columns"),
            (8, "Limits reached"),
            (16, "Invalid character"),
            (32, "Invalid contig"),
            (64, "Invalid tag"),
        ];
        let mut out = String::new();
        let mut code = self.bits();
        let mut ok = true;
        let add = |out: &mut String, desc: &str| -> bool {
            let rem = maxbuffer - out.len();
            let need = desc.len() + usize::from(rem != maxbuffer);
            if rem > need {
                if rem != maxbuffer {
                    out.push(',');
                }
                out.push_str(desc);
                true
            } else {
                let pos = if rem <= 4 { maxbuffer - 4 } else { out.len() };
                out.truncate(pos);
                out.push_str("...");
                false
            }
        };
        for (bit, desc) in DESCRIPTIONS {
            if code & bit != 0 {
                if !add(&mut out, desc) {
                    ok = false;
                    break;
                }
                code &= !bit;
            }
        }
        if code != 0 && ok {
            add(&mut out, "Unknown error");
        }
        out
    }
}

/// Set in `unpacked` when a value needs 64 bits (HTSlib's `BCF_IS_64BIT`).
pub(crate) const IS_64BIT: u32 = 1 << 30;

pub(crate) const DIRTY_ID: u8 = 1;
pub(crate) const DIRTY_ALS: u8 = 2;
pub(crate) const DIRTY_FLT: u8 = 4;
pub(crate) const DIRTY_INF: u8 = 8;

/// Where a decoded INFO or FORMAT field's bytes live (HTSlib's `vptr`/`p` and `*_free`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) enum FieldData {
    /// Marked for removal (a NULL `vptr`).
    #[default]
    Removed,
    /// In the record's block, starting at this offset (the start of the key encoding).
    InBlock(usize),
    /// In its own buffer (the key encoding, descriptor, and values).
    Owned(Vec<u8>),
}

/// A decoded INFO field (HTSlib's `bcf_info_t`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct InfoField {
    pub(crate) key: i32,
    pub(crate) ty: u8,
    pub(crate) len: i32,
    /// HTSlib's `v1` union: integer values, or float bits in the low 32 bits.
    pub(crate) v1: i64,
    pub(crate) data: FieldData,
    pub(crate) vptr_off: u32,
    pub(crate) vptr_len: u32,
}

impl InfoField {
    /// Returns the tag's dictionary ID.
    #[inline]
    #[must_use]
    pub fn key(&self) -> i32 {
        self.key
    }

    /// Returns the BCF type of the values.
    #[inline]
    #[must_use]
    pub fn bcf_type(&self) -> Option<BcfType> {
        BcfType::from_code(self.ty)
    }

    /// Returns the raw BCF type code of the values.
    #[inline]
    #[must_use]
    pub fn type_code(&self) -> u8 {
        self.ty
    }

    /// Returns the number of values (1 for scalars, 0 for flags).
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.len.max(0) as usize
    }

    /// Returns whether the field has no values (a flag).
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len <= 0
    }

    /// Returns whether the field is marked for removal.
    #[inline]
    #[must_use]
    pub fn is_removed(&self) -> bool {
        self.data == FieldData::Removed
    }

    /// Returns the integer value of a single-value field (HTSlib's `v1.i`).
    #[inline]
    #[must_use]
    pub fn int_value(&self) -> i64 {
        self.v1
    }

    /// Returns the float value of a single-value float field (HTSlib's `v1.f`).
    #[inline]
    #[must_use]
    pub fn float_value(&self) -> f32 {
        f32::from_bits(self.v1 as u32)
    }
}

/// A decoded FORMAT field (HTSlib's `bcf_fmt_t`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FormatField {
    pub(crate) id: i32,
    pub(crate) n: i32,
    pub(crate) size: i32,
    pub(crate) ty: u8,
    pub(crate) data: FieldData,
    pub(crate) p_off: u32,
    pub(crate) p_len: u32,
}

impl FormatField {
    /// Returns the tag's dictionary ID.
    #[inline]
    #[must_use]
    pub fn id(&self) -> i32 {
        self.id
    }

    /// Returns the number of values per sample.
    #[inline]
    #[must_use]
    pub fn n(&self) -> usize {
        self.n.max(0) as usize
    }

    /// Returns the number of bytes per sample.
    #[inline]
    #[must_use]
    pub fn size(&self) -> usize {
        self.size.max(0) as usize
    }

    /// Returns the BCF type of the values.
    #[inline]
    #[must_use]
    pub fn bcf_type(&self) -> Option<BcfType> {
        BcfType::from_code(self.ty)
    }

    /// Returns the raw BCF type code.
    #[inline]
    #[must_use]
    pub fn type_code(&self) -> u8 {
        self.ty
    }

    /// Returns whether the field is marked for removal.
    #[inline]
    #[must_use]
    pub fn is_removed(&self) -> bool {
        self.data == FieldData::Removed
    }
}

/// The decoded parts of a record (HTSlib's `bcf_dec_t`).
#[derive(Clone, Debug, Default)]
pub(crate) struct Decoded {
    /// ID text (`.` when missing).
    pub(crate) id: Vec<u8>,
    /// REF and ALT texts, each followed by a NUL.
    pub(crate) als: Vec<u8>,
    /// Start of each allele in `als`.
    pub(crate) allele_starts: Vec<u32>,
    pub(crate) flt: Vec<i32>,
    pub(crate) info: Vec<InfoField>,
    pub(crate) fmt: Vec<FormatField>,
    pub(crate) var: Vec<Variant>,
    /// Variant types of all alleles, or -1 when not computed.
    pub(crate) var_type: i32,
    pub(crate) shared_dirty: u8,
    pub(crate) indiv_dirty: bool,
    /// Emptied buffers of edited fields, reused by later edits instead of allocating.
    spare: Vec<Vec<u8>>,
}

/// How many emptied field buffers a record keeps for reuse.
const MAX_SPARE_BUFFERS: usize = 32;

impl Decoded {
    /// An empty buffer for an edited field: a spare one when there is one.
    #[inline]
    pub(crate) fn take_buf(&mut self) -> Vec<u8> {
        let mut buf = self.spare.pop().unwrap_or_default();
        buf.clear();
        buf
    }

    /// Keeps a field's own buffer for reuse.
    #[inline]
    pub(crate) fn recycle(&mut self, data: FieldData) {
        if let FieldData::Owned(buf) = data {
            self.put_buf(buf);
        }
    }

    /// Keeps a buffer for reuse.
    #[inline]
    pub(crate) fn put_buf(&mut self, buf: Vec<u8>) {
        keep_spare(&mut self.spare, buf);
    }

    /// Clears the INFO fields, keeping their buffers for reuse.
    fn clear_info(&mut self) {
        for f in &mut self.info {
            if let FieldData::Owned(buf) = &mut f.data {
                keep_spare(&mut self.spare, std::mem::take(buf));
            }
        }
        self.info.clear();
    }

    /// Clears the FORMAT fields, keeping their buffers for reuse.
    fn clear_fmt(&mut self) {
        for f in &mut self.fmt {
            if let FieldData::Owned(buf) = &mut f.data {
                keep_spare(&mut self.spare, std::mem::take(buf));
            }
        }
        self.fmt.clear();
    }
}

/// Adds a buffer to the spare ones, up to [`MAX_SPARE_BUFFERS`].
#[inline]
fn keep_spare(spare: &mut Vec<Vec<u8>>, buf: Vec<u8>) {
    if spare.len() < MAX_SPARE_BUFFERS && buf.capacity() > 0 {
        spare.push(buf);
    }
}

/// One VCF/BCF record (HTSlib's `bcf1_t`).
#[derive(Clone, Debug)]
pub struct Record {
    pub(crate) rid: i32,
    pub(crate) pos: i64,
    pub(crate) rlen: i64,
    pub(crate) qual: f32,
    pub(crate) n_info: u32,
    pub(crate) n_allele: u32,
    pub(crate) n_fmt: u32,
    pub(crate) n_sample: u32,
    pub(crate) shared: Vec<u8>,
    pub(crate) indiv: Vec<u8>,
    pub(crate) d: Decoded,
    pub(crate) max_unpack: Unpack,
    pub(crate) unpacked: u32,
    pub(crate) unpack_size: [usize; 3],
    pub(crate) errcode: RecordError,
}

impl Default for Record {
    fn default() -> Self {
        Self::new()
    }
}

impl Record {
    /// Creates an empty record (HTSlib's `bcf_init`).
    #[must_use]
    pub fn new() -> Self {
        Self {
            rid: 0,
            pos: 0,
            rlen: 0,
            qual: 0.0,
            n_info: 0,
            n_allele: 0,
            n_fmt: 0,
            n_sample: 0,
            shared: Vec::new(),
            indiv: Vec::new(),
            d: Decoded {
                var_type: -1,
                ..Decoded::default()
            },
            max_unpack: Unpack::empty(),
            unpacked: 0,
            unpack_size: [0; 3],
            errcode: RecordError::empty(),
        }
    }

    /// Resets the record for reuse, keeping its buffers (HTSlib's `bcf_clear`).
    pub fn clear(&mut self) {
        self.rid = 0;
        self.pos = 0;
        self.rlen = 0;
        self.unpacked = 0;
        self.qual = encoding::float_missing();
        self.n_info = 0;
        self.n_allele = 0;
        self.n_fmt = 0;
        self.n_sample = 0;
        self.shared.clear();
        self.indiv.clear();
        self.d.var_type = -1;
        self.d.shared_dirty = 0;
        self.d.indiv_dirty = false;
        self.d.flt.clear();
        self.d.clear_info();
        self.d.clear_fmt();
        self.errcode = RecordError::empty();
        self.d.als.clear();
        self.d.allele_starts.clear();
        self.d.id.clear();
    }

    /// Limits how much of a VCF line is parsed (HTSlib's `max_unpack`): fields after the given
    /// level are skipped. Empty parses everything.
    pub fn set_max_unpack(&mut self, which: Unpack) {
        self.max_unpack = which;
    }

    // ---- fixed fields ----

    /// Returns the contig ID (HTSlib's `rid`).
    #[inline]
    #[must_use]
    pub fn rid(&self) -> i32 {
        self.rid
    }

    /// Sets the contig ID.
    #[inline]
    pub fn set_rid(&mut self, rid: i32) {
        self.rid = rid;
    }

    /// Returns the contig name (HTSlib's `bcf_seqname`).
    #[inline]
    #[must_use]
    pub fn chrom<'h>(&self, header: &'h Header) -> Option<&'h [u8]> {
        header.contig_name(self.rid)
    }

    /// Returns the contig name or `(unknown)` (HTSlib's `bcf_seqname_safe`).
    #[inline]
    #[must_use]
    pub fn chrom_safe<'h>(&self, header: &'h Header) -> &'h [u8] {
        self.chrom(header).unwrap_or(b"(unknown)")
    }

    /// Returns the 0-based position (HTSlib's `pos`).
    #[inline]
    #[must_use]
    pub fn pos(&self) -> i64 {
        self.pos
    }

    /// Sets the 0-based position.
    #[inline]
    pub fn set_pos(&mut self, pos: i64) {
        self.pos = pos;
    }

    /// Returns the length on the reference (HTSlib's `rlen`): from END, SVLEN, FORMAT/LEN, or
    /// the REF allele.
    #[inline]
    #[must_use]
    pub fn rlen(&self) -> i64 {
        self.rlen
    }

    /// Sets the reference length.
    #[inline]
    pub fn set_rlen(&mut self, rlen: i64) {
        self.rlen = rlen;
    }

    /// Returns the 0-based end position (exclusive): `pos + rlen`.
    #[inline]
    #[must_use]
    pub fn end(&self) -> i64 {
        self.pos + self.rlen
    }

    /// Returns QUAL, or `None` when missing.
    #[inline]
    #[must_use]
    pub fn qual(&self) -> Option<f32> {
        (!encoding::float_is_missing(self.qual)).then_some(self.qual)
    }

    /// Returns QUAL including the missing sentinel.
    #[inline]
    #[must_use]
    pub fn qual_raw(&self) -> f32 {
        self.qual
    }

    /// Sets QUAL (`None` for missing).
    #[inline]
    pub fn set_qual(&mut self, qual: Option<f32>) {
        self.qual = qual.unwrap_or_else(encoding::float_missing);
    }

    /// Returns the number of alleles, REF included.
    #[inline]
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.n_allele as usize
    }

    /// Returns the number of INFO fields.
    #[inline]
    #[must_use]
    pub fn info_count(&self) -> usize {
        self.n_info as usize
    }

    /// Returns the number of FORMAT fields.
    #[inline]
    #[must_use]
    pub fn format_count(&self) -> usize {
        self.n_fmt as usize
    }

    /// Returns the number of samples.
    #[inline]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.n_sample as usize
    }

    /// Sets the number of samples of a record built from scratch (HTSlib's `n_sample`); the
    /// FORMAT updates set it from the header.
    #[inline]
    pub fn set_sample_count(&mut self, n: usize) {
        self.n_sample = n as u32;
    }

    /// Returns the problems found while parsing or reading (HTSlib's `errcode`).
    #[inline]
    #[must_use]
    pub fn errors(&self) -> RecordError {
        self.errcode
    }

    /// Clears the recorded problems.
    #[inline]
    pub fn clear_errors(&mut self) {
        self.errcode = RecordError::empty();
    }

    /// Returns which parts are decoded.
    #[inline]
    #[must_use]
    pub fn unpacked(&self) -> Unpack {
        Unpack::from_bits_truncate(self.unpacked)
    }

    /// Returns whether a value needs 64 bits, which BCF cannot store.
    #[inline]
    #[must_use]
    pub fn is_64bit(&self) -> bool {
        self.unpacked & IS_64BIT != 0
    }

    /// Clears the flag behind [`is_64bit`](Self::is_64bit), as HTSlib's `bcf_dup` does. The
    /// record is then refused as BCF only for values that BCF's 32-bit fields cannot store.
    #[inline]
    pub fn clear_64bit(&mut self) {
        self.unpacked &= !IS_64BIT;
    }

    /// Returns the encoded site block (HTSlib's `shared`).
    #[inline]
    #[must_use]
    pub fn shared_bytes(&self) -> &[u8] {
        &self.shared
    }

    /// Returns the encoded sample block (HTSlib's `indiv`).
    #[inline]
    #[must_use]
    pub fn indiv_bytes(&self) -> &[u8] {
        &self.indiv
    }

    /// Returns the size of REF and ALT in the encoded site block (HTSlib's `unpack_size[1]`),
    /// known once the alleles are unpacked or the record is synced.
    #[inline]
    #[must_use]
    pub fn encoded_alleles_len(&self) -> usize {
        self.unpack_size[1]
    }

    /// Returns the bytes the record's buffers hold allocated: the encoded blocks, the decoded
    /// fields, and the spare field buffers kept for reuse. A record keeps its buffers when it
    /// is cleared, so this is what keeping it for reuse costs.
    #[must_use]
    pub fn allocated_bytes(&self) -> usize {
        fn bytes<T>(v: &Vec<T>) -> usize {
            v.capacity() * std::mem::size_of::<T>()
        }
        fn owned(data: &FieldData) -> usize {
            match data {
                FieldData::Owned(buf) => buf.capacity(),
                FieldData::Removed | FieldData::InBlock(_) => 0,
            }
        }
        let d = &self.d;
        let info: usize = d.info.iter().map(|f| owned(&f.data)).sum();
        let fmt: usize = d.fmt.iter().map(|f| owned(&f.data)).sum();
        let spare: usize = d.spare.iter().map(Vec::capacity).sum();
        bytes(&self.shared)
            + bytes(&self.indiv)
            + bytes(&d.id)
            + bytes(&d.als)
            + bytes(&d.allele_starts)
            + bytes(&d.flt)
            + bytes(&d.info)
            + bytes(&d.fmt)
            + bytes(&d.var)
            + bytes(&d.spare)
            + info
            + fmt
            + spare
    }

    // ---- decoded fields ----

    /// Returns the ID text, `.` when missing. Needs [`Unpack::STR`].
    ///
    /// # Panics
    ///
    /// Panics if the record has encoded data that is not unpacked.
    #[inline]
    #[must_use]
    pub fn id(&self) -> &[u8] {
        self.assert_unpacked(Unpack::STR);
        &self.d.id
    }

    /// Returns allele `i` (0 is REF). Needs [`Unpack::STR`].
    ///
    /// # Panics
    ///
    /// Panics if `i` is out of range or the alleles are not unpacked.
    #[inline]
    #[must_use]
    pub fn allele(&self, i: usize) -> &[u8] {
        self.assert_unpacked(Unpack::STR);
        let start = self.d.allele_starts[i] as usize;
        let rest = &self.d.als[start..];
        &rest[..rest.iter().position(|&c| c == 0).unwrap_or(rest.len())]
    }

    /// Iterates over the alleles, REF first. Needs [`Unpack::STR`].
    #[inline]
    pub fn alleles(&self) -> impl ExactSizeIterator<Item = &[u8]> + '_ {
        (0..self.n_allele as usize).map(move |i| self.allele(i))
    }

    /// Returns the decoded alleles, REF first, each followed by a NUL, or `None` until
    /// [`Unpack::STR`] decodes them. Records whose alleles start with the same ones share a
    /// prefix.
    #[inline]
    #[must_use]
    pub fn decoded_alleles(&self) -> Option<&[u8]> {
        (self.unpacked & Unpack::STR.bits() != 0 || self.shared.is_empty()).then_some(&self.d.als)
    }

    /// Iterates over the alleles, REF first, without decoding them: from the encoded block
    /// until [`Unpack::STR`] decodes them, as [`Self::alleles`] gives them afterwards.
    pub fn peek_alleles(&self) -> impl ExactSizeIterator<Item = &[u8]> + '_ {
        let decoded = self.unpacked & Unpack::STR.bits() != 0 || self.shared.is_empty();
        PeekAlleles {
            rec: self,
            decoded,
            i: 0,
            p: if decoded {
                0
            } else {
                skip_sized(&self.shared, 0)
            },
        }
    }

    /// Returns the FILTER IDs (empty when `.`). Needs [`Unpack::FLT`].
    ///
    /// # Panics
    ///
    /// Panics if FILTER is not unpacked.
    #[inline]
    #[must_use]
    pub fn filters(&self) -> &[i32] {
        self.assert_unpacked(Unpack::FLT);
        &self.d.flt
    }

    /// Returns the decoded INFO fields, including ones marked for removal. Needs
    /// [`Unpack::INFO`].
    ///
    /// # Panics
    ///
    /// Panics if INFO is not unpacked.
    #[inline]
    #[must_use]
    pub fn info_fields(&self) -> &[InfoField] {
        self.assert_unpacked(Unpack::INFO);
        &self.d.info[..self.n_info as usize]
    }

    /// Returns the decoded FORMAT fields, including ones marked for removal. Needs
    /// [`Unpack::FMT`].
    ///
    /// # Panics
    ///
    /// Panics if FORMAT is not unpacked.
    #[inline]
    #[must_use]
    pub fn format_fields(&self) -> &[FormatField] {
        if self.n_sample > 0 {
            self.assert_unpacked(Unpack::FMT);
        }
        &self.d.fmt[..(self.n_fmt as usize).min(self.d.fmt.len())]
    }

    /// Returns the value bytes of an INFO field (HTSlib's `vptr`, `vptr_len` bytes).
    #[inline]
    #[must_use]
    pub fn info_data<'a>(&'a self, field: &'a InfoField) -> &'a [u8] {
        let len = field.vptr_len as usize;
        match &field.data {
            FieldData::Removed => &[],
            FieldData::InBlock(start) => {
                let s = start + field.vptr_off as usize;
                &self.shared[s..s + len]
            }
            FieldData::Owned(v) => &v[field.vptr_off as usize..field.vptr_off as usize + len],
        }
    }

    /// Returns the value bytes of a FORMAT field for all samples (HTSlib's `p`, `p_len` bytes).
    #[inline]
    #[must_use]
    pub fn format_data<'a>(&'a self, field: &'a FormatField) -> &'a [u8] {
        let len = field.p_len as usize;
        match &field.data {
            FieldData::Removed => &[],
            FieldData::InBlock(start) => {
                let s = start + field.p_off as usize;
                &self.indiv[s..s + len]
            }
            FieldData::Owned(v) => &v[field.p_off as usize..field.p_off as usize + len],
        }
    }

    /// Returns the value bytes of one sample of a FORMAT field.
    #[inline]
    #[must_use]
    pub fn format_sample_data<'a>(&'a self, field: &'a FormatField, sample: usize) -> &'a [u8] {
        let size = field.size();
        &self.format_data(field)[sample * size..(sample + 1) * size]
    }

    /// The whole encoded block of an INFO field: key, descriptor, and values.
    pub(crate) fn info_block<'a>(&'a self, field: &'a InfoField) -> &'a [u8] {
        let len = (field.vptr_off + field.vptr_len) as usize;
        match &field.data {
            FieldData::Removed => &[],
            FieldData::InBlock(start) => &self.shared[*start..start + len],
            FieldData::Owned(v) => &v[..len],
        }
    }

    pub(crate) fn format_block<'a>(&'a self, field: &'a FormatField) -> &'a [u8] {
        let len = (field.p_off + field.p_len) as usize;
        match &field.data {
            FieldData::Removed => &[],
            FieldData::InBlock(start) => &self.indiv[*start..start + len],
            FieldData::Owned(v) => &v[..len],
        }
    }

    #[inline]
    #[track_caller]
    fn assert_unpacked(&self, which: Unpack) {
        assert!(
            self.shared.is_empty() || self.unpacked & which.bits() == which.bits(),
            "record fields used before unpack({which:?})"
        );
    }

    // ---- unpacking ----

    /// Decodes parts of the record (HTSlib's `bcf_unpack`). Unpacking FILTER implies the
    /// alleles, and INFO all site fields. Already decoded parts are not redone.
    pub fn unpack(&mut self, which: Unpack) {
        if self.shared.is_empty() {
            return;
        }
        let mut which = which.bits();
        if which & Unpack::FLT.bits() != 0 {
            which |= Unpack::STR.bits();
        }
        if which & Unpack::INFO.bits() != 0 {
            which |= Unpack::SHR.bits();
        }
        let mut ptr = 0usize;
        if which & Unpack::STR.bits() != 0 && self.unpacked & Unpack::STR.bits() == 0 {
            self.d.id.clear();
            ptr = fmt_sized_array_into(&mut self.d.id, &self.shared, ptr);
            self.unpack_size[0] = ptr;
            self.d.als.clear();
            self.d.allele_starts.clear();
            let start = ptr;
            for _ in 0..self.n_allele {
                self.d.allele_starts.push(self.d.als.len() as u32);
                ptr = fmt_sized_array_into(&mut self.d.als, &self.shared, ptr);
                self.d.als.push(0);
            }
            self.unpack_size[1] = ptr - start;
            self.unpacked |= Unpack::STR.bits();
        }
        if which & Unpack::FLT.bits() != 0 && self.unpacked & Unpack::FLT.bits() == 0 {
            let mut p = self.unpack_size[0] + self.unpack_size[1];
            let start = p;
            self.d.flt.clear();
            if self.shared[p] >> 4 != 0 {
                let (n, ty, used) = encoding::dec_size(&self.shared[p..]);
                p += used;
                for _ in 0..n.max(0) {
                    let (v, used) = encoding::dec_int1(&self.shared[p..], ty);
                    self.d.flt.push(v as i32);
                    p += used;
                }
            } else {
                p += 1;
            }
            self.unpack_size[2] = p - start;
            self.unpacked |= Unpack::FLT.bits();
        }
        if which & Unpack::INFO.bits() != 0 && self.unpacked & Unpack::INFO.bits() == 0 {
            let mut p = self.unpack_size[0] + self.unpack_size[1] + self.unpack_size[2];
            // Decoded in place, into fields built where they live (edited fields' buffers are
            // kept for reuse first).
            self.d.clear_info();
            self.d
                .info
                .resize_with(self.n_info as usize, InfoField::default);
            for field in &mut self.d.info {
                p = unpack_info_core1(&self.shared, p, field);
            }
            self.unpacked |= Unpack::INFO.bits();
        }
        if which & Unpack::FMT.bits() != 0
            && self.n_sample > 0
            && self.unpacked & Unpack::FMT.bits() == 0
        {
            let mut p = 0;
            self.d.clear_fmt();
            self.d
                .fmt
                .resize_with(self.n_fmt as usize, FormatField::default);
            for field in &mut self.d.fmt {
                p = unpack_fmt_core1(&self.indiv, p, self.n_sample as usize, field);
            }
            self.unpacked |= Unpack::FMT.bits();
        }
    }

    /// The offset of FILTER in the encoded shared block, past ID and the alleles.
    fn encoded_filter_offset(&self) -> usize {
        let mut p = 0;
        for _ in 0..=self.n_allele {
            p = skip_sized(&self.shared, p);
        }
        p
    }

    /// The offset of INFO in the encoded shared block, past FILTER.
    pub(crate) fn encoded_info_offset(&self) -> usize {
        skip_sized(&self.shared, self.encoded_filter_offset())
    }

    /// Whether one of the FILTER IDs, or -1 for `.`, is in `ids` (the FILTER test of HTSlib's
    /// `_reader_fill_buffer`). Undecoded FILTER is read from the encoded block, leaving ID,
    /// REF and ALT encoded.
    pub(crate) fn has_filter_in(&mut self, ids: &[i32]) -> bool {
        if self.shared.is_empty() || self.unpacked & Unpack::FLT.bits() != 0 {
            self.unpack(Unpack::FLT);
            let flt = self.filters();
            return if flt.is_empty() {
                ids.contains(&-1)
            } else {
                flt.iter().any(|f| ids.contains(f))
            };
        }
        let shared = &self.shared;
        let mut p = self.encoded_filter_offset();
        let mut n = 0;
        if shared[p] >> 4 != 0 {
            let (len, ty, used) = encoding::dec_size(&shared[p..]);
            p += used;
            for _ in 0..len.max(0) {
                let (v, used) = encoding::dec_int1(&shared[p..], ty);
                p += used;
                if ids.contains(&(v as i32)) {
                    return true;
                }
                n += 1;
            }
        }
        n == 0 && ids.contains(&-1)
    }

    // ---- syncing ----

    /// Re-encodes edited parts into the BCF blocks (HTSlib's `bcf1_sync`), dropping fields
    /// marked for removal. Called before writing BCF and by [`Self::copy_from`].
    pub fn sync(&mut self) {
        let shared_rebuilt = if self.shared.is_empty() {
            // A record built through the API: encode everything.
            let mut tmp = std::mem::take(&mut self.shared);
            tmp.clear();
            self.sync_id(&mut tmp);
            self.unpack_size[0] = tmp.len();
            let prev = tmp.len();
            self.sync_alleles(&mut tmp);
            self.unpack_size[1] = tmp.len() - prev;
            let prev = tmp.len();
            self.sync_filter(&mut tmp);
            self.unpack_size[2] = tmp.len() - prev;
            self.sync_info(&mut tmp);
            self.shared = tmp;
            true
        } else if self.d.shared_dirty != 0 {
            if self.unpacked & Unpack::STR.bits() == 0 {
                self.unpack(Unpack::STR);
            }
            let mut tmp = Vec::with_capacity(self.shared.len() + 16);
            let mut ori = 0usize;
            let dirty = self.d.shared_dirty;
            // ID
            if dirty & DIRTY_ID != 0 {
                self.sync_id(&mut tmp);
            } else {
                tmp.extend_from_slice(&self.shared[ori..ori + self.unpack_size[0]]);
            }
            ori += self.unpack_size[0];
            self.unpack_size[0] = tmp.len();
            let mut prev = tmp.len();
            // REF+ALT
            if dirty & DIRTY_ALS != 0 {
                self.sync_alleles(&mut tmp);
            } else {
                tmp.extend_from_slice(&self.shared[ori..ori + self.unpack_size[1]]);
                if self.rlen == 0 && self.n_allele > 0 {
                    self.rlen = self.allele(0).len() as i64;
                }
            }
            ori += self.unpack_size[1];
            self.unpack_size[1] = tmp.len() - prev;
            prev = tmp.len();
            if self.unpacked & Unpack::FLT.bits() != 0 {
                if dirty & DIRTY_FLT != 0 {
                    self.sync_filter(&mut tmp);
                } else if !self.d.flt.is_empty() {
                    tmp.extend_from_slice(&self.shared[ori..ori + self.unpack_size[2]]);
                } else {
                    enc_vint(&mut tmp, &[], None);
                }
                ori += self.unpack_size[2];
                self.unpack_size[2] = tmp.len() - prev;
                if self.unpacked & Unpack::INFO.bits() != 0 && dirty & DIRTY_INF != 0 {
                    self.sync_info(&mut tmp);
                    ori = self.shared.len();
                }
            }
            if ori < self.shared.len() {
                tmp.extend_from_slice(&self.shared[ori..]);
            }
            self.shared = tmp;
            true
        } else {
            false
        };
        if shared_rebuilt && self.unpacked & Unpack::INFO.bits() != 0 {
            let mut off = self.unpack_size[0] + self.unpack_size[1] + self.unpack_size[2];
            let n = self.n_info as usize;
            let Decoded { info, spare, .. } = &mut self.d;
            for info in &mut info[..n] {
                // The field's own buffer is kept for later edits.
                if let FieldData::Owned(buf) =
                    std::mem::replace(&mut info.data, FieldData::InBlock(off))
                {
                    keep_spare(spare, buf);
                }
                off += (info.vptr_len + info.vptr_off) as usize;
            }
        }

        if self.n_sample > 0 && self.n_fmt > 0 && (self.indiv.is_empty() || self.d.indiv_dirty) {
            // A record built through the API has no block to read from: reuse its buffer.
            let mut tmp = if self.indiv.is_empty() {
                std::mem::take(&mut self.indiv)
            } else {
                Vec::with_capacity(self.indiv.len())
            };
            let mut irm: Option<usize> = None;
            let n_fmt = (self.n_fmt as usize).min(self.d.fmt.len());
            for i in 0..n_fmt {
                if self.d.fmt[i].is_removed() {
                    irm.get_or_insert(i);
                    continue;
                }
                tmp.extend_from_slice(self.format_block(&self.d.fmt[i]));
                if let Some(r) = irm.as_mut() {
                    self.d.fmt.swap(*r, i);
                    while *r <= i && !self.d.fmt[*r].is_removed() {
                        *r += 1;
                    }
                }
            }
            if let Some(r) = irm {
                self.n_fmt = r as u32;
            }
            self.indiv = tmp;
            let mut off = 0usize;
            let n = self.n_fmt as usize;
            let Decoded { fmt, spare, .. } = &mut self.d;
            for f in &mut fmt[..n] {
                if let FieldData::Owned(buf) =
                    std::mem::replace(&mut f.data, FieldData::InBlock(off))
                {
                    keep_spare(spare, buf);
                }
                off += (f.p_len + f.p_off) as usize;
            }
        }
        if self.n_sample == 0 {
            self.n_fmt = 0;
        }
        self.d.shared_dirty = 0;
        self.d.indiv_dirty = false;
    }

    fn sync_id(&self, out: &mut Vec<u8>) {
        if !self.d.id.is_empty() && self.d.id != b"." {
            enc_vchar(out, &self.d.id);
        } else {
            enc_size(out, 0, BcfType::Char as u8);
        }
    }

    fn sync_alleles(&mut self, out: &mut Vec<u8>) {
        for i in 0..self.n_allele as usize {
            let start = self.d.allele_starts[i] as usize;
            let rest = &self.d.als[start..];
            let a = &rest[..rest.iter().position(|&c| c == 0).unwrap_or(rest.len())];
            enc_vchar(out, a);
        }
        if self.rlen == 0 && self.n_allele > 0 {
            self.rlen = self.allele(0).len() as i64;
        }
    }

    fn sync_filter(&self, out: &mut Vec<u8>) {
        enc_vint(out, &self.d.flt, None);
    }

    fn sync_info(&mut self, out: &mut Vec<u8>) {
        let mut irm: Option<usize> = None;
        let n = self.n_info as usize;
        for i in 0..n {
            if self.d.info[i].is_removed() {
                irm.get_or_insert(i);
                continue;
            }
            out.extend_from_slice(self.info_block(&self.d.info[i]));
            if let Some(r) = irm.as_mut() {
                self.d.info.swap(*r, i);
                while *r <= i && !self.d.info[*r].is_removed() {
                    *r += 1;
                }
            }
        }
        if let Some(r) = irm {
            self.n_info = r as u32;
        }
    }

    /// Copies `src` into `self`, syncing `src` first (HTSlib's `bcf_copy`). Decoded fields are
    /// not copied: unpack the copy as if it had been read. Unlike `bcf_copy`, the 64-bit flag
    /// is kept, so the copy is still refused as BCF.
    pub fn copy_from(&mut self, src: &mut Self) {
        src.sync();
        self.clear();
        self.rid = src.rid;
        self.pos = src.pos;
        self.rlen = src.rlen;
        self.unpacked = src.unpacked & IS_64BIT;
        self.qual = src.qual;
        self.n_info = src.n_info;
        self.n_allele = src.n_allele;
        self.n_fmt = src.n_fmt;
        self.n_sample = src.n_sample;
        self.shared.extend_from_slice(&src.shared);
        self.indiv.extend_from_slice(&src.indiv);
    }

    /// Keeps only the samples in `imap` (indices into the current samples, -1 to skip), in
    /// that order (HTSlib's `bcf_subset`). An empty `imap` drops all sample data.
    pub fn subset(&mut self, imap: &[i32]) {
        let mut ind = Vec::new();
        if imap.is_empty() {
            self.n_sample = 0;
        } else {
            let mut p = 0;
            let n_sample = self.n_sample as usize;
            for _ in 0..self.n_fmt {
                let mut f = FormatField::default();
                p = unpack_fmt_core1(&self.indiv, p, n_sample, &mut f);
                let FieldData::InBlock(start) = f.data else {
                    unreachable!("unpacking reads from the block");
                };
                let data_start = start + f.p_off as usize;
                enc_int1(&mut ind, f.id);
                enc_size(&mut ind, f.n.max(0) as usize, f.ty);
                let size = f.size();
                for &j in imap {
                    if j >= 0 {
                        let s = data_start + j as usize * size;
                        ind.extend_from_slice(&self.indiv[s..s + size]);
                    }
                }
            }
            self.n_sample = imap.iter().filter(|&&j| j >= 0).count() as u32;
        }
        if self.n_sample == 0 {
            self.n_fmt = 0;
        }
        self.indiv = ind;
        self.unpacked &= !Unpack::FMT.bits();
    }

    /// Keeps the samples selected with [`Header::set_samples`] (HTSlib's `bcf_subset_format`).
    pub fn subset_format(&mut self, header: &Header) {
        let Some(keep) = &header.keep_samples else {
            return;
        };
        if header.sample_count() == 0 {
            self.indiv.clear();
            self.n_sample = 0;
            return;
        }
        let mut out = Vec::with_capacity(self.indiv.len());
        let mut p = 0;
        self.d.fmt.clear();
        for _ in 0..self.n_fmt {
            let mut f = FormatField::default();
            p = unpack_fmt_core1(&self.indiv, p, self.n_sample as usize, &mut f);
            let FieldData::InBlock(start) = f.data else {
                unreachable!("unpacking reads from the block");
            };
            let new_start = out.len();
            out.extend_from_slice(&self.indiv[start..start + f.p_off as usize]);
            let size = f.size();
            let data = start + f.p_off as usize;
            for j in 0..keep.n_original {
                if keep.contains(j) {
                    out.extend_from_slice(&self.indiv[data + j * size..data + (j + 1) * size]);
                }
            }
            f.p_len = (out.len() - new_start) as u32 - f.p_off;
            f.data = FieldData::InBlock(new_start);
            self.d.fmt.push(f);
        }
        self.indiv = out;
        self.unpacked |= Unpack::FMT.bits();
        self.n_sample = header.sample_count() as u32;
    }
}

/// Appends a sized typed array as text, as HTSlib's `bcf_fmt_sized_array` does for ID and
/// alleles, returning the position after it.
pub(crate) fn fmt_sized_array_into(out: &mut Vec<u8>, data: &[u8], p: usize) -> usize {
    let (n, ty, used) = encoding::dec_size(&data[p..]);
    let start = p + used;
    let n = n.max(0) as usize;
    crate::variant::vcf::format::fmt_array(out, n, ty, &data[start..]);
    start + (n << type_shift(ty))
}

/// The iterator of [`Record::peek_alleles`].
struct PeekAlleles<'r> {
    rec: &'r Record,
    decoded: bool,
    /// The next allele.
    i: usize,
    /// The offset of the next allele in the encoded block.
    p: usize,
}

impl<'r> Iterator for PeekAlleles<'r> {
    type Item = &'r [u8];

    #[inline]
    fn next(&mut self) -> Option<&'r [u8]> {
        if self.i >= self.rec.n_allele as usize {
            return None;
        }
        self.i += 1;
        if self.decoded {
            return Some(self.rec.allele(self.i - 1));
        }
        // As `bcf_fmt_sized_array` decodes it: up to the first NUL.
        let data = &self.rec.shared;
        let (n, ty, used) = encoding::dec_size(&data[self.p..]);
        let start = self.p + used;
        self.p = start + ((n.max(0) as usize) << type_shift(ty));
        let a = &data[start..self.p];
        Some(&a[..a.iter().position(|&c| c == 0).unwrap_or(a.len())])
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = (self.rec.n_allele as usize).saturating_sub(self.i);
        (n, Some(n))
    }
}

impl ExactSizeIterator for PeekAlleles<'_> {}

/// Returns the offset after the sized typed array at `p` (a type descriptor and its values).
#[inline]
fn skip_sized(data: &[u8], p: usize) -> usize {
    let (n, ty, used) = encoding::dec_size(&data[p..]);
    p + used + ((n.max(0) as usize) << type_shift(ty))
}

/// Decodes one INFO field at `p` (HTSlib's `bcf_unpack_info_core1`).
pub(crate) fn unpack_info_core1(data: &[u8], p: usize, info: &mut InfoField) -> usize {
    let start = p;
    let (key, used) = encoding::dec_typed_int1(&data[p..]);
    let mut q = p + used;
    let (len, ty, used) = encoding::dec_size(&data[q..]);
    q += used;
    info.key = key as i32;
    info.len = len as i32;
    info.ty = ty;
    info.vptr_off = (q - start) as u32;
    info.data = FieldData::InBlock(start);
    info.v1 = 0;
    let mut bytes = len.max(0) as usize;
    if len == 1 {
        match ty {
            1 | 7 => info.v1 = i64::from(data[q] as i8),
            2 => {
                info.v1 = i64::from(encoding::le_i16(&data[q..]));
                bytes <<= 1;
            }
            3 => {
                info.v1 = i64::from(encoding::le_i32(&data[q..]));
                bytes <<= 2;
            }
            5 => {
                info.v1 = i64::from(encoding::le_u32(&data[q..]));
                bytes <<= 2;
            }
            4 => {
                info.v1 = encoding::le_i64(&data[q..]);
                bytes <<= 3;
            }
            _ => {}
        }
    } else {
        bytes <<= type_shift(ty);
    }
    info.vptr_len = bytes as u32;
    q + bytes
}

/// Decodes one FORMAT field at `p` (HTSlib's `bcf_unpack_fmt_core1`).
pub(crate) fn unpack_fmt_core1(
    data: &[u8],
    p: usize,
    n_sample: usize,
    fmt: &mut FormatField,
) -> usize {
    let start = p;
    let (id, used) = encoding::dec_typed_int1(&data[p..]);
    let mut q = p + used;
    let (n, ty, used) = encoding::dec_size(&data[q..]);
    q += used;
    fmt.id = id as i32;
    fmt.n = n as i32;
    fmt.ty = ty;
    fmt.size = (n << type_shift(ty)) as i32;
    fmt.p_off = (q - start) as u32;
    fmt.data = FieldData::InBlock(start);
    let len = n_sample * fmt.size.max(0) as usize;
    fmt.p_len = len as u32;
    q + len
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peek_alleles_match_decoded_alleles() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .expect("the header parses");
        for alts in [".", "T", "<*>", "AT,ACGTACGTACGTAC,<DEL>"] {
            let mut line = format!("1\t5\trs7\tACGTACGTACGTA\t{alts}\t.\t.\t.").into_bytes();
            let mut parsed = Record::new();
            crate::variant::vcf::parse(&mut header, &mut line, &mut parsed).expect("parses");
            let mut packed = Vec::new();
            crate::variant::bcf::pack_record(&mut parsed, &mut packed);
            let mut rec = Record::new();
            crate::variant::bcf::unpack_record(&packed, &mut rec).expect("unpacks");
            let peeked: Vec<Vec<u8>> = rec.peek_alleles().map(<[u8]>::to_vec).collect();
            assert_eq!(rec.peek_alleles().len(), rec.allele_count());
            rec.unpack(Unpack::STR);
            let decoded: Vec<Vec<u8>> = rec.alleles().map(<[u8]>::to_vec).collect();
            assert_eq!(peeked, decoded, "{alts}");
            let again: Vec<Vec<u8>> = rec.peek_alleles().map(<[u8]>::to_vec).collect();
            assert_eq!(again, decoded, "{alts}");
        }
    }

    #[test]
    fn allocated_bytes_count_blocks_and_edited_fields() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
              ##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        )
        .expect("the header parses");
        let mut rec = Record::new();
        assert_eq!(rec.allocated_bytes(), 0);
        let mut line = b"1\t5\t.\tA\tT\t.\t.\t.\tDP\t1\t2".to_vec();
        crate::variant::vcf::parse(&mut header, &mut line, &mut rec).expect("parses");
        let parsed = rec.allocated_bytes();
        assert!(
            parsed >= rec.shared_bytes().len() + rec.indiv_bytes().len(),
            "{parsed} bytes"
        );
        rec.unpack(Unpack::ALL);
        rec.update_format_i32(&header, b"DP", &[300_000, 400_000])
            .expect("DP is a FORMAT tag");
        let edited = rec.allocated_bytes();
        assert!(edited > parsed, "{edited} after editing, {parsed} before");
        // Cleared records keep their buffers for reuse.
        rec.clear();
        let cleared = rec.allocated_bytes();
        assert!(
            cleared >= parsed,
            "{cleared} after clearing, {parsed} before"
        );
    }

    #[test]
    fn remove_all_info_cut_matches_decoding() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n##FILTER=<ID=q10,Description=\"q\">\n\
              ##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
              ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
              ##INFO=<ID=S,Number=.,Type=String,Description=\"s\">\n\
              ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n##contig=<ID=1>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n",
        )
        .expect("the header parses");
        let lines = [
            "1\t5\trs1\tACGT\tA,ACGTTT\t30\tq10\tDP=300000;S=abc,defghijk\tGT\t0/1",
            "1\t9\t.\tA\t<DEL>\t.\t.\tEND=500;DP=3\tGT\t1|1",
            "1\t12\t.\tA\tT\t.\tPASS\t.\tGT\t./.",
        ];
        for line in lines {
            let mut text = line.as_bytes().to_vec();
            let mut parsed = Record::new();
            crate::variant::vcf::parse(&mut header, &mut text, &mut parsed).expect("parses");
            let mut packed = Vec::new();
            crate::variant::bcf::pack_record(&mut parsed, &mut packed);
            let mut out = Vec::new();
            for decode_first in [false, true] {
                let mut rec = Record::new();
                crate::variant::bcf::unpack_record(&packed, &mut rec).expect("unpacks");
                if decode_first {
                    rec.unpack(Unpack::INFO);
                }
                rec.remove_all_info(&header);
                let mut vcf = Vec::new();
                crate::variant::vcf::format_record(&header, &rec, &mut vcf).expect("formats");
                let mut bcf = Vec::new();
                crate::variant::bcf::encode_record(&header, &mut rec, &mut bcf).expect("encodes");
                out.push((vcf, bcf, rec.rlen));
            }
            assert_eq!(out[0], out[1], "{line}");
        }
    }

    #[test]
    fn has_filter_in_reads_the_encoded_block_like_decoding() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n##FILTER=<ID=q10,Description=\"q\">\n\
              ##FILTER=<ID=s50,Description=\"s\">\n##contig=<ID=1>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .expect("the header parses");
        let (q10, s50) = (
            header.id(b"q10").expect("q10") as i32,
            header.id(b"s50").expect("s50") as i32,
        );
        let sets: [&[i32]; 6] = [&[-1], &[0], &[q10], &[s50], &[q10, -1], &[0, s50]];
        for filter in [".", "PASS", "q10", "s50", "q10;s50"] {
            let mut line =
                format!("1\t5\trs1;rs2\tACGTACGTACG\tA,ACGTTTTTTTTTT\t.\t{filter}\t.").into_bytes();
            let mut parsed = Record::new();
            crate::variant::vcf::parse(&mut header, &mut line, &mut parsed).expect("parses");
            let mut packed = Vec::new();
            crate::variant::bcf::pack_record(&mut parsed, &mut packed);
            for ids in sets {
                let mut rec = Record::new();
                crate::variant::bcf::unpack_record(&packed, &mut rec).expect("unpacks");
                let mut decoded = rec.clone();
                decoded.unpack(Unpack::FLT);
                let want = if decoded.filters().is_empty() {
                    ids.contains(&-1)
                } else {
                    decoded.filters().iter().any(|f| ids.contains(f))
                };
                assert_eq!(rec.has_filter_in(ids), want, "FILTER {filter}, ids {ids:?}");
                assert_eq!(rec.unpacked & Unpack::STR.bits(), 0);
            }
        }
    }

    #[test]
    fn strerror_matches_htslib() {
        let e = RecordError::CTG_UNDEF | RecordError::TAG_UNDEF;
        assert_eq!(
            e.describe(1024),
            "Contig not defined in header,Tag not defined in header"
        );
        assert_eq!(RecordError::empty().describe(1024), "");
    }

    #[test]
    fn add_id_finds_whole_lists() {
        let mut rec = Record::new();
        rec.update_id(Some(b"rs1;rs2"));
        rec.add_id(b"rs1;rs2");
        assert_eq!(rec.id(), b"rs1;rs2");
        rec.add_id(b"rs2");
        assert_eq!(rec.id(), b"rs1;rs2");
        rec.add_id(b"rs");
        assert_eq!(rec.id(), b"rs1;rs2;rs");
        rec.update_id(None);
        rec.add_id(b"x");
        assert_eq!(rec.id(), b"x");
    }

    #[test]
    fn format_lookup_in_the_encoded_block_matches_decoding() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n\
              ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
              ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
              ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
              ##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"v\">\n\
              ##FORMAT=<ID=S,Number=1,Type=String,Description=\"s\">\n\
              ##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"x\">\n##contig=<ID=1>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        )
        .expect("the header parses");
        let lines = [
            "1\t5\t.\tA\tT,G\t.\t.\t.\tGT:DP:AD:VAF:S\t0/1:300:1,2,3:0.5,.:ab\t1|1:.:4,.:.:cdef",
            "1\t9\t.\tA\tT\t.\t.\t.\tGT:S:DP\t0/0:x:7\t./.:yz:40000",
        ];
        for line in lines {
            let mut text = line.as_bytes().to_vec();
            let mut parsed = Record::new();
            crate::variant::vcf::parse(&mut header, &mut text, &mut parsed).expect("parses");
            let mut packed = Vec::new();
            crate::variant::bcf::pack_record(&mut parsed, &mut packed);
            let mut rec = Record::new();
            crate::variant::bcf::unpack_record(&packed, &mut rec).expect("unpacks");
            let mut decoded = rec.clone();
            decoded.unpack(Unpack::FMT);
            for tag in [&b"GT"[..], b"DP", b"AD", b"X", b"VAF", b"S", b"NONE"] {
                let (mut a, mut b) = (Vec::new(), Vec::new());
                let got = rec.get_format_i32(&header, tag, &mut a);
                assert_eq!(
                    got,
                    decoded.get_format_i32(&header, tag, &mut b),
                    "{line} {tag:?}"
                );
                assert_eq!(a, b, "{line} {tag:?}");
                let (mut a, mut b) = (Vec::new(), Vec::new());
                let got = rec.get_format_f32(&header, tag, &mut a);
                assert_eq!(
                    got,
                    decoded.get_format_f32(&header, tag, &mut b),
                    "{line} {tag:?}"
                );
                let bits = |v: &[f32]| v.iter().map(|f| f.to_bits()).collect::<Vec<_>>();
                assert_eq!(bits(&a), bits(&b), "{line} {tag:?}");
                assert_eq!(rec.unpacked & Unpack::FMT.bits(), 0, "{line} {tag:?}");
            }
            assert_eq!(
                rec.clone().get_format_chars(&header, b"S"),
                decoded.get_format_chars(&header, b"S")
            );
        }
    }

    #[test]
    fn format_reads_of_a_subset_record_with_the_full_header_fail() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n\
              ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
              ##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"v\">\n##contig=<ID=1>\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        )
        .expect("the header parses");
        let mut text = b"1\t5\t.\tA\tT,G\t.\t.\t.\tAD:VAF\t1,2,3:0.5,.\t4,5,6:.,.".to_vec();
        let mut parsed = Record::new();
        crate::variant::vcf::parse(&mut header, &mut text, &mut parsed).expect("parses");
        let mut packed = Vec::new();
        crate::variant::bcf::pack_record(&mut parsed, &mut packed);
        let mut rec = Record::new();
        crate::variant::bcf::unpack_record(&packed, &mut rec).expect("unpacks");
        rec.subset(&[1]);
        let mut decoded = rec.clone();
        decoded.unpack(Unpack::FMT);
        for r in [&mut rec, &mut decoded] {
            let (mut ints, mut floats) = (Vec::new(), Vec::new());
            assert_eq!(
                r.get_format_i32(&header, b"AD", &mut ints),
                Err(GetError::TooFewSamples)
            );
            assert_eq!(
                r.get_format_f32(&header, b"VAF", &mut floats),
                Err(GetError::TooFewSamples)
            );
        }
    }
}
