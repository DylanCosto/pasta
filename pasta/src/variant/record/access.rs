//! Reading INFO and FORMAT values (HTSlib's `bcf_get_info_values`, `bcf_get_format_values`,
//! `bcf_get_format_string`) and typed views of encoded values.

use super::{Record, Unpack};
use crate::variant::{
    Header, HeaderRecordType, ValueType,
    encoding::{
        self, FLOAT_VECTOR_END, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING, INT16_VECTOR_END,
        INT32_MISSING, INT32_VECTOR_END, INT64_MISSING,
    },
    header::SpecialTag,
    log,
};

/// Why a value could not be read (HTSlib's negative return codes).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum GetError {
    /// The tag is not defined in the header for this kind of field (-1).
    NotInHeader,
    /// The tag is defined with a different type (-2).
    WrongType,
    /// The tag is not present in this record, or is marked for removal (-3).
    Absent,
    /// The record holds fewer samples than the header names, as when it was subset without
    /// the header. HTSlib reads past the record's values.
    TooFewSamples,
}

impl std::fmt::Display for GetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NotInHeader => "the tag is not defined in the header",
            Self::WrongType => "the tag has a different type",
            Self::Absent => "the tag is not present in the record",
            Self::TooFewSamples => "the record has fewer samples than the header",
        })
    }
}

impl std::error::Error for GetError {}

/// The values of a FORMAT field for all samples.
struct FormatView<'a> {
    /// The encoded type.
    ty: u8,
    /// Values per sample.
    n: usize,
    /// Bytes per sample.
    size: usize,
    data: &'a [u8],
}

/// A genotype allele as BCF encodes it: `(allele + 1) << 1 | phased`, 0 for a missing allele
/// (HTSlib's `bcf_gt_*` macros).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Genotype(pub i32);

impl Genotype {
    /// A missing allele (`.`).
    pub const MISSING: Self = Self(0);

    /// An unphased allele (`bcf_gt_unphased`).
    #[inline]
    #[must_use]
    pub fn unphased(allele: i32) -> Self {
        Self((allele + 1) << 1)
    }

    /// A phased allele (`bcf_gt_phased`).
    #[inline]
    #[must_use]
    pub fn phased(allele: i32) -> Self {
        Self(((allele + 1) << 1) | 1)
    }

    /// Returns whether the allele is missing (`bcf_gt_is_missing`).
    #[inline]
    #[must_use]
    pub fn is_missing(self) -> bool {
        self.0 >> 1 == 0
    }

    /// Returns whether the allele is phased (`bcf_gt_is_phased`).
    #[inline]
    #[must_use]
    pub fn is_phased(self) -> bool {
        self.0 & 1 != 0
    }

    /// Returns the allele index, -1 for missing (`bcf_gt_allele`).
    #[inline]
    #[must_use]
    pub fn allele(self) -> i32 {
        (self.0 >> 1) - 1
    }

    /// Returns whether this is the vector-end padding of a shorter ploidy.
    #[inline]
    #[must_use]
    pub fn is_vector_end(self) -> bool {
        self.0 == INT32_VECTOR_END
    }
}

/// A typed view of encoded BCF values (little-endian bytes of one BCF type).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TypedValues<'a> {
    /// No values.
    Null,
    /// 8-bit integers.
    Int8(&'a [u8]),
    /// 16-bit integers.
    Int16(&'a [u8]),
    /// 32-bit integers.
    Int32(&'a [u8]),
    /// 64-bit integers.
    Int64(&'a [u8]),
    /// 32-bit floats.
    Float(&'a [u8]),
    /// Characters.
    Char(&'a [u8]),
}

impl<'a> TypedValues<'a> {
    /// Views `data` as values of the raw BCF type code `ty`.
    #[inline]
    #[must_use]
    pub fn new(ty: u8, data: &'a [u8]) -> Self {
        match ty {
            1 => Self::Int8(data),
            2 => Self::Int16(data),
            3 => Self::Int32(data),
            4 => Self::Int64(data),
            5 => Self::Float(data),
            7 => Self::Char(data),
            _ => Self::Null,
        }
    }

    /// Returns the number of values.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Null => 0,
            Self::Int8(d) | Self::Char(d) => d.len(),
            Self::Int16(d) => d.len() / 2,
            Self::Int32(d) | Self::Float(d) => d.len() / 4,
            Self::Int64(d) => d.len() / 8,
        }
    }

    /// Returns whether there are no values.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns integer value `i` widened to `i32`, with the missing and vector-end sentinels
    /// mapped to [`INT32_MISSING`] and [`INT32_VECTOR_END`]; `None` for non-integer types.
    #[inline]
    #[must_use]
    pub fn int(&self, i: usize) -> Option<i32> {
        Some(match self {
            Self::Int8(d) => {
                let v = d[i] as i8;
                if v == INT8_MISSING {
                    INT32_MISSING
                } else if v == INT8_VECTOR_END {
                    INT32_VECTOR_END
                } else {
                    i32::from(v)
                }
            }
            Self::Int16(d) => {
                let v = encoding::le_i16(&d[i * 2..]);
                if v == INT16_MISSING {
                    INT32_MISSING
                } else if v == INT16_VECTOR_END {
                    INT32_VECTOR_END
                } else {
                    i32::from(v)
                }
            }
            Self::Int32(d) => encoding::le_i32(&d[i * 4..]),
            _ => return None,
        })
    }

    /// Returns float value `i` (missing and vector end keep their NaN payloads).
    #[inline]
    #[must_use]
    pub fn float(&self, i: usize) -> Option<f32> {
        match self {
            Self::Float(d) => Some(f32::from_bits(encoding::le_u32(&d[i * 4..]))),
            _ => None,
        }
    }

    /// Returns the raw bytes.
    #[inline]
    #[must_use]
    pub fn bytes(&self) -> &'a [u8] {
        match self {
            Self::Null => &[],
            Self::Int8(d)
            | Self::Int16(d)
            | Self::Int32(d)
            | Self::Int64(d)
            | Self::Float(d)
            | Self::Char(d) => d,
        }
    }
}

/// Appends the values of each sample, `n` per sample (`size` bytes of `width`-byte values), up
/// to and including the vector end `end`, padded with vector ends (the loop of
/// `bcf_get_format_values`). `dec` decodes a value to the 32-bit word HTSlib stores, and `word`
/// makes an output value of it.
#[inline(always)]
#[expect(
    clippy::too_many_arguments,
    reason = "the layout of an encoded FORMAT field and the two conversions"
)]
fn fill_words<T: Copy>(
    out: &mut Vec<T>,
    data: &[u8],
    nsmpl: usize,
    n: usize,
    size: usize,
    width: usize,
    end: u32,
    dec: impl Fn(&[u8]) -> u32,
    word: impl Fn(u32) -> T,
) {
    // Filled with vector ends first, so a sample's values stop at its first vector end.
    let start = out.len();
    out.resize(start + nsmpl * n, word(end));
    let dst = &mut out[start..];
    for s in 0..nsmpl {
        let sample = &data[s * size..s * size + n * width];
        let row = &mut dst[s * n..(s + 1) * n];
        for j in 0..n {
            let v = dec(&sample[j * width..(j + 1) * width]);
            row[j] = word(v);
            if v == end {
                break;
            }
        }
    }
}

/// Reads a FORMAT field into `out` as `bcf_get_format_values` reads it into an array of 32-bit
/// values: `n` per sample, padded with vector ends. Integers are widened to `i32` (with the
/// `i32` missing and vector-end values) and floats kept as they are; HTSlib stores the bits
/// whatever the header says, so integers read as floats, and floats read as integers, keep
/// their bits. `word` makes a value of each 32-bit word. Returns the number of values.
///
/// # Errors
///
/// Returns [`GetError::WrongType`] for another encoded type (where HTSlib exits), leaving
/// `out` as it is.
fn format_words<T: Copy>(
    field: &FormatView<'_>,
    nsmpl: usize,
    out: &mut Vec<T>,
    word: impl Fn(u32) -> T,
) -> Result<usize, GetError> {
    let (n, size, data) = (field.n, field.size, field.data);
    if data.len() < nsmpl * size {
        return Err(GetError::TooFewSamples);
    }
    if n == 0 || nsmpl == 0 {
        out.clear();
        return Ok(0);
    }
    let widen = |v: i32, missing: i32, end: i32| {
        if v == missing {
            INT32_MISSING
        } else if v == end {
            INT32_VECTOR_END
        } else {
            v
        }
    };
    let int_end = INT32_VECTOR_END as u32;
    // One loop per encoded width, so that the type is not matched for every value.
    match TypedValues::new(field.ty, &[]) {
        TypedValues::Int8(_) => {
            out.clear();
            let dec = |b: &[u8]| {
                let v = i32::from(b[0] as i8);
                widen(v, i32::from(INT8_MISSING), i32::from(INT8_VECTOR_END)) as u32
            };
            fill_words(out, data, nsmpl, n, size, 1, int_end, dec, word);
        }
        TypedValues::Int16(_) => {
            out.clear();
            let dec = |b: &[u8]| {
                let v = i32::from(encoding::le_i16(b));
                widen(v, i32::from(INT16_MISSING), i32::from(INT16_VECTOR_END)) as u32
            };
            fill_words(out, data, nsmpl, n, size, 2, int_end, dec, word);
        }
        TypedValues::Int32(_) => {
            out.clear();
            let dec = |b: &[u8]| encoding::le_i32(b) as u32;
            fill_words(out, data, nsmpl, n, size, 4, int_end, dec, word);
        }
        TypedValues::Float(_) => {
            out.clear();
            let end = FLOAT_VECTOR_END;
            fill_words(out, data, nsmpl, n, size, 4, end, encoding::le_u32, word);
        }
        TypedValues::Null | TypedValues::Int64(_) | TypedValues::Char(_) => {
            return Err(GetError::WrongType);
        }
    }
    Ok(nsmpl * n)
}

impl Record {
    fn info_index(
        &mut self,
        header: &Header,
        tag: &[u8],
        ty: ValueType,
    ) -> Result<Option<usize>, GetError> {
        let id = header.id(tag).ok_or(GetError::NotInHeader)?;
        self.info_index_id(header, id, ty)
    }

    /// [`Self::info_index`] for a dictionary ID (as [`Header::id`] returns).
    fn info_index_id(
        &mut self,
        header: &Header,
        id: u32,
        ty: ValueType,
    ) -> Result<Option<usize>, GetError> {
        if !header.id_exists(HeaderRecordType::Info, id) {
            return Err(GetError::NotInHeader);
        }
        if header.id_type(HeaderRecordType::Info, id) != Some(ty) {
            return Err(GetError::WrongType);
        }
        if self.unpacked & Unpack::INFO.bits() == 0 {
            self.unpack(Unpack::INFO);
        }
        Ok(self.d.info[..self.n_info as usize]
            .iter()
            .position(|f| f.key == id as i32))
    }

    /// Returns whether an INFO flag is set (HTSlib's `bcf_get_info_flag`).
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an INFO flag of the header.
    pub fn get_info_flag(&mut self, header: &Header, tag: &[u8]) -> Result<bool, GetError> {
        Ok(self.info_index(header, tag, ValueType::Flag)?.is_some())
    }

    /// Reads an integer INFO field into `out` (HTSlib's `bcf_get_info_int32`), with missing
    /// values as [`INT32_MISSING`]; stops at a vector end. Returns the number of values.
    /// Values stored as floats (a BCF whose header was changed) come out as their bits, as in
    /// HTSlib.
    ///
    /// # Errors
    ///
    /// Returns an error, leaving `out` as it is, when the tag is not an integer INFO tag, is
    /// absent, or is stored with no values or as characters (HTSlib's message is printed).
    pub fn get_info_i32(
        &mut self,
        header: &Header,
        tag: &[u8],
        out: &mut Vec<i32>,
    ) -> Result<usize, GetError> {
        let id = header.id(tag).ok_or(GetError::NotInHeader)?;
        self.get_info_i32_id(header, id, out)
    }

    /// [`Self::get_info_i32`] for a tag given by its dictionary ID ([`Header::id`]), which
    /// saves the name lookup when the same tag is read from many records.
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::get_info_i32`] does.
    pub fn get_info_i32_id(
        &mut self,
        header: &Header,
        id: u32,
        out: &mut Vec<i32>,
    ) -> Result<usize, GetError> {
        let i = self
            .info_index_id(header, id, ValueType::Integer)?
            .ok_or(GetError::Absent)?;
        self.info_words(header, i, out, |w| w as i32)
    }

    /// Reads an integer INFO field as 64-bit values (HTSlib's `bcf_get_info_int64`).
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an integer INFO tag or is absent.
    pub fn get_info_i64(
        &mut self,
        header: &Header,
        tag: &[u8],
        out: &mut Vec<i64>,
    ) -> Result<usize, GetError> {
        let mut tmp = Vec::new();
        self.get_info_i32(header, tag, &mut tmp)?;
        out.clear();
        out.extend(tmp.iter().map(|&v| {
            if v == INT32_MISSING {
                INT64_MISSING
            } else {
                i64::from(v)
            }
        }));
        Ok(out.len())
    }

    /// Reads a float INFO field into `out` (HTSlib's `bcf_get_info_float`); stops at a vector
    /// end. Values stored as integers (a BCF whose header was changed) come out as the bits
    /// of the `i32` values, as in HTSlib: `5` reads as a denormal, a missing value as `-0.0`.
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::get_info_i32`] does.
    pub fn get_info_f32(
        &mut self,
        header: &Header,
        tag: &[u8],
        out: &mut Vec<f32>,
    ) -> Result<usize, GetError> {
        let id = header.id(tag).ok_or(GetError::NotInHeader)?;
        self.get_info_f32_id(header, id, out)
    }

    /// [`Self::get_info_f32`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::get_info_f32`] does.
    pub fn get_info_f32_id(
        &mut self,
        header: &Header,
        id: u32,
        out: &mut Vec<f32>,
    ) -> Result<usize, GetError> {
        let i = self
            .info_index_id(header, id, ValueType::Float)?
            .ok_or(GetError::Absent)?;
        self.info_words(header, i, out, f32::from_bits)
    }

    /// Reads INFO field `i` into `out` as `bcf_get_info_values` reads a field into an array of
    /// 32-bit values, up to the first vector end: integers widened to `i32` (missing values
    /// as [`INT32_MISSING`]), floats as they are. HTSlib stores the bits whatever the header
    /// says, so integers read as floats, and floats read as integers, keep their bits:
    /// `word` makes a value of each 32-bit word. Returns the number of values.
    ///
    /// # Errors
    ///
    /// Returns [`GetError::Absent`] for a field marked for removal, and
    /// [`GetError::WrongType`] (with HTSlib's message) for one of another encoded type. `out`
    /// is left as it is on error.
    fn info_words<T>(
        &self,
        header: &Header,
        i: usize,
        out: &mut Vec<T>,
        word: impl Fn(u32) -> T,
    ) -> Result<usize, GetError> {
        let field = &self.d.info[i];
        if field.is_removed() {
            return Err(GetError::Absent);
        }
        let values = TypedValues::new(field.ty, self.info_data(field));
        let n = field.len().min(values.len());
        match values {
            TypedValues::Int8(_) | TypedValues::Int16(_) | TypedValues::Int32(_) => {
                out.clear();
                let ints = (0..n).map_while(|j| values.int(j));
                out.extend(
                    ints.take_while(|&v| v != INT32_VECTOR_END)
                        .map(|v| word(v as u32)),
                );
            }
            TypedValues::Float(data) => {
                out.clear();
                let bits = data.as_chunks::<4>().0[..n]
                    .iter()
                    .map(|c| u32::from_le_bytes(*c));
                out.extend(bits.take_while(|&b| b != FLOAT_VECTOR_END).map(word));
            }
            TypedValues::Null | TypedValues::Int64(_) | TypedValues::Char(_) => {
                log::error(
                    "bcf_get_info_values",
                    format_args!(
                        "Unexpected type {} at {}:{}",
                        field.ty,
                        String::from_utf8_lossy(self.chrom_safe(header)),
                        self.pos + 1
                    ),
                );
                return Err(GetError::WrongType);
            }
        }
        Ok(out.len())
    }

    /// Reads a string INFO field (HTSlib's `bcf_get_info_string`), returning its bytes.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not a string INFO tag or is absent.
    pub fn get_info_string(&mut self, header: &Header, tag: &[u8]) -> Result<&[u8], GetError> {
        let i = self
            .info_index(header, tag, ValueType::String)?
            .ok_or(GetError::Absent)?;
        let field = &self.d.info[i];
        if field.is_removed() {
            return Err(GetError::Absent);
        }
        let n = field.len();
        Ok(&self.info_data(field)[..n.min(field.vptr_len as usize)])
    }

    fn format_view(
        &self,
        header: &Header,
        tag: &[u8],
        ty: ValueType,
    ) -> Result<FormatView<'_>, GetError> {
        let id = header.id(tag).ok_or(GetError::NotInHeader)?;
        self.format_view_id(header, id, ty)
    }

    /// [`Self::format_view`] for a dictionary ID (as [`Header::id`] returns). Before FORMAT is
    /// decoded, finds the field in the encoded block without decoding the others.
    fn format_view_id(
        &self,
        header: &Header,
        id: u32,
        ty: ValueType,
    ) -> Result<FormatView<'_>, GetError> {
        if !header.id_exists(HeaderRecordType::Format, id) {
            return Err(GetError::NotInHeader);
        }
        let want = if header.special_id(SpecialTag::Gt) == Some(id) {
            ValueType::String
        } else {
            ty
        };
        if header.id_type(HeaderRecordType::Format, id) != Some(want) {
            return Err(GetError::WrongType);
        }
        if self.unpacked & Unpack::FMT.bits() == 0 && !self.shared.is_empty() {
            if self.n_sample == 0 {
                return Err(GetError::Absent);
            }
            // Read as `unpack_fmt_core1` decodes each field.
            let mut p = 0;
            for _ in 0..self.n_fmt {
                let (fid, used) = encoding::dec_typed_int1(&self.indiv[p..]);
                let mut q = p + used;
                let (n, ty, used) = encoding::dec_size(&self.indiv[q..]);
                q += used;
                let size = ((n << encoding::type_shift(ty)) as i32).max(0) as usize;
                let len = self.n_sample as usize * size;
                if fid as i32 == id as i32 {
                    return Ok(FormatView {
                        ty,
                        n: (n as i32).max(0) as usize,
                        size,
                        data: &self.indiv[q..q + len],
                    });
                }
                p = q + len;
            }
            return Err(GetError::Absent);
        }
        let n = (self.n_fmt as usize).min(self.d.fmt.len());
        let field = self.d.fmt[..n]
            .iter()
            .find(|f| f.id == id as i32)
            .ok_or(GetError::Absent)?;
        if field.is_removed() {
            return Err(GetError::Absent);
        }
        Ok(FormatView {
            ty: field.ty,
            n: field.n(),
            size: field.size(),
            data: self.format_data(field),
        })
    }

    /// Reads an integer FORMAT field for all samples (HTSlib's `bcf_get_format_int32`): the
    /// same number of values per sample, padded with [`INT32_VECTOR_END`]. Genotypes come out
    /// as [`Genotype`] encodings. Returns the number of values. Values stored as floats (a
    /// BCF whose header was changed) come out as their bits, as in HTSlib.
    ///
    /// # Errors
    ///
    /// Returns an error, leaving `out` as it is, when the tag is not an integer FORMAT tag
    /// (or GT), is absent, or is stored as characters, or when the record has fewer samples
    /// than `header`.
    pub fn get_format_i32(
        &mut self,
        header: &Header,
        tag: &[u8],
        out: &mut Vec<i32>,
    ) -> Result<usize, GetError> {
        let id = header.id(tag).ok_or(GetError::NotInHeader)?;
        self.get_format_i32_id(header, id, out)
    }

    /// [`Self::get_format_i32`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::get_format_i32`] does.
    pub fn get_format_i32_id(
        &mut self,
        header: &Header,
        id: u32,
        out: &mut Vec<i32>,
    ) -> Result<usize, GetError> {
        let field = self.format_view_id(header, id, ValueType::Integer)?;
        format_words(&field, header.sample_count(), out, |w| w as i32)
    }

    /// Reads the genotypes of all samples (HTSlib's `bcf_get_genotypes`).
    ///
    /// # Errors
    ///
    /// Returns an error when GT is not defined or absent.
    pub fn get_genotypes(
        &mut self,
        header: &Header,
        out: &mut Vec<i32>,
    ) -> Result<usize, GetError> {
        self.get_format_i32(header, b"GT", out)
    }

    /// Reads a float FORMAT field for all samples (HTSlib's `bcf_get_format_float`), padded
    /// with vector ends. Values stored as integers (a BCF whose header was changed) come out
    /// as the bits of the `i32` values, as in HTSlib.
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::get_format_i32`] does.
    pub fn get_format_f32(
        &mut self,
        header: &Header,
        tag: &[u8],
        out: &mut Vec<f32>,
    ) -> Result<usize, GetError> {
        let id = header.id(tag).ok_or(GetError::NotInHeader)?;
        self.get_format_f32_id(header, id, out)
    }

    /// [`Self::get_format_f32`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::get_format_f32`] does.
    pub fn get_format_f32_id(
        &mut self,
        header: &Header,
        id: u32,
        out: &mut Vec<f32>,
    ) -> Result<usize, GetError> {
        let field = self.format_view_id(header, id, ValueType::Float)?;
        format_words(&field, header.sample_count(), out, f32::from_bits)
    }

    /// Returns the raw bytes of a string FORMAT field for all samples (HTSlib's
    /// `bcf_get_format_char`): `n` bytes per sample, NUL-padded. Use
    /// [`Self::format_strings`] for one string per sample.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not a string FORMAT tag or is absent.
    pub fn get_format_chars(
        &mut self,
        header: &Header,
        tag: &[u8],
    ) -> Result<(&[u8], usize), GetError> {
        let field = self.format_view(header, tag, ValueType::String)?;
        let n = field.n * header.sample_count();
        Ok((&field.data[..n.min(field.data.len())], field.n))
    }

    /// Returns one string per sample of a string FORMAT field, each ending at its first NUL
    /// (HTSlib's `bcf_get_format_string`).
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not a string FORMAT tag or is absent.
    pub fn format_strings(&mut self, header: &Header, tag: &[u8]) -> Result<Vec<&[u8]>, GetError> {
        let (data, n) = self.get_format_chars(header, tag)?;
        Ok(data
            .chunks(n.max(1))
            .map(|c| &c[..c.iter().position(|&b| b == 0).unwrap_or(c.len())])
            .collect())
    }

    /// Returns a typed view of an INFO field's values.
    #[must_use]
    pub fn info_values<'a>(&'a self, field: &'a super::InfoField) -> TypedValues<'a> {
        TypedValues::new(field.ty, self.info_data(field))
    }

    /// Returns a typed view of one sample's values of a FORMAT field.
    #[must_use]
    pub fn format_values<'a>(
        &'a self,
        field: &'a super::FormatField,
        sample: usize,
    ) -> TypedValues<'a> {
        TypedValues::new(field.ty, self.format_sample_data(field, sample))
    }

    /// Finds a decoded INFO field by tag name (HTSlib's `bcf_get_info`), unpacking INFO.
    pub fn info_field(&mut self, header: &Header, tag: &[u8]) -> Option<&super::InfoField> {
        let id = header
            .id(tag)
            .filter(|&id| header.id_exists(HeaderRecordType::Info, id))?;
        self.info_field_by_id(id as i32)
    }

    /// Finds a decoded INFO field by ID (HTSlib's `bcf_get_info_id`), unpacking INFO.
    pub fn info_field_by_id(&mut self, id: i32) -> Option<&super::InfoField> {
        if self.unpacked & Unpack::INFO.bits() == 0 {
            self.unpack(Unpack::INFO);
        }
        self.d.info[..self.n_info as usize]
            .iter()
            .find(|f| f.key == id)
    }

    /// Finds a decoded FORMAT field by tag name (HTSlib's `bcf_get_fmt`), unpacking FORMAT.
    pub fn format_field(&mut self, header: &Header, tag: &[u8]) -> Option<&super::FormatField> {
        let id = header
            .id(tag)
            .filter(|&id| header.id_exists(HeaderRecordType::Format, id))?;
        self.format_field_by_id(id as i32)
    }

    /// Finds a decoded FORMAT field by ID (HTSlib's `bcf_get_fmt_id`), unpacking FORMAT.
    pub fn format_field_by_id(&mut self, id: i32) -> Option<&super::FormatField> {
        if self.unpacked & Unpack::FMT.bits() == 0 {
            self.unpack(Unpack::FMT);
        }
        let n = (self.n_fmt as usize).min(self.d.fmt.len());
        self.d.fmt[..n].iter().find(|f| f.id == id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::{bcf, encoding::FLOAT_MISSING, vcf};

    /// A header with tags `I` and `F` (INFO and FORMAT) and INFO `S` of these types.
    fn header(i: &str, f: &str, s: &str) -> Header {
        let text = format!(
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=I,Number=.,Type={i},Description=\"i\">\n\
             ##INFO=<ID=F,Number=.,Type={f},Description=\"f\">\n\
             ##INFO=<ID=S,Number=1,Type={s},Description=\"s\">\n\
             ##FORMAT=<ID=I,Number=.,Type={i},Description=\"i\">\n\
             ##FORMAT=<ID=F,Number=.,Type={f},Description=\"f\">\n\
             ##contig=<ID=1>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        );
        Header::parse(text.as_bytes()).expect("the header parses")
    }

    /// `line` written with `I` as integers, `F` as floats, and `S` as a string, then read
    /// back with the header types changed (as `bcftools reheader` can): as decoded from VCF,
    /// and as unpacked from BCF.
    fn records(line: &str) -> (Header, [Record; 2]) {
        let mut written = header("Integer", "Float", "String");
        let mut text = line.as_bytes().to_vec();
        let mut parsed = Record::new();
        vcf::parse(&mut written, &mut text, &mut parsed).expect("the line parses");
        let mut packed = Vec::new();
        bcf::pack_record(&mut parsed.clone(), &mut packed);
        let mut unpacked = Record::new();
        bcf::unpack_record(&packed, &mut unpacked).expect("the record unpacks");
        let read = header("Float", "Integer", "Integer");
        for tag in [&b"I"[..], b"F", b"S"] {
            assert_eq!(read.id(tag), written.id(tag), "{tag:?}");
        }
        (read, [parsed, unpacked])
    }

    fn bits(v: &[f32]) -> Vec<u32> {
        v.iter().map(|f| f.to_bits()).collect()
    }

    #[test]
    fn values_stored_as_another_type_keep_their_bits() {
        let (h, recs) = records("1\t5\t.\tA\tT\t.\t.\tI=5,.,-3;F=1.5,.\tI:F\t1,.:2.5\t7:.");
        let word = |v: i32| v as u32;
        for mut rec in recs {
            let mut f = Vec::new();
            assert_eq!(rec.get_info_f32(&h, b"I", &mut f), Ok(3));
            assert_eq!(bits(&f), [5, word(INT32_MISSING), word(-3)]);
            let mut i = Vec::new();
            assert_eq!(rec.get_info_i32(&h, b"F", &mut i), Ok(2));
            assert_eq!(i, [1.5f32.to_bits() as i32, FLOAT_MISSING as i32]);
            assert_eq!(rec.get_format_f32(&h, b"I", &mut f), Ok(4));
            let ends = [1, word(INT32_MISSING), 7, word(INT32_VECTOR_END)];
            assert_eq!(bits(&f), ends);
            assert_eq!(rec.get_format_i32(&h, b"F", &mut i), Ok(2));
            assert_eq!(i, [2.5f32.to_bits() as i32, FLOAT_MISSING as i32]);
        }
    }

    #[test]
    fn values_that_cannot_be_read_leave_the_buffer() {
        // `I` without a value is stored with no type, which HTSlib reports.
        let (h, recs) = records("1\t5\t.\tA\tT\t.\t.\tI;S=ab\tF\t1\t2");
        for mut rec in recs {
            let mut f = vec![9.0];
            let (got, messages) = log::capture(|| rec.get_info_f32(&h, b"I", &mut f));
            assert_eq!(got, Err(GetError::WrongType));
            assert_eq!(messages.len(), 1);
            assert_eq!(bits(&f), [9f32.to_bits()]);
            let mut i = vec![9];
            assert_eq!(rec.get_info_i32(&h, b"F", &mut i), Err(GetError::Absent));
            let (got, _) = log::capture(|| rec.get_info_i32(&h, b"S", &mut i));
            assert_eq!(got, Err(GetError::WrongType));
            assert_eq!(i, [9]);
            assert_eq!(rec.get_format_f32(&h, b"I", &mut f), Err(GetError::Absent));
            assert_eq!(bits(&f), [9f32.to_bits()]);
        }
    }
}
