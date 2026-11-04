//! Editing records: INFO, FORMAT, FILTER, alleles, and ID (HTSlib's `bcf_update_*`,
//! `bcf_add_filter`, `bcf_remove_filter`, `bcf_add_id`).

use std::sync::atomic::AtomicBool;

use super::{
    DIRTY_ALS, DIRTY_FLT, DIRTY_ID, DIRTY_INF, FieldData, FormatField, InfoField, Record,
    RecordError, Unpack, rlen::get_rlen, unpack_fmt_core1, unpack_info_core1,
};
use crate::variant::{
    Header, HeaderRecordType,
    encoding::{self, BcfType},
    header::SpecialTag,
    log,
};

/// An error editing a record.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum UpdateError {
    /// The tag is not defined in the header for this kind of field.
    NotInHeader,
    /// END must be a single integer.
    InvalidEnd,
    /// The number of values is not a multiple of the number of samples.
    NotPerSample,
}

impl std::fmt::Display for UpdateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NotInHeader => "the tag is not defined in the header",
            Self::InvalidEnd => "END must be a single integer",
            Self::NotPerSample => "the number of values is not a multiple of the sample count",
        })
    }
}

impl std::error::Error for UpdateError {}

static NEGATIVE_RLEN_WARNED: AtomicBool = AtomicBool::new(false);

/// INFO values to store; see [`Record::update_info`].
#[derive(Clone, Copy, Debug)]
pub(crate) enum InfoUpdate<'a> {
    Remove,
    Int(&'a [i32]),
    Float(&'a [f32]),
    Str(&'a [u8]),
    Flag,
}

impl Record {
    /// Sets an integer INFO field (HTSlib's `bcf_update_info_int32`); an empty slice removes it.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an INFO tag of the header, or for an END that is
    /// not a single value.
    pub fn update_info_i32(
        &mut self,
        header: &Header,
        key: &[u8],
        values: &[i32],
    ) -> Result<(), UpdateError> {
        let update = if values.is_empty() {
            InfoUpdate::Remove
        } else {
            InfoUpdate::Int(values)
        };
        self.update_info(header, key, update)
    }

    /// Sets a float INFO field (HTSlib's `bcf_update_info_float`); an empty slice removes it.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an INFO tag of the header, or for END.
    pub fn update_info_f32(
        &mut self,
        header: &Header,
        key: &[u8],
        values: &[f32],
    ) -> Result<(), UpdateError> {
        let update = if values.is_empty() {
            InfoUpdate::Remove
        } else {
            InfoUpdate::Float(values)
        };
        self.update_info(header, key, update)
    }

    /// Sets a string INFO field (HTSlib's `bcf_update_info_string`); `None` removes it.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an INFO tag of the header, or for END.
    pub fn update_info_string(
        &mut self,
        header: &Header,
        key: &[u8],
        value: Option<&[u8]>,
    ) -> Result<(), UpdateError> {
        let update = value.map_or(InfoUpdate::Remove, InfoUpdate::Str);
        self.update_info(header, key, update)
    }

    /// Sets or removes an INFO flag (HTSlib's `bcf_update_info_flag`).
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an INFO tag of the header, or for END.
    pub fn update_info_flag(
        &mut self,
        header: &Header,
        key: &[u8],
        set: bool,
    ) -> Result<(), UpdateError> {
        let update = if set {
            InfoUpdate::Flag
        } else {
            InfoUpdate::Remove
        };
        self.update_info(header, key, update)
    }

    /// Removes an INFO field.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not an INFO tag of the header.
    pub fn remove_info(&mut self, header: &Header, key: &[u8]) -> Result<(), UpdateError> {
        self.update_info(header, key, InfoUpdate::Remove)
    }

    /// [`Self::update_info_i32`] for a tag given by its dictionary ID ([`Header::id`]), which
    /// saves the name lookup when the same tag is set in many records.
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_info_i32`] does.
    pub fn update_info_i32_id(
        &mut self,
        header: &Header,
        id: u32,
        values: &[i32],
    ) -> Result<(), UpdateError> {
        let update = if values.is_empty() {
            InfoUpdate::Remove
        } else {
            InfoUpdate::Int(values)
        };
        self.update_info_id(header, id, update)
    }

    /// [`Self::update_info_f32`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_info_f32`] does.
    pub fn update_info_f32_id(
        &mut self,
        header: &Header,
        id: u32,
        values: &[f32],
    ) -> Result<(), UpdateError> {
        let update = if values.is_empty() {
            InfoUpdate::Remove
        } else {
            InfoUpdate::Float(values)
        };
        self.update_info_id(header, id, update)
    }

    /// [`Self::remove_info`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Returns an error when the ID is not an INFO tag of the header.
    pub fn remove_info_id(&mut self, header: &Header, id: u32) -> Result<(), UpdateError> {
        self.update_info_id(header, id, InfoUpdate::Remove)
    }

    /// Appends a copy of another record's INFO field, without decoding its values, under
    /// dictionary ID `key` of `header` (bcftools merge's copy of a `bcf_info_t`). Unlike
    /// [`Self::update_info_i32`], an existing field with the same key is not replaced. END and
    /// SVLEN update the record's length.
    pub fn push_info_copy(&mut self, header: &Header, src: &Record, field: &InfoField, key: u32) {
        if self.unpacked & Unpack::INFO.bits() == 0 {
            self.unpack(Unpack::INFO);
        }
        let mut s = self.d.take_buf();
        encoding::enc_int1(&mut s, key as i32);
        encoding::enc_size(&mut s, field.len(), field.ty);
        s.extend_from_slice(src.info_data(field));
        let mut f = InfoField::default();
        unpack_info_core1(&s, 0, &mut f);
        f.data = FieldData::Owned(s);
        self.d.info.truncate(self.n_info as usize);
        self.d.info.push(f);
        self.n_info += 1;
        self.d.shared_dirty |= DIRTY_INF;
        self.unpacked |= Unpack::INFO.bits();
        if header.special_id(SpecialTag::End) == Some(key)
            || header.special_id(SpecialTag::Svlen) == Some(key)
        {
            self.rlen = get_rlen(header, self);
        }
    }

    /// Marks every INFO field for removal in one pass (what bcftools annotate `-x INFO` does),
    /// updating `rlen` when END or SVLEN was among them.
    pub fn remove_all_info(&mut self, header: &Header) {
        let end = header.special_id(SpecialTag::End).map(|i| i as i32);
        let svlen = header.special_id(SpecialTag::Svlen).map(|i| i as i32);
        let mut rlen_tag = false;
        if self.unpacked & Unpack::INFO.bits() == 0
            && self.d.shared_dirty == 0
            && !self.shared.is_empty()
        {
            // Cut INFO off the encoded block rather than decode it, reading only the keys.
            let start = self.encoded_info_offset();
            let mut p = start;
            for _ in 0..self.n_info {
                let (key, used) = encoding::dec_typed_int1(&self.shared[p..]);
                rlen_tag |= Some(key as i32) == end || Some(key as i32) == svlen;
                let (n, ty, size) = encoding::dec_size(&self.shared[p + used..]);
                p += used + size + ((n.max(0) as usize) << encoding::type_shift(ty));
            }
            self.shared.truncate(start);
            self.n_info = 0;
            if rlen_tag {
                self.rlen = get_rlen(header, self);
            }
            return;
        }
        if self.unpacked & Unpack::INFO.bits() == 0 {
            self.unpack(Unpack::INFO);
        }
        for i in 0..self.n_info as usize {
            let inf = &mut self.d.info[i];
            if inf.is_removed() {
                continue;
            }
            rlen_tag |= Some(inf.key) == end || Some(inf.key) == svlen;
            let old = std::mem::replace(&mut inf.data, FieldData::Removed);
            inf.vptr_off = 0;
            inf.vptr_len = 0;
            self.d.recycle(old);
            self.d.shared_dirty |= DIRTY_INF;
        }
        if rlen_tag {
            self.rlen = get_rlen(header, self);
        }
    }

    /// Marks every FORMAT field but `keep` (a dictionary ID, such as GT's) for removal in one
    /// pass (what bcftools annotate `-x FORMAT` does).
    pub fn remove_all_format(&mut self, header: &Header, keep: Option<u32>) {
        if self.unpacked & Unpack::FMT.bits() == 0 {
            self.unpack(Unpack::FMT);
        }
        let len = header.special_id(SpecialTag::Len).map(|i| i as i32);
        let keep = keep.map(|i| i as i32);
        let mut rlen_tag = false;
        let n_fmt = (self.n_fmt as usize).min(self.d.fmt.len());
        for i in 0..n_fmt {
            let f = &mut self.d.fmt[i];
            if f.is_removed() || Some(f.id) == keep {
                continue;
            }
            rlen_tag |= Some(f.id) == len;
            let old = std::mem::replace(&mut f.data, FieldData::Removed);
            self.d.recycle(old);
            self.d.indiv_dirty = true;
        }
        if rlen_tag {
            self.rlen = get_rlen(header, self);
        }
    }

    /// HTSlib's `bcf_update_info`.
    pub(crate) fn update_info(
        &mut self,
        header: &Header,
        key: &[u8],
        update: InfoUpdate<'_>,
    ) -> Result<(), UpdateError> {
        let id = header.id(key).ok_or(UpdateError::NotInHeader)?;
        self.update_info_id(header, id, update)
    }

    /// [`Self::update_info`] for a dictionary ID.
    fn update_info_id(
        &mut self,
        header: &Header,
        inf_id: u32,
        update: InfoUpdate<'_>,
    ) -> Result<(), UpdateError> {
        if !header.id_exists(HeaderRecordType::Info, inf_id) {
            return Err(UpdateError::NotInHeader);
        }
        if self.unpacked & Unpack::INFO.bits() == 0 {
            self.unpack(Unpack::INFO);
        }
        let is_end = header.special_id(SpecialTag::End) == Some(inf_id);
        let is_svlen = header.special_id(SpecialTag::Svlen) == Some(inf_id);
        let found = self.d.info[..self.n_info as usize]
            .iter()
            .position(|f| f.key == inf_id as i32);

        if matches!(update, InfoUpdate::Remove) {
            if let Some(i) = found {
                let inf = &mut self.d.info[i];
                let old = std::mem::replace(&mut inf.data, FieldData::Removed);
                inf.vptr_off = 0;
                inf.vptr_len = 0;
                self.d.recycle(old);
                self.d.shared_dirty |= DIRTY_INF;
            }
            if is_end || is_svlen {
                self.rlen = get_rlen(header, self);
            }
            return Ok(());
        }

        if is_end {
            let ok_type = matches!(update, InfoUpdate::Int(_));
            let n = match update {
                InfoUpdate::Int(v) => v.len(),
                InfoUpdate::Float(v) => v.len(),
                _ => 1,
            };
            if n != 1 {
                log::error(
                    "bcf_update_info",
                    format_args!(
                        "END info tag should only have one value at {}:{}",
                        String::from_utf8_lossy(self.chrom_safe(header)),
                        self.pos + 1
                    ),
                );
                self.errcode |= RecordError::TAG_INVALID;
                return Err(UpdateError::InvalidEnd);
            }
            if !ok_type {
                let ty = match update {
                    InfoUpdate::Float(_) => 2,
                    InfoUpdate::Flag => 0,
                    _ => 3,
                };
                log::error(
                    "bcf_update_info",
                    format_args!(
                        "Wrong type ({ty}) for END info tag at {}:{}",
                        String::from_utf8_lossy(self.chrom_safe(header)),
                        self.pos + 1
                    ),
                );
                self.errcode |= RecordError::TAG_INVALID;
                return Err(UpdateError::InvalidEnd);
            }
        }

        let mut s = self.d.take_buf();
        encoding::enc_int1(&mut s, inf_id as i32);
        match update {
            InfoUpdate::Int(v) => encoding::enc_vint(&mut s, v, None),
            InfoUpdate::Float(v) => encoding::enc_vfloat(&mut s, v),
            InfoUpdate::Str(v) => {
                let v = &v[..v.iter().position(|&c| c == 0).unwrap_or(v.len())];
                encoding::enc_vchar(&mut s, v);
            }
            InfoUpdate::Flag => encoding::enc_size(&mut s, 0, BcfType::Null as u8),
            InfoUpdate::Remove => unreachable!("removal is handled above"),
        }

        match found {
            Some(i) => {
                let fits = {
                    let inf = &self.d.info[i];
                    !inf.is_removed() && s.len() <= (inf.vptr_len + inf.vptr_off) as usize
                };
                if fits {
                    let inf = &self.d.info[i];
                    if s.len() != (inf.vptr_len + inf.vptr_off) as usize {
                        self.d.shared_dirty |= DIRTY_INF;
                    }
                    let data = std::mem::take(&mut self.d.info[i].data);
                    match data {
                        FieldData::InBlock(start) => {
                            self.shared[start..start + s.len()].copy_from_slice(&s);
                            let mut f = InfoField::default();
                            unpack_info_core1(&self.shared, start, &mut f);
                            self.d.info[i] = f;
                        }
                        FieldData::Owned(mut buf) => {
                            buf[..s.len()].copy_from_slice(&s);
                            let mut f = InfoField::default();
                            unpack_info_core1(&buf, 0, &mut f);
                            f.data = FieldData::Owned(buf);
                            self.d.info[i] = f;
                        }
                        FieldData::Removed => unreachable!("checked above"),
                    }
                    self.d.put_buf(s);
                } else {
                    let mut f = InfoField::default();
                    unpack_info_core1(&s, 0, &mut f);
                    f.data = FieldData::Owned(s);
                    let old = std::mem::replace(&mut self.d.info[i], f);
                    self.d.recycle(old.data);
                    self.d.shared_dirty |= DIRTY_INF;
                }
            }
            None => {
                let mut f = InfoField::default();
                unpack_info_core1(&s, 0, &mut f);
                f.data = FieldData::Owned(s);
                self.d.info.truncate(self.n_info as usize);
                self.d.info.push(f);
                self.n_info += 1;
                self.d.shared_dirty |= DIRTY_INF;
            }
        }
        self.unpacked |= Unpack::INFO.bits();

        if is_end
            && let InfoUpdate::Int(&[end]) = update
            && end != encoding::INT32_MISSING
            && i64::from(end) <= self.pos
        {
            log::warn_once(
                &NEGATIVE_RLEN_WARNED,
                "bcf_update_info",
                format_args!(
                    "INFO/END={end} is smaller than POS at {}:{}",
                    String::from_utf8_lossy(self.chrom_safe(header)),
                    self.pos + 1
                ),
            );
        }
        if is_end || is_svlen {
            self.rlen = get_rlen(header, self);
        }
        Ok(())
    }

    /// Sets an integer FORMAT field for all samples (HTSlib's `bcf_update_format_int32`);
    /// `values` holds the same number of values for each sample, padded with
    /// [`encoding::INT32_VECTOR_END`]. An empty slice removes the field.
    ///
    /// # Errors
    ///
    /// Returns an error when the tag is not a FORMAT tag of the header, or the number of
    /// values is not a multiple of the number of samples.
    pub fn update_format_i32(
        &mut self,
        header: &Header,
        key: &[u8],
        values: &[i32],
    ) -> Result<(), UpdateError> {
        self.update_format(header, key, FormatUpdate::Int(values))
    }

    /// Sets a float FORMAT field for all samples (HTSlib's `bcf_update_format_float`).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_format_i32`] does.
    pub fn update_format_f32(
        &mut self,
        header: &Header,
        key: &[u8],
        values: &[f32],
    ) -> Result<(), UpdateError> {
        self.update_format(header, key, FormatUpdate::Float(values))
    }

    /// Sets a fixed-width character FORMAT field: `values` holds the same number of bytes
    /// for each sample, NUL-padded (HTSlib's `bcf_update_format_char`).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_format_i32`] does.
    pub fn update_format_chars(
        &mut self,
        header: &Header,
        key: &[u8],
        values: &[u8],
    ) -> Result<(), UpdateError> {
        self.update_format(header, key, FormatUpdate::Char(values))
    }

    /// Sets a string FORMAT field from one string per sample, padding them to the longest
    /// (HTSlib's `bcf_update_format_string`). An empty slice removes the field.
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_format_i32`] does.
    pub fn update_format_strings(
        &mut self,
        header: &Header,
        key: &[u8],
        values: &[&[u8]],
    ) -> Result<(), UpdateError> {
        if values.is_empty() {
            return self.update_format(header, key, FormatUpdate::Char(&[]));
        }
        let max_len = values.iter().map(|v| v.len()).max().unwrap_or(0);
        let mut out = vec![0u8; max_len * values.len()];
        for (i, v) in values.iter().enumerate() {
            let v = &v[..v.iter().position(|&c| c == 0).unwrap_or(v.len())];
            out[i * max_len..i * max_len + v.len()].copy_from_slice(v);
        }
        self.update_format(header, key, FormatUpdate::Char(&out))
    }

    /// Sets genotypes, encoded as [`crate::variant::Genotype`] values with vector-end padding
    /// (HTSlib's `bcf_update_genotypes`).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_format_i32`] does.
    pub fn update_genotypes(&mut self, header: &Header, values: &[i32]) -> Result<(), UpdateError> {
        match header.special_id(SpecialTag::Gt) {
            Some(id) => self.update_format_id(header, id, FormatUpdate::Int(values)),
            None => self.update_format(header, b"GT", FormatUpdate::Int(values)),
        }
    }

    /// Removes a FORMAT field.
    ///
    /// # Errors
    ///
    /// Never fails for an undefined tag (as in HTSlib).
    pub fn remove_format(&mut self, header: &Header, key: &[u8]) -> Result<(), UpdateError> {
        self.update_format(header, key, FormatUpdate::Int(&[]))
    }

    /// [`Self::update_format_i32`] for a tag given by its dictionary ID ([`Header::id`]),
    /// which saves the name lookup when the same tag is set in many records.
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_format_i32`] does.
    pub fn update_format_i32_id(
        &mut self,
        header: &Header,
        id: u32,
        values: &[i32],
    ) -> Result<(), UpdateError> {
        self.update_format_id(header, id, FormatUpdate::Int(values))
    }

    /// [`Self::update_format_f32`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Returns an error as [`Self::update_format_i32`] does.
    pub fn update_format_f32_id(
        &mut self,
        header: &Header,
        id: u32,
        values: &[f32],
    ) -> Result<(), UpdateError> {
        self.update_format_id(header, id, FormatUpdate::Float(values))
    }

    /// [`Self::remove_format`] for a tag given by its dictionary ID ([`Header::id`]).
    ///
    /// # Errors
    ///
    /// Never fails for an undefined tag (as in HTSlib).
    pub fn remove_format_id(&mut self, header: &Header, id: u32) -> Result<(), UpdateError> {
        self.update_format_id(header, id, FormatUpdate::Int(&[]))
    }

    /// HTSlib's `bcf_update_format`.
    fn update_format(
        &mut self,
        header: &Header,
        key: &[u8],
        update: FormatUpdate<'_>,
    ) -> Result<(), UpdateError> {
        match header.id(key) {
            Some(id) => self.update_format_id(header, id, update),
            None if update.len() == 0 => Ok(()),
            None => Err(UpdateError::NotInHeader),
        }
    }

    /// [`Self::update_format`] for a dictionary ID.
    fn update_format_id(
        &mut self,
        header: &Header,
        fmt_id: u32,
        update: FormatUpdate<'_>,
    ) -> Result<(), UpdateError> {
        let n = update.len();
        if !header.id_exists(HeaderRecordType::Format, fmt_id) {
            return if n == 0 {
                Ok(())
            } else {
                Err(UpdateError::NotInHeader)
            };
        }
        if self.unpacked & Unpack::FMT.bits() == 0 {
            self.unpack(Unpack::FMT);
        }
        let n_fmt = (self.n_fmt as usize).min(self.d.fmt.len());
        let found = self.d.fmt[..n_fmt]
            .iter()
            .position(|f| f.id == fmt_id as i32);
        let is_len = header.special_id(SpecialTag::Len) == Some(fmt_id);
        if n == 0 {
            if let Some(i) = found {
                let old = std::mem::replace(&mut self.d.fmt[i].data, FieldData::Removed);
                self.d.recycle(old);
                self.d.indiv_dirty = true;
            }
            if is_len {
                self.rlen = get_rlen(header, self);
            }
            return Ok(());
        }
        self.n_sample = header.sample_count() as u32;
        let n_sample = self.n_sample as usize;
        if n_sample == 0 || !n.is_multiple_of(n_sample) {
            return Err(UpdateError::NotPerSample);
        }
        let nps = n / n_sample;
        let mut s = self.d.take_buf();
        s.reserve(n * 4 + 8);
        encoding::enc_int1(&mut s, fmt_id as i32);
        match update {
            FormatUpdate::Int(v) => encoding::enc_vint(&mut s, v, Some(nps)),
            FormatUpdate::Float(v) => {
                encoding::enc_size(&mut s, nps, BcfType::Float as u8);
                encoding::push_floats(&mut s, v);
            }
            FormatUpdate::Char(v) => {
                encoding::enc_size(&mut s, nps, BcfType::Char as u8);
                s.extend_from_slice(v);
            }
        }
        match found {
            None => {
                let mut f = FormatField::default();
                unpack_fmt_core1(&s, 0, n_sample, &mut f);
                f.data = FieldData::Owned(s);
                self.d.fmt.truncate(n_fmt);
                if n_fmt >= 1 && header.special_id(SpecialTag::Gt) == Some(fmt_id) {
                    self.d.fmt.insert(0, f);
                } else {
                    self.d.fmt.push(f);
                }
                self.n_fmt += 1;
                self.d.indiv_dirty = true;
            }
            Some(i) => {
                let fits = {
                    let f = &self.d.fmt[i];
                    !f.is_removed() && s.len() <= (f.p_len + f.p_off) as usize
                };
                if fits {
                    let f = &self.d.fmt[i];
                    if s.len() != (f.p_len + f.p_off) as usize {
                        self.d.indiv_dirty = true;
                    }
                    let data = std::mem::take(&mut self.d.fmt[i].data);
                    match data {
                        FieldData::InBlock(start) => {
                            self.indiv[start..start + s.len()].copy_from_slice(&s);
                            let mut f = FormatField::default();
                            unpack_fmt_core1(&self.indiv, start, n_sample, &mut f);
                            self.d.fmt[i] = f;
                        }
                        FieldData::Owned(mut buf) => {
                            buf[..s.len()].copy_from_slice(&s);
                            let mut f = FormatField::default();
                            unpack_fmt_core1(&buf, 0, n_sample, &mut f);
                            f.data = FieldData::Owned(buf);
                            self.d.fmt[i] = f;
                        }
                        FieldData::Removed => unreachable!("checked above"),
                    }
                    self.d.put_buf(s);
                } else {
                    let mut f = FormatField::default();
                    unpack_fmt_core1(&s, 0, n_sample, &mut f);
                    f.data = FieldData::Owned(s);
                    let old = std::mem::replace(&mut self.d.fmt[i], f);
                    self.d.recycle(old.data);
                    self.d.indiv_dirty = true;
                }
            }
        }
        self.unpacked |= Unpack::FMT.bits();
        if is_len {
            self.rlen = get_rlen(header, self);
        }
        Ok(())
    }

    /// Replaces the FILTER IDs (HTSlib's `bcf_update_filter`); empty means `.`.
    pub fn update_filter(&mut self, ids: &[i32]) {
        if self.unpacked & Unpack::FLT.bits() == 0 {
            self.unpack(Unpack::FLT);
        }
        self.d.shared_dirty |= DIRTY_FLT;
        self.d.flt.clear();
        self.d.flt.extend_from_slice(ids);
    }

    /// Adds a FILTER ID (HTSlib's `bcf_add_filter`): `PASS` (0) replaces all others, and
    /// another filter replaces `PASS`. Returns whether the filter was added.
    pub fn add_filter(&mut self, id: i32) -> bool {
        if self.unpacked & Unpack::FLT.bits() == 0 {
            self.unpack(Unpack::FLT);
        }
        if self.d.flt.contains(&id) {
            return false;
        }
        self.d.shared_dirty |= DIRTY_FLT;
        if id == 0 || (self.d.flt.len() == 1 && self.d.flt[0] == 0) {
            self.d.flt.clear();
        }
        self.d.flt.push(id);
        true
    }

    /// Removes a FILTER ID (HTSlib's `bcf_remove_filter`); with [`PassOnEmpty::Pass`], an
    /// emptied FILTER becomes `PASS`.
    pub fn remove_filter(&mut self, id: i32, pass: PassOnEmpty) {
        if self.unpacked & Unpack::FLT.bits() == 0 {
            self.unpack(Unpack::FLT);
        }
        let Some(i) = self.d.flt.iter().position(|&f| f == id) else {
            return;
        };
        self.d.shared_dirty |= DIRTY_FLT;
        self.d.flt.remove(i);
        if self.d.flt.is_empty() && pass == PassOnEmpty::Pass {
            self.add_filter(0);
        }
    }

    /// Returns whether the record has a filter (`.` means `PASS`, which an empty FILTER also
    /// counts as), or `None` when the filter is not defined in the header (HTSlib's
    /// `bcf_has_filter`).
    pub fn has_filter(&mut self, header: &Header, name: &[u8]) -> Option<bool> {
        let name: &[u8] = if name == b"." { b"PASS" } else { name };
        let id = header
            .id(name)
            .filter(|&id| header.id_exists(HeaderRecordType::Filter, id))?;
        if self.unpacked & Unpack::FLT.bits() == 0 {
            self.unpack(Unpack::FLT);
        }
        if id == 0 && self.d.flt.is_empty() {
            return Some(true);
        }
        Some(self.d.flt.contains(&(id as i32)))
    }

    /// Replaces the alleles, REF first (HTSlib's `bcf_update_alleles`), and updates `rlen`.
    pub fn update_alleles<A: AsRef<[u8]>>(&mut self, header: &Header, alleles: &[A]) {
        if self.unpacked & Unpack::STR.bits() == 0 {
            self.unpack(Unpack::STR);
        }
        self.d.als.clear();
        self.d.allele_starts.clear();
        for a in alleles {
            let a = a.as_ref();
            let a = &a[..a.iter().position(|&c| c == 0).unwrap_or(a.len())];
            self.d.allele_starts.push(self.d.als.len() as u32);
            self.d.als.extend_from_slice(a);
            self.d.als.push(0);
        }
        self.n_allele = alleles.len() as u32;
        self.after_alleles_change(header);
    }

    /// Replaces the alleles from a comma-separated string (HTSlib's `bcf_update_alleles_str`).
    pub fn update_alleles_str(&mut self, header: &Header, alleles: &[u8]) {
        if self.unpacked & Unpack::STR.bits() == 0 {
            self.unpack(Unpack::STR);
        }
        let alleles = &alleles[..alleles
            .iter()
            .position(|&c| c == 0)
            .unwrap_or(alleles.len())];
        self.d.als.clear();
        self.d.allele_starts.clear();
        self.d.allele_starts.push(0);
        for &c in alleles {
            if c == b',' {
                self.d.als.push(0);
                self.d.allele_starts.push(self.d.als.len() as u32);
            } else {
                self.d.als.push(c);
            }
        }
        self.d.als.push(0);
        self.n_allele = self.d.allele_starts.len() as u32;
        self.after_alleles_change(header);
    }

    fn after_alleles_change(&mut self, header: &Header) {
        self.d.shared_dirty |= DIRTY_ALS;
        self.d.var_type = -1;
        self.rlen = get_rlen(header, self);
    }

    /// Replaces the ID (HTSlib's `bcf_update_id`); `None` sets `.`.
    pub fn update_id(&mut self, id: Option<&[u8]>) {
        if self.unpacked & Unpack::STR.bits() == 0 {
            self.unpack(Unpack::STR);
        }
        self.d.id.clear();
        let id = id.unwrap_or(b".");
        self.d
            .id
            .extend_from_slice(&id[..id.iter().position(|&c| c == 0).unwrap_or(id.len())]);
        self.d.shared_dirty |= DIRTY_ID;
    }

    /// Adds an ID to the `;`-separated list unless present (HTSlib's `bcf_add_id`).
    pub fn add_id(&mut self, id: &[u8]) {
        if self.unpacked & Unpack::STR.bits() == 0 {
            self.unpack(Unpack::STR);
        }
        // HTSlib looks for `id` as a whole `;`-delimited run, so an ID list like `a;b` is
        // found too.
        let hay = &self.d.id;
        let present = !id.is_empty()
            && memchr::memmem::find_iter(hay, id).any(|at| {
                let end = at + id.len();
                (at == 0 || hay[at - 1] == b';') && (end == hay.len() || hay[end] == b';')
            });
        if present {
            return;
        }
        if self.d.id.as_slice() == b"." || self.d.id.is_empty() {
            self.d.id.clear();
        } else {
            self.d.id.push(b';');
        }
        self.d.id.extend_from_slice(id);
        self.d.shared_dirty |= DIRTY_ID;
    }
}

/// Whether removing the last filter sets `PASS`.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum PassOnEmpty {
    /// Leave FILTER empty (`.`).
    Empty,
    /// Set `PASS`.
    Pass,
}

#[derive(Clone, Copy, Debug)]
enum FormatUpdate<'a> {
    Int(&'a [i32]),
    Float(&'a [f32]),
    Char(&'a [u8]),
}

impl FormatUpdate<'_> {
    fn len(&self) -> usize {
        match self {
            Self::Int(v) => v.len(),
            Self::Float(v) => v.len(),
            Self::Char(v) => v.len(),
        }
    }
}
