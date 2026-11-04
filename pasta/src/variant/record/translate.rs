//! Translating a record's IDs between headers (HTSlib's `bcf_translate`).

use super::{DIRTY_FLT, DIRTY_INF, FieldData, Record, RecordError, Unpack};
use crate::variant::{
    Header,
    encoding::{enc_int1, enc_size},
    header::Translation,
    log,
};

/// The encoded width class of an ID, as `bcf_translate` decides whether IDs can be rewritten
/// in place.
fn id_width(id: i32) -> u8 {
    if id >> 7 == 0 {
        1
    } else if id >> 15 == 0 {
        2
    } else {
        4
    }
}

/// Rewrites the ID of an encoded field block in place (the key starts at byte 1).
fn write_id(block: &mut [u8], id: i32, width: u8) {
    match width {
        1 => block[1] = id as u8,
        2 => block[1..3].copy_from_slice(&(id as i16).to_le_bytes()),
        _ => block[1..5].copy_from_slice(&id.to_le_bytes()),
    }
}

impl Record {
    /// Rewrites contig, FILTER, INFO, and FORMAT IDs from `src`'s dictionaries to `dst`'s
    /// (HTSlib's `bcf_translate`). The translation tables are computed once per destination
    /// and cached in `src`; another destination, or a change to either header's dictionaries,
    /// builds them again (HTSlib keeps those of the first destination).
    ///
    /// HTSlib writes INFO IDs of 128 and more at the wrong offset when rewriting them in place,
    /// corrupting BCF output; this writes them after the type byte, as it does for FORMAT.
    ///
    /// # Errors
    ///
    /// Returns the record's unchecked parse errors, which HTSlib treats as fatal.
    pub fn translate(&mut self, dst: &Header, src: &mut Header) -> Result<(), RecordError> {
        if !self.errcode.is_empty() {
            log::error(
                "bcf_translate",
                format_args!(
                    "Unchecked error ({} {}) at {}:{}, exiting",
                    self.errcode.bits(),
                    self.errcode.describe(1024),
                    String::from_utf8_lossy(self.chrom_safe(src)),
                    self.pos + 1
                ),
            );
            return Err(self.errcode);
        }
        let dst_stamp = dst.dictionaries.stamp();
        if !src.translation.is_for(dst_stamp) {
            let mut changed = 0;
            let ids: Vec<i32> = (0..src.id_count() as u32)
                .map(|i| match src.id_name(i) {
                    None => -1,
                    Some(name) => {
                        let t = dst.id(name).map_or(-1, |t| t as i32);
                        if t != -1 && t != i as i32 {
                            changed += 1;
                        }
                        t
                    }
                })
                .collect();
            let ctgs: Vec<i32> = (0..src.contig_count() as i32)
                .map(|i| match src.contig_name(i) {
                    None => -1,
                    Some(name) => {
                        let t = dst.contig_id(name).map_or(-1, |t| t as i32);
                        if t != -1 && t != i {
                            changed += 1;
                        }
                        t
                    }
                })
                .collect();
            src.translation = if changed == 0 {
                Translation::Identity { dst: dst_stamp }
            } else {
                Translation::Map {
                    dst: dst_stamp,
                    ids,
                    ctgs,
                }
            };
        }
        let Translation::Map { ids, ctgs, .. } = &src.translation else {
            return Ok(());
        };
        self.unpack(Unpack::ALL);
        let map = |table: &[i32], id: i32| table.get(id as usize).copied().unwrap_or(-1);
        let t = map(ctgs, self.rid);
        if t >= 0 {
            self.rid = t;
        }
        for i in 0..self.d.flt.len() {
            let t = map(ids, self.d.flt[i]);
            if t >= 0 {
                self.d.flt[i] = t;
            }
            self.d.shared_dirty |= DIRTY_FLT;
        }
        for i in 0..self.n_info as usize {
            let dst_id = map(ids, self.d.info[i].key);
            if dst_id >= 0 {
                self.rekey_info(i, dst_id);
            }
        }
        for i in 0..(self.n_fmt as usize).min(self.d.fmt.len()) {
            let dst_id = map(ids, self.d.fmt[i].id);
            if dst_id >= 0 {
                self.rekey_format(i, dst_id);
            }
        }
        Ok(())
    }

    /// Changes the dictionary ID of every INFO field with ID `old` to `new`, re-encoding the
    /// key (bcftools annotate's `rename_annots_info`). The header must define `new`.
    pub fn rename_info_id(&mut self, old: i32, new: i32) {
        if self.unpacked & Unpack::INFO.bits() == 0 {
            self.unpack(Unpack::INFO);
        }
        for i in 0..self.n_info as usize {
            if self.d.info[i].key == old && !self.d.info[i].is_removed() {
                self.rekey_info(i, new);
            }
        }
    }

    /// Changes the dictionary ID of every FORMAT field with ID `old` to `new` (bcftools
    /// annotate's `rename_annots_format`).
    pub fn rename_format_id(&mut self, old: i32, new: i32) {
        if self.unpacked & Unpack::FMT.bits() == 0 {
            self.unpack(Unpack::FMT);
        }
        for i in 0..(self.n_fmt as usize).min(self.d.fmt.len()) {
            if self.d.fmt[i].id == old && !self.d.fmt[i].is_removed() {
                self.rekey_format(i, new);
            }
        }
    }

    /// Changes FILTER ID `old` to `new` (bcftools annotate's renaming of filters).
    pub fn rename_filter_id(&mut self, old: i32, new: i32) {
        if self.unpacked & Unpack::FLT.bits() == 0 {
            self.unpack(Unpack::FLT);
        }
        for f in &mut self.d.flt {
            if *f == old {
                *f = new;
            }
        }
        self.d.shared_dirty |= DIRTY_FLT;
    }

    /// Sets the ID of INFO field `i`, rewriting its encoded key.
    fn rekey_info(&mut self, i: usize, dst_id: i32) {
        let src_id = self.d.info[i].key;
        self.d.info[i].key = dst_id;
        if self.d.info[i].is_removed() {
            return;
        }
        let (sw, dw) = (id_width(src_id), id_width(dst_id));
        if sw == dw {
            match &mut self.d.info[i].data {
                FieldData::InBlock(start) => write_id(&mut self.shared[*start..], dst_id, dw),
                FieldData::Owned(buf) => write_id(buf, dst_id, dw),
                FieldData::Removed => {}
            }
        } else {
            let info = &self.d.info[i];
            let mut s = Vec::new();
            enc_int1(&mut s, dst_id);
            enc_size(&mut s, info.len.max(0) as usize, info.ty);
            let off = s.len() as u32;
            s.extend_from_slice(self.info_data(info));
            let info = &mut self.d.info[i];
            info.vptr_off = off;
            info.data = FieldData::Owned(s);
            self.d.shared_dirty |= DIRTY_INF;
        }
    }

    /// Sets the ID of FORMAT field `i`, rewriting its encoded key.
    fn rekey_format(&mut self, i: usize, dst_id: i32) {
        let src_id = self.d.fmt[i].id;
        self.d.fmt[i].id = dst_id;
        if self.d.fmt[i].is_removed() {
            return;
        }
        let (sw, dw) = (id_width(src_id), id_width(dst_id));
        if sw == dw {
            match &mut self.d.fmt[i].data {
                FieldData::InBlock(start) => write_id(&mut self.indiv[*start..], dst_id, dw),
                FieldData::Owned(buf) => write_id(buf, dst_id, dw),
                FieldData::Removed => {}
            }
        } else {
            let f = &self.d.fmt[i];
            let mut s = Vec::new();
            enc_int1(&mut s, dst_id);
            enc_size(&mut s, f.n.max(0) as usize, f.ty);
            let off = s.len() as u32;
            s.extend_from_slice(self.format_data(f));
            let f = &mut self.d.fmt[i];
            f.p_off = off;
            f.data = FieldData::Owned(s);
            self.d.indiv_dirty = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::variant::{Header, Record, Unpack, vcf};

    /// A header whose ID and contig dictionaries list `ids` and `contigs` in that order.
    fn header(ids: &[&str], contigs: &[&str]) -> Header {
        let mut text = String::from("##fileformat=VCFv4.2\n");
        for id in ids {
            text += &format!("##INFO=<ID={id},Number=1,Type=Integer,Description=\"{id}\">\n");
        }
        for contig in contigs {
            text += &format!("##contig=<ID={contig}>\n");
        }
        text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        Header::parse(text.as_bytes()).expect("the header parses")
    }

    fn record(src: &mut Header) -> Record {
        let mut line = b"c2\t5\t.\tA\tG\t.\t.\tA=1;B=2".to_vec();
        let mut rec = Record::new();
        vcf::parse(src, &mut line, &mut rec).expect("the line parses");
        rec
    }

    /// The record's contig and INFO keys, by name in `header`.
    fn names(rec: &mut Record, header: &Header) -> (Vec<u8>, Vec<Vec<u8>>) {
        rec.unpack(Unpack::ALL);
        let keys = rec
            .info_fields()
            .iter()
            .map(|f| header.id_name(f.key as u32).expect("a known ID").to_vec())
            .collect();
        (rec.chrom_safe(header).to_vec(), keys)
    }

    #[test]
    fn tables_are_rebuilt_for_another_destination_and_after_a_clone() {
        let mut src = header(&["A", "B"], &["c1", "c2"]);
        let one = header(&["B", "A"], &["c2", "c1"]);
        let two = header(&["X", "A", "Y", "B"], &["c0", "c1", "c2"]);
        let want = (b"c2".to_vec(), vec![b"A".to_vec(), b"B".to_vec()]);

        let mut rec = record(&mut src);
        rec.translate(&one, &mut src).expect("translate");
        assert_eq!(names(&mut rec, &one), want);
        // HTSlib would use the first destination's tables here.
        let mut rec = record(&mut src);
        rec.translate(&two, &mut src).expect("translate");
        assert_eq!(names(&mut rec, &two), want);

        // A clone of the destination is another destination, and a clone of the source keeps
        // nothing it could misuse.
        let mut cloned = src.clone();
        let one_again = one.clone();
        let mut rec = record(&mut cloned);
        rec.translate(&one_again, &mut cloned).expect("translate");
        assert_eq!(names(&mut rec, &one_again), want);
        let mut rec = record(&mut cloned);
        rec.translate(&two, &mut cloned).expect("translate");
        assert_eq!(names(&mut rec, &two), want);
    }

    #[test]
    fn tables_are_rebuilt_when_the_destination_changes() {
        let mut src = header(&["A", "B"], &["c1", "c2"]);
        let mut dst = header(&["A"], &["c1"]);
        let mut rec = record(&mut src);
        rec.translate(&dst, &mut src).expect("translate");
        dst.append(b"##INFO=<ID=B,Number=1,Type=Integer,Description=\"B\">")
            .expect("append");
        dst.append(b"##contig=<ID=c2>").expect("append");
        let mut rec = record(&mut src);
        rec.translate(&dst, &mut src).expect("translate");
        assert_eq!(
            names(&mut rec, &dst),
            (b"c2".to_vec(), vec![b"A".to_vec(), b"B".to_vec()])
        );
    }
}
