//! The reference length of a record (HTSlib's `get_rlen`): the longest of REF, INFO/END,
//! INFO/SVLEN of `<DEL>`-like alleles, and FORMAT/LEN of gVCF `<*>` blocks.

use super::{
    FieldData, FormatField, InfoField, Record, Unpack, unpack_fmt_core1, unpack_info_core1,
};
use crate::variant::{
    Header, HeaderRecordType,
    encoding::{self, type_shift},
    header::SpecialTag,
};

/// Returns whether SVLEN gives the reference length for an ALT allele: a symbolic `<CNV>`,
/// `<DEL>`, `<DUP>`, or `<INV>` (or a subtype such as `<DEL:ME>`) (HTSlib's
/// `svlen_on_ref_for_vcf_alt`).
#[must_use]
pub fn svlen_on_ref_for_alt(alt: &[u8]) -> bool {
    if alt.first() != Some(&b'<') || alt.len() < 5 {
        return false;
    }
    if alt[4] != b'>' && alt[4] != b':' {
        return false;
    }
    if !matches!(&alt[..4], b"<CNV" | b"<DEL" | b"<DUP" | b"<INV") {
        return false;
    }
    alt[alt.len() - 1] == b'>'
}

/// Finds an INFO field by ID among decoded fields, as `bcf_get_info` does (the first match,
/// even one marked for removal), when the tag is defined as INFO.
fn decoded_info<'a>(header: &Header, rec: &'a Record, name: &[u8]) -> Option<&'a InfoField> {
    let id = header.id(name)?;
    if !header.id_exists(HeaderRecordType::Info, id) {
        return None;
    }
    rec.d.info[..rec.n_info as usize]
        .iter()
        .find(|f| f.key == id as i32)
}

fn decoded_fmt<'a>(header: &Header, rec: &'a Record, name: &[u8]) -> Option<&'a FormatField> {
    let id = header.id(name)?;
    if !header.id_exists(HeaderRecordType::Format, id) {
        return None;
    }
    rec.d.fmt[..(rec.n_fmt as usize).min(rec.d.fmt.len())]
        .iter()
        .find(|f| f.id == id as i32)
}

/// Reads value `i` of an integer field of raw type `ty`, mapping missing to 0; `None` for a
/// non-integer type.
fn int_at(data: &[u8], ty: u8, i: usize) -> Option<i64> {
    Some(match ty {
        1 => {
            let v = data[i] as i8;
            if v == encoding::INT8_MISSING {
                0
            } else {
                i64::from(v)
            }
        }
        2 => {
            let v = encoding::le_i16(&data[i * 2..]);
            if v == encoding::INT16_MISSING {
                0
            } else {
                i64::from(v)
            }
        }
        3 => {
            let v = encoding::le_i32(&data[i * 4..]);
            if v == encoding::INT32_MISSING {
                0
            } else {
                i64::from(v)
            }
        }
        4 => {
            let v = encoding::le_i64(&data[i * 8..]);
            if v == encoding::INT64_MISSING { 0 } else { v }
        }
        _ => return None,
    })
}

/// Computes the reference length (HTSlib's `get_rlen`).
pub(crate) fn get_rlen(header: &Header, v: &Record) -> i64 {
    let shared = &v.shared;
    let e = shared.len();
    let mut f = 0usize;
    let n_allele = v.n_allele as usize;
    let mut svlenals = vec![0u8; 1 + n_allele / 8];
    let mut use_svlen = false;
    let mut gvcf = false;
    let mut len_ref: i64 = 0;
    let endid = header.special_id(SpecialTag::End).map_or(-1, i64::from);

    if v.unpacked & Unpack::STR.bits() != 0 || v.d.shared_dirty & super::DIRTY_ALS != 0 {
        for i in 1..n_allele {
            let a = v.allele(i);
            if a.first() != Some(&b'<') {
                continue;
            }
            if svlen_on_ref_for_alt(a) {
                svlenals[i >> 3] |= 1 << (i & 7);
                use_svlen = true;
            } else if a == b"<*>" || a == b"<NON_REF>" {
                gvcf = true;
            }
        }
        f += v.unpack_size[0] + v.unpack_size[1];
        len_ref = if n_allele > 0 {
            v.allele(0).len() as i64
        } else {
            0
        };
    } else if f < e {
        let (size, ty, used) = encoding::dec_size(&shared[f..]);
        f += used + ((size.max(0) as usize) << type_shift(ty));
        for i in 0..n_allele {
            let (size, ty, used) = encoding::dec_size(&shared[f..]);
            f += used;
            let size = size.max(0) as usize;
            if i == 0 {
                len_ref = size as i64;
            } else if size > 0 && shared[f] == b'<' {
                let a = &shared[f..f + size];
                if svlen_on_ref_for_alt(a) {
                    svlenals[i >> 3] |= 1 << (i & 7);
                    use_svlen = true;
                } else if a == b"<*>" || a == b"<NON_REF>" {
                    gvcf = true;
                }
            }
            f += size << type_shift(ty);
        }
    }
    // FILTER
    if v.unpacked & Unpack::FLT.bits() != 0 {
        f += v.unpack_size[2];
    } else if f < e {
        let (size, ty, used) = encoding::dec_size(&shared[f..]);
        f += used + ((size.max(0) as usize) << type_shift(ty));
    }

    let svlenid = if use_svlen {
        header.special_id(SpecialTag::Svlen).map_or(-1, i64::from)
    } else {
        -1
    };

    // INFO: owned copies of the END and SVLEN fields (with their bytes) when read packed.
    let mut end_lcl: Option<(InfoField, &[u8])> = None;
    let mut svlen_lcl: Option<(InfoField, &[u8])> = None;
    if svlenid >= 0 || endid >= 0 {
        if v.unpacked & Unpack::INFO.bits() != 0 || v.d.shared_dirty & super::DIRTY_INF != 0 {
            end_lcl = decoded_info(header, v, b"END").map(|i| (i.clone(), v.info_data(i)));
            svlen_lcl = decoded_info(header, v, b"SVLEN").map(|i| (i.clone(), v.info_data(i)));
        } else if f < e {
            for _ in 0..v.n_info {
                let (id, used) = encoding::dec_typed_int1(&shared[f..]);
                let t;
                if id == endid {
                    let mut field = InfoField::default();
                    t = unpack_info_core1(shared, f, &mut field);
                    let data = &shared[f + field.vptr_off as usize..t];
                    end_lcl = Some((field, data));
                    if svlen_lcl.is_some() || svlenid < 0 {
                        break;
                    }
                } else if id == svlenid {
                    let mut field = InfoField::default();
                    t = unpack_info_core1(shared, f, &mut field);
                    let data = &shared[f + field.vptr_off as usize..t];
                    svlen_lcl = Some((field, data));
                    if end_lcl.is_some() || endid < 0 {
                        break;
                    }
                } else {
                    let g = f + used;
                    let (size, ty, used) = encoding::dec_size(&shared[g..]);
                    t = g + used + ((size.max(0) as usize) << type_shift(ty));
                }
                f = t;
            }
        }
    }

    let lenid = if gvcf {
        header.special_id(SpecialTag::Len).map_or(-1, i64::from)
    } else {
        -1
    };
    let mut len_lcl: Option<(FormatField, &[u8])> = None;
    if lenid >= 0 {
        let indiv = &v.indiv;
        if v.unpacked & Unpack::FMT.bits() != 0 || v.d.indiv_dirty {
            len_lcl = decoded_fmt(header, v, b"LEN").map(|f| (f.clone(), v.format_data(f)));
        } else if !indiv.is_empty() {
            let mut f = 0usize;
            for _ in 0..v.n_fmt {
                let (id, used) = encoding::dec_typed_int1(&indiv[f..]);
                if id == lenid {
                    let mut field = FormatField::default();
                    let t = unpack_fmt_core1(indiv, f, v.n_sample as usize, &mut field);
                    let data = &indiv[f + field.p_off as usize..t];
                    len_lcl = Some((field, data));
                    break;
                }
                let g = f + used;
                let (size, ty, used) = encoding::dec_size(&indiv[g..]);
                f = g + used + ((size.max(0) as usize * v.n_sample as usize) << type_shift(ty));
            }
        }
    }

    let pos = v.pos;
    let mut end: i64 = 0;
    if let Some((info, _)) = &end_lcl
        && info.data != FieldData::Removed
    {
        end = info.v1;
        end = match info.ty {
            1 if end == i64::from(encoding::INT8_MISSING) => 0,
            2 if end == i64::from(encoding::INT16_MISSING) => 0,
            3 if end == i64::from(encoding::INT32_MISSING) => 0,
            4 if end == encoding::INT64_MISSING => 0,
            1..=4 => end,
            _ => 0,
        };
    }

    let mut len: i64 = 0;
    if let Some((info, data)) = &svlen_lcl
        && info.data != FieldData::Removed
    {
        let mut i = 0usize;
        while (i as i32) < info.len && i + 1 < n_allele {
            // SVLEN value i belongs to allele i + 1. HTSlib reads the byte of allele i
            // (`svlenals[i >> 3]`), the wrong byte at ALT 8, 16, ...; vtools reads allele
            // i + 1's own bit.
            if svlenals[(i + 1) >> 3] & (1 << ((i + 1) & 7)) == 0 {
                i += 1;
                continue;
            }
            let Some(tmp) = int_at(data, info.ty, i) else {
                len = 0;
                break;
            };
            let tmp = tmp.abs();
            if len < tmp {
                len = tmp;
            }
            i += 1;
        }
    }
    if (svlen_lcl.is_none() || len == 0) && end != 0 {
        len = if end > pos { end - pos - 1 } else { 0 };
    }
    let end_svlen = pos + len + 1;

    len = 0;
    if let Some((fmt, data)) = &len_lcl
        && fmt.data != FieldData::Removed
    {
        let mut offset = 0usize;
        let n = fmt.n.max(0) as usize;
        'samples: for _ in 0..v.n_sample {
            let mut j = 0;
            while j < n {
                let Some(tmp) = int_at(&data[offset..], fmt.ty, j) else {
                    len = 0;
                    break 'samples;
                };
                if len < tmp {
                    len = tmp;
                }
                j += 1;
            }
            offset += j << type_shift(fmt.ty);
        }
    }
    if (len_lcl.is_none() || len == 0) && end != 0 {
        len = if end > pos { end - pos } else { 0 };
    }
    let end_fmtlen = pos + len;

    let hpos = end.max(end_svlen).max(end_fmtlen);
    let len = hpos - pos;
    len.max(len_ref)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parses `line` and returns its rlen as parsed from VCF and as recomputed from packed BCF.
    fn rlens(header: &mut Header, line: &str) -> (i64, i64) {
        let mut text = line.as_bytes().to_vec();
        let mut parsed = Record::new();
        crate::variant::vcf::parse(header, &mut text, &mut parsed).expect("the line parses");
        let mut packed = Vec::new();
        crate::variant::bcf::pack_record(&mut parsed, &mut packed);
        let mut rec = Record::new();
        crate::variant::bcf::unpack_record(&packed, &mut rec).expect("the record unpacks");
        (parsed.rlen, get_rlen(header, &rec))
    }

    /// Builds a record at POS 100 with `n` ALT alleles: `svs` places symbolic alleles (1-based
    /// ALT index, allele, SVLEN value) and the rest are distinct short alleles with SVLEN `.`.
    fn sv_line(n: usize, svs: &[(usize, &str, &str)]) -> String {
        const FILL: [&str; 17] = [
            "C", "G", "T", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT",
            "TA", "TC",
        ];
        let (mut alt, mut svlen) = (Vec::new(), Vec::new());
        for k in 1..=n {
            let (a, l) = svs
                .iter()
                .find(|sv| sv.0 == k)
                .map_or((FILL[k - 1], "."), |sv| (sv.1, sv.2));
            alt.push(a);
            svlen.push(l);
        }
        format!(
            "1\t100\t.\tA\t{}\t.\t.\tSVLEN={}",
            alt.join(","),
            svlen.join(",")
        )
    }

    #[test]
    fn svlen_counts_only_svlen_on_ref_alleles_at_every_alt_index() {
        let mut header = Header::parse(
            b"##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n\
              ##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"l\">\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .expect("the header parses");
        let mut lines = Vec::new();
        // A 50-base `<DEL>` as the last of n ALTs, across the byte boundaries of the bitset.
        for n in [1, 7, 8, 9, 15, 16, 17] {
            lines.push(sv_line(n, &[(n, "<DEL>", "-50")]));
        }
        // The SVLEN of an `<INS>` at ALT 16 is not on REF, even with a `<DEL>` at ALT 8.
        lines.push(sv_line(16, &[(8, "<DEL>", "-50"), (16, "<INS>", "1000")]));
        for line in lines {
            assert_eq!(rlens(&mut header, &line), (51, 51), "{line}");
        }
    }
}
