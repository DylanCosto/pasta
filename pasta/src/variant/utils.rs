//! Record utilities from HTSlib's `vcfutils.c`: allele counts, genotype types, and allele
//! removal.

mod remove_alleles;

pub use self::remove_alleles::{remove_allele_set, remove_alleles, trim_alleles};
use crate::variant::{
    Header, HeaderRecordType, Record, Unpack,
    encoding::{self, INT8_VECTOR_END, INT16_VECTOR_END, INT32_VECTOR_END},
    log,
    record::FormatField,
};

/// Where [`calc_ac`] may take allele counts from.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum AcSource {
    /// INFO/AC and INFO/AN when present, else the genotypes (`BCF_UN_INFO|BCF_UN_FMT`).
    InfoThenGenotypes,
    /// Only the genotypes (`BCF_UN_FMT`).
    Genotypes,
    /// Only INFO/AC and INFO/AN (`BCF_UN_INFO`).
    Info,
}

/// The type of a genotype (HTSlib's `GT_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum GenotypeType {
    /// Homozygous reference (`GT_HOM_RR`).
    HomRef,
    /// Heterozygous reference/alternate (`GT_HET_RA`).
    HetRefAlt,
    /// Homozygous alternate (`GT_HOM_AA`).
    HomAlt,
    /// Heterozygous between two alternate alleles (`GT_HET_AA`).
    HetAltAlt,
    /// Haploid reference (`GT_HAPL_R`).
    HaploidRef,
    /// Haploid alternate (`GT_HAPL_A`).
    HaploidAlt,
    /// Missing (`GT_UNKN`).
    Unknown,
}

/// An error from [`calc_ac`] that HTSlib treats as fatal.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CalcAcError(pub String);

impl std::fmt::Display for CalcAcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CalcAcError {}

/// Reads value `i` of an integer field as `i32`, keeping the vector-end sentinel distinct.
#[inline]
fn gt_value(data: &[u8], ty: u8, i: usize) -> Option<i32> {
    Some(match ty {
        1 => {
            let v = data[i] as i8;
            if v == INT8_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(v)
            }
        }
        2 => {
            let v = encoding::le_i16(&data[i * 2..]);
            if v == INT16_VECTOR_END {
                INT32_VECTOR_END
            } else {
                i32::from(v)
            }
        }
        3 => encoding::le_i32(&data[i * 4..]),
        _ => return None,
    })
}

/// The most alleles [`calc_ac`] counts in a pass each; past it, values are counted one by one.
const MAX_COUNTED_ALLELES: usize = 16;

/// Counts alleles (HTSlib's `bcf_calc_ac`): `ac[i]` for allele `i`, from INFO/AC and AN or
/// from the genotypes. Returns whether counts were found.
///
/// # Errors
///
/// Returns the errors HTSlib exits on: AN smaller than the sum of AC, an allele index out of
/// range, or an unexpected type.
pub fn calc_ac(
    header: &Header,
    line: &mut Record,
    ac: &mut Vec<i32>,
    which: AcSource,
) -> Result<bool, CalcAcError> {
    let n_allele = line.allele_count();
    ac.clear();
    ac.resize(n_allele.max(1), 0);
    let rid = line.rid();
    let chrom =
        || String::from_utf8_lossy(header.contig_name(rid).unwrap_or(b"(null)")).into_owned();
    if which != AcSource::Genotypes {
        line.unpack(Unpack::INFO);
        let an_id = header.id(b"AN");
        let ac_id = header.id(b"AC");
        let mut an: i64 = -1;
        let mut ac_field = None;
        if let (Some(an_id), Some(ac_id)) = (an_id, ac_id) {
            for f in line.info_fields() {
                if f.key() == an_id as i32 {
                    an = f.int_value();
                } else if f.key() == ac_id as i32 {
                    ac_field = Some(f.clone());
                }
            }
        }
        if an >= 0
            && let Some(f) = ac_field
        {
            if f.len() != n_allele.saturating_sub(1) {
                static WARNED: std::sync::atomic::AtomicBool =
                    std::sync::atomic::AtomicBool::new(false);
                log::warn_once(
                    &WARNED,
                    "bcf_calc_ac",
                    format_args!(
                        "Incorrect number of AC fields at {}:{}. (This message is printed only once.)\n",
                        chrom(),
                        line.pos() + 1
                    ),
                );
                return Ok(false);
            }
            let data = line.info_data(&f);
            let mut nac: i64 = 0;
            for i in 0..f.len() {
                let v = match f.type_code() {
                    1 => i64::from(data[i] as i8),
                    2 => i64::from(encoding::le_i16(&data[i * 2..])),
                    3 => i64::from(encoding::le_i32(&data[i * 4..])),
                    t => {
                        log::error(
                            "bcf_calc_ac",
                            format_args!("Unexpected type {t} at {}:{}", chrom(), line.pos() + 1),
                        );
                        return Err(CalcAcError(format!("unexpected type {t}")));
                    }
                };
                ac[i + 1] = v as i32;
                nac += v;
            }
            if an < nac {
                log::error(
                    "bcf_calc_ac",
                    format_args!("Incorrect AN/AC counts at {}:{}", chrom(), line.pos() + 1),
                );
                return Err(CalcAcError("incorrect AN/AC counts".into()));
            }
            ac[0] = (an - nac) as i32;
            return Ok(true);
        }
    }
    if which != AcSource::Info {
        let Some(gt_id) = header.id(b"GT") else {
            return Ok(false);
        };
        line.unpack(Unpack::FMT);
        let Some(gt) = line
            .format_fields()
            .iter()
            .find(|f| f.id() == gt_id as i32)
            .cloned()
        else {
            return Ok(false);
        };
        let n = gt.n();
        let incorrect_allele = |i: usize, allele: i32| {
            let sample = String::from_utf8_lossy(&header.samples()[i]).into_owned();
            log::error(
                "bcf_calc_ac",
                format_args!(
                    "Incorrect allele (\"{}\") in {sample} at {}:{}",
                    allele - 1,
                    chrom(),
                    line.pos() + 1
                ),
            );
            CalcAcError("incorrect allele".into())
        };
        if gt.type_code() == 1 && n > 0 {
            // The usual 8-bit genotypes, read straight from the field's block.
            let data = line.format_data(&gt);
            let data = &data[..data.len().min(n * line.sample_count())];
            // When every value is an allele or missing (no vector end, nothing out of range),
            // count each allele in a pass of its own, which the compiler vectorizes; adding to
            // `ac` value by value waits on the previous addition to the same count.
            let max_allele = i32::try_from(n_allele).unwrap_or(i32::MAX);
            let plain = data.iter().fold(true, |ok, &b| {
                let allele = i32::from(b as i8) >> 1;
                ok & (allele >= 0) & (allele <= max_allele)
            });
            if plain && n_allele <= MAX_COUNTED_ALLELES {
                for (k, count) in (1..=max_allele).zip(ac.iter_mut()) {
                    *count = data
                        .iter()
                        .map(|&b| i32::from(i32::from(b as i8) >> 1 == k))
                        .sum();
                }
                return Ok(true);
            }
            for (i, p) in data.chunks_exact(n).take(line.sample_count()).enumerate() {
                for &b in p {
                    let v = b as i8;
                    if v == INT8_VECTOR_END {
                        break;
                    }
                    let allele = i32::from(v) >> 1;
                    if allele == 0 {
                        continue;
                    }
                    if allele as usize > n_allele {
                        return Err(incorrect_allele(i, allele));
                    }
                    ac[(allele - 1) as usize] += 1;
                }
            }
            return Ok(true);
        }
        for i in 0..line.sample_count() {
            let p = line.format_sample_data(&gt, i);
            for ial in 0..n {
                let Some(val) = gt_value(p, gt.type_code(), ial) else {
                    log::error(
                        "bcf_calc_ac",
                        format_args!(
                            "Unexpected type {} at {}:{}",
                            gt.type_code(),
                            chrom(),
                            line.pos() + 1
                        ),
                    );
                    return Err(CalcAcError("unexpected GT type".into()));
                };
                if val == INT32_VECTOR_END {
                    break;
                }
                if val >> 1 == 0 {
                    continue;
                }
                if (val >> 1) as usize > n_allele {
                    return Err(incorrect_allele(i, val >> 1));
                }
                ac[((val >> 1) - 1) as usize] += 1;
            }
        }
        return Ok(true);
    }
    Ok(false)
}

/// Classifies one sample's genotype (HTSlib's `bcf_gt_type`), also returning the smallest and
/// largest alternate allele seen (0-based; 0 when none).
#[must_use]
pub fn gt_type(rec: &Record, gt: &FormatField, sample: usize) -> (GenotypeType, i32, i32) {
    let p = rec.format_sample_data(gt, sample);
    let mut nals = 0;
    let mut has_ref = false;
    let mut has_alt = 0;
    let mut ial = 0i32;
    let mut jal = 0i32;
    for i in 0..gt.n() {
        let Some(val) = gt_value(p, gt.type_code(), i) else {
            break;
        };
        if val == INT32_VECTOR_END {
            break;
        }
        if val >> 1 == 0 {
            return (GenotypeType::Unknown, 0, 0);
        }
        let tmp = val >> 1;
        if tmp > 1 {
            if ial == 0 {
                ial = tmp;
                has_alt = 1;
            } else if tmp != ial {
                if tmp < ial {
                    jal = ial;
                    ial = tmp;
                } else {
                    jal = tmp;
                }
                has_alt = 2;
            }
        } else {
            has_ref = true;
        }
        nals += 1;
    }
    let ial_out = if ial > 0 { ial - 1 } else { ial };
    let jal_out = if jal > 0 { jal - 1 } else { jal };
    let ty = if nals == 0 {
        GenotypeType::Unknown
    } else if nals == 1 {
        if has_ref {
            GenotypeType::HaploidRef
        } else {
            GenotypeType::HaploidAlt
        }
    } else if !has_ref {
        if has_alt == 1 {
            GenotypeType::HomAlt
        } else {
            GenotypeType::HetAltAlt
        }
    } else if has_alt == 0 {
        GenotypeType::HomRef
    } else {
        GenotypeType::HetRefAlt
    };
    (ty, ial_out, jal_out)
}

/// Returns the GT field of a record, unpacking FORMAT (HTSlib's `bcf_get_fmt(hdr, rec,
/// "GT")`).
pub fn gt_field(header: &Header, rec: &mut Record) -> Option<FormatField> {
    let id = header
        .id(b"GT")
        .filter(|&id| header.id_exists(HeaderRecordType::Format, id))?;
    rec.format_field_by_id(id as i32).cloned()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::vcf;

    /// `calc_ac` from the genotypes of each record of `text`, or `None` for an error.
    fn genotype_counts(text: &str) -> Vec<Option<Vec<i32>>> {
        let split = text.find("\n#CHROM").expect("a #CHROM line");
        let body_start = text[split + 1..].find('\n').expect("records") + split + 2;
        let mut header = Header::parse(&text.as_bytes()[..body_start]).expect("a header");
        let (counts, _) = log::capture(|| {
            text[body_start..]
                .lines()
                .map(|line| {
                    let mut rec = Record::new();
                    let mut buf = line.as_bytes().to_vec();
                    vcf::parse(&mut header, &mut buf, &mut rec).expect("a record");
                    let mut ac = Vec::new();
                    calc_ac(&header, &mut rec, &mut ac, AcSource::Genotypes)
                        .ok()
                        .map(|_| ac)
                })
                .collect()
        });
        counts
    }

    #[test]
    fn genotype_counts_take_missing_haploid_and_bad_alleles_into_account() {
        let text = "##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0|0\t0/1\t1|1\t./.
1\t2\t.\tA\tC,G\t.\t.\t.\tGT\t0|2\t2/2\t1\t.|1
1\t3\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0\t1|.\t1/1
1\t4\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/2\t1/1\t0/1
";
        assert_eq!(
            genotype_counts(text),
            [
                Some(vec![3, 3]),
                Some(vec![1, 2, 3]),
                Some(vec![3, 3]),
                // Allele 2 of a site with two alleles.
                None,
            ]
        );
    }
}
