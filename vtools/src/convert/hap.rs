//! The hidden fields behind `bcftools convert`'s IMPUTE2 and SHAPEIT outputs: genotype
//! probabilities (`%_GT_TO_PROB3`, `%_PL_TO_PROB3`, `%_GP_TO_PROB3`) and haplotypes
//! (`%_GT_TO_HAP`, `%_GT_TO_HAP2`). They print all samples at once, in header order.

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Record,
    encoding::{
        INT8_MISSING, INT8_VECTOR_END, INT32_VECTOR_END, float_is_missing, float_is_vector_end,
    },
    number::push_int,
};

use super::{BT_INT8, State, location, push_fixed6};

/// HTSlib's `bcf_gt_is_missing`.
fn gt_is_missing(v: i32) -> bool {
    v >> 1 == 0
}

/// HTSlib's `bcf_gt_allele`.
fn gt_allele(v: i32) -> i32 {
    (v >> 1) - 1
}

/// `%_GT_TO_PROB3`: three genotype probabilities per sample from GT.
pub(super) fn gt_to_prob3(
    st: &mut State,
    header: &Header,
    rec: &mut Record,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let n = rec.get_genotypes(header, &mut st.ibuf).unwrap_or(0);
    let nsmpl = st.samples.len();
    if n == 0 || nsmpl == 0 {
        bail!("Error parsing GT tag at {}", location(header, rec));
    }
    let n = n / nsmpl;
    for i in 0..nsmpl {
        let ptr = &st.ibuf[i * n..(i + 1) * n];
        let ploidy = ptr.iter().take_while(|&&v| v != INT32_VECTOR_END).count();
        match ploidy {
            2 => {
                if gt_is_missing(ptr[0]) {
                    out.extend_from_slice(b" 0.33 0.33 0.33");
                } else if gt_allele(ptr[0]) != gt_allele(ptr[1]) {
                    out.extend_from_slice(b" 0 1 0");
                } else if gt_allele(ptr[0]) == 1 {
                    out.extend_from_slice(b" 0 0 1");
                } else {
                    out.extend_from_slice(b" 1 0 0");
                }
            }
            1 => {
                if gt_is_missing(ptr[0]) {
                    out.extend_from_slice(b" 0.5 0.0 0.5");
                } else if gt_allele(ptr[0]) == 1 {
                    out.extend_from_slice(b" 0 0 1");
                } else {
                    out.extend_from_slice(b" 1 0 0");
                }
            }
            j => bail!("FIXME: not ready for ploidy {j}"),
        }
    }
    Ok(())
}

/// `%_PL_TO_PROB3`: genotype probabilities from phred-scaled likelihoods.
pub(super) fn pl_to_prob3(
    st: &mut State,
    header: &Header,
    rec: &mut Record,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let n = rec.get_format_i32(header, b"PL", &mut st.ibuf).unwrap_or(0);
    let nsmpl = st.samples.len();
    if n == 0 || nsmpl == 0 {
        bail!("Error parsing PL tag at {}", location(header, rec));
    }
    let n = n / nsmpl;
    let prob = |pl: i32| 10f64.powf(-0.1 * f64::from(pl));
    for i in 0..nsmpl {
        let ptr = &st.ibuf[i * n..(i + 1) * n];
        let mut sum = 0f32;
        let mut j = 0;
        while j < n && ptr[j] != INT32_VECTOR_END {
            sum = (f64::from(sum) + prob(ptr[j])) as f32;
            j += 1;
        }
        let sum = f64::from(sum);
        let at = |k: usize| ptr.get(k).copied().unwrap_or(INT32_VECTOR_END);
        out.push(b' ');
        push_fixed6(out, prob(at(0)) / sum);
        if j == rec.allele_count() {
            // Haploid.
            out.extend_from_slice(b" 0 ");
            push_fixed6(out, prob(at(1)) / sum);
        } else {
            out.push(b' ');
            push_fixed6(out, prob(at(1)) / sum);
            out.push(b' ');
            push_fixed6(out, prob(at(2)) / sum);
        }
    }
    Ok(())
}

/// `%_GP_TO_PROB3`: genotype probabilities from FORMAT/GP.
pub(super) fn gp_to_prob3(
    st: &mut State,
    header: &Header,
    rec: &mut Record,
    out: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let n = rec.get_format_f32(header, b"GP", &mut st.fbuf).unwrap_or(0);
    let nsmpl = st.samples.len();
    if n == 0 || nsmpl == 0 {
        bail!("Error parsing GP tag at {}", location(header, rec));
    }
    let n = n / nsmpl;
    for i in 0..nsmpl {
        let ptr = &mut st.fbuf[i * n..(i + 1) * n];
        let mut j = 0;
        while j < n {
            if float_is_vector_end(ptr[j]) {
                break;
            }
            if float_is_missing(ptr[j]) {
                ptr[j] = 0.0;
                j += 1;
                continue;
            }
            if ptr[j] < 0.0 || ptr[j] > 1.0 {
                let mut v = Vec::new();
                push_fixed6(&mut v, f64::from(ptr[j]));
                bail!(
                    "[{}:{}] GP value outside range [0,1]; vtools convert expects the VCF4.3+ spec for the GP field encoding genotype posterior probabilities",
                    location(header, rec),
                    String::from_utf8_lossy(&v)
                );
            }
            j += 1;
        }
        let at = |k: usize| f64::from(ptr.get(k).copied().unwrap_or(0.0));
        let (a, b, c) = if j == rec.allele_count() {
            (at(0), 0.0, at(1))
        } else {
            (at(0), at(1), at(2))
        };
        for v in [a, b, c] {
            out.push(b' ');
            push_fixed6(out, v);
        }
    }
    Ok(())
}

/// How `%_GT_TO_HAP2` differs from `%_GT_TO_HAP`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Haploid {
    /// Haploid genotypes print as `A -` (`%_GT_TO_HAP`).
    Keep,
    /// Haploid genotypes print as diploid homozygotes (`%_GT_TO_HAP2`).
    ToDiploid,
}

/// `%_GT_TO_HAP` and `%_GT_TO_HAP2`: two haplotype columns per sample (`0`, `1`, `?`, `-`),
/// with `*` marking unphased alleles (bcftools' `process_gt_to_hap` and `process_gt_to_hap2`).
pub(super) fn gt_to_hap(
    st: &State,
    header: &Header,
    rec: &Record,
    out: &mut Vec<u8>,
    haploid: Haploid,
) -> anyhow::Result<()> {
    let gt_id = header.id(b"GT").map_or(-1, |i| i as i32);
    if !super::id_exists(header, HeaderRecordType::Format, gt_id) {
        bail!("FORMAT/GT tag not present at {}", location(header, rec));
    }
    let fmts = rec.format_fields();
    let Some(gt) = fmts.iter().find(|f| f.id() == gt_id) else {
        bail!("FORMAT/GT tag not present at {}", location(header, rec));
    };
    if rec.allele_count() > 100 {
        bail!(
            "Too many alleles ({}) at {}",
            rec.allele_count(),
            location(header, rec)
        );
    }
    if gt.type_code() != BT_INT8 {
        bail!(
            "Uh, too many alleles ({}) or redundant BCF representation at {}",
            rec.allele_count(),
            location(header, rec)
        );
    }
    let n = gt.n();
    if haploid == Haploid::Keep && n != 1 && n != 2 {
        bail!(
            "Uh, ploidy of {n} not supported, see {}",
            location(header, rec)
        );
    }
    let data = rec.format_data(gt);
    let push_allele = |out: &mut Vec<u8>, v: i8| push_int(out, i64::from(gt_allele(i32::from(v))));
    for i in 0..st.samples.len() {
        let p0 = data.get(i * n).map_or(INT8_VECTOR_END, |&b| b as i8);
        // With one value per sample, the second allele is absent (bcftools reads the next
        // sample's first allele here).
        let p1 = if n >= 2 {
            data.get(i * n + 1).map_or(INT8_VECTOR_END, |&b| b as i8)
        } else {
            INT8_VECTOR_END
        };
        if haploid == Haploid::Keep && n == 1 {
            match p0 >> 1 {
                1 => out.extend_from_slice(b"0 - "),
                2 => out.extend_from_slice(b"1 - "),
                _ if p0 == INT8_MISSING => out.extend_from_slice(b"? ? "),
                _ => {
                    push_allele(out, p0);
                    out.extend_from_slice(b" - ");
                }
            }
            continue;
        }
        let first = p0 >> 1;
        if first == 1 || first == 2 {
            let a = if first == 1 { b'0' } else { b'1' };
            match p1 {
                3 => out.extend_from_slice(&[a, b' ', b'0', b' ']),
                5 => out.extend_from_slice(&[a, b' ', b'1', b' ']),
                INT8_VECTOR_END => match haploid {
                    Haploid::Keep => out.extend_from_slice(&[a, b' ', b'-', b' ']),
                    Haploid::ToDiploid => out.extend_from_slice(&[a, b' ', a, b' ']),
                },
                2 => out.extend_from_slice(&[a, b'*', b' ', b'0', b'*', b' ']),
                4 => out.extend_from_slice(&[a, b'*', b' ', b'1', b'*', b' ']),
                _ if gt_is_missing(i32::from(p1)) => out.extend_from_slice(b"? ? "),
                _ if p1 & 1 != 0 => {
                    out.extend_from_slice(&[a, b' ']);
                    push_allele(out, p1);
                    out.push(b' ');
                }
                _ => {
                    out.extend_from_slice(&[a, b'*', b' ']);
                    push_allele(out, p1);
                    out.extend_from_slice(b"* ");
                }
            }
        } else if gt_is_missing(i32::from(p0)) {
            if haploid == Haploid::Keep && p1 == INT8_VECTOR_END {
                out.extend_from_slice(b"? - ");
            } else {
                out.extend_from_slice(b"? ? ");
            }
        } else if p1 == INT8_VECTOR_END {
            // Something else than the first ALT: use REF.
            match haploid {
                Haploid::Keep => out.extend_from_slice(b"0 - "),
                Haploid::ToDiploid => out.extend_from_slice(b"0 0 "),
            }
        } else {
            let unphased = p1 & 1 == 0;
            push_allele(out, p0);
            if unphased {
                out.push(b'*');
            }
            out.push(b' ');
            push_allele(out, p1);
            if unphased {
                out.push(b'*');
            }
            out.push(b' ');
        }
    }
    // Delete the last space.
    out.pop();
    Ok(())
}
