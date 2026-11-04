//! Choosing the allele to apply from FORMAT/GT, and IUPAC codes (the genotype parts of
//! bcftools' `apply_variant`, `iupac_add_gt`, and `iupac_set_allele`).

/// BCF's `int8` type code (`BCF_BT_INT8`).
pub const BT_INT8: u8 = 1;

/// The `int8` end-of-vector marker (`bcf_int8_vector_end`).
pub const INT8_VECTOR_END: u8 = 0x81;

/// `bcf_gt_is_missing` on a raw `int8` genotype byte.
#[inline]
pub fn gt_is_missing(v: u8) -> bool {
    v >> 1 == 0
}

/// `bcf_gt_allele` on a raw `int8` genotype byte.
#[inline]
pub fn gt_allele(v: u8) -> i64 {
    i64::from(v >> 1) - 1
}

/// `bcf_gt_is_phased` on a raw `int8` genotype byte.
#[inline]
pub fn gt_is_phased(v: u8) -> bool {
    v & 1 != 0
}

/// bcftools' `iupac2bitmask`: the bases an IUPAC code stands for, or `None` for other
/// characters.
pub fn iupac_to_bitmask(c: u8) -> Option<u8> {
    const A: u8 = 1;
    const C: u8 = 2;
    const G: u8 = 4;
    const T: u8 = 8;
    let c = if c >= 97 { c - 32 } else { c };
    Some(match c {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T,
        b'M' => A | C,
        b'R' => A | G,
        b'W' => A | T,
        b'S' => C | G,
        b'Y' => C | T,
        b'K' => G | T,
        b'V' => A | C | G,
        b'H' => A | C | T,
        b'D' => A | G | T,
        b'B' => C | G | T,
        b'N' => A | C | G | T,
        _ => return None,
    })
}

/// bcftools' `bitmask2iupac`.
pub fn bitmask_to_iupac(mask: u8) -> u8 {
    const CODES: &[u8; 16] = b".ACMGRSVTWYHKDBN";
    if mask == 0 || mask > 15 {
        0
    } else {
        CODES[mask as usize]
    }
}

/// Whether every base of an allele is an IUPAC code (not symbolic, a breakend, or `*`).
fn is_iupac_allele(allele: &[u8]) -> bool {
    allele.iter().all(|&c| iupac_to_bitmask(c).is_some())
}

/// Merges the alleles marked in `used` (all with `None`) into an IUPAC-coded copy of the
/// longest ALT, written over that allele. Returns its index, or the first used allele (an ALT
/// in preference to REF) when no ALT is a plain sequence (bcftools' `iupac_set_allele`; with
/// `None`, the REF,ALT branch of `apply_variant`).
pub fn set_iupac_allele(
    als: &mut [Vec<u8>],
    used: Option<&[bool]>,
    bitmask: &mut Vec<u8>,
) -> Option<usize> {
    bitmask.clear();
    let mut alt_len = 0usize;
    let mut ialt = None;
    let mut fallback: Option<usize> = None;
    for (i, allele) in als.iter().enumerate() {
        if used.is_some_and(|u| !u.get(i).copied().unwrap_or(false)) {
            continue;
        }
        if fallback.is_none_or(|f| f == 0) {
            fallback = Some(i);
        }
        if !is_iupac_allele(allele) {
            continue;
        }
        if allele.len() > bitmask.len() {
            bitmask.resize(allele.len(), 0);
        }
        if i > 0 && allele.len() > alt_len {
            alt_len = allele.len();
            ialt = Some(i);
        }
        for (m, &c) in bitmask.iter_mut().zip(allele) {
            *m |= iupac_to_bitmask(c).unwrap_or(0);
        }
    }
    if let Some(i) = ialt {
        for (dst, &m) in als[i].iter_mut().zip(bitmask.iter()) {
            *dst = bitmask_to_iupac(m);
        }
        return Some(i);
    }
    fallback
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_round_trip() {
        for c in b"ACGTMRWSYKVHDBN" {
            assert_eq!(bitmask_to_iupac(iupac_to_bitmask(*c).unwrap()), *c);
        }
        assert_eq!(iupac_to_bitmask(b'a'), Some(1));
        assert_eq!(iupac_to_bitmask(b'<'), None);
    }

    #[test]
    fn longest_alt_gets_codes() {
        let mut als = vec![b"A".to_vec(), b"C".to_vec(), b"AT".to_vec()];
        let mut mask = Vec::new();
        assert_eq!(set_iupac_allele(&mut als, None, &mut mask), Some(2));
        assert_eq!(als[2], b"MT");
        let mut als = vec![b"A".to_vec(), b"<DEL>".to_vec()];
        let used = [true, true];
        assert_eq!(set_iupac_allele(&mut als, Some(&used), &mut mask), Some(1));
        let used = [true, false];
        assert_eq!(set_iupac_allele(&mut als, Some(&used), &mut mask), Some(0));
    }
}
