//! Matching alleles of records whose REF alleles may differ in length (bcftools' `vcmp.c`).
//!
//! Two records at one position, say `CA>C` and `C>T`, describe the same ALT only when the
//! longer REF's extra bases are appended to the shorter record's ALT. [`Vcmp::set_ref`] finds
//! those extra bases; [`Vcmp::find_allele`] then compares ALTs taking them into account.

/// A list of alleles, REF first: a record's, or an annotation line's.
pub(super) trait AlleleList {
    /// The number of alleles.
    fn count(&self) -> usize;
    /// Allele `i`.
    fn allele_at(&self, i: usize) -> &[u8];
}

impl AlleleList for pasta::variant::Record {
    fn count(&self) -> usize {
        self.allele_count()
    }

    fn allele_at(&self, i: usize) -> &[u8] {
        self.allele(i)
    }
}

impl<A: AsRef<[u8]>> AlleleList for [A] {
    fn count(&self) -> usize {
        self.len()
    }

    fn allele_at(&self, i: usize) -> &[u8] {
        self[i].as_ref()
    }
}

/// Allele matching state (bcftools' `vcmp_t`).
#[derive(Debug, Default)]
pub(super) struct Vcmp {
    /// The extra REF bases of the longer REF, upper-cased.
    dref: Vec<u8>,
    /// Which REF is longer: positive when the first, negative when the second, 0 when equal.
    ndref: isize,
    /// The allele map of the last [`Self::map_ar_values`].
    map: Vec<i32>,
    /// The diploid genotype map of the last [`Self::map_dip_g_values`].
    map_dip: Vec<i32>,
}

/// The length of the case-insensitive common prefix of two alleles.
fn common_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b)
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count()
}

/// Whether `tail` is the extra REF bases `dref` (already upper-cased).
fn is_dref(tail: &[u8], dref: &[u8]) -> bool {
    tail.len() == dref.len()
        && tail
            .iter()
            .zip(dref)
            .all(|(x, d)| x.to_ascii_uppercase() == *d)
}

impl Vcmp {
    /// Records how the REF alleles differ (bcftools' `vcmp_set_ref`); `false` when they are
    /// not compatible (neither is a prefix of the other).
    pub(super) fn set_ref(&mut self, ref1: &[u8], ref2: &[u8]) -> bool {
        self.ndref = 0;
        let n = common_prefix(ref1, ref2);
        if n == ref1.len() && n == ref2.len() {
            return true;
        }
        if n < ref1.len() && n < ref2.len() {
            return false;
        }
        self.dref.clear();
        if n < ref1.len() {
            self.dref
                .extend(ref1[n..].iter().map(u8::to_ascii_uppercase));
            self.ndref = self.dref.len() as isize;
        } else {
            self.dref
                .extend(ref2[n..].iter().map(u8::to_ascii_uppercase));
            self.ndref = -(self.dref.len() as isize);
        }
        true
    }

    /// Returns the index (counted from `from`) of the allele of `als1[from..]` matching `al2`
    /// (bcftools' `vcmp_find_allele`), after [`Self::set_ref`] with the REF of `als1`'s record
    /// first.
    pub(super) fn find_allele<L: AlleleList + ?Sized>(
        &self,
        als1: &L,
        from: usize,
        al2: &[u8],
    ) -> Option<usize> {
        (from..als1.count()).position(|i| self.allele_matches(als1.allele_at(i), al2))
    }

    fn allele_matches(&self, a: &[u8], b: &[u8]) -> bool {
        let n = common_prefix(a, b);
        let (ra, rb) = (&a[n..], &b[n..]);
        if !ra.is_empty() && !rb.is_empty() {
            return false;
        }
        if self.ndref == 0 {
            return ra.is_empty() && rb.is_empty();
        }
        if !ra.is_empty() {
            // The first allele is longer: its tail must be the extra REF bases.
            self.ndref > 0 && is_dref(ra, &self.dref)
        } else {
            self.ndref < 0 && is_dref(rb, &self.dref)
        }
    }

    /// Maps the `n` values of a Number=A (`n` = ALTs of `als2`) or Number=R (`n` = alleles of
    /// `als2`) field of the record with alleles `als2` to the alleles of `als1` (bcftools'
    /// `vcmp_map_ARvalues`): entry `i` is the index in `als1` (ALT-based for A) of allele `i`,
    /// or -1. `None` when the REF alleles are not compatible.
    pub(super) fn map_ar_values<A: AlleleList + ?Sized, B: AlleleList + ?Sized>(
        &mut self,
        n: usize,
        als1: &A,
        als2: &B,
    ) -> Option<&[i32]> {
        if !self.set_ref(als1.allele_at(0), als2.allele_at(0)) {
            return None;
        }
        let ifrom = usize::from(n != als2.count());
        self.map.clear();
        self.map.resize(n, -1);
        for i in ifrom..als2.count() {
            let found = self.find_allele(als1, ifrom, als2.allele_at(i));
            if let Some(slot) = self.map.get_mut(i - ifrom) {
                *slot = found.map_or(-1, |k| k as i32);
            }
        }
        Some(&self.map)
    }

    /// The allele map of the last [`Self::map_ar_values`].
    pub(super) fn map(&self) -> &[i32] {
        &self.map
    }

    /// The diploid map of the last [`Self::map_dip_g_values`].
    pub(super) fn map_dip(&self) -> &[i32] {
        &self.map_dip
    }

    /// Extends the last allele map to diploid genotypes (bcftools' `vcmp_map_dipGvalues`).
    pub(super) fn map_dip_g_values(&mut self) -> &[i32] {
        let n = self.map.len();
        self.map_dip.clear();
        for i in 0..n {
            for j in 0..=i {
                let (a, b) = (self.map[i], self.map[j]);
                self.map_dip.push(if a >= 0 && b >= 0 {
                    alleles_to_gt(a, b)
                } else {
                    -1
                });
            }
        }
        &self.map_dip
    }
}

/// The index of a diploid genotype (HTSlib's `bcf_alleles2gt`).
fn alleles_to_gt(a: i32, b: i32) -> i32 {
    if a > b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alleles_match_across_ref_lengths() {
        let mut v = Vcmp::default();
        assert!(v.set_ref(b"CA", b"C"));
        // C>T in the second record is CA>TA in the first.
        let als = [b"CA".as_slice(), b"TA", b"C"];
        assert_eq!(v.find_allele(&als[..], 1, b"T"), Some(0));
        assert_eq!(v.find_allele(&als[..], 1, b"G"), None);
        assert!(!v.set_ref(b"CA", b"G"));
        assert!(v.set_ref(b"c", b"C"));
        assert_eq!(
            v.find_allele(&[b"c".as_slice(), b"t"][..], 1, b"T"),
            Some(0)
        );
    }

    #[test]
    fn ar_maps_follow_vcmp() {
        let mut v = Vcmp::default();
        let als1 = [b"A".as_slice(), b"G", b"T"];
        let als2 = [b"A".as_slice(), b"T", b"C"];
        let map = v
            .map_ar_values(2, &als1[..], &als2[..])
            .expect("compatible");
        assert_eq!(map, [1, -1]);
        let map = v
            .map_ar_values(3, &als1[..], &als2[..])
            .expect("compatible");
        assert_eq!(map, [0, 2, -1]);
        assert_eq!(v.map_dip_g_values(), [0, 3, 5, -1, -1, -1]);
    }
}
