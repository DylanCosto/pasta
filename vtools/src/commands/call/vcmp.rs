//! Comparing alleles across records with different REF lengths (bcftools' `vcmp.c`).

/// bcftools' `vcmp_t`: the REF difference of two records being compared.
#[derive(Clone, Debug, Default)]
pub struct Vcmp {
    /// The uppercased extra REF suffix of the longer REF.
    dref: Vec<u8>,
    /// Positive when the first REF is longer, negative when the second is.
    ndref: i32,
}

impl Vcmp {
    /// Compares two REF alleles (`vcmp_set_ref`); `false` when they are incompatible.
    pub fn set_ref(&mut self, ref1: &[u8], ref2: &[u8]) -> bool {
        self.ndref = 0;
        let n = ref1
            .iter()
            .zip(ref2)
            .take_while(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        let (a_left, b_left) = (n < ref1.len(), n < ref2.len());
        if !a_left && !b_left {
            return true;
        }
        if a_left && b_left {
            return false;
        }
        let longer = if a_left { ref1 } else { ref2 };
        self.dref.clear();
        self.dref
            .extend(longer[n..].iter().map(u8::to_ascii_uppercase));
        self.ndref = self.dref.len() as i32;
        if b_left {
            self.ndref = -self.ndref;
        }
        true
    }

    /// Finds `al2` among `als1` (`vcmp_find_allele`), allowing for the REF difference.
    pub fn find_allele<A: AsRef<[u8]>>(&self, als1: &[A], al2: &[u8]) -> Option<usize> {
        als1.iter().position(|a1| {
            let a1 = a1.as_ref();
            let n = a1
                .iter()
                .zip(al2)
                .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
                .count();
            let (a, b) = (&a1[n..], &al2[n..]);
            if !a.is_empty() && !b.is_empty() {
                return false;
            }
            if self.ndref == 0 {
                return a.is_empty() && b.is_empty();
            }
            let tail = if !a.is_empty() {
                if self.ndref < 0 {
                    return false;
                }
                a
            } else {
                if self.ndref > 0 {
                    return false;
                }
                b
            };
            tail.len() == self.dref.len()
                && tail
                    .iter()
                    .zip(&self.dref)
                    .all(|(x, d)| x.to_ascii_uppercase() == *d)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_alleles_across_ref_lengths() {
        let mut v = Vcmp::default();
        assert!(v.set_ref(b"A", b"a"));
        assert_eq!(v.find_allele(&[b"C".as_slice(), b"G"], b"g"), Some(1));
        assert!(!v.set_ref(b"A", b"C"));
        // REF ACG vs A: the first record's alleles carry the extra CG.
        assert!(v.set_ref(b"ACG", b"A"));
        assert_eq!(v.find_allele(&[b"TCG".as_slice(), b"A"], b"T"), Some(0));
        assert_eq!(v.find_allele(&[b"A".as_slice()], b"T"), None);
    }
}
