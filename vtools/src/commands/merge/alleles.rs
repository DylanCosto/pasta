//! Allele lists and their merging (bcftools merge's `merge_alleles` and `normalize_alleles`).

/// A list of alleles, REF first, whose buffers are kept for reuse when it is cleared.
#[derive(Debug, Default)]
pub(super) struct Alleles {
    v: Vec<Vec<u8>>,
    n: usize,
}

impl Alleles {
    /// Empties the list, keeping the buffers.
    #[inline]
    pub(super) fn clear(&mut self) {
        self.n = 0;
    }

    /// Returns the number of alleles.
    #[inline]
    pub(super) fn len(&self) -> usize {
        self.n
    }

    /// Appends an allele.
    #[inline]
    pub(super) fn push(&mut self, a: &[u8]) {
        if self.n == self.v.len() {
            self.v.push(Vec::new());
        }
        let slot = &mut self.v[self.n];
        slot.clear();
        slot.extend_from_slice(a);
        self.n += 1;
    }

    /// Returns allele `i`.
    #[inline]
    pub(super) fn get(&self, i: usize) -> &[u8] {
        &self.v[..self.n][i]
    }

    /// Returns allele `i` for changes.
    #[inline]
    pub(super) fn get_mut(&mut self, i: usize) -> &mut Vec<u8> {
        &mut self.v[..self.n][i]
    }

    /// Returns the alleles.
    #[inline]
    pub(super) fn as_slice(&self) -> &[Vec<u8>] {
        &self.v[..self.n]
    }

    /// Replaces the list with the alleles of a record.
    pub(super) fn assign<'a>(&mut self, als: impl Iterator<Item = &'a [u8]>) {
        self.clear();
        for a in als {
            self.push(a);
        }
    }
}

/// Trims the bases shared by the ends of all alleles, keeping at least one (bcftools'
/// `normalize_alleles`): REF=GTTT ALT=GTT becomes REF=GT ALT=G.
pub(super) fn normalize_alleles(als: &mut Alleles) {
    let nals = als.len();
    if als.get(0).len() <= 1 {
        return;
    }
    let len0 = als.get(0).len();
    let mut i = 1;
    let mut done = false;
    while i < len0 {
        for j in 1..nals {
            let a = als.get(j);
            if i >= a.len() {
                done = true;
            }
            let differ = match a.len().checked_sub(i) {
                Some(k) => a[k] != als.get(0)[len0 - i],
                None => true,
            };
            if differ {
                done = true;
                break;
            }
        }
        if done {
            break;
        }
        i += 1;
    }
    if i > 1 {
        let i = i - 1;
        for j in 0..nals {
            let a = als.get_mut(j);
            let l = a.len();
            a.truncate(l - i);
        }
    }
}

/// REF length as bcftools' `merge_alleles` takes it (1 for one-base or empty alleles).
#[inline]
fn ref_len(a: &[u8]) -> usize {
    if a.len() <= 1 { 1 } else { a.len() }
}

/// Whether an allele is symbolic or an overlapping deletion, which REF expansion leaves alone.
#[inline]
fn is_symbolic_or_star(a: &[u8]) -> bool {
    matches!(a.first(), Some(b'<' | b'*'))
}

/// Why two allele lists could not be merged.
#[derive(Debug)]
pub(super) struct RefPrefixesDiffer;

/// Merges alleles `a` (REF first) into `b`, mapping each `a` allele to its index in `b` (bcftools'
/// `merge_alleles`). REFs of different lengths are reconciled by extending the shorter one's
/// alleles. When the REFs differ only in case, both lists are uppercased; returns whether `a`
/// was.
///
/// For example, `a` = ACG,AC,A and `b` = ACGT,A give `b` = ACGT,A,ACT,AT and map 0,2,3.
#[expect(
    clippy::needless_range_loop,
    reason = "a port of merge_alleles, indexing the alleles and map together"
)]
pub(super) fn merge_alleles(
    a: &mut Alleles,
    map: &mut Vec<usize>,
    b: &mut Alleles,
    tmp: &mut Vec<u8>,
) -> Result<bool, RefPrefixesDiffer> {
    let na = a.len();
    map.clear();
    map.resize(na, 0);
    let rla = ref_len(a.get(0));
    let rlb = ref_len(b.get(0));

    // The most common case: the same SNP.
    if na == 2
        && b.len() == 2
        && rla == 1
        && rlb == 1
        && a.get(1).len() == 1
        && b.get(1).len() == 1
        && a.get(1)[0] == b.get(1)[0]
    {
        map[1] = 1;
        return Ok(false);
    }

    // The REF prefixes must be identical.
    let n = rla.min(rlb);
    let (ra, rb) = (a.get(0), b.get(0));
    let (pa, pb) = (&ra[..n.min(ra.len())], &rb[..n.min(rb.len())]);
    let mut uppercased = false;
    if pa != pb {
        if !pa.eq_ignore_ascii_case(pb) {
            eprintln!(
                "The REF prefixes differ: {} vs {} ({rla},{rlb})",
                String::from_utf8_lossy(ra),
                String::from_utf8_lossy(rb)
            );
            return Err(RefPrefixesDiffer);
        }
        for i in 0..na {
            a.get_mut(i).make_ascii_uppercase();
        }
        for i in 0..b.len() {
            b.get_mut(i).make_ascii_uppercase();
        }
        uppercased = true;
    }

    // Extend the `b` alleles.
    if rla > rlb {
        for i in 0..b.len() {
            if is_symbolic_or_star(b.get(i)) {
                continue;
            }
            b.get_mut(i).extend_from_slice(&a.get(0)[rlb..]);
        }
    }

    // Add the `a` alleles that `b` lacks.
    for i in 1..na {
        tmp.clear();
        tmp.extend_from_slice(a.get(i));
        if rlb > rla && !is_symbolic_or_star(a.get(i)) {
            tmp.extend_from_slice(&b.get(0)[rla..]);
        }
        let nb = b.len();
        match (1..nb).find(|&j| tmp.eq_ignore_ascii_case(b.get(j))) {
            Some(j) => map[i] = j,
            None => {
                map[i] = nb;
                b.push(tmp);
            }
        }
    }
    Ok(uppercased)
}

/// Returns the index of the first unobserved allele (`<*>`, `<NON_REF>`, `<X>`), 0 when there
/// is none (bcftools' `get_unseen_allele`).
pub(super) fn unseen_allele<'a>(als: impl Iterator<Item = &'a [u8]>) -> usize {
    als.enumerate()
        .skip(1)
        .find(|(_, a)| is_unobserved(a))
        .map_or(0, |(i, _)| i)
}

/// Whether an allele stands for unobserved alleles (`<*>`, `<NON_REF>`, `<X>`).
#[inline]
pub(super) fn is_unobserved(a: &[u8]) -> bool {
    matches!(a, b"<*>" | b"<NON_REF>" | b"<X>")
}

/// The index of genotype `a/b` in Number=G order (HTSlib's `bcf_alleles2gt`).
#[inline]
pub(super) fn alleles2gt(a: usize, b: usize) -> usize {
    if a > b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn list(als: &[&str]) -> Alleles {
        let mut l = Alleles::default();
        for a in als {
            l.push(a.as_bytes());
        }
        l
    }

    fn strings(l: &Alleles) -> Vec<String> {
        l.as_slice()
            .iter()
            .map(|a| String::from_utf8_lossy(a).into_owned())
            .collect()
    }

    #[test]
    fn merges_deletions_of_different_ref_lengths() {
        let mut a = list(&["ACG", "AC", "A"]);
        let mut b = list(&["ACGT", "A"]);
        let mut map = Vec::new();
        let mut tmp = Vec::new();
        assert!(merge_alleles(&mut a, &mut map, &mut b, &mut tmp).is_ok());
        assert_eq!(strings(&b), ["ACGT", "A", "ACT", "AT"]);
        assert_eq!(map, [0, 2, 3]);
    }

    #[test]
    fn extends_b_and_keeps_symbolic_alleles() {
        let mut a = list(&["ACGT", "A"]);
        let mut b = list(&["A", "G", "<*>"]);
        let mut map = Vec::new();
        let mut tmp = Vec::new();
        assert!(merge_alleles(&mut a, &mut map, &mut b, &mut tmp).is_ok());
        assert_eq!(strings(&b), ["ACGT", "GCGT", "<*>", "A"]);
        assert_eq!(map, [0, 3]);
    }

    #[test]
    fn rejects_different_refs_and_uppercases_case_differences() {
        let mut map = Vec::new();
        let mut tmp = Vec::new();
        let mut a = list(&["C", "T"]);
        let mut b = list(&["A", "G"]);
        assert!(merge_alleles(&mut a, &mut map, &mut b, &mut tmp).is_err());
        let mut a = list(&["a", "t"]);
        let mut b = list(&["A", "G"]);
        assert!(merge_alleles(&mut a, &mut map, &mut b, &mut tmp).unwrap());
        assert_eq!(strings(&b), ["A", "G", "T"]);
    }

    #[test]
    fn normalizes_shared_suffixes() {
        let mut l = list(&["GTTT", "GTT"]);
        normalize_alleles(&mut l);
        assert_eq!(strings(&l), ["GT", "G"]);
        let mut l = list(&["ACT", "AT", "GCT"]);
        normalize_alleles(&mut l);
        assert_eq!(strings(&l), ["AC", "A", "GC"]);
        let mut l = list(&["A", "G"]);
        normalize_alleles(&mut l);
        assert_eq!(strings(&l), ["A", "G"]);
    }

    #[test]
    fn genotype_indexes() {
        assert_eq!(alleles2gt(0, 0), 0);
        assert_eq!(alleles2gt(1, 0), 1);
        assert_eq!(alleles2gt(1, 1), 2);
        assert_eq!(alleles2gt(0, 2), 3);
        assert_eq!(unseen_allele([&b"A"[..], b"G", b"<*>"].into_iter()), 2);
    }
}
