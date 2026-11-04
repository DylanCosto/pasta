//! Number=G values of any ploidy. VCF orders the genotypes of ploidy `p` (alleles
//! `a1 <= a2 <= ... <= ap`) by the sum of `C(a_i + i - 1, i)`, which for two alleles is HTSlib's
//! `bcf_alleles2gt`. bcftools `norm` handles one and two alleles only; these give the other
//! ploidies the same treatment.

/// The binomial coefficient `C(n, k)`, or `None` when it does not fit in `usize`.
fn binomial(n: usize, k: usize) -> Option<usize> {
    if k > n {
        return Some(0);
    }
    let k = k.min(n - k);
    let mut c: u128 = 1;
    for i in 1..=k {
        // C(n - k + i, i) from C(n - k + i - 1, i - 1): exact at every step.
        c = c * (n - k + i) as u128 / i as u128;
        if c > usize::MAX as u128 {
            return None;
        }
    }
    usize::try_from(c).ok()
}

/// The number of genotypes of `ploidy` alleles drawn from `nals` alleles, `None` when it
/// does not fit in `usize`.
pub(super) fn genotype_count(nals: usize, ploidy: usize) -> Option<usize> {
    if nals == 0 {
        return Some(usize::from(ploidy == 0));
    }
    binomial(nals + ploidy - 1, ploidy)
}

/// The ploidy whose genotypes of `nals` alleles number `count`, for two alleles or more.
pub(super) fn ploidy_of(nals: usize, count: usize) -> Option<usize> {
    if nals < 2 {
        return None;
    }
    (1..)
        .map_while(|p| {
            genotype_count(nals, p)
                .filter(|&c| c <= count)
                .map(|c| (p, c))
        })
        .find_map(|(p, c)| (c == count).then_some(p))
}

/// The index of a genotype, its alleles in ascending order, in a Number=G list.
pub(super) fn genotype_index(alleles: &[usize]) -> usize {
    alleles
        .iter()
        .enumerate()
        .map(|(i, &a)| binomial(a + i, i + 1).expect("an index below the genotype count"))
        .sum()
}

/// The index of the genotype with `nalt` copies of allele `ialt` and REF otherwise, for
/// `ploidy` alleles.
pub(super) fn ref_alt_index(ialt: usize, ploidy: usize, nalt: usize) -> usize {
    (ploidy - nalt..ploidy)
        .map(|i| binomial(ialt + i, i + 1).expect("an index below the genotype count"))
        .sum()
}

/// For each genotype of ploidy `ploidy` of a record whose alleles go to `map` (the joined
/// alleles' indexes), in Number=G order, its index among the joined record's genotypes.
pub(super) fn genotype_map(map: &[usize], ploidy: usize, out: &mut Vec<usize>) {
    out.clear();
    let n = map.len();
    if n == 0 || ploidy == 0 {
        return;
    }
    let mut als = vec![0; ploidy];
    let mut mapped = vec![0; ploidy];
    loop {
        for (m, &a) in mapped.iter_mut().zip(&als) {
            *m = map[a];
        }
        mapped.sort_unstable();
        out.push(genotype_index(&mapped));
        // The next genotype: raise the first allele that can be raised (to stay at most the
        // next one, or the last allele) and reset the ones before it.
        let Some(i) = (0..ploidy).find(|&i| als[i] < als.get(i + 1).copied().unwrap_or(n - 1))
        else {
            return;
        };
        als[i] += 1;
        als[..i].fill(0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::norm::alleles2gt;

    #[test]
    fn counts_and_ploidies() {
        assert_eq!(genotype_count(2, 1), Some(2));
        assert_eq!(genotype_count(3, 2), Some(6));
        assert_eq!(genotype_count(3, 3), Some(10));
        assert_eq!(genotype_count(4, 3), Some(20));
        assert_eq!(genotype_count(2, 4), Some(5));
        assert_eq!(genotype_count(1000, 1000), None);
        assert_eq!(ploidy_of(2, 4), Some(3));
        assert_eq!(ploidy_of(3, 15), Some(4));
        assert_eq!(ploidy_of(3, 3), Some(1));
        assert_eq!(ploidy_of(3, 6), Some(2));
        assert_eq!(ploidy_of(3, 7), None);
        assert_eq!(ploidy_of(3, 1), None);
        assert_eq!(ploidy_of(1, 1), None);
    }

    #[test]
    fn genotypes_come_in_vcf_order() {
        // Triploid, three alleles: 000 001 011 111 002 012 112 022 122 222.
        let order = [
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
            [0, 0, 2],
            [0, 1, 2],
            [1, 1, 2],
            [0, 2, 2],
            [1, 2, 2],
            [2, 2, 2],
        ];
        for (i, g) in order.iter().enumerate() {
            assert_eq!(genotype_index(g), i, "{g:?}");
        }
        let mut map = Vec::new();
        genotype_map(&[0, 1, 2], 3, &mut map);
        assert_eq!(map, (0..10).collect::<Vec<_>>());
        assert_eq!(ref_alt_index(2, 3, 0), 0);
        assert_eq!(ref_alt_index(2, 3, 1), 4);
        assert_eq!(ref_alt_index(2, 3, 2), 7);
        assert_eq!(ref_alt_index(2, 3, 3), 9);
    }

    #[test]
    fn maps_match_the_diploid_and_haploid_ones() {
        let mut out = Vec::new();
        for map in [&[0, 1][..], &[0, 2], &[0, 3, 1], &[0, 2, 1, 4]] {
            genotype_map(map, 2, &mut out);
            let mut want = Vec::new();
            for i in 0..map.len() {
                for j in 0..=i {
                    want.push(alleles2gt(map[i], map[j]));
                }
            }
            assert_eq!(out, want, "{map:?}");
            genotype_map(map, 1, &mut out);
            assert_eq!(out, map);
            for (a, &m) in map.iter().enumerate().skip(1) {
                assert_eq!(ref_alt_index(a, 2, 1), alleles2gt(0, a));
                assert_eq!(ref_alt_index(m, 2, 2), alleles2gt(m, m));
            }
        }
        // A biallelic triploid joined as the second ALT: 000 002 022 222.
        genotype_map(&[0, 2], 3, &mut out);
        assert_eq!(out, [0, 4, 7, 9]);
    }
}
