//! Byte histograms for the entropy coders.
//!
//! Counting into one table, a run of a byte (or of a pair of bytes) makes each increment wait for
//! the last one to be stored. These count alternate positions into separate tables and add them
//! up, as htscodecs' `hist8` and `hist1_4` do.

const ALPHABET_SIZE: usize = 256;

/// Inputs at least this long are counted into four order-1 tables (htscodecs' threshold for its
/// two); below it, the other tables would cost more than they save.
const MIN_SPLIT_ORDER_1_LEN: usize = 500_000;

/// Counts the bytes of `src`.
pub(super) fn order_0(src: &[u8]) -> [u32; ALPHABET_SIZE] {
    let mut tables = [[0; ALPHABET_SIZE]; 4];
    let (chunks, remainder) = src.as_chunks::<4>();

    for chunk in chunks {
        for (table, &b) in tables.iter_mut().zip(chunk) {
            table[usize::from(b)] += 1;
        }
    }

    for &b in remainder {
        tables[0][usize::from(b)] += 1;
    }

    let [mut frequencies, t1, t2, t3] = tables;

    for (((f, a), b), c) in frequencies.iter_mut().zip(t1).zip(t2).zip(t3) {
        *f += a + b + c;
    }

    frequencies
}

/// Adds the count of each pair of adjacent bytes of `src` to `frequencies`, by the first byte,
/// then the second.
pub(super) fn add_order_1(src: &[u8], frequencies: &mut [[u32; ALPHABET_SIZE]; ALPHABET_SIZE]) {
    let Some((&first, rest)) = src.split_first() else {
        return;
    };

    if src.len() < MIN_SPLIT_ORDER_1_LEN {
        let mut prev = first;

        for &b in rest {
            frequencies[usize::from(prev)][usize::from(b)] += 1;
            prev = b;
        }

        return;
    }

    // Pairs go to four tables in turn.
    let mut others = vec![[[0; ALPHABET_SIZE]; ALPHABET_SIZE]; 3];
    let [t1, t2, t3] = &mut others[..] else {
        unreachable!();
    };

    let mut prev = first;
    let (quads, remainder) = rest.as_chunks::<4>();

    for &[a, b, c, d] in quads {
        frequencies[usize::from(prev)][usize::from(a)] += 1;
        t1[usize::from(a)][usize::from(b)] += 1;
        t2[usize::from(b)][usize::from(c)] += 1;
        t3[usize::from(c)][usize::from(d)] += 1;
        prev = d;
    }

    for &b in remainder {
        frequencies[usize::from(prev)][usize::from(b)] += 1;
        prev = b;
    }

    for table in &others {
        for (row, other_row) in frequencies.iter_mut().zip(table) {
            for (f, g) in row.iter_mut().zip(other_row) {
                *f += g;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn src(len: usize) -> Vec<u8> {
        // Runs and changes, as in quality scores.
        (0..len).map(|i| b"FFFF:FF,F"[(i * 7 / 3) % 9]).collect()
    }

    #[test]
    fn test_order_0() {
        for len in [0, 1, 3, 4, 5, 1001] {
            let src = src(len);
            let mut expected = [0; ALPHABET_SIZE];

            for &b in &src {
                expected[usize::from(b)] += 1;
            }

            assert_eq!(order_0(&src), expected, "len {len}");
        }
    }

    #[test]
    fn test_add_order_1() {
        for len in [0, 1, 2, 3]
            .into_iter()
            .chain((0..5).map(|n| MIN_SPLIT_ORDER_1_LEN + n))
        {
            let src = src(len);
            let mut expected = vec![[0; ALPHABET_SIZE]; ALPHABET_SIZE];

            for pair in src.windows(2) {
                expected[usize::from(pair[0])][usize::from(pair[1])] += 1;
            }

            let mut actual = [[0; ALPHABET_SIZE]; ALPHABET_SIZE];
            add_order_1(&src, &mut actual);
            assert!(actual[..] == expected[..], "len {len}");
        }
    }
}
