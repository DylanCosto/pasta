//! The bucket layout of HTSlib's bin tables (`KHASH_MAP_INIT_INT(bin, bins_t)` in `hts.c`),
//! which sets the order `hts_idx_save` writes bins in.

/// khash's `__ac_HASH_UPPER`.
const HASH_UPPER: f64 = 0.77;

/// An insert-only khash table of bin numbers: keys are placed, and the table grows and
/// rehashes in place, exactly as `kh_put` and `kh_resize` do, so [`Self::keys`] visits them
/// in khash's order.
#[derive(Debug, Default)]
pub(super) struct BinTable {
    buckets: Vec<Option<u32>>,
    size: usize,
    upper_bound: usize,
}

impl BinTable {
    /// `kh_put`: grows a full table first (even when `key` is present), then adds `key` to the
    /// first empty bucket along its probe sequence unless it is there already.
    pub(super) fn put(&mut self, key: u32) {
        if self.size >= self.upper_bound {
            self.grow();
        }
        let mask = self.buckets.len() - 1;
        // `kh_int_hash_func` is the key itself.
        let mut i = key as usize & mask;
        let mut step = 0;
        while let Some(k) = self.buckets[i] {
            if k == key {
                return;
            }
            step += 1;
            i = (i + step) & mask;
        }
        self.buckets[i] = Some(key);
        self.size += 1;
    }

    /// `kh_resize` to the next power of two (at least 4) of a table without deleted buckets,
    /// with its in-place kick-out rehash: each key, in bucket order, takes the first free new
    /// bucket along its probe sequence, and a key not rehashed yet that sits there moves on in
    /// its place.
    fn grow(&mut self) {
        let new_n = (self.buckets.len() + 1).next_power_of_two().max(4);
        let mask = new_n - 1;
        let mut old = std::mem::take(&mut self.buckets);
        let mut buckets = vec![None; new_n];
        let mut taken = vec![false; new_n];
        for j in 0..old.len() {
            let Some(mut key) = old[j].take() else {
                continue;
            };
            loop {
                let mut i = key as usize & mask;
                let mut step = 0;
                while taken[i] {
                    step += 1;
                    i = (i + step) & mask;
                }
                taken[i] = true;
                let kicked = old.get_mut(i).and_then(Option::take);
                buckets[i] = Some(key);
                match kicked {
                    Some(k) => key = k,
                    None => break,
                }
            }
        }
        self.buckets = buckets;
        self.upper_bound = (new_n as f64 * HASH_UPPER + 0.5) as usize;
    }

    /// The keys in bucket order (`kh_begin` to `kh_end`).
    pub(super) fn keys(&self) -> impl Iterator<Item = u32> + '_ {
        self.buckets.iter().flatten().copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order(keys: &[u32]) -> Vec<u32> {
        let mut table = BinTable::default();
        for &k in keys {
            table.put(k);
        }
        table.keys().collect()
    }

    #[test]
    fn keys_sit_at_their_hash_modulo_the_table_size() {
        // Four buckets: 4681 & 3 = 1, 37450 & 3 = 2, 4682 & 3 = 2 (probes to 3).
        assert_eq!(order(&[4681, 37450, 4682]), [4681, 37450, 4682]);
        assert_eq!(order(&[4682, 4681]), [4681, 4682]);
    }

    #[test]
    fn a_put_of_a_present_key_grows_a_full_table() {
        // Three keys fill four buckets (upper bound 3): 11 probes past 3 to bucket 0. The
        // fourth put adds nothing but rehashes them into eight buckets.
        assert_eq!(order(&[3, 11, 1]), [11, 1, 3]);
        assert_eq!(order(&[3, 11, 1, 1]), [1, 11, 3]);
    }

    #[test]
    fn growing_kicks_out_keys_not_rehashed_yet() {
        // In four buckets 11 sits at 0, 14 at 2, and 6 at 3. Rehashed first, 11 takes bucket 3
        // and kicks 6 out to its own bucket, 6, before 14 is rehashed (which then probes to 7).
        assert_eq!(order(&[14, 6, 11, 11]), [11, 6, 14]);
        // The layout bcftools 1.24 writes for the bins of a VCF's second contig, the pseudo-bin
        // 37450 last and put twice.
        assert_eq!(
            order(&[4685, 4686, 4687, 4688, 37450, 37450]),
            [4688, 37450, 4685, 4686, 4687]
        );
    }
}
