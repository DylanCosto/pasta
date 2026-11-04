use std::{collections::HashMap, io};

use crate::core::Region;

use super::Record;

/// A FASTA index.
///
/// Records are looked up by name through a map built once, so lookups take constant time
/// whatever the number of sequences (e.g., a transcriptome's hundreds of thousands). As in
/// HTSlib, the first record with a name is the one found.
#[derive(Clone, Debug, Default)]
pub struct Index {
    records: Vec<Record>,
    /// The position in `records` of the first record with each name.
    positions: HashMap<Vec<u8>, usize>,
}

impl Index {
    /// Returns start position of the given region.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use pasta::core::Region;
    /// use pasta::fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let index = fai::Index::from(vec![
    ///     fai::Record::new("sq0", 10946, 4, line_base_count, line_width),
    /// ]);
    ///
    /// let region = Region::new("sq0", ..);
    /// assert_eq!(index.query(&region)?, 4);
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn query(&self, region: &Region) -> io::Result<u64> {
        self.get(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", region.name()),
                )
            })
            .and_then(|record| record.query(region.interval()))
    }

    /// Returns the record of the sequence named `name`, if any.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use pasta::fasta::fai;
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let index = fai::Index::from(vec![
    ///     fai::Record::new("sq0", 10946, 4, line_base_count, line_width),
    /// ]);
    ///
    /// assert_eq!(index.get(b"sq0").map(|record| record.length()), Some(10946));
    /// assert!(index.get(b"sq1").is_none());
    /// ```
    pub fn get(&self, name: &[u8]) -> Option<&Record> {
        self.positions.get(name).map(|&i| &self.records[i])
    }
}

impl PartialEq for Index {
    fn eq(&self, other: &Self) -> bool {
        // The positions follow from the records.
        self.records == other.records
    }
}

impl Eq for Index {}

impl AsRef<[Record]> for Index {
    fn as_ref(&self) -> &[Record] {
        &self.records
    }
}

impl From<Vec<Record>> for Index {
    fn from(records: Vec<Record>) -> Self {
        let mut positions = HashMap::with_capacity(records.len());

        for (i, record) in records.iter().enumerate() {
            positions.entry(record.name().to_vec()).or_insert(i);
        }

        Self { records, positions }
    }
}

impl From<Index> for Vec<Record> {
    fn from(index: Index) -> Self {
        index.records
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use super::*;

    #[test]
    fn test_get_finds_the_first_record_with_a_name_among_many() {
        let line_base_count = const { NonZero::new(80).unwrap() };
        let line_width = const { NonZero::new(81).unwrap() };

        let mut records: Vec<_> = (0..100_000u64)
            .map(|i| Record::new(format!("ctg{i}"), 100, i * 200, line_base_count, line_width))
            .collect();
        // A duplicate name: the first record is the one found, as with a scan.
        records.push(Record::new("ctg7", 1, 0, line_base_count, line_width));

        let index = Index::from(records);

        assert_eq!(
            index.get(b"ctg99999").map(Record::position),
            Some(99_999 * 200)
        );
        assert_eq!(index.get(b"ctg7").map(Record::length), Some(100));
        assert!(index.get(b"ctg100000").is_none());

        let region = Region::new("ctg12345", ..);
        assert_eq!(index.query(&region).ok(), Some(12_345 * 200));
    }
}
