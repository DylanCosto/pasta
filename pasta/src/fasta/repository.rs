//! Sequence repository and adapters.

mod adapter;
pub mod adapters;

pub use self::adapter::Adapter;

use std::{
    collections::HashMap,
    fmt, fs, io,
    path::{Path, PathBuf},
    sync::{Arc, Mutex, RwLock, Weak},
};

use crate::core::region::Interval;

use super::record::Sequence;

/// A normalized MD5 checksum, once calculated, locked while it is.
type ChecksumSlot = Arc<Mutex<Option<[u8; 16]>>>;

// The adapter and the cache have separate locks, so that reading from the adapter (e.g., a range
// read for one slice) does not block cache hits. To avoid deadlocks, the cache lock is never held
// while taking the adapter lock.
struct AdapterCache {
    adapter: Mutex<Box<dyn Adapter>>,
    cache: RwLock<HashMap<Vec<u8>, Arc<Sequence>>>,
    /// Sequences released from the cache, still returned while something else holds them.
    released: Mutex<Released>,
    /// Sequence lengths, by name, once known.
    lengths: RwLock<HashMap<Vec<u8>, usize>>,
    /// Normalized MD5 checksums of whole sequences, by name, once calculated. Each has its own
    /// lock, so that a checksum is calculated once while others are looked up.
    checksums: Mutex<HashMap<Vec<u8>, ChecksumSlot>>,
    /// The `M5` checksums of a sequence dictionary (see [`Repository::load_dictionary`]), by name.
    dictionary: RwLock<HashMap<Vec<u8>, [u8; 16]>>,
}

impl AdapterCache {
    fn cached(&self, name: &[u8]) -> Option<Arc<Sequence>> {
        if let Some(sequence) = self.cache.read().unwrap().get(name) {
            return Some(Arc::clone(sequence));
        }

        self.released
            .lock()
            .unwrap()
            .sequences
            .get(name)
            .and_then(Weak::upgrade)
    }
}

/// The sequences released from the cache (see [`Repository::release`]).
#[derive(Default)]
struct Released {
    sequences: HashMap<Vec<u8>, Weak<Sequence>>,
    /// The number of sequences at which those dropped since are removed. Removing them only
    /// once the map has doubled keeps releases constant time on average, even while many
    /// released sequences are held (e.g., by a container of records on many small contigs).
    prune_at: usize,
}

impl Released {
    /// The fewest sequences that are pruned.
    const MIN_PRUNE_AT: usize = 64;

    fn insert(&mut self, name: &[u8], sequence: Weak<Sequence>) {
        self.sequences.insert(name.into(), sequence);

        if self.sequences.len() >= self.prune_at {
            self.sequences
                .retain(|_, sequence| sequence.strong_count() > 0);
            self.prune_at = (2 * self.sequences.len()).max(Self::MIN_PRUNE_AT);
        }
    }
}

/// A caching sequence repository.
pub struct Repository(Arc<AdapterCache>);

impl Repository {
    /// Creates a sequence repository.
    pub fn new<A>(adapter: A) -> Self
    where
        A: Adapter + 'static,
    {
        Self(Arc::new(AdapterCache {
            adapter: Mutex::new(Box::new(adapter)),
            cache: RwLock::new(HashMap::new()),
            released: Mutex::new(Released::default()),
            lengths: RwLock::new(HashMap::new()),
            checksums: Mutex::new(HashMap::new()),
            dictionary: RwLock::new(HashMap::new()),
        }))
    }

    /// Reads the `M5` checksums of the sequence dictionary of the FASTA file `src`, if it has one,
    /// and returns the dictionary's path.
    ///
    /// The dictionary is `src` with its extension (after any `.gz` or `.bgz`) replaced by `.dict`,
    /// as Picard's `CreateSequenceDictionary` and GATK name it, or else `src` followed by `.dict`
    /// (`samtools dict -o`). Its `@SQ` lines give the checksums [`Self::dictionary_md5`] returns,
    /// trusted as GATK trusts them: they are not checked against the sequences.
    ///
    /// # Errors
    ///
    /// Returns an error if the dictionary exists but cannot be read, or has a malformed `M5`.
    pub fn load_dictionary(&self, src: &Path) -> io::Result<Option<PathBuf>> {
        use crate::sam::header::record::value::map::reference_sequence::Md5Checksum;

        let Some(path) = dictionary_paths(src)
            .into_iter()
            .find(|path| path.is_file())
        else {
            return Ok(None);
        };

        let text = fs::read_to_string(&path)?;
        let mut checksums = HashMap::new();

        for line in text.lines().filter(|line| line.starts_with("@SQ\t")) {
            let field = |tag: &str| {
                line.split('\t')
                    .find_map(|field| field.strip_prefix(tag)?.strip_prefix(':'))
            };

            if let (Some(name), Some(md5)) = (field("SN"), field("M5")) {
                let md5: Md5Checksum = md5.parse().map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("{}: invalid M5 for {name}: {md5}", path.display()),
                    )
                })?;
                checksums.insert(name.as_bytes().to_vec(), <[u8; 16]>::from(md5));
            }
        }

        *self.0.dictionary.write().unwrap() = checksums;

        Ok(Some(path))
    }

    /// Returns the checksum of the sequence with the given name in the sequence dictionary read
    /// by [`Self::load_dictionary`], if it has one. This never reads the sequence.
    pub fn dictionary_md5(&self, name: &[u8]) -> Option<[u8; 16]> {
        self.0.dictionary.read().unwrap().get(name).copied()
    }

    /// Returns the sequence of the given name.
    pub fn get(&self, name: &[u8]) -> Option<io::Result<Arc<Sequence>>> {
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(sequence));
        }

        let mut adapter = self.0.adapter.lock().unwrap();

        // A concurrent miss may have populated the cache while this call waited.
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(sequence));
        }

        let record = match adapter.get(name)? {
            Ok(record) => record,
            Err(e) => return Some(Err(e)),
        };

        self.0.released.lock().unwrap().sequences.remove(name);

        let mut cache = self.0.cache.write().unwrap();

        let s = cache
            .entry(name.into())
            .or_insert_with(|| Arc::new(record.into_sequence()));

        Some(Ok(Arc::clone(s)))
    }

    /// Returns whether the sequence of the given name is cached.
    pub(crate) fn is_cached(&self, name: &[u8]) -> bool {
        self.0.cache.read().unwrap().contains_key(name)
    }

    /// Removes the sequence of the given name from the cache. Until the last copy of it is
    /// dropped, [`Self::get`] still returns it; after that, it is loaded again.
    pub(crate) fn release(&self, name: &[u8]) {
        let sequence = self.0.cache.write().unwrap().remove(name);

        if let Some(sequence) = sequence {
            let mut released = self.0.released.lock().unwrap();
            released.insert(name, Arc::downgrade(&sequence));
        }
    }

    /// Returns the subsequence in the given interval of the sequence with the given name.
    ///
    /// If the whole sequence is not cached, this reads only the interval, when the adapter
    /// supports it (see [`Adapter::get_range`]), and does not cache it. The subsequence ends
    /// early if the interval extends past the end of the sequence.
    ///
    /// This returns `None` if the sequence is not cached and the adapter does not support
    /// ranges or does not have the sequence. Use [`Self::get`] in that case.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::core::Position;
    /// use pasta::fasta::{self as fasta, record::{Definition, Sequence}, Repository};
    ///
    /// let repository = Repository::new(vec![fasta::Record::new(
    ///     Definition::new("sq0", None),
    ///     Sequence::from(b"ACGT".to_vec()),
    /// )]);
    ///
    /// // The `Vec<Record>` adapter does not read ranges.
    /// let start = Position::try_from(2)?;
    /// let end = Position::try_from(3)?;
    /// assert!(repository.get_range(b"sq0", (start..=end).into()).is_none());
    ///
    /// // Once the whole sequence is cached, its subsequences are, too.
    /// repository.get(b"sq0").transpose()?;
    /// let sequence = repository.get_range(b"sq0", (start..=end).into()).transpose()?;
    /// assert_eq!(sequence.as_deref().map(|s| s.as_ref()), Some(&b"CG"[..]));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get_range(&self, name: &[u8], interval: Interval) -> Option<io::Result<Arc<Sequence>>> {
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(Arc::new(slice_to_end(&sequence, interval))));
        }

        let mut adapter = self.0.adapter.lock().unwrap();

        // The whole sequence may have been loaded while this call waited.
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(Arc::new(slice_to_end(&sequence, interval))));
        }

        // Ranges are not cached: a sorted file's slices read mostly disjoint spans, so a later
        // slice would rarely find its span already read.
        Some(adapter.get_range(name, interval)?.map(Arc::new))
    }

    /// Returns the length of the sequence with the given name.
    ///
    /// This does not read the sequence if it is cached or the adapter knows its length (see
    /// [`Adapter::len`]), e.g., from a FASTA index; otherwise, it reads and caches it, as
    /// [`Self::get`] does. Lengths are remembered.
    ///
    /// This returns `None` if the repository does not have the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::fasta::{self as fasta, record::{Definition, Sequence}, Repository};
    ///
    /// let repository = Repository::new(vec![fasta::Record::new(
    ///     Definition::new("sq0", None),
    ///     Sequence::from(b"ACGT".to_vec()),
    /// )]);
    ///
    /// assert_eq!(repository.sequence_len(b"sq0").transpose()?, Some(4));
    /// assert!(repository.sequence_len(b"sq1").is_none());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn sequence_len(&self, name: &[u8]) -> Option<io::Result<usize>> {
        if let Some(&len) = self.0.lengths.read().unwrap().get(name) {
            return Some(Ok(len));
        }

        let len = match self.0.cached(name) {
            Some(sequence) => sequence.len(),
            None => {
                let len = self.0.adapter.lock().unwrap().len(name);

                match len {
                    Some(Ok(len)) => len,
                    Some(Err(e)) => return Some(Err(e)),
                    None => match self.get(name)? {
                        Ok(sequence) => sequence.len(),
                        Err(e) => return Some(Err(e)),
                    },
                }
            }
        };

        self.0.lengths.write().unwrap().insert(name.into(), len);

        Some(Ok(len))
    }

    /// Returns the MD5 checksum of the whole sequence with the given name, normalized as `@SQ`
    /// `M5` fields hold it (see [`crate::cram::calculate_normalized_sequence_digest`]).
    ///
    /// The checksum is calculated once. A sequence that is not cached is read to calculate it but
    /// is not cached, so that checking a reference sequence does not keep all of it in memory.
    ///
    /// This returns `None` if the repository does not have the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::calculate_normalized_sequence_digest;
    /// use pasta::fasta::{self as fasta, record::{Definition, Sequence}, Repository};
    ///
    /// let repository = Repository::new(vec![fasta::Record::new(
    ///     Definition::new("sq0", None),
    ///     Sequence::from(b"acgt".to_vec()),
    /// )]);
    ///
    /// let checksum = repository.normalized_md5(b"sq0").transpose()?;
    /// assert_eq!(checksum, Some(calculate_normalized_sequence_digest(b"ACGT")));
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn normalized_md5(&self, name: &[u8]) -> Option<io::Result<[u8; 16]>> {
        use crate::cram::calculate_normalized_sequence_digest;

        let slot = Arc::clone(
            self.0
                .checksums
                .lock()
                .unwrap()
                .entry(name.into())
                .or_default(),
        );
        // Another call calculating this checksum holds the slot until it is done.
        let mut slot = slot.lock().unwrap();

        if let Some(checksum) = *slot {
            return Some(Ok(checksum));
        }

        let checksum = match self.0.cached(name) {
            Some(sequence) => calculate_normalized_sequence_digest(&sequence[..]),
            None => {
                let mut adapter = self.0.adapter.lock().unwrap();

                // The whole sequence may have been loaded while this call waited.
                if let Some(sequence) = self.0.cached(name) {
                    drop(adapter);
                    calculate_normalized_sequence_digest(&sequence[..])
                } else {
                    let result = adapter.get(name)?;
                    // Other sequences are read while this one is hashed.
                    drop(adapter);

                    match result {
                        Ok(record) => {
                            calculate_normalized_sequence_digest(record.sequence().as_ref())
                        }
                        Err(e) => return Some(Err(e)),
                    }
                }
            }
        };

        *slot = Some(checksum);

        Some(Ok(checksum))
    }

    /// Returns the number of cached sequences.
    pub fn len(&self) -> usize {
        self.0.cache.read().unwrap().len()
    }

    /// Returns whether any sequences are cached.
    pub fn is_empty(&self) -> bool {
        self.0.cache.read().unwrap().is_empty()
    }

    /// Clears the sequence cache, and the lengths and checksums remembered.
    pub fn clear(&self) {
        self.0.cache.write().unwrap().clear();
        *self.0.released.lock().unwrap() = Released::default();
        self.0.lengths.write().unwrap().clear();
        self.0.checksums.lock().unwrap().clear();
    }
}

/// The paths a sequence dictionary of the FASTA file `src` may have, in the order they are tried.
fn dictionary_paths(src: &Path) -> Vec<PathBuf> {
    let uncompressed = match src.extension().and_then(|ext| ext.to_str()) {
        Some("gz" | "bgz") => src.with_extension(""),
        _ => src.to_path_buf(),
    };

    let mut appended = src.as_os_str().to_owned();
    appended.push(".dict");

    vec![uncompressed.with_extension("dict"), PathBuf::from(appended)]
}

/// Returns the subsequence in the interval, ending early at the end of the sequence.
fn slice_to_end(sequence: &Sequence, interval: Interval) -> Sequence {
    let src = sequence.as_ref();

    let start = interval.start().map(|p| usize::from(p) - 1).unwrap_or(0);
    let end = interval.end().map(usize::from).unwrap_or(src.len());

    let start = start.min(src.len());
    let end = end.clamp(start, src.len());

    Sequence::from(src[start..end].to_vec())
}

impl Clone for Repository {
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl fmt::Debug for Repository {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Repository")
            .field("cache", &*self.0.cache.read().unwrap())
            .finish()
    }
}

impl Default for Repository {
    fn default() -> Self {
        Self::new(adapters::Empty::new())
    }
}

#[cfg(test)]
mod tests {
    use std::{
        sync::{Barrier, mpsc},
        time::Duration,
    };

    use super::*;
    use crate::fasta::{
        Record,
        record::{Definition, Sequence},
    };

    #[test]
    fn test_get() -> io::Result<()> {
        let sq0 = Record::new(
            Definition::new("sq0", None),
            Sequence::from(b"ACGT".to_vec()),
        );
        let repository = Repository::new(vec![sq0.clone()]);

        assert_eq!(
            repository.get(b"sq0").transpose()?,
            Some(Arc::new(sq0.sequence().clone()))
        );
        assert_eq!(repository.get(b"sq1").transpose()?, None);

        Ok(())
    }

    // Loads every name as "ACGT", counting loads. Range reads wait at the barrier when they start
    // and before they load.
    struct TestAdapter(Arc<Mutex<usize>>, Arc<Barrier>);

    impl Adapter for TestAdapter {
        fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>> {
            *self.0.lock().unwrap() += 1;
            std::thread::sleep(Duration::from_millis(10));
            let sequence = Sequence::from(b"ACGT".to_vec());
            Some(Ok(Record::new(Definition::new(name, None), sequence)))
        }

        fn get_range(&mut self, name: &[u8], _: Interval) -> Option<io::Result<Sequence>> {
            self.1.wait();
            self.1.wait();
            Some(self.get(name)?.map(Record::into_sequence))
        }
    }

    #[test]
    fn test_release() -> io::Result<()> {
        let adapter = TestAdapter(Arc::default(), Arc::new(Barrier::new(1)));
        let loads = Arc::clone(&adapter.0);
        let repository = Repository::new(adapter);

        let sequence = repository.get(b"sq0").transpose()?;
        repository.release(b"sq0");
        assert!(!repository.is_cached(b"sq0"));

        // Still held, so it is not loaded again.
        let held = repository.get(b"sq0").transpose()?;
        assert!(Arc::ptr_eq(
            sequence.as_ref().unwrap(),
            held.as_ref().unwrap()
        ));
        assert_eq!(*loads.lock().unwrap(), 1);

        drop((sequence, held));
        repository.get(b"sq0").transpose()?;
        assert_eq!(*loads.lock().unwrap(), 2);
        assert!(repository.is_cached(b"sq0"));

        Ok(())
    }

    #[test]
    fn test_concurrent_cache_misses_load_once() {
        let adapter = TestAdapter(Arc::default(), Arc::new(Barrier::new(2)));
        let loads = Arc::clone(&adapter.0);
        let repository = Repository::new(adapter);
        std::thread::scope(|scope| {
            for _ in 0..4 {
                scope.spawn(|| repository.get(b"sq0").unwrap().unwrap());
            }
        });
        assert_eq!(*loads.lock().unwrap(), 1);
    }

    #[test]
    fn test_cache_hits_do_not_wait_for_adapter_reads() {
        let barrier = Arc::new(Barrier::new(2));
        let repository = &Repository::new(TestAdapter(Arc::default(), Arc::clone(&barrier)));
        repository.get(b"sq0").unwrap().unwrap();
        std::thread::scope(|scope| {
            let range_read = scope.spawn(move || repository.get_range(b"sq1", (..).into()));
            barrier.wait();
            let (tx, rx) = mpsc::channel();
            scope.spawn(move || tx.send(repository.get(b"sq0").is_some()).unwrap());
            let result = rx.recv_timeout(Duration::from_secs(10));
            barrier.wait();
            assert_eq!(result, Ok(true), "a cache hit waited for an adapter read");
            assert!(range_read.join().unwrap().is_some());
        });
    }

    #[test]
    fn test_get_range_matches_whole_sequence() -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Cursor;

        use crate::core::Region;

        use crate::fasta::{fai, io::IndexedReader};

        let data = Cursor::new(b">sq0\nACGT\nACGT\nAC\n>sq1\nNN\n");
        let mut index = fai::io::Reader::new(&b"sq0\t10\t5\t4\t5\nsq1\t2\t23\t2\t3\n"[..]);
        let reader = IndexedReader::new(data, index.read_index()?);
        let ranges = Repository::new(adapters::IndexedReader::new(reader));
        let sequence = Sequence::from(b"ACGTACGTAC".to_vec());
        let whole = Repository::new(vec![Record::new(Definition::new("sq0", None), sequence)]);
        whole.get(b"sq0").transpose()?;

        // Regions and their sequences. Ranges end early at the end of the sequence.
        for case in "sq0:3-6=GTAC sq0:8-12=TAC sq0:6=CGTAC sq0:11= sq0:13-14=".split(' ') {
            let (region, expected) = case.split_once('=').unwrap();
            let region: Region = region.parse()?;
            for r in [&ranges, &whole] {
                let actual = r.get_range(region.name(), region.interval()).transpose()?;
                assert_eq!(
                    actual.as_deref().map(|s| s.as_ref()),
                    Some(expected.as_bytes())
                );
            }
        }

        assert!(ranges.is_empty(), "ranges are not cached");
        // A missing sequence is `None` for both adapters.
        assert!(ranges.get(b"sqx").is_none());
        assert!(ranges.get_range(b"sqx", (..).into()).is_none());
        assert!(whole.get_range(b"sqx", (..).into()).is_none());

        Ok(())
    }

    #[test]
    fn test_sequence_len_and_normalized_md5_do_not_cache_sequences()
    -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Cursor;

        use crate::cram::calculate_normalized_sequence_digest;
        use crate::fasta::{fai, io::IndexedReader};

        let data = Cursor::new(b">sq0\nACGT\nacgt\nAC\n>sq1\nNN\n");
        let mut index = fai::io::Reader::new(&b"sq0\t10\t5\t4\t5\nsq1\t2\t23\t2\t3\n"[..]);
        let reader = IndexedReader::new(data, index.read_index()?);
        let indexed = Repository::new(adapters::IndexedReader::new(reader));
        let records = Repository::new(vec![Record::new(
            Definition::new("sq0", None),
            Sequence::from(b"ACGTacgtAC".to_vec()),
        )]);

        let md5 = calculate_normalized_sequence_digest(b"ACGTACGTAC");

        for repository in [&indexed, &records] {
            assert_eq!(repository.sequence_len(b"sq0").transpose()?, Some(10));
            assert_eq!(repository.normalized_md5(b"sq0").transpose()?, Some(md5));
            assert!(repository.sequence_len(b"sqx").is_none());
            assert!(repository.normalized_md5(b"sqx").is_none());
        }

        // The index gives the length, and the checksum is calculated without keeping the sequence.
        assert!(indexed.is_empty());
        assert_eq!(indexed.sequence_len(b"sq1").transpose()?, Some(2));
        assert!(indexed.is_empty());

        Ok(())
    }
}
