//! The reference sequence (`-f`): whole sequences read through the `.fai` index, as
//! `faidx_fetch_seq(fai, name, 0, INT_MAX)` returns them, and cached for the pileup.

use std::{
    cell::RefCell,
    collections::HashMap,
    io,
    path::Path,
    sync::{
        Arc, Mutex,
        atomic::{AtomicU64, Ordering},
    },
};

use pasta::{fasta::fai, variant::log};

use crate::faidx::{Faidx, LoadError};

/// An indexed FASTA file, plain or BGZF-compressed.
pub struct Fasta {
    faidx: Faidx,
    records: HashMap<Vec<u8>, fai::Record>,
}

impl Fasta {
    /// Opens a FASTA file, building and saving its index if it has none, as `fai_load` does.
    ///
    /// # Errors
    ///
    /// Returns an error, after printing HTSlib's messages, if the file or its index cannot
    /// be read or built.
    pub fn open(path: &Path) -> Result<Self, LoadError> {
        let faidx = Faidx::load(path)?;
        let records = faidx
            .index()
            .as_ref()
            .iter()
            .map(|r| (r.name().to_vec(), r.clone()))
            .collect();
        Ok(Self { faidx, records })
    }

    /// The whole sequence `name`, without line breaks, or `None` if it is not in the index.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read, or ends before the sequence does.
    pub fn fetch(&mut self, name: &[u8]) -> io::Result<Option<Vec<u8>>> {
        let Some(r) = self.records.get(name) else {
            return Ok(None);
        };
        let len = r.length();
        if len == 0 {
            return Ok(Some(Vec::new()));
        }
        let line_bases = r.line_base_count().get();
        let line_width = r.line_width().get();
        let offset = |p: u64| r.position() + p / line_bases * line_width + p % line_bases;
        let first = offset(0);
        let last = offset(len - 1);
        let mut raw = vec![0u8; (last - first + 1) as usize];
        self.faidx.read_at(first, &mut raw)?;
        // fai_retrieve keeps the printable characters. In a well-formed file those are the
        // first `line_bases` bytes of each line, so move whole lines together and check them.
        let (lb, lw) = (line_bases as usize, line_width as usize);
        let (mut out, mut at) = (0, 0);
        while at < raw.len() {
            let n = lb.min(raw.len() - at);
            raw.copy_within(at..at + n, out);
            out += n;
            at += lw;
        }
        raw.truncate(out);
        // Branch-free per block, so that the check vectorizes.
        let graphic = |block: &[u8]| {
            block
                .iter()
                .fold(0u8, |bad, &c| bad | u8::from(c.wrapping_sub(0x21) > 0x5d))
                == 0
        };
        if raw.len() as u64 != len || !raw.chunks(1 << 12).all(graphic) {
            // An irregular file: read it again and keep the printable characters.
            raw.resize((last - first + 1) as usize, 0);
            self.faidx.read_at(first, &mut raw)?;
            raw.retain(u8::is_ascii_graphic);
        }
        if raw.len() as u64 != len {
            return Err(io::ErrorKind::UnexpectedEof.into());
        }
        Ok(Some(raw))
    }
}

/// A cached sequence by contig ID.
type CachedSequence = (usize, Arc<Vec<u8>>);

/// Reference sequences by contig ID, loaded on demand and shared between threads.
pub struct References {
    /// Tells this set of references from others in [`LAST`].
    id: u64,
    names: Vec<Vec<u8>>,
    /// The FASTA file and the cached sequences, most recently used first; `None` without
    /// `-f`.
    cache: Option<Mutex<Cache>>,
}

struct Cache {
    fasta: Fasta,
    seqs: Vec<CachedSequence>,
}

/// How many sequences the cache keeps.
const SLOTS: usize = 4;

/// The source of [`References::id`].
static NEXT_ID: AtomicU64 = AtomicU64::new(0);

/// A sequence a thread looked up.
struct LastSequence {
    /// The [`References::id`] it was looked up in.
    refs: u64,
    tid: i32,
    seq: Arc<Vec<u8>>,
}

thread_local! {
    /// The sequence this thread looked up last. The pileup looks a sequence up for every read
    /// and column, so on many threads the shared cache's lock would be contended.
    static LAST: RefCell<Option<LastSequence>> = const { RefCell::new(None) };
}

impl References {
    pub fn new(fasta: Option<Fasta>, names: Vec<Vec<u8>>) -> Self {
        Self {
            id: NEXT_ID.fetch_add(1, Ordering::Relaxed),
            names,
            cache: fasta.map(|fasta| {
                Mutex::new(Cache {
                    fasta,
                    seqs: Vec::new(),
                })
            }),
        }
    }

    /// Whether a FASTA file was given.
    pub fn has_fasta(&self) -> bool {
        self.cache.is_some()
    }

    /// `mplp_get_ref`: the sequence of contig `tid`, or `None`. A sequence that is missing
    /// or cannot be read is reported, as HTSlib does, and looked up again the next time, as
    /// bcftools caches only the sequences it read.
    pub fn get(&self, tid: i32) -> Option<Arc<Vec<u8>>> {
        let last = LAST.with_borrow(|last| match last {
            Some(last) if last.refs == self.id && last.tid == tid => Some(Arc::clone(&last.seq)),
            _ => None,
        });
        if last.is_some() {
            return last;
        }
        let seq = self.get_shared(tid)?;
        LAST.set(Some(LastSequence {
            refs: self.id,
            tid,
            seq: Arc::clone(&seq),
        }));
        Some(seq)
    }

    /// [`Self::get`] from the cache shared by the threads.
    fn get_shared(&self, tid: i32) -> Option<Arc<Vec<u8>>> {
        let cache = self.cache.as_ref()?;
        let name = usize::try_from(tid).ok().and_then(|i| self.names.get(i))?;
        let tid = tid as usize;
        let mut cache = cache.lock().expect("lock poisoned");
        if let Some(i) = cache.seqs.iter().position(|&(t, _)| t == tid) {
            let entry = cache.seqs.remove(i);
            let seq = Arc::clone(&entry.1);
            cache.seqs.insert(0, entry);
            return Some(seq);
        }
        let seq = match cache.fasta.fetch(name) {
            Ok(Some(seq)) => Arc::new(seq),
            Ok(None) => {
                log::error(
                    "faidx_adjust_position",
                    format_args!(
                        "The sequence \"{}\" was not found",
                        String::from_utf8_lossy(name)
                    ),
                );
                return None;
            }
            Err(e) => {
                let reason = if e.kind() == io::ErrorKind::UnexpectedEof {
                    "unexpected end of file"
                } else {
                    "error reading file"
                };
                log::error(
                    "fai_retrieve",
                    format_args!("Failed to retrieve block: {reason}"),
                );
                return None;
            }
        };
        cache.seqs.insert(0, (tid, Arc::clone(&seq)));
        cache.seqs.truncate(SLOTS);
        Some(seq)
    }
}

#[cfg(test)]
mod tests {
    use std::{io::Write, path::PathBuf};

    use pasta::bgzf;

    use super::*;

    /// A new, empty directory for one test, which the test removes.
    fn test_dir(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("vtools-mpileup-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("creating the test directory");
        dir
    }

    /// Two sequences of 60 bases a line, the second spanning several BGZF blocks.
    fn fasta() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let a: Vec<u8> = (0..150).map(|i| b"ACGT"[i % 4]).collect();
        let b: Vec<u8> = (0..5000).map(|i| b"TTGCA"[i % 5]).collect();
        let mut text = Vec::new();
        for (name, seq) in [(&b"a"[..], &a), (b"b", &b)] {
            text.extend_from_slice(b">");
            text.extend_from_slice(name);
            text.push(b'\n');
            for line in seq.chunks(60) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
        }
        (text, a, b)
    }

    fn bgzip(data: &[u8]) -> Vec<u8> {
        let mut w = bgzf::io::Writer::new(Vec::new());
        for chunk in data.chunks(1000) {
            w.write_all(chunk).expect("writing to memory");
            w.flush().expect("writing to memory");
        }
        w.finish().expect("writing to memory")
    }

    #[test]
    fn reads_a_bgzipped_reference_with_or_without_its_indexes() {
        let dir = test_dir("bgzf-ref");
        let (text, a, b) = fasta();
        let gz = dir.join("ref.fa.gz");
        std::fs::write(&gz, bgzip(&text)).expect("writing the reference");
        // The first open builds .fai and .gzi; the second reads them.
        for round in 0..2 {
            let mut fasta = Fasta::open(&gz).expect("opening the reference");
            assert_eq!(
                fasta.fetch(b"a").expect("reading a"),
                Some(a.clone()),
                "{round}"
            );
            assert_eq!(
                fasta.fetch(b"b").expect("reading b"),
                Some(b.clone()),
                "{round}"
            );
            assert_eq!(fasta.fetch(b"c").expect("looking up c"), None, "{round}");
        }
        for index in ["ref.fa.gz.fai", "ref.fa.gz.gzi"] {
            assert!(dir.join(index).is_file(), "{index} is written");
        }
        // The CRAM decoder's repository reads the same file.
        let repo =
            super::super::input::reference_repository(Some(&gz)).expect("opening the repository");
        let seq = repo
            .get(b"b")
            .expect("b is in the repository")
            .expect("b reads");
        assert_eq!(AsRef::<[u8]>::as_ref(&*seq), &b[..]);
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }

    #[test]
    fn a_missing_or_unreadable_sequence_is_no_reference() {
        let dir = test_dir("short-ref");
        let (text, a, _) = fasta();
        let path = dir.join("ref.fa");
        std::fs::write(&path, &text).expect("writing the reference");
        let fasta = Fasta::open(&path).expect("opening the reference");
        // The file loses its end after indexing: reading b fails.
        std::fs::write(&path, &text[..text.len() / 2]).expect("truncating the reference");
        let names = vec![b"a".to_vec(), b"b".to_vec(), b"c".to_vec()];
        let refs = References::new(Some(fasta), names);
        assert_eq!(refs.get(0).as_deref(), Some(&a));
        assert_eq!(refs.get(1), None);
        assert_eq!(refs.get(2), None);
        assert_eq!(refs.get(3), None);
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }

    #[test]
    fn refuses_a_fifo_reference_without_waiting_for_a_writer() {
        let (dir, fifo) = crate::raw_input::make_fifo("mpileup-ref");
        let path = fifo.clone();
        let opened =
            crate::raw_input::within_seconds("opening a FASTA FIFO with no writer", move || {
                Fasta::open(&path).is_ok()
            });
        assert!(!opened);
        let mut fai = fifo.into_os_string();
        fai.push(".fai");
        assert!(
            !Path::new(&fai).exists(),
            "no index is written next to the FIFO"
        );
        std::fs::remove_dir_all(&dir).expect("removing the test directory");
    }
}
