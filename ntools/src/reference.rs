//! Random access to an indexed (`.fai`) FASTA file.

use std::{
    collections::HashMap,
    ffi::OsString,
    path::{Path, PathBuf},
};

use anyhow::Context;
use pasta::fasta::{self as fasta, fai};

use crate::commands::faidx::{Source, build_index, read_bases};

/// Builds and saves a FASTA file's `.fai` index if it has none, as HTSlib's `fai_load` does.
pub fn ensure_index(path: &Path) -> anyhow::Result<()> {
    let mut index = OsString::from(path);
    index.push(".fai");

    if !Path::new(&index).exists() {
        build_index(path).with_context(|| format!("Could not build faidx: {}", path.display()))?;
    }

    Ok(())
}

/// An indexed FASTA file, plain or BGZF compressed (with its `.gzi`), read through its index.
pub struct IndexedFasta {
    source: Source,
    records: HashMap<Vec<u8>, fai::Record>,
}

impl IndexedFasta {
    pub fn open(path: &Path) -> anyhow::Result<Self> {
        ensure_index(path)?;
        let index = fasta::fs::read_associated_index(path)
            .with_context(|| format!("Could not load faidx: {}", path.display()))?;
        let records = index
            .as_ref()
            .iter()
            .map(|record| (record.name().to_vec(), record.clone()))
            .collect();
        let mut gzi = OsString::from(path);
        gzi.push(".gzi");
        let gzi = PathBuf::from(gzi);

        let mut source = Source::open(path, Some(&gzi))?;

        // A compressed file without its `.gzi` gets both indexes built, as HTSlib's `FAI_CREATE`
        // does; the `.fai` offsets are of the uncompressed text either way.
        if source.is_compressed() && !gzi.exists() {
            build_index(path)
                .with_context(|| format!("Could not build faidx: {}", path.display()))?;
            source = Source::open(path, Some(&gzi))?;
        }

        Ok(Self { source, records })
    }

    /// `faidx_fetch_seq64` of the 0-based closed interval `[beg, end]`, clamped to the sequence
    /// as htslib does. Returns `None` if the sequence is not in the index.
    pub fn fetch(
        &self,
        name: &[u8],
        beg: i64,
        end: i64,
        dst: &mut Vec<u8>,
    ) -> anyhow::Result<bool> {
        let Some(record) = self.records.get(name) else {
            return Ok(false);
        };

        let length = i64::try_from(record.length())?;
        let mut beg = beg;
        let mut end = end;

        if end < beg {
            beg = end;
        }

        beg = beg.clamp(0, length);
        end = if end < 0 { 0 } else { end.min(length - 1) };

        dst.clear();

        if end >= beg {
            let mut raw = Vec::new();
            read_bases(
                &self.source,
                record,
                usize::try_from(beg + 1)?,
                usize::try_from(end + 1)?,
                &mut raw,
                dst,
            )?;
        }

        Ok(true)
    }
}
