//! Whole-header operations: merging, combining, sample subsets, and read-time sample selection.

use super::{Header, HeaderError, HeaderRecordType, HeaderTextFlavor, KeepSamples};
use crate::variant::{log, readlist};

impl Header {
    /// Adds `src`'s lines that `self` lacks (HTSlib's `bcf_hdr_combine`). Generic lines are
    /// compared by key only; definitions by ID. Returns whether any definition had a different
    /// Number or Type (which is warned about).
    ///
    /// # Errors
    ///
    /// Returns an error when an added line's IDX conflicts.
    pub fn combine(&mut self, src: &Self) -> Result<bool, HeaderError> {
        self.merge_lines(src, MergeKind::Combine)
    }

    /// Merges `src` into `self` (HTSlib's `bcf_hdr_merge` with a destination): like
    /// [`Self::combine`], and the `##fileformat` becomes the newer of the two versions.
    ///
    /// # Errors
    ///
    /// Returns an error when an added line's IDX conflicts.
    pub fn merge(&mut self, src: &Self) -> Result<(), HeaderError> {
        self.merge_lines(src, MergeKind::Merge).map(|_| ())
    }

    /// Starts a merge from `src` (HTSlib's `bcf_hdr_merge(NULL, src)`): `src` reparsed without
    /// its `IDX` fields.
    ///
    /// # Errors
    ///
    /// Returns an error when the text does not reparse.
    pub fn merge_new(src: &Self) -> Result<Self, HeaderError> {
        let mut text = Vec::new();
        src.format_into(&mut text, HeaderTextFlavor::Vcf);
        Self::parse(&text)
    }

    fn merge_lines(&mut self, src: &Self, kind: MergeKind) -> Result<bool, HeaderError> {
        let ndst_ori: Vec<u32> = self.order.clone();
        let mut conflict = false;
        for rec in src.records() {
            match rec.kind {
                HeaderRecordType::Generic if rec.value.is_some() => {
                    let found = ndst_ori.iter().find_map(|&s| {
                        let d = self.record_by_serial(s)?;
                        (d.kind == HeaderRecordType::Generic && d.key == rec.key).then_some(s)
                    });
                    match found {
                        None => {
                            self.add_record(rec.dup())?;
                        }
                        Some(_) if kind == MergeKind::Merge && rec.key == b"fileformat" => {
                            if src.version() > self.version() {
                                let v = rec.value.clone().unwrap_or_default();
                                self.set_version(&v);
                            }
                        }
                        Some(_) => {}
                    }
                }
                HeaderRecordType::Structured => {
                    if let Some(j) = rec.find_key(b"ID")
                        && self
                            .get_record(
                                HeaderRecordType::Structured,
                                b"ID",
                                Some(&rec.vals[j]),
                                Some(&rec.key),
                            )
                            .is_none()
                    {
                        self.add_record(rec.dup())?;
                    }
                }
                HeaderRecordType::Generic => {}
                _ => {
                    let Some(j) = rec.find_key(b"ID") else {
                        continue;
                    };
                    let id = rec.vals[j].clone();
                    if self.get_record(rec.kind, b"ID", Some(&id), None).is_none() {
                        self.add_record(rec.dup())?;
                    } else if matches!(rec.kind, HeaderRecordType::Info | HeaderRecordType::Format)
                    {
                        let slot = rec.kind as usize;
                        let info_of = |h: &Self| {
                            h.id(&id)
                                .and_then(|i| h.ids[i as usize].as_ref())
                                .map_or(15, |e| e.info[slot])
                        };
                        let (s, d) = (info_of(src), info_of(self));
                        let name = String::from_utf8_lossy(&id).into_owned();
                        let func = match kind {
                            MergeKind::Combine => "bcf_hdr_combine",
                            MergeKind::Merge => "bcf_hdr_merge",
                        };
                        if (s >> 8) & 0xf != (d >> 8) & 0xf {
                            log::warn(
                                func,
                                format_args!(
                                    "Trying to combine \"{name}\" tag definitions of different lengths"
                                ),
                            );
                            conflict = true;
                        }
                        if (s >> 4) & 0xf != (d >> 4) & 0xf {
                            log::warn(
                                func,
                                format_args!(
                                    "Trying to combine \"{name}\" tag definitions of different types"
                                ),
                            );
                            conflict = true;
                        }
                    }
                }
            }
        }
        Ok(conflict)
    }

    /// Returns a header with only the given samples, in the given order (HTSlib's
    /// `bcf_hdr_subset`), and for each requested sample its index in `self` (or -1 when absent,
    /// in which case it is left out).
    ///
    /// # Errors
    ///
    /// Returns an error for a duplicate sample in `samples` or a malformed `#CHROM` line.
    pub fn subset(&self, samples: &[&[u8]]) -> Result<(Self, Vec<i32>), HeaderError> {
        let mut text = Vec::new();
        self.format_into(&mut text, HeaderTextFlavor::Bcf);
        let mut h = Self::new_for_writing();
        h.set_version(&self.version_string());
        let mut imap = vec![-1i32; samples.len()];
        let mut out: Vec<u8>;
        if self.sample_count() > 0 {
            let start = if text.starts_with(b"#CHROM\t") {
                0
            } else {
                text.windows(8)
                    .position(|w| w == b"\n#CHROM\t")
                    .map(|p| p + 1)
                    .ok_or(HeaderError::InvalidSampleLine)?
            };
            let end = if samples.is_empty() { 7 } else { 8 };
            let mut p = start;
            let mut i = 0;
            while let Some(t) = text[p..].iter().position(|&c| c == b'\t') {
                p += t;
                if i >= end {
                    break;
                }
                i += 1;
                p += 1;
            }
            if i != end {
                log::error(
                    "bcf_hdr_subset",
                    format_args!("Wrong number of columns in header #CHROM line"),
                );
                return Err(HeaderError::InvalidSampleLine);
            }
            out = text[..p].to_vec();
            let mut seen = rustc_hash::FxHashSet::default();
            for (k, &s) in samples.iter().enumerate() {
                if !seen.insert(s) {
                    let s = String::from_utf8_lossy(s).into_owned();
                    log::error(
                        "bcf_hdr_subset",
                        format_args!("Duplicate sample name \"{s}\""),
                    );
                    return Err(HeaderError::DuplicateSample(s));
                }
                imap[k] = self.sample_id(s).map_or(-1, |i| i as i32);
                if imap[k] < 0 {
                    continue;
                }
                out.push(b'\t');
                out.extend_from_slice(s);
            }
        } else {
            out = text;
        }
        while matches!(out.last(), Some(0 | b'\n')) {
            out.pop();
        }
        out.push(b'\n');
        h.parse_into(&out)?;
        Ok((h, imap))
    }

    /// Selects the samples to keep while reading (HTSlib's `bcf_hdr_set_samples`): a
    /// comma-separated list or a file of names (`^` to exclude them), `-` for all samples, or
    /// `None` for none. Records read afterwards hold only the kept samples.
    ///
    /// Returns 0, or the 1-based position in the list of the first name not in the header.
    ///
    /// # Errors
    ///
    /// Returns an error when the list file cannot be read.
    pub fn set_samples(
        &mut self,
        samples: Option<&[u8]>,
        source: readlist::ListSource,
    ) -> Result<usize, HeaderError> {
        if samples == Some(b"-") {
            return Ok(0);
        }
        let n_original = self.sample_count();
        let mut bits = vec![0u8; n_original / 8 + 1];
        let Some(samples) = samples else {
            self.keep_samples = Some(KeepSamples { bits, n_original });
            self.samples.clear();
            self.sample_map.clear();
            return Ok(0);
        };
        let exclude = samples.first() == Some(&b'^');
        if exclude {
            for i in 0..n_original {
                bits[i / 8] |= 1 << (i % 8);
            }
        }
        let list = readlist::read_list(if exclude { &samples[1..] } else { samples }, source)
            .map_err(|e| HeaderError::SampleList(e.to_string()))?;
        let mut ret = 0;
        for (i, name) in list.iter().enumerate() {
            let Some(idx) = self.sample_id(name) else {
                if ret == 0 {
                    ret = i + 1;
                }
                continue;
            };
            let idx = idx as usize;
            if exclude {
                bits[idx / 8] &= !(1 << (idx % 8));
            } else {
                bits[idx / 8] |= 1 << (idx % 8);
            }
        }
        let keep = KeepSamples { bits, n_original };
        let kept: Vec<Box<[u8]>> = (0..n_original)
            .filter(|&i| keep.contains(i))
            .map(|i| self.samples[i].clone())
            .collect();
        self.samples.clear();
        self.sample_map.clear();
        if kept.is_empty() {
            self.keep_samples = None;
        } else {
            for s in kept {
                let id = self.samples.len() as u32;
                self.sample_map.insert(s.clone(), id);
                self.samples.push(s);
            }
            self.keep_samples = Some(keep);
        }
        Ok(ret)
    }

    /// Returns whether read-time sample selection is active.
    #[must_use]
    pub fn has_sample_selection(&self) -> bool {
        self.keep_samples.is_some()
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MergeKind {
    Combine,
    Merge,
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    #[test]
    fn subset_reorders() {
        let h = Header::parse(TEXT.as_bytes()).unwrap();
        let (s, imap) = h.subset(&[b"C", b"X", b"A"]).unwrap();
        assert_eq!(imap, [2, -1, 0]);
        let names: Vec<&[u8]> = s.samples().iter().map(|b| &**b).collect();
        assert_eq!(names, [&b"C"[..], b"A"]);
        let (s, _) = h.subset(&[]).unwrap();
        assert!(
            String::from_utf8(s.to_vcf_text())
                .unwrap()
                .ends_with("\tINFO\n")
        );
    }

    #[test]
    fn combine_adds_missing() {
        let mut a = Header::parse(TEXT.as_bytes()).unwrap();
        let b = Header::parse(
            b"##fileformat=VCFv4.3\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        )
        .unwrap();
        a.merge(&b).unwrap();
        assert!(a.id(b"AF").is_some());
        assert_eq!(a.version(), 4_003_000);
    }
}
