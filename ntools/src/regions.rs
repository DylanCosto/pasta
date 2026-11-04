//! Region strings and BED files as samtools reads them, and merged genomic intervals for `-L`
//! and `-M`.

use std::{
    fs,
    io::{self, Read},
    path::Path,
};

use anyhow::{Context, bail};
use flate2::read::MultiGzDecoder;
use pasta::core::{Position, Region};
use pasta::sam;

use crate::hts_region;

/// Sorted, merged, 1-based closed intervals for each reference sequence (by header index).
#[derive(Clone, Debug, Default)]
pub struct RegionSet {
    by_reference: Vec<Vec<(usize, usize)>>,
}

impl RegionSet {
    pub fn new(header: &sam::Header) -> Self {
        Self {
            by_reference: vec![Vec::new(); header.reference_sequences().len()],
        }
    }

    /// Adds the intervals of a BED or position list file, read as samtools reads them (see
    /// [`read_bed`]). Names not in the header are ignored, as in samtools.
    pub fn add_bed(&mut self, header: &sam::Header, path: &Path) -> anyhow::Result<()> {
        for interval in read_bed(path)? {
            if let Some(id) = header
                .reference_sequences()
                .get_index_of(&interval.name[..])
                && interval.end > interval.start
            {
                let start = usize::try_from(interval.start)?;
                let end = usize::try_from(interval.end)?;
                self.by_reference[id].push((start + 1, end));
            }
        }

        Ok(())
    }

    /// Adds a samtools region string (see [`parse_region`]).
    pub fn add_region(&mut self, header: &sam::Header, raw: &str) -> anyhow::Result<()> {
        let Some((id, start, end)) = parse_region(header, raw) else {
            bail!("invalid region or unknown reference sequence: {raw}");
        };

        let end = end.unwrap_or_else(|| reference_sequence_length(header, id));

        if end >= start {
            self.by_reference[id].push((start, end));
        }

        Ok(())
    }

    /// Adds a 1-based closed interval of a reference sequence (by header index); call
    /// [`Self::merge`] afterwards.
    pub fn add_interval(&mut self, reference_sequence_id: usize, start: usize, end: usize) {
        if let Some(intervals) = self.by_reference.get_mut(reference_sequence_id) {
            intervals.push((start, end));
        }
    }

    /// The overlaps of this set's intervals with `other`'s, as samtools intersects `-L` with
    /// regions under `-M` (`bed_filter`).
    pub fn intersection(&self, other: &Self) -> Self {
        let by_reference = self
            .by_reference
            .iter()
            .zip(&other.by_reference)
            .map(|(a, b)| {
                let mut overlaps = Vec::new();

                for &(start, end) in a {
                    for &(s, e) in b {
                        if s <= end && e >= start {
                            overlaps.push((start.max(s), end.min(e)));
                        }
                    }
                }

                overlaps
            })
            .collect();

        let mut set = Self { by_reference };
        set.merge();
        set
    }

    /// Sorts and merges overlapping or adjacent intervals.
    pub fn merge(&mut self) {
        for intervals in &mut self.by_reference {
            intervals.sort_unstable();

            let mut merged: Vec<(usize, usize)> = Vec::with_capacity(intervals.len());

            for &(start, end) in intervals.iter() {
                match merged.last_mut() {
                    Some(last) if start <= last.1.saturating_add(1) => last.1 = last.1.max(end),
                    _ => merged.push((start, end)),
                }
            }

            *intervals = merged;
        }
    }

    /// Merged intervals of one reference sequence.
    pub fn intervals(&self, reference_sequence_id: usize) -> &[(usize, usize)] {
        self.by_reference
            .get(reference_sequence_id)
            .map_or(&[], Vec::as_slice)
    }

    /// Returns the index of the first merged interval that `[start, end]` overlaps.
    pub fn first_overlap(
        &self,
        reference_sequence_id: usize,
        start: usize,
        end: usize,
    ) -> Option<usize> {
        let intervals = self.intervals(reference_sequence_id);
        let i = intervals.partition_point(|&(_, e)| e < start);
        intervals.get(i).filter(|&&(s, _)| s <= end).map(|_| i)
    }
}

/// Returns a record's reference sequence ID and the 1-based closed interval samtools uses for
/// overlap tests: the reference span of a mapped record (at least one base), or one base at the
/// position of an unmapped one. Returns `None` for records with no reference or position.
pub fn record_interval(
    record: &dyn sam::alignment::Record,
    header: &sam::Header,
) -> io::Result<Option<(usize, usize, usize)>> {
    let (Some(id), Some(start)) = (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };

    let span = if record.flags()?.is_unmapped() {
        0
    } else {
        record.cigar().alignment_span()?
    };

    let start = usize::from(start);
    let end = start + span.max(1) - 1;

    Ok(Some((id, start, end)))
}

/// Formats a 1-based closed interval as a region query string.
pub fn region_string(name: &[u8], start: usize, end: usize) -> anyhow::Result<Region> {
    let name = std::str::from_utf8(name).context("reference sequence name is not UTF-8")?;
    let start = Position::try_from(start)?;
    let end = Position::try_from(end)?;
    Ok(Region::new(name, start..=end))
}

/// Returns the header index of the reference sequence named `name`, as HTSlib looks names up:
/// an exact `SN` match first, else the first reference sequence listing `name` in its `AN`
/// alternative names.
pub fn reference_sequence_index(header: &sam::Header, name: &[u8]) -> Option<usize> {
    use sam::header::record::value::map::reference_sequence::tag::ALTERNATIVE_NAMES;

    let reference_sequences = header.reference_sequences();

    reference_sequences.get_index_of(name).or_else(|| {
        reference_sequences.values().position(|reference_sequence| {
            reference_sequence
                .other_fields()
                .get(&ALTERNATIVE_NAMES)
                .is_some_and(|names| names.split(|&b| b == b',').any(|alias| alias == name))
        })
    })
}

/// Parses a samtools region string against a header as `hts_parse_region` does: `chr1`,
/// `chr1:100` (to the end), `chr1:100-`, `chr1:-200`, `chr1:1,000-2k`, and `{chr1:A}:5-6` for
/// names containing colons. Returns the reference sequence ID and the 1-based closed interval,
/// with `None` as the end for the end of the reference sequence, or `None` if the region is
/// invalid or names no reference sequence of the header.
pub fn parse_region(header: &sam::Header, raw: &str) -> Option<(usize, usize, Option<usize>)> {
    let (id, start, end) =
        hts_region::parse_region(raw, |name| reference_sequence_index(header, name))?;

    let start = usize::try_from(start).ok()? + 1;
    let end = if end == hts_region::POS_MAX {
        None
    } else {
        Some(usize::try_from(end).ok()?)
    };

    Some((id, start, end))
}

/// Parses a samtools region string against a header (see [`parse_region`]) into a query region
/// naming the header's reference sequence.
pub fn parse_query_region(header: &sam::Header, raw: &str) -> anyhow::Result<Region> {
    let (id, start, end) =
        parse_region(header, raw).with_context(|| format!("invalid region: {raw}"))?;

    let name = header
        .reference_sequences()
        .get_index(id)
        .map(|(name, _)| name)
        .expect("parse_region returns IDs of the header");
    let name = std::str::from_utf8(name).context("reference sequence name is not UTF-8")?;
    let start = Position::try_from(start)?;

    Ok(match end {
        Some(end) => Region::new(name, start..=Position::try_from(end)?),
        None if start == Position::MIN => Region::new(name, ..),
        None => Region::new(name, start..),
    })
}

fn reference_sequence_length(header: &sam::Header, id: usize) -> usize {
    header
        .reference_sequences()
        .get_index(id)
        .map_or(0, |(_, reference_sequence)| {
            reference_sequence.length().get()
        })
}

/// An interval of a BED or position list file: 0-based and half-open.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BedInterval {
    pub name: Vec<u8>,
    pub start: u64,
    pub end: u64,
}

/// Reads a BED file as samtools' `bed_read` does, for `view -L`, `depth -b`, `mpileup -l` and
/// the like: plain or gzip/BGZF-compressed (`-` is stdin), columns separated by any whitespace,
/// and two-column lines (`chr pos`) taken as 1-based positions. Blank lines and `#` comments
/// are skipped, as are `track` and `browser` lines that do not parse as intervals. Intervals
/// are returned in file order.
pub fn read_bed(path: &Path) -> anyhow::Result<Vec<BedInterval>> {
    let text = if path == Path::new("-") {
        let mut raw = Vec::new();
        io::stdin().lock().read_to_end(&mut raw)?;
        decompress(raw, path)?
    } else {
        read_text_file(path)?
    };

    let mut intervals = Vec::new();

    // As in `bed_read`, line numbers count the nonempty lines.
    for (n, line) in text
        .split(|&c| c == b'\n')
        .filter(|line| !line.is_empty())
        .enumerate()
    {
        let line = line.trim_ascii_start();

        if line.is_empty() || line[0] == b'#' {
            continue;
        }

        let name_len = line
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(line.len());
        let name = &line[..name_len];

        let (mut start, mut end, num) = if name_len < line.len() {
            scan_two_u64(&line[name_len + 1..])
        } else {
            (0, 0, 0)
        };

        if num == 1 {
            // A position, counting from 1.
            end = start;
            start = start.wrapping_sub(1);
        }

        if num < 1 || end < start {
            // Checked here rather than first, in case a reference sequence is called so.
            if name == b"browser" || name == b"track" {
                continue;
            }

            if num < 1 {
                bail!(
                    "[bed_read] Parse error reading \"{}\" at line {}",
                    path.display(),
                    n + 1
                );
            }

            bail!(
                "[bed_read] Parse error reading \"{}\" at line {} : end ({end}) must not be less \
                 than start ({start})",
                path.display(),
                n + 1
            );
        }

        intervals.push(BedInterval {
            name: name.to_vec(),
            start,
            end,
        });
    }

    Ok(intervals)
}

/// Reads a plain or gzip/BGZF-compressed text file whole, as HTSlib's `gzopen` and `hopen`
/// read BED files.
pub fn read_text_file(path: &Path) -> anyhow::Result<Vec<u8>> {
    let raw = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    decompress(raw, path)
}

/// Decompresses `raw` if it is gzip (BGZF included) and returns it as it is otherwise.
fn decompress(raw: Vec<u8>, path: &Path) -> anyhow::Result<Vec<u8>> {
    if !raw.starts_with(&[0x1f, 0x8b]) {
        return Ok(raw);
    }

    let mut text = Vec::new();
    MultiGzDecoder::new(&raw[..])
        .read_to_end(&mut text)
        .with_context(|| format!("reading {}", path.display()))?;
    Ok(text)
}

/// `sscanf(s, "%" SCNu64 " %" SCNu64, ...)`: returns the values and how many were read.
fn scan_two_u64(s: &[u8]) -> (u64, u64, u32) {
    fn scan(s: &[u8]) -> Option<(u64, &[u8])> {
        let s = s.trim_ascii_start();
        let (neg, s) = match s.first() {
            Some(b'-') => (true, &s[1..]),
            Some(b'+') => (false, &s[1..]),
            _ => (false, s),
        };

        let digits = s.iter().take_while(|c| c.is_ascii_digit()).count();

        if digits == 0 {
            return None;
        }

        let n = std::str::from_utf8(&s[..digits])
            .ok()
            .and_then(|d| d.parse::<u64>().ok())
            .unwrap_or(u64::MAX);

        Some((if neg { n.wrapping_neg() } else { n }, &s[digits..]))
    }

    match scan(s) {
        None => (0, 0, 0),
        Some((start, rest)) => match scan(rest) {
            None => (start, 0, 1),
            Some((end, _)) => (start, end, 2),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> sam::Header {
        "@SQ\tSN:sq0\tLN:1000\n@SQ\tSN:sq1\tLN:500\n@SQ\tSN:HLA:1\tLN:300\n@SQ\tSN:X:1\tLN:300\n@SQ\tSN:X\tLN:300\n"
            .parse()
            .expect("valid header")
    }

    fn bed_file(contents: &[u8]) -> io::Result<tempfile::NamedTempFile> {
        use std::io::Write;

        let mut file = tempfile::NamedTempFile::new()?;
        file.write_all(contents)?;
        Ok(file)
    }

    fn interval(name: &str, start: u64, end: u64) -> BedInterval {
        BedInterval {
            name: name.as_bytes().to_vec(),
            start,
            end,
        }
    }

    #[test]
    fn test_parse_region() {
        let header = header();

        assert_eq!(parse_region(&header, "sq0"), Some((0, 1, None)));
        assert_eq!(
            parse_region(&header, "sq0:1,000-1,001"),
            Some((0, 1000, Some(1001)))
        );
        assert_eq!(
            parse_region(&header, "sq0:1k-1001"),
            Some((0, 1000, Some(1001)))
        );
        assert_eq!(
            parse_region(&header, "{sq0}:100-200"),
            Some((0, 100, Some(200)))
        );
        assert_eq!(parse_region(&header, "sq1:100"), Some((1, 100, None)));
        assert_eq!(parse_region(&header, "sq1:100-"), Some((1, 100, None)));
        assert_eq!(parse_region(&header, "sq1:-100"), Some((1, 1, Some(100))));
        assert_eq!(parse_region(&header, "HLA:1"), Some((2, 1, None)));
        assert_eq!(parse_region(&header, "HLA:1:5-6"), Some((2, 5, Some(6))));
        assert_eq!(parse_region(&header, "{HLA:1}:5-6"), Some((2, 5, Some(6))));
        assert_eq!(parse_region(&header, "{X}:1"), Some((4, 1, None)));
        assert_eq!(parse_region(&header, "{X:1}"), Some((3, 1, None)));
        // Both `X:1` and `X` name reference sequences.
        assert_eq!(parse_region(&header, "X:1"), None);
        assert_eq!(parse_region(&header, "sq0:200-100"), None);
        assert_eq!(parse_region(&header, "sq2"), None);
    }

    #[test]
    fn test_parse_region_resolves_alternative_names() -> anyhow::Result<()> {
        // `AN` aliases resolve to their `SN`, and an exact `SN` takes precedence.
        let header: sam::Header =
            "@SQ\tSN:chr1\tLN:50\tAN:1,one,chr2\n@SQ\tSN:chr2\tLN:50\n".parse()?;

        assert_eq!(parse_region(&header, "1:1-4"), Some((0, 1, Some(4))));
        assert_eq!(parse_region(&header, "one:1-4"), Some((0, 1, Some(4))));
        assert_eq!(parse_region(&header, "chr2:1-4"), Some((1, 1, Some(4))));
        assert_eq!(parse_region(&header, "two:1-4"), None);
        Ok(())
    }

    #[test]
    fn test_parse_query_region() -> anyhow::Result<()> {
        let header = header();

        let region = parse_query_region(&header, "{HLA:1}:1,000-2k")?;
        assert_eq!(region.name(), b"HLA:1");
        assert_eq!(region.interval().start(), Some(Position::try_from(1000)?));
        assert_eq!(region.interval().end(), Some(Position::try_from(2000)?));

        assert!(parse_query_region(&header, "{sq0}:1x").is_err());

        Ok(())
    }

    #[test]
    fn test_read_bed() -> anyhow::Result<()> {
        let file = bed_file(
            b"# comment\ntrack name=x\nbrowser position sq0\n\n  sq0  10\t20 name\nsq1\t5\r\nsq0 0 0\n",
        )?;

        assert_eq!(
            read_bed(file.path())?,
            [
                interval("sq0", 10, 20),
                interval("sq1", 4, 5),
                interval("sq0", 0, 0)
            ]
        );

        Ok(())
    }

    #[test]
    fn test_read_bed_compressed() -> anyhow::Result<()> {
        use std::io::Write;

        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(b"sq0\t10\t20\n")?;
        let file = bed_file(&encoder.finish()?)?;

        assert_eq!(read_bed(file.path())?, [interval("sq0", 10, 20)]);

        Ok(())
    }

    #[test]
    fn test_read_text_file() -> anyhow::Result<()> {
        use std::io::Write;

        let text = b"sq0\t10\t20\tname\t\xff\n# comment\nsq1 5 6\n";

        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip.write_all(text)?;

        // Several BGZF blocks.
        let mut bgzf = pasta::bgzf::io::Writer::new(Vec::new());
        for line in text.split_inclusive(|&c| c == b'\n') {
            bgzf.write_all(line)?;
            bgzf.flush()?;
        }

        for contents in [text.to_vec(), gzip.finish()?, bgzf.finish()?] {
            let file = bed_file(&contents)?;
            assert_eq!(read_text_file(file.path())?, text);
        }

        Ok(())
    }

    #[test]
    fn test_read_bed_with_errors() -> anyhow::Result<()> {
        let file = bed_file(b"sq0\t10\t20\n\nsq0\t20\t10\n")?;
        let message = read_bed(file.path()).unwrap_err().to_string();
        assert!(
            message.ends_with("at line 2 : end (10) must not be less than start (20)"),
            "{message}"
        );

        let file = bed_file(b"sq0\tx\n")?;
        let message = read_bed(file.path()).unwrap_err().to_string();
        assert!(message.ends_with("at line 1"), "{message}");

        // A position of 0 is before the start.
        assert!(read_bed(bed_file(b"sq0\t0\n")?.path()).is_err());

        Ok(())
    }

    #[test]
    fn test_add_bed() -> anyhow::Result<()> {
        let header = header();
        let mut set = RegionSet::new(&header);

        let file = bed_file(b"sq0 99 200\nsq0 150\nsq3\t1\t2\nsq1\t10\t10\n")?;
        set.add_bed(&header, file.path())?;
        set.merge();

        assert_eq!(set.intervals(0), [(100, 200)]);
        assert!(set.intervals(1).is_empty());

        Ok(())
    }

    #[test]
    fn test_merge_and_first_overlap() -> anyhow::Result<()> {
        let header = header();
        let mut set = RegionSet::new(&header);

        set.add_region(&header, "sq0:100-200")?;
        set.add_region(&header, "sq0:150-250")?;
        set.add_region(&header, "sq0:251-260")?;
        set.add_region(&header, "sq0:400-500")?;
        set.add_region(&header, "sq1")?;
        set.add_region(&header, "{HLA:1}:1,000-2k")?;
        set.merge();

        assert_eq!(set.intervals(0), [(100, 260), (400, 500)]);
        assert_eq!(set.intervals(1), [(1, 500)]);
        assert_eq!(set.intervals(2), [(1000, 2000)]);

        assert_eq!(set.first_overlap(0, 1, 99), None);
        assert_eq!(set.first_overlap(0, 90, 100), Some(0));
        assert_eq!(set.first_overlap(0, 250, 450), Some(0));
        assert_eq!(set.first_overlap(0, 261, 399), None);
        assert_eq!(set.first_overlap(0, 300, 400), Some(1));
        assert_eq!(set.first_overlap(1, 500, 600), Some(0));

        Ok(())
    }
}
