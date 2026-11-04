//! The genetic map and the transition hooks of `roh` (`load_genmap`, `get_genmap_rate`,
//! `set_tprob_genmap`, and `set_tprob_rrate` in `vcfroh.c`).

use std::path::Path;

use anyhow::bail;
use pasta::variant::{index::split_idx_name, log, number, readlist};

use super::hmm::{Matrix, TransitionHook};
use crate::raw_input::strerror;

/// The header line of an IMPUTE2 genetic map.
const GENMAP_HEADER: &[u8] = b"position COMBINED_rate(cM/Mb) Genetic_Map(cM)";

/// `vcfroh.c`'s `STATE_HW`, the Hardy-Weinberg state.
pub const HW: usize = 0;
/// `vcfroh.c`'s `STATE_AZ`, the autozygous state.
pub const AZ: usize = 1;

/// One genetic map point (`genmap_t`).
#[derive(Clone, Copy, Debug)]
pub struct MapPoint {
    /// 0-based position.
    pos: i32,
    /// Genetic distance from the start, in Morgans (cM scaled by 0.01).
    rate: f64,
}

/// The file name of a genetic map, `{CHROM}` replaced by the sequence name.
fn map_file_name(pattern: &str, chr: &[u8]) -> String {
    match pattern.find("{CHROM}") {
        Some(i) => format!(
            "{}{}{}",
            &pattern[..i],
            String::from_utf8_lossy(chr),
            &pattern[i + 7..]
        ),
        None => pattern.to_owned(),
    }
}

/// Loads the genetic map of a sequence (`load_genmap`). Returns `None` when the file cannot be
/// opened, which skips the sequence. `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
///
/// # Errors
///
/// Returns bcftools' errors for a file with the wrong header, unparsable lines, or no points,
/// and an error for a file that cannot be read (a truncated or corrupt gzip member), which
/// bcftools reports when closing it (`Close failed`).
pub fn load(pattern: &str, chr: &[u8]) -> anyhow::Result<Option<Vec<MapPoint>>> {
    let fname = map_file_name(pattern, chr);
    let data = split_idx_name(&fname).0;
    let mut reader = match readlist::open_text(Path::new(data)) {
        Ok(reader) => reader,
        Err(e) => {
            log::error(
                "hts_open_format",
                format_args!("Failed to open file \"{data}\" : {}", strerror(&e)),
            );
            return Ok(None);
        }
    };
    let lines = readlist::read_all_lines_from(&mut reader)
        .map_err(|e| anyhow::anyhow!("Error: failed to read {fname}: {e}"))?;
    let mut lines = lines.into_iter();
    let first = lines.next().unwrap_or_default();
    if first != GENMAP_HEADER {
        bail!(
            "Unexpected header in {fname}, found:\n\t[{}], but expected:\n\t[{}]",
            String::from_utf8_lossy(&first),
            String::from_utf8_lossy(GENMAP_HEADER)
        );
    }
    let mut map = Vec::new();
    // hts_getline stops at the first empty line.
    for line in lines.take_while(|l| !l.is_empty()) {
        let bad = || {
            anyhow::anyhow!(
                "Could not parse {fname}: {}",
                String::from_utf8_lossy(&line)
            )
        };
        let (pos, mut i) = number::strtol(&line, 0);
        if i == 0 {
            return Err(bad());
        }
        // Skip the separator and the second column.
        i += 1;
        while i < line.len() && !line[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= line.len() {
            return Err(bad());
        }
        let (rate, end) = number::strtod(&line, i + 1);
        if end == i + 1 {
            return Err(bad());
        }
        map.push(MapPoint {
            pos: (pos as i32).wrapping_sub(1),
            rate: rate * 0.01,
        });
    }
    if map.is_empty() {
        bail!("Genetic map empty?");
    }
    Ok(Some(map))
}

/// A genetic map with the search position that `get_genmap_rate` starts from
/// (`args->igenmap`).
#[derive(Clone, Copy, Debug)]
pub struct MapCursor<'a> {
    map: &'a [MapPoint],
    pub i: usize,
}

impl<'a> MapCursor<'a> {
    pub fn new(map: &'a [MapPoint], i: usize) -> Self {
        Self { map, i }
    }

    /// The probability of a recombination between two positions (`get_genmap_rate`).
    fn rate(&mut self, mut start: i32, mut end: i32) -> f64 {
        let gm = self.map;
        if gm.is_empty() {
            // bcftools keeps the old points of a sequence whose map is missing, but with no
            // points counted the search stays put and the rate is 0.
            return 0.0;
        }
        // i: the point at or before start.
        let mut i = self.i;
        if gm[i].pos > start {
            while i > 0 && gm[i].pos > start {
                i -= 1;
            }
        } else {
            while i + 1 < gm.len() && gm[i + 1].pos < start {
                i += 1;
            }
        }
        // j: the point at or after end.
        let mut j = i;
        while j + 1 < gm.len() && gm[j].pos < end {
            j += 1;
        }
        if i == j {
            self.i = i;
            return 0.0;
        }
        if start < gm[i].pos {
            start = gm[i].pos;
        }
        if end > gm[j].pos {
            end = gm[j].pos;
        }
        let rate =
            (gm[j].rate - gm[i].rate) / f64::from(gm[j].pos - gm[i].pos) * f64::from(end - start);
        self.i = j;
        rate
    }
}

/// Scales the switching probabilities by the probability of a crossover `ci`.
#[inline]
fn apply_crossover(ci: f64, tprob: &mut Matrix<2>) {
    let ci = if ci > 1.0 { 1.0 } else { ci };
    tprob[HW][AZ] *= ci;
    tprob[AZ][HW] *= ci;
    tprob[AZ][AZ] = 1.0 - tprob[HW][AZ];
    tprob[HW][HW] = 1.0 - tprob[AZ][HW];
}

/// How the transitions depend on the distance between sites.
#[derive(Clone, Copy, Debug)]
pub enum Transitions<'a> {
    /// Only through the precalculated matrix powers.
    Constant,
    /// Scaled by the genetic map, and by the recombination rate when it is not 0
    /// (`set_tprob_genmap`).
    GenMap {
        cursor: MapCursor<'a>,
        rec_rate: f64,
    },
    /// Scaled by a constant recombination rate per base (`set_tprob_rrate`).
    RecRate(f64),
}

impl Transitions<'_> {
    /// The genetic map's search position, 0 without a map.
    pub fn map_index(&self) -> usize {
        match self {
            Self::GenMap { cursor, .. } => cursor.i,
            _ => 0,
        }
    }

    /// Sets the genetic map's search position (`args->igenmap = ...`).
    pub fn set_map_index(&mut self, i: usize) {
        if let Self::GenMap { cursor, .. } = self {
            cursor.i = i;
        }
    }
}

impl TransitionHook<2> for Transitions<'_> {
    #[inline]
    fn set_tprob(&mut self, prev_pos: u32, pos: u32, tprob: &mut Matrix<2>) {
        match self {
            Self::Constant => {}
            Self::GenMap { cursor, rec_rate } => {
                let mut ci = cursor.rate(prev_pos as i32, pos as i32);
                if *rec_rate != 0.0 {
                    ci *= *rec_rate;
                }
                apply_crossover(ci, tprob);
            }
            Self::RecRate(rate) => {
                let ci = f64::from(pos.wrapping_sub(prev_pos)) * *rate;
                apply_crossover(ci, tprob);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn map(points: &[(i32, f64)]) -> Vec<MapPoint> {
        points
            .iter()
            .map(|&(pos, rate)| MapPoint { pos, rate })
            .collect()
    }

    #[test]
    fn idx_names_read_the_data_file() {
        let (dir, name) = crate::raw_input::write_idx_named(
            "genmap-idx",
            "gm.22.txt",
            "position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n1 1.0 0.0\n1000 1.0 1.0\n",
        );
        let pattern = name.replace("gm.22.txt", "gm.{CHROM}.txt");
        let points = load(&pattern, b"22");
        let _ = std::fs::remove_dir_all(&dir);
        let points = points
            .expect("the lines parse")
            .expect("the data file is read");
        assert_eq!(points.len(), 2);
    }

    #[test]
    fn file_names_substitute_the_sequence() {
        assert_eq!(map_file_name("gm.{CHROM}.txt", b"22"), "gm.22.txt");
        assert_eq!(map_file_name("gm.txt", b"22"), "gm.txt");
    }

    #[test]
    fn rate_interpolates_between_points() {
        let m = map(&[(0, 0.0), (100, 0.01), (200, 0.03)]);
        let mut c = MapCursor::new(&m, 0);
        assert_eq!(c.rate(50, 150), 0.03 / 200.0 * 100.0);
        assert_eq!(c.i, 2);
        // Within one interval of the map, nothing is gained.
        let mut c = MapCursor::new(&m, 2);
        assert_eq!(c.rate(250, 300), 0.0);
    }
}
