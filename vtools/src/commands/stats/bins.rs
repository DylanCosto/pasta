//! Histograms of bcftools stats: the log-scaled `dist_t` (`dist.c`), the fixed-width `idist_t`,
//! and the allele-frequency bins of `--af-bins` (`bin.c`).

use anyhow::bail;
use pasta::variant::{
    number::strtod,
    readlist::{self, ListSource},
};

use crate::zeroed::{self, AllocError};

/// A histogram that is exact up to `10^npow` and log-scaled above (bcftools' `dist_t`).
#[derive(Clone, Debug)]
pub struct Dist {
    bins: Vec<u64>,
    npow: i32,
    nexact: u32,
    nlevel: u32,
}

impl Dist {
    /// Creates a histogram exact up to `10^npow` (`dist_init`).
    pub fn new(npow: i32) -> Self {
        let nexact = 10f64.powi(npow) as u32;
        let nlevel = (f64::from(nexact) - 10f64.powi(npow - 1)) as u32;
        Self {
            bins: Vec::new(),
            npow,
            nexact,
            nlevel,
        }
    }

    /// Counts a value (`dist_insert`).
    pub fn insert(&mut self, value: u32) {
        let ibin = if value <= self.nexact {
            value as usize
        } else {
            let npow = f64::from(value).log10() as i32;
            let level = npow - self.npow + 1;
            let step = f64::from(10f64.powi(level) as u32);
            let fine = (f64::from(value) - 10f64.powi(npow)) / step;
            (f64::from(self.nexact) + f64::from(self.nlevel) * f64::from(level - 1) + fine) as usize
        };
        if ibin >= self.bins.len() {
            self.bins.resize(ibin + 1, 0);
        }
        self.bins[ibin] += 1;
    }

    /// Returns the number of bins (`dist_nbins`).
    pub fn len(&self) -> usize {
        self.bins.len()
    }

    /// Returns the count of bin `idx` and the smallest value it holds (`dist_get`).
    pub fn get(&self, idx: usize) -> (u64, u32) {
        let beg = if (idx as u64) < u64::from(self.nexact) {
            idx as u32
        } else {
            let idx = idx as u32;
            let level = (idx - self.nexact) / self.nlevel + 1;
            let bin = idx - self.nexact - self.nlevel * (level - 1);
            let step = 10f64.powi(level as i32) as u32;
            (10f64.powi(level as i32 + self.npow - 1) + f64::from(step) * f64::from(bin)) as u32
        };
        (self.bins[idx], beg)
    }

    /// Adds another histogram's counts.
    pub fn merge(&mut self, other: &Self) {
        if other.bins.len() > self.bins.len() {
            self.bins.resize(other.bins.len(), 0);
        }
        for (a, b) in self.bins.iter_mut().zip(&other.bins) {
            *a += b;
        }
    }
}

/// The `-d/--depth` binning: `min,max,step`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DepthBins {
    /// Values below go to the first bin.
    pub min: i32,
    /// Values above go to the last bin.
    pub max: i32,
    /// The bin width.
    pub step: i32,
}

impl Default for DepthBins {
    fn default() -> Self {
        Self {
            min: 0,
            max: 500,
            step: 1,
        }
    }
}

/// A fixed-width histogram with underflow and overflow bins (bcftools' `idist_t`).
#[derive(Clone, Debug)]
pub struct IDist {
    /// The binning.
    pub bins: DepthBins,
    /// The counts: underflow, `(max - min) / step + 2` bins, and overflow.
    pub vals: Vec<u64>,
}

impl IDist {
    /// Creates an empty histogram (`idist_init`).
    ///
    /// # Errors
    ///
    /// Fails when the system cannot allocate the bins (a huge `-d` range).
    pub fn try_new(bins: DepthBins) -> Result<Self, AllocError> {
        let m_vals = 4 + (bins.max - bins.min) / bins.step;
        Ok(Self {
            bins,
            vals: zeroed::counts(m_vals as usize)?,
        })
    }

    /// Counts a value (`idist`).
    #[inline]
    pub fn add(&mut self, val: i32) {
        let i = if val < self.bins.min {
            0
        } else if val > self.bins.max {
            self.vals.len() - 1
        } else {
            1 + ((val - self.bins.min) / self.bins.step) as usize
        };
        self.vals[i] += 1;
    }

    /// Returns the smallest value of bin `i`. bcftools' `idist_i2bin` returns `i - 1 + min`,
    /// which is wrong when the step is not 1.
    pub fn bin_start(&self, i: usize) -> i64 {
        i64::from(self.bins.min) + (i as i64 - 1) * i64::from(self.bins.step)
    }

    /// Adds another histogram's counts.
    pub fn merge(&mut self, other: &Self) {
        for (a, b) in self.vals.iter_mut().zip(&other.vals) {
            *a += b;
        }
    }
}

/// Allele-frequency bin boundaries from `--af-bins` (bcftools' `bin_t`).
#[derive(Clone, Debug)]
pub struct AfBins {
    bins: Vec<f32>,
}

impl AfBins {
    /// Reads boundaries from a comma-separated list or, without a comma, a file, and adds 0
    /// and 1 when missing (`bin_init(list, 0, 1)`).
    pub fn parse(list_def: &str) -> anyhow::Result<Self> {
        let (min, max) = (0f32, 1f32);
        let source = if list_def.contains(',') {
            ListSource::Inline
        } else {
            ListSource::File
        };
        let Ok(list) = readlist::read_list(list_def.as_bytes(), source) else {
            bail!("Error: failed to read {list_def}");
        };
        let mut bins = Vec::with_capacity(list.len() + 2);
        for item in &list {
            let (v, end) = strtod(item, 0);
            if end != item.len() {
                bail!(
                    "Could not parse {list_def}: {}",
                    String::from_utf8_lossy(item)
                );
            }
            let v = v as f32;
            if v < min || v > max {
                bail!(
                    "Expected values from the interval [{min:.6},{max:.6}], found {}",
                    String::from_utf8_lossy(item)
                );
            }
            bins.push(v);
        }
        if bins.len() < 2 {
            // bcftools fails an assertion here.
            bail!("Error: --af-bins needs at least two values: {list_def}");
        }
        let max_err = (bins[1] - bins[0]) * 1e-6;
        if (bins[0] - min).abs() > max_err {
            bins.insert(0, min);
        }
        if (bins[bins.len() - 1] - max).abs() > max_err {
            bins.push(max);
        }
        Ok(Self { bins })
    }

    /// Returns the number of boundaries (`bin_get_size`).
    pub fn len(&self) -> usize {
        self.bins.len()
    }

    /// Returns boundary `i` (`bin_get_value`).
    pub fn value(&self, i: usize) -> f32 {
        self.bins[i]
    }

    /// Returns the bin of a value, with half-open intervals `[a, b)` (`bin_get_idx`).
    pub fn index(&self, value: f32) -> usize {
        let n = self.bins.len();
        if self.bins[n - 1] < value {
            return n - 1;
        }
        let mut imin: isize = 0;
        let mut imax: isize = n as isize - 2;
        while imin < imax {
            let i = (imin + imax) / 2;
            let b = self.bins[i as usize];
            if value < b {
                imax = i - 1;
            } else if value > b {
                imin = i + 1;
            } else {
                return i as usize;
            }
        }
        if imax >= 0 && self.bins[imax as usize] <= value {
            return imax as usize;
        }
        if imin > 0 { (imin - 1) as usize } else { 0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dist_is_exact_then_log_scaled() {
        let mut d = Dist::new(5);
        d.insert(0);
        d.insert(100_000);
        d.insert(100_001);
        d.insert(1_234_567);
        assert_eq!(d.get(0), (1, 0));
        assert_eq!(d.get(100_000), (2, 100_000));
        // 1,234,567 is at level 2 (step 100) in the bin starting at 1,234,500.
        let ibin = 100_000 + 90_000 + (1_234_567 - 1_000_000) / 100;
        assert_eq!(d.len(), ibin + 1);
        assert_eq!(d.get(ibin), (1, 1_234_500));
    }

    #[test]
    fn idist_bins() {
        let mut d = IDist::try_new(DepthBins {
            min: 10,
            max: 20,
            step: 5,
        })
        .expect("allocate six bins");
        for v in [5, 10, 14, 15, 20, 21] {
            d.add(v);
        }
        assert_eq!(d.vals, [1, 2, 1, 1, 0, 1]);
        assert_eq!(d.bin_start(2), 15);
    }

    #[test]
    fn af_bins_search() {
        let b = AfBins::parse("0.1,0.5").expect("the bins parse");
        assert_eq!(b.len(), 4);
        assert_eq!(b.index(0.0), 0);
        assert_eq!(b.index(0.05), 0);
        assert_eq!(b.index(0.1), 1);
        assert_eq!(b.index(0.3), 1);
        assert_eq!(b.index(0.5), 2);
        assert_eq!(b.index(1.0), 2);
    }
}
