use std::{error, fmt};

const MIN: u8 = 0;
const MAX: u8 = 12;

/// A DEFLATE compression level, 0 to 12, on libdeflate's scale, which zflate keeps: HTSlib maps
/// zlib's levels 1–9 onto its 1, 2, 3, 5, 6, 7, 8, 10 and 12.
///
/// The default is 5, zflate's fast level, chosen for pipelines: it writes BAM 1.3–1.7 times as
/// fast as libdeflate's 6, under 1% larger. VCF and BCF compress worse at 5 (11–12% larger than
/// HTSlib's default, libdeflate's 7), so variant output uses [`Self::VARIANTS`].
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct CompressionLevel(u8);

impl CompressionLevel {
    /// No compression.
    pub const NONE: Self = Self(0);

    /// A compression level optimized for speed.
    pub const FAST: Self = Self(1);

    /// A compression level optimized for compression rate.
    pub const BEST: Self = Self(MAX);

    /// The default for VCF and BCF output, 6: HTSlib's zlib-style level 5. It writes them 1.3–1.4
    /// times as fast as HTSlib's default (libdeflate's 7), with files 5% larger.
    pub const VARIANTS: Self = Self(6);

    /// Creates a compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::bgzf::io::writer::CompressionLevel;
    /// assert_eq!(CompressionLevel::new(0), Some(CompressionLevel::NONE));
    /// assert!(CompressionLevel::new(255).is_none());
    /// ```
    pub const fn new(n: u8) -> Option<Self> {
        if n <= MAX { Some(Self(n)) } else { None }
    }

    /// Returns the inner value.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::bgzf::io::writer::CompressionLevel;
    /// assert_eq!(CompressionLevel::NONE.get(), 0);
    /// ```
    pub const fn get(&self) -> u8 {
        self.0
    }
}

impl Default for CompressionLevel {
    fn default() -> Self {
        Self(5)
    }
}

/// An error returned when a raw DEFLATE compression level fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromU8Error {
    Invalid(u8),
}

impl error::Error for TryFromU8Error {}

impl fmt::Display for TryFromU8Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(n) => write!(f, "invalid input: {n}"),
        }
    }
}

impl TryFrom<u8> for CompressionLevel {
    type Error = TryFromU8Error;

    fn try_from(n: u8) -> Result<Self, Self::Error> {
        match n {
            MIN..=MAX => Ok(Self(n)),
            _ => Err(TryFromU8Error::Invalid(n)),
        }
    }
}

impl From<CompressionLevel> for u8 {
    fn from(compression_level: CompressionLevel) -> Self {
        compression_level.0
    }
}

impl CompressionLevel {
    /// The same level, as zflate takes it.
    pub(crate) fn to_zflate(self) -> zflate::CompressionLevel {
        zflate::CompressionLevel::new(self.0).expect("zflate takes levels 0 to 12")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(CompressionLevel::default(), CompressionLevel(5));
        assert_eq!(CompressionLevel::BEST, CompressionLevel(12));
    }
}
