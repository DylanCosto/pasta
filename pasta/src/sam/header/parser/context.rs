use crate::sam::header::record::value::map::header::Version;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Context {
    allow_duplicate_tags: bool,
    lenient: bool,
}

impl Context {
    /// Returns a context that accepts what htslib 1.24 accepts (see [`super::Parser::lenient`]).
    pub fn lenient() -> Self {
        Self {
            allow_duplicate_tags: true,
            lenient: true,
        }
    }

    pub fn allow_duplicate_tags(&self) -> bool {
        self.allow_duplicate_tags
    }

    /// Whether the header is parsed as htslib parses it rather than as the SAM specification
    /// requires.
    pub fn is_lenient(&self) -> bool {
        self.lenient
    }
}

impl Default for Context {
    fn default() -> Self {
        Self::from(Version::default())
    }
}

impl From<Version> for Context {
    fn from(version: Version) -> Self {
        Self {
            // SAM 1.6 no longer allows duplicate tags. See § 1.3 "The header section" (2025-08-12):
            // "Within each (non-`@CO`) header line, no field tag may appear more than once..."
            allow_duplicate_tags: version < Version::new(1, 6),
            lenient: false,
        }
    }
}
