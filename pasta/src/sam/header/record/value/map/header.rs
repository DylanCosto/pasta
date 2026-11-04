//! SAM header record header map value.

mod builder;
pub mod group_order;
pub mod sort_order;
pub mod tag;
pub mod version;

pub use self::{tag::Tag, version::Version};

use bstr::BString;

use self::builder::Builder;
use super::{Inner, Map, OtherFields};

/// A SAM header record header map value.
///
/// The header describes file-level metadata. The format version is guaranteed to be set.
///
/// A header read by a lenient parser (see [`crate::sam::header::Parser::lenient`]) may have no
/// `VN` field, or one that is not `<major>.<minor>`, as htslib allows. Its [`Map::version`] is
/// then the default version, and the original field (or its absence) is kept so that writing the
/// header reproduces it, until the version is changed with [`Map::version_mut`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    pub(crate) version: Version,
    pub(crate) raw_version: Option<RawVersion>,
}

/// A `VN` field that a lenient parser could not read as a [`Version`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum RawVersion {
    /// The record has no `VN` field.
    Missing,
    /// The `VN` value is not `<major>.<minor>`.
    Invalid(BString),
}

impl Inner for Header {
    type StandardTag = tag::Standard;
    type Builder = Builder;
}

impl Map<Header> {
    /// Creates a SAM header record header map value with a format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::header::record::value::{
    ///     map::{self, header::Version},
    ///     Map,
    /// };
    ///
    /// let header = Map::<map::Header>::new(Version::new(1, 6));
    /// ```
    pub fn new(version: Version) -> Self {
        Self {
            inner: Header {
                version,
                raw_version: None,
            },
            other_fields: OtherFields::new(),
        }
    }

    /// Returns the format version.
    ///
    /// For a header read leniently without a valid `VN` field, this is [`Version::default`].
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::header::record::value::{
    ///     map::{self, header::Version},
    ///     Map,
    /// };
    ///
    /// let header = Map::<map::Header>::new(Version::new(1, 6));
    /// assert_eq!(header.version(), Version::new(1, 6));
    /// ```
    pub fn version(&self) -> Version {
        self.inner.version
    }

    /// Returns a mutable reference to the format version.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::header::record::value::{map::{self, header::Version}, Map};
    /// let mut header = Map::<map::Header>::default();
    /// *header.version_mut() = Version::new(1, 5);
    /// assert_eq!(header.version(), Version::new(1, 5));
    /// ```
    pub fn version_mut(&mut self) -> &mut Version {
        // The caller may set a version, which must then be written instead of the field as read.
        self.inner.raw_version = None;
        &mut self.inner.version
    }

    pub(crate) fn raw_version(&self) -> Option<&RawVersion> {
        self.inner.raw_version.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Map::<Header>::default();
        assert_eq!(header.version(), Version::default());
    }
}
