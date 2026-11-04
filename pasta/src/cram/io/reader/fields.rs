use bitflags::bitflags;

use crate::sam::alignment::record::data::field::Tag;

bitflags! {
    /// The record fields a CRAM reader decodes, as HTSlib's `CRAM_OPT_REQUIRED_FIELDS`.
    ///
    /// Decoding a subset skips the data series only the other fields use, and the blocks that
    /// hold only those data series are not decompressed. The skipped fields of decoded records
    /// are left empty or unset (a skipped sequence reads as missing). Fields that share a block
    /// with a required field are decoded too, so a subset is a lower bound.
    ///
    /// [`Self::DATA`] decodes the data fields stored as tags, all of them or only those
    /// [`Tags`] selects. The read group is a separate data series, decoded with
    /// [`Self::READ_GROUP`].
    ///
    /// Reference sequences are read only for [`Self::SEQUENCE`], and missing `MD` and `NM` tags
    /// are added only with both [`Self::SEQUENCE`] and [`Self::DATA`] (and only those of the two
    /// that [`Tags`] selects).
    ///
    /// The default is every field.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::reader::Fields;
    ///
    /// let fields = Fields::FLAGS | Fields::MAPPING_QUALITY;
    /// assert!(!fields.contains(Fields::SEQUENCE));
    /// assert_eq!(Fields::default(), Fields::all());
    /// ```
    #[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
    pub struct Fields: u16 {
        /// The read name (`QNAME`).
        const NAME = 1 << 0;
        /// The flags (`FLAG`), including the mate flags of attached mates.
        const FLAGS = 1 << 1;
        /// The reference sequence ID (`RNAME`).
        const REFERENCE_SEQUENCE_ID = 1 << 2;
        /// The alignment start (`POS`).
        const ALIGNMENT_START = 1 << 3;
        /// The mapping quality (`MAPQ`).
        const MAPPING_QUALITY = 1 << 4;
        /// The CIGAR operations (`CIGAR`), and with them the alignment end.
        const CIGAR = 1 << 5;
        /// The mate reference sequence ID (`RNEXT`).
        const MATE_REFERENCE_SEQUENCE_ID = 1 << 6;
        /// The mate alignment start (`PNEXT`).
        const MATE_ALIGNMENT_START = 1 << 7;
        /// The template length (`TLEN`).
        const TEMPLATE_LENGTH = 1 << 8;
        /// The sequence (`SEQ`).
        const SEQUENCE = 1 << 9;
        /// The quality scores (`QUAL`).
        const QUALITY_SCORES = 1 << 10;
        /// The data fields stored as tags; see [`Tags`].
        const DATA = 1 << 11;
        /// The read group (`RG`), which CRAM stores as a data series rather than a tag.
        const READ_GROUP = 1 << 12;
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self::all()
    }
}

/// The data fields (tags) a CRAM reader decodes with [`Fields::DATA`].
///
/// Each tag is usually stored in an external block of its own, which is then not decompressed
/// when the tag is left out (a tag that shares a block with a decoded field is read but still
/// left out). A decoded record holds the selected tags in the order a full decode gives them,
/// regenerated `MD` and `NM` tags included, as if the others were removed.
///
/// The read group is not a tag in CRAM; [`Fields::READ_GROUP`] decodes it. Only an `RG` tag that
/// a writer stored as a tag (because the header has no such read group) is selected here.
///
/// The default is every tag.
///
/// # Examples
///
/// ```
/// use pasta::cram::io::reader::Tags;
/// use pasta::sam::alignment::record::data::field::Tag;
///
/// let tags = Tags::Except(vec![Tag::BASE_MODIFICATIONS, Tag::BASE_MODIFICATION_PROBABILITIES]);
/// assert!(tags.contains(Tag::EDIT_DISTANCE));
/// assert!(!tags.contains(Tag::BASE_MODIFICATIONS));
/// ```
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub enum Tags {
    /// Every tag.
    #[default]
    All,
    /// Only these tags.
    Only(Vec<Tag>),
    /// Every tag but these.
    Except(Vec<Tag>),
}

impl Tags {
    /// Returns whether `tag` is selected.
    pub fn contains(&self, tag: Tag) -> bool {
        match self {
            Self::All => true,
            Self::Only(tags) => tags.contains(&tag),
            Self::Except(tags) => !tags.contains(&tag),
        }
    }

    /// Returns whether every tag is selected.
    pub fn is_all(&self) -> bool {
        match self {
            Self::All => true,
            Self::Only(_) => false,
            Self::Except(tags) => tags.is_empty(),
        }
    }
}
