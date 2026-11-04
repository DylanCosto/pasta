//! The data series, tags, and blocks a slice decodes for a set of required fields and tags.
//!
//! This follows HTSlib's `cram_dependent_data_series`: the fields name data series, and any data
//! series that shares a block with a decoded one (or uses the core block when a decoded one does)
//! is decoded too, since the reads of a shared block must stay in order. Tags are handled the
//! same way, one by one: a tag that is not required but shares a block with a decoded data series
//! or tag is decoded and then left out of the record.

use crate::cram::{
    container::{
        CompressionHeader, block,
        compression_header::{
            Encoding,
            data_series_encodings::DataSeries,
            encoding::codec::{Byte, ByteArray, Integer},
        },
    },
    io::reader::{Fields, MdNm, Tags},
};
use crate::sam::alignment::record::data::field::Tag;

/// A set of data series, and whether tags are decoded.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct DataSeriesSet(u32);

/// The bit for tags, after those of the data series.
const TAGS: u32 = 1 << 31;

impl DataSeriesSet {
    pub(crate) const ALL: Self = Self(u32::MAX);

    const fn of(data_series: &[DataSeries]) -> Self {
        let mut bits = 0;
        let mut i = 0;

        while i < data_series.len() {
            bits |= bit(data_series[i]);
            i += 1;
        }

        Self(bits)
    }

    pub(crate) fn contains(self, data_series: DataSeries) -> bool {
        self.0 & bit(data_series) != 0
    }

    pub(crate) fn contains_tags(self) -> bool {
        self.0 & TAGS != 0
    }

    fn intersects(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }

    fn insert(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

const fn bit(data_series: DataSeries) -> u32 {
    1 << data_series as u32
}

/// The data series that set CIGAR operations: the features that consume read or reference
/// bases, and the read length. Substitutions and stretches of bases or scores only change bases
/// or quality scores.
const CIGAR: DataSeriesSet = DataSeriesSet::of(&[
    DataSeries::ReadLengths,
    DataSeries::FeatureCounts,
    DataSeries::FeatureCodes,
    DataSeries::FeaturePositionDeltas,
    DataSeries::DeletionLengths,
    DataSeries::InsertionBases,
    DataSeries::SoftClipBases,
    DataSeries::HardClipLengths,
    DataSeries::PaddingLengths,
    DataSeries::ReferenceSkipLengths,
]);

/// The data series read per feature.
const FEATURE_VALUES: DataSeriesSet = DataSeriesSet::of(&[
    DataSeries::FeatureCodes,
    DataSeries::FeaturePositionDeltas,
    DataSeries::DeletionLengths,
    DataSeries::StretchesOfBases,
    DataSeries::StretchesOfQualityScores,
    DataSeries::BaseSubstitutionCodes,
    DataSeries::InsertionBases,
    DataSeries::ReferenceSkipLengths,
    DataSeries::PaddingLengths,
    DataSeries::HardClipLengths,
    DataSeries::SoftClipBases,
    // Read-base and quality-score features read these too, so reading them for any record keeps
    // their blocks in order only if every feature is read.
    DataSeries::Bases,
    DataSeries::QualityScores,
]);

const FEATURES: DataSeriesSet = DataSeriesSet::of(&[
    DataSeries::FeatureCounts,
    DataSeries::FeatureCodes,
    DataSeries::FeaturePositionDeltas,
]);

/// What a slice decodes for a set of required fields and tags.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Selection {
    fields: Fields,
    series: DataSeriesSet,
    /// The external blocks the decoded data series and tags read, or `None` for all of them.
    blocks: Option<Vec<block::ContentId>>,
    /// The tags (by their keys' block content IDs) that are decoded, or `None` for all of them.
    decoded_tags: Option<Vec<block::ContentId>>,
    /// The tags decoded records keep.
    kept_tags: Tags,
    /// Whether the CRAM flags tag is read, for `MD` and `NM` regeneration.
    reads_cram_flags: bool,
    md_nm: MdNmTags,
}

/// What happens to a tag's values.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TagAction {
    /// Decoded and kept.
    Keep,
    /// Read, to keep a shared block in order, and left out.
    Discard,
    /// Not read.
    Skip,
}

/// Which of the `MD` and `NM` tags a record lacks are regenerated.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct MdNmTags {
    pub(crate) md: bool,
    pub(crate) nm: bool,
}

impl MdNmTags {
    #[cfg(test)]
    pub(crate) const ALL: Self = Self { md: true, nm: true };

    pub(crate) fn any(self) -> bool {
        self.md || self.nm
    }
}

/// The CRAM flags tag, which says whether `MD` and `NM` are regenerated.
const CRAM_FLAGS: Tag = Tag::new(b'c', b'F');

impl Selection {
    /// Every field and tag, with missing `MD` and `NM` tags left out.
    pub(crate) fn all() -> Self {
        Self {
            fields: Fields::all(),
            series: DataSeriesSet::ALL,
            blocks: None,
            decoded_tags: None,
            kept_tags: Tags::All,
            reads_cram_flags: true,
            md_nm: MdNmTags::default(),
        }
    }

    pub(crate) fn new(
        compression_header: &CompressionHeader,
        fields: Fields,
        tags: &Tags,
        md_nm: MdNm,
    ) -> Self {
        let md_nm = md_nm_tags(fields, tags, md_nm);

        if fields.is_all() && tags.is_all() {
            return Self {
                md_nm,
                ..Self::all()
            };
        }

        let kept_tags = if fields.contains(Fields::DATA) {
            tags.clone()
        } else {
            Tags::Only(Vec::new())
        };

        // The CRAM flags tag is always left out of records, but it is read when MD or NM may be
        // regenerated. This does not depend on `md_nm`, so that the blocks a slice decodes do not.
        let reads_cram_flags = md_nm_tags(fields, tags, MdNm::Regenerate).any();
        let is_read = |tag: Tag| kept_tags.contains(tag) || (reads_cram_flags && tag == CRAM_FLAGS);

        let mut decoded_tags: Vec<block::ContentId> = compression_header
            .tag_encodings()
            .keys()
            .copied()
            .filter(|&id| is_read(key_tag(id)))
            .collect();

        let mut series = required_data_series(fields);

        let blocks = loop {
            let prev = (series, decoded_tags.len());

            if !decoded_tags.is_empty() {
                series.insert(DataSeriesSet(TAGS | bit(DataSeries::TagSetIds)));
            }

            if series.intersects(FEATURE_VALUES) {
                series.insert(FEATURES);
            }

            if series.contains(DataSeries::Bases) || series.contains(DataSeries::QualityScores) {
                series.insert(DataSeriesSet::of(&[DataSeries::ReadLengths]));
            }

            let mut uses = Uses::default();

            for_each_encoding(compression_header, |item, encoding_uses| {
                if is_selected(series, &decoded_tags, item) {
                    uses.add(&encoding_uses);
                }
            });

            for_each_encoding(compression_header, |item, encoding_uses| {
                if !is_selected(series, &decoded_tags, item) && uses.overlaps(&encoding_uses) {
                    match item {
                        Item::DataSeries(data_series) => {
                            series.insert(DataSeriesSet::of(&[data_series]));
                        }
                        Item::Tag(id) => decoded_tags.push(id),
                    }
                }
            });

            if (series, decoded_tags.len()) == prev {
                break uses.blocks;
            }
        };

        Self {
            fields,
            series,
            blocks: Some(blocks),
            decoded_tags: Some(decoded_tags),
            kept_tags,
            reads_cram_flags,
            md_nm,
        }
    }

    pub(crate) fn data_series(&self) -> DataSeriesSet {
        self.series
    }

    /// Returns what happens to the values of the tag whose key has this block content ID.
    pub(crate) fn tag_action(&self, id: block::ContentId) -> TagAction {
        let tag = key_tag(id);

        if self
            .decoded_tags
            .as_ref()
            .is_some_and(|decoded_tags| !decoded_tags.contains(&id))
        {
            TagAction::Skip
        } else if self.kept_tags.contains(tag) || (self.reads_cram_flags && tag == CRAM_FLAGS) {
            // The CRAM flags tag is kept for MD and NM regeneration, which removes it.
            TagAction::Keep
        } else {
            TagAction::Discard
        }
    }

    /// Returns whether the records read the external block with this content ID.
    pub(crate) fn reads_block(&self, block_content_id: block::ContentId) -> bool {
        self.blocks
            .as_ref()
            .is_none_or(|blocks| blocks.contains(&block_content_id))
    }

    /// Returns whether the records' bases are rebuilt from reference sequences.
    pub(crate) fn reads_reference_sequences(&self) -> bool {
        self.fields.contains(Fields::SEQUENCE)
    }

    /// Returns which missing `MD` and `NM` tags are regenerated.
    pub(crate) fn md_nm(&self) -> MdNmTags {
        self.md_nm
    }
}

/// Missing `MD` and `NM` tags are regenerated only when both the sequence and the data are
/// decoded, and only those of the two that are selected.
fn md_nm_tags(fields: Fields, tags: &Tags, md_nm: MdNm) -> MdNmTags {
    if md_nm == MdNm::Regenerate && fields.contains(Fields::SEQUENCE | Fields::DATA) {
        MdNmTags {
            md: tags.contains(Tag::MISMATCHED_POSITIONS),
            nm: tags.contains(Tag::EDIT_DISTANCE),
        }
    } else {
        MdNmTags::default()
    }
}

/// The tag of a tag key's block content ID (its two tag bytes and its type byte).
fn key_tag(id: block::ContentId) -> Tag {
    let [_, l, r, _] = id.to_be_bytes();
    Tag::new(l, r)
}

/// A data series or a tag (by its key's block content ID).
#[derive(Clone, Copy, Debug)]
enum Item {
    DataSeries(DataSeries),
    Tag(block::ContentId),
}

fn is_selected(series: DataSeriesSet, decoded_tags: &[block::ContentId], item: Item) -> bool {
    match item {
        Item::DataSeries(data_series) => series.contains(data_series),
        Item::Tag(id) => decoded_tags.contains(&id),
    }
}

/// The data series that `fields` name, as HTSlib sets them.
fn required_data_series(fields: Fields) -> DataSeriesSet {
    use DataSeries::{
        AlignmentStarts, BamFlags, BaseSubstitutionCodes, Bases, CramFlags, MappingQualities,
        MateAlignmentStarts, MateDistances, MateFlags, MateReferenceSequenceIds, Names,
        QualityScores, ReadGroupIds, ReferenceSequenceIds, StretchesOfBases,
        StretchesOfQualityScores, TemplateLengths,
    };

    // Every record reads its flags; CRAM flags say which other data series it has.
    let mut series = DataSeriesSet::of(&[BamFlags, CramFlags]);

    let required: [(Fields, DataSeriesSet); 12] = [
        // Without stored names, mates take the name of the first record of their template.
        (Fields::NAME, DataSeriesSet::of(&[Names, MateDistances])),
        // Mate flags come from detached records' mate flags and attached records' mates.
        (
            Fields::FLAGS,
            DataSeriesSet::of(&[MateFlags, MateDistances]),
        ),
        (
            Fields::REFERENCE_SEQUENCE_ID,
            DataSeriesSet::of(&[ReferenceSequenceIds]),
        ),
        (
            Fields::ALIGNMENT_START,
            DataSeriesSet::of(&[AlignmentStarts]),
        ),
        (
            Fields::MAPPING_QUALITY,
            DataSeriesSet::of(&[MappingQualities]),
        ),
        (Fields::CIGAR, CIGAR),
        (
            Fields::MATE_REFERENCE_SEQUENCE_ID,
            DataSeriesSet::of(&[
                MateFlags,
                MateDistances,
                MateReferenceSequenceIds,
                ReferenceSequenceIds,
            ]),
        ),
        (
            Fields::MATE_ALIGNMENT_START,
            DataSeriesSet::of(&[
                MateFlags,
                MateDistances,
                MateAlignmentStarts,
                AlignmentStarts,
            ]),
        ),
        (
            Fields::TEMPLATE_LENGTH,
            DataSeriesSet(
                CIGAR.0
                    | DataSeriesSet::of(&[
                        MateFlags,
                        MateDistances,
                        TemplateLengths,
                        AlignmentStarts,
                        ReferenceSequenceIds,
                    ])
                    .0,
            ),
        ),
        (
            Fields::SEQUENCE,
            DataSeriesSet(
                CIGAR.0
                    | DataSeriesSet::of(&[
                        Bases,
                        BaseSubstitutionCodes,
                        StretchesOfBases,
                        AlignmentStarts,
                        ReferenceSequenceIds,
                    ])
                    .0,
            ),
        ),
        (
            Fields::QUALITY_SCORES,
            DataSeriesSet(
                CIGAR.0 | DataSeriesSet::of(&[QualityScores, StretchesOfQualityScores]).0,
            ),
        ),
        // Tags are selected one by one (`Selection::new`).
        (Fields::READ_GROUP, DataSeriesSet::of(&[ReadGroupIds])),
    ];

    for (field, data_series) in required {
        if fields.contains(field) {
            series.insert(data_series);
        }
    }

    series
}

/// The blocks an encoding reads.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct Uses {
    core: bool,
    blocks: Vec<block::ContentId>,
}

impl Uses {
    fn add(&mut self, other: &Self) {
        self.core |= other.core;

        for &id in &other.blocks {
            if !self.blocks.contains(&id) {
                self.blocks.push(id);
            }
        }
    }

    fn overlaps(&self, other: &Self) -> bool {
        (self.core && other.core) || other.blocks.iter().any(|id| self.blocks.contains(id))
    }

    fn of_integer(encoding: &Encoding<Integer>) -> Self {
        match encoding.get() {
            Integer::External { block_content_id } => Self::external(*block_content_id),
            Integer::Huffman { alphabet, .. } => Self::huffman(alphabet),
            Integer::Golomb { .. }
            | Integer::Beta { .. }
            | Integer::Subexp { .. }
            | Integer::GolombRice { .. }
            | Integer::Gamma { .. } => Self::core(),
        }
    }

    fn of_byte(encoding: &Encoding<Byte>) -> Self {
        match encoding.get() {
            Byte::External { block_content_id } => Self::external(*block_content_id),
            Byte::Huffman { alphabet, .. } => Self::huffman(alphabet),
            Byte::Beta { .. } => Self::core(),
        }
    }

    fn of_byte_array(encoding: &Encoding<ByteArray>) -> Self {
        match encoding.get() {
            ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => {
                let mut uses = Self::of_integer(len_encoding);
                uses.add(&Self::of_byte(value_encoding));
                uses
            }
            ByteArray::ByteArrayStop {
                block_content_id, ..
            } => Self::external(*block_content_id),
        }
    }

    fn external(block_content_id: block::ContentId) -> Self {
        Self {
            core: false,
            blocks: vec![block_content_id],
        }
    }

    fn core() -> Self {
        Self {
            core: true,
            blocks: Vec::new(),
        }
    }

    /// A single-symbol Huffman code reads no bits.
    fn huffman(alphabet: &[i32]) -> Self {
        if alphabet.len() == 1 {
            Self::default()
        } else {
            Self::core()
        }
    }
}

/// Calls `f` with each data series and tag and the blocks its encoding reads.
fn for_each_encoding<F>(compression_header: &CompressionHeader, mut f: F)
where
    F: FnMut(Item, Uses),
{
    use DataSeries as Ds;

    let encodings = compression_header.data_series_encodings();

    let integers = [
        (Ds::BamFlags, encodings.bam_flags()),
        (Ds::CramFlags, encodings.cram_flags()),
        (Ds::ReferenceSequenceIds, encodings.reference_sequence_ids()),
        (Ds::ReadLengths, encodings.read_lengths()),
        (Ds::AlignmentStarts, encodings.alignment_starts()),
        (Ds::ReadGroupIds, encodings.read_group_ids()),
        (Ds::MateFlags, encodings.mate_flags()),
        (
            Ds::MateReferenceSequenceIds,
            encodings.mate_reference_sequence_ids(),
        ),
        (Ds::MateAlignmentStarts, encodings.mate_alignment_starts()),
        (Ds::TemplateLengths, encodings.template_lengths()),
        (Ds::MateDistances, encodings.mate_distances()),
        (Ds::TagSetIds, encodings.tag_set_ids()),
        (Ds::FeatureCounts, encodings.feature_counts()),
        (
            Ds::FeaturePositionDeltas,
            encodings.feature_position_deltas(),
        ),
        (Ds::DeletionLengths, encodings.deletion_lengths()),
        (Ds::ReferenceSkipLengths, encodings.reference_skip_lengths()),
        (Ds::PaddingLengths, encodings.padding_lengths()),
        (Ds::HardClipLengths, encodings.hard_clip_lengths()),
        (Ds::MappingQualities, encodings.mapping_qualities()),
    ];

    for (data_series, encoding) in integers {
        if let Some(encoding) = encoding {
            f(Item::DataSeries(data_series), Uses::of_integer(encoding));
        }
    }

    let bytes = [
        (Ds::FeatureCodes, encodings.feature_codes()),
        (
            Ds::BaseSubstitutionCodes,
            encodings.base_substitution_codes(),
        ),
        (Ds::Bases, encodings.bases()),
        (Ds::QualityScores, encodings.quality_scores()),
    ];

    for (data_series, encoding) in bytes {
        if let Some(encoding) = encoding {
            f(Item::DataSeries(data_series), Uses::of_byte(encoding));
        }
    }

    let byte_arrays = [
        (Ds::Names, encodings.names()),
        (Ds::StretchesOfBases, encodings.stretches_of_bases()),
        (
            Ds::StretchesOfQualityScores,
            encodings.stretches_of_quality_scores(),
        ),
        (Ds::InsertionBases, encodings.insertion_bases()),
        (Ds::SoftClipBases, encodings.soft_clip_bases()),
    ];

    for (data_series, encoding) in byte_arrays {
        if let Some(encoding) = encoding {
            f(Item::DataSeries(data_series), Uses::of_byte_array(encoding));
        }
    }

    for (&id, encoding) in compression_header.tag_encodings() {
        f(Item::Tag(id), Uses::of_byte_array(encoding));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compression_header() -> CompressionHeader {
        use crate::cram::container::compression_header::{
            DataSeriesEncodings, PreservationMap, TagEncodings,
        };

        CompressionHeader::new(
            PreservationMap::default(),
            DataSeriesEncodings::init(),
            TagEncodings::default(),
        )
    }

    #[test]
    fn test_new_with_all_fields() {
        let selection = Selection::new(
            &compression_header(),
            Fields::all(),
            &Tags::All,
            MdNm::Regenerate,
        );
        assert_eq!(selection.data_series(), DataSeriesSet::ALL);
        assert!(selection.reads_block(42));
        assert!(selection.reads_reference_sequences());
        assert_eq!(selection.md_nm(), MdNmTags::ALL);
    }

    #[test]
    fn test_new_with_flags() {
        let selection = Selection::new(
            &compression_header(),
            Fields::FLAGS,
            &Tags::All,
            MdNm::Regenerate,
        );
        let series = selection.data_series();

        for data_series in [
            DataSeries::BamFlags,
            DataSeries::CramFlags,
            DataSeries::MateFlags,
            DataSeries::MateDistances,
        ] {
            assert!(series.contains(data_series), "{data_series:?}");
            assert!(selection.reads_block(block::ContentId::from(data_series)));
        }

        for data_series in [
            DataSeries::Names,
            DataSeries::FeatureCounts,
            DataSeries::Bases,
            DataSeries::QualityScores,
            DataSeries::TagSetIds,
            DataSeries::ReadGroupIds,
        ] {
            assert!(!series.contains(data_series), "{data_series:?}");
            assert!(!selection.reads_block(block::ContentId::from(data_series)));
        }

        assert!(!series.contains_tags());
        assert!(!selection.reads_reference_sequences());
        assert_eq!(selection.md_nm(), MdNmTags::default());
    }

    #[test]
    fn test_new_with_cigar() {
        let selection = Selection::new(
            &compression_header(),
            Fields::CIGAR,
            &Tags::All,
            MdNm::Regenerate,
        );
        let series = selection.data_series();

        assert!(series.contains(DataSeries::FeatureCodes));
        assert!(series.contains(DataSeries::SoftClipBases));
        assert!(!series.contains(DataSeries::Bases));
        assert!(!series.contains(DataSeries::QualityScores));
        assert!(!series.contains(DataSeries::BaseSubstitutionCodes));
    }

    #[test]
    fn test_new_with_a_shared_block() {
        use crate::cram::container::compression_header::{
            DataSeriesEncodings, PreservationMap, TagEncodings,
        };

        // Mapping qualities share the quality scores' block, so reading either reads both, and
        // quality scores are read in features too.
        let mut encodings = DataSeriesEncodings::init();
        encodings.mapping_qualities = Some(Encoding::new(Integer::External {
            block_content_id: block::ContentId::from(DataSeries::QualityScores),
        }));

        let compression_header = CompressionHeader::new(
            PreservationMap::default(),
            encodings,
            TagEncodings::default(),
        );

        let selection = Selection::new(
            &compression_header,
            Fields::MAPPING_QUALITY,
            &Tags::All,
            MdNm::Regenerate,
        );
        let series = selection.data_series();

        assert!(series.contains(DataSeries::QualityScores));
        assert!(series.contains(DataSeries::FeatureCodes));
        assert!(series.contains(DataSeries::ReadLengths));
        assert!(!series.contains(DataSeries::Names));
    }

    #[test]
    fn test_new_with_the_core_block() {
        use crate::cram::container::compression_header::{
            DataSeriesEncodings, PreservationMap, TagEncodings,
        };

        // Every data series in the core block is read when one is.
        let mut encodings = DataSeriesEncodings::init();
        encodings.bam_flags = Some(Encoding::new(Integer::Beta { offset: 0, len: 12 }));
        encodings.names = Some(Encoding::new(ByteArray::ByteArrayLength {
            len_encoding: Encoding::new(Integer::Gamma { offset: 0 }),
            value_encoding: Encoding::new(Byte::External {
                block_content_id: block::ContentId::from(DataSeries::Names),
            }),
        }));

        let compression_header = CompressionHeader::new(
            PreservationMap::default(),
            encodings,
            TagEncodings::default(),
        );

        let selection = Selection::new(
            &compression_header,
            Fields::FLAGS,
            &Tags::All,
            MdNm::Regenerate,
        );
        assert!(selection.data_series().contains(DataSeries::Names));
        assert!(selection.reads_block(block::ContentId::from(DataSeries::Names)));
        assert!(!selection.reads_block(block::ContentId::from(DataSeries::Bases)));
    }

    /// A compression header whose tags `MM`, `ML`, and `HP` each have a block of their own, and
    /// whose `XS` shares `HP`'s block.
    fn compression_header_with_tags() -> (CompressionHeader, [block::ContentId; 4]) {
        use crate::cram::container::compression_header::{
            DataSeriesEncodings, PreservationMap, TagEncodings, preservation_map::tag_sets::Key,
        };
        use crate::sam::alignment::record::data::field::Type;

        let keys = [
            Key::new(Tag::BASE_MODIFICATIONS, Type::String),
            Key::new(Tag::BASE_MODIFICATION_PROBABILITIES, Type::Array),
            Key::new(Tag::new(b'H', b'P'), Type::UInt8),
            Key::new(Tag::new(b'X', b'S'), Type::UInt8),
        ];
        let ids = keys.map(block::ContentId::from);

        let mut tag_encodings = TagEncodings::default();

        for (i, &id) in ids.iter().enumerate() {
            let block_content_id = if i == 3 { ids[2] } else { id };

            tag_encodings.insert(
                id,
                Encoding::new(ByteArray::ByteArrayStop {
                    stop_byte: 0,
                    block_content_id,
                }),
            );
        }

        let compression_header = CompressionHeader::new(
            PreservationMap::default(),
            DataSeriesEncodings::init(),
            tag_encodings,
        );

        (compression_header, ids)
    }

    #[test]
    fn test_new_with_only_some_tags() {
        let (compression_header, [mm, ml, hp, xs]) = compression_header_with_tags();
        let tags = Tags::Only(vec![Tag::new(b'H', b'P')]);
        let selection = Selection::new(&compression_header, Fields::DATA, &tags, MdNm::Regenerate);

        assert!(selection.data_series().contains(DataSeries::TagSetIds));
        assert!(!selection.data_series().contains(DataSeries::ReadGroupIds));
        assert_eq!(selection.tag_action(hp), TagAction::Keep);
        // XS shares HP's block, so it is read but left out.
        assert_eq!(selection.tag_action(xs), TagAction::Discard);
        assert_eq!(selection.tag_action(mm), TagAction::Skip);
        assert_eq!(selection.tag_action(ml), TagAction::Skip);
        assert!(selection.reads_block(hp));
        assert!(!selection.reads_block(mm));
        assert!(!selection.reads_block(ml));
    }

    #[test]
    fn test_new_with_all_fields_except_some_tags() {
        let (compression_header, [mm, ml, hp, xs]) = compression_header_with_tags();
        let tags = Tags::Except(vec![
            Tag::BASE_MODIFICATIONS,
            Tag::BASE_MODIFICATION_PROBABILITIES,
        ]);
        let selection = Selection::new(&compression_header, Fields::all(), &tags, MdNm::Regenerate);

        assert_eq!(selection.tag_action(hp), TagAction::Keep);
        assert_eq!(selection.tag_action(xs), TagAction::Keep);
        assert_eq!(selection.tag_action(mm), TagAction::Skip);
        assert_eq!(selection.tag_action(ml), TagAction::Skip);
        assert!(!selection.reads_block(mm));
        assert!(selection.reads_block(block::ContentId::from(DataSeries::QualityScores)));
        assert_eq!(selection.md_nm(), MdNmTags::ALL);
    }

    #[test]
    fn test_new_with_the_read_group() {
        let (compression_header, ids) = compression_header_with_tags();
        let selection = Selection::new(
            &compression_header,
            Fields::READ_GROUP,
            &Tags::All,
            MdNm::Regenerate,
        );

        assert!(selection.data_series().contains(DataSeries::ReadGroupIds));
        assert!(!selection.data_series().contains(DataSeries::TagSetIds));
        assert!(!selection.data_series().contains_tags());

        for id in ids {
            assert_eq!(selection.tag_action(id), TagAction::Skip);
            assert!(!selection.reads_block(id));
        }
    }

    #[test]
    fn test_md_nm_tags() {
        let fields = Fields::SEQUENCE | Fields::DATA;
        let nm = Tags::Only(vec![Tag::EDIT_DISTANCE]);

        assert_eq!(
            md_nm_tags(fields, &Tags::All, MdNm::Regenerate),
            MdNmTags::ALL
        );
        assert_eq!(
            md_nm_tags(fields, &nm, MdNm::Regenerate),
            MdNmTags {
                md: false,
                nm: true
            }
        );
        assert_eq!(md_nm_tags(fields, &nm, MdNm::AsStored), MdNmTags::default());
        assert_eq!(
            md_nm_tags(Fields::DATA, &Tags::All, MdNm::Regenerate),
            MdNmTags::default()
        );
    }

    #[test]
    fn test_key_tag() {
        use crate::cram::container::compression_header::preservation_map::tag_sets::Key;
        use crate::sam::alignment::record::data::field::Type;

        let key = Key::new(Tag::BASE_MODIFICATION_PROBABILITIES, Type::Array);
        assert_eq!(
            key_tag(block::ContentId::from(key)),
            Tag::BASE_MODIFICATION_PROBABILITIES
        );
    }
}
