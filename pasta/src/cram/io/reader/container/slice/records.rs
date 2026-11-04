mod data;
mod external_data_readers;
mod selection;

pub use external_data_readers::ExternalDataReaders;
use selection::TagAction;
pub(crate) use selection::{DataSeriesSet, MdNmTags, Selection};

use std::{borrow::Cow, error, fmt, io};

use crate::core::Position;
use crate::sam;

use crate::cram::{
    Record,
    container::{
        CompressionHeader, ReferenceSequenceContext, block,
        compression_header::{
            Encoding, data_series_encodings::DataSeries, encoding::codec::ByteArray,
            preservation_map::tag_sets,
        },
    },
    io::BitReader,
    record::{Feature, Flags, MateFlags, feature},
};

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReadRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag_sets::Key),
}

impl error::Error for ReadRecordError {}

impl fmt::Display for ReadRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDataSeriesEncoding(data_series) => {
                write!(f, "missing data series encoding: {data_series:?}")
            }
            Self::MissingTagEncoding(key) => write!(f, "missing tag encoding: {key:?}"),
        }
    }
}

pub struct Records<'c, 'ch: 'c> {
    compression_header: &'ch CompressionHeader,
    core_data_reader: BitReader<'c>,
    external_data_readers: ExternalDataReaders<'c>,
    reference_sequence_context: ReferenceSequenceContext,
    id: u64,
    prev_alignment_start: Option<Position>,
    /// Each tag set's tags that are read, and whether each is kept.
    resolved_tag_sets: Vec<Vec<ResolvedTag<'ch>>>,
    /// The data series read; the fields of the others are left unset.
    series: DataSeriesSet,
}

impl<'c, 'ch: 'c> Records<'c, 'ch> {
    pub fn new(
        compression_header: &'ch CompressionHeader,
        core_data_reader: BitReader<'c>,
        external_data_readers: ExternalDataReaders<'c>,
        reference_sequence_context: ReferenceSequenceContext,
        initial_id: u64,
        selection: &Selection,
    ) -> io::Result<Self> {
        let series = selection.data_series();

        let initial_alignment_start = match reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.alignment_start()),
            _ => None,
        };

        let resolved_tag_sets = if series.contains_tags() {
            resolve_tag_sets(compression_header)?
        } else {
            Vec::new()
        };

        // Skipped tags read nothing, so they are left out of the tag sets.
        let resolved_tag_sets = resolved_tag_sets
            .into_iter()
            .map(|tag_set| {
                tag_set
                    .into_iter()
                    .filter_map(|(key, encoding)| {
                        let is_kept = match selection.tag_action(block::ContentId::from(key)) {
                            TagAction::Keep => true,
                            TagAction::Discard => false,
                            TagAction::Skip => return None,
                        };

                        Some(ResolvedTag {
                            key,
                            decoder: TagDecoder::new(encoding, &external_data_readers),
                            is_kept,
                        })
                    })
                    .collect()
            })
            .collect();

        Ok(Self {
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
            id: initial_id,
            prev_alignment_start: initial_alignment_start,
            resolved_tag_sets,
            series,
        })
    }

    fn reads(&self, data_series: DataSeries) -> bool {
        self.series.contains(data_series)
    }

    pub fn read_record(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        record.id = self.id;

        record.bam_flags = self.read_bam_flags()?;
        record.cram_flags = self.read_cram_flags()?;

        self.read_positions(record)?;
        self.read_names(record)?;
        self.read_mate(record)?;
        self.read_data(record)?;

        if record.bam_flags.is_unmapped() {
            self.read_unmapped_read(record)?;
        } else {
            self.read_mapped_read(record)?;
        }

        self.id += 1;
        self.prev_alignment_start = record.alignment_start;

        Ok(())
    }

    fn read_bam_flags(&mut self) -> io::Result<sam::alignment::record::Flags> {
        self.compression_header
            .data_series_encodings()
            .bam_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::BamFlags))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::alignment::record::Flags::from)
    }

    fn read_cram_flags(&mut self) -> io::Result<Flags> {
        self.compression_header
            .data_series_encodings()
            .cram_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::CramFlags))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Flags::from)
    }

    fn read_positions(&mut self, record: &mut Record) -> io::Result<()> {
        record.reference_sequence_id = match self.reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.reference_sequence_id()),
            ReferenceSequenceContext::None => None,
            ReferenceSequenceContext::Many if self.reads(DataSeries::ReferenceSequenceIds) => {
                self.read_reference_sequence_id()?
            }
            ReferenceSequenceContext::Many => None,
        };

        if self.reads(DataSeries::ReadLengths) {
            record.read_length = self.read_read_length()?;
        }

        if self.reads(DataSeries::AlignmentStarts) {
            record.alignment_start = self.read_alignment_start()?;
        }

        if self.reads(DataSeries::ReadGroupIds) {
            record.read_group_id = self.read_read_group_id()?;
        }

        Ok(())
    }

    fn read_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .reference_sequence_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReferenceSequenceIds))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                UNMAPPED => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_read_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .read_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReadLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_alignment_start(&mut self) -> io::Result<Option<Position>> {
        let alignment_starts_are_deltas = self
            .compression_header
            .preservation_map()
            .alignment_starts_are_deltas();

        let alignment_start_or_delta = self
            .compression_header
            .data_series_encodings()
            .alignment_starts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::AlignmentStarts))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

        let alignment_start = if alignment_starts_are_deltas {
            let prev_alignment_start = i32::try_from(
                self.prev_alignment_start
                    .map(usize::from)
                    .unwrap_or_default(),
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            prev_alignment_start + alignment_start_or_delta
        } else {
            alignment_start_or_delta
        };

        usize::try_from(alignment_start)
            .map(Position::new)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    fn read_read_group_id(&mut self) -> io::Result<Option<usize>> {
        // § 10.2 "CRAM positional data" (2021-10-15): "-1 for no group".
        const MISSING: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .read_group_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReadGroupIds))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                MISSING => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_names(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        let preservation_map = self.compression_header.preservation_map();

        // Missing read names are generated when resolving mates.
        if preservation_map.records_have_names() && self.reads(DataSeries::Names) {
            record.name = self.read_name()?;
        }

        Ok(())
    }

    fn read_name(&mut self) -> io::Result<Option<Cow<'c, [u8]>>> {
        const MISSING: &[u8] = b"*";

        self.compression_header
            .data_series_encodings()
            .names()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Names))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .map(|buf| match &buf[..] {
                MISSING => None,
                _ => Some(buf),
            })
    }

    fn read_mate(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        if record.cram_flags.is_detached() {
            if self.reads(DataSeries::MateFlags) {
                record.mate_flags = self.read_mate_flags()?;
            }

            if record.mate_flags.is_on_negative_strand() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if record.mate_flags.is_unmapped() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
            }

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.records_have_names() && self.reads(DataSeries::Names) {
                record.name = self.read_name()?;
            }

            if self.reads(DataSeries::MateReferenceSequenceIds) {
                record.mate_reference_sequence_id = self.read_mate_reference_sequence_id()?;
            }

            if self.reads(DataSeries::MateAlignmentStarts) {
                record.mate_alignment_start = self.read_mate_alignment_start()?;
            }

            if self.reads(DataSeries::TemplateLengths) {
                record.template_length = self.read_template_length()?;
            }
        } else if record.cram_flags.mate_is_downstream() && self.reads(DataSeries::MateDistances) {
            record.mate_distance = self.read_mate_distance().map(Some)?;
        }

        Ok(())
    }

    fn read_mate_flags(&mut self) -> io::Result<MateFlags> {
        self.compression_header
            .data_series_encodings()
            .mate_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateFlags))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(MateFlags::from)
    }

    fn read_mate_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .mate_reference_sequence_ids()
            .ok_or_else(|| {
                missing_data_series_encoding_error(DataSeries::MateReferenceSequenceIds)
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| match id {
                UNMAPPED => Ok(None),
                _ => usize::try_from(id)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_mate_alignment_start(&mut self) -> io::Result<Option<Position>> {
        self.compression_header
            .data_series_encodings()
            .mate_alignment_starts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateAlignmentStarts))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Position::new)
    }

    fn read_template_length(&mut self) -> io::Result<i32> {
        self.compression_header
            .data_series_encodings()
            .template_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::TemplateLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_mate_distance(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .mate_distances()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateDistances))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_data(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        if !self.reads(DataSeries::TagSetIds) {
            return Ok(());
        }

        let tag_set_id = self.read_tag_set_id()?;

        if !self.series.contains_tags() {
            return Ok(());
        }

        let tag_set = self
            .resolved_tag_sets
            .get(tag_set_id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing tag set"))?;

        record.data.reserve(tag_set.len());

        for tag in tag_set {
            let src = tag
                .decoder
                .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

            if tag.is_kept {
                let value = self::data::read_value(src, tag.key.ty())?;
                record.data.push((tag.key.tag(), value));
            }
        }

        Ok(())
    }

    fn read_tag_set_id(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .tag_set_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::TagSetIds))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapped_read(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        if self.reads(DataSeries::FeatureCounts) {
            let feature_count = self.read_feature_count()?;
            record.features.reserve(feature_count);

            let mut prev_position = 0;

            for _ in 0..feature_count {
                let feature = self.read_feature(prev_position)?;
                prev_position = usize::from(feature.position());
                record.features.push(feature);
            }
        }

        if self.reads(DataSeries::MappingQualities) {
            record.mapping_quality = self.read_mapping_quality()?;
        }

        if record.cram_flags.quality_scores_are_stored_as_array()
            && self.reads(DataSeries::QualityScores)
        {
            record.quality_scores = self.read_quality_scores(record.read_length)?;
        }

        Ok(())
    }

    fn read_feature_count(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .feature_counts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeatureCounts))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_feature(&mut self, prev_position: usize) -> io::Result<Feature<'c>> {
        use feature::Code;

        let code = self.read_feature_code()?;

        let delta = self.read_feature_position_delta()?;
        let position = Position::try_from(prev_position + delta)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        match code {
            Code::Bases => {
                let bases = self
                    .read_or_default(DataSeries::StretchesOfBases, Self::read_stretches_of_bases)?;
                Ok(Feature::Bases { position, bases })
            }
            Code::Scores => {
                let quality_scores = self.read_or_default(
                    DataSeries::StretchesOfQualityScores,
                    Self::read_stretches_of_quality_scores,
                )?;

                Ok(Feature::Scores {
                    position,
                    quality_scores,
                })
            }
            Code::ReadBase => {
                let base = self.read_or_default(DataSeries::Bases, Self::read_base)?;
                let quality_score =
                    self.read_or_default(DataSeries::QualityScores, Self::read_quality_score)?;

                Ok(Feature::ReadBase {
                    position,
                    base,
                    quality_score,
                })
            }
            Code::Substitution => {
                let code = self.read_or_default(
                    DataSeries::BaseSubstitutionCodes,
                    Self::read_base_substitution_code,
                )?;
                Ok(Feature::Substitution { position, code })
            }
            Code::Insertion => {
                let bases =
                    self.read_or_default(DataSeries::InsertionBases, Self::read_insertion_bases)?;
                Ok(Feature::Insertion { position, bases })
            }
            Code::Deletion => {
                let len =
                    self.read_or_default(DataSeries::DeletionLengths, Self::read_deletion_length)?;
                Ok(Feature::Deletion { position, len })
            }
            Code::InsertBase => {
                let base = self.read_or_default(DataSeries::Bases, Self::read_base)?;
                Ok(Feature::InsertBase { position, base })
            }
            Code::QualityScore => {
                let quality_score =
                    self.read_or_default(DataSeries::QualityScores, Self::read_quality_score)?;

                Ok(Feature::QualityScore {
                    position,
                    quality_score,
                })
            }
            Code::ReferenceSkip => {
                let len = self.read_or_default(
                    DataSeries::ReferenceSkipLengths,
                    Self::read_reference_skip_length,
                )?;
                Ok(Feature::ReferenceSkip { position, len })
            }
            Code::SoftClip => {
                let bases =
                    self.read_or_default(DataSeries::SoftClipBases, Self::read_soft_clip_bases)?;
                Ok(Feature::SoftClip { position, bases })
            }
            Code::Padding => {
                let len =
                    self.read_or_default(DataSeries::PaddingLengths, Self::read_padding_length)?;
                Ok(Feature::Padding { position, len })
            }
            Code::HardClip => {
                let len =
                    self.read_or_default(DataSeries::HardClipLengths, Self::read_hard_clip_length)?;
                Ok(Feature::HardClip { position, len })
            }
        }
    }

    /// Reads a feature's value from `data_series`, or returns a placeholder if it is not read.
    fn read_or_default<T, F>(&mut self, data_series: DataSeries, read: F) -> io::Result<T>
    where
        T: Default,
        F: FnOnce(&mut Self) -> io::Result<T>,
    {
        if self.reads(data_series) {
            read(self)
        } else {
            Ok(T::default())
        }
    }

    fn read_feature_code(&mut self) -> io::Result<feature::Code> {
        self.compression_header
            .data_series_encodings()
            .feature_codes()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeatureCodes))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| {
                feature::Code::try_from(id)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_feature_position_delta(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .feature_position_deltas()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeaturePositionDeltas))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_stretches_of_bases(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::StretchesOfBases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_stretches_of_quality_scores(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_quality_scores()
            .ok_or_else(|| {
                missing_data_series_encoding_error(DataSeries::StretchesOfQualityScores)
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Bases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_quality_score(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::QualityScores))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base_substitution_code(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .base_substitution_codes()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::BaseSubstitutionCodes))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_insertion_bases(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .insertion_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::InsertionBases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_deletion_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .deletion_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::DeletionLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_reference_skip_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .reference_skip_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReferenceSkipLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_soft_clip_bases(&mut self) -> io::Result<Cow<'c, [u8]>> {
        self.compression_header
            .data_series_encodings()
            .soft_clip_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::SoftClipBases))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_padding_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .padding_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::PaddingLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_hard_clip_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .hard_clip_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::HardClipLengths))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapping_quality(
        &mut self,
    ) -> io::Result<Option<sam::alignment::record::MappingQuality>> {
        self.compression_header
            .data_series_encodings()
            .mapping_qualities()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MappingQualities))?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::alignment::record::MappingQuality::new)
    }

    fn read_unmapped_read(&mut self, record: &mut Record<'c>) -> io::Result<()> {
        // Unmapped records have no mapping quality data series. HTSlib decodes them as 0.
        record.mapping_quality = sam::alignment::record::MappingQuality::new(0);

        if self.reads(DataSeries::Bases) {
            record.sequence = self.read_sequence(record.read_length)?;
        }

        if record.cram_flags.quality_scores_are_stored_as_array()
            && self.reads(DataSeries::QualityScores)
        {
            record.quality_scores = self.read_quality_scores(record.read_length)?;
        }

        Ok(())
    }

    fn read_sequence(&mut self, read_length: usize) -> io::Result<Cow<'c, [u8]>> {
        let encoding = self
            .compression_header
            .data_series_encodings()
            .bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Bases))?;

        encoding.get().decode_take(
            &mut self.core_data_reader,
            &mut self.external_data_readers,
            read_length,
        )
    }

    fn read_quality_scores(&mut self, read_length: usize) -> io::Result<Cow<'c, [u8]>> {
        const MISSING: u8 = 0xff;

        let encoding = self
            .compression_header
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::QualityScores))?;

        let src = encoding.get().decode_take(
            &mut self.core_data_reader,
            &mut self.external_data_readers,
            read_length,
        )?;

        if src.iter().all(|&n| n == MISSING) {
            Ok(Cow::from(&[]))
        } else {
            Ok(src)
        }
    }
}

/// A tag of a tag set that is read.
struct ResolvedTag<'ch> {
    key: tag_sets::Key,
    decoder: TagDecoder<'ch>,
    /// Whether the value is kept, or only read to keep a shared block in order.
    is_kept: bool,
}

/// How a tag's values are read: directly from an external block for the encodings htslib writes,
/// otherwise through the generic byte array codec.
enum TagDecoder<'ch> {
    /// Values of a fixed length (a single-symbol Huffman length) from one external block.
    FixedLength {
        len: usize,
        reader: usize,
    },
    /// Values ending at a stop byte in one external block.
    Stop {
        stop_byte: u8,
        reader: usize,
    },
    Generic(&'ch Encoding<ByteArray>),
}

impl<'ch> TagDecoder<'ch> {
    fn new(encoding: &'ch Encoding<ByteArray>, readers: &ExternalDataReaders<'_>) -> Self {
        use crate::cram::container::compression_header::encoding::codec::{Byte, Integer};

        match encoding.get() {
            ByteArray::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => {
                if let (Integer::Huffman { alphabet, .. }, Byte::External { block_content_id }) =
                    (len_encoding.get(), value_encoding.get())
                    && let [len] = alphabet[..]
                    && let Ok(len) = usize::try_from(len)
                    && let Some(reader) = readers.index_of(block_content_id)
                {
                    return Self::FixedLength { len, reader };
                }
            }
            ByteArray::ByteArrayStop {
                stop_byte,
                block_content_id,
            } => {
                if let Some(reader) = readers.index_of(block_content_id) {
                    return Self::Stop {
                        stop_byte: *stop_byte,
                        reader,
                    };
                }
            }
        }

        Self::Generic(encoding)
    }

    /// Reads a value as [`ByteArray`]'s `decode` does.
    #[inline]
    fn decode<'de>(
        &self,
        core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
    ) -> io::Result<Cow<'de, [u8]>> {
        match *self {
            Self::FixedLength { len, reader } => {
                let src = external_data_readers.get_mut_by_index(reader);
                let buf = src
                    .split_off(..len)
                    .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;
                Ok(Cow::from(buf))
            }
            Self::Stop { stop_byte, reader } => {
                let src = external_data_readers.get_mut_by_index(reader);

                let Some(i) = src.iter().position(|&b| b == stop_byte) else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "missing byte array stop byte",
                    ));
                };

                let (buf, rest) = src.split_at(i);
                *src = &rest[1..];

                Ok(Cow::from(buf))
            }
            Self::Generic(encoding) => encoding.decode(core_data_reader, external_data_readers),
        }
    }
}

#[expect(clippy::type_complexity)]
fn resolve_tag_sets(
    compression_header: &CompressionHeader,
) -> io::Result<Vec<Vec<(tag_sets::Key, &Encoding<ByteArray>)>>> {
    compression_header
        .preservation_map()
        .tag_sets()
        .iter()
        .map(|tag_set| {
            tag_set
                .iter()
                .map(|&key| {
                    let id = block::ContentId::from(key);

                    let encoding =
                        compression_header.tag_encodings().get(&id).ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                ReadRecordError::MissingTagEncoding(key),
                            )
                        })?;

                    Ok((key, encoding))
                })
                .collect()
        })
        .collect()
}

fn missing_data_series_encoding_error(data_series: DataSeries) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        ReadRecordError::MissingDataSeriesEncoding(data_series),
    )
}
