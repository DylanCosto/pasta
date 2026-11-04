use crate::fasta;

use super::compression;
use crate::cram::{container::BlockContentEncoderMap, file_definition::Version};

const DEFAULT_RECORDS_PER_SLICE: usize = 10240;

#[derive(Clone, Debug)]
pub struct Context {
    pub reference_sequence_repository: fasta::Repository,
    pub records_per_slice: usize,
    /// A slice also ends once its records hold this many bases.
    pub bases_per_slice: usize,
    pub preserve_read_names: bool,
    pub encode_alignment_start_positions_as_deltas: bool,
    pub version: Version,
    pub block_content_encoder_map: BlockContentEncoderMap,
    /// How blocks are compressed, with methods chosen per block by the writer's
    /// [`compression::Compressor`]; without it, blocks use [`Self::block_content_encoder_map`].
    pub compression: Option<compression::Options>,
}

impl Context {
    /// Whether records are stored against reference sequences whose `@SQ` `M5` does not match
    /// them ([`compression::Options::ignore_md5`]).
    pub fn ignores_md5(&self) -> bool {
        self.compression.is_some_and(|options| options.ignore_md5)
    }

    /// Whether `MD` and `NM` tags are kept even when readers would regenerate them
    /// ([`compression::Options::store_md`], [`compression::Options::store_nm`]).
    pub fn stores_md_nm(&self) -> (bool, bool) {
        self.compression.map_or((false, false), |options| {
            (options.store_md, options.store_nm)
        })
    }
}

impl Default for Context {
    fn default() -> Self {
        Self {
            reference_sequence_repository: fasta::Repository::default(),
            records_per_slice: DEFAULT_RECORDS_PER_SLICE,
            bases_per_slice: usize::MAX,
            preserve_read_names: true,
            encode_alignment_start_positions_as_deltas: true,
            version: Version::default(),
            block_content_encoder_map: BlockContentEncoderMap::default(),
            compression: None,
        }
    }
}
