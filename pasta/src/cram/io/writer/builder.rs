use std::{
    fs::File,
    io::{self, Write},
    num::NonZero,
    path::Path,
};

use crate::fasta;

use super::{Batch, Context, DEFAULT_SLICES_PER_CONTAINER, Writer, compression, parallel::Workers};
use crate::cram::{codecs::Encoder, container::BlockContentEncoderMap, file_definition::Version};

// § 7 "Container header structure" (2025-04-07): "record counter: 0-based sequential index of
// records in the file/stream."
const MIN_RECORD_COUNTER: u64 = 0;

/// Records per slice with a fixed block content-encoder map.
const FIXED_RECORDS_PER_SLICE: usize = 10240;

/// A CRAM writer builder.
#[derive(Default)]
pub struct Builder {
    context: Context,
    compression: Option<compression::Options>,
    block_content_encoder_map: Option<BlockContentEncoderMap>,
    worker_count: Option<NonZero<usize>>,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::writer::Builder;
    /// use pasta::fasta;
    ///
    /// let repository = fasta::Repository::default();
    /// let builder = Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.context.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets whether to preserve read names.
    ///
    /// If `false`, read names are discarded.
    ///
    /// The default is `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::writer::Builder;
    /// let builder = Builder::default().preserve_read_names(false);
    /// ```
    pub fn preserve_read_names(mut self, value: bool) -> Self {
        self.context.preserve_read_names = value;
        self
    }

    /// Sets whether to encode alignment start positions as deltas.
    ///
    /// If `false`, record alignment start positions are written with their actual values.
    ///
    /// The default is `true`.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::writer::Builder;
    /// let builder = Builder::default()
    ///     .encode_alignment_start_positions_as_deltas(false);
    /// ```
    pub fn encode_alignment_start_positions_as_deltas(mut self, value: bool) -> Self {
        self.context.encode_alignment_start_positions_as_deltas = value;
        self
    }

    /// Sets how blocks are compressed: the CRAM version, level, and methods allowed. Each block's
    /// method is then chosen as HTSlib does. This is the default (CRAM 3.1, normal profile), unless
    /// a block content-encoder map is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{
    ///     file_definition::Version,
    ///     io::writer::{
    ///         Builder,
    ///         compression::{Options, Profile},
    ///     },
    /// };
    ///
    /// let builder = Builder::default()
    ///     .set_compression_options(Options::new(Version::new(3, 0), Profile::Small));
    /// ```
    pub fn set_compression_options(mut self, options: compression::Options) -> Self {
        self.compression = Some(options);
        self
    }

    /// Sets a fixed encoder for each block content ID instead of choosing methods per block.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::{container::BlockContentEncoderMap, io::writer::Builder};
    ///
    /// let block_content_encoder_map = BlockContentEncoderMap::default();
    /// let builder = Builder::default()
    ///     .set_block_content_encoder_map(block_content_encoder_map);
    /// ```
    pub fn set_block_content_encoder_map(
        mut self,
        block_content_encoder_map: BlockContentEncoderMap,
    ) -> Self {
        self.block_content_encoder_map = Some(block_content_encoder_map);
        self
    }

    /// Sets the number of worker threads that encode and compress containers.
    ///
    /// By default, containers are encoded on the thread that writes the records. With workers,
    /// records are converted (when written with [`Writer::write_bam_record`]), encoded, and
    /// compressed on the workers, several containers at a time, and written in order. The output
    /// is the same either way. When the system cannot start that many threads, fewer are used, or
    /// none.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    /// use pasta::cram::io::writer::Builder;
    /// let builder = Builder::default().set_worker_count(NonZero::new(4).unwrap());
    /// ```
    pub fn set_worker_count(mut self, worker_count: NonZero<usize>) -> Self {
        self.worker_count = Some(worker_count);
        self
    }

    /// Builds a CRAM writer from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use pasta::cram::io::writer::Builder;
    /// let writer = Builder::default().build_from_path("out.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, dst: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        File::create(dst).map(|file| self.build_from_writer(file))
    }

    /// Builds a CRAM writer from a writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::cram::io::writer::Builder;
    /// let writer = Builder::default().build_from_writer(Vec::new());
    /// ```
    pub fn build_from_writer<W>(mut self, writer: W) -> Writer<W>
    where
        W: Write,
    {
        let mut compressor = None;

        let records_per_slice = match self.block_content_encoder_map {
            // A fixed encoder per block.
            Some(block_content_encoder_map) => {
                if uses_cram_3_1_codecs(&block_content_encoder_map) {
                    self.context.version = Version::new(3, 1);
                }

                self.context.block_content_encoder_map = block_content_encoder_map;
                FIXED_RECORDS_PER_SLICE
            }
            // Methods chosen per block, as HTSlib does.
            None => {
                let options = self.compression.unwrap_or_default();
                self.context.version = options.version;
                self.context.compression = Some(options);
                compressor = Some(compression::Compressor::new(options));
                self.context.bases_per_slice = options.bases_per_slice;
                options.records_per_slice
            }
        };

        let records_per_container = DEFAULT_SLICES_PER_CONTAINER * records_per_slice;
        self.context.records_per_slice = records_per_slice;

        Writer {
            inner: Some(writer),
            unfinished_header: None,
            context: self.context,
            compressor,
            batch: Batch::default(),
            records_per_container,
            record_counter: MIN_RECORD_COUNTER,
            bases: 0,
            missing_reference_sequence_names: Default::default(),
            pending_header: None,
            last_reference_sequence: None,
            read_groups: Default::default(),
            reference_sequences: Default::default(),
            // Without threads, containers are encoded on this one.
            workers: self.worker_count.and_then(|n| Workers::start(n).ok()),
            failed: false,
            many_references: false,
        }
    }
}

pub fn uses_cram_3_1_codecs(block_content_encoder_map: &BlockContentEncoderMap) -> bool {
    fn is_cram_3_1_codec(encoder: &Encoder) -> bool {
        matches!(
            encoder,
            Encoder::RansNx16(_) | Encoder::AdaptiveArithmeticCoding(_) | Encoder::NameTokenizer
        )
    }

    if let Some(encoder) = block_content_encoder_map.core_data_encoder()
        && is_cram_3_1_codec(encoder)
    {
        return true;
    }

    block_content_encoder_map
        .data_series_encoders()
        .iter()
        .chain(block_content_encoder_map.tag_values_encoders().values())
        .flatten()
        .any(is_cram_3_1_codec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uses_cram_3_1_codecs() {
        use crate::cram::codecs::rans_nx16::Flags;

        let block_content_encoder_map = BlockContentEncoderMap::default();
        assert!(!uses_cram_3_1_codecs(&block_content_encoder_map));

        let block_content_encoder_map = BlockContentEncoderMap::builder()
            .set_core_data_encoder(Some(Encoder::RansNx16(Flags::empty())))
            .build();
        assert!(uses_cram_3_1_codecs(&block_content_encoder_map));
    }
}
