//! Alignment writer builder.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use crate::bam;
use crate::bgzf;
use crate::cram::{self as cram, container::BlockContentEncoderMap};
use crate::fasta;
use crate::sam;

use super::Writer;
use crate::util::alignment::io::{CompressionMethod, Format};

/// An alignment writer builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
    block_content_encoder_map: BlockContentEncoderMap,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// If not set, a default compression method is selected depending on the format.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::util::alignment::{self, io::CompressionMethod};
    /// let builder = alignment::io::writer::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the output.
    ///
    /// If not set, a default format is used.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::util::alignment::{self, io::Format};
    /// let builder = alignment::io::writer::Builder::default()
    ///     .set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// This is only used when the output format is CRAM.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::fasta;
    /// use pasta::util::alignment::{self, io::Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::io::writer::Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets the block content-encoder map.
    ///
    /// This is only used when the output format is CRAM.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::cram::container::BlockContentEncoderMap;
    /// use pasta::util::alignment::{self, io::Format};
    ///
    /// let builder = alignment::io::writer::Builder::default()
    ///     .set_block_content_encoder_map(BlockContentEncoderMap::default());
    /// ```
    pub fn set_block_content_encoder_map(
        mut self,
        block_content_encoder_map: BlockContentEncoderMap,
    ) -> Self {
        self.block_content_encoder_map = block_content_encoder_map;
        self
    }

    /// Builds an alignment writer from a path.
    ///
    /// If the format is not set, it is detected from the path extension. If the compression
    /// method is not set, it is the format's standard one: BGZF for BAM, none for CRAM, and for
    /// SAM, BGZF when the path ends in `.gz` or `.bgz`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use pasta::util::alignment::{self, io::Format};
    /// let writer = alignment::io::writer::Builder::default()
    ///     .build_from_path("out.sam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer<File>>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src);
        }

        // BAM and CRAM have one standard compression each, which `build_from_writer` picks; the
        // extension only chooses between plain and bgzipped SAM. A format that was set explicitly
        // thus keeps its compression whatever the path is called (e.g., `out.bam.tmp`).
        if self.compression_method.is_none() && matches!(self.format, None | Some(Format::Sam)) {
            self.compression_method = Some(detect_compression_method_from_path_extension(src));
        }

        File::create(src).and_then(|writer| self.build_from_writer(writer))
    }

    /// Builds an alignment writer from a writer.
    ///
    /// If the format is not set, a default format is used. If the compression method is not set, a
    /// default one is determined by the format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use pasta::util::alignment::{self, io::Format};
    /// let writer = alignment::io::writer::Builder::default()
    ///     .build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> io::Result<Writer<W>>
    where
        W: Write,
    {
        use super::Inner;

        let format = self.format.unwrap_or(Format::Sam);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => match format {
                Format::Sam | Format::Cram => None,
                Format::Bam => Some(CompressionMethod::Bgzf),
            },
        };

        let inner = match (format, compression_method) {
            (Format::Sam, None) => Inner::Sam(sam::io::Writer::new(BufWriter::new(writer))),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                Inner::SamGz(sam::io::Writer::new(bgzf::io::Writer::new(writer)))
            }
            (Format::Bam, None) => Inner::BamRaw(bam::io::Writer::from(BufWriter::new(writer))),
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                Inner::Bam(bam::io::Writer::new(writer))
            }
            (Format::Cram, None) => Inner::Cram(Box::new(
                cram::io::writer::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .set_block_content_encoder_map(self.block_content_encoder_map)
                    .build_from_writer(writer),
            )),
            (Format::Cram, Some(CompressionMethod::Bgzf)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM cannot be bgzip-compressed",
                ));
            }
        };

        Ok(Writer(inner))
    }
}

pub(crate) fn detect_compression_method_from_path_extension<P>(path: P) -> Option<CompressionMethod>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("bam" | "gz" | "bgz") => Some(CompressionMethod::Bgzf),
        _ => None,
    }
}

pub(crate) fn detect_format_from_path_extension<P>(path: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();

    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        Some("gz" | "bgz") => {
            let file_stem = path.file_stem().and_then(|stem| stem.to_str())?;
            file_stem.ends_with("sam").then_some(Format::Sam)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_compression_method_from_path_extension() {
        assert_eq!(
            detect_compression_method_from_path_extension("out.sam"),
            None
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.sam.gz"),
            Some(CompressionMethod::Bgzf)
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.sam.bgz"),
            Some(CompressionMethod::Bgzf)
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.bam"),
            Some(CompressionMethod::Bgzf)
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.cram"),
            None
        );

        assert!(detect_format_from_path_extension("out.fa").is_none());
    }

    #[test]
    fn an_explicit_format_keeps_its_compression_whatever_the_extension() -> io::Result<()> {
        use crate::bgzf::io::writer::BGZF_EOF;

        let dir = std::env::temp_dir().join(format!("pasta-writer-builder-{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let header = sam::Header::default();

        let write = |format, name: &str| -> io::Result<Vec<u8>> {
            let path = dir.join(name);
            let mut writer = Builder::default()
                .set_format(format)
                .build_from_path(&path)?;
            writer.write_header(&header)?;
            writer.finish(&header)?;
            drop(writer);
            std::fs::read(&path)
        };

        // BAM is BGZF-compressed, whatever the path's extension.
        for name in ["out.bam", "out.bam.tmp", "out"] {
            let data = write(Format::Bam, name)?;
            assert!(data.starts_with(&[0x1f, 0x8b]), "{name} is BGZF");
            assert!(
                data.ends_with(&BGZF_EOF),
                "{name} ends with the BGZF EOF block"
            );
        }

        // SAM is bgzipped only when the path says so.
        assert!(write(Format::Sam, "out.sam.gz")?.ends_with(&BGZF_EOF));
        assert!(!write(Format::Sam, "out.tmp")?.starts_with(&[0x1f, 0x8b]));

        // CRAM is never bgzipped, even when the path ends in `.gz`.
        assert!(write(Format::Cram, "out.cram.gz")?.starts_with(b"CRAM"));

        std::fs::remove_dir_all(&dir)?;

        Ok(())
    }

    #[test]
    fn test_detect_format_from_path_extension() {
        assert_eq!(
            detect_format_from_path_extension("out.sam"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.sam.gz"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.sam.bgz"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.bam"),
            Some(Format::Bam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.cram"),
            Some(Format::Cram)
        );

        assert!(detect_format_from_path_extension("out.fa").is_none());
    }
}
