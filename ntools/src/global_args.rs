//! samtools' shared long options (`sam_opts.c`), for the commands whose samtools counterparts
//! accept them.

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use pasta::cram::{
    self as cram,
    io::{
        reader::MdNm,
        writer::compression::{self, Profile},
    },
};

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    commands::index::{IndexKind, index_file},
    output::{OutputFormat, split_index_path},
};

/// `--verbosity`, HTSlib's logging level. It changes only diagnostics, so it is accepted and
/// ignored.
#[derive(Debug, clap::Args)]
pub struct Verbosity {
    /// Ignored (HTSlib's logging level).
    #[arg(long = "verbosity", value_name = "INT", allow_negative_numbers = true)]
    verbosity: Option<i32>,
}

/// `--write-index` for commands that, like their samtools counterparts, write no output to
/// index: samtools accepts and ignores it.
#[derive(Debug, clap::Args)]
pub struct IgnoredWriteIndex {
    /// Ignored, as in samtools: there is no alignment output to index.
    #[arg(long = "write-index")]
    write_index: bool,
}

/// `--input-fmt-option OPT[=VAL]`, HTSlib's input format options.
#[derive(Debug, Default, clap::Args)]
pub struct InputFormatOptions {
    /// Input format option: `decode_md=0` keeps CRAM MD/NM tags as stored, `reference=FILE` is
    /// --reference, and `nthreads=N` is ignored.
    #[arg(long = "input-fmt-option", value_name = "OPT[=VAL]")]
    input_format_options: Vec<String>,
}

impl InputFormatOptions {
    /// Options given by a command's own argument.
    pub fn new(input_format_options: Vec<String>) -> Self {
        Self {
            input_format_options,
        }
    }

    /// These options, after the ones a `--input-fmt FORMAT,OPT[=VAL]...` value adds; the format
    /// itself is detected, so its name is ignored.
    pub fn with_input_format(&self, input_format: Option<&str>) -> Self {
        let mut input_format_options: Vec<String> = input_format
            .map(|s| s.split(',').skip(1).map(str::to_owned).collect())
            .unwrap_or_default();

        input_format_options.extend(self.input_format_options.iter().cloned());

        Self {
            input_format_options,
        }
    }

    /// Checks the options, for commands that read no alignments they would change.
    pub fn validate(&self) -> anyhow::Result<()> {
        self.parse().map(|_| ())
    }

    /// Takes `reference=FILE` as the reference when `--reference` was not given.
    pub fn configure(&self, read_options: &mut ReadOptions) -> anyhow::Result<()> {
        for option in self.parse()? {
            if let InputOption::Reference(path) = option
                && read_options.reference.is_none()
            {
                read_options.reference = Some(path);
            }
        }

        Ok(())
    }

    /// Applies the options that change how an open reader decodes (`decode_md`).
    pub fn apply(&self, reader: &mut AlignmentReader) -> anyhow::Result<()> {
        for option in self.parse()? {
            if let InputOption::DecodeMd(md_nm) = option {
                reader.set_md_nm(md_nm);
            }
        }

        Ok(())
    }

    fn parse(&self) -> anyhow::Result<Vec<InputOption>> {
        self.input_format_options
            .iter()
            .map(|s| InputOption::parse(s))
            .collect()
    }
}

/// A supported HTSlib input format option.
enum InputOption {
    DecodeMd(MdNm),
    Reference(PathBuf),
    /// A performance option, ignored.
    Ignored,
}

impl InputOption {
    fn parse(s: &str) -> anyhow::Result<Self> {
        let (name, value) = s.split_once('=').unwrap_or((s, ""));

        match name {
            "decode_md" => {
                // HTSlib reads a missing value as 1 and the value with `strtol`.
                let n: i64 = if value.is_empty() {
                    1
                } else {
                    value
                        .parse()
                        .with_context(|| format!("invalid input format option: {s}"))?
                };

                Ok(Self::DecodeMd(if n == 0 {
                    MdNm::AsStored
                } else {
                    MdNm::Regenerate
                }))
            }
            "reference" if !value.is_empty() => Ok(Self::Reference(PathBuf::from(value))),
            "nthreads" => Ok(Self::Ignored),
            _ => bail!("unsupported input format option: {s}"),
        }
    }
}

/// `--output-fmt-option OPT[=VAL]`, HTSlib's output format options.
#[derive(Debug, Default, clap::Args)]
pub struct OutputFormatOptions {
    /// Output format option: `level=N` (compression level); for CRAM, also `version=3.0|3.1`,
    /// a profile (`fast`, `normal`, `small`, `archive`), `use_bzip2`, `use_lzma`, `use_rans`,
    /// `use_tok`, `use_fqz`, `use_arith`, `store_md`, `store_nm`, `ignore_md5` (`=0` or `=1`),
    /// `seqs_per_slice=N`, and `bases_per_slice=N`; `nthreads=N` is ignored. `-O FORMAT,OPT,...`
    /// also sets them.
    #[arg(long = "output-fmt-option", value_name = "OPT=VAL")]
    output_format_options: Vec<String>,
}

/// A CRAM codec that `use_*` options turn on or off.
#[derive(Clone, Copy)]
enum Codec {
    Bzip2,
    Lzma,
    Rans,
    Tok,
    Fqz,
    Arith,
}

/// One output format option.
enum OutputOption {
    Level(u8),
    Version(cram::file_definition::Version),
    Profile(Profile),
    UseCodec(Codec, bool),
    SeqsPerSlice(usize),
    BasesPerSlice(usize),
    StoreMd(bool),
    StoreNm(bool),
    IgnoreMd5(bool),
    Ignored,
}

impl OutputOption {
    fn parse(s: &str) -> anyhow::Result<Self> {
        let (key, value) = s.split_once('=').unwrap_or((s, ""));

        let flag = |value: &str| -> anyhow::Result<bool> {
            match value {
                "" | "1" => Ok(true),
                "0" => Ok(false),
                _ => bail!("invalid output format option: {s}"),
            }
        };

        Ok(match key {
            "level" => Self::Level(
                value
                    .parse()
                    .with_context(|| format!("invalid compression level: {value}"))?,
            ),
            "version" => Self::Version(match value {
                "3.0" | "3" => cram::file_definition::Version::new(3, 0),
                "3.1" => cram::file_definition::Version::new(3, 1),
                _ => bail!("unsupported CRAM version: {value}"),
            }),
            "fast" if value.is_empty() => Self::Profile(Profile::Fast),
            "normal" if value.is_empty() => Self::Profile(Profile::Normal),
            "small" if value.is_empty() => Self::Profile(Profile::Small),
            "archive" if value.is_empty() => Self::Profile(Profile::Archive),
            "use_bzip2" => Self::UseCodec(Codec::Bzip2, flag(value)?),
            "use_lzma" => Self::UseCodec(Codec::Lzma, flag(value)?),
            "use_rans" => Self::UseCodec(Codec::Rans, flag(value)?),
            "use_tok" => Self::UseCodec(Codec::Tok, flag(value)?),
            "use_fqz" => Self::UseCodec(Codec::Fqz, flag(value)?),
            "use_arith" => Self::UseCodec(Codec::Arith, flag(value)?),
            "seqs_per_slice" => Self::SeqsPerSlice(
                value
                    .parse()
                    .ok()
                    .filter(|&n| n > 0)
                    .with_context(|| format!("invalid output format option: {s}"))?,
            ),
            "bases_per_slice" => Self::BasesPerSlice(
                value
                    .parse()
                    .ok()
                    .filter(|&n| n > 0)
                    .with_context(|| format!("invalid output format option: {s}"))?,
            ),
            "store_md" => Self::StoreMd(flag(value)?),
            "store_nm" => Self::StoreNm(flag(value)?),
            "ignore_md5" => Self::IgnoreMd5(flag(value)?),
            "nthreads" => Self::Ignored,
            _ => bail!("unsupported output format option: {s}"),
        })
    }
}

impl OutputFormatOptions {
    fn parse(&self) -> anyhow::Result<Vec<OutputOption>> {
        self.output_format_options
            .iter()
            .map(|s| OutputOption::parse(s))
            .collect()
    }

    /// The compression level from `level=N`, if given.
    pub fn explicit_level(&self) -> anyhow::Result<Option<u8>> {
        let mut level = None;

        for option in self.parse()? {
            if let OutputOption::Level(n) = option {
                level = Some(n);
            }
        }

        Ok(level)
    }

    /// The compression level from `level=N`, or `default`.
    pub fn level(&self, default: u8) -> anyhow::Result<u8> {
        let mut level = default;

        for option in self.parse()? {
            if let OutputOption::Level(n) = option {
                level = n;
            }
        }

        Ok(level)
    }

    /// The CRAM compression options: CRAM 3.1 with the normal profile, as in samtools, changed by
    /// the options given. As in HTSlib, a profile sets the level only when none is given.
    pub fn cram(&self) -> anyhow::Result<compression::Options> {
        let options = self.parse()?;

        let version = options
            .iter()
            .filter_map(|option| match option {
                OutputOption::Version(version) => Some(*version),
                _ => None,
            })
            .next_back()
            .unwrap_or(cram::file_definition::Version::new(3, 1));

        let profile = options
            .iter()
            .filter_map(|option| match option {
                OutputOption::Profile(profile) => Some(*profile),
                _ => None,
            })
            .next_back()
            .unwrap_or_default();

        let mut cram = compression::Options::new(version, profile);
        let bases_given = options
            .iter()
            .any(|option| matches!(option, OutputOption::BasesPerSlice(_)));

        for option in options {
            match option {
                OutputOption::Level(level) => cram.level = level,
                OutputOption::UseCodec(codec, value) => {
                    let use_codec = match codec {
                        Codec::Bzip2 => &mut cram.use_bzip2,
                        Codec::Lzma => &mut cram.use_lzma,
                        Codec::Rans => &mut cram.use_rans,
                        Codec::Tok => &mut cram.use_tok,
                        Codec::Fqz => &mut cram.use_fqz,
                        Codec::Arith => &mut cram.use_arith,
                    };
                    *use_codec = value;
                }
                // As in HTSlib, records per slice also set the bases per slice unless given.
                OutputOption::SeqsPerSlice(n) => {
                    cram.records_per_slice = n;

                    if !bases_given {
                        cram.bases_per_slice = n * compression::BASES_PER_RECORD;
                    }
                }
                OutputOption::BasesPerSlice(n) => cram.bases_per_slice = n,
                OutputOption::StoreMd(value) => cram.store_md = value,
                OutputOption::StoreNm(value) => cram.store_nm = value,
                OutputOption::IgnoreMd5(value) => cram.ignore_md5 = value,
                OutputOption::Version(_) | OutputOption::Profile(_) | OutputOption::Ignored => {}
            }
        }

        Ok(cram)
    }
}

/// An index to build for an alignment output file once it is written (`--write-index`,
/// samtools' `auto_index`).
pub struct IndexRequest {
    output: PathBuf,
    index: Option<PathBuf>,
    kind: IndexKind,
}

impl IndexRequest {
    /// Checks that `output` (`None` or `-` for stdout), without any `##idx##INDEX` suffix, can be
    /// indexed. The index goes to `index`, or next to the output: CSI for BAM unless `index` ends
    /// in `.bai`, and CRAI for CRAM.
    pub fn new(
        output: Option<&Path>,
        index: Option<PathBuf>,
        format: OutputFormat,
    ) -> anyhow::Result<Self> {
        let Some(output) = output.filter(|p| *p != Path::new("-")) else {
            bail!("--write-index needs an output file");
        };

        let is_bai = index
            .as_deref()
            .is_some_and(|p| p.as_os_str().as_encoded_bytes().ends_with(b".bai"));

        let kind = match format {
            OutputFormat::Bam if is_bai => IndexKind::Default,
            OutputFormat::Bam => IndexKind::Csi,
            OutputFormat::Cram => IndexKind::Default,
            OutputFormat::Sam | OutputFormat::SamGz => {
                bail!("--write-index needs BAM or CRAM output")
            }
        };

        Ok(Self {
            output: output.to_path_buf(),
            index,
            kind,
        })
    }

    /// Checks that an output given as `-o` names it (`FILE` or `FILE##idx##INDEX`, or `None` or
    /// `-` for stdout) can be indexed; see [`Self::new`].
    pub fn for_output(output: Option<&Path>, format: OutputFormat) -> anyhow::Result<Self> {
        let (output, index) = match output {
            Some(path) => {
                let (output, index) = split_index_path(path);
                (Some(output), index)
            }
            None => (None, None),
        };

        Self::new(output.as_deref(), index, format)
    }

    /// Indexes the written output.
    pub fn write(self, threads: usize) -> anyhow::Result<()> {
        index_file(&self.output, self.index, self.kind, threads)
    }
}

/// Indexes an alignment output file, given as `-o` names it (`FILE` or `FILE##idx##INDEX`), as
/// samtools' `--write-index` does; see [`IndexRequest::new`].
pub fn write_index(
    output: Option<&Path>,
    format: OutputFormat,
    threads: usize,
) -> anyhow::Result<()> {
    IndexRequest::for_output(output, format)?.write(threads)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request(output: &str, format: OutputFormat) -> anyhow::Result<IndexRequest> {
        IndexRequest::for_output(Some(Path::new(output)), format)
    }

    #[test]
    fn test_index_request() -> anyhow::Result<()> {
        let bam = request("out.bam", OutputFormat::Bam)?;
        assert_eq!(
            (bam.output, bam.index, bam.kind),
            (PathBuf::from("out.bam"), None, IndexKind::Csi)
        );

        let bai = request("out.bam##idx##custom.bai", OutputFormat::Bam)?;
        assert_eq!(
            (bai.output, bai.index, bai.kind),
            (
                PathBuf::from("out.bam"),
                Some(PathBuf::from("custom.bai")),
                IndexKind::Default
            )
        );

        let csi = request("out.bam##idx##custom.idx", OutputFormat::Bam)?;
        assert_eq!(csi.kind, IndexKind::Csi);

        let crai = request("out.cram##idx##custom.crai", OutputFormat::Cram)?;
        assert_eq!(
            (crai.output, crai.kind),
            (PathBuf::from("out.cram"), IndexKind::Default)
        );

        assert!(request("out.sam", OutputFormat::Sam).is_err());
        assert!(request("-", OutputFormat::Bam).is_err());
        assert!(IndexRequest::for_output(None, OutputFormat::Bam).is_err());

        Ok(())
    }
}
