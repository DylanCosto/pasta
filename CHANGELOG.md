# Changelog

Notable changes to pasta, ntools, and vtools. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The three crates share one version.

## [Unreleased]

### Added

- pasta: `bgzf::io::writer::CompressionLevel::VARIANTS`, level 6, the default for VCF and BCF.

### Changed

- pasta: CRAM 3.1 read names are compressed as HTSlib compresses them. The name tokenizer tries
  the one or two rANS settings htscodecs tries for each kind of token at the compression level,
  instead of eight on every token stream, and stripes the streams of numbers. On 2 million GIAB
  HG003 NovaSeq reads, single-threaded, `ntools view -C` takes 11.4 s instead of 15.1 s
  (samtools: 7.2 s), and the file is 4.7% smaller, 0.3% smaller than samtools'. The compressed
  bytes of CRAM output change. `Encoder::NameTokenizer` uses the settings of HTSlib's default
  level.
- pasta compresses and decompresses DEFLATE (BGZF blocks, CRAM's gzip blocks) and computes CRC32
  with zflate 1.0, a Rust port of libdeflate, in place of libdeflate: the build has no C code, and
  the `libdeflate` feature is gone. On x86-64 (Xeon Gold 5418Y), reading BAM is 7–13% faster end
  to end (`view`, `flagstat`, `index`).
- BAM and compressed SAM are written at level 5 unless a level is given: 1.3–1.7 times as fast as
  by 0.2.0, which wrote at libdeflate's level 6, and under 1% larger. VCF.gz and BCF, and
  `ntools bgzip`'s output, are written at level 6 (`CompressionLevel::VARIANTS`, HTSlib's
  zlib-style 5): 1.3–1.4 times as fast as at HTSlib's default, libdeflate's level 7, and 5%
  larger. Given levels keep their meaning; `-O z6`, `-O b6` and `bgzip -l 6` give HTSlib's size.
  `CompressionLevel::default()` is 5.
- ntools: on Linux, `sort` and `collate` free the disk space of their temporary files while
  merging them, as far as the merge has read. A pipeline of commands that spill then needs
  temporary space for about one copy of the reads instead of one per command: on GIAB HG003
  chr20 (1.1 GB BAM), `collate | fixmate | sort | markdup` peaks at 2.0 GB of temporary files
  instead of 3.7 GB (samtools: 2.4 GB), in the same time.
- Library: `bgzf::io::writer::CompressionLevel` is on libdeflate's scale, 0 to 12, in every build.
  Builds without the `libdeflate` feature took zlib's 0 to 9, so `CompressionLevel::new(9)`
  there was the highest level and is now 9 of 12 (`CompressionLevel::BEST` is 12). The
  `From<CompressionLevel>` conversions to `i32` and to libdeflater's level are removed; use
  `CompressionLevel::get`.
- A BGZF block that inflates to fewer bytes than its ISIZE, or a CRAM gzip block to fewer than
  its size, is an error (`InvalidData`), as in HTSlib; such blocks were accepted when the CRC
  matched (BGZF) or always (CRAM gzip).
- Compressed bytes are zflate's, so `ntools bgzip` and saved CSI/TBI indexes no longer match
  HTSlib's byte for byte; they decompress to the same data, in the same blocks.
- On 64-bit Arm, the multithreaded BGZF reader's jobs each decode up to 8 queued blocks
  together, and vtools' frame reader its chunks, with zflate's decoder that decodes two at a
  time on one thread (`bgzf::io::raw::inflate_frames`); on an Apple M4 Pro, long-read BAMs
  inflate about 25% faster than one block at a time. A block is queued as soon as it is read, so
  blocks are batched only while the workers are behind, and none waits for the next to be read.
  Elsewhere a job still takes one block. When the reader pauses (a seek, `get_mut`), the reader
  thread stops at the blocks it has taken, rather than filling every free buffer first.

## [0.2.1] - 2026-11-05

### Added

- pasta: `bam::record::codec::decoder::decode_bases`, which decodes 4-bit BAM bases to their
  characters 32 at a time with SSSE3 (x86-64, chosen at run time) or NEON (aarch64).

### Changed

- pasta decodes BAM sequences with `decode_bases`: `bam::record::Sequence`'s `append_to` and
  conversion to a `RecordBuf` sequence, decoding to `RecordBuf`, and the SAM writer. So do
  ntools' SAM output (`view`), `fastq`/`fasta` and filter expressions: `ntools fastq -@ 8` takes
  19% less time on 35x Illumina reads, and `view` to SAM up to 3% less on long reads.

## [0.2.0] - 2026-11-04

The first public release.

- ntools: a samtools 1.24 replacement (35 commands, and HTSlib's `bgzip`), reading and writing
  SAM, BAM, and CRAM 3.0/3.1, with `-@` threads in every command that reads or writes BGZF or
  CRAM, including CRAM encoding and decoding.
- vtools: a bcftools 1.24 replacement (23 commands and four built-in plugins), reading and
  writing VCF and BCF, with `--threads` workers that parse, decode, format, and encode records.
- pasta: the library both tools are built on, a fork of noodles with a port of HTSlib 1.24's
  `vcf.c` (`pasta::variant`) and pileup engine.
- Correct output over agreement with samtools and bcftools: where those tools crash, read unset
  memory, or give a wrong result, ntools and vtools do the right thing. Each tool's README lists
  the intentional differences. For example, BAQ keeps the quality of a base it is certain of
  (x86-64 samtools sets it to 0), and CRAM slices without a reference checksum are checked
  against the reference's sequence dictionary (`ref.dict`) or `@SQ` `LN`.
- Prebuilt binaries for Linux x86_64 (glibc 2.28 or later) and macOS on Apple silicon.
