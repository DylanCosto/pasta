# Changelog

Notable changes to pasta, ntools, and vtools. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The three crates share one version.

## [Unreleased]

### Added

- pasta: `cram::NormalizedSequenceDigest`, which calculates the `M5` checksum of a reference
  sequence (`calculate_normalized_sequence_digest`) over its parts, given in turn.

### Changed

- pasta: CRAM is written as fast as samtools writes it single-threaded, faster with threads, and
  no larger. On GIAB HG003 chr20 (1.1 GB BAM, GRCh38), `ntools view -C` writes CRAM 3.1 in 56.7 s
  against samtools' 63.9 s single-threaded (0.9.0: 85.5 s), and in 3.9 s against 12.3 s at
  `-@ 16` (0.9.0: 16.7 s); CRAM 3.0 in 57.0 s against 56.8 s, and 4.0 s against 12.2 s. On the
  whole 35x genome (44 GB BAM) at `-@ 16`, CRAM 3.1 takes 160.9 s against 223.9 s (0.9.0:
  270 s), with 8% less CPU time than samtools. The files are within 0.05% of samtools' size
  (0.9.0's CRAM 3.1 was 5% larger), and samtools reads them back to the same records.
  - The read-name tokenizer is htscodecs' `tok3` encoder: each name is encoded against the name
    htscodecs picks (found with a trie, as htscodecs does), the prefix of a format it recognizes
    (Illumina, Ion Torrent, Oxford Nanopore, PacBio) is one token, and the token streams are laid
    out as htscodecs lays them out. rANS 4x8 and Nx16 normalize frequencies as htscodecs does,
    store data that does not compress as it is, and code the parts of a striped stream with the
    method that suits each. Given the same data, rANS and read-name blocks come out byte for byte
    as HTSlib 1.24 writes them.
  - Tags are stored as HTSlib stores them, from the BAM record's own bytes: fixed-size values
    without a length each, strings ended by a tab.
  - The `M5` checksums of a header's `@SQ` lines that lack them are calculated on all `-@` threads
    while records are encoded, reading each sequence a part at a time where the reference is
    indexed. MD5 is pasta's own, about 10% faster than the `md-5` crate, which pasta no longer
    depends on.
  - The rANS coders, the record conversion, and the slice encoder allocate and branch less.
    Records borrow their reference sequence instead of each counting a reference to it, on which
    the threads contended.
  The compressed bytes of CRAM output change.
- ntools: writing CRAM no longer calculates the reference's `M5` checksums before the first
  record; the CRAM writer does, alongside the records.
- ntools: `collate -m` is the memory for buffering records in all, not per thread: by default
  768M whatever `-@`, instead of 768M for each thread. On GIAB HG003 chr20 (1.1 GB BAM) at
  `-@ 16`, `collate -Ou` peaks at 910 MB instead of 7.7 GB and takes 17.9 s instead of 15.8 s
  (samtools: 34.2 s, 348 MB); at `-@ 4`, 860 MB instead of 4.1 GB, in 31.5 s instead of 21.6 s
  (samtools: 39.8 s). Inputs whose records fit in 768M are collated in memory, as before.
- ntools: `mpileup` sends records from its reader in batches of at most 4 MB as well as 4,096
  reads, so long reads no longer fill batches of hundreds of megabytes. On 9,784 nanopore reads,
  `mpileup -M` peaks at 70 MB instead of 520 MB (samtools: 44 MB), and `mpileup` at 57 MB
  instead of 430 MB, a little faster and with the same output.
- ntools: merging sorted runs (the temporary files of `sort` and `collate`, and `merge`'s inputs)
  compares records by a key computed once per record, as sorting in memory does, instead of, for
  `collate`, hashing both read names in every comparison. `collate -Ou -@ 4 -m 640M` of GIAB
  HG003 chr20 (1.1 GB BAM) to `/dev/null`, which spills, takes 27.0 s instead of 33.3 s
  (samtools: 27.9 s), with the same output.

## [0.9.0] - 2026-11-08

### Added

- pasta: `bgzf::io::writer::CompressionLevel::VARIANTS`, level 6, the default for VCF and BCF.

### Changed

- pasta: CRAM 3.1 read names are compressed as HTSlib compresses them. The name tokenizer tries
  the one or two rANS settings htscodecs tries for each kind of token at the compression level,
  instead of eight on every token stream, and stripes the streams of numbers. It also no longer
  allocates for each token, and, like htscodecs, codes a repeat of a slice's first name as a
  duplicate. On 2 million GIAB HG003 NovaSeq reads, single-threaded, `ntools view -C` takes 8.9 s
  instead of 15.1 s, and the file is 4.7% smaller; on the whole 35x genome at `-@ 16`, it takes
  270 s instead of 473 s (samtools: 227 s). The compressed bytes of CRAM output change.
  `Encoder::NameTokenizer` uses the settings of HTSlib's default level.
- pasta compresses and decompresses DEFLATE (BGZF blocks, CRAM's gzip blocks) and computes CRC32
  with zflate 1.0, a Rust port of libdeflate, in place of libdeflate: the build has no C code, and
  the `libdeflate` feature is gone. On x86-64 (Xeon Gold 5418Y), reading BAM is 7–13% faster end
  to end (`view`, `flagstat`, `index`).
- BAM and compressed SAM are written at level 5 unless a level is given: commands that write BAM
  are 1.2–1.7 times as fast as in 0.2.0, which wrote at libdeflate's level 6, and their files
  under 1% larger. VCF.gz and BCF, and `ntools bgzip`'s output, are written at level 6
  (`CompressionLevel::VARIANTS`, HTSlib's zlib-style 5): 1.1–1.6 times as fast as at HTSlib's
  default, libdeflate's level 7, and 5% larger. Given levels keep their meaning; `-O z6`, `-O b6`
  and `bgzip -l 6` give HTSlib's size. `CompressionLevel::default()` is 5.
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
