# Changelog

Notable changes to pasta, ntools, and vtools. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The three crates share one version.

## [Unreleased]

### Added

- pasta: `pileup::Read::quality_score`, the quality of one base: 0xff (missing) past the stored
  qualities, as for every base of a read whose `SEQ` is `*`.

### Changed

- pasta: `pileup::Read::base` gives 15 (`N`) past the stored sequence instead of panicking, so
  a read whose `SEQ` is `*` and whose CIGAR consumes query bases piles up as unknown bases.
- ntools: `calmd` and `mpileup` now fail with "CIGAR and query sequence lengths differ" on a
  mapped read whose CIGAR and sequence differ in length, as samtools does, where they used to
  write output for it. `ntools view` still reads such a read.

### Fixed

- pasta: a BAI, CSI, or tabix index with a bin ID past the last bin of its depth (other than the
  metadata pseudo-bin) fails to load with `InvalidData` ("invalid bin ID"), instead of loading
  and then panicking on the first query of that reference sequence.
- pasta: `ReferenceSequence::query` skips bins past the last bin of the depth instead of
  panicking on them.
- pasta: reading a BAI, CSI, or tabix index no longer reserves memory for the reference
  sequence or bin count it states before reading them (at most 65,536 up front), so a corrupt
  count fails with an error instead of aborting on a failed allocation of up to hundreds of GB.
- pasta: `variant::index::Index::parse` fails with `InvalidData` on a target, bin, chunk, or
  linear index count that the remaining data cannot hold, instead of reserving memory for it
  and aborting.
- pasta: reading a BAI, CSI, tabix, or VCF index (`variant::index::Index::parse`) with more
  reference sequences or bins than fit in memory, such as a small file of compressed zeros that
  states 100 million empty reference sequences, fails with `OutOfMemory` ("could not allocate
  N reference sequences: out of memory"; from `csi::io::Reader`, an `InvalidData` error caused
  by it) instead of aborting on a failed allocation. The lists grow as the entries are read.
- ntools, vtools: such an index makes the command fail with that error (vtools, as bcftools
  does: `Could not load local index file ... : Cannot allocate memory`) instead of aborting.
- pasta: `variant::index::report_load_error` reports an index too large to fit in memory as
  HTSlib does, `Cannot allocate memory`, not `Invalid argument`.
- ntools: a BAI, CSI, or tabix index with a bin ID past the last bin of its depth no longer
  loads, where samtools 1.24 ignores the bin. A region query (`view`, `mpileup -r`, and the
  others) fails with "invalid bin ID" instead of panicking, where samtools returns records; and
  `idxstats` warns that it cannot load the index and counts the records by reading the whole
  file, where samtools prints the counts in the index's metadata.
- pasta: a BAM record `block_size` below 32, the length of a record's fixed fields, is an error
  (`InvalidData`), as in HTSlib; a `block_size` of 0 ended the records as if at EOF.
- pasta: a corrupt BAM record `block_size` no longer allocates and zero-fills up to 4 GiB before
  the record is read; the record buffer grows with the bytes read, so a truncated record is an
  `UnexpectedEof` error at the cost of the bytes present. `ntools merge` reads its BAM inputs
  the same way.
- pasta: a corrupt BAM header `n_ref` or reference sequence `l_name` no longer reserves memory
  for the count it states (up to tens of GB, which aborted the process); a header that ends
  early is an `UnexpectedEof` error.
- pasta: a BAM record `block_size` above `i32::MAX` is an error (`InvalidData`), as in HTSlib; a
  few megabytes of BGZF could decompress to the 4 GiB it states and the reader read it all before
  failing. `ntools merge` rejects it the same way.
- pasta: a BAM record whose fixed fields give an `l_read_name` of 0, or lengths of name, CIGAR,
  sequence and qualities that do not fit in its `block_size`, is an `InvalidData` error found
  before the rest of the record is read, as in HTSlib; it was read in full (up to 2 GiB) first,
  then the overrun was an `UnexpectedEof` error and an `l_read_name` of 0 was accepted.
  `ntools merge` checks its BAM inputs the same way.
- pasta: the pileup no longer panics on a read whose CIGAR is a single deletion or skip (`5D`,
  `5N`), which HTSlib reads from before its CIGAR; it covers its positions as a deletion or skip
  at query position 0, as one that starts a longer CIGAR does.
- pasta: overlapping mates of which one has no sequence (`SEQ` `*`) are left as they are,
  instead of failing the pileup ("invalid CIGAR in overlapping mates") or panicking.
- pasta: `pileup::baq::cap_mapq` (`mpileup -C`) no longer panics on a read without its
  sequence; its bases are unknown and their qualities missing.
- pasta: `pileup::Read::new` refuses a mapped read whose CIGAR and sequence differ in length
  ("CIGAR and query sequence lengths differ", as HTSlib refuses to read it), on which
  `vtools mpileup` panicked; a read without its sequence is still accepted.
- vtools: `mpileup` no longer panics on a read with a zero-length match (`2M0M2D`, `0M2D` with
  `SEQ` `*`), which puts it in the pileup just past its last base; the base there is an `N` of
  quality 0, as samtools shows it, so the default `-Q 1` leaves it out.

## [0.9.0] - 2026-11-08

### Added

- pasta: `cram::NormalizedSequenceDigest`, which calculates the `M5` checksum of a reference
  sequence (`calculate_normalized_sequence_digest`) over its parts, given in turn.
- pasta: `bgzf::io::writer::CompressionLevel::VARIANTS`, level 6, the default for VCF and BCF.

### Changed

- pasta: CRAM is written as fast as samtools writes it single-threaded, faster with threads, and
  no larger. On GIAB HG003 chr20 (1.1 GB BAM, GRCh38), `ntools view -C` writes CRAM 3.1 in 56.7 s
  against samtools' 63.9 s single-threaded, and in 3.9 s against 12.3 s at `-@ 16`; CRAM 3.0 in
  57.0 s against 56.8 s, and 4.0 s against 12.2 s. On the whole 35x genome (44 GB BAM) at
  `-@ 16`, CRAM 3.1 takes 160.9 s against samtools' 223.9 s (0.2.1: 473 s), with 8% less CPU
  time than samtools. The files are within 0.05% of samtools' size, and samtools reads them back
  to the same records.
  - The read-name tokenizer is htscodecs' `tok3` encoder: each name is encoded against the name
    htscodecs picks (found with a trie, as htscodecs does), the prefix of a format it recognizes
    (Illumina, Ion Torrent, Oxford Nanopore, PacBio) is one token, and the token streams are laid
    out as htscodecs lays them out. rANS 4x8 and Nx16 normalize frequencies as htscodecs does,
    store data that does not compress as it is, and code the parts of a striped stream with the
    method that suits each. Given the same data, rANS and read-name blocks come out byte for byte
    as HTSlib 1.24 writes them. `Encoder::NameTokenizer` uses the settings of HTSlib's default
    level.
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

### Fixed

- pasta: read names that take more than 127 tokens are no longer coded with the name tokenizer,
  whose blocks HTSlib then could not read (it accepts 128 token positions); the method trials
  pick another method, and an `Encoder::NameTokenizer` chosen for the block gzips them.
- ntools: `mpileup --help` no longer says that `-M` is unsupported.

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
