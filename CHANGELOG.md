# Changelog

Notable changes to pasta, ntools, and vtools. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The three crates share one version.

## [Unreleased]

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
