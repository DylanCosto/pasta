# pasta

Rust-native SAM/BAM/CRAM/VCF/BCF library, with ntools and vtools: samtools- and
bcftools-compatible command-line tools that put correct output first.

- [`pasta`](pasta): the library. SAM, BAM, and CRAM; VCF and BCF (`pasta::variant`); BGZF; BAI,
  CSI, CRAI, and TBI indexes; FASTA; and an HTSlib-compatible pileup engine. It began as a fork
  of [noodles](https://github.com/zaeleus/noodles) and is developed independently.
- [`ntools`](ntools): a samtools 1.24 replacement, with HTSlib's `bgzip`.
- [`vtools`](vtools): a bcftools 1.24 replacement.

## Install

Each [release](https://github.com/DylanCosto/pasta/releases) has prebuilt `ntools` and `vtools`
binaries, with SHA-256 checksums, for Linux x86_64 (`x86_64-linux`; glibc 2.28 or newer: RHEL 8,
Ubuntu 18.10 and later) and macOS on Apple silicon (`aarch64-macos`):

```sh
tar -xzf pasta-v0.2.0-x86_64-linux.tar.gz
pasta-v0.2.0-x86_64-linux/ntools --version
```

The macOS binaries are not signed: after downloading the archive with a browser, run
`xattr -d com.apple.quarantine pasta-v0.2.0-aarch64-macos.tar.gz` before unpacking it.

To build from source, with Rust 1.91 or newer and a C compiler:

```sh
cargo install --locked --path ntools
cargo install --locked --path vtools
```

The default build links libdeflate, which is C. `cargo build --release --no-default-features -p
ntools -p vtools` builds without any C code; deflate then uses zlib-rs, which gives the same
decompressed data but other compressed bytes.

## Policy

- Output must be right by the format specifications (SAM, BAM, CRAM, VCF, BCF, BGZF, and the
  index formats) and by what each command documents. A wrong result is a bug, whatever samtools
  or bcftools gives.
- Commands, options, output formats, and exit statuses follow samtools, bcftools, and HTSlib
  1.24, so that scripts keep working. Where those tools have a bug (a crash, a read of unset
  memory, a wrong result), ntools and vtools do the right thing rather than copy it; each tool's
  README lists these and the other intentional differences.
- Compressed bytes can differ from samtools' and bcftools', since they depend on the deflate
  library and on how blocks are cut.

## Commands

Each tool's README lists every option: [ntools](ntools/README.md#commands) and
[vtools](vtools/README.md#commands).

| Tool | Status | Commands |
|---|---|---|
| ntools | Supported | `view`, `index`, `merge`, `fqidx`, `flagstat`, `idxstats`, `depth`, `coverage`, `fastq`, `fasta`, `quickcheck`, `reheader`, `cat`, `split`, `collate`, `fixmate`, `dict`, `bedcov`, `addreplacerg`, `markdup`, `stats`, `mpileup`, `calmd`, `ampliconclip`, `import`, `reset`, `checksum`, `head`, `flags`, `samples`, `consensus`, `cram-size`, `reference`, `rmdup`, `help`, and `bgzip` (HTSlib's) |
| ntools | Partial | `sort`: no minimiser order (`-M`, `-I`, `-K`, `-w`). `faidx`: no `--write-index` with BGZF output |
| ntools | Not provided | `ampliconstats`, `phase`, `targetcut`, `depad`, `tview`, `version` (use `ntools --version`) |
| vtools | Supported | `view`, `filter`, `head`, `annotate`, `query`, `convert`, `index`, `isec`, `sort`, `reheader`, `concat`, `tabix`, `norm`, `consensus`, `stats`, `call`, `mpileup`, `gtcheck`, `merge`, `roh`, `csq`, `help`, `version` |
| vtools | Partial | Plugins: the built-in `fill-tags`, `setGT`, `split-vep`, and `counts`; no shared-library plugins |
| vtools | Not provided | `cnv`, `polysomy`, `som`; Perl functions in filter expressions |

## Formats

- Alignments: SAM (plain, gzip, or BGZF), BAM, and CRAM 3.0 and 3.1 (the default for CRAM
  output), read and written. CRAM 2.1 and 4.x are not read. CRAM references come from
  `--reference`/`-T`, or from the `@SQ` `UR` field when it names a local FASTA with a `.fai`;
  a sequence dictionary beside the FASTA (`ref.dict`) lets ntools check every slice against it.
  Lookup by MD5 (`REF_PATH`, `REF_CACHE`, reference servers) is not implemented. Every CRAM
  codec HTSlib 1.24 uses is read, and the writer chooses each block's method as HTSlib does.
- Variants: VCF (plain, gzip, or BGZF) and BCFv2.2 (compressed or not), read and written.
- Indexes: BAI, CSI, CRAI, TBI, GZI, and FAI, read and written. An index under another name can
  be given as `FILE##idx##INDEX`.
- Inputs: local files, standard input (`-`), pipes, and FIFOs. Remote names (`http://`, `s3://`,
  ...) are not opened; download such files first.

## Platforms

| Platform | Status |
|---|---|
| Linux x86_64 | Tested in CI; release binaries need glibc 2.28 or later |
| macOS aarch64 (Apple silicon) | Tested in CI; release binaries |
| Linux aarch64, other Unix systems | Not tested |
| Windows | Not supported (ntools uses Unix-only APIs) |

SIMD code (AVX2 for CRAM quality decoding and `vtools gtcheck`) is chosen at run time, with a
portable fallback, and results do not depend on the platform: `vtools call -c` computes its
80-bit sums in software where there is no x87 unit, and BAQ keeps the quality of a base it is
certain of everywhere (x86-64 samtools sets it to 0, arm64 samtools keeps it).

## Library

`pasta` shares its version with ntools and vtools and is not published to crates.io. Its
[README](pasta/README.md) lists the modules and the changes from noodles.

## Contributing

Bug reports and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md), and
[SECURITY.md](SECURITY.md) for vulnerabilities. Changes are listed in [CHANGELOG.md](CHANGELOG.md).

## License

MIT; see [LICENSE](LICENSE) and [NOTICE](NOTICE), which credits the projects pasta derives from.
