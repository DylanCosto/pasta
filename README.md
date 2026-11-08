# pasta

Rust-native SAM/BAM/CRAM/VCF/BCF library, with ntools and vtools: samtools- and
bcftools-compatible command-line tools, inspired by [noodles](https://github.com/zaeleus/noodles).

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
tar -xzf pasta-v0.9.0-x86_64-linux.tar.gz
pasta-v0.9.0-x86_64-linux/ntools --version
```

The macOS binaries are not signed: after downloading the archive with a browser, run
`xattr -d com.apple.quarantine pasta-v0.9.0-aarch64-macos.tar.gz` before unpacking it.

To build from source, with Rust 1.91 or newer:

```sh
cargo install --locked --path ntools
cargo install --locked --path vtools
```

The build is pure Rust: DEFLATE (BGZF and CRAM's gzip blocks) and CRC32 go through
[zflate](https://crates.io/crates/zflate), a Rust port of libdeflate ([Compression](#compression)).

## Usage

`ntools` takes samtools' commands and options, and `vtools` takes bcftools', so existing
commands and scripts work with the name changed. `ntools --help` lists the commands and
`ntools COMMAND --help` gives a command's options; `vtools` alone lists its commands, and
`vtools COMMAND` without arguments prints the command's usage, as bcftools does.

Alignments:

```sh
ntools faidx ref.fa                                      # index the reference
ntools sort -@ 8 -o sample.bam sample.sam                # sort by coordinate, to BAM
ntools index sample.bam
ntools view -C -T ref.fa -@ 8 -o sample.cram sample.bam  # BAM to CRAM
ntools view -b -T ref.fa -@ 8 -o copy.bam sample.cram    # CRAM to BAM
ntools view -c -q 30 sample.bam chr1:1000-2000           # count reads in a region
ntools flagstat -@ 8 sample.cram
ntools stats -@ 8 sample.bam > sample.stats
ntools collate -O -u sample.bam | ntools fixmate -m -u - - |
    ntools sort -u - | ntools markdup - marked.bam       # mark duplicates
```

Variants:

```sh
vtools mpileup -f ref.fa -Ou sample.bam |
    vtools call -mv -Oz -o calls.vcf.gz -W               # call variants, with an index
vtools view -i 'QUAL>=30 && INFO/DP>=10' calls.vcf.gz
vtools filter -s LowQual -e 'QUAL<30' -Oz -o filtered.vcf.gz calls.vcf.gz
vtools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -r chr1:1-5000 calls.vcf.gz
vtools norm -f ref.fa -m -any -Ob -o norm.bcf -W calls.vcf.gz
vtools view --threads 4 -Ob -o calls.bcf calls.vcf.gz    # VCF to BCF
vtools stats calls.vcf.gz
```

Reading CRAM, ntools adds the `MD` and `NM` tags that CRAM writers leave out, as samtools does
(`--input-fmt-option decode_md=0` turns this off). The [Formats](#formats) section says where
CRAM references come from.

## Performance

The largest speedups, as wall time on a 24-vCPU Intel Xeon VM against samtools and bcftools 1.24,
with the same threads in both tools. The ntools rows were measured with ntools 0.2.0, which
compressed with libdeflate as samtools does, and the vtools rows with vtools 0.9.0, which
compresses with zflate, at a lower level by default ([Compression](#compression)):

| Tool | Command | Input | samtools/bcftools | ntools/vtools | Speedup |
|---|---|---|---:|---:|---:|
| ntools | `stats -@ 16` | HG002 PacBio Revio, whole genome (63 GB BAM) | 2,170 s | 294 s | 7.38x |
| ntools | `view -c -@ 16` | HG003 35x Illumina, whole genome (14 GB CRAM 3.1) | 23.5 s | 5.1 s | 4.64x |
| ntools | `flagstat -@ 16` | HG003 35x Illumina, whole genome (14 GB CRAM 3.1) | 30.4 s | 7.6 s | 4.02x |
| ntools | `stats -@ 16` | HG003 35x Illumina, whole genome (44 GB BAM) | 1,441 s | 554 s | 2.60x |
| ntools | `view` CRAM to BAM `-@ 16` | HG003 ONT R10.4.1, whole genome (69 GB CRAM) | 746 s | 445 s | 1.67x |
| vtools | `stats -s - --threads 4` | 1000 Genomes chr22, 500 samples | 16.6 s | 1.3 s | 13.19x |
| vtools | `view` VCF.gz to VCF `--threads 4` | 1000 Genomes chr22, 500 samples | 19.5 s | 1.6 s | 11.97x |
| vtools | `concat` 4 chunks `-Ob` | 1000 Genomes chr22, 500 samples | 8.3 s | 0.7 s | 11.62x |
| vtools | `view` VCF.gz to BCF `--threads 4` | 1000 Genomes chr22, 3,202 samples | 66.9 s | 8.4 s | 7.93x |
| vtools | `merge -Ob --threads 4` | 1000 Genomes chr22, 3,202 samples (sample halves) | 72.6 s | 29.8 s | 2.44x |

More rows, and the cases where samtools or bcftools is faster, are in the
[ntools](ntools/README.md#performance) and [vtools](vtools/README.md#performance) READMEs.

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

## Compression

BGZF blocks and CRAM's gzip blocks are compressed and decompressed, and their CRC32s computed, by
[zflate](https://github.com/DylanCosto/zflate), a Rust port of libdeflate. On an Intel Xeon Gold
5418Y:

- BAM (and compressed SAM) is written at level 5 unless a level is given, under 1% larger than by
  ntools 0.2.0, which wrote at libdeflate's level 6: `view -b` of long reads is 1.3 times as fast
  as in 0.2.0 single-threaded and 1.2 times at `-@ 8`, and `sort -@ 4` and `markdup` of short
  reads 1.5 and 1.7 times.
- VCF.gz and BCF, and `ntools bgzip`'s output, are written at level 6 (HTSlib's zlib-style 5), 5%
  larger than at HTSlib's default, libdeflate's level 7: `ntools bgzip` is 1.3 times as fast as
  HTSlib's single-threaded and 1.6 times at `-@ 8`, and vtools writes BCF and VCF.gz with
  threads 1.1–1.4 times as fast as 0.2.0. `-O z6`, `-O b6` and `bgzip -l 6` give HTSlib's size.
- Reading BAM is 7–13% faster than in 0.2.0, end to end (`view`, `flagstat`, `index`).

Given levels keep their meaning: they map onto libdeflate's scale, which zflate keeps, as HTSlib
maps them, and CRAM's gzip blocks keep HTSlib's levels. The compressed bytes are zflate's, not
HTSlib's, so `ntools bgzip` output and index files differ from HTSlib's byte for byte; they
decompress to the same data.

## Platforms

| Platform | Status |
|---|---|
| Linux x86_64 | Tested in CI; release binaries need glibc 2.28 or later |
| macOS aarch64 (Apple silicon) | Tested in CI; release binaries |
| Linux aarch64, other Unix systems | Not tested |
| Windows | Not supported (ntools uses Unix-only APIs) |

SIMD code (AVX2 for CRAM quality decoding and `vtools gtcheck`) is chosen at run time, with a
portable fallback; so is zflate's (BMI2, PCLMULQDQ, VPCLMULQDQ and AVX-512 on x86_64; the CRC and
PMULL instructions on arm64, where the multithreaded BGZF reader also decodes blocks two at a
time). Results do not depend on the platform: `vtools call -c` computes its
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
