# pasta

A Rust library for high-throughput sequencing data: SAM, BAM, and CRAM alignments (`sam`, `bam`,
`cram`), VCF and BCF variant data (`variant`, a port of HTSlib 1.24's `vcf.c` with its tabix and
CSI index building), BGZF compression (`bgzf`), CSI/BAI/CRAI/tabix indexes (`csi`, `bam::bai`,
`cram::crai`, `tabix`), FASTA references (`fasta`), positions and regions (`core`), a unified
alignment reader and writer (`util::alignment`), and an HTSlib-compatible pileup engine
(`pileup`).

It is one crate so that the compiler can inline across formats; both tools in this repository,
`ntools` and `vtools`, are built on it.

`variant` is the only variant API. The VCF/BCF modules that came from noodles (`vcf`, `bcf`, and
`util::variant`) silently corrupted genotype and INFO data on ordinary writes and were removed.

## Performance

Against noodles 0.119.0 (noodles-bam 0.98, -cram 0.102, -vcf 0.96, -bcf 0.94, -bgzf 0.52) on a
24-vCPU Intel Xeon Gold 5418Y VM, both built with fat LTO and one codegen unit into one benchmark
program. Its alignment code is the same source for both libraries; its variant code uses each
library's own API. Wall time, medians of three runs (one for runs over 30 s), pinned to one core;
the two `-@ 8` rows give CPU time, as their wall time varied with other load on the machine.
noodles was timed on 2026-11-08 and pasta 0.9.1 on 2026-11-20, with two noodles rows re-timed then
as a check (within 1%). `+ libdeflate` is noodles with noodles-bgzf's optional `libdeflate`
feature; without it, noodles uses zlib-rs. Every workload gives the same records and checksums
with both libraries.

Inputs: chr20 of GIAB HG003 Illumina NovaSeq PCR-free 35x (`ilmn`, 19.5 M records, 1.1 GB BAM)
and of HG002 PacBio Revio HiFi (`hifi`, 221,258 records, 1.5 GB BAM); 2 M HG003 Illumina reads on
chr20 (`2M`, 109 MB BAM), also as CRAM 3.0 and 3.1 written by samtools 1.24; and the first
200,000 records of 1000 Genomes chr22 with 3,202 samples, as VCF.gz (79 MB) and BCF (59 MB).

| Workload | noodles | + libdeflate | pasta | Speedup | vs + libdeflate |
|---|---:|---:|---:|---:|---:|
| Read BAM, fields on access (`ilmn`) | 10.31 s | 6.67 s | 6.52 s | 1.58x | 1.02x |
| Read BAM, fields on access (`hifi`) | 14.77 s | 9.62 s | 8.64 s | 1.71x | 1.11x |
| Read BAM, `-@ 8` (`ilmn`), CPU | 12.00 s | 8.42 s | 7.50 s | 1.60x | 1.12x |
| Read BAM, every field decoded (`ilmn`) | 39.22 s | | 35.14 s | 1.12x | |
| Read BAM, every field decoded (`hifi`) | 34.10 s | | 27.74 s | 1.23x | |
| BAM to SAM (`2M`) | 2.07 s | | 1.56 s | 1.33x | |
| BAM to BAM (`2M`) | 11.35 s | 8.33 s | 5.25 s | 2.16x | 1.59x |
| BAM to BAM, `-@ 8` (`2M`), CPU | 12.53 s | 9.49 s | 5.96 s | 2.10x | 1.59x |
| Read CRAM 3.1 (`2M`) | 11.16 s | | 3.66 s | 3.05x | |
| Read CRAM 3.0 (`2M`) | 6.89 s | | 3.71 s | 1.86x | |
| BAM to CRAM 3.0 (`2M`) | 15.36 s | | 6.47 s | 2.37x | |
| Read BCF, positions and alleles | 0.82 s | | 0.59 s | 1.39x | |
| Read BCF, all genotypes | 20.28 s | | 3.69 s | 5.50x | |
| Read VCF.gz, positions and alleles | 11.69 s | 10.57 s | 3.54 s | 3.30x | 2.99x |
| Read VCF.gz, all genotypes | 58.63 s | | 6.57 s | 8.92x | |
| VCF.gz to BCF | 124.10 s | | 9.06 s | 13.70x | |
| BCF to VCF | 81.35 s | | 5.58 s | 14.58x | |

- Output sizes: BAM 109.0 MB (noodles, level 6), 110.5 MB (+ libdeflate, level 6), and 114.4 MB
  (pasta, level 5); CRAM 3.0 53.4 MB (noodles) and 38.7 MB (pasta; samtools' is 38.6 MB); BCF
  57.2 MB (noodles) and 59.9 MB (pasta).
- noodles writes only CRAM 3.0. pasta's default, CRAM 3.1, takes 6.82 s and gives 33.8 MB
  (samtools': 33.7 MB).
- Memory: under 20 MB in both libraries for every BAM, VCF, and BCF workload. Reading CRAM takes
  27 MB in pasta against 131 MB (3.0) and 159 MB (3.1) in noodles; writing CRAM 3.0, 98 MB
  against 127 MB.
- pasta's VCF reader parses every field of each line, genotypes included, as HTSlib does; noodles'
  splits the line and decodes fields when they are read. The "positions and alleles" rows
  therefore do more work in pasta.

## Changes from noodles

pasta started from the noodles crates that ntools used (`noodles-core`, `-bgzf`, `-csi`, `-fasta`,
`-sam`, `-bam`, `-cram`, `-util`, and the new `-pileup`) as of the `DylanCosto/noodles` fork at
commit `0b0f56e01` (upstream zaeleus/noodles as of 2026-10-15 plus the changes below), merged
into one crate without the async API and the VCF/BCF utilities. Changes made before the move:

- CRAM: correctness fixes, table-driven rANS decoding, direct tok3 name decoding, and
  reference-range retrieval (submitted upstream as zaeleus/noodles#440–#443).
- CRAM: unmapped records decode with a mapping quality of 0, as in HTSlib.
- CRAM: sequences are reconstructed in bulk (`Record::with_sequence`), and tag blocks are found
  without hashing.
- CRAM: `Record::encode_bam_record` encodes a decoded record as raw BAM directly, so CRAM input
  takes the same raw BAM paths as BAM input.
- CRAM: 32-state order-1 rANS Nx16 streams (quality scores) decode with AVX2 when the CPU has
  it, falling back to the scalar decoder otherwise.
- BAM: `Record` exposes its encoded bytes (`AsRef<[u8]>`, `TryFrom<Vec<u8>>`,
  `From<Record> for Vec<u8>`) for copying records without re-encoding.
- SAM: integer arrays are formatted through a buffer, which speeds up long-read tags (`MM`/`ML`).
- CRAM: decoding regenerates `MD`/`NM` tags that HTSlib's writer drops (as HTSlib's reader
  does, unless turned off with `io::reader::MdNm::AsStored`, HTSlib's `decode_md=0`), and fills
  bases past the end of the reference with `N` instead of panicking.
- CRAM: the writer stores records on sequences missing from the reference against a consensus
  embedded in the slice, or verbatim in slices on many sequences (their containers are not
  "reference required"), instead of failing, and `@SQ` lines may lack `M5`.
  Records in multi-reference slices that do not require the reference no longer get one, like
  single-reference slices. FASTA: the indexed reader adapter returns `None` for a name missing
  from the index, like the other adapters.
- CRAM: a reference sequence missing from the repository is an error rather than a panic, and
  `Reader::set_reference_sequence_repository` lets a reference be chosen after reading the
  header. (A corrupt container whose records need reference bases that their slice does not
  provide can still panic when their sequences are read.)
- CRAM: `io::reader::read_block`, `Block`, and `Container::src` expose a container's raw blocks
  (method as stored, content type and ID, sizes, compressed bytes) for block-level tools such
  as `cram-size`.
- BAM: a record whose CIGAR is held in a `CG` field (over 65,535 operations) no longer lists
  that field in its data, as in HTSlib; `Record::cigar` already returned its operations.
- `pileup` (new): HTSlib-compatible pileup (`bam_plp`), BAQ, mapping quality capping,
  MD/NM calculation, and MM/ML base modifications over encoded BAM records, which `mpileup` and
  `calmd` use. Positions inside a CIGAR operation are resolved from a cached
  cursor without reading the CIGAR, the per-position scan reads a compact array rather than the
  buffered records, and overlapping mates' qualities are combined in lockstep along shared match
  operations.
- CRAM writer: each block's compression method is chosen as HTSlib 1.24 does
  (`io::writer::compression`: CRAM 3.0/3.1, levels, the `fast`/`normal`/`small`/`archive`
  profiles), and slices also end at a base limit. The default is CRAM 3.1 with the normal profile.
- CRAM codecs: the name tokenizer encoder splits names into letters, digits, and punctuation and
  encodes numbers as deltas from the previous name as the decoder resolves it (about 0.03
  bytes per sequential name, from 3.6); rANS Nx16 encodes with precomputed reciprocals and
  only touches the order-1 contexts present; frequency scaling no longer overflows 32 bits for
  symbols seen more than 2^20 times in a block; the arithmetic coder accepts empty streams.
- CRAM reader: read names stored as `*` stay missing; names are generated only for slices that
  store none. Indexing counts records with a reference but no position from position 0, as
  HTSlib does.
