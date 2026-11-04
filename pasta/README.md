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
