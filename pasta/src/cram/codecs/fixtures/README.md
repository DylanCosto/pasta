# Codec fixtures

These are streams encoded by htscodecs 1.6.7 (as bundled with htslib 1.24), each
checked to decode to its input with the matching htscodecs decoder.

  * `r4x8-<input>-o<order>.bin`: `rans_compress(in, len, &out_len, order)`
    (rANS 4x8), with `order` 0 or 1.
  * `r4x16-<input>-<flags>.bin`: `rans_compress_4x16(in, len, &out_len, order)`
    (rANS Nx16), named by the flags byte written. The `order` arguments are the
    flags, except for `qual-64` (0x44; htscodecs adds CAT) and `qual-49`
    (`0x49 | 3 << 8`, i.e., STRIPE with 3 ways).

The inputs are

  * `qual`: `qual.bin`, 1027 bytes of quality score-like data (8 symbols, with
    runs);
  * `mid`: `sparse(8000, 1024, 8)`, for which htscodecs uses 12-bit order-1
    frequencies;
  * `big`: `sparse(70001, 2048, 4)`; and
  * `wide`: `sparse(4099, 16, 32)`.

`sparse` is defined in `../test_utils.rs`. The fixtures were written by a short
C program linked against htslib that builds the same inputs with the same LCG
(`Rng`).

`names.bin` holds 35 NUL-terminated Illumina HiSeq read names
(`HSQ1004:134:C0D8DACXX:…`). `tok3-rans.bin` and `tok3-arith.bin` are those names
encoded by htscodecs 1.6.7's `tok3_encode_names` at level 9, with rANS Nx16 and
with the arithmetic coder.

`names-formats.bin` holds 26 NUL-terminated read names of the formats htscodecs
recognizes (Illumina, Ion Torrent, Oxford Nanopore, PacBio) and others,
including duplicates, a name that begins an earlier one, padded numbers and an
empty name. `tok3-rans-formats-<level>.bin` are those names encoded by
`tok3_encode_names` at levels 1 and 3 (the level HTSlib uses for names at its
default level 5), with rANS Nx16.
