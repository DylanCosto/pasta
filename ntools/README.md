# ntools

A samtools-compatible command-line toolkit built on the [pasta](../pasta) library. It aims to
produce the same results as samtools 1.24 while being faster, especially with threads.

## Install

```sh
cargo install --path ntools
```

## Commands

| Command | samtools equivalent | Notes |
|---|---|---|
| `view` | `view` | SAM/BAM/CRAM conversion; `-c`, `-f/-F/-G`, `--rf`, `-q`, `-m`, `-e` (filter expressions), `-L`, `--region-file`, `-M`, `-P`, `-s`, `--subsample/--subsample-seed` (with samtools' `@CO` line and header-derived `auto` seed), `-N`, `-r/-R/-n`, `-d/-D`, `-x` (and `^TAGS`), `--keep-tag`, `--add-flags/--remove-flags`, `-U`, `-p`, `-B`, `-t`, `-X` (index after the input), `-l` (library), `-h/-H`, `--no-header`, `--save-counts`, `-b/-C/-u/-1`, `-O FMT,OPT,...` and `--output-fmt-option` (CRAM version, level, profiles, codecs), `--write-index`, `-?`, `-S`, regions and `*` |
| `sort` | `sort` | coordinate, `-n` (natural), `-N` (bytewise), `-t TAG` (then by position, or by name with `-n`/`-N`), `--template-coordinate` (needs `MC` tags from `fixmate -m`; the last of `-n`/`-N`/`--template-coordinate` wins and `-t` overrides it, as in samtools); `-m`, `-T`, `-@`; `--write-index` (coordinate order only; `-o FILE##idx##INDEX`); stable, like samtools. No minimiser order (`-M`, `-I`, `-K`, `-w`); `-H` and `-R`, which only modify it, are accepted and ignored |
| `merge` | `merge` | SAM, BAM, and CRAM inputs; headers merged as samtools does, renaming colliding `@RG`/`@PG` IDs (`-s` seed) or combining them (`-c`, `-p`), with `RG`/`PG` tags translated; `-n/-N`, `-t`, `--template-coordinate`, `-r`, `-R`/`-L` (indexed inputs, `-X` index files), `-b`, `-h`, `-1/-l/-u`, `-f`, `-o`, `--write-index` (output `FILE##idx##INDEX` names the index) |
| `index` | `index` | BAI (default) or CSI (`-c`, `-m` minimum shift, with HTSlib's level count) for BAM; CRAI for CRAM; `-M` (several inputs), `-o` or a second argument naming the index; like samtools, checks the records' order rather than the header's `SO` |
| `faidx`, `fqidx` | `faidx`, `fqidx` | FASTA and FASTQ, plain or BGZF (`.fai` and `.gzi` as HTSlib builds them); regions as HTSlib parses them (`{name}` quoting, `1,000`, `2k`), `-r`, `-o` (BGZF for `.gz`), `-n`, `-c`, `-i`, `--mark-strand`, `--fai-idx/--gzi-idx`, `--write-index`, `-f` |
| `flagstat` | `flagstat` | default, `-O tsv`, `-O json` |
| `idxstats` | `idxstats` | from BAM index metadata (`-X` names the index); scans CRAM/SAM |
| `depth` | `depth` | `-a/-aa`, `-r`, `-b` (with `-r`, the positions in both), `-J`, `-s`, `-q/-Q`, `-g/-G`, `--incl-flags`, `--require-flags`, `-l`, `-H`, `-o`, several inputs, `-f` (input list), `-X` (index files after the inputs); `-d/-m` are ignored, as in samtools |
| `coverage` | `coverage` | per reference or `-r` region; `-q/-Q`, `-l`, `--rf/--incl-flags`, `--ff/--excl-flags`, `-d` (the pileup's limit, applied as HTSlib does), `--min-depth`, `-m/-D/-A/-w` histograms (sized to the terminal like samtools), `-H`, `-o`, several inputs or `-b` list |
| `fastq`, `fasta` | `fastq`, `fasta` | `-0/-1/-2/-s/-o`, `-n/-N`, `-t/-T` (`*` for all tags), `-F/-f/-G/--rf`, `-d/-D` (tag filters), `-O` (OQ qualities), `-i` (CASAVA), `--i1/--i2`, `--index-format`, `--barcode-tag/--quality-tag`, `-v`, `-U/--UMI-tag`, `--no-sc/--sc-aux/--no-sc-bkp`, `-c`, `.gz` output |
| `quickcheck` | `quickcheck` | `-v`, `-q`, `-u`; stdin, pipes and FIFOs (whose EOF block, as in HTSlib, cannot be checked) |
| `reheader` | `reheader` | BAM and CRAM (2.1 and 3.x); `-P`, `-c CMD`, `-i` (CRAM in place); copies compressed blocks and containers, checking CRAM container header CRC32s; reads stdin, pipes and FIFOs |
| `cat` | `cat` | BAM and CRAM (2.1 and 3.x); `-b`, `-h`, `-o`, later inputs' `@RG` lines added; CRAM `-r` (region, `*`, or `#:start-end` containers), `-p N/M`, `-f`, `-q`; copies compressed blocks and containers, checking CRAM container header CRC32s and re-encoding only CRAM containers that must be filtered or get new read group numbers; reads stdin, pipes and FIFOs (`-r` and `-p` need an indexed file) |
| `split` | `split` | by read group or `-d` tag; `-f`, `-p/--zero-pad`, `-u`, `-M`, `--write-index` |
| `collate` | `collate` | groups templates by name; `-f`, `-O/-o`, `-u`, `-l`, `-m`, `-T`; `-n` and `-r` are accepted and ignored (they size samtools' temporary files) |
| `fixmate` | `fixmate` | `-m`, `-c`, `-r`, `-p`, `-z`, `-M` (base modification tags); parallel with `-@` |
| `dict` | `dict` | plain, gzip, or BGZF FASTA; `-a`, `-s`, `-u`, `-H`, `-o`, `-A` (`AN` aliases), `-l` (`AH` for alternate loci) |
| `bedcov` | `bedcov` | `-c`, `-d`, `-j`, `-Q`, `-g/-G`, `-H`, `--max-depth` (the pileup's limit, applied as HTSlib does); several inputs; `-X` (index files after the inputs) |
| `addreplacerg` | `addreplacerg` | `-r` (repeatable, `\t` escapes), `-R`, `-m overwrite_all\|orphan_only`, `-w`, `--write-index`, `--input-fmt` |
| `markdup` | `markdup` | `-r`, `-l`, `-m t\|s`, `-c`, `-S`, `-t`, `-d` (optical duplicates), `-T`, `-s/-f`, `--json`, `--no-multi-dup`, `--include-fails`, `--use-read-groups`, `--duplicate-count`, `--read-coords`/`--coords-order`, `--barcode-tag`/`-name`/`-rgx`, `--move-umi-to-tag`, `--umi-separator`, `--write-index` |
| `stats` | `stats` | all sections; `-r`, `-t`, regions, `-f/-F/-d`, `-c`, `-i`, `-m`, `-q`, `-l`, `-I`, `-S/-P`, `-p`, `-g`, `-x`, `--GC-depth`, `--ref-stats`; per-cycle counts use `-@` threads |
| `mpileup` | `mpileup` | text pileup with BAQ (`-B`, `-E`), overlap removal (`-x`), `-f` (or `--reference`), `-r`, `-l`, `-a/-aa`, `-q/-Q`, `-d`, `-A`, `-C`, `-G`, `-R`, `-6`, `--rf/--ff`, `-b`, `-X`, `-s`, `-O`, `--output-BP-5`, `--output-QNAME`, `--output-extra`, `--output-sep/-empty`, `--reverse-del`, `--no-output-ins/-del/-ends`, `-M` (base modifications from MM/ML, also in insertions unless `--no-output-ins-mods`); reads and BAQ use `-@` threads, and indexed regular files with several sequences pile up one chromosome per thread (stdin, pipes and FIFOs are read through) |
| `calmd` | `calmd` | `-e`, `-r` (BAQ), `-A`, `-E`, `-n`, `-C`, `-d`, `-q`, `-N`, `-Q`, `-b/-u`, `--output-fmt`; records are processed on `-@` threads |
| `ampliconclip` | `ampliconclip` | soft or `--hard-clip`, `--both-ends`, `--strand`, `--tolerance`, `--original`, `--keep-tag`, `--clipped`, `--fail`, `--filter-len/--fail-len/--unmap-len`, `--no-excluded`, `--rejects-file`, `--primer-counts`, `-f` |
| `import` | `import` | FASTQ/FASTA (plain or gzipped) to unaligned SAM/BAM/CRAM; `-0/-1/-2/-s`, `--i1/--i2` (index reads to `BC`/`QT`), `-b`, `-i` (CASAVA), `-U` (UMIs), `-N` (SRA names), `-T` (tags from comments, such as `MM`/`ML`), `-r/-R`, `--order`, `-O`, `-u`, `-o` |
| `reset` | `reset` | `-x`/`--keep-tag` (including `^` lists), `--reject-PG`, `--no-RG`, `--dupflag`, `-O`, `-o` |
| `checksum` | `checksum` | `-F/-f/-b`, `-c`, `-t` (including `*` negation), `-O`, `-P/-C/-M`, `-z`, `-N`, `-q`, `-v`, `-a`, `-T`, `-B`, `-m` (merge reports), `-o` |
| `head` | `head` | `-h` (header lines), `-n` (records) |
| `flags` | `flags` | numbers and flag names, both ways |
| `samples` | `samples` | `-h`, `-i`, `-T`, `-f/-F` (matching references), `-X`, `-o`, paths from stdin |
| `consensus` | `consensus` | FASTA/FASTQ/pileup; `simple` and Bayesian modes (`bayesian`, `_r`, `_p`, `_m`, `_116`), `-X` presets and `-t` calibration files, `-a/-aa`, `--show-ins/--show-del/--mark-ins`, `-A`, `-r`, `--regions-file`, `-T`, every tuning option; `-@` on an indexed input computes `-Z` blocks in parallel, with the same output as samtools run the same way (which, like samtools, can differ from a serial run at block edges with `--mark-ins` or `-A`) |
| `cram-size` | `cram-size` | sizes by content ID and expanded compression method (gzip level, rANS order and transforms, ...), with the data series using each block; `-v`, `-e` (encodings), `-o` |
| `reference` | `reference` | from `MD` tags (any input) or embedded CRAM references (`-e`); `-r` (indexed), `-q`, `-o` |
| `rmdup` | `rmdup` | paired-end, `-s`, `-S`; libraries from `@RG` `LB`; single-end mode holds only the reads that can still be replaced |
| `bgzip` | HTSlib's `bgzip` | compress (`-c`, `-k`, `-f`, `-o`, `-l`, `-@`, `--binary`), `-d`, `-t`, `-b`/`-s`, `-i`/`-I`, `-r`, `-g`, several files; blocks cut where HTSlib cuts them (65,280 bytes, ending at line ends for text, with header lines in blocks of their own), so `.gzi` indexes list the same uncompressed offsets; `-l` levels mapped as HTSlib maps them (without one, level 6; the bytes are zflate's, not HTSlib's); output names, overwrite prompts, timestamps, and exit statuses as in HTSlib |

Region arguments (`view`, `depth -r`, `mpileup -r`, `coverage -r`, `stats`, `merge -R`, and the
rest) are parsed as HTSlib's `hts_parse_region` parses them: `chr1:1,000-2,000`, `chr1:1k-2k`,
`chr1:100` and `chr1:100-` (to the end), `chr1:-200`, and `{name}` quoting for names containing
colons, which also resolve unquoted when that is unambiguous. Interval files (`view -L`,
`depth -b`, `mpileup -l`, `merge -L`, `consensus --regions-file`) are read as samtools'
`bed_read` reads them: plain, gzip, or BGZF; columns separated by any whitespace; two-column
`chr pos` lines as 1-based positions; and `#`, `track`, and `browser` lines skipped.

CRAM input without `--reference` reads each reference sequence from the FASTA file its `@SQ`
line's `UR` field names (a path or a `file:` URL), when that file exists locally with a `.fai`
index, as HTSlib does; a sequence whose line names no usable file is looked up in the other named
files. A sequence found nowhere fails, naming the file tried and why it did not serve (`ntools:
x.cram: no reference for chrS (tried UR:/path/ref.fa, not found); give one with -T/--reference`),
where samtools prints HTSlib's warnings and a link to its reference documentation. ntools does
not search `REF_PATH` or `REF_CACHE`. CRAM output without `--reference` takes the same files;
records on sequences none has are stored as with a reference that lacks them (below), with one
warning, as in samtools (75 KB against samtools' 74 KB for 49,484 short reads).

Commands take samtools' shared long options where samtools does: `--input-fmt-option`
(`decode_md=0|1`, `reference=FILE`; `nthreads` is ignored), `--output-fmt-option` (`level=N`,
the CRAM options below; `nthreads` is ignored), and `--write-index` (on `view`, which indexes its `-U` output too, `split`,
`addreplacerg`, `markdup`, and `faidx`; accepted and ignored elsewhere, as in samtools). As in
HTSlib, an alignment output named `FILE##idx##INDEX` is written to `FILE` and indexed to `INDEX`
(BAI if it ends in `.bai`, otherwise CSI for BAM). Other format options are errors rather
than silently ignored. `--verbosity` is accepted and ignored.

`-@ N` adds worker threads throughout. Beyond BGZF compression and decompression, which is what
samtools mostly uses its threads for, ntools also formats SAM text, decodes CRAM containers
(region queries a slice per task, so even one region uses every thread), encodes them, sorts,
and fixes mates in parallel.

A typical duplicate-marking pipeline is unchanged apart from the program name:

```sh
ntools collate -Ou -@ 8 in.bam | ntools fixmate -m -@ 8 - - | ntools sort -@ 8 - \
  | ntools markdup -@ 8 - out.bam
```

On the short-read corpus this pipeline takes 1.30 s, against 3.24 s for the same samtools
pipeline, and marks the same 26,749 duplicates in identical records.

## Compatibility

ntools has samtools 1.24's commands, options, and output formats, and HTSlib's `bgzip`. Where
samtools or HTSlib gives a wrong result (a crash, wrong output, or output another tool rejects),
ntools gives the right one instead of copying it; those cases are below, with the other
intentional differences.

Intentional differences:

- BGZF output (BAM, compressed SAM) is written at level 5 unless a level is given, and `bgzip`'s
  at 6: faster than samtools' default, and a little larger
  ([Compression](../README.md#compression)). The compressed bytes are zflate's, not HTSlib's;
  they decompress to the same data.
- `@PG` lines name ntools. Use `--no-PG` (or `-P` for `reheader`) for identical headers.
- SAM text is read as HTSlib's `sam_parse1` reads it: an RNAME or RNEXT missing from the header
  (by `SN` or `AN`) is read as `*`, a record with no reference sequence, POS 0, or CIGAR `*`
  gets the unmapped flag, and RNEXT with PNEXT 0 is read as `*`, each with HTSlib's warning.
  Where HTSlib's result is wrong, ntools differs:
  - The position of a reference sequence dropped because it is missing from the header goes
    with it: `x 0 chrX 10` reads as `x 4 * 0`, and RNEXT `chrY` PNEXT 40 as `* 0`. HTSlib keeps
    POS and PNEXT (`x 4 * 10`), placing a read at a position on no reference sequence, which
    htsjdk rejects ("Alignment start should be 0 because reference name = *").
  - RNEXT `=` is the RNAME as written: `g 1 chr1 0 … = 40` keeps its mate at `chr1:40`. HTSlib
    drops RNAME for POS 0 first, so `=` becomes `*` (with PNEXT 40 kept).
  - HTSlib prints a warning for every record; ntools prints each distinct one once per input.
- A write to a closed pipe (`| head`) kills ntools with SIGPIPE, as it kills samtools, but only
  once the command has removed its temporary files; samtools leaves `sort` and `collate` spill
  files behind.
- `view` on a CRAM `*` region returns every unplaced unmapped read; samtools 1.24 returns only
  some of them from the corpus CRAM. ntools matches samtools' output for the equivalent BAM.
- `view -e` evaluates HTSlib's filter expression language with HTSlib's semantics (undefined
  values, tag presence as truth, no short-circuiting), with these fixes: an undefined value
  `%` a number stays undefined (HTSlib turns `[XX] % 2` into 0 when `XX` is absent, so
  `[XX] % 2 == 0` selects every record without the tag), and `%` by a value between -1 and 1 is
  undefined (HTSlib crashes). A syntax error is reported before any record is read; samtools
  reports it at the first record, and not at all for an empty input. A regular expression
  computed from record fields is compiled for each record; HTSlib reuses the first pattern it
  saw.
- `view -e` regular expressions are POSIX extended expressions run by the Rust `regex` crate on
  bytes (matching for `=~`/`!~` only needs a yes or no, so leftmost-longest rules do not
  matter). Back-references (`\1`) are not supported. As in glibc, `\w`, `\W`, `\s`, `\S`,
  `\b`, `\B`, `\<`, `\>` keep their GNU meaning, a backslash is literal inside brackets, and a
  backslash before another letter or digit is that character (`\d` is `d`, not a digit).
- `view -c -n` on CRAM counts the records with a read group, as on BAM. samtools 1.24 requests
  no CRAM fields for `-n` alone, so its records have no read group and it counts 0 (unless
  another option, such as `-r`, `-R`, `-l`, or `-d`, decodes read groups).
- `view -R ^FILE` removes the `@RG` lines of the excluded read groups. samtools instead keeps
  only those lines and removes the lines of the read groups whose records it outputs.
- `view -P` applies `-p` and `-U` to the records it fetches (the reads in the regions and their
  mates); samtools ignores both options with `-P`.
- `view -B` on a record without base qualities leaves them all missing. HTSlib marks only the
  first one missing and keeps merged values in the others (SAM output shows `*` either way).
  pasta cannot parse `B` operations in SAM text, so `-B` needs BAM input.
- `view -T ref.fa` on SAM without `@SQ` lines uses `ref.fa.fai` only if it exists; HTSlib
  builds it.
- `merge` renames clashing `@RG`/`@PG` IDs with the same random suffixes as samtools for a
  given `-s` seed; without `-s`, both seed from the clock, so suffixes differ between runs.
- `merge` refuses to overwrite an existing output given as `FILE##idx##INDEX` without `-f`;
  samtools checks for a file with the whole name, so it overwrites.
- `merge --write-index` indexes the output after writing it, and checks while writing that it is
  sorted by coordinate (samtools fails the same way, from its on-the-fly indexer).
- `view -L` and `merge -L` ignore zero-length BED intervals (`chr1 1000 1000`). samtools'
  `bed_overlap` treats one as overlapping the reads that span the point between the two bases.
- `sort --template-coordinate` and `merge --template-coordinate` compare reference sequence IDs
  directly; samtools subtracts them, which overflows for a record flagged as mapped without a
  reference sequence. With `-N --template-coordinate -n`, ntools sorts by name naturally;
  samtools sorts bytewise, since `-N` anywhere switches it to bytewise names.
- BAI/CSI bins are written in insertion order, not HTSlib's hash order. samtools reads them and
  returns identical query results.
- Writing CRAM, ntools leaves out the `MD` and `NM` tags that readers regenerate with the same
  values, as HTSlib does, and keeps them for alignments past the end of the reference or with a
  read `N` on a reference `N`. Reading CRAM, ntools regenerates `MD`/`NM`
  as HTSlib does: for mapped records with a reference, missing tags are added after the stored
  ones, `MD` before `NM` (a `cF` tag can turn this off, and is removed). Like samtools, `fastq`,
  `mpileup`, `depth`, and `coverage` do not regenerate them, so `fastq -T` and
  `mpileup --output-extra` show only stored tags.
- Reading CRAM records that overhang the end of the reference, ntools fills the bases past the
  end with `N`, as HTSlib does when it can decode them; HTSlib fails on such slices when it
  must regenerate `MD`.
- CRAM output is CRAM 3.1 by default, as in samtools, with each block's compression method
  chosen the way HTSlib chooses it (trials of the allowed methods, then the best on average,
  re-tried periodically). `-O cram,OPTS` and `--output-fmt-option` take `version=3.0|3.1`,
  `level=N`, the `fast`/`normal`/`small`/`archive` profiles, `use_bzip2`, `use_lzma`,
  `use_rans`, `use_tok`, `use_fqz`, `use_arith`, `seqs_per_slice`, `bases_per_slice`,
  `store_md` and `store_nm` (keep the `MD`/`NM` tags readers regenerate), and `ignore_md5` (see
  below). Files are within a few percent of samtools' size (the short-read corpus: 0.87 MB
  against 0.86 MB by default, 0.58 MB against 0.65 MB with `archive`). Gzip blocks use the default strategy, not
  HTSlib's `Z_FILTERED`/`Z_RLE`; fqzcomp uses one parameter set; and read names use the rANS
  name tokenizer, never its arithmetic-coded form. Containers are encoded and compressed on the
  `-@` threads, with the same output for any thread count. Single-threaded, writing is slower
  than samtools (2.6 s against 2.0 s for the short-read corpus); with `-@ 16` it is faster
  (0.35 s against 0.49 s). ntools rejects HTSlib's other CRAM options rather than ignore them:
  `embed_ref` and `no_ref` (sequences missing from the reference are handled as below),
  `lossy_names`, output `reference=` (use `--reference`), and on input `ignore_md5` and
  `required_fields` (with which samtools prints `*` for the fields left out).
- CRAM read names stored as `*` stay `*`. Files written without read names (`lossy_names`) get
  names from the record counter (`1`, `2`, ...); HTSlib prefixes the file name (`x.cram:2`).
- `view -1` and `-u` write BAM unless `-b`, `-C`, `-O`, or the output name chooses another format;
  `--output-fmt-option level=N` overrides them, as in HTSlib.
- CRAM output fails, as in samtools, on records on a sequence whose `@SQ` `M5` is not the
  checksum of that sequence in `--reference` (e.g., reads aligned to UCSC's hg38 analysis set,
  written with the GCA no-alt FASTA): such a file would not decode with the reference its
  records are stored against. `@SQ` `LN` fields that differ from the reference get HTSlib's
  warning. With `--output-fmt-option ignore_md5=1`, the records are written, and the `M5` fields
  of the sequences the reference has are replaced with its checksums, so that ntools and
  samtools read the file back with that reference; samtools keeps the input's `M5` fields, which
  then describe another reference (and suggests `embed_ref=2`, which ntools does not take).
- CRAM input checks the reference for slices without a reference checksum (slices with many
  reference sequences, which samtools writes for long reads); HTSlib decodes them against
  whatever the reference holds. Where the reference has a sequence dictionary beside it
  (`ref.dict` or `ref.fa.dict`, from `ntools dict`, `samtools dict`, or Picard), its `M5` must
  match the header's, which also catches a reference of the right lengths with other bases (such
  as UCSC's hg38 analysis set for GRCh38). Without one, the length must match `@SQ` `LN` (a
  reference shorter than `LN` is an error); `ntools dict ref.fa -o ref.dict` makes one.
- CRAM output with a reference that lacks some `@SQ` sequences warns as HTSlib does. Containers
  holding records on those sequences do not require the reference: slices on one reference
  sequence embed a consensus of their records to store the bases against, as HTSlib's
  `embed_ref=2` does, and other slices store the bases verbatim. HTSlib fills its consensus from
  `MD` tags where reads have them (so it can drop matching `MD`/`NM`), and for unsorted input
  drops the reference for the whole file (52 MB against 1.0 MB for 200,000 short reads, half on a
  missing sequence). Records against a consensus decode without regenerated `MD`/`NM`, as in
  HTSlib.
- `reheader` and `cat` recompress the block holding the end of the header with pasta's encoder,
  so files can differ in size by a few bytes; records are identical.
- `reheader` and `cat` on CRAM fail (exit 1) at a data container whose header fails its CRC32
  check, naming its offset. samtools 1.24 reports the failure, then stops as if the file had
  ended and exits 0, leaving the output without the rest of the records. As in HTSlib, the
  copied blocks' own CRC32s are not checked, since they are not decompressed.
- `reheader` and `cat` on BAM check the BGZF framing of every block they copy and fail (exit 1)
  on a BAM cut inside a block (an interrupted download, say). samtools 1.24 warns, copies the
  partial block and exits 0, so its output passes `quickcheck` but cannot be decoded. A BAM that
  ends between blocks without an EOF marker is copied, with samtools' "EOF marker is absent"
  warning.
- `cat` on CRAM re-encodes the containers it filters to a region (`-r` without `-f`) and those
  whose single read group gets a new number with pasta's encoder (which needs the reference,
  from `--reference` or the `@SQ` `UR` fields, and does not write CRAM 2.1); HTSlib edits the
  compression header for the latter. Decoded records are the same, and `MD`/`NM` tags are dropped
  where HTSlib would regenerate them. `cat -p N/M` without `-r` divides all of a file's
  containers; samtools 1.24 divides a single container, so part 1 copies the whole file and the
  others nothing. Containers are numbered (`-q`, `-p`, `#:`) by the index's entries for each
  reference in turn, then unplaced ones, which is HTSlib's order unless slices nest.
- `collate` groups each template's records together, READ1 before READ2, as samtools does, but
  templates come out in a different (hash) order, and records with the same READ1/READ2 flags
  (a primary and its secondary or supplementary alignments, or reads with neither flag) keep
  their input order, which samtools' unstable sort does not. `collate -f` keeps every primary
  record, as its usage says ("only primary alignments"); samtools 1.24 also drops primaries
  without exactly one of READ1 and READ2, so single-end input comes out empty. ntools writes the
  records of a template together even when they do not pair up, where samtools writes the
  leftovers at the end. Output named from a prefix in BGZF SAM (`--output-fmt sam.gz`, and
  `split`'s `%.`) ends in `.sam.gz`; samtools names it `.sam`, which hides its compression.
- `ampliconclip` fixes samtools bugs: `OA` tags follow the SAM specification
  (`RNAME,POS,strand,CIGAR,MAPQ,NM;`), where samtools writes the read name as RNAME and drops
  the comma before MAPQ without an `NM` tag; clipping the 3' end of a read that ends in a
  deletion keeps its other CIGAR operations, which samtools drops (leaving, say, `1D` for 18
  bases, so the read is then unmapped); and clipped records get a recalculated BAM bin.
- `import` parses `-T` tags that follow a CASAVA barcode (`-i`), which samtools loses. It reads
  an empty input file as having no records; samtools fails on it ("truncated
  file"), as it does on the empty `-0` file `samtools fastq` often writes. An `-r` read group
  line without an ID, or with a malformed field, is an error; samtools writes nothing and exits
  successfully.
- `markdup`: with `-d`, duplicates whose names have no coordinates are not compared with each
  other; samtools compares them using uninitialized memory, so a chain mixing parseable and
  unparseable names can come out differently. `--read-coords` and `--barcode-rgx` take Rust
  `regex` syntax rather than POSIX extended regular expressions (the same for ordinary
  patterns). Without `-S`, `--move-umi-to-tag` moves the UMI of every read; samtools 1.24 moves
  only those of reads still buffered at the end of the input.
- `consensus` fixes samtools bugs instead of reproducing them. Pileup output with `-a` lists each
  column once; samtools repeats growing runs of depth-0 lines after positions whose consensus is
  a deletion (and prints such positions as depth 0). `-@` gives the serial output, and on an
  unindexed input writes pileup lines, which samtools drops. A read's qualities never come from
  other reads (samtools takes a deleted or padded position's quality from stale bytes of earlier
  records), and qualities above 100 use the tables' quality-100 values instead of reading past
  them. With `-@` and a `-T` FASTA that lacks a needed sequence, ntools stops with an error where
  samtools writes truncated or garbage sequence. Positions past the end of a `-T` sequence (reads
  placed beyond `@SQ LN`) are filled as without `-T`, with `N` and quality 0; samtools reads past
  the end of its copy of the reference there and prints whatever bytes follow it, which can end
  a sequence line early and leave FASTQ sequence and quality of different lengths.
- `coverage` averages base qualities over the bases that have one: a read with `QUAL` `*`
  (stored as `0xff`) is left out of `meanbaseq`, as one without a sequence is (with the same
  warning), where samtools counts its bases as quality 255. samtools also reports a `meanbaseq` (column 8) of
  255 for CRAM input unless `-Q` is above 0, as it does not decode CRAM base qualities otherwise;
  ntools reports the mean base quality, as for BAM. Histogram bins count in 64 bits, so `-D` cannot overflow.
- `depth -X` with an odd number of files is an error (samtools' check never fires, and it
  pairs the files up wrongly), and `-H` labels the columns with the input names; samtools
  prints the index file names.
- `view` applies `--add-flags`/`--remove-flags` with `-P` too; samtools skips them there. An
  out-of-range `--subsample` is an error with a nonzero exit status (samtools prints the error but
  exits 0).
- `view -c -o -` prints the count to stdout (samtools writes it to a file named `-`), and
  `view -c` ignores `--write-index` (samtools 1.24 crashes). With `--write-index`, `view` checks
  that `-o` and any `-U` output are BAM or CRAM files before writing anything; samtools writes
  the outputs and then fails.
- With a compression level in `--output-fmt-option`, samtools writes SAM output BGZF-compressed;
  ntools writes plain SAM (the same text).
- `bedcov -d` adds its column without changing the depth sum, which counts deletions and
  reference skips unless `-j`, as documented; samtools also drops them from the sum when `-d`
  is given. The `-d` column counts aligned bases, as in samtools.
- BAQ (`mpileup`, `calmd -r`) keeps the quality of a base it is certain of, such as in a long
  read that matches the reference for thousands of bases. samtools on x86-64 converts that
  base's infinite BAQ phred value to `INT_MIN` and sets its quality to 0; on arm64 it keeps the
  quality, as ntools does.
- `mpileup` stops with an error when an input's reference sequences differ from the first
  input's (positions are compared by reference index). samtools mixes up such inputs, or fails
  partway when one has more sequences.
- `depth`, `coverage`, and `cat` with several inputs stop with an error when an input's
  reference sequences (names and lengths, in order) are not the first input's, or the first of
  them (an input with fewer sequences, at the same indexes, is fine). samtools
  labels every input's reads with the first input's sequence names, so it reports reads on the
  wrong sequence, drops them, crashes, or (`cat`) writes a BAM whose records name the wrong
  sequences.
- `reference` advances the reference position over `N` (reference skip) CIGAR operations;
  samtools does not, so it writes the bases after a splice at the wrong positions (76 wrong
  bases on the tagged corpus, against none for ntools). A region with no data is labelled with
  its range (`>chr:1-100`), as one with data is; samtools labels it as the whole sequence. An
  `MD` tag that runs past the CIGAR is reported as an error; samtools fails without a message.
- `fastq`/`fasta` group records by name. samtools compares each record with a name it can
  overwrite: after a copy of a read without qualities is replaced by one with them, the next
  template is merged into the group, losing reads or pairing reads of different templates.
  Index reads (`--i1/--i2`) end with the barcode tag: samtools reads on past its end (as with
  `i*i*` and a barcode without a separator) and writes the following tag bytes as an index
  read. Barcodes that are not strings give no index reads; samtools reads their binary values
  as text. As in samtools, reads without qualities are written with `B` qualities, and `-v`
  applies only to index reads.
- `fastq`/`fasta --no-sc` remaps base modifications (`MM`, `ML`, `MN`) to the bases left after
  removing soft clips, in the read's original orientation; samtools keeps the tags of the
  untrimmed sequence, so its calls land on the wrong bases or past the read. Tags that do not
  match the sequence are removed, with a warning.
- `fixmate -M` trims the stale base modification tags of a hard-clipped secondary or
  supplementary alignment correctly where samtools does not. A call carries one `ML` value per
  modification code (`C+mh` has two): samtools drops one per call, so it writes a corrupt record
  and fails ("Corrupted aux data"). It counts the clipped bases on the primary's stored sequence
  rather than the sequenced one, so for a reverse-strand primary it drops calls or keeps the
  wrong ones (on long reads with synthetic hard-clipped supplementaries, 336 of 724 lost or
  changed calls, against none for ntools), and an `N` entry counts only `N` bases rather than
  every base. Tags that do not fit the primary's sequence (or whose `ML` is not a `B:C` array)
  are invalid and removed, as invalid tags are on any record; samtools drops `ML` and every tag
  after it, and writes corrupt records when `ML` has fewer values than `MM` calls. Tags that
  HTSlib accepts but that are outside the `MM` grammar (such as a delta written `+1`) are kept
  unchanged with a warning, with `MN` set to the primary's length so that readers skip them.
- `fixmate` soft clips only read bases when it trims an alignment that runs past the end of
  its reference: samtools also counts deletions and reference skips, so `3M10D2M` three bases
  before the end becomes `3M3D9S` for a 5-base read, a corrupt record (ntools writes `3M3D2S`).
  `checksum -a`, which only hashes the CIGAR, keeps samtools' result so checksums stay comparable.
- `checksum -a` and `checksum -z` trim a long CIGAR (more than 65,535 operations, stored in a `CG`
  tag) that runs past the end of its reference as they trim any other CIGAR; samtools 1.24 aborts
  with "double free or corruption".
- Only CRAM 3.0 and 3.1 can be read; CRAM 2.1 input is an error ("unsupported CRAM version"),
  except for `reheader`, which copies its containers. `quickcheck` reports a CRAM 2.1 header as
  unreadable. CRAM 4.x (HTSlib's experimental draft, which samtools reads) is an error too, rather
  than being decoded as CRAM 3, whose integer encoding differs.
- Headers that break the SAM specification are read as samtools reads them, lookups resolve the
  same way (the first `@HD` and `@RG` with an ID, the last `@PG` with an ID), and they are
  written back as samtools writes them: every line in its place with its tags as written, a BAM
  header byte for byte until something changes it (`--no-PG`). The differences:
  - samtools parses a BAM header only when it needs to (to add an `@PG` line, or to change or
    look up a header line), and fails then on what HTSlib's parser rejects; `view -c` or
    `view --no-PG` read such a file. ntools always parses the header, so it rejects these files
    for every command: a duplicate `@SQ` name, a duplicate name in the binary reference
    dictionary, and a malformed line (e.g., `@CO` without a tab).
  - An `@SQ` line in a BAM whose `LN` is zero, missing or invalid takes the binary dictionary's
    length; samtools fails when it parses it ("invalid or no LN").
  - A BAM header text padded with NULs ends at the first NUL. samtools writes the NULs into SAM
    output with `--no-PG` (`view -h`), which neither tool can then read.
  - ntools warns about a duplicate `@RG` ID when it reads the header; samtools only once it
    parses it (see above).
  - `merge` merges each input's read groups and programs, one per ID; samtools copies (and
    renames, if their IDs collide) every `@RG` and `@PG` line, so an input with duplicate IDs
    gives it more lines.
- An empty (0-byte) SAM file reads as a file with an empty header and no records; samtools fails
  to read its header.
- Reading a CRAM detached record with the mate-reverse flag but not the paired flag (flag 52),
  ntools keeps the flags as stored; HTSlib adds the paired flag (53).
- `bgzip -i -o OUT` with several inputs indexes the whole output. HTSlib restarts the index at
  each input, so its `.gzi` covers only the last one, with offsets that point at the wrong data.
- `bgzip -b N` at exactly the end of the data writes nothing and succeeds; HTSlib fails ("BGZF
  block offset set beyond block size").
- `bgzip --version` names ntools.
- `faidx` and `fqidx` write `.fai` and `.gzi` files (and the `.fai` of `--write-index`) to a
  temporary file renamed over the destination once complete, so a failure writing them leaves
  any existing index as it was; samtools leaves a partial `.fai` (both exit 1). As in samtools,
  a BGZF input truncated inside a block is an error that writes no index, and one ending exactly
  at a block boundary only warns that its EOF marker is absent.
- `dict` reads plain, gzip, and BGZF FASTA (and stdin, as `-` or when no input is given), as
  samtools does, but a truncated or corrupt gzip input is an error; samtools writes a dictionary
  of the sequences it read before the damage (shortening the last one) and exits 0. An input
  with no sequences gets a warning (samtools writes the `@HD` line silently); both exit 0.
- Writing CRAM, ntools attaches a pair as mates (leaving their mate fields and template lengths
  to the decoder) only when HTSlib's decoder reproduces them; HTSlib's encoder also attaches
  pairs that start at the same position with the second segment first and the shorter, whose
  template length signs its decoder then swaps.
- An output (`-o`) that is one of the inputs (the same device and inode) is refused before it is
  opened, with the input left untouched; samtools truncates the input first.

## Performance

Measured on an Intel Xeon Gold 5418Y VM (24 vCPUs, 122 GB) against samtools/HTSlib 1.24, which
compresses with libdeflate, as wall time. `-@ N` adds N threads in both tools. Corpus rows
("short" is 1M synthetic 150 bp paired reads, "long" 9,784 nanopore reads with MM/ML tags) are
medians of 5 warm runs of ntools 0.9.0 on 2026-11-07. Whole-genome rows (GIAB HG003 35x Illumina,
HG002 PacBio Revio, HG003 ONT R10.4.1) are one run each, at `-@ 16`, of ntools 0.2.0 on
2026-11-04, which compressed with libdeflate, except `view BAM->CRAM 3.1`, re-timed with ntools
0.9.0 on 2026-11-08. ntools now compresses with zflate ([Compression](../README.md#compression)),
which writes BAM faster.

| Operation | samtools (s) | ntools (s) | Speedup |
|---|---:|---:|---:|
| view BAM->SAM, long | 2.120 | 1.560 | 1.36x |
| view BAM->BAM, long, -@4 | 3.350 | 1.790 | 1.87x |
| view CRAM->SAM, short | 0.870 | 0.930 | 0.94x |
| view CRAM->SAM, long, -@8 | 0.550 | 0.330 | 1.67x |
| sort, short, -@4 | 0.970 | 0.390 | 2.49x |
| index, long, -@4 | 0.320 | 0.270 | 1.19x |
| markdup, short | 3.800 | 1.260 | 3.02x |
| stats, short | 1.650 | 1.010 | 1.63x |
| mpileup -f (BAQ), short, -@8 (samtools mpileup has no threads) | 27.730 | 3.870 | 7.17x |
| flagstat, 35x HG003 short reads (BAM) | 92.150 | 70.360 | 1.31x |
| view -c, 35x HG003 short reads (CRAM 3.1) | 23.460 | 5.060 | 4.64x |
| flagstat, 35x HG003 short reads (CRAM 3.1) | 30.420 | 7.560 | 4.02x |
| stats, 35x HG003 short reads (BAM) | 1440.650 | 553.840 | 2.60x |
| view BAM->CRAM 3.1, 35x HG003 short reads | 227.240 | 270.010 | 0.84x |
| collate, fixmate -m, sort, markdup (piped), 35x HG003 short reads | 2869.200 | 2909.910 | 0.99x |
| stats, HG002 Revio (BAM) | 2169.500 | 293.970 | 7.38x |
| sort -n, HG002 Revio (BAM) | 1011.740 | 646.610 | 1.56x |
| view CRAM->BAM, HG003 ONT (CRAM) | 745.710 | 445.340 | 1.67x |

Counting CRAM records (`view -c`, `flagstat`, `idxstats`) decodes only the fields the command
reads, as samtools does, and also selects tags one by one (`-d TAG`, the read group filters),
which HTSlib does not: speedups of 1.06x to 6.49x on the corpus as CRAM 3.1. `bgzip` of the GIAB
HG002 VCF (2.8 GB) takes 12.4 s against HTSlib's 16.1 s single-threaded, and 1.5 s against 2.4 s
at `-@ 8`, at its default level 6, for a file 5% larger than HTSlib's default.

### Where samtools is faster

- Writing CRAM (`view -C`) takes more time than HTSlib: on 2 million of the 35x short reads,
  single-threaded, 1.25x samtools' time (8.9 s against 7.1 s), for a file 0.3% smaller than
  samtools'; on the whole 35x genome at `-@ 16`, 1.19x (270 s against 227 s), with 1.6 times its
  CPU time. ntools 0.2.0 took 2.1x on both, and 1.05x on the Revio long reads; since then, read
  names are compressed as HTSlib compresses them.
- `view` CRAM to SAM of the short reads is slightly slower single-threaded (0.94x). On the 69 GB
  ONT CRAM, `view -c` and `flagstat` (0.86x and 0.89x in 0.2.0) mostly wait for the file to be
  read and use less CPU time than samtools.
- The whole-genome `collate | fixmate -m | sort | markdup` pipeline was level in 0.2.0 (0.99x),
  and its temporary files took 155.5 GB at peak against samtools' 79.8 GB. ntools now frees the
  space of its temporary files as it merges them: on GIAB HG003 chr20 at `-@ 2`, the pipeline
  peaks at 2.0 GB of temporary files against samtools' 2.4 GB (3.7 GB before), and takes 81 s
  against 153 s. `markdup` stops speeding up past `-@ 4`, limited by the thread that matches
  duplicates.
- Memory: `collate -Ou` takes 350 MB against samtools' 26 MB on the corpus, `mpileup -M` on the
  long reads 509 MB against 43 MB, `stats` on the long reads at `-@ 8` 843 MB against 660 MB,
  and writing CRAM 113 MB against 28 MB on the 2 million short reads.
