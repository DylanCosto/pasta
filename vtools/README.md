# vtools

A bcftools-compatible command-line toolkit built on the [pasta](https://github.com/DylanCosto/pasta/tree/main/pasta) library. It aims to
produce the same results as bcftools 1.24 while being faster, especially with threads.

Its engine is `pasta::variant`, a Rust port of HTSlib's VCF/BCF data model (`bcf_hdr_t`,
`bcf1_t`), VCF parser and formatter, BCF codec, CSI/TBI indexing, and the synced reader that
reads several files in step, written from the HTSlib 1.24 sources so that output matches byte
for byte.

## Install

```sh
cargo install --path vtools
```

## Commands

| Command | bcftools equivalent | Notes |
|---|---|---|
| `view` | `view` | VCF/BCF conversion (`-O v\|z\|b\|u[0-9]`, `-l`, `-o` with the type taken from the name as bcftools does, `--threads`, `-W/--write-index[=csi\|tbi]`, `--no-version`); `-h/-H/--with-header`; regions through the index (`-r`, `-R`, positional regions, `--regions-overlap`); streamed targets (`-t`, `-T`, `^` exclusion, `--targets-overlap`); samples (`-s`, `-S`, `^`, `--force-samples`), `-G`, `-a`, `-A/-AA`, `-I`; filters `-c/-C`, `-q/-Q` (`:nref\|alt1\|minor\|major\|nonmajor`), `-f`, `-g`, `-k/-n`, `-m/-M`, `-p/-P`, `-u/-U`, `-v/-V`, `-x/-X`, `-i/-e` (filter expressions, as in `filter`); `--verbosity`; `-o/--output-file`; option values may start with `-` (`--threads -1`, `-C -1`), as getopt takes them; without an input on a terminal, bcftools' usage (exit 1) |
| `filter` | `filter` | `-i/-e` with bcftools' whole expression language (no Perl): INFO and FORMAT tags with vector and sample subscripts (`TAG[1]`, `TAG[0,2-]`, `FMT/AD[0:1]`, `AD[*:GT]`, `GT[@samples.txt]`), GT classes (`het`, `hom`, `hap`, `mis`, `ref`, `alt`, `RR`, `RA`, `AA`, `aA`, `R`, `A`) and genotype strings, `QUAL`, `FILTER` (`=`, `!=`, `~`, `!~`), `ID` (also `@file`), `CHROM`, `POS`, `REF`, `ALT`, `TYPE`, `ILEN`, `N_ALT`, `N_SAMPLES`, `N_MISSING`, `F_MISSING`, computed `AC`, `AN`, `AF`, `MAC`, `MAF`, functions `MAX`, `MIN`, `AVG/MEAN`, `MEDIAN`, `STDEV`, `SUM`, `ABS`, `COUNT`, `STRLEN`, `PHRED`, `binom`, `fisher`, `N_PASS`, `F_PASS`, `SMPL_*`/`s*`, arithmetic (`+ - * / %`), regexes (`~`, `!~`, `/i`), `&`/`&&` and `\|`/`\|\|` per-sample logic, `"."` as missing; `-s/--soft-filter` (`+` for `FilterN`), `-m/--mode +x`, `-g/--SnpGap INT[:TYPE]`, `-G/--IndelGap`, `--mask`, `-M/--mask-file` (`^`, `--mask-overlap`), `-S/--set-GTs .\|0`; `-r`, `-R`, positional regions, `--regions-overlap`, `-t`, `-T`, `--targets-overlap`; `-o`, `-O`, `--threads`, `-W/--write-index`, `--no-version`, `-v/--verbosity` |
| `head` | `head` | `-h/--headers`, `-n/--records`, `-s/--samples` (the later of `-n` and `-s` sets the count), `-v/--verbosity`; numbers as `strtoull` reads them (`0x10`, `010`, `-1`) |
| `annotate` | `annotate` | Annotations from a VCF/BCF (read in step through the indexes, `--pair-logic`/`--collapse snps\|indels\|both\|all\|any\|some\|exact\|none\|id`) or a tab-delimited file (tabix-indexed and streamed, or held in memory when it has END/TO columns, plain or bgzipped, BED by extension; `--min-overlap`, `--single-overlaps`); `-c` with every column form (`CHROM`, `POS`, `FROM`/`BEG`, `TO`/`END`, `REF`, `ALT`, `ID`, `QUAL`, `FILTER`, `INFO`, `FORMAT`, `[INFO/]TAG`, `FMT/TAG`, `-` to skip, `^` exclusions, `+`/`-`/`=`/`.` modes, `~ID`, `~INFO/END`, `~TAG`, `-POS`, `NEW:=OLD`, `ID:=TAG`, `INFO/TAG:=ID\|FILTER`, sample columns), `-C`, `-l` (`first`, `append`, `append-missing`, `unique`, `sum`, `avg`, `min`, `max`), `-x` (`ID`, `QUAL`, `FILTER`, `INFO`, `FORMAT`, `INFO/TAG`, `FMT/TAG`, `FILTER/NAME`, `^` keep lists, `##key` and structured header lines), `-k`, `-h`, `-H`, `-I [+]FORMAT`, `-m [+-]TAG`, `-i/-e` with `{}` values from the annotation line, `-s/-S`, `--rename-chrs`, `--rename-annots`, `--force`, `-r/-R/--regions-overlap`, `-o/-O/--threads/-W/--no-version/-v`. With `--threads`, a single VCF input is parsed and the output formatted or encoded on worker threads |
| `query` | `query` | `-f` with every field of bcftools' format strings (`%CHROM`, `%POS`, `%POS0`, `%END`, `%END0`, `%ID`, `%REF`, `%ALT{N}`, `%FIRST_ALT`, `%QUAL`, `%FILTER`, `%TYPE`, `%IS_TS`, `%MASK`, `%LINE`, `%INFO`, `%INFO/TAG{N}`, `%TAG`, `%FORMAT`, and in `[ ]` `%SAMPLE`, `%GT`, `%TGT`, `%IUPACGT`, `%TAG{N}`, `%PBINOM(TAG)`, `%TBCSQ`; `%N_PASS(EXPR)` and filter functions such as `%MAX(FMT/DP)`; the hidden `%VKX`, `%RSX`, `%_CHROM_POS_ID`), `\n`/`\t` escapes; `-H/-HH`; `-i/-e` with per-sample results and `-F`; `-l`; `-s/-S` (`^`, `--force-samples`); `-r/-R/-t/-T`, `--regions-overlap`, `--targets-overlap`; `-u`; `-N`; `-o`; `-v`; several input files; the hidden `-a`; `--verbosity`; `--threads` (decompression, a vtools addition) |
| `convert` | `convert` | VCF/BCF output (`-O`, `-o`, `--threads`, `-W`, `--no-version`); `-i/-e`, `-s/-S`, `-r/-R/-t/-T` and overlap options; `--gvcf2vcf -f`; `-g/--gensample` (`--tag GT\|PL\|GP`, `--3N6`, `--vcf-ids`, `--sex`, `--keep-duplicates`) and `-G/--gensample2vcf`; `-h/--haplegendsample` and `-H/--haplegendsample2vcf`; `--hapsample` and `--hapsample2vcf` (`--haploid2diploid`, `--vcf-ids`, `--sex`); `--tsv2vcf -f` with `-c` columns (`ID`, `CHROM`, `POS`, `AA`, `REF`, `ALT`, `-`) and `-s/-S`; `-v/--verbosity` |
| `index` | `index` | CSI (`-c`, default, `-m` min shift) for BCF and bgzipped VCF, TBI (`-t`) for VCF (`-c`, `-t`, and `-m` apply in command-line order, as in bcftools: `-t -c` still writes TBI); `-o/--output-file`, `-f`, `--threads`, `-v`; the BGZF and EOF-marker checks on standard input (a pipe cannot be checked); `-n`, `-s`, `-a` statistics from an existing index (`FILE##idx##INDEX` names it, or a `.csi`/`.tbi` name alone) |
| `isec` | `isec` | Two or more bgzipped and indexed files (`-l` file list); `-n [+-=~]N` (`~` bitmask), `-C`, and the four-way split of two files; `-p DIR` (`0000.*` files, `sites.txt`, `README.txt`, the indexes bcftools makes) with `-w LIST`, or one file's records with `-w N` (`-o`, `-O v\|z\|b\|u[0-9]`, `-W[=csi\|tbi]`), or the sites list; `-c none\|snps\|indels\|both\|all\|any\|some\|id` (cumulative); `-f`; `-r`, `-R`, `--regions-overlap`, `-t`, `-T` (`^`), `--targets-overlap`; `--threads`, `--no-version`, `-v`. `-i/-e` (one, or one per file, `-` for none) |
| `sort` | `sort` | `-m` (k/M/G), `-T`, `-o`, `-O`, `-W/--write-index`, `-v`. Sorts in memory when the input fits in `-m`; otherwise sorted blocks spill to BGZF temporary files merged as bcftools merges them (32 files per layer), with the last block kept in memory. Memory is counted as bcftools counts its packed records, so a given `-m` spills the same blocks and prints the same `Merging N temporary files` |
| `reheader` | `reheader` | `-h`, `-s/-N` (names or `old new` pairs, with `\` escapes), `-n`, `-f` (`.fai`: update lengths, drop and add contigs), `-o`, `--threads`, `-T` (ignored), `-v`. Bgzipped VCF and BCF keep their compressed record blocks; BCF records are still read and checked against the new header |
| `concat` | `concat` | Plain concatenation (bgzipped VCF to `-Oz` and BCF to `-Ob` copy the compressed blocks of unchanged records), `-n/--naive`, `--naive-force`, `-a` (HTSlib's synced reader with `-d snps\|indels\|both\|all\|exact\|none`, `-D`, `-r`, `-R`, `--regions-overlap`), `-l/--ligate` with `--ligate-force`, `--ligate-warn`, `-c`, `-q`; `-f`, `-G`, `-o`, `-O`, `--threads`, `-W`, `--no-version`, `-v 0\|1` |
| `tabix` | `tabix` (bcftools' minimal subcommand) | Index (`-p gff\|bed\|sam\|vcf` or by file name, `-s`, `-b`, `-e`, `-0`, `-S`, `-c`, `-m` for CSI, `-f`), query regions, `-a` to print all lines |
| `norm` | `norm` | Left-alignment against an indexed FASTA (`-f`, plain or bgzipped; `.fai` and `.gzi` built and saved when missing), `-c e\|w\|x\|s` with bcftools' messages, `-N[NUM]`/`--do-not-normalize`, `-g/--gff-annot` (right-align in forward-strand transcripts); `-m +\|-[snps\|indels\|both\|any]` with INFO/FORMAT Number=A/R/G splitting and joining, `--multi-overlaps`, `--keep-sum AD`, `-s/--strict-filter`; `-a/--atomize`, `--atom-overlaps`; `-d snps\|indels\|both\|all\|any\|exact\|none`, `-D`; `--old-rec-tag`, `-S pos\|lex`, `-w`, `--force`, `-i/-e` (through the filter hook), `-r/-R/-t/-T` with overlap modes, and the shared output options; `-h/--help` prints bcftools' usage. With `--threads`, VCF lines are parsed and output records formatted or encoded on worker threads, written one record at a time as a serial run writes them (so a fatal error loses the same buffered output) |
| `consensus` | `consensus` | Every 1.24 option: `-f` (plain, gzip or bgzip FASTA, or standard input; `>chr:from-to` headers as regions), `-H 1\|2\|N\|R\|A\|L\|S\|LR\|LA\|SR\|SA\|I\|NpIu`, `-s`/`-S` (`^`, `-` for REF,ALT), IUPAC codes from the GTs of several samples, `-I`, `-a`, `-M`, `-c` chain files, `-m` (BED or CHROM POS END) with `--mask-with CHAR\|uc\|lc` per mask, `--mark-del`, `--mark-ins`, `--mark-snv`, `-p`, `-o`, `-i/-e`, `--regions-overlap`, `-v`; `<DEL>`, `<*>`, `<NON_REF>`, and the overlap warnings and errors as bcftools words them. `--threads` (decompression) is a vtools addition |
| `stats` | `stats` | One or two inputs, every section bcftools prints (SN, TSTV, SiS, AF, QUAL, IDD, ST, DP, PSC, PSI, HWE, VAF; FS with `-E`, ICS/ICL with `-F`, USR with `-u`): `--af-bins`, `--af-tag`, `-1`, `-d`, `-E`, `-f`, `-F`, `-I`, `-r`, `-R`, `--regions-overlap`, `-s`, `-S`, `-t`, `-T`, `--targets-overlap`, `-u`, `--threads`, `-v`, `--verbosity=INT`; `-i/-e` through the filter engine (applied to each input). Two indexed inputs are read with the synced reader: sets for records only in A, only in B, and in both, paired with `-c snps\|indels\|both\|all\|any\|some\|none`; with `-s/-S`, the genotype concordance of the samples both inputs have (GCsAF, GCiAF, NRDs, NRDi, GCsS, GCiS, GCTs, GCTi) and with `-v` the per-site PSD and DBG lines |
| `call` | `call` | The multiallelic (`-m`) and consensus (`-c`) callers with every 1.24 option: `-A`, `-*`, `-a`/`-f` (`GQ`, `GP`, `PV4`), `-F AN,AC`, `-G FILE\|-`, `--group-samples-tag`, `-g` gVCF blocks, `-i`, `-M`/`-N`, `-V snps\|indels`, `-v`, `-C alleles\|trio`, `-n`, `-p`, `-P`, `--ploidy` (and `-X`/`-Y`), `--ploidy-file`, `-s`/`-S` (lists, sexes or ploidies, PED files, `^`), `-r`/`-R`/`--regions-overlap`, `-t`/`-T`, `-o`/`-O`/`--threads`/`-W`/`--no-version`/`--verbosity`, `-h`/`-?` (usage) |
| `mpileup` | `mpileup` | SAM/BAM/CRAM inputs (`-b`, several files, read groups as samples: `-s`, `-S`, `-G`, `--ignore-RG`); `-f` (plain or bgzipped FASTA; `.fai` and `.gzi` built and saved when missing)/`--no-reference`; regions through the index (`-r`, `-R`) and streamed targets (`-t`, `-T`, `^`); read filters `-q`, `-A`, `--ls/--ns/--lu/--nu` (`--rf`, `--ff`; numbers or flag names), `-x`, `-d`, `-6`, `-C`; BAQ (`-B`, `-D`, `-E`, `-M`); base qualities `-Q`, `--max-BQ`, `--delta-BQ`; every `-a` annotation (`?` lists them); `-g` gVCF blocks; indels: the default model, `--indels-2.0`, `--indels-cns` (edlib), `-I`, `-L`, `-m`, `-F`, `-p`, `-o`/`--open-prob`, `-e`, `-h`, `--indel-bias`, `--del-bias`, `--score-vs-ref`, `--indel-size`, `--seqq-offset`, `--poly-mqual`, `--ambig-reads`, `-X` profiles (`-X list`); `--seed`; `-P` (accepted and ignored, as in bcftools 1.24); `-O`, `-o` (a number is `--open-prob`), `--threads`, `-W`, `--no-version`, `-v`; usage without arguments (`-U` and `--help` are rejected, as in bcftools) |
| `gtcheck` | `gtcheck` | The cross-check of every pair of samples in one file, or query samples against `-g FILE` (read in step through the indexes); `-E INT` (genotype error model; `0` counts mismatches), `-u TAG1[,TAG2]` (`GT`, `PL`), `--no-HWE-prob`, `-H`, `--keep-refs`, `--n-matches INT` (negative: by HWE score), `-s/-S` with `qry:`/`gt:` (and `-`), `-p/-P` (`--plot`), `--distinctive-sites NUM[,MEM[,TMP]]`, `-i/-e` with `qry:`/`gt:` (and the old `-e INT`), `-r/-R/-t/-T`, `--regions-overlap`, `--targets-overlap`, `--dry-run`, `-o`, `-O t\|z[0-9]`, `-v`; every report section (the header, timing and INFO lines, DCv2, DS). `-a`, `-c`, `-G`, and `--target-sample` stop with bcftools' messages (unimplemented or deprecated in 1.24). `--threads INT` (a vtools addition) decompresses on threads and shares the pair sums out to them when there are enough pairs |
| `merge` | `merge` | Two or more files (`-l` file list, `--force-single` for one), indexed (`-r`, `-R`, `--regions-overlap`) or streamed (`--no-index`/`--force-no-index`); `-m none\|snps\|indels\|both\|all\|any\|exact\|id\|snp-ins-del` (prefixes, a trailing `,`, `*`/`**` to trim unobserved alleles); alleles normalized across files with Number=A/R/G INFO and FORMAT values remapped (unknown `<*>`/`<NON_REF>` values filled in); `-i TAG:sum\|avg\|min\|max\|join` (and `-`, with bcftools' defaults); `-M TAG:.\|max\|NUM`; `-F x\|+`; `-0`; `-g -\|REF.FA` gVCF blocks; `-L INT` (LAA and `L*` fields); `--force-samples` (`2:NAME`), `--use-header`, `--print-header`; `-f`; `-o`, `-O`, `--threads`, `-W`, `--no-version`, `-v`; `-h` prints bcftools' usage. With `--threads`, VCF output lines are formatted on worker threads, and a fatal error loses the same buffered output as a serial run |
| `roh` | `roh` | The HMM of `vcfroh.c` and `HMM.c` with bcftools' numbers to the last digit: PL likelihoods or `-G/--GTs-only`; allele frequencies from INFO/AC and AN, `--AF-tag`, `--AF-file` (tabix-indexed `CHROM POS REF,ALT AF`), `--AF-dflt`, `-e/--estimate-AF [GT,\|PL,]FILE\|-`; `-m/--genetic-map` (IMPUTE2 format, `{CHROM}` in the name, plain or gzipped), `-M/--rec-rate`; `-a/--hw-to-az`, `-H/--az-to-hw`, `-V/--viterbi-training`; `-b/--buffer-size INT[,INT]` (negative for megabytes); `-i/--ignore-homref`, `-I/--skip-indels`, `--include-noalt`; `--include/--exclude`; `-s/-S` (`^`); `-r/-R/-t/-T`, `--regions-overlap`, `--targets-overlap`; `-o`, `-O s\|r\|z`; `--threads` (decompression, and the samples' HMMs on workers); `-v/--verbosity` |
| `csq` | `csq` | Haplotype-aware consequences with every 1.24 option: `-f`, `-g` (GFF3, plain or gzipped; the GFF reader is shared with `norm -g`), `-p a\|m\|r\|R\|s`, `-l/--local-csq`, `-s/-S` (`-` to ignore samples, `^`), `-n/--ncsq`, `-B/--trim-protein-seq`, `-b/--brief-predictions`, `-c/--custom-tag`, `-C/--genetic-code` (`l` lists the tables), `-G/--greedy`, `--force`, `--unify-chr-names`, `--dump-gff`, `-i/-e`, `-r/-R/-t/-T` with overlap modes, `-o`, `-O b\|u\|z\|v\|t[0-9]` (`t`: the tab-delimited text output), `--threads`, `-W`, `--no-version`, `-v/--verbosity`; INFO/BCSQ and the FORMAT/BCSQ bitmask (`-q` fails as in bcftools) |
| `plugin`, `+NAME` | `plugin`, `+NAME` | Built-in plugins (no shared libraries): `+fill-tags`, `+setGT`, `+split-vep`, `+counts`. The shared options of `vcfplugin.c` before `--`: `-i/-e`, `-r/-R`, `--regions-overlap`, `-t/-T` (`^`), `--targets-overlap`, `-o`, `-O`, `--threads`, `-W[=csi\|tbi]`, `--no-version`, `-l` (lists the four), `-v`, `-V`, `-h` (the plugin's usage); the input file is the first argument that is not an option (or standard input). `+fill-tags -t` with every tag (`AN`, `AC`, `AC_Hom`, `AC_Het`, `AC_Hemi`, `AF` (from AC/AN without samples), `MAF`, `NS`, `HWE`, `ExcHet`, `END`, `TYPE`, `F_MISSING`, `VAF`, `VAF1`, `all`, `INFO/`/`FORMAT/` prefixes) and `[INFO/\|FORMAT/]TAG[:N]=[int(\|float(]EXPR[)]` expressions through the filter engine, per population with `-S`; `-d`, `-l`. `+setGT -t ./.\|./x\|.\|a\|q\|?\|b:TAG<CMP>VAL\|r:FLOAT` (combinable) and `-n .\|0\|c:GT\|i\|m\|M\|X\|p\|u` (and combinations such as `0p`, `pM`), `-i/-e` with per-sample results, the binomial test, `-s` (glibc's `drand48`). `+split-vep -a`, `-A`, `-c` (names, indexes, ranges, `-`, `:TYPE`), `--columns-types`, `--csq-field`/`--consequence-field`, `-d`, `-f` (the `query` format engine), `-g [+]FILE`, `--gene-list-fields`, `-H/-HH`, `-l`, `-p`, `-s TR:CSQ:PRN` (`all`, `worst`, `primary`, `pick`, `mane`, `FIELD=\|!=\|~\|!~VALUE`; severity ranges `CSQ`, `CSQ+`, `CSQ-`; `all\|worst`), `-S`, `-u`, `-x/-X`, `-i/-e` (undefined subfields become columns), `-r/-R/-t/-T` and overlaps, `-o`, `-O` (gzip text with `-f`), `-W`, `--no-version`, `-v`, and `--threads` (a vtools addition). `+counts`. With `--threads`, VCF lines are parsed and output records formatted or encoded on worker threads |

Commands take bcftools' shared options as bcftools does: options may be repeated (the last one
wins), `-O` and `-l` combine in command-line order, an output name ending in `.vcf`, `.vcf.gz`,
`.vcf.bgz`, or `.bcf` overrides `-O`, and `--write-index` writes `OUT.csi` (or `OUT.tbi`), or the
index named by `-o OUT##idx##INDEX`. Where bcftools writes `##bcftools_viewVersion` and
`##bcftools_viewCommand` header lines, vtools writes `##vtools_viewVersion` and
`##vtools_viewCommand` (neither with `--no-version`).

Command lines are read as bcftools' `getopt_long` reads them, for every command and plugin
(`src/shortopt.rs`, in front of each command's clap arguments): a long option may be shortened to
any unique prefix (`--mask-w`, `--no-vers`); a required value is the next word even when it
starts with `-` (`--threads -1`, `-s -x`, `stats -E -u`) or is attached (`-i'QUAL>10'`,
`--format=...`); an optional value is taken only attached (`-N1000`, `-Wtbi`,
`--write-index=tbi`, `--verbosity=4`), never from the next word; short options cluster
(`-HGv snps`, `-Ov -o out`); options and operands mix in any order (with `POSIXLY_CORRECT` set,
as in glibc, the first operand ends the options), and `--` ends the options.
An unknown option, a missing value, a value given to a flag (`--no-header=x`), or an ambiguous
prefix prints getopt's message (`view: unrecognized option '--x'`, `view: invalid option --
'Z'`, `view: option '--min' is ambiguous; possibilities: ...`) and then the command's usage, and
exits as the command's usage does (1; 255 for `call` and `csq`). `-h`/`--help` print the usage
to standard error and exit the same way, where bcftools does (where `-h` is a real option,
such as `view -h` or `head -h INT`, it keeps its meaning). `tabix` uses plain `getopt`, as
bcftools' does: an unknown option is a warning and is ignored, and `--name` is a cluster of
short options. A unit test checks every short and long option of the 1.24 sources (the table in
`src/shortopt/bcftools-1.24.tsv`) against each command's definition, including how every prefix
of every long name resolves.

`vtools` itself follows bcftools' `main`: no arguments print the command list to standard error
(exit 1); `help`, `--help`, or `-h` alone print it to standard output (exit 0); `help COMMAND`
runs `COMMAND` without arguments (its usage); `--version`, `-v`, or `version` print the version,
and `--version-only` prints `VERSION+pasta-VERSION`; `+NAME` is `plugin NAME`; anything else is
`[E::main] unrecognized command 'NAME'` (exit 1), including the dropped `cnv`, `polysomy`, and
`som`.

Exit statuses follow bcftools: 255 for errors bcftools reports with `error()`, and HTSlib's
synced-reader codes for read errors (8 for a malformed VCF line, 9 for a malformed BCF record, 6
when a BGZF file lacks its end-of-file marker, after printing every record). Where bcftools
takes damaged input for the end of the file and exits 0, or only warns that a BGZF input it read
to the end lacks its EOF marker block, vtools fails with 255 (see Damaged input under
Compatibility). Warnings and errors on stderr are worded as HTSlib and bcftools word them
(`[W::vcf_parse] ...`), including the warnings HTSlib prints only once.

`--threads N` adds worker threads throughout. Beyond BGZF compression and decompression, which
is what bcftools uses its threads for, vtools parses VCF lines, decodes BCF records, and formats
or encodes output on the workers, while keeping the output (and stderr) identical to a serial
run: a line that needs a header change (an undefined contig or tag) is handed back to the main
thread, which updates the header as HTSlib does. A failed write (a full disk) stops the
compression threads at once and is reported as bcftools reports it: with the command's error for
a record it cannot write (`[main_vcfview] Error: cannot write to FILE`, exit status 255), also
when the threads notice it only while the output is closed. Only a failure on the last block or
the end-of-file marker gets the command's close message (`Error: close failed FILE`). A write to
a pipe whose reader went away (`| head`, a FIFO), on standard output, an `-o` FIFO or standard
error and from whichever thread, kills vtools with SIGPIPE and no message, as it kills bcftools
(status 141 in the shell): Rust programs ignore SIGPIPE, so vtools restores the default action
at startup.

Inputs named by a pipe are read as standard input is, without seeking: process substitution
(`<(...)`), `/dev/stdin`, and FIFOs work wherever bcftools reads them, including `mpileup`
alignment inputs. Options that need an index (`-r`, `merge` without `--no-index`) fail on them
as in bcftools, and commands that read an input twice (`concat`, `annotate -a`) read a pipe on
from where the first reading left it, as bcftools does.

## Compatibility

vtools has bcftools 1.24's commands, options, and output formats. Where bcftools or HTSlib gives a
wrong result (a crash, wrong output, or output another tool rejects), vtools gives the right one
instead of copying it; those cases are below, with the other intentional differences.

BGZF blocks are laid out as HTSlib lays them out: a block of compressed VCF ends at a record, as
`vcf_write` ends it, unless the record alone is longer than a block, and BCF blocks hold 0xff00
bytes. So output cut short by a fatal error ends at a record, and index files (CSI and TBI, from
`index`, `tabix`, and `--write-index`) have the layout of HTSlib's.

Intentional differences:

- BCF and VCF.gz are written at level 6 unless a level is given (`-O b5`/`-O z5` in HTSlib's
  terms): 1.3–1.4 times as fast as bcftools' default, with files 5% larger; `-O b6`/`-O z6` gives
  bcftools' size ([Compression](../README.md#compression)). Index files are written at level 5.
  The compressed bytes are zflate's, not HTSlib's; they decompress to the same data.
- `query -N` (`--disable-automatic-newline`) takes no value. bcftools declares the long form
  with a required argument, so `--disable-automatic-newline` swallows the next argument there.
- `query -H` with `%LINE` prints no header, as bcftools intends; bcftools fails writing the empty
  header ("cannot write to standard output").
- `%MASK` prints nothing, as in bcftools, which never gives the format engine its readers.
  (`query` with several input files reads them in step through their indexes and prints the
  first file's records, ordered as bcftools' synced reader pairs them.)
- `query --threads` decompresses on worker threads; bcftools query has no such option.
- Positions from 2^31 on (possible in VCF text): `query` (`%POS`, `%POS0`, `%END`, `%END0`) and
  `convert -g/--hapsample` (the position and the `CHROM:POS_REF_ALT` IDs) print them as they
  are; bcftools prints them through 32-bit integers, as negative or wrapped numbers.
- `query -s/-S` removing every sample, on VCF input with an undeclared contig or tag: HTSlib
  keeps the whole sample dictionary when no sample is left, so the header update for the
  undeclared name brings the samples back, and bcftools evaluates `-i/-e` with the samples'
  values in the records after it. vtools keeps the samples removed.
- `filter`, `query`, and `convert` with both `-r` and `-R` (or `-t` and `-T`, `-s` and `-S`)
  use the later one, inline list or file; bcftools takes the later value but reads it as a file
  name once the file form was given anywhere, so `-R FILE -r REGION` fails there.
- `query -v` applies `--regions-overlap` and `--targets-overlap` to every listed file;
  bcftools ignores them with `-v`.
- `%N_PASS(EXPR)` with an expression on site values only (`%N_PASS(QUAL>3)`) prints 0: no sample
  is tested. bcftools reads the per-sample results the expression never set, and crashes.
- `convert -g/-h/--hapsample` with an output (sample, hap, or legend) file that cannot be created
  stops with `Error writing FILE: ...`; bcftools crashes.
- `%RSX` of an ID shorter than two characters (`.`) prints `00000000`; bcftools reads past the
  end of the ID.
- `convert --haploid2diploid` (`%_GT_TO_HAP2`) at sites where every genotype is haploid (one GT
  value per sample) prints each haploid allele as a homozygote; bcftools reads the next
  sample's allele as the second one.
- `convert -g --tag PL/GP` (`%_PL_TO_PROB3`, `%_GP_TO_PROB3`) for a sample with fewer than
  three values (a missing PL, or a site where no sample has more than one value) prints the
  missing probabilities as `nan` (PL) or `0` (GP); bcftools reads past the sample's values, into
  the next sample's or a previous record's.
- `convert --tsv2vcf` with genotype letters (`AA`) lists the ALT alleles in the order their
  allele numbers were given; bcftools lists them in A, C, G, T order, which mislabels the
  genotypes of sites where the second ALT seen sorts first.
- `convert --tsv2vcf` at a position past the end of its reference sequence writes REF as `.`
  in VCF output, as bcftools does in BCF output; bcftools' VCF output has an empty REF there.
- The gzipped outputs of `convert -g/-h/--hapsample` are written by a different compressor, so
  the compressed bytes differ (with `--threads`, the hap file is BGZF in both, compressed on the
  workers, or on the main thread when none can start).
- `view -g` with `-s`/`-S` classifies the genotypes of the selected samples. bcftools loops over
  the samples of the original header, reading past the end of the subset data.
- `view -A/-AA` with `-s`/`-S` trims the unobserved allele (`<*>`, `<NON_REF>`, `<X>`) of the
  subset record, giving the output of `view -s | view -A`. bcftools trims it with the original
  header after subsetting the samples; at a site where it trims the allele, it stops with
  `[E::bcf_remove_allele_set] Problem updating genotypes` and `[subset_vcf] Error: failed to
  trim the unobserved allele` (or `Broken VCF record` on a record with only GT), exit 255.
- When translating record IDs between headers (`Record::translate`, HTSlib's `bcf_translate`),
  INFO IDs of 128 or more that are rewritten in place are written after the type byte, as HTSlib
  does for FORMAT IDs; HTSlib writes them over the type byte, corrupting BCF output.
- Removing alleles (`view -a`, `view -s` with alleles left unused, `-A`) keeps each Number=G
  field's own width when a Number=LG field (local alleles, as `LPL`) comes before it in FORMAT.
  HTSlib's `bcf_remove_allele_set` keeps one width for both, so a `PL` after an `LPL` is cut to
  the `LPL`'s width: four alleles cut to three leave three of `PL`'s six values, with no
  warning.
- BCF output refuses a record with more than 65535 alleles or INFO fields, more than 255 FORMAT
  fields, or a shared or per-sample block of 2 GiB or more
  (`[E::bcf_write] Too many ALT alleles at CHROM:POS`, `Too many INFO entries`,
  `Too many FORMAT fields`, `The record at .. is too long`, `The FORMAT at .. is too long`), and
  the command fails as it does for any record it cannot write. The VCF parser refuses that many
  alleles, INFO fields, or FORMAT fields in both tools, but edits can make them: `norm -m+`
  joining onto a 65535-allele record, or `annotate` adding an INFO field to a record with 65535
  or a FORMAT field to one with 255. HTSlib keeps these counts in 16- and 8-bit fields, so
  bcftools crashes or writes a corrupt record there; vtools writes such records as VCF. HTSlib
  writes blocks up to 4 GiB, but the size of a single value of 2^31 elements or more wraps.
- A record's length on the reference (HTSlib's `get_rlen`: `%END`, the overlap of `-r`/`-t` with
  `--regions-overlap`/`--targets-overlap` 1 or 2, and the length written to BCF) takes the SVLEN
  of each `<DEL>`, `<DUP>`, `<INV>`, or `<CNV>` allele from that allele. HTSlib tests the flag
  of the allele before it at ALT 8, 16, 24, ...: there the SVLEN is ignored (a 50-base `<DEL>`
  at ALT 8 at POS 2000 ends at 2000, not 2050, and overlapping regions miss it), and at ALT 16,
  24, ... any SVLEN counts when ALT 8, 16, ... is such an allele (an `<INS>` with SVLEN 1000 at
  ALT 16 makes END 6000 at POS 5000). BCF input keeps the length stored in the file, as in
  HTSlib, so BCF written by bcftools keeps its wrong length.
- On truncated BGZF input, the `[E::bgzf_read_block]` message is printed only by `head` and
  `tabix` (without the block offset); the other commands report the read error ("BCF read
  error") with bcftools' exit status, or with 255 where bcftools exits 0 (see Damaged input).
- With `--threads`, a truncated or corrupt BGZF input gives the serial run's output, message,
  and exit status (the records before the damaged block, then "BCF read error", 9). bcftools'
  threaded reader decodes ahead and may instead fail while reading the header ("could not parse
  header", 255) and print nothing.
- Damaged input: a BGZF block with a bad CRC32 or a corrupt deflate stream, a block cut short,
  or a truncated or malformed BCF record. bcftools' `call`, plain `concat` (every output mode),
  and `head` take the failed read for the end of the input and exit 0, the records after it
  lost, printing only HTSlib's checksum or decompression message. vtools writes the same output
  (the records before the damaged block, in full) and then fails with status 255:
  `Error: BCF read error` (`call`), `Error: BCF read error: FILE` (`concat`),
  `[main_vcfhead] Error reading "FILE"` (`head`). Side files read as they go fail the same way
  where bcftools takes a failed read for their end: a tabix-indexed targets file (`-T`, and
  `roh --AF-file`) gives `Error: failed to read the targets: FILE: REASON` from every command,
  an `annotate -a` tab-delimited file `Error: failed to read the annotations: FILE: REASON`, and
  `stats -E` `Error: failed to read the exons: FILE: REASON`; `call -T` read as a list fails
  with `Error: failed to read the targets file FILE: REASON` (bcftools crashes, as it does with
  an `annotate -a` tab-delimited file and `-c CHROM,FROM,TO`, and in `call -c` on a damaged BGZF
  input) and `roh -m` with `Error: failed to read FILE: REASON` (bcftools fails with "Close
  failed"). The commands that read through the synced reader already fail as bcftools does. An
  unparsable VCF line still ends `call`, plain `concat`, and `head` quietly, as in bcftools
  (HTSlib reports the line), unless the input is plain gzip (not BGZF): vtools then reads the
  rest of the gzip stream first, whose CRC32 and length it checks for every member, and a
  corrupt deflate stream or a checksum that does not match makes the line a failed read: status
  255 (`Error: BCF read error`, `[main_vcfhead] Error reading "FILE"`), or 9 in `view`
  (bcftools' status for a failed read), where bcftools exits 0 or reports a VCF parse error (8
  in `view`). Garbage before such an error that happens to parse is caught at the member's end,
  as before; garbage that fails a command on its content fails it as in bcftools.
  `concat --naive` (bgzipped VCF and BCF alike) and `reheader` of bgzipped VCF copy compressed
  blocks without decompressing them, as bcftools does, checking only each block's header and
  length: a block with a bad CRC32 or corrupt deflate data is copied as it is, with status 0 in
  both. A BGZF input read to its end without the EOF marker block (all its blocks intact, as a
  writer killed between blocks or a tool that never writes the marker leaves it) gets HTSlib's
  warning `EOF marker is absent. The input may be truncated`, and then, once the command has
  written its output, fails with `Error: no BGZF EOF marker; the input may be truncated: FILE`
  (status 255), where bcftools exits 0 after the warning: `concat` (also `--naive`, which
  bcftools copies without a word), `sort`, `head`, `reheader`, `tabix -a` and tabix indexing,
  every command reading a pipe (which cannot be checked up front), and side files read to their
  end (targets and regions files, `annotate -a` tab-delimited files, `stats -E`, `roh -m`). The
  synced-reader commands still fail up front on such a file with status 6, as in bcftools; an
  index query that stops before the end is not affected; `index -f`, whose point is to index
  such a file anyway, keeps bcftools' warning and status 0.
- An empty output type (`-O ''`) is rejected as not recognised by every command that takes `-O`.
  bcftools reads the compression level from past the end of the empty argument (the text of the
  next argument).
- HTSlib's informational and debugging messages (`[I::...]`, `[D::...]`, printed with
  `--verbosity 4` or more, such as `[I::hts_idx_check_local] Using alignment file`) are not
  printed; warnings and errors are.
- `-r`/`-R`, `-t`/`-T`, and `-s`/`-S` given together: the later option wins, value and kind.
  bcftools keeps the "is a file" flag of an earlier `-R`, `-T`, or `-S`, so `-T FILE -t REGION`
  opens REGION as a file and fails.
- A regions or targets file with a tabix index is read as HTSlib streams it: regions in file
  order, without merging overlaps. HTSlib queries an indexed targets file per sequence instead,
  which visits the same lines for sorted files.
- `FILE##idx##INDEX` names (an index stored under another name, as HTSlib reads them) work
  wherever bcftools opens a file with `hts_open` or loads an index: inputs; `annotate -a` (VCF,
  BCF, tab-delimited with or without END/TO columns, BED; TBI or CSI) and `-h`; `-R`/`-T` region
  and target files; `filter --mask-file`; `consensus -m`; `call --ploidy-file` and `-T`;
  `csq -g` and `norm --gff-annot`; `roh -m`; the `convert` text inputs; `reheader -h`; `concat`;
  `mpileup` alignment files (BAI, CSI, or CRAI) and `-R`/`-T` files; `+fill-tags -S`; and
  `+split-vep -S`. A region, target, annotation, or mask file read whole rather than through its
  tabix index is read as BED when `FILE` ends in `.bed`, `.bed.gz`, or `.bed.bgz`. bcftools
  tests the whole name, which ends with the index's name, so it reads such a BED file as 1-based
  tab-delimited regions, one base off at each start.
- `norm -a`: bcftools reuses its atom structures without resetting the length of the REF/ALT
  prefix they share, so whether two atoms at one position conflict (and get `*` and genotype
  allele 2) can depend on earlier records. vtools resets it, giving the result bcftools gives
  for the record on its own (with `-r`, say). This changes 237 output records, at 155
  positions, of the 1000 Genomes chr22 sample of 200.
- `norm -m+` with a Number=G FORMAT field that has a single value per sample (`.`): bcftools
  reads values past the end of the record's array for the last samples; vtools reads missing
  values there.
- `norm -m-` with a missing Number=G string FORMAT value (`.`) followed by other samples:
  bcftools' `split_format_string` skips the missing value without moving on to the next sample,
  so every later sample keeps its unsplit value, cut to the longest value split before it, and
  the field is dropped when the first sample's value is missing (with `--force`, a count check
  on those unsplit values can drop it too). vtools splits every sample.
- `norm -m+` and `-m-` with Number=G FORMAT values (PL, GL, GP, and strings) of a sample of
  ploidy three or more (`0/0/1` with four PL values for one ALT): bcftools sizes these values
  for haploid and diploid samples only, so `-m+` stops with `could not merge FORMAT tag` (or,
  when the polyploid record comes first, writes values at the wrong places and drops the rest: a
  triploid with three alleles needs 10 values), `-m-` stops with `wrong number of fields` (with
  `--force` it drops the field), and string fields stop with `neither diploid nor haploid`.
  vtools takes a sample's ploidy p from its number of values (C(n+p-1, p) for n alleles) and
  joins and splits the values in VCF's genotype order: `-m+` puts each record's values at its
  genotypes' places among the joined alleles, a genotype no record has is missing, and where
  records share a genotype a later record's number replaces an earlier one's while a string
  keeps the first record's value, as bcftools does for diploid samples; `-m-` gives each ALT the
  p+1 genotypes of REF and that ALT. Samples of each ploidy keep their own number of values;
  haploid and diploid fields are joined and split as before, and joined GTs of ploidy three or
  more are what bcftools writes. `-m+` stops with `could not merge FORMAT tag` when, in a field
  where any sample has values of ploidy three or more in any of the joined records, a sample's
  number of values in some record does not fit its ploidy in the first record (numbers and
  strings, including a haploid/diploid difference in another sample; where bcftools exits 0 with
  values at the wrong places, when the polyploid record comes first, and for strings also when
  the polyploid values only appear in a later record, this is a new error), and when the joined
  values of such a field would take more than 1 GiB across the samples (512 MiB for strings;
  ploidy 150 over six alleles, say), where bcftools stops anyway.
- `norm -m` on BCF records with fewer samples than the header: bcftools reads past the end of
  the record's FORMAT data for the missing samples and fails on whatever it finds there (the
  count in its 'wrong number of fields' message varies from run to run); when nothing else fails
  first, vtools fails as `bcf_write` does on such records
  (`[E::bcf_write] Broken VCF record, the number of columns ... does not match the number of samples`,
  then `[flush_buffer] Error: cannot write to ...`), exit status 255.
- `norm` right-alignment (`-g` forward-strand transcripts, `--right-align`): bcftools pads each
  window with `aln_win`+1 bases but moves on by `aln_win`, repeating one base when a repeat
  spans more than 100 bp, and does not stop at the end of the sequence; vtools pads with the
  bases that follow and stops there. A one-base `<DEL>` (REF and the expanded ALT the same base)
  right-aligns to the end of the sequence: bcftools reads past it (garbage REF, or a crash),
  vtools puts it at the last base.
- `norm -c` (and `-f` generally) at positions above 2^31: bcftools passes the 0-based position
  to the FASTA reader as a 32-bit `int`, which wraps (a record at 2^32 is checked against the
  first base of its sequence, so `-c x` keeps it when that base matches); vtools reads the
  reference at the position itself (past the end of the sequence, a mismatch).
- `norm -d` at positions above 2^31: bcftools compares each record's position with the previous
  one kept in a 32-bit `int` (`prev_pos` in `flush_buffer`), which wraps, so it misses the
  duplicates there; vtools removes them.
- `norm -g` reads the GFF with the GFF3 reader `csq` uses (bcftools' `gff.c`), so it prints the
  same warnings and notes and fails on the same files. `--force` does not relax the CDS phase
  checks, as in bcftools `norm`.
- `csq` writes `##vtools_csqVersion` and `##vtools_csqCommand` where bcftools writes
  `##bcftools/csqVersion` and `##bcftools/csqCommand`, and its `-Ot` text output starts with
  `# This file was produced by: vtools +csq(...)` and `# The command line was:\tvtools +csq ...`.
- `csq`: where bcftools reads memory it never wrote (translated sequences compared past their
  end, the bytes after a REF shorter than its END-based length, bases other than A, C, G, T,
  a, c, g, t, n read through its 117-entry base table), vtools reads what an earlier translation
  left there or zeros, the next alleles' bytes, and non-ACGT codes respectively. bcftools
  release builds skip its assertions; where what such a build does is well defined, vtools
  does the same, its reads being bounds-checked: a triploid genotype sets bits of the next
  consequence, a start-lost frameshift node that other haplotypes extend reports their
  consequences, and a haplotype shorter than its transcript (overlapping CDS) gets the same
  wrapped padding and amino-acid positions as in bcftools' unsigned arithmetic. Where the
  skipped assertion is all that stops a read outside an array, vtools stops with the
  assertion's message instead: a GT allele index past the ALT alleles fails with `csq:
  assertion failed in test_cds: ial < rec->n_allele` (exit status 255), where bcftools reads
  past the alleles (garbage output, or a crash).
- `csq` on a VCF without `##contig` lines: bcftools reads the first name of the empty contig
  list to suggest `--unify-chr-names` and crashes; vtools goes on, names the contigs the VCF
  parser adds, and stops (or warns with `--force`) when one is not in the reference or the
  GFF, as bcftools does for declared contigs.
- `csq` fixes bcftools' use of recycled records: a `@POS` reference to a compound consequence
  printed at a record that was already written names that record (bcftools names whichever
  later record reuses its buffer), and a REF compared at an offset past its end (a record
  starting more than 10 bases before a transcript, with END making it overlap the CDS) is not
  compared with whatever memory follows. With `-s - -Ot`, bcftools drops the samples while
  reading, so `-i/-e` expressions on FORMAT fields see no samples (and `-i` crashes); vtools
  evaluates them on the samples. Duplicate names in `-s` are used once (bcftools crashes). A GFF
  line cut off right after the strand column is read as having no attributes (bcftools reads
  the rest of the previous, longer line's buffer), and translating a transcript with a CDS
  that ends before it starts is an error (bcftools wraps its length around and crashes).
- `csq` with VCF output, when the first record of a new sequence is one it does not annotate (no
  ALT, a `*` gVCF allele, or excluded by `-i/-e`) at the position of a buffered record of the
  previous sequence other than the last: bcftools flushes the previous sequence only before a
  record it annotates, so the new record takes over that position in its buffer index and
  bcftools stops with `This should not happen.. CHROM:POS` (exit status 255); vtools flushes the
  previous sequence first.
- Every `stats` check passes (also with `--threads`, and with two inputs); the checks of the
  differences below compare with the equivalent bcftools command.
- `stats A B` with `-s/-S` collects per-sample VAFs only when A's header defines FORMAT/AD, as
  bcftools does; when only B defines it, bcftools crashes on the first VAF of a record only in
  B, and vtools counts no VAF for it.
- `stats A B` with `-s/-S` compares the two inputs' genotypes by class, as bcftools does:
  hom-ref, het, hom-alt, and het of two ALTs. At a multiallelic site, genotypes of the same
  class with different alleles count as a match (`C/C` against `G/G`, both hom-alt), so the
  concordance tables (GCsAF, GCiAF, NRDs, NRDi, GCsS, GCiS, GCTs, GCTi) do not show such allele
  mismatches.
- `stats` fixes bcftools crashes and misparses: `-s ^LIST` and `-S ^FILE` exclude the listed
  samples (bcftools keeps the `^` on the first name and fails); `-u TAG:MIN` and
  `-u TAG:MIN:MAX` use the default for what is left out (bcftools reads the next argument as the
  bin count); `-I` with `-u` works (bcftools crashes); a missing `--af-tag` value puts the
  allele in the unknown-frequency bin and a missing `-u` Float value is skipped (bcftools
  crashes on both); `-d` with a step below 1 or a range too large for an int, `-u` bounds that
  are infinite or NaN, and `--af-bins` with fewer than two values are errors (bcftools divides
  by zero, crashes, or fails an assertion); `-d` and `-u` bins that the memory cannot hold
  (`-d 0,2000000000,1` under `ulimit -v 8000000`) stop with
  `Error: could not allocate N bytes of memory for --depth MIN,MAX,STEP` (or `--user-tstv TAG`)
  and status 255, where bcftools crashes. With `--threads`, each worker needs its own copy of
  those bins: when the memory holds them for fewer workers, `stats` runs with fewer, or without
  threads, as bcftools runs.
- `stats` IDD genotype counts and mean VAF: a symbolic or breakend allele has length 0, so its
  VAF goes to no printed bin. HTSlib leaves the length of such alleles unset, so bcftools reuses
  whatever an earlier record left in that slot, which can add genotypes to an indel length.
  Similarly, for an indel within 50 bases of the end of a contig, bcftools' `-F` repeat search
  reads past the fetched sequence; vtools reads it as the end of the sequence.
- `stats` with `-s` and `-S` together (or `-r`/`-R`, `-t`/`-T`): the later one is used, as a list
  or a file as its option says; bcftools uses the later value but reads it as a file if the file
  form was given anywhere.
- `stats` writes its report in one write, so `stats FILE | head -c 10` usually finishes before
  `head` closes the pipe, with status 0; bcftools writes its report 4 KB at a time with work in
  between, so it is killed by SIGPIPE (141). When the reader is already gone, both are killed by
  SIGPIPE.
- `stats -s` VAF of a genotype whose allele index is past the FORMAT/AD values of every sample
  (all AD values `.`, say, at a multiallelic or indel site): bcftools reads the next sample's
  values (as vtools does) and, for the last sample, whatever memory follows the AD data, which can
  add a VAF; vtools reads the missing value there and adds none.
- `stats -d MIN,MAX,STEP` labels each DP bin with its first depth; bcftools labels bin `i` as
  `i-1+MIN`, which is wrong when STEP is above 1.
- `stats --threads` adds the sums behind IDD's mean VAF in a different order, which could change
  a rounded value in the last printed digit in rare cases.
- FASTA references (`-f` of `mpileup`, `norm`, `csq`, `convert`, `merge -g`, and `stats -F`)
  named by a pipe (process substitution, `/dev/stdin` fed by a pipe, a FIFO) fail before they
  are read, with HTSlib's message for an index that cannot be saved next to them
  (`[E::fai_build3_core] Failed to open FASTA index NAME.fai : No such file or directory`).
  bcftools reads the pipe to index it and then fails to save the index (with `Permission denied`
  for `/dev/stdin`, and `Failed to make bgzf index NAME.gzi` for a bgzipped pipe), or, for a
  FIFO in a writable directory, saves `FIFO.fai` and waits forever to read the FIFO again. Where
  it can write `/dev/stdin.fai` (as root), bcftools saves the index of the pipe, then fails
  every read of the reference and still writes its output and exits 0.
- `mpileup -t ^REGIONS`/`-T ^FILE` (excluded targets) keeps the reads that are not entirely
  inside an excluded region, as bcftools' code comment says it should. bcftools 1.24 instead
  keeps only the reads that overlap an excluded region, so the positions it reports next to an
  excluded region have lower depth and positions further away are missing.
- `mpileup` writes `##vtoolsVersion` and `##vtoolsCommand` where bcftools writes
  `##bcftoolsVersion` and `##bcftoolsCommand`.
- `mpileup` reproduces what release builds of bcftools do where its code relies on undefined
  behaviour: `int` conversions of out-of-range doubles (scores of failed alignments) give
  `INT_MIN` as on x86 (except for BAQ, below), `--indels-2.0` alignments that fail bcftools'
  (disabled) assertions with empty windows score as probaln_glocal scores them, and the
  consensus candidate count and heap that bcftools never resets between positions are carried
  over. A `--indels-2.0` window that would read outside a read is skipped with a message
  (bcftools reads past the read).
- `mpileup` on CRAM checks the `-f` reference for slices without a reference checksum (slices
  with many reference sequences, as in long-read CRAMs), which bcftools decodes against whatever
  the reference holds: against the `M5` of a sequence dictionary beside it (`ref.dict`), or
  else against `@SQ` `LN`.
- `mpileup` BAQ keeps the quality of a base it is certain of, such as in a long read that
  matches the reference for thousands of bases. bcftools on x86-64 converts that base's
  infinite BAQ phred value to `INT_MIN` and sets its quality to 0; on arm64 it keeps the
  quality, as vtools does.
- `mpileup` reads SAM text as HTSlib's `sam_parse1` does (an RNAME missing from the header is
  unmapped, and so on), with HTSlib's warnings; HTSlib prints one for every record, vtools each
  distinct one once per input. As in ntools, a position whose reference sequence is dropped
  goes with it, and RNEXT `=` is the RNAME as written; neither changes the pileup.
- `mpileup` skips mapped reads without a sequence (SEQ `*`) whose CIGAR consumes query bases;
  bcftools reads their bases and qualities from past the end of the record (its output for
  them is garbage, or it crashes).
- `mpileup -S FILE` with a samples file that cannot be read is an error; bcftools ignores it
  and keeps every sample.
- `mpileup` decodes CRAM slices with reads that overhang the end of the reference (the bases
  past the end are `N`, as in ntools); HTSlib fails on such slices ("CRAM CIGAR extends beyond
  slice reference extents") and bcftools stops with "failed to read from input file".
- `mpileup` with several inputs stops with an error, before writing anything, when a file's
  reference sequences (names and lengths, in order) are not the first file's, or the first of
  them. bcftools takes every file's contig names from the first header, so it piles reads up on
  the wrong contig and exits 0, or, for a file with more contigs than the first, stops with
  "failed to read from input file" after writing the records before it.
- `mpileup --threads N` piles up on `N` workers with the output of a serial run. Regions
  (`-r`/`-R`), or the chromosomes of indexed inputs, are cut into chunks of 20-200 kb (longer
  with a large `-M`), planned from the inputs' indexes: a chromosome or region that no index
  query reads from is skipped, and a chunk with no reads is read by the chunk before it. So
  sparse inputs (panels, slices, empty files) open each input once per stretch of reads rather
  than once per 200 kb of the header's sequences, and are not slower than a serial run. A chunk
  starts reading a margin before its first position to rebuild the
  pileup there (buffered reads, BAQ already applied, waiting mates, and the state bcftools
  carries between positions: the BQBZ and MQSBZ values indel records repeat from the last SNP,
  and `--indels-2.0`'s candidate count), and its output is used only if a fingerprint of that
  state equals the state the previous chunk ended with; otherwise the previous chunk's pileup
  carries on serially. Deep pileups (more than 255 reads of a sample) are subsampled with
  `drand48` as in HTSlib; each chunk takes the generator state after the draws of the chunks
  before it. Positions with many costly indel realignments (long reads over long repeats)
  realign their reads on threads borrowed from a budget of `2N` that the workers share (each
  holds one while it piles up a chunk), so at most `2N` threads compute. Each chunk in flight
  reads every input on its own, so the chunks piled up ahead (at most `2N`) are as many as the
  soft open-file limit (`ulimit -n`) leaves room for; when not even one fits, the run is serial
  and keeps one reader of each input open, as bcftools does (the limit is not raised). gVCF
  output (`-g`) and unindexed inputs run serially, as does a run whose workers cannot start
  (`ulimit -u`); the serial run shares the `N` decompression threads between the inputs, and
  reads each input on a thread of its own, or on its main thread when that thread cannot start.
- `mpileup --threads N` with fewer open files allowed (`ulimit -n`) than its inputs need:
  bcftools opens each input once and keeps it open (loading its index right after it when there
  are regions), and fails on the first input it cannot open or index, before it prints anything.
  Unless the output is gVCF, vtools reads every header first, one file at a time, loading the
  index and closing the file, then opens the inputs again for the pileup in the same order. So
  it prints the `[mpileup] N samples in M input files` lines before the error, which names the
  first input it cannot open again. It never needs more open files than bcftools (with regions,
  one fewer), and the exit status is the same. Without `--threads`, vtools keeps each input open
  as bcftools does and fails at the same input.
- On a fatal error, bcftools exits without flushing its output buffers; vtools emulates HTSlib's
  file buffer sizes so the same output is lost. Only closing an output completes it: an output
  left on an error keeps the BGZF blocks already handed to the compression threads, but not the
  block being filled nor the end-of-file marker. So with `--threads`, an error found while
  reading or processing (an input that fails to open, a malformed record, an undefined tag in an
  expression, a failed write to the one output) leaves the serial run's output, byte for byte. A
  record out of order for `-W/--write-index` (unsorted, or a BCF contig that comes back) is
  found as it is written, so the output also ends where the serial run and bcftools, with or
  without threads, end it. A write error on one of several outputs (`isec -p`) is the exception:
  the other outputs' threads keep writing until the error is seen, so their threaded output can
  differ from the serial run's, as bcftools' threaded outputs also do.
- `isec -p` writes `vtools isec ...` as the command line in `README.txt`. Its indexes (`.tbi` for
  `-Oz`, `.csi` for `-Ob`, or as `-W` asks) are built while writing rather than by rereading each
  file after closing it; they are the same indexes. `--threads` also decompresses the inputs.
- `isec -l FILE` with several files on the command line reads them all, then the list; bcftools
  1.24 overwrites list entries when shifting them behind more than one command-line file.
- `isec -R FILE -r LIST` (or `-T` then `-t`): the later option wins and is read as what it is;
  bcftools keeps reading a later `-r`/`-t` string as a file name.
- `isec` does not parse the FORMAT columns of VCF lines whose records it never writes: those of
  every input for the sites list, and of inputs left out by `-w` or `-C`, unless `-i/-e` applies
  to them. So their sample data draws no warnings (undefined FORMAT tags) or errors (malformed
  values, where bcftools stops with a VCF parse error).
- The synced reader (`pasta::variant::synced`): when a record queued for output at a position is
  no longer buffered (HTSlib asserts), that reader has no record there; `bcf_sr_set_samples`
  with a `^` exclusion list keeps the first reader's other samples (HTSlib reads the list with
  the `^` and keeps no sample);
  `seek` returns the number of readers positioned on the sequence (HTSlib's returns 0). In
  `regidx`, merging overlapping regions keeps each payload with its region (bcftools moves the
  regions but not the payloads).
- An output (`-o`) that is one of the inputs (the same device and inode) is refused before it is
  opened, with status 255 and the input left untouched; bcftools truncates the input first.
- `concat -l` fails (status 255) on a chunk whose first record cannot be read (a corrupt or
  truncated BGZF block); bcftools skips it as an empty chunk and exits 0.
- Inputs read in step (`concat -a`/`-l`, `isec`, `merge`, `annotate -a`) share the `--threads`
  decompression threads, each getting `threads / N`; bcftools shares one pool of `threads`.

Filter expressions (`-i/-e`) and `filter` differ from bcftools where it has bugs, reads memory
it should not, or crashes:

- `ID!~"regex"` is true when no ID matches; bcftools does not negate `!~` for ID (it tests
  `!=` twice), so `!~` behaves as `~`.
- `ALT[0,1]` gives exactly the listed alleles; bcftools takes the number of listed indices as
  the last index, adding or dropping alleles when the list has gaps.
- A subscript past the end of a value reads nothing (or a missing value): `STR[0,2]` of a
  two-item string list (bcftools reads a previous record's string past the end),
  `FMT/TAG[*:3]` with fewer values per sample, and `binom()`/`fisher()` with alleles past the
  values (bcftools reads the next sample's values).
- Missing integer values (`.`) are missing: `AC<0`, `SUM(AC)<0`, and `AC[*]<0` are false and
  `AC="."` is true for `AC=.` or `AC=.,1`, `FMT/AD[*:GT]<0` is false for `AD=.`, and `query -f
  '%SUM(AC)'` prints `.` for `AC=.`, as for floats and the other FORMAT subscripts. bcftools
  copies them as the numbers -128 or -2147483648 for INFO tags that are not Number=1 and for
  FORMAT tags with `[GT]` subscripts, also into `binom()` and `fisher()`. An index past the only
  value of such an INFO tag reads nothing, as past several (`AC[1]` of `AC=5`, `AC[0]` of
  `AC=.`); bcftools reads the one value at any index without the missing check.
- `SMPL_*()`, `binom()`, and `fisher()` of a FORMAT tag absent from the record (or of
  arithmetic on it that has no values) give missing values, and `INFO/STR=@file` compares an
  empty string; bcftools reuses the previous record's values. Arithmetic between a tag with a
  sample subscript (`FMT/AD[1:*]`) and one over all samples reads 0 for the unselected samples,
  where bcftools reads values it has not set.
- `ABS()` and `PHRED()` of `N_MISSING`, `F_MISSING`, or `PHRED()` are computed; bcftools
  leaves them unset.
- `TAG[*:GT]` subscripts use the record's own genotypes. bcftools caches them by the record's
  address in memory and does not refresh them at records without the tag, so a later record at
  the same address can use another record's genotypes.
- `5!=AD[*:GT]` gives the same result as `AD[*:GT]!=5`. With the `[GT]` values on the right,
  bcftools checks the left side for their early end and compares the padding of samples with
  fewer alleles as missing values.
- `filter -s +` names the new filter `Filter2` when `Filter1` is taken; bcftools appends the
  names it tries (`Filter1Filter2`).
- A modulo by zero, and the computed `AC`, `AF`, `MAC`, and `MAF` of a record without INFO/AC
  and AN in a file without samples, are missing, and a function without arguments (`MAX()`) or
  a sample subscript (`FMT/DP[1]`) in a file without samples is an error; bcftools crashes.
- Regexes are POSIX basic regular expressions as in bcftools (with GNU's `\+`, `\?`, `\|`), run
  by the `regex` crate. Back-references (`\1` to `\9`) match as glibc's `regexec` matches them,
  on a backtracking matcher that tries each state (a place in the pattern, a position, and the
  groups a later reference reads) once, so its time grows polynomially with the text and its
  stack is on the heap; a reference to a group not closed before it is refused, as in bcftools.
  glibc has bugs with back-references, where vtools matches what the pattern says: it misses
  some matches of a reference to a repeated group (`^\(a\)\{0,2\}\1` on `aa`), it matches a
  reference to a group that took no part in the match where the match starts (`\(\(a*\)b\)*\2$`
  matches every text, at its end), it crashes on some nested repetitions
  (`\(\(a*\)\2\{1,2\}\2\)*.*`), and its time grows with the cube of the length on others or
  worse (`^\(a*\)*\1x$` on 1,000 `a`s and `bx`: 38 s, vtools 0.2 s). On such patterns vtools'
  memory grows with the square of the length (about 430 MB for a 3 kb text). glibc also ignores
  some anchors (`^`, `$`, `\b`, `\B`, `\<`, `\>`, `` \` ``, `\'`) of the copies it makes of a
  group it repeats with `\{m,n\}` or `\+`: one inside the second or a later copy that does not
  end it, and, with `\{m,n\}`, most of those that end a copy another copy follows (`\(^.\)\{2\}`
  matches `ab`, `\(.\>\)\{2\}$` matches `ab`, `\(\>.\)\{2\}$` matches `a.b`). Its answer for a
  pattern with back-references can also depend on the texts it was given before (`\(.\)\1\>`
  matches `**` when `aa` came first), so bcftools' results can depend on record order. vtools
  matches what the pattern says, for each text alone.
- `filter -g/-G` holds back as many records as the gap filters need. bcftools holds at most 64
  and then silently drops the oldest (for example, a long deletion with many SNPs inside it).
- `filter -S .` and `-S 0` keep a missing INFO/AN, and each missing INFO/AC value (`AC=.,2`),
  missing when they take away the alleles of the genotypes they reset. bcftools subtracts from
  the missing value's integer, so `AC=.` becomes `AC=2147483646` and `AN=.` a number such as
  2147483644 or -2147483646.
- `concat` checks that each chromosome block is sorted, as its "not sorted, consider running
  with -a" error intends; bcftools never updates the previous position, so it concatenates
  overlapping or out-of-order files on one chromosome into unsorted output.
- `concat -l` compares the phase of the overlap by the alleles the genotype indexes name, so two
  chunks that list a site's ALT alleles in different orders (`C,G` with `1|2`, `G,C` with `2|1`)
  agree. bcftools compares the indexes, counts such sites as phase mismatches, and can flip the
  phase of everything after the overlap.
- Plain `concat` (without `-a`, `-l`, or `--naive`) fails with `Error: BCF read error: FILE`
  (status 255) at a truncated or malformed BCF record, or at a truncated or corrupt BGZF block
  of BCF input, after writing the records before it in full (see Damaged input). bcftools takes
  it for the end of that file and goes on with the next one, dropping the rest of the file with
  status 0 (with only HTSlib's `[E::bgzf_read]` or `[W::bcf_record_check]` message, if any).
  `-a` and `-l` fail with "BCF read error" in both; `--naive` copies such records as they are in
  both.
- `concat --naive` decompresses the last BGZF block of each input before copying it. bcftools
  checks only each block's header and size, so a file cut inside its last block and then given
  an EOF marker whose 28 bytes make up the missing length is copied with its corrupt block, with
  status 0; vtools fails with `Error: the last bgzf block of FILE is truncated or corrupt`
  (status 255). Other cut blocks fail in both ("Could not read N bytes", "Could not parse the
  header of a bgzf block"). An input without the EOF marker block is copied and then fails (see
  Damaged input).
- `concat` copies VCF lines as they are when the output is VCF (bcftools' fast path) only when the
  output really is VCF. bcftools decides by `-O` alone, so `-o out.bcf` without `-O` writes VCF
  text into a BCF stream.
- `concat` and `reheader` pass BGZF blocks through without recompressing them where the data is
  unchanged (bcftools does this only for `reheader` of bgzipped VCF and `concat --naive`), so
  `-O b0`-`b9` levels apply only to recompressed parts, and `-Ob` output keeps the input's
  compression. The decompressed output is the same. A write that fails while `concat` copies
  blocks is reported as `Failed to write 0 bytes`, without HTSlib's `[E::bgzf_flush]` line;
  bcftools names the length of the VCF line it was writing.
- `concat` (every mode: plain, `--naive`, `-l`, `-a`) and `annotate -a` read their inputs twice.
  A FIFO is opened the second time without waiting for a writer: once its writer is gone it
  reads as empty ("Input is not detected as bcf or vcf format" and "Failed to parse header",
  "unknown file type", or the `--naive` error), where bcftools waits forever. A pipe is read on
  from where the first reading left it, as in bcftools, but the first readings take different
  amounts: a bgzipped pipe of roughly 64 to 200 KB is drained by bcftools (the empty-input
  messages) and not by vtools (`[E::hts_hopen]` and `[E::hts_open_format] ... Exec format
  error`). A tab-delimited `annotate -a` pipe or FIFO is not read again: without END columns it
  fails with "Could not initialize the annotation file", and with them nothing is annotated, as
  bcftools does when the type detection drained the pipe; bcftools reads what a larger pipe left
  ("Expected tabix-indexed annotation file", "Failed to parse", or annotations from the
  remaining lines).
- `tabix` indexes a FIFO, which it reads once; HTSlib's `tabix` opens it a second time and waits
  forever.
- `sort` orders records with the same position by REF and ALT as bytes in every case. bcftools
  compares them case-insensitively when merging temporary files, so records whose alleles differ
  only in case could come out in a different order depending on `-m`.
- `sort -m` is read and converted as bcftools does (`-m 1Mb` is one byte, a negative value
  wraps to a size that cannot be allocated, which fails with bcftools' message). bcftools
  allocates 90% of `-m` up front and fails when that is more than the system can give
  (`-m 1000G`); vtools uses memory as records arrive, so such a run sorts in memory.
- `sort -Ob`/`-Ou` of positions of 2^32 or more fails as HTSlib's `bcf_write` does ("64-bit
  values not representable in BCF"). bcftools sort reads its packed temporary records back
  without the VCF parser's 64-bit flag, so it writes BCF with such positions truncated to 32
  bits (2^32 + 4 becomes 4). Below 2^32 both write the records, which BCF's unsigned 32-bit POS
  holds.
- `norm -m-`, `-m+`, `-a`, and `-d` with `-Ob`/`-Ou`, at positions of 2^32 or more, fail as
  HTSlib's `bcf_write` does ("64-bit values not representable in BCF", exit 255). bcftools norm
  buffers copies of its records, which lose the VCF parser's 64-bit flag, so it writes BCF with
  such positions truncated to 32 bits (2^32 + 4 becomes 4). From 2^31 to 2^32 - 1 both write the
  records, which BCF's unsigned 32-bit POS holds.
- `isec` prints positions of 2^31 or more as they are in the sites list (standard output and
  `sites.txt`); bcftools prints them as 32-bit integers, so 2^31 becomes -2147483648 and 2^32
  becomes 0.
- `reheader -s` with "old new" pairs: a line that starts with whitespace has an empty old name.
  bcftools reuses the previous line's old name (a stale buffer), or crashes on the first line.
- `reheader -n/-s` of a header whose `#CHROM` line ends with FORMAT (no samples) adds the new
  samples after FORMAT; bcftools reads past the end of the header text looking for the samples
  and crashes (VCF), writes garbage (bgzipped VCF), or leaves the samples out (VCF with `-f`).
- `reheader -h` of a BCF file whose records are on a contig the new header lacks, while it has a
  contig with a later IDX, stops with "The CHROM is not defined"; bcftools compares the contig
  name with the missing one (a null pointer) and crashes.
- `reheader -o FILE` of bgzipped VCF input fails with the reason when FILE cannot be created;
  bcftools does not check that it opened the output and crashes.
- `reheader` of BCF, bgzipped or not, with a truncated or malformed record (cut short, a length
  that runs past the data or is shorter than the fixed fields, or a record HTSlib's
  `bcf_record_check` refuses) writes the records before it and then fails with
  `[reheader_bcf] Error: close failed .. FILE` (status 255), as it does for a truncated BGZF
  block; bcftools takes the bad record for the end of the file and exits with status 0. Bgzipped
  VCF is copied in whole blocks: a truncated block or one whose header does not parse (or
  garbage after the EOF block) is not copied, the blocks before it are, and `reheader` fails
  with HTSlib's `bgzf_read_block` error and `Error reading FILE` (status 255); bcftools copies
  the bytes as they are, which writes corrupt BGZF, with status 0. A block whose data is corrupt
  (a bad CRC32) is copied as it is by both, with status 0. Bgzipped VCF or BCF without the EOF
  marker block gets HTSlib's warning (bcftools prints it only for BCF) and then fails with
  status 255 after the output (see Damaged input).
- When the reader of the output goes away (`| head`, or a FIFO whose reader exits), vtools is
  killed by SIGPIPE with nothing more on standard error (status 141 in the shell), as bcftools
  is. `sort` removes its temporary directory and files before it dies; bcftools leaves them
  behind. When bcftools sort exits with `error()` rather than its clean-up (an index that cannot
  be made or written), vtools removes the temporary files too, without the "Cleaning" message.
  Started with SIGPIPE ignored or blocked by its parent (for example from Python's `os.system`),
  bcftools gets EPIPE, prints a write error (`Error: close failed stdout`) and exits with 255;
  vtools is still killed by SIGPIPE then.
- With standard output closed (`vtools view x.bcf >&-`), a command that writes to it fails with
  its write error and status 255, as bcftools does. vtools keeps descriptor 1 taken (by the root
  directory, open only for reading) so that no file it opens later is given it. bcftools gives
  descriptor 1 to its first input, so with `-o /dev/stdout` (or `/dev/fd/1`) it truncates its
  input and writes its output over it. vtools cannot open that output (`Is a directory`), fails
  with status 255 and leaves the input as it was.
- `head`, `stats`, `index -s/-n`, `query -l` and `+counts` report a failed write to standard
  output (a closed standard output, or a full disk) and exit with status 255. bcftools ignores
  it and exits with status 0, having written nothing. Usage, version and list output
  (`--version`, `-h`, `+split-vep -l`, plugin lists) still ignores a closed standard output, as
  in bcftools.
- `tabix` reads from the first indexed record for the region `.`; HTSlib does not seek when that
  record is at virtual offset 0 (a file without header lines) and continues from wherever an
  earlier region left the file.
- `merge` fills a Number=A FORMAT value of an allele a file lacks with the value of its unknown
  allele (`<*>`, `<NON_REF>`), as it does for Number=R and G. bcftools reads one value further,
  past the sample's values when the unknown allele is the last, and writes whatever is there.
- `merge` sets the length of a line whose INFO/END it copies from the input (symbolic alleles),
  so BCF output and `--write-index` indexes cover the whole record. bcftools keeps the REF
  length (the VCF text is the same).
- `merge --use-header` with a header whose sample count differs from the inputs' stops with an
  error; bcftools fails an assertion (and aborts).
- `merge` where bcftools reads memory it never set: a Number=A FORMAT field at a line without
  ALT alleles copies nothing (bcftools writes past its buffer), and `-L` at a site without
  FORMAT/PL where a file has more ALT alleles than allowed gives up localizing with bcftools'
  warning (bcftools reads a stray pointer, and warns when it is not NULL). With `--threads`,
  `merge` formats VCF output lines on the workers.
- `merge` never writes the same line forever. Where every record left at a position is kept
  from the line of the open gVCF blocks (`-m exact` with a different number of alleles, or
  `-m none`/`indels`/`exact` when the most frequent allele is a block's), bcftools writes
  the block line again and again; vtools merges the blocks with the first record left, in
  file order, as the other `-m` modes do (or stops with bcftools' error when the REFs
  differ, as with a block whose REF is N).
- `merge -i` rules over a Number=G tag given per allele (a REF-only site) in one file and per
  genotype in another read zeros where bcftools reads past the values (zeros the first time,
  values of earlier sites after that). A Number=A/G/R INFO tag merged into the buffer of
  another one (the tag is Number=. in the header of the file it was first copied from) is
  only written within that buffer, where bcftools writes past it.
- `merge` keeps records at positions above 2^31; bcftools keeps the current position in a 32-bit
  `int`, which wraps there (to -1, the value for no position, at 2^32), so the record is dropped
  (with `--force-single`, indexed or `--no-index` inputs alike). With `-Ob`/`-Ou`, a record at
  2^32 or more fails as HTSlib's `bcf_write` does ("64-bit values not representable in BCF",
  exit 255) rather than being written with its position truncated to 32 bits.
- `merge --no-index` with contigs that are not in the headers (declared on the fly while
  reading): bcftools looks the contig up in the output header, which lacks it, and crashes (a
  build whose check catches the missing contig first stops with `Uh` and status 255, as vtools
  does); vtools stops with bcftools' check message (`Uh`). Indexed inputs get the contigs from
  their indexes, as in bcftools.
- `consensus` reads only the samples it uses from VCF text input (none in REF,ALT mode) when no
  `-i/-e` expression is given and the `-H` mode does not depend on the other samples (`-H N`
  with N > 1 and `-H NpIu` look at the padding of genotypes shorter than the longest at the
  site, so they read every sample); bcftools parses every sample. So malformed columns of the
  unused samples are not reported, and a site whose GT bcftools stores in 16 bits because an
  unused sample has an allele index above 62 is applied where bcftools stops with its "Todo:
  GT field represented with BCF_BT_INT8" error. `--threads` (BGZF decompression) is a vtools
  addition.
- `consensus -a`: bcftools refills the whole unwritten part of the sequence with the absent
  character after every FASTA line, which is quadratic over long stretches without variants
  (37 s on the 8 Mb short-read reference); vtools fills each base once, with the same output.
- `consensus` writes the error for a REF that does not match the FASTA with the sequence up to
  the end of its buffer; bcftools prints past the end (stale or uninitialized bytes).
- `consensus` writes a stretch of over 61,440 bases with no record pending as it reads it, so
  when it stops with an error (such as a REF that does not match the FASTA) it has written more
  of the sequence than bcftools: bcftools' output is a prefix of vtools'. The exit status and
  messages are the same.
- `consensus` buffers 64 KB of output, so when the reader of a short output has already gone, it
  prints its warnings and summary before the write that kills it with SIGPIPE; bcftools writes
  its output in 4 KB pieces and dies at its first write, before printing them. The exit status
  is the same.
- `consensus` never reads or writes before the start of its sequence buffer. bcftools can write
  there for a variant whose first base is not shared by REF and ALT, that starts before a
  `>chr:from-to` region, and whose ALT is shorter than the overlap (undefined behaviour in C).
  Nor does it read outside the buffer for a missing (`-M`) reference block (a gVCF `./.` block):
  bcftools copies every base of the block from its buffer, from before it when the block starts
  before a `>chr:from-to` region and past it when the block ends after the region or the
  sequence, and crashes (segfault) or fails on what it finds there ("Requested allele outside
  valid range"). vtools clips the block to the part of the region or sequence it covers and
  masks just that, without bcftools' "trimming variant" warning, which would print the whole
  block. It copies only those bases, so a long block uses less memory than in bcftools (125 MB
  against 206 MB for a 40 Mb block).
- `consensus` with a `>chr:from-to` header with a negative coordinate treats it as 0.
- `consensus` with a variant that starts past the end of its FASTA sequence (a position past
  the reference, or past 2^31): bcftools computes a negative length, prints "trimming
  variant", writes before the allele's buffer, and aborts on the heap corruption ("free():
  invalid next size") without writing the rest of its output; vtools prints the same warning,
  skips the variant, and goes on.
- `consensus` with a symbolic allele such as `<INS>` at the position where an applied variant
  ends: HTSlib does not set the length of a symbolic allele, so bcftools decides whether the
  site overlaps with a length left over from an earlier record, and either skips it or stops
  with its error about unsupported symbolic alleles; vtools takes the length as 0 and skips the
  site as overlapping ("The site ... overlaps with another variant, skipping...").

`call` matches bcftools' numbers exactly: the same floating-point operations in the same order
(including `float` intermediates), glibc's `lgamma` for Fisher's test, and the 80-bit x87
`long double` sums of `prob1.c` behind `-c`'s QUAL, AF1, and FQ (on the x87 unit on x86-64, in a
software emulation elsewhere). I16 values are converted to integers as bcftools' C casts convert
them on x86: a `.` (missing) or a value of 2^31 or more becomes `INT_MIN`, HTSlib's missing
integer, and sums of such counts wrap as C's `int` does. Such a read count is written as `.` in
DP4. PV4 is left out when the read counts hold a `.`, and written, from the wrapped counts, when
they hold a value past 2^31. With `-c`, MQ is `.` when the read counts hold either, or the
squared mapping-quality sums (I16[9], I16[11]) a `.`; with `-m`, MQ is 0 when the read counts
hold either, and `.` when the mapping-quality sums (I16[8], I16[10]) hold a `.`. Its deliberate
differences:

- `-t REGION` that is not a file: bcftools opens the value as a file, fails, and calls every site;
  vtools calls the sites that start in the regions (`-t ^REGION` calls the others). With
  `-C alleles`, bcftools crashes on such a value; vtools stops with "Expected -t or -T". (`-T
  FILE` and `-t FILE` read `CHROM POS [REF,ALT]` lines as bcftools does.)
- `-G FILE`: bcftools keys the groups on their names without the first character, so groups
  `A` and `B` are one group; vtools uses the whole name.
- `-s`/`-S` lists with lines that are not sample names as they stand (names missing from the
  input, `#` comment lines, leading spaces): bcftools subsets the header with the first lines of
  the list as they are and the records with the samples found, so it stops with `No such sample`
  or crashes when no sample is left; vtools uses the samples found, and reports `No matching
  sample found` when there is none.
- `-c` on records without INFO/I16 (such as `call` output): bcftools reuses the previous
  record's I16 values, and a record with fewer than 16 values replaces only the first ones; it
  crashes when the first record has none, and reads memory it never set when the first has fewer
  than 16. vtools does the same (also with `--threads`, which carries the values from one batch
  of records to the next) but uses zeros for values no record has set: zero counts (DP4=0,0,0,0,
  MQ=0) for a first record without I16.
- `-c` with PL values above 255 (`mpileup` caps them at 255), and `-m` filling a missing PL from
  an unseen-allele PL above 255: bcftools reads past the end of its 256-entry probability table;
  vtools converts them as `-m` converts other values.
- `-c` with missing PLs among a sample's first three values: bcftools uses them as indexes into
  its probability table and crashes; vtools reads them as `-m` does: a sample whose first PL is
  missing, or whose PLs end early, has no data (as `mpileup` writes 0,0,0), and another missing
  PL is 255.
- PLs with one value per allele for every sample (haploid only, which `mpileup` never writes):
  bcftools reads them with the diploid stride, past the end of the record's values into what
  earlier records left in its buffer (and with `-g` it can stop with an error about another
  record); vtools reads the values past the end as missing (`-m`) or as 0 (`-c`).
- `-m` at a site whose best alleles include the unseen allele `<*>` that is not kept: bcftools
  writes genotypes that point at it although it is no longer listed (`1/1` with ALT `.`, `1/2`
  with one ALT), with GP from values it never set (NaN, or for a haploid sample other samples'
  values). vtools writes that allele missing (`./.`, `0/.`, `1/.`, `.`), counts only the alleles
  it writes in AC and AN, and writes GQ and GP as for a sample without data.
- `-g` with `-A`, when a block starts at a site whose PL was removed and continues with sites
  that keep it: bcftools reads a buffer it never filled (and crashes); the block keeps no PL.
- Number=R Float fields of trimmed alleles stay Float; bcftools writes their bits as integers
  (`mpileup` writes no such fields).
- `-C trio` is disabled in bcftools 1.24 (it stops with `todo: constrained trio calling
  temporarily disabled` at the first variant site); vtools does the same.
- `-C alleles` at a site without FORMAT/PL (or without samples) whose alleles change: bcftools
  reads the PL values it could not get (and crashes, or divides by zero samples); vtools goes on
  to `-m`'s checks, which stop with its error for such sites.
- `-C alleles -i` with a record at POS 1 (the first record of a sequence, or another one there):
  bcftools writes the targets with no record up to the position before each record, passed as an
  unsigned `POS-1`, which wraps at POS 1, so it writes every later target of the sequence as
  missing, before that record, and then drops the calls at those targets; vtools writes nothing
  before POS 1.
- A sample with a ploidy neither caller genotypes (above 2, or negative; ploidies are bytes, so
  258 is diploid as in bcftools), from `--ploidy-file`, `--ploidy`, or `-S`: bcftools `-m`
  writes nonsense there (AN counts the ploidy, the genotypes are haploid, PL is missing) and,
  built with assertions, aborts at the first variant site; `-c` computes nonsense. vtools writes
  the calls before it and stops with an error (exit status 255) at the first record where a
  sample has such a ploidy, with the ploidies set as bcftools sets them: samples given a ploidy
  in `-S` have the largest one only until the ploidy of a sex first changes, and a region on a
  sequence the input lacks never applies, so a file whose unsupported ploidy no sample gets at a
  called record is called as bcftools calls it. `-c` sizes its tables by the largest ploidy,
  which must be 1 or 2: bcftools crashes otherwise (also when no sample has it), and vtools
  stops before calling.
- With `--threads`, vtools parses, calls, and formats on the workers (bcftools only compresses
  output); `-C alleles`, `-i`, `-r`/`-R`, `-t REGION`, and per-sample ploidies given as numbers in
  `-S` run serially. An unreadable record ends the input quietly as in a serial run: the records
  HTSlib's synced reader has read ahead without returning (those at the last position read, and
  one more at a position shared by two sequences) are not called or written.

`gtcheck` sums every pair's scores in the same order as bcftools, so the report is the same
(the sums of a batch of sites are added a block of pairs at a time, with AVX2 where the CPU has
it, and with `--threads` on several threads, without changing any pair's order of additions).
The first lines name vtools and its command line where bcftools names itself. Its deliberate
differences:

- An unprefixed `-e EXPR` (not a number) also sets bcftools' error probability to 0, as if
  `-E 0` were given: bcftools parses it as the old `-e INT` first and keeps the result. vtools
  keeps `-E` (40 by default) and uses EXPR as the filter.
- `-p/-P` with `--no-HWE-prob`, and `--n-matches` with a negative value and `--no-HWE-prob`:
  bcftools reads the match counts it did not allocate (and crashes); vtools prints the HWE
  columns as 0, as in the other modes, and ranks the matches of a negative `--n-matches` as tied.
- `-g FILE` whose header lists the contigs in another order than the query's: bcftools compares
  the two headers' contig numbers, so it skips sites both files have as `POS+ALT mismatch` and
  can compare a query site with the last `-g` record of another contig that has the same number
  in the other header. vtools matches the contigs by name.
- The DS section of `--distinctive-sites` starts with its comment lines; bcftools builds them
  but clears its buffer before writing them.
- `--distinctive-sites NUM,MEM,TMP` holds the sites in memory up to MEM (counted as bcftools
  counts it, a pointer and a record per site) and past it writes them in sorted runs to a
  temporary file named from TMP. bcftools writes each run to a file of its own and merges them
  all at once, holding a record and an open file for each, so its memory grows with the number
  of runs; vtools merges at most as many runs at a time as MEM holds records (as many as 1 MB
  holds when MEM is less, and from 2 to 1024 runs), in passes through a second temporary file
  when there are more, so the merge takes about MEM however many runs there are. The report is
  the same: sites with the same number of pairs are ordered by the same `lrand48` sequence (seed
  0) as in bcftools, and then in the order they came, as bcftools orders them when they fit in
  MEM (past MEM, bcftools' merge can swap two sites that also tie on the random number). vtools
  makes the files only past MEM and removes them as soon as they are made; bcftools makes one
  even when every site fits in MEM, so it fails with a TMP it cannot write. When TMP is given,
  bcftools reads MEM without its unit (`100M,/scratch` as 100 bytes); vtools keeps the unit.
- With one file, a site where the query's tag falls back to the other one (PL asked for but
  absent from the record, say) reads the second sample of each pair with that tag too (the INFO
  counters still say what bcftools says); bcftools reads GT values as PLs there, past the end of
  the array for the last sample.
- `-O z` writes gzip (compressed differently); the text is the same.

`annotate`'s deliberate differences:

- `-i/-e` with typed external values (`{int:NAME}`, `{str:NAME}`, `{float:NAME}`): bcftools takes
  the type as the column name and stops with `No such column: int`; vtools uses NAME.
- `-i/-e` with `{NAME}` values when `-c` lists no annotation columns (only `CHROM,POS`, say):
  bcftools reads the column index from uninitialised memory (and may crash); vtools stops with
  `No such column: NAME`.
- `-i/-e` with a quoted string before a `{}` value: bcftools rescans the closing quote as an
  opening one and stops with `Could not parse the expression`; vtools skips quoted text only.
- Number=A/R INFO values given as a single `.` for a multiallelic annotation (for Number=R
  strings, also a biallelic one): for integers and floats bcftools' count check reads past the
  values or can never pass, so it usually stops with `Incorrect number of values`; for strings
  it copies the values allele by allele, finds no second value in `.`, and stops with
  `Failed to copy a string field`, from a VCF, a BCF, or a tab file. vtools takes the `.` as
  missing for every allele, which gives what bcftools writes for the `.` spelled out per allele
  (`.,.`).
- FORMAT Number=A/R/G values given as a single `.` at a site whose alleles differ from the
  record's: in every sample of the annotation VCF (`-c FMT/AD` with `AD=.`, say), or, for
  Number=G, in a sample next to haploid ones (`PL=.` beside `PL=1,2,3`), which is then read as
  diploid from a row of haploid width. bcftools reads one value per source allele or genotype
  past the sample's values and writes whatever memory holds (the next sample's values, or
  garbage after the last sample); vtools writes missing values for those alleles or genotypes.
  In the second case the vector ends of the row still cut the sample's values short, in both
  tools.
- Number=G FORMAT values (`-c FMT/PL`) onto a record whose values are of haploid width (a
  haploid sample next to one with `PL=.`) from a diploid source sample: bcftools writes past the
  record's values and aborts (`double free or corruption`, exit 134); vtools keeps the record's
  width and writes the genotypes that fit in it.
- FORMAT Number=A/R/G transfers from a VCF whose samples are ordered differently or only partly
  shared: bcftools swaps the source and destination sample indexes (and sizes the output by the
  source's sample count); vtools maps each destination sample to its source sample. Whether a
  Number=A field's alleles are unchanged is decided from the ALT map alone; bcftools reads one
  entry past it.
- `-l TAG:sum` (and other merges) on a Number=A/R tag: bcftools aborts on an assertion; vtools
  reports an error.
- `-x FILTER` with FILTER definitions added back (`-a file.vcf -c FILTER`, `-h`, `-H`): bcftools
  frees the removed `##FILTER` lines but keeps pointers to them, so a definition with the same
  ID is silently not added and the output has FILTER values without definitions; vtools adds
  them.
- `-x FILTER` together with `-x FILTER/NAME` or `^FILTER/NAME`: bcftools sets FILTER from a null
  list of one ID and crashes; vtools empties FILTER. An empty `-x` item (`-x ,ID`) is ignored;
  bcftools calls an uninitialised handler.
- `--set-id` with `%POS` at positions from 2^31 on: bcftools prints the position as a 32-bit
  `int` (negative, or 0 at 2^32); vtools prints it as it is.
- `-a FILE.tab` at a record at POS 0: HTSlib's region reader reports an overlap before it has
  read an annotation line, and bcftools parses that empty line and crashes; vtools reads the
  annotation lines as for any other position (none overlap POS 0 in practice).
- `-a FILE.tab` at positions above 2^31: bcftools reads the annotation's position into a 32-bit
  `int`, which wraps, so it does not annotate those records, and at POS 2^32 it crashes; vtools
  annotates them.
- Columns past the end of a tab-delimited line read as empty strings, and `-c` with `REF` but no
  `ALT` matches no allele; bcftools reads unrelated memory in both cases.
- FORMAT/GT copied from an annotation VCF or BCF (`-c FMT/GT`, `-c FORMAT`, `FMT/NEW:=GT`) whose
  record lists the alleles in another order, or has others: bcftools copies the allele indexes
  as they are, so they name other alleles (ALT `C,G` with GT `1/1`, C/C, onto ALT `G,C` gives
  `1/1`, G/G, and an index past the record's ALTs can be written). vtools rewrites each index to
  the same allele among the alleles the record has once every column is set. Alleles are matched
  as for Number=A/R values, so a longer REF on either side still matches. An allele the record
  lacks becomes missing (`.`), as a Number=A/R value of an allele the annotation lacks does.
  Every allele does when the REF alleles are not compatible (`--pair-logic any`, say). Phasing,
  ploidy and missing alleles are kept. When `-c` also takes the ALT alleles from the annotation
  (`ALT`, whether or not `+ALT` is listed too, or `+ALT` onto a record without ALTs), the record
  ends up with the annotation's ALTs, so the indexes are copied as they are, as bcftools copies
  them. This holds whichever column comes first in the list. With `REF` but not `ALT`, the
  genotypes are mapped onto the annotation's REF and the record's ALTs. Where the annotation's
  alleles are the record's, or its first ones, nothing changes.

`roh` matches bcftools' numbers exactly (the HMM's operations in `HMM.c`'s order, with its
precalculated powers of the transition matrix). Its output names vtools in the first two lines.
`--threads` runs the samples' HMMs on workers, with the same output. Deliberate differences:

- A malformed record (an unparsable VCF line, a malformed BCF record) or damaged input (see
  Damaged input) stops `roh` with the reader's error and status 255, before it writes any
  results. bcftools takes the failed read for the end of the input, writes the results of the
  sites before it, and exits 0.
- `-b` with samples whose sites differ (as with `-i`): bcftools' model keeps one snapshot
  pointer for all samples, so a sample's last run of a sequence can overwrite the snapshot that
  another sample's next buffer starts from. vtools keeps each sample's snapshot apart, which
  gives the results bcftools gives for each sample run alone.
- `-V` with a genetic map per chromosome (`-m` with `{CHROM}`) and more than one chromosome:
  bcftools keeps one map, replaced at each new chromosome, and with `-V` runs the HMM only at
  the end of the input, so it trains and decodes every chromosome with the map of the last one
  (or with none, when the last chromosome's map is missing and it is skipped). vtools trains and
  decodes each chromosome with its own map.
- A missing INFO value (`.`) for `--AF-tag`: bcftools computes with it (its NaN), which turns
  every later state of the sequence into `nan`; vtools takes the frequency as unknown, as for a
  site without the tag (`--AF-dflt` applies, or the site counts as no AF).
- NaN for `-a/--hw-to-az` or `-H/--az-to-hw`, or as an `--AF-file` frequency: bcftools accepts
  it (and prints `nan` qualities); vtools rejects it as out of range, with bcftools' messages.
- Positions from 2^31 on: bcftools keeps the previous position in a 32-bit `int`, so it misses
  duplicate positions there, and prints positions and region lengths with `%d` (negative);
  vtools finds the duplicates and prints the positions unsigned. As in bcftools, positions are
  kept in 32 bits.
- `-e -` with `-s/-S` on VCF input: bcftools drops the samples not in `-s` while parsing (a
  speed-up), so it estimates AF from the `-s` samples only; on BCF input it uses all samples, as
  documented. vtools uses all samples for both.
- `-s/-S` and `-e FILE` naming the same sample, on VCF input: bcftools crashes (its sample list
  holds the duplicate); vtools runs.
- `--include/--exclude` reading FORMAT fields with `-s/-S` on VCF input: bcftools compiles the
  expression for all samples and evaluates it on records holding only the `-s` samples (wrong
  results); vtools keeps all samples then, giving bcftools' results for BCF input.
- `-G` at a site without diploid-width GTs (every sample haploid, a sample of ploidy above 2, or
  no GT field): bcftools crashes; vtools skips the site.
- `-G` and `-e GT` (the default for `-e`) at a site where haploid GTs sit next to diploid ones
  (male chrX, say): a haploid GT is padded with the vector end, which bcftools reads as an
  allele, so `-G` takes a haploid `0` or `1` for a het (with likelihoods favouring HW) and `-e
  GT` counts the padding as an ALT allele, which shifts the frequency used for every sample at
  the site. vtools skips the haploid samples for `-G` (a haploid call says nothing about
  autozygosity; without `-G` their PLs are skipped too) and counts a haploid call as one allele
  for `-e GT`, as AC/AN do.
- `--AF-tag` with fewer values than the index of the ALT allele, and `--AF-tag` with
  `--include-noalt` at a site without ALT: bcftools reads outside the values; vtools takes the
  frequency as unknown.
- A malformed record: bcftools ends its loop there, prints the results so far, and exits with 0;
  vtools stops with the error (exit status 255).
- `-s` and `-S` together (or `-r`/`-R`, `-t`/`-T`): the later one is used; bcftools uses the
  later value but reads it as a file if the file form was given anywhere.
- The usage (`-h`) leaves out bcftools' `roh-viz` example: vtools has no such script.

Command lines (all commands; see the getopt paragraph under Commands). Deliberate differences:

- The usage a command prints after a bad option, and for `-h`/`--help`, is vtools' option list
  (for commands whose own usage text vtools has, that text); the exit status and stream are
  bcftools'. `vtools` with no arguments or `help` lists vtools' commands and its version.
- Every command takes `--help`. Where bcftools has no `--help` (it is an unknown option there,
  which also ends in the usage), vtools matches it only in full, so `view --he` is still
  `--header-only` as in bcftools.
- vtools' `concat` also takes `--output-file` (as bcftools' `view`, `index`, `query`, and `sort`
  do), so `--output-` is ambiguous there; bcftools takes it as `--output-type`.
- An argument that clap still rejects after the getopt rewrite (a value that is not UTF-8)
  stops with clap's message and exit status 255, not clap's 2.
- `--threads` with a negative count in `view`, `norm` ("Failed to create threads") and `convert`
  ("Could not initialize --threads N") fails as in bcftools; `filter`, and the vtools additions
  `query --threads` and `consensus --threads`, run without threads, as HTSlib's
  `hts_set_threads` ignores it.
- `reheader -o` into a directory that does not exist fails with the error (exit status 255);
  bcftools 1.24 crashes (segmentation fault) after failing to open the output.
- vtools reads the whole command line before it checks any option's value, so when a command line
  has two problems, getopt's (an unknown option, a missing value, an ambiguous prefix) is
  reported first. bcftools checks each value as its getopt loop reads the option, so a bad value
  before the getopt problem wins there (`stats -i A -e B --nosuch` fails with "only one -i or -e
  expression" and status 255 in bcftools, with the usage and status 1 in vtools). Likewise an
  option given more than once is checked only for its last value (`sort -v x -v 1` runs; bcftools
  rejects `-v x`).
- `--threads` is read as bcftools reads it (`int n_threads = strtol(optarg, 0, 0)`), so a count
  that does not fit an `int` wraps as in C (`4294967297` is one thread). Each pool then starts
  at most 1024 threads, and under an address-space limit (`ulimit -v`) at most one per 64 MB of
  it; under such a limit glibc's malloc arenas are also capped at a quarter of it (at most eight
  per online CPU, glibc's own cap), unless `MALLOC_ARENA_MAX` or `glibc.malloc.arena_max` in
  `GLIBC_TUNABLES` sets them. The output does not depend on the count. When the system cannot
  start all the threads of a pool (`ulimit -u`, `ulimit -v`), vtools runs it with half of those
  it could start, leaving room for the next; a reader whose own thread cannot start reads on the
  calling thread. HTSlib has no limit: bcftools tries to start every thread asked for, so a huge
  count (`--threads 99999999999`, 1215752191 as an `int`) hangs it until the system runs out of
  threads or memory, or, under a memory limit, fails with "Failed to create threads" (status
  255; `convert`: "Could not initialize --threads N"; `csq`: "Failed to create N extra threads";
  `reheader` of BCF: "Could not initialize threading"). vtools runs those commands with the
  threads it can start. When not even one thread can start, vtools fails with those messages
  where bcftools does (`view`, `annotate`, `norm`, `merge`, `stats`, `roh`, `concat`, `csq`,
  `convert`, `reheader`), and `filter`, `call`, `isec`, `mpileup` and the plugins run without
  threads, as in bcftools (`mpileup` reads each input on a thread of its own, or on its main
  thread when that thread cannot start). When a thread can start but a command's workers,
  started after its reader's and writer's threads, then cannot (those threads, or other
  processes under `ulimit -u`, took the room), the command runs without them, with the same
  output (`stats`, `call`, `mpileup`, and the hap file of `convert --hapsample`, compressed on
  the main thread); `concat` checks that a thread starts before its output's compression threads
  start, as bcftools creates its thread pool before giving it to the output. Counts between the
  number of cores and the cap run, but slowly: vtools' worker pools are not meant for more
  threads than the machine has.
- A list file (regions, targets, masks) whose name ends in `.fai` and whose first line has a name
  and four integer columns is read as text. HTSlib detects such a file as a FASTA index and
  refuses to open it ("Exec format error"), so bcftools fails there.

Plugins (`plugin`, `+NAME`) are built in: `-l` lists `counts`, `fill-tags`, `setGT`, and
`split-vep`, an unknown name fails as a plugin that cannot be loaded (status 255), and `-V`
prints vtools' version. `--threads` is also accepted by `+split-vep` (bcftools' has no such
option). Deliberate differences, where bcftools has bugs:

- `-v`/`--verbosity` only sets the verbosity. bcftools' option handler falls through to `-o`, so
  `-v2` writes the output to a file named `2` (and `-v` makes the output name null).
- `+fill-tags -t INFO/AC_Hemi` is accepted; bcftools accepts `INFO_Hemi` (a typo) instead, which
  vtools also accepts.
- `+fill-tags -t FORMAT/TAG...=EXPR` for a FORMAT tag that already exists compares the
  definition with the FORMAT tag's type; bcftools compares it with the type of the INFO tag of
  the same name, so it fails ("Type definition is different") unless an INFO tag of that name
  and type exists.
- `+fill-tags -S`: a sample name of one or two characters is accepted; bcftools rejects the line
  ("Could not parse the file").
- A sample without a GT value (a sample column with fewer fields than FORMAT, such as `4` for
  `DP:GT`, which HTSlib fills with the integer missing value): `+fill-tags` counts it as a
  missing genotype; bcftools reads the value as allele -65 and writes outside its bit set
  (usually crashing). `+setGT -n m/M` stops with HTSlib's `Incorrect allele` error there, where
  bcftools' allele counting writes outside its array.
- `+setGT -t q -n i` swaps the alleles of each selected sample; bcftools swaps those of the first
  sample once per selected sample.
- `+setGT -n X` (and `c:` with `X`) when a sample's first FORMAT/AD value is missing: the allele
  with the most reads among the others is used; bcftools compares with the value before the
  sample's (the previous sample's last value, or memory before the array).
- `+setGT -n c:GT` on records without FORMAT/GT gives every sample the custom genotype (as
  bcftools does for one sample; with more samples bcftools reads and writes outside its array).
  `+counts` and `+setGT` run without an input from a terminal print the usage; bcftools crashes.
- When both `-r` and `-R` (or `-t` and `-T`) are given, the later one is used, inline list or
  file, as in vtools' other commands (also in `+split-vep`); bcftools takes the later value but
  reads it as a file name once `-R` (`-T`) was given anywhere, so `-R FILE -r REGION` fails.
- `+counts` reads only CHROM to ALT of each VCF line (the rest is not parsed unless `-i/-e`
  needs it), so a malformed QUAL, FILTER, INFO, or sample column no longer stops it; bcftools
  parses whole lines. This makes it several times faster on VCF. (With or without `--threads`.)
- `+split-vep -o` with a name that cannot be created: bcftools does not check the open and
  crashes; vtools stops with the error.
- `+setGT -n c:GT` with more alleles than the genotypes of the record (such as `c:0/0/1` on
  diploid genotypes): bcftools widens the genotype array in place, copying each sample forward
  over values it has not moved yet, so every sample but the first reads its first allele twice
  (`0/1` becomes `0/0`, `1/2` becomes `1/1`) before the targets are tested; vtools widens the
  genotypes intact. A custom allele index that is negative or does not fit in 32 bits is set
  as a missing allele; bcftools reads -1 to -4 as its own codes for `m`, `M`, `X`, and `.`
  (without working out the minor or major allele), wraps the larger values, and writes other
  negative indexes as invalid alleles such as `-5`.
- `+split-vep -c` with an index range that ends before it starts (`3-1`, or `2-` with the end
  left out): bcftools counts a negative number of columns (no columns, or a crash for `-2`);
  vtools stops with `No such column`.
- `+fill-tags -t TAG:N=EXPR` with a negative N: bcftools declares `Number=-1` and writes values
  from memory it never set; vtools declares the same and writes no value.
- `+fill-tags` with polyploid genotypes (more than two alleles, such as `0/0/1`) counts every
  called allele in AN, AC, AF, MAF, AC_Het, AC_Hom, AC_Hemi, and the `-d` half calls, as
  bcftools' own `view -s` recalculation of AC and AN does; bcftools counts each distinct allele
  of a genotype once (twice for a hom), so `0/0/1` counts as `0/1` and `1/1/1` as `1/1`. HWE and
  ExcHet, diploid tests, leave polyploid genotypes out, as bcftools already does for haploid and
  partly missing diploid ones; bcftools takes them in as those diploid genotypes and stops with
  "Expected diploid genotypes" when that makes the allele total odd. Diploid and haploid
  genotypes are counted as in bcftools.

## Performance

Measured on an Intel Xeon Gold 5418Y VM (24 vCPUs, 122 GB) on 2026-11-04 against bcftools/HTSlib
1.24, both with libdeflate at the same compression levels (these are vtools 0.2.0's timings; vtools
now compresses with zflate, at level 6 unless a level is given), as wall time: medians of 3 warm
runs, and one run each on the cohort. Commands are single-threaded unless the row says `--threads`.
`kg` is 1000 Genomes chr22 with 500 samples (1.1 M records), `giab` the one-sample GIAB HG002 calls
(4.0 M records), and the cohort rows use chr22 of the 3,202-sample 1000 Genomes 30x panel.

| Operation | bcftools (s) | vtools (s) | Speedup |
|---|---:|---:|---:|
| view VCF.gz to VCF, kg | 20.520 | 7.710 | 2.66x |
| view VCF.gz to VCF, kg, --threads 4 | 19.570 | 1.600 | 12.23x |
| view VCF.gz to BCF, kg, --threads 4 | 11.390 | 1.900 | 5.99x |
| view -f PASS -i 'FMT/DP>40', giab | 4.290 | 4.500 | 0.95x |
| query -f '%CHROM %POS[ %GT]', kg | 13.810 | 6.870 | 2.01x |
| filter -s LowAF -e 'AF<0.01', kg | 9.160 | 4.710 | 1.94x |
| index -t, kg vcf.gz | 3.950 | 1.590 | 2.48x |
| concat 4 chunks -Ob, kg | 8.240 | 0.690 | 11.94x |
| sort -Ob, giab | 27.890 | 25.540 | 1.09x |
| norm -m-any -Ob, kg | 9.150 | 8.270 | 1.11x |
| stats -s -, kg | 6.320 | 3.040 | 2.08x |
| mpileup -f -a AD,DP, short reads 1 Mb, --threads 4 | 3.310 | 0.880 | 3.76x |
| call -m, one-sample mpileup (1 M records) | 2.050 | 1.960 | 1.05x |
| view VCF.gz to BCF, cohort, --threads 4 | 67.290 | 13.220 | 5.09x |
| merge sample halves -Ob, cohort, --threads 4 | 72.440 | 29.850 | 2.43x |
| concat --ligate 2 chunks -Ob, cohort | 88.810 | 89.400 | 0.99x |

### Where bcftools is faster

- `view -f PASS -i 'FMT/DP>40'` on giab (0.95x): with one sample each record carries little work,
  and vtools' filter evaluation costs about 50 ns per record more than bcftools'.
- `concat --ligate` on the cohort (0.99x), where most of the time is BGZF compression on one
  thread in both tools.
- Within 5% of bcftools, though not slower: on giab, `view` BCF to VCF,
  `query -f '[%GT %DP %AD]'`, and `norm -d exact`; `call` on small inputs; and `+split-vep -d`.
- Memory: `csq` takes about twice bcftools' (287 MB against 139 MB; 748 MB against 461 MB with
  `-Ot`). `norm -m-any`, `annotate -a` with a tab file, `roh -s` with one sample,
  `query -f '%LINE'` on giab, and `mpileup -f` take 7 to 18 MB more. With `--threads`, memory
  grows with the thread count: at `--threads 16`, `mpileup | call` takes 192 MB (bcftools:
  50 MB), `view` BCF to VCF.gz 175 MB (18 MB), and `concat` 118 MB (19 MB).
