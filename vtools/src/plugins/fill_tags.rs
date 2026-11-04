//! `+fill-tags`: set INFO/AN, AC, AC_Hom, AC_Het, AC_Hemi, AF, MAF, NS, HWE, ExcHet, END,
//! TYPE, F_MISSING, FORMAT/VAF and VAF1, and `TAG=func(...)` expressions, for all samples and
//! for populations of samples (bcftools' `plugins/fill-tags.c`).

use std::path::Path;

use anyhow::bail;
use clap::{CommandFactory, FromArgMatches, Parser};
use pasta::variant::{
    Header, HeaderRecordType, Number, Record, Unpack, ValueType,
    encoding::{FLOAT_MISSING, FLOAT_VECTOR_END, INT32_MISSING, INT32_VECTOR_END},
    index::split_idx_name,
    readlist,
};

use super::{Exit, InitHeaders, Loaded, Output, Plugin};
use crate::{filter::Filter, raw_input::log_open_error, shortopt};

/// The description `plugin -lv` prints.
pub const ABOUT: &str = "Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS; FORMAT/VAF and more.\n";

/// The plugin's usage text.
pub fn usage() -> String {
    "\n\
About: Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS\n\
\x20      FORMAT tag VAF, custom INFO/TAG=func(FMT/TAG).\n\
\x20      See examples below, run with -l for detailed description.\n\
Usage: vtools +fill-tags [General Options] -- [Plugin Options]\n\
Options:\n\
\x20  run \"vtools plugin\" for a list of common options\n\
\n\
Plugin options:\n\
\x20  -d, --drop-missing          do not count half-missing genotypes \"./1\" as hemizygous\n\
\x20  -l, --list-tags             list available tags with description\n\
\x20  -t, --tags LIST             list of output tags, \"all\" for all tags\n\
\x20  -S, --samples-file FILE     list of samples (first column) and comma-separated list of populations (second column)\n\
\n\
Example:\n\
\x20  # Print a detailed list of available tags\n\
\x20  vtools +fill-tags -- -l\n\
\n\
\x20  # Fill INFO/AN and INFO/AC\n\
\x20  vtools +fill-tags in.bcf -Ob -o out.bcf -- -t AN,AC\n\
\n\
\x20  # Fill (almost) all available tags\n\
\x20  vtools +fill-tags in.bcf -Ob -o out.bcf -- -t all\n\
\n\
\x20  # Calculate HWE for sample groups (possibly multiple) read from a file\n\
\x20  vtools +fill-tags in.bcf -Ob -o out.bcf -- -S sample-group.txt -t HWE\n\
\n\
\x20  # Calculate total read depth (INFO/DP) from per-sample depths (FORMAT/DP)\n\
\x20  vtools +fill-tags in.bcf -Ob -o out.bcf -- -t 'DP:1=int(sum(FORMAT/DP))'\n\
\n\
\x20  # Calculate per-sample read depth (FORMAT/DP) from per-sample allelic depths (FORMAT/AD)\n\
\x20  vtools +fill-tags in.bcf -Ob -o out.bcf -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))'\n\
\n\
\x20  # Add number of samples which pass (INFO/good) and fail (INFO/bad) a binomial test\n\
\x20  vtools +fill-tags in.bcf -- -t 'good=N_PASS(binom(FMT/AD[:0],FMT/AD[:1])>=1e-5)','bad=N_PASS(binom(FMT/AD[:0],FMT/AD[:1])<1e-5)'\n\
\n\
\x20  # Annotate with phred-scaled p-value of fisher exact test, use the DP4 or ADF,ADR tags\n\
\x20  vtools +fill-tags in.bcf -- -t 'FMT/FT:1=phred(fisher(FMT/DP4))'\n\
\x20  vtools +fill-tags in.bcf -- -t 'FMT/FT:1=phred(fisher(FMT/ADF,FMT/ADR))'\n\
\n\
\x20  # Annotate with allelic fraction\n\
\x20  vtools +fill-tags in.bcf -Ob -o out.bcf -- -t FORMAT/VAF\n\
\n"
    .to_owned()
}

/// The `-l` text (bcftools' `list_tags`, printed through `error()`).
const LIST_TAGS: &str = "\
INFO/AC        Number:A  Type:Integer  ..  Allele count in genotypes\n\
INFO/AC_Hom    Number:A  Type:Integer  ..  Allele counts in homozygous genotypes\n\
INFO/AC_Het    Number:A  Type:Integer  ..  Allele counts in heterozygous genotypes\n\
INFO/AC_Hemi   Number:A  Type:Integer  ..  Allele counts in hemizygous genotypes\n\
INFO/AF        Number:A  Type:Float    ..  Allele frequency from FMT/GT or AC,AN if FMT/GT is not present\n\
INFO/AN        Number:1  Type:Integer  ..  Total number of alleles in called genotypes\n\
INFO/ExcHet    Number:A  Type:Float    ..  Test excess heterozygosity; 1=good, 0=bad\n\
INFO/END       Number:1  Type:Integer  ..  End position of the variant\n\
INFO/F_MISSING Number:1  Type:Float    ..  Fraction of missing genotypes, synonymous with 'F_MISSING=F_PASS(GT=\"mis\")'\n\
INFO/HWE       Number:A  Type:Float    ..  HWE test (PMID:15789306); 1=good, 0=bad\n\
INFO/MAF       Number:1  Type:Float    ..  Frequency of the second most common allele\n\
INFO/NS        Number:1  Type:Integer  ..  Number of samples with data\n\
INFO/TYPE      Number:.  Type:String   ..  The record type (REF,SNP,MNP,INDEL,etc)\n\
FORMAT/VAF     Number:A  Type:Float    ..  The fraction of reads with the alternate allele, requires FORMAT/AD or ADF+ADR\n\
FORMAT/VAF1    Number:1  Type:Float    ..  The same as FORMAT/VAF but for all alternate alleles cumulatively\n\
TAG:Number=Type(EXPR)                  ..  Experimental support for user expressions such as DP:1=int(sum(DP))\n\
\x20              If Number and Type are not given (e.g. DP=sum(DP)), variable number (Number=.) of floating point\n\
\x20              values (Type=Float) will be used.\n";

/// The plugin options (bcftools' getopt string `?ht:dS:l`).
#[derive(Debug, Parser)]
#[command(
    no_binary_name = true,
    disable_help_flag = true,
    disable_version_flag = true,
    args_override_self = true
)]
struct Args {
    /// Arguments left over (an error).
    rest: Vec<String>,
    /// Do not count half-missing genotypes "./1" as hemizygous.
    #[arg(short = 'd', long = "drop-missing")]
    drop_missing: bool,
    /// List available tags with description.
    #[arg(short = 'l', long = "list-tags", action = clap::ArgAction::Count)]
    list_tags: u8,
    /// List of output tags, "all" for all tags.
    #[arg(short = 't', long = "tags", allow_hyphen_values = true)]
    tags: Option<String>,
    /// List of samples (first column) and comma-separated list of populations (second column).
    #[arg(short = 'S', long = "samples-file")]
    samples_file: Option<String>,
    /// Print the usage.
    #[arg(short = 'h', long = "help", short_alias = '?', action = clap::ArgAction::Count)]
    help: u8,
}

/// The options that end bcftools' getopt loop when it reads them.
#[derive(Clone, Copy, Debug)]
enum Stop {
    /// `-l`: the tag list.
    List,
    /// `-h`: the usage.
    Help,
}

/// The tags to set (bcftools' `SET_*` bits).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Tags(u32);

impl Tags {
    const AN: u32 = 1 << 0;
    const AC: u32 = 1 << 1;
    const AC_HOM: u32 = 1 << 2;
    const AC_HET: u32 = 1 << 3;
    const AC_HEMI: u32 = 1 << 4;
    const AF: u32 = 1 << 5;
    const NS: u32 = 1 << 6;
    const MAF: u32 = 1 << 7;
    const HWE: u32 = 1 << 8;
    const EXC_HET: u32 = 1 << 9;
    const FUNC: u32 = 1 << 10;
    const END: u32 = 1 << 11;
    const TYPE: u32 = 1 << 12;
    const VAF: u32 = 1 << 13;
    const VAF1: u32 = 1 << 14;

    fn has(self, bits: u32) -> bool {
        self.0 & bits != 0
    }
}

/// Genotype counts of one allele in one population (bcftools' `counts_t`).
#[derive(Clone, Copy, Debug, Default)]
struct Counts {
    nhom: i32,
    nhet: i32,
    nhemi: i32,
    nac: i32,
    /// The part of `nhom` from polyploid genotypes, which HWE and ExcHet leave out.
    poly_nhom: i32,
    /// The part of `nhet` from polyploid genotypes, which HWE and ExcHet leave out.
    poly_nhet: i32,
}

impl Counts {
    fn total(self) -> i32 {
        self.nhet + self.nhom + self.nhemi + self.nac
    }

    /// Allele count in diploid homozygous genotypes, the input of HWE and ExcHet.
    fn diploid_nhom(self) -> i32 {
        self.nhom - self.poly_nhom
    }

    /// Allele count in diploid heterozygous genotypes, the input of HWE and ExcHet.
    fn diploid_nhet(self) -> i32 {
        self.nhet - self.poly_nhet
    }
}

/// The value type of an expression tag.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FtfType {
    Int,
    Float,
}

/// A `TAG=func(...)` expression for one population (bcftools' `ftf_t`).
struct Ftf {
    /// The tag with the population suffix.
    tag: Vec<u8>,
    ty: FtfType,
    /// `Some(n)` for `TAG:n=...`, `None` for a variable number of values.
    count: Option<i64>,
    info: bool,
    filter: Filter,
    fval: Vec<f32>,
    ival: Vec<i32>,
}

/// A population: a list of samples and its counts (bcftools' `pop_t`).
struct Pop {
    name: String,
    suffix: String,
    /// Header sample indices.
    samples: Vec<usize>,
    ns: i32,
    counts: Vec<Counts>,
    ftf: Vec<Ftf>,
    /// Output tag IDs, by tag.
    ids: TagIds,
}

/// The output header IDs of a population's INFO tags.
#[derive(Default)]
struct TagIds {
    an: Option<u32>,
    ac: Option<u32>,
    ns: Option<u32>,
    ac_hom: Option<u32>,
    ac_het: Option<u32>,
    ac_hemi: Option<u32>,
    af: Option<u32>,
    maf: Option<u32>,
    hwe: Option<u32>,
    exc_het: Option<u32>,
}

/// The allele fraction tag to compute from FORMAT/AD.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum VafTag {
    /// FORMAT/VAF: the fraction of each alternate allele.
    Vaf,
    /// FORMAT/VAF1: the fraction of all alternate alleles together.
    Vaf1,
}

/// The plugin state (bcftools' `args_t`).
struct FillTags {
    tags: Tags,
    warned: u32,
    drop_missing: bool,
    /// Whether FORMAT is needed (bcftools' `unpack & BCF_UN_FMT`).
    fmt: bool,
    gt_id: Option<u32>,
    /// The populations; the last is all samples.
    pops: Vec<Pop>,
    /// For each sample, the populations it is in, the all-samples one first.
    smpl2pop: Vec<Vec<usize>>,
    farr: Vec<f32>,
    iarr: Vec<i32>,
    hwe_probs: Vec<f64>,
    /// Distinct alleles of the current genotype (bcftools' `kbitset`).
    alleles: Vec<usize>,
    vaf_id: Option<u32>,
    vaf1_id: Option<u32>,
}

/// Initializes the plugin (bcftools' `init`).
pub fn init(argv: &[String], hdrs: InitHeaders<'_>) -> anyhow::Result<Loaded> {
    // argv[0] is the input file, which getopt takes as the program name (for its messages).
    let prog = argv.first().map_or("fill-tags", String::as_str);
    let argv = shortopt::rewrite_strings(&Args::command(), prog, argv.get(1..).unwrap_or_default())
        .map_err(|msg| Exit::error(format!("{msg}{}", usage())))?;
    let matches = Args::command()
        .try_get_matches_from(argv)
        .map_err(|_| Exit::error(usage()))?;
    let args = Args::from_arg_matches(&matches).map_err(|_| Exit::error(usage()))?;
    let stops = super::in_order(
        &matches,
        &[],
        &[
            ("list_tags", args.list_tags, Stop::List),
            ("help", args.help, Stop::Help),
        ],
    );
    match stops.first() {
        Some(Stop::List) => return Err(Exit::error(LIST_TAGS)),
        Some(Stop::Help) => return Err(Exit::error(usage())),
        None => {}
    }
    if !args.rest.is_empty() {
        return Err(Exit::error(usage()));
    }
    let Some((hdr, out)) = hdrs else {
        return Err(Exit::error(usage()));
    };
    let mut ft = FillTags {
        tags: Tags::default(),
        warned: 0,
        drop_missing: args.drop_missing,
        fmt: false,
        gt_id: hdr.id(b"GT"),
        pops: Vec::new(),
        smpl2pop: Vec::new(),
        farr: Vec::new(),
        iarr: Vec::new(),
        hwe_probs: Vec::new(),
        alleles: Vec::new(),
        vaf_id: None,
        vaf1_id: None,
    };
    if let Some(f) = &args.samples_file {
        ft.parse_samples(hdr, f)?;
    }
    ft.init_pops(hdr);
    let tags_str = args.tags.as_deref().unwrap_or("all");
    let tags = ft.parse_tags(hdr, out, tags_str)?;
    ft.tags.0 |= tags;
    ft.add_header_lines(out);
    ft.resolve_ids(out);
    Ok(Loaded {
        plugin: Box::new(ft),
        output: Output::Records,
        max_unpack: Unpack::empty(),
    })
}

/// C's `isspace`.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// Splits at `delim`, skipping empty fields (HTSlib's `ksplit_core`).
fn ksplit(s: &[u8], delim: u8) -> impl Iterator<Item = &[u8]> {
    s.split(move |&c| c == delim).filter(|f| !f.is_empty())
}

/// Splits a tag list at commas outside brackets (bcftools' `parse_tag_list`).
fn parse_tag_list(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut start = 0;
    for (i, c) in s.bytes().enumerate() {
        if c == b',' && depth == 0 {
            out.push(&s[start..i]);
            start = i + 1;
        } else if c == b'[' || c == b'(' {
            depth += 1;
        } else if (c == b']' || c == b')') && depth > 0 {
            depth -= 1;
        }
    }
    out.push(&s[start..]);
    out
}

/// Whether `s` starts with `prefix`, ignoring case (`strncasecmp`).
fn starts_with_nocase(s: &str, prefix: &str) -> bool {
    s.len() >= prefix.len() && s.as_bytes()[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
}

/// bcftools' `float_set_from_double`.
fn float_from_double(v: f64) -> f32 {
    if is_missing_or_end(v) {
        f32::from_bits(FLOAT_MISSING)
    } else {
        v as f32
    }
}

/// bcftools' `int32_from_double`: `(int)round(v)`, which on x86 (`cvttsd2si`) gives
/// `INT_MIN` (written as missing) for NaN and values out of range, where Rust's cast
/// saturates.
fn int32_from_double(v: f64) -> i32 {
    if is_missing_or_end(v) {
        return INT32_MISSING;
    }
    let r = v.round();
    if r.is_nan() || !(-2_147_483_648.0..2_147_483_648.0).contains(&r) {
        i32::MIN
    } else {
        r as i32
    }
}

/// HTSlib's `bcf_double_is_missing_or_vector_end`.
fn is_missing_or_end(v: f64) -> bool {
    let b = v.to_bits();
    b == 0x7ff0_0000_0000_0001 || b == 0x7ff0_0000_0000_0002
}

impl FillTags {
    /// Reads the `-S` populations (bcftools' `parse_samples`); `FILE##idx##INDEX` reads
    /// `FILE`, as `hts_open` does.
    fn parse_samples(&mut self, hdr: &Header, fname: &str) -> anyhow::Result<()> {
        let data = split_idx_name(fname).0;
        let lines = readlist::read_all_lines(Path::new(data)).map_err(|e| {
            log_open_error(data, &e);
            anyhow::anyhow!("Could not read: {fname}")
        })?;
        let mut seen = std::collections::HashSet::new();
        let mut nsmpl = 0;
        for line in &lines {
            let bad = || {
                anyhow::anyhow!(
                    "Could not parse the file: {}",
                    String::from_utf8_lossy(line)
                )
            };
            // The populations are the last whitespace-separated word, the sample the rest.
            let mut end = line.len();
            while end > 0 && is_space(line[end - 1]) {
                end -= 1;
            }
            if end <= 1 {
                return Err(bad());
            }
            let mut beg = end;
            while beg > 0 && !is_space(line[beg - 1]) {
                beg -= 1;
            }
            if beg <= 1 {
                return Err(bad());
            }
            let pop_names = &line[beg..end];
            let mut send = beg - 1;
            while send > 0 && is_space(line[send - 1]) {
                send -= 1;
            }
            if send == 0 {
                return Err(bad());
            }
            let smpl = &line[..send];
            let smpl_str = String::from_utf8_lossy(smpl);
            let Some(ismpl) = hdr.sample_id(smpl) else {
                eprintln!("Warning: The sample not present in the VCF: {smpl_str}");
                continue;
            };
            if !seen.insert(smpl.to_vec()) {
                eprintln!("Warning: The sample is listed twice in {fname}: {smpl_str}");
                continue;
            }
            for pop_name in ksplit(pop_names, b',') {
                let pop_name = String::from_utf8_lossy(pop_name).into_owned();
                let ipop = match self.pops.iter().position(|p| p.name == pop_name) {
                    Some(i) => i,
                    None => {
                        self.pops.push(Pop::new(format!("_{pop_name}"), pop_name));
                        self.pops.len() - 1
                    }
                };
                self.pops[ipop].samples.push(ismpl as usize);
            }
            nsmpl += 1;
        }
        if nsmpl != hdr.sample_count() {
            eprintln!(
                "Warning: {nsmpl} samples in the list, {} samples in the VCF.",
                hdr.sample_count()
            );
        }
        if self.pops.is_empty() {
            bail!("No populations given?");
        }
        Ok(())
    }

    /// Adds the population of all samples and maps samples to populations (bcftools'
    /// `init_pops`).
    fn init_pops(&mut self, hdr: &Header) {
        self.pops.push(Pop::new(String::new(), String::new()));
        let all = self.pops.len() - 1;
        self.smpl2pop = vec![vec![all]; hdr.sample_count()];
        for (i, pop) in self.pops.iter().enumerate().take(all) {
            for &s in &pop.samples {
                self.smpl2pop[s].push(i);
            }
        }
    }

    /// Parses `-t` (bcftools' `parse_tags`); returns the tags to set.
    fn parse_tags(&mut self, hdr: &Header, out: &mut Header, s: &str) -> anyhow::Result<u32> {
        self.warned = 0;
        let mut flag = 0;
        for tag in parse_tag_list(s) {
            let is = |names: &[&str]| names.iter().any(|n| tag.eq_ignore_ascii_case(n));
            if is(&["all"]) {
                flag |= !(Tags::END | Tags::TYPE);
                flag |= self.parse_func(hdr, out, "F_MISSING:1=F_MISSING", "F_MISSING")?;
                self.warned = !(Tags::END | Tags::TYPE);
                self.fmt = true;
                continue;
            }
            let fixed = [
                (&["AN", "INFO/AN"][..], Tags::AN),
                (&["AC", "INFO/AC"], Tags::AC),
                (&["NS", "INFO/NS"], Tags::NS),
                (&["AC_Hom", "INFO/AC_Hom"], Tags::AC_HOM),
                (&["AC_Het", "INFO/AC_Het"], Tags::AC_HET),
                // bcftools spells the second name INFO_Hemi.
                (&["AC_Hemi", "INFO/AC_Hemi", "INFO_Hemi"], Tags::AC_HEMI),
                (&["AF", "INFO/AF"], Tags::AF),
                (&["MAF", "INFO/MAF"], Tags::MAF),
                (&["HWE", "INFO/HWE"], Tags::HWE),
                (&["ExcHet", "INFO/ExcHet"], Tags::EXC_HET),
                (&["VAF", "FORMAT/VAF"], Tags::VAF),
                (&["VAF1", "FORMAT/VAF1"], Tags::VAF1),
            ];
            if let Some((_, bit)) = fixed.iter().find(|(names, _)| is(names)) {
                flag |= bit;
                self.fmt = true;
            } else if is(&["END", "INFO/END"]) {
                flag |= Tags::END;
            } else if is(&["TYPE", "INFO/TYPE"]) {
                flag |= Tags::TYPE;
            } else if is(&["F_MISSING", "INFO/F_MISSING"]) {
                flag |= self.parse_func(hdr, out, "F_MISSING:1=F_MISSING", "F_MISSING")?;
                self.fmt = true;
            } else if let Some(eq) = tag.find('=') {
                flag |= self.parse_func(hdr, out, tag, &tag[eq + 1..])?;
                self.fmt = true;
            } else {
                return Err(Exit::usage(format!(
                    "Error parsing \"--tags {s}\": the tag \"{tag}\" is not supported\n"
                )));
            }
        }
        Ok(flag)
    }

    /// Adds an expression tag for every population (bcftools' `parse_func`).
    fn parse_func(
        &mut self,
        hdr: &Header,
        out: &mut Header,
        tag_expr: &str,
        expr: &str,
    ) -> anyhow::Result<u32> {
        for i in 0..self.pops.len() {
            let ftf = parse_func_pop(hdr, out, &self.pops[i], tag_expr, expr)?;
            self.pops[i].ftf.push(ftf);
        }
        Ok(Tags::FUNC)
    }

    /// Adds the header lines of the fixed tags, for each population.
    fn add_header_lines(&self, out: &mut Header) {
        let per_pop = |out: &mut Header, f: &dyn Fn(&str, &str) -> String| {
            for p in &self.pops {
                let name = if p.name.is_empty() {
                    String::new()
                } else {
                    format!(" in {}", p.name)
                };
                let _ = out.append(f(&p.suffix, &name).as_bytes());
            }
        };
        let t = self.tags;
        if t.has(Tags::AN) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=AN{s},Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes{n}\">"
                )
            });
        }
        if t.has(Tags::AC) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=AC{s},Number=A,Type=Integer,Description=\"Allele count in genotypes{n}\">"
                )
            });
        }
        if t.has(Tags::NS) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=NS{s},Number=1,Type=Integer,Description=\"Number of samples with data{n}\">"
                )
            });
        }
        if t.has(Tags::AC_HOM) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=AC_Hom{s},Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes{n}\">"
                )
            });
        }
        if t.has(Tags::AC_HET) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=AC_Het{s},Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes{n}\">"
                )
            });
        }
        if t.has(Tags::AC_HEMI) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=AC_Hemi{s},Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes{n}\">"
                )
            });
        }
        if t.has(Tags::AF) {
            per_pop(out, &|s, n| {
                format!("##INFO=<ID=AF{s},Number=A,Type=Float,Description=\"Allele frequency{n}\">")
            });
        }
        if t.has(Tags::MAF) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=MAF{s},Number=1,Type=Float,Description=\"Frequency of the second most common allele{n}\">"
                )
            });
        }
        if t.has(Tags::HWE) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=HWE{s},Number=A,Type=Float,Description=\"HWE test{n} (PMID:15789306); 1=good, 0=bad\">"
                )
            });
        }
        if t.has(Tags::END) {
            let _ = out.append(b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">");
        }
        if t.has(Tags::TYPE) {
            let _ =
                out.append(b"##INFO=<ID=TYPE,Number=.,Type=String,Description=\"Variant type\">");
        }
        if t.has(Tags::EXC_HET) {
            per_pop(out, &|s, n| {
                format!(
                    "##INFO=<ID=ExcHet{s},Number=A,Type=Float,Description=\"Test excess heterozygosity{n}; 1=good, 0=bad\">"
                )
            });
        }
        if t.has(Tags::VAF) {
            let _ = out.append(b"##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"The fraction of reads with alternate allele (nALT/nSumAll)\">");
        }
        if t.has(Tags::VAF1) {
            let _ = out.append(b"##FORMAT=<ID=VAF1,Number=1,Type=Float,Description=\"The fraction of reads with alternate alleles (nSumALT/nSumAll)\">");
        }
    }

    /// Looks up the output IDs of the tags.
    fn resolve_ids(&mut self, out: &Header) {
        let id = |name: &str, suffix: &str| out.id(format!("{name}{suffix}").as_bytes());
        for p in &mut self.pops {
            let s = &p.suffix;
            let ids = TagIds {
                an: id("AN", s),
                ac: id("AC", s),
                ns: id("NS", s),
                ac_hom: id("AC_Hom", s),
                ac_het: id("AC_Het", s),
                ac_hemi: id("AC_Hemi", s),
                af: id("AF", s),
                maf: id("MAF", s),
                hwe: id("HWE", s),
                exc_het: id("ExcHet", s),
            };
            p.ids = ids;
        }
        self.vaf_id = out.id(b"VAF");
        self.vaf1_id = out.id(b"VAF1");
    }
}

impl Pop {
    fn new(suffix: String, name: String) -> Self {
        Self {
            name,
            suffix,
            samples: Vec::new(),
            ns: 0,
            counts: Vec::new(),
            ftf: Vec::new(),
            ids: TagIds::default(),
        }
    }
}

/// Parses `TAG[:N]=EXPR` for one population and adds its header line (bcftools'
/// `parse_func_pop`).
fn parse_func_pop(
    hdr: &Header,
    out: &mut Header,
    pop: &Pop,
    tag_expr: &str,
    expr: &str,
) -> anyhow::Result<Ftf> {
    let mut ty = FtfType::Float;
    let mut info = true;
    let eq = tag_expr.find('=').unwrap_or(tag_expr.len());
    let mut beg = 0;
    if starts_with_nocase(tag_expr, "info/") {
        beg = 5;
    } else if starts_with_nocase(tag_expr, "format/") {
        beg = 7;
        info = false;
    } else if starts_with_nocase(tag_expr, "fmt/") {
        beg = 4;
        info = false;
    }
    let mut dst = &tag_expr[beg.min(eq)..eq];
    let mut count = None;
    if let Some(colon) = dst.find(':') {
        let digits = &dst[colon + 1..];
        let (n, end) = pasta::variant::number::strtol(digits.as_bytes(), 0);
        if end != digits.len() {
            bail!("Could not parse the expression: {tag_expr}");
        }
        count = Some(n);
        dst = &dst[..colon];
    }

    let mut filter_expr = expr;
    if expr.ends_with(')') {
        for (prefix, t) in [
            ("int(", FtfType::Int),
            ("integer(", FtfType::Int),
            ("float(", FtfType::Float),
        ] {
            if starts_with_nocase(expr, prefix) {
                filter_expr = &expr[prefix.len()..expr.len() - 1];
                ty = t;
                break;
            }
        }
    }

    let tag = format!("{dst}{}", pop.suffix);
    let kind = if info {
        HeaderRecordType::Info
    } else {
        HeaderRecordType::Format
    };
    let col = if info { "INFO" } else { "FORMAT" };
    let existing = hdr.id(tag.as_bytes()).filter(|&id| hdr.id_exists(kind, id));
    if let Some(id) = existing {
        let number = hdr.id_number(kind, id);
        let want_code = if count.is_some() {
            Number::Fixed(0).code()
        } else {
            Number::Variable.code()
        };
        if number.map(Number::code) != Some(want_code) {
            bail!(
                "Error: the field {col}/{tag} already exists and its Number definition is different"
            );
        }
        if let (Some(n), Some(Number::Fixed(have))) = (count, number)
            && i64::from(have) != n
        {
            bail!(
                "Error: the field {col}/{tag} already exists and its Number definition is different from {n}"
            );
        }
        let want_type = match ty {
            FtfType::Int => ValueType::Integer,
            FtfType::Float => ValueType::Float,
        };
        if hdr.id_type(kind, id) != Some(want_type) {
            bail!(
                "Error: the field {col}/{tag} already exists and its Type definition is different"
            );
        }
    } else {
        let mut line = format!("##{col}=<ID={tag},Number=");
        match count {
            Some(n) => line.push_str(&n.to_string()),
            None => line.push('.'),
        }
        line.push_str(match ty {
            FtfType::Int => ",Type=Integer",
            FtfType::Float => ",Type=Float",
        });
        line.push_str(",Description=\"Added by +fill-tags expression ");
        for c in tag_expr.chars() {
            if c == '"' {
                line.push('\\');
            }
            line.push(c);
        }
        if !pop.name.is_empty() {
            line.push_str(" in ");
            line.push_str(&pop.name);
        }
        line.push_str("\">");
        let _ = out.append(line.as_bytes());
    }
    let mut filter = Filter::init(hdr, filter_expr)?;
    if !pop.name.is_empty() {
        let mut samples = vec![false; hdr.sample_count()];
        for &s in &pop.samples {
            samples[s] = true;
        }
        filter.set_samples(&samples);
    }
    Ok(Ftf {
        tag: tag.into_bytes(),
        ty,
        count,
        info,
        filter,
        fval: Vec::new(),
        ival: Vec::new(),
    })
}

/// The error of a failed update.
fn update_error(tag: &[u8], hdr: &Header, rec: &Record) -> anyhow::Error {
    anyhow::anyhow!(
        "Error occurred while updating {} at {}:{}",
        String::from_utf8_lossy(tag),
        String::from_utf8_lossy(rec.chrom_safe(hdr)),
        rec.pos() + 1
    )
}

/// Evaluates an expression tag and sets its values (bcftools' `ftf_filter_expr`).
fn ftf_filter_expr(
    ftf: &mut Ftf,
    hdr: &Header,
    out: &Header,
    rec: &mut Record,
) -> anyhow::Result<()> {
    ftf.filter.test(hdr, rec)?;
    let (vals, nval1): (&[f64], usize) = ftf.filter.doubles().unwrap_or((&[], 0));
    let nval = vals.len();
    let ret = if ftf.info {
        let nfill = ftf.count.map_or(nval, |c| usize::try_from(c).unwrap_or(0));
        let ncopy = nval.min(nfill);
        match ftf.ty {
            FtfType::Float => {
                ftf.fval.clear();
                ftf.fval
                    .extend(vals[..ncopy].iter().map(|&v| float_from_double(v)));
                ftf.fval.resize(nfill, f32::from_bits(FLOAT_MISSING));
                rec.update_info_f32(out, &ftf.tag, &ftf.fval)
            }
            FtfType::Int => {
                ftf.ival.clear();
                ftf.ival
                    .extend(vals[..ncopy].iter().map(|&v| int32_from_double(v)));
                ftf.ival.resize(nfill, INT32_MISSING);
                rec.update_info_i32(out, &ftf.tag, &ftf.ival)
            }
        }
    } else {
        let nsmpl = rec.sample_count();
        if nval1 * nsmpl > nval {
            bail!(
                "Error: The expression returned {nval} value(s) total but {nval1} value(s) per sample are required; there are {nsmpl} samples"
            );
        }
        let nfill = ftf.count.map_or(nval1, |c| usize::try_from(c).unwrap_or(0));
        let ncopy = nval1.min(nfill);
        match ftf.ty {
            FtfType::Float => {
                ftf.fval.clear();
                for k in 0..nsmpl {
                    let src = &vals[k * nval1..k * nval1 + ncopy];
                    ftf.fval.extend(src.iter().map(|&v| float_from_double(v)));
                    for _ in ncopy..nfill {
                        ftf.fval.push(f32::from_bits(FLOAT_MISSING));
                    }
                }
                rec.update_format_f32(out, &ftf.tag, &ftf.fval)
            }
            FtfType::Int => {
                ftf.ival.clear();
                for k in 0..nsmpl {
                    let src = &vals[k * nval1..k * nval1 + ncopy];
                    ftf.ival.extend(src.iter().map(|&v| int32_from_double(v)));
                    for _ in ncopy..nfill {
                        ftf.ival.push(INT32_MISSING);
                    }
                }
                rec.update_format_i32(out, &ftf.tag, &ftf.ival)
            }
        }
    };
    ret.map_err(|_| update_error(&ftf.tag, hdr, rec))
}

/// Hardy-Weinberg and excess heterozygosity p-values (Wigginton 2005, PMID: 15789306;
/// bcftools' `calc_hwe`).
fn calc_hwe(probs: &mut Vec<f64>, nref: i32, nalt: i32, nhet: i32) -> anyhow::Result<(f32, f32)> {
    let ngt = (nref + nalt) / 2;
    let nrare = nref.min(nalt);
    if (nrare & 1) ^ (nhet & 1) != 0 {
        bail!(
            "nrare/nhet should be both odd or even: nrare={nrare} nref={nref} nalt={nalt} nhet={nhet}"
        );
    }
    if nrare < nhet {
        bail!("Fewer rare alleles than hets? nrare={nrare} nref={nref} nalt={nalt} nhet={nhet}");
    }
    if (nref + nalt) & 1 != 0 {
        bail!("Expected diploid genotypes: nref={nref} nalt={nalt}");
    }
    let nrare_u = nrare as usize;
    probs.clear();
    probs.resize(nrare_u + 1, 0.0);

    let mut mid =
        (f64::from(nrare) * f64::from(nref + nalt - nrare) / f64::from(nref + nalt)) as i32;
    if (nrare & 1) ^ (mid & 1) != 0 {
        mid += 1;
    }
    let mut hom_r = (nrare - mid) / 2;
    let mut hom_c = ngt - mid - hom_r;
    probs[mid as usize] = 1.0;
    let mut sum = 1.0;
    let mut het = mid;
    while het > 1 {
        let h = het as usize;
        probs[h - 2] = probs[h] * f64::from(het) * (f64::from(het) - 1.0)
            / (4.0 * (f64::from(hom_r) + 1.0) * (f64::from(hom_c) + 1.0));
        sum += probs[h - 2];
        hom_r += 1;
        hom_c += 1;
        het -= 2;
    }
    het = mid;
    hom_r = (nrare - mid) / 2;
    hom_c = ngt - het - hom_r;
    while het <= nrare - 2 {
        let h = het as usize;
        probs[h + 2] = probs[h] * 4.0 * f64::from(hom_r) * f64::from(hom_c)
            / ((f64::from(het) + 2.0) * (f64::from(het) + 1.0));
        sum += probs[h + 2];
        hom_r -= 1;
        hom_c -= 1;
        het += 2;
    }
    for p in probs.iter_mut() {
        *p /= sum;
    }
    let nhet_u = nhet as usize;
    let mut prob = probs[nhet_u];
    for p in &probs[nhet_u + 1..=nrare_u] {
        prob += p;
    }
    let exc_het = prob;
    prob = 0.0;
    for &p in probs.iter() {
        if p > probs[nhet_u] {
            continue;
        }
        prob += p;
    }
    if prob > 1.0 {
        prob = 1.0;
    }
    Ok((prob as f32, exc_het as f32))
}

/// Reads value `i` of a typed integer vector as a 32-bit value, mapping the narrower
/// missing and vector-end values.
#[inline]
fn gt_value(data: &[u8], width: usize, i: usize) -> i32 {
    match width {
        1 => {
            let v = data[i] as i8;
            match v {
                i8::MIN => INT32_MISSING,
                -127 => INT32_VECTOR_END,
                _ => i32::from(v),
            }
        }
        2 => {
            let v = i16::from_le_bytes([data[2 * i], data[2 * i + 1]]);
            match v {
                i16::MIN => INT32_MISSING,
                -32767 => INT32_VECTOR_END,
                _ => i32::from(v),
            }
        }
        _ => i32::from_le_bytes([
            data[4 * i],
            data[4 * i + 1],
            data[4 * i + 2],
            data[4 * i + 3],
        ]),
    }
}

/// How [`count_polyploid`] counts the alleles of a genotype.
#[derive(Clone, Copy)]
enum GtKind {
    /// Partly missing, with `-d`.
    Half,
    Het,
    Hom,
    /// Partly missing, without `-d`.
    Hemi,
}

/// Adds each called allele of a polyploid genotype (`width`-byte values, all checked against
/// the record's alleles) to its counts, keeping hom and het counts apart for HWE and ExcHet.
#[cold]
#[inline(never)]
fn count_polyploid(counts: &mut [Counts], gt: &[u8], width: usize, kind: GtKind) {
    for k in 0..gt.len() / width {
        let v = gt_value(gt, width, k);
        if v >> 1 == 0 || v == INT32_MISSING {
            continue;
        }
        let c = &mut counts[((v >> 1) - 1) as usize];
        match kind {
            GtKind::Half => c.nac += 1,
            GtKind::Het => {
                c.nhet += 1;
                c.poly_nhet += 1;
            }
            GtKind::Hom => {
                c.nhom += 1;
                c.poly_nhom += 1;
            }
            GtKind::Hemi => c.nhemi += 1,
        }
    }
}

impl FillTags {
    /// Counts genotypes and sets the INFO tags derived from them (bcftools' `process_fmt`).
    fn process_fmt(&mut self, hdr: &Header, out: &Header, rec: &mut Record) -> anyhow::Result<()> {
        let Some(gt_id) = self.gt_id else {
            return Ok(());
        };
        rec.unpack(Unpack::FMT);
        let Some(fi) = rec
            .format_fields()
            .iter()
            .position(|f| f.id() == gt_id as i32)
        else {
            return Ok(());
        };
        let nals = rec.allele_count();
        if nals == 0 {
            // A BCF record without even a REF allele has nothing to count.
            return Ok(());
        }
        self.iarr.resize(nals.max(self.iarr.len()), 0);
        self.farr.resize((2 * nals).max(self.farr.len()), 0.0);
        for pop in &mut self.pops {
            pop.ns = 0;
            pop.counts.clear();
            pop.counts.resize(nals, Counts::default());
        }
        match rec.format_fields()[fi].type_code() {
            1 => self.count_gts::<1>(hdr, rec, fi, nals)?,
            2 => self.count_gts::<2>(hdr, rec, fi, nals)?,
            3 => self.count_gts::<4>(hdr, rec, fi, nals)?,
            t => bail!(
                "The GT type is not recognised: {t} at {}:{}",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1
            ),
        }
        self.set_count_tags(hdr, out, rec, nals)
    }

    /// Counts the genotypes of every sample in the record's FORMAT field `fi`, GT values of
    /// `W` bytes (bcftools' `BRANCH_INT`, one instance per integer type).
    fn count_gts<const W: usize>(
        &mut self,
        hdr: &Header,
        rec: &Record,
        fi: usize,
        nals: usize,
    ) -> anyhow::Result<()> {
        let field = &rec.format_fields()[fi];
        let n = field.n();
        let nsmpl = hdr.sample_count();
        let data = rec.format_data(field);
        let size = field.size();
        for i in 0..nsmpl {
            let p = &data[i * size..(i + 1) * size];
            self.alleles.clear();
            let mut nals_gt = 0;
            let mut ial = 0;
            while ial < n {
                let v = gt_value(p, W, ial);
                if v == INT32_VECTOR_END {
                    break;
                }
                ial += 1;
                // HTSlib's bcf_gt_is_missing; a missing sample also reads as missing.
                if v >> 1 == 0 || v == INT32_MISSING {
                    continue;
                }
                let idx = (v >> 1) - 1;
                nals_gt += 1;
                let idx = match usize::try_from(idx) {
                    Ok(a) if a < nals => a,
                    _ => bail!(
                        "Incorrect allele (\"{idx}\") in {} at {}:{}",
                        String::from_utf8_lossy(&hdr.samples()[i]),
                        String::from_utf8_lossy(rec.chrom_safe(hdr)),
                        rec.pos() + 1
                    ),
                };
                if !self.alleles.contains(&idx) {
                    self.alleles.push(idx);
                }
            }
            if nals_gt == 0 {
                continue;
            }
            let is_hom = self.alleles.len() == 1;
            let (is_hemi, is_half) = if nals_gt != ial {
                if self.drop_missing {
                    (false, true)
                } else {
                    (true, false)
                }
            } else {
                (nals_gt == 1, false)
            };
            // bcftools counts each distinct allele once (twice for a hom), whatever the
            // ploidy, so 0/0/1 adds one REF and 1/1/1 two ALTs. vtools counts every allele
            // of a polyploid genotype, as bcftools' own bcf_calc_ac does, and keeps those
            // hom and het counts apart for HWE. Diploid and haploid genotypes take
            // bcftools' path, which counts them right.
            let is_poly = ial > 2;
            for &ipop in &self.smpl2pop[i] {
                let pop = &mut self.pops[ipop];
                if is_poly {
                    let kind = if is_half {
                        GtKind::Half
                    } else if !is_hom {
                        GtKind::Het
                    } else if !is_hemi {
                        GtKind::Hom
                    } else {
                        GtKind::Hemi
                    };
                    count_polyploid(&mut pop.counts, &p[..ial * W], W, kind);
                } else {
                    for &a in &self.alleles {
                        let c = &mut pop.counts[a];
                        if is_half {
                            c.nac += 1;
                        } else if !is_hom {
                            c.nhet += 1;
                        } else if !is_hemi {
                            c.nhom += 2;
                        } else {
                            c.nhemi += 1;
                        }
                    }
                }
                pop.ns += 1;
            }
        }
        Ok(())
    }

    /// Sets NS, AN, AF, MAF, AC, AC_Het, AC_Hom, AC_Hemi, HWE, and ExcHet from the counts.
    fn set_count_tags(
        &mut self,
        hdr: &Header,
        out: &Header,
        rec: &mut Record,
        nals: usize,
    ) -> anyhow::Result<()> {
        let t = self.tags;
        let set_i = |rec: &mut Record, id: Option<u32>, name: &str, sfx: &str, v: &[i32]| {
            let r = match id {
                Some(id) => rec.update_info_i32_id(out, id, v),
                None => Err(pasta::variant::record::UpdateError::NotInHeader),
            };
            r.map_err(|_| update_error(format!("{name}{sfx}").as_bytes(), hdr, rec))
        };
        let set_f = |rec: &mut Record, id: Option<u32>, name: &str, sfx: &str, v: &[f32]| {
            let r = match id {
                Some(id) => rec.update_info_f32_id(out, id, v),
                None => Err(pasta::variant::record::UpdateError::NotInHeader),
            };
            r.map_err(|_| update_error(format!("{name}{sfx}").as_bytes(), hdr, rec))
        };
        if t.has(Tags::NS) {
            for p in &self.pops {
                set_i(rec, p.ids.ns, "NS", &p.suffix, &[p.ns])?;
            }
        }
        if t.has(Tags::AN) {
            for p in &self.pops {
                let an: i32 = p.counts.iter().map(|c| c.total()).sum();
                set_i(rec, p.ids.an, "AN", &p.suffix, &[an])?;
            }
        }
        if t.has(Tags::AF | Tags::MAF) {
            for p in &self.pops {
                let mut an: i32 = 0;
                if nals > 1 {
                    for (j, c) in p.counts.iter().enumerate() {
                        self.farr[j] = c.total() as f32;
                        // bcftools adds the float to an int.
                        an = (an as f32 + self.farr[j]) as i32;
                    }
                    if an != 0 {
                        for v in &mut self.farr[..nals] {
                            *v /= an as f32;
                        }
                    } else {
                        for v in &mut self.farr[..nals] {
                            *v = f32::from_bits(FLOAT_MISSING);
                        }
                    }
                }
                if t.has(Tags::AF) {
                    set_f(rec, p.ids.af, "AF", &p.suffix, &self.farr[1..nals.max(1)])?;
                }
                if nals > 1 && t.has(Tags::MAF) {
                    if an != 0 {
                        self.farr[..nals]
                            .sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                    }
                    set_f(rec, p.ids.maf, "MAF", &p.suffix, &self.farr[1..2])?;
                }
            }
        }
        if t.has(Tags::AC) {
            for p in &self.pops {
                if nals > 1 {
                    for (j, c) in p.counts.iter().enumerate() {
                        self.iarr[j] = c.total();
                    }
                }
                set_i(rec, p.ids.ac, "AC", &p.suffix, &self.iarr[1..nals.max(1)])?;
            }
        }
        let per_alt = |iarr: &mut Vec<i32>, p: &Pop, f: fn(&Counts) -> i32| {
            if nals > 1 {
                for j in 1..nals {
                    iarr[j - 1] = f(&p.counts[j]);
                }
            }
        };
        if t.has(Tags::AC_HET) {
            for p in &self.pops {
                per_alt(&mut self.iarr, p, |c| c.nhet);
                set_i(
                    rec,
                    p.ids.ac_het,
                    "AC_Het",
                    &p.suffix,
                    &self.iarr[..nals - 1],
                )?;
            }
        }
        if t.has(Tags::AC_HOM) {
            for p in &self.pops {
                per_alt(&mut self.iarr, p, |c| c.nhom);
                set_i(
                    rec,
                    p.ids.ac_hom,
                    "AC_Hom",
                    &p.suffix,
                    &self.iarr[..nals - 1],
                )?;
            }
        }
        if t.has(Tags::AC_HEMI) && nals > 1 {
            for p in &self.pops {
                per_alt(&mut self.iarr, p, |c| c.nhemi);
                set_i(
                    rec,
                    p.ids.ac_hemi,
                    "AC_Hemi",
                    &p.suffix,
                    &self.iarr[..nals - 1],
                )?;
            }
        }
        if t.has(Tags::HWE | Tags::EXC_HET) {
            for p in &self.pops {
                let (fhwe, fexc) = self.farr.split_at_mut(nals);
                if nals > 1 {
                    fhwe.fill(0.0);
                    fexc[..nals].fill(0.0);
                    // HWE and ExcHet are diploid tests (calc_hwe fails on an odd allele
                    // count), so they use diploid genotypes only, as bcftools already does
                    // for haploid and partly missing ones; bcftools takes polyploid ones in
                    // with their distinct alleles, as if 0/0/1 were 0/1.
                    let mut nref_tot = p.counts[0].diploid_nhom();
                    for c in &p.counts {
                        nref_tot += c.diploid_nhet();
                    }
                    for j in 1..nals {
                        let c = p.counts[j];
                        let nhet = c.diploid_nhet();
                        let nref = nref_tot - nhet;
                        let nalt = nhet + c.diploid_nhom();
                        (fhwe[j - 1], fexc[j - 1]) = if nref > 0 && nalt > 0 {
                            calc_hwe(&mut self.hwe_probs, nref, nalt, nhet)?
                        } else {
                            (1.0, 1.0)
                        };
                    }
                }
                if t.has(Tags::HWE) {
                    set_f(rec, p.ids.hwe, "HWE", &p.suffix, &fhwe[..nals - 1])?;
                }
                if t.has(Tags::EXC_HET) {
                    set_f(rec, p.ids.exc_het, "ExcHet", &p.suffix, &fexc[..nals - 1])?;
                }
            }
        }
        Ok(())
    }

    /// Sets AF from INFO/AC and AN in a file without samples (bcftools' `process_info_af`).
    fn process_info_af(
        &mut self,
        hdr: &Header,
        out: &Header,
        rec: &mut Record,
    ) -> anyhow::Result<()> {
        if !self.tags.has(Tags::AF) || hdr.sample_count() > 0 {
            return Ok(());
        }
        if !matches!(rec.get_info_i32(hdr, b"AN", &mut self.iarr), Ok(1)) {
            return Ok(());
        }
        let an = self.iarr[0];
        if an == 0 {
            return Ok(());
        }
        let n = match rec.get_info_i32(hdr, b"AC", &mut self.iarr) {
            Ok(n) if n + 1 == rec.allele_count() => n,
            _ => return Ok(()),
        };
        self.farr.clear();
        self.farr.extend(
            self.iarr[..n]
                .iter()
                .map(|&ac| (f64::from(ac) / f64::from(an)) as f32),
        );
        rec.update_info_f32(out, b"AF", &self.farr)
            .map_err(|_| update_error(b"AF", hdr, rec))
    }

    /// Sets FORMAT/VAF and VAF1 from FORMAT/AD (bcftools' `process_vaf_vaf1`).
    fn process_vaf_vaf1(
        &mut self,
        hdr: &Header,
        out: &Header,
        rec: &mut Record,
    ) -> anyhow::Result<()> {
        if !self.tags.has(Tags::VAF | Tags::VAF1) || rec.allele_count() <= 1 {
            return Ok(());
        }
        let niarr = rec.get_format_i32(hdr, b"AD", &mut self.iarr).unwrap_or(0);
        if niarr == 0 {
            if self.warned & (Tags::VAF | Tags::VAF1) == 0 {
                eprintln!(
                    "Warning: cannot add the VAF/VAF1 annotations, the required FORMAT/AD tag is missing at {}:{}.\n         (This warning is printed only once.)",
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1
                );
            }
            self.warned |= Tags::VAF | Tags::VAF1;
            return Ok(());
        }
        let nsmpl = hdr.sample_count();
        if niarr != nsmpl * rec.allele_count() {
            return Ok(());
        }
        if self.tags.has(Tags::VAF) {
            self.process_vaf(hdr, out, rec, VafTag::Vaf)?;
        }
        if self.tags.has(Tags::VAF1) {
            self.process_vaf(hdr, out, rec, VafTag::Vaf1)?;
        }
        Ok(())
    }

    /// Computes VAF or VAF1 for each sample (bcftools' `process_vaf`).
    fn process_vaf(
        &mut self,
        hdr: &Header,
        out: &Header,
        rec: &mut Record,
        tag: VafTag,
    ) -> anyhow::Result<()> {
        let per_alt = tag == VafTag::Vaf;
        let nsmpl = hdr.sample_count();
        let nval = self.iarr.len() / nsmpl;
        let nval1 = if per_alt { rec.allele_count() - 1 } else { 1 };
        self.farr.clear();
        self.farr.resize(nval1 * nsmpl, 0.0);
        for i in 0..nsmpl {
            let src = &self.iarr[i * nval..(i + 1) * nval];
            let dst = &mut self.farr[i * nval1..(i + 1) * nval1];
            let mut sum = 0f32;
            let mut complete = true;
            for &v in src {
                if v == INT32_MISSING || v == INT32_VECTOR_END {
                    complete = false;
                    break;
                }
                sum += v as f32;
            }
            if !complete {
                dst[0] = f32::from_bits(FLOAT_MISSING);
                for d in &mut dst[1..] {
                    *d = f32::from_bits(FLOAT_VECTOR_END);
                }
                continue;
            }
            if per_alt {
                for (j, d) in dst.iter_mut().enumerate() {
                    *d = if sum != 0.0 {
                        src[j + 1] as f32 / sum
                    } else {
                        0.0
                    };
                }
            } else {
                dst[0] = if sum != 0.0 {
                    (sum - src[0] as f32) / sum
                } else {
                    0.0
                };
            }
        }
        let (id, name): (Option<u32>, &[u8]) = match tag {
            VafTag::Vaf => (self.vaf_id, b"VAF"),
            VafTag::Vaf1 => (self.vaf1_id, b"VAF1"),
        };
        let r = match id {
            Some(id) => rec.update_format_f32_id(out, id, &self.farr),
            None => Err(pasta::variant::record::UpdateError::NotInHeader),
        };
        r.map_err(|_| update_error(name, hdr, rec))
    }
}

impl Plugin for FillTags {
    fn process(&mut self, hdr: &Header, out: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        for pop in &mut self.pops {
            for ftf in &mut pop.ftf {
                ftf_filter_expr(ftf, hdr, out, rec)?;
            }
        }
        if self.fmt {
            self.process_fmt(hdr, out, rec)?;
            self.process_info_af(hdr, out, rec)?;
            self.process_vaf_vaf1(hdr, out, rec)?;
        }
        if self.tags.has(Tags::END) {
            // C stores the end in an int32_t (modulo 2^32 on the platforms HTSlib runs on).
            let end = (rec.pos() + rec.rlen()) as i32;
            rec.update_info_i32(out, b"END", &[end]).map_err(|_| {
                anyhow::anyhow!(
                    "Error occurred while updating INFO/END at {}:{}",
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1
                )
            })?;
        }
        if self.tags.has(Tags::TYPE) {
            set_type(hdr, out, rec)?;
        }
        Ok(true)
    }
}

/// Sets INFO/TYPE to the record's variant types.
fn set_type(hdr: &Header, out: &Header, rec: &mut Record) -> anyhow::Result<()> {
    use pasta::variant::VariantType as V;
    rec.unpack(Unpack::STR);
    let ty = rec.variant_types();
    let mut s = String::new();
    if ty.is_empty() {
        s.push_str("REF,");
    }
    for (bit, name) in [
        (V::SNP, "SNP,"),
        (V::MNP, "MNP,"),
        (V::INDEL, "INDEL,"),
        (V::OTHER, "OTHER,"),
        (V::BND, "BND,"),
        (V::OVERLAP, "OVERLAP,"),
    ] {
        if ty.contains(bit) {
            s.push_str(name);
        }
    }
    if s.is_empty() {
        return Ok(());
    }
    s.pop();
    rec.update_info_string(out, b"TYPE", Some(s.as_bytes()))
        .map_err(|_| {
            anyhow::anyhow!(
                "Error occurred while updating INFO/TYPE at {}:{}",
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1
            )
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tag_lists_split_outside_brackets() {
        assert_eq!(
            parse_tag_list("AN,X=sum(FMT/AD[0,1]),AC"),
            ["AN", "X=sum(FMT/AD[0,1])", "AC"]
        );
        assert_eq!(parse_tag_list(""), [""]);
    }

    #[test]
    fn hwe_of_a_balanced_sample() {
        let mut probs = Vec::new();
        // 10 samples: 4 0/0, 4 0/1, 2 1/1 -> nref = 12, nalt = 8, nhet = 4.
        let (hwe, exc) = calc_hwe(&mut probs, 12, 8, 4).expect("valid counts");
        assert!((0.0..=1.0).contains(&hwe), "{hwe}");
        assert!((0.0..=1.0).contains(&exc), "{exc}");
        assert!(calc_hwe(&mut probs, 3, 2, 2).is_err());
    }

    /// Runs `+fill-tags` with `opts` over records given as `ALT GT1 GT2`, returning each
    /// record's INFO `tags` as `TAG=values` (`TAG=` when unset).
    fn fill(opts: &[&str], tags: &[&str], recs: &[&str]) -> Vec<Vec<String>> {
        let mut hdr = Header::parse(
            b"##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n",
        )
        .expect("the header parses");
        let mut out = hdr.dup().expect("the header copies");
        let argv: Vec<String> = std::iter::once("in.vcf")
            .chain(opts.iter().copied())
            .map(String::from)
            .collect();
        let mut plugin = init(&argv, Some((&hdr, &mut out)))
            .expect("the plugin starts")
            .plugin;
        let mut rec = Record::new();
        let (mut ivals, mut fvals) = (Vec::new(), Vec::new());
        let mut values = |rec: &mut Record, tag: &str| -> String {
            let vals: Vec<String> =
                if let Ok(n) = rec.get_info_i32(&out, tag.as_bytes(), &mut ivals) {
                    ivals[..n].iter().map(i32::to_string).collect()
                } else if let Ok(n) = rec.get_info_f32(&out, tag.as_bytes(), &mut fvals) {
                    fvals[..n].iter().map(|v| format!("{v:.4}")).collect()
                } else {
                    Vec::new()
                };
            format!("{tag}={}", vals.join(","))
        };
        recs.iter()
            .map(|r| {
                let [alt, gt1, gt2] = r.split(' ').collect::<Vec<_>>()[..] else {
                    panic!("bad record {r:?}");
                };
                let mut line = format!("1\t5\t.\tA\t{alt}\t.\t.\t.\tGT\t{gt1}\t{gt2}").into_bytes();
                pasta::variant::vcf::parse(&mut hdr, &mut line, &mut rec).expect("the line parses");
                plugin
                    .process(&hdr, &out, &mut rec)
                    .expect("the record is processed");
                tags.iter().map(|t| values(&mut rec, t)).collect()
            })
            .collect()
    }

    #[test]
    fn polyploid_genotypes_count_every_allele() {
        let tags = ["AN", "AC", "AF", "AC_Het", "AC_Hom", "AC_Hemi"];
        let opts = ["-t", "AN,AC,AF,AC_Het,AC_Hom,AC_Hemi"];
        let recs = [
            "G 0/0/1 0/0",
            "G 1/1/1 0/0",
            "G 0/0/1/1 0/1",
            "G 0/1/. 0/1",
            "G 0/0/1 ./.",
            "G 1/1/. 1",
        ];
        let want = [
            "AN=5 AC=1 AF=0.2000 AC_Het=1 AC_Hom=0 AC_Hemi=0",
            "AN=5 AC=3 AF=0.6000 AC_Het=0 AC_Hom=3 AC_Hemi=0",
            "AN=6 AC=3 AF=0.5000 AC_Het=3 AC_Hom=0 AC_Hemi=0",
            "AN=4 AC=2 AF=0.5000 AC_Het=2 AC_Hom=0 AC_Hemi=0",
            "AN=3 AC=1 AF=0.3333 AC_Het=1 AC_Hom=0 AC_Hemi=0",
            "AN=3 AC=3 AF=1.0000 AC_Het=0 AC_Hom=0 AC_Hemi=3",
        ];
        let got: Vec<String> = fill(&opts, &tags, &recs)
            .iter()
            .map(|r| r.join(" "))
            .collect();
        assert_eq!(got, want);
        // -d counts the called alleles of a partly missing genotype as half calls.
        let got: Vec<String> = fill(&["-d", "-t", "AN,AC,AC_Hemi"], &tags[..2], &recs[3..])
            .iter()
            .map(|r| r.join(" "))
            .collect();
        assert_eq!(got, ["AN=4 AC=2", "AN=3 AC=1", "AN=3 AC=3"]);
    }

    #[test]
    fn diploid_and_haploid_counts_are_unchanged() {
        // bcftools 1.24 gives the same values.
        let tags = ["AN", "AC", "AC_Het", "AC_Hom", "AC_Hemi", "HWE", "ExcHet"];
        let opts = ["-t", "AN,AC,AC_Het,AC_Hom,AC_Hemi,HWE,ExcHet"];
        let recs = [
            "G,T 1/2 1/1",
            "G 0/1 1",
            "G ./1 0/0",
            "G 0/1 0/1",
            "G ./. .",
        ];
        let want = [
            "AN=4 AC=3,1 AC_Het=1,1 AC_Hom=2,0 AC_Hemi=0,0 HWE=1.0000,1.0000 ExcHet=1.0000,1.0000",
            "AN=3 AC=2 AC_Het=1 AC_Hom=0 AC_Hemi=1 HWE=1.0000 ExcHet=1.0000",
            "AN=3 AC=1 AC_Het=0 AC_Hom=0 AC_Hemi=1 HWE=1.0000 ExcHet=1.0000",
            "AN=4 AC=2 AC_Het=2 AC_Hom=0 AC_Hemi=0 HWE=1.0000 ExcHet=0.6667",
            "AN=0 AC=0 AC_Het=0 AC_Hom=0 AC_Hemi=0 HWE=1.0000 ExcHet=1.0000",
        ];
        let got: Vec<String> = fill(&opts, &tags, &recs)
            .iter()
            .map(|r| r.join(" "))
            .collect();
        assert_eq!(got, want);
    }

    #[test]
    fn hwe_leaves_polyploid_genotypes_out() {
        // Counting 1/1/1 as 1/1 next to a 0/0 gives HWE 1/3; with the diploid 0/0 alone
        // there is no ALT allele to test. Next to 0/1, 0/0/1 is not a second 0/1 either.
        let got = fill(
            &["-t", "HWE,ExcHet"],
            &["HWE", "ExcHet"],
            &["G 1/1/1 0/0", "G 0/0/1 0/1"],
        );
        assert_eq!(
            got,
            [
                ["HWE=1.0000", "ExcHet=1.0000"],
                ["HWE=1.0000", "ExcHet=1.0000"]
            ]
        );
    }

    #[test]
    fn int_conversion_is_x86_truncation() {
        assert_eq!(int32_from_double(2.5), 3);
        assert_eq!(int32_from_double(-2.5), -3);
        assert_eq!(int32_from_double(2_147_483_647.0), i32::MAX);
        assert_eq!(int32_from_double(2_147_483_647.5), i32::MIN);
        assert_eq!(int32_from_double(f64::INFINITY), i32::MIN);
        assert_eq!(int32_from_double(f64::NEG_INFINITY), i32::MIN);
        assert_eq!(int32_from_double(f64::NAN), i32::MIN);
        assert_eq!(int32_from_double(-2_147_483_648.4), i32::MIN);
    }
}
