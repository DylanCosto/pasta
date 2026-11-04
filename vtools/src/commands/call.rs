//! `vtools call`: SNP/indel calling from `mpileup` output (bcftools call, `vcfcall.c`).

mod ccall;
mod em;
mod f80;
mod gvcf;
mod kfunc;
mod mcall;
mod parallel;
mod ploidy;
mod prob1;
mod regidx;
mod samples;
mod targets;
pub(crate) mod vcmp;

use crate::options::last_index;
use std::{
    io::{self, IsTerminal},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    GetError, Header, HeaderRecordType, Number, Record, Unpack, UpdateError, ValueType,
    encoding::{INT32_VECTOR_END, float_vector_end},
    io::{OutputType, WriteError, Writer},
    readlist::ListSource,
    synced::{Overlap, Regions, SyncedError, SyncedOptions, SyncedReader},
};

use self::{
    ccall::{CCall, I16Batch},
    gvcf::{CalledSite, Gvcf, GvcfItem, GvcfRanges},
    mcall::MCall,
    ploidy::Ploidy,
    samples::SampleSet,
    targets::{TargetKind, Targets},
};
use crate::{commands::head::apply_verbosity, output, version};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF file from `bcftools mpileup` (`-` or none for stdin); further
    /// arguments are ignored, as bcftools does.
    inputs: Vec<String>,

    /// Do not append version and command line to the header.
    #[arg(long = "no-version")]
    no_version: bool,

    /// Write output to a file [standard output].
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,

    /// u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "u|b|v|z[0-9]", action = clap::ArgAction::Append)]
    output_type: Vec<String>,

    /// Predefined ploidy, 'list' to print available settings, append '?' for details [2].
    #[arg(long = "ploidy", value_name = "ASSEMBLY[?]", action = clap::ArgAction::Append)]
    ploidy: Vec<String>,

    /// Space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY.
    #[arg(long = "ploidy-file", value_name = "FILE")]
    ploidy_file: Option<String>,

    /// Restrict to comma-separated list of regions.
    #[arg(short = 'r', long = "regions", value_name = "REGION", action = clap::ArgAction::Append)]
    regions: Vec<String>,

    /// Restrict to regions listed in a file.
    #[arg(short = 'R', long = "regions-file", value_name = "FILE", action = clap::ArgAction::Append)]
    regions_file: Vec<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(
        long = "regions-overlap",
        value_name = "0|1|2",
        action = clap::ArgAction::Append
    )]
    regions_overlap: Vec<String>,

    /// List of samples to include [all samples].
    #[arg(short = 's', long = "samples", value_name = "LIST", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    samples: Vec<String>,

    /// PED file or a file with an optional column with sex [all samples].
    #[arg(short = 'S', long = "samples-file", value_name = "FILE", action = clap::ArgAction::Append)]
    samples_file: Vec<String>,

    /// Similar to -r but streams rather than index-jumps.
    #[arg(short = 't', long = "targets", value_name = "REGION", action = clap::ArgAction::Append)]
    targets: Vec<String>,

    /// Similar to -R but streams rather than index-jumps.
    #[arg(short = 'T', long = "targets-file", value_name = "FILE", action = clap::ArgAction::Append)]
    targets_file: Vec<String>,

    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Keep all possible alternate alleles at variant sites.
    #[arg(short = 'A', long = "keep-alts")]
    keep_alts: bool,

    /// Keep the unobserved allele <*> or <NON_REF>.
    #[arg(short = '*', long = "keep-unseen-allele")]
    keep_unseen_allele: bool,

    /// Optional tags to output (lowercase allowed); '?' to list available tags.
    #[arg(short = 'a', long = "annotate", value_name = "LIST", action = clap::ArgAction::Append)]
    annotate: Vec<String>,

    /// Deprecated alias of -a.
    #[arg(short = 'f', long = "format-fields", value_name = "LIST", action = clap::ArgAction::Append)]
    format_fields: Vec<String>,

    /// Use prior allele frequencies, determined from these pre-filled tags.
    #[arg(
        short = 'F',
        long = "prior-freqs",
        value_name = "AN,AC",
        action = clap::ArgAction::Append
    )]
    prior_freqs: Vec<String>,

    /// Group samples by population (file with "sample\tgroup") or "-" for single-sample
    /// calling.
    #[arg(
        short = 'G',
        long = "group-samples",
        value_name = "FILE|-",
        allow_hyphen_values = true
    )]
    group_samples: Option<String>,

    /// The tag to use with -G, by default FORMAT/QS and FORMAT/AD are checked automatically.
    #[arg(long = "group-samples-tag", value_name = "TAG")]
    group_samples_tag: Option<String>,

    /// Group non-variant sites into gVCF blocks by minimum per-sample DP.
    #[arg(
        short = 'g',
        long = "gvcf",
        value_name = "INT,[...]",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    gvcf: Vec<String>,

    /// Output also sites missed by mpileup but present in -T.
    #[arg(short = 'i', long = "insert-missed")]
    insert_missed: bool,

    /// Skip sites where REF is N (the default).
    #[arg(short = 'N', long = "skip-Ns", action = clap::ArgAction::Count)]
    skip_ns: u8,

    /// Keep sites with masked reference allele (REF=N).
    #[arg(short = 'M', long = "keep-masked-refs", visible_alias = "keep-masked-ref", action = clap::ArgAction::Count)]
    keep_masked_ref: u8,

    /// Skip indels/snps.
    #[arg(short = 'V', long = "skip-variants", value_name = "TYPE", action = clap::ArgAction::Append)]
    skip_variants: Vec<String>,

    /// Output variant sites only.
    #[arg(short = 'v', long = "variants-only")]
    variants_only: bool,

    /// Verbosity level.
    #[arg(
        long = "verbosity",
        value_name = "INT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    verbosity: Vec<String>,

    /// Automatically index the output files [off].
    #[arg(
        short = 'W',
        long = "write-index",
        value_name = "FMT",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "",
        action = clap::ArgAction::Append
    )]
    write_index: Vec<String>,

    /// The original calling method (conflicts with -m).
    #[arg(short = 'c', long = "consensus-caller")]
    consensus_caller: bool,

    /// One of: alleles, trio (see manual).
    #[arg(short = 'C', long = "constrain", value_name = "STR", action = clap::ArgAction::Append)]
    constrain: Vec<String>,

    /// Alternative model for multiallelic and rare-variant calling (conflicts with -c).
    #[arg(short = 'm', long = "multiallelic-caller")]
    multiallelic_caller: bool,

    /// Likelihood of novel mutation for constrained trio calling [1e-8,1e-9,1e-9].
    #[arg(
        short = 'n',
        long = "novel-rate",
        value_name = "FLOAT,[...]",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    novel_rate: Vec<String>,

    /// Variant if P(ref|D)<FLOAT with -c [0.5].
    #[arg(
        short = 'p',
        long = "pval-threshold",
        value_name = "FLOAT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    pval_threshold: Vec<String>,

    /// Mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3].
    #[arg(
        short = 'P',
        long = "prior",
        value_name = "FLOAT",
        allow_hyphen_values = true,
        action = clap::ArgAction::Append
    )]
    prior: Vec<String>,

    /// Deprecated: --ploidy X.
    #[arg(short = 'X', long = "chromosome-X", action = clap::ArgAction::Count)]
    chromosome_x: u8,

    /// Deprecated: --ploidy Y.
    #[arg(short = 'Y', long = "chromosome-Y", action = clap::ArgAction::Count)]
    chromosome_y: u8,

    /// Print help.
    #[arg(short = 'h', long = "help", action = clap::ArgAction::Count)]
    help: u8,

    /// Print help.
    #[arg(short = '?', hide = true, action = clap::ArgAction::Count)]
    help_q: u8,
}

/// The `CALL_*` flags that change how sites are called.
#[derive(Clone, Copy, Debug, Default)]
pub struct CallFlags {
    /// `-A`: keep all alternate alleles at variant sites.
    pub keep_alts: bool,
    /// `-v`: variant sites only.
    pub var_only: bool,
    /// `-C trio`.
    pub constrain_trio: bool,
    /// `-C alleles`.
    pub constrain_alleles: bool,
    /// `-*`: keep the unseen allele.
    pub keep_unseen: bool,
}

/// The optional annotations of `-a`.
#[derive(Clone, Copy, Debug, Default)]
pub struct OutputTags {
    pub gq: bool,
    pub gp: bool,
    pub pv4: bool,
}

/// Settings shared by the callers (the read-only part of bcftools' `call_t`).
#[derive(Clone, Debug, Default)]
pub struct Shared {
    pub flags: CallFlags,
    pub tags: OutputTags,
    /// `-p`: the P(ref|D) threshold of `-c`.
    pub pref: f64,
    /// `-F AN,AC`.
    pub prior_an: Option<Vec<u8>>,
    pub prior_ac: Option<Vec<u8>>,
}

/// The calling model.
#[derive(Clone, Debug)]
pub enum Caller {
    /// `-m`.
    Multi(Box<MCall>),
    /// `-c`.
    Consensus(Box<CCall>),
}

impl Caller {
    /// Calls one site (bcftools' `mcall` or `ccall`).
    ///
    /// # Errors
    ///
    /// Returns bcftools' fatal errors (missing PL or QS, incompatible alleles, ...).
    pub fn call(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        shared: &Shared,
        ploidy: &[u8],
        unseen: usize,
        tgt_als: Option<&[Vec<u8>]>,
    ) -> anyhow::Result<i32> {
        match self {
            Self::Multi(m) => m.call(hdr, rec, shared, ploidy, unseen, tgt_als),
            Self::Consensus(c) => c.call(hdr, rec, shared, ploidy, unseen),
        }
    }

    /// Starts a batch of records on a worker thread: with `-c`, records without INFO/I16
    /// take `anno16` until one reads its own (`-m` keeps nothing from record to record).
    pub fn start_batch(&mut self, anno16: &[f32]) {
        if let Self::Consensus(c) = self {
            c.start_batch(anno16);
        }
    }

    /// What the batch since [`Self::start_batch`] did with the I16 values of `-c`.
    pub fn end_batch(&self) -> Option<I16Batch> {
        match self {
            Self::Multi(_) => None,
            Self::Consensus(c) => Some(c.end_batch()),
        }
    }
}

/// Maps a getter result to HTSlib's return value: the count, or -1 (tag not defined), -2
/// (wrong type), -3 (absent).
pub fn get_code(r: Result<usize, GetError>) -> i32 {
    match r {
        Ok(n) => n as i32,
        Err(GetError::NotInHeader) => -1,
        Err(GetError::WrongType) => -2,
        Err(_) => -3,
    }
}

/// A tag of the output header resolved to its dictionary ID once, so that per-record reads
/// and writes skip the name lookup. Each method returns what its name-based counterpart
/// returns, including for a tag the header does not define.
#[derive(Clone, Copy, Debug)]
pub struct Tag(Option<u32>);

impl Tag {
    /// Resolves `name` in `hdr`.
    pub fn new(hdr: &Header, name: &[u8]) -> Self {
        Self(hdr.id(name))
    }

    /// `bcf_get_info_int32`, as [`get_code`] reports it.
    pub fn get_info_i32(self, hdr: &Header, rec: &mut Record, out: &mut Vec<i32>) -> i32 {
        self.0
            .map_or(-1, |id| get_code(rec.get_info_i32_id(hdr, id, out)))
    }

    /// `bcf_get_info_float`, as [`get_code`] reports it.
    pub fn get_info_f32(self, hdr: &Header, rec: &mut Record, out: &mut Vec<f32>) -> i32 {
        self.0
            .map_or(-1, |id| get_code(rec.get_info_f32_id(hdr, id, out)))
    }

    /// `bcf_get_format_int32`, as [`get_code`] reports it.
    pub fn get_format_i32(self, hdr: &Header, rec: &mut Record, out: &mut Vec<i32>) -> i32 {
        self.0
            .map_or(-1, |id| get_code(rec.get_format_i32_id(hdr, id, out)))
    }

    /// `bcf_update_info_int32` (an empty slice removes the field).
    pub fn update_info_i32(
        self,
        hdr: &Header,
        rec: &mut Record,
        values: &[i32],
    ) -> Result<(), UpdateError> {
        let id = self.0.ok_or(UpdateError::NotInHeader)?;
        rec.update_info_i32_id(hdr, id, values)
    }

    /// `bcf_update_info_float` (an empty slice removes the field).
    pub fn update_info_f32(
        self,
        hdr: &Header,
        rec: &mut Record,
        values: &[f32],
    ) -> Result<(), UpdateError> {
        let id = self.0.ok_or(UpdateError::NotInHeader)?;
        rec.update_info_f32_id(hdr, id, values)
    }

    /// Removes an INFO field (`bcf_update_info_int32` with no values); a tag the header
    /// does not define is ignored.
    pub fn remove_info(self, hdr: &Header, rec: &mut Record) {
        if let Some(id) = self.0 {
            // Fails only for a tag that is not an INFO tag, which leaves nothing to remove.
            let _ = rec.remove_info_id(hdr, id);
        }
    }

    /// `bcf_update_format_int32` (an empty slice removes the field).
    pub fn update_format_i32(
        self,
        hdr: &Header,
        rec: &mut Record,
        values: &[i32],
    ) -> Result<(), UpdateError> {
        match self.0 {
            Some(id) => rec.update_format_i32_id(hdr, id, values),
            None if values.is_empty() => Ok(()),
            None => Err(UpdateError::NotInHeader),
        }
    }

    /// `bcf_update_format_float` (an empty slice removes the field).
    pub fn update_format_f32(
        self,
        hdr: &Header,
        rec: &mut Record,
        values: &[f32],
    ) -> Result<(), UpdateError> {
        match self.0 {
            Some(id) => rec.update_format_f32_id(hdr, id, values),
            None if values.is_empty() => Ok(()),
            None => Err(UpdateError::NotInHeader),
        }
    }

    /// Removes a FORMAT field; a tag the header does not define is ignored.
    pub fn remove_format(self, hdr: &Header, rec: &mut Record) {
        if let Some(id) = self.0 {
            // Never fails for removal.
            let _ = rec.remove_format_id(hdr, id);
        }
    }
}

/// HTSlib's `bcf_is_snp`: every allele is one base (not `*`), or mpileup's `<*>`/`<X>`.
pub fn is_snp(rec: &Record) -> bool {
    rec.peek_alleles().all(is_snp_allele)
}

/// Whether an allele passes `bcf_is_snp`'s test.
fn is_snp_allele(a: &[u8]) -> bool {
    (a.len() == 1 && a[0] != b'*') || a.starts_with(b"<X>") || a.starts_with(b"<*>")
}

/// The Number=R INFO and FORMAT tags of the output header, and buffers for trimming them
/// (bcftools' `mcall_trim_and_update_numberR`).
#[derive(Clone, Debug, Default)]
pub struct NumberR {
    /// The dictionary IDs and types of the Number=R INFO tags.
    info: Vec<(u32, ValueType)>,
    /// The same for FORMAT tags.
    fmt: Vec<(u32, ValueType)>,
    selected: Vec<usize>,
    ints: Vec<i32>,
    floats: Vec<f32>,
    out_ints: Vec<i32>,
    out_floats: Vec<f32>,
}

impl NumberR {
    /// Finds the Number=R Integer and Float tags of a header.
    pub fn new(hdr: &Header) -> Self {
        let tags = |kind: HeaderRecordType| -> Vec<(u32, ValueType)> {
            let mut out = Vec::new();
            for id in 0..hdr.id_count() as u32 {
                if hdr.id_number(kind, id) != Some(Number::R) || hdr.id_name(id).is_none() {
                    continue;
                }
                if let Some(ty @ (ValueType::Integer | ValueType::Float)) = hdr.id_type(kind, id) {
                    out.push((id, ty));
                }
            }
            out
        };
        Self {
            info: tags(HeaderRecordType::Info),
            fmt: tags(HeaderRecordType::Format),
            ..Self::default()
        }
    }

    /// Drops the values of removed alleles from the Number=R fields, given the old-to-new
    /// allele map.
    ///
    /// # Errors
    ///
    /// Returns an error when a field cannot be updated.
    pub fn trim(
        &mut self,
        hdr: &Header,
        rec: &mut Record,
        als_map: &[i32],
        nals_ori: usize,
        nals_new: usize,
    ) -> anyhow::Result<()> {
        if nals_ori == nals_new {
            return Ok(());
        }
        rec.unpack(Unpack::INFO | Unpack::FMT);
        let mut selected = std::mem::take(&mut self.selected);
        selected.clear();
        for f in rec.info_fields().iter().filter(|f| !f.is_removed()) {
            if let Some(k) = self.info.iter().position(|t| t.0 == f.key() as u32) {
                selected.push(k);
            }
        }
        for &k in &selected {
            let (id, ty) = self.info[k];
            if ty == ValueType::Integer {
                if get_code(rec.get_info_i32_id(hdr, id, &mut self.ints)) <= 0 {
                    continue;
                }
                remap_site(
                    &self.ints,
                    als_map,
                    nals_new,
                    INT32_VECTOR_END,
                    &mut self.out_ints,
                );
                rec.update_info_i32_id(hdr, id, &self.out_ints)?;
            } else {
                if get_code(rec.get_info_f32_id(hdr, id, &mut self.floats)) <= 0 {
                    continue;
                }
                let pad = float_vector_end();
                remap_site(&self.floats, als_map, nals_new, pad, &mut self.out_floats);
                rec.update_info_f32_id(hdr, id, &self.out_floats)?;
            }
        }
        selected.clear();
        for f in rec.format_fields().iter().filter(|f| !f.is_removed()) {
            if let Some(k) = self.fmt.iter().position(|t| t.0 == f.id() as u32) {
                selected.push(k);
            }
        }
        let nsmpl = hdr.sample_count();
        for &k in &selected {
            let (id, ty) = self.fmt[k];
            if ty == ValueType::Integer {
                let n = get_code(rec.get_format_i32_id(hdr, id, &mut self.ints));
                if n <= 0 {
                    continue;
                }
                if n as usize != nals_ori * nsmpl {
                    bail!("assertion failed: nret==nals_ori*nsmpl");
                }
                let pad = INT32_VECTOR_END;
                remap_samples(
                    &self.ints,
                    als_map,
                    nals_ori,
                    nals_new,
                    pad,
                    &mut self.out_ints,
                );
                rec.update_format_i32_id(hdr, id, &self.out_ints)?;
            } else {
                let n = get_code(rec.get_format_f32_id(hdr, id, &mut self.floats));
                if n <= 0 {
                    continue;
                }
                if n as usize != nals_ori * nsmpl {
                    bail!("assertion failed: nret==nals_ori*nsmpl");
                }
                let pad = float_vector_end();
                remap_samples(
                    &self.floats,
                    als_map,
                    nals_ori,
                    nals_new,
                    pad,
                    &mut self.out_floats,
                );
                rec.update_format_f32_id(hdr, id, &self.out_floats)?;
            }
        }
        self.selected = selected;
        Ok(())
    }
}

/// Keeps the values of the kept alleles of a site field (a single REF value when only REF
/// is left).
fn remap_site<T: Copy>(values: &[T], als_map: &[i32], nals_new: usize, pad: T, out: &mut Vec<T>) {
    out.clear();
    if nals_new == 1 {
        out.extend_from_slice(&values[..1]);
        return;
    }
    out.resize(nals_new, pad);
    for (j, &k) in als_map.iter().enumerate() {
        // An index past the kept alleles (the unseen allele when it is among the best alleles
        // but not on output) is a write past the values bcftools updates the field with.
        if let (Ok(k), Some(&v)) = (usize::try_from(k), values.get(j))
            && let Some(o) = out.get_mut(k)
        {
            *o = v;
        }
    }
}

/// Keeps the values of the kept alleles of each sample.
fn remap_samples<T: Copy>(
    values: &[T],
    als_map: &[i32],
    nals_ori: usize,
    nals_new: usize,
    pad: T,
    out: &mut Vec<T>,
) {
    let nsmpl = values.len() / nals_ori;
    out.clear();
    out.resize(nsmpl * nals_new, pad);
    for s in 0..nsmpl {
        for (k, &l) in als_map.iter().enumerate().take(nals_ori) {
            // An index past the kept alleles (see `remap_site`) writes over the next sample's
            // first value in bcftools, which that sample then sets, or past the end.
            if let Ok(l) = usize::try_from(l)
                && l < nals_new
            {
                out[s * nals_new + l] = values[s * nals_ori + k];
            }
        }
    }
}

/// bcftools' `parse_output_tags`.
fn parse_output_tags(s: &str, tags: &mut OutputTags) -> anyhow::Result<()> {
    // The C loop stops at the end of the string, so an empty string or a trailing comma
    // adds no item, while an empty item before a comma matches GQ (a zero-length prefix).
    let mut rest = s;
    while !rest.is_empty() {
        let (item, tail) = rest.split_once(',').unwrap_or((rest, ""));
        rest = tail;
        // strncasecmp(ss, NAME, len): the item is a case-insensitive prefix of the name.
        let is =
            |name: &str| item.len() <= name.len() && name[..item.len()].eq_ignore_ascii_case(item);
        if is("GQ") || is("FORMAT/GQ") || is("FMT/GQ") {
            tags.gq = true;
        } else if is("GP") || is("FORMAT/GP") || is("FMT/GP") {
            tags.gp = true;
        } else if is("PV4") || is("INFO/PV4") {
            tags.pv4 = true;
        } else {
            eprintln!("Could not parse \"{s}\"");
            std::process::exit(1);
        }
    }
    Ok(())
}

fn list_annotations() {
    eprint!(
        "\n\
         Optional INFO annotations available with -m (\"INFO/\" prefix is optional):\n  \
         INFO/PV4   .. P-values for strand bias, baseQ bias, mapQ bias and tail distance bias (Number=4,Type=Float)\n\
         \n\
         Optional FORMAT annotations available with -m (\"FORMAT/\" prefix is optional):\n  \
         FORMAT/GQ  .. Phred-scaled genotype quality (Number=1,Type=Integer)\n  \
         FORMAT/GP  .. Phred-scaled genotype posterior probabilities (Number=G,Type=Float)\n\
         \n"
    );
}

/// bcftools call's usage text.
const USAGE: &str = r#"
About:   SNP/indel variant calling from VCF/BCF. To be used in conjunction with vtools mpileup.
         This command replaces the former "bcftools view" caller. Some of the original
         functionality has been temporarily lost in the process of transition to htslib,
         but will be added back on popular demand. The original calling model can be
         invoked with the -c option.
Usage:   vtools call [options] <in.vcf.gz>

File format options:
       --no-version                Do not append version and command line to the header
   -o, --output FILE               Write output to a file [standard output]
   -O, --output-type b|u|z|v       Output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
       --ploidy ASSEMBLY[?]        Predefined ploidy, 'list' to print available settings, append '?' for details [2]
       --ploidy-file FILE          Space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY
   -r, --regions REGION            Restrict to comma-separated list of regions
   -R, --regions-file FILE         Restrict to regions listed in a file
       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
   -s, --samples LIST              List of samples to include [all samples]
   -S, --samples-file FILE         PED file or a file with an optional column with sex (see man page for details) [all samples]
   -t, --targets REGION            Similar to -r but streams rather than index-jumps
   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps
       --threads INT               Use multithreading with INT worker threads [0]

Input/output options:
   -A, --keep-alts                 Keep all possible alternate alleles at variant sites
   -*, --keep-unseen-allele        Keep the unobserved allele <*> or <NON_REF>
   -a, --annotate LIST             Optional tags to output (lowercase allowed); '?' to list available tags
   -F, --prior-freqs AN,AC         Use prior allele frequencies, determined from these pre-filled tags
   -G, --group-samples FILE|-      Group samples by population (file with "sample\tgroup") or "-" for single-sample calling.
                                   This requires FORMAT/QS or other Number=R,Type=Integer tag such as FORMAT/AD
       --group-samples-tag TAG     The tag to use with -G, by default FORMAT/QS and FORMAT/AD are checked automatically
   -g, --gvcf INT,[...]            Group non-variant sites into gVCF blocks by minimum per-sample DP
   -i, --insert-missed             Output also sites missed by mpileup but present in -T
   -M, --keep-masked-ref           Keep sites with masked reference allele (REF=N)
   -V, --skip-variants TYPE        Skip indels/snps
   -v, --variants-only             Output variant sites only
       --verbosity INT             Verbosity level
   -W, --write-index[=FMT]         Automatically index the output files [off]

Consensus/variant calling options:
   -c, --consensus-caller          The original calling method (conflicts with -m)
   -C, --constrain STR             One of: alleles, trio (see manual)
   -m, --multiallelic-caller       Alternative model for multiallelic and rare-variant calling (conflicts with -c)
   -n, --novel-rate FLOAT,[...]    Likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]
   -p, --pval-threshold FLOAT      Variant if P(ref|D)<FLOAT with -c [0.5]
   -P, --prior FLOAT               Mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]

Example:
   # See also http://samtools.github.io/bcftools/howtos/variant-calling.html
   vtools mpileup -Ou -f reference.fa alignments.bam | vtools call -mv -Ob -o calls.bcf

"#;

/// Prints bcftools call's usage and exits as it does.
fn usage() -> ! {
    eprint!("{USAGE}");
    std::process::exit(255);
}

/// C's `strtod` over the whole string: `None` when anything is left over.
fn strtod_all(s: &str) -> Option<f64> {
    let (v, end) = pasta::variant::number::strtod(s.as_bytes(), 0);
    (end == s.len()).then_some(v)
}

/// bcftools' `parse_novel_rate`: only validated, as trio calling is disabled in bcftools.
fn parse_novel_rate(s: &str) -> anyhow::Result<()> {
    let (_, end) = pasta::variant::number::strtod(s.as_bytes(), 0);
    if end == 0 {
        bail!("Could not parse --novel-rate {s}");
    }
    Ok(())
}

/// The values of a repeated option with their command-line positions.
fn occurrences<'a>(matches: &ArgMatches, id: &str, values: &'a [String]) -> Vec<(usize, &'a str)> {
    matches
        .indices_of(id)
        .map(|idx| idx.zip(values.iter().map(String::as_str)).collect())
        .unwrap_or_default()
}

/// A `-x VALUE` / `-X FILE` list: the value, and whether it names a file.
type ListSpec = (String, ListSource);

/// The last of a `-x VALUE` / `-X FILE` pair by position, read as a file when any `-X` was
/// given (bcftools sets its `*_is_file` flag for good once it sees the file option).
fn last_of_pair(
    matches: &ArgMatches,
    inline: (&str, &[String]),
    file: (&str, &[String]),
) -> Option<ListSpec> {
    let mut all = occurrences(matches, inline.0, inline.1);
    let files = occurrences(matches, file.0, file.1);
    let source = if files.is_empty() {
        ListSource::Inline
    } else {
        ListSource::File
    };
    all.extend(files);
    all.sort_by_key(|(i, _)| *i);
    all.last().map(|&(_, v)| (v.to_owned(), source))
}

/// Whether streamed `-t` regions are the sites to call or the sites to skip.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TargetsLogic {
    Include,
    Exclude,
}

/// `-t` regions the reader streams, to call or to skip.
type StreamTargets = (Regions, TargetsLogic);

/// What the main loop needs from the options.
struct Config {
    shared: Shared,
    acgt_only: bool,
    indel_only: bool,
    no_indel: bool,
    insert_missed: bool,
    consensus: bool,
    threads: usize,
}

/// An option that bcftools call's option loop checks or reports on as it reads it.
#[derive(Clone, Copy, Debug)]
enum LoopOption<'a> {
    Help,
    OutputType(&'a str),
    Verbosity(&'a str),
    /// `-a`, or the deprecated `-f`, which warns first.
    Annotate(&'a str, Deprecated),
    PriorFreqs(&'a str),
    Gvcf(&'a str),
    Constrain(&'a str),
    SkipVariants(&'a str),
    PvalThreshold(&'a str),
    Prior(&'a str),
    NovelRate(&'a str),
    RegionsOverlap(&'a str),
    WriteIndex(&'a str),
    /// `--ploidy`, or the deprecated `-X` and `-Y`, which warn.
    Ploidy(&'a str, Deprecated),
}

/// Whether an option is the deprecated spelling of another.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Deprecated {
    Yes,
    No,
}

/// What the checked options of bcftools call's option loop set.
struct LoopOptions {
    shared: Shared,
    /// `-P`.
    theta: f64,
    gvcf: Option<Gvcf>,
    /// `-V snps`.
    indel_only: bool,
    /// `-V indels`.
    no_indel: bool,
    regions_overlap: Overlap,
    /// The last of --ploidy, -X and -Y.
    ploidy_alias: Option<String>,
}

impl LoopOptions {
    /// The checked options in command-line order.
    fn ordered<'a>(args: &'a Args, matches: &ArgMatches) -> Vec<(usize, LoopOption<'a>)> {
        let mut out = Vec::new();
        let mut add = |id: &str, values: &'a [String], make: fn(&'a str) -> LoopOption<'a>| {
            out.extend(
                occurrences(matches, id, values)
                    .into_iter()
                    .map(|(i, v)| (i, make(v))),
            );
        };
        add("output_type", &args.output_type, LoopOption::OutputType);
        add("verbosity", &args.verbosity, LoopOption::Verbosity);
        add("annotate", &args.annotate, |v| {
            LoopOption::Annotate(v, Deprecated::No)
        });
        add("format_fields", &args.format_fields, |v| {
            LoopOption::Annotate(v, Deprecated::Yes)
        });
        add("prior_freqs", &args.prior_freqs, LoopOption::PriorFreqs);
        add("gvcf", &args.gvcf, LoopOption::Gvcf);
        add("constrain", &args.constrain, LoopOption::Constrain);
        add(
            "skip_variants",
            &args.skip_variants,
            LoopOption::SkipVariants,
        );
        add(
            "pval_threshold",
            &args.pval_threshold,
            LoopOption::PvalThreshold,
        );
        add("prior", &args.prior, LoopOption::Prior);
        add("novel_rate", &args.novel_rate, LoopOption::NovelRate);
        add(
            "regions_overlap",
            &args.regions_overlap,
            LoopOption::RegionsOverlap,
        );
        add("write_index", &args.write_index, LoopOption::WriteIndex);
        add("ploidy", &args.ploidy, |v| {
            LoopOption::Ploidy(v, Deprecated::No)
        });
        for (id, option) in [
            ("chromosome_x", LoopOption::Ploidy("X", Deprecated::Yes)),
            ("chromosome_y", LoopOption::Ploidy("Y", Deprecated::Yes)),
            ("help", LoopOption::Help),
            ("help_q", LoopOption::Help),
        ] {
            if matches.get_count(id) > 0
                && let Some(idx) = matches.indices_of(id)
            {
                out.extend(idx.map(|i| (i, option)));
            }
        }
        out.sort_by_key(|o| o.0);
        out
    }

    /// Reads the checked options in command-line order, stopping at the first error as
    /// bcftools does; `None` when `-a ?` listed the annotations (exit status 1).
    fn parse(args: &Args, matches: &ArgMatches) -> anyhow::Result<Option<Self>> {
        let mut o = Self {
            shared: Shared {
                pref: 0.5,
                ..Shared::default()
            },
            theta: 1.1e-3,
            gvcf: None,
            indel_only: false,
            no_indel: false,
            regions_overlap: Overlap::Record,
            ploidy_alias: None,
        };
        for (_, option) in Self::ordered(args, matches) {
            match option {
                LoopOption::Help => usage(),
                LoopOption::OutputType(v) => {
                    output::requested_type(&[output::OutputSetting::Type(v.to_owned())])?;
                }
                LoopOption::Verbosity(v) => apply_verbosity(v)?,
                LoopOption::Annotate(v, deprecated) => {
                    if deprecated == Deprecated::Yes {
                        eprintln!(
                            "Warning: -f, --format-fields will be deprecated, please use -a, --annotate instead."
                        );
                    }
                    if v.starts_with('?') {
                        list_annotations();
                        return Ok(None);
                    }
                    parse_output_tags(v, &mut o.shared.tags)?;
                }
                LoopOption::PriorFreqs(f) => {
                    let Some((an, ac)) = f.split_once(',') else {
                        bail!("Expected two tags with -F (e.g. AN,AC), got \"{f}\"");
                    };
                    o.shared.prior_an = Some(an.as_bytes().to_vec());
                    o.shared.prior_ac = Some(ac.as_bytes().to_vec());
                }
                LoopOption::Gvcf(g) => match GvcfRanges::new(g) {
                    Some(g) => o.gvcf = Some(Gvcf::new(g)),
                    None => bail!("Could not parse: --gvcf {g}"),
                },
                LoopOption::Constrain(c) => {
                    if c.eq_ignore_ascii_case("alleles") {
                        o.shared.flags.constrain_alleles = true;
                    } else if c.eq_ignore_ascii_case("trio") {
                        o.shared.flags.constrain_trio = true;
                    } else {
                        bail!("Unknown argument to -C: \"{c}\"");
                    }
                }
                LoopOption::SkipVariants(v) => {
                    if v.eq_ignore_ascii_case("snps") {
                        o.indel_only = true;
                    } else if v.eq_ignore_ascii_case("indels") {
                        o.no_indel = true;
                    } else {
                        bail!(
                            "Unknown skip category \"{v}\" (-S argument must be \"snps\" or \"indels\")"
                        );
                    }
                }
                LoopOption::PvalThreshold(p) => {
                    o.shared.pref = strtod_all(p)
                        .with_context(|| format!("Could not parse: --pval-threshold {p}"))?;
                }
                LoopOption::Prior(p) => {
                    o.theta = strtod_all(p).with_context(|| {
                        format!("Could not parse, expected float argument: -P {p}")
                    })?;
                }
                LoopOption::NovelRate(n) => parse_novel_rate(n)?,
                LoopOption::RegionsOverlap(s) => {
                    o.regions_overlap = Overlap::parse(s)
                        .with_context(|| format!("Could not parse: --regions-overlap {s}"))?;
                }
                LoopOption::WriteIndex(w) => {
                    output::parse_write_index(w)?;
                }
                LoopOption::Ploidy(alias, deprecated) => {
                    if deprecated == Deprecated::Yes {
                        eprintln!(
                            "Warning: -{alias} will be deprecated, please use --ploidy instead."
                        );
                    }
                    o.ploidy_alias = Some(alias.to_owned());
                }
            }
        }
        Ok(Some(o))
    }
}

pub fn run(args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let LoopOptions {
        mut shared,
        theta,
        gvcf,
        indel_only,
        no_indel,
        regions_overlap,
        ploidy_alias,
    } = match LoopOptions::parse(&args, matches)? {
        Some(o) => o,
        None => return Ok(ExitCode::from(1)),
    };
    // -M and -N: the last one wins.
    let acgt_only = {
        let m = last_index(matches, "keep_masked_ref");
        let n = last_index(matches, "skip_ns");
        match (args.keep_masked_ref > 0, args.skip_ns > 0) {
            (true, true) => n > m,
            (true, false) => false,
            _ => true,
        }
    };
    shared.flags.keep_alts = args.keep_alts;
    shared.flags.keep_unseen = args.keep_unseen_allele;
    shared.flags.var_only = args.variants_only;
    // `int n_threads = strtol(optarg, 0, 0)`; HTSlib ignores a count below one.
    let threads = output::threads(args.threads.as_deref()).unwrap_or(0);

    let ploidy = match (&args.ploidy_file, &ploidy_alias) {
        (Some(f), _) => Some(Ploidy::from_file(f, 2)),
        (None, Some(a)) => Some(Ploidy::predefined(a)),
        (None, None) => None,
    };

    let input = match args.inputs.first() {
        Some(i) => i.clone(),
        None => {
            if io::stdin().is_terminal() {
                usage();
            }
            "-".to_owned()
        }
    };

    if gvcf.is_none() && !shared.flags.var_only && args.multiallelic_caller {
        eprintln!(
            "Warning: `vtools call -m` emits all callable sites, not just variants; use -v for variants only"
        );
    }
    let samples_spec = last_of_pair(
        matches,
        ("samples", &args.samples),
        ("samples_file", &args.samples_file),
    );
    let ploidy = match ploidy {
        Some(p) => p?,
        None => {
            if !samples_spec
                .as_ref()
                .is_some_and(|s| s.1 == ListSource::File)
            {
                eprintln!(
                    "Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid"
                );
            }
            Ploidy::from_str(ploidy::DEFAULT_PLOIDY, 2)?
        }
    };

    let consensus = args.consensus_caller;
    if consensus && args.multiallelic_caller {
        bail!("Only one of -c or -m options can be given");
    }
    if !consensus && !args.multiallelic_caller {
        bail!("Expected -c or -m option");
    }
    if consensus && gvcf.is_some() {
        bail!("gvcf -g option not functional with -c calling mode yet");
    }
    let targets_spec = last_of_pair(
        matches,
        ("targets", &args.targets),
        ("targets_file", &args.targets_file),
    );
    if shared.flags.constrain_alleles {
        if targets_spec.is_none() {
            bail!("Expected -t or -T with \"-C alleles\"");
        }
        if !args.multiallelic_caller {
            bail!("The \"-C alleles\" mode requires -m");
        }
    }
    if args.insert_missed && !shared.flags.constrain_alleles {
        bail!("The -i option requires -C alleles");
    }
    if shared.flags.var_only && gvcf.is_some() {
        bail!("The two options cannot be combined: --variants-only and --gvcf");
    }
    if args.group_samples.is_some() && !args.multiallelic_caller {
        bail!("The -G feature is supported only with the -m calling mode");
    }

    let config = Config {
        shared,
        acgt_only,
        indel_only,
        no_indel,
        insert_missed: args.insert_missed,
        consensus,
        threads,
    };
    let regions_spec = last_of_pair(
        matches,
        ("regions", &args.regions),
        ("regions_file", &args.regions_file),
    );
    init_and_run(
        &args,
        matches,
        config,
        &input,
        ploidy,
        gvcf,
        theta,
        samples_spec,
        targets_spec,
        regions_spec.map(|r| (r, regions_overlap)),
    )
}

/// bcftools' `init_data`, the main loop, and `destroy_data`.
#[expect(
    clippy::too_many_arguments,
    reason = "the state of bcftools call's main"
)]
fn init_and_run(
    args: &Args,
    matches: &ArgMatches,
    config: Config,
    input: &str,
    mut ploidy: Ploidy,
    mut gvcf: Option<Gvcf>,
    theta: f64,
    samples_spec: Option<ListSpec>,
    targets_spec: Option<ListSpec>,
    regions_spec: Option<(ListSpec, Overlap)>,
) -> anyhow::Result<ExitCode> {
    let shared = &config.shared;
    let (mut targets, stream_targets) = match &targets_spec {
        Some((spec, source)) => load_targets(spec, *source, shared.flags)?,
        None => (None, None),
    };
    let mut reader = open_reader(input, regions_spec.as_ref(), stream_targets, config.threads)?;

    // HTSlib leaves `errnum` set for a missing BGZF EOF marker, so bcftools call stops at
    // the first record.
    let eof = if reader.take_pending_error().is_some() {
        EofMarker::Missing
    } else {
        EofMarker::Present
    };

    let set = select_samples(
        reader.header(),
        &mut ploidy,
        samples_spec.as_ref(),
        shared.flags,
    )?;
    let nsex = ploidy.nsex();
    let pmax = ploidy.max();
    // Ploidies are bytes, as in bcftools, which takes them modulo 256.
    let mut sample_ploidy: Vec<u8> = vec![pmax as u8; set.nsamples];
    // -c sizes its tables by every sample's ploidy at the start, the largest one: bcftools
    // reads past them when it is 0 and crashes when it is above 2.
    if args.consensus_caller && !(1..=2).contains(&(pmax as u8)) {
        bail!(
            "Unsupported ploidy {pmax} in --ploidy-file: -c calls with the largest ploidy, which must be 1 or 2"
        );
    }

    if gvcf.is_some() {
        add_gvcf_header(reader.header_mut())?;
    }
    let mut hdr = output_header(reader.header(), &set)?;
    let (mut writer, out_type) = open_writer(args, matches, config.threads)?;
    let mut caller = init_caller(args, &mut hdr, &config, &set, &sample_ploidy, theta)?;
    hdr.remove(HeaderRecordType::Info, Some(b"QS"));
    hdr.remove(HeaderRecordType::Info, Some(b"I16"));
    if !args.no_version {
        version::append_version(&mut hdr, "call")?;
    }
    let out_name = output_name(args);
    write_header(args, &mut writer, &hdr, out_type, &out_name)?;

    if let Some(g) = gvcf.as_mut() {
        g.ranges.resolve(&hdr);
    }
    if let Some(t) = targets.as_mut() {
        t.out_name.clone_from(&out_name);
    }
    let mut state = MainLoop {
        hdr: &hdr,
        config: &config,
        shared,
        set: &set,
        ploidy: &ploidy,
        sex2ploidy: vec![0; nsex],
        sex2ploidy_prev: vec![pmax; nsex],
        sample_ploidy: &mut sample_ploidy,
        ploidy_checked: false,
        targets: targets.as_mut(),
        gvcf: gvcf.as_mut(),
        out_name: &out_name,
        eof,
    };
    // Without a thread to start, the records are called here, as bcftools goes on without its
    // output threads.
    let called = if config.threads > 0
        && parallel::supported(&state, &reader)
        && crate::output::threads_start(config.threads)
    {
        parallel::run(&mut state, &mut reader, &mut writer, &caller)
    } else {
        Ok(parallel::Called::NoThreads)
    };
    let result = match called {
        Ok(parallel::Called::Done) => Ok(()),
        Ok(parallel::Called::NoThreads) => state.run(&mut reader, &mut writer, &mut caller),
        Err(e) => Err(e),
    };
    // After a failed read the output is closed, keeping the calls before it as bcftools
    // does, and the read error is returned unless closing fails. An unsupported ploidy stops
    // the same way, so that the output ends at a record, with or without threads.
    let read_failure = match result {
        Ok(()) => None,
        Err(e) if e.is::<output::InputReadFailure>() => Some(e),
        Err(e) => {
            writer.abandon();
            return Err(e);
        }
    };
    match writer.finish() {
        Ok(()) => {}
        Err(WriteError::Index(_)) => bail!("Error: cannot write to index"),
        Err(WriteError::Deferred(e)) => {
            return Err(output::write_failure(e, || {
                anyhow::anyhow!("[main_vcfcall] Error: failed to write to {out_name}")
            }));
        }
        Err(WriteError::Io(e)) => {
            return Err(output::write_failure(e, || {
                anyhow::anyhow!("[destroy_data] Error: close failed .. {out_name}")
            }));
        }
        Err(e) => {
            return Err(anyhow::Error::new(e)
                .context(format!("[destroy_data] Error: close failed .. {out_name}")));
        }
    }
    if let Some(e) = read_failure {
        return Err(e);
    }
    Ok(ExitCode::SUCCESS)
}

/// Reads `-t`/`-T`: a file of CHROM POS [REF,ALT] targets, or regions the reader streams
/// (bcftools reads -t as a file name too; a -t value that is not a file is taken as a region
/// list here).
fn load_targets(
    spec: &str,
    source: ListSource,
    flags: CallFlags,
) -> anyhow::Result<(Option<Targets>, Option<StreamTargets>)> {
    if source == ListSource::File || Path::new(spec).exists() {
        let kind = if flags.constrain_alleles {
            TargetKind::Alleles
        } else {
            TargetKind::Positions
        };
        return Ok((Targets::load(spec, kind)?, None));
    }
    if flags.constrain_alleles {
        bail!("Expected -t or -T with \"-C alleles\"");
    }
    // `^REGIONS` excludes the regions, as in the other commands.
    let (list, logic) = match spec.strip_prefix('^') {
        Some(list) => (list, TargetsLogic::Exclude),
        None => (spec, TargetsLogic::Include),
    };
    let t = Regions::new(list.as_bytes(), ListSource::Inline)
        .map_err(|_| anyhow::anyhow!("Failed to read the targets: {spec}"))?;
    Ok((None, Some((t, logic))))
}

/// Opens the input with the `-r`/`-R` regions and the streamed `-t` regions.
fn open_reader(
    input: &str,
    regions_spec: Option<&(ListSpec, Overlap)>,
    stream_targets: Option<StreamTargets>,
    threads: usize,
) -> anyhow::Result<SyncedReader> {
    let regions = match regions_spec {
        Some(((spec, source), overlap)) => {
            let mut r = Regions::new(spec.as_bytes(), *source)
                .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
            r.overlap = *overlap;
            Some(r)
        }
        None => None,
    };
    let options = SyncedOptions {
        regions,
        targets_exclude: stream_targets
            .as_ref()
            .is_some_and(|t| t.1 == TargetsLogic::Exclude),
        targets: stream_targets.map(|t| t.0),
        apply_filters: None,
        threads,
    };
    match SyncedReader::open(input, options) {
        Ok(r) => Ok(r),
        Err(e) => {
            let name = if input == "-" {
                "standard input"
            } else {
                input
            };
            bail!("Failed to read from {name}: {e}");
        }
    }
}

/// Reads `-s`/`-S`, or takes all the samples of the input, and limits each sample's sex to
/// those the ploidy declares.
fn select_samples(
    hdr: &Header,
    ploidy: &mut Ploidy,
    samples_spec: Option<&ListSpec>,
    flags: CallFlags,
) -> anyhow::Result<SampleSet> {
    let mut set = SampleSet::default();
    if let Some((spec, source)) = samples_spec {
        set = samples::set_samples(hdr, ploidy, spec, *source)?;
        if flags.constrain_trio {
            if 3 * set.families.len() != set.nsamples {
                bail!("Expected only trios in {spec}, sorry!");
            }
            eprintln!(
                "Detected {} samples in {} trio families",
                set.nsamples,
                set.families.len()
            );
        }
    }
    let nsex = ploidy.nsex();
    if set.nsamples == 0 {
        set.nsamples = hdr.sample_count();
        set.sample2sex = vec![nsex as i32 - 1; set.nsamples];
    }
    for s in &mut set.sample2sex {
        if *s >= nsex as i32 {
            *s = nsex as i32 - 1;
        }
    }
    Ok(set)
}

/// Checks that the input defines FORMAT/DP, which `-g` needs, and adds the gVCF header lines.
fn add_gvcf_header(hdr: &mut Header) -> anyhow::Result<()> {
    let ok = hdr
        .id(b"DP")
        .is_some_and(|id| hdr.id_exists(HeaderRecordType::Format, id));
    if !ok {
        bail!(
            "--gvcf output mode requires FORMAT/DP tag, which is not present in the input header"
        );
    }
    Gvcf::update_header(hdr)
}

/// The input header, subset to the `-s`/`-S` samples when they are not all the samples in
/// order.
fn output_header(input: &Header, set: &SampleSet) -> anyhow::Result<Header> {
    if set.samples_map.is_some() {
        let names: Vec<&[u8]> = set.names.iter().map(Vec::as_slice).collect();
        let (h, _) = input
            .subset(&names)
            .context("Error occurred while subsetting samples")?;
        if h.sample_count() == 0 {
            bail!("No matching sample found");
        }
        Ok(h)
    } else {
        let h = input.dup().context("Could not copy the header")?;
        for n in &set.names {
            if h.sample_id(n).is_none() {
                bail!("No such sample: {}", String::from_utf8_lossy(n));
            }
        }
        Ok(h)
    }
}

/// The output name in messages: `-o`, or `-` for standard output.
fn output_name(args: &Args) -> String {
    args.output
        .as_deref()
        .map_or_else(|| "-".into(), |p| p.display().to_string())
}

/// Opens the output of `-o` and `-O`, returning the writer and its output type.
fn open_writer(
    args: &Args,
    matches: &ArgMatches,
    threads: usize,
) -> anyhow::Result<(Writer, OutputType)> {
    let (out_type, level) = output::resolve(
        &output_settings(args, matches),
        Some(args.output.as_deref().unwrap_or(Path::new("-"))),
    )?;
    let out_path = args.output.as_deref().map(output::data_path);
    match Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, threads),
    ) {
        Ok(w) => Ok((w, out_type)),
        Err(e) => {
            let name = output_name(args);
            let e = crate::raw_input::strerror(&e);
            bail!(
                "[E::hts_open_format] Failed to open file \"{name}\" : {e}\nError: cannot write to \"{name}\": {e}"
            );
        }
    }
}

/// Sets up the caller of `-c` or `-m` and adds its header lines (`ccall_init` and
/// `mcall_init`).
fn init_caller(
    args: &Args,
    hdr: &mut Header,
    config: &Config,
    set: &SampleSet,
    sample_ploidy: &[u8],
    theta: f64,
) -> anyhow::Result<Caller> {
    let shared = &config.shared;
    let nsmpl = hdr.sample_count();
    if config.consensus {
        append_ccall_header(hdr, shared)?;
        return Ok(Caller::Consensus(Box::new(CCall::new(
            hdr,
            nsmpl,
            sample_ploidy,
        ))));
    }
    let (groups, tag) = samples::sample_groups(
        hdr,
        args.group_samples.as_deref(),
        args.group_samples_tag.as_deref(),
    )?;
    if shared.flags.constrain_trio {
        samples::check_trios(hdr, shared, set, sample_ploidy)?;
    }
    append_mcall_header(hdr, shared)?;
    Ok(Caller::Multi(Box::new(MCall::new(
        hdr,
        groups,
        tag,
        theta,
        sample_ploidy,
        nsmpl,
    ))))
}

/// Writes the header and starts the `-W` index.
fn write_header(
    args: &Args,
    writer: &mut Writer,
    hdr: &Header,
    out_type: OutputType,
    out_name: &str,
) -> anyhow::Result<()> {
    if writer.write_header(hdr).is_err() {
        bail!("[init_data] Error: cannot write the header to {out_name}");
    }
    if let Some(wi) = args.write_index.last() {
        let fmt = output::parse_write_index(wi)?;
        let ok = match output::index_request(fmt, out_type, args.output.as_deref()) {
            Some(req) => writer.request_index(req).is_ok(),
            None => false,
        };
        if !ok {
            bail!("Error: failed to initialise index for {out_name}");
        }
    }
    writer.start_index(hdr)?;
    Ok(())
}

/// Collects `-O` in command-line order.
fn output_settings(args: &Args, matches: &ArgMatches) -> Vec<output::OutputSetting> {
    occurrences(matches, "output_type", &args.output_type)
        .into_iter()
        .map(|(_, v)| output::OutputSetting::Type(v.to_owned()))
        .collect()
}

/// The header lines of `mcall_init`.
fn append_mcall_header(hdr: &mut Header, shared: &Shared) -> anyhow::Result<()> {
    let mut lines: Vec<&[u8]> = Vec::new();
    if shared.flags.constrain_trio {
        lines.push(b"##FORMAT=<ID=CGT,Number=1,Type=Integer,Description=\"Constrained Genotype (0-based index to Number=G ordering).\">");
        lines.push(b"##FORMAT=<ID=UGT,Number=1,Type=Integer,Description=\"Unconstrained Genotype (0-based index to Number=G ordering).\">");
    }
    lines.push(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if shared.tags.gq {
        lines.push(
            b"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled Genotype Quality\">",
        );
    }
    if shared.tags.gp {
        lines.push(b"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities in the range 0 to 1\">");
    }
    lines.push(b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">");
    lines.push(b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    lines.push(b"##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">");
    lines.push(b"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">");
    if shared.tags.pv4 {
        lines.push(b"##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">");
    }
    for l in lines {
        hdr.append(l)?;
    }
    Ok(())
}

/// The header lines of `ccall_init`.
fn append_ccall_header(hdr: &mut Header, shared: &Shared) -> anyhow::Result<()> {
    hdr.append(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    if shared.tags.gq {
        hdr.append(b"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">")?;
    }
    if shared.tags.gp {
        bail!("Sorry, -f GP is not supported with -c");
    }
    for l in [
        b"##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele frequency (assuming HWE)\">".as_slice(),
        b"##INFO=<ID=AF2,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)\">",
        b"##INFO=<ID=AC1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele count (no HWE assumption)\">",
        b"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">",
        b"##INFO=<ID=FQ,Number=1,Type=Float,Description=\"Phred probability of all samples being the same\">",
        b"##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">",
        b"##INFO=<ID=G3,Number=3,Type=Float,Description=\"ML estimate of genotype frequencies\">",
        b"##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Chi^2 based HWE test P-value based on G3\">",
        b"##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">",
    ] {
        hdr.append(l)?;
    }
    Ok(())
}

/// What a record became before calling.
pub enum Prepared {
    /// Skipped (filtered out).
    Skip,
    /// To be called, with the unseen allele's index (0 for none).
    Call { unseen: usize },
}

/// The main loop's state (the parts of bcftools' `args_t` it uses).
pub struct MainLoop<'a> {
    hdr: &'a Header,
    config: &'a Config,
    shared: &'a Shared,
    set: &'a SampleSet,
    ploidy: &'a Ploidy,
    sex2ploidy: Vec<i32>,
    sex2ploidy_prev: Vec<i32>,
    sample_ploidy: &'a mut Vec<u8>,
    /// Whether [`check_ploidy`] has passed `sample_ploidy` as it stands.
    ploidy_checked: bool,
    targets: Option<&'a mut Targets>,
    gvcf: Option<&'a mut Gvcf>,
    out_name: &'a str,
    eof: EofMarker,
}

/// Fails when a sample has a ploidy that neither caller genotypes, at the first record where
/// it has it. bcftools -m counts it in AN at reference sites and aborts at the first variant
/// site (where assertions are compiled in), and -c calls nonsense.
fn check_ploidy(hdr: &Header, chr: &[u8], pos: i64, sample_ploidy: &[u8]) -> anyhow::Result<()> {
    if let Some((i, p)) = sample_ploidy.iter().enumerate().find(|(_, p)| **p > 2) {
        let sample = hdr.samples().get(i).map(|s| String::from_utf8_lossy(s));
        // The calls before the record are kept, as after a failed read.
        return Err(output::InputReadFailure(format!(
            "Unsupported ploidy {p} of sample {} at {}:{}: only ploidy 0, 1 and 2 can be called",
            sample.unwrap_or_default(),
            String::from_utf8_lossy(chr),
            pos + 1
        ))
        .into());
    }
    Ok(())
}

/// The checks of `main_vcfcall`'s loop before calling: variant types, REF=N, the unseen
/// allele, and `-v` at reference sites.
fn prepare(config: &Config, shared: &Shared, rec: &Record) -> Prepared {
    // One pass over the alleles: whether the record is a SNP, and its first unseen allele.
    let mut alleles = rec.peek_alleles();
    let ref_allele = alleles.next().unwrap_or_default();
    let mut snp = rec.allele_count() == 0 || is_snp_allele(ref_allele);
    let mut unseen = 0;
    for (i, a) in (1..).zip(alleles) {
        snp &= is_snp_allele(a);
        if unseen == 0
            && (a.first() == Some(&b'X')
                || (a.first() == Some(&b'<')
                    && (a.get(1) == Some(&b'X') || a.get(1) == Some(&b'*'))
                    && a.get(2) == Some(&b'>')))
        {
            unseen = i;
        }
    }
    let is_indel = !snp;
    if config.indel_only && !is_indel {
        return Prepared::Skip;
    }
    if config.no_indel && is_indel {
        return Prepared::Skip;
    }
    if config.acgt_only && matches!(ref_allele.first(), Some(b'N' | b'n')) {
        return Prepared::Skip;
    }
    let n = rec.allele_count();
    let is_ref = n == 1 || (n == 2 && unseen > 0);
    if is_ref && shared.flags.var_only {
        return Prepared::Skip;
    }
    Prepared::Call { unseen }
}

impl MainLoop<'_> {
    /// bcftools' `set_ploidy`, returning whether the samples' ploidies changed.
    pub fn set_ploidy(&mut self, rec: &Record) -> bool {
        let chr = rec.chrom_safe(self.hdr);
        self.ploidy
            .query(chr, rec.pos() as u32, &mut self.sex2ploidy);
        if self.sex2ploidy == self.sex2ploidy_prev {
            return false;
        }
        for (p, &sex) in self.sample_ploidy.iter_mut().zip(&self.set.sample2sex) {
            *p = if sex < 0 {
                (-sex) as u8
            } else {
                self.sex2ploidy[sex as usize] as u8
            };
        }
        std::mem::swap(&mut self.sex2ploidy, &mut self.sex2ploidy_prev);
        true
    }

    /// Reads the next record to call (bcftools' `next_line`), with its `-C alleles`
    /// target alleles.
    fn next_line(
        &mut self,
        reader: &mut SyncedReader,
        rec: &mut Box<Record>,
    ) -> anyhow::Result<Option<Option<Vec<Vec<u8>>>>> {
        let subset = self.set.samples_map.as_deref();
        match self.targets.as_deref_mut() {
            Some(t) if t.has_alleles() => t.next_with_alleles(reader, subset, rec, self.eof),
            t => loop {
                if !read_record(reader, rec, self.eof)? {
                    return Ok(None);
                }
                if let Some(t) = &t
                    && !t.matches(reader.header(), rec)
                {
                    continue;
                }
                if let Some(map) = subset {
                    rec.subset(map);
                }
                return Ok(Some(None));
            },
        }
    }

    /// The serial main loop of `main_vcfcall`.
    fn run(
        &mut self,
        reader: &mut SyncedReader,
        writer: &mut Writer,
        caller: &mut Caller,
    ) -> anyhow::Result<()> {
        let mut rec = Box::new(Record::new());
        let mut dp = Vec::new();
        let mut spare = None;
        let out_name = self.out_name.to_owned();
        let hdr = self.hdr;
        while let Some(tgt_als) = self.next_line(reader, &mut rec)? {
            if self.shared.flags.constrain_alleles && tgt_als.is_none() {
                continue;
            }
            let Prepared::Call { unseen } = prepare(self.config, self.shared, &rec) else {
                continue;
            };
            rec.unpack(Unpack::ALL);
            if self.ploidy.nsex() > 0 && self.set_ploidy(&rec) {
                self.ploidy_checked = false;
            }
            if !self.ploidy_checked {
                check_ploidy(hdr, rec.chrom_safe(hdr), rec.pos(), self.sample_ploidy)?;
                self.ploidy_checked = true;
            }
            if self.config.insert_missed
                && let Some(t) = self.targets.as_deref_mut()
            {
                t.flush(hdr, Some(&rec), &mut |line| {
                    write_one(writer, hdr, line, &out_name)
                })?;
            }
            let ret = caller.call(
                hdr,
                &mut rec,
                self.shared,
                self.sample_ploidy,
                unseen,
                tgt_als.as_deref(),
            )?;
            if self.shared.flags.var_only && ret == 0 && self.gvcf.is_none() {
                continue;
            }
            if let Some(g) = self.gvcf.as_deref_mut() {
                let item = g.ranges.classify(
                    hdr,
                    &mut rec,
                    CalledSite::from_call(ret),
                    &mut dp,
                    &mut spare,
                )?;
                let write_rec = matches!(item, GvcfItem::Write { .. });
                let freed = g.advance(hdr, Some(item), |line| {
                    write_one(writer, hdr, line, &out_name)
                })?;
                if freed.is_some() {
                    spare = freed;
                }
                if !write_rec {
                    continue;
                }
            }
            write_one(writer, hdr, &mut rec, &out_name)?;
        }
        self.finish(writer)
    }

    /// Flushes the gVCF block and the targets not seen (end of `main_vcfcall`).
    fn finish(&mut self, writer: &mut Writer) -> anyhow::Result<()> {
        let hdr = self.hdr;
        let out_name = self.out_name.to_owned();
        if let Some(g) = self.gvcf.as_deref_mut() {
            g.advance(hdr, None, |line| write_one(writer, hdr, line, &out_name))?;
        }
        if self.config.insert_missed
            && let Some(t) = self.targets.as_deref_mut()
        {
            t.flush(hdr, None, &mut |line| {
                write_one(writer, hdr, line, &out_name)
            })?;
        }
        Ok(())
    }
}

/// Whether the input has its BGZF EOF marker. HTSlib reads past a missing one but leaves
/// `errnum` set, which `bcftools call` reports at the first record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum EofMarker {
    /// Present, or not a BGZF input.
    Present,
    /// Missing.
    Missing,
}

/// Reads a record from the synced reader. An unparsable VCF line ends the input quietly, as
/// `bcf_sr_next_line` does for `bcftools call` (HTSlib reports the line), but records with
/// errors (or any record of an input without its EOF marker) are fatal. So is a failed read (a
/// truncated or corrupt BGZF block, a malformed BCF record), which bcftools also takes for the
/// end of the input, dropping the records after it and the one read ahead, with status 0.
fn read_record(
    reader: &mut SyncedReader,
    rec: &mut Box<Record>,
    eof: EofMarker,
) -> anyhow::Result<bool> {
    match reader.next_record_boxed(rec) {
        Ok(true) => {}
        Ok(false) | Err(SyncedError::VcfParse) => return Ok(false),
        Err(e) => return Err(output::InputReadFailure(format!("Error: {e}")).into()),
    }
    if eof == EofMarker::Missing || !rec.errors().is_empty() {
        bail!("Error: could not parse the input VCF");
    }
    Ok(true)
}

/// Writes one record (`bcf_write1`).
fn write_one(
    writer: &mut Writer,
    hdr: &Header,
    rec: &mut Record,
    out_name: &str,
) -> anyhow::Result<()> {
    match writer.write_record(hdr, rec) {
        Ok(()) => Ok(()),
        Err(WriteError::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => {
            Err(anyhow::Error::new(e))
        }
        Err(_) => bail!("[main_vcfcall] Error: failed to write to {out_name}"),
    }
}
