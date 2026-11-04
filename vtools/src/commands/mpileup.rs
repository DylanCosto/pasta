//! `vtools mpileup`: genotype likelihoods from alignments (bcftools mpileup).
//!
//! A port of bcftools 1.24's `mpileup.c`, `bam2bcf.c`, `bam2bcf_indel.c`, `bam2bcf_iaux.c`
//! (`--indels-2.0`), `bam2bcf_edlib.c` and `edlib.c` (`--indels-cns`), `read_consensus.c`,
//! `bam_sample.c`, and `gvcf.c` on top of pasta's HTSlib-compatible pileup engine and BAQ.
//!
//! With `--threads`, regions (`-r/-R`), or whole chromosomes of indexed inputs, are piled up
//! on worker threads and written in order. The random numbers HTSlib draws to subsample deep
//! pileups are handed to each region in sequence, so the output is that of a serial run.
//! The workers and the realignments of heavy indel positions share one budget of 2N threads.

mod budget;
mod call;
mod consensus;
mod edlib;
mod errmod;
mod gvcf;
mod iaux;
mod indel;
mod input;
mod parallel;
mod pileup;
mod reference;
mod regidx;
mod sample;

use std::{
    io::Write as _,
    path::{Path, PathBuf},
    process::ExitCode,
    thread,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, Record,
    io::{OutputType, WriteError, Writer, reader::is_stream},
    readlist::ListSource,
};

use self::{
    budget::ThreadBudget,
    call::{AmbigReads, CallSettings},
    errmod::Rng48,
    gvcf::Gvcf,
    input::AlnReader,
    parallel::{Encoding, Piled},
    pileup::{Conf, Engine, InputFile, RngSource, Shared, Sink, Span},
    reference::{Fasta, References},
    regidx::{LineFormat, REGIDX_MAX, RegIdx},
    sample::Samples,
};
use crate::output::{self, OutputSetting};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input SAM/BAM/CRAM files.
    inputs: Vec<String>,

    /// Quality is in the Illumina-1.3+ encoding.
    #[arg(short = '6', long = "illumina1.3+", action = clap::ArgAction::Count)]
    illumina13: u8,
    /// Include anomalous read pairs, with flag PAIRED but not PROPER_PAIR set.
    #[arg(short = 'A', long = "count-orphans", action = clap::ArgAction::Count)]
    count_orphans: u8,
    /// List of input BAM filenames, one per line.
    #[arg(short = 'b', long = "bam-list", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    bam_list: Vec<String>,
    /// Disable BAQ (per-Base Alignment Quality).
    #[arg(short = 'B', long = "no-BAQ", visible_alias = "no-baq", action = clap::ArgAction::Count)]
    no_baq: u8,
    /// Adjust mapping quality [0].
    #[arg(short = 'C', long = "adjust-MQ", visible_alias = "adjust-mq", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    adjust_mq: Vec<String>,
    /// Apply BAQ everywhere, not just in problematic regions.
    #[arg(short = 'D', long = "full-BAQ", visible_alias = "full-baq", action = clap::ArgAction::Count)]
    full_baq: u8,
    /// Max raw per-file depth; avoids excessive memory usage [250].
    #[arg(short = 'd', long = "max-depth", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    max_depth: Vec<String>,
    /// Recalculate BAQ on the fly, ignore existing BQs.
    #[arg(short = 'E', long = "redo-BAQ", visible_alias = "redo-baq", action = clap::ArgAction::Count)]
    redo_baq: u8,
    /// Faidx indexed reference sequence file.
    #[arg(short = 'f', long = "fasta-ref", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    fasta_ref: Vec<String>,
    /// Do not require fasta reference file.
    #[arg(long = "no-reference", action = clap::ArgAction::Count)]
    no_reference: u8,
    /// Select or exclude read groups listed in the file.
    #[arg(short = 'G', long = "read-groups", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    read_groups: Vec<String>,
    /// Skip alignments with mapQ smaller than INT [0].
    #[arg(short = 'q', long = "min-MQ", visible_alias = "min-mq", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    min_mq: Vec<String>,
    /// Skip bases with baseQ/BAQ smaller than INT [1].
    #[arg(short = 'Q', long = "min-BQ", visible_alias = "min-bq", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    min_bq: Vec<String>,
    /// Limit baseQ/BAQ to no more than INT [60].
    #[arg(long = "max-BQ", visible_alias = "max-bq", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    max_bq: Vec<String>,
    /// Use neighbour_qual + INT if less than qual [30].
    #[arg(long = "delta-BQ", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    delta_bq: Vec<String>,
    /// Comma separated list of regions in which pileup is generated.
    #[arg(short = 'r', long = "regions", visible_alias = "region", value_name = "REG", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions: Vec<String>,
    /// Restrict to regions listed in a file.
    #[arg(short = 'R', long = "regions-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    regions_file: Vec<String>,
    /// Ignore RG tags (one BAM = one sample).
    #[arg(long = "ignore-RG", visible_alias = "ignore-rg", action = clap::ArgAction::Count)]
    ignore_rg: u8,
    /// Skip reads with all of the bits set.
    #[arg(long = "skip-all-set", visible_alias = "ls", value_name = "STR|INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    skip_all_set: Vec<String>,
    /// Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP].
    #[arg(long = "skip-any-set", visible_aliases = ["ns", "ff"], value_name = "STR|INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    skip_any_set: Vec<String>,
    /// Skip reads with all of the bits unset.
    #[arg(long = "skip-all-unset", visible_aliases = ["lu", "rf"], value_name = "STR|INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    skip_all_unset: Vec<String>,
    /// Skip reads with any of the bits unset.
    #[arg(long = "skip-any-unset", visible_alias = "nu", value_name = "STR|INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    skip_any_unset: Vec<String>,
    /// Comma separated list of samples to include.
    #[arg(short = 's', long = "samples", value_name = "LIST", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    samples: Vec<String>,
    /// File of samples to include.
    #[arg(short = 'S', long = "samples-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    samples_file: Vec<String>,
    /// Similar to -r but streams rather than index-jumps.
    #[arg(short = 't', long = "targets", value_name = "REG", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    targets: Vec<String>,
    /// Similar to -R but streams rather than index-jumps.
    #[arg(short = 'T', long = "targets-file", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    targets_file: Vec<String>,
    /// Disable read-pair overlap detection.
    #[arg(short = 'x', long = "ignore-overlaps", action = clap::ArgAction::Count)]
    ignore_overlaps: u8,
    /// Random number seed used for sampling deep regions [0].
    #[arg(long = "seed", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    seed: Vec<String>,
    /// Optional tags to output; '?' to list available tags.
    #[arg(short = 'a', long = "annotate", value_name = "LIST", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    annotate: Vec<String>,
    /// Group non-variant sites into gVCF blocks by minimum per-sample DP.
    #[arg(short = 'g', long = "gvcf", value_name = "INT[,...]", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    gvcf: Vec<String>,
    /// Do not append version and command line to the header.
    #[arg(long = "no-version", action = clap::ArgAction::Count)]
    no_version: u8,
    /// Write output to FILE [standard output].
    #[arg(long = "output", value_name = "FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output: Vec<String>,
    /// Phred-scaled gap open seq error probability [40], or with a file name, the output.
    #[arg(short = 'o', value_name = "INT|FILE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    o_short: Vec<String>,
    /// Phred-scaled gap open seq error probability [40].
    #[arg(long = "open-prob", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    open_prob: Vec<String>,
    /// 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF;
    /// 0-9 compression level [v].
    #[arg(short = 'O', long = "output-type", value_name = "TYPE", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    output_type: Vec<String>,
    /// Use multithreading with INT worker threads [0].
    #[arg(long = "threads", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    threads: Vec<String>,
    /// Verbosity level.
    #[arg(short = 'v', long = "verbosity", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    verbosity: Vec<String>,
    /// Automatically index the output files [off].
    #[arg(short = 'W', long = "write-index", value_name = "FMT", num_args = 0..=1, require_equals = true, default_missing_value = "", action = clap::ArgAction::Append)]
    write_index: Vec<String>,
    /// Specify platform profile (use "-X list" for details).
    #[arg(short = 'X', long = "config", value_name = "STR", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    config: Vec<String>,
    /// Phred-scaled gap extension seq error probability [20].
    #[arg(short = 'e', long = "ext-prob", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    ext_prob: Vec<String>,
    /// Minimum fraction of gapped reads [0.05].
    #[arg(short = 'F', long = "gap-frac", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    gap_frac: Vec<String>,
    /// Coefficient for homopolymer errors [500].
    #[arg(short = 'h', long = "tandem-qual", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    tandem_qual: Vec<String>,
    /// Do not perform indel calling.
    #[arg(short = 'I', long = "skip-indels", action = clap::ArgAction::Count)]
    skip_indels: u8,
    /// Maximum per-file depth for INDEL calling [250].
    #[arg(short = 'L', long = "max-idepth", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    max_idepth: Vec<String>,
    /// Minimum number gapped reads for indel candidates [2].
    #[arg(short = 'm', long = "min-ireads", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    min_ireads: Vec<String>,
    /// Maximum length of read to pass to BAQ algorithm [500].
    #[arg(short = 'M', long = "max-read-len", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    max_read_len: Vec<String>,
    /// Apply -m and -F per-sample for increased sensitivity.
    #[arg(short = 'p', long = "per-sample-mF", visible_alias = "per-sample-mf", action = clap::ArgAction::Count)]
    per_sample: u8,
    /// Comma separated list of platforms for indels [all] (accepted and ignored, as in
    /// bcftools 1.24).
    #[arg(short = 'P', long = "platforms", value_name = "STR", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    platforms: Vec<String>,
    /// What to do with ambiguous indel reads: drop,incAD,incAD0 [drop].
    #[arg(long = "ambig-reads", visible_alias = "ar", value_name = "STR", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    ambig_reads: Vec<String>,
    /// Raise to favour recall over precision [1.00].
    #[arg(long = "indel-bias", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    indel_bias: Vec<String>,
    /// Relative likelihood of insertion to deletion [0.00].
    #[arg(long = "del-bias", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    del_bias: Vec<String>,
    /// Ratio of score vs ref (1) or 2nd-best allele (0) [0.00].
    #[arg(long = "score-vs-ref", value_name = "FLOAT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    score_vs_ref: Vec<String>,
    /// Approximate maximum indel size considered [110].
    #[arg(long = "indel-size", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    indel_size: Vec<String>,
    /// New EXPERIMENTAL indel calling model (diploid reference consensus).
    #[arg(long = "indels-2.0", action = clap::ArgAction::Count)]
    indels_20: u8,
    /// New EXPERIMENTAL indel calling model with edlib.
    #[arg(long = "indels-cns", action = clap::ArgAction::Count)]
    indels_cns: u8,
    /// Indel-cns tuning for indel seq-qual scores [120].
    #[arg(long = "seqq-offset", value_name = "INT", action = clap::ArgAction::Append, allow_hyphen_values = true)]
    seqq_offset: Vec<String>,
    /// Disable CNS mode, to use after a -X profile.
    #[arg(long = "no-indels-cns", action = clap::ArgAction::Count)]
    no_indels_cns: u8,
    /// (Edlib mode) Use minimum quality within homopolymers.
    #[arg(long = "poly-mqual", action = clap::ArgAction::Count)]
    poly_mqual: u8,
    /// Disable --poly-mqual.
    #[arg(long = "no-poly-mqual", action = clap::ArgAction::Count)]
    no_poly_mqual: u8,
    /// (Not an option of bcftools 1.24; rejected as bcftools rejects it.)
    #[arg(short = 'U', hide = true, action = clap::ArgAction::Count)]
    removed_u: u8,
}

/// An error bcftools reports and exits from with status 1 (`goto err`, `exit(EXIT_FAILURE)`).
#[derive(Debug)]
struct Exit1(String);

impl std::fmt::Display for Exit1 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Exit1 {}

fn exit1(msg: impl Into<String>) -> anyhow::Error {
    anyhow::Error::new(Exit1(msg.into()))
}

/// C's `atoi`.
fn atoi(s: &str) -> i32 {
    pasta::variant::number::atoi(s.as_bytes()) as i32
}

/// C's `atof`.
fn atof(s: &str) -> f64 {
    pasta::variant::number::atof(s.as_bytes())
}

/// The options, as bcftools' option loop leaves them (`mplp_conf_t` and locals).
struct Options {
    min_mq: i32,
    min_baseq: i32,
    max_baseq: i32,
    delta_baseq: i32,
    capq_thres: i32,
    max_depth: i32,
    max_indel_depth: i32,
    openq: i32,
    extq: i32,
    tandemq: i32,
    min_frac: f64,
    indel_bias: f64,
    min_support: i32,
    vs_ref: f64,
    flag: u32,
    skip_any_unset: i32,
    skip_all_unset: i32,
    skip_any_set: i32,
    skip_all_set: i32,
    output_fname: Option<String>,
    output_settings: Vec<OutputSetting>,
    record_cmd_line: bool,
    n_threads: usize,
    fmt_flag: u32,
    max_read_len: i32,
    ambig_reads: AmbigReads,
    indel_win_size: i32,
    poly_mqual: bool,
    seqq_offset: i32,
    del_bias: f64,
    seed: i64,
    gvcf: Option<Gvcf>,
    fai_fname: Option<String>,
    /// The `-f` reference, loaded as bcftools loads it in the option loop.
    fasta: Option<Fasta>,
    noref: bool,
    /// `-r/-R`: the regions, a list or a file.
    reg: Option<(String, ListSource)>,
    bed: Option<RegIdx>,
    bed_logic: bool,
    file_list: Option<String>,
    write_index: Option<output::WriteIndex>,
    use_orphan: bool,
    indels_v20: bool,
    edlib: bool,
    samples: Samples,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            min_mq: 0,
            min_baseq: 1,
            max_baseq: 60,
            delta_baseq: 30,
            capq_thres: 0,
            max_depth: 250,
            max_indel_depth: 250,
            openq: 40,
            extq: 20,
            tandemq: 500,
            min_frac: 0.05,
            indel_bias: 1.0,
            min_support: 2,
            vs_ref: 0.0,
            flag: pileup::NO_ORPHAN
                | pileup::REALN
                | pileup::REALN_PARTIAL
                | pileup::SMART_OVERLAPS,
            skip_any_unset: 0,
            skip_all_unset: 0,
            skip_any_set: 0x4 | 0x100 | 0x200 | 0x400,
            skip_all_set: 0,
            output_fname: None,
            output_settings: Vec::new(),
            record_cmd_line: true,
            n_threads: 0,
            fmt_flag: call::INFO_BQBZ
                | call::INFO_IDV
                | call::INFO_IMF
                | call::INFO_MQ0F
                | call::INFO_MQBZ
                | call::INFO_MQSBZ
                | call::INFO_RPBZ
                | call::INFO_SCBZ
                | call::INFO_SGB
                | call::INFO_VDB
                | call::FMT_AD,
            max_read_len: 500,
            ambig_reads: AmbigReads::Drop,
            indel_win_size: 110,
            poly_mqual: false,
            seqq_offset: 120,
            del_bias: 0.0,
            seed: 0,
            gvcf: None,
            fai_fname: None,
            fasta: None,
            noref: false,
            reg: None,
            bed: None,
            bed_logic: false,
            file_list: None,
            write_index: None,
            use_orphan: false,
            indels_v20: false,
            edlib: false,
            samples: Samples::default(),
        }
    }
}

/// HTSlib's `bam_str2flag`: a number (C `strtol` with base 0, so `0x` and leading-`0` octal,
/// and trailing text ignored), or comma-separated flag names in any case; -1 on error.
fn str2flag(s: &str) -> i32 {
    if let Some(n) = strtol_base0(s) {
        // `long` to `int`, as the C assignment truncates it.
        return n as i32;
    }
    let mut flag = 0;
    for name in s.split(',') {
        flag |= match name.to_ascii_uppercase().as_str() {
            "PAIRED" => 0x1,
            "PROPER_PAIR" => 0x2,
            "UNMAP" => 0x4,
            "MUNMAP" => 0x8,
            "REVERSE" => 0x10,
            "MREVERSE" => 0x20,
            "READ1" => 0x40,
            "READ2" => 0x80,
            "SECONDARY" => 0x100,
            "QCFAIL" => 0x200,
            "DUP" => 0x400,
            "SUPPLEMENTARY" => 0x800,
            // The empty string has no names; an empty name between commas is an error.
            "" if s.is_empty() => 0,
            _ => return -1,
        };
    }
    flag
}

/// C's `strtol(s, &end, 0)`, or `None` when it converts nothing (`end == s`).
fn strtol_base0(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    let mut i = b
        .iter()
        .take_while(|c| matches!(c, b' ' | b'\t'..=b'\r'))
        .count();
    let neg = match b.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let hex = b.get(i) == Some(&b'0')
        && matches!(b.get(i + 1), Some(b'x' | b'X'))
        && b.get(i + 2).is_some_and(u8::is_ascii_hexdigit);
    let (radix, start) = if hex {
        (16, i + 2)
    } else if b.get(i) == Some(&b'0') {
        (8, i)
    } else {
        (10, i)
    };
    let digits: Vec<u32> = b[start.min(b.len())..]
        .iter()
        .map_while(|&c| char::from(c).to_digit(radix))
        .collect();
    if digits.is_empty() {
        return None;
    }
    let magnitude = digits
        .iter()
        .try_fold(0i128, |acc, &d| {
            let v = acc * i128::from(radix) + i128::from(d);
            (v <= i128::from(i64::MAX) + 1).then_some(v)
        })
        .unwrap_or(i128::from(i64::MAX) + 1);
    let value = if neg { -magnitude } else { magnitude };
    Some(value.clamp(i128::from(i64::MIN), i128::from(i64::MAX)) as i64)
}

/// Checks `-O`'s value as bcftools mpileup's option loop does (it has its own messages).
fn check_output_type(v: &str) -> anyhow::Result<()> {
    let level = |s: &str| pasta::variant::number::strtol(s.as_bytes(), 0);
    match v.as_bytes().first() {
        // bcftools reads the level from past the end of an empty value (README).
        None => bail!("The output type \"\" not recognised"),
        Some(b'b' | b'u' | b'z' | b'v') => {}
        _ => {
            let (l, end) = level(v);
            if end != v.len() || !(0..=9).contains(&l) {
                bail!("The output type \"{v}\" not recognised");
            }
        }
    }
    if let Some(rest) = v.get(1..).filter(|r| !r.is_empty()) {
        let (l, end) = level(rest);
        if end != rest.len() || !(0..=9).contains(&l) {
            bail!("Could not parse argument: --output-type {rest}");
        }
    }
    Ok(())
}

/// The FORMAT tags of `-a`: name, flag, and the deprecation warning printed when it is set.
const FORMAT_TAGS: [(&str, u32, &str); 12] = [
    ("AD", call::FMT_AD, ""),
    ("ADF", call::FMT_ADF, ""),
    ("ADR", call::FMT_ADR, ""),
    ("DP", call::FMT_DP, ""),
    (
        "DP4",
        call::FMT_DP4,
        "[warning] tag DP4 functional, but deprecated. Please switch to `ADF` and `ADR` in future.\n",
    ),
    (
        "DPR",
        call::FMT_DPR,
        "[warning] tag DPR functional, but deprecated. Please switch to `AD` in future.\n",
    ),
    (
        "DV",
        call::FMT_DV,
        "[warning] tag DV functional, but deprecated. Please switch to `AD` in future.\n",
    ),
    ("NMBZ", call::FMT_NMBZ, ""),
    ("QM", call::FMT_QM, ""),
    ("QS", call::FMT_QS, ""),
    ("SP", call::FMT_SP, ""),
    ("SCR", call::FMT_SCR, ""),
];

/// The INFO tags of `-a`. bcftools' `SET_INFO_FLAG` prints no warning for the deprecated
/// `INFO/DPR`.
const INFO_TAGS: [(&str, u32); 17] = [
    ("DPR", call::INFO_DPR),
    ("AD", call::INFO_AD),
    ("ADF", call::INFO_ADF),
    ("ADR", call::INFO_ADR),
    ("BQBZ", call::INFO_BQBZ),
    ("FS", call::INFO_FS),
    ("IDV", call::INFO_IDV),
    ("IMF", call::INFO_IMF),
    ("MQ0F", call::INFO_MQ0F),
    ("MQBZ", call::INFO_MQBZ),
    ("NM", call::INFO_NM),
    ("NMBZ", call::INFO_NMBZ),
    ("RPBZ", call::INFO_RPBZ),
    ("SCBZ", call::INFO_SCBZ),
    ("SCR", call::INFO_SCR),
    ("SGB", call::INFO_SGB),
    ("VDB", call::INFO_VDB),
];

/// `parse_format_flag` (`-a`).
fn parse_format_flag(flag: &mut u32, s: &str) -> anyhow::Result<()> {
    let tags = pasta::variant::readlist::read_list(s.as_bytes(), ListSource::Inline)?;
    for raw in tags {
        let t = String::from_utf8_lossy(&raw).into_owned();
        let (exclude, tag) = match t.strip_prefix('-') {
            Some(rest) => (true, rest),
            None => (false, t.as_str()),
        };
        let set = |flag: &mut u32, bit: u32| {
            if exclude {
                *flag &= !bit;
            } else {
                *flag |= bit;
            }
        };
        let format_tag = FORMAT_TAGS.iter().find(|(name, ..)| {
            tag.eq_ignore_ascii_case(name)
                || tag.eq_ignore_ascii_case(&format!("FMT/{name}"))
                || tag.eq_ignore_ascii_case(&format!("FORMAT/{name}"))
        });
        if let Some(&(_, bit, msg)) = format_tag {
            if !msg.is_empty() {
                eprint!("{msg}");
            }
            set(flag, bit);
            continue;
        }
        if tag.len() > 5 && tag[..5].eq_ignore_ascii_case("INFO/") {
            let rest = &tag[5..];
            if let Some(&(_, bit)) = INFO_TAGS.iter().find(|(n, _)| rest.eq_ignore_ascii_case(n)) {
                set(flag, bit);
                continue;
            }
        }
        return Err(exit1(format!("Could not parse tag \"{tag}\" in \"{s}\"")));
    }
    Ok(())
}

/// `print_usage`, with the default values.
const USAGE: &str = r#"
Usage: vtools mpileup [options] in1.bam [in2.bam [...]]

Input options:
  -6, --illumina1.3+      Quality is in the Illumina-1.3+ encoding
  -A, --count-orphans     Include anomalous read pairs, with flag PAIRED but not PROPER_PAIR set
  -b, --bam-list FILE     List of input BAM filenames, one per line
  -B, --no-BAQ            Disable BAQ (per-Base Alignment Quality)
  -C, --adjust-MQ INT     Adjust mapping quality [0]
  -D, --full-BAQ          Apply BAQ everywhere, not just in problematic regions
  -d, --max-depth INT     Max raw per-file depth; avoids excessive memory usage [250]
  -E, --redo-BAQ          Recalculate BAQ on the fly, ignore existing BQs
  -f, --fasta-ref FILE    Faidx indexed reference sequence file
      --no-reference      Do not require fasta reference file
  -G, --read-groups FILE  Select or exclude read groups listed in the file
  -q, --min-MQ INT        Skip alignments with mapQ smaller than INT [0]
  -Q, --min-BQ INT        Skip bases with baseQ/BAQ smaller than INT [1]
      --max-BQ INT        Limit baseQ/BAQ to no more than INT [60]
      --delta-BQ INT      Use neighbour_qual + INT if less than qual [30]
  -r, --regions REG[,...] Comma separated list of regions in which pileup is generated
  -R, --regions-file FILE Restrict to regions listed in a file
      --ignore-RG         Ignore RG tags (one BAM = one sample)
  --ls, --skip-all-set STR|INT  Skip reads with all of the bits set []
  --ns, --skip-any-set STR|INT  Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP]
  --lu, --skip-all-unset STR|INT  Skip reads with all of the bits unset []
  --nu, --skip-any-unset STR|INT  Skip reads with any of the bits unset []
  -s, --samples LIST      Comma separated list of samples to include
  -S, --samples-file FILE File of samples to include
  -t, --targets REG[,...] Similar to -r but streams rather than index-jumps
  -T, --targets-file FILE Similar to -R but streams rather than index-jumps
  -x, --ignore-overlaps   Disable read-pair overlap detection
      --seed INT          Random number seed used for sampling deep regions [0]

Output options:
  -a, --annotate LIST     Optional tags to output; '\?' to list available tags []
  -g, --gvcf INT[,...]    Group non-variant sites into gVCF blocks according
                          To minimum per-sample DP
      --no-version        Do not append version and command line to the header
  -o, --output FILE       Write output to FILE [standard output]
  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;
                          'z' compressed VCF; 'v' uncompressed VCF; 0-9 compression level [v]
      --threads INT       Use multithreading with INT worker threads [0]
  -v, --verbosity INT     Verbosity level
  -W, --write-index[=FMT] Automatically index the output files [off]

SNP/INDEL genotype likelihoods options:
  -X, --config STR        Specify platform profile (use "-X list" for details)
  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [20]
  -F, --gap-frac FLOAT    Minimum fraction of gapped reads [0.05]
  -h, --tandem-qual INT   Coefficient for homopolymer errors [500]
  -I, --skip-indels       Do not perform indel calling
  -L, --max-idepth INT    Maximum per-file depth for INDEL calling [250]
  -m, --min-ireads INT    Minimum number gapped reads for indel candidates [2]
  -M, --max-read-len INT  Maximum length of read to pass to BAQ algorithm [500]
  -o, --open-prob INT     Phred-scaled gap open seq error probability [40]
  -p, --per-sample-mF     Apply -m and -F per-sample for increased sensitivity
  -P, --platforms STR     Comma separated list of platforms for indels [all]
  --ar, --ambig-reads STR   What to do with ambiguous indel reads: drop,incAD,incAD0 [drop]
      --indel-bias FLOAT  Raise to favour recall over precision [1.00]
      --del-bias FLOAT    Relative likelihood of insertion to deletion [0.00]
      --score-vs-ref FLOAT
                          Ratio of score vs ref (1) or 2nd-best allele (0) [0.00]
      --indel-size INT    Approximate maximum indel size considered [110]
      --indels-2.0        New EXPERIMENTAL indel calling model (diploid reference consensus)
      --indels-cns        New EXPERIMENTAL indel calling model with edlib
      --seqq-offset       Indel-cns tuning for indel seq-qual scores [120]
      --no-indels-cns     Disable CNS mode, to use after a -X profile
      --poly-mqual        (Edlib mode) Use minimum quality within homopolymers

Notes: Assuming diploid individuals.

Example:
   # See also http://samtools.github.io/bcftools/howtos/variant-calling.html
   vtools mpileup -Ou -f reference.fa alignments.bam | vtools call -mv -Ob -o calls.bcf

"#;

/// `list_annotations`
const ANNOTATIONS: &str = "Annotations added by default are in this list prefixed with \"*\". To suppress their output, run with
e.g. \"-a -FORMAT/AD\".

FORMAT annotation tags available (\"FORMAT/\" prefix is optional):

* FORMAT/AD   .. Allelic depth (Number=R,Type=Integer)
  FORMAT/ADF  .. Allelic depths on the forward strand (Number=R,Type=Integer)
  FORMAT/ADR  .. Allelic depths on the reverse strand (Number=R,Type=Integer)
  FORMAT/DP   .. Number of high-quality bases (Number=1,Type=Integer)
  FORMAT/NMBZ .. Mann-Whitney U-z test of Number of Mismatches within supporting reads (Number=1,Type=Float)
  FORMAT/QM   .. Allele phred-score quality mean for use with +trio-dnm (Number=R,Type=Integer)
  FORMAT/QS   .. Allele phred-score quality sum for use with `call -mG` and +trio-dnm (Number=R,Type=Integer)
  FORMAT/SP   .. Phred-scaled strand bias P-value (Number=1,Type=Integer)
  FORMAT/SCR  .. Number of soft-clipped reads (Number=1,Type=Integer)

INFO annotation tags available:

  INFO/AD    .. Total allelic depth (Number=R,Type=Integer)
  INFO/ADF   .. Total allelic depths on the forward strand (Number=R,Type=Integer)
  INFO/ADR   .. Total allelic depths on the reverse strand (Number=R,Type=Integer)
* INFO/BQBZ  .. Mann-Whitney U test of Base Quality Bias (Number=1,Type=Float)
  INFO/FS    .. Fisher's exact test P-value to detect strand bias (Number=1,Type=Float)
* INFO/IDV   .. Maximum number of raw reads supporting an indel (Number=1,Type=Integer)
* INFO/IMF   .. Maximum fraction of raw reads supporting an indel (Number=1,Type=Float)
* INFO/MQ0F  .. Fraction of reads with zero mapping quality (Number=1,Type=Float)
* INFO/MQBZ  .. Mann-Whitney U test of Mapping Quality Bias (Number=1,Type=Float)
* INFO/MQSBZ .. Mann-Whitney U-z test of Mapping Quality vs Strand Bias (Number=1,Type=Float)
  INFO/NM    .. Approximate average number of mismatches in ref and alt reads, experimental (Number=2,Type=Float)
  INFO/NMBZ  .. Mann-Whitney U-z test of Number of Mismatches within supporting reads (Number=1,Type=Float)
* INFO/RPBZ  .. Mann-Whitney U test of Read Position Bias (Number=1,Type=Float)
* INFO/SCBZ  .. Mann-Whitney U-z test of Soft-Clip Length Bias (Number=1,Type=Float)
  INFO/SCR   .. Number of soft-clipped reads (Number=1,Type=Integer)
* INFO/SGB   .. Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf (Number=1,Type=Float)
* INFO/VDB   .. Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (Number=1,Type=Float)

";

/// `print_profiles`
const PROFILES: &str = "Configuration profiles activated with -X, --config:

1.12
    -Q13 -h100 -m1 -F0.002

bgi, bgi-1.20
    --indels-cns -B --indel-size 80 -F0.1 --indel-bias 0.9 --seqq-offset 120

illumina-1.18
    --indel-size 110

illumina
illumina-1.20
    --indels-cns --indel-size 110

ont
    -B -Q5 --max-BQ 30 -I

ont-sup, ont-sup-1.20
    --indels-cns -B -Q1 --max-BQ 35 -F0.2 -o15 -e1 -h110 --delta-BQ 99\\
    --del-bias 0.4 --indel-bias 0.7 --poly-mqual --seqq-offset 130\\
    --indel-size 80

pacbio-ccs-1.18
    -D -Q5 --max-BQ 50 -F0.1 -o25 -e1 --delta-BQ 10 \\
    -M99999 --indel-size 110

pacbio-ccs, pacbio-ccs-1.20
    --indels-cns -B -Q5 --max-BQ 50 -F0.1 -o25 -e1 -h300 --delta-BQ 10 \\
    --del-bias 0.4 --poly-mqual --indel-bias 0.9 --seqq-offset 118\\
    --indel-size 80 --score-vs-ref 0.7

ultima, ultima-1.20
    --indels-cns -B -Q1 --max-BQ 30 -F0.15 -o20 -e10 -h250 --delta-BQ 10 \\
    --del-bias 0.3 --indel-bias 0.7 --poly-mqual --seqq-offset 140 \\
    --indel-size 80 --score-vs-ref 0.3


";

/// `-X`: applies a configuration profile.
fn apply_profile(o: &mut Options, name: &str) -> anyhow::Result<()> {
    let n = name.to_ascii_lowercase();
    match n.as_str() {
        "pacbio-ccs-1.18" => {
            o.min_frac = 0.1;
            o.min_baseq = 5;
            o.max_baseq = 50;
            o.delta_baseq = 10;
            o.openq = 25;
            o.extq = 1;
            o.flag |= pileup::REALN_PARTIAL;
            o.max_read_len = 99999;
        }
        "pacbio-ccs" | "pacbio-ccs-1.20" => {
            o.min_frac = 0.1;
            o.min_baseq = 5;
            o.max_baseq = 50;
            o.delta_baseq = 10;
            o.tandemq = 300;
            o.openq = 25;
            o.extq = 1;
            o.flag &= !pileup::REALN;
            o.del_bias = 0.4;
            o.indel_bias = 1.0 / 0.9;
            o.seqq_offset = 118;
            o.poly_mqual = true;
            o.edlib = true;
            o.vs_ref = 0.7;
            o.indel_win_size = 80;
        }
        "ont" => {
            eprintln!(
                "With old ONT data may be beneficial to also run vtools call with a higher -P, eg -P0.01 or -P 0.1"
            );
            o.min_baseq = 5;
            o.max_baseq = 30;
            o.flag &= !pileup::REALN;
            o.flag |= pileup::NO_INDEL;
        }
        "ont-sup" | "ont-sup-1.20" => {
            o.min_frac = 0.2;
            o.min_baseq = 1;
            o.max_baseq = 35;
            o.delta_baseq = 99;
            o.openq = 15;
            o.extq = 1;
            o.flag &= !pileup::REALN;
            o.max_read_len = 9_999_999;
            o.del_bias = 0.4;
            o.poly_mqual = true;
            o.edlib = true;
            o.tandemq = 110;
            o.indel_bias = 1.0 / 0.7;
            o.seqq_offset = 130;
            o.indel_win_size = 80;
        }
        "ultima" | "ultima-1.20" => {
            o.min_frac = 0.15;
            o.min_baseq = 1;
            o.max_baseq = 30;
            o.delta_baseq = 10;
            o.openq = 20;
            o.extq = 10;
            o.tandemq = 250;
            o.flag &= !pileup::REALN;
            o.del_bias = 0.3;
            o.poly_mqual = true;
            o.edlib = true;
            o.indel_bias = 1.0 / 0.7;
            o.seqq_offset = 140;
            o.vs_ref = 0.3;
            o.indel_win_size = 80;
        }
        "1.12" => {
            o.min_frac = 0.002;
            o.min_support = 1;
            o.min_baseq = 13;
            o.tandemq = 100;
            o.flag &= !pileup::REALN_PARTIAL;
            o.flag |= pileup::REALN;
        }
        "illumina-1.18" => {
            o.indel_win_size = 110;
            o.flag |= pileup::REALN_PARTIAL;
        }
        "illumina" | "illumina-1.20" => {
            o.edlib = true;
            o.indel_win_size = 110;
            o.flag |= pileup::REALN_PARTIAL;
            o.indel_bias = 1.0;
            o.seqq_offset = 125;
        }
        "bgi" | "bgi-1.20" => {
            o.min_frac = 0.1;
            o.edlib = true;
            o.indel_bias = 1.0;
            o.seqq_offset = 120;
            o.flag |= pileup::REALN_PARTIAL;
            o.indel_win_size = 80;
        }
        "list" | "help" => {
            print!("{PROFILES}");
            let _ = std::io::stdout().flush();
            return Err(exit1(""));
        }
        _ => {
            return Err(exit1(format!(
                "Unknown configuration name '{name}'\nPlease use '-X list' to show available choices."
            )));
        }
    }
    Ok(())
}

/// The options in command-line order, as `(index, id, value)`; flags have an empty value.
fn ordered_options(matches: &ArgMatches) -> Vec<(usize, &'static str, &str)> {
    const VALUED: &[&str] = &[
        "bam_list",
        "adjust_mq",
        "max_depth",
        "fasta_ref",
        "read_groups",
        "min_mq",
        "min_bq",
        "max_bq",
        "delta_bq",
        "regions",
        "regions_file",
        "skip_all_set",
        "skip_any_set",
        "skip_all_unset",
        "skip_any_unset",
        "samples",
        "samples_file",
        "targets",
        "targets_file",
        "seed",
        "annotate",
        "gvcf",
        "output",
        "o_short",
        "open_prob",
        "output_type",
        "threads",
        "verbosity",
        "write_index",
        "config",
        "ext_prob",
        "gap_frac",
        "tandem_qual",
        "max_idepth",
        "min_ireads",
        "max_read_len",
        "platforms",
        "ambig_reads",
        "indel_bias",
        "del_bias",
        "score_vs_ref",
        "indel_size",
        "seqq_offset",
    ];
    const FLAGS: &[&str] = &[
        "illumina13",
        "count_orphans",
        "no_baq",
        "full_baq",
        "redo_baq",
        "no_reference",
        "ignore_rg",
        "ignore_overlaps",
        "no_version",
        "skip_indels",
        "per_sample",
        "indels_20",
        "indels_cns",
        "no_indels_cns",
        "poly_mqual",
        "no_poly_mqual",
        "removed_u",
    ];
    let mut out = Vec::new();
    for &id in VALUED {
        if let (Some(idx), Some(vals)) = (matches.indices_of(id), matches.get_many::<String>(id)) {
            for (i, v) in idx.zip(vals) {
                out.push((i, id, v.as_str()));
            }
        }
    }
    for &id in FLAGS {
        if matches.get_count(id) > 0
            && let Some(idx) = matches.indices_of(id)
        {
            for i in idx {
                out.push((i, id, ""));
            }
        }
    }
    out.sort_by_key(|(i, _, _)| *i);
    out
}

/// The option loop of `main_mpileup`.
fn parse_options(matches: &ArgMatches) -> anyhow::Result<Options> {
    let mut o = Options::default();
    for (_, id, v) in ordered_options(matches) {
        match id {
            "verbosity" => crate::commands::head::apply_verbosity(v)?,
            "ignore_overlaps" => o.flag &= !pileup::SMART_OVERLAPS,
            "skip_any_unset" => o.skip_any_unset = parse_skip_flags(v, "nf")?,
            "skip_all_unset" => o.skip_all_unset = parse_skip_flags(v, "if")?,
            "skip_any_set" => o.skip_any_set = parse_skip_flags(v, "ef")?,
            "skip_all_set" => o.skip_all_set = parse_skip_flags(v, "df")?,
            "output" => o.output_fname = Some(v.to_owned()),
            "open_prob" => o.openq = atoi(v),
            "ignore_rg" => o.samples.ignore_read_groups(),
            "gvcf" => o.gvcf = Some(Gvcf::new(v)?),
            "fasta_ref" => {
                // HTSlib has printed why the file cannot be loaded.
                o.fasta = Some(Fasta::open(Path::new(v)).map_err(|_| exit1(""))?);
                o.fai_fname = Some(v.to_owned());
            }
            "no_reference" => o.noref = true,
            "no_version" => o.record_cmd_line = false,
            "threads" => {
                // `strtol` into an int: a negative count means no threads.
                o.n_threads = usize::try_from(crate::output::thread_count(v)).unwrap_or(0);
            }
            "max_depth" => o.max_depth = atoi(v),
            // `reg_is_file` stays set once -R is given.
            "regions" => {
                let source = match o.reg {
                    Some((_, ListSource::File)) => ListSource::File,
                    _ => ListSource::Inline,
                };
                o.reg = Some((v.to_owned(), source));
            }
            "regions_file" => o.reg = Some((v.to_owned(), ListSource::File)),
            "targets" | "targets_file" => {
                let source = if id == "targets" {
                    ListSource::Inline
                } else {
                    ListSource::File
                };
                let (include, idx) = read_targets(v, source)?;
                // A later `^` list does not turn the logic back.
                if include {
                    o.bed_logic = true;
                }
                o.bed = Some(idx);
            }
            "platforms" => {}
            "per_sample" => o.flag |= pileup::PER_SAMPLE,
            "no_baq" => o.flag &= !pileup::REALN,
            "full_baq" => o.flag &= !pileup::REALN_PARTIAL,
            "skip_indels" => o.flag |= pileup::NO_INDEL,
            "redo_baq" => o.flag |= pileup::REDO_BAQ,
            "illumina13" => o.flag |= pileup::ILLUMINA13,
            "samples" => o
                .samples
                .add_samples(v, ListSource::Inline)
                .map_err(|_| anyhow::anyhow!("Could not read samples: {v}"))?,
            "samples_file" => o
                .samples
                .add_samples(v, ListSource::File)
                .map_err(|_| anyhow::anyhow!("Could not read samples: {v}"))?,
            "output_type" => {
                check_output_type(v)?;
                o.output_settings.push(OutputSetting::Type(v.to_owned()));
            }
            "adjust_mq" => o.capq_thres = atoi(v),
            "min_mq" => o.min_mq = atoi(v),
            "min_bq" => o.min_baseq = atoi(v),
            "max_bq" => o.max_baseq = atoi(v),
            "delta_bq" => o.delta_baseq = atoi(v),
            "bam_list" => o.file_list = Some(v.to_owned()),
            "o_short" => {
                let (value, end) = pasta::variant::number::strtol(v.as_bytes(), 0);
                if end == v.len() {
                    o.openq = value as i32;
                } else {
                    o.output_fname = Some(v.to_owned());
                }
            }
            "ext_prob" => o.extq = atoi(v),
            "tandem_qual" => o.tandemq = atoi(v),
            "indel_bias" => {
                let f = atof(v);
                o.indel_bias = if f < 1e-2 { 1.0 / 1e2 } else { 1.0 / f };
            }
            "score_vs_ref" => {
                o.vs_ref = atof(v);
                if o.vs_ref > 1.0 {
                    o.vs_ref = 1.0;
                }
            }
            "indel_size" => o.indel_win_size = parse_indel_size(v)?,
            "indels_20" => {
                o.indels_v20 = true;
                o.edlib = false;
            }
            "write_index" => o.write_index = Some(output::parse_write_index(v)?),
            "indels_cns" => {
                o.edlib = true;
                o.indels_v20 = false;
            }
            "no_indels_cns" => o.edlib = false,
            "seqq_offset" => o.seqq_offset = atoi(v).clamp(100, 200),
            "del_bias" => o.del_bias = atof(v),
            "poly_mqual" => o.poly_mqual = true,
            "no_poly_mqual" => o.poly_mqual = false,
            "count_orphans" => o.use_orphan = true,
            "gap_frac" => o.min_frac = atof(v),
            "min_ireads" => o.min_support = atoi(v),
            "max_idepth" => o.max_indel_depth = atoi(v),
            "read_groups" => o.samples.add_read_groups(v)?,
            "annotate" => {
                if v.starts_with('?') {
                    eprint!("{ANNOTATIONS}");
                    return Err(exit1(""));
                }
                parse_format_flag(&mut o.fmt_flag, v)?;
            }
            "max_read_len" => o.max_read_len = atoi(v),
            "config" => apply_profile(&mut o, v)?,
            "seed" => o.seed = i64::from(atoi(v)),
            "ambig_reads" => o.ambig_reads = parse_ambig_reads(v)?,
            "removed_u" => return Err(exit1("Invalid option: 'U'")),
            _ => unreachable!("unhandled option {id}"),
        }
    }
    Ok(o)
}

/// `--skip-*` (bcftools' `--nf`, `--if`, `--ef` and `--df`): the flag bits.
fn parse_skip_flags(v: &str, option: &str) -> anyhow::Result<i32> {
    let flag = str2flag(v);
    if flag < 0 {
        return Err(exit1(format!("Could not parse --{option} {v}")));
    }
    Ok(flag)
}

/// `-t/-T`: whether reads in the targets are included (no `^`), and the targets.
fn read_targets(v: &str, source: ListSource) -> anyhow::Result<(bool, RegIdx)> {
    let (include, body) = match v.strip_prefix('^') {
        Some(rest) => (false, rest),
        None => (true, v),
    };
    let idx = match source {
        ListSource::Inline => RegIdx::from_list(body, LineFormat::Region)
            .map_err(|e| regidx_error(&e, format!("Could not parse the targets: {body}")))?,
        ListSource::File => RegIdx::from_file(body, LineFormat::for_file(body)).map_err(|e| {
            regidx_error(
                &e,
                format!("vtools mpileup: Could not read file \"{body}\""),
            )
        })?,
    };
    Ok((include, idx))
}

/// The smallest `--indel-size`.
const MIN_INDEL_SIZE: i32 = 20;

/// `--indel-size`, raised to [`MIN_INDEL_SIZE`] with a warning.
fn parse_indel_size(v: &str) -> anyhow::Result<i32> {
    let (value, end) = pasta::variant::number::strtol(v.as_bytes(), 0);
    if end != v.len() {
        bail!("Could not parse argument: --indel-size {v}");
    }
    let size = value as i32;
    if size < MIN_INDEL_SIZE {
        eprintln!(
            "Warning: running with --indel-size {MIN_INDEL_SIZE}, the requested value is too small"
        );
        return Ok(MIN_INDEL_SIZE);
    }
    Ok(size)
}

/// `--ambig-reads`
fn parse_ambig_reads(v: &str) -> anyhow::Result<AmbigReads> {
    Ok(if v.eq_ignore_ascii_case("drop") {
        AmbigReads::Drop
    } else if v.eq_ignore_ascii_case("incAD") {
        AmbigReads::IncAd
    } else if v.eq_ignore_ascii_case("incAD0") {
        AmbigReads::IncAd0
    } else {
        bail!("The option to --ambig-reads not recognised: {v}");
    })
}

/// Reports a region list that `regidx_init` rejects: its own message (a line it cannot parse,
/// or a file it cannot open), then mpileup's.
fn regidx_error(e: &anyhow::Error, msg: String) -> anyhow::Error {
    eprintln!("{e:#}");
    exit1(msg)
}

/// `read_file_list`
fn read_file_list(path: &str) -> anyhow::Result<Vec<String>> {
    let text = std::fs::read(path)
        .map_err(|e| exit1(format!("{path}: {}", crate::raw_input::strerror(&e))))?;
    let mut files = Vec::new();
    for line in text.split(|&c| c == b'\n') {
        let mut len = line.len();
        while len > 0 && matches!(line[len - 1], b' ' | b'\t'..=b'\r') {
            len -= 1;
        }
        if len == 0 {
            continue;
        }
        let name = String::from_utf8_lossy(&line[..len]).into_owned();
        let is_url = name.find(':').is_some_and(|i| {
            name[..i]
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"+.-".contains(&c))
        });
        if !is_url && !Path::new(&name).exists() {
            let printable = line[..len].iter().all(|&c| (b' '..=b'~').contains(&c));
            return Err(exit1(if printable {
                format!("The file list \"{path}\" appears broken, could not locate: {name}")
            } else {
                format!("Does the file \"{path}\" really contain a list of files and do all exist?")
            }));
        }
        files.push(name);
    }
    if files.is_empty() {
        return Err(exit1(format!("No files read from {path}")));
    }
    Ok(files)
}

/// The header lines of the annotations, each with the `-a` flag that adds it (`None` for
/// lines always added), in bcftools' order.
const ANNOTATION_HEADER_LINES: &[(Option<u32>, &str)] = &[
    (
        Some(call::INFO_IDV),
        "##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of raw reads supporting an indel\">",
    ),
    (
        Some(call::INFO_IMF),
        "##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of raw reads supporting an indel\">",
    ),
    (
        None,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">",
    ),
    (
        Some(call::INFO_VDB),
        "##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">",
    ),
    (
        Some(call::INFO_RPBZ),
        "##INFO=<ID=RPBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)\">",
    ),
    (
        Some(call::INFO_MQBZ),
        "##INFO=<ID=MQBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)\">",
    ),
    (
        Some(call::INFO_BQBZ),
        "##INFO=<ID=BQBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)\">",
    ),
    (
        Some(call::INFO_MQSBZ),
        "##INFO=<ID=MQSBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)\">",
    ),
    (
        Some(call::INFO_NM),
        "##INFO=<ID=NM,Number=2,Type=Float,Description=\"Average number of mismatches in ref and alt reads (approximate, experimental, make me localized?)\">",
    ),
    (
        Some(call::INFO_NMBZ),
        "##INFO=<ID=NMBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better; approximate, experimental, make me localized?)\">",
    ),
    (
        Some(call::FMT_NMBZ),
        "##FORMAT=<ID=NMBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better)\">",
    ),
    (
        Some(call::INFO_SCBZ),
        "##INFO=<ID=SCBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)\">",
    ),
    (
        Some(call::INFO_FS),
        "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher's exact test P-value to detect strand bias\">",
    ),
    (
        Some(call::INFO_SGB),
        "##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf\">",
    ),
    (
        Some(call::INFO_MQ0F),
        "##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">",
    ),
    (
        None,
        "##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">",
    ),
    (
        None,
        "##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">",
    ),
    (
        None,
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">",
    ),
    (
        Some(call::FMT_DP),
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">",
    ),
    (
        Some(call::FMT_DV),
        "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">",
    ),
    (
        Some(call::FMT_DPR),
        "##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">",
    ),
    (
        Some(call::INFO_DPR),
        "##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">",
    ),
    (
        Some(call::FMT_DP4),
        "##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">",
    ),
    (
        Some(call::FMT_SP),
        "##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">",
    ),
    (
        Some(call::FMT_AD),
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (high-quality bases)\">",
    ),
    (
        Some(call::FMT_ADF),
        "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand (high-quality bases)\">",
    ),
    (
        Some(call::FMT_ADR),
        "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand (high-quality bases)\">",
    ),
    (
        Some(call::FMT_QM),
        "##FORMAT=<ID=QM,Number=R,Type=Integer,Description=\"Phred-score allele quality mean used by `+trio-dnm`\">",
    ),
    (
        Some(call::FMT_QS),
        "##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Phred-score allele quality sum used by `call -mG` and `+trio-dnm`\">",
    ),
    (
        Some(call::INFO_AD),
        "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths (high-quality bases)\">",
    ),
    (
        Some(call::INFO_ADF),
        "##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand (high-quality bases)\">",
    ),
    (
        Some(call::INFO_SCR),
        "##INFO=<ID=SCR,Number=1,Type=Integer,Description=\"Number of soft-clipped reads (at high-quality bases)\">",
    ),
    (
        Some(call::FMT_SCR),
        "##FORMAT=<ID=SCR,Number=1,Type=Integer,Description=\"Per-sample number of soft-clipped reads (at high-quality bases)\">",
    ),
    (
        Some(call::INFO_ADR),
        "##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand (high-quality bases)\">",
    ),
];

/// Builds the VCF header (`mpileup`'s header creation).
fn build_header(
    o: &Options,
    contigs: &[(Vec<u8>, i64)],
    samples: &[Vec<u8>],
) -> anyhow::Result<Header> {
    let mut h = Header::new_for_writing();
    let mut add = |line: &str| -> anyhow::Result<()> {
        h.append(line.as_bytes())
            .with_context(|| format!("[mpileup] failed to update the header: {line}"))?;
        Ok(())
    };
    if o.record_cmd_line {
        add(&format!(
            "##vtoolsVersion={}+pasta-{}",
            env!("CARGO_PKG_VERSION"),
            env!("CARGO_PKG_VERSION")
        ))?;
        let mut cmd = String::from("##vtoolsCommand=mpileup");
        for a in std::env::args().skip(2) {
            cmd.push(' ');
            cmd.push_str(&a);
        }
        add(&cmd)?;
    }
    if let Some(f) = &o.fai_fname {
        add(&format!("##reference=file://{f}"))?;
    }
    for (name, len) in contigs {
        add(&format!(
            "##contig=<ID={},length={len}>",
            String::from_utf8_lossy(name)
        ))?;
    }
    let ff = o.fmt_flag;
    add("##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">")?;
    add(
        "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">",
    )?;
    for &(bit, line) in ANNOTATION_HEADER_LINES {
        if bit.is_none_or(|bit| ff & bit != 0) {
            add(line)?;
        }
    }
    if o.gvcf.is_some() {
        Gvcf::update_header(&mut h)?;
    }
    for s in samples {
        h.add_sample(s)?;
    }
    Ok(h)
}

/// Writes records directly, through the gVCF blocks if any (`flush_bcf_records`).
struct WriterSink<'w> {
    writer: &'w mut Writer,
    gvcf: Option<Gvcf>,
}

impl WriterSink<'_> {
    fn finish(&mut self, hdr: &Header) -> anyhow::Result<()> {
        if let Some(g) = &mut self.gvcf {
            let writer = &mut *self.writer;
            g.flush(hdr, |line| {
                writer.write_record(hdr, line)?;
                Ok(())
            })?;
        }
        Ok(())
    }
}

impl Sink for WriterSink<'_> {
    fn emit(&mut self, hdr: &Header, rec: &mut Record) -> anyhow::Result<()> {
        let Some(g) = &mut self.gvcf else {
            self.writer.write_record(hdr, rec)?;
            return Ok(());
        };
        let writer = &mut *self.writer;
        let write_rec = g.write(hdr, rec, |line| {
            writer.write_record(hdr, line)?;
            Ok(())
        })?;
        if write_rec {
            self.writer.write_record(hdr, rec)?;
        }
        Ok(())
    }
}

pub fn run(_args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    match run_inner(matches) {
        Ok(code) => Ok(code),
        Err(e) => match e.downcast_ref::<Exit1>() {
            Some(Exit1(msg)) => {
                if !msg.is_empty() {
                    eprintln!("{msg}");
                }
                Ok(ExitCode::from(1))
            }
            None => Err(e),
        },
    }
}

/// Whether the command line gives no option or input, bcftools' `argc == 1`. This reads the
/// parsed arguments, not the process's, which differ when vtools runs inside another program
/// (its tests).
fn no_arguments(matches: &ArgMatches) -> bool {
    matches
        .ids()
        .all(|id| matches.value_source(id.as_str()) != Some(clap::parser::ValueSource::CommandLine))
}

fn run_inner(matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    let mut o = parse_options(matches)?;
    if o.gvcf.is_some() && o.fmt_flag & call::FMT_DP == 0 {
        eprintln!("[warning] The -a DP option is required with --gvcf, switching on.");
        o.fmt_flag |= call::FMT_DP;
    }
    if o.flag & pileup::REALN == 0 && o.flag & pileup::REDO_BAQ != 0 {
        return Err(exit1("Error: The -B option cannot be combined with -E"));
    }
    if o.use_orphan {
        o.flag &= !pileup::NO_ORPHAN;
    }
    if no_arguments(matches) {
        eprint!("{USAGE}");
        return Err(exit1(""));
    }
    if o.fai_fname.is_none() && !o.noref {
        return Err(exit1(
            "Error: mpileup requires the --fasta-ref option by default; use --no-reference to run without a fasta reference",
        ));
    }
    let files: Vec<String> = match &o.file_list {
        Some(list) => read_file_list(list)?,
        None => matches
            .get_many::<String>("inputs")
            .map(|v| v.cloned().collect())
            .unwrap_or_default(),
    };
    if files.is_empty() {
        return Err(exit1("[mpileup] no input file/data given"));
    }
    mpileup(o, files)
}

/// The error for a region whose sequence is not in the header. bcftools queries
/// `SEQ:BEG-END` and then `SEQ` alone: a sequence name that itself reads as a region of a
/// known contig gets a different message.
fn region_not_found(
    header: &input::AlnHeader,
    seq: &[u8],
    beg: u32,
    end: u32,
    file: &str,
) -> anyhow::Error {
    let seq = String::from_utf8_lossy(seq);
    let parses_as_region = seq
        .rfind(':')
        .is_some_and(|i| header.tid(&seq.as_bytes()[..i]).is_some());
    if parses_as_region {
        exit1(format!(
            "[E::mpileup] fail to parse region '{seq}:{}-{}'",
            u64::from(beg) + 1,
            u64::from(end) + 1
        ))
    } else {
        exit1(format!(
            "[E::mpileup] the sequence \"{seq}\" not found: {file}"
        ))
    }
}

/// A region of `-r/-R`: the sequence name and the 0-based start and inclusive end.
type Region = (Vec<u8>, u32, u32);

/// A span to pile up, with its 0-based start and inclusive end.
type SpanRange = (Span, u32, u32);

/// `-r/-R`: the regions in order, overlapping ones merged.
fn read_regions(o: &Options) -> anyhow::Result<Option<Vec<Region>>> {
    let Some((spec, source)) = &o.reg else {
        return Ok(None);
    };
    let idx = match source {
        ListSource::File => RegIdx::from_file(spec, LineFormat::for_file(spec)),
        ListSource::Inline => RegIdx::from_list(spec, LineFormat::Region),
    };
    let mut idx =
        idx.map_err(|e| regidx_error(&e, format!("Could not parse the regions: {spec}")))?;
    idx.set_merge_overlaps();
    Ok(Some(idx.all_regions()))
}

/// The inputs with usable read groups.
struct Inputs {
    files: Vec<InputFile>,
    /// The first file's header: contig names for every file.
    first: input::AlnHeader,
    /// The number of leading regions whose sequences every file has.
    contigs_ok_upto: usize,
    /// Whether every file is indexed and has the first file's contigs, so chromosomes can
    /// be piled up separately (never for a serial run, which does not look for the indexes
    /// without regions).
    all_indexed: bool,
    /// The file descriptors that readers of every file use.
    readers: input::Descriptors,
}

/// Opens the inputs one after another and registers their samples. With regions, each file
/// must be indexed and have the first region's sequence, and its index is loaded once it is
/// opened, as bcftools does.
///
/// For a serial run (`serial_io_threads` gives its decompression threads), each file stays
/// open for the pileup, as in bcftools. A run that may be on threads opens the inputs again
/// for each chunk; it loads the indexes also without regions, to pile up the chromosomes
/// separately, and closes each file once read, so that it holds no more descriptors than
/// bcftools, which loads no index then.
fn open_inputs(
    samples: &mut Samples,
    files: &[String],
    regions: Option<&[Region]>,
    serial_io_threads: Option<usize>,
    cram_reference: Option<&Path>,
) -> anyhow::Result<Inputs> {
    let mut inputs: Vec<InputFile> = Vec::new();
    let mut first: Option<input::AlnHeader> = None;
    let mut contigs_ok_upto = usize::MAX;
    let mut all_indexed = true;
    let mut readers = input::Descriptors::default();
    for name in files {
        let (reader, header) =
            AlnReader::open(name, cram_reference, serial_io_threads.unwrap_or(0))
                .map_err(|e| exit1(format!("{e:#}")))?;
        let Some(bam_id) = samples.add_bam(Some(&header.text), name)? else {
            continue;
        };
        // Reads are piled up by reference sequence index, named from the first file's header,
        // so every file's sequences must be the first file's, or the first of them. (bcftools
        // mixes up files whose sequences differ.)
        if let Some(f) = &first
            && (!f.names.starts_with(&header.names) || !f.lengths.starts_with(&header.lengths))
        {
            return Err(exit1(format!(
                "[mpileup] the reference sequences of {name} differ from those of {} \
                 (they must be the first file's, or the first of them, in the same order)",
                inputs[0].name
            )));
        }
        let index = if regions.is_some() || serial_io_threads.is_none() {
            reader.load_index()
        } else {
            None
        };
        let has_index = index.is_some();
        readers += reader.descriptors();
        // Chromosomes are piled up separately only when every file has the same contigs.
        all_indexed &= has_index && first.as_ref().is_none_or(|f| f.names == header.names);
        if let Some(regs) = regions
            && !regs.is_empty()
        {
            if !has_index {
                return Err(exit1(format!("[mpileup] fail to load index for {name}")));
            }
            let h = first.as_ref().unwrap_or(&header);
            for (k, (seq, beg, end)) in regs.iter().enumerate() {
                if h.tid(seq).is_none() && k < contigs_ok_upto {
                    if k == 0 {
                        return Err(region_not_found(h, seq, *beg, *end, name));
                    }
                    contigs_ok_upto = k;
                }
            }
        }
        // A serial run keeps every reader for the pileup. Standard input, a pipe, and a FIFO
        // cannot be opened again (the bytes are gone, or `open` waits for a writer): keep
        // their reader for any run.
        let keep = serial_io_threads.is_some() || is_stream(name);
        let opened = keep.then(|| (reader, header.clone()));
        inputs.push(InputFile {
            name: name.clone(),
            bam_id,
            opened: std::sync::Mutex::new(opened),
            index,
        });
        if first.is_none() {
            first = Some(header);
        }
    }
    if regions.is_none() && !all_indexed {
        // Read whole: no chromosome is queried.
        for file in &mut inputs {
            file.index = None;
        }
    }
    let Some(first) = first else {
        return Err(exit1(
            "[mpileup] failed to find a file header with usable read groups",
        ));
    };
    Ok(Inputs {
        files: inputs,
        first,
        contigs_ok_upto,
        all_indexed,
        readers,
    })
}

/// Opens the output and writes the header to it, starting the index if `-W` asks for one.
fn open_output(
    o: &Options,
    first: &input::AlnHeader,
) -> anyhow::Result<(Writer, Header, OutputType)> {
    let (out_type, level) =
        output::resolve(&o.output_settings, o.output_fname.as_deref().map(Path::new))?;
    let out_path = o
        .output_fname
        .as_deref()
        .map(|p| output::data_path(Path::new(p)));
    let mut writer = Writer::create(
        out_path.as_deref(),
        &output::write_options(out_type, level, o.n_threads),
    )
    .map_err(|e| {
        exit1(format!(
            "[mpileup] failed to write to {}: {e}",
            o.output_fname.as_deref().unwrap_or("standard output")
        ))
    })?;
    let contigs: Vec<(Vec<u8>, i64)> = first
        .names
        .iter()
        .cloned()
        .zip(first.lengths.iter().copied())
        .collect();
    let hdr = build_header(o, &contigs, o.samples.names())?;
    writer.write_header(&hdr)?;
    // `init_index2`, after the header as in bcftools: standard output and uncompressed VCF
    // cannot be indexed.
    if let Some(fmt) = o.write_index {
        let indexed =
            output::index_request(fmt, out_type, o.output_fname.as_deref().map(Path::new))
                .is_some_and(|req| {
                    writer.request_index(req).is_ok() && writer.start_index(&hdr).is_ok()
                });
        if !indexed {
            writer.abandon();
            bail!(
                "Error: failed to initialise index for {}",
                o.output_fname.as_deref().unwrap_or("(null)")
            );
        }
    }
    Ok((writer, hdr, out_type))
}

/// The per-sample depth limit below which bcftools notes what `-d` allows per sample.
const NOTE_PER_SAMPLE_DEPTH: f64 = 250.0;

/// The pileup's depth limit per input (`-d`, 0 for none), after bcftools' notes about it.
fn max_depth(o: &Options, n_files: usize, nsmpl: usize) -> usize {
    if o.max_depth == 0 {
        eprintln!("[mpileup] Max depth set to maximum value ({})", i32::MAX);
        return i32::MAX as usize;
    }
    eprintln!(
        "[mpileup] maximum number of reads per input file set to -d {}",
        o.max_depth
    );
    let total = f64::from(o.max_depth) * n_files as f64;
    if total > f64::from(1 << 20) {
        eprintln!(
            "Warning: Potential memory hog, up to {:.0}M reads in the pileup!",
            total
        );
    }
    if total / (nsmpl as f64) < NOTE_PER_SAMPLE_DEPTH {
        eprintln!(
            "Note: The maximum per-sample depth with -d {} is {:.1}x",
            o.max_depth,
            total / nsmpl as f64
        );
    }
    o.max_depth.max(0) as usize
}

/// The pileup's and the caller's settings (`mplp_conf_t` and `bcf_call_init`), with a budget
/// of `threads` CPU threads. Takes the targets out of `o`.
fn pileup_conf(o: &mut Options, max_depth: usize, nsmpl: usize, threads: usize) -> Conf {
    let call = CallSettings {
        fmt_flag: o.fmt_flag,
        ambig_reads: o.ambig_reads,
        min_baseq: o.min_baseq,
        max_baseq: o.max_baseq,
        delta_baseq: o.delta_baseq,
        openq: o.openq,
        extq: o.extq,
        tandemq: o.tandemq,
        min_support: o.min_support as u32,
        min_frac: o.min_frac,
        per_sample_flt: o.flag & pileup::PER_SAMPLE != 0,
        indel_win_size: o.indel_win_size,
        indels_v20: o.indels_v20,
        edlib: o.edlib,
        seqq_offset: o.seqq_offset,
        poly_mqual: o.poly_mqual,
        indel_bias: o.indel_bias as f32,
        del_bias: o.del_bias as f32,
        vs_ref: o.vs_ref as f32,
        threads: ThreadBudget::new(threads),
    };
    let mut bed = o.bed.take();
    if let Some(b) = &mut bed {
        b.build();
    }
    Conf {
        flag: o.flag,
        min_mq: o.min_mq,
        capq_thres: o.capq_thres,
        max_depth,
        max_indel_depth: o.max_indel_depth.saturating_mul(nsmpl as i32),
        max_read_len: o.max_read_len,
        skip_any_unset: o.skip_any_unset as u16,
        skip_all_unset: o.skip_all_unset as u16,
        skip_any_set: o.skip_any_set as u16,
        skip_all_set: o.skip_all_set as u16,
        bed,
        bed_logic: o.bed_logic,
        call,
    }
}

/// The spans to pile up: each region, or the whole input (by chromosome with threads). A
/// region whose sequence is missing ends the list, and its error is returned to report
/// after the regions before it are written.
fn plan_spans(
    o: &Options,
    regions: Option<&[Region]>,
    inputs: &Inputs,
) -> (Vec<SpanRange>, Option<anyhow::Error>) {
    let first = &inputs.first;
    let mut spans: Vec<SpanRange> = Vec::new();
    let Some(regs) = regions else {
        let by_chromosome =
            o.n_threads > 0 && o.gvcf.is_none() && !first.names.is_empty() && inputs.all_indexed;
        if by_chromosome {
            spans.extend((0..first.names.len()).map(|tid| {
                (
                    Span::Region {
                        tid,
                        beg: 0,
                        end: REGIDX_MAX,
                    },
                    0u32,
                    u32::MAX,
                )
            }));
        } else {
            spans.push((Span::All, 0, u32::MAX));
        }
        return (spans, None);
    };
    for (k, (seq, beg, end)) in regs.iter().enumerate() {
        if k >= inputs.contigs_ok_upto {
            let e = region_not_found(first, seq, *beg, *end, &inputs.files[0].name);
            return (spans, Some(e));
        }
        let tid = first.tid(seq).expect("checked above");
        spans.push((
            Span::Region {
                tid,
                beg: *beg,
                end: *end,
            },
            *beg,
            *end,
        ));
    }
    (spans, None)
}

/// `mpileup()`
fn mpileup(mut o: Options, files: Vec<String>) -> anyhow::Result<ExitCode> {
    // The error model's tables take a while; build them while the inputs are opened.
    // Without a thread to start, the tables are built when first used.
    let _ = thread::Builder::new().spawn(|| {
        std::sync::LazyLock::force(&errmod::ERROR_MODEL);
    });
    let regions = read_regions(&o)?;
    let cram_reference = o.fai_fname.as_ref().map(PathBuf::from);
    // Without threads, or with gVCF output, the run is serial: the inputs stay open from
    // their headers on, with the decompression threads shared out between them (each is read
    // on a thread of its own, or on this one when that thread cannot start).
    let serial = o.n_threads == 0 || o.gvcf.is_some();
    let serial_io_threads = serial.then(|| o.n_threads / files.len().max(1));
    let inputs = open_inputs(
        &mut o.samples,
        &files,
        regions.as_deref(),
        serial_io_threads,
        cram_reference.as_deref(),
    )?;
    let nsmpl = o.samples.names().len();
    eprintln!(
        "[mpileup] {nsmpl} samples in {} input files",
        inputs.files.len()
    );
    let (mut writer, hdr, out_type) = open_output(&o, &inputs.first)?;
    let max_depth = max_depth(&o, inputs.files.len(), nsmpl);
    let refs = References::new(o.fasta.take(), inputs.first.names.clone());
    let (spans, pending_error) = plan_spans(&o, regions.as_deref(), &inputs);
    // Piling up on threads opens every input once per chunk in flight; where the file
    // descriptors do not allow even one chunk, the serial run opens each input once, as
    // bcftools does.
    let window = if o.n_threads > 0 && !spans.is_empty() && o.gvcf.is_none() {
        parallel::window(o.n_threads, inputs.readers, parallel::free_descriptors())
    } else {
        None
    };
    // On threads, each of the N workers holds a thread of the budget while it piles up a
    // chunk, and heavy indel positions borrow the other N (and those of idle workers); the
    // serial run holds one and lends the rest.
    let budget = if window.is_some() {
        2 * o.n_threads
    } else {
        o.n_threads
    };
    let conf = pileup_conf(&mut o, max_depth, nsmpl, budget);

    let shared = Shared {
        conf: &conf,
        files: &inputs.files,
        samples: &o.samples,
        header: &inputs.first,
        refs: &refs,
        out_header: &hdr,
        cram_reference: cram_reference.as_deref(),
        io_threads: 0,
        nsmpl,
    };
    let piled = if let Some(window) = window {
        let encoding = if out_type.is_vcf() {
            Encoding::Vcf
        } else {
            Encoding::Bcf
        };
        parallel::run(
            &shared,
            &spans,
            o.n_threads,
            window,
            o.seed,
            &mut writer,
            encoding,
        )
    } else {
        Ok(Piled::NoThreads)
    };
    let result = if let Ok(Piled::NoThreads) = piled {
        // The inputs not kept open (where no chunk fits under the open-file limit, or no
        // worker thread starts) are opened again with the decompression threads shared out
        // between them.
        let shared = Shared {
            io_threads: o.n_threads / inputs.files.len(),
            ..shared
        };
        let _busy = conf.call.threads.hold();
        let mut engine = Engine::new(&conf.call, nsmpl, Rng48::new(o.seed), RngSource::Known);
        let mut sink = WriterSink {
            writer: &mut writer,
            gvcf: o.gvcf.take(),
        };
        pileup::run(&shared, &spans, &mut engine, &mut sink).and_then(|()| sink.finish(&hdr))
    } else {
        piled.map(drop)
    };
    // bcftools' `flush_bcf_records` error for a record it cannot write.
    let name = o.output_fname.as_deref().unwrap_or("standard output");
    let record_failure =
        || anyhow::anyhow!("[flush_bcf_records] Error: failed to write the record to {name}");
    let result = result.and_then(|()| match pending_error {
        Some(e) => Err(e),
        None => Ok(()),
    });
    if let Err(e) = result {
        writer.abandon();
        return Err(match e.downcast::<WriteError>() {
            Ok(e) => output::write_failure(e, record_failure),
            Err(e) => e,
        });
    }
    let close_failed = || {
        let name = o.output_fname.as_deref().unwrap_or("(null)");
        anyhow::anyhow!("[mpileup] Error: close failed .. {name}")
    };
    match writer.finish() {
        Ok(()) => Ok(ExitCode::SUCCESS),
        Err(e @ WriteError::Deferred(_)) => Err(output::write_failure(e, record_failure)),
        Err(e @ WriteError::Io(_)) => Err(output::write_failure(e, close_failed)),
        Err(e) => Err(e.into()),
    }
}

#[cfg(test)]
mod tests {
    use std::thread;

    use pasta::variant::io::reader::open_again;

    use super::*;
    use crate::raw_input::{make_fifo, within_seconds};

    #[test]
    fn keeps_the_reader_of_a_fifo_for_the_pileup() {
        let (dir, fifo) = make_fifo("mpileup-fifo");
        let name = fifo.to_string_lossy().into_owned();
        let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n\
                   r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t*\n";
        // Opening the FIFO to write waits for the reader below, then the writer leaves:
        // opening it again would wait forever.
        let feeder = {
            let fifo = fifo.clone();
            thread::spawn(move || std::fs::write(&fifo, sam))
        };
        let kept = within_seconds("opening a SAM FIFO", move || {
            let inputs = open_inputs(&mut Samples::default(), &[name], None, None, None)?;
            let kept = inputs.files[0]
                .opened
                .lock()
                .expect("input lock poisoned")
                .is_some();
            anyhow::Ok(kept)
        });
        // A reader for the writer, should the FIFO not have been opened.
        let partner = open_again(&fifo);
        let fed = feeder.join().expect("the writer thread does not panic");
        drop(partner);
        let _ = std::fs::remove_dir_all(&dir);
        fed.expect("write the SAM text");
        let kept = kept.expect("the FIFO opens");
        assert!(kept, "the pileup would open the FIFO again");
    }

    #[test]
    fn checks_output_types_as_mpileup_does() {
        for ok in ["v", "z", "b", "u", "z5", "b0", "3"] {
            assert!(check_output_type(ok).is_ok(), "-O {ok}");
        }
        let msg = |v: &str| check_output_type(v).map_err(|e| e.to_string()).unwrap_err();
        assert!(check_output_type("").is_err(), "-O ''");
        assert_eq!(msg("q"), "The output type \"q\" not recognised");
        assert_eq!(msg("12"), "The output type \"12\" not recognised");
        assert_eq!(msg("vx"), "Could not parse argument: --output-type x");
        assert_eq!(msg("z10"), "Could not parse argument: --output-type 10");
    }

    /// Where the helper test below finds its inputs and writes its outputs, and the limit of
    /// open files it runs under; set only in the child process that runs it.
    const LOW_NOFILE_DIR: &str = "VTOOLS_TEST_LOW_NOFILE_DIR";
    const LOW_NOFILE_LIMIT: &str = "VTOOLS_TEST_LOW_NOFILE_LIMIT";

    /// Runs `vtools ARGS` in this process.
    fn run_vtools(args: &[&str]) -> anyhow::Result<ExitCode> {
        use clap::{CommandFactory, FromArgMatches};

        let args = args.iter().map(std::ffi::OsString::from).collect();
        let crate::TopLevel::Run(args) = crate::top_level(args) else {
            panic!("vtools exits before running the subcommand");
        };
        let matches = crate::Cli::command().try_get_matches_from(&args)?;
        let cli = crate::Cli::from_arg_matches(&matches)?;
        let sub = matches
            .subcommand()
            .map(|(_, sub)| sub)
            .expect("clap requires a subcommand");
        crate::run(cli.command, sub)
    }

    /// Writes a coordinate-sorted BAM file of 50 bp reads every 1000 bp of a 500 kb contig, and
    /// its index.
    fn write_indexed_bam(path: &Path) {
        let positions: Vec<_> = (1..499_900).step_by(1000).map(|pos| ("c1", pos)).collect();
        write_bam(path, &[("c1", 500_000)], &positions);
    }

    /// Writes a BAM file of 50 bp reads at `reads` (sequence and 1-based position, sorted)
    /// with `sequences` in its header, and its index.
    fn write_bam(path: &Path, sequences: &[(&str, u64)], reads: &[(&str, usize)]) {
        use pasta::{bam, sam, sam::alignment::io::Write as _};

        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
        for (name, len) in sequences {
            text += &format!("@SQ\tSN:{name}\tLN:{len}\n");
        }
        let seq = "ACGTTGCA".repeat(7);
        for (i, (name, pos)) in reads.iter().enumerate() {
            text += &format!(
                "r{i}\t0\t{name}\t{pos}\t60\t50M\t*\t0\t0\t{}\t{}\n",
                &seq[..50],
                "I".repeat(50)
            );
        }
        let mut reader = sam::io::Reader::new(text.as_bytes());
        let header = reader.read_header().expect("parse the SAM header");
        let file = std::fs::File::create(path).expect("create the BAM file");
        let mut writer = bam::io::Writer::new(file);
        writer.write_header(&header).expect("write the BAM header");
        for record in reader.record_bufs(&header) {
            let record = record.expect("parse a SAM record");
            writer
                .write_alignment_record(&header, &record)
                .expect("write a BAM record");
        }
        writer.try_finish().expect("finish the BAM file");
        let bam::Index::Bai(index) = bam::fs::index(path).expect("index the BAM file") else {
            panic!("a short contig gets a BAI index");
        };
        let mut bai = path.as_os_str().to_owned();
        bai.push(".bai");
        bam::bai::fs::write(bai, &index).expect("write the BAI index");
    }

    /// Some of GRCh38's sequences.
    const GRCH38: [(&str, u64); 4] = [
        ("chr1", 248_956_422),
        ("chr2", 242_193_529),
        ("chr3", 198_295_559),
        ("chrM", 16_569),
    ];

    #[test]
    fn sparse_inputs_pile_up_on_threads_in_a_chunk_per_window_of_records() {
        let dir =
            std::env::temp_dir().join(format!("vtools-mpileup-sparse-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        // Two windows of reads in each of 3 inputs, and an input with none. Each window is in
        // one 16 kb bin of the index, and in one 200 kb chunk.
        let reads: Vec<_> = (1_100_001..1_110_000)
            .step_by(40)
            .map(|pos| ("chr1", pos))
            .chain(
                (150_100_001..150_105_000)
                    .step_by(40)
                    .map(|pos| ("chr3", pos)),
            )
            .collect();
        let mut files = Vec::new();
        for i in 0..4 {
            let path = dir.join(format!("in{i}.bam"));
            write_bam(&path, &GRCH38, if i < 3 { &reads } else { &[] });
            files.push(path.display().to_string());
        }
        let plan = |files: &[String]| {
            let inputs = open_inputs(&mut Samples::default(), files, None, None, None)
                .expect("open the inputs");
            assert!(inputs.all_indexed);
            let spans: Vec<SpanRange> = (0..GRCH38.len())
                .map(|tid| {
                    let span = Span::Region {
                        tid,
                        beg: 0,
                        end: REGIDX_MAX,
                    };
                    (span, 0, u32::MAX)
                })
                .collect();
            parallel::plan_inputs(&inputs.files, &inputs.first.lengths, &spans, 4, 2_000).len()
        };
        // Each chunk's walker opens every input (and the main thread's, if it carries on the
        // pileup through the chunk): a chunk per window of records, rather than per 200 kb
        // of the sequences (2,446 chunks), and none for the empty input.
        let chunks = plan(&files);
        let empty_chunks = plan(&files[3..]);
        let run = |inputs: &[String], threads: &str| {
            let out = dir.join(format!("out-{}-{threads}.vcf", inputs.len()));
            let out = out.display().to_string();
            let mut args = vec![
                "vtools",
                "mpileup",
                "--no-version",
                "--no-reference",
                "--threads",
                threads,
                "-o",
                &out,
            ];
            args.extend(inputs.iter().map(String::as_str));
            let code = run_vtools(&args).expect("mpileup succeeds");
            assert_eq!(code, ExitCode::SUCCESS, "--threads {threads}");
            std::fs::read(&out).expect("read the output")
        };
        let outputs: Vec<_> = [&files[..], &files[3..]]
            .iter()
            .map(|inputs| {
                let serial = run(inputs, "0");
                let threaded = ["1", "4"].map(|t| run(inputs, t) == serial);
                (serial, threaded)
            })
            .collect();
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(chunks, 2);
        assert_eq!(empty_chunks, 0);
        let records = |vcf: &[u8]| {
            vcf.split(|&b| b == b'\n')
                .filter(|l| l.first().is_some_and(|&b| b != b'#'))
                .count()
        };
        assert!(
            records(&outputs[0].0) > 10_000,
            "the serial run writes records"
        );
        assert_eq!(records(&outputs[1].0), 0, "no records from the empty input");
        for (serial, threaded) in &outputs {
            assert!(!serial.is_empty());
            assert_eq!(
                threaded,
                &[true, true],
                "--threads 1 and 4 write the serial output"
            );
        }
    }

    #[test]
    fn many_inputs_pile_up_on_threads_under_a_low_open_file_limit() {
        let dir =
            std::env::temp_dir().join(format!("vtools-mpileup-nofile-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        for i in 0..8 {
            write_indexed_bam(&dir.join(format!("in{i}.bam")));
        }
        // The contig makes 25 chunks at --threads 4, which could keep up to 10 sets of the 8
        // inputs open at once. 30 descriptors leave no room for a second set, so the inputs
        // run serially; 80 leave room for a few chunks at once, as they are written.
        let exe = std::env::current_exe().expect("the test binary's path");
        let mut failures = Vec::new();
        for limit in [30, 80] {
            let status = std::process::Command::new(&exe)
                .args([
                    "commands::mpileup::tests::piles_up_under_the_open_file_limit",
                    "--exact",
                    "--ignored",
                    "--quiet",
                ])
                .env(LOW_NOFILE_DIR, &dir)
                .env(LOW_NOFILE_LIMIT, limit.to_string())
                .stdout(std::process::Stdio::null())
                .stderr(std::process::Stdio::null())
                .status()
                .expect("run the test binary");
            if !status.success() {
                failures.push(format!("limit {limit}: {status}"));
            }
        }
        let read = |name: &str| std::fs::read(dir.join(name)).unwrap_or_default();
        let serial = read("out-unlimited-0.vcf");
        let outputs: Vec<_> = ["30-0", "30-4", "80-4"]
            .iter()
            .map(|run| (run, read(&format!("out-{run}.vcf")) == serial))
            .collect();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(failures.is_empty(), "mpileup failed: {failures:?}");
        assert!(serial.len() > 1000, "the serial run writes records");
        for (run, same) in outputs {
            assert!(same, "limit-threads {run} differs from the serial run");
        }
    }

    /// Run by the test above in a child process.
    #[test]
    #[ignore = "run in a child process by many_inputs_pile_up_on_threads_under_a_low_open_file_limit"]
    fn piles_up_under_the_open_file_limit() {
        let (Some(dir), Some(limit)) = (
            std::env::var_os(LOW_NOFILE_DIR),
            std::env::var(LOW_NOFILE_LIMIT).ok(),
        ) else {
            return;
        };
        let dir = PathBuf::from(dir);
        let inputs: Vec<String> = (0..8)
            .map(|i| dir.join(format!("in{i}.bam")).display().to_string())
            .collect();
        let mpileup = |tag: &str, threads: &str| {
            let out = dir
                .join(format!("out-{tag}-{threads}.vcf"))
                .display()
                .to_string();
            let mut args = vec![
                "vtools",
                "mpileup",
                "--no-version",
                "--no-reference",
                "--threads",
                threads,
                "-o",
                &out,
            ];
            args.extend(inputs.iter().map(String::as_str));
            let code = run_vtools(&args).expect("mpileup succeeds");
            assert_eq!(code, ExitCode::SUCCESS, "--threads {threads}");
        };
        if limit == "80" {
            mpileup("unlimited", "0");
        }
        let soft: libc::rlim_t = limit.parse().expect("a number of files");
        let mut lim = libc::rlimit {
            rlim_cur: 0,
            rlim_max: 0,
        };
        // SAFETY: getrlimit only writes the limits into the struct it is given.
        let got = unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut lim) };
        assert_eq!(got, 0, "getrlimit: {}", std::io::Error::last_os_error());
        lim.rlim_cur = soft.min(lim.rlim_max);
        // SAFETY: setrlimit only reads the struct it is given.
        let set = unsafe { libc::setrlimit(libc::RLIMIT_NOFILE, &lim) };
        assert_eq!(set, 0, "setrlimit: {}", std::io::Error::last_os_error());
        if limit == "30" {
            mpileup(&limit, "0");
        }
        mpileup(&limit, "4");
    }

    /// Where the helper test below finds its inputs and writes its outputs; set only in the
    /// child process that runs it.
    const TIGHT_NOFILE_DIR: &str = "VTOOLS_TEST_TIGHT_NOFILE_DIR";

    /// The runs of the helper test below: a tag, `--threads`, and the regions.
    const TIGHT_RUNS: [(&str, &str, Option<&str>); 4] = [
        ("r-0", "0", Some("c1")),
        ("r-4", "4", Some("c1")),
        ("all-0", "0", None),
        ("all-4", "4", None),
    ];

    #[test]
    fn many_inputs_pile_up_with_the_open_files_bcftools_needs() {
        let dir = std::env::temp_dir().join(format!("vtools-mpileup-tight-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create the test directory");
        for i in 0..8 {
            write_indexed_bam(&dir.join(format!("in{i}.bam")));
        }
        let exe = std::env::current_exe().expect("the test binary's path");
        let status = std::process::Command::new(&exe)
            .args([
                "commands::mpileup::tests::piles_up_with_the_open_files_bcftools_needs",
                "--exact",
                "--ignored",
                "--quiet",
            ])
            .env(TIGHT_NOFILE_DIR, &dir)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .expect("run the test binary");
        let read = |name: &str| std::fs::read(dir.join(name)).unwrap_or_default();
        let outputs: Vec<_> = TIGHT_RUNS
            .iter()
            .map(|(tag, ..)| {
                let want = read(&format!("want-{tag}.vcf"));
                (tag, want.len(), read(&format!("got-{tag}.vcf")) == want)
            })
            .collect();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(status.success(), "mpileup failed: {status}");
        for (tag, len, same) in outputs {
            assert!(len > 1000, "the unlimited run {tag} writes records");
            assert!(same, "run {tag} differs from the unlimited one");
        }
    }

    /// Run by the test above in a child process: piles up 8 indexed BAM files, with and
    /// without regions and threads, with no more open files than bcftools needs (the 8
    /// inputs, and the output or, while the inputs are opened, an index).
    #[test]
    #[ignore = "run in a child process by many_inputs_pile_up_with_the_open_files_bcftools_needs"]
    fn piles_up_with_the_open_files_bcftools_needs() {
        let Some(dir) = std::env::var_os(TIGHT_NOFILE_DIR) else {
            return;
        };
        let dir = PathBuf::from(dir);
        let inputs: Vec<String> = (0..8)
            .map(|i| dir.join(format!("in{i}.bam")).display().to_string())
            .collect();
        let mpileup = |out: &str, threads: &str, region: Option<&str>| {
            let out = dir.join(out).display().to_string();
            let mut args = vec![
                "vtools",
                "mpileup",
                "--no-version",
                "--no-reference",
                "--threads",
                threads,
                "-o",
                &out,
            ];
            if let Some(region) = region {
                args.extend(["-r", region]);
            }
            args.extend(inputs.iter().map(String::as_str));
            let code = run_vtools(&args).expect("mpileup succeeds");
            assert_eq!(code, ExitCode::SUCCESS, "--threads {threads} -r {region:?}");
        };
        for (tag, _, region) in TIGHT_RUNS {
            mpileup(&format!("want-{tag}.vcf"), "0", region);
        }
        // Listing the open descriptors opens one more.
        let open = std::fs::read_dir("/dev/fd").map_or(0, Iterator::count) - 1;
        let mut lim = libc::rlimit {
            rlim_cur: 0,
            rlim_max: 0,
        };
        // SAFETY: getrlimit only writes the limits into the struct it is given.
        let got = unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut lim) };
        assert_eq!(got, 0, "getrlimit: {}", std::io::Error::last_os_error());
        lim.rlim_cur = (open as libc::rlim_t + 9).min(lim.rlim_max);
        // SAFETY: setrlimit only reads the struct it is given.
        let set = unsafe { libc::setrlimit(libc::RLIMIT_NOFILE, &lim) };
        assert_eq!(set, 0, "setrlimit: {}", std::io::Error::last_os_error());
        for (tag, threads, region) in TIGHT_RUNS {
            mpileup(&format!("got-{tag}.vcf"), threads, region);
        }
    }

    #[test]
    fn parses_flags() {
        assert_eq!(str2flag("UNMAP,DUP"), 0x404);
        assert_eq!(str2flag("1024"), 1024);
        assert_eq!(str2flag("0x10"), 16);
        assert_eq!(str2flag("010"), 8);
        assert_eq!(str2flag("12abc"), 12);
        assert_eq!(str2flag("unmap,Dup"), 0x404);
        assert_eq!(str2flag("UNMAP,"), -1);
        assert_eq!(str2flag(""), 0);
        assert_eq!(str2flag("NOSUCH"), -1);
    }
}
