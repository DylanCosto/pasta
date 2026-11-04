//! The text report of bcftools stats (`print_header` and `print_stats` in `vcfstats.c`),
//! with C's `printf` number formatting.

use std::io::{self, Write};

use super::{
    bins::AfBins,
    collect::{Config, M_INDEL, NAF_HWE, Stats, VAF_BINS},
    compare::{GtCompare, GtCounts},
    indel_context::REP_LEN,
};

/// What the report describes besides the statistics.
#[derive(Debug)]
pub struct ReportInfo<'a> {
    /// The input names, `-` for standard input.
    pub files: &'a [String],
    /// The number of samples in each input's header.
    pub nsamples: &'a [usize],
    /// The names of the samples of `-s/-S`.
    pub sample_names: &'a [Vec<u8>],
    /// `-I/--split-by-ID`.
    pub split_by_id: bool,
    /// `-v`: per-site discordance output (two files).
    pub verbose_sites: bool,
    /// `-E/--exons` was given.
    pub exons: bool,
    /// `-F` loaded a reference.
    pub indel_ctx: bool,
}

/// How glibc's `printf` writes a value that is not finite (`-nan` when the sign bit is set, as
/// it is for x86's default NaN from `0/0` or `0*inf`), or `None` for a finite value.
fn c_nonfinite(x: f64) -> Option<&'static str> {
    match (x.is_nan(), x.is_infinite(), x.is_sign_negative()) {
        (true, _, true) => Some("-nan"),
        (true, _, false) => Some("nan"),
        (_, true, true) => Some("-inf"),
        (_, true, false) => Some("inf"),
        _ => None,
    }
}

/// C's `%.0f`.
fn fmt_f0(x: f64) -> String {
    match c_nonfinite(x) {
        Some(s) => s.to_owned(),
        None => format!("{x:.0}"),
    }
}

/// C's `%e`: one digit, six decimals, and an exponent of at least two digits.
pub fn fmt_e(x: f64) -> String {
    if let Some(s) = c_nonfinite(x) {
        return s.to_owned();
    }
    let s = format!("{x:.6e}");
    match s.split_once('e') {
        Some((mantissa, exp)) => {
            let (sign, digits) = match exp.strip_prefix('-') {
                Some(d) => ('-', d),
                None => ('+', exp),
            };
            format!("{mantissa}e{sign}{digits:0>2}")
        }
        None => s,
    }
}

/// The name of an input in the report.
fn display_name(f: &str) -> &str {
    if f == "-" { "<STDIN>" } else { f }
}

/// Prints the header of the report (`print_header`).
pub fn print_header(out: &mut impl Write, info: &ReportInfo<'_>) -> io::Result<()> {
    writeln!(
        out,
        "# This file was produced by vtools stats ({}+pasta-{}) and can be plotted using plot-vcfstats.",
        env!("CARGO_PKG_VERSION"),
        env!("CARGO_PKG_VERSION")
    )?;
    write!(out, "# The command line was:\tvtools stats ")?;
    for a in std::env::args().skip(2) {
        write!(out, " {a}")?;
    }
    write!(out, "\n#\n")?;
    writeln!(
        out,
        "# Definition of sets:\n# ID\t[2]id\t[3]tab-separated file names"
    )?;
    if info.files.len() == 1 {
        let fname = display_name(&info.files[0]);
        if info.split_by_id {
            writeln!(
                out,
                "ID\t0\t{fname}:known (sites with ID different from \".\")"
            )?;
            writeln!(out, "ID\t1\t{fname}:novel (sites where ID column is \".\")")?;
        } else {
            writeln!(out, "ID\t0\t{fname}")?;
        }
    } else {
        let f0 = display_name(&info.files[0]);
        let f1 = display_name(&info.files[1]);
        writeln!(out, "ID\t0\t{f0}")?;
        writeln!(out, "ID\t1\t{f1}")?;
        writeln!(out, "ID\t2\t{f0}\t{f1}")?;
        if info.verbose_sites {
            write!(
                out,
                "# Verbose per-site discordance output.\n\
                 # PSD\t[2]CHROM\t[3]POS\t[4]Number of matches\t[5]Number of mismatches\t[6]NRD\n"
            )?;
            write!(
                out,
                "# Verbose per-site and per-sample output. Genotype codes: 0:HomRefRef, 1:HomAltAlt, 2:HetAltRef, 3:HetAltAlt, 4:haploidRef, 5:haploidAlt\n\
                 # DBG\t[2]CHROM\t[3]POS\t[4]Sample\t[5]GT in {f0}\t[6]GT in {f1}\n"
            )?;
        }
    }
    Ok(())
}

/// The frequency printed for allele-frequency bin `i`.
fn af_value(cfg: &Config, bins: Option<&AfBins>, i: usize) -> f64 {
    match bins {
        Some(b) => f64::from(b.value(i) + b.value(i - 1)) * 0.5,
        None => (i as f64 - 1.0) / (cfg.m_af as f64 - 1.0),
    }
}

/// Prints the statistics (`print_stats`). Moves the singleton counts into the first
/// frequency bin, as bcftools does while printing.
pub fn print_stats(
    out: &mut impl Write,
    cfg: &Config,
    info: &ReportInfo<'_>,
    sets: &mut [Stats],
    gtcmp: Option<&mut GtCompare>,
) -> io::Result<()> {
    print_summary(out, info, sets)?;
    print_tstv(out, sets)?;
    if info.exons {
        print_frameshifts(out, sets)?;
    }
    if info.indel_ctx {
        print_indel_context(out, sets)?;
    }
    print_singletons(out, sets)?;
    let mut gtcmp = gtcmp;
    if let Some(g) = gtcmp.as_deref_mut() {
        g.fold_singletons();
    }
    print_af(out, cfg, sets)?;
    print_qual(out, sets)?;
    print_user(out, cfg, sets)?;
    print_idd(out, sets)?;
    print_subst(out, sets)?;
    if info.files.len() > 1
        && !info.sample_names.is_empty()
        && let Some(g) = gtcmp
    {
        print_concordance(out, cfg, info, g)?;
    }
    print_depth(out, sets)?;
    if !info.sample_names.is_empty() {
        print_per_sample(out, cfg, info, sets)?;
    }
    if sets[0].smpl_vaf.is_some() {
        print_vaf(out, info, sets)?;
    }
    Ok(())
}

fn print_summary(out: &mut impl Write, info: &ReportInfo<'_>, sets: &[Stats]) -> io::Result<()> {
    out.write_all(
        b"# SN, Summary numbers:\n\
#   number of records   .. number of data rows in the VCF\n\
#   number of no-ALTs   .. reference-only sites, ALT is either \".\" or identical to REF\n\
#   number of SNPs      .. number of rows with a SNP\n\
#   number of MNPs      .. number of rows with a MNP, such as CC>TT\n\
#   number of indels    .. number of rows with an indel\n\
#   number of others    .. number of rows with other type, for example a symbolic allele or\n\
#                          a complex substitution, such as ACT>TCGA\n\
#   number of multiallelic sites     .. number of rows with multiple alternate alleles\n\
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs\n\
# \n\
#   Note that rows containing multiple types will be counted multiple times, in each\n\
#   counter. For example, a row with a SNP and an indel increments both the SNP and\n\
#   the indel counter.\n\
# \n\
# SN\t[2]id\t[3]key\t[4]value\n",
    )?;
    for (id, n) in info.nsamples.iter().enumerate() {
        writeln!(out, "SN\t{id}\tnumber of samples:\t{n}")?;
    }
    for (id, s) in sets.iter().enumerate() {
        writeln!(out, "SN\t{id}\tnumber of records:\t{}", s.n_records)?;
        writeln!(out, "SN\t{id}\tnumber of no-ALTs:\t{}", s.n_noalts)?;
        writeln!(out, "SN\t{id}\tnumber of SNPs:\t{}", s.n_snps)?;
        writeln!(out, "SN\t{id}\tnumber of MNPs:\t{}", s.n_mnps)?;
        writeln!(out, "SN\t{id}\tnumber of indels:\t{}", s.n_indels)?;
        writeln!(out, "SN\t{id}\tnumber of others:\t{}", s.n_others)?;
        writeln!(out, "SN\t{id}\tnumber of multiallelic sites:\t{}", s.n_mals)?;
        writeln!(
            out,
            "SN\t{id}\tnumber of multiallelic SNP sites:\t{}",
            s.n_snp_mals
        )?;
    }
    Ok(())
}

/// `x / y` as a C float, or 0 when `y` is 0.
fn ratio(x: u64, y: u64) -> f64 {
    if y == 0 {
        0.0
    } else {
        f64::from(x as f32 / y as f32)
    }
}

fn print_tstv(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(
        b"# TSTV, transitions/transversions\n\
#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)\n\
#   - transversions, see https://en.wikipedia.org/wiki/Transversion\n\
# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)\n",
    )?;
    for (id, s) in sets.iter().enumerate() {
        let ts: u64 = s.af_ts.iter().sum();
        let tv: u64 = s.af_tv.iter().sum();
        writeln!(
            out,
            "TSTV\t{id}\t{ts}\t{tv}\t{:.2}\t{}\t{}\t{:.2}",
            ratio(ts, tv),
            s.ts_alt1,
            s.tv_alt1,
            ratio(s.ts_alt1, s.tv_alt1)
        )?;
    }
    Ok(())
}

fn print_frameshifts(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(b"# FS, Indel frameshifts:\n# FS\t[2]id\t[3]in-frame\t[4]out-frame\t[5]not applicable\t[6]out/(in+out) ratio\t[7]in-frame (1st ALT)\t[8]out-frame (1st ALT)\t[9]not applicable (1st ALT)\t[10]out/(in+out) ratio (1st ALT)\n")?;
    let fs_ratio = |out: i64, inn: i64| -> f64 {
        if out == 0 {
            0.0
        } else {
            f64::from(out as f32 / (inn + out) as f32)
        }
    };
    for (id, s) in sets.iter().enumerate() {
        writeln!(
            out,
            "FS\t{id}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{:.2}",
            s.in_frame,
            s.out_frame,
            s.na_frame,
            fs_ratio(s.out_frame, s.in_frame),
            s.in_frame_alt1,
            s.out_frame_alt1,
            s.na_frame_alt1,
            fs_ratio(s.out_frame_alt1, s.in_frame_alt1)
        )?;
    }
    Ok(())
}

fn print_indel_context(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(
        b"# ICS, Indel context:\n\
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]\n\
# ICS\t[2]id\t[3]repeat-consistent\t[4]repeat-inconsistent\t[5]not applicable\t[6]c/(c+i) ratio\n",
    )?;
    for (id, s) in sets.iter().enumerate() {
        let (mut nc, mut ni) = (0u64, 0u64);
        for r in &s.n_repeat {
            nc += r[0] + r[2];
            ni += r[1] + r[3];
        }
        writeln!(
            out,
            "ICS\t{id}\t{nc}\t{ni}\t{}\t{:.4}",
            s.n_repeat_na,
            ratio(nc, nc + ni)
        )?;
    }
    out.write_all(
        b"# ICL, Indel context by length:\n\
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]\n\
# ICL\t[2]id\t[3]length of repeat element\t[4]repeat-consistent deletions)\t[5]repeat-inconsistent deletions\t[6]consistent insertions\t[7]inconsistent insertions\t[8]c/(c+i) ratio\n",
    )?;
    for (id, s) in sets.iter().enumerate() {
        for i in 1..REP_LEN {
            let r = s.n_repeat[i];
            let nc = r[0] + r[2];
            let ni = r[1] + r[3];
            writeln!(
                out,
                "ICL\t{id}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                i + 1,
                r[0],
                r[1],
                r[2],
                r[3],
                ratio(nc, nc + ni)
            )?;
        }
    }
    Ok(())
}

fn print_singletons(out: &mut impl Write, sets: &mut [Stats]) -> io::Result<()> {
    out.write_all(
        b"# SiS, Singleton stats:\n\
#   - allele count, i.e. the number of singleton genotypes (AC=1)\n\
#   - number of transitions, see above\n\
#   - number of transversions, see above\n\
#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]\n\
# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n",
    )?;
    for (id, s) in sets.iter_mut().enumerate() {
        let r = [s.af_repeats[0][0], s.af_repeats[1][0], s.af_repeats[2][0]];
        writeln!(
            out,
            "SiS\t{id}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.af_snps[0],
            s.af_ts[0],
            s.af_tv[0],
            r[0] + r[1] + r[2],
            r[0],
            r[1],
            r[2]
        )?;
        // The singletons go into the first frequency bin from here on.
        s.af_snps[1] += s.af_snps[0];
        s.af_ts[1] += s.af_ts[0];
        s.af_tv[1] += s.af_tv[0];
        for rep in &mut s.af_repeats {
            rep[1] += rep[0];
        }
    }
    Ok(())
}

fn print_af(out: &mut impl Write, cfg: &Config, sets: &[Stats]) -> io::Result<()> {
    out.write_all(b"# AF, Stats by non-reference allele frequency:\n# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n")?;
    for (id, s) in sets.iter().enumerate() {
        for i in 1..cfg.m_af {
            let r = [s.af_repeats[0][i], s.af_repeats[1][i], s.af_repeats[2][i]];
            if s.af_snps[i] + s.af_ts[i] + s.af_tv[i] + r[0] + r[1] + r[2] == 0 {
                continue;
            }
            let af = af_value(cfg, cfg.af_bins.as_ref(), i);
            writeln!(
                out,
                "AF\t{id}\t{af:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.af_snps[i],
                s.af_ts[i],
                s.af_tv[i],
                r[0] + r[1] + r[2],
                r[0],
                r[1],
                r[2]
            )?;
        }
    }
    Ok(())
}

fn print_qual(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(b"# QUAL, Stats by quality\n# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels\n")?;
    let qval_of = |beg: u32| -> f32 {
        if beg > 0 {
            (0.1 * f64::from(beg - 1)) as f32
        } else {
            -1.0
        }
    };
    for (id, s) in sets.iter().enumerate() {
        let (n_ts, n_tv, n_in) = (s.qual_ts.len(), s.qual_tv.len(), s.qual_indels.len());
        let n_max = n_ts.max(n_tv).max(n_in);
        for i in 0..n_max {
            let (mut nts, mut ntv, mut nin) = (0u64, 0u64, 0u64);
            let mut qval = -1f32;
            if i < n_ts {
                let (n, beg) = s.qual_ts.get(i);
                nts = n;
                qval = qval_of(beg);
            }
            if i < n_tv {
                let (n, beg) = s.qual_tv.get(i);
                ntv = n;
                if qval == -1.0 {
                    qval = qval_of(beg);
                }
            }
            if i < n_in {
                let (n, beg) = s.qual_indels.get(i);
                nin = n;
                if qval == -1.0 {
                    qval = qval_of(beg);
                }
            }
            if nts + ntv + nin == 0 {
                continue;
            }
            write!(out, "QUAL\t{id}\t")?;
            if qval == -1.0 {
                write!(out, ".")?;
            } else {
                write!(out, "{:.1}", f64::from(qval))?;
            }
            writeln!(out, "\t{}\t{nts}\t{ntv}\t{nin}", nts + ntv)?;
        }
    }
    Ok(())
}

fn print_user(out: &mut impl Write, cfg: &Config, sets: &[Stats]) -> io::Result<()> {
    for (i, def) in cfg.user.iter().enumerate() {
        let (tag, idx, min, max) = (String::from_utf8_lossy(&def.tag), def.idx, def.min, def.max);
        writeln!(
            out,
            "# USR:{tag}/{idx}\t[2]id\t[3]{tag}/{idx}\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)"
        )?;
        for (id, s) in sets.iter().enumerate() {
            let usr = &s.usr[i];
            let nbins = usr.vals_ts.len();
            for j in 0..nbins {
                let (ts, tv) = (usr.vals_ts[j], usr.vals_tv[j]);
                if ts + tv == 0 {
                    continue;
                }
                // With one bin this is 0/0, a NaN.
                let val = min + (max - min) * j as f32 / (nbins as f32 - 1.0);
                let val = f64::from(val);
                let v = if usr.is_float {
                    fmt_e(val)
                } else {
                    fmt_f0(val)
                };
                writeln!(out, "USR:{tag}/{idx}\t{id}\t{v}\t{}\t{ts}\t{tv}", ts + tv)?;
            }
        }
    }
    Ok(())
}

fn print_idd(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(b"# IDD, InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF\n")?;
    let vaf = |out: &mut dyn Write, s: &Stats, bin: usize| -> io::Result<()> {
        match s.nvaf.get(bin) {
            Some(&n) if n > 0 => writeln!(out, "{n}\t{:.2}", s.dvaf[bin] / n as f64),
            _ => writeln!(out, "0\t."),
        }
    };
    for (id, s) in sets.iter().enumerate() {
        for i in (0..M_INDEL).rev() {
            if s.deletions[i] == 0 {
                continue;
            }
            write!(out, "IDD\t{id}\t{}\t{}\t", -(i as i64) - 1, s.deletions[i])?;
            vaf(out, s, M_INDEL - i - 1)?;
        }
        for i in 0..M_INDEL {
            if s.insertions[i] == 0 {
                continue;
            }
            write!(out, "IDD\t{id}\t{}\t{}\t", i + 1, s.insertions[i])?;
            vaf(out, s, M_INDEL + i + 1)?;
        }
    }
    Ok(())
}

fn print_subst(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(b"# ST, Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n")?;
    const ACGT: &[u8; 4] = b"ACGT";
    for (id, s) in sets.iter().enumerate() {
        for t in 0..15 {
            if t >> 2 == (t & 3) {
                continue;
            }
            writeln!(
                out,
                "ST\t{id}\t{}>{}\t{}",
                char::from(ACGT[t >> 2]),
                char::from(ACGT[t & 3]),
                s.subst[t]
            )?;
        }
    }
    Ok(())
}

fn print_depth(out: &mut impl Write, sets: &[Stats]) -> io::Result<()> {
    out.write_all(
        b"# DP, depth:\n\
#   - set id, see above\n\
#   - the depth bin, corresponds to the depth (unless --depth was given)\n\
#   - number of genotypes with this depth (zero unless -s/-S was given)\n\
#   - fraction of genotypes with this depth (zero unless -s/-S was given)\n\
#   - number of sites with this depth\n\
#   - fraction of sites with this depth\n\
# DP, Depth distribution\n# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)\n",
    )?;
    for (id, s) in sets.iter().enumerate() {
        let sum: u64 = s.dp.vals.iter().sum();
        let sum_sites: u64 = s.dp_sites.vals.iter().sum();
        let m = s.dp.vals.len();
        let pct = |v: u64, sum: u64| {
            if sum == 0 {
                0.0
            } else {
                v as f64 * 100.0 / sum as f64
            }
        };
        for i in 0..m {
            let (g, st) = (s.dp.vals[i], s.dp_sites.vals[i]);
            if g == 0 && st == 0 {
                continue;
            }
            write!(out, "DP\t{id}\t")?;
            if i == 0 {
                write!(out, "<{}", s.dp.bins.min)?;
            } else if i + 1 == m {
                write!(out, ">{}", s.dp.bins.max)?;
            } else {
                write!(out, "{}", s.dp.bin_start(i))?;
            }
            writeln!(
                out,
                "\t{g}\t{:.6}\t{st}\t{:.6}",
                pct(g, sum),
                pct(st, sum_sites)
            )?;
        }
    }
    Ok(())
}

fn print_per_sample(
    out: &mut impl Write,
    cfg: &Config,
    info: &ReportInfo<'_>,
    sets: &mut [Stats],
) -> io::Result<()> {
    out.write_all(b"# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.\n")?;
    out.write_all(b"# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing\n")?;
    for (id, s) in sets.iter().enumerate() {
        for (c, name) in s.smpl.iter().zip(info.sample_names) {
            let dp = if c.ndp == 0 {
                0.0
            } else {
                c.dp as f32 / c.ndp as f32
            };
            out.write_all(format!("PSC\t{id}\t").as_bytes())?;
            out.write_all(name)?;
            writeln!(
                out,
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
                c.hom_rr,
                c.hom_aa,
                c.hets,
                c.ts,
                c.tv,
                c.indels,
                f64::from(dp),
                c.sngl,
                c.hap_ref,
                c.hap_alt,
                c.missing
            )?;
        }
    }
    out.write_all(b"# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.\n")?;
    out.write_all(b"# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms\n")?;
    for (id, s) in sets.iter().enumerate() {
        for (c, name) in s.smpl.iter().zip(info.sample_names) {
            let [na, inn, outf] = if cfg.exons { c.frm_shifts } else { [0; 3] };
            let r = if inn + outf == 0 {
                0.0
            } else {
                outf as f64 / (inn + outf) as f64
            };
            out.write_all(format!("PSI\t{id}\t").as_bytes())?;
            out.write_all(name)?;
            writeln!(
                out,
                "\t{inn}\t{outf}\t{na}\t{r:.2}\t{}\t{}\t{}\t{}",
                c.ins_hets, c.del_hets, c.ins_homs, c.del_homs
            )?;
        }
    }

    out.write_all(b"# HWE\n# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile\n")?;
    for (id, s) in sets.iter_mut().enumerate() {
        for i in 0..NAF_HWE {
            s.af_hwe[i + NAF_HWE] += s.af_hwe[i];
        }
        for i in 1..cfg.m_af {
            let row = &s.af_hwe[i * NAF_HWE..(i + 1) * NAF_HWE];
            let sum_tot: u64 = row.iter().sum();
            if sum_tot == 0 {
                continue;
            }
            let af = af_value(cfg, cfg.af_bins.as_ref(), i);
            write!(out, "HWE\t{id}\t{af:.6}\t{sum_tot}")?;
            let mut nprn = 3;
            let mut sum_tmp = 0u64;
            let q = |j: usize| f64::from(j as f32 / NAF_HWE as f32);
            for (j, &v) in row.iter().enumerate() {
                sum_tmp += v;
                let frac = sum_tmp as f32 / sum_tot as f32;
                if frac >= 0.75 {
                    while nprn > 0 {
                        write!(out, "\t{:.6}", q(j))?;
                        nprn -= 1;
                    }
                    break;
                }
                if frac >= 0.5 {
                    while nprn > 1 {
                        write!(out, "\t{:.6}", q(j))?;
                        nprn -= 1;
                    }
                    continue;
                }
                if frac >= 0.25 {
                    while nprn > 2 {
                        write!(out, "\t{:.6}", q(j))?;
                        nprn -= 1;
                    }
                }
            }
            writeln!(out)?;
        }
    }
    Ok(())
}

fn print_vaf(out: &mut impl Write, info: &ReportInfo<'_>, sets: &[Stats]) -> io::Result<()> {
    out.write_all(b"# VAF, Variant Allele Frequency determined as fraction of alternate reads in FORMAT/AD\n# VAF\t[2]id\t[3]sample\t[4]SNV VAF distribution\t[5]indel VAF distribution\n")?;
    for (id, s) in sets.iter().enumerate() {
        let Some(vafs) = &s.smpl_vaf else { continue };
        for (v, name) in vafs.iter().zip(info.sample_names) {
            out.write_all(format!("VAF\t{id}\t").as_bytes())?;
            out.write_all(name)?;
            out.write_all(b"\t")?;
            for j in 0..VAF_BINS {
                write!(out, "{}{}", if j > 0 { "," } else { "" }, v.snv[j])?;
            }
            out.write_all(b"\t")?;
            for j in 0..VAF_BINS {
                write!(out, "{}{}", if j > 0 { "," } else { "" }, v.indel[j])?;
            }
            out.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Pearson's r² of the dosages, as bcftools computes it.
fn r2(c: &GtCounts) -> f64 {
    if c.n == 0.0 {
        return 0.0;
    }
    let mut r = c.yx - c.x * c.y / c.n;
    r /= ((c.xx - c.x * c.x / c.n) * (c.yy - c.y * c.y / c.n)).sqrt();
    r * r
}

// Genotype classes (bcftools' `type2stats`): RR, RA, AA hom, AA het, missing.
const RR: usize = 0;
const RA: usize = 1;
const AA: usize = 2;
const AAHET: usize = 3;

/// Prints the genotype concordance of two files (the `nreaders > 1` part of `print_stats`).
fn print_concordance(
    out: &mut impl Write,
    cfg: &Config,
    info: &ReportInfo<'_>,
    g: &GtCompare,
) -> io::Result<()> {
    writeln!(
        out,
        "SN\t2\tnumber of samples:\t{}",
        info.sample_names.len()
    )?;
    print_concordance_by_af(out, cfg, g)?;
    print_concordance_by_sample(out, info, g)?;
    print_concordance_table(out, info, g)
}

/// Prints the GCsAF/GCiAF concordance by allele frequency and the NRDs/NRDi discordance.
fn print_concordance_by_af(out: &mut impl Write, cfg: &Config, g: &GtCompare) -> io::Result<()> {
    for x in 0..2 {
        let (stats, c) = if x == 0 {
            out.write_all(b"# GCsAF, Genotype concordance by non-reference allele frequency (SNPs)\n# GCsAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of genotypes\n")?;
            (&g.af_snps, 's')
        } else {
            out.write_all(b"# GCiAF, Genotype concordance by non-reference allele frequency (indels)\n# GCiAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of genotypes\n")?;
            (&g.af_indels, 'i')
        };
        let mut nrd_m = [0u64; 4];
        let mut nrd_mm = [0u64; 4];
        for (i, s) in stats.iter().enumerate() {
            let mut n = 0u64;
            let mut m = [0u64; 4];
            let mut mm = [0u64; 4];
            for j in 0..4 {
                for k in 0..4 {
                    let v = s.gt2gt[j][k];
                    n += v;
                    if j == k {
                        nrd_m[j] += v;
                        m[j] += v;
                    } else {
                        nrd_mm[j] += v;
                        mm[j] += v;
                    }
                }
            }
            if i == 0 || n == 0 {
                continue;
            }
            let r = r2(s);
            let af = af_value(cfg, cfg.af_bins.as_ref(), i);
            write!(out, "GC{c}AF\t2\t{af:.6}")?;
            write!(out, "\t{}\t{}\t{}", m[RR], m[RA], m[AA])?;
            write!(out, "\t{}\t{}\t{}", mm[RR], mm[RA], mm[AA])?;
            if s.n != 0.0 && !r.is_nan() {
                write!(out, "\t{r:.6}")?;
            } else {
                write!(out, "\t0")?;
            }
            writeln!(out, "\t{:.0}", s.n)?;
        }
        if x == 0 {
            out.write_all(
                b"# NRD and discordance is calculated as follows:\n\
#   m .. number of matches\n\
#   x .. number of mismatches\n\
#   NRD = 100 * (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)\n\
#   RR discordance = 100 * xRR / (xRR + mRR)\n\
#   RA discordance = 100 * xRA / (xRA + mRA)\n\
#   AA discordance = 100 * xAA / (xAA + mAA)\n\
# Non-Reference Discordance (NRD), SNPs\n# NRDs\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n",
            )?;
        } else {
            out.write_all(b"# Non-Reference Discordance (NRD), indels\n# NRDi\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n")?;
        }
        let m = nrd_m[RA] + nrd_m[AA] + nrd_m[AAHET];
        let mm = nrd_mm[RR] + nrd_mm[RA] + nrd_mm[AA] + nrd_mm[AAHET];
        let pct = |a: u64, b: u64| {
            if a + b == 0 {
                0.0
            } else {
                a as f64 * 100.0 / (a + b) as f64
            }
        };
        writeln!(
            out,
            "NRD{c}\t2\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            pct(mm, m),
            pct(nrd_mm[RR], nrd_m[RR]),
            pct(nrd_mm[RA], nrd_m[RA]),
            pct(nrd_mm[AA], nrd_m[AA])
        )?;
    }
    Ok(())
}

/// Prints the GCsS/GCiS concordance by sample.
fn print_concordance_by_sample(
    out: &mut impl Write,
    info: &ReportInfo<'_>,
    g: &GtCompare,
) -> io::Result<()> {
    for x in 0..2 {
        let (stats, c) = if x == 0 {
            out.write_all(b"# GCsS, Genotype concordance by sample (SNPs)\n# GCsS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n")?;
            (&g.smpl_snps, 's')
        } else {
            out.write_all(b"# GCiS, Genotype concordance by sample (indels)\n# GCiS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n")?;
            (&g.smpl_indels, 'i')
        };
        for (s, name) in stats.iter().zip(info.sample_names) {
            let m = s.gt2gt[RA][RA] + s.gt2gt[AA][AA];
            let mut mm = 0u64;
            for j in 0..3 {
                for k in 0..3 {
                    if j != k {
                        mm += s.gt2gt[j][k];
                    }
                }
            }
            let r = r2(s);
            let nrd = if m + mm == 0 {
                0.0
            } else {
                mm as f64 * 100.0 / (m + mm) as f64
            };
            out.write_all(format!("GC{c}S\t2\t").as_bytes())?;
            out.write_all(name)?;
            write!(out, "\t{nrd:.3}")?;
            write!(
                out,
                "\t{}\t{}\t{}",
                s.gt2gt[RR][RR], s.gt2gt[RA][RA], s.gt2gt[AA][AA]
            )?;
            write!(
                out,
                "\t{}\t{}\t{}",
                s.gt2gt[RR][RA] + s.gt2gt[RR][AA],
                s.gt2gt[RA][RR] + s.gt2gt[RA][AA],
                s.gt2gt[AA][RR] + s.gt2gt[AA][RA]
            )?;
            if s.n != 0.0 && !r.is_nan() {
                writeln!(out, "\t{r:.6}")?;
            } else {
                writeln!(out, "\t0")?;
            }
        }
    }
    Ok(())
}

/// Prints the GCTs/GCTi tables of genotype classes in one file against the other, by sample.
fn print_concordance_table(
    out: &mut impl Write,
    info: &ReportInfo<'_>,
    g: &GtCompare,
) -> io::Result<()> {
    const CLASSES: [&str; 5] = ["RR Hom", "RA Het", "AA Hom", "AA Het", "missing"];
    for x in 0..2 {
        let (stats, c) = if x == 0 {
            out.write_all(b"# GCTs, Genotype concordance table (SNPs)\n# GCTs")?;
            (&g.smpl_snps, 's')
        } else {
            out.write_all(b"# GCTi, Genotype concordance table (indels)\n# GCTi")?;
            (&g.smpl_indels, 'i')
        };
        let mut col = 2;
        write!(out, "\t[{col}]sample")?;
        for a in CLASSES {
            for b in CLASSES {
                col += 1;
                write!(out, "\t[{col}]{a} -> {b}")?;
            }
        }
        writeln!(out)?;
        for (s, name) in stats.iter().zip(info.sample_names) {
            out.write_all(format!("GCT{c}\t").as_bytes())?;
            out.write_all(name)?;
            for row in &s.gt2gt {
                for v in row {
                    write!(out, "\t{v}")?;
                }
            }
            writeln!(out)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e_format_matches_c() {
        assert_eq!(fmt_e(0.000_123_456_789), "1.234568e-04");
        assert_eq!(fmt_e(0.0), "0.000000e+00");
        assert_eq!(fmt_e(12345.0), "1.234500e+04");
        assert_eq!(fmt_e(1e100), "1.000000e+100");
        assert_eq!(fmt_e(-f64::NAN), "-nan");
        assert_eq!(fmt_e(f64::INFINITY), "inf");
        assert_eq!(fmt_f0(-f64::NAN), "-nan");
        assert_eq!(fmt_f0(2.5), "2");
    }
}
