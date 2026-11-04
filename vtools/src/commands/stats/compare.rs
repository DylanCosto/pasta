//! Genotype concordance between two files (the `matched==3` part of `do_sample_stats` in
//! `vcfstats.c`), for bcftools stats with two inputs and `-s/-S`.
//!
//! The pairing of records between the files is the multi-file synced reader's job; this module
//! takes a pair of records at the same site and counts how their genotypes agree.

use std::io::{self, Write};

use pasta::variant::{
    FormatField, Header, HeaderRecordType, Record, VariantType,
    utils::{self, GenotypeType},
};

/// Finds a FORMAT field of an unpacked record by tag (`bcf_get_fmt`).
fn find_format<'r>(header: &Header, rec: &'r Record, tag: &[u8]) -> Option<&'r FormatField> {
    let id = header
        .id(tag)
        .filter(|&id| header.id_exists(HeaderRecordType::Format, id))?;
    rec.format_fields().iter().find(|f| f.id() == id as i32)
}

/// HTSlib's numeric `GT_*` code of a genotype type.
fn gt_code(gt: GenotypeType) -> u8 {
    match gt {
        GenotypeType::HomRef => 0,
        GenotypeType::HomAlt => 1,
        GenotypeType::HetRefAlt => 2,
        GenotypeType::HetAltAlt => 3,
        GenotypeType::HaploidRef => 4,
        GenotypeType::HaploidAlt => 5,
        GenotypeType::Unknown => 6,
    }
}

/// The concordance class of a genotype (bcftools' `type2stats`): RR, RA, AA, AA het,
/// missing; haploid genotypes count as homozygous.
fn class(gt: GenotypeType) -> usize {
    match gt {
        GenotypeType::HomRef | GenotypeType::HaploidRef => 0,
        GenotypeType::HetRefAlt => 1,
        GenotypeType::HomAlt | GenotypeType::HaploidAlt => 2,
        GenotypeType::HetAltAlt => 3,
        GenotypeType::Unknown => 4,
    }
}

/// The alternate-allele dosage of a genotype (bcftools' `type2dosage`).
fn dosage(gt: GenotypeType) -> f32 {
    match gt {
        GenotypeType::HetRefAlt | GenotypeType::HaploidAlt => 1.0,
        GenotypeType::HomAlt | GenotypeType::HetAltAlt => 2.0,
        _ => 0.0,
    }
}

/// Whether a genotype is haploid (bcftools' `type2ploidy` of -1).
fn is_haploid(gt: GenotypeType) -> bool {
    matches!(gt, GenotypeType::HaploidRef | GenotypeType::HaploidAlt)
}

/// Concordance counts and dosage sums (bcftools' `gtcmp_t`).
#[derive(Clone, Copy, Debug, Default)]
pub struct GtCounts {
    /// Genotype class in the first file by class in the second.
    pub gt2gt: [[u64; 5]; 5],
    pub y: f64,
    pub yy: f64,
    pub x: f64,
    pub xx: f64,
    pub yx: f64,
    pub n: f64,
}

/// Genotype concordance by allele frequency and by sample, for SNPs and indels.
#[derive(Clone, Debug)]
pub struct GtCompare {
    pub af_snps: Vec<GtCounts>,
    pub af_indels: Vec<GtCounts>,
    pub smpl_snps: Vec<GtCounts>,
    pub smpl_indels: Vec<GtCounts>,
}

/// One side of a pair of records at the same site; the record's FORMAT fields are unpacked.
#[derive(Debug)]
pub struct Side<'a> {
    pub header: &'a Header,
    pub rec: &'a Record,
    /// The header columns of the compared samples.
    pub samples: &'a [usize],
}

impl GtCompare {
    /// Creates empty counts for `m_af` frequency bins and `n_smpl` samples.
    pub fn new(m_af: usize, n_smpl: usize) -> Self {
        Self {
            af_snps: vec![GtCounts::default(); m_af],
            af_indels: vec![GtCounts::default(); m_af],
            smpl_snps: vec![GtCounts::default(); n_smpl],
            smpl_indels: vec![GtCounts::default(); n_smpl],
        }
    }

    /// Moves the singleton dosage sums into the first frequency bin, as `print_stats` does
    /// (without `x` or the genotype table, as bcftools).
    pub fn fold_singletons(&mut self) {
        for v in [&mut self.af_snps, &mut self.af_indels] {
            if v.len() < 2 {
                continue;
            }
            let s = v[0];
            let d = &mut v[1];
            d.y += s.y;
            d.yy += s.yy;
            d.xx += s.xx;
            d.yx += s.yx;
            d.n += s.n;
        }
    }

    /// Compares the genotypes of a record pair. `iaf1` is the frequency bin of the first
    /// file's first ALT allele and `line_type` its record's variant types; `sample_names`
    /// and `verbose` give the `-v` per-site output.
    ///
    /// # Errors
    ///
    /// Returns write errors of the verbose output.
    pub fn add(
        &mut self,
        a: &Side<'_>,
        b: &Side<'_>,
        iaf1: usize,
        line_type: VariantType,
        verbose: Option<(&mut dyn Write, &[Vec<u8>])>,
    ) -> io::Result<()> {
        let (Some(fmt0), Some(fmt1)) = (
            find_format(a.header, a.rec, b"GT"),
            find_format(b.header, b.rec, b"GT"),
        ) else {
            return Ok(());
        };
        let snp = line_type.contains(VariantType::SNP);
        let (af_stats, smpl_stats) = if snp {
            (&mut self.af_snps, &mut self.smpl_snps)
        } else {
            (&mut self.af_indels, &mut self.smpl_indels)
        };
        for (is, (&s0, &s1)) in a.samples.iter().zip(b.samples).enumerate() {
            let gt0 = utils::gt_type(a.rec, fmt0, s0).0;
            let gt1 = utils::gt_type(b.rec, fmt1, s1).0;
            let (i0, i1) = (class(gt0), class(gt1));
            af_stats[iaf1].gt2gt[i0][i1] += 1;
            smpl_stats[is].gt2gt[i0][i1] += 1;
            if gt0 == GenotypeType::Unknown || gt1 == GenotypeType::Unknown {
                continue;
            }
            if is_haploid(gt0) != is_haploid(gt1) {
                continue;
            }
            let y = f64::from(dosage(gt0));
            let x = f64::from(dosage(gt1));
            for c in [&mut smpl_stats[is], &mut af_stats[iaf1]] {
                c.yx += y * x;
                c.x += x;
                c.xx += x * x;
                c.y += y;
                c.yy += y * y;
                c.n += 1.0;
            }
        }
        if let Some((out, names)) = verbose {
            let chrom = a.header.contig_name(a.rec.rid()).unwrap_or(b"");
            let pos = a.rec.pos() + 1;
            let (mut nm, mut nmm, mut nrefm) = (0u32, 0u32, 0u32);
            for (is, (&s0, &s1)) in a.samples.iter().zip(b.samples).enumerate() {
                let gt = utils::gt_type(a.rec, fmt0, s0).0;
                if gt == GenotypeType::Unknown {
                    continue;
                }
                let gt2 = utils::gt_type(b.rec, fmt1, s1).0;
                if gt2 == GenotypeType::Unknown {
                    continue;
                }
                if gt != gt2 {
                    nmm += 1;
                    out.write_all(b"DBG\t")?;
                    out.write_all(chrom)?;
                    write!(out, "\t{pos}\t")?;
                    out.write_all(&names[is])?;
                    writeln!(out, "\t{}\t{}", gt_code(gt), gt_code(gt2))?;
                } else {
                    if gt != GenotypeType::HomRef {
                        nrefm += 1;
                    }
                    nm += 1;
                }
            }
            let nrd = if nrefm + nmm == 0 {
                0.0
            } else {
                (100.0 * f64::from(nmm) / f64::from(nrefm + nmm)) as f32
            };
            out.write_all(b"PSD\t")?;
            out.write_all(chrom)?;
            writeln!(out, "\t{pos}\t{nm}\t{nmm}\t{:.6}", f64::from(nrd))?;
        }
        Ok(())
    }
}
