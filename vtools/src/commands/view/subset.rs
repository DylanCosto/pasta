//! The per-record selection and subsetting of bcftools view (`subset_vcf`).

use anyhow::bail;
use pasta::variant::{
    Header, Record, Unpack, VariantMatch, VariantType,
    encoding::{INT8_VECTOR_END, INT16_VECTOR_END, INT32_VECTOR_END},
    utils::{self, AcSource, GenotypeType},
};

use super::{GtRequirement, Selection};
use crate::filter::{FilterLogic, RecordFilter};

/// Which allele count a `-c/-C/-q/-Q` threshold applies to (bcftools' `ALLELE_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlleleKind {
    /// Non-reference alleles (`nref`).
    NonRef,
    /// The least frequent allele (`minor`).
    Minor,
    /// The first alternate allele (`alt1`).
    Alt1,
    /// The most frequent allele (`major`).
    Major,
    /// All but the most frequent allele (`nonmajor`).
    NonMajor,
}

impl AlleleKind {
    /// Parses an allele type (bcftools' `set_allele_type`).
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        Ok(match s {
            "minor" => Self::Minor,
            "alt1" => Self::Alt1,
            "nref" => Self::NonRef,
            "major" => Self::Major,
            "nonmajor" => Self::NonMajor,
            _ => bail!(
                "Error: allele type not recognised. Expected one of nref|alt1|minor|major|nonmajor, got \"{s}\"."
            ),
        })
    }
}

/// Variant types to include or exclude (bcftools' `view_vtype_t`).
#[derive(Clone, Copy, Debug, Default)]
pub struct VariantTypes {
    mask: VariantType,
    reference: bool,
}

impl VariantTypes {
    /// Parses a `-v/-V` list (bcftools' `parse_vtype`).
    pub fn parse(list: &str) -> anyhow::Result<Self> {
        let mut t = Self::default();
        // An empty list is one empty type, as `hts_readlist` splits it.
        for ty in list.split(',') {
            match ty {
                "snps" | "snp" => t.mask |= VariantType::SNP,
                "mnps" | "mnp" => t.mask |= VariantType::MNP,
                "indels" | "indel" => t.mask |= VariantType::INDEL,
                "other" => t.mask |= VariantType::OTHER,
                "bnd" => t.mask |= VariantType::BND,
                "ref" => t.reference = true,
                "overlap" => t.mask |= VariantType::OVERLAP,
                _ => bail!(
                    "[E::{ty}] unknown type. Accepted types are snps, indels, mnps, ref, bnd, other, overlap"
                ),
            }
        }
        Ok(t)
    }

    fn is_set(&self) -> bool {
        !self.mask.is_empty() || self.reference
    }

    /// bcftools' `has_vtype`.
    fn matches(&self, rec: &mut Record) -> bool {
        if self.reference
            && rec
                .has_variant_types(VariantType::empty(), VariantMatch::Exact)
                .0
        {
            return true;
        }
        !self.mask.is_empty() && rec.has_variant_types(self.mask, VariantMatch::Overlap).0
    }
}

/// bcftools view's record selection state.
#[derive(Default)]
pub struct Subsetter {
    pub min_alleles: i64,
    pub max_alleles: i64,
    pub known: bool,
    pub novel: bool,
    pub include: VariantTypes,
    pub exclude: VariantTypes,
    pub calc_ac: bool,
    pub update_info: bool,
    pub sample_subset: bool,
    pub imap: Vec<i32>,
    pub private_vars: Option<Selection>,
    pub uncalled: Option<Selection>,
    pub phased: Option<Selection>,
    pub gt_type: Option<GtRequirement>,
    pub min_ac: Option<(i32, AlleleKind)>,
    pub max_ac: Option<(i32, AlleleKind)>,
    pub min_af: Option<(f32, AlleleKind)>,
    pub max_af: Option<(f32, AlleleKind)>,
    pub trim_alts: bool,
    pub trim_star_allele: u8,
    pub sites_only: bool,
    ac: Vec<i32>,
    ac_sub: Vec<i32>,
}

impl Subsetter {
    /// Returns whether any selection or subsetting applies (otherwise every record passes
    /// unchanged).
    pub fn is_active(&self) -> bool {
        self.min_alleles != 0
            || self.max_alleles != 0
            || self.known
            || self.novel
            || self.include.is_set()
            || self.exclude.is_set()
            || self.calc_ac
            || self.sample_subset
            || self.private_vars.is_some()
            || self.uncalled.is_some()
            || self.phased.is_some()
            || self.gt_type.is_some()
            || self.trim_alts
            || self.trim_star_allele > 0
            || self.sites_only
    }

    /// Decides whether to output a record, subsetting and updating it (bcftools'
    /// `subset_vcf`).
    pub fn keep(
        &mut self,
        hdr: &Header,
        hsub: Option<&Header>,
        line: &mut Record,
        filter: &mut Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    ) -> anyhow::Result<bool> {
        let n_allele = line.allele_count() as i64;
        if self.min_alleles != 0 && n_allele < self.min_alleles {
            return Ok(false);
        }
        if self.max_alleles != 0 && n_allele > self.max_alleles {
            return Ok(false);
        }
        if self.novel || self.known {
            line.unpack(Unpack::STR);
            let is_novel = line.id() == b".";
            if self.novel && !is_novel {
                return Ok(false);
            }
            if self.known && is_novel {
                return Ok(false);
            }
        }
        if self.include.is_set() && !self.include.matches(line) {
            return Ok(false);
        }
        if self.exclude.is_set() && self.exclude.matches(line) {
            return Ok(false);
        }
        if let Some((f, logic)) = filter {
            let ret = f.pass(hdr, line)?;
            match logic {
                FilterLogic::Include if !ret => return Ok(false),
                FilterLogic::Exclude if ret => return Ok(false),
                _ => {}
            }
        }

        let mut an: i64 = 0;
        let mut non_ref_ac: i64 = 0;
        if self.calc_ac {
            utils::calc_ac(hdr, line, &mut self.ac, AcSource::InfoThenGenotypes)
                .map_err(anyhow::Error::new)?;
            non_ref_ac = self.ac[1..].iter().map(|&v| i64::from(v)).sum();
            an = self.ac.iter().map(|&v| i64::from(v)).sum();
        } else {
            self.ac.clear();
            self.ac.resize(line.allele_count().max(1), 0);
        }

        let mut update_ac = self.calc_ac;
        if self.sample_subset {
            line.subset(&self.imap);
            if self.calc_ac && utils::gt_field(hdr, line).is_none() {
                update_ac = false;
            }
            if update_ac {
                let hsub = hsub.expect("a subset header exists with sample subsetting");
                utils::calc_ac(hsub, line, &mut self.ac_sub, AcSource::Genotypes)
                    .map_err(anyhow::Error::new)?;
                an = self.ac_sub.iter().map(|&v| i64::from(v)).sum();
                self.ac.clone_from(&self.ac_sub);
                let non_ref_ac_sub: i64 = self.ac_sub[1..].iter().map(|&v| i64::from(v)).sum();
                match self.private_vars {
                    Some(Selection::Include)
                        if !(non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub) =>
                    {
                        return Ok(false);
                    }
                    Some(Selection::Exclude)
                        if non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub =>
                    {
                        return Ok(false);
                    }
                    _ => {}
                }
                non_ref_ac = non_ref_ac_sub;
            }
        }

        if let Some(req) = self.gt_type
            && let Some(gt) = utils::gt_field(hdr, line)
        {
            let (mut nhet, mut nhom, mut nmiss) = (false, false, false);
            // bcftools loops over the original header's samples, reading past the data after
            // -s subsetting; this loops over the record's samples.
            for i in 0..line.sample_count() {
                let (ty, _, _) = utils::gt_type(line, &gt, i);
                match ty {
                    GenotypeType::HetRefAlt | GenotypeType::HetAltAlt => {
                        if req == GtRequirement::NoHet {
                            return Ok(false);
                        }
                        nhet = true;
                    }
                    GenotypeType::Unknown => {
                        if req == GtRequirement::NoMissing {
                            return Ok(false);
                        }
                        nmiss = true;
                    }
                    _ => {
                        if req == GtRequirement::NoHom {
                            return Ok(false);
                        }
                        nhom = true;
                    }
                }
            }
            match req {
                GtRequirement::NeedHom if !nhom => return Ok(false),
                GtRequirement::NeedHet if !nhet => return Ok(false),
                GtRequirement::NeedMissing if !nmiss => return Ok(false),
                _ => {}
            }
        }

        let (mut minor_ac, mut major_ac) = (0i64, 0i64);
        if self.calc_ac {
            minor_ac = i64::from(self.ac[0]);
            major_ac = i64::from(self.ac[0]);
            for &v in &self.ac[1..line.allele_count().max(1)] {
                minor_ac = minor_ac.min(i64::from(v));
                major_ac = major_ac.max(i64::from(v));
            }
        }
        let n_allele = line.allele_count();
        let alt1 = if n_allele >= 2 {
            i64::from(self.ac[1])
        } else {
            0
        };
        if let Some((min, kind)) = self.min_ac.filter(|(v, _)| *v != -1) {
            let min = i64::from(min);
            let fail = match kind {
                AlleleKind::NonRef => min > non_ref_ac,
                AlleleKind::Minor => min > minor_ac,
                AlleleKind::Alt1 => n_allele < 2 || min > alt1,
                AlleleKind::Major => min > major_ac,
                AlleleKind::NonMajor => min > an - major_ac,
            };
            if fail {
                return Ok(false);
            }
        }
        if let Some((max, kind)) = self.max_ac.filter(|(v, _)| *v != -1) {
            let max = i64::from(max);
            let fail = match kind {
                AlleleKind::NonRef => max < non_ref_ac,
                AlleleKind::Minor => max < minor_ac,
                AlleleKind::Alt1 => n_allele >= 2 && max < alt1,
                AlleleKind::Major => max < major_ac,
                AlleleKind::NonMajor => max < an - major_ac,
            };
            if fail {
                return Ok(false);
            }
        }
        let freq = |x: i64| x as f64 / an as f64;
        if let Some((min, kind)) = self.min_af.filter(|(v, _)| *v != -1.0) {
            if an == 0 {
                return Ok(false);
            }
            let min = f64::from(min);
            let fail = match kind {
                AlleleKind::NonRef => min > freq(non_ref_ac),
                AlleleKind::Minor => min > freq(minor_ac),
                AlleleKind::Alt1 => n_allele < 2 || min > freq(alt1),
                AlleleKind::Major => min > freq(major_ac),
                AlleleKind::NonMajor => min > freq(an - major_ac),
            };
            if fail {
                return Ok(false);
            }
        }
        if let Some((max, kind)) = self.max_af.filter(|(v, _)| *v != -1.0) {
            if an == 0 {
                return Ok(false);
            }
            let max = f64::from(max);
            let fail = match kind {
                AlleleKind::NonRef => max < freq(non_ref_ac),
                AlleleKind::Minor => max < freq(minor_ac),
                AlleleKind::Alt1 => n_allele >= 2 && max < freq(alt1),
                AlleleKind::Major => max < freq(major_ac),
                AlleleKind::NonMajor => max < freq(an - major_ac),
            };
            if fail {
                return Ok(false);
            }
        }
        match self.uncalled {
            Some(Selection::Include) if an > 0 => return Ok(false),
            Some(Selection::Exclude) if an == 0 => return Ok(false),
            _ => {}
        }
        if update_ac && self.update_info {
            let acs = if n_allele > 1 {
                &self.ac[1..n_allele]
            } else {
                &[][..]
            };
            line.update_info_i32(hdr, b"AC", acs)
                .map_err(anyhow::Error::new)?;
            line.update_info_i32(hdr, b"AN", &[an as i32])
                .map_err(anyhow::Error::new)?;
        }
        if self.trim_alts {
            let h = hsub.unwrap_or(hdr);
            if utils::trim_alleles(h, line).is_err() {
                bail!(
                    "Error: Could not trim alleles at {}:{}",
                    String::from_utf8_lossy(line.chrom_safe(h)),
                    line.pos() + 1
                );
            }
        }
        if self.trim_star_allele > 0 && (line.allele_count() > 2 || self.trim_star_allele > 1) {
            line.unpack(Unpack::STR);
            let unseen = (1..line.allele_count())
                .find(|&i| matches!(line.allele(i), b"<*>" | b"<NON_REF>" | b"<X>"));
            if let Some(i) = unseen {
                let mut rm = vec![false; line.allele_count()];
                rm[i] = true;
                // The record has the subset's samples now. bcftools passes the input header
                // here and fails on any -s subset; the subset header trims the allele.
                let h = hsub.unwrap_or(hdr);
                if utils::remove_allele_set(h, line, &rm).is_err() {
                    bail!(
                        "[subset_vcf] Error: failed to trim the unobserved allele at {}:{}",
                        String::from_utf8_lossy(line.chrom_safe(h)),
                        line.pos() + 1
                    );
                }
            }
        }
        if let Some(sel) = self.phased {
            let phased = all_phased(hdr, line);
            match sel {
                Selection::Include if !phased => return Ok(false),
                Selection::Exclude if phased => return Ok(false),
                _ => {}
            }
        }
        if self.sites_only {
            line.subset(&[]);
        }
        Ok(true)
    }
}

/// Returns whether every sample's genotype is phased; haploid genotypes count as phased
/// (bcftools' `bcf_all_phased`).
fn all_phased(header: &Header, line: &mut Record) -> bool {
    line.unpack(Unpack::FMT);
    let Some(gt) = utils::gt_field(header, line) else {
        return true;
    };
    let n = gt.n();
    for s in 0..line.sample_count() {
        let p = line.format_sample_data(&gt, s);
        let value = |i: usize| -> i32 {
            match gt.type_code() {
                1 => {
                    let v = p[i] as i8;
                    if v == INT8_VECTOR_END {
                        INT32_VECTOR_END
                    } else {
                        i32::from(v)
                    }
                }
                2 => {
                    let v = i16::from_le_bytes([p[2 * i], p[2 * i + 1]]);
                    if v == INT16_VECTOR_END {
                        INT32_VECTOR_END
                    } else {
                        i32::from(v)
                    }
                }
                _ => i32::from_le_bytes([p[4 * i], p[4 * i + 1], p[4 * i + 2], p[4 * i + 3]]),
            }
        };
        let mut sample_phased = false;
        for i in 0..n {
            let v = value(i);
            if n == 1 || (v == INT32_VECTOR_END && i == 1) {
                sample_phased = true;
                break;
            }
            if v == INT32_VECTOR_END {
                break;
            }
            if v >> 1 == 0 {
                continue;
            }
            if v & 1 != 0 {
                sample_phased = true;
                break;
            }
        }
        if !sample_phased {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use pasta::variant::vcf;

    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">\n\
##FORMAT=<ID=XA,Number=A,Type=Integer,Description=\"x\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

    /// Runs `line` through `sub` with the samples subset to `sample`, as `view -s -I` does,
    /// and returns the record as VCF text.
    fn keep_subset(mut sub: Subsetter, sample: &str, line: &str) -> String {
        let mut hdr = Header::parse(HEADER.as_bytes()).expect("the header parses");
        let (hsub, imap) = hdr
            .subset(&[sample.as_bytes()])
            .expect("the sample is in the header");
        sub.sample_subset = true;
        sub.imap = imap;
        let mut rec = Record::new();
        vcf::parse(&mut hdr, &mut line.as_bytes().to_vec(), &mut rec).expect("the line parses");
        let kept = sub
            .keep(&hdr, Some(&hsub), &mut rec, &mut None)
            .expect("the record is processed");
        assert!(kept, "the record is kept");
        let mut out = Vec::new();
        vcf::format_record(&hsub, &rec, &mut out).expect("the record formats");
        String::from_utf8(out).expect("the record is text")
    }

    #[test]
    fn trim_unseen_allele_after_sample_subset() {
        let sub = || Subsetter {
            trim_star_allele: 1,
            ..Subsetter::default()
        };
        for alt in ["<NON_REF>", "<*>"] {
            let line = format!(
                "1\t10\t.\tA\tC,{alt}\t50\tPASS\tDP=5\tGT:AD:PL:XA\t\
0/1:3,2,0:10,0,20,30,40,50:1,2\t1/1:0,5,0:60,30,0,70,80,90:3,4"
            );
            assert_eq!(
                keep_subset(sub(), "s1", &line),
                "1\t10\t.\tA\tC\t50\tPASS\tDP=5\tGT:AD:PL:XA\t0/1:3,2:10,0,20:1\n",
                "ALT {alt}"
            );
            assert_eq!(
                keep_subset(sub(), "s2", &line),
                "1\t10\t.\tA\tC\t50\tPASS\tDP=5\tGT:AD:PL:XA\t1/1:0,5:60,30,0:3\n",
                "ALT {alt}"
            );
        }
        let line = "1\t10\t.\tA\tC,<NON_REF>\t50\tPASS\tDP=5\tGT\t0/1\t1/1";
        assert_eq!(
            keep_subset(sub(), "s2", line),
            "1\t10\t.\tA\tC\t50\tPASS\tDP=5\tGT\t1/1\n"
        );
    }
}
