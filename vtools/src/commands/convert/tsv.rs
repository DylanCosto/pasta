//! Text to VCF: `--tsv2vcf`, `--gensample2vcf`, `--hapsample2vcf`, and
//! `--haplegendsample2vcf`. Lines are split into whitespace-separated columns, each filled into
//! the record by a setter (bcftools' `tsv2vcf.c` and the setters of `vcfconvert.c`).

use std::path::Path;

use anyhow::bail;
use pasta::variant::{
    Genotype, Header, Record,
    encoding::INT32_VECTOR_END,
    index::split_idx_name,
    io::Writer,
    number::{strtod, strtol},
    readlist::{self, ListSource},
};

use super::{
    CloseMessage, Settings, finish_output, open_output,
    text::{LineReader, Reference},
    write_record,
};
use crate::version;

/// Writes a record built from text. Its fields exist only decoded, so they are encoded first
/// (HTSlib formats VCF from the decoded fields).
fn write_built(
    writer: &mut Writer,
    header: &Header,
    rec: &mut Record,
    s: &Settings,
    func: &str,
) -> anyhow::Result<()> {
    rec.sync();
    write_record(writer, header, rec, s, func)
}

/// The byte at `i`, or NUL past the end (the C string's terminator).
fn at(s: &[u8], i: usize) -> u8 {
    s.get(i).copied().unwrap_or(0)
}

/// C's `isspace`.
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// The text from `i` to the end of the line (what bcftools prints with `%s`).
fn rest(s: &[u8], i: usize) -> String {
    String::from_utf8_lossy(s.get(i..).unwrap_or(b"")).into_owned()
}

/// A column setter (bcftools' `tsv_setter_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Setter {
    Chrom,
    Pos,
    Id,
    Ref,
    Alt,
    Aa,
    ChromPosRefAlt,
    ChromPosRefAltOrChrom,
    ChromPosRefAltOrId,
    ChromPosRefAltIdOrDie,
    VerifyPos,
    VerifyRefAlt,
    GtGp,
    Haps,
}

/// Columns and the current field (bcftools' `tsv_t`): `ss..se` is the field being set.
struct Tsv {
    cols: Vec<(Option<Vec<u8>>, Option<Setter>)>,
    ss: usize,
    se: usize,
}

impl Tsv {
    /// Parses a comma-separated column list; `-` columns are skipped (`tsv_init`).
    fn new(spec: &str) -> Self {
        let mut names: Vec<&str> = spec.split(',').collect();
        if names.last() == Some(&"") {
            names.pop();
        }
        let cols = names
            .into_iter()
            .map(|n| {
                (
                    (!n.eq_ignore_ascii_case("-")).then(|| n.as_bytes().to_vec()),
                    None,
                )
            })
            .collect();
        Self { cols, ss: 0, se: 0 }
    }

    /// Sets the setter of the first column named `id`; `false` if there is none
    /// (`tsv_register`).
    fn register(&mut self, id: &str, setter: Setter) -> bool {
        for col in &mut self.cols {
            if col
                .0
                .as_ref()
                .is_some_and(|n| n.eq_ignore_ascii_case(id.as_bytes()))
            {
                col.1 = Some(setter);
                return true;
            }
        }
        false
    }

    /// Moves to the next field (`tsv_next`); `false` when there is none.
    fn next(&mut self, s: &[u8]) -> bool {
        if at(s, self.se) == 0 {
            return false;
        }
        if self.ss == self.se {
            while at(s, self.se) != 0 && !is_space(at(s, self.se)) {
                self.se += 1;
            }
            return true;
        }
        while at(s, self.se) != 0 && is_space(at(s, self.se)) {
            self.se += 1;
        }
        self.ss = self.se;
        while at(s, self.se) != 0 && !is_space(at(s, self.se)) {
            self.se += 1;
        }
        true
    }

    /// Fills the record from a line (`tsv_parse`): `true` on success, `false` when a setter
    /// rejects the line or no column has a setter.
    fn parse(&mut self, ctx: &mut Ctx<'_>, rec: &mut Record, s: &[u8]) -> anyhow::Result<bool> {
        let mut status = 0;
        self.ss = 0;
        self.se = 0;
        let mut icol = 0;
        while at(s, self.ss) != 0 && icol < self.cols.len() {
            while at(s, self.se) != 0 && !is_space(at(s, self.se)) {
                self.se += 1;
            }
            if let Some(setter) = self.cols[icol].1 {
                if !ctx.set(setter, self, rec, s)? {
                    return Ok(false);
                }
                status += 1;
            }
            while at(s, self.se) != 0 && is_space(at(s, self.se)) {
                self.se += 1;
            }
            self.ss = self.se;
            icol += 1;
        }
        Ok(status > 0)
    }
}

/// Genotype counts of `--tsv2vcf`.
#[derive(Default)]
struct Counts {
    missing: u64,
    hom_rr: u64,
    het_ra: u64,
    hom_aa: u64,
    het_aa: u64,
}

/// What the setters share (the relevant part of bcftools' `args_t`).
struct Ctx<'a> {
    header: &'a Header,
    vcf_ids: bool,
    rev_als: bool,
    gts: Vec<i32>,
    gps: Vec<f32>,
    tsv_ref: Vec<u8>,
    tsv_alt: Vec<u8>,
    alleles: Vec<u8>,
    reference: Option<&'a mut Reference>,
    base: Vec<u8>,
    counts: Counts,
}

impl<'a> Ctx<'a> {
    fn new(header: &'a Header, settings: &Settings) -> Self {
        Self {
            header,
            vcf_ids: settings.vcf_ids,
            rev_als: false,
            gts: Vec::new(),
            gps: Vec::new(),
            tsv_ref: Vec::new(),
            tsv_alt: Vec::new(),
            alleles: Vec::new(),
            reference: None,
            base: Vec::new(),
            counts: Counts::default(),
        }
    }

    /// Runs a setter; `false` when it rejects the line.
    fn set(
        &mut self,
        setter: Setter,
        tsv: &mut Tsv,
        rec: &mut Record,
        s: &[u8],
    ) -> anyhow::Result<bool> {
        let field = &s[tsv.ss..tsv.se.min(s.len())];
        Ok(match setter {
            Setter::Chrom => self.set_chrom(rec, field),
            Setter::Pos => {
                let (v, end) = strtol(s, tsv.ss);
                rec.set_pos(v - 1);
                end != tsv.ss
            }
            Setter::Id => {
                rec.update_id(Some(field));
                true
            }
            Setter::Ref => {
                self.tsv_ref.extend_from_slice(field);
                if !self.tsv_alt.is_empty() {
                    self.set_ref_alt(rec);
                }
                true
            }
            Setter::Alt => {
                self.tsv_alt.extend_from_slice(field);
                if !self.tsv_ref.is_empty() {
                    self.set_ref_alt(rec);
                }
                true
            }
            Setter::Aa => self.set_aa(tsv, rec, s)?,
            Setter::ChromPosRefAlt => {
                if !self.set_chrom_pos_ref_alt(tsv, rec, s)? {
                    bail!(
                        "Could not parse the CHROM:POS_REF_ALT[_END] string: {}",
                        rest(s, tsv.ss)
                    );
                }
                true
            }
            Setter::ChromPosRefAltOrChrom => {
                self.set_chrom_pos_ref_alt(tsv, rec, s)? || self.set_chrom(rec, field)
            }
            Setter::ChromPosRefAltOrId => {
                if self.set_chrom_pos_ref_alt(tsv, rec, s)? {
                    return Ok(true);
                }
                // Mark the record as unset.
                rec.set_pos(-1);
                if self.vcf_ids {
                    rec.update_id(Some(field));
                }
                true
            }
            Setter::ChromPosRefAltIdOrDie => {
                if rec.pos() != -1 {
                    if self.vcf_ids {
                        rec.update_id(Some(field));
                    }
                    return Ok(true);
                }
                if !self.set_chrom_pos_ref_alt(tsv, rec, s)? {
                    bail!(
                        "Could not parse the CHROM:POS_REF_ALT[_END] string: {}",
                        rest(s, tsv.ss)
                    );
                }
                true
            }
            Setter::VerifyPos => {
                let (v, end) = strtol(s, tsv.ss);
                if end == tsv.ss {
                    bail!("Could not parse POS: {}", rest(s, tsv.ss));
                }
                if rec.pos() != i64::from(v as i32) - 1 {
                    bail!("POS mismatch: {}", rest(s, tsv.ss));
                }
                true
            }
            Setter::VerifyRefAlt => {
                self.verify_ref_alt(tsv, rec, s)?;
                true
            }
            Setter::GtGp => self.set_gt_gp(tsv, rec, s)?,
            Setter::Haps => self.set_haps(tsv, rec, s)?,
        })
    }

    fn set_chrom(&self, rec: &mut Record, field: &[u8]) -> bool {
        match self.header.contig_id(field) {
            Some(rid) => {
                rec.set_rid(rid as i32);
                true
            }
            None => {
                rec.set_rid(-1);
                false
            }
        }
    }

    /// Sets REF and ALT once both columns are read (`_set_ref_alt`).
    fn set_ref_alt(&mut self, rec: &mut Record) {
        self.alleles.clear();
        self.alleles.extend_from_slice(&self.tsv_ref);
        if self.tsv_alt != b"." && self.tsv_alt != self.tsv_ref {
            self.alleles.push(b',');
            self.alleles.extend_from_slice(&self.tsv_alt);
        }
        rec.update_alleles_str(self.header, &self.alleles);
        self.tsv_ref.clear();
        self.tsv_alt.clear();
        self.alleles.clear();
    }

    /// Parses `CHROM:POS_REF_ALT[_END]` (`_set_chrom_pos_ref_alt`); `false` when it does not
    /// parse.
    fn set_chrom_pos_ref_alt(
        &mut self,
        tsv: &Tsv,
        rec: &mut Record,
        s: &[u8],
    ) -> anyhow::Result<bool> {
        let fe = tsv.se;
        let mut se = tsv.ss;
        while se < fe && at(s, se) != b':' {
            se += 1;
        }
        if at(s, se) != b':' {
            return Ok(false);
        }
        let Some(rid) = self.header.contig_id(&s[tsv.ss..se]) else {
            return Ok(false);
        };
        let (pos, ss) = strtol(s, se + 1);
        if ss == se + 1 {
            return Ok(false);
        }
        let pos = pos - 1;

        self.alleles.clear();
        let mut ss = ss + 1;
        let mut se = ss;
        while se < fe && at(s, se) != b'_' {
            se += 1;
        }
        if at(s, se) != b'_' {
            return Ok(false);
        }
        self.alleles.extend_from_slice(&s[ss..se]);
        se += 1;
        ss = se;
        while se < fe && at(s, se) != b'_' && !is_space(at(s, se)) {
            se += 1;
        }
        self.alleles.push(b',');
        self.alleles.extend_from_slice(&s[ss..se.min(s.len())]);

        // END, optional.
        if at(s, se) == b'_' {
            let (end, ss) = strtol(s, se + 1);
            if ss == se + 1 || !(1..=i64::from(i32::MAX)).contains(&end) {
                return Ok(false);
            }
            rec.update_info_i32(self.header, b"END", &[end as i32])
                .map_err(|e| anyhow::anyhow!("{e}"))?;
        }
        rec.set_rid(rid as i32);
        rec.set_pos(pos);
        let alleles = std::mem::take(&mut self.alleles);
        rec.update_alleles_str(self.header, &alleles);
        self.alleles = alleles;
        Ok(true)
    }

    /// Checks REF and ALT against the record, noting swapped alleles
    /// (`tsv_setter_verify_ref_alt`); consumes the ALT column too.
    fn verify_ref_alt(&mut self, tsv: &mut Tsv, rec: &Record, s: &[u8]) -> anyhow::Result<()> {
        let allele = |i: usize| -> &[u8] {
            if i < rec.allele_count() {
                rec.allele(i)
            } else {
                b""
            }
        };
        self.rev_als = false;
        let field = &s[tsv.ss..tsv.se];
        if field != allele(0) {
            if field != allele(1) {
                bail!(
                    "REF/ALT mismatch: [{}][{}]",
                    rest(s, tsv.ss),
                    String::from_utf8_lossy(allele(1))
                );
            }
            self.rev_als = true;
        }
        while at(s, tsv.se) != 0 && is_space(at(s, tsv.se)) {
            tsv.se += 1;
        }
        tsv.ss = tsv.se;
        while at(s, tsv.se) != 0 && !is_space(at(s, tsv.se)) {
            tsv.se += 1;
        }
        let field = &s[tsv.ss..tsv.se];
        let expected = if self.rev_als { allele(0) } else { allele(1) };
        if field != expected {
            bail!(
                "REF/ALT mismatch: [{}][{}]",
                rest(s, tsv.ss),
                String::from_utf8_lossy(expected)
            );
        }
        Ok(())
    }

    /// Reads three genotype probabilities per sample into GT and GP (`tsv_setter_gt_gp`).
    fn set_gt_gp(&mut self, tsv: &mut Tsv, rec: &mut Record, s: &[u8]) -> anyhow::Result<bool> {
        let nsamples = self.header.sample_count();
        self.gts.resize(nsamples * 2, 0);
        self.gps.resize(nsamples * 3, 0.0);
        for i in 0..nsamples {
            let mut vals = [0f32; 3];
            for (k, what) in ["first", "second", "third"].into_iter().enumerate() {
                let (v, end) = strtod(s, tsv.ss);
                tsv.se = end;
                if end == tsv.ss {
                    eprintln!("Could not parse {what} value of {}-th sample", i + 1);
                    return Ok(false);
                }
                tsv.ss = end + 1;
                vals[k] = v as f32;
            }
            let [mut aa, ab, mut bb] = vals;
            if self.rev_als {
                std::mem::swap(&mut aa, &mut bb);
            }
            self.gps[3 * i..3 * i + 3].copy_from_slice(&[aa, ab, bb]);
            let (a, b) = if aa >= ab {
                if aa >= bb { (0, 0) } else { (1, 1) }
            } else if ab >= bb {
                (0, 1)
            } else {
                (1, 1)
            };
            self.gts[2 * i] = Genotype::unphased(a).0;
            self.gts[2 * i + 1] = Genotype::unphased(b).0;
        }
        if at(s, tsv.se) != 0 {
            bail!("Could not parse: {}", rest(s, tsv.ss));
        }
        if rec.update_genotypes(self.header, &self.gts).is_err() {
            bail!("Could not update GT field");
        }
        if rec
            .update_format_f32(self.header, b"GP", &self.gps)
            .is_err()
        {
            bail!("Could not update GP field");
        }
        Ok(true)
    }

    /// Reads two haplotype columns per sample (`0`, `1`, `?`, `-`, with `*` for unphased)
    /// into GT (`tsv_setter_haps`).
    fn set_haps(&mut self, tsv: &Tsv, rec: &mut Record, s: &[u8]) -> anyhow::Result<bool> {
        let nsamples = self.header.sample_count();
        self.gts.resize(nsamples * 2, 0);
        let (a0, a1) = if self.rev_als {
            (Genotype::phased(1).0, Genotype::phased(0).0)
        } else {
            (Genotype::phased(0).0, Genotype::phased(1).0)
        };
        let mut nup = 0;
        for i in 0..nsamples {
            let ss = tsv.ss + 4 * i + nup;
            let c = |k: usize| at(s, ss + k);
            let mut up = 0;
            for all in 0..2 {
                // Checking for premature ending.
                if c(0) == 0 || c(1) == 0 || c(2) == 0 || (up > 0 && (c(3) == 0 || c(4) == 0)) {
                    eprint!(
                        "Wrong number of fields at {}-th sample ([{}][{}][{}]). ",
                        i + 1,
                        char::from(c(0)),
                        char::from(c(1)),
                        char::from(c(2))
                    );
                    return Ok(false);
                }
                let ch = c(all * 2 + up);
                self.gts[2 * i + all] = match ch {
                    b'0' => a0,
                    b'1' => a1,
                    // A phased missing allele.
                    b'?' => Genotype::phased(-1).0,
                    b'-' => INT32_VECTOR_END,
                    _ => {
                        eprintln!("Could not parse: [{}][{}]", char::from(ch), rest(s, tsv.ss));
                        return Ok(false);
                    }
                };
                if c(all * 2 + up + 1) == b'*' {
                    up += 1;
                }
            }
            if up > 0 && up != 2 {
                eprint!(
                    "Missing unphased marker '*': [{}][{}]",
                    char::from(c(2 + up)),
                    rest(s, tsv.ss)
                );
                return Ok(false);
            }
            // Unphased alleles.
            if up > 0 {
                for g in &mut self.gts[2 * i..2 * i + 2] {
                    *g = Genotype::unphased(Genotype(*g).allele()).0;
                }
            }
            nup += up;
        }
        let last = (nsamples.max(1) - 1) * 4 + 3 + nup;
        if at(s, tsv.ss + last) != 0 {
            eprint!("nup: {nup}");
            eprint!(
                "Wrong number of fields ({}-th column = [{}]). ",
                nsamples * 2,
                char::from(at(s, tsv.ss + (nsamples.max(1) - 1) * 4 + nup))
            );
            return Ok(false);
        }
        if rec.update_genotypes(self.header, &self.gts).is_err() {
            bail!("Could not update GT field");
        }
        Ok(true)
    }

    /// Reads one genotype column per sample as letters (`AG`, `T`, `--`) against the reference
    /// base (`tsv_setter_aa`); `false` for indels, which are skipped.
    fn set_aa(&mut self, tsv: &mut Tsv, rec: &mut Record, s: &[u8]) -> anyhow::Result<bool> {
        let header = self.header;
        let chrom = rec.chrom_safe(header).to_vec();
        let location = || format!("{}:{}", String::from_utf8_lossy(&chrom), rec.pos() + 1);
        let reference = self
            .reference
            .as_deref_mut()
            .expect("--tsv2vcf loads the reference");
        if !reference.fetch(&chrom, rec.pos(), &mut self.base)? {
            bail!("faidx_fetch_seq failed at {}", location());
        }
        let ref_base = self.base.first().copied().unwrap_or(0).to_ascii_uppercase();
        let iref = acgt_to_5(ref_base);
        // The allele number of each of A, C, G, T, N, in order of appearance.
        let mut alleles = [-1i32; 5];
        let mut nals = 1;
        alleles[iref] = 0;
        let nsamples = header.sample_count();
        self.gts.resize(nsamples * 2, 0);
        for i in 0..nsamples {
            if i > 0 && !tsv.next(s) {
                bail!("Too few columns for {nsamples} samples at {}", location());
            }
            if tsv.se - tsv.ss > 2 {
                bail!(
                    "Error parsing the site {}, expected two characters",
                    location()
                );
            }
            let (c0, c1) = (at(s, tsv.ss), at(s, tsv.ss + 1));
            if c0 == b'-' || c0 == b'.' {
                self.gts[2 * i] = Genotype::MISSING.0;
                self.gts[2 * i + 1] = Genotype::MISSING.0;
                self.counts.missing += 1;
                continue;
            }
            if c0 == b'I' || c0 == b'D' {
                // Something else than a SNP.
                return Ok(false);
            }
            let a0 = acgt_to_5(c0.to_ascii_uppercase());
            // Like bcftools, the byte after a one-letter genotype is read as its second
            // allele unless the line ends there.
            let a1 = if c1 != 0 {
                acgt_to_5(c1.to_ascii_uppercase())
            } else {
                a0
            };
            for a in [a0, a1] {
                if alleles[a] < 0 {
                    alleles[a] = nals;
                    nals += 1;
                }
            }
            self.gts[2 * i] = Genotype::unphased(alleles[a0]).0;
            self.gts[2 * i + 1] = if c1 != 0 {
                Genotype::unphased(alleles[a1]).0
            } else {
                INT32_VECTOR_END
            };
            if iref == a0 && iref == a1 {
                self.counts.hom_rr += 1;
            } else if iref == a0 || iref == a1 {
                self.counts.het_ra += 1;
            } else if a0 == a1 {
                self.counts.hom_aa += 1;
            } else {
                self.counts.het_aa += 1;
            }
        }
        // REF, then the ALT bases in the order their allele numbers were given (bcftools
        // lists them in A, C, G, T, N order, which mislabels the genotypes of sites with
        // several ALTs).
        self.alleles.clear();
        self.alleles.push(ref_base);
        for n in 1..nals {
            if let Some(k) = alleles.iter().position(|&a| a == n) {
                self.alleles.push(b',');
                self.alleles.push(b"ACGTN"[k]);
            }
        }
        let alleles = std::mem::take(&mut self.alleles);
        rec.update_alleles_str(header, &alleles);
        self.alleles = alleles;
        if rec
            .update_genotypes(header, &self.gts[..nsamples * 2])
            .is_err()
        {
            bail!("Could not update the GT field");
        }
        Ok(true)
    }
}

/// A, C, G, T as 0-3, anything else 4 (bcftools' `acgt_to_5`).
fn acgt_to_5(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// Splits `FILE1,FILE2` or derives the names from a prefix and suffixes.
fn input_names(infname: &str, suffixes: &[&str]) -> Vec<String> {
    match infname.split_once(',') {
        Some((a, b)) => vec![a.to_owned(), b.to_owned()],
        None => suffixes.iter().map(|s| format!("{infname}{s}")).collect(),
    }
}

/// Opens a text input; `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
fn open_lines(name: &str) -> anyhow::Result<LineReader> {
    LineReader::open(split_idx_name(name).0).map_err(|_| anyhow::anyhow!("Could not read: {name}"))
}

/// Reads a sample file, skipping `skip` header lines; names end at the first space.
fn add_samples(header: &mut Header, sample_fname: &str, skip: usize) -> anyhow::Result<usize> {
    let Ok(lines) = readlist::read_list(sample_fname.as_bytes(), ListSource::File) else {
        bail!("Could not read {sample_fname}");
    };
    for line in lines.iter().skip(skip) {
        let end = line.iter().position(|&c| is_space(c)).unwrap_or(line.len());
        // bcftools ignores failures (an empty or duplicate name is reported and skipped).
        let _ = header.add_sample(&line[..end]);
    }
    Ok(lines.len())
}

/// Whether a converted file has genotype probabilities (`FORMAT/GP`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Gp {
    Declared,
    Omitted,
}

/// The header lines every `*2vcf` conversion starts with.
fn new_header(s: &Settings, chrom: &[u8], gp: Gp) -> anyhow::Result<Header> {
    let mut header = Header::new_for_writing();
    header.append(
        b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
    )?;
    header.append(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    if gp == Gp::Declared {
        header.append(
            b"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">",
        )?;
    }
    let mut contig = b"##contig=<ID=".to_vec();
    contig.extend_from_slice(chrom);
    contig.extend_from_slice(b",length=2147483647>");
    header.append(&contig)?;
    if s.record_cmd_line {
        version::append_version(&mut header, "convert")?;
    }
    Ok(header)
}

/// The first whitespace-delimited token of `s` from `i`: its end.
fn token_end(s: &[u8], mut i: usize) -> usize {
    while i < s.len() && !is_space(s[i]) {
        i += 1;
    }
    i
}

/// IMPUTE2 gen/sample to VCF with GT and GP (bcftools' `gensample_to_vcf`).
pub(super) fn gensample_to_vcf(s: &Settings) -> anyhow::Result<()> {
    let names = input_names(&s.infname, &[".gen.gz", ".samples"]);
    let (gen_fname, sample_fname) = (&names[0], &names[1]);
    let mut gen_fh = open_lines(gen_fname)?;
    let mut line = Vec::new();
    if !gen_fh.next_line(&mut line)? {
        bail!("Empty file: {gen_fname}");
    }
    // The chromosome: the first column with --3N6, else from CHROM:POS_REF_ALT in the first or
    // second column.
    let fail = || {
        anyhow::anyhow!(
            "Could not determine CHROM in {gen_fname}: {}",
            rest(&line, 0)
        )
    };
    let mut se = token_end(&line, 0);
    if se == line.len() {
        return Err(fail());
    }
    let chrom = if s.gen_3n6 {
        line[..se].to_vec()
    } else {
        let mut sb = 0;
        let mut sc = line.iter().position(|&c| c == b':');
        if sc.is_none_or(|c| c > se) {
            se = token_end(&line, se);
            if se == line.len() {
                return Err(fail());
            }
            se += 1;
            sb = se;
            sc = line[sb..].iter().position(|&c| c == b':').map(|i| i + sb);
            if sc.is_none() {
                return Err(fail());
            }
        }
        line[sb..sc.expect("checked above")].to_vec()
    };

    let mut header = new_header(s, &chrom, Gp::Declared)?;
    let mut tsv = if s.gen_3n6 {
        let mut t = Tsv::new("CHROM,CHROM_POS_REF_ALT,ID,POS,REF_ALT,GT_GP");
        t.register("CHROM", Setter::Chrom);
        t
    } else {
        Tsv::new("CHROM_POS_REF_ALT,ID,POS,REF_ALT,GT_GP")
    };
    tsv.register("CHROM_POS_REF_ALT", Setter::ChromPosRefAltOrId);
    tsv.register("ID", Setter::ChromPosRefAltIdOrDie);
    tsv.register("POS", Setter::VerifyPos);
    tsv.register("REF_ALT", Setter::VerifyRefAlt);
    tsv.register("GT_GP", Setter::GtGp);

    add_samples(&mut header, sample_fname, 2)?;
    let mut writer = open_output(s, &header, "gensample_to_vcf")?;
    let mut ctx = Ctx::new(&header, s);
    let mut rec = Record::new();
    let mut total = 0u64;
    let result = (|| -> anyhow::Result<()> {
        loop {
            rec.clear();
            total += 1;
            if !tsv.parse(&mut ctx, &mut rec, &line)? {
                bail!("Error occurred while parsing: {}", rest(&line, 0));
            }
            write_built(&mut writer, &header, &mut rec, s, "gensample_to_vcf")?;
            if !gen_fh.next_line(&mut line)? {
                return Ok(());
            }
        }
    })();
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    finish_output(writer, s, "gensample_to_vcf", CloseMessage::Short)?;
    eprintln!("Number of processed rows: \t{total}");
    Ok(())
}

/// IMPUTE2 hap/legend/sample to VCF (bcftools' `haplegendsample_to_vcf`).
pub(super) fn haplegendsample_to_vcf(s: &Settings) -> anyhow::Result<()> {
    if s.vcf_ids {
        bail!(
            "The option --haplegendsample2vcf cannot be combined with --vcf-ids. This is because the\n\
             ID column must be formatted as \"CHROM:POS_REF_ALT\" to check sanity of the operation"
        );
    }
    let (hap_fname, leg_fname, sample_fname) = match s.infname.split_once(',') {
        Some((a, rest)) => match rest.split_once(',') {
            Some((b, c)) => (a.to_owned(), b.to_owned(), c.to_owned()),
            None => bail!(
                "Could not parse hap/legend/sample file names: {}",
                s.infname
            ),
        },
        None => (
            format!("{}.hap.gz", s.infname),
            format!("{}.legend.gz", s.infname),
            format!("{}.samples", s.infname),
        ),
    };
    let mut hap_fh = open_lines(&hap_fname)?;
    let mut leg_fh = open_lines(&leg_fname)?;
    let mut line = Vec::new();
    // Skip the legend's header line.
    for _ in 0..2 {
        if !leg_fh.next_line(&mut line)? {
            bail!("Empty file: {leg_fname}");
        }
    }
    let Some(colon) = line.iter().position(|&c| c == b':') else {
        bail!("Expected CHROM:POS_REF_ALT in first column of {leg_fname}");
    };
    let chrom = line[..colon].to_vec();

    let mut leg_tsv = Tsv::new("CHROM_POS_REF_ALT,POS,REF_ALT");
    leg_tsv.register("CHROM_POS_REF_ALT", Setter::ChromPosRefAlt);
    leg_tsv.register("POS", Setter::VerifyPos);
    leg_tsv.register("REF_ALT", Setter::VerifyRefAlt);
    let mut hap_tsv = Tsv::new("HAPS");
    hap_tsv.register("HAPS", Setter::Haps);

    let mut header = new_header(s, &chrom, Gp::Omitted)?;
    add_samples(&mut header, &sample_fname, 1)?;
    let mut writer = open_output(s, &header, "haplegendsample_to_vcf")?;
    let mut ctx = Ctx::new(&header, s);
    let mut rec = Record::new();
    let mut total = 0u64;
    let result = (|| -> anyhow::Result<()> {
        loop {
            rec.clear();
            total += 1;
            if !leg_tsv.parse(&mut ctx, &mut rec, &line)? {
                bail!(
                    "Error occurred while parsing {leg_fname}: {}",
                    rest(&line, 0)
                );
            }
            if !hap_fh.next_line(&mut line)? {
                bail!("Different number of records in {leg_fname} and {hap_fname}?");
            }
            if !hap_tsv.parse(&mut ctx, &mut rec, &line)? {
                bail!(
                    "Error occurred while parsing {hap_fname}: {}",
                    rest(&line, 0)
                );
            }
            write_built(&mut writer, &header, &mut rec, s, "haplegendsample_to_vcf")?;
            if !leg_fh.next_line(&mut line)? {
                if hap_fh.next_line(&mut line)? {
                    bail!("Different number of records in {leg_fname} and {hap_fname}?");
                }
                return Ok(());
            }
        }
    })();
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    finish_output(writer, s, "haplegendsample_to_vcf", CloseMessage::Short)?;
    eprintln!("Number of processed rows: \t{total}");
    Ok(())
}

/// SHAPEIT hap/sample to VCF (bcftools' `hapsample_to_vcf`).
pub(super) fn hapsample_to_vcf(s: &Settings) -> anyhow::Result<()> {
    let names = input_names(&s.infname, &[".hap.gz", ".samples"]);
    let (hap_fname, sample_fname) = (&names[0], &names[1]);
    let mut hap_fh = open_lines(hap_fname)?;
    let mut line = Vec::new();
    if !hap_fh.next_line(&mut line)? {
        bail!("Empty file: {hap_fname}");
    }
    // The chromosome, from CHROM:POS_REF_ALT in the first column (--vcf-ids) or the second.
    let fail = || {
        anyhow::anyhow!(
            "Could not determine CHROM in {hap_fname}: {}",
            rest(&line, 0)
        )
    };
    let mut sb = 0;
    let mut se = token_end(&line, 0);
    if se == line.len() {
        return Err(fail());
    }
    if !s.vcf_ids {
        se += 1;
        sb = se;
        se = token_end(&line, se);
        if se == line.len() {
            return Err(fail());
        }
        if !line[sb..].contains(&b':') {
            bail!(
                "Could not determine CHROM in the second column of {hap_fname}: {}",
                rest(&line, 0)
            );
        }
    }
    let Some(sc) = line[sb..]
        .iter()
        .position(|&c| c == b':')
        .map(|i| i + sb)
        .filter(|&c| c <= se)
    else {
        return Err(fail());
    };
    let chrom = line[sb..sc].to_vec();

    let mut header = new_header(s, &chrom, Gp::Omitted)?;
    let mut tsv = if s.vcf_ids {
        let mut t = Tsv::new("CHROM_POS_REF_ALT,ID,POS,REF_ALT,HAPS");
        t.register("ID", Setter::Id);
        t
    } else {
        let mut t = Tsv::new("CHROM,CHROM_POS_REF_ALT,POS,REF_ALT,HAPS");
        t.register("CHROM", Setter::ChromPosRefAltOrChrom);
        t
    };
    tsv.register("CHROM_POS_REF_ALT", Setter::ChromPosRefAlt);
    tsv.register("POS", Setter::VerifyPos);
    tsv.register("REF_ALT", Setter::VerifyRefAlt);
    tsv.register("HAPS", Setter::Haps);

    add_samples(&mut header, sample_fname, 2)?;
    let mut writer = open_output(s, &header, "hapsample_to_vcf")?;
    let mut ctx = Ctx::new(&header, s);
    let mut rec = Record::new();
    let mut total = 0u64;
    let result = (|| -> anyhow::Result<()> {
        loop {
            rec.clear();
            total += 1;
            if !tsv.parse(&mut ctx, &mut rec, &line)? {
                bail!("Error occurred while parsing: {}", rest(&line, 0));
            }
            write_built(&mut writer, &header, &mut rec, s, "hapsample_to_vcf")?;
            if !hap_fh.next_line(&mut line)? {
                return Ok(());
            }
        }
    })();
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    finish_output(writer, s, "hapsample_to_vcf", CloseMessage::Short)?;
    eprintln!("Number of processed rows: \t{total}");
    Ok(())
}

/// Tab- or space-separated sites (and optionally genotypes) to VCF (bcftools' `tsv_to_vcf`).
pub(super) fn tsv_to_vcf(s: &Settings) -> anyhow::Result<()> {
    let Some(ref_fname) = &s.fasta_ref else {
        bail!("--tsv2vcf requires the --fasta-ref option");
    };
    let mut reference = Reference::open(Path::new(ref_fname))
        .map_err(|_| anyhow::anyhow!("Could not load the reference {ref_fname}"))?;

    let mut header = Header::new_for_writing();
    for (name, len) in reference.sequences() {
        let mut contig = b"##contig=<ID=".to_vec();
        contig.extend_from_slice(&name);
        contig.extend_from_slice(format!(",length={}>", len as i32).as_bytes());
        header.append(&contig)?;
    }
    header.append(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    if s.record_cmd_line {
        version::append_version(&mut header, "convert")?;
    }
    if let Some((list, source)) = &s.samples {
        let Ok(names) = readlist::read_list(list.as_bytes(), *source) else {
            bail!("Could not parse {list}");
        };
        for name in names {
            // bcftools ignores failures (an empty or duplicate name is reported and skipped).
            let _ = header.add_sample(&name);
        }
    }

    let mut writer = open_output(s, &header, "tsv_to_vcf")?;

    let mut tsv = Tsv::new(s.columns.as_deref().unwrap_or("ID,CHROM,POS,AA"));
    if !tsv.register("CHROM", Setter::Chrom) {
        bail!("Expected CHROM column");
    }
    if !tsv.register("POS", Setter::Pos) {
        bail!("Expected POS column");
    }
    if !tsv.register("ID", Setter::Id) && s.columns.is_none() {
        bail!("Expected ID column");
    }
    if !tsv.register("AA", Setter::Aa) {
        if s.samples.is_some() {
            bail!("Expected AA column with -s/-S");
        }
        if !tsv.register("REF", Setter::Ref) || !tsv.register("ALT", Setter::Alt) {
            bail!("Expected REF and ALT columns when AA was not given");
        }
    }

    let mut ctx = Ctx::new(&header, s);
    ctx.reference = Some(&mut reference);
    let mut rec = Record::new();
    let (mut total, mut skipped, mut written) = (0u64, 0u64, 0u64);
    let mut in_fh = open_lines(&s.infname)?;
    let mut line = Vec::new();
    let result = (|| -> anyhow::Result<()> {
        while in_fh.next_line(&mut line)? {
            if line[0] == b'#' {
                continue;
            }
            rec.clear();
            total += 1;
            if tsv.parse(&mut ctx, &mut rec, &line)? {
                write_built(&mut writer, &header, &mut rec, s, "tsv_to_vcf")?;
                written += 1;
            } else {
                skipped += 1;
            }
        }
        Ok(())
    })();
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    finish_output(writer, s, "tsv_to_vcf", CloseMessage::Named)?;
    eprintln!("Rows total: \t{total}");
    eprintln!("Rows skipped: \t{skipped}");
    eprintln!("Sites written: \t{written}");
    if s.samples.is_some() {
        let c = &ctx.counts;
        eprintln!("Missing GTs: \t{}", c.missing);
        eprintln!("Hom RR: \t{}", c.hom_rr);
        eprintln!("Het RA: \t{}", c.het_ra);
        eprintln!("Hom AA: \t{}", c.hom_aa);
        eprintln!("Het AA: \t{}", c.het_aa);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn columns_parse_as_bcftools() {
        let mut t = Tsv::new("ID,-,pos,AA,");
        assert_eq!(t.cols.len(), 4);
        assert!(t.cols[1].0.is_none());
        assert!(t.register("POS", Setter::Pos));
        assert!(!t.register("REF", Setter::Ref));
        assert_eq!(t.cols[2].1, Some(Setter::Pos));
    }

    #[test]
    fn next_moves_over_whitespace() {
        let line = b"a  bb\tc";
        let mut t = Tsv::new("X");
        assert!(t.next(line));
        assert_eq!(&line[t.ss..t.se], b"a");
        assert!(t.next(line));
        assert_eq!(&line[t.ss..t.se], b"bb");
        assert!(t.next(line));
        assert_eq!(&line[t.ss..t.se], b"c");
        assert!(!t.next(line));
    }
}
