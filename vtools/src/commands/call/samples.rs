//! Samples, sexes, PED families, and `-G` groups (bcftools' `set_samples`,
//! `parse_ped_samples`, and `init_sample_groups`).

use std::collections::HashMap;

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType,
    readlist::{self, ListSource},
};

use super::{Shared, mcall::SampleGroup, ploidy::Ploidy};

/// The samples to call and their sexes.
#[derive(Clone, Debug, Default)]
pub struct SampleSet {
    /// The number of samples called.
    pub nsamples: usize,
    /// The names of the samples given with `-s`/`-S` (empty without them).
    pub names: Vec<Vec<u8>>,
    /// For each called sample, its index in the input, when a subset or reordering is needed.
    pub samples_map: Option<Vec<i32>>,
    /// Each sample's sex ID, or minus its ploidy for an undeclared sex 1 or 2.
    pub sample2sex: Vec<i32>,
    /// Trio families: father, mother, child.
    pub families: Vec<[i32; 3]>,
}

/// One PED-derived line: `name sex`, and the families.
struct Ped {
    lines: Vec<Vec<u8>>,
    families: Vec<[i32; 3]>,
}

/// bcftools' `parse_ped_samples`: `None` when the lines are not a PED file.
fn parse_ped(ploidy: &Ploidy, vals: &[Vec<u8>]) -> anyhow::Result<Option<Ped>> {
    let mut name2idx: HashMap<Vec<u8>, i32> = HashMap::new();
    let mut lines: Vec<Vec<u8>> = Vec::new();
    let mut families = Vec::new();
    let mut add = |lines: &mut Vec<Vec<u8>>, name: &[u8], sex: &[u8]| -> i32 {
        if let Some(&i) = name2idx.get(name) {
            return i;
        }
        let mut l = name.to_vec();
        l.push(b' ');
        l.extend_from_slice(sex);
        lines.push(l);
        let i = lines.len() as i32 - 1;
        name2idx.insert(name.to_vec(), i);
        i
    };
    for val in vals {
        // Columns: family, sample, father, mother, sex (and more). Each ends at the first
        // whitespace after it.
        let ws = |c: u8| c.is_ascii_whitespace();
        let mut starts = Vec::with_capacity(5);
        let mut i = 0;
        while i < val.len() && starts.len() < 5 {
            if ws(val[i]) {
                i += 1;
                while i < val.len() && ws(val[i]) {
                    i += 1;
                }
                starts.push(i);
                continue;
            }
            i += 1;
        }
        if starts.len() < 4 {
            return Ok(None);
        }
        let col = |k: usize| {
            let s = starts[k];
            let e = val[s..]
                .iter()
                .position(|&c| ws(c))
                .map_or(val.len(), |p| s + p);
            &val[s..e]
        };
        let (name, father, mother) = (col(0), col(1), col(2));
        let mut sex = col(3);
        if ploidy.sex_id(sex).is_none() {
            if sex == b"1" && ploidy.sex_id(b"M").is_some() {
                sex = b"M";
            } else if sex == b"2" && ploidy.sex_id(b"F").is_some() {
                sex = b"F";
            } else {
                bail!(
                    "[E::parse_ped_samples] The sex \"{}\" has not been declared in --ploidy/--ploidy-file",
                    String::from_utf8_lossy(sex)
                );
            }
        }
        let child = add(&mut lines, name, sex);
        if father != b"0" && mother != b"0" {
            let f = add(&mut lines, father, b"M");
            let m = add(&mut lines, mother, b"F");
            families.push([f, m, child]);
        }
    }
    Ok(Some(Ped { lines, families }))
}

/// bcftools' `set_samples`: reads `-s`/`-S` (a sample list, a list with sexes or ploidies,
/// or a PED file; `^` excludes).
///
/// # Errors
///
/// Returns bcftools' errors for unreadable lists and undeclared sexes.
pub fn set_samples(
    hdr: &Header,
    ploidy: &mut Ploidy,
    spec: &str,
    source: ListSource,
) -> anyhow::Result<SampleSet> {
    let (negate, name) = match spec.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, spec),
    };
    let Ok(mut lines) = readlist::read_list(name.as_bytes(), source) else {
        bail!("Could not read the file: {name}");
    };
    let mut families = Vec::new();
    match parse_ped(ploidy, &lines)? {
        Some(ped) => {
            lines = ped.lines;
            families = ped.families;
        }
        None if source == ListSource::File => eprintln!("Note: could not parse as PED: {name}"),
        None => {}
    }

    let mut sel = Selection::new(hdr.sample_count(), ploidy.nsex() as i32 - 1);
    if negate {
        sel.exclude(hdr, &lines)?;
    } else {
        sel.include(hdr, ploidy, &lines)?;
    }

    // bcftools maps the PED-line indices of the families as if they were header indices.
    // Its `assert(nmiss==0 || nmiss==3)` on families with members missing from the VCF is
    // compiled out in release builds; such families only matter to -C trio, which then
    // stops on the trio count.
    for fam in &mut families {
        for s in fam.iter_mut() {
            *s = sel.old2new.get(*s as usize).copied().unwrap_or(-1);
        }
    }
    let nsmpl = sel.names.len();
    sel.samples_map.truncate(nsmpl);
    Ok(SampleSet {
        nsamples: nsmpl,
        names: sel.names,
        samples_map: sel.map_needed.then_some(sel.samples_map),
        sample2sex: sel.sample2sex,
        families,
    })
}

/// The first word of a `-s`/`-S` line and the index where it ends; `None` for a comment.
fn first_word(line: &[u8]) -> anyhow::Result<Option<(&[u8], usize)>> {
    let ws = |c: u8| c.is_ascii_whitespace();
    let s = line.iter().position(|&c| !ws(c)).unwrap_or(line.len());
    if s == line.len() {
        bail!("Could not parse: {}", String::from_utf8_lossy(line));
    }
    if line[s] == b'#' {
        return Ok(None);
    }
    let e = line[s..]
        .iter()
        .position(|&c| ws(c))
        .map_or(line.len(), |p| s + p);
    Ok(Some((&line[s..e], e)))
}

/// The samples `-s`/`-S` keep, as `set_samples` works them out.
struct Selection {
    /// The names of the kept samples, in output order.
    names: Vec<Vec<u8>>,
    /// For each kept sample, its index in the header (room for all of them).
    samples_map: Vec<i32>,
    /// For each header sample, its index on output, or -1.
    old2new: Vec<i32>,
    /// Each sample's sex, as [`SampleSet::sample2sex`] (room for all of them).
    sample2sex: Vec<i32>,
    /// Whether the kept samples are not all the header's samples in order.
    map_needed: bool,
}

impl Selection {
    /// Nothing kept yet, of `nhdr` samples of default sex `dflt_sex`.
    fn new(nhdr: usize, dflt_sex: i32) -> Self {
        Self {
            names: Vec::new(),
            samples_map: vec![0; nhdr],
            old2new: vec![-1; nhdr],
            sample2sex: vec![dflt_sex; nhdr],
            map_needed: false,
        }
    }

    /// Keeps the listed samples in the order listed, with the sex or ploidy after each name.
    fn include(
        &mut self,
        hdr: &Header,
        ploidy: &mut Ploidy,
        lines: &[Vec<u8>],
    ) -> anyhow::Result<()> {
        for line in lines {
            let Some((sample, end)) = first_word(line)? else {
                continue;
            };
            let Some(ismpl) = hdr.sample_id(sample) else {
                eprintln!(
                    "Warning: No such sample in the VCF: {}",
                    String::from_utf8_lossy(sample)
                );
                continue;
            };
            let ismpl = ismpl as usize;
            if self.old2new[ismpl] != -1 {
                eprintln!(
                    "Warning: The sample is listed multiple times: {}",
                    String::from_utf8_lossy(sample)
                );
                continue;
            }
            let nsmpl = self.names.len();
            let mut p = (end + 1).min(line.len());
            while p < line.len() && line[p].is_ascii_whitespace() {
                p += 1;
            }
            // The sex is the rest of the line (bcftools does not cut it at whitespace).
            let sex: &[u8] = if p == line.len() { b"2" } else { &line[p..] };
            self.sample2sex[nsmpl] = match ploidy.sex_id(sex) {
                Some(_) => ploidy.add_sex(sex) as i32,
                None if sex.len() == 1 && matches!(sex[0], b'0' | b'1' | b'2') => {
                    -i32::from(sex[0] - b'0')
                }
                None => bail!(
                    "[E::set_samples] The sex \"{}\" has not been declared in --ploidy/--ploidy-file",
                    String::from_utf8_lossy(sex)
                ),
            };
            if ismpl != nsmpl {
                self.map_needed = true;
            }
            self.samples_map[nsmpl] = ismpl as i32;
            self.old2new[ismpl] = nsmpl as i32;
            self.names.push(sample.to_vec());
        }
        if self.names.len() != hdr.sample_count() {
            self.map_needed = true;
        }
        Ok(())
    }

    /// Keeps the samples not listed (`^`), in header order.
    fn exclude(&mut self, hdr: &Header, lines: &[Vec<u8>]) -> anyhow::Result<()> {
        let mut keep = vec![true; hdr.sample_count()];
        for line in lines {
            let Some((sample, _)) = first_word(line)? else {
                continue;
            };
            match hdr.sample_id(sample) {
                Some(i) => keep[i as usize] = false,
                None => eprintln!(
                    "Warning: No such sample in the VCF: {}",
                    String::from_utf8_lossy(sample)
                ),
            }
        }
        for (i, k) in keep.iter().enumerate() {
            if !k {
                continue;
            }
            let nsmpl = self.names.len();
            self.names.push(hdr.samples()[i].to_vec());
            self.samples_map[nsmpl] = i as i32;
            self.old2new[i] = nsmpl as i32;
        }
        self.map_needed = true;
        Ok(())
    }
}

/// bcftools' `init_sample_groups`: one group of all samples, one group per sample (`-G -`),
/// or the groups of a file; and the FORMAT tag giving their allele frequencies.
///
/// # Errors
///
/// Returns bcftools' errors for missing tags and malformed or incomplete group files.
pub fn sample_groups(
    hdr: &Header,
    groups: Option<&str>,
    tag: Option<&str>,
) -> anyhow::Result<(Vec<SampleGroup>, Option<Vec<u8>>)> {
    let nsmpl = hdr.sample_count();
    let Some(groups) = groups else {
        let all = (0..nsmpl as u32).collect();
        return Ok((vec![SampleGroup::new(all)], None));
    };
    let is_fmt = |t: &[u8]| {
        hdr.id(t)
            .is_some_and(|id| hdr.id_exists(HeaderRecordType::Format, id))
    };
    let tag: Vec<u8> = match tag {
        Some(t) => {
            if hdr.id(t.as_bytes()).is_none() {
                bail!("No such tag \"{t}\"");
            }
            if !is_fmt(t.as_bytes()) {
                bail!("No such FORMAT tag \"{t}\"");
            }
            t.as_bytes().to_vec()
        }
        None if is_fmt(b"QS") => b"QS".to_vec(),
        None if is_fmt(b"AD") => b"AD".to_vec(),
        None => bail!(
            "Error: neither \"AD\" nor \"QS\" FORMAT tag exists and no alternative given with -G"
        ),
    };
    if groups == "-" {
        let g = (0..nsmpl as u32)
            .map(|i| SampleGroup::new(vec![i]))
            .collect();
        return Ok((g, Some(tag)));
    }
    let Ok(lines) = readlist::read_list(groups.as_bytes(), ListSource::File) else {
        bail!("Could not read the file: {groups}");
    };
    let mut smpl2grp = vec![0usize; nsmpl];
    let mut grp2idx: HashMap<Vec<u8>, usize> = HashMap::new();
    let ws = |c: u8| c.is_ascii_whitespace();
    for line in &lines {
        let err = || {
            anyhow::anyhow!(
                "Could not parse the line in {groups}, expected a sample name followed by tab and a population name: {}",
                String::from_utf8_lossy(line)
            )
        };
        let e = line.iter().position(|&c| ws(c)).ok_or_else(err)?;
        let g = line[e..].iter().position(|&c| !ws(c)).ok_or_else(err)? + e;
        let name = &line[..e];
        let Some(ismpl) = hdr.sample_id(name) else {
            continue;
        };
        let ismpl = ismpl as usize;
        if smpl2grp[ismpl] != 0 {
            bail!(
                "Error: the sample \"{}\" is listed twice in {groups}",
                String::from_utf8_lossy(name)
            );
        }
        // bcftools keys groups on the name without its first character, merging groups
        // such as "A" and "B"; the whole name is used here.
        let key = line[g..].to_vec();
        let n = grp2idx.len();
        let igrp = *grp2idx.entry(key).or_insert(n);
        smpl2grp[ismpl] = igrp + 1;
    }
    if grp2idx.is_empty() {
        bail!("Could not parse the file, no matching samples found: {groups}");
    }
    let mut out: Vec<SampleGroup> = (0..grp2idx.len())
        .map(|_| SampleGroup::new(Vec::new()))
        .collect();
    for (i, &g) in smpl2grp.iter().enumerate() {
        if g == 0 {
            bail!(
                "Error: The sample \"{}\" is not listed in {groups}",
                String::from_utf8_lossy(&hdr.samples()[i])
            );
        }
        out[g - 1].smpl.push(i as u32);
    }
    Ok((out, Some(tag)))
}

/// The checks of `mcall_init_trios`: the `-F` tags, and family ploidies.
///
/// # Errors
///
/// Returns bcftools' errors for undefined `-F` tags and unsupported family ploidies.
pub fn check_trios(
    hdr: &Header,
    shared: &Shared,
    set: &SampleSet,
    ploidy: &[u8],
) -> anyhow::Result<()> {
    if let (Some(an), Some(ac)) = (&shared.prior_an, &shared.prior_ac) {
        for t in [an, ac] {
            let name = String::from_utf8_lossy(t);
            let Some(id) = hdr.id(t) else {
                bail!("No such tag \"{name}\"");
            };
            if !hdr.id_exists(HeaderRecordType::Format, id) {
                bail!("No such FORMAT tag \"{name}\"");
            }
        }
    }
    for fam in &set.families {
        let p: Vec<i32> = fam
            .iter()
            .map(|&s| ploidy.get(s as usize).map_or(0, |&p| i32::from(p)))
            .collect();
        let incorrect = || anyhow::anyhow!("Incorrect ploidy: {} {} {}", p[0], p[1], p[2]);
        if p[0] == 2 {
            if p[1] != 2 || p[2] != 2 {
                return Err(incorrect());
            }
            continue;
        }
        if p[0] != 1 || p[1] == 1 {
            return Err(incorrect());
        }
        if p[1] == 2 {
            if p[2] == 0 {
                return Err(incorrect());
            }
        } else if p[2] == 2 {
            return Err(incorrect());
        }
    }
    Ok(())
}
