//! bcftools' assignment of reads to output samples (`bam_sample.c`): by read group (`RG`)
//! and the header's `@RG SM` fields, restricted and renamed by `-s/-S` and `-G`, or one sample
//! per file with `--ignore-RG`.

use std::collections::HashMap;

use anyhow::bail;
use pasta::{pileup::Read, variant::readlist};

/// One input file's read groups (`file_t`).
#[derive(Debug, Default)]
struct FileSamples {
    name: String,
    /// Read group ID to output sample index, `-1` for a read group that is ignored.
    rg2idx: HashMap<Vec<u8>, i32>,
    has_rg_map: bool,
    /// The sample of every read, when all read groups are one sample.
    default_idx: i32,
}

/// `bam_smpl_t`
#[derive(Debug, Default)]
pub struct Samples {
    files: Vec<FileSamples>,
    ignore_rg: bool,
    /// Output sample names, in order.
    names: Vec<Vec<u8>>,
    name2idx: HashMap<Vec<u8>, i32>,
    /// `-s/-S`: input sample name to output name.
    sample_list: Option<HashMap<Vec<u8>, Vec<u8>>>,
    sample_logic: bool,
    /// `-G`: read group (or `RG\tFILE`, `*\tFILE`) to sample name, `\t` to keep the name.
    rg_list: Option<HashMap<Vec<u8>, Vec<u8>>>,
    rg_logic: bool,
}

/// Splits a `-s/-G` list row into whitespace-separated fields, honouring `\` escapes, as
/// bcftools reads them (at most `n` fields).
fn split_escaped(row: &[u8], n: usize) -> Vec<Vec<u8>> {
    let mut fields = Vec::new();
    let mut i = 0;
    for f in 0..n {
        if f > 0 {
            if i >= row.len() {
                break;
            }
            while i < row.len() && row[i].is_ascii_whitespace() {
                i += 1;
            }
        }
        let mut field = Vec::new();
        let mut escaped = false;
        while i < row.len() {
            let c = row[i];
            if c == b'\\' && !escaped {
                escaped = true;
                i += 1;
                continue;
            }
            if c.is_ascii_whitespace() && !escaped {
                break;
            }
            field.push(c);
            escaped = false;
            i += 1;
        }
        fields.push(field);
    }
    fields
}

impl Samples {
    /// `bam_smpl_ignore_readgroups`
    pub fn ignore_read_groups(&mut self) {
        self.ignore_rg = true;
    }

    /// `bam_smpl_add_samples` (`-s` or, with `is_file`, `-S`).
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read.
    pub fn add_samples(&mut self, list: &str, source: readlist::ListSource) -> anyhow::Result<()> {
        let list = match list.strip_prefix('^') {
            Some(rest) => rest,
            None => {
                self.sample_logic = true;
                list
            }
        };
        let rows = readlist::read_list(list.as_bytes(), source)
            .map_err(|_| anyhow::anyhow!("Could not read samples: {list}"))?;
        if rows.is_empty() {
            return Ok(());
        }
        let mut map = HashMap::new();
        for row in rows {
            let mut f = split_escaped(&row, 2);
            let ori = f.remove(0);
            let ren = f
                .pop()
                .filter(|r| !r.is_empty())
                .unwrap_or_else(|| ori.clone());
            map.insert(ori, ren);
        }
        self.sample_list = Some(map);
        Ok(())
    }

    /// `bam_smpl_add_readgroups` (`-G`).
    ///
    /// # Errors
    ///
    /// Returns an error if a read group is given two samples.
    pub fn add_read_groups(&mut self, list: &str) -> anyhow::Result<()> {
        let list = match list.strip_prefix('^') {
            Some(rest) => rest,
            None => {
                self.rg_logic = true;
                list
            }
        };
        let Ok(rows) = readlist::read_list(list.as_bytes(), readlist::ListSource::File) else {
            return Ok(());
        };
        if rows.is_empty() {
            return Ok(());
        }
        let mut map: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        for row in rows {
            let f = split_escaped(&row, 3);
            let mut fld1 = f.first().cloned().unwrap_or_default();
            let mut fld2 = f.get(1).cloned().unwrap_or_default();
            if let Some(fld3) = f.get(2).filter(|s| !s.is_empty()) {
                // ID FILE SAMPLE
                fld1.push(b'\t');
                fld1.extend_from_slice(&fld2);
                fld2 = fld3.clone();
            }
            let value = if fld2.is_empty() {
                b"\t".to_vec()
            } else {
                fld2
            };
            match map.get(&fld1) {
                None => {
                    map.insert(fld1, value);
                }
                Some(old) if *old != value => bail!(
                    "Error: The read group \"{}\" was assigned to two different samples: \"{}\" and \"{}\"",
                    String::from_utf8_lossy(&fld1),
                    String::from_utf8_lossy(old),
                    String::from_utf8_lossy(&value)
                ),
                Some(_) => {}
            }
        }
        self.rg_list = Some(map);
        Ok(())
    }

    /// `bsmpl_add_readgroup`
    fn add_read_group(&mut self, file: usize, rg_id: &[u8], smpl_name: Option<&[u8]>) {
        let mut ismpl = -1;
        if let Some(name) = smpl_name {
            ismpl = match self.name2idx.get(name) {
                Some(&i) => i,
                None => {
                    self.names.push(name.to_vec());
                    let i = self.names.len() as i32 - 1;
                    self.name2idx.insert(name.to_vec(), i);
                    i
                }
            };
        }
        let file = &mut self.files[file];
        if rg_id == b"*" {
            file.default_idx = ismpl;
            return;
        }
        file.has_rg_map = true;
        file.rg2idx.entry(rg_id.to_vec()).or_insert(ismpl);
    }

    /// `bsmpl_keep_readgroup`: whether to keep a read group, and the sample to rename it to.
    fn keep_read_group(&self, file: usize, rg_id: &[u8]) -> (bool, Option<Vec<u8>>) {
        let rg_list = self.rg_list.as_ref().expect("called with -G");
        let fname = self.files[file].name.as_bytes();
        let mut key = rg_id.to_vec();
        let mut rg_smpl = rg_list.get(&key);
        if rg_smpl.is_none() {
            key.push(b'\t');
            key.extend_from_slice(fname);
            rg_smpl = rg_list.get(&key);
        }
        if rg_smpl.is_none() {
            let mut key = b"*\t".to_vec();
            key.extend_from_slice(fname);
            rg_smpl = rg_list.get(&key);
        }
        if rg_smpl.is_none() && self.rg_logic {
            return (false, None);
        }
        if rg_smpl.is_some() && !self.rg_logic {
            return (false, None);
        }
        let rename = rg_smpl.filter(|s| s.first() != Some(&b'\t')).cloned();
        (true, rename)
    }

    /// `bam_smpl_add_bam`: registers a file and its read groups; `None` when the file has no
    /// usable read group and is skipped.
    ///
    /// # Errors
    ///
    /// Returns an error for the read group IDs `*` and `?`.
    pub fn add_bam(
        &mut self,
        header_text: Option<&[u8]>,
        fname: &str,
    ) -> anyhow::Result<Option<usize>> {
        self.files.push(FileSamples {
            name: fname.to_owned(),
            default_idx: -1,
            ..FileSamples::default()
        });
        let ifile = self.files.len() - 1;
        let hdr = match header_text {
            Some(h) if !self.ignore_rg => h,
            _ => {
                self.add_read_group(ifile, b"*", Some(fname.as_bytes()));
                return Ok(Some(ifile));
            }
        };

        let mut bam_smpls: Vec<Vec<u8>> = Vec::new();
        let mut first_smpl = -1;
        let mut nskipped = 0;
        let find = |from: usize, pat: &[u8]| -> Option<usize> {
            memchr::memmem::find(&hdr[from..], pat).map(|i| i + from)
        };
        let mut p = Some(0usize);
        while let Some(pp) = p {
            let Some(q) = find(pp, b"@RG") else { break };
            let eol = memchr::memchr(b'\n', &hdr[q + 3..]).map(|i| i + q + 3);
            if q > 0 && hdr[q - 1] != b'\n' {
                p = eol;
                continue;
            }
            let pp = q + 3;
            let id = find(pp, b"\tID:").map(|i| i + 4);
            let sm = find(pp, b"\tSM:").map(|i| i + 4);
            let (Some(qi), Some(ri)) = (id, sm) else {
                break;
            };
            let field_end = |s: usize| {
                hdr[s..]
                    .iter()
                    .position(|&c| c == b'\t' || c == b'\n' || c == 0)
                    .map_or(hdr.len(), |i| i + s)
            };
            let rg_id = &hdr[qi..field_end(qi)];
            let mut smpl: Vec<u8> = hdr[ri..field_end(ri)].to_vec();
            if rg_id == b"*" || rg_id == b"?" {
                bail!(
                    "Error: the read group IDs \"*\" and \"?\" have a special meaning in the mpileup code. Please fix the code or the bam: {fname}"
                );
            }
            let mut accept = true;
            if let Some(list) = &self.sample_list {
                let name = list.get(&smpl);
                if !self.sample_logic {
                    accept = name.is_none();
                } else if let Some(name) = name {
                    smpl = name.clone();
                } else {
                    accept = false;
                }
            }
            if accept && self.rg_list.is_some() {
                let (keep, rename) = self.keep_read_group(ifile, rg_id);
                accept = keep;
                if let Some(r) = rename {
                    smpl = r;
                }
            }
            if accept {
                self.add_read_group(ifile, rg_id, Some(&smpl));
            } else {
                self.add_read_group(ifile, rg_id, None);
                nskipped += 1;
            }
            if first_smpl < 0
                && let Some(&i) = self.name2idx.get(&smpl)
            {
                first_smpl = i;
            }
            if !bam_smpls.contains(&smpl) {
                bam_smpls.push(smpl);
            }
            p = eol;
        }
        let nsmpls = bam_smpls.len();

        let mut smpl_name: Option<Vec<u8>> = None;
        let mut accept_null_rg = true;
        if self.rg_list.is_some() {
            let (keep, rename) = self.keep_read_group(ifile, b"?");
            if !keep {
                accept_null_rg = false;
            }
            smpl_name = rename;
        }
        if self.sample_list.is_some() && first_smpl == -1 {
            accept_null_rg = false;
        }
        if !accept_null_rg && first_smpl == -1 {
            self.files.pop();
            return Ok(None);
        }
        if !accept_null_rg {
            return Ok(Some(ifile));
        }
        if nsmpls == 1 && nskipped == 0 {
            self.files[ifile].default_idx = first_smpl;
            return Ok(Some(ifile));
        }
        let name = smpl_name.unwrap_or_else(|| {
            if first_smpl == -1 {
                fname.as_bytes().to_vec()
            } else {
                self.names[first_smpl as usize].clone()
            }
        });
        self.add_read_group(ifile, b"?", Some(&name));
        Ok(Some(ifile))
    }

    /// The output sample names.
    pub fn names(&self) -> &[Vec<u8>] {
        &self.names
    }

    /// `bam_smpl_get_sample_id`: the output sample of a read, or -1 to skip it.
    pub fn sample_id(&self, file: usize, read: &Read) -> i32 {
        let file = &self.files[file];
        if file.default_idx >= 0 {
            return file.default_idx;
        }
        if !file.has_rg_map {
            // khash_str2int_get on a NULL hash fails.
            return -1;
        }
        let rg: &[u8] = match read.tag(*b"RG") {
            Some((_, value)) => &value[..memchr::memchr(0, value).unwrap_or(value.len())],
            None => b"?",
        };
        if let Some(&i) = file.rg2idx.get(rg) {
            return i;
        }
        file.rg2idx.get(&b"?"[..]).copied().unwrap_or(-1)
    }
}
