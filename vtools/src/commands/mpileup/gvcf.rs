//! gVCF blocks (`-g/--gvcf`): a port of bcftools' `gvcf.c`, which merges runs of reference
//! sites with similar per-sample depth into one record with `INFO/END` and `INFO/MIN_DP`.

use anyhow::bail;
use pasta::variant::{Header, Record};

/// `gvcf_t`
#[derive(Debug)]
pub struct Gvcf {
    dp_range: Vec<i32>,
    prev_range: usize,
    dp: Vec<i32>,
    pl: Vec<i32>,
    npl: i32,
    tmp: Vec<i32>,
    gts: Vec<i32>,
    qsum: Vec<f32>,
    rid: i32,
    start: i64,
    end: i64,
    min_dp: i32,
    als: Vec<u8>,
    line: Record,
}

impl Gvcf {
    /// `gvcf_init`: parses the comma-separated depth thresholds.
    ///
    /// # Errors
    ///
    /// Returns an error if the list cannot be parsed.
    pub fn new(dp_ranges: &str) -> anyhow::Result<Self> {
        let err = || anyhow::anyhow!("Could not parse: --gvcf {dp_ranges}");
        let bytes = dp_ranges.as_bytes();
        let mut dp_range = Vec::new();
        let mut ss = 0;
        while ss < bytes.len() {
            let (v, se) = pasta::variant::number::strtol(bytes, ss);
            if se == ss {
                return Err(err());
            }
            dp_range.push(v as i32);
            if se < bytes.len() && bytes[se] == b',' && se + 1 < bytes.len() {
                ss = se + 1;
                continue;
            } else if se == bytes.len() {
                break;
            }
            return Err(err());
        }
        Ok(Self {
            dp_range,
            prev_range: 0,
            dp: Vec::new(),
            pl: Vec::new(),
            npl: 0,
            tmp: Vec::new(),
            gts: Vec::new(),
            qsum: Vec::new(),
            rid: -1,
            start: 0,
            end: 0,
            min_dp: 0,
            als: Vec::new(),
            line: Record::new(),
        })
    }

    /// `gvcf_update_header`
    ///
    /// # Errors
    ///
    /// Returns an error if the header lines cannot be added.
    pub fn update_header(hdr: &mut Header) -> anyhow::Result<()> {
        hdr.append(b"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">")?;
        hdr.append(b"##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum per-sample depth in this gVCF block\">")?;
        Ok(())
    }

    /// Adds a record, calling `emit` with each record to write (a block it ends). Returns
    /// whether `rec` itself is to be written after them.
    ///
    /// # Errors
    ///
    /// Returns `emit`'s errors, or an error for a PL count that does not match.
    pub fn write<F>(&mut self, hdr: &Header, rec: &mut Record, emit: F) -> anyhow::Result<bool>
    where
        F: FnMut(&mut Record) -> anyhow::Result<()>,
    {
        // `flush_bcf_records`: a site without an alternate allele, or with only `<*>`.
        let n = rec.allele_count();
        let is_ref = n == 1 || (n == 2 && rec.allele(1).starts_with(b"<*>"));
        self.gvcf_write(hdr, Some(rec), is_ref, emit)
    }

    /// Writes the pending block, if any, through `emit`.
    ///
    /// # Errors
    ///
    /// Returns `emit`'s errors.
    pub fn flush<F>(&mut self, hdr: &Header, emit: F) -> anyhow::Result<()>
    where
        F: FnMut(&mut Record) -> anyhow::Result<()>,
    {
        self.gvcf_write(hdr, None, false, emit)?;
        Ok(())
    }

    /// `gvcf_write`: adds a record (`None` to flush), calling `emit` with each record to write.
    /// Returns whether `rec` itself is to be written after the flushed block.
    fn gvcf_write<F>(
        &mut self,
        hdr: &Header,
        rec: Option<&mut Record>,
        is_ref: bool,
        mut emit: F,
    ) -> anyhow::Result<bool>
    where
        F: FnMut(&mut Record) -> anyhow::Result<()>,
    {
        let nsmpl = hdr.sample_count();
        let mut can_collapse = is_ref;
        let mut dp_range = 0usize;
        let mut min_dp = 0i32;
        if rec.is_none() && self.prev_range == 0 {
            return Ok(false);
        }
        let mut needs_flush = !can_collapse;
        let mut rec = rec;
        if let Some(r) = rec.as_deref_mut()
            && can_collapse
        {
            match r.get_format_i32(hdr, b"DP", &mut self.tmp) {
                Ok(n) if n == nsmpl => {
                    min_dp = *self.tmp[..nsmpl].iter().min().unwrap_or(&0);
                    let mut i = 0;
                    while i < self.dp_range.len() {
                        if min_dp < self.dp_range[i] {
                            break;
                        }
                        i += 1;
                    }
                    dp_range = i;
                    if dp_range == 0 {
                        needs_flush = true;
                        can_collapse = false;
                    }
                }
                _ => needs_flush = true,
            }
        }
        if self.prev_range != 0 && self.prev_range != dp_range {
            needs_flush = true;
        }
        match rec.as_deref() {
            None => needs_flush = true,
            Some(r) => {
                if self.rid != r.rid() || r.pos() > self.end + 1 {
                    needs_flush = true;
                }
            }
        }
        if self.prev_range != 0 && needs_flush {
            if let Some(r) = rec.as_deref()
                && r.rid() == self.rid
                && r.pos() == self.end
            {
                self.end -= 1;
            }
            self.end += 1;
            let line = &mut self.line;
            line.clear();
            line.set_rid(self.rid);
            line.set_pos(self.start);
            line.update_alleles_str(hdr, &self.als);
            line.set_rlen(self.end - self.start);
            if self.start + 1 < self.end {
                line.update_info_i32(hdr, b"END", &[self.end as i32])?;
            }
            line.update_info_i32(hdr, b"MIN_DP", &[self.min_dp])?;
            if !self.qsum.is_empty() {
                line.update_info_f32(hdr, b"QS", &self.qsum)?;
            }
            if !self.gts.is_empty() {
                line.update_genotypes(hdr, &self.gts)?;
            }
            if self.npl > 0 {
                line.update_format_i32(hdr, b"PL", &self.pl[..self.npl as usize])?;
            }
            line.update_format_i32(hdr, b"DP", &self.dp[..nsmpl])?;
            emit(line)?;
            self.prev_range = 0;
            self.rid = -1;
            self.npl = 0;
            self.qsum.clear();
            self.gts.clear();
            if rec.is_none() {
                return Ok(false);
            }
        }
        let rec = rec.expect("only flushing returns above");
        if can_collapse {
            if self.prev_range == 0 {
                self.dp.clear();
                self.dp.extend_from_slice(&self.tmp[..nsmpl]);
                self.npl = match rec.get_format_i32(hdr, b"PL", &mut self.pl) {
                    Ok(n) => n as i32,
                    Err(_) => -1,
                };
                let mut q = Vec::new();
                self.qsum = match rec.get_info_f32(hdr, b"QS", &mut q) {
                    Ok(_) => q,
                    Err(_) => Vec::new(),
                };
                let mut g = Vec::new();
                self.gts = match rec.get_format_i32(hdr, b"GT", &mut g) {
                    Ok(_) => g,
                    Err(_) => Vec::new(),
                };
                self.rid = rec.rid();
                self.start = rec.pos();
                self.als.clear();
                for i in 0..rec.allele_count() {
                    if i > 0 {
                        self.als.push(b',');
                    }
                    self.als.extend_from_slice(rec.allele(i));
                }
                self.min_dp = min_dp;
            } else {
                if self.min_dp > min_dp {
                    self.min_dp = min_dp;
                }
                for i in 0..nsmpl {
                    if self.dp[i] > self.tmp[i] {
                        self.dp[i] = self.tmp[i];
                    }
                }
                match rec.get_format_i32(hdr, b"PL", &mut self.tmp) {
                    Ok(n) => {
                        if n != nsmpl * 3 {
                            bail!("Unexpected number of PL fields");
                        }
                        for i in 0..nsmpl {
                            if self.pl[3 * i + 1] > self.tmp[3 * i + 1] {
                                self.pl[3 * i + 1] = self.tmp[3 * i + 1];
                                self.pl[3 * i + 2] = self.tmp[3 * i + 2];
                            } else if self.pl[3 * i + 1] == self.tmp[3 * i + 1]
                                && self.pl[3 * i + 2] > self.tmp[3 * i + 2]
                            {
                                self.pl[3 * i + 2] = self.tmp[3 * i + 2];
                            }
                        }
                    }
                    Err(_) => self.npl = 0,
                }
            }
            self.prev_range = dp_range;
            let mut e = Vec::new();
            self.end = match rec.get_info_i32(hdr, b"END", &mut e) {
                Ok(1) => i64::from(e[0]) - 1,
                _ => rec.pos(),
            };
            return Ok(false);
        }
        if is_ref && min_dp != 0 {
            rec.update_info_i32(hdr, b"MIN_DP", &[min_dp])?;
        }
        Ok(true)
    }
}
