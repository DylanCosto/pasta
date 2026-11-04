//! The report: the INFO counters, the DCv2 lines of every pair (or the best matches of each
//! sample), and the DS lines of `--distinctive-sites` (bcftools' `report` and
//! `report_distinctive_sites`).

use std::{
    cmp::{Ordering, Reverse},
    collections::BinaryHeap,
    fmt::Write as _,
    fs::{File, OpenOptions},
    io::{self, BufWriter, Seek, Write},
    mem,
    ops::Range,
    os::unix::fs::{FileExt, OpenOptionsExt},
};

use anyhow::Context;

use super::{Gtcheck, kernels::Discordance};

/// C's `%e`: one digit, six decimals, and an exponent of at least two digits.
fn push_e(out: &mut Vec<u8>, tmp: &mut String, x: f64) {
    if x.is_nan() {
        out.extend_from_slice(if x.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
        return;
    }
    if x.is_infinite() {
        out.extend_from_slice(if x < 0.0 { b"-inf" } else { b"inf" });
        return;
    }
    tmp.clear();
    write!(tmp, "{x:.6e}").expect("writing to a String cannot fail");
    let (mantissa, exp) = tmp.split_once('e').expect("{:e} has an exponent");
    out.extend_from_slice(mantissa.as_bytes());
    out.push(b'e');
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => (b'-', d),
        None => (b'+', exp),
    };
    out.push(sign);
    if digits.len() < 2 {
        out.push(b'0');
    }
    out.extend_from_slice(digits.as_bytes());
}

/// Appends formatted text to a line of the report.
fn push_fmt(out: &mut Vec<u8>, args: std::fmt::Arguments<'_>) {
    out.write_fmt(args).expect("writing to a Vec cannot fail");
}

const DCV2_HEADER: &str = "# DCv2, discordance version 2:
#     - Query sample
#     - Genotyped sample
#     - Discordance, given either as an abstract score or number of mismatches, see the options -E/-u
#       in man page for details. Note that samples with high missingness have fewer sites compared,
#       which results in lower overall discordance. Therefore it is advisable to use the average score
#       per site rather than the absolute value, i.e. divide the value by the number of sites compared
#       (smaller value = better match)
#     - Average negative log of HWE probability at matching sites, attempts to quantify the following
#       intuition: rare genotype matches are more informative than common genotype matches, hence two
#       samples with similar discordance can be further stratified by the HWE score (bigger value = better
#       match, the observed concordance was less likely to occur by chance)
#     - Number of sites compared for this pair of samples (bigger = more informative)
#     - Number of matching genotypes
#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes
";

/// The scratch space of the DCv2 lines, each built here before it is written.
#[derive(Debug, Default)]
struct LineBuf {
    line: Vec<u8>,
    /// A number formatted for [`push_e`].
    num: String,
}

/// One candidate match of a sample when only the best `--n-matches` are printed (bcftools'
/// `idbl_t`).
#[derive(Clone, Copy, Debug)]
struct Candidate {
    /// The other sample, as an index into the sample list.
    ism: usize,
    /// The pair.
    idx: usize,
    val: f64,
}

/// bcftools' `cmp_idbl`: by value, then by sample.
fn cmp_candidate(a: &Candidate, b: &Candidate) -> Ordering {
    if a.val < b.val {
        Ordering::Less
    } else if a.val > b.val {
        Ordering::Greater
    } else {
        a.ism.cmp(&b.ism)
    }
}

impl Gtcheck {
    /// Writes the counters and the DCv2 section (bcftools' `report`).
    pub(super) fn report(&self, out: &mut impl Write) -> io::Result<()> {
        let c = &self.counts;
        let info = [
            ("sites-compared", c.ncmp),
            ("sites-skipped-no-match", c.no_match),
            ("sites-skipped-multiallelic", c.not_biallelic),
            ("sites-skipped-monoallelic", c.mono),
            ("sites-skipped-no-data", c.no_data),
            ("sites-skipped-GT-not-diploid", c.dip_gt),
            ("sites-skipped-PL-not-diploid", c.dip_pl),
            ("sites-skipped-filtering-expression", c.filter),
            ("sites-used-PL-vs-PL", c.nused[0][0]),
            ("sites-used-PL-vs-GT", c.nused[0][1]),
            ("sites-used-GT-vs-PL", c.nused[1][0]),
            ("sites-used-GT-vs-GT", c.nused[1][1]),
        ];
        for (name, n) in info {
            writeln!(out, "INFO\t{name}\t{n}")?;
        }
        out.write_all(DCV2_HEADER.as_bytes())?;

        let mut buf = LineBuf::default();
        if let Some(pairs) = &self.pairs {
            for (i, p) in pairs.iter().enumerate() {
                self.dcv2_line(out, &mut buf, p.qry, p.gt, i)?;
            }
            return Ok(());
        }
        let ntop = self.ntop;
        let trim = ntop != 0 && self.ngt > ntop;
        if !trim {
            let mut idx = 0;
            for i in 0..self.nqry {
                let iqry = self.qry_sample(i);
                let ngt = if self.cross_check { i } else { self.ngt };
                for j in 0..ngt {
                    let k = if self.cross_check {
                        idx
                    } else {
                        self.rect_index(i, j)
                    };
                    self.dcv2_line(out, &mut buf, iqry, self.gt_sample(j), k)?;
                    idx += 1;
                }
            }
        } else if !self.cross_check {
            let mut arr = Vec::with_capacity(self.ngt);
            for i in 0..self.nqry {
                arr.clear();
                arr.extend((0..self.ngt).map(|j| {
                    let idx = self.rect_index(i, j);
                    Candidate {
                        ism: j,
                        idx,
                        val: self.sort_value(idx),
                    }
                }));
                arr.sort_unstable_by(cmp_candidate);
                let iqry = self.qry_sample(i);
                for cand in &arr[..ntop] {
                    self.dcv2_line(out, &mut buf, iqry, self.gt_sample(cand.ism), cand.idx)?;
                }
            }
        } else {
            let narr = self.nqry - 1;
            let mut arr = Vec::with_capacity(narr);
            for i in 0..self.nqry {
                arr.clear();
                let start = i * i.saturating_sub(1) / 2;
                arr.extend((0..i).map(|j| Candidate {
                    ism: j,
                    idx: start + j,
                    val: self.sort_value(start + j),
                }));
                arr.extend((i..narr).map(|j| {
                    let idx = j * (j + 1) / 2 + i;
                    Candidate {
                        ism: j + 1,
                        idx,
                        val: self.sort_value(idx),
                    }
                }));
                arr.sort_unstable_by(cmp_candidate);
                let iqry = self.qry_sample(i);
                // Each pair is printed once, by the later of its two samples.
                for cand in arr[..ntop].iter().filter(|c| c.ism < i) {
                    self.dcv2_line(out, &mut buf, iqry, self.qry_sample(cand.ism), cand.idx)?;
                }
            }
        }
        Ok(())
    }

    /// The value that ranks pair `idx` for `--n-matches`: the average discordance, or the
    /// negated average HWE score.
    fn sort_value(&self, idx: usize) -> f64 {
        let s = &self.sums;
        if self.sort_by_hwe {
            // bcftools dereferences the missing nmatch array with --no-HWE-prob.
            let n = s.nmatch(idx);
            return match &s.hwe {
                Some((h, _)) if n != 0 => -h[idx] / f64::from(n),
                _ => 0.0,
            };
        }
        let ncnt = s.ncnt[idx];
        if ncnt == 0 {
            return 0.0;
        }
        match &s.disc {
            Discordance::Count(d) => f64::from(d[idx]) / f64::from(ncnt),
            Discordance::Score(d) => d[idx] / f64::from(ncnt),
        }
    }

    /// Writes the DCv2 line of pair `idx` between header samples `iqry` and `igt`.
    fn dcv2_line(
        &self,
        out: &mut impl Write,
        buf: &mut LineBuf,
        iqry: usize,
        igt: usize,
        idx: usize,
    ) -> io::Result<()> {
        let s = &self.sums;
        let LineBuf { line, num } = buf;
        line.clear();
        line.extend_from_slice(b"DCv2\t");
        line.extend_from_slice(&self.qry_names[iqry]);
        line.push(b'\t');
        line.extend_from_slice(&self.gt_names()[igt]);
        line.push(b'\t');
        match &s.disc {
            Discordance::Count(d) => push_fmt(line, format_args!("{}", d[idx])),
            Discordance::Score(d) => push_e(line, num, d[idx]),
        }
        line.push(b'\t');
        push_e(line, num, s.hwe_avg(idx));
        push_fmt(line, format_args!("\t{}\t{}\n", s.ncnt[idx], s.nmatch(idx)));
        out.write_all(line)
    }

    /// Writes the DS section of `--distinctive-sites` (bcftools'
    /// `report_distinctive_sites`).
    pub(super) fn report_distinctive_sites(&mut self, out: &mut impl Write) -> anyhow::Result<()> {
        let Some(ds) = self.distinctive.as_mut() else {
            return Ok(());
        };
        out.write_all(
            b"# DS, distinctive sites:\n\
              #     - chromosome\n\
              #     - position\n\
              #     - cumulative number of pairs distinguished by this block\n\
              #     - block id\n\
              #DS\t[2]Chromosome\t[3]Position\t[4]Cumulative number of distinct pairs\t[5]Block id\n",
        )?;
        let ndiff_min = ds.ndiff_min(self.npairs);
        let mut block = vec![0u64; ds.words];
        let (mut ndiff_tot, mut iblock) = (0i64, 0i64);
        let header = self.sr.header(0);
        ds.for_each_sorted(|site, bits| {
            let mut ndiff_new = 0;
            for (b, &w) in block.iter_mut().zip(bits) {
                ndiff_new += (w & !*b).count_ones();
                *b |= w;
            }
            if ndiff_new == 0 {
                return Ok(());
            }
            ndiff_tot += i64::from(ndiff_new);
            let chrom = header.contig_name(site.rid).unwrap_or(b"(null)");
            out.write_all(b"DS\t")?;
            out.write_all(chrom)?;
            writeln!(out, "\t{}\t{ndiff_tot}\t{iblock}", i64::from(site.pos) + 1)?;
            if ndiff_tot >= ndiff_min {
                iblock += 1;
                ndiff_tot = 0;
                block.fill(0);
            }
            Ok(())
        })
    }

    /// The names of the header the `-g` samples come from.
    fn gt_names(&self) -> &[Box<[u8]>] {
        if self.two_files {
            self.sr.header(1).samples()
        } else {
            &self.qry_names
        }
    }
}

/// The sites of `--distinctive-sites`: those where some pairs differ, with the set of pairs.
/// They are held in memory up to MEM, counted as bcftools counts them (a pointer and a record
/// per site); beyond it they go in sorted runs to a temporary file, merged for the report, as
/// bcftools' `extsort` does. Unlike `extsort`, which merges all the runs at once, at most
/// [`fan_in`] runs are merged at a time, in passes, so the merge's memory does not grow with
/// the number of runs.
#[derive(Debug)]
pub(super) struct Distinctive {
    /// The requested number of pairs (bcftools' updated `distinctive_sites`).
    target: f64,
    words: usize,
    sites: Vec<DiffSite>,
    /// The differing pairs of each site, `words` words per site.
    bits: Vec<u64>,
    /// The differing pairs of the current site.
    pub(super) current: Vec<u64>,
    rng: Rand48,
    /// MEM, in bytes.
    max_mem: usize,
    /// The most sites held in memory, from MEM: each run of the temporary file has this many
    /// sites but the last.
    run_cap: usize,
    /// TMP: the temporary file is named from it and `XXXXXX`.
    tmp_prefix: String,
    spilled: Option<Spilled>,
}

#[derive(Clone, Copy, Debug)]
struct DiffSite {
    ndiff: u32,
    rid: i32,
    pos: i32,
    /// Shuffles sites with the same count (bcftools' `rand`).
    rand: u32,
}

impl DiffSite {
    /// The report's order: by decreasing number of pairs, ties in random order and then in
    /// the order the sites came (`diff_sites_cmp`, with a stable sort).
    fn key(&self, index: u64) -> (Reverse<u32>, u32, u64) {
        (Reverse(self.ndiff), self.rand, index)
    }
}

/// The sites written to the temporary file, in sorted runs, each record a [`DiffSite`], its
/// index among all the sites (for ties), and its pairs.
#[derive(Debug)]
struct Spilled {
    file: File,
    /// The number of records of each run; the last run may have fewer.
    run_len: u64,
    /// The number of records written.
    len: u64,
}

impl Spilled {
    fn runs(&self) -> u64 {
        self.len.div_ceil(self.run_len)
    }

    /// Merges groups of `fan_in` runs into longer runs, in a second temporary file named from
    /// `tmp_prefix` and then back and forth, until there are at most `fan_in` runs. The records
    /// keep their indices, so the order of ties does not change.
    ///
    /// # Errors
    ///
    /// Fails when a temporary file cannot be made, read, or written.
    fn merge_down(
        &mut self,
        fan_in: usize,
        record: usize,
        budget: usize,
        tmp_prefix: &str,
    ) -> anyhow::Result<()> {
        let fan_in = fan_in as u64;
        if self.runs() <= fan_in {
            return Ok(());
        }
        let mut other = temp_file(tmp_prefix)?;
        while self.runs() > fan_in {
            other.set_len(0).context(WRITE_ERR)?;
            other.rewind().context(WRITE_ERR)?;
            let group = self.run_len * fan_in;
            let mut w = BufWriter::new(&other);
            let mut first = 0;
            while first < self.len {
                let end = self.len.min(first + group);
                merge_runs(
                    &self.file,
                    first..end,
                    self.run_len,
                    record,
                    budget,
                    |rec| w.write_all(rec).context(WRITE_ERR),
                )?;
                first = end;
            }
            w.flush().context(WRITE_ERR)?;
            drop(w);
            mem::swap(&mut self.file, &mut other);
            self.run_len = group;
        }
        Ok(())
    }
}

/// The bytes of a [`DiffSite`] and its index in the temporary file, before its pairs.
const RECORD_HEAD: usize = 24;

/// The most bytes read from one run of the temporary file at a time.
const RUN_BUFFER: usize = 1 << 20;

/// The least memory the read buffers of a merge share, whatever MEM is: with a MEM of a few
/// records, merging only as many runs at a time would take many more passes over the file.
const MERGE_MEM_MIN: usize = 1 << 20;

/// The fewest runs merged at a time, even when MEM holds fewer records.
const MIN_FAN_IN: usize = 2;

/// The most runs merged at a time, however many records MEM holds.
const MAX_FAN_IN: usize = 1024;

const READ_ERR: &str = "Error: failed to read from the temporary file";
const WRITE_ERR: &str = "Error: failed to write to the temporary file";

/// The memory the read buffers of a merge share: MEM, or [`MERGE_MEM_MIN`] when MEM is less.
fn merge_budget(max_mem: usize) -> usize {
    max_mem.max(MERGE_MEM_MIN)
}

/// The most runs merged at a time: as many records as the merge's memory holds, from
/// [`MIN_FAN_IN`] to [`MAX_FAN_IN`].
fn fan_in(max_mem: usize, record: usize) -> usize {
    (merge_budget(max_mem) / record).clamp(MIN_FAN_IN, MAX_FAN_IN)
}

/// Merges the sorted runs of `run_len` records that make up records `range` of `file`,
/// calling `emit` with each record in the report's order. The runs' read buffers share
/// `budget` bytes, with one record each at least.
fn merge_runs(
    file: &File,
    range: Range<u64>,
    run_len: u64,
    record: usize,
    budget: usize,
    mut emit: impl FnMut(&[u8]) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let nruns = (range.end - range.start).div_ceil(run_len);
    let share = budget / usize::try_from(nruns).unwrap_or(usize::MAX).max(1);
    let per_read = (share.min(RUN_BUFFER) / record).max(1) as u64;
    let mut runs: Vec<RunReader> = (0..nruns)
        .map(|i| {
            let next = range.start + i * run_len;
            RunReader {
                next,
                end: range.end.min(next + run_len),
                buf: Vec::new(),
                at: 0,
            }
        })
        .collect();
    let mut heap = BinaryHeap::with_capacity(runs.len());
    for (r, run) in runs.iter_mut().enumerate() {
        if run.fill(file, record, per_read).context(READ_ERR)? {
            heap.push(Reverse((run.key(), r)));
        }
    }
    while let Some(Reverse((_, r))) = heap.pop() {
        let run = &mut runs[r];
        emit(&run.buf[run.at..run.at + record])?;
        run.at += record;
        if run.fill(file, record, per_read).context(READ_ERR)? {
            heap.push(Reverse((run.key(), r)));
        }
    }
    Ok(())
}

impl Distinctive {
    /// bcftools' `diff_sites_init`: resolves the number of pairs (a fraction of them when at
    /// most 1).
    ///
    /// # Errors
    ///
    /// Fails when that number is not positive.
    pub(super) fn new(value: f64, npairs: usize) -> anyhow::Result<Self> {
        let n = if value <= 1.0 {
            npairs as f64 * value
        } else {
            value
        };
        let nsites = c_double_to_int(n);
        if nsites <= 0 {
            anyhow::bail!("The value for --distinctive-sites was set too low: {nsites}");
        }
        let target = if nsites as usize > npairs {
            eprintln!(
                "Warning: The value for --distinctive-sites is bigger than is the number of pairs, all discordant sites be printed."
            );
            npairs as f64 + 1.0
        } else {
            f64::from(nsites)
        };
        let words = npairs.div_ceil(64);
        Ok(Self {
            target,
            words,
            sites: Vec::new(),
            bits: Vec::new(),
            current: vec![0; words],
            rng: Rand48::new(0),
            max_mem: usize::MAX,
            run_cap: usize::MAX,
            tmp_prefix: String::new(),
            spilled: None,
        })
    }

    /// Holds at most `max_mem` bytes of sites, writing the others to a temporary file named
    /// from `tmp_prefix`.
    pub(super) fn spill_to(&mut self, max_mem: usize, tmp_prefix: &str) {
        // bcftools' `diff_sites_t` (three ints, a long, and the words) and a pointer to it.
        let site_mem = 32 + 8 * self.words;
        self.max_mem = max_mem;
        self.run_cap = (max_mem / site_mem).max(1);
        tmp_prefix.clone_into(&mut self.tmp_prefix);
    }

    /// Clears the current site's pairs (`diff_sites_reset`).
    pub(super) fn reset(&mut self) {
        self.current.fill(0);
    }

    /// Marks pair `i` as differing at the current site (`kbs_insert`).
    #[inline]
    pub(super) fn insert(&mut self, i: usize) {
        self.current[i / 64] |= 1 << (i % 64);
    }

    /// Keeps the current site (`diff_sites_push` and `extsort_push`).
    ///
    /// # Errors
    ///
    /// Fails when the sites held cannot be written to the temporary file.
    pub(super) fn push(&mut self, ndiff: u32, rid: i32, pos: i64) -> anyhow::Result<()> {
        if self.sites.len() >= self.run_cap {
            self.spill()?;
        }
        self.sites.push(DiffSite {
            ndiff,
            rid,
            pos: pos as i32,
            rand: self.rng.next(),
        });
        self.bits.extend_from_slice(&self.current);
        Ok(())
    }

    /// The indices of the sites held, in the report's order.
    fn sorted(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.sites.len()).collect();
        order.sort_by_key(|&k| self.sites[k].key(0));
        order
    }

    /// Writes the sites held to the temporary file as a sorted run (`_buf_flush`).
    fn spill(&mut self) -> anyhow::Result<()> {
        let order = self.sorted();
        let spilled = match &mut self.spilled {
            Some(s) => s,
            None => self.spilled.insert(Spilled {
                file: temp_file(&self.tmp_prefix)?,
                run_len: self.run_cap as u64,
                len: 0,
            }),
        };
        let first = spilled.len;
        let mut w = BufWriter::new(&spilled.file);
        for k in order {
            let site = &self.sites[k];
            for v in [site.ndiff, site.rid as u32, site.pos as u32, site.rand] {
                w.write_all(&v.to_ne_bytes())?;
            }
            w.write_all(&(first + k as u64).to_ne_bytes())?;
            for word in &self.bits[k * self.words..(k + 1) * self.words] {
                w.write_all(&word.to_ne_bytes())?;
            }
        }
        w.flush()
            .context("Error: failed to write to the temporary file")?;
        spilled.len += self.sites.len() as u64;
        self.sites.clear();
        self.bits.clear();
        Ok(())
    }

    /// Calls `f` with each site and its pairs in the report's order, merging the runs of the
    /// temporary file when there are any (`extsort_sort` and `extsort_shift`).
    fn for_each_sorted(
        &mut self,
        f: impl FnMut(&DiffSite, &[u64]) -> io::Result<()>,
    ) -> anyhow::Result<()> {
        let fan_in = fan_in(self.max_mem, RECORD_HEAD + 8 * self.words);
        self.merge_sorted(fan_in, f)
    }

    /// [`Self::for_each_sorted`], merging at most `fan_in` runs at a time.
    fn merge_sorted(
        &mut self,
        fan_in: usize,
        mut f: impl FnMut(&DiffSite, &[u64]) -> io::Result<()>,
    ) -> anyhow::Result<()> {
        let words = self.words;
        if self.spilled.is_none() {
            for k in self.sorted() {
                f(&self.sites[k], &self.bits[k * words..(k + 1) * words])?;
            }
            return Ok(());
        }
        if !self.sites.is_empty() {
            self.spill()?;
        }
        // The merge's buffers take the place of the sites'.
        self.sites = Vec::new();
        self.bits = Vec::new();
        let Some(spilled) = &mut self.spilled else {
            unreachable!("the sites were written out");
        };
        let record = RECORD_HEAD + 8 * words;
        let budget = merge_budget(self.max_mem);
        spilled.merge_down(fan_in, record, budget, &self.tmp_prefix)?;
        let mut bits = vec![0u64; words];
        merge_runs(
            &spilled.file,
            0..spilled.len,
            spilled.run_len,
            record,
            budget,
            |rec| {
                let (site, _) = record_site(rec);
                record_bits(rec, &mut bits);
                Ok(f(&site, &bits)?)
            },
        )
    }

    fn ndiff_min(&self, npairs: usize) -> i64 {
        if self.target <= npairs as f64 {
            i64::from(c_double_to_int(self.target))
        } else {
            npairs as i64
        }
    }
}

/// A run of the temporary file being merged, read a buffer at a time.
struct RunReader {
    /// The next record to read into the buffer, and the end of the run.
    next: u64,
    end: u64,
    buf: Vec<u8>,
    /// Where the current record starts in the buffer.
    at: usize,
}

impl RunReader {
    /// Reads the next records when the buffer is used up; returns whether there is a current
    /// record.
    fn fill(&mut self, file: &File, record: usize, per_read: u64) -> io::Result<bool> {
        if self.at < self.buf.len() {
            return Ok(true);
        }
        let n = per_read.min(self.end - self.next);
        if n == 0 {
            return Ok(false);
        }
        self.buf.resize(n as usize * record, 0);
        file.read_exact_at(&mut self.buf, self.next * record as u64)?;
        self.next += n;
        self.at = 0;
        Ok(true)
    }

    /// The current record's place in the report's order.
    fn key(&self) -> (Reverse<u32>, u32, u64) {
        let (site, index) = record_site(&self.buf[self.at..]);
        site.key(index)
    }
}

/// The site and index of a record of the temporary file.
fn record_site(rec: &[u8]) -> (DiffSite, u64) {
    let u32_at = |i: usize| u32::from_ne_bytes(rec[i..i + 4].try_into().expect("4 bytes"));
    let site = DiffSite {
        ndiff: u32_at(0),
        rid: u32_at(4) as i32,
        pos: u32_at(8) as i32,
        rand: u32_at(12),
    };
    let index = u64::from_ne_bytes(rec[16..RECORD_HEAD].try_into().expect("8 bytes"));
    (site, index)
}

/// The pairs of a record of the temporary file.
fn record_bits(rec: &[u8], bits: &mut [u64]) {
    for (w, chunk) in bits.iter_mut().zip(rec[RECORD_HEAD..].as_chunks::<8>().0) {
        *w = u64::from_ne_bytes(*chunk);
    }
}

/// Creates a temporary file named from `prefix` and six more characters, as `mkstemp`
/// does, and removes its name: the file goes when it is closed.
fn temp_file(prefix: &str) -> anyhow::Result<File> {
    const CHARS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| u64::from(d.subsec_nanos()));
    let mut seed = (u64::from(std::process::id()) << 32) ^ nanos;
    for _ in 0..1000 {
        let mut name = prefix.to_owned();
        for _ in 0..6 {
            seed = seed.wrapping_mul(0x5_deec_e66d).wrapping_add(0xb);
            name.push(char::from(CHARS[(seed >> 33) as usize % CHARS.len()]));
        }
        match OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&name)
        {
            Ok(file) => {
                let _ = std::fs::remove_file(&name);
                return Ok(file);
            }
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
            Err(_) => break,
        }
    }
    anyhow::bail!("Error: failed to open a temporary file {prefix}XXXXXX")
}

/// C's conversion of a double to `int` as x86-64 does it: truncation, with `INT_MIN` for
/// values out of range and NaN.
fn c_double_to_int(x: f64) -> i32 {
    if x.is_nan() || x >= 2_147_483_648.0 || x <= -2_147_483_649.0 {
        i32::MIN
    } else {
        x as i32
    }
}

/// glibc's `srand48`/`lrand48` generator, which HTSlib's `hts_lrand48` uses on Linux.
#[derive(Clone, Copy, Debug)]
struct Rand48(u64);

impl Rand48 {
    const MASK: u64 = (1 << 48) - 1;

    fn new(seed: i64) -> Self {
        Self((((seed as u64) & 0xffff_ffff) << 16) | 0x330e)
    }

    fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(0x5_deec_e66d).wrapping_add(0xb) & Self::MASK;
        (self.0 >> 17) as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c_style_exponents() {
        let mut out = Vec::new();
        let mut tmp = String::new();
        for x in [0.0, 1.5, 123_456.789, 1e-7, 2.5e100, f64::INFINITY] {
            push_e(&mut out, &mut tmp, x);
            out.push(b' ');
        }
        assert_eq!(
            String::from_utf8(out).expect("ascii"),
            "0.000000e+00 1.500000e+00 1.234568e+05 1.000000e-07 2.500000e+100 inf "
        );
    }

    /// The sites `for_each_sorted` gives, as (position, ndiff, pairs).
    fn report_order(ds: &mut Distinctive) -> Vec<(i32, u32, Vec<u64>)> {
        let mut got = Vec::new();
        ds.for_each_sorted(|site, bits| {
            got.push((site.pos, site.ndiff, bits.to_vec()));
            Ok(())
        })
        .expect("the sites are read back");
        got
    }

    /// A directory for one test's temporary files, removed when dropped.
    struct TestDir(std::path::PathBuf);

    impl TestDir {
        fn new(name: &str) -> Self {
            let dir = std::env::temp_dir().join(format!("vtools-ds-{name}-{}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("create the test directory");
            Self(dir)
        }

        fn prefix(&self) -> String {
            format!("{}/ds", self.0.display())
        }

        fn files(&self) -> usize {
            std::fs::read_dir(&self.0)
                .expect("read the test directory")
                .count()
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// 500 sites of 130 pairs, holding at most `max_mem` bytes of them in memory.
    fn sites(max_mem: usize, prefix: &str) -> Distinctive {
        let npairs = 130;
        let mut ds = Distinctive::new(0.5, npairs).expect("a valid count");
        ds.spill_to(max_mem, prefix);
        let mut x = 1u64;
        for pos in 0..500 {
            ds.reset();
            let mut nd = 0;
            for i in 0..npairs {
                x = x.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                // Few pairs per site, so that many sites tie on the count.
                if x >> 60 == 0 {
                    ds.insert(i);
                    nd += 1;
                }
            }
            ds.push(nd, 0, pos).expect("the sites are written out");
        }
        ds
    }

    #[test]
    fn sites_past_mem_are_merged_back_in_order() {
        let dir = TestDir::new("spill");
        let prefix = dir.prefix();
        let mut runs = Vec::new();
        for max_mem in [usize::MAX, 1, 1000, 10_000] {
            let mut ds = sites(max_mem, &prefix);
            assert_eq!(ds.spilled.is_some(), max_mem != usize::MAX, "MEM {max_mem}");
            runs.push(report_order(&mut ds));
        }
        assert_eq!(dir.files(), 0, "the temporary file is removed at once");
        for (got, mem) in runs[1..].iter().zip([1, 1000, 10_000]) {
            assert!(got == &runs[0], "MEM {mem} gives the in-memory order");
        }
        assert!(runs[0].windows(2).all(|w| w[0].1 >= w[1].1));
    }

    #[test]
    fn many_runs_are_merged_a_few_at_a_time() {
        let dir = TestDir::new("passes");
        let prefix = dir.prefix();
        let want = report_order(&mut sites(usize::MAX, &prefix));
        // MEM 1 and 120 bytes: runs of 1 and 2 sites (bcftools counts 56 bytes per site).
        for (max_mem, nruns) in [(1, 500), (120, 250)] {
            for fan_in in [2, 3, 7, 1000] {
                let mut ds = sites(max_mem, &prefix);
                let spilled = ds.spilled.as_ref().expect("the sites are written out");
                // The last run is still in memory.
                assert_eq!(spilled.runs() + 1, nruns, "MEM {max_mem}");
                let mut got = Vec::new();
                ds.merge_sorted(fan_in, |site, bits| {
                    got.push((site.pos, site.ndiff, bits.to_vec()));
                    Ok(())
                })
                .expect("the sites are read back");
                let spilled = ds.spilled.as_ref().expect("the sites are written out");
                assert!(
                    spilled.runs() <= fan_in as u64,
                    "MEM {max_mem}: {} runs merged at once, at most {fan_in} wanted",
                    spilled.runs()
                );
                assert!(got == want, "MEM {max_mem}, {fan_in} runs at a time");
                assert_eq!(dir.files(), 0, "the temporary files are removed at once");
            }
        }
    }

    #[test]
    fn fan_in_follows_mem_within_bounds() {
        // MEM below the merge's minimum memory counts as that minimum.
        assert_eq!(fan_in(4096, 2512), MERGE_MEM_MIN / 2512);
        assert_eq!(fan_in(10_000_000, 100_000), 100);
        // Records too big for MEM still merge two runs at a time.
        assert_eq!(fan_in(4096, 4 << 20), MIN_FAN_IN);
        assert_eq!(fan_in(500_000_000, 32), MAX_FAN_IN);
    }

    #[test]
    fn rand48_matches_glibc() {
        // glibc: srand48(0); lrand48() x3.
        let mut r = Rand48::new(0);
        assert_eq!(
            [r.next(), r.next(), r.next()],
            [366_850_414, 1_610_402_240, 206_956_554]
        );
    }
}
