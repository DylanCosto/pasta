//! `concat -l`: ligate overlapping phased chunks, matching the phase of the overlap and
//! breaking phase sets where the match is poor (bcftools' `phased_push`, `phased_flush`, and
//! the phased part of `concat`).

use anyhow::bail;
use pasta::variant::{
    Genotype, Header, Record, Unpack, encoding::INT32_MISSING, io::Writer, synced::SyncedReaders,
};

use super::{Concat, ExitOne, LastWrite, active_index, next_line};
use crate::output;

/// Whether a site lies where two chunks overlap (bcftools' `is_overlap`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Site {
    /// Only one chunk covers the site: its record is written at once.
    Single,
    /// Both chunks cover the site: the pair is buffered until the overlap ends.
    Overlap,
}

/// The chunk of an overlapping pair a buffered record comes from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Chunk {
    First,
    Second,
}

/// The ligation options.
pub struct Options {
    /// `-c`: write PS only where a phase set starts.
    pub compact_ps: bool,
    /// `-q`: break the phase set below this phasing quality.
    pub min_pq: i32,
    /// `--ligate-force`.
    pub force: bool,
    /// `--ligate-warn`.
    pub warn: bool,
}

/// The ligation state (the phasing fields of bcftools' `args_t`).
pub struct Ligate {
    opts: Options,
    start_pos: Vec<i64>,
    ifname: usize,
    /// Record pairs of the overlap: `buf[2k]` from the first file, `buf[2k+1]` from the
    /// second, present per `buf_mask[k]` (1, 2, or 3).
    buf: Vec<Record>,
    buf_mask: Vec<u8>,
    nbuf: usize,
    swap_phase: Vec<bool>,
    nswap: usize,
    nmatch: Vec<i32>,
    nmism: Vec<i32>,
    phase_qual: Vec<i32>,
    phase_set: Vec<i32>,
    phase_set_changed: bool,
    prev_chr: i64,
    prev_pos_check: i64,
    gta: Vec<i32>,
    gtb: Vec<i32>,
    site_drop_warned: bool,
    gt_absent_warned: bool,
}

/// The site of a record, as `CHROM:POS` for messages.
fn site(hdr: &Header, rec: &Record) -> String {
    format!(
        "{}:{}",
        String::from_utf8_lossy(rec.chrom_safe(hdr)),
        rec.pos() + 1
    )
}

impl Ligate {
    /// Removes empty files and checks the start positions (the phased part of bcftools'
    /// `init_data`).
    pub fn new(state: &mut Concat, start_pos: Vec<i64>, opts: Options) -> anyhow::Result<Self> {
        let mut kept_names = Vec::new();
        let mut kept_pos = Vec::new();
        for (name, &pos) in state.fnames.iter().zip(&start_pos) {
            if pos != -2 {
                kept_names.push(name.clone());
                kept_pos.push(pos);
            }
        }
        state.fnames = kept_names;
        let start_pos = kept_pos;
        for i in 1..start_pos.len() {
            if start_pos[i - 1] != -1 && start_pos[i] != -1 && start_pos[i] < start_pos[i - 1] {
                bail!(
                    "The files not in ascending order: {} in {}, {} in {}",
                    start_pos[i - 1] + 1,
                    state.fnames[i - 1],
                    start_pos[i] + 1,
                    state.fnames[i]
                );
            }
        }
        let nsmpl = state.out_hdr.sample_count();
        Ok(Self {
            opts,
            start_pos,
            ifname: 0,
            buf: Vec::new(),
            buf_mask: Vec::new(),
            nbuf: 0,
            swap_phase: vec![false; nsmpl],
            nswap: 0,
            nmatch: vec![0; nsmpl],
            nmism: vec![0; nsmpl],
            phase_qual: vec![0; nsmpl],
            phase_set: vec![0; nsmpl],
            phase_set_changed: false,
            prev_chr: -1,
            prev_pos_check: -1,
            gta: Vec::new(),
            gtb: Vec::new(),
            site_drop_warned: false,
            gt_absent_warned: false,
        })
    }

    fn warn_drop(&mut self, chr: &str, pos: i64) {
        if self.site_drop_warned {
            return;
        }
        eprintln!(
            "Warning: Dropping the site {chr}:{}. The --ligate option is intended for VCFs with perfect\n         overlap, sites in overlapping regions present in one but missing in other are dropped.\n         This warning is printed only once.",
            pos + 1
        );
        self.site_drop_warned = true;
    }

    /// The phased part of bcftools' `concat`: keeps at most two files open and feeds their
    /// records to [`Self::push`].
    pub fn run(
        &mut self,
        state: &mut Concat,
        files: &mut SyncedReaders,
        writer: &mut Writer,
    ) -> anyhow::Result<()> {
        let nfnames = state.fnames.len();
        while self.ifname < nfnames {
            let mut new_file = false;
            while files.reader_count() < 2 && self.ifname < nfnames {
                let f = &state.fnames[self.ifname];
                if let Err(e) = files.add_reader_again(f) {
                    bail!("Failed to open {f}: {e}");
                }
                new_file = true;
                self.ifname += 1;
                if self.start_pos[self.ifname - 1] == -1 {
                    break;
                }
                if self.ifname < nfnames && self.start_pos[self.ifname] == -1 {
                    break;
                }
            }
            let mut seek: Option<(i64, Option<u32>)> = None;
            if files.has_line(0) {
                let chr = files_line_chrom(files, 0).into_bytes();
                let pos = files.line(0).expect("reader 0 has a line").pos();
                seek_to(files, &chr, pos)?;
                seek = Some((pos, state.out_hdr.contig_id(&chr)));
            } else if new_file {
                files.seek_start();
            }
            while next_line(files)? {
                let mut site = if files.reader_count() > 1 {
                    Site::Overlap
                } else {
                    Site::Single
                };
                if !files.has_line(0) {
                    if files.region_done(0) {
                        self.flush(state, files, writer)?;
                        files.remove_reader(0);
                        site = Site::Single;
                    } else if self.opts.warn {
                        let ir = active_index(files).expect("a reader has a line");
                        let pos = files.line(ir).expect("active").pos();
                        let chr = files_line_chrom(files, ir);
                        self.warn_drop(&chr, pos);
                        continue;
                    } else if !self.opts.force {
                        let ir = active_index(files).expect("a reader has a line");
                        let pos = files.line(ir).expect("active").pos();
                        bail!(
                            "Error: The --ligate option is intended for VCFs with perfect overlap, the site {}:{} breaks the assumption",
                            files_line_chrom(files, ir),
                            pos + 1
                        );
                    }
                }
                let ir = active_index(files).expect("a reader has a line");
                let chr = files_line_chrom(files, ir).into_bytes();
                let pos = files.line(ir).expect("active").pos();
                // This can happen after a seek: an indel may start before the position sought.
                if let Some((seek_pos, Some(seek_chr))) = seek
                    && seek_pos > pos
                    && state.out_hdr.contig_id(&chr) == Some(seek_chr)
                {
                    continue;
                }
                seek = None;
                let mut must_seek = false;
                while self.ifname < nfnames
                    && self.start_pos[self.ifname] != -1
                    && pos >= self.start_pos[self.ifname]
                {
                    must_seek = true;
                    let f = &state.fnames[self.ifname];
                    if let Err(e) = files.add_reader_again(f) {
                        bail!("Failed to open {f}: {e}");
                    }
                    self.ifname += 1;
                }
                if must_seek {
                    seek_to(files, &chr, pos)?;
                    seek = Some((pos, state.out_hdr.contig_id(&chr)));
                    continue;
                }
                if files.reader_count() > 1
                    && !files.has_line(1)
                    && !files.region_done(1)
                    && !self.opts.force
                {
                    if self.opts.warn {
                        self.warn_drop(&String::from_utf8_lossy(&chr), pos);
                    } else {
                        bail!(
                            "Error: The --ligate option is intended for VCFs with perfect overlap, the site {}:{} breaks the assumption",
                            String::from_utf8_lossy(&chr),
                            pos + 1
                        );
                    }
                    continue;
                }
                self.push(state, files, writer, site)?;
            }
            if let Some(e) = files.take_pending_error() {
                bail!("Error: {e}");
            }
            if files.reader_count() > 0 {
                self.flush(state, files, writer)?;
                while files.reader_count() > 0 {
                    files.remove_reader(0);
                }
            }
        }
        Ok(())
    }

    /// bcftools' `phase_update`: swaps the haplotypes of the samples marked for swapping.
    fn phase_update(&mut self, hdr: &Header, rec: &mut Record) {
        let Ok(n) = rec.get_genotypes(hdr, &mut self.gta) else {
            return;
        };
        if n == 0 {
            return;
        }
        for (i, &swap) in self.swap_phase.iter().enumerate() {
            if !swap || 2 * i + 1 >= n {
                continue;
            }
            let (g0, g1) = (Genotype(self.gta[2 * i]), Genotype(self.gta[2 * i + 1]));
            if g0.is_missing() || g1.is_vector_end() || !g1.is_phased() {
                continue;
            }
            self.gta[2 * i] = g1.0;
            self.gta[2 * i + 1] = g0.0 | 1;
        }
        let _ = rec.update_genotypes(hdr, &self.gta[..n]);
    }

    /// Writes PS (unless compact and unchanged) and the record, checking the order.
    fn write(
        &mut self,
        state: &Concat,
        writer: &mut Writer,
        rec: &mut Record,
        tag: &str,
        func: &'static str,
    ) -> anyhow::Result<()> {
        if !self.opts.compact_ps || self.phase_set_changed {
            let _ = rec.update_format_i32(&state.out_hdr, b"PS", &self.phase_set);
            self.phase_set_changed = false;
        }
        let written = LastWrite::Named(func);
        state.last_write.set(written);
        if let Err(e) = writer.write_record(&state.out_hdr, rec) {
            return Err(output::write_failure(e, || {
                written.failure(&state.output_name)
            }));
        }
        if rec.pos() < self.prev_pos_check {
            bail!(
                "FIXME, disorder: {} vs {}  [{tag}]",
                site(&state.out_hdr, rec),
                self.prev_pos_check + 1
            );
        }
        self.prev_pos_check = rec.pos();
        Ok(())
    }

    /// bcftools' `phased_flush`: decides per sample whether the second file's phase is
    /// swapped relative to the first from the overlap, writes the first half of the overlap
    /// from the first file and the second half from the second.
    fn flush(
        &mut self,
        state: &mut Concat,
        files: &mut SyncedReaders,
        writer: &mut Writer,
    ) -> anyhow::Result<()> {
        if self.nbuf == 0 {
            return Ok(());
        }
        let nsmpl = state.out_hdr.sample_count();
        if files.reader_count() < 2 {
            bail!("[phased_flush] Error: the overlapping files were closed");
        }
        let (ahdr, bhdr) = (files.header(0), files.header(1));
        let mut i = 0;
        while i < self.nbuf {
            if self.buf_mask[i / 2] != 3 {
                i += 2;
                continue;
            }
            let (a, b) = self.buf.split_at_mut(i + 1);
            let arec = &mut a[i];
            let brec = &mut b[0];
            let Ok(na) = arec.get_genotypes(ahdr, &mut self.gta) else {
                warn_gt_absent(&mut self.gt_absent_warned, ahdr, arec);
                i += 2;
                continue;
            };
            if na != 2 * nsmpl {
                i += 2;
                continue;
            }
            let Ok(nb) = brec.get_genotypes(bhdr, &mut self.gtb) else {
                warn_gt_absent(&mut self.gt_absent_warned, bhdr, brec);
                i += 2;
                continue;
            };
            if nb != 2 * nsmpl {
                i += 2;
                continue;
            }
            // The phase is compared by the alleles the indexes name, so that `C,G` with `1|2`
            // and `G,C` with `2|1` agree; bcftools compares the indexes. An index past the
            // alleles is compared as itself.
            arec.unpack(Unpack::STR);
            brec.unpack(Unpack::STR);
            let (arec, brec) = (&*arec, &*brec);
            fn allele(rec: &Record, g: Genotype) -> Result<&[u8], i32> {
                usize::try_from(g.allele())
                    .ok()
                    .filter(|&k| k < rec.allele_count())
                    .map(|k| rec.allele(k))
                    .ok_or(g.allele())
            }
            for j in 0..nsmpl {
                let ga = [Genotype(self.gta[2 * j]), Genotype(self.gta[2 * j + 1])];
                let gb = [Genotype(self.gtb[2 * j]), Genotype(self.gtb[2 * j + 1])];
                if ga[1].is_vector_end() || gb[1].is_vector_end() {
                    continue;
                }
                if ga.iter().chain(&gb).any(|g| g.is_missing()) {
                    continue;
                }
                if !ga[1].is_phased() || !gb[1].is_phased() {
                    continue;
                }
                let ([a0, a1], [b0, b1]) =
                    (ga.map(|g| allele(arec, g)), gb.map(|g| allele(brec, g)));
                if a0 == a1 || b0 == b1 {
                    continue;
                }
                if a0 == b0 && a1 == b1 {
                    if self.swap_phase[j] {
                        self.nmism[j] += 1;
                    } else {
                        self.nmatch[j] += 1;
                    }
                }
                if a0 == b1 && a1 == b0 {
                    if self.swap_phase[j] {
                        self.nmatch[j] += 1;
                    } else {
                        self.nmism[j] += 1;
                    }
                }
            }
            i += 2;
        }
        let mut i = 0;
        while i < self.nbuf / 2 {
            let mask = self.buf_mask[i / 2];
            let (k, hdr) = if mask & 1 != 0 {
                (i, files.header_mut(0))
            } else {
                (i + 1, files.header_mut(1))
            };
            let mut rec = std::mem::take(&mut self.buf[k]);
            let chunk = if mask & 1 != 0 {
                Chunk::First
            } else {
                Chunk::Second
            };
            let result = self.write_buffered(state, writer, &mut rec, hdr, chunk, "1");
            self.buf[k] = rec;
            result?;
            i += 2;
        }
        self.nswap = 0;
        for j in 0..nsmpl {
            if self.nmatch[j] >= self.nmism[j] {
                self.swap_phase[j] = false;
            } else {
                self.swap_phase[j] = true;
                self.nswap += 1;
            }
            if self.nmatch[j] != 0 && self.nmism[j] != 0 {
                // Entropy-inspired quality. The factor 0.7 shifts and scales to (0,1).
                let f = f64::from(self.nmatch[j]) / f64::from(self.nmatch[j] + self.nmism[j]);
                self.phase_qual[j] =
                    (99.0 * (0.7 + f * f.ln() + (1.0 - f) * (1.0 - f).ln()) / 0.7) as i32;
            } else {
                self.phase_qual[j] = 99;
            }
            self.nmatch[j] = 0;
            self.nmism[j] = 0;
        }
        let mut pq_printed = false;
        while i < self.nbuf {
            let mask = self.buf_mask[i / 2];
            let (k, hdr) = if mask & 2 != 0 {
                (i + 1, files.header_mut(1))
            } else {
                (i, files.header_mut(0))
            };
            let mut rec = std::mem::take(&mut self.buf[k]);
            let result = (|| -> anyhow::Result<()> {
                state.translate(&mut rec, hdr)?;
                if !pq_printed && mask == 3 {
                    let _ = rec.update_format_i32(&state.out_hdr, b"PQ", &self.phase_qual);
                    pq_printed = true;
                    for j in 0..nsmpl {
                        if self.phase_qual[j] < self.opts.min_pq {
                            self.phase_set[j] = (rec.pos() + 1) as i32;
                            self.phase_set_changed = true;
                        } else if self.opts.compact_ps {
                            self.phase_set[j] = INT32_MISSING;
                        }
                    }
                }
                if self.nswap > 0 {
                    self.phase_update(&state.out_hdr, &mut rec);
                }
                self.write(state, writer, &mut rec, "2", "phased_flush")
            })();
            self.buf[k] = rec;
            result?;
            i += 2;
        }
        self.nbuf = 0;
        Ok(())
    }

    /// Writes a record of the first half of the overlap.
    fn write_buffered(
        &mut self,
        state: &Concat,
        writer: &mut Writer,
        rec: &mut Record,
        hdr: &mut Header,
        chunk: Chunk,
        tag: &str,
    ) -> anyhow::Result<()> {
        state.translate(rec, hdr)?;
        if self.nswap > 0 && chunk == Chunk::First {
            self.phase_update(&state.out_hdr, rec);
        }
        self.write(state, writer, rec, tag, "phased_flush")
    }

    /// bcftools' `phased_push`: writes a record outside overlaps, or buffers the pair of an
    /// overlap.
    fn push(
        &mut self,
        state: &mut Concat,
        files: &mut SyncedReaders,
        writer: &mut Writer,
        kind: Site,
    ) -> anyhow::Result<()> {
        let has_a = files.has_line(0);
        let has_b = files.reader_count() > 1 && files.has_line(1);
        for (i, has) in [(0, has_a), (1, has_b)] {
            if !has {
                continue;
            }
            let fname = files.file_name(i).to_owned();
            let (rec, hdr) = files.line_with_header(i).expect("the reader has a line");
            if !rec.errors().is_empty() {
                bail!("Parse error at {}, cannot proceed: {fname}", site(hdr, rec));
            }
        }
        let first = if has_a { 0 } else { 1 };
        let chr = files_line_chrom(files, first).into_bytes();
        let pos = files.line(first).expect("the reader has a line").pos();
        let chr_id = state.out_hdr.contig_id(&chr).map_or(-1, i64::from);
        if self.prev_chr < 0 || self.prev_chr != chr_id {
            if self.prev_chr >= 0 {
                self.flush(state, files, writer)?;
            }
            self.phase_set
                .iter_mut()
                .for_each(|p| *p = (pos + 1) as i32);
            self.phase_set_changed = true;
            let seen = usize::try_from(chr_id)
                .ok()
                .and_then(|c| state.seen_seq.get_mut(c));
            let Some(seen) = seen else {
                return Err(ExitOne.into());
            };
            if *seen {
                bail!(
                    "The chromosome block {} is not contiguous",
                    String::from_utf8_lossy(&chr)
                );
            }
            *seen = true;
            self.prev_chr = chr_id;
            self.prev_pos_check = -1;
        }
        if kind == Site::Single {
            let fname = files.file_name(0).to_owned();
            let (rec, hdr) = files
                .line_with_header_mut(0)
                .expect("the first reader has a line");
            state.translate(rec, hdr)?;
            if self.nswap > 0 {
                self.phase_update(&state.out_hdr, rec);
            }
            if !self.opts.compact_ps || self.phase_set_changed {
                let _ = rec.update_format_i32(&state.out_hdr, b"PS", &self.phase_set);
                self.phase_set_changed = false;
            }
            let written = LastWrite::Named("phased_push");
            state.last_write.set(written);
            if let Err(e) = writer.write_record(&state.out_hdr, rec) {
                return Err(output::write_failure(e, || {
                    written.failure(&state.output_name)
                }));
            }
            if rec.pos() < self.prev_pos_check {
                bail!(
                    "FIXME, disorder: {} in {fname} vs {} written  [3]",
                    site(&state.out_hdr, rec),
                    self.prev_pos_check + 1
                );
            }
            self.prev_pos_check = rec.pos();
            return Ok(());
        }
        self.nbuf += 2;
        while self.buf.len() < self.nbuf {
            self.buf.push(Record::new());
        }
        if self.buf_mask.len() < self.nbuf / 2 {
            self.buf_mask.resize(self.nbuf / 2, 0);
        }
        if has_a {
            files.swap_line(0, &mut self.buf[self.nbuf - 2]);
        }
        if has_b {
            files.swap_line(1, &mut self.buf[self.nbuf - 1]);
        }
        self.buf_mask[self.nbuf / 2 - 1] = u8::from(has_a) | (u8::from(has_b) << 1);
        Ok(())
    }
}

/// Warns once that a record has no GT.
fn warn_gt_absent(warned: &mut bool, hdr: &Header, rec: &Record) {
    if !*warned {
        eprintln!(
            "GT is not present at {}. (This warning is printed only once.)",
            site(hdr, rec)
        );
        *warned = true;
    }
}

/// The contig name of reader `i`'s current line.
fn files_line_chrom(files: &SyncedReaders, i: usize) -> String {
    match files.line(i) {
        Some(rec) => String::from_utf8_lossy(rec.chrom_safe(files.header(i))).into_owned(),
        None => String::new(),
    }
}

/// Positions the readers at `chr:pos` (`bcf_sr_seek`).
fn seek_to(files: &mut SyncedReaders, chr: &[u8], pos: i64) -> anyhow::Result<()> {
    if let Err(e) = files.seek(chr, pos) {
        bail!("Error: {e}");
    }
    Ok(())
}
