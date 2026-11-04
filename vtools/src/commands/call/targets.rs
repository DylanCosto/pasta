//! `-T` targets of `bcftools call`: exact positions, with alleles for `-C alleles`, and
//! the lines of `-i/--insert-missed` (bcftools' `tgt_parse`, `next_line`, and `tgt_flush`).

use std::collections::VecDeque;

use anyhow::bail;
use pasta::variant::{
    Header, Record, Unpack, encoding::INT32_VECTOR_END, index::split_idx_name, readlist,
    synced::SyncedReader,
};

use super::{
    EofMarker, read_record,
    regidx::{REGIDX_MAX, RegIdx, strtod_u32},
    vcmp::Vcmp,
};

/// The alleles of a target (`tgt_als_t`).
#[derive(Clone, Debug, Default)]
pub struct TargetAlleles {
    alleles: Vec<Vec<u8>>,
    used: bool,
}

/// What a targets file lists.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TargetKind {
    /// `CHROM POS` (`-T` alone).
    Positions,
    /// `CHROM POS REF,ALT` (`-C alleles`).
    Alleles,
}

/// The targets and the state of `-C alleles` and `-i`.
#[derive(Debug)]
pub struct Targets {
    idx: RegIdx<TargetAlleles>,
    kind: TargetKind,
    /// Records read ahead at one position (bcftools' `vcfbuf`).
    buf: VecDeque<Box<Record>>,
    #[expect(
        clippy::vec_box,
        reason = "records are swapped with the reader's boxed buffer, a pointer swap"
    )]
    spare: Vec<Box<Record>>,
    /// The sequence and position of the last record matched (`tgt_itr`).
    cur: Option<(Vec<u8>, u32)>,
    /// The same for the previous record called (`tgt_itr_prev`).
    prev: Option<(Vec<u8>, u32)>,
    /// The line written for targets missing from the input (`missed_line`).
    missed: Record,
    /// The output name, for errors.
    pub out_name: String,
}

/// A target: sequence, 0-based position, alleles.
type Target = (Vec<u8>, u32, TargetAlleles);

/// One line of a targets file.
enum TargetLine {
    Target(Target),
    Comment,
    /// A line `tgt_parse` rejects (the message is printed).
    Malformed,
}

/// A `-T` line: `CHROM POS` and, with `-C alleles`, `REF,ALT` (`tgt_parse`).
fn parse_line(line: &[u8], kind: TargetKind) -> anyhow::Result<TargetLine> {
    let ws = |c: u8| c.is_ascii_whitespace();
    let text = String::from_utf8_lossy(line);
    let s = line.iter().position(|&c| !ws(c)).unwrap_or(line.len());
    if s == line.len() {
        eprintln!("Could not parse the line: {text}");
        return Ok(TargetLine::Malformed);
    }
    if line[s] == b'#' {
        return Ok(TargetLine::Comment);
    }
    let e = line[s..]
        .iter()
        .position(|&c| ws(c))
        .map_or(line.len(), |p| s + p);
    let chr = line[s..e].to_vec();
    if e == line.len() {
        eprintln!("Could not parse the line: {text}");
        return Ok(TargetLine::Malformed);
    }
    let (beg, used) = strtod_u32(&line[e + 1..]);
    if used == 0 {
        eprintln!("Could not parse tab line: {text}");
        return Ok(TargetLine::Malformed);
    }
    if beg == 0 {
        eprintln!("Could not parse tab line, expected 1-based coordinate: {text}");
        return Ok(TargetLine::Malformed);
    }
    let beg = beg - 1;
    let mut als = TargetAlleles::default();
    if kind == TargetKind::Alleles {
        let mut p = e + 1 + used + 1;
        while p < line.len() {
            let start = p;
            while p < line.len() && line[p] != b',' && !ws(line[p]) {
                p += 1;
            }
            als.alleles.push(line[start..p].to_vec());
            let stop = p >= line.len() || ws(line[p]);
            p += 1;
            if stop {
                break;
            }
        }
        if als.alleles.len() < 2 {
            bail!(
                "Unable to parse the -T file; expected CHROM\\tPOS\\tREF,ALT with -C alleles but found instead:\n\t{text}"
            );
        }
    }
    Ok(TargetLine::Target((chr, beg, als)))
}

/// mpileup output has no MNPs, and `<*>` only at non-indel sites (`is_indel` in
/// `vcfcall.c`).
fn is_indel<A: AsRef<[u8]>>(als: &[A]) -> bool {
    if als
        .get(1)
        .is_some_and(|a| a.as_ref().first() == Some(&b'<'))
    {
        return false;
    }
    als.iter().any(|a| {
        let a = a.as_ref();
        a.first() != Some(&b'<') && a.len() > 1
    })
}

impl Targets {
    /// Reads a targets file (`regidx_init` with `tgt_parse`). `None` when the file cannot be
    /// read or has malformed lines: bcftools then calls every site. `FILE##idx##INDEX` reads
    /// `FILE`, as `hts_open` does.
    ///
    /// # Errors
    ///
    /// Returns an error for `-C alleles` lines without alleles.
    pub fn load(path: &str, kind: TargetKind) -> anyhow::Result<Option<Self>> {
        let data = split_idx_name(path).0;
        let Ok(mut reader) = readlist::open_text(std::path::Path::new(data)) else {
            eprintln!(
                "[E::hts_open_format] Failed to open file \"{data}\" : No such file or directory"
            );
            return Ok(None);
        };
        // A failed read (a truncated or corrupt gzip member) fails bcftools' `regidx_init`,
        // and call then crashes on the missing index: an error here.
        let lines = readlist::read_all_lines_from(&mut reader)
            .map_err(|e| anyhow::anyhow!("Error: failed to read the targets file {data}: {e}"))?;
        let mut idx = RegIdx::default();
        // regidx stops at the first empty line (hts_getline returns 0 for it).
        for line in lines.iter().take_while(|l| !l.is_empty()) {
            match parse_line(line, kind)? {
                TargetLine::Target((chr, beg, als)) => idx.push(&chr, beg, beg, als),
                TargetLine::Comment => {}
                TargetLine::Malformed => return Ok(None),
            }
        }
        idx.finish();
        Ok(Some(Self {
            idx,
            kind,
            buf: VecDeque::new(),
            spare: Vec::new(),
            cur: None,
            prev: None,
            missed: Record::new(),
            out_name: "-".to_owned(),
        }))
    }

    /// Whether the targets carry alleles (`-C alleles`).
    pub fn has_alleles(&self) -> bool {
        self.kind == TargetKind::Alleles
    }

    /// Whether a target starts exactly at the record's position.
    pub fn matches(&self, rhdr: &Header, rec: &Record) -> bool {
        let chr = rec.chrom_safe(rhdr);
        let pos = rec.pos() as u32;
        self.idx.overlaps(chr, pos, pos).any(|h| h.beg == pos)
    }

    /// `next_line` with `-C alleles`: buffers the records at one position and returns the
    /// one whose alleles best match an unused target, with that target's alleles (`None`
    /// when no target is left for the record).
    pub fn next_with_alleles(
        &mut self,
        reader: &mut SyncedReader,
        subset: Option<&[i32]>,
        out: &mut Box<Record>,
        eof: EofMarker,
    ) -> anyhow::Result<Option<Option<Vec<Vec<u8>>>>> {
        let same = |a: &Record, b: &Record| a.rid() == b.rid() && a.pos() == b.pos();
        let full = match self.buf.len() {
            0 | 1 => false,
            n => !same(&self.buf[0], &self.buf[n - 1]),
        };
        if !full {
            loop {
                let mut rec = self.spare.pop().unwrap_or_default();
                if !read_record(reader, &mut rec, eof)? {
                    self.spare.push(rec);
                    break;
                }
                if !self.matches(reader.header(), &rec) {
                    self.spare.push(rec);
                    continue;
                }
                if let Some(map) = subset {
                    rec.subset(map);
                }
                rec.unpack(Unpack::STR);
                self.buf.push_back(rec);
                let n = self.buf.len();
                if !same(&self.buf[0], &self.buf[n - 1]) {
                    break;
                }
            }
        }
        let nbuf = self.buf.len();
        if nbuf == 0 {
            return Ok(None);
        }
        let mut n = nbuf;
        while n > 1 {
            if same(&self.buf[0], &self.buf[n - 1]) {
                break;
            }
            n -= 1;
        }

        let rec0 = &self.buf[0];
        let chr = rec0.chrom_safe(reader.header()).to_vec();
        let pos = rec0.pos() as u32;
        let hits: Vec<_> = self.idx.overlaps(&chr, pos, pos).collect();
        self.cur = hits.first().map(|h| (chr.clone(), h.beg));
        let mut best: Option<(i32, usize, usize)> = None;
        let mut vcmp = Vcmp::default();
        for i in 0..n {
            let rec = &self.buf[i];
            let rec_als: Vec<&[u8]> = rec.alleles().collect();
            let rec_indel = if is_indel(&rec_als) { 1 } else { -1 };
            for (hi, h) in hits.iter().enumerate() {
                if i64::from(h.beg) != rec.pos() {
                    continue;
                }
                let als = self.idx.payload(h);
                if als.used {
                    continue;
                }
                let mut nmatch = 0;
                if vcmp.set_ref(rec_als[0], &als.alleles[0]) {
                    nmatch += 1;
                    if rec_als.len() > 1 && als.alleles.len() > 1 {
                        for a in &als.alleles[1..] {
                            if vcmp.find_allele(&rec_als[1..], a).is_some() {
                                nmatch += 1;
                            }
                        }
                    }
                }
                let als_indel = if is_indel(&als.alleles) { 1 } else { -1 };
                nmatch *= rec_indel * als_indel;
                if best.is_none_or(|b| nmatch > b.0) {
                    best = Some((nmatch, hi, i));
                }
            }
        }
        let (tgt, ibuf) = match best {
            Some((_, hi, i)) => {
                let t = self.idx.payload_mut(&hits[hi]);
                t.used = true;
                (Some(t.alleles.clone()), i)
            }
            None => (None, 0),
        };
        let rec = self.buf.remove(ibuf).expect("the buffer has the record");
        let old = std::mem::replace(out, rec);
        self.spare.push(old);
        Ok(Some(tgt))
    }

    /// `tgt_flush`: writes the unused targets before the record (all of them at the end).
    pub fn flush(
        &mut self,
        hdr: &Header,
        rec: Option<&Record>,
        write: &mut dyn FnMut(&mut Record) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        match rec {
            Some(rec) => {
                let chr = rec.chrom_safe(hdr).to_vec();
                // bcftools passes `rec->pos-1` as a `uint32_t`, which wraps for a record at
                // POS 1 and writes every later target of the sequence as missing before it
                // (their calls are then dropped); nothing comes before POS 1.
                let end = (rec.pos() > 0).then(|| (rec.pos() - 1) as u32);
                // `prev` is replaced below.
                match self.prev.take() {
                    None => {
                        if let Some(end) = end {
                            self.flush_region(hdr, &chr, 0, end, write)?;
                        }
                    }
                    Some((pchr, pbeg)) if pchr != chr => {
                        self.flush_region(hdr, &pchr, pbeg.wrapping_add(1), REGIDX_MAX, write)?;
                        if let Some(end) = end {
                            self.flush_region(hdr, &chr, 0, end, write)?;
                        }
                    }
                    Some((pchr, pbeg)) => {
                        if let Some(end) = end {
                            self.flush_region(hdr, &pchr, pbeg, end, write)?;
                        }
                    }
                }
                self.prev = self.cur.clone();
            }
            None => {
                if let Some((pchr, pbeg)) = self.prev.take() {
                    self.flush_region(hdr, &pchr, pbeg, REGIDX_MAX, write)?;
                }
                let names: Vec<Vec<u8>> = self.idx.seq_names().to_vec();
                for chr in names {
                    self.flush_region(hdr, &chr, 0, REGIDX_MAX, write)?;
                }
            }
        }
        Ok(())
    }

    /// `tgt_flush_region`.
    fn flush_region(
        &mut self,
        hdr: &Header,
        chr: &[u8],
        beg: u32,
        end: u32,
        write: &mut dyn FnMut(&mut Record) -> anyhow::Result<()>,
    ) -> anyhow::Result<()> {
        let hits: Vec<_> = self
            .idx
            .overlaps(chr, beg, end)
            .filter(|h| h.beg >= beg)
            .collect();
        for h in hits {
            if self.idx.payload(&h).used {
                continue;
            }
            // A target on a sequence missing from the header cannot be written.
            let Some(rid) = hdr.contig_id(chr) else {
                eprintln!("[E::vcf_format] Invalid BCF, CONTIG id=-1 not present in the header");
                bail!(
                    "[tgt_flush_region] Error: failed to write to {}",
                    self.out_name
                );
            };
            let nsmpl = hdr.sample_count();
            let line = &mut self.missed;
            if line.allele_count() == 0 {
                // `init_missed_line`: missing haploid genotypes, missing QUAL.
                if nsmpl > 0 {
                    let mut gts = vec![0i32; nsmpl * 2];
                    for g in gts.as_chunks_mut::<2>().0 {
                        g[1] = INT32_VECTOR_END;
                    }
                    line.update_genotypes(hdr, &gts)?;
                }
                line.set_qual(None);
            }
            line.set_rid(rid as i32);
            line.set_pos(i64::from(h.beg));
            let t = self.idx.payload_mut(&h);
            line.update_alleles(hdr, &t.alleles);
            t.used = true;
            write(line)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    use clap::{Args as _, Command, FromArgMatches as _};

    /// Runs `call -m -i -C alleles` on `vcf` with targets `tgt`; returns CHROM to QUAL of
    /// each record.
    fn call_targets(name: &str, vcf: &str, tgt: &str) -> anyhow::Result<Vec<String>> {
        let dir = std::env::temp_dir().join(format!("vtools-call-{name}-{}", std::process::id()));
        fs::create_dir_all(&dir)?;
        let input = dir.join("in.vcf");
        let targets = dir.join("tgt.tsv");
        let out = dir.join("out.vcf");
        fs::write(&input, vcf)?;
        fs::write(&targets, tgt)?;
        let argv = [
            "call",
            "--no-version",
            "-mi",
            "-C",
            "alleles",
            "-T",
            targets.to_str().expect("UTF-8 path"),
            "-o",
            out.to_str().expect("UTF-8 path"),
            input.to_str().expect("UTF-8 path"),
        ];
        let matches =
            super::super::Args::augment_args(Command::new("call")).try_get_matches_from(argv)?;
        super::super::run(super::super::Args::from_arg_matches(&matches)?, &matches)?;
        let text = fs::read_to_string(&out)?;
        for f in [&input, &targets, &out] {
            fs::remove_file(f)?;
        }
        fs::remove_dir(&dir)?;
        Ok(text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').take(6).collect::<Vec<_>>().join(" "))
            .collect())
    }

    #[test]
    fn a_record_at_pos_1_keeps_the_later_targets() -> anyhow::Result<()> {
        // bcftools flushes the targets up to POS-1 as a `uint32_t`, so a record at POS 1
        // wrote every later target of its sequence as missing and dropped their calls.
        let vcf = "##fileformat=VCFv4.2\n##contig=<ID=1>\n##contig=<ID=2>\n\
                   ##INFO=<ID=I16,Number=16,Type=Float,Description=\"I16\">\n\
                   ##INFO=<ID=QS,Number=R,Type=Float,Description=\"QS\">\n\
                   ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let site = |chr: u8, pos: u32, r: char, a: char| {
            format!(
                "{chr}\t{pos}\t.\t{r}\t{a},<*>\t0\t.\tI16=5,5,5,5,200,4000,200,4000,600,18000,\
                 600,18000,100,1000,100,1000;QS=0.5,0.5,0\tPL\t200,0,200,255,255,255\n"
            )
        };
        // The first record, the first record of a later sequence, and a second record at
        // POS 1; a target with no record (2:5) is still written as missing, in order.
        let input = [
            site(1, 1, 'A', 'C'),
            site(1, 1, 'A', 'G'),
            site(1, 10, 'G', 'T'),
            site(2, 1, 'A', 'C'),
            site(2, 10, 'G', 'T'),
        ]
        .concat();
        let tgt = "1\t1\tA,C\n1\t1\tA,G\n1\t10\tG,T\n2\t1\tA,C\n2\t5\tC,T\n2\t10\tG,T\n";
        let lines = call_targets("pos1", &format!("{vcf}{input}"), tgt)?;
        assert_eq!(
            lines,
            [
                "1 1 . A C 167.406",
                "1 1 . A G 167.406",
                "1 10 . G T 167.406",
                "2 1 . A C 167.406",
                "2 5 . C T .",
                "2 10 . G T 167.406"
            ]
        );
        Ok(())
    }

    #[test]
    fn idx_names_read_the_data_file() {
        let (dir, name) =
            crate::raw_input::write_idx_named("targets-idx", "t.tsv", "chr1\t100\tA,G\n");
        let t = Targets::load(&name, TargetKind::Alleles);
        let _ = std::fs::remove_dir_all(&dir);
        let t = t.expect("the lines parse").expect("the data file is read");
        assert_eq!(t.idx.overlaps(b"chr1", 99, 99).count(), 1);
    }
}
