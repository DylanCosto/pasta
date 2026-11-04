//! Processing records for `+split-vep`: splitting the annotation into transcripts and
//! subfields, selecting transcripts, setting the INFO tags, and writing VCF/BCF or text
//! (bcftools' `process_record`, `filter_and_output`, and their helpers).

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use anyhow::bail;
use flate2::{Compression, write::GzEncoder};
use pasta::variant::{
    Header, Record,
    encoding::{FLOAT_MISSING, INT32_MISSING},
    io::OutputType,
    number,
    synced::SyncedReader,
};

use super::{Annot, Args, ColType, SelectTr, Severity, SplitVep, TrOp};
use crate::{
    filter::FilterLogic,
    output::{self, OutputSetting},
    plugins::{Sink, Source},
    raw_input::strerror,
    version,
};

/// `-s :any`: no consequence filtering (bcftools' `SELECT_CSQ_ANY`).
pub(super) const SELECT_CSQ_ANY: i32 = -1;

/// The annotation of the current record split into transcripts and subfields (bcftools'
/// `cols_tr` and `cols_csq`).
#[derive(Debug, Default)]
pub(super) struct Csq {
    /// The annotation text.
    buf: Vec<u8>,
    /// Transcript ranges into `buf`.
    trs: Vec<(usize, usize)>,
    /// Subfield ranges of each transcript.
    fields: Vec<Vec<(usize, usize)>>,
    /// The transcripts in use (after `-g`).
    ncols: usize,
}

impl Csq {
    /// Splits `buf` (bcftools' `cols_split` by `,` then `|`); returns whether a transcript
    /// has too few subfields to hold the consequence.
    fn split(&mut self, csq_idx: usize) -> bool {
        self.trs.clear();
        let mut start = 0;
        for i in memchr::memchr_iter(b',', &self.buf) {
            self.trs.push((start, i));
            start = i + 1;
        }
        self.trs.push((start, self.buf.len()));
        if self.fields.len() < self.trs.len() {
            self.fields.resize_with(self.trs.len(), Vec::new);
        }
        let mut need_fix = false;
        for (k, &(b, e)) in self.trs.iter().enumerate() {
            let f = &mut self.fields[k];
            f.clear();
            let mut s = b;
            for i in memchr::memchr_iter(b'|', &self.buf[b..e]) {
                f.push((s, b + i));
                s = b + i + 1;
            }
            f.push((s, e));
            need_fix |= csq_idx >= f.len();
        }
        self.ncols = self.trs.len();
        need_fix
    }

    fn field(&self, tr: usize, idx: usize) -> &[u8] {
        let (b, e) = self.fields[tr][idx];
        &self.buf[b..e]
    }
}

impl Severity {
    /// The severity of one consequence term, adding unknown terms to the scale (bcftools'
    /// `csq_term_severity`).
    fn term(&mut self, term: &[u8]) -> i32 {
        self.key.clear();
        self.key.extend(term.iter().map(u8::to_ascii_lowercase));
        if let Some(&s) = self.map.get(self.key.as_slice()) {
            return s;
        }
        let key = self.key.clone();
        let known = self
            .scale
            .iter()
            .find(|s| memchr::memmem::find(&key, s).is_some());
        let severity = match known {
            Some(s) => *self.map.get(s).expect("scale terms are mapped"),
            None => self.scale.len() as i32 + 1,
        };
        if known.is_none() {
            eprintln!(
                "Note: assigning a (high) severity score to a new consequence, use -S to override: {} -> {severity}",
                String::from_utf8_lossy(&key)
            );
        }
        self.map.insert(key.clone(), severity);
        self.scale.push(key);
        severity
    }

    /// The smallest and largest severity of `a&b&c` consequences, or the one equal to
    /// `exact` (bcftools' `csq_to_severity`).
    fn of(&mut self, csq: &[u8], exact: Option<i32>) -> (i32, i32) {
        let (mut min, mut max) = (i32::MAX, -1);
        // bcftools' loop stops at the end of the string: an empty string, or the empty term
        // after a trailing `&`, is not looked up.
        if csq.is_empty() {
            return (min, max);
        }
        let nterms = memchr::memchr_iter(b'&', csq).count() + 1 - usize::from(csq.ends_with(b"&"));
        for term in csq.split(|&c| c == b'&').take(nterms) {
            let s = self.term(term);
            match exact {
                None => {
                    min = min.min(s);
                    max = max.max(s);
                }
                Some(e) if s == e => return (s, s),
                Some(_) => {}
            }
        }
        (min, max)
    }
}

/// Parses comma-separated floats (bcftools' `parse_array_real`).
fn parse_array_real(tag: &str, s: &[u8], out: &mut Vec<f32>, warned: &mut bool) {
    out.clear();
    let mut bp = 0;
    while bp < s.len() {
        let (v, ep) = number::strtod(s, bp);
        let mut ep = ep;
        if ep == bp {
            let missing_dot =
                s.get(ep) == Some(&b'.') && matches!(s.get(ep + 1), None | Some(b','));
            if !*warned && !missing_dot {
                eprintln!(
                    "Warning: Could not parse, not a numeric list {tag}=\"{}\", check the -c and --columns-types options.\n         This message is printed only once.",
                    String::from_utf8_lossy(s)
                );
                *warned = true;
            }
            out.push(f32::from_bits(FLOAT_MISSING));
        } else {
            out.push(v as f32);
        }
        while ep < s.len() && s[ep] != b',' {
            ep += 1;
        }
        bp = if ep < s.len() { ep + 1 } else { ep };
    }
}

/// Parses comma-separated integers (bcftools' `parse_array_int32`).
fn parse_array_int32(tag: &str, s: &[u8], out: &mut Vec<i32>, warned: &mut bool) {
    out.clear();
    let mut bp = 0;
    while bp < s.len() {
        let (v, ep) = number::strtol(s, bp);
        let mut ep = ep;
        if ep == bp {
            let missing_dot =
                s.get(ep) == Some(&b'.') && matches!(s.get(ep + 1), None | Some(b','));
            if !*warned && !missing_dot {
                eprintln!(
                    "Warning: Could not parse, not a numeric list {tag}=\"{}\", check the -c and --columns-types options.\n         This message is printed only once.",
                    String::from_utf8_lossy(s)
                );
                *warned = true;
            }
            out.push(INT32_MISSING);
        } else {
            // C assigns the long to an int32_t, keeping the low 32 bits.
            out.push(v as i32);
        }
        while ep < s.len() && s[ep] != b',' {
            ep += 1;
        }
        bp = if ep < s.len() { ep + 1 } else { ep };
    }
}

/// Where records or text go.
enum Out {
    /// `-f` text, plain or gzip-compressed.
    Text {
        writer: Box<dyn Write>,
        name: String,
    },
    /// VCF/BCF records.
    Vcf { sink: Box<Sink>, name: String },
}

/// Opens the `-f` text output: gzip with `-Oz`/`-Ob`, plain otherwise (bcftools'
/// `bgzf_open(..., "wg"/"wu")`).
fn open_text(args: &Args, name: &str) -> anyhow::Result<Box<dyn Write>> {
    let settings: Vec<OutputSetting> = args
        .output_type
        .iter()
        .map(|t| OutputSetting::Type(t.clone()))
        .collect();
    let compressed = matches!(
        output::requested_type(&settings)?,
        OutputType::VcfBgzf | OutputType::Bcf
    );
    let inner: Box<dyn Write> = if name == "-" {
        Box::new(output::stdout())
    } else {
        Box::new(
            File::create(Path::new(name))
                .map_err(|e| anyhow::anyhow!("Could not open {name}: {}", strerror(&e)))?,
        )
    };
    Ok(if compressed {
        Box::new(BufWriter::with_capacity(
            1 << 16,
            GzEncoder::new(inner, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(1 << 16, inner))
    })
}

/// The main loop of `run`: opens the output and processes every record.
pub(super) fn run(
    mut sv: SplitVep,
    args: &Args,
    mut reader: SyncedReader,
    mut hdr_out: Header,
    threads: usize,
) -> anyhow::Result<()> {
    let name = args
        .output
        .as_deref()
        .map_or_else(|| "-".to_owned(), |p| p.display().to_string());
    let mut out = if sv.has_format {
        let mut writer = open_text(args, &name)?;
        if args.print_header > 0 {
            let mut line = Vec::new();
            if let Some(c) = &mut sv.convert {
                c.header_line(&hdr_out, &mut line);
            }
            writer.write_all(&line).map_err(|e| write_error(e, &name))?;
        }
        Out::Text { writer, name }
    } else {
        if !args.no_version {
            version::append_version(&mut hdr_out, "split-vep")?;
        }
        let settings: Vec<OutputSetting> = args
            .output_type
            .iter()
            .map(|t| OutputSetting::Type(t.clone()))
            .collect();
        let sink = Sink::create(
            &settings,
            args.output.as_deref(),
            args.write_index.last().map(String::as_str),
            &hdr_out,
            &name,
            threads,
        )?
        .with_write_error(format!("Failed to write to {name}"));
        Out::Vcf {
            sink: Box::new(sink),
            name,
        }
    };
    let mut source = Source::new(&mut reader, threads, sv.max_unpack);
    let mut rec = source.take();
    let result = loop {
        match source.next(&mut reader, &mut rec) {
            Ok(true) => {}
            Ok(false) => break Ok(()),
            Err(e) => break Err(anyhow::anyhow!("Error: {e}")),
        }
        if let Err(e) =
            sv.process_record(reader.header(), &hdr_out, &mut rec, &mut out, &mut source)
        {
            break Err(e);
        }
    };
    let result = result.and_then(|()| match reader.take_pending_error() {
        Some(e) => Err(anyhow::anyhow!("Error: {e}")),
        None => Ok(()),
    });
    match out {
        Out::Text { mut writer, name } => {
            result?;
            writer.flush().map_err(|e| write_error(e, &name))?;
        }
        Out::Vcf { sink, name } => {
            if let Err(e) = result {
                sink.abandon(&hdr_out);
                return Err(e);
            }
            sink.finish(&hdr_out, &name)?;
        }
    }
    Ok(())
}

/// A write error: a closed pipe passes through, for `main` to die of SIGPIPE.
fn write_error(e: io::Error, name: &str) -> anyhow::Error {
    if e.kind() == io::ErrorKind::BrokenPipe {
        return anyhow::Error::new(e);
    }
    anyhow::anyhow!("Failed to write to {name}")
}

impl SplitVep {
    /// Splits the annotation, fixing commas inside subfields when the counts do not add up
    /// (bcftools' `split_csq_fields`).
    fn split_csq_fields(&mut self, hdr: &Header, rec: &Record) {
        if !self.csq.split(self.csq_idx) {
            return;
        }
        let nfield = self.fields.len() as i64;
        if !self.warned.nfields {
            let mut nfield_diff: i64 = 0;
            for f in &self.csq.fields[..self.csq.trs.len()] {
                let d = f.len() as i64 - nfield;
                if nfield_diff < d.abs() {
                    nfield_diff = d;
                }
            }
            eprintln!(
                "Warning: The number of INFO/{} subfields at {}:{} does not match the header definition,\n         expected {nfield} subfields, found as {} as {}. (This warning is printed only once.)",
                self.vep_tag,
                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                rec.pos() + 1,
                if nfield_diff > 0 { "much" } else { "few" },
                nfield + nfield_diff
            );
        }
        self.warned.nfields = true;
        // LoF_info subfields can contain commas: rely on the number of pipes matching the
        // header definition and replace the other commas with slashes.
        let mut npipe: i64 = 0;
        for c in self.csq.buf.iter_mut().rev() {
            if *c == b'|' {
                npipe += 1;
                continue;
            }
            if *c != b',' {
                continue;
            }
            if npipe != 0 && nfield > 1 && npipe % (nfield - 1) == 0 {
                npipe = 0;
                continue;
            }
            *c = b'/';
        }
        self.csq.split(self.csq_idx);
    }

    /// Keeps only the transcripts of `-g` genes (bcftools' `restrict_csqs_to_genes`).
    fn restrict_to_genes(&mut self) {
        let Some(genes) = &self.genes else {
            return;
        };
        let n = self.csq.trs.len();
        let mut hits = std::mem::take(&mut self.list_tr);
        hits.clear();
        let mut nhit = 0;
        for i in 0..n {
            let f = &self.csq.fields[i];
            let hit = genes
                .fields
                .iter()
                .any(|&idx| idx < f.len() && genes.names.contains(self.csq.field(i, idx)));
            hits.push(usize::from(hit));
            nhit += usize::from(hit);
        }
        if nhit == 0 {
            if !genes.prioritize {
                self.csq.ncols = 0;
            }
            self.list_tr = hits;
            return;
        }
        // Move the hits first, swapping as bcftools does.
        let (mut i, mut j) = (0usize, n as isize - 1);
        while (i as isize) < j {
            if hits[i] == 1 {
                i += 1;
                continue;
            }
            if hits[j as usize] == 0 {
                j -= 1;
                continue;
            }
            self.csq.fields.swap(i, j as usize);
            self.csq.trs.swap(i, j as usize);
            i += 1;
            j -= 1;
        }
        self.csq.ncols = nhit;
        self.list_tr = hits;
    }

    /// The transcripts `-s` selects, into `list_tr`.
    fn select_transcripts(&mut self, hdr: &Header, rec: &Record) -> anyhow::Result<()> {
        self.list_tr.clear();
        let n = self.csq.ncols;
        if n == 0 {
            return Ok(());
        }
        match &self.select_tr {
            SelectTr::Expr(e) => {
                for i in 0..n {
                    let nf = self.csq.fields[i].len();
                    if e.idx >= nf {
                        bail!(
                            "Too few columns at {}:{} .. {} (field {}) >= {nf}",
                            String::from_utf8_lossy(rec.chrom_safe(hdr)),
                            rec.pos() + 1,
                            e.idx,
                            e.field
                        );
                    }
                    let v = self.csq.field(i, e.idx);
                    let keep = match &e.op {
                        TrOp::Eq(x) => v == x.as_slice(),
                        TrOp::Ne(x) => v != x.as_slice(),
                        TrOp::Re(re) => re.is_match(v),
                        TrOp::NotRe(re) => !re.is_match(v),
                    };
                    if keep {
                        self.list_tr.push(i);
                    }
                }
            }
            SelectTr::Worst => {
                let (mut max_sev, mut imax) = (-1, 0);
                for i in 0..n {
                    let nf = self.csq.fields[i].len();
                    if self.csq_idx >= nf {
                        bail!(
                            "Too few columns at {}:{} .. {} ({}) >= {nf}",
                            String::from_utf8_lossy(rec.chrom_safe(hdr)),
                            rec.pos() + 1,
                            self.csq_idx,
                            self.csq_field
                        );
                    }
                    let (b, e) = self.csq.fields[i][self.csq_idx];
                    let (_, max) = self.severity.of(&self.csq.buf[b..e], None);
                    if max_sev < max {
                        imax = i;
                        max_sev = max;
                    }
                }
                self.list_tr.push(imax);
            }
            SelectTr::All => self.list_tr.extend(0..n),
        }
        Ok(())
    }

    /// Whether a consequence passes the `-s :CSQ` severity range (bcftools'
    /// `csq_severity_pass`).
    fn severity_pass(&mut self, csq: (usize, usize)) -> bool {
        if self.min_severity == SELECT_CSQ_ANY && self.max_severity == SELECT_CSQ_ANY {
            return true;
        }
        let exact = (self.min_severity == self.max_severity).then_some(self.min_severity);
        let (min, max) = self.severity.of(&self.csq.buf[csq.0..csq.1], exact);
        !(max < self.min_severity || min > self.max_severity)
    }

    /// The most severe of `a&b&c` (bcftools' `csq_rewrite_worst`).
    fn rewrite_worst(&mut self, range: (usize, usize)) -> (usize, usize) {
        let s = &self.csq.buf[range.0..range.1];
        if memchr::memchr(b'&', s).is_none() {
            return range;
        }
        let (mut worst, mut smax) = (range, -1);
        let mut start = range.0;
        let ends = memchr::memchr_iter(b'&', s)
            .map(|i| range.0 + i)
            .chain(std::iter::once(range.1));
        for end in ends {
            let sev = self.severity.term(&self.csq.buf[start..end]);
            if smax < sev {
                smax = sev;
                worst = (start, end);
            }
            start = end + 1;
        }
        worst
    }

    /// Processes one record (bcftools' `process_record`).
    fn process_record(
        &mut self,
        hdr: &Header,
        hdr_out: &Header,
        rec: &mut Record,
        out: &mut Out,
        source: &mut Source,
    ) -> anyhow::Result<()> {
        let got = match rec.get_info_string(hdr, self.vep_tag.as_bytes()) {
            Ok(s) => {
                let s = memchr::memchr(0, s).map_or(s, |i| &s[..i]);
                self.csq.buf.clear();
                self.csq.buf.extend_from_slice(s);
                !s.is_empty()
            }
            Err(_) => false,
        };
        if !got {
            if !self.drop_sites {
                self.reset_annots();
                let matched = Matched {
                    severity_pass: true,
                    all_missing: true,
                };
                self.filter_and_output(hdr_out, rec, matched, out, source)?;
            }
            return Ok(());
        }
        self.split_csq_fields(hdr, rec);
        self.restrict_to_genes();
        self.select_transcripts(hdr, rec)?;
        self.reset_annots();
        let mut matched = Matched::NONE;
        let list = std::mem::take(&mut self.list_tr);
        for &i in &list {
            let nf = self.csq.fields[i].len();
            if self.csq_idx >= nf {
                self.list_tr = list;
                bail!(
                    "Too few columns at {}:{} .. {} ({}) >= {nf}",
                    String::from_utf8_lossy(rec.chrom_safe(hdr)),
                    rec.pos() + 1,
                    self.csq_idx,
                    self.csq_field
                );
            }
            let csq = self.csq.fields[i][self.csq_idx];
            if !self.severity_pass(csq) {
                continue;
            }
            matched.severity_pass = true;
            for j in 0..self.annots.len() {
                let idx = self.annots[j].idx;
                let range = match idx {
                    None => Some(self.csq.trs[i]),
                    Some(idx) if idx >= nf => {
                        if !self.warned.too_few_fields {
                            eprintln!(
                                "Warning: fewer {} fields than expected at {}:{}, filling with dots. This warning is printed only once.",
                                self.vep_tag,
                                String::from_utf8_lossy(rec.chrom_safe(hdr)),
                                rec.pos() + 1
                            );
                            self.warned.too_few_fields = true;
                        }
                        append(&mut self.annots[j], b".");
                        continue;
                    }
                    Some(idx) => {
                        let r = self.csq.fields[i][idx];
                        if r.0 == r.1 {
                            None
                        } else if idx == self.csq_idx && self.prn_csq_worst {
                            Some(self.rewrite_worst(r))
                        } else {
                            Some(r)
                        }
                    }
                };
                match range {
                    Some((b, e)) => {
                        let v = &self.csq.buf[b..e];
                        append(&mut self.annots[j], v);
                        matched.all_missing = false;
                    }
                    None => append(&mut self.annots[j], b"."),
                }
            }
            if self.duplicate {
                self.filter_and_output(hdr_out, rec, matched, out, source)?;
                self.reset_annots();
                matched = Matched::NONE;
            }
        }
        self.list_tr = list;
        if !matched.severity_pass && self.drop_sites {
            return Ok(());
        }
        if !self.duplicate {
            self.filter_and_output(hdr_out, rec, matched, out, source)?;
        }
        Ok(())
    }

    fn reset_annots(&mut self) {
        for a in &mut self.annots {
            a.value.clear();
        }
    }

    /// Sets the INFO tags, applies `-i/-e`, and writes the record or its text (bcftools'
    /// `filter_and_output`).
    fn filter_and_output(
        &mut self,
        hdr_out: &Header,
        rec: &mut Record,
        matched: Matched,
        out: &mut Out,
        source: &mut Source,
    ) -> anyhow::Result<()> {
        let mut updated = 0;
        for a in &self.annots {
            if a.value.is_empty() {
                continue;
            }
            if let Some(id) = a.id {
                // bcftools ignores the result of these updates.
                let _ = match a.ty {
                    ColType::Real => {
                        parse_array_real(&a.tag, &a.value, &mut self.farr, &mut self.warned.real);
                        rec.update_info_f32_id(hdr_out, id, &self.farr)
                    }
                    ColType::Int => {
                        parse_array_int32(&a.tag, &a.value, &mut self.iarr, &mut self.warned.int);
                        rec.update_info_i32_id(hdr_out, id, &self.iarr)
                    }
                    ColType::Str | ColType::Flag => {
                        rec.update_info_string(hdr_out, a.tag.as_bytes(), Some(&a.value))
                    }
                };
            }
            updated += 1;
        }
        let mut smpl_pass = None;
        if let Some((filter, logic)) = &mut self.filter {
            let pass = filter.test(hdr_out, rec)?;
            if pass != (*logic == FilterLogic::Include) {
                return Ok(());
            }
            smpl_pass = filter.sample_pass();
        }
        match out {
            Out::Text { writer, name } => {
                if self.annots.is_empty() {
                    if !matched.severity_pass {
                        return Ok(());
                    }
                } else if self.drop_sites && (updated == 0 || matched.all_missing) {
                    return Ok(());
                }
                let convert = self.convert.as_mut().expect("-f has a format");
                convert.line(hdr_out, rec, smpl_pass, &mut self.line)?;
                if !self.line.is_empty() {
                    writer
                        .write_all(&self.line)
                        .map_err(|e| write_error(e, name))?;
                }
            }
            Out::Vcf { sink, .. } => {
                let spare = sink.spare().unwrap_or_else(|| source.take());
                // Without -d a record is written once, at the end of its processing.
                let done = if self.duplicate {
                    let mut copy = spare;
                    copy.copy_from(rec);
                    copy
                } else {
                    std::mem::replace(rec, spare)
                };
                if let Some(back) = sink.write(hdr_out, done)? {
                    source.recycle(back);
                }
            }
        }
        Ok(())
    }
}

/// What the transcripts of a record that go into one output line matched.
#[derive(Clone, Copy, Debug)]
struct Matched {
    /// Whether a transcript passed the `-s` severity selection.
    severity_pass: bool,
    /// Whether all the extracted annotations were empty.
    all_missing: bool,
}

impl Matched {
    /// Before any transcript is looked at.
    const NONE: Self = Self {
        severity_pass: false,
        all_missing: true,
    };
}

/// Appends a value to an annotation, comma-separated (bcftools' `annot_append`).
fn append(a: &mut Annot, value: &[u8]) {
    if !a.value.is_empty() {
        a.value.push(b',');
    }
    a.value.extend_from_slice(value);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_lists() {
        let mut f = Vec::new();
        let mut warned = false;
        parse_array_real("X", b"0.5,.,1e-3,", &mut f, &mut warned);
        assert_eq!(f.len(), 3);
        assert_eq!(f[0], 0.5);
        assert_eq!(f[1].to_bits(), FLOAT_MISSING);
        assert!(!warned);
        let mut i = Vec::new();
        parse_array_int32("X", b"3,abc,-2", &mut i, &mut warned);
        assert_eq!(i, [3, INT32_MISSING, -2]);
        assert!(warned);
    }

    #[test]
    fn transcripts_and_subfields() {
        let mut c = Csq {
            buf: b"A|missense|G1,C|intron".to_vec(),
            ..Csq::default()
        };
        assert!(!c.split(1));
        assert_eq!(c.trs.len(), 2);
        assert_eq!(c.field(0, 1), b"missense");
        assert_eq!(c.field(1, 1), b"intron");
        assert!(c.split(2));
    }
}
