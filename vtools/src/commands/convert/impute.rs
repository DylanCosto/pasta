//! VCF to IMPUTE2 gen/sample and hap/legend/sample, and to SHAPEIT hap/sample (bcftools'
//! `vcf_to_gensample`, `vcf_to_haplegendsample`, and `vcf_to_hapsample`).

use anyhow::bail;
use pasta::variant::{
    Header, Record,
    number::push_int,
    readlist::{self, ListSource},
};

use super::{
    Input, Settings, open_vcf, split_names,
    text::{TextCompression, TextWriter, has_gz_suffix},
};
use crate::raw_input::strerror;

/// Reads `--sex`: one `SAMPLE M|F` per line, giving `1` or `2` per header sample (bcftools'
/// `init_sample2sex`).
fn init_sample2sex(header: &Header, sex_fname: &str) -> anyhow::Result<Vec<u8>> {
    let mut sex = vec![0u8; header.sample_count()];
    let Ok(lines) = readlist::read_list(sex_fname.as_bytes(), ListSource::File) else {
        bail!("Could not read {sex_fname}");
    };
    for line in &lines {
        let end = line
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(line.len());
        let Some(id) = header.sample_id(&line[..end]) else {
            continue;
        };
        let rest = &line[end..];
        let first = rest
            .iter()
            .copied()
            .find(|c| !c.is_ascii_whitespace())
            .unwrap_or(0);
        sex[id as usize] = match first {
            b'M' => b'1',
            b'F' => b'2',
            _ => bail!(
                "Could not parse {sex_fname}: {}",
                String::from_utf8_lossy(line)
            ),
        };
    }
    if let Some(i) = sex.iter().position(|&c| c == 0) {
        bail!(
            "Missing sex for sample {} in {sex_fname}",
            String::from_utf8_lossy(&header.samples()[i])
        );
    }
    Ok(sex)
}

/// The names given as `-` or `.`-suppressed list items: `None` for an empty or `.` item.
fn named(item: &str) -> Option<String> {
    (!item.is_empty() && item != ".").then(|| item.to_owned())
}

/// Writes the `.samples` file of gen/sample and hap/sample outputs.
fn write_oxford_samples(s: &Settings, header: &Header, fname: &str) -> anyhow::Result<()> {
    let sex = match &s.sex {
        Some(f) => Some(init_sample2sex(header, f)?),
        None => None,
    };
    let compression = if fname.len() > 3 && has_gz_suffix(fname) {
        TextCompression::Gzip
    } else {
        TextCompression::Plain
    };
    let mut out = TextWriter::create(fname, compression)
        .map_err(|e| anyhow::anyhow!("Error writing {fname}: {}", strerror(&e)))?;
    out.write(if sex.is_some() {
        b"ID_1 ID_2 missing sex\n0 0 0 0\n"
    } else {
        b"ID_1 ID_2 missing\n0 0 0\n"
    })?;
    let mut line = Vec::new();
    for (i, name) in header.samples().iter().enumerate() {
        line.clear();
        line.extend_from_slice(name);
        line.push(b' ');
        line.extend_from_slice(name);
        line.extend_from_slice(b" 0");
        if let Some(sex) = &sex {
            line.push(b' ');
            line.push(sex[i]);
        }
        line.push(b'\n');
        out.write(&line)?;
    }
    out.close()
}

/// How a data output (gen, hap, legend) is compressed: gzip unless the name lacks `.gz`.
fn data_compression(fname: &str) -> TextCompression {
    if has_gz_suffix(fname) {
        TextCompression::Gzip
    } else {
        TextCompression::Plain
    }
}

/// Why a record was skipped.
#[derive(Default)]
struct Skipped {
    no_alt: u64,
    non_biallelic: u64,
    filtered: u64,
}

impl Skipped {
    /// Applies the filter and the biallelic requirement; `true` when the record is kept.
    fn keep(&mut self, input: &mut Input, rec: &mut Record) -> anyhow::Result<bool> {
        if !input.pass(rec)? {
            self.filtered += 1;
            return Ok(false);
        }
        // An ALT allele is required.
        if rec.allele_count() < 2 {
            self.no_alt += 1;
            return Ok(false);
        }
        // Biallelic sites only.
        if rec.allele_count() > 2 {
            if self.non_biallelic == 0 {
                eprintln!(
                    "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'vtools norm -m-'."
                );
            }
            self.non_biallelic += 1;
            return Ok(false);
        }
        Ok(true)
    }

    fn total(&self) -> u64 {
        self.no_alt + self.non_biallelic + self.filtered
    }
}

/// Reads the next record; a read error is fatal.
fn next(input: &mut Input, rec: &mut Record) -> anyhow::Result<bool> {
    match input.reader.next_record(rec) {
        Ok(b) => Ok(b),
        Err(e) => bail!("Error: {e}"),
    }
}

/// Formats a record through the input's format string.
fn format_line(input: &mut Input, rec: &mut Record, out: &mut Vec<u8>) -> anyhow::Result<()> {
    let convert = input
        .convert
        .as_mut()
        .expect("the conversion opened the input with a format");
    convert.line(input.reader.header(), rec, None, out)
}

/// VCF to IMPUTE2 gen/sample (bcftools' `vcf_to_gensample`).
pub(super) fn vcf_to_gensample(s: &Settings) -> anyhow::Result<()> {
    let mut fmt = String::new();
    if s.gen_3n6 {
        fmt.push_str("%CHROM ");
    }
    fmt.push_str("%CHROM:%POS\\_%REF\\_%FIRST_ALT ");
    if s.vcf_ids {
        fmt.push_str("%ID ");
    } else {
        fmt.push_str("%CHROM:%POS\\_%REF\\_%FIRST_ALT ");
    }
    fmt.push_str("%POS %REF %FIRST_ALT");
    match s.tag.as_deref() {
        None | Some("GT") => fmt.push_str("%_GT_TO_PROB3"),
        Some("PL") => fmt.push_str("%_PL_TO_PROB3"),
        Some("GP") => fmt.push_str("%_GP_TO_PROB3"),
        Some(t) => bail!("todo: --tag {t}"),
    }
    fmt.push('\n');
    let mut input = open_vcf(s, Some(&fmt))?;

    let files = split_names(&s.outfname);
    let (gen_fname, sample_fname) = match files.len() {
        1 => (
            Some(format!("{}.gen.gz", s.outfname)),
            Some(format!("{}.samples", s.outfname)),
        ),
        2 => (named(&files[0]), named(&files[1])),
        _ => bail!("Error parsing --gensample filenames: {}", s.outfname),
    };
    if let Some(f) = &gen_fname {
        eprintln!("Gen file: {f}");
    }
    if let Some(f) = &sample_fname {
        eprintln!("Sample file: {f}");
    }
    if let Some(f) = &sample_fname {
        write_oxford_samples(s, input.reader.header(), f)?;
    }
    let Some(gen_fname) = gen_fname else {
        return Ok(());
    };

    let mut gout = TextWriter::create(&gen_fname, data_compression(&gen_fname))
        .map_err(|e| anyhow::anyhow!("Error writing {gen_fname}: {}", strerror(&e)))?;
    let mut skipped = Skipped::default();
    let (mut ndup, mut nok) = (0u64, 0u64);
    let mut prev: Option<(i32, i64)> = None;
    let mut rec = Record::new();
    let mut str = Vec::new();
    while next(&mut input, &mut rec)? {
        if !skipped.keep(&mut input, &mut rec)? {
            continue;
        }
        // Skip duplicate positions, or else SHAPEIT complains.
        let here = (rec.rid(), rec.pos());
        if !s.keep_duplicates && prev == Some(here) {
            ndup += 1;
            continue;
        }
        prev = Some(here);
        format_line(&mut input, &mut rec, &mut str)?;
        if !str.is_empty() {
            gout.write(&str)?;
            nok += 1;
        }
    }
    if let Some(e) = input.reader.take_pending_error() {
        bail!("Error: {e}");
    }
    eprintln!(
        "{nok} records written, {} skipped: {}/{}/{}/{ndup} no-ALT/non-biallelic/filtered/duplicated",
        skipped.total() + ndup,
        skipped.no_alt,
        skipped.non_biallelic,
        skipped.filtered
    );
    gout.close()
}

/// VCF to IMPUTE2 hap/legend/sample (bcftools' `vcf_to_haplegendsample`).
pub(super) fn vcf_to_haplegendsample(s: &Settings) -> anyhow::Result<()> {
    let fmt = if s.hap2dip {
        "%_GT_TO_HAP2\n"
    } else {
        "%_GT_TO_HAP\n"
    };
    let mut input = open_vcf(s, Some(fmt))?;

    let files = split_names(&s.outfname);
    let (hap_fname, legend_fname, sample_fname) = match files.len() {
        1 => (
            Some(format!("{}.hap.gz", s.outfname)),
            Some(format!("{}.legend.gz", s.outfname)),
            Some(format!("{}.samples", s.outfname)),
        ),
        3 => (named(&files[0]), named(&files[1]), named(&files[2])),
        _ => bail!("Error parsing --hapslegendsample filenames: {}", s.outfname),
    };
    if let Some(f) = &hap_fname {
        eprintln!("Hap file: {f}");
    }
    if let Some(f) = &legend_fname {
        eprintln!("Legend file: {f}");
    }
    if let Some(f) = &sample_fname {
        eprintln!("Sample file: {f}");
    }
    if let Some(f) = &sample_fname {
        let header = input.reader.header();
        let sex = match &s.sex {
            Some(sf) => Some(init_sample2sex(header, sf)?),
            None => None,
        };
        let compression = if f.len() > 3 && has_gz_suffix(f) {
            TextCompression::Gzip
        } else {
            TextCompression::Plain
        };
        let mut out = TextWriter::create(f, compression)
            .map_err(|e| anyhow::anyhow!("Error writing {f}: {}", strerror(&e)))?;
        out.write(b"sample population group sex\n")?;
        let mut line = Vec::new();
        for (i, name) in header.samples().iter().enumerate() {
            line.clear();
            for _ in 0..3 {
                line.extend_from_slice(name);
                line.push(b' ');
            }
            line.push(sex.as_ref().map_or(b'2', |x| x[i]));
            line.push(b'\n');
            out.write(&line)?;
        }
        out.close()?;
    }
    if hap_fname.is_none() && legend_fname.is_none() {
        return Ok(());
    }

    let mut hout = match &hap_fname {
        Some(f) => Some(
            TextWriter::create_threaded(f, data_compression(f), s.threads)
                .map_err(|e| anyhow::anyhow!("Error writing {f}: {}", strerror(&e)))?,
        ),
        None => None,
    };
    let mut lout = match &legend_fname {
        Some(f) => {
            let mut w = TextWriter::create(f, data_compression(f))
                .map_err(|e| anyhow::anyhow!("Error writing {f}: {}", strerror(&e)))?;
            w.write(b"id position a0 a1\n")?;
            Some(w)
        }
        None => None,
    };

    let mut skipped = Skipped::default();
    let mut nok = 0u64;
    let mut rec = Record::new();
    let mut str = Vec::new();
    let mut legend = Vec::new();
    while next(&mut input, &mut rec)? {
        if !skipped.keep(&mut input, &mut rec)? {
            continue;
        }
        format_line(&mut input, &mut rec, &mut str)?;
        if str.is_empty() {
            continue;
        }
        if let Some(h) = &mut hout {
            h.write(&str)?;
        }
        if let Some(l) = &mut lout {
            legend.clear();
            let id = rec.id();
            if s.vcf_ids && id != b"." {
                legend.extend_from_slice(id);
            } else {
                legend.extend_from_slice(rec.chrom_safe(input.reader.header()));
                legend.push(b':');
                push_int(&mut legend, rec.pos() + 1);
                legend.push(b'_');
                legend.extend_from_slice(rec.allele(0));
                legend.push(b'_');
                legend.extend_from_slice(rec.allele(1));
            }
            legend.push(b' ');
            push_int(&mut legend, rec.pos() + 1);
            legend.push(b' ');
            legend.extend_from_slice(rec.allele(0));
            legend.push(b' ');
            legend.extend_from_slice(rec.allele(1));
            legend.push(b'\n');
            l.write(&legend)?;
        }
        nok += 1;
    }
    if let Some(e) = input.reader.take_pending_error() {
        bail!("Error: {e}");
    }
    eprintln!(
        "{nok} records written, {} skipped: {}/{}/{} no-ALT/non-biallelic/filtered",
        skipped.total(),
        skipped.no_alt,
        skipped.non_biallelic,
        skipped.filtered
    );
    if let Some(h) = hout {
        h.close()?;
    }
    if let Some(l) = lout {
        l.close()?;
    }
    Ok(())
}

/// VCF to SHAPEIT hap/sample (bcftools' `vcf_to_hapsample`).
pub(super) fn vcf_to_hapsample(s: &Settings) -> anyhow::Result<()> {
    let mut fmt = String::from(if s.vcf_ids {
        "%CHROM:%POS\\_%REF\\_%FIRST_ALT %ID %POS %REF %FIRST_ALT "
    } else {
        "%CHROM %CHROM:%POS\\_%REF\\_%FIRST_ALT %POS %REF %FIRST_ALT "
    });
    fmt.push_str(if s.hap2dip {
        "%_GT_TO_HAP2\n"
    } else {
        "%_GT_TO_HAP\n"
    });
    let mut input = open_vcf(s, Some(&fmt))?;

    let files = split_names(&s.outfname);
    let (hap_fname, sample_fname) = match files.len() {
        1 => (
            Some(format!("{}.hap.gz", s.outfname)),
            Some(format!("{}.samples", s.outfname)),
        ),
        2 => (named(&files[0]), named(&files[1])),
        _ => bail!("Error parsing --hapsample filenames: {}", s.outfname),
    };
    if let Some(f) = &hap_fname {
        eprintln!("Hap file: {f}");
    }
    if let Some(f) = &sample_fname {
        eprintln!("Sample file: {f}");
    }
    if let Some(f) = &sample_fname {
        write_oxford_samples(s, input.reader.header(), f)?;
    }
    let Some(hap_fname) = hap_fname else {
        return Ok(());
    };

    let mut hout = TextWriter::create_threaded(&hap_fname, data_compression(&hap_fname), s.threads)
        .map_err(|e| anyhow::anyhow!("Error writing {hap_fname}: {}", strerror(&e)))?;
    let mut skipped = Skipped::default();
    let mut nok = 0u64;
    let mut rec = Record::new();
    let mut str = Vec::new();
    while next(&mut input, &mut rec)? {
        if !skipped.keep(&mut input, &mut rec)? {
            continue;
        }
        format_line(&mut input, &mut rec, &mut str)?;
        if str.is_empty() {
            continue;
        }
        hout.write(&str)?;
        nok += 1;
    }
    if let Some(e) = input.reader.take_pending_error() {
        bail!("Error: {e}");
    }
    eprintln!(
        "{nok} records written, {} skipped: {}/{}/{} no-ALT/non-biallelic/filtered",
        skipped.total(),
        skipped.no_alt,
        skipped.non_biallelic,
        skipped.filtered
    );
    hout.close()
}
