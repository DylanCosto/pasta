//! The output header (bcftools merge's `merge_headers` and `hdr_add_localized_tags`, and
//! HTSlib's `bcf_hdr_set` for `--use-header`).

use anyhow::bail;
use pasta::variant::{Header, HeaderRecordType, Number, readlist};

/// What happens to a sample name that is already in the output header.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum DuplicateSamples {
    /// It is an error.
    Fail,
    /// It gets a prefix (`--force-samples`).
    Rename,
}

/// Adds a file's header lines and samples to the output header (bcftools' `merge_headers`).
/// A sample name already present gets `clash_prefix:` prepended, as often as needed, with
/// [`DuplicateSamples::Rename`]; otherwise it is an error.
pub(super) fn merge_headers(
    hw: &mut Header,
    hr: &Header,
    clash_prefix: &str,
    duplicates: DuplicateSamples,
) -> anyhow::Result<()> {
    hw.merge(hr)?;
    for name in hr.samples() {
        let mut name: Vec<u8> = name.to_vec();
        while hw.sample_id(&name).is_some() {
            if duplicates == DuplicateSamples::Fail {
                bail!(
                    "Error: Duplicate sample names ({}), use --force-samples to proceed anyway.",
                    String::from_utf8_lossy(&name)
                );
            }
            let mut renamed = clash_prefix.as_bytes().to_vec();
            renamed.push(b':');
            renamed.extend_from_slice(&name);
            name = renamed;
        }
        hw.add_sample(&name)?;
    }
    Ok(())
}

/// Adds FORMAT/LAA and an `L<tag>` definition for each Number=A/G/R FORMAT tag, for `-L`
/// (bcftools' `hdr_add_localized_tags`).
pub(super) fn add_localized_tags(hdr: &mut Header) -> anyhow::Result<()> {
    let mut lines: Vec<Vec<u8>> = Vec::new();
    for rec in hdr.records() {
        if rec.kind() != HeaderRecordType::Format {
            continue;
        }
        let Some(id) = rec.get(b"ID") else {
            continue;
        };
        let Some(dict_id) = hdr.id(id) else {
            continue;
        };
        let localize = matches!(
            hdr.id_number(HeaderRecordType::Format, dict_id),
            Some(Number::G | Number::A | Number::R)
        );
        if !localize {
            continue;
        }
        let mut line = b"##".to_vec();
        line.extend_from_slice(rec.key());
        line.extend_from_slice(b"=<");
        let mut nout = 0;
        for (k, v) in rec.fields() {
            if k == b"IDX" {
                continue;
            }
            if nout > 0 {
                line.push(b',');
            }
            if k == b"ID" {
                line.extend_from_slice(b"ID=L");
                line.extend_from_slice(v);
            } else if k == b"Number" {
                line.extend_from_slice(b"Number=.");
            } else if k == b"Description" && v.first() == Some(&b'"') {
                line.extend_from_slice(b"Description=\"Localized field: ");
                line.extend_from_slice(&v[1..]);
            } else {
                line.extend_from_slice(k);
                line.push(b'=');
                line.extend_from_slice(v);
            }
            nout += 1;
        }
        line.push(b'>');
        lines.push(line);
    }
    if lines.is_empty() {
        return Ok(());
    }
    hdr.append(
        b"##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"Localized alleles: subset of alternate alleles relevant for each sample\">",
    )?;
    for l in lines {
        hdr.append(&l)?;
    }
    Ok(())
}

/// Adds the header lines and samples of a text file to `hdr` (HTSlib's `bcf_hdr_set`): every
/// line but the last must be a `##` line, and the last the `#CHROM` line.
pub(super) fn read_header_file(hdr: &mut Header, path: &str) -> anyhow::Result<()> {
    let lines = readlist::read_lines_or_inline(path.as_bytes())?;
    let Some((sample_line, meta)) = lines.split_last() else {
        bail!("no header lines");
    };
    for l in meta {
        if !l.starts_with(b"##") {
            bail!("not a header line");
        }
        hdr.append(l)?;
    }
    const MANDATORY: &[u8] = b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let Some(rest) = sample_line.strip_prefix(MANDATORY) else {
        bail!("invalid #CHROM line");
    };
    if rest.is_empty() {
        return Ok(());
    }
    let Some(samples) = rest.strip_prefix(b"\tFORMAT\t") else {
        bail!("invalid #CHROM line");
    };
    for s in samples.split(|&c| c == b'\t') {
        hdr.add_sample(s)?;
    }
    Ok(())
}
