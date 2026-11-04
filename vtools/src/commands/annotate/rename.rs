//! `--rename-chrs` and `--rename-annots` (and `-c NEW:=OLD` without `-a`): renaming
//! sequences and tags (bcftools' `rename_chrs`, `rename_annots_init1`, and `rename_annots`).

use anyhow::bail;
use pasta::variant::{
    Header, HeaderRecordType, Record,
    readlist::{self, ListSource},
};

use super::columns::format_line;

fn starts_with_ci(s: &[u8], prefix: &[u8]) -> bool {
    s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}

fn lossy(s: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(s)
}

/// Splits `old<whitespace>new[<whitespace>...]` into its first two words.
fn two_words(line: &[u8]) -> Option<(&[u8], &[u8])> {
    let end = line.iter().position(u8::is_ascii_whitespace)?;
    let rest = line[end..].trim_ascii_start();
    let new_end = rest
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(rest.len());
    Some((&line[..end], &rest[..new_end]))
}

/// bcftools' `rename_chrs`: renames the contigs listed in `fname` (`old new` per line) in the
/// output header. Names are looked up as the input header has them.
pub(super) fn rename_chrs(fname: &str, hdr: &Header, hdr_out: &mut Header) -> anyhow::Result<()> {
    let Ok(lines) = readlist::read_list(fname.as_bytes(), ListSource::File) else {
        bail!("Could not read: {fname}");
    };
    for line in &lines {
        if !line.iter().any(u8::is_ascii_whitespace) {
            bail!("Could not parse: {fname}");
        }
        let (old, new) = two_words(line).expect("the line has whitespace");
        let Some(rid) = hdr.contig_id(old) else {
            continue;
        };
        hdr_out.rename_contig(rid, new);
    }
    Ok(())
}

/// One tag rename (bcftools' `rename_tag_t`).
#[derive(Debug)]
pub(super) struct Rename {
    kind: HeaderRecordType,
    old_id: i32,
    new_id: i32,
}

/// Reads `--rename-annots FILE` lines (`TYPE/old new`) after the `-c` renames.
pub(super) fn read_rename_file(
    fname: &str,
    renames: &mut Vec<(Vec<u8>, Vec<u8>)>,
) -> anyhow::Result<()> {
    let Ok(lines) = readlist::read_list(fname.as_bytes(), ListSource::File) else {
        bail!("Could not read: {fname}");
    };
    for line in &lines {
        match two_words(line) {
            Some((old, new)) if !new.is_empty() => renames.push((old.to_vec(), new.to_vec())),
            _ => bail!("Could not parse: {}", lossy(line)),
        }
    }
    Ok(())
}

/// bcftools' `rename_annots_init1`: renames one tag in `hdr_out`, in place when it can.
fn init_one(
    old_tag: &[u8],
    new_tag: &[u8],
    hdr_out: &mut Header,
) -> anyhow::Result<Option<Rename>> {
    let (kind, ori) = if starts_with_ci(old_tag, b"info/") {
        (HeaderRecordType::Info, &old_tag[5..])
    } else if starts_with_ci(old_tag, b"format/") {
        (HeaderRecordType::Format, &old_tag[7..])
    } else if starts_with_ci(old_tag, b"fmt/") {
        (HeaderRecordType::Format, &old_tag[4..])
    } else if starts_with_ci(old_tag, b"filter/") {
        (HeaderRecordType::Filter, &old_tag[7..])
    } else {
        bail!(
            "Cannot rename \"{}\" to \"{}\"",
            lossy(old_tag),
            lossy(new_tag)
        );
    };
    let new = if starts_with_ci(new_tag, b"info/") {
        if kind != HeaderRecordType::Info {
            bail!("Cannot transfer {} to INFO", lossy(ori));
        }
        &new_tag[5..]
    } else if starts_with_ci(new_tag, b"format/") || starts_with_ci(new_tag, b"fmt/") {
        if kind != HeaderRecordType::Format {
            bail!("Cannot transfer {} to FORMAT", lossy(ori));
        }
        &new_tag[if starts_with_ci(new_tag, b"fmt/") {
            4
        } else {
            7
        }..]
    } else if starts_with_ci(new_tag, b"filter/") {
        if kind != HeaderRecordType::Filter {
            bail!("Cannot transfer {} to FILTER", lossy(ori));
        }
        &new_tag[7..]
    } else {
        new_tag
    };
    let fail = || anyhow::anyhow!("Cannot rename \"{}\" to \"{}\"", lossy(ori), lossy(new));
    let old_id = hdr_out.id(ori).ok_or_else(fail)?;
    let old_line = match hdr_out.definition(kind, ori) {
        Some(r) => format_line(r),
        None => return Err(fail()),
    };
    let shared = [
        HeaderRecordType::Info,
        HeaderRecordType::Format,
        HeaderRecordType::Filter,
    ]
    .into_iter()
    .any(|k| k != kind && hdr_out.id_exists(k, old_id));
    let new_existing = hdr_out.id(new);
    if !shared && new_existing.is_none_or(|i| i == 0) {
        if !hdr_out.rename_id(old_id, new) {
            bail!("The tag \"{}\" already exists in the header", lossy(new));
        }
        return Ok(None);
    }
    if hdr_out.definition(kind, new).is_some() {
        let class = match kind {
            HeaderRecordType::Info => "INFO",
            HeaderRecordType::Format => "FORMAT",
            _ => "FILTER",
        };
        bail!(
            "The tag \"{class}/{}\" already exists in the header",
            lossy(new)
        );
    }
    // Replace ID=old with ID=new in the formatted line.
    let Some(at) = old_line.windows(3).position(|w| w == b"ID=") else {
        bail!(
            "Malformed header record for {}: {}",
            lossy(ori),
            lossy(&old_line)
        );
    };
    let idp = at + 3;
    let mut idend = idp;
    while idend < old_line.len() && old_line[idend] != b',' && old_line[idend] != b'>' {
        idend += 1;
    }
    let mut line = old_line[..idp].to_vec();
    line.extend_from_slice(new);
    line.extend_from_slice(&old_line[idend..]);
    if hdr_out.append(&line).is_err() {
        bail!("Failed to append renamed header line: {}", lossy(&line));
    }
    let Some(new_id) = hdr_out.id(new) else {
        bail!("Failed to obtain new header ID for {}", lossy(new));
    };
    hdr_out.remove(kind, Some(ori));
    Ok(Some(Rename {
        kind,
        old_id: old_id as i32,
        new_id: new_id as i32,
    }))
}

/// bcftools' `rename_annots_init`: applies the renames to `hdr_out` and returns those that
/// records need (when the numeric ID changes).
pub(super) fn init_renames(
    renames: &[(Vec<u8>, Vec<u8>)],
    hdr_out: &mut Header,
) -> anyhow::Result<Vec<Rename>> {
    let mut out = Vec::new();
    for (old, new) in renames {
        if let Some(r) = init_one(old, new, hdr_out)?
            && r.old_id != r.new_id
        {
            out.push(r);
        }
    }
    Ok(out)
}

/// bcftools' `rename_annots`: moves a record's renamed tags to their new IDs.
pub(super) fn apply(renames: &[Rename], rec: &mut Record) {
    for r in renames {
        match r.kind {
            HeaderRecordType::Info => rec.rename_info_id(r.old_id, r.new_id),
            HeaderRecordType::Format => rec.rename_format_id(r.old_id, r.new_id),
            _ => rec.rename_filter_id(r.old_id, r.new_id),
        }
    }
}
