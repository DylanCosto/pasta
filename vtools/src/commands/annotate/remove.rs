//! `-x/--remove`: annotations to remove from records and the header (bcftools'
//! `init_remove_annots` and the `remove_*` handlers).

use std::collections::HashMap;

use anyhow::bail;
use pasta::variant::{Header, HeaderRecordType, Record, record::PassOnEmpty};

/// One removal (bcftools' `rm_tag_t` with its handler).
#[derive(Debug)]
pub(super) enum Removal {
    Id,
    /// All filters (`None`), or one, by its input-header ID (`None` when undefined).
    Filter(Option<(Vec<u8>, Option<i32>)>),
    Qual,
    /// All INFO fields.
    Info,
    InfoTag(Vec<u8>),
    /// All FORMAT fields but GT.
    Format,
    FormatTag(Vec<u8>),
}

/// The removals and whether removing the last filter sets PASS (bcftools' `flt_keep_pass`).
#[derive(Debug)]
pub(super) struct Removals {
    pub(super) tags: Vec<Removal>,
    pub(super) keep_pass: PassOnEmpty,
}

fn starts_with_ci(s: &[u8], prefix: &[u8]) -> bool {
    s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix)
}

fn lossy(s: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(s)
}

fn defined(h: &Header, kind: HeaderRecordType, name: &[u8]) -> bool {
    h.id(name).is_some_and(|id| h.id_exists(kind, id))
}

/// The bit of a class in bcftools' keep set (`1 << BCF_HL_*`).
fn class_bit(kind: HeaderRecordType) -> u8 {
    match kind {
        HeaderRecordType::Filter => 1,
        HeaderRecordType::Info => 2,
        _ => 4,
    }
}

/// The `BCF_HL_*` number of a class.
fn class_num(kind: HeaderRecordType) -> u8 {
    match kind {
        HeaderRecordType::Filter => 0,
        HeaderRecordType::Info => 1,
        _ => 2,
    }
}

/// Removes all lines of a class from the output header, except FORMAT/GT (bcftools'
/// `remove_hdr_lines`): INFO and FORMAT tags become undefined, FILTER tags stay defined.
fn remove_hdr_lines(hdr_out: &mut Header, kind: HeaderRecordType) {
    let ids: Vec<Vec<u8>> = hdr_out
        .records()
        .filter(|r| r.kind() == kind)
        .filter_map(|r| r.id().map(<[u8]>::to_vec))
        .collect();
    for id in ids {
        if kind == HeaderRecordType::Format && id == b"GT" {
            continue;
        }
        hdr_out.remove(kind, Some(&id));
        if kind != HeaderRecordType::Filter {
            hdr_out.undefine_id(kind, &id);
        }
    }
}

/// What `init_remove_annots` needs besides the headers.
pub(super) struct RemoveSpec<'a> {
    pub(super) list: &'a str,
    pub(super) keep_sites: bool,
    /// Whether `--set-id` reads INFO tags (which removal would run before).
    pub(super) set_id_reads_info: bool,
    pub(super) set_id_fmt: &'a str,
}

/// bcftools' `init_remove_annots`: parses `-x` and removes the definitions from `hdr_out`
/// (unless `-k`).
#[expect(
    clippy::too_many_lines,
    reason = "a line-for-line port of init_remove_annots"
)]
pub(super) fn init_remove_annots(
    spec: &RemoveSpec<'_>,
    hdr: &Header,
    hdr_out: &mut Header,
) -> anyhow::Result<Removals> {
    let mut out = Removals {
        tags: Vec::new(),
        keep_pass: PassOnEmpty::Empty,
    };
    let mut keep: HashMap<Vec<u8>, u8> = HashMap::new();
    let (mut keep_info, mut keep_fmt, mut keep_flt) = (false, false, false);
    let keep_sites = spec.keep_sites;
    let mut items: Vec<&[u8]> = spec.list.as_bytes().split(|&c| c == b',').collect();
    // bcftools stops at the end of the string: a trailing comma adds no item.
    if spec.list.ends_with(',') {
        items.pop();
    }
    if spec.list.is_empty() {
        items.clear();
    }
    let info_conflict = |tag: &[u8]| -> anyhow::Result<()> {
        if spec.set_id_reads_info {
            bail!(
                "Error: `--remove INFO/{}` is executed first, cannot combine with `--set-id {}`",
                lossy(tag),
                spec.set_id_fmt
            );
        }
        Ok(())
    };
    for item in items {
        let mut kind = None;
        let mut s = item;
        for (prefix, k, is_keep) in [
            (&b"INFO/"[..], HeaderRecordType::Info, false),
            (b"INF/", HeaderRecordType::Info, false),
            (b"FORMAT/", HeaderRecordType::Format, false),
            (b"FMT/", HeaderRecordType::Format, false),
            (b"FILTER/", HeaderRecordType::Filter, false),
            (b"^INFO/", HeaderRecordType::Info, true),
            (b"^INF/", HeaderRecordType::Info, true),
            (b"^FORMAT/", HeaderRecordType::Format, true),
            (b"^FMT/", HeaderRecordType::Format, true),
            (b"^FILTER/", HeaderRecordType::Filter, true),
        ] {
            if starts_with_ci(item, prefix) {
                kind = Some(k);
                s = &item[prefix.len()..];
                if is_keep {
                    match k {
                        HeaderRecordType::Info => keep_info = true,
                        HeaderRecordType::Format => keep_fmt = true,
                        _ => keep_flt = true,
                    }
                }
                break;
            }
        }
        match kind {
            Some(HeaderRecordType::Filter) => {
                if !keep_flt {
                    out.keep_pass = PassOnEmpty::Pass;
                    let id = hdr.id(s);
                    if id.is_some_and(|id| hdr.id_exists(HeaderRecordType::Filter, id)) {
                        if !keep_sites {
                            hdr_out.remove(HeaderRecordType::Filter, Some(s));
                        }
                    } else if keep_sites {
                        bail!(
                            "Error: The filter \"{}\" is not defined in the header, cannot use the -k option",
                            lossy(s)
                        );
                    } else {
                        eprintln!(
                            "Warning: The filter \"{}\" is not defined in the header",
                            lossy(s)
                        );
                    }
                    out.tags
                        .push(Removal::Filter(Some((s.to_vec(), id.map(|i| i as i32)))));
                } else {
                    *keep.entry(s.to_vec()).or_default() |= class_bit(HeaderRecordType::Filter);
                }
            }
            Some(k) => {
                let is_keep = (k == HeaderRecordType::Format && keep_fmt)
                    || (k == HeaderRecordType::Info && keep_info);
                if !defined(hdr, k, s) {
                    if keep_sites {
                        bail!(
                            "Error: The tag \"{}\" is not defined in the header, cannot use the -k option",
                            lossy(s)
                        );
                    }
                    eprintln!(
                        "Warning: The tag \"{}\" not defined in the header",
                        lossy(s)
                    );
                    if k == HeaderRecordType::Info {
                        info_conflict(s)?;
                        out.tags.push(Removal::InfoTag(s.to_vec()));
                    } else {
                        out.tags.push(Removal::FormatTag(s.to_vec()));
                    }
                } else if is_keep {
                    *keep.entry(s.to_vec()).or_default() |= class_bit(k);
                } else {
                    if k == HeaderRecordType::Info {
                        info_conflict(s)?;
                        out.tags.push(Removal::InfoTag(s.to_vec()));
                    } else {
                        out.tags.push(Removal::FormatTag(s.to_vec()));
                    }
                    if !keep_sites {
                        hdr_out.remove(k, Some(s));
                    }
                }
            }
            None => {
                if s.eq_ignore_ascii_case(b"ID") {
                    out.tags.push(Removal::Id);
                } else if s.eq_ignore_ascii_case(b"FILTER") {
                    out.tags.push(Removal::Filter(None));
                    if !keep_sites {
                        remove_hdr_lines(hdr_out, HeaderRecordType::Filter);
                    }
                } else if s.eq_ignore_ascii_case(b"QUAL") {
                    out.tags.push(Removal::Qual);
                } else if s.eq_ignore_ascii_case(b"INFO") {
                    if spec.set_id_reads_info {
                        bail!(
                            "Error: `--remove INFO` is executed first, cannot combine with `--set-id {}`",
                            spec.set_id_fmt
                        );
                    }
                    out.tags.push(Removal::Info);
                    if !keep_sites {
                        remove_hdr_lines(hdr_out, HeaderRecordType::Info);
                    }
                } else if s.eq_ignore_ascii_case(b"FMT") || s.eq_ignore_ascii_case(b"FORMAT") {
                    out.tags.push(Removal::Format);
                    if !keep_sites {
                        remove_hdr_lines(hdr_out, HeaderRecordType::Format);
                    }
                } else if !s.is_empty() {
                    if defined(hdr, HeaderRecordType::Info, s) {
                        bail!("Error: did you mean INFO/{}?", lossy(s));
                    }
                    if defined(hdr, HeaderRecordType::Format, s) {
                        bail!("Error: did you mean FORMAT/{}?", lossy(s));
                    }
                    if !keep_sites {
                        match s.strip_prefix(b"##") {
                            Some(key) => hdr_out.remove(HeaderRecordType::Generic, Some(key)),
                            None => hdr_out.remove(HeaderRecordType::Structured, Some(s)),
                        }
                    }
                }
            }
        }
    }
    if keep_flt || keep_info || keep_fmt {
        let lines: Vec<(HeaderRecordType, Vec<u8>)> = hdr
            .records()
            .filter(|r| {
                matches!(
                    r.kind(),
                    HeaderRecordType::Filter | HeaderRecordType::Info | HeaderRecordType::Format
                )
            })
            .filter_map(|r| Some((r.kind(), r.id()?.to_vec())))
            .collect();
        for (kind, id) in lines {
            let wanted = match kind {
                HeaderRecordType::Filter => keep_flt,
                HeaderRecordType::Info => keep_info,
                _ => keep_fmt,
            };
            if !wanted {
                continue;
            }
            // bcftools tests `value >> type`: a tag kept as INFO also keeps a FILTER of the
            // same name, and a FORMAT tag both.
            if keep.get(&id).is_some_and(|&v| v >> class_num(kind) != 0) {
                if kind == HeaderRecordType::Filter && id == b"PASS" {
                    out.keep_pass = PassOnEmpty::Pass;
                }
                continue;
            }
            if !keep_sites {
                hdr_out.remove(kind, Some(&id));
            }
            out.tags.push(match kind {
                HeaderRecordType::Info => Removal::InfoTag(id),
                HeaderRecordType::Format => Removal::FormatTag(id),
                _ => {
                    let hid = hdr.id(&id).map(|i| i as i32);
                    Removal::Filter(Some((id, hid)))
                }
            });
        }
    }
    if out.tags.is_empty() {
        bail!("No matching tag in -x {}", spec.list);
    }
    Ok(out)
}

impl Removals {
    /// Applies the removals to a record (the `rm[i].handler` calls of `annotate_line`).
    pub(super) fn apply(&self, hdr: &Header, rec: &mut Record) -> anyhow::Result<()> {
        for tag in &self.tags {
            match tag {
                Removal::Id => rec.update_id(None),
                Removal::Filter(None) => {
                    // bcftools passes the PASS flag as the filter count; with no IDs that only
                    // empties FILTER.
                    rec.update_filter(&[]);
                }
                Removal::Filter(Some((_, None))) => bail!(
                    "Error: Cannot proceed, not even with the --force option, bad things could happen.\n       Note that \"vtools annotate -x FILTER\" can be used to remove ALL filters.\n       Even better, use \"vtools view -h\" and \"vtools reheader\" to fix the header!"
                ),
                Removal::Filter(Some((_, Some(id)))) => rec.remove_filter(*id, self.keep_pass),
                Removal::Qual => rec.set_qual(None),
                Removal::Info => rec.remove_all_info(hdr),
                Removal::InfoTag(key) => {
                    let _ = rec.remove_info(hdr, key);
                }
                Removal::Format => rec.remove_all_format(hdr, hdr.id(b"GT")),
                Removal::FormatTag(key) => {
                    let _ = rec.remove_format(hdr, key);
                }
            }
        }
        Ok(())
    }
}
