//! Which samples receive annotations, and from which source samples (`-s`/`-S`; bcftools'
//! `init_sample_map` with `smpl_ilist_init`).

use anyhow::bail;
use pasta::variant::{
    Header,
    readlist::{self, ListSource},
};

/// The samples to annotate.
#[derive(Debug, Default)]
pub(super) struct SampleMap {
    /// Destination sample -> source sample (-1 for none), when the samples do not correspond
    /// one to one (bcftools' `sample_map`).
    pub(super) map: Option<Vec<i32>>,
    /// The number of annotated samples (bcftools' `nsmpl_annot`).
    pub(super) nsmpl_annot: usize,
    /// Whether a destination sample is in the source at all (`init_sample_map` returning -1
    /// otherwise, which fails a FORMAT transfer).
    pub(super) any_match: bool,
}

/// A parsed sample list (bcftools' `smpl_ilist_t` with `idx` and `pair`).
struct SampleList {
    idx: Vec<usize>,
    pair: Vec<Option<Vec<u8>>>,
}

/// Whether the character at `i` is whitespace not escaped by an odd number of backslashes.
fn is_space_or_escaped(s: &[u8], i: usize) -> bool {
    if !s[i].is_ascii_whitespace() {
        return false;
    }
    s[..i].iter().rev().take_while(|&&c| c == b'\\').count() % 2 == 0
}

/// How the sample list is read: in list order (`SMPL_REORDER`, tab-delimited annotations)
/// or in header order with `src dst` pairs (`SMPL_PAIR2`, VCF annotations).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ListKind {
    Tab,
    Vcf,
}

/// bcftools' `smpl_ilist_init(dst, list, is_file, SMPL_STRICT|SMPL_SINGLE|...)`.
fn sample_list(
    dst: &Header,
    list: &str,
    source: ListSource,
    kind: ListKind,
) -> anyhow::Result<SampleList> {
    let negate = list.starts_with('^');
    let body = if negate { &list[1..] } else { list };
    let Ok(names) = readlist::read_list(body.as_bytes(), source) else {
        bail!("Could not parse {list}");
    };
    let reorder = kind == ListKind::Tab && !negate;
    let nsmpl = dst.sample_count();
    let mut listed = Vec::new();
    let mut tmp = vec![false; nsmpl];
    let mut pair: Vec<Option<Vec<u8>>> = vec![None; nsmpl];
    for name in &names {
        let split = (0..name.len()).find(|&i| is_space_or_escaped(name, i));
        let (smpl1, smpl2) = match split {
            Some(i) => (&name[..i], Some(&name[i + 1..])),
            None => (name.as_slice(), None),
        };
        let lookup = match (kind, smpl2) {
            (ListKind::Vcf, Some(s2)) => s2,
            _ => smpl1,
        };
        let Some(idx) = dst.sample_id(lookup) else {
            bail!("No such sample: \"{}\"", String::from_utf8_lossy(lookup));
        };
        let idx = idx as usize;
        if reorder {
            listed.push(idx);
        } else {
            tmp[idx] = true;
            if smpl2.is_some() && kind == ListKind::Vcf {
                pair[idx] = Some(smpl1.to_vec());
            }
        }
    }
    if reorder {
        let n = listed.len();
        return Ok(SampleList {
            idx: listed,
            pair: vec![None; n],
        });
    }
    let mut out = SampleList {
        idx: Vec::new(),
        pair: Vec::new(),
    };
    for (i, &t) in tmp.iter().enumerate() {
        if t != negate {
            out.idx.push(i);
            out.pair.push(if negate { None } else { pair[i].take() });
        }
    }
    Ok(out)
}

/// bcftools' `init_sample_map`: `src` is the annotation VCF's header (`None` for a
/// tab-delimited file), `dst` the input's.
pub(super) fn init_sample_map(
    src: Option<&Header>,
    dst: &Header,
    samples: Option<(&str, ListSource)>,
) -> anyhow::Result<SampleMap> {
    let nsmpl = dst.sample_count();
    let Some((list, source)) = samples else {
        let mut sm = SampleMap {
            map: None,
            nsmpl_annot: nsmpl,
            any_match: true,
        };
        let Some(src) = src else {
            return Ok(sm);
        };
        sm.any_match = src.samples().iter().any(|s| dst.sample_id(s).is_some());
        if sm.any_match {
            sm.map = Some(
                dst.samples()
                    .iter()
                    .map(|s| src.sample_id(s).map_or(-1, |i| i as i32))
                    .collect(),
            );
        }
        return Ok(sm);
    };
    let kind = if src.is_none() {
        ListKind::Tab
    } else {
        ListKind::Vcf
    };
    let ilist = sample_list(dst, list, source, kind)?;
    if ilist.idx.is_empty() {
        bail!("Could not parse the samples: {list}");
    }
    let mut map = vec![-1i32; nsmpl];
    for (i, &idst) in ilist.idx.iter().enumerate() {
        let src_name: &[u8] = match &ilist.pair[i] {
            Some(p) => p,
            None => &dst.samples()[idst],
        };
        let mut isrc = i as i32;
        if let Some(src) = src {
            let Some(id) = src.sample_id(src_name) else {
                bail!(
                    "Sample \"{}\" not found in the annotation file",
                    String::from_utf8_lossy(src_name)
                );
            };
            isrc = id as i32;
        }
        map[idst] = isrc;
    }
    Ok(SampleMap {
        map: Some(map),
        nsmpl_annot: ilist.idx.len(),
        any_match: true,
    })
}
