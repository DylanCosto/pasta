//! Removing alleles from records (HTSlib's `bcf_remove_allele_set`, `bcf_remove_alleles`, and
//! `bcf_trim_alleles`), with the INFO and FORMAT fields that depend on them.

use crate::variant::{
    Header, HeaderRecordType, Record, Unpack,
    encoding::{
        FLOAT_MISSING, FLOAT_VECTOR_END, INT8_MISSING, INT8_VECTOR_END, INT16_MISSING,
        INT16_VECTOR_END, INT32_MISSING, INT32_VECTOR_END, enc_int1, enc_size, le_i16, le_i32,
        le_i64, le_u32,
    },
    log,
    record::{DIRTY_INF, FieldData, GetError, InfoField},
};

/// An error removing alleles.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RemoveAllelesError(pub String);

impl std::fmt::Display for RemoveAllelesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for RemoveAllelesError {}

/// The HTSlib function the messages are attributed to.
const FUNC: &str = "bcf_remove_allele_set";

// `BCF_VL_*` length codes.
const VL_VAR: u32 = 1;
const VL_A: u32 = 2;
const VL_G: u32 = 3;
const VL_R: u32 = 4;
const VL_LA: u32 = 6;
const VL_LG: u32 = 7;
const VL_LR: u32 = 8;

const VL_A_G_R: u32 = 1 << VL_A | 1 << VL_G | 1 << VL_R;
const VL_LA_LG_LR: u32 = 1 << VL_LA | 1 << VL_LG | 1 << VL_LR;
const VL_A_G_R_LA_LG_LR: u32 = VL_A_G_R | VL_LA_LG_LR;
const VL_A_R_LA_LR: u32 = 1 << VL_A | 1 << VL_R | 1 << VL_LA | 1 << VL_LR;

// `BCF_HT_*` value types.
const HT_FLAG: u32 = 0;
const HT_INT: u32 = 1;
const HT_REAL: u32 = 2;
const HT_STR: u32 = 3;

// Raw BCF type codes.
const BT_INT8: u8 = 1;
const BT_INT16: u8 = 2;
const BT_INT32: u8 = 3;
const BT_INT64: u8 = 4;
const BT_FLOAT: u8 = 5;
const BT_CHAR: u8 = 7;

/// The `Number=` names of the `BCF_VL_*` codes, for messages.
const CARDINALITIES: [&str; 10] = ["fixed", ".", "A", "G", "R", "P", "LA", "LG", "LR", "M"];

fn cardinality(vlen: u32) -> &'static str {
    CARDINALITIES.get(vlen as usize).copied().unwrap_or("?")
}

/// HTSlib's `bcf_hdr_id2length`: the `BCF_VL_*` code of an ID.
fn id_length(header: &Header, kind: HeaderRecordType, id: i32) -> u32 {
    u32::try_from(id)
        .ok()
        .and_then(|id| header.id_info(kind, id))
        .map_or(0, |info| (info >> 8) & 0xf)
}

/// HTSlib's `bcf_hdr_id2type`: the `BCF_HT_*` code of an ID.
fn id_type(header: &Header, kind: HeaderRecordType, id: i32) -> u32 {
    u32::try_from(id).map_or(15, |id| header.id_type_code(kind, id))
}

/// HTSlib's `bcf_hdr_int2id(BCF_DT_ID)`.
fn tag_name(header: &Header, id: i32) -> &[u8] {
    u32::try_from(id)
        .ok()
        .and_then(|id| header.id_name(id))
        .unwrap_or(b"")
}

/// HTSlib's negative return codes of `bcf_get_*_values`.
fn get_code(e: GetError) -> i64 {
    match e {
        GetError::NotInHeader => -1,
        GetError::WrongType => -2,
        GetError::Absent => -3,
        // HTSlib has no code for it (it reads past the record's values); -4 is its other
        // failure.
        GetError::TooFewSamples => -4,
    }
}

/// A value of an integer or float field, with its sentinels.
trait Value: Copy {
    const MISSING: Self;
    const VECTOR_END: Self;
    fn is_vector_end(self) -> bool;
}

impl Value for i32 {
    const MISSING: Self = INT32_MISSING;
    const VECTOR_END: Self = INT32_VECTOR_END;
    fn is_vector_end(self) -> bool {
        self == INT32_VECTOR_END
    }
}

impl Value for f32 {
    const MISSING: Self = f32::from_bits(FLOAT_MISSING);
    const VECTOR_END: Self = f32::from_bits(FLOAT_VECTOR_END);
    fn is_vector_end(self) -> bool {
        self.to_bits() == FLOAT_VECTOR_END
    }
}

/// Local-allele bookkeeping built from FORMAT/LAA.
struct LocalAlleles {
    /// Per sample (`stride` entries each): the new position of each original local allele
    /// (REF first), or -1 when removed.
    map: Vec<i32>,
    stride: usize,
    /// Per sample: the original number of local alleles, REF included.
    lr_orig: Vec<usize>,
}

/// Which alleles a FORMAT field's values follow.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum AlleleScope {
    /// The record's alleles (Number=A/R/G).
    Record,
    /// Each sample's local alleles, listed in FORMAT/LAA (Number=LA/LR/LG).
    Local,
}

impl AlleleScope {
    /// The scope of a field with `BCF_VL_*` length code `vlen`.
    fn of(vlen: u32) -> Self {
        if VL_LA_LG_LR & (1 << vlen) != 0 {
            Self::Local
        } else {
            Self::Record
        }
    }
}

/// What `bcf_remove_allele_set` works with once the alleles are renumbered.
struct Ctx<'a> {
    rm: &'a [bool],
    /// Old allele index to new, -1 when removed.
    map: Vec<i32>,
    chrom: String,
    pos: i64,
    n_r_ori: usize,
    n_r_new: usize,
    n_a_ori: usize,
    n_a_new: usize,
    n_g_ori: usize,
    n_sample: usize,
    local: Option<LocalAlleles>,
}

impl Ctx<'_> {
    /// HTSlib's `kbs_exists(rm_set, i)`.
    fn is_rm(&self, i: usize) -> bool {
        is_rm(self.rm, i)
    }

    /// The allele map for sample `j`: the local one for Number=LA/LG/LR fields.
    fn sample_map(&self, scope: AlleleScope, j: usize) -> &[i32] {
        match (&self.local, scope) {
            (Some(l), AlleleScope::Local) => &l.map[j * l.stride..(j + 1) * l.stride],
            _ => &self.map,
        }
    }

    fn lr_orig(&self, j: usize) -> usize {
        self.local.as_ref().map_or(0, |l| l.lr_orig[j])
    }

    /// Sample `j`'s original numbers of alleles (REF included) and of genotypes, for a
    /// Number=G or LG field.
    fn sample_n_ori(&self, scope: AlleleScope, j: usize) -> (usize, usize) {
        match scope {
            AlleleScope::Local => {
                let n_r = self.lr_orig(j);
                (n_r, n_r * (n_r + 1) / 2)
            }
            AlleleScope::Record => (self.n_r_ori, self.n_g_ori),
        }
    }

    /// Logs an error and returns it.
    fn fail(&self, msg: std::fmt::Arguments<'_>) -> RemoveAllelesError {
        let msg = msg.to_string();
        log::error(FUNC, format_args!("{msg}"));
        RemoveAllelesError(msg)
    }
}

fn is_rm(rm: &[bool], i: usize) -> bool {
    rm.get(i).copied().unwrap_or(false)
}

/// Reads entry `i` of an allele map; out-of-range entries count as removed.
fn map_at(map: &[i32], i: usize) -> i32 {
    map.get(i).copied().unwrap_or(-1)
}

/// Removes the alleles whose indices are set in `rm` (HTSlib's `bcf_remove_allele_set`),
/// updating Number=A/R/G INFO and FORMAT fields (and the local-allele LAA, LA/LR/LG fields),
/// the `<CNV:TR>` tandem-repeat INFO tags, and the genotypes. Index 0 is REF and should not
/// be set; indices past `rm` are kept.
///
/// # Errors
///
/// Returns an error (after logging it as HTSlib does) when a field has an unexpected number
/// of values or cannot be read or updated.
pub fn remove_allele_set(
    header: &Header,
    line: &mut Record,
    rm: &[bool],
) -> Result<(), RemoveAllelesError> {
    line.unpack(Unpack::ALL);

    // Map the old ALT numbering to the new one and build the new alleles.
    let n_allele = line.allele_count();
    if n_allele == 0 {
        return Ok(());
    }
    let mut map = vec![0i32; n_allele];
    let mut alleles: Vec<Vec<u8>> = Vec::with_capacity(n_allele);
    alleles.push(line.allele(0).to_vec());
    let mut have_cnv_tr = false;
    let mut nrm = 0;
    let mut j = 1;
    for (i, slot) in map.iter_mut().enumerate().skip(1) {
        let allele = line.allele(i);
        if allele == b"<CNV:TR>" {
            have_cnv_tr = true;
        }
        if is_rm(rm, i) {
            *slot = -1;
            nrm += 1;
            continue;
        }
        alleles.push(allele.to_vec());
        *slot = j;
        j += 1;
    }
    if nrm == 0 {
        return Ok(());
    }

    // REF is never removed, so at least one allele remains.
    let n_r_ori = n_allele;
    let n_r_new = n_allele - nrm;
    let n_a_ori = n_r_ori - 1;
    let n_a_new = n_r_new - 1;
    let n_g_ori = n_r_ori * (n_r_ori + 1) / 2;
    let n_g_new = n_r_new * (n_r_new + 1) / 2;

    line.update_alleles(header, &alleles);

    if have_cnv_tr {
        fixup_cnv_tr_info_tags(header, line, n_a_ori, rm);
    }

    let mut ctx = Ctx {
        rm,
        map,
        chrom: String::from_utf8_lossy(line.chrom_safe(header)).into_owned(),
        pos: line.pos() + 1,
        n_r_ori,
        n_r_new,
        n_a_ori,
        n_a_new,
        n_g_ori,
        n_sample: line.sample_count(),
        local: None,
    };

    remove_from_info(header, line, &ctx)?;
    update_genotypes(header, line, &ctx)?;
    ctx.local = update_local_alleles(header, line, &ctx)?;
    remove_from_format(header, line, &ctx, n_g_new)
}

/// Removes the alleles from Number=A/R/G INFO fields (and the Number=. fields that hold a
/// multiple of Number=A).
fn remove_from_info(
    header: &Header,
    line: &mut Record,
    ctx: &Ctx<'_>,
) -> Result<(), RemoveAllelesError> {
    let mut idat: Vec<i32> = Vec::new();
    let mut fdat: Vec<f32> = Vec::new();
    let mut sdat: Vec<u8> = Vec::new();
    let mut out: Vec<u8> = Vec::new();
    let mut i = 0;
    while i < line.info_count() {
        let key = line.d.info[i].key;
        i += 1;
        let mut vlen = id_length(header, HeaderRecordType::Info, key);
        let name = tag_name(header, key);
        let tag = String::from_utf8_lossy(name);

        // INFO types with Number=. that really hold a multiple of Number=A.
        let multiple = if vlen == VL_VAR {
            is_special_info_type(name)
        } else {
            1
        };
        if multiple > 1 {
            vlen = VL_A;
        }
        if VL_A_G_R & (1 << vlen) == 0 {
            continue;
        }

        let ty = id_type(header, HeaderRecordType::Info, key);
        if ty == HT_FLAG {
            continue;
        }

        let nret = match ty {
            HT_STR => line.get_info_string(header, name).map(|v| {
                sdat.clear();
                sdat.extend_from_slice(v);
                sdat.push(0);
                v.len()
            }),
            HT_INT => line.get_info_i32(header, name, &mut idat),
            HT_REAL => line.get_info_f32(header, name, &mut fdat),
            _ => Err(GetError::WrongType),
        };
        let nret = match nret {
            Ok(n) => n,
            Err(e) => {
                return Err(ctx.fail(format_args!(
                    "Could not access INFO/{tag} at {}:{} [{}]",
                    ctx.chrom,
                    ctx.pos,
                    get_code(e)
                )));
            }
        };
        if nret == 0 {
            continue; // no data for this tag
        }

        if ty == HT_STR {
            out.clear();
            if !remove_from_info_string(ctx, &sdat, &mut out, vlen, multiple, &tag)? {
                continue; // missing
            }
            let value = (!out.is_empty()).then_some(out.as_slice());
            if line.update_info_string(header, name, value).is_err() {
                return Err(ctx.fail(format_args!(
                    "Could not update INFO/{tag} at {}:{} [-1]",
                    ctx.chrom, ctx.pos
                )));
            }
            continue;
        }

        if nret == 1 {
            // The value could be missing: check.
            let field = &line.d.info[i - 1];
            let data = line.info_data(field);
            let missing = match field.ty {
                BT_INT8 => data.first().is_some_and(|&v| v as i8 == INT8_MISSING),
                BT_INT16 => data.len() >= 2 && le_i16(data) == INT16_MISSING,
                BT_INT32 => data.len() >= 4 && le_i32(data) == INT32_MISSING,
                BT_FLOAT => data.len() >= 4 && le_u32(data) == FLOAT_MISSING,
                other => return Err(ctx.fail(format_args!("Unexpected type {other}"))),
            };
            if missing {
                continue;
            }
        }

        let res = if ty == HT_INT {
            let ndat = remove_from_info_values(ctx, &mut idat, nret, vlen, multiple, &tag)?;
            line.update_info_i32(header, name, &idat[..ndat.min(idat.len())])
        } else {
            let ndat = remove_from_info_values(ctx, &mut fdat, nret, vlen, multiple, &tag)?;
            line.update_info_f32(header, name, &fdat[..ndat.min(fdat.len())])
        };
        if res.is_err() {
            return Err(ctx.fail(format_args!(
                "Could not update INFO/{tag} at {}:{} [-1]",
                ctx.chrom, ctx.pos
            )));
        }
    }
    Ok(())
}

/// Removes the alleles from a comma-separated string INFO value (`dat`, NUL-terminated) into
/// `out`. Returns `false` for a missing value, which is left alone.
fn remove_from_info_string(
    ctx: &Ctx<'_>,
    dat: &[u8],
    out: &mut Vec<u8>,
    vlen: u32,
    multiple: usize,
    tag: &str,
) -> Result<bool, RemoveAllelesError> {
    let at = |p: usize| dat.get(p).copied().unwrap_or(0);
    let s = at(0);
    let (mut ss, mut se) = (0usize, 0usize);
    if vlen == VL_A || vlen == VL_R {
        let (nexp, inc) = if vlen == VL_A {
            (ctx.n_a_ori * multiple, 1)
        } else {
            (ctx.n_r_ori, 0)
        };
        let mut j = 0;
        while j < nexp {
            if at(se) == 0 {
                break;
            }
            while at(se) != 0 && at(se) != b',' {
                se += 1;
            }
            if ctx.is_rm(j / multiple + inc) {
                if at(se) != 0 {
                    se += 1;
                }
                ss = se;
                j += 1;
                continue;
            }
            if !out.is_empty() {
                out.push(b',');
            }
            out.extend_from_slice(&dat[ss..se]);
            if at(se) != 0 {
                se += 1;
            }
            ss = se;
            j += 1;
        }
        if j == 1 && s == b'.' {
            return Ok(false);
        }
        if j != nexp {
            return Err(ctx.fail(format_args!(
                "Unexpected number of values in INFO/{tag} at {}:{}; expected Number={}={nexp}, but found {j}",
                ctx.chrom,
                ctx.pos,
                cardinality(vlen)
            )));
        }
    } else {
        // Number=G, assuming diploid genotypes.
        let mut n = 0;
        'outer: for j in 0..ctx.n_r_ori {
            for k in 0..=j {
                if at(se) == 0 {
                    break 'outer;
                }
                while at(se) != 0 && at(se) != b',' {
                    se += 1;
                }
                n += 1;
                if ctx.is_rm(j) || ctx.is_rm(k) {
                    if at(se) != 0 {
                        se += 1;
                    }
                    ss = se;
                    continue;
                }
                if !out.is_empty() {
                    out.push(b',');
                }
                out.extend_from_slice(&dat[ss..se]);
                if at(se) != 0 {
                    se += 1;
                }
                ss = se;
            }
            if at(se) == 0 {
                break;
            }
        }
        if n == 1 && s == b'.' {
            return Ok(false);
        }
        if n != ctx.n_g_ori {
            return Err(ctx.fail(format_args!(
                "Unexpected number of values in INFO/{tag} at {}:{}; expected Number=G={}, but found {n}",
                ctx.chrom, ctx.pos, ctx.n_g_ori
            )));
        }
    }
    Ok(true)
}

/// Removes the alleles from the `nret` values of a numeric INFO field in place, returning the
/// new number of values.
fn remove_from_info_values<T: Value>(
    ctx: &Ctx<'_>,
    dat: &mut [T],
    nret: usize,
    vlen: u32,
    multiple: usize,
    tag: &str,
) -> Result<usize, RemoveAllelesError> {
    if vlen == VL_A || vlen == VL_R {
        let (ntop, ndat, inc) = if vlen == VL_A {
            if nret != ctx.n_a_ori * multiple {
                return Err(ctx.fail(format_args!(
                    "Unexpected number of values in INFO/{tag} at {}:{}; expected Number={}={}, but found {nret}",
                    ctx.chrom,
                    ctx.pos,
                    if multiple == 1 { "A" } else { "." },
                    ctx.n_a_ori
                )));
            }
            (ctx.n_a_ori * multiple, ctx.n_a_new * multiple, 1)
        } else {
            if nret != ctx.n_r_ori {
                return Err(ctx.fail(format_args!(
                    "Unexpected number of values in INFO/{tag} at {}:{}; expected Number=R={}, but found {nret}",
                    ctx.chrom, ctx.pos, ctx.n_r_ori
                )));
            }
            (ctx.n_r_ori, ctx.n_r_new, 0)
        };
        let mut k = 0;
        for j in 0..ntop.min(dat.len()) {
            if dat[j].is_vector_end() {
                dat[k] = dat[j];
                break;
            }
            if ctx.is_rm(j / multiple + inc) {
                continue;
            }
            dat[k] = dat[j];
            k += 1;
        }
        Ok(ndat)
    } else {
        // Number=G
        if nret != ctx.n_g_ori {
            return Err(ctx.fail(format_args!(
                "Unexpected number of values in INFO/{tag} at {}:{}; expected Number=G={}, but found {nret}",
                ctx.chrom, ctx.pos, ctx.n_g_ori
            )));
        }
        let ndat = ctx.n_r_new * (ctx.n_r_new + 1) / 2;
        let mut l_new = 0;
        let mut l_ori = 0;
        for j in 0..ctx.n_r_ori {
            for k in 0..=j {
                let Some(&v) = dat.get(l_ori) else {
                    break;
                };
                l_ori += 1;
                if v.is_vector_end() {
                    dat[l_new] = v;
                    break;
                }
                if ctx.is_rm(j) || ctx.is_rm(k) {
                    continue;
                }
                dat[l_new] = v;
                l_new += 1;
            }
        }
        Ok(ndat)
    }
}

/// Renumbers the alleles in GT, turning removed ones into missing alleles.
fn update_genotypes(
    header: &Header,
    line: &mut Record,
    ctx: &Ctx<'_>,
) -> Result<(), RemoveAllelesError> {
    let renumbered = (1..ctx.n_r_ori).any(|i| ctx.map[i] != i as i32);
    if !renumbered {
        return Ok(());
    }
    let mut dat: Vec<i32> = Vec::new();
    let Ok(nret) = line.get_genotypes(header, &mut dat) else {
        return Ok(());
    };
    if nret == 0 || ctx.n_sample == 0 {
        return Ok(());
    }
    let nret = nret / ctx.n_sample;
    let samples = if nret > 0 { ctx.n_sample } else { 0 };
    for sample in dat.chunks_mut(nret.max(1)).take(samples) {
        for v in sample.iter_mut() {
            if *v >> 1 == 0 {
                continue; // missing
            }
            if *v == INT32_VECTOR_END {
                break;
            }
            let al = (*v >> 1) - 1;
            let mapped = usize::try_from(al)
                .ok()
                .filter(|&a| a < ctx.n_r_ori)
                .map(|a| ctx.map[a]);
            let Some(mapped) = mapped.filter(|&m| m >= -1) else {
                return Err(ctx.fail(format_args!(
                    "Problem updating genotypes at {}:{} [ al<nR_ori && map[al]>=-1 :: al={al},nR_ori={},map[al]={} ]",
                    ctx.chrom,
                    ctx.pos,
                    ctx.n_r_ori,
                    mapped.unwrap_or(-1)
                )));
            };
            // A removed allele becomes missing, keeping the phasing bit.
            *v = (if mapped < 0 { 0 } else { (mapped + 1) << 1 }) | (*v & 1);
        }
    }
    if line
        .update_genotypes(header, &dat[..nret * ctx.n_sample])
        .is_err()
    {
        return Err(ctx.fail(format_args!(
            "Could not update FORMAT/GT at {}:{} [-1]",
            ctx.chrom, ctx.pos
        )));
    }
    Ok(())
}

/// Removes the alleles from FORMAT/LAA and returns how the local alleles were renumbered, or
/// `None` without LAA.
fn update_local_alleles(
    header: &Header,
    line: &mut Record,
    ctx: &Ctx<'_>,
) -> Result<Option<LocalAlleles>, RemoveAllelesError> {
    let mut laa: Vec<i32> = Vec::new();
    let num_laa = match line.get_format_i32(header, b"LAA", &mut laa) {
        Ok(n) => n,
        Err(GetError::NotInHeader | GetError::Absent) => return Ok(None),
        Err(GetError::WrongType) => {
            return Err(RemoveAllelesError(
                "FORMAT/LAA has an unexpected type".into(),
            ));
        }
        Err(e @ GetError::TooFewSamples) => {
            return Err(RemoveAllelesError(format!("FORMAT/LAA: {e}")));
        }
    };
    if num_laa == 0 || ctx.n_sample == 0 {
        return Ok(None);
    }

    // Remove the LAA values of removed alleles, recording where the remaining ones went.
    let n_sample = ctx.n_sample;
    let num_laa_vals = num_laa / n_sample;
    let stride = num_laa_vals + 1;
    let mut map = vec![0i32; stride * n_sample];
    let mut lr_orig = vec![0usize; n_sample];
    let mut max_k = 0;
    let mut laa_changed = false;
    for i in 0..n_sample {
        let sample_laa = &mut laa[i * num_laa_vals..(i + 1) * num_laa_vals];
        let sample_map = &mut map[i * stride..(i + 1) * stride];
        sample_map[0] = 0;
        let mut j = 0;
        let mut k = 0;
        while j < num_laa_vals {
            let v = sample_laa[j];
            if v == INT32_VECTOR_END || v == INT32_MISSING {
                break;
            }
            let allele = usize::try_from(v)
                .ok()
                .filter(|&a| a > 0 && a < ctx.n_r_ori)
                .unwrap_or(0);
            if allele == 0 || ctx.map[allele] < 0 {
                sample_map[j + 1] = -1;
                laa_changed = true;
                j += 1;
                continue;
            }
            if allele as i32 != ctx.map[allele] {
                laa_changed = true;
            }
            sample_laa[k] = ctx.map[allele];
            k += 1;
            sample_map[j + 1] = k as i32;
            j += 1;
        }
        lr_orig[i] = j + 1;
        max_k = max_k.max(k);
        for m in &mut sample_map[j + 1..] {
            *m = -1;
        }
        for (kk, v) in sample_laa.iter_mut().enumerate().skip(k) {
            *v = if kk > 0 {
                INT32_VECTOR_END
            } else {
                INT32_MISSING
            };
        }
    }
    if laa_changed {
        let mut num_laa = num_laa;
        if max_k < num_laa_vals {
            // The maximum number of values has shrunk: consolidate.
            if max_k > 0 {
                for i in 1..n_sample {
                    laa.copy_within(i * num_laa_vals..i * num_laa_vals + max_k, i * max_k);
                }
                num_laa = n_sample * max_k;
            } else {
                // No values left: all the referenced alleles were removed. Store missing values
                // so LAA stays (removing it would invalidate LAD, LPL, etc.).
                laa[..n_sample].fill(INT32_MISSING);
                num_laa = n_sample;
            }
        }
        if line
            .update_format_i32(header, b"LAA", &laa[..num_laa])
            .is_err()
        {
            return Err(RemoveAllelesError("could not update FORMAT/LAA".into()));
        }
    }
    Ok(Some(LocalAlleles {
        map,
        stride,
        lr_orig,
    }))
}

/// Removes the alleles from Number=A/R/G and LA/LR/LG FORMAT fields, assuming haploid or
/// diploid genotypes.
fn remove_from_format(
    header: &Header,
    line: &mut Record,
    ctx: &Ctx<'_>,
    n_g_new: usize,
) -> Result<(), RemoveAllelesError> {
    let mut idat: Vec<i32> = Vec::new();
    let mut fdat: Vec<f32> = Vec::new();
    let mut sdat: Vec<u8> = Vec::new();
    let mut out: Vec<u8> = Vec::new();
    let mut i = 0;
    while i < line.format_count().min(line.d.fmt.len()) {
        let id = line.d.fmt[i].id;
        i += 1;
        let vlen = id_length(header, HeaderRecordType::Format, id);
        if VL_A_G_R_LA_LG_LR & (1 << vlen) == 0 {
            continue;
        }
        let scope = AlleleScope::of(vlen);

        let ty = id_type(header, HeaderRecordType::Format, id);
        if ty == HT_FLAG {
            continue;
        }
        let name = tag_name(header, id);
        let tag = String::from_utf8_lossy(name);

        let nret = match ty {
            HT_STR => line.get_format_chars(header, name).map(|(v, _)| {
                sdat.clear();
                sdat.extend_from_slice(v);
                v.len()
            }),
            HT_INT => line.get_format_i32(header, name, &mut idat),
            HT_REAL => line.get_format_f32(header, name, &mut fdat),
            _ => Err(GetError::WrongType),
        };
        let nret = match nret {
            Ok(n) => n,
            Err(e) => {
                return Err(ctx.fail(format_args!(
                    "Could not access FORMAT/{tag} at {}:{} [{}]",
                    ctx.chrom,
                    ctx.pos,
                    get_code(e)
                )));
            }
        };
        if nret == 0 || ctx.n_sample == 0 {
            continue; // no data for this tag
        }

        if ty == HT_STR {
            out.clear();
            remove_from_format_strings(
                ctx,
                &sdat,
                &mut out,
                nret / ctx.n_sample,
                vlen,
                scope,
                &tag,
            )?;
            if line.update_format_chars(header, name, &out).is_err() {
                return Err(ctx.fail(format_args!(
                    "Could not update FORMAT/{tag} at {}:{} [-1]",
                    ctx.chrom, ctx.pos
                )));
            }
            continue;
        }

        let nori = nret / ctx.n_sample;
        if nori == 1 && !(vlen == VL_A && nori == ctx.n_a_ori) {
            // All values may be missing: check.
            let field = &line.d.fmt[i - 1];
            let data = line.format_data(field);
            let size = field.size();
            let sample = |j: usize| data.get(j * size..).unwrap_or(&[]);
            let is_missing: fn(&[u8]) -> bool = match field.ty {
                BT_INT8 => |d| d.first().is_some_and(|&v| v as i8 == INT8_MISSING),
                BT_INT16 => |d| d.len() >= 2 && le_i16(d) == INT16_MISSING,
                BT_INT32 => |d| d.len() >= 4 && le_i32(d) == INT32_MISSING,
                BT_FLOAT => |d| d.len() >= 4 && le_u32(d) == FLOAT_MISSING,
                other => return Err(ctx.fail(format_args!("Unexpected type {other}"))),
            };
            if (0..ctx.n_sample).all(|j| is_missing(sample(j))) {
                continue;
            }
        }

        let res = if ty == HT_INT {
            let ndat = remove_from_format_values(ctx, &mut idat, nori, vlen, scope, n_g_new, &tag)?;
            line.update_format_i32(header, name, &idat[..ndat])
        } else {
            let ndat = remove_from_format_values(ctx, &mut fdat, nori, vlen, scope, n_g_new, &tag)?;
            line.update_format_f32(header, name, &fdat[..ndat])
        };
        if res.is_err() {
            return Err(ctx.fail(format_args!(
                "Could not update FORMAT/{tag} at {}:{} [-1]",
                ctx.chrom, ctx.pos
            )));
        }
    }
    Ok(())
}

/// Removes the alleles from a string FORMAT field (`size` bytes per sample in `dat`) into
/// `out`, keeping `size` bytes per sample.
fn remove_from_format_strings(
    ctx: &Ctx<'_>,
    dat: &[u8],
    out: &mut Vec<u8>,
    size: usize,
    vlen: u32,
    scope: AlleleScope,
    tag: &str,
) -> Result<(), RemoveAllelesError> {
    if VL_A_R_LA_LR & (1 << vlen) != 0 {
        remove_from_allele_strings(ctx, dat, out, size, vlen, scope, tag)
    } else {
        remove_from_genotype_strings(ctx, dat, out, size, vlen, scope, tag)
    }
}

/// [`remove_from_format_strings`] for Number=A, R, LA, and LR fields.
fn remove_from_allele_strings(
    ctx: &Ctx<'_>,
    dat: &[u8],
    out: &mut Vec<u8>,
    size: usize,
    vlen: u32,
    scope: AlleleScope,
    tag: &str,
) -> Result<(), RemoveAllelesError> {
    let at = |p: usize| dat.get(p).copied().unwrap_or(0);
    let (mut nexp, inc) = match vlen {
        VL_A => (ctx.n_a_ori, 1),
        VL_R => (ctx.n_r_ori, 0),
        _ => {
            if ctx.local.is_none() {
                return Err(no_laa(ctx, tag, vlen));
            }
            (0, usize::from(vlen == VL_LA))
        }
    };
    for j in 0..ctx.n_sample {
        let mut ss = j * size;
        let se = ss + size;
        let mut ptr = ss;
        let s = at(ss);
        let l = out.len();
        let sample_map = ctx.sample_map(scope, j);
        if scope == AlleleScope::Local {
            nexp = ctx.lr_orig(j).saturating_sub(inc);
        }
        let mut k_src = 0;
        let mut k_dst = 0;
        while k_src < nexp {
            if ptr >= se || at(ptr) == 0 {
                break;
            }
            while ptr < se && at(ptr) != 0 && at(ptr) != b',' {
                ptr += 1;
            }
            if map_at(sample_map, k_src + inc) < 0 {
                ptr += 1;
                ss = ptr;
                k_src += 1;
                continue;
            }
            if k_dst > 0 {
                out.push(b',');
            }
            out.extend_from_slice(&dat[ss..ptr]);
            ptr += 1;
            ss = ptr;
            k_dst += 1;
            k_src += 1;
        }
        if k_src != nexp && !(k_src == 1 && s == b'.') {
            return Err(ctx.fail(format_args!(
                "Unexpected number of values in FORMAT/{tag} at {}:{}; expected Number={}={nexp}, but found {k_src}",
                ctx.chrom,
                ctx.pos,
                cardinality(vlen)
            )));
        }
        for p in out.len() - l..size {
            out.push(if p == 0 { b'.' } else { 0 });
        }
    }
    Ok(())
}

/// [`remove_from_format_strings`] for Number=G and LG fields, diploid or haploid.
fn remove_from_genotype_strings(
    ctx: &Ctx<'_>,
    dat: &[u8],
    out: &mut Vec<u8>,
    size: usize,
    vlen: u32,
    scope: AlleleScope,
    tag: &str,
) -> Result<(), RemoveAllelesError> {
    let at = |p: usize| dat.get(p).copied().unwrap_or(0);
    for j in 0..ctx.n_sample {
        let mut ss = j * size;
        let se = ss + size;
        let mut ptr = ss;
        let s = at(ss);
        let l = out.len();
        let mut k_dst = 0;
        let (sample_n_r_ori, sample_n_g_ori) = ctx.sample_n_ori(scope, j);
        let sample_map = ctx.sample_map(scope, j);
        let mut nexp = 0;
        while ptr < se {
            if at(ptr) == 0 {
                break;
            }
            if at(ptr) == b',' {
                nexp += 1;
            }
            ptr += 1;
        }
        if ptr != ss {
            nexp += 1;
        }
        if nexp != sample_n_g_ori && nexp != sample_n_r_ori && !(nexp == 1 && s == b'.') {
            return Err(ctx.fail(format_args!(
                "Unexpected number of values in FORMAT/{tag} at {}:{}; expected Number={}={sample_n_g_ori}(diploid) or {sample_n_r_ori}(haploid), but found {nexp}",
                ctx.chrom,
                ctx.pos,
                cardinality(vlen)
            )));
        }
        ptr = ss;
        if nexp == 1 && s == b'.' {
            out.push(b'.');
        } else if nexp == sample_n_g_ori {
            // diploid
            'outer: for ia in 0..sample_n_r_ori {
                for ib in 0..=ia {
                    if ptr >= se || at(ptr) == 0 {
                        break 'outer;
                    }
                    while ptr < se && at(ptr) != 0 && at(ptr) != b',' {
                        ptr += 1;
                    }
                    if map_at(sample_map, ia) < 0 || map_at(sample_map, ib) < 0 {
                        ptr += 1;
                        ss = ptr;
                        continue;
                    }
                    if k_dst > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(&dat[ss..ptr]);
                    ptr += 1;
                    ss = ptr;
                    k_dst += 1;
                }
                if ptr >= se || at(ptr) == 0 {
                    break;
                }
            }
        } else {
            // haploid
            let mut k_src = 0;
            while k_src < sample_n_r_ori {
                if ptr >= se || at(ptr) == 0 {
                    break;
                }
                while ptr < se && at(ptr) != 0 && at(ptr) != b',' {
                    ptr += 1;
                }
                if map_at(sample_map, k_src) < 0 {
                    ptr += 1;
                    ss = ptr;
                    k_src += 1;
                    continue;
                }
                if k_dst > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(&dat[ss..ptr]);
                ptr += 1;
                ss = ptr;
                k_dst += 1;
                k_src += 1;
            }
            if k_src != ctx.n_r_ori {
                return Err(ctx.fail(format_args!(
                    "Unexpected number of values in FORMAT/{tag} at {}:{}; expected Number={}={}(haploid), but found {k_src}",
                    ctx.chrom,
                    ctx.pos,
                    cardinality(vlen),
                    ctx.n_r_ori
                )));
            }
        }
        let written = out.len() - l;
        if written < size {
            out.resize(l + size, 0);
        }
    }
    Ok(())
}

fn no_laa(ctx: &Ctx<'_>, tag: &str, vlen: u32) -> RemoveAllelesError {
    ctx.fail(format_args!(
        "No LAA data at {}:{}; required by FORMAT/{tag} with Number={}",
        ctx.chrom,
        ctx.pos,
        cardinality(vlen)
    ))
}

/// Removes the alleles from the values of a numeric FORMAT field (`nori` per sample) in
/// place, returning the new total number of values.
fn remove_from_format_values<T: Value>(
    ctx: &Ctx<'_>,
    dat: &mut Vec<T>,
    nori: usize,
    vlen: u32,
    scope: AlleleScope,
    n_g_new: usize,
    tag: &str,
) -> Result<usize, RemoveAllelesError> {
    if VL_A_R_LA_LR & (1 << vlen) != 0 || (vlen == VL_G && nori == ctx.n_r_ori) {
        remove_from_allele_values(ctx, dat, nori, vlen, scope, tag)
    } else {
        remove_from_genotype_values(ctx, dat, nori, scope, n_g_new, tag)
    }
}

/// [`remove_from_format_values`] for Number=A, R, LA, LR, and haploid Number=G fields.
fn remove_from_allele_values<T: Value>(
    ctx: &Ctx<'_>,
    dat: &mut [T],
    nori: usize,
    vlen: u32,
    scope: AlleleScope,
    tag: &str,
) -> Result<usize, RemoveAllelesError> {
    let (inc, nnew) = match vlen {
        VL_A => {
            if nori != ctx.n_a_ori {
                return Err(ctx.fail(format_args!(
                    "Unexpected number of values in FORMAT/{tag} at {}:{}; expected Number=A={}, but found {nori}",
                    ctx.chrom, ctx.pos, ctx.n_a_ori
                )));
            }
            (1, ctx.n_a_new)
        }
        VL_R | VL_G => {
            if vlen == VL_R && nori != ctx.n_r_ori {
                return Err(ctx.fail(format_args!(
                    "Unexpected number of values in FORMAT/{tag} at {}:{}; expected Number=R={}, but found {nori}",
                    ctx.chrom, ctx.pos, ctx.n_r_ori
                )));
            }
            (0, ctx.n_r_new)
        }
        _ => {
            if ctx.local.is_none() {
                return Err(no_laa(ctx, tag, vlen));
            }
            (usize::from(vlen == VL_LA), nori)
        }
    };
    for j in 0..ctx.n_sample {
        let src = j * nori;
        let dst = j * nnew;
        let sample_map = ctx.sample_map(scope, j);
        let sample_nori = match scope {
            AlleleScope::Local => ctx.lr_orig(j).saturating_sub(inc).min(nori),
            AlleleScope::Record => nori,
        };
        let mut k_dst = 0;
        for k_src in 0..sample_nori {
            let v = dat[src + k_src];
            if v.is_vector_end() {
                break;
            }
            if map_at(sample_map, k_src + inc) < 0 {
                continue;
            }
            dat[dst + k_dst] = v;
            k_dst += 1;
        }
        if k_dst == 0 {
            dat[dst] = T::MISSING;
            k_dst += 1;
        }
        while k_dst < nnew {
            dat[dst + k_dst] = T::VECTOR_END;
            k_dst += 1;
        }
    }
    Ok(nnew * ctx.n_sample)
}

/// [`remove_from_format_values`] for Number=G and LG fields, diploid or a mixture of haploid
/// and diploid.
fn remove_from_genotype_values<T: Value>(
    ctx: &Ctx<'_>,
    dat: &mut Vec<T>,
    nori: usize,
    scope: AlleleScope,
    n_g_new: usize,
    tag: &str,
) -> Result<usize, RemoveAllelesError> {
    let n_sample = ctx.n_sample;
    if scope == AlleleScope::Record && nori != ctx.n_g_ori {
        return Err(ctx.fail(format_args!(
            "Unexpected number of values in FORMAT/{tag} at {}:{}; expected Number=G={}, but found {nori}",
            ctx.chrom, ctx.pos, ctx.n_g_ori
        )));
    }
    // A Number=LG field keeps its own width. HTSlib stores it in the width the Number=G fields
    // that follow use too, so a PL after an LPL is cut to the LPL's width.
    let n_g_new = if scope == AlleleScope::Local {
        nori
    } else {
        n_g_new
    };
    let ndat = n_g_new * n_sample;
    if dat.len() < ndat {
        dat.resize(ndat, T::VECTOR_END);
    }
    let get = |dat: &[T], i: usize| dat.get(i).copied().unwrap_or(T::VECTOR_END);
    let set = |dat: &mut [T], i: usize, v: T| {
        if let Some(slot) = dat.get_mut(i) {
            *slot = v;
        }
    };
    for j in 0..n_sample {
        let src = j * nori;
        let dst = j * n_g_new;
        let sample_map = ctx.sample_map(scope, j);
        let (sample_n_r_ori, sample_n_g_ori) = ctx.sample_n_ori(scope, j);
        let nset = (0..sample_n_g_ori)
            .take_while(|&k| !get(dat, src + k).is_vector_end())
            .count();
        let mut k_dst = 0;
        if nset == sample_n_r_ori {
            // haploid
            for k_src in 0..sample_n_r_ori {
                if map_at(sample_map, k_src) < 0 {
                    continue;
                }
                let v = get(dat, src + k_src);
                set(dat, dst + k_dst, v);
                k_dst += 1;
            }
        } else {
            // diploid
            let mut k_src = 0;
            let mut ia = 0;
            while ia < sample_n_r_ori {
                for ib in 0..=ia {
                    let v = get(dat, src + k_src);
                    k_src += 1;
                    if v.is_vector_end() {
                        set(dat, dst + k_dst, v);
                        ia = ctx.n_r_ori;
                        break;
                    }
                    if map_at(sample_map, ia) < 0 || map_at(sample_map, ib) < 0 {
                        continue;
                    }
                    set(dat, dst + k_dst, v);
                    k_dst += 1;
                }
                ia += 1;
            }
        }
        while k_dst < n_g_new {
            set(dat, dst + k_dst, T::VECTOR_END);
            k_dst += 1;
        }
    }
    Ok(ndat)
}

/// Returns how many values per ALT allele an INFO tag with Number=. really holds (HTSlib's
/// `is_special_info_type`): 2 for the confidence intervals, 4 for the mobile elements.
fn is_special_info_type(name: &[u8]) -> usize {
    match name {
        b"CICN" | b"CIEND" | b"CILEN" | b"CIPOS" => 2,
        b"MEINFO" | b"METRANS" => 4,
        _ => 1,
    }
}

// ---- <CNV:TR> tags ----

/// The parts of an INFO field the `<CNV:TR>` helpers read (HTSlib keeps a pointer to it).
struct InfoSnapshot {
    ty: u8,
    len: i32,
    data: Vec<u8>,
}

/// Snapshots an INFO field by tag name (HTSlib's `bcf_get_info`).
fn info_snapshot(header: &Header, line: &Record, tag: &[u8]) -> Option<(usize, InfoSnapshot)> {
    let id = header
        .id(tag)
        .filter(|&id| header.id_exists(HeaderRecordType::Info, id))?;
    let i = line.d.info[..line.info_count()]
        .iter()
        .position(|f| f.key == id as i32)?;
    let field = &line.d.info[i];
    Some((
        i,
        InfoSnapshot {
            ty: field.ty,
            len: field.len,
            data: line.info_data(field).to_vec(),
        },
    ))
}

/// Reads value `index` of an integer or float INFO field as `i32` (HTSlib's
/// `get_int32_info_value`).
fn get_int32_info_value(info: &InfoSnapshot, index: usize) -> i32 {
    let len = info.len.max(0) as usize;
    if index >= len {
        return INT32_MISSING;
    }
    let d = &info.data;
    match info.ty {
        BT_INT8 => match d.get(index) {
            Some(&v) => {
                let val = v as i8;
                if val > INT8_VECTOR_END {
                    i32::from(val)
                } else {
                    INT32_VECTOR_END - (i32::from(INT8_VECTOR_END) - i32::from(val))
                }
            }
            None => INT32_MISSING,
        },
        BT_INT16 => match d.get(index * 2..index * 2 + 2) {
            Some(p) => {
                let val = le_i16(p);
                if val > INT16_VECTOR_END {
                    i32::from(val)
                } else {
                    INT32_VECTOR_END - (i32::from(INT16_VECTOR_END) - i32::from(val))
                }
            }
            None => INT32_MISSING,
        },
        BT_INT32 => d
            .get(index * 4..index * 4 + 4)
            .map_or(INT32_MISSING, le_i32),
        BT_FLOAT => match d.get(index * 4..index * 4 + 4) {
            Some(p) => {
                let bits = le_u32(p);
                if bits == FLOAT_MISSING {
                    INT32_MISSING
                } else if bits == FLOAT_VECTOR_END {
                    INT32_VECTOR_END
                } else {
                    f32::from_bits(bits) as i32
                }
            }
            None => INT32_MISSING,
        },
        _ => INT32_MISSING,
    }
}

/// The repeat count of ALT allele `index` from RN, one without RN (HTSlib's `get_rn_value`).
fn get_rn_value(rn: Option<&InfoSnapshot>, index: usize) -> i32 {
    let val = rn.map_or(1, |rn| get_int32_info_value(rn, index));
    val.max(0)
}

/// Updates `v1` after a field shrinks to one value (HTSlib's `set_info_v1`); `vptr` holds the
/// values.
fn set_info_v1(info: &mut InfoField, vptr: &[u8]) {
    match info.ty {
        BT_INT8 => info.v1 = i64::from(vptr[0] as i8),
        BT_INT16 => info.v1 = i64::from(le_i16(vptr)),
        BT_INT32 => info.v1 = i64::from(le_i32(vptr)),
        BT_INT64 => info.v1 = le_i64(vptr),
        // Only the float half of HTSlib's union is written.
        BT_FLOAT => info.v1 = (info.v1 & !0xffff_ffff) | i64::from(le_u32(vptr)),
        _ => {}
    }
}

/// What a `<CNV:TR>` helper did to a field.
enum Trim {
    /// Left unchanged (an unexpected type or count).
    Skipped,
    /// Trimmed in place.
    Trimmed,
    /// Every value went: mark the field for removal (HTSlib's `mark_for_removal`).
    Remove,
    /// Trimmed into a new buffer (key encoding, descriptor, values).
    Moved(Vec<u8>),
}

/// Rewrites the key and type descriptor before the values after `info.len` changed (HTSlib's
/// `fixup_info_length_code`). `block` starts at the key encoding. Returns a new buffer when
/// the descriptor grew.
fn fixup_info_length_code(info: &mut InfoField, block: &mut [u8]) -> Trim {
    let mut buf = Vec::with_capacity(24);
    enc_int1(&mut buf, info.key);
    enc_size(&mut buf, info.len.max(0) as usize, info.ty);
    let new_len = buf.len();
    let vptr_off = info.vptr_off as usize;
    let vptr_len = info.vptr_len as usize;
    if new_len <= vptr_off {
        // The length code kept its size or shrank: move the values down after it.
        block[..new_len].copy_from_slice(&buf);
        if new_len < vptr_off {
            block.copy_within(vptr_off..vptr_off + vptr_len, new_len);
            info.vptr_off = new_len as u32;
        }
        Trim::Trimmed
    } else {
        // Grown: this should not happen as values are only removed.
        buf.extend_from_slice(&block[vptr_off..vptr_off + vptr_len]);
        info.vptr_off = new_len as u32;
        Trim::Moved(buf)
    }
}

/// Runs a `<CNV:TR>` helper on INFO field `idx` with its bytes (from the key encoding on) and
/// applies the outcome. Returns whether the field moved to a new buffer.
fn with_info_block(
    line: &mut Record,
    idx: usize,
    f: impl FnOnce(&mut InfoField, &mut [u8]) -> Trim,
) -> bool {
    let Record { shared, d, .. } = line;
    let info = &mut d.info[idx];
    let mut data = std::mem::take(&mut info.data);
    let res = match &mut data {
        FieldData::InBlock(start) => {
            let start = *start;
            f(info, &mut shared[start..])
        }
        FieldData::Owned(buf) => f(info, buf),
        FieldData::Removed => Trim::Skipped,
    };
    match res {
        Trim::Remove => {
            info.data = FieldData::Removed;
            info.vptr_off = 0;
            info.vptr_len = 0;
            false
        }
        Trim::Moved(buf) => {
            info.data = FieldData::Owned(buf);
            true
        }
        Trim::Skipped | Trim::Trimmed => {
            info.data = data;
            false
        }
    }
}

/// Removes the values of removed alleles from an integer or float `<CNV:TR>` tag with
/// Number=. (HTSlib's `trim_int_cnv_tr_int_tags`): RB, RUL, CIRB, CIRUC, RUC, or RUB (with
/// `ruc`, whose counts give the items per repeat). Returns whether the field moved.
#[expect(clippy::too_many_arguments, reason = "mirrors HTSlib's helper")]
fn trim_int_cnv_tr_int_tags(
    header: &Header,
    line: &mut Record,
    idx: usize,
    rm: &[bool],
    id: &[u8],
    rn: Option<&InfoSnapshot>,
    ruc: Option<&InfoSnapshot>,
    num_alt_orig: usize,
    orig_total: i64,
) -> bool {
    let count: i64 = if id.first() == Some(&b'C') { 2 } else { 1 };
    let field = &line.d.info[idx];
    let ty = id_type(header, HeaderRecordType::Info, field.key);
    let vlen = id_length(header, HeaderRecordType::Info, field.key);
    const ELEMENT_SIZES: [i64; 8] = [0, 1, 2, 4, 0, 4, 0, 0];
    let element_size = ELEMENT_SIZES[usize::from(field.ty & 0x7)];

    // Give up on an unexpected type, Number, or count, or data that might be too short.
    if (ty != HT_INT && ty != HT_REAL)
        || element_size == 0
        || vlen != VL_VAR
        || i64::from(field.len) != orig_total
        || i64::from(field.vptr_len) < orig_total * element_size * count
    {
        return false;
    }

    with_info_block(line, idx, |info, block| {
        let vptr = info.vptr_off as usize;
        let mut unit = 0;
        let (mut orig_pos, mut new_pos) = (0usize, 0usize);
        let mut new_total: i64 = 0;
        for allele in 0..num_alt_orig {
            let mut n_repeats = get_rn_value(rn, allele);
            let mut n_items = i64::from(n_repeats);
            if let Some(ruc) = ruc {
                // For RUB, add up the units in RUC.
                n_items = 0;
                while n_repeats > 0 {
                    n_repeats -= 1;
                    let n_units = get_int32_info_value(ruc, unit);
                    unit += 1;
                    n_items += i64::from(n_units.max(0));
                }
            }
            let byte_len = (n_items * element_size * count) as usize;
            if is_rm(rm, allele + 1) {
                orig_pos += byte_len;
                continue;
            }
            if new_pos < orig_pos {
                block.copy_within(vptr + orig_pos..vptr + orig_pos + byte_len, vptr + new_pos);
            }
            orig_pos += byte_len;
            new_pos += byte_len;
            new_total += n_items;
        }
        if new_total == 0 {
            return Trim::Remove;
        }
        info.vptr_len = new_pos as u32;
        info.len = new_total as i32;
        if info.len == 1 {
            set_info_v1(info, &block[vptr..]);
        }
        fixup_info_length_code(info, block)
    })
}

/// Removes the values of removed alleles from the string `<CNV:TR>` tag RUS (HTSlib's
/// `trim_int_cnv_tr_str_tags`). Returns whether the field moved.
fn trim_int_cnv_tr_str_tags(
    header: &Header,
    line: &mut Record,
    idx: usize,
    rm: &[bool],
    rn: Option<&InfoSnapshot>,
    num_alt_orig: usize,
) -> bool {
    let field = &line.d.info[idx];
    let ty = id_type(header, HeaderRecordType::Info, field.key);
    let vlen = id_length(header, HeaderRecordType::Info, field.key);
    if ty != HT_STR || field.ty != BT_CHAR || vlen != VL_VAR {
        return false;
    }

    with_info_block(line, idx, |info, block| {
        let vptr = info.vptr_off as usize;
        let lim = info.vptr_len as usize;
        let (mut orig_pos, mut new_pos) = (0usize, 0usize);
        {
            let data = &mut block[vptr..vptr + lim];
            for allele in 0..num_alt_orig {
                let mut n_items = get_rn_value(rn, allele);
                let start = orig_pos;
                let mut end = start;
                while n_items > 0 {
                    n_items -= 1;
                    while end < lim && data[end] != 0 && data[end] != b',' {
                        end += 1;
                    }
                    if end == lim || data[end] == 0 {
                        break;
                    }
                    end += 1;
                }
                let len = end - start;
                if is_rm(rm, allele + 1) {
                    orig_pos += len;
                    continue;
                }
                if new_pos < orig_pos {
                    data.copy_within(orig_pos..orig_pos + len, new_pos);
                }
                orig_pos += len;
                new_pos += len;
            }
            if new_pos == 0 {
                return Trim::Remove;
            }
            if new_pos >= orig_pos {
                return Trim::Trimmed;
            }
            data[new_pos] = 0;
            // Dropping items at the end can leave a trailing comma: remove it.
            if data[new_pos - 1] == b',' {
                new_pos -= 1;
                data[new_pos] = 0;
            }
        }
        info.len = new_pos as i32;
        info.vptr_len = new_pos as u32;
        fixup_info_length_code(info, block)
    })
}

/// Removes the values of removed alleles from the INFO tags of `<CNV:TR>` (tandem repeat)
/// alleles (HTSlib's `fixup_cnv_tr_info_tags`):
///
/// - RN: repeat number, Number=A, handled as other tags of this type;
/// - RUS, RUL, RB, RUC: repeat unit sequence, unit length, sequence length, and unit count,
///   with RN values per allele;
/// - CIRB, CIRUC: confidence intervals around RB and RUC, two items per value;
/// - RUB: number of bases in each repeat unit, with RUC values per repeat.
fn fixup_cnv_tr_info_tags(header: &Header, line: &mut Record, num_alt_orig: usize, rm: &[bool]) {
    let rn = info_snapshot(header, line, b"RN").map(|(_, s)| s);
    let ruc = info_snapshot(header, line, b"RUC");

    // Total numbers of values of RUS etc. and of RUB, to check the counts found later.
    let mut orig_total_repeats: i64 = 0;
    let mut orig_total_units: i64 = 0;
    let mut unit = 0;
    for allele in 0..num_alt_orig {
        let mut n_repeats = get_rn_value(rn.as_ref(), allele);
        orig_total_repeats += i64::from(n_repeats);
        if let Some((_, ruc)) = &ruc {
            while n_repeats > 0 {
                n_repeats -= 1;
                let n_units = get_int32_info_value(ruc, unit);
                unit += 1;
                orig_total_units += i64::from(n_units.max(0));
            }
        }
    }

    for i in 0..line.info_count() {
        let id = tag_name(header, line.d.info[i].key);
        if !matches!(id.first(), Some(b'C' | b'R')) {
            continue;
        }
        // RUC is left for last, as RUB needs it intact. Tags with an unexpected type or count
        // are left unchanged.
        let moved = match id {
            b"RB" | b"RUL" | b"CIRB" | b"CIRUC" => trim_int_cnv_tr_int_tags(
                header,
                line,
                i,
                rm,
                id,
                rn.as_ref(),
                None,
                num_alt_orig,
                orig_total_repeats,
            ),
            b"RUS" => trim_int_cnv_tr_str_tags(header, line, i, rm, rn.as_ref(), num_alt_orig),
            b"RUB" if ruc.is_some() => trim_int_cnv_tr_int_tags(
                header,
                line,
                i,
                rm,
                id,
                rn.as_ref(),
                ruc.as_ref().map(|(_, s)| s),
                num_alt_orig,
                orig_total_units,
            ),
            _ => false,
        };
        // A removed field, or one moved to a new buffer, needs the INFO block rebuilt.
        if moved || line.d.info[i].is_removed() {
            line.d.shared_dirty |= DIRTY_INF;
        }
    }

    if let Some((idx, _)) = ruc {
        trim_int_cnv_tr_int_tags(
            header,
            line,
            idx,
            rm,
            b"RUC",
            rn.as_ref(),
            None,
            num_alt_orig,
            orig_total_repeats,
        );
    }
}

/// Removes the alleles in a bit mask (HTSlib's `bcf_remove_alleles`).
///
/// # Errors
///
/// Returns an error when a field cannot be updated.
pub fn remove_alleles(
    header: &Header,
    line: &mut Record,
    rm_mask: u64,
) -> Result<(), RemoveAllelesError> {
    let rm: Vec<bool> = (0..line.allele_count())
        .map(|i| i > 0 && i < 64 && rm_mask & (1 << i) != 0)
        .collect();
    remove_allele_set(header, line, &rm)
}

/// Removes ALT alleles no genotype uses (HTSlib's `bcf_trim_alleles`); returns how many were
/// removed (0 without GT).
///
/// # Errors
///
/// Returns an error for a genotype allele out of range or a failed removal.
pub fn trim_alleles(header: &Header, line: &mut Record) -> Result<usize, RemoveAllelesError> {
    let Some(gt) = super::gt_field(header, line) else {
        return Ok(0);
    };
    let n_allele = line.allele_count();
    let mut ac = vec![0u32; n_allele];
    for i in 0..line.sample_count() {
        let values = line.format_values(&gt, i);
        for ial in 0..gt.n() {
            let Some(val) = values.int(ial) else {
                crate::variant::log::error(
                    "bcf_trim_alleles",
                    format_args!(
                        "Unexpected GT {} at {}:{}",
                        gt.type_code(),
                        String::from_utf8_lossy(line.chrom_safe(header)),
                        line.pos() + 1
                    ),
                );
                return Ok(0);
            };
            if val == crate::variant::encoding::INT32_VECTOR_END {
                break;
            }
            if val >> 1 == 0 {
                continue;
            }
            let a = ((val >> 1) - 1) as usize;
            if a >= n_allele {
                crate::variant::log::error(
                    "bcf_trim_alleles",
                    format_args!(
                        "Allele index is out of bounds at {}:{}",
                        String::from_utf8_lossy(line.chrom_safe(header)),
                        line.pos() + 1
                    ),
                );
                return Err(RemoveAllelesError("allele index out of bounds".into()));
            }
            ac[a] += 1;
        }
    }
    let rm: Vec<bool> = (0..n_allele).map(|i| i > 0 && ac[i] == 0).collect();
    let nrm = rm.iter().filter(|&&r| r).count();
    if nrm > 0 {
        remove_allele_set(header, line, &rm)?;
    }
    Ok(nrm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::vcf;

    /// Parses `text` (header and records), trims each record's unused alleles, and returns
    /// the records formatted directly and after a BCF round trip, with the messages logged.
    fn trim_vcf(text: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
        let split = text.find("\n#CHROM").unwrap();
        let body_start = text[split + 1..].find('\n').unwrap() + split + 2;
        let mut header = Header::parse(&text.as_bytes()[..body_start]).unwrap();
        let (mut direct, mut synced) = (Vec::new(), Vec::new());
        let (_, messages) = log::capture(|| {
            for line in text[body_start..].lines() {
                let mut rec = Record::new();
                let mut buf = line.as_bytes().to_vec();
                vcf::parse(&mut header, &mut buf, &mut rec).unwrap();
                if trim_alleles(&header, &mut rec).is_err() {
                    continue;
                }
                let mut out = Vec::new();
                vcf::format_record(&header, &rec, &mut out).unwrap();
                direct.push(String::from_utf8(out).unwrap().trim_end().to_owned());
                let mut copy = Record::new();
                copy.copy_from(&mut rec);
                let mut out = Vec::new();
                vcf::format_record(&header, &copy, &mut out).unwrap();
                synced.push(String::from_utf8(out).unwrap().trim_end().to_owned());
            }
        });
        let messages = messages.iter().map(|m| format!("{m:?}")).collect();
        (direct, synced, messages)
    }

    const BASIC: &str = "##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">
##INFO=<ID=FL,Number=0,Type=Flag,Description=\"f\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"a\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">
##INFO=<ID=RC,Number=R,Type=Integer,Description=\"r\">
##INFO=<ID=GI,Number=G,Type=Integer,Description=\"g\">
##INFO=<ID=GF,Number=G,Type=Float,Description=\"g\">
##INFO=<ID=SA,Number=A,Type=String,Description=\"s\">
##INFO=<ID=SR,Number=R,Type=String,Description=\"s\">
##INFO=<ID=SG,Number=G,Type=String,Description=\"s\">
##INFO=<ID=CIPOS,Number=.,Type=Integer,Description=\"c\">
##INFO=<ID=MEINFO,Number=.,Type=String,Description=\"m\">
##INFO=<ID=VV,Number=.,Type=Integer,Description=\"v\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">
##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"p\">
##FORMAT=<ID=FA,Number=A,Type=Float,Description=\"a\">
##FORMAT=<ID=IA,Number=A,Type=Integer,Description=\"a\">
##FORMAT=<ID=XA,Number=A,Type=String,Description=\"s\">
##FORMAT=<ID=XR,Number=R,Type=String,Description=\"s\">
##FORMAT=<ID=XG,Number=G,Type=String,Description=\"s\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
";

    #[test]
    fn remove_alleles_basic() {
        let text = format!(
            "{BASIC}{}\n{}\n{}\n{}\n{}\n",
            "1\t100\tr1\tA\tC,G,T\t50\tPASS\tDP=10;FL;AC=1,2,3;AF=0.1,0.2,0.3;RC=5,6,7,8;\
             GI=0,1,2,3,4,5,6,7,8,9;GF=0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5;SA=a,b,c;\
             SR=r,a,b,c;SG=g0,g1,g2,g3,g4,g5,g6,g7,g8,g9;CIPOS=-1,1,-2,2,-3,3;\
             MEINFO=A,1,2,+,B,3,4,-,C,5,6,+;VV=1,2\tGT:AD:PL:GP:FA:IA:XA:XR:XG\t\
             0/1:1,2,3,4:0,1,2,3,4,5,6,7,8,9:0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9:0.1,0.2,0.3:\
             1,2,3:x,y,z:r,x,y,z:a,b,c,d,e,f,g,h,i,j\t\
             1|3:5,6,7,8:10,11,12,13,14,15,16,17,18,19:1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9:\
             .:.:.:.:.\t3:9,10,11,12:20,21,22,23:.:0.5,0.6,0.7:4,5,6:p,q,r:R,P,Q,S:h0,h1,h2,h3",
            "1\t200\tr2\tAC\tA,ACC\t.\t.\tAC=.;AF=.;RC=.;GI=.;SA=.;SR=.;SG=.\tGT:AD:PL:FA:XA\t\
             0/2:.:.:.:.\t./.:.:.:.:.\t0|.:.:.:.:.",
            "1\t300\tr3\tG\tT,C\t.\t.\tAC=0,1;AF=0,0.5;RC=3,0,1;GI=1,2,3,4,5,6;SA=u,v;SR=x,y,z\t\
             GT:AD:PL:FA\t0/0:3,0,0:0,10,20,30,40,50:0,0\t0/0:3,0,0:0,10,20,30,40,50:.,.\t.:.:.:.",
            "1\t400\tr4\tT\tA,C,G\t.\t.\tAC=1,1,1\tGT:PL\t0/1:0,1,2,3,4,5,6,7,8,9\t\
             2/3:0,1,2,3,4,5,6,7,8,9\t1/3:0,1,2,3,4,5,6,7,8,9",
            "1\t600\tr6\tT\tA,C\t.\t.\tAC=1,2\tGT:PL:AD\t0/2:0,1,2,3,4,5:1,2,3\t2:7,8,9:4,5,6\t\
             0/2:0,1,2,3:1,2",
        );
        // From `bcftools view -a -I` (HTSlib 1.24).
        let want = [
            "1\t100\tr1\tA\tC,T\t50\tPASS\tDP=10;FL;AC=1,3;AF=0.1,0.3;RC=5,6,8;GI=0,1,2,6,7,9;\
             GF=0.5,1.5,2.5,6.5,7.5,9.5;SA=a,c;SR=r,a,c;SG=g0,g1,g2,g6,g7,g9;CIPOS=-1,1,-3,3;\
             MEINFO=A,1,2,+,C,5,6,+;VV=1,2\tGT:AD:PL:GP:FA:IA:XA:XR:XG\t0/1:1,2,4:0,1,2,6,7,9:\
             0,0.1,0.2,0.6,0.7,0.9:0.1,0.3:1,3:x,z:r,x,z:a,b,c,g,h,j\t1|2:5,6,8:10,11,12,16,17,19:\
             1,1.1,1.2,1.6,1.7,1.9:.:.:.:.:.\t2:9,10,12:20,21,23:.:0.5,0.7:4,6:p,r:R,P,S:h0,h1,h3",
            "1\t200\tr2\tAC\tACC\t.\t.\tAC=.;AF=.;RC=.;GI=.;SA=.;SR=.;SG=.\tGT:AD:PL:FA:XA\t\
             0/1:.:.:.:.\t./.:.:.:.:.\t0|.:.:.:.:.",
            "1\t300\tr3\tG\t.\t.\t.\tRC=3;GI=1;SR=x\tGT:AD:PL\t0/0:3:0\t0/0:3:0\t.:.:.",
            "1\t400\tr4\tT\tA,C,G\t.\t.\tAC=1,1,1\tGT:PL\t0/1:0,1,2,3,4,5,6,7,8,9\t\
             2/3:0,1,2,3,4,5,6,7,8,9\t1/3:0,1,2,3,4,5,6,7,8,9",
            "1\t600\tr6\tT\tC\t.\t.\tAC=2\tGT:PL:AD\t0/1:0,3,5:1,3\t1:7,9:4,6\t0/1:0,3:1",
        ];
        let (direct, synced, messages) = trim_vcf(&text);
        assert_eq!(direct, want);
        assert_eq!(synced, want);
        assert!(messages.is_empty(), "{messages:?}");
    }

    #[test]
    fn unexpected_number_of_strings() {
        let text = format!(
            "{BASIC}1\t500\tr5\tT\tA,C\t.\t.\tSA=a,,c;AC=1,2\tGT:XA:XR\t0/2:,q:p,,r\t2:a,:,b,\t0:,:,,\n"
        );
        let (direct, _, messages) = trim_vcf(&text);
        assert!(direct.is_empty());
        assert_eq!(messages.len(), 1);
        assert!(
            messages[0].contains(
                "[E::bcf_remove_allele_set] Unexpected number of values in FORMAT/XA at 1:500; \
                 expected Number=A=2, but found 1"
            ),
            "{messages:?}"
        );
    }

    const CNV_TR: &str = "##fileformat=VCFv4.5
##contig=<ID=1,length=100000>
##INFO=<ID=RN,Number=A,Type=Integer,Description=\"Repeat number\">
##INFO=<ID=RUS,Number=.,Type=String,Description=\"Repeat unit sequence\">
##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"Repeat unit length\">
##INFO=<ID=RB,Number=.,Type=Integer,Description=\"Repeat bases\">
##INFO=<ID=RUC,Number=.,Type=Float,Description=\"Repeat unit count\">
##INFO=<ID=RUB,Number=.,Type=Integer,Description=\"Repeat unit bases\">
##INFO=<ID=CIRB,Number=.,Type=Integer,Description=\"CI RB\">
##INFO=<ID=CIRUC,Number=.,Type=Float,Description=\"CI RUC\">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"len\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
";

    #[test]
    fn remove_alleles_cnv_tr() {
        let tr3 = "<CNV:TR>,<CNV:TR>,<CNV:TR>\t.\t.\tRN=1,2,1;RUS=CAG,CAG,TT,A;RUL=3,3,2,1;\
                   RB=30,15,8,5;RUC=2,1,3,2;RUB=3,3,3,3,3,3,2,1;CIRB=-1,1,-2,2,-3,3,-4,4;\
                   CIRUC=-0.5,0.5,-1,1,-1.5,1.5,-2,2;SVLEN=30,23,5\tGT:AD";
        let text = format!(
            "{CNV_TR}\
             1\t100\tt1\tA\t{tr3}\t0/1:1,2,3,4\t3/3:5,6,7,8\n\
             1\t200\tt2\tA\t{tr3}\t0/1:1,2,3,4\t2/2:5,6,7,8\n\
             1\t300\tt3\tA\t{tr3}\t0/3:1,2,3,4\t3/3:5,6,7,8\n\
             1\t400\tt4\tA\t<CNV:TR>,<CNV:TR>\t.\t.\tRN=1,1;RUS=CAG,TT;RUL=3,2;RB=30,8;RUC=10,4;\
             RUB=3,3,3,3,3,3,3,3,3,3,2,2,2,2;CIRB=-1,1,-2,2;CIRUC=-1,1,-2,2;SVLEN=30,8\tGT:AD\t\
             0/1:1,2,3\t1/1:5,6,7\n\
             1\t500\tt5\tA\t<CNV:TR>,<CNV:TR>\t.\t.\tRN=1,1;RUS=CAG,T;RUL=3,1;RB=30,8;RUC=10,8;\
             CIRB=-1,1,-2,2\tGT\t0/2\t2/2\n\
             1\t600\tt6\tA\t<CNV:TR>,<DEL>\t.\t.\tRN=2,.;RUS=CAG,CAGG;RUL=3,4;RB=30,40\tGT\t0/2\t2/2\n\
             1\t700\tt7\tA\t<CNV:TR>,<CNV:TR>\t.\t.\tRN=1,1;RUS=CAG,T;RB=30,8,9\tGT\t0/1\t1/1\n"
        );
        // From `bcftools view -a -I` (HTSlib 1.24). A string shortened to one character is
        // printed from the stale `v1` (a NUL) until the record is re-read from BCF.
        let mut want = vec![
            "1\t100\tt1\tA\t<CNV:TR>,<CNV:TR>\t.\t.\tRN=1,1;RUS=CAG,A;RUL=3,1;RB=30,5;RUC=2,2;\
             RUB=3,3,2,1;CIRB=-1,1,-2,2,-3,3,-4,4;CIRUC=-0.5,0.5,-1,1,-1.5,1.5,-2,2;SVLEN=30,5\t\
             GT:AD\t0/1:1,2,4\t2/2:5,6,8"
                .to_owned(),
            "1\t200\tt2\tA\t<CNV:TR>,<CNV:TR>\t.\t.\tRN=1,2;RUS=CAG,CAG,TT;RUL=3,3,2;RB=30,15,8;\
             RUC=2,1,3;RUB=3,3,3,3,3,3;CIRB=-1,1,-2,2,-3,3,-4,4;CIRUC=-0.5,0.5,-1,1,-1.5,1.5,-2,2;\
             SVLEN=30,23\tGT:AD\t0/1:1,2,3\t2/2:5,6,7"
                .to_owned(),
            "1\t300\tt3\tA\t<CNV:TR>\t.\t.\tRN=1;RUS=\0;RUL=1;RB=5;RUC=2;RUB=2,1;\
             CIRB=-1,1,-2,2,-3,3,-4,4;CIRUC=-0.5,0.5,-1,1,-1.5,1.5,-2,2;SVLEN=5\tGT:AD\t0/1:1,4\t\
             1/1:5,8"
                .to_owned(),
            "1\t400\tt4\tA\t<CNV:TR>\t.\t.\tRN=1;RUS=CAG;RUL=3;RB=30;RUC=10;\
             RUB=3,3,3,3,3,3,3,3,3,3;CIRB=-1,1,-2,2;CIRUC=-1,1,-2,2;SVLEN=30\tGT:AD\t0/1:1,2\t\
             1/1:5,6"
                .to_owned(),
            "1\t500\tt5\tA\t<CNV:TR>\t.\t.\tRN=1;RUS=\0;RUL=1;RB=8;RUC=8;CIRB=-1,1,-2,2\tGT\t\
             0/1\t1/1"
                .to_owned(),
            "1\t600\tt6\tA\t<DEL>\t.\t.\tRN=.\tGT\t0/1\t1/1".to_owned(),
            "1\t700\tt7\tA\t<CNV:TR>\t.\t.\tRN=1;RUS=CAG;RB=30,8,9\tGT\t0/1\t1/1".to_owned(),
        ];
        let (direct, synced, messages) = trim_vcf(&text);
        assert_eq!(direct, want);
        want[2] = want[2].replace("RUS=\0", "RUS=A");
        want[4] = want[4].replace("RUS=\0", "RUS=T");
        assert_eq!(synced, want);
        assert!(messages.is_empty(), "{messages:?}");
    }

    const LOCAL: &str = "##fileformat=VCFv4.5
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"local alleles\">
##FORMAT=<ID=LAD,Number=LR,Type=Integer,Description=\"a\">
##FORMAT=<ID=LPL,Number=LG,Type=Integer,Description=\"p\">
##FORMAT=<ID=LAF,Number=LA,Type=Float,Description=\"f\">
##FORMAT=<ID=LXS,Number=LR,Type=String,Description=\"s\">
##FORMAT=<ID=LGS,Number=LG,Type=String,Description=\"s\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
";

    #[test]
    fn remove_alleles_local() {
        let text = format!(
            "{LOCAL}\
             1\t100\tl1\tA\tC,G,T\t.\t.\t.\tGT:LAA:LAD:LPL:LAF:LXS:LGS:AD\t\
             0/1:1,3:5,6,7:0,1,2,3,4,5:0.1,0.3:a,b,c:p,q,r,s,t,u:1,2,3,4\t\
             0/2:2:5,6:0,1,2:0.2:a,b:p,q,r:1,2,3,4\t0/0:.:9:0:.:a:p:1,2,3,4\n\
             1\t200\tl2\tA\tC,G,T\t.\t.\t.\tGT:LAA:LAD:LPL:AD\t0/1:1,2:5,6,7:0,1,2,3,4,5:1,2,3,4\t\
             0/1:1:5,6:0,1,2:1,2,3,4\t1/1:1,3:9,8,7:0,1,2,3,4,5:1,2,3,4\n\
             1\t300\tl3\tA\tC,G,T\t.\t.\t.\tGT:LAA:LAD:LPL\t0/0:2:5,6:0,1,2\t\
             0/0:3,2:5,6,7:0,1,2,3,4,5\t0/0:.:9:0\n\
             1\t400\tl4\tA\tC,G\t.\t.\t.\tGT:LAD\t0/1:1,2\t0/0:3,4\t0/0:1,1\n"
        );
        // From `bcftools view -a -I` (HTSlib 1.24).
        let want = [
            "1\t100\tl1\tA\tC,G\t.\t.\t.\tGT:LAA:LAD:LPL:LAF:LXS:LGS:AD\t\
             0/1:1:5,6:0,1,2:0.1:a,b:p,q,r:1,2,3\t0/2:2:5,6:0,1,2:0.2:a,b:p,q,r:1,2,3\t\
             0/0:.:9:0:.:a:p:1,2,3",
            "1\t200\tl2\tA\tC\t.\t.\t.\tGT:LAA:LAD:LPL:AD\t0/1:1:5,6:0,1,2:1,2\t0/1:1:5,6:0,1,2:1,2\t\
             1/1:1:9,8:0,1,2:1,2",
            "1\t300\tl3\tA\t.\t.\t.\t.\tGT:LAA:LAD:LPL\t0/0:.:5:0\t0/0:.:5:0\t0/0:.:9:0",
        ];
        let (direct, synced, messages) = trim_vcf(&text);
        assert_eq!(direct, want);
        assert_eq!(synced, want);
        assert_eq!(messages.len(), 1);
        assert!(
            messages[0].contains(
                "[E::bcf_remove_allele_set] No LAA data at 1:400; required by FORMAT/LAD with \
                 Number=LR"
            ),
            "{messages:?}"
        );
    }

    #[test]
    fn remove_alleles_keeps_global_genotype_values_after_local_ones() {
        // HTSlib sizes the Number=G fields after a Number=LG field by the LG field's width,
        // so PL lost three of its six values in either order.
        let header = "##fileformat=VCFv4.5
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"local alleles\">
##FORMAT=<ID=LPL,Number=LG,Type=Integer,Description=\"p\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
";
        let pl = "0,1,2,3,4,5,6,7,8,9";
        for format in ["GT:LAA:LPL:PL", "GT:PL:LAA:LPL"] {
            let values = |gt: &str| match format {
                "GT:LAA:LPL:PL" => format!("{gt}:1:0,1,2:{pl}"),
                _ => format!("{gt}:{pl}:1:0,1,2"),
            };
            let text = format!(
                "{header}1\t100\tg1\tA\tC,G,T\t.\t.\t.\t{format}\t{}\t{}\n",
                values("0/1"),
                values("0/2")
            );
            let (direct, synced, _) = trim_vcf(&text);
            for line in [&direct[0], &synced[0]] {
                let fields: Vec<&str> = line.split('\t').collect();
                assert_eq!(fields[4], "C,G", "{format}: {line}");
                let keys: Vec<&str> = fields[8].split(':').collect();
                let at = keys.iter().position(|&k| k == "PL").expect("PL");
                for sample in &fields[9..] {
                    let pl: Vec<&str> = sample
                        .split(':')
                        .nth(at)
                        .expect("PL value")
                        .split(',')
                        .collect();
                    assert_eq!(pl, ["0", "1", "2", "3", "4", "5"], "{format}: {line}");
                }
            }
        }
    }

    #[test]
    fn remove_alleles_by_mask() {
        let text = format!(
            "{BASIC}1\t100\tr1\tA\tC,G,T\t50\tPASS\tAC=1,2,3;RC=5,6,7,8\tGT:PL\t\
             0/1:0,1,2,3,4,5,6,7,8,9\t1|3:10,11,12,13,14,15,16,17,18,19\t\
             3:20,21,22,23\n"
        );
        let split = text.find("\n#CHROM").unwrap();
        let body_start = text[split + 1..].find('\n').unwrap() + split + 2;
        let mut header = Header::parse(&text.as_bytes()[..body_start]).unwrap();
        let mut rec = Record::new();
        let mut buf = text[body_start..].trim_end().as_bytes().to_vec();
        vcf::parse(&mut header, &mut buf, &mut rec).unwrap();
        // Removing nothing changes nothing.
        remove_alleles(&header, &mut rec, 1).unwrap();
        // Remove C and T.
        remove_alleles(&header, &mut rec, 0b1010).unwrap();
        let mut out = Vec::new();
        vcf::format_record(&header, &rec, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "1\t100\tr1\tA\tG\t50\tPASS\tAC=2;RC=5,7\tGT:PL\t0/.:0,3,5\t.|.:10,13,15\t\
             .:20,22\n"
        );
    }
}
