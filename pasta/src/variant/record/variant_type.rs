//! Variant classification of ALT alleles (HTSlib's `bcf_set_variant_type(s)` and the
//! `bcf_has_variant_type(s)` family).

use bitflags::bitflags;

use super::{Record, Unpack};

bitflags! {
    /// Variant types (HTSlib's `VCF_*`). The empty set is `VCF_REF`.
    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
    pub struct VariantType: u32 {
        /// A single-base substitution (`VCF_SNP`).
        const SNP = 1;
        /// A multi-base substitution (`VCF_MNP`).
        const MNP = 1 << 1;
        /// An insertion or deletion (`VCF_INDEL`).
        const INDEL = 1 << 2;
        /// Anything else, including symbolic alleles (`VCF_OTHER`).
        const OTHER = 1 << 3;
        /// A breakend (`VCF_BND`).
        const BND = 1 << 4;
        /// An overlapping deletion, ALT `*` (`VCF_OVERLAP`).
        const OVERLAP = 1 << 5;
        /// An insertion; implies INDEL (`VCF_INS`).
        const INS = 1 << 6;
        /// A deletion; implies INDEL (`VCF_DEL`).
        const DEL = 1 << 7;
    }
}

impl VariantType {
    /// The types HTSlib's older `bcf_get_variant_type(s)` report (without INS and DEL).
    pub const ORIGINAL: Self = Self::SNP
        .union(Self::MNP)
        .union(Self::INDEL)
        .union(Self::OTHER)
        .union(Self::BND)
        .union(Self::OVERLAP);
}

/// The classification of one allele (HTSlib's `bcf_variant_t`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Hash)]
pub struct Variant {
    /// The variant type; empty for a reference allele.
    pub ty: VariantType,
    /// The number of bases affected, negative for deletions.
    pub n: i32,
}

/// How [`Record::has_variant_types`] compares types (HTSlib's `bcf_variant_match`).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum VariantMatch {
    /// All types present must be in the mask (`bcf_match_subset`).
    Subset,
    /// The types present must equal the mask (`bcf_match_exact`).
    Exact,
    /// Any type present is in the mask (`bcf_match_overlap`).
    Overlap,
}

/// Classifies `alt` against `ref_` (HTSlib's `bcf_set_variant_type`).
#[must_use]
pub fn classify(ref_: &[u8], alt: &[u8]) -> Variant {
    let r = |i: usize| ref_.get(i).copied().unwrap_or(0);
    let a = |i: usize| alt.get(i).copied().unwrap_or(0);
    let refv = |n| Variant {
        ty: VariantType::empty(),
        n,
    };
    let make = |ty, n| Variant { ty, n };
    if a(0) == b'*' && a(1) == 0 {
        return make(VariantType::OVERLAP, 0);
    }
    if r(1) == 0 && a(1) == 0 {
        if a(0) == b'.' || r(0) == a(0) {
            return refv(0);
        }
        if a(0) == b'X' {
            return refv(0);
        }
        return make(VariantType::SNP, 1);
    }
    if a(0) == b'<' {
        if a(1) == b'X' && a(2) == b'>' {
            return refv(0);
        }
        if a(1) == b'*' && a(2) == b'>' {
            return refv(0);
        }
        if &alt[1..] == b"NON_REF>" {
            return refv(0);
        }
        return make(VariantType::OTHER, 0);
    }
    if a(0) == b']' || a(0) == b'[' {
        return make(VariantType::BND, 0);
    }
    let up = |c: u8| c.to_ascii_uppercase();
    let mut ri = 0;
    let mut ai = 0;
    while r(ri) != 0 && a(ai) != 0 && up(r(ri)) == up(a(ai)) {
        ri += 1;
        ai += 1;
    }
    if a(ai) != 0 && r(ri) == 0 {
        while a(ai) != 0 {
            ai += 1;
        }
        if a(ai - 1) == b']' || a(ai - 1) == b'[' {
            return make(VariantType::BND, 0);
        }
        return make(VariantType::INDEL | VariantType::INS, ai as i32 - ri as i32);
    } else if r(ri) != 0 && a(ai) == 0 {
        while r(ri) != 0 {
            ri += 1;
        }
        return make(VariantType::INDEL | VariantType::DEL, ai as i32 - ri as i32);
    } else if r(ri) == 0 && a(ai) == 0 {
        return refv(0);
    }
    let (r0, a0) = (ri, ai);
    let mut re = ri;
    let mut ae = ai;
    while r(re + 1) != 0 {
        re += 1;
    }
    while a(ae + 1) != 0 {
        ae += 1;
    }
    if a(ae) == b']' || a(ae) == b'[' {
        return make(VariantType::BND, 0);
    }
    while re > r0 && ae > a0 && up(r(re)) == up(a(ae)) {
        re -= 1;
        ae -= 1;
    }
    if ae == a0 {
        if re == r0 {
            return make(VariantType::SNP, 1);
        }
        let n = -((re - r0) as i32);
        if up(r(re)) == up(a(ae)) {
            return make(VariantType::INDEL | VariantType::DEL, n);
        }
        return make(VariantType::OTHER, n);
    } else if re == r0 {
        let n = (ae - a0) as i32;
        if up(r(re)) == up(a(ae)) {
            return make(VariantType::INDEL | VariantType::INS, n);
        }
        return make(VariantType::OTHER, n);
    }
    let (rl, al) = ((re - r0) as i32, (ae - a0) as i32);
    let ty = if rl == al {
        VariantType::MNP
    } else {
        VariantType::OTHER
    };
    let n = if rl > al { -(rl + 1) } else { al + 1 };
    make(ty, n)
}

impl Record {
    /// Classifies all alleles (HTSlib's `bcf_set_variant_types`), if not done since the
    /// alleles last changed.
    #[inline]
    fn ensure_variant_types(&mut self) {
        if self.d.var_type == -1 {
            self.set_variant_types();
        }
    }

    /// Classifies all alleles (HTSlib's `bcf_set_variant_types`).
    fn set_variant_types(&mut self) {
        if self.unpacked & Unpack::STR.bits() == 0 {
            self.unpack(Unpack::STR);
        }
        let n = self.n_allele as usize;
        self.d.var.clear();
        self.d.var.push(Variant::default());
        let mut all = VariantType::empty();
        for i in 1..n {
            let v = classify(self.allele(0), self.allele(i));
            all |= v.ty;
            self.d.var.push(v);
        }
        self.d.var_type = all.bits() as i32;
    }

    /// Returns the union of the alleles' variant types, without INS/DEL (HTSlib's
    /// `bcf_get_variant_types`).
    #[inline]
    pub fn variant_types(&mut self) -> VariantType {
        self.ensure_variant_types();
        VariantType::from_bits_truncate(self.d.var_type as u32) & VariantType::ORIGINAL
    }

    /// Returns the union of the alleles' variant types, including INS/DEL.
    pub fn variant_types_all(&mut self) -> VariantType {
        self.ensure_variant_types();
        VariantType::from_bits_truncate(self.d.var_type as u32)
    }

    /// Returns the classification of allele `i`, or `None` when out of range (HTSlib's
    /// `bcf_get_variant_type` and `bcf_variant_length`).
    #[inline]
    pub fn variant(&mut self, i: usize) -> Option<Variant> {
        self.ensure_variant_types();
        self.d.var.get(i).copied()
    }

    /// Returns whether allele `i` has a type in `mask`, or is a reference allele when `mask` is
    /// empty (HTSlib's `bcf_has_variant_type`).
    #[inline]
    pub fn has_variant_type(&mut self, i: usize, mask: VariantType) -> Option<bool> {
        let v = self.variant(i)?;
        Some(if mask.is_empty() {
            v.ty.is_empty()
        } else {
            v.ty.intersects(mask)
        })
    }

    /// Tests the record's variant types against `mask` (HTSlib's `bcf_has_variant_types`). An
    /// empty `mask` with [`VariantMatch::Exact`] tests for a reference-only site.
    pub fn has_variant_types(&mut self, mask: VariantType, mode: VariantMatch) -> HasVariantTypes {
        self.ensure_variant_types();
        let mut ty = VariantType::from_bits_truncate(self.d.var_type as u32);
        if mode == VariantMatch::Overlap {
            return HasVariantTypes(!(mask & ty).is_empty());
        }
        let indel_parts = VariantType::INS | VariantType::DEL;
        if mask.intersects(indel_parts) && !mask.contains(VariantType::INDEL) {
            ty &= !VariantType::INDEL;
        } else if mask.contains(VariantType::INDEL) && !mask.intersects(indel_parts) {
            ty &= !indel_parts;
        }
        if mode == VariantMatch::Subset {
            if !(ty - mask).is_empty() {
                return HasVariantTypes(false);
            }
            return HasVariantTypes(!(mask & ty).is_empty());
        }
        if mask.is_empty() {
            return HasVariantTypes(ty.is_empty());
        }
        HasVariantTypes(ty == mask && !ty.is_empty())
    }
}

/// The result of [`Record::has_variant_types`]: whether the record matched.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub struct HasVariantTypes(pub bool);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classification() {
        assert_eq!(classify(b"A", b"C").ty, VariantType::SNP);
        assert_eq!(
            classify(b"A", b"AT").ty,
            VariantType::INDEL | VariantType::INS
        );
        assert_eq!(
            classify(b"AT", b"A"),
            Variant {
                ty: VariantType::INDEL | VariantType::DEL,
                n: -1
            }
        );
        assert_eq!(classify(b"AC", b"GT").ty, VariantType::MNP);
        assert_eq!(classify(b"A", b"<DEL>").ty, VariantType::OTHER);
        assert_eq!(classify(b"A", b"<*>").ty, VariantType::empty());
        assert_eq!(classify(b"A", b"*").ty, VariantType::OVERLAP);
        assert_eq!(classify(b"A", b"A[1:10[").ty, VariantType::BND);
    }
}
