//! `%VKX`: the 64-bit VariantKey of a variant, from bcftools' `variantkey.h` (MIT licensed,
//! by Nicola Asuni and GENOMICS plc).

/// Bit offset of the chromosome code.
const SHIFT_CHROM: u32 = 59;
/// Bit offset of the position.
const SHIFT_POS: u32 = 31;

/// Encodes a numeric chromosome name (`encode_numeric_chrom`); 0 when not a number.
fn encode_numeric_chrom(chrom: &[u8]) -> u8 {
    let mut v = chrom[0].wrapping_sub(b'0');
    for &c in &chrom[1..] {
        if !c.is_ascii_digit() {
            return 0;
        }
        v = v.wrapping_mul(10).wrapping_add(c - b'0');
    }
    v
}

/// Encodes a chromosome name (`encode_chrom`): 1-22, X=23, Y=24, M/MT=25, 0 otherwise.
fn encode_chrom(mut chrom: &[u8]) -> u8 {
    if chrom.len() > 3 && chrom[..3].eq_ignore_ascii_case(b"chr") {
        chrom = &chrom[3..];
    }
    let Some(&c0) = chrom.first() else {
        return 0;
    };
    if c0.is_ascii_digit() {
        return encode_numeric_chrom(chrom);
    }
    if chrom.len() == 1 || (chrom.len() == 2 && (chrom[1] == b'T' || chrom[1] == b't')) {
        return match c0 {
            b'X' | b'x' => 23,
            b'Y' | b'y' => 24,
            b'M' | b'm' => 25,
            _ => 0,
        };
    }
    0
}

fn encode_base(c: u8) -> u32 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

fn encode_allele(h: &mut u32, bitpos: &mut u8, s: &[u8]) -> bool {
    for &c in s {
        let v = encode_base(c);
        if v > 3 {
            return false;
        }
        *bitpos -= 2;
        *h |= v << *bitpos;
    }
    true
}

/// The reversible encoding of short alleles (`encode_refalt_rev`), or `None`.
fn encode_refalt_rev(r: &[u8], a: &[u8]) -> Option<u32> {
    let mut h = ((r.len() as u32) << 27) | ((a.len() as u32) << 23);
    let mut bitpos = 23u8;
    (encode_allele(&mut h, &mut bitpos, r) && encode_allele(&mut h, &mut bitpos, a)).then_some(h)
}

/// Mixes two 32-bit hashes as MurmurHash3 does (`muxhash`).
fn muxhash(mut k: u32, mut h: u32) -> u32 {
    k = k.wrapping_mul(0xcc9e_2d51);
    k = k.rotate_left(15);
    k = k.wrapping_mul(0x1b87_3593);
    h ^= k;
    h = h.rotate_left(13);
    h.wrapping_mul(5).wrapping_add(0xe654_6b64)
}

fn encode_packchar(c: u8) -> u32 {
    if c < b'A' {
        27
    } else if c >= b'a' {
        u32::from(c - b'a' + 1)
    } else {
        u32::from(c - b'A' + 1)
    }
}

/// Packs up to six characters into 32 bits (`pack_chars` and `pack_chars_tail`).
fn pack_chars(s: &[u8]) -> u32 {
    s.iter().enumerate().fold(0, |h, (i, &c)| {
        h ^ (encode_packchar(c) << (1 + 5 * (5 - i as u32)))
    })
}

fn hash32(mut s: &[u8]) -> u32 {
    let mut h = 0;
    while s.len() >= 6 {
        h = muxhash(pack_chars(&s[..6]), h);
        s = &s[6..];
    }
    if !s.is_empty() {
        h = muxhash(pack_chars(s), h);
    }
    h
}

fn encode_refalt_hash(r: &[u8], a: &[u8]) -> u32 {
    let mut h = muxhash(hash32(a), muxhash(0x3, hash32(r)));
    h ^= h >> 16;
    h = h.wrapping_mul(0x85eb_ca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2_ae35);
    h ^= h >> 16;
    (h >> 1) | 1
}

fn encode_refalt(r: &[u8], a: &[u8]) -> u32 {
    if r.len() + a.len() <= 11
        && let Some(h) = encode_refalt_rev(r, a)
    {
        return h;
    }
    encode_refalt_hash(r, a)
}

/// The VariantKey of CHROM, the 0-based POS, REF, and the first ALT (`variantkey`).
pub(super) fn variantkey(chrom: &[u8], pos: u32, r: &[u8], a: &[u8]) -> u64 {
    (u64::from(encode_chrom(chrom)) << SHIFT_CHROM)
        | (u64::from(pos) << SHIFT_POS)
        | u64::from(encode_refalt(r, a))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_match_the_reference_implementation() {
        // CHROM 1 in the top five bits, REF and ALT lengths, then A=0 and C=1.
        assert_eq!(variantkey(b"1", 0, b"A", b"C"), 0x0800_0000_0888_0000);
        assert_eq!(encode_chrom(b"chrX"), 23);
        assert_eq!(encode_chrom(b"MT"), 25);
        assert_eq!(encode_chrom(b"chr22"), 22);
        assert_eq!(encode_chrom(b"HLA"), 0);
    }
}
