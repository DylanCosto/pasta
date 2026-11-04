//! Translation of spliced coding sequences (bcftools' `cds_translate`), and the `kstring_t`
//! buffers it fills.

use anyhow::bail;

use super::{N_REF_PAD, gencode::Gencode};
use crate::gff::Strand;

/// A `kstring_t`: bytes and a length. bcftools truncates these strings by writing a NUL and
/// sometimes reads past the length, so the bytes beyond it stay (reads past the buffer give
/// 0).
#[derive(Clone, Debug, Default)]
pub(super) struct KStr {
    buf: Vec<u8>,
    len: usize,
}

impl KStr {
    /// The length (`l`).
    #[inline]
    pub(super) fn len(&self) -> usize {
        self.len
    }

    /// Sets the length to zero, keeping the bytes (`str.l = 0`).
    #[inline]
    pub(super) fn clear(&mut self) {
        self.len = 0;
    }

    /// Whether nothing was ever written (C's NULL `s`: `kputsn` allocates even for no bytes).
    #[inline]
    pub(super) fn is_null(&self) -> bool {
        self.buf.is_empty()
    }

    /// The byte at `i`, also past the length; 0 past the buffer.
    #[inline]
    pub(super) fn at(&self, i: usize) -> u8 {
        self.buf.get(i).copied().unwrap_or(0)
    }

    /// The bytes up to the length.
    #[inline]
    pub(super) fn as_bytes(&self) -> &[u8] {
        &self.buf[..self.len.min(self.buf.len())]
    }

    #[inline]
    fn set(&mut self, i: usize, c: u8) {
        if i < self.buf.len() {
            self.buf[i] = c;
        } else {
            self.buf.resize(i, 0);
            self.buf.push(c);
        }
    }

    /// Appends a byte without terminating the string (`kputc_`).
    #[inline]
    pub(super) fn push(&mut self, c: u8) {
        self.set(self.len, c);
        self.len += 1;
    }

    /// Writes a NUL after the string (`kputc_(0, s); s->l--`).
    #[inline]
    pub(super) fn terminate(&mut self) {
        self.set(self.len, 0);
    }

    /// Appends bytes and a NUL (`kputsn`).
    pub(super) fn put(&mut self, bytes: impl IntoIterator<Item = u8>) {
        for c in bytes {
            self.push(c);
        }
        self.terminate();
    }

    /// Sets the length, keeping the bytes (`s->l -= n`).
    #[inline]
    pub(super) fn set_len(&mut self, n: usize) {
        self.len = n;
    }

    /// Cuts the string to `n` bytes with a NUL (`s->l = n; s->s[n] = 0`).
    #[inline]
    pub(super) fn truncate(&mut self, n: usize) {
        self.len = n;
        self.set(n, 0);
    }
}

/// C's `strncmp(a + ao, b + bo, n) == 0` on two buffers.
pub(super) fn strncmp_eq(a: &KStr, ao: usize, b: &KStr, bo: usize, n: usize) -> bool {
    for i in 0..n {
        let (x, y) = (a.at(ao + i), b.at(bo + i));
        if x != y {
            return false;
        }
        if x == 0 {
            return true;
        }
    }
    true
}

/// The byte at `i` of `s`, 0 outside it.
#[inline]
fn byte(s: &[u8], i: isize) -> u8 {
    usize::try_from(i)
        .ok()
        .and_then(|i| s.get(i))
        .copied()
        .unwrap_or(0)
}

/// A codon starting at `i` of `s`.
#[inline]
fn codon_at(s: &[u8], i: isize) -> [u8; 3] {
    [byte(s, i), byte(s, i + 1), byte(s, i + 2)]
}

/// Where a translation stops (bcftools' `fill` argument).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Fill {
    /// Only the codons of the sequence.
    No,
    /// A frameshift: continue to the end (forward strand) or the start (reverse strand) of
    /// the transcript.
    ToEnd,
}

/// The part of a spliced sequence to translate (bcftools' `seq` argument of
/// `cds_translate`).
pub(super) struct SplicedSeq<'a> {
    /// The bases (`seq.s`, `seq.l`).
    pub bases: &'a [u8],
    /// The length of the whole spliced sequence (`seq.m`).
    pub total: i64,
}

/// Translates a part of a spliced transcript, padding it with reference bases to whole
/// codons (bcftools' `cds_translate`).
///
/// `sref` is the spliced reference with [`N_REF_PAD`] bases on both ends; `sbeg` the offset
/// of `seq` in the spliced query; `rbeg` and `rend` the reference range `seq` replaces.
#[expect(
    clippy::too_many_arguments,
    reason = "line-for-line port of cds_translate"
)]
pub(super) fn cds_translate(
    gc: &Gencode,
    sref: &[u8],
    seq: &SplicedSeq<'_>,
    sbeg: u32,
    rbeg: u32,
    rend: u32,
    strand: Strand,
    tseq: &mut KStr,
    tseq_stop: &mut KStr,
    fill: Fill,
) -> anyhow::Result<()> {
    tseq.clear();
    tseq_stop.clear();
    let s = seq.bases;
    if s.is_empty() {
        tseq.push(b'?');
        tseq.terminate();
        tseq_stop.push(b'?');
        tseq_stop.terminate();
        return Ok(());
    }
    let pad = N_REF_PAD as isize;
    let mut tmp = [0u8; 3];
    match strand {
        Strand::Fwd => {
            // left padding
            // bcftools asserts npad<=rbeg; release builds go on, and so does this port, as
            // its reads are bounds-checked.
            let npad = (sbeg % 3) as isize;
            let mut i: isize = 0;
            while i < npad {
                tmp[i as usize] = byte(sref, rbeg as isize + i - npad + pad);
                i += 1;
            }
            while i < 3 && i - npad < s.len() as isize {
                tmp[i as usize] = s[(i - npad) as usize];
                i += 1;
            }
            let len = s.len() as isize - i + npad;
            if i == 3 {
                tseq.push(gc.dna2aa(tmp));
                tseq_stop.push(gc.dna2stop(tmp));
                let mut codon = 3 - npad;
                let end = codon + len - 1 - (len % 3);
                while codon < end {
                    let c = codon_at(s, codon);
                    tseq.push(gc.dna2aa(c));
                    tseq_stop.push(gc.dna2stop(c));
                    codon += 3;
                }
                let end = s.len() as isize - 1;
                i = 0;
                while codon + i <= end {
                    tmp[i as usize] = byte(s, codon + i);
                    i += 1;
                }
            }
            // right padding
            let mut codon = rend as isize + pad;
            if i > 0 {
                while i < 3 {
                    tmp[i as usize] = byte(sref, codon);
                    codon += 1;
                    i += 1;
                }
                tseq.push(gc.dna2aa(tmp));
                tseq_stop.push(gc.dna2stop(tmp));
            }
            if fill == Fill::ToEnd {
                let end = sref.len() as isize - pad;
                while codon + 3 <= end {
                    let c = codon_at(sref, codon);
                    tseq.push(gc.dna2aa(c));
                    tseq_stop.push(gc.dna2stop(c));
                    codon += 3;
                }
            }
        }
        Strand::Rev => {
            // right padding: the number of bases to take from the reference
            // In C this is size_t arithmetic, which wraps when the sequence is shorter than
            // expected; bcftools then prints a message and asserts, but its release builds
            // go on, as does this port (its reads are bounds-checked).
            let total = seq.total as u64;
            let used = u64::from(sbeg) + s.len() as u64;
            let npad = (total.wrapping_sub(used) % 3) as i64;
            if used + npad as u64 > total {
                eprintln!(
                    "sbeg={sbeg}  seq.l={} seq.m={} npad={npad}",
                    s.len() as i32,
                    seq.total as i32
                );
            }
            let mut i: isize = match npad {
                2 => {
                    tmp[1] = byte(sref, rend as isize + pad);
                    tmp[2] = byte(sref, rend as isize + pad + 1);
                    0
                }
                1 => {
                    tmp[2] = byte(sref, rend as isize + pad);
                    1
                }
                _ => 2,
            };
            let mut end = s.len() as isize;
            while i >= 0 && end > 0 {
                end -= 1;
                tmp[i as usize] = s[end as usize];
                i -= 1;
            }
            if i == -1 {
                tseq.push(gc.cdna2aa(tmp));
                tseq_stop.push(gc.cdna2stop(tmp));
                let mut codon = end - 3;
                while codon >= 0 {
                    let c = codon_at(s, codon);
                    tseq.push(gc.cdna2aa(c));
                    tseq_stop.push(gc.cdna2stop(c));
                    codon -= 3;
                }
                if codon == -2 {
                    tmp[2] = s[0];
                    i = 1;
                } else if codon == -1 {
                    tmp[1] = s[0];
                    tmp[2] = byte(s, 1);
                    i = 0;
                } else {
                    i = -1;
                }
            }
            // left padding
            let mut end = pad + rbeg as isize;
            if i >= 0 {
                while i >= 0 && end >= 0 {
                    end -= 1;
                    tmp[i as usize] = byte(sref, end);
                    i -= 1;
                }
                tseq.push(gc.cdna2aa(tmp));
                tseq_stop.push(gc.cdna2stop(tmp));
            }
            if fill == Fill::ToEnd {
                let mut codon = end - 3;
                while codon >= pad {
                    let c = codon_at(sref, codon);
                    tseq.push(gc.cdna2aa(c));
                    tseq_stop.push(gc.cdna2stop(c));
                    codon -= 3;
                }
            }
        }
        Strand::Unknown => bail!("Should not happen: {}", strand.code()),
    }
    tseq.terminate();
    tseq_stop.terminate();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::csq::gencode::{GencodeChoice, choose};

    fn standard() -> &'static Gencode {
        match choose(None).expect("table 0 exists") {
            GencodeChoice::Table(g) => g,
            GencodeChoice::List => panic!("table 0 is not a listing"),
        }
    }

    /// Pads a spliced sequence with `N_REF_PAD` bases on both ends.
    fn padded(seq: &[u8]) -> Vec<u8> {
        let mut s = vec![b'N'; N_REF_PAD];
        s.extend_from_slice(seq);
        s.extend(std::iter::repeat_n(b'N', N_REF_PAD));
        s
    }

    fn translate(
        sref: &[u8],
        seq: &[u8],
        total: i64,
        (sbeg, rbeg, rend): (u32, u32, u32),
        strand: Strand,
    ) -> (Vec<u8>, Vec<u8>) {
        let (mut aa, mut stop) = (KStr::default(), KStr::default());
        cds_translate(
            standard(),
            sref,
            &SplicedSeq { bases: seq, total },
            sbeg,
            rbeg,
            rend,
            strand,
            &mut aa,
            &mut stop,
            Fill::No,
        )
        .expect("the translation succeeds");
        (aa.as_bytes().to_vec(), stop.as_bytes().to_vec())
    }

    #[test]
    fn translates_whole_codons_on_both_strands() {
        let cds = b"ATGAAATAG";
        let sref = padded(cds);
        let (aa, stop) = translate(&sref, cds, 9, (0, 0, 9), Strand::Fwd);
        assert_eq!((aa.as_slice(), stop.as_slice()), (&b"MK*"[..], &b"M-*"[..]));
        // The same protein on the reverse strand: CTA TTT CAT read backwards.
        let rc = b"CTATTTCAT";
        let sref = padded(rc);
        let (aa, _) = translate(&sref, rc, 9, (0, 0, 9), Strand::Rev);
        assert_eq!(aa, b"MK*");
    }

    #[test]
    fn pads_partial_codons_with_reference_bases() {
        // A substitution of the second codon's middle base: AAA>AGA (K>R).
        let sref = padded(b"ATGAAATAG");
        let (aa, _) = translate(&sref, b"G", 9, (4, 4, 5), Strand::Fwd);
        assert_eq!(aa, b"R");
    }

    #[test]
    fn empty_sequences_translate_to_unknown() {
        let sref = padded(b"ATG");
        let (aa, stop) = translate(&sref, b"", 3, (0, 0, 0), Strand::Fwd);
        assert_eq!((aa.as_slice(), stop.as_slice()), (&b"?"[..], &b"?"[..]));
    }

    #[test]
    fn reverse_padding_wraps_like_size_t() {
        // A sequence ending past the spliced total (as with a duplicated CDS): C computes
        // (total - (sbeg + len)) % 3 in size_t, (2^64 - 21) % 3 = 1, and goes on.
        let sref = padded(b"CTATTTCATCTATTTCAT");
        let (aa, _) = translate(&sref, b"A", 59, (79, 3, 4), Strand::Rev);
        let (aa1, _) = translate(&sref, b"A", 5, (3, 3, 4), Strand::Rev);
        // npad 1 either way: the base after the sequence completes the codon.
        assert_eq!(aa, aa1);
    }

    #[test]
    fn kstr_keeps_bytes_past_its_length() {
        let mut k = KStr::default();
        k.put(*b"MKV");
        k.truncate(1);
        assert_eq!(k.as_bytes(), b"M");
        // bcftools compares past the length: the NUL, then what an earlier string left.
        assert_eq!((k.at(1), k.at(2), k.at(99)), (0, b'V', 0));
        assert!(strncmp_eq(&k, 0, &k, 0, 8));
    }

    #[test]
    fn kstr_is_null_until_written() {
        let mut k = KStr::default();
        assert!(k.is_null());
        // kputsn allocates even for no bytes.
        k.put([]);
        assert!(!k.is_null());
        assert_eq!(k.len(), 0);
        k.clear();
        assert!(
            !k.is_null(),
            "clearing keeps the buffer, as `s->l = 0` does"
        );
    }
}
