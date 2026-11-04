//! GFF3 gene models (bcftools' `gff.c`): genes, transcripts, and their CDS, exon, and UTR
//! features, indexed by position for `csq` and `norm -g`.
//!
//! The parsing logic follows bcftools: a gene line (`ID=gene:...`), transcript lines with a
//! supported biotype (`ID=transcript:...;Parent=gene:...`), and CDS, exon, and UTR lines
//! (`Parent=transcript:...`). Features of transcripts that are not kept are dropped. The
//! orders of transcripts and genes follow HTSlib's khash tables, so that features are indexed
//! (and listed) in bcftools' order.

use std::{
    collections::HashMap,
    io::{BufRead, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use pasta::variant::{index::data_file, readlist, regidx::RegIdx};

use crate::raw_input::{log_open_error, strerror};

/// The last two intron bases at the 5' end (bcftools' `N_SPLICE_DONOR`).
pub const N_SPLICE_DONOR: u32 = 2;
/// Intron bases of the splice region (bcftools' `N_SPLICE_REGION_INTRON`).
pub const N_SPLICE_REGION_INTRON: u32 = 8;

/// A strand (bcftools' `STRAND_*`; the discriminants are bcftools' values).
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Strand {
    Rev = 0,
    Fwd = 1,
    Unknown = 2,
}

impl Strand {
    /// bcftools' numeric value.
    #[must_use]
    pub fn code(self) -> u32 {
        self as u32
    }

    /// The GFF strand column character.
    fn gff_char(self) -> char {
        match self {
            Self::Fwd => '+',
            Self::Rev => '-',
            Self::Unknown => '.',
        }
    }
}

/// The 5' end of a transcript's CDS is incomplete (bcftools' `TRIM_5PRIME`).
pub const TRIM_5PRIME: u32 = 1;
/// The 3' end of a transcript's CDS is incomplete (bcftools' `TRIM_3PRIME`).
pub const TRIM_3PRIME: u32 = 2;

/// The bit that marks coding biotypes (bcftools' `GF_coding_bit`).
const GF_CODING_BIT: u32 = 6;

/// Biotypes and feature types (bcftools' `GF_*`).
#[expect(dead_code, reason = "the complete list of bcftools' GF_* values")]
pub mod gf {
    use super::GF_CODING_BIT;

    pub const MT_RRNA: u32 = 1;
    pub const MT_TRNA: u32 = 2;
    pub const LINCRNA: u32 = 3;
    pub const MIRNA: u32 = 4;
    pub const MISC_RNA: u32 = 5;
    pub const RRNA: u32 = 6;
    pub const SNRNA: u32 = 7;
    pub const SNORNA: u32 = 8;
    pub const PROCESSED_TRANSCRIPT: u32 = 9;
    pub const ANTISENSE: u32 = 10;
    pub const MACRO_LNCRNA: u32 = 11;
    pub const RIBOZYME_LC: u32 = 12;
    pub const SRNA: u32 = 13;
    pub const SCRNA: u32 = 14;
    pub const SCARNA: u32 = 15;
    pub const SENSE_INTRONIC: u32 = 16;
    pub const SENSE_OVERLAPPING: u32 = 17;
    pub const PSEUDOGENE: u32 = 18;
    pub const PROCESSED_PSEUDOGENE: u32 = 19;
    pub const ARTIFACT: u32 = 20;
    pub const IG_PSEUDOGENE: u32 = 21;
    pub const IG_C_PSEUDOGENE: u32 = 22;
    pub const IG_J_PSEUDOGENE: u32 = 23;
    pub const IG_V_PSEUDOGENE: u32 = 24;
    pub const TR_V_PSEUDOGENE: u32 = 25;
    pub const TR_J_PSEUDOGENE: u32 = 26;
    pub const MT_TRNA_PSEUDOGENE: u32 = 27;
    pub const MISC_RNA_PSEUDOGENE: u32 = 28;
    pub const MIRNA_PSEUDOGENE: u32 = 29;
    pub const RIBOZYME: u32 = 30;
    pub const RETAINED_INTRON: u32 = 31;
    pub const RETROTRANSPOSED: u32 = 32;
    pub const TRNA_PSEUDOGENE: u32 = 33;
    pub const TRANSCRIBED_PROCESSED_PSEUDOGENE: u32 = 34;
    pub const TRANSCRIBED_UNPROCESSED_PSEUDOGENE: u32 = 35;
    pub const TRANSCRIBED_UNITARY_PSEUDOGENE: u32 = 36;
    pub const TRANSLATED_UNPROCESSED_PSEUDOGENE: u32 = 37;
    pub const TRANSLATED_PROCESSED_PSEUDOGENE: u32 = 38;
    pub const KNOWN_NCRNA: u32 = 39;
    pub const UNITARY_PSEUDOGENE: u32 = 40;
    pub const UNPROCESSED_PSEUDOGENE: u32 = 41;
    pub const LRG_GENE: u32 = 42;
    pub const PRIME3_OVERLAPPING_NCRNA: u32 = 43;
    pub const DISRUPTED_DOMAIN: u32 = 44;
    pub const VAULTRNA: u32 = 45;
    pub const BIDIRECTIONAL_PROMOTER_LNCRNA: u32 = 46;
    pub const AMBIGUOUS_ORF: u32 = 47;
    pub const LNCRNA: u32 = 48;

    const C: u32 = 1 << GF_CODING_BIT;
    pub const PROTEIN_CODING: u32 = 1 | C;
    pub const POLYMORPHIC_PSEUDOGENE: u32 = 2 | C;
    pub const IG_C: u32 = 3 | C;
    pub const IG_D: u32 = 4 | C;
    pub const IG_J: u32 = 5 | C;
    pub const IG_LV: u32 = 6 | C;
    pub const IG_V: u32 = 7 | C;
    pub const TR_C: u32 = 8 | C;
    pub const TR_D: u32 = 9 | C;
    pub const TR_J: u32 = 10 | C;
    pub const TR_V: u32 = 11 | C;
    pub const NMD: u32 = 12 | C;
    pub const NON_STOP_DECAY: u32 = 13 | C;

    const S: u32 = 1 << (GF_CODING_BIT + 1);
    pub const CDS: u32 = S + 1;
    pub const EXON: u32 = S + 2;
    pub const UTR3: u32 = S + 3;
    pub const UTR5: u32 = S + 4;
}

/// Whether a biotype is coding (bcftools' `GF_is_coding`).
#[inline]
#[must_use]
pub fn is_coding(ty: u32) -> bool {
    ty & (1 << GF_CODING_BIT) != 0
}

const GF_STRINGS_NONCODING: [&str; 48] = [
    "MT_rRNA",
    "MT_tRNA",
    "lincRNA",
    "miRNA",
    "misc_RNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "processed_transcript",
    "antisense",
    "macro_lncRNA",
    "ribozyme",
    "sRNA",
    "scRNA",
    "scaRNA",
    "sense_intronic",
    "sense_overlapping",
    "pseudogene",
    "processed_pseudogene",
    "artifact",
    "IG_pseudogene",
    "IG_C_pseudogene",
    "IG_J_pseudogene",
    "IG_V_pseudogene",
    "TR_V_pseudogene",
    "TR_J_pseudogene",
    "MT_tRNA_pseudogene",
    "misc_RNA_pseudogene",
    "miRNA_pseudogene",
    "ribozyme",
    "retained_intron",
    "retrotransposed",
    "Trna_pseudogene",
    "transcribed_processed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "transcribed_unitary_pseudogene",
    "translated_unprocessed_pseudogene",
    "translated_processed_pseudogene",
    "known_ncRNA",
    "unitary_pseudogene",
    "unprocessed_pseudogene",
    "LRG_gene",
    "3_prime_overlapping_ncRNA",
    "disrupted_domain",
    "vaultRNA",
    "bidirectional_promoter_lncRNA",
    "ambiguous_orf",
    "lncRNA",
];
const GF_STRINGS_CODING: [&str; 13] = [
    "protein_coding",
    "polymorphic_pseudogene",
    "IG_C",
    "IG_D",
    "IG_J",
    "IG_LV",
    "IG_V",
    "TR_C",
    "TR_D",
    "TR_J",
    "TR_V",
    "NMD",
    "non_stop_decay",
];
const GF_STRINGS_SPECIAL: [&str; 4] = ["CDS", "exon", "3_prime_UTR", "5_prime_UTR"];

/// The GFF name of a biotype or feature type (bcftools' `gf_type2gff_string`).
#[must_use]
pub fn type2gff_string(ty: u32) -> &'static str {
    if !is_coding(ty) {
        if ty < (1 << GF_CODING_BIT) {
            return GF_STRINGS_NONCODING
                .get((ty as usize).wrapping_sub(1))
                .copied()
                .unwrap_or("");
        }
        let t = ty & ((1 << (GF_CODING_BIT + 1)) - 1);
        return GF_STRINGS_SPECIAL
            .get((t as usize).wrapping_sub(1))
            .copied()
            .unwrap_or("");
    }
    let t = ty & ((1 << GF_CODING_BIT) - 1);
    GF_STRINGS_CODING
        .get((t as usize).wrapping_sub(1))
        .copied()
        .unwrap_or("")
}

/// bcftools' `CDS_PHASE_UNKN`.
pub const CDS_PHASE_UNKN: u32 = 3;

/// A coding exon of a transcript (bcftools' `gf_cds_t`).
#[derive(Clone, Copy, Debug)]
pub struct Cds {
    /// The transcript (an index into [`Gff::transcripts`]).
    pub tr: u32,
    /// The first base on the reference strand (0-based).
    pub beg: u32,
    /// The offset of the first base within the spliced coding sequence.
    pub pos: u32,
    /// The length.
    pub len: u32,
    /// The index of this CDS within the transcript.
    pub icds: u32,
    /// 0, 1, 2, or [`CDS_PHASE_UNKN`].
    pub phase: u32,
}

/// A gene (bcftools' `gf_gene_t`).
#[derive(Clone, Debug, Default)]
pub struct Gene {
    /// The `Name=` attribute, or the ID; `None` for genes that only transcripts refer to.
    pub name: Option<Vec<u8>>,
    pub iseq: u32,
    id: u32,
    beg: u32,
    end: u32,
    strand: u32,
    used: bool,
}

/// An exon (bcftools' `gf_exon_t`).
#[derive(Clone, Copy, Debug)]
pub struct Exon {
    pub beg: u32,
    pub end: u32,
    pub tr: u32,
}

/// Which UTR (bcftools' `utr_t`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum UtrKind {
    Prime3,
    Prime5,
}

/// A UTR (bcftools' `gf_utr_t`).
#[derive(Clone, Copy, Debug)]
pub struct Utr {
    pub which: UtrKind,
    pub beg: u32,
    pub end: u32,
    pub tr: u32,
}

/// A transcript (bcftools' `gf_tscript_t`).
#[derive(Clone, Debug)]
pub struct Transcript {
    /// The transcript ID number ([`Gff::transcript_name`]).
    pub id: u32,
    /// The extent on the reference strand (0-based, inclusive).
    pub beg: u32,
    pub end: u32,
    pub strand: Strand,
    used: bool,
    /// The CDSs ordered by position (indexes into [`Gff::cds`]).
    pub cds: Vec<u32>,
    /// [`TRIM_5PRIME`] and [`TRIM_3PRIME`] bits.
    pub trim: u32,
    /// The biotype ([`gf`]).
    pub ty: u32,
    /// The gene (an index into [`Gff::genes`]).
    pub gene: u32,
}

/// Options of the GFF reader (bcftools' `gff_set`).
#[derive(Clone, Debug, Default)]
pub struct GffOptions {
    pub verbosity: i32,
    /// Skip transcripts whose CDS phases are inconsistent instead of failing.
    pub force: bool,
    /// Write the parsed features to this file (`--dump-gff`).
    pub dump: Option<PathBuf>,
}

/// String IDs mapped to numbers (bcftools' `id_tbl_t`).
#[derive(Default)]
struct IdTable {
    ids: HashMap<Vec<u8>, u32>,
    names: Vec<Vec<u8>>,
}

impl IdTable {
    fn register(&mut self, s: &[u8]) -> u32 {
        if let Some(&id) = self.ids.get(s) {
            return id;
        }
        let id = self.names.len() as u32;
        self.names.push(s.to_vec());
        self.ids.insert(s.to_vec(), id);
        id
    }
}

/// A key of a khash table, with HTSlib's hash function for its type.
pub(crate) trait KhashKey: Eq {
    fn khash(&self) -> u32;
}

impl KhashKey for Vec<u8> {
    /// HTSlib's `kh_str_hash_func` (`__ac_FNV1a_hash_string`).
    fn khash(&self) -> u32 {
        self.iter().fold(2_166_136_261, |h: u32, &c| {
            (h ^ u32::from(c)).wrapping_mul(16_777_619)
        })
    }
}

impl KhashKey for u32 {
    /// HTSlib's `kh_int_hash_func`.
    fn khash(&self) -> u32 {
        *self
    }
}

/// Keys with counts, kept in HTSlib's khash bucket order so that listing them gives
/// bcftools' order (a table that is only ever inserted into, like `khash_str2int` for the
/// ignored biotypes and `kh_int2tscript` for the transcripts).
pub(crate) struct KhashOrder<K> {
    /// The buckets: a key and how many times it was put.
    buckets: Vec<Option<(K, usize)>>,
    size: usize,
    upper_bound: usize,
}

impl<K> Default for KhashOrder<K> {
    fn default() -> Self {
        Self {
            buckets: Vec::new(),
            size: 0,
            upper_bound: 0,
        }
    }
}

impl<K: KhashKey> KhashOrder<K> {
    /// The bucket for `key` in a table whose length is a power of two: the one holding
    /// `key`, or the first empty one along khash's probe sequence.
    fn probe(buckets: &[Option<(K, usize)>], key: &K) -> usize {
        let mask = buckets.len() - 1;
        let mut i = key.khash() as usize & mask;
        let mut step = 0;
        while let Some((k, _)) = &buckets[i] {
            if k == key {
                break;
            }
            step += 1;
            i = (i + step) & mask;
        }
        i
    }

    /// khash's `kh_resize` growing a table without deleted entries: the same placement as
    /// its in-place kick-out rehash.
    fn grow(&mut self) {
        let old_n = self.buckets.len();
        let new_n = (old_n + 1).next_power_of_two().max(4);
        let mut old = std::mem::take(&mut self.buckets);
        let mut taken = vec![false; new_n];
        let mut out: Vec<Option<(K, usize)>> = Vec::new();
        out.resize_with(new_n, || None);
        let mask = new_n - 1;
        for j in 0..old_n {
            let Some(mut entry) = old[j].take() else {
                continue;
            };
            loop {
                let mut i = entry.0.khash() as usize & mask;
                let mut step = 0;
                while taken[i] {
                    step += 1;
                    i = (i + step) & mask;
                }
                taken[i] = true;
                // An entry not rehashed yet is kicked out to be placed next.
                let kicked = old.get_mut(i).and_then(Option::take);
                out[i] = Some(entry);
                match kicked {
                    Some(k) => entry = k,
                    None => break,
                }
            }
        }
        self.buckets = out;
        self.upper_bound = (new_n as f64 * 0.77 + 0.5) as usize;
    }

    /// Puts `key` (`kh_put`), adding one to its count.
    pub(crate) fn put(&mut self, key: K) {
        if self.size >= self.upper_bound {
            self.grow();
        }
        let i = Self::probe(&self.buckets, &key);
        match &mut self.buckets[i] {
            Some((_, n)) => *n += 1,
            slot @ None => {
                *slot = Some((key, 1));
                self.size += 1;
            }
        }
    }

    /// The keys and their counts in bucket order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = &(K, usize)> {
        self.buckets.iter().flatten()
    }

    /// The number of distinct keys.
    pub(crate) fn len(&self) -> usize {
        self.size
    }
}

/// The kind of a feature line that belongs to a transcript (bcftools' `ftr_t.type`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FeatureKind {
    Cds,
    Exon,
    Utr3,
    Utr5,
}

/// An exon, CDS, or UTR line (bcftools' `ftr_t`).
#[derive(Clone, Copy, Debug)]
struct Feature {
    kind: FeatureKind,
    beg: u32,
    end: u32,
    trid: u32,
    strand: Strand,
    phase: u32,
}

/// How many times each warning fired (bcftools' `gff_t.warned`).
#[derive(Default)]
struct Warned {
    unknown_chr: usize,
    unknown_tscript_biotype: usize,
    unknown_strand: usize,
    unknown_phase: usize,
    duplicate_id: usize,
    unknown_cds_phase: usize,
    incomplete_cds: usize,
    wrong_phase: usize,
    overlapping_cds: usize,
    ftr_out_of_bounds: usize,
}

/// A parsed GFF3 file (bcftools' `gff_t`).
pub struct Gff {
    /// The sequence names, in order of appearance.
    seqs: Vec<Vec<u8>>,
    seq2int: HashMap<Vec<u8>, u32>,
    tscript_ids: IdTable,
    pub genes: Vec<Gene>,
    pub transcripts: Vec<Transcript>,
    pub cds: Vec<Cds>,
    pub utrs: Vec<Utr>,
    pub exons: Vec<Exon>,
    /// CDSs by position (payload: an index into [`Self::cds`]).
    pub idx_cds: RegIdx<u32>,
    /// UTRs by position (payload: an index into [`Self::utrs`]).
    pub idx_utr: RegIdx<u32>,
    /// Exons extended by [`N_SPLICE_REGION_INTRON`] on both sides (payload: an index into
    /// [`Self::exons`]).
    pub idx_exon: RegIdx<u32>,
    /// Transcripts by position (payload: an index into [`Self::transcripts`]).
    pub idx_tscript: RegIdx<u32>,
}

/// The parsing state (bcftools' `aux_t` and the warning counters of `gff_t`).
struct Parser {
    opts: GffOptions,
    gff: Gff,
    features: Vec<Feature>,
    /// Gene ID number to gene (bcftools' `gid2gene`), and its khash order.
    gid2gene: HashMap<u32, u32>,
    gene_order: KhashOrder<u32>,
    /// Transcript ID number to transcript (bcftools' `id2tr`), and its khash order.
    id2tr: HashMap<u32, u32>,
    tscript_order: KhashOrder<u32>,
    ignored_biotypes: KhashOrder<Vec<u8>>,
    gene_ids: IdTable,
    warned: Warned,
}

/// The columns of a GFF line (bcftools' `gff_parse_line` locals and `aux_t` pointers).
struct Line<'a> {
    chr: &'a [u8],
    /// The third column and everything after it (bcftools compares prefixes).
    ty: &'a [u8],
    kind: Option<FeatureKind>,
    is_gene_line: bool,
    beg: u32,
    end: u32,
    strand: Strand,
    phase: u32,
    id: Option<&'a [u8]>,
    parent: Option<&'a [u8]>,
    biotype: Option<&'a [u8]>,
    name: Option<&'a [u8]>,
}

/// C's `strtol(s, &end, 10)` from `at`: the value and the index after it (`at` when nothing
/// parses).
fn strtol_at(s: &[u8], at: usize) -> (i64, usize) {
    if at >= s.len() {
        return (0, at);
    }
    let (v, n) = pasta::variant::number::strtol(&s[at..], 0);
    (v, at + n)
}

/// How a biotype prefix is compared.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Case {
    Exact,
    Ignore,
}

/// bcftools' `gff_parse_biotype` as a table: the cases of its switch on the first character,
/// in order, with the prefix each tests and the biotype it gives.
const BIOTYPES: &[(u8, &[u8], Case, u32)] = &[
    (b'p', b"protein_coding", Case::Exact, gf::PROTEIN_CODING),
    (b'p', b"pseudogene", Case::Exact, gf::PSEUDOGENE),
    (
        b'p',
        b"processed_transcript",
        Case::Exact,
        gf::PROCESSED_TRANSCRIPT,
    ),
    (
        b'p',
        b"processed_pseudogene",
        Case::Exact,
        gf::PROCESSED_PSEUDOGENE,
    ),
    (
        b'p',
        b"polymorphic_pseudogene",
        Case::Exact,
        gf::POLYMORPHIC_PSEUDOGENE,
    ),
    (b'a', b"artifact", Case::Exact, gf::ARTIFACT),
    (b'a', b"antisense", Case::Exact, gf::ANTISENSE),
    (b'a', b"ambiguous_orf", Case::Exact, gf::AMBIGUOUS_ORF),
    (b'I', b"IG_pseudogene", Case::Exact, gf::IG_PSEUDOGENE),
    (b'I', b"IG_C_pseudogene", Case::Exact, gf::IG_C_PSEUDOGENE),
    (b'I', b"IG_J_pseudogene", Case::Exact, gf::IG_J_PSEUDOGENE),
    (b'I', b"IG_V_pseudogene", Case::Exact, gf::IG_V_PSEUDOGENE),
    (b'I', b"IG_C", Case::Exact, gf::IG_C),
    (b'I', b"IG_D", Case::Exact, gf::IG_D),
    (b'I', b"IG_J", Case::Exact, gf::IG_J),
    (b'I', b"IG_V", Case::Exact, gf::IG_V),
    (b'I', b"IG_LV", Case::Exact, gf::IG_LV),
    (b'T', b"TR_V_pseudogene", Case::Exact, gf::TR_V_PSEUDOGENE),
    (b'T', b"TR_J_pseudogene", Case::Exact, gf::TR_J_PSEUDOGENE),
    (b'T', b"TR_C", Case::Exact, gf::TR_C),
    (b'T', b"TR_D", Case::Exact, gf::TR_D),
    (b'T', b"TR_J", Case::Exact, gf::TR_J),
    (b'T', b"TR_V", Case::Exact, gf::TR_V),
    (
        b'M',
        b"Mt_tRNA_pseudogene",
        Case::Exact,
        gf::MT_TRNA_PSEUDOGENE,
    ),
    (b'M', b"Mt_tRNA", Case::Ignore, gf::MT_TRNA),
    (b'M', b"Mt_rRNA", Case::Ignore, gf::MT_RRNA),
    (b'M', b"MRNA", Case::Ignore, gf::PROTEIN_CODING),
    (b'l', b"lincRNA", Case::Exact, gf::LINCRNA),
    (b'l', b"lncRNA", Case::Exact, gf::LNCRNA),
    (b'm', b"macro_lncRNA", Case::Exact, gf::MACRO_LNCRNA),
    (
        b'm',
        b"misc_RNA_pseudogene",
        Case::Exact,
        gf::MISC_RNA_PSEUDOGENE,
    ),
    (b'm', b"miRNA_pseudogene", Case::Exact, gf::MIRNA_PSEUDOGENE),
    (b'm', b"miRNA", Case::Exact, gf::MIRNA),
    (b'm', b"misc_RNA", Case::Exact, gf::MISC_RNA),
    (b'm', b"mRNA", Case::Ignore, gf::PROTEIN_CODING),
    (b'r', b"rRNA", Case::Exact, gf::RRNA),
    (b'r', b"ribozyme", Case::Exact, gf::RIBOZYME),
    (b'r', b"retained_intron", Case::Exact, gf::RETAINED_INTRON),
    (b'r', b"retrotransposed", Case::Exact, gf::RETROTRANSPOSED),
    (b's', b"snRNA", Case::Exact, gf::SNRNA),
    (b's', b"sRNA", Case::Exact, gf::SRNA),
    (b's', b"scRNA", Case::Exact, gf::SCRNA),
    (b's', b"scaRNA", Case::Exact, gf::SCARNA),
    (b's', b"snoRNA", Case::Exact, gf::SNORNA),
    (b's', b"sense_intronic", Case::Exact, gf::SENSE_INTRONIC),
    (
        b's',
        b"sense_overlapping",
        Case::Exact,
        gf::SENSE_OVERLAPPING,
    ),
    (b't', b"tRNA_pseudogene", Case::Exact, gf::TRNA_PSEUDOGENE),
    (
        b't',
        b"transcribed_processed_pseudogene",
        Case::Exact,
        gf::TRANSCRIBED_PROCESSED_PSEUDOGENE,
    ),
    (
        b't',
        b"transcribed_unprocessed_pseudogene",
        Case::Exact,
        gf::TRANSCRIBED_UNPROCESSED_PSEUDOGENE,
    ),
    (
        b't',
        b"transcribed_unitary_pseudogene",
        Case::Exact,
        gf::TRANSCRIBED_UNITARY_PSEUDOGENE,
    ),
    (
        b't',
        b"translated_unprocessed_pseudogene",
        Case::Exact,
        gf::TRANSLATED_UNPROCESSED_PSEUDOGENE,
    ),
    (
        b't',
        b"translated_processed_pseudogene",
        Case::Exact,
        gf::TRANSLATED_PROCESSED_PSEUDOGENE,
    ),
    (b'n', b"nonsense_mediated_decay", Case::Exact, gf::NMD),
    (b'n', b"non_stop_decay", Case::Exact, gf::NON_STOP_DECAY),
    (b'N', b"NMD", Case::Exact, gf::NMD),
    (b'k', b"known_ncrna", Case::Exact, gf::KNOWN_NCRNA),
    (
        b'u',
        b"unitary_pseudogene",
        Case::Exact,
        gf::UNITARY_PSEUDOGENE,
    ),
    (
        b'u',
        b"unprocessed_pseudogene",
        Case::Exact,
        gf::UNPROCESSED_PSEUDOGENE,
    ),
    (b'L', b"LRG_gene", Case::Exact, gf::LRG_GENE),
    (
        b'3',
        b"3prime_overlapping_ncRNA",
        Case::Ignore,
        gf::PRIME3_OVERLAPPING_NCRNA,
    ),
    (
        b'3',
        b"3_prime_overlapping_ncRNA",
        Case::Ignore,
        gf::PRIME3_OVERLAPPING_NCRNA,
    ),
    (b'd', b"disrupted_domain", Case::Exact, gf::DISRUPTED_DOMAIN),
    (b'v', b"vaultRNA", Case::Exact, gf::VAULTRNA),
    (
        b'b',
        b"bidirectional_promoter_lncRNA",
        Case::Exact,
        gf::BIDIRECTIONAL_PROMOTER_LNCRNA,
    ),
];

/// bcftools' `gff_parse_biotype`: the biotype of a string starting with a biotype name, or 0.
fn parse_biotype(s: &[u8]) -> u32 {
    let Some(&first) = s.first() else { return 0 };
    BIOTYPES
        .iter()
        .filter(|&&(c, ..)| c == first)
        .find(|&&(_, p, case, _)| match case {
            Case::Exact => s.starts_with(p),
            Case::Ignore => s.len() >= p.len() && s[..p.len()].eq_ignore_ascii_case(p),
        })
        .map_or(0, |&(.., ty)| ty)
}

/// The index of the next tab at or after `at`, plus one (bcftools' `gff_skip`).
fn skip_column(line: &[u8], at: usize) -> anyhow::Result<usize> {
    match line.get(at..).and_then(|s| memchr::memchr(b'\t', s)) {
        Some(i) => Ok(at + i + 1),
        None => bail!(
            "[gff.c:229 gff_skip] Could not parse the line: {}",
            String::from_utf8_lossy(line)
        ),
    }
}

/// Splits a GFF line into its columns; `None` for lines to skip (bcftools' `gff_parse_line`
/// up to the dispatch on the line type).
fn split_line(line: &[u8]) -> anyhow::Result<Option<Line<'_>>> {
    if line.is_empty() || line[0] == b'#' {
        return Ok(None);
    }
    let text = || String::from_utf8_lossy(line).into_owned();
    // 1st column: the sequence
    let Some(tab) = memchr::memchr(b'\t', line) else {
        bail!(
            "[gff.c:236 gff_parse_chr] Could not parse the line: {}",
            text()
        );
    };
    let chr = &line[..tab];
    // 2nd column: skipped
    let mut ss = skip_column(line, tab + 1)?;

    // 3rd column: the feature type
    let ty = &line[ss..];
    let mut is_gene_line = false;
    let mut kind = None;
    if ty.starts_with(b"exon\t") {
        kind = Some(FeatureKind::Exon);
        ss += 5;
    } else if ty.starts_with(b"CDS\t") {
        kind = Some(FeatureKind::Cds);
        ss += 4;
    } else if ty.starts_with(b"three_prime_UTR\t") {
        kind = Some(FeatureKind::Utr3);
        ss += 16;
    } else if ty.starts_with(b"five_prime_UTR\t") {
        kind = Some(FeatureKind::Utr5);
        ss += 15;
    } else if ty.starts_with(b"biological_region\t") {
        return Ok(None);
    } else if ty.starts_with(b"gene\t") {
        is_gene_line = true;
        ss += 5;
    } else {
        ss = skip_column(line, ss)?;
    }

    // 4th and 5th columns: beg, end
    let (beg, se) = strtol_at(line, ss);
    if se == ss {
        bail!(
            "[gff.c:244 gff_parse_beg_end] Could not parse the line:\n\t{}\n\t{}",
            text(),
            String::from_utf8_lossy(&line[ss.min(line.len())..])
        );
    }
    ss = se + 1;
    let (end, se) = strtol_at(line, ss);
    if se == ss {
        bail!(
            "[gff.c:247 gff_parse_beg_end] Could not parse the line: {}",
            text()
        );
    }
    ss = se + 1;

    // 6th column: skipped
    ss = skip_column(line, ss)?;

    // 7th column: strand
    let strand = match line.get(ss) {
        Some(b'+') => Strand::Fwd,
        Some(b'-') => Strand::Rev,
        _ => Strand::Unknown,
    };
    ss += 2;

    // 8th column: phase
    let phase = match line.get(ss) {
        Some(c @ b'0'..=b'2') => u32::from(c - b'0'),
        _ => CDS_PHASE_UNKN,
    };
    ss += 2;

    // 9th column: id, parent, name, biotype
    let (mut id, mut parent, mut biotype, mut name) = (None, None, None, None);
    let attrs = line.get(ss..).unwrap_or_default();
    for attr in attrs.split(|&c| c == b';') {
        if let Some(v) = attr.strip_prefix(b"ID=") {
            if let Some(v) = v.strip_prefix(b"gene:") {
                id = Some(v);
                is_gene_line = true;
            } else {
                id = Some(v.strip_prefix(b"transcript:").unwrap_or(v));
            }
        } else if let Some(v) = attr.strip_prefix(b"Name=") {
            name = Some(v);
        } else if let Some(v) = attr.strip_prefix(b"Parent=") {
            parent = Some(
                v.strip_prefix(b"gene:")
                    .or_else(|| v.strip_prefix(b"transcript:"))
                    .unwrap_or(v),
            );
        } else if let Some(v) = attr.strip_prefix(b"biotype=") {
            biotype = Some(v);
        } else if let Some(v) = attr.strip_prefix(b"gene_biotype=") {
            biotype = Some(v);
        }
    }
    Ok(Some(Line {
        chr,
        ty,
        kind,
        is_gene_line,
        beg: (beg - 1) as u32,
        end: (end - 1) as u32,
        strand,
        phase,
        id,
        parent,
        biotype,
        name,
    }))
}

impl Parser {
    /// Prints a warning the first time (every time with verbosity above 1).
    fn warn(&self, count: usize, msg: impl FnOnce() -> String) {
        if count == 0 || self.opts.verbosity > 1 {
            eprintln!("{}", msg());
        }
    }

    /// bcftools' `feature_set_seq`.
    fn set_seq(&mut self, chr: &[u8]) -> u32 {
        if let Some(&i) = self.gff.seq2int.get(chr) {
            return i;
        }
        let i = self.gff.seqs.len() as u32;
        self.gff.seqs.push(chr.to_vec());
        self.gff.seq2int.insert(chr.to_vec(), i);
        i
    }

    /// bcftools' `gene_init`.
    fn gene_init(&mut self, gene_id: u32) -> u32 {
        if let Some(&g) = self.gid2gene.get(&gene_id) {
            return g;
        }
        let g = self.gff.genes.len() as u32;
        self.gff.genes.push(Gene::default());
        self.gid2gene.insert(gene_id, g);
        self.gene_order.put(gene_id);
        g
    }

    /// Registers one GFF line (bcftools' `gff_parse_line`).
    fn parse_line(&mut self, raw: &[u8]) -> anyhow::Result<()> {
        let Some(line) = split_line(raw)? else {
            return Ok(());
        };
        let text = || String::from_utf8_lossy(raw).into_owned();
        let parent = match line.parent {
            Some(p) if !line.is_gene_line => p,
            _ => {
                self.parse_gene(&line, text);
                return Ok(());
            }
        };
        if let Some(kind) = line.kind {
            // bcftools' `gff_parse_exon`
            let trid = self.gff.tscript_ids.register(parent);
            if line.strand == Strand::Unknown && self.opts.verbosity > 0 {
                self.warn(self.warned.unknown_strand, || {
                    format!(
                        "Warning: Ignoring GFF feature with unknown strand .. {}",
                        text()
                    )
                });
                self.warned.unknown_strand += 1;
            }
            if line.phase == CDS_PHASE_UNKN && self.opts.verbosity > 0 {
                self.warn(self.warned.unknown_phase, || {
                    format!(
                        "Warning: Ignoring GFF feature with unknown phase .. {}",
                        text()
                    )
                });
                self.warned.unknown_phase += 1;
            }
            self.set_seq(line.chr);
            self.features.push(Feature {
                kind,
                beg: line.beg,
                end: line.end,
                trid,
                strand: line.strand,
                phase: line.phase,
            });
            return Ok(());
        }
        self.parse_transcript(&line, parent, text)
    }

    /// bcftools' `gff_parse_transcript`.
    fn parse_transcript(
        &mut self,
        line: &Line<'_>,
        parent: &[u8],
        text: impl Fn() -> String,
    ) -> anyhow::Result<()> {
        let mut ty = line.biotype.map_or(0, parse_biotype);
        if ty == 0 {
            ty = parse_biotype(line.ty);
        }
        if ty == 0 {
            if let Some(b) = line.biotype {
                self.ignored_biotypes.put(b.to_vec());
            } else if self.opts.verbosity > 0 {
                self.warn(self.warned.unknown_tscript_biotype, || {
                    format!(
                        "Warning: Ignoring transcript with unknown biotype .. {}",
                        text()
                    )
                });
                self.warned.unknown_tscript_biotype += 1;
            }
            return Ok(());
        }
        let Some(id) = line.id else {
            bail!(
                "[gff.c:441 gff_parse_transcript] Could not parse the line, neither \"ID=transcript:\" nor \"ID=\" substring is present: {}",
                text()
            );
        };
        let trid = self.gff.tscript_ids.register(id);
        let gene_id = self.gene_ids.register(parent);
        let gene = self.gene_init(gene_id);
        let tr = self.gff.transcripts.len() as u32;
        self.gff.transcripts.push(Transcript {
            id: trid,
            beg: line.beg,
            end: line.end,
            strand: line.strand,
            used: false,
            cds: Vec::new(),
            trim: 0,
            ty,
            gene,
        });
        self.id2tr.insert(trid, tr);
        self.tscript_order.put(trid);
        Ok(())
    }

    /// Registers a gene line (bcftools' `gff_parse_gene`).
    fn parse_gene(&mut self, line: &Line<'_>, text: impl Fn() -> String) {
        let Some(id) = line.id else { return };
        let gene_id = self.gene_ids.register(id);
        let g = self.gene_init(gene_id);
        if self.gff.genes[g as usize].name.is_some() {
            self.warn(self.warned.duplicate_id, || {
                format!(
                    "Warning: The GFF contains features with duplicate id .. {}",
                    text()
                )
            });
            self.warned.duplicate_id += 1;
            return;
        }
        let iseq = self.set_seq(line.chr);
        let name = match line.name {
            Some(n) => n.to_vec(),
            None => self.gene_ids.names[gene_id as usize].clone(),
        };
        let gene = &mut self.gff.genes[g as usize];
        gene.iseq = iseq;
        gene.beg = line.beg;
        gene.end = line.end;
        gene.strand = line.strand.code();
        gene.id = gene_id;
        gene.name = Some(name);
    }

    /// The sequence name of a transcript's gene (bcftools' `chr_beg_end`).
    fn tr_seq(&self, tr: u32) -> Vec<u8> {
        let g = self.gff.transcripts[tr as usize].gene;
        let iseq = self.gff.genes[g as usize].iseq;
        self.gff
            .seqs
            .get(iseq as usize)
            .cloned()
            .unwrap_or_default()
    }

    /// Attaches exons, CDSs, and UTRs to their transcripts (the feature loop of bcftools'
    /// `gff_parse`).
    fn attach_features(&mut self) -> anyhow::Result<()> {
        let features = std::mem::take(&mut self.features);
        for f in &features {
            let Some(&tr) = self.id2tr.get(&f.trid) else {
                continue;
            };
            let t = &mut self.gff.transcripts[tr as usize];
            if f.beg < t.beg || f.end > t.end {
                if self.warned.ftr_out_of_bounds == 0 || self.opts.verbosity > 1 {
                    eprintln!(
                        "Warning: The GFF contains features outside the transcript boundaries .. {}",
                        String::from_utf8_lossy(&self.gff.tscript_ids.names[t.id as usize])
                    );
                }
                self.warned.ftr_out_of_bounds += 1;
                t.beg = t.beg.min(f.beg);
                t.end = t.end.max(f.end);
            }
            t.used = true;
            let g = t.gene;
            self.gff.genes[g as usize].used = true;
            match f.kind {
                FeatureKind::Cds => {
                    // bcftools' `register_cds`
                    let t = &mut self.gff.transcripts[tr as usize];
                    if t.strand != f.strand {
                        bail!(
                            "Conflicting strand in transcript {} .. {} vs {}",
                            f.trid,
                            t.strand.code(),
                            f.strand.code()
                        );
                    }
                    let c = self.gff.cds.len() as u32;
                    self.gff.cds.push(Cds {
                        tr,
                        beg: f.beg,
                        pos: 0,
                        len: f.end.wrapping_sub(f.beg).wrapping_add(1),
                        icds: 0,
                        phase: f.phase,
                    });
                    t.cds.push(c);
                }
                FeatureKind::Exon => {
                    // bcftools' `register_exon`
                    let seq = self.tr_seq(tr);
                    let e = self.gff.exons.len() as u32;
                    self.gff.exons.push(Exon {
                        beg: f.beg,
                        end: f.end,
                        tr,
                    });
                    self.gff.idx_exon.push(
                        &seq,
                        f.beg.wrapping_sub(N_SPLICE_REGION_INTRON),
                        f.end.wrapping_add(N_SPLICE_REGION_INTRON),
                        e,
                    );
                }
                FeatureKind::Utr3 | FeatureKind::Utr5 => {
                    // bcftools' `register_utr`
                    let seq = self.tr_seq(tr);
                    let u = self.gff.utrs.len() as u32;
                    self.gff.utrs.push(Utr {
                        which: if f.kind == FeatureKind::Utr3 {
                            UtrKind::Prime3
                        } else {
                            UtrKind::Prime5
                        },
                        beg: f.beg,
                        end: f.end,
                        tr,
                    });
                    self.gff.idx_utr.push(&seq, f.beg, f.end, u);
                }
            }
        }
        Ok(())
    }

    fn tr_name(&self, tr: u32) -> String {
        let id = self.gff.transcripts[tr as usize].id;
        String::from_utf8_lossy(&self.gff.tscript_ids.names[id as usize]).into_owned()
    }

    /// Checks a transcript's CDS phases, walking `order` (bcftools' sanity check loop of
    /// `tscript_init_cds`). Returns `false` when the transcript is skipped with `--force`.
    fn check_phases(
        &mut self,
        tr: u32,
        order: impl Iterator<Item = usize>,
    ) -> anyhow::Result<bool> {
        let mut len: i32 = 0;
        for i in order {
            let c = self.gff.cds[self.gff.transcripts[tr as usize].cds[i] as usize];
            if c.phase == CDS_PHASE_UNKN {
                if self.opts.verbosity > 0 {
                    self.warn(self.warned.unknown_cds_phase, || {
                        format!(
                            "Warning: CDS with unknown phase, could not verify reading frame in transcript {}",
                            self.tr_name(tr)
                        )
                    });
                    self.warned.unknown_cds_phase += 1;
                }
                len = len.wrapping_add(c.len as i32);
                continue;
            }
            let phase = if c.phase > 0 { 3 - c.phase as i32 } else { 0 };
            if phase != len % 3 {
                if !self.opts.force {
                    bail!(
                        "Error: GFF3 assumption failed for transcript {}, CDS={}: phase!=len%3 (phase={phase}, len={len}). Use the --force option to proceed anyway (at your own risk).",
                        self.tr_name(tr),
                        c.beg.wrapping_add(1)
                    );
                }
                if self.opts.verbosity > 0 {
                    self.warn(self.warned.wrong_phase, || {
                        format!(
                            "Warning: The GFF has inconsistent phase column in transcript {}, skipping. CDS pos={}: phase!=len%3 (phase={phase}, len={len})",
                            self.tr_name(tr),
                            c.beg.wrapping_add(1)
                        )
                    });
                    self.warned.wrong_phase += 1;
                }
                return Ok(false);
            }
            len = len.wrapping_add(c.len as i32);
        }
        Ok(true)
    }

    /// Sorts, trims, and checks the CDSs of every transcript and indexes transcripts and
    /// CDSs (bcftools' `tscript_init_cds`).
    fn init_cds(&mut self) -> anyhow::Result<()> {
        let order: Vec<u32> = self.tscript_order.iter().map(|&(trid, _)| trid).collect();
        for trid in order {
            let tr = self.id2tr[&trid];
            let seq = self.tr_seq(tr);
            let (tbeg, tend) = {
                let t = &self.gff.transcripts[tr as usize];
                (t.beg, t.end)
            };
            self.gff.idx_tscript.push(&seq, tbeg, tend, tr);
            if self.gff.transcripts[tr as usize].cds.is_empty() {
                continue;
            }
            let mut cds = std::mem::take(&mut self.gff.transcripts[tr as usize].cds);
            cds.sort_by_key(|&c| self.gff.cds[c as usize].beg);
            self.gff.transcripts[tr as usize].cds = cds;
            let ncds = self.gff.transcripts[tr as usize].cds.len();
            let ok = match self.gff.transcripts[tr as usize].strand {
                Strand::Fwd => {
                    let c0 = self.gff.transcripts[tr as usize].cds[0] as usize;
                    let c = &mut self.gff.cds[c0];
                    if c.phase != CDS_PHASE_UNKN {
                        if c.phase != 0 {
                            self.gff.transcripts[tr as usize].trim |= TRIM_5PRIME;
                        }
                        c.beg = c.beg.wrapping_add(c.phase);
                        c.len = c.len.wrapping_sub(c.phase);
                        c.phase = 0;
                    }
                    self.check_phases(tr, 0..ncds)?
                }
                Strand::Rev => {
                    self.trim_rev_5prime(tr);
                    self.check_phases(tr, (0..ncds).rev())?
                }
                Strand::Unknown => continue,
            };
            if !ok {
                continue;
            }
            self.finish_cds(tr, &seq);
        }
        Ok(())
    }

    /// Trims the 5' end of a reverse-strand transcript by the phase of its last CDS.
    fn trim_rev_5prime(&mut self, tr: u32) {
        let t = &self.gff.transcripts[tr as usize];
        let last = t.cds.len() - 1;
        let cds_ids = t.cds.clone();
        if self.gff.cds[cds_ids[last] as usize].phase == CDS_PHASE_UNKN {
            return;
        }
        // The phase can be bigger than the CDS (Mus_musculus.GRCm38.85 ENSMUST00000163141).
        let mut i = last as isize;
        let mut phase = self.gff.cds[cds_ids[last] as usize].phase;
        if phase != 0 {
            self.gff.transcripts[tr as usize].trim |= TRIM_5PRIME;
        }
        while i >= 0 && phase > self.gff.cds[cds_ids[i as usize] as usize].len {
            let c = &mut self.gff.cds[cds_ids[i as usize] as usize];
            phase -= c.len;
            c.phase = 0;
            c.len = 0;
            i -= 1;
        }
        // bcftools reads before the first CDS here when all were consumed.
        let Some(&k) = usize::try_from(i).ok().and_then(|i| cds_ids.get(i)) else {
            return;
        };
        if self.opts.verbosity > 0 && self.gff.cds[k as usize].phase != 0 {
            self.warn(self.warned.incomplete_cds, || truncated(&self.tr_name(tr)));
            self.warned.incomplete_cds += 1;
        }
        let c = &mut self.gff.cds[k as usize];
        c.len = c.len.wrapping_sub(c.phase);
        c.phase = 0;
    }

    /// Numbers the CDSs, warns about overlaps, trims an incomplete 3' end, sets the offsets,
    /// and indexes the CDSs (the end of bcftools' `tscript_init_cds`).
    fn finish_cds(&mut self, tr: u32, seq: &[u8]) {
        let cds_ids = self.gff.transcripts[tr as usize].cds.clone();
        let ncds = cds_ids.len();
        let mut len: i32 = 0;
        for i in 0..ncds {
            self.gff.cds[cds_ids[i] as usize].icds = i as u32;
            len = len.wrapping_add(self.gff.cds[cds_ids[i] as usize].len as i32);
            if i == 0 {
                continue;
            }
            let a = self.gff.cds[cds_ids[i - 1] as usize];
            let b = self.gff.cds[cds_ids[i] as usize];
            if a.beg.wrapping_add(a.len).wrapping_sub(1) >= b.beg && self.opts.verbosity > 0 {
                self.warn(self.warned.overlapping_cds, || {
                    format!(
                        "Warning: GFF contains overlapping CDS {}, {}-{} and {}-{} (ribosomal slippage?)",
                        self.tr_name(tr),
                        a.beg.wrapping_add(1),
                        a.beg.wrapping_add(a.len),
                        b.beg.wrapping_add(1),
                        b.beg.wrapping_add(b.len)
                    )
                });
                self.warned.overlapping_cds += 1;
            }
        }
        if len % 3 != 0 {
            // Very common (13k transcripts in GRCh37 with an incomplete 3' CDS).
            if self.opts.verbosity > 0 {
                self.warn(self.warned.incomplete_cds, || truncated(&self.tr_name(tr)));
                self.warned.incomplete_cds += 1;
            }
            self.gff.transcripts[tr as usize].trim |= TRIM_3PRIME;
            match self.gff.transcripts[tr as usize].strand {
                Strand::Fwd => {
                    let mut i = ncds as isize - 1;
                    while i >= 0 && len % 3 != 0 {
                        let c = &mut self.gff.cds[cds_ids[i as usize] as usize];
                        let r = (len % 3) as u32;
                        let dlen = if c.len >= r { r } else { c.len };
                        c.len -= dlen;
                        len -= dlen as i32;
                        i -= 1;
                    }
                }
                Strand::Rev => {
                    let mut i = 0;
                    while i < ncds && len % 3 != 0 {
                        let c = &mut self.gff.cds[cds_ids[i] as usize];
                        let r = (len % 3) as u32;
                        let dlen = if c.len >= r { r } else { c.len };
                        c.len -= dlen;
                        c.beg = c.beg.wrapping_add(dlen);
                        len -= dlen as i32;
                        i += 1;
                    }
                }
                Strand::Unknown => {}
            }
        }
        let mut len: u32 = 0;
        for &c in &cds_ids {
            let cds = &mut self.gff.cds[c as usize];
            cds.pos = len;
            len = len.wrapping_add(cds.len);
            let (beg, end) = (cds.beg, cds.beg.wrapping_add(cds.len).wrapping_sub(1));
            self.gff.idx_cds.push(seq, beg, end, c);
        }
    }

    /// Writes the parsed features (bcftools' `gff_dump`).
    fn dump(&mut self, path: &Path) -> anyhow::Result<()> {
        let file = std::fs::File::create(path)
            .map_err(|e| anyhow::anyhow!("Failed to open {}: {}", path.display(), strerror(&e)))?;
        let mut out = flate2::write::GzEncoder::new(
            std::io::BufWriter::new(file),
            flate2::Compression::default(),
        );
        let werr = || format!("Error writing {}", path.display());
        let strand_char = |s: Strand| s.gff_char();
        let gff = &mut self.gff;
        for &(gid, _) in self.gene_order.iter() {
            let gene = &gff.genes[self.gid2gene[&gid] as usize];
            let seq = gff
                .seqs
                .get(gene.iseq as usize)
                .map_or("", |s| std::str::from_utf8(s).unwrap_or(""));
            let strand = match gene.strand {
                1 => '+',
                0 => '-',
                _ => '.',
            };
            let name = gene
                .name
                .as_deref()
                .map_or_else(|| "(null)".into(), String::from_utf8_lossy);
            writeln!(
                out,
                "{seq}\t.\tgene\t{}\t{}\t.\t{strand}\t.\tID={};Name={name};used={}",
                gene.beg.wrapping_add(1),
                gene.end.wrapping_add(1),
                String::from_utf8_lossy(&self.gene_ids.names[gene.id as usize]),
                i32::from(gene.used)
            )
            .with_context(werr)?;
        }
        let names = &gff.tscript_ids.names;
        for r in gff.idx_tscript.iter() {
            let tr = &gff.transcripts[*r.payload as usize];
            let gene = &gff.genes[tr.gene as usize];
            let ty = if tr.ty == gf::PROTEIN_CODING {
                "mRNA"
            } else {
                type2gff_string(tr.ty)
            };
            writeln!(
                out,
                "{}\t.\t{ty}\t{}\t{}\t.\t{}\t.\tID={};Parent={};biotype={};used={}",
                String::from_utf8_lossy(r.seq),
                r.beg.wrapping_add(1),
                r.end.wrapping_add(1),
                strand_char(tr.strand),
                String::from_utf8_lossy(&names[tr.id as usize]),
                String::from_utf8_lossy(&self.gene_ids.names[gene.id as usize]),
                type2gff_string(tr.ty),
                i32::from(tr.used)
            )
            .with_context(werr)?;
        }
        for r in gff.idx_cds.iter() {
            let cds = &gff.cds[*r.payload as usize];
            let tr = &gff.transcripts[cds.tr as usize];
            let phase = if cds.phase == 3 {
                '.'
            } else {
                char::from(b'0' + cds.phase as u8)
            };
            writeln!(
                out,
                "{}\t.\tCDS\t{}\t{}\t.\t{}\t{phase}\tParent={}",
                String::from_utf8_lossy(r.seq),
                cds.beg.wrapping_add(1),
                cds.beg.wrapping_add(cds.len),
                strand_char(tr.strand),
                String::from_utf8_lossy(&names[tr.id as usize])
            )
            .with_context(werr)?;
        }
        for r in gff.idx_utr.iter() {
            let utr = &gff.utrs[*r.payload as usize];
            let tr = &gff.transcripts[utr.tr as usize];
            writeln!(
                out,
                "{}\t.\t{}_prime_UTR\t{}\t{}\t.\t{}\t.\tParent={}",
                String::from_utf8_lossy(r.seq),
                if utr.which == UtrKind::Prime3 {
                    "three"
                } else {
                    "five"
                },
                utr.beg.wrapping_add(1),
                utr.end.wrapping_add(1),
                strand_char(tr.strand),
                String::from_utf8_lossy(&names[tr.id as usize])
            )
            .with_context(werr)?;
        }
        for r in gff.idx_exon.iter() {
            let exon = &gff.exons[*r.payload as usize];
            let tr = &gff.transcripts[exon.tr as usize];
            writeln!(
                out,
                "{}\t.\texon\t{}\t{}\t.\t{}\t.\tParent={}",
                String::from_utf8_lossy(r.seq),
                exon.beg.wrapping_add(1),
                exon.end.wrapping_add(1),
                strand_char(tr.strand),
                String::from_utf8_lossy(&names[tr.id as usize])
            )
            .with_context(werr)?;
        }
        out.finish()
            .and_then(|mut w| w.flush())
            .with_context(|| format!("Error: close failed .. {}", path.display()))?;
        Ok(())
    }
}

/// bcftools' note about a transcript whose CDS is not a whole number of codons.
fn truncated(name: &str) -> String {
    format!("Note: truncated transcript {name} with incomplete CDS (this is very common)")
}

impl Gff {
    /// Reads and indexes a GFF3 file, plain or compressed (bcftools' `gff_parse`).
    /// `FILE##idx##INDEX` reads `FILE`, as `hts_open` does.
    ///
    /// # Errors
    ///
    /// Returns bcftools' errors for unreadable or unparsable files, conflicting strands,
    /// inconsistent CDS phases (without `force`), and files without usable transcripts.
    pub fn parse(path: &Path, opts: GffOptions) -> anyhow::Result<Self> {
        if opts.verbosity > 0 {
            eprintln!("Parsing {} ...", path.display());
        }
        let fail = || format!("Failed to read {}", path.display());
        let data = data_file(path);
        let mut reader = readlist::open_text(data).map_err(|e| {
            log_open_error(&data.display().to_string(), &e);
            anyhow::anyhow!(fail())
        })?;
        let mut p = Parser {
            opts,
            gff: Gff {
                seqs: Vec::new(),
                seq2int: HashMap::new(),
                tscript_ids: IdTable::default(),
                genes: Vec::new(),
                transcripts: Vec::new(),
                cds: Vec::new(),
                utrs: Vec::new(),
                exons: Vec::new(),
                idx_cds: RegIdx::new(),
                idx_utr: RegIdx::new(),
                idx_exon: RegIdx::new(),
                idx_tscript: RegIdx::new(),
            },
            features: Vec::new(),
            gid2gene: HashMap::new(),
            gene_order: KhashOrder::default(),
            id2tr: HashMap::new(),
            tscript_order: KhashOrder::default(),
            ignored_biotypes: KhashOrder::default(),
            gene_ids: IdTable::default(),
            warned: Warned::default(),
        };
        let mut line = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line).with_context(fail)? == 0 {
                break;
            }
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            if line.last() == Some(&b'\r') {
                line.pop();
            }
            let end = memchr::memchr(0, &line).unwrap_or(line.len());
            p.parse_line(&line[..end])?;
        }
        p.attach_features()?;
        p.init_cds()?;
        let verbosity = p.opts.verbosity;
        if verbosity > 0 {
            eprintln!(
                "Indexed {} transcripts, {} exons, {} CDSs, {} UTRs",
                p.gff.idx_tscript.region_count(),
                p.gff.idx_exon.region_count(),
                p.gff.idx_cds.region_count(),
                p.gff.idx_utr.region_count()
            );
        }
        if verbosity > 0 && p.ignored_biotypes.len() > 0 {
            eprintln!("Ignored the following biotypes:");
            for (biotype, n) in p.ignored_biotypes.iter() {
                let b = String::from_utf8_lossy(biotype);
                if b == "TCE" {
                    eprintln!("\t{n}x\t.. TCE (\"To be Experimentally Confirmed\")");
                } else {
                    eprintln!("\t{n}x\t.. {b}");
                }
            }
        }
        if verbosity > 0 {
            let w = &p.warned;
            let suppressed: usize = [
                w.unknown_chr,
                w.unknown_tscript_biotype,
                w.unknown_strand,
                w.unknown_phase,
                w.duplicate_id,
                w.unknown_cds_phase,
                w.incomplete_cds,
                w.wrong_phase,
                w.overlapping_cds,
            ]
            .iter()
            .map(|&n| {
                if verbosity > 1 {
                    0
                } else {
                    n.saturating_sub(1)
                }
            })
            .sum();
            if suppressed > 0 {
                eprintln!(
                    "Warning: {suppressed} warnings were suppressed, increase verbosity to see them all"
                );
            }
        }
        if let Some(dump) = p.opts.dump.take() {
            p.dump(&dump)?;
        }
        if p.gff.idx_tscript.region_count() == 0 {
            bail!(
                "Error: No usable transcripts found, likely a failure to parse a non-standard GFF file. Please check if the misc/gff2gff\n       or misc/gff2gff.py script can fix the problem (both do different things). See also the man page for the description\n       of the expected format http://samtools.github.io/bcftools/bcftools-man.html#csq"
            );
        }
        Ok(p.gff)
    }

    /// Whether features were read for a sequence (bcftools' `gff_has_seq`).
    #[must_use]
    pub fn has_seq(&self, seq: &[u8]) -> bool {
        self.seq2int.contains_key(seq)
    }

    /// The name of the `i`-th sequence (bcftools' `gff_iseq`).
    #[must_use]
    pub fn seq_name(&self, i: usize) -> Option<&[u8]> {
        self.seqs.get(i).map(Vec::as_slice)
    }

    /// The name of a transcript ID (bcftools' `gff_id2string`).
    #[must_use]
    pub fn transcript_name(&self, id: u32) -> &[u8] {
        &self.tscript_ids.names[id as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order<K: KhashKey + Clone>(keys: &[K]) -> Vec<(K, usize)> {
        let mut t = KhashOrder::default();
        for k in keys {
            t.put(k.clone());
        }
        t.iter().cloned().collect()
    }

    #[test]
    fn khash_order_matches_bcftools_listing() {
        // bcftools 1.24 lists these ignored biotypes as foo_rna, bar, TCE.
        let keys: Vec<Vec<u8>> = ["TCE", "foo_rna", "foo_rna", "bar"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect();
        let got: Vec<(String, usize)> = order(&keys)
            .into_iter()
            .map(|(k, n)| (String::from_utf8(k).expect("ASCII keys"), n))
            .collect();
        assert_eq!(
            got,
            [("foo_rna".into(), 2), ("bar".into(), 1), ("TCE".into(), 1)]
        );
    }

    #[test]
    fn khash_int_order_wraps_around_the_buckets() {
        // Integer keys hash to themselves: in 4 buckets 6 lands in bucket 2, before 3.
        let got: Vec<u32> = order(&[3u32, 6]).into_iter().map(|(k, _)| k).collect();
        assert_eq!(got, [6, 3]);
        // Growing rehashes everything: 4 buckets become 8 at the fourth key, and the kicked
        // out 7 moves to bucket 7.
        let got: Vec<u32> = order(&[7u32, 3, 1, 6])
            .into_iter()
            .map(|(k, _)| k)
            .collect();
        assert_eq!(got, [1, 3, 6, 7]);
    }

    #[test]
    fn biotypes_match_by_prefix() {
        assert_eq!(parse_biotype(b"protein_coding"), gf::PROTEIN_CODING);
        assert_eq!(parse_biotype(b"mRNA\t1\t2"), gf::PROTEIN_CODING);
        assert_eq!(parse_biotype(b"IG_LV_gene"), gf::IG_LV);
        assert_eq!(parse_biotype(b"lnc_RNA"), 0);
        assert_eq!(type2gff_string(gf::PROTEIN_CODING), "protein_coding");
        assert_eq!(type2gff_string(gf::LNCRNA), "lncRNA");
        assert_eq!(type2gff_string(gf::CDS), "CDS");
    }
}
