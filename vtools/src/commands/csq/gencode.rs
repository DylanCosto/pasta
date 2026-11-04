//! Genetic code tables and codon translation (bcftools' `gencode_tables`, `dna2aa`,
//! `cdna2aa`, `dna2stop`, `cdna2stop`, and `init_gencode`).

use std::io::Write;

use anyhow::bail;

/// A genetic code table (bcftools' `gencode_t`): amino acids and start/stop marks of the 64
/// codons ordered AAA, AAC, ..., TTT.
#[derive(Debug)]
pub(super) struct Gencode {
    id: i32,
    name: &'static str,
    code: [u8; 64],
    stop: [u8; 64],
}

/// The tables, generated with bcftools' `misc/gencode-tables`.
#[rustfmt::skip]
static GENCODE_TABLES: [Gencode; 27] = [
    Gencode { id: 0, name: "Standard simplified", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", stop: *b"--------------M---------------------------------*-*-----*-------" },
    Gencode { id: 1, name: "Standard", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", stop: *b"--------------M---------------M-----------------*-*-----*-----M-" },
    Gencode { id: 2, name: "Vertebrate Mitochondrial", code: *b"KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"--------*-*-MMMM------------------------------M-*-*-------------" },
    Gencode { id: 3, name: "Yeast Mitochondrial", code: *b"KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"------------M-M-------------------------------M-*-*-------------" },
    Gencode { id: 4, name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"------------MMMM--------------M---------------M-*-*---------M-M-" },
    Gencode { id: 5, name: "Invertebrate Mitochondrial", code: *b"KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"------------MMMM------------------------------M-*-*-----------M-" },
    Gencode { id: 6, name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF", stop: *b"--------------M-----------------------------------------*-------" },
    Gencode { id: 9, name: "Echinoderm Mitochondrial; Flatworm Mitochondrial", code: *b"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"--------------M-------------------------------M-*-*-------------" },
    Gencode { id: 10, name: "Euplotid Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF", stop: *b"--------------M---------------------------------*-*-------------" },
    Gencode { id: 11, name: "Bacterial, Archaeal and Plant Plastid", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", stop: *b"------------MMMM--------------M---------------M-*-*-----*-----M-" },
    Gencode { id: 12, name: "Alternative Yeast Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", stop: *b"--------------M---------------M-----------------*-*-----*-------" },
    Gencode { id: 13, name: "Ascidian Mitochondrial", code: *b"KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"------------M-M-------------------------------M-*-*-----------M-" },
    Gencode { id: 14, name: "Alternative Flatworm Mitochondrial", code: *b"NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF", stop: *b"--------------M-----------------------------------*-------------" },
    Gencode { id: 15, name: "Blepharisma Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF", stop: *b"--------------M---------------------------------*-------*-------" },
    Gencode { id: 16, name: "Chlorophycean Mitochondrial", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF", stop: *b"--------------M---------------------------------*-------*-------" },
    Gencode { id: 21, name: "Trematode Mitochondrial", code: *b"NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"--------------M-------------------------------M-*-*-------------" },
    Gencode { id: 22, name: "Scenedesmus obliquus Mitochondrial Code", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF", stop: *b"--------------M---------------------------------*---*---*-------" },
    Gencode { id: 23, name: "Thraustochytrium mitochondrial code", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF", stop: *b"--------------MM------------------------------M-*-*-----*---*---" },
    Gencode { id: 24, name: "Pterobranchia Mitochondrial", code: *b"KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", stop: *b"--------------M---------------M---------------M-*-*-----------M-" },
    Gencode { id: 25, name: "Candidate Division SR1 and Gracilibacteria", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF", stop: *b"--------------M-------------------------------M-*-*-----------M-" },
    Gencode { id: 26, name: "Pachysolen tannophilus Nuclear Code", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLALEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", stop: *b"--------------M---------------M-----------------*-*-----*-------" },
    Gencode { id: 27, name: "Karyorelict Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLF", stop: *b"--------------M-----------------------------------------*-------" },
    Gencode { id: 28, name: "Condylostoma Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLF", stop: *b"--------------M---------------------------------*-*-----*-------" },
    Gencode { id: 29, name: "Mesodinium Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYYYSSSS*CWCLFLF", stop: *b"--------------M-----------------------------------------*-------" },
    Gencode { id: 30, name: "Peritrich Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSS*CWCLFLF", stop: *b"--------------M-----------------------------------------*-------" },
    Gencode { id: 31, name: "Blastocrithidia Nuclear", code: *b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSSWCWCLFLF", stop: *b"--------------M---------------------------------*-*-------------" },
    Gencode { id: 33, name: "Cephalodiscidae Mitochondrial UAA-Tyr", code: *b"KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF", stop: *b"--------------M---------------M---------------M---*-----------M-" },
];

/// Base codes on the forward strand (bcftools' `nt4`): A, C, G, T (either case) are 0..=3,
/// everything else 4. bcftools reads past its 117-entry table for bytes above `t`; those
/// count as 4 here.
const fn nt4_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b'a' as usize] = 0;
    t[b'c' as usize] = 1;
    t[b'g' as usize] = 2;
    t[b't' as usize] = 3;
    t
}

/// Complemented base codes (bcftools' `cnt4`).
const fn cnt4_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 3;
    t[b'C' as usize] = 2;
    t[b'G' as usize] = 1;
    t[b'T' as usize] = 0;
    t[b'a' as usize] = 3;
    t[b'c' as usize] = 2;
    t[b'g' as usize] = 1;
    t[b't' as usize] = 0;
    t
}

static NT4: [u8; 256] = nt4_table();
static CNT4: [u8; 256] = cnt4_table();

impl Gencode {
    /// The codon index of three forward-strand bases; above 63 with a non-ACGT first base.
    /// As in bcftools, a non-ACGT second or third base can still give an index below 64.
    #[inline]
    fn fwd_index(codon: [u8; 3]) -> usize {
        (usize::from(NT4[usize::from(codon[0])]) << 4)
            | (usize::from(NT4[usize::from(codon[1])]) << 2)
            | usize::from(NT4[usize::from(codon[2])])
    }

    /// The codon index of three reference-strand bases read on the reverse strand.
    #[inline]
    fn rev_index(codon: [u8; 3]) -> usize {
        (usize::from(CNT4[usize::from(codon[2])]) << 4)
            | (usize::from(CNT4[usize::from(codon[1])]) << 2)
            | usize::from(CNT4[usize::from(codon[0])])
    }

    /// The amino acid of a forward-strand codon (bcftools' `dna2aa`).
    #[inline]
    pub(super) fn dna2aa(&self, codon: [u8; 3]) -> u8 {
        self.code
            .get(Self::fwd_index(codon))
            .copied()
            .unwrap_or(b'X')
    }

    /// The amino acid of a reverse-strand codon given on the reference strand (`cdna2aa`).
    #[inline]
    pub(super) fn cdna2aa(&self, codon: [u8; 3]) -> u8 {
        self.code
            .get(Self::rev_index(codon))
            .copied()
            .unwrap_or(b'X')
    }

    /// The start/stop mark of a forward-strand codon: `M`, `*`, or `-` (`dna2stop`).
    #[inline]
    pub(super) fn dna2stop(&self, codon: [u8; 3]) -> u8 {
        self.stop.get(Self::fwd_index(codon)).copied().unwrap_or(0)
    }

    /// The start/stop mark of a reverse-strand codon (`cdna2stop`).
    #[inline]
    pub(super) fn cdna2stop(&self, codon: [u8; 3]) -> u8 {
        self.stop.get(Self::rev_index(codon)).copied().unwrap_or(0)
    }
}

/// What `-C/--genetic-code` asks for.
pub(super) enum GencodeChoice {
    Table(&'static Gencode),
    /// `l`: print the tables and exit.
    List,
}

/// Picks the table of `-C` (bcftools' `init_gencode`); `None` is the default table 0.
pub(super) fn choose(arg: Option<&str>) -> anyhow::Result<GencodeChoice> {
    let s = arg.unwrap_or("0");
    if s.eq_ignore_ascii_case("l") {
        return Ok(GencodeChoice::List);
    }
    let b = s.as_bytes();
    let (id, end) = pasta::variant::number::strtol(b, 0);
    if end != b.len() {
        bail!("Could not parse argument: --genetic-code {s}");
    }
    match GENCODE_TABLES.iter().find(|g| i64::from(g.id) == id) {
        Some(g) => Ok(GencodeChoice::Table(g)),
        None => bail!("Could not parse argument: --genetic-code {s}, no such table"),
    }
}

/// Prints the tables (`-C l`).
pub(super) fn print_tables(out: &mut impl Write) -> std::io::Result<()> {
    const ACGT: &[u8; 4] = b"ACGT";
    writeln!(
        out,
        "# The tables are ordered by codon as AAA,AAC,AAG,...,ACA,ACC,...,TTT:"
    )?;
    write!(out, "#\n#\t")?;
    for &b in ACGT {
        out.write_all(&[b; 16])?;
    }
    write!(out, "\n#\t")?;
    for _ in 0..4 {
        for &b in ACGT {
            out.write_all(&[b; 4])?;
        }
    }
    write!(out, "\n#\t")?;
    for _ in 0..16 {
        out.write_all(ACGT)?;
    }
    write!(out, "\n#\n\n")?;
    for g in &GENCODE_TABLES {
        writeln!(
            out,
            "{}\t{}\n\t{}\n\t{}\n",
            g.id,
            g.name,
            String::from_utf8_lossy(&g.code),
            String::from_utf8_lossy(&g.stop)
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_code_translates_both_strands() {
        let GencodeChoice::Table(g) = choose(None).expect("table 0 exists") else {
            panic!("table 0 is not a listing");
        };
        assert_eq!(g.dna2aa(*b"ATG"), b'M');
        assert_eq!(g.dna2stop(*b"ATG"), b'M');
        assert_eq!(g.dna2stop(*b"TAA"), b'*');
        // CAT on the reference strand is ATG on the reverse strand.
        assert_eq!(g.cdna2aa(*b"CAT"), b'M');
        assert_eq!(g.dna2aa(*b"NAA"), b'X');
        // bcftools' quirk: a non-ACGT middle base still indexes the table.
        assert_eq!(g.dna2aa(*b"ANA"), g.code[16]);
    }
}
