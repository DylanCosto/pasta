//! Ploidy by region and sex (bcftools' `ploidy.c`) and the predefined ploidy sets of
//! `--ploidy` (from `vcfcall.c`).

use std::{collections::HashMap, io::Write};

use anyhow::bail;
use pasta::variant::index::split_idx_name;

use super::regidx::{RegIdx, parse_tab};

/// The ploidy of one sex in a region (`sex_ploidy_t`).
#[derive(Clone, Copy, Debug)]
struct SexPloidy {
    sex: usize,
    ploidy: i32,
}

/// bcftools' `ploidy_t`.
#[derive(Clone, Debug, Default)]
pub struct Ploidy {
    nsex: usize,
    dflt: i32,
    min: i32,
    max: i32,
    sex2dflt: Vec<i32>,
    idx: RegIdx<SexPloidy>,
    sex2id: HashMap<Vec<u8>, usize>,
    id2sex: Vec<Vec<u8>>,
}

/// A predefined ploidy set (`ploidy_predef_t`).
struct Predefined {
    alias: &'static str,
    about: &'static str,
    ploidy: &'static str,
}

const PREDEFINED: &[Predefined] = &[
    Predefined {
        alias: "GRCh37",
        about: "Human Genome reference assembly GRCh37 / hg19",
        ploidy: "X 1 60000 M 1\n\
                 X 2699521 154931043 M 1\n\
                 Y 1 59373566 M 1\n\
                 Y 1 59373566 F 0\n\
                 MT 1 16569 M 1\n\
                 MT 1 16569 F 1\n\
                 chrX 1 60000 M 1\n\
                 chrX 2699521 154931043 M 1\n\
                 chrY 1 59373566 M 1\n\
                 chrY 1 59373566 F 0\n\
                 chrM 1 16569 M 1\n\
                 chrM 1 16569 F 1\n\
                 *  * *     M 2\n\
                 *  * *     F 2\n",
    },
    Predefined {
        alias: "GRCh38",
        about: "Human Genome reference assembly GRCh38 / hg38",
        ploidy: "X 1 9999 M 1\n\
                 X 2781480 155701381 M 1\n\
                 Y 1 57227415 M 1\n\
                 Y 1 57227415 F 0\n\
                 MT 1 16569 M 1\n\
                 MT 1 16569 F 1\n\
                 chrX 1 9999 M 1\n\
                 chrX 2781480 155701381 M 1\n\
                 chrY 1 57227415 M 1\n\
                 chrY 1 57227415 F 0\n\
                 chrM 1 16569 M 1\n\
                 chrM 1 16569 F 1\n\
                 *  * *     M 2\n\
                 *  * *     F 2\n",
    },
    Predefined {
        alias: "X",
        about: "Treat male samples as haploid and female as diploid regardless of the chromosome name",
        ploidy: "*  * *     M 1\n\
                 *  * *     F 2\n",
    },
    Predefined {
        alias: "Y",
        about: "Treat male samples as haploid and female as no-copy, regardless of the chromosome name",
        ploidy: "*  * *     M 1\n\
                 *  * *     F 0\n",
    },
    Predefined {
        alias: "1",
        about: "Treat all samples as haploid",
        ploidy: "*  * *     * 1\n",
    },
    Predefined {
        alias: "2",
        about: "Treat all samples as diploid",
        ploidy: "*  * *     * 2\n",
    },
];

/// The ploidy definition used when none is given: sexes "0", "1", and "2" with that ploidy.
pub const DEFAULT_PLOIDY: &str = "* * * 0 0\n* * * 1 1\n* * * 2 2\n";

impl Ploidy {
    fn empty() -> Self {
        Self {
            min: -1,
            max: -1,
            ..Self::default()
        }
    }

    /// Reads a ploidy file (`ploidy_init`); `FILE##idx##INDEX` reads `FILE`, as `hts_open`
    /// does.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be read or a line cannot be parsed.
    pub fn from_file(path: &str, dflt: i32) -> anyhow::Result<Self> {
        let data = split_idx_name(path).0;
        let Ok(text) = pasta::variant::readlist::read_all_lines(std::path::Path::new(data)) else {
            eprintln!(
                "[E::hts_open_format] Failed to open file \"{data}\" : No such file or directory"
            );
            bail!("Could not initialize ploidy");
        };
        let mut p = Self::empty();
        // regidx stops at the first empty line (hts_getline returns 0 for it).
        for line in text.iter().take_while(|l| !l.is_empty()) {
            if !p.parse_line(line)? {
                bail!("Could not initialize ploidy");
            }
        }
        p.finish(dflt);
        Ok(p)
    }

    /// Parses ploidy lines from a string (`ploidy_init_string`).
    ///
    /// # Errors
    ///
    /// Returns an error when a line cannot be parsed.
    pub fn from_str(text: &str, dflt: i32) -> anyhow::Result<Self> {
        let mut p = Self::empty();
        for line in text.split(['\r', '\n']) {
            let line = line.trim_start();
            if line.is_empty() {
                continue;
            }
            // Malformed lines of a string are skipped, as with bcftools.
            p.parse_line(line.as_bytes())?;
        }
        p.finish(dflt);
        Ok(p)
    }

    /// The predefined set `--ploidy ALIAS` names (`init_ploidy`), printing the list (or the
    /// definition, after `?`) and exiting as bcftools does for unknown aliases.
    ///
    /// # Errors
    ///
    /// Returns an error when the predefined lines cannot be parsed.
    pub fn predefined(alias: &str) -> anyhow::Result<Self> {
        let (alias, detailed) = match alias.strip_suffix('?') {
            Some(a) => (a, true),
            None => (alias, false),
        };
        let found = PREDEFINED
            .iter()
            .find(|p| p.alias.eq_ignore_ascii_case(alias));
        let mut err = std::io::stderr().lock();
        match found {
            None => {
                let _ = writeln!(err, "\nPRE-DEFINED PLOIDY FILES\n");
                let _ = writeln!(err, " * Columns are: CHROM,FROM,TO,SEX,PLOIDY");
                let _ = writeln!(err, " * Coordinates are 1-based inclusive.");
                let _ = writeln!(err, " * A '*' means any value not otherwise defined.\n");
                for p in PREDEFINED {
                    let _ = writeln!(err, "{}\n   .. {}\n", p.alias, p.about);
                    if detailed {
                        let _ = writeln!(err, "{}", p.ploidy);
                    }
                }
                let _ = writeln!(err, "Run as --ploidy <alias> (e.g. --ploidy GRCh37).");
                let _ = writeln!(
                    err,
                    "To see the detailed ploidy definition, append a question mark (e.g. --ploidy GRCh37?).\n"
                );
                std::process::exit(255);
            }
            Some(p) if detailed => {
                let _ = write!(err, "{}", p.ploidy);
                std::process::exit(255);
            }
            Some(p) => Self::from_str(p.ploidy, 2),
        }
    }

    /// `ploidy_parse`: one line `CHROM FROM TO SEX PLOIDY` (or `* * * SEX PLOIDY`). Returns
    /// `false` for a line with too few fields, which fails a ploidy file.
    fn parse_line(&mut self, line: &[u8]) -> anyhow::Result<bool> {
        let ws = |c: u8| c.is_ascii_whitespace();
        let text = String::from_utf8_lossy(line).into_owned();
        let mut i = line.iter().position(|&c| !ws(c)).unwrap_or(line.len());
        let is_default = line.get(i) == Some(&b'*') && line.get(i + 1).is_none_or(|&c| ws(c));
        let region = if is_default {
            None
        } else {
            match parse_tab(line) {
                Ok(Some(r)) => Some((r.0.to_vec(), r.1, r.2)),
                Ok(None) => return Ok(true),
                Err(msg) => {
                    eprintln!("{msg}");
                    return Ok(false);
                }
            }
        };
        for _ in 0..3 {
            while i < line.len() && !ws(line[i]) {
                i += 1;
            }
            if i == line.len() {
                return Ok(false);
            }
            while i < line.len() && ws(line[i]) {
                i += 1;
            }
        }
        if i == line.len() {
            return Ok(false);
        }
        let start = i;
        while i < line.len() && !ws(line[i]) {
            i += 1;
        }
        if i == line.len() || i == start {
            bail!("Could not parse: {text}");
        }
        let sex_name = &line[start..i];
        let sex = match self.sex2id.get(sex_name) {
            Some(&s) => s,
            None => {
                self.nsex += 1;
                self.id2sex.push(sex_name.to_vec());
                self.sex2id.insert(sex_name.to_vec(), self.nsex - 1);
                self.sex2dflt.push(-1);
                self.nsex - 1
            }
        };
        let (value, end) = pasta::variant::number::strtol(line, i);
        if end == i {
            bail!("Could not parse: {text}");
        }
        let ploidy = value as i32;
        if self.min < 0 || ploidy < self.min {
            self.min = ploidy;
        }
        if self.max < 0 || ploidy > self.max {
            self.max = ploidy;
        }
        match region {
            None => self.sex2dflt[sex] = ploidy,
            Some((chr, beg, end)) => self.idx.push(&chr, beg, end, SexPloidy { sex, ploidy }),
        }
        Ok(true)
    }

    /// `_set_defaults`.
    fn finish(&mut self, dflt: i32) {
        self.idx.finish();
        let dflt = match self.sex2id.get(b"*".as_slice()) {
            Some(&i) => self.sex2dflt[i],
            None => dflt,
        };
        for d in &mut self.sex2dflt {
            if *d == -1 {
                *d = dflt;
            }
        }
        self.dflt = dflt;
        if self.min < 0 || dflt < self.min {
            self.min = dflt;
        }
        if self.max < 0 || dflt > self.max {
            self.max = dflt;
        }
    }

    /// The number of sexes (`ploidy_nsex`).
    pub fn nsex(&self) -> usize {
        self.nsex
    }

    /// The ID of a sex (`ploidy_sex2id`).
    pub fn sex_id(&self, sex: &[u8]) -> Option<usize> {
        self.sex2id.get(sex).copied()
    }

    /// Adds a sex with the default ploidy (`ploidy_add_sex`).
    pub fn add_sex(&mut self, sex: &[u8]) -> usize {
        if let Some(&i) = self.sex2id.get(sex) {
            return i;
        }
        self.nsex += 1;
        self.id2sex.push(sex.to_vec());
        self.sex2dflt.push(self.dflt);
        self.sex2id.insert(sex.to_vec(), self.nsex - 1);
        self.nsex - 1
    }

    /// The largest ploidy (`ploidy_max`).
    pub fn max(&self) -> i32 {
        self.dflt.max(self.max)
    }

    /// The ploidy of each sex at a position (`ploidy_query`).
    pub fn query(&self, seq: &[u8], pos: u32, sex2ploidy: &mut [i32]) {
        let mut hits = self.idx.overlaps(seq, pos, pos).peekable();
        if hits.peek().is_none() {
            sex2ploidy.copy_from_slice(&self.sex2dflt);
            return;
        }
        sex2ploidy.fill(self.dflt);
        for h in hits {
            let sp = self.idx.payload(&h);
            if sp.ploidy != self.dflt {
                sex2ploidy[sp.sex] = sp.ploidy;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn idx_names_read_the_data_file() {
        let (dir, name) = crate::raw_input::write_idx_named("ploidy-idx", "p.txt", "* * * M 3\n");
        let p = Ploidy::from_file(&name, 2);
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(p.expect("the data file is read").max(), 3);
    }
}
