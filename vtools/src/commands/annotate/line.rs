//! A line of a tab-delimited annotation file split into columns (bcftools' `annot_line_t`).

use anyhow::bail;

use super::{columns::TabColumns, vcmp::AlleleList};

/// A tab-delimited annotation line, its columns, and (with `-c ...,REF,ALT`) its alleles.
#[derive(Debug, Default)]
pub(super) struct AnnotLine {
    line: Vec<u8>,
    /// Column ranges in `line`.
    cols: Vec<(usize, usize)>,
    /// Allele ranges in `line`: REF, then the comma-separated ALTs.
    als: Vec<(usize, usize)>,
    pub(super) rid: i32,
    /// The 0-based inclusive start and end.
    pub(super) start: i64,
    pub(super) end: i64,
}

impl AnnotLine {
    /// A line of `ncols` missing values (`.`), used to append missing values for ALTs without
    /// an annotation (bcftools' `aline_missing`).
    pub(super) fn missing(ncols: usize) -> Self {
        let mut l = Self::default();
        for i in 0..ncols {
            if i > 0 {
                l.line.push(b'\t');
            }
            l.cols.push((l.line.len(), l.line.len() + 1));
            l.line.push(b'.');
        }
        l
    }

    /// Splits a line into columns and alleles (bcftools' `parse_annot_line`).
    pub(super) fn parse(&mut self, text: &[u8], tab: &TabColumns) -> anyhow::Result<()> {
        self.line.clear();
        self.line.extend_from_slice(text);
        self.cols.clear();
        let mut beg = 0;
        for (i, &c) in self.line.iter().enumerate() {
            if c == b'\t' {
                self.cols.push((beg, i));
                beg = i + 1;
            }
        }
        self.cols.push((beg, self.line.len()));
        self.als.clear();
        if tab.ref_ != -1 {
            let ncols = self.cols.len();
            for idx in [tab.ref_, tab.alt] {
                if idx >= ncols as i32 {
                    bail!(
                        "Could not parse the line, expected {}+ columns, found {ncols}:\n\t{}",
                        idx + 1,
                        String::from_utf8_lossy(text)
                    );
                }
            }
            self.als.push(self.cols[tab.ref_ as usize]);
            if tab.alt >= 0 {
                let (b, e) = self.cols[tab.alt as usize];
                let mut s = b;
                for i in b..e {
                    if self.line[i] == b',' {
                        self.als.push((s, i));
                        s = i + 1;
                    }
                }
                self.als.push((s, e));
            }
        }
        Ok(())
    }

    /// Column `i`, empty when the line has fewer columns.
    #[inline]
    pub(super) fn col(&self, i: i32) -> &[u8] {
        usize::try_from(i)
            .ok()
            .and_then(|i| self.cols.get(i))
            .map_or(&[], |&(b, e)| &self.line[b..e])
    }

    /// The number of columns.
    #[inline]
    pub(super) fn ncols(&self) -> usize {
        self.cols.len()
    }

    /// The number of alleles (0 without REF and ALT columns).
    #[inline]
    pub(super) fn nals(&self) -> usize {
        self.als.len()
    }
}

impl AlleleList for AnnotLine {
    fn count(&self) -> usize {
        self.als.len()
    }

    fn allele_at(&self, i: usize) -> &[u8] {
        let (b, e) = self.als[i];
        &self.line[b..e]
    }
}
