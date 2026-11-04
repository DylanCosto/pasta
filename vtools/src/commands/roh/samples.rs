//! Sample lists in header order (bcftools' `smpl_ilist_init` with `SMPL_NONE`).

use anyhow::bail;
use pasta::variant::{
    Header,
    readlist::{self, ListSource},
};

/// Whether names missing from the header are reported (`SMPL_VERBOSE`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Verbosity {
    Verbose,
    Quiet,
}

/// Sample indices in header order (`smpl_ilist_t`).
#[derive(Clone, Debug)]
pub struct SampleList(Vec<usize>);

/// The sample name of a list entry: up to the first unescaped whitespace.
fn entry_name(entry: &[u8]) -> &[u8] {
    let mut i = 0;
    while i < entry.len() {
        if entry[i].is_ascii_whitespace() {
            let backslashes = entry[..i].iter().rev().take_while(|&&c| c == b'\\').count();
            if backslashes % 2 == 0 {
                break;
            }
        }
        i += 1;
    }
    &entry[..i]
}

impl SampleList {
    /// Reads a comma-separated list or a file of sample names, `^` to exclude them, or all
    /// samples for `None`. Names not in the header are skipped (and reported with
    /// [`Verbosity::Verbose`]).
    ///
    /// # Errors
    ///
    /// Returns bcftools' error when the list cannot be read.
    pub fn new(
        hdr: &Header,
        list: Option<(&str, ListSource)>,
        verbosity: Verbosity,
    ) -> anyhow::Result<Self> {
        let nsamples = hdr.sample_count();
        let Some((list, source)) = list else {
            return Ok(Self((0..nsamples).collect()));
        };
        let (negate, body) = match list.strip_prefix('^') {
            Some(rest) => (true, rest),
            None => (false, list),
        };
        let Ok(entries) = readlist::read_list(body.as_bytes(), source) else {
            bail!("Could not parse {list}");
        };
        let mut listed = vec![false; nsamples];
        for entry in &entries {
            let name = entry_name(entry);
            match hdr.sample_id(name) {
                Some(i) => listed[i as usize] = true,
                None if verbosity == Verbosity::Verbose => {
                    eprintln!("No such sample: \"{}\"", String::from_utf8_lossy(name));
                }
                None => {}
            }
        }
        Ok(Self(
            (0..nsamples).filter(|&i| listed[i] != negate).collect(),
        ))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &usize> {
        self.0.iter()
    }

    pub fn into_vec(self) -> Vec<usize> {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_end_at_unescaped_whitespace() {
        assert_eq!(entry_name(b"NA1 extra"), b"NA1");
        assert_eq!(entry_name(b"N\\ A1\tx"), b"N\\ A1");
        assert_eq!(entry_name(b"NA1"), b"NA1");
    }
}
