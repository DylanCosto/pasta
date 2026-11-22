mod flags;
pub mod parameter;

use std::{io, num::NonZero};

pub use self::{flags::Flags, parameter::Parameter};
use self::{flags::read_flags, parameter::fqz_decode_single_param};
use crate::cram::io::reader::num::read_u8;

// § 6.2 "FQZComp Data Stream" (2023-03-15).
const SELECTOR_TABLE_SIZE: usize = 256;

pub struct Parameters {
    pub gflags: Flags,
    selector_table: Option<Vec<u8>>,
    pub params: Vec<Parameter>,
    pub max_symbol_count: NonZero<usize>,
    pub selector_count: Option<NonZero<usize>>,
}

impl Parameters {
    pub(super) fn selector_table(&self) -> Option<&[u8]> {
        self.selector_table.as_deref()
    }
}

pub fn fqz_decode_params(src: &mut &[u8]) -> io::Result<Parameters> {
    read_version(src)?;

    let flags = read_flags(src)?;

    let mut selector_count = None;

    let parameter_count = if flags.has_parameter_count() {
        let n = read_parameter_count(src)?;
        selector_count = Some(n.checked_add(1).expect("attempt to add with overflow"));
        n
    } else {
        NonZero::<usize>::MIN
    };

    let selector_table = if flags.has_selector_table() {
        selector_count = read_selector_count(src).map(Some)?;
        read_selector_table(src).map(Some)?
    } else {
        None
    };

    let params: Vec<_> = (0..parameter_count.get())
        .map(|_| fqz_decode_single_param(src))
        .collect::<io::Result<_>>()?;

    // SAFETY: `params` is nonempty.
    let max_symbol_count = params
        .iter()
        .map(|param| param.symbol_count())
        .max()
        .unwrap();

    Ok(Parameters {
        gflags: flags,
        selector_table,
        params,
        max_symbol_count,
        selector_count,
    })
}

fn read_version(src: &mut &[u8]) -> io::Result<()> {
    const VERSION: u8 = 5;

    match read_u8(src)? {
        VERSION => Ok(()),
        n => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid version: expected {VERSION}, got {n}"),
        )),
    }
}

fn read_parameter_count(src: &mut &[u8]) -> io::Result<NonZero<usize>> {
    read_u8(src).and_then(|n| {
        NonZero::try_from(usize::from(n)).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_max_selector(src: &mut &[u8]) -> io::Result<NonZero<u8>> {
    read_u8(src).and_then(|n| {
        NonZero::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_selector_count(src: &mut &[u8]) -> io::Result<NonZero<usize>> {
    read_max_selector(src).map(|n| {
        NonZero::<usize>::from(n)
            .checked_add(1)
            .expect("attempt to add with overflow")
    })
}

fn read_selector_table(src: &mut &[u8]) -> io::Result<Vec<u8>> {
    read_array(src, SELECTOR_TABLE_SIZE)
}

/// The most run lengths an array is read from, as in htscodecs.
const MAX_RUN_COUNT: usize = 1024;

/// Reads an array of `n` values that do not fall, as the number of times each value from 0 is
/// repeated.
///
/// A run of 255 or more continues in the next length, and a length repeated is followed by how
/// many more times it repeats. A run past the end of the array is cut short.
pub fn read_array(src: &mut &[u8], n: usize) -> io::Result<Vec<u8>> {
    fn invalid_runs() -> io::Error {
        io::Error::new(io::ErrorKind::InvalidData, "invalid fqzcomp array")
    }

    let mut runs = Vec::new();
    let mut len = 0;
    let mut last = None;

    while len < n {
        let run = read_u8(src)?;

        runs.push(run);
        len += usize::from(run);

        if last == Some(run) {
            let copies = read_u8(src).map(usize::from)?;

            runs.extend(std::iter::repeat_n(run, copies));
            len += usize::from(run) * copies;
        }

        if runs.len() > MAX_RUN_COUNT {
            return Err(invalid_runs());
        }

        last = Some(run);
    }

    let mut array = Vec::with_capacity(n);
    let mut runs = runs.into_iter();
    let mut value = 0;

    while array.len() < n {
        let mut run_len = 0;

        loop {
            let part = runs.next().ok_or_else(invalid_runs)?;
            run_len += usize::from(part);

            if part != u8::MAX {
                break;
            }
        }

        if run_len > 0 {
            let value = u8::try_from(value).map_err(|_| invalid_runs())?;
            array.extend(std::iter::repeat_n(value, run_len.min(n - array.len())));
        }

        value += 1;
    }

    Ok(array)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_array() -> io::Result<()> {
        // Values 0 to 9, in runs of 1, 2, 3, 3, 3, 3, 3, 4, 4 and 5: a run length repeated is
        // followed by how many more times it repeats.
        let src = [1, 2, 3, 3, 3, 4, 4, 0, 5];
        let expected: Vec<u8> = [1, 2, 3, 3, 3, 3, 3, 4, 4, 5]
            .into_iter()
            .zip(0..)
            .flat_map(|(len, value)| std::iter::repeat_n(value, len))
            .collect();
        assert_eq!(read_array(&mut &src[..], 31)?, expected);

        // A run of 255 or more continues in the next length.
        assert_eq!(read_array(&mut &[255, 45][..], 300)?, [0; 300]);

        let mut expected = vec![0; 255];
        expected.push(1);
        assert_eq!(read_array(&mut &[255, 0, 1][..], 256)?, expected);

        // What follows the array is left.
        let mut src = &[2, 2, 0, 9][..];
        assert_eq!(read_array(&mut src, 4)?, [0, 0, 1, 1]);
        assert_eq!(src, [9]);

        Ok(())
    }

    #[test]
    fn test_read_array_without_the_first_value() -> io::Result<()> {
        // No 0, then two 1s and two 2s: the first length repeats nothing before it.
        assert_eq!(read_array(&mut &[0, 2, 2, 0][..], 4)?, [1, 1, 2, 2]);
        assert_eq!(read_array(&mut &[0, 0, 0, 3][..], 3)?, [2, 2, 2]);

        Ok(())
    }

    #[test]
    fn test_read_array_with_a_run_past_the_end() -> io::Result<()> {
        // As in htscodecs, the run is cut short.
        assert_eq!(read_array(&mut &[5][..], 3)?, [0, 0, 0]);
        assert_eq!(read_array(&mut &[1, 200][..], 4)?, [0, 1, 1, 1]);

        Ok(())
    }

    #[test]
    fn test_read_array_with_invalid_runs() {
        fn kind(src: &[u8], n: usize) -> io::ErrorKind {
            read_array(&mut &src[..], n).unwrap_err().kind()
        }

        // The lengths end early.
        assert_eq!(kind(&[], 4), io::ErrorKind::UnexpectedEof);
        assert_eq!(kind(&[2, 1], 4), io::ErrorKind::UnexpectedEof);
        assert_eq!(kind(&[1, 1], 4), io::ErrorKind::UnexpectedEof);

        // A run of 255 or more without its last length.
        assert_eq!(kind(&[1, 255], 4), io::ErrorKind::InvalidData);

        // 257 runs of 1, and so a value past 255.
        assert_eq!(
            kind(&[1, 1, 255, 1, 1, 41], 300),
            io::ErrorKind::InvalidData
        );

        // More runs of nothing than any array has.
        let mut src = vec![0];
        src.extend([0, 255].repeat(8));
        src.push(4);
        assert_eq!(kind(&src, 4), io::ErrorKind::InvalidData);
    }
}
