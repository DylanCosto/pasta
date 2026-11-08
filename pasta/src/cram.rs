//! Handles the reading and writing of the CRAM format.

pub mod codecs;
pub mod container;
pub mod crai;
pub mod file_definition;
pub mod fs;
mod huffman;
pub mod io;
mod md5;
pub mod record;

use self::md5::Md5;

pub use self::{file_definition::FileDefinition, record::Record};

const MAGIC_NUMBER: [u8; 4] = *b"CRAM";
const MD5_OUTPUT_SIZE: usize = md5::OUTPUT_SIZE;

/// Calculates the MD5 checksum of a reference sequence as `@SQ` `M5` fields hold it.
///
/// _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.3.2 "Reference MD5
/// calculation": lowercase bases are uppercased, and characters outside `!`..`~` are skipped.
pub fn calculate_normalized_sequence_digest(sequence: &[u8]) -> [u8; MD5_OUTPUT_SIZE] {
    let mut digest = NormalizedSequenceDigest::default();
    digest.update(sequence);
    digest.finalize()
}

/// The checksum of [`calculate_normalized_sequence_digest`], calculated over the parts of a
/// sequence given in turn, so that the whole sequence need not be in memory at once.
///
/// # Examples
///
/// ```
/// use pasta::cram::{NormalizedSequenceDigest, calculate_normalized_sequence_digest};
///
/// let mut digest = NormalizedSequenceDigest::new();
/// digest.update(b"ACG");
/// digest.update(b"t\n");
/// assert_eq!(digest.finalize(), calculate_normalized_sequence_digest(b"ACGT"));
/// ```
pub struct NormalizedSequenceDigest {
    hasher: Md5,
    buf: Box<[u8; Self::CHUNK_SIZE]>,
}

impl std::fmt::Debug for NormalizedSequenceDigest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NormalizedSequenceDigest")
            .finish_non_exhaustive()
    }
}

impl Default for NormalizedSequenceDigest {
    fn default() -> Self {
        Self::new()
    }
}

impl NormalizedSequenceDigest {
    const CHUNK_SIZE: usize = 1 << 16;

    /// Starts a checksum.
    pub fn new() -> Self {
        Self {
            hasher: Md5::new(),
            buf: Box::new([0; Self::CHUNK_SIZE]),
        }
    }

    /// Adds the next part of the sequence.
    pub fn update(&mut self, sequence: &[u8]) {
        for chunk in sequence.chunks(Self::CHUNK_SIZE) {
            // Whether the chunk is all graphic characters, and whether any is lowercase, in one
            // pass that does not stop early (so that it is vectorized).
            let (is_graphic, has_lowercase) =
                chunk
                    .iter()
                    .fold((true, false), |(is_graphic, has_lowercase), b| {
                        (
                            is_graphic & b.is_ascii_graphic(),
                            has_lowercase | b.is_ascii_lowercase(),
                        )
                    });

            if is_graphic && !has_lowercase {
                self.hasher.update(chunk);
            } else if is_graphic {
                // "All lowercase characters are converted to uppercase."
                for (dst, src) in self.buf.iter_mut().zip(chunk) {
                    *dst = src.to_ascii_uppercase();
                }

                self.hasher.update(&self.buf[..chunk.len()]);
            } else {
                // "All characters outside of the inclusive range 33 ('!') to 126 ('~') are
                // stripped out."
                let mut n = 0;

                for &b in chunk.iter().filter(|b| b.is_ascii_graphic()) {
                    self.buf[n] = b.to_ascii_uppercase();
                    n += 1;
                }

                self.hasher.update(&self.buf[..n]);
            }
        }
    }

    /// Returns the checksum of the parts given.
    pub fn finalize(self) -> [u8; MD5_OUTPUT_SIZE] {
        self.hasher.finalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalized_sequence_digest_in_parts() {
        // Over many chunks, with lowercase and stripped characters in some, given in uneven parts.
        let sequence: Vec<u8> = (0..300_000u32)
            .map(|i| match i % 7919 {
                0 => b'\n',
                1..=1000 => b"acgt"[(i % 4) as usize],
                _ => b"ACGTN"[(i % 5) as usize],
            })
            .collect();

        let expected = {
            let normalized: Vec<u8> = sequence
                .iter()
                .filter(|b| b.is_ascii_graphic())
                .map(u8::to_ascii_uppercase)
                .collect();
            Md5::digest(&normalized)
        };

        assert_eq!(calculate_normalized_sequence_digest(&sequence), expected);

        let mut digest = NormalizedSequenceDigest::default();

        for part in sequence.chunks(12_345) {
            digest.update(part);
        }

        assert_eq!(digest.finalize(), expected);
    }

    #[test]
    fn test_calculate_normalized_sequence_digest() {
        assert_eq!(
            calculate_normalized_sequence_digest(b"ACGT"),
            [
                0xf1, 0xf8, 0xf4, 0xbf, 0x41, 0x3b, 0x16, 0xad, 0x13, 0x57, 0x22, 0xaa, 0x45, 0x91,
                0x04, 0x3e
            ]
        );

        assert_eq!(
            calculate_normalized_sequence_digest(b" AC\tgt\n"),
            [
                0xf1, 0xf8, 0xf4, 0xbf, 0x41, 0x3b, 0x16, 0xad, 0x13, 0x57, 0x22, 0xaa, 0x45, 0x91,
                0x04, 0x3e
            ]
        );

        // _Sequence Alignment/Map Format Specification_ (2024-11-06) § 1.3.2 "Reference MD5
        // calculation"
        assert_eq!(
            calculate_normalized_sequence_digest(b"ACGT ACGT ACGT\nacgt acgt acgt\n... 12345 !!!"),
            [
                0xdf, 0xab, 0xdb, 0xb3, 0x6e, 0x23, 0x9a, 0x6d, 0xa8, 0x89, 0x57, 0x84, 0x1f, 0x32,
                0xb8, 0xe4
            ]
        );
    }
}
