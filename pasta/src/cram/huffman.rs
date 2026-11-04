use std::{collections::HashMap, io};

use crate::cram::io::BitReader;

type CodeBook = HashMap<i32, (i32, u32)>;

pub struct CanonicalHuffmanDecoder {
    code_book_by_len: HashMap<u32, Vec<(i32, i32)>>,
    sorted_lens: Vec<u32>,
}

impl CanonicalHuffmanDecoder {
    /// Builds a decoder for the canonical code given by `alphabet` and `bit_lens`.
    ///
    /// # Errors
    ///
    /// Returns an error if the code lengths do not describe a canonical code; see [`validate`].
    pub fn new(alphabet: &[i32], bit_lens: &[u32]) -> io::Result<Self> {
        let code_book = build_canonical_code_book(alphabet, bit_lens)?;

        let mut code_book_by_len = HashMap::new();

        for (symbol, (code, len)) in code_book {
            let entries = code_book_by_len.entry(len).or_insert_with(Vec::new);
            entries.push((symbol, code));
        }

        let sorted_lens = {
            let mut lens: Vec<_> = code_book_by_len.keys().copied().collect();
            lens.sort_unstable();
            lens
        };

        Ok(Self {
            code_book_by_len,
            sorted_lens,
        })
    }

    pub fn decode(&self, reader: &mut BitReader<'_>) -> io::Result<i32> {
        if self.sorted_lens.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "cannot decode with an empty Huffman alphabet",
            ));
        }

        let mut prev_len = 0;
        let mut input_code = 0;

        for &len in &self.sorted_lens {
            input_code <<= len - prev_len;

            let b = reader.read_i32(len - prev_len)?;
            input_code |= b;

            let entry = self.code_book_by_len[&len]
                .iter()
                .find(|(_, code)| input_code == *code);

            if let Some((symbol, _)) = entry {
                return Ok(*symbol);
            }

            prev_len = len;
        }

        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "could not find symbol",
        ))
    }
}

/// The longest code HTSlib accepts (`max_code_bits`), so that every code fits in an `i32`.
const MAX_CODE_LEN: u32 = 31;

/// Checks that `alphabet` and `bit_lens` describe a canonical Huffman code, as HTSlib's
/// `cram_huffman_decode_init` does.
///
/// An empty alphabet is valid here, but decoding with it fails.
///
/// # Errors
///
/// Returns an error if the two lists differ in length, a code is longer than 31 bits, or there
/// are more codes of some length than the code space allows.
pub fn validate(alphabet: &[i32], bit_lens: &[u32]) -> io::Result<()> {
    build_canonical_code_book(alphabet, bit_lens).map(|_| ())
}

fn build_canonical_code_book(alphabet: &[i32], bit_lens: &[u32]) -> io::Result<CodeBook> {
    if alphabet.len() != bit_lens.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Huffman alphabet has {} symbols but {} bit lengths",
                alphabet.len(),
                bit_lens.len()
            ),
        ));
    }

    if let Some(&len) = bit_lens.iter().find(|&&len| len > MAX_CODE_LEN) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Huffman code length ({len}) is greater than maximum supported ({MAX_CODE_LEN})"
            ),
        ));
    }

    let sorted_alphabet = {
        let mut pairs: Vec<_> = alphabet.iter().zip(bit_lens.iter()).collect();
        pairs.sort_by_key(|&(symbol, bit_len)| (bit_len, symbol));
        pairs
    };

    let mut code_book = CodeBook::with_capacity(sorted_alphabet.len());

    // As in HTSlib, each code is one more than the previous one, shifted left to the new length.
    // Lengths are at most 31 bits, so the codes and their bounds fit in an `i64`.
    let mut code: i64 = -1;
    let mut max_code: i64 = 0;
    let mut prev_bit_len = 0;

    for (&symbol, &bit_len) in sorted_alphabet {
        code += 1;

        if code > max_code {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Huffman code lengths are oversubscribed",
            ));
        }

        if bit_len > prev_bit_len {
            code <<= bit_len - prev_bit_len;
            max_code = (1 << bit_len) - 1;
            prev_bit_len = bit_len;
        }

        let code =
            i32::try_from(code).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        code_book.insert(symbol, (code, bit_len));
    }

    Ok(code_book)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() -> io::Result<()> {
        let symbols = [0x4e, 0x44, 0x4c];
        let bit_lens = [1, 2, 2];
        let decoder = CanonicalHuffmanDecoder::new(&symbols, &bit_lens)?;

        let data = [0b01011000];
        let mut reader = BitReader::new(&data[..]);

        assert_eq!(decoder.decode(&mut reader)?, 0x4e);
        assert_eq!(decoder.decode(&mut reader)?, 0x44);
        assert_eq!(decoder.decode(&mut reader)?, 0x4c);
        assert_eq!(decoder.decode(&mut reader)?, 0x4e);

        Ok(())
    }

    #[test]
    fn test_build_canonical_code_book() -> io::Result<()> {
        let symbols = [65, 66, 67, 68, 69, 70];
        let bit_lens = [1, 3, 3, 3, 4, 4];

        let code_book = build_canonical_code_book(&symbols, &bit_lens)?;

        assert_eq!(code_book.len(), 6);

        assert_eq!(code_book[&65], (0b0, 1));
        assert_eq!(code_book[&66], (0b100, 3));
        assert_eq!(code_book[&67], (0b101, 3));
        assert_eq!(code_book[&68], (0b110, 3));
        assert_eq!(code_book[&69], (0b1110, 4));
        assert_eq!(code_book[&70], (0b1111, 4));

        Ok(())
    }

    #[test]
    fn test_build_canonical_code_book_with_invalid_lengths() {
        fn t(alphabet: &[i32], bit_lens: &[u32]) {
            assert!(
                matches!(
                    build_canonical_code_book(alphabet, bit_lens),
                    Err(e) if e.kind() == io::ErrorKind::InvalidData
                ),
                "alphabet {alphabet:?}, bit_lens {bit_lens:?}"
            );
        }

        t(&[65, 66], &[1]);
        t(&[65], &[32]);
        t(&[65, 66], &[0, 0]);
        t(&[65, 66, 67], &[1, 1, 1]);
        t(&[65, 66, 67], &[1, 1, 31]);
    }

    #[test]
    fn test_decode_with_empty_alphabet() -> io::Result<()> {
        assert!(validate(&[], &[]).is_ok());

        let decoder = CanonicalHuffmanDecoder::new(&[], &[])?;
        let mut reader = BitReader::new(&[0xff][..]);

        assert!(matches!(
            decoder.decode(&mut reader),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
