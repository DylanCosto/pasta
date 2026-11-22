use std::{
    io::{self, Write},
    num::NonZero,
};

use super::RangeCoder;

const MAX_SYMBOL_COUNT: usize = 1 << 8;

/// What coding a symbol adds to its frequency.
const STEP: u16 = 16;

/// The total above which the frequencies are halved.
const MAX_TOTAL_FREQUENCY: u32 = (1 << 16) - 17;

/// A symbol and its frequency.
///
/// A model is a slice of these. The first holds the total of the frequencies, and no symbol; the
/// rest hold the symbols, kept roughly in order of falling frequency. A frequency is never above
/// the total, so the comparison that moves a symbol ahead of its neighbor leaves the first
/// symbol where it is without a test for it.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Cell {
    frequency: u16,
    symbol: u8,
}

/// The cells of a model in which each of `symbol_count` symbols has been seen once.
fn initial_cells(symbol_count: NonZero<usize>) -> Vec<Cell> {
    let len = symbol_count.get();

    assert!(len <= MAX_SYMBOL_COUNT);

    let total = Cell {
        // `len` <= 256.
        frequency: len as u16,
        symbol: 0,
    };

    let symbols = (0..=u8::MAX).take(len).map(|symbol| Cell {
        frequency: 1,
        symbol,
    });

    std::iter::once(total).chain(symbols).collect()
}

fn decode(cells: &mut [Cell], src: &mut &[u8], range_coder: &mut RangeCoder) -> io::Result<u8> {
    let freq = range_coder.range_get_freq(u32::from(cells[0].frequency));

    let mut acc = 0;
    let mut x = 1;

    let frequency = loop {
        // A valid stream gives a frequency below the total, which some symbol's range holds.
        let Some(cell) = cells.get(x) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid arithmetic-coded symbol",
            ));
        };

        let frequency = u32::from(cell.frequency);

        if acc + frequency > freq {
            break frequency;
        }

        acc += frequency;
        x += 1;
    };

    range_coder.range_decode(src, acc, frequency)?;

    Ok(update(cells, x))
}

fn encode<W>(
    cells: &mut [Cell],
    writer: &mut W,
    range_coder: &mut RangeCoder,
    sym: u8,
) -> io::Result<()>
where
    W: Write,
{
    let mut acc = 0;
    let mut x = 1;

    while cells[x].symbol != sym {
        acc += u32::from(cells[x].frequency);
        x += 1;
    }

    range_coder.range_encode(
        writer,
        acc,
        u32::from(cells[x].frequency),
        u32::from(cells[0].frequency),
    )?;

    update(cells, x);

    Ok(())
}

/// Counts the symbol of cell `x`, which is returned.
fn update(cells: &mut [Cell], x: usize) -> u8 {
    let symbol = cells[x].symbol;

    // A frequency is at most the total, which is at most `MAX_TOTAL_FREQUENCY` between symbols.
    cells[x].frequency += STEP;
    let total = u32::from(cells[0].frequency) + u32::from(STEP);

    if total > MAX_TOTAL_FREQUENCY {
        let mut total = 0;

        for cell in &mut cells[1..] {
            cell.frequency -= cell.frequency / 2;
            total += u32::from(cell.frequency);
        }

        cells[0].frequency = total as u16;
    } else {
        cells[0].frequency = total as u16;
    }

    // Before the first symbol is the total, which it never exceeds.
    if cells[x].frequency > cells[x - 1].frequency {
        cells.swap(x, x - 1);
    }

    symbol
}

/// An adaptive model of up to 256 symbols.
#[derive(Clone, Debug)]
pub struct Model {
    cells: Vec<Cell>,
}

impl Model {
    pub fn new(symbol_count: NonZero<usize>) -> Self {
        Self {
            cells: initial_cells(symbol_count),
        }
    }

    pub fn decode(&mut self, src: &mut &[u8], range_coder: &mut RangeCoder) -> io::Result<u8> {
        decode(&mut self.cells, src, range_coder)
    }

    pub fn encode<W>(
        &mut self,
        writer: &mut W,
        range_coder: &mut RangeCoder,
        sym: u8,
    ) -> io::Result<()>
    where
        W: Write,
    {
        encode(&mut self.cells, writer, range_coder, sym)
    }
}

/// Models of the same symbols, chosen by a context, in one allocation.
#[derive(Clone, Debug)]
pub struct ModelBank {
    cells: Vec<Cell>,
    /// The cells of a model.
    stride: usize,
}

impl ModelBank {
    pub fn new(symbol_count: NonZero<usize>, model_count: usize) -> Self {
        let model = initial_cells(symbol_count);

        Self {
            stride: model.len(),
            cells: model.repeat(model_count),
        }
    }

    /// Decodes a symbol with model `i`.
    ///
    /// # Panics
    ///
    /// Panics if there is no model `i`.
    #[inline]
    pub fn decode(
        &mut self,
        i: usize,
        src: &mut &[u8],
        range_coder: &mut RangeCoder,
    ) -> io::Result<u8> {
        decode(self.model_mut(i), src, range_coder)
    }

    /// Encodes a symbol with model `i`.
    ///
    /// # Panics
    ///
    /// Panics if there is no model `i`, or it does not have the symbol.
    #[inline]
    pub fn encode<W>(
        &mut self,
        i: usize,
        writer: &mut W,
        range_coder: &mut RangeCoder,
        sym: u8,
    ) -> io::Result<()>
    where
        W: Write,
    {
        encode(self.model_mut(i), writer, range_coder, sym)
    }

    fn model_mut(&mut self, i: usize) -> &mut [Cell] {
        &mut self.cells[i * self.stride..][..self.stride]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn symbol_count(n: usize) -> NonZero<usize> {
        NonZero::new(n).unwrap()
    }

    fn frequencies(model: &Model) -> Vec<(u8, u16)> {
        model.cells[1..]
            .iter()
            .map(|cell| (cell.symbol, cell.frequency))
            .collect()
    }

    #[test]
    fn test_encode_orders_symbols_by_frequency() -> io::Result<()> {
        let mut model = Model::new(symbol_count(3));
        let mut range_coder = RangeCoder::default();
        let mut dst = Vec::new();

        assert_eq!(model.cells[0].frequency, 3);
        assert_eq!(frequencies(&model), [(0, 1), (1, 1), (2, 1)]);

        model.encode(&mut dst, &mut range_coder, 2)?;
        assert_eq!(model.cells[0].frequency, 19);
        assert_eq!(frequencies(&model), [(0, 1), (2, 17), (1, 1)]);

        model.encode(&mut dst, &mut range_coder, 2)?;
        assert_eq!(frequencies(&model), [(2, 33), (0, 1), (1, 1)]);

        // The first symbol stays first.
        model.encode(&mut dst, &mut range_coder, 2)?;
        assert_eq!(model.cells[0].frequency, 51);
        assert_eq!(frequencies(&model), [(2, 49), (0, 1), (1, 1)]);

        Ok(())
    }

    #[test]
    fn test_encode_halves_the_frequencies_past_the_maximum_total() -> io::Result<()> {
        let mut model = Model::new(symbol_count(2));
        let mut range_coder = RangeCoder::default();
        let mut dst = Vec::new();

        // 2 + 16 * 4094 = 65506, and one more symbol passes 65519.
        for _ in 0..4094 {
            model.encode(&mut dst, &mut range_coder, 1)?;
        }

        assert_eq!(model.cells[0].frequency, 65506);
        assert_eq!(frequencies(&model), [(1, 65505), (0, 1)]);

        model.encode(&mut dst, &mut range_coder, 1)?;

        // 65521 - 32760 = 32761, and 1 - 0 = 1.
        assert_eq!(model.cells[0].frequency, 32762);
        assert_eq!(frequencies(&model), [(1, 32761), (0, 1)]);

        Ok(())
    }

    #[test]
    fn test_decode_follows_encode() -> io::Result<()> {
        let symbols: Vec<u8> = (0..20_000u32)
            .map(|i| (i.wrapping_mul(2_654_435_761) >> 24) as u8 % 7 % 5)
            .collect();

        let mut encoder = Model::new(symbol_count(5));
        let mut range_coder = RangeCoder::default();
        let mut dst = Vec::new();

        for &sym in &symbols {
            encoder.encode(&mut dst, &mut range_coder, sym)?;
        }

        range_coder.range_encode_end(&mut dst)?;

        let mut src = &dst[..];
        let mut decoder = Model::new(symbol_count(5));
        let mut range_coder = RangeCoder::new(&mut src)?;

        for &sym in &symbols {
            assert_eq!(decoder.decode(&mut src, &mut range_coder)?, sym);
        }

        assert_eq!(decoder.cells, encoder.cells);

        Ok(())
    }

    #[test]
    fn test_decode_rejects_a_frequency_past_the_total() -> io::Result<()> {
        // A code of all ones is past the range of the last of 3 symbols.
        let mut src = &[0x00, 0xff, 0xff, 0xff, 0xff][..];
        let mut range_coder = RangeCoder::new(&mut src)?;
        let mut model = Model::new(symbol_count(3));

        assert!(matches!(
            model.decode(&mut src, &mut range_coder),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_model_bank_keeps_a_model_per_context() -> io::Result<()> {
        let contexts = [0, 2, 2, 1, 0, 2, 2, 2, 1, 0, 0, 2];
        let symbols = [3, 1, 1, 0, 3, 1, 2, 1, 0, 3, 3, 1];

        let mut bank = ModelBank::new(symbol_count(4), 3);
        let mut models = vec![Model::new(symbol_count(4)); 3];

        let mut bank_dst = Vec::new();
        let mut bank_coder = RangeCoder::default();
        let mut dst = Vec::new();
        let mut range_coder = RangeCoder::default();

        for (&i, &sym) in contexts.iter().zip(&symbols) {
            bank.encode(i, &mut bank_dst, &mut bank_coder, sym)?;
            models[i].encode(&mut dst, &mut range_coder, sym)?;
        }

        bank_coder.range_encode_end(&mut bank_dst)?;
        range_coder.range_encode_end(&mut dst)?;

        assert_eq!(bank_dst, dst);

        for (i, model) in models.iter().enumerate() {
            assert_eq!(bank.model_mut(i), &model.cells[..]);
        }

        let mut src = &bank_dst[..];
        let mut bank = ModelBank::new(symbol_count(4), 3);
        let mut range_coder = RangeCoder::new(&mut src)?;

        for (&i, &sym) in contexts.iter().zip(&symbols) {
            assert_eq!(bank.decode(i, &mut src, &mut range_coder)?, sym);
        }

        Ok(())
    }
}
