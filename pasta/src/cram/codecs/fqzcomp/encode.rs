use std::{
    io::{self, Write},
    num::NonZero,
};

use super::{
    Models,
    parameters::{self, parameter},
};
use crate::cram::{
    codecs::aac::RangeCoder,
    io::writer::num::{write_u8, write_u16_le, write_uint7},
};

/// Compresses quality scores, the records of `lens` scores each, choosing the parameters as
/// htscodecs' `fqz_compress` does with its basic strategy.
pub fn encode(lens: &[usize], src: &[u8]) -> io::Result<Vec<u8>> {
    encode_with(Strategy::Basic, lens, src)
}

/// The preset a block's parameters start from, before the scores themselves adjust them
/// (htscodecs' `strat`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "the writer does not yet try the other strategies, as HTSlib does"
    )
)]
pub(super) enum Strategy {
    /// What HTSlib uses below level 7.
    Basic,
    /// For reads like the HiSeq 2000's.
    HiSeq,
    /// For reads like the MiSeq's.
    MiSeq,
    /// For reads like Ion Torrent's: the context is two scores and no more.
    IonTorrent,
}

/// A strategy's starting bit widths and where each sits in the 16-bit context.
struct Preset {
    q_bits: u32,
    q_shift: u32,
    p_bits: u32,
    /// `None` chooses the shift from the first record's length.
    p_shift: Option<u32>,
    d_bits: u32,
    d_shift: u32,
    q_loc: u32,
    s_loc: u32,
    p_loc: u32,
    d_loc: u32,
    /// Whether a selector by average quality is always tried (`do_qa` of -1) or never (0).
    selects_by_average: bool,
}

impl Strategy {
    fn preset(self) -> Preset {
        let (q, p, d, loc, selects_by_average) = match self {
            Self::Basic => ((10, 5), (4, None), (2, 1), (0, 14, 10, 14), true),
            Self::HiSeq => ((8, 5), (7, Some(0)), (0, 0), (0, 14, 8, 14), true),
            Self::MiSeq => ((12, 6), (2, Some(0)), (2, 3), (0, 9, 12, 14), false),
            Self::IonTorrent => ((12, 6), (0, Some(0)), (0, 0), (0, 12, 0, 0), false),
        };

        Preset {
            q_bits: q.0,
            q_shift: q.1,
            p_bits: p.0,
            p_shift: p.1,
            d_bits: d.0,
            d_shift: d.1,
            q_loc: loc.0,
            s_loc: loc.1,
            p_loc: loc.2,
            d_loc: loc.3,
            selects_by_average,
        }
    }
}

const POSITIONS_TABLE_SIZE: usize = 1024;
const DELTAS_TABLE_SIZE: usize = 256;

/// The context mask: there are 2^16 quality models.
const CONTEXT_MASK: u32 = 0xffff;

/// The parameters of a block, as htscodecs' `fqz_pick_parameters` chooses them, and each
/// record's selector.
struct Parameters {
    gflags: parameters::Flags,
    /// The largest selector, 0 with none.
    max_sel: u8,
    /// One past the largest symbol the quality models code.
    symbol_count: NonZero<usize>,
    flags: parameter::Flags,
    /// The largest symbol stored, or the number of symbols with a quality map.
    max_sym: u8,
    q_bits: u32,
    q_shift: u32,
    q_loc: u32,
    s_loc: u32,
    p_loc: u32,
    d_loc: u32,
    /// The symbol each score is coded as: its rank among the scores present with a quality
    /// map, else itself.
    q_map: [u8; 256],
    /// The scores present, in order, with a quality map.
    q_map_symbols: Option<Vec<u8>>,
    p_tab: Option<[u8; POSITIONS_TABLE_SIZE]>,
    d_tab: Option<[u8; DELTAS_TABLE_SIZE]>,
    /// Each record's selector, with [`parameter::Flags::DO_SEL`].
    selectors: Vec<u8>,
}

pub(super) fn encode_with(strategy: Strategy, lens: &[usize], src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    let len =
        u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(&mut dst, len)?;

    let lens = fit_lens(lens, src.len())?;
    let parameters = pick_parameters(strategy, &lens, src);
    write_parameters(&mut dst, &parameters)?;

    let selector_count = NonZero::new(usize::from(parameters.max_sel)).map(|n| n.saturating_add(1));
    let mut models = Models::new(parameters.symbol_count, selector_count);
    let mut range_coder = RangeCoder::default();

    let has_selector = parameters.flags.has_selector();
    let has_duplicates = parameters.flags.has_duplicates();
    let is_fixed_length = parameters.flags.is_fixed_length();

    let q_mask = (1 << parameters.q_bits) - 1;

    // The tables as the context takes them: shifted into place, and zero if absent.
    let p_tab: [u32; POSITIONS_TABLE_SIZE] = std::array::from_fn(|i| {
        parameters
            .p_tab
            .map_or(0, |tab| u32::from(tab[i]) << parameters.p_loc)
    });
    let d_tab: [u32; DELTAS_TABLE_SIZE] = std::array::from_fn(|i| {
        parameters
            .d_tab
            .map_or(0, |tab| u32::from(tab[i]) << parameters.d_loc)
    });

    let mut i = 0;
    let mut last_len = 0;

    for (rec, &len) in lens.iter().enumerate() {
        if len == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "a record has no quality scores",
            ));
        }

        let selector = if has_selector {
            let selector = parameters.selectors[rec];

            if let Some(model) = models.sel.as_mut() {
                model.encode(&mut dst, &mut range_coder, selector)?;
            }

            u32::from(selector)
        } else {
            0
        };

        if !is_fixed_length || rec == 0 {
            encode_length(&mut dst, &mut range_coder, &mut models, len)?;
        }

        let scores = &src[i..i + len];

        if has_duplicates {
            let is_duplicate = i > 0 && len == last_len && src[i - len..i] == *scores;
            models
                .dup
                .encode(&mut dst, &mut range_coder, u8::from(is_duplicate))?;

            last_len = len;

            if is_duplicate {
                i += len;
                continue;
            }
        }

        let mut ctx = 0;
        let mut q_ctx: u32 = 0;
        let mut delta = 0;
        let mut prev_q = 0;

        for (&score, pos) in scores.iter().zip((1..=len).rev()) {
            let q = parameters.q_map[usize::from(score)];
            models
                .qual
                .encode(ctx as usize, &mut dst, &mut range_coder, q)?;

            // The qualities table is the identity: htscodecs never writes another.
            q_ctx = (q_ctx << parameters.q_shift).wrapping_add(u32::from(q));

            ctx = (q_ctx & q_mask) << parameters.q_loc;
            ctx += p_tab[pos.min(POSITIONS_TABLE_SIZE - 1)];
            ctx += d_tab[delta.min(DELTAS_TABLE_SIZE - 1)];
            ctx += selector << parameters.s_loc;
            ctx &= CONTEXT_MASK;

            delta += usize::from(prev_q != q);
            prev_q = q;
        }

        i += len;
    }

    range_coder.range_encode_end(&mut dst)?;

    Ok(dst)
}

/// The lengths of the records, which must add up to `len`: as htscodecs does, a record that
/// runs past the end is cut short, and the last takes what is left over.
fn fit_lens(lens: &[usize], len: usize) -> io::Result<Vec<usize>> {
    if lens.is_empty() && len > 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "quality scores without records",
        ));
    }

    let mut lens = lens.to_vec();
    let mut total = 0;

    for n in &mut lens {
        *n = (*n).min(len - total);
        total += *n;
    }

    if let Some(last) = lens.last_mut() {
        *last += len - total;
    }

    Ok(lens)
}

/// htscodecs' `fast_log`: a natural logarithm from the bits of a double, close enough to rank
/// entropies by.
fn fast_log(a: f64) -> f64 {
    (a.to_bits() as i64).wrapping_sub(4_606_921_278_410_026_770) as f64 * 1.539_095_918_623_324e-16
}

/// What the scores of a block say about how to code them (htscodecs' `fqz_qual_stats`).
struct Stats {
    /// The scores present.
    is_present: [bool; 256],
    /// Whether enough records repeat the one before them to say so for each.
    has_duplicates: bool,
    /// Each record's selector by its average quality, and the largest, if that codes smaller.
    selectors: Option<(Vec<u8>, u8)>,
}

/// The positions in a record that statistics are kept apart for, by how far from its end.
const POSITION_BINS: usize = 32;

/// Ten times the largest average quality told apart.
const AVERAGE_BINS: usize = 2560;

fn gather_stats(lens: &[usize], src: &[u8], tries_selector: bool) -> Stats {
    let mut is_present = [false; 256];
    let mut symbol_count = 0;

    // Ten times each record's average quality, rounded, and how many records have each.
    let mut averages = Vec::with_capacity(lens.len());
    let mut average_counts = [0u32; AVERAGE_BINS];

    let mut duplicates = 0usize;
    let mut i = 0;
    let mut last_len = 0;

    for &len in lens {
        let scores = &src[i..i + len];

        if i > 0 && len == last_len && src[i - len..i] == *scores {
            duplicates += 1;
        }

        let mut total = 0u32;

        for &q in scores {
            total = total.wrapping_add(u32::from(q));

            if !is_present[usize::from(q)] {
                is_present[usize::from(q)] = true;
                symbol_count += 1;
            }
        }

        let average = if len == 0 {
            0
        } else {
            (f64::from(total) * 10.0 / len as f64 + 0.5) as u32
        };

        averages.push(average);
        average_counts[(average as usize).min(AVERAGE_BINS - 1)] += 1;

        last_len = len;
        i += len;
    }

    let has_duplicates = (lens.len() + 1) / (duplicates + 1) < 500;

    let selectors = tries_selector
        .then(|| select_by_average(lens, src, symbol_count, &averages, average_counts))
        .flatten();

    Stats {
        is_present,
        has_duplicates,
        selectors,
    }
}

/// Splits the records into four classes by average quality, and returns each record's class
/// (or pair of classes) as its selector if the scores have less entropy by class than
/// together, after what the selectors cost.
fn select_by_average(
    lens: &[usize],
    src: &[u8],
    symbol_count: usize,
    averages: &[u32],
    average_counts: [u32; AVERAGE_BINS],
) -> Option<(Vec<u8>, u8)> {
    // With few symbols the scores compress well and a selector costs more of what is left, so
    // the classes are not of equal size.
    let quantiles = if symbol_count > 8 {
        [0.2, 0.5, 0.8]
    } else {
        [0.05, 0.22, 0.60]
    };

    let record_count = lens.len();

    // The class of each average: counts are turned into classes up to each quantile.
    let mut classes = [0u8; AVERAGE_BINS];
    let mut total = 0u32;
    let mut bin = 0;

    for (class, quantile) in (0..).zip(quantiles) {
        while bin < AVERAGE_BINS {
            total += average_counts[bin];

            if f64::from(total) > quantile * record_count as f64 {
                break;
            }

            classes[bin] = class;
            bin += 1;
        }
    }

    classes[bin..].fill(3);

    let class_of = |average: u32| classes[(average as usize).min(AVERAGE_BINS - 1)];

    // Counts of each score at each position, of every eighth record: all together, by pair of
    // classes, and by class.
    let mut scores_1 = vec![[0u32; 256]; POSITION_BINS];
    let mut scores_2 = vec![[[0u32; 256]; POSITION_BINS]; 2];
    let mut scores_4 = vec![[[0u32; 256]; POSITION_BINS]; 4];
    let mut totals_1 = [0u32; POSITION_BINS];
    let mut totals_2 = [[0u32; POSITION_BINS]; 2];
    let mut totals_4 = [[0u32; POSITION_BINS]; 4];

    let mut i = 0;

    for (rec, (&len, &average)) in lens.iter().zip(averages).enumerate() {
        if rec % 8 == 0 {
            let class_4 = usize::from(class_of(average));
            let class_2 = class_4 / 2;

            for (&q, pos) in src[i..i + len].iter().zip((1..=len).rev()) {
                let (x, q) = (pos % POSITION_BINS, usize::from(q));

                scores_4[class_4][x][q] += 1;
                totals_4[class_4][x] += 1;
                scores_2[class_2][x][q] += 1;
                totals_2[class_2][x] += 1;
                scores_1[x][q] += 1;
                totals_1[x] += 1;
            }
        }

        i += len;
    }

    // The sums are taken in htscodecs' order, so that they round as its do.
    let term = |count: u32, total: u32| {
        if count == 0 {
            0.0
        } else {
            f64::from(count) * fast_log(f64::from(count) / f64::from(total))
        }
    };

    let (mut e1, mut e2, mut e4) = (0.0, 0.0, 0.0);

    for x in 0..POSITION_BINS {
        for q in 0..256 {
            e1 += term(scores_1[x][q], totals_1[x]);

            for class in 0..2 {
                e2 += term(scores_2[class][x][q], totals_2[class][x]);
            }

            for class in 0..4 {
                e4 += term(scores_4[class][x][q], totals_4[class][x]);
            }
        }
    }

    // In bytes.
    let scale = -std::f64::consts::LN_2 / 8.0;
    let (e1, e2, e4) = (e1 / scale, e2 / scale, e4 / scale);

    // A selector takes context bits from the rest, by how much is a guess: it must save 2%.
    const MARGIN: f64 = 0.98;

    let cost_4 = (record_count / 4) as f64;
    let cost_2 = (record_count / 8) as f64;

    let selectors: Vec<u8> = if e4 + cost_4 < e2 * MARGIN + cost_2 && e4 + cost_4 < e1 * MARGIN {
        averages.iter().map(|&average| class_of(average)).collect()
    } else if e2 + cost_2 < e1 * MARGIN {
        averages
            .iter()
            .map(|&average| class_of(average) >> 1)
            .collect()
    } else {
        return None;
    };

    // The largest selector in use is stored, and sizes the selector model. With every record
    // in class 0 there is nothing to select.
    let max_sel = selectors.iter().copied().max().filter(|&max| max > 0)?;

    Some((selectors, max_sel))
}

/// Chooses a block's parameters as htscodecs' `fqz_pick_parameters` does.
fn pick_parameters(strategy: Strategy, lens: &[usize], src: &[u8]) -> Parameters {
    // About the square root of a count of changes.
    const DELTA_SQRT: [u8; 64] = [
        0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 7,
    ];

    /// Blocks shorter than this have fewer scores than a wide context can learn from.
    const SMALL_BLOCK_LEN: usize = 5_000_000;
    const TINY_BLOCK_LEN: usize = 300_000;

    let Preset {
        mut q_bits,
        mut q_shift,
        mut p_bits,
        mut p_shift,
        mut d_bits,
        d_shift,
        q_loc,
        mut s_loc,
        mut p_loc,
        mut d_loc,
        selects_by_average,
    } = strategy.preset();

    let stats = gather_stats(lens, src, selects_by_average);

    // A selector's bits come out of the positions and deltas, or the deltas, or the scores,
    // whether or not the selector is used.
    if selects_by_average {
        if p_bits > 0 && d_bits > 0 {
            s_loc = d_loc - 1;
            p_bits -= 1;
            d_bits -= 1;
            d_loc += 1;
        } else if d_bits >= 2 {
            s_loc = d_loc;
            d_bits -= 2;
            d_loc += 2;
        } else if q_bits >= 2 {
            q_bits -= 2;
            p_loc -= 2;
            // htscodecs takes a further bit here for a selector by read number, which a
            // strategy whose reads are told apart sets.
            s_loc = 16 - 2 - u32::from(strategy == Strategy::HiSeq);

            if q_bits == 6 && q_shift == 5 {
                q_bits -= 1;
            }
        }
    }

    let symbol_count = stats.is_present.iter().filter(|&&p| p).count();
    let max_symbol = stats.is_present.iter().rposition(|&p| p).unwrap_or(0);

    let stores_q_map = symbol_count <= 8 && symbol_count * 2 < max_symbol;
    let is_fixed_length = lens.windows(2).all(|w| w[0] == w[1]);

    let p_shift = *p_shift.get_or_insert_with(|| {
        let first_len = lens.first().copied().unwrap_or(0);
        let shift =
            (first_len as f64 / f64::from(1u32 << p_bits)).ln() / std::f64::consts::LN_2 + 0.5;

        if 0.0 > shift { 0 } else { shift as u32 }
    });
    let mut p_shift = p_shift;

    if symbol_count <= 4 {
        // Like the NovaSeq's.
        q_shift = 2;

        if src.len() < SMALL_BLOCK_LEN {
            p_bits = 2;
            p_shift = 5;
        }
    } else if symbol_count <= 8 {
        // Like the HiSeq X's.
        q_bits = q_bits.min(9);
        q_shift = 3;

        if src.len() < SMALL_BLOCK_LEN {
            q_bits = 6;
        }
    }

    if src.len() < TINY_BLOCK_LEN {
        q_bits = q_shift;
        d_bits = 2;
    }

    let mut q_map: [u8; 256] = std::array::from_fn(|i| i as u8);
    let mut q_map_symbols = None;
    let mut max_sym = max_symbol as u8;

    if stores_q_map {
        let symbols: Vec<u8> = (0..=u8::MAX)
            .filter(|&q| stats.is_present[usize::from(q)])
            .collect();

        for (rank, &q) in (0..).zip(&symbols) {
            q_map[usize::from(q)] = rank;
        }

        // At most 8.
        max_sym = symbols.len() as u8;
        q_map_symbols = Some(symbols);
    }

    let p_tab = (p_bits > 0)
        .then(|| std::array::from_fn(|i| ((1usize << p_bits) - 1).min(i >> p_shift.min(31)) as u8));

    let d_tab = (d_bits > 0).then(|| {
        let max = (1u8 << d_bits) - 1;
        std::array::from_fn(|i| DELTA_SQRT[(i >> d_shift).min(DELTA_SQRT.len() - 1)].min(max))
    });

    let mut flags = parameter::Flags::empty();
    flags.set(parameter::Flags::HAVE_DTAB, d_tab.is_some());
    flags.set(parameter::Flags::HAVE_PTAB, p_tab.is_some());
    flags.set(parameter::Flags::DO_SEL, stats.selectors.is_some());
    flags.set(parameter::Flags::DO_LEN, is_fixed_length);
    flags.set(parameter::Flags::DO_DEDUP, stats.has_duplicates);
    flags.set(parameter::Flags::HAVE_QMAP, stores_q_map);

    let mut gflags = parameters::Flags::empty();
    let (selectors, max_sel) = stats.selectors.unwrap_or_default();

    // The selectors pick no parameters, only context bits, so their table is all zeros.
    gflags.set(parameters::Flags::HAVE_S_TAB, flags.has_selector());

    Parameters {
        gflags,
        max_sel,
        symbol_count: NonZero::new(usize::from(max_sym) + 1).expect("at least one symbol"),
        flags,
        max_sym,
        q_bits,
        q_shift,
        q_loc,
        s_loc,
        p_loc,
        d_loc,
        q_map,
        q_map_symbols,
        p_tab,
        d_tab,
        selectors,
    }
}

/// Writes the parameters as htscodecs' `fqz_store_parameters` does.
fn write_parameters(dst: &mut Vec<u8>, parameters: &Parameters) -> io::Result<()> {
    const VERSION: u8 = 5;

    /// A parameter set's starting context.
    const CONTEXT: u16 = 0;

    write_u8(dst, VERSION)?;
    write_u8(dst, u8::from(parameters.gflags))?;

    if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
        write_u8(dst, parameters.max_sel)?;
        write_array(dst, &[0; 256])?;
    }

    write_u16_le(dst, CONTEXT)?;
    write_u8(dst, u8::from(parameters.flags))?;
    write_u8(dst, parameters.max_sym)?;

    // Each is below 16.
    let nibbles = |high: u32, low: u32| ((high << 4) | low) as u8;
    write_u8(dst, nibbles(parameters.q_bits, parameters.q_shift))?;
    write_u8(dst, nibbles(parameters.q_loc, parameters.s_loc))?;
    write_u8(dst, nibbles(parameters.p_loc, parameters.d_loc))?;

    if let Some(symbols) = &parameters.q_map_symbols {
        dst.write_all(symbols)?;
    }

    if let Some(tab) = &parameters.p_tab {
        write_array(dst, tab)?;
    }

    if let Some(tab) = &parameters.d_tab {
        write_array(dst, tab)?;
    }

    Ok(())
}

fn write_array<W>(writer: &mut W, data: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let mut rle1 = Vec::new();

    let mut i = 0;
    let mut j = 0;

    while j < data.len() {
        let start = j;

        while j < data.len() && usize::from(data[j]) == i {
            j += 1;
        }

        let mut len = j - start;

        loop {
            let rle = len.min(255);
            rle1.push(rle as u8);

            len -= rle;

            if rle != 255 {
                break;
            }
        }

        i += 1;
    }

    let mut rle2 = Vec::new();
    j = 0;
    let mut last = -1;

    while j < rle1.len() {
        let curr = rle1[j];
        j += 1;

        rle2.push(curr);

        if i32::from(curr) == last {
            let start = j;
            let mut len = 0;

            while j < rle1.len() && i32::from(rle1[j]) == last && len < 255 {
                j += 1;
                len = j - start;
            }

            rle2.push(len as u8);
        } else {
            last = i32::from(curr);
        }
    }

    writer.write_all(&rle2)?;

    Ok(())
}

fn encode_length<W>(
    writer: &mut W,
    range_coder: &mut RangeCoder,
    models: &mut Models,
    len: usize,
) -> io::Result<()>
where
    W: Write,
{
    let n = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    models.len[0].encode(writer, range_coder, (n & 0xff) as u8)?;
    models.len[1].encode(writer, range_coder, ((n >> 8) & 0xff) as u8)?;
    models.len[2].encode(writer, range_coder, ((n >> 16) & 0xff) as u8)?;
    models.len[3].encode(writer, range_coder, ((n >> 24) & 0xff) as u8)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cram::codecs::fqzcomp::decode;

    const STRATEGIES: [Strategy; 4] = [
        Strategy::Basic,
        Strategy::HiSeq,
        Strategy::MiSeq,
        Strategy::IonTorrent,
    ];

    /// The records of a fixture: a line of scores each, as ASCII from `!`.
    fn records(text: &[u8]) -> (Vec<usize>, Vec<u8>) {
        let lines = text.split(|&b| b == b'\n').filter(|line| !line.is_empty());

        let mut lens = Vec::new();
        let mut scores = Vec::new();

        for line in lines {
            lens.push(line.len());
            scores.extend(line.iter().map(|b| b - b'!'));
        }

        (lens, scores)
    }

    // The expected bytes are htscodecs 1.6.7's, from `fqz_compress(4, ...)` with strategy 0.
    #[test]
    fn test_encode() -> io::Result<()> {
        let data = [
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0],
        ];

        let lens: Vec<_> = data.iter().map(|scores| scores.len()).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let actual = encode(&lens, &src)?;

        let expected = [
            0x19, 0x05, 0x02, 0x03, 0xff, 0x01, 0x00, 0x00, 0x6a, 0x03, 0x22, 0x0d, 0xaf, 0x20,
            0x20, 0x01, 0xff, 0xff, 0x01, 0xa3, 0x02, 0x06, 0x0a, 0xee, 0x00, 0xc2, 0x7f, 0xff,
            0xf3, 0x01, 0x94, 0xbb, 0x08, 0xe6, 0x4a, 0xff, 0xc6, 0x24, 0x39, 0x66, 0x70, 0xc7,
            0xae,
        ];

        assert_eq!(actual, expected);
        assert_eq!(decode(&actual)?, src);

        Ok(())
    }

    #[test]
    fn test_encode_with_do_len() -> io::Result<()> {
        let data = [
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0, 0, 0, 0, 1, 1],
        ];

        let lens: Vec<_> = data.iter().map(|scores| scores.len()).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let actual = encode(&lens, &src)?;

        let expected = [
            0x1e, 0x05, 0x02, 0x03, 0xff, 0x01, 0x00, 0x00, 0x6e, 0x03, 0x22, 0x0d, 0xaf, 0x20,
            0x20, 0x01, 0xff, 0xff, 0x01, 0xa3, 0x02, 0x06, 0x0a, 0xee, 0x00, 0xc2, 0x7f, 0xff,
            0xf3, 0x01, 0x94, 0xbb, 0x0d, 0xa8, 0xd1, 0x27, 0x83, 0xbb, 0x79, 0x00,
        ];

        assert_eq!(actual, expected);
        assert_eq!(decode(&actual)?, src);

        Ok(())
    }

    /// The output is htscodecs' byte for byte, which chooses the same parameters and selectors
    /// from the scores, at each of its strategies.
    #[test]
    fn test_encode_matches_htscodecs() -> io::Result<()> {
        fn check(name: &str, text: &[u8], outputs: [&[u8]; 4]) -> io::Result<()> {
            let (lens, src) = records(text);

            for (strategy, output) in STRATEGIES.into_iter().zip(outputs) {
                let actual = encode_with(strategy, &lens, &src)?;
                assert_eq!(actual, output, "{name} records, {strategy:?}");
                assert_eq!(decode(&actual)?, src, "{name} records, {strategy:?}");
            }

            Ok(())
        }

        check(
            "mixed",
            include_bytes!("../fixtures/qual-records-mixed.txt"),
            [
                include_bytes!("../fixtures/fqz-mixed-0.bin"),
                include_bytes!("../fixtures/fqz-mixed-1.bin"),
                include_bytes!("../fixtures/fqz-mixed-2.bin"),
                include_bytes!("../fixtures/fqz-mixed-3.bin"),
            ],
        )?;

        check(
            "binned",
            include_bytes!("../fixtures/qual-records-binned.txt"),
            [
                include_bytes!("../fixtures/fqz-binned-0.bin"),
                include_bytes!("../fixtures/fqz-binned-1.bin"),
                include_bytes!("../fixtures/fqz-binned-2.bin"),
                include_bytes!("../fixtures/fqz-binned-3.bin"),
            ],
        )
    }

    /// Whatever the parameters chosen, the scores decode: with many values or few, records
    /// of one length or several, repeated records, and averages that differ enough for a
    /// selector.
    #[test]
    fn test_encode_round_trips() -> io::Result<()> {
        let mut x = 0x9e37_79b9_u32;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            x
        };

        let mut blocks: Vec<(Vec<usize>, Vec<u8>)> = vec![
            (Vec::new(), Vec::new()),
            (vec![1], vec![40]),
            (vec![3, 3], vec![255, 0, 255, 255, 0, 255]),
        ];

        for (record_count, symbols, fixed_len) in [
            (400, &[2, 11, 25, 37][..], Some(150)),
            (400, &[0, 1, 2, 3, 4, 5, 6][..], None),
            (300, &(0..94).collect::<Vec<u8>>()[..], None),
            (3000, &(0..42).collect::<Vec<u8>>()[..], Some(100)),
        ] {
            let mut lens = Vec::new();
            let mut src = Vec::new();

            for _ in 0..record_count {
                let len = fixed_len.unwrap_or_else(|| 20 + next() as usize % 300);

                if !lens.is_empty() && lens.last() == Some(&len) && next() % 16 == 0 {
                    let start = src.len() - len;
                    src.extend_from_within(start..);
                } else {
                    // Good reads and poor ones.
                    let bias = if next() % 3 == 0 {
                        0
                    } else {
                        symbols.len() / 2
                    };

                    for _ in 0..len {
                        let i = bias + next() as usize % (symbols.len() - bias);
                        src.push(symbols[i]);
                    }
                }

                lens.push(len);
            }

            blocks.push((lens, src));
        }

        for (lens, src) in &blocks {
            for strategy in STRATEGIES {
                let encoded = encode_with(strategy, lens, src)?;
                assert_eq!(
                    decode(&encoded)?,
                    *src,
                    "{} records, {strategy:?}",
                    lens.len()
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_fit_lens() -> io::Result<()> {
        assert_eq!(fit_lens(&[3, 4], 7)?, [3, 4]);
        // The last record takes what is left over, and one past the end is cut short.
        assert_eq!(fit_lens(&[3, 4], 9)?, [3, 6]);
        assert_eq!(fit_lens(&[3, 4, 5], 5)?, [3, 2, 0]);
        assert_eq!(fit_lens(&[], 0)?, []);
        assert!(fit_lens(&[], 1).is_err());

        Ok(())
    }

    #[test]
    fn test_encode_with_an_empty_record() {
        assert!(encode(&[2, 0, 2], &[1, 2, 3, 4]).is_err());
    }
}
