use std::io;

use super::{
    Models,
    parameters::{Parameters, fqz_decode_params},
};
use crate::cram::{codecs::aac::RangeCoder, io::reader::num::read_uint7_as};

/// The contexts: the quality models are chosen by the low 16 bits.
const CONTEXT_MASK: u32 = 0xffff;

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let uncompressed_size = read_uncompressed_size(&mut src)?;

    let params = fqz_decode_params(&mut src)?;

    let mut models = Models::new(params.max_symbol_count, params.selector_count);
    let mut range_coder = RangeCoder::new(&mut src)?;

    let mut dst = vec![0; uncompressed_size];
    let mut i = 0;

    let mut rec_no = 0;
    let mut last_len = 0;
    let mut rev_len = Vec::new();

    while i < uncompressed_size {
        let record = read_record(
            &mut src,
            &params,
            &mut range_coder,
            &mut models,
            rec_no,
            last_len,
        )?;

        rec_no += 1;
        last_len = record.len;

        if record.is_duplicate {
            let len = record.len as usize;
            copy_record(&mut dst, i, len)?;
            rev_len.push((record.is_reversed, len));
            i += len;
            continue;
        }

        let remaining = uncompressed_size - i;

        // As htscodecs does, the position counts down from the length and wraps, so a record
        // of no scores takes the rest of the block.
        let len = match record.len as usize {
            0 => remaining,
            n => n.min(remaining),
        };

        rev_len.push((record.is_reversed, len));

        let param = &params.params[record.parameter];

        let selector = if param.flags().has_selector() {
            u32::from(record.selector) << param.s_loc
        } else {
            0
        };

        let base = u32::from(param.context) + selector;

        let mut ctx = usize::from(param.context);
        let mut q_ctx: u32 = 0;
        let mut pos = record.len;
        let mut delta: u32 = 0;
        let mut prev_q = 0;

        for d in &mut dst[i..i + len] {
            let q = models.qual.decode(ctx, &mut src, &mut range_coder)?;
            let j = usize::from(q);

            *d = param.quality_map[j];

            q_ctx = (q_ctx << param.q_shift).wrapping_add(param.qualities_table[j]);

            let mut next = base + ((q_ctx & param.q_mask) << param.q_loc);
            next += param.positions_table[pos.min(1023) as usize];
            next += param.deltas_table[delta.min(255) as usize];
            ctx = (next & CONTEXT_MASK) as usize;

            delta += u32::from(prev_q != q);
            prev_q = q;
            pos = pos.wrapping_sub(1);
        }

        i += len;
    }

    if params.gflags.has_reversed_values() {
        reverse_qualities(&mut dst, &rev_len);
    }

    Ok(dst)
}

fn read_uncompressed_size(src: &mut &[u8]) -> io::Result<usize> {
    read_uint7_as(src)
}

/// What precedes a record's scores.
#[derive(Debug, Default)]
struct Record {
    /// The index of the parameter set its scores are coded with.
    parameter: usize,
    selector: u8,
    len: u32,
    is_reversed: bool,
    is_duplicate: bool,
}

fn read_record(
    src: &mut &[u8],
    parameters: &Parameters,
    range_coder: &mut RangeCoder,
    models: &mut Models,
    rec_no: usize,
    last_len: u32,
) -> io::Result<Record> {
    let mut record = Record::default();

    if let Some(model) = models.sel.as_mut() {
        record.selector = model.decode(src, range_coder)?;

        if let Some(table) = parameters.selector_table() {
            record.parameter = usize::from(table[usize::from(record.selector)]);
        }
    }

    let param = parameters.params.get(record.parameter).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "invalid fqzcomp parameter set")
    })?;

    record.len = if !param.flags().is_fixed_length() || rec_no == 0 {
        read_length(src, range_coder, models)?
    } else {
        last_len
    };

    if parameters.gflags.has_reversed_values() {
        record.is_reversed = models.rev.decode(src, range_coder).map(decode_bool)?;
    }

    if param.flags().has_duplicates() {
        record.is_duplicate = models.dup.decode(src, range_coder).map(decode_bool)?;
    }

    Ok(record)
}

fn read_length(
    src: &mut &[u8],
    range_coder: &mut RangeCoder,
    models: &mut Models,
) -> io::Result<u32> {
    let mut buf = [0; 4];

    for (d, model) in buf.iter_mut().zip(&mut models.len) {
        *d = model.decode(src, range_coder)?;
    }

    Ok(u32::from_le_bytes(buf))
}

/// Reverses the records marked so, given each record's mark and length.
fn reverse_qualities(qual: &mut [u8], rev_len: &[(bool, usize)]) {
    let mut records = qual;

    for &(is_reversed, len) in rev_len {
        let Some((record, rest)) = records.split_at_mut_checked(len) else {
            break;
        };

        if is_reversed {
            record.reverse();
        }

        records = rest;
    }
}

fn decode_bool(n: u8) -> bool {
    n != 0
}

/// Copies the `len` scores before `pos` to `pos`: a record that repeats the one before it.
fn copy_record(buf: &mut [u8], pos: usize, len: usize) -> io::Result<()> {
    let start = pos.checked_sub(len);
    let end = pos.checked_add(len).filter(|&end| end <= buf.len());

    let (Some(start), Some(_)) = (start, end) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid fqzcomp duplicate record",
        ));
    };

    buf.copy_within(start..pos, pos);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() -> io::Result<()> {
        let src = [
            0x19, 0x05, 0x00, 0x00, 0x00, 0x20, 0x03, 0x82, 0x7f, 0x0f, 0x01, 0x01, 0x7d, 0xff,
            0xff, 0x01, 0x84, 0x00, 0x09, 0xff, 0xff, 0xf6, 0x01, 0x65, 0x00, 0x86, 0x2e, 0x98,
            0xea, 0xca, 0x71, 0x6f, 0x22, 0xcd, 0xd8, 0x40,
        ];

        let actual = decode(&src)?;

        let expected = [
            0, 0, 0, 1, 1, 2, 1, 1, 0, 0, // 1
            0, 1, 2, 3, 3, 3, 3, 3, 3, 3, // 2
            2, 1, 1, 0, 0, // 3
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_reverse_qualities() {
        let mut data = Vec::from(b"ndlsndlsndls");
        reverse_qualities(&mut data, &[(false, 4), (true, 4), (false, 4)]);
        assert_eq!(data, b"ndlssldnndls");

        // A record of no scores, and one past the end.
        let mut data = Vec::from(b"ndlsndls");
        reverse_qualities(&mut data, &[(true, 0), (true, 4), (true, 5)]);
        assert_eq!(data, b"sldnndls");
    }

    #[test]
    fn test_decode_bool() {
        assert!(!decode_bool(0));
        assert!(decode_bool(1));
        assert!(decode_bool(2));
        assert!(decode_bool(255));
    }

    #[test]
    fn test_copy_record() -> io::Result<()> {
        let mut buf = *b"ndls\0\0\0\0";
        copy_record(&mut buf[..], 4, 2)?;
        assert_eq!(buf, *b"ndlsls\0\0");

        copy_record(&mut buf[..], 6, 0)?;
        assert_eq!(buf, *b"ndlsls\0\0");

        // No record before it, and no room for it.
        assert!(copy_record(&mut buf[..], 2, 3).is_err());
        assert!(copy_record(&mut buf[..], 6, 3).is_err());

        Ok(())
    }

    /// A payload that is not what an encoder wrote is an error or other scores, whatever
    /// lengths, duplicates and symbols it gives.
    #[test]
    fn test_decode_with_a_corrupt_payload() -> io::Result<()> {
        let streams: [&[u8]; 4] = [
            include_bytes!("../fixtures/fqz-mixed-0.bin"),
            include_bytes!("../fixtures/fqz-mixed-2.bin"),
            include_bytes!("../fixtures/fqz-binned-0.bin"),
            include_bytes!("../fixtures/fqz-binned-3.bin"),
        ];

        for stream in streams {
            let mut payload = stream;
            let len = read_uncompressed_size(&mut payload)?;
            fqz_decode_params(&mut payload)?;
            let start = stream.len() - payload.len();

            assert_eq!(decode(stream)?.len(), len);

            for i in (start..stream.len()).step_by(11).take(64) {
                for mask in [0x01, 0x80, 0xff] {
                    let mut corrupt = stream.to_vec();
                    corrupt[i] ^= mask;

                    if let Ok(scores) = decode(&corrupt) {
                        assert_eq!(scores.len(), len);
                    }
                }
            }

            // The payload ends early.
            for end in (start..stream.len()).step_by(97) {
                assert!(decode(&stream[..end]).is_err(), "{end} of {}", stream.len());
            }
        }

        Ok(())
    }
}
