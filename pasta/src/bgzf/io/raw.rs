//! Raw BGZF frames, for tools that copy compressed blocks without recompressing them (bcftools
//! `reheader` and `concat --naive`, through HTSlib's `bgzf_read_block` and `bgzf_raw_read`) and
//! for tools that choose block boundaries themselves (`bgzip`).

use std::io::{self, Read};

use super::{
    reader::frame,
    writer::{BGZF_EOF, CompressionLevel, write_frame},
};

/// The BGZF end-of-file marker: an empty block.
pub const EOF_MARKER: [u8; 28] = BGZF_EOF;

/// Reads the next complete frame (header, compressed data, and trailer) into `buf`. Returns
/// `false` at the end of the stream.
///
/// # Errors
///
/// Returns an error for a truncated or malformed frame header.
pub fn read_frame<R: Read>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<bool> {
    Ok(frame::read_frame_into(reader, buf)?.is_some())
}

/// Decompresses a complete frame into `out`, replacing its contents.
///
/// # Errors
///
/// Returns an error for a malformed frame or a checksum mismatch.
pub fn inflate_frame(frame: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    frame::inflate_frame_into(frame, out)
}

/// Decompresses complete frames, each into its `out`, replacing its contents: the result of
/// [`inflate_frame`] for each, decompressed together, which on 64-bit Arm overlaps two frames'
/// decoding on one thread.
///
/// # Panics
///
/// If `frames` and `outs` differ in length.
pub fn inflate_frames(frames: &[&[u8]], outs: &mut [&mut Vec<u8>]) -> Vec<io::Result<()>> {
    assert_eq!(frames.len(), outs.len(), "one output per frame");
    frame::inflate_frames_into(frames, outs)
}

/// Compresses `data` into a complete frame in `out`, replacing its contents.
///
/// # Errors
///
/// Returns an error if the frame would exceed the BGZF maximum block size of 64 KiB (`data` is
/// at most 65280 bytes in HTSlib's blocks, which always fit).
pub fn deflate_frame(data: &[u8], level: CompressionLevel, out: &mut Vec<u8>) -> io::Result<()> {
    use crate::bgzf::deflate;

    let mut compressed = Vec::new();
    let crc32 = deflate::encode(data, level, &mut compressed)?;

    out.clear();
    write_frame(out, &compressed, crc32, data.len())?;

    Ok(())
}

/// Returns whether a frame is the BGZF end-of-file marker block.
#[must_use]
pub fn is_eof_marker(frame: &[u8]) -> bool {
    frame == BGZF_EOF
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    #[test]
    fn frames_round_trip() -> io::Result<()> {
        let mut w = crate::bgzf::io::Writer::new(Vec::new());
        w.write_all(b"noodles")?;
        let data = w.finish()?;
        let mut src = &data[..];
        let mut frame = Vec::new();
        assert!(read_frame(&mut src, &mut frame)?);
        let mut out = Vec::new();
        inflate_frame(&frame, &mut out)?;
        assert_eq!(out, b"noodles");
        assert!(read_frame(&mut src, &mut frame)?);
        assert!(is_eof_marker(&frame));
        assert!(!read_frame(&mut src, &mut frame)?);
        Ok(())
    }

    #[test]
    fn deflated_frames_inflate() -> io::Result<()> {
        let data: Vec<u8> = (0..65280u32).map(|i| (i * 7 % 251) as u8).collect();
        let mut frame = Vec::new();
        let mut out = Vec::new();

        for level in [CompressionLevel::NONE, CompressionLevel::default()] {
            deflate_frame(&data, level, &mut frame)?;
            let mut src = &frame[..];
            let mut read = Vec::new();
            assert!(read_frame(&mut src, &mut read)?);
            assert_eq!(read, frame);
            inflate_frame(&frame, &mut out)?;
            assert_eq!(out, data);
        }

        deflate_frame(b"", CompressionLevel::default(), &mut frame)?;
        inflate_frame(&frame, &mut out)?;
        assert!(out.is_empty());
        assert!(is_eof_marker(&EOF_MARKER));

        Ok(())
    }

    #[test]
    fn frames_inflate_together_as_one_by_one() -> io::Result<()> {
        let mut frames = Vec::new();
        for (i, level) in [0, 1, 5, 6, 9, 12].into_iter().enumerate() {
            let data: Vec<u8> = (0..65280u32)
                .map(|j| ((j * (i as u32 + 3)) % 61) as u8)
                .collect();
            let mut frame = Vec::new();
            deflate_frame(
                &data,
                CompressionLevel::new(level).expect("level"),
                &mut frame,
            )?;
            frames.push(frame);
        }
        frames.push(EOF_MARKER.to_vec());
        // A frame whose data no longer matches its CRC, one whose DEFLATE data is damaged, and
        // one cut short.
        let mut bad_crc = frames[2].clone();
        let at = bad_crc.len() - 8;
        bad_crc[at] ^= 1;
        let mut bad_data = frames[3].clone();
        bad_data[20] ^= 0xff;
        let short = frames[4][..frames[4].len() - 3].to_vec();
        frames.extend([bad_crc, bad_data, short, frames[1].clone()]);

        let srcs: Vec<&[u8]> = frames.iter().map(|f| &f[..]).collect();
        let mut bufs = vec![Vec::new(); srcs.len()];
        let mut outs: Vec<&mut Vec<u8>> = bufs.iter_mut().collect();
        let results = inflate_frames(&srcs, &mut outs);
        for (i, (result, src)) in results.into_iter().zip(&srcs).enumerate() {
            let mut want = Vec::new();
            match (result, inflate_frame(src, &mut want)) {
                (Ok(()), Ok(())) => assert_eq!(bufs[i], want, "frame {i}"),
                (Err(got), Err(want)) => assert_eq!(got.kind(), want.kind(), "frame {i}"),
                (got, want) => panic!("frame {i}: {got:?} against {want:?}"),
            }
        }
        Ok(())
    }
}
