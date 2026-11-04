//! `concat -n/--naive`: concatenates bgzipped VCF or BCF files by copying their compressed
//! blocks, writing only the first file's header (bcftools' `naive_concat`).

use std::{
    io::{self, Read},
    path::Path,
    time::Instant,
};

use anyhow::bail;
use pasta::{
    bgzf::io::raw,
    variant::{
        Header, HeaderRecordType,
        index::split_idx_name,
        io::{
            Compression as ReaderCompression, OutputType, ReadError, Reader, WriteOptions, Writer,
            reader::report_missing_eof,
        },
    },
};

use crate::{
    output::{self, OutputSetting},
    raw_input::{Compression, Format, FrameReader, RawInput, report_exec_format_error, strerror},
};

/// The options of a naive concatenation.
pub struct Options<'a> {
    pub fnames: &'a [String],
    pub output: Option<&'a Path>,
    pub settings: &'a [OutputSetting],
    pub explicit_output_type: bool,
    pub check_headers: bool,
    pub verbose: bool,
}

/// The BGZF end-of-file marker block, which is not copied.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0,
    0, 0, 0, 0, 0, 0, 0,
];

/// The largest BGZF block (HTSlib's `BGZF_MAX_BLOCK_SIZE`).
const MAX_BLOCK: usize = 0x10000;

/// Bytes of frames passed to the writer at once.
const COPY_BATCH: usize = 1 << 20;

/// HTSlib's BGZF error code for a failed write (`BGZF_ERR_IO`), which bcftools prints.
const BGZF_ERR_IO: i32 = 4;

/// bcftools' error for a failed `bgzf_write` or `bgzf_flush` of the first block of a file.
fn block_failure(e: io::Error) -> anyhow::Error {
    output::write_failure(e, || anyhow::anyhow!("\nError: {BGZF_ERR_IO}"))
}

/// bcftools' `naive_concat_check_headers`: the same samples, and the same IDX for every
/// FILTER, INFO, FORMAT, and contig, in every file.
fn check_headers(fnames: &[String]) -> anyhow::Result<()> {
    eprintln!("Checking the headers of {} files.", fnames.len());
    let mut hdr0: Option<Header> = None;
    for (i, f) in fnames.iter().enumerate() {
        let mut reader = match Reader::open(f, 0) {
            Ok(reader) => reader,
            Err(ReadError::UnknownBinary) => {
                report_exec_format_error(f);
                bail!("Failed to open: {f}");
            }
            Err(_) => bail!("Failed to open: {f}"),
        };
        let Ok(hdr) = reader.read_header() else {
            bail!("Failed to parse header: {f}");
        };
        let Some(h0) = &hdr0 else {
            hdr0 = Some(hdr);
            continue;
        };
        let _ = i;
        if h0.sample_count() != hdr.sample_count() {
            bail!(
                "Cannot concatenate, different number of samples: {} vs {} in {} vs {f}",
                h0.sample_count(),
                hdr.sample_count(),
                fnames[0]
            );
        }
        if h0.samples() != hdr.samples() {
            bail!(
                "Cannot concatenate, different samples in {} vs {f}",
                fnames[0]
            );
        }
        if reader.compression() != ReaderCompression::Bgzf {
            bail!("The --naive option works only for compressed BCFs or VCFs, sorry :-/");
        }
        check_hrecs(h0, &hdr, &fnames[0], f)?;
        check_hrecs(&hdr, h0, f, &fnames[0])?;
    }
    eprintln!("Done, the headers are compatible.");
    Ok(())
}

/// bcftools' `_check_hrecs`.
fn check_hrecs(hdr0: &Header, hdr: &Header, fname0: &str, fname: &str) -> anyhow::Result<()> {
    for r0 in hdr0.records() {
        let ty = match r0.kind() {
            HeaderRecordType::Filter => "FILTER",
            HeaderRecordType::Info => "INFO",
            HeaderRecordType::Format => "FORMAT",
            HeaderRecordType::Contig => "contig",
            _ => continue,
        };
        let id = r0.get(b"ID").unwrap_or_default();
        let ids = String::from_utf8_lossy(id);
        let Some(r) = hdr.get_record(r0.kind(), b"ID", Some(id), None) else {
            bail!(
                "Cannot use --naive, incompatible headers, the tag {ty}/{ids} not present in {fname}"
            );
        };
        let (Some(idx0), Some(idx)) = (r0.get(b"IDX"), r.get(b"IDX")) else {
            bail!("fixme: unexpected IDX<0 for {ty}/{ids} in {fname0} or {fname}");
        };
        if idx0 != idx {
            bail!(
                "Cannot use --naive, use --naive-force instead: different order the tag {ty}/{ids} in {fname0} vs {fname}"
            );
        }
    }
    Ok(())
}

/// Reads decompressed bytes across blocks (HTSlib's `bgzf_read` after `bgzf_read_block`).
struct BlockCursor<R: Read> {
    fr: FrameReader<R>,
    offset: usize,
}

impl<R: Read> BlockCursor<R> {
    /// Fills `out`; `false` when the input ends first.
    fn read(&mut self, out: &mut Vec<u8>, n: usize) -> io::Result<bool> {
        out.clear();
        while out.len() < n {
            if self.offset >= self.fr.data.len() {
                self.fr.read_block()?;
                self.offset = 0;
                if self.fr.data.is_empty() {
                    return Ok(false);
                }
            }
            let k = (n - out.len()).min(self.fr.data.len() - self.offset);
            out.extend_from_slice(&self.fr.data[self.offset..self.offset + k]);
            self.offset += k;
        }
        Ok(true)
    }
}

/// bcftools' `print_vcf_gz_header`: reads the header lines from the decompressed blocks and
/// returns the offset of the first record in the current block.
fn read_vcf_gz_header<R: Read>(
    fr: &mut FrameReader<R>,
    hdr: &mut Vec<u8>,
    fname: &str,
) -> anyhow::Result<usize> {
    if fr.data[0] != b'#' {
        bail!(
            "Could not parse the header, expected '#', found '{}'",
            fr.data[0] as char
        );
    }
    let at = |d: &[u8], i: usize| d.get(i).copied().unwrap_or(0);
    let next = |fr: &mut FrameReader<R>| -> anyhow::Result<()> {
        fr.read_block()
            .map_err(|_| anyhow::anyhow!("\nError reading {fname}"))
    };
    let mut nskip = 1usize;
    loop {
        if at(&fr.data, nskip) == b'\n' {
            nskip += 1;
            if nskip >= fr.data.len() {
                hdr.extend_from_slice(&fr.data[..nskip.min(fr.data.len())]);
                next(fr)?;
                if fr.data.is_empty() {
                    break;
                }
                nskip = 0;
            }
            if at(&fr.data, nskip) != b'#' {
                hdr.extend_from_slice(&fr.data[..nskip]);
                break;
            }
        }
        nskip += 1;
        if nskip >= fr.data.len() {
            hdr.extend_from_slice(&fr.data);
            next(fr)?;
            if fr.data.is_empty() {
                break;
            }
            nskip = 0;
        }
    }
    Ok(nskip)
}

/// bcftools' `naive_concat`.
pub fn naive_concat(o: &Options<'_>) -> anyhow::Result<()> {
    if o.check_headers {
        check_headers(o.fnames)?;
    }
    // The format requested with -O alone (bcftools does not look at the file name here).
    let (requested, _) = output::resolve(o.settings, None)?;
    // Without compression threads, as bcftools writes: the blocks are copied as they are, and
    // a failed write is reported (and logged by HTSlib) when it happens.
    let mut writer = Writer::create(
        o.output.map(output::data_path).as_deref(),
        &WriteOptions {
            output_type: OutputType::VcfBgzf,
            level: None,
            threads: 0,
        },
    )
    .map_err(|e| anyhow::anyhow!("Error: {}", strerror(&e)))?;
    let result = copy_files(o, &mut writer, requested);
    if let Err(e) = result {
        writer.abandon();
        return Err(e);
    }
    writer.finish().map_err(|e| {
        let message = format!("Error: {}", strerror(&e));
        output::write_failure(e, || anyhow::anyhow!(message))
    })?;
    Ok(())
}

fn copy_files(o: &Options<'_>, writer: &mut Writer, requested: OutputType) -> anyhow::Result<()> {
    let out_vcf = requested.is_vcf();
    let out_bgzf = matches!(requested, OutputType::VcfBgzf | OutputType::Bcf);
    let mut file_types = 0u8;
    let mut tmp = Vec::new();
    let mut batch: Vec<u8> = Vec::with_capacity(COPY_BATCH + 2 * MAX_BLOCK);
    let mut inflated = Vec::new();
    let out_name = o
        .output
        .map_or_else(|| "-".into(), |p| p.display().to_string());
    for (i, f) in o.fnames.iter().enumerate() {
        let t0 = Instant::now();
        if o.verbose {
            eprint!("Concatenating {f}");
        }
        // Checking the headers read the file once: a stream goes on from where that left it,
        // without waiting for the writer of a FIFO (bcftools waits forever). `hts_open` reads
        // `FILE` of `FILE##idx##INDEX`.
        let data = split_idx_name(f).0;
        let opened = if o.check_headers {
            RawInput::open_again(data)
        } else {
            RawInput::open(data)
        };
        let Ok(input) = opened else {
            bail!("\nFailed to open: {f}");
        };
        if input.is_unknown_binary() {
            report_exec_format_error(data);
            bail!("\nFailed to open: {f}");
        }
        if input.compression != Compression::Bgzf {
            bail!("\nThe --naive option works only for compressed BCFs or VCFs");
        }
        let is_vcf = input.format == Format::Vcf;
        file_types |= if is_vcf { 1 } else { 2 };
        if file_types == 3 {
            bail!(
                "\nThe --naive option works only for compressed files of the same type, all BCFs or all VCFs"
            );
        }
        if o.explicit_output_type {
            if out_vcf != is_vcf {
                bail!(
                    "\nThe --naive option works only for the output of the same type, all BCFs or all VCFs"
                );
            }
            if !out_bgzf {
                bail!(
                    "\nThe --naive option works only for the output of the same compression type"
                );
            }
        }
        let mut fr = FrameReader::new(input.source, data);
        if fr.read_block().is_err() || fr.data.is_empty() {
            bail!("\nFailed to read {f}: Success");
        }
        let nskip = if is_vcf {
            tmp.clear();
            let n = read_vcf_gz_header(&mut fr, &mut tmp, f)?;
            if i == 0 {
                writer.write_raw(&tmp).map_err(|e| {
                    let n = tmp.len();
                    output::write_failure(e, || anyhow::anyhow!("Failed to write {n} bytes"))
                })?;
            }
            n
        } else {
            let mut cur = BlockCursor { fr, offset: 0 };
            let mut magic = Vec::new();
            if !cur.read(&mut magic, 9).unwrap_or(false) {
                bail!("\nFailed to read the BCF header in {f}");
            }
            if &magic[..5] != b"BCF\x02\x02" {
                bail!("\nInvalid BCF magic string in {f}");
            }
            let l_text = u32::from_le_bytes(magic[5..9].try_into().expect("4 bytes")) as usize;
            if !cur.read(&mut tmp, l_text).unwrap_or(false) {
                bail!("\nFailed to read the BCF header in {f}");
            }
            if i == 0 {
                for part in [&magic, &tmp] {
                    writer.write_raw(part).map_err(|e| {
                        let n = part.len();
                        output::write_failure(e, || {
                            anyhow::anyhow!("\nFailed to write {n} bytes to {out_name}")
                        })
                    })?;
                }
            }
            fr = cur.fr;
            cur.offset
        };
        if fr.data.len() > nskip {
            writer.write_raw(&fr.data[nskip..]).map_err(block_failure)?;
        }
        writer.write_bgzf_frames(&[]).map_err(block_failure)?;
        // Copy the remaining blocks as they are, without the EOF markers. The last block is
        // held back until the input ends, to be checked.
        // Whether the input has ended with the EOF marker block so far: reading the header
        // may have reached its end.
        let mut ends_with_marker = fr.last_block_empty();
        let src = fr.inner_mut();
        batch.clear();
        let mut last = 0;
        loop {
            let start = batch.len();
            batch.resize(start + 18, 0);
            let n = read_up_to(src, &mut batch[start..])?;
            if n == 0 {
                batch.truncate(start);
                break;
            }
            let head = &batch[start..start + 18];
            let ok = n == 18
                && head[0] == 31
                && head[1] == 139
                && head[2] == 8
                && head[3] & 4 != 0
                && u16::from_le_bytes([head[10], head[11]]) == 6
                && head[12] == b'B'
                && head[13] == b'C'
                && u16::from_le_bytes([head[14], head[15]]) == 2;
            let nblock = usize::from(u16::from_le_bytes([head[16], head[17]])) + 1;
            let bad = if !ok {
                Some(format!("\nCould not parse the header of a bgzf block: {f}"))
            } else if !(18..=MAX_BLOCK).contains(&nblock) {
                Some(format!(
                    "\nError: malformed bgzf block in {f}: block size {nblock} out of valid range [18,{MAX_BLOCK}]"
                ))
            } else {
                batch.resize(start + nblock, 0);
                (read_up_to(src, &mut batch[start + 18..])? != nblock - 18)
                    .then(|| format!("\nCould not read {nblock} bytes: {f}"))
            };
            if let Some(msg) = bad {
                // bcftools has copied the blocks before this one.
                write_frames(writer, &batch[..start])?;
                bail!(msg);
            }
            ends_with_marker = batch[start..] == BGZF_EOF;
            if ends_with_marker {
                batch.truncate(start);
                continue;
            }
            last = start;
            if start >= COPY_BATCH {
                write_frames(writer, &batch[..start])?;
                batch.drain(..start);
                last = 0;
            }
        }
        // bcftools checks only the size of each block, so a file cut inside its last block and
        // then given an EOF marker whose bytes make up the missing length passes as a whole
        // block, and is copied corrupt. The last block is decompressed to catch that.
        if batch.len() > last && raw::inflate_frame(&batch[last..], &mut inflated).is_err() {
            write_frames(writer, &batch[..last])?;
            bail!("\nError: the last bgzf block of {f} is truncated or corrupt");
        }
        write_frames(writer, &batch)?;
        // bcftools copies a file without the EOF marker block (as a writer killed between two
        // blocks leaves it) without a word; vtools fails when it is done (README Damaged input).
        if !ends_with_marker {
            report_missing_eof(data);
        }
        if o.verbose {
            eprintln!("\t{:.6} seconds", t0.elapsed().as_secs_f64());
        }
    }
    Ok(())
}

/// Appends BGZF frames to the output. bcftools copies one frame at a time, so a failure names
/// the size of the first (`bgzf_raw_write` returning -1).
fn write_frames(writer: &mut Writer, frames: &[u8]) -> anyhow::Result<()> {
    writer.write_bgzf_frames(frames).map_err(|e| {
        let len = frames
            .get(16..18)
            .map_or(0, |n| usize::from(u16::from_le_bytes([n[0], n[1]])) + 1);
        output::write_failure(e, || {
            anyhow::anyhow!("\nWrite failed, wrote -1 instead of {len} bytes.")
        })
    })
}

fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use pasta::bgzf::io::writer::CompressionLevel;

    use super::*;
    use crate::raw_input::TestDir;

    /// A bgzipped VCF in three blocks: the header with a record, more records, and the EOF
    /// marker.
    fn vcf_gz() -> Vec<u8> {
        let mut out = Vec::new();
        let mut frame = Vec::new();
        for text in [
            &b"##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tC\t.\t.\t.\n"[..],
            b"1\t2\t.\tA\tG\t.\t.\t.\n1\t3\t.\tA\tT\t.\t.\t.\n",
        ] {
            raw::deflate_frame(text, CompressionLevel::default(), &mut frame)
                .expect("the block compresses");
            out.extend_from_slice(&frame);
        }
        out.extend_from_slice(&BGZF_EOF);
        out
    }

    fn concat(dir: &TestDir, name: &str, data: &[u8]) -> anyhow::Result<()> {
        let path = dir.join(name);
        std::fs::write(&path, data).expect("write the input");
        let fnames = [path.display().to_string()];
        let out = dir.join(format!("{name}.out"));
        naive_concat(&Options {
            fnames: &fnames,
            output: Some(&out),
            settings: &[],
            explicit_output_type: false,
            check_headers: false,
            verbose: false,
        })
    }

    #[test]
    fn a_cut_last_block_fails_even_when_an_eof_marker_fills_it_up() {
        let dir = TestDir::new("naive-cut-block");
        let good = vcf_gz();
        concat(&dir, "good.vcf.gz", &good).expect("an intact file is copied");
        // The records block cut short by 28 bytes, which the EOF marker after it makes up: the
        // block has its full size but does not decompress.
        let mut filled = good[..good.len() - 2 * BGZF_EOF.len()].to_vec();
        filled.extend_from_slice(&BGZF_EOF);
        let e = concat(&dir, "filled.vcf.gz", &filled).expect_err("the corrupt block fails");
        assert_eq!(
            e.to_string(),
            format!(
                "\nError: the last bgzf block of {} is truncated or corrupt",
                dir.join("filled.vcf.gz").display()
            )
        );
        // Cut short without the marker, the block is short of its size.
        let cut = &good[..good.len() - BGZF_EOF.len() - 10];
        let e = concat(&dir, "cut.vcf.gz", cut).expect_err("the short block fails");
        assert!(e.to_string().starts_with("\nCould not read "), "{e}");
    }

    #[test]
    fn a_broken_block_fails_after_copying_the_blocks_before_it() {
        use std::io::Read as _;

        let dir = TestDir::new("naive-broken-block");
        let good = vcf_gz();
        // The header block, then uncompressed blocks of records (about 1 MB in all, more than
        // the output buffers), then bytes that are no BGZF block.
        let mut data = good[..good.len() - BGZF_EOF.len()].to_vec();
        let mut frame = Vec::new();
        let mut text = Vec::new();
        for i in 4..40_000 {
            text.extend_from_slice(format!("1\t{i}\t.\tA\tC\t.\t.\t.\n").as_bytes());
        }
        for chunk in text.chunks(0xff00) {
            raw::deflate_frame(chunk, CompressionLevel::NONE, &mut frame)
                .expect("the block compresses");
            data.extend_from_slice(&frame);
        }
        data.extend_from_slice(&[0u8; 18]);
        let e = concat(&dir, "broken.vcf.gz", &data).expect_err("the broken block fails");
        assert!(
            e.to_string()
                .starts_with("\nCould not parse the header of a bgzf block"),
            "{e}"
        );
        // bcftools has copied the blocks before it, which reach the file but for what its
        // buffer holds.
        let out = std::fs::read(dir.join("broken.vcf.gz.out")).expect("the output is there");
        let mut copied = Vec::new();
        let _ = flate2::read::MultiGzDecoder::new(&out[..]).read_to_end(&mut copied);
        assert!(
            copied.len() > text.len() / 2,
            "{} of {} bytes of records copied",
            copied.len(),
            text.len()
        );
    }
}
