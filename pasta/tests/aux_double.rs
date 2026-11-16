//! HTSlib's non-standard `d` (double) tag type, which samtools writes to BAM: read and written
//! by SAM and BAM, refused by the CRAM writer as HTSlib refuses it.

use std::io;

use pasta::bam;
use pasta::cram;
use pasta::sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::data::field::{Tag, Value},
        record_buf::data::field::Value as ValueBuf,
    },
};

const XD: Tag = Tag::new(b'X', b'D');
const HP: Tag = Tag::new(b'H', b'P');

const SAM: &[u8] = b"@HD\tVN:1.6\n\
r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXD:d:0.125\tHP:i:1\n\
r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXD:d:-1e+300\tHP:i:2\n";

fn read_sam() -> io::Result<(sam::Header, Vec<RecordBuf>)> {
    let mut reader = sam::io::Reader::new(SAM);
    let header = reader.read_header()?;
    let records = reader
        .record_bufs(&header)
        .collect::<io::Result<Vec<_>>>()?;
    Ok((header, records))
}

fn write_bam(header: &sam::Header, records: &[RecordBuf]) -> io::Result<Vec<u8>> {
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_header(header)?;

    for record in records {
        writer.write_alignment_record(header, record)?;
    }

    writer.try_finish()?;
    Ok(writer.into_inner().into_inner())
}

#[test]
fn test_sam_to_bam_to_sam() -> io::Result<()> {
    let (header, records) = read_sam()?;

    assert_eq!(records[0].data().get(&XD), Some(&ValueBuf::Double(0.125)));
    assert_eq!(records[1].data().get(&XD), Some(&ValueBuf::Double(-1e300)));

    let src = write_bam(&header, &records)?;

    let mut reader = bam::io::Reader::new(&src[..]);
    let header = reader.read_header()?;
    let mut writer = sam::io::Writer::new(Vec::new());

    for (i, result) in reader.records().enumerate() {
        let record = result?;
        let data = record.data();

        // The tag after the double is found.
        assert!(matches!(
            data.get(&HP).transpose()?,
            Some(Value::UInt8(n)) if usize::from(n) == i + 1
        ));
        assert!(matches!(data.get(&XD).transpose()?, Some(Value::Double(_))));

        // The encoded double is 8 bytes, as HTSlib stores it.
        let expected = if i == 0 { 0.125f64 } else { -1e300 };
        let mut field = b"XDd".to_vec();
        field.extend(expected.to_le_bytes());
        assert!(data.as_bytes().starts_with(&field));

        writer.write_alignment_record(&header, &record)?;
    }

    // As HTSlib prints a double: `d:%g`.
    assert_eq!(
        writer.get_ref().as_slice(),
        &b"r0\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXD:d:0.125\tHP:i:1\n\
r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXD:d:-1e+300\tHP:i:2\n"[..]
    );

    Ok(())
}

#[test]
fn test_cram_writer_refuses_a_double() -> io::Result<()> {
    let (header, records) = read_sam()?;

    let mut writer = cram::io::Writer::new(Vec::new());
    writer.write_header(&header)?;

    let result = writer
        .write_alignment_record(&header, &records[0])
        .and_then(|()| writer.finish(&header));

    assert!(matches!(
        result,
        Err(e) if e.kind() == io::ErrorKind::InvalidInput
    ));

    Ok(())
}
