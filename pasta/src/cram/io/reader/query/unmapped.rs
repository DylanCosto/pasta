use std::{
    io::{self, Read},
    mem,
};

use crate::bam;
use crate::sam;

use super::{
    super::{Container, Reader},
    Batch,
};

/// A reader over unplaced unmapped records (`*`).
///
/// This is created by calling [`Reader::query_unmapped`]. It reads every container from the
/// first one the index lists unplaced records in to the end of the stream and, like HTSlib,
/// returns only the unplaced records: a slice with many reference sequences can also hold placed
/// records.
pub struct UnmappedQuery<'r, 'h: 'r, R>
where
    R: Read,
{
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
    container: Container,
    records: Batch,
    record: bam::Record,
    is_done: bool,
}

impl<'r, 'h: 'r, R> UnmappedQuery<'r, 'h, R>
where
    R: Read,
{
    /// `reader` must be positioned at the first container with unplaced records; when
    /// `is_done` is set, there are none.
    pub(crate) fn new(reader: &'r mut Reader<R>, header: &'h sam::Header, is_done: bool) -> Self {
        Self {
            reader,
            header,
            container: Container::default(),
            records: Batch::default(),
            record: bam::Record::default(),
            is_done,
        }
    }

    /// Reads a record as a BAM record, as [`super::Query::read_record`] does.
    ///
    /// If successful, 1 is returned. If 0 is returned, the stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::bam;
    /// use pasta::cram::{self as cram, crai};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let mut query = reader.query_unmapped(&header, &index)?;
    ///
    /// let mut record = bam::Record::default();
    ///
    /// while query.read_record(&mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut bam::Record) -> io::Result<usize> {
        loop {
            if self.records.pop_into(record) {
                return Ok(1);
            }

            if self.is_done || !self.read_next_container()? {
                self.is_done = true;
                return Ok(0);
            }
        }
    }

    /// Decodes the next container, keeping its unplaced records. Returns `false` at the end of
    /// the stream.
    fn read_next_container(&mut self) -> io::Result<bool> {
        if self.reader.read_container(&mut self.container)? == 0 {
            return Ok(false);
        }

        let compression_header = self.container.compression_header()?;
        self.records.clear();

        for result in self.container.slices() {
            let slice = result?;
            let fields = self.reader.required_fields();
            let tags = self.reader.required_tags();
            let (core_data_src, external_data_srcs) =
                slice.decode_blocks_for(&compression_header, fields, tags)?;

            let records = slice.records_with_fields(
                self.reader.reference_sequence_repository.clone(),
                self.header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                self.reader.md_nm(),
                fields,
                tags,
            )?;

            for record in records
                .iter()
                .filter(|record| record.reference_sequence_id.is_none())
            {
                self.records.push_if(self.header, record, |_| Ok(true))?;
            }
        }

        Ok(true)
    }
}

impl<R> Iterator for UnmappedQuery<'_, '_, R>
where
    R: Read,
{
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut raw = mem::take(&mut self.record);

        let result = match self.read_record(&mut raw) {
            Ok(0) => return None,
            Ok(_) => sam::alignment::RecordBuf::try_from_alignment_record(self.header, &raw),
            Err(e) => Err(e),
        };

        self.record = raw;

        Some(result)
    }
}
