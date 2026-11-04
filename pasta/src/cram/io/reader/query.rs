pub(crate) mod parallel;
mod records;
mod unmapped;

use std::{
    io::{self, Read, Seek, SeekFrom},
    mem, slice,
};

use crate::bam;
use crate::core::{Position, region::Interval};
use crate::sam::{self as sam, alignment::Record as _};

use self::records::Records;
pub use self::unmapped::UnmappedQuery;
use super::{Container, Reader};
use crate::cram::{self as cram, crai};

/// A reader over records that intersect a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h: 'r, 'i: 'r, R>
where
    R: Read + Seek,
{
    reader: &'r mut Reader<R>,

    header: &'h sam::Header,

    index: slice::Iter<'i, crai::Record>,

    reference_sequence_id: usize,
    interval: Interval,

    container: Container,
    // The offset `container` was read from, so that a container is read once for all of the
    // index records (slices) that point into it.
    container_offset: Option<u64>,

    records: Batch,
    record: bam::Record,
}

impl<'r, 'h: 'r, 'i: 'r, R> Query<'r, 'h, 'i, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'r mut Reader<R>,
        header: &'h sam::Header,
        index: &'i crai::Index,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader,

            header,

            index: index.iter(),

            reference_sequence_id,
            interval,

            container: Container::default(),
            container_offset: None,

            records: Batch::default(),
            record: bam::Record::default(),
        }
    }

    /// Reads a record as a BAM record.
    ///
    /// The BAM record is encoded from the decoded CRAM fields directly (see
    /// [`cram::Record::encode_bam_record`]), which is much cheaper than building a
    /// [`sam::alignment::RecordBuf`]. `record`'s buffer is reused.
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
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// let mut record = bam::Record::default();
    ///
    /// while query.read_record(&mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_record(&mut self, record: &mut bam::Record) -> io::Result<usize> {
        loop {
            if self.records.pop_into(record) {
                return Ok(1);
            }

            if !self.read_next_slice()? {
                return Ok(0);
            }
        }
    }

    /// Reads a record into an alignment record buffer.
    ///
    /// If successful, 1 is returned. If 0 is returned, the stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram::{self as cram, crai};
    /// use pasta::fasta;
    /// use pasta::sam;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    ///
    /// while query.read_record_buf(&mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_record_buf(&mut self, record: &mut sam::alignment::RecordBuf) -> io::Result<usize> {
        let mut raw = mem::take(&mut self.record);
        let n = self.read_record(&mut raw)?;

        if n != 0 {
            *record = sam::alignment::RecordBuf::try_from_alignment_record(self.header, &raw)?;
        }

        self.record = raw;

        Ok(n)
    }

    /// Returns an iterator over records.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram::{self as cram, crai};
    /// use pasta::fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn records(self) -> Records<'r, 'h, 'i, R> {
        Records::new(self)
    }

    /// Calls `f` with each record that intersects the region, in file order, without encoding
    /// records as BAM.
    ///
    /// Records are decoded with the reader's required fields
    /// ([`Reader::set_required_fields`]), so this suits records with only some fields decoded,
    /// whose BAM encoding would be incomplete.
    ///
    /// # Errors
    ///
    /// Returns an error if a container cannot be read or decoded, or the first error `f`
    /// returns.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use pasta::cram::{self as cram, crai};
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// let mut n = 0;
    /// query.for_each_record(|_| {
    ///     n += 1;
    ///     Ok::<_, io::Error>(())
    /// })?;
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn for_each_record<F, E>(&mut self, mut f: F) -> Result<(), E>
    where
        F: FnMut(&cram::Record<'_>) -> Result<(), E>,
        E: From<io::Error>,
    {
        let (reference_sequence_id, interval) = (self.reference_sequence_id, self.interval);

        loop {
            let Some(index_record) = self.index.next() else {
                return Ok(());
            };

            if !index_record_intersects(index_record, reference_sequence_id, interval) {
                continue;
            }

            if !self.read_container_of(index_record)? {
                return Ok(());
            }

            let compression_header = self.container.compression_header()?;
            let slice = self.container.slice_at_landmark(index_record.landmark())?;
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

            for record in &records {
                if record_intersects(record, reference_sequence_id, interval) {
                    f(record)?;
                }
            }
        }
    }

    /// Sets the query aside, releasing the reader, header, and index: [`Self::resume`] carries it
    /// on with the same ones, the reader not read or moved in between. This lets a query be read
    /// a part at a time by a type that owns the reader.
    pub fn suspend(self) -> SuspendedQuery {
        SuspendedQuery {
            index_left: self.index.len(),
            reference_sequence_id: self.reference_sequence_id,
            interval: self.interval,
            container: self.container,
            container_offset: self.container_offset,
            records: self.records,
            record: self.record,
        }
    }

    /// Carries on a query set aside by [`Self::suspend`] with the reader, header, and index it
    /// was reading.
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
    /// let region = "sq0:8-13".parse()?;
    /// let mut suspended = reader.query(&header, &index, &region)?.suspend();
    ///
    /// let mut record = bam::Record::default();
    ///
    /// loop {
    ///     let mut query = cram::io::reader::Query::resume(&mut reader, &header, &index, suspended);
    ///     if query.read_record(&mut record)? == 0 {
    ///         break;
    ///     }
    ///     suspended = query.suspend();
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn resume(
        reader: &'r mut Reader<R>,
        header: &'h sam::Header,
        index: &'i crai::Index,
        suspended: SuspendedQuery,
    ) -> Self {
        let left = index.len().saturating_sub(suspended.index_left);
        Self {
            reader,
            header,
            index: index[left..].iter(),
            reference_sequence_id: suspended.reference_sequence_id,
            interval: suspended.interval,
            container: suspended.container,
            container_offset: suspended.container_offset,
            records: suspended.records,
            record: suspended.record,
        }
    }

    /// Reads the container that `index_record` points into, unless it was the last one read.
    /// Returns `false` at the end of the stream.
    fn read_container_of(&mut self, index_record: &crai::Record) -> io::Result<bool> {
        let offset = index_record.offset();

        if self.container_offset != Some(offset) {
            self.container_offset = None;
            self.reader.seek(SeekFrom::Start(offset))?;

            if self.reader.read_container(&mut self.container)? == 0 {
                return Ok(false);
            }

            self.container_offset = Some(offset);
        }

        Ok(true)
    }

    /// Decodes the slice of the next index record, keeping the records that intersect the
    /// region. Returns `false` when there are no more slices.
    fn read_next_slice(&mut self) -> io::Result<bool> {
        let Some(index_record) = self.index.next() else {
            return Ok(false);
        };

        if !index_record_intersects(index_record, self.reference_sequence_id, self.interval) {
            return Ok(true);
        }

        if !self.read_container_of(index_record)? {
            return Ok(false);
        }

        let compression_header = self.container.compression_header()?;
        let slice = self.container.slice_at_landmark(index_record.landmark())?;
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

        let (reference_sequence_id, interval) = (self.reference_sequence_id, self.interval);
        self.records.clear();

        for record in &records {
            if may_intersect(record, reference_sequence_id, interval) {
                self.records.push_if(self.header, record, |r| {
                    intersects(r, reference_sequence_id, interval)
                })?;
            }
        }

        Ok(true)
    }
}

/// A query set aside without its reader, header, and index ([`Query::suspend`]).
pub struct SuspendedQuery {
    /// The index records not yet read: the last ones.
    index_left: usize,
    reference_sequence_id: usize,
    interval: Interval,
    container: Container,
    container_offset: Option<u64>,
    records: Batch,
    record: bam::Record,
}

/// The BAM records encoded from one CRAM slice or container, with their buffers kept for reuse.
#[derive(Default)]
struct Batch {
    records: Vec<bam::Record>,
    len: usize,
    next: usize,
}

impl Batch {
    fn clear(&mut self) {
        self.len = 0;
        self.next = 0;
    }

    /// Encodes `record` as BAM and keeps it if `keep` accepts it.
    fn push_if<F>(
        &mut self,
        header: &sam::Header,
        record: &cram::Record<'_>,
        keep: F,
    ) -> io::Result<()>
    where
        F: FnOnce(&bam::Record) -> io::Result<bool>,
    {
        if self.len == self.records.len() {
            self.records.push(bam::Record::default());
        }

        let dst = &mut self.records[self.len];
        let mut buf = Vec::from(mem::take(dst));
        record.encode_bam_record(header, &mut buf)?;
        *dst = bam::Record::try_from(buf)?;

        if keep(dst)? {
            self.len += 1;
        }

        Ok(())
    }

    /// Swaps the next kept record into `record`, which hands `record`'s buffer to the batch.
    fn pop_into(&mut self, record: &mut bam::Record) -> bool {
        if self.next < self.len {
            mem::swap(record, &mut self.records[self.next]);
            self.next += 1;
            true
        } else {
            false
        }
    }
}

/// Whether a CRAM record can intersect the region, judged from its reference sequence and
/// start only, to skip encoding records that cannot. [`intersects`] decides exactly.
fn may_intersect(
    record: &cram::Record<'_>,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    if record.reference_sequence_id != Some(reference_sequence_id) {
        return false;
    }

    match (record.alignment_start, region_interval.end()) {
        (Some(start), Some(end)) => start <= end,
        (Some(_), None) => true,
        (None, _) => false,
    }
}

/// Whether an encoded record intersects the region, with its end calculated from its CIGAR, as
/// [`sam::alignment::RecordBuf::alignment_end`] does.
fn intersects(
    record: &bam::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    if record.reference_sequence_id().transpose()? != Some(reference_sequence_id) {
        return Ok(false);
    }

    if interval_is_unbounded(region_interval) {
        return Ok(true);
    }

    match (
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            Ok(region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}

/// Whether a decoded CRAM record intersects the region, as [`intersects`] decides for its BAM
/// encoding.
fn record_intersects(
    record: &cram::Record<'_>,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    if record.reference_sequence_id != Some(reference_sequence_id) {
        return false;
    }

    if interval_is_unbounded(region_interval) {
        return true;
    }

    let Some(start) = record.alignment_start else {
        return false;
    };

    // A BAM record's span is its CIGAR's, which an unmapped record lacks; without a span, the
    // record ends where it starts.
    let span = if record.bam_flags.is_unmapped() {
        0
    } else {
        cram::record::calculate_alignment_span(record.read_length, &record.features)
    };

    let end = match span {
        0 => Some(start),
        _ => Position::new(usize::from(start) + span - 1),
    };

    end.is_some_and(|end| region_interval.intersects((start..=end).into()))
}

pub(crate) fn index_record_intersects(
    record: &crai::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> bool {
    let Some(id) = record.reference_sequence_id() else {
        return false;
    };

    if id != reference_sequence_id {
        return false;
    }

    if interval_is_unbounded(region_interval) {
        true
    } else {
        let Some(start) = record.alignment_start() else {
            return false;
        };

        let span = record.alignment_span();

        if span == 0 {
            return false;
        }

        let end = start
            .checked_add(span - 1)
            .expect("attempt to add with overflow");

        let alignment_interval = (start..=end).into();
        region_interval.intersects(alignment_interval)
    }
}

fn interval_is_unbounded(interval: Interval) -> bool {
    interval.start().is_none() && interval.end().is_none()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_record_intersects() {
        let start = const { Position::new(5).unwrap() };
        let end = const { Position::new(13).unwrap() };
        let interval = Interval::from(start..=end);

        let record = crai::Record::new(Some(0), Position::new(8), 3, 0, 0, 0);
        assert!(index_record_intersects(&record, 0, interval));

        let record = crai::Record::new(Some(0), Position::new(8), 3, 0, 0, 0);
        assert!(index_record_intersects(&record, 0, Interval::from(..)));

        let record = crai::Record::new(None, None, 0, 0, 0, 0);
        assert!(!index_record_intersects(&record, 0, interval));

        let record = crai::Record::new(Some(1), None, 0, 0, 0, 0);
        assert!(!index_record_intersects(&record, 0, interval));

        let record = crai::Record::new(Some(0), None, 0, 0, 0, 0);
        assert!(!index_record_intersects(&record, 0, interval));

        let record = crai::Record::new(Some(0), None, 0, 0, 0, 0);
        assert!(!index_record_intersects(&record, 0, interval));
    }
}
