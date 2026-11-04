//! `ntools fixmate`: fill in mate information from name-grouped alignments (samtools fixmate).
//!
//! This follows `bam_mating_core` and `bam_sanitize` in samtools 1.24 step by step, including
//! the order in which tags are replaced.

use std::{
    collections::BTreeMap,
    panic::{self, AssertUnwindSafe},
    path::PathBuf,
    thread,
};

use anyhow::{Context, bail};
use bstr::BString;
use pasta::core::Position;
use pasta::sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::data::field::Value,
    },
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};

use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, OutputFormatOptions, Verbosity};
use crate::{
    alignment::{AlignmentReader, ReadOptions},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    sanitize::Sanitize,
};

/// Base qualities below this are not counted in mate scores (`MD_MIN_QUALITY`).
const MIN_SCORED_QUALITY: u8 = 15;

/// Largest CIGAR operation length (28 bits).
const MAX_OP_LEN: usize = (1 << 28) - 1;

const MATE_CIGAR: Tag = Tag::new(b'M', b'C');
const MATE_MAPPING_QUALITY: Tag = Tag::new(b'M', b'Q');
const MATE_SCORE: Tag = Tag::new(b'm', b's');
const TEMPLATE_CIGAR: Tag = Tag::new(b'c', b't');

/// Tags removed from unmapped records when sanitizing auxiliary data.
const UNMAPPED_REMOVED_TAGS: [Tag; 4] = [
    Tag::EDIT_DISTANCE,
    Tag::MISMATCHED_POSITIONS,
    Tag::new(b'C', b'G'),
    Tag::new(b'S', b'M'),
];

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Name-grouped input SAM, BAM, or CRAM file.
    input: String,

    /// Output file.
    output: PathBuf,

    /// Remove unmapped reads and secondary alignments.
    #[arg(short = 'r')]
    remove_reads: bool,

    /// Disable the FR proper pair check.
    #[arg(short = 'p')]
    no_proper_pair_check: bool,

    /// Add template CIGAR (ct) tags.
    #[arg(short = 'c')]
    template_cigar: bool,

    /// Add mate score (ms) tags.
    #[arg(short = 'm')]
    mate_score: bool,

    /// Fix base modification tags (MM/ML/MN).
    #[arg(short = 'M')]
    base_mods: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Sanitize alignment fields [all].
    #[arg(short = 'z', long = "sanitize", default_value = "all")]
    sanitize: String,

    /// Output format [default: from the output extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(
        short = '@',
        long = "threads",
        default_value_t = 0,
        value_parser = crate::getopt::negative_as_0::<usize>
    )]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    output_format_options: OutputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

#[derive(Clone, Copy, Debug)]
struct Options {
    remove_reads: bool,
    proper_pair_check: bool,
    template_cigar: bool,
    mate_score: bool,
    base_mods: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let sanitize = Sanitize::parse(&args.sanitize)?;

    let mut read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };
    args.input_format_options.configure(&mut read_options)?;

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    args.input_format_options.apply(&mut reader)?;
    let mut header = reader.read_header()?;

    let is_coordinate_sorted = header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| so == COORDINATE);

    if is_coordinate_sorted {
        bail!("[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname.");
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let header = header;

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&args.output))
        .unwrap_or(OutputFormat::Bam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed {
            0
        } else {
            args.output_format_options.level(6)?
        })?,
        threads: args.threads,
        reference: args.reference,
        cram: args.output_format_options.cram()?,
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header(&header)?;

    let options = Options {
        remove_reads: args.remove_reads,
        proper_pair_check: !args.no_proper_pair_check,
        template_cigar: args.template_cigar,
        mate_score: args.mate_score,
        base_mods: args.base_mods,
    };

    let lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length().get())
        .collect();

    if args.threads > 0
        && matches!(reader, AlignmentReader::Bam(_))
        && writer.preformatted_kind() == Some(OutputFormat::Bam)
    {
        let config = BatchConfig {
            header: &header,
            lengths: &lengths,
            sanitize,
            options,
        };

        fix_in_parallel(&mut reader, &config, args.threads, &mut writer)?;
        writer.finish(&header)?;
        return Ok(());
    }

    let mut template: Vec<RecordBuf> = Vec::new();
    let mut framed = Vec::new();

    let flush =
        |template: &mut Vec<RecordBuf>, writer: &mut AlignmentWriter| -> anyhow::Result<()> {
            fix_template(template, options, &header)?;

            for record in template.iter() {
                let flags = record.flags();

                if options.remove_reads && flags.intersects(Flags::SECONDARY | Flags::UNMAPPED) {
                    continue;
                }

                writer.write(&header, &crate::alignment::Rec::Buf(record))?;
            }

            template.clear();
            Ok(())
        };

    reader.for_each(&header, |record| {
        let mut record = match &record {
            // The BAM decoder builds record buffers much faster than the generic conversion.
            crate::alignment::Rec::Bam(r) => {
                framed.clear();
                crate::runs::write_raw_record(&mut framed, r.as_ref())?;
                let mut buf = RecordBuf::default();
                pasta::bam::io::Reader::from(&framed[..]).read_record_buf(&header, &mut buf)?;
                buf
            }
            _ => RecordBuf::try_from_alignment_record(&header, record.as_dyn())?,
        };
        sanitize_record(&mut record, &lengths, sanitize)?;

        if template
            .first()
            .is_some_and(|first| first.name() != record.name())
        {
            flush(&mut template, &mut writer)?;
        }

        template.push(record);
        Ok(())
    })?;

    if !template.is_empty() {
        flush(&mut template, &mut writer)?;
    }

    writer.finish(&header)?;

    Ok(())
}

/// Raw record bytes per batch sent to a worker.
const BATCH_CAPACITY: usize = 1 << 20;

/// Batches in flight per worker: sent and not yet written, whether queued, being fixed, or fixed
/// and waiting for an earlier batch. This bounds memory, even behind a slow batch, while keeping
/// the workers busy.
const IN_FLIGHT_PER_WORKER: usize = 4;

/// What a worker needs to fix a batch of templates.
struct BatchConfig<'a> {
    header: &'a sam::Header,
    lengths: &'a [usize],
    sanitize: Sanitize,
    options: Options,
}

/// Fixes templates on `threads` workers: this thread batches raw records at template
/// boundaries, workers decode, fix, and re-encode them, and a writer thread writes batches in
/// input order.
fn fix_in_parallel(
    reader: &mut AlignmentReader,
    config: &BatchConfig<'_>,
    threads: usize,
    writer: &mut AlignmentWriter,
) -> anyhow::Result<()> {
    let header = config.header;
    let in_flight = threads * IN_FLIGHT_PER_WORKER;

    thread::scope(|scope| {
        // One permit per batch that may be in flight; the writer thread returns each once the
        // batch is written. Neither channel can then hold more than `in_flight` batches.
        let (permits_tx, permits_rx) = crossbeam_channel::bounded(in_flight);

        for _ in 0..in_flight {
            permits_tx
                .send(())
                .expect("the permit channel has room for every permit");
        }

        let (work_tx, work_rx) = crossbeam_channel::bounded::<(u64, Vec<u8>)>(in_flight);
        let (done_tx, done_rx) =
            crossbeam_channel::bounded::<(u64, anyhow::Result<Vec<u8>>)>(in_flight);

        for _ in 0..threads {
            let work_rx = work_rx.clone();
            let done_tx = done_tx.clone();

            scope.spawn(move || {
                for (id, batch) in work_rx {
                    // A panic still sends a result: a batch that never arrives would leave the
                    // writer waiting for it, and this thread's reader for permits, forever.
                    let result =
                        panic::catch_unwind(AssertUnwindSafe(|| fix_batch(&batch, config)))
                            .unwrap_or_else(|_| Err(anyhow::anyhow!("fixmate thread panicked")));

                    if done_tx.send((id, result)).is_err() {
                        break;
                    }
                }
            });
        }

        drop(work_rx);
        drop(done_tx);

        let write = scope.spawn(move || -> anyhow::Result<()> {
            let mut pending = BTreeMap::new();
            let mut next = 0;

            for (id, result) in done_rx {
                pending.insert(id, result);

                // Errors surface in batch order, after the batches before them are written.
                while let Some(result) = pending.remove(&next) {
                    writer.write_bam_records(&result?)?;
                    next += 1;

                    // The batch is written: let another one in. The reader may have finished.
                    let _ = permits_tx.send(());
                }
            }

            Ok(())
        });

        // Sends a batch once fewer than `in_flight` are unwritten. The permits run out (the
        // channel closes) only when the writer thread stopped.
        let send = |id: u64, batch: Vec<u8>| -> anyhow::Result<()> {
            permits_rx
                .recv()
                .ok()
                .and_then(|()| work_tx.send((id, batch)).ok())
                .ok_or_else(|| anyhow::anyhow!("fixmate writer stopped"))
        };

        let mut batch = Vec::with_capacity(BATCH_CAPACITY);
        let mut last_name = Vec::new();
        let mut id = 0;

        let produced = reader.for_each(header, |record| {
            let crate::alignment::Rec::Bam(record) = record else {
                bail!("parallel fixmate needs BAM input");
            };

            let raw = record.as_ref();
            let name = crate::sort_order::name(raw);

            // Batches end only between templates.
            if batch.len() >= BATCH_CAPACITY && name != last_name.as_slice() {
                let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_CAPACITY));
                send(id, full)?;
                id += 1;
            }

            last_name.clear();
            last_name.extend_from_slice(name);
            crate::runs::write_raw_record(&mut batch, raw)?;

            Ok(())
        });

        let sent = match produced {
            Ok(()) if !batch.is_empty() => send(id, batch),
            other => other,
        };

        drop(work_tx);

        let written = write
            .join()
            .map_err(|_| anyhow::anyhow!("fixmate writer thread panicked"))?;

        // A writer error explains a failed send, so report it first.
        written.and(sent)
    })
}

/// Decodes, fixes, and re-encodes the templates of one batch.
fn fix_batch(batch: &[u8], config: &BatchConfig<'_>) -> anyhow::Result<Vec<u8>> {
    let mut reader = pasta::bam::io::Reader::from(batch);
    let mut encoder = pasta::bam::io::Writer::from(Vec::with_capacity(batch.len()));
    let mut template: Vec<RecordBuf> = Vec::new();

    let mut flush = |template: &mut Vec<RecordBuf>| -> anyhow::Result<()> {
        fix_template(template, config.options, config.header)?;

        for record in template.iter() {
            let removed = config.options.remove_reads
                && record
                    .flags()
                    .intersects(Flags::SECONDARY | Flags::UNMAPPED);

            if !removed {
                encoder.write_alignment_record(config.header, record)?;
            }
        }

        template.clear();
        Ok(())
    };

    loop {
        let mut record = RecordBuf::default();

        if reader.read_record_buf(config.header, &mut record)? == 0 {
            break;
        }

        sanitize_record(&mut record, config.lengths, config.sanitize)?;

        if template
            .first()
            .is_some_and(|first| first.name() != record.name())
        {
            flush(&mut template)?;
        }

        template.push(record);
    }

    if !template.is_empty() {
        flush(&mut template)?;
    }

    Ok(encoder.into_inner())
}

/// The 0-based position, or -1.
fn pos(record: &RecordBuf) -> i64 {
    record.alignment_start().map_or(-1, |p| {
        i64::try_from(usize::from(p)).unwrap_or(i64::MAX) - 1
    })
}

fn set_pos(record: &mut RecordBuf, pos: i64) {
    *record.alignment_start_mut() = usize::try_from(pos + 1).ok().and_then(Position::new);
}

fn tid(record: &RecordBuf) -> Option<usize> {
    record.reference_sequence_id()
}

fn reference_span(ops: &[Op]) -> usize {
    ops.iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum()
}

/// HTSlib's `bam_endpos`: the 0-based exclusive end, covering at least one base.
fn end_pos(record: &RecordBuf) -> i64 {
    let span = if record.flags().is_unmapped() {
        0
    } else {
        reference_span(record.cigar().as_ref())
    };

    pos(record) + i64::try_from(span.max(1)).unwrap_or(i64::MAX)
}

/// HTSlib's `bam_sanitize`.
fn sanitize_record(
    record: &mut RecordBuf,
    lengths: &[usize],
    sanitize: Sanitize,
) -> anyhow::Result<()> {
    if sanitize.has(Sanitize::POS) && tid(record).is_none() {
        set_pos(record, -1);

        if sanitize.has(Sanitize::UNMAP) {
            *record.flags_mut() |= Flags::UNMAPPED;
        }
    }

    if sanitize.has(Sanitize::CIGAR) && !record.flags().is_unmapped() {
        if pos(record) < 0 && sanitize.has(Sanitize::UNMAP) {
            *record.flags_mut() |= Flags::UNMAPPED;
        } else {
            let length = tid(record)
                .and_then(|id| lengths.get(id))
                .map_or(0, |&n| i64::try_from(n).unwrap_or(i64::MAX));

            if pos(record) >= length && sanitize.has(Sanitize::UNMAP) {
                *record.flags_mut() |= Flags::UNMAPPED;

                if sanitize.has(Sanitize::POS) {
                    *record.reference_sequence_id_mut() = None;
                    set_pos(record, -1);
                }
            } else if end_pos(record) > length {
                trim(record, length);
            }
        }
    }

    if record.flags().is_unmapped() {
        if sanitize.has(Sanitize::CIGAR) && !record.cigar().as_ref().is_empty() {
            record.cigar_mut().as_mut().clear();
        }

        if sanitize.has(Sanitize::MQUAL) {
            *record.mapping_quality_mut() = MappingQuality::new(0);
        }

        if sanitize.has(Sanitize::AUX) {
            for tag in UNMAPPED_REMOVED_TAGS {
                while record.data_mut().shift_remove(&tag).is_some() {}
            }
        }
    }

    if sanitize.has(Sanitize::CIGARX) && !record.flags().is_unmapped() {
        for op in record.cigar_mut().as_mut() {
            if matches!(op.kind(), Kind::SequenceMatch | Kind::SequenceMismatch) {
                *op = Op::new(Kind::Match, op.len());
            }
        }
    }

    if sanitize.has(Sanitize::CIGDUP) && !record.flags().is_unmapped() {
        merge_cigar_ops(record.cigar_mut().as_mut());
    }

    Ok(())
}

/// Merges adjacent operations of the same kind and drops empty ones (`FIX_CIGDUP`).
fn merge_cigar_ops(ops: &mut Vec<Op>) {
    let mut merged: Vec<Op> = Vec::with_capacity(ops.len());

    for &op in ops.iter() {
        if op.is_empty() {
            continue;
        }

        match merged.last_mut() {
            Some(last) if last.kind() == op.kind() => {
                let total = last.len() + op.len();

                if total > MAX_OP_LEN {
                    *last = Op::new(op.kind(), MAX_OP_LEN);
                    merged.push(Op::new(op.kind(), total - MAX_OP_LEN));
                } else {
                    *last = Op::new(op.kind(), total);
                }
            }
            _ => merged.push(op),
        }
    }

    *ops = merged;
}

/// HTSlib's `bam_trim`: soft clips the alignment past reference position `end`.
fn trim(record: &mut RecordBuf, end: i64) {
    let ops: Vec<Op> = record.cigar().as_ref().to_vec();
    let mut position = pos(record);
    let mut i = 0;
    let mut kind = Kind::Match;
    let mut len = 0;

    while i < ops.len() {
        kind = ops[i].kind();
        len = i64::try_from(ops[i].len()).unwrap_or(i64::MAX);

        if kind.consumes_reference() {
            position += len;

            if position > end {
                break;
            }
        }

        i += 1;
    }

    if i == ops.len() {
        return;
    }

    let mut kept: Vec<Op> = ops[..i].to_vec();

    let to_len = |n: i64| usize::try_from(n).unwrap_or(0);

    // The soft clip counts only read bases. HTSlib also adds the lengths of deletions and
    // reference skips, which gives a CIGAR longer than the sequence (`3M10D2M` 3 bases from the
    // end becomes `3M3D9S` for 5 bases).
    let read_len = |kind: Kind, len: usize| if kind.consumes_read() { len } else { 0 };

    let mut clip = if position - len < end {
        // The operation straddles the end: split it.
        kept.push(Op::new(kind, to_len(end - (position - len))));
        read_len(kind, to_len(position - end))
    } else if position > end {
        *record.flags_mut() |= Flags::UNMAPPED;
        *record.flags_mut() &= !Flags::PROPERLY_SEGMENTED;
        return;
    } else {
        read_len(kind, ops[i].len())
    };

    let mut hard_clips = Vec::new();

    for op in &ops[i + 1..] {
        if op.kind() == Kind::HardClip {
            hard_clips.push(*op);
        } else {
            clip += read_len(op.kind(), op.len());
        }
    }

    if clip > 0 {
        kept.push(Op::new(Kind::SoftClip, clip));
    }

    kept.extend(hard_clips);
    *record.cigar_mut() = kept.into_iter().collect();
}

fn cigar_string(record: &RecordBuf) -> String {
    let ops = record.cigar().as_ref();

    if ops.is_empty() {
        return "*".to_owned();
    }

    ops.iter()
        .map(|op| format!("{}{}", op.len(), char::from(kind_char(op.kind()))))
        .collect()
}

fn kind_char(kind: Kind) -> u8 {
    match kind {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    }
}

/// Removes a tag (keeping the order of the others) and appends a new value, as samtools'
/// `bam_aux_del` followed by `bam_aux_append` does.
fn replace_tag(record: &mut RecordBuf, tag: Tag, value: Value) {
    let data = record.data_mut();
    data.shift_remove(&tag);
    data.insert(tag, value);
}

/// `sync_mate_inner`
fn sync_mate_inner(src: &RecordBuf, dest: &mut RecordBuf) {
    *dest.mate_reference_sequence_id_mut() = tid(src);
    *dest.mate_alignment_start_mut() = src.alignment_start();

    let flags = dest.flags_mut();
    flags.set(
        Flags::MATE_REVERSE_COMPLEMENTED,
        src.flags().is_reverse_complemented(),
    );

    if src.flags().is_unmapped() {
        flags.insert(Flags::MATE_UNMAPPED);
    }
}

/// `sync_mq_mc`
fn sync_mq_mc(src: &RecordBuf, dest: &mut RecordBuf) {
    if !src.flags().is_unmapped() {
        let mq = src.mapping_quality().map_or(255, u8::from);
        replace_tag(dest, MATE_MAPPING_QUALITY, Value::Int32(i32::from(mq)));
    }

    if !src.flags().is_unmapped() || !dest.flags().is_unmapped() {
        replace_tag(
            dest,
            MATE_CIGAR,
            Value::String(BString::from(cigar_string(src))),
        );
    }
}

/// `sync_unmapped_pos_inner`
fn sync_unmapped_pos_inner(src: &RecordBuf, dest: &mut RecordBuf) {
    if dest.flags().is_unmapped() && !src.flags().is_unmapped() {
        *dest.reference_sequence_id_mut() = tid(src);
        *dest.alignment_start_mut() = src.alignment_start();
    }
}

/// `calc_mate_score`: the sum of base qualities of at least 15. Missing qualities count as 255,
/// as they are stored in BAM.
fn mate_score(record: &RecordBuf) -> u32 {
    let scores = record.quality_scores().as_ref();

    if scores.is_empty() {
        let n = u32::try_from(record.sequence().len()).unwrap_or(u32::MAX);
        return n.saturating_mul(255);
    }

    scores
        .iter()
        .filter(|&&q| q >= MIN_SCORED_QUALITY)
        .map(|&q| u32::from(q))
        .sum()
}

/// `plausibly_properly_paired`
fn plausibly_properly_paired(a: &RecordBuf, b: &RecordBuf) -> bool {
    if a.flags().is_unmapped() || b.flags().is_unmapped() || tid(a) != tid(b) {
        return false;
    }

    let five_prime = |r: &RecordBuf| {
        if r.flags().is_reverse_complemented() {
            end_pos(r)
        } else {
            pos(r)
        }
    };

    let (first, second) = if five_prime(a) > five_prime(b) {
        (b, a)
    } else {
        (a, b)
    };

    !first.flags().is_reverse_complemented() && second.flags().is_reverse_complemented()
}

/// `bam_template_cigar`: sets `ct` on the lower-positioned record.
fn template_cigar(a: &mut RecordBuf, b: &mut RecordBuf) {
    let unusable = |r: &RecordBuf| tid(r).is_none() || pos(r) < 0 || r.flags().is_unmapped();

    if tid(a) != tid(b) || unusable(a) || unusable(b) {
        return;
    }

    let swapped = pos(a) > pos(b);
    let (b1, b2): (&RecordBuf, &RecordBuf) = if swapped { (b, a) } else { (a, b) };

    let segment = |r: &RecordBuf| {
        if r.flags().is_first_segment() {
            '1'
        } else {
            '2'
        }
    };
    let strand = |r: &RecordBuf| {
        if r.flags().is_reverse_complemented() {
            'R'
        } else {
            'F'
        }
    };

    let mut s = String::new();
    s.push(segment(b1));
    s.push(strand(b1));
    for op in b1.cigar().as_ref() {
        s.push_str(&format!("{}{}", op.len(), char::from(kind_char(op.kind()))));
    }
    s.push_str(&(pos(b2) - end_pos(b1)).to_string());
    s.push('T');
    s.push(segment(b2));
    s.push(strand(b2));
    for op in b2.cigar().as_ref() {
        s.push_str(&format!("{}{}", op.len(), char::from(kind_char(op.kind()))));
    }

    a.data_mut().shift_remove(&TEMPLATE_CIGAR);
    b.data_mut().shift_remove(&TEMPLATE_CIGAR);

    let target = if swapped { b } else { a };
    target
        .data_mut()
        .insert(TEMPLATE_CIGAR, Value::String(BString::from(s)));
}

/// Returns mutable references to two distinct elements.
fn pair_mut(records: &mut [RecordBuf], i: usize, j: usize) -> (&mut RecordBuf, &mut RecordBuf) {
    assert_ne!(i, j, "a record cannot be its own mate");

    if i < j {
        let (left, right) = records.split_at_mut(j);
        (&mut left[i], &mut right[0])
    } else {
        let (left, right) = records.split_at_mut(i);
        (&mut right[0], &mut left[j])
    }
}

/// Fixes one template (records sharing a name), as `bam_mating_core` does.
fn fix_template(
    records: &mut [RecordBuf],
    options: Options,
    header: &sam::Header,
) -> anyhow::Result<()> {
    let not_primary = Flags::SECONDARY | Flags::SUPPLEMENTARY;

    let mut pre: Option<usize> = None;
    let mut cur: Option<usize> = None;
    let mut by_segment: [Option<usize>; 2] = [None, None];
    let mut pre_end = 0;

    for n in 0..records.len() {
        let flags = records[n].flags();

        if flags.intersects(not_primary) {
            continue;
        }

        let is_r2 = usize::from(flags.is_last_segment());

        if options.base_mods {
            base_mods::fix_primary(&mut records[n], header)?;
        }

        let Some(p) = pre else {
            pre = Some(n);
            by_segment[is_r2] = Some(n);
            pre_end = if flags.is_unmapped() {
                0
            } else {
                end_pos(&records[n])
            };
            continue;
        };

        cur = Some(n);
        by_segment[is_r2] = Some(n);
        let cur_end = if flags.is_unmapped() {
            0
        } else {
            end_pos(&records[n])
        };

        let (a, b) = pair_mut(records, p, n);

        *a.flags_mut() |= Flags::SEGMENTED;
        *b.flags_mut() |= Flags::SEGMENTED;

        // sync_mate
        sync_unmapped_pos_inner(a, b);
        sync_unmapped_pos_inner(b, a);
        sync_mate_inner(a, b);
        sync_mate_inner(b, a);
        sync_mq_mc(a, b);
        sync_mq_mc(b, a);

        let unmapped = Flags::UNMAPPED | Flags::MATE_UNMAPPED;

        if tid(a) == tid(b) && !b.flags().intersects(unmapped) && !a.flags().intersects(unmapped) {
            let cur5 = if b.flags().is_reverse_complemented() {
                cur_end
            } else {
                pos(b)
            };
            let pre5 = if a.flags().is_reverse_complemented() {
                pre_end
            } else {
                pos(a)
            };

            *b.template_length_mut() = i32::try_from(pre5 - cur5).context("TLEN overflows")?;
            *a.template_length_mut() = i32::try_from(cur5 - pre5).context("TLEN overflows")?;
        } else {
            *b.template_length_mut() = 0;
            *a.template_length_mut() = 0;
        }

        if options.template_cigar {
            template_cigar(a, b);
        }

        if options.proper_pair_check && !plausibly_properly_paired(a, b) {
            *a.flags_mut() &= !Flags::PROPERLY_SEGMENTED;
            *b.flags_mut() &= !Flags::PROPERLY_SEGMENTED;
        }

        if options.mate_score {
            let (score_a, score_b) = (mate_score(a), mate_score(b));
            replace_tag(b, MATE_SCORE, score_value(score_a));
            replace_tag(a, MATE_SCORE, score_value(score_b));
        }

        if options.remove_reads {
            let orphaned = Flags::MATE_REVERSE_COMPLEMENTED | Flags::PROPERLY_SEGMENTED;

            if a.flags().is_unmapped() {
                *b.flags_mut() &= !orphaned;
            }

            if b.flags().is_unmapped() {
                *a.flags_mut() &= !orphaned;
            }
        }
    }

    if let (None, Some(p)) = (cur, pre) {
        let record = &mut records[p];
        *record.mate_reference_sequence_id_mut() = None;
        *record.mate_alignment_start_mut() = None;
        *record.template_length_mut() = 0;
        *record.flags_mut() &= !(Flags::MATE_REVERSE_COMPLEMENTED | Flags::PROPERLY_SEGMENTED);
    }

    for n in 0..records.len() {
        let flags = records[n].flags();

        if !flags.intersects(not_primary) {
            continue;
        }

        // Supplementary alignments get their mate's primary information.
        if flags.is_supplementary() && flags.is_segmented() {
            let mate = by_segment[usize::from(!flags.is_last_segment())];

            if let Some(m) = mate {
                let (mate, record) = pair_mut(records, m, n);
                sync_mate_inner(mate, record);
                sync_mq_mc(mate, record);
            }
        }

        if options.base_mods {
            match by_segment[usize::from(flags.is_last_segment())] {
                Some(p) => {
                    let (primary, record) = pair_mut(records, p, n);
                    base_mods::fix_secondary(primary, record, header)?;
                }
                // Records with base modifications but no known primary are only validated.
                None => base_mods::validate(&mut records[n], header)?,
            }
        }
    }

    Ok(())
}

/// `fix_MM` and its helpers: base modification tags are checked against the sequence, and those
/// of secondary and supplementary alignments trimmed where hard clips removed bases.
mod base_mods {
    use anyhow::bail;
    use pasta::pileup::{Read, baq::SEQ_NT16_TABLE, mods::BaseMods};
    use pasta::sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{cigar::op::Kind, data::field::Tag},
            record_buf::data::field::{Value, value::Array},
        },
    };

    const MM: Tag = Tag::new(b'M', b'M');
    const ML: Tag = Tag::new(b'M', b'L');
    const MN: Tag = Tag::new(b'M', b'N');

    /// The tags `delete_mod_tags` removes.
    const MOD_TAGS: [Tag; 5] = [MM, Tag::new(b'M', b'm'), ML, Tag::new(b'M', b'l'), MN];

    /// `fix_MM` for a primary alignment.
    pub(super) fn fix_primary(record: &mut RecordBuf, header: &sam::Header) -> anyhow::Result<()> {
        rename_draft_tags(record);

        if record.data().get(&MM).is_none() {
            delete_tags(record);
            return Ok(());
        }

        let (end5, end3) = hard_clips(record);
        let mn = mn_value(record);
        let len = record.sequence().len();

        if end5 == 0 && end3 == 0 && mn <= 0 {
            // No MN tag, but also no clipping: the MM tag is assumed to be valid.
            if len > 0 {
                set_mn(record, len)?;
            }
        } else if (end5 > 0 || end3 > 0) && i64::try_from(len).ok() != Some(mn) {
            // The bases were hard clipped after the modifications were called.
            delete_tags(record);
        }

        validate(record, header)
    }

    /// `fix_MM` for a secondary or supplementary alignment with a known primary.
    pub(super) fn fix_secondary(
        primary: &RecordBuf,
        record: &mut RecordBuf,
        header: &sam::Header,
    ) -> anyhow::Result<()> {
        rename_draft_tags(record);

        if record.data().get(&MM).is_none() {
            delete_tags(record);
            return Ok(());
        }

        let (end5, end3) = hard_clips(record);
        let mn = mn_value(record);
        let len = record.sequence().len();
        let primary_len = primary.sequence().len();

        // An MN matching the sequence length means the tags are already valid.
        if i64::try_from(len).ok() == Some(mn) {
            return validate(record, header);
        }

        if primary_len != len + end3 + end5 {
            delete_tags(record);
            return Ok(());
        }

        if (end5 > 0 || end3 > 0) && (mn < 0 || i64::try_from(primary_len).ok() == Some(mn)) {
            match trim(primary, record, end5, end3, header)? {
                Trimmed::Yes => {}
                // The tags were called on the primary's sequence; ones that do not fit it are
                // removed, as `validate` removes invalid tags.
                Trimmed::Invalid => delete_tags(record),
                Trimmed::No if !describes(record, primary, header)? => delete_tags(record),
                Trimmed::No => {
                    let name = record
                        .name()
                        .map(|name| name.to_string())
                        .unwrap_or_default();
                    eprintln!(
                        "[W::fix_MM] cannot trim the base modification tags of {name} to its \
                         hard clips; keeping them for the {primary_len}-base primary sequence"
                    );
                    record.data_mut().shift_remove(&MN);
                    set_mn(record, primary_len)?;
                }
            }

            return Ok(());
        }

        // samtools ignores this failure (a non-integer MN) for these records, skipping the
        // validation.
        if len > 0 && set_mn(record, len).is_err() {
            return Ok(());
        }

        validate(record, header)
    }

    /// Removes the base modification tags if they do not parse or do not fit the sequence.
    pub(super) fn validate(record: &mut RecordBuf, header: &sam::Header) -> anyhow::Result<()> {
        if !is_valid(record, header)? {
            delete_tags(record);
        }

        Ok(())
    }

    /// Whether the base modification tags parse and fit the sequence, as HTSlib checks them.
    fn is_valid(record: &RecordBuf, header: &sam::Header) -> anyhow::Result<bool> {
        let mut writer = pasta::bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(header, record)?;
        let mut raw = writer.into_inner();
        raw.drain(..size_of::<u32>());

        let read = Read::new(raw)?;

        Ok(BaseMods::parse(&read)
            .and_then(|mut mods| mods.validate(&read))
            .is_ok())
    }

    /// Whether the MM and ML tags of `record` are valid for the sequence of `primary`.
    fn describes(
        record: &RecordBuf,
        primary: &RecordBuf,
        header: &sam::Header,
    ) -> anyhow::Result<bool> {
        let mut probe = primary.clone();
        let data = probe.data_mut();

        for tag in MOD_TAGS {
            data.shift_remove(&tag);
        }

        for tag in [MM, ML] {
            if let Some(value) = record.data().get(&tag) {
                data.insert(tag, value.clone());
            }
        }

        is_valid(&probe, header)
    }

    /// `find_tags`: the draft tags Mm, Ml, and Mn become MM, ML, and MN, in place.
    fn rename_draft_tags(record: &mut RecordBuf) {
        let upgrade = |tag: Tag| match *tag.as_ref() {
            [b'M', c @ (b'm' | b'l' | b'n')] => Tag::new(b'M', c.to_ascii_uppercase()),
            _ => tag,
        };

        let data = record.data_mut();

        if data.keys().all(|tag| upgrade(tag) == tag) {
            return;
        }

        let fields: Vec<(Tag, Value)> = data
            .iter()
            .map(|(tag, value)| (upgrade(tag), value.clone()))
            .collect();

        *data = fields.into_iter().collect();
    }

    /// `delete_mod_tags`
    fn delete_tags(record: &mut RecordBuf) {
        let data = record.data_mut();

        for tag in MOD_TAGS {
            data.shift_remove(&tag);
        }
    }

    /// The MN value, or -1 if it is absent or not an integer.
    fn mn_value(record: &RecordBuf) -> i64 {
        record.data().get(&MN).and_then(Value::as_int).unwrap_or(-1)
    }

    /// `bam_aux_update_int`: sets MN in place, or appends it, with HTSlib's smallest type.
    fn set_mn(record: &mut RecordBuf, len: usize) -> anyhow::Result<()> {
        let value = match u32::try_from(len)? {
            n if n < u32::from(u8::MAX) => Value::UInt8(u8::try_from(n)?),
            n if n < u32::from(u16::MAX) => Value::UInt16(u16::try_from(n)?),
            n => Value::UInt32(n),
        };

        let data = record.data_mut();

        match data.get_mut(&MN) {
            Some(old) if old.is_int() => *old = value,
            Some(_) => bail!("MN tag is not an integer"),
            None => {
                data.insert(MN, value);
            }
        }

        Ok(())
    }

    /// The 5' and 3' hard clip lengths.
    fn hard_clips(record: &RecordBuf) -> (usize, usize) {
        let ops = record.cigar().as_ref();

        let left = match ops.first() {
            Some(op) if op.kind() == Kind::HardClip => op.len(),
            _ => 0,
        };

        let right = match ops {
            [_, .., last] if last.kind() == Kind::HardClip => last.len(),
            _ => 0,
        };

        if record.flags().is_reverse_complemented() {
            (right, left)
        } else {
            (left, right)
        }
    }

    /// The result of [`trim`].
    enum Trimmed {
        /// The tags now describe the clipped sequence.
        Yes,
        /// The tags do not fit the primary's sequence they were called on.
        Invalid,
        /// The tags are left unchanged: they are outside the MM grammar, or trimming them did
        /// not give valid tags.
        No,
    }

    /// `trim_MM`: drops the MM deltas and ML likelihoods of the bases hard clipped from
    /// `record`, using its primary's sequence to count them, and sets MN to the new length.
    fn trim(
        primary: &RecordBuf,
        record: &mut RecordBuf,
        end5: usize,
        end3: usize,
        header: &sam::Header,
    ) -> anyhow::Result<Trimmed> {
        let counts = BaseCounts::new(primary, end5, end3);

        let data = record.data();

        let Some(Value::String(mm)) = data.get(&MM) else {
            return Ok(Trimmed::Invalid);
        };

        let ml = match data.get(&ML) {
            Some(Value::Array(Array::UInt8(ml))) => Some(&ml[..]),
            None => None,
            Some(_) => return Ok(Trimmed::Invalid),
        };

        let (new_mm, new_ml) = match trim_values(mm, ml, &counts) {
            Ok(values) => values,
            Err(TrimError::Invalid) => return Ok(Trimmed::Invalid),
            Err(TrimError::Unparsed) => return Ok(Trimmed::No),
        };

        let saved: Vec<(Tag, Value)> = [MM, ML, MN]
            .into_iter()
            .filter_map(|tag| data.get(&tag).map(|value| (tag, value.clone())))
            .collect();

        let data = record.data_mut();

        if let Some(value) = data.get_mut(&MM) {
            *value = Value::String(new_mm.into());
        }

        if let (Some(ml), Some(value)) = (new_ml, data.get_mut(&ML)) {
            *value = Value::Array(Array::UInt8(ml));
        }

        data.shift_remove(&MN);

        let len = record.sequence().len();

        if len > 0 {
            set_mn(record, len)?;
        }

        if is_valid(record, header)? {
            return Ok(Trimmed::Yes);
        }

        let data = record.data_mut();

        for (tag, value) in saved {
            data.insert(tag, value);
        }

        Ok(Trimmed::No)
    }

    /// The 4-bit code of the complement of each 4-bit base code (HTSlib's `seqi_rc`).
    const SEQI_RC: [usize; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

    /// The 4-bit code of `N`, which MM entries for `N` match on every base.
    const ANY_BASE: usize = 15;

    /// The number of bases of each 4-bit code in the primary's sequence, in its original
    /// (sequenced) orientation, before the 5' clip point, before the 3' one, and in all. The `N`
    /// counts are of every base, as in HTSlib.
    #[derive(Debug)]
    struct BaseCounts {
        before5: [i64; 16],
        before3: [i64; 16],
        total: [i64; 16],
    }

    impl BaseCounts {
        fn new(primary: &RecordBuf, end5: usize, end3: usize) -> Self {
            let bases = primary.sequence().as_ref();
            let reverse = primary.flags().is_reverse_complemented();

            let code = |i: usize| {
                if reverse {
                    SEQI_RC[usize::from(SEQ_NT16_TABLE[usize::from(bases[bases.len() - 1 - i])])]
                } else {
                    usize::from(SEQ_NT16_TABLE[usize::from(bases[i])])
                }
            };

            let count = |range: std::ops::Range<usize>, mut counts: [i64; 16]| {
                for i in range {
                    counts[code(i)] += 1;
                }

                counts
            };

            let end3 = bases.len().saturating_sub(end3).max(end5);
            let before5 = count(0..end5, [0; 16]);
            let before3 = count(end5..end3, before5);
            let total = count(end3..bases.len(), before3);

            let with_any = |mut counts: [i64; 16]| {
                counts[ANY_BASE] = counts.iter().sum();
                counts
            };

            Self {
                before5: with_any(before5),
                before3: with_any(before3),
                total: with_any(total),
            }
        }
    }

    /// One `MM` entry: its header (base, strand, codes, and `.`/`?`), the number of `ML` values
    /// per call (one per code, or one for a ChEBI number), and its deltas.
    struct Entry<'a> {
        header: &'a [u8],
        base: usize,
        stride: usize,
        deltas: Vec<i64>,
    }

    /// Parses one `MM` entry (without its `;`), or `None` if it is malformed.
    fn parse_entry(entry: &[u8]) -> Option<Entry<'_>> {
        let (header, deltas) = match entry.iter().position(|&b| b == b',') {
            Some(i) => (&entry[..i], Some(&entry[i + 1..])),
            None => (entry, None),
        };

        let [base, strand, codes @ ..] = header else {
            return None;
        };

        let base = match base {
            b'A' | b'C' | b'G' | b'T' | b'N' => usize::from(SEQ_NT16_TABLE[usize::from(*base)]),
            b'U' => usize::from(SEQ_NT16_TABLE[usize::from(b'T')]),
            _ => return None,
        };

        if !matches!(strand, b'+' | b'-') {
            return None;
        }

        let codes = match codes {
            [codes @ .., b'.' | b'?'] => codes,
            codes => codes,
        };

        let stride = if !codes.is_empty() && codes.iter().all(u8::is_ascii_digit) {
            1
        } else if !codes.is_empty() && codes.iter().all(u8::is_ascii_alphabetic) {
            codes.len()
        } else {
            return None;
        };

        let deltas = match deltas {
            Some(deltas) => deltas
                .split(|&b| b == b',')
                .map(|delta| {
                    if delta.is_empty() || !delta.iter().all(u8::is_ascii_digit) {
                        return None;
                    }

                    std::str::from_utf8(delta).ok()?.parse().ok()
                })
                .collect::<Option<Vec<i64>>>()?,
            None => Vec::new(),
        };

        Some(Entry {
            header,
            base,
            stride,
            deltas,
        })
    }

    /// Why MM and ML could not be trimmed.
    #[derive(Debug, PartialEq, Eq)]
    enum TrimError {
        /// MM does not follow the grammar of the SAM tags specification.
        Unparsed,
        /// ML does not have one value per modification code for each call, or a call is past
        /// the last base of its type in the primary's sequence (both errors for HTSlib).
        Invalid,
    }

    /// The MM and ML values after trimming.
    type TrimmedValues = (Vec<u8>, Option<Vec<u8>>);

    /// The trimming of `trim_MM` on the MM and ML values: the calls on bases outside the clip
    /// points are dropped with their `ML` values (one per modification code), and the first kept
    /// delta of each entry counts from the 5' clip point.
    fn trim_values(
        mm: &[u8],
        ml: Option<&[u8]>,
        counts: &BaseCounts,
    ) -> Result<TrimmedValues, TrimError> {
        let entries = match mm {
            [] => &[][..],
            [entries @ .., b';'] => entries,
            entries => entries,
        };

        let mut out = Vec::with_capacity(mm.len());
        let mut ml_out = Vec::new();
        let mut ml_pos = 0;

        for entry in entries.split(|&b| b == b';').filter(|_| !mm.is_empty()) {
            let entry = parse_entry(entry).ok_or(TrimError::Unparsed)?;
            let start = counts.before5[entry.base];
            let end = counts.before3[entry.base];

            out.extend_from_slice(entry.header);

            // The index among the bases of this type of the base each call is on.
            let mut index = -1;
            // The index of the last base before the kept calls.
            let mut last = start - 1;

            for delta in entry.deltas {
                index += delta + 1;

                if index >= counts.total[entry.base] {
                    return Err(TrimError::Invalid);
                }

                let values = match ml {
                    Some(ml) => Some(
                        ml.get(ml_pos..ml_pos + entry.stride)
                            .ok_or(TrimError::Invalid)?,
                    ),
                    None => None,
                };

                ml_pos += entry.stride;

                if (start..end).contains(&index) {
                    out.push(b',');
                    out.extend_from_slice((index - last - 1).to_string().as_bytes());
                    last = index;
                    ml_out.extend_from_slice(values.unwrap_or_default());
                }
            }

            out.push(b';');
        }

        if ml.is_some_and(|ml| ml.len() != ml_pos) {
            return Err(TrimError::Invalid);
        }

        Ok((out, ml.map(|_| ml_out)))
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// Counts of C bases (code 2) and of all bases (`N`): before the 5' clip point, before
        /// the 3' one, and in all.
        fn counts(c: [i64; 3], all: [i64; 3]) -> BaseCounts {
            let with = |i: usize| {
                let mut counts = [0; 16];
                counts[2] = c[i];
                counts[ANY_BASE] = all[i];
                counts
            };

            BaseCounts {
                before5: with(0),
                before3: with(1),
                total: with(2),
            }
        }

        #[test]
        fn test_trim_values() {
            let ml = Some(&[10, 20, 30][..]);
            let trim = |counts| trim_values(b"C+m,1,3,2;", ml, &counts);

            // Two C bases clipped at the 5' end, ten kept: the first call is dropped.
            let (mm, ml_out) = trim(counts([2, 12, 12], [2, 12, 12])).expect("trimmed");
            assert_eq!(mm, b"C+m,3,2;");
            assert_eq!(ml_out, Some(vec![20, 30]));

            // A clip point within a gap shortens it.
            let (mm, ml_out) = trim(counts([3, 12, 12], [3, 12, 12])).expect("trimmed");
            assert_eq!(mm, b"C+m,2,2;");
            assert_eq!(ml_out, Some(vec![20, 30]));

            // Bases clipped at the 3' end.
            let (mm, ml_out) = trim(counts([0, 5, 12], [0, 5, 12])).expect("trimmed");
            assert_eq!(mm, b"C+m,1;");
            assert_eq!(ml_out, Some(vec![10]));

            // A call past the primary's last C.
            assert_eq!(trim(counts([0, 5, 8], [0, 5, 8])), Err(TrimError::Invalid));

            // Without ML, and with every call clipped.
            let (mm, ml_out) =
                trim_values(b"C+m?,0,0;", None, &counts([2, 2, 2], [2, 2, 2])).expect("trimmed");
            assert_eq!(mm, b"C+m?;");
            assert_eq!(ml_out, None);
        }

        #[test]
        fn test_trim_values_with_combined_codes() {
            let counts = counts([2, 4, 4], [2, 4, 4]);

            // C+mh: ML holds the m and h likelihoods of each call, interleaved.
            let ml = [10, 20, 30, 40, 50, 60, 70, 80];
            let (mm, ml_out) = trim_values(b"C+mh,0,0,0,0;", Some(&ml), &counts).expect("trimmed");
            assert_eq!(mm, b"C+mh,0,0;");
            assert_eq!(ml_out, Some(vec![50, 60, 70, 80]));

            // A ChEBI code has one likelihood per call.
            let (mm, ml_out) =
                trim_values(b"C+76792.,2,0;", Some(&[10, 20]), &counts).expect("trimmed");
            assert_eq!(mm, b"C+76792.,0,0;");
            assert_eq!(ml_out, Some(vec![10, 20]));

            // Too few likelihoods for two codes per call.
            let short = [10, 20, 30, 40];
            assert_eq!(
                trim_values(b"C+mh,0,0,0,0;", Some(&short), &counts),
                Err(TrimError::Invalid)
            );
        }

        #[test]
        fn test_trim_values_with_several_entries() {
            // C+m? and C+h? each have their own deltas and likelihoods, in order; G-o has no
            // calls, and N+a counts every base. One C of four and three bases of ten are
            // clipped at the 5' end, and two Cs and six bases at the 3' end.
            let counts = counts([1, 4, 6], [3, 10, 16]);

            let mm = b"C+m?,0,1,1;C+h?,1,0;G-o;N+a,2,3;";
            let ml = [1, 2, 3, 4, 5, 6, 7];
            let (out, ml_out) = trim_values(mm, Some(&ml), &counts).expect("trimmed");
            assert_eq!(out, b"C+m?,1;C+h?,0,0;G-o;N+a,3;");
            assert_eq!(ml_out, Some(vec![2, 4, 5, 7]));

            // The likelihoods must cover every entry's calls, and no more.
            for ml in [&ml[..6], &[0; 8]] {
                assert_eq!(trim_values(mm, Some(ml), &counts), Err(TrimError::Invalid));
            }
        }

        /// Runs `fix_secondary` on the second record of a SAM file, after the first (its
        /// primary), and returns the second's MM, ML, and MN tags as SAM text.
        fn fix_pair(sam: &str) -> anyhow::Result<String> {
            let mut reader = sam::io::Reader::new(sam.as_bytes());
            let header = reader.read_header()?;
            let mut records: Vec<RecordBuf> =
                reader.record_bufs(&header).collect::<Result<_, _>>()?;
            let [primary, record] = &mut records[..] else {
                anyhow::bail!("expected two records");
            };

            fix_secondary(primary, record, &header)?;

            let mut writer = sam::io::Writer::new(Vec::new());
            writer.write_alignment_record(&header, record)?;
            let line = String::from_utf8(writer.into_inner())?;
            let tags: Vec<&str> = line.trim_end().split('\t').skip(11).collect();
            Ok(tags.join(" "))
        }

        const HEADER: &str = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n";

        #[test]
        fn test_fix_secondary_trims_hard_clipped_tags() -> anyhow::Result<()> {
            // One modification code: as samtools.
            let sam = format!(
                "{HEADER}r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\t\
                 MM:Z:C+m,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:4\n\
                 r\t2048\tchr1\t11\t60\t2H2M\t*\t0\t0\tCC\tII\t\
                 MM:Z:C+m,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:4\n"
            );
            assert_eq!(fix_pair(&sam)?, "MM:Z:C+m,0,0; ML:B:C,30,40 MN:i:2");

            // Two codes per call (samtools writes a corrupt record).
            let sam = format!(
                "{HEADER}r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\t\
                 MM:Z:C+mh,0,0,0,0;\tML:B:C,10,20,30,40,50,60,70,80\tMN:i:4\n\
                 r\t2048\tchr1\t11\t60\t2H2M\t*\t0\t0\tCC\tII\t\
                 MM:Z:C+mh,0,0,0,0;\tML:B:C,10,20,30,40,50,60,70,80\tMN:i:4\n"
            );
            assert_eq!(fix_pair(&sam)?, "MM:Z:C+mh,0,0; ML:B:C,50,60,70,80 MN:i:2");

            // Separate entries for each code.
            let sam = format!(
                "{HEADER}r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\t\
                 MM:Z:C+m?,0,0,0,0;C+h?,1,1;\tML:B:C,10,20,30,40,50,60\tMN:i:4\n\
                 r\t2048\tchr1\t11\t60\t2H2M\t*\t0\t0\tCC\tII\t\
                 MM:Z:C+m?,0,0,0,0;C+h?,1,1;\tML:B:C,10,20,30,40,50,60\tMN:i:4\n"
            );
            assert_eq!(
                fix_pair(&sam)?,
                "MM:Z:C+m?,0,0;C+h?,1; ML:B:C,30,40,60 MN:i:2"
            );

            // Reverse strand: the 5' clip is on the right, and the bases are counted in the
            // sequenced orientation (CCAACC; samtools counts the stored GGTTGG and drops every
            // call).
            let sam = format!(
                "{HEADER}r\t16\tchr1\t1\t60\t6M\t*\t0\t0\tGGTTGG\tIIIIII\t\
                 MM:Z:C+m,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:6\n\
                 r\t2064\tchr1\t11\t60\t4M2H\t*\t0\t0\tGGTT\tIIII\t\
                 MM:Z:C+m,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:6\n"
            );
            assert_eq!(fix_pair(&sam)?, "MM:Z:C+m,0,0; ML:B:C,30,40 MN:i:4");

            // An N entry counts every base (samtools counts only N bases and drops every call).
            let sam = format!(
                "{HEADER}r\t0\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tIIIIII\t\
                 MM:Z:N+n,0,0,0,0,0,0;\tML:B:C,1,2,3,4,5,6\tMN:i:6\n\
                 r\t2048\tchr1\t11\t60\t2H4M\t*\t0\t0\tGTAC\tIIII\t\
                 MM:Z:N+n,0,0,0,0,0,0;\tML:B:C,1,2,3,4,5,6\tMN:i:6\n"
            );
            assert_eq!(fix_pair(&sam)?, "MM:Z:N+n,0,0,0,0; ML:B:C,3,4,5,6 MN:i:4");

            // A forward supplementary of a reverse primary.
            let sam = format!(
                "{HEADER}r\t16\tchr1\t1\t60\t6M\t*\t0\t0\tGGTTGG\tIIIIII\t\
                 MM:Z:C+m,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:6\n\
                 r\t2048\tchr1\t11\t60\t2H4M\t*\t0\t0\tAACC\tIIII\t\
                 MM:Z:C+m,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:6\n"
            );
            assert_eq!(fix_pair(&sam)?, "MM:Z:C+m,0,0; ML:B:C,30,40 MN:i:4");

            Ok(())
        }

        #[test]
        fn test_fix_secondary_keeps_tags_it_cannot_trim() -> anyhow::Result<()> {
            // HTSlib accepts a delta written with a plus sign, so these tags are valid for the
            // primary, but they are outside the MM grammar and are not trimmed: they are kept,
            // with MN marking them as describing the primary's sequence.
            let sam = format!(
                "{HEADER}r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\n\
                 r\t2048\tchr1\t11\t60\t2H2M\t*\t0\t0\tCC\tII\t\
                 MM:Z:C+m,0,0,+0,0;\tML:B:C,10,20,30,40\n"
            );
            assert_eq!(
                fix_pair(&sam)?,
                "MM:Z:C+m,0,0,+0,0; ML:B:C,10,20,30,40 MN:i:4"
            );

            // Tags that do not fit the primary are invalid and removed.
            let sam = format!(
                "{HEADER}r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tCCCC\tIIII\n\
                 r\t2048\tchr1\t11\t60\t2H2M\t*\t0\t0\tCC\tII\t\
                 MM:Z:C+mh,0,0,0,0;\tML:B:C,10,20,30,40\tMN:i:4\n"
            );
            assert_eq!(fix_pair(&sam)?, "");

            Ok(())
        }

        #[test]
        fn test_trim_values_rejects_malformed_tags() {
            for mm in [
                &b"C+,0;"[..],
                b"C*m,0;",
                b"X+m,0;",
                b"C+m,a;",
                b"C+m,0,;",
                b"C+m,+1;",
            ] {
                assert_eq!(
                    trim_values(mm, None, &counts([0, 4, 4], [0, 4, 4])),
                    Err(TrimError::Unparsed),
                    "{:?}",
                    String::from_utf8_lossy(mm)
                );
            }
        }
    }
}

/// samtools appends mate scores as `i` (32-bit) values.
fn score_value(score: u32) -> Value {
    i32::try_from(score).map_or(Value::UInt32(score), Value::Int32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_cigar_ops() {
        let mut ops = vec![
            Op::new(Kind::Match, 5),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Insertion, 0),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Match, 4),
        ];
        merge_cigar_ops(&mut ops);
        assert_eq!(
            ops,
            [
                Op::new(Kind::Match, 8),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 4)
            ]
        );
    }

    #[test]
    fn test_trim() {
        use pasta::core::Position;

        let record = |start: usize, ops: &[Op]| {
            RecordBuf::builder()
                .set_alignment_start(Position::new(start).unwrap())
                .set_cigar(ops.iter().copied().collect())
                .build()
        };

        // 3M10D2M at 115 on a 120-base reference: the deletion crosses the end, and the soft
        // clip holds only the 2 remaining read bases.
        let mut r = record(
            115,
            &[
                Op::new(Kind::Match, 3),
                Op::new(Kind::Deletion, 10),
                Op::new(Kind::Match, 2),
            ],
        );
        trim(&mut r, 120);
        assert_eq!(
            r.cigar().as_ref(),
            [
                Op::new(Kind::Match, 3),
                Op::new(Kind::Deletion, 3),
                Op::new(Kind::SoftClip, 2),
            ]
        );

        // 10M at 195, then a hard clip.
        let mut r = record(195, &[Op::new(Kind::Match, 10), Op::new(Kind::HardClip, 4)]);
        trim(&mut r, 200);
        assert_eq!(
            r.cigar().as_ref(),
            [
                Op::new(Kind::Match, 6),
                Op::new(Kind::SoftClip, 4),
                Op::new(Kind::HardClip, 4),
            ]
        );
    }
}
